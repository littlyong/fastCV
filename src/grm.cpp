#include "fastcv/grm.h"
#include "fastcv/memory_estimate.h"
#include "fastcv/misc_util.h"

#ifdef FASTCV_HAS_OPENMP
#include <omp.h>
#endif

#ifdef FASTCV_USE_MKL
#include <mkl.h>
#elif defined(FASTCV_HAS_BLAS)
#include <cblas.h>
#endif

#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>

namespace fastcv {

// ============================================================
// Helper: directly call cblas_dgemm for GEMM, bypassing
// Eigen's BLAS wrapper which has overhead for dynamic matrices.
// C = alpha * A * B + beta * C  (row-major: C = alpha * A * B')
// ============================================================
static void my_dgemm(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE TransA,
                     CBLAS_TRANSPOSE TransB,
                     int M, int N, int K,
                     double alpha,
                     const double* A, int lda,
                     const double* B, int ldb,
                     double beta,
                     double* C, int ldc) {
#if defined(FASTCV_USE_MKL)
    cblas_dgemm(layout, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
#elif defined(FASTCV_HAS_BLAS)
    cblas_dgemm(layout, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
#else
    Eigen::Map<const Eigen::MatrixXd> eigA(A, (TransA == CblasNoTrans) ? M : K,
                                              (TransA == CblasNoTrans) ? K : M);
    Eigen::Map<const Eigen::MatrixXd> eigB(B, (TransB == CblasNoTrans) ? K : N,
                                              (TransB == CblasNoTrans) ? N : K);
    Eigen::Map<Eigen::MatrixXd> eigC(C, M, N);
    if (TransA == CblasNoTrans && TransB == CblasTrans)
        eigC.noalias() = alpha * eigA * eigB.transpose() + beta * eigC;
    else
        eigC.noalias() = alpha * eigA * eigB + beta * eigC;
#endif
}

// ============================================================
// compute_grm: streaming SNP-by-SNP with block multiplication
// GCTA-style: G = (1/Σ 2p(1-p)) * Σ_k (g-2p)(g-2p)'
// When n is small enough, uses OpenMP parallel SNP processing
// with per-thread G_local accumulators for speed.
// ============================================================
Eigen::MatrixXd compute_grm(const PlinkReader& reader,
                            const std::vector<int>& sample_idx,
                            int block_size,
                            int n_threads,
                            Eigen::VectorXd* maf_out,
                            const std::vector<int>* snp_idx) {
    int n = static_cast<int>(sample_idx.size());
    int m = reader.n_snps();

    int num_snps;
    if (snp_idx && !snp_idx->empty()) {
        num_snps = static_cast<int>(snp_idx->size());
    } else {
        num_snps = m;
    }

    set_threads(n_threads);

    // Decide parallel vs serial: use parallel when per-thread G_local is small
    int n_omp = get_effective_threads(n_threads);
    bool use_parallel = (num_snps == m) && (n_omp > 1) &&
                        (static_cast<long>(n_omp) * n * n * 8 < 8L * 1024 * 1024 * 1024);

    // For parallel mode: cap block_size to limit per-thread memory
    // Per-thread temp = n * bs * (8+4) bytes (Z_block + geno_flat), target < 200 MB
    int par_bs = block_size;
    if (use_parallel) {
        long max_bs = 200L * 1024 * 1024 / (static_cast<long>(n) * 12);
        par_bs = std::max(1000, std::min(par_bs, static_cast<int>(max_bs)));
    }

    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(n, n);
    if (maf_out) maf_out->resize(m);
    double sigma2_sum = 0.0;

    auto t0 = std::chrono::high_resolution_clock::now();
    int total_blocks = (num_snps + block_size - 1) / block_size;
    int par_total_blocks = use_parallel ? (num_snps + par_bs - 1) / par_bs : total_blocks;

    if (use_parallel) {
        // === Parallel SNP processing ===
        // G_local allocated inside parallel region (NUMA first-touch: each thread
        // initializes its own accumulator on the local NUMA node)
        std::vector<Eigen::MatrixXd> G_local(n_omp);
        std::vector<double> sigma2_local(n_omp, 0.0);

        std::cout << "      [GRM] Parallel mode: " << n_omp << " threads, bs="
                  << par_bs << ", " << par_total_blocks << " blocks\n" << std::flush;

        int blocks_done = 0;
        int print_interval = std::max(1, par_total_blocks / 20);

        #pragma omp parallel num_threads(n_omp)
        {
            int tid = omp_get_thread_num();
            // NUMA first-touch: each thread initializes its own accumulator
            G_local[tid] = Eigen::MatrixXd::Zero(n, n);
            sigma2_local[tid] = 0.0;
#ifdef FASTCV_USE_MKL
            int prev_mkl = mkl_set_num_threads_local(1);
#endif
            #pragma omp for schedule(dynamic, 1)
            for (int s = 0; s < num_snps; s += par_bs) {
                int end = std::min(s + par_bs, num_snps);
                int bs = end - s;

                Eigen::MatrixXd Z_block;
                Eigen::VectorXd maf_block;
                reader.read_snp_block(s, end, sample_idx, Z_block, maf_block, nullptr, true);

                for (int j = 0; j < bs; ++j)
                    sigma2_local[tid] += 2.0 * maf_block(j) * (1.0 - maf_block(j));

#if defined(FASTCV_USE_MKL) || defined(FASTCV_HAS_BLAS)
                cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans,
                            n, bs, 1.0, Z_block.data(), n,
                            1.0, G_local[tid].data(), n);
#else
                my_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                         n, n, bs, 1.0, Z_block.data(), n,
                         Z_block.data(), n, 1.0, G_local[tid].data(), n);
#endif

                if (maf_out)
                    for (int j = 0; j < bs; ++j)
                        (*maf_out)(s + j) = maf_block(j);

                int done;
                #pragma omp atomic capture
                done = blocks_done++;
                ++done;
                if (done % print_interval == 0 || done == par_total_blocks) {
                    #pragma omp critical(progress_bar)
                    {
                        detail::print_progress_bar(done * 100 / par_total_blocks);
                    }
                }
            }

#ifdef FASTCV_USE_MKL
            mkl_set_num_threads_local(prev_mkl);
#endif
        }

        // Sequential reduction
        for (int t = 0; t < n_omp; ++t) {
            G += G_local[t];
            sigma2_sum += sigma2_local[t];
        }
        // Release per-thread memory immediately
        G_local.clear();

    } else if (snp_idx && !snp_idx->empty()) {
        // Use LD-pruned SNP subset (serial)
        for (int b = 0; b < num_snps; b += block_size) {
            int bend = std::min(b + block_size, num_snps);
            int bs = bend - b;

            Eigen::MatrixXd Z_block;
            Eigen::VectorXd maf_block;
            Z_block.resize(n, bs);
            maf_block.resize(bs);
            for (int j = 0; j < bs; ++j) {
                Eigen::MatrixXd Z_col;
                Eigen::VectorXd maf_j;
                reader.read_snp_block((*snp_idx)[b + j], (*snp_idx)[b + j] + 1,
                                      sample_idx, Z_col, maf_j, nullptr, true);
                Z_block.col(j) = Z_col.col(0);
                maf_block(j) = maf_j(0);
            }

            for (int j = 0; j < bs; ++j) {
                double p = maf_block(j);
                sigma2_sum += 2.0 * p * (1.0 - p);
            }

#if defined(FASTCV_USE_MKL) || defined(FASTCV_HAS_BLAS)
            cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans,
                        n, bs, 1.0, Z_block.data(), n,
                        1.0, G.data(), n);
#else
            my_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                     n, n, bs, 1.0, Z_block.data(), n,
                     Z_block.data(), n, 1.0, G.data(), n);
#endif

            if (maf_out)
                for (int j = 0; j < bs; ++j)
                    (*maf_out)((*snp_idx)[b + j]) = maf_block(j);

            int cur_block = b / block_size + 1;
            detail::print_progress_bar(cur_block * 100 / total_blocks);
        }
    } else {
        // Serial mode (n too large for per-thread G_local)
        for (int s = 0; s < num_snps; s += block_size) {
            int end = std::min(s + block_size, num_snps);

            Eigen::MatrixXd Z_block;
            Eigen::VectorXd maf_block;
            reader.read_snp_block(s, end, sample_idx, Z_block, maf_block, nullptr, true);

            int bs = end - s;

            for (int j = 0; j < bs; ++j) {
                double p = maf_block(j);
                sigma2_sum += 2.0 * p * (1.0 - p);
            }

#if defined(FASTCV_USE_MKL) || defined(FASTCV_HAS_BLAS)
            cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans,
                        n, bs, 1.0, Z_block.data(), n,
                        1.0, G.data(), n);
#else
            my_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                     n, n, bs, 1.0, Z_block.data(), n,
                     Z_block.data(), n, 1.0, G.data(), n);
#endif

            int cur_block = s / block_size + 1;
            detail::print_progress_bar(cur_block * 100 / total_blocks);

            if (maf_out) {
                for (int j = 0; j < bs; ++j) (*maf_out)(s + j) = maf_block(j);
            }
        }
    }

    // Normalize by Σ 2p(1-p)
    G /= sigma2_sum;

#ifdef FASTCV_USE_MKL
    mkl_free_buffers();
#endif

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "\n      [GRM blocks] " << (use_parallel ? par_total_blocks : total_blocks)
              << " blocks done in " << std::chrono::duration<double>(t1 - t0).count() << "s\n"
              << std::flush;

    return G;
}

// ============================================================
// compute_grm_cross: G_cross = (1/Σ 2p(1-p)) * Σ D_new * D_train'
// Uses read_snp_block_dual: single fread per block, MAF from
// training set, prevents data leakage.
// ============================================================
Eigen::MatrixXd compute_grm_cross(const PlinkReader& reader,
                                  const std::vector<int>& sample_idx_new,
                                  const std::vector<int>& sample_idx_train,
                                  int block_size,
                                  int n_threads,
                                  double sigma2_sum) {
    int n_new = static_cast<int>(sample_idx_new.size());
    int n_train = static_cast<int>(sample_idx_train.size());
    int m = reader.n_snps();

    set_threads(n_threads);

    Eigen::MatrixXd G_cross = Eigen::MatrixXd::Zero(n_new, n_train);

    for (int s = 0; s < m; s += block_size) {
        int end = std::min(s + block_size, m);

        // Single I/O: read both train and new samples together (center only)
        Eigen::MatrixXd Z_train, Z_new;
        Eigen::VectorXd maf_block;
        reader.read_snp_block_dual(s, end, sample_idx_train, sample_idx_new,
                                   Z_train, Z_new, maf_block, true);

        int bs = end - s;

        // G_cross += Z_new_scaled * Z_train_scaled'
        my_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                 n_new, n_train, bs,
                 1.0,
                 Z_new.data(), n_new,
                 Z_train.data(), n_train,
                 1.0, G_cross.data(), n_new);
    }

    // Normalize by Σ 2p(1-p)
    G_cross /= sigma2_sum;

#ifdef FASTCV_USE_MKL
    mkl_free_buffers();
#endif

    return G_cross;
}

// ============================================================
// compute_grm_with_cross: single-pass GRM + G_cross + MAF
// OpenMP parallel when n is small enough for per-thread accumulators.
// ============================================================
GrmAndCrossResult compute_grm_with_cross(
    const PlinkReader& reader,
    const std::vector<int>& train_idx,
    const std::vector<int>& test_idx,
    int block_size,
    int n_threads,
    const std::vector<int>* snp_idx)
{
    int n_train = static_cast<int>(train_idx.size());
    int n_test = static_cast<int>(test_idx.size());
    int m = reader.n_snps();
    bool has_test = (n_test > 0);

    int num_snps;
    if (snp_idx && !snp_idx->empty()) {
        num_snps = static_cast<int>(snp_idx->size());
    } else {
        num_snps = m;
    }

    set_threads(n_threads);

    // Decide parallel vs serial: parallel when per-thread accumulators fit in memory
    int n_omp = get_effective_threads(n_threads);
    long per_thread_bytes = static_cast<long>(n_train) * n_train * 8;
    if (has_test) per_thread_bytes += static_cast<long>(n_test) * n_train * 8;
    bool use_parallel = (num_snps == m) && (n_omp > 1) &&
                        (static_cast<long>(n_omp) * per_thread_bytes < 8L * 1024 * 1024 * 1024);

    // Cap block_size for parallel mode to limit per-thread temp memory
    int par_bs = block_size;
    if (use_parallel) {
        long max_bs = 200L * 1024 * 1024 / (static_cast<long>(n_train + n_test) * 12);
        par_bs = std::max(1000, std::min(par_bs, static_cast<int>(max_bs)));
    }

    GrmAndCrossResult result;
    result.G = Eigen::MatrixXd::Zero(n_train, n_train);
    if (has_test) result.G_cross = Eigen::MatrixXd::Zero(n_test, n_train);
    result.maf.resize(m);
    double sigma2_sum = 0.0;

    auto t0 = std::chrono::high_resolution_clock::now();
    int total_blocks = (num_snps + block_size - 1) / block_size;
    int par_total_blocks = use_parallel ? (num_snps + par_bs - 1) / par_bs : total_blocks;

    if (use_parallel) {
        // === Parallel SNP processing ===
        // NUMA first-touch: each thread initializes its own accumulator inside
        // the parallel region so memory is allocated on the local NUMA node
        std::vector<Eigen::MatrixXd> G_local(n_omp);
        std::vector<Eigen::MatrixXd> Gx_local(has_test ? n_omp : 0);
        std::vector<double> sigma2_local(n_omp, 0.0);

        std::cout << "      [GRM+Gx] Parallel mode: " << n_omp << " threads, bs="
                  << par_bs << ", " << par_total_blocks << " blocks\n" << std::flush;

        int blocks_done = 0;
        int print_interval = std::max(1, par_total_blocks / 20);

        #pragma omp parallel num_threads(n_omp)
        {
            int tid = omp_get_thread_num();
            G_local[tid] = Eigen::MatrixXd::Zero(n_train, n_train);
            if (has_test) Gx_local[tid] = Eigen::MatrixXd::Zero(n_test, n_train);
            sigma2_local[tid] = 0.0;
#ifdef FASTCV_USE_MKL
            int prev_mkl = mkl_set_num_threads_local(1);
#endif
            #pragma omp for schedule(dynamic, 1)
            for (int s = 0; s < num_snps; s += par_bs) {
                int end = std::min(s + par_bs, num_snps);
                int bs = end - s;

                Eigen::MatrixXd Z_train, Z_test;
                Eigen::VectorXd maf_block;

                if (has_test) {
                    reader.read_snp_block_dual(s, end, train_idx, test_idx,
                                               Z_train, Z_test, maf_block, true);
                } else {
                    reader.read_snp_block(s, end, train_idx, Z_train, maf_block, nullptr, true);
                }

                for (int j = 0; j < bs; ++j)
                    sigma2_local[tid] += 2.0 * maf_block(j) * (1.0 - maf_block(j));

#if defined(FASTCV_USE_MKL) || defined(FASTCV_HAS_BLAS)
                cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans,
                            n_train, bs, 1.0,
                            Z_train.data(), n_train,
                            1.0, G_local[tid].data(), n_train);
#else
                my_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                         n_train, n_train, bs, 1.0,
                         Z_train.data(), n_train, Z_train.data(), n_train,
                         1.0, G_local[tid].data(), n_train);
#endif

                if (has_test) {
                    my_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                             n_test, n_train, bs, 1.0,
                             Z_test.data(), n_test, Z_train.data(), n_train,
                             1.0, Gx_local[tid].data(), n_test);
                }

                // MAF: each thread writes to disjoint indices — no race
                for (int j = 0; j < bs; ++j)
                    result.maf(s + j) = maf_block(j);

                int done;
                #pragma omp atomic capture
                done = blocks_done++;
                ++done;
                if (done % print_interval == 0 || done == par_total_blocks) {
                    #pragma omp critical(progress_bar)
                    {
                        detail::print_progress_bar(done * 100 / par_total_blocks);
                    }
                }
            }

#ifdef FASTCV_USE_MKL
            mkl_set_num_threads_local(prev_mkl);
#endif
        }

        // Sequential reduction
        for (int t = 0; t < n_omp; ++t) {
            result.G += G_local[t];
            if (has_test) result.G_cross += Gx_local[t];
            sigma2_sum += sigma2_local[t];
        }
        G_local.clear();
        Gx_local.clear();

    } else {
        // === Serial mode ===
        for (int s = 0; s < num_snps; s += block_size) {
            int end = std::min(s + block_size, num_snps);
            int bs = end - s;

            Eigen::MatrixXd Z_train, Z_test;
            Eigen::VectorXd maf_block;

            if (has_test) {
                reader.read_snp_block_dual(s, end, train_idx, test_idx,
                                           Z_train, Z_test, maf_block, true);
            } else {
                reader.read_snp_block(s, end, train_idx, Z_train, maf_block, nullptr, true);
            }

            for (int j = 0; j < bs; ++j) {
                double p = maf_block(j);
                sigma2_sum += 2.0 * p * (1.0 - p);
            }

#if defined(FASTCV_USE_MKL) || defined(FASTCV_HAS_BLAS)
            cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans,
                        n_train, bs, 1.0,
                        Z_train.data(), n_train,
                        1.0, result.G.data(), n_train);
#else
            my_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                     n_train, n_train, bs, 1.0,
                     Z_train.data(), n_train, Z_train.data(), n_train,
                     1.0, result.G.data(), n_train);
#endif

            if (has_test) {
                my_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                         n_test, n_train, bs, 1.0,
                         Z_test.data(), n_test, Z_train.data(), n_train,
                         1.0, result.G_cross.data(), n_test);
            }

            for (int j = 0; j < bs; ++j) result.maf(s + j) = maf_block(j);

            int cur_block = s / block_size + 1;
            detail::print_progress_bar(cur_block * 100 / total_blocks);
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "\n      [GRM+Gx] " << (use_parallel ? par_total_blocks : total_blocks)
              << " blocks done in " << std::chrono::duration<double>(t1 - t0).count() << "s\n"
              << std::flush;

    // Normalize by Σ 2p(1-p)
    result.G /= sigma2_sum;
    if (has_test) result.G_cross /= sigma2_sum;

#ifdef FASTCV_USE_MKL
    mkl_free_buffers();
#endif

    return result;
}

// ============================================================
// project_test_pcs: V-based streaming test projection
// test_pcs = (1/σ²) · Σ_j z_test_j · (z_train_j' · W)
// where W = U · diag(1/λ)
// Uses OpenMP parallel when n_train is small enough.
// ============================================================
Eigen::MatrixXd project_test_pcs(const PlinkReader& reader,
                                 const std::vector<int>& test_idx,
                                 const std::vector<int>& train_idx,
                                 const Eigen::MatrixXd& U,
                                 const Eigen::VectorXd& lambda,
                                 double sigma2_sum,
                                 int block_size,
                                 int n_threads,
                                 const std::vector<int>* snp_idx) {
    int n_test = static_cast<int>(test_idx.size());
    int n_train = static_cast<int>(train_idx.size());
    int k = static_cast<int>(lambda.size());
    int m = reader.n_snps();

    int num_snps;
    if (snp_idx && !snp_idx->empty()) num_snps = static_cast<int>(snp_idx->size());
    else num_snps = m;

    if (block_size <= 0) block_size = 1000;

    set_threads(n_threads);

    // W = U · diag(1/λ)   (n_train × k)
    Eigen::VectorXd inv_lambda = lambda.array().inverse();
    Eigen::MatrixXd W = U * inv_lambda.asDiagonal();

    Eigen::MatrixXd test_pcs = Eigen::MatrixXd::Zero(n_test, k);

    auto t0 = std::chrono::high_resolution_clock::now();
    int total_blocks = (num_snps + block_size - 1) / block_size;

    // Parallel when no snp_idx and enough threads
    int n_omp = get_effective_threads(n_threads);
    bool use_parallel = (num_snps == m) && (n_omp > 1) && !(snp_idx && !snp_idx->empty());

    // Cap block_size for parallel mode to limit per-thread memory
    // Per-thread persistent: par_bs * (bps + (n_train+n_test)*12 + k*8 + 8)
    // bps=raw_buf, (n_train+n_test)*12=geno_a/b(int4)+Z_a/b(double8), k*8=V_block
    // Target ~50 MB/thread — V-proj compute is light (dgemm K=k≈5), no need for large blocks
    int par_bs = block_size;
    if (use_parallel) {
        int bps = reader.bytes_per_snp();
        long per_snp_mem = static_cast<long>(bps)
                         + static_cast<long>(n_train + n_test) * 12
                         + static_cast<long>(k) * 8 + 8;
        long max_bs = 50L * 1024 * 1024 / per_snp_mem;
        par_bs = std::max(1000, std::min(par_bs, static_cast<int>(max_bs)));
    }

    if (use_parallel) {
        int par_total_blocks = (num_snps + par_bs - 1) / par_bs;
        // Per-thread accumulators (test_pcs_local is tiny: n_test × k)
        std::vector<Eigen::MatrixXd> pcs_local(n_omp);

        std::cout << "      [V-proj] Parallel mode: " << n_omp << " threads, bs="
                  << par_bs << ", " << par_total_blocks << " blocks\n" << std::flush;

        int blocks_done = 0;
        int print_interval = std::max(1, par_total_blocks / 20);

        #pragma omp parallel num_threads(n_omp)
        {
            int tid = omp_get_thread_num();
            // NUMA first-touch: each thread initializes its own accumulator
            pcs_local[tid] = Eigen::MatrixXd::Zero(n_test, k);
#ifdef FASTCV_USE_MKL
            int prev_mkl = mkl_set_num_threads_local(1);
#endif

            // Pre-allocate per-thread buffers (reused across all blocks, zero malloc in loop)
            int bps = reader.bytes_per_snp();
            DualReadBuf buf;
            buf.raw_buf.resize(static_cast<size_t>(par_bs) * bps);
            buf.geno_a.resize(static_cast<size_t>(n_train) * par_bs);
            buf.geno_b.resize(static_cast<size_t>(n_test) * par_bs);
            buf.loc_a.resize(n_train);
            buf.loc_b.resize(n_test);
            for (int i = 0; i < n_train; ++i) {
                int si = train_idx[i];
                buf.loc_a[i] = {si / 4, (si % 4) * 2};
            }
            for (int i = 0; i < n_test; ++i) {
                int si = test_idx[i];
                buf.loc_b[i] = {si / 4, (si % 4) * 2};
            }

            Eigen::MatrixXd Z_train(n_train, par_bs);
            Eigen::MatrixXd Z_test(n_test, par_bs);
            Eigen::VectorXd maf_block(par_bs);
            Eigen::MatrixXd V_block(par_bs, k);

            #pragma omp for schedule(dynamic, 1)
            for (int s = 0; s < num_snps; s += par_bs) {
                int end = std::min(s + par_bs, num_snps);
                int bs = end - s;

                reader.read_snp_block_dual_buf(s, end, train_idx, test_idx,
                                               Z_train, Z_test, maf_block, buf, true);

                // V_block = Z_train' · W   (bs × k)
                my_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                         bs, k, n_train,
                         1.0,
                         Z_train.data(), n_train,
                         W.data(), n_train,
                         0.0, V_block.data(), par_bs);

                // test_pcs_local += Z_test · V_block
                my_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                         n_test, k, bs,
                         1.0,
                         Z_test.data(), n_test,
                         V_block.data(), par_bs,
                         1.0, pcs_local[tid].data(), n_test);

                int done;
                #pragma omp atomic capture
                done = blocks_done++;
                ++done;
                if (done % print_interval == 0 || done == par_total_blocks) {
                    #pragma omp critical(progress_bar)
                    {
                        detail::print_progress_bar(done * 100 / par_total_blocks);
                    }
                }
            }

#ifdef FASTCV_USE_MKL
            mkl_set_num_threads_local(prev_mkl);
#endif
        }

        for (int t = 0; t < n_omp; ++t)
            test_pcs += pcs_local[t];
        pcs_local.clear();
    } else {
        for (int s = 0; s < num_snps; s += block_size) {
            int end = std::min(s + block_size, num_snps);
            int bs = end - s;

            Eigen::MatrixXd Z_train, Z_test;
            Eigen::VectorXd maf_block;
            reader.read_snp_block_dual(s, end, train_idx, test_idx,
                                       Z_train, Z_test, maf_block, true);

            Eigen::MatrixXd V_block(bs, k);
            my_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                     bs, k, n_train,
                     1.0,
                     Z_train.data(), n_train,
                     W.data(), n_train,
                     0.0, V_block.data(), bs);

            my_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                     n_test, k, bs,
                     1.0,
                     Z_test.data(), n_test,
                     V_block.data(), bs,
                     1.0, test_pcs.data(), n_test);

            int cur_block = s / block_size + 1;
            detail::print_progress_bar(cur_block * 100 / total_blocks);
        }
    }

    // Normalize by σ²
    test_pcs /= sigma2_sum;

#ifdef FASTCV_USE_MKL
    mkl_free_buffers();
#endif

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "\n      [V-proj] " << total_blocks
              << " blocks done in " << std::chrono::duration<double>(t1 - t0).count() << "s\n"
              << std::flush;

    return test_pcs;
}

} // namespace fastcv
