#include "fastcv/dim_reduction.h"
#include "fastcv/grm.h"
#include "fastcv/memory_estimate.h"
#include "fastcv/misc_util.h"

#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/Util/SelectionRule.h>

#ifdef FASTCV_USE_MKL
#include <mkl_lapack.h>
#include <mkl_types.h>
#endif

#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <numeric>
#include <vector>
#include <chrono>
#include <iostream>

#ifdef FASTCV_HAS_OPENMP
#include <omp.h>
#endif

namespace fastcv {

// ============================================================
// Private helpers
// ============================================================
void GrmPcaDimReduction::eigendecompose(Eigen::MatrixXd& G, int n, int k,
                                         Eigen::MatrixXd& U,
                                         Eigen::VectorXd& lambda) {
    auto do_eigen_sign_fix = [&]() {
        for (int c = 0; c < k; ++c) {
            int max_idx;
            U.col(c).cwiseAbs().maxCoeff(&max_idx);
            if (U(max_idx, c) < 0) U.col(c) = -U.col(c);
        }
    };

    int ncv = std::min(std::max(4 * k + 1, 20), n);
    std::cout << "    [EIG] Spectra k=" << k << " ncv=" << ncv << " n=" << n << "\n" << std::flush;
    {
        Spectra::DenseSymMatProd<double, Eigen::Upper> op(G);
        Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double, Eigen::Upper>> solver(op, k, ncv);
        solver.init();
        int nconv = solver.compute(Spectra::SortRule::LargestAlge, 1000);
        if (nconv >= k) {
            U = solver.eigenvectors(k);
            lambda = solver.eigenvalues();
            do_eigen_sign_fix();
        } else {
#ifdef FASTCV_USE_MKL
            std::cout << "    [EIG] Spectra converged only " << nconv
                      << ", falling back to MKL dsyevr\n" << std::flush;
            double vl = 0.0, vu = 0.0, abstol = 0.0;
            MKL_INT il = n - k + 1, iu = n;
            MKL_INT m_out = 0;
            double work_query;
            MKL_INT lwork = -1, liwork = -1, iwork_query, info = 0;
            MKL_INT mkl_n = n, mkl_ldz = n;
            Eigen::VectorXd sup_ev(n);
            Eigen::MatrixXd sup_z(n, k);
            std::vector<MKL_INT> isuppz(2 * n);

            dsyevr("V", "I", "U", &mkl_n, G.data(), &mkl_ldz,
                    &vl, &vu, &il, &iu, &abstol, &m_out,
                    sup_ev.data(), sup_z.data(), &mkl_ldz,
                    isuppz.data(), &work_query, &lwork, &iwork_query, &liwork, &info);
            if (info != 0) throw std::runtime_error("dsyevr workspace query failed (info=" + std::to_string(info) + ")");
            lwork = static_cast<MKL_INT>(work_query);
            liwork = iwork_query;
            std::vector<double> work(lwork);
            std::vector<MKL_INT> iwork(liwork);
            dsyevr("V", "I", "U", &mkl_n, G.data(), &mkl_ldz,
                    &vl, &vu, &il, &iu, &abstol, &m_out,
                    sup_ev.data(), sup_z.data(), &mkl_ldz,
                    isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
            if (info != 0) throw std::runtime_error("dsyevr failed (info=" + std::to_string(info) + ")");
            for (int i = 0; i < k; ++i) { lambda(i) = sup_ev(k-1-i); U.col(i) = sup_z.col(k-1-i); }
            do_eigen_sign_fix();
#else
            throw std::runtime_error("Spectra converged only " + std::to_string(nconv) + "/" + std::to_string(k));
#endif
        }
    }
    std::cout << "    [EIG] Eigendecomposition done\n" << std::flush;
}

void GrmPcaDimReduction::fill_maf_params(const Eigen::VectorXd& maf,
                                          DimReductionResult& result) {
    result.maf = maf;
    result.centers = 2.0 * maf;
    result.scales = (2.0 * maf.array() * (1.0 - maf.array())).sqrt();
    result.sigma2_sum = 2.0 * (maf.array() * (1.0 - maf.array())).sum();
    for (int j = 0; j < result.scales.size(); ++j) {
        if (result.scales(j) < 1e-10) result.scales(j) = 1.0;
    }
}

// ============================================================
// GrmPcaDimReduction: set_block_size
// ============================================================
void GrmPcaDimReduction::set_block_size(int block_size, int n_samples) {
    if (block_size > 0) {
        grm_block_size_ = block_size;
        std::cout << "    [GRM] Block size: " << grm_block_size_ << " (manual)\n" << std::flush;
        return;
    }

    long l3 = detect_l3_per_numa();
    grm_block_size_ = auto_detect_block_size(l3, n_samples);

    if (l3 > 0) {
        std::cout << "    [GRM] Block size: " << grm_block_size_
                  << " (auto, L3/NUMA=" << (l3 / (1024 * 1024)) << " MiB, n=" << n_samples << ")\n"
                  << std::flush;
    } else {
        std::cout << "    [GRM] Block size: " << grm_block_size_ << " (default, L3 detection failed)\n" << std::flush;
    }
}

// ============================================================
// GrmPcaDimReduction::fit()
// ============================================================
DimReductionResult GrmPcaDimReduction::fit(const PlinkReader& reader,
                                           const std::vector<int>& sample_idx,
                                           int n_components,
                                           int n_threads) {
    int n = static_cast<int>(sample_idx.size());
    int m = reader.n_snps();

    set_threads(n_threads);

    // Auto-detect block size on first call (grm_block_size_ == 0)
    if (grm_block_size_ == 0) {
        set_block_size(0, n);
    }

    // 1. Compute GRM (n×n), streaming SNP-by-SNP; also collect MAF
    Eigen::VectorXd maf;
    std::cout << "    [GRM] Computing " << n << "x" << n << " GRM from " << m << " SNPs...\n" << std::flush;
    auto t0 = std::chrono::high_resolution_clock::now();
    auto grm_res = compute_grm(reader, sample_idx, grm_block_size_, n_threads, &maf);
    auto t1 = std::chrono::high_resolution_clock::now();
    double grm_time = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "    [GRM] Done in " << grm_time << "s\n" << std::flush;

    // 2. Top-k eigendecomposition via Spectra (Arnoldi iteration)
    //    Spectra is O(n^2 * k * n_iter) — efficient when GRM doesn't fit in cache.
    //    MKL dsyevr (O(n^3) tridiagonal reduction) is used as fallback.
    int k = std::min(n_components, n - 1);
    if (k <= 0) {
        throw std::runtime_error("n_components must be > 0 and < n_samples");
    }

    Eigen::MatrixXd U(n, k);
    Eigen::VectorXd lambda(k);
    eigendecompose(grm_res.G, n, k, U, lambda);

    // 3. Build result (MAF already computed during GRM pass)
    DimReductionResult result;
    result.method_name = "grm_pca";
    result.n_components = k;
    result.eigenvectors = std::move(U);
    result.eigenvalues = std::move(lambda);
    fill_maf_params(maf, result);
    result.sigma2_sum = grm_res.sigma2_sum;

    return result;
}

// ============================================================
// GrmPcaDimReduction::transform()
// V-based streaming: test_pcs = (1/σ²) · Z_test · Z_train' · U · diag(1/√λ)
// ============================================================
Eigen::MatrixXd GrmPcaDimReduction::transform(const PlinkReader& reader,
                                               const std::vector<int>& new_idx,
                                               const DimReductionResult& params,
                                               const std::vector<int>& train_idx,
                                               int n_threads) const {
    set_threads(n_threads);

    return project_test_pcs(reader, new_idx, train_idx,
                            params.eigenvectors, params.eigenvalues,
                            params.sigma2_sum, grm_block_size_, n_threads);
}

// ============================================================
// GrmPcaDimReduction::fit_and_transform()
// Two-pass approach:
//   Pass 1: GRM + MAF (train only) → eigendecompose → U, λ
//   Pass 2: V-based streaming projection for test samples
// ============================================================
GrmPcaDimReduction::FitTransformResult GrmPcaDimReduction::fit_and_transform(
    const PlinkReader& reader,
    const std::vector<int>& train_idx,
    const std::vector<int>& test_idx,
    int n_components,
    int n_threads) const {

    int n_train = static_cast<int>(train_idx.size());
    int m = reader.n_snps();
    int k = std::min(n_components, n_train - 1);
    if (k <= 0) {
        throw std::runtime_error("n_components must be > 0 and < n_train");
    }

    set_threads(n_threads);

    // Auto-detect block size if not set
    int bs = grm_block_size_;
    if (bs <= 0) {
        long l3 = detect_l3_per_numa();
        bs = auto_detect_block_size(l3, n_train);
    }

    // === Pass 1: Compute GRM + MAF (train only) ===
    Eigen::VectorXd maf;
    std::cout << "    [GRM] Computing " << n_train << "x" << n_train
              << " GRM from " << m << " SNPs...\n" << std::flush;
    auto t0 = std::chrono::high_resolution_clock::now();
    auto grm_res = compute_grm(reader, train_idx, bs, n_threads, &maf);
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "    [GRM] Done in " << std::chrono::duration<double>(t1 - t0).count() << "s\n" << std::flush;

    double sigma2_sum = grm_res.sigma2_sum;

    // === Eigendecomposition ===
    Eigen::MatrixXd U(n_train, k);
    Eigen::VectorXd lambda(k);
    eigendecompose(grm_res.G, n_train, k, U, lambda);

    // Release G immediately — no longer needed after eigendecomposition
    grm_res.G = Eigen::MatrixXd();

    // === Build result ===
    FitTransformResult result;
    result.params.method_name = "grm_pca";
    result.params.n_components = k;
    result.params.eigenvectors = std::move(U);
    result.params.eigenvalues = std::move(lambda);
    result.params.sigma2_sum = sigma2_sum;
    result.train_pcs = result.params.eigenvectors;  // raw eigenvectors

    // === Pass 2: V-based test projection ===
    if (!test_idx.empty()) {
        std::cout << "    [V-proj] Projecting " << test_idx.size()
                  << " test samples via V-based streaming...\n" << std::flush;
        result.test_pcs = project_test_pcs(reader, test_idx, train_idx,
                                           result.params.eigenvectors,
                                           result.params.eigenvalues,
                                           sigma2_sum, bs, n_threads);
    }

    fill_maf_params(maf, result.params);
    return result;
}

// ============================================================
// GrmPcaDimReduction::fit_and_transform_gcross()
// Single-pass approach:
//   1. compute_grm_with_cross: G + G_cross + MAF in one SNP pass
//   2. Eigendecompose G → U, λ
//   3. test_pcs = G_cross * U * diag(1/λ)
// ============================================================
GrmPcaDimReduction::FitTransformResult GrmPcaDimReduction::fit_and_transform_gcross(
    const PlinkReader& reader,
    const std::vector<int>& train_idx,
    const std::vector<int>& test_idx,
    int n_components,
    int n_threads) const {

    int n_train = static_cast<int>(train_idx.size());
    int n_test = static_cast<int>(test_idx.size());
    int m = reader.n_snps();
    int k = std::min(n_components, n_train - 1);
    if (k <= 0) {
        throw std::runtime_error("n_components must be > 0 and < n_train");
    }

    set_threads(n_threads);

    // Auto-detect block size if not set
    int bs = grm_block_size_;
    if (bs <= 0) {
        long l3 = detect_l3_per_numa();
        bs = auto_detect_block_size(l3, n_train);
    }

    // === Single pass: GRM + G_cross + MAF ===
    std::cout << "    [GRM+Gx] Computing " << n_train << "x" << n_train << " GRM + "
              << n_test << "x" << n_train << " G_cross from " << m
              << " SNPs (single pass)...\n" << std::flush;
    auto t0 = std::chrono::high_resolution_clock::now();
    auto grm_result = compute_grm_with_cross(reader, train_idx, test_idx, bs, n_threads);
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "    [GRM+Gx] Single pass done in "
              << std::chrono::duration<double>(t1 - t0).count() << "s\n" << std::flush;

    // σ² from GRM pass (already computed correctly)
    double sigma2_sum = grm_result.sigma2_sum;

    // === Eigendecomposition of G ===
    Eigen::MatrixXd U(n_train, k);
    Eigen::VectorXd lambda(k);
    eigendecompose(grm_result.G, n_train, k, U, lambda);

    // Release G immediately
    grm_result.G = Eigen::MatrixXd();

    // === Compute test PCs: test_pcs = G_cross * U * diag(1/λ) ===
    Eigen::VectorXd inv_lambda = lambda.array().inverse();
    Eigen::MatrixXd test_pcs = grm_result.G_cross * U * inv_lambda.asDiagonal();

    // Release G_cross
    grm_result.G_cross = Eigen::MatrixXd();

    // === Build result ===
    FitTransformResult result;
    result.params.method_name = "grm_pca";
    result.params.n_components = k;
    result.params.eigenvectors = std::move(U);
    result.params.eigenvalues = std::move(lambda);
    result.params.sigma2_sum = sigma2_sum;
    result.train_pcs = result.params.eigenvectors;
    result.test_pcs = std::move(test_pcs);

    fill_maf_params(grm_result.maf, result.params);

    return result;
}

// ============================================================
// standard_pca_fit(): for small in-memory matrices (e.g., phenotype)
// ============================================================
StandardPcaResult standard_pca_fit(const Eigen::MatrixXd& X,
                                   int n_components,
                                   int n_threads) {
    int n = X.rows();
    int p = X.cols();
    int k = std::min(n_components, std::min(n, p));

    set_threads(n_threads);

    // Center
    Eigen::VectorXd centers(p);
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < p; ++j) {
        centers(j) = X.col(j).mean();
    }

    Eigen::MatrixXd Xc = X.rowwise() - centers.transpose();

    // Scale
    Eigen::VectorXd scales(p);
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < p; ++j) {
        double var = Xc.col(j).squaredNorm() / (n - 1);
        scales(j) = std::sqrt(var);
        if (scales(j) < 1e-10) scales(j) = 1.0;
        Xc.col(j) /= scales(j);
    }

    // SVD: Xc = U S V'
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Xc, Eigen::ComputeThinV);
    Eigen::MatrixXd rotation = svd.matrixV().leftCols(k);  // (p × k)

    StandardPcaResult result;
    result.centers = centers;
    result.scales = scales;
    result.rotation = rotation;
    return result;
}

// ============================================================
// standard_pca_transform()
// ============================================================
Eigen::MatrixXd standard_pca_transform(const Eigen::MatrixXd& X_new,
                                       const StandardPcaResult& params,
                                       int n_threads) {
    int n = X_new.rows();
    int p = X_new.cols();

    set_threads(n_threads);

    // Center and scale using stored parameters
    Eigen::MatrixXd Xc(n, p);
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < p; ++j) {
        Xc.col(j) = (X_new.col(j).array() - params.centers(j)) / params.scales(j);
    }

    // Project: Xc * V
    return Xc * params.rotation;
}

// ============================================================
// Factory
// ============================================================
std::unique_ptr<DimReduction> create_dim_reduction(const std::string& method) {
    if (method == "grm_pca" || method == "pca") {
        return std::make_unique<GrmPcaDimReduction>();
    }
    throw std::runtime_error("Unknown dimensionality reduction method: " + method +
                             ". Supported: grm_pca");
}

} // namespace fastcv
