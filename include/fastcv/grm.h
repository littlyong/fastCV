#pragma once

#include "fastcv/types.h"

namespace fastcv {

/// Compute GRM G = Z Z' / m for given sample indices, streaming SNP-by-SNP.
/// Uses block matrix multiplication for efficiency.
/// @param reader    PLINK reader (must be open)
/// @param sample_idx  Sample indices to include
/// @param block_size  Number of SNPs per block (default 1000)
/// @param n_threads   OpenMP threads (0 = auto)
/// @param[out] maf_out  If non-null, filled with per-SNP MAF (m,)
/// @param snp_idx  If non-null, only use these SNP indices (for LD-pruned GRM)
/// @return n × n GRM matrix
Eigen::MatrixXd compute_grm(const PlinkReader& reader,
                            const std::vector<int>& sample_idx,
                            int block_size = 1000,
                            int n_threads = 0,
                            Eigen::VectorXd* maf_out = nullptr,
                            const std::vector<int>* snp_idx = nullptr);

/// Compute cross-GRM G_cross = Z_new Z_train' / m for two sample sets.
/// Used for projecting test samples via G_cross · U · diag(sqrt(m/lambda)).
/// @param reader       PLINK reader (must be open)
/// @param sample_idx_new    New (test) sample indices
/// @param sample_idx_train  Training sample indices
/// @param block_size   Number of SNPs per block (default 1000)
/// @param n_threads    OpenMP threads (0 = auto)
/// @param sigma2_sum   Σ 2p(1-p) from training set, used as normalizer
/// @return n_new × n_train cross-GRM matrix
Eigen::MatrixXd compute_grm_cross(const PlinkReader& reader,
                                  const std::vector<int>& sample_idx_new,
                                  const std::vector<int>& sample_idx_train,
                                  int block_size = 1000,
                                  int n_threads = 0,
                                  double sigma2_sum = 0.0);

/// Combined GRM + cross-GRM in a single SNP-by-SNP pass.
/// Reads SNPs once, computes both G = Z_train Z_train' / m and
/// G_cross = Z_test Z_train' / m simultaneously. MAF from training set.
/// @param reader       PLINK reader (must be open)
/// @param train_idx    Training sample indices
/// @param test_idx     Test sample indices (empty = no cross-GRM)
/// @param block_size   Number of SNPs per block
/// @param n_threads    OpenMP threads (0 = auto)
struct GrmAndCrossResult {
    Eigen::MatrixXd G;        // n_train × n_train
    Eigen::MatrixXd G_cross;  // n_test × n_train (empty if test_idx empty)
    Eigen::VectorXd maf;      // per-SNP MAF (m,)
};

GrmAndCrossResult compute_grm_with_cross(
    const PlinkReader& reader,
    const std::vector<int>& train_idx,
    const std::vector<int>& test_idx,
    int block_size,
    int n_threads,
    const std::vector<int>* snp_idx = nullptr);

/// Project test samples onto training PCs via V-based streaming.
/// Avoids forming the n_test × n_train G_cross matrix entirely.
///
/// Computes: test_pcs = (1/σ²) · Z_test · Z_train' · U · diag(1/√λ)
/// Equivalently, for each SNP block:
///   W = U · diag(1/√λ)           (n_train × k, precomputed)
///   V_block = Z_train' · W        (bs × k)
///   test_pcs += Z_test · V_block  (n_test × k)
/// test_pcs /= σ²
///
/// Complexity per SNP: O(n_train·k + n_test·k) vs O(n_test·n_train) for G_cross.
/// Memory: O(n_test·k) vs O(n_test·n_train) for G_cross.
/// @param reader       PLINK reader (must be open)
/// @param test_idx     Test sample indices
/// @param train_idx    Training sample indices (MAF computed from these)
/// @param U            Eigenvectors from GRM eigendecomposition (n_train × k)
/// @param lambda       Eigenvalues (k)
/// @param sigma2_sum   Σ 2p(1-p) from training set
/// @param block_size   SNPs per block (0 → default 1000)
/// @param n_threads    OpenMP threads (0 = auto)
/// @return n_test × k projected PCs
Eigen::MatrixXd project_test_pcs(const PlinkReader& reader,
                                 const std::vector<int>& test_idx,
                                 const std::vector<int>& train_idx,
                                 const Eigen::MatrixXd& U,
                                 const Eigen::VectorXd& lambda,
                                 double sigma2_sum,
                                 int block_size = 1000,
                                 int n_threads = 0,
                                 const std::vector<int>* snp_idx = nullptr);

} // namespace fastcv
