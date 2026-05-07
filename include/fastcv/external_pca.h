#pragma once

#include "fastcv/types.h"
#include <string>
#include <vector>

namespace fastcv {

/// Result from external PCA tool (hiblup/plink).
/// PCs are in "natural" scale: train_pcs = U * diag(S),
/// matching the scale used internally by fastCV.
struct ExternalPcaResult {
    Eigen::MatrixXd train_pcs;    // n_train × k  (scaled PC scores)
    Eigen::VectorXd singular_values;  // k  (singular values of Z)
};

/// Run PCA via an external tool (hiblup or plink).
/// Writes a --keep file for training samples, calls the tool as a subprocess,
/// parses the output, and returns scaled PCs + singular values.
///
/// @param tool         "hiblup" or "plink"
/// @param tool_path    Path to the binary (empty = search PATH)
/// @param bed_prefix   PLINK file prefix (e.g., "/path/to/plink_test")
/// @param train_ids    Training sample FIDs (one per line in --keep file)
/// @param n_pcs        Number of PCs to compute
/// @param n_threads    Threads for the external tool
/// @param work_dir     Directory for temporary files
/// @param fold_idx     Fold index (for naming temp files)
ExternalPcaResult run_external_pca(
    const std::string& tool,
    const std::string& tool_path,
    const std::string& bed_prefix,
    const std::vector<std::string>& train_ids,
    int n_pcs,
    int n_threads,
    const std::string& work_dir,
    int fold_idx);

/// Compute test PCs using G_cross and external PCA results.
/// test_PCs = m * G_cross * P_unit * diag(1/S)
/// where P_unit = U (unit-norm eigenvectors from external tool).
///
/// @param G_cross        n_test × n_train cross-GRM (computed by our streaming code)
/// @param train_pcs      n_train × k (scaled: U*diag(S), from run_external_pca)
/// @param singular_values  k (S_i, from run_external_pca)
/// @param m_snps         Number of SNPs
/// @return               n_test × k test PCs in natural scale
Eigen::MatrixXd project_test_external(
    const Eigen::MatrixXd& G_cross,
    const Eigen::MatrixXd& train_pcs,
    const Eigen::VectorXd& singular_values,
    int m_snps);

} // namespace fastcv
