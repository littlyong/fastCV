#pragma once

#include "fastcv/types.h"

namespace fastcv {

/// Lightweight context for per-fold export (avoids accumulating full CvResult).
struct FoldExportContext {
    std::string output_dir;
    std::vector<std::string> sample_ids;
    std::vector<int> fold_ids;
    std::vector<std::string> snp_ids;
    std::vector<std::string> cov_col_names;
    std::vector<std::string> trait_names;
    int n_geno_pcs = 0;
    int n_pheno_pcs = 0;
    CorrectionMode correction_mode = CorrectionMode::GENO_ONLY;
    std::string separator = " ";
};

/// Export global metadata: fold_assignments.csv + config.json.
/// Call once before the fold loop.
void export_cv_meta(const FoldExportContext& ctx);

/// Export a single fold's results immediately.
/// Call inside the fold loop, then release fold memory.
/// If save_details=false, only exports y_train_residual.txt and y_test_residual.txt.
void export_single_fold(const FoldExportContext& ctx,
                        int fold_idx,
                        const SharedFoldData& shared,
                        const std::vector<FoldData>& trait_folds,
                        bool save_details = true);

/// Export outer test residuals for nested CV mode.
void export_test_residuals(const FoldExportContext& ctx,
                            int fold_idx,
                            const std::vector<int>& test_idx,
                            const std::vector<FoldData>& trait_folds);

/// Export nested CV results for one outer fold.
/// Always outputs train_val_residuals.txt and nested_fold_ids.txt.
/// If save_details=true, also outputs nested/nested_XX/ subdirectories with
/// per-nested-fold MAF, beta, PCs, eigenvalues.
void export_nested_cv_results(
    const FoldExportContext& ctx,
    int outer_fold_idx,
    const std::vector<int>& train_pool_idx,  // indices into ctx.sample_ids
    const std::vector<int>& nested_fold_ids,  // per-sample nested fold (1-based)
    const Eigen::MatrixXd& residual_matrix,   // (n_train_pool × n_nested_folds)
    bool save_details,
    const Eigen::VectorXd& outer_eigenvalues,
    const Eigen::VectorXd& outer_maf,
    // per-nested-fold details (save_details=true):
    const std::vector<SharedFoldData>& nested_shared,
    const std::vector<std::vector<FoldData>>& nested_trait_folds);

/// Export a single nested fold's detail files (immediate per-fold export).
/// Outputs nested/nested_XX/ subdirectory with MAF, beta, PCs, eigenvalues.
void export_nested_fold_details(
    const FoldExportContext& ctx,
    int outer_fold_idx,
    int nested_fold_idx,
    const SharedFoldData& nested_shared,
    const std::vector<FoldData>& nested_trait_folds);

/// Export nested_fold_ids.txt (call once before nested loop).
void export_nested_fold_ids(
    const FoldExportContext& ctx,
    int outer_fold_idx,
    const std::vector<int>& train_pool_idx,
    const std::vector<int>& nested_fold_ids);

/// Initialize train_val_residuals.txt with header (call before nested loop).
/// Returns the file path for subsequent appends.
std::string init_train_val_residuals(
    const FoldExportContext& ctx,
    int outer_fold_idx,
    int n_nested_folds);

/// Append one column of residuals to train_val_residuals.txt.
/// trait_suffix is appended to the column name (e.g. "_trait1" or "" for single-trait).
void append_train_val_residual_column(
    const std::string& path,
    const std::vector<std::string>& sample_ids,
    const std::vector<int>& train_pool_idx,
    int nested_fold_idx,
    const Eigen::VectorXd& residuals,
    const std::string& trait_suffix = "");

} // namespace fastcv
