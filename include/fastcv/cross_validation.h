#pragma once

#include "fastcv/types.h"
#include <vector>

namespace fastcv {

/// Create fold assignments (1 to n_folds), single stratification variable.
/// @param n            Number of samples
/// @param n_folds      Number of folds
/// @param stratify_by  Group labels for stratification
/// @param seed         Random seed
/// @return Vector of fold assignments (1-indexed)
std::vector<int> create_folds(int n, int n_folds,
                              const std::vector<int>& stratify_by,
                              int seed = 42);

/// Create fold assignments with multi-factor stratification.
///
/// Strategy: form composite cells from all factor combinations,
/// then round-robin assign within each cell. This ensures every
/// factor level is proportionally represented in every fold.
///
/// @param stratify_cols  Vector of factor columns (each is n-length int vector)
/// @param n_folds       Number of folds
/// @param seed          Random seed
/// @return Vector of fold assignments (1-indexed)
std::vector<int> create_folds_multi(const std::vector<std::vector<int>>& stratify_cols,
                                   int n_folds,
                                   int seed = 42);

/// Detect which columns in a covariate matrix are factors.
/// A factor column has <= max_categories unique values (integer-coded).
/// @param cov       Covariate matrix (n_samples x n_covars), int-coded
/// @param max_categories  Columns with <= this many unique values are factors
/// @return Vector of factor column indices
std::vector<int> detect_factor_cols(const Eigen::MatrixXi& cov, int max_categories = 10);

/// Split data into train and test based on fold assignment
struct TrainTestSplit {
    std::vector<int> train;
    std::vector<int> test;
};
TrainTestSplit train_test_split(const std::vector<int>& fold_ids, int test_fold);

} // namespace fastcv
