#include "fastcv/cross_validation.h"

#include <random>
#include <map>
#include <set>
#include <iostream>
#include <stdexcept>

namespace fastcv {

// ============================================================
// create_folds() — single stratification
// ============================================================
std::vector<int> create_folds(int n, int n_folds,
                              const std::vector<int>& stratify_by,
                              int seed) {
    if (n < n_folds) {
        throw std::runtime_error("Number of samples (" + std::to_string(n) +
                                 ") less than n_folds (" + std::to_string(n_folds) + ")");
    }

    std::vector<int> folds(n);
    std::mt19937 rng(seed);

    if (stratify_by.empty()) {
        // No stratification: random shuffle + round-robin
        std::vector<int> indices(n);
        for (int i = 0; i < n; ++i) indices[i] = i;
        std::shuffle(indices.begin(), indices.end(), rng);
        for (int i = 0; i < n; ++i) {
            folds[indices[i]] = (i % n_folds) + 1;
        }
    } else {
        if (static_cast<int>(stratify_by.size()) != n) {
            throw std::runtime_error("stratify_by size (" +
                                     std::to_string(stratify_by.size()) +
                                     ") != n (" + std::to_string(n) + ")");
        }
        // Collect indices per group, shuffle, round-robin assign
        std::map<int, std::vector<int>> groups;
        for (int i = 0; i < n; ++i) {
            groups[stratify_by[i]].push_back(i);
        }
        for (auto& [val, indices] : groups) {
            std::shuffle(indices.begin(), indices.end(), rng);
            for (size_t i = 0; i < indices.size(); ++i) {
                folds[indices[i]] = (static_cast<int>(i) % n_folds) + 1;
            }
        }
    }

    return folds;
}

// ============================================================
// detect_factor_cols() — auto-detect factor-type columns
// ============================================================
std::vector<int> detect_factor_cols(const Eigen::MatrixXi& cov, int max_categories) {
    int n_cols = cov.cols();
    std::vector<int> factor_cols;

    for (int c = 0; c < n_cols; ++c) {
        std::set<int> unique_vals;
        for (int i = 0; i < cov.rows(); ++i) {
            unique_vals.insert(cov(i, c));
            if (static_cast<int>(unique_vals.size()) > max_categories) {
                break;  // not a factor
            }
        }
        if (static_cast<int>(unique_vals.size()) <= max_categories && unique_vals.size() > 1) {
            factor_cols.push_back(c);
        }
    }

    return factor_cols;
}

// ============================================================
// create_folds_multi() — multi-factor stratified folds
//
// Algorithm:
//   1. For each sample, compute a composite key from all factor columns.
//      key = (f1[0], f2[0], ...), (f1[1], f2[1], ...), ...
//   2. Group samples by composite key (each group is a "cell").
//   3. Within each cell, shuffle and round-robin assign to folds.
//
// This guarantees proportional representation of every factor level
// in every fold, as long as n_folds <= min_cell_size.
// ============================================================
std::vector<int> create_folds_multi(const std::vector<std::vector<int>>& stratify_cols,
                                   int n_folds,
                                   int seed) {
    if (stratify_cols.empty()) {
        throw std::runtime_error("No stratification columns provided");
    }

    int n = stratify_cols[0].size();
    for (const auto& col : stratify_cols) {
        if (static_cast<int>(col.size()) != n) {
            throw std::runtime_error("All stratification columns must have size n");
        }
    }

    std::vector<int> folds(n);

    // 1. Build composite key for each sample
    //    key = f0 * stride1 + f1 * stride2 + ...
    std::vector<long long> composite_key(n, 0);
    {
        long long stride = 1;
        for (int c = 0; c < static_cast<int>(stratify_cols.size()); ++c) {
            for (int i = 0; i < n; ++i) {
                composite_key[i] += static_cast<long long>(stratify_cols[c][i]) * stride;
            }
            // Shift to avoid collisions with large values
            stride *= 1000;
        }
    }

    // 2. Group samples by composite key
    std::map<long long, std::vector<int>> cells;
    for (int i = 0; i < n; ++i) {
        cells[composite_key[i]].push_back(i);
    }

    // Log info
    int n_cells = static_cast<int>(cells.size());
    int min_cell = n, max_cell = 0;
    for (const auto& [key, indices] : cells) {
        min_cell = std::min(min_cell, static_cast<int>(indices.size()));
        max_cell = std::max(max_cell, static_cast<int>(indices.size()));
    }
    std::cout << "[INFO] Multi-factor stratification: " << stratify_cols.size()
              << " factors, " << n_cells << " cells (min_cell=" << min_cell
              << ", max_cell=" << max_cell << ", n_folds=" << n_folds << ")\n";
    if (min_cell < n_folds) {
        std::cout << "[WARN] Some cells have fewer samples than folds ("
                  << min_cell << " < " << n_folds << "). "
                  << "Stratification may not be perfectly balanced.\n";
    }

    // 3. Within each cell, shuffle and round-robin assign
    std::mt19937 rng(seed);
    for (auto& [key, indices] : cells) {
        std::shuffle(indices.begin(), indices.end(), rng);
        for (size_t i = 0; i < indices.size(); ++i) {
            folds[indices[i]] = (static_cast<int>(i) % n_folds) + 1;
        }
    }

    return folds;
}

// ============================================================
// train_test_split()
// ============================================================
TrainTestSplit train_test_split(const std::vector<int>& fold_ids, int test_fold) {
    TrainTestSplit split;
    for (size_t i = 0; i < fold_ids.size(); ++i) {
        if (fold_ids[i] == test_fold) {
            split.test.push_back(static_cast<int>(i));
        } else {
            split.train.push_back(static_cast<int>(i));
        }
    }
    return split;
}

} // namespace fastcv
