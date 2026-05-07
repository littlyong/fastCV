#include "fastcv/data_export.h"
#include "fastcv/misc_util.h"
#include "json.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>

namespace fastcv {

using json = nlohmann::json;

// ============================================================
// export_cv_meta(): fold_assignments.csv + config.json
// ============================================================
void export_cv_meta(const FoldExportContext& ctx) {
    mkdir_p(ctx.output_dir);

    // --- fold_assignments.txt ---
    {
        std::string path = ctx.output_dir + "/fold_assignments.txt";
        std::ofstream out(path);
        if (!out.is_open()) throw std::runtime_error("Cannot write: " + path);
        out << "sample_id fold\n";
        for (size_t i = 0; i < ctx.sample_ids.size(); ++i)
            out << ctx.sample_ids[i] << " " << ctx.fold_ids[i] << "\n";
        std::cout << "[INFO] Saved: " << path << "\n";
    }

    // --- datainfo.json ---
    {
        std::string path = ctx.output_dir + "/datainfo.json";
        json cfg;
        cfg["n_geno_pcs"] = ctx.n_geno_pcs;
        cfg["n_pheno_pcs"] = ctx.n_pheno_pcs;
        cfg["correction_mode"] = correction_mode_to_string(ctx.correction_mode);
        cfg["n_traits"] = static_cast<int>(ctx.trait_names.size());
        json traits_json = json::array();
        for (const auto& t : ctx.trait_names) traits_json.push_back(t);
        cfg["traits"] = traits_json;

        std::ofstream out(path);
        out << cfg.dump(2) << "\n";
        std::cout << "[INFO] Saved: " << path << "\n";
    }
}

// ============================================================
// export_single_fold(): export one fold's data and free memory
// ============================================================
void export_single_fold(const FoldExportContext& ctx,
                        int fold_idx,
                        const SharedFoldData& shared,
                        const std::vector<FoldData>& trait_folds,
                        bool save_details) {
    std::string fold_dir = make_fold_dir(ctx.output_dir, fold_idx);
    mkdir_p(fold_dir);

    int n_train = static_cast<int>(shared.train_idx.size());
    int n_test = static_cast<int>(shared.test_idx.size());
    int n_traits = static_cast<int>(trait_folds.size());

    // --- y_train_residual.txt (always exported) ---
    {
        std::string path = fold_dir + "/y_train_residual.txt";
        std::ofstream out(path);
        out << "sample_id";
        for (const auto& tn : ctx.trait_names) out << " " << tn;
        out << "\n";
        for (int i = 0; i < n_train; ++i) {
            out << ctx.sample_ids[shared.train_idx[i]];
            for (int t = 0; t < n_traits; ++t)
                out << " " << trait_folds[t].y_train_residual(i);
            out << "\n";
        }
    }

    // --- y_test_residual.txt (always exported) ---
    {
        std::string path = fold_dir + "/y_test_residual.txt";
        std::ofstream out(path);
        out << "sample_id";
        for (const auto& tn : ctx.trait_names) out << " " << tn;
        out << "\n";
        for (int i = 0; i < n_test; ++i) {
            out << ctx.sample_ids[shared.test_idx[i]];
            for (int t = 0; t < n_traits; ++t)
                out << " " << trait_folds[t].y_test_residual(i);
            out << "\n";
        }
    }

    // --- Intermediate detail files (gated by save_details) ---
    if (save_details) {
        // --- train_geno_pcs.txt ---
        {
            std::string path = fold_dir + "/train_geno_pcs.txt";
            std::ofstream out(path);
            out << "sample_id";
            for (int k = 0; k < ctx.n_geno_pcs; ++k)
                out << " geno_pc" << (k + 1);
            out << "\n";
            for (int i = 0; i < n_train; ++i) {
                out << ctx.sample_ids[shared.train_idx[i]];
                for (int k = 0; k < ctx.n_geno_pcs; ++k)
                    out << " " << shared.geno_pcs_train(i, k);
                out << "\n";
            }
        }

        // --- geno_eigenvalues.txt ---
        {
            std::string path = fold_dir + "/geno_eigenvalues.txt";
            std::ofstream out(path);
            const auto& ev = shared.geno_dim.eigenvalues;
            out << "pc eigenvalue\n";
            for (int k = 0; k < ev.size(); ++k)
                out << "pc" << (k + 1) << " " << ev(k) << "\n";
        }

        // --- test_geno_pcs.txt ---
        {
            std::string path = fold_dir + "/test_geno_pcs.txt";
            std::ofstream out(path);
            out << "sample_id";
            for (int k = 0; k < ctx.n_geno_pcs; ++k)
                out << " geno_pc" << (k + 1);
            out << "\n";
            for (int i = 0; i < n_test; ++i) {
                out << ctx.sample_ids[shared.test_idx[i]];
                for (int k = 0; k < ctx.n_geno_pcs; ++k)
                    out << " " << shared.geno_pcs_test(i, k);
                out << "\n";
            }
        }

        // --- train_pheno_pcs.txt (only in geno_pheno mode) ---
        if (ctx.correction_mode == CorrectionMode::GENO_PHENO &&
            shared.pheno_pcs_train.cols() > 0) {
            // Train
            {
                std::string path = fold_dir + "/train_pheno_pcs.txt";
                std::ofstream out(path);
                int n_ph = static_cast<int>(shared.pheno_pcs_train.cols());
                out << "sample_id";
                for (int k = 0; k < n_ph; ++k)
                    out << " pheno_pc" << (k + 1);
                out << "\n";
                for (int i = 0; i < n_train; ++i) {
                    out << ctx.sample_ids[shared.train_idx[i]];
                    for (int k = 0; k < n_ph; ++k)
                        out << " " << shared.pheno_pcs_train(i, k);
                    out << "\n";
                }
            }
            // Test
            {
                std::string path = fold_dir + "/test_pheno_pcs.txt";
                std::ofstream out(path);
                int n_ph = static_cast<int>(shared.pheno_pcs_test.cols());
                out << "sample_id";
                for (int k = 0; k < n_ph; ++k)
                    out << " pheno_pc" << (k + 1);
                out << "\n";
                for (int i = 0; i < n_test; ++i) {
                    out << ctx.sample_ids[shared.test_idx[i]];
                    for (int k = 0; k < n_ph; ++k)
                        out << " " << shared.pheno_pcs_test(i, k);
                    out << "\n";
                }
            }
        }

        // --- train_beta.txt ---
        {
            std::string path = fold_dir + "/train_beta.txt";
            std::ofstream out(path);
            out << "coefficient";
            for (const auto& tn : ctx.trait_names) out << " " << tn;
            out << "\n";

            int n_beta = static_cast<int>(trait_folds[0].beta.size());
            std::vector<std::string> row_names;
            row_names.push_back("intercept");
            for (const auto& name : ctx.cov_col_names) row_names.push_back(name);
            for (int i = 0; i < ctx.n_geno_pcs; ++i)
                row_names.push_back("geno_pc" + std::to_string(i + 1));
            if (ctx.correction_mode == CorrectionMode::GENO_PHENO)
                for (int i = 0; i < ctx.n_pheno_pcs; ++i)
                    row_names.push_back("pheno_pc" + std::to_string(i + 1));

            for (int r = 0; r < n_beta; ++r) {
                out << row_names[r];
                for (int t = 0; t < n_traits; ++t)
                    out << " " << trait_folds[t].beta(r);
                out << "\n";
            }
        }

        // --- train_maf.txt ---
        {
            std::string path = fold_dir + "/train_maf.txt";
            std::ofstream out(path);
            out << "snp_id maf\n";
            for (int j = 0; j < shared.geno_dim.maf.size(); ++j)
                out << ctx.snp_ids[j] << " " << shared.geno_dim.maf(j) << "\n";
        }
    }

    std::cout << "[INFO] Fold " << fold_idx << " saved to: " << fold_dir << "\n";
}

// ============================================================
// export_test_residuals(): outer test residuals for nested mode
// ============================================================
void export_test_residuals(const FoldExportContext& ctx,
                            int fold_idx,
                            const std::vector<int>& test_idx,
                            const std::vector<FoldData>& trait_folds) {
    std::string fold_dir = make_fold_dir(ctx.output_dir, fold_idx);
    mkdir_p(fold_dir);

    int n_test = static_cast<int>(test_idx.size());
    int n_traits = static_cast<int>(trait_folds.size());

    std::string path = fold_dir + "/y_test_residual.txt";
    std::ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("Cannot write: " + path);
    out << "sample_id";
    for (const auto& tn : ctx.trait_names) out << " " << tn;
    out << "\n";
    for (int i = 0; i < n_test; ++i) {
        out << ctx.sample_ids[test_idx[i]];
        for (int t = 0; t < n_traits; ++t)
            out << " " << trait_folds[t].y_test_residual(i);
        out << "\n";
    }
    std::cout << "[INFO] Saved: " << path << "\n";
}

// ============================================================
// export_nested_cv_results()
// ============================================================
void export_nested_cv_results(
    const FoldExportContext& ctx,
    int outer_fold_idx,
    const std::vector<int>& train_pool_idx,
    const std::vector<int>& nested_fold_ids,
    const Eigen::MatrixXd& residual_matrix,
    bool save_details,
    const Eigen::VectorXd& outer_eigenvalues,
    const Eigen::VectorXd& outer_maf,
    const std::vector<SharedFoldData>& nested_shared,
    const std::vector<std::vector<FoldData>>& nested_trait_folds)
{
    std::string fold_dir = make_fold_dir(ctx.output_dir, outer_fold_idx);
    mkdir_p(fold_dir);

    int n_pool = static_cast<int>(train_pool_idx.size());
    int n_nested = static_cast<int>(residual_matrix.cols());
    int n_traits = static_cast<int>(ctx.trait_names.size());

    // --- nested_fold_ids.txt ---
    {
        std::string path = fold_dir + "/nested_fold_ids.txt";
        std::ofstream out(path);
        if (!out.is_open()) throw std::runtime_error("Cannot write: " + path);
        out << "sample_id nested_fold\n";
        for (int i = 0; i < n_pool; ++i)
            out << ctx.sample_ids[train_pool_idx[i]] << " " << nested_fold_ids[i] << "\n";
        std::cout << "[INFO] Saved: " << path << "\n";
    }

    // --- train_val_residuals.txt ---
    {
        std::string path = fold_dir + "/train_val_residuals.txt";
        std::ofstream out(path);
        if (!out.is_open()) throw std::runtime_error("Cannot write: " + path);
        out << "sample_id";
        for (int j = 0; j < n_nested; ++j)
            out << " nested_" << (j + 1);
        out << "\n";
        for (int i = 0; i < n_pool; ++i) {
            out << ctx.sample_ids[train_pool_idx[i]];
            for (int j = 0; j < n_nested; ++j)
                out << " " << residual_matrix(i, j);
            out << "\n";
        }
        std::cout << "[INFO] Saved: " << path << "\n";
    }

    // --- Optional: nested details ---
    if (save_details && !nested_shared.empty()) {
        for (size_t nf = 0; nf < nested_shared.size(); ++nf) {
            export_nested_fold_details(ctx, outer_fold_idx, nf + 1,
                                       nested_shared[nf], nested_trait_folds[nf]);
        }
    }

    std::cout << "[INFO] Nested CV results saved to: " << fold_dir << "\n";
}

// ============================================================
// export_nested_fold_details(): one nested fold's detail files
// ============================================================
void export_nested_fold_details(
    const FoldExportContext& ctx,
    int outer_fold_idx,
    int nested_fold_idx,
    const SharedFoldData& nested_shared,
    const std::vector<FoldData>& nested_trait_folds)
{
    std::string fold_dir = make_fold_dir(ctx.output_dir, outer_fold_idx);
    std::string nested_dir = fold_dir + "/nested";
    mkdir_p(nested_dir);

    const auto& shared = nested_shared;
    const auto& trait_folds = nested_trait_folds;
    int n_nt = static_cast<int>(shared.train_idx.size());
    int n_nv = static_cast<int>(shared.test_idx.size());
    int n_traits = static_cast<int>(ctx.trait_names.size());

    std::string nd = nested_dir + "/nested_" +
                     (nested_fold_idx < 10 ? "0" : "") +
                     std::to_string(nested_fold_idx);
    mkdir_p(nd);

    // geno_eigenvalues.txt
    {
        std::string path = nd + "/geno_eigenvalues.txt";
        std::ofstream out(path);
        const auto& ev = shared.geno_dim.eigenvalues;
        out << "pc eigenvalue\n";
        for (int k = 0; k < ev.size(); ++k)
            out << "pc" << (k + 1) << " " << ev(k) << "\n";
    }

    // train_geno_pcs.txt
    {
        std::string path = nd + "/train_geno_pcs.txt";
        std::ofstream out(path);
        out << "sample_id";
        for (int k = 0; k < ctx.n_geno_pcs; ++k)
            out << " geno_pc" << (k + 1);
        out << "\n";
        for (int i = 0; i < n_nt; ++i) {
            out << ctx.sample_ids[shared.train_idx[i]];
            for (int k = 0; k < ctx.n_geno_pcs; ++k)
                out << " " << shared.geno_pcs_train(i, k);
            out << "\n";
        }
    }

    // val_geno_pcs.txt
    {
        std::string path = nd + "/val_geno_pcs.txt";
        std::ofstream out(path);
        out << "sample_id";
        for (int k = 0; k < ctx.n_geno_pcs; ++k)
            out << " geno_pc" << (k + 1);
        out << "\n";
        for (int i = 0; i < n_nv; ++i) {
            out << ctx.sample_ids[shared.test_idx[i]];
            for (int k = 0; k < ctx.n_geno_pcs; ++k)
                out << " " << shared.geno_pcs_test(i, k);
            out << "\n";
        }
    }

    // train_beta.txt
    {
        std::string path = nd + "/train_beta.txt";
        std::ofstream out(path);
        out << "coefficient";
        for (const auto& tn : ctx.trait_names) out << " " << tn;
        out << "\n";
        int n_beta = static_cast<int>(trait_folds[0].beta.size());
        std::vector<std::string> row_names;
        row_names.push_back("intercept");
        for (const auto& name : ctx.cov_col_names) row_names.push_back(name);
        for (int i = 0; i < ctx.n_geno_pcs; ++i)
            row_names.push_back("geno_pc" + std::to_string(i + 1));
        if (ctx.correction_mode == CorrectionMode::GENO_PHENO)
            for (int i = 0; i < ctx.n_pheno_pcs; ++i)
                row_names.push_back("pheno_pc" + std::to_string(i + 1));
        for (int r = 0; r < n_beta; ++r) {
            out << row_names[r];
            for (int t = 0; t < n_traits; ++t)
                out << " " << trait_folds[t].beta(r);
            out << "\n";
        }
    }

    // train_maf.txt
    {
        std::string path = nd + "/train_maf.txt";
        std::ofstream out(path);
        out << "snp_id maf\n";
        for (int j = 0; j < shared.geno_dim.maf.size(); ++j)
            out << ctx.snp_ids[j] << " " << shared.geno_dim.maf(j) << "\n";
    }

    std::cout << "[INFO] Saved: " << nd << "\n";
}

// ============================================================
// export_nested_fold_ids(): nested fold assignments
// ============================================================
void export_nested_fold_ids(
    const FoldExportContext& ctx,
    int outer_fold_idx,
    const std::vector<int>& train_pool_idx,
    const std::vector<int>& nested_fold_ids)
{
    std::string fold_dir = make_fold_dir(ctx.output_dir, outer_fold_idx);
    mkdir_p(fold_dir);

    std::string path = fold_dir + "/nested_fold_ids.txt";
    std::ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("Cannot write: " + path);
    out << "sample_id nested_fold\n";
    for (size_t i = 0; i < train_pool_idx.size(); ++i)
        out << ctx.sample_ids[train_pool_idx[i]] << " " << nested_fold_ids[i] << "\n";
    std::cout << "[INFO] Saved: " << path << "\n";
}

// ============================================================
// init_train_val_residuals(): create file with header only
// ============================================================
std::string init_train_val_residuals(
    const FoldExportContext& ctx,
    int outer_fold_idx,
    int n_nested_folds)
{
    std::string fold_dir = make_fold_dir(ctx.output_dir, outer_fold_idx);
    mkdir_p(fold_dir);

    std::string path = fold_dir + "/train_val_residuals.txt";
    std::ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("Cannot write: " + path);
    out << "sample_id";
    out << "\n";
    return path;
}

// ============================================================
// append_train_val_residual_column(): rewrite file with new column
// ============================================================
void append_train_val_residual_column(
    const std::string& path,
    const std::vector<std::string>& sample_ids,
    const std::vector<int>& train_pool_idx,
    int nested_fold_idx,
    const Eigen::VectorXd& residuals,
    const std::string& trait_suffix)
{
    // Read existing content (header + previous columns)
    std::vector<std::string> lines;
    {
        std::ifstream in(path);
        if (!in.is_open()) throw std::runtime_error("Cannot read: " + path);
        std::string line;
        while (std::getline(in, line)) lines.push_back(line);
    }

    // Append new column to header
    lines[0] += " nested_" + std::to_string(nested_fold_idx) + trait_suffix;

    // Append residual values to each data line (or add new lines if first column)
    int n_pool = static_cast<int>(train_pool_idx.size());
    int existing_cols = static_cast<int>(lines.size()) - 1;  // minus header

    if (existing_cols == 0) {
        // First column: write data lines
        for (int i = 0; i < n_pool; ++i)
            lines.push_back(sample_ids[train_pool_idx[i]] + " " +
                            std::to_string(residuals(i)));
    } else {
        // Append to existing data lines
        for (int i = 0; i < n_pool && i < existing_cols; ++i)
            lines[i + 1] += " " + std::to_string(residuals(i));
    }

    // Rewrite file
    std::ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("Cannot write: " + path);
    for (const auto& line : lines) out << line << "\n";
}

} // namespace fastcv
