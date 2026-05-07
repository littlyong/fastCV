#include "fastcv/pc_corrector.h"
#include "fastcv/dim_reduction.h"
#include "fastcv/cross_validation.h"
#include "fastcv/data_export.h"
#include "fastcv/external_pca.h"
#include "fastcv/grm.h"
#include "fastcv/ld_prune.h"
#include "fastcv/misc_util.h"

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <numeric>
#include <unordered_map>
#include <memory>

#ifdef FASTCV_HAS_OPENMP
#include <omp.h>
#endif

#ifdef FASTCV_USE_MKL
#include <mkl.h>
#endif

namespace fastcv {

// ============================================================
// TeeStream: writes to both std::cout and a log file
// ============================================================
class TeeStream : public std::streambuf {
public:
    TeeStream(std::streambuf* primary, const std::string& log_path)
        : primary_(primary) {
        log_file_.open(log_path);
        if (!log_file_.is_open()) {
            std::cerr << "[WARN] Cannot open log file: " << log_path << "\n";
        }
    }
    ~TeeStream() {
        std::ostream(primary_).flush();
        log_file_.close();
    }

    std::streambuf* original_buf() const { return primary_; }

protected:
    int overflow(int c) override {
        if (c != EOF) {
            primary_->sputc(c);
            if (log_file_.is_open()) log_file_.put(c);
        }
        return c;
    }

    int sync() override {
        primary_->pubsync();
        if (log_file_.is_open()) log_file_.flush();
        return 0;
    }

private:
    std::streambuf* primary_;
    std::ofstream log_file_;
};

// ============================================================
// Timing helper: returns elapsed seconds and pretty-prints
// ============================================================
struct Timer {
    std::chrono::high_resolution_clock::time_point start;
    const char* label;
    Timer(const char* l) : start(std::chrono::high_resolution_clock::now()), label(l) {}
    double elapsed() const {
        return std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - start).count();
    }
    ~Timer() {
        double sec = elapsed();
        std::cout << "    [" << label << "] " << sec << "s\n" << std::flush;
    }
};

// ============================================================
// Internal helper: compute OLS hat matrix (X'X)^{-1} X'
// LLT with SVD fallback for near-singular design matrices
// ============================================================
static Eigen::MatrixXd compute_ols_hat_matrix(const Eigen::MatrixXd& X) {
    Eigen::MatrixXd XtX = X.transpose() * X;
    Eigen::MatrixXd XtX_inv;
    Eigen::LLT<Eigen::MatrixXd> llt(XtX);
    if (llt.info() == Eigen::Success) {
        XtX_inv = llt.solve(Eigen::MatrixXd::Identity(XtX.rows(), XtX.cols()));
    } else {
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(XtX, Eigen::ComputeThinU | Eigen::ComputeThinV);
        XtX_inv = svd.solve(Eigen::MatrixXd::Identity(XtX.rows(), XtX.cols()));
    }
    return XtX_inv * X.transpose();
}

// ============================================================
// Internal helper: build design matrix X = [1, cov, geno_pcs, pheno_pcs]
// ============================================================
static Eigen::MatrixXd build_design_matrix(
    const Eigen::MatrixXd& geno_pcs,
    const Eigen::MatrixXd* pheno_pcs,
    const Eigen::MatrixXd* cov = nullptr)
{
    int n = geno_pcs.rows();
    int n_c = cov ? cov->cols() : 0;
    int n_g = geno_pcs.cols();
    int n_p = pheno_pcs ? pheno_pcs->cols() : 0;
    int p = 1 + n_c + n_g + n_p;

    Eigen::MatrixXd X(n, p);
    X.col(0).setOnes();  // intercept
    int offset = 1;
    if (n_c > 0) {
        X.block(0, offset, n, n_c) = *cov;
        offset += n_c;
    }
    X.block(0, offset, n, n_g) = geno_pcs;
    offset += n_g;
    if (n_p > 0) {
        X.block(0, offset, n, n_p) = *pheno_pcs;
    }
    return X;
}

// ============================================================
// prepare_cv_data() — main pipeline (streaming per-fold output)
// ============================================================
CvResult prepare_cv_data(const CvConfig& config) {
    auto pipeline_start = std::chrono::high_resolution_clock::now();

    // Open genotype data
    PlinkReader reader;
    reader.open(config.bed_file);

    // Load phenotype data (all traits)
    PhenoData pheno;
    {
        std::ifstream fin(config.pheno_file);
        if (!fin.is_open()) {
            throw std::runtime_error("Cannot open phenotype file: " + config.pheno_file);
        }
        std::string header_line;
        std::getline(fin, header_line);
        std::istringstream hss(header_line);
        std::vector<std::string> col_names;
        std::string col;
        while (hss >> col) col_names.push_back(col);

        int id_col_idx = -1;
        for (int i = 0; i < static_cast<int>(col_names.size()); ++i) {
            if (col_names[i] == config.pheno_id_col) { id_col_idx = i; break; }
        }
        if (id_col_idx < 0) {
            throw std::runtime_error("ID column '" + config.pheno_id_col +
                                     "' not found in phenotype file");
        }

        std::vector<int> trait_cols;
        for (int i = 0; i < static_cast<int>(col_names.size()); ++i) {
            if (i != id_col_idx) trait_cols.push_back(i);
        }
        pheno.trait_names.clear();
        for (int tc : trait_cols) pheno.trait_names.push_back(col_names[tc]);

        std::vector<std::string> pheno_ids;
        std::vector<std::vector<double>> pheno_values;
        std::string line;
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::vector<std::string> fields;
            std::string field;
            while (iss >> field) fields.push_back(field);
            if (static_cast<int>(fields.size()) <= id_col_idx) continue;
            pheno_ids.push_back(fields[id_col_idx]);
            std::vector<double> vals;
            for (int tc : trait_cols) vals.push_back(std::stod(fields[tc]));
            pheno_values.push_back(vals);
        }

        int n_pheno = static_cast<int>(pheno_ids.size());
        int n_traits = static_cast<int>(trait_cols.size());
        pheno.pheno_all.resize(n_pheno, n_traits);
        for (int i = 0; i < n_pheno; ++i)
            for (int j = 0; j < n_traits; ++j)
                pheno.pheno_all(i, j) = pheno_values[i][j];
        pheno.sample_ids = pheno_ids;
        pheno.is_multi_pheno = (n_traits > 1);
    }

    // Match samples between genotype and phenotype
    int n = reader.n_samples();
    std::vector<std::string> geno_ids = reader.sample_ids();
    std::unordered_map<std::string, int> pheno_id_map;
    for (int i = 0; i < static_cast<int>(pheno.sample_ids.size()); ++i)
        pheno_id_map[pheno.sample_ids[i]] = i;

    std::vector<int> matched_geno_idx, matched_pheno_idx;
    std::vector<std::string> matched_ids;
    for (int i = 0; i < n; ++i) {
        auto it = pheno_id_map.find(geno_ids[i]);
        if (it != pheno_id_map.end()) {
            matched_geno_idx.push_back(i);
            matched_pheno_idx.push_back(it->second);
            matched_ids.push_back(geno_ids[i]);
        }
    }
    if (matched_geno_idx.empty()) {
        throw std::runtime_error("No matching samples between genotype and phenotype files");
    }

    int n_geno_total = n;
    int n_pheno_total = static_cast<int>(pheno.sample_ids.size());

    // Build matched phenotype matrix (all traits)
    int n_matched = static_cast<int>(matched_geno_idx.size());
    int n_traits = static_cast<int>(pheno.trait_names.size());
    Eigen::MatrixXd pheno_all_matched(n_matched, n_traits);
    for (int i = 0; i < n_matched; ++i)
        pheno_all_matched.row(i) = pheno.pheno_all.row(matched_pheno_idx[i]);

    // Parse trait expression
    std::vector<int> selected_traits;
    if (!config.trait_name.empty()) {
        // Select trait by column name
        for (int i = 0; i < n_traits; ++i) {
            if (pheno.trait_names[i] == config.trait_name) {
                selected_traits.push_back(i);
                break;
            }
        }
        if (selected_traits.empty()) {
            throw std::runtime_error("Trait name '" + config.trait_name +
                                     "' not found in phenotype file");
        }
    } else {
        selected_traits = parse_trait_expr(config.trait_expr, n_traits);
    }
    if (selected_traits.empty()) {
        throw std::runtime_error("No valid traits selected: '" + config.trait_expr +
                                 "' (available: 1-" + std::to_string(n_traits) + ")");
    }

    // Set up log file (tee std::cout to output_dir/fastcv.log)
    // Must live until function returns
    std::string log_path;
    std::unique_ptr<TeeStream> log_tee;
    {
        // Ensure output directory exists
        mkdir_p(config.output_dir);

        log_path = config.output_dir + "/fastcv.log";
        log_tee = std::make_unique<TeeStream>(std::cout.rdbuf(), log_path);
        std::cout.rdbuf(log_tee.get());
    }

    // Print configuration
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "fastCV: Cross-Validation Data Preparation\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << "Traits (" << (config.trait_name.empty() ? config.trait_expr : config.trait_name) << "): ";
    for (size_t i = 0; i < selected_traits.size(); ++i) {
        if (i > 0 && i == selected_traits.size() - 1) std::cout << " and ";
        else if (i > 0) std::cout << ", ";
        std::cout << pheno.trait_names[selected_traits[i]];
    }
    std::cout << " (" << selected_traits.size() << " of " << n_traits << ")\n";
    std::cout << "Sample matching:\n";
    std::cout << "  Genotype:  " << n_geno_total << " samples\n";
    std::cout << "  Phenotype: " << n_pheno_total << " samples\n";
    std::cout << "  Geno-Pheno intersection: " << n_matched << " samples";
    if (n_geno_total > n_matched || n_pheno_total > n_matched)
        std::cout << " (" << (n_geno_total - n_matched) << " geno-only, "
                  << (n_pheno_total - n_matched) << " pheno-only excluded)";
    std::cout << "\n";
    std::cout << "  SNPs: " << reader.n_snps() << "\n";
    std::cout << "Folds: " << config.n_folds << " | Seed: " << config.seed << "\n";
    std::cout << "Correction: " << correction_mode_to_string(config.correction_mode) << "\n";
    std::cout << "  Genotype PCs: " << config.n_geno_pcs << "\n";
    if (config.correction_mode == CorrectionMode::GENO_PHENO)
        std::cout << "  Phenotype PCs: " << config.n_pheno_pcs << " (from all " << n_traits << " traits)\n";
    std::cout << "DimReduction: " << config.dim_reduction_method << "\n";
    if (config.n_nested_folds > 0) {
        int eff_nested_seed = (config.nested_seed > 0) ? config.nested_seed : config.seed;
        std::cout << "Nested CV: " << config.n_nested_folds
                  << " folds (seed=" << eff_nested_seed << ")\n";
    }
    std::cout << std::string(60, '=') << "\n\n";

    // Create folds (with multi-factor stratification if covariates provided)
    std::vector<int> fold_ids;
    Eigen::MatrixXd cov_matched;
    std::vector<std::string> cov_col_names_stored;
    if (!config.cov_file.empty()) {
        std::cout << "Loading covariates from: " << config.cov_file << "\n";
        std::ifstream cov_fin(config.cov_file);
        if (!cov_fin.is_open()) {
            throw std::runtime_error("Cannot open covariate file: " + config.cov_file);
        }
        std::string cov_header;
        std::getline(cov_fin, cov_header);
        std::istringstream cov_hss(cov_header);
        std::vector<std::string> cov_col_names;
        std::string cov_col;
        while (cov_hss >> cov_col) cov_col_names.push_back(cov_col);

        int cov_id_idx = -1;
        for (int i = 0; i < static_cast<int>(cov_col_names.size()); ++i) {
            if (cov_col_names[i] == config.cov_id_col) { cov_id_idx = i; break; }
        }
        if (cov_id_idx < 0) {
            throw std::runtime_error("ID column '" + config.cov_id_col +
                                     "' not found in covariate file");
        }

        // Determine which covariate columns to use
        std::vector<int> cov_use_idx;  // indices into cov_col_names (excluding ID col)
        if (config.cov_select == "all") {
            for (int i = 0; i < static_cast<int>(cov_col_names.size()); ++i)
                if (i != cov_id_idx) cov_use_idx.push_back(i);
        } else {
            // Parse comma-separated names or indices
            std::istringstream ss(config.cov_select);
            std::string part;
            while (std::getline(ss, part, ',')) {
                // Try as integer first
                try {
                    int idx = std::stoi(part);
                    if (idx >= 1 && idx <= static_cast<int>(cov_col_names.size()))
                        cov_use_idx.push_back(idx - 1);
                } catch (...) {
                    // Treat as column name
                    for (int i = 0; i < static_cast<int>(cov_col_names.size()); ++i) {
                        if (cov_col_names[i] == part) {
                            cov_use_idx.push_back(i);
                            break;
                        }
                    }
                }
            }
        }

        std::vector<std::vector<int>> cov_values;
        std::vector<std::string> cov_ids;
        std::string cov_line;
        while (std::getline(cov_fin, cov_line)) {
            if (cov_line.empty()) continue;
            std::istringstream iss(cov_line);
            std::vector<std::string> fields;
            std::string field;
            while (iss >> field) fields.push_back(field);
            if (static_cast<int>(fields.size()) <= cov_id_idx) continue;
            cov_ids.push_back(fields[cov_id_idx]);
            std::vector<int> vals;
            for (int ci : cov_use_idx) {
                try { vals.push_back(std::stoi(fields[ci])); }
                catch (...) { vals.push_back(0); }
            }
            cov_values.push_back(vals);
        }
        int n_cov_samples = static_cast<int>(cov_ids.size());
        int n_cov_cols = static_cast<int>(cov_use_idx.size());

        cov_col_names_stored.clear();
        for (int ci : cov_use_idx) cov_col_names_stored.push_back(cov_col_names[ci]);

        std::cout << "  Covariate samples: " << n_cov_samples
                  << ", columns: " << n_cov_cols << "\n";

        std::unordered_map<std::string, int> cov_id_map;
        for (int i = 0; i < n_cov_samples; ++i) cov_id_map[cov_ids[i]] = i;

        // Three-way intersection: genotype ∩ phenotype ∩ covariate
        std::vector<int> new_geno_idx, new_pheno_idx;
        std::vector<std::string> new_ids, dropped_from_cov;
        for (int i = 0; i < n_matched; ++i) {
            auto it = cov_id_map.find(matched_ids[i]);
            if (it != cov_id_map.end()) {
                new_geno_idx.push_back(matched_geno_idx[i]);
                new_pheno_idx.push_back(matched_pheno_idx[i]);
                new_ids.push_back(matched_ids[i]);
            } else {
                dropped_from_cov.push_back(matched_ids[i]);
            }
        }
        int n_three_way = static_cast<int>(new_ids.size());
        std::cout << "  Covariate:  " << n_cov_samples << " samples\n";
        std::cout << "  Three-way intersection: " << n_three_way << " samples\n";

        if (!dropped_from_cov.empty()) {
            std::cout << "  WARNING: " << dropped_from_cov.size()
                      << " samples excluded (in geno+pheno but not covariate):\n";
            int show = std::min(static_cast<int>(dropped_from_cov.size()), 10);
            for (int i = 0; i < show; ++i)
                std::cout << "    " << dropped_from_cov[i] << "\n";
            if (static_cast<int>(dropped_from_cov.size()) > 10)
                std::cout << "    ... and " << (dropped_from_cov.size() - 10) << " more\n";

            matched_geno_idx = std::move(new_geno_idx);
            matched_pheno_idx = std::move(new_pheno_idx);
            matched_ids = std::move(new_ids);
            n_matched = n_three_way;

            // Rebuild phenotype matrix with filtered samples
            pheno_all_matched.resize(n_matched, n_traits);
            for (int i = 0; i < n_matched; ++i)
                pheno_all_matched.row(i) = pheno.pheno_all.row(matched_pheno_idx[i]);
        }

        // Build covariate matrix (all samples guaranteed to have covariates)
        cov_matched.setZero(n_matched, n_cov_cols);
        for (int i = 0; i < n_matched; ++i) {
            auto it = cov_id_map.find(matched_ids[i]);
            for (int j = 0; j < n_cov_cols; ++j)
                cov_matched(i, j) = static_cast<double>(cov_values[it->second][j]);
        }

        // Stratification: use specified column or auto-detect factor columns
        if (!config.stratified_cov_name.empty()) {
            // Find the stratification column by name
            int strat_col = -1;
            for (int i = 0; i < n_cov_cols; ++i) {
                if (cov_col_names_stored[i] == config.stratified_cov_name) {
                    strat_col = i;
                    break;
                }
            }
            if (strat_col < 0) {
                throw std::runtime_error("Stratified covariate column '" +
                                         config.stratified_cov_name + "' not found in covariate file");
            }
            std::vector<int> strat_vals(n_matched);
            Eigen::MatrixXi cov_matched_int = cov_matched.cast<int>();
            for (int i = 0; i < n_matched; ++i)
                strat_vals[i] = cov_matched_int(i, strat_col);
            fold_ids = create_folds(n_matched, config.n_folds, strat_vals, config.seed);
            std::cout << "  Stratifying by column: " << config.stratified_cov_name << "\n";
        } else {
            // Auto-detect factor columns for stratification
            Eigen::MatrixXi cov_matched_int = cov_matched.cast<int>();
            std::vector<int> factor_cols = detect_factor_cols(cov_matched_int, 10);
            if (factor_cols.empty()) {
                std::cout << "  No factor columns detected, using random folds.\n";
                std::vector<int> empty_strat;
                fold_ids = create_folds(n_matched, config.n_folds, empty_strat, config.seed);
            } else {
                std::cout << "  Detected " << factor_cols.size() << " factor column(s):";
                for (int fc : factor_cols) {
                    std::string col_name = (fc < n_cov_cols)
                        ? cov_col_names_stored[fc]
                        : "cov_" + std::to_string(fc);
                    std::cout << " " << col_name;
                }
                std::cout << "\n";
                std::vector<std::vector<int>> stratify_cols(factor_cols.size());
                for (size_t i = 0; i < factor_cols.size(); ++i) {
                    stratify_cols[i].resize(n_matched);
                    for (int j = 0; j < n_matched; ++j)
                        stratify_cols[i][j] = cov_matched_int(j, factor_cols[i]);
                }
                fold_ids = create_folds_multi(stratify_cols, config.n_folds, config.seed);
            }
        }
    } else {
        std::cout << "No covariate file specified, using random folds.\n";
        std::vector<int> empty_strat;
        fold_ids = create_folds(n_matched, config.n_folds, empty_strat, config.seed);
    }

    // Create dimensionality reduction method
    auto dim_method = create_dim_reduction(config.dim_reduction_method);

    // Set GRM block size if user specified (otherwise auto-detect on first fit)
    if (config.grm_block_size > 0) {
        auto* grm_pca = dynamic_cast<GrmPcaDimReduction*>(dim_method.get());
        if (grm_pca) grm_pca->set_block_size(config.grm_block_size, n_matched);
    }

    // Optional LD pruning before GRM computation
    std::vector<int> pruned_snp_idx;
    const std::vector<int>* snp_idx_ptr = nullptr;
    if (config.ld_prune) {
        LdPruneParams ld_params;
        ld_params.r2_threshold = config.ld_r2_threshold;
        ld_params.window_size = config.ld_window_size;
        ld_params.maf_min = config.ld_maf_min;
        std::cout << "  [LD prune] Pruning SNPs (r² < " << ld_params.r2_threshold
                  << ", window=" << ld_params.window_size
                  << ", MAF > " << ld_params.maf_min << ")...\n" << std::flush;
        pruned_snp_idx = ld_prune_snps(reader, matched_geno_idx, ld_params);
        snp_idx_ptr = &pruned_snp_idx;
        auto* grm_pca = dynamic_cast<GrmPcaDimReduction*>(dim_method.get());
        if (grm_pca) grm_pca->set_snp_idx(snp_idx_ptr);
    }

    // Collect SNP IDs for export
    std::vector<std::string> snp_ids;
    snp_ids.reserve(reader.n_snps());
    for (const auto& si : reader.snp_info()) snp_ids.push_back(si.snp_id);

    // Build trait names for export
    std::vector<std::string> trait_names_export;
    for (int t : selected_traits) trait_names_export.push_back(pheno.trait_names[t]);

    // Set up export context (lightweight, reused for every fold)
    FoldExportContext export_ctx;
    export_ctx.output_dir = config.output_dir;
    export_ctx.sample_ids = matched_ids;
    export_ctx.fold_ids = fold_ids;
    export_ctx.snp_ids = snp_ids;
    export_ctx.cov_col_names = cov_col_names_stored;
    export_ctx.trait_names = trait_names_export;
    export_ctx.n_geno_pcs = config.n_geno_pcs;
    export_ctx.n_pheno_pcs = config.n_pheno_pcs;
    export_ctx.correction_mode = config.correction_mode;
    export_ctx.separator = config.separator;

    // Export global metadata (fold_assignments.csv + config.json)
    {
        Timer t("META");
        export_cv_meta(export_ctx);
    }

    // Build minimal CvResult for return (only metadata, no fold data)
    CvResult result;
    result.cov_col_names = cov_col_names_stored;
    result.config.dim_reduction_method = config.dim_reduction_method;
    result.config.n_samples = n_matched;
    result.config.n_snps = reader.n_snps();
    result.config.n_folds = config.n_folds;
    result.config.n_geno_pcs = config.n_geno_pcs;
    result.config.n_pheno_pcs = config.n_pheno_pcs;
    result.config.correction_mode = config.correction_mode;
    result.config.seed = config.seed;
    result.config.stratify_col = config.stratified_cov_name;
    {
        auto now = std::chrono::system_clock::now();
        auto time = std::chrono::system_clock::to_time_t(now);
        std::ostringstream oss;
        oss << std::put_time(std::localtime(&time), "%Y%m%d_%H%M%S");
        result.config.timestamp = oss.str();
    }
    result.fold_ids = fold_ids;
    result.sample_ids = matched_ids;
    result.snp_ids = snp_ids;

    // ========================================================
    // Main loop: for each fold, compute, export, release memory
    // ========================================================
    double total_fold_time = 0.0;

    for (int f = 0; f < config.n_folds; ++f) {
        auto fold_start = std::chrono::high_resolution_clock::now();
        std::cout << "\n--- Fold " << (f + 1) << "/" << config.n_folds << " ---\n";

        auto split = train_test_split(fold_ids, f + 1);
        std::vector<int> train_idx, test_idx;
        for (int idx : split.train) train_idx.push_back(matched_geno_idx[idx]);
        for (int idx : split.test) test_idx.push_back(matched_geno_idx[idx]);

        std::vector<int> train_pheno_idx, test_pheno_idx;
        for (int idx : split.train) train_pheno_idx.push_back(idx);
        for (int idx : split.test) test_pheno_idx.push_back(idx);

        std::cout << "  Train: " << train_idx.size() << " | Test: " << test_idx.size() << "\n";

        // Prepare covariate sub-matrices
        Eigen::MatrixXd cov_train, cov_test;
        const Eigen::MatrixXd* cov_train_ptr = nullptr;
        const Eigen::MatrixXd* cov_test_ptr = nullptr;
        if (cov_matched.cols() > 0) {
            Eigen::VectorXi ti(train_pheno_idx.size()), te(test_pheno_idx.size());
            for (size_t i = 0; i < train_pheno_idx.size(); ++i) ti(i) = train_pheno_idx[i];
            for (size_t i = 0; i < test_pheno_idx.size(); ++i) te(i) = test_pheno_idx[i];
            cov_train = cov_matched(ti, Eigen::all).cast<double>();
            cov_test = cov_matched(te, Eigen::all).cast<double>();
            cov_train_ptr = &cov_train;
            cov_test_ptr = &cov_test;
        }

        // Prepare phenotype sub-matrices (all traits, for PCA)
        Eigen::MatrixXd pheno_train, pheno_test;
        const Eigen::MatrixXd* pheno_train_ptr = nullptr;
        const Eigen::MatrixXd* pheno_test_ptr = nullptr;
        if (config.correction_mode == CorrectionMode::GENO_PHENO && n_traits > 1) {
            Eigen::VectorXi ti(train_pheno_idx.size()), te(test_pheno_idx.size());
            for (size_t i = 0; i < train_pheno_idx.size(); ++i) ti(i) = train_pheno_idx[i];
            for (size_t i = 0; i < test_pheno_idx.size(); ++i) te(i) = test_pheno_idx[i];
            pheno_train = pheno_all_matched(ti, Eigen::all);
            pheno_test = pheno_all_matched(te, Eigen::all);
            pheno_train_ptr = &pheno_train;
            pheno_test_ptr = &pheno_test;
        }

        // --- Shared computation: genotype PCA (once per fold) ---
        std::cout << "  Computing genotype PCA...\n";
        SharedFoldData shared;
        shared.fold_idx = f + 1;
        shared.train_idx = split.train;
        shared.test_idx = split.test;

        bool use_external = (config.pca_backend != "internal");
        Eigen::VectorXd ext_singular_values;  // only used for external path

        if (use_external) {
            // External tool path: hiblup/plink computes GRM + eigendecomposition
            std::vector<std::string> train_fids;
            for (int idx : split.train) train_fids.push_back(matched_ids[idx]);

            // Get PLINK prefix (bed_file without .bed extension)
            std::string bed_prefix = config.bed_file;
            auto dot = bed_prefix.rfind('.');
            if (dot != std::string::npos) bed_prefix = bed_prefix.substr(0, dot);

            ExternalPcaResult ext_result = run_external_pca(
                config.pca_backend, config.pca_tool_path,
                bed_prefix, train_fids,
                config.n_geno_pcs, config.n_threads,
                config.output_dir, f + 1);

            shared.geno_pcs_train = ext_result.train_pcs;
            ext_singular_values = ext_result.singular_values;

            // Test projection via G_cross + external PCA results
            {
                Timer t("GENO proj test (external)");
                Eigen::MatrixXd G_cross = compute_grm_cross(reader, test_idx, train_idx,
                                                            5000, config.n_threads);
                shared.geno_pcs_test = project_test_external(
                    G_cross, shared.geno_pcs_train, ext_singular_values,
                    reader.n_snps());
            }
        } else if (config.use_gcross_proj && !test_idx.empty()) {
            // G_cross single-pass: read SNPs once, compute G + G_cross simultaneously
            // Avoids second SNP pass at the cost of holding G_cross (n_test × n_train) in memory
            auto* grm_pca = dynamic_cast<GrmPcaDimReduction*>(dim_method.get());
            if (grm_pca) {
                Timer t("GENO PCA fit+proj (G_cross)");
                auto ft_result = grm_pca->fit_and_transform_gcross(
                    reader, train_idx, test_idx,
                    config.n_geno_pcs, config.n_threads);
                shared.geno_dim = std::move(ft_result.params);
                shared.geno_pcs_train = std::move(ft_result.train_pcs);
                shared.geno_pcs_test = std::move(ft_result.test_pcs);
            } else {
                // Fallback (shouldn't happen with default method)
                Timer t("GENO PCA fit");
                shared.geno_dim = dim_method->fit(reader, train_idx,
                                                   config.n_geno_pcs, config.n_threads);
                shared.geno_pcs_train = shared.geno_dim.eigenvectors;
                Timer t2("GENO proj test");
                shared.geno_pcs_test = dim_method->transform(reader, test_idx,
                                                              shared.geno_dim, train_idx, config.n_threads);
            }
        } else {
            // Internal path: 2-pass (fit then transform) for lower peak memory
            // G and G_cross are never alive simultaneously
            {
                Timer t("GENO PCA fit");
                shared.geno_dim = dim_method->fit(reader, train_idx,
                                                   config.n_geno_pcs, config.n_threads);
            }
            // Training PCs = U (raw eigenvectors)
            {
                shared.geno_pcs_train = shared.geno_dim.eigenvectors;
            }
            // Test PCs via cross-GRM (second SNP pass, G_cross ~328 MB, G already released)
            {
                Timer t("GENO proj test");
                shared.geno_pcs_test = dim_method->transform(reader, test_idx,
                                                              shared.geno_dim, train_idx, config.n_threads);
            }
        }

        // --- Shared computation: phenotype PCA (once per fold, if geno_pheno) ---
        const Eigen::MatrixXd* pheno_pcs_train_ptr = nullptr;
        const Eigen::MatrixXd* pheno_pcs_test_ptr = nullptr;
        if (config.correction_mode == CorrectionMode::GENO_PHENO && pheno_train_ptr) {
            std::cout << "  Computing phenotype PCA...\n";
            StandardPcaResult pheno_pca_result;
            {
                Timer t("PHENO PCA fit");
                pheno_pca_result = standard_pca_fit(*pheno_train_ptr,
                                                    config.n_pheno_pcs, config.n_threads);
            }
            {
                Timer t("PHENO proj train");
                shared.pheno_pcs_train = standard_pca_transform(*pheno_train_ptr, pheno_pca_result, config.n_threads);
            }
            {
                Timer t("PHENO proj test");
                shared.pheno_pcs_test = standard_pca_transform(*pheno_test_ptr, pheno_pca_result, config.n_threads);
            }
            shared.pheno_dim.method_name = "standard_pca";
            shared.pheno_dim.n_components = config.n_pheno_pcs;
            shared.pheno_dim.centers = pheno_pca_result.centers;
            shared.pheno_dim.scales = pheno_pca_result.scales;
            shared.pheno_dim.rotation = pheno_pca_result.rotation;
            pheno_pcs_train_ptr = &shared.pheno_pcs_train;
            pheno_pcs_test_ptr = &shared.pheno_pcs_test;
        }

        // --- Build design matrices ---
        {
            Timer t("DESIGN+OLS");
            Eigen::MatrixXd X_train = build_design_matrix(shared.geno_pcs_train, pheno_pcs_train_ptr, cov_train_ptr);
            Eigen::MatrixXd X_test = build_design_matrix(shared.geno_pcs_test, pheno_pcs_test_ptr, cov_test_ptr);

            // Pre-compute (X'X)^{-1} X' for all traits at once
            Eigen::MatrixXd XtX_inv_Xt = compute_ols_hat_matrix(X_train);

            // Collect per-trait fold results for this fold
            std::vector<FoldData> trait_folds(selected_traits.size());
            for (size_t t = 0; t < selected_traits.size(); ++t) {
                int trait = selected_traits[t];

                // Extract trait values
                Eigen::VectorXd y_train(train_pheno_idx.size());
                Eigen::VectorXd y_test(test_pheno_idx.size());
                for (size_t i = 0; i < train_pheno_idx.size(); ++i)
                    y_train(i) = pheno_all_matched(train_pheno_idx[i], trait);
                for (size_t i = 0; i < test_pheno_idx.size(); ++i)
                    y_test(i) = pheno_all_matched(test_pheno_idx[i], trait);

                // OLS: beta = (X'X)^{-1} X'y  (reuse pre-computed)
                Eigen::VectorXd beta = XtX_inv_Xt * y_train;

                // Residuals
                trait_folds[t].fold_idx = f + 1;
                trait_folds[t].y_train_residual = y_train - X_train * beta;
                trait_folds[t].y_test_residual = y_test - X_test * beta;
                trait_folds[t].beta = beta;
            }

            // --- Export this fold ---
            if (config.n_nested_folds <= 0) {
                // 2-way mode
                Timer t2("EXPORT");
                export_single_fold(export_ctx, f + 1, shared, trait_folds, config.save_details_cv);
            } else {
                // ===== Nested CV mode =====
                Timer t2("NESTED CV");

                // Export outer test residuals first
                {
                    Timer t3("TEST EXPORT");
                    export_test_residuals(export_ctx, f + 1, split.test, trait_folds);
                }

                // Export outer fold details (MAF, beta, PCs, eigenvalues) if requested
                if (config.save_details_cv) {
                    Timer t3("DETAIL EXPORT");
                    export_single_fold(export_ctx, f + 1, shared, trait_folds, true);
                }

                int n_nested = config.n_nested_folds;
                int n_pool = static_cast<int>(split.train.size());
                int eff_nested_seed = (config.nested_seed > 0) ? config.nested_seed : config.seed;

                // 1. Create nested folds within the training pool
                std::vector<int> nested_fold_ids;
                if (cov_matched.cols() > 0) {
                    // Stratified nested folds using covariates of training pool
                    Eigen::MatrixXi cov_pool_int(n_pool, cov_matched.cols());
                    for (int i = 0; i < n_pool; ++i)
                        for (int c = 0; c < cov_matched.cols(); ++c)
                            cov_pool_int(i, c) = static_cast<int>(cov_matched(train_pheno_idx[i], c));
                    std::vector<int> factor_cols = detect_factor_cols(cov_pool_int, 10);
                    if (factor_cols.empty()) {
                        std::vector<int> empty_strat;
                        nested_fold_ids = create_folds(n_pool, n_nested, empty_strat, eff_nested_seed);
                    } else {
                        std::vector<std::vector<int>> stratify_cols(factor_cols.size());
                        for (size_t i = 0; i < factor_cols.size(); ++i) {
                            stratify_cols[i].resize(n_pool);
                            for (int j = 0; j < n_pool; ++j)
                                stratify_cols[i][j] = cov_pool_int(j, factor_cols[i]);
                        }
                        nested_fold_ids = create_folds_multi(stratify_cols, n_nested, eff_nested_seed);
                    }
                } else {
                    std::vector<int> empty_strat;
                    nested_fold_ids = create_folds(n_pool, n_nested, empty_strat, eff_nested_seed);
                }

                // 2. Export nested_fold_ids.txt immediately
                export_nested_fold_ids(export_ctx, f + 1, split.train, nested_fold_ids);

                // 3. Initialize train_val_residuals.txt (header only)
                std::string residuals_path = init_train_val_residuals(export_ctx, f + 1, n_nested);

                // 4. Nested loop: compute, export immediately, release memory
                for (int nf = 0; nf < n_nested; ++nf) {
                    std::cout << "    Nested fold " << (nf + 1) << "/" << n_nested << "\n" << std::flush;

                    // Split training pool into nested_train / nested_val
                    std::vector<int> nested_train_local, nested_val_local;
                    for (int i = 0; i < n_pool; ++i) {
                        if (nested_fold_ids[i] == nf + 1)
                            nested_val_local.push_back(i);
                        else
                            nested_train_local.push_back(i);
                    }

                    int n_nt = static_cast<int>(nested_train_local.size());
                    int n_nv = static_cast<int>(nested_val_local.size());

                    // Convert to geno indices and pheno indices
                    std::vector<int> nested_train_geno, nested_val_geno;
                    std::vector<int> nested_train_pheno, nested_val_pheno;
                    for (int idx : nested_train_local) {
                        nested_train_geno.push_back(matched_geno_idx[split.train[idx]]);
                        nested_train_pheno.push_back(split.train[idx]);
                    }
                    for (int idx : nested_val_local) {
                        nested_val_geno.push_back(matched_geno_idx[split.train[idx]]);
                        nested_val_pheno.push_back(split.train[idx]);
                    }

                    // Compute GRM from nested_train
                    DimReductionResult nested_dim;
                    {
                        Timer t3("NEST GRM fit");
                        nested_dim = dim_method->fit(reader, nested_train_geno,
                                                      config.n_geno_pcs, config.n_threads);
                    }

                    // Training PCs = U (raw eigenvectors)
                    Eigen::MatrixXd nested_train_pcs = nested_dim.eigenvectors;

                    // Validation PCs via V-based projection
                    Eigen::MatrixXd nested_val_pcs;
                    {
                        Timer t3("NEST proj val");
                        nested_val_pcs = dim_method->transform(reader, nested_val_geno,
                                                                nested_dim, nested_train_geno,
                                                                config.n_threads);
                    }

                    // Build design matrices for nested train and val
                    Eigen::MatrixXd X_nt, X_nv;

                    // Covariates for nested subsets
                    Eigen::MatrixXd cov_nt, cov_nv;
                    const Eigen::MatrixXd* cov_nt_ptr = nullptr;
                    const Eigen::MatrixXd* cov_nv_ptr = nullptr;
                    if (cov_matched.cols() > 0) {
                        Eigen::VectorXi nti(n_nt), nvi(n_nv);
                        for (int i = 0; i < n_nt; ++i) nti(i) = nested_train_pheno[i];
                        for (int i = 0; i < n_nv; ++i) nvi(i) = nested_val_pheno[i];
                        cov_nt = cov_matched(nti, Eigen::all).cast<double>();
                        cov_nv = cov_matched(nvi, Eigen::all).cast<double>();
                        cov_nt_ptr = &cov_nt;
                        cov_nv_ptr = &cov_nv;
                    }

                    // Phenotype PCA for nested subsets (if geno_pheno mode)
                    Eigen::MatrixXd pheno_nt, pheno_nv;
                    if (config.correction_mode == CorrectionMode::GENO_PHENO && pheno_train_ptr) {
                        Eigen::VectorXi nti(n_nt), nvi(n_nv);
                        for (int i = 0; i < n_nt; ++i) nti(i) = nested_train_pheno[i];
                        for (int i = 0; i < n_nv; ++i) nvi(i) = nested_val_pheno[i];
                        pheno_nt = pheno_all_matched(nti, Eigen::all);
                        pheno_nv = pheno_all_matched(nvi, Eigen::all);
                        StandardPcaResult pheno_pca_nt;
                        pheno_pca_nt = standard_pca_fit(pheno_nt, config.n_pheno_pcs, config.n_threads);
                        Eigen::MatrixXd pp_nt = standard_pca_transform(pheno_nt, pheno_pca_nt, config.n_threads);
                        Eigen::MatrixXd pp_nv = standard_pca_transform(pheno_nv, pheno_pca_nt, config.n_threads);
                        X_nt = build_design_matrix(nested_train_pcs, &pp_nt, cov_nt_ptr);
                        X_nv = build_design_matrix(nested_val_pcs, &pp_nv, cov_nv_ptr);
                        pheno_nt = Eigen::MatrixXd();  // release immediately
                        pheno_nv = Eigen::MatrixXd();
                    } else {
                        X_nt = build_design_matrix(nested_train_pcs, nullptr, cov_nt_ptr);
                        X_nv = build_design_matrix(nested_val_pcs, nullptr, cov_nv_ptr);
                    }

                    // OLS for each trait: beta = (X'X)^{-1} X'y
                    Eigen::MatrixXd XtX_inv_Xt_nt = compute_ols_hat_matrix(X_nt);

                    // Assemble PCs for all pool samples in pool order
                    Eigen::MatrixXd all_pcs(n_pool, X_nt.cols());
                    {
                        int row = 0;
                        for (int local_idx : nested_train_local)
                            all_pcs.row(local_idx) = X_nt.row(row++);
                        row = 0;
                        for (int local_idx : nested_val_local)
                            all_pcs.row(local_idx) = X_nv.row(row++);
                    }

                    // Compute residuals for each trait and append to train_val_residuals.txt
                    {
                        Eigen::VectorXd y_pool(n_pool);
                        Eigen::VectorXd y_nt(n_nt);
                        for (size_t t = 0; t < selected_traits.size(); ++t) {
                            int trait = selected_traits[t];
                            for (int i = 0; i < n_pool; ++i)
                                y_pool(i) = pheno_all_matched(split.train[i], trait);
                            for (int i = 0; i < n_nt; ++i)
                                y_nt(i) = pheno_all_matched(nested_train_pheno[i], trait);

                            Eigen::VectorXd beta_nt = XtX_inv_Xt_nt * y_nt;
                            Eigen::VectorXd residuals = y_pool - all_pcs * beta_nt;

                            std::string suffix = selected_traits.size() > 1
                                ? ("_" + pheno.trait_names[trait]) : "";
                            append_train_val_residual_column(
                                residuals_path, export_ctx.sample_ids, split.train,
                                nf + 1, residuals, suffix);
                        }
                    }

                    // Export detail files immediately if requested
                    if (config.save_details_nested) {
                        SharedFoldData ns;
                        ns.fold_idx = nf + 1;
                        ns.train_idx = nested_train_local;
                        ns.test_idx = nested_val_local;
                        ns.geno_pcs_train = nested_train_pcs;
                        ns.geno_pcs_test = nested_val_pcs;
                        ns.geno_dim = nested_dim;

                        std::vector<FoldData> ntf(selected_traits.size());
                        for (size_t t = 0; t < selected_traits.size(); ++t) {
                            int trait = selected_traits[t];
                            Eigen::VectorXd y_nt(n_nt);
                            for (int i = 0; i < n_nt; ++i)
                                y_nt(i) = pheno_all_matched(nested_train_pheno[i], trait);
                            Eigen::VectorXd beta_nt = XtX_inv_Xt_nt * y_nt;
                            ntf[t].fold_idx = nf + 1;
                            ntf[t].beta = beta_nt;
                        }

                        export_nested_fold_details(export_ctx, f + 1, nf + 1, ns, ntf);
                        // Detail data is now written to disk; local copies freed below
                    }

                    // Release nested fold memory immediately
                    nested_train_pcs = Eigen::MatrixXd();
                    nested_val_pcs = Eigen::MatrixXd();
                    X_nt = Eigen::MatrixXd();
                    X_nv = Eigen::MatrixXd();
                    cov_nt = Eigen::MatrixXd();
                    cov_nv = Eigen::MatrixXd();
                    all_pcs = Eigen::MatrixXd();
                    nested_dim.eigenvectors = Eigen::MatrixXd();
                    nested_dim.eigenvalues = Eigen::VectorXd();
                    nested_dim.maf = Eigen::VectorXd();
                    XtX_inv_Xt_nt = Eigen::MatrixXd();
                }

                std::cout << "[INFO] Nested CV results saved to: fold_" << (f + 1) << "\n";
            }
        }

        // --- Release fold-specific memory ---
        shared.geno_pcs_train = Eigen::MatrixXd();
        shared.geno_pcs_test = Eigen::MatrixXd();
        shared.pheno_pcs_train = Eigen::MatrixXd();
        shared.pheno_pcs_test = Eigen::MatrixXd();
        shared.geno_dim.eigenvectors = Eigen::MatrixXd();
        shared.geno_dim.eigenvalues = Eigen::VectorXd();
        shared.geno_dim.maf = Eigen::VectorXd();
        shared.geno_dim.centers = Eigen::VectorXd();
        shared.geno_dim.scales = Eigen::VectorXd();
        shared.geno_dim.rotation = Eigen::MatrixXd();
        shared.pheno_dim = DimReductionResult();
        cov_train = Eigen::MatrixXd();
        cov_test = Eigen::MatrixXd();
        pheno_train = Eigen::MatrixXd();
        pheno_test = Eigen::MatrixXd();

#ifdef FASTCV_USE_MKL
        mkl_free_buffers();
#endif

        auto fold_end = std::chrono::high_resolution_clock::now();
        double fold_sec = std::chrono::duration<double>(fold_end - fold_start).count();
        total_fold_time += fold_sec;
        std::cout << "  Fold " << (f + 1) << " total: " << fold_sec << "s"
                  << " (cumulative: " << total_fold_time << "s)\n" << std::flush;
    }

    // ========================================================
    // Summary
    // ========================================================
    auto pipeline_end = std::chrono::high_resolution_clock::now();
    double total_sec = std::chrono::duration<double>(pipeline_end - pipeline_start).count();
    double data_load_sec = std::chrono::duration<double>(
        // data load ended when fold loop started
        total_sec - total_fold_time
    ).count();

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "DATA PREPARATION COMPLETE\n";
    std::cout << "  Traits: " << selected_traits.size() << " | Folds: " << config.n_folds << "\n";
    std::cout << "  Fold processing:  " << total_fold_time << "s\n";
    std::cout << "  Total wall time:   " << total_sec << "s\n";
    std::cout << "  Avg time/fold:     " << (total_fold_time / config.n_folds) << "s\n";
    std::cout << "  Output: " << config.output_dir << "\n";
    std::cout << std::string(60, '=') << "\n\n";

    reader.close();

    // Restore cout before TeeStream is destroyed (prevents use-after-free segfault)
    if (log_tee) {
        std::cout.rdbuf(log_tee->original_buf());
    }

    return result;
}

} // namespace fastcv
