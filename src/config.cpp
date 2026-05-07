#include "fastcv/config.h"
#include "json.hpp"
#include <fstream>
#include <stdexcept>

namespace fastcv {

using json = nlohmann::json;

CvConfig load_config_from_string(const std::string& json_str) {
    json j = json::parse(json_str);
    CvConfig config;

    // ============================================================
    // Helper: resolve path relative to data_dir
    // ============================================================
    std::string data_dir;
    if (j.contains("paths") && j["paths"].contains("data_dir")) {
        data_dir = j["paths"]["data_dir"];
    }

    auto resolve_path = [&data_dir](const std::string& p) -> std::string {
        if (p.empty()) return p;
        if (p[0] == '/' || data_dir.empty()) return p;
        return data_dir + "/" + p;
    };

    // ============================================================
    // Data paths (new format: data.genotype / data.phenotype / data.covariate)
    // ============================================================
    if (j.contains("data")) {
        auto& data = j["data"];

        // Genotype
        if (data.contains("genotype")) {
            auto& geno = data["genotype"];
            if (geno.contains("bed_file"))
                config.bed_file = resolve_path(geno["bed_file"]);
        }

        // Phenotype (new field names: pheno_file, id_col_name, trait_idx/trait_name)
        if (data.contains("phenotype")) {
            auto& ph = data["phenotype"];
            if (ph.contains("pheno_file"))
                config.pheno_file = resolve_path(ph["pheno_file"]);
            else if (ph.contains("file"))
                config.pheno_file = resolve_path(ph["file"]);  // backward compat

            if (ph.contains("id_col_name"))
                config.pheno_id_col = ph["id_col_name"];
            else if (ph.contains("id_col"))
                config.pheno_id_col = ph["id_col"];  // backward compat

            if (ph.contains("trait_name"))
                config.trait_name = ph["trait_name"];
            else if (ph.contains("trait_idx")) {
                // trait_idx can be int or string
                if (ph["trait_idx"].is_number())
                    config.trait_expr = std::to_string(ph["trait_idx"].get<int>());
                else
                    config.trait_expr = ph["trait_idx"];
            } else if (ph.contains("traits")) {
                config.trait_expr = ph["traits"];  // backward compat
            }
        }

        // Covariate (new: cov_file, id_col_name, cov_name/cov_idx)
        if (data.contains("covariate")) {
            auto& cov = data["covariate"];
            if (cov.contains("cov_file"))
                config.cov_file = resolve_path(cov["cov_file"]);
            else if (cov.contains("file"))
                config.cov_file = resolve_path(cov["file"]);  // backward compat

            if (cov.contains("id_col_name"))
                config.cov_id_col = cov["id_col_name"];

            if (cov.contains("cov_name"))
                config.cov_select = cov["cov_name"];
            else if (cov.contains("cov_idx"))
                config.cov_select = cov["cov_idx"];
        }
    }

    // ============================================================
    // Cross-validation
    // ============================================================
    if (j.contains("cross_validation")) {
        auto& cv = j["cross_validation"];

        // cv sub-object (new format)
        if (cv.contains("cv")) {
            auto& cv_main = cv["cv"];
            if (cv_main.contains("n_folds")) config.n_folds = cv_main["n_folds"];
            if (cv_main.contains("seed")) config.seed = cv_main["seed"];
        } else {
            // Backward compat: n_folds and seed at cross_validation level
            if (cv.contains("n_folds")) config.n_folds = cv["n_folds"];
            if (cv.contains("seed")) config.seed = cv["seed"];
        }

        // nested_cv sub-object (new format)
        if (cv.contains("nested_cv")) {
            auto& ncv = cv["nested_cv"];
            if (ncv.contains("n_nested_folds")) config.n_nested_folds = ncv["n_nested_folds"];
            if (ncv.contains("nested_seed")) config.nested_seed = ncv["nested_seed"];
        } else {
            // Backward compat: at cross_validation level
            if (cv.contains("n_nested_folds")) config.n_nested_folds = cv["n_nested_folds"];
            if (cv.contains("nested_seed")) config.nested_seed = cv["nested_seed"];
        }

        // save_details sub-object (new format)
        if (cv.contains("save_details")) {
            auto& sd = cv["save_details"];
            if (sd.contains("cv")) config.save_details_cv = sd["cv"];
            if (sd.contains("nested_cv")) config.save_details_nested = sd["nested_cv"];
        } else if (cv.contains("save_nested_details")) {
            // Backward compat: save_nested_details → save_details_nested
            config.save_details_nested = cv["save_nested_details"];
        }

        // Stratified covariate name
        if (cv.contains("stratified_cov_name"))
            config.stratified_cov_name = cv["stratified_cov_name"];

        // Correction
        if (cv.contains("correction")) {
            auto& corr = cv["correction"];
            if (corr.contains("mode"))
                config.correction_mode = parse_correction_mode(corr["mode"]);
            if (corr.contains("n_geno_pcs")) config.n_geno_pcs = corr["n_geno_pcs"];
            if (corr.contains("n_pheno_pcs")) config.n_pheno_pcs = corr["n_pheno_pcs"];
        }

        // Dimensionality reduction
        if (cv.contains("dim_reduction")) {
            auto& dr = cv["dim_reduction"];
            if (dr.contains("method"))
                config.dim_reduction_method = dr["method"];
            if (dr.contains("pca_backend"))
                config.pca_backend = dr["pca_backend"];
            if (dr.contains("pca_tool_path"))
                config.pca_tool_path = dr["pca_tool_path"];
            if (dr.contains("grm_block_size"))
                config.grm_block_size = dr["grm_block_size"];
            if (dr.contains("ld_prune"))
                config.ld_prune = dr["ld_prune"];
            if (dr.contains("ld_r2_threshold"))
                config.ld_r2_threshold = dr["ld_r2_threshold"];
            if (dr.contains("ld_window_size"))
                config.ld_window_size = dr["ld_window_size"];
            if (dr.contains("ld_maf_min"))
                config.ld_maf_min = dr["ld_maf_min"];
            if (dr.contains("use_gcross_proj"))
                config.use_gcross_proj = dr["use_gcross_proj"];
        }
    }

    // ============================================================
    // Other settings (new format)
    // ============================================================
    bool output_dir_explicit = false;
    if (j.contains("other")) {
        auto& other = j["other"];
        if (other.contains("separator")) config.separator = other["separator"];
        if (other.contains("n_threads")) config.n_threads = other["n_threads"];
        if (other.contains("output_dir")) {
            config.output_dir = other["output_dir"];
            output_dir_explicit = true;
        }
    }

    // ============================================================
    // Top-level fallbacks (backward compat)
    // ============================================================
    if (config.n_threads == 0 && j.contains("n_threads"))
        config.n_threads = j["n_threads"];

    // paths.output_dir: apply if "other" section didn't specify output_dir
    if (!output_dir_explicit && j.contains("paths") && j["paths"].contains("output_dir"))
        config.output_dir = j["paths"]["output_dir"];

    // Old "output" section → ignore (format was csv-only anyway)

    return config;
}

CvConfig load_config(const std::string& json_file) {
    std::ifstream fin(json_file);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open config file: " + json_file);
    }
    std::string content((std::istreambuf_iterator<char>(fin)),
                         std::istreambuf_iterator<char>());
    auto config = load_config_from_string(content);
    validate_config(config);
    return config;
}

void validate_config(const CvConfig& config) {
    std::vector<std::string> errors;

    if (config.bed_file.empty()) {
        errors.push_back("genotype bed_file is required");
    }
    if (config.pheno_file.empty()) {
        errors.push_back("phenotype file is required");
    }
    if (config.n_folds < 2) {
        errors.push_back("n_folds must be >= 2");
    }
    if (config.n_geno_pcs < 1) {
        errors.push_back("n_geno_pcs must be >= 1");
    }
    if (config.n_nested_folds < 0) {
        errors.push_back("n_nested_folds must be >= 0 (0 = no nested CV)");
    }
    if (config.correction_mode == CorrectionMode::GENO_PHENO && config.n_pheno_pcs < 1) {
        errors.push_back("n_pheno_pcs must be >= 1 in geno_pheno mode");
    }
    if (config.dim_reduction_method != "grm_pca" && config.dim_reduction_method != "pca") {
        errors.push_back("Unsupported dim_reduction method: " + config.dim_reduction_method);
    }

    if (!errors.empty()) {
        std::string msg = "Configuration validation failed:\n";
        for (const auto& e : errors) msg += "  - " + e + "\n";
        throw std::runtime_error(msg);
    }
}

std::string generate_template_config() {
    json j;

    // Top-level descriptions (ignored by load_config_from_string)
    j["_comment"] = {
        {"data.genotype.bed_file",               "PLINK .bed file path (required)"},
        {"data.phenotype.pheno_file",            "Phenotype file path, whitespace-delimited (required)"},
        {"data.phenotype.id_col_name",           "Sample ID column name in phenotype file"},
        {"data.phenotype.trait_name",            "Trait column name (takes priority over trait_idx, empty = use trait_idx)"},
        {"data.phenotype.trait_idx",             "Trait selection: '1', 'all', '2-10', '1,3,5-7'"},
        {"data.covariate.cov_file",              "Covariate file (optional, empty = no covariates)"},
        {"data.covariate.id_col_name",           "Sample ID column name in covariate file"},
        {"data.covariate.cov_name",              "Covariate selection: 'all' or comma-separated names/indices"},
        {"cross_validation.cv.n_folds",          "Number of CV folds (>= 2)"},
        {"cross_validation.cv.seed",             "Random seed for fold assignment"},
        {"cross_validation.nested_cv.n_nested_folds", "Nested CV folds (0 = no nested CV)"},
        {"cross_validation.nested_cv.nested_seed",    "Nested CV seed (0 = same as seed)"},
        {"cross_validation.correction.mode",     "Correction mode: 'geno_only' or 'geno_pheno'"},
        {"cross_validation.correction.n_geno_pcs",    "Number of genotype PCs for correction (>= 1)"},
        {"cross_validation.correction.n_pheno_pcs",   "Number of phenotype PCs (geno_pheno mode, >= 1)"},
        {"cross_validation.dim_reduction.method",      "Dim reduction: 'pca' (default, GRM-based) or 'grm_pca' (legacy)"},
        {"cross_validation.dim_reduction.pca_backend", "PCA backend: 'internal' (default), 'hiblup', 'plink'"},
        {"cross_validation.dim_reduction.pca_tool_path", "Path to hiblup/plink binary (empty = search PATH)"},
        {"cross_validation.dim_reduction.grm_block_size", "GRM block size in SNPs (0 = auto-detect from L3 cache)"},
        {"cross_validation.dim_reduction.ld_prune",         "Enable LD pruning before GRM (true/false)"},
        {"cross_validation.dim_reduction.ld_r2_threshold",  "Max r-squared for LD pruning"},
        {"cross_validation.dim_reduction.ld_window_size",   "Sliding window size in SNPs for LD pruning"},
        {"cross_validation.dim_reduction.ld_maf_min",       "Minimum allele frequency for LD pruning"},
        {"cross_validation.dim_reduction.use_gcross_proj",  "G_cross single-pass (reads SNPs once, higher peak RAM)"},
        {"cross_validation.save_details.cv",           "Save MAF/beta/PCs per outer CV fold"},
        {"cross_validation.save_details.nested_cv",    "Save details per nested CV fold"},
        {"cross_validation.stratified_cov_name",       "Covariate column for stratified fold assignment (empty = none)"},
        {"other.n_threads",    "OpenMP threads (0 = auto-detect)"},
        {"other.output_dir",   "Output directory for CV results"},
        {"other.separator",    "Output file field separator"}
    };

    j["data"] = {
        {"genotype", {{"bed_file", "genotype.bed"}}},
        {"phenotype", {
            {"pheno_file", "phenotype.txt"},
            {"id_col_name", "sample"},
            {"trait_name", ""},
            {"trait_idx", "1"}
        }},
        {"covariate", {
            {"cov_file", ""},
            {"id_col_name", "sample"},
            {"cov_name", "all"}
        }}
    };
    j["cross_validation"] = {
        {"cv", {
            {"n_folds", 5},
            {"seed", 42}
        }},
        {"nested_cv", {
            {"n_nested_folds", 0},
            {"nested_seed", 0}
        }},
        {"correction", {
            {"mode", "geno_only"},
            {"n_geno_pcs", 5},
            {"n_pheno_pcs", 10}
        }},
        {"dim_reduction", {
            {"method", "pca"},
            {"pca_backend", "internal"},
            {"pca_tool_path", ""},
            {"grm_block_size", 0},
            {"ld_prune", false},
            {"ld_r2_threshold", 0.2},
            {"ld_window_size", 50},
            {"ld_maf_min", 0.01},
            {"use_gcross_proj", false}
        }},
        {"save_details", {
            {"cv", false},
            {"nested_cv", false}
        }},
        {"stratified_cov_name", ""}
    };
    j["other"] = {
        {"separator", " "},
        {"n_threads", 0},
        {"output_dir", "cv_output"}
    };
    return j.dump(2);
}

} // namespace fastcv
