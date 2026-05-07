#include "fastcv/cli.h"
#include "fastcv/config.h"
#include "fastcv/pc_corrector.h"
#include "fastcv/memory_estimate.h"
#include "fastcv/types.h"

#include "fastcv/cross_validation.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

namespace fastcv {

namespace {

void print_help() {
    std::cout << R"(
fastcv - Fast Cross-Validation Data Preparation for Genomic Prediction
=====================================================================

Usage: fastcv <command> [options]

Commands:
  prepare    Prepare CV data with fold-wise PC correction
  validate   Validate configuration file
  config     Generate configuration template
  help       Show this help message

Options:
  -c, --config <file>     Configuration file (JSON)
  -o, --output <file>     Output file for template

Examples:
  fastcv prepare -c config.json
  fastcv validate -c config.json
  fastcv config -o my_config.json
  fastcv help

)";
}

int cmd_prepare(const std::string& config_file) {
    // NUMA-aware thread binding to prevent cross-NUMA memory access
    // These MUST be set before the first OpenMP parallel region.
    // Use overwrite=1 to override any existing (possibly wrong) settings.
    setenv("OMP_PROC_BIND", "close", 1);
    setenv("OMP_PLACES", "cores", 1);
    // Intel MKL affinity (only effective with Intel OpenMP runtime)
    setenv("KMP_AFFINITY", "compact", 0);

    std::cout << "[INFO] Loading configuration: " << config_file << "\n";
    auto config = load_config(config_file);
    prepare_cv_data(config);  // exports each fold inline

    return 0;
}

int cmd_validate(const std::string& config_file) {
    std::cout << "[INFO] Validating configuration: " << config_file << "\n";
    try {
        auto config = load_config(config_file);
        std::cout << "\n  Configuration is valid!\n";

        // === Dataset info ===
        int n_samples = 0, n_snps = 0;
        bool bed_available = false;
        PlinkReader reader;
        if (!config.bed_file.empty()) {
            try {
                reader.open(config.bed_file);
                n_samples = reader.n_samples();
                n_snps = reader.n_snps();
                bed_available = true;
            } catch (const std::exception& e) {
                std::cout << "  [WARN] Cannot open BED: " << e.what() << "\n";
            }
        }

        std::cout << "\n=== Dataset ===\n";
        if (bed_available) {
            std::cout << "  Bed file:       " << config.bed_file << "\n";
            std::cout << "  Samples:        " << n_samples << "\n";
            std::cout << "  SNPs:           " << n_snps << "\n";
        } else {
            std::cout << "  Bed file:       " << config.bed_file << " (NOT FOUND)\n";
            std::cout << "  Samples:        unknown\n";
            std::cout << "  SNPs:           unknown\n";
        }
        std::cout << "  Pheno file:     " << config.pheno_file << "\n";
        if (!config.cov_file.empty())
            std::cout << "  Cov file:       " << config.cov_file << "\n";

        // === CV config ===
        std::cout << "\n=== Cross-Validation ===\n";
        std::cout << "  Folds:          " << config.n_folds << "\n";
        std::cout << "  Seed:           " << config.seed << "\n";
        std::cout << "  Correction:     " << correction_mode_to_string(config.correction_mode) << "\n";
        std::cout << "  Genotype PCs:   " << config.n_geno_pcs << "\n";
        if (config.correction_mode == CorrectionMode::GENO_PHENO)
            std::cout << "  Phenotype PCs:  " << config.n_pheno_pcs << "\n";
        std::cout << "  DimReduction:   " << config.dim_reduction_method << "\n";
        if (config.n_nested_folds > 0)
            std::cout << "  Nested folds:   " << config.n_nested_folds << "\n";
        if (config.ld_prune)
            std::cout << "  LD prune:       ON (r2<" << config.ld_r2_threshold
                      << ", window=" << config.ld_window_size << ")\n";

        // === Memory prediction ===
        MemoryEstimate est = estimate_memory(
            n_samples, n_snps,
            config.n_folds, config.n_geno_pcs,
            config.n_threads, config.grm_block_size,
            config.use_gcross_proj);

        print_memory_report(est);

        if (bed_available) reader.close();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nValidation FAILED: " << e.what() << "\n";
        return 1;
    }
}

int cmd_config(const std::string& output_file) {
    std::string template_json = generate_template_config();
    if (output_file.empty()) {
        std::cout << template_json << "\n";
    } else {
        std::ofstream out(output_file);
        if (!out.is_open()) {
            std::cerr << "Cannot write to: " << output_file << "\n";
            return 1;
        }
        out << template_json << "\n";
        std::cout << "Template configuration written to: " << output_file << "\n";
    }
    return 0;
}

} // anonymous namespace

// ============================================================
// run_cli()
// ============================================================
int run_cli(const std::vector<std::string>& args) {
    if (args.empty()) {
        print_help();
        return 0;
    }

    const std::string& cmd = args[0];

    if (cmd == "help" || cmd == "--help" || cmd == "-h") {
        print_help();
        return 0;
    }

    if (cmd == "prepare" || cmd == "run") {
        std::string config_file;
        for (size_t i = 1; i < args.size(); ++i) {
            if ((args[i] == "-c" || args[i] == "--config") && i + 1 < args.size()) {
                config_file = args[++i];
            }
        }
        if (config_file.empty()) {
            std::cerr << "Error: configuration file required. Use -c <file>\n";
            return 1;
        }
        return cmd_prepare(config_file);
    }

    if (cmd == "validate") {
        std::string config_file;
        for (size_t i = 1; i < args.size(); ++i) {
            if ((args[i] == "-c" || args[i] == "--config") && i + 1 < args.size()) {
                config_file = args[++i];
            }
        }
        if (config_file.empty()) {
            std::cerr << "Error: configuration file required. Use -c <file>\n";
            return 1;
        }
        return cmd_validate(config_file);
    }

    if (cmd == "config") {
        std::string output_file;
        for (size_t i = 1; i < args.size(); ++i) {
            if ((args[i] == "-o" || args[i] == "--output") && i + 1 < args.size()) {
                output_file = args[++i];
            }
        }
        return cmd_config(output_file);
    }

    std::cerr << "Unknown command: " << cmd << "\n";
    std::cerr << "Run 'fastcv help' for usage information.\n";
    return 1;
}

} // namespace fastcv
