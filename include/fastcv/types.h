#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

namespace fastcv {

// ============================================================
// Forward declarations
// ============================================================
class PlinkReader;

// ============================================================
// Enums
// ============================================================
enum class CorrectionMode { GENO_ONLY, GENO_PHENO };

// ============================================================
// SNP metadata
// ============================================================
struct SnpInfo {
    std::string chr;
    std::string snp_id;
    double cm = 0.0;
    double pos = 0.0;
    std::string a1;  // minor allele
    std::string a2;  // major allele
};

// ============================================================
// Reusable I/O buffers for dual-block reads (avoids malloc per block)
// ============================================================
struct DualReadBuf {
    std::vector<unsigned char> raw_buf;
    std::vector<int> geno_a;
    std::vector<int> geno_b;
    std::vector<std::pair<int,int>> loc_a;
    std::vector<std::pair<int,int>> loc_b;
};

// ============================================================
// PLINK streaming reader (never loads full n×m matrix)
// ============================================================
class PlinkReader {
public:
    PlinkReader() = default;
    ~PlinkReader();

    // Non-copyable, movable
    PlinkReader(const PlinkReader&) = delete;
    PlinkReader& operator=(const PlinkReader&) = delete;
    PlinkReader(PlinkReader&&) noexcept;
    PlinkReader& operator=(PlinkReader&&) noexcept;

    /// Open PLINK .bed/.bim/.fam files (infer .bim/.fam from .bed path)
    void open(const std::string& bed_file);

    /// Check if files are loaded
    bool is_open() const { return bed_fd_ >= 0; }

    int n_samples() const { return n_samples_; }
    int n_snps() const { return n_snps_; }
    const std::vector<std::string>& sample_ids() const { return sample_ids_; }
    const std::vector<SnpInfo>& snp_info() const { return snp_info_; }

    /// Read a single SNP for specified sample indices.
    /// Returns genotype vector (n_selected), coded as 0/1/2/-9 (missing).
    Eigen::VectorXi read_snp(int snp_idx,
                             const std::vector<int>& sample_idx) const;

    /// Read a block of SNPs for specified sample indices.
    /// Z_block: (n_selected × actual_block_size), center/scaled
    /// maf_block: (actual_block_size), allele frequencies (actual freq, not MAF)
    /// If external_maf is provided, use it for centering/scaling instead of
    /// computing from the data. This is needed for cross-GRM where test
    /// samples must be standardized using training-set allele frequencies.
    /// If center_only is true, only center (g - 2p) without dividing by sqrt(2p(1-p)).
    void read_snp_block(int snp_start, int snp_end,
                        const std::vector<int>& sample_idx,
                        Eigen::MatrixXd& Z_block,
                        Eigen::VectorXd& maf_block,
                        const Eigen::VectorXd* external_maf = nullptr,
                        bool center_only = false) const;

    /// Read a block of SNPs for TWO sample sets in a single I/O pass.
    /// MAF is computed from sample_idx_a only; both Z_a and Z_b are
    /// center/scaled using that MAF (prevents data leakage).
    /// If center_only is true, only center (g - 2p) without dividing by sqrt(2p(1-p)).
    void read_snp_block_dual(int snp_start, int snp_end,
                             const std::vector<int>& sample_idx_a,
                             const std::vector<int>& sample_idx_b,
                             Eigen::MatrixXd& Z_a,
                             Eigen::MatrixXd& Z_b,
                             Eigen::VectorXd& maf_block,
                             bool center_only = false) const;

    /// Compute MAF for all SNPs given sample indices
    Eigen::VectorXd compute_maf(const std::vector<int>& sample_idx) const;

    /// Buffered version of read_snp_block_dual: reuses pre-allocated DualReadBuf.
    /// Z_a, Z_b, maf_block must be pre-allocated at max block size; only the
    /// first (snp_end - snp_start) columns/elements are written.
    void read_snp_block_dual_buf(int snp_start, int snp_end,
                                 const std::vector<int>& sample_idx_a,
                                 const std::vector<int>& sample_idx_b,
                                 Eigen::MatrixXd& Z_a,
                                 Eigen::MatrixXd& Z_b,
                                 Eigen::VectorXd& maf_block,
                                 DualReadBuf& buf,
                                 bool center_only = false) const;

    void close();

    /// Read raw SNP bytes via pread (thread-safe, no mmap).
    /// Returns bytes_per_snp bytes into user-provided buffer.
    void read_snp_raw(int snp_idx, unsigned char* buf) const;

    /// BED file path
    const std::string& bed_file_path() const { return bed_file_path_; }
    long bed_data_offset() const { return bed_data_offset_; }
    int bytes_per_snp() const { return (n_samples_ + 3) / 4; }

private:
    void read_fam_(const std::string& fam_file);
    void read_bim_(const std::string& bim_file);
    void validate_bed_();
    void decode_snp_(int snp_idx, const std::vector<int>& sample_idx,
                     Eigen::VectorXi& geno) const;

    int bed_fd_ = -1;                        // file descriptor for pread
    long bed_file_size_ = 0;                 // total file size
    long bed_data_offset_ = 0;               // byte offset where genotype data starts (3)
    int n_samples_ = 0;
    int n_snps_ = 0;
    std::string bed_file_path_;
    std::vector<std::string> sample_ids_;
    std::vector<SnpInfo> snp_info_;
};

// ============================================================
// Phenotype data
// ============================================================
struct PhenoData {
    Eigen::VectorXd y;                 // target trait
    Eigen::MatrixXd pheno_all;         // full phenotype matrix (multi-trait mode)
    std::vector<std::string> sample_ids;
    std::vector<std::string> trait_names;
    bool is_multi_pheno = false;
};

// ============================================================
// Dimensionality reduction result
// ============================================================
struct DimReductionResult {
    std::string method_name;
    int n_components = 0;

    // GRM-PCA: store eigenvectors (U) and eigenvalues (lambda), not V
    Eigen::MatrixXd eigenvectors;      // U (n × k)
    Eigen::VectorXd eigenvalues;       // lambda (k)

    // Per-feature statistics for centering/scaling new data
    Eigen::VectorXd centers;           // per-column means
    Eigen::VectorXd scales;            // per-column std devs
    Eigen::VectorXd maf;               // allele frequencies (actual freq, not MAF)

    // Sum of per-SNP variances: Σ 2p(1-p), used as GRM normalizer
    double sigma2_sum = 0.0;

    // For phenotype PCA: store rotation matrix V directly (pheno matrix is small)
    Eigen::MatrixXd rotation;          // V (p × k), only for phenotype PCA
};

// ============================================================
// Per-fold shared data (identical across traits)
// ============================================================
struct SharedFoldData {
    int fold_idx = 0;
    std::vector<int> train_idx;
    std::vector<int> test_idx;

    // PCA scores (n × k)
    Eigen::MatrixXd geno_pcs_train;
    Eigen::MatrixXd geno_pcs_test;
    Eigen::MatrixXd pheno_pcs_train;   // empty if not geno_pheno
    Eigen::MatrixXd pheno_pcs_test;

    // DimReduction parameters
    DimReductionResult geno_dim;       // contains MAF
    DimReductionResult pheno_dim;
};

// ============================================================
// Per-trait per-fold data
// ============================================================
struct FoldData {
    int fold_idx = 0;

    // Corrected residuals
    Eigen::VectorXd y_train_residual;
    Eigen::VectorXd y_test_residual;

    // Correction coefficients
    Eigen::VectorXd beta;
};

// ============================================================
// CV configuration
// ============================================================
struct CvConfig {
    // Trait selection: "1", "all", "2-10", "1,3,5-7" (used if trait_name is empty)
    std::string trait_expr = "1";
    // Trait selection by column name (takes priority over trait_expr if non-empty)
    std::string trait_name;

    // CV settings
    int n_folds = 5;
    int seed = 42;
    // Stratification: covariate column name used for stratified fold assignment.
    // Empty = no stratification (or auto-detect factor columns from covariate file).
    std::string stratified_cov_name;

    // PC correction
    int n_geno_pcs = 5;
    int n_pheno_pcs = 10;
    CorrectionMode correction_mode = CorrectionMode::GENO_ONLY;

    // Dimensionality reduction method: "pca" (default, GRM-based) or "grm_pca" (legacy)
    std::string dim_reduction_method = "pca";

    // External PCA backend: "internal" (default), "hiblup", "plink"
    std::string pca_backend = "internal";
    std::string pca_tool_path;           // path to hiblup/plink binary (empty = search PATH)

    // GRM block size (SNPs per block). 0 = auto-detect based on L3 cache.
    int grm_block_size = 0;

    // LD pruning (optional, disabled by default)
    bool ld_prune = false;
    double ld_r2_threshold = 0.2;   // maximum allowed r²
    int ld_window_size = 50;        // sliding window size (SNPs)
    double ld_maf_min = 0.01;       // minimum allele frequency

    // Nested cross-validation
    int n_nested_folds = 0;           // 0 = no nested CV (2-way), >0 = nested folds
    int nested_seed = 0;              // 0 = same as seed

    // Save intermediate details (MAF/beta/PCs per fold)
    bool save_details_cv = false;      // save details for outer CV folds
    bool save_details_nested = false;  // save details for nested CV folds

    // G_cross single-pass projection (reads SNPs once instead of twice)
    bool use_gcross_proj = false;      // true = single-pass G_cross, false = 2-pass V-proj

    // Performance
    int n_threads = 0;                 // 0 = auto-detect

    // I/O
    std::string output_dir = "cv_output";
    std::string separator = " ";       // output file field separator

    // Data paths
    std::string bed_file;
    std::string pheno_file;
    std::string pheno_id_col = "sample";
    std::string cov_file;              // covariate file (optional)
    std::string cov_id_col = "sample"; // ID column in covariate file
    // Covariate column selection: "all" (default) or comma-separated names/indices
    std::string cov_select = "all";
};

// ============================================================
// Per-trait result
// ============================================================
struct TraitResult {
    std::string trait_name;
    std::vector<FoldData> folds;
};

// ============================================================
// CV result
// ============================================================
struct CvResult {
    struct Meta {
        std::string dim_reduction_method;
        int n_samples = 0;
        int n_snps = 0;
        int n_folds = 0;
        int n_geno_pcs = 0;
        int n_pheno_pcs = 0;
        CorrectionMode correction_mode = CorrectionMode::GENO_ONLY;
        int seed = 0;
        std::string stratify_col;
        std::string timestamp;
    } config;

    std::vector<int> fold_ids;
    std::vector<std::string> sample_ids;
    std::vector<std::string> snp_ids;
    std::vector<std::string> cov_col_names;
    std::vector<SharedFoldData> shared_folds;  // per-fold shared (PCs, MAF)
    std::vector<TraitResult> traits;             // per-trait fold data
};

// ============================================================
// Correction fit result
// ============================================================
struct CorrectionFitResult {
    Eigen::VectorXd y_residual;
    DimReductionResult geno_dim;
    DimReductionResult pheno_dim;
    Eigen::VectorXd beta;
};

// ============================================================
// Utility: parse correction mode from string
// ============================================================
inline CorrectionMode parse_correction_mode(const std::string& s) {
    if (s == "geno_pheno") return CorrectionMode::GENO_PHENO;
    return CorrectionMode::GENO_ONLY;  // default
}

inline std::string correction_mode_to_string(CorrectionMode m) {
    return (m == CorrectionMode::GENO_PHENO) ? "geno_pheno" : "geno_only";
}

/// Parse trait expression → 0-indexed trait indices.
/// Supports: "1", "all", "2-10", "1,3,5-7"
inline std::vector<int> parse_trait_expr(const std::string& expr, int n_traits) {
    std::vector<int> result;
    if (expr == "all") {
        for (int i = 0; i < n_traits; ++i) result.push_back(i);
        return result;
    }
    // Parse comma-separated parts: "1,3,5-7"
    std::istringstream ss(expr);
    std::string part;
    while (std::getline(ss, part, ',')) {
        auto dash = part.find('-');
        if (dash != std::string::npos) {
            int a = std::stoi(part.substr(0, dash));
            int b = std::stoi(part.substr(dash + 1));
            for (int i = a; i <= b; ++i) {
                if (i >= 1 && i <= n_traits) result.push_back(i - 1);
            }
        } else {
            int v = std::stoi(part);
            if (v >= 1 && v <= n_traits) result.push_back(v - 1);
        }
    }
    return result;
}

} // namespace fastcv
