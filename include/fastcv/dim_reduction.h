#pragma once

#include "fastcv/types.h"
#include <memory>

namespace fastcv {

// ============================================================
// Abstract interface for dimensionality reduction methods
// ============================================================
class DimReduction {
public:
    virtual ~DimReduction() = default;

    /// Fit on training data (streaming via PlinkReader)
    /// @param reader      PLINK reader (must be open)
    /// @param sample_idx  Training sample indices
    /// @param n_components  Number of components to extract
    /// @param n_threads   OpenMP threads (0 = auto)
    virtual DimReductionResult fit(const PlinkReader& reader,
                                   const std::vector<int>& sample_idx,
                                   int n_components,
                                   int n_threads = 0) = 0;

    /// Project new data using fitted parameters (streaming via PlinkReader)
    /// @param reader      PLINK reader (must be open)
    /// @param new_idx     New sample indices
    /// @param params      Fitted parameters from fit()
    /// @param train_idx   Training sample indices (needed for cross-GRM)
    /// @param n_threads   OpenMP threads (0 = auto)
    /// @return n_new × n_components projection matrix
    virtual Eigen::MatrixXd transform(const PlinkReader& reader,
                                      const std::vector<int>& new_idx,
                                      const DimReductionResult& params,
                                      const std::vector<int>& train_idx,
                                      int n_threads = 0) const = 0;

    virtual std::string name() const = 0;
};

// ============================================================
// GRM-based PCA (default implementation)
//
// Algorithm:
//   1. Compute G = Z Z' / σ² (n×n GRM, streaming, GCTA-style)
//   2. Top-k eigendecomposition of G via Spectra → U (n×k), λ (k)
//   3. Store U, λ, MAF, σ²
//   4. Test projection (V-based streaming):
//      W = U · diag(1/√λ)
//      For each SNP j: test_pcs += z_test_j · (z_train_j' · W) / σ²
//      Avoids forming n_test × n_train G_cross matrix.
// ============================================================
class GrmPcaDimReduction : public DimReduction {
public:
    DimReductionResult fit(const PlinkReader& reader,
                           const std::vector<int>& sample_idx,
                           int n_components,
                           int n_threads = 0) override;

    Eigen::MatrixXd transform(const PlinkReader& reader,
                              const std::vector<int>& new_idx,
                              const DimReductionResult& params,
                              const std::vector<int>& train_idx,
                              int n_threads = 0) const override;

    std::string name() const override { return "grm_pca"; }

    /// Set GRM block size. Pass 0 for auto-detect (based on L3 cache).
    void set_block_size(int block_size, int n_samples);

    /// Set LD-pruned SNP indices. Pass nullptr or empty to use all SNPs.
    /// The pointer must remain valid for the lifetime of fit()/transform() calls.
    void set_snp_idx(const std::vector<int>* snp_idx) { snp_idx_ = snp_idx; }

    /// Combined fit + transform in a single SNP pass.
    /// Computes GRM, G_cross, eigendecomposition, and projections together,
    /// reading SNPs only once instead of twice (separate fit + transform).
    struct FitTransformResult {
        DimReductionResult params;   // eigenvectors, eigenvalues, maf, etc.
        Eigen::MatrixXd train_pcs;   // n_train × k
        Eigen::MatrixXd test_pcs;    // n_test × k
    };

    FitTransformResult fit_and_transform(
        const PlinkReader& reader,
        const std::vector<int>& train_idx,
        const std::vector<int>& test_idx,
        int n_components,
        int n_threads = 0) const;

    /// Single-pass G_cross approach: reads SNPs once, computes G + G_cross,
    /// eigendecomposes G, then projects test via G_cross * U * diag(1/λ).
    /// Avoids the second SNP pass at the cost of holding G_cross in memory.
    FitTransformResult fit_and_transform_gcross(
        const PlinkReader& reader,
        const std::vector<int>& train_idx,
        const std::vector<int>& test_idx,
        int n_components,
        int n_threads = 0) const;

private:
    int grm_block_size_ = 0;  // 0 = auto-detect
    const std::vector<int>* snp_idx_ = nullptr;  // LD-pruned SNP indices (optional)

    static void eigendecompose(Eigen::MatrixXd& G, int n, int k,
                               Eigen::MatrixXd& U, Eigen::VectorXd& lambda);
    static void fill_maf_params(const Eigen::VectorXd& maf, DimReductionResult& result);
};

// ============================================================
// Standard PCA for small matrices (e.g., phenotype matrix)
// ============================================================
struct StandardPcaResult {
    Eigen::VectorXd centers;
    Eigen::VectorXd scales;
    Eigen::MatrixXd rotation;  // V (p × k)
};

/// Standard PCA via center → scale → SVD (for in-memory matrices)
StandardPcaResult standard_pca_fit(const Eigen::MatrixXd& X,
                                   int n_components,
                                   int n_threads = 0);

/// Apply standard PCA transform to new data
Eigen::MatrixXd standard_pca_transform(const Eigen::MatrixXd& X_new,
                                       const StandardPcaResult& params,
                                       int n_threads = 0);

// ============================================================
// Factory function
// ============================================================
std::unique_ptr<DimReduction> create_dim_reduction(const std::string& method);

} // namespace fastcv
