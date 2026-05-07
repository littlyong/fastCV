#pragma once

#include "fastcv/types.h"
#include <vector>

namespace fastcv {

struct LdPruneParams {
    double r2_threshold = 0.2;    // maximum allowed r²
    int window_size = 50;         // SNPs to consider in sliding window
    double maf_min = 0.01;        // minimum allele frequency to keep
};

/// Compute LD-pruned SNP indices (greedy, chromosome-aware).
/// Returns sorted vector of SNP indices to KEEP.
/// MAF is computed from the given sample set; both MAF filter and
/// pairwise LD check are applied.
std::vector<int> ld_prune_snps(const PlinkReader& reader,
                               const std::vector<int>& sample_idx,
                               const LdPruneParams& params = {});

} // namespace fastcv
