#pragma once

#include <string>
#include <cstdint>

namespace fastcv {

struct MemoryEstimate {
    int n_samples = 0;
    int n_snps = 0;
    int n_train = 0;
    int n_test = 0;
    int n_omp = 0;
    int block_size = 0;       // base block size (auto or manual)
    int par_bs_grm = 0;       // parallel block size for GRM phase
    int par_bs_vproj = 0;     // parallel block size for V-proj phase
    long l3_per_numa = 0;     // detected L3 cache bytes per NUMA node

    bool use_gcross_proj = false;
    bool grm_parallel = false;
    bool vproj_parallel = false;

    long grm_phase_peak = 0;
    long vproj_phase_peak = 0;
    long total_peak = 0;
};

/// Detect L3 cache size (bytes) per NUMA node via lscpu or /proc/cpuinfo.
long detect_l3_per_numa();

/// Get effective thread count (0 = auto-detect, capped at 64).
int get_effective_threads(int n_threads);

/// Set both OpenMP and MKL thread count (avoids oversubscription).
void set_threads(int n_threads);

/// Auto-detect GRM block size from L3 cache and sample count.
int auto_detect_block_size(long l3_per_numa, int n_samples);

/// Compute peak memory estimate for a CV run.
MemoryEstimate estimate_memory(int n_samples, int n_snps, int n_folds,
                               int n_geno_pcs, int n_threads_cfg,
                               int grm_block_size_cfg, bool use_gcross_proj);

/// Format bytes as human-readable string ("1.5 GB", "256 MB").
std::string format_bytes(long bytes);

/// Print formatted memory prediction report to stdout.
void print_memory_report(const MemoryEstimate& est);

} // namespace fastcv
