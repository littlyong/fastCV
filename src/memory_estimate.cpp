#include "fastcv/memory_estimate.h"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>

#ifdef FASTCV_HAS_OPENMP
#include <omp.h>
#endif

#ifdef FASTCV_USE_MKL
#include <mkl.h>
#endif

namespace fastcv {

// ============================================================
// detect_l3_per_numa(): parse lscpu / /proc/cpuinfo
// ============================================================
long detect_l3_per_numa() {
    // Method 1: parse lscpu
    FILE* pipe = popen("lscpu 2>/dev/null", "r");
    if (pipe) {
        long total_l3 = 0;
        int numa_nodes = 1;
        char buf[256];
        while (fgets(buf, sizeof(buf), pipe)) {
            std::string line(buf);
            if (line.find("L3 cache:") != std::string::npos) {
                long val;
                if (sscanf(line.c_str(), " L3 cache: %ld", &val) == 1) {
                    if (line.find("KiB") != std::string::npos) total_l3 = val * 1024;
                    else if (line.find("MiB") != std::string::npos) total_l3 = val * 1024 * 1024;
                    else total_l3 = val;
                }
            }
            if (line.find("NUMA node(s):") != std::string::npos) {
                int n;
                if (sscanf(line.c_str(), " NUMA node(s): %d", &n) == 1)
                    numa_nodes = n;
            }
        }
        pclose(pipe);
        if (total_l3 > 0 && numa_nodes > 0)
            return total_l3 / numa_nodes;
    }

    // Method 2: parse /proc/cpuinfo
    std::ifstream cpuinfo("/proc/cpuinfo");
    if (cpuinfo.is_open()) {
        long total_l3 = 0;
        std::string line;
        while (std::getline(cpuinfo, line)) {
            if (line.find("cache size") != std::string::npos) {
                long val;
                if (sscanf(line.c_str(), " cache size : %ld", &val) == 1)
                    total_l3 = val * 1024;
            }
        }
        std::ifstream numa("/sys/devices/system/node/possible");
        if (numa.is_open()) {
            std::string nline;
            std::getline(numa, nline);
            auto dash = nline.find('-');
            if (dash != std::string::npos) {
                int nn = std::stoi(nline.substr(dash + 1)) + 1;
                if (total_l3 > 0) return total_l3 / nn;
            }
        }
    }
    return 0;
}

// ============================================================
// get_effective_threads()
// ============================================================
int get_effective_threads(int n_threads) {
    if (n_threads > 0) return n_threads;
#ifdef FASTCV_HAS_OPENMP
    return std::min(omp_get_max_threads(), 64);
#else
    return 1;
#endif
}

// ============================================================
// set_threads(): set both OpenMP and MKL thread count
// ============================================================
void set_threads(int n_threads) {
    int nth = get_effective_threads(n_threads);
#ifdef FASTCV_HAS_OPENMP
    omp_set_num_threads(nth);
#endif
#ifdef FASTCV_USE_MKL
    mkl_set_num_threads(nth);
#endif
}

// ============================================================
// auto_detect_block_size()
// ============================================================
int auto_detect_block_size(long l3_per_numa, int n_samples) {
    if (n_samples <= 0) return 2000;
    if (l3_per_numa <= 0) return 2000;
    long target = std::min(static_cast<long>(1024LL * 1024 * 1024), l3_per_numa * 10);
    long bs = target / static_cast<long>(n_samples * 8);
    return static_cast<int>(std::max(500L, std::min(50000L, bs)));
}

// ============================================================
// format_bytes()
// ============================================================
std::string format_bytes(long bytes) {
    const long GB = 1024L * 1024 * 1024;
    const long MB = 1024L * 1024;
    const long KB = 1024L;
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1);
    if (bytes >= GB)       oss << static_cast<double>(bytes) / GB << " GB";
    else if (bytes >= MB)  oss << static_cast<double>(bytes) / MB << " MB";
    else if (bytes >= KB)  oss << static_cast<double>(bytes) / KB << " KB";
    else                   oss << bytes << " B";
    return oss.str();
}

// ============================================================
// estimate_memory()
// ============================================================
MemoryEstimate estimate_memory(int n_samples, int n_snps, int n_folds,
                               int n_geno_pcs, int n_threads_cfg,
                               int grm_block_size_cfg, bool use_gcross_proj) {
    MemoryEstimate est;
    est.n_samples = n_samples;
    est.n_snps = n_snps;
    est.use_gcross_proj = use_gcross_proj;

    if (n_samples <= 0) return est;

    // Fold split
    est.n_train = n_samples * (n_folds - 1) / n_folds;
    est.n_test = n_samples - est.n_train;
    if (est.n_train <= 0) est.n_train = 1;

    // Hardware
    est.l3_per_numa = detect_l3_per_numa();
    est.n_omp = get_effective_threads(n_threads_cfg);

    // Block size
    if (grm_block_size_cfg > 0) {
        est.block_size = grm_block_size_cfg;
    } else {
        est.block_size = auto_detect_block_size(est.l3_per_numa, est.n_train);
    }

    int k = n_geno_pcs;
    int n = est.n_train;
    int nt = est.n_test;
    int nth = est.n_omp;
    int bps = (n_samples + 3) / 4;
    const long GB8 = 8L * 1024 * 1024 * 1024;

    // === GRM phase ===
    {
        long g_local_bytes = static_cast<long>(n) * n * 8;
        est.grm_parallel = (nth > 1) && (static_cast<long>(nth) * g_local_bytes < GB8);

        if (est.grm_parallel) {
            long max_bs = 200L * 1024 * 1024 / (static_cast<long>(n) * 12);
            est.par_bs_grm = std::max(1000, std::min(est.block_size, static_cast<int>(max_bs)));
            long per_thread_accum = g_local_bytes;
            long per_thread_temp = static_cast<long>(n) * est.par_bs_grm * 12;
            est.grm_phase_peak = static_cast<long>(nth) * (per_thread_accum + per_thread_temp);
        } else {
            est.par_bs_grm = est.block_size;
            long g_bytes = g_local_bytes;
            long z_temp = static_cast<long>(n) * est.block_size * 12;
            est.grm_phase_peak = g_bytes + z_temp;
        }
    }

    // === Second phase (V-proj or G_cross) ===
    if (use_gcross_proj) {
        // GRM + G_cross single pass
        long accum = static_cast<long>(n) * n * 8 + static_cast<long>(nt) * n * 8;
        bool par = (nth > 1) && (static_cast<long>(nth) * accum < GB8);

        if (par) {
            long max_bs = 200L * 1024 * 1024 / (static_cast<long>(n + nt) * 12);
            int par_bs = std::max(1000, std::min(est.block_size, static_cast<int>(max_bs)));
            long per_thread_temp = static_cast<long>(n + nt) * par_bs * 12;
            est.vproj_phase_peak = static_cast<long>(nth) * (accum + per_thread_temp);
        } else {
            long temp = static_cast<long>(n + nt) * est.block_size * 12;
            est.vproj_phase_peak = accum + temp;
        }
        est.vproj_parallel = par;
        est.par_bs_vproj = 0;  // N/A for gcross
    } else {
        // V-proj (2-pass)
        long per_snp_mem = static_cast<long>(bps)
                         + static_cast<long>(n + nt) * 12
                         + static_cast<long>(k) * 8 + 8;
        long max_bs = 50L * 1024 * 1024 / per_snp_mem;
        est.par_bs_vproj = std::max(1000, std::min(est.block_size, static_cast<int>(max_bs)));

        est.vproj_parallel = (nth > 1);

        if (est.vproj_parallel) {
            long per_thread = static_cast<long>(est.par_bs_vproj) * per_snp_mem
                            + static_cast<long>(nt) * k * 8;
            long w_matrix = static_cast<long>(n) * k * 8;
            est.vproj_phase_peak = static_cast<long>(nth) * per_thread + w_matrix;
        } else {
            long temp = static_cast<long>(n + nt) * est.block_size * 8;
            long w_matrix = static_cast<long>(n) * k * 8;
            est.vproj_phase_peak = temp + w_matrix;
        }
    }

    // === Total peak: max of all phases ===
    est.total_peak = std::max(est.grm_phase_peak, est.vproj_phase_peak);

    return est;
}

// ============================================================
// print_memory_report()
// ============================================================
void print_memory_report(const MemoryEstimate& est) {
    auto W = [](const std::string& s, int w) {
        return s.size() >= static_cast<size_t>(w) ? s : s + std::string(w - s.size(), ' ');
    };

    if (est.n_samples <= 0) {
        std::cout << "\n=== Memory Prediction ===\n"
                  << "  (Precise estimates require a valid BED file.)\n";
        return;
    }

    std::cout << "\n=== Hardware ===\n";
    if (est.l3_per_numa > 0)
        std::cout << "  L3 cache/NUMA:  " << format_bytes(est.l3_per_numa) << "\n";
    else
        std::cout << "  L3 cache/NUMA:  (not detected)\n";
    std::cout << "  Threads:        " << est.n_omp << "\n";
    std::cout << "  Base block sz:  " << est.block_size << "\n";

    std::cout << "\n=== Fold Split ===\n";
    std::cout << "  Total samples:  " << est.n_samples << "\n";
    std::cout << "  Train/fold:     " << est.n_train << "\n";
    std::cout << "  Test/fold:      " << est.n_test << "\n";

    std::cout << "\n=== Memory Prediction ===\n";
    std::cout << "  " << W("Phase", 20) << W("Parallel", 12)
              << W("BlockSz", 10) << W("Peak RAM", 12) << "\n";
    std::cout << "  " << std::string(54, '-') << "\n";

    // GRM row
    {
        std::string par = est.grm_parallel
                          ? ("YES (" + std::to_string(est.n_omp) + "t)")
                          : "NO";
        std::cout << "  " << W("GRM", 20) << W(par, 12)
                  << W(std::to_string(est.par_bs_grm), 10)
                  << W(format_bytes(est.grm_phase_peak), 12) << "\n";
    }

    // Second phase row
    if (est.use_gcross_proj) {
        std::string par = est.vproj_parallel
                          ? ("YES (" + std::to_string(est.n_omp) + "t)")
                          : "NO";
        std::cout << "  " << W("GRM+Gx (single)", 20) << W(par, 12)
                  << W("-", 10)
                  << W(format_bytes(est.vproj_phase_peak), 12) << "\n";
    } else {
        std::string par = est.vproj_parallel
                          ? ("YES (" + std::to_string(est.n_omp) + "t)")
                          : "NO";
        std::cout << "  " << W("V-proj (2-pass)", 20) << W(par, 12)
                  << W(std::to_string(est.par_bs_vproj), 10)
                  << W(format_bytes(est.vproj_phase_peak), 12) << "\n";
    }

    std::cout << "  " << std::string(54, '-') << "\n";
    std::cout << "  " << W("TOTAL PEAK", 20) << W("", 12) << W("", 10)
              << format_bytes(est.total_peak) << "\n";
}

} // namespace fastcv
