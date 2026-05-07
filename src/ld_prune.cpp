#include "fastcv/ld_prune.h"
#include "fastcv/misc_util.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <chrono>

namespace fastcv {

std::vector<int> ld_prune_snps(const PlinkReader& reader,
                               const std::vector<int>& sample_idx,
                               const LdPruneParams& params) {
    int m = reader.n_snps();
    int n = static_cast<int>(sample_idx.size());
    const auto& snp_info = reader.snp_info();
    int bytes_per_snp = (reader.n_samples() + 3) / 4;

    std::vector<int> keep_idx;
    int removed_maf = 0;
    int removed_ld = 0;

    auto t0 = std::chrono::high_resolution_clock::now();

    // Process SNPs in blocks, keep a sliding window of recently kept SNPs
    // (centered/scaled genotype vectors) for LD checking.
    // Window resets at chromosome boundaries.
    const int block_size = 1000;
    std::vector<Eigen::VectorXd> kept_z;  // centered/scaled genotypes of kept SNPs in current chr window
    std::string prev_chr;

    // Pre-compute sample byte/bit lookups (shared across all blocks)
    std::vector<std::pair<int,int>> sample_loc(n);
    for (int i = 0; i < n; ++i) {
        int si = sample_idx[i];
        sample_loc[i] = {si / 4, (si % 4) * 2};
    }

    int total_blocks = (m + block_size - 1) / block_size;

    for (int s = 0; s < m; s += block_size) {
        int end = std::min(s + block_size, m);
        int bs = end - s;

        // Read block via pread into local buffer
        long block_bytes = static_cast<long>(bs) * bytes_per_snp;
        std::vector<unsigned char> raw_buf(block_bytes);
        // Read each SNP individually via read_snp_raw
        for (int j = 0; j < bs; ++j) {
            reader.read_snp_raw(s + j, raw_buf.data() + static_cast<long>(j) * bytes_per_snp);
        }
        const unsigned char* block_start = raw_buf.data();

        // Decode + center/scale
        for (int j = 0; j < bs; ++j) {
            int global_j = s + j;
            const unsigned char* snp_bytes = block_start + j * bytes_per_snp;

            // Decode genotypes + compute MAF in one pass
            std::vector<int> geno_col(n);
            int sum = 0, valid = 0;
            static const int decode_tbl[4] = {0, -9, 1, 2};
            for (int i = 0; i < n; ++i) {
                auto [byte_idx, bit_offset] = sample_loc[i];
                int g = decode_tbl[(snp_bytes[byte_idx] >> bit_offset) & 0x03];
                geno_col[i] = g;
                if (g >= 0) { sum += g; ++valid; }
            }
            double freq = (valid > 0) ? sum / (2.0 * valid) : 0.5;

            // MAF filter
            if (freq < params.maf_min || freq > (1.0 - params.maf_min)) {
                removed_maf++;
                continue;
            }

            // Chromosome boundary: reset window
            const std::string& chr = snp_info[global_j].chr;
            if (chr != prev_chr) {
                kept_z.clear();
                prev_chr = chr;
            }

            // Center/scale
            double scale = std::sqrt(2.0 * freq * (1.0 - freq));
            if (scale < 1e-10) continue;
            double inv_scale = 1.0 / scale;
            double center = 2.0 * freq;
            Eigen::VectorXd z(n);
            for (int i = 0; i < n; ++i) {
                z(i) = (geno_col[i] >= 0) ? (geno_col[i] - center) * inv_scale : 0.0;
            }

            // Check LD with all kept SNPs in the window
            bool prune = false;
            for (const auto& kz : kept_z) {
                double r = z.dot(kz) / static_cast<double>(n);
                if (r * r > params.r2_threshold) {
                    prune = true;
                    break;
                }
            }

            if (!prune) {
                keep_idx.push_back(global_j);
                kept_z.push_back(std::move(z));
                // Limit window size
                if (static_cast<int>(kept_z.size()) > params.window_size) {
                    kept_z.erase(kept_z.begin());
                }
            } else {
                removed_ld++;
            }
        }

        // Progress
        int cur_block = s / block_size + 1;
        detail::print_progress_bar(cur_block * 100 / total_blocks, "LD prune");
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "\n      [LD prune] " << keep_idx.size() << "/" << m
              << " SNPs kept (MAF removed: " << removed_maf
              << ", LD removed: " << removed_ld << ")"
              << " in " << std::chrono::duration<double>(t1 - t0).count() << "s\n"
              << std::flush;

    return keep_idx;
}

} // namespace fastcv
