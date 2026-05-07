#pragma once

#include <cstdint>
#include <cstddef>
#include <array>

namespace fastcv {

// ============================================================
// Compile-time byte → 4-genotype lookup table
// BED encoding: 00→0(hom major), 01→-9(missing), 10→1(het), 11→2(hom minor)
// Indexed as kGenoLut[byte_value][geno_index], where geno_index = bit_offset / 2
// ============================================================
namespace detail {

constexpr int8_t decode_geno(int bits) {
    switch (bits) {
        case 0: return 0;
        case 1: return -9;
        case 2: return 1;
        case 3: return 2;
        default: return 0;
    }
}

constexpr std::array<std::array<int8_t, 4>, 256> make_geno_lut() {
    std::array<std::array<int8_t, 4>, 256> lut{};
    for (int b = 0; b < 256; ++b)
        for (int k = 0; k < 4; ++k)
            lut[b][k] = decode_geno((b >> (k * 2)) & 0x03);
    return lut;
}

} // namespace detail

constexpr auto kGenoLut = detail::make_geno_lut();

// ============================================================
// AVX2-accelerated center/scale
// Converts n int32 genotypes to double, applying (g - center) * inv_scale
// for g >= 0, else 0. Missing values (g = -9) are zeroed.
// ============================================================
#ifdef __AVX2__
#include <immintrin.h>

static inline void center_scale_avx2(const int* geno, int n,
                                       double* z_out,
                                       double center, double inv_scale) {
    const __m256d center_v = _mm256_set1_pd(center);
    const __m256d inv_scale_v = _mm256_set1_pd(inv_scale);
    const __m256d zero_v = _mm256_setzero_pd();

    int i = 0;
    for (; i + 3 < n; i += 4) {
        __m128i gi = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&geno[i]));
        __m256d gd = _mm256_cvtepi32_pd(gi);
        // Branchless mask: all-ones where g >= 0, all-zeros where g = -9
        __m256d valid = _mm256_cmp_pd(gd, zero_v, _CMP_GE_OQ);
        __m256d scaled = _mm256_mul_pd(_mm256_sub_pd(gd, center_v), inv_scale_v);
        _mm256_storeu_pd(&z_out[i], _mm256_and_pd(scaled, valid));
    }
    // Scalar remainder
    for (; i < n; ++i)
        z_out[i] = (geno[i] >= 0) ? (geno[i] - center) * inv_scale : 0.0;
}

#endif // __AVX2__

} // namespace fastcv
