#include "fastcv/types.h"
#include "fastcv/simd_util.h"
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

namespace fastcv {

// ============================================================
// Destructor / Move
// ============================================================
PlinkReader::~PlinkReader() { close(); }

PlinkReader::PlinkReader(PlinkReader&& o) noexcept
    : bed_fd_(o.bed_fd_), bed_file_size_(o.bed_file_size_),
      bed_data_offset_(o.bed_data_offset_),
      n_samples_(o.n_samples_), n_snps_(o.n_snps_),
      bed_file_path_(std::move(o.bed_file_path_)),
      sample_ids_(std::move(o.sample_ids_)),
      snp_info_(std::move(o.snp_info_))
{
    o.bed_fd_ = -1;
    o.bed_file_size_ = 0;
    o.bed_data_offset_ = 0;
    o.n_samples_ = 0;
    o.n_snps_ = 0;
}

PlinkReader& PlinkReader::operator=(PlinkReader&& o) noexcept {
    if (this != &o) {
        close();
        bed_fd_ = o.bed_fd_;
        bed_file_size_ = o.bed_file_size_;
        bed_data_offset_ = o.bed_data_offset_;
        n_samples_ = o.n_samples_;
        n_snps_ = o.n_snps_;
        bed_file_path_ = std::move(o.bed_file_path_);
        sample_ids_ = std::move(o.sample_ids_);
        snp_info_ = std::move(o.snp_info_);
        o.bed_fd_ = -1;
        o.bed_file_size_ = 0;
        o.bed_data_offset_ = 0;
        o.n_samples_ = 0;
        o.n_snps_ = 0;
    }
    return *this;
}

// ============================================================
// open()
// ============================================================
void PlinkReader::open(const std::string& bed_file) {
    close();

    std::string bim_file = bed_file;
    std::string fam_file = bed_file;
    auto replace_ext = [](std::string& s, const std::string& ext) {
        auto pos = s.rfind('.');
        if (pos != std::string::npos) s = s.substr(0, pos) + ext;
        else s += ext;
    };
    replace_ext(bim_file, ".bim");
    replace_ext(fam_file, ".fam");

    read_fam_(fam_file);
    read_bim_(bim_file);

    bed_fd_ = ::open(bed_file.c_str(), O_RDONLY);
    if (bed_fd_ < 0) {
        throw std::runtime_error("Cannot open BED file: " + bed_file);
    }

    struct stat st;
    if (fstat(bed_fd_, &st) < 0) {
        ::close(bed_fd_);
        bed_fd_ = -1;
        throw std::runtime_error("Cannot stat BED file: " + bed_file);
    }
    bed_file_size_ = static_cast<long>(st.st_size);

    if (bed_file_size_ < 3) {
        ::close(bed_fd_);
        bed_fd_ = -1;
        throw std::runtime_error("Invalid BED file: too short to read header");
    }
    int bps = (n_samples_ + 3) / 4;
    long expected_min = 3 + static_cast<long>(n_snps_) * bps;
    if (bed_file_size_ < expected_min) {
        ::close(bed_fd_);
        bed_fd_ = -1;
        throw std::runtime_error("BED file too small: expected at least " +
                                 std::to_string(expected_min) + " bytes, got " +
                                 std::to_string(bed_file_size_));
    }

    validate_bed_();

    bed_data_offset_ = 3;
    bed_file_path_ = bed_file;
}

void PlinkReader::close() {
    if (bed_fd_ >= 0) {
        ::close(bed_fd_);
        bed_fd_ = -1;
    }
    bed_file_size_ = 0;
    bed_data_offset_ = 0;
}

// ============================================================
// read_snp_raw(): pread one SNP's raw bytes (thread-safe)
// ============================================================
void PlinkReader::read_snp_raw(int snp_idx, unsigned char* buf) const {
    int bps = (n_samples_ + 3) / 4;
    long offset = bed_data_offset_ + static_cast<long>(snp_idx) * bps;
    ssize_t got = ::pread(bed_fd_, buf, bps, offset);
    if (got != bps) {
        throw std::runtime_error("pread failed for SNP " + std::to_string(snp_idx));
    }
}

// ============================================================
// Private: read FAM file
// ============================================================
void PlinkReader::read_fam_(const std::string& fam_file) {
    std::ifstream fin(fam_file);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open FAM file: " + fam_file);
    }
    sample_ids_.clear();
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string fid, iid;
        iss >> fid >> iid;
        if (fid == iid) sample_ids_.push_back(fid);
        else sample_ids_.push_back(fid + "_" + iid);
    }
    n_samples_ = static_cast<int>(sample_ids_.size());
}

// ============================================================
// Private: read BIM file
// ============================================================
void PlinkReader::read_bim_(const std::string& bim_file) {
    std::ifstream fin(bim_file);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open BIM file: " + bim_file);
    }
    snp_info_.clear();
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        SnpInfo info;
        iss >> info.chr >> info.snp_id >> info.cm >> info.pos >> info.a1 >> info.a2;
        snp_info_.push_back(info);
    }
    n_snps_ = static_cast<int>(snp_info_.size());
}

// ============================================================
// Private: validate BED file header via pread
// ============================================================
void PlinkReader::validate_bed_() {
    unsigned char hdr[3];
    ssize_t got = ::pread(bed_fd_, hdr, 3, 0);
    if (got != 3) {
        throw std::runtime_error("Cannot read BED file header");
    }
    if (hdr[0] != 0x6c || hdr[1] != 0x1b) {
        throw std::runtime_error("Invalid BED file: wrong magic number");
    }
    if (hdr[2] != 0x01) {
        throw std::runtime_error("Only SNP-major mode (0x01) is supported, got: " +
                                 std::to_string(hdr[2]));
    }
}

// ============================================================
// Private: decode a single SNP from BED via pread
// ============================================================
void PlinkReader::decode_snp_(int snp_idx, const std::vector<int>& sample_idx,
                              Eigen::VectorXi& geno) const {
    int n = static_cast<int>(sample_idx.size());
    geno.resize(n);

    int bps = (n_samples_ + 3) / 4;
    std::vector<unsigned char> snp_buf(bps);
    read_snp_raw(snp_idx, snp_buf.data());

    for (int i = 0; i < n; ++i) {
        int si = sample_idx[i];
        int byte_idx = si / 4;
        int bit_offset = (si % 4) * 2;
        int bits = (snp_buf[byte_idx] >> bit_offset) & 0x03;
        switch (bits) {
            case 0: geno(i) = 0;  break;
            case 1: geno(i) = -9; break;
            case 2: geno(i) = 1;  break;
            case 3: geno(i) = 2;  break;
        }
    }
}

// ============================================================
// read_snp(): read one SNP, raw genotypes (not centered/scaled)
// ============================================================
Eigen::VectorXi PlinkReader::read_snp(int snp_idx,
                                      const std::vector<int>& sample_idx) const {
    Eigen::VectorXi geno;
    decode_snp_(snp_idx, sample_idx, geno);
    return geno;
}

// ============================================================
// read_snp_block(): read a block of SNPs via pread, center/scale
// ============================================================
void PlinkReader::read_snp_block(int snp_start, int snp_end,
                                 const std::vector<int>& sample_idx,
                                 Eigen::MatrixXd& Z_block,
                                 Eigen::VectorXd& maf_block,
                                 const Eigen::VectorXd* external_maf,
                                 bool center_only) const {
    int block_size = snp_end - snp_start;
    int n = static_cast<int>(sample_idx.size());
    int bps = (n_samples_ + 3) / 4;

    // Read the entire block in one pread call
    long block_bytes = static_cast<long>(block_size) * bps;
    long block_offset = bed_data_offset_ + static_cast<long>(snp_start) * bps;
    std::vector<unsigned char> raw_buf(block_bytes);
    ssize_t got = ::pread(bed_fd_, raw_buf.data(), block_bytes, block_offset);
    if (got != block_bytes) {
        throw std::runtime_error("pread failed for SNP block " +
                                 std::to_string(snp_start) + "-" + std::to_string(snp_end));
    }
    const unsigned char* block_start = raw_buf.data();

    // Pre-compute sample byte indices and bit offsets
    std::vector<std::pair<int,int>> sample_loc(n);
    for (int i = 0; i < n; ++i) {
        int si = sample_idx[i];
        sample_loc[i] = {si / 4, (si % 4) * 2};
    }

    std::vector<int> geno_flat(n * block_size);
    maf_block.resize(block_size);

    if (external_maf) {
        maf_block = *external_maf;
    }

    for (int j = 0; j < block_size; ++j) {
        const unsigned char* snp_bytes = block_start + j * bps;
        int* geno_col = geno_flat.data() + j * n;
        int sum = 0;
        int valid = 0;
        for (int i = 0; i < n; ++i) {
            auto [byte_idx, bit_offset] = sample_loc[i];
            int g = kGenoLut[snp_bytes[byte_idx]][bit_offset >> 1];
            geno_col[i] = g;
            if (g >= 0) { sum += g; ++valid; }
        }
        if (!external_maf) {
            double freq = (valid > 0) ? sum / (2.0 * valid) : 0.5;
            maf_block(j) = freq;
        }
    }

    // Center (and optionally scale) into Z_block
    Z_block.resize(n, block_size);
    if (center_only) {
#ifdef __AVX2__
        for (int j = 0; j < block_size; ++j)
            center_scale_avx2(geno_flat.data() + j * n, n, &Z_block(0, j),
                              2.0 * maf_block(j), 1.0);
#else
        for (int j = 0; j < block_size; ++j) {
            double center = 2.0 * maf_block(j);
            const int* geno_col = geno_flat.data() + j * n;
            double* z_col = &Z_block(0, j);
            for (int i = 0; i < n; ++i) {
                int g = geno_col[i];
                z_col[i] = (g >= 0) ? (g - center) : 0.0;
            }
        }
#endif
    } else {
#ifdef __AVX2__
        for (int j = 0; j < block_size; ++j) {
            double p = maf_block(j);
            double scale = std::sqrt(2.0 * p * (1.0 - p));
            if (scale < 1e-10) scale = 1.0;
            const int* geno_col = geno_flat.data() + j * n;
            double* z_col = &Z_block(0, j);
            center_scale_avx2(geno_col, n, z_col, 2.0 * p, 1.0 / scale);
        }
#else
        for (int j = 0; j < block_size; ++j) {
            double p = maf_block(j);
            double scale = std::sqrt(2.0 * p * (1.0 - p));
            if (scale < 1e-10) scale = 1.0;
            double inv_scale = 1.0 / scale;
            double center = 2.0 * p;
            const int* geno_col = geno_flat.data() + j * n;
            double* z_col = &Z_block(0, j);
            for (int i = 0; i < n; ++i) {
                int g = geno_col[i];
                z_col[i] = (g >= 0) ? (g - center) * inv_scale : 0.0;
            }
        }
#endif
    }
}

// ============================================================
// read_snp_block_dual(): read one SNP block for two sample sets
// via single pread. MAF from sample_idx_a used for both.
// ============================================================
void PlinkReader::read_snp_block_dual(int snp_start, int snp_end,
                                      const std::vector<int>& sample_idx_a,
                                      const std::vector<int>& sample_idx_b,
                                      Eigen::MatrixXd& Z_a,
                                      Eigen::MatrixXd& Z_b,
                                      Eigen::VectorXd& maf_block,
                                      bool center_only) const {
    int block_size = snp_end - snp_start;
    int n_a = static_cast<int>(sample_idx_a.size());
    int n_b = static_cast<int>(sample_idx_b.size());
    int bps = (n_samples_ + 3) / 4;

    // Read the entire block in one pread call
    long block_bytes = static_cast<long>(block_size) * bps;
    long block_offset = bed_data_offset_ + static_cast<long>(snp_start) * bps;
    std::vector<unsigned char> raw_buf(block_bytes);
    ssize_t got = ::pread(bed_fd_, raw_buf.data(), block_bytes, block_offset);
    if (got != block_bytes) {
        throw std::runtime_error("pread failed for SNP block (dual) " +
                                 std::to_string(snp_start) + "-" + std::to_string(snp_end));
    }
    const unsigned char* block_start = raw_buf.data();

    std::vector<std::pair<int,int>> loc_a(n_a), loc_b(n_b);
    for (int i = 0; i < n_a; ++i) {
        int si = sample_idx_a[i];
        loc_a[i] = {si / 4, (si % 4) * 2};
    }
    for (int i = 0; i < n_b; ++i) {
        int si = sample_idx_b[i];
        loc_b[i] = {si / 4, (si % 4) * 2};
    }

    std::vector<int> geno_a(n_a * block_size);
    std::vector<int> geno_b(n_b * block_size);
    maf_block.resize(block_size);

    for (int j = 0; j < block_size; ++j) {
        const unsigned char* snp_bytes = block_start + j * bps;
        int* col_a = geno_a.data() + j * n_a;
        int* col_b = geno_b.data() + j * n_b;

        int sum = 0, valid = 0;
        for (int i = 0; i < n_a; ++i) {
            auto [byte_idx, bit_offset] = loc_a[i];
            int g = kGenoLut[snp_bytes[byte_idx]][bit_offset >> 1];
            col_a[i] = g;
            if (g >= 0) { sum += g; ++valid; }
        }
        for (int i = 0; i < n_b; ++i) {
            auto [byte_idx, bit_offset] = loc_b[i];
            col_b[i] = kGenoLut[snp_bytes[byte_idx]][bit_offset >> 1];
        }

        double freq = (valid > 0) ? sum / (2.0 * valid) : 0.5;
        maf_block(j) = freq;
    }

    Z_a.resize(n_a, block_size);
    Z_b.resize(n_b, block_size);

    if (center_only) {
#ifdef __AVX2__
        for (int j = 0; j < block_size; ++j) {
            center_scale_avx2(geno_a.data() + j * n_a, n_a, &Z_a(0, j),
                              2.0 * maf_block(j), 1.0);
            center_scale_avx2(geno_b.data() + j * n_b, n_b, &Z_b(0, j),
                              2.0 * maf_block(j), 1.0);
        }
#else
        for (int j = 0; j < block_size; ++j) {
            double center = 2.0 * maf_block(j);
            const int* col_a = geno_a.data() + j * n_a;
            double* z_a = &Z_a(0, j);
            for (int i = 0; i < n_a; ++i)
                z_a[i] = (col_a[i] >= 0) ? (col_a[i] - center) : 0.0;
            const int* col_b = geno_b.data() + j * n_b;
            double* z_b = &Z_b(0, j);
            for (int i = 0; i < n_b; ++i)
                z_b[i] = (col_b[i] >= 0) ? (col_b[i] - center) : 0.0;
        }
#endif
    } else {
#ifdef __AVX2__
        for (int j = 0; j < block_size; ++j) {
            double p = maf_block(j);
            double scale = std::sqrt(2.0 * p * (1.0 - p));
            if (scale < 1e-10) scale = 1.0;

            center_scale_avx2(geno_a.data() + j * n_a, n_a, &Z_a(0, j), 2.0 * p, 1.0 / scale);
            center_scale_avx2(geno_b.data() + j * n_b, n_b, &Z_b(0, j), 2.0 * p, 1.0 / scale);
        }
#else
        for (int j = 0; j < block_size; ++j) {
            double p = maf_block(j);
            double scale = std::sqrt(2.0 * p * (1.0 - p));
            if (scale < 1e-10) scale = 1.0;
            double inv_scale = 1.0 / scale;
            double center = 2.0 * p;

            const int* col_a = geno_a.data() + j * n_a;
            double* z_a = &Z_a(0, j);
            for (int i = 0; i < n_a; ++i)
                z_a[i] = (col_a[i] >= 0) ? (col_a[i] - center) * inv_scale : 0.0;

            const int* col_b = geno_b.data() + j * n_b;
            double* z_b = &Z_b(0, j);
            for (int i = 0; i < n_b; ++i)
                z_b[i] = (col_b[i] >= 0) ? (col_b[i] - center) * inv_scale : 0.0;
        }
#endif
    }
}

// ============================================================
// read_snp_block_dual_buf(): pre-allocated buffer version
// ============================================================
void PlinkReader::read_snp_block_dual_buf(
    int snp_start, int snp_end,
    const std::vector<int>& sample_idx_a,
    const std::vector<int>& sample_idx_b,
    Eigen::MatrixXd& Z_a,
    Eigen::MatrixXd& Z_b,
    Eigen::VectorXd& maf_block,
    DualReadBuf& buf,
    bool center_only) const
{
    int block_size = snp_end - snp_start;
    int n_a = static_cast<int>(sample_idx_a.size());
    int n_b = static_cast<int>(sample_idx_b.size());
    int bps = (n_samples_ + 3) / 4;

    long block_bytes = static_cast<long>(block_size) * bps;
    long block_offset = bed_data_offset_ + static_cast<long>(snp_start) * bps;
    ssize_t got = ::pread(bed_fd_, buf.raw_buf.data(), block_bytes, block_offset);
    if (got != block_bytes) {
        throw std::runtime_error("pread failed for SNP block (dual buf) " +
                                 std::to_string(snp_start) + "-" + std::to_string(snp_end));
    }

    for (int j = 0; j < block_size; ++j) {
        const unsigned char* snp_bytes = buf.raw_buf.data() + static_cast<long>(j) * bps;
        int* col_a = buf.geno_a.data() + static_cast<long>(j) * n_a;
        int* col_b = buf.geno_b.data() + static_cast<long>(j) * n_b;

        int sum = 0, valid = 0;
        for (int i = 0; i < n_a; ++i) {
            auto [byte_idx, bit_offset] = buf.loc_a[i];
            int g = kGenoLut[snp_bytes[byte_idx]][bit_offset >> 1];
            col_a[i] = g;
            if (g >= 0) { sum += g; ++valid; }
        }
        for (int i = 0; i < n_b; ++i) {
            auto [byte_idx, bit_offset] = buf.loc_b[i];
            col_b[i] = kGenoLut[snp_bytes[byte_idx]][bit_offset >> 1];
        }

        double freq = (valid > 0) ? sum / (2.0 * valid) : 0.5;
        maf_block(j) = freq;
    }

    if (center_only) {
#ifdef __AVX2__
        for (int j = 0; j < block_size; ++j) {
            center_scale_avx2(buf.geno_a.data() + static_cast<long>(j) * n_a, n_a, &Z_a(0, j),
                              2.0 * maf_block(j), 1.0);
            center_scale_avx2(buf.geno_b.data() + static_cast<long>(j) * n_b, n_b, &Z_b(0, j),
                              2.0 * maf_block(j), 1.0);
        }
#else
        for (int j = 0; j < block_size; ++j) {
            double center = 2.0 * maf_block(j);
            const int* ca = buf.geno_a.data() + static_cast<long>(j) * n_a;
            double* za = &Z_a(0, j);
            for (int i = 0; i < n_a; ++i)
                za[i] = (ca[i] >= 0) ? (ca[i] - center) : 0.0;
            const int* cb = buf.geno_b.data() + static_cast<long>(j) * n_b;
            double* zb = &Z_b(0, j);
            for (int i = 0; i < n_b; ++i)
                zb[i] = (cb[i] >= 0) ? (cb[i] - center) : 0.0;
        }
#endif
    } else {
#ifdef __AVX2__
        for (int j = 0; j < block_size; ++j) {
            double p = maf_block(j);
            double scale = std::sqrt(2.0 * p * (1.0 - p));
            if (scale < 1e-10) scale = 1.0;
            center_scale_avx2(buf.geno_a.data() + static_cast<long>(j) * n_a, n_a, &Z_a(0, j), 2.0 * p, 1.0 / scale);
            center_scale_avx2(buf.geno_b.data() + static_cast<long>(j) * n_b, n_b, &Z_b(0, j), 2.0 * p, 1.0 / scale);
        }
#else
        for (int j = 0; j < block_size; ++j) {
            double p = maf_block(j);
            double scale = std::sqrt(2.0 * p * (1.0 - p));
            if (scale < 1e-10) scale = 1.0;
            double inv_scale = 1.0 / scale;
            double center = 2.0 * p;
            const int* ca = buf.geno_a.data() + static_cast<long>(j) * n_a;
            double* za = &Z_a(0, j);
            for (int i = 0; i < n_a; ++i)
                za[i] = (ca[i] >= 0) ? (ca[i] - center) * inv_scale : 0.0;
            const int* cb = buf.geno_b.data() + static_cast<long>(j) * n_b;
            double* zb = &Z_b(0, j);
            for (int i = 0; i < n_b; ++i)
                zb[i] = (cb[i] >= 0) ? (cb[i] - center) * inv_scale : 0.0;
        }
#endif
    }
}

// ============================================================
// compute_maf(): compute MAF for all SNPs
// ============================================================
Eigen::VectorXd PlinkReader::compute_maf(const std::vector<int>& sample_idx) const {
    int n = static_cast<int>(sample_idx.size());
    Eigen::VectorXd maf(n_snps_);

    const int block_size = 1000;
    Eigen::VectorXi snp_buf(n);

    for (int s = 0; s < n_snps_; s += block_size) {
        int end = std::min(s + block_size, n_snps_);
        for (int j = s; j < end; ++j) {
            decode_snp_(j, sample_idx, snp_buf);
            double sum = 0.0;
            int valid = 0;
            for (int i = 0; i < n; ++i) {
                int g = snp_buf(i);
                if (g >= 0) { sum += g; ++valid; }
            }
            double freq = (valid > 0) ? sum / (2.0 * valid) : 0.5;
            maf(j) = freq;
        }
    }
    return maf;
}

} // namespace fastcv
