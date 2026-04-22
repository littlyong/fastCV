# fastCV

**fastCV** is a high-performance C++17 tool for cross-validation data preparation in genomic prediction. It performs fold-wise principal component (PC) correction using genotype and/or phenotype covariates, producing corrected residuals for downstream genomic analyses (e.g., genomic BLUP, Bayesian models).

Designed for large-scale datasets (millions of SNPs, tens of thousands of samples), fastCV uses streaming PLINK I/O (`pread`, no mmap) and GRM-based dimensionality reduction, keeping peak memory at O(n²) instead of O(nm).

## Features

- **Memory-efficient streaming I/O** — Reads PLINK `.bed` via `pread()` (thread-safe, no mmap); never loads the full n×m genotype matrix
- **GRM-PCA dimensionality reduction** — GCTA-style GRM with block-wise `dsyrk`, top-k eigenvectors via Spectra
- **Leak-free fold-wise correction** — PCs computed exclusively from training samples; test samples projected via V-based streaming (no G_cross matrix)
- **OpenMP + MKL parallelism** — Multi-threaded GRM, eigendecomposition, and fold-level parallelism with NUMA-aware thread binding
- **Multi-trait support** — Correct single or multiple traits: `"1"`, `"all"`, `"2-10"`, `"1,3,5-7"`
- **Geno + Pheno correction** — Optional phenotype PCA correction alongside genotype PCs
- **LD pruning** — Built-in LD-based SNP pruning before GRM computation
- **Nested cross-validation** — Automatic hyperparameter tuning via nested folds
- **External PCA backends** — Integrate with HIBLUP or PLINK for PCA
- **Memory prediction** — `fastcv validate` estimates peak RAM from L3 cache and dataset size
- **AVX2-accelerated decoding** — Byte lookup table + SIMD center/scale for PLINK 2-bit genotype decoding

## Algorithm

### GRM-PCA (default path, 2-pass V-projection)

```
Pass 1 (train only):
  G = Z_train · Z_train' / Σ(2p(1-p))     GCTA-style GRM (streaming dsyrk)
  G = U Λ U'                               top-k eigendecomposition (Spectra)

Pass 2 (test projection, no G_cross matrix):
  W = U · diag(1/λ)                        (n_train × k)
  V_block = Z_train' · W                    per-block (bs × k)
  test_pcs += Z_test · V_block             streaming accumulation
  test_pcs /= σ²                            normalize
```

**Memory:** Peak = n_train² × 8 bytes (GRM), released before V-proj.

### G_cross single-pass (`use_gcross_proj: true`)

Computes G and G_cross simultaneously in one SNP pass. Higher peak memory (G + G_cross coexist) but avoids the second SNP scan.

## Installation

### Prerequisites

| Requirement           | Version                             |
| --------------------- | ----------------------------------- |
| CMake                 | >= 3.16                             |
| C++ compiler          | C++17 (GCC >= 8, Clang >= 7)        |
| OpenMP                | Optional (recommended)              |
| Intel MKL or OpenBLAS | Optional (Eigen fallback if absent) |

### Build

```bash
git clone https://github.com/<user>/fastCV.git
cd fastCV
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

The binary is `build/fastcv`.

### Build Options

| CMake Option          | Default | Description                              |
| --------------------- | ------- | ---------------------------------------- |
| `FASTCV_BUILD_CLI`    | ON      | Build CLI executable                     |
| `FASTCV_BUILD_SHARED` | OFF     | Build shared library (static by default) |

### BLAS Backend Selection

fastCV auto-detects the best available BLAS:

1. **Intel MKL** — highest performance (detected via `mkl.h`)
2. **OpenBLAS** — vendored source build in `third_party/`
3. **Eigen** — pure C++ fallback (no external dependency)

To use MKL, place headers and static libraries under `third_party/mkl_env/`.

## Quick Start

### 1. Generate configuration template

```bash
./fastcv config -o my_config.json
```

### 2. Edit configuration

Open `my_config.json` and set your data paths:

```json
{
  "data": {
    "genotype": {
      "bed_file": "/path/to/genotype.bed"
    },
    "phenotype": {
      "pheno_file": "/path/to/phenotype.txt",
      "id_col_name": "sample",
      "trait_idx": "1"
    },
    "covariate": {
      "cov_file": "/path/to/covariate.txt",
      "id_col_name": "sample",
      "cov_name": "all"
    }
  },
  "cross_validation": {
    "cv": { "n_folds": 5, "seed": 42 },
    "correction": { "mode": "geno_only", "n_geno_pcs": 5 },
    "dim_reduction": { "method": "pca" }
  },
  "other": { "n_threads": 0, "output_dir": "cv_output" }
}
```

### 3. Validate and predict memory

```bash
./fastcv validate -c my_config.json
```

Output:

```
  Configuration is valid!

=== Dataset ===
  Bed file:       /path/to/genotype.bed
  Samples:        3,284
  SNPs:           45,678
  Pheno file:     /path/to/phenotype.txt

=== Cross-Validation ===
  Folds:          5
  Seed:           42
  Correction:     geno_only
  Genotype PCs:   5

=== Hardware ===
  L3 cache/NUMA:  256 MiB
  Threads:        32 (auto)
  Base block sz:  4096

=== Fold Split ===
  Total samples:  3,284
  Train/fold:     2,627
  Test/fold:      657

=== Memory Prediction ===
  Phase               Parallel    BlockSz   Peak RAM
  ------------------------------------------------------
  GRM                 YES (32t)   4096      2.15 GB
  V-proj (2-pass)     YES (32t)   1528      0.52 GB
  ------------------------------------------------------
  TOTAL PEAK                                2.15 GB
```

### 4. Run cross-validation

```bash
./fastcv prepare -c my_config.json
```

## CLI Reference

```
Usage: fastcv <command> [options]

Commands:
  prepare    Run cross-validation with fold-wise PC correction
  validate   Validate config and predict peak memory usage
  config     Generate a configuration template with all options
  help       Show help message

Options:
  -c, --config <file>     Configuration file (for prepare / validate)
  -o, --output <file>     Output file (for config)
```

## Configuration Reference

`fastcv config -o my_config.json` generates a template with inline `_comment` descriptions. All fields:

### Data Paths

| JSON Path                    | Type   | Default    | Description                                            |
| ---------------------------- | ------ | ---------- | ------------------------------------------------------ |
| `data.genotype.bed_file`     | string | —          | PLINK `.bed` file (required)                           |
| `data.phenotype.pheno_file`  | string | —          | Phenotype file (required)                              |
| `data.phenotype.id_col_name` | string | `"sample"` | Sample ID column in phenotype file                     |
| `data.phenotype.trait_name`  | string | `""`       | Trait column name (overrides `trait_idx`)              |
| `data.phenotype.trait_idx`   | string | `"1"`      | Trait selection: `"1"`, `"all"`, `"2-10"`, `"1,3,5-7"` |
| `data.covariate.cov_file`    | string | `""`       | Covariate file (optional)                              |
| `data.covariate.id_col_name` | string | `"sample"` | Sample ID column in covariate file                     |
| `data.covariate.cov_name`    | string | `"all"`    | Covariate selection: `"all"` or comma-separated names  |

### Cross-Validation

| JSON Path                                 | Type   | Default       | Description                               |
| ----------------------------------------- | ------ | ------------- | ----------------------------------------- |
| `cross_validation.cv.n_folds`             | int    | 5             | Number of CV folds (≥ 2)                  |
| `cross_validation.cv.seed`                | int    | 42            | Random seed                               |
| `cross_validation.correction.mode`        | string | `"geno_only"` | `"geno_only"` or `"geno_pheno"`           |
| `cross_validation.correction.n_geno_pcs`  | int    | 5             | Number of genotype PCs (≥ 1)              |
| `cross_validation.correction.n_pheno_pcs` | int    | 10            | Number of phenotype PCs (geno_pheno mode) |
| `cross_validation.stratified_cov_name`    | string | `""`          | Covariate column for stratified folds     |

### Dimensionality Reduction

| JSON Path                                        | Type   | Default      | Description                                    |
| ------------------------------------------------ | ------ | ------------ | ---------------------------------------------- |
| `cross_validation.dim_reduction.method`          | string | `"pca"`      | `"pca"` (GRM-based) or `"grm_pca"` (legacy)    |
| `cross_validation.dim_reduction.pca_backend`     | string | `"internal"` | `"internal"`, `"hiblup"`, or `"plink"`         |
| `cross_validation.dim_reduction.pca_tool_path`   | string | `""`         | Path to external binary                        |
| `cross_validation.dim_reduction.grm_block_size`  | int    | 0            | SNPs per block (0 = auto-detect from L3)       |
| `cross_validation.dim_reduction.use_gcross_proj` | bool   | false        | Single-pass G_cross (higher RAM, one SNP scan) |

### LD Pruning

| JSON Path                                        | Type   | Default | Description                  |
| ------------------------------------------------ | ------ | ------- | ---------------------------- |
| `cross_validation.dim_reduction.ld_prune`        | bool   | false   | Enable LD pruning before GRM |
| `cross_validation.dim_reduction.ld_r2_threshold` | double | 0.2     | Maximum r² for LD pruning    |
| `cross_validation.dim_reduction.ld_window_size`  | int    | 50      | Sliding window size (SNPs)   |
| `cross_validation.dim_reduction.ld_maf_min`      | double | 0.01    | Minimum allele frequency     |

### Nested Cross-Validation

| JSON Path                                   | Type | Default | Description                                |
| ------------------------------------------- | ---- | ------- | ------------------------------------------ |
| `cross_validation.nested_cv.n_nested_folds` | int  | 0       | Nested CV folds (0 = disabled)             |
| `cross_validation.nested_cv.nested_seed`    | int  | 0       | Seed for nested folds (0 = same as `seed`) |

### Output & Performance

| JSON Path                                 | Type   | Default       | Description                      |
| ----------------------------------------- | ------ | ------------- | -------------------------------- |
| `cross_validation.save_details.cv`        | bool   | false         | Save MAF/beta/PCs per outer fold |
| `cross_validation.save_details.nested_cv` | bool   | false         | Save details per nested fold     |
| `other.n_threads`                         | int    | 0             | OpenMP threads (0 = auto-detect) |
| `other.output_dir`                        | string | `"cv_output"` | Output directory                 |
| `other.separator`                         | string | `" "`         | Output field separator           |

## Input File Formats

### PLINK Binary Files

Standard PLINK 1 `.bed` / `.bim` / `.fam` triple (SNP-major mode). All three files must share the same prefix.

### Phenotype File

Whitespace-delimited text file with an ID column and trait columns:

```
sample  height  weight  yield
ID001   165.3   72.1    4.5
ID002   152.8   58.4    3.8
ID003   NA      65.0    5.1
```

- Missing values: `NA` or empty
- ID column name matches `id_col_name` (default: `"sample"`)
- Sample IDs must match FAM file

### Covariate File (optional)

Same format as phenotype file, with numeric covariate columns:

```
sample  age  location
ID001   3.5  1
ID002   2.1  2
```

## Output Structure

```
cv_output/
├── config.json                    # Run configuration (echoed)
├── fold_assignments.csv           # sample_id, fold_id
├── fold_01/
│   ├── train_samples.txt
│   ├── test_samples.txt
│   ├── geno_pcs_train.csv         # Genotype PCs (training)
│   ├── geno_pcs_test.csv          # Genotype PCs (test, projected)
│   ├── pheno_pcs_train.csv        # Phenotype PCs (geno_pheno mode)
│   ├── pheno_pcs_test.csv
│   ├── y_train_residual.csv       # sample_id, y_original, y_residual
│   ├── y_test_residual.csv
│   └── correction_params.json     # beta, eigenvectors, eigenvalues, MAF
├── fold_02/
│   └── ...
```

## C++ API

### Full Pipeline

```cpp
#include <fastcv/fastcv.h>

fastcv::CvConfig config = fastcv::load_config("config.json");
fastcv::CvResult result = fastcv::prepare_cv_data(config);
fastcv::export_cv_result(result, config.output_dir);
```

### GRM Computation

```cpp
#include <fastcv/grm.h>

// Streaming GRM (GCTA-style)
Eigen::VectorXd maf;
Eigen::MatrixXd G = fastcv::compute_grm(reader, sample_idx, 1000, 8, &maf);

// Combined GRM + cross-GRM in single pass
auto [G, G_cross, maf] = fastcv::compute_grm_with_cross(reader, train_idx, test_idx, 1000, 8);
```

### PLINK Streaming Reader

```cpp
#include <fastcv/types.h>

fastcv::PlinkReader reader;
reader.open("genotype.bed");
// reader.n_samples(), reader.n_snps()

// Read raw genotypes (0/1/2/-9)
Eigen::VectorXi geno = reader.read_snp(snp_idx, sample_idx);

// Read block, center & scale
Eigen::MatrixXd Z_block;
Eigen::VectorXd maf_block;
reader.read_snp_block(0, 1000, sample_idx, Z_block, maf_block);
```

## Performance

### Benchmarks (n=2797, m=1.2M SNPs, k=5, 16 threads, Intel Xeon)

| Phase                                    | Time |
| ---------------------------------------- | ---- |
| GRM (parallel dsyrk)                     | ~30s |
| Eigendecomposition (Spectra, k=5)        | <1s  |
| V-proj (parallel, pre-allocated buffers) | ~3s  |
| PC correction + export (5 folds)         | <1s  |

### Memory (n=2797, m=1.2M, 5-fold, k=5, 16 threads)

| Phase                           | Peak RAM    |
| ------------------------------- | ----------- |
| GRM (per-thread G_local)        | ~640 MB     |
| V-proj (per-thread DualReadBuf) | ~835 MB     |
| **Total peak**                  | **~835 MB** |

## Project Structure

```
fastCV/
├── CMakeLists.txt
├── include/fastcv/
│   ├── fastcv.h                # Umbrella header
│   ├── types.h                 # Core types (CvConfig, PlinkReader, etc.)
│   ├── grm.h                   # GRM functions
│   ├── dim_reduction.h         # DimReduction interface
│   ├── pc_corrector.h          # PC correction
│   ├── cross_validation.h      # Fold assignment
│   ├── data_export.h           # CSV export
│   ├── ld_prune.h              # LD pruning
│   ├── simd_util.h             # AVX2 genotype decoding
│   ├── memory_estimate.h       # Memory prediction
│   ├── config.h                # JSON config loader
│   └── cli.h                   # CLI entry point
├── src/                        # Implementation
├── configs/                    # Example configurations
└── third_party/                # Vendored dependencies
```

## Dependencies

| Library                                           | Purpose                  | License                  |
| ------------------------------------------------- | ------------------------ | ------------------------ |
| [Eigen](https://eigen.tuxfamily.org/) 3.4+        | Linear algebra           | MPL2                     |
| [Spectra](https://spectralib.org/)                | Top-k eigendecomposition | MPL2                     |
| [nlohmann/json](https://github.com/nlohmann/json) | JSON parsing             | MIT                      |
| Intel MKL (optional)                              | Optimized BLAS/LAPACK    | Intel Simplified License |
| OpenBLAS (optional)                               | Optimized BLAS           | BSD-3                    |

## License

MIT
