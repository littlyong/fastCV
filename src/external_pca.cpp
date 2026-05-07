#include "fastcv/external_pca.h"
#include "fastcv/grm.h"
#include "fastcv/misc_util.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <chrono>

namespace fastcv {

namespace {

static std::string find_tool(const std::string& name, const std::string& user_path) {
    if (!user_path.empty()) return user_path;
    // Try system PATH
    std::string cmd = "which " + name + " 2>/dev/null";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return "";
    char buf[256];
    std::string result;
    if (fgets(buf, sizeof(buf), pipe)) {
        result = buf;
        while (!result.empty() && (result.back() == '\n' || result.back() == '\r'))
            result.pop_back();
    }
    pclose(pipe);
    return result;
}

static void write_keep_file(const std::string& path,
                            const std::vector<std::string>& ids) {
    std::ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("Cannot write: " + path);
    for (const auto& id : ids) out << id << "\n";
}

// Parse hiblup .pc file: "id\tPC1\tPC2\t..."
static Eigen::MatrixXd parse_hiblup_pc(const std::string& path, int n_pcs) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Cannot read: " + path);
    std::string header;
    std::getline(fin, header);

    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string id;
        iss >> id;
        std::vector<double> vals(n_pcs);
        for (int k = 0; k < n_pcs; ++k) iss >> vals[k];
        rows.push_back(vals);
    }

    int n = static_cast<int>(rows.size());
    Eigen::MatrixXd pcs(n, n_pcs);
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n_pcs; ++k)
            pcs(i, k) = rows[i][k];
    return pcs;
}

// Parse hiblup .pcp file: extract "Standard deviation" row
static Eigen::VectorXd parse_hiblup_pcp(const std::string& path, int n_pcs) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Cannot read: " + path);

    std::string line;
    while (std::getline(fin, line)) {
        if (line.find("Standard deviation") != std::string::npos) {
            std::istringstream iss(line);
            std::string label;
            iss >> label >> label;  // skip "Standard deviation"
            Eigen::VectorXd sd(n_pcs);
            for (int k = 0; k < n_pcs; ++k) iss >> sd(k);
            return sd;
        }
    }
    throw std::runtime_error("Cannot find 'Standard deviation' in: " + path);
}

// Parse plink .eigenvec file: "FID IID PC1 PC2 ..."
static Eigen::MatrixXd parse_plink_eigenvec(const std::string& path, int n_pcs) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Cannot read: " + path);

    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string fid, iid;
        iss >> fid >> iid;
        std::vector<double> vals(n_pcs);
        for (int k = 0; k < n_pcs; ++k) iss >> vals[k];
        rows.push_back(vals);
    }

    int n = static_cast<int>(rows.size());
    Eigen::MatrixXd pcs(n, n_pcs);
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n_pcs; ++k)
            pcs(i, k) = rows[i][k];
    return pcs;
}

// Parse plink .eigenval file: one eigenvalue per line
static Eigen::VectorXd parse_plink_eigenval(const std::string& path, int n_pcs) {
    std::ifstream fin(path);
    if (!fin.is_open()) throw std::runtime_error("Cannot read: " + path);

    Eigen::VectorXd evals(n_pcs);
    for (int k = 0; k < n_pcs; ++k) {
        if (!(fin >> evals(k))) {
            throw std::runtime_error("Failed to read eigenvalue " + std::to_string(k+1) +
                                     " from: " + path);
        }
    }
    return evals;
}

} // anonymous namespace

// ============================================================
// run_external_pca()
// ============================================================
ExternalPcaResult run_external_pca(
    const std::string& tool,
    const std::string& tool_path,
    const std::string& bed_prefix,
    const std::vector<std::string>& train_ids,
    int n_pcs,
    int n_threads,
    const std::string& work_dir,
    int fold_idx)
{
    std::string binary = find_tool(tool, tool_path);
    if (binary.empty()) {
        throw std::runtime_error(tool + " not found. Install it or set pca_tool_path in config.");
    }

    std::string fold_tag = "fold_" + std::string(fold_idx < 10 ? "0" : "") + std::to_string(fold_idx);
    std::string keep_file = work_dir + "/" + fold_tag + "_train_keep.txt";
    std::string out_prefix = work_dir + "/" + fold_tag + "_pca";

    // Write --keep file
    write_keep_file(keep_file, train_ids);

    // Build and execute command
    std::string cmd;
    if (tool == "hiblup") {
        cmd = binary + " --pca --bfile " + bed_prefix +
              " --keep " + keep_file +
              " --npc " + std::to_string(n_pcs) +
              " --out " + out_prefix +
              " --threads " + std::to_string(n_threads) +
              " 2>&1";
    } else if (tool == "plink") {
        cmd = binary + " --bfile " + bed_prefix +
              " --keep " + keep_file +
              " --pca " + std::to_string(n_pcs) +
              " --out " + out_prefix +
              " --threads " + std::to_string(n_threads) +
              " 2>&1";
    } else {
        throw std::runtime_error("Unsupported external PCA tool: " + tool);
    }

    std::cout << "    [EXT] Running " << tool << " for " << fold_tag << "...\n" << std::flush;
    auto t0 = std::chrono::high_resolution_clock::now();

    int ret = std::system(cmd.c_str());
    if (ret != 0) {
        throw std::runtime_error(tool + " failed (exit code " + std::to_string(ret) +
                                 ") for " + fold_tag);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double ext_sec = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "    [EXT] " << tool << " done in " << ext_sec << "s\n" << std::flush;

    // Parse output
    ExternalPcaResult result;
    int n_train = static_cast<int>(train_ids.size());

    if (tool == "hiblup") {
        Eigen::MatrixXd P_unit = parse_hiblup_pc(out_prefix + ".pc", n_pcs);
        Eigen::VectorXd std_dev = parse_hiblup_pcp(out_prefix + ".pcp", n_pcs);

        if (P_unit.rows() != n_train) {
            throw std::runtime_error("hiblup returned " + std::to_string(P_unit.rows()) +
                                     " PCs but expected " + std::to_string(n_train));
        }

        // hiblup outputs unit-norm eigenvectors U (||P||^2=1).
        // std_dev = S / sqrt(n), so S = std_dev * sqrt(n)
        // Natural scale PCs: U * diag(S)
        double sqrt_n = std::sqrt(static_cast<double>(n_train));
        result.singular_values = std_dev * sqrt_n;
        result.train_pcs = P_unit * result.singular_values.asDiagonal();

    } else if (tool == "plink") {
        Eigen::MatrixXd P_unit = parse_plink_eigenvec(out_prefix + ".eigenvec", n_pcs);
        Eigen::VectorXd eigenval = parse_plink_eigenval(out_prefix + ".eigenval", n_pcs);

        if (P_unit.rows() != n_train) {
            throw std::runtime_error("plink returned " + std::to_string(P_unit.rows()) +
                                     " PCs but expected " + std::to_string(n_train));
        }

        // plink eigenval = eigenvalues of covariance matrix Z*Z'/(n-1)
        // S^2 = eigenval * (n-1)
        int n_minus_1 = n_train - 1;
        result.singular_values = (eigenval.array() * n_minus_1).sqrt();
        result.train_pcs = P_unit * result.singular_values.asDiagonal();
    }

    return result;
}

// ============================================================
// project_test_external()
// test_PCs = m * G_cross * U * diag(1/S)
//          = m * G_cross * train_pcs * diag(1/S^2)
// ============================================================
Eigen::MatrixXd project_test_external(
    const Eigen::MatrixXd& G_cross,
    const Eigen::MatrixXd& train_pcs,
    const Eigen::VectorXd& singular_values,
    int m_snps)
{
    // scale = m / S^2
    Eigen::VectorXd scale = (static_cast<double>(m_snps) /
                             singular_values.array().square()).matrix();
    // test_PCs = m * G_cross * U * diag(1/S)
    //          = G_cross * train_pcs * diag(m / S^2)
    return G_cross * train_pcs * scale.asDiagonal();
}

} // namespace fastcv
