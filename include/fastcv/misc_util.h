#pragma once

#include <string>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <sys/stat.h>

namespace fastcv {
namespace detail {

inline void print_progress_bar(int pct, const char* label = nullptr) {
    int bar_width = 30;
    int filled = pct * bar_width / 100;
    std::string bar(bar_width, ' ');
    std::fill(bar.begin(), bar.begin() + filled, '=');
    std::cerr << "\r        ";
    if (label) std::cerr << "[" << label << "] ";
    std::cerr << "[" << bar << "] " << std::setw(3) << pct << "%" << std::flush;
}

} // namespace detail

inline void mkdir_p(const std::string& path) {
    for (size_t i = 1; i < path.size(); ++i) {
        if (path[i] == '/') {
            std::string sub = path.substr(0, i);
            ::mkdir(sub.c_str(), 0755);
        }
    }
    ::mkdir(path.c_str(), 0755);
}

inline std::string make_fold_dir(const std::string& base, int fold_idx) {
    return base + "/fold_" + (fold_idx < 10 ? "0" : "") + std::to_string(fold_idx);
}

} // namespace fastcv
