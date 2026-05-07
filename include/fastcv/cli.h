#pragma once

#include "fastcv/types.h"
#include <string>
#include <vector>

namespace fastcv {

/// Run the CLI with given arguments
/// @param args  Command-line arguments (excluding program name)
/// @return Exit code (0 = success)
int run_cli(const std::vector<std::string>& args);

} // namespace fastcv
