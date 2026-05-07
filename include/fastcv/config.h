#pragma once

#include "fastcv/types.h"
#include <string>

namespace fastcv {

/// Load CV configuration from JSON file
CvConfig load_config(const std::string& json_file);

/// Load CV configuration from JSON string
CvConfig load_config_from_string(const std::string& json_str);

/// Validate configuration and throw on error
void validate_config(const CvConfig& config);

/// Generate a template config JSON string
std::string generate_template_config();

} // namespace fastcv
