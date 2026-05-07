#pragma once

#include "fastcv/types.h"
#include "fastcv/dim_reduction.h"

namespace fastcv {

/// Main pipeline: prepare CV data with fold-wise PC correction.
/// Orchestrates fold creation, correction, and returns all results.
///
/// @param config  CV configuration (paths, parameters, etc.)
/// @return Complete CV result with all fold data
CvResult prepare_cv_data(const CvConfig& config);

} // namespace fastcv
