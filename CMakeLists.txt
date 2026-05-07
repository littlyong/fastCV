cmake_minimum_required(VERSION 3.16)
project(fastCV VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ============================================================
# Options
# ============================================================
option(FASTCV_BUILD_CLI "Build the fastcv CLI executable" ON)
option(FASTCV_BUILD_SHARED "Build shared library (OFF = static)" OFF)

# ============================================================
# All dependencies are in third_party/ — no system libs needed
# ============================================================
set(THIRD_PARTY "${CMAKE_CURRENT_SOURCE_DIR}/third_party")

# --- Eigen (header-only, vendored) ---
set(EIGEN_INCLUDE_DIR "${THIRD_PARTY}/eigen")
if(NOT EXISTS "${EIGEN_INCLUDE_DIR}/Eigen/Dense")
    message(FATAL_ERROR "Eigen not found at ${EIGEN_INCLUDE_DIR}")
endif()
message(STATUS "Eigen:       ${EIGEN_INCLUDE_DIR}")

# --- Spectra (header-only, vendored) ---
set(SPECTRA_INCLUDE_DIR "${THIRD_PARTY}/Spectra")
if(NOT EXISTS "${SPECTRA_INCLUDE_DIR}/SymEigsSolver.h")
    message(FATAL_ERROR "Spectra not found at ${SPECTRA_INCLUDE_DIR}")
endif()
message(STATUS "Spectra:     ${SPECTRA_INCLUDE_DIR}")

# --- nlohmann/json (header-only, vendored) ---
set(JSON_INCLUDE_DIR "${THIRD_PARTY}")
if(NOT EXISTS "${JSON_INCLUDE_DIR}/json.hpp")
    message(FATAL_ERROR "nlohmann/json not found at ${JSON_INCLUDE_DIR}")
endif()
message(STATUS "JSON:        ${JSON_INCLUDE_DIR}/json.hpp")

# --- OpenBLAS (compiled from source) ---
set(OPENBLAS_DIR "${THIRD_PARTY}/OpenBLAS")
if(EXISTS "${OPENBLAS_DIR}/libopenblas.a")
    message(STATUS "OpenBLAS:    ${OPENBLAS_DIR}/libopenblas.a")
else()
    message(WARNING "OpenBLAS not found — run: cd third_party/OpenBLAS && make NO_FORTRAN=1 USE_OPENMP=1 -j16")
endif()

# --- Intel MKL (via conda env) ---
set(MKL_DIR "${THIRD_PARTY}/mkl_env")
if(EXISTS "${MKL_DIR}/lib/libmkl_gf_lp64.a")
    message(STATUS "MKL:         ${MKL_DIR}/lib/")
    set(FASTCV_HAS_MKL 1)
else()
    message(STATUS "MKL:         not found")
    set(FASTCV_HAS_MKL 0)
endif()

# --- OpenMP: always use find_package to get -fopenmp compiler flag ---
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    set(GOMP_LIB OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP:      ${OpenMP_CXX_LIBRARIES} (compiler flags: ${OpenMP_CXX_FLAGS})")
else()
    message(WARNING "OpenMP not found — will build single-threaded")
endif()

# ============================================================
# BLAS selection: MKL > OpenBLAS > none
# ============================================================
if(FASTCV_HAS_MKL)
    # MKL GNU thread + gomp for internal OpenMP parallelism
    # Use --start-group to resolve circular dependencies between MKL static libs
    set(BLAS_LIBRARIES
        -Wl,--start-group
        ${MKL_DIR}/lib/libmkl_gf_lp64.a
        ${MKL_DIR}/lib/libmkl_gnu_thread.a
        ${MKL_DIR}/lib/libmkl_core.a
        -Wl,--end-group
        ${GOMP_LIB}
        -lpthread -lm -ldl)
    set(BLAS_INCLUDE_DIR "${MKL_DIR}/include")
    set(BLAS_TYPE "Intel MKL")
elseif(EXISTS "${OPENBLAS_DIR}/libopenblas.a")
    set(BLAS_LIBRARIES
        ${OPENBLAS_DIR}/libopenblas.a
        ${GOMP_LIB}
        -lpthread -lm)
    set(BLAS_INCLUDE_DIR "${OPENBLAS_DIR}")
    set(BLAS_TYPE "OpenBLAS (local)")
else()
    set(BLAS_TYPE "none (Eigen built-in)")
endif()
message(STATUS "BLAS:        ${BLAS_TYPE}")

# ============================================================
# fastcv library
# ============================================================
set(FASTCV_SOURCES
    src/plink_reader.cpp
    src/grm.cpp
    src/dim_reduction.cpp
    src/cross_validation.cpp
    src/pc_corrector.cpp
    src/data_export.cpp
    src/config.cpp
    src/external_pca.cpp
    src/ld_prune.cpp
    src/memory_estimate.cpp
)

add_library(fastcv ${FASTCV_SOURCES})
target_include_directories(fastcv PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${EIGEN_INCLUDE_DIR}>
)
target_include_directories(fastcv PRIVATE
    ${SPECTRA_INCLUDE_DIR}
    ${JSON_INCLUDE_DIR}
)
target_compile_features(fastcv PUBLIC cxx_std_17)
target_compile_options(fastcv PUBLIC -march=native)

# OpenMP
if(GOMP_LIB)
    target_link_libraries(fastcv PUBLIC ${GOMP_LIB})
    target_compile_definitions(fastcv PUBLIC FASTCV_HAS_OPENMP)
    # Ensure MKL conda env lib dir is in rpath for runtime
    if(FASTCV_HAS_MKL)
        set_target_properties(fastcv PROPERTIES
            BUILD_RPATH "${MKL_DIR}/lib"
            INSTALL_RPATH "${MKL_DIR}/lib"
        )
    endif()
endif()

# BLAS
if(BLAS_LIBRARIES)
    target_link_libraries(fastcv PUBLIC ${BLAS_LIBRARIES})
    target_compile_definitions(fastcv PUBLIC FASTCV_HAS_BLAS)
    if(BLAS_INCLUDE_DIR)
        target_include_directories(fastcv PRIVATE ${BLAS_INCLUDE_DIR})
    endif()
    if(FASTCV_HAS_MKL)
        target_compile_definitions(fastcv PUBLIC FASTCV_USE_MKL)
    endif()
endif()

# ============================================================
# fastcv CLI executable
# ============================================================
if(FASTCV_BUILD_CLI)
    add_executable(fastcv_cli
        src/cli.cpp
        src/main.cpp
    )
    target_link_libraries(fastcv_cli PRIVATE fastcv)
    if(FASTCV_HAS_MKL)
        set_target_properties(fastcv_cli PROPERTIES
            BUILD_RPATH "${MKL_DIR}/lib"
            INSTALL_RPATH "${MKL_DIR}/lib"
        )
    endif()
    set_target_properties(fastcv_cli PROPERTIES OUTPUT_NAME fastcv)
endif()

# ============================================================
# Install rules
# ============================================================
include(GNUInstallDirs)
install(TARGETS fastcv
    EXPORT fastcv-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/fastcv DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT fastcv-targets
    FILE fastcv-targets.cmake
    NAMESPACE fastcv::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fastcv
)
