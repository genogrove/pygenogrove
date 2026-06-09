cmake_minimum_required(VERSION 3.15)
project(pygenogrove)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Single source of truth for the package version: pyproject.toml. scikit-build-core
# passes it in as SKBUILD_PROJECT_VERSION, which we forward to the module as a
# compile definition (see target_compile_definitions below) so the version lives
# in exactly one place. A plain `cmake` dev build (no scikit-build-core) gets a
# clearly-non-release marker instead of a stale hardcoded number.
if(DEFINED SKBUILD_PROJECT_VERSION)
    set(PYGENOGROVE_VERSION "${SKBUILD_PROJECT_VERSION}")
else()
    set(PYGENOGROVE_VERSION "0.0.0+local")
endif()

# Fetch pybind11
include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.13.6
)
FetchContent_MakeAvailable(pybind11)

# Add genogrove library — build it as a library only. Its test/CLI/benchmark
# targets need generated headers (e.g. config/version.hpp) and extra deps we
# don't use here, and would otherwise break the wheel build.
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_CLI OFF CACHE BOOL "" FORCE)
set(BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
# genogrove's static lib must be position-independent so it can be linked
# into the pybind11 shared module (.so); otherwise ld rejects the relocation.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_subdirectory(external/genogrove)

# Create Python module
pybind11_add_module(pygenogrove src/bindings.cpp)
target_link_libraries(pygenogrove PRIVATE genogrove)
target_include_directories(pygenogrove PRIVATE
    ${CMAKE_SOURCE_DIR}/external/genogrove/include
    # genogrove's version macros live in a header it generates at configure time
    # into its own binary dir; expose it so bindings.cpp can report which
    # genogrove it was built against (pygenogrove.__genogrove_version__).
    ${genogrove_BINARY_DIR}/include
)

# PYGENOGROVE_VERSION -> pygenogrove.__version__ (single-sourced from pyproject.toml).
target_compile_definitions(pygenogrove PRIVATE
    PYGENOGROVE_VERSION="${PYGENOGROVE_VERSION}"
)

# Set output directory
set_target_properties(pygenogrove PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Install the extension into the wheel root so `import pygenogrove` resolves.
# scikit-build-core populates the wheel from the CMake install tree; without
# this rule the module is built but never packaged, leaving it unimportable.
install(TARGETS pygenogrove LIBRARY DESTINATION .)