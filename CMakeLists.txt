cmake_minimum_required(VERSION 3.15)
project(pygenogrove)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fetch pybind11
include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.13.6
)
FetchContent_MakeAvailable(pybind11)

# Add genogrove library
add_subdirectory(external/genogrove)

# Create Python module
pybind11_add_module(pygenogrove src/bindings.cpp)
target_link_libraries(pygenogrove PRIVATE genogrove)
target_include_directories(pygenogrove PRIVATE
    ${CMAKE_SOURCE_DIR}/external/genogrove/include
)

# Set output directory
set_target_properties(pygenogrove PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)