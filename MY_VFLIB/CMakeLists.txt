cmake_minimum_required(VERSION 3.18)
project(vflib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(VFLIB_SOURCES
    src/argedit.cpp
    src/argloader.cpp
    src/argraph.cpp
    src/error.cpp
    src/gene.cpp
    src/gene_mesh.cpp
    src/match.cpp
    src/sd_state.cpp
    src/sortnodes.cpp
    src/ull_state.cpp
    src/ull_sub_state.cpp
    src/vf2_mono_state.cpp
    src/vf2_state.cpp
    src/vf2_sub_state.cpp
    src/vf_mono_state.cpp
    src/vf_state.cpp
    src/vf_sub_state.cpp
    src/xsubgraph.cpp
)

add_library(vflib STATIC ${VFLIB_SOURCES})

target_include_directories(vflib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
