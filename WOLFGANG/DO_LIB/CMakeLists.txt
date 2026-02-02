cmake_minimum_required(VERSION 3.18)
project(do_lib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(DOLIB_SOURCES
    blas_implement.cpp
    broydn.cpp
    determinant.cpp
    diagonalize.cpp
    error.cpp
    inverse.cpp
    stringc.cpp
    svd.cpp
    xtime.cpp
)

add_library(dolib STATIC ${DOLIB_SOURCES})

target_include_directories(dolib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/template
)

# Try to find LAPACK/BLAS for optimized routines
find_package(LAPACK QUIET)
if(LAPACK_FOUND)
    target_link_libraries(dolib PUBLIC ${LAPACK_LIBRARIES})
    target_compile_definitions(dolib PRIVATE USE_LAPACK)
endif()
