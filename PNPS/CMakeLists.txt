cmake_minimum_required(VERSION 3.18)
project(pnps_lib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find Python for headers
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

set(PNPS_SOURCES
    src/buildworldni.cpp
    src/contworld.cpp
    src/mapio.cpp
    src/nernstplanksolver.cpp
    src/pbwithljsolver.cpp
    # src/pmfcalculation.cpp  # Missing header pmfcalculation.h
    src/pnpdebug.cpp
    src/pnpinterfaces.cpp
    src/pnpsapp.cpp
    src/pnps.cpp
    src/pnpstructs.cpp
    src/pnputil.cpp
    src/poissonboltzmannsolver.cpp
    src/poissonnernstplancksolver.cpp
    src/poissonsolver.cpp
    src/sissgausps.cpp
)

add_library(pnps STATIC ${PNPS_SOURCES})

target_include_directories(pnps PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/tinyxml
    ${CMAKE_SOURCE_DIR}/MOLSET/sources
)

# Link Python for header access
target_link_libraries(pnps PUBLIC Python::Module)

# Optional CUDA support
option(PNPS_ENABLE_CUDA "Enable CUDA support for PNPS" OFF)
if(PNPS_ENABLE_CUDA)
    enable_language(CUDA)
    target_sources(pnps PRIVATE
        src/cudabldworld.cu
        src/cudasolver.cu
        src/cudasolverdouble.cu
    )
endif()
