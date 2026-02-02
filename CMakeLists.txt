cmake_minimum_required(VERSION 3.18)
project(molset_bindings LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find required packages
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(SWIG REQUIRED)
find_package(ZLIB REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem system)
include(UseSWIG)

# --- Build dependency libraries via subdirectories ---
add_subdirectory(MORT_LIB)
add_subdirectory(MY_VFLIB)
add_subdirectory(JUMNA_LIB)
add_subdirectory(GEOBALL)
add_subdirectory(DELPHI_LIB)
add_subdirectory(PNPS)
# WOLFGANG/DO_LIB and WOLFGANG/IPACK have legacy C++ template code
# that requires source changes for modern compilers. Disabled for now.
# add_subdirectory(WOLFGANG)

# --- Directory definitions for MOLSET ---
set(MOLSET_DIR "${CMAKE_SOURCE_DIR}/MOLSET")

# --- MOLSET sources (C/C++ and Fortran) ---
file(GLOB MOLSET_CXX_SOURCES
    "${MOLSET_DIR}/sources/*.c"
    "${MOLSET_DIR}/sources/*.cc"
    "${MOLSET_DIR}/sources/*.cxx"
    "${MOLSET_DIR}/sources/*.cpp"
)

file(GLOB MOLSET_FORTRAN_SOURCES
    "${MOLSET_DIR}/sources/*.f"
    "${MOLSET_DIR}/sources/*.F"
    "${MOLSET_DIR}/sources/*.f90"
    "${MOLSET_DIR}/sources/*.F90"
    "${MOLSET_DIR}/sources/math/*.f"
    "${MOLSET_DIR}/sources/math/*.F"
)

set(MOLSET_SOURCES ${MOLSET_CXX_SOURCES} ${MOLSET_FORTRAN_SOURCES})

# --- SWIG interface ---
# The SWIG interface file defines %module molsetc
set(SWIG_I_FILE "${MOLSET_DIR}/sources/hamolset.i")

# Where SWIG writes the generated Python proxy (molsetc.py)
set(CMAKE_SWIG_OUTDIR "${CMAKE_CURRENT_BINARY_DIR}/swig_py")

# SWIG include paths for %include directives
set(CMAKE_SWIG_FLAGS
    "-I${MOLSET_DIR}/sources"
    "-I${CMAKE_SOURCE_DIR}/MORT_LIB"
    "-I${CMAKE_SOURCE_DIR}/MY_VFLIB/src"
    "-I${CMAKE_SOURCE_DIR}/JUMNA_LIB/sources"
    "-I${CMAKE_SOURCE_DIR}/GEOBALL"
    "-I${CMAKE_SOURCE_DIR}/DELPHI_LIB/source"
    "-I${CMAKE_SOURCE_DIR}/PNPS/src"
)

set_source_files_properties(${SWIG_I_FILE} PROPERTIES
    CPLUSPLUS ON
    USE_TARGET_INCLUDE_DIRECTORIES ON
)

# Create the SWIG Python module (_molsetc)
# Module name matches %module molsetc in the .i file
swig_add_library(molsetc
    TYPE MODULE
    LANGUAGE python
    SOURCES ${SWIG_I_FILE} ${MOLSET_SOURCES}
)

target_include_directories(molsetc PRIVATE
    "${MOLSET_DIR}/sources"
    "${MOLSET_DIR}/sources/math"
    "${CMAKE_SOURCE_DIR}/WOLFGANG/DO_LIB/include"
    "${CMAKE_SOURCE_DIR}/WOLFGANG/DO_LIB/template"
    "${CMAKE_SOURCE_DIR}/WOLFGANG/IPACK"
    "${CMAKE_SOURCE_DIR}/WOLFGANG/generic"
)

# Define HARLEM_MOD to avoid constant redefinitions with PNPS
# Define CLK_TCK for legacy code (deprecated in modern systems)
target_compile_definitions(molsetc PRIVATE
    HARLEM_MOD
    CLK_TCK=CLOCKS_PER_SEC
)

# Link all dependency libraries
set(MOLSET_LINK_LIBS
    mort
    vflib
    geoball
    delphi
    pnps
    Python::Module
    ZLIB::ZLIB
    ${LAPACK_LIBRARIES}
    Boost::filesystem
    Boost::system
)

# Add jumna if it was built (conditional in JUMNA_LIB/CMakeLists.txt)
if(TARGET jumna)
    list(APPEND MOLSET_LINK_LIBS jumna)
endif()

# Add Fortran runtime if we have Fortran sources
if(MOLSET_FORTRAN_SOURCES OR TARGET jumna OR TARGET geoball OR TARGET delphi)
    # Link gfortran runtime on GCC
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        list(APPEND MOLSET_LINK_LIBS gfortran)
    endif()
endif()

swig_link_libraries(molsetc ${MOLSET_LINK_LIBS})

# --- Install into the wheel ---
# scikit-build-core installs into the wheel
# Package directory matches wheel.packages in pyproject.toml
install(TARGETS molsetc
    LIBRARY DESTINATION molset
    RUNTIME DESTINATION molset
    ARCHIVE DESTINATION molset
)

# Install the SWIG-generated molsetc.py proxy
install(FILES "${CMAKE_SWIG_OUTDIR}/molsetc.py"
    DESTINATION molset
)
