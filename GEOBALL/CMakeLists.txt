cmake_minimum_required(VERSION 3.18)
project(geoball_lib LANGUAGES C Fortran)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Legacy Fortran code - requires fixed-form with extended line length
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-form -ffixed-line-length-none -std=legacy")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fixed -extend-source")
endif()

set(GEOBALL_C_SOURCES
    GeoBall.c
    alf_gmp.c
    sos_minor_gmp.c
)

set(GEOBALL_FORTRAN_SOURCES
    adjust.f
    alfcx.f
    compute_vol.f
    delcx.f
    prepare_deriv.f
    spacefill_vol.f
    truncate_real.f
    vector.f
)

add_library(geoball STATIC ${GEOBALL_C_SOURCES} ${GEOBALL_FORTRAN_SOURCES})

target_include_directories(geoball PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link GMP library if available
find_library(GMP_LIBRARY gmp)
if(GMP_LIBRARY)
    target_link_libraries(geoball PUBLIC ${GMP_LIBRARY})
endif()
