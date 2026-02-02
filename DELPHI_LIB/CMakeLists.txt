cmake_minimum_required(VERSION 3.18)
project(delphi_lib LANGUAGES Fortran)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Legacy Fortran 77 code - requires fixed-form with extended line length
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-form -ffixed-line-length-none -std=legacy")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fixed -extend-source")
endif()

set(DELPHI_SOURCES
    source/itit4j.f
    source/mem4.f
    source/mkeps4f.f
    source/non.f
    source/phintp4.f
    source/qdiff4s.f
    source/relfac4b.f
    source/scaler4.f
    source/setbc4.f
    source/setcrg4.f
    source/setin4.f
    source/setout4.f
)

add_library(delphi STATIC ${DELPHI_SOURCES})

target_include_directories(delphi PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/source
)
