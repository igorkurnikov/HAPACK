cmake_minimum_required(VERSION 3.18)
project(jumna_lib LANGUAGES Fortran)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Legacy Fortran code - requires fixed-form with extended line length
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-form -ffixed-line-length-none -std=legacy")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fixed -extend-source")
endif()

file(GLOB JUMNA_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/*.f
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/*.F
)

if(JUMNA_SOURCES)
    add_library(jumna STATIC ${JUMNA_SOURCES})
    target_include_directories(jumna PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/sources
    )
else()
    message(STATUS "No JUMNA Fortran sources found")
endif()
