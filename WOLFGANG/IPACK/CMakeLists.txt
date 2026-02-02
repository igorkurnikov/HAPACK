cmake_minimum_required(VERSION 3.18)
project(ipack_lib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(IPACK_SOURCES
    alloc.cpp
    angular.cpp
    basis.cpp
    coef_array.cpp
    coef_set.cpp
    const.cpp
    contract.cpp
    e_field.cpp
    element.cpp
    f.cpp
    fix.cpp
    flexfile.cpp
    four_stor.cpp
    function.cpp
    hrr1.cpp
    hrr2.cpp
    hrr_set.cpp
    input.cpp
    integ_array.cpp
    integ.cpp
    integ_file.cpp
    ipack.cpp
    kinetic.cpp
    memmgr.cpp
    moment.cpp
    nuc_e_field.cpp
    nuclear.cpp
    nuc_mom.cpp
    nuc_rep.cpp
    one_el.cpp
    operators.cpp
    orbinfo.cpp
    output1.cpp
    output2.cpp
    output3.cpp
    output.cpp
    overlap.cpp
    print.cpp
    qc_utilities.cpp
    quad.cpp
    rec.cpp
    rec_stor.cpp
    sphere.cpp
    spin_orb.cpp
    storage.cpp
    symdesig.cpp
    symlogic.cpp
    symop.cpp
    tpair.cpp
    tquadruple.cpp
    transpose.cpp
    twoel0.cpp
    twoel.cpp
    types.cpp
    vrr.cpp
)

add_library(ipack STATIC ${IPACK_SOURCES})

target_include_directories(ipack PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# IPACK depends on DO_LIB
if(TARGET dolib)
    target_link_libraries(ipack PUBLIC dolib)
endif()
