cmake_minimum_required(VERSION 3.18)
project(mort_lib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find Boost (required for algorithm/string and numeric/ublas)
# Check system location first
set(BOOST_ROOT "/usr")
find_package(Boost REQUIRED)

set(MORT_SOURCES
    ambfmt/ambfrc-amoeba.cpp
    ambfmt/amboff.cpp
    ambfmt/chain.cpp
    ambfmt/energee.cpp
    ambfmt/exclude.cpp
    ambfmt/extent.cpp
    ambfmt/les.cpp
    ambfmt/parm.cpp
    ambfmt/parmset.cpp
    ambfmt/prep.cpp
    ambfmt/prmtop-amoeba.cpp
    ambfmt/prmtop.cpp
    atmask/atom.cpp
    atmask/condition.cpp
    atmask/dist.cpp
    atmask/logic.cpp
    atmask/mask.cpp
    atmask/paren.cpp
    atmask/resd.cpp
    capbox/addions.cpp
    capbox/octree.cpp
    capbox/shaper.cpp
    capbox/solute.cpp
    capbox/solvate.cpp
    capbox/solvent.cpp
    capbox/region.cpp
    common/numvec.cpp
    common/matrix.cpp
    common/pertab.cpp
    common/fortran.cpp
    common/stralgo.cpp
    common/crdalgo.cpp
    common/funstack.cpp
    common/geometry.cpp
    common/hashcode.cpp
    enefrc/bonding.cpp
    enefrc/ctrlparm.cpp
    enefrc/nonbond-lst.cpp
    enefrc/nonbond-egb.cpp
    enefrc/nonbond-pbc.cpp
    format/mol.cpp
    format/mdb.cpp
    format/ext.cpp
    format/frc.cpp
    format/top.cpp
    format/sdf.cpp
    format/xyz.cpp
    guilib/command.cpp
    guilib/console.cpp
    guilib/control.cpp
    guilib/grammar.cpp
    guilib/mainwin.cpp
    guilib/leaplog.cpp
    object/atom.cpp
    object/angl.cpp
    object/bond.cpp
    object/dihe.cpp
    object/impr.cpp
    object/resd.cpp
    object/atmv.cpp
    object/rela.cpp
    object/data.cpp
    object/comp.cpp
    object/root.cpp
    object/mole.cpp
    object/morf.cpp
    objfun/wrapper.cpp
    objfun/atomfun.cpp
    objfun/bondfun.cpp
    objfun/fixbond.cpp
    objfun/geomfun.cpp
    objfun/hydrfun.cpp
    objfun/merge.cpp
    objfun/parmfun.cpp
    pdbent/atom.cpp
    pdbent/bond.cpp
    pdbent/resd.cpp
    pdbent/mdlize.cpp
    pdbent/second.cpp
    pdbent/pdbent.cpp
    pdbent/namemap.cpp
    smarts/parse.cpp
    smarts/reader.cpp
    smarts/status.cpp
    smarts/finder.cpp
    smarts/atomexpr.cpp
    smarts/tyrule.cpp
    smarts/ring.cpp
    smarts/subst.cpp
    smarts/smarts.cpp
    tripos/atom.cpp
    tripos/bond.cpp
    tripos/head.cpp
    tripos/resd.cpp
    tripos/sect.cpp
    tripos/mol2.cpp
)

add_library(mort STATIC ${MORT_SOURCES})

target_link_libraries(mort PUBLIC Boost::boost)

target_include_directories(mort PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/ambfmt
    ${CMAKE_CURRENT_SOURCE_DIR}/atmask
    ${CMAKE_CURRENT_SOURCE_DIR}/capbox
    ${CMAKE_CURRENT_SOURCE_DIR}/common
    ${CMAKE_CURRENT_SOURCE_DIR}/enefrc
    ${CMAKE_CURRENT_SOURCE_DIR}/format
    ${CMAKE_CURRENT_SOURCE_DIR}/guilib
    ${CMAKE_CURRENT_SOURCE_DIR}/object
    ${CMAKE_CURRENT_SOURCE_DIR}/objfun
    ${CMAKE_CURRENT_SOURCE_DIR}/pdbent
    ${CMAKE_CURRENT_SOURCE_DIR}/smarts
    ${CMAKE_CURRENT_SOURCE_DIR}/tripos
)
