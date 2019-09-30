PROJECT(GLEAP)

cmake_minimum_required(VERSION 2.0)

SET(CMAKE_INSTALL_PREFIX $ENV{AMBERHOME} CACHE INTERNAL "prefix")

FIND_PATH(PYTHONINC Python.h $ENV{AMBERHOME}/include /usr/include/ /usr/local/include/ )
FIND_PATH(BOOSTPINC python.h $ENV{AMBERHOME}/include/boost-1_34_1 /usr/include/ /usr/local/include/ )
FIND_LIBRARY(PYTHONLIB NAMES python $ENV{AMBERHOME}/lib/ python /usr/lib/ /usr/local/lib/ )
FIND_LIBRARY(BOOSTPLIB NAMES boost_python-1_34_1 $ENV{AMBERHOME}/lib/ /usr/lib/ /usr/local/lib )

INCLUDE_DIRECTORIES(
    ${GLEAP_SOURCE_DIR}/mortsrc
    ${GLEAP_SOURCE_DIR}/freelib
    ${PYTHONINC}
    ${BOOSTPINC}
)

OPTION(CMAKE_VERBOSE_MAKEFILE "if all commands will be echoed to the console during the make" ON)

SUBDIRS(mortsrc)

