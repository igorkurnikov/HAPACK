SET(LIB_PYCOMMON_SOURCE pycommon.cpp)
ADD_LIBRARY(pycommon MODULE ${LIB_PYCOMMON_SOURCE})
INSTALL_TARGETS(/lib pycommon)
TARGET_LINK_LIBRARIES(pycommon /usr/lib/libboost_python.so)
TARGET_LINK_LIBRARIES(pycommon mort)

SET(LIB_PYOBJECT_SOURCE pyobject.cpp)
ADD_LIBRARY(pyobject MODULE ${LIB_PYOBJECT_SOURCE})
INSTALL_TARGETS(/lib pyobject)
TARGET_LINK_LIBRARIES(pyobject /usr/lib/libboost_python.so)
TARGET_LINK_LIBRARIES(pyobject mort)

SET(LIB_PYMORTGL_SOURCE pymortgl.cpp)
ADD_LIBRARY(pymortgl MODULE ${LIB_PYMORTGL_SOURCE})
INSTALL_TARGETS(/lib pymortgl)
TARGET_LINK_LIBRARIES(pymortgl /usr/lib/libboost_python.so)
TARGET_LINK_LIBRARIES(pymortgl mort mortgl)

INSTALL_FILES(/lib .py hashcode pymort)
