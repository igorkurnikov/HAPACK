SET(LIB_MORTGL_SOURCE 
glutils.cpp
glstyle.cpp
graphic.cpp
drawing.cpp
ribbon.cpp
mouse.cpp
ortho.cpp
)

ADD_LIBRARY(mortgl SHARED ${LIB_MORTGL_SOURCE})
INSTALL_TARGETS(/lib mortgl)
TARGET_LINK_LIBRARIES(mortgl /usr/lib/libGL.so /usr/lib/libGLU.so mort)
