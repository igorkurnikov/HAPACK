noinst_LIBRARIES = libdo.a
AM_CPPFLAGS = -DGCC_NEW -DIPACK -fPIC \
 -I$(top_srcdir)/WOLFGANG/DO_LIB \
 -I$(top_srcdir)/WOLFGANG/DO_LIB/include \
 -I$(top_srcdir)/WOLFGANG/DO_LIB/template \
 $(MPI_INCLUDE)
libdo_a_SOURCES=stringc.cpp inverse.cpp determinant.cpp diagonalize.cpp xtime.cpp svd.cpp
