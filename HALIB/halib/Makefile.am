noinst_LIBRARIES = libhalib.a
pkgscript_DATA = halib.py
AM_CPPFLAGS = \
 -I$(top_srcdir)/HARLEM/sources \
 -I$(top_srcdir)/HARLEM/sources/math \
 -I$(top_srcdir)/HARLEM/sources/gl \
 -I$(top_srcdir)/HARLEM/sources/mVue \
 -I$(top_srcdir)/HALIB/include \
  $(PYTHON_CPPFLAGS) \
  @HARLEM_MOD_DEFS@ -fPIC

libhalib_a_SOURCES=./../src/haobject.cpp   ./../src/haxml.cpp   ./../src/haconsts.cpp ./../src/tinyxml.cpp ./../src/tinyxmlparser.cpp ./../src/tinyxmlerror.cpp halib_wrap.cpp

# halib_wrap_dep = \
#  ./../src/haobject.h \
#  ./../src/haconsts.h \
#  ./../src/tinyxml.h

halib_wrap.cpp:$(srcdir)/../include/halib.i
	$(SWIG) -python -shadow -c++ @HARLEM_MOD_DEFS@ -o halib_wrap.cpp $(srcdir)/../include/halib.i

halib.py:halib_wrap.cpp

CLEANFILES = halib_wrap.cpp halib.py
