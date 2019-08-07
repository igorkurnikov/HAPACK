noinst_LIBRARIES = libpnp.a
pkgscript_DATA = pnpsll.py
AM_CPPFLAGS = \
 -I$(top_srcdir)/HARLEM/sources \
 -I$(top_srcdir)/HARLEM/sources/math \
 -I$(top_srcdir)/HARLEM/sources/gl \
 -I$(top_srcdir)/HARLEM/sources/mVue \
 -I$(top_srcdir)/HALIB/include \
 -I$(top_srcdir)/MY_VFLIB/src \
 -I$(top_srcdir)/WOLFGANG/DO_LIB/include \
 -I$(top_srcdir)/WOLFGANG/DO_LIB/template \
 -I$(top_srcdir)/WOLFGANG/IPACK  \
 -I$(top_srcdir)/WOLFGANG/generic \
 -I$(top_srcdir)/PNPS/src \
 $(OPEN_SCENE_GRAPH_VIEW_INCLUDE) \
 $(BOOST_INCLUDE) \
 $(LIBSSH_INCLUDE)
# -I$(top_srcdir)/PNPS/src 

if WITH_CUDA
AM_CPPFLAGS += -DWITH_CUDA $(NVCCHOSTFLAGS) -D DBG0 -D DBG1 -D DBG2 @HARLEM_MOD_DEFS@  $(PYTHON_CPPFLAGS) $(WX_CPPFLAGS) $(OPEN_SCENE_GRAPH_VIEW_FLAG) $(MPI_INCLUDE)  -fPIC
else
AM_CPPFLAGS += -D DBG0 -D DBG1 -D DBG2 @HARLEM_MOD_DEFS@ $(PYTHON_CPPFLAGS) $(WX_CPPFLAGS) $(OPEN_SCENE_GRAPH_VIEW_FLAG) $(MPI_INCLUDE) -fPIC
endif

#  -Wno-non-template-friend -Wno-deprecated
if WITH_CUDA
libpnp_a_SOURCES = \
  ./../src/pnps.cpp ./../src/pnpinterfaces.cpp ./../src/pnpsapp.cpp ./../src/contworld.cpp ./../src/buildworldni.cpp ./../src/mapio.cpp ./../src/poissonsolver.cpp ./../src/poissonboltzmannsolver.cpp ./../src/nernstplanksolver.cpp ./../src/poissonnernstplancksolver.cpp ./../src/pnputil.cpp ./../src/pbwithljsolver.cpp pnpsll_wrap.cpp ./../src/sissgausps.cpp ./../src/pnpstructs.cpp ./../src/cudasolver.cu  ./../src/cudasolverdouble.cu ./../src/cudabldworld.cu
else
libpnp_a_SOURCES = \
  ./../src/pnps.cpp ./../src/pnpinterfaces.cpp ./../src/pnpsapp.cpp ./../src/contworld.cpp ./../src/buildworldni.cpp ./../src/mapio.cpp ./../src/poissonsolver.cpp ./../src/poissonboltzmannsolver.cpp ./../src/nernstplanksolver.cpp ./../src/poissonnernstplancksolver.cpp  ./../src/pnputil.cpp ./../src/pbwithljsolver.cpp pnpsll_wrap.cpp ./../src/sissgausps.cpp ./../src/pnpstructs.cpp ./../src/pnpdebug.cpp
endif


cudasolver_dep = ./../src/pnpstructs.cpp

.cu.o:
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<


pnpsll_wrap_dep = \
 ./../src/pnpmod.h \
 ./../src/mapio.h \
 ./../src/contworld.h \
 ./../src/calcctrl.h \
 ./../src/pmfcalculation.h \
 ./../src/pnpsapp.h \
 ./../src/calcctrl.h \
 ./../src/hasurface.h

pnpsll_wrap.cpp:$(top_srcdir)/PNPS/src/pnpsll.i
	$(SWIG) -python -shadow -c++ -py3 $(INCLUDES) @HARLEM_MOD_DEFS@ -o pnpsll_wrap.cpp $(top_srcdir)/PNPS/src/pnpsll.i

pnpsll.py:pnpsll_wrap.cpp


