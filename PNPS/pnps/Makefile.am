bin_PROGRAMS = pnps

if WITH_CUDA
AM_CPPFLAGS = -DWITH_CUDA $(NVCCHOSTFLAGS) -I$(top_srcdir)/HALIB/include  -I$(top_srcdir)/PNPS/src  @ZLIB_CPPFLAGS@ \
			-D HALIB_BASE -D PNPONLY  -D DBG0 -D DBG1 -D DBG2 
else
AM_CPPFLAGS = -I$(top_srcdir)/HALIB/include  -I$(top_srcdir)/PNPS/src  @ZLIB_CPPFLAGS@ \
			-D HALIB_BASE -D PNPONLY  -D DBG0 -D DBG1 -D DBG2
endif
#CXX=nvcc
#NVCCFLAGS = -g --ptxas-options=-v --device-emulation

AM_LDFLAGS = 
LDADD = ../../HALIB/halib_base/libhalib_base.a @ZLIB_LIBS@ -lcudart -L/usr/local/cuda/lib
if WITH_CUDA
pnps_SOURCES = ./../src/pnps.cpp ./../src/pnpinterfaces.cpp ./../src/pnpsapp.cpp ./../src/contworld.cpp ./../src/buildworldni.cpp ./../src/mapio.cpp ./../src/poissonsolver.cpp ./../src/poissonboltzmannsolver.cpp ./../src/nernstplanksolver.cpp ./../src/poissonnernstplancksolver.cpp ./../src/calcctrl.cpp ./../src/pmfcalculation.cpp ./../src/pnputil.cpp ./../src/pbwithljsolver.cpp ./../src/sissgausps.cpp ./../src/pnpstructs.cpp ./../src/cudasolver.cu ./../src/cudasolverdouble.cu ./../src/cudabldworld.cu
else
pnps_SOURCES = ./../src/pnps.cpp ./../src/pnpinterfaces.cpp ./../src/pnpsapp.cpp ./../src/contworld.cpp ./../src/buildworldni.cpp ./../src/mapio.cpp ./../src/poissonsolver.cpp ./../src/poissonboltzmannsolver.cpp ./../src/nernstplanksolver.cpp ./../src/poissonnernstplancksolver.cpp ./../src/calcctrl.cpp ./../src/pmfcalculation.cpp ./../src/pnputil.cpp ./../src/pbwithljsolver.cpp ./../src/sissgausps.cpp ./../src/pnpstructs.cpp
endif


cudasolver_dep = ./../src/pnpstructs.cpp

.cu.o:
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<

