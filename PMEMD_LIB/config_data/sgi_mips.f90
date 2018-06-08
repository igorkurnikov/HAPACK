DIRFRC_DEFINES = -DDIRFRC_EFS -DDIRFRC_NOVEC
CPP = /lib/cpp
CPPFLAGS = -P
F90_DEFINES =

F90 = f90
MODULE_SUFFIX = mod
CPUFLAGS = $(MIPS_PROCESSOR_TYPE) -mips4 -n32
F90FLAGS = $(CPUFLAGS) -nocpp -old_rl -c

F90_OPT_DBG = -g
F90_OPT_LO =  -O1 -OPT:roundoff=3:IEEE_arithmetic=3:fast_bit_intrinsics=ON
F90_OPT_MED = -O2 -OPT:roundoff=3:IEEE_arithmetic=3:fast_bit_intrinsics=ON
F90_OPT_HI =  -O3 -OPT:roundoff=3:IEEE_arithmetic=3:fast_bit_intrinsics=ON
F90_OPT_DFLT =  $(F90_OPT_HI)

CC = cc
CFLAGS = $(CPUFLAGS)

LOAD = f90
LOADFLAGS = $(CPUFLAGS)
LOADLIBS = -lfastm
