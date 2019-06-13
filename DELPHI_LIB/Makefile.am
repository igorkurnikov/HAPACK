noinst_LIBRARIES = libdelphi.a

if WITH_GFORTRAN
AM_FFLAGS = -ffixed-line-length-132 -fPIC
else
AM_FFLAGS = -132 -fPIC
endif


AM_CPPFLAGS = -I$(top_srcdir)/DELPHI_LIB/source
libdelphi_a_SOURCES=\
source/itit4j.f source/non.f source/relfac4b.f source/setcrg4.f \
source/mem4.f source/phintp4.f source/scaler4.f source/setin4.f \
source/mkeps4f.f source/qdiff4s.f source/setbc4.f source/setout4.f
