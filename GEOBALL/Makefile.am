noinst_LIBRARIES = libgeoball.a

AM_FFLAGS = -ffixed-line-length-132 -fPIC 
AM_CPPFLAGS = -fPIC 
#AM_FFLAGS = -132
INCLUDES = -I$(top_srcdir)/GEOBALL

libgeoball_a_SOURCES=\
GeoBall.c \
sos_minor_gmp.c alf_gmp.c \
adjust.f delcx.f alfcx.f truncate_real.f \
compute_vol.f prepare_deriv.f spacefill_vol.f vector.f
