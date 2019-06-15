lib_LTLIBRARIES = libgeoball.la

AM_FFLAGS = -ffixed-line-length-132 -fPIC 
AM_CPPFLAGS = -fPIC -I$(top_srcdir)/GEOBALL
#AM_FFLAGS = -132

libgeoball_la_SOURCES=\
GeoBall.c \
alf_gmp.c sos_minor_gmp.c \
adjust.f delcx.f alfcx.f truncate_real.f \
compute_vol.f prepare_deriv.f spacefill_vol.f vector.f

