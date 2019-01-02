noinst_LIBRARIES = libjumna.a

if WITH_GFORTRAN
AM_FFLAGS = -ffixed-line-length-132 -fPIC
else
AM_FFLAGS = -132 -fPIC
endif

INCLUDES = -I$(top_srcdir)/JUMNA_LIB/sources

libjumna_a_SOURCES=\
sources/jumcall.f sources/ampar91.f sources/ampar94.f \
sources/ang.f sources/assemb.f sources/axeout.f   \
sources/backbo.f sources/build.f sources/cataxe.f \
sources/close.f sources/combi.f sources/cosdir.f sources/deltor.f   \
sources/delval.f sources/dis.f sources/disth.f \
sources/ecomp.f sources/ecomp91.f sources/ecomp94.f   \
sources/enecyl.f sources/energ91.f sources/energ94.f \
sources/energy.f sources/equim.f sources/equiv.f sources/fitsug.f   \
sources/flex.f sources/grads.f sources/gradt.f \
sources/helix.f sources/helloc.f sources/joie.f   \
sources/kapa.f sources/kline.f sources/ligaxe.f \
sources/ligput.f sources/liner.f sources/lmoout.f   \
sources/loops.f sources/mc11a.f sources/microb.f \
sources/minfor.f sources/minim.f   \
sources/molin.f sources/move.f sources/openb.f \
sources/pairc.f sources/pairs.f sources/penalty.f   \
sources/present.f sources/putbac.f sources/renoe.f \
sources/reset.f sources/setbac.f sources/setd.f   \
sources/setgeo.f sources/setvar.f sources/tbasl.f \
sources/tjunl.f sources/torp.f $(ADD_SRC)


