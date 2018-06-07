bin_PROGRAMS = harlem_nogui_bin
bin_SCRIPTS  = harlem_nogui

#AM_LDFLAGS = -static-intel -lz -lm -lpthread -lutil -ldl
AM_LDFLAGS = -lz  -lutil 


AM_CPPFLAGS = \
  @HARLEM_MOD_DEFS@ -DHA_NOGUI \
 -I$(top_srcdir)/HARLEM/sources \
 -I$(top_srcdir)/HARLEM/sources/math \
 -I$(top_srcdir)/HARLEM/sources/gl \
 -I$(top_srcdir)/HALIB/include \
 -I$(top_srcdir)/MY_VFLIB/src \
 -I$(top_srcdir)/WOLFGANG/DO_LIB/include \
 -I$(top_srcdir)/WOLFGANG/DO_LIB/template \
 -I$(top_srcdir)/WOLFGANG/IPACK  \
 -I$(top_srcdir)/WOLFGANG/generic \
 -I$(top_srcdir)/PNPS/src \
 -I$(top_srcdir)/PNPS/pnpmod \
 -I$(top_srcdir)/../MYPROG/AMBER11_IGOR/AmberTools/src/gleap/mortsrc \
 $(PYTHON_CPPFLAGS) \
 $(BOOST_INCLUDE) \
 $(MPI_INCLUDE) \
  $(WX_CPPFLAGS) \
  -Wno-non-template-friend -Wno-deprecated -fpermissive 

harlem_nogui_bin_LDADD = \
 $(top_builddir)/HALIB/halib/halib.a \
 $(top_builddir)/DELPHI_LIB/libdelphi.a \
 $(top_builddir)/MY_VFLIB/libvf.a \
 $(top_builddir)/WOLFGANG/DO_LIB/libdo.a \
 $(top_builddir)/WOLFGANG/IPACK/libipack.a \
 $(top_builddir)/PNPS/libpnp/libpnp.a \
 $(top_builddir)/GEOBALL/libgeoball.a  -lgmp \
 $(top_builddir)/JUMNA_LIB/libjumna.a \
 $(top_srcdir)/../PROG_SRC/g98/util.a \
 -L$(top_srcdir)/../MYPROG/AMBER11_IGOR/AmberTools/src/gleap/mortsrc -lmort
 
if WITH_LIB_PMEMD
harlem_nogui_bin_LDADD += $(top_builddir)/PMEMD_LIB/libpmemd.a
endif
 
harlem_nogui_bin_LDADD +=  $(MPI_LDFLAGS) $(MPI_LIBS) $(BOOST_LIBS)\
 $(WX_LIBS_BASE) \
 $(FLIBS) $(BLAS_LIBS) $(LAPACK_LIBS)  \
 $(RESOURCES_BIN) \
 -L$(PYTHON_PREFIX)/lib $(PYTHON_LDFLAGS)
     
VPATH = $(top_srcdir)/HARLEM/sources $(top_srcdir)/HARLEM/sources/math
     
harlem_main_files = \
 abstree.cpp \
 apbsmod.cpp  \
 atom_mapping.cpp \
 canvas3d.cpp \
 command.cpp \
 electrostmod.cpp \
 elmod.cpp \
 etcoupl.cpp \
 fminbr.cpp \
 gaufile.cpp \
 infile.cpp \
 haatbasdb.cpp \
 haatgroup.cpp \
 haatom.cpp \
 haatombasis.cpp \
 habond.cpp \
 hacolors.cpp \
 hacompmod.cpp \
 hacoord.cpp \
 hadalton.cpp \
 hagaussian.cpp \
 haempirical.cpp \
 haflexmod.cpp  \
 hahbhp.cpp \
 haintcrd.cpp \
 haintengine.cpp \
 haintermol.cpp \
 haio.cpp \
 halinalg.cpp \
 halocexcit.cpp \
 halocorb.cpp \
 hamatdb.cpp \
 hamolecule.cpp \
 hamolmech.cpp \
 hamolset.cpp \
 hamolview.cpp \
 ha_mort_mm.cpp \
 ha_mort_mol.cpp \
 hamultipole.cpp \
 harlemapp.cpp \
 harlemapp_wx.cpp \
 hapseudopot.cpp \
 haresdb.cpp \
 haproteined.cpp \
 haqchem.cpp \
 hampi.cpp \
 harpaham.cpp \
 harpavec.cpp \
 hascattermod.cpp \
 hasimulator.cpp \
 hastl.cpp \
 hastring.cpp \
 hasurface.cpp \
 hatests.cpp \
 hatypes.cpp \
 hazindo.cpp \
 march_cube.cpp \
 mm_driver_amber.cpp \
 mm_elements.cpp \
 mm_force_field.cpp \
 mm_model.cpp \
 mm_params.cpp \
 mm_traj_anal.cpp \
 moleditor.cpp \
 molmech_evt_handler.cpp \
 molset_evt_handler.cpp \
 ndo.cpp \
 nuclacidmod.cpp \
 object3d.cpp \
 outfile.cpp \
 protonredox.cpp \
 qc_params.cpp \
 repres.cpp \
 rigidbodycoord.cpp \
 script.cpp \
 stmmod.cpp \
 trajanal.cpp \
 transfor.cpp \
 vec3d.cpp \
 axx_fit.f \
 delphutil.f \
 gextr.f \
 gmath.f \
 gsurf.f \
 fd05.f \
 id05.f \
 mc11.f \
 va13.f \
 za02.f \
 hawx_add.cpp  \
 $(ADD_SRC)

harlem_nogui_bin_SOURCES = $(harlem_main_files) \
 hamolset_wrap.cpp \
 hasvnrev.cpp

hamolset_wrap.cpp:$(top_srcdir)/HARLEM/sources/hamolset.i
	$(SWIG) -python -shadow -c++ $(INCLUDES) $(harlem_DEFVAR) -DHA_NOGUI -o hamolset_wrap.cpp $(top_srcdir)/HARLEM/sources/hamolset.i

molset.py:hamolset_wrap.cpp

hasvnrev.cpp: $(harlem_main_files) $(harlem_wx_files) $(top_srcdir)/HARLEM/sources/hamolset.i
	echo '#include "hasvnrev.h"'> hasvnrev.cpp
	echo -n 'const char* HaSVNRevision(){const char* SVN_Version = "' >> hasvnrev.cpp
	echo -n `svn info $(top_srcdir)/|grep "Revision"` >> hasvnrev.cpp
	echo '"; return SVN_Version; }'   >> hasvnrev.cpp
	echo -n 'const char* HaSVNDate(){const char* SVN_Date = "' >> hasvnrev.cpp
	echo -n `svn info $(top_srcdir)/|grep "Last Changed Date"` >> hasvnrev.cpp
	echo '"; return SVN_Date; }'   >> hasvnrev.cpp
#svnversion -n $(top_srcdir)/ >> hasvnrev.cpp
FORCE:

pkgscript_DATA = molset.py 

CLEANFILES = hamolset_wrap.cpp hasvnrev.cpp molset.py harlemwin.res \
	harlembin.lib harlembin.exp 

######################################
wx_libs_dirs := $(filter -L%,${WX_LIBS})

harlem_nogui:
	@echo "#!/bin/bash">harlem_nogui
	@echo "#!/bin/bash">harlem_nogui
	@echo "#">>harlem_nogui
	@echo "# Script to run HARLEM">>harlem_nogui
	@echo "#">>harlem_nogui
	@echo "export HARLEM_ROOTDIR=${prefix}">>harlem_nogui
	@echo "export HARLEM_HOME=\$${HARLEM_ROOTDIR}/">>harlem_nogui
	@echo "export PATH=\"\$${HARLEM_ROOTDIR}/bin:\$$PATH\"">>harlem_nogui
#paths to python interfaces
	@echo "export PYTHONPATH=\"\$${PYTHONPATH}:${pkgscriptdir}\"">>harlem_nogui
#Paths to shared libraries
	@echo "export LD_LIBRARY_PATH=\"\$${LD_LIBRARY_PATH}\\">>harlem_nogui
#mkl
	@echo ":${MKL_LIBS}\\"|sed 's,-L,,g'|sed 's,  *,:,g'|sed 's,\t\t*,:,g'>>harlem_nogui
#wxWidgets
	@echo ":${wx_libs_dirs}\\"|sed 's,-L,,g'|sed 's,  *,:,g'|sed 's,\t\t*,:,g'>>harlem_nogui       
#harlem's shared library, not yet done
	@echo ":\$${HARLEM_ROOTDIR}/lib\\">>harlem_nogui
	@echo ":\$${HARLEM_ROOTDIR}/lib/HARLEM\\">>harlem_nogui
#LIBSSH
	@echo ":${LIBSSH_DIR}/libssh\"">>harlem_nogui
	@echo "echo \$$LD_LIBRARY_PATH">>harlem_nogui
	p="harlem_nogui_bin";f=`echo "$$p" | sed 's,^.*/,,;s/$(EXEEXT)$$//;$(transform);s/$$/$(EXEEXT)/'`;\
	echo "nice -n 0 $(bindir)/$$f \$$1 \$$2 \$$3 \$$4 \$$5 \$$6 \$$7 \$$8 \$$9">>harlem_nogui
	@chmod 755 harlem_nogui
	@cat harlem_nogui

#harlemtest now is shell which allowed to run Harlem from build directory
#, without need to install
noinst_SCRIPTS = harlemtestnow
harlemtestnow:
	@echo "#!/bin/bash">harlemtestnow
	@echo "#!/bin/bash">harlemtestnow
	@echo "#">>harlemtestnow
	@echo "# Script to run HARLEM">>harlemtestnow
	@echo "#">>harlemtestnow
	@echo "export HARLEM_ROOTDIR=${abs_top_builddir}/HARLEM">>harlemtestnow
	@echo "export HARLEM_HOME=${abs_top_srcdir}/HARLEM/">>harlemtestnow
	@echo "export PATH=\"\$$PATH\\">>harlemtestnow
	@echo ":${abs_top_builddir}/HARLEM\\">>harlemtestnow
	@echo "\"">>harlemtestnow
#paths to python interfaces
	@echo "export PYTHONPATH=\"\$${PYTHONPATH}\\">>harlemtestnow
	@echo ":${abs_top_builddir}/HARLEM\\">>harlemtestnow
	@echo ":${abs_top_builddir}/HALIB/halib\\">>harlemtestnow
	@echo ":${abs_top_builddir}/PNPS/libpnp\\">>harlemtestnow
	@echo "\"">>harlemtestnow
#Paths to shared libraries
	@echo "export LD_LIBRARY_PATH=\"\$${LD_LIBRARY_PATH}\\">>harlemtestnow
#mkl
	@echo ":${MKL_LIBS}\\"|sed 's,-L,,g'|sed 's,  *,:,g'|sed 's,\t\t*,:,g'>>harlemtestnow
#harlem's shared library, not yet done
#	@echo ":\$${HARLEM_ROOTDIR}/lib\\">>harlemtestnow
#	@echo ":\$${HARLEM_ROOTDIR}/lib/HARLEM\\">>harlemtestnow
#LIBSSH
	@echo ":${LIBSSH_DIR}/libssh\"">>harlemtestnow
	@echo "echo \$$LD_LIBRARY_PATH">>harlemtestnow
	@echo "nice -n 0 ${abs_top_builddir}/HARLEM/harlembin \$$1 \$$2 \$$3 \$$4 \$$5 \$$6 \$$7 \$$8 \$$9">>harlemtestnow
	@chmod 755 harlemtestnow
	@cat harlemtestnow
# 	chmod 755 harlem



