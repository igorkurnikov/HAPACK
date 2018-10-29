
bin_PROGRAMS = harlembin
bin_SCRIPTS  = harlem

#AM_LDFLAGS = -all-static -lz -lm -lpthread -lutil -ldl
AM_LDFLAGS = -lz -lutil 

AM_CPPFLAGS = \
  @HARLEM_MOD_DEFS@ \
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
 $(MPI_INCLUDE) \
 $(BOOST_INCLUDE) \
  $(WX_CPPFLAGS) $(OPEN_SCENE_GRAPH_VIEW_FLAG) \
  -Wno-non-template-friend -Wno-deprecated -fpermissive
harlembin_LDADD = \
 $(top_builddir)/HALIB/halib/halib.a \
 $(top_builddir)/DELPHI_LIB/libdelphi.a \
 $(top_builddir)/MY_VFLIB/libvf.a \
 $(top_builddir)/WOLFGANG/DO_LIB/libdo.a \
 $(top_builddir)/WOLFGANG/IPACK/libipack.a \
 $(top_builddir)/PNPS/libpnp/libpnp.a \
 $(top_builddir)/GEOBALL/libgeoball.a -lgmp \
 $(top_builddir)/JUMNA_LIB/libjumna.a  \
 $(top_srcdir)/../PROG_SRC/g98/util.a \
 -L$(top_srcdir)/../MYPROG/AMBER11_IGOR/AmberTools/src/gleap/mortsrc -lmort


#$(top_builddir)/SANDER_LIB/libsander.a

if WITH_LIB_PMEMD
harlembin_LDADD += $(top_builddir)/PMEMD_LIB/libpmemd.a
endif

if WITH_CUDA
harlembin_LDADD += $(NVCCLIBS)
endif

harlembin_LDADD += $(MPI_LDFLAGS) $(MPI_LIBS) $(BOOST_LIBS)\
 $(WX_LIBS) \
 $(FLIBS) $(BLAS_LIBS) $(LAPACK_LIBS)  \
 $(RESOURCES_BIN) \
 -L$(PYTHON_PREFIX)/lib $(PYTHON_LDFLAGS) \
 -lplplotwxwidgets -lplplot -lplplotcxx -lplplotwxwidgetsd -lplplotcxxd 
 

harlem_wx_files = \
 sources/dialogs_wx_1.cpp \
 sources/dialogs_wx_2.cpp \
 sources/mm_dialogs_wx.cpp \
 sources/qc_dialogs_wx.cpp \
 sources/wx_prot_redox_dlg.cpp \
 sources/hamainframe_wx.cpp \
 sources/harlemapp_wx.cpp \
 sources/ha_wx_aux_1.cpp \
 sources/ha_wx_ed_wdr.cpp \
 sources/ha_wx_res_wdr.cpp \
 sources/ha_wx_res_2_wdr.cpp \
 sources/ha_wx_res_mikola_wdr.cpp \
 sources/ha_wx_res_molflex_wdr.cpp \
 sources/wxMolED.cpp \
 sources/wxMolFlex.cpp \
 sources/hachart.cpp \
 sources/hachart_wx.cpp \
 sources/wxpymod.cpp 
 
harlem_main_files = \
 sources/apbsmod.cpp  \
 sources/abstree.cpp \
 sources/atom_mapping.cpp \
 sources/canvas3d.cpp \
 sources/command.cpp \
 sources/electrostmod.cpp \
 sources/elmod.cpp \
 sources/etcoupl.cpp \
 sources/fminbr.cpp \
 sources/gaufile.cpp \
 sources/infile.cpp \
 sources/haatbasdb.cpp \
 sources/haatgroup.cpp \
 sources/haatom.cpp \
 sources/haatombasis.cpp \
 sources/habond.cpp \
 sources/hacolors.cpp \
 sources/hacoord.cpp \
 sources/hacompmod.cpp \
 sources/hadalton.cpp \
 sources/hagaussian.cpp \
 sources/haempirical.cpp \
 sources/haflexmod.cpp  \
 sources/hahbhp.cpp \
 sources/haintcrd.cpp \
 sources/haintengine.cpp \
 sources/haintermol.cpp \
 sources/haio.cpp \
 sources/halinalg.cpp \
 sources/halocexcit.cpp \
 sources/halocorb.cpp \
 sources/hamatdb.cpp \
 sources/hamolecule.cpp \
 sources/hamolmech.cpp \
 sources/hamolset.cpp \
 sources/hamolview.cpp \
 sources/ha_mort_mm.cpp \
 sources/ha_mort_mol.cpp \
 sources/hamultipole.cpp \
 sources/harlemapp.cpp \
 sources/hapseudopot.cpp \
 sources/haresdb.cpp \
 sources/haproteined.cpp \
 sources/haqchem.cpp \
 sources/hampi.cpp \
 sources/harpaham.cpp \
 sources/harpavec.cpp \
 sources/hascattermod.cpp \
 sources/hasimulator.cpp \
 sources/hastl.cpp \
 sources/hastring.cpp \
 sources/hasurface.cpp \
 sources/hatests.cpp \
 sources/hatypes.cpp \
 sources/hazindo.cpp \
 sources/march_cube.cpp \
 sources/mm_driver_amber.cpp \
 sources/mm_driver_gromacs.cpp \
 sources/mm_elements.cpp \
 sources/mm_force_field.cpp \
 sources/mm_model.cpp \
 sources/mm_params.cpp \
 sources/mm_traj_anal.cpp \
 sources/moleditor.cpp \
 sources/molmech_evt_handler.cpp \
 sources/molset_evt_handler.cpp \
 sources/ndo.cpp \
 sources/nuclacidmod.cpp \
 sources/object3d.cpp \
 sources/outfile.cpp \
 sources/protonredox.cpp \
 sources/repres.cpp \
 sources/qc_params.cpp \
 sources/rigidbodycoord.cpp \
 sources/script.cpp \
 sources/stmmod.cpp \
 sources/trajanal.cpp \
 sources/transfor.cpp \
 sources/vec3d.cpp \
 sources/axx_fit.f \
 sources/delphutil.f \
 sources/gextr.f \
 sources/gmath.f \
 sources/gsurf.f \
 sources/math/fd05.f \
 sources/math/id05.f \
 sources/math/mc11.f \
 sources/math/va13.f \
 sources/math/za02.f \
 $(ADD_SRC)

 harlembin_SOURCES = $(harlem_main_files) $(harlem_wx_files)\
 hamolset_wrap.cpp \
 hasvnrev.cpp

#if WITH_LIB_PMEMD
#harlembin_SOURCES += sources/mm_driver_amber.cpp
#endif


hamolset_wrap.cpp:$(top_srcdir)/HARLEM/sources/hamolset.i
	$(SWIG) -python -shadow -c++ $(INCLUDES) $(harlem_DEFVAR) -o hamolset_wrap.cpp $(top_srcdir)/HARLEM/sources/hamolset.i


molset.py:hamolset_wrap.cpp

hasvnrev.cpp: $(harlem_main_files) $(harlem_wx_files) sources/hamolset.i
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

CLEANFILES = hamolset_wrap.cpp hasvnrev.cpp molset.py harlemwin.res harlembin.exe.manifest \
	harlembin.lib harlembin.exp harlem.exe 

######################################
wx_libs_dirs := $(filter -L%,${WX_LIBS})

harlem:
	@echo "#!/bin/bash">harlem
	@echo "#!/bin/bash">harlem
	@echo "#">>harlem
	@echo "# Script to run HARLEM">>harlem
	@echo "#">>harlem
	@echo "export HARLEM_ROOTDIR=${prefix}">>harlem
	@echo "export HARLEM_HOME=\$${HARLEM_ROOTDIR}/">>harlem
	@echo "export PATH=\"\$${HARLEM_ROOTDIR}/bin:\$$PATH\"">>harlem
#paths to python interfaces
	@echo "export PYTHONPATH=\"\$${PYTHONPATH}:${pkgscriptdir}\"">>harlem
#Paths to shared libraries
	@echo "export LD_LIBRARY_PATH=\"\$${LD_LIBRARY_PATH}\\">>harlem
#mkl
	@echo ":${MKL_LIBS}\\"|sed 's,-L,,g'|sed 's,  *,:,g'|sed 's,\t\t*,:,g'>>harlem
#wxWidgets
	@echo ":${wx_libs_dirs}\\"|sed 's,-L,,g'|sed 's,  *,:,g'|sed 's,\t\t*,:,g'>>harlem       
#harlem's shared library, not yet done
	@echo ":\$${HARLEM_ROOTDIR}/lib\\">>harlem
	@echo ":\$${HARLEM_ROOTDIR}/lib/HARLEM\\">>harlem
#LIBSSH
	@echo ":${LIBSSH_DIR}/libssh\"">>harlem
	@echo "echo \$$LD_LIBRARY_PATH">>harlem
	p="harlembin";f=`echo "$$p" | sed 's,^.*/,,;s/$(EXEEXT)$$//;$(transform);s/$$/$(EXEEXT)/'`;\
	echo "nice -n 0 $(bindir)/$$f \$$1 \$$2 \$$3 \$$4 \$$5 \$$6 \$$7 \$$8 \$$9">>harlem
	@chmod 755 harlem
	@cat harlem

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



