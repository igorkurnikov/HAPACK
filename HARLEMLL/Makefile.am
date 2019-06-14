lib_LTLIBRARIES = lib_harlem.la

AM_LDFLAGS = -lz -lutil -shared 

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
 -I$(top_srcdir)/../AMBER11_IGOR/AmberTools/src/gleap/mortsrc \
 $(PYTHON_CPPFLAGS) \
 $(MPI_INCLUDE) \
 $(BOOST_INCLUDE) \
  $(WX_CPPFLAGS) $(OPEN_SCENE_GRAPH_VIEW_FLAG) \
  -Wno-non-template-friend -Wno-deprecated -fpermissive 

lib_harlem_la_LIBADD = \
 $(top_builddir)/HALIB/halib/libhalib.a \
 $(top_builddir)/DELPHI_LIB/libdelphi.a \
 $(top_builddir)/MY_VFLIB/libvf.a \
 $(top_builddir)/WOLFGANG/DO_LIB/libdo.a \
 $(top_builddir)/WOLFGANG/IPACK/libipack.a \
 $(top_builddir)/PNPS/libpnp/libpnp.a \
 $(top_builddir)/GEOBALL/libgeoball.a -lgmp \
 $(top_builddir)/JUMNA_LIB/libjumna.a  \
 $(top_srcdir)/../../PROG_SRC/g98/util.a \
 -L$(top_srcdir)/../AMBER11_IGOR/AmberTools/src/gleap/mortsrc -lmort

if WITH_LIB_PMEMD
lib_harlem_la_LIBADD += $(top_builddir)/PMEMD_LIB/libpmemd.a
endif

if WITH_CUDA
lib_harlem_la_LIBADD += $(NVCCLIBS)
endif

lib_harlem_la_LIBADD += $(MPI_LDFLAGS) $(MPI_LIBS) $(BOOST_LIBS)\
 $(WX_LIBS) \
 $(FLIBS) $(BLAS_LIBS) $(LAPACK_LIBS)  \
 $(RESOURCES_BIN) \
 -L$(PYTHON_PREFIX)/lib $(PYTHON_LDFLAGS) \
 -lplplotwxwidgets -lplplot -lplplotcxx -lplplotwxwidgets -lplplotcxx 
 

harlem_wx_files = \
 $(top_srcdir)/HARLEM/sources/dialogs_wx_1.cpp \
 $(top_srcdir)/HARLEM/sources/dialogs_wx_2.cpp \
 $(top_srcdir)/HARLEM/sources/mm_dialogs_wx.cpp \
 $(top_srcdir)/HARLEM/sources/qc_dialogs_wx.cpp \
 $(top_srcdir)/HARLEM/sources/wx_prot_redox_dlg.cpp \
 $(top_srcdir)/HARLEM/sources/hamainframe_wx.cpp \
 $(top_srcdir)/HARLEM/sources/harlemapp_wx.cpp \
 $(top_srcdir)/HARLEM/sources/ha_wx_aux_1.cpp \
 $(top_srcdir)/HARLEM/sources/ha_wx_ed_wdr.cpp \
 $(top_srcdir)/HARLEM/sources/ha_wx_res_wdr.cpp \
 $(top_srcdir)/HARLEM/sources/ha_wx_res_2_wdr.cpp \
 $(top_srcdir)/HARLEM/sources/ha_wx_res_mikola_wdr.cpp \
 $(top_srcdir)/HARLEM/sources/ha_wx_res_molflex_wdr.cpp \
 $(top_srcdir)/HARLEM/sources/wxMolED.cpp \
 $(top_srcdir)/HARLEM/sources/wxMolFlex.cpp \
 $(top_srcdir)/HARLEM/sources/hachart.cpp \
 $(top_srcdir)/HARLEM/sources/hachart_wx.cpp \
 $(top_srcdir)/HARLEM/sources/wxpymod.cpp 
 
harlem_main_files = \
 $(top_srcdir)/HARLEM/sources/apbsmod.cpp  \
 $(top_srcdir)/HARLEM/sources/abstree.cpp \
 $(top_srcdir)/HARLEM/sources/atom_mapping.cpp \
 $(top_srcdir)/HARLEM/sources/canvas3d.cpp \
 $(top_srcdir)/HARLEM/sources/command.cpp \
 $(top_srcdir)/HARLEM/sources/electrostmod.cpp \
 $(top_srcdir)/HARLEM/sources/elmod.cpp \
 $(top_srcdir)/HARLEM/sources/etcoupl.cpp \
 $(top_srcdir)/HARLEM/sources/fminbr.cpp \
 $(top_srcdir)/HARLEM/sources/gaufile.cpp \
 $(top_srcdir)/HARLEM/sources/infile.cpp \
 $(top_srcdir)/HARLEM/sources/haatbasdb.cpp \
 $(top_srcdir)/HARLEM/sources/haatgroup.cpp \
 $(top_srcdir)/HARLEM/sources/haatom.cpp \
 $(top_srcdir)/HARLEM/sources/haatombasis.cpp \
 $(top_srcdir)/HARLEM/sources/habond.cpp \
 $(top_srcdir)/HARLEM/sources/hacolors.cpp \
 $(top_srcdir)/HARLEM/sources/hacoord.cpp \
 $(top_srcdir)/HARLEM/sources/hacompmod.cpp \
 $(top_srcdir)/HARLEM/sources/hadalton.cpp \
 $(top_srcdir)/HARLEM/sources/hagaussian.cpp \
 $(top_srcdir)/HARLEM/sources/haempirical.cpp \
 $(top_srcdir)/HARLEM/sources/haflexmod.cpp  \
 $(top_srcdir)/HARLEM/sources/hahbhp.cpp \
 $(top_srcdir)/HARLEM/sources/haintcrd.cpp \
 $(top_srcdir)/HARLEM/sources/haintengine.cpp \
 $(top_srcdir)/HARLEM/sources/haintermol.cpp \
 $(top_srcdir)/HARLEM/sources/haio.cpp \
 $(top_srcdir)/HARLEM/sources/halinalg.cpp \
 $(top_srcdir)/HARLEM/sources/halocexcit.cpp \
 $(top_srcdir)/HARLEM/sources/halocorb.cpp \
 $(top_srcdir)/HARLEM/sources/hamatdb.cpp \
 $(top_srcdir)/HARLEM/sources/hamolecule.cpp \
 $(top_srcdir)/HARLEM/sources/hamolmech.cpp \
 $(top_srcdir)/HARLEM/sources/hamolset.cpp \
 $(top_srcdir)/HARLEM/sources/hamolview.cpp \
 $(top_srcdir)/HARLEM/sources/ha_mort_mm.cpp \
 $(top_srcdir)/HARLEM/sources/ha_mort_mol.cpp \
 $(top_srcdir)/HARLEM/sources/hamultipole.cpp \
 $(top_srcdir)/HARLEM/sources/harlemapp.cpp \
 $(top_srcdir)/HARLEM/sources/hapseudopot.cpp \
 $(top_srcdir)/HARLEM/sources/haresdb.cpp \
 $(top_srcdir)/HARLEM/sources/haproteined.cpp \
 $(top_srcdir)/HARLEM/sources/haqchem.cpp \
 $(top_srcdir)/HARLEM/sources/hampi.cpp \
 $(top_srcdir)/HARLEM/sources/harpaham.cpp \
 $(top_srcdir)/HARLEM/sources/harpavec.cpp \
 $(top_srcdir)/HARLEM/sources/hascattermod.cpp \
 $(top_srcdir)/HARLEM/sources/hasimulator.cpp \
 $(top_srcdir)/HARLEM/sources/hastl.cpp \
 $(top_srcdir)/HARLEM/sources/hastring.cpp \
 $(top_srcdir)/HARLEM/sources/hasurface.cpp \
 $(top_srcdir)/HARLEM/sources/hatests.cpp \
 $(top_srcdir)/HARLEM/sources/hatypes.cpp \
 $(top_srcdir)/HARLEM/sources/hazindo.cpp \
 $(top_srcdir)/HARLEM/sources/march_cube.cpp \
 $(top_srcdir)/HARLEM/sources/mm_driver_amber.cpp \
 $(top_srcdir)/HARLEM/sources/mm_driver_gromacs.cpp \
 $(top_srcdir)/HARLEM/sources/mm_elements.cpp \
 $(top_srcdir)/HARLEM/sources/mm_force_field.cpp \
 $(top_srcdir)/HARLEM/sources/mm_model.cpp \
 $(top_srcdir)/HARLEM/sources/mm_params.cpp \
 $(top_srcdir)/HARLEM/sources/mm_traj_anal.cpp \
 $(top_srcdir)/HARLEM/sources/moleditor.cpp \
 $(top_srcdir)/HARLEM/sources/molmech_evt_handler.cpp \
 $(top_srcdir)/HARLEM/sources/molset_evt_handler.cpp \
 $(top_srcdir)/HARLEM/sources/ndo.cpp \
 $(top_srcdir)/HARLEM/sources/nuclacidmod.cpp \
 $(top_srcdir)/HARLEM/sources/object3d.cpp \
 $(top_srcdir)/HARLEM/sources/outfile.cpp \
 $(top_srcdir)/HARLEM/sources/protonredox.cpp \
 $(top_srcdir)/HARLEM/sources/repres.cpp \
 $(top_srcdir)/HARLEM/sources/qc_params.cpp \
 $(top_srcdir)/HARLEM/sources/rigidbodycoord.cpp \
 $(top_srcdir)/HARLEM/sources/script.cpp \
 $(top_srcdir)/HARLEM/sources/stmmod.cpp \
 $(top_srcdir)/HARLEM/sources/trajanal.cpp \
 $(top_srcdir)/HARLEM/sources/transfor.cpp \
 $(top_srcdir)/HARLEM/sources/vec3d.cpp \
 $(top_srcdir)/HARLEM/sources/axx_fit.f \
 $(top_srcdir)/HARLEM/sources/delphutil.f \
 $(top_srcdir)/HARLEM/sources/gextr.f \
 $(top_srcdir)/HARLEM/sources/gmath.f \
 $(top_srcdir)/HARLEM/sources/gsurf.f \
 $(top_srcdir)/HARLEM/sources/math/fd05.f \
 $(top_srcdir)/HARLEM/sources/math/id05.f \
 $(top_srcdir)/HARLEM/sources/math/mc11.f \
 $(top_srcdir)/HARLEM/sources/math/va13.f \
 $(top_srcdir)/HARLEM/sources/math/za02.f \
 $(ADD_SRC)

lib_harlem_la_SOURCES = $(harlem_main_files) $(harlem_wx_files)\
 hamolset_wrap.cpp 

#if WITH_LIB_PMEMD
#lib_harlem_la_SOURCES += $(top_srcdir)/HARLEM/sources/mm_driver_amber.cpp
#endif


hamolset_wrap.cpp:$(top_srcdir)/HARLEM/sources/hamolset.i
	$(SWIG) -python -shadow -py3 -c++ $(INCLUDES) $(harlem_DEFVAR) -o hamolset_wrap.cpp $(top_srcdir)/HARLEM/sources/hamolset.i


molset.py:hamolset_wrap.cpp


pkgscript_DATA = molset.py 

CLEANFILES = hamolset_wrap.cpp molset.py lib_harlem.so 








