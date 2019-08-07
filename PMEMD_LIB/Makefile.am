noinst_LIBRARIES = libpmemd.a
if WITH_GFORTRAN
AM_FCFLAGS = -ffree-form -DDIRFRC_COMTRANS -DFFTLOADBAL_2PROC -DINCLUDE_MPIF $(MPI_INCLUDE) -I/opt/mpich/include/f90base  -fPIC
else
AM_FCFLAGS = -fpe0  -DDIRFRC_COMTRANS -DFFTLOADBAL_2PROC -DINCLUDE_MPIF $(MPI_INCLUDE)   -fPIC
endif

AM_CPPFLAGS = -fPIC \
              -I$(top_srcdir)/PMEMD_LIB/src -I/opt/intel/mkl/10.0.1.014/include 

VPATH = $(top_srcdir)/PMEMD_LIB/src

CLEANFILES = *.i *.mod


libpmemd_a_SOURCES =   alltasks_setup.F90 amoeba_adjust.F90 amoeba_bonds.F90 \
        amoeba_direct.F90 amoeba_flags.F90 amoeba_force.F90 amoeba_induced.F90 \
	amoeba_interface.F90 amoeba_multipoles.F90 amoeba_opbend_angles.F90 \
	amoeba_pitorsions.F90 amoeba_recip.F90 amoeba_reg_angles.F90 amoeba_self.F90 \
	amoeba_stretch_bend.F90 amoeba_stretch_torsions.F90 amoeba_torsion_torsion.F90 \
	amoeba_torsions.F90 amoeba_trig_angles.F90 amoeba_ureyb.F90  amoeba_valence.F90  \
	amoeba_vdw.F90 angles.F90  bonds.F90  bspline.F90  cit.F90 constraints.F90  \
	dihedrals.F90 dynamics.F90 dynamics_dat.F90 ene_frc_splines.F90 erfcfun.F90  \
	file_io.F90 file_io_add.F90 file_io_dat.F90  fft1d.F90 force_proxy.F90 \
        gb_ene.F90 gb_force.F90 gbl_constants.F90 gbl_datatypes.F90 img.F90  loadbal.F90 \
	mdin_amoeba_dat.F90 mdin_ctrl_dat.F90 mdin_ewald_dat.F90  \
	nb_pairlist.F90 nonbond_constr.F90 \
	parallel.F90 parallel_dat.F90 pbc.F90 pme_direct.F90 pme_fft.F90 pme_force.F90  \
	pme_recip.F90 pme_setup.F90  pmemd_clib.c prmtop_dat.F90 \
	random.F90 runfiles.F90 shake.F90  state_info.F90 timers.F90 veclib.F90

include $(top_srcdir)/PMEMD_LIB/src/f90depends
