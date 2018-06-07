!*******************************************************************************
!
! Module: mdin_ctrl_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module mdin_ctrl_dat_mod

use file_io_dat_mod
use mdin_amoeba_dat_mod

  implicit none

! The mdin ctrl input options are grouped in the order of input, with the
! common blocks containing valid options and the private data containing
! either deprecated options or sander 9 options not supported by pmemd.
! Some variables here actually are not &ctrl namelist parameters, but values
! derived from them that are needed during mdin &ctrl printout, etc.

  integer,save          ::      ntf, ibelly, &
                                nrespa, ntp, &
                                ntc, jfastw, igb, alpb, rbornstat, &
                                nrespai, iamoeba 

  character(4),save     ::      ihwtnm1,ihwtnm2,iowtnm, iwtnm


  ! New variables for pmemd.  These are set based on actual input of the values,
  ! or alternatively on the value of the 'cut' input variable, or based on
  ! default.  They are what actually gets used to determine cutoffs for pme:
  ! es_cutoff - electrostatic cutoff.
  ! vdw_cutoff - vdw cutoff.

  double precision, save     :: dielc, es_cutoff, vdw_cutoff, &                    
                                scnb, scee, intdiel, extdiel, saltcon, &          
                                rgbmax, offset, surften, cut_inner, gb_cutoff, &      
                                gb_alpha, gb_beta, gb_gamma, gb_fs_max, &        
                                gb_kappa, gb_neckscale, arad, &                         
                                bbox_xmin, bbox_ymin, bbox_zmin, &                     
                                bbox_xmax, bbox_ymax, bbox_zmax                       
                              

  ! Note - gb_alpha, gb_beta, gb_gamma, gb_fs_max, gb_kappa and gb_neckscale
  !        would be better factored elsewhere, but as usual the sander order
  !        of operations makes this difficult...

  ! Undocumented control variables, for dev use:

  integer, save :: dbg_atom_redistribution      ! force more frequent atom
                                                ! redistribution; actually get
                                                ! better testing now in assoc
                                                ! with fft slab redistribution.
  integer, save :: loadbal_verbose              ! loadbalancing verbosity.

  ! Logical variables used to indicate the type of potential to use.
  ! Currently the options include PME and Generalized Born.  These are
  ! set in the broadcast module for mpi, and are based on the igb values.
  ! Direct use of igb is undesirable because it designates multiple gb
  ! implementations.

  logical, save :: using_gb_potential
  logical, save :: using_pme_potential
  
  logical, save :: dirfrc_efs = .true.    ! flag to compute direct part of pme sum using splines
 
end module mdin_ctrl_dat_mod
 
subroutine set_mdin_ctrl_int(vals)
! Set values of mdin_ctrl_int  structure
	use mdin_ctrl_dat_mod
	use mdin_amoeba_dat_mod
	
	implicit none
	integer, intent(in)    :: vals(*)

    integer ixx
	character(4)          watdef(4)
	! Define default water residue name and the names of water oxygen & hydrogens
    data watdef/'WAT ', 'O   ', 'H1  ', 'H2  '/

	ixx    = vals(1)
	ixx    = vals(2) 
	ixx    = vals(3)
	ixx    = vals(4) 
	ixx    = vals(5)
	ixx    = vals(6)  
    ixx    = vals(7) 
    ixx    = vals(8) 
    ixx    = vals(9) 
    ixx    = vals(10) 
    ixx    = vals(11)  
    ixx    = vals(12) 
    ixx    = vals(13) 
    ixx    = vals(14) 
    ixx    = vals(15) 
    ntf    = vals(16) 
    ixx    = vals(17) 
    ixx    = vals(18) 
    ixx    = vals(19) 
    ibelly = vals(20) 
    ixx    = vals(21)  
    ixx    = vals(22)  
    ixx    = vals(23) 
    ixx    = vals(24) 
    ixx    = vals(25) 
    ixx    = vals(26) 
    nrespa = vals(27)  
    ixx    = vals(28) 
    ixx    = vals(29) 
    ixx    = vals(30) 
    ntp    = vals(31) 
    ntc    = vals(32) 
    jfastw = vals(33) 
    ixx    = vals(34)  
    igb    = vals(35) 
    alpb   = vals(36) 
    rbornstat = vals(37) 
    ixx    = vals(38)  
    nrespai   = vals(39)  
    iamoeba   = vals(40)
    
    do_amoeba_valence = vals(41)  
	do_amoeba_nonbond = vals(42)  
	do_bond           = vals(43) 
	do_ureyb          = vals(44) 
	do_reg_angle      = vals(45)  
	do_trig_angle     = vals(46)  
    do_opbend         = vals(47)  
    do_torsion        = vals(48) 
    do_pi_torsion     = vals(49)  
	do_strbend        = vals(50)  
	do_torsion_torsion = vals(51) 
	do_str_torsion     = vals(52)
	do_recip           = vals(53)
	do_adjust          = vals(54) 
	do_direct          = vals(55) 
	do_self            = vals(56) 
	do_vdw             = vals(57) 
	do_induced         = vals(58) 
	do_vdw_taper       = vals(59) 
	do_vdw_longrange   = vals(60) 
	beeman_integrator  = vals(61)
	dipole_scf_iter_max = vals(62) 
	amoeba_verbose      = vals(63) 
        
    read(watdef(1), '(A4)') iwtnm
	read(watdef(2), '(A4)') iowtnm
	read(watdef(3), '(A4)') ihwtnm1
	read(watdef(4), '(A4)') ihwtnm2 
	
end subroutine set_mdin_ctrl_int

subroutine set_mdin_ctrl_dbl(vals)
! Set values of mdin_ctrl_dbl  structure
	use mdin_ctrl_dat_mod
	use mdin_amoeba_dat_mod
	
	implicit none
	double precision, intent(in)    :: vals(*)
	
	double precision ::  xx
	
	dielc      = vals(1) 
	es_cutoff  = vals(2)  
	vdw_cutoff = vals(3) 
    scnb       = vals(4)  
    scee       = vals(5) 
    xx         = vals(6) 
    xx         = vals(7) 
    xx         = vals(8) 
    xx         = vals(9) 
    xx         = vals(10) 
    xx         = vals(11) 
    xx         = vals(12) 
    xx         = vals(13) 
    xx         = vals(14) 
    xx         = vals(15) 
    xx         = vals(16) 
    xx         = vals(17) 
    xx         = vals(18) 
    xx         = vals(19) 
    xx         = vals(20)  
    intdiel    = vals(21) 
    extdiel    = vals(22) 
    saltcon    = vals(23) 
    rgbmax     = vals(24) 
    offset     = vals(25) 
    surften    = vals(26) 
    cut_inner  = vals(27) 
    gb_cutoff  = vals(28) 
    gb_alpha   = vals(29) 
    gb_beta    = vals(30)  
    gb_gamma   = vals(31) 
    gb_fs_max  = vals(32) 
    gb_kappa   = vals(33) 
    gb_neckscale = vals(34) 
    arad       = vals(35) 
    bbox_xmin  = vals(36) 
    bbox_ymin  = vals(37) 
    bbox_zmin  = vals(38)  
    bbox_xmax  = vals(39) 
    bbox_ymax  = vals(40) 
    bbox_zmax  = vals(41) 
    
    compress   = vals(42) 
	dipole_scf_tol = vals(43) 
	ee_dsum_cut    = vals(44) 
	ee_damped_cut  = vals(45) 
	sor_coefficient = vals(46) 
	thole_expon_coeff = vals(47)  
	vdw_taper         = vals(48)  
	
end subroutine set_mdin_ctrl_dbl



