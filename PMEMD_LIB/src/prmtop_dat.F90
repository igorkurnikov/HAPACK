#include "copyright.i"

!*******************************************************************************
!
! Module:  prmtop_dat_mod
!
! Description: <TBS>s
!              
!*******************************************************************************

module prmtop_dat_mod

  use gbl_datatypes_mod
  use file_io_dat_mod

  implicit none
  
  integer, parameter    :: max_dihed_dups = 10000

! Global data.  This stuff should be broadcast:

  ! Starting with nttyp, the integer values are derived rather than read.
  ! nttyp is the number of 6-12 vdw parameters.

  integer, parameter    :: prmtop_int_cnt = 27

  integer        ::  natom, ntypes, nbonh, ntheth, nphih, next, nres, &  ! 7
                        nbona, ntheta, nphia, numbnd, numang, nptra, nphb, &      ! 14
                        ifbox_dummy, ifcap_dummy, nspm, numextra_dummy, ncopy_dummy, nttyp, &                    ! 20
                        bonda_idx, anglea_idx, diheda_idx, &                                  ! 23
                        gbl_bond_allocsize, gbl_angle_allocsize, &                          ! 25
                        gbl_dihed_allocsize, num_dist_constr                                ! 27

  common / prmtop_int / natom, ntypes, nbonh, ntheth, nphih, next, nres, &
                        nbona, ntheta, nphia, numbnd, numang, nptra, nphb, &
                        ifbox_dummy, ifcap_dummy, nspm, numextra_dummy, ncopy_dummy, nttyp, &
                        bonda_idx, anglea_idx, diheda_idx, &
                        gbl_bond_allocsize, gbl_angle_allocsize, &
                        gbl_dihed_allocsize, num_dist_constr

  save  :: / prmtop_int /

! Cap data not currently supported...

 integer, parameter    :: prmtop_dbl_cnt = 4

 double precision      cutcap, xcap, ycap, zcap

 common / prmtop_dbl / cutcap, xcap, ycap, zcap


  ! atm_qterm = the atom partial charge array.
  ! atm_iac = TBS
  ! typ_ico = TBS
  ! atm_igraph = the atom name array.
  ! gbl_res_atms = The residue atoms index array.

  ! Atom and residue physical consts, names, indices into parameters,
  ! parameters, etc.

  double precision,     allocatable, save       :: atm_qterm(:)   ! Amoeba
  integer,              allocatable, save       :: atm_iac(:)
  integer,              allocatable, save       :: typ_ico(:)
  double precision,     allocatable, save       :: gbl_cn1(:)
  double precision,     allocatable, save       :: gbl_cn2(:)
  double precision,     allocatable, save       :: gbl_rk(:)
  double precision,     allocatable, save       :: gbl_req(:)
  double precision,     allocatable, save       :: gbl_tk(:)
  double precision,     allocatable, save       :: gbl_teq(:)
  double precision,     allocatable, save       :: gbl_asol(:)
  double precision,     allocatable, save       :: gbl_bsol(:)
  double precision,     allocatable, save       :: gbl_pk(:)
  double precision,     allocatable, save       :: gbl_pn(:)
  double precision,     allocatable, save       :: gbl_phase(:)

  ! Derived from atom paramaters:

  double precision,     allocatable, save       :: gbl_gamc(:)
  double precision,     allocatable, save       :: gbl_gams(:)
  double precision,     allocatable, save       :: gbl_fmn(:) 
  integer,              allocatable, save       :: gbl_ipn(:)

  character(4),         allocatable, save       :: atm_igraph(:)
  integer,              allocatable, save       :: gbl_res_atms(:)

  ! For Generalized Born:

  double precision,     allocatable, save       :: atm_gb_fs(:)
  double precision,     allocatable, save       :: atm_gb_radii(:)

  ! This is the original atom bond, angle, dihedral 
  ! and non-bonded constraints information.
  ! It is used at initialization and when atoms are reassigned to make
  ! the arrays that are actually used in the bonds, angles, and dihedrals
  ! modules.

  type(bond_rec),       allocatable, save       :: gbl_bond(:)
  type(angle_rec),      allocatable, save       :: gbl_angle(:)
  type(dihed_rec),      allocatable, save       :: gbl_dihed(:)
  type(dist_constr_rec),  allocatable, save       :: gbl_dist_constr(:)

  ! NOTE! These are deallocated after setup, unless we are using
  !       Generalized Born:

  ! atm_numex = TBS
  ! gbl_natex = TBS

  integer,              allocatable, save       :: atm_numex(:)
  integer,              allocatable, save       :: gbl_natex(:)

  ! Main title from prmtop.  No need to broadcast, since only master does i/o.
  
  character(80), save   :: prmtop_ititl

! Old prmtop variables that must be printed out but that are otherwise unused:

  integer, private      :: mbona
  integer, private      :: mphia
  integer, private      :: mtheta
  integer, private      :: natyp
  integer, private      :: nhparm
  integer, private      :: nmxrs
  integer, private      :: nparm

! Copies of original values of prmtop variables, the original values of which
! need to be printed out, but which are modified prior to printout:

  integer, private      :: orig_nphih
  integer, private      :: orig_nphia
  
! PRMTOP type definition:

  ! Enumerations for prmtop_type

  integer, parameter                    :: prmtop_type_undefined = 0
  integer, parameter                    :: prmtop_type_nonamoeba = 1
  integer, parameter                    :: prmtop_type_amoeba = 2

  integer, save                         :: prmtop_type = prmtop_type_undefined


! Hide internal routines:

!  private        alloc_amoeba_prmtop_mem

!  private        alloc_prmtop_mem, duplicate_dihedrals, calc_dihedral_parms

contains

!*******************************************************************************
!
! Subroutine:  alloc_prmtop_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_prmtop_mem()

  use mdin_ctrl_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

! Local variables:

  integer                       :: alloc_failed

  ! Generic allocations required for just about everything:

  if( allocated(atm_qterm)) deallocate(atm_qterm)
  if( allocated(atm_iac)) deallocate(atm_iac)
  if( allocated(typ_ico)) deallocate(typ_ico)
  if( allocated(gbl_cn1)) deallocate(gbl_cn1)
  if( allocated(gbl_cn2)) deallocate(gbl_cn2)
  if( allocated(gbl_rk)) deallocate(gbl_rk)
  if( allocated(gbl_req)) deallocate(gbl_req)
  if( allocated(gbl_tk)) deallocate(gbl_tk)
  if( allocated(gbl_teq)) deallocate(gbl_teq)
  if( allocated(gbl_res_atms)) deallocate(gbl_res_atms)
  if( allocated(atm_igraph)) deallocate(atm_igraph)

  allocate(atm_qterm(natom), &
           atm_iac(natom), &
           typ_ico(ntypes * ntypes), &
           gbl_cn1(nttyp), &
           gbl_cn2(nttyp), &
           gbl_rk(numbnd), &
           gbl_req(numbnd), &
           gbl_tk(numang), &
           gbl_teq(numang), &
           gbl_res_atms(nres + 1), &
           atm_igraph(natom), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  if (nphb .ne. 0) then
  
    if( allocated(gbl_asol)) deallocate(gbl_asol)
	if( allocated(gbl_bsol)) deallocate(gbl_bsol)

    allocate(gbl_asol(nphb), &
             gbl_bsol(nphb), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

  end if

  ! Allocate space for atom bond, angle and dihedral information:
  if( allocated(gbl_bond)) deallocate(gbl_bond)
  if( allocated(gbl_angle)) deallocate(gbl_angle)
  if( allocated(gbl_dihed)) deallocate(gbl_dihed)
  if( allocated(gbl_dist_constr)) deallocate(gbl_dist_constr)
  
  allocate(gbl_bond(gbl_bond_allocsize), &
           gbl_angle(gbl_angle_allocsize), &
           gbl_dihed(gbl_dihed_allocsize), &
           gbl_dist_constr(num_dist_constr), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  if (nptra > 0) then
  
    if( allocated(gbl_pk)) deallocate(gbl_pk)
	if( allocated(gbl_pn)) deallocate(gbl_pn)
	if( allocated(gbl_phase)) deallocate(gbl_phase)
	if( allocated(gbl_gamc)) deallocate(gbl_gamc)
	if( allocated(gbl_gams)) deallocate(gbl_gams)
	if( allocated(gbl_fmn)) deallocate(gbl_fmn)
	if( allocated(gbl_ipn)) deallocate(gbl_ipn)
  
    allocate(gbl_pk(nptra), &
             gbl_pn(nptra), &
             gbl_phase(nptra), &
             gbl_gamc(nptra), &
             gbl_gams(nptra), &
             gbl_fmn(nptra), &
             gbl_ipn(nptra), &
             stat = alloc_failed)
             
    if (alloc_failed .ne. 0) call setup_alloc_error

    ! We know these get initialized properly, so no need to zero.

  end if

  ! Allocations for Generalized Born:

  if (using_gb_potential) then

	if( allocated(atm_gb_fs)) deallocate(atm_gb_fs)
	if( allocated(atm_gb_radii)) deallocate(atm_gb_radii)

    allocate(atm_gb_fs(natom), &
             atm_gb_radii(natom), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

  end if

  ! Polarizabilities would need to have space alloc'd if they are supported.

  ! Allocation of stuff that will later be deallocated:

  if( allocated(atm_numex)) deallocate(atm_numex)
  if( allocated(gbl_natex)) deallocate(gbl_natex)

  allocate(atm_numex(natom), &
           gbl_natex(next), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error
  	
  return

end subroutine alloc_prmtop_mem

subroutine alloc_amoeba_prmtop_mem()

!  use pmemd_lib_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

! Local variables:

  integer                       :: alloc_failed

  ! Generic allocations required for just about everything:

  allocate(atm_qterm(natom), &
!           atm_mass(natom), &
           atm_iac(natom), &
           typ_ico(ntypes * ntypes), &
           gbl_res_atms(nres + 1), &
           atm_igraph(natom), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  atm_qterm(:) = 0.d0   ! Safety measure...

  ! Allocation of stuff that will later be deallocated:

  allocate(atm_numex(natom), &
           gbl_natex(next), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  ! Allocation of stuff that will later be deallocated that is only used in
  ! the master:

!  if (master) then
!
!   allocate(atm_isymbl(natom), &
!             atm_itree(natom), &
!             stat = alloc_failed)
!
!    if (alloc_failed .ne. 0) call setup_alloc_error
!
!  end if

  return

end subroutine alloc_amoeba_prmtop_mem
!*******************************************************************************
!
! Subroutine:  calc_dihedral_parms
!
! Description: Routine to calc additional parameters for the vectorised ephi.
!
!*******************************************************************************

subroutine calc_dihedral_parms(numphi, pk, pn, phase, gamc, gams, ipn, fmn)

  use gbl_constants_mod

  implicit none

! Formal arguments:

  integer           numphi
  double precision  pk(*)
  double precision  pn(*)
  double precision  phase(*)
  double precision  gamc(*)
  double precision  gams(*)
  integer           ipn(*)
  double precision  fmn(*)

! Local variables:

  double precision  dum
  double precision  dumc
  double precision  dums
  double precision  four
  integer           i
  double precision  one
  double precision  pim
  double precision  tenm3
  double precision  tenm10
  double precision  zero

  data zero, one, tenm3, tenm10/0.0d+00, 1.0d+00, 1.0d-03, 1.0d-06/
  data four /4.0d+00/

  pim = four*atan(one)

  do i = 1, numphi
    dum = phase(i)
    if (abs(dum - PI) .le. tenm3) dum = sign(pim, dum)
    dumc = cos(dum)
    dums = sin(dum)
    if (abs(dumc) .le. tenm10) dumc = zero
    if (abs(dums) .le. tenm10) dums = zero
    gamc(i) = dumc*pk(i)
    gams(i) = dums*pk(i)
    fmn(i) = one
    if (pn(i) .le. zero) fmn(i) = zero
    pn(i) = abs(pn(i))
    ipn(i) = int(pn(i) + tenm3)
  end do

  return

end subroutine calc_dihedral_parms

end module prmtop_dat_mod

subroutine set_atm_qterm(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(atm_qterm)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(atm_qterm)
	   endif
	   allocate(atm_qterm(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	atm_qterm(1:n) = vals(1:n) 
end subroutine set_atm_qterm

subroutine set_atm_iac(vals,n) 
	
	use prmtop_dat_mod
	
	implicit none
	integer, intent(in)         :: vals(*)
	integer, intent(in)         :: n

	integer          ::   nold
	integer          ::   alloc_failed

	if( n .eq. 0) return
	nold = size(atm_iac)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(atm_iac)
	   endif
	   allocate(atm_iac(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	atm_iac(1:n) = vals(1:n)
	 
end subroutine set_atm_iac

subroutine set_typ_ico(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	integer, intent(in)         :: vals(*)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed

	if( n .eq. 0) return
	nold = size(typ_ico)
	if( nold .ne. n ) then
	    if( allocated(typ_ico)) deallocate(typ_ico)
	   allocate(typ_ico(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	typ_ico(1:n) = vals(1:n)
	 
end subroutine set_typ_ico

subroutine set_gbl_cn1(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_cn1)
	if( nold .ne. n ) then
	    if( allocated(gbl_cn1)) deallocate(gbl_cn1)
	   allocate(gbl_cn1(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_cn1(1:n) = vals(1:n) 
     
end subroutine set_gbl_cn1

subroutine set_gbl_cn2(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_cn2)
	if( nold .ne. n ) then
	    if( allocated(gbl_cn2)) deallocate(gbl_cn2)
	   allocate(gbl_cn2(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_cn2(1:n) = vals(1:n) 
   
end subroutine set_gbl_cn2

subroutine set_gbl_rk(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_rk)
	if( nold .ne. n ) then
	    if( allocated(gbl_rk)) deallocate(gbl_rk)
	   allocate(gbl_rk(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_rk(1:n) = vals(1:n) 

end subroutine set_gbl_rk

subroutine set_gbl_req(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_req)
	if( nold .ne. n ) then
	   if( allocated(gbl_req)) deallocate(gbl_req)
	   allocate(gbl_req(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_req(1:n) = vals(1:n) 
end subroutine set_gbl_req

subroutine set_gbl_tk(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_tk)
	if( nold .ne. n ) then
	   if( allocated(gbl_tk) ) deallocate(gbl_tk)
	   allocate(gbl_tk(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif

	gbl_tk(1:n) = vals(1:n)
	 
end subroutine set_gbl_tk

subroutine set_gbl_teq(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_teq)
	if( nold .ne. n ) then
	   if( allocated(gbl_teq) ) deallocate(gbl_teq)
	   allocate(gbl_teq(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif

	gbl_teq(1:n) = vals(1:n)
	 
end subroutine set_gbl_teq

subroutine set_gbl_asol(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_asol)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(gbl_asol)
	   endif
	   allocate(gbl_asol(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_asol(1:n) = vals(1:n)
	 
end subroutine set_gbl_asol

subroutine set_gbl_bsol(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_bsol)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(gbl_bsol)
	   endif
	   allocate(gbl_bsol(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_bsol(1:n) = vals(1:n)
	 
end subroutine set_gbl_bsol

subroutine set_gbl_pk(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_pk)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(gbl_pk)
	   endif
	   allocate(gbl_pk(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_pk(1:n) = vals(1:n) 
	
end subroutine set_gbl_pk

subroutine set_gbl_pn(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_pn)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(gbl_pn)
	   endif
	   allocate(gbl_pn(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_pn(1:n) = vals(1:n)
	 
end subroutine set_gbl_pn

subroutine set_gbl_phase(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_phase)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(gbl_phase)
	   endif
	   allocate(gbl_phase(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_phase(1:n) = vals(1:n)
	 
end subroutine set_gbl_phase

subroutine set_gbl_gamc(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_gamc)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(gbl_gamc)
	   endif
	   allocate(gbl_gamc(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_gamc(1:n) = vals(1:n)
	 
end subroutine set_gbl_gamc

subroutine set_gbl_gams(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_gams)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(gbl_gams)
	   endif
	   allocate(gbl_gams(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_gams(1:n) = vals(1:n)
	 
end subroutine set_gbl_gams

subroutine set_gbl_fmn(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_fmn)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(gbl_fmn)
	   endif
	   allocate(gbl_fmn(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_fmn(1:n) = vals(1:n)
	 
end subroutine set_gbl_fmn


subroutine int_alloc_arr_set(this,ivals,n)
	implicit none
	integer                 :: this(*)
	integer, intent(in)     :: ivals(*)
	integer, intent(in)     :: n
	
	this(1:n) = ivals(1:n) 
end subroutine int_alloc_arr_set

subroutine set_gbl_ipn(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	integer, intent(in)         :: vals(*)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_ipn)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(gbl_ipn)
	   endif
	   allocate(gbl_ipn(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_ipn(1:n) = vals(1:n)
	 
end subroutine set_gbl_ipn

subroutine set1_atm_igraph(i,vals) 
! Set particular value of atm_igraph Fortran 1-based idx
	use prmtop_dat_mod
	
	implicit none
	character*4, intent(in)    :: vals
	integer, intent(in)         :: i
	
	atm_igraph(i) = vals
end subroutine set1_atm_igraph

subroutine set_gbl_res_atms(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	integer, intent(in)         :: vals(*)
	integer, intent(in)         :: n

	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_res_atms)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(gbl_res_atms)
	   endif
	   allocate(gbl_res_atms(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_res_atms(1:n) = vals(1:n)
	 
end subroutine set_gbl_res_atms

subroutine set_atm_gb_fs(vals,n)
	use prmtop_dat_mod
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(atm_gb_fs)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(atm_gb_fs)
	   endif
	   allocate(atm_gb_fs(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	atm_gb_fs(1:n) = vals(1:n) 
	
end subroutine set_atm_gb_fs

subroutine set_atm_gb_radii(vals,n)
	
	use prmtop_dat_mod
	
	implicit none
	double precision, intent(in)         :: vals(*)
	integer, intent(in)                  :: n

	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(atm_gb_radii)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(atm_gb_radii)
	   endif
	   allocate(atm_gb_radii(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	atm_gb_radii(1:n) = vals(1:n)
	 
end subroutine set_atm_gb_radii

subroutine set1_gbl_bond(i,atm_i_n,atm_j_n,parm_idx_n) 
! Set particular value of gbl_bond Fortran 1-based idx
	use prmtop_dat_mod
	
	implicit none
	integer, intent(in)         :: i
	integer, intent(in)         :: atm_i_n
	integer, intent(in)         :: atm_j_n
	integer, intent(in)         :: parm_idx_n
	
	gbl_bond(i)%atm_i = atm_i_n
	gbl_bond(i)%atm_j = atm_j_n
	gbl_bond(i)%parm_idx = parm_idx_n
	
end subroutine set1_gbl_bond

subroutine set1_gbl_angle(i,atm_i_n,atm_j_n,atm_k_n,parm_idx_n) 
! Set particular value of gbl_angle Fortran 1-based idx
	use prmtop_dat_mod
	
	implicit none
	integer, intent(in)         :: i
	integer, intent(in)         :: atm_i_n
	integer, intent(in)         :: atm_j_n
	integer, intent(in)         :: atm_k_n
	integer, intent(in)         :: parm_idx_n
	
	gbl_angle(i)%atm_i = atm_i_n
	gbl_angle(i)%atm_j = atm_j_n
	gbl_angle(i)%atm_k = atm_k_n
	gbl_angle(i)%parm_idx = parm_idx_n
	
end subroutine set1_gbl_angle

subroutine set1_gbl_dihed(i,atm_i_n,atm_j_n,atm_k_n,atm_l_n,parm_idx_n) 
! Set particular value of gbl_dihed Fortran 1-based idx
	use prmtop_dat_mod
	
	implicit none
	integer, intent(in)         :: i
	integer, intent(in)         :: atm_i_n
	integer, intent(in)         :: atm_j_n
	integer, intent(in)         :: atm_k_n
	integer, intent(in)         :: atm_l_n
	integer, intent(in)         :: parm_idx_n
	
	gbl_dihed(i)%atm_i = atm_i_n
	gbl_dihed(i)%atm_j = atm_j_n
	gbl_dihed(i)%atm_k = atm_k_n
	gbl_dihed(i)%atm_l = atm_l_n
	gbl_dihed(i)%parm_idx = parm_idx_n
	
end subroutine set1_gbl_dihed

subroutine set1_gbl_dist_constr(i,atm_i_n,atm_j_n,acoef_n,bcoef_n,cnt_type_n) 
! Set particular value of gbl_dist_constr Fortran 1-based idx
	use prmtop_dat_mod
	
	implicit none
	integer, intent(in)         :: i
	integer, intent(in)         :: atm_i_n
	integer, intent(in)         :: atm_j_n
	double precision, intent(in):: acoef_n
	double precision, intent(in):: bcoef_n
	integer, intent(in)         :: cnt_type_n
	
!	write(*,*) "  "
!	write(*,*) " i=  ", i
!	write(*,*) " set1_gbl_dist_constr() atm_i_n = ", atm_i_n
!	write(*,*) " set1_gbl_dist_constr() atm_j_n = ", atm_j_n
!	write(*,*) " set1_gbl_dist_constr() acoef_n = ", acoef_n
!	write(*,*) " set1_gbl_dist_constr() bcoef_n = ", bcoef_n
!	write(*,*) " set1_gbl_dist_constr() cnt_type_n = ", cnt_type_n
!	write(*,*) "  "
	 
	gbl_dist_constr(i)%atm_i = atm_i_n
	gbl_dist_constr(i)%atm_j = atm_j_n
	gbl_dist_constr(i)%acoef = acoef_n
	gbl_dist_constr(i)%bcoef = bcoef_n
	gbl_dist_constr(i)%cnt_type = cnt_type_n
	
end subroutine set1_gbl_dist_constr

subroutine set_atm_numex(vals,n)
	use prmtop_dat_mod
	implicit none
	integer, intent(in)         :: vals(*)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(atm_numex)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(atm_numex)
	   endif
	   allocate(atm_numex(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	atm_numex(1:n) = vals(1:n)
	 
end subroutine set_atm_numex

subroutine set_gbl_natex(vals,n)
	use prmtop_dat_mod
	implicit none
	integer, intent(in)         :: vals(*)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	if( n .eq. 0) return
	nold = size(gbl_natex)
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(gbl_natex)
	   endif
	   allocate(gbl_natex(n), stat = alloc_failed)
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	gbl_natex(1:n) = vals(1:n)
	 
end subroutine set_gbl_natex

subroutine set_prmtop_ititl(vals) 
! Set value of prmtop title
	use prmtop_dat_mod
	
	implicit none
	character(80), intent(in)    :: vals
	integer i;
	
	prmtop_ititl(1:80) = vals(1:80)
	
end subroutine set_prmtop_ititl

subroutine set_prmtop_int(vals)
! Set values of /prmtop_int/ common	
	use prmtop_dat_mod
	
	implicit none
	integer, intent(in)    :: vals(*)

	natom   = vals(1)
	ntypes  = vals(2)
	nbonh   = vals(3)
	ntheth  = vals(4)
	nphih   = vals(5)
	next    = vals(6)
	nres    = vals(7)
    nbona   = vals(8)
    ntheta  = vals(9)
    nphia   = vals(10)
    numbnd  = vals(11)
    numang  = vals(12)
    nptra   = vals(13)
    nphb    = vals(14)
    ifbox_dummy   = vals(15)
    ifcap_dummy   = vals(16)
    nspm    = vals(17)
    numextra_dummy = vals(18)
    ncopy_dummy    = vals(19)
    nttyp   = vals(20)
    bonda_idx = vals(21)
    anglea_idx = vals(22)
    diheda_idx = vals(23)
    gbl_bond_allocsize = vals(24)
    gbl_angle_allocsize = vals(25)
    gbl_dihed_allocsize = vals(26)
    num_dist_constr   = vals(27)
    
end subroutine set_prmtop_int


