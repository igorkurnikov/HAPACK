#include "copyright.i"

!*******************************************************************************
!
! Module: shake_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module shake_mod

  implicit none

  type shake_bond_rec
    integer             :: atm_i
    integer             :: atm_j 
    double precision    :: parm
  end type shake_bond_rec

  ! Note: fastwat_res_lst and my_fastwat_res_lst are actually lists of the
  ! first atm_id in each water residue and each water residue processed by
  ! this processor, respectively.

  integer, allocatable, save, private   :: nonfastwat_bond_lst(:)    ! used only in MPI mode

  integer, save, private        :: nonfastwat_bond_cnt

  type(shake_bond_rec), allocatable, save, private :: my_nonfastwat_bond_dat(:)
  integer, save, private                           :: my_nonfastwat_bond_cnt

  double precision, save, private                  :: box_half(3)

  private       wrap_crds

contains

!*******************************************************************************
!
! Subroutine:   claim_my_nonfastwat_bonds
!
! Description:  <TBS>
!               Used only in MPI
!*******************************************************************************

subroutine claim_my_nonfastwat_bonds(atm_owner_map)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none
  
  integer               :: atm_owner_map(*)

! Local variables:

  integer               :: alloc_failed
  integer               :: atm_i, atm_j
  integer               :: bond_idx        
  integer               :: bonddat_idx        
  integer               :: bondlst_idx

  ! Determine how many nonfastwater bonds you actually own.

  my_nonfastwat_bond_cnt = 0

  do bondlst_idx = 1, nonfastwat_bond_cnt

    bond_idx = nonfastwat_bond_lst(bondlst_idx)

    atm_i = gbl_bond(bond_idx)%atm_i
    atm_j = gbl_bond(bond_idx)%atm_j

    if (atm_owner_map(atm_i) .eq. mytaskid) then

      my_nonfastwat_bond_cnt = my_nonfastwat_bond_cnt + 1
          
      if (atm_owner_map(atm_j) .ne. mytaskid) then
        write(error_msg, *) 'Partition error in shake, task ', mytaskid
        call mol_mech_error
      end if

    else

      if (atm_owner_map(atm_j) .eq. mytaskid) then
        write(error_msg, *) 'Partition error in shake, task ', mytaskid
        call mol_mech_error
      end if

    end if

  end do

  ! Allocate space. Note it is reallocatable.

  if (allocated(my_nonfastwat_bond_dat)) then

    if (size(my_nonfastwat_bond_dat) .lt. my_nonfastwat_bond_cnt) then

      deallocate(my_nonfastwat_bond_dat)
      allocate(my_nonfastwat_bond_dat(my_nonfastwat_bond_cnt), &
               stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error

    end if

  else

    allocate(my_nonfastwat_bond_dat(my_nonfastwat_bond_cnt), &
             stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error

  end if

  bonddat_idx = 0

  ! Set up your bond data.

  do bondlst_idx = 1, nonfastwat_bond_cnt

    bond_idx = nonfastwat_bond_lst(bondlst_idx)

    atm_i = gbl_bond(bond_idx)%atm_i

    if (atm_owner_map(atm_i) .eq. mytaskid) then

      bonddat_idx = bonddat_idx + 1

      my_nonfastwat_bond_dat(bonddat_idx)%atm_i = gbl_bond(bond_idx)%atm_i
      my_nonfastwat_bond_dat(bonddat_idx)%atm_j = gbl_bond(bond_idx)%atm_j
      my_nonfastwat_bond_dat(bonddat_idx)%parm = &
        gbl_req(gbl_bond(bond_idx)%parm_idx)**2 

    end if

  end do

  return

end subroutine claim_my_nonfastwat_bonds

!*******************************************************************************
!
! Subroutine:   shake_setup
!
! Description:
!
!*******************************************************************************

subroutine shake_setup(atm_owner_map)

  use constraints_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none
  
  integer               :: atm_owner_map(*)

! Local variables:

  integer               :: alloc_failed
  integer               :: bond_cnt
  integer               :: bond_idx        
  integer               :: bonddat_idx        
  integer               :: bondlst_idx

  if (ntc .eq. 2) then

    bond_cnt = nbonh

  else if (ntc .eq. 3) then

    if (master .and. numtasks .gt. 1) write(error_msg, *) 'this parallel version only works for ntc < 3'
    call mol_mech_error

    bond_cnt = nbonh + nbona

  else

    ! shake not done on anything.

    nonfastwat_bond_cnt = 0
    my_nonfastwat_bond_cnt = 0

    return

  end if

  ! Determine how many nonfastwater bonds there are.

  nonfastwat_bond_cnt = bond_cnt

  ! If you need it, set up the nonfastwater bonds list.

  if (numtasks .gt. 1) then         ! Molecular Dynamics

    ! Set up the nonfastwater bonds list.  This is only used for MPI
    ! molecular dynamics, not minimizations.

	if( allocated(nonfastwat_bond_lst)) deallocate(nonfastwat_bond_lst)

    allocate(nonfastwat_bond_lst(nonfastwat_bond_cnt), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error

    bondlst_idx = 0

    do bond_idx = 1, bond_cnt
        bondlst_idx = bondlst_idx + 1
        nonfastwat_bond_lst(bondlst_idx) = bond_idx
    end do

  end if

  ! Claim all the nonfastwater bonds.  This is not true for mpi molecular
  ! dynamics, and will be corrected in claim_nonfastwater_bonds.

  my_nonfastwat_bond_cnt = nonfastwat_bond_cnt

  if ( numtasks .gt. 1) then     

     call  claim_my_nonfastwat_bonds(atm_owner_map)
  
  else

    ! Allocate space.

	if( allocated(my_nonfastwat_bond_dat)) deallocate(my_nonfastwat_bond_dat)

    allocate(my_nonfastwat_bond_dat(nonfastwat_bond_cnt), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error

    bonddat_idx = 0

    do bond_idx = 1, bond_cnt
        bonddat_idx = bonddat_idx + 1
        my_nonfastwat_bond_dat(bonddat_idx)%atm_i = gbl_bond(bond_idx)%atm_i
        my_nonfastwat_bond_dat(bonddat_idx)%atm_j = gbl_bond(bond_idx)%atm_j
        my_nonfastwat_bond_dat(bonddat_idx)%parm = &
        gbl_req(gbl_bond(bond_idx)%parm_idx)**2
    end do
    
  end if ! if( numtasks .gt. 1) 

  return

end subroutine shake_setup

!*******************************************************************************
!
! Subroutine:  shake
!
! Description:
!
! All the bonds involving hydrogen atoms are loaded first in struct gbl_bond
! followed by those involving non-hydrogen atoms and the perturbed atoms.
!
! Mods for version 4.1:
!   - Add fastwat_bonds(i), so that waters shaken by fast analytic 3-point
!     shake will not be shake-n here -- dap
!
!
! The only MPI code specific to this routine is wrappers for
! all the I/O to allow only the master to write.
!
!
!*******************************************************************************

subroutine shake(x, xp, ntb_par, igroup, box, mass_inv,tol_lcl)

  use constraints_mod
  use dynamics_dat_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pbc_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  double precision       :: x(3, *)
  double precision       :: xp(3, *)
  integer                :: ntb_par
  integer                :: igroup(*)
  double precision       :: box(3)
  double precision       :: mass_inv(*)
  double precision       :: tol_lcl
  

! Local variables:

  double precision       :: acor
  double precision       :: box_half_min
  double precision       :: box_half_min_sq
  double precision       :: diff
  integer                :: i, j, m
  integer                :: iter_cnt
  integer                :: my_bond_idx
  integer                :: ns
  integer                :: ibelly_lcl
  logical                :: done

  double precision       :: rpij2
  double precision       :: rrpr
  
  double precision       :: toler
  double precision       :: winvi, winvj
  double precision       :: xh
  double precision       :: xij(3)
  double precision       :: xpij(3)

  double precision, save :: zero = 0.0d0

  integer                :: skips(natom)

  ! Calc size of max vector that is guaranteed to fit in box w/out imaging:

  box_half(:) = box(:) * 0.5d0
  box_half_min = min(box_half(1), box_half(2), box_half(3))

  box_half_min_sq = box_half_min * box_half_min

  ibelly_lcl = ibelly

  skips(:) = 1

  do iter_cnt = 1, 3000

    done = .true.       ! until proven otherwise.

! Loop over all the bonds that are not 3-point waters:

    do my_bond_idx = 1, my_nonfastwat_bond_cnt

      i = my_nonfastwat_bond_dat(my_bond_idx)%atm_i
      j = my_nonfastwat_bond_dat(my_bond_idx)%atm_j

      if (skips(i) .lt. iter_cnt .and. skips(j) .lt. iter_cnt) cycle

! Calc nominal distance squared:

      rpij2 = zero

      do  m = 1, 3
        xpij(m) = xp(m, i) - xp(m, j)
        rpij2 = rpij2 + xpij(m)**2
      end do

! BUGBUG - The following boundary check is no longer used in sander 8+, and
!          is probably unnecessary since we are dealing with atom coords (not
!          image crds) here and this stuff is bonded...  We leave it in for
!          now since it is probably not doing any harm in terms of results and
!          minimal harm in terms of performance...

! If boundary condition is not present skip it:

      ns = 0 

      if (ntb_par .gt. 0) then

        if (rpij2 .ge. box_half_min_sq) then

! Apply the boundary & recalc the distance squared:

          ns = 1
          call wrap_crds(xpij,box)
          rpij2 = zero

          do m = 1, 3
            rpij2 = rpij2 + xpij(m)**2
          end do

        end if

      end if

! Apply the correction:

      toler = my_nonfastwat_bond_dat(my_bond_idx)%parm

      diff = toler - rpij2

      if (abs(diff) .lt. toler * tol_lcl) cycle

      do  m = 1, 3
        xij(m) = x(m, i) - x(m, j)
      end do

      if (ns .ne. 0) then
        call wrap_crds(xij,box)
      end if

! Shake resetting of coordinate is done here:

      rrpr = zero

      do  m = 1, 3
        rrpr = rrpr + xij(m) * xpij(m)
      end do

      if (rrpr .lt. toler * 1.0d-06) then ! Deviation too large.  Kill PMEMD.
        if (master) then
          write(error_msg, 321) iter_cnt, my_bond_idx, i, j
        end if
        call mol_mech_error
      end if

      winvi = mass_inv(i)
      winvj = mass_inv(j)

! If belly option is on then resetting of the frozen atom is to be prevented:

      if (ibelly_lcl .gt. 0) then
        if (igroup(i) .le. 0) winvi = zero
        if (igroup(j) .le. 0) winvj = zero
      end if

      acor = diff / (rrpr * (mass_inv(i) + mass_inv(j) + &
             mass_inv(i) + mass_inv(j)))

      do m = 1, 3
        xh = xij(m) * acor
        xp(m, i) = xp(m, i) + xh * winvi
        xp(m, j) = xp(m, j) - xh * winvj
      end do

      skips(i) = iter_cnt + 1 ! Process in this and the next iteration.
      skips(j) = iter_cnt + 1 ! Process in this and the next iteration.

      done = .false.

    end do

    if (done) return    ! Normal exit

  end do

  ! We failed to accomplish coordinate resetting.  Kill PMEMD.

  if (master) write(mdout, 311)

  call mol_mech_error

311 format(/5x, 'Coordinate resetting (shake) was not accomplished', &
           /5x, 'within 3000 iterations')
321 format(/5x, 'Coordinate resetting cannot be accomplished,', &
           /5x, 'deviation is too large', &
           /5x, 'iter_cnt, my_bond_idx, i and j are :', 7i5)

end subroutine shake

!*******************************************************************************
!
! Internal Subroutine:  wrap_crds
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine wrap_crds(xpij,box)

  use pbc_mod

  implicit none

! Formal arguments:

  double precision  xpij(3)
  double precision  box(3)

! Local variables:

  integer           m

  do m = 1, 3

    if (xpij(m) .ge. box_half(m)) then

      xpij(m) = xpij(m) - box(m)

    else if (xpij(m) .lt. -box_half(m)) then

      xpij(m) = xpij(m) + box(m)

    end if

  end do

  return

end subroutine wrap_crds

end module shake_mod
