#include "copyright.i"

!*******************************************************************************
!
! Module:  angles_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module angles_mod

  use gbl_datatypes_mod

  implicit none

! The following are derived from prmtop angle info:

  integer, save                         :: cit_ntheth, cit_ntheta

  type(angle_rec), allocatable, save    :: cit_angle(:)

contains

!*******************************************************************************
!
! Subroutine:  angles_setup
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine angles_setup(use_atm_map,atm_owner_map)

  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer                       :: use_atm_map(natom)
  integer                       :: atm_owner_map(natom)

! Local variables:

  integer               :: alloc_failed
  type(angle_rec)       :: angles_copy(ntheth + ntheta)
  integer               :: atm_i, atm_j, atm_k, angles_idx
  integer               :: my_angle_cnt

! This routine can handle reallocation, and thus can be called multiple
! times.

! Find all angles for which this process owns either atom:

  my_angle_cnt = 0

  do angles_idx = 1, ntheth

    atm_i = gbl_angle(angles_idx)%atm_i
    atm_j = gbl_angle(angles_idx)%atm_j
    atm_k = gbl_angle(angles_idx)%atm_k

    if (atm_owner_map(atm_i) .eq. mytaskid) then
      my_angle_cnt = my_angle_cnt + 1
      angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
    else if (atm_owner_map(atm_j) .eq. mytaskid) then
      my_angle_cnt = my_angle_cnt + 1
      angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
    else if (atm_owner_map(atm_k) .eq. mytaskid) then
      my_angle_cnt = my_angle_cnt + 1
      angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
    end if

  end do

  cit_ntheth = my_angle_cnt

  do angles_idx = anglea_idx, anglea_idx + ntheta - 1

    atm_i = gbl_angle(angles_idx)%atm_i
    atm_j = gbl_angle(angles_idx)%atm_j
    atm_k = gbl_angle(angles_idx)%atm_k

    if (atm_owner_map(atm_i) .eq. mytaskid) then
      my_angle_cnt = my_angle_cnt + 1
      angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
    else if (atm_owner_map(atm_j) .eq. mytaskid) then
      my_angle_cnt = my_angle_cnt + 1
      angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
    else if (atm_owner_map(atm_k) .eq. mytaskid) then
      my_angle_cnt = my_angle_cnt + 1
      angles_copy(my_angle_cnt) = gbl_angle(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
    end if

  end do

  cit_ntheta = my_angle_cnt - cit_ntheth

  if (my_angle_cnt .gt. 0) then
    if (allocated(cit_angle)) then
      if (size(cit_angle) .lt. my_angle_cnt) then
        deallocate(cit_angle)
        allocate(cit_angle(my_angle_cnt), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
      end if
    else
      allocate(cit_angle(my_angle_cnt), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
    end if
    cit_angle(1:my_angle_cnt) = angles_copy(1:my_angle_cnt)
  end if

!BEGIN DBG
! write(0,*)'task,cit_ntheth,cit_ntheta,', mytaskid, &
!            cit_ntheth, cit_ntheta
! write(0,*)'task,    ntheth,    ntheta', mytaskid, &
!                ntheth,     ntheta
!END DBG
  return

end subroutine angles_setup

!*******************************************************************************
!
! Subroutine:  get_angle_energy
!
! Description:  Routine to get the bond energies and forces for potentials of
!               the type ct*(t-t0)**2.
!
!*******************************************************************************

subroutine get_angle_energy(angle_cnt, angle, x, frc, angle_energy, atm_owner_map)

  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: angle_cnt
  type(angle_rec)       :: angle(*)
  double precision      :: x(3, *)
  double precision      :: frc(3, *)
  double precision      :: angle_energy
  integer               :: atm_owner_map(*)

! Local variables:

  double precision, parameter   :: pt999 = 0.9990d0

  double precision      :: ant
  double precision      :: cst
  double precision      :: dfw
  double precision      :: eaw
  double precision      :: rij, rik, rkj
  double precision      :: xij, yij, zij
  double precision      :: xkj, ykj, zkj
  double precision      :: cii, cik, ckk
  double precision      :: da
  double precision      :: df
  double precision      :: dt1, dt2, dt3, dt4, dt5, dt6
  double precision      :: sth
  integer               :: i, j, k, ic, jn

  angle_energy = 0.0d0

! Grand loop for the angle stuff:

  do jn = 1, angle_cnt

    i = angle(jn)%atm_i
    j = angle(jn)%atm_j
    k = angle(jn)%atm_k
    ic = angle(jn)%parm_idx

! Calculation of the angle:

    xij = x(1, i) - x(1, j)
    xkj = x(1, k) - x(1, j)

    yij = x(2, i) - x(2, j)
    ykj = x(2, k) - x(2, j)

    zij = x(3, i) - x(3, j)
    zkj = x(3, k) - x(3, j)

    rij = xij * xij + yij * yij + zij * zij
    rkj = xkj * xkj + ykj * ykj + zkj * zkj

    rik = sqrt(rij * rkj)

    cst = min(pt999, max(-pt999, (xij * xkj + yij * ykj + zij * zkj) / rik))

    ant = acos(cst)

! Calculation of the energy and deriv:

    da = ant - gbl_teq(ic)

    df = gbl_tk(ic) * da
    eaw = df * da
    dfw = -(df + df) / sin(ant)

! Calculation of the force: 

    cik = dfw / rik
    sth = dfw * cst
    cii = sth / rij
    ckk = sth / rkj

    dt1 = cik * xkj - cii * xij
    dt2 = cik * ykj - cii * yij
    dt3 = cik * zkj - cii * zij
    dt4 = cik * xij - ckk * xkj
    dt5 = cik * yij - ckk * ykj
    dt6 = cik * zij - ckk * zkj
    
    ! We use atm_i to determine who sums up the energy...
    if (atm_owner_map(i) .eq. mytaskid) then

      angle_energy = angle_energy + eaw

      frc(1, i) = frc(1, i) - dt1
      frc(2, i) = frc(2, i) - dt2
      frc(3, i) = frc(3, i) - dt3

    end if

    if (atm_owner_map(j) .eq. mytaskid) then

      frc(1, j) = frc(1, j) + dt1 + dt4
      frc(2, j) = frc(2, j) + dt2 + dt5
      frc(3, j) = frc(3, j) + dt3 + dt6

    end if

    if (atm_owner_map(k) .eq. mytaskid) then

      frc(1, k) = frc(1, k) - dt4
      frc(2, k) = frc(2, k) - dt5
      frc(3, k) = frc(3, k) - dt6

    end if
        
  end do

  return

end subroutine get_angle_energy

end module angles_mod
