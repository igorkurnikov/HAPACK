#include "copyright.i"

!*******************************************************************************
!
! Module:  constraints_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module constraints_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  get_crd_constraint_energy
!
! Description: Routine to put harmonic constraints for position.
!
! Mods for Rev A by GLS.
!
!*******************************************************************************

subroutine get_crd_constraint_energy(natc, econ, jrc, x, frc, xc, weit, atm_owner_map )

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: natc
  double precision      :: econ
  integer               :: jrc(*)
  double precision      :: x(3, *)
  double precision      :: frc(3, *)
  double precision      :: xc(3, *)
  double precision      :: weit(*)
  integer               :: atm_owner_map(*) 

! Local variables:

  double precision      :: ax, ay, az
  double precision      :: eadd
  integer               :: i, j
  double precision      :: wt
  double precision      :: wx, wy, wz

  econ = 0.0d+00

! BUGBUG - A more efficient implementation would modify jrc to contain
!          only atoms owned by this processor, but this constraint stuff
!          is probably not that much used...

  do j = 1, natc
    i = jrc(j)
    if (atm_owner_map(i) .eq. mytaskid) then
      wt = weit(j)
      ax = x(1, i) - xc(1, i)
      ay = x(2, i) - xc(2, i)
      az = x(3, i) - xc(3, i)
      wx = wt * ax
      wy = wt * ay
      wz = wt * az
      eadd = wx * ax + wy * ay + wz * az
      econ = econ + eadd
      frc(1, i) = frc(1, i) - (wx + wx)
      frc(2, i) = frc(2, i) - (wy + wy)
      frc(3, i) = frc(3, i) - (wz + wz)
    end if
  end do

  return

end subroutine get_crd_constraint_energy

end module constraints_mod
