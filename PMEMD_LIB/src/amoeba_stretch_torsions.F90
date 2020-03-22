#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_stretch_torsions_mod
!
! Description: <TBS>
!
!*******************************************************************************

! BUGBUG - Current test cases don't include coverage for this module...

module amoeba_stretch_torsions_mod

  implicit none

  private

  integer, save :: do_amoeba_stretch_torsions_flag, num_list, num_params

  integer, allocatable, save            :: list(:,:)
  double precision, allocatable, save   :: force_constant(:)
  double precision, allocatable, save   :: periodicity(:)
  double precision, allocatable, save   :: cosphase(:)
  double precision, allocatable, save   :: sinphase(:)
  double precision, allocatable, save   :: bond_equil_value(:)

  ! BUGBUG - move to local storage:

  double precision, save                :: energy
  double precision, save                :: virial(3,3)

  double precision, parameter           :: stretch_tor_unit = 1.d0

  public        am_stretch_torsions_zero_flag
  public        am_stretch_torsions_suser_bit
  public        am_stretch_torsions_eval

contains

!*******************************************************************************!
! Subroutine:  am_stretch_torsions_zero_flag
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_stretch_torsions_zero_flag
  implicit none
  do_amoeba_stretch_torsions_flag = 0
  return
end subroutine am_stretch_torsions_zero_flag

!*******************************************************************************!
! Subroutine:  am_stretch_torsions_suser_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_stretch_torsions_suser_bit(do_this)

  use amoeba_flags_mod

  implicit none

! Formal arguments:

  integer, intent(in)   :: do_this

  call set_user_bit(do_this, do_amoeba_stretch_torsions_flag)

  return

end subroutine am_stretch_torsions_suser_bit

!*******************************************************************************!
! Subroutine:  am_stretch_torsions_eval
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_stretch_torsions_eval(crd, frc, ene, vir)

  use amoeba_flags_mod

  implicit none

! Formal arguments:

  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: frc(3, *)
  double precision, intent(out)         :: ene
  double precision, intent(in out)      :: vir(3, 3)

! Local variables:

  double precision                      :: fn(num_list)
  double precision                      :: dfn_darg(num_list)
  double precision                      :: cosarg(num_list)
  double precision                      :: sinarg(num_list)
  double precision                      :: darg1_dcrd(12, num_list)
  double precision                      :: darg2_dcrd(6, num_list)
  double precision                      :: arg2(num_list)

  energy = 0.d0
  virial(:,:) = 0.d0

  if (do_amoeba_stretch_torsions_flag .ne. proceed) return

  call am_stretch_torsions_get_tor_arg(crd, list, cosarg, sinarg, darg1_dcrd)

  call am_stretch_torsions_bond_arg(crd, list, bond_equil_value, &
                                    arg2, darg2_dcrd)

  call am_stretch_torsions_func(list, force_constant, periodicity, cosphase, &
                                sinphase, cosarg, sinarg, fn, dfn_darg)

  call am_stretch_torsions_get_ene_frc(list, fn, dfn_darg, darg1_dcrd, &
                                       arg2, darg2_dcrd, crd, frc)
  ene = energy
  vir(:,:) = vir(:, :) + virial(:, :)

  return

end subroutine am_stretch_torsions_eval

!*******************************************************************************!
! Subroutine:  am_stretch_torsions_get_tor_arg
!
! Description:
!
! This routine calculates torsion function argument and its derivatives with
! respect to atomic positions of atom i, j, k, l. 
!
! INPUT variables:
!    nangles:  number of angles in list
!    crd the atomic coord array
!    list: 4 x ntorsions array giving for each torsion the index of the first
!           atom, index of the second atom, index of the third, fourth,
!           and index into the torsion parameter table giving the 
!           18 terms of 6th order fourier expansion
! OUTPUT variables:
!    cosarg, array of cosines of phi angles 
!    sinarg, array of sines of phi angles 
!    darg_dcrdijkl   derivs of arg wrt crds of atom i, j, k, l
!
!*******************************************************************************

subroutine am_stretch_torsions_get_tor_arg(crd, list, cosarg, sinarg, &
                                           darg_dcrdijkl)

  implicit none

! Formal arguments:

  double precision, intent(in)  :: crd(3, *)
  integer, intent(in)           :: list(5, *)
  double precision, intent(out) :: cosarg(*)
  double precision, intent(out) :: sinarg(*)
  double precision, intent(out) :: darg_dcrdijkl(12, *)

! Local variables:

  integer                       :: i, j, k, l, n, m
  double precision              :: crd_abcd(12)
  double precision              :: gradphi_abcd(12)
  double precision              :: cosphi
  double precision              :: sinphi

  do n = 1, num_list

    i = list(1, n)
    j = list(2, n)
    k = list(3, n)
    l = list(4, n)

    do m = 1, 3
      crd_abcd(m) = crd(m, i)
      crd_abcd(m + 3) = crd(m, j)
      crd_abcd(m + 6) = crd(m, k)
      crd_abcd(m + 9) = crd(m, l)
    end do

    call am_val_geom_torsion(crd_abcd, gradphi_abcd, cosphi, sinphi)

    cosarg(n) = cosphi
    sinarg(n) = sinphi

    do m = 1, 12
      darg_dcrdijkl(m, n) = gradphi_abcd(m)
    end do

  end do

  return

end subroutine am_stretch_torsions_get_tor_arg

!*******************************************************************************!
! Subroutine:  am_stretch_torsions_bond_arg
!
! Description:
!
! This routine calculates bond function argument and its derivatives with
! respect to atomic positions of atom j, k. (tor list atoms are i, j, k, l
! This gets used in a stretch-torsion energy routine 
! Note that derivatives with respect 
! to atom k are just the negative of those with respect to atom i
! NOTE that periodic boundary conditions are not used here
!      i.e. imaging is done on a per molecule basis
!
! INPUT variables:
!    crd the atomic coord array
!    list: 3 x nbond array giving for each bond the index of the first
!           atom, index of the second atom, and index into the bond
!           parameter table giving the equilibrium length
!    equil_value the list of ref bond lengths
! OUTPUT variables:
!    arg, array of bond function args
!    darg_dcrd   derivs of arg wrt crds of atom j, k
!
!*******************************************************************************

subroutine am_stretch_torsions_bond_arg(crd, list, equil_value, arg, darg_dcrd)

  implicit none

! Formal arguments:

  double precision, intent(in)  :: crd(3, *)
  integer, intent(in)           :: list(5, *)
  double precision, intent(in)  :: equil_value(*)
  double precision, intent(out) :: arg(*)
  double precision, intent(out) :: darg_dcrd(6, *)

! Local variables:

  integer                       :: j, k, n, it
  double precision              :: bl, dx, dy, dz

  do n = 1, num_list

    j = list(2, n)
    k = list(3, n)
    it = list(5, n)
    dx = crd(1, j) - crd(1, k)
    dy = crd(2, j) - crd(2, k)
    dz = crd(3, j) - crd(3, k)
    bl = sqrt(dx * dx + dy * dy + dz * dz)

! Recall that the reference bond length is 2nd parameter.

    arg(n) = bl - equil_value(it)

! Differentiate bl to get darg_dcrd.

    darg_dcrd(1, n) = dx / bl
    darg_dcrd(2, n) = dy / bl
    darg_dcrd(3, n) = dz / bl
    darg_dcrd(4, n) = -dx / bl
    darg_dcrd(5, n) = -dy / bl
    darg_dcrd(6, n) = -dz / bl

  end do

  return

end subroutine am_stretch_torsions_bond_arg

!*******************************************************************************!
! Subroutine:  am_stretch_torsions_func
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_stretch_torsions_func(list, force_constant, periodicity, &
                                    cosphase, sinphase, cosarg, sinarg, &
                                    func, dfunc_darg)
  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: list(5, *)
  double precision, intent(in)  :: force_constant(*)
  double precision, intent(in)  :: periodicity(*)
  double precision, intent(in)  :: cosphase(*)
  double precision, intent(in)  :: sinphase(*)
  double precision, intent(in)  :: cosarg(*)
  double precision, intent(in)  :: sinarg(*)
  double precision, intent(out) :: func(*)
  double precision, intent(out) :: dfunc_darg(*)

! Local variables:

  integer                       :: n, m, it, degree
  double precision              :: cosine(10), sine(10) ! max of degree is 10
  double precision              :: amplitude, period
  double precision              :: cos_phase, sin_phase

  do n = 1, num_list

    cosine(1) = cosarg(n)
    sine(1) = sinarg(n)
    it = list(5, n)
    amplitude = force_constant(it)
    period = periodicity(it)
    cos_phase = cosphase(it)
    sin_phase = sinphase(it)
    degree = nint(period)

    if (degree .eq. 2) then
      cosine(2) = cosine(1) * cosine(1) - sine(1) * sine(1)
      sine(2) = sine(1) * cosine(1) + cosine(1) * sine(1)
    else if (degree .eq. 3) then
      cosine(2) = cosine(1) * cosine(1) - sine(1) * sine(1)
      sine(2) = sine(1) * cosine(1) + cosine(1) * sine(1)
      cosine(3) = cosine(2) * cosine(1) - sine(2) * sine(1)
      sine(3) = sine(2) * cosine(1) + cosine(2) * sine(1)
    else if (degree .gt. 3) then
      if (degree .gt. 10) then
        write(error_msg, *)'AM_TORSIONS_func: degree too big: ', degree
        call mol_mech_error
      end if
      do m = 2, degree
        cosine(m) = cosine(m - 1) * cosine(1) - sine(m - 1) * sine(1)
        sine(m) = sine(m - 1) * cosine(1) + cosine(m - 1) * sine(1)
      end do
    end if

    func(n) =  amplitude * (1.d0 + cos_phase * cosine(degree) + &
                                     sin_phase * sine(degree))
    dfunc_darg(n) = amplitude * period * (sin_phase * cosine(degree) - &
                                            cos_phase * sine(degree))
  end do 

  return

end subroutine am_stretch_torsions_func

!*******************************************************************************!
! Subroutine:  am_stretch_torsions_get_ene_frc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_stretch_torsions_get_ene_frc(list, fn, dfn_darg1, darg1_dcrd, &
                                           arg2, darg2_dcrd, crd, frc)

  implicit none

! Formal arguments:

  integer, intent(in)                   :: list(5, *) ! 4 atoms plus param ptr
  double precision, intent(in)          :: fn(*)
  double precision, intent(in)          :: dfn_darg1(*)
  double precision, intent(in)          :: darg1_dcrd(12, *)
  double precision, intent(in)          :: arg2(*)
  double precision, intent(in)          :: darg2_dcrd(6, *)
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: frc(3, *)

! Local variables:

  integer                               :: i, j, k, l, n, m, p, q
  double precision                      :: term, f(12)

  do n = 1, num_list

    i = list(1, n)
    j = list(2, n)
    k = list(3, n)
    l = list(4, n)
    energy = energy + stretch_tor_unit * arg2(n) * fn(n)

! Apply product rule and then chain rule to get deriv of energy with respect 
! to crds of i, j, k, l.
! dfn_darg1 holds the deriv of fn with respect to its arg
! while darg1_dcrd holds the derivs of arg with respect to crds of i, j, k, l.
! Recall force is negative of grad.

    term = stretch_tor_unit * arg2(n) * dfn_darg1(n)

    do m = 1, 12
      f(m) = term * darg1_dcrd(m, n)
    end do

    ! 2nd term from product rule

    term = stretch_tor_unit * fn(n)

    do m = 1, 6
      f(m + 3) = f(m + 3) + term * darg2_dcrd(m, n)
    end do

    frc(1, i) = frc(1, i) - f(1)
    frc(2, i) = frc(2, i) - f(2)
    frc(3, i) = frc(3, i) - f(3)
    frc(1, j) = frc(1, j) - f(4)
    frc(2, j) = frc(2, j) - f(5)
    frc(3, j) = frc(3, j) - f(6)
    frc(1, k) = frc(1, k) - f(7)
    frc(2, k) = frc(2, k) - f(8)
    frc(3, k) = frc(3, k) - f(9)
    frc(1, l) = frc(1, l) - f(10)
    frc(2, l) = frc(2, l) - f(11)
    frc(3, l) = frc(3, l) - f(12)

! Now get virial.

    do q = 1, 3
      do p = 1, 3
        virial(p, q) = virial(p, q) + f(p) * crd(q, i) + &
                                      f(p + 3) * crd(q, j) + &
                                      f(p + 6) * crd(q, k) + &
                                      f(p + 9) * crd(q, l)
      end do
    end do
  end do

  return

end subroutine am_stretch_torsions_get_ene_frc

end module amoeba_stretch_torsions_mod
