#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_pitorsions_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_pitorsions_mod

  implicit none

  private

  integer, save         :: do_amoeba_pitorsions_flag, num_list, num_params

  integer, allocatable, save            :: list(:,:)
  double precision, allocatable, save   :: force_constant(:)
  double precision, allocatable, save   :: periodicity(:)
  double precision, allocatable, save   :: cosphase(:)
  double precision, allocatable, save   :: sinphase(:)

  ! BUGBUG - move to local storage:

  double precision, save                :: energy
  double precision, save                :: virial(3,3)

  public        am_pitorsions_zero_flag
  public        am_pitorsions_set_user_bit
  public        am_pitorsions_eval


contains

!*******************************************************************************!
! Subroutine:  am_pitorsions_zero_flag
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_pitorsions_zero_flag
  implicit none
  do_amoeba_pitorsions_flag = 0
  return
end subroutine am_pitorsions_zero_flag

!*******************************************************************************!
! Subroutine:  am_pitorsions_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_pitorsions_set_user_bit(do_this)

  use amoeba_flags_mod

  implicit none

! Formal arguments:

  integer, intent(in)   :: do_this

  call set_user_bit(do_this, do_amoeba_pitorsions_flag)

  return

end subroutine am_pitorsions_set_user_bit

!*******************************************************************************!
! Subroutine:  am_pitorsions_eval
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_pitorsions_eval(crd, frc, ene, vir)

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
  double precision                      :: darg_dcrd(18, num_list)

  energy = 0.d0
  virial(:,:) = 0.d0

  if (do_amoeba_pitorsions_flag .ne. proceed) return

  call am_pitorsions_get_args(crd, list, cosarg, sinarg, darg_dcrd)

  call am_pitorsions_func(list, force_constant, periodicity,  &
                          cosphase, sinphase, cosarg, sinarg, fn, dfn_darg)

  call am_pitorsions_get_ene_frc(list, fn, dfn_darg, darg_dcrd, crd, frc)

  ene = energy
  vir(:, :) = vir(:, :) + virial(:, :)

  return

end subroutine am_pitorsions_eval

!*******************************************************************************!
! Subroutine:  am_pitorsions_get_args
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_pitorsions_get_args(crd, list, cosarg, sinarg, darg_dcrd_ijklmn)

  implicit none

! Formal arguments:

  double precision, intent(in)  :: crd(3, *)
  integer, intent(in)           :: list(7, *) !6 atoms plus parm ptr
  double precision, intent(out) :: cosarg(*)
  double precision, intent(out) :: sinarg(*)
  double precision, intent(out) :: darg_dcrd_ijklmn(18, *)

! Local variables:

  integer                       :: h, i, j, k, l, m, n, p
  double precision              :: ril(3), rjl(3), rlk(3), rmk(3), rnk(3)
  double precision              :: p_ijl(3), p_mnk(3)
  double precision              :: dp_ijl_dril_p(3), dp_ijl_drjl_p(3)
  double precision              :: dp_mnk_drmk_p(3), dp_mnk_drnk_p(3)
  double precision              :: crd_abcd(12), gradphi_abcd(12)
  double precision              :: cosphi, sinphi
  double precision              :: termi, termj, termm, termn, VEC3D_dotprod

  do h = 1, num_list

    i = list(1, h)
    j = list(2, h)
    k = list(3, h)
    l = list(4, h)
    m = list(5, h)
    n = list(6, h)

    do p = 1, 3
      ril(p) = crd(p, i) - crd(p, l)
      rjl(p) = crd(p, j) - crd(p, l)
      rlk(p) = crd(p, l) - crd(p, k)
      rmk(p) = crd(p, m) - crd(p, k)
      rnk(p) = crd(p, n) - crd(p, k)
    end do

    call am_val_vec3d_get_perp_to_vecs(ril, rjl, p_ijl, &
                                       dp_ijl_dril_p, dp_ijl_drjl_p)

    call am_val_vec3d_get_perp_to_vecs(rmk, rnk, p_mnk, &
                                       dp_mnk_drmk_p, dp_mnk_drnk_p)

    ! Use the regular torsion results for artificial sites a, b, c, d.

    do p = 1, 3
      crd_abcd(p) = crd(p, k) + p_ijl(p)
      crd_abcd(p + 3) = crd(p, k)
      crd_abcd(p + 6) = crd(p, l)
      crd_abcd(p + 9) = crd(p, l) + p_mnk(p)
    end do

    call am_val_geom_torsion(crd_abcd, gradphi_abcd, cosphi, sinphi)

    cosarg(h) = cosphi
    sinarg(h) = sinphi

    ! Transfer grad wrt point a to i, j, l and grad wrt point d to m, n, k.
    ! Recall that movement of i parallel to p_ijk causes change in p_ijk,
    ! hence a in direction dp_ijl_dril_p

    termi = dp_ijl_dril_p(1) * gradphi_abcd(1) + &
            dp_ijl_dril_p(2) * gradphi_abcd(2) + &
            dp_ijl_dril_p(3) * gradphi_abcd(3)
    termj = dp_ijl_drjl_p(1) * gradphi_abcd(1) + &
            dp_ijl_drjl_p(2) * gradphi_abcd(2) + &
            dp_ijl_drjl_p(3) * gradphi_abcd(3)
    termm = dp_mnk_drmk_p(1) * gradphi_abcd(10) + &
            dp_mnk_drmk_p(2) * gradphi_abcd(11) + &
            dp_mnk_drmk_p(3) * gradphi_abcd(12)
    termn = dp_mnk_drnk_p(1) * gradphi_abcd(10) + &
            dp_mnk_drnk_p(2) * gradphi_abcd(11) + &
            dp_mnk_drnk_p(3) * gradphi_abcd(12)

    do p = 1, 3

      darg_dcrd_ijklmn(p, h) = termi * p_ijl(p)  ! for i

      darg_dcrd_ijklmn(3 + p, h) = termj * p_ijl(p) ! for j

      darg_dcrd_ijklmn(6 + p, h) = gradphi_abcd(p) + gradphi_abcd(p + 3) - &
                                     (termm + termn) * p_mnk(p) ! for k

      darg_dcrd_ijklmn(9 + p, h) = gradphi_abcd(p + 6) + &
                                     gradphi_abcd(p + 9) - &
                                     (termi + termj) * p_ijl(p) ! for l
                                     
      darg_dcrd_ijklmn(12 + p, h) = termm * p_mnk(p)  ! for m
      
      darg_dcrd_ijklmn(15 + p, h) = termn * p_mnk(p)  ! for n

    end do
    
  end do

  return

end subroutine am_pitorsions_get_args

!*******************************************************************************!
! Subroutine:  am_pitorsions_func
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_pitorsions_func(list, force_constant, periodicity, cosphase, &
                              sinphase, cosarg, sinarg, func, dfunc_darg)

  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: list(7, *)           ! 6 atoms plus parm ptr
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
  double precision              :: amplitude, period, cos_phase, sin_phase

  do n = 1, num_list

    cosine(1) = cosarg(n)
    sine(1) = sinarg(n)
    it = list(7, n)
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
        write(error_msg,*) 'AM_PITORSIONS_func: degree too big: ', degree
        call mol_mech_error
      end if
      do m = 2, degree
        cosine(m) = cosine(m-1) * cosine(1) - sine(m-1) * sine(1)
        sine(m) = sine(m-1) * cosine(1) + cosine(m-1) * sine(1)
      end do
    end if

    func(n) =  amplitude * (1.d0 + cos_phase * cosine(degree) + &
                                     sin_phase * sine(degree))

    dfunc_darg(n) = amplitude * period * (sin_phase * cosine(degree) - &
                                            cos_phase * sine(degree))
  end do 

  return

end subroutine am_pitorsions_func

!*******************************************************************************!
! Subroutine:  am_pitorsions_get_ene_frc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_pitorsions_get_ene_frc(list, func, dfunc, darg_dcrd_ijklmn, &
                                     crd, frc)

  implicit none

! Formal arguments:

  integer, intent(in)                   :: list(7, *) ! 6 atoms plus param ptr
  double precision, intent(in)          :: func(*)
  double precision, intent(in)          :: dfunc(*)
  double precision, intent(in)          :: darg_dcrd_ijklmn(18, *)
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: frc(3, *)

! Local variables:

  integer                               :: i, j, k, l, n, m, h, p, q
  double precision                      :: term, f(18)

  do h = 1, num_list

    i = list(1, h)
    j = list(2, h)
    k = list(3, h)
    l = list(4, h)
    m = list(5, h)
    n = list(6, h)

    energy = energy + func(h)

! Apply chain rule to get deriv of energy with respect to crds of i, j, k, l.
! df holds the deriv of f with respect to its arg, while darg_dcrdijkl holds
! the derivs of arg with respect to crds of i, j, k, l.
! Recall force is negative of grad

    term = dfunc(h)

    do p = 1, 18
      f(p) = term * darg_dcrd_ijklmn(p, h)
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
    frc(1, m) = frc(1, m) - f(13)
    frc(2, m) = frc(2, m) - f(14)
    frc(3, m) = frc(3, m) - f(15)
    frc(1, n) = frc(1, n) - f(16)
    frc(2, n) = frc(2, n) - f(17)
    frc(3, n) = frc(3, n) - f(18)

! Now get virial.

    do q = 1, 3
      do p = 1, 3
        virial(p, q) = virial(p, q) + f(p) * crd(q, i) + &
                                      f(p + 3) * crd(q, j) + &
                                      f(p + 6) * crd(q, k) + &
                                      f(p + 9) * crd(q, l) + &
                                      f(p + 12) * crd(q, m) + &
                                      f(p + 15) * crd(q, n)
      end do
    end do

  end do

  return

end subroutine am_pitorsions_get_ene_frc

subroutine set_pitorsions_list(list_new,n)
	
	implicit none
	integer, intent(in)         :: list_new(7,*)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	nold = size(list,2)
	num_list = n
	
    if( nold .ne. n ) then
        if( nold .gt. 0) then 
            deallocate(list)
        endif
        allocate( list(7,num_list), stat = alloc_failed )
        if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( num_list .eq. 0) return 
	    
	list(:,1:n) = list_new(:,1:n) 
	
end subroutine set_pitorsions_list

subroutine set_pitorsions_params(fc,per,phase,n)
	
	implicit none
	double precision, intent(in)         :: fc(n)
	double precision, intent(in)         :: per(n)
	double precision, intent(in)         :: phase(n)
	integer, intent(in)         :: n
	
	integer          ::   i
	integer          ::   nold
	integer          ::   alloc_failed
	
	nold = size(force_constant)
	num_params = n
	
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(force_constant)
	      deallocate(periodicity)
	      deallocate(cosphase)
	      deallocate(sinphase)
	   endif
	   allocate( force_constant(num_params), periodicity(num_params), cosphase(num_params), sinphase(num_params), stat = alloc_failed )
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( num_params .eq. 0) return
	    
	force_constant(1:num_params) = fc(1:num_params) 
	periodicity(1:num_params)    = per(1:num_params)
	do i = 1, num_params
        cosphase(i) = cos(phase(i))
        sinphase(i) = sin(phase(i))
    end do 
	
end subroutine set_pitorsions_params

subroutine set_valid_bit(ival)

    use amoeba_flags_mod
    
    implicit none
	integer, intent(in)         :: ival
    
    if(ival .eq. 1) then
        do_amoeba_pitorsions_flag = ibset(do_amoeba_pitorsions_flag, valid_bit)
    else
        do_amoeba_pitorsions_flag = ibclr(do_amoeba_pitorsions_flag, valid_bit)
    endif
        
end subroutine set_valid_bit

end module amoeba_pitorsions_mod
