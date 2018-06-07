#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_opbend_angles_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_opbend_angles_mod

  implicit none

  private

! Data that should be broadcast to slave processes from the master:

  integer, save         :: do_amoeba_opbend_angles_flag, &
                                        num_list, num_params, ftable_degree

  integer, allocatable, save            :: list(:,:)
  double precision, allocatable, save   :: force_constant(:)
  double precision, allocatable, save   :: ftable_coeff(:)

  ! BUGBUG - move to local storage:

  double precision, save                :: energy
  double precision, save                :: virial(3,3)

  double precision, parameter   :: pt999999 = 0.999999d0
  double precision, parameter   :: pi = 3.14159265358979323846d0
  double precision, parameter   :: radians_to_degrees = 180.d0 / pi
  double precision, parameter   :: degrees_to_radians = pi / 180.d0
  double precision, parameter   :: opbend_unit = 0.02191418d0

  public        am_opbend_angles_zero_flag
  public        am_opbend_angles_set_user_bit
  public        am_opbend_angles_eval

contains

!*******************************************************************************!
! Subroutine:  am_opbend_angles_zero_flag
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_opbend_angles_zero_flag
  implicit none
  do_amoeba_opbend_angles_flag = 0
  return
end subroutine am_opbend_angles_zero_flag

!*******************************************************************************
!
! Subroutine:  alloc_amoeba_opbend_angles_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_amoeba_opbend_angles_mem(num_ints, num_reals)

  use parallel_dat_mod
  use amoeba_flags_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals     

! Local variables:

  integer                       :: alloc_failed

  ! List has 4 atoms plus param ptr...

  allocate(list(5, num_list), &
           force_constant(num_params), &
           ftable_coeff(0:ftable_degree), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(force_constant) + size(ftable_coeff)
  num_ints = num_ints + size(list)

  return

end subroutine alloc_amoeba_opbend_angles_mem

!*******************************************************************************!
! Subroutine:  am_opbend_angles_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_opbend_angles_set_user_bit(do_this)

  use amoeba_flags_mod

  implicit none

! Formal arguments:

  integer, intent(in)   :: do_this

  call set_user_bit(do_this, do_amoeba_opbend_angles_flag)

  return

end subroutine am_opbend_angles_set_user_bit

!*******************************************************************************!
! Subroutine:  am_opbend_angles_eval
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_opbend_angles_eval(crd, frc, ene, vir)

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
  double precision                      :: arg(num_list)
  double precision                      :: darg_dcrd(12,num_list)

  energy = 0.d0
  virial(:,:) = 0.d0

  if (do_amoeba_opbend_angles_flag .ne. proceed) return

  call am_opbend_angles_get_args(crd, list, arg, darg_dcrd)

  call am_val_ftab_eval_f_df(num_list, ftable_degree, ftable_coeff,  &
                             arg, fn, dfn_darg)

  call am_opbend_angles_get_ene_frc(list, force_constant,  &
                                    fn, dfn_darg, darg_dcrd, crd, frc)
  ene = energy

  vir(:,:) = vir(:, :) + virial(:, :)

  return

end subroutine am_opbend_angles_eval

!*******************************************************************************!
! Subroutine:  am_opbend_angles_get_args
!
! Description:
!
! This routine calculates angle function argument and its derivatives with
! respect to atomic positions of atom i, j, k, l. 
!
! INPUT variables:
!    nangles:  number of angles in list
!    crd the atomic coord array
!    alist: 5 x nangles array giving for each angle the index of the first
!           atom, index of the second atom, index of the third,
!           index of the 4th for calc of projection
!           and parameter pointer
! OUTPUT variables:
!    arg, array of angle function args
!    darg_dcrdijkl   derivs of arg wrt crds of atom i, j, k, l
!
!*******************************************************************************

subroutine am_opbend_angles_get_args(crd, alist, arg, darg_dcrd_ijkl)

  implicit none

! Formal arguments:

  double precision, intent(in)  :: crd(3, *)
  integer, intent(in)           :: alist(5, *)
  double precision, intent(out) :: arg(*)
  double precision, intent(out) :: darg_dcrd_ijkl(12, *)

! Local variables

  integer                       :: i, j, k, l, n, m, ind
  double precision              :: p(3), dp_dv1_p(3), dp_dv2_p(3)
  double precision              :: v1(3), v2(3), v3(3)
  double precision              :: dotp, lenjl2, lenjl
  double precision              :: ang, darg_dang, dang_dsinang
  double precision              :: dsinang_dv1(3)
  double precision              :: dsinang_dv2(3)
  double precision              :: dsinang_dv3(3)
  double precision              :: dsinang_dp(3)
  double precision              :: term1, term2, sinang

! note units are degrees not radians...possibly change back later

  do n = 1, num_list

    i = alist(1, n)
    j = alist(2, n)
    k = alist(3, n)
    l = alist(4, n)

    do m = 1, 3
      v1(m) = crd(m, i) - crd(m, l)
      v2(m) = crd(m, k) - crd(m, l)
      v3(m) = crd(m, j) - crd(m, l)
    end do

    call am_val_vec3d_get_perp_to_vecs(v1, v2, p, dp_dv1_p, dp_dv2_p)

    ! use arcsin formulation since angle is near 0

    dotp = p(1) * v3(1) + p(2) * v3(2) + p(3) * v3(3)
    lenjl2 = v3(1) * v3(1) + v3(2) * v3(2) + v3(3) * v3(3)
    lenjl = sqrt(lenjl2)
    sinang = dotp / lenjl
    sinang = max(-pt999999, sinang)
    sinang = min(pt999999, sinang)
    ang = radians_to_degrees * asin(sinang)
    dang_dsinang = radians_to_degrees / sqrt(1.d0 - sinang**2)
    arg(n) = dabs(ang)  ! reference angle is zero
    darg_dang = 1.d0

    if (ang < 0) darg_dang = -1.d0

    do m = 1, 3
      dsinang_dv3(m) = p(m) / lenjl - (dotp * v3(m)) / (lenjl * lenjl2)
      dsinang_dp(m) = v3(m) / lenjl - (dotp * p(m)) / lenjl
    end do

    term1 = dp_dv1_p(1) * dsinang_dp(1) + dp_dv1_p(2) * dsinang_dp(2)+ &
            dp_dv1_p(3) * dsinang_dp(3)
    term2 = dp_dv2_p(1) * dsinang_dp(1) + dp_dv2_p(2) * dsinang_dp(2)+ &
            dp_dv2_p(3) * dsinang_dp(3)

    do m = 1, 3
      dsinang_dv1(m) = term1 * p(m)
      dsinang_dv2(m) = term2 * p(m)
    end do

    do m = 1, 3
      darg_dcrd_ijkl(m, n) = darg_dang * dang_dsinang * dsinang_dv1(m)
      darg_dcrd_ijkl(m + 3, n) = darg_dang * dang_dsinang * dsinang_dv3(m)
      darg_dcrd_ijkl(m + 6, n) = darg_dang * dang_dsinang * dsinang_dv2(m)
      darg_dcrd_ijkl(m + 9, n) = -(darg_dcrd_ijkl(m, n) + &
                                     darg_dcrd_ijkl(m + 3, n) + &
                                     darg_dcrd_ijkl(m + 6, n))
    end do
    
  end do

  return

end subroutine am_opbend_angles_get_args

!*******************************************************************************!
! Subroutine:  am_opbend_angles_get_ene_frc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_opbend_angles_get_ene_frc(alist, force_constant, func, dfunc, &
                                      darg_dcrdijkl, crd, frc)

  implicit none

! Formal arguments:

  integer, intent(in)                   :: alist(5, *) ! 4 atoms plus param ptr
  double precision, intent(in)          :: force_constant(*)
  double precision, intent(in)          :: func(*)
  double precision, intent(in)          :: dfunc(*)
  double precision, intent(in)          :: darg_dcrdijkl(12, *)
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: frc(3, *)

! Local variables:

  integer                               :: i, j, k, l, n, m, p, q, ind, it
  double precision                      :: term, Frc_K, f(12)

  do n = 1, num_list

    i = alist(1, n)
    j = alist(2, n)
    k = alist(3, n)
    l = alist(4, n)
    it = alist(5, n)
    Frc_K = opbend_unit * force_constant(it)
    energy = energy + Frc_K * func(n)

! Apply chain rule to get deriv of energy with respect to crds of i, j, k, l.
! df holds the deriv of f with respect to its arg
! while darg_dcrdijkl holds the derivs of arg with respect to crds of
! i, j, k, l.  Recall force is negative of grad.

    term = Frc_K * dfunc(n)

    do m = 1, 12
      f(m) = term * darg_dcrdijkl(m, n)
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

! Now get virial:

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

end subroutine am_opbend_angles_get_ene_frc

subroutine set_opbend_angles_list(list_new,n)
	
	implicit none
	integer, intent(in)         :: list_new(5,*)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	nold = size(list,2)
	num_list = n
	
    if( nold .ne. n ) then
        if( nold .gt. 0) then 
            deallocate(list)
        endif
        allocate( list(5,num_list), stat = alloc_failed )
        if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( num_list .eq. 0) return 
	    
	list(:,1:n) = list_new(:,1:n) 
	
end subroutine set_opbend_angles_list

subroutine set_opbend_angles_params(fc,n)
	
	implicit none
	double precision, intent(in)         :: fc(n)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	nold = size(force_constant)
	num_params = n
	
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(force_constant)
	   endif
	   allocate( force_constant(num_params), stat = alloc_failed )
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( num_params .eq. 0) return
	    
	force_constant(1:num_params) = fc(1:num_params) 
	
end subroutine set_opbend_angles_params

subroutine set_ftable(ftab,n)

	implicit none
	double precision, intent(in)         :: ftab(0:n)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	nold = size(ftable_coeff) - 1
	ftable_degree = n
	
	if( nold .ne. n ) then
	   if( nold .gt. -1) then
	      deallocate(ftable_coeff)
	   endif
	   allocate( ftable_coeff(0:ftable_degree), stat = alloc_failed )
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( ftable_degree .lt. 1) return
	        
	ftable_coeff(0:n) = ftab(0:n) 
	
end subroutine set_ftable

subroutine set_valid_bit(ival)

    use amoeba_flags_mod
    
    implicit none
	integer, intent(in)         :: ival
    
    if(ival .eq. 1) then
        do_amoeba_opbend_angles_flag = ibset(do_amoeba_opbend_angles_flag, valid_bit)
    else
        do_amoeba_opbend_angles_flag = ibclr(do_amoeba_opbend_angles_flag, valid_bit)
    endif
        
end subroutine set_valid_bit

end module amoeba_opbend_angles_mod
