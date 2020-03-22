#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_trig_angles_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_trig_angles_mod

  implicit none

  private

  integer,save ::  do_amoeba_trig_angles_flag, num_list, num_params, ftable_degree

  integer, allocatable, save            :: list(:,:)
  double precision, allocatable, save   :: force_constant(:)
  double precision, allocatable, save   :: equil_value(:)
  double precision, allocatable, save   :: ftable_coeff(:)

  ! BUGBUG - move to local storage:

  double precision, save                :: energy
  double precision, save                :: virial(3,3)

  double precision                      :: xic, yic, zic
  double precision                      :: xkc, ykc, zkc
  double precision                      :: ang
  double precision                      :: dang_dxic, dang_dyic, dang_dzic
  double precision                      :: dang_dxkc, dang_dykc, dang_dzkc

  double precision, parameter           :: pt999999 = 0.999999d0
  double precision, parameter           :: pi= 3.14159265358979323846d0
  double precision, parameter           :: radians_to_degrees = 180.d0 / pi
  double precision, parameter           :: degrees_to_radians = pi / 180.d0

  public        am_trig_angles_zero_flag
  public        am_trig_angles_set_user_bit
  public        am_trig_angles_eval
  public        set_trig_angles_list
  public        set_trig_angles_params
  public        set_ftable
  public        set_valid_bit
 
contains

!*******************************************************************************!
! Subroutine:  am_trig_angles_zero_flag
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_trig_angles_zero_flag
  implicit none
  do_amoeba_trig_angles_flag = 0
  return
end subroutine am_trig_angles_zero_flag

!*******************************************************************************!
! Subroutine:  am_trig_angles_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_trig_angles_set_user_bit(do_this)

  use amoeba_flags_mod

  implicit none

! Formal arguments:

  integer, intent(in)   :: do_this

  call set_user_bit(do_this, do_amoeba_trig_angles_flag)

  return

end subroutine am_trig_angles_set_user_bit

!*******************************************************************************!
! Subroutine:  am_trig_angles_eval
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_trig_angles_eval(crd, frc, ene, vir)

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
  double precision                      :: darg_dcrd(12, num_list)

  energy = 0.d0

  virial(:,:) = 0.d0

  if (do_amoeba_trig_angles_flag .ne. proceed) return

  call am_trig_angles_get_args(crd, list, equil_value, arg, darg_dcrd)

  call am_val_ftab_eval_f_df(num_list, ftable_degree, ftable_coeff,  &
                             arg, fn, dfn_darg)

  call am_trig_angles_get_ene_frc(list, force_constant,  &
                                  fn, dfn_darg, darg_dcrd, crd, frc)
  ene = energy

  vir(:,:) = vir(:, :) + virial(:, :)

  return

end subroutine am_trig_angles_eval

!*******************************************************************************!
! Subroutine:  am_trig_angles_get_args
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
!           and then the param pointer
!    equil_value the list of ref angles
! OUTPUT variables:
!    arg, array of angle function args
!    darg_dcrdijkl   derivs of arg wrt crds of atom i, j, k, l
!
!*******************************************************************************

subroutine am_trig_angles_get_args(crd, alist, equil_value, arg, darg_dcrdijkl)

  implicit none

! Formal arguments:

  double precision, intent(in)  :: crd(3, *)
  integer, intent(in)           :: alist(5, *)
  double precision, intent(in)  :: equil_value(*)
  double precision, intent(out) :: arg(*)
  double precision, intent(out) :: darg_dcrdijkl(12, *)

! Local variables:

  integer                       :: i, j, k, l, n, m, it

  double precision              :: xil, yil, zil
  double precision              :: xjl, yjl, zjl
  double precision              :: xkl, ykl, zkl
  double precision              :: px, py, pz
  double precision              :: siz
  double precision              :: dotp
  double precision              :: xcen, ycen, zcen
  double precision              :: ang0
  double precision              :: dang_dxcen, dang_dycen, dang_dzcen

  double precision              :: v1(3), v2(3), v3(3)
  double precision              :: p(3), dp_dv1_p(3), dp_dv2_p(3)
  double precision              :: dotp1, dotp2
  double precision              :: dxcen_dv1_p, dycen_dv1_p, dzcen_dv1_p
  double precision              :: dxcen_dv2_p, dycen_dv2_p, dzcen_dv2_p
  double precision              :: dang_dv1_p, dang_dv2_p

! note units are degrees not radians...possibly change back later.

  do n = 1, num_list

    i = alist(1, n)
    j = alist(2, n)
    k = alist(3, n)
    l = alist(4, n)
    it = alist(5, n)
    ang0 = equil_value(it)

!------first get the projected center xcen, ycen, zcen
    do m = 1, 3

      v1(m) = crd(m, i) - crd(m, l)
      v2(m) = crd(m, k) - crd(m, l)
      v3(m) = crd(m, j) - crd(m, l)
    end do

    call am_val_vec3d_get_perp_to_vecs(v1, v2, p, dp_dv1_p, dp_dv2_p)

    dotp = v3(1) * p(1) + v3(2) * p(2) + v3(3) * p(3)
    xcen = crd(1, j) - dotp * p(1)
    ycen = crd(2, j) - dotp * p(2)
    zcen = crd(3, j) - dotp * p(3)

! dcen_dv1_p = -(v3.dp_dv1_p)p - (v3.p)dp_dv1_p
! dcen_dv2_p = -(v3.dp_dv2_p)p - (v3.p)dp_dv2_p

    dotp1 = v3(1) * dp_dv1_p(1) + v3(2) * dp_dv1_p(2) + v3(3) * dp_dv1_p(3)
    dotp2 = v3(1) * dp_dv2_p(1) + v3(2) * dp_dv2_p(2) + v3(3) * dp_dv2_p(3)

    dxcen_dv1_p = -(dotp1 * p(1) + dotp * dp_dv1_p(1))
    dycen_dv1_p = -(dotp1 * p(2) + dotp * dp_dv1_p(2))
    dzcen_dv1_p = -(dotp1 * p(3) + dotp * dp_dv1_p(3))

    dxcen_dv2_p = -(dotp2 * p(1) + dotp * dp_dv2_p(1))
    dycen_dv2_p = -(dotp2 * p(2) + dotp * dp_dv2_p(2))
    dzcen_dv2_p = -(dotp2 * p(3) + dotp * dp_dv2_p(3))

    xic = crd(1, i) - xcen
    yic = crd(2, i) - ycen
    zic = crd(3, i) - zcen

    xkc = crd(1, k) - xcen
    ykc = crd(2, k) - ycen
    zkc = crd(3, k) - zcen

!--------next get angle as in a normal angle function

    call  am_trig_angles_getang_cenproj()

!  derivative of angle wrt position of projected center (base of ic & kc)

    dang_dxcen = -(dang_dxic + dang_dxkc)
    dang_dycen = -(dang_dyic + dang_dykc)
    dang_dzcen = -(dang_dzic + dang_dzkc)

! extra deriv component due to effect of v1 on ang through cen

    dang_dv1_p = dang_dxcen * dxcen_dv1_p + dang_dycen * dycen_dv1_p + &
                 dang_dzcen * dzcen_dv1_p
    dang_dv2_p = dang_dxcen * dxcen_dv2_p + dang_dycen * dycen_dv2_p + &
                 dang_dzcen * dzcen_dv2_p

    arg(n) = ang - ang0

! coordiindate v1_p given by dot product v1 and p ==> dv1_p_dv1 = p

    darg_dcrdijkl(1, n) = dang_dxic + dang_dv1_p * p(1)
    darg_dcrdijkl(2, n) = dang_dyic + dang_dv1_p * p(2)
    darg_dcrdijkl(3, n) = dang_dzic + dang_dv1_p * p(3)

! indext for atom j

    darg_dcrdijkl(4, n) = dang_dxcen
    darg_dcrdijkl(5, n) = dang_dycen
    darg_dcrdijkl(6, n) = dang_dzcen

! indext for atom k

    darg_dcrdijkl(7, n) = dang_dxkc + dang_dv2_p * p(1)
    darg_dcrdijkl(8, n) = dang_dykc + dang_dv2_p * p(2)
    darg_dcrdijkl(9, n) = dang_dzkc + dang_dv2_p * p(3)

! deriv wrt crds of l due entirely to effect of v1, v2 on ang through cen

    darg_dcrdijkl(10, n) = -(dang_dv1_p + dang_dv2_p) * p(1)
    darg_dcrdijkl(11, n) = -(dang_dv1_p + dang_dv2_p) * p(2)
    darg_dcrdijkl(12, n) = -(dang_dv1_p + dang_dv2_p) * p(3)

  end do

  return

end subroutine am_trig_angles_get_args

!*******************************************************************************!
! Subroutine:  am_trig_angles_getang_cenproj
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_trig_angles_getang_cenproj()

  implicit none

  ! Local variables:

  double precision :: dotp, lenic2, lenkc2, lenp, cosang, dang_dcosang
  double precision :: dcosang_dxic, dcosang_dyic, dcosang_dzic
  double precision :: dcosang_dxkc, dcosang_dykc, dcosang_dzkc

! Cosine of angle is given by dot product of rij and rkj
! divided by the product of the lengths of rij and rkj.

  dotp = xic * xkc + yic * ykc + zic * zkc
  lenic2 = xic**2 + yic**2 + zic**2
  lenkc2 = xkc**2 + ykc**2 + zkc**2
  lenp = sqrt(lenic2 * lenkc2)
  cosang = dotp / lenp

! Avoid angle of pi and 0; really you need to use cosangle formulation
! near those points. however due to severe strain you could get bad values
! even for reference angle not near 0 or pi.

  cosang = max(-pt999999, cosang)
  cosang = min(pt999999, cosang)
  ang = radians_to_degrees * acos(cosang)
  dang_dcosang = -radians_to_degrees / sqrt(1.d0-cosang**2) 

! deriv of dotp wrt xic is xkc; deriv of lenp^-1 wrt xic is
! lenkc^-1 * (-lenic^-2) * (xic / lenic) = -xic / (lenp * lenic2).
! Similar for others.

  dcosang_dxic = xkc / lenp - (dotp * xic) / (lenp * lenic2)
  dcosang_dyic = ykc / lenp - (dotp * yic) / (lenp * lenic2)
  dcosang_dzic = zkc / lenp - (dotp * zic) / (lenp * lenic2)
  dcosang_dxkc = xic / lenp - (dotp * xkc) / (lenp * lenkc2)
  dcosang_dykc = yic / lenp - (dotp * ykc) / (lenp * lenkc2)
  dcosang_dzkc = zic / lenp - (dotp * zkc) / (lenp * lenkc2)

! Now use the chain rule.

  dang_dxic = dang_dcosang * dcosang_dxic
  dang_dyic = dang_dcosang * dcosang_dyic
  dang_dzic = dang_dcosang * dcosang_dzic
  dang_dxkc = dang_dcosang * dcosang_dxkc
  dang_dykc = dang_dcosang * dcosang_dykc
  dang_dzkc = dang_dcosang * dcosang_dzkc

  return

end subroutine am_trig_angles_getang_cenproj

!*******************************************************************************!
! Subroutine:  am_trig_angles_get_ene_frc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_trig_angles_get_ene_frc(alist, force_constant, func, dfunc,  & 
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

  integer                               :: i, j, k, l, n, m, p, q, it
  double precision                      :: term, Frc_K, f(12)
  double precision                      :: factor

  ! correct bending constant for units in degrees instead of radians

  factor = degrees_to_radians**2

  do n = 1, num_list

    i = alist(1, n)
    j = alist(2, n)
    k = alist(3, n)
    l = alist(4, n)
    it = alist(5, n)
    Frc_K = factor * force_constant(it)
    energy = energy + Frc_K * func(n)

! Apply chain rule to get deriv of energy with respect to crds of i, j, k, l.
! df holds the deriv of f with respect to its arg while darg_dcrdijkl holds
! the derivs of arg with respect to crds of i, j, k, l.
! Recall force is negative of grad.

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

! now get virial

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

end subroutine am_trig_angles_get_ene_frc

subroutine set_trig_angles_list(list_new,n)
	
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
	
end subroutine set_trig_angles_list

subroutine set_trig_angles_params(fc,eq,n)
	
	implicit none
	double precision, intent(in)         :: fc(n)
	double precision, intent(in)         :: eq(n)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	nold = size(force_constant)
	num_params = n
	
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(force_constant)
	      deallocate(equil_value)
	   endif
	   allocate( force_constant(num_params), equil_value(num_params), stat = alloc_failed )
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( num_params .eq. 0) return
	    
	force_constant(1:num_params) = fc(1:num_params) 
	equil_value(1:num_params)    = eq(1:num_params)
	
end subroutine set_trig_angles_params

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
        do_amoeba_trig_angles_flag = ibset(do_amoeba_trig_angles_flag, valid_bit)
    else
        do_amoeba_trig_angles_flag = ibclr(do_amoeba_trig_angles_flag, valid_bit)
    endif
        
end subroutine set_valid_bit

end module amoeba_trig_angles_mod
