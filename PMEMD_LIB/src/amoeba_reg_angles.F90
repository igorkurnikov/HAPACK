#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_reg_angles_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_reg_angles_mod

  implicit none

  private

! Data that should be broadcast to slave processes from the master:

  integer, parameter    :: amoeba_reg_angles_int_cnt = 4

  integer, save         :: do_amoeba_reg_angles_flag, num_list, &
                                        num_params, ftable_degree

  integer, allocatable, save            :: list(:,:)
  double precision, allocatable, save   :: force_constant(:)
  double precision, allocatable, save   :: equil_value(:)
  double precision, allocatable, save   :: ftable_coeff(:)

  ! BUGBUG - move to local storage:

  double precision, save                :: energy
  double precision, save                :: virial(3,3)

  double precision, parameter           :: pt999999 = 0.999999d0
  double precision, parameter           :: pi= 3.14159265358979323846d0
  double precision, parameter           :: radians_to_degrees = 180.d0 / pi
  double precision, parameter           :: degrees_to_radians = pi / 180.d0

  public        am_reg_angles_zero_flag
  public        am_reg_angles_set_user_bit
  public        am_reg_angles_eval
  public        set_reg_angles_list
  public        set_reg_angles_params
  public        set_ftable
  public        set_valid_bit
  
contains

!*******************************************************************************!
! Subroutine:  am_reg_angles_zero_flag
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_reg_angles_zero_flag
  implicit none
  do_amoeba_reg_angles_flag = 0
  return
end subroutine am_reg_angles_zero_flag

!*******************************************************************************!
! Subroutine:  am_reg_angles_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_reg_angles_set_user_bit(do_this)

  use amoeba_flags_mod

  implicit none

! Formal arguments:

  integer, intent(in)   :: do_this

  call set_user_bit(do_this, do_amoeba_reg_angles_flag)

  return

end subroutine am_reg_angles_set_user_bit

!*******************************************************************************!
! Subroutine:  am_reg_angles_eval
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_reg_angles_eval(crd, frc, ene, vir)

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
  double precision                      :: darg_dcrd(9, num_list)

  energy = 0.d0
  virial(:,:) = 0.d0

  if (do_amoeba_reg_angles_flag .ne. proceed) return

  call am_reg_angles_get_args(crd, list, equil_value, arg, darg_dcrd)

  call am_val_ftab_eval_f_df(num_list, ftable_degree, ftable_coeff,  &
                             arg, fn, dfn_darg)

  call am_reg_angles_get_ene_frc(list, force_constant,  &
                                 fn, dfn_darg, darg_dcrd, crd, frc)
  ene = energy

  vir(:,:) = vir(:, :) + virial(:, :)

  return

end subroutine am_reg_angles_eval

!*******************************************************************************!
! Subroutine:  am_reg_angles_get_args
!
! Description:
!
! This routine calculates angle function argument and its derivatives with
! respect to atomic positions of atom i, j, k. 
!
! INPUT variables:
!    crd the atomic coord array
!    alist: 4 x nangles array giving for each angle the index of the first
!           atom, index of the second atom, index of the thirs,
!           and the param table pointer
!    equil_value the list of angle ref angle
! OUTPUT variables:
!    arg, array of angle function args
!    darg_dcrdi   derivs of arg wrt crds of atom i
!    darg_dcrdj   derivs of arg wrt crds of atom j
!    darg_dcrdk   derivs of arg wrt crds of atom k
!
!*******************************************************************************

subroutine am_reg_angles_get_args(crd, alist, equil_value, arg, darg_dcrdijk)

  implicit none

! Formal arguments:

  double precision, intent(in)  :: crd(3, *)
  integer, intent(in)           :: alist(4, *)
  double precision, intent(in)  :: equil_value(*)
  double precision, intent(out) :: arg(*)
  double precision, intent(out) :: darg_dcrdijk(9, *)

! Local variables:

  integer                       :: i, j, k, n, it

  double precision              :: xij, yij, zij, xkj, ykj, zkj
  double precision              :: cosang    
  double precision              :: dotp
  double precision              :: lenij2, lenkj2
  double precision              :: lenp
  double precision              :: ang, ang0
  double precision              :: dang_dcosang
  double precision              :: dcosang_dxij, dcosang_dyij, dcosang_dzij
  double precision              :: dcosang_dxkj, dcosang_dykj, dcosang_dzkj

! note units are degrees not radians...possibly change back later

  do n = 1, num_list

    i = alist(1, n)
    j = alist(2, n)
    k = alist(3, n)
    it = alist(4, n)
    ang0 = equil_value(it)
    xij = crd(1, i) - crd(1, j)
    yij = crd(2, i) - crd(2, j)
    zij = crd(3, i) - crd(3, j)
    xkj = crd(1, k) - crd(1, j)
    ykj = crd(2, k) - crd(2, j)
    zkj = crd(3, k) - crd(3, j)

! Cosine of angle is given by dot product of rij and rkj
! divided by the product of the lengths of rij and rkj.

    dotp = xij * xkj + yij * ykj + zij * zkj
    lenij2 = xij**2 + yij**2 + zij**2
    lenkj2 = xkj**2 + ykj**2 + zkj**2
    lenp = sqrt(lenij2 * lenkj2)
    cosang = dotp / lenp

! Avoid angle of pi and 0; really you need to use arcsin formulation
! near those points. however due to severe strain you could get bad values
! even for reference angle not near 0 or pi.

    cosang = max(-pt999999, cosang)
    cosang = min(pt999999, cosang)
    ang = radians_to_degrees * acos(cosang)
    dang_dcosang = -radians_to_degrees / sqrt(1.d0-cosang**2) 

! Again note units are degrees for now
! ang = acos(cosang)
! dang_dcosang = -1.d0 / sqrt(1.d0-cosang**2) 
! Angle function argument is angle ijk minus reference angle.
! Reference angle is aparm(2, at).

    arg(n) = ang - ang0

! Deriv of dotp wrt xij is xkj; deriv of lenp^-1 wrt xij is
! lenkj^-1 * (-lenij^-2) * (xij / lenij) = -xij / (lenp * lenij2).
! Similar for others.

    dcosang_dxij = xkj / lenp - (dotp * xij) / (lenp * lenij2)
    dcosang_dyij = ykj / lenp - (dotp * yij) / (lenp * lenij2)
    dcosang_dzij = zkj / lenp - (dotp * zij) / (lenp * lenij2)
    dcosang_dxkj = xij / lenp - (dotp * xkj) / (lenp * lenkj2)
    dcosang_dykj = yij / lenp - (dotp * ykj) / (lenp * lenkj2)
    dcosang_dzkj = zij / lenp - (dotp * zkj) / (lenp * lenkj2)

! Now use the chain rule

! First the i crds

    darg_dcrdijk(1, n) = dang_dcosang * dcosang_dxij
    darg_dcrdijk(2, n) = dang_dcosang * dcosang_dyij
    darg_dcrdijk(3, n) = dang_dcosang * dcosang_dzij

! Next the k crds

    darg_dcrdijk(7, n) = dang_dcosang * dcosang_dxkj
    darg_dcrdijk(8, n) = dang_dcosang * dcosang_dykj
    darg_dcrdijk(9, n) = dang_dcosang * dcosang_dzkj

! Finally the j crds

    darg_dcrdijk(4, n) = -dang_dcosang * (dcosang_dxij + dcosang_dxkj)
    darg_dcrdijk(5, n) = -dang_dcosang * (dcosang_dyij + dcosang_dykj)
    darg_dcrdijk(6, n) = -dang_dcosang * (dcosang_dzij + dcosang_dzkj)

  end do

  return

end subroutine am_reg_angles_get_args

!*******************************************************************************!
! Subroutine:  am_reg_angles_get_ene_frc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_reg_angles_get_ene_frc(alist, force_constant, func, dfunc,  & 
                                     darg_dcrdijk, crd, frc)

  implicit none

! Formal arguments:

  integer, intent(in)                   :: alist(4, *) ! 3 atoms plus param ptr
  double precision, intent(in)          :: force_constant(*)
  double precision, intent(in)          :: func(*)
  double precision, intent(in)          :: dfunc(*)
  double precision, intent(in)          :: darg_dcrdijk(9, *)
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: frc(3, *)

! Local variables:

  integer                               :: i, j, k, it, n, p, q
  double precision                      :: term, Frc_K, f(9)
  double precision                      :: factor

  ! correct bending constant for units in degrees instead of radians

  factor = degrees_to_radians**2

  do n = 1, num_list

    i = alist(1, n)
    j = alist(2, n)
    k = alist(3, n)
    it = alist(4, n)
    Frc_K = factor * force_constant(it)
    energy = energy + Frc_K * func(n)

! Apply chain rule to get deriv of energy with respect to crds of i, j, k.
! df holds the deriv of f with respect to its arg (i.e. a-a0)
! while darg_dcrdijk holds the derivs of arg with respect to crds of i, j, k.
! Recall force is negative of grad.

    term = Frc_K * dfunc(n)
    f(1) = term * darg_dcrdijk(1, n)
    f(2) = term * darg_dcrdijk(2, n) 
    f(3) = term * darg_dcrdijk(3, n) 
    f(4) = term * darg_dcrdijk(4, n)
    f(5) = term * darg_dcrdijk(5, n) 
    f(6) = term * darg_dcrdijk(6, n) 
    f(7) = term * darg_dcrdijk(7, n)
    f(8) = term * darg_dcrdijk(8, n) 
    f(9) = term * darg_dcrdijk(9, n) 
    frc(1, i) = frc(1, i) - f(1)
    frc(2, i) = frc(2, i) - f(2)
    frc(3, i) = frc(3, i) - f(3)
    frc(1, j) = frc(1, j) - f(4)
    frc(2, j) = frc(2, j) - f(5)
    frc(3, j) = frc(3, j) - f(6)
    frc(1, k) = frc(1, k) - f(7)
    frc(2, k) = frc(2, k) - f(8)
    frc(3, k) = frc(3, k) - f(9)

! now get virial

    do q = 1, 3
      do p = 1, 3
        virial(p, q) = virial(p, q) + f(p) * crd(q, i) + &
                                      f(p + 3) * crd(q, j) + &
                                      f(p + 6) * crd(q, k)
      end do
    end do

  end do

  return

end subroutine am_reg_angles_get_ene_frc

subroutine set_reg_angles_list(list_new,n)
	
	implicit none
	integer, intent(in)         :: list_new(4,*)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	nold = size(list,2)
	num_list = n
	
    if( nold .ne. n ) then
        if( nold .gt. 0) then 
            deallocate(list)
        endif
        allocate( list(4,num_list), stat = alloc_failed )
        if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( num_list .eq. 0) return 
	    
	list(:,1:n) = list_new(:,1:n) 
	
end subroutine set_reg_angles_list

subroutine set_reg_angles_params(fc,eq,n)
	
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
	
end subroutine set_reg_angles_params

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
        do_amoeba_reg_angles_flag = ibset(do_amoeba_reg_angles_flag, valid_bit)
    else
        do_amoeba_reg_angles_flag = ibclr(do_amoeba_reg_angles_flag, valid_bit)
    endif
        
end subroutine set_valid_bit

end module amoeba_reg_angles_mod
