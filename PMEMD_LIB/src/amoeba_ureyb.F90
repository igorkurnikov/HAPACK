#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_ureyb_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_ureyb_mod

  implicit none

  private

  ! Data that should be broadcast to slave processes from the master:

  integer, save         :: do_amoeba_ureyb_flag, num_list, num_params, &
                                ftable_degree

  integer, allocatable, save            :: list(:,:)
  double precision, allocatable, save   :: force_constant(:)
  double precision, allocatable, save   :: equil_value(:)
  double precision, allocatable, save   :: ftable_coeff(:)

  ! BUGBUG - move to local storage:

  double precision, save                :: energy
  double precision, save                :: virial(3,3)

  public        am_ureyb_zero_flag
  public        am_ureyb_set_user_bit
  public        am_ureyb_eval

contains

!*******************************************************************************!
! Subroutine:  am_ureyb_zero_flag
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_ureyb_zero_flag
  implicit none
  do_amoeba_ureyb_flag = 0
  return
end subroutine am_ureyb_zero_flag


!*******************************************************************************!
! Subroutine:  am_ureyb_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_ureyb_set_user_bit(do_this)

  use amoeba_flags_mod

  implicit none

! Formal arguments:

  integer, intent(in)   :: do_this

  call set_user_bit(do_this, do_amoeba_ureyb_flag)

  return

end subroutine am_ureyb_set_user_bit


!*******************************************************************************!
! Subroutine:  am_ureyb_eval
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_ureyb_eval(crd, frc, ene, vir)

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
  double precision                      :: darg_dcrd(3, num_list)

! Initialize.

  energy = 0.d0
  virial(:, :) = 0.d0

  if (do_amoeba_ureyb_flag .ne. proceed) return

  call am_ureyb_get_args(crd, list, equil_value, arg, darg_dcrd)

  call am_val_ftab_eval_f_df(num_list, ftable_degree, ftable_coeff,  &
                             arg, fn, dfn_darg)

  call am_ureyb_get_ene_frc(list, force_constant, fn, dfn_darg, darg_dcrd, &
                            crd, frc)
  ene = energy

  vir(:,:) = vir(:, :) + virial(:, :)

  return

end subroutine am_ureyb_eval

!*******************************************************************************!
! Subroutine:  am_ureyb_get_args
!
! Description:
!
! These urey bradley routines are essentially identical to the regular bond 
! routines---
!
!*******************************************************************************

subroutine am_ureyb_get_args(crd, blist, equil_value, arg, darg_dcrdi)

  implicit none

  ! Formal arguments:

  double precision, intent(in)          :: crd(3, *)
  integer, intent(in)                   :: blist(3, *)
  double precision, intent(in)          :: equil_value(*)
  double precision, intent(out)         :: arg(*)
  double precision, intent(out)         :: darg_dcrdi(3, *)

! Local variables:

  integer                               :: i, j, bpar, n, it
  double precision                      :: bl, dx, dy, dz

  do n = 1, num_list

    i = blist(1, n)
    j = blist(2, n)
    it = blist(3, n)
    dx = crd(1, i) - crd(1, j)
    dy = crd(2, i) - crd(2, j)
    dz = crd(3, i) - crd(3, j)
    bl = sqrt(dx * dx + dy * dy + dz * dz)

! Recall the reference bond length is 2nd parameter.

    arg(n) = bl - equil_value(it)

! Differentiate bl to get darg_dcrdi.

    darg_dcrdi(1, n) = dx / bl
    darg_dcrdi(2, n) = dy / bl
    darg_dcrdi(3, n) = dz / bl

  end do

  return

end subroutine am_ureyb_get_args

!*******************************************************************************!
! Subroutine:  am_ureyb_get_ene_frc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_ureyb_get_ene_frc(blist, force_constant, fn, dfn, darg_dcrdi, &
                                crd, frc)

  implicit none

! Formal arguments:

  integer, intent(in)                   :: blist(3, *)
  double precision, intent(in)          :: force_constant(*)
  double precision, intent(in)          :: fn(*)
  double precision, intent(in)          :: dfn(*)
  double precision, intent(in)          :: darg_dcrdi(3, *)
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: frc(3, *)

! Local variables:

  integer                               :: i, j, bt, n, p, q, it
  double precision                      :: f(6), frc_K

  do n = 1, num_list

    i = blist(1, n)
    j = blist(2, n)
    it = blist(3, n)
    frc_K = force_constant(it)
    energy = energy + frc_K * fn(n)

! Apply chain rule to get deriv of energy with respect to crds of i.
! df holds the deriv of f with respect to its arg (i.e. b-b0)
! while darg_dcrdi holds the derivs of b-b0 with respect to crds of i

    f(1) = frc_K * dfn(n) * darg_dcrdi(1, n)
    f(2) = frc_K * dfn(n) * darg_dcrdi(2, n)
    f(3) = frc_K * dfn(n) * darg_dcrdi(3, n)

! Deriv wrt j is opposite to that wrt i.

    f(4) = -f(1) 
    f(5) = -f(2) 
    f(6) = -f(3)

! Recall force is negative of grad.

    frc(1, i) = frc(1, i) - f(1)
    frc(2, i) = frc(2, i) - f(2)
    frc(3, i) = frc(3, i) - f(3)
    frc(1, j) = frc(1, j) - f(4)
    frc(2, j) = frc(2, j) - f(5)
    frc(3, j) = frc(3, j) - f(6)

! Update the virial.

    do q = 1, 3
      do p = 1, 3
        virial(p, q) = virial(p, q) + f(p) * crd(q, i) + f(p + 3) * crd(q, j)
      end do
    end do
  end do

  return

end subroutine am_ureyb_get_ene_frc

subroutine set_ureyb_list(list_new,n)
	
	implicit none
	integer, intent(in)         :: list_new(3,*)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	nold = size(list,2)
	num_list = n
	
    if( nold .ne. n ) then
        if( nold .gt. 0) then 
            deallocate(list)
        endif
        allocate( list(3,num_list), stat = alloc_failed )
        if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( num_list .eq. 0) return 
	    
	list(:,1:n) = list_new(:,1:n) 
	
end subroutine set_ureyb_list

subroutine set_ureyb_params(fc,eq,n)
	
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
	
end subroutine set_ureyb_params

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
        do_amoeba_ureyb_flag = ibset(do_amoeba_ureyb_flag, valid_bit)
    else
        do_amoeba_ureyb_flag = ibclr(do_amoeba_ureyb_flag, valid_bit)
    endif
        
end subroutine set_valid_bit

end module amoeba_ureyb_mod
