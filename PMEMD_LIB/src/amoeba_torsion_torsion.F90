#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_torsion_torsion_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_torsion_torsion_mod

  implicit none

  private

  integer,save ::   do_amoeba_tor_tor_flag, num_list, num_params

  integer, allocatable, save            :: list(:,:)

  type angle_angle_functable
    integer                     :: dim1 = 0
    integer                     :: dim2 = 0
    double precision, pointer   :: angle1(:) => null()
    double precision, pointer   :: angle2(:) => null()
    double precision, pointer   :: func(:,:) => null()
    double precision, pointer   :: dfunc_dangle1(:,:) => null()
    double precision, pointer   :: dfunc_dangle2(:,:) => null()
    double precision, pointer   :: d2func_dangle1_dangle2(:,:) => null()
  end type  angle_angle_functable

  type(angle_angle_functable), allocatable, save  :: torsion_torsion_tbl(:) 

  double precision, parameter   :: pi = 3.14159265358979323846d0
  double precision, parameter   :: radians_to_degrees = 180.d0 / pi
  double precision, parameter   :: degrees_to_radians = pi / 180.d0

  ! BUGBUG - move to local storage:

  double precision, save                :: energy
  double precision, save                :: virial(3,3)

  public        am_tor_tor_zero_flag
  public        am_tor_tor_set_user_bit
  public        am_tor_tor_eval

contains

!*******************************************************************************!
! Subroutine:  am_tor_tor_zero_flag
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_tor_tor_zero_flag
  implicit none
  do_amoeba_tor_tor_flag = 0
  return
end subroutine am_tor_tor_zero_flag

!*******************************************************************************!
! Subroutine:  am_tor_tor_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_tor_tor_set_user_bit(do_this)

  use amoeba_flags_mod

  implicit none

! Formal arguments:

  integer, intent(in)   :: do_this

  call set_user_bit(do_this, do_amoeba_tor_tor_flag)

  return

end subroutine am_tor_tor_set_user_bit

!*******************************************************************************!
! Subroutine:  am_tor_tor_eval
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_tor_tor_eval(crd, frc, ene, vir)

  use amoeba_flags_mod

  implicit none

! Formal arguments:

  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: frc(3, *)
  double precision, intent(out)         :: ene
  double precision, intent(in out)      :: vir(3, 3)

! Local variables:

  double precision                      :: arg1(num_list)
  double precision                      :: arg2(num_list)
  double precision                      :: darg1_dcrd(12, num_list)
  double precision                      :: darg2_dcrd(12, num_list)
  double precision                      :: func(num_list)
  double precision                      :: dfunc_darg1(num_list)
  double precision                      :: dfunc_darg2(num_list)

  energy = 0.d0
  virial(:, :) = 0.d0

  if (do_amoeba_tor_tor_flag .ne. proceed) return

  call am_tor_tor_get_args(list, crd, arg1, arg2, darg1_dcrd, darg2_dcrd)

  call am_tor_tor_func(list, arg1, arg2, torsion_torsion_tbl, &
                       func, dfunc_darg1, dfunc_darg2)

  call am_tor_tor_get_ene_frc(list, crd, func, dfunc_darg1, dfunc_darg2, &
                              darg1_dcrd, darg2_dcrd, frc)
  ene = energy

  vir(:, :) = vir(:, :) + virial(:, :)

  return

end subroutine am_tor_tor_eval

!*******************************************************************************!
! Subroutine:  am_tor_tor_get_args
!
! Description:
!
! This routine calculates torsion-torsion function arguments and their
!  derivatives with  respect to atomic positions of atom i, j, k, l, m.
!
! INPUT variables:
!    crd the atomic coord array
!    list: 6 x ntortor array giving for each torsion-torsion 
!           the index of the first, second, third, fourth and fifth atoms
!           and the param table pointer
! OUTPUT variables:
!    for each torsion-torsion in list
!    arg1--the torsion angle of i, j, k, l
!    arg2--the torsion angle of j, k, l, m
!    darg1_dcrdijkl   derivs of arg1 wrt crds of atom i, j, k, l
!    darg2_dcrdjklm   derivs of arg2 wrt crds of atom j, k, l, m
!
!*******************************************************************************

subroutine am_tor_tor_get_args(list, crd, arg1, arg2,   &
                               darg1_dcrdijkl, darg2_dcrdjklm)

  implicit none

! Formal arguments:

  integer, intent(in)           :: list(6, *) !5 atoms plus param ptr
  double precision, intent(in)  :: crd(3, *)
  double precision, intent(out) :: arg1(*)
  double precision, intent(out) :: arg2(*)
  double precision, intent(out) :: darg1_dcrdijkl(12, *)
  double precision, intent(out) :: darg2_dcrdjklm(12, *)

! Local variables:

  integer                       :: i, j, k, l, m, n, p
  double precision              :: crd_abcd(12), gradphi_abcd(12)
  double precision              :: cosphi, sinphi
  double precision              :: phi

  do n = 1, num_list

    i = list(1, n)
    j = list(2, n)
    k = list(3, n)
    l = list(4, n)
    m = list(5, n)

    do p = 1, 3
      crd_abcd(p) = crd(p, i)
      crd_abcd(p + 3) = crd(p, j)
      crd_abcd(p + 6) = crd(p, k)
      crd_abcd(p + 9) = crd(p, l)
    end do

    call am_val_geom_torsion(crd_abcd, gradphi_abcd, cosphi, sinphi)

    phi = radians_to_degrees * acos(cosphi)

    if (sinphi .ge. 0.d0) then
      arg1(n) = phi
    else
      arg1(n) = -phi
    end if

    do p = 1, 12
      darg1_dcrdijkl(p, n) = radians_to_degrees * gradphi_abcd(p)
    end do

    do p = 1, 3
      crd_abcd(p) =   crd(p, j)
      crd_abcd(p + 3) = crd(p, k)
      crd_abcd(p + 6) = crd(p, l)
      crd_abcd(p + 9) = crd(p, m)
    end do

    call am_val_geom_torsion(crd_abcd, gradphi_abcd, cosphi, sinphi)

    phi = radians_to_degrees * acos(cosphi)

    if (sinphi .ge. 0.d0) then
      arg2(n) = phi
    else
      arg2(n) = -phi
    end if

    do p = 1, 12
      darg2_dcrdjklm(p, n) = radians_to_degrees * gradphi_abcd(p)
    end do

  end do

  return

end subroutine am_tor_tor_get_args

!*******************************************************************************!
! Subroutine:  am_tor_tor_func
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_tor_tor_func(list, arg1, arg2, tortor_table,  &
                           func, dfunc_darg1, dfunc_darg2)

  implicit none

! Formal arguments:

  integer, intent(in)                     :: list(6, *) !5 atoms plus param ptr

  double precision, intent(in)            :: arg1(*)
  double precision, intent(in)            :: arg2(*)
  type(angle_angle_functable), intent(in) :: tortor_table(*)
  double precision, intent(out)           :: func(*)
  double precision, intent(out)           :: dfunc_darg1(*)
  double precision, intent(out)           :: dfunc_darg2(*)

! Local variables:

  integer                                 :: am_val_real_array_index
  integer                                 :: n, it, ind1, ind2
  double precision                        :: ang1, ang2
  double precision                        :: ang1_lo, ang1_hi
  double precision                        :: ang2_lo, ang2_hi
  double precision                        :: f(4), e
  double precision                        :: df_da1(4), df_da2(4)
  double precision                        :: d2f_da1_da2(4)
  double precision                        :: de_dang1, de_dang2

  do n = 1, num_list

    it = list(6, n)
    ang1 = arg1(n)
    ang2 = arg2(n)

    ind1 = am_val_real_array_index(ang1, tortor_table(it)%angle1, &
                                  tortor_table(it)%dim1) 

    ind2 = am_val_real_array_index(ang2, tortor_table(it)%angle1, &
                                  tortor_table(it)%dim1) 

    ang1_lo = tortor_table(it)%angle1(ind1)
    ang1_hi = tortor_table(it)%angle1(ind1 + 1)
    ang2_lo = tortor_table(it)%angle2(ind2)
    ang2_hi = tortor_table(it)%angle2(ind2 + 1)

    ! counter-clockwise order around surrounding table vertices

    f(1) = tortor_table(it)%func(ind1, ind2)
    f(2) = tortor_table(it)%func(ind1 + 1, ind2)
    f(3) = tortor_table(it)%func(ind1 + 1, ind2 + 1)
    f(4) = tortor_table(it)%func(ind1, ind2 + 1)

    df_da1(1) = tortor_table(it)%dfunc_dangle1(ind1, ind2)
    df_da1(2) = tortor_table(it)%dfunc_dangle1(ind1 + 1, ind2)
    df_da1(3) = tortor_table(it)%dfunc_dangle1(ind1 + 1, ind2 + 1)
    df_da1(4) = tortor_table(it)%dfunc_dangle1(ind1, ind2 + 1)

    df_da2(1) = tortor_table(it)%dfunc_dangle2(ind1, ind2)
    df_da2(2) = tortor_table(it)%dfunc_dangle2(ind1 + 1, ind2)
    df_da2(3) = tortor_table(it)%dfunc_dangle2(ind1 + 1, ind2 + 1)
    df_da2(4) = tortor_table(it)%dfunc_dangle2(ind1, ind2 + 1)

    d2f_da1_da2(1) = tortor_table(it)%d2func_dangle1_dangle2(ind1, ind2)
    d2f_da1_da2(2) = tortor_table(it)%d2func_dangle1_dangle2(ind1 + 1, ind2)
    d2f_da1_da2(3) = tortor_table(it)%d2func_dangle1_dangle2(ind1 + 1, ind2 + 1)
    d2f_da1_da2(4) = tortor_table(it)%d2func_dangle1_dangle2(ind1, ind2 + 1)

    call am_val_bcuint1(f, df_da1, df_da2, d2f_da1_da2, ang1_lo, ang1_hi, &
                        ang2_lo, ang2_hi, ang1, ang2, e, de_dang1, de_dang2)
    func(n) = e
    dfunc_darg1(n) = de_dang1
    dfunc_darg2(n) = de_dang2

  end do

  return

end subroutine am_tor_tor_func

!*******************************************************************************!
! Subroutine:  am_tor_tor_get_ene_frc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_tor_tor_get_ene_frc(list, crd, fn, dfn_darg1, dfn_darg2, &
                                  darg1_dcrdijkl, darg2_dcrdjklm, frc)

  implicit none

! Formal arguments:

  integer, intent(in)                   :: list(6, *) !5 atoms plus param ptr
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in)          :: fn(*)
  double precision, intent(in)          :: dfn_darg1(*)
  double precision, intent(in)          :: dfn_darg2(*)
  double precision, intent(in)          :: darg1_dcrdijkl(12, *)
  double precision, intent(in)          :: darg2_dcrdjklm(12, *)
  double precision, intent(in out)      :: frc(3, *)

! Local variables:

  integer                               ::  i, j, k, l, m, n, p, q
  double precision                      :: f(12), g(12), term

  do n = 1, num_list

    i = list(1, n)
    j = list(2, n)
    k = list(3, n)
    l = list(4, n)
    m = list(5, n)

    energy = energy + fn(n)

! Apply chain rule to get deriv of energy with respect to crds of i, j, k, l, m.
! dfn_darg1 holds the deriv of fn with respect to arg1.
! dfn_darg2 holds the deriv of fn with respect to arg2
! while darg1_dcrdijkl holds the derivs of arg1 with respect to crds of
! i, j, k, l and darg2_dcrdjklm holds the derivs of arg2 with respect to crds of
! j, k, l, m.
! Recall force is negative of grad.

    term = dfn_darg1(n)

    do p = 1, 12
      f(p) = term * darg1_dcrdijkl(p, n)
    end do

    term = dfn_darg2(n)

    do p = 1, 12
      g(p) = term * darg2_dcrdjklm(p, n)
    end do

    do p = 1, 3
      frc(p, i) = frc(p, i) - f(p)
      frc(p, j) = frc(p, j) - f(p + 3) - g(p)
      frc(p, k) = frc(p, k) - f(p + 6) - g(p + 3)
      frc(p, l) = frc(p, l) - f(p + 9) - g(p + 6)
      frc(p, m) = frc(p, m) - g(p + 9)
    end do

! Now get virial.

    do q = 1, 3
      do p = 1, 3
        virial(p, q) = virial(p, q) + f(p) * crd(q, i) + &
                                      f(p + 3) * crd(q, j) + &
                                      f(p + 6) * crd(q, k) + &
                                      f(p + 9) * crd(q, l) + &
                                      g(p) * crd(q, j) + &
                                      g(p + 3) * crd(q, k) + &
                                      g(p + 6) * crd(q, l) + &
                                      g(p + 9) * crd(q, m)
      end do
    end do

  end do

  return

end subroutine am_tor_tor_get_ene_frc

subroutine set_torsion_torsion_list(list_new,n)
	
	implicit none
	integer, intent(in)         :: list_new(6,*)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	nold = size(list,2)
	num_list = n
	
    if( nold .ne. n ) then
        if( nold .gt. 0) then 
            deallocate(list)
        endif
        allocate( list(6,num_list), stat = alloc_failed )
        if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( num_list .eq. 0) return 
	    
	list(:,1:n) = list_new(:,1:n) 
	
end subroutine set_torsion_torsion_list

subroutine set_torsion_torsion_num_params(n)
	
	implicit none
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	integer          ::   i
	integer          ::   dim1 
	integer          ::   dim2 
	
	nold = size(torsion_torsion_tbl) 
	num_params = n
	
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      do i = 1, nold
	        if( associated(torsion_torsion_tbl(i)%angle1) ) deallocate(torsion_torsion_tbl(i)%angle1)
	        if( associated(torsion_torsion_tbl(i)%angle2) ) deallocate(torsion_torsion_tbl(i)%angle2)
	        if( associated(torsion_torsion_tbl(i)%func)   ) deallocate(torsion_torsion_tbl(i)%func)
	        if( associated(torsion_torsion_tbl(i)%dfunc_dangle1) ) deallocate(torsion_torsion_tbl(i)%dfunc_dangle1)
	        if( associated(torsion_torsion_tbl(i)%dfunc_dangle2) ) deallocate(torsion_torsion_tbl(i)%dfunc_dangle2)
	        if( associated(torsion_torsion_tbl(i)%d2func_dangle1_dangle2) ) deallocate(torsion_torsion_tbl(i)%d2func_dangle1_dangle2) 
	      enddo
	      deallocate(torsion_torsion_tbl)
	   endif
	   allocate( torsion_torsion_tbl(num_params), stat = alloc_failed )
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	return
	   
end subroutine set_torsion_torsion_num_params

subroutine set_torsion_torsion_params(i,dim1,dim2,angle1,angle2,func,dfunc_dangle1,dfunc_dangle2,d2func_dangle1_dangle2)
	
	implicit none
	integer, intent(in)          :: dim1
	integer, intent(in)          :: dim2
	integer, intent(in)          :: i
	double precision, intent(in) :: angle1(dim1)
	double precision, intent(in) :: angle2(dim2)
	double precision, intent(in) :: func(dim1,dim2)
	double precision, intent(in) :: dfunc_dangle1(dim1,dim2)
	double precision, intent(in) :: dfunc_dangle2(dim1,dim2)
	double precision, intent(in) :: d2func_dangle1_dangle2(dim1,dim2)
	
	integer  :: dim1_old
	integer  :: dim2_old
	integer  ::   alloc_failed
	
	dim1_old = size(torsion_torsion_tbl(i)%angle1);
	dim2_old = size(torsion_torsion_tbl(i)%angle2);
	
	if( (dim1_old .gt. 0 .and. dim1_old .ne. dim1) .or. (dim2_old .gt. 0 .and. dim2_old .ne. dim2)) then
	    deallocate(torsion_torsion_tbl(i)%angle1)
	    deallocate(torsion_torsion_tbl(i)%angle2)
	    deallocate(torsion_torsion_tbl(i)%func)
	    deallocate(torsion_torsion_tbl(i)%dfunc_dangle1)
	    deallocate(torsion_torsion_tbl(i)%dfunc_dangle2)
	    deallocate(torsion_torsion_tbl(i)%d2func_dangle1_dangle2)
	endif
	
	if( dim1 .ne. dim1_old .or. dim2 .ne. dim2_old ) then 
        torsion_torsion_tbl(i)%dim1 = dim1
        torsion_torsion_tbl(i)%dim2 = dim2
        allocate(torsion_torsion_tbl(i)%angle1(dim1), &
                 torsion_torsion_tbl(i)%angle2(dim2), &
                 torsion_torsion_tbl(i)%func(dim1, dim2), &
                 torsion_torsion_tbl(i)%dfunc_dangle1(dim1, dim2), &
                 torsion_torsion_tbl(i)%dfunc_dangle2(dim1, dim2), &
                 torsion_torsion_tbl(i)%d2func_dangle1_dangle2(dim1, dim2), &
                 stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
    endif
	    
	torsion_torsion_tbl(i)%angle1(1:dim1) = angle1(1:dim1)
	torsion_torsion_tbl(i)%angle2(1:dim2) = angle1(1:dim2)
	torsion_torsion_tbl(i)%func(:,:) = func(:,:)
	torsion_torsion_tbl(i)%dfunc_dangle1(:,:) = dfunc_dangle1(:,:)
	torsion_torsion_tbl(i)%dfunc_dangle2(:,:) = dfunc_dangle2(:,:)
	torsion_torsion_tbl(i)%d2func_dangle1_dangle2(:,:) = d2func_dangle1_dangle2(:,:)
	
end subroutine set_torsion_torsion_params

subroutine set_valid_bit(ival)

    use amoeba_flags_mod
    
    implicit none
	integer, intent(in)         :: ival
    
    if(ival .eq. 1) then
        do_amoeba_tor_tor_flag = ibset(do_amoeba_tor_tor_flag, valid_bit)
    else
        do_amoeba_tor_tor_flag = ibclr(do_amoeba_tor_tor_flag, valid_bit)
    endif
        
end subroutine set_valid_bit

end module amoeba_torsion_torsion_mod
