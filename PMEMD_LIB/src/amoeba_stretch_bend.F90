#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_stretch_bend_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_stretch_bend_mod

  implicit none

  private

  integer,save        :: do_amoeba_stretch_bend_flag, num_list, num_params

  integer, allocatable, save            :: list(:,:)
  double precision, allocatable, save   :: force_constant(:)
  double precision, allocatable, save   :: angle_equil_value(:)
  double precision, allocatable, save   :: bond1_equil_value(:)
  double precision, allocatable, save   :: bond2_equil_value(:)

  ! BUGBUG - move to local storage:

  double precision, save                :: energy
  double precision, save                :: virial(3,3)

  double precision, parameter           :: pt999999 = 0.999999d0
  double precision, parameter           :: pi= 3.14159265358979323846d0
  double precision, parameter           :: radians_to_degrees = 180.d0 / pi
  double precision, parameter           :: degrees_to_radians = pi / 180.d0

  public        am_stretch_bend_zero_flag
  public        am_stretch_bend_set_user_bit
  public        am_stretch_bend_eval

contains

!*******************************************************************************!
! Subroutine:  am_stretch_bend_zero_flag
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_stretch_bend_zero_flag
  implicit none
  do_amoeba_stretch_bend_flag = 0
  return
end subroutine am_stretch_bend_zero_flag

!*******************************************************************************!
! Subroutine:  am_stretch_bend_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_stretch_bend_set_user_bit(do_this)

  use amoeba_flags_mod

  implicit none

! Formal arguments:

  integer, intent(in)   :: do_this

  call set_user_bit(do_this, do_amoeba_stretch_bend_flag)

  return

end subroutine am_stretch_bend_set_user_bit

!*******************************************************************************!
! Subroutine:  am_stretch_bend_eval
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_stretch_bend_eval(crd, frc, ene, vir)

  use amoeba_flags_mod

  implicit none

! Formal arguments:

  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: frc(3, *)
  double precision, intent(out)         :: ene
  double precision, intent(in out)      :: vir(3, 3)

! Local variables:

  double precision                      ::  arg1(num_list)
  double precision                      ::  arg2(num_list)
  double precision                      ::  darg1_dcrd(9, num_list)
  double precision                      ::  darg2_dcrd(9, num_list)

  energy = 0.d0
  virial(:, :) = 0.d0

  if (do_amoeba_stretch_bend_flag .ne. proceed) return

! First get the geometry terms.
! First the stretch contribution.

  call am_stretch_bend_get_stretchargs(crd, list, bond1_equil_value, &
                                       bond2_equil_value, arg1, darg1_dcrd)

! Next the bend contribution. Can use standard angle function.

  call am_stretch_bend_get_bendargs(crd, list, angle_equil_value, &
                                    arg2, darg2_dcrd)

! Finally call the energy force evaluator. Energy is given by a constant
! times arg1 x arg2 --- so doesn't fit into the usual functable methodology
! which assumes functions of one argument.

  call am_stretch_bend_get_ene_frc(list, force_constant, arg1, arg2, &
                                   darg1_dcrd, darg2_dcrd, crd, frc)
  ene = energy
  vir(:, :) = vir(:, :) + virial(:, :)

  return

end subroutine am_stretch_bend_eval

!*******************************************************************************!
! Subroutine:  am_stretch_bend_get_stretchargs
!
! Description:
!
! This routine calculates ij & kj bond function arguments and 
! their derivatives with
! respect to atomic positions of atom i, j, k. This gets used
! in a stretch bend energy routine 
! NOTE that periodic boundary conditions are not used here
!      i.e. imaging is done on a per molecule basis
!
! INPUT variables:
!    mstrbend: the size of the stretchbend list
!    crd the atomic coord array
!    sblist: 4 x mstrbend array giving for each strbend the index of the first
!           atom, index of the second atom, index of the thirs,
!           and index into the stretch bend parameter table giving the 
!           force constant and equilibrium angle (as well as ideal bond dists)
! OUTPUT variables:
!    arg, array of stretch bend function args
!    darg_dcrdi   derivs of arg wrt crds of atom i
!    darg_dcrdj   derivs of arg wrt crds of atom j
!    darg_dcrdk   derivs of arg wrt crds of atom k
!
!*******************************************************************************

subroutine am_stretch_bend_get_stretchargs(crd, sblist, bond1_equil_value, &
                                           bond2_equil_value, arg, &
                                           darg_dcrdijk)

  implicit none

! Formal arguments:

  double precision, intent(in)  :: crd(3, *)
  integer, intent(in)           :: sblist(4, *)
  double precision, intent(in)  :: bond1_equil_value(*)
  double precision, intent(in)  :: bond2_equil_value(*)
  double precision, intent(out) :: arg(*)
  double precision, intent(out) :: darg_dcrdijk(9, *)

! Local variables:

  integer                       :: i, j, k, n, it
  double precision              :: xij, yij, zij, xkj, ykj, zkj, bl

  do n = 1, num_list

    i = sblist(1, n)
    j = sblist(2, n)
    k = sblist(3, n)
    it = sblist(4, n)

! First the ij bond.

    xij = crd(1, i) - crd(1, j)
    yij = crd(2, i) - crd(2, j)
    zij = crd(3, i) - crd(3, j)
    bl = sqrt(xij**2 + yij**2 + zij**2)
    arg(n) = bl - bond1_equil_value(it)

! Differentiate bl to get darg_dcrd of i.

    darg_dcrdijk(1, n) = xij / bl
    darg_dcrdijk(2, n) = yij / bl
    darg_dcrdijk(3, n) = zij / bl

! Next add the kj bond contributions.

    xkj = crd(1, k) - crd(1, j)
    ykj = crd(2, k) - crd(2, j)
    zkj = crd(3, k) - crd(3, j)
    bl = sqrt(xkj**2 + ykj**2 + zkj**2)
    arg(n) = arg(n) + bl - bond2_equil_value(it)
    darg_dcrdijk(7, n) = xkj / bl
    darg_dcrdijk(8, n) = ykj / bl
    darg_dcrdijk(9, n) = zkj / bl

! j get the negative of the forces on i, k.

    darg_dcrdijk(4, n) = -(darg_dcrdijk(1, n) + darg_dcrdijk(7, n))
    darg_dcrdijk(5, n) = -(darg_dcrdijk(2, n) + darg_dcrdijk(8, n))
    darg_dcrdijk(6, n) = -(darg_dcrdijk(3, n) + darg_dcrdijk(9, n))

  end do

  return

end subroutine am_stretch_bend_get_stretchargs

!*******************************************************************************!
! Subroutine:  am_stretch_bend_get_bendargs
!
! Description:
!
! This routine calculates angle function argument and its derivatives with
! respect to atomic positions of atom i, j, k. 
!
! INPUT variables:
!    nangles:  number of angles in list
!    crd the atomic coord array
!    alist: 3 x nangles array giving for each angle the index of the first
!           atom, index of the second atom, index of the thirs,
! OUTPUT variables:
!    arg, array of angle function args
!    darg_dcrdijk   derivs of arg wrt crds of atom i, j, k
!
!*******************************************************************************

subroutine am_stretch_bend_get_bendargs(crd, sblist, angle_equil_value,  &
                                        arg, darg_dcrdijk)

  implicit none

! Formal arguments:

  double precision, intent(in)  :: crd(3, *)
  integer, intent(in)           :: sblist(4, *)
  double precision, intent(in)  :: angle_equil_value(*)
  double precision, intent(out) :: arg(*)
  double precision, intent(out) :: darg_dcrdijk(9, *)

! Local variables:

  integer                       :: i, j, k, n, it
  double precision              :: xij, yij, zij
  double precision              :: xkj, ykj, zkj
  double precision              :: cosang
  double precision              :: dotp
  double precision              :: lenij2, lenkj2
  double precision              :: lenp
  double precision              :: ang
  double precision              :: dang_dcosang
  double precision              :: dcosang_dxij, dcosang_dyij, dcosang_dzij
  double precision              :: dcosang_dxkj, dcosang_dykj, dcosang_dzkj
  double precision              :: pi
  double precision              :: conv

! Note units are degrees not radians...possibly change back later.

  do n = 1, num_list

    i = sblist(1, n)
    j = sblist(2, n)
    k = sblist(3, n)
    it = sblist(4, n)
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

! Avoid angle of pi and 0; really you need to use cosangle formulation
! near those points. However, due to severe strain you could get bad values
! even for reference angle not near 0 or pi.

    cosang = max(-pt999999, cosang)
    cosang = min(pt999999, cosang)
    ang = radians_to_degrees * acos(cosang)
    dang_dcosang = -radians_to_degrees / sqrt(1.d0-cosang**2) 

! Again note units are degrees for now
!   ang = acos(cosang)
!   dang_dcosang = -1.d0 / sqrt(1.d0-cosang**2) 
! Angle function argument is angle ijk minus reference angle.
! Reference angle is aparm(2, at).

    arg(n) = ang - angle_equil_value(it)

! Deriv of dotp wrt xij is xkj; deriv of lenp^-1 wrt xij is
! lenkj^-1 * (-lenij^-2) * (xij / lenij) = -xij / (lenp * lenij2).
! Similar for others.

    dcosang_dxij = xkj / lenp - (dotp * xij) / (lenp * lenij2)
    dcosang_dyij = ykj / lenp - (dotp * yij) / (lenp * lenij2)
    dcosang_dzij = zkj / lenp - (dotp * zij) / (lenp * lenij2)

    dcosang_dxkj = xij / lenp - (dotp * xkj) / (lenp * lenkj2)
    dcosang_dykj = yij / lenp - (dotp * ykj) / (lenp * lenkj2)
    dcosang_dzkj = zij / lenp - (dotp * zkj) / (lenp * lenkj2)

! Now use the chain rule.
! First the i crds.

    darg_dcrdijk(1, n) = dang_dcosang * dcosang_dxij
    darg_dcrdijk(2, n) = dang_dcosang * dcosang_dyij
    darg_dcrdijk(3, n) = dang_dcosang * dcosang_dzij

! Next the k crds.

    darg_dcrdijk(7, n) = dang_dcosang * dcosang_dxkj
    darg_dcrdijk(8, n) = dang_dcosang * dcosang_dykj
    darg_dcrdijk(9, n) = dang_dcosang * dcosang_dzkj

! Finally the j crds.

    darg_dcrdijk(4, n) = -dang_dcosang * (dcosang_dxij + dcosang_dxkj)
    darg_dcrdijk(5, n) = -dang_dcosang * (dcosang_dyij + dcosang_dykj)
    darg_dcrdijk(6, n) = -dang_dcosang * (dcosang_dzij + dcosang_dzkj)

  end do

  return

end subroutine am_stretch_bend_get_bendargs

!*******************************************************************************!
! Subroutine:  am_stretch_bend_get_ene_frc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_stretch_bend_get_ene_frc(sblist, force_constant, arg1, arg2,   &
                                       darg1_dcrdijk, darg2_dcrdijk, crd, frc)

  implicit none

! Formal arguments:

  integer, intent(in)                   :: sblist(4, *)
  double precision, intent(in)          :: force_constant(*)
  double precision, intent(in)          :: arg1(*)
  double precision, intent(in)          :: arg2(*)
  double precision, intent(in)          :: darg1_dcrdijk(9, *)
  double precision, intent(in)          :: darg2_dcrdijk(9, *)
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: frc(3, *)

! Local variables:

  integer                               :: i, j, k, n, p, q, it
  double precision                      :: Frc_K, term1, term2, f(9)

  do n = 1, num_list

    i = sblist(1, n)
    j = sblist(2, n)
    k = sblist(3, n)
    it = sblist(4, n)
    Frc_K = degrees_to_radians * force_constant(it)

    ! prod of stretch, bend terms:
    energy = energy + Frc_K * arg1(n) * arg2(n)

! Apply product rule to get deriv of energy with respect to crds of i, j, k
! while e.g. darg1_dcrdi holds the derivs of arg1 with respect to crds of i.
! Recall force is negative of grad.

    term1 = Frc_K * arg1(n)
    term2 = Frc_K * arg2(n)

    f(1) = term2 * darg1_dcrdijk(1, n) + term1 * darg2_dcrdijk(1, n)
    f(2) = term2 * darg1_dcrdijk(2, n) + term1 * darg2_dcrdijk(2, n)
    f(3) = term2 * darg1_dcrdijk(3, n) + term1 * darg2_dcrdijk(3, n)
    f(4) = term2 * darg1_dcrdijk(4, n) + term1 * darg2_dcrdijk(4, n)
    f(5) = term2 * darg1_dcrdijk(5, n) + term1 * darg2_dcrdijk(5, n)
    f(6) = term2 * darg1_dcrdijk(6, n) + term1 * darg2_dcrdijk(6, n)
    f(7) = term2 * darg1_dcrdijk(7, n) + term1 * darg2_dcrdijk(7, n)
    f(8) = term2 * darg1_dcrdijk(8, n) + term1 * darg2_dcrdijk(8, n)
    f(9) = term2 * darg1_dcrdijk(9, n) + term1 * darg2_dcrdijk(9, n)

    frc(1, i) = frc(1, i) - f(1)
    frc(2, i) = frc(2, i) - f(2)
    frc(3, i) = frc(3, i) - f(3)
    frc(1, j) = frc(1, j) - f(4)
    frc(2, j) = frc(2, j) - f(5)
    frc(3, j) = frc(3, j) - f(6)
    frc(1, k) = frc(1, k) - f(7)
    frc(2, k) = frc(2, k) - f(8)
    frc(3, k) = frc(3, k) - f(9)

! Now get virial.

    do q = 1, 3
      do p = 1, 3
        virial(p, q) = virial(p, q) + f(p) * crd(q, i) + &
                                      f(p + 3) * crd(q, j) + &
                                      f(p + 6) * crd(q, k)
      end do
    end do
  end do

  return

end subroutine am_stretch_bend_get_ene_frc

subroutine set_stretch_bend_list(list_new,n)
	
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
	
end subroutine set_stretch_bend_list

subroutine set_stretch_bend_params(fc,aeq,beq1,beq2,n)
	
	implicit none
	double precision, intent(in)         :: fc(n)
	double precision, intent(in)         :: aeq(n)
	double precision, intent(in)         :: beq1(n)
	double precision, intent(in)         :: beq2(n)
	integer, intent(in)         :: n
	
	integer          ::   i
	integer          ::   nold
	integer          ::   alloc_failed
	
	nold = size(force_constant)
	num_params = n
	
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(force_constant)
	      deallocate(angle_equil_value)
	      deallocate(bond1_equil_value)
	      deallocate(bond2_equil_value)
	   endif
	   allocate( force_constant(num_params), angle_equil_value(num_params), bond1_equil_value(num_params), &
                     bond2_equil_value(num_params), stat = alloc_failed )
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( num_params .eq. 0) return
	    
	force_constant(1:num_params)     = fc(1:num_params) 
	angle_equil_value(1:num_params)  = aeq(1:num_params)
	bond1_equil_value(1:num_params)  = beq1(1:num_params)
	bond2_equil_value(1:num_params)  = beq2(1:num_params)
	
end subroutine set_stretch_bend_params

subroutine set_valid_bit(ival)

    use amoeba_flags_mod
    
    implicit none
	integer, intent(in)         :: ival
    
    if(ival .eq. 1) then
        do_amoeba_stretch_bend_flag = ibset(do_amoeba_stretch_bend_flag, valid_bit)
    else
        do_amoeba_stretch_bend_flag = ibclr(do_amoeba_stretch_bend_flag, valid_bit)
    endif
        
end subroutine set_valid_bit

end module amoeba_stretch_bend_mod
