#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_vdw_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_vdw_mod

implicit none

private

  integer,save ::    do_amoeba_vdw_flag, vdw_param_cnt, vdw_atom_cnt

 
  ! halgren's vdw buffered by these:

  double precision, public, save                :: vdw_buf_gamma
  double precision, public, save                :: vdw_buf_delta 

  integer, public, save, allocatable            :: vdw_atom_type(:)
  integer, public, save, allocatable            :: vdw_atom_parent(:)
  double precision, public, save, allocatable   :: vdw_atom_parent_crd_wt(:)
  double precision, public, save, allocatable   :: vdw_epsilon(:,:)
  double precision, public, save, allocatable   :: vdw_radius(:,:)

  ! Data that is initialized in final amoeba setup:

  double precision, save              :: ene_vdw_longrange_factor
  double precision, save              :: vir_vdw_longrange_factor
  double precision, save              :: vdw_switch_on
  double precision, public, save      :: vdw_switch_on_2
  double precision, save              :: vdw_switch_off
  double precision, public, save      :: vdw_switch_off_2
  double precision, public, save      :: c0, c1, c2, c3, c4, c5

  ! This is allocated in initialization, but initialized in final setup:

  integer, public, save, allocatable            :: vdw_type_count(:)

  public        do_amoeba_vdw_flag
  public        am_vdw_zero_flag
  public        am_vdw_set_user_bit
  public        am_vdw_adjust_ene_frc
  public        am_vdw_longrange_ene

contains

!*******************************************************************************!
! Subroutine:  am_vdw_zero_flag
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_vdw_zero_flag
  implicit none
  do_amoeba_vdw_flag = 0
  return
end subroutine am_vdw_zero_flag

!*******************************************************************************!
! Subroutine:  am_vdw_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_vdw_set_user_bit(do_this)

  use amoeba_flags_mod

  implicit none

! Formal arguments:

  integer, intent(in) :: do_this

  if (do_this .eq. 1) then
    do_amoeba_vdw_flag = ibset(do_amoeba_vdw_flag, user_bit)
  else
    do_amoeba_vdw_flag = ibclr(do_amoeba_vdw_flag, user_bit)
  end if

  return

end subroutine am_vdw_set_user_bit

!*******************************************************************************
!
! Subroutine:  am_vdw_switch
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_vdw_switch()

  use mdin_amoeba_dat_mod, only : vdw_taper
  use mdin_ctrl_dat_mod, only : vdw_cutoff
  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Local variables:

  double precision      :: denom

  if (vdw_taper .ge. 1.d0 .or. vdw_taper .le. 0.d0) then
     write(error_msg, *) 'Bad value for vdw_taper: ', vdw_taper, 'should be fraction'
     call mol_mech_error
  end if

  vdw_switch_on = vdw_taper * vdw_cutoff
  vdw_switch_on_2 = vdw_switch_on * vdw_switch_on
  vdw_switch_off = vdw_cutoff
  vdw_switch_off_2 = vdw_cutoff * vdw_cutoff
  
  denom = (vdw_switch_off - vdw_switch_on)**5

  c0 = vdw_switch_off * vdw_switch_off_2 * &
       (vdw_switch_off_2 - 5.d0 * vdw_switch_off * vdw_switch_on + &
        10.d0 * vdw_switch_on_2) / denom

  c1 = -30.d0 * vdw_switch_on_2 * vdw_switch_off_2 / denom

  c2 = 30.d0 * (vdw_switch_off_2 * vdw_switch_on + &
                vdw_switch_off * vdw_switch_on_2) / denom

  c3 = -10.d0 * (vdw_switch_off_2 + 4.d0 * vdw_switch_off * vdw_switch_on + &
                 vdw_switch_on_2) / denom

  c4 = 15.d0 * (vdw_switch_off + vdw_switch_on) / denom

  c5 = -6.d0 / denom

  return

end subroutine am_vdw_switch

!*******************************************************************************!
! Subroutine:  am_vdw_longrange_factor
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_vdw_longrange_factor(atm_cnt)

  use mdin_amoeba_dat_mod, only : do_vdw_taper
  use gbl_constants_mod, only : PI

  implicit none

! Formal arguments:

  integer, intent(in)   :: atm_cnt

! Local variables:

  integer               :: ier, n, nt, kdel, ndel, i, j
  double precision      :: r, r1, r2, r3, r4, r5
  double precision      :: req, eps, f, sume, sumv, t1, t2, rho, switch, delr
  double precision      :: dt1drho, dt2drho, drhodr, dfdr, dswitch_dr, g1, g2

  do n = 1, vdw_param_cnt
    vdw_type_count(n) = 0
  end do

  do n = 1, atm_cnt
    nt = vdw_atom_type(n)
    vdw_type_count(nt) = vdw_type_count(nt) + 1
  end do

! Choose r2 to be 60 Angstroms or ~20 sigma.
! Choose delr to be 0.05 angstroms, or ndel = 20 * (r2 - r1).

  if (do_vdw_taper .eq. 1) then
     r1 = vdw_switch_on
  else
     r1 = vdw_switch_off
  end if

  r2 = 60.d0
  ndel = nint(100.d0 * (r2 - r1))
  delr = (r2 - r1) / ndel
  ene_vdw_longrange_factor = 0.d0
  vir_vdw_longrange_factor = 0.d0

! Use trapezoidal rule for each integral.

  do j = 1, vdw_param_cnt
    do i = 1, vdw_param_cnt
      req = vdw_radius(i, j)
      eps = vdw_epsilon(i, j)
      sume = 0.d0
      sumv = 0.d0

      do kdel = 1, ndel

        r = r1 + (kdel - 1) * delr
        rho = r / req
        t1 = ((1.d0 + vdw_buf_delta) / (rho + vdw_buf_delta))**7
        t2 = (1.d0 + vdw_buf_gamma) / (rho**7 + vdw_buf_gamma)
        dt1drho = -7.d0 * t1 / (rho + vdw_buf_delta)
        dt2drho = -7.d0 * t2 * (rho**6 / (rho**7 + vdw_buf_gamma))
        drhodr = 1.d0 / req
        f = eps * t1 * (t2 - 2.d0)
        dfdr = eps * (dt1drho * (t2 - 2.d0) + t1 * dt2drho) * drhodr

        if (r < vdw_switch_off) then
          g1 = f
          g2 = dfdr
          r2 = r * r
          r3 = r2 * r
          r4 = r3 * r
          r5 = r4 * r
          switch = c5 * r5 + c4 * r4 + c3 * r3 + c2 * r2 + c1 * r + c0
          dswitch_dr = 5.d0 * c5 * r4 + 4.d0 * c4 * r3 + 3.d0 * c3 * r2 + &
                       2.d0 * c2 * r + c1
          dfdr = switch * dfdr + f * dswitch_dr
          f = switch * f
          f = g1 - f            ! need what's not done explicitly.
          dfdr = g2 - dfdr      ! need what's not done explicitly.
        end if

        if (kdel .eq. 1 .or. kdel .eq. ndel) then
          sume = sume + 0.5d0 * r**2 * f * delr
          sumv = sumv + 0.5d0 * r**3 * dfdr * delr
        else
          sume = sume + r**2 * f * delr
          sumv = sumv + r**3 * dfdr * delr
        end if

      end do

      ! Note the 2 * pi below not 4 * pi---since we do each i, j pair 2x.

      ene_vdw_longrange_factor = ene_vdw_longrange_factor + 2.d0 * PI * &
                                 vdw_type_count(i) * vdw_type_count(j) * sume

      vir_vdw_longrange_factor = vir_vdw_longrange_factor + 2.d0 * PI * &
                                 vdw_type_count(i) * vdw_type_count(j) * sumv
    end do
  end do
  
  return

end subroutine am_vdw_longrange_factor

!*******************************************************************************!
! Subroutine:  am_vdw_longrange_ene
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_vdw_longrange_ene(ene_vdw, vir_tensor)

  use pbc_mod, only: uc_volume
  use mdin_amoeba_dat_mod, only : do_vdw_longrange
  use amoeba_flags_mod

  implicit none

! Formal arguments:

  double precision, intent(out)         :: ene_vdw
  double precision, intent(in out)      :: vir_tensor(3, 3)

   ene_vdw = 0.d0

   if (do_amoeba_vdw_flag .ne. proceed) return
   if (do_vdw_longrange .ne. 1) return

  ene_vdw = ene_vdw_longrange_factor / uc_volume

  vir_tensor(1, 1) = vir_tensor(1, 1) + vir_vdw_longrange_factor / &
                     (3.d0 * uc_volume)

  vir_tensor(2, 2) = vir_tensor(2, 2) + vir_vdw_longrange_factor / &
                     (3.d0 * uc_volume)

  vir_tensor(3, 3) = vir_tensor(3, 3) + vir_vdw_longrange_factor / &
                     (3.d0 * uc_volume)
  return

end subroutine am_vdw_longrange_ene

!*******************************************************************************!
! Subroutine:  am_vdw_adjust_ene_frc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_vdw_adjust_ene_frc(crd, atm_owner_map, num_adjust_list, &
                                 adjust_list, vdw_weight, ene_vdw, frc, virial)

  use mdin_amoeba_dat_mod, only : do_vdw_taper
  use amoeba_flags_mod
  use parallel_dat_mod
  use img_mod

  implicit none

! Formal arguments:

  double precision, intent(in)          :: crd(3, *)
  integer, intent(in)                   :: atm_owner_map(*)
  integer, intent(in)                   :: num_adjust_list
  integer, intent(in)                   :: adjust_list(3, *)
  double precision, intent(in)          :: vdw_weight(9)
  double precision, intent(in out)      :: ene_vdw
  double precision, intent(in out)      :: frc(3, *)
  double precision, intent(in out)      :: virial(3, 3)

! Local variables:

  integer                               :: i, j, k, ih, jh, it, jt, n
  integer                               :: img_id

  double precision                      :: wi, wj
  double precision                      :: xi, yi, zi
  double precision                      :: xj, yj, zj
  double precision                      :: delx, dely, delz
  double precision                      :: delr, delr2
  double precision                      :: eps
  double precision                      :: rad
  double precision                      :: rho
  double precision                      :: t1, t2
  double precision                      :: dt1drho, dt2drho, drhodr
  double precision                      :: term
  double precision                      :: dfx, dfy, dfz
  double precision                      :: vxx, vxy, vxz, vyy, vyz, vzz
  double precision                      :: prefac
  double precision                      :: switch, dswitch_dr
  double precision                      :: f, dfdr
  double precision                      :: delr3, delr4, delr5

  ene_vdw = 0.d0

  if (do_amoeba_vdw_flag .ne. proceed) return

  vxx = 0.d0
  vxy = 0.d0
  vxz = 0.d0
  vyy = 0.d0
  vyz = 0.d0
  vzz = 0.d0

  do n = 1, num_adjust_list
    if( numtasks .gt. 1) then
        img_id = gbl_atm_img_map(adjust_list(1, n))
        if (img_id .lt. my_img_lo) cycle
        if (img_id .gt. my_img_hi) cycle
    endif 

    i = adjust_list(1, n)
    j = adjust_list(2, n)
    k = adjust_list(3, n)
    ih = vdw_atom_parent(i)
    jh = vdw_atom_parent(j)
    it = vdw_atom_type(i)
    jt = vdw_atom_type(j)
    eps = vdw_epsilon(it, jt)
    
!    eps = 0.0d0 ! IGOR TMP DEBUG 
    
    rad = vdw_radius(it, jt)
    wi = vdw_atom_parent_crd_wt(i)
    wj = vdw_atom_parent_crd_wt(j)
    xi = wi * crd(1, ih) + (1.d0 - wi) * crd(1, i)
    yi = wi * crd(2, ih) + (1.d0 - wi) * crd(2, i)
    zi = wi * crd(3, ih) + (1.d0 - wi) * crd(3, i)
    xj = wj * crd(1, jh) + (1.d0 - wj) * crd(1, j)
    yj = wj * crd(2, jh) + (1.d0 - wj) * crd(2, j)
    zj = wj * crd(3, jh) + (1.d0 - wj) * crd(3, j)
    delx = xj - xi
    dely = yj - yi
    delz = zj - zi
    delr2 = delx * delx + dely * dely + delz * delz

    if (delr2 < vdw_switch_off_2) then
      delr = sqrt(delr2)
      rho = delr / rad
      t1 = ((1.d0 + vdw_buf_delta) / (rho + vdw_buf_delta))**7
      t2 = (1.d0 + vdw_buf_gamma) / (rho**7 + vdw_buf_gamma)
      dt1drho = -7.d0 * t1 / (rho + vdw_buf_delta)
      dt2drho = -7.d0 * t2 * (rho**6 / (rho**7 + vdw_buf_gamma))
      drhodr = 1.d0 / rad
      prefac = vdw_weight(k) * eps
      f = prefac * t1 * (t2 - 2.d0)
      dfdr = prefac * (dt1drho * (t2 - 2.d0) + t1 * dt2drho) * drhodr
      if (do_vdw_taper .eq. 1 .and. delr2 .gt. vdw_switch_on_2) then
        delr3 = delr2 * delr
        delr4 = delr3 * delr
        delr5 = delr4 * delr
        switch = c5 * delr5 + c4 * delr4 + c3 * delr3 + c2 * delr2 + &
                 c1 * delr + c0
        dswitch_dr = 5.d0 * c5 * delr4 + 4.d0 * c4 * delr3 + &
                     3.d0 * c3 * delr2 + 2.d0 * c2 * delr + c1
        f = switch * f
        dfdr = switch * dfdr + f * dswitch_dr
      end if
      ene_vdw = ene_vdw + f
      term = dfdr / delr
      dfx = term * delx
      dfy = term * dely
      dfz = term * delz
      
      ! Recall ddelx_dxi = -1.

      frc(1, i) = frc(1, i) + (1.d0 - wi) * dfx
      frc(2, i) = frc(2, i) + (1.d0 - wi) * dfy
      frc(3, i) = frc(3, i) + (1.d0 - wi) * dfz
      frc(1, ih) = frc(1, ih) + wi * dfx
      frc(2, ih) = frc(2, ih) + wi * dfy
      frc(3, ih) = frc(3, ih) + wi * dfz
      frc(1, j) = frc(1, j) - (1.d0 - wj) * dfx
      frc(2, j) = frc(2, j) - (1.d0 - wj) * dfy
      frc(3, j) = frc(3, j) - (1.d0 - wj) * dfz
      frc(1, jh) = frc(1, jh) - wj * dfx
      frc(2, jh) = frc(2, jh) - wj * dfy
      frc(3, jh) = frc(3, jh) - wj * dfz
      vxx = vxx + dfx * delx
      vxy = vxy + dfx * dely
      vxz = vxz + dfx * delz
      vyy = vyy + dfy * dely
      vyz = vyz + dfy * delz
      vzz = vzz + dfz * delz
    end if
  end do !n = 1, num_adjust_list

  virial(1, 1) = virial(1, 1) + vxx
  virial(1, 2) = virial(1, 2) + vxy
  virial(1, 3) = virial(1, 3) + vxz
  virial(2, 1) = virial(2, 1) + vxy
  virial(2, 2) = virial(2, 2) + vyy
  virial(2, 3) = virial(2, 3) + vyz
  virial(3, 1) = virial(3, 1) + vxz
  virial(3, 2) = virial(3, 2) + vyz
  virial(3, 3) = virial(3, 3) + vzz

  return

end subroutine am_vdw_adjust_ene_frc

subroutine set_vdw_types_list(vdw_atom_type_new, vdw_atom_parent_new, vdw_atom_parent_crd_wt_new, n)
	
	implicit none
	integer, intent(in)          :: vdw_atom_type_new(n)
	integer, intent(in)          :: vdw_atom_parent_new(n)
	double precision, intent(in) :: vdw_atom_parent_crd_wt_new(n)
	integer, intent(in)          :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	nold = size(vdw_atom_type)
	vdw_atom_cnt = n
	
    if( nold .ne. n ) then
        if( nold .gt. 0) then 
            deallocate(vdw_atom_type)
            deallocate(vdw_atom_parent)
            deallocate(vdw_atom_parent_crd_wt)
        endif
        allocate( vdw_atom_type(vdw_atom_cnt),   &
                  vdw_atom_parent(vdw_atom_cnt), &
                  vdw_atom_parent_crd_wt(vdw_atom_cnt), &
                  stat = alloc_failed )
        if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( vdw_atom_cnt .eq. 0) return 
	    
	vdw_atom_type(1:n)          = vdw_atom_type_new(1:n)
	vdw_atom_parent(1:n)        = vdw_atom_parent_new(1:n) 
	vdw_atom_parent_crd_wt(1:n) = vdw_atom_parent_crd_wt_new(1:n) 
	
	 call am_vdw_switch()
	
end subroutine set_vdw_types_list

subroutine set_vdw_params(vdw_epsilon_new,vdw_radius_new,vdw_buf_delta_new,vdw_buf_gamma_new,n)
	
	implicit none
	double precision, intent(in)         :: vdw_epsilon_new(n,n)
	double precision, intent(in)         :: vdw_radius_new(n,n)
	double precision, intent(in)         :: vdw_buf_delta_new
	double precision, intent(in)         :: vdw_buf_gamma_new
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	vdw_buf_delta = vdw_buf_delta_new
	vdw_buf_gamma = vdw_buf_gamma_new
	
	nold = size(vdw_type_count)
	vdw_param_cnt = n
	
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(vdw_type_count)
	      deallocate(vdw_epsilon)
	      deallocate(vdw_radius)
	   endif
	   allocate( vdw_type_count(vdw_param_cnt), &
	             vdw_epsilon(vdw_param_cnt,vdw_param_cnt), & 
	             vdw_radius(vdw_param_cnt,vdw_param_cnt),  &
	             stat = alloc_failed )
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( vdw_param_cnt .eq. 0) return
	    
	vdw_epsilon(1:vdw_param_cnt,1:vdw_param_cnt)    = vdw_epsilon_new(1:vdw_param_cnt,1:vdw_param_cnt) 
	vdw_radius(1:vdw_param_cnt,1:vdw_param_cnt)     = vdw_radius_new(1:vdw_param_cnt,1:vdw_param_cnt)
	
	call am_vdw_longrange_factor(vdw_atom_cnt) 
	
end subroutine set_vdw_params

subroutine set_valid_bit(ival)

    use amoeba_flags_mod
    
    implicit none
	integer, intent(in)         :: ival
    
    if(ival .eq. 1) then
        do_amoeba_vdw_flag = ibset(do_amoeba_vdw_flag, valid_bit)
    else
        do_amoeba_vdw_flag = ibclr(do_amoeba_vdw_flag, valid_bit)
    endif
        
end subroutine set_valid_bit

end module amoeba_vdw_mod
