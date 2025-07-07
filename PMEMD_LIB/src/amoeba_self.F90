#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_self_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_self_mod

  implicit none

  private

  integer, save ::         do_amoeba_self_flag

#include "amoeba_mpole_index.i"

  public        am_self_zero_flag
  public        am_self_set_user_bit
  public        am_self_permfield
  public        am_self_dipole_field
  public        am_self_ene_torque

contains

!*******************************************************************************!
! Subroutine:  am_self_zero_flag
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_self_zero_flag
  implicit none
  do_amoeba_self_flag = 0
  return
end subroutine am_self_zero_flag

!*******************************************************************************!
! Subroutine:  am_self_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_self_set_user_bit(do_this)

  use amoeba_flags_mod
  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in) :: do_this

  ! Set the valid bit---this part always since no parmread needed.

  do_amoeba_self_flag = ibset(do_amoeba_self_flag, valid_bit)

  if (do_this .eq. 1) then ! do in all cases
    do_amoeba_self_flag = ibset(do_amoeba_self_flag, user_induce_bit)
    do_amoeba_self_flag = ibset(do_amoeba_self_flag, user_postinduce_bit)
  else if (do_this .eq. 2) then ! do the induction, not the post-induction
    do_amoeba_self_flag = ibset(do_amoeba_self_flag, user_induce_bit)
    do_amoeba_self_flag = ibclr(do_amoeba_self_flag, user_postinduce_bit)
  else if (do_this .eq. 3) then ! do the post-induction, not the induction
    do_amoeba_self_flag = ibclr(do_amoeba_self_flag, user_induce_bit)
    do_amoeba_self_flag = ibset(do_amoeba_self_flag, user_postinduce_bit)
  else if (do_this .eq. 0) then 
    do_amoeba_self_flag = ibclr(do_amoeba_self_flag, user_induce_bit)
    do_amoeba_self_flag = ibclr(do_amoeba_self_flag, user_postinduce_bit)
  else
    error_msg = 'am_self_set_user_bit: bad value of user do_this'
    call mol_mech_error
  end if

  return

end subroutine am_self_set_user_bit

!*******************************************************************************!
! Subroutine:  am_self_permfield
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_self_permfield(img_cnt, direct_field, polar_field)

  use amoeba_multipoles_mod, only : global_multipole
  use amoeba_flags_mod
  use gbl_constants_mod, only : PI
  use mdin_ewald_dat_mod, only : ew_coeff
  use timers_mod
  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: img_cnt
  double precision, intent(in out)      :: direct_field(3, *)
  double precision, intent(in out)      :: polar_field(3, *)

! Local variables:

  double precision                      :: factor
  integer                               :: n
  integer                               :: atm_id, img_id

  if (iand(do_amoeba_self_flag, proceed_induce) .ne. proceed_induce) return

  factor = 4.d0 * ew_coeff**3 / (3.d0 * sqrt(PI))

  do img_id = my_img_lo, my_img_hi
  
    atm_id = gbl_img_atm_map(img_id)

    direct_field(1,atm_id) = direct_field(1,atm_id) - factor * &
                             global_multipole(Ind_100,atm_id)
    direct_field(2,atm_id) = direct_field(2,atm_id) - factor * &
                             global_multipole(Ind_010,atm_id)
    direct_field(3,atm_id) = direct_field(3,atm_id) - factor * &
                             global_multipole(Ind_001,atm_id)
    polar_field(1,atm_id) = polar_field(1,atm_id) - factor * &
                            global_multipole(Ind_100,atm_id)
    polar_field(2,atm_id) = polar_field(2,atm_id) - factor * &
                            global_multipole(Ind_010,atm_id)
    polar_field(3,atm_id) = polar_field(3,atm_id) - factor * &
                            global_multipole(Ind_001,atm_id)
  end do

  call update_pme_time(pme_misc_timer)

  return

end subroutine am_self_permfield

!*******************************************************************************!
! Subroutine:  am_self_dipole_field
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_self_dipole_field(img_cnt, ind_dip1, ind_dip2, &
                                dip_field1, dip_field2)

  use amoeba_flags_mod
  use gbl_constants_mod, only : PI
  use mdin_ewald_dat_mod, only : ew_coeff
  use timers_mod
  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: img_cnt
  double precision, intent(in)          :: ind_dip1(3, *)
  double precision, intent(in)          :: ind_dip2(3, *)
  double precision, intent(in out)      :: dip_field1(3, *)
  double precision, intent(in out)      :: dip_field2(3, *)

! Local variables:

  double precision                      :: factor
  integer                               :: atm_id, img_id

  if (iand(do_amoeba_self_flag, proceed_induce) .ne. proceed_induce) return

  factor = 4.d0 * ew_coeff**3 / (3.d0 * sqrt(PI))

  do img_id = my_img_lo, my_img_hi
  
    atm_id = gbl_img_atm_map(img_id)

    dip_field1(1, atm_id) = dip_field1(1, atm_id) - factor * ind_dip1(1, atm_id)
    dip_field1(2, atm_id) = dip_field1(2, atm_id) - factor * ind_dip1(2, atm_id)
    dip_field1(3, atm_id) = dip_field1(3, atm_id) - factor * ind_dip1(3, atm_id)
    dip_field2(1, atm_id) = dip_field2(1, atm_id) - factor * ind_dip2(1, atm_id)
    dip_field2(2, atm_id) = dip_field2(2, atm_id) - factor * ind_dip2(2, atm_id)
    dip_field2(3, atm_id) = dip_field2(3, atm_id) - factor * ind_dip2(3, atm_id)
  end do

  call update_pme_time(pme_misc_timer)

  return

end subroutine am_self_dipole_field

!*******************************************************************************!
! Subroutine:  am_self_ene_torque
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_self_ene_torque(img_cnt, ind_dip_d, ind_dip_p, &
                              ene_perm, ene_ind)

  use amoeba_multipoles_mod, only : global_multipole, &
                                    coulomb_const_kcal_per_mole, &
                                    torque_field

  use amoeba_flags_mod
  use gbl_constants_mod, only : PI
  use mdin_ewald_dat_mod, only : ew_coeff
  use timers_mod
  use parallel_dat_mod
  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: img_cnt
  double precision, intent(in)          :: ind_dip_d(3, *)
  double precision, intent(in)          :: ind_dip_p(3, *)
  double precision, intent(in out)      :: ene_perm
  double precision, intent(in out)      :: ene_ind

! Local variables:

  double precision                      :: delx, dely, delz
  double precision                      :: B(0:4)
  double precision                      :: fac
  double precision                      :: fact
  double precision                      :: gmi(10)
  double precision                      :: phi(10)
  double precision                      :: i_di(3)
  double precision                      :: i_mi(3)
  double precision                      :: e_pp
  double precision                      :: e_ind
  double precision                      :: e_add
  double precision                      :: Rn(1)
  double precision                      :: Rn_1(4)
  double precision                      :: Rn_2(10)
  double precision                      :: Rn_3(20)
  double precision                      :: Rn_4(35)
  integer                               :: j, n
  integer                               :: atm_id, img_id
  
  double precision                      :: fScale,fScaleDip,fScaleQdp

  ene_perm = 0.d0
  ene_ind = 0.d0

  if (iand(do_amoeba_self_flag, proceed_postinduce) .ne. proceed_postinduce) &
    return

  fact = 2.d0 * ew_coeff / sqrt(PI)
  fac = -2.d0 * ew_coeff * ew_coeff
  do j = 0, 4
    B(j) = fact / (2.d0 * j + 1)
    fact = fac * fact
  end do
  
  delx = 0.d0
  dely = 0.d0
  delz = 0.d0

  n = 4

  Rn(Ind_000) = B(n)
  Rn_1(Ind_000) = B(n-1)
  Rn_1(Ind_100) = 0.d0   ! Rn_1(Ind_100) = delx * Rn(Ind_000)
  Rn_1(Ind_010) = 0.d0   ! Rn_1(Ind_010) = dely * Rn(Ind_000)
  Rn_1(Ind_010) = 0.d0   ! Rn_1(Ind_001) = delz * Rn(Ind_000)
  Rn_2(Ind_000) = B(n-2)
  Rn_2(Ind_100) = 0.0d0  ! Rn_2(Ind_100) = delx * Rn_1(Ind_000)
  Rn_2(Ind_010) = 0.0d0  ! Rn_2(Ind_010) = dely * Rn_1(Ind_000)
  Rn_2(Ind_001) = 0.0d0  ! Rn_2(Ind_001) = delz * Rn_1(Ind_000)
  Rn_2(Ind_200) = B(n-1) ! Rn_2(Ind_200) = Rn_1(Ind_000) + delx * Rn_1(Ind_100)
  Rn_2(Ind_020) = B(n-1) ! Rn_2(Ind_020) = Rn_1(Ind_000) + dely * Rn_1(Ind_010)
  Rn_2(Ind_002) = B(n-1) ! Rn_2(Ind_002) = Rn_1(Ind_000) + delz * Rn_1(Ind_001)
  Rn_2(Ind_110) = 0.0d0  ! Rn_2(Ind_110) = delx * Rn_1(Ind_010)
  Rn_2(Ind_101) = 0.0d0  ! Rn_2(Ind_101) = delx * Rn_1(Ind_001)
  Rn_2(Ind_011) = 0.0d0  ! Rn_2(Ind_011) = dely * Rn_1(Ind_001)
  Rn_3(Ind_000) = B(n-3) 
  Rn_3(Ind_100) = 0.0d0  ! Rn_3(Ind_100) = delx * Rn_2(Ind_000)
  Rn_3(Ind_010) = 0.0d0  ! Rn_3(Ind_010) = dely * Rn_2(Ind_000)
  Rn_3(Ind_001) = 0.0d0  ! Rn_3(Ind_001) = delz * Rn_2(Ind_000)
  Rn_3(Ind_200) = B(n-2) ! Rn_3(Ind_200) = Rn_2(Ind_000) + delx * Rn_2(Ind_100)
  Rn_3(Ind_020) = B(n-2) ! Rn_3(Ind_020) = Rn_2(Ind_000) + dely * Rn_2(Ind_010)
  Rn_3(Ind_002) = B(n-2) ! Rn_3(Ind_002) = Rn_2(Ind_000) + delz * Rn_2(Ind_001)
  Rn_3(Ind_110) = 0.0d0 ! Rn_3(Ind_110) = delx * Rn_2(Ind_010)
  Rn_3(Ind_101) = 0.0d0 ! Rn_3(Ind_101) = delx * Rn_2(Ind_001)
  Rn_3(Ind_011) = 0.0d0 ! Rn_3(Ind_011) = dely * Rn_2(Ind_001)
  Rn_3(Ind_300) = 0.0d0 ! Rn_3(Ind_300) = 2.d0 * Rn_2(Ind_100) + delx * Rn_2(Ind_200)
  Rn_3(Ind_030) = 0.0d0 ! Rn_3(Ind_030) = 2.d0 * Rn_2(Ind_010) + dely * Rn_2(Ind_020)
  Rn_3(Ind_003) = 0.0d0 ! Rn_3(Ind_003) = 2.d0 * Rn_2(Ind_001) + delz * Rn_2(Ind_002)
  Rn_3(Ind_210) = 0.0d0 ! Rn_3(Ind_210) = dely * Rn_2(Ind_200)
  Rn_3(Ind_201) = 0.0d0 ! Rn_3(Ind_201) = delz * Rn_2(Ind_200)
  Rn_3(Ind_120) = 0.0d0 ! Rn_3(Ind_120) = delx * Rn_2(Ind_020)
  Rn_3(Ind_021) = 0.0d0 ! Rn_3(Ind_021) = delz * Rn_2(Ind_020)
  Rn_3(Ind_102) = 0.0d0 ! Rn_3(Ind_102) = delx * Rn_2(Ind_002)
  Rn_3(Ind_012) = 0.0d0 ! Rn_3(Ind_012) = dely * Rn_2(Ind_002)
  Rn_3(Ind_111) = 0.0d0 ! Rn_3(Ind_111) = delx * Rn_2(Ind_011)
  Rn_4(Ind_000) = B(n-4) 
  Rn_4(Ind_100) = 0.0d0 ! Rn_4(Ind_100) = delx * Rn_3(Ind_000)
  Rn_4(Ind_010) = 0.0d0 ! Rn_4(Ind_010) = dely * Rn_3(Ind_000)
  Rn_4(Ind_001) = 0.0d0 ! Rn_4(Ind_001) = delz * Rn_3(Ind_000)`
  Rn_4(Ind_200) = B(n-3)  ! Rn_4(Ind_200) = Rn_3(Ind_000) + delx * Rn_3(Ind_100)
  Rn_4(Ind_020) = B(n-3)  ! Rn_4(Ind_020) = Rn_3(Ind_000) + dely * Rn_3(Ind_010)
  Rn_4(Ind_002) = B(n-3)  ! Rn_4(Ind_002) = Rn_3(Ind_000) + delz * Rn_3(Ind_001)
  Rn_4(Ind_110) = 0.0d0 ! Rn_4(Ind_110) = delx * Rn_3(Ind_010)
  Rn_4(Ind_101) = 0.0d0 ! Rn_4(Ind_101) = delx * Rn_3(Ind_001)
  Rn_4(Ind_011) = 0.0d0 ! Rn_4(Ind_011) = dely * Rn_3(Ind_001)
  Rn_4(Ind_300) = 0.0d0 ! Rn_4(Ind_300) = 2.d0 * Rn_3(Ind_100) + delx * Rn_3(Ind_200)
  Rn_4(Ind_030) = 0.0d0 ! Rn_4(Ind_030) = 2.d0 * Rn_3(Ind_010) + dely * Rn_3(Ind_020)
  Rn_4(Ind_003) = 0.0d0 ! Rn_4(Ind_003) = 2.d0 * Rn_3(Ind_001) + delz * Rn_3(Ind_002)
  Rn_4(Ind_210) = 0.0d0 ! Rn_4(Ind_210) = dely * Rn_3(Ind_200)
  Rn_4(Ind_201) = 0.0d0 ! Rn_4(Ind_201) = delz * Rn_3(Ind_200)
  Rn_4(Ind_120) = 0.0d0 ! Rn_4(Ind_120) = delx * Rn_3(Ind_020)
  Rn_4(Ind_021) = 0.0d0 ! Rn_4(Ind_021) = delz * Rn_3(Ind_020)
  Rn_4(Ind_102) = 0.0d0 ! Rn_4(Ind_102) = delx * Rn_3(Ind_002)
  Rn_4(Ind_012) = 0.0d0 ! Rn_4(Ind_012) = dely * Rn_3(Ind_002)
  Rn_4(Ind_111) = 0.0d0 ! Rn_4(Ind_111) = delx * Rn_3(Ind_011)
  Rn_4(Ind_400) = 3.d0 * B(n-2)  ! Rn_4(Ind_400) = 3.d0 * Rn_3(Ind_200) + delx * Rn_3(Ind_300)
  Rn_4(Ind_040) = 3.d0 * B(n-2)  ! Rn_4(Ind_040) = 3.d0 * Rn_3(Ind_020) + dely * Rn_3(Ind_030)
  Rn_4(Ind_004) = 3.d0 * B(n-2)  ! Rn_4(Ind_004) = 3.d0 * Rn_3(Ind_002) + delz * Rn_3(Ind_003)
  Rn_4(Ind_310) = 0.0d0 ! Rn_4(Ind_310) = dely * Rn_3(Ind_300)
  Rn_4(Ind_301) = 0.0d0 ! Rn_4(Ind_301) = delz * Rn_3(Ind_300)
  Rn_4(Ind_130) = 0.0d0 ! Rn_4(Ind_130) = delx * Rn_3(Ind_030)
  Rn_4(Ind_031) = 0.0d0 ! Rn_4(Ind_031) = delz * Rn_3(Ind_030) 
  Rn_4(Ind_103) = 0.0d0 ! Rn_4(Ind_103) = delx * Rn_3(Ind_003)
  Rn_4(Ind_013) = 0.0d0 ! Rn_4(Ind_013) = dely * Rn_3(Ind_003)
  Rn_4(Ind_220) = B(n-2)  ! Rn_4(Ind_220) = Rn_3(Ind_020) + delx * Rn_3(Ind_120)
  Rn_4(Ind_202) = B(n-2)  ! Rn_4(Ind_202) = Rn_3(Ind_002) + delx * Rn_3(Ind_102)
  Rn_4(Ind_022) = B(n-2)  ! Rn_4(Ind_022) = Rn_3(Ind_002) + dely * Rn_3(Ind_012)
  Rn_4(Ind_211) = 0.0d0 ! Rn_4(Ind_211) = dely * Rn_3(Ind_201)
  Rn_4(Ind_121) = 0.0d0 ! Rn_4(Ind_121) = delx * Rn_3(Ind_021)
  Rn_4(Ind_112) = 0.0d0 ! Rn_4(Ind_112) = delx * Rn_3(Ind_012)

!  B(0) = 2.d0 * ew_coeff / sqrt(PI)
!  B(1) = (4.d0/3.0d0) * ew_coeff * ew_coeff * ew_coeff / sqrt(PI)
!  B(2) = (8.d0/5.0d0) * ew_coeff * ew_coeff * ew_coeff * ew_coeff * ew_coeff / sqrt(PI)
  
  fScale = 2.d0 * ew_coeff / sqrt(PI)
  B(0) = fScale
  fScaleDip = (2.d0/3.0d0) * ew_coeff * ew_coeff
  B(1) = fScale * fScaleDip
  fScaleQdp = (4.d0/5.0d0) * ew_coeff * ew_coeff * ew_coeff * ew_coeff
  B(2) = fScale * fScaleQdp

  do img_id = my_img_lo, my_img_hi
  
    atm_id = gbl_img_atm_map(img_id)

    do j = 1, 10
       gmi(j) = global_multipole(j, atm_id)
    end do

    do j = 1, 3
       i_di(j) = ind_dip_d(j, atm_id)
       i_mi(j) = ind_dip_d(j, atm_id) + ind_dip_p(j, atm_id)
    end do

    ! self-field due to permanent mpoles at atm_id and derivs wrt r_j-r_i (at 0)
    
    phi(Ind_000) = B(0) * gmi(Ind_000)+ &
                   B(1)*( gmi(Ind_200) + gmi(Ind_020) + gmi(Ind_002))
    phi(Ind_100) = -B(1)*gmi(Ind_100)
    phi(Ind_010) = -B(1)*gmi(Ind_010)
    phi(Ind_001) = -B(1)*gmi(Ind_001)
    phi(Ind_200) = B(1)*gmi(Ind_000)+ &
                   B(2)* ( 3.0d0* gmi(Ind_200) + gmi(Ind_020) + gmi(Ind_002))
    phi(Ind_020) = B(1)*gmi(Ind_000)+ &
                   B(2)*( gmi(Ind_200)+3.d0*gmi(Ind_020)+ gmi(Ind_002)) 
    phi(Ind_002) = B(1)*gmi(Ind_000)+ &
                   B(2)*( gmi(Ind_200)+gmi(Ind_020)+ 3.d0*gmi(Ind_002)) 
    phi(Ind_110) = B(2)*gmi(Ind_110) 
    phi(Ind_101) = B(2)*gmi(Ind_101)
    phi(Ind_011) = B(2)*gmi(Ind_011)
    
    e_pp = phi(Ind_000)*gmi(Ind_000) + phi(Ind_100)*gmi(Ind_100) + &
           phi(Ind_010)*gmi(Ind_010) + phi(Ind_001)*gmi(Ind_001) + &
           phi(Ind_200)*gmi(Ind_200) + phi(Ind_020)*gmi(Ind_020) + &
           phi(Ind_002)*gmi(Ind_002) + phi(Ind_110)*gmi(Ind_110) + &
           phi(Ind_101)*gmi(Ind_101) + phi(Ind_011)*gmi(Ind_011)

    e_ind = phi(Ind_100)*i_di(1) + phi(Ind_010)*i_di(2) + &
            phi(Ind_001)*i_di(3)
    
    ene_perm = ene_perm - e_pp  ! subtract self-term.

    ene_ind = ene_ind - e_ind   ! subtract self-term.

    ! Torque due to permanent mpoles (should be zero).

    do j = 1, 10
       torque_field(j, atm_id) = torque_field(j, atm_id) - phi(j)
    end do

    ! Field due to induced dipoles.

    phi(Ind_000) = Rn_4(Ind_100)*i_mi(1) + Rn_4(Ind_010)*i_mi(2) + &
                   Rn_4(Ind_001)*i_mi(3)
    phi(Ind_100) = -(Rn_4(Ind_200)*i_mi(1) + Rn_4(Ind_110)*i_mi(2) + &
                     Rn_4(Ind_101)*i_mi(3))
    phi(Ind_010) = -(Rn_4(Ind_110)*i_mi(1) + Rn_4(Ind_020)*i_mi(2) + &
                     Rn_4(Ind_011)*i_mi(3))
    phi(Ind_001) = -(Rn_4(Ind_101)*i_mi(1) + Rn_4(Ind_011)*i_mi(2) + &
                     Rn_4(Ind_002)*i_mi(3))
    phi(Ind_200) = Rn_4(Ind_300)*i_mi(1) + Rn_4(Ind_210)*i_mi(2) + &
                   Rn_4(Ind_201)*i_mi(3)
    phi(Ind_020) = Rn_4(Ind_120)*i_mi(1) + Rn_4(Ind_030)*i_mi(2) + &
                   Rn_4(Ind_021)*i_mi(3)
    phi(Ind_002) = Rn_4(Ind_102)*i_mi(1) + Rn_4(Ind_012)*i_mi(2) + &
                   Rn_4(Ind_003)*i_mi(3)
    phi(Ind_110) = Rn_4(Ind_210)*i_mi(1) + Rn_4(Ind_120)*i_mi(2) + &
                   Rn_4(Ind_111)*i_mi(3)
    phi(Ind_101) = Rn_4(Ind_201)*i_mi(1) + Rn_4(Ind_111)*i_mi(2) + &
                   Rn_4(Ind_102)*i_mi(3)
    phi(Ind_011) = Rn_4(Ind_111)*i_mi(1) + Rn_4(Ind_021)*i_mi(2) + &
                   Rn_4(Ind_012)*i_mi(3)
    do j = 1, 10
       torque_field(j, atm_id) = torque_field(j, atm_id) - 0.5d0 * phi(j)
    end do

  end do

  ene_perm = 0.5d0 * coulomb_const_kcal_per_mole * ene_perm
  ene_ind = 0.5d0 * coulomb_const_kcal_per_mole * ene_ind

  call update_pme_time(pme_misc_timer)

  return

end subroutine am_self_ene_torque

end module amoeba_self_mod
