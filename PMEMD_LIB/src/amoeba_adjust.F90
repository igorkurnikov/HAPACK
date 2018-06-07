#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_adjust_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_adjust_mod

  implicit none

  private

#include "amoeba_mpole_index.i"

  integer, save ::      do_amoeba_adjust_flag, num_adjust_list

  integer, save, allocatable            :: adjust_list(:,:)

  ! This is a list of sublists analogous to atm_nb_mask(:) in the main
  ! pme code, but rather than containing atm id's, it contains the appropriate
  ! mutual_weight(:) index.

  integer, save, allocatable            :: adjust_weight_idxs(:)

! adjust list has 3rd component describing the relationship of the pair
! to each other: 9 possible values. Value is 1,2,3,4 if they are
! 1-2,1-3,1-4 or 1-5 pair; 5,6,7,8 if they are 1-2,1-3,1-4 or 1-5 pair
! within a polar group, or 9 if they are in a polar group but beyond 1-5

  double precision, save                :: vdw_weight(9)
  double precision, save                :: mpole_weight(9)
  double precision, save                :: direct_weight(9)
  double precision, save                :: polar_weight(9)
  double precision, save                :: mutual_weight(9)

  public        num_adjust_list
  public        am_adjust_zero_flag
  public        am_adjust_set_user_bit
  public        init_adjust_weight_idxs
  public        am_adjust_permfield
  public        am_adjust_dip_dip_fields
  public        am_adjust_ene_frc

contains

!*******************************************************************************!
! Subroutine:  am_adjust_zero_flag
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_adjust_zero_flag
  implicit none
  do_amoeba_adjust_flag = 0
  return
end subroutine am_adjust_zero_flag

!*******************************************************************************!
! Subroutine:  am_adjust_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_adjust_set_user_bit(do_this)

  use amoeba_flags_mod
  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in) :: do_this

  if (do_this .eq. 1) then      ! do in all cases
    do_amoeba_adjust_flag = ibset(do_amoeba_adjust_flag, user_induce_bit)
    do_amoeba_adjust_flag = ibset(do_amoeba_adjust_flag, user_postinduce_bit)
  else if (do_this .eq. 2) then ! do the induction, not the post-induction
    do_amoeba_adjust_flag = ibset(do_amoeba_adjust_flag, user_induce_bit)
    do_amoeba_adjust_flag = ibclr(do_amoeba_adjust_flag, user_postinduce_bit)
  else if (do_this .eq. 3) then ! do the post-induction, not the induction
    do_amoeba_adjust_flag = ibclr(do_amoeba_adjust_flag, user_induce_bit)
    do_amoeba_adjust_flag = ibset(do_amoeba_adjust_flag, user_postinduce_bit)
  else if (do_this .eq. 0) then 
    do_amoeba_adjust_flag = ibclr(do_amoeba_adjust_flag, user_induce_bit)
    do_amoeba_adjust_flag = ibclr(do_amoeba_adjust_flag, user_postinduce_bit)
  else
    error_msg = 'am_adjust_set_user_bit: bad value of user do_this'
    call mol_mech_error 
  end if

  return

end subroutine am_adjust_set_user_bit

!*******************************************************************************
!
! Subroutine:  init_adjust_weight_idxs
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine init_adjust_weight_idxs

  use prmtop_dat_mod, only : natom
  use pme_force_mod, only : atm_nb_maskdata, atm_nb_mask
  use file_io_dat_mod, only : mdout, error_msg
  use gbl_constants_mod, only : error_hdr

  implicit none

! Local variables:

  integer       :: adjust_lst_idx
  integer       :: atm_i, atm_j
  integer       :: found_cnt_i, found_cnt_j
  integer       :: wt_idx
  integer       :: mask_idx, mask_ptr, num_mask

  if (num_adjust_list .le. 0) return

  adjust_weight_idxs(:) = 0

  found_cnt_i = 0
  found_cnt_j = 0

  do adjust_lst_idx = 1, num_adjust_list

    atm_i = adjust_list(1, adjust_lst_idx)
    atm_j = adjust_list(2, adjust_lst_idx)
    wt_idx = adjust_list(3, adjust_lst_idx)

    num_mask = atm_nb_maskdata(atm_i)%nummask
    if (num_mask .gt. 0) then
      mask_ptr = atm_nb_maskdata(atm_i)%maskptr
      do mask_idx = mask_ptr + 1, mask_ptr + num_mask
        if (atm_nb_mask(mask_idx) .eq. atm_j) then
          if (adjust_weight_idxs(mask_idx) .ne. 0) then
            error_msg = 'Amoeba and Amber nonbonded adjust lists are inconsistent 1!'
            call mol_mech_error 
          end if
          adjust_weight_idxs(mask_idx) = wt_idx
          found_cnt_i = found_cnt_i + 1
          exit
        end if
      end do
    end if

    num_mask = atm_nb_maskdata(atm_j)%nummask
    if (num_mask .gt. 0) then
      mask_ptr = atm_nb_maskdata(atm_j)%maskptr
      do mask_idx = mask_ptr + 1, mask_ptr + num_mask
        if (atm_nb_mask(mask_idx) .eq. atm_i) then
          if (adjust_weight_idxs(mask_idx) .ne. 0) then
            error_msg = 'Amoeba and Amber nonbonded adjust lists are inconsistent 2!'
            call mol_mech_error 
          end if
          adjust_weight_idxs(mask_idx) = wt_idx
          found_cnt_j = found_cnt_j + 1
          exit
        end if
      end do
    end if

  end do

  if (found_cnt_i .ne. num_adjust_list .or. &
      found_cnt_j .ne. num_adjust_list) then
         error_msg = 'Amoeba and Amber nonbonded adjust lists are inconsistent 3!'
         call mol_mech_error 
  end if


  do mask_idx = 1, 2 * num_adjust_list
    if (adjust_weight_idxs(mask_idx) .eq. 0) then
        error_msg = 'Amoeba and Amber nonbonded adjust lists are inconsistent 5!'
        call mol_mech_error 
    end if
  end do

  return

end subroutine init_adjust_weight_idxs

!*******************************************************************************!
! Subroutine:  am_adjust_permfield
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_adjust_permfield(crd, direct_gradphi, polar_gradphi, &
                               dipole_dipole_tensor)

  use amoeba_multipoles_mod, only : global_multipole
  use amoeba_induced_mod, only : is_polarizable, screen_polar
  use prmtop_dat_mod, only : atm_qterm          ! is sq_polinv
  use mdin_amoeba_dat_mod, only : thole_expon_coeff
  use mdin_ewald_dat_mod
  use pme_direct_mod
  use timers_mod

  implicit none

! Formal Arguments:

  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: direct_gradphi(3, *)  ! corresponds to dip_field_d
  double precision, intent(in out)      :: polar_gradphi(3, *)   ! corresponds to dip_field_p
  double precision, intent(out)         :: dipole_dipole_tensor(6, *)

  call am_adjust_calc_permfield(crd, num_adjust_list, adjust_list, ew_coeff, &
                                eedtbdns, gbl_eed_cub, thole_expon_coeff, screen_polar, &
                                atm_qterm, is_polarizable, global_multipole, &
                                direct_weight, polar_weight, mutual_weight, &
                                direct_gradphi, polar_gradphi, &
                                dipole_dipole_tensor)

  call update_pme_time(adjust_masked_timer)

  return

end subroutine am_adjust_permfield

!*******************************************************************************!
! Subroutine:  am_adjust_calc_permfield
!
! Description: <TBS>
!
! here the difference between direct_gradphi ( dip_field_d )
! and                         polar_gradphi ( dip_field_p ) is originated 
!*******************************************************************************

subroutine  am_adjust_calc_permfield(crd, num_adjust_list, adjust_list, &
                                    ewaldcof, eedtbdns, eed_cub, &
                                    thole_expon_coeff, screen_polar, sq_polinv, &
                                    is_polarizable, global_multipole, &
                                    direct_weight, polar_weight, &
                                    mutual_weight, direct_gradphi, &
                                    polar_gradphi, dipole_dipole_tensor)

  use amoeba_flags_mod
  use file_io_dat_mod
  use parallel_dat_mod
  use img_mod

  implicit none

! Formal arguments:

  double precision, intent(in)          :: crd(3, *)
  integer, intent(in)                   :: num_adjust_list
  integer, intent(in)                   :: adjust_list(3, *)
  double precision, intent(in)          :: ewaldcof
  double precision, intent(in)          :: eedtbdns
  double precision, intent(in)          :: eed_cub(4, *)
  double precision, intent(in)          :: thole_expon_coeff
  double precision, intent(in)          :: screen_polar(*)
  double precision, intent(in)          :: sq_polinv(*)
  logical, intent(in)                   :: is_polarizable(*)
  double precision, intent(in)          :: global_multipole(10, *)
  double precision, intent(in)          :: direct_weight(*)
  double precision, intent(in)          :: polar_weight(*)
  double precision, intent(in)          :: mutual_weight(*)
  double precision, intent(in out)      :: direct_gradphi(3, *) ! corresponds to dip_field_d
  double precision, intent(in out)      :: polar_gradphi(3, *)  ! corresponds to dip_field_p
  double precision, intent(out)         :: dipole_dipole_tensor(6, *)

! Local variables:

  double precision                      :: delx, dely, delz
  double precision                      :: delr2, delr, delr2inv
  double precision                      :: x, dx, switch, d_switch_dx
  double precision                      :: dxdr
  double precision                      :: D(3), A(0:3), B(0:3), BD(3)
  double precision                      :: fac, fact, del, gphi_i(3), gphi_j(3)
  double precision                      :: expon, expo
  double precision                      :: clam3, clam5, clam7
  double precision                      :: delr3inv, delr5inv, delr7inv
  double precision                      :: Rn(1), Rn_1(4), Rn_2(10), Rn_3(20)
  integer                               :: i, j, k, n, n_adj, ind
  integer                               :: img_id

  double precision, parameter           :: third = 1.d0/3.d0
  logical                               :: single_thole_val

  if (iand(do_amoeba_adjust_flag, proceed_induce) .ne. proceed_induce) return

  fac = 2.d0 * ewaldcof * ewaldcof
  del = 1.d0 / eedtbdns               !  eedtbdns = 5000.0 - is it slope of switch function ?
  dxdr = ewaldcof
  
  if( thole_expon_coeff .gt. 0.000001) then 
    single_thole_val = .true.  
  else
    single_thole_val = .false.
  endif  

  do n_adj = 1, num_adjust_list
    if( numtasks .gt. 1) then
        img_id = gbl_atm_img_map(adjust_list(1, n_adj))
        if (img_id .lt. my_img_lo) cycle
        if (img_id .gt. my_img_hi) cycle
    endif 

    i = adjust_list(1, n_adj)
    j = adjust_list(2, n_adj)

    if (is_polarizable(i) .or. is_polarizable(j)) then
    
      k = adjust_list(3, n_adj)
      delx = crd(1, j) - crd(1, i)
      dely = crd(2, j) - crd(2, i)
      delz = crd(3, j) - crd(3, i)
      delr2 = delx * delx + dely * dely + delz * delz
      delr = sqrt(delr2)
      delr2inv = 1.d0 / delr2
      x = dxdr * delr

      ! Only supported option is cubic spline on switch...

      ind = eedtbdns * x + 1
      dx = x - (ind - 1.d0) * del
      switch = eed_cub(1, ind) + dx * (eed_cub(2, ind) + &
               dx * (eed_cub(3, ind) + dx * eed_cub(4, ind) * third) * 0.5d0)

      d_switch_dx = eed_cub(2, ind) + dx * (eed_cub(3, ind) + &
                      dx * eed_cub(4, ind) * 0.5d0)

      !------------------------------------------------------------
      ! McMurchie-Davidson recursion holds for any smooth function of r
      ! that is, to get the higher order derivs wrt x, y, z of g(r)
      ! define R(0, 0, 0, 0) = g(r)
      ! next  R(0, 0, 0, n+1) = -(1/r)d/dr R(0, 0, 0, n)
      ! then denote  R(t, u, v, n) = (d/dx)^t (d/dy)^u (d/dz)^v R(0, 0, 0, n)
      ! quantities of interest obtained by setting n = 0
      ! McMurchie-Davidson says that 
      !  R(t+1, u, v, n) = t * R(t-1, u, v, n+1) + x * R(t, u, v, n)
      !  R(t, u+1, v, n) = u * R(t, u-1, v, n+1) + y * R(t, u, v, n)
      !  R(t, u, v+1, n) = v * R(t, u, v-1, n+1) + z * R(t, u, v, n)
      ! use similar data structures as in direct sum
      !------------------------------------------------------------- 
      ! calc the contributions for cancelling reciprocal sum for ij
      ! -erf = erfc - 1.0 ---first get the erfc part as in direct sum

      B(0) = switch * delr * delr2inv
      fact = d_switch_dx * dxdr
      B(1) = (B(0) - fact) * delr2inv
      fact = fac * fact
      B(2) = (3.d0 * B(1) - fact) * delr2inv
      fact = fac * fact
      B(3) = (5.d0 * B(2) - fact) * delr2inv

      ! damping factors

      delr3inv = delr2inv / delr
      delr5inv = delr3inv * delr2inv
      delr7inv = delr5inv * delr2inv
      if( single_thole_val ) then
         expon = thole_expon_coeff * delr2 * delr * sq_polinv(i) * sq_polinv(j)
!         write(*,*) "am_adjust_calc_permfield()  thole_expon_coeff  = ", thole_expon_coeff 
      else                       
         expon = screen_polar(i)*screen_polar(j) * delr2 * delr * sq_polinv(i) * sq_polinv(j)
!         write(*,*) "am_adjust_calc_permfield()  screen_polar(i)*screen_polar(j) = ", screen_polar(i)*screen_polar(j)
      endif  
      expo = exp(-expon)

      ! clam3 = 1.d0-lam3, clam5 = 1.d0-lam5 etc. where 
      ! lam is from ponder's paper

      clam3 = expo
      clam5 = (1.d0 + expon) * expo
      clam7 = (1.d0 + expon + 0.6d0 * expon**2) * expo
      BD(1) = clam3 * delr3inv
      BD(2) = 3.d0 * clam5 * delr5inv
      BD(3) = 15.d0 * clam7 * delr7inv

      ! get the 1/r part
      ! note that damped coulomb part given by A(j) - BD(j) while 
      ! the erf part is given by A(j) - B(j)

      A(0) = 1.d0 * delr * delr2inv
      A(1) = A(0) * delr2inv
      A(2) = 3.d0 * A(1) * delr2inv
      A(3) = 5.d0 * A(2) * delr2inv

      ! first do the mutual field to store dipole-dipole tensor
      ! smaller recursion than for permanent field

      D(1) = -(mutual_weight(k) * (A(1) - BD(1)) + (B(1) - A(1)))
      D(2) = mutual_weight(k) * (A(2) - BD(2))  + (B(2) - A(2))
      n = 2
      Rn(Ind_000) = D(n)
      Rn_1(Ind_000) = D(n-1)
      Rn_1(Ind_100) = delx * Rn(Ind_000)
      Rn_1(Ind_010) = dely * Rn(Ind_000)
      Rn_1(Ind_001) = delz * Rn(Ind_000)
      !Rn_2(Ind_000) = D(n-2) no need for this
      !Rn_2(Ind_100) = delx*Rn_1(Ind_000) ditto
      !Rn_2(Ind_010) = dely*Rn_1(Ind_000) ditto
      !Rn_2(Ind_001) = delz*Rn_1(Ind_000) ditto
      Rn_2(Ind_200) = Rn_1(Ind_000) + delx * Rn_1(Ind_100)
      Rn_2(Ind_020) = Rn_1(Ind_000) + dely * Rn_1(Ind_010)
      Rn_2(Ind_002) = Rn_1(Ind_000) + delz * Rn_1(Ind_001)
      Rn_2(Ind_110) = delx * Rn_1(Ind_010)
      Rn_2(Ind_101) = delx * Rn_1(Ind_001)
      Rn_2(Ind_011) = dely * Rn_1(Ind_001)

      ! store dipole_dipole tensor 

      dipole_dipole_tensor(1,n_adj) = Rn_2(Ind_200)
      dipole_dipole_tensor(2,n_adj) = Rn_2(Ind_110)
      dipole_dipole_tensor(3,n_adj) = Rn_2(Ind_101)
      dipole_dipole_tensor(4,n_adj) = Rn_2(Ind_020)
      dipole_dipole_tensor(5,n_adj) = Rn_2(Ind_011)
      dipole_dipole_tensor(6,n_adj) = Rn_2(Ind_002)

      ! next do the direct field 

      D(1) = -(direct_weight(k) * (A(1) - BD(1)) + (B(1) - A(1)))
      D(2) = direct_weight(k) * (A(2) - BD(2))  + (B(2) - A(2))
      D(3) = -(direct_weight(k) * (A(3) - BD(3)) + (B(3) - A(3)))
      n = 3
      Rn(Ind_000) = D(n)
      Rn_1(Ind_000) = D(n-1)
      Rn_1(Ind_100) = delx * Rn(Ind_000)
      Rn_1(Ind_010) = dely * Rn(Ind_000)
      Rn_1(Ind_001) = delz * Rn(Ind_000)
      Rn_2(Ind_000) = D(n-2)
      Rn_2(Ind_100) = delx * Rn_1(Ind_000)
      Rn_2(Ind_010) = dely * Rn_1(Ind_000)
      Rn_2(Ind_001) = delz * Rn_1(Ind_000)
      Rn_2(Ind_200) = Rn_1(Ind_000) + delx * Rn_1(Ind_100)
      Rn_2(Ind_020) = Rn_1(Ind_000) + dely * Rn_1(Ind_010)
      Rn_2(Ind_002) = Rn_1(Ind_000) + delz * Rn_1(Ind_001)
      Rn_2(Ind_110) = delx * Rn_1(Ind_010)
      Rn_2(Ind_101) = delx * Rn_1(Ind_001)
      Rn_2(Ind_011) = dely * Rn_1(Ind_001)
      Rn_3(Ind_100) = delx * Rn_2(Ind_000)
      Rn_3(Ind_010) = dely * Rn_2(Ind_000)
      Rn_3(Ind_001) = delz * Rn_2(Ind_000)
      Rn_3(Ind_200) = Rn_2(Ind_000) + delx * Rn_2(Ind_100)
      Rn_3(Ind_020) = Rn_2(Ind_000) + dely * Rn_2(Ind_010)
      Rn_3(Ind_002) = Rn_2(Ind_000) + delz * Rn_2(Ind_001)
      Rn_3(Ind_110) = delx * Rn_2(Ind_010)
      Rn_3(Ind_101) = delx * Rn_2(Ind_001)
      Rn_3(Ind_011) = dely * Rn_2(Ind_001)
      Rn_3(Ind_300) = 2.d0 * Rn_2(Ind_100) + delx * Rn_2(Ind_200)
      Rn_3(Ind_030) = 2.d0 * Rn_2(Ind_010) + dely * Rn_2(Ind_020)
      Rn_3(Ind_003) = 2.d0 * Rn_2(Ind_001) + delz * Rn_2(Ind_002)
      Rn_3(Ind_210) = dely * Rn_2(Ind_200)
      Rn_3(Ind_201) = delz * Rn_2(Ind_200)
      Rn_3(Ind_120) = delx * Rn_2(Ind_020)
      Rn_3(Ind_021) = delz * Rn_2(Ind_020)
      Rn_3(Ind_102) = delx * Rn_2(Ind_002)
      Rn_3(Ind_012) = dely * Rn_2(Ind_002)
      Rn_3(Ind_111) = delx * Rn_2(Ind_011)

      ! finally get field components

      if (is_polarizable(i)) then

        gphi_i(1) = Rn_3(Ind_100) * global_multipole(Ind_000,j) + &
                    Rn_3(Ind_200) * global_multipole(Ind_100,j) + &
                    Rn_3(Ind_110) * global_multipole(Ind_010,j) + &
                    Rn_3(Ind_101) * global_multipole(Ind_001,j) + &
                    Rn_3(Ind_300) * global_multipole(Ind_200,j) + &
                    Rn_3(Ind_120) * global_multipole(Ind_020,j) + &
                    Rn_3(Ind_102) * global_multipole(Ind_002,j) + &
                    Rn_3(Ind_210) * global_multipole(Ind_110,j) + &
                    Rn_3(Ind_201) * global_multipole(Ind_101,j) + &
                    Rn_3(Ind_111) * global_multipole(Ind_011,j)
        gphi_i(2) = Rn_3(Ind_010) * global_multipole(Ind_000,j) + &
                    Rn_3(Ind_110) * global_multipole(Ind_100,j) + &
                    Rn_3(Ind_020) * global_multipole(Ind_010,j) + &
                    Rn_3(Ind_011) * global_multipole(Ind_001,j) + &
                    Rn_3(Ind_210) * global_multipole(Ind_200,j) + &
                    Rn_3(Ind_030) * global_multipole(Ind_020,j) + &
                    Rn_3(Ind_012) * global_multipole(Ind_002,j) + &
                    Rn_3(Ind_120) * global_multipole(Ind_110,j) + &
                    Rn_3(Ind_111) * global_multipole(Ind_101,j) + &
                    Rn_3(Ind_021) * global_multipole(Ind_011,j)
        gphi_i(3) = Rn_3(Ind_001) * global_multipole(Ind_000,j) + &
                    Rn_3(Ind_101) * global_multipole(Ind_100,j) + &
                    Rn_3(Ind_011) * global_multipole(Ind_010,j) + &
                    Rn_3(Ind_002) * global_multipole(Ind_001,j) + &
                    Rn_3(Ind_201) * global_multipole(Ind_200,j) + &
                    Rn_3(Ind_021) * global_multipole(Ind_020,j) + &
                    Rn_3(Ind_003) * global_multipole(Ind_002,j) + &
                    Rn_3(Ind_111) * global_multipole(Ind_110,j) + &
                    Rn_3(Ind_102) * global_multipole(Ind_101,j) + &
                    Rn_3(Ind_012) * global_multipole(Ind_011,j)

        ! negative for i since d/dx_i = -d/delx

        direct_gradphi(1,i) = direct_gradphi(1,i) - gphi_i(1)
        direct_gradphi(2,i) = direct_gradphi(2,i) - gphi_i(2)
        direct_gradphi(3,i) = direct_gradphi(3,i) - gphi_i(3)

      end if ! is_polarizable(i)

      if (is_polarizable(j)) then

        ! note negative contribs for dipoles of i since d/dx_i = -d/delx)
        ! (look at contributions to electrostatic potential at j -these will
        !  contain negative contributions due to dipolar contribs at i--the
        !  extra deriv in tensor due to grad at j has no negative signs)

        gphi_j(1) = Rn_3(Ind_100) * global_multipole(Ind_000,i) - &
                    Rn_3(Ind_200) * global_multipole(Ind_100,i) - &
                    Rn_3(Ind_110) * global_multipole(Ind_010,i) - &
                    Rn_3(Ind_101) * global_multipole(Ind_001,i) + &
                    Rn_3(Ind_300) * global_multipole(Ind_200,i) + &
                    Rn_3(Ind_120) * global_multipole(Ind_020,i) + &
                    Rn_3(Ind_102) * global_multipole(Ind_002,i) + &
                    Rn_3(Ind_210) * global_multipole(Ind_110,i) + &
                    Rn_3(Ind_201) * global_multipole(Ind_101,i) + &
                    Rn_3(Ind_111) * global_multipole(Ind_011,i)
        gphi_j(2) = Rn_3(Ind_010) * global_multipole(Ind_000,i) - &
                    Rn_3(Ind_110) * global_multipole(Ind_100,i) - &
                    Rn_3(Ind_020) * global_multipole(Ind_010,i) - &
                    Rn_3(Ind_011) * global_multipole(Ind_001,i) + &
                    Rn_3(Ind_210) * global_multipole(Ind_200,i) + &
                    Rn_3(Ind_030) * global_multipole(Ind_020,i) + &
                    Rn_3(Ind_012) * global_multipole(Ind_002,i) + &
                    Rn_3(Ind_120) * global_multipole(Ind_110,i) + &
                    Rn_3(Ind_111) * global_multipole(Ind_101,i) + &
                    Rn_3(Ind_021) * global_multipole(Ind_011,i)
        gphi_j(3) = Rn_3(Ind_001) * global_multipole(Ind_000,i) - &
                    Rn_3(Ind_101) * global_multipole(Ind_100,i) - &
                    Rn_3(Ind_011) * global_multipole(Ind_010,i) - &
                    Rn_3(Ind_002) * global_multipole(Ind_001,i) + &
                    Rn_3(Ind_201) * global_multipole(Ind_200,i) + &
                    Rn_3(Ind_021) * global_multipole(Ind_020,i) + &
                    Rn_3(Ind_003) * global_multipole(Ind_002,i) + &
                    Rn_3(Ind_111) * global_multipole(Ind_110,i) + &
                    Rn_3(Ind_102) * global_multipole(Ind_101,i) + &
                    Rn_3(Ind_012) * global_multipole(Ind_011,i)

        direct_gradphi(1,j) = direct_gradphi(1,j) + gphi_j(1)
        direct_gradphi(2,j) = direct_gradphi(2,j) + gphi_j(2)
        direct_gradphi(3,j) = direct_gradphi(3,j) + gphi_j(3)

      end if !(is_polarizable(j)) then

      ! next do the polar field --negate the odd order terms

      D(1) = -(polar_weight(k) * (A(1) - BD(1)) + (B(1) - A(1)))
      D(2) = polar_weight(k) * (A(2) - BD(2))  + (B(2) - A(2))
      D(3) = -(polar_weight(k) * (A(3) - BD(3)) + (B(3) - A(3)))
      n = 3
      Rn(Ind_000) = D(n)
      Rn_1(Ind_000) = D(n-1)
      Rn_1(Ind_100) = delx * Rn(Ind_000)
      Rn_1(Ind_010) = dely * Rn(Ind_000)
      Rn_1(Ind_001) = delz * Rn(Ind_000)
      Rn_2(Ind_000) = D(n-2)
      Rn_2(Ind_100) = delx * Rn_1(Ind_000)
      Rn_2(Ind_010) = dely * Rn_1(Ind_000)
      Rn_2(Ind_001) = delz * Rn_1(Ind_000)
      Rn_2(Ind_200) = Rn_1(Ind_000) + delx * Rn_1(Ind_100)
      Rn_2(Ind_020) = Rn_1(Ind_000) + dely * Rn_1(Ind_010)
      Rn_2(Ind_002) = Rn_1(Ind_000) + delz * Rn_1(Ind_001)
      Rn_2(Ind_110) = delx * Rn_1(Ind_010)
      Rn_2(Ind_101) = delx * Rn_1(Ind_001)
      Rn_2(Ind_011) = dely * Rn_1(Ind_001)
      Rn_3(Ind_100) = delx * Rn_2(Ind_000)
      Rn_3(Ind_010) = dely * Rn_2(Ind_000)
      Rn_3(Ind_001) = delz * Rn_2(Ind_000)
      Rn_3(Ind_200) = Rn_2(Ind_000) + delx * Rn_2(Ind_100)
      Rn_3(Ind_020) = Rn_2(Ind_000) + dely * Rn_2(Ind_010)
      Rn_3(Ind_002) = Rn_2(Ind_000) + delz * Rn_2(Ind_001)
      Rn_3(Ind_110) = delx * Rn_2(Ind_010)
      Rn_3(Ind_101) = delx * Rn_2(Ind_001)
      Rn_3(Ind_011) = dely * Rn_2(Ind_001)
      Rn_3(Ind_300) = 2.d0 * Rn_2(Ind_100) + delx * Rn_2(Ind_200)
      Rn_3(Ind_030) = 2.d0 * Rn_2(Ind_010) + dely * Rn_2(Ind_020)
      Rn_3(Ind_003) = 2.d0 * Rn_2(Ind_001) + delz * Rn_2(Ind_002)
      Rn_3(Ind_210) = dely * Rn_2(Ind_200)
      Rn_3(Ind_201) = delz * Rn_2(Ind_200)
      Rn_3(Ind_120) = delx * Rn_2(Ind_020)
      Rn_3(Ind_021) = delz * Rn_2(Ind_020)
      Rn_3(Ind_102) = delx * Rn_2(Ind_002)
      Rn_3(Ind_012) = dely * Rn_2(Ind_002)
      Rn_3(Ind_111) = delx * Rn_2(Ind_011)

      ! finally get field components 

      if (is_polarizable(i)) then

        gphi_i(1) = Rn_3(Ind_100) * global_multipole(Ind_000,j) + &
                    Rn_3(Ind_200) * global_multipole(Ind_100,j) + &
                    Rn_3(Ind_110) * global_multipole(Ind_010,j) + &
                    Rn_3(Ind_101) * global_multipole(Ind_001,j) + &
                    Rn_3(Ind_300) * global_multipole(Ind_200,j) + &
                    Rn_3(Ind_120) * global_multipole(Ind_020,j) + &
                    Rn_3(Ind_102) * global_multipole(Ind_002,j) + &
                    Rn_3(Ind_210) * global_multipole(Ind_110,j) + &
                    Rn_3(Ind_201) * global_multipole(Ind_101,j) + &
                    Rn_3(Ind_111) * global_multipole(Ind_011,j)
        gphi_i(2) = Rn_3(Ind_010) * global_multipole(Ind_000,j) + &
                    Rn_3(Ind_110) * global_multipole(Ind_100,j) + &
                    Rn_3(Ind_020) * global_multipole(Ind_010,j) + &
                    Rn_3(Ind_011) * global_multipole(Ind_001,j) + &
                    Rn_3(Ind_210) * global_multipole(Ind_200,j) + &
                    Rn_3(Ind_030) * global_multipole(Ind_020,j) + &
                    Rn_3(Ind_012) * global_multipole(Ind_002,j) + &
                    Rn_3(Ind_120) * global_multipole(Ind_110,j) + &
                    Rn_3(Ind_111) * global_multipole(Ind_101,j) + &
                    Rn_3(Ind_021) * global_multipole(Ind_011,j)
        gphi_i(3) = Rn_3(Ind_001) * global_multipole(Ind_000,j) + &
                    Rn_3(Ind_101) * global_multipole(Ind_100,j) + &
                    Rn_3(Ind_011) * global_multipole(Ind_010,j) + &
                    Rn_3(Ind_002) * global_multipole(Ind_001,j) + &
                    Rn_3(Ind_201) * global_multipole(Ind_200,j) + &
                    Rn_3(Ind_021) * global_multipole(Ind_020,j) + &
                    Rn_3(Ind_003) * global_multipole(Ind_002,j) + &
                    Rn_3(Ind_111) * global_multipole(Ind_110,j) + &
                    Rn_3(Ind_102) * global_multipole(Ind_101,j) + &
                    Rn_3(Ind_012) * global_multipole(Ind_011,j)

        ! negative for i since d/dx_i = -d/delx

        polar_gradphi(1,i) = polar_gradphi(1,i) - gphi_i(1)
        polar_gradphi(2,i) = polar_gradphi(2,i) - gphi_i(2)
        polar_gradphi(3,i) = polar_gradphi(3,i) - gphi_i(3)

      end if ! is_polarizable(i)

      if (is_polarizable(j)) then

         ! note negative contribs for dipoles of i since d/dx_i = -d/delx)
         ! (look at contributions to electrostatic potential at j -these will
         !  contain negative contributions due to dipolar contribs at i--the
         !  extra deriv in tensor due to grad at j has no negative signs)
         gphi_j(1) = Rn_3(Ind_100) * global_multipole(Ind_000,i) - &
                     Rn_3(Ind_200) * global_multipole(Ind_100,i) - &
                     Rn_3(Ind_110) * global_multipole(Ind_010,i) - &
                     Rn_3(Ind_101) * global_multipole(Ind_001,i) + &
                     Rn_3(Ind_300) * global_multipole(Ind_200,i) + &
                     Rn_3(Ind_120) * global_multipole(Ind_020,i) + &
                     Rn_3(Ind_102) * global_multipole(Ind_002,i) + &
                     Rn_3(Ind_210) * global_multipole(Ind_110,i) + &
                     Rn_3(Ind_201) * global_multipole(Ind_101,i) + &
                     Rn_3(Ind_111) * global_multipole(Ind_011,i)
         gphi_j(2) = Rn_3(Ind_010) * global_multipole(Ind_000,i) - &
         
                     Rn_3(Ind_110) * global_multipole(Ind_100,i) - &
                     Rn_3(Ind_020) * global_multipole(Ind_010,i) - &
                     Rn_3(Ind_011) * global_multipole(Ind_001,i) + &
                     Rn_3(Ind_210) * global_multipole(Ind_200,i) + &
                     Rn_3(Ind_030) * global_multipole(Ind_020,i) + &
                     Rn_3(Ind_012) * global_multipole(Ind_002,i) + &
                     Rn_3(Ind_120) * global_multipole(Ind_110,i) + &
                     Rn_3(Ind_111) * global_multipole(Ind_101,i) + &
                     Rn_3(Ind_021) * global_multipole(Ind_011,i)
         gphi_j(3) = Rn_3(Ind_001) * global_multipole(Ind_000,i) - &
                     Rn_3(Ind_101) * global_multipole(Ind_100,i) - &
                     Rn_3(Ind_011) * global_multipole(Ind_010,i) - &
                     Rn_3(Ind_002) * global_multipole(Ind_001,i) + &
                     Rn_3(Ind_201) * global_multipole(Ind_200,i) + &
                     Rn_3(Ind_021) * global_multipole(Ind_020,i) + &
                     Rn_3(Ind_003) * global_multipole(Ind_002,i) + &
                     Rn_3(Ind_111) * global_multipole(Ind_110,i) + &
                     Rn_3(Ind_102) * global_multipole(Ind_101,i) + &
                     Rn_3(Ind_012) * global_multipole(Ind_011,i)
         polar_gradphi(1,j) = polar_gradphi(1,j) + gphi_j(1)
         polar_gradphi(2,j) = polar_gradphi(2,j) + gphi_j(2)
         polar_gradphi(3,j) = polar_gradphi(3,j) + gphi_j(3)

      end if !(is_polarizable(j)) then

    end if ! (is_polarizable(i) .or. is_polarizable(j)) then

  end do

  return

end subroutine am_adjust_calc_permfield

!*******************************************************************************!
! Subroutine:  am_adjust_dip_dip_fields
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_adjust_dip_dip_fields(dipole_dipole_tensor, &
                                    ind_dip_d, ind_dip_p, &
                                    dip_field_d, dip_field_p)

  use amoeba_flags_mod
  use timers_mod

  implicit none

! Formal arguments:

  double precision, intent(in)          :: dipole_dipole_tensor(3, *)
  double precision, intent(in)          :: ind_dip_d(3, *)
  double precision, intent(in)          :: ind_dip_p(3, *)
  double precision, intent(in out)      :: dip_field_d(3, *)
  double precision, intent(in out)      :: dip_field_p(3, *) 

  if (iand(do_amoeba_adjust_flag, proceed_induce) .ne. proceed_induce) return

  call am_adjust_calc_dip_dip_fields(num_adjust_list, adjust_list, &
                                     dipole_dipole_tensor, ind_dip_d, &
                                     ind_dip_p, dip_field_d, dip_field_p)

  call update_pme_time(adjust_masked_timer)

  return

end subroutine am_adjust_dip_dip_fields

!*******************************************************************************!
! Subroutine:  am_adjust_calc_dip_dip_fields
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_adjust_calc_dip_dip_fields(num_adjust_list, adjust_list, &
                                         dipole_dipole_tensor, ind_dip_d, &
                                         ind_dip_p, dip_field_d, dip_field_p)

  use img_mod
  use parallel_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: num_adjust_list
  integer, intent(in)                   :: adjust_list(3, *)
  double precision, intent(in)          :: dipole_dipole_tensor(6, *)
  double precision, intent(in)          :: ind_dip_d(3, *)
  double precision, intent(in)          :: ind_dip_p(3, *)
  double precision, intent(in out)      :: dip_field_d(3, *)
  double precision, intent(in out)      :: dip_field_p(3, *) 

! Local variables:

  integer                               :: i, j, n
  integer                               :: img_id


  do n = 1, num_adjust_list

    if(numtasks .gt. 1) then
        img_id = gbl_atm_img_map(adjust_list(1, n))
        if (img_id .lt. my_img_lo) cycle
        if (img_id .gt. my_img_hi) cycle
    endif

    i = adjust_list(1, n)
    j = adjust_list(2, n)

    ! minus signs due to deriv wrt position of i

    dip_field_d(1, i) = dip_field_d(1, i) -  &
      (dipole_dipole_tensor(1, n) * ind_dip_d(1, j) + &
       dipole_dipole_tensor(2, n) * ind_dip_d(2, j) + &
       dipole_dipole_tensor(3, n) * ind_dip_d(3, j))

    dip_field_d(2, i) = dip_field_d(2, i) -  &
      (dipole_dipole_tensor(2, n) * ind_dip_d(1, j) + &
       dipole_dipole_tensor(4, n) * ind_dip_d(2, j) + &
       dipole_dipole_tensor(5, n) * ind_dip_d(3, j))

    dip_field_d(3, i) = dip_field_d(3, i) -  &
      (dipole_dipole_tensor(3, n) * ind_dip_d(1, j) + &
       dipole_dipole_tensor(5, n) * ind_dip_d(2, j) + &
       dipole_dipole_tensor(6, n) * ind_dip_d(3, j))

    dip_field_d(1, j) = dip_field_d(1, j) -  &
      (dipole_dipole_tensor(1, n) * ind_dip_d(1, i) + &
       dipole_dipole_tensor(2, n) * ind_dip_d(2, i) + &
       dipole_dipole_tensor(3, n) * ind_dip_d(3, i))

    dip_field_d(2, j) = dip_field_d(2, j) -  &
      (dipole_dipole_tensor(2, n) * ind_dip_d(1, i) + &
       dipole_dipole_tensor(4, n) * ind_dip_d(2, i) + &
       dipole_dipole_tensor(5, n) * ind_dip_d(3, i))

    dip_field_d(3, j) = dip_field_d(3, j) -  &
      (dipole_dipole_tensor(3, n) * ind_dip_d(1, i) + &
       dipole_dipole_tensor(5, n) * ind_dip_d(2, i) + &
       dipole_dipole_tensor(6, n) * ind_dip_d(3, i))

    ! now do other set of dipoles

    dip_field_p(1, i) = dip_field_p(1, i) -  &
      (dipole_dipole_tensor(1, n) * ind_dip_p(1, j) + &
       dipole_dipole_tensor(2, n) * ind_dip_p(2, j) + &
       dipole_dipole_tensor(3, n) * ind_dip_p(3, j))

    dip_field_p(2, i) = dip_field_p(2, i) -  &
      (dipole_dipole_tensor(2, n) * ind_dip_p(1, j) + &
       dipole_dipole_tensor(4, n) * ind_dip_p(2, j) + &
       dipole_dipole_tensor(5, n) * ind_dip_p(3, j))

    dip_field_p(3, i) = dip_field_p(3, i) -  &
      (dipole_dipole_tensor(3, n) * ind_dip_p(1, j) + &
       dipole_dipole_tensor(5, n) * ind_dip_p(2, j) + &
       dipole_dipole_tensor(6, n) * ind_dip_p(3, j))

    dip_field_p(1, j) = dip_field_p(1, j) -  &
      (dipole_dipole_tensor(1, n) * ind_dip_p(1, i) + &
       dipole_dipole_tensor(2, n) * ind_dip_p(2, i) + &
       dipole_dipole_tensor(3, n) * ind_dip_p(3, i))

    dip_field_p(2, j) = dip_field_p(2, j) -  &
      (dipole_dipole_tensor(2, n) * ind_dip_p(1, i) + &
       dipole_dipole_tensor(4, n) * ind_dip_p(2, i) + &
       dipole_dipole_tensor(5, n) * ind_dip_p(3, i))

    dip_field_p(3, j) = dip_field_p(3, j) -  &
      (dipole_dipole_tensor(3, n) * ind_dip_p(1, i) + &
       dipole_dipole_tensor(5, n) * ind_dip_p(2, i) + &
       dipole_dipole_tensor(6, n) * ind_dip_p(3, i))

  end do

  return

end subroutine am_adjust_calc_dip_dip_fields

!*******************************************************************************!
! Subroutine:  am_adjust_ene_frc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_adjust_ene_frc(crd, ind_dip_d, ind_dip_p, &
                             ene_perm, ene_ind, ene_vdw, frc, virial, atm_owner_map )

  use amoeba_multipoles_mod, only : global_multipole
  use amoeba_induced_mod, only : is_polarizable, screen_polar
  use prmtop_dat_mod, only : atm_qterm          ! is sq_polinv
  use mdin_amoeba_dat_mod, only : thole_expon_coeff
  use amoeba_vdw_mod, only : am_vdw_adjust_ene_frc
  use mdin_ewald_dat_mod
  use pme_direct_mod
  use timers_mod
  
  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in)          :: ind_dip_d(3, *)
  double precision, intent(in)          :: ind_dip_p(3, *)
  double precision, intent(in out)      :: ene_perm
  double precision, intent(in out)      :: ene_ind
  double precision, intent(in out)      :: ene_vdw
  double precision, intent(in out)      :: frc(3, *)
  double precision, intent(in out)      :: virial(3, 3)
  integer, intent(in out)               :: atm_owner_map(*)

  ene_perm = 0.d0
  ene_ind = 0.d0
  ene_vdw = 0.d0

  call am_vdw_adjust_ene_frc(crd, atm_owner_map, num_adjust_list, &
                             adjust_list, vdw_weight, ene_vdw, frc, virial)

  call am_adjust_calc_ene_frc(crd, atm_owner_map, num_adjust_list, &
                              adjust_list, ew_coeff, eedtbdns, gbl_eed_cub, &
                              thole_expon_coeff, screen_polar, atm_qterm, is_polarizable, &
                              global_multipole, ind_dip_d, ind_dip_p, &
                              mpole_weight, polar_weight, direct_weight, &
                              mutual_weight, ene_perm, ene_ind, frc, virial)

  call update_pme_time(adjust_masked_timer)

  return

end subroutine am_adjust_ene_frc

!*******************************************************************************!
! Subroutine:  am_adjust_calc_ene_frc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_adjust_calc_ene_frc(crd, atm_owner_map, num_adjust_list, &
                                  adjust_list, ewaldcof, eedtbdns, eed_cub, &
                                  thole_expon_coeff, screen_polar, sq_polinv, is_polarizable,&
                                  global_multipole, ind_dip_d, ind_dip_p, &
                                  mpole_weight, polar_weight, direct_weight, &
                                  mutual_weight, ene_perm, ene_ind, frc, virial)

  use amoeba_multipoles_mod, only : coulomb_const_kcal_per_mole, torque_field
  use amoeba_flags_mod
  use file_io_dat_mod
  use parallel_dat_mod
  use img_mod

  implicit none

! Formal arguments:

  double precision, intent(in)          :: crd(3, *)
  integer, intent(in)                   :: atm_owner_map(*)
  integer, intent(in)                   :: num_adjust_list
  integer, intent(in)                   :: adjust_list(3, *)
  double precision, intent(in)          :: ewaldcof
  double precision, intent(in)          :: eedtbdns
  double precision, intent(in)          :: eed_cub(4, *)
  double precision, intent(in)          :: thole_expon_coeff
  double precision, intent(in)          :: screen_polar(*)
  double precision, intent(in)          :: sq_polinv(*)
  logical, intent(in)                   :: is_polarizable(*)
  double precision, intent(in)          :: global_multipole(10, *)
  double precision, intent(in)          :: ind_dip_d(3, *)
  double precision, intent(in)          :: ind_dip_p(3, *)
  double precision, intent(in)          :: mpole_weight(*)
  double precision, intent(in)          :: polar_weight(*)
  double precision, intent(in)          :: direct_weight(*)
  double precision, intent(in)          :: mutual_weight(*)
  double precision, intent(in out)      :: ene_perm
  double precision, intent(in out)      :: ene_ind
  double precision, intent(in out)      :: frc(3, *)
  double precision, intent(in out)      :: virial(3, 3)

! Local variables:

  double precision                      :: delx, dely, delz
  double precision                      :: delr, delr2, delr2inv
  double precision                      :: x, dx, switch, d_switch_dx
  double precision                      :: dxdr
  double precision                      :: vxx, vxy, vxz
  double precision                      :: vyx, vyy, vyz
  double precision                      :: vzx, vzy, vzz
  double precision                      :: D(0:5), A(0:5), B(0:5), C(0:5), BD(4)
  double precision                      :: fac, fact, del
  double precision                      :: expon, expo
  double precision                      :: clam3, clam5, clam7, clam9
  double precision                      :: delr3inv, delr5inv
  double precision                      :: delr7inv, delr9inv
  double precision                      :: Rn(1), Rn_1(4), Rn_2(10)
  double precision                      :: Rn_3(20), Rn_4(35), Rn_5(56)
  double precision                      :: gmi(10), gmj(10)
  double precision                      :: phi(20)
  double precision                      :: e_pp
  double precision                      :: g_pp(3)
  double precision                      :: e_ind
  double precision                      :: e_add
  double precision                      :: g_ind(3)
  double precision                      :: i_di(3), i_dj(3)
  double precision                      :: i_pi(3), i_pj(3)

  double precision, parameter           :: third = 1.d0 / 3.d0
  double precision, parameter           :: const1 = 0.6d0
  double precision, parameter           :: const2 = 18.d0 / 35.d0
  double precision, parameter           :: const3 = 9.d0 / 35.d0
          
  integer                               :: i, j, k, n, n_adj, ind, jj
  integer                               :: img_id
  logical                               :: single_thole_val

  if (iand(do_amoeba_adjust_flag, proceed_postinduce) .ne. &
      proceed_postinduce) return

  fac = 2.d0 * ewaldcof * ewaldcof
  del = 1.d0 / eedtbdns
  dxdr = ewaldcof
  
  if( thole_expon_coeff .gt. 0.000001) then 
    single_thole_val = .true.
  else
    single_thole_val = .false.
  endif  
  
  vxx = 0.d0
  vxy = 0.d0
  vxz = 0.d0
  vyx = 0.d0
  vyy = 0.d0
  vyz = 0.d0
  vzx = 0.d0
  vzy = 0.d0
  vzz = 0.d0
  

  do n_adj = 1, num_adjust_list
    if( numtasks .gt. 1) then
        img_id = gbl_atm_img_map(adjust_list(1, n_adj))
        if (img_id .lt. my_img_lo) cycle
        if (img_id .gt. my_img_hi) cycle
    endif 
    i = adjust_list(1, n_adj)
    j = adjust_list(2, n_adj)
    k = adjust_list(3, n_adj)
    delx = crd(1, j) - crd(1, i)
    dely = crd(2, j) - crd(2, i)
    delz = crd(3, j) - crd(3, i)
    delr2 = delx * delx + dely * dely + delz * delz
    delr = sqrt(delr2)
    delr2inv = 1.d0 / delr2
    x = dxdr * delr

    ! Only supported option is cubic spline on switch...

    ind = eedtbdns * x + 1
    dx = x - (ind - 1.d0) * del

    switch = eed_cub(1, ind) + dx * (eed_cub(2, ind) + &
             dx * (eed_cub(3, ind) + dx * eed_cub(4, ind) * third) * 0.5d0)

    d_switch_dx = eed_cub(2, ind) + dx * (eed_cub(3, ind) + &
                  dx * eed_cub(4, ind) * 0.5d0)

    !------------------------------------------------------------
    ! McMurchie-Davidson recursion holds for any smooth function of r
    ! that is, to get the higher order derivs wrt x, y, z of g(r)
    ! define R(0, 0, 0, 0) = g(r)
    ! next  R(0, 0, 0, n+1) = -(1/r)d/dr R(0, 0, 0, n)
    ! then denote  R(t, u, v, n) = (d/dx)^t (d/dy)^u (d/dz)^v R(0, 0, 0, n)
    ! quantities of interest obtained by setting n = 0
    ! McMurchie-Davidson says that 
    !  R(t+1, u, v, n) = t * R(t-1, u, v, n+1) + x * R(t, u, v, n)
    !  R(t, u+1, v, n) = u * R(t, u-1, v, n+1) + y * R(t, u, v, n)
    !  R(t, u, v+1, n) = v * R(t, u, v-1, n+1) + z * R(t, u, v, n)
    ! use similar data structures as in direct sum
    !------------------------------------------------------------- 
    ! calc the contributions for cancelling reciprocal sum for ij
    ! -erf = erfc - 1.0 ---first get the erfc part as in direct sum

    B(0) = switch * delr * delr2inv
    fact = d_switch_dx * dxdr
    B(1) = (B(0) - fact) * delr2inv
    fact = fac * fact
    B(2) = (3.d0 * B(1) - fact) * delr2inv
    fact = fac * fact
    B(3) = (5.d0 * B(2) - fact) * delr2inv
    fact = fac * fact
    B(4) = (7.d0 * B(3) - fact) * delr2inv
    fact = fac * fact
    B(5) = (9.d0 * B(4) - fact) * delr2inv

    ! the -erf part is given by B(j) - A(j); weighted direct by
    ! mpole_weight * A(j)

    A(0) = 1.d0 * delr * delr2inv
    A(1) = A(0) * delr2inv
    A(2) = 3.d0 * A(1) * delr2inv
    A(3) = 5.d0 * A(2) * delr2inv
    A(4) = 7.d0 * A(3) * delr2inv
    A(5) = 9.d0 * A(4) * delr2inv

    ! damped part for permanent-induced or induced-induced interactions

    delr3inv = delr2inv / delr
    delr5inv = delr3inv * delr2inv
    delr7inv = delr5inv * delr2inv
    delr9inv = delr7inv * delr2inv
    if( single_thole_val ) then
         expon = thole_expon_coeff * delr2 * delr * sq_polinv(i) * sq_polinv(j)
!         write(*,*) " am_adjust_calc_ene_frc  thole_expon_coeff = ", thole_expon_coeff
    else                       
         expon = screen_polar(i)*screen_polar(j) * delr2 * delr * sq_polinv(i) * sq_polinv(j)
!         write(*,*) " am_adjust_calc_ene_frc  screen_polar(i)*screen_polar(j) = ", screen_polar(i)*screen_polar(j)
    endif  
    expo = exp(-expon)

    ! clam3 = 1.d0-lam3, clam5 = 1.d0-lam5 etc. where 
    ! lam is from ponder's paper

    clam3 = expo
    clam5 = (1.d0 + expon) * expo
    clam7 = (1.d0 + expon + const1 * expon**2) * expo
    clam9 = (1.d0 + expon + const2 * expon**2 + const3 * expon**3) * expo
    BD(1) = clam3 * delr3inv
    BD(2) = 3.d0 * clam5 * delr5inv
    BD(3) = 15.d0 * clam7 * delr7inv
    BD(4) = 105.d0 * clam9 * delr9inv

    ! first handle the permanent-permanent interactions

    do jj = 0, 5
       C(jj) = (B(jj) - A(jj)) + mpole_weight(k) * A(jj)
    end do

    ! negate the odd order to get sign right

    C(1) = -C(1)
    C(3) = -C(3)
    C(5) = -C(5)

    ! get the interaction tensor 

    n = 5
    Rn(Ind_000) = C(n)
    Rn_1(Ind_000) = C(n-1)
    Rn_1(Ind_100) = delx * Rn(Ind_000)
    Rn_1(Ind_010) = dely * Rn(Ind_000)
    Rn_1(Ind_001) = delz * Rn(Ind_000)
    Rn_2(Ind_000) = C(n-2)
    Rn_2(Ind_100) = delx * Rn_1(Ind_000)
    Rn_2(Ind_010) = dely * Rn_1(Ind_000)
    Rn_2(Ind_001) = delz * Rn_1(Ind_000)
    Rn_2(Ind_200) = Rn_1(Ind_000) + delx * Rn_1(Ind_100)
    Rn_2(Ind_020) = Rn_1(Ind_000) + dely * Rn_1(Ind_010)
    Rn_2(Ind_002) = Rn_1(Ind_000) + delz * Rn_1(Ind_001)
    Rn_2(Ind_110) = delx * Rn_1(Ind_010)
    Rn_2(Ind_101) = delx * Rn_1(Ind_001)
    Rn_2(Ind_011) = dely * Rn_1(Ind_001)
    Rn_3(Ind_000) = C(n-3) 
    Rn_3(Ind_100) = delx * Rn_2(Ind_000)
    Rn_3(Ind_010) = dely * Rn_2(Ind_000)
    Rn_3(Ind_001) = delz * Rn_2(Ind_000)
    Rn_3(Ind_200) = Rn_2(Ind_000) + delx * Rn_2(Ind_100)
    Rn_3(Ind_020) = Rn_2(Ind_000) + dely * Rn_2(Ind_010)
    Rn_3(Ind_002) = Rn_2(Ind_000) + delz * Rn_2(Ind_001)
    Rn_3(Ind_110) = delx * Rn_2(Ind_010)
    Rn_3(Ind_101) = delx * Rn_2(Ind_001)
    Rn_3(Ind_011) = dely * Rn_2(Ind_001)
    Rn_3(Ind_300) = 2.d0 * Rn_2(Ind_100) + delx * Rn_2(Ind_200)
    Rn_3(Ind_030) = 2.d0 * Rn_2(Ind_010) + dely * Rn_2(Ind_020)
    Rn_3(Ind_003) = 2.d0 * Rn_2(Ind_001) + delz * Rn_2(Ind_002)
    Rn_3(Ind_210) = dely * Rn_2(Ind_200)
    Rn_3(Ind_201) = delz * Rn_2(Ind_200)
    Rn_3(Ind_120) = delx * Rn_2(Ind_020)
    Rn_3(Ind_021) = delz * Rn_2(Ind_020)
    Rn_3(Ind_102) = delx * Rn_2(Ind_002)
    Rn_3(Ind_012) = dely * Rn_2(Ind_002)
    Rn_3(Ind_111) = delx * Rn_2(Ind_011)
    Rn_4(Ind_000) = C(n-4) 
    Rn_4(Ind_100) = delx * Rn_3(Ind_000)
    Rn_4(Ind_010) = dely * Rn_3(Ind_000)
    Rn_4(Ind_001) = delz * Rn_3(Ind_000)
    Rn_4(Ind_200) = Rn_3(Ind_000) + delx * Rn_3(Ind_100)
    Rn_4(Ind_020) = Rn_3(Ind_000) + dely * Rn_3(Ind_010)
    Rn_4(Ind_002) = Rn_3(Ind_000) + delz * Rn_3(Ind_001)
    Rn_4(Ind_110) = delx * Rn_3(Ind_010)
    Rn_4(Ind_101) = delx * Rn_3(Ind_001)
    Rn_4(Ind_011) = dely * Rn_3(Ind_001)
    Rn_4(Ind_300) = 2.d0 * Rn_3(Ind_100) + delx * Rn_3(Ind_200)
    Rn_4(Ind_030) = 2.d0 * Rn_3(Ind_010) + dely * Rn_3(Ind_020)
    Rn_4(Ind_003) = 2.d0 * Rn_3(Ind_001) + delz * Rn_3(Ind_002)
    Rn_4(Ind_210) = dely * Rn_3(Ind_200)
    Rn_4(Ind_201) = delz * Rn_3(Ind_200)
    Rn_4(Ind_120) = delx * Rn_3(Ind_020)
    Rn_4(Ind_021) = delz * Rn_3(Ind_020)
    Rn_4(Ind_102) = delx * Rn_3(Ind_002)
    Rn_4(Ind_012) = dely * Rn_3(Ind_002)
    Rn_4(Ind_111) = delx * Rn_3(Ind_011)
    Rn_4(Ind_400) = 3.d0 * Rn_3(Ind_200) + delx * Rn_3(Ind_300)
    Rn_4(Ind_040) = 3.d0 * Rn_3(Ind_020) + dely * Rn_3(Ind_030)
    Rn_4(Ind_004) = 3.d0 * Rn_3(Ind_002) + delz * Rn_3(Ind_003)
    Rn_4(Ind_310) = dely * Rn_3(Ind_300)
    Rn_4(Ind_301) = delz * Rn_3(Ind_300)
    Rn_4(Ind_130) = delx * Rn_3(Ind_030)
    Rn_4(Ind_031) = delz * Rn_3(Ind_030)
    Rn_4(Ind_103) = delx * Rn_3(Ind_003)
    Rn_4(Ind_013) = dely * Rn_3(Ind_003)
    Rn_4(Ind_220) = Rn_3(Ind_020) + delx * Rn_3(Ind_120)
    Rn_4(Ind_202) = Rn_3(Ind_002) + delx * Rn_3(Ind_102)
    Rn_4(Ind_022) = Rn_3(Ind_002) + dely * Rn_3(Ind_012)
    Rn_4(Ind_211) = dely * Rn_3(Ind_201)
    Rn_4(Ind_121) = delx * Rn_3(Ind_021)
    Rn_4(Ind_112) = delx * Rn_3(Ind_012)
    Rn_5(Ind_000) = C(n-5) 
    Rn_5(Ind_100) = delx * Rn_4(Ind_000)
    Rn_5(Ind_010) = dely * Rn_4(Ind_000)
    Rn_5(Ind_001) = delz * Rn_4(Ind_000)
    Rn_5(Ind_200) = Rn_4(Ind_000) + delx * Rn_4(Ind_100)
    Rn_5(Ind_020) = Rn_4(Ind_000) + dely * Rn_4(Ind_010)
    Rn_5(Ind_002) = Rn_4(Ind_000) + delz * Rn_4(Ind_001)
    Rn_5(Ind_110) = delx * Rn_4(Ind_010)
    Rn_5(Ind_101) = delx * Rn_4(Ind_001)
    Rn_5(Ind_011) = dely * Rn_4(Ind_001)
    Rn_5(Ind_300) = 2.d0 * Rn_4(Ind_100) + delx * Rn_4(Ind_200)
    Rn_5(Ind_030) = 2.d0 * Rn_4(Ind_010) + dely * Rn_4(Ind_020)
    Rn_5(Ind_003) = 2.d0 * Rn_4(Ind_001) + delz * Rn_4(Ind_002)
    Rn_5(Ind_210) = dely * Rn_4(Ind_200)
    Rn_5(Ind_201) = delz * Rn_4(Ind_200)
    Rn_5(Ind_120) = delx * Rn_4(Ind_020)
    Rn_5(Ind_021) = delz * Rn_4(Ind_020)
    Rn_5(Ind_102) = delx * Rn_4(Ind_002)
    Rn_5(Ind_012) = dely * Rn_4(Ind_002)
    Rn_5(Ind_111) = delx * Rn_4(Ind_011)
    Rn_5(Ind_400) = 3.d0 * Rn_4(Ind_200) + delx * Rn_4(Ind_300)
    Rn_5(Ind_040) = 3.d0 * Rn_4(Ind_020) + dely * Rn_4(Ind_030)
    Rn_5(Ind_004) = 3.d0 * Rn_4(Ind_002) + delz * Rn_4(Ind_003)
    Rn_5(Ind_310) = dely * Rn_4(Ind_300)
    Rn_5(Ind_301) = delz * Rn_4(Ind_300)
    Rn_5(Ind_130) = delx * Rn_4(Ind_030)
    Rn_5(Ind_031) = delz * Rn_4(Ind_030)
    Rn_5(Ind_103) = delx * Rn_4(Ind_003)
    Rn_5(Ind_013) = dely * Rn_4(Ind_003)
    Rn_5(Ind_220) = Rn_4(Ind_020) + delx * Rn_4(Ind_120)
    Rn_5(Ind_202) = Rn_4(Ind_002) + delx * Rn_4(Ind_102)
    Rn_5(Ind_022) = Rn_4(Ind_002) + dely * Rn_4(Ind_012)
    Rn_5(Ind_211) = dely * Rn_4(Ind_201)
    Rn_5(Ind_121) = delx * Rn_4(Ind_021)
    Rn_5(Ind_112) = delx * Rn_4(Ind_012)
    Rn_5(Ind_500) = 4.d0 * Rn_4(Ind_300) + delx * Rn_4(Ind_400)
    Rn_5(Ind_050) = 4.d0 * Rn_4(Ind_030) + dely * Rn_4(Ind_040)
    Rn_5(Ind_005) = 4.d0 * Rn_4(Ind_003) + delz * Rn_4(Ind_004)
    Rn_5(Ind_410) = dely * Rn_4(Ind_400)
    Rn_5(Ind_401) = delz * Rn_4(Ind_400)
    Rn_5(Ind_140) = delx * Rn_4(Ind_040)
    Rn_5(Ind_041) = delz * Rn_4(Ind_040)
    Rn_5(Ind_104) = delx * Rn_4(Ind_004)
    Rn_5(Ind_014) = dely * Rn_4(Ind_004)
    Rn_5(Ind_320) = Rn_4(Ind_300) + dely * Rn_4(Ind_310)
    Rn_5(Ind_302) = Rn_4(Ind_300) + delz * Rn_4(Ind_301)
    Rn_5(Ind_230) = Rn_4(Ind_030) + delx * Rn_4(Ind_130)
    Rn_5(Ind_032) = Rn_4(Ind_030) + delz * Rn_4(Ind_031)
    Rn_5(Ind_203) = Rn_4(Ind_003) + delx * Rn_4(Ind_103)
    Rn_5(Ind_023) = Rn_4(Ind_003) + dely * Rn_4(Ind_013)
    Rn_5(Ind_311) = dely * Rn_4(Ind_301)
    Rn_5(Ind_131) = delx * Rn_4(Ind_031)
    Rn_5(Ind_113) = delx * Rn_4(Ind_013)
    Rn_5(Ind_221) = delz * Rn_4(Ind_220)
    Rn_5(Ind_212) = dely * Rn_4(Ind_202)
    Rn_5(Ind_122) = delx * Rn_4(Ind_022)

    ! phi array (electrostatic potential at i due to j permanent moments
    ! and derivs of that esp wrt r_i)
    ! minus signs arise due to derivs of r_j - r_i wrt r_i

    do jj = 1, 10
       gmj(jj) = global_multipole(jj, j)
       gmi(jj) = global_multipole(jj, i)
    end do

    phi(Ind_000) = Rn_5(Ind_000) * gmj(Ind_000)+Rn_5(Ind_100) * gmj(Ind_100)+ &
                   Rn_5(Ind_010) * gmj(Ind_010)+Rn_5(Ind_001) * gmj(Ind_001)+ &
                   Rn_5(Ind_200) * gmj(Ind_200)+Rn_5(Ind_020) * gmj(Ind_020)+ &
                   Rn_5(Ind_002) * gmj(Ind_002)+Rn_5(Ind_110) * gmj(Ind_110)+ &
                   Rn_5(Ind_101) * gmj(Ind_101)+Rn_5(Ind_011) * gmj(Ind_011)
    phi(Ind_100) = -(Rn_5(Ind_100) * gmj(Ind_000)+Rn_5(Ind_200) * gmj(Ind_100)+&
                   Rn_5(Ind_110) * gmj(Ind_010)+Rn_5(Ind_101) * gmj(Ind_001)+ &
                   Rn_5(Ind_300) * gmj(Ind_200)+Rn_5(Ind_120) * gmj(Ind_020)+ &
                   Rn_5(Ind_102) * gmj(Ind_002)+Rn_5(Ind_210) * gmj(Ind_110)+ &
                   Rn_5(Ind_201) * gmj(Ind_101)+Rn_5(Ind_111) * gmj(Ind_011))
    phi(Ind_010) = -(Rn_5(Ind_010) * gmj(Ind_000)+Rn_5(Ind_110) * gmj(Ind_100)+&
                   Rn_5(Ind_020) * gmj(Ind_010)+Rn_5(Ind_011) * gmj(Ind_001)+ &
                   Rn_5(Ind_210) * gmj(Ind_200)+Rn_5(Ind_030) * gmj(Ind_020)+ &
                   Rn_5(Ind_012) * gmj(Ind_002)+Rn_5(Ind_120) * gmj(Ind_110)+ &
                   Rn_5(Ind_111) * gmj(Ind_101)+Rn_5(Ind_021) * gmj(Ind_011))
    phi(Ind_001) = -(Rn_5(Ind_001) * gmj(Ind_000)+Rn_5(Ind_101) * gmj(Ind_100)+&
                   Rn_5(Ind_011) * gmj(Ind_010)+Rn_5(Ind_002) * gmj(Ind_001)+ &
                   Rn_5(Ind_201) * gmj(Ind_200)+Rn_5(Ind_021) * gmj(Ind_020)+ &
                   Rn_5(Ind_003) * gmj(Ind_002)+Rn_5(Ind_111) * gmj(Ind_110)+ &
                   Rn_5(Ind_102) * gmj(Ind_101)+Rn_5(Ind_012) * gmj(Ind_011))
    phi(Ind_200) = Rn_5(Ind_200) * gmj(Ind_000)+Rn_5(Ind_300) * gmj(Ind_100)+ &
                   Rn_5(Ind_210) * gmj(Ind_010)+Rn_5(Ind_201) * gmj(Ind_001)+ &
                   Rn_5(Ind_400) * gmj(Ind_200)+Rn_5(Ind_220) * gmj(Ind_020)+ &
                   Rn_5(Ind_202) * gmj(Ind_002)+Rn_5(Ind_310) * gmj(Ind_110)+ &
                   Rn_5(Ind_301) * gmj(Ind_101)+Rn_5(Ind_211) * gmj(Ind_011)
    phi(Ind_020) = Rn_5(Ind_020) * gmj(Ind_000)+Rn_5(Ind_120) * gmj(Ind_100)+ &
                   Rn_5(Ind_030) * gmj(Ind_010)+Rn_5(Ind_021) * gmj(Ind_001)+ &
                   Rn_5(Ind_220) * gmj(Ind_200)+Rn_5(Ind_040) * gmj(Ind_020)+ &
                   Rn_5(Ind_022) * gmj(Ind_002)+Rn_5(Ind_130) * gmj(Ind_110)+ &
                   Rn_5(Ind_121) * gmj(Ind_101)+Rn_5(Ind_031) * gmj(Ind_011)
    phi(Ind_002) = Rn_5(Ind_002) * gmj(Ind_000)+Rn_5(Ind_102) * gmj(Ind_100)+ &
                   Rn_5(Ind_012) * gmj(Ind_010)+Rn_5(Ind_003) * gmj(Ind_001)+ &
                   Rn_5(Ind_202) * gmj(Ind_200)+Rn_5(Ind_022) * gmj(Ind_020)+ &
                   Rn_5(Ind_004) * gmj(Ind_002)+Rn_5(Ind_112) * gmj(Ind_110)+ &
                   Rn_5(Ind_103) * gmj(Ind_101)+Rn_5(Ind_013) * gmj(Ind_011)
    phi(Ind_110) = Rn_5(Ind_110) * gmj(Ind_000)+Rn_5(Ind_210) * gmj(Ind_100)+ &
                   Rn_5(Ind_120) * gmj(Ind_010)+Rn_5(Ind_111) * gmj(Ind_001)+ &
                   Rn_5(Ind_310) * gmj(Ind_200)+Rn_5(Ind_130) * gmj(Ind_020)+ &
                   Rn_5(Ind_112) * gmj(Ind_002)+Rn_5(Ind_220) * gmj(Ind_110)+ &
                   Rn_5(Ind_211) * gmj(Ind_101)+Rn_5(Ind_121) * gmj(Ind_011)
    phi(Ind_101) = Rn_5(Ind_101) * gmj(Ind_000)+Rn_5(Ind_201) * gmj(Ind_100)+ &
                   Rn_5(Ind_111) * gmj(Ind_010)+Rn_5(Ind_102) * gmj(Ind_001)+ &
                   Rn_5(Ind_301) * gmj(Ind_200)+Rn_5(Ind_121) * gmj(Ind_020)+ &
                   Rn_5(Ind_103) * gmj(Ind_002)+Rn_5(Ind_211) * gmj(Ind_110)+ &
                   Rn_5(Ind_202) * gmj(Ind_101)+Rn_5(Ind_112) * gmj(Ind_011)
    phi(Ind_011) = Rn_5(Ind_011) * gmj(Ind_000)+Rn_5(Ind_111) * gmj(Ind_100)+ &
                   Rn_5(Ind_021) * gmj(Ind_010)+Rn_5(Ind_012) * gmj(Ind_001)+ &
                   Rn_5(Ind_211) * gmj(Ind_200)+Rn_5(Ind_031) * gmj(Ind_020)+ &
                   Rn_5(Ind_013) * gmj(Ind_002)+Rn_5(Ind_121) * gmj(Ind_110)+ &
                   Rn_5(Ind_112) * gmj(Ind_101)+Rn_5(Ind_022) * gmj(Ind_011)
    phi(Ind_300) = -(Rn_5(Ind_300) * gmj(Ind_000)+Rn_5(Ind_400) * gmj(Ind_100)+&
                   Rn_5(Ind_310) * gmj(Ind_010)+Rn_5(Ind_301) * gmj(Ind_001)+ &
                   Rn_5(Ind_500) * gmj(Ind_200)+Rn_5(Ind_320) * gmj(Ind_020)+ &
                   Rn_5(Ind_302) * gmj(Ind_002)+Rn_5(Ind_410) * gmj(Ind_110)+ &
                   Rn_5(Ind_401) * gmj(Ind_101)+Rn_5(Ind_311) * gmj(Ind_011))
    phi(Ind_030) = -(Rn_5(Ind_030) * gmj(Ind_000)+Rn_5(Ind_130) * gmj(Ind_100)+&
                   Rn_5(Ind_040) * gmj(Ind_010)+Rn_5(Ind_031) * gmj(Ind_001)+ &
                   Rn_5(Ind_230) * gmj(Ind_200)+Rn_5(Ind_050) * gmj(Ind_020)+ &
                   Rn_5(Ind_032) * gmj(Ind_002)+Rn_5(Ind_140) * gmj(Ind_110)+ &
                   Rn_5(Ind_131) * gmj(Ind_101)+Rn_5(Ind_041) * gmj(Ind_011))
    phi(Ind_003) = -(Rn_5(Ind_003) * gmj(Ind_000)+Rn_5(Ind_103) * gmj(Ind_100)+&
                   Rn_5(Ind_013) * gmj(Ind_010)+Rn_5(Ind_004) * gmj(Ind_001)+ &
                   Rn_5(Ind_203) * gmj(Ind_200)+Rn_5(Ind_023) * gmj(Ind_020)+ &
                   Rn_5(Ind_005) * gmj(Ind_002)+Rn_5(Ind_113) * gmj(Ind_110)+ &
                   Rn_5(Ind_104) * gmj(Ind_101)+Rn_5(Ind_014) * gmj(Ind_011))
    phi(Ind_210) = -(Rn_5(Ind_210) * gmj(Ind_000)+Rn_5(Ind_310) * gmj(Ind_100)+&
                   Rn_5(Ind_220) * gmj(Ind_010)+Rn_5(Ind_211) * gmj(Ind_001)+ &
                   Rn_5(Ind_410) * gmj(Ind_200)+Rn_5(Ind_230) * gmj(Ind_020)+ &
                   Rn_5(Ind_212) * gmj(Ind_002)+Rn_5(Ind_320) * gmj(Ind_110)+ &
                   Rn_5(Ind_311) * gmj(Ind_101)+Rn_5(Ind_221) * gmj(Ind_011))
    phi(Ind_201) = -(Rn_5(Ind_201) * gmj(Ind_000)+Rn_5(Ind_301) * gmj(Ind_100)+&
                   Rn_5(Ind_211) * gmj(Ind_010)+Rn_5(Ind_202) * gmj(Ind_001)+ &
                   Rn_5(Ind_401) * gmj(Ind_200)+Rn_5(Ind_221) * gmj(Ind_020)+ &
                   Rn_5(Ind_203) * gmj(Ind_002)+Rn_5(Ind_311) * gmj(Ind_110)+ &
                   Rn_5(Ind_302) * gmj(Ind_101)+Rn_5(Ind_212) * gmj(Ind_011))
    phi(Ind_120) = -(Rn_5(Ind_120) * gmj(Ind_000)+Rn_5(Ind_220) * gmj(Ind_100)+&
                   Rn_5(Ind_130) * gmj(Ind_010)+Rn_5(Ind_121) * gmj(Ind_001)+ &
                   Rn_5(Ind_320) * gmj(Ind_200)+Rn_5(Ind_140) * gmj(Ind_020)+ &
                   Rn_5(Ind_122) * gmj(Ind_002)+Rn_5(Ind_230) * gmj(Ind_110)+ &
                   Rn_5(Ind_221) * gmj(Ind_101)+Rn_5(Ind_131) * gmj(Ind_011))
    phi(Ind_021) = -(Rn_5(Ind_021) * gmj(Ind_000)+Rn_5(Ind_121) * gmj(Ind_100)+&
                   Rn_5(Ind_031) * gmj(Ind_010)+Rn_5(Ind_022) * gmj(Ind_001)+ &
                   Rn_5(Ind_221) * gmj(Ind_200)+Rn_5(Ind_041) * gmj(Ind_020)+ &
                   Rn_5(Ind_023) * gmj(Ind_002)+Rn_5(Ind_131) * gmj(Ind_110)+ &
                   Rn_5(Ind_122) * gmj(Ind_101)+Rn_5(Ind_032) * gmj(Ind_011))
    phi(Ind_102) = -(Rn_5(Ind_102) * gmj(Ind_000)+Rn_5(Ind_202) * gmj(Ind_100)+&
                   Rn_5(Ind_112) * gmj(Ind_010)+Rn_5(Ind_103) * gmj(Ind_001)+ &
                   Rn_5(Ind_302) * gmj(Ind_200)+Rn_5(Ind_122) * gmj(Ind_020)+ &
                   Rn_5(Ind_104) * gmj(Ind_002)+Rn_5(Ind_212) * gmj(Ind_110)+ &
                   Rn_5(Ind_203) * gmj(Ind_101)+Rn_5(Ind_113) * gmj(Ind_011))
    phi(Ind_012) = -(Rn_5(Ind_012) * gmj(Ind_000)+Rn_5(Ind_112) * gmj(Ind_100)+&
                   Rn_5(Ind_022) * gmj(Ind_010)+Rn_5(Ind_013) * gmj(Ind_001)+ &
                   Rn_5(Ind_212) * gmj(Ind_200)+Rn_5(Ind_032) * gmj(Ind_020)+ &
                   Rn_5(Ind_014) * gmj(Ind_002)+Rn_5(Ind_122) * gmj(Ind_110)+ &
                   Rn_5(Ind_113) * gmj(Ind_101)+Rn_5(Ind_023) * gmj(Ind_011))
    phi(Ind_111) = -(Rn_5(Ind_111) * gmj(Ind_000)+Rn_5(Ind_211) * gmj(Ind_100)+&
                   Rn_5(Ind_121) * gmj(Ind_010)+Rn_5(Ind_112) * gmj(Ind_001)+ &
                   Rn_5(Ind_311) * gmj(Ind_200)+Rn_5(Ind_131) * gmj(Ind_020)+ &
                   Rn_5(Ind_113) * gmj(Ind_002)+Rn_5(Ind_221) * gmj(Ind_110)+ &
                   Rn_5(Ind_212) * gmj(Ind_101)+Rn_5(Ind_122) * gmj(Ind_011))

    e_pp = phi(Ind_000) * gmi(Ind_000) + phi(Ind_100) * gmi(Ind_100) + &
           phi(Ind_010) * gmi(Ind_010) + phi(Ind_001) * gmi(Ind_001) + &
           phi(Ind_200) * gmi(Ind_200) + phi(Ind_020) * gmi(Ind_020) + &
           phi(Ind_002) * gmi(Ind_002) + phi(Ind_110) * gmi(Ind_110) + &
           phi(Ind_101) * gmi(Ind_101) + phi(Ind_011) * gmi(Ind_011)

    ! gradient of e_pp wrt r_i

    g_pp(1) = phi(Ind_100) * gmi(Ind_000) + phi(Ind_200) * gmi(Ind_100) + &
              phi(Ind_110) * gmi(Ind_010) + phi(Ind_101) * gmi(Ind_001) + &
              phi(Ind_300) * gmi(Ind_200) + phi(Ind_120) * gmi(Ind_020) + &
              phi(Ind_102) * gmi(Ind_002) + phi(Ind_210) * gmi(Ind_110) + &
              phi(Ind_201) * gmi(Ind_101) + phi(Ind_111) * gmi(Ind_011)
    g_pp(2) = phi(Ind_010) * gmi(Ind_000) + phi(Ind_110) * gmi(Ind_100) + &
              phi(Ind_020) * gmi(Ind_010) + phi(Ind_011) * gmi(Ind_001) + &
              phi(Ind_210) * gmi(Ind_200) + phi(Ind_030) * gmi(Ind_020) + &
              phi(Ind_012) * gmi(Ind_002) + phi(Ind_120) * gmi(Ind_110) + &
              phi(Ind_111) * gmi(Ind_101) + phi(Ind_021) * gmi(Ind_011)
    g_pp(3) = phi(Ind_001) * gmi(Ind_000) + phi(Ind_101) * gmi(Ind_100) + &
              phi(Ind_011) * gmi(Ind_010) + phi(Ind_002) * gmi(Ind_001) + &
              phi(Ind_201) * gmi(Ind_200) + phi(Ind_021) * gmi(Ind_020) + &
              phi(Ind_003) * gmi(Ind_002) + phi(Ind_111) * gmi(Ind_110) + &
              phi(Ind_102) * gmi(Ind_101) + phi(Ind_012) * gmi(Ind_011)

    ! torque field at i due to permanent mpoles of j

    do jj = 1, 10
       torque_field(jj, i) = torque_field(jj, i) + phi(jj)
    end do

    ! next do field at j due to the permanent mpoles of i to get torque at j
    ! electrostatic potential at j due to permanent mpoles at i
    ! and derivatives of that with respect to r_j
    ! minus signs due to deriv of r_j-r_i wrt r_i entering into potential

    phi(Ind_000) = Rn_5(Ind_000) * gmi(Ind_000)-Rn_5(Ind_100) * gmi(Ind_100)- &
                   Rn_5(Ind_010) * gmi(Ind_010)-Rn_5(Ind_001) * gmi(Ind_001)+ &
                   Rn_5(Ind_200) * gmi(Ind_200)+Rn_5(Ind_020) * gmi(Ind_020)+ &
                   Rn_5(Ind_002) * gmi(Ind_002)+Rn_5(Ind_110) * gmi(Ind_110)+ &
                   Rn_5(Ind_101) * gmi(Ind_101)+Rn_5(Ind_011) * gmi(Ind_011)
    phi(Ind_100) = Rn_5(Ind_100) * gmi(Ind_000)-Rn_5(Ind_200) * gmi(Ind_100)- &
                   Rn_5(Ind_110) * gmi(Ind_010)-Rn_5(Ind_101) * gmi(Ind_001)+ &
                   Rn_5(Ind_300) * gmi(Ind_200)+Rn_5(Ind_120) * gmi(Ind_020)+ &
                   Rn_5(Ind_102) * gmi(Ind_002)+Rn_5(Ind_210) * gmi(Ind_110)+ &
                   Rn_5(Ind_201) * gmi(Ind_101)+Rn_5(Ind_111) * gmi(Ind_011)
    phi(Ind_010) = Rn_5(Ind_010) * gmi(Ind_000)-Rn_5(Ind_110) * gmi(Ind_100)- &
                   Rn_5(Ind_020) * gmi(Ind_010)-Rn_5(Ind_011) * gmi(Ind_001)+ &
                   Rn_5(Ind_210) * gmi(Ind_200)+Rn_5(Ind_030) * gmi(Ind_020)+ &
                   Rn_5(Ind_012) * gmi(Ind_002)+Rn_5(Ind_120) * gmi(Ind_110)+ &
                   Rn_5(Ind_111) * gmi(Ind_101)+Rn_5(Ind_021) * gmi(Ind_011)
    phi(Ind_001) = Rn_5(Ind_001) * gmi(Ind_000)-Rn_5(Ind_101) * gmi(Ind_100)- &
                   Rn_5(Ind_011) * gmi(Ind_010)-Rn_5(Ind_002) * gmi(Ind_001)+ &
                   Rn_5(Ind_201) * gmi(Ind_200)+Rn_5(Ind_021) * gmi(Ind_020)+ &
                   Rn_5(Ind_003) * gmi(Ind_002)+Rn_5(Ind_111) * gmi(Ind_110)+ &
                   Rn_5(Ind_102) * gmi(Ind_101)+Rn_5(Ind_012) * gmi(Ind_011)
    phi(Ind_200) = Rn_5(Ind_200) * gmi(Ind_000)-Rn_5(Ind_300) * gmi(Ind_100)- &
                   Rn_5(Ind_210) * gmi(Ind_010)-Rn_5(Ind_201) * gmi(Ind_001)+ &
                   Rn_5(Ind_400) * gmi(Ind_200)+Rn_5(Ind_220) * gmi(Ind_020)+ &
                   Rn_5(Ind_202) * gmi(Ind_002)+Rn_5(Ind_310) * gmi(Ind_110)+ &
                   Rn_5(Ind_301) * gmi(Ind_101)+Rn_5(Ind_211) * gmi(Ind_011)
    phi(Ind_020) = Rn_5(Ind_020) * gmi(Ind_000)-Rn_5(Ind_120) * gmi(Ind_100)- &
                   Rn_5(Ind_030) * gmi(Ind_010)-Rn_5(Ind_021) * gmi(Ind_001)+ &
                   Rn_5(Ind_220) * gmi(Ind_200)+Rn_5(Ind_040) * gmi(Ind_020)+ &
                   Rn_5(Ind_022) * gmi(Ind_002)+Rn_5(Ind_130) * gmi(Ind_110)+ &
                   Rn_5(Ind_121) * gmi(Ind_101)+Rn_5(Ind_031) * gmi(Ind_011)
    phi(Ind_002) = Rn_5(Ind_002) * gmi(Ind_000)-Rn_5(Ind_102) * gmi(Ind_100)- &
                   Rn_5(Ind_012) * gmi(Ind_010)-Rn_5(Ind_003) * gmi(Ind_001)+ &
                   Rn_5(Ind_202) * gmi(Ind_200)+Rn_5(Ind_022) * gmi(Ind_020)+ &
                   Rn_5(Ind_004) * gmi(Ind_002)+Rn_5(Ind_112) * gmi(Ind_110)+ &
                   Rn_5(Ind_103) * gmi(Ind_101)+Rn_5(Ind_013) * gmi(Ind_011)
    phi(Ind_110) = Rn_5(Ind_110) * gmi(Ind_000)-Rn_5(Ind_210) * gmi(Ind_100)- &
                   Rn_5(Ind_120) * gmi(Ind_010)-Rn_5(Ind_111) * gmi(Ind_001)+ &
                   Rn_5(Ind_310) * gmi(Ind_200)+Rn_5(Ind_130) * gmi(Ind_020)+ &
                   Rn_5(Ind_112) * gmi(Ind_002)+Rn_5(Ind_220) * gmi(Ind_110)+ &
                   Rn_5(Ind_211) * gmi(Ind_101)+Rn_5(Ind_121) * gmi(Ind_011)
    phi(Ind_101) = Rn_5(Ind_101) * gmi(Ind_000)-Rn_5(Ind_201) * gmi(Ind_100)- &
                   Rn_5(Ind_111) * gmi(Ind_010)-Rn_5(Ind_102) * gmi(Ind_001)+ &
                   Rn_5(Ind_301) * gmi(Ind_200)+Rn_5(Ind_121) * gmi(Ind_020)+ &
                   Rn_5(Ind_103) * gmi(Ind_002)+Rn_5(Ind_211) * gmi(Ind_110)+ &
                   Rn_5(Ind_202) * gmi(Ind_101)+Rn_5(Ind_112) * gmi(Ind_011)
    phi(Ind_011) = Rn_5(Ind_011) * gmi(Ind_000)-Rn_5(Ind_111) * gmi(Ind_100)- &
                   Rn_5(Ind_021) * gmi(Ind_010)-Rn_5(Ind_012) * gmi(Ind_001)+ &
                   Rn_5(Ind_211) * gmi(Ind_200)+Rn_5(Ind_031) * gmi(Ind_020)+ &
                   Rn_5(Ind_013) * gmi(Ind_002)+Rn_5(Ind_121) * gmi(Ind_110)+ &
                   Rn_5(Ind_112) * gmi(Ind_101)+Rn_5(Ind_022) * gmi(Ind_011)

    ! torque field at j due to permanent mpoles of i

    do jj = 1, 10
      torque_field(jj, j) = torque_field(jj, j) + phi(jj)
    end do

    e_ind = 0.d0
    g_ind(1) = 0.d0
    g_ind(2) = 0.d0
    g_ind(3) = 0.d0

    if (is_polarizable(i) .or. is_polarizable(j)) then
      do jj = 1, 3
        i_di(jj) = ind_dip_d(jj, i)
        i_dj(jj) = ind_dip_d(jj, j)
        i_pi(jj) = ind_dip_p(jj, i)
        i_pj(jj) = ind_dip_p(jj, j)
      end do

      ! first the direct dipoles interact with polar field
      ! recall A gives coulomb, BD damped and B the erfc * 1/r contributions

      do jj = 1, 4
        C(jj) = (B(jj) - A(jj)) + polar_weight(k) * (A(jj) - BD(jj))
      end do

      ! negate the odd order to get sign right

      C(1) = -C(1)
      C(3) = -C(3)

      ! get the interaction tensor 
       
      n = 4
      Rn(Ind_000) = C(n)
      Rn_1(Ind_000) = C(n-1)
      Rn_1(Ind_100) = delx * Rn(Ind_000)
      Rn_1(Ind_010) = dely * Rn(Ind_000)
      Rn_1(Ind_001) = delz * Rn(Ind_000)
      Rn_2(Ind_000) = C(n-2)
      Rn_2(Ind_100) = delx * Rn_1(Ind_000)
      Rn_2(Ind_010) = dely * Rn_1(Ind_000)
      Rn_2(Ind_001) = delz * Rn_1(Ind_000)
      Rn_2(Ind_200) = Rn_1(Ind_000) + delx * Rn_1(Ind_100)
      Rn_2(Ind_020) = Rn_1(Ind_000) + dely * Rn_1(Ind_010)
      Rn_2(Ind_002) = Rn_1(Ind_000) + delz * Rn_1(Ind_001)
      Rn_2(Ind_110) = delx * Rn_1(Ind_010)
      Rn_2(Ind_101) = delx * Rn_1(Ind_001)
      Rn_2(Ind_011) = dely * Rn_1(Ind_001)
      Rn_3(Ind_000) = C(n-3) 
      Rn_3(Ind_100) = delx * Rn_2(Ind_000)
      Rn_3(Ind_010) = dely * Rn_2(Ind_000)
      Rn_3(Ind_001) = delz * Rn_2(Ind_000)
      Rn_3(Ind_200) = Rn_2(Ind_000) + delx * Rn_2(Ind_100)
      Rn_3(Ind_020) = Rn_2(Ind_000) + dely * Rn_2(Ind_010)
      Rn_3(Ind_002) = Rn_2(Ind_000) + delz * Rn_2(Ind_001)
      Rn_3(Ind_110) = delx * Rn_2(Ind_010)
      Rn_3(Ind_101) = delx * Rn_2(Ind_001)
      Rn_3(Ind_011) = dely * Rn_2(Ind_001)
      Rn_3(Ind_300) = 2.d0 * Rn_2(Ind_100) + delx * Rn_2(Ind_200)
      Rn_3(Ind_030) = 2.d0 * Rn_2(Ind_010) + dely * Rn_2(Ind_020)
      Rn_3(Ind_003) = 2.d0 * Rn_2(Ind_001) + delz * Rn_2(Ind_002)
      Rn_3(Ind_210) = dely * Rn_2(Ind_200)
      Rn_3(Ind_201) = delz * Rn_2(Ind_200)
      Rn_3(Ind_120) = delx * Rn_2(Ind_020)
      Rn_3(Ind_021) = delz * Rn_2(Ind_020)
      Rn_3(Ind_102) = delx * Rn_2(Ind_002)
      Rn_3(Ind_012) = dely * Rn_2(Ind_002)
      Rn_3(Ind_111) = delx * Rn_2(Ind_011)
      !Rn_4(Ind_000) = C(n-4) NOT NEEDED
      Rn_4(Ind_100) = delx * Rn_3(Ind_000)
      Rn_4(Ind_010) = dely * Rn_3(Ind_000)
      Rn_4(Ind_001) = delz * Rn_3(Ind_000)
      Rn_4(Ind_200) = Rn_3(Ind_000) + delx * Rn_3(Ind_100)
      Rn_4(Ind_020) = Rn_3(Ind_000) + dely * Rn_3(Ind_010)
      Rn_4(Ind_002) = Rn_3(Ind_000) + delz * Rn_3(Ind_001)
      Rn_4(Ind_110) = delx * Rn_3(Ind_010)
      Rn_4(Ind_101) = delx * Rn_3(Ind_001)
      Rn_4(Ind_011) = dely * Rn_3(Ind_001)
      Rn_4(Ind_300) = 2.d0 * Rn_3(Ind_100) + delx * Rn_3(Ind_200)
      Rn_4(Ind_030) = 2.d0 * Rn_3(Ind_010) + dely * Rn_3(Ind_020)
      Rn_4(Ind_003) = 2.d0 * Rn_3(Ind_001) + delz * Rn_3(Ind_002)
      Rn_4(Ind_210) = dely * Rn_3(Ind_200)
      Rn_4(Ind_201) = delz * Rn_3(Ind_200)
      Rn_4(Ind_120) = delx * Rn_3(Ind_020)
      Rn_4(Ind_021) = delz * Rn_3(Ind_020)
      Rn_4(Ind_102) = delx * Rn_3(Ind_002)
      Rn_4(Ind_012) = dely * Rn_3(Ind_002)
      Rn_4(Ind_111) = delx * Rn_3(Ind_011)
      Rn_4(Ind_400) = 3.d0 * Rn_3(Ind_200) + delx * Rn_3(Ind_300)
      Rn_4(Ind_040) = 3.d0 * Rn_3(Ind_020) + dely * Rn_3(Ind_030)
      Rn_4(Ind_004) = 3.d0 * Rn_3(Ind_002) + delz * Rn_3(Ind_003)
      Rn_4(Ind_310) = dely * Rn_3(Ind_300)
      Rn_4(Ind_301) = delz * Rn_3(Ind_300)
      Rn_4(Ind_130) = delx * Rn_3(Ind_030)
      Rn_4(Ind_031) = delz * Rn_3(Ind_030)
      Rn_4(Ind_103) = delx * Rn_3(Ind_003)
      Rn_4(Ind_013) = dely * Rn_3(Ind_003)
      Rn_4(Ind_220) = Rn_3(Ind_020) + delx * Rn_3(Ind_120)
      Rn_4(Ind_202) = Rn_3(Ind_002) + delx * Rn_3(Ind_102)
      Rn_4(Ind_022) = Rn_3(Ind_002) + dely * Rn_3(Ind_012)
      Rn_4(Ind_211) = dely * Rn_3(Ind_201)
      Rn_4(Ind_121) = delx * Rn_3(Ind_021)
      Rn_4(Ind_112) = delx * Rn_3(Ind_012)

      if (is_polarizable(i)) then

        ! induced direct dipoles at i interact with permanent at j

        !phi(Ind_000)  NOT NEEDED
        !phi(Ind_000)= &
                !Rn_4(Ind_000) * gmj(Ind_000)+Rn_4(Ind_100) * gmj(Ind_100)+ &
                !Rn_4(Ind_010) * gmj(Ind_010)+Rn_4(Ind_001) * gmj(Ind_001)+ &
                !Rn_4(Ind_200) * gmj(Ind_200)+Rn_4(Ind_020) * gmj(Ind_020)+ &
                !Rn_4(Ind_002) * gmj(Ind_002)+Rn_4(Ind_110) * gmj(Ind_110)+ &
                !Rn_4(Ind_101) * gmj(Ind_101)+Rn_4(Ind_011) * gmj(Ind_011)
        phi(Ind_100)= &
              -(Rn_4(Ind_100) * gmj(Ind_000)+Rn_4(Ind_200) * gmj(Ind_100)+ &
                Rn_4(Ind_110) * gmj(Ind_010)+Rn_4(Ind_101) * gmj(Ind_001)+ &
                Rn_4(Ind_300) * gmj(Ind_200)+Rn_4(Ind_120) * gmj(Ind_020)+ &
                Rn_4(Ind_102) * gmj(Ind_002)+Rn_4(Ind_210) * gmj(Ind_110)+ &
                Rn_4(Ind_201) * gmj(Ind_101)+Rn_4(Ind_111) * gmj(Ind_011))
        phi(Ind_010)= &
              -(Rn_4(Ind_010) * gmj(Ind_000)+Rn_4(Ind_110) * gmj(Ind_100)+ &
                Rn_4(Ind_020) * gmj(Ind_010)+Rn_4(Ind_011) * gmj(Ind_001)+ &
                Rn_4(Ind_210) * gmj(Ind_200)+Rn_4(Ind_030) * gmj(Ind_020)+ &
                Rn_4(Ind_012) * gmj(Ind_002)+Rn_4(Ind_120) * gmj(Ind_110)+ &
                Rn_4(Ind_111) * gmj(Ind_101)+Rn_4(Ind_021) * gmj(Ind_011))
        phi(Ind_001)= &
              -(Rn_4(Ind_001) * gmj(Ind_000)+Rn_4(Ind_101) * gmj(Ind_100)+ &
                Rn_4(Ind_011) * gmj(Ind_010)+Rn_4(Ind_002) * gmj(Ind_001)+ &
                Rn_4(Ind_201) * gmj(Ind_200)+Rn_4(Ind_021) * gmj(Ind_020)+ &
                Rn_4(Ind_003) * gmj(Ind_002)+Rn_4(Ind_111) * gmj(Ind_110)+ &
                Rn_4(Ind_102) * gmj(Ind_101)+Rn_4(Ind_012) * gmj(Ind_011))
        phi(Ind_200)=  &
                Rn_4(Ind_200) * gmj(Ind_000)+Rn_4(Ind_300) * gmj(Ind_100)+ &
                Rn_4(Ind_210) * gmj(Ind_010)+Rn_4(Ind_201) * gmj(Ind_001)+ &
                Rn_4(Ind_400) * gmj(Ind_200)+Rn_4(Ind_220) * gmj(Ind_020)+ &
                Rn_4(Ind_202) * gmj(Ind_002)+Rn_4(Ind_310) * gmj(Ind_110)+ &
                Rn_4(Ind_301) * gmj(Ind_101)+Rn_4(Ind_211) * gmj(Ind_011)
        phi(Ind_020)= &
                Rn_4(Ind_020) * gmj(Ind_000)+Rn_4(Ind_120) * gmj(Ind_100)+ &
                Rn_4(Ind_030) * gmj(Ind_010)+Rn_4(Ind_021) * gmj(Ind_001)+ &
                Rn_4(Ind_220) * gmj(Ind_200)+Rn_4(Ind_040) * gmj(Ind_020)+ &
                Rn_4(Ind_022) * gmj(Ind_002)+Rn_4(Ind_130) * gmj(Ind_110)+ &
                Rn_4(Ind_121) * gmj(Ind_101)+Rn_4(Ind_031) * gmj(Ind_011)
        phi(Ind_002)= &
                Rn_4(Ind_002) * gmj(Ind_000)+Rn_4(Ind_102) * gmj(Ind_100)+ &
                Rn_4(Ind_012) * gmj(Ind_010)+Rn_4(Ind_003) * gmj(Ind_001)+ &
                Rn_4(Ind_202) * gmj(Ind_200)+Rn_4(Ind_022) * gmj(Ind_020)+ &
                Rn_4(Ind_004) * gmj(Ind_002)+Rn_4(Ind_112) * gmj(Ind_110)+ &
                Rn_4(Ind_103) * gmj(Ind_101)+Rn_4(Ind_013) * gmj(Ind_011)
        phi(Ind_110)= &
                Rn_4(Ind_110) * gmj(Ind_000)+Rn_4(Ind_210) * gmj(Ind_100)+ &
                Rn_4(Ind_120) * gmj(Ind_010)+Rn_4(Ind_111) * gmj(Ind_001)+ &
                Rn_4(Ind_310) * gmj(Ind_200)+Rn_4(Ind_130) * gmj(Ind_020)+ &
                Rn_4(Ind_112) * gmj(Ind_002)+Rn_4(Ind_220) * gmj(Ind_110)+ &
                Rn_4(Ind_211) * gmj(Ind_101)+Rn_4(Ind_121) * gmj(Ind_011)
        phi(Ind_101)= &
                Rn_4(Ind_101) * gmj(Ind_000)+Rn_4(Ind_201) * gmj(Ind_100)+ &
                Rn_4(Ind_111) * gmj(Ind_010)+Rn_4(Ind_102) * gmj(Ind_001)+ &
                Rn_4(Ind_301) * gmj(Ind_200)+Rn_4(Ind_121) * gmj(Ind_020)+ &
                Rn_4(Ind_103) * gmj(Ind_002)+Rn_4(Ind_211) * gmj(Ind_110)+ &
                Rn_4(Ind_202) * gmj(Ind_101)+Rn_4(Ind_112) * gmj(Ind_011)
        phi(Ind_011)=  &
                Rn_4(Ind_011) * gmj(Ind_000)+Rn_4(Ind_111) * gmj(Ind_100)+ &
                Rn_4(Ind_021) * gmj(Ind_010)+Rn_4(Ind_012) * gmj(Ind_001)+ &
                Rn_4(Ind_211) * gmj(Ind_200)+Rn_4(Ind_031) * gmj(Ind_020)+ &
                Rn_4(Ind_013) * gmj(Ind_002)+Rn_4(Ind_121) * gmj(Ind_110)+ &
                Rn_4(Ind_112) * gmj(Ind_101)+Rn_4(Ind_022) * gmj(Ind_011)
        e_ind = e_ind + 0.5d0 *  &
                (phi(Ind_100) * i_di(1) + phi(Ind_010) * i_di(2) + &
                 phi(Ind_001) * i_di(3))

        g_ind(1) = g_ind(1) + 0.5d0 *  &
                   (phi(Ind_200) * i_di(1) + phi(Ind_110) * i_di(2) + &
                    phi(Ind_101) * i_di(3))
        g_ind(2) = g_ind(2) + 0.5d0 *  &
                   (phi(Ind_110) * i_di(1) + phi(Ind_020) * i_di(2) + &
                    phi(Ind_011) * i_di(3))
        g_ind(3) = g_ind(3) + 0.5d0 *  &
                   (phi(Ind_101) * i_di(1) + phi(Ind_011) * i_di(2) + &
                    phi(Ind_002) * i_di(3))

        ! calculate the potential at j due to direct induced at i
        ! and derivatives of that with respect to r_j
        ! minus signs due to deriv of r_j-r_i wrt r_i entering into potential

        phi(Ind_000) = -(Rn_4(Ind_100) * i_di(1)+Rn_4(Ind_010) * i_di(2) + &
                         Rn_4(Ind_001) * i_di(3))
        phi(Ind_100) = -(Rn_4(Ind_200) * i_di(1)+Rn_4(Ind_110) * i_di(2) + &
                         Rn_4(Ind_101) * i_di(3))
        phi(Ind_010) = -(Rn_4(Ind_110) * i_di(1)+Rn_4(Ind_020) * i_di(2) + &
                         Rn_4(Ind_011) * i_di(3))
        phi(Ind_001) = -(Rn_4(Ind_101) * i_di(1)+Rn_4(Ind_011) * i_di(2) + &
                         Rn_4(Ind_002) * i_di(3))
        phi(Ind_200) = -(Rn_4(Ind_300) * i_di(1)+Rn_4(Ind_210) * i_di(2) + &
                         Rn_4(Ind_201) * i_di(3))
        phi(Ind_020) = -(Rn_4(Ind_120) * i_di(1)+Rn_4(Ind_030) * i_di(2) + &
                         Rn_4(Ind_021) * i_di(3))
        phi(Ind_002) = -(Rn_4(Ind_102) * i_di(1)+Rn_4(Ind_012) * i_di(2) + &
                         Rn_4(Ind_003) * i_di(3))
        phi(Ind_110) = -(Rn_4(Ind_210) * i_di(1)+Rn_4(Ind_120) * i_di(2) + &
                         Rn_4(Ind_111) * i_di(3))
        phi(Ind_101) = -(Rn_4(Ind_201) * i_di(1)+Rn_4(Ind_111) * i_di(2) + &
                         Rn_4(Ind_102) * i_di(3))
        phi(Ind_011) = -(Rn_4(Ind_111) * i_di(1)+Rn_4(Ind_021) * i_di(2) + &
                         Rn_4(Ind_012) * i_di(3))

        ! torque field at j due to direct induced mpoles of i
        ! note factor of 1/2

        do jj = 1, 10
          torque_field(jj, j) = torque_field(jj, j) + 0.5d0 * phi(jj)
        end do

      end if !is_polarizable(i)

      if (is_polarizable(j)) then

        ! induced direct dipoles at j interact with permanent at i
        ! phi array (electrostatic potential at i due to j induced moments
        ! and derivs of that esp wrt r_i)
        ! first that due to ind_dip_d for energy contribution
        ! minus signs arise due to derivs of r_j - r_i wrt r_i

        phi(Ind_000) = Rn_4(Ind_100) * i_dj(1)+Rn_4(Ind_010) * i_dj(2) + &
                       Rn_4(Ind_001) * i_dj(3)
        phi(Ind_100) = -(Rn_4(Ind_200) * i_dj(1)+Rn_4(Ind_110) * i_dj(2) + &
                         Rn_4(Ind_101) * i_dj(3))
        phi(Ind_010) = -(Rn_4(Ind_110) * i_dj(1)+Rn_4(Ind_020) * i_dj(2) + &
                         Rn_4(Ind_011) * i_dj(3))
        phi(Ind_001) = -(Rn_4(Ind_101) * i_dj(1)+Rn_4(Ind_011) * i_dj(2) + &
                         Rn_4(Ind_002) * i_dj(3))
        phi(Ind_200) = Rn_4(Ind_300) * i_dj(1)+Rn_4(Ind_210) * i_dj(2) + &
                       Rn_4(Ind_201) * i_dj(3)
        phi(Ind_020) = Rn_4(Ind_120) * i_dj(1)+Rn_4(Ind_030) * i_dj(2) + &
                       Rn_4(Ind_021) * i_dj(3)
        phi(Ind_002) = Rn_4(Ind_102) * i_dj(1)+Rn_4(Ind_012) * i_dj(2) + &
                       Rn_4(Ind_003) * i_dj(3)
        phi(Ind_110) = Rn_4(Ind_210) * i_dj(1)+Rn_4(Ind_120) * i_dj(2) + &
                       Rn_4(Ind_111) * i_dj(3)
        phi(Ind_101) = Rn_4(Ind_201) * i_dj(1)+Rn_4(Ind_111) * i_dj(2) + &
                       Rn_4(Ind_102) * i_dj(3)
        phi(Ind_011) = Rn_4(Ind_111) * i_dj(1)+Rn_4(Ind_021) * i_dj(2) + &
                       Rn_4(Ind_012) * i_dj(3)
        phi(Ind_300) = -(Rn_4(Ind_400) * i_dj(1)+Rn_4(Ind_310) * i_dj(2) + &
                         Rn_4(Ind_301) * i_dj(3))
        phi(Ind_030) = -(Rn_4(Ind_130) * i_dj(1)+Rn_4(Ind_040) * i_dj(2) + &
                         Rn_4(Ind_031) * i_dj(3))
        phi(Ind_003) = -(Rn_4(Ind_103) * i_dj(1)+Rn_4(Ind_013) * i_dj(2) + &
                         Rn_4(Ind_004) * i_dj(3))
        phi(Ind_210) = -(Rn_4(Ind_310) * i_dj(1)+Rn_4(Ind_220) * i_dj(2) + &
                         Rn_4(Ind_211) * i_dj(3))
        phi(Ind_201) = -(Rn_4(Ind_301) * i_dj(1)+Rn_4(Ind_211) * i_dj(2) + &
                         Rn_4(Ind_202) * i_dj(3))
        phi(Ind_120) = -(Rn_4(Ind_220) * i_dj(1)+Rn_4(Ind_130) * i_dj(2) + &
                         Rn_4(Ind_121) * i_dj(3))
        phi(Ind_021) = -(Rn_4(Ind_121) * i_dj(1)+Rn_4(Ind_031) * i_dj(2) + &
                         Rn_4(Ind_022) * i_dj(3))
        phi(Ind_102) = -(Rn_4(Ind_202) * i_dj(1)+Rn_4(Ind_112) * i_dj(2) + &
                         Rn_4(Ind_103) * i_dj(3))
        phi(Ind_012) = -(Rn_4(Ind_112) * i_dj(1)+Rn_4(Ind_022) * i_dj(2) + &
                         Rn_4(Ind_013) * i_dj(3))
        phi(Ind_111) = -(Rn_4(Ind_211) * i_dj(1)+Rn_4(Ind_121) * i_dj(2) + &
                         Rn_4(Ind_112) * i_dj(3))
                         
        e_add = 0.5d0 * &
                (phi(Ind_000) * gmi(Ind_000)+phi(Ind_100) * gmi(Ind_100)+ &
                 phi(Ind_010) * gmi(Ind_010)+phi(Ind_001) * gmi(Ind_001)+ &
                 phi(Ind_200) * gmi(Ind_200)+phi(Ind_020) * gmi(Ind_020)+ &
                 phi(Ind_002) * gmi(Ind_002)+phi(Ind_110) * gmi(Ind_110)+ & 
                 phi(Ind_101) * gmi(Ind_101)+phi(Ind_011) * gmi(Ind_011))
        
        e_ind = e_ind + e_add
        
        g_ind(1) = g_ind(1) + 0.5d0 * &
                (phi(Ind_100) * gmi(Ind_000)+phi(Ind_200) * gmi(Ind_100)+ &
                 phi(Ind_110) * gmi(Ind_010)+phi(Ind_101) * gmi(Ind_001)+ &
                 phi(Ind_300) * gmi(Ind_200)+phi(Ind_120) * gmi(Ind_020)+ &
                 phi(Ind_102) * gmi(Ind_002)+phi(Ind_210) * gmi(Ind_110)+ &
                 phi(Ind_201) * gmi(Ind_101)+phi(Ind_111) * gmi(Ind_011))
        g_ind(2) = g_ind(2) + 0.5d0 * &
                (phi(Ind_010) * gmi(Ind_000)+phi(Ind_110) * gmi(Ind_100)+ &
                 phi(Ind_020) * gmi(Ind_010)+phi(Ind_011) * gmi(Ind_001)+ &
                 phi(Ind_210) * gmi(Ind_200)+phi(Ind_030) * gmi(Ind_020)+ &
                 phi(Ind_012) * gmi(Ind_002)+phi(Ind_120) * gmi(Ind_110)+ &
                 phi(Ind_111) * gmi(Ind_101)+phi(Ind_021) * gmi(Ind_011))
        g_ind(3) = g_ind(3) + 0.5d0 * &
                (phi(Ind_001) * gmi(Ind_000)+phi(Ind_101) * gmi(Ind_100)+ &
                 phi(Ind_011) * gmi(Ind_010)+phi(Ind_002) * gmi(Ind_001)+ &
                 phi(Ind_201) * gmi(Ind_200)+phi(Ind_021) * gmi(Ind_020)+ &
                 phi(Ind_003) * gmi(Ind_002)+phi(Ind_111) * gmi(Ind_110)+ &
                 phi(Ind_102) * gmi(Ind_101)+phi(Ind_012) * gmi(Ind_011))

        ! torque field at i due to direct dipoles of j
        ! note factor of 1/2

        do jj = 1, 10
          torque_field(jj, i) = torque_field(jj, i) + 0.5d0 * phi(jj)
        end do

      end if ! is_polarizable(j)) then

      ! next the polar dipoles interact with direct field in grad
      ! recall A gives coulomb, BD damped and B the erfc * 1/r contributions

      do jj = 1, 4
        C(jj) = (B(jj) - A(jj)) + direct_weight(k) * (A(jj) - BD(jj))
      end do
      
      ! negate the odd order to get sign right

      C(1) = -C(1)
      C(3) = -C(3)

      ! get the interaction tensor 

      n = 4
      Rn(Ind_000) = C(n)
      Rn_1(Ind_000) = C(n-1)
      Rn_1(Ind_100) = delx * Rn(Ind_000)
      Rn_1(Ind_010) = dely * Rn(Ind_000)
      Rn_1(Ind_001) = delz * Rn(Ind_000)
      Rn_2(Ind_000) = C(n-2)
      Rn_2(Ind_100) = delx * Rn_1(Ind_000)
      Rn_2(Ind_010) = dely * Rn_1(Ind_000)
      Rn_2(Ind_001) = delz * Rn_1(Ind_000)
      Rn_2(Ind_200) = Rn_1(Ind_000) + delx * Rn_1(Ind_100)
      Rn_2(Ind_020) = Rn_1(Ind_000) + dely * Rn_1(Ind_010)
      Rn_2(Ind_002) = Rn_1(Ind_000) + delz * Rn_1(Ind_001)
      Rn_2(Ind_110) = delx * Rn_1(Ind_010)
      Rn_2(Ind_101) = delx * Rn_1(Ind_001)
      Rn_2(Ind_011) = dely * Rn_1(Ind_001)
      Rn_3(Ind_000) = C(n-3) 
      Rn_3(Ind_100) = delx * Rn_2(Ind_000)
      Rn_3(Ind_010) = dely * Rn_2(Ind_000)
      Rn_3(Ind_001) = delz * Rn_2(Ind_000)
      Rn_3(Ind_200) = Rn_2(Ind_000) + delx * Rn_2(Ind_100)
      Rn_3(Ind_020) = Rn_2(Ind_000) + dely * Rn_2(Ind_010)
      Rn_3(Ind_002) = Rn_2(Ind_000) + delz * Rn_2(Ind_001)
      Rn_3(Ind_110) = delx * Rn_2(Ind_010)
      Rn_3(Ind_101) = delx * Rn_2(Ind_001)
      Rn_3(Ind_011) = dely * Rn_2(Ind_001)
      Rn_3(Ind_300) = 2.d0 * Rn_2(Ind_100) + delx * Rn_2(Ind_200)
      Rn_3(Ind_030) = 2.d0 * Rn_2(Ind_010) + dely * Rn_2(Ind_020)
      Rn_3(Ind_003) = 2.d0 * Rn_2(Ind_001) + delz * Rn_2(Ind_002)
      Rn_3(Ind_210) = dely * Rn_2(Ind_200)
      Rn_3(Ind_201) = delz * Rn_2(Ind_200)
      Rn_3(Ind_120) = delx * Rn_2(Ind_020)
      Rn_3(Ind_021) = delz * Rn_2(Ind_020)
      Rn_3(Ind_102) = delx * Rn_2(Ind_002)
      Rn_3(Ind_012) = dely * Rn_2(Ind_002)
      Rn_3(Ind_111) = delx * Rn_2(Ind_011)
      !Rn_4(Ind_000) = C(n-4) NOT NEEDED
      Rn_4(Ind_100) = delx * Rn_3(Ind_000)
      Rn_4(Ind_010) = dely * Rn_3(Ind_000)
      Rn_4(Ind_001) = delz * Rn_3(Ind_000)
      Rn_4(Ind_200) = Rn_3(Ind_000) + delx * Rn_3(Ind_100)
      Rn_4(Ind_020) = Rn_3(Ind_000) + dely * Rn_3(Ind_010)
      Rn_4(Ind_002) = Rn_3(Ind_000) + delz * Rn_3(Ind_001)
      Rn_4(Ind_110) = delx * Rn_3(Ind_010)
      Rn_4(Ind_101) = delx * Rn_3(Ind_001)
      Rn_4(Ind_011) = dely * Rn_3(Ind_001)
      Rn_4(Ind_300) = 2.d0 * Rn_3(Ind_100) + delx * Rn_3(Ind_200)
      Rn_4(Ind_030) = 2.d0 * Rn_3(Ind_010) + dely * Rn_3(Ind_020)
      Rn_4(Ind_003) = 2.d0 * Rn_3(Ind_001) + delz * Rn_3(Ind_002)
      Rn_4(Ind_210) = dely * Rn_3(Ind_200)
      Rn_4(Ind_201) = delz * Rn_3(Ind_200)
      Rn_4(Ind_120) = delx * Rn_3(Ind_020)
      Rn_4(Ind_021) = delz * Rn_3(Ind_020)
      Rn_4(Ind_102) = delx * Rn_3(Ind_002)
      Rn_4(Ind_012) = dely * Rn_3(Ind_002)
      Rn_4(Ind_111) = delx * Rn_3(Ind_011)
      Rn_4(Ind_400) = 3.d0 * Rn_3(Ind_200) + delx * Rn_3(Ind_300)
      Rn_4(Ind_040) = 3.d0 * Rn_3(Ind_020) + dely * Rn_3(Ind_030)
      Rn_4(Ind_004) = 3.d0 * Rn_3(Ind_002) + delz * Rn_3(Ind_003)
      Rn_4(Ind_310) = dely * Rn_3(Ind_300)
      Rn_4(Ind_301) = delz * Rn_3(Ind_300)
      Rn_4(Ind_130) = delx * Rn_3(Ind_030)
      Rn_4(Ind_031) = delz * Rn_3(Ind_030)
      Rn_4(Ind_103) = delx * Rn_3(Ind_003)
      Rn_4(Ind_013) = dely * Rn_3(Ind_003)
      Rn_4(Ind_220) = Rn_3(Ind_020) + delx * Rn_3(Ind_120)
      Rn_4(Ind_202) = Rn_3(Ind_002) + delx * Rn_3(Ind_102)
      Rn_4(Ind_022) = Rn_3(Ind_002) + dely * Rn_3(Ind_012)
      Rn_4(Ind_211) = dely * Rn_3(Ind_201)
      Rn_4(Ind_121) = delx * Rn_3(Ind_021)
      Rn_4(Ind_112) = delx * Rn_3(Ind_012)

      if (is_polarizable(i)) then
        ! induced direct dipoles at i interact with permanent at j
        !phi(Ind_000)  NOT NEEDED
        !phi(Ind_000)= &
                !Rn_4(Ind_000) * gmj(Ind_000)+Rn_4(Ind_100) * gmj(Ind_100)+ &
                !Rn_4(Ind_010) * gmj(Ind_010)+Rn_4(Ind_001) * gmj(Ind_001)+ &
                !Rn_4(Ind_200) * gmj(Ind_200)+Rn_4(Ind_020) * gmj(Ind_020)+ &
                !Rn_4(Ind_002) * gmj(Ind_002)+Rn_4(Ind_110) * gmj(Ind_110)+ &
                !Rn_4(Ind_101) * gmj(Ind_101)+Rn_4(Ind_011) * gmj(Ind_011)
        phi(Ind_100)= &
              -(Rn_4(Ind_100) * gmj(Ind_000)+Rn_4(Ind_200) * gmj(Ind_100)+ &
                Rn_4(Ind_110) * gmj(Ind_010)+Rn_4(Ind_101) * gmj(Ind_001)+ &
                Rn_4(Ind_300) * gmj(Ind_200)+Rn_4(Ind_120) * gmj(Ind_020)+ &
                Rn_4(Ind_102) * gmj(Ind_002)+Rn_4(Ind_210) * gmj(Ind_110)+ &
                Rn_4(Ind_201) * gmj(Ind_101)+Rn_4(Ind_111) * gmj(Ind_011))
        phi(Ind_010)= &
              -(Rn_4(Ind_010) * gmj(Ind_000)+Rn_4(Ind_110) * gmj(Ind_100)+ &
                Rn_4(Ind_020) * gmj(Ind_010)+Rn_4(Ind_011) * gmj(Ind_001)+ &
                Rn_4(Ind_210) * gmj(Ind_200)+Rn_4(Ind_030) * gmj(Ind_020)+ &
                Rn_4(Ind_012) * gmj(Ind_002)+Rn_4(Ind_120) * gmj(Ind_110)+ &
                Rn_4(Ind_111) * gmj(Ind_101)+Rn_4(Ind_021) * gmj(Ind_011))
        phi(Ind_001)= &
              -(Rn_4(Ind_001) * gmj(Ind_000)+Rn_4(Ind_101) * gmj(Ind_100)+ &
                Rn_4(Ind_011) * gmj(Ind_010)+Rn_4(Ind_002) * gmj(Ind_001)+ &
                Rn_4(Ind_201) * gmj(Ind_200)+Rn_4(Ind_021) * gmj(Ind_020)+ &
                Rn_4(Ind_003) * gmj(Ind_002)+Rn_4(Ind_111) * gmj(Ind_110)+ &
                Rn_4(Ind_102) * gmj(Ind_101)+Rn_4(Ind_012) * gmj(Ind_011))
        phi(Ind_200)=  &
                Rn_4(Ind_200) * gmj(Ind_000)+Rn_4(Ind_300) * gmj(Ind_100)+ &
                Rn_4(Ind_210) * gmj(Ind_010)+Rn_4(Ind_201) * gmj(Ind_001)+ &
                Rn_4(Ind_400) * gmj(Ind_200)+Rn_4(Ind_220) * gmj(Ind_020)+ &
                Rn_4(Ind_202) * gmj(Ind_002)+Rn_4(Ind_310) * gmj(Ind_110)+ &
                Rn_4(Ind_301) * gmj(Ind_101)+Rn_4(Ind_211) * gmj(Ind_011)
        phi(Ind_020)= &
                Rn_4(Ind_020) * gmj(Ind_000)+Rn_4(Ind_120) * gmj(Ind_100)+ &
                Rn_4(Ind_030) * gmj(Ind_010)+Rn_4(Ind_021) * gmj(Ind_001)+ &
                Rn_4(Ind_220) * gmj(Ind_200)+Rn_4(Ind_040) * gmj(Ind_020)+ &
                Rn_4(Ind_022) * gmj(Ind_002)+Rn_4(Ind_130) * gmj(Ind_110)+ &
                Rn_4(Ind_121) * gmj(Ind_101)+Rn_4(Ind_031) * gmj(Ind_011)
        phi(Ind_002)= &
                Rn_4(Ind_002) * gmj(Ind_000)+Rn_4(Ind_102) * gmj(Ind_100)+ &
                Rn_4(Ind_012) * gmj(Ind_010)+Rn_4(Ind_003) * gmj(Ind_001)+ &
                Rn_4(Ind_202) * gmj(Ind_200)+Rn_4(Ind_022) * gmj(Ind_020)+ &
                Rn_4(Ind_004) * gmj(Ind_002)+Rn_4(Ind_112) * gmj(Ind_110)+ &
                Rn_4(Ind_103) * gmj(Ind_101)+Rn_4(Ind_013) * gmj(Ind_011)
        phi(Ind_110)= &
                Rn_4(Ind_110) * gmj(Ind_000)+Rn_4(Ind_210) * gmj(Ind_100)+ &
                Rn_4(Ind_120) * gmj(Ind_010)+Rn_4(Ind_111) * gmj(Ind_001)+ &
                Rn_4(Ind_310) * gmj(Ind_200)+Rn_4(Ind_130) * gmj(Ind_020)+ &
                Rn_4(Ind_112) * gmj(Ind_002)+Rn_4(Ind_220) * gmj(Ind_110)+ &
                Rn_4(Ind_211) * gmj(Ind_101)+Rn_4(Ind_121) * gmj(Ind_011)
        phi(Ind_101)= &
                Rn_4(Ind_101) * gmj(Ind_000)+Rn_4(Ind_201) * gmj(Ind_100)+ &
                Rn_4(Ind_111) * gmj(Ind_010)+Rn_4(Ind_102) * gmj(Ind_001)+ &
                Rn_4(Ind_301) * gmj(Ind_200)+Rn_4(Ind_121) * gmj(Ind_020)+ &
                Rn_4(Ind_103) * gmj(Ind_002)+Rn_4(Ind_211) * gmj(Ind_110)+ &
                Rn_4(Ind_202) * gmj(Ind_101)+Rn_4(Ind_112) * gmj(Ind_011)
        phi(Ind_011)=  &
                Rn_4(Ind_011) * gmj(Ind_000)+Rn_4(Ind_111) * gmj(Ind_100)+ &
                Rn_4(Ind_021) * gmj(Ind_010)+Rn_4(Ind_012) * gmj(Ind_001)+ &
                Rn_4(Ind_211) * gmj(Ind_200)+Rn_4(Ind_031) * gmj(Ind_020)+ &
                Rn_4(Ind_013) * gmj(Ind_002)+Rn_4(Ind_121) * gmj(Ind_110)+ &
                Rn_4(Ind_112) * gmj(Ind_101)+Rn_4(Ind_022) * gmj(Ind_011)

        g_ind(1) = g_ind(1) + 0.5d0 *  &
                (phi(Ind_200) * i_pi(1) + phi(Ind_110) * i_pi(2) + &
                   phi(Ind_101) * i_pi(3))
        g_ind(2) = g_ind(2) + 0.5d0 *  &
                (phi(Ind_110) * i_pi(1) + phi(Ind_020) * i_pi(2) + &
                   phi(Ind_011) * i_pi(3))
        g_ind(3) = g_ind(3) + 0.5d0 *  &
                (phi(Ind_101) * i_pi(1) + phi(Ind_011) * i_pi(2) + &
                   phi(Ind_002) * i_pi(3))

        ! calculate the potential at j due to polar induced at i
        ! and derivatives of that with respect to r_j
        ! minus signs due to deriv of r_j-r_i wrt r_i entering into potential

        phi(Ind_000) = -(Rn_4(Ind_100) * i_pi(1)+Rn_4(Ind_010) * i_pi(2) + &
                         Rn_4(Ind_001) * i_pi(3))
        phi(Ind_100) = -(Rn_4(Ind_200) * i_pi(1)+Rn_4(Ind_110) * i_pi(2) + &
                         Rn_4(Ind_101) * i_pi(3))
        phi(Ind_010) = -(Rn_4(Ind_110) * i_pi(1)+Rn_4(Ind_020) * i_pi(2) + &
                         Rn_4(Ind_011) * i_pi(3))
        phi(Ind_001) = -(Rn_4(Ind_101) * i_pi(1)+Rn_4(Ind_011) * i_pi(2) + &
                         Rn_4(Ind_002) * i_pi(3))
        phi(Ind_200) = -(Rn_4(Ind_300) * i_pi(1)+Rn_4(Ind_210) * i_pi(2) + &
                         Rn_4(Ind_201) * i_pi(3))
        phi(Ind_020) = -(Rn_4(Ind_120) * i_pi(1)+Rn_4(Ind_030) * i_pi(2) + &
                         Rn_4(Ind_021) * i_pi(3))
        phi(Ind_002) = -(Rn_4(Ind_102) * i_pi(1)+Rn_4(Ind_012) * i_pi(2) + &
                         Rn_4(Ind_003) * i_pi(3))
        phi(Ind_110) = -(Rn_4(Ind_210) * i_pi(1)+Rn_4(Ind_120) * i_pi(2) + &
                         Rn_4(Ind_111) * i_pi(3))
        phi(Ind_101) = -(Rn_4(Ind_201) * i_pi(1)+Rn_4(Ind_111) * i_pi(2) + &
                         Rn_4(Ind_102) * i_pi(3))
        phi(Ind_011) = -(Rn_4(Ind_111) * i_pi(1)+Rn_4(Ind_021) * i_pi(2) + &
                         Rn_4(Ind_012) * i_pi(3))

        ! torque field at j due to direct induced mpoles of i
        ! note factor of 1/2

        do jj = 1, 10
           torque_field(jj, j) = torque_field(jj, j) + 0.5d0 * phi(jj)
        end do

      end if !is_polarizable(i)

      if (is_polarizable(j)) then

        ! induced polar dipoles at j interact with permanent at i
        ! phi array (electrostatic potential at i due to j induced moments
        ! and derivs of that esp wrt r_i)
        ! minus signs arise due to derivs of r_j - r_i wrt r_i

        phi(Ind_000) = Rn_4(Ind_100) * i_pj(1)+Rn_4(Ind_010) * i_pj(2) + &
                       Rn_4(Ind_001) * i_pj(3)
        phi(Ind_100) = -(Rn_4(Ind_200) * i_pj(1)+Rn_4(Ind_110) * i_pj(2) + &
                         Rn_4(Ind_101) * i_pj(3))
        phi(Ind_010) = -(Rn_4(Ind_110) * i_pj(1)+Rn_4(Ind_020) * i_pj(2) + &
                         Rn_4(Ind_011) * i_pj(3))
        phi(Ind_001) = -(Rn_4(Ind_101) * i_pj(1)+Rn_4(Ind_011) * i_pj(2) + &
                         Rn_4(Ind_002) * i_pj(3))
        phi(Ind_200) = Rn_4(Ind_300) * i_pj(1)+Rn_4(Ind_210) * i_pj(2) + &
                       Rn_4(Ind_201) * i_pj(3)
        phi(Ind_020) = Rn_4(Ind_120) * i_pj(1)+Rn_4(Ind_030) * i_pj(2) + &
                       Rn_4(Ind_021) * i_pj(3)
        phi(Ind_002) = Rn_4(Ind_102) * i_pj(1)+Rn_4(Ind_012) * i_pj(2) + &
                       Rn_4(Ind_003) * i_pj(3)
        phi(Ind_110) = Rn_4(Ind_210) * i_pj(1)+Rn_4(Ind_120) * i_pj(2) + &
                       Rn_4(Ind_111) * i_pj(3)
        phi(Ind_101) = Rn_4(Ind_201) * i_pj(1)+Rn_4(Ind_111) * i_pj(2) + &
                       Rn_4(Ind_102) * i_pj(3)
        phi(Ind_011) = Rn_4(Ind_111) * i_pj(1)+Rn_4(Ind_021) * i_pj(2) + &
                       Rn_4(Ind_012) * i_pj(3)
        phi(Ind_300) = -(Rn_4(Ind_400) * i_pj(1)+Rn_4(Ind_310) * i_pj(2) + &
                         Rn_4(Ind_301) * i_pj(3))
        phi(Ind_030) = -(Rn_4(Ind_130) * i_pj(1)+Rn_4(Ind_040) * i_pj(2) + &
                         Rn_4(Ind_031) * i_pj(3))
        phi(Ind_003) = -(Rn_4(Ind_103) * i_pj(1)+Rn_4(Ind_013) * i_pj(2) + &
                         Rn_4(Ind_004) * i_pj(3))
        phi(Ind_210) = -(Rn_4(Ind_310) * i_pj(1)+Rn_4(Ind_220) * i_pj(2) + &
                         Rn_4(Ind_211) * i_pj(3))
        phi(Ind_201) = -(Rn_4(Ind_301) * i_pj(1)+Rn_4(Ind_211) * i_pj(2) + &
                         Rn_4(Ind_202) * i_pj(3))
        phi(Ind_120) = -(Rn_4(Ind_220) * i_pj(1)+Rn_4(Ind_130) * i_pj(2) + &
                         Rn_4(Ind_121) * i_pj(3))
        phi(Ind_021) = -(Rn_4(Ind_121) * i_pj(1)+Rn_4(Ind_031) * i_pj(2) + &
                         Rn_4(Ind_022) * i_pj(3))
        phi(Ind_102) = -(Rn_4(Ind_202) * i_pj(1)+Rn_4(Ind_112) * i_pj(2) + &
                         Rn_4(Ind_103) * i_pj(3))
        phi(Ind_012) = -(Rn_4(Ind_112) * i_pj(1)+Rn_4(Ind_022) * i_pj(2) + &
                         Rn_4(Ind_013) * i_pj(3))
        phi(Ind_111) = -(Rn_4(Ind_211) * i_pj(1)+Rn_4(Ind_121) * i_pj(2) + &
                         Rn_4(Ind_112) * i_pj(3))

        g_ind(1) = g_ind(1) + 0.5d0 *  &
                (phi(Ind_100) * gmi(Ind_000)+phi(Ind_200) * gmi(Ind_100)+ &
                 phi(Ind_110) * gmi(Ind_010)+phi(Ind_101) * gmi(Ind_001)+ &
                 phi(Ind_300) * gmi(Ind_200)+phi(Ind_120) * gmi(Ind_020)+ &
                 phi(Ind_102) * gmi(Ind_002)+phi(Ind_210) * gmi(Ind_110)+ &
                 phi(Ind_201) * gmi(Ind_101)+phi(Ind_111) * gmi(Ind_011))
        g_ind(2) = g_ind(2) + 0.5d0 *  &
                (phi(Ind_010) * gmi(Ind_000)+phi(Ind_110) * gmi(Ind_100)+ &
                 phi(Ind_020) * gmi(Ind_010)+phi(Ind_011) * gmi(Ind_001)+ &
                 phi(Ind_210) * gmi(Ind_200)+phi(Ind_030) * gmi(Ind_020)+ &
                 phi(Ind_012) * gmi(Ind_002)+phi(Ind_120) * gmi(Ind_110)+ &
                 phi(Ind_111) * gmi(Ind_101)+phi(Ind_021) * gmi(Ind_011))
        g_ind(3) = g_ind(3) + 0.5d0 *  &
                (phi(Ind_001) * gmi(Ind_000)+phi(Ind_101) * gmi(Ind_100)+ &
                 phi(Ind_011) * gmi(Ind_010)+phi(Ind_002) * gmi(Ind_001)+ &
                 phi(Ind_201) * gmi(Ind_200)+phi(Ind_021) * gmi(Ind_020)+ &
                 phi(Ind_003) * gmi(Ind_002)+phi(Ind_111) * gmi(Ind_110)+ &
                 phi(Ind_102) * gmi(Ind_101)+phi(Ind_012) * gmi(Ind_011))

        ! torque field at i due to polar dipoles of j
        ! note factor of 1/2

        do jj = 1, 10
          torque_field(jj, i) = torque_field(jj, i) + 0.5d0 * phi(jj)
        end do

      end if ! is_polarizable(j)) then

      ! finally the induced-induced interactions

      if (is_polarizable(i) .and. is_polarizable(j)) then

        ! recall A gives coulomb, BD damped and B the erfc * 1/r contributions

        do jj = 1, 3
          C(jj) = (B(jj) - A(jj)) + mutual_weight(k) * (A(jj) - BD(jj))
        end do

        ! negate the odd order to get sign right

        C(1) = -C(1)
        C(3) = -C(3)

        ! get the interaction tensor 

        n = 3
        Rn(Ind_000) = C(n)
        Rn_1(Ind_000) = C(n-1)
        Rn_1(Ind_100) = delx * Rn(Ind_000)
        Rn_1(Ind_010) = dely * Rn(Ind_000)
        Rn_1(Ind_001) = delz * Rn(Ind_000)
        Rn_2(Ind_000) = C(n-2)
        Rn_2(Ind_100) = delx * Rn_1(Ind_000)
        Rn_2(Ind_010) = dely * Rn_1(Ind_000)
        Rn_2(Ind_001) = delz * Rn_1(Ind_000)
        Rn_2(Ind_200) = Rn_1(Ind_000) + delx * Rn_1(Ind_100)
        Rn_2(Ind_020) = Rn_1(Ind_000) + dely * Rn_1(Ind_010)
        Rn_2(Ind_002) = Rn_1(Ind_000) + delz * Rn_1(Ind_001)
        Rn_2(Ind_110) = delx * Rn_1(Ind_010)
        Rn_2(Ind_101) = delx * Rn_1(Ind_001)
        Rn_2(Ind_011) = dely * Rn_1(Ind_001)
        !Rn_3(Ind_000) = C(n-3) NOT NEEDED
        Rn_3(Ind_100) = delx * Rn_2(Ind_000)
        Rn_3(Ind_010) = dely * Rn_2(Ind_000)
        Rn_3(Ind_001) = delz * Rn_2(Ind_000)
        Rn_3(Ind_200) = Rn_2(Ind_000) + delx * Rn_2(Ind_100)
        Rn_3(Ind_020) = Rn_2(Ind_000) + dely * Rn_2(Ind_010)
        Rn_3(Ind_002) = Rn_2(Ind_000) + delz * Rn_2(Ind_001)
        Rn_3(Ind_110) = delx * Rn_2(Ind_010)
        Rn_3(Ind_101) = delx * Rn_2(Ind_001)
        Rn_3(Ind_011) = dely * Rn_2(Ind_001)
        Rn_3(Ind_300) = 2.d0 * Rn_2(Ind_100) + delx * Rn_2(Ind_200)
        Rn_3(Ind_030) = 2.d0 * Rn_2(Ind_010) + dely * Rn_2(Ind_020)
        Rn_3(Ind_003) = 2.d0 * Rn_2(Ind_001) + delz * Rn_2(Ind_002)
        Rn_3(Ind_210) = dely * Rn_2(Ind_200)
        Rn_3(Ind_201) = delz * Rn_2(Ind_200)
        Rn_3(Ind_120) = delx * Rn_2(Ind_020)
        Rn_3(Ind_021) = delz * Rn_2(Ind_020)
        Rn_3(Ind_102) = delx * Rn_2(Ind_002)
        Rn_3(Ind_012) = dely * Rn_2(Ind_002)
        Rn_3(Ind_111) = delx * Rn_2(Ind_011)
        ! phi array (electrostatic pot at i due to j induced direct
        ! and derivs of that esp wrt r_i)
        phi(Ind_200) = Rn_3(Ind_300) * i_dj(1)+Rn_3(Ind_210) * i_dj(2)+ &
                       Rn_3(Ind_201) * i_dj(3)
        phi(Ind_020) = Rn_3(Ind_120) * i_dj(1)+Rn_3(Ind_030) * i_dj(2)+ &
                       Rn_3(Ind_021) * i_dj(3)
        phi(Ind_002) = Rn_3(Ind_102) * i_dj(1)+Rn_3(Ind_012) * i_dj(2)+ &
                       Rn_3(Ind_003) * i_dj(3)
        phi(Ind_110) = Rn_3(Ind_210) * i_dj(1)+Rn_3(Ind_120) * i_dj(2)+ &
                       Rn_3(Ind_111) * i_dj(3)
        phi(Ind_101) = Rn_3(Ind_201) * i_dj(1)+Rn_3(Ind_111) * i_dj(2)+ &
                       Rn_3(Ind_102) * i_dj(3)
        phi(Ind_011) = Rn_3(Ind_111) * i_dj(1)+Rn_3(Ind_021) * i_dj(2)+ &
                       Rn_3(Ind_012) * i_dj(3)
        g_ind(1) = g_ind(1) + 0.5d0* &
                      (phi(Ind_200) * i_pi(1)+phi(Ind_110) * i_pi(2) + &
                       phi(Ind_101) * i_pi(3))
        g_ind(2) = g_ind(2) + 0.5d0* &
                      (phi(Ind_110) * i_pi(1)+phi(Ind_020) * i_pi(2) + &
                       phi(Ind_011) * i_pi(3))
        g_ind(3) = g_ind(3) + 0.5d0* &
                      (phi(Ind_101) * i_pi(1)+phi(Ind_011) * i_pi(2) + &
                       phi(Ind_002) * i_pi(3))

        ! phi array (electrostatic pot at i due to j induced polar
        ! and derivs of that esp wrt r_i)

        phi(Ind_200) = Rn_3(Ind_300) * i_pj(1)+Rn_3(Ind_210) * i_pj(2)+ &
                       Rn_3(Ind_201) * i_pj(3)
        phi(Ind_020) = Rn_3(Ind_120) * i_pj(1)+Rn_3(Ind_030) * i_pj(2)+ &
                       Rn_3(Ind_021) * i_pj(3)
        phi(Ind_002) = Rn_3(Ind_102) * i_pj(1)+Rn_3(Ind_012) * i_pj(2)+ &
                       Rn_3(Ind_003) * i_pj(3)
        phi(Ind_110) = Rn_3(Ind_210) * i_pj(1)+Rn_3(Ind_120) * i_pj(2)+ &
                       Rn_3(Ind_111) * i_pj(3)
        phi(Ind_101) = Rn_3(Ind_201) * i_pj(1)+Rn_3(Ind_111) * i_pj(2)+ &
                       Rn_3(Ind_102) * i_pj(3)
        phi(Ind_011) = Rn_3(Ind_111) * i_pj(1)+Rn_3(Ind_021) * i_pj(2)+ &
                       Rn_3(Ind_012) * i_pj(3)
        g_ind(1) = g_ind(1) + 0.5d0* &
                      (phi(Ind_200) * i_di(1)+phi(Ind_110) * i_di(2) + &
                       phi(Ind_101) * i_di(3))
        g_ind(2) = g_ind(2) + 0.5d0* &
                      (phi(Ind_110) * i_di(1)+phi(Ind_020) * i_di(2) + &
                       phi(Ind_011) * i_di(3))
        g_ind(3) = g_ind(3) + 0.5d0* &
                      (phi(Ind_101) * i_di(1)+phi(Ind_011) * i_di(2) + &
                       phi(Ind_002) * i_di(3))

      end if !is_polarizable(i) .and. is_polarizable(j) 

    end if !(is_polarizable(i) .or. is_polarizable(j)) then

    ! frc is negative gradient
    
    ene_perm = ene_perm + coulomb_const_kcal_per_mole * e_pp
    ene_ind  = ene_ind + coulomb_const_kcal_per_mole * e_ind
    frc(1, j) = frc(1, j) + coulomb_const_kcal_per_mole * (g_pp(1) + g_ind(1))
    frc(2, j) = frc(2, j) + coulomb_const_kcal_per_mole * (g_pp(2) + g_ind(2))
    frc(3, j) = frc(3, j) + coulomb_const_kcal_per_mole * (g_pp(3) + g_ind(3))
    frc(1, i) = frc(1, i) - coulomb_const_kcal_per_mole * (g_pp(1) + g_ind(1))
    frc(2, i) = frc(2, i) - coulomb_const_kcal_per_mole * (g_pp(2) + g_ind(2))
    frc(3, i) = frc(3, i) - coulomb_const_kcal_per_mole * (g_pp(3) + g_ind(3))
    vxx = vxx - coulomb_const_kcal_per_mole * delx * (g_pp(1) + g_ind(1))
    vxy = vxy - coulomb_const_kcal_per_mole * delx * (g_pp(2) + g_ind(2))
    vxz = vxz - coulomb_const_kcal_per_mole * delx * (g_pp(3) + g_ind(3))
    vyx = vyx - coulomb_const_kcal_per_mole * dely * (g_pp(1) + g_ind(1))
    vyy = vyy - coulomb_const_kcal_per_mole * dely * (g_pp(2) + g_ind(2))
    vyz = vyz - coulomb_const_kcal_per_mole * dely * (g_pp(3) + g_ind(3))
    vzx = vzx - coulomb_const_kcal_per_mole * delz * (g_pp(1) + g_ind(1))
    vzy = vzy - coulomb_const_kcal_per_mole * delz * (g_pp(2) + g_ind(2))
    vzz = vzz - coulomb_const_kcal_per_mole * delz * (g_pp(3) + g_ind(3))

  end do !n_adj = 1, num_adjust_list

  virial(1, 1) = virial(1, 1) + vxx
  virial(1, 2) = virial(1, 2) + 0.5d0 * (vxy + vyx)
  virial(1, 3) = virial(1, 3) + 0.5d0 * (vxz + vzx)
  virial(2, 1) = virial(2, 1) + 0.5d0 * (vxy + vyx)
  virial(2, 2) = virial(2, 2) + vyy
  virial(2, 3) = virial(2, 3) + 0.5d0 * (vyz + vzy)
  virial(3, 1) = virial(3, 1) + 0.5d0 * (vxz + vzx)
  virial(3, 2) = virial(3, 2) + 0.5d0 * (vyz + vzy)
  virial(3, 3) = virial(3, 3) + vzz

  return

end subroutine am_adjust_calc_ene_frc

subroutine set_amoeba_adjust_list(list_new,n)
	
	implicit none
	integer, intent(in)         :: list_new(3,n)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	
	nold = size(adjust_list,2)
	num_adjust_list = n
	
    if( nold .ne. n ) then
        if( nold .gt. 0) then 
            deallocate(adjust_list)
            deallocate(adjust_weight_idxs)
        endif
        allocate( adjust_list(3,num_adjust_list), &
                  adjust_weight_idxs(2 * num_adjust_list), &
                  stat = alloc_failed )
        if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( num_adjust_list .eq. 0) return 
	    
	adjust_list(:,1:n) = list_new(:,1:n) 
	
end subroutine set_amoeba_adjust_list

subroutine set_adjust_weights( vdw_weight_new, mpole_weight_new, direct_weight_new, polar_weight_new, mutual_weight_new )
	
	implicit none
	double precision, intent(in) :: vdw_weight_new(9)
	double precision, intent(in) :: mpole_weight_new(9)
	double precision, intent(in) :: direct_weight_new(9)
	double precision, intent(in) :: polar_weight_new(9)
	double precision, intent(in) :: mutual_weight_new(9)
		    
	vdw_weight(1:9)    = vdw_weight_new(1:9)
	mpole_weight(1:9)  = mpole_weight_new(1:9)
	direct_weight(1:9) = direct_weight_new(1:9)
	polar_weight(1:9)  = polar_weight_new(1:9)
	mutual_weight(1:9) = mutual_weight_new(1:9) 
	
end subroutine set_adjust_weights

subroutine set_valid_bit(ival)

    use amoeba_flags_mod
    
    implicit none
	integer, intent(in)         :: ival
    
    if(ival .eq. 1) then
        do_amoeba_adjust_flag = ibset(do_amoeba_adjust_flag, valid_bit)
    else
        do_amoeba_adjust_flag = ibclr(do_amoeba_adjust_flag, valid_bit)
    endif
        
end subroutine set_valid_bit

end module amoeba_adjust_mod
