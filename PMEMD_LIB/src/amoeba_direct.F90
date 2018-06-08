#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_direct_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_direct_mod

  implicit none

  private

  integer, save ::  do_amoeba_direct_flag

#include "amoeba_mpole_index.i" 

  ! Process-local storage:

  integer, save                         :: num_pairs_in_ee_cut = 0
  integer, save                         :: size_dipole_dipole_list = -1
  integer, save                         :: num_tensor

  double precision, save, allocatable   :: dipole_dipole_tensor(:,:)
  integer, save, allocatable            :: dipole_dipole_list(:,:)

  double precision, parameter           :: safety = 1.5d0   !  factor to predict maximum number of dipole-dipole interactions increased from 1.25
  double precision, parameter           :: checklist = 0.9d0

  public        am_direct_zero_flag
  public        am_direct_set_user_bit
  public        am_direct_permfield
  public        am_direct_dip_dip_field
  public        am_direct_ene_frc

contains

!*******************************************************************************!
! Subroutine:  am_direct_zero_flag
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_direct_zero_flag
  implicit none
  do_amoeba_direct_flag = 0
  return
end subroutine am_direct_zero_flag

!*******************************************************************************!
! Subroutine:  am_direct_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_direct_set_user_bit(do_this)

  use amoeba_flags_mod
  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in) :: do_this

! Set the valid bit---this part always since no parmread needed.

  do_amoeba_direct_flag = ibset(do_amoeba_direct_flag, valid_bit)

  if (do_this .eq. 1) then        ! do in all cases
    do_amoeba_direct_flag = ibset(do_amoeba_direct_flag, user_induce_bit)
    do_amoeba_direct_flag = ibset(do_amoeba_direct_flag, user_postinduce_bit)
  else if (do_this .eq. 2) then   ! do the induction, not the post-induction
    do_amoeba_direct_flag = ibset(do_amoeba_direct_flag, user_induce_bit)
    do_amoeba_direct_flag = ibclr(do_amoeba_direct_flag, user_postinduce_bit)
  else if (do_this .eq. 3) then   ! do the post-induction, not the induction
    do_amoeba_direct_flag = ibclr(do_amoeba_direct_flag, user_induce_bit)
    do_amoeba_direct_flag = ibset(do_amoeba_direct_flag, user_postinduce_bit)
  else if (do_this .eq. 0) then 
    do_amoeba_direct_flag = ibclr(do_amoeba_direct_flag, user_induce_bit)
    do_amoeba_direct_flag = ibclr(do_amoeba_direct_flag, user_postinduce_bit)
  else
    write(error_msg,*)'am_direct_set_user_bit: bad value of user do_this = ',do_this
    call mol_mech_error 
  end if

  return

end subroutine am_direct_set_user_bit

!*******************************************************************************!
! Subroutine:  am_direct_permfield
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_direct_permfield(cart_dipole_field, ipairs, img, tranvec, &
                               img_atm_map)

  use amoeba_multipoles_mod, only : global_multipole
  use amoeba_induced_mod, only : is_polarizable, screen_polar, damp_polar_strength, damp_polar_sensitivity, &
                         damp_polar_rad, polarizability_corr
  use mdin_amoeba_dat_mod, only : ee_dsum_cut, ee_damped_cut, thole_expon_coeff
  use file_io_dat_mod
  use mdin_ewald_dat_mod
  use pme_direct_mod
  use img_mod
  use parallel_dat_mod
  use timers_mod

  implicit none

! Formal arguments:

  double precision, intent(in out)      :: cart_dipole_field(3, *)
  integer                               :: ipairs(*)
  type(img_rec), intent(in)             :: img(*)
  double precision, intent(in)          :: tranvec(1:3, 0:17)
  integer, intent(in)                   :: img_atm_map(*)

! Local variables:

  double precision                      :: x_i, y_i, z_i
  double precision                      :: x_tran(1:3, 0:17)
  double precision                      :: ee_dsum_cut2
  double precision                      :: ee_damped_cut2
  integer                               :: ier
  integer                               :: atm_i
  integer                               :: i
  integer                               :: i_loc
  integer                               :: j_loc
  integer                               :: ipairs_idx
  integer                               :: img_i
  integer                               :: pair_cnt
  integer                               :: alloc_failed
#ifdef DIRFRC_COMTRANS
  ! flag - 1 if translation not needed
  integer                               :: common_tran
#endif /* DIRFRC_COMTRANS */

  ! check if dipole_dipole_list big enough

  if (num_pairs_in_ee_cut .gt. checklist * size_dipole_dipole_list) then

    if (allocated(dipole_dipole_tensor)) deallocate(dipole_dipole_tensor)
    if (allocated(dipole_dipole_list)) deallocate(dipole_dipole_list)

    call am_direct_count_num_ee_pairs(ipairs, img, tranvec, &
                                      ee_dsum_cut, num_pairs_in_ee_cut)

    size_dipole_dipole_list = safety * num_pairs_in_ee_cut

    allocate(dipole_dipole_tensor(6, size_dipole_dipole_list), &
             dipole_dipole_list(2, size_dipole_dipole_list), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

  end if

  ee_dsum_cut2 = ee_dsum_cut * ee_dsum_cut
  ee_damped_cut2 = ee_damped_cut * ee_damped_cut

  ipairs_idx = 1
  num_tensor = 0
  
  do img_i = my_img_lo, my_img_hi

#ifdef DIRFRC_COMTRANS
    ! Common translation (ie. no translation) flag is packed at
    ! the front of each sublist followed by the count(s) of sublist
    ! image pair entries.

    common_tran = ipairs(ipairs_idx)
    ipairs_idx = ipairs_idx + 1
#endif /* DIRFRC_COMTRANS */
    
    ! Electrostatic evaluation-only count followed by
    ! full evaluation count packed at the front of each pair sublist.

    pair_cnt = ipairs(ipairs_idx) + ipairs(ipairs_idx + 1)
    ipairs_idx = ipairs_idx + 2

    if (pair_cnt .gt. 0) then

      x_i = img(img_i)%x
      y_i = img(img_i)%y
      z_i = img(img_i)%z

#ifdef DIRFRC_COMTRANS
      if (common_tran .eq. 0) then
#endif /* DIRFRC_COMTRANS */
        ! We need all the translation vectors:
        do i = 0, 17
          x_tran(1, i) = tranvec(1, i) - x_i
          x_tran(2, i) = tranvec(2, i) - y_i
          x_tran(3, i) = tranvec(3, i) - z_i
        end do
#ifdef DIRFRC_COMTRANS
      else
        ! Just put the x,y,z values in the middle cell
        x_tran(1, 13) = - x_i 
        x_tran(2, 13) = - y_i
        x_tran(3, 13) = - z_i
      end if
#endif /* DIRFRC_COMTRANS */

      atm_i = img_atm_map(img_i)

      call am_direct_permfield_i(atm_i, img, ipairs(ipairs_idx), x_tran, &
                                 pair_cnt, ew_coeff, eedtbdns, gbl_eed_cub, &
                                 ee_dsum_cut2, ee_damped_cut2, &
                                 dipole_dipole_tensor, dipole_dipole_list, &
                                 size_dipole_dipole_list, num_tensor, &
                                 global_multipole, cart_dipole_field, &
                                 img_atm_map)  
    
      ipairs_idx = ipairs_idx + pair_cnt

    end if
  end do
  
  call update_pme_time(dir_frc_sum_timer)

  return

contains

!*******************************************************************************!
! Subroutine:  am_direct_permfield_i
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_direct_permfield_i(atm_i, img, ipairs_sublst, x_tran, &
                                 pair_cnt, ewaldcof, eedtbdns, eed_cub, &
                                 ee_dsum_cut2, ee_damped_cut2, &
                                 dipole_dipole_tensor, dipole_dipole_list, &
                                 size_dipole_dipole_list, num_tensor, &
                                 global_multipole, gradphi, &
                                 img_atm_map)

  use amoeba_flags_mod
  use parallel_dat_mod
  use img_mod
  use amoeba_induced_mod, only : ndamp

  implicit none

! Formal arguments:

  integer, intent(in)                   :: atm_i
  type(img_rec), intent(in)             :: img(*)
  integer                               :: ipairs_sublst(*)
  double precision, intent(in)          :: x_tran(1:3, 0:17)
  integer, intent(in)                   :: pair_cnt
  double precision, intent(in)          :: ewaldcof
  double precision, intent(in)          :: eedtbdns
  double precision, intent(in)          :: eed_cub(4, *)
  double precision, intent(in)          :: ee_dsum_cut2
  double precision, intent(in)          :: ee_damped_cut2
  double precision, intent(out)         :: dipole_dipole_tensor(6, *)
  integer, intent(out)                  :: dipole_dipole_list(2, *)
  integer, intent(in)                   :: size_dipole_dipole_list
  integer, intent(in out)               :: num_tensor
  double precision, intent(in)          :: global_multipole(10, *)
  double precision, intent(in out)      :: gradphi(3, *)
  integer, intent(in)                   :: img_atm_map(*)

! Local variables:

  integer, parameter                    :: vs = 128      ! vector size
  integer, parameter                    :: mask27 = Z"07FFFFFF"
  double precision, parameter           :: third = 1.d0/3.d0

  integer                               :: atm_j, img_j, enc_img
  integer                               :: itran
  integer                               :: sublst_head
  integer                               :: vec_idx, vec_cnt, vec_max
  integer                               :: n, ind
  double precision                      :: delx, dely, delz
  double precision                      :: delr, delr2, delr2inv
  double precision                      :: x, dx, switch, d_switch_dx
  double precision                      :: dxdr
  double precision                      :: B(0:3), BD(3)
  double precision                      :: fac, fact
  double precision                      :: del
  double precision                      :: del_damp
  double precision                      :: del_damp2
  double precision                      :: gphi_i(3), gphi_j(3)
  double precision                      :: asq
  double precision                      :: expon, expo
  double precision                      :: clam3, clam5, clam7
  double precision                      :: delr3inv, delr5inv, delr7inv
  double precision                      :: Rn(1), Rn_1(4), Rn_2(10), Rn_3(20)

  integer                               :: img_j_vec(vs)
  double precision                      :: del_vec(3, vs)
  double precision                      :: delr2_vec(vs)
  
  logical                               :: single_thole_val

  if (iand(do_amoeba_direct_flag, proceed_induce) .ne. proceed_induce) return

  fac = 2.d0 * ewaldcof * ewaldcof
  del = 1.d0 / eedtbdns
  dxdr = ewaldcof
  
  if( thole_expon_coeff .gt. 0.000001) then 
    single_thole_val = .true.
  else
    single_thole_val = .false.
  endif  
    
  if( single_thole_val ) asq = thole_expon_coeff * img(img_i)%qterm
   
  sublst_head = 0

  do while (sublst_head .lt. pair_cnt)

    vec_max = min(pair_cnt - sublst_head, vs)

    vec_cnt = 0

#ifdef DIRFRC_COMTRANS
    if (common_tran .eq. 1) then

      do vec_idx = 1, vec_max
        img_j = ipairs_sublst(sublst_head + vec_idx)
        delx = img(img_j)%x + x_tran(1, 13)
        dely = img(img_j)%y + x_tran(2, 13)
        delz = img(img_j)%z + x_tran(3, 13)
        delr2 = delx * delx + dely * dely + delz * delz
        if (delr2 .lt. ee_dsum_cut2) then
          vec_cnt = vec_cnt + 1
          img_j_vec(vec_cnt) = img_j
          del_vec(1, vec_cnt) = delx
          del_vec(2, vec_cnt) = dely
          del_vec(3, vec_cnt) = delz
          delr2_vec(vec_cnt) = delr2
        end if
      end do

    else
#endif /* DIRFRC_COMTRANS */

      do vec_idx = 1, vec_max
        enc_img = ipairs_sublst(sublst_head + vec_idx)
        img_j = iand(enc_img, mask27)
        itran = ishft(enc_img, -27)
        delx = img(img_j)%x + x_tran(1, itran)
        dely = img(img_j)%y + x_tran(2, itran)
        delz = img(img_j)%z + x_tran(3, itran)
        delr2 = delx * delx + dely * dely + delz * delz
        if (delr2 .lt. ee_dsum_cut2) then
          vec_cnt = vec_cnt + 1
          img_j_vec(vec_cnt) = img_j
          del_vec(1, vec_cnt) = delx
          del_vec(2, vec_cnt) = dely
          del_vec(3, vec_cnt) = delz
          delr2_vec(vec_cnt) = delr2
        end if
      end do

#ifdef DIRFRC_COMTRANS
    end if
#endif /* DIRFRC_COMTRANS */

    do vec_idx = 1, vec_cnt  ! some over non-bonded list neighbours within electrostatic cutoff r_ij < ee_dsum_cut

      img_j = img_j_vec(vec_idx)
      delx = del_vec(1, vec_idx)
      dely = del_vec(2, vec_idx)
      delz = del_vec(3, vec_idx)
      delr2 = delr2_vec(vec_idx) ! r_ij^2
      atm_j = img_atm_map(img_j)

      if ((is_polarizable(atm_i) .or. is_polarizable(atm_j))) then
        !---------------------------------------------------------
        ! McMurchie-Davidson recursion for interaction tensor:
        ! interaction at point charge level given by complementary Boys
        ! B(0) = int_0^1 exp(-pr^2^t^2)dt
        ! complementary Boys is BC(0) = 1/r - B(0)
        ! (d/dr) B(0) = (-2p)*r int_0^1 t^2*exp(-pr^2t^2)dt = (-2p)*r B(1)
        ! and so if R(0,0,0,n) = (-2p)^n B(n) then we have
        ! (1/r)*(d/dr) R(0,0,0,n) = R(0,0,0,n+1)
        ! Now let R(t,u,v,n) = (d/dx)^t (d/dy)^u (d/dz)^v R(0,0,0,n)
        ! Then e.g. R(t+1,u,v,n) = t*R(t-1,u,v,n+1) + x*R(t,u,v,n)
        ! proof:  
        ! R(t+1,u,v,n) = (d/dx)^t+1 (d/dy)^u (d/dz)^v R(0,0,0,n)
        !              = (d/dx)^t (d/dy)^u (d/dz)^v (d/dx) R(0,0,0,n)
        !              = (d/dx)^t (d/dy)^u (d/dz)^v x*(1/r)(d/dr)R(0,0,0,n)
        !              = (d/dx)^t (d/dy)^u (d/dz)^v [x*R(0,0,0,n+1)]
        !              = (d/dx)^t [x*R(0,u,v,n+1)]
        !              = t*R(t-1,u,v,n+1) + x*R(t,u,v,n+1) (Leibniz)
        ! similar recursions hold for R(t,u+1,v,n),R(t,u,v+1,n)
        ! R(t,u+1,v,n) = u*R(t,u-1,v,n+1) + y*R(t,u,v,n+1)
        ! R(t,u,v+1,n) = v*R(t,u,v-1,n+1) + z*R(t,u,v,n+1)
        ! below array is packed---hence use of Ind_tuv
        ! Rn(Ind_tuv) denotes R(t,u,v,n)
        ! Due to its form--we recur downwards in n
        !---------------------------------------------------------
        ! top n is 3 for dipole fields
        ! get boys and R(0,0,0,n), n = 0,1,2,3

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
         
        ! TD Got the idea for B_l from Walter Smith's CCP5 article 1982
        ! Ewald for point multipoles
        ! B_l satisfies grad_i B_l(|r_j - r_i|) = (r_j - r_i)B_{l+1}(|r_j-r_i|)
        ! grad_j B_l(|r_j - r_i|) = -grad_i B_l(|r_j - r_i|)
         
        B(0) = switch * delr * delr2inv
        fact = d_switch_dx * dxdr
        B(1) = (B(0) - fact) * delr2inv
        fact = fac * fact
        B(2) = (3.d0 * B(1) - fact) * delr2inv
        fact = fac * fact
        B(3) = (5.d0 * B(2) - fact) * delr2inv

        if (delr2 .lt. ee_damped_cut2) then
          !-------------------------------------------------------
          ! McMurchie-Davidson holds for damped tensor as well---in fact,
          ! BD below satisfies BD(n+1) = (1/r)(d/dr)BD(n)
          ! RD(0,0,0,n) = BD(n), n = 1,2,3
          ! RD(t,u,v,n) = (d/dx)^t (d/dy)^u (d/dz)^v R(0,0,0,n)
          !-------------------------------------------------------
          delr3inv = delr2inv / delr
          delr5inv = delr3inv * delr2inv
          delr7inv = delr5inv * delr2inv
          if( single_thole_val ) then
             expon = asq * delr2 * delr * img(img_j)%qterm
!             write(*,*) "am_direct_permfield_i  thole_expon_coeff = ", thole_expon_coeff
          else                       
             expon = screen_polar(atm_i)*screen_polar(atm_j)*delr2*delr * img(img_i)%qterm * img(img_j)%qterm
!             write(*,*) "am_direct_permfield_i  screen_polar(atm_i)*screen_polar(atm_j) = ", screen_polar(atm_i)*screen_polar(atm_j)
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
          ! correct the boys factors by damped factors
          ! ewald dsum tensor and damped tensor both satisfy 
          ! McMurchie-Davidson recursion-thus their sum does as well 
          ! use recur for sum--starting with summed factors
          B(1) = B(1) - BD(1)
          B(2) = B(2) - BD(2)
          B(3) = B(3) - BD(3)
        end if !(delr2 < damped_cut2) then
        ! negate the odd order boys factors
        B(1) = -B(1)
        B(3) = -B(3)
        n = 3
        Rn(Ind_000) = B(n)
        Rn_1(Ind_000) = B(n-1)
        Rn_1(Ind_100) = delx*Rn(Ind_000)
        Rn_1(Ind_010) = dely*Rn(Ind_000)
        Rn_1(Ind_001) = delz*Rn(Ind_000)
        Rn_2(Ind_000) = B(n-2)
        Rn_2(Ind_100) = delx*Rn_1(Ind_000)
        Rn_2(Ind_010) = dely*Rn_1(Ind_000)
        Rn_2(Ind_001) = delz*Rn_1(Ind_000)
        Rn_2(Ind_200) = Rn_1(Ind_000) + delx*Rn_1(Ind_100)
        Rn_2(Ind_020) = Rn_1(Ind_000) + dely*Rn_1(Ind_010)
        Rn_2(Ind_002) = Rn_1(Ind_000) + delz*Rn_1(Ind_001)
        Rn_2(Ind_110) = delx*Rn_1(Ind_010)
        Rn_2(Ind_101) = delx*Rn_1(Ind_001)
        Rn_2(Ind_011) = dely*Rn_1(Ind_001)
        !Rn_3(Ind_000) = B(n-3) --don't need this one in any field calc
        Rn_3(Ind_100) = delx*Rn_2(Ind_000)
        Rn_3(Ind_010) = dely*Rn_2(Ind_000)
        Rn_3(Ind_001) = delz*Rn_2(Ind_000)
        Rn_3(Ind_200) = Rn_2(Ind_000) + delx*Rn_2(Ind_100)
        Rn_3(Ind_020) = Rn_2(Ind_000) + dely*Rn_2(Ind_010)
        Rn_3(Ind_002) = Rn_2(Ind_000) + delz*Rn_2(Ind_001)
        Rn_3(Ind_110) = delx*Rn_2(Ind_010)
        Rn_3(Ind_101) = delx*Rn_2(Ind_001)
        Rn_3(Ind_011) = dely*Rn_2(Ind_001)
        Rn_3(Ind_300) = 2.d0*Rn_2(Ind_100) + delx*Rn_2(Ind_200)
        Rn_3(Ind_030) = 2.d0*Rn_2(Ind_010) + dely*Rn_2(Ind_020)
        Rn_3(Ind_003) = 2.d0*Rn_2(Ind_001) + delz*Rn_2(Ind_002)
        Rn_3(Ind_210) = dely*Rn_2(Ind_200)
        Rn_3(Ind_201) = delz*Rn_2(Ind_200)
        Rn_3(Ind_120) = delx*Rn_2(Ind_020)
        Rn_3(Ind_021) = delz*Rn_2(Ind_020)
        Rn_3(Ind_102) = delx*Rn_2(Ind_002)
        Rn_3(Ind_012) = dely*Rn_2(Ind_002)
        Rn_3(Ind_111) = delx*Rn_2(Ind_011)
        ! UNDER MPI ONLY tensors for which you own img_i (img of atm_i) are
        ! stored!
        ! store the dipole-dipole component
        num_tensor = num_tensor + 1
        if (num_tensor .gt. size_dipole_dipole_list) then
          write( error_msg,*)'num_tensor = ',num_tensor, " > ", "size_dipole_dipole_list = ",size_dipole_dipole_list
          call mol_mech_error 
        end if
        !if( num_tensor .eq. 1 ) then
        !    write(*,*) " am_direct_permfield_i() pt 1, num_tensor, atm_i, atm_j ", num_tensor, atm_i, atm_j  
        !    write(*,'(" img(img_i)%qterm, img(img_j)%qterm =", 2(F14.8,1X) )') img(img_i)%qterm,img(img_j)%qterm
        !    write(*,'(" switch,eedtbdns,dxdr =", 3(F14.8,1X) )') switch,eedtbdns,dxdr
        !    write(*,'(" B(1),B(2),B(3) =", 3(F14.8,1X) )') B(1),B(2),B(3)
        !    write(*,'(" delx, dely, delz =", 3(F14.8,1X) )') delx, dely, delz
        !    write(*,'(" Rn_2(Ind_000), Rn_2(Ind_100) =", 2(F14.8,1X) )') Rn_2(Ind_000), Rn_2(Ind_100)
        !    write(*,'(" Rn_2(Ind_010), Rn_2(Ind_001) =", 2(F14.8,1X) )') Rn_2(Ind_010), Rn_2(Ind_001)
        !    write(*,'(" Rn_3(Ind_200), Rn_3(Ind_020), Rn_3(Ind_002) =", 3(F14.8,1X) )') Rn_3(Ind_200), Rn_3(Ind_020), Rn_3(Ind_002)
        !endif
        dipole_dipole_list(1, num_tensor) = atm_i
        dipole_dipole_list(2, num_tensor) = atm_j
        dipole_dipole_tensor(1, num_tensor) = Rn_3(Ind_200)
        dipole_dipole_tensor(2, num_tensor) = Rn_3(Ind_110)
        dipole_dipole_tensor(3, num_tensor) = Rn_3(Ind_101)
        dipole_dipole_tensor(4, num_tensor) = Rn_3(Ind_020)
        dipole_dipole_tensor(5, num_tensor) = Rn_3(Ind_011)
        dipole_dipole_tensor(6, num_tensor) = Rn_3(Ind_002)
        ! next do the field components
        if (is_polarizable(atm_i)) then
          gphi_i(1) = Rn_3(Ind_100)*global_multipole(Ind_000,atm_j) + &
                      Rn_3(Ind_200)*global_multipole(Ind_100,atm_j) + &
                      Rn_3(Ind_110)*global_multipole(Ind_010,atm_j) + &
                      Rn_3(Ind_101)*global_multipole(Ind_001,atm_j) + &
                      Rn_3(Ind_300)*global_multipole(Ind_200,atm_j) + &
                      Rn_3(Ind_120)*global_multipole(Ind_020,atm_j) + &
                      Rn_3(Ind_102)*global_multipole(Ind_002,atm_j) + &
                      Rn_3(Ind_210)*global_multipole(Ind_110,atm_j) + &
                      Rn_3(Ind_201)*global_multipole(Ind_101,atm_j) + &
                      Rn_3(Ind_111)*global_multipole(Ind_011,atm_j)
          gphi_i(2) = Rn_3(Ind_010)*global_multipole(Ind_000,atm_j) + &
                      Rn_3(Ind_110)*global_multipole(Ind_100,atm_j) + &
                      Rn_3(Ind_020)*global_multipole(Ind_010,atm_j) + &
                      Rn_3(Ind_011)*global_multipole(Ind_001,atm_j) + &
                      Rn_3(Ind_210)*global_multipole(Ind_200,atm_j) + &
                      Rn_3(Ind_030)*global_multipole(Ind_020,atm_j) + &
                      Rn_3(Ind_012)*global_multipole(Ind_002,atm_j) + &
                      Rn_3(Ind_120)*global_multipole(Ind_110,atm_j) + &
                      Rn_3(Ind_111)*global_multipole(Ind_101,atm_j) + &
                      Rn_3(Ind_021)*global_multipole(Ind_011,atm_j)
          gphi_i(3) = Rn_3(Ind_001)*global_multipole(Ind_000,atm_j) + &
                      Rn_3(Ind_101)*global_multipole(Ind_100,atm_j) + &
                      Rn_3(Ind_011)*global_multipole(Ind_010,atm_j) + &
                      Rn_3(Ind_002)*global_multipole(Ind_001,atm_j) + &
                      Rn_3(Ind_201)*global_multipole(Ind_200,atm_j) + &
                      Rn_3(Ind_021)*global_multipole(Ind_020,atm_j) + &
                      Rn_3(Ind_003)*global_multipole(Ind_002,atm_j) + &
                      Rn_3(Ind_111)*global_multipole(Ind_110,atm_j) + &
                      Rn_3(Ind_102)*global_multipole(Ind_101,atm_j) + &
                      Rn_3(Ind_012)*global_multipole(Ind_011,atm_j)
          ! minus sign due to derivs wrt crds of i
          gradphi(1,atm_i) = gradphi(1,atm_i) - gphi_i(1)
          gradphi(2,atm_i) = gradphi(2,atm_i) - gphi_i(2)
          gradphi(3,atm_i) = gradphi(3,atm_i) - gphi_i(3)
          
#ifdef VAR_POLAR         
          if( damp_polar_strength(atm_j) > 0.0 .and. damp_polar_sensitivity(atm_i) > 0.0 ) then
              del_damp = delr - (damp_polar_rad(atm_i) + damp_polar_rad(atm_j)) 
!              del_damp2 = del_damp ** 6  ! IGOR CHANGING_DAMPING POWER
               del_damp2 = del_damp ** ndamp
!              write(*,'("am_direct_permfield_i pt1 ",2I5," strength= ",F6.4," sensitivity= ",F6.4," del_damp= ",F8.4)') &
!                    atm_j, atm_i, damp_polar_strength(atm_j), damp_polar_sensitivity(atm_i),del_damp
              polarizability_corr(atm_i) = polarizability_corr(atm_i) + damp_polar_strength(atm_j) * damp_polar_sensitivity(atm_i)/del_damp2
          endif 
#endif
        end if ! ( is_polarizable(atm_i) ) then
  
        if (is_polarizable(atm_j)) then
        ! note negative contribs for dipoles of i since d/dx_i = -d/delx)
        ! (look at contributions to electrostatic potential at atm_j -these will
        !  contain negative contributions due to dipolar contribs at i--the
        !  extra deriv in tensor due to grad at atm_j has no negative signs)
          gphi_j(1) = Rn_3(Ind_100)*global_multipole(Ind_000,atm_i) - &
                      Rn_3(Ind_200)*global_multipole(Ind_100,atm_i) - &
                      Rn_3(Ind_110)*global_multipole(Ind_010,atm_i) - &
                      Rn_3(Ind_101)*global_multipole(Ind_001,atm_i) + &
                      Rn_3(Ind_300)*global_multipole(Ind_200,atm_i) + &
                      Rn_3(Ind_120)*global_multipole(Ind_020,atm_i) + &
                      Rn_3(Ind_102)*global_multipole(Ind_002,atm_i) + &
                      Rn_3(Ind_210)*global_multipole(Ind_110,atm_i) + &
                      Rn_3(Ind_201)*global_multipole(Ind_101,atm_i) + &
                      Rn_3(Ind_111)*global_multipole(Ind_011,atm_i)
          gphi_j(2) = Rn_3(Ind_010)*global_multipole(Ind_000,atm_i) - &
                      Rn_3(Ind_110)*global_multipole(Ind_100,atm_i) - &
                      Rn_3(Ind_020)*global_multipole(Ind_010,atm_i) - &
                      Rn_3(Ind_011)*global_multipole(Ind_001,atm_i) + &
                      Rn_3(Ind_210)*global_multipole(Ind_200,atm_i) + &
                      Rn_3(Ind_030)*global_multipole(Ind_020,atm_i) + &
                      Rn_3(Ind_012)*global_multipole(Ind_002,atm_i) + &
                      Rn_3(Ind_120)*global_multipole(Ind_110,atm_i) + &
                      Rn_3(Ind_111)*global_multipole(Ind_101,atm_i) + &
                      Rn_3(Ind_021)*global_multipole(Ind_011,atm_i)
          gphi_j(3) = Rn_3(Ind_001)*global_multipole(Ind_000,atm_i) - &
                      Rn_3(Ind_101)*global_multipole(Ind_100,atm_i) - &
                      Rn_3(Ind_011)*global_multipole(Ind_010,atm_i) - &
                      Rn_3(Ind_002)*global_multipole(Ind_001,atm_i) + &
                      Rn_3(Ind_201)*global_multipole(Ind_200,atm_i) + &
                      Rn_3(Ind_021)*global_multipole(Ind_020,atm_i) + &
                      Rn_3(Ind_003)*global_multipole(Ind_002,atm_i) + &
                      Rn_3(Ind_111)*global_multipole(Ind_110,atm_i) + &
                      Rn_3(Ind_102)*global_multipole(Ind_101,atm_i) + &
                      Rn_3(Ind_012)*global_multipole(Ind_011,atm_i)
          gradphi(1, atm_j) = gradphi(1, atm_j) + gphi_j(1)
          gradphi(2, atm_j) = gradphi(2, atm_j) + gphi_j(2)
          gradphi(3, atm_j) = gradphi(3, atm_j) + gphi_j(3)
          
#ifdef VAR_POLAR
          if( damp_polar_strength(atm_i) > 0.0 .and. damp_polar_sensitivity(atm_j) > 0.0 ) then
              del_damp = delr - (damp_polar_rad(atm_i) + damp_polar_rad(atm_j)) 
!              del_damp2 = del_damp ** 6  ! IGOR CHANGING_DAMPING POWER
               del_damp2 = del_damp ** ndamp
!              write(*,'("am_direct_permfield_i pt 2 ",2I5," strength= ",F6.4," sensitivity= ",F6.4," del_damp= ",F8.4)') &
!                       atm_i, atm_j, damp_polar_strength(atm_i), damp_polar_sensitivity(atm_j), del_damp
              polarizability_corr(atm_j) = polarizability_corr(atm_j) + damp_polar_strength(atm_i) * damp_polar_sensitivity(atm_j)/del_damp2
          endif 
#endif

        end if ! (is_polarizable(atm_j)) then
  
      end if

    end do

    sublst_head = sublst_head + vec_max

  end do

  return

end subroutine am_direct_permfield_i

end subroutine am_direct_permfield

!*******************************************************************************!
! Subroutine:  am_direct_ene_frc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_direct_ene_frc(ipairs, img, tranvec, crd, ind_dip_d, ind_dip_p, &
                             ene_perm, ene_ind, ene_vdw, frc, img_frc, virial, &
                             img_atm_map)

  use amoeba_multipoles_mod, only : global_multipole
  use amoeba_induced_mod, only : is_polarizable, screen_polar
  use prmtop_dat_mod
  use mdin_amoeba_dat_mod, only : ee_dsum_cut, ee_damped_cut, thole_expon_coeff
  use amoeba_vdw_mod
  use mdin_ewald_dat_mod
  use pme_direct_mod
  use img_mod
  use timers_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer                               :: ipairs(*)
  type(img_rec), intent(in)             :: img(*)
  double precision, intent(in)          :: tranvec(1:3, 0:17)
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in)          :: ind_dip_d(3, *)
  double precision, intent(in)          :: ind_dip_p(3, *)
  double precision, intent(in out)      :: ene_perm
  double precision, intent(in out)      :: ene_ind
  double precision, intent(in out)      :: ene_vdw
  double precision, intent(in out)      :: frc(3, *)
  double precision, intent(in out)      :: img_frc(3, *)
  double precision, intent(in out)      :: virial(3, 3)
  integer, intent(in)                   :: img_atm_map(*)

! Local variables:

  double precision                      :: x_i, y_i, z_i
  double precision                      :: x_tran(1:3, 0:17)
  double precision                      :: ee_dsum_cut2
  double precision                      :: ee_damped_cut2
  integer                               :: i
  integer                               :: ipairs_idx
  integer                               :: img_i
  integer                               :: pair_cnt
  integer                               :: atm_i
#ifdef DIRFRC_COMTRANS
  ! flag - 1 if translation not needed
  integer                               :: common_tran
#endif /* DIRFRC_COMTRANS */

  ee_dsum_cut2 = ee_dsum_cut * ee_dsum_cut
  ee_damped_cut2 = ee_damped_cut * ee_damped_cut

  ene_perm = 0.d0
  ene_ind = 0.d0
  ene_vdw = 0.d0

  ipairs_idx = 1

  do img_i = my_img_lo, my_img_hi

#ifdef DIRFRC_COMTRANS
    ! Common translation (ie. no translation) flag is packed at
    ! the front of each sublist followed by the count(s) of sublist
    ! image pair entries.

    common_tran = ipairs(ipairs_idx)
    ipairs_idx = ipairs_idx + 1
#endif /* DIRFRC_COMTRANS */
    
    ! Electrostatic evaluation-only count followed by
    ! full evaluation count packed at the front of each pair sublist.

    pair_cnt = ipairs(ipairs_idx) + ipairs(ipairs_idx + 1)
    ipairs_idx = ipairs_idx + 2

    if (pair_cnt .gt. 0) then

      x_i = img(img_i)%x
      y_i = img(img_i)%y
      z_i = img(img_i)%z

#ifdef DIRFRC_COMTRANS
      if (common_tran .eq. 0) then
#endif /* DIRFRC_COMTRANS */
        ! We need all the translation vectors:
        do i = 0, 17
          x_tran(1, i) = tranvec(1, i) - x_i
          x_tran(2, i) = tranvec(2, i) - y_i
          x_tran(3, i) = tranvec(3, i) - z_i
        end do
#ifdef DIRFRC_COMTRANS
      else
        ! Just put the x,y,z values in the middle cell
        x_tran(1, 13) = - x_i
        x_tran(2, 13) = - y_i
        x_tran(3, 13) = - z_i
      end if
#endif /* DIRFRC_COMTRANS */

      atm_i = img_atm_map(img_i)

      call am_direct_ene_frc_i(atm_i, img, ipairs(ipairs_idx), x_tran, &
                               pair_cnt, ew_coeff, eedtbdns, gbl_eed_cub, &
                               ee_dsum_cut2, ee_damped_cut2, &
                               thole_expon_coeff, screen_polar, &
                               is_polarizable, ind_dip_d, ind_dip_p, &
                               global_multipole, ene_perm, &
                               ene_ind, img_frc, virial, img_atm_map)

      call am_vdw_direct_ene_frc_i(atm_i, img, ipairs(ipairs_idx), x_tran, &
                                   pair_cnt, crd, ene_vdw, frc, img_frc, &
                                   virial, img_atm_map)

      ipairs_idx = ipairs_idx + pair_cnt

    end if
  end do

  call update_pme_time(dir_frc_sum_timer)

  return

contains

!*******************************************************************************!
! Subroutine:  am_direct_ene_frc_i
!
! Description: part electrostatic energy for direct sum of pme expansion for atom 
!              for consequitive atom images in the pairlist
!
!*******************************************************************************

subroutine am_direct_ene_frc_i(atm_i, img, ipairs_sublst, x_tran, &
                               pair_cnt, ewaldcof, eedtbdns, eed_cub, &
                               ee_dsum_cut2, ee_damped_cut2, &
                               thole_expon_coeff, screen_polar, &
                               is_polarizable, ind_dip_d, ind_dip_p, &
                               global_multipole, ene_perm, &
                               ene_ind, img_frc, virial, img_atm_map)

  use amoeba_multipoles_mod, only : coulomb_const_kcal_per_mole, torque_field
  use amoeba_flags_mod
  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: atm_i
  type(img_rec), intent(in)             :: img(*)
  integer                               :: ipairs_sublst(*)
  double precision, intent(in)          :: x_tran(1:3, 0:17)
  integer, intent(in)                   :: pair_cnt
  double precision, intent(in)          :: ewaldcof
  double precision, intent(in)          :: eedtbdns
  double precision, intent(in)          :: eed_cub(4, *)
  double precision, intent(in)          :: ee_dsum_cut2
  double precision, intent(in)          :: ee_damped_cut2
  double precision, intent(in)          :: thole_expon_coeff
  double precision, intent(in)          :: screen_polar(*)
  logical, intent(in)                   :: is_polarizable(*)
  double precision, intent(in)          :: ind_dip_d(3, *)
  double precision, intent(in)          :: ind_dip_p(3, *)
  double precision, intent(in)          :: global_multipole(10, *)
  double precision, intent(in out)      :: ene_perm       ! energy of permanent multipoles ??
  double precision, intent(in out)      :: ene_ind        ! polarization energy ??
  double precision, intent(in out)      :: img_frc(3, *)  ! forces ??? 
  double precision, intent(in out)      :: virial(3, 3)
  integer, intent(in)                   :: img_atm_map(*)

! Local variables:

  integer, parameter            :: vs = 128      ! vector size
  integer, parameter            :: mask27 = Z"07FFFFFF"
  double precision, parameter   :: const1 = 0.6d0
  double precision, parameter   :: const2 = 18.d0 / 35.d0
  double precision, parameter   :: const3 = 9.d0 / 35.d0
  double precision, parameter   :: third = 1.d0 / 3.d0

  integer                       :: atm_j, img_j, enc_img
  integer                       :: ind, jj, n
  integer                       :: itran
  integer                       :: sublst_head
  integer                       :: vec_idx, vec_cnt, vec_max
  double precision              :: delx, dely, delz
  double precision              :: delr, delr2, delr2inv
  double precision              :: x, dx, switch, d_switch_dx
  double precision              :: dxdr
  double precision              :: vxx, vxy, vxz, vyx, vyy, vyz, vzx, vzy, vzz
  double precision              :: B(0:5), BD(5)
  double precision              :: fac, fact
  double precision              :: del
  double precision              :: gm(10), phi(20), gmj(10), gmi(10), tmi(10)
  double precision              :: e_pp, e_ind, g_pp(3), g_ind(3)
  double precision              :: e_add
  double precision              :: i_di(3), i_pi(3), i_mi(3)
  double precision              :: i_dj(3), i_pj(3), i_mj(3)
  double precision              :: asq, expon, expo, clam3, clam5, clam7, clam9
  double precision              :: delr3inv, delr5inv, delr7inv, delr9inv
  double precision              :: Rn(1), Rn_1(4), Rn_2(10), Rn_3(20), Rn_4(35)
  double precision              :: Rn_5(56)

  integer                       :: img_j_vec(vs)
  double precision              :: del_vec(3, vs)
  double precision              :: delr2_vec(vs)
  
  logical                       :: single_thole_val

  if (iand(do_amoeba_direct_flag, proceed_postinduce) .ne. proceed_postinduce) &
    return

  fac = 2.d0 * ewaldcof * ewaldcof
  del = 1.d0 / eedtbdns
  dxdr = ewaldcof
  
  if( thole_expon_coeff .gt. 0.000001) then 
    single_thole_val = .true.
  else
    single_thole_val = .false.
  endif  
    
  if( single_thole_val ) asq = thole_expon_coeff * img(img_i)%qterm
  
  
  vxx = 0.d0
  vxy = 0.d0
  vxz = 0.d0
  vyx = 0.d0
  vyy = 0.d0
  vyz = 0.d0
  vzx = 0.d0
  vzy = 0.d0
  vzz = 0.d0

  sublst_head = 0

  do while (sublst_head .lt. pair_cnt)

    vec_max = min(pair_cnt - sublst_head, vs)

    vec_cnt = 0

#ifdef DIRFRC_COMTRANS
    if (common_tran .eq. 1) then

      do vec_idx = 1, vec_max
        img_j = ipairs_sublst(sublst_head + vec_idx)
        delx = img(img_j)%x + x_tran(1, 13)
        dely = img(img_j)%y + x_tran(2, 13)
        delz = img(img_j)%z + x_tran(3, 13)
        delr2 = delx * delx + dely * dely + delz * delz
        if (delr2 .lt. ee_dsum_cut2) then
          vec_cnt = vec_cnt + 1
          img_j_vec(vec_cnt) = img_j
          del_vec(1, vec_cnt) = delx
          del_vec(2, vec_cnt) = dely
          del_vec(3, vec_cnt) = delz
          delr2_vec(vec_cnt) = delr2
        end if
      end do

    else
#endif /* DIRFRC_COMTRANS */

      do vec_idx = 1, vec_max
        enc_img = ipairs_sublst(sublst_head + vec_idx)
        img_j = iand(enc_img, mask27)
        itran = ishft(enc_img, -27)
        delx = img(img_j)%x + x_tran(1, itran)
        dely = img(img_j)%y + x_tran(2, itran)
        delz = img(img_j)%z + x_tran(3, itran)
        delr2 = delx * delx + dely * dely + delz * delz
        if (delr2 .lt. ee_dsum_cut2) then
          vec_cnt = vec_cnt + 1
          img_j_vec(vec_cnt) = img_j
          del_vec(1, vec_cnt) = delx
          del_vec(2, vec_cnt) = dely
          del_vec(3, vec_cnt) = delz
          delr2_vec(vec_cnt) = delr2
        end if
      end do

#ifdef DIRFRC_COMTRANS
    end if
#endif /* DIRFRC_COMTRANS */

    do vec_idx = 1, vec_cnt

      img_j = img_j_vec(vec_idx)
      delx = del_vec(1, vec_idx)
      dely = del_vec(2, vec_idx)
      delz = del_vec(3, vec_idx)
      delr2 = delr2_vec(vec_idx)
  
      !---------------------------------------------------------
      ! McMurchie-Davidson recursion for interaction tensor:
      ! interaction at point charge level given by complementary Boys
      ! B(0) = int_0^1 exp(-pr^2^t^2)dt
      ! complementary Boys is BC(0) = 1/r - B(0)
      ! (d/dr) B(0) = (-2p)*r int_0^1 t^2*exp(-pr^2t^2)dt = (-2p)*r B(1)
      ! and so if R(0,0,0,n) = (-2p)^n B(n) then we have
      ! (1/r)*(d/dr) R(0,0,0,n) = R(0,0,0,n+1)
      ! Now let R(t,u,v,n) = (d/dx)^t (d/dy)^u (d/dz)^v R(0,0,0,n)
      ! Then e.g. R(t+1,u,v,n) = t*R(t-1,u,v,n+1) + x*R(t,u,v,n)
      ! proof:  
      ! R(t+1,u,v,n) = (d/dx)^t+1 (d/dy)^u (d/dz)^v R(0,0,0,n)
      !              = (d/dx)^t (d/dy)^u (d/dz)^v (d/dx) R(0,0,0,n)
      !              = (d/dx)^t (d/dy)^u (d/dz)^v x*(1/r)(d/dr)R(0,0,0,n)
      !              = (d/dx)^t (d/dy)^u (d/dz)^v [x*R(0,0,0,n+1)]
      !              = (d/dx)^t [x*R(0,u,v,n+1)]
      !              = t*R(t-1,u,v,n+1) + x*R(t,u,v,n+1) (Leibniz)
      ! similar recursions hold for R(t,u+1,v,n),R(t,u,v+1,n)
      ! R(t,u+1,v,n) = u*R(t,u-1,v,n+1) + y*R(t,u,v,n+1)
      ! R(t,u,v+1,n) = v*R(t,u,v-1,n+1) + z*R(t,u,v,n+1)
      ! below array is packed---hence use of Ind_tuv
      ! Rn(Ind_tuv) denotes R(t,u,v,n)
      ! Due to its form--we recur downwards in n
      !---------------------------------------------------------
      ! top n is 5 for energy and forces
      ! get boys and R(0,0,0,n), n = 0,...,5

      atm_j = img_atm_map(img_j)

      delr = sqrt(delr2)
      delr2inv = 1.d0 / delr2
      x = dxdr * delr         

      ! Only supported option is cubic spline on switch...

      ind = eedtbdns * x + 1
      dx = x - (ind - 1.d0) * del
      switch = eed_cub(1, ind) + dx * (eed_cub(2, ind)+ &
               dx * (eed_cub(3, ind) + dx * eed_cub(4, ind) * third) * 0.5d0)

      d_switch_dx = eed_cub(2, ind) + dx * (eed_cub(3, ind)+ &
                    dx * eed_cub(4, ind) * 0.5d0)
        
      ! TD Got the idea for B_l from Walter Smith's CCP5 article 1982
      ! Ewald for point multipoles
      ! B_l satisfies grad_i B_l(|r_j - r_i|) = (r_j - r_i)B_{l+1}(|r_j-r_i|)
      ! grad_j B_l(|r_j - r_i|) = -grad_i B_l(|r_j - r_i|)
        
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
      ! negate the odd order boys factors
      B(1) = -B(1)
      B(3) = -B(3)
      B(5) = -B(5)
      n = 5
      Rn(Ind_000) = B(n)
      Rn_1(Ind_000) = B(n-1)
      Rn_1(Ind_100) = delx*Rn(Ind_000)
      Rn_1(Ind_010) = dely*Rn(Ind_000)
      Rn_1(Ind_001) = delz*Rn(Ind_000)
      Rn_2(Ind_000) = B(n-2)
      Rn_2(Ind_100) = delx*Rn_1(Ind_000)
      Rn_2(Ind_010) = dely*Rn_1(Ind_000)
      Rn_2(Ind_001) = delz*Rn_1(Ind_000)
      Rn_2(Ind_200) = Rn_1(Ind_000) + delx*Rn_1(Ind_100)
      Rn_2(Ind_020) = Rn_1(Ind_000) + dely*Rn_1(Ind_010)
      Rn_2(Ind_002) = Rn_1(Ind_000) + delz*Rn_1(Ind_001)
      Rn_2(Ind_110) = delx*Rn_1(Ind_010)
      Rn_2(Ind_101) = delx*Rn_1(Ind_001)
      Rn_2(Ind_011) = dely*Rn_1(Ind_001)
      Rn_3(Ind_000) = B(n-3) 
      Rn_3(Ind_100) = delx*Rn_2(Ind_000)
      Rn_3(Ind_010) = dely*Rn_2(Ind_000)
      Rn_3(Ind_001) = delz*Rn_2(Ind_000)
      Rn_3(Ind_200) = Rn_2(Ind_000) + delx*Rn_2(Ind_100)
      Rn_3(Ind_020) = Rn_2(Ind_000) + dely*Rn_2(Ind_010)
      Rn_3(Ind_002) = Rn_2(Ind_000) + delz*Rn_2(Ind_001)
      Rn_3(Ind_110) = delx*Rn_2(Ind_010)
      Rn_3(Ind_101) = delx*Rn_2(Ind_001)
      Rn_3(Ind_011) = dely*Rn_2(Ind_001)
      Rn_3(Ind_300) = 2.d0*Rn_2(Ind_100) + delx*Rn_2(Ind_200)
      Rn_3(Ind_030) = 2.d0*Rn_2(Ind_010) + dely*Rn_2(Ind_020)
      Rn_3(Ind_003) = 2.d0*Rn_2(Ind_001) + delz*Rn_2(Ind_002)
      Rn_3(Ind_210) = dely*Rn_2(Ind_200)
      Rn_3(Ind_201) = delz*Rn_2(Ind_200)
      Rn_3(Ind_120) = delx*Rn_2(Ind_020)
      Rn_3(Ind_021) = delz*Rn_2(Ind_020)
      Rn_3(Ind_102) = delx*Rn_2(Ind_002)
      Rn_3(Ind_012) = dely*Rn_2(Ind_002)
      Rn_3(Ind_111) = delx*Rn_2(Ind_011)
      Rn_4(Ind_000) = B(n-4) 
      Rn_4(Ind_100) = delx*Rn_3(Ind_000)
      Rn_4(Ind_010) = dely*Rn_3(Ind_000)
      Rn_4(Ind_001) = delz*Rn_3(Ind_000)
      Rn_4(Ind_200) = Rn_3(Ind_000) + delx*Rn_3(Ind_100)
      Rn_4(Ind_020) = Rn_3(Ind_000) + dely*Rn_3(Ind_010)
      Rn_4(Ind_002) = Rn_3(Ind_000) + delz*Rn_3(Ind_001)
      Rn_4(Ind_110) = delx*Rn_3(Ind_010)
      Rn_4(Ind_101) = delx*Rn_3(Ind_001)
      Rn_4(Ind_011) = dely*Rn_3(Ind_001)
      Rn_4(Ind_300) = 2.d0*Rn_3(Ind_100) + delx*Rn_3(Ind_200)
      Rn_4(Ind_030) = 2.d0*Rn_3(Ind_010) + dely*Rn_3(Ind_020)
      Rn_4(Ind_003) = 2.d0*Rn_3(Ind_001) + delz*Rn_3(Ind_002)
      Rn_4(Ind_210) = dely*Rn_3(Ind_200)
      Rn_4(Ind_201) = delz*Rn_3(Ind_200)
      Rn_4(Ind_120) = delx*Rn_3(Ind_020)
      Rn_4(Ind_021) = delz*Rn_3(Ind_020)
      Rn_4(Ind_102) = delx*Rn_3(Ind_002)
      Rn_4(Ind_012) = dely*Rn_3(Ind_002)
      Rn_4(Ind_111) = delx*Rn_3(Ind_011)
      Rn_4(Ind_400) = 3.d0*Rn_3(Ind_200) + delx*Rn_3(Ind_300)
      Rn_4(Ind_040) = 3.d0*Rn_3(Ind_020) + dely*Rn_3(Ind_030)
      Rn_4(Ind_004) = 3.d0*Rn_3(Ind_002) + delz*Rn_3(Ind_003)
      Rn_4(Ind_310) = dely*Rn_3(Ind_300)
      Rn_4(Ind_301) = delz*Rn_3(Ind_300)
      Rn_4(Ind_130) = delx*Rn_3(Ind_030)
      Rn_4(Ind_031) = delz*Rn_3(Ind_030)
      Rn_4(Ind_103) = delx*Rn_3(Ind_003)
      Rn_4(Ind_013) = dely*Rn_3(Ind_003)
      Rn_4(Ind_220) = Rn_3(Ind_020) + delx*Rn_3(Ind_120)
      Rn_4(Ind_202) = Rn_3(Ind_002) + delx*Rn_3(Ind_102)
      Rn_4(Ind_022) = Rn_3(Ind_002) + dely*Rn_3(Ind_012)
      Rn_4(Ind_211) = dely*Rn_3(Ind_201)
      Rn_4(Ind_121) = delx*Rn_3(Ind_021)
      Rn_4(Ind_112) = delx*Rn_3(Ind_012)
      Rn_5(Ind_000) = B(n-5) 
      Rn_5(Ind_100) = delx*Rn_4(Ind_000)
      Rn_5(Ind_010) = dely*Rn_4(Ind_000)
      Rn_5(Ind_001) = delz*Rn_4(Ind_000)
      Rn_5(Ind_200) = Rn_4(Ind_000) + delx*Rn_4(Ind_100)
      Rn_5(Ind_020) = Rn_4(Ind_000) + dely*Rn_4(Ind_010)
      Rn_5(Ind_002) = Rn_4(Ind_000) + delz*Rn_4(Ind_001)
      Rn_5(Ind_110) = delx*Rn_4(Ind_010)
      Rn_5(Ind_101) = delx*Rn_4(Ind_001)
      Rn_5(Ind_011) = dely*Rn_4(Ind_001)
      Rn_5(Ind_300) = 2.d0*Rn_4(Ind_100) + delx*Rn_4(Ind_200)
      Rn_5(Ind_030) = 2.d0*Rn_4(Ind_010) + dely*Rn_4(Ind_020)
      Rn_5(Ind_003) = 2.d0*Rn_4(Ind_001) + delz*Rn_4(Ind_002)
      Rn_5(Ind_210) = dely*Rn_4(Ind_200)
      Rn_5(Ind_201) = delz*Rn_4(Ind_200)
      Rn_5(Ind_120) = delx*Rn_4(Ind_020)
      Rn_5(Ind_021) = delz*Rn_4(Ind_020)
      Rn_5(Ind_102) = delx*Rn_4(Ind_002)
      Rn_5(Ind_012) = dely*Rn_4(Ind_002)
      Rn_5(Ind_111) = delx*Rn_4(Ind_011)
      Rn_5(Ind_400) = 3.d0*Rn_4(Ind_200) + delx*Rn_4(Ind_300)
      Rn_5(Ind_040) = 3.d0*Rn_4(Ind_020) + dely*Rn_4(Ind_030)
      Rn_5(Ind_004) = 3.d0*Rn_4(Ind_002) + delz*Rn_4(Ind_003)
      Rn_5(Ind_310) = dely*Rn_4(Ind_300)
      Rn_5(Ind_301) = delz*Rn_4(Ind_300)
      Rn_5(Ind_130) = delx*Rn_4(Ind_030)
      Rn_5(Ind_031) = delz*Rn_4(Ind_030)
      Rn_5(Ind_103) = delx*Rn_4(Ind_003)
      Rn_5(Ind_013) = dely*Rn_4(Ind_003)
      Rn_5(Ind_220) = Rn_4(Ind_020) + delx*Rn_4(Ind_120)
      Rn_5(Ind_202) = Rn_4(Ind_002) + delx*Rn_4(Ind_102)
      Rn_5(Ind_022) = Rn_4(Ind_002) + dely*Rn_4(Ind_012)
      Rn_5(Ind_211) = dely*Rn_4(Ind_201)
      Rn_5(Ind_121) = delx*Rn_4(Ind_021)
      Rn_5(Ind_112) = delx*Rn_4(Ind_012)
      Rn_5(Ind_500) = 4.d0*Rn_4(Ind_300) + delx*Rn_4(Ind_400)
      Rn_5(Ind_050) = 4.d0*Rn_4(Ind_030) + dely*Rn_4(Ind_040)
      Rn_5(Ind_005) = 4.d0*Rn_4(Ind_003) + delz*Rn_4(Ind_004)
      Rn_5(Ind_410) = dely*Rn_4(Ind_400)
      Rn_5(Ind_401) = delz*Rn_4(Ind_400)
      Rn_5(Ind_140) = delx*Rn_4(Ind_040)
      Rn_5(Ind_041) = delz*Rn_4(Ind_040)
      Rn_5(Ind_104) = delx*Rn_4(Ind_004)
      Rn_5(Ind_014) = dely*Rn_4(Ind_004)
      Rn_5(Ind_320) = Rn_4(Ind_300) + dely*Rn_4(Ind_310)
      Rn_5(Ind_302) = Rn_4(Ind_300) + delz*Rn_4(Ind_301)
      Rn_5(Ind_230) = Rn_4(Ind_030) + delx*Rn_4(Ind_130)
      Rn_5(Ind_032) = Rn_4(Ind_030) + delz*Rn_4(Ind_031)
      Rn_5(Ind_203) = Rn_4(Ind_003) + delx*Rn_4(Ind_103)
      Rn_5(Ind_023) = Rn_4(Ind_003) + dely*Rn_4(Ind_013)
      Rn_5(Ind_311) = dely*Rn_4(Ind_301)
      Rn_5(Ind_131) = delx*Rn_4(Ind_031)
      Rn_5(Ind_113) = delx*Rn_4(Ind_013)
      Rn_5(Ind_221) = delz*Rn_4(Ind_220)
      Rn_5(Ind_212) = dely*Rn_4(Ind_202)
      Rn_5(Ind_122) = delx*Rn_4(Ind_022)
      ! phi array (electrostatic potential at i due to j permanent mpoles
      ! and derivs of that esp wrt r_i)
      ! minus signs arise due to derivs of r_j - r_i wrt r_i
      do jj = 1, 10
        gmi(jj) = global_multipole(jj, atm_i)
        tmi(jj) = global_multipole(jj, atm_i) ! used for torque contrib
        gmj(jj) = global_multipole(jj, atm_j)
      end do
      do jj = 1, 3
        g_ind(jj) = 0.d0
        i_dj(jj) = ind_dip_d(jj, atm_j)
        i_di(jj) = ind_dip_d(jj, atm_i)
        i_pj(jj) = ind_dip_p(jj, atm_j)
        i_pi(jj) = ind_dip_p(jj, atm_i)
        i_mj(jj) = ind_dip_d(jj, atm_j) + ind_dip_p(jj, atm_j)
        i_mi(jj) = ind_dip_d(jj, atm_i) + ind_dip_p(jj, atm_i)
        tmi(jj + 1) = tmi(jj + 1) + 0.5d0 * i_mi(jj) !used for torque contrib
      end do
      ! initialize induction contributions
      e_ind = 0.d0
      phi(Ind_000) = Rn_5(Ind_000)*gmj(Ind_000)+Rn_5(Ind_100)*gmj(Ind_100)+ &   ! electrostatic potential at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_010)*gmj(Ind_010)+Rn_5(Ind_001)*gmj(Ind_001)+ &
                     Rn_5(Ind_200)*gmj(Ind_200)+Rn_5(Ind_020)*gmj(Ind_020)+ &
                     Rn_5(Ind_002)*gmj(Ind_002)+Rn_5(Ind_110)*gmj(Ind_110)+ &
                     Rn_5(Ind_101)*gmj(Ind_101)+Rn_5(Ind_011)*gmj(Ind_011)
      phi(Ind_100) = -(Rn_5(Ind_100)*gmj(Ind_000)+Rn_5(Ind_200)*gmj(Ind_100)+ &  ! dipole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_110)*gmj(Ind_010)+Rn_5(Ind_101)*gmj(Ind_001)+ &
                     Rn_5(Ind_300)*gmj(Ind_200)+Rn_5(Ind_120)*gmj(Ind_020)+ &
                     Rn_5(Ind_102)*gmj(Ind_002)+Rn_5(Ind_210)*gmj(Ind_110)+ &
                     Rn_5(Ind_201)*gmj(Ind_101)+Rn_5(Ind_111)*gmj(Ind_011))
      phi(Ind_010) = -(Rn_5(Ind_010)*gmj(Ind_000)+Rn_5(Ind_110)*gmj(Ind_100)+ &  ! dipole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_020)*gmj(Ind_010)+Rn_5(Ind_011)*gmj(Ind_001)+ &
                     Rn_5(Ind_210)*gmj(Ind_200)+Rn_5(Ind_030)*gmj(Ind_020)+ &
                     Rn_5(Ind_012)*gmj(Ind_002)+Rn_5(Ind_120)*gmj(Ind_110)+ &
                     Rn_5(Ind_111)*gmj(Ind_101)+Rn_5(Ind_021)*gmj(Ind_011))
      phi(Ind_001) = -(Rn_5(Ind_001)*gmj(Ind_000)+Rn_5(Ind_101)*gmj(Ind_100)+ &  ! dipole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_011)*gmj(Ind_010)+Rn_5(Ind_002)*gmj(Ind_001)+ &
                     Rn_5(Ind_201)*gmj(Ind_200)+Rn_5(Ind_021)*gmj(Ind_020)+ &
                     Rn_5(Ind_003)*gmj(Ind_002)+Rn_5(Ind_111)*gmj(Ind_110)+ &
                     Rn_5(Ind_102)*gmj(Ind_101)+Rn_5(Ind_012)*gmj(Ind_011))
      phi(Ind_200) = Rn_5(Ind_200)*gmj(Ind_000)+Rn_5(Ind_300)*gmj(Ind_100)+ &   ! quadrupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_210)*gmj(Ind_010)+Rn_5(Ind_201)*gmj(Ind_001)+ &
                     Rn_5(Ind_400)*gmj(Ind_200)+Rn_5(Ind_220)*gmj(Ind_020)+ &
                     Rn_5(Ind_202)*gmj(Ind_002)+Rn_5(Ind_310)*gmj(Ind_110)+ &
                     Rn_5(Ind_301)*gmj(Ind_101)+Rn_5(Ind_211)*gmj(Ind_011)
      phi(Ind_020) = Rn_5(Ind_020)*gmj(Ind_000)+Rn_5(Ind_120)*gmj(Ind_100)+ &   ! quadrupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_030)*gmj(Ind_010)+Rn_5(Ind_021)*gmj(Ind_001)+ &
                     Rn_5(Ind_220)*gmj(Ind_200)+Rn_5(Ind_040)*gmj(Ind_020)+ &
                     Rn_5(Ind_022)*gmj(Ind_002)+Rn_5(Ind_130)*gmj(Ind_110)+ &
                     Rn_5(Ind_121)*gmj(Ind_101)+Rn_5(Ind_031)*gmj(Ind_011)
      phi(Ind_002) = Rn_5(Ind_002)*gmj(Ind_000)+Rn_5(Ind_102)*gmj(Ind_100)+ &   ! quadrupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_012)*gmj(Ind_010)+Rn_5(Ind_003)*gmj(Ind_001)+ &
                     Rn_5(Ind_202)*gmj(Ind_200)+Rn_5(Ind_022)*gmj(Ind_020)+ &
                     Rn_5(Ind_004)*gmj(Ind_002)+Rn_5(Ind_112)*gmj(Ind_110)+ &
                     Rn_5(Ind_103)*gmj(Ind_101)+Rn_5(Ind_013)*gmj(Ind_011)
      phi(Ind_110) = Rn_5(Ind_110)*gmj(Ind_000)+Rn_5(Ind_210)*gmj(Ind_100)+ &   ! quadrupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_120)*gmj(Ind_010)+Rn_5(Ind_111)*gmj(Ind_001)+ &
                     Rn_5(Ind_310)*gmj(Ind_200)+Rn_5(Ind_130)*gmj(Ind_020)+ &
                     Rn_5(Ind_112)*gmj(Ind_002)+Rn_5(Ind_220)*gmj(Ind_110)+ &
                     Rn_5(Ind_211)*gmj(Ind_101)+Rn_5(Ind_121)*gmj(Ind_011)
      phi(Ind_101) = Rn_5(Ind_101)*gmj(Ind_000)+Rn_5(Ind_201)*gmj(Ind_100)+ &   ! quadrupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_111)*gmj(Ind_010)+Rn_5(Ind_102)*gmj(Ind_001)+ &
                     Rn_5(Ind_301)*gmj(Ind_200)+Rn_5(Ind_121)*gmj(Ind_020)+ &
                     Rn_5(Ind_103)*gmj(Ind_002)+Rn_5(Ind_211)*gmj(Ind_110)+ &
                     Rn_5(Ind_202)*gmj(Ind_101)+Rn_5(Ind_112)*gmj(Ind_011)
      phi(Ind_011) = Rn_5(Ind_011)*gmj(Ind_000)+Rn_5(Ind_111)*gmj(Ind_100)+ &   ! quadrupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_021)*gmj(Ind_010)+Rn_5(Ind_012)*gmj(Ind_001)+ &
                     Rn_5(Ind_211)*gmj(Ind_200)+Rn_5(Ind_031)*gmj(Ind_020)+ &
                     Rn_5(Ind_013)*gmj(Ind_002)+Rn_5(Ind_121)*gmj(Ind_110)+ &
                     Rn_5(Ind_112)*gmj(Ind_101)+Rn_5(Ind_022)*gmj(Ind_011)
      phi(Ind_300) = -(Rn_5(Ind_300)*gmj(Ind_000)+Rn_5(Ind_400)*gmj(Ind_100)+ &  ! octupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_310)*gmj(Ind_010)+Rn_5(Ind_301)*gmj(Ind_001)+ &
                     Rn_5(Ind_500)*gmj(Ind_200)+Rn_5(Ind_320)*gmj(Ind_020)+ &
                     Rn_5(Ind_302)*gmj(Ind_002)+Rn_5(Ind_410)*gmj(Ind_110)+ &
                     Rn_5(Ind_401)*gmj(Ind_101)+Rn_5(Ind_311)*gmj(Ind_011))
      phi(Ind_030) = -(Rn_5(Ind_030)*gmj(Ind_000)+Rn_5(Ind_130)*gmj(Ind_100)+ & ! octupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_040)*gmj(Ind_010)+Rn_5(Ind_031)*gmj(Ind_001)+ &
                     Rn_5(Ind_230)*gmj(Ind_200)+Rn_5(Ind_050)*gmj(Ind_020)+ &
                     Rn_5(Ind_032)*gmj(Ind_002)+Rn_5(Ind_140)*gmj(Ind_110)+ &
                     Rn_5(Ind_131)*gmj(Ind_101)+Rn_5(Ind_041)*gmj(Ind_011))
      phi(Ind_003) = -(Rn_5(Ind_003)*gmj(Ind_000)+Rn_5(Ind_103)*gmj(Ind_100)+ & ! octupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_013)*gmj(Ind_010)+Rn_5(Ind_004)*gmj(Ind_001)+ &
                     Rn_5(Ind_203)*gmj(Ind_200)+Rn_5(Ind_023)*gmj(Ind_020)+ &
                     Rn_5(Ind_005)*gmj(Ind_002)+Rn_5(Ind_113)*gmj(Ind_110)+ &
                     Rn_5(Ind_104)*gmj(Ind_101)+Rn_5(Ind_014)*gmj(Ind_011))
      phi(Ind_210) = -(Rn_5(Ind_210)*gmj(Ind_000)+Rn_5(Ind_310)*gmj(Ind_100)+ &  ! octupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_220)*gmj(Ind_010)+Rn_5(Ind_211)*gmj(Ind_001)+ &
                     Rn_5(Ind_410)*gmj(Ind_200)+Rn_5(Ind_230)*gmj(Ind_020)+ &
                     Rn_5(Ind_212)*gmj(Ind_002)+Rn_5(Ind_320)*gmj(Ind_110)+ &
                     Rn_5(Ind_311)*gmj(Ind_101)+Rn_5(Ind_221)*gmj(Ind_011))
      phi(Ind_201) = -(Rn_5(Ind_201)*gmj(Ind_000)+Rn_5(Ind_301)*gmj(Ind_100)+ &  ! octupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_211)*gmj(Ind_010)+Rn_5(Ind_202)*gmj(Ind_001)+ &
                     Rn_5(Ind_401)*gmj(Ind_200)+Rn_5(Ind_221)*gmj(Ind_020)+ &
                     Rn_5(Ind_203)*gmj(Ind_002)+Rn_5(Ind_311)*gmj(Ind_110)+ &
                     Rn_5(Ind_302)*gmj(Ind_101)+Rn_5(Ind_212)*gmj(Ind_011))
      phi(Ind_120) = -(Rn_5(Ind_120)*gmj(Ind_000)+Rn_5(Ind_220)*gmj(Ind_100)+ &  ! octupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_130)*gmj(Ind_010)+Rn_5(Ind_121)*gmj(Ind_001)+ &
                     Rn_5(Ind_320)*gmj(Ind_200)+Rn_5(Ind_140)*gmj(Ind_020)+ &
                     Rn_5(Ind_122)*gmj(Ind_002)+Rn_5(Ind_230)*gmj(Ind_110)+ &
                     Rn_5(Ind_221)*gmj(Ind_101)+Rn_5(Ind_131)*gmj(Ind_011))
      phi(Ind_021) = -(Rn_5(Ind_021)*gmj(Ind_000)+Rn_5(Ind_121)*gmj(Ind_100)+ &  ! octupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_031)*gmj(Ind_010)+Rn_5(Ind_022)*gmj(Ind_001)+ &
                     Rn_5(Ind_221)*gmj(Ind_200)+Rn_5(Ind_041)*gmj(Ind_020)+ &
                     Rn_5(Ind_023)*gmj(Ind_002)+Rn_5(Ind_131)*gmj(Ind_110)+ &
                     Rn_5(Ind_122)*gmj(Ind_101)+Rn_5(Ind_032)*gmj(Ind_011))
      phi(Ind_102) = -(Rn_5(Ind_102)*gmj(Ind_000)+Rn_5(Ind_202)*gmj(Ind_100)+ &  ! octupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_112)*gmj(Ind_010)+Rn_5(Ind_103)*gmj(Ind_001)+ &
                     Rn_5(Ind_302)*gmj(Ind_200)+Rn_5(Ind_122)*gmj(Ind_020)+ &
                     Rn_5(Ind_104)*gmj(Ind_002)+Rn_5(Ind_212)*gmj(Ind_110)+ &
                     Rn_5(Ind_203)*gmj(Ind_101)+Rn_5(Ind_113)*gmj(Ind_011))
      phi(Ind_012) = -(Rn_5(Ind_012)*gmj(Ind_000)+Rn_5(Ind_112)*gmj(Ind_100)+ &  ! octupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_022)*gmj(Ind_010)+Rn_5(Ind_013)*gmj(Ind_001)+ &
                     Rn_5(Ind_212)*gmj(Ind_200)+Rn_5(Ind_032)*gmj(Ind_020)+ &
                     Rn_5(Ind_014)*gmj(Ind_002)+Rn_5(Ind_122)*gmj(Ind_110)+ &
                     Rn_5(Ind_113)*gmj(Ind_101)+Rn_5(Ind_023)*gmj(Ind_011))
      phi(Ind_111) = -(Rn_5(Ind_111)*gmj(Ind_000)+Rn_5(Ind_211)*gmj(Ind_100)+ &  ! octupole electrostatic field at position of atm_i from multipoles of atm_j
                     Rn_5(Ind_121)*gmj(Ind_010)+Rn_5(Ind_112)*gmj(Ind_001)+ &
                     Rn_5(Ind_311)*gmj(Ind_200)+Rn_5(Ind_131)*gmj(Ind_020)+ &
                     Rn_5(Ind_113)*gmj(Ind_002)+Rn_5(Ind_221)*gmj(Ind_110)+ &
                     Rn_5(Ind_212)*gmj(Ind_101)+Rn_5(Ind_122)*gmj(Ind_011))
      e_pp = phi(Ind_000)*gmi(Ind_000) + phi(Ind_100)*gmi(Ind_100) + &          ! energy of interaction of multipoles of atm_i with the multipole field of atm_j
             phi(Ind_010)*gmi(Ind_010) + phi(Ind_001)*gmi(Ind_001) + &
             phi(Ind_200)*gmi(Ind_200) + phi(Ind_020)*gmi(Ind_020) + &
             phi(Ind_002)*gmi(Ind_002) + phi(Ind_110)*gmi(Ind_110) + &
             phi(Ind_101)*gmi(Ind_101) + phi(Ind_011)*gmi(Ind_011)

      ! gradient of e_pp wrt r_i

      g_pp(1) = phi(Ind_100)*gmi(Ind_000) + phi(Ind_200)*gmi(Ind_100) + &
                phi(Ind_110)*gmi(Ind_010) + phi(Ind_101)*gmi(Ind_001) + &
                phi(Ind_300)*gmi(Ind_200) + phi(Ind_120)*gmi(Ind_020) + &
                phi(Ind_102)*gmi(Ind_002) + phi(Ind_210)*gmi(Ind_110) + &
                phi(Ind_201)*gmi(Ind_101) + phi(Ind_111)*gmi(Ind_011)
      g_pp(2) = phi(Ind_010)*gmi(Ind_000) + phi(Ind_110)*gmi(Ind_100) + &
                phi(Ind_020)*gmi(Ind_010) + phi(Ind_011)*gmi(Ind_001) + &
                phi(Ind_210)*gmi(Ind_200) + phi(Ind_030)*gmi(Ind_020) + &
                phi(Ind_012)*gmi(Ind_002) + phi(Ind_120)*gmi(Ind_110) + &
                phi(Ind_111)*gmi(Ind_101) + phi(Ind_021)*gmi(Ind_011)
      g_pp(3) = phi(Ind_001)*gmi(Ind_000) + phi(Ind_101)*gmi(Ind_100) + &
                phi(Ind_011)*gmi(Ind_010) + phi(Ind_002)*gmi(Ind_001) + &
                phi(Ind_201)*gmi(Ind_200) + phi(Ind_021)*gmi(Ind_020) + &
                phi(Ind_003)*gmi(Ind_002) + phi(Ind_111)*gmi(Ind_110) + &
                phi(Ind_102)*gmi(Ind_101) + phi(Ind_012)*gmi(Ind_011)

      ! torque field at i due to permanent mpoles of j

      do jj = 1, 10
        torque_field(jj, atm_i) = torque_field(jj, atm_i) + phi(jj)
      end do

      if (is_polarizable(atm_i)) then

        e_ind = e_ind + 0.5d0 * &           ! interaction of polarizable dipole of atm_i with dipole field created by multipoles of atm_j       
                (phi(Ind_100) * i_di(1) + phi(Ind_010) * i_di(2) + &
                 phi(Ind_001) * i_di(3))
        
        g_ind(1) = g_ind(1) + 0.5d0 * &
                   (phi(Ind_200) * i_mi(1) + phi(Ind_110) * i_mi(2) + &
                    phi(Ind_101) * i_mi(3))

        g_ind(2) = g_ind(2) + 0.5d0 * &
                   (phi(Ind_110) * i_mi(1) + phi(Ind_020) * i_mi(2) + &
                    phi(Ind_011) * i_mi(3))

        g_ind(3) = g_ind(3) + 0.5d0 * &
                   (phi(Ind_101) * i_mi(1) + phi(Ind_011) * i_mi(2) + &
                    phi(Ind_002) * i_mi(3))

      end if ! is_polarizable(atm_i) 

      ! get the field at j due to permanent + induced mpoles of i 
      ! for torque at j
      ! electrostatic potential at j due to permanent + induced mpoles at i
      ! and derivatives of that with respect to r_j
      ! minus signs due to deriv of r_j-r_i wrt r_i entering into potential

      phi(Ind_000) = Rn_5(Ind_000)*tmi(Ind_000)-Rn_5(Ind_100)*tmi(Ind_100)- &
                     Rn_5(Ind_010)*tmi(Ind_010)-Rn_5(Ind_001)*tmi(Ind_001)+ &
                     Rn_5(Ind_200)*tmi(Ind_200)+Rn_5(Ind_020)*tmi(Ind_020)+ &
                     Rn_5(Ind_002)*tmi(Ind_002)+Rn_5(Ind_110)*tmi(Ind_110)+ &
                     Rn_5(Ind_101)*tmi(Ind_101)+Rn_5(Ind_011)*tmi(Ind_011)
      phi(Ind_100) = Rn_5(Ind_100)*tmi(Ind_000)-Rn_5(Ind_200)*tmi(Ind_100)- &
                     Rn_5(Ind_110)*tmi(Ind_010)-Rn_5(Ind_101)*tmi(Ind_001)+ &
                     Rn_5(Ind_300)*tmi(Ind_200)+Rn_5(Ind_120)*tmi(Ind_020)+ &
                     Rn_5(Ind_102)*tmi(Ind_002)+Rn_5(Ind_210)*tmi(Ind_110)+ &
                     Rn_5(Ind_201)*tmi(Ind_101)+Rn_5(Ind_111)*tmi(Ind_011)
      phi(Ind_010) = Rn_5(Ind_010)*tmi(Ind_000)-Rn_5(Ind_110)*tmi(Ind_100)- &
                     Rn_5(Ind_020)*tmi(Ind_010)-Rn_5(Ind_011)*tmi(Ind_001)+ &
                     Rn_5(Ind_210)*tmi(Ind_200)+Rn_5(Ind_030)*tmi(Ind_020)+ &
                     Rn_5(Ind_012)*tmi(Ind_002)+Rn_5(Ind_120)*tmi(Ind_110)+ &
                     Rn_5(Ind_111)*tmi(Ind_101)+Rn_5(Ind_021)*tmi(Ind_011)
      phi(Ind_001) = Rn_5(Ind_001)*tmi(Ind_000)-Rn_5(Ind_101)*tmi(Ind_100)- &
                     Rn_5(Ind_011)*tmi(Ind_010)-Rn_5(Ind_002)*tmi(Ind_001)+ &
                     Rn_5(Ind_201)*tmi(Ind_200)+Rn_5(Ind_021)*tmi(Ind_020)+ &
                     Rn_5(Ind_003)*tmi(Ind_002)+Rn_5(Ind_111)*tmi(Ind_110)+ &
                     Rn_5(Ind_102)*tmi(Ind_101)+Rn_5(Ind_012)*tmi(Ind_011)
      phi(Ind_200) = Rn_5(Ind_200)*tmi(Ind_000)-Rn_5(Ind_300)*tmi(Ind_100)- &
                     Rn_5(Ind_210)*tmi(Ind_010)-Rn_5(Ind_201)*tmi(Ind_001)+ &
                     Rn_5(Ind_400)*tmi(Ind_200)+Rn_5(Ind_220)*tmi(Ind_020)+ &
                     Rn_5(Ind_202)*tmi(Ind_002)+Rn_5(Ind_310)*tmi(Ind_110)+ &
                     Rn_5(Ind_301)*tmi(Ind_101)+Rn_5(Ind_211)*tmi(Ind_011)
      phi(Ind_020) = Rn_5(Ind_020)*tmi(Ind_000)-Rn_5(Ind_120)*tmi(Ind_100)- &
                     Rn_5(Ind_030)*tmi(Ind_010)-Rn_5(Ind_021)*tmi(Ind_001)+ &
                     Rn_5(Ind_220)*tmi(Ind_200)+Rn_5(Ind_040)*tmi(Ind_020)+ &
                     Rn_5(Ind_022)*tmi(Ind_002)+Rn_5(Ind_130)*tmi(Ind_110)+ &
                     Rn_5(Ind_121)*tmi(Ind_101)+Rn_5(Ind_031)*tmi(Ind_011)
      phi(Ind_002) = Rn_5(Ind_002)*tmi(Ind_000)-Rn_5(Ind_102)*tmi(Ind_100)- &
                     Rn_5(Ind_012)*tmi(Ind_010)-Rn_5(Ind_003)*tmi(Ind_001)+ &
                     Rn_5(Ind_202)*tmi(Ind_200)+Rn_5(Ind_022)*tmi(Ind_020)+ &
                     Rn_5(Ind_004)*tmi(Ind_002)+Rn_5(Ind_112)*tmi(Ind_110)+ &
                     Rn_5(Ind_103)*tmi(Ind_101)+Rn_5(Ind_013)*tmi(Ind_011)
      phi(Ind_110) = Rn_5(Ind_110)*tmi(Ind_000)-Rn_5(Ind_210)*tmi(Ind_100)- &
                     Rn_5(Ind_120)*tmi(Ind_010)-Rn_5(Ind_111)*tmi(Ind_001)+ &
                     Rn_5(Ind_310)*tmi(Ind_200)+Rn_5(Ind_130)*tmi(Ind_020)+ &
                     Rn_5(Ind_112)*tmi(Ind_002)+Rn_5(Ind_220)*tmi(Ind_110)+ &
                     Rn_5(Ind_211)*tmi(Ind_101)+Rn_5(Ind_121)*tmi(Ind_011)
      phi(Ind_101) = Rn_5(Ind_101)*tmi(Ind_000)-Rn_5(Ind_201)*tmi(Ind_100)- &
                     Rn_5(Ind_111)*tmi(Ind_010)-Rn_5(Ind_102)*tmi(Ind_001)+ &
                     Rn_5(Ind_301)*tmi(Ind_200)+Rn_5(Ind_121)*tmi(Ind_020)+ &
                     Rn_5(Ind_103)*tmi(Ind_002)+Rn_5(Ind_211)*tmi(Ind_110)+ &
                     Rn_5(Ind_202)*tmi(Ind_101)+Rn_5(Ind_112)*tmi(Ind_011)
      phi(Ind_011) = Rn_5(Ind_011)*tmi(Ind_000)-Rn_5(Ind_111)*tmi(Ind_100)- &
                     Rn_5(Ind_021)*tmi(Ind_010)-Rn_5(Ind_012)*tmi(Ind_001)+ &
                     Rn_5(Ind_211)*tmi(Ind_200)+Rn_5(Ind_031)*tmi(Ind_020)+ &
                     Rn_5(Ind_013)*tmi(Ind_002)+Rn_5(Ind_121)*tmi(Ind_110)+ &
                     Rn_5(Ind_112)*tmi(Ind_101)+Rn_5(Ind_022)*tmi(Ind_011)

      ! torque field at j due to permanent + induced mpoles of i

      do jj = 1, 10
        torque_field(jj, atm_j) = torque_field(jj, atm_j) + phi(jj)
      end do

      if (is_polarizable(atm_j)) then

      ! phi array (electrostatic potential at i due to j induced moments
      ! and derivs of that esp wrt r_i)
      ! first that due to ind_dip_d for energy contribution
      ! minus signs arise due to derivs of r_j - r_i wrt r_i
        phi(Ind_000) = Rn_5(Ind_100)*i_dj(1) + Rn_5(Ind_010)*i_dj(2) + &   ! electrostatic potential at position of atm_i from polarizable dipoles of atm_j
                       Rn_5(Ind_001)*i_dj(3)

        phi(Ind_100) = -(Rn_5(Ind_200)*i_dj(1) + Rn_5(Ind_110)*i_dj(2) + &  ! dipole electrostatic field at position of atm_i from polarizable dipoles of atm_j
                         Rn_5(Ind_101)*i_dj(3))

        phi(Ind_010) = -(Rn_5(Ind_110)*i_dj(1) + Rn_5(Ind_020)*i_dj(2) + & ! dipole electrostatic field at position of atm_i from polarizable dipoles of atm_j
                         Rn_5(Ind_011)*i_dj(3))

        phi(Ind_001) = -(Rn_5(Ind_101)*i_dj(1) + Rn_5(Ind_011)*i_dj(2) + & ! dipole electrostatic field at position of atm_i from polarizable dipoles of atm_j
                         Rn_5(Ind_002)*i_dj(3))

        phi(Ind_200) = Rn_5(Ind_300)*i_dj(1) + Rn_5(Ind_210)*i_dj(2) + &  ! quadrupole electrostatic field at position of atm_i from polarizable dipoles of atm_j
                       Rn_5(Ind_201)*i_dj(3)

        phi(Ind_020) = Rn_5(Ind_120)*i_dj(1) + Rn_5(Ind_030)*i_dj(2) + &  ! quadrupole electrostatic field at position of atm_i from polarizable dipoles of atm_j
                       Rn_5(Ind_021)*i_dj(3)

        phi(Ind_002) = Rn_5(Ind_102)*i_dj(1) + Rn_5(Ind_012)*i_dj(2) + &  ! quadrupole electrostatic field at position of atm_i from polarizable dipoles of atm_j
                       Rn_5(Ind_003)*i_dj(3)

        phi(Ind_110) = Rn_5(Ind_210)*i_dj(1) + Rn_5(Ind_120)*i_dj(2) + &  ! quadrupole electrostatic field at position of atm_i from polarizable dipoles of atm_j
                       Rn_5(Ind_111)*i_dj(3)

        phi(Ind_101) = Rn_5(Ind_201)*i_dj(1) + Rn_5(Ind_111)*i_dj(2) + &  ! quadrupole electrostatic field at position of atm_i from polarizable dipoles of atm_j
                       Rn_5(Ind_102)*i_dj(3)

        phi(Ind_011) = Rn_5(Ind_111)*i_dj(1) + Rn_5(Ind_021)*i_dj(2) + &  ! quadrupole electrostatic field at position of atm_i from polarizable dipoles of atm_j
                       Rn_5(Ind_012)*i_dj(3)

        e_add = 0.5d0 * &                                                    ! electrostatic interactions of multipoles of atm_i with polarizable dipoles of atm_j
                (phi(Ind_000)*gmi(Ind_000)+phi(Ind_100)*gmi(Ind_100)+ &
                 phi(Ind_010)*gmi(Ind_010)+phi(Ind_001)*gmi(Ind_001)+ &
                 phi(Ind_200)*gmi(Ind_200)+phi(Ind_020)*gmi(Ind_020)+ &
                 phi(Ind_002)*gmi(Ind_002)+phi(Ind_110)*gmi(Ind_110)+ &
                 phi(Ind_101)*gmi(Ind_101)+phi(Ind_011)*gmi(Ind_011))
        
        e_ind = e_ind + e_add

        ! next that due to ind_dip_d+ind_dip_p for force contribution

        phi(Ind_000) = Rn_5(Ind_100)*i_mj(1) + Rn_5(Ind_010)*i_mj(2) + &
                       Rn_5(Ind_001)*i_mj(3)

        phi(Ind_100) = -(Rn_5(Ind_200)*i_mj(1) + Rn_5(Ind_110)*i_mj(2) + &
                         Rn_5(Ind_101)*i_mj(3))

        phi(Ind_010) = -(Rn_5(Ind_110)*i_mj(1) + Rn_5(Ind_020)*i_mj(2) + &
                         Rn_5(Ind_011)*i_mj(3))

        phi(Ind_001) = -(Rn_5(Ind_101)*i_mj(1) + Rn_5(Ind_011)*i_mj(2) + &
                         Rn_5(Ind_002)*i_mj(3))

        phi(Ind_200) = Rn_5(Ind_300)*i_mj(1) + Rn_5(Ind_210)*i_mj(2) + &
                       Rn_5(Ind_201)*i_mj(3)

        phi(Ind_020) = Rn_5(Ind_120)*i_mj(1) + Rn_5(Ind_030)*i_mj(2) + &
                       Rn_5(Ind_021)*i_mj(3)

        phi(Ind_002) = Rn_5(Ind_102)*i_mj(1) + Rn_5(Ind_012)*i_mj(2) + &
                       Rn_5(Ind_003)*i_mj(3)

        phi(Ind_110) = Rn_5(Ind_210)*i_mj(1) + Rn_5(Ind_120)*i_mj(2) + &
                       Rn_5(Ind_111)*i_mj(3)

        phi(Ind_101) = Rn_5(Ind_201)*i_mj(1) + Rn_5(Ind_111)*i_mj(2) + &
                       Rn_5(Ind_102)*i_mj(3)

        phi(Ind_011) = Rn_5(Ind_111)*i_mj(1) + Rn_5(Ind_021)*i_mj(2) + &
                       Rn_5(Ind_012)*i_mj(3)

        phi(Ind_300) = -(Rn_5(Ind_400)*i_mj(1) + Rn_5(Ind_310)*i_mj(2) + &
                         Rn_5(Ind_301)*i_mj(3))

        phi(Ind_030) = -(Rn_5(Ind_130)*i_mj(1) + Rn_5(Ind_040)*i_mj(2) + &
                         Rn_5(Ind_031)*i_mj(3))

        phi(Ind_003) = -(Rn_5(Ind_103)*i_mj(1) + Rn_5(Ind_013)*i_mj(2) + &
                         Rn_5(Ind_004)*i_mj(3))

        phi(Ind_210) = -(Rn_5(Ind_310)*i_mj(1) + Rn_5(Ind_220)*i_mj(2) + &
                         Rn_5(Ind_211)*i_mj(3))

        phi(Ind_201) = -(Rn_5(Ind_301)*i_mj(1) + Rn_5(Ind_211)*i_mj(2) + &
                         Rn_5(Ind_202)*i_mj(3))

        phi(Ind_120) = -(Rn_5(Ind_220)*i_mj(1) + Rn_5(Ind_130)*i_mj(2) + &
                         Rn_5(Ind_121)*i_mj(3))

        phi(Ind_021) = -(Rn_5(Ind_121)*i_mj(1) + Rn_5(Ind_031)*i_mj(2) + &
                         Rn_5(Ind_022)*i_mj(3))

        phi(Ind_102) = -(Rn_5(Ind_202)*i_mj(1) + Rn_5(Ind_112)*i_mj(2) + &
                         Rn_5(Ind_103)*i_mj(3))

        phi(Ind_012) = -(Rn_5(Ind_112)*i_mj(1) + Rn_5(Ind_022)*i_mj(2) + &
                         Rn_5(Ind_013)*i_mj(3))

        phi(Ind_111) = -(Rn_5(Ind_211)*i_mj(1) + Rn_5(Ind_121)*i_mj(2) + &
                         Rn_5(Ind_112)*i_mj(3))

        g_ind(1) = g_ind(1) + 0.5d0* &
                   (phi(Ind_100)*gmi(Ind_000)+phi(Ind_200)*gmi(Ind_100)+ &
                    phi(Ind_110)*gmi(Ind_010)+phi(Ind_101)*gmi(Ind_001)+ &
                    phi(Ind_300)*gmi(Ind_200)+phi(Ind_120)*gmi(Ind_020)+ &
                    phi(Ind_102)*gmi(Ind_002)+phi(Ind_210)*gmi(Ind_110)+ &
                    phi(Ind_201)*gmi(Ind_101)+phi(Ind_111)*gmi(Ind_011))

        g_ind(2) = g_ind(2) + 0.5d0* &
                   (phi(Ind_010)*gmi(Ind_000)+phi(Ind_110)*gmi(Ind_100)+ &
                    phi(Ind_020)*gmi(Ind_010)+phi(Ind_011)*gmi(Ind_001)+ &
                    phi(Ind_210)*gmi(Ind_200)+phi(Ind_030)*gmi(Ind_020)+ &
                    phi(Ind_012)*gmi(Ind_002)+phi(Ind_120)*gmi(Ind_110)+ &
                    phi(Ind_111)*gmi(Ind_101)+phi(Ind_021)*gmi(Ind_011))

        g_ind(3) = g_ind(3) + 0.5d0* &
                   (phi(Ind_001)*gmi(Ind_000)+phi(Ind_101)*gmi(Ind_100)+ &
                    phi(Ind_011)*gmi(Ind_010)+phi(Ind_002)*gmi(Ind_001)+ &
                    phi(Ind_201)*gmi(Ind_200)+phi(Ind_021)*gmi(Ind_020)+ &
                    phi(Ind_003)*gmi(Ind_002)+phi(Ind_111)*gmi(Ind_110)+ &
                    phi(Ind_102)*gmi(Ind_101)+phi(Ind_012)*gmi(Ind_011))

        ! torque field at i due to induced moments of j
        !  note the factor of 1/2

        do jj = 1, 10
          torque_field(jj, atm_i) = torque_field(jj, atm_i) + 0.5d0 * phi(jj)
        end do

        if (is_polarizable(atm_i)) then ! i, j both polarizable

        ! phi array (electrostatic potential at i due to j induced moments
        ! and derivs of that esp wrt r_i)

          phi(Ind_200) = Rn_5(Ind_300)*i_dj(1) + Rn_5(Ind_210)*i_dj(2) + &
                         Rn_5(Ind_201)*i_dj(3)

          phi(Ind_020) = Rn_5(Ind_120)*i_dj(1) + Rn_5(Ind_030)*i_dj(2) + &
                         Rn_5(Ind_021)*i_dj(3)

          phi(Ind_002) = Rn_5(Ind_102)*i_dj(1) + Rn_5(Ind_012)*i_dj(2) + &
                         Rn_5(Ind_003)*i_dj(3)

          phi(Ind_110) = Rn_5(Ind_210)*i_dj(1) + Rn_5(Ind_120)*i_dj(2) + &
                         Rn_5(Ind_111)*i_dj(3)

          phi(Ind_101) = Rn_5(Ind_201)*i_dj(1) + Rn_5(Ind_111)*i_dj(2) + &
                         Rn_5(Ind_102)*i_dj(3)

          phi(Ind_011) = Rn_5(Ind_111)*i_dj(1) + Rn_5(Ind_021)*i_dj(2) + &
                         Rn_5(Ind_012)*i_dj(3)

          g_ind(1) = g_ind(1) + 0.5d0 * &
                     (phi(Ind_200)*i_pi(1)+phi(Ind_110)*i_pi(2) + &
                      phi(Ind_101)*i_pi(3))

          g_ind(2) = g_ind(2) + 0.5d0 * &
                     (phi(Ind_110)*i_pi(1)+phi(Ind_020)*i_pi(2) + &
                      phi(Ind_011)*i_pi(3))

          g_ind(3) = g_ind(3) + 0.5d0 * &
                     (phi(Ind_101)*i_pi(1)+phi(Ind_011)*i_pi(2) + &
                      phi(Ind_002)*i_pi(3))

          phi(Ind_200) = Rn_5(Ind_300)*i_pj(1) + Rn_5(Ind_210)*i_pj(2) + &
                         Rn_5(Ind_201)*i_pj(3)

          phi(Ind_020) = Rn_5(Ind_120)*i_pj(1) + Rn_5(Ind_030)*i_pj(2) + &
                         Rn_5(Ind_021)*i_pj(3)

          phi(Ind_002) = Rn_5(Ind_102)*i_pj(1) + Rn_5(Ind_012)*i_pj(2) + &
                         Rn_5(Ind_003)*i_pj(3)

          phi(Ind_110) = Rn_5(Ind_210)*i_pj(1) + Rn_5(Ind_120)*i_pj(2) + &
                         Rn_5(Ind_111)*i_pj(3)

          phi(Ind_101) = Rn_5(Ind_201)*i_pj(1) + Rn_5(Ind_111)*i_pj(2) + &
                         Rn_5(Ind_102)*i_pj(3)

          phi(Ind_011) = Rn_5(Ind_111)*i_pj(1) + Rn_5(Ind_021)*i_pj(2) + &
                         Rn_5(Ind_012)*i_pj(3)

          g_ind(1) = g_ind(1) + 0.5d0 * &
                     (phi(Ind_200)*i_di(1)+phi(Ind_110)*i_di(2) + &
                      phi(Ind_101)*i_di(3))

          g_ind(2) = g_ind(2) + 0.5d0 * &
                     (phi(Ind_110)*i_di(1)+phi(Ind_020)*i_di(2) + &
                      phi(Ind_011)*i_di(3))

          g_ind(3) = g_ind(3) + 0.5d0 * &
                     (phi(Ind_101)*i_di(1)+phi(Ind_011)*i_di(2) + &
                      phi(Ind_002)*i_di(3))
 
        end if !is_polarizable(atm_i)
      end if !is_polarizable(atm_j) 
      
      if ((delr2 .lt. ee_damped_cut2) .and. &
           (is_polarizable(atm_i) .or. is_polarizable(atm_j))) then

        !-------------------------------------------------------
        ! McMurchie-Davidson holds for damped tensor as well---in fact,
        ! BD below satisfies BD(n+1) = (1/r)(d/dr)BD(n)
        ! RD(0,0,0,n) = BD(n), n = 1,2,3
        ! RD(t,u,v,n) = (d/dx)^t (d/dy)^u (d/dz)^v R(0,0,0,n)
        !-------------------------------------------------------

        delr3inv = delr2inv / delr
        delr5inv = delr3inv * delr2inv
        delr7inv = delr5inv * delr2inv
        delr9inv = delr7inv * delr2inv
        if( single_thole_val ) then
          expon = asq * delr2 * delr * img(img_j)%qterm
        else      
          expon = screen_polar(atm_i)*screen_polar(atm_j)*delr2*delr* img(img_i)%qterm * img(img_j)%qterm      
!          write(*,*) "am_direct_ene_frc_i  screen_polar(atm_i)*screen_polar(atm_j) = ", screen_polar(atm_i)*screen_polar(atm_j)
        endif  
        expo = exp(-expon)

        ! clam3 = 1.d0-lam3, clam5 = 1.d0-lam5 etc. where 
        ! lam is from ponder's paper

        clam3 = expo
        clam5 = (1.d0 + expon) * expo
        clam7 = (1.d0 + expon + const1 * expon**2) * expo
        clam9 = (1.d0 + expon + const2 * expon**2 + const3 * expon**3) * expo
        BD(1) = -clam3 * delr3inv
        BD(2) = 3.d0 * clam5 * delr5inv
        BD(3) = -15.d0 * clam7 * delr7inv
        BD(4) = 105.d0 * clam9 * delr9inv
        n = 4
        Rn(Ind_000) = BD(n)
        Rn_1(Ind_000) = BD(n-1)
        Rn_1(Ind_100) = delx*Rn(Ind_000)
        Rn_1(Ind_010) = dely*Rn(Ind_000)
        Rn_1(Ind_001) = delz*Rn(Ind_000)
        Rn_2(Ind_000) = BD(n-2)
        Rn_2(Ind_100) = delx*Rn_1(Ind_000)
        Rn_2(Ind_010) = dely*Rn_1(Ind_000)
        Rn_2(Ind_001) = delz*Rn_1(Ind_000)
        Rn_2(Ind_200) = Rn_1(Ind_000) + delx*Rn_1(Ind_100)
        Rn_2(Ind_020) = Rn_1(Ind_000) + dely*Rn_1(Ind_010)
        Rn_2(Ind_002) = Rn_1(Ind_000) + delz*Rn_1(Ind_001)
        Rn_2(Ind_110) = delx*Rn_1(Ind_010)
        Rn_2(Ind_101) = delx*Rn_1(Ind_001)
        Rn_2(Ind_011) = dely*Rn_1(Ind_001)
        Rn_3(Ind_000) = BD(n-3) 
        Rn_3(Ind_100) = delx*Rn_2(Ind_000)
        Rn_3(Ind_010) = dely*Rn_2(Ind_000)
        Rn_3(Ind_001) = delz*Rn_2(Ind_000)
        Rn_3(Ind_200) = Rn_2(Ind_000) + delx*Rn_2(Ind_100)
        Rn_3(Ind_020) = Rn_2(Ind_000) + dely*Rn_2(Ind_010)
        Rn_3(Ind_002) = Rn_2(Ind_000) + delz*Rn_2(Ind_001)
        Rn_3(Ind_110) = delx*Rn_2(Ind_010)
        Rn_3(Ind_101) = delx*Rn_2(Ind_001)
        Rn_3(Ind_011) = dely*Rn_2(Ind_001)
        Rn_3(Ind_300) = 2.d0*Rn_2(Ind_100) + delx*Rn_2(Ind_200)
        Rn_3(Ind_030) = 2.d0*Rn_2(Ind_010) + dely*Rn_2(Ind_020)
        Rn_3(Ind_003) = 2.d0*Rn_2(Ind_001) + delz*Rn_2(Ind_002)
        Rn_3(Ind_210) = dely*Rn_2(Ind_200)
        Rn_3(Ind_201) = delz*Rn_2(Ind_200)
        Rn_3(Ind_120) = delx*Rn_2(Ind_020)
        Rn_3(Ind_021) = delz*Rn_2(Ind_020)
        Rn_3(Ind_102) = delx*Rn_2(Ind_002)
        Rn_3(Ind_012) = dely*Rn_2(Ind_002)
        Rn_3(Ind_111) = delx*Rn_2(Ind_011)
        !Rn_4(Ind_000) = BD(n-4) NOT NEEDED
        Rn_4(Ind_100) = delx*Rn_3(Ind_000)
        Rn_4(Ind_010) = dely*Rn_3(Ind_000)
        Rn_4(Ind_001) = delz*Rn_3(Ind_000)
        Rn_4(Ind_200) = Rn_3(Ind_000) + delx*Rn_3(Ind_100)
        Rn_4(Ind_020) = Rn_3(Ind_000) + dely*Rn_3(Ind_010)
        Rn_4(Ind_002) = Rn_3(Ind_000) + delz*Rn_3(Ind_001)
        Rn_4(Ind_110) = delx*Rn_3(Ind_010)
        Rn_4(Ind_101) = delx*Rn_3(Ind_001)
        Rn_4(Ind_011) = dely*Rn_3(Ind_001)
        Rn_4(Ind_300) = 2.d0*Rn_3(Ind_100) + delx*Rn_3(Ind_200)
        Rn_4(Ind_030) = 2.d0*Rn_3(Ind_010) + dely*Rn_3(Ind_020)
        Rn_4(Ind_003) = 2.d0*Rn_3(Ind_001) + delz*Rn_3(Ind_002)
        Rn_4(Ind_210) = dely*Rn_3(Ind_200)
        Rn_4(Ind_201) = delz*Rn_3(Ind_200)
        Rn_4(Ind_120) = delx*Rn_3(Ind_020)
        Rn_4(Ind_021) = delz*Rn_3(Ind_020)
        Rn_4(Ind_102) = delx*Rn_3(Ind_002)
        Rn_4(Ind_012) = dely*Rn_3(Ind_002)
        Rn_4(Ind_111) = delx*Rn_3(Ind_011)
        Rn_4(Ind_400) = 3.d0*Rn_3(Ind_200) + delx*Rn_3(Ind_300)
        Rn_4(Ind_040) = 3.d0*Rn_3(Ind_020) + dely*Rn_3(Ind_030)
        Rn_4(Ind_004) = 3.d0*Rn_3(Ind_002) + delz*Rn_3(Ind_003)
        Rn_4(Ind_310) = dely*Rn_3(Ind_300)
        Rn_4(Ind_301) = delz*Rn_3(Ind_300)
        Rn_4(Ind_130) = delx*Rn_3(Ind_030)
        Rn_4(Ind_031) = delz*Rn_3(Ind_030)
        Rn_4(Ind_103) = delx*Rn_3(Ind_003)
        Rn_4(Ind_013) = dely*Rn_3(Ind_003)
        Rn_4(Ind_220) = Rn_3(Ind_020) + delx*Rn_3(Ind_120)
        Rn_4(Ind_202) = Rn_3(Ind_002) + delx*Rn_3(Ind_102)
        Rn_4(Ind_022) = Rn_3(Ind_002) + dely*Rn_3(Ind_012)
        Rn_4(Ind_211) = dely*Rn_3(Ind_201)
        Rn_4(Ind_121) = delx*Rn_3(Ind_021)
        Rn_4(Ind_112) = delx*Rn_3(Ind_012)

        if (is_polarizable(atm_i)) then

          ! phi(Ind_000) NOT NEEDED
          !phi(Ind_000)= &
                   !Rn_4(Ind_000)*gmj(Ind_000)+Rn_4(Ind_100)*gmj(Ind_100)+ &
                   !Rn_4(Ind_010)*gmj(Ind_010)+Rn_4(Ind_001)*gmj(Ind_001)+ &
                   !Rn_4(Ind_200)*gmj(Ind_200)+Rn_4(Ind_020)*gmj(Ind_020)+ &
                   !Rn_4(Ind_002)*gmj(Ind_002)+Rn_4(Ind_110)*gmj(Ind_110)+ &
                   !Rn_4(Ind_101)*gmj(Ind_101)+Rn_4(Ind_011)*gmj(Ind_011)

          phi(Ind_100)= &
                 -(Rn_4(Ind_100)*gmj(Ind_000)+Rn_4(Ind_200)*gmj(Ind_100)+ &
                   Rn_4(Ind_110)*gmj(Ind_010)+Rn_4(Ind_101)*gmj(Ind_001)+ &
                   Rn_4(Ind_300)*gmj(Ind_200)+Rn_4(Ind_120)*gmj(Ind_020)+ &
                   Rn_4(Ind_102)*gmj(Ind_002)+Rn_4(Ind_210)*gmj(Ind_110)+ &
                   Rn_4(Ind_201)*gmj(Ind_101)+Rn_4(Ind_111)*gmj(Ind_011))

          phi(Ind_010)= &
                 -(Rn_4(Ind_010)*gmj(Ind_000)+Rn_4(Ind_110)*gmj(Ind_100)+ &
                   Rn_4(Ind_020)*gmj(Ind_010)+Rn_4(Ind_011)*gmj(Ind_001)+ &
                   Rn_4(Ind_210)*gmj(Ind_200)+Rn_4(Ind_030)*gmj(Ind_020)+ &
                   Rn_4(Ind_012)*gmj(Ind_002)+Rn_4(Ind_120)*gmj(Ind_110)+ &
                   Rn_4(Ind_111)*gmj(Ind_101)+Rn_4(Ind_021)*gmj(Ind_011))

          phi(Ind_001)= &
                  -(Rn_4(Ind_001)*gmj(Ind_000)+Rn_4(Ind_101)*gmj(Ind_100)+ &
                   Rn_4(Ind_011)*gmj(Ind_010)+Rn_4(Ind_002)*gmj(Ind_001)+ &
                   Rn_4(Ind_201)*gmj(Ind_200)+Rn_4(Ind_021)*gmj(Ind_020)+ &
                   Rn_4(Ind_003)*gmj(Ind_002)+Rn_4(Ind_111)*gmj(Ind_110)+ &
                   Rn_4(Ind_102)*gmj(Ind_101)+Rn_4(Ind_012)*gmj(Ind_011))

          phi(Ind_200)=  &
                   Rn_4(Ind_200)*gmj(Ind_000)+Rn_4(Ind_300)*gmj(Ind_100)+ &
                   Rn_4(Ind_210)*gmj(Ind_010)+Rn_4(Ind_201)*gmj(Ind_001)+ &
                   Rn_4(Ind_400)*gmj(Ind_200)+Rn_4(Ind_220)*gmj(Ind_020)+ &
                   Rn_4(Ind_202)*gmj(Ind_002)+Rn_4(Ind_310)*gmj(Ind_110)+ &
                   Rn_4(Ind_301)*gmj(Ind_101)+Rn_4(Ind_211)*gmj(Ind_011)

          phi(Ind_020)= &
                   Rn_4(Ind_020)*gmj(Ind_000)+Rn_4(Ind_120)*gmj(Ind_100)+ &
                   Rn_4(Ind_030)*gmj(Ind_010)+Rn_4(Ind_021)*gmj(Ind_001)+ &
                   Rn_4(Ind_220)*gmj(Ind_200)+Rn_4(Ind_040)*gmj(Ind_020)+ &
                   Rn_4(Ind_022)*gmj(Ind_002)+Rn_4(Ind_130)*gmj(Ind_110)+ &
                   Rn_4(Ind_121)*gmj(Ind_101)+Rn_4(Ind_031)*gmj(Ind_011)

          phi(Ind_002)= &
                   Rn_4(Ind_002)*gmj(Ind_000)+Rn_4(Ind_102)*gmj(Ind_100)+ &
                   Rn_4(Ind_012)*gmj(Ind_010)+Rn_4(Ind_003)*gmj(Ind_001)+ &
                   Rn_4(Ind_202)*gmj(Ind_200)+Rn_4(Ind_022)*gmj(Ind_020)+ &
                   Rn_4(Ind_004)*gmj(Ind_002)+Rn_4(Ind_112)*gmj(Ind_110)+ &
                   Rn_4(Ind_103)*gmj(Ind_101)+Rn_4(Ind_013)*gmj(Ind_011)

          phi(Ind_110)= &
                   Rn_4(Ind_110)*gmj(Ind_000)+Rn_4(Ind_210)*gmj(Ind_100)+ &
                   Rn_4(Ind_120)*gmj(Ind_010)+Rn_4(Ind_111)*gmj(Ind_001)+ &
                   Rn_4(Ind_310)*gmj(Ind_200)+Rn_4(Ind_130)*gmj(Ind_020)+ &
                   Rn_4(Ind_112)*gmj(Ind_002)+Rn_4(Ind_220)*gmj(Ind_110)+ &
                   Rn_4(Ind_211)*gmj(Ind_101)+Rn_4(Ind_121)*gmj(Ind_011)

          phi(Ind_101)= &
                   Rn_4(Ind_101)*gmj(Ind_000)+Rn_4(Ind_201)*gmj(Ind_100)+ &
                   Rn_4(Ind_111)*gmj(Ind_010)+Rn_4(Ind_102)*gmj(Ind_001)+ &
                   Rn_4(Ind_301)*gmj(Ind_200)+Rn_4(Ind_121)*gmj(Ind_020)+ &
                   Rn_4(Ind_103)*gmj(Ind_002)+Rn_4(Ind_211)*gmj(Ind_110)+ &
                   Rn_4(Ind_202)*gmj(Ind_101)+Rn_4(Ind_112)*gmj(Ind_011)

          phi(Ind_011)=  &
                   Rn_4(Ind_011)*gmj(Ind_000)+Rn_4(Ind_111)*gmj(Ind_100)+ &
                   Rn_4(Ind_021)*gmj(Ind_010)+Rn_4(Ind_012)*gmj(Ind_001)+ &
                   Rn_4(Ind_211)*gmj(Ind_200)+Rn_4(Ind_031)*gmj(Ind_020)+ &
                   Rn_4(Ind_013)*gmj(Ind_002)+Rn_4(Ind_121)*gmj(Ind_110)+ &
                   Rn_4(Ind_112)*gmj(Ind_101)+Rn_4(Ind_022)*gmj(Ind_011)

          ! minus sign since we remove damped correction

          e_ind = e_ind - 0.5d0 * &
                  (phi(Ind_100) * i_di(1) + phi(Ind_010) * i_di(2) + &
                   phi(Ind_001) * i_di(3))
    
          g_ind(1) = g_ind(1) - 0.5d0 * &
                     (phi(Ind_200) * i_mi(1) + phi(Ind_110) * i_mi(2) + &
                      phi(Ind_101) * i_mi(3))

          g_ind(2) = g_ind(2) - 0.5d0 * &
                     (phi(Ind_110) * i_mi(1) + phi(Ind_020) * i_mi(2) + &
                      phi(Ind_011) * i_mi(3))

          g_ind(3) = g_ind(3) - 0.5d0 * &
                     (phi(Ind_101) * i_mi(1) + phi(Ind_011) * i_mi(2) + &
                      phi(Ind_002) * i_mi(3))

          ! next do torque field at j due to induced at i
          ! potential at j do to dipoles negative due to derivs wrt r_i
          ! higher order are derivs of pot wrt r_j so no sign change

          phi(Ind_000) = -(Rn_4(Ind_100)*i_mi(1)+Rn_4(Ind_010)*i_mi(2) + &
                           Rn_4(Ind_001)*i_mi(3))

          phi(Ind_100) = -(Rn_4(Ind_200)*i_mi(1)+Rn_4(Ind_110)*i_mi(2) + &
                           Rn_4(Ind_101)*i_mi(3))

          phi(Ind_010) = -(Rn_4(Ind_110)*i_mi(1)+Rn_4(Ind_020)*i_mi(2) + &
                           Rn_4(Ind_011)*i_mi(3))

          phi(Ind_001) = -(Rn_4(Ind_101)*i_mi(1)+Rn_4(Ind_011)*i_mi(2) + &
                           Rn_4(Ind_002)*i_mi(3))

          phi(Ind_200) = -(Rn_4(Ind_300)*i_mi(1)+Rn_4(Ind_210)*i_mi(2) + &
                           Rn_4(Ind_201)*i_mi(3))

          phi(Ind_020) = -(Rn_4(Ind_120)*i_mi(1)+Rn_4(Ind_030)*i_mi(2) + &
                           Rn_4(Ind_021)*i_mi(3))

          phi(Ind_002) = -(Rn_4(Ind_102)*i_mi(1)+Rn_4(Ind_012)*i_mi(2) + &
                           Rn_4(Ind_003)*i_mi(3))

          phi(Ind_110) = -(Rn_4(Ind_210)*i_mi(1)+Rn_4(Ind_120)*i_mi(2) + &
                           Rn_4(Ind_111)*i_mi(3))

          phi(Ind_101) = -(Rn_4(Ind_201)*i_mi(1)+Rn_4(Ind_111)*i_mi(2) + &
                           Rn_4(Ind_102)*i_mi(3))

          phi(Ind_011) = -(Rn_4(Ind_111)*i_mi(1)+Rn_4(Ind_021)*i_mi(2) + &
                           Rn_4(Ind_012)*i_mi(3))

          ! torque field at j due to induced moments of i
          !  note the factor of 1/2 (term ~ 1/2 the total induced moment)
          ! the minus sign is since we remove damped contributions

          do jj = 1, 10
            torque_field(jj, atm_j) = torque_field(jj, atm_j) - 0.5d0 * phi(jj)
          end do

        end if !is_polarizable(atm_i)

        if (is_polarizable(atm_j)) then

          ! phi array (electrostatic potential at i due to j induced moments
          ! and derivs of that esp wrt r_i)
          ! first that due to ind_dip_d for energy contribution
          ! minus signs arise due to derivs of r_j - r_i wrt r_i

          phi(Ind_000) = Rn_4(Ind_100)*i_dj(1)+Rn_4(Ind_010)*i_dj(2) + &
                         Rn_4(Ind_001)*i_dj(3)

          phi(Ind_100) = -(Rn_4(Ind_200)*i_dj(1)+Rn_4(Ind_110)*i_dj(2) + &
                           Rn_4(Ind_101)*i_dj(3))

          phi(Ind_010) = -(Rn_4(Ind_110)*i_dj(1)+Rn_4(Ind_020)*i_dj(2) + &
                           Rn_4(Ind_011)*i_dj(3))

          phi(Ind_001) = -(Rn_4(Ind_101)*i_dj(1)+Rn_4(Ind_011)*i_dj(2) + &
                           Rn_4(Ind_002)*i_dj(3))

          phi(Ind_200) = Rn_4(Ind_300)*i_dj(1)+Rn_4(Ind_210)*i_dj(2) + &
                         Rn_4(Ind_201)*i_dj(3)

          phi(Ind_020) = Rn_4(Ind_120)*i_dj(1)+Rn_4(Ind_030)*i_dj(2) + &
                         Rn_4(Ind_021)*i_dj(3)

          phi(Ind_002) = Rn_4(Ind_102)*i_dj(1)+Rn_4(Ind_012)*i_dj(2) + &
                         Rn_4(Ind_003)*i_dj(3)

          phi(Ind_110) = Rn_4(Ind_210)*i_dj(1)+Rn_4(Ind_120)*i_dj(2) + &
                         Rn_4(Ind_111)*i_dj(3)

          phi(Ind_101) = Rn_4(Ind_201)*i_dj(1)+Rn_4(Ind_111)*i_dj(2) + &
                         Rn_4(Ind_102)*i_dj(3)

          phi(Ind_011) = Rn_4(Ind_111)*i_dj(1)+Rn_4(Ind_021)*i_dj(2) + &
                         Rn_4(Ind_012)*i_dj(3)

          ! minus sign since we remove damped contributions

          e_add = 0.5d0 * &
                  (phi(Ind_000)*gmi(Ind_000)+phi(Ind_100)*gmi(Ind_100)+ &
                   phi(Ind_010)*gmi(Ind_010)+phi(Ind_001)*gmi(Ind_001)+ &
                   phi(Ind_200)*gmi(Ind_200)+phi(Ind_020)*gmi(Ind_020)+ &
                   phi(Ind_002)*gmi(Ind_002)+phi(Ind_110)*gmi(Ind_110)+ &
                   phi(Ind_101)*gmi(Ind_101)+phi(Ind_011)*gmi(Ind_011))
          
          e_ind = e_ind - e_add

          ! next that due to ind_dip_d+ind_dip_p for force contribution

          phi(Ind_000) = Rn_4(Ind_100)*i_mj(1)+Rn_4(Ind_010)*i_mj(2) + &
                         Rn_4(Ind_001)*i_mj(3)

          phi(Ind_100) = -(Rn_4(Ind_200)*i_mj(1)+Rn_4(Ind_110)*i_mj(2) + &
                           Rn_4(Ind_101)*i_mj(3))

          phi(Ind_010) = -(Rn_4(Ind_110)*i_mj(1)+Rn_4(Ind_020)*i_mj(2) + &
                           Rn_4(Ind_011)*i_mj(3))

          phi(Ind_001) = -(Rn_4(Ind_101)*i_mj(1)+Rn_4(Ind_011)*i_mj(2) + &
                           Rn_4(Ind_002)*i_mj(3))

          phi(Ind_200) = Rn_4(Ind_300)*i_mj(1)+Rn_4(Ind_210)*i_mj(2) + &
                         Rn_4(Ind_201)*i_mj(3)

          phi(Ind_020) = Rn_4(Ind_120)*i_mj(1)+Rn_4(Ind_030)*i_mj(2) + &
                         Rn_4(Ind_021)*i_mj(3)

          phi(Ind_002) = Rn_4(Ind_102)*i_mj(1)+Rn_4(Ind_012)*i_mj(2) + &
                         Rn_4(Ind_003)*i_mj(3)

          phi(Ind_110) = Rn_4(Ind_210)*i_mj(1)+Rn_4(Ind_120)*i_mj(2) + &
                         Rn_4(Ind_111)*i_mj(3)

          phi(Ind_101) = Rn_4(Ind_201)*i_mj(1)+Rn_4(Ind_111)*i_mj(2) + &
                         Rn_4(Ind_102)*i_mj(3)

          phi(Ind_011) = Rn_4(Ind_111)*i_mj(1)+Rn_4(Ind_021)*i_mj(2) + &
                         Rn_4(Ind_012)*i_mj(3)

          phi(Ind_300) = -(Rn_4(Ind_400)*i_mj(1)+Rn_4(Ind_310)*i_mj(2) + &
                           Rn_4(Ind_301)*i_mj(3))

          phi(Ind_030) = -(Rn_4(Ind_130)*i_mj(1)+Rn_4(Ind_040)*i_mj(2) + &
                           Rn_4(Ind_031)*i_mj(3))

          phi(Ind_003) = -(Rn_4(Ind_103)*i_mj(1)+Rn_4(Ind_013)*i_mj(2) + &
                           Rn_4(Ind_004)*i_mj(3))

          phi(Ind_210) = -(Rn_4(Ind_310)*i_mj(1)+Rn_4(Ind_220)*i_mj(2) + &
                           Rn_4(Ind_211)*i_mj(3))

          phi(Ind_201) = -(Rn_4(Ind_301)*i_mj(1)+Rn_4(Ind_211)*i_mj(2) + &
                           Rn_4(Ind_202)*i_mj(3))

          phi(Ind_120) = -(Rn_4(Ind_220)*i_mj(1)+Rn_4(Ind_130)*i_mj(2) + &
                           Rn_4(Ind_121)*i_mj(3))

          phi(Ind_021) = -(Rn_4(Ind_121)*i_mj(1)+Rn_4(Ind_031)*i_mj(2) + &
                           Rn_4(Ind_022)*i_mj(3))

          phi(Ind_102) = -(Rn_4(Ind_202)*i_mj(1)+Rn_4(Ind_112)*i_mj(2) + &
                           Rn_4(Ind_103)*i_mj(3))

          phi(Ind_012) = -(Rn_4(Ind_112)*i_mj(1)+Rn_4(Ind_022)*i_mj(2) + &
                           Rn_4(Ind_013)*i_mj(3))

          phi(Ind_111) = -(Rn_4(Ind_211)*i_mj(1)+Rn_4(Ind_121)*i_mj(2) + &
                           Rn_4(Ind_112)*i_mj(3))

          ! minus sign since we remove damped contributions

          g_ind(1) = g_ind(1) - 0.5d0 * &
                     (phi(Ind_100)*gmi(Ind_000)+phi(Ind_200)*gmi(Ind_100)+ &
                      phi(Ind_110)*gmi(Ind_010)+phi(Ind_101)*gmi(Ind_001)+ &
                      phi(Ind_300)*gmi(Ind_200)+phi(Ind_120)*gmi(Ind_020)+ &
                      phi(Ind_102)*gmi(Ind_002)+phi(Ind_210)*gmi(Ind_110)+ &
                      phi(Ind_201)*gmi(Ind_101)+phi(Ind_111)*gmi(Ind_011))

          g_ind(2) = g_ind(2) - 0.5d0 * &
                     (phi(Ind_010)*gmi(Ind_000)+phi(Ind_110)*gmi(Ind_100)+ &
                      phi(Ind_020)*gmi(Ind_010)+phi(Ind_011)*gmi(Ind_001)+ &
                      phi(Ind_210)*gmi(Ind_200)+phi(Ind_030)*gmi(Ind_020)+ &
                      phi(Ind_012)*gmi(Ind_002)+phi(Ind_120)*gmi(Ind_110)+ &
                      phi(Ind_111)*gmi(Ind_101)+phi(Ind_021)*gmi(Ind_011))

          g_ind(3) = g_ind(3) - 0.5d0 * &
                     (phi(Ind_001)*gmi(Ind_000)+phi(Ind_101)*gmi(Ind_100)+ &
                      phi(Ind_011)*gmi(Ind_010)+phi(Ind_002)*gmi(Ind_001)+ &
                      phi(Ind_201)*gmi(Ind_200)+phi(Ind_021)*gmi(Ind_020)+ &
                      phi(Ind_003)*gmi(Ind_002)+phi(Ind_111)*gmi(Ind_110)+ &
                      phi(Ind_102)*gmi(Ind_101)+phi(Ind_012)*gmi(Ind_011))

          ! torque field at i due to induced moments of j
          !  note the factor of 1/2 (term ~ 1/2 the total induced moment)
          ! the minus sign is since we remove damped contributions

          do jj = 1, 10
            torque_field(jj, atm_i) = torque_field(jj, atm_i) - 0.5d0 * phi(jj)
          end do

          if (is_polarizable(atm_i)) then ! i, j both polarizable

            ! phi array (electrostatic pot at i due to j induced moments
            ! and derivs of that esp wrt r_i)

            phi(Ind_200) = Rn_4(Ind_300)*i_dj(1)+Rn_4(Ind_210)*i_dj(2)+ &
                           Rn_4(Ind_201)*i_dj(3)

            phi(Ind_020) = Rn_4(Ind_120)*i_dj(1)+Rn_4(Ind_030)*i_dj(2)+ &
                           Rn_4(Ind_021)*i_dj(3)

            phi(Ind_002) = Rn_4(Ind_102)*i_dj(1)+Rn_4(Ind_012)*i_dj(2)+ &
                           Rn_4(Ind_003)*i_dj(3)

            phi(Ind_110) = Rn_4(Ind_210)*i_dj(1)+Rn_4(Ind_120)*i_dj(2)+ &
                           Rn_4(Ind_111)*i_dj(3)

            phi(Ind_101) = Rn_4(Ind_201)*i_dj(1)+Rn_4(Ind_111)*i_dj(2)+ &
                           Rn_4(Ind_102)*i_dj(3)

            phi(Ind_011) = Rn_4(Ind_111)*i_dj(1)+Rn_4(Ind_021)*i_dj(2)+ &
                           Rn_4(Ind_012)*i_dj(3)

            ! minus sign since we remove damped contributions

            g_ind(1) = g_ind(1) - 0.5d0 * &
                       (phi(Ind_200)*i_pi(1)+phi(Ind_110)*i_pi(2) + &
                        phi(Ind_101)*i_pi(3))

            g_ind(2) = g_ind(2) - 0.5d0 * &
                       (phi(Ind_110)*i_pi(1)+phi(Ind_020)*i_pi(2) + &
                        phi(Ind_011)*i_pi(3))

            g_ind(3) = g_ind(3) - 0.5d0 * &
                       (phi(Ind_101)*i_pi(1)+phi(Ind_011)*i_pi(2) + &
                        phi(Ind_002)*i_pi(3))

            phi(Ind_200) = Rn_4(Ind_300)*i_pj(1)+Rn_4(Ind_210)*i_pj(2)+ &
                           Rn_4(Ind_201)*i_pj(3)

            phi(Ind_020) = Rn_4(Ind_120)*i_pj(1)+Rn_4(Ind_030)*i_pj(2)+ &
                           Rn_4(Ind_021)*i_pj(3)

            phi(Ind_002) = Rn_4(Ind_102)*i_pj(1)+Rn_4(Ind_012)*i_pj(2)+ &
                           Rn_4(Ind_003)*i_pj(3)

            phi(Ind_110) = Rn_4(Ind_210)*i_pj(1)+Rn_4(Ind_120)*i_pj(2)+ &
                           Rn_4(Ind_111)*i_pj(3)

            phi(Ind_101) = Rn_4(Ind_201)*i_pj(1)+Rn_4(Ind_111)*i_pj(2)+ &
                           Rn_4(Ind_102)*i_pj(3)

            phi(Ind_011) = Rn_4(Ind_111)*i_pj(1)+Rn_4(Ind_021)*i_pj(2)+ &
                           Rn_4(Ind_012)*i_pj(3)

            ! minus sign since we remove damped contributions

            g_ind(1) = g_ind(1) - 0.5d0 * &
                       (phi(Ind_200)*i_di(1)+phi(Ind_110)*i_di(2) + &
                        phi(Ind_101)*i_di(3))

            g_ind(2) = g_ind(2) - 0.5d0 * &
                       (phi(Ind_110)*i_di(1)+phi(Ind_020)*i_di(2) + &
                        phi(Ind_011)*i_di(3))

            g_ind(3) = g_ind(3) - 0.5d0 * &
                       (phi(Ind_101)*i_di(1)+phi(Ind_011)*i_di(2) + &
                        phi(Ind_002)*i_di(3))

          end if !is_polarizable(atm_i)
        end if !is_polarizable(atm_j) 
      end if !(delr2 < ee_damped_cut2) .and. &
             !(is_polarizable(atm_i) .or. is_polarizable(atm_j))

      ! img_frc is negative gradient

      ene_perm = ene_perm + coulomb_const_kcal_per_mole * e_pp
      ene_ind  = ene_ind + coulomb_const_kcal_per_mole * e_ind

      img_frc(1, img_j) = &
        img_frc(1, img_j) + coulomb_const_kcal_per_mole * (g_pp(1) + g_ind(1))
      img_frc(2, img_j) = &
        img_frc(2, img_j) + coulomb_const_kcal_per_mole * (g_pp(2) + g_ind(2))
      img_frc(3, img_j) = &
        img_frc(3, img_j) + coulomb_const_kcal_per_mole * (g_pp(3) + g_ind(3))
      img_frc(1, img_i) = &
        img_frc(1, img_i) - coulomb_const_kcal_per_mole * (g_pp(1) + g_ind(1))
      img_frc(2, img_i) = &
        img_frc(2, img_i) - coulomb_const_kcal_per_mole * (g_pp(2) + g_ind(2))
      img_frc(3, img_i) = &
        img_frc(3, img_i) - coulomb_const_kcal_per_mole * (g_pp(3) + g_ind(3))

      vxx = vxx - delx * (g_pp(1) + g_ind(1))
      vxy = vxy - delx * (g_pp(2) + g_ind(2))
      vxz = vxz - delx * (g_pp(3) + g_ind(3))
      vyx = vyx - dely * (g_pp(1) + g_ind(1))
      vyy = vyy - dely * (g_pp(2) + g_ind(2))
      vyz = vyz - dely * (g_pp(3) + g_ind(3))
      vzx = vzx - delz * (g_pp(1) + g_ind(1))
      vzy = vzy - delz * (g_pp(2) + g_ind(2))
      vzz = vzz - delz * (g_pp(3) + g_ind(3))

    end do

    sublst_head = sublst_head + vec_max

  end do

  vxx = vxx * coulomb_const_kcal_per_mole
  vxy = vxy * coulomb_const_kcal_per_mole
  vxz = vxz * coulomb_const_kcal_per_mole
  vyx = vyx * coulomb_const_kcal_per_mole
  vyy = vyy * coulomb_const_kcal_per_mole
  vyz = vyz * coulomb_const_kcal_per_mole
  vzx = vzx * coulomb_const_kcal_per_mole
  vzy = vzy * coulomb_const_kcal_per_mole
  vzz = vzz * coulomb_const_kcal_per_mole

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

end subroutine am_direct_ene_frc_i

!*******************************************************************************!
! Subroutine:  am_vdw_direct_ene_frc_i
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_vdw_direct_ene_frc_i(atm_i, img, ipairs_sublst, x_tran, &
                                   pair_cnt, crd, ene_vdw, frc, img_frc, &
                                   virial, img_atm_map)

  use mdin_amoeba_dat_mod, only : do_vdw_taper
  use amoeba_flags_mod
  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: atm_i 
  type(img_rec), intent(in)             :: img(*)
  integer, intent(in)                   :: ipairs_sublst(*)
  double precision                      :: x_tran(1:3, 0:17)
  integer, intent(in)                   :: pair_cnt
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: ene_vdw
  double precision, intent(in out)      :: frc(3, *)
  double precision, intent(in out)      :: img_frc(3, *)
  double precision, intent(in out)      :: virial(3, 3)
  integer, intent(in)                   :: img_atm_map(*)

! Local variables:

  integer               :: itran, it, jt, ih, jh, idx
  integer               :: sublst_idx
  integer               :: atm_j, img_j
  double precision      :: wi, wj
  double precision      :: delx, dely, delz
  double precision      :: delr, delr2
  double precision      :: eps
  double precision      :: rad
  double precision      :: rho
  double precision      :: t1, t2
  double precision      :: dt1drho, dt2drho, drhodr
  double precision      :: term
  double precision      :: dfx, dfy, dfz
  double precision      :: rho6, rho7
  double precision      :: vxx, vxy, vxz, vyy, vyz, vzz
  double precision      :: switch, dswitch_dr
  double precision      :: f, dfdr
  double precision      :: delr3, delr4, delr5
  integer, parameter    :: mask27 = Z"07FFFFFF"

  if (do_amoeba_vdw_flag .ne. proceed) return

  ih = vdw_atom_parent(atm_i)
  wi = vdw_atom_parent_crd_wt(atm_i)

  ! We have to set up the translation vectors again due to weighting...

  x_i = x_i + wi * (crd(1, ih) - crd(1, atm_i))
  y_i = y_i + wi * (crd(2, ih) - crd(2, atm_i))
  z_i = z_i + wi * (crd(3, ih) - crd(3, atm_i))

#ifdef DIRFRC_COMTRANS
  if (common_tran .eq. 0) then
#endif /* DIRFRC_COMTRANS */
    ! We need all the translation vectors:
    do idx = 0, 17
      x_tran(1, idx) = tranvec(1, idx) - x_i
      x_tran(2, idx) = tranvec(2, idx) - y_i
      x_tran(3, idx) = tranvec(3, idx) - z_i
    end do
#ifdef DIRFRC_COMTRANS
  else
    ! Just put the x,y,z values in the middle cell
    x_tran(1, 13) = - x_i
    x_tran(2, 13) = - y_i
    x_tran(3, 13) = - z_i
  end if
#endif /* DIRFRC_COMTRANS */

  it = vdw_atom_type(atm_i)

  vxx = 0.d0
  vxy = 0.d0
  vxz = 0.d0
  vyy = 0.d0
  vyz = 0.d0
  vzz = 0.d0

  do sublst_idx = 1, pair_cnt

#ifdef DIRFRC_COMTRANS
    if (common_tran .eq. 1) then
      img_j = ipairs_sublst(sublst_idx)
      itran = 13
    else
#endif /* DIRFRC_COMTRANS */
      img_j = iand(ipairs_sublst(sublst_idx), mask27)
      itran = ishft(ipairs_sublst(sublst_idx), -27)
#ifdef DIRFRC_COMTRANS
    end if
#endif /* DIRFRC_COMTRANS */

    atm_j = img_atm_map(img_j)
    jh = vdw_atom_parent(atm_j)
    wj = vdw_atom_parent_crd_wt(atm_j)

    ! Calculate delx etc between vdw points, not the atoms.

    delx = wj * (crd(1, jh) - crd(1, atm_j)) + img(img_j)%x + x_tran(1, itran)
    dely = wj * (crd(2, jh) - crd(2, atm_j)) + img(img_j)%y + x_tran(2, itran)
    delz = wj * (crd(3, jh) - crd(3, atm_j)) + img(img_j)%z + x_tran(3, itran)

    delr2 = delx * delx + dely * dely + delz * delz

    if (delr2 .lt. vdw_switch_off_2) then

      jt = vdw_atom_type(atm_j)
      eps = vdw_epsilon(jt, it)
      
!      eps = 0.0d0 ! IGOR TMP DEBUG
      
      rad = vdw_radius(jt, it)
      delr = sqrt(delr2)
      rho = delr / rad
      rho6 = rho**6
      rho7 = rho6 * rho
      t1 = ((1.d0 + vdw_buf_delta) / (rho + vdw_buf_delta))**7
      t2 = (1.d0 + vdw_buf_gamma) / (rho7 + vdw_buf_gamma)
      dt1drho = -7.d0 * t1 / (rho + vdw_buf_delta)
      dt2drho = -7.d0 * t2 * (rho6 / (rho7 + vdw_buf_gamma))
      drhodr = 1.d0 / rad
      f = eps * t1 * (t2 - 2.d0)
      dfdr = eps * (dt1drho * (t2 - 2.d0) + t1 * dt2drho) * drhodr

      if (do_vdw_taper .eq. 1 .and. delr2 .gt. vdw_switch_on_2) then
        delr3 = delr2 * delr
        delr4 = delr3 * delr
        delr5 = delr4 * delr
        switch = c5 * delr5 + c4 * delr4 + c3 * delr3 + c2 * delr2 + &
                 c1 * delr + c0
        dswitch_dr = 5.d0 * c5 * delr4 + 4.d0 * c4 * delr3 + &
                     3.d0 * c3 * delr2 + 2.d0 * c2 * delr + c1
        dfdr = switch * dfdr + f * dswitch_dr
        f = switch * f
      end if

      ene_vdw = ene_vdw + f
      term = dfdr / delr
      dfx = term * delx
      dfy = term * dely
      dfz = term * delz

      ! Recall ddelx_dxi = -1.

      img_frc(1, img_i) = img_frc(1, img_i) + (1.d0 - wi) * dfx
      img_frc(2, img_i) = img_frc(2, img_i) + (1.d0 - wi) * dfy
      img_frc(3, img_i) = img_frc(3, img_i) + (1.d0 - wi) * dfz
      frc(1, ih) = frc(1, ih) + wi * dfx
      frc(2, ih) = frc(2, ih) + wi * dfy
      frc(3, ih) = frc(3, ih) + wi * dfz
      img_frc(1, img_j) = img_frc(1, img_j) - (1.d0 - wj) * dfx
      img_frc(2, img_j) = img_frc(2, img_j) - (1.d0 - wj) * dfy
      img_frc(3, img_j) = img_frc(3, img_j) - (1.d0 - wj) * dfz
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

  end do

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

end subroutine am_vdw_direct_ene_frc_i

end subroutine am_direct_ene_frc

!*******************************************************************************!
! Subroutine:  am_direct_count_num_ee_pairs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_direct_count_num_ee_pairs(ipairs, img, tranvec, ee_dsum_cut, &
                                        num_pairs_in_ee_cut)

  use img_mod

  implicit none

! Formal arguments:

  integer                       :: ipairs(*)
  type(img_rec), intent(in)     :: img(*)
  double precision, intent(in)  :: tranvec(1:3, 0:17)
  double precision, intent(in)  :: ee_dsum_cut
  integer, intent(out)          :: num_pairs_in_ee_cut

! Local variables:

  double precision      x_i, y_i, z_i
  double precision      x_tran(1:3, 0:17)
  double precision      ee_dsum_cut2
  integer               i
  integer               ipairs_idx
  integer               img_i
  integer               pair_cnt
#ifdef DIRFRC_COMTRANS
  integer               common_tran    ! flag - 1 if translation not needed
#endif /* DIRFRC_COMTRANS */

  ee_dsum_cut2 = ee_dsum_cut * ee_dsum_cut

  ipairs_idx = 1

  do img_i = my_img_lo, my_img_hi

#ifdef DIRFRC_COMTRANS
    ! Common translation (ie. no translation) flag is packed at
    ! the front of each sublist followed by the count(s) of sublist
    ! image pair entries.

    common_tran = ipairs(ipairs_idx)
    ipairs_idx = ipairs_idx + 1
#endif /* DIRFRC_COMTRANS */
    
    ! Electrostatic evaluation-only count followed by
    ! full evaluation count packed at the front of each pair sublist.

    pair_cnt = ipairs(ipairs_idx) + ipairs(ipairs_idx + 1)
    ipairs_idx = ipairs_idx + 2

    if (pair_cnt .gt. 0) then

      x_i = img(img_i)%x
      y_i = img(img_i)%y
      z_i = img(img_i)%z

#ifdef DIRFRC_COMTRANS
      if (common_tran .eq. 0) then
#endif /* DIRFRC_COMTRANS */
        ! We need all the translation vectors:
        do i = 0, 17
          x_tran(1, i) = tranvec(1, i) - x_i
          x_tran(2, i) = tranvec(2, i) - y_i
          x_tran(3, i) = tranvec(3, i) - z_i
        end do
#ifdef DIRFRC_COMTRANS
      else
        ! Just put the x,y,z values in the middle cell
        x_tran(1, 13) = - x_i
        x_tran(2, 13) = - y_i
        x_tran(3, 13) = - z_i
      end if
#endif /* DIRFRC_COMTRANS */

      call am_direct_increment_ee_pairs(img, ipairs(ipairs_idx), x_tran, &
                                        ee_dsum_cut2, pair_cnt, &
                                        num_pairs_in_ee_cut)

      ipairs_idx = ipairs_idx + pair_cnt

    end if
  end do

  return

contains

!*******************************************************************************!
! Subroutine:  am_direct_increment_ee_pairs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_direct_increment_ee_pairs(img, ipairs_sublst, x_tran, &
                                        ee_dsum_cut2, pair_cnt, &
                                        num_pairs_in_ee_cut)

  implicit none

! Formal arguments:

  type(img_rec), intent(in)     :: img(*)
  integer                       :: ipairs_sublst(*)
  double precision, intent(in)  :: x_tran(1:3, 0:17)
  double precision, intent(in)  :: ee_dsum_cut2
  integer, intent(in)           :: pair_cnt
  integer, intent(in out)       :: num_pairs_in_ee_cut

! Local variables:

  integer                       :: nxt_img_j, img_j
  integer                       :: itran
  integer                       :: sublst_idx
  integer                       :: saved_pairlist_val
  double precision              :: nxt_delx, nxt_dely, nxt_delz
  double precision              :: delx, dely, delz, delr2
  integer, parameter            :: mask27 = Z"07FFFFFF"

  ! The pairlist must have one dummy end entry to cover reading past the
  ! end of the list...

  saved_pairlist_val = ipairs_sublst(pair_cnt + 1)

  ipairs_sublst(pair_cnt + 1) = ipairs_sublst(pair_cnt)

#ifdef DIRFRC_COMTRANS
  if (common_tran .eq. 1) then
    nxt_img_j = ipairs_sublst(1)
    itran = 13
  else
#endif /* DIRFRC_COMTRANS */
    nxt_img_j = iand(ipairs_sublst(1), mask27)
    itran = ishft(ipairs_sublst(1), -27)
#ifdef DIRFRC_COMTRANS
  end if
#endif /* DIRFRC_COMTRANS */

  nxt_delx = img(nxt_img_j)%x + x_tran(1, itran)
  nxt_dely = img(nxt_img_j)%y + x_tran(2, itran)
  nxt_delz = img(nxt_img_j)%z + x_tran(3, itran)

  do sublst_idx = 2, pair_cnt + 1

    img_j = nxt_img_j
    delx = nxt_delx
    dely = nxt_dely
    delz = nxt_delz

#ifdef DIRFRC_COMTRANS
    if (common_tran .eq. 1) then
      nxt_img_j = ipairs_sublst(sublst_idx)
      itran = 13
    else
#endif /* DIRFRC_COMTRANS */
      nxt_img_j = iand(ipairs_sublst(sublst_idx), mask27)
      itran = ishft(ipairs_sublst(sublst_idx), -27)
#ifdef DIRFRC_COMTRANS
    end if
#endif /* DIRFRC_COMTRANS */

    nxt_delx = img(nxt_img_j)%x + x_tran(1, itran)
    nxt_dely = img(nxt_img_j)%y + x_tran(2, itran)
    nxt_delz = img(nxt_img_j)%z + x_tran(3, itran)

    delr2 = delx * delx + dely * dely + delz * delz

    if (delr2 .lt. ee_dsum_cut2) num_pairs_in_ee_cut = num_pairs_in_ee_cut + 1

  end do

  ipairs_sublst(pair_cnt + 1) = saved_pairlist_val

  return

end subroutine am_direct_increment_ee_pairs

end subroutine am_direct_count_num_ee_pairs

!*******************************************************************************!
! Subroutine:  am_direct_dip_dip_field
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_direct_dip_dip_field(ind_dip_d, ind_dip_p, dip_field_d, &
                                   dip_field_p)

  use amoeba_flags_mod
  use timers_mod

  implicit none

! Formal arguments:

  double precision, intent(in)          :: ind_dip_d(3, *)
  double precision, intent(in)          :: ind_dip_p(3, *)
  double precision, intent(in out)      :: dip_field_d(3, *)
  double precision, intent(in out)      :: dip_field_p(3, *)

  if (iand(do_amoeba_direct_flag, proceed_induce) .ne. proceed_induce) return

  call am_direct_calc_dipdip_field(num_tensor, dipole_dipole_list, &
                                   dipole_dipole_tensor, ind_dip_d, &
                                   ind_dip_p, dip_field_d, dip_field_p)

  call update_pme_time(dir_frc_sum_timer)

  return

end subroutine am_direct_dip_dip_field

!*******************************************************************************!
! Subroutine:  am_direct_calc_dipdip_field
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_direct_calc_dipdip_field(num_tensor, dipole_dipole_list, &
                                       dipole_dipole_tensor, ind_dip_d, &
                                       ind_dip_p, dip_field_d, dip_field_p)

  implicit none

! Formal arguments:

  integer, intent(in)                   :: num_tensor
  integer, intent(in)                   :: dipole_dipole_list(2, *)
  double precision, intent(in)          :: dipole_dipole_tensor(6, *)
  double precision, intent(in)          :: ind_dip_d(3, *)
  double precision, intent(in)          :: ind_dip_p(3, *)
  double precision, intent(in out)      :: dip_field_d(3, *)
  double precision, intent(in out)      :: dip_field_p(3, *)

! Local variables:

  integer                               :: i, j, n

  ! UNDER MPI ONLY tensors for which you own img_i (img of atm_i) are
  ! stored!
  do n = 1, num_tensor
    i = dipole_dipole_list(1, n)
    j = dipole_dipole_list(2, n)

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

    ! other set of dipoles, fields

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

end subroutine am_direct_calc_dipdip_field

end module amoeba_direct_mod
