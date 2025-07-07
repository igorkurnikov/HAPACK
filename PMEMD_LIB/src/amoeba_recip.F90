#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_recip_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_recip_mod

  use pme_fft_mod
  use parallel_dat_mod
  use pme_recip_mod

  implicit none

  private

  integer, save ::   do_amoeba_recip_flag

  integer, save                         :: amba_fft_alloc_size
  logical, save                         :: perm_field_done = .false.
  integer, save                         :: my_chg_cnt
  integer, parameter                    :: dr_order = 3

  ! This storage should be allocated on entry to and deallocated on exit
  ! from a nonbond evaluation step.  This will allow for resizing the data
  ! as needed.

  double precision, allocatable, save   :: gbl_G_func(:)
  double precision, allocatable, save   :: gbl_Qperm(:)
  double precision, allocatable, save   :: gbl_Q1(:)  ! ind_dip_d spread on the grid using splines 
  double precision, allocatable, save   :: gbl_Q2(:)  ! ind_dip_p spread on the grid using splines 
  double precision, allocatable, save   :: gbl_theta1(:,:,:)  ! spline coef of multipoles 1st dir ??? indexed by deriv, grid spline idx and charge num
  double precision, allocatable, save   :: gbl_theta2(:,:,:)  ! spline coef of multipoles 2nd dir ??? 
  double precision, allocatable, save   :: gbl_theta3(:,:,:)  ! spline coef of multipoles 3rd dir ??? 
  integer, allocatable, save            :: gbl_init_grid_ind(:,:) ! initial grid indexes for spline coef of multipole expansions 

  ! BUGBUG - Must eventually img-base these:
  double precision, allocatable, save   :: gbl_fractional_multipole(:,:)
  double precision, allocatable, save   :: gbl_perm_F_field(:,:)
  integer, allocatable, save            :: gbl_my_chgs(:)

  public gbl_fractional_multipole
  public gbl_perm_F_field
  
#include "amoeba_mpole_index.i"

  public        am_recip_zero_flag
  public        am_recip_set_user_bit
  public        am_recip_perm_field
  public        am_recip_dipole_field
  public        am_recip_ene_frc
  public        amoeba_recip_final_setup
  public        amoeba_recip_step_alloc_dat
  public        amoeba_recip_step_dealloc_dat
         
contains

!*******************************************************************************!
! Subroutine:  am_recip_zero_flag
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_zero_flag
  implicit none
  do_amoeba_recip_flag = 0
  return
end subroutine am_recip_zero_flag

!*******************************************************************************!
! Subroutine:  am_recip_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_set_user_bit(do_this)

  use amoeba_flags_mod
  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in) :: do_this

  if (do_this .eq. 1) then      ! do in all cases
    do_amoeba_recip_flag = ibset(do_amoeba_recip_flag, user_induce_bit)
    do_amoeba_recip_flag = ibset(do_amoeba_recip_flag, user_postinduce_bit)
  else if (do_this .eq. 2) then ! do the induction, not the post-induction
    do_amoeba_recip_flag = ibset(do_amoeba_recip_flag, user_induce_bit)
    do_amoeba_recip_flag = ibclr(do_amoeba_recip_flag, user_postinduce_bit)
  else if (do_this .eq. 3) then ! do the post-induction, not the induction
    do_amoeba_recip_flag = ibclr(do_amoeba_recip_flag, user_induce_bit)
    do_amoeba_recip_flag = ibset(do_amoeba_recip_flag, user_postinduce_bit)
  else if (do_this .eq. 0) then 
    do_amoeba_recip_flag = ibclr(do_amoeba_recip_flag, user_induce_bit)
    do_amoeba_recip_flag = ibclr(do_amoeba_recip_flag, user_postinduce_bit)
  else
    error_msg =  'am_recip_set_user_bit: bad value of user do_this'
    call mol_mech_error
  end if

  return

end subroutine am_recip_set_user_bit

!*******************************************************************************!
! Subroutine:  amoeba_recip_final_setup
!
! Description: <TBS>
!
!*******************************************************************************

subroutine amoeba_recip_final_setup

  use amoeba_flags_mod
  use file_io_dat_mod
  use parallel_dat_mod
  use mdin_ewald_dat_mod
  use prmtop_dat_mod

  implicit none

! Local variables:

  integer, parameter    :: min_bspline_order = dr_order + 2
  integer, parameter    :: max_bspline_order = 25

  if (bspl_order .lt. min_bspline_order) then
    write(error_msg,*) 'Spline order too small! Min = ', min_bspline_order
    call mol_mech_error
  else if (bspl_order .gt. max_bspline_order) then
    write(error_msg,*) 'bspline_order too big! Max = ', max_bspline_order
    call mol_mech_error
  end if

  if( numtasks .gt. 1) then
     amba_fft_alloc_size = max(xy_slab_dbl_cnt * max_xy_slab_cnt, &
                               zx_slab_dbl_cnt * max_zx_slab_cnt)
  else
     amba_fft_alloc_size = 2 * fft_x_dim * fft_y_dim * fft_z_dim
     max_recip_imgs = natom
     recip_img_lo = 1
     recip_img_hi = natom
     recip_img_range_wraps = .false.
  endif

  do_amoeba_recip_flag = ibset(do_amoeba_recip_flag, valid_bit)

  return

end subroutine amoeba_recip_final_setup

!*******************************************************************************!
! Subroutine:  amoeba_recip_step_alloc_dat
!
! Description: <TBS>
!
!*******************************************************************************

subroutine amoeba_recip_step_alloc_dat(num_ints, num_reals)

  use mdin_ewald_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed
  
  allocate(gbl_Qperm(amba_fft_alloc_size), &
           gbl_Q1(amba_fft_alloc_size), &
           gbl_Q2(amba_fft_alloc_size), &
           gbl_G_func(amba_fft_alloc_size), &
           gbl_theta1(0:dr_order, bspl_order, max_recip_imgs), &
           gbl_theta2(0:dr_order, bspl_order, max_recip_imgs), &
           gbl_theta3(0:dr_order, bspl_order, max_recip_imgs), &
           gbl_init_grid_ind(3, max_recip_imgs), &
           gbl_fractional_multipole(10, max_recip_imgs), &
           gbl_perm_F_field(20, max_recip_imgs), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(gbl_init_grid_ind)

  num_reals = num_reals + size(gbl_Qperm) + &
                          size(gbl_Q1) + &
                          size(gbl_Q2) + &
                          size(gbl_G_func) + &
                          size(gbl_theta1) + &
                          size(gbl_theta2) + &
                          size(gbl_theta3) + &
                          size(gbl_fractional_multipole) + &
                          size(gbl_perm_F_field)

    allocate(gbl_my_chgs(max_recip_imgs), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
    num_ints = num_ints + size(gbl_my_chgs)

  return

end subroutine amoeba_recip_step_alloc_dat

!*******************************************************************************!
! Subroutine:  amoeba_recip_step_dealloc_dat
!
! Description: <TBS>
!
!*******************************************************************************

subroutine amoeba_recip_step_dealloc_dat(num_ints, num_reals)

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

  num_ints = num_ints - size(gbl_init_grid_ind)

  num_reals = num_reals - size(gbl_Qperm) - &
                          size(gbl_Q1) - &
                          size(gbl_Q2) - &
                          size(gbl_G_func) - &
                          size(gbl_theta1) - &
                          size(gbl_theta2) - &
                          size(gbl_theta3) - &
                          size(gbl_fractional_multipole) - &
                          size(gbl_perm_F_field)

  deallocate(gbl_Qperm, &
             gbl_Q1, gbl_Q2, &
             gbl_G_func, &
             gbl_theta1, gbl_theta2, gbl_theta3, &
             gbl_init_grid_ind, &
             gbl_fractional_multipole, &
             gbl_perm_F_field)

    num_ints = num_ints - size(gbl_my_chgs)
    deallocate(gbl_my_chgs)

  return

end subroutine amoeba_recip_step_dealloc_dat

!*******************************************************************************!
! Subroutine:  am_recip_perm_field
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_perm_field(atm_cnt, crd, cart_dipole_field)

  use amoeba_multipoles_mod, only : global_multipole
  use amoeba_induced_mod, only : is_polarizable
  use amoeba_flags_mod
  use mdin_ewald_dat_mod
  use pme_recip_mod
  use pbc_mod
  use img_mod
  use parallel_dat_mod
  use timers_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: atm_cnt
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: cart_dipole_field(3, *)

! Local variables:

  double precision      :: FdipF(3 * max_recip_imgs)
  double precision      :: Qtmp(amba_fft_alloc_size)

  if (iand(do_amoeba_recip_flag, proceed_induce) .ne. proceed_induce) return

  ! Next get the bspline coeffs--later for MPI fill the select array
  ! for slabs i.e. those atoms landing in this proc's slab

  ! fill theta1-3; these are saved for use in induction scf & final energy, frc
  
  call am_recip_bspline_fill(atm_cnt, crd, gbl_img_atm_map, &
                             dr_order, bspl_order, gbl_init_grid_ind, &
                             gbl_theta1, gbl_theta2, gbl_theta3, &
                             gbl_my_chgs, my_chg_cnt)  ! gbl_my_chgs, my_chg_cnt are set here
                             
  call am_recip_global_to_fractional(atm_cnt, recip, nfft1, nfft2, nfft3,  &
                                     global_multipole, &
                                     gbl_fractional_multipole, &
                                     gbl_my_chgs, my_chg_cnt)

  call update_pme_time(bspline_timer)

  call am_recip_perm_fill_grid(atm_cnt, bspl_order, dr_order, &
                               gbl_fractional_multipole, &
                               gbl_theta1, gbl_theta2, gbl_theta3, &
                               gbl_init_grid_ind, gbl_Q1, &
                               gbl_my_chgs, my_chg_cnt)

  call update_pme_time(grid_charges_timer)

  call fft3drc_forward(gbl_Q1, Qtmp, fft_x_dim, fft_y_dim, fft_z_dim) !  FFT transform b-splined charge distribution  gbl_Q1

! Result saved in Qperm will be used for final virial, force calculation.

  call array_copy(gbl_Qperm, Qtmp, 2 * nfft3 * fft_x_dim * my_zx_slab_cnt)

  call update_pme_time(fft_timer)

  call am_recip_get_recip_gfunc(gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                                recip, uc_volume, ew_coeff,  &
                                nfft1, fft_x_dim, nfft2, nfft3, gbl_G_func)  ! find G-functions to multiply charge distribution FFT transform to find 
                                                                             ! fft transform of electric field   formula (2.46) of Darden_2004

  call am_recip_g_times_q(nfft1, fft_x_dim, nfft2, nfft3, gbl_G_func, Qtmp) ! multiply FFT transformed charge distribution by G -function to find
                                                                            ! FFT transform of electric field

  call update_pme_time(scalar_sum_timer)

  call fft3drc_back(Qtmp, gbl_Q1, fft_x_dim, fft_y_dim, fft_z_dim) ! back FFT transform to find electric field on the grid in normal space 

  call update_pme_time(fft_timer)

  call am_recip_get_perm_f_field(atm_cnt, is_polarizable, gbl_init_grid_ind, &
                                 dr_order, bspl_order, &
                                 gbl_theta1, gbl_theta2, gbl_theta3, & 
                                 gbl_Q1, gbl_perm_F_field, FdipF, &
                                 gbl_my_chgs, my_chg_cnt)  ! find electric field (monopole, dipole, quadrupole) on the atomic centers (fractional coordinates) 
                                 
  call am_recip_fdip_to_cdip_field(atm_cnt, is_polarizable, &
                                   nfft1, nfft2, nfft3, &
                                   recip, FdipF, cart_dipole_field, &
                                   gbl_my_chgs, my_chg_cnt)  ! find electric field on polarizable atomic centers ( cartesian coordinates ) 

  perm_field_done = .true.

  call update_pme_time(grad_sum_timer)

  return

end subroutine am_recip_perm_field

!*******************************************************************************!
! Subroutine:  am_recip_dipole_field
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_dipole_field(atm_cnt, crd, ind_dip1, ind_dip2, &
                                 dip_field1, dip_field2)

  use amoeba_flags_mod
  use amoeba_induced_mod, only : is_polarizable
  use mdin_ewald_dat_mod
  use pbc_mod
  use timers_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  double precision, intent(in)  :: crd(3, *)
  double precision, intent(in)  :: ind_dip1(3, *), ind_dip2(3, *)
  double precision, intent(out) :: dip_field1(3, *), dip_field2(3, *)

! Local variables:

  double precision      :: fdip1(3 * max_recip_imgs)
  double precision      :: fdip2(3 * max_recip_imgs)
  double precision      :: frac_field1(3 * max_recip_imgs)
  double precision      :: frac_field2(3 * max_recip_imgs)
  double precision      :: Qtmp(amba_fft_alloc_size)

  if (iand(do_amoeba_recip_flag, proceed_induce) .ne. proceed_induce) return

  call am_recip_dip_to_frac_dip(atm_cnt, is_polarizable, ind_dip1, ind_dip2, &
                                fdip1, fdip2, gbl_my_chgs, my_chg_cnt)

  call update_pme_time(bspline_timer)

  call am_recip_dipole_fill_grids(atm_cnt, bspl_order, dr_order, &
                                  is_polarizable, fdip1, fdip2, &
                                  gbl_theta1, gbl_theta2, gbl_theta3, &
                                  gbl_init_grid_ind, gbl_Q1, gbl_Q2, &
                                  gbl_my_chgs, my_chg_cnt)

  call update_pme_time(grid_charges_timer)

  call fft3drc_forward(gbl_Q1, Qtmp, fft_x_dim, fft_y_dim, fft_z_dim)

  call update_pme_time(fft_timer)

  call am_recip_g_times_q(nfft1, fft_x_dim, nfft2, nfft3, gbl_G_func, Qtmp)

  call update_pme_time(scalar_sum_timer)

  call fft3drc_back(Qtmp, gbl_Q1, fft_x_dim, fft_y_dim, fft_z_dim)

  call fft3drc_forward(gbl_Q2, Qtmp, fft_x_dim, fft_y_dim, fft_z_dim)

  call update_pme_time(fft_timer)

  call am_recip_g_times_q(nfft1, fft_x_dim, nfft2, nfft3, gbl_G_func, Qtmp)

  call update_pme_time(scalar_sum_timer)

  call fft3drc_back(Qtmp, gbl_Q2, fft_x_dim, fft_y_dim, fft_z_dim)

  call update_pme_time(fft_timer)

  call am_recip_get_2dip_f_fields(atm_cnt, is_polarizable, gbl_init_grid_ind, &
                                  dr_order, bspl_order, &
                                  gbl_theta1, gbl_theta2, gbl_theta3, &
                                  gbl_Q1, gbl_Q2, frac_field1, frac_field2, &
                                  gbl_my_chgs, my_chg_cnt)

  call am_recip_fdip_to_cdip_field(atm_cnt, is_polarizable, &
                                   nfft1, nfft2, nfft3, &
                                   recip, frac_field1, dip_field1, &
                                   gbl_my_chgs, my_chg_cnt)
  
  call am_recip_fdip_to_cdip_field(atm_cnt, is_polarizable, &
                                   nfft1, nfft2, nfft3, &
                                   recip, frac_field2, dip_field2, &
                                   gbl_my_chgs, my_chg_cnt)

  call update_pme_time(grad_sum_timer)

  return

end subroutine am_recip_dipole_field

!*******************************************************************************!
! Subroutine:  am_recip_ene_frc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_ene_frc(atm_cnt, crd, ind_dip1, ind_dip2, ene_perm, &
                            ene_ind, img_frc, virial, netfrcs)

  use amoeba_multipoles_mod, only : torque_field, global_multipole
  use amoeba_induced_mod, only : is_polarizable
  use amoeba_multipoles_mod, only : coulomb_const_kcal_per_mole
  use amoeba_flags_mod
  use img_mod
  use mdin_ewald_dat_mod
  use pbc_mod
  use timers_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: atm_cnt
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in)          :: ind_dip1(3, *)
  double precision, intent(in)          :: ind_dip2(3, *)
  double precision, intent(out)         :: ene_perm
  double precision, intent(out)         :: ene_ind
  double precision, intent(in out)      :: img_frc(3, *)
  double precision, intent(in out)      :: virial(3, 3)
  double precision, intent(out)         :: netfrcs(3)

! Local variables:

  integer                       :: alloc_failed
  double precision, allocatable :: cdf(:)
  double precision              :: fdip1(3 * max_recip_imgs)
  double precision              :: fdip2(3 * max_recip_imgs)
  double precision              :: frac_field1(20 * max_recip_imgs)
  double precision              :: frac_field2(20 * max_recip_imgs)
  double precision              :: Qtmp1(amba_fft_alloc_size)
  double precision              :: Qtmp2(amba_fft_alloc_size)

  ene_perm = 0.d0
  ene_ind = 0.d0

  if (iand(do_amoeba_recip_flag, proceed_postinduce) .ne. proceed_postinduce) &
    return

  ! This occurs if am_recip_perm_field was not called since last call to here,
  ! i.e. no induced dipoles we need permanent field for ene_perm.

  if (.not. perm_field_done) then
    allocate(cdf(3 * atm_cnt), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
    cdf(:) = 0.d0
    call am_recip_perm_field(atm_cnt, crd, cdf)

    deallocate(cdf)  ! do not use electric field at polariazable atom positions
  end if

  perm_field_done = .false.  ! reset for next force call
  
!  write(iunit_debug,*) " am_recip_ene_frc() pt 1  my_chg_cnt = ",my_chg_cnt

  call am_recip_dip_to_frac_dip(atm_cnt, is_polarizable, ind_dip1, ind_dip2, &
                                fdip1, fdip2, gbl_my_chgs, my_chg_cnt)
 
  call update_pme_time(bspline_timer)

  call am_recip_dipole_fill_grids(atm_cnt, bspl_order, dr_order, &
                                  is_polarizable, fdip1, fdip2,  &
                                  gbl_theta1, gbl_theta2, gbl_theta3, &
                                  gbl_init_grid_ind, gbl_Q1, gbl_Q2, &
                                  gbl_my_chgs, my_chg_cnt)                       ! gbl_Q1, gbl_Q2 - ind_dip1, ind_dip2  distributed on the grid
                                  
  call update_pme_time(grid_charges_timer)

  call fft3drc_forward(gbl_Q1, Qtmp1, fft_x_dim, fft_y_dim, fft_z_dim)  ! Qtmp1 - Fourie transform of ind_dip1 
  call fft3drc_forward(gbl_Q2, Qtmp2, fft_x_dim, fft_y_dim, fft_z_dim)  ! Qtmp1 - Fourie transform of ind_dip2

  call update_pme_time(fft_timer)

  ! recall gbl_Qperm was saved during perm field calc

  call am_recip_scalar_sum(recip, ew_coeff, nfft1, fft_x_dim, nfft2, nfft3, &
                           gbl_G_func, gbl_Qperm, Qtmp1, Qtmp2, virial)  ! find electric field of polarized charges (Qtmp1, Qtmp2) in k-space  
                                                                         ! k-space charge destribution of permanent charges gbl_Qperm is used for virial calculations

  call update_pme_time(scalar_sum_timer)

  call fft3drc_back(Qtmp1, gbl_Q1, fft_x_dim, fft_y_dim, fft_z_dim)
  call fft3drc_back(Qtmp2, gbl_Q2, fft_x_dim, fft_y_dim, fft_z_dim)  ! transform electric field of polarized dipoles into real space ( reciprocal space grid )
 
  call update_pme_time(fft_timer)

  ! get the field due to the two sets of dipoles
  ! (up to grad of quadrupole terms i.e. octupolar order)

  call am_recip_get_ind_f_fields(atm_cnt, gbl_init_grid_ind, dr_order, &
                                 bspl_order, &
                                 gbl_theta1, gbl_theta2, gbl_theta3,  &
                                 gbl_Q1, gbl_Q2, frac_field1, frac_field2, &    ! transform electric field from polarized dipoles on the grid ( gbl_Q1, gbl_Q2 ) 
                                 gbl_my_chgs, my_chg_cnt)                       ! to multipole fields at charge centers in fractional coordinates ( frac_field1, frac_field2 ) 

  call am_recip_perm_ene_grad(atm_cnt, nfft1, nfft2, nfft3, recip, &
                              gbl_perm_F_field, gbl_fractional_multipole, & 
                              ene_perm, img_frc, gbl_my_chgs, my_chg_cnt)
                              
  call am_recip_ind_ene_grad(atm_cnt, nfft1, nfft2, nfft3, &
                             is_polarizable, recip, &
                             gbl_perm_F_field, gbl_fractional_multipole, &
                             frac_field1, frac_field2, &                      ! frac_field1, frac_field2 - multipole fields from induced charges at charge centers in fractional coordinates ( frac_field1, frac_field2 ) 
                             fdip1, fdip2, ene_ind, img_frc, &   ! fdip1, fdip2 - induced dipoles in fractional coordinates
                             gbl_my_chgs, my_chg_cnt)

  call am_recip_field_torque_virial(atm_cnt, is_polarizable, &
                                    nfft1, nfft2, nfft3, recip, &
                                    gbl_perm_F_field, &
                                    frac_field1, frac_field2, &
                                    torque_field, virial, global_multipole, &
                                    ind_dip1, ind_dip2, gbl_my_chgs, my_chg_cnt)

  call update_pme_time(grad_sum_timer)

  virial(:, :) = coulomb_const_kcal_per_mole * virial(:, :)

  call am_recip_calc_netfrcs(atm_cnt, img_frc, netfrcs, gbl_my_chgs, my_chg_cnt)

  return

end subroutine am_recip_ene_frc

!*******************************************************************************!
! Subroutine:  am_recip_perm_fill_grid
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_perm_fill_grid(img_cnt, bspline_order, dr_order, &
                                   Fmpole, theta1, theta2, theta3, &
                                   init_grid_ind, Q, &
                                   my_chgs, my_chg_cnt)

  use mdin_ewald_dat_mod
  use pme_fft_mod
  use img_mod

  implicit none

  ! Note hardwired for quadrupoles--

  ! Formal arguments:

  integer, intent(in)           :: img_cnt
  integer, intent(in)           :: bspline_order
  integer, intent(in)           :: dr_order
  double precision, intent(in)  :: Fmpole(10, *)
  double precision, intent(in)  :: theta1(0:dr_order, bspline_order, *)
  double precision, intent(in)  :: theta2(0:dr_order, bspline_order, *)
  double precision, intent(in)  :: theta3(0:dr_order, bspline_order, *)
  integer, intent(in)           :: init_grid_ind(3, *)
  double precision, intent(out) :: Q(2 * fft_x_dim, fft_y_dim, my_xy_slab_cnt)
  integer, intent(in)           :: my_chgs(*), my_chg_cnt
                       
  integer                       :: img_id, my_imgs_idx
  integer                       :: igrd0, jgrd0, kgrd0
  integer                       :: ith1, ith2, ith3
  integer                       :: i, i0, j, j0, k, k0, kq
  integer                       :: kbot0, kbot, ktop
  double precision              :: t0, t1, t2
  double precision              :: u0, u1, u2
  double precision              :: v0, v1, v2
  double precision              :: term0, term1, term2

  Q(:,:,:) = 0.d0

  kbot0 = my_xy_slab_start
  kbot = kbot0 + 1
  ktop = kbot0 + my_xy_slab_cnt
  
  do my_imgs_idx = 1, my_chg_cnt
    img_id = my_chgs(my_imgs_idx)

    igrd0 = init_grid_ind(1, my_imgs_idx) ! begin index in 1st direction
    jgrd0 = init_grid_ind(2, my_imgs_idx) ! begin index in 2nd direction
    kgrd0 = init_grid_ind(3, my_imgs_idx) ! begin index in 3rd direction
    k0 = kgrd0
    do ith3 = 1, bspline_order
      k0 = k0 + 1
      k = k0 + 1 + (nfft3 - isign(nfft3, k0)) / 2

      if ( k .ge. kbot .and. k .le. ktop ) then
        kq = k - kbot0
        
        j0 = jgrd0
        v0 = theta3(0, ith3, my_imgs_idx) ! theta3
        v1 = theta3(1, ith3, my_imgs_idx) ! 1st deriv of theta3
        v2 = theta3(2, ith3, my_imgs_idx) ! 2nd deriv of theta3
        do ith2 = 1, bspline_order
          j0 = j0 + 1
          j = j0 + 1 + (nfft2 - isign(nfft2, j0)) / 2
          u0 = theta2(0, ith2, my_imgs_idx) ! theta2
          u1 = theta2(1, ith2, my_imgs_idx) ! 1st deriv of theta2
          u2 = theta2(2, ith2, my_imgs_idx) ! 2nd deriv of theta2
          i0 = igrd0

! hardwire our knowledge of layout of theta1, 2, 3 to pre-assemble factors

          term0 = Fmpole(Ind_000, my_imgs_idx) * u0 * v0 + &
                  Fmpole(Ind_010, my_imgs_idx) * u1 * v0 + &
                  Fmpole(Ind_001, my_imgs_idx) * u0 * v1 + &
                  Fmpole(Ind_020, my_imgs_idx) * u2 * v0 + &
                  Fmpole(Ind_002, my_imgs_idx) * u0 * v2 + &
                  Fmpole(Ind_011, my_imgs_idx) * u1 * v1

          term1 = Fmpole(Ind_100, my_imgs_idx) * u0 * v0 + &
                  Fmpole(Ind_110, my_imgs_idx) * u1 * v0 + &
                  Fmpole(Ind_101, my_imgs_idx) * u0 * v1

          term2 = Fmpole(Ind_200, my_imgs_idx) * u0 * v0

          do ith1 = 1, bspline_order
            i0 = i0 + 1
            i = i0 + 1 + (nfft1 - isign(nfft1, i0)) / 2
            t0 = theta1(0, ith1, my_imgs_idx) ! theta1
            t1 = theta1(1, ith1, my_imgs_idx) ! 1st deriv of theta1
            t2 = theta1(2, ith1, my_imgs_idx) ! 2nd deriv of theta1
            Q(i, j, kq) = Q(i, j, kq) + term0 * t0 + term1 * t1 + term2 * t2
          end do
        end do ! ith2 = 1, bspline_order
      end if
    end do ! ith3 = 1, bspline_order
  end do ! n = 1, img_cnt

  return

end subroutine am_recip_perm_fill_grid

!*******************************************************************************!
! Subroutine:  am_recip_dipole_fill_grids
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_dipole_fill_grids(img_cnt, bspline_order, dr_order, &
                                      is_polarizable, fdip1, fdip2,  &
                                      theta1, theta2, theta3, &
                                      init_grid_ind, Q1, Q2, &
                                      my_chgs, my_chg_cnt)

  use mdin_ewald_dat_mod
  use pme_fft_mod
  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: img_cnt
  integer, intent(in)           :: bspline_order
  integer, intent(in)           :: dr_order
  logical, intent(in)           :: is_polarizable(*)
  double precision, intent(in)  :: fdip1(3, *) ! induced dipoles  (ind_dip_p ) in reciprocal space  
  double precision, intent(in)  :: fdip2(3, *) ! induced dipoles  (ind_dip_d ) in reciprocal space  
  double precision, intent(in)  :: theta1(0:dr_order, bspline_order, *)  ! coefficients to convert to convert dipoles to grid charges??? using splines  
  double precision, intent(in)  :: theta2(0:dr_order, bspline_order, *)  ! for each charge in 3 directions
  double precision, intent(in)  :: theta3(0:dr_order, bspline_order, *)  !
  integer, intent(in)           :: init_grid_ind(3, *)  !  begin indexes in 3 directions for specific charges ?? 
  double precision, intent(out) :: Q1(2 * fft_x_dim, fft_y_dim, my_xy_slab_cnt) ! b-spline charge expansions on the grid corresponding to induced dipoles ind_dip_p (ifip1)
  double precision, intent(out) :: Q2(2 * fft_x_dim, fft_y_dim, my_xy_slab_cnt) ! b-spalne charge exapnsions on the grid corresponding to induced dipoles ind_dip_d (ifip2)
  integer, intent(in)           :: my_chgs(*), my_chg_cnt  ! charges that belong to a given MPI processor ??

! Local variables:

  integer                       :: img_id, atm_id, my_imgs_idx
  integer                       :: igrd0, jgrd0, kgrd0
  integer                       :: ith1, ith2, ith3
  integer                       :: i, i0, j, j0, k, k0, kq
  integer                       :: kbot0, kbot, ktop
  double precision              :: t0, t1
  double precision              :: u0, u1
  double precision              :: v0, v1
  double precision              :: term1_0, term1_1
  double precision              :: term2_0, term2_1

  Q1(:,:,:) = 0.d0
  Q2(:,:,:) = 0.d0

  kbot0 = my_xy_slab_start
  kbot = kbot0 + 1
  ktop = kbot0 + my_xy_slab_cnt

  do my_imgs_idx = 1, my_chg_cnt
    img_id = my_chgs(my_imgs_idx)

    atm_id = gbl_img_atm_map(img_id)

    if (is_polarizable(atm_id)) then
      igrd0 = init_grid_ind(1, my_imgs_idx) ! begin index in 1st direction
      jgrd0 = init_grid_ind(2, my_imgs_idx) ! begin index in 2nd direction
      kgrd0 = init_grid_ind(3, my_imgs_idx) ! begin index in 3rd direction
      k0 = kgrd0
      do ith3 = 1, bspline_order
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3, k0)) / 2

        if ( k .ge. kbot .and. k .le. ktop ) then
          kq = k - kbot0
          
          j0 = jgrd0
          v0 = theta3(0, ith3, my_imgs_idx) ! theta3
          v1 = theta3(1, ith3, my_imgs_idx) ! 1st deriv of theta3
          do ith2 = 1, bspline_order
            j0 = j0 + 1
            j = j0 + 1 + (nfft2 - isign(nfft2, j0)) / 2
            u0 = theta2(0, ith2, my_imgs_idx) ! theta2
            u1 = theta2(1, ith2, my_imgs_idx) ! 1st deriv of theta2
            i0 = igrd0
            term1_0 = fdip1(2, my_imgs_idx) * u1 * v0 + &
                      fdip1(3, my_imgs_idx) * u0 * v1
            term2_0 = fdip2(2, my_imgs_idx) * u1 * v0 + &
                      fdip2(3, my_imgs_idx) * u0 * v1
            term1_1 = fdip1(1, my_imgs_idx) * u0 * v0
            term2_1 = fdip2(1, my_imgs_idx) * u0 * v0
            do ith1 = 1, bspline_order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1, i0)) / 2  ! if i0 < 0  then i = i0 + 1 + nfft1 
              t0 = theta1(0, ith1, my_imgs_idx) ! theta1
              t1 = theta1(1, ith1, my_imgs_idx) ! 1st deriv of theta1
              Q1(i, j, kq) = Q1(i, j, kq) + term1_0 * t0 + term1_1 * t1
              Q2(i, j, kq) = Q2(i, j, kq) + term2_0 * t0 + term2_1 * t1
            end do
          end do ! ith2 = 1, bspline_order
        end if
      end do ! ith3 = 1, bspline_order
    end if ! is_polarizable(atm_id)
  end do ! n = 1, img_cnt

  return

end subroutine am_recip_dipole_fill_grids

!*******************************************************************************
!
! Subroutine:  am_recip_bspline_fill
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine am_recip_bspline_fill(img_cnt, crd, img_atm_map, &
                                 dr_order, bspline_order, init_grid_ind, &
                                 theta1, theta2, theta3, &
                                 my_chgs, my_chg_cnt)

  use pme_fft_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: img_cnt
  double precision      :: crd(3, *)
  integer               :: img_atm_map(*)
  integer               :: dr_order
  integer               :: bspline_order
  integer               :: init_grid_ind(3, *)
  double precision      :: theta1(0:dr_order, bspline_order, *)
  double precision      :: theta2(0:dr_order, bspline_order, *)
  double precision      :: theta3(0:dr_order, bspline_order, *)
  integer               :: my_chgs(*), my_chg_cnt

! Local variables:

  double precision      :: crd_x, crd_y, crd_z
  double precision      :: fract(3)
  double precision      :: weight(3)
  double precision      :: recip_stk(3, 3)
  integer               :: atm_id
  integer               :: i, k

  integer               :: img_lo, img_hi, img_lo2, img_hi2
  integer               :: range_cnt, range_ctr
  integer               :: kbot0, ktop1
  logical               :: kwraps

  my_chg_cnt = 0

  if( numtasks .gt. 1) then  
    if (recip_img_range_wraps) then
        range_cnt = 2
        img_lo = recip_img_lo
        img_hi = img_cnt
        img_lo2 = 1
        img_hi2 = recip_img_hi
    else
        range_cnt = 1
        img_lo = recip_img_lo
        img_hi = recip_img_hi
    end if

    kbot0 = xy_slab_start(mytaskid) 
    ktop1 = kbot0 + my_xy_slab_cnt + bspline_order - 2
    kwraps = (ktop1 .ge. nfft3)
  else
    range_cnt = 1
    img_lo = 1
    img_hi = img_cnt
  endif
 
  recip_stk(:,:) = recip(:,:)
  
!  write(iunit_debug,*)" am_recip_bspline_fill() pt 1  recip_img_range_wraps = ", &
!                        recip_img_range_wraps
  
!  write(iunit_debug,*)" am_recip_bspline_fill() pt 1  kwraps, kbot0,ktop1,nfft3 = ", &
!                            kwraps, kbot0,ktop1,nfft3
                            
!  write(iunit_debug,*)" am_recip_bspline_fill() pt 1  range_ctr, range_cnt,img_lo,img_hi = ", &
!                            range_ctr, range_cnt,img_lo,img_hi                            

  do range_ctr = 1, range_cnt
  do i = img_lo, img_hi
  
  if (img_atm_map(i) .lt. 0) cycle    ! Skip unused images.

! Unfortunately we need fractional coords that are not available without going
! back to atm_crd().  This is the case because the coordinates in img_crd
! were based on fractional coords at list build time, but as movement occurs
! you can't reconstruct the fractional.

    atm_id = img_atm_map(i)
    crd_x = crd(1, atm_id)
    crd_y = crd(2, atm_id)
    crd_z = crd(3, atm_id)

    fract(1) = crd_x * recip_stk(1,1) + &
               crd_y * recip_stk(2,1) +  &
               crd_z * recip_stk(3,1)

    fract(1) = dble(nfft1) * (fract(1) - dnint(fract(1)) + 0.5d0)

    fract(2) = crd_x * recip_stk(1,2) + &
               crd_y * recip_stk(2,2) + &
               crd_z * recip_stk(3,2)

    fract(2) = dble(nfft2) * (fract(2) - dnint(fract(2)) + 0.5d0)

    fract(3) = crd_x * recip_stk(1,3) + &
               crd_y * recip_stk(2,3) + &
               crd_z * recip_stk(3,3)

    fract(3) = dble(nfft3) * (fract(3) - dnint(fract(3)) + 0.5d0)

    if( numtasks .gt. 1) then
        k = int(fract(3))
        
!        write(iunit_debug,*)" am_recip_bspline_fill() pt 3  k  = ", k

        if (kwraps) then
!            write(iunit_debug,*)" am_recip_bspline_fill() pt 4 "
            if (k .lt. kbot0 .and. k .gt. (ktop1 - nfft3) ) cycle
!            write(iunit_debug,*)" am_recip_bspline_fill() pt 5 "
        else
!            write(iunit_debug,*)" am_recip_bspline_fill() pt 6 "
            if (k .lt. kbot0 .or. k .gt. ktop1) cycle
!            write(iunit_debug,*)" am_recip_bspline_fill() pt 7 "
        end if

        my_chg_cnt = my_chg_cnt + 1
        my_chgs(my_chg_cnt) = i
    else
        my_chg_cnt = my_chg_cnt + 1
        my_chgs(my_chg_cnt) = i
    endif

    weight(:) = fract(:) - int(fract(:))

    if( numtasks .gt. 1) then
        init_grid_ind(:, my_chg_cnt) = int(fract(:)) - bspline_order

        call am_recip_bspline_fill_gen(weight(1), bspline_order, dr_order, &
                                       theta1(0, 1, my_chg_cnt))
        call am_recip_bspline_fill_gen(weight(2), bspline_order, dr_order, &
                                       theta2(0, 1, my_chg_cnt))
        call am_recip_bspline_fill_gen(weight(3), bspline_order, dr_order, &
                                       theta3(0, 1, my_chg_cnt))

    else
        init_grid_ind(:, i) = int(fract(:)) - bspline_order

        call am_recip_bspline_fill_gen(weight(1), bspline_order, dr_order, &
                                       theta1(0, 1, i))
        call am_recip_bspline_fill_gen(weight(2), bspline_order, dr_order, &
                                       theta2(0, 1, i))
        call am_recip_bspline_fill_gen(weight(3), bspline_order, dr_order, &
                                       theta3(0, 1, i))
    endif

  end do

  if (numtasks .gt. 1 .and. recip_img_range_wraps) then
    img_lo = img_lo2
    img_hi = img_hi2
  end if

  end do

  return

end subroutine am_recip_bspline_fill

!*******************************************************************************!
! Subroutine:  am_recip_bspline_fill_gen
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_bspline_fill_gen(w, spline_order, dr_order, new_array)

  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision      :: w
  integer               :: spline_order
  integer               :: dr_order
  double precision      :: new_array(0:dr_order, spline_order)

! Local variables:

  double precision      :: array(spline_order, spline_order)
  integer               :: k, j

! init spline_order 2

  array(2, 2) = w      
  array(1, 2) = 1.d0 - w      

! one pass to spline_order 3

  array(3, 3) = 0.5d0 * w * array(2, 2)
  array(2, 3) = 0.5d0 * ((w + 1.d0) * array(1, 2) + (2.d0-w) * array(2, 2))
  array(1, 3) = 0.5d0 * (1.d0-w) * array(1, 2)

! compute standard b-spline recursion

  do k = 4, spline_order
    call am_recip_bspline_one_pass_recur(w, k, array(1, k-1), array(1, k)) 
  end do

! do derivatives

  if (dr_order .gt. 0) then
    call am_recip_bspline_diff(array(1, spline_order-1), spline_order)
    if (dr_order .gt. 1) then
      call am_recip_bspline_diff(array(1, spline_order-2), spline_order-1)
      call am_recip_bspline_diff(array(1, spline_order-2), spline_order)
      if (dr_order .gt. 2) then
        call am_recip_bspline_diff(array(1, spline_order-3), spline_order-2)
        call am_recip_bspline_diff(array(1, spline_order-3), spline_order-1)
        call am_recip_bspline_diff(array(1, spline_order-3), spline_order)
        if (dr_order .gt. 3) then
          call am_recip_bspline_diff(array(1, spline_order-4), spline_order-3)
          call am_recip_bspline_diff(array(1, spline_order-4), spline_order-2)
          call am_recip_bspline_diff(array(1, spline_order-4), spline_order-1)
          call am_recip_bspline_diff(array(1, spline_order-4), spline_order)
          if (dr_order .gt. 4) then
            call am_recip_bspline_diff(array(1, spline_order-5), spline_order-4)
            call am_recip_bspline_diff(array(1, spline_order-5), spline_order-3)
            call am_recip_bspline_diff(array(1, spline_order-5), spline_order-2)
            call am_recip_bspline_diff(array(1, spline_order-5), spline_order-1)
            call am_recip_bspline_diff(array(1, spline_order-5), spline_order)
            if (dr_order .gt. 5) then
               error_msg =  'spline derivs of order > 5 not implemented!'
               call mol_mech_error
            end if !(dr_order .gt. 5) then
          end if !(dr_order .gt. 4) then
        end if !(dr_order .gt. 3) then
      end if !(dr_order .gt. 2) then
    end if !(dr_order .gt. 1) then
  end if !(dr_order .gt. 0) then

! re-arrange array

  do k = 1, spline_order
    do j = 0, dr_order
      new_array(j, k) = array(k, spline_order-j)
    end do
  end do

  return

end subroutine am_recip_bspline_fill_gen

!*******************************************************************************!
! Subroutine:  am_recip_bspline_one_pass_recur
!
! Description:
!
! Using notation from Essmann et al; w = u-[u] and
! array(j) = M_n(w + order - j)  where n is order
! RECURSION:  M_n(w) = (w/(n-1)) * M_n-1(w)+((n-w)/(n-1)) * M_n-1(w-1)
! so M_n(w+n-j) = ((w+n-j)/(n-1)) * M_n-1(w+n-j)+((j-w)/(n-1)) * M_n-1(w+n-j-1)
! and new(j) = ((w+n-j)/(n-1)) * old(j-1) + ((j-w)/(n-1)) * old(j)
! where old is array before one_pass (thus n->n-1) and new is array afterwards
!
!*******************************************************************************

subroutine am_recip_bspline_one_pass_recur(w, n, old, new)

  implicit none

! Formal arguments:

  double precision      :: w
  integer               :: n
  double precision      :: old(*)
  double precision      :: new(*)

! Local variables:

  double precision      :: div
  integer               :: j

  div = 1.d0 / (n - 1)
  new(n) = div * w * old(n - 1)
  do j = 1, n - 2
    new(n - j) = div * ((w+j) * old(n - j - 1) + (n - j - w) * old(n - j))
  end do
  new(1) = div * (1 - w) * old(1)

  return

end subroutine am_recip_bspline_one_pass_recur

!*******************************************************************************!
! Subroutine:  am_recip_bspline_diff
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_bspline_diff(c, n)

  implicit none

! Formal arguments:

  double precision      :: c(*)
  integer               :: n

! Using notation from Essmann et al; w = u-[u] and
! array(j) = M_n(w + order - j)  where n is order
! DERIVATIVE:    d/dw M_n(w) = M_n-1(w) - M_n-1(w-1)
! i.e.   d/dw M_n(w+n-j) = M_n-1(w+n-j) - M_n-1(w+n-j-1)
! i.e.   new(j) = old(j-1) - old(j)
! where old is array before one_pass (thus n->n-1) and new is array afterwards
! do backwards to do in place

! Local variables:

  integer               :: j

  c(n) = c(n - 1)
  do j = n - 1, 2, -1
    c(j) = c(j - 1) - c(j)
  end do
  c(1) = -c(1)

  return

end subroutine am_recip_bspline_diff

!*******************************************************************************!
! Subroutine:  am_recip_get_recip_gfunc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_get_recip_gfunc(prefac1, prefac2, prefac3,  &
                                    recip, volume, ewald_coeff,  &
                                    nfft1, fft_x_dim, nfft2, nfft3, G)

  implicit none

! Formal arguments:

  integer, intent(in)           :: nfft1, fft_x_dim, nfft2, nfft3
  double precision, intent(in)  :: prefac1(nfft1)
  double precision, intent(in)  :: prefac2(nfft2)
  double precision, intent(in)  :: prefac3(nfft3)
  double precision, intent(in)  :: recip(3, 3)
  double precision, intent(in)  :: volume
  double precision, intent(in)  :: ewald_coeff
  double precision, intent(out) :: G(nfft3, fft_x_dim, my_zx_slab_cnt)

! Local variables:

  integer                       :: k1, k2, k2q, k3, k10
  integer                       :: m1, m2, m3, nf1, nf2, nf3
  double precision              :: pi, fac, mhat1, mhat2, mhat3, msq, denom

  pi = 3.14159265358979323846d0
  fac = pi**2 / ewald_coeff**2
  nf1 = nfft1 / 2
  if (2 * nf1 .lt. nfft1) nf1 = nf1 + 1
  nf2 = nfft2 / 2
  if (2 * nf2 .lt. nfft2) nf2 = nf2 + 1
  nf3 = nfft3 / 2
  if (2 * nf3 .lt. nfft3) nf3 = nf3 + 1

  do k2q = 1, my_zx_slab_cnt
    k2 = k2q + my_zx_slab_start

    m2 = k2 - 1
    if (k2 .gt. nf2) m2 = k2 - 1 - nfft2
    do k3 = 1, nfft3
      m3 = k3 - 1
      if (k3 .gt. nf3) m3 = k3 - 1 - nfft3
      k10 = 1
      if (my_zx_slab_start .eq. 0) then
        if (k3 + k2 .eq. 2) k10 = 2
      end if
      do k1 = k10, nf1 + 1
        m1 = k1 - 1
        if (k1 .gt. nf1) m1 = k1 - 1 - nfft1
        mhat1 = recip(1, 1) * m1 + recip(1, 2) * m2 + recip(1, 3) * m3
        mhat2 = recip(2, 1) * m1 + recip(2, 2) * m2 + recip(2, 3) * m3
        mhat3 = recip(3, 1) * m1 + recip(3, 2) * m2 + recip(3, 3) * m3
        msq = mhat1 * mhat1 + mhat2 * mhat2 + mhat3 * mhat3              ! msq = m^2  - 
        denom = pi * volume * msq                                        ! formula (2.45) of JCP_04_120_00073 ( Darden_2004 )
        G(k3, k1, k2q) = exp(-fac * msq)* prefac1(k1) * &                ! G(k3,k1,k2) =  1/(pi * V * m^2) ) * exp( -(pi^2 * m^2 / beta^2) ) * prefac1(k1) * prefac1(k2) * prefac1(k3)
                        prefac2(k2) * prefac3(k3) / denom                ! prefac1(k1) = (b1(m1/K1))^2  prefac1(k2) = (b2(m2/K1))^2  prefac1(k3) = (b3(m3/K1))^2
      end do ! k1 = k10, nf1 + 1
    end do ! k3 = 1, nfft3
  end do ! k2 = 1, nfft2

  return

end subroutine am_recip_get_recip_gfunc

!*******************************************************************************!
! Subroutine:  am_recip_g_times_q
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_g_times_q(nfft1, fft_x_dim, nfft2, nfft3, G, Q)

  implicit none

! Formal arguments:

  integer, intent(in)           :: nfft1, fft_x_dim, nfft2, nfft3
  double precision, intent(in)  :: G(nfft3, fft_x_dim, my_zx_slab_cnt)
  double precision, intent(out) :: Q(2, nfft3, fft_x_dim, my_zx_slab_cnt)

! Local variables:

  integer                       :: k1, k2, k2q, k3, k10, nf1

  nf1 = nfft1 / 2
  if (2 * nf1 .lt. nfft1) nf1 = nf1 + 1

!........Insist that Q(1, 1, 1, 1) is set to 0 (true already for neutral)

  if (my_zx_slab_start .eq. 0) then
    Q(1, 1, 1, 1) = 0.d0
    Q(2, 1, 1, 1) = 0.d0
  end if

  do k2q = 1, my_zx_slab_cnt
    k2 = k2q + my_zx_slab_start

    do k3 = 1, nfft3
      k10 = 1
      if (my_zx_slab_start .eq. 0) then
        if (k3 + k2 .eq. 2) k10 = 2
      end if
      do k1 = k10, nf1 + 1
        Q(1, k3, k1, k2q) = G(k3, k1, k2q) * Q(1, k3, k1, k2q)
        Q(2, k3, k1, k2q) = G(k3, k1, k2q) * Q(2, k3, k1, k2q)
      end do ! k1 = k10, nf1 + 1
    end do ! k3 = 1, nfft3
  end do ! k2 = 1, nfft2

  return

end subroutine am_recip_g_times_q

!*******************************************************************************!
! Subroutine:  am_recip_scalar_sum
!
! Description: <TBS>
!
!*******************************************************************************

subroutine  am_recip_scalar_sum(recip, ewald_coeff, &
                               nfft1, fft_x_dim, nfft2, nfft3, &
                               G, Q, Q1, Q2, virial)                       

  implicit none

! Formal arguments:

  double precision, intent(in)          :: recip(3, 3)
  double precision, intent(in)          :: ewald_coeff
  integer, intent(in)                   :: nfft1, fft_x_dim, nfft2, nfft3
  double precision, intent(in)          :: G(nfft3,fft_x_dim,my_zx_slab_cnt)
  double precision, intent(in)          :: Q(2,nfft3,fft_x_dim,my_zx_slab_cnt)      !  Fourie transform of permanent charge distribution 
  double precision, intent(in out)      :: Q1(2,nfft3,fft_x_dim,my_zx_slab_cnt)     !  Fourie transform of ind_dip_d (induced charge distribution 1)
  double precision, intent(in out)      :: Q2(2,nfft3,fft_x_dim,my_zx_slab_cnt)     !  Fourie transform of ind_dip_p (induced charge distribution 1)
  double precision, intent(out)         :: virial(3, 3)

! Local variables:

  integer                               :: k1, k2, k2q, k3, k10
  integer                               :: m1, m2, m3
  integer                               :: nf1, nf2, nf3
  integer                               :: k1s, k2s, k3s
  integer                               :: m1s, m2s, m3s
  double precision                      :: pi
  double precision                      :: fac
  double precision                      :: mhat1, mhat2, mhat3
  double precision                      :: msq
  double precision                      :: struc2
  double precision                      :: eterm
  double precision                      :: vterm
  double precision                      :: q11, q12, q21, q22
  double precision                      :: tmp1, tmp2
  double precision                      :: mult
  double precision                      :: vxx, vxy, vxz, vyy, vyz, vzz

  pi = 3.14159265358979323846d0

  fac = pi**2 / ewald_coeff**2

  nf1 = nfft1 / 2
  if (2 * nf1 .lt. nfft1) nf1 = nf1 + 1
  nf2 = nfft2 / 2
  if (2 * nf2 .lt. nfft2) nf2 = nf2 + 1
  nf3 = nfft3 / 2
  if (2 * nf3 .lt. nfft3) nf3 = nf3 + 1

  vxx = 0.d0
  vxy = 0.d0
  vxz = 0.d0
  vyy = 0.d0
  vyz = 0.d0
  vzz = 0.d0

  do k2q = 1, my_zx_slab_cnt
    k2 = k2q + my_zx_slab_start

    m2 = k2 - 1
    if (k2 .gt. nf2) m2 = k2 - 1 - nfft2
    
    do k3 = 1, nfft3
      m3 = k3 - 1
      if (k3 .gt. nf3) m3 = k3 - 1 - nfft3
      k10 = 1
      if (my_zx_slab_start .eq. 0) then
        if (k3 + k2 .eq. 2) k10 = 2
      end if
      do k1 = k10, nf1 + 1
        if (k1 .gt. 1) then
          mult = 2.d0
        else
          mult = 1.d0
        end if
        m1 = k1 - 1
        if (k1 .gt. nf1) m1 = k1 - 1 - nfft1
        mhat1 = recip(1, 1) * m1 + recip(1, 2) * m2 + recip(1, 3) * m3
        mhat2 = recip(2, 1) * m1 + recip(2, 2) * m2 + recip(2, 3) * m3
        mhat3 = recip(3, 1) * m1 + recip(3, 2) * m2 + recip(3, 3) * m3
        msq = mhat1 * mhat1 + mhat2 * mhat2 + mhat3 * mhat3
        eterm = mult * G(k3, k1, k2q)
        vterm = 2.d0 * (fac * msq + 1.d0) / msq            ! 2*(1 + pi^2 * m^2/beta^2)/m^2 
        q11 = Q(1, k3, k1, k2q) + Q1(1, k3, k1, k2q)
        q12 = Q(2, k3, k1, k2q) + Q1(2, k3, k1, k2q)
        q21 = Q(1, k3, k1, k2q) + Q2(1, k3, k1, k2q)
        q22 = Q(2, k3, k1, k2q) + Q2(2, k3, k1, k2q)
        struc2 = q11 * q21 + q12 * q22                     ! S(m)* S(-m) 
        tmp1 = eterm * struc2                              ! G(k3, k1, k2q)*S(m)* S(-m) =  1/(2*pi*V) * exp( - pi^2 * m^2/beta^2 ) /m^2 
        tmp2 = tmp1 * vterm
        vxx = vxx + tmp2 * mhat1 * mhat1 - tmp1
        vxy = vxy + tmp2 * mhat1 * mhat2
        vxz = vxz + tmp2 * mhat1 * mhat3
        vyy = vyy + tmp2 * mhat2 * mhat2 - tmp1
        vyz = vyz + tmp2 * mhat2 * mhat3
        vzz = vzz + tmp2 * mhat3 * mhat3 - tmp1
        ! transform Q1, Q2. No need to do Q--done in permanent field routine
        Q1(1, k3, k1, k2q) = G(k3, k1, k2q) * Q1(1, k3, k1, k2q)
        Q1(2, k3, k1, k2q) = G(k3, k1, k2q) * Q1(2, k3, k1, k2q)
        Q2(1, k3, k1, k2q) = G(k3, k1, k2q) * Q2(1, k3, k1, k2q)
        Q2(2, k3, k1, k2q) = G(k3, k1, k2q) * Q2(2, k3, k1, k2q)
      end do ! k1 = k10, nf1 + 1
    end do ! k3 = 1, nfft3
  end do ! k2 = 1, my_zx_slab_cnt
  virial(1, 1) = virial(1, 1) + 0.5d0 * vxx
  virial(1, 2) = virial(1, 2) + 0.5d0 * vxy
  virial(2, 1) = virial(2, 1) + 0.5d0 * vxy
  virial(1, 3) = virial(1, 3) + 0.5d0 * vxz
  virial(3, 1) = virial(3, 1) + 0.5d0 * vxz
  virial(2, 2) = virial(2, 2) + 0.5d0 * vyy
  virial(2, 3) = virial(2, 3) + 0.5d0 * vyz
  virial(3, 2) = virial(3, 2) + 0.5d0 * vyz
  virial(3, 3) = virial(3, 3) + 0.5d0 * vzz

  return

end subroutine am_recip_scalar_sum

!*******************************************************************************!
! Subroutine:  am_recip_get_perm_f_field
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_get_perm_f_field(img_cnt, is_polarizable, &
                                     init_grid_ind, dr_order, bspline_order, &
                                     theta1, theta2, theta3,  &
                                     Q_p, Fperm_field, Fperm_dipfield, &
                                     my_chgs, my_chg_cnt)

  use img_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pme_fft_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: img_cnt
  logical, intent(in)           :: is_polarizable(*)
  integer, intent(in)           :: init_grid_ind(3, *)
  integer, intent(in)           :: dr_order
  integer, intent(in)           :: bspline_order
  double precision, intent(in)  :: theta1(0:dr_order, bspline_order, *)
  double precision, intent(in)  :: theta2(0:dr_order, bspline_order, *)
  double precision, intent(in)  :: theta3(0:dr_order, bspline_order, *)
  double precision, intent(in)  :: Q_p(2 * fft_x_dim, fft_y_dim, my_xy_slab_cnt)
  double precision, intent(out) :: Fperm_field(20, *)
  double precision, intent(out) :: Fperm_dipfield(3, *)
  integer, intent(in)           :: my_chgs(*), my_chg_cnt

! Note order of field elements
! 1 Field   corresponds to 000 in spline indices
! 2-4 -> F_x, F_y, F_z respectively or to 100, 010, 001 in indices
! 5-10 -> F_xx, F_yy, F_zz, F_xy, F_xz, F_yz resp. or to 
!         200, 020, 002, 110, 101, 011 in indices
! 11-20 -> F_xxx, F_yyy, F_zzz, F_xxy, F_xxz, F_xyy, F_yyz, F_xzz, F_yzz,
!          F_xyz or to 300, 030, 003, 210, 201, 120, 021, 102, 012, 111
!          in indices

! Local variables:

  double precision              :: tq_p, t_p(0:3)
  double precision              :: u(0:3), v(0:3)
  double precision              :: tu_p(10), tuv_p(20)
  integer                       :: img_id, atm_id, my_imgs_idx
  integer                       :: m, i, i0, j, j0, k, k0, kq
  integer                       :: kbot0, kbot, ktop
  integer                       :: igrd0, jgrd0, kgrd0
  integer                       :: ith1, ith2, ith3

  kbot0 = my_xy_slab_start
  kbot = kbot0 + 1
  ktop = kbot0 + my_xy_slab_cnt

  do my_imgs_idx = 1, my_chg_cnt
    img_id = my_chgs(my_imgs_idx)

    igrd0 = init_grid_ind(1, my_imgs_idx) ! begin index in 1st direction
    jgrd0 = init_grid_ind(2, my_imgs_idx) ! begin index in 2nd direction
    kgrd0 = init_grid_ind(3, my_imgs_idx) ! begin index in 3rd direction

    do m = 1, 20
      tuv_p(m) = 0.d0
    end do

    k0 = kgrd0

    do ith3 = 1, bspline_order
      do m = 1, 10
        tu_p(m) = 0.d0
      end do
      do m = 0, 3
        v(m) = theta3(m, ith3, my_imgs_idx)
      end do
      k0 = k0 + 1
      k = k0 + 1 + (nfft3 - isign(nfft3, k0)) / 2
      
      if ( k .ge. kbot .and. k .le. ktop ) then
        kq = k - kbot0
        j0 = jgrd0
        do ith2 = 1, bspline_order
          j0 = j0 + 1
          j = j0 + 1 + (nfft2 - isign(nfft2, j0)) / 2
          i0 = igrd0
          do m = 0, 3
            u(m) = theta2(m, ith2, my_imgs_idx)
            t_p(m) = 0.d0
          end do
          do ith1 = 1, bspline_order
            i0 = i0 + 1
            i = i0 + 1 + (nfft1 - isign(nfft1, i0)) / 2
            tq_p = Q_p(i, j, kq)
            t_p(0) = t_p(0) + tq_p * theta1(0, ith1, my_imgs_idx)
            t_p(1) = t_p(1) + tq_p * theta1(1, ith1, my_imgs_idx)
            t_p(2) = t_p(2) + tq_p * theta1(2, ith1, my_imgs_idx)
            t_p(3) = t_p(3) + tq_p * theta1(3, ith1, my_imgs_idx)
          end do ! ith1 = 1, bspline_order
          tu_p(Ind_00) = tu_p(Ind_00) + t_p(0) * u(0)
          tu_p(Ind_10) = tu_p(Ind_10) + t_p(1) * u(0)
          tu_p(Ind_01) = tu_p(Ind_01) + t_p(0) * u(1)
          tu_p(Ind_20) = tu_p(Ind_20) + t_p(2) * u(0)
          tu_p(Ind_02) = tu_p(Ind_02) + t_p(0) * u(2)
          tu_p(Ind_11) = tu_p(Ind_11) + t_p(1) * u(1)
          tu_p(Ind_30) = tu_p(Ind_30) + t_p(3) * u(0)
          tu_p(Ind_03) = tu_p(Ind_03) + t_p(0) * u(3)
          tu_p(Ind_21) = tu_p(Ind_21) + t_p(2) * u(1)
          tu_p(Ind_12) = tu_p(Ind_12) + t_p(1) * u(2)
        end do ! ith2 = 1, Spline_order
      end if
      tuv_p(Ind_000) = tuv_p(Ind_000) + tu_p(Ind_00) * v(0)
      tuv_p(Ind_100) = tuv_p(Ind_100) + tu_p(Ind_10) * v(0)
      tuv_p(Ind_010) = tuv_p(Ind_010) + tu_p(Ind_01) * v(0)
      tuv_p(Ind_001) = tuv_p(Ind_001) + tu_p(Ind_00) * v(1)
      tuv_p(Ind_200) = tuv_p(Ind_200) + tu_p(Ind_20) * v(0)
      tuv_p(Ind_020) = tuv_p(Ind_020) + tu_p(Ind_02) * v(0)
      tuv_p(Ind_002) = tuv_p(Ind_002) + tu_p(Ind_00) * v(2)
      tuv_p(Ind_110) = tuv_p(Ind_110) + tu_p(Ind_11) * v(0)
      tuv_p(Ind_101) = tuv_p(Ind_101) + tu_p(Ind_10) * v(1)
      tuv_p(Ind_011) = tuv_p(Ind_011) + tu_p(Ind_01) * v(1)
      tuv_p(Ind_300) = tuv_p(Ind_300) + tu_p(Ind_30) * v(0)
      tuv_p(Ind_030) = tuv_p(Ind_030) + tu_p(Ind_03) * v(0)
      tuv_p(Ind_003) = tuv_p(Ind_003) + tu_p(Ind_00) * v(3)
      tuv_p(Ind_210) = tuv_p(Ind_210) + tu_p(Ind_21) * v(0)
      tuv_p(Ind_201) = tuv_p(Ind_201) + tu_p(Ind_20) * v(1)
      tuv_p(Ind_120) = tuv_p(Ind_120) + tu_p(Ind_12) * v(0)
      tuv_p(Ind_021) = tuv_p(Ind_021) + tu_p(Ind_02) * v(1)
      tuv_p(Ind_102) = tuv_p(Ind_102) + tu_p(Ind_10) * v(2)
      tuv_p(Ind_012) = tuv_p(Ind_012) + tu_p(Ind_01) * v(2)
      tuv_p(Ind_111) = tuv_p(Ind_111) + tu_p(Ind_11) * v(1)
    end do ! ith3 = 1, Spline_order

    atm_id = gbl_img_atm_map(img_id)

    do m = 1, 20
      Fperm_field(m, my_imgs_idx) = tuv_p(m)
    end do

    if (is_polarizable(atm_id)) then
      Fperm_dipfield(1, my_imgs_idx) = tuv_p(Ind_100)
      Fperm_dipfield(2, my_imgs_idx) = tuv_p(Ind_010)
      Fperm_dipfield(3, my_imgs_idx) = tuv_p(Ind_001)
    else
      Fperm_dipfield(1, my_imgs_idx) = 0.d0
      Fperm_dipfield(2, my_imgs_idx) = 0.d0
      Fperm_dipfield(3, my_imgs_idx) = 0.d0
    end if

  end do

  return

end subroutine am_recip_get_perm_f_field

!*******************************************************************************!
! Subroutine:  am_recip_get_2dip_f_fields
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_get_2dip_f_fields(img_cnt, is_polarizable, &
                                      init_grid_ind, dr_order, bspline_order, &
                                      theta1, theta2, theta3,  &
                                      Q1, Q2, Fdip_field1, Fdip_field2, &
                                      my_chgs, my_chg_cnt)

  use mdin_ewald_dat_mod
  use pme_fft_mod
  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: img_cnt
  logical, intent(in)           :: is_polarizable(*)
  integer, intent(in)           :: init_grid_ind(3, *)
  integer, intent(in)           :: dr_order
  integer, intent(in)           :: bspline_order
  double precision, intent(in)  :: theta1(0:dr_order, bspline_order, *)
  double precision, intent(in)  :: theta2(0:dr_order, bspline_order, *)
  double precision, intent(in)  :: theta3(0:dr_order, bspline_order, *)
  double precision, intent(in)  :: Q1(2 * fft_x_dim, fft_y_dim, my_xy_slab_cnt)
  double precision, intent(in)  :: Q2(2 * fft_x_dim, fft_y_dim, my_xy_slab_cnt)
  double precision, intent(out) :: Fdip_field1(3, *)
  double precision, intent(out) :: Fdip_field2(3, *)
  integer, intent(in)           :: my_chgs(*), my_chg_cnt

! Local variables:

  double precision              :: tq1, tq2
  double precision              :: t0_1, t0_2
  double precision              :: t1_1, t1_2
  double precision              :: u0, u1
  double precision              :: v0, v1
  double precision              :: tu00_1, tu00_2
  double precision              :: tu10_1, tu10_2
  double precision              :: tu01_1, tu01_2
  double precision              :: tuv001_1, tuv001_2
  double precision              :: tuv010_1, tuv010_2
  double precision              :: tuv100_1, tuv100_2
  integer                       :: img_id, atm_id, my_imgs_idx
  integer                       :: i, i0, j, j0, k, k0, kq
  integer                       :: kbot0, kbot, ktop
  integer                       :: igrd0, jgrd0, kgrd0
  integer                       :: ith1, ith2, ith3

  kbot0 = my_xy_slab_start
  kbot = kbot0 + 1
  ktop = kbot0 + my_xy_slab_cnt

  do my_imgs_idx = 1, my_chg_cnt
    img_id = my_chgs(my_imgs_idx)

    atm_id = gbl_img_atm_map(img_id)
    if (is_polarizable(atm_id)) then
      igrd0 = init_grid_ind(1, my_imgs_idx) ! begin index in 1st direction
      jgrd0 = init_grid_ind(2, my_imgs_idx) ! begin index in 2nd direction
      kgrd0 = init_grid_ind(3, my_imgs_idx) ! begin index in 3rd direction
      tuv001_1 = 0.d0
      tuv010_1 = 0.d0
      tuv100_1 = 0.d0
      tuv001_2 = 0.d0
      tuv010_2 = 0.d0
      tuv100_2 = 0.d0
      k0 = kgrd0
      do ith3 = 1, bspline_order
        v0 = theta3(0, ith3, my_imgs_idx) ! theta3
        v1 = theta3(1, ith3, my_imgs_idx) ! 1st deriv of theta3
        tu00_1 = 0.d0
        tu10_1 = 0.d0
        tu01_1 = 0.d0
        tu00_2 = 0.d0
        tu10_2 = 0.d0
        tu01_2 = 0.d0
        k0 = k0 + 1
        k = k0 + 1 + (nfft3 - isign(nfft3, k0)) / 2

        if ( k .ge. kbot .and. k .le. ktop ) then
          kq = k - kbot0

          j0 = jgrd0
          do ith2 = 1, bspline_order
            j0 = j0 + 1
            j = j0 + 1 + (nfft2 - isign(nfft2, j0)) / 2
            i0 = igrd0
            u0 = theta2(0, ith2, my_imgs_idx) ! theta2
            u1 = theta2(1, ith2, my_imgs_idx) ! 1st deriv of theta2
            t0_1 = 0.d0
            t1_1 = 0.d0
            t0_2 = 0.d0
            t1_2 = 0.d0
            do ith1 = 1, bspline_order
              i0 = i0 + 1
              i = i0 + 1 + (nfft1 - isign(nfft1, i0)) / 2
              tq1 = Q1(i, j, kq)
              tq2 = Q2(i, j, kq)
              t0_1 = t0_1 + tq1 * theta1(0, ith1, my_imgs_idx)
              t1_1 = t1_1 + tq1 * theta1(1, ith1, my_imgs_idx)
              t0_2 = t0_2 + tq2 * theta1(0, ith1, my_imgs_idx)
              t1_2 = t1_2 + tq2 * theta1(1, ith1, my_imgs_idx)
            end do ! ith1 = 1, bspline_order
            tu00_1 = tu00_1 + t0_1 * u0
            tu10_1 = tu10_1 + t1_1 * u0
            tu01_1 = tu01_1 + t0_1 * u1
            tu00_2 = tu00_2 + t0_2 * u0
            tu10_2 = tu10_2 + t1_2 * u0
            tu01_2 = tu01_2 + t0_2 * u1
          end do ! ith2 = 1, Spline_order

        end if

        tuv100_1 = tuv100_1 + tu10_1 * v0
        tuv010_1 = tuv010_1 + tu01_1 * v0
        tuv001_1 = tuv001_1 + tu00_1 * v1
        tuv100_2 = tuv100_2 + tu10_2 * v0
        tuv010_2 = tuv010_2 + tu01_2 * v0
        tuv001_2 = tuv001_2 + tu00_2 * v1
      end do ! ith3 = 1, Spline_order
      Fdip_field1(1, my_imgs_idx) = tuv100_1
      Fdip_field1(2, my_imgs_idx) = tuv010_1
      Fdip_field1(3, my_imgs_idx) = tuv001_1
      Fdip_field2(1, my_imgs_idx) = tuv100_2
      Fdip_field2(2, my_imgs_idx) = tuv010_2
      Fdip_field2(3, my_imgs_idx) = tuv001_2
    end if ! is_polarizable(atm_id)
  end do ! my_imgs_idx = 1, img_cnt

  return

end subroutine am_recip_get_2dip_f_fields

!*******************************************************************************!
! Subroutine:  am_recip_get_ind_f_fields
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_get_ind_f_fields(img_cnt, init_grid_ind, dr_order, &
                                     bspline_order, theta1, theta2, theta3,  &
                                     Q1, Q2, Fdip_field1, Fdip_field2, &
                                     my_chgs, my_chg_cnt)

  use mdin_ewald_dat_mod
  use pme_fft_mod
  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: img_cnt
  integer, intent(in)           :: init_grid_ind(3, *)
  integer, intent(in)           :: dr_order
  integer, intent(in)           :: bspline_order
  double precision, intent(in)  :: theta1(0:dr_order, bspline_order, *)
  double precision, intent(in)  :: theta2(0:dr_order, bspline_order, *)
  double precision, intent(in)  :: theta3(0:dr_order, bspline_order, *)
  double precision, intent(in)  :: Q1(2 * fft_x_dim, fft_y_dim, my_xy_slab_cnt)
  double precision, intent(in)  :: Q2(2 * fft_x_dim, fft_y_dim, my_xy_slab_cnt)
  double precision, intent(out) :: Fdip_field1(20, *)
  double precision, intent(out) :: Fdip_field2(20, *)
  integer, intent(in)           :: my_chgs(*)
  integer, intent(in)           :: my_chg_cnt

! Note order of field elements
! 1 Field   corresponds to 000 in spline indices
! 2-4 -> F_x, F_y, F_z respectively or to 100, 010, 001 in indices
! 5-10 -> F_xx, F_yy, F_zz, F_xy, F_xz, F_yz resp. or to 
!         200, 020, 002, 110, 101, 011 in indices
! 11-20 -> F_xxx, F_yyy, F_zzz, F_xxy, F_xxz, F_xyy, F_yyz, F_xzz, F_yzz,
!          F_xyz or to
!          300, 030, 003, 210, 201, 120, 021, 102, 012, 111  in indices

! Local variables:

  double precision              :: tq_1, tq_2
  double precision              :: t_1(0:3), t_2(0:3)
  double precision              :: u(0:3), v(0:3)
  double precision              :: tu_1(10), tu_2(10)
  double precision              :: tuv_1(20), tuv_2(20)
  integer                       :: img_id, my_imgs_idx
  integer                       :: m, i, i0, j, j0, k, k0, kq
  integer                       :: kbot0, kbot, ktop
  integer                       :: igrd0, jgrd0, kgrd0
  integer                       :: ith1, ith2, ith3

  kbot0 = my_xy_slab_start
  kbot = kbot0 + 1
  ktop = kbot0 + my_xy_slab_cnt

  do my_imgs_idx = 1, my_chg_cnt
    img_id = my_chgs(my_imgs_idx)

    igrd0 = init_grid_ind(1, my_imgs_idx) ! begin index in 1st direction
    jgrd0 = init_grid_ind(2, my_imgs_idx) ! begin index in 2nd direction
    kgrd0 = init_grid_ind(3, my_imgs_idx) ! begin index in 3rd direction
    do m = 1, 20
      tuv_1(m) = 0.d0
      tuv_2(m) = 0.d0
    end do
    k0 = kgrd0
    do ith3 = 1, bspline_order
      do m = 1, 10
        tu_1(m) = 0.d0
        tu_2(m) = 0.d0
      end do
      do m = 0, 3
        v(m) = theta3(m, ith3, my_imgs_idx)
      end do
      k0 = k0 + 1
      k = k0 + 1 + (nfft3 - isign(nfft3, k0)) / 2

      if ( k .ge. kbot .and. k .le. ktop ) then
        kq = k - kbot0

        j0 = jgrd0
        do ith2 = 1, bspline_order
          j0 = j0 + 1
          j = j0 + 1 + (nfft2 - isign(nfft2, j0)) / 2
          i0 = igrd0
          do m = 0, 3
            u(m) = theta2(m, ith2, my_imgs_idx)
            t_1(m) = 0.d0
            t_2(m) = 0.d0
          end do
          do ith1 = 1, bspline_order
            i0 = i0 + 1
            i = i0 + 1 + (nfft1 - isign(nfft1, i0)) / 2
            tq_1 = Q1(i, j, kq)
            tq_2 = Q2(i, j, kq)
            t_1(0) = t_1(0) + tq_1 * theta1(0, ith1, my_imgs_idx)
            t_1(1) = t_1(1) + tq_1 * theta1(1, ith1, my_imgs_idx)
            t_1(2) = t_1(2) + tq_1 * theta1(2, ith1, my_imgs_idx)
            t_1(3) = t_1(3) + tq_1 * theta1(3, ith1, my_imgs_idx)
            t_2(0) = t_2(0) + tq_2 * theta1(0, ith1, my_imgs_idx)
            t_2(1) = t_2(1) + tq_2 * theta1(1, ith1, my_imgs_idx)
            t_2(2) = t_2(2) + tq_2 * theta1(2, ith1, my_imgs_idx)
            t_2(3) = t_2(3) + tq_2 * theta1(3, ith1, my_imgs_idx)
          end do ! ith1 = 1, bspline_order
          tu_1(Ind_00) = tu_1(Ind_00) + t_1(0) * u(0)
          tu_1(Ind_10) = tu_1(Ind_10) + t_1(1) * u(0)
          tu_1(Ind_01) = tu_1(Ind_01) + t_1(0) * u(1)
          tu_1(Ind_20) = tu_1(Ind_20) + t_1(2) * u(0)
          tu_1(Ind_02) = tu_1(Ind_02) + t_1(0) * u(2)
          tu_1(Ind_11) = tu_1(Ind_11) + t_1(1) * u(1)
          tu_1(Ind_30) = tu_1(Ind_30) + t_1(3) * u(0)
          tu_1(Ind_03) = tu_1(Ind_03) + t_1(0) * u(3)
          tu_1(Ind_21) = tu_1(Ind_21) + t_1(2) * u(1)
          tu_1(Ind_12) = tu_1(Ind_12) + t_1(1) * u(2)
 
          tu_2(Ind_00) = tu_2(Ind_00) + t_2(0) * u(0)
          tu_2(Ind_10) = tu_2(Ind_10) + t_2(1) * u(0)
          tu_2(Ind_01) = tu_2(Ind_01) + t_2(0) * u(1)
          tu_2(Ind_20) = tu_2(Ind_20) + t_2(2) * u(0)
          tu_2(Ind_02) = tu_2(Ind_02) + t_2(0) * u(2)
          tu_2(Ind_11) = tu_2(Ind_11) + t_2(1) * u(1)
          tu_2(Ind_30) = tu_2(Ind_30) + t_2(3) * u(0)
          tu_2(Ind_03) = tu_2(Ind_03) + t_2(0) * u(3)
          tu_2(Ind_21) = tu_2(Ind_21) + t_2(2) * u(1)
          tu_2(Ind_12) = tu_2(Ind_12) + t_2(1) * u(2)
 
        end do ! ith2 = 1, Spline_order
        
      end if

      tuv_1(Ind_000) = tuv_1(Ind_000) + tu_1(Ind_00) * v(0)
      tuv_1(Ind_100) = tuv_1(Ind_100) + tu_1(Ind_10) * v(0)
      tuv_1(Ind_010) = tuv_1(Ind_010) + tu_1(Ind_01) * v(0)
      tuv_1(Ind_001) = tuv_1(Ind_001) + tu_1(Ind_00) * v(1)
      tuv_1(Ind_200) = tuv_1(Ind_200) + tu_1(Ind_20) * v(0)
      tuv_1(Ind_020) = tuv_1(Ind_020) + tu_1(Ind_02) * v(0)
      tuv_1(Ind_002) = tuv_1(Ind_002) + tu_1(Ind_00) * v(2)
      tuv_1(Ind_110) = tuv_1(Ind_110) + tu_1(Ind_11) * v(0)
      tuv_1(Ind_101) = tuv_1(Ind_101) + tu_1(Ind_10) * v(1)
      tuv_1(Ind_011) = tuv_1(Ind_011) + tu_1(Ind_01) * v(1)
      tuv_1(Ind_300) = tuv_1(Ind_300) + tu_1(Ind_30) * v(0)
      tuv_1(Ind_030) = tuv_1(Ind_030) + tu_1(Ind_03) * v(0)
      tuv_1(Ind_003) = tuv_1(Ind_003) + tu_1(Ind_00) * v(3)
      tuv_1(Ind_210) = tuv_1(Ind_210) + tu_1(Ind_21) * v(0)
      tuv_1(Ind_201) = tuv_1(Ind_201) + tu_1(Ind_20) * v(1)
      tuv_1(Ind_120) = tuv_1(Ind_120) + tu_1(Ind_12) * v(0)
      tuv_1(Ind_021) = tuv_1(Ind_021) + tu_1(Ind_02) * v(1)
      tuv_1(Ind_102) = tuv_1(Ind_102) + tu_1(Ind_10) * v(2)
      tuv_1(Ind_012) = tuv_1(Ind_012) + tu_1(Ind_01) * v(2)
      tuv_1(Ind_111) = tuv_1(Ind_111) + tu_1(Ind_11) * v(1)

      tuv_2(Ind_000) = tuv_2(Ind_000) + tu_2(Ind_00) * v(0)
      tuv_2(Ind_100) = tuv_2(Ind_100) + tu_2(Ind_10) * v(0)
      tuv_2(Ind_010) = tuv_2(Ind_010) + tu_2(Ind_01) * v(0)
      tuv_2(Ind_001) = tuv_2(Ind_001) + tu_2(Ind_00) * v(1)
      tuv_2(Ind_200) = tuv_2(Ind_200) + tu_2(Ind_20) * v(0)
      tuv_2(Ind_020) = tuv_2(Ind_020) + tu_2(Ind_02) * v(0)
      tuv_2(Ind_002) = tuv_2(Ind_002) + tu_2(Ind_00) * v(2)
      tuv_2(Ind_110) = tuv_2(Ind_110) + tu_2(Ind_11) * v(0)
      tuv_2(Ind_101) = tuv_2(Ind_101) + tu_2(Ind_10) * v(1)
      tuv_2(Ind_011) = tuv_2(Ind_011) + tu_2(Ind_01) * v(1)
      tuv_2(Ind_300) = tuv_2(Ind_300) + tu_2(Ind_30) * v(0)
      tuv_2(Ind_030) = tuv_2(Ind_030) + tu_2(Ind_03) * v(0)
      tuv_2(Ind_003) = tuv_2(Ind_003) + tu_2(Ind_00) * v(3)
      tuv_2(Ind_210) = tuv_2(Ind_210) + tu_2(Ind_21) * v(0)
      tuv_2(Ind_201) = tuv_2(Ind_201) + tu_2(Ind_20) * v(1)
      tuv_2(Ind_120) = tuv_2(Ind_120) + tu_2(Ind_12) * v(0)
      tuv_2(Ind_021) = tuv_2(Ind_021) + tu_2(Ind_02) * v(1)
      tuv_2(Ind_102) = tuv_2(Ind_102) + tu_2(Ind_10) * v(2)
      tuv_2(Ind_012) = tuv_2(Ind_012) + tu_2(Ind_01) * v(2)
      tuv_2(Ind_111) = tuv_2(Ind_111) + tu_2(Ind_11) * v(1)

    end do ! ith3 = 1, Spline_order

    do m = 1, 20
      Fdip_field1(m, my_imgs_idx) = tuv_1(m)
      Fdip_field2(m, my_imgs_idx) = tuv_2(m)
    end do

  end do

  return

end subroutine am_recip_get_ind_f_fields

!*******************************************************************************!
! Subroutine:  am_recip_perm_ene_grad
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_perm_ene_grad(img_cnt, nfft1, nfft2, nfft3, recip, &
                                  F_field, F_mpole, ene_perm, img_frc, &
                                  my_chgs, my_chg_cnt)

  use amoeba_multipoles_mod, only : coulomb_const_kcal_per_mole
  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: img_cnt
  integer, intent(in)                   :: nfft1, nfft2, nfft3
  double precision, intent(in)          :: recip(3, 3)
  double precision, intent(in)          :: F_field(20, *)
  double precision, intent(in)          :: F_mpole(10, *)
  double precision, intent(out)         :: ene_perm
  double precision, intent(in out)      :: img_frc(3, *)
  integer, intent(in)           :: my_chgs(*), my_chg_cnt

! Local variables:

  double precision                      :: f1, f2, f3, dfx, dfy, dfz
  double precision                      :: mpole_xform_3x3(3, 3)
  double precision                      :: field_xform_3x3(3, 3)
  integer                               :: img_id, my_imgs_idx
  integer                               :: k, j1, j2, j3
 
  call am_recip_xform_matrices(nfft1, nfft2, nfft3, recip, &
                               mpole_xform_3x3, field_xform_3x3)

  do my_imgs_idx = 1, my_chg_cnt
    img_id = my_chgs(my_imgs_idx)

    f1 = 0.d0
    f2 = 0.d0
    f3 = 0.d0
    do k = 1, 10
      ene_perm = ene_perm + F_mpole(k, my_imgs_idx) * F_field(k, my_imgs_idx)
      j1 = deriv1(k) ! indexes for correspoding derivative of the field 
      j2 = deriv2(k)
      j3 = deriv3(k)
      f1 = f1 + F_mpole(k, my_imgs_idx) * F_field(j1, my_imgs_idx)
      f2 = f2 + F_mpole(k, my_imgs_idx) * F_field(j2, my_imgs_idx)
      f3 = f3 + F_mpole(k, my_imgs_idx) * F_field(j3, my_imgs_idx)
    end do

    ! force is negative of gradient
    ! transform from scaled fractional to cartesian--same as fields

    dfx = field_xform_3x3(1, 1) * f1 + field_xform_3x3(1, 2) * f2 + &
          field_xform_3x3(1, 3) * f3
    dfy = field_xform_3x3(2, 1) * f1 + field_xform_3x3(2, 2) * f2 + &
          field_xform_3x3(2, 3) * f3
    dfz = field_xform_3x3(3, 1) * f1 + field_xform_3x3(3, 2) * f2 + &
          field_xform_3x3(3, 3) * f3

    img_frc(1, img_id) = img_frc(1, img_id) - coulomb_const_kcal_per_mole * dfx 
    img_frc(2, img_id) = img_frc(2, img_id) - coulomb_const_kcal_per_mole * dfy 
    img_frc(3, img_id) = img_frc(3, img_id) - coulomb_const_kcal_per_mole * dfz

  end do ! n = 1, img_cnt

  ene_perm = 0.5d0 * coulomb_const_kcal_per_mole * ene_perm

  return

end subroutine am_recip_perm_ene_grad

!*******************************************************************************!
! Subroutine:  am_recip_ind_ene_grad
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_ind_ene_grad(img_cnt, nfft1, nfft2, nfft3, &
                                 is_polarizable, recip, &
                                 F_field, F_mpole, F_dip1_field, F_dip2_field, &
                                 F_dip1, F_dip2, ene_ind, img_frc, &
                                 my_chgs, my_chg_cnt)

  use amoeba_multipoles_mod, only : coulomb_const_kcal_per_mole
  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: img_cnt
  integer, intent(in)                   :: nfft1, nfft2, nfft3
  logical, intent(in)                   :: is_polarizable(*)
  double precision, intent(in)          :: recip(3, 3)
  double precision, intent(in)          :: F_field(20, *)
  double precision, intent(in)          :: F_mpole(10, *)
  double precision, intent(in)          :: F_dip1_field(20, *)  ! multipole field from induced dipoles ind_dip_d 
  double precision, intent(in)          :: F_dip2_field(20, *)  ! multipole field from induced dipoles ind_dip_p
  double precision, intent(in)          :: F_dip1(3, *)         ! induced dipoles (ind_dip_d) in fractional coordinates
  double precision, intent(in)          :: F_dip2(3, *)         ! induced dipoles (ind_dip_p) in fractional coordinates
  double precision, intent(out)         :: ene_ind
  double precision, intent(in out)      :: img_frc(3, *)
  integer, intent(in)                   :: my_chgs(*), my_chg_cnt


! Local variables:

  double precision                      :: f1, f2, f3, dfx, dfy, dfz
  double precision                      :: mpole_xform_3x3(3, 3)
  double precision                      :: field_xform_3x3(3, 3)
  integer                               :: img_id, atm_id, my_imgs_idx
  integer                               :: k, j1, j2, j3
  double precision                      :: eadd
 
  call am_recip_xform_matrices(nfft1, nfft2, nfft3, recip, &
                               mpole_xform_3x3, field_xform_3x3)

  do my_imgs_idx = 1, my_chg_cnt
    img_id = my_chgs(my_imgs_idx)

    atm_id = gbl_img_atm_map(img_id)
    f1 = 0.d0
    f2 = 0.d0
    f3 = 0.d0
    ! add the force on mpoles due to ave of dip fields
    do k = 1, 10
      j1 = deriv1(k)
      j2 = deriv2(k)
      j3 = deriv3(k)
      f1 = f1 + 0.5d0 * F_mpole(k, my_imgs_idx) * &
                        (F_dip1_field(j1, my_imgs_idx) + &
                         F_dip2_field(j1, my_imgs_idx))

      f2 = f2 + 0.5d0 * F_mpole(k, my_imgs_idx) * &
                        (F_dip1_field(j2, my_imgs_idx) + &
                         F_dip2_field(j2, my_imgs_idx))

      f3 = f3 + 0.5d0 * F_mpole(k, my_imgs_idx) * &
                        (F_dip1_field(j3, my_imgs_idx) + &
                         F_dip2_field(j3, my_imgs_idx))
    end do
    
    if (is_polarizable(atm_id)) then

      ! induced energy is interaction of direct dipoles with perm field

      eadd = F_dip1(1, my_imgs_idx) * F_field(Ind_100, my_imgs_idx)
      ene_ind = ene_ind + eadd
    
      eadd = F_dip1(2, my_imgs_idx) * F_field(Ind_010, my_imgs_idx)
      ene_ind = ene_ind + eadd
      
      eadd = F_dip1(3, my_imgs_idx) * F_field(Ind_001, my_imgs_idx)
      ene_ind = ene_ind + eadd

      do k = 2, 4

        j1 = deriv1(k)
        j2 = deriv2(k)
        j3 = deriv3(k)

        ! add the force on ave of dipoles due to perm field

        f1 = f1 + 0.5d0 * (F_dip1(k - 1, my_imgs_idx) + &
                           F_dip2(k - 1, my_imgs_idx)) * &
                           F_field(j1, my_imgs_idx)
        f2 = f2 + 0.5d0 * (F_dip1(k - 1, my_imgs_idx) + &
                           F_dip2(k - 1, my_imgs_idx)) * &
                           F_field(j2, my_imgs_idx)
        f3 = f3 + 0.5d0 * (F_dip1(k - 1, my_imgs_idx) + &
                           F_dip2(k - 1, my_imgs_idx)) * &
                           F_field(j3, my_imgs_idx)

        ! next the forces of dips on each other

        f1 = f1 + 0.5d0 * (F_dip1(k-1, my_imgs_idx) * &
                           F_dip2_field(j1, my_imgs_idx) + &
                           F_dip2(k-1, my_imgs_idx) * &
                           F_dip1_field(j1, my_imgs_idx))
        f2 = f2 + 0.5d0 * (F_dip1(k-1, my_imgs_idx) * &
                           F_dip2_field(j2, my_imgs_idx) + &
                           F_dip2(k-1, my_imgs_idx) * &
                           F_dip1_field(j2, my_imgs_idx))
        f3 = f3 + 0.5d0 * (F_dip1(k-1, my_imgs_idx) * &
                           F_dip2_field(j3, my_imgs_idx) + &
                           F_dip2(k-1, my_imgs_idx) * &
                           F_dip1_field(j3, my_imgs_idx))

      end do
    end if !is_polarizable(atm_id)

! force is negative of gradient
! transform from scaled fractional to cartesian--same as fields

    dfx = field_xform_3x3(1, 1) * f1 + field_xform_3x3(1, 2) * f2 + &
          field_xform_3x3(1, 3) * f3
    dfy = field_xform_3x3(2, 1) * f1 + field_xform_3x3(2, 2) * f2 + &
          field_xform_3x3(2, 3) * f3
    dfz = field_xform_3x3(3, 1) * f1 + field_xform_3x3(3, 2) * f2 + &
          field_xform_3x3(3, 3) * f3

    img_frc(1, img_id) = img_frc(1, img_id) - coulomb_const_kcal_per_mole * dfx 
    img_frc(2, img_id) = img_frc(2, img_id) - coulomb_const_kcal_per_mole * dfy 
    img_frc(3, img_id) = img_frc(3, img_id) - coulomb_const_kcal_per_mole * dfz

  end do !my_imgs_idx = 1, img_cnt

  ene_ind = 0.5d0 * coulomb_const_kcal_per_mole * ene_ind

  return

end subroutine am_recip_ind_ene_grad

!*******************************************************************************!
! Subroutine:  am_recip_global_to_fractional
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_global_to_fractional(img_cnt, recip, &
                                         nfft1, nfft2, nfft3, &
                                         glob_mpole, frac_mpole, &
                                         my_chgs, my_chg_cnt)

  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: img_cnt
  double precision, intent(in)  :: recip(3, 3)
  integer, intent(in)           :: nfft1, nfft2, nfft3
  double precision, intent(in)  :: glob_mpole(10, *)
  double precision, intent(out) :: frac_mpole(10, *)
  integer, intent(in)           :: my_chgs(*), my_chg_cnt

! Local variables:

  integer                       :: order, dimxy
  integer                       :: img_id, atm_id, my_imgs_idx
  double precision              :: mpole_xform_3x3(3, 3), field_xform_3x3(3, 3)
  double precision              :: Mpole_xy(MAXMP * MAXMP)

  order = 10
  dimxy = 10

! first get mpole_xform_3x3

  call am_recip_xform_matrices(nfft1, nfft2, nfft3, recip, &
                               mpole_xform_3x3, field_xform_3x3)

  call xform_mpole_matrix(mpole_xform_3x3, Mpole_xy, order)

  do my_imgs_idx = 1, my_chg_cnt
    
    img_id = my_chgs(my_imgs_idx)
    
    atm_id = gbl_img_atm_map(img_id)
    call xform_mpole(Mpole_xy, dimxy, glob_mpole(:, atm_id), &
                     frac_mpole(:, my_imgs_idx), order)
  end do

  return

end subroutine am_recip_global_to_fractional

!*******************************************************************************!
! Subroutine:  am_recip_field_torque_virial
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_field_torque_virial(img_cnt, is_polarizable, &
                                        nfft1, nfft2, nfft3, recip, &
                                        F_field, F_dip1_field, F_dip2_field, &
                                        torque_field, virial, mpole_p, &
                                        dipole1, dipole2, my_chgs, my_chg_cnt)

  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: img_cnt
  logical, intent(in)                   :: is_polarizable(*)
  integer, intent(in)                   :: nfft1, nfft2, nfft3
  double precision, intent(in)          :: recip(3, 3)
  double precision, intent(in)          :: F_field(20, *)
  double precision, intent(in)          :: F_dip1_field(20, *)
  double precision, intent(in)          :: F_dip2_field(20, *)
  double precision, intent(in out)      :: torque_field(10, *)
  double precision, intent(in out)      :: virial(3, 3)
  double precision, intent(in)          :: mpole_p(10, *)
  double precision, intent(in)          :: dipole1(3, *)
  double precision, intent(in)          :: dipole2(3, *)

  integer, intent(in)                   :: my_chgs(*), my_chg_cnt


! Local variables:

  integer                               :: i, j, k
  integer                               :: img_id, atm_id, my_imgs_idx
  integer                               :: order, dimxy
  double precision                      :: field_p(10)
  double precision                      :: field_d1(10)
  double precision                      :: field_d2(10)
  double precision                      :: mpole_xform_3x3(3, 3)
  double precision                      :: field_xform_3x3(3, 3)
  double precision                      :: Field_xy(MAXMP * MAXMP)

  order = 10
  dimxy = 10

  ! transform fields from fractional to cartesian
  
  call am_recip_xform_matrices(nfft1, nfft2, nfft3, recip, &
                               mpole_xform_3x3, field_xform_3x3)

  ! get the higher order terms

  call xform_mpole_field_matrix(field_xform_3x3, Field_xy, order)

  do my_imgs_idx = 1, my_chg_cnt
    img_id = my_chgs(my_imgs_idx)

    atm_id = gbl_img_atm_map(img_id)

    call xform_field(Field_xy, dimxy, F_field(:, my_imgs_idx), &
                     field_p, order)

    call xform_field(Field_xy, dimxy, F_dip1_field(:, my_imgs_idx), &
                     field_d1, order)

    call xform_field(Field_xy, dimxy, F_dip2_field(:, my_imgs_idx), &
                     field_d2, order)

    ! torque is due to permanent field plus ave of dip fields

    do i = 1, 10
      torque_field(i, atm_id) = torque_field(i, atm_id) +  &
                           field_p(i) + 0.5d0 * (field_d1(i) + field_d2(i))
    end do

    ! extra virial terms due to dipolar contributions

    do j = 1, 3
      do k = 1, 3
        virial(j, k) = virial(j, k) - &
          (torque_field(j + 1, atm_id) * mpole_p(k + 1, atm_id) + &
          0.5d0 * field_p(j + 1) * (dipole1(k, atm_id) &
          + dipole2(k, atm_id)) + &
          0.5d0 * (field_d1(j + 1) * dipole2(k, atm_id) + &
          field_d2(j + 1) * dipole1(k, atm_id)))
      end do
    end do

    ! quadrupolar contribs

    virial(1, 1) = virial(1, 1) - &
                   (2.d0 * mpole_p(5, atm_id) * torque_field(5, atm_id) + &
                    mpole_p(8, atm_id) * torque_field(8, atm_id) + &
                    mpole_p(9, atm_id) * torque_field(9, atm_id))

    virial(2, 2) = virial(2, 2) - &
                   (2.d0 * mpole_p(6, atm_id) * torque_field(6, atm_id) + &
                    mpole_p(8, atm_id) * torque_field(8, atm_id) + &
                    mpole_p(10,  atm_id) * torque_field(10, atm_id))

    virial(3, 3) = virial(3, 3) - &
                   (2.d0 * mpole_p(7, atm_id) * torque_field(7, atm_id) + &
                    mpole_p(9, atm_id) * torque_field(9, atm_id) + &
                    mpole_p(10, atm_id) * torque_field(10, atm_id))

    virial(1, 2) = virial(1, 2) - &
                   (2.d0 * mpole_p(6, atm_id) * torque_field(8, atm_id) + &
                    mpole_p(8, atm_id) * torque_field(5, atm_id) + &
                    mpole_p(10, atm_id) * torque_field(9, atm_id))

    virial(2, 1) = virial(2, 1) - &
                   (2.d0 * mpole_p(5, atm_id) * torque_field(8, atm_id) + &
                    mpole_p(8, atm_id) * torque_field(6, atm_id) + &
                    mpole_p(9, atm_id) * torque_field(10, atm_id))

    virial(1, 3) = virial(1, 3) - &
                   (2.d0 * mpole_p(7, atm_id) * torque_field(9, atm_id) + &
                    mpole_p(9, atm_id) * torque_field(5, atm_id) + &
                    mpole_p(10, atm_id) * torque_field(8, atm_id))

    virial(3, 1) = virial(3, 1) - &
                   (2.d0 * mpole_p(5, atm_id) * torque_field(9, atm_id) + &
                    mpole_p(8, atm_id) * torque_field(10, atm_id) + &
                    mpole_p(9, atm_id) * torque_field(7, atm_id))
           
    virial(2, 3) = virial(2, 3) - &
                   (2.d0 * mpole_p(7, atm_id) * torque_field(10, atm_id) + &
                    mpole_p(9, atm_id) * torque_field(8, atm_id) + &
                    mpole_p(10, atm_id) * torque_field(6, atm_id))
           
    virial(3, 2) = virial(3, 2) - &
                   (2.d0 * mpole_p(6, atm_id) * torque_field(10, atm_id) + &
                    mpole_p(8, atm_id) * torque_field(9, atm_id) + &
                    mpole_p(10, atm_id) * torque_field(7, atm_id))

  end do

  !symmetrize

  virial(1, 2) = 0.5d0 * (virial(1, 2) + virial(2, 1))
  virial(1, 3) = 0.5d0 * (virial(1, 3) + virial(3, 1))
  virial(2, 3) = 0.5d0 * (virial(2, 3) + virial(3, 2))
  virial(2, 1) = virial(1, 2)
  virial(3, 1) = virial(1, 3)
  virial(3, 2) = virial(2, 3)
  
  return

end subroutine am_recip_field_torque_virial

!*******************************************************************************!
! Subroutine:  am_recip_fdip_to_cdip_fiel
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_fdip_to_cdip_field(img_cnt, is_polarizable, &
                                       nfft1, nfft2, nfft3, recip, &
                                       frac_dipole_field, cart_dipole_field, &
                                       my_chgs, my_chg_cnt)

  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: img_cnt
  logical, intent(in)           :: is_polarizable(*)
  integer, intent(in)           :: nfft1, nfft2, nfft3
  double precision, intent(in)  :: recip(3, 3)
  double precision, intent(in)  :: frac_dipole_field(3, *)
  double precision, intent(out) :: cart_dipole_field(3, *)
  integer, intent(in)           :: my_chgs(*), my_chg_cnt

! Local variables:

  integer                       :: img_id, atm_id
  integer                       :: my_imgs_idx
  integer                       :: i, j
  double precision              :: mpole_xform_3x3(3, 3), field_xform_3x3(3, 3)

  call am_recip_xform_matrices(nfft1, nfft2, nfft3, recip, &
                               mpole_xform_3x3, field_xform_3x3)
  do my_imgs_idx = 1, my_chg_cnt
    img_id = my_chgs(my_imgs_idx)

    atm_id = gbl_img_atm_map(img_id)
    if (is_polarizable(atm_id)) then
      do i = 1, 3
        ! add recip contribution to cart_dipole_field
        do j = 1, 3
          cart_dipole_field(i, atm_id) = cart_dipole_field(i, atm_id) + &
            field_xform_3x3(i, j) * frac_dipole_field(j, my_imgs_idx)
        end do
      end do
    end if
  end do

  return

end subroutine am_recip_fdip_to_cdip_field

!*******************************************************************************!
! Subroutine:  am_recip_dip_to_frac_dip
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_dip_to_frac_dip(img_cnt, is_polarizable, &
                                    ind_dip1, ind_dip2, f_dipole1, f_dipole2, &
                                    my_chgs, my_chg_cnt)

  use img_mod
  use pbc_mod
  use mdin_ewald_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: img_cnt
  logical, intent(in)           :: is_polarizable(*)
  double precision, intent(in)  :: ind_dip1(3, *)
  double precision, intent(in)  :: ind_dip2(3, *)
  double precision, intent(out) :: f_dipole1(3, *)
  double precision, intent(out) :: f_dipole2(3, *)
  integer, intent(in)           :: my_chgs(*), my_chg_cnt

! Local variables:

  integer                       :: img_id, atm_id
  integer                       :: my_imgs_idx
  integer                       :: i, j
  double precision              :: mpole_xform_3x3(3, 3), field_xform_3x3(3, 3)

  call am_recip_xform_matrices(nfft1, nfft2, nfft3, recip, &
                               mpole_xform_3x3, field_xform_3x3)

  do my_imgs_idx = 1, my_chg_cnt
    img_id = my_chgs(my_imgs_idx)

    atm_id = gbl_img_atm_map(img_id)
    if (is_polarizable(atm_id)) then
      do i = 1, 3
        f_dipole1(i, my_imgs_idx) = 0.d0
        f_dipole2(i, my_imgs_idx) = 0.d0
        do j = 1, 3
          f_dipole1(i, my_imgs_idx) = f_dipole1(i, my_imgs_idx) + &
                            mpole_xform_3x3(i, j) * ind_dip1(j, atm_id)
          f_dipole2(i, my_imgs_idx) = f_dipole2(i, my_imgs_idx) + &
                            mpole_xform_3x3(i, j) * ind_dip2(j, atm_id)
        end do
      end do
    end if
  end do

  return

end subroutine am_recip_dip_to_frac_dip

!*******************************************************************************!
! Subroutine:  am_recip_xform_matrices
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_xform_matrices(nfft1, nfft2, nfft3, recip, &
                                   mpole_xform_3x3, field_xform_3x3)

  implicit none

! Formal arguments:

  integer, intent(in)           :: nfft1, nfft2, nfft3
  double precision, intent(in)  :: recip(3, 3)
  double precision, intent(out) :: mpole_xform_3x3(3, 3)  
  double precision, intent(out) :: field_xform_3x3(3, 3)

! Local variables:

  double precision              :: du1_dx, du1_dy, du1_dz
  double precision              :: du2_dx, du2_dy, du2_dz
  double precision              :: du3_dx, du3_dy, du3_dz

  du1_dx = nfft1 * recip(1, 1)  !  = 1/(a/nfft1) for orthog  = 1/dx 
  du1_dy = nfft1 * recip(2, 1)
  du1_dz = nfft1 * recip(3, 1)
  du2_dx = nfft2 * recip(1, 2) 
  du2_dy = nfft2 * recip(2, 2)  ! = 1/(b/nfft2) for orthog  = 1/dy
  du2_dz = nfft2 * recip(3, 2)
  du3_dx = nfft3 * recip(1, 3)
  du3_dy = nfft3 * recip(2, 3)
  du3_dz = nfft3 * recip(3, 3)  ! = 1/(c/nfft3) for orthog  = 1/dz

  field_xform_3x3(1, 1) = du1_dx
  field_xform_3x3(1, 2) = du2_dx
  field_xform_3x3(1, 3) = du3_dx
  field_xform_3x3(2, 1) = du1_dy
  field_xform_3x3(2, 2) = du2_dy
  field_xform_3x3(2, 3) = du3_dy
  field_xform_3x3(3, 1) = du1_dz
  field_xform_3x3(3, 2) = du2_dz
  field_xform_3x3(3, 3) = du3_dz

  mpole_xform_3x3(1, 1) = du1_dx
  mpole_xform_3x3(1, 2) = du1_dy
  mpole_xform_3x3(1, 3) = du1_dz
  mpole_xform_3x3(2, 1) = du2_dx
  mpole_xform_3x3(2, 2) = du2_dy
  mpole_xform_3x3(2, 3) = du2_dz
  mpole_xform_3x3(3, 1) = du3_dx
  mpole_xform_3x3(3, 2) = du3_dy
  mpole_xform_3x3(3, 3) = du3_dz

  return

end subroutine am_recip_xform_matrices

!*******************************************************************************!
! Subroutine:  am_test_frac_cart_ene
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_test_frac_cart_ene(img_cnt, nfft1, nfft2, nfft3, recip, &
                                 F_mpole, F_field, C_mpole, my_chgs, my_chg_cnt)

  use amoeba_multipoles_mod, only : coulomb_const_kcal_per_mole
  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: img_cnt
  integer, intent(in)           :: nfft1, nfft2, nfft3
  double precision, intent(in)  :: recip(3, 3)
  double precision, intent(in)  :: F_mpole(10, *)
  double precision, intent(in)  :: F_field(20, *)
  double precision, intent(in)  :: C_mpole(10, *)
  integer, intent(in)           :: my_chgs(*), my_chg_cnt


! Local variables:

  double precision              :: ene_F, ene_C
  double precision              :: C_field(10)
  double precision              :: mpole_xform_3x3(3, 3)
  double precision              :: field_xform_3x3(3, 3)
  double precision              :: Field_xy(MAXMP * MAXMP)
  integer                       :: k, dimxy, order
  integer                       :: img_id, atm_id, my_imgs_idx

  ene_F = 0.d0
  ene_C = 0.d0

  call am_recip_xform_matrices(nfft1, nfft2, nfft3, recip, &
                               mpole_xform_3x3, field_xform_3x3)

  call xform_mpole_field_matrix(field_xform_3x3, Field_xy, order)

  do my_imgs_idx = 1, my_chg_cnt
    img_id = my_chgs(my_imgs_idx)

  atm_id = gbl_img_atm_map(img_id)
    call xform_field(Field_xy, dimxy, F_field(:, my_imgs_idx), C_field, order)
    do k = 1, 10
      ene_F = ene_F + F_mpole(k, my_imgs_idx) * F_field(k, my_imgs_idx)
      ene_C = ene_C + C_mpole(k, atm_id) * C_field(k)
    end do
  end do

  ene_F = 0.5d0 * coulomb_const_kcal_per_mole * ene_F
  ene_C = 0.5d0 * coulomb_const_kcal_per_mole * ene_C

  write(6, '(a, 2(1x, f14.4))') 'am_test_frac_cart:  ene_F, ene_C, = ', &
                                ene_F, ene_C

  return

end subroutine am_test_frac_cart_ene

!*******************************************************************************!
! Subroutine:   am_recip_calc_netfrcs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_recip_calc_netfrcs(img_cnt, img_frc, netfrcs, my_chgs, my_chg_cnt)

  use parallel_dat_mod
  use mdin_ewald_dat_mod, only : netfrc
  use timers_mod
  use img_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: img_cnt
  double precision, intent(in out)      :: img_frc(3, *)
  double precision, intent(out)         :: netfrcs(3)
  integer, intent(in)                   :: my_chgs(*), my_chg_cnt

! Local variables:

  integer                               :: img_id, my_imgs_idx

  netfrcs(:) = 0.d0

  if (netfrc .eq. 0) return

  do my_imgs_idx = 1, my_chg_cnt
    img_id = my_chgs(my_imgs_idx)
    netfrcs(:) = netfrcs(:) + img_frc(:, img_id)
  end do

  call update_pme_time(pme_misc_timer)

  return

end subroutine am_recip_calc_netfrcs

end module amoeba_recip_mod
