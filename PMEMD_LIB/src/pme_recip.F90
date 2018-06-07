#include "copyright.i"

!*******************************************************************************
!
! Module: pme_recip_mod
!
! Description: <TBS>
!
! NOTE NOTE NOTE: This code assumes frc_int .eq. 0 and should only be used
!                 under these conditions!!!
!*******************************************************************************

module pme_recip_mod

  implicit none

  double precision, allocatable, save, private  :: m1_exp_tbl(:)
  double precision, allocatable, save, private  :: m2_exp_tbl(:)
  double precision, allocatable, save, private  :: m3_exp_tbl(:)

  double precision, allocatable, save           :: gbl_prefac1(:)
  double precision, allocatable, save           :: gbl_prefac2(:)
  double precision, allocatable, save           :: gbl_prefac3(:)

  integer, save                       :: max_recip_imgs
  integer, save                       :: recip_img_lo, recip_img_hi
  logical, save                       :: recip_img_range_wraps

contains

!*******************************************************************************
!
! Subroutine:  do_pmesh_kspace
!
! Description:  <TBS>
!
! INPUT:
!
! OUTPUT:
!
! img_frc:      Forces incremented by k-space sum.
! virial:       Virial due to k-space sum (valid for atomic scaling;
!               rigid molecule virial needs a correction term not computed here.
!*******************************************************************************

subroutine do_pmesh_kspace(crd, img, img_frc, eer, virial, frcx, frcy, frcz)

  use img_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pbc_mod
  use pme_fft_mod
  use prmtop_dat_mod
  use timers_mod

  implicit none

! Formal arguments:

  double precision      :: crd(3,*)
  type(img_rec)         :: img(*)
  double precision      :: img_frc(3, *)
  double precision      :: eer
  double precision      :: virial(3, 3)
  double precision      :: frcx, frcy, frcz

! Local variables:

! NOTE - For FFTW, *_q need to be 16 byte aligned if they are moved off
!        the stack (and the stack must be aligned if they are on the stack).

  integer               :: my_chg_cnt  ! Used only in MPI

  double precision      :: theta(bspl_order * 3 * max_recip_imgs)
  double precision      :: dtheta(bspl_order * 3 * max_recip_imgs)

  double precision      :: xyz_q(max(xy_slab_dbl_cnt * max_xy_slab_cnt, &
                                     zx_slab_dbl_cnt * max_zx_slab_cnt))
  double precision      :: zxy_q(max(xy_slab_dbl_cnt * max_xy_slab_cnt, &
                                     zx_slab_dbl_cnt * max_zx_slab_cnt))
  integer               :: my_chgs(max_recip_imgs)  ! used only for MPI

!  double precision      :: xyz_q(2 * fft_x_dim * fft_y_dim * fft_z_dim)  ! not MPI - reduced to this for numtasks == 1
!  double precision      :: zxy_q(2 * fft_x_dim * fft_y_dim * fft_z_dim)

  integer               :: ifracts(3, max_recip_imgs)

  ! Exponential tables are only practical to use for orthogonal unit cells, and
  ! only need updating when constant pressure runs are being done.

  if (ntp .gt. 0 .and. is_orthog .ne. 0) call load_m_exp_tbls

  if (is_orthog .ne. 0) then
      call get_grid_weights(natom, img, ifracts, theta, dtheta, bspl_order, &
                            gbl_img_atm_map, my_chgs, my_chg_cnt)
  else
      call get_grid_weights_nonorthog(natom, crd, gbl_img_atm_map, ifracts, &
                                      theta, dtheta, bspl_order, my_chgs, &
                                      my_chg_cnt)
  end if
 
  call update_pme_time(bspline_timer)

! Fill Charge Grid.  Charges are approximated on an even grid.

  call fill_charge_grid(xyz_q, theta, bspl_order, img, ifracts, my_chgs, &
                        my_chg_cnt)

  call update_pme_time(grid_charges_timer)

  call fft3drc_forward(xyz_q, zxy_q, fft_x_dim, fft_y_dim, fft_z_dim)

  call update_pme_time(fft_timer)

  if (is_orthog .ne. 0) then
    call scalar_sumrc(zxy_q, ew_coeff, uc_volume, &
                      gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                      nfft1, nfft2, nfft3, &
                      fft_x_dim, fft_y_dim, fft_z_dim, eer, virial)
  else
    call scalar_sumrc_nonorthog(zxy_q, ew_coeff, uc_volume, &
                                gbl_prefac1, gbl_prefac2, gbl_prefac3, &
                                nfft1, nfft2, nfft3, &
                                fft_x_dim, fft_y_dim, fft_z_dim, &
                                eer, virial)
  end if

  call update_pme_time(scalar_sum_timer)

  call fft3drc_back(zxy_q, xyz_q, fft_x_dim, fft_y_dim, fft_z_dim)

  call update_pme_time(fft_timer)

  call grad_sum(xyz_q, theta, dtheta, bspl_order, &
                img, img_frc, ifracts, my_chgs, my_chg_cnt, &
                fft_x_dim, fft_y_dim, fft_z_dim, &
                frcx, frcy, frcz)

  call update_pme_time(grad_sum_timer)

  return

end subroutine do_pmesh_kspace

!*******************************************************************************
!
! Subroutine:  fill_charge_grid
!
! Description: <TBS>
!
! INPUT:
!
! theta1, theta2, theta3:       Spline coeff arrays.
! ifracts:                      int(scaled fractional coords).
! nfft1, nfft2, nfft3:          Logical charge grid dimensions.
!
! fft_x_dim, fft_y_dim, fft_z_dim: Physical charge grid dims.
!
! order:                    Order of spline interpolation.
!
! OUTPUT:
!
! q:                            Charge grid.
!              
!*******************************************************************************

subroutine fill_charge_grid(q, theta, order, img, ifracts, my_chgs, my_chg_cnt)

  use img_mod
  use mdin_ewald_dat_mod
  use pme_fft_mod
  use prmtop_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: order
  double precision      :: theta(order, 3, *)
  type(img_rec)         :: img(*)
  integer               :: ifracts(3, *)
  double precision      :: q(2 * fft_x_dim, fft_y_dim, my_xy_slab_cnt)

  integer               :: my_chgs(*), my_chg_cnt

! Local variables:

  integer               :: i, j, k
  integer               :: i_base, j_base, k_base
  integer               :: ith1, ith2, ith3
  integer               :: my_chg_idx
  double precision      :: charge
  double precision      :: k_term, j_term

! integer               :: kbot0
  integer               :: kbot, ktop

  integer               :: i_tbl(-order : nfft1)
  integer               :: j_tbl(-order : nfft2)
  integer               :: k_tbl(-order : nfft3)

#ifdef RANGE_DBG
  integer               :: dbg_img_map(natom)
  integer               :: dbg_img_map_res(natom)

  dbg_img_map(:) = 0
#endif /* RANGE_DBG */

! Old usages; kept for documentation during code conversion...
! kbot0 = my_xy_slab_start
! kbot = kbot0 + 1
! ktop = kbot0 + my_xy_slab_cnt

! Zero the Charge grids:

  call zero_q_array(q, 2 * fft_x_dim * fft_y_dim * my_xy_slab_cnt)

  ! Initialize the indexing tables.  It actually produces faster code to
  ! do this here rather than caching the info remotely.

  do i = 0, nfft1
    i_tbl(i) = i + 1
  end do
  do i = -order, -1
    i_tbl(i) = i + nfft1 + 1
  end do

  do j = 0, nfft2
    j_tbl(j) = j + 1
  end do
  do j = -order, -1
    j_tbl(j) = j + nfft2 + 1
  end do

  if(numtasks .gt. 1) then
      kbot = 1
      ktop = my_xy_slab_cnt
      do k = 0, nfft3
        k_tbl(k) = k - my_xy_slab_start + 1
        if (k_tbl(k) .lt. kbot .or. k_tbl(k) .gt. ktop) k_tbl(k) = 0
      end do
      do k = -order, -1
        k_tbl(k) = k + nfft3 - my_xy_slab_start + 1
        if (k_tbl(k) .lt. kbot .or. k_tbl(k) .gt. ktop) k_tbl(k) = 0
      end do
  else
      do k = 0, nfft3
        k_tbl(k) = k + 1
      end do
      do k = -order, -1
        k_tbl(k) = k + nfft3 + 1
      end do
  endif

  ! We special-case order 4, the default.
  ! We special-case order 6, for future amoeba use.

  if (order .eq. 4) then
    do my_chg_idx = 1, my_chg_cnt
      charge = img(my_chgs(my_chg_idx))%qterm

      i_base = ifracts(1, my_chg_idx)
      j_base = ifracts(2, my_chg_idx)
      k_base = ifracts(3, my_chg_idx)

      do ith3 = 1, 4

        k = k_tbl(k_base + ith3)

        if (k .ne. 0) then

#ifdef RANGE_DBG
    if(numtasks .gt. 1) then
          dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
    endif
#endif /* RANGE_DBG */
          k_term =  theta(ith3, 3, my_chg_idx) * charge
          do ith2 = 1, 4
            j = j_tbl(j_base + ith2)
            j_term = k_term * theta(ith2, 2, my_chg_idx)
            q(i_tbl(i_base+1), j, k) = q(i_tbl(i_base+1), j, k) + &
                                       theta(1, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+2), j, k) = q(i_tbl(i_base+2), j, k) + &
                                       theta(2, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+3), j, k) = q(i_tbl(i_base+3), j, k) + &
                                       theta(3, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+4), j, k) = q(i_tbl(i_base+4), j, k) + &
                                       theta(4, 1, my_chg_idx) * j_term
          end do !  (ith2)
        end if !  (k != 0)
      end do !  (ith3)  
    end do !  (my_chg_idx)
  else if (order .eq. 6) then
    do my_chg_idx = 1, my_chg_cnt
      charge = img(my_chgs(my_chg_idx))%qterm

      i_base = ifracts(1, my_chg_idx)
      j_base = ifracts(2, my_chg_idx)
      k_base = ifracts(3, my_chg_idx)

      do ith3 = 1, 6
        k = k_tbl(k_base + ith3)

        if (k .ne. 0) then
#ifdef RANGE_DBG
        if( numtasks .gt. 1) then
          dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
        endif
#endif /* RANGE_DBG */
          k_term =  theta(ith3, 3, my_chg_idx) * charge
          do ith2 = 1, 6
            j = j_tbl(j_base + ith2)
            j_term = k_term * theta(ith2, 2, my_chg_idx)
            q(i_tbl(i_base+1), j, k) = q(i_tbl(i_base+1), j, k) + &
                                       theta(1, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+2), j, k) = q(i_tbl(i_base+2), j, k) + &
                                       theta(2, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+3), j, k) = q(i_tbl(i_base+3), j, k) + &
                                       theta(3, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+4), j, k) = q(i_tbl(i_base+4), j, k) + &
                                       theta(4, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+5), j, k) = q(i_tbl(i_base+5), j, k) + &
                                       theta(5, 1, my_chg_idx) * j_term
            q(i_tbl(i_base+6), j, k) = q(i_tbl(i_base+6), j, k) + &
                                       theta(6, 1, my_chg_idx) * j_term
          end do ! (ith2)
        end if ! ( k != 0)
      end do ! (ith3)
    end do ! (my_chg_idx)
  else ! (order != 4, != 6)
    do my_chg_idx = 1, my_chg_cnt
      charge = img(my_chgs(my_chg_idx))%qterm

      i_base = ifracts(1, my_chg_idx)
      j_base = ifracts(2, my_chg_idx)
      k_base = ifracts(3, my_chg_idx)

      do ith3 = 1, order
      
        k = k_tbl(k_base + ith3)
        
        if (k .ne. 0) then

#ifdef RANGE_DBG
        if(numtasks .gt. 1) then
          dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
        endif
#endif /* RANGE_DBG */
          k_term =  theta(ith3, 3, my_chg_idx) * charge
          do ith2 = 1, order
            j = j_tbl(j_base + ith2)
            j_term = k_term * theta(ith2, 2, my_chg_idx)
            do ith1 = 1, order
              i = i_tbl(i_base + ith1)
              q(i, j, k) = q(i, j, k) + theta(ith1, 1, my_chg_idx) * j_term
            end do ! (ith1)
          end do  ! (ith2 )
        end if  ! ( k != 0) 
      end do  ! (ith3)
    end do ! (my_chg_idx)
  end if ! ( order == 4)

#ifdef RANGE_DBG
  if( numtasks .gt. 1) then
    call mpi_reduce(dbg_img_map, dbg_img_map_res, size(dbg_img_map), &
                      mpi_integer, mpi_sum, 0, lib_mpi_comm, err_code_mpi)
    if (master) then
        do i = 1, natom
            if (dbg_img_map_res(i) .ne. order) then
               write(0,*)'DBG: IMAGE ', i, 'processed ', &
                      dbg_img_map_res(i), 'times in fill_charge_grid code!!!'
            end if
        end do
    end if
  endif
#endif /* RANGE_DBG */

  return

end subroutine fill_charge_grid

!*******************************************************************************
!
! Subroutine:  zero_q_array
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine zero_q_array(array, num)

  implicit none

  double precision      array(*)
  integer               num

  integer               i

  do i = 1, num
    array(i) = 0.d0
  end do

  return

end subroutine zero_q_array

!*******************************************************************************
!
! Subroutine:  grad_sum
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine grad_sum(q, theta, dtheta, order, img, img_frc, ifracts, &
                    my_chgs, my_chg_cnt, x_dim, y_dim, z_dim, &
                    frcx, frcy, frcz)

  use img_mod
  use mdin_ewald_dat_mod
  use prmtop_dat_mod
  use pbc_mod
  use pme_fft_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: x_dim, y_dim, z_dim
  double precision      :: q(2 * x_dim, y_dim, z_dim)
  integer               :: order
  double precision      :: theta(order, 3, *)
  double precision      :: dtheta(order, 3, *)
  type(img_rec)         :: img(*)
  double precision      :: img_frc(3, *)
  integer               :: ifracts(3, *)
  double precision      :: frcx, frcy, frcz

  integer               :: my_chgs(*), my_chg_cnt

! Local variables:

  integer               :: i, j, k
  integer               :: i_base, j_base, k_base
  integer               :: ith1, ith2, ith3
  integer               :: my_chg_idx
  double precision      :: f1, f2, f3
  double precision      :: charge, qterm
  double precision      :: f1_term, f2_term, f3_term
  double precision      :: dfx, dfy, dfz
  double precision      :: recip_11, recip_22, recip_33
  double precision      :: dnfft1, dnfft2, dnfft3

  integer               :: img_i
! integer               :: kbot0
  integer               :: kbot, ktop

  integer               :: i_tbl(-order : nfft1)
  integer               :: j_tbl(-order : nfft2)
  integer               :: k_tbl(-order : nfft3)

#ifdef RANGE_DBG
  integer               :: dbg_img_map(natom)
  integer               :: dbg_img_map_res(natom)

  dbg_img_map(:) = 0
#endif /* RANGE_DBG */


! Old usages; kept for documentation during code conversion...
! kbot0 = my_xy_slab_start
! kbot = kbot0 + 1
! ktop = kbot0 + my_xy_slab_cnt

  ! Initialize the indexing tables.  It actually produces faster code to
  ! do this here rather than caching the info remotely.

  do i = 0, nfft1
    i_tbl(i) = i + 1
  end do
  do i = -order, -1
    i_tbl(i) = i + nfft1 + 1
  end do

  do j = 0, nfft2
    j_tbl(j) = j + 1
  end do
  do j = -order, -1
    j_tbl(j) = j + nfft2 + 1
  end do

  if( numtasks .gt. 1) then
      kbot = 1
      ktop = my_xy_slab_cnt
      do k = 0, nfft3
        k_tbl(k) = k - my_xy_slab_start + 1
        if (k_tbl(k) .lt. kbot .or. k_tbl(k) .gt. ktop) k_tbl(k) = 0
      end do
      do k = -order, -1
        k_tbl(k) = k + nfft3 - my_xy_slab_start + 1
        if (k_tbl(k) .lt. kbot .or. k_tbl(k) .gt. ktop) k_tbl(k) = 0
      end do
  else
      do k = 0, nfft3
        k_tbl(k) = k + 1
      end do
      do k = -order, -1
        k_tbl(k) = k + nfft3 + 1
      end do
  endif

  recip_11 = recip(1, 1)
  recip_22 = recip(2, 2)
  recip_33 = recip(3, 3)

  dnfft1 = dble(nfft1)
  dnfft2 = dble(nfft2)
  dnfft3 = dble(nfft3)

  do my_chg_idx = 1, my_chg_cnt
    img_i = my_chgs(my_chg_idx)
    charge = img(img_i)%qterm

    i_base = ifracts(1, my_chg_idx)
    j_base = ifracts(2, my_chg_idx)
    k_base = ifracts(3, my_chg_idx)

    f1 = 0.d0
    f2 = 0.d0
    f3 = 0.d0

    ! We special-case order 4, the default.
    ! We special-case order 6, for future amoeba use.

    if (order .eq. 4) then
      do ith3 = 1, 4

        k = k_tbl(k_base + ith3)

        if (k .ne. 0) then
#ifdef RANGE_DBG
          if(numtasks .gt. 1) then
              dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
          endif
#endif /* RANGE_DBG */
          do ith2 = 1, 4

            j = j_tbl(j_base + ith2)

            f1_term = theta(ith2, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
            f2_term = dtheta(ith2, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
            f3_term = theta(ith2, 2, my_chg_idx) * dtheta(ith3, 3, my_chg_idx)

! Force is negative of grad:

            qterm = q(i_tbl(i_base+1), j, k)
            f1 = f1 - qterm * dtheta(1, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(1, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(1, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+2), j, k)
            f1 = f1 - qterm * dtheta(2, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(2, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(2, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+3), j, k)
            f1 = f1 - qterm * dtheta(3, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(3, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(3, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+4), j, k)
            f1 = f1 - qterm * dtheta(4, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(4, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(4, 1, my_chg_idx) * f3_term

          end do ! (ith2)
        end if ! ( k != 0)
      end do ! (ith3)
    else if (order .eq. 6) then
      do ith3 = 1, 6

        k = k_tbl(k_base + ith3)

        if (k .ne. 0) then
#ifdef RANGE_DBG
          if(numtasks .gt. 1) then
            dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
          endif
#endif /* RANGE_DBG */
          do ith2 = 1, 6

            j = j_tbl(j_base + ith2)

            f1_term = theta(ith2, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
            f2_term = dtheta(ith2, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
            f3_term = theta(ith2, 2, my_chg_idx) * dtheta(ith3, 3, my_chg_idx)

! Force is negative of grad:

            qterm = q(i_tbl(i_base+1), j, k)
            f1 = f1 - qterm * dtheta(1, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(1, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(1, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+2), j, k)
            f1 = f1 - qterm * dtheta(2, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(2, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(2, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+3), j, k)
            f1 = f1 - qterm * dtheta(3, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(3, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(3, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+4), j, k)
            f1 = f1 - qterm * dtheta(4, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(4, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(4, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+5), j, k)
            f1 = f1 - qterm * dtheta(5, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(5, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(5, 1, my_chg_idx) * f3_term
            qterm = q(i_tbl(i_base+6), j, k)
            f1 = f1 - qterm * dtheta(6, 1, my_chg_idx) * f1_term
            f2 = f2 - qterm * theta(6, 1, my_chg_idx) * f2_term
            f3 = f3 - qterm * theta(6, 1, my_chg_idx) * f3_term

          end do ! (ith2)
        end if ! ( k != 0)
      end do ! (ith3)
    else ! (order != 4 && != 6)
      do ith3 = 1, order

        k = k_tbl(k_base + ith3)

        if (k .ne. 0) then
#ifdef RANGE_DBG
          if( numtasks .gt. 1) then
            dbg_img_map(my_chgs(my_chg_idx)) = dbg_img_map(my_chgs(my_chg_idx))+1
          endif
#endif /* RANGE_DBG */
          do ith2 = 1, order

            j = j_tbl(j_base + ith2)

            f1_term = theta(ith2, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
            f2_term = dtheta(ith2, 2, my_chg_idx) * theta(ith3, 3, my_chg_idx)
            f3_term = theta(ith2, 2, my_chg_idx) * dtheta(ith3, 3, my_chg_idx)

            do ith1 = 1, order

              qterm = q(i_tbl(i_base+ith1), j, k)

! Force is negative of grad:

              f1 = f1 - qterm * dtheta(ith1, 1, my_chg_idx) * f1_term
              f2 = f2 - qterm * theta(ith1, 1, my_chg_idx) * f2_term
              f3 = f3 - qterm * theta(ith1, 1, my_chg_idx) * f3_term
  
            end do ! (ith1)
          end do ! (ith2)
        end if ! ( k != 0)
      end do ! (ith3)
    end if ! (order == 4,6 or other)

    f1 = f1 * dnfft1 * charge
    f2 = f2 * dnfft2 * charge
    f3 = f3 * dnfft3 * charge

    if (is_orthog .ne. 0) then
      dfx = recip_11 * f1
      dfy = recip_22 * f2
      dfz = recip_33 * f3
    else
      dfx = recip(1, 1) * f1 + recip(1, 2) * f2 + recip(1, 3) * f3
      dfy = recip(2, 1) * f1 + recip(2, 2) * f2 + recip(2, 3) * f3
      dfz = recip(3, 1) * f1 + recip(3, 2) * f2 + recip(3, 3) * f3
    end if

    img_frc(1, img_i) = img_frc(1, img_i) + dfx
    img_frc(2, img_i) = img_frc(2, img_i) + dfy
    img_frc(3, img_i) = img_frc(3, img_i) + dfz

    frcx = frcx + dfx
    frcy = frcy + dfy
    frcz = frcz + dfz

  end do ! (my_chg_idx) 

#ifdef RANGE_DBG
  if( numtasks .gt. 1) then
    call mpi_reduce(dbg_img_map, dbg_img_map_res, size(dbg_img_map), &
                      mpi_integer, mpi_sum, 0, lib_mpi_comm, err_code_mpi)
    if (master) then
        do i = 1, natom
            if (dbg_img_map_res(i) .ne. order) then
                write(0,*)'DBG: IMAGE ', i, 'processed ', &
                          dbg_img_map_res(i), 'times in grad_sum code!!!'
            end if
        end do
    end if
  endif
#endif /* RANGE_DBG */

  return

end subroutine grad_sum


!*******************************************************************************
!
! Subroutine:  scalar_sumrc
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine scalar_sumrc(q, ewaldcof, vol, prefac1, prefac2, prefac3, &
                        nfft1, nfft2, nfft3, x_dim, y_dim, z_dim, &
                        eer, vir)

!  use parallel_dat_mod
  use pbc_mod
  use pme_fft_mod

  implicit none

! Formal arguments:

  integer               :: nfft1, nfft2, nfft3
  integer               :: x_dim, y_dim, z_dim
  double precision      :: q(2,z_dim,x_dim,y_dim)
  double precision      :: ewaldcof
  double precision      :: vol
  double precision      :: prefac1(nfft1), prefac2(nfft2), prefac3(nfft3)
  double precision      :: eer
  double precision      :: vir(3,3)

! Local variables:

  double precision      :: energy, fac_2, pi_vol_inv, qterm
  double precision      :: eterm, eterm12, eterm_struc2, vterm
  double precision      :: eterms, eterms12, eterms_struc2s, vterms
  double precision      :: mhat1, mhat2, mhat3, msq_inv
  double precision      :: mhat1s, mhat2s, mhat3s, msqs_inv
  double precision      :: recip_11, recip_22, recip_33
  double precision      :: vir_11, vir_22, vir_33, vir_21, vir_31, vir_32
  integer               :: k1, k2, k2q, k3, k3_start, m1, m2, m3
  integer               :: k1s, k2s, k3s, m1s, m2s, m3s
  integer               :: nf1, nf2, nf3
  integer               :: m3_tbl(nfft3)
  integer               :: k3s_tbl(nfft3)
  integer               :: m3s_tbl(nfft3)

  recip_11 = recip(1, 1)
  recip_22 = recip(2, 2)
  recip_33 = recip(3, 3)

  fac_2 = (2.d0 * PI * PI) / (ewaldcof * ewaldcof)

  pi_vol_inv = 1.d0 / (PI * vol)

  nf1 = nfft1 / 2
  ! There is an assumption that nfft1 must be even!
  nf2 = nfft2 / 2
  if (2 * nf2 .lt. nfft2) nf2 = nf2 + 1
  nf3 = nfft3 / 2
  if (2 * nf3 .lt. nfft3) nf3 = nf3 + 1

  energy = 0.d0

  vir_11 = 0.d0
  vir_22 = 0.d0
  vir_33 = 0.d0

  vir_21 = 0.d0
  vir_31 = 0.d0
  vir_32 = 0.d0

! Tables used to avoid branching in the innermost loop:

  do k3 = 1, nfft3

    if (k3 .le. nf3) then
      m3_tbl(k3) = k3 - 1
    else
      m3_tbl(k3)= k3 - 1 - nfft3
    end if

    k3s = mod(nfft3 - k3 + 1, nfft3) + 1
    k3s_tbl(k3) = k3s

    if (k3s .le. nf3) then
      m3s_tbl(k3) = k3s - 1
    else
      m3s_tbl(k3) = k3s - 1 - nfft3
    end if

  end do

! Insist that q(1,1,1,1) is set to 0.d0 (true already for neutral)
! All results using these elements are calculated, but add 0.d0, so
! it is like they are not used.

  if(my_zx_slab_start .eq. 0) then
    q(1, 1, 1, 1) = 0.d0
    q(2, 1, 1, 1) = 0.d0
  endif

! Big loop:

! k2q is the z index into the actual q array in this process;
! k2 is the index that would be used if the entire q array, which only exists
! for the uniprocessor case.

  do k2q = 1, my_zx_slab_cnt

    k2 = k2q + my_zx_slab_start

    if (k2 .le. nf2) then
      m2 = k2 - 1
    else
      m2 = k2 - 1 - nfft2
    end if

    mhat2 = recip_22 * m2

    k2s = mod(nfft2 - k2 + 1, nfft2) + 1
    
    if (k2s .le. nf2) then
      m2s = k2s - 1
    else
      m2s = k2s - 1 - nfft2
    end if

    mhat2s = recip_22 * m2s

    do k1 = 1, nf1 + 1

      if (k1 .le. nf1) then
        m1 = k1 - 1
      else
        m1 = k1 - 1 - nfft1
      end if
      mhat1 = recip_11 * m1
      eterm12 = m1_exp_tbl(m1) * m2_exp_tbl(m2) * &
                prefac1(k1) * prefac2(k2) * pi_vol_inv

      k3_start = 1
      if (my_zx_slab_start .eq. 0) then
        if (k2 .eq. 1) then
          if (k1 .eq. 1) then
            k3_start = 2
          end if
        end if
      end if

      if (k1 .gt. 1 .and. k1 .le. nfft1) then

        k1s = nfft1 - k1 + 2
        if (k1s .le. nf1) then
          m1s = k1s - 1
        else
          m1s = k1s - 1 - nfft1
        end if
        mhat1s = recip_11 * m1s
        eterms12 = m1_exp_tbl(m1s) * m2_exp_tbl(m2s) * &
                   prefac1(k1s) * prefac2(k2s) * pi_vol_inv

        do k3 = k3_start, nfft3
          m3 = m3_tbl(k3)
          mhat3 = recip_33 * m3
          k3s = k3s_tbl(k3)
          m3s = m3s_tbl(k3)
          mhat3s = recip_33 * m3s
          msq_inv = 1.d0 / (mhat1 * mhat1 + mhat2 * mhat2 + mhat3 * mhat3)
          msqs_inv = 1.d0 / (mhat1s*mhat1s + mhat2s*mhat2s + mhat3s*mhat3s)
          ! The product of the following 3 table lookups is exp(-fac * msq):
          ! (two of the lookups occurred in calculating eterm12)
          eterm = eterm12 * m3_exp_tbl(m3) * prefac3(k3) * msq_inv
          ! The product of the following 3 table lookups is exp(-fac * msqs):
          ! (two of the lookups occurred in calculating eterms12)
          eterms = eterms12 * m3_exp_tbl(m3s) * prefac3(k3s) * msqs_inv
          qterm = (q(1, k3, k1, k2q) * q(1, k3, k1, k2q) + &
                   q(2, k3, k1, k2q) * q(2, k3, k1, k2q))
          q(1, k3, k1, k2q) = eterm * q(1, k3, k1, k2q)
          q(2, k3, k1, k2q) = eterm * q(2, k3, k1, k2q)
          eterm_struc2 = eterm * qterm
          eterms_struc2s = eterms * qterm
          energy = energy + eterm_struc2 + eterms_struc2s
          vterm = (fac_2 + 2.d0 * msq_inv) * eterm_struc2
          vterms = (fac_2 + 2.d0 * msqs_inv) * eterms_struc2s
          vir_21 = vir_21 + vterm * mhat1 * mhat2 + vterms * mhat1s * mhat2s
          vir_31 = vir_31 + vterm * mhat1 * mhat3 + vterms * mhat1s * mhat3s
          vir_32 = vir_32 + vterm * mhat2 * mhat3 + vterms * mhat2s * mhat3s
          vir_11 = vir_11 + vterm * mhat1 * mhat1 - eterm_struc2 + &
                            vterms * mhat1s * mhat1s - eterms_struc2s
          vir_22 = vir_22 + vterm * mhat2 * mhat2 - eterm_struc2 + &
                            vterms * mhat2s * mhat2s - eterms_struc2s
          vir_33 = vir_33 + vterm * mhat3 * mhat3 - eterm_struc2 + &
                            vterms * mhat3s * mhat3s - eterms_struc2s
        end do

      else

        do k3 = k3_start, nfft3
          m3 = m3_tbl(k3)
          mhat3 = recip_33 * m3
          msq_inv = 1.d0 / (mhat1 * mhat1 + mhat2 * mhat2 + mhat3 * mhat3)
          ! The product of the following 3 table lookups is exp(-fac * msq):
          ! (two of the lookups occurred in calculating eterm12)
          eterm = eterm12 * m3_exp_tbl(m3) * prefac3(k3) * msq_inv
          qterm = (q(1, k3, k1, k2q) * q(1, k3, k1, k2q) + &
                   q(2, k3, k1, k2q) * q(2, k3, k1, k2q))
          q(1, k3, k1, k2q) = eterm * q(1, k3, k1, k2q)
          q(2, k3, k1, k2q) = eterm * q(2, k3, k1, k2q)
          eterm_struc2 = eterm * qterm
          energy = energy + eterm_struc2
          vterm = (fac_2 + 2.d0 * msq_inv) * eterm_struc2
          vir_21 = vir_21 + vterm * mhat1 * mhat2
          vir_31 = vir_31 + vterm * mhat1 * mhat3
          vir_32 = vir_32 + vterm * mhat2 * mhat3
          vir_11 = vir_11 + vterm * mhat1 * mhat1 - eterm_struc2
          vir_22 = vir_22 + vterm * mhat2 * mhat2 - eterm_struc2
          vir_33 = vir_33 + vterm * mhat3 * mhat3 - eterm_struc2
        end do

      end if

    end do
  end do

  eer = 0.5d0 * energy

  vir(1, 1) = 0.5d0 * vir_11
  vir(2, 1) = 0.5d0 * vir_21
  vir(3, 1) = 0.5d0 * vir_31

  vir(1, 2) = 0.5d0 * vir_21
  vir(2, 2) = 0.5d0 * vir_22
  vir(3, 2) = 0.5d0 * vir_32

  vir(1, 3) = 0.5d0 * vir_31
  vir(2, 3) = 0.5d0 * vir_32
  vir(3, 3) = 0.5d0 * vir_33

  return

end subroutine scalar_sumrc

!*******************************************************************************
!
! Subroutine:  scalar_sumrc_nonorthog
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine scalar_sumrc_nonorthog(q, ewaldcof, vol, &
                                  prefac1, prefac2, prefac3, &
                                  nfft1, nfft2, nfft3, &
                                  x_dim, y_dim, z_dim, &
                                  eer, vir)

  use parallel_dat_mod
  use pbc_mod
  use pme_fft_mod

  implicit none

! Formal arguments:

  integer               :: nfft1, nfft2, nfft3
  integer               :: x_dim, y_dim, z_dim
  double precision      :: q(2,z_dim,x_dim,y_dim)
  double precision      :: ewaldcof
  double precision      :: vol
  double precision      :: prefac1(nfft1), prefac2(nfft2), prefac3(nfft3)
  double precision      :: eer
  double precision      :: vir(3,3)

! Local variables:

  double precision      :: energy, fac, fac_2, pi_vol_inv
  double precision      :: eterm, eterm_struc2, vterm
  double precision      :: eterms, eterms_struc2s, vterms
  double precision      :: mhat1, mhat2, mhat3, msq, msq_inv
  double precision      :: mhat1s, mhat2s, mhat3s, msqs, msqs_inv
  double precision      :: recip_stk(3,3)
  double precision      :: vir_11, vir_22, vir_33, vir_21, vir_31, vir_32
  integer               :: k, k1, k2, k2q, k3, k1_start, m1, m2, m3
  integer               :: k1s, k2s, k3s, m1s, m2s, m3s
  integer               :: nf1, nf2, nf3

  recip_stk(:,:) = recip(:,:)

  fac = (PI * PI) / (ewaldcof * ewaldcof)

  fac_2 = 2.d0 * fac

  pi_vol_inv = 1.d0 / (PI * vol)

  nf1 = nfft1 / 2
  if (2 * nf1 .lt. nfft1) nf1 = nf1 + 1
  nf2 = nfft2 / 2
  if (2 * nf2 .lt. nfft2) nf2 = nf2 + 1
  nf3 = nfft3 / 2
  if (2 * nf3 .lt. nfft3) nf3 = nf3 + 1

  energy = 0.d0

  vir_11 = 0.d0
  vir_22 = 0.d0
  vir_33 = 0.d0

  vir_21 = 0.d0
  vir_31 = 0.d0
  vir_32 = 0.d0

! Insist that q(1,1,1,1) is set to 0 (true already for neutral)

  if(my_zx_slab_start .eq. 0) then
    q(1, 1, 1, 1) = 0.d0
    q(2, 1, 1, 1) = 0.d0
  endif

! Big loop:

  do k2q = 1, my_zx_slab_cnt

    k2 = k2q + my_zx_slab_start

    if (k2 .le. nf2) then
      m2 = k2 - 1
    else
      m2 = k2 - 1 - nfft2
    end if

    k2s = mod(nfft2 - k2 + 1, nfft2) + 1
    
    if (k2s .le. nf2) then
      m2s = k2s - 1
    else
      m2s = k2s - 1 - nfft2
    end if

    do k3 = 1, nfft3

      k1_start = 1
      if (my_zx_slab_start .eq. 0) then
        if (k3 + k2 .eq. 2) k1_start = 2
      endif

      if (k3 .le. nf3) then
        m3 = k3 - 1
      else
        m3 = k3 - 1 - nfft3
      end if

      k3s = mod(nfft3 - k3 + 1, nfft3) + 1

      if (k3s .le. nf3) then
        m3s = k3s - 1
      else
        m3s = k3s - 1 - nfft3
      end if

      do k1 = k1_start, nf1 + 1

        if (k1 .le. nf1) then
          m1 = k1 - 1
        else
          m1 = k1 - 1 - nfft1
        end if

        mhat1 = recip_stk(1,1) * m1 + recip_stk(1,2) * m2 + recip_stk(1,3) * m3
        mhat2 = recip_stk(2,1) * m1 + recip_stk(2,2) * m2 + recip_stk(2,3) * m3
        mhat3 = recip_stk(3,1) * m1 + recip_stk(3,2) * m2 + recip_stk(3,3) * m3

        msq = mhat1 * mhat1 + mhat2 * mhat2 + mhat3 * mhat3

        msq_inv = 1.d0 / msq

        ! Getting the exponential via table lookup is currently not done
        ! for nonorthogonal unit cells.

        eterm = exp(-fac * msq) * prefac1(k1) * prefac2(k2) * prefac3(k3) * &
                pi_vol_inv * msq_inv

        vterm = fac_2 + 2.d0 * msq_inv

        eterm_struc2 = eterm * (q(1, k3, k1, k2q) * q(1, k3, k1, k2q) + &
                                q(2, k3, k1, k2q) * q(2, k3, k1, k2q))

        energy = energy + eterm_struc2

        vir_21 = vir_21 + eterm_struc2 * (vterm * mhat1 * mhat2)
        vir_31 = vir_31 + eterm_struc2 * (vterm * mhat1 * mhat3)
        vir_32 = vir_32 + eterm_struc2 * (vterm * mhat2 * mhat3)

        vir_11 = vir_11 + eterm_struc2 * (vterm * mhat1 * mhat1 - 1.d0)
        vir_22 = vir_22 + eterm_struc2 * (vterm * mhat2 * mhat2 - 1.d0)
        vir_33 = vir_33 + eterm_struc2 * (vterm * mhat3 * mhat3 - 1.d0)

        if (k1 .gt. 1 .and. k1 .le. nfft1) then

          k1s = nfft1 - k1 + 2

          if (k1s .le. nf1) then
            m1s = k1s - 1
          else
            m1s = k1s - 1 - nfft1
          end if

          mhat1s = recip_stk(1,1) * m1s + &
                   recip_stk(1,2) * m2s + &
                   recip_stk(1,3) * m3s

          mhat2s = recip_stk(2,1) * m1s + &
                   recip_stk(2,2) * m2s + &
                   recip_stk(2,3) * m3s

          mhat3s = recip_stk(3,1) * m1s + &
                   recip_stk(3,2) * m2s + &
                   recip_stk(3,3) * m3s

          msqs = mhat1s * mhat1s + mhat2s * mhat2s + mhat3s * mhat3s

          msqs_inv = 1.d0 / msqs

          ! Getting the exponential via table lookup is currently not done
          ! for nonorthogonal unit cells.

          eterms = exp(-fac * msqs) * &
                   prefac1(k1s) * prefac2(k2s) * prefac3(k3s) * &
                   pi_vol_inv * msqs_inv

          vterms = fac_2 + 2.d0 * msqs_inv

          eterms_struc2s = eterms * (q(1, k3, k1, k2q) * q(1, k3, k1, k2q) + &
                                     q(2, k3, k1, k2q) * q(2, k3, k1, k2q))

          energy = energy + eterms_struc2s

          vir_21 = vir_21 + eterms_struc2s * (vterms * mhat1s * mhat2s)
          vir_31 = vir_31 + eterms_struc2s * (vterms * mhat1s * mhat3s)
          vir_32 = vir_32 + eterms_struc2s * (vterms * mhat2s * mhat3s)

          vir_11 = vir_11 + eterms_struc2s * (vterms * mhat1s * mhat1s - 1.d0)
          vir_22 = vir_22 + eterms_struc2s * (vterms * mhat2s * mhat2s - 1.d0)
          vir_33 = vir_33 + eterms_struc2s * (vterms * mhat3s * mhat3s - 1.d0)

        endif

        q(1, k3, k1, k2q) = eterm * q(1, k3, k1, k2q)
        q(2, k3, k1, k2q) = eterm * q(2, k3, k1, k2q)

      end do
    end do
  end do

  eer = 0.5d0 * energy

  vir(1, 1) = 0.5d0 * vir_11
  vir(2, 1) = 0.5d0 * vir_21
  vir(3, 1) = 0.5d0 * vir_31

  vir(1, 2) = 0.5d0 * vir_21
  vir(2, 2) = 0.5d0 * vir_22
  vir(3, 2) = 0.5d0 * vir_32

  vir(1, 3) = 0.5d0 * vir_31
  vir(2, 3) = 0.5d0 * vir_32
  vir(3, 3) = 0.5d0 * vir_33

  return

end subroutine scalar_sumrc_nonorthog

!*******************************************************************************
!
! Subroutine:  get_grid_weights
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine get_grid_weights(img_cnt, img, ifracts, theta, dtheta, order, &
                            img_atm_map, my_chgs, my_chg_cnt)

  use bspline_mod
  use pme_fft_mod
  use img_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: img_cnt
  type(img_rec)         :: img(*)
  integer               :: ifracts(3, *)
  integer               :: order
  double precision      :: theta(order, 3, *)
  double precision      :: dtheta(order, 3, *)

  integer               :: img_atm_map(*)
  integer               :: my_chgs(*), my_chg_cnt


! Local variables:

  double precision      :: box_x, box_y, box_z
  double precision      :: crd_x, crd_y, crd_z
  double precision      :: factor1, factor2, factor3
  double precision      :: fract(3)
  double precision      :: weight(3)
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
      ktop1 = kbot0 + my_xy_slab_cnt + order - 2

      kwraps = (ktop1 .ge. nfft3)
  else  ! ( numtasks == 1)
      range_cnt = 1
      img_lo = 1
      img_hi = img_cnt
  end if 
  
  box_x = ucell(1, 1)
  box_y = ucell(2, 2)
  box_z = ucell(3, 3)

  ! Scaling factors to get from cit table indexes to nfft indexes:

  factor1 = dble(nfft1) * recip(1, 1)
  factor2 = dble(nfft2) * recip(2, 2)
  factor3 = dble(nfft3) * recip(3, 3)

  do range_ctr = 1, range_cnt

    do i = img_lo, img_hi

      if (img_atm_map(i) .lt. 0) cycle    ! Skip unused images.

      crd_x = img(i)%x
      crd_y = img(i)%y
      crd_z = img(i)%z

      if (crd_x .ge. box_x) then
        crd_x = crd_x - box_x
      else if (crd_x .lt. 0.d0) then
        crd_x = crd_x + box_x
      end if

      if (crd_y .ge. box_y) then
        crd_y = crd_y - box_y
      else if (crd_y .lt. 0.d0) then
        crd_y = crd_y + box_y
      end if

      if (crd_z .ge. box_z) then
        crd_z = crd_z - box_z
      else if (crd_z .lt. 0.d0) then
        crd_z = crd_z + box_z
      end if

      fract(1) = factor1 * crd_x
      fract(2) = factor2 * crd_y
      fract(3) = factor3 * crd_z

      if( numtasks .gt. 1) then
          k = int(fract(3))
          if (kwraps) then
            if (k .lt. kbot0 .and. k .gt. ktop1 - nfft3) cycle
          else
            if (k .lt. kbot0 .or. k .gt. ktop1) cycle
          end if
      endif

      my_chg_cnt = my_chg_cnt + 1
      my_chgs(my_chg_cnt) = i

      weight(:) = fract(:) - int(fract(:))

      ifracts(:, my_chg_cnt) = int(fract(:)) - order
      call fill_bspline_1_3d(weight, order, &
                             theta(1, 1, my_chg_cnt), dtheta(1, 1, my_chg_cnt))
    end do ! (i)

    if (recip_img_range_wraps) then
      img_lo = img_lo2
      img_hi = img_hi2
    end if

  end do ! (range_ctr)

  return

end subroutine get_grid_weights

!*******************************************************************************
!
! Subroutine:  claim_recip_imgs
!
! Description: <TBS>
!              Used only in MPI              
!*******************************************************************************

subroutine claim_recip_imgs(img_cnt, fraction, box, crd_idx_tbl, img, &
                            img_iac, img_atm_map)

  use cit_mod
  use pme_fft_mod
  use img_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: img_cnt
  double precision      :: fraction(3, img_cnt)
  double precision      :: box(3)
  type(cit_tbl_rec)     :: crd_idx_tbl(0 : cit_tbl_x_dim - 1, &
                                       0 : cit_tbl_y_dim - 1, &
                                       0 : cit_tbl_z_dim - 1)
  type(img_rec)         :: img(*)
  integer               :: img_iac(img_cnt)
  integer               :: img_atm_map(img_cnt)

! Local variables:

  integer               :: k_bot0, k_top1
  integer               :: i, j, k
  integer               :: atm_id
  integer               :: z_bkt_lo, z_bkt_hi
  logical               :: z_wraps
  double precision      :: x_box, y_box, z_box
  double precision      :: z_lo, z_hi, z_cur

  double precision      :: crd_x, crd_y, crd_z

  integer               :: img_lo, img_hi
  
  if (.not. i_do_recip) then    ! Map null range; not actually ever used...
    recip_img_lo = 1
    recip_img_hi = 0
    recip_img_range_wraps = .false.
    call get_max_recip_imgs(img_cnt)
    return
  end if

  x_box = box(1)
  y_box = box(2)
  z_box = box(3)

  ! If there really is only one task doing reciprocal force calcs, we have to
  ! catch it here, because the code further below assumes that one process does
  ! not need to map the entire range of images.  Besides, this is more
  ! efficient for the special case...

  if (my_xy_slab_cnt .eq. fft_z_dim) then

    recip_img_lo = 1
    recip_img_hi = img_cnt
    recip_img_range_wraps = .false.
    call get_max_recip_imgs(img_cnt)

    do i = recip_img_lo, recip_img_hi
      atm_id = img_atm_map(i)
      if (atm_id .lt. 0) then
        atm_id = - atm_id
        img_atm_map(i) = atm_id
        img(i)%x = fraction(1, atm_id) * x_box
        img(i)%y = fraction(2, atm_id) * y_box
        img(i)%z = fraction(3, atm_id) * z_box
        img(i)%qterm = atm_qterm(atm_id)
        img_iac(i) = atm_iac(atm_id)
      end if
    end do

    return

  end if

  k_bot0 = xy_slab_start(mytaskid)
  k_top1 = k_bot0 + my_xy_slab_cnt + bspl_order - 2
  
  ! We have to allow for motion of skinnb, plus we have to allow for fact that
  ! atom z crds are rounded down to the grid point.  Hence, the + 1 for z_hi

  z_lo = dble(k_bot0) * z_box / dble(nfft3) - (0.5d0 * skinnb)
  z_hi = dble(k_top1 + 1) * z_box / dble(nfft3) + (0.5d0 * skinnb)

  if (z_lo .lt. 0.d0) then
    z_wraps = .true.
    z_lo = z_lo + z_box
  else if (z_hi .ge. z_box) then
    z_wraps = .true.
    z_hi = z_hi - z_box
  else
    z_wraps = .false.
  end if

  z_bkt_lo = int(z_lo * dble(cit_tbl_z_dim) / z_box)
  z_bkt_hi = int(z_hi * dble(cit_tbl_z_dim) / z_box)

!  if (z_bkt_lo .eq. z_bkt_hi) z_wraps = .false.  ! TEMPORAL IGOR FIX 
  
  if (z_wraps) then

    img_hi = 0

    find_img_hi_loop1: &
    do k = z_bkt_hi, 0, -1
      do j = cit_tbl_y_dim - 1, 0, -1
        do i = cit_tbl_x_dim - 1, 0, -1

          if (crd_idx_tbl(i,j,k)%img_lo .ne. 0) then
            img_hi = crd_idx_tbl(i,j,k)%img_hi
            z_bkt_hi = k
            exit find_img_hi_loop1
          end if

        end do
      end do
    end do find_img_hi_loop1

    if (img_hi .eq. 0) then

      z_bkt_hi = cit_tbl_z_dim - 1
      z_wraps = .false.

    else

      find_img_lo_loop1: &
      do k = z_bkt_lo, cit_tbl_z_dim - 1
        do j = 0, cit_tbl_y_dim - 1
          do i = 0, cit_tbl_x_dim - 1
            img_lo = crd_idx_tbl(i,j,k)%img_lo
            if (crd_idx_tbl(i,j,k)%img_lo .ne. 0) then
              exit find_img_lo_loop1
            end if
          end do
        end do
      end do find_img_lo_loop1

      if (img_lo .eq. 0) then

        z_bkt_lo = 0
        z_wraps = .false.
    
      end if

    end if

  end if
  
!  write(iunit_debug,*) " claim_recip_imgs() pt 3   z_wraps =",z_wraps

  if (.not. z_wraps) then

    img_hi = 0

    find_img_hi_loop2: &
    do k = z_bkt_hi, z_bkt_lo, -1
      do j = cit_tbl_y_dim - 1, 0, -1
        do i = cit_tbl_x_dim - 1, 0, -1

          if (crd_idx_tbl(i,j,k)%img_lo .ne. 0) then
            img_hi = crd_idx_tbl(i,j,k)%img_hi
            z_bkt_hi = k
            exit find_img_hi_loop2
          end if

        end do
      end do
    end do find_img_hi_loop2

    if (img_hi .eq. 0) then   ! Empty range.  Highly unlikely.
      recip_img_lo = 0
      recip_img_hi = -1
      recip_img_range_wraps = .false.
      write(0,*)'DBG: Empty range return!!'
      call get_max_recip_imgs(img_cnt)
      return
    end if

    find_img_lo_loop2: &
    do k = z_bkt_lo, z_bkt_hi
      do j = 0, cit_tbl_y_dim - 1
        do i = 0, cit_tbl_x_dim - 1
          img_lo = crd_idx_tbl(i,j,k)%img_lo
          if (crd_idx_tbl(i,j,k)%img_lo .ne. 0) then
            exit find_img_lo_loop2
          end if
        end do
      end do
    end do find_img_lo_loop2

    ! Not possible for the above img_lo find to fail...

  end if

  ! NOTE - The range recip_img_lo to recip_img_hi is not all image-mapped; thus
  !        all users of this range must check for mapping; the range purely
  !        establishes outer bounds over which you have to look for mapped
  !        images!!!

  recip_img_lo = img_lo
  recip_img_hi = img_hi
  recip_img_range_wraps = z_wraps
  
  ! BUGBUG:
  ! The following is not as efficient as possible, but relatively foolproof.
  ! It can be set up to skip over the owned image range later...

  if (recip_img_range_wraps) then

    do i = recip_img_lo, img_cnt

      atm_id = img_atm_map(i)

      if (atm_id .lt. 0) then

        atm_id = - atm_id

        z_cur = fraction(3, atm_id) * z_box

        if (z_cur .ge. z_lo) then
          img_atm_map(i) = atm_id
          img(i)%x = fraction(1, atm_id) * x_box
          img(i)%y = fraction(2, atm_id) * y_box
          img(i)%z = z_cur
          img(i)%qterm = atm_qterm(atm_id)
          img_iac(i) = atm_iac(atm_id)
        end if

      end if

    end do

    do i = 1, recip_img_hi

      atm_id = img_atm_map(i)

      if (atm_id .lt. 0) then

        atm_id = - atm_id

        z_cur = fraction(3, atm_id) * z_box

        if (z_cur .lt. z_hi) then
          img_atm_map(i) = atm_id
          img(i)%x = fraction(1, atm_id) * x_box
          img(i)%y = fraction(2, atm_id) * y_box
          img(i)%z = z_cur
          img(i)%qterm = atm_qterm(atm_id)
          img_iac(i) = atm_iac(atm_id)
        end if

      end if

    end do

  else  ! Range does not wrap...

    do i = recip_img_lo, recip_img_hi

      atm_id = img_atm_map(i)

      if (atm_id .lt. 0) then

        atm_id = - atm_id

        z_cur = fraction(3, atm_id) * z_box

        if (z_cur .ge. z_lo .and. z_cur .lt. z_hi) then
          img_atm_map(i) = atm_id
          img(i)%x = fraction(1, atm_id) * x_box
          img(i)%y = fraction(2, atm_id) * y_box
          img(i)%z = z_cur
          img(i)%qterm = atm_qterm(atm_id)
          img_iac(i) = atm_iac(atm_id)
        end if

      end if

    end do

  end if

  call get_max_recip_imgs(img_cnt)

! BEGIN DBG
! write(0,*)'DBG: task, recip_img_lo,hi=', mytaskid, recip_img_lo, recip_img_hi
! write(0,*)'DBG: task, recip_img_range_wraps=', mytaskid, recip_img_range_wraps
! END DBG

  return

end subroutine claim_recip_imgs

!*******************************************************************************
!
! Subroutine:  claim_recip_imgs_nonorthog
!
! Description: <TBS>
!              Used only in MPI              
!*******************************************************************************

subroutine claim_recip_imgs_nonorthog(img_cnt, fraction, crd_idx_tbl, img, &
                                      img_iac, img_atm_map)

  use cit_mod
  use pme_fft_mod
  use img_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: img_cnt
  double precision      :: fraction(3, img_cnt)
  type(cit_tbl_rec)     :: crd_idx_tbl(0 : cit_tbl_x_dim - 1, &
                                       0 : cit_tbl_y_dim - 1, &
                                       0 : cit_tbl_z_dim - 1)
  type(img_rec)         :: img(*)
  integer               :: img_iac(img_cnt)
  integer               :: img_atm_map(img_cnt)

! Local variables:

  integer               :: k_bot0, k_top1
  integer               :: i, j, k
  integer               :: atm_id
  integer               :: z_bkt_lo, z_bkt_hi
  logical               :: z_wraps
  double precision      :: ucell_stk(3, 3)
  double precision      :: f1, f2, f3
! double precision      :: z_lo, z_hi, z_cur
  double precision      :: z_lo, z_hi

  double precision      :: crd_x, crd_y, crd_z

  integer               :: img_lo, img_hi

  ucell_stk(:,:) = ucell(:,:)

  if (.not. i_do_recip) then    ! Map null range; not actually ever used...
    recip_img_lo = 1
    recip_img_hi = 0
    recip_img_range_wraps = .false.
    call get_max_recip_imgs(img_cnt)
    return
  end if

  ! If there really is only one task doing reciprocal force calcs, we have to
  ! catch it here, because the code further below assumes that one process does
  ! not need to map the entire range of images.  Besides, this is more
  ! efficient for the special case...

  if (my_xy_slab_cnt .eq. fft_z_dim) then

    recip_img_lo = 1
    recip_img_hi = img_cnt
    recip_img_range_wraps = .false.
    call get_max_recip_imgs(img_cnt)

    do i = recip_img_lo, recip_img_hi
      atm_id = img_atm_map(i)
      if (atm_id .lt. 0) then
        atm_id = - atm_id
        img_atm_map(i) = atm_id
        f1 = fraction(1, atm_id)
        f2 = fraction(2, atm_id)
        f3 = fraction(3, atm_id)
        ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
        ! we can simplify the expression in this critical inner loop
        img(i)%x = f1 * ucell_stk(1, 1) + &
                   f2 * ucell_stk(1, 2) + &
                   f3 * ucell_stk(1, 3)
        img(i)%y = f2 * ucell_stk(2, 2) + &
                   f3 * ucell_stk(2, 3)
        img(i)%z = f3 * ucell_stk(3, 3)
        img(i)%qterm = atm_qterm(atm_id)
        img_iac(i) = atm_iac(atm_id)
      end if
    end do

    return

  end if

  k_bot0 = xy_slab_start(mytaskid)
  k_top1 = k_bot0 + my_xy_slab_cnt + bspl_order - 2

  ! We have to allow for motion of skinnb, plus we have to allow for fact that
  ! atom z crds are rounded down to the grid point.  Hence, the + 1 for z_hi

  z_lo = dble(k_bot0) / dble(nfft3) - &
         (0.5d0 * cut_factor(3) * skinnb) / ucell_stk(3,3)      ! a fractional

  z_hi = dble(k_top1 + 1) / dble(nfft3) + &
         (0.5d0 * cut_factor(3) * skinnb) / ucell_stk(3,3)      ! a fractional

  if (z_lo .lt. 0.d0) then
    z_wraps = .true.
    z_lo = z_lo + 1.d0
  else if (z_hi .ge. 1.d0) then
    z_wraps = .true.
    z_hi = z_hi - 1.d0
  else
    z_wraps = .false.
  end if

  z_bkt_lo = int(z_lo * dble(cit_tbl_z_dim))
  z_bkt_hi = int(z_hi * dble(cit_tbl_z_dim))

  if (z_bkt_lo .eq. z_bkt_hi) z_wraps = .false.

  if (z_wraps) then

    img_hi = 0

    find_img_hi_loop1: &
    do k = z_bkt_hi, 0, -1
      do j = cit_tbl_y_dim - 1, 0, -1
        do i = cit_tbl_x_dim - 1, 0, -1

          if (crd_idx_tbl(i,j,k)%img_lo .ne. 0) then
            img_hi = crd_idx_tbl(i,j,k)%img_hi
            z_bkt_hi = k
            exit find_img_hi_loop1
          end if

        end do
      end do
    end do find_img_hi_loop1

    if (img_hi .eq. 0) then

      z_bkt_hi = cit_tbl_z_dim - 1
      z_wraps = .false.

    else

      find_img_lo_loop1: &
      do k = z_bkt_lo, cit_tbl_z_dim - 1
        do j = 0, cit_tbl_y_dim - 1
          do i = 0, cit_tbl_x_dim - 1
            img_lo = crd_idx_tbl(i,j,k)%img_lo
            if (crd_idx_tbl(i,j,k)%img_lo .ne. 0) then
              exit find_img_lo_loop1
            end if
          end do
        end do
      end do find_img_lo_loop1

      if (img_lo .eq. 0) then

        z_bkt_lo = 0
        z_wraps = .false.
    
      end if

    end if

  end if

  if (.not. z_wraps) then

    img_hi = 0

    find_img_hi_loop2: &
    do k = z_bkt_hi, z_bkt_lo, -1
      do j = cit_tbl_y_dim - 1, 0, -1
        do i = cit_tbl_x_dim - 1, 0, -1

          if (crd_idx_tbl(i,j,k)%img_lo .ne. 0) then
            img_hi = crd_idx_tbl(i,j,k)%img_hi
            z_bkt_hi = k
            exit find_img_hi_loop2
          end if

        end do
      end do
    end do find_img_hi_loop2

    if (img_hi .eq. 0) then   ! Empty range.  Highly unlikely.
      recip_img_lo = 0
      recip_img_hi = -1
      recip_img_range_wraps = .false.
      write(0,*)'DBG: Empty range return!!'
      call get_max_recip_imgs(img_cnt)
      return
    end if

    find_img_lo_loop2: &
    do k = z_bkt_lo, z_bkt_hi
      do j = 0, cit_tbl_y_dim - 1
        do i = 0, cit_tbl_x_dim - 1
          img_lo = crd_idx_tbl(i,j,k)%img_lo
          if (crd_idx_tbl(i,j,k)%img_lo .ne. 0) then
            exit find_img_lo_loop2
          end if
        end do
      end do
    end do find_img_lo_loop2

    ! Not possible for the above img_lo find to fail...

  end if

  ! NOTE - The range recip_img_lo to recip_img_hi is not all image-mapped; thus
  !        all users of this range must check for mapping; the range purely
  !        establishes outer bounds over which you have to look for mapped
  !        images!!!

  recip_img_lo = img_lo
  recip_img_hi = img_hi
  recip_img_range_wraps = z_wraps

  ! BUGBUG:
  ! The following is not as efficient as possible, but relatively foolproof.
  ! It can be set up to skip over the owned image range later...

  if (recip_img_range_wraps) then

    do i = recip_img_lo, img_cnt

      atm_id = img_atm_map(i)

      if (atm_id .lt. 0) then

        atm_id = - atm_id

        f3 = fraction(3, atm_id)
!       z_cur = f3 * ucell_stk(3, 3)
        
!       if (z_cur .ge. z_lo) then
        if (f3 .ge. z_lo) then
          img_atm_map(i) = atm_id

          f1 = fraction(1, atm_id)
          f2 = fraction(2, atm_id)

          ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
          ! we can simplify the expression in this critical inner loop:

          img(i)%x = f1 * ucell_stk(1, 1) + &
                     f2 * ucell_stk(1, 2) + &
                     f3 * ucell_stk(1, 3)

          img(i)%y = f2 * ucell_stk(2, 2) + &
                     f3 * ucell_stk(2, 3)

          img(i)%z = f3 * ucell_stk(3, 3)

          img(i)%qterm = atm_qterm(atm_id)
          img_iac(i) = atm_iac(atm_id)
        end if

      end if

    end do

    do i = 1, recip_img_hi

      atm_id = img_atm_map(i)

      if (atm_id .lt. 0) then

        atm_id = - atm_id

        f3 = fraction(3, atm_id)
!       z_cur = f3 * ucell_stk(3, 3)
        
!       if (z_cur .lt. z_hi) then
        if (f3 .lt. z_hi) then
          img_atm_map(i) = atm_id

          f1 = fraction(1, atm_id)
          f2 = fraction(2, atm_id)

          ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
          ! we can simplify the expression in this critical inner loop

          img(i)%x = f1 * ucell_stk(1, 1) + &
                     f2 * ucell_stk(1, 2) + &
                     f3 * ucell_stk(1, 3)

          img(i)%y = f2 * ucell_stk(2, 2) + &
                     f3 * ucell_stk(2, 3)

          img(i)%z = f3 * ucell_stk(3, 3)

          img(i)%qterm = atm_qterm(atm_id)
          img_iac(i) = atm_iac(atm_id)
        end if

      end if

    end do

  else  ! Range does not wrap...

    do i = recip_img_lo, recip_img_hi

      atm_id = img_atm_map(i)

      if (atm_id .lt. 0) then

        atm_id = - atm_id

        f3 = fraction(3, atm_id)
!       z_cur = f3 * ucell_stk(3, 3)
        
!       if (z_cur .ge. z_lo .and. z_cur .lt. z_hi) then
        if (f3 .ge. z_lo .and. f3 .lt. z_hi) then
          img_atm_map(i) = atm_id

          f1 = fraction(1, atm_id)
          f2 = fraction(2, atm_id)

          ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
          ! we can simplify the expression in this critical inner loop

          img(i)%x = f1 * ucell_stk(1, 1) + &
                     f2 * ucell_stk(1, 2) + &
                     f3 * ucell_stk(1, 3)

          img(i)%y = f2 * ucell_stk(2, 2) + &
                     f3 * ucell_stk(2, 3)

          img(i)%z = f3 * ucell_stk(3, 3)

          img(i)%qterm = atm_qterm(atm_id)
          img_iac(i) = atm_iac(atm_id)
        end if

      end if

    end do

  end if

  call get_max_recip_imgs(img_cnt)

! BEGIN DBG
! write(0,*)'DBG: task, recip_img_lo,hi=', mytaskid, recip_img_lo, recip_img_hi
! write(0,*)'DBG: task, recip_img_range_wraps=', mytaskid, recip_img_range_wraps
! END DBG

  return

end subroutine claim_recip_imgs_nonorthog

!*******************************************************************************
!
! Subroutine:   get_max_recip_imgs
!
! Description: <TBS>
!              Used only in MPI
!*******************************************************************************

subroutine get_max_recip_imgs(img_cnt)

  implicit none

! Formal arguments:

  integer               :: img_cnt

  if (recip_img_range_wraps) then
    max_recip_imgs = img_cnt - recip_img_lo + 1
    max_recip_imgs = max_recip_imgs + recip_img_hi
  else
    max_recip_imgs = recip_img_hi - recip_img_lo + 1
  end if

  return

end subroutine get_max_recip_imgs

!*******************************************************************************
!
! Subroutine:  get_grid_weights_nonorthog
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine get_grid_weights_nonorthog(img_cnt, crd, img_atm_map, ifracts, &
                                      theta, dtheta, order, my_chgs, &
                                      my_chg_cnt)

  use bspline_mod
  use pme_fft_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: img_cnt
  double precision      :: crd(3, *)
  integer               :: img_atm_map(*)
  integer               :: ifracts(3, *)
  integer               :: order
  double precision      :: theta(order, 3, *)
  double precision      :: dtheta(order, 3, *)

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
      ktop1 = kbot0 + my_xy_slab_cnt + order - 2

      kwraps = (ktop1 .ge. nfft3)
  else ! (numtasks == 1)
      range_cnt = 1
      img_lo = 1
      img_hi = img_cnt
  endif

  recip_stk(:,:) = recip(:,:)


  do range_ctr = 1, range_cnt
    do i = img_lo, img_hi

      if (img_atm_map(i) .lt. 0) cycle    ! Skip unused images.

! Unfortunately we need fractional coords that are not available without going
! back to atom coordinates crd().  This is the case because the coordinates in img_crd
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

          if (kwraps) then
              if (k .lt. kbot0 .and. k .gt. ktop1 - nfft3) cycle
          else
              if (k .lt. kbot0 .or. k .gt. ktop1) cycle
          end if
      endif 

      my_chg_cnt = my_chg_cnt + 1
      my_chgs(my_chg_cnt) = i

      weight(:) = fract(:) - int(fract(:))

      ifracts(:, my_chg_cnt) = int(fract(:)) - order

      call fill_bspline_1_3d(weight, order, &
                             theta(1, 1, my_chg_cnt), dtheta(1, 1, my_chg_cnt))

    end do ! (i) 
    if (recip_img_range_wraps) then
      img_lo = img_lo2
      img_hi = img_hi2
    end if
  end do ! (range_ctr)

  return

end subroutine get_grid_weights_nonorthog

!*******************************************************************************
!
! Subroutine:  allocate_m_exp_tbls
!
! Description: <TBS>
!
!*******************************************************************************

subroutine allocate_m_exp_tbls

  use mdin_ewald_dat_mod

  implicit none

! Local variables:

  integer               :: alloc_failed

  if( allocated(m1_exp_tbl)) deallocate(m1_exp_tbl)
  if( allocated(m2_exp_tbl)) deallocate(m2_exp_tbl)
  if( allocated(m3_exp_tbl)) deallocate(m3_exp_tbl)

  allocate(m1_exp_tbl(-(nfft1/2 + 1) : nfft1/2 + 1), &
           m2_exp_tbl(-(nfft2/2 + 1) : nfft2/2 + 1), &
           m3_exp_tbl(-(nfft3/2 + 1) : nfft3/2 + 1), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  return

end subroutine allocate_m_exp_tbls

!*******************************************************************************
!
! Subroutine:  load_m_exp_tbls
!
! Description: <TBS>
!
!*******************************************************************************

subroutine load_m_exp_tbls

  use mdin_ewald_dat_mod
  use pbc_mod

  implicit none

! Local variables:

  double precision      :: fac
  double precision      :: recip_11_2, recip_22_2, recip_33_2
  integer               :: i

  fac = -(PI * PI)/(ew_coeff * ew_coeff)

  recip_11_2 = recip(1,1) * recip(1,1)
  recip_22_2 = recip(2,2) * recip(2,2)
  recip_33_2 = recip(3,3) * recip(3,3)

  do i = -(nfft1/2 + 1), nfft1/2 + 1
    m1_exp_tbl(i) = exp(fac * dble(i) * dble(i) * recip_11_2)
  end do

  do i = -(nfft2/2 + 1), nfft2/2 + 1
    m2_exp_tbl(i) = exp(fac * dble(i) * dble(i) * recip_22_2)
  end do

  do i = -(nfft3/2 + 1), nfft3/2 + 1
    m3_exp_tbl(i) = exp(fac * dble(i) * dble(i) * recip_33_2)
  end do

  return

end subroutine load_m_exp_tbls

!*******************************************************************************
!
! Subroutine:  pme_recip_setup
!
! Description:
!              
!*******************************************************************************

subroutine pme_recip_setup(imin_par)

  use pme_fft_mod
  use loadbal_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer imin_par

! Local variables:

  integer               :: alloc_failed
  
  if( allocated(gbl_prefac1)) deallocate(gbl_prefac1)
  if( allocated(gbl_prefac2)) deallocate(gbl_prefac2)
  if( allocated(gbl_prefac3)) deallocate(gbl_prefac3)

  allocate(gbl_prefac1(nfft1), &
           gbl_prefac2(nfft2), &
           gbl_prefac3(nfft3), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  gbl_prefac1(:) = 0.d0
  gbl_prefac2(:) = 0.d0
  gbl_prefac3(:) = 0.d0

  call load_prefacs(gbl_prefac1, gbl_prefac2, gbl_prefac3)

  if( numtasks .gt. 1 ) then
  
  ! For mpi we do fft loadbalancing unless nrespa is in effect, this is
  ! a minimization, or we are using Amoeba pme.  We also only do loadbalancing
  ! if the number of slabs per task goes below a threshhold, or if we
  ! have 2 procs and have defined the manifest const FFTLOADBAL_2PROC.
  
      if (nrespa .eq. 1 .and. &
          imin_par .eq. 0 .and. iamoeba .eq. 0 .and. &
#ifdef FFTLOADBAL_2PROC
          (numtasks * bspl_order .gt. nfft3 .or. numtasks .eq. 2)) then
#else
           numtasks * bspl_order .gt. nfft3) then
#endif
        fft_slab_redist_enabled = .true.
        fft_slab_redist_needed = .true.
      end if   
  endif
      
  call fft_setup(recip_numtasks, fft_slab_redist_enabled)
  
  if( numtasks .eq. 1 ) then
     max_recip_imgs = natom
  endif
  
  if (is_orthog .ne. 0) then
    call allocate_m_exp_tbls
    call load_m_exp_tbls
  end if

  return

contains

!*******************************************************************************
!
! Internal Subroutine:  load_prefacs
!
! Description:  This routine loads the moduli of the inverse DFT of the
!               B splines.  bsp_mod1-3 hold these values, and nfft1-3 are the
!               grid dimensions.  Order is the order of the B spline approx.
!
!*******************************************************************************

subroutine load_prefacs(prefac1, prefac2, prefac3)

  use bspline_mod
  use mdin_ewald_dat_mod

  implicit none

  ! Formal arguments:

  double precision      :: prefac1(*), prefac2(*), prefac3(*)

  ! Local variables:

  double precision      :: array(bspl_order)
  double precision      :: bsp_arr(max(nfft1,nfft2,nfft3))
  double precision      :: bsp_mod(max(nfft1,nfft2,nfft3))
  double precision      :: w
  integer               :: i

  w = 0.d0
  call fill_bspline_0(w, bspl_order, array)

  bsp_arr(:) = 0.d0

  do i = 1, bspl_order - 1
   bsp_arr(i) = array(bspl_order - i)
  end do

  call dftmod(bsp_mod, bsp_arr, nfft1)

  call factor_lambda(bsp_mod, nfft1, prefac1)

  call dftmod(bsp_mod, bsp_arr, nfft2)

  call factor_lambda(bsp_mod, nfft2, prefac2)

  call dftmod(bsp_mod, bsp_arr, nfft3)

  call factor_lambda(bsp_mod, nfft3, prefac3)

  return

end subroutine load_prefacs

!*******************************************************************************
!
! Internal Subroutine:  dftmod
!
! Description:  Computes the modulus of the discrete fourier transform of
!               bsp_arr, storing it into bsp_mod.
!
!*******************************************************************************

subroutine dftmod(bsp_mod, bsp_arr, nfft)

  implicit none

  integer               nfft

  double precision      bsp_mod(nfft), bsp_arr(nfft)

  integer               j, k

  double precision      sum1, sum2, twopi, arg, tiny

  twopi = 2.d0 * PI

  tiny = 1.d-7

  do k = 1, nfft
    sum1 = 0.d0
    sum2 = 0.d0

    do j = 1, nfft
      arg = twopi * (k - 1) * (j - 1)/nfft
      sum1 = sum1 + bsp_arr(j) * cos(arg)  ! Re( M_p(k+1) *exp(2*pi*i*m_i/K_i)
      sum2 = sum2 + bsp_arr(j) * sin(arg)
    end do

    bsp_mod(k) = sum1**2 + sum2**2 

  end do

! Fix the ONE case where exponential Euler spline interpolation fails.
! This arbitrary assignment to avoid division by zero is okay
! since it happens with highest frequency reciprocal vectors out in tail
! of reciprocal sum.

  do k = 1, nfft
    if (bsp_mod(k) .lt. tiny) bsp_mod(k) = 0.5d0 * (bsp_mod(k-1) + bsp_mod(k+1))
  end do

  return

end subroutine dftmod

!*******************************************************************************
!
! Internal Subroutine:  factor_lambda
!
! Description:  Factor in optimal lambda coefficient for bspline approximant
!               of complex exponentials, thus modify influence function.
!
!*******************************************************************************

subroutine factor_lambda(bsp_mod, nfft, prefac)
! construct full prefac(k) = |b_i(m_i/K_i)|^2  by (lambda_p(m_i/K_i))^2/bsp_mod(m_i/K_i)    
! Multiply lambda_p(m_i/K_i) - optimizing coef to fit complex exponential in the least-square sense formulas (2.31 - 2.32) of Darden_2000
  use mdin_ewald_dat_mod

  implicit none

  integer               nfft

  double precision      bsp_mod(nfft), prefac(nfft)

  integer               k, nf, m, order2

  double precision      lambda, gsum, gsum2

  integer, parameter    :: kcut = 50

  nf = nfft / 2
  order2 = 2 * bspl_order

  do k = 1, nfft
    m = k - 1
    if (k .gt. nf)m = k - 1 - nfft
    order2 = 2 * bspl_order
    if (m .eq. 0) then
      lambda = 1.d0
    else
      call gamma_sum(gsum, m, nfft, bspl_order, kcut)
      call gamma_sum(gsum2, m, nfft, order2, kcut)
      lambda = gsum/gsum2
    end if
    prefac(k) = lambda**2/bsp_mod(k)
  end do

  return

end subroutine factor_lambda

!*******************************************************************************
!
! Internal Subroutine:  gamma_sum
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine gamma_sum(gsum, m, nfft, order, kcut)
!
! gsum = m_i/K_i * sum_k ( 1/( m_i/K_i + k )^p  k = -kcut, + kcut   formula (2.32)  jcp_00_113_10913 (Darden 2000)
  implicit none

  double precision      gsum

  integer               m, nfft, order, kcut

  double precision      frac, x

  integer               k

  if (m .eq. 0) then
    gsum = 1.d0
    return
  end if

  frac = dble(m)/nfft   !  frac = m_i/K_i      
  x = PI * frac         !  x = pi * m_i/K_i 
  gsum = 1.d0

  do k = 1, kcut
    gsum = gsum + (x/(x + PI * k))**order        ! m_i/K_i * sum_k ( 1/( m_i/K_i + k )^p  k = - Inf, + Inf
  end do 

  do k = 1, kcut
    gsum = gsum + (x/(x - PI * k))**order
  end do 

  return

end subroutine gamma_sum

end subroutine pme_recip_setup

end module pme_recip_mod
