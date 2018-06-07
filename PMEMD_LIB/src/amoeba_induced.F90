#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_induced_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_induced_mod

  implicit none

  private

  integer, save ::   do_amoeba_induced_flag, polar_atm_cnt

  double precision, save, allocatable   :: polarizability(:)  
  double precision, save, allocatable   :: screen_polar(:)
  double precision, save, allocatable   :: damp_polar_strength(:)  ! Damping polarization coefficent 1, if > 0 atom is affecting polarization of other atoms
  double precision, save, allocatable   :: damp_polar_sensitivity(:)  ! Damping polarization coefficent 2, if > 0 atom polarizabilty is affected by the presence of other atoms
  double precision, save, allocatable   :: damp_polar_rad(:)    ! atom radius to determine range of polarization damping interactions
  double precision, save, allocatable   :: polarizability_corr(:) ! corrections to polarizabilities due to presence of other atoms  
  double precision, save, allocatable   :: hpolar(:)   ! atom hyperpolarizabilities
   
  ! Data that is allocated and initialized in final amoeba setup.  Some
  ! of this may not need to be global (BUGBUG - check it out).

  logical, save, allocatable            :: is_polarizable(:)
  double precision, save, allocatable   :: ind_dip_d(:,:)        !  induced dipoles due to ?
  double precision, save, allocatable   :: ind_dip_p(:,:)        !  induced dipoles due to ?
  

  ! Only relevant in master...

  logical, save, public                 :: print_amoeba_dip_info

  public        polarizability
  public        polarizability_corr
  public        damp_polar_strength
  public        damp_polar_sensitivity
  public        damp_polar_rad
  public        screen_polar
  public        is_polarizable
  public        ind_dip_d
  public        ind_dip_p
  public        am_induced_zero_flag
  public        am_induced_set_user_bit
  public        am_induced_eval_mpi
  public        am_induced_eval_uni
  public        am_var_polar_frc
  
  integer, parameter, public :: ndamp = 6

contains

!*******************************************************************************!
! Subroutine:  am_induced_zero_flag
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_induced_zero_flag
  implicit none
  do_amoeba_induced_flag = 0
  return
end subroutine am_induced_zero_flag

!*******************************************************************************!
! Subroutine:  am_induced_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_induced_set_user_bit(do_this)

  use amoeba_flags_mod
  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in) :: do_this

  if (do_this .eq. 1) then
    do_amoeba_induced_flag = ibset(do_amoeba_induced_flag, user_bit)
    print_amoeba_dip_info = .true.
  else
    do_amoeba_induced_flag = ibclr(do_amoeba_induced_flag, user_bit)
    print_amoeba_dip_info = .false.
  end if

  return

end subroutine am_induced_set_user_bit

!*******************************************************************************!
! Subroutine:  am_induced_eval
!
! Description: <TBS>
!
!*******************************************************************************

! MPI IMPLEMENTATION
subroutine  am_induced_eval_mpi(atm_cnt, crd, diprms, dipiter, num_adjust_list, &
                               tranvec )

  use mdin_amoeba_dat_mod, only : sor_coefficient, dipole_scf_tol, &
                                  dipole_scf_iter_max, amoeba_verbose

  use amoeba_flags_mod
  use file_io_dat_mod
  use parallel_dat_mod
  use timers_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  double precision, intent(in)  :: crd(3, *)
  double precision, intent(out) :: diprms
  double precision, intent(out) :: dipiter
  integer, intent(in)           :: num_adjust_list
  double precision, intent(in)  :: tranvec(*)

! Local variables:

  integer                       :: iter
  integer, save                 :: done
  double precision              :: rms1, rms2, rms, oldrms
  double precision              :: dip_field_d(3, atm_cnt)
  double precision              :: dip_field_p(3, atm_cnt)
  double precision              :: dip_d_perm(3, atm_cnt)
  double precision              :: dip_p_perm(3, atm_cnt)
  double precision              :: old_dip_d(3, atm_cnt)
  double precision              :: old_dip_p(3, atm_cnt)
  double precision              :: reduced_dip_d(3, atm_cnt)
  double precision              :: reduced_dip_p(3, atm_cnt)
  double precision              :: adj_dip_dip_tensor(6, num_adjust_list)

  diprms = 0.d0
  dipiter = 0.d0

  if (do_amoeba_induced_flag .ne. proceed) return

  dip_field_d(:,:) = 0.d0
  dip_field_p(:,:) = 0.d0
  polarizability_corr(:) = 0.d0
  
  call am_nonbond_perm_fields(atm_cnt, is_polarizable, crd, &
                              dip_field_d, dip_field_p, adj_dip_dip_tensor, &
                              tranvec)

  call update_time(nonbond_time)

  call mpi_allreduce(dip_field_d, dbl_mpi_recv_buf, 3 * atm_cnt, &
                     mpi_double_precision, mpi_sum, lib_mpi_comm, &
                     err_code_mpi)

  call array_copy(reduced_dip_d, dbl_mpi_recv_buf, 3 * atm_cnt)

  call mpi_allreduce(dip_field_p, dbl_mpi_recv_buf, 3 * atm_cnt, &
                     mpi_double_precision, mpi_sum, lib_mpi_comm, &
                     err_code_mpi)

  call array_copy(reduced_dip_p, dbl_mpi_recv_buf, 3 * atm_cnt)
  
  call mpi_allreduce( polarizability_corr, dbl_mpi_recv_buf, atm_cnt, &
                     mpi_double_precision, mpi_sum, lib_mpi_comm, &
                     err_code_mpi)

  call array_copy(polarizability_corr, dbl_mpi_recv_buf, atm_cnt)


  call update_time(fcve_dist_time)
  call zero_pme_time()

  call am_induced_fields_to_ind_dips(atm_cnt, is_polarizable, &
                                     dip_field_d, dip_field_p, &
                                     polarizability, hpolar, polarizability_corr, ind_dip_d, ind_dip_p)

  ! Do an all reduce on the ind_dip_[dp] arrays to get input for
  ! am_nonbond_dip_dip_fields

  call update_time(nonbond_time)

  call mpi_allreduce(ind_dip_d, dbl_mpi_recv_buf, 3 * atm_cnt, &
                     mpi_double_precision, mpi_sum, lib_mpi_comm, &
                     err_code_mpi)

  call array_copy(reduced_dip_d, dbl_mpi_recv_buf, 3 * atm_cnt)

  call mpi_allreduce(ind_dip_p, dbl_mpi_recv_buf, 3 * atm_cnt, &
                     mpi_double_precision, mpi_sum, lib_mpi_comm, &
                     err_code_mpi)

  call array_copy(reduced_dip_p, dbl_mpi_recv_buf, 3 * atm_cnt)

  ! Save the dips due to permanent fields

  dip_d_perm(:,1:atm_cnt) = reduced_dip_d(:,1:atm_cnt)
  dip_p_perm(:,1:atm_cnt) = reduced_dip_p(:,1:atm_cnt)

  call update_time(fcve_dist_time)
  call zero_pme_time()

  iter = 0
  rms = 1.d0

  done = 0
  
  do 

    ! copy dipoles to old

    old_dip_d(:,1:atm_cnt) = reduced_dip_d(:,1:atm_cnt)
    old_dip_p(:,1:atm_cnt) = reduced_dip_p(:,1:atm_cnt)

    call am_nonbond_dip_dip_fields(atm_cnt, crd, &
                                   reduced_dip_d, reduced_dip_p,  &
                                   adj_dip_dip_tensor, dip_field_d, dip_field_p)

    ! get dipoles due to dipole fields

    call am_induced_fields_to_ind_dips(atm_cnt, is_polarizable, &
                                       dip_field_d, dip_field_p, &
                                       polarizability, hpolar, polarizability_corr, ind_dip_d, ind_dip_p)

    call update_time(nonbond_time)

    ! Get the next reduced set of values for ind_dip_[dp]:

    call mpi_allreduce(ind_dip_d, dbl_mpi_recv_buf, 3 * atm_cnt, &
                       mpi_double_precision, mpi_sum, lib_mpi_comm, &
                       err_code_mpi)

    call array_copy(reduced_dip_d, dbl_mpi_recv_buf, 3 * atm_cnt)


    call mpi_allreduce(ind_dip_p, dbl_mpi_recv_buf, 3 * atm_cnt, &
                       mpi_double_precision, mpi_sum, lib_mpi_comm, &
                       err_code_mpi)

    call array_copy(reduced_dip_p, dbl_mpi_recv_buf, 3 * atm_cnt)

    call update_time(fcve_dist_time)
    call zero_pme_time()

    ! add dipoles due to permfields

    call array_add(reduced_dip_d, dip_d_perm, 3 * atm_cnt)
    call array_add(reduced_dip_p, dip_p_perm, 3 * atm_cnt)

    call am_induced_sor_update(atm_cnt, is_polarizable, sor_coefficient, &
                               reduced_dip_d, reduced_dip_p,  &
                               old_dip_d, old_dip_p)

    ! In theory, everyone would get the same result if they executed the
    ! following code.  In reality, there may be rounding errors in mpi
    ! dependent on operation order, so we just have the master do the next
    ! step.  Note however, that we have already gotten the result to all
    ! processes once we converge...

    if (master) then

      oldrms = rms

      call am_induced_rms_diff(atm_cnt, is_polarizable, rms1, rms2, &
                               reduced_dip_d, reduced_dip_p,  &
                               old_dip_d, old_dip_p)
      rms = max(rms1, rms2)
      iter = iter + 1

      call update_pme_time(pme_misc_timer)
      call update_time(nonbond_time)
      if (rms .lt. dipole_scf_tol) then
        done = 1
        call mpi_bcast(done, 1, mpi_integer, 0, lib_mpi_comm, err_code_mpi)
      else if (rms .gt. oldrms) then
        error_msg = 'dipoles failed to converge: rms increasing!'
        done = -1 ! Error flag...
        call mpi_bcast(done, 1, mpi_integer, 0, lib_mpi_comm, err_code_mpi)
        call mol_mech_error
      else if (iter .gt. dipole_scf_iter_max) then
        done = -1 ! Error flag...
        error_msg = 'dipoles failed to converge: iter max exceeded!'
        call mpi_bcast(done, 1, mpi_integer, 0, lib_mpi_comm, err_code_mpi)
        call mol_mech_error
      else
        call mpi_bcast(done, 1, mpi_integer, 0, lib_mpi_comm, err_code_mpi)
      end if
    else
      call update_pme_time(pme_misc_timer)
      call update_time(nonbond_time)
      call mpi_bcast(done, 1, mpi_integer, 0, lib_mpi_comm, err_code_mpi)
    end if
    call update_time(fcve_dist_time)
    call zero_pme_time()


    if (done .gt. 0) then
      exit
    else if (done .lt. 0) then
      call mol_mech_error
    end if

  end do

  ind_dip_d(:,:) = reduced_dip_d(:,:)
  ind_dip_p(:,:) = reduced_dip_p(:,:)

  if (master .and. amoeba_verbose .ne. 0) &
     write(mdout, '(a, i5, 1x, e12.4)') 'out of loop, iter,rms = ', iter, rms

  if (master) then
    diprms = rms
    dipiter = iter
  end if

  call update_pme_time(pme_misc_timer)

  return

end subroutine am_induced_eval_mpi
! UNIPROCESSOR IMPLEMENTATION
subroutine am_induced_eval_uni(atm_cnt, crd, diprms, dipiter, num_adjust_list, &
                               tranvec )

  use mdin_amoeba_dat_mod, only : sor_coefficient, dipole_scf_tol, &
                                  dipole_scf_iter_max, amoeba_verbose

  use amoeba_flags_mod
  use file_io_dat_mod
  use parallel_dat_mod
  use timers_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  double precision, intent(in)  :: crd(3, *)
  double precision, intent(out) :: diprms
  double precision, intent(out) :: dipiter
  integer, intent(in)           :: num_adjust_list
  double precision, intent(in)  :: tranvec(*)

! Local variables:

  integer                       :: iter
  integer                       :: i
  double precision              :: rms1, rms2, oldrms
  double precision, save        :: rms
  double precision              :: dip_field_d(3, atm_cnt)
  double precision              :: dip_field_p(3, atm_cnt)
  double precision              :: dip_d_perm(3, atm_cnt)
  double precision              :: dip_p_perm(3, atm_cnt)
  double precision              :: old_dip_d(3, atm_cnt)
  double precision              :: old_dip_p(3, atm_cnt)
  double precision              :: adj_dip_dip_tensor(6, num_adjust_list)

  if (do_amoeba_induced_flag .ne. proceed) return

  dip_field_d(:,:) = 0.d0
  dip_field_p(:,:) = 0.d0
  polarizability_corr(:) = 0.d0

  call am_nonbond_perm_fields(atm_cnt, is_polarizable, crd, &
                              dip_field_d, dip_field_p, adj_dip_dip_tensor, &
                              tranvec)
  
  call am_induced_fields_to_ind_dips(atm_cnt, is_polarizable, &
                                     dip_field_d, dip_field_p, &
                                     polarizability, hpolar, polarizability_corr, ind_dip_d, ind_dip_p)

  ! Save the dips due to permanent fields

  dip_d_perm(:,1:atm_cnt) = ind_dip_d(:,1:atm_cnt)
  dip_p_perm(:,1:atm_cnt) = ind_dip_p(:,1:atm_cnt)

  iter = 0
  rms = 1.d0
 
  call update_pme_time(pme_misc_timer)

  do 

    ! copy dipoles to old

    old_dip_d(:,1:atm_cnt) = ind_dip_d(:,1:atm_cnt)
    old_dip_p(:,1:atm_cnt) = ind_dip_p(:,1:atm_cnt)

    call am_nonbond_dip_dip_fields(atm_cnt, crd, ind_dip_d, ind_dip_p,  &
                                   adj_dip_dip_tensor, dip_field_d, dip_field_p) ! dip_field_d, dip_field_p - dipole fields due to induced dipoles

    ! get dipoles due to dipole fields

    call am_induced_fields_to_ind_dips(atm_cnt, is_polarizable, &
                                       dip_field_d, dip_field_p, &
                                       polarizability, hpolar, polarizability_corr, ind_dip_d, ind_dip_p)

    ! add dipoles due to permfields

    call array_add(ind_dip_d, dip_d_perm, 3 * atm_cnt)
    call array_add(ind_dip_p, dip_p_perm, 3 * atm_cnt)

    call am_induced_sor_update(atm_cnt, is_polarizable, sor_coefficient, &
                               ind_dip_d, ind_dip_p,  old_dip_d, old_dip_p)

    call am_induced_rms_diff(atm_cnt, is_polarizable, rms1, rms2, &
                             ind_dip_d, ind_dip_p,  old_dip_d, old_dip_p)

    oldrms = rms
    rms = max(rms1, rms2)
    iter = iter + 1

    if (rms .lt. dipole_scf_tol) exit

    if (rms .gt. oldrms) then
        error_msg = 'dipoles failed to converge: rms increasing!'
        call mol_mech_error
    end if

    if (iter .gt. dipole_scf_iter_max) then
        error_msg =  'dipoles failed to converge: iter max exceeded!'
        call mol_mech_error
    end if

  end do

  call update_pme_time(pme_misc_timer)

  if (master .and. amoeba_verbose .ne. 0) &
     write(mdout, '(a, i5, 1x, e12.4)') 'out of loop, iter,rms = ', iter, rms

  dipiter = iter
  diprms = rms

  return

end subroutine am_induced_eval_uni
                               
subroutine am_var_polar_frc(ipairs, img, tranvec, crd, e_hpol_ind,  &
                            frc, img_frc, virial, img_atm_map)

  use prmtop_dat_mod
  use mdin_ewald_dat_mod
  use mdin_amoeba_dat_mod, only : ee_dsum_cut
  use amoeba_multipoles_mod, only : coulomb_const_kcal_per_mole
  use img_mod
  use timers_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer                               :: ipairs(*)
  type(img_rec), intent(in)             :: img(*)
  double precision, intent(in)          :: tranvec(1:3, 0:17)
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: e_hpol_ind
  double precision, intent(in out)      :: frc(3, *)
  double precision, intent(in out)      :: img_frc(3, *)
  double precision, intent(in out)      :: virial(3, 3)
  integer, intent(in)                   :: img_atm_map(*)

! Local variables:

  double precision                      :: x_i, y_i, z_i
  double precision                      :: x_tran(1:3, 0:17)
  double precision                      :: ee_dsum_cut2
  double precision                      :: hpol, pol, hpol_pref 
  integer                               :: i
  integer                               :: itran
  integer                               :: ipairs_idx
  integer                               :: sublst_idx
  integer                               :: img_i,img_j
  integer                               :: pair_cnt
  integer                               :: atm_i,atm_j
  double precision                      :: delx, dely, delz
  double precision                      :: delr, delr2
  double precision                      :: pref, ind_dip_sq,ind_dip_3, cf
  double precision                      :: del_damp,del_damp_n1
  
  double precision, parameter   :: small = 1.d-24
  
#ifdef DIRFRC_COMTRANS
  ! flag - 1 if translation not needed
  integer                               :: common_tran
#endif /* DIRFRC_COMTRANS */
  integer, parameter    :: mask27 = Z"07FFFFFF"

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
      
      if( is_polarizable(atm_i) ) then
         pol =  polarizability(atm_i)
         hpol = hpolar(atm_i)
         if( pol .gt. small .and. hpol .gt. small ) then
             
             ind_dip_sq = ind_dip_p(1,atm_i)**2 + ind_dip_p(2,atm_i)**2 + ind_dip_p(3,atm_i)**2
             ind_dip_3 = ind_dip_sq * sqrt( ind_dip_sq )
             e_hpol_ind = e_hpol_ind - coulomb_const_kcal_per_mole * hpol * ind_dip_3 /( 3.0*pol*pol*pol) 
             
         endif
      endif  
    
      do sublst_idx = 1, pair_cnt

#ifdef DIRFRC_COMTRANS
        if (common_tran .eq. 1) then
            img_j = ipairs(ipairs_idx + sublst_idx)
            itran = 13
        else
#endif /* DIRFRC_COMTRANS */
            img_j = iand(ipairs(ipairs_idx + sublst_idx), mask27)
            itran = ishft(ipairs(ipairs_idx + sublst_idx), -27)
#ifdef DIRFRC_COMTRANS
        end if
#endif /* DIRFRC_COMTRANS */
        atm_j = img_atm_map(img_j)
        
         delx = img(img_j)%x + x_tran(1, itran)
         dely = img(img_j)%y + x_tran(2, itran)
         delz = img(img_j)%z + x_tran(3, itran)

         delr2 = delx * delx + dely * dely + delz * delz
         
         if (delr2 .lt. ee_dsum_cut2) then
           
             if( damp_polar_strength(atm_j) > 0.0 .and. damp_polar_sensitivity(atm_i) > 0.0 ) then
                 
                 delr = sqrt(delr2)
                 
                 del_damp = delr - (damp_polar_rad(atm_i) + damp_polar_rad(atm_j)) 
!                 del_damp7 = del_damp**7   ! IGOR CHANGING_DAMPING POWER
                 del_damp_n1 = del_damp ** (ndamp + 1)
                 cf = damp_polar_strength(atm_j) * damp_polar_sensitivity(atm_i) 
                  
!
! force corrections due to position dependence of polarizability
!

                 ind_dip_sq = ind_dip_p(1,atm_i)**2 + ind_dip_p(2,atm_i)**2 + ind_dip_p(3,atm_i)**2
                 
                 pol  = polarizability(atm_i)
                 hpol = hpolar(atm_i) 
                 
!                 pref = 6.0d0 * cf/del_damp7 /polarizability(atm_i) ! IGOR CHANGING_DAMPING POWER
!                 pref = 2.0d0 * cf/del_damp3 /polarizability(atm_i)
                  pref = ndamp * cf/del_damp_n1 /pol 
                  
                  hpol_pref = 0.0
                  if( pol .gt. small .and. hpol .gt. small ) then 
                      hpol_pref = hpol * sqrt(ind_dip_sq)/(pol * pol) 
                  endif
                  
                 pref = -pref * ind_dip_sq * ( 0.5d0 - hpol_pref )
                 
                 pref = pref * coulomb_const_kcal_per_mole
                 
                 img_frc(1, img_j) = img_frc(1, img_j) - pref * delx/delr
                 img_frc(2, img_j) = img_frc(2, img_j) - pref * dely/delr
                 img_frc(3, img_j) = img_frc(3, img_j) - pref * delz/delr
                 
                 img_frc(1, img_i) = img_frc(1, img_i) + pref * delx/delr
                 img_frc(2, img_i) = img_frc(2, img_i) + pref * dely/delr
                 img_frc(3, img_i) = img_frc(3, img_i) + pref * delz/delr
                
             endif
             
             if( damp_polar_strength(atm_i) .gt. 0.0 .and. damp_polar_sensitivity(atm_j) .gt. 0.0 ) then
                 
                 delr = sqrt(delr2)
                 
                 del_damp = delr - (damp_polar_rad(atm_i) + damp_polar_rad(atm_j)) 
!                 del_damp7 = del_damp**7  ! IGOR CHANGING_DAMPING POWER
                 del_damp_n1 = del_damp**(ndamp + 1)
                 cf = damp_polar_strength(atm_j) * damp_polar_sensitivity(atm_i) 
                 
!
! force corrections due to position dependence of polarizability
!

                 ind_dip_sq = ind_dip_p(1,atm_j)**2 + ind_dip_p(2,atm_j)**2 + ind_dip_p(3,atm_j)**2 
                 
                 pol  = polarizability(atm_j)
                 hpol = hpolar(atm_j) 
                 
!                 pref = 6.0d0 * cf/del_damp7 /polarizability(atm_j) ! IGOR CHANGING_DAMPING POWER
!                 pref = 2.0d0 * cf/del_damp3 /polarizability(atm_j)
                  pref = ndamp * cf/del_damp_n1 /pol 
                  
                  hpol_pref = 0.0
                  if( pol .gt. small .and. hpol .gt. small ) then 
                      hpol_pref = hpol*sqrt(ind_dip_sq)/(pol*pol) 
                  endif
                  
                 pref = -pref * ind_dip_sq * ( 0.5d0 - hpol_pref )
                 
                 pref = pref * coulomb_const_kcal_per_mole
                 
                 img_frc(1, img_i) = img_frc(1, img_i) - pref * delx/delr
                 img_frc(2, img_i) = img_frc(2, img_i) - pref * dely/delr
                 img_frc(3, img_i) = img_frc(3, img_i) - pref * delz/delr
                 
                 img_frc(1, img_j) = img_frc(1, img_j) + pref * delx/delr
                 img_frc(2, img_j) = img_frc(2, img_j) + pref * dely/delr
                 img_frc(3, img_j) = img_frc(3, img_j) + pref * delz/delr

              endif 
         endif
      end do

      ipairs_idx = ipairs_idx + pair_cnt
    end if
  end do  
  
  call update_pme_time(pme_misc_timer)
  
  return      
end subroutine am_var_polar_frc                           

subroutine set_atm_polar(polarizability_new,is_polarizable_new, screen_polar_new,  &
                         damp_polar_coef1_new, damp_polar_coef2_new, damp_polar_rad_new, hpolar_new,n)
		
	implicit none
	double precision, intent(in)         :: polarizability_new(n)
	integer, intent(in)                  :: is_polarizable_new(n)
    double precision, intent(in)         :: screen_polar_new(n)
	double precision, intent(in)         :: damp_polar_coef1_new(n)
    double precision, intent(in)         :: damp_polar_coef2_new(n)
    double precision, intent(in)         :: damp_polar_rad_new(n)
    double precision, intent(in)         :: hpolar_new(n)
	integer, intent(in)         :: n
	
	integer          ::   nold
	integer          ::   alloc_failed
	integer          ::   i
	
	nold = size(polarizability)
	polar_atm_cnt = n
	
	if( nold .ne. n ) then
	   if( nold .gt. 0) then
	      deallocate(polarizability)
          deallocate(screen_polar)
          deallocate(damp_polar_strength)
          deallocate(damp_polar_sensitivity)
          deallocate(damp_polar_rad)
          deallocate(hpolar)
          deallocate(is_polarizable)
	      deallocate(ind_dip_d)
	      deallocate(ind_dip_p)    
          deallocate(polarizability_corr)
	   endif
	   allocate( polarizability(polar_atm_cnt), &
                 screen_polar(polar_atm_cnt), &
                 damp_polar_strength(polar_atm_cnt), &
                 damp_polar_sensitivity(polar_atm_cnt), &
                 damp_polar_rad(polar_atm_cnt), &
                 hpolar(polar_atm_cnt), &
                 is_polarizable(polar_atm_cnt), &
                 ind_dip_d(3, polar_atm_cnt), &
                 ind_dip_p(3, polar_atm_cnt), &
                 polarizability_corr(polar_atm_cnt), &             
	             stat = alloc_failed )
	   if( alloc_failed .ne. 0 ) call setup_alloc_error
	endif
	
	if( polar_atm_cnt .eq. 0) return
	    
	polarizability(1:polar_atm_cnt)   = polarizability_new(1:polar_atm_cnt)
    screen_polar(1:polar_atm_cnt)     = screen_polar_new(1:polar_atm_cnt)
	damp_polar_strength(1:polar_atm_cnt) = damp_polar_coef1_new(1:polar_atm_cnt)
    damp_polar_sensitivity(1:polar_atm_cnt) = damp_polar_coef2_new(1:polar_atm_cnt)
    damp_polar_rad(1:polar_atm_cnt)   = damp_polar_rad_new(1:polar_atm_cnt)
    hpolar(1:polar_atm_cnt)   = hpolar_new(1:polar_atm_cnt)
	do i = 1, polar_atm_cnt
	    if( is_polarizable_new(i) .eq. 1) then 
	        is_polarizable(i) = .true. 
	    else
	        is_polarizable(i) = .false.
	    endif
	    ind_dip_d(:, i) = 0.d0
        ind_dip_p(:, i) = 0.d0
    enddo 
	
end subroutine set_atm_polar

subroutine set_valid_bit(ival)

    use amoeba_flags_mod
    
    implicit none
	integer, intent(in)         :: ival
    
    if(ival .eq. 1) then
        do_amoeba_induced_flag = ibset(do_amoeba_induced_flag, valid_bit)
    else
        do_amoeba_induced_flag = ibclr(do_amoeba_induced_flag, valid_bit)
    endif
        
end subroutine set_valid_bit

subroutine get_ind_dip( natom, ind_dip_d_out, ind_dip_p_out )

implicit none

    integer,          intent(in)   :: natom
    double precision, intent(out)  :: ind_dip_d_out(*)
    double precision, intent(out)  :: ind_dip_p_out(*)

    integer          ::   i

    if ( 3*natom > size(ind_dip_d) ) return
    if ( 3*natom > size(ind_dip_p) ) return

	do i = 1, natom
	    if( is_polarizable(i)) then 
	        ind_dip_d_out(3*i - 2) = ind_dip_d(1, i)
	        ind_dip_d_out(3*i - 1) = ind_dip_d(2, i)
	        ind_dip_d_out(3*i)     = ind_dip_d(3, i)
	        ind_dip_p_out(3*i - 2) = ind_dip_p(1, i)
	        ind_dip_p_out(3*i - 1) = ind_dip_p(2, i)
	        ind_dip_p_out(3*i)     = ind_dip_p(3, i)
	    else
	        ind_dip_d_out(3*i - 2) = 0.d0
	        ind_dip_d_out(3*i - 1) = 0.d0
	        ind_dip_d_out(3*i)     = 0.d0
	        ind_dip_p_out(3*i - 2) = 0.d0
	        ind_dip_p_out(3*i - 1) = 0.d0
	        ind_dip_p_out(3*i)     = 0.d0	        
	    endif
    enddo 

end subroutine get_ind_dip

end module amoeba_induced_mod

!*******************************************************************************!
! Subroutine:  am_induced_add_cart_to_dip
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_induced_add_cart_to_dip(atm_cnt, is_polarizable, &
                                      cart_dipole_field, &
                                      dip_field_d, dip_field_p)

  use timers_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: atm_cnt
  logical, intent(in)                   :: is_polarizable(*)
  double precision, intent(in)          :: cart_dipole_field(3, *)
  double precision, intent(in out)      :: dip_field_d(3, *)
  double precision, intent(in out)      :: dip_field_p(3, *)

! Local variables:

  integer                               :: n

  ! NOTE - This is intentionally over entire atom range at the moment.

  do n = 1, atm_cnt
    if (is_polarizable(n)) then
      dip_field_d(:, n) = dip_field_d(:, n) + cart_dipole_field(:, n)
      dip_field_p(:, n) = dip_field_p(:, n) + cart_dipole_field(:, n)
    end if
  end do

  call update_pme_time(pme_misc_timer)

  return

end subroutine am_induced_add_cart_to_dip

!*******************************************************************************!
! Subroutine:  am_induced_fields_to_ind_dips
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_induced_fields_to_ind_dips(atm_cnt, is_polarizable, &
                                         dip_field_d, dip_field_p, &
                                         polarizability, hpolar, polarizability_corr, &
                                         ind_dip_d, ind_dip_p)

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  logical, intent(in)           :: is_polarizable(*)
  double precision, intent(in)  :: dip_field_d(3, *)
  double precision, intent(in)  :: dip_field_p(3, *)
  double precision, intent(in)  :: polarizability(*)
  double precision, intent(in)  :: hpolar(*)
  double precision, intent(in)  :: polarizability_corr(*)
  double precision, intent(out) :: ind_dip_d(3, *)
  double precision, intent(out) :: ind_dip_p(3, *)

! Local variables:

  integer                       :: n
  double precision              :: pol
  double precision              :: hpol
  double precision              :: rf,hcf_d, hcf_p, df
  
  double precision, parameter   :: small = 1.d-24

  ! NOTE - This is intentionally over entire atom range at the moment.

  !write(*,*) "current polarizabilities:"
  do n = 1, atm_cnt

    if (is_polarizable(n)) then
      pol =  polarizability(n)
      hpol = hpolar(n)
      hcf_d = 1.0d0 ! hyperpolarizabilty correction for direct induced polarizability
      hcf_p = 1.0d0 ! hyperpolarizabilty correction for mutual induced polarizability
#ifdef VAR_POLAR
      if( polarizability_corr(n) .gt. 0.0 ) then 
          rf = 1.0d0/( 1.0d0 +  polarizability_corr(n) )
          pol = pol * rf
          hpol = hpol * rf * rf
      end if
      
      if( pol .gt. small .and. hpol .gt. small ) then
          df = sqrt( dip_field_d(1, n)*dip_field_d(1, n) + dip_field_d(2, n)*dip_field_d(2, n) + dip_field_d(3, n)*dip_field_d(3, n))
 !         hcf_d = 1 + df * hpol/pol
          if( 4.0*hpol*df < pol ) then
             hcf_d = 2.0/( 1.0 + sqrt(1.0 - 4.0*hpol*df/pol ))
          else
             hcf_d = 2.0
          endif
          ! write(*,'("hpol correction: = ",I3,3F12.6)') n, hpol, df, hcf_d
          df = sqrt( dip_field_p(1, n)*dip_field_p(1, n) + dip_field_p(2, n)*dip_field_p(2, n) + dip_field_p(3, n)*dip_field_p(3, n))
 !         hcf_p = 1 + df * hpol/pol
          if( 4.0*hpol*df < pol ) then 
             hcf_p = 2.0/( 1.0 + sqrt(1.0 - 4.0*hpol*df/pol ))
          else
             hcf_p = 2.0
          endif
      end if
      
#endif
      !write(*,*) n, " corr = ", polarizability_corr(n), " polar = ", pol
      !write(*,'("dip_field_d = ",I3,3(F12.6,2X))') n, dip_field_d(1, n), dip_field_d(2, n), dip_field_d(3, n)
      ! minus sign since our field is actually gradient of potential

      ind_dip_d(:, n) = -pol * dip_field_d(:, n) * hcf_d
      ind_dip_p(:, n) = -pol * dip_field_p(:, n) * hcf_p
      
      !write(*,'("hpolar = ",I3,3(F12.6,2X))') n, hpolar(n)
      !write(*,'("ind_dip_d = ",I3,3(F12.6,2X))') n, ind_dip_d(1, n), ind_dip_d(2, n), ind_dip_d(3, n)
    else
      ind_dip_d(:, n) = 0.d0
      ind_dip_p(:, n) = 0.d0
    end if
  end do

  return

end subroutine am_induced_fields_to_ind_dips

!*******************************************************************************!
! Subroutine:  am_induced_sor_update
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_induced_sor_update(atm_cnt, is_polarizable, sor_coeff, &
                                 ind_dip_d, ind_dip_p, old_dip_d, old_dip_p)

  implicit none

! Formal arguments:

  integer, intent(in)                   :: atm_cnt
  logical, intent(in)                   :: is_polarizable(*)
  double precision, intent(in)          :: sor_coeff
  double precision, intent(in out)      :: ind_dip_d(3, *)
  double precision, intent(in out)      :: ind_dip_p(3, *)
  double precision, intent(in)          :: old_dip_d(3, *)
  double precision, intent(in)          :: old_dip_p(3, *)

! Local variables:

  double precision                      :: c_sor
  integer                               :: i, n

  ! NOTE - This is intentionally over entire atom range at the moment.

  c_sor = 1.d0 - sor_coeff

  do n = 1, atm_cnt
    if (is_polarizable(n)) then
      ind_dip_d(:, n) = c_sor * old_dip_d(:, n) + sor_coeff * ind_dip_d(:, n)
      ind_dip_p(:, n) = c_sor * old_dip_p(:, n) + sor_coeff * ind_dip_p(:, n)
    end if
  end do

  return

end subroutine am_induced_sor_update

!*******************************************************************************!
! Subroutine:  am_induced_rms_diff
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_induced_rms_diff(atm_cnt, is_polarizable, rms1, rms2, &
                               ind_dip_d, ind_dip_p, old_dip_d, old_dip_p)

  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  logical, intent(in)           :: is_polarizable(*)
  double precision, intent(out) :: rms1
  double precision, intent(out) :: rms2
  double precision, intent(in)  :: ind_dip_d(3, *)
  double precision, intent(in)  :: ind_dip_p(3, *)
  double precision, intent(in)  :: old_dip_d(3, *)
  double precision, intent(in)  :: old_dip_p(3, *)

! Local variables:

  double precision              :: debye
  integer                       :: n, num

  debye = 4.8033324d0
  num = 0
  rms1 = 0.d0
  rms2 = 0.d0

  do n = 1, atm_cnt

    if (is_polarizable(n)) then
      num = num + 1

      rms1 = rms1 + (ind_dip_d(1, n)-old_dip_d(1, n))**2 + &
                    (ind_dip_d(2, n)-old_dip_d(2, n))**2 + &
                    (ind_dip_d(3, n)-old_dip_d(3, n))**2

      rms2 = rms2 + (ind_dip_p(1, n)-old_dip_p(1, n))**2 + &
                    (ind_dip_p(2, n)-old_dip_p(2, n))**2 + &
                    (ind_dip_p(3, n)-old_dip_p(3, n))**2
    end if
  end do

  if (num .eq. 0) then
    error_msg =  'am_induced_rms_diff: num polarizable = 0!!'
    call mol_mech_error
  end if

  rms1 = debye * sqrt(rms1 / num)
  rms2 = debye * sqrt(rms2 / num)

  return

end subroutine am_induced_rms_diff


