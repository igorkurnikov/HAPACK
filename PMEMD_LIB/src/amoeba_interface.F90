#include "copyright.i"

!*******************************************************************************!
! Module: amoeba_interface_mod
!
! Description: <TBS>
!
!*******************************************************************************

module amoeba_interface_mod

  implicit none

  private

  public        am_val_eval
  public        am_nonbond_eval
 
contains

!*******************************************************************************!
! Subroutine:  am_val_eval
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_val_eval(crd, frc, sander_vir, ebond, eangle, etors)

  use amoeba_bonds_mod, only : am_bonds_eval
  use amoeba_ureyb_mod, only : am_ureyb_eval
  use amoeba_reg_angles_mod, only : am_reg_angles_eval
  use amoeba_trig_angles_mod, only : am_trig_angles_eval
  use amoeba_opbend_angles_mod, only : am_opbend_angles_eval
  use amoeba_torsions_mod, only : am_torsions_eval
  use amoeba_stretch_torsions_mod, only : am_stretch_torsions_eval
  use amoeba_pitorsions_mod, only : am_pitorsions_eval
  use amoeba_stretch_bend_mod, only : am_stretch_bend_eval
  use amoeba_torsion_torsion_mod, only : am_tor_tor_eval
  use mdin_amoeba_dat_mod, only : do_amoeba_valence, amoeba_verbose
  use mdin_ctrl_dat_mod, only : iamoeba
  use file_io_dat_mod
  use timers_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: frc(3, *)
  double precision, intent(in out)      :: sander_vir(3)
  double precision, intent(in out)      :: ebond
  double precision, intent(in out)      :: eangle
  double precision, intent(in out)      :: etors

! Local variables:

  double precision                      :: ene(10)
  double precision                      :: vir_tensor(3, 3)

  ene = 0.d0
  vir_tensor(:, :) = 0.d0

  if (master) then      ! temporary conditional!
  if (iamoeba .eq. 1 .and. do_amoeba_valence .eq. 1) then
    call am_bonds_eval(crd, frc, ene(1), vir_tensor)
    call am_ureyb_eval(crd, frc, ene(2), vir_tensor)
   
    call update_time(bond_time)

    call am_reg_angles_eval(crd, frc, ene(3), vir_tensor)
    call am_trig_angles_eval(crd, frc, ene(4), vir_tensor)
    call am_opbend_angles_eval(crd, frc, ene(5), vir_tensor)

    call update_time(angle_time)

    call am_torsions_eval(crd, frc, ene(6), vir_tensor)
    call am_pitorsions_eval(crd, frc, ene(7), vir_tensor)
    call am_stretch_bend_eval(crd, frc, ene(8), vir_tensor)
    call am_tor_tor_eval(crd, frc, ene(9), vir_tensor)
    call am_stretch_torsions_eval(crd, frc, ene(10), vir_tensor)

    call update_time(dihedral_time)
  end if
  end if                ! temporary conditional!
  
  if (amoeba_verbose .ne. 0) then
     if( numtasks .gt. 1) then
        dbl_mpi_send_buf(1:10) = ene(1:10)
        dbl_mpi_send_buf(11:13) = vir_tensor(:,1)
        dbl_mpi_send_buf(14:16) = vir_tensor(:,2)
        dbl_mpi_send_buf(17:19) = vir_tensor(:,3)

        call mpi_reduce(dbl_mpi_send_buf, dbl_mpi_recv_buf, 19, &
                         mpi_double_precision, mpi_sum, 0, lib_mpi_comm, &
                         err_code_mpi)
     endif
     
     if (master) then
       if( numtasks .gt. 1)then
         ene(1:10) = dbl_mpi_recv_buf(1:10)
         vir_tensor(:,1) = dbl_mpi_recv_buf(11:13)
         vir_tensor(:,2) = dbl_mpi_recv_buf(14:16)
         vir_tensor(:,3) = dbl_mpi_recv_buf(17:19)
       endif 
       write(mdout, '(a, 3(1x, e16.8))') &
         'valence energies: bond,ureyb,angle ', ene(1), ene(2), ene(3)

       write(mdout, '(a, 3(1x, f14.4))') &
         'valence energies: trangle,opbend,tor ', ene(4), ene(5), ene(6)

       write(mdout, '(a, 3(1x, f14.4))') &
         'valence energies: pitors,strbend,tortor ', ene(7), ene(8), ene(9)

       write(mdout, '(a, 1x, f14.4)') &
         'valence energies: strtor ', ene(10)

       write(mdout, '(a, 3(1x, f14.4))') &
         'valence vir = ', vir_tensor(1, 1), vir_tensor(1, 2), vir_tensor(1, 3)

       write(mdout, '(a, 3(1x, f14.4))') &
         'valence vir = ', vir_tensor(2, 1), vir_tensor(2, 2), vir_tensor(2, 3)

       write(mdout, '(a, 3(1x, f14.4))') &
         'valence vir = ', vir_tensor(3, 1), vir_tensor(3, 2), vir_tensor(3, 3)

      if( numtasks .gt. 1) then
      ! Restore original single process values in master.  This will
      ! be reduced again later.  This is amoeba_verbose-only code (debug).
          ene(1:10) = dbl_mpi_send_buf(1:10)
          vir_tensor(:,1) = dbl_mpi_send_buf(11:13)
          vir_tensor(:,2) = dbl_mpi_send_buf(14:16)
          vir_tensor(:,3) = dbl_mpi_send_buf(17:19)
       endif
    end if
  end if

  ebond = ene(1) ! bonds energy
  eangle = ene(2) + ene(3) + ene(4) + ene(5) + ene(8) ! angle energy
  etors = ene(6) + ene(7) + ene(9) + ene(10) ! torsion energy
  sander_vir(1) = sander_vir(1) + vir_tensor(1, 1)
  sander_vir(2) = sander_vir(2) + vir_tensor(2, 2)
  sander_vir(3) = sander_vir(3) + vir_tensor(3, 3)
  
  return

end subroutine am_val_eval

!*******************************************************************************!
! Subroutine:  am_nonbond_eval
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_nonbond_eval(atm_cnt, crd, frc, img_frc, sander_vir, evdw, eelt, &
                           epolar, evdw_14, diprms, dipiter, netfrcs, atm_owner_map, &
                           tranvec )

  use mdin_amoeba_dat_mod, only : do_amoeba_nonbond, amoeba_verbose
  use mdin_ctrl_dat_mod, only : iamoeba

  use amoeba_multipoles_mod, only : am_mpole_local_to_global, &
                                    torque_field, &
                                    global_multipole, &
                                    am_mpole_torque_to_force

  use amoeba_induced_mod, only : am_induced_eval_uni, am_induced_eval_mpi
  use amoeba_adjust_mod, only : num_adjust_list
  use amoeba_recip_mod, only : amoeba_recip_step_alloc_dat, &
                               amoeba_recip_step_dealloc_dat
  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: atm_cnt
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in out)      :: frc(3, *)
  double precision, intent(in out)      :: img_frc(3, *)
  double precision, intent(in out)      :: sander_vir(3)
  double precision, intent(out)         :: evdw
  double precision, intent(out)         :: eelt
  double precision, intent(out)         :: epolar
  double precision, intent(out)         :: evdw_14
  double precision, intent(out)         :: dipiter
  double precision, intent(out)         :: diprms
  double precision, intent(out)         :: netfrcs(3)
  integer, intent(in out)               :: atm_owner_map(atm_cnt)
  double precision, intent(in)          :: tranvec(*)       ! translational vectors for atom images in different partition (cit) buckets (17 different variants)  

! Local variables:

  integer                               :: i,j, k
  integer                               :: num_ints, num_reals  ! not yet used.
  double precision                      :: vir_tensor(3, 3)

  num_ints = 0
  num_reals = 0

  evdw = 0.d0
  eelt = 0.d0
  epolar = 0.d0
  evdw_14 = 0.d0
  vir_tensor(:,:) = 0.d0
  netfrcs(:) = 0.d0             ! In case the code below does not execute...

  if (iamoeba .eq. 1 .and. do_amoeba_nonbond .eq. 1) then

    torque_field(1:10,1:atm_cnt) = 0.d0
    global_multipole(1:10,1:atm_cnt) = 0.d0

    call am_mpole_local_to_global(crd)

    call amoeba_recip_step_alloc_dat(num_ints, num_reals)

    if( numtasks .gt. 1) then
        call am_induced_eval_mpi(atm_cnt, crd, diprms, dipiter, num_adjust_list, & 
                                 tranvec)
    else
        call am_induced_eval_uni(atm_cnt, crd, diprms, dipiter, num_adjust_list, &
                                 tranvec)
    endif

    ! The frc array must be zero'd prior to entry into am_nonbond_ene_frc(),
    ! or the netfrcs calculation will be incorrect.
    call am_nonbond_ene_frc(atm_cnt, crd, eelt, epolar, evdw, &
                            evdw_14, frc, img_frc, vir_tensor, netfrcs, atm_owner_map, &
                            tranvec )

    ! Add the torque contributions
    
!    do i = 1,atm_cnt
!        frc(1,i) = 0.0d0   ! TMP DEBUG IGOR
!        frc(2,i) = 0.0d0   ! TMP DEBUG IGOR 
!        frc(3,i) = 0.0d0   ! TMP DEBUG IGOR
!    enddo

    call am_mpole_torque_to_force(atm_cnt, crd, frc, vir_tensor)
    
!    write(mdout,*)" Electric torque Forces on atoms "  !  TMP DEBUG IGOR
!    do i = 1,atm_cnt
!         write(mdout, '(i5,f20.10,f20.10,f20.10)')i,frc(1,i),frc(2,i),frc(3,i)
!    enddo

    ! call  dump_dipoles(frc, atm_cnt, 40)
    ! if (atm_cnt > 0) stop

    if (amoeba_verbose .ne. 0) then
      if( numtasks .gt. 1) then
        dbl_mpi_send_buf(1:3) = vir_tensor(:,1)
        dbl_mpi_send_buf(4:6) = vir_tensor(:,2)
        dbl_mpi_send_buf(7:9) = vir_tensor(:,3)

        call mpi_reduce(dbl_mpi_send_buf, dbl_mpi_recv_buf, 9, &
                        mpi_double_precision, mpi_sum, 0, lib_mpi_comm, &
                        err_code_mpi)
      endif 

      if (master) then
        if( numtasks .gt. 1) then
          vir_tensor(:,1) = dbl_mpi_recv_buf(1:3)
          vir_tensor(:,2) = dbl_mpi_recv_buf(4:6)
          vir_tensor(:,3) = dbl_mpi_recv_buf(7:9)
        endif 
        write(mdout, '(a, 3(1x, g16.8))') &
          'nonbond vir = ', vir_tensor(1, 1), vir_tensor(1, 2), vir_tensor(1, 3)

        write(mdout, '(a, 3(1x, g16.8))') &
          'nonbond vir = ', vir_tensor(2, 1), vir_tensor(2, 2), vir_tensor(2, 3)

        write(mdout, '(a, 3(1x, g16.8))') &
          'nonbond vir = ', vir_tensor(3, 1), vir_tensor(3, 2), vir_tensor(3, 3)

        if( numtasks .gt. 1) then
        ! Restore original single process values in master.  This will
        ! be reduced again later.  This is amoeba_verbose-only code (debug).
            vir_tensor(:,1) = dbl_mpi_send_buf(1:3)
            vir_tensor(:,2) = dbl_mpi_send_buf(4:6)
            vir_tensor(:,3) = dbl_mpi_send_buf(7:9)
        endif
      end if
    end if
    
    call amoeba_recip_step_dealloc_dat(num_ints, num_reals)

  end if

  sander_vir(1) = sander_vir(1) + vir_tensor(1, 1)
  sander_vir(2) = sander_vir(2) + vir_tensor(2, 2)
  sander_vir(3) = sander_vir(3) + vir_tensor(3, 3)
   
  return

end subroutine am_nonbond_eval

end module amoeba_interface_mod

!*******************************************************************************!
! Subroutine:  am_nonbond_perm_fields
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_nonbond_perm_fields(atm_cnt, is_polarizable, crd, &
                                  dip_field_d, dip_field_p, adj_dip_dip_tensor, &
                                  tranvec )

  use amoeba_recip_mod, only : am_recip_perm_field
  use amoeba_direct_mod, only : am_direct_permfield
  use amoeba_adjust_mod, only : am_adjust_permfield
  use amoeba_self_mod, only : am_self_permfield
  use img_mod
  use nb_pairlist_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  logical, intent(in)           :: is_polarizable(*)
  double precision, intent(in)  :: crd(3, *)
  double precision, intent(out) :: dip_field_d(3, atm_cnt)
  double precision, intent(out) :: dip_field_p(3, atm_cnt)
  double precision, intent(out) :: adj_dip_dip_tensor(6, *)
  double precision, intent(in)  :: tranvec(*)

! Local variables:

  double precision              :: cart_dipole_field(3, atm_cnt) ! electrical field on polarizable centers (cartesian coordinates) from reciprocal space part of PME sum  
  integer                       :: i
  
  cart_dipole_field(:,:) = 0.d0

  call am_recip_perm_field(atm_cnt, crd, cart_dipole_field)
    
  call am_direct_permfield(cart_dipole_field, gbl_ipairs, gbl_img, &
                           tranvec, gbl_img_atm_map)
                           
  call am_adjust_permfield(crd, dip_field_d, dip_field_p, adj_dip_dip_tensor)
  
!  write(*,*) " "
!  write(*,*) " am_nonbond_perm_fields() pt 1  dip_field_d: "
!  do i=1, atm_cnt
!     write(*,'(I3,3F14.6)') i, dip_field_d(1,i),dip_field_d(2,i),dip_field_d(3,i)
!  enddo
!  write(*,*) " "
!  write(*,*) " am_nonbond_perm_fields() pt 1  dip_field_p: "
!  do i=1, atm_cnt
!     write(*,'(I3,3F14.6)') i, dip_field_p(1,i),dip_field_p(2,i),dip_field_p(3,i)
!  enddo
  
  call am_self_permfield(atm_cnt, dip_field_d, dip_field_p)
  
  call am_induced_add_cart_to_dip(atm_cnt, is_polarizable, cart_dipole_field, &
                                  dip_field_d, dip_field_p)
  
!  write(*,*) " "
!  write(*,*) " am_nonbond_perm_fields() pt end  dip_field_d: "
!  do i=1, atm_cnt
!     write(*,'(I3,3F14.6)') i, dip_field_d(1,i),dip_field_d(2,i),dip_field_d(3,i)
!  enddo

!  write(*,*) " "
!  write(*,*) " am_nonbond_perm_fields() pt end  dip_field_p: "
!  do i=1, atm_cnt
!     write(*,'(I3,3F14.6)') i, dip_field_p(1,i),dip_field_p(2,i),dip_field_p(3,i)
!  enddo
  return

end subroutine am_nonbond_perm_fields

!*******************************************************************************!
! Subroutine:  am_nonbond_dip_dip_fields
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_nonbond_dip_dip_fields(atm_cnt, crd, ind_dip_d, ind_dip_p, &
                                     adj_dip_dip_tensor, &
                                     dip_field_d, dip_field_p)

  use amoeba_recip_mod, only : am_recip_dipole_field
  use amoeba_direct_mod, only : am_direct_dip_dip_field
  use amoeba_adjust_mod, only : am_adjust_dip_dip_fields
  use amoeba_self_mod, only : am_self_dipole_field

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  double precision, intent(in)  :: crd(3, *)
  double precision, intent(in)  :: ind_dip_d(3, *)
  double precision, intent(in)  :: ind_dip_p(3, *)
  double precision, intent(in)  :: adj_dip_dip_tensor(6, *)
  double precision, intent(out) :: dip_field_d(3, atm_cnt)
  double precision, intent(out) :: dip_field_p(3, atm_cnt)

  dip_field_d(:,:) = 0.d0
  dip_field_p(:,:) = 0.d0

  call am_recip_dipole_field(atm_cnt, crd, ind_dip_d, ind_dip_p, &
                             dip_field_d, dip_field_p)

  call am_direct_dip_dip_field(ind_dip_d, ind_dip_p, dip_field_d, dip_field_p)

  call am_adjust_dip_dip_fields(adj_dip_dip_tensor, ind_dip_d, ind_dip_p, &
                                dip_field_d, dip_field_p)

  call am_self_dipole_field(atm_cnt, ind_dip_d, ind_dip_p, dip_field_d, &
                            dip_field_p)

  return

end subroutine am_nonbond_dip_dip_fields

!*******************************************************************************!
! Subroutine:  am_nonbond_ene_frc
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_nonbond_ene_frc(atm_cnt, crd, ene_perm, ene_ind, ene_vdw, &
                              ene_vdw_14, frc, img_frc, vir_tensor, netfrcs, atm_owner_map, &
                              tranvec )

  use amoeba_recip_mod, only : am_recip_ene_frc, gbl_fractional_multipole, gbl_perm_F_field
  use amoeba_direct_mod, only : am_direct_ene_frc
  use amoeba_adjust_mod, only : am_adjust_ene_frc
  use amoeba_self_mod, only : am_self_ene_torque
  use amoeba_vdw_mod, only : am_vdw_longrange_ene
  use amoeba_induced_mod, only : ind_dip_d, ind_dip_p, am_var_polar_frc
  use amoeba_multipoles_mod, only : coulomb_const_kcal_per_mole,global_multipole,torque_field
  use mdin_amoeba_dat_mod, only : amoeba_verbose
  use file_io_dat_mod
  use img_mod
  use nb_pairlist_mod
  use parallel_dat_mod
  use pme_fft_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: atm_cnt
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(out)         :: ene_perm
  double precision, intent(out)         :: ene_ind
  double precision, intent(out)         :: ene_vdw
  double precision, intent(out)         :: ene_vdw_14
  double precision, intent(in out)      :: frc(3, *)
  double precision, intent(in out)      :: img_frc(3, *)
  double precision, intent(in out)      :: vir_tensor(3, 3)
  double precision, intent(out)         :: netfrcs(3)
  integer, intent(in out)               :: atm_owner_map(atm_cnt)
  double precision, intent(in)          :: tranvec(*)

! Local variables:

  double precision                      :: e_rec_perm
  double precision                      :: e_rec_ind
  double precision                      :: e_dir_perm
  double precision                      :: e_dir_ind
  double precision                      :: e_adj_perm
  double precision                      :: e_adj_ind
  double precision                      :: e_self_perm
  double precision                      :: e_self_ind
  double precision                      :: e_dir_vdw
  double precision                      :: e_adj_vdw
  double precision                      :: e_rec_vdw
  double precision                      :: e_hpol_ind  ! hyperpolarizable induced electrostatic energy
  
  integer                               :: i,j

  ! NOTE: There is an assumption here that on entry the img_frc array is zero'd.
  !       If it is not, the net force calculation will be incorrect!
  
!  do i = 1,atm_cnt
!    img_frc(:, i) = 0.d0 ! IGOR TMP DEBUG
!    frc(:, i) = 0.d0 ! IGOR TMP DEBUG 
!  end do
  
  call am_recip_ene_frc(atm_cnt, crd, ind_dip_d, ind_dip_p, e_rec_perm, &
                        e_rec_ind, img_frc, vir_tensor, netfrcs)
  
  call am_direct_ene_frc(gbl_ipairs, gbl_img, tranvec, crd, &
                         ind_dip_d, ind_dip_p, e_dir_perm, e_dir_ind, &
                         e_dir_vdw, frc, img_frc, vir_tensor, gbl_img_atm_map)
  
  call am_adjust_ene_frc(crd, ind_dip_d, ind_dip_p, &
                         e_adj_perm, e_adj_ind, e_adj_vdw, frc, vir_tensor, atm_owner_map )
  
  call am_self_ene_torque(atm_cnt, ind_dip_d, ind_dip_p, &  ! IGOR TMP DEBUG
                          e_self_perm, e_self_ind)          ! IGOR TMP DEBUG
  
  if (master) then
    call am_vdw_longrange_ene(e_rec_vdw, vir_tensor)    ! IGOR TMP DEBUG
  else
    e_rec_vdw = 0.d0
  end if

  e_hpol_ind = 0.0d0
#ifdef VAR_POLAR
!  call am_var_polar_frc(gbl_ipairs, gbl_img, tranvec, crd, e_hpol_ind, &   ! IGOR TMP DEBUG
!                       frc, img_frc, vir_tensor, gbl_img_atm_map)          ! IGOR TMP DEBUG
#endif

!  amoeba_verbose = 1         ! IGOR TMP DEBUG
  if (amoeba_verbose .ne. 0) then
    if( numtasks .gt. 1) then
        dbl_mpi_send_buf(1) =  e_rec_perm
        dbl_mpi_send_buf(2) =  e_dir_perm
        dbl_mpi_send_buf(3) =  e_adj_perm
        dbl_mpi_send_buf(4) =  e_self_perm
        dbl_mpi_send_buf(5) =  e_rec_ind
        dbl_mpi_send_buf(6) =  e_dir_ind
        dbl_mpi_send_buf(7) =  e_adj_ind
        dbl_mpi_send_buf(8) =  e_self_ind
        dbl_mpi_send_buf(9) =  e_dir_vdw
        dbl_mpi_send_buf(10) = e_adj_vdw
        dbl_mpi_send_buf(11) = e_rec_vdw
        dbl_mpi_send_buf(12) = e_hpol_ind

        call mpi_reduce(dbl_mpi_send_buf, dbl_mpi_recv_buf, 12, &
                        mpi_double_precision, mpi_sum, 0, lib_mpi_comm, &
                        err_code_mpi)
    endif

    if (master) then
      if( numtasks .gt. 1) then
        e_rec_perm = dbl_mpi_recv_buf(1)
        e_dir_perm = dbl_mpi_recv_buf(2)
        e_adj_perm = dbl_mpi_recv_buf(3)
        e_self_perm = dbl_mpi_recv_buf(4)
        e_rec_ind = dbl_mpi_recv_buf(5)
        e_dir_ind = dbl_mpi_recv_buf(6)
        e_adj_ind = dbl_mpi_recv_buf(7)
        e_self_ind = dbl_mpi_recv_buf(8)
        e_dir_vdw = dbl_mpi_recv_buf(9)
        e_adj_vdw = dbl_mpi_recv_buf(10)
        e_rec_vdw = dbl_mpi_recv_buf(11)
        e_hpol_ind = dbl_mpi_recv_buf(12)
      endif
      write(mdout, '(a, /, 4(1x, f20.10))') &
        'e_rec_perm,e_dir_perm,e_adj_perm,e_self_perm = ', &
        e_rec_perm, e_dir_perm, e_adj_perm, e_self_perm

      write(mdout, '(a, /, 5(1x, f20.10))') &
        'e_rec_ind,e_dir_ind,e_adj_ind,e_self_ind,e_hpol_ind  = ', &
        e_rec_ind, e_dir_ind, e_adj_ind, e_self_ind, e_hpol_ind

      write(mdout, '(a, /, 3(1x, f20.10))') &
        'e_dir_vdw,e_adj_vdw,e_rec_vdw = ', &
        e_dir_vdw, e_adj_vdw, e_rec_vdw
 
! IGOR TMP DEBUG:
!      write(mdout,'(a,i5)') "atm_cnt = ",atm_cnt
!      write(mdout,'(a,i5,i5,i5)')" FFT dimensions: ",fft_x_dim,fft_y_dim,fft_z_dim
      
!      write(mdout,*)" Multipoles of oxygen "
!      do j = 4,4
!        do i = 1,10
!         write(mdout, '(i5,i5,f22.12)')i,j,global_multipole(i,j)
!        enddo
!      enddo
!      write(mdout,*)" Fractional Multipoles of oxygen  (image 5) "
!      do j = 5,5
!        do i = 1,10
!         write(mdout, '(i5,i5,f22.12)')i,j,gbl_fractional_multipole(i,j)
!        enddo
!      enddo
!      write(mdout,*)" Components of the electric field on oxygen (image 5)"
!      do j = 5,5
!        do i = 1,10
!         write(mdout, '(i5,i5,f22.12)')i,j,gbl_perm_F_field(i,j)
!        enddo
!      enddo
!      write(mdout,*)" Forces on atoms "
!      do i = 1,6
!         write(mdout, '(i5,f20.10,f20.10,f20.10)')i,frc(1,i),frc(2,i),frc(3,i)
!      enddo
!      write(mdout,*)" Forces on atoms (through its images) "
!      do i = 1,6
!         j = gbl_atm_img_map(i)
!         write(mdout, '(i5,f20.10,f20.10,f20.10)')i,img_frc(1,j),img_frc(2,j),img_frc(3,j)
!      enddo
!      write(mdout,*)" Electric Forces on atoms (total) "
!      do i = 1,6
!         j = gbl_atm_img_map(i)
!         write(mdout, '(i5,f20.10,f20.10,f20.10)')i,(img_frc(1,j)+frc(1,i)),(img_frc(2,j)+frc(2,i)),(img_frc(3,j)+frc(3,i))
!      enddo
!      write(mdout,*)"Torque field zero on atoms 0-order"
!      do i = 1,6
!         j = gbl_atm_img_map(i)
!         write(mdout, '(i5,f20.10,f20.10,f20.10)')i,torque_field(1,i)
!      enddo
!      write(mdout,*)"Torque field zero on atoms 1-order"
!      do i = 1,6
!         j = gbl_atm_img_map(i)
!         write(mdout, '(i5,f20.10,f20.10,f20.10)')i,torque_field(2,i),torque_field(3,i),torque_field(4,i)
!      enddo
!      write(mdout,*)"Torque field zero on atoms 2-order"
!      do i = 1,6
!         j = gbl_atm_img_map(i)
!         write(mdout, '(i5,f20.10,f20.10,f20.10,f20.10,f20.10,f20.10)')i,coulomb_const_kcal_per_mole*torque_field(5,i),coulomb_const_kcal_per_mole*torque_field(6,i),coulomb_const_kcal_per_mole*torque_field(7,i),coulomb_const_kcal_per_mole*torque_field(8,i),coulomb_const_kcal_per_mole*torque_field(9,i),coulomb_const_kcal_per_mole*torque_field(10,i)
!      enddo
      
      if( numtasks .gt. 1) then
        ! Restore original single process values in master.  This will
        ! be reduced again later.  This is amoeba_verbose-only code (debug).
        e_rec_perm = dbl_mpi_send_buf(1)
        e_dir_perm = dbl_mpi_send_buf(2)
        e_adj_perm = dbl_mpi_send_buf(3)
        e_self_perm = dbl_mpi_send_buf(4)
        e_rec_ind = dbl_mpi_send_buf(5)
        e_dir_ind = dbl_mpi_send_buf(6)
        e_adj_ind = dbl_mpi_send_buf(7)
        e_self_ind = dbl_mpi_send_buf(8)
        e_dir_vdw = dbl_mpi_send_buf(9)
        e_adj_vdw = dbl_mpi_send_buf(10)
        e_rec_vdw = dbl_mpi_send_buf(11)
        e_hpol_ind = dbl_mpi_send_buf(12)
      endif
    end if
  end if

  ene_perm = e_rec_perm + e_dir_perm + e_adj_perm + e_self_perm
  ene_ind = e_rec_ind + e_dir_ind + e_adj_ind + e_self_ind + e_hpol_ind

  ene_vdw = e_dir_vdw + e_rec_vdw
  ene_vdw_14 = e_adj_vdw

  return

end subroutine am_nonbond_ene_frc

!*******************************************************************************!
! Subroutine:  am_nonbond_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_nonbond_set_user_bit(do_recip, do_adjust, do_direct, do_self, &
                                   do_vdw, do_induce)

  use amoeba_recip_mod, only : am_recip_zero_flag, am_recip_set_user_bit
  use amoeba_adjust_mod, only : am_adjust_zero_flag, am_adjust_set_user_bit
  use amoeba_direct_mod, only : am_direct_zero_flag, am_direct_set_user_bit
  use amoeba_self_mod, only : am_self_zero_flag, am_self_set_user_bit
  use amoeba_vdw_mod, only : am_vdw_zero_flag, am_vdw_set_user_bit
  use amoeba_induced_mod, only : am_induced_zero_flag, am_induced_set_user_bit

  implicit none

! Formal arguments:

  integer, intent(in)  :: do_recip
  integer, intent(in)  :: do_adjust
  integer, intent(in)  :: do_direct
  integer, intent(in)  :: do_self
  integer, intent(in)  :: do_vdw
  integer, intent(in)  :: do_induce
                          
  call am_recip_zero_flag
  call am_adjust_zero_flag
  call am_direct_zero_flag
  call am_self_zero_flag
  call am_vdw_zero_flag
  call am_induced_zero_flag

  call am_recip_set_user_bit(do_recip)
  call am_adjust_set_user_bit(do_adjust)
  call am_direct_set_user_bit(do_direct)
  call am_self_set_user_bit(do_self)
  call am_vdw_set_user_bit(do_vdw)
  call am_induced_set_user_bit(do_induce)

  return

end subroutine am_nonbond_set_user_bit

!*******************************************************************************!
! Subroutine:  am_val_set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************

subroutine am_val_set_user_bit(do_bond, do_ureyb, do_reg_angle, &
                               do_trig_angle, do_opbend_angle, do_torsions, &
                               do_str_torsions, do_pitorsions,  &
                               do_stretch_bend, do_torsion_torsion)

  use amoeba_bonds_mod, only : am_bonds_zero_flag, &
                               am_bonds_set_user_bit
  use amoeba_ureyb_mod, only : am_ureyb_zero_flag, &
                               am_ureyb_set_user_bit
  use amoeba_reg_angles_mod, only : am_reg_angles_zero_flag, &
                                    am_reg_angles_set_user_bit
  use amoeba_trig_angles_mod, only : am_trig_angles_zero_flag, &
                                     am_trig_angles_set_user_bit
  use amoeba_opbend_angles_mod, only : am_opbend_angles_zero_flag, &
                                       am_opbend_angles_set_user_bit
  use amoeba_torsions_mod, only : am_torsions_zero_flag, &
                                  am_torsions_set_user_bit
  use amoeba_stretch_torsions_mod, only : am_stretch_torsions_zero_flag, &
                                          am_stretch_torsions_suser_bit
  use amoeba_pitorsions_mod, only : am_pitorsions_zero_flag, &
                                    am_pitorsions_set_user_bit
  use amoeba_stretch_bend_mod, only : am_stretch_bend_zero_flag, &
                                      am_stretch_bend_set_user_bit
  use amoeba_torsion_torsion_mod, only : am_tor_tor_zero_flag, &
                                         am_tor_tor_set_user_bit

  implicit none

! Formal arguments:

  integer, intent(in)   :: do_bond
  integer, intent(in)   :: do_ureyb
  integer, intent(in)   :: do_reg_angle
  integer, intent(in)   :: do_trig_angle
  integer, intent(in)   :: do_opbend_angle
  integer, intent(in)   :: do_torsions
  integer, intent(in)   :: do_pitorsions
  integer, intent(in)   :: do_str_torsions
  integer, intent(in)   :: do_stretch_bend
  integer, intent(in)   :: do_torsion_torsion

  call am_bonds_zero_flag
  call am_ureyb_zero_flag
  call am_reg_angles_zero_flag
  call am_trig_angles_zero_flag
  call am_opbend_angles_zero_flag
  call am_torsions_zero_flag
  call am_stretch_torsions_zero_flag
  call am_pitorsions_zero_flag
  call am_stretch_bend_zero_flag
  call am_tor_tor_zero_flag

  call am_bonds_set_user_bit(do_bond)
  call am_ureyb_set_user_bit(do_ureyb)
  call am_reg_angles_set_user_bit(do_reg_angle)
  call am_trig_angles_set_user_bit(do_trig_angle)
  call am_opbend_angles_set_user_bit(do_opbend_angle)
  call am_torsions_set_user_bit(do_torsions)
  call am_stretch_torsions_suser_bit(do_str_torsions)
  call am_pitorsions_set_user_bit(do_pitorsions)
  call am_stretch_bend_set_user_bit(do_stretch_bend)
  call am_tor_tor_set_user_bit(do_torsion_torsion)

  return

end subroutine am_val_set_user_bit


!*******************************************************************************!
! Subroutine:  array_add
!
! Description: <TBS>
!
!*******************************************************************************

subroutine array_add(a, b, num)

  implicit none

! Formal arguments:

  double precision, intent(in out)      :: a(*)
  double precision, intent(in)          :: b(*)
  integer, intent(in)                   :: num

! Local variables:

  integer                               :: n

  do n = 1, num
    a(n) = a(n) + b(n)
  end do

  return

end subroutine array_add

!*******************************************************************************!
! Subroutine:  array_copy
!
! Description: <TBS>
!
!*******************************************************************************

subroutine array_copy(a, b, num)

  implicit none

! Formal arguments:

  double precision, intent(in out)      :: a(1:num)
  double precision, intent(in)          :: b(1:num)
  integer, intent(in)                   :: num

  a(:) = b(:)

  return

end subroutine array_copy

!*******************************************************************************!
! Subroutine:  dump_3d_array
!
! Description: <TBS>
!
!*******************************************************************************

subroutine dump_3d_array(array, cnt, nf)

  implicit none

! Formal arguments:

  double precision      :: array(3, *)
  integer               :: cnt, nf

! Local variables:

  integer               :: j, n

  do n = 1, cnt
    write(nf, '(3f20.12)') (array(j, n), j = 1, 3)
  end do

  return

end subroutine dump_3d_array

!*******************************************************************************!
! Subroutine:  reduce_dump_3d_array
!
! Description: <TBS>
!
!*******************************************************************************

subroutine reduce_dump_3d_array(array, cnt, nf)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision      :: array(3, *)
  integer               :: cnt, nf

! Local variables:

  integer               :: j, n
  double precision      :: outbuf(3, cnt)

  call mpi_allreduce(array, dbl_mpi_recv_buf, 3 * cnt, &
                     mpi_double_precision, mpi_sum, lib_mpi_comm, &
                     err_code_mpi)
  call array_copy(outbuf, dbl_mpi_recv_buf, 3 * cnt)

  if (master) then
    do n = 1, cnt
      write(nf, '(3f20.12)') (outbuf(j, n), j = 1, 3)
    end do
  end if

  return

end subroutine reduce_dump_3d_array


