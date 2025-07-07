#include "copyright.i"

!*******************************************************************************
!
! Module: amoeba_force_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module amoeba_force_mod

  use gbl_datatypes_mod

  implicit none

  ! Potential energies, with breakdown, from pme.  This is intended to be the
  ! external interface to potential energy data produced by this module, in
  ! particular by subroutine amoeba_force()

  type amba_pot_ene_rec
    sequence
    double precision    :: total        ! EPtot
    double precision    :: vdw          ! VDWAALS
    double precision    :: elec         ! EELEC
    double precision    :: hbond        ! EHBOND, always 0.d0
    double precision    :: bond         ! BOND
    double precision    :: angle        ! ANGLE
    double precision    :: dihedral     ! DIHED
    double precision    :: vdw_14       ! 1-4 NB
    double precision    :: elec_14      ! 1-4 EEL, always 0.d0
    double precision    :: restraint    ! RESTRAINT
    double precision    :: polar        ! EPOLZ
  end type amba_pot_ene_rec

  integer, parameter    :: amba_pot_ene_rec_size = 11

  type(amba_pot_ene_rec), parameter      :: null_amba_pot_ene_rec = &
    amba_pot_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)

contains

!*******************************************************************************
!
! Subroutine:  amoeba_force
!
! Description: <TBS>
!              
!*******************************************************************************

! MPI IMPLEMENTATION!
subroutine amoeba_force_mpi(atm_cnt, crd, saved_crd, box, saved_box, frc, mass, img_atm_map, atm_img_map, &
                            atm_owner_map, my_atm_lst, my_atm_cnt, new_list, pot_ene, diprms, dipiter, virial, &
                            atm_jrc, atm_xc, atm_weight, igroup, natc, tranvec, imin_par )

  use constraints_mod
  use dist_constr_mod
  use dynamics_dat_mod
  use dynamics_mod
  use pme_direct_mod
  use pme_recip_mod
  use img_mod
!  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_pairlist_mod
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use pme_force_mod
  use prmtop_dat_mod
  use runfiles_mod
  use timers_mod
  use amoeba_interface_mod
  use mdin_amoeba_dat_mod
  use mdin_ewald_dat_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3, atm_cnt)
  double precision              :: saved_crd(3, atm_cnt)
  double precision              :: box(3)
  double precision              :: saved_box(3)
  double precision              :: frc(3, atm_cnt)
  double precision              :: mass(atm_cnt)
  integer                       :: img_atm_map(atm_cnt)
  integer                       :: atm_img_map(atm_cnt)
  integer                       :: atm_owner_map(atm_cnt)
  integer                       :: my_atm_lst(*)        ! used only in MPI
  integer                       :: my_atm_cnt           ! used only in MPI
  logical                       :: new_list
  type(amba_pot_ene_rec)        :: pot_ene
  double precision              :: diprms
  double precision              :: dipiter
  double precision, optional    :: virial(3)            ! Only used for MD
  integer                       :: atm_jrc(*)
  double precision              :: atm_xc(3,*)
  double precision              :: atm_weight(*)
  integer                       :: igroup(*)
  integer                       :: natc
  double precision              :: nb_constr_ene
  double precision              :: tranvec(*)
  integer                       :: imin_par

! Local variables:

  double precision              :: netfrcs(3)
  double precision              :: virial_lcl(3)
  integer                       :: atm_lst_idx
  integer                       :: alloc_failed
  integer                       :: i, j

  double precision, allocatable :: img_frc(:,:)
#ifdef UNDEF
  double precision, allocatable :: nb_frc(:,:)
#else
  double precision, allocatable :: summed_frcs(:,:)
#endif

  call zero_time()
  call zero_pme_time()

! Zero energies that are stack or call parameters:

  pot_ene = null_amba_pot_ene_rec

  virial_lcl(:) = 0.d0

! Zero internal energies, virials, etc.

#ifdef UNDEF
  ! BUGBUG: nb_frc, img_frc not actually used yet...
  allocate(nb_frc(3, atm_cnt), &
           stat = alloc_failed)
#else
  allocate(img_frc(3, atm_cnt), &
           summed_frcs(3, atm_cnt), &
           stat = alloc_failed)
#endif /* UNDEF */
  if (alloc_failed .ne. 0) call setup_alloc_error
  
! Calculate the non-bonded contributions:

! Direct part of ewald plus vdw, hbond, pairlist setup and image claiming:

  if (ntp .gt. 0) call fill_tranvec(tranvec)

  call update_pme_time(pme_misc_timer)

  ! The following encapsulation (save_imgcrds) seems to be necessary to
  ! prevent an optimization bug with the SGI f90 compiler.  Sigh...


  if (new_list) then
    call pme_list(atm_cnt, crd, saved_crd, box, saved_box, igroup, gbl_excl_img_pairlst, atm_nb_maskdata, &
                  atm_nb_mask, atm_owner_map, my_atm_lst, my_atm_cnt, tranvec, imin_par )
    call save_imgcrds(atm_cnt, used_img_lo, used_img_hi, gbl_img, &
                        gbl_saved_imgcrd) 
  else
    call adjust_imgcrds(atm_cnt, used_img_lo, used_img_hi, gbl_img, &
                        img_atm_map, gbl_saved_imgcrd, crd, &
                        saved_crd, box, saved_box, ntp)
  end if
  
  if (used_img_range_wraps) then
    do i = used_img_lo, atm_cnt
      img_frc(1, i) = 0.d0
      img_frc(2, i) = 0.d0
      img_frc(3, i) = 0.d0
    end do
    do i = 1, used_img_hi
      img_frc(1, i) = 0.d0
      img_frc(2, i) = 0.d0
      img_frc(3, i) = 0.d0
    end do
  else
    do i = used_img_lo, used_img_hi
      img_frc(1, i) = 0.d0
      img_frc(2, i) = 0.d0
      img_frc(3, i) = 0.d0
    end do
  end if

  ! We also need to zero anything on the extra used atom list since nonbonded
  ! forces for it will get sent to the atom owner.  We don't calc any such
  ! forces for these atoms, but they are on the send atom list in order to
  ! get their coordinates updated.

  do j = 1, extra_used_atm_cnt
    i = atm_img_map(gbl_extra_used_atms(j))
    img_frc(:, i) = 0.d0
  end do
  
! BUGBUG - Because we are accumulating into frc instead of img_frc here,
!          we have to clear the frc array first.  Because we are using an
!          allreduce to send frc info around, we have to clear the entire
!          array...

  frc(:,:) = 0.d0       ! Temporary...  REMOVE LATER...
                        ! BUT NOTE - there is currently an assumption that the
                        !            img_frc array is entering am_nonbond_eval
                        !            zero'd - necessary for netfrcs calc to
                        !            be correct...  The netfrcs(:) array itself
                        !            will be zero'd in the call before being
                        !            summed, but the returned value is a frc
                        !            sum at this point.

  call am_nonbond_eval(atm_cnt, crd, frc, img_frc, virial_lcl, &
                       pot_ene%vdw, pot_ene%elec, pot_ene%polar, &
                       pot_ene%vdw_14, diprms, dipiter, netfrcs, atm_owner_map, &
                       tranvec )
                       
! BUGBUG - netfrcs must be accumulated globally and corrected for in the
!          parallel implementation! TBD!!!

! Clear the force array. We delay copying the nonbonded forces into the force
! array in order to be able to schedule i/o later, and batch it up.
! BUGBUG - NOT DONE THAT WAY YET...

#ifdef UNDEF
  do atm_lst_idx = 1, my_atm_cnt
    i = my_atm_lst(atm_lst_idx)
    frc(:,i) = 0.d0
  end do
#endif

! Calculate the other contributions:

  call update_time(nonbond_time)

  ! BUGBUG - For now, only the master does valence terms calcs; this is the
  !          case because we have to actually assign the workload based on
  !          individual valence term list structures if we are going to
  !          distribute this load.  The process selection must occur inside
  !          am_val_eval() to prevent hanging on verbose printout...

  call am_val_eval(crd, frc, virial_lcl, &
                   pot_ene%bond, pot_ene%angle, pot_ene%dihedral)
                   
!  write(iunit_debug,*) " amoeba_force_mpi() pt 1  pot_ene%bond = ",pot_ene%bond
!  write(iunit_debug,*) " amoeba_force_mpi() pt 1  pot_ene%polar = ",pot_ene%polar

! Calculate the position constraint energy:

  if (natc .gt. 0) then
    call get_crd_constraint_energy(natc, pot_ene%restraint, atm_jrc, &
                                   crd, frc, atm_xc, atm_weight, atm_owner_map )
  endif
  
  if( cit_num_dist_constr .gt. 0 ) then	 
	 call get_dist_constr_energy(cit_num_dist_constr, cit_dist_constr, crd, frc, nb_constr_ene, atm_owner_map)
	 pot_ene%restraint = pot_ene%restraint + nb_constr_ene
  endif

  ! Sum up total potential energy for this task:

  pot_ene%total = pot_ene%vdw + &
                  pot_ene%elec + &
                  pot_ene%hbond + &
                  pot_ene%bond + &
                  pot_ene%angle + &
                  pot_ene%dihedral + &
                  pot_ene%vdw_14 + &
                  pot_ene%elec_14 + &
                  pot_ene%restraint + &
                  pot_ene%polar
                               
  ! BUGBUG - The standard pmemd force distribution scheme is not yet ready for
  !          prime time in an amoeba context.

  call update_time(nonbond_time)
#ifdef UNDEF
  if (new_list) then
#ifdef SLOW_NONBLOCKING_MPI
    call get_img_frc_distribution(atm_cnt, gbl_atm_offsets, gbl_taskmap, &
                                  gbl_inv_taskmap, &
                                  gbl_send_atm_lst, gbl_send_atm_cnts, &
                                  gbl_recv_atm_lsts, gbl_recv_atm_cnts)
#else
    call get_img_frc_distribution(atm_cnt, gbl_atm_offsets, gbl_taskmap, &
                                  gbl_send_atm_lst, gbl_send_atm_cnts, &
                                  gbl_recv_atm_lsts, gbl_recv_atm_cnts)
#endif
  end if
  
  call distribute_img_frcs(atm_cnt, img_frc, nb_frc, atm_img_map, &
                           gbl_atm_offsets, gbl_taskmap, &
                           gbl_inv_taskmap, &
                           gbl_send_atm_lst, gbl_send_atm_cnts, &
                           gbl_recv_atm_lsts, gbl_recv_atm_cnts, &
                           dbl_mpi_send_buf, dbl_mpi_recv_buf)
#else
  ! BUGBUG - Temporary amoeba force distribution code.  Copy the img_frc
  ! values to frc and do an all_reduce of frc.
  
  if (used_img_range_wraps) then
    do i = used_img_lo, atm_cnt
      j = img_atm_map(i)
      if( j .gt. 0) frc(:, j) = frc(:, j) + img_frc(:, i)  ! IGOR TEMP FIX    do we just add img_frc to frc??
    end do
    do i = 1, used_img_hi
      j = img_atm_map(i)
      if( j .gt. 0) frc(:, j) = frc(:, j) + img_frc(:, i)  ! IGOR TEMP FIX 
    end do
  else
    do i = used_img_lo, used_img_hi
      j = img_atm_map(i)
      if( j .gt. 0) frc(:, j) = frc(:, j) + img_frc(:, i)
    end do
  end if

  call mpi_allreduce(frc, summed_frcs, 3 * atm_cnt, mpi_double_precision, &
                     mpi_sum, lib_mpi_comm, err_code_mpi)
  frc(:,:) = summed_frcs(:,:)
#endif /* UNDEF */

! Add pot_ene, ekcmt(1:3) (md only), and the bulk of the ew_ene_vir common
! block together from all nodes.

  call distribute_amoeba_enes_virs_netfrcs(pot_ene, virial_lcl, netfrcs )

  call update_time(fcve_dist_time)
  call zero_pme_time()

  ! BUGBUG - This IS the correct thing to do, given the current context, where
  !          forces are globally reduced as a force (not image force) array,
  !          the reduction has already been done, and now we are correcting
  !          the forces for the atoms we own.  Once this is in an image force
  !          context, it should change...

  if (netfrc .ne. 0) then
    netfrcs(:) = netfrcs(:) / atm_cnt
    ! Do net force corrections to your own atoms:
    do atm_lst_idx = 1, my_atm_cnt
      i = my_atm_lst(atm_lst_idx)
      frc(:, i) = frc(:, i) - netfrcs(:)
    end do
  end if

  ! Copy image forces to atom forces: 

#ifdef UNDEF
  ! BUGBUG - Not currently needed, but it will be needed...
  do atm_lst_idx = 1, my_atm_cnt
    i = my_atm_lst(atm_lst_idx)
    frc(1, i) = frc(1, i) + nb_frc(1, i)
    frc(2, i) = frc(2, i) + nb_frc(2, i)
    frc(3, i) = frc(3, i) + nb_frc(3, i)
  end do
#endif /* UNDEF */

! If belly is on then set the belly atom forces to zero:

!  if (ibelly .gt. 0) call bellyf(atm_cnt, igroup, frc)

#ifdef UNDEF
  ! BUGBUG - Not needed yet...
  deallocate(nb_frc)
#else
  deallocate(img_frc, summed_frcs)
#endif /* UNDEF */

  virial(:) = virial_lcl(:)

  call update_time(nonbond_time)
  call update_pme_time(pme_misc_timer)
  
!  write(iunit_debug,*) " amoeba_force_mpi() pt end  pot_ene%bond = ",pot_ene%bond
!  write(iunit_debug,*) " amoeba_force_mpi() pt end  pot_ene%polar = ",pot_ene%polar

  return

end subroutine amoeba_force_mpi
!*******************************************************************************
!
! Subroutine:  amoeba_force_uni
!
! Description: <TBS>
!              
!*******************************************************************************

! UNIPROCESSOR IMPLEMENTATION!
subroutine amoeba_force_uni(atm_cnt, crd, saved_crd, box, saved_box, frc, mass, img_atm_map, atm_img_map, &
                            atm_owner_map, my_atm_lst, my_atm_cnt, new_list, pot_ene, diprms, dipiter, virial, &
                            atm_jrc, atm_xc, atm_weight, igroup, natc, tranvec, imin_par )

  use constraints_mod
  use dist_constr_mod
  use dynamics_dat_mod
  use dynamics_mod
  use pme_direct_mod
  use pme_recip_mod
  use img_mod
!  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_pairlist_mod
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use pme_force_mod
  use prmtop_dat_mod
  use runfiles_mod
  use timers_mod
  use amoeba_interface_mod
  use mdin_amoeba_dat_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3, atm_cnt)
  double precision              :: saved_crd(3, atm_cnt)
  double precision              :: box(3)
  double precision              :: saved_box(3)
  double precision              :: frc(3, atm_cnt)
  double precision              :: mass(atm_cnt)
  integer                       :: img_atm_map(atm_cnt)
  integer                       :: atm_img_map(atm_cnt)
  integer                       :: atm_owner_map(atm_cnt)
  integer                       :: my_atm_lst(*)        ! used only in MPI
  integer                       :: my_atm_cnt           ! used only in MPI
  logical                       :: new_list
  type(amba_pot_ene_rec)        :: pot_ene
  double precision              :: diprms
  double precision              :: dipiter
  double precision, optional    :: virial(3)            ! Only used for MD
  integer                       :: atm_jrc(*)
  double precision              :: atm_xc(3,*)
  double precision              :: atm_weight(*)
  integer                       :: igroup(*)
  integer                       :: natc
  double precision              :: nb_constr_ene
  double precision              :: tranvec(*)
  integer                       :: imin_par

! Local variables:

  double precision              :: netfrcs(3)
  double precision              :: virial_lcl(3)
  integer                       :: alloc_failed
  integer                       :: i, j

  double precision, allocatable :: img_frc(:,:)

  call zero_time()
  call zero_pme_time()

! Zero energies that are stack or call parameters:

  pot_ene = null_amba_pot_ene_rec

  virial_lcl(:) = 0.d0

! Zero internal energies, virials, etc.

  allocate(img_frc(3, atm_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

! Calculate the non-bonded contributions:

! Direct part of ewald plus vdw, hbond, pairlist setup and image claiming:

  if (ntp .gt. 0) call fill_tranvec(tranvec)

  call update_pme_time(pme_misc_timer)
  
!  write(mdout,*) "Hello I am in amoeba_force()  pt 1"
!  open(11,file="ene_fort.dat",FORM='formatted',status='unknown')
!  write(11,*) "Hello I am in amoeba_force()  pt 2"
!  close(11)
!  call flush_mdout()

  ! The following encapsulation (save_imgcrds) seems to be necessary to
  ! prevent an optimization bug with the SGI f90 compiler.  Sigh...


  if (new_list) then
    call pme_list(atm_cnt, crd, saved_crd, box, saved_box, igroup, gbl_excl_img_pairlst, atm_nb_maskdata, &
                  atm_nb_mask, atm_owner_map, my_atm_lst, my_atm_cnt, tranvec, imin_par )
    call save_imgcrds(atm_cnt, used_img_lo, used_img_hi, gbl_img, &
                      gbl_saved_imgcrd)
  else
    call adjust_imgcrds(atm_cnt, used_img_lo, used_img_hi, gbl_img, &
                        img_atm_map, gbl_saved_imgcrd, crd, &
                        saved_crd, box, saved_box, ntp)
  end if

  ! Time for this gets lumped in am_nonbond_eval()...

  img_frc(:,:) = 0.d0

  frc(:,:) = 0.d0       ! Note there is an assumption that the img_frc array is
                        ! entering am_nonbond_eval zero'd - necessary for
                        ! netfrcs calc to be correct...  The netfrcs(:) array
                        ! itself will be zero'd in the call before being
                        ! summed, but the returned value is a frc sum at this
                        ! point.

  call am_nonbond_eval(atm_cnt, crd, frc, img_frc, virial_lcl, &
                       pot_ene%vdw, pot_ene%elec, pot_ene%polar, &
                       pot_ene%vdw_14, diprms, dipiter, netfrcs, atm_owner_map, & 
                       tranvec )

  if (netfrc .ne. 0) then
    netfrcs(:) = netfrcs(:) / atm_cnt

    do i = 1, atm_cnt
      frc(:, i) = frc(:, i) - netfrcs(:)
    end do
  end if

  do i = 1, atm_cnt
    j = atm_img_map(i)
    frc(1, i) = frc(1, i) + img_frc(1, j)
    frc(2, i) = frc(2, i) + img_frc(2, j)
    frc(3, i) = frc(3, i) + img_frc(3, j)
  end do

! Calculate the other contributions:

  call update_time(nonbond_time)
  call update_pme_time(pme_misc_timer)

  call am_val_eval(crd, frc, virial_lcl, &
                   pot_ene%bond, pot_ene%angle, pot_ene%dihedral)

! Calculate the position constraint energy:

  if (natc .gt. 0) then
    call get_crd_constraint_energy(natc, pot_ene%restraint, atm_jrc, &
                                   crd, frc, atm_xc, atm_weight, atm_owner_map )
  endif
  
  if( cit_num_dist_constr .gt. 0 ) then	 
	 call get_dist_constr_energy(cit_num_dist_constr, cit_dist_constr, crd, frc, nb_constr_ene, atm_owner_map)
	 pot_ene%restraint = pot_ene%restraint + nb_constr_ene
  endif

  ! Sum up total potential energy for this task:

  pot_ene%total = pot_ene%vdw + &
                  pot_ene%elec + &
                  pot_ene%hbond + &
                  pot_ene%bond + &
                  pot_ene%angle + &
                  pot_ene%dihedral + &
                  pot_ene%vdw_14 + &
                  pot_ene%elec_14 + &
                  pot_ene%restraint + &
                  pot_ene%polar
               
! If belly is on then set the belly atom forces to zero:

!  if (ibelly .gt. 0) call bellyf(atm_cnt, igroup, frc)

  deallocate(img_frc)

  call update_time(nonbond_time)
  call update_pme_time(pme_misc_timer)

  virial(:) = virial_lcl(:)
  
  return

end subroutine amoeba_force_uni

!*******************************************************************************
!
! Subroutine:  distribute_amoeba_enes_virs_netfrcs
!
! Description: We reduce the appropriate subset of values in the ene array,
!              and the pme_ene_vir common block
!*******************************************************************************

subroutine distribute_amoeba_enes_virs_netfrcs(pot_ene, virial, netfrcs )

  use parallel_dat_mod

  implicit none

! Formal arguments:

  type(amba_pot_ene_rec), intent(in out)        :: pot_ene
  double precision, intent(in out)              :: virial(3)
  double precision, intent(in out)              :: netfrcs(3)

! Local variables:

  type amba_dat
    sequence
    type(amba_pot_ene_rec)      :: pot_ene
    double precision            :: virial(3)
    double precision            :: netfrcs(3)
  end type amba_dat

  type(amba_dat), save          :: dat_in, dat_out
  integer                       :: buf_size

  dat_in%pot_ene = pot_ene
  dat_in%virial(:) = virial(:)
  dat_in%netfrcs(:) = netfrcs(:)

  buf_size = amba_pot_ene_rec_size + 3 + 3

  call mpi_allreduce(dat_in%pot_ene%total, dat_out%pot_ene%total, &
                     buf_size, mpi_double_precision, &
                     mpi_sum, lib_mpi_comm, err_code_mpi)

  pot_ene = dat_out%pot_ene
  virial(:) = dat_out%virial(:)
  netfrcs(:) = dat_out%netfrcs(:)

  return

end subroutine distribute_amoeba_enes_virs_netfrcs

end module amoeba_force_mod
