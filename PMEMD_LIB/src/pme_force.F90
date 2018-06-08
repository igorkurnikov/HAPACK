#include "copyright.i"

!*******************************************************************************
!
! Module: pme_force_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module pme_force_mod

  use gbl_datatypes_mod

  implicit none

  ! Potential energies, with breakdown, from pme.  This is intended to be the
  ! external interface to potential energy data produced by this module, in
  ! particular by subroutine pme_force().

  type pme_pot_ene_rec
    sequence
    double precision    :: total
    double precision    :: vdw_tot      ! total of dir, recip
    double precision    :: vdw_dir
    double precision    :: vdw_recip
    double precision    :: elec_tot     ! total of dir, recip, nb_adjust, self
    double precision    :: elec_dir
    double precision    :: elec_recip
    double precision    :: elec_nb_adjust
    double precision    :: elec_self
    double precision    :: hbond
    double precision    :: bond
    double precision    :: angle
    double precision    :: dihedral
    double precision    :: vdw_14
    double precision    :: elec_14
    double precision    :: restraint
  end type pme_pot_ene_rec

  integer, parameter    :: pme_pot_ene_rec_size = 16

  type(pme_pot_ene_rec), parameter      :: null_pme_pot_ene_rec = &
    pme_pot_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                    0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)

  ! Virials, with breakdown, from pme.  This is intended to be the
  ! structure for storing virials; it currently is not exported but could be.

  type pme_virial_rec
    sequence
    double precision    :: molecular(3, 3)
    double precision    :: atomic(3, 3)
    double precision    :: elec_direct(3, 3)
    double precision    :: elec_nb_adjust(3, 3)
    double precision    :: elec_recip(3, 3)
    double precision    :: elec_recip_vdw_corr(3, 3)
    double precision    :: elec_recip_self(3, 3)
    double precision    :: elec_14(3, 3)
    double precision    :: eedvir               ! used in Darden's error est.
  end type pme_virial_rec

  integer, parameter    :: pme_virial_rec_size = 73

  type(pme_virial_rec), parameter      :: null_pme_virial_rec = &
    pme_virial_rec(9*0.d0,9*0.d0,9*0.d0,9*0.d0,9*0.d0,9*0.d0,9*0.d0,9*0.d0,0.d0)

  integer, allocatable, save            :: gbl_excl_img_pairlst(:)
  integer, allocatable, save            :: gbl_nvdwcls(:)
  type(maskdata_rec), allocatable, save :: atm_nb_maskdata(:)
  integer, allocatable, save            :: atm_nb_mask(:)
  
  logical, save         :: setup_not_done = .true.

  ! The following variables don't need to be broadcast:

  integer, save         :: irespa = 0

  ! Net force and molecular virial correction factor; both used internally.

  double precision, save, private       :: frcx, frcy, frcz
  double precision, save, private       :: molvir_netfrc_corr(3, 3)
  
  ! private used in pme_force
  
  double precision, allocatable, save, private   :: nb_frc(:,:) ! used only in MPI
  double precision, allocatable, save, private   :: img_frc(:,:)

contains

!*******************************************************************************
!
! Subroutine:  alloc_pme_force_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_pme_force_mem(atm_cnt, ext_cnt, ntypes)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  integer                       :: ext_cnt
  integer                       :: ntypes

! Local variables:

  integer               :: alloc_failed

  if( allocated(gbl_excl_img_pairlst)) deallocate(gbl_excl_img_pairlst)
  if( allocated(gbl_nvdwcls))          deallocate(gbl_nvdwcls)
  if( allocated(atm_nb_maskdata)) deallocate(atm_nb_maskdata)
  if( allocated(atm_nb_mask))     deallocate(atm_nb_mask)

  allocate(gbl_excl_img_pairlst(atm_cnt + ext_cnt + ext_cnt), &
           gbl_nvdwcls(ntypes), &
           atm_nb_maskdata(atm_cnt), &
           atm_nb_mask(ext_cnt + ext_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  gbl_excl_img_pairlst(:) = 0
  gbl_nvdwcls(:) = 0
  atm_nb_maskdata(:) = maskdata_rec(0,0)
  atm_nb_mask(:) = 0
  
  if( numtasks > 1) then
      if( allocated(nb_frc))  deallocate(nb_frc)
      allocate(nb_frc(3,atm_cnt),stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
  endif
  
  if( allocated(img_frc)) deallocate(img_frc)
  allocate(img_frc(3,atm_cnt),stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error
  
  return

end subroutine alloc_pme_force_mem

!*******************************************************************************
!
! Subroutine:  bcast_pme_force_dat
!
! Description: <TBS>
!              
!        Used only with MPI
!*******************************************************************************

subroutine bcast_pme_force_dat(atm_cnt, ext_cnt, ntypes)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer       :: atm_cnt
  integer       :: ext_cnt
  integer       :: ntypes

! Local variables:

  ! Nothing to broadcast.  We just allocate storage in the non-master nodes.

  if (.not. master) then
    call alloc_pme_force_mem(atm_cnt, ext_cnt, ntypes)
  end if

  ! The allocated data is not initialized from the master node.

  return

end subroutine bcast_pme_force_dat

!*******************************************************************************
!
! Subroutine:  pme_force
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine pme_force(atm_cnt, crd, saved_crd, box, saved_box, vel, frc, mass, img_atm_map, atm_img_map, atm_owner_map, &
                     my_atm_lst, my_atm_cnt, new_list, pot_ene, virial,&
                     ekcmt, pme_err_est, atm_jrc, atm_xc, atm_weight, igroup, natc, tranvec, imin_par )

  use angles_mod
  use bonds_mod
  use constraints_mod
  use dist_constr_mod
  use dihedrals_mod
  use dynamics_dat_mod
  use dynamics_mod
  use pme_direct_mod
  use pme_recip_mod
  use loadbal_mod
  use img_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_pairlist_mod
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use prmtop_dat_mod
  use runfiles_mod
  use timers_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3, atm_cnt)
  double precision              :: saved_crd(3, atm_cnt)
  double precision              :: box(3)
  double precision              :: saved_box(3)
  double precision              :: vel(3, atm_cnt)
  double precision              :: frc(3, atm_cnt)
  double precision              :: mass(atm_cnt)
  integer                       :: img_atm_map(atm_cnt)
  integer                       :: atm_img_map(atm_cnt)

  integer                       :: atm_owner_map(atm_cnt)
  integer                       :: my_atm_lst(atm_cnt)  ! used only in MPI
  integer                       :: my_atm_cnt           ! used only in MPI

  logical                       :: new_list
  type(pme_pot_ene_rec)         :: pot_ene
  double precision, optional    :: virial(3)            ! Only used for MD
  double precision, optional    :: ekcmt(3)             ! Only used for MD
  double precision, optional    :: pme_err_est          ! Only used for MD
  
  integer                       :: atm_jrc(*)
  double precision              :: atm_xc(3,*)
  double precision              :: atm_weight(*)
  integer                       :: igroup(*)
  integer                       :: natc
  
  double precision              :: tranvec(*)
  integer                       :: imin_par

! Local variables:

  type(pme_virial_rec)          :: vir
  double precision              :: vir_vs_ene
  integer                       :: atm_lst_idx   ! used only in MPI
  integer                       :: alloc_failed
  integer                       :: i, j
  
  integer,save                  :: print_info

  call zero_time()
  call zero_pme_time()
  
!  write(mdout,*) "Hello I am in pme_force  pt 1"
!  open(11,file="ene_fort.dat",FORM='formatted',status='unknown')
!  write(11,*) "Hello I am in pme_force  pt 2"
!  close(11)
!  call flush_mdout()

! Zero energies that are stack or call parameters:
    
  pot_ene = null_pme_pot_ene_rec
  vir = null_pme_virial_rec

  virial(:) = 0.d0
  ekcmt(:) = 0.d0
  pme_err_est = 0.d0

! Zero internal energies, virials, etc.

  vir_vs_ene = 0.d0

  frcx = 0.d0
  frcy = 0.d0
  frcz = 0.d0

  molvir_netfrc_corr(:,:) = 0.d0

! If no force calcs are to be done, clear the frc array and bag out now.

  if (ntf .eq. 8) then
    frc(:,:) = 0.d0
    return
  end if
  
! Calculate the non-bonded contributions:

! Direct part of ewald plus vdw, hbond, pairlist setup and image claiming:

  if (ntp .gt. 0) call fill_tranvec(tranvec)

  ! The following encapsulation (save_imgcrds) seems to be necessary to
  ! prevent an optimization bug with the SGI f90 compiler.  Sigh...
  

  if (new_list) then 
    
    call pme_list(atm_cnt, crd, saved_crd, box, saved_box, igroup, gbl_excl_img_pairlst, atm_nb_maskdata, &
                  atm_nb_mask, atm_owner_map, my_atm_lst, my_atm_cnt, tranvec, imin_par)
                  
    if( numtasks .gt. 1) then
       call start_loadbal_timer
    endif
    call save_imgcrds(atm_cnt, used_img_lo, used_img_hi, gbl_img, &
                      gbl_saved_imgcrd) 
  else
    if( numtasks .gt. 1) then
      call start_loadbal_timer
    endif
    call adjust_imgcrds(atm_cnt, used_img_lo, used_img_hi, gbl_img, &
                        img_atm_map, gbl_saved_imgcrd, crd, &
                        saved_crd, box, saved_box, ntp)
  end if

! BEGIN DBG
!if( numtasks .gt. 1) then
! if (used_img_range_wraps) then
!   do i = used_img_hi + 1, used_img_lo - 1
!     if (img_atm_map(i) .gt. 0) write(mdout,*)'DBG: GBL_IMG RNG VIOLATION! IMG: ',i
!   end do
! else
!   do i = 1, used_img_lo - 1
!     if (img_atm_map(i) .gt. 0) write(mdout,*)'DBG: GBL_IMG RNG VIOLATION! IMG: ',i
!   end do
!   do i = used_img_hi + 1, atm_cnt
!     if (img_atm_map(i) .gt. 0) write(mdout,*)'DBG: GBL_IMG RNG VIOLATION! IMG: ',i
!   end do
! end if
!end if
! END DBG

  if( numtasks .gt. 1) then

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

  else ! (numtasks == 1)
    do i = 1, atm_cnt
        img_frc(1, i) = 0.d0
        img_frc(2, i) = 0.d0
        img_frc(3, i) = 0.d0
    end do
  end if ! (numtasks > 1)
  
! Don't do recip if PME is not invoked. Don't do it this step unless
! mod(irespa,nrepsa) = 0

  if (mod(irespa, nrespa) .eq. 0) then

    ! Self energy:

    ! The small amount of time used here gets lumped with the recip stuff...

    if (master) then
      call self(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self )
    end if

    ! Reciprocal energy:

    if(numtasks .gt. 1) then
    ! We intentionally keep the load balance counter running through the
    ! fft dist transposes in the recip code; synchronization will mess up the
    ! times a bit, but when not all nodes are doing recip calcs, it will
    ! really help with load balancing.
        if (i_do_recip) then
            call update_pme_time(pme_misc_timer)
            call update_loadbal_timer(elapsed_100usec_other)
            call do_pmesh_kspace(crd, gbl_img, img_frc, pot_ene%elec_recip, &
                                 vir%elec_recip, frcx, frcy, frcz)
            call update_loadbal_timer(elapsed_100usec_recipfrc)
            if (nrespa .gt. 1) call respa_scale(atm_cnt, img_frc, nrespa)
        end if
    else
        call update_pme_time(pme_misc_timer) 
    
        call do_pmesh_kspace(crd, gbl_img, img_frc, pot_ene%elec_recip, &
                             vir%elec_recip, frcx, frcy, frcz)
        if (nrespa .gt. 1) call respa_scale(atm_cnt, img_frc, nrespa)
    end if

! Long range dispersion contributions:

! Continuum method:
    
    if (vdwmeth .eq. 1 .and. master) then
      call vdw_correction(pot_ene%vdw_recip, vir%elec_recip_vdw_corr )
    end if

  end if      ! respa

  if( numtasks .gt. 1) then
    call update_loadbal_timer(elapsed_100usec_other)
  endif

! Direct part of ewald plus vdw, hbond, force and energy calculations:

  call update_pme_time(pme_misc_timer)
  call get_nb_energy(img_frc, gbl_img, gbl_eed_cub, gbl_ipairs, tranvec, &
                     pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                     vir%eedvir, vir%elec_direct)
  call update_pme_time(dir_frc_sum_timer)

! Adjust energies, forces for masked out pairs:

  call nb_adjust(gbl_img, crd, img_frc, img_atm_map, gbl_excl_img_pairlst, &
                 gbl_eed_cub, pot_ene%elec_nb_adjust, vir%elec_nb_adjust)
  if( numtasks .gt. 1) then
    call update_loadbal_timer(elapsed_100usec_dirfrc)
  endif
  call update_pme_time(adjust_masked_timer)

! Calculate total nonbonded energy components.

  pot_ene%vdw_tot = pot_ene%vdw_dir + &
                    pot_ene%vdw_recip
                    
  pot_ene%elec_tot = pot_ene%elec_dir + &
                     pot_ene%elec_recip + &
                     pot_ene%elec_nb_adjust + &
                     pot_ene%elec_self
  
!  write(36,*) " Hello I am in pme_force  pt 2"
                     
  write(iunit_debug,*) "pme_force: "
  write(iunit_debug,*) "pot_ene%elec_dir = ", pot_ene%elec_dir
  write(iunit_debug,*) "pot_ene%elec_recip = ", pot_ene%elec_recip
  write(iunit_debug,*) "pot_ene%elec_nb_adjust = ", pot_ene%elec_nb_adjust
  write(iunit_debug,*) "pot_ene%elec_self = ", pot_ene%elec_self
  write(iunit_debug,*) "pot_ene%elec_tot = ", pot_ene%elec_tot
  write(iunit_debug,*) "  "

  call get_atm_rel_crd(my_mol_cnt, gbl_mol_atms, gbl_mol_com, crd, &
                         atm_rel_crd, gbl_my_mol_lst)

  if( numtasks .gt. 1) then

! Clear the force array. We delay copying the nonbonded forces into the force
! array in order to be able to schedule i/o later, and batch it up.

      do atm_lst_idx = 1, my_atm_cnt
         i = my_atm_lst(atm_lst_idx)
         frc(:,i) = 0.d0
      end do

      call update_pme_time(pme_misc_timer)
      call update_time(nonbond_time)

! Calculate the other contributions:

      call pme_bonded_force(crd, frc, pot_ene, vir%molecular, vir%elec_14, &
                            atm_jrc,atm_xc,atm_weight,natc,atm_owner_map )

  ! Sum up total potential energy for this task:

      pot_ene%total = pot_ene%vdw_tot + &
                      pot_ene%elec_tot + &
                      pot_ene%hbond + &
                      pot_ene%bond + &
                      pot_ene%angle + &
                      pot_ene%dihedral + &
                      pot_ene%vdw_14 + &
                      pot_ene%elec_14 + &
                      pot_ene%restraint
                               
  ! Adjustment of total energy for constraint energies does not seem
  ! consistent, but it matches sander...

! The above stuff gets lumped as "other time"...

      call update_loadbal_timer(elapsed_100usec_other)

      if (new_list) then
#ifdef SLOW_NONBLOCKING_MPI
        call get_img_frc_distribution(atm_cnt, my_atm_cnt, gbl_atm_offsets, gbl_taskmap, &
                                      gbl_inv_taskmap, &
                                      gbl_send_atm_lst, gbl_send_atm_cnts, &
                                      gbl_recv_atm_lsts, gbl_recv_atm_cnts)
#else
        call get_img_frc_distribution(atm_cnt, my_atm_cnt, gbl_atm_offsets, gbl_taskmap, &
                                      gbl_send_atm_lst, gbl_send_atm_cnts, &
                                      gbl_recv_atm_lsts, gbl_recv_atm_cnts)
#endif
      end if

      call distribute_img_frcs(atm_cnt, img_frc, nb_frc, atm_img_map, &
                               gbl_atm_offsets, gbl_taskmap, &
                               gbl_inv_taskmap, &
                               gbl_send_atm_lst, gbl_send_atm_cnts, &
                               gbl_recv_atm_lsts, gbl_recv_atm_cnts, &
                               dbl_mpi_send_buf, dbl_mpi_recv_buf, my_atm_lst, my_atm_cnt )

      call update_time(fcve_dist_time)
      call zero_pme_time()

  ! First phase of virial work.  We need just the nonbonded forces at this
  ! stage.

      do atm_lst_idx = 1, my_atm_cnt
          i = my_atm_lst(atm_lst_idx)
          vir%molecular(:,1) = vir%molecular(:,1) + nb_frc(:,i) * atm_rel_crd(1,i)
          vir%molecular(:,2) = vir%molecular(:,2) + nb_frc(:,i) * atm_rel_crd(2,i)
          vir%molecular(:,3) = vir%molecular(:,3) + nb_frc(:,i) * atm_rel_crd(3,i)
          molvir_netfrc_corr(:, 1) = molvir_netfrc_corr(:, 1) + atm_rel_crd(1, i)
          molvir_netfrc_corr(:, 2) = molvir_netfrc_corr(:, 2) + atm_rel_crd(2, i)
          molvir_netfrc_corr(:, 3) = molvir_netfrc_corr(:, 3) + atm_rel_crd(3, i)
      end do

  ! Finish up virial work; Timing is inconsequential...

          vir%atomic(:,:) = vir%elec_recip(:,:) + &
                            vir%elec_direct(:,:) + &
                            vir%elec_nb_adjust(:,:) + &
                            vir%elec_recip_vdw_corr(:,:) + &
                            vir%elec_recip_self(:,:) + &
                            vir%elec_14(:,:)

          vir%molecular(:,:) = vir%molecular(:,:) + vir%atomic(:,:)
    
      call get_ekcom(my_mol_cnt, gbl_mol_atms, gbl_mol_mass_inv, &
                       ekcmt, vel, mass, gbl_my_mol_lst)

      call update_time(nonbond_time)
      call update_pme_time(pme_misc_timer)

! Add pot_ene, ekcmt(1:3) (md only), and the bulk of the ew_ene_vir common
! block together from all nodes.

      call distribute_enes_virs_netfrcs(pot_ene, vir, ekcmt)
    
      call update_time(fcve_dist_time)
      call zero_pme_time()

! Adjust forces for net force (comes from recip work):

      frcx = frcx / atm_cnt
      frcy = frcy / atm_cnt
      frcz = frcz / atm_cnt

  ! Copy image forces to atom forces, correcting for netfrc as you go, if
  ! appropriate.  Do not remove net force if netfrc = 0; e.g. in minimization.

      if (netfrc .gt. 0) then

          do atm_lst_idx = 1, my_atm_cnt
            i = my_atm_lst(atm_lst_idx)
            frc(1, i) = frc(1, i) + nb_frc(1, i) - frcx
            frc(2, i) = frc(2, i) + nb_frc(2, i) - frcy
            frc(3, i) = frc(3, i) + nb_frc(3, i) - frcz
          end do

    ! Correct the molecular virial for netfrc:

            vir%molecular(1,:) = vir%molecular(1,:) - molvir_netfrc_corr(1,:) * frcx
            vir%molecular(2,:) = vir%molecular(2,:) - molvir_netfrc_corr(2,:) * frcy
            vir%molecular(3,:) = vir%molecular(3,:) - molvir_netfrc_corr(3,:) * frcz
  
      else  ! ( netfrc == 0 )

          do atm_lst_idx = 1, my_atm_cnt
            i = my_atm_lst(atm_lst_idx)
            frc(1, i) = frc(1, i) + nb_frc(1, i)
            frc(2, i) = frc(2, i) + nb_frc(2, i)
            frc(3, i) = frc(3, i) + nb_frc(3, i)
          end do
          
      end if  ! ( netfrc .gt. 0 )

  else  ! ( numtasks == 1 )

! Copy nonbonded forces into force array, adjusting for net force if 
! appropriate as you go.

      frcx = frcx / atm_cnt
      frcy = frcy / atm_cnt
      frcz = frcz / atm_cnt

  ! Copy image forces to atom forces, correcting for netfrc as you go, if
  ! appropriate.  Do not remove net force if netfrc = 0; e.g. in minimization.

      if (netfrc .gt. 0) then
         do i = 1, atm_cnt
           j = atm_img_map(i)
           frc(1, i) = img_frc(1, j) - frcx
           frc(2, i) = img_frc(2, j) - frcy
           frc(3, i) = img_frc(3, j) - frcz
         end do
      else
         do i = 1, atm_cnt
           j = atm_img_map(i)
           frc(1, i) = img_frc(1, j)
           frc(2, i) = img_frc(2, j)
           frc(3, i) = img_frc(3, j)
         end do
      end if
  
  ! First phase of virial work.  We need just the nonbonded forces at this
  ! stage, though there will be a further correction in get_dihed.

        do i = 1, atm_cnt
          vir%molecular(:,1) = vir%molecular(:,1) + frc(:,i) * atm_rel_crd(1,i)
          vir%molecular(:,2) = vir%molecular(:,2) + frc(:,i) * atm_rel_crd(2,i)
          vir%molecular(:,3) = vir%molecular(:,3) + frc(:,i) * atm_rel_crd(3,i)
        end do

        call get_ekcom(my_mol_cnt, gbl_mol_atms, gbl_mol_mass_inv, &
                       ekcmt, vel, mass, gbl_my_mol_lst)

      call update_time(nonbond_time)
      call update_pme_time(pme_misc_timer)

! Calculate the other contributions:

      call pme_bonded_force(crd, frc, pot_ene, vir%molecular, vir%elec_14, &
                            atm_jrc,atm_xc,atm_weight,natc, atm_owner_map )

  ! Sum up total potential energy for this task:

      pot_ene%total = pot_ene%vdw_tot + &
                      pot_ene%elec_tot + &
                      pot_ene%hbond + &
                      pot_ene%bond + &
                      pot_ene%angle + &
                      pot_ene%dihedral + &
                      pot_ene%vdw_14 + &
                      pot_ene%elec_14 + &
                      pot_ene%restraint
                               
  ! Adjustment of total energy for constraint energies does not seem
  ! consistent, but it matches sander...

      call zero_pme_time()
  
  ! Finish up virial work; Timing is inconsequential...

        vir%atomic(:,:) = vir%elec_recip(:,:) + &
                          vir%elec_direct(:,:) + &
                          vir%elec_nb_adjust(:,:) + &
                          vir%elec_recip_vdw_corr(:,:) + &
                          vir%elec_recip_self(:,:) + &
                          vir%elec_14(:,:)

        vir%molecular(:,:) = vir%molecular(:,:) + vir%atomic(:,:)

  end if  ! ( numtasks > 1)

! Calculate vir_vs_ene in master.

  if (master) then

    vir_vs_ene = vir%elec_recip(1, 1) + &
                 vir%elec_recip(2, 2) + &
                 vir%elec_recip(3, 3) + &
                 vir%eedvir + &
                 vir%elec_nb_adjust(1, 1) + &
                 vir%elec_nb_adjust(2, 2) + &
                 vir%elec_nb_adjust(3, 3)

    ! Avoid divide-by-zero for pure neutral systems (l-j spheres).

    if (pot_ene%elec_tot .ne. 0.0d0) then
      vir_vs_ene = abs(vir_vs_ene + pot_ene%elec_tot)/abs(pot_ene%elec_tot)
    else
      vir_vs_ene = 0.0d0
    end if

  end if

! Save virials in form used in runmd:

  virial(1) = 0.5d0 * vir%molecular(1, 1)
  virial(2) = 0.5d0 * vir%molecular(2, 2)
  virial(3) = 0.5d0 * vir%molecular(3, 3)
  pme_err_est = vir_vs_ene

  if (master .and. verbose_pme .gt. 0) then
    call write_netfrc(frcx, frcy, frcz)
    call pme_verbose_print(pot_ene, vir, vir_vs_ene)
  end if
  
  setup_not_done = .false.

  call update_time(nonbond_time)
  call update_pme_time(pme_misc_timer)

  return

end subroutine pme_force

!*******************************************************************************
!
! Subroutine:  nb_adjust
!
! Description:  The part of ewald due to gaussian counterion about an atom you
!               are bonded to or otherwise for which you do not compute the
!               nonbond pair force. NECESSARY since you ARE computing this pair
!               in the reciprocal sum.
!*******************************************************************************

subroutine nb_adjust(img, crd, img_frc, img_atm_map, excl_img_pairlst, &
                     eed_cub, ene, virial)

  use img_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  type(img_rec)         :: img(natom)
  double precision      :: crd(3, natom)
  double precision      :: img_frc(3, natom)
  integer               :: img_atm_map(*)
  integer               :: excl_img_pairlst(*)
  double precision      :: eed_cub(*)
  double precision      :: ene
  double precision      :: virial(3, 3)

! Local variables:

  integer               :: atm_i, atm_j
  integer               :: img_i, img_j
  integer               :: excl_img_j_cnt
  integer               :: nxt_sublst
  integer               :: sublst_idx
  double precision      :: ene_stk
  double precision      :: ew_coeff_stk
  double precision      :: eedtbdns_stk
  double precision      :: delx, dely, delz, delr, delr2, delrinv
  double precision      :: cgi, cgi_cgj
  double precision      :: del, dx, x
  double precision      :: erfc, derfc, d0, d1
  double precision      :: df, dfx, dfy, dfz
  double precision      :: vxx, vxy, vxz, vyy, vyz, vzz
  double precision      :: x_i, y_i, z_i
  integer               :: ind

  double precision, parameter   :: half = 1.d0/2.d0
  double precision, parameter   :: third = 1.d0/3.d0

  ene_stk = 0.d0
  ew_coeff_stk = ew_coeff
  eedtbdns_stk = eedtbdns

  vxx = 0.d0
  vxy = 0.d0
  vxz = 0.d0
  vyy = 0.d0
  vyz = 0.d0
  vzz = 0.d0

  del = 1.d0 / eedtbdns

  nxt_sublst = 0

  do img_i = my_img_lo, my_img_hi

    cgi = img(img_i)%qterm

    nxt_sublst = nxt_sublst + 1
    excl_img_j_cnt = excl_img_pairlst(nxt_sublst)

    atm_i = img_atm_map(img_i)

    x_i = crd(1, atm_i)
    y_i = crd(2, atm_i)
    z_i = crd(3, atm_i)

    do sublst_idx = nxt_sublst + 1, nxt_sublst + excl_img_j_cnt

      img_j = excl_img_pairlst(sublst_idx)
      atm_j = img_atm_map(img_j)

      delx = crd(1, atm_j) - x_i
      dely = crd(2, atm_j) - y_i
      delz = crd(3, atm_j) - z_i

      delr2 = delx * delx + dely * dely + delz * delz

      ! Similar code to that in short_ene; however the only valid option
      ! here is the erfc switch, and dxdr = ewaldcof.

      delr = sqrt(delr2)
      delrinv = 1.d0 / delr

      x = ew_coeff_stk * delr
      ind = int(eedtbdns_stk * x)
      dx = x - dble(ind) * del
      ind = ishft(ind, 2)             ! 4 * ind
      
      ! Cubic spline on erfc derfc:

      erfc = eed_cub(1 + ind) + dx * (eed_cub(2 + ind) + &
             dx * (eed_cub(3 + ind) + dx * eed_cub(4 + ind) * third) * half)

      derfc = eed_cub(2 + ind) + dx * (eed_cub(3 + ind) + &
              dx * eed_cub(4 + ind) * half)

      d0 = (erfc - 1.d0) * delrinv
      d1 = (d0 - ew_coeff_stk * derfc) * delrinv * delrinv

      cgi_cgj = cgi * img(img_j)%qterm

      ene_stk = ene_stk + cgi_cgj * d0

      df = cgi_cgj * d1

      dfx = delx * df
      dfy = dely * df
      dfz = delz * df

      vxx = vxx - dfx * delx
      vxy = vxy - dfx * dely
      vxz = vxz - dfx * delz
      vyy = vyy - dfy * dely
      vyz = vyz - dfy * delz
      vzz = vzz - dfz * delz

      img_frc(1, img_j) = img_frc(1, img_j) + dfx
      img_frc(2, img_j) = img_frc(2, img_j) + dfy
      img_frc(3, img_j) = img_frc(3, img_j) + dfz
      img_frc(1, img_i) = img_frc(1, img_i) - dfx
      img_frc(2, img_i) = img_frc(2, img_i) - dfy
      img_frc(3, img_i) = img_frc(3, img_i) - dfz

    end do

    nxt_sublst = nxt_sublst + excl_img_j_cnt

  end do

  virial(1, 1) = vxx
  virial(1, 2) = vxy
  virial(2, 1) = vxy
  virial(1, 3) = vxz
  virial(3, 1) = vxz
  virial(2, 2) = vyy
  virial(2, 3) = vyz
  virial(3, 2) = vyz
  virial(3, 3) = vzz

  ene = ene_stk

  return

end subroutine nb_adjust

!*******************************************************************************
!
! Subroutine:  distribute_enes_virs_netfrcs
!
! Description: We reduce the appropriate subset of values in the ene array,
!              the ekcmt array, and the pme_ene_vir common block
!
!              Used only with MPI
!*******************************************************************************

subroutine distribute_enes_virs_netfrcs(pot_ene, vir, ekcmt)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  type(pme_pot_ene_rec) :: pot_ene
  type(pme_virial_rec)  :: vir
  double precision      :: ekcmt(3)

! Local variables:

  type pme_dat
    sequence
    type(pme_pot_ene_rec)       :: pot_ene
    type(pme_virial_rec)        :: vir
    double precision            :: ekcmt(3)
    double precision            :: molvir_netfrc_corr(3,3)
    double precision            :: frcx, frcy, frcz
  end type pme_dat

  type(pme_dat), save           :: dat_in, dat_out
  integer                       :: buf_size

  dat_in%pot_ene = pot_ene
  dat_in%vir = vir
  
  dat_in%ekcmt(:) = ekcmt(:)
  
  dat_in%molvir_netfrc_corr(:,:) = molvir_netfrc_corr(:,:)
  dat_in%frcx = frcx
  dat_in%frcy = frcy
  dat_in%frcz = frcz

  buf_size = pme_pot_ene_rec_size + pme_virial_rec_size + 3 + 9 + 3

  call mpi_allreduce(dat_in%pot_ene%total, dat_out%pot_ene%total, &
                     buf_size, mpi_double_precision, &
                     mpi_sum, lib_mpi_comm, err_code_mpi)

  pot_ene = dat_out%pot_ene
  vir = dat_out%vir
  ekcmt(:) = dat_out%ekcmt(:)
  molvir_netfrc_corr(:,:) = dat_out%molvir_netfrc_corr(:,:)
  frcx = dat_out%frcx
  frcy = dat_out%frcy
  frcz = dat_out%frcz

  return

end subroutine distribute_enes_virs_netfrcs

!*******************************************************************************
!
! Subroutine:  self
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine self(ene, ewaldcof, vol, vir)

  use gbl_constants_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  double precision      :: ene, ewaldcof, vol, vir(3, 3)

! Local variables:

  integer                       :: i
  double precision              :: ee_plasma
  double precision, save        :: factor
  double precision, save        :: sqrt_pi
  double precision, save        :: sumq
  double precision, save        :: sumq2

! Only compute sumq and sumq2 at beginning. They don't change. This code is
! only executed by the master, so we precalc anything we can...

  if (setup_not_done ) then

    factor = -0.5d0 * PI / (ewaldcof * ewaldcof)

    sqrt_pi = sqrt(PI)

    sumq = 0.d0
    sumq2 = 0.d0

    do i = 1, natom
      sumq = sumq + atm_qterm(i)
      sumq2 = sumq2 + atm_qterm(i) * atm_qterm(i)
    end do

  end if

  ee_plasma = factor * sumq * sumq / vol

  ene = - sumq2 * ewaldcof / sqrt_pi + ee_plasma

  ! The off-diagonal elements are already zero.

  vir(1,1) = -ee_plasma
  vir(2,2) = -ee_plasma
  vir(3,3) = -ee_plasma

  return

end subroutine self

!*******************************************************************************
!
! Subroutine:  vdw_correction
!
! Description:  Get analytic estimate of energy and virial corrections due to
!               dispersion interactions beyond the cutoff.
!*******************************************************************************

subroutine vdw_correction(ene, virial)

  use mdin_ctrl_dat_mod
  use pbc_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  double precision      :: ene, virial(3, 3)

! Local variables:

  integer                       :: i, j, ic, iaci
  double precision, save        :: ene_factor            ! Result of precalc.
  double precision              :: prefac, term

! Only compute ene_factor at the beginning. It doesn't change. This code is
! only executed by the master, so we precalc anything we can...

  if (setup_not_done ) then

    term = 0.d0

    ! Will later divide by volume, which is all that could change:

    prefac = 2.d0 * PI / (3.d0 * vdw_cutoff**3)

    do i = 1, ntypes
      iaci = ntypes * (i - 1)
      do j = 1, ntypes
        ic = typ_ico(iaci + j)
        if (ic .gt. 0) term = term + &
                              gbl_nvdwcls(i) * gbl_nvdwcls(j) * gbl_cn2(ic)
      end do
    end do
    
    ene_factor = -prefac * term

  end if

  ene = ene_factor / uc_volume

  ! The off-diagonal elements are already zero.

  virial(1, 1) = - 2.d0 * ene
  virial(2, 2) = - 2.d0 * ene
  virial(3, 3) = - 2.d0 * ene

  return

end subroutine vdw_correction

!*******************************************************************************
!
! Subroutine:  respa_scale
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine respa_scale(atm_cnt, img_frc, nrespa)

  use img_mod
  use prmtop_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: img_frc(3, atm_cnt)
  integer               :: nrespa

! Local variables:

  integer               :: i

  if( numtasks .gt. 1) then
    if (used_img_range_wraps) then
        do i = used_img_lo, atm_cnt
            img_frc(1, i) = nrespa * img_frc(1, i)
            img_frc(2, i) = nrespa * img_frc(2, i)
            img_frc(3, i) = nrespa * img_frc(3, i)
        end do
        do i = 1, used_img_hi
            img_frc(1, i) = nrespa * img_frc(1, i)
            img_frc(2, i) = nrespa * img_frc(2, i)
            img_frc(3, i) = nrespa * img_frc(3, i)
        end do
    else
        do i = used_img_lo, used_img_hi
            img_frc(1, i) = nrespa * img_frc(1, i)
            img_frc(2, i) = nrespa * img_frc(2, i)
            img_frc(3, i) = nrespa * img_frc(3, i)
        end do
    end if
  else
    do i = 1, atm_cnt
        img_frc(1, i) = nrespa * img_frc(1, i)
        img_frc(2, i) = nrespa * img_frc(2, i)
        img_frc(3, i) = nrespa * img_frc(3, i)
    end do
  endif

  return

end subroutine respa_scale

!*******************************************************************************
!
! Subroutine:  pme_bonded_force
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine pme_bonded_force(crd, frc, pot_ene, molvir, e14_vir, &
                            atm_jrc,atm_xc,atm_weight,natc, atm_owner_map)

  use angles_mod
  use bonds_mod
  use constraints_mod
  use dist_constr_mod
  use dihedrals_mod
  use dynamics_mod
  use dynamics_dat_mod
  use mdin_ctrl_dat_mod
  use timers_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision              :: crd(3, *)
  double precision              :: frc(3, *)
  type(pme_pot_ene_rec)         :: pot_ene
  double precision              :: molvir(3, 3)
  double precision              :: e14_vir(3, 3)
  
  integer                       :: atm_jrc(*)
  double precision              :: atm_xc(3, *)
  double precision              :: atm_weight(*)
  integer                       :: natc
  integer                       :: atm_owner_map(*)
  

  ! These energy variables are temporaries, for summing. DON'T use otherwise!

  double precision              :: bond_ene
  double precision              :: angle_ene
  double precision              :: dihedral_ene
  double precision              :: vdw_14_ene
  double precision              :: elec_14_ene
  double precision              :: nb_constr_ene

! Bond energy contribution:

  if (ntf .le. 1) then
    if (cit_nbonh .gt. 0) then
      call get_bond_energy(cit_nbonh, cit_h_bond, crd, frc, bond_ene, atm_owner_map)
      pot_ene%bond = bond_ene
    end if
  end if

  if (ntf .le. 2) then
    if (cit_nbona .gt. 0) then
      call get_bond_energy(cit_nbona, cit_a_bond, crd, frc, bond_ene, atm_owner_map)
      pot_ene%bond = pot_ene%bond + bond_ene
    end if
  end if

  call update_time(bond_time)

! Angle energy contribution:

  if (ntf .le. 3) then
    if (cit_ntheth .gt. 0) then
      call get_angle_energy(cit_ntheth, cit_angle, crd, frc, angle_ene, atm_owner_map)
      pot_ene%angle = angle_ene
    end if
  end if

  if (ntf .le. 4) then
    if (cit_ntheta .gt. 0) then
      call get_angle_energy(cit_ntheta, cit_angle(cit_ntheth+1), &
                            crd, frc, angle_ene, atm_owner_map)
      pot_ene%angle = pot_ene%angle + angle_ene
    end if
  end if

  call update_time(angle_time)

! Dihedral energy contribution:

  if (ntf .le. 5) then
    if (cit_nphih .gt. 0) then
      call get_dihed_energy(atm_owner_map, cit_nphih, cit_dihed, crd, frc, &
                            dihedral_ene, vdw_14_ene, elec_14_ene, &
                            atm_rel_crd, molvir, e14_vir)
      pot_ene%dihedral = dihedral_ene
      pot_ene%vdw_14 = vdw_14_ene
      pot_ene%elec_14 = elec_14_ene
    end if
  end if

  if (ntf .le. 6) then
    if (cit_nphia .gt. 0) then
      call get_dihed_energy(atm_owner_map, cit_nphia, cit_dihed(cit_nphih + 1), crd, frc, &
                            dihedral_ene, vdw_14_ene, elec_14_ene, &
                            atm_rel_crd, molvir, e14_vir)
      pot_ene%dihedral = pot_ene%dihedral + dihedral_ene
      pot_ene%vdw_14 = pot_ene%vdw_14 + vdw_14_ene
      pot_ene%elec_14 = pot_ene%elec_14 + elec_14_ene
    end if
  end if

  call update_time(dihedral_time)

! Calculate the position constraints energy:

  if (natc .gt. 0) then
     call get_crd_constraint_energy(natc, pot_ene%restraint, atm_jrc, &
                                    crd, frc, atm_xc, atm_weight, atm_owner_map )
  endif

! Calculate atom-atom distance constraints energy

  if( cit_num_dist_constr .gt. 0 ) then
	 
	 call get_dist_constr_energy(cit_num_dist_constr, cit_dist_constr, crd, frc, nb_constr_ene, atm_owner_map )
	 pot_ene%restraint = pot_ene%restraint + nb_constr_ene
  
  endif
!  write(iunit_debug,*) " pme_force() pot_ene%restraint = ",pot_ene%restraint

  return

end subroutine pme_bonded_force


!*******************************************************************************
!
! Subroutine:  write_netfrc
!
! Description:  Get the netfrc's back into the external axis order and print
!               them out.  We do this all in a separate subroutine just to
!               keep from cluttering up critical code.
!              
!*******************************************************************************

subroutine write_netfrc(frcx, frcy, frcz)

  use file_io_dat_mod

  implicit none

! Formal arguments:

  double precision      :: frcx, frcy, frcz

! Local variables:

  double precision      :: netfrc_out(3)

  netfrc_out(1) = frcx
  netfrc_out(2) = frcy
  netfrc_out(3) = frcz

  write(mdout, 33) netfrc_out(1), netfrc_out(2), netfrc_out(3)

  return

33     format(1x, 'NET FORCE PER ATOM: ', 3(1x, e12.4))

end subroutine write_netfrc

!*******************************************************************************
!
! Subroutine:  pme_verbose_print
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine pme_verbose_print(pot_ene, vir, vir_vs_ene)

  use file_io_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  type(pme_pot_ene_rec) :: pot_ene
  type(pme_virial_rec)  :: vir 
  double precision      :: vir_vs_ene

! Local variables:

  if (.not. master) return

  if (verbose_pme .ge. 1) then
    write(mdout, '(4(/,5x,a,f22.12))') &
          'Evdw                   = ', pot_ene%vdw_tot, &
          'Ehbond                 = ', pot_ene%hbond, &
          'Ecoulomb               = ', pot_ene%elec_tot
    write(mdout, '(2(/,5x,a,f22.12))') &
          'Iso virial             = ',  &
          vir%molecular(1, 1) + vir%molecular(2, 2) + vir%molecular(3, 3), &
          'Eevir vs. Ecoulomb     = ', vir_vs_ene
  end if

  if (verbose_pme .ge. 2) then
    write(mdout, '(4(/,5x,a,f22.12),/)') &
          'E electrostatic (self) = ', pot_ene%elec_self, &
          '                (rec)  = ', pot_ene%elec_recip, &
          '                (dir)  = ', pot_ene%elec_dir, &
          '                (adj)  = ', pot_ene%elec_nb_adjust
    write(mdout, 30) vir%molecular(1, 1), &
                     vir%molecular(1, 2), &
                     vir%molecular(1, 3)
    write(mdout, 30) vir%molecular(2, 1), &
                     vir%molecular(2, 2), &
                     vir%molecular(2, 3)
    write(mdout, 30) vir%molecular(3, 1), &
                     vir%molecular(3, 2), &
                     vir%molecular(3, 3)
30                   format(5x, 'MOLECULAR VIRIAL: ', 3(1x, e14.8))
   
  end if

  if (verbose_pme .eq. 3) then
    write(mdout, *) '--------------------------------------------'
    write(mdout, 31) vir%elec_recip(1, 1), &
                     vir%elec_recip(1, 2), &
                     vir%elec_recip(1, 3)
    write(mdout, 31) vir%elec_recip(2, 1), &
                     vir%elec_recip(2, 2), &
                     vir%elec_recip(2, 3)
    write(mdout, 31) vir%elec_recip(3, 1), &
                     vir%elec_recip(3, 2), &
                     vir%elec_recip(3, 3)
    write(mdout, *) '..................'
31     format(5x, 'Reciprocal VIRIAL: ', 3(1x, e14.8))
    write(mdout, 32) vir%elec_direct(1, 1), &
                     vir%elec_direct(1, 2), &
                     vir%elec_direct(1, 3)
    write(mdout, 32) vir%elec_direct(2, 1), &
                     vir%elec_direct(2, 2), &
                     vir%elec_direct(2, 3)
    write(mdout, 32) vir%elec_direct(3, 1), &
                     vir%elec_direct(3, 2), &
                     vir%elec_direct(3, 3)
    write(mdout, *) '..................'
32     format(5x, 'Direct VIRIAL: ', 3(1x, e14.8))
    write(mdout, 38) vir%eedvir
    write(mdout, *) '..................'
38     format(5x, 'Dir Sum EE vir trace: ', e14.8)
    write(mdout, 33) vir%elec_nb_adjust(1, 1), &
                     vir%elec_nb_adjust(1, 2), &
                     vir%elec_nb_adjust(1, 3)
    write(mdout, 33) vir%elec_nb_adjust(2, 1), &
                     vir%elec_nb_adjust(2, 2), &
                     vir%elec_nb_adjust(2, 3)
    write(mdout, 33) vir%elec_nb_adjust(3, 1), &
                     vir%elec_nb_adjust(3, 2), &
                     vir%elec_nb_adjust(3, 3)
    write(mdout, *) '..................'
33     format(5x, 'Adjust VIRIAL: ', 3(1x, e14.8))
    write(mdout, 34) vir%elec_recip_vdw_corr(1, 1), &
                     vir%elec_recip_vdw_corr(1, 2), &
                     vir%elec_recip_vdw_corr(1, 3)
    write(mdout, 34) vir%elec_recip_vdw_corr(2, 1), &
                     vir%elec_recip_vdw_corr(2, 2), &
                     vir%elec_recip_vdw_corr(2, 3)
    write(mdout, 34) vir%elec_recip_vdw_corr(3, 1), &
                     vir%elec_recip_vdw_corr(3, 2), &
                     vir%elec_recip_vdw_corr(3, 3)
    write(mdout, *) '..................'
34     format(5x, 'Recip Disp. VIRIAL: ', 3(1x, e14.8))
    write(mdout, 35) vir%elec_recip_self(1, 1), &
                     vir%elec_recip_self(1, 2), &
                     vir%elec_recip_self(1, 3)
    write(mdout, 35) vir%elec_recip_self(2, 1), &
                     vir%elec_recip_self(2, 2), &
                     vir%elec_recip_self(2, 3)
    write(mdout, 35) vir%elec_recip_self(3, 1), &
                     vir%elec_recip_self(3, 2), &
                     vir%elec_recip_self(3, 3)
    write(mdout, *) '..................'
35     format(5x, 'Self VIRIAL: ', 3(1x, e14.8))
    write(mdout, 36) vir%elec_14(1, 1), &
                     vir%elec_14(1, 2), &
                     vir%elec_14(1, 3)
    write(mdout, 36) vir%elec_14(2, 1), &
                     vir%elec_14(2, 2), &
                     vir%elec_14(2, 3)
    write(mdout, 36) vir%elec_14(3, 1), &
                     vir%elec_14(3, 2), &
                     vir%elec_14(3, 3)
    write(mdout, *) '..................'
36     format(5x, 'E14 VIRIAL: ', 3(1x, e14.8))
    write(mdout, 37) vir%atomic(1, 1), &
                     vir%atomic(1, 2), &
                     vir%atomic(1, 3)
    write(mdout, 37) vir%atomic(2, 1), &
                     vir%atomic(2, 2), &
                     vir%atomic(2, 3)
    write(mdout, 37) vir%atomic(3, 1), &
                     vir%atomic(3, 2), &
                     vir%atomic(3, 3)
37     format(5x, 'Atomic VIRIAL: ', 3(1x, e14.8))
    write(mdout, *)'--------------------------------------------'
  end if

  return

end subroutine pme_verbose_print

end module pme_force_mod

subroutine set_irespa(irespa_new)

    use pme_force_mod
    
    implicit none
    
    integer :: irespa_new
    
    irespa = irespa_new

end subroutine set_irespa 
