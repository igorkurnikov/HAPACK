#include "copyright.i"

!*******************************************************************************
!
! Module: alltasks_setup_mod
!
! Description:  Setup of data structures for uniprocessor code as well as
!               mpi master and slave processors.
!              
!*******************************************************************************

module alltasks_setup_mod

use file_io_dat_mod

  implicit none

  private       do_initial_pme_atom_division  ! used only in MPI
  private       check_atm_division            ! used only in MPI

contains

!*******************************************************************************
!
! Subroutine:  alltasks_setup
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alltasks_setup(crd,box,vel,igroup,imin_par,ntb_par,atm_owner_map,my_atm_lst,my_atm_cnt,mass_inv)

  use angles_mod
  use bonds_mod
  use cit_mod
  use constraints_mod
  use dihedrals_mod
  use dist_constr_mod
  use dynamics_mod
  use dynamics_dat_mod
  use ene_frc_splines_mod
  use pme_direct_mod
  use pme_force_mod
  use pme_recip_mod
  use img_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_pairlist_mod
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use prmtop_dat_mod
  use random_mod
  use mdin_amoeba_dat_mod
  use amoeba_recip_mod
  use amoeba_interface_mod

  implicit none
 
 ! formal variables
 
  double precision      :: crd(3,*)
  double precision      :: box(3)
  double precision      :: vel(3,*)
  integer               :: igroup(*) 
  integer               :: imin_par
  integer               :: ntb_par
  integer               :: atm_owner_map(*)
  integer               :: my_atm_lst(*)
  integer               :: my_atm_cnt 
  double precision      :: mass_inv(*)
  
! Local variables:

  integer               :: alloc_failed

  integer               :: i                   ! used only in MPI
  integer               :: node                ! used only in MPI
  integer               :: taskmap_idx         ! used only in MPI
  integer, allocatable  :: extra_used_atms(:)  ! used only in MPI
  integer,allocatable   :: use_atm_map(:)
   
  call gen_atom_maps_setup(atm_owner_map,my_atm_lst,my_atm_cnt)

!  write(iunit_debug,*)"alltasks_setup pt 1 numtasks = ",numtasks

  if(numtasks .gt. 1) then
     
    call bcast_dynamics_dat(natom, nspm) ! dynamics_mod
 
    if (.not. master .and. using_pme_potential) then
        call set_cit_tbl_dims(box, vdw_cutoff + skinnb)
    end if
    
    if (using_pme_potential) then
        call bcast_img_dat(natom)                              ! img_mod
        call bcast_nb_pairlist_dat(natom, vdw_cutoff + skinnb) ! nb_pairlist_mod
        call bcast_pme_force_dat(natom, next, ntypes)          ! pme_force_mod
        call bcast_pme_direct_dat                              ! pme_direct_mod
        
        call bcast_ene_frc_splines_dat                         ! ene_frc_splines_mod

        ! Set up data structures used in n to n mpi exchanges:

        ! The gbl_taskmap is an ordering of tasks other than this task ordered from
        ! mytaskid + 1 to mytaskid - 1, with wraparound at numtasks.

        if(allocated(gbl_taskmap)) deallocate(gbl_taskmap)
        if(allocated(gbl_inv_taskmap)) deallocate(gbl_inv_taskmap)
    
        allocate(gbl_taskmap(numtasks - 1), &
                gbl_inv_taskmap(numtasks - 1), &
                stat = alloc_failed)

        if (alloc_failed .ne. 0) call setup_alloc_error

        ! The following taskmaps set up the order for what is essentially a
        ! "synchronous shuffle exchange" (Tam and Wang) which we discovered
        ! independently.  We use it for our asynchronous comm.
        
        taskmap_idx = 0
        node = mytaskid + 1

        do
          if (node .ge. numtasks) node = 0
          if (node .eq. mytaskid) exit
          taskmap_idx = taskmap_idx + 1
          gbl_taskmap(taskmap_idx) = node
          node = node + 1
        end do

        ! The order of gbl_inv_taskmap is inverted.  If you receive and send
        ! sequentially, it works best to recv via gbl_taskmap while you
        ! simultaneously send via gbl_inv_taskmap (or vice versa).

        do taskmap_idx = 1, numtasks - 1
            gbl_inv_taskmap(taskmap_idx) = gbl_taskmap(numtasks - taskmap_idx)
        end do

        ! We must know about fft slab allocations before we do atom division:

        call pme_recip_setup(imin_par)

        ! Divide atoms up among the processors.  The atom division is redone
        ! periodically under cit, and is either residue or molecule-based, with
        ! locality.  In other words, under cit a contiguous block of atoms owned by
        ! each process is a thing of the past.

        call do_initial_pme_atom_division(natom,crd, atm_owner_map, my_atm_lst, my_atm_cnt, imin_par)
       
       ! We can now check that atom division is along residue boundaries.  We do
       ! this for constant pressure runs, where division is molecule-based. There
       ! appears to be a bug in leap whereby users can somehow create residues
       ! split between two molecules.  Ideally, we would find this problem for
       ! single and multiprocessor code, constant pressure and constant volume.
       ! Pragmatically, we just check in pmemd in the case where it will wreak
       ! havoc due to improper atom division among processors.

        if (ntb_par .eq. 2) call check_atm_division(atm_owner_map)

    else if (using_gb_potential) then

        ! For the first GB implementation, we attempt to evenly divide the
        ! atom workload on residue boundaries.  We will keep track of atom
        ! ownership, but will keep all coordinates updated in all processes,
        ! assuming that large cutoff size and small atom count make it somewhere
        ! between unnecessary and perhaps even detrimental to do otherwise
        ! (in other words, attaining processor spatial locality is not practical
        ! with a small irregularly shaped population of atoms, given the large
        ! nonbonded cutoffs typical of GB).

        call do_initial_gb_atom_division(atm_owner_map, my_atm_lst, my_atm_cnt)
        
    end if  !  if (using_pme_potential)
  end if ! if (numtasks .gt. 1)

  ! Bond-angle-dihedral ownership needs to be established.  We use the
  ! use_atm_map for both pme and GB here; in reality the atom usage info
  ! will not be kept for GB though.

  allocate(use_atm_map(natom), stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error
  use_atm_map(:) = 0
 
  call bonds_setup(use_atm_map, atm_owner_map)
  call angles_setup(use_atm_map, atm_owner_map)  
  call dihedrals_setup(use_atm_map, atm_owner_map)
  call dist_constr_setup(use_atm_map, atm_owner_map)
  
  if( numtasks .gt. 1) then

  ! Extra used atoms setup is now possible.  This is only done for pme:

      if (using_pme_potential) then

	      if(allocated(extra_used_atms)) deallocate(extra_used_atms)
          allocate(extra_used_atms(natom), stat = alloc_failed) ! temporary
          if (alloc_failed .ne. 0) call setup_alloc_error
          extra_used_atms(:) = 0

          extra_used_atm_cnt = 0

          do i = 1, natom
              if (use_atm_map(i) .ne. 0) then
                  if (atm_owner_map(i) .ne. mytaskid) then
                      extra_used_atm_cnt = extra_used_atm_cnt + 1
                      extra_used_atms(extra_used_atm_cnt) = i
                  end if
              end if
          end do

          if (extra_used_atm_cnt .gt. 0) then

	          if( allocated(gbl_extra_used_atms) ) deallocate(gbl_extra_used_atms)
              allocate(gbl_extra_used_atms(extra_used_atm_cnt), stat = alloc_failed)
              if (alloc_failed .ne. 0) call setup_alloc_error

              gbl_extra_used_atms(1:extra_used_atm_cnt) = &
              extra_used_atms(1:extra_used_atm_cnt)

          end if

          deallocate(extra_used_atms) ! Release temporary buffer.

      else if (using_gb_potential) then

          extra_used_atm_cnt = 0      ! not used...

      end if ! end: if (using_pme_potential)  

  ! Send/recv list allocs can now be done because my_atm_cnt is known. Only
  ! used for pme...

      if (using_pme_potential) then

	      if(allocated(gbl_send_atm_lst))  deallocate(gbl_send_atm_lst)
	      if(allocated(gbl_send_atm_cnts)) deallocate(gbl_send_atm_cnts)
	      if(allocated(gbl_recv_atm_lsts)) deallocate(gbl_recv_atm_lsts)
	      if(allocated(gbl_recv_atm_cnts)) deallocate(gbl_recv_atm_cnts)
	
          allocate(gbl_send_atm_lst(natom), &
                  gbl_send_atm_cnts(0 : numtasks - 1), &
                  gbl_recv_atm_lsts(my_atm_cnt, numtasks - 1), &
                  gbl_recv_atm_cnts(0 : numtasks - 1), &
                  stat = alloc_failed)

          if (alloc_failed .ne. 0) call setup_alloc_error

          gbl_send_atm_lst(:) = 0
          gbl_send_atm_cnts(:) = 0
          gbl_recv_atm_lsts(:,:) = 0
          gbl_recv_atm_cnts(:) = 0

      end if ! end: if (using_pme_potential) 

  ! Allocate buffers for mpi i/o that will remain allocated throughout the
  ! run.  This is done because some mpi implementations (myrinet in particular)
  ! mmap the buffers, and using stack buffers could have negative performance
  ! impacts in result (ifc also mmap/munmaps the dynamic stack space).
  ! These buffers are not used when other static allocations are available, say
  ! when broadcasting initialization data from the master.  All i/o using
  ! these buffers must be completed before routine exit (ie., waits MUST be
  ! done on nonblocking i/o).
  ! We must be sure that the buffer size variables have been updated before
  ! this call...

      if (using_pme_potential) then
#ifdef SLOW_NONBLOCKING_MPI
          call set_minimum_mpi_bufs_size(3 * natom)
#else
          call set_minimum_mpi_bufs_size(max(3*natom, 3*my_atm_cnt*(numtasks-1)))
#endif
      else if (using_gb_potential) then
          call set_minimum_mpi_bufs_size(3 * natom)
      end if

  else !  ( numtasks == 1) /* begin non-MPI code */

  ! We set up reciprocal force data structures here in parallel with
  ! where it has to be done for mpi code:

    if (using_pme_potential) call pme_recip_setup(imin_par)

  ! gbl_bond is still needed for shake setup and resetup, but gbl_angle and
  ! gbl_dihed can be deallocated .

  !  deallocate(gbl_angle, gbl_dihed)

  ! gbl_dist_constr can be deallocated too, I assume (IGOR)

  !  deallocate(gbl_dist_constr)
   
  endif !  ( numtasks > 1)

  if (allocated(use_atm_map)) deallocate(use_atm_map)
  
!  write(iunit_debug,*)"alltasks_setup pt end "

  return

end subroutine alltasks_setup

!*******************************************************************************
!
! Internal Subroutine:  do_initial_pme_atom_division
!
! Description:  Determine atom distribution for parallel processing.
!
!*******************************************************************************

subroutine do_initial_pme_atom_division(atm_cnt,crd, atm_owner_map, my_atm_lst, my_atm_cnt, imin_par)

  use cit_mod
  use pme_fft_mod
  use gbl_constants_mod
  use loadbal_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use pbc_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3,atm_cnt)
  integer                       :: atm_owner_map(atm_cnt)
  integer                       :: my_atm_lst(atm_cnt)
  integer                       :: my_atm_cnt
  integer                       :: imin_par

! Local variables:

  integer               :: alloc_failed
  integer               :: node
  double precision      :: fraction(3, atm_cnt) ! in range 0.0 - +0.999...
  integer               :: crd_idx_lst_tbl(0 : cit_tbl_x_dim - 1, &
                                           0 : cit_tbl_y_dim - 1, &
                                           0 : cit_tbl_z_dim - 1)
  type(atm_lst_rec)     :: atm_lst(0:atm_cnt)

  
! The following checks for minimum atoms, residues, and molecules per
! processor are intended to avoid hitting conditions we may not
! have considered.  The use of more processors than atoms, residues or
! molecules would be a bit ridiculous anyway...

! Check for sufficient atoms for the number of processors:

  if (natom .lt. 10 * numtasks) then
    error_msg = 'ERROR: Must have 10x more atoms than processors!'
    call mol_mech_error
  end if

! Check for sufficient residues for the number of processors:

  if (nres .lt. 4 * numtasks) then
    error_msg = 'ERROR: Must have 4x more residues than processors!'
    call mol_mech_error
  end if

! For CP MD, check for sufficient molecules for the number of processors:

  if (ntp .ne. 0 .and. nspm .lt. 4 * numtasks) then
    error_msg = 'ERROR: Must have 4x more molecules than processors!'
    call mol_mech_error
  end if
  
  ! Set up image division for force calcs.  Images ARE assigned in contiguous
  ! blocks, without wrap around.  For respa runs and minimizations, asymmetric
  ! fft slab load balancing (basically, trying to use as few tasks doing
  ! fft's/recip force) will not be done for various technical reasons.

  if (nrespa .eq. 1 .and. imin_par .eq. 0) then
    call divide_images_recip_biased(fft_workload_estimate, gbl_img_div_tbl )
  else
    call divide_images_evenly(gbl_img_div_tbl) ! Equal image division
  end if                                               ! res boundaries ignored.
  
  ! Master needs space to receive data for load balancing.  The 5 values
  ! for each task are the "direct force time" due to image nonbonded calcs,
  ! the "reciprocal force time" due to pme nonbonded reciprocal force calcs,
  ! the "bad force time" due to bond-angle-dihedral force calcs, the cit setup
  ! time, and the send_atm_cnts total for each node, which can be used to
  ! determine when to do redistribution of the atom workload.

  if (master) then
	if(allocated(gbl_loadbal_node_dat)) deallocate(gbl_loadbal_node_dat)
    allocate(gbl_loadbal_node_dat(5, 0:numtasks - 1), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
    gbl_loadbal_node_dat(:,:) = 0
  end if

  ! All nodes divide the work in runmd. It is a residue-based division
  ! for constant volume simulations and a molecule-based division for
  ! constant pressure simulations.  The division is map- and list-based, in
  ! order to improve locality.  The division will be redone periodically to
  ! maintain reasonable locality as atoms move.
  
  ! Make the atom to image map:
  
  call get_fract_crds(atm_cnt, crd, fraction)
  call setup_crd_idx_tbl(atm_cnt, fraction, crd_idx_lst_tbl, atm_lst)
  call divide_atoms(atm_cnt, fraction, crd_idx_lst_tbl, atm_lst, atm_owner_map, my_atm_lst, my_atm_cnt)

  return

end subroutine do_initial_pme_atom_division

!*******************************************************************************
!
! Internal Subroutine:  do_initial_gb_atom_division
!
! Description:  Determine atom distribution for parallel processing.
!
!               Used Only in MPI
!*******************************************************************************

subroutine do_initial_gb_atom_division(atm_owner_map,my_atm_lst, my_atm_cnt)

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none
  
! Formal Arguments

  integer               :: atm_owner_map(*)
  integer               :: my_atm_lst(*)
  integer               :: my_atm_cnt

! Local variables:

  integer               :: alloc_failed
  integer               :: i
  integer               :: task_id 

! The following checks for minimum atoms and residues per
! processor are intended to avoid hitting conditions we may not
! have considered.  The use of more processors than atoms would be a bit
! ridiculous; with GB though, we may actually get close to 1 processor per
! residue.  We thus allow a lower limit there of 1.01 residue per processor.

! Check for sufficient atoms for the number of processors:

  if (natom .lt. 10 * numtasks) then
	error_msg = 'ERROR: Must have 10x more atoms than processors!'
    call mol_mech_error
  end if

! Check for sufficient residues for the number of processors:

  if (dble(nres) .lt. dble(numtasks) * 1.01d0) then
    error_msg = 'ERROR: Must have 1.01x more residues than processors!'
    call mol_mech_error
  end if
	
  ! Do residue-based atom division.

  call gb_atom_division(atm_owner_map)

  ! All the other data structures can be set up using gbl_vec_rcvcnts.

  my_atm_cnt = gbl_vec_rcvcnts(mytaskid)

!BEGIN DBG
! if (master) write(6,*)'Node atom cnts =', gbl_vec_rcvcnts(:)
!END DBG

  gbl_atm_offsets(0) = 0

  do task_id = 0, numtasks - 1
    gbl_atm_offsets(task_id + 1) = gbl_atm_offsets(task_id) + &
                                   gbl_vec_rcvcnts(task_id)
  end do

  ! The atom list is not strictly necessary in a GB context, as atoms
  ! are owned in contiguous blocks.  However, the guts of pmemd is set up
  ! without this assumption, so...
  
  do i = 1, my_atm_cnt
     my_atm_lst(i) = gbl_atm_offsets(mytaskid) + i
  end do

  do task_id = 0, numtasks - 1
    gbl_vec_rcvcnts(task_id) = 3 * gbl_vec_rcvcnts(task_id)
  end do

  gbl_vec_rcvcnts(numtasks) = 0

  gbl_vec_offsets(0) = 0

  do task_id = 0, numtasks - 1
    gbl_vec_offsets(task_id + 1) = gbl_vec_offsets(task_id) + &
                                   gbl_vec_rcvcnts(task_id)
  end do

  return

contains

!*******************************************************************************
!
! Subroutine:  gb_atom_division
!
! Description:  <TBS>
!
!               Used Only in MPI
!*******************************************************************************

subroutine gb_atom_division(atm_owner_map)

  implicit none
  
  integer               :: atm_owner_map(natom)

  integer               :: task_id
  integer               :: first_atm_id, last_atm_id
  integer               :: first_res_id, last_res_id
  integer               :: res_assigned
  integer               :: task_res_cnt
  integer               :: task_res_cnts(0:numtasks - 1)  ! Number of residue per MPI process, task_id(rank) starts with 0
  double precision      :: per_task_target
  double precision      :: total_target

  per_task_target = dble(nres) / dble(numtasks)

  total_target = 0.d0
  res_assigned = 0

  do task_id = 0, numtasks - 2
    total_target = total_target + per_task_target
    task_res_cnt = int(total_target) - res_assigned
    res_assigned = res_assigned + task_res_cnt
    task_res_cnts(task_id) = task_res_cnt              ! Number of residue per MPI process, task_id (or rank) starts with 0
  end do

  task_res_cnts(numtasks - 1) = nres - res_assigned

  ! Initialize the various data structures that will reflect the atom division:

  gbl_vec_rcvcnts(:) = 0
  atm_owner_map(:) = -1

  ! Use the task_res_cnts() array to do the division...

  first_res_id = 1

  do task_id = 0, numtasks - 1
    last_res_id = first_res_id + task_res_cnts(task_id) - 1
    first_atm_id = gbl_res_atms(first_res_id)
    last_atm_id = gbl_res_atms(last_res_id + 1) - 1
    gbl_vec_rcvcnts(task_id) = last_atm_id - first_atm_id + 1   !
    atm_owner_map(first_atm_id:last_atm_id) = task_id
    first_res_id = last_res_id + 1
  end do
  
  return

end subroutine gb_atom_division

end subroutine do_initial_gb_atom_division

!*******************************************************************************
!
! Subroutine:  check_atm_division
!
! Description:  Confirm that no residues are split between two molecules.
!               This is intended for use under constant pressure, mpi.
!
!               Used Only in MPI
!*******************************************************************************

subroutine check_atm_division(atm_owner_map)

  use gbl_constants_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer       :: atm_owner_map(natom)

! Local variables:

  integer       :: atm_idx
  integer       :: first_res_atm
  integer       :: last_res_atm
  integer       :: res_idx
  integer       :: res_owner
  
  do res_idx = 1, nres

    first_res_atm = gbl_res_atms(res_idx)
    last_res_atm = gbl_res_atms(res_idx + 1) - 1

    res_owner = atm_owner_map(first_res_atm)

    do atm_idx = first_res_atm + 1, last_res_atm

      if (atm_owner_map(atm_idx) .ne. res_owner) then

        if (master) then
          write(error_msg, '(a,a,a,i6,a,i7,a,i7,a)') &
            'Bad residue/molecule data !', extra_line_hdr, &
            'Residue ', res_idx, '(atoms ', first_res_atm, '-', &
            last_res_atm, ') belong to multiple MPI nodes.'
          call mol_mech_error
        else
          error_msg = ""
          call mol_mech_error
        end if

      end if

    end do

  end do

  return

end subroutine check_atm_division

subroutine gen_atom_maps_setup(atm_owner_map,my_atm_lst,my_atm_cnt)
! Allocate atom mapping structures and set them for serial run 
!
    use parallel_dat_mod
    use file_io_dat_mod
    use prmtop_dat_mod
    use dynamics_dat_mod
    
    implicit none
    
! Formal Parameters:

    integer  :: atm_owner_map(natom)
    integer  :: my_atm_lst(natom)
    integer  :: my_atm_cnt
    
! Local variables
    
    integer :: i
    integer :: alloc_failed
    
! Allocate storage for various structures associated with keeping track of the
! division of the atom and image workload.  The gbl_img_div_tbl has one extra
! elements at the end for passing a boolean back to slave nodes indicating that
! atom workload redistribution is needed.

  if(allocated(gbl_img_div_tbl)) deallocate(gbl_img_div_tbl)   ! used only in MPI PME 
  if(allocated(gbl_vec_offsets)) deallocate(gbl_vec_offsets)
  if(allocated(gbl_atm_offsets)) deallocate(gbl_atm_offsets)
  if(allocated(gbl_vec_rcvcnts)) deallocate(gbl_vec_rcvcnts)
 
  allocate(gbl_img_div_tbl(0:numtasks + 1), &
           gbl_vec_offsets(0:numtasks), &
           gbl_atm_offsets(0:numtasks), &
           gbl_vec_rcvcnts(0:numtasks), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error    
  
  if( numtasks .eq. 1) then
    master = .true.
    mytaskid = 0
    my_atm_cnt = natom
    my_mol_cnt = nspm  
    gbl_vec_offsets(0) = 0
    gbl_vec_rcvcnts(0) = natom
    do i = 1,natom
       atm_owner_map(i) = 0
       my_atm_lst(i)    = i 
    enddo
    do i = 1,nspm
       gbl_my_mol_lst(i) = i
    enddo
  endif   
    
end subroutine gen_atom_maps_setup

end module alltasks_setup_mod

subroutine open_iunit_debug(irank)

    use parallel_dat_mod

    implicit none
    
! Formal argumnets:

    integer  :: irank
    
! Local variables

    logical :: is_opened

    iunit_debug = 50 + irank    
!    idebug_flag = 1
!   iunit_debug = 0
    
    inquire( iunit_debug, opened = is_opened )
    if( .not. is_opened) then
        open(iunit_debug,FORM='formatted',status='unknown')
    endif

end subroutine open_iunit_debug
