#include "copyright.i"

!*******************************************************************************
!
! Module: master_setup_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module master_setup_mod

use file_io_dat_mod

  implicit none 

contains

!*******************************************************************************
!
! Subroutine:  master_setup
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine master_setup_2(crd,box,belly_atm_cnt,igroup,mass)

  use cit_mod
  use constraints_mod
  use dynamics_mod
  use dynamics_dat_mod
#ifdef DIRFRC_EFS
  use ene_frc_splines_mod
#endif /* DIRFRC_EFS */
  use pme_direct_mod
  use pme_force_mod
  use file_io_mod
  use gbl_constants_mod
  use img_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_pairlist_mod
  use nmr_calls_mod
  use pbc_mod
  use prmtop_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision      :: crd(3,*)
  double precision      :: box(3)
  integer               :: belly_atm_cnt
  integer               :: igroup(*)
  double precision      :: mass(*)

! Local variables:

  integer               :: i
  integer               :: ifind
  integer               :: inerr
  character(80)         :: mdin_title

! Here stuff that needed to be called only once after the model is set

  inerr = 0
   
 ! Check if ifbox variable from prmtop file matches actual angles. 

  if (ifbox .eq. 1) then
     if (abs(pbc_alpha - 90.0d0 ) .gt. 1.d-5 .or. &
          abs(pbc_beta - 90.0d0) .gt. 1.d-5 .or. &
          abs(pbc_gamma - 90.0d0) .gt. 1.d-5) then
          ifbox = 3
          write(*,'(a)') '     Setting ifbox to 3 for non-orthogonal unit cell'
     end if
  end if

  if (ifbox .eq. 2) then
    if (abs(pbc_alpha - 109.4712190d0) .gt. 1.d-5 .or. &
        abs(pbc_beta - 109.4712190d0) .gt. 1.d-5 .or. &
        abs(pbc_gamma - 109.4712190d0) .gt. 1.d-5) then
      write(mdout,'(/2x,a)') &
        'Error: ifbox=2 in prmtop but angles are not correct'
      inerr = 1
    end if
  end if
  
  if (using_pme_potential) then

! Set up coordinate index table dimensions; 
! for partitioning of atoms into regions probably ??
! these dimensions will not change on small changes of box dimenstions
! so no need to recompute for small changes of periodical box 

    call set_cit_tbl_dims(box, vdw_cutoff + skinnb) 

  end if

! If the user has requested NMR restraints, do a cursory read of the
! restraints file(s) now to determine the amount of memory necessary
! for these restraints, and allocate the memory in the master.

  if (nmropt .ne. 0) call init_nmr_dat()

! Init dynamics data; 

  call init_dynamics_dat(natom, nres, nspm, gbl_res_atms, atm_nsp, &
                         mass, crd)

  if (using_pme_potential) then

    ! Set up image dynamic memory:

    call alloc_img_mem(natom)

    ! Set up pairlist memory:

    call alloc_nb_pairlist_mem(natom, vdw_cutoff + skinnb)
    
    ! Set up ewald variables and memory:

    call alloc_pme_force_mem(natom, next, ntypes)
    call init_pme_direct_dat()
#ifdef DIRFRC_EFS
    call init_ene_frc_splines_dat()
#endif /* DIRFRC_EFS */

  end if

! Consistency checking:

  if (ntb .ne. 0 .and. ntp .ne. 0 .and. ifbox .eq. 0) then
    write(*, '(a,a)') error_hdr, &
      'the combination ntb != 0, ntp != 0, ifbox == 0 is not supported!'
    inerr = 1
  end if

  if (using_pme_potential) then

    if (vdw_cutoff .ge. box(1) * 0.5d0 .or. &
        vdw_cutoff .ge. box(2) * 0.5d0 .or. &
        vdw_cutoff .ge. box(3) * 0.5d0) then
      write(*, '(a,a)') error_hdr, &
                            'max cut must be < half smallest box dimension!'
      write(*, '(a,a)') extra_line_hdr, 'max cut=', vdw_cutoff
      write(*, '(a,a)') extra_line_hdr, 'box(1)=', box(1)
      write(*, '(a,a)') extra_line_hdr, 'box(2)=', box(2)
      write(*, '(a,a)') extra_line_hdr, 'box(3)=', box(3)
      inerr = 1
    end if

  end if

! Warnings:

  if (using_pme_potential .and. ibelly .gt. 0) then
     write(*, '(a,/,a,/)') 'Warning: Although EWALD will work with belly', &
           '(for equilibration), it is not strictly correct!'
  end if
  
  if (inerr .eq. 1) then
    write(error_msg, '(a,/)') ' Input errors occurred. Terminating execution.' 
    call mol_mech_error
  end if

! Load the constrained (or belly) atoms. these are read as groups:

  belly_atm_cnt = 0

! All the bond, angle, and dihedral parameters may be changed here as the
! bond, angle, and dihedral arrays are repacked! Note in particular that
! diheda_idx may also be changed.  We also count atoms in the "belly" here,
! which is probably redundant (also done in rgroup() above).

  if (ibelly .gt. 0) then

    call remove_nonbelly_bnd_ang_dihed( igroup )
    call count_belly_atoms(natom, belly_atm_cnt, igroup)

    if (.not. using_pme_potential) then

      ! The only allowable belly here has just the first belly_atm_cnt atoms
      ! in the moving part.  Confirm this.

      do i = belly_atm_cnt + 1, natom
        if (igroup(i) .ne. 0) then
          write(error_msg, '(a,a,a)')'When ibelly != 0 and igb != 0, the moving part must', &
                                ' be at the start of the molecule, which seems to', &
                                '  not be the case!'
          call mol_mech_error
        end if
      end do

    end if
  end if

! If we are reading NMR restraints/weight changes, read them, and then determine
! how many of the torsional parameters are improper:

  if (nmropt .ne. 0) then
    call nmr_read(crd, mass, mdin, mdout)
    call set_num_improp_dihed(nphih, gbl_dihed, &
                              nphia, gbl_dihed(diheda_idx), nptra)
  endif

  return

end subroutine master_setup_2

end module master_setup_mod


!*******************************************************************************
!
! Subroutine:   open_pmemd_output_files
!
! Description:  Routine to open the dumping and restart files.
!*******************************************************************************

subroutine start_mdout_log(box)
	
  use file_io_mod
  use prmtop_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use cit_mod
  use pbc_mod

  implicit none
  
  ! Formal arguments
  
  double precision      :: box(3)

  character(8)          :: date
  character(10)         :: time
  integer               :: itmp(3)
  double precision      :: rtmp(3)
  
  call amopen(mdout, mdout_name, 'U', 'F', 'W')
  write(mdout, 1000)
  write(mdout, '(a, /)') '| PMEMD implementation of SANDER, Release 9'
  
  call date_and_time(DATE=date, TIME=time) 

  write(mdout,'(12(a),/)') '| Run on ', date(5:6), '/', date(7:8), '/',  &
        date(1:4), ' at ', time(1:2), ':', time(3:4), ':', time(5:6)

  if (using_pme_potential) then

    itmp(1) = cit_tbl_x_dim
    itmp(2) = cit_tbl_y_dim
    itmp(3) = cit_tbl_z_dim

    write(mdout, '(a, 3i5)')'| Coordinate Index Table dimensions: ', &
                            itmp(1), itmp(2), itmp(3)

    rtmp(1) = box(1)/cit_tbl_x_dim
    rtmp(2) = box(2)/cit_tbl_y_dim
    rtmp(3) = box(3)/cit_tbl_z_dim

    write(mdout,'(a, 3f10.4, /)')'| Direct force subcell size = ', &
                                 rtmp(1), rtmp(2), rtmp(3)

  end if
  
  return
  
1000 format(/10x, 55('-'), /10x, &
            'Amber 9  SANDER                              2006', &
             /10x, 55('-')/)

end subroutine start_mdout_log


subroutine open_pmemd_output_files

  use bintraj_mod
  use file_io_mod
  use mdin_ctrl_dat_mod
  use prmtop_dat_mod
  use parallel_dat_mod

  implicit none

  character(10), parameter      :: file_version = '9.00'
  integer                       :: box_flag
  character(100)                :: bin4_title

! Echo the file assignments to the user:

  refc_name = "refc"

  write(mdout,*) "  "
  write(mdout, '("File Assignments:", /, 6("|", a7, ": ", a, /))' ) & 
                    'MDOUT',  mdout_name(1:70),  &
                    'RESTRT', restrt_name(1:70), 'REFC',   refc_name(1:70),   &
                    'MDVEL',  mdvel_name(1:70),  'MDEN',   mden_name(1:70),   &
                     'MDCRD',  mdcrd_name(1:70)

! ioutfm .ne. 0 selects binary output, theoretically for all files below. In
! reality though, we never open mden for binary output.

  if (ioutfm .le. 0) then               ! Formatted dumping:

    if (ntwx .gt. 0) then
      call amopen(mdcrd, mdcrd_name, 'U', 'F', 'W')
      write(mdcrd, 1000) prmtop_ititl
    end if

    if (ntwv .gt. 0) then
      call amopen(mdvel, mdvel_name, 'U', 'F', 'W')
      write(mdvel, 1000) prmtop_ititl
    end if

  else if (ioutfm .eq. 1) then

    call open_binary_files

  else if (ioutfm .eq. 2) then  ! The new "bin4" efficiency format...

    if (ntwx .gt. 0) then

      bin4_title = trim(mdcrd_name) // '.bin4'
      call amopen(mdcrd, bin4_title, 'U', 'U', 'W') 
      write(mdcrd) file_version
      write(mdcrd) prmtop_ititl

      if (ntb .gt. 0) then
        box_flag = 1
      else
        box_flag = 0
      end if

      if (ntwprt .ne. 0) then
        write(mdcrd) ntwprt, box_flag
      else
        write(mdcrd) natom, box_flag
      end if

    end if 

    if (ntwv .gt. 0) then

      bin4_title = trim(mdvel_name) // '.bin4'
      call amopen(mdvel, bin4_title, 'U', 'U', 'W')   
      write(mdvel) file_version
      write(mdvel) prmtop_ititl

      box_flag = 0

      if (ntwprt .ne. 0) then
        write(mdvel) ntwprt, box_flag
      else
        write(mdvel) natom, box_flag
      end if

    end if

  end if

! Open the energies file:

  if (ntwe .gt. 0) then
    call amopen(mden, mden_name, 'U', 'F', 'W')
  end if

! Open the restart file:

  if (ntxo .le. 0) then
    call amopen(restrt, restrt_name, 'U', 'U', 'W') 
  else
    call amopen(restrt, restrt_name, 'U', 'F', 'W') 
  end if

  return

1000 format(a80)

end subroutine open_pmemd_output_files

