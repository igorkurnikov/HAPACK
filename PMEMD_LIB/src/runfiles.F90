#include "copyright.i"
!*******************************************************************************
!
! Module: runfiles_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module runfiles_mod

use file_io_dat_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  corpac
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine corpac(iend, crd, istart, nf)

  implicit none

! Formal arguments:

  integer               :: iend
  double precision      :: crd(*)
  integer               :: istart
  integer               :: nf

! Local variables:

  integer, save         :: imax = 0
  integer               :: i
  integer               :: iobuf_cnt
  double precision      :: iobuf(3, 150)

  double precision, parameter   :: rmin = -999.99d0
  double precision, parameter   :: rmax = 9999.99d0
  integer, parameter            :: iobuf_siz = size(iobuf, 2)

  ! It is ASSUMED that iend is a multiple of 3; anything else would be
  ! invalid!  Also istart is a multiple of 3 (including 0) + 1.

! NOTE - corpac input may actually be crds, vels, or box data.

  ! We only check for out-of-range values as long as they have not occurred
  ! before.  We assume that once they have occurred, it is unlikely that
  ! we will get reverse diffusion.

  ! BUGBUG - This check should be unnecessary for vels, as any system with
  !          vel velocities past f8.3 format is well past exploding.

  if (imax .eq. 0) then
    do i = istart, iend
      if (crd(i) .gt. rmax .or. crd(i) .lt. rmin) then
        imax = 1
        exit
      end if
    end do
  end if

  i = istart

  do

    iobuf_cnt = 0

    do

      if (iobuf_cnt .ge. iobuf_siz .or. i .gt. iend) exit

      iobuf_cnt = iobuf_cnt + 1

      iobuf(1, iobuf_cnt) = crd(i )
      iobuf(2, iobuf_cnt) = crd(i + 1)
      iobuf(3, iobuf_cnt) = crd(i + 2)

      i = i + 3

    end do

    if (imax .eq. 0) then
      write(nf, 1000) iobuf(:, 1:iobuf_cnt)
    else
      write(nf, 1001) iobuf(:, 1:iobuf_cnt)
    end if

    if (i .gt. iend) exit
      
  end do

 1000 format(10f8.3)
 1001 format(10f8.2)

  return

end subroutine corpac

!*******************************************************************************
!
! Internal Subroutine:  mdwrit
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine mdwrit(imin_par, ntwr_par, ntxo_par, ntb_par, nstep, atm_cnt, crd, box, vel, tt)

  use file_io_mod
  use mdin_ctrl_dat_mod

  implicit none

! Formal arguments:

  integer               :: imin_par
  integer               :: ntwr_par
  integer               :: ntxo_par
  integer               :: ntb_par
  integer               :: nstep
  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: box(3)
  double precision      :: vel(3, atm_cnt)
  double precision      :: tt

! Local variables:

  integer               :: istart, iend
  character(12)         :: num
  character(89)         :: restrt2_name

! Write/rewind the restrt:

  call write_restart(restrt, imin_par, ntxo_par, ntb_par, atm_cnt, crd, box, vel, tt)

  rewind(restrt)

! Consider whether to save 2ndary restrt:

  if (ntwr_par .ge. 0) return

  do iend = 1, 80
    if (restrt_name(iend:iend) .le. ' ') exit
  end do

  iend = iend - 1

  write(num, '(i12)') nstep

  do istart = 1, 12
    if (num(istart:istart) .ne. ' ') exit
  end do

  write(restrt2_name, '(a,a,a)') restrt_name(1:iend), '_', num(istart:12)

  write(mdout, '(a,a)') ' writing ', restrt2_name

  if (ntxo_par .eq. 0) then
    call amopen(restrt2, restrt2_name, 'U', 'U', 'W')
  else
    call amopen(restrt2, restrt2_name, 'U', 'F', 'W')
  end if

  call write_restart(restrt2, imin_par, ntxo_par, ntb_par, atm_cnt, crd, box, vel, tt)

  close(restrt2)

  return

end subroutine mdwrit

!*******************************************************************************
!
! Subroutine:  mdeng
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine mdeng(nstep, time, si, box)

  use mdin_ctrl_dat_mod
  use pbc_mod
  use state_info_mod

  implicit none

! Formal arguments:

  integer               :: nstep
  double precision      :: time
  double precision      :: si(*)        ! State information.
  double precision      :: box(3)

! Local variables:

  integer               :: i
  logical, save         :: first = .true.
  character(16), save   :: labs(41)

  data labs/'Nsteps  ', 'time(ps)  ', 'Etot  ', 'EKinetic  ', &
            'Temp  ', 'T_solute ', 'T_solv  ', 'Pres_scal_solu ', &
            'Pres_scal_solv ', 'BoxX  ', 'BoxY  ', 'BoxZ  ', &
            'volume  ', 'pres_X  ', 'pres_Y  ', 'pres_Z  ', &
            'Pressure ', 'EKCoM_x ', 'EKCoM_y ', 'EKCoM_z', &
            'EKComTot ', 'VIRIAL_x ', 'VIRIAL_y ', 'VIRIAL_z ', &
            'VIRIAL_tot ', 'E_pot  ', 'E_vdw  ', 'E_el  ', &
            'E_hbon  ', 'E_bon  ', 'E_angle  ', 'E_dih  ', &
            'E_14vdw  ', 'E_14el  ', 'E_const  ', 'E_pol  ', &
            'AV_permMoment ', 'AV_indMoment ', 'AV_totMoment ', &
            'Density', 'dV/dlambda'/

! Define various terms:

  if (first) then
    ! up to Ekinetic:
    write(mden, 1) 'L0 ', (labs(i), i = 1, 4)
    ! up to Pres_scal_solu:
    write(mden, 1) 'L1 ', (labs(i), i = 5, 8)
    ! up to boxZ:
    write(mden, 1) 'L2 ', (labs(i), i = 9, 12)
    ! up to pres_Z:
    write(mden, 1) 'L3 ', (labs(i), i = 13, 16)
    ! up to EKCoM_z:
    write(mden, 1) 'L4 ', (labs(i), i = 17, 20)
    ! up to VIRIAL_z:
    write(mden, 1) 'L5 ', (labs(i), i = 21, 24)
    ! up to E_el:
    write(mden, 1) 'L6 ', (labs(i), i = 25, 28)
    ! up to E_dih:
    write(mden, 1) 'L7 ', (labs(i), i = 29, 32)
    ! up to E_pol:
    write(mden, 1) 'L8 ', (labs(i), i = 33, 36)
    ! up to Density or dV/dlambda:
    write(mden, 1) 'L9 ', (labs(i), i = 37, 41)
1   format(a, 10(1x, a))
    first = .false.
  end if

! Write values for this step:

  ! Pres_scal_solu and Pres_scal_solv are not supported anymore; output values
  ! are fixed at 1.d0

  ! up to Ekinetic:
  write(mden, 2) 'L0 ', nstep, time, si(si_tot_ene), si(si_kin_ene)

  ! up to Pres_scal_solu:
  write(mden, 3) 'L1 ', si(si_temp), &
                        si(si_temp_solute), &
                        si(si_temp_solvent), 1.d0

  ! up to boxZ:
  write(mden, 3) 'L2 ', 1.d0, box(1), box(2), box(3)

  ! up to pres_Z:
  write(mden, 3) 'L3 ', si(si_volume), si(si_press_0),si(si_press_1),si(si_press_2)

  ! up to EKCoM_z:
  write(mden, 3) 'L4 ', si(si_tot_press), si(si_ekcmt_0),si(si_ekcmt_1),si(si_ekcmt_2)

  ! up to VIRIAL_z:
  write(mden, 3) 'L5 ', si(si_tot_ekcmt), si(si_vir_0),si(si_vir_1),si(si_vir_2)

  ! up to E_el:
  write(mden, 3) 'L6 ', si(si_tot_virial), &
                        si(si_pot_ene), si(si_vdw_ene), si(si_elect_ene)

  ! up to E_dih:
  write(mden, 3) 'L7 ', si(si_hbond_ene), &
                        si(si_bond_ene), si(si_angle_ene), si(si_dihedral_ene)

  ! up to E_pol :
  write(mden, 3) 'L8 ', si(si_vdw_14_ene), si(si_elect_14_ene), &
                        si(si_restraint_ene), si(si_polar)

  ! up to dV/dlambda, includes 3 0.d0 fields:
  write(mden, 3) 'L9 ', 0.d0, 0.d0, 0.d0, si(si_density), si(si_dv_dlambda)

2 format(a, i8, 20(2x, e17.10))
3 format(a, 20(e17.10, 2x))

  return

end subroutine mdeng

!*******************************************************************************
!
! Subroutine:  write_restart
!
! Description: Routine to write final coordinates and velocities.
!
! EWALD: dump ewald specific box information to the restrt files
!        (may only be necessary with ntx=7).
!*******************************************************************************

subroutine write_restart(nf, imin_par, ntxo_par, ntb_par, atm_cnt, crd, box, vel, tt)

  use mdin_ctrl_dat_mod
  use pbc_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: ntxo_par
  integer               :: imin_par
  integer               :: ntb_par
  integer               :: i
  integer               :: nf
  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: box(3)
  double precision      :: vel(3, atm_cnt)
  double precision      :: tt           ! Not used if minimization
  logical               :: odd_atm_cnt

  odd_atm_cnt = mod(atm_cnt, 2) .ne. 0

  if (ntxo_par .eq. 0) then                 ! Formatted writing:

    write(nf, 9008) prmtop_ititl  
    if (imin_par .eq. 0) then
      if (atm_cnt .lt. 100000) then
        write(nf, 9018) atm_cnt, tt
      else
        write(nf, 9019) atm_cnt, tt
      end if
      do i = 1, atm_cnt - 1, 2
        write(nf, 9028) crd(1, i), crd(2, i), crd(3, i), &
                        crd(1, i+1), crd(2, i+1), crd(3, i+1)
      end do
      if (odd_atm_cnt) then
        write(nf, 9028) crd(1, atm_cnt), &
                        crd(2, atm_cnt), &
                        crd(3, atm_cnt)
      end if
      do i = 1, atm_cnt - 1, 2
        write(nf, 9028) vel(1, i), vel(2, i), vel(3, i), &
                        vel(1, i+1), vel(2, i+1), vel(3, i+1)
      end do
      if (odd_atm_cnt) then
        write(nf, 9028) vel(1, atm_cnt), &
                        vel(2, atm_cnt), &
                        vel(3, atm_cnt)
      end if
    else
      if (atm_cnt .lt. 100000) then
        write(nf, 9018) atm_cnt
      else
        write(nf, 9019) atm_cnt
      end if
      do i = 1, atm_cnt - 1, 2
        write(nf, 9028) crd(1, i), crd(2, i), crd(3, i), &
                        crd(1, i+1), crd(2, i+1), crd(3, i+1)
      end do
      if (odd_atm_cnt) then
        write(nf, 9028) crd(1, atm_cnt), &
                        crd(2, atm_cnt), &
                        crd(3, atm_cnt)
      end if
    end if
    if (ntb_par .ne. 0) write(nf, 9028) box(1),box(2),box(3), &
                                    pbc_alpha,pbc_beta,pbc_gamma

  else                                  ! Binary writing:

    write(nf) prmtop_ititl
    if (imin_par .eq. 0) then
      write(nf) atm_cnt, tt
      write(nf) (crd(1, i), crd(2, i), crd(3, i), i = 1, atm_cnt)
      write(nf) (vel(1, i), vel(2, i), vel(3, i), i = 1, atm_cnt)
    else
      write(nf) atm_cnt, 0.d0
      write(nf) (crd(1, i), crd(2, i), crd(3, i), i = 1, atm_cnt)
    end if

    ! Sander does not provide box/angle info here, which strikes me as a bug.

    if (ntb_par .ne. 0) write(nf) box(1), box(2), box(3), &
                              pbc_alpha, pbc_beta, pbc_gamma
  end if

 9008 format(a80)
 9018 format(i5, 5e15.7)
 9019 format(i6, 5e15.7)
 9028 format(6f12.7)

  return

end subroutine write_restart

subroutine flush_mdout()

   use mdin_ctrl_dat_mod
     
   close(mdout)
   open(unit=mdout, file=mdout_name, status='OLD', position='APPEND')
  
end subroutine 

subroutine close_mdout()

   use mdin_ctrl_dat_mod
     
   close(mdout)
  
end subroutine 

end module runfiles_mod
