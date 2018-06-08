#include "copyright.i"

!*******************************************************************************
!
! Module: timers_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module timers_mod

use file_io_dat_mod

  implicit none

! The main routine calls second() or walltime() to set these:

  double precision, save        :: run_start_cputime
  double precision, save        :: run_setup_end_cputime
  double precision, save        :: run_end_cputime
  integer, save                 :: run_start_walltime
  integer, save                 :: run_setup_end_walltime
  integer, save                 :: run_end_walltime

! Generic timer constants; use with the time_stats array.

! NOTE: fcve is frc-crds-vel-ene

  integer, parameter    :: fcve_dist_time               =  1
  integer, parameter    :: nonbond_time                 =  2
  integer, parameter    :: bond_time                    =  3
  integer, parameter    :: angle_time                   =  4
  integer, parameter    :: dihedral_time                =  5
  integer, parameter    :: shake_time                   =  6
  integer, parameter    :: runmd_time                   =  7
  integer, parameter    :: other_time                   =  8
  integer, parameter    :: nonsetup_time                =  9
  integer, parameter    :: max_generic_timer            =  9

! Fields dropped from usage (in the old timsts array):

! time_stats(7) - was Caldis, Pol. der
! time_stats(8) - was Calrate, Corf.1, Efielde
! time_stats(9) - was Dspev, Corf.2, Indip
! time_stats(10) - was Matmul, Drates
! time_stats(11) - was Kmat
! time_stats(12) - was Force, whatever that means...
! time_stats(13) - was Dinten
! time_stats(14) - was Remarc
! time_stats(15) - was RingCurr
! time_stats(16) - was Electro.
! time_stats(17) - was Anisotr.
! time_stats(18) - was ShiftDer
! time_stats(19) - was Noecalc1 - tbd
! time_stats(20) - was Noecalc2
! time_stats(21) - was GBrad
! time_stats(22) - was GBraddist

! PME timer constants; use these in update_pme_time().

  integer, parameter    :: cit_setup_timer         = max_generic_timer + 1
  integer, parameter    :: build_list_timer        = max_generic_timer + 2
  integer, parameter    :: bspline_timer           = max_generic_timer + 3
  integer, parameter    :: grid_charges_timer      = max_generic_timer + 4
  integer, parameter    :: scalar_sum_timer        = max_generic_timer + 5
  integer, parameter    :: grad_sum_timer          = max_generic_timer + 6
  integer, parameter    :: fft_timer               = max_generic_timer + 7
  integer, parameter    :: dir_frc_sum_timer       = max_generic_timer + 8
  integer, parameter    :: adjust_masked_timer     = max_generic_timer + 9
  integer, parameter    :: pme_misc_timer          = max_generic_timer + 10
  integer, parameter    :: atm_reassign_timer      = max_generic_timer + 11
  integer, parameter    :: img_reassign_timer      = max_generic_timer + 12
  integer, parameter    :: fft_slab_reassign_timer = max_generic_timer + 13
  integer, parameter    :: max_pme_timer           = max_generic_timer + 13

! Amoeba-specific timer constants; use these in update_pme_time().  For
! Amoeba, we use the pme timers plus these specific Amoeba timers where
! appropriate.

  ! BUGBUG - We are "lumping" rather than "splitting" in the early
  !          stages of development, which means we have some reclassification
  !          work to do later, potentially.  In this phase of dev, though,
  !          amoeba time categorization will look much like that of standard
  !          pme ff's, with anything that does not fit thrown into
  !          pme_misc_timer.

  integer, parameter    :: max_amoeba_timer        = max_pme_timer + 0

! GB timer constants; use these in update_gb_time().

  integer, parameter    :: calc_gb_rad_timer       = max_generic_timer + 1
  integer, parameter    :: calc_gb_diag_timer      = max_generic_timer + 2
  integer, parameter    :: calc_gb_offdiag_timer   = max_generic_timer + 3
  integer, parameter    :: dist_gb_rad_timer       = max_generic_timer + 4
  integer, parameter    :: max_gb_timer            = max_generic_timer + 4

! Maximum timer field must be kept in sync...

  integer, parameter    :: max_timer               = max_pme_timer

! Data for generic and pme timing routines:

  double precision, save        :: time_stats(max_timer) = 0.d0

  double precision, private, save     :: generic_time1, generic_time2
  double precision, private, save     :: pme_time1, pme_time2
  double precision, private, save     :: gb_time1, gb_time2

! The TIME_TEST timing facility is for high resolution performance test
! timing, intended to be used temporarily in development.  There is provision
! for labeling the timepoints and including i/o metrics (or whatever other
! numbers you want to record).  Change time_test_max_cnt as required, depending
! on how many calls you need to monitor.

  logical, save                 :: test_timers_enabled = .false.
  integer, parameter            :: time_test_max_cnt = 50
  integer, save                 :: call_cnt(time_test_max_cnt)
  double precision, save        :: io_bytes(time_test_max_cnt)
  double precision, save        :: start_cpu(time_test_max_cnt)
  double precision, save        :: elapsed_cpu(time_test_max_cnt)
  integer, save                 :: start_wall_sec(time_test_max_cnt)
  integer, save                 :: elapsed_wall_sec(time_test_max_cnt)
  integer, save                 :: start_wall_usec(time_test_max_cnt)
  integer, save                 :: elapsed_wall_usec(time_test_max_cnt)
  character(len=50), save       :: tt_id_str(time_test_max_cnt)

contains

!*******************************************************************************
!
! Subroutine:  second
!
! Description:  Get CPU time in seconds, as a real (returns double precision).
!               This is an approximation of CPU time, probably only accurate
!               in the msec range, based on the f90 cpu_time intrinsic. Timings
!               should be assumed to be relative instead of absolute.
!*******************************************************************************

subroutine second(cpu_sec)

  implicit none

! Formal arguments:

  double precision, intent(out) :: cpu_sec

! Local variables:

  real          :: cpu_time_result

  call cpu_time(cpu_time_result)
  cpu_sec = dble(cpu_time_result)

  return

end subroutine second

!*******************************************************************************
!
! Subroutine:  wall
!
! Description:  Get local time in seconds since 2000 (for relative wallclock
!               timings).
!
!*******************************************************************************

subroutine wall(wall_sec)

  implicit none

! Formal arguments:

  integer, intent(out)  :: wall_sec

! Local variables:

  character(10) :: char_date
  character(10) :: char_time
  character(5)  :: char_zone

  integer       :: dt(8)
  integer       :: year, month, day, hour, minute, sec, msec
  integer       :: feb_days
  integer       :: i

  call date_and_time(char_date, char_time, char_zone, dt)

  year = dt(1)
  month = dt(2)
  day = dt(3)
  hour = dt(5)
  minute = dt(6)
  sec = dt(7)
  msec = dt(8)

! years to days:

  wall_sec = 0

! Count the number of days to the current year. You have to determine which
! years are leap years.  The proper algorithm is:
! (year % 4 == 0 && year % 100 != 0) || year % 400 == 0.  This simplifies to
! year % 4 for the 2000-2099. Then we need to fix the code ;-)

  do i = 2000, year - 1
    if (modulo(i,4) .ne. 0) then
      wall_sec = wall_sec + 365
    else
      wall_sec = wall_sec + 366
    end if
  end do

! Is the current year a leap year:

  if (modulo(year,4) .ne. 0) then
    feb_days = 28
  else
    feb_days = 29
  endif

! Months to days

  if (month .eq. 2) then
    wall_sec = wall_sec + 31
  else if (month .eq. 3) then
    wall_sec = wall_sec + 31 + feb_days
  else if (month .eq. 4) then
    wall_sec = wall_sec + 31 + feb_days + 31
  else if (month .eq. 5) then
    wall_sec = wall_sec + 31 + feb_days + 31 + 30
  else if (month .eq. 6) then
    wall_sec = wall_sec + 31 + feb_days + 31 + 30 + 31
  else if (month .eq. 7) then
    wall_sec = wall_sec + 31 + feb_days + 31 + 30 + 31 + 30
  else if (month .eq. 8) then
    wall_sec = wall_sec + 31 + feb_days + 31 + 30 + 31 + 30 + 31
  else if (month .eq. 9) then
    wall_sec = wall_sec + 31 + feb_days + 31 + 30 + 31 + 30 + 31 + 31
  else if (month .eq. 10) then
    wall_sec = wall_sec + 31 + feb_days + 31 + 30 + 31 + 30 + 31 + 31 + 30
  else if (month .eq. 11) then
    wall_sec = wall_sec + 31 + feb_days + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31
  else if (month .eq. 12) then
    wall_sec = wall_sec + 31 + feb_days + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30
  endif

  wall_sec = wall_sec + day

! Days to hours:

  wall_sec = wall_sec * 24 + hour

! Hours to minutes:

  wall_sec = wall_sec * 60 + minute

! Minutes to seconds:

  wall_sec = wall_sec * 60 + sec

! Round on msec:

  if (msec .ge. 500) wall_sec = wall_sec + 1

  return

end subroutine wall

!*******************************************************************************
!
! Subroutine:  profile_cpu
!
! Description:  Produce a nonsetup cpu usage profile for all processes.
!
! Author: George Seibel
!*******************************************************************************

subroutine profile_cpu(imin_par, igb, iamoeba)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: imin_par
  integer               :: igb
  integer               :: iamoeba

! Local variables:

  integer               :: i
  double precision      :: percent(max_timer) ! % of total time for each 
                                              ! generic time stat.
  double precision      :: time_sum           ! sum of counted routines.
  double precision      :: percent_factor

  double precision      :: avg_time(max_timer) ! used only in MPI

  time_stats(nonsetup_time) = run_end_cputime - run_setup_end_cputime

  time_sum = 0.d0
  do i = 1, runmd_time
    time_sum = time_sum + time_stats(i)
  end do
  time_stats(other_time) = time_stats(nonsetup_time) - time_sum

  if( numtasks .gt. 1) then
    call profile_parallel_cpu(avg_time, igb, iamoeba)

    if (.not. master) return

    time_stats(1:max_timer) = avg_time(1:max_timer)
  endif

  percent_factor = 100.d0 / time_stats(nonsetup_time)

  do i = 1, max_timer
    percent(i) = time_stats(i) * percent_factor
  end do

! The | characters are for test/dacdif

  if( numtasks .gt. 1) then
     write(mdout, 10) &
      '|  NonSetup CPU Time in Major Routines, Average for All Tasks:'
  else
     write(mdout, 10) &
       '|  NonSetup CPU Time in Major Routines:'
  endif
  write(mdout, 10) '|'
  write(mdout, 10) '|     Routine           Sec        %'
  write(mdout, 10) '|     ------------------------------'
  if( numtasks .gt. 1) then
    write(mdout, 20) 'DataDistrib', &
                   time_stats(fcve_dist_time), percent(fcve_dist_time)
  endif
  write(mdout, 20) 'Nonbond    ', &
                   time_stats(nonbond_time), percent(nonbond_time)
  write(mdout, 20) 'Bond       ', &
                   time_stats(bond_time), percent(bond_time)
  write(mdout, 20) 'Angle      ', &
                   time_stats(angle_time), percent(angle_time)
  write(mdout, 20) 'Dihedral   ', &
                   time_stats(dihedral_time), percent(dihedral_time)
  write(mdout, 20) 'Shake      ', &
                   time_stats(shake_time), percent(shake_time)

  if (imin_par .eq. 0) then
    write(mdout, 20) 'RunMD      ', time_stats(runmd_time), percent(runmd_time)
  end if
  write(mdout, 20) 'Other      ', time_stats(other_time), percent(other_time)
  write(mdout, 10) '|     ------------------------------'
  write(mdout, 20) 'Total      ', time_stats(nonsetup_time)

  if (igb .eq. 0) then
     if (iamoeba .eq. 0) then
         call profile_pme_time(avg_time, percent)
     else
         call profile_amoeba_time(avg_time, percent)
     endif
  else if (igb .eq. 1 .or. igb .eq. 2 .or. igb .eq. 5 .or. igb .eq. 7 .or. igb .eq. 20) then
     call profile_gb_time(avg_time, percent)
  end if

  return

10 format(a)
15 format(/, a)
20 format('|', 5x, a, f11.2, f8.2)
30 format('|', 5x, a, f11.2, f8.2, a)

end subroutine profile_cpu

!*******************************************************************************
!
! Subroutine:  profile_pme_time
!
! Description:  PME time profiling.
!              
!*******************************************************************************

subroutine profile_pme_time(avg_time, percent)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision      :: avg_time(max_timer)
  double precision      :: percent(max_timer)

  if (.not. master) return

  if( numtasks .gt. 1) then
    write(mdout, 15) &
            '|  PME Nonbond Pairlist CPU Time, Average for All Tasks:'
  else
    write(mdout, 15) &
            '|  PME Nonbond Pairlist CPU Time:'
  endif
  write(mdout, 10) '|'
  write(mdout, 10) '|     Routine              Sec        %'
  write(mdout, 10) '|     ---------------------------------'
  write(mdout, 20) 'Set Up Cit    ', &
                   time_stats(cit_setup_timer), percent(cit_setup_timer)
  write(mdout, 20) 'Build List    ', &
                   time_stats(build_list_timer), percent(build_list_timer)
  write(mdout, 10) '|     ---------------------------------'
  write(mdout, 20) 'Total         ', time_stats(cit_setup_timer) + &
                                     time_stats(build_list_timer), &
                                     percent(cit_setup_timer) + &
                                     percent(build_list_timer)

  if( numtasks .gt. 1) then
    write(mdout, 15) &
            '|  PME Direct Force CPU Time, Average for All Tasks:'
  else
    write(mdout, 15) &
            '|  PME Direct Force CPU Time:'
  endif
  write(mdout, 10) '|'
  write(mdout, 10) '|     Routine              Sec        %'
  write(mdout, 10) '|     ---------------------------------'
  write(mdout, 20) 'NonBonded Calc', &
                   time_stats(dir_frc_sum_timer), percent(dir_frc_sum_timer)
  write(mdout, 20) 'Exclude Masked', &
                   time_stats(adjust_masked_timer), percent(adjust_masked_timer)
  write(mdout, 20) 'Other         ', &
                   time_stats(pme_misc_timer), percent(pme_misc_timer)
  write(mdout, 10) '|     ---------------------------------'
  write(mdout, 20) 'Total         ', time_stats(dir_frc_sum_timer) + &
                                     time_stats(adjust_masked_timer) + &
                                     time_stats(pme_misc_timer), &
                                     percent(dir_frc_sum_timer) + &
                                     percent(adjust_masked_timer) + &
                                     percent(pme_misc_timer)

  if( numtasks .gt. 1) then
    write(mdout, 15) &
            '|  PME Reciprocal Force CPU Time, Average for All Tasks:'
  else
    write(mdout, 15) &
            '|  PME Reciprocal Force CPU Time:'
  endif
  write(mdout, 10) '|'
  write(mdout, 10) '|     Routine              Sec        %'
  write(mdout, 10) '|     ---------------------------------'
  write(mdout, 20) '1D bspline    ', &
                   time_stats(bspline_timer), percent(bspline_timer)
  write(mdout, 20) 'Grid Charges  ', &
                   time_stats(grid_charges_timer), percent(grid_charges_timer)
  write(mdout, 20) 'Scalar Sum    ', &
                   time_stats(scalar_sum_timer), percent(scalar_sum_timer)
  write(mdout, 20) 'Gradient Sum  ', &
                   time_stats(grad_sum_timer), percent(grad_sum_timer)
  write(mdout, 20) 'FFT           ', &
                   time_stats(fft_timer), percent(fft_timer)
  write(mdout, 10) '|     ---------------------------------'
  write(mdout, 20) 'Total         ', time_stats(bspline_timer) + &
                                     time_stats(grid_charges_timer) + &
                                     time_stats(scalar_sum_timer) + &
                                     time_stats(grad_sum_timer) + &
                                     time_stats(fft_timer), &
                                     percent(bspline_timer) + &
                                     percent(grid_charges_timer) + &
                                     percent(scalar_sum_timer) + &
                                     percent(grad_sum_timer) + &
                                     percent(fft_timer)

  if( numtasks .gt. 1) then
    write(mdout, 15) &
            '|  PME Load Balancing CPU Time, Average for All Tasks:'
    write(mdout, 10) '|'
    write(mdout, 10) '|     Routine                 Sec        %'
    write(mdout, 10) '|     ------------------------------------'
    write(mdout, 20) 'Atom Reassign    ', &
                       time_stats(atm_reassign_timer), percent(atm_reassign_timer)
    write(mdout, 20) 'Image Reassign   ', &
                       time_stats(img_reassign_timer), percent(img_reassign_timer)
    write(mdout, 20) 'FFT Slab Reassign', &
                       time_stats(fft_slab_reassign_timer), &
                       percent(fft_slab_reassign_timer)
    write(mdout, 10) '|     ------------------------------------'
    write(mdout, 20) 'Total            ', time_stats(atm_reassign_timer) + &
                                          time_stats(img_reassign_timer) + &
                                          time_stats(fft_slab_reassign_timer), &
                                          percent(atm_reassign_timer) + &
                                          percent(img_reassign_timer) + &
                                          percent(fft_slab_reassign_timer)
  endif

10 format(a)
15 format(/, a)
20 format('|', 5x, a, f11.2, f8.2)
30 format('|', 5x, a, f11.2, f8.2, a)

  return

end subroutine profile_pme_time

!*******************************************************************************
!
! Subroutine:  profile_amoeba_time
!
! Description:  Amoeba time profiling.
!              
!*******************************************************************************

subroutine profile_amoeba_time(avg_time, percent)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision      :: avg_time(max_timer)
  double precision      :: percent(max_timer)

  if (.not. master) return

  if( numtasks .gt. 1) then 
      write(mdout, 15) &
             '|  Amoeba Nonbond Pairlist CPU Time, Average for All Tasks:'
  else
      write(mdout, 15) &
             '|  Amoeba Nonbond Pairlist CPU Time:'
  endif
  write(mdout, 10) '|'
  write(mdout, 10) '|     Routine              Sec        %'
  write(mdout, 10) '|     ---------------------------------'
  write(mdout, 20) 'Set Up Cit    ', &
                   time_stats(cit_setup_timer), percent(cit_setup_timer)
  write(mdout, 20) 'Build List    ', &
                   time_stats(build_list_timer), percent(build_list_timer)
  write(mdout, 10) '|     ---------------------------------'
  write(mdout, 20) 'Total         ', time_stats(cit_setup_timer) + &
                                     time_stats(build_list_timer), &
                                     percent(cit_setup_timer) + &
                                     percent(build_list_timer)

  if(numtasks .gt. 1) then
       write(mdout, 15) &
            '|  Amoeba Direct Force CPU Time, Average for All Tasks:'
  else
       write(mdout, 15) &
            '|  Amoeba Direct Force CPU Time:'
  endif
  write(mdout, 10) '|'
  write(mdout, 10) '|     Routine              Sec        %'
  write(mdout, 10) '|     ---------------------------------'
  write(mdout, 20) 'NonBonded Calc', &
                   time_stats(dir_frc_sum_timer), percent(dir_frc_sum_timer)
  write(mdout, 20) 'Exclude Masked', &
                   time_stats(adjust_masked_timer), percent(adjust_masked_timer)
  write(mdout, 20) 'Other         ', &
                   time_stats(pme_misc_timer), percent(pme_misc_timer)
  write(mdout, 10) '|     ---------------------------------'
  write(mdout, 20) 'Total         ', time_stats(dir_frc_sum_timer) + &
                                     time_stats(adjust_masked_timer) + &
                                     time_stats(pme_misc_timer), &
                                     percent(dir_frc_sum_timer) + &
                                     percent(adjust_masked_timer) + &
                                     percent(pme_misc_timer)

  if( numtasks .gt. 1) then
      write(mdout, 15) &
            '|  Amoeba Reciprocal Force CPU Time, Average for All Tasks:'
  else
      write(mdout, 15) &
            '|  Amoeba Reciprocal Force CPU Time:'
  endif
  write(mdout, 10) '|'
  write(mdout, 10) '|     Routine              Sec        %'
  write(mdout, 10) '|     ---------------------------------'
  write(mdout, 20) '1D bspline    ', &
                   time_stats(bspline_timer), percent(bspline_timer)
  write(mdout, 20) 'Grid Charges  ', &
                   time_stats(grid_charges_timer), percent(grid_charges_timer)
  write(mdout, 20) 'Scalar Sum    ', &
                   time_stats(scalar_sum_timer), percent(scalar_sum_timer)
  write(mdout, 20) 'Gradient Sum  ', &
                   time_stats(grad_sum_timer), percent(grad_sum_timer)
  write(mdout, 20) 'FFT           ', &
                   time_stats(fft_timer), percent(fft_timer)
  write(mdout, 10) '|     ---------------------------------'
  write(mdout, 20) 'Total         ', time_stats(bspline_timer) + &
                                     time_stats(grid_charges_timer) + &
                                     time_stats(scalar_sum_timer) + &
                                     time_stats(grad_sum_timer) + &
                                     time_stats(fft_timer), &
                                     percent(bspline_timer) + &
                                     percent(grid_charges_timer) + &
                                     percent(scalar_sum_timer) + &
                                     percent(grad_sum_timer) + &
                                     percent(fft_timer)

  if(numtasks .gt. 1) then
     write(mdout, 15) &
           '|  Amoeba Load Balancing CPU Time, Average for All Tasks:'
     write(mdout, 10) '|'
     write(mdout, 10) '|     Routine                 Sec        %'
     write(mdout, 10) '|     ------------------------------------'
     write(mdout, 20) 'Atom Reassign    ', &
                       time_stats(atm_reassign_timer), percent(atm_reassign_timer)
     write(mdout, 20) 'Image Reassign   ', &
                       time_stats(img_reassign_timer), percent(img_reassign_timer)
     write(mdout, 20) 'FFT Slab Reassign', &
                       time_stats(fft_slab_reassign_timer), &
                       percent(fft_slab_reassign_timer)
     write(mdout, 10) '|     ------------------------------------'
     write(mdout, 20) 'Total            ', time_stats(atm_reassign_timer) + &
                                           time_stats(img_reassign_timer) + &
                                           time_stats(fft_slab_reassign_timer), &
                                           percent(atm_reassign_timer) + &
                                           percent(img_reassign_timer) + &
                                           percent(fft_slab_reassign_timer)
  endif

10 format(a)
15 format(/, a)
20 format('|', 5x, a, f11.2, f8.2)
30 format('|', 5x, a, f11.2, f8.2, a)

  return

end subroutine profile_amoeba_time

!*******************************************************************************
!
! Subroutine:  profile_gb_time
!
! Description:  GB time profiling.
!              
!*******************************************************************************

subroutine profile_gb_time(avg_time, percent)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision      :: avg_time(max_timer)
  double precision      :: percent(max_timer)

  if (.not. master) return

  if( numtasks .gt. 1) then
    write(mdout, 15) &
            '|  Generalized Born CPU Time, Average for All Tasks:'
  else
    write(mdout, 15) &
            '|  Generalized Born CPU Time:'
  endif
  write(mdout, 10) '|'
  write(mdout, 10) '|     Routine                 Sec        %'
  write(mdout, 10) '|     ------------------------------------'
  write(mdout, 20) 'Radii Calc       ', &
                   time_stats(calc_gb_rad_timer), percent(calc_gb_rad_timer)
  write(mdout, 20) 'Diagonal Calc    ', &
                   time_stats(calc_gb_diag_timer), percent(calc_gb_diag_timer)
  write(mdout, 20) 'Off Diagonal Calc', &
                   time_stats(calc_gb_offdiag_timer), &
                   percent(calc_gb_offdiag_timer)
  if( numtasks .gt. 1 ) then
    write(mdout, 20) 'Radii Distrib    ', &
                       time_stats(dist_gb_rad_timer), &
                       percent(dist_gb_rad_timer)
    write(mdout, 10) '|     ---------------------------------'
    write(mdout, 20) 'Total            ', time_stats(calc_gb_rad_timer) + &
                                          time_stats(calc_gb_diag_timer) + &
                                          time_stats(calc_gb_offdiag_timer) + &
                                          time_stats(dist_gb_rad_timer), &
                                          percent(calc_gb_rad_timer) + &
                                          percent(calc_gb_diag_timer) + &
                                          percent(calc_gb_offdiag_timer) + &
                                          percent(dist_gb_rad_timer)
  else
    write(mdout, 10) '|     ---------------------------------'
    write(mdout, 20) 'Total            ', time_stats(calc_gb_rad_timer) + &
                                          time_stats(calc_gb_diag_timer) + &
                                          time_stats(calc_gb_offdiag_timer), &
                                          percent(calc_gb_rad_timer) + &
                                          percent(calc_gb_diag_timer) + &
                                          percent(calc_gb_offdiag_timer)
  endif

10 format(a)
15 format(/, a)
20 format('|', 5x, a, f11.2, f8.2)

  return

end subroutine profile_gb_time

!*******************************************************************************
!
! Subroutine:  profile_parallel_cpu
!
! Description:
!
!*******************************************************************************

subroutine profile_parallel_cpu(avg_time, igb, iamoeba)

  use parallel_dat_mod

  implicit none

! Formal Arguments:

  double precision              :: avg_time(max_timer)
  integer                       :: igb
  integer                       :: iamoeba

! Local variables:

  character(30)                 :: fmt1
  integer                       :: i, j
  integer                       :: stat_mpi(mpi_status_size)
  double precision              :: max_time(max_generic_timer)
  double precision              :: min_time(max_generic_timer)
  double precision              :: std_dev(max_generic_timer)
  double precision, save        :: time_stats_mpi_buf(max_timer)

50 format(i4, f8.1, f10.1, 6f7.1, f10.1)
51 format(a4, f8.1, f10.1, 6f7.1, f10.1)

  if (master) then

    fmt1 = '(a4, a8, a10, 6a7, a10)'

    do i = 1, max_generic_timer
      avg_time(i) = time_stats(i)
      std_dev(i) = time_stats(i) * time_stats(i)
      min_time(i) = time_stats(i)
      max_time(i) = time_stats(i)
    end do

    ! We leave RunMD time in place, even for minimizations, where 
    ! it will always be 0.d0

    write(iunit_debug, '(/, a, /)') &
      'Major Routine Parallel Profiling - NonSetup CPU Seconds:'
    write(iunit_debug, fmt1)' ', 'D', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '
    write(iunit_debug, fmt1)' ', 'a', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '
    write(iunit_debug, fmt1)' ', 't', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '
    write(iunit_debug, fmt1)' ', 'a', ' ', ' ', ' ', 'D', ' ', ' ', ' ', ' '
    write(iunit_debug, fmt1)' ', 'D', 'N', ' ', ' ', 'i', ' ', ' ', ' ', ' '
    write(iunit_debug, fmt1)' ', 'i', 'o', ' ', ' ', 'h', ' ', ' ', ' ', ' '
    write(iunit_debug, fmt1)' ', 's', 'n', ' ', 'A', 'e', 'S', 'R', 'O', 'T'
    write(iunit_debug, fmt1)'T', 't', 'B', 'B', 'n', 'd', 'h', 'u', 't', 'o'
    write(iunit_debug, fmt1)'a', 'r', 'o', 'o', 'g', 'r', 'a', 'n', 'h', 't'
    write(iunit_debug, fmt1)'s', 'i', 'n', 'n', 'l', 'a', 'k', 'M', 'e', 'a'
    write(iunit_debug, fmt1)'k', 'b', 'd', 'd', 'e', 'l', 'e', 'D', 'r', 'l'
    write(iunit_debug, '(78a)') ('-', i = 1, 78)
    write(iunit_debug, 50) mytaskid, &
                       time_stats(fcve_dist_time), &
                       time_stats(nonbond_time), &
                       time_stats(bond_time), &
                       time_stats(angle_time), &
                       time_stats(dihedral_time), &
                       time_stats(shake_time), &
                       time_stats(runmd_time), &
                       time_stats(other_time), &
                       time_stats(nonsetup_time)

    do i = 1, numtasks - 1

! Receive time stat from each node in turn:

      call mpi_recv(time_stats_mpi_buf, max_generic_timer, &
                    mpi_double_precision, i, i, lib_mpi_comm, stat_mpi, &
                    err_code_mpi)

! Compute total time, write times to iunit_debug:

      do j = 1, max_generic_timer
        avg_time(j) = avg_time(j) + time_stats_mpi_buf(j)
        std_dev(j) = std_dev(j) + time_stats_mpi_buf(j) * time_stats_mpi_buf(j)
        if (time_stats_mpi_buf(j) .lt. min_time(j)) &
          min_time(j) = time_stats_mpi_buf(j)
        if (time_stats_mpi_buf(j) .gt. max_time(j)) &
          max_time(j) = time_stats_mpi_buf(j)
      end do

      write(iunit_debug, 50) i, &
        time_stats_mpi_buf(fcve_dist_time), &
        time_stats_mpi_buf(nonbond_time), &
        time_stats_mpi_buf(bond_time), &
        time_stats_mpi_buf(angle_time), &
        time_stats_mpi_buf(dihedral_time), &
        time_stats_mpi_buf(shake_time), &
        time_stats_mpi_buf(runmd_time), &
        time_stats_mpi_buf(other_time), &
        time_stats_mpi_buf(nonsetup_time)

    end do

    write(iunit_debug, '(78a)') ('-', i = 1, 78)

! Print statistics:

    do i = 1, max_generic_timer
      avg_time(i) = avg_time(i)/numtasks
      std_dev(i) = sqrt(abs(std_dev(i)/numtasks - avg_time(i) * avg_time(i)))
    end do

    write(iunit_debug, 51) 'avg', avg_time(1:max_generic_timer)
    write(iunit_debug, 51) 'min', min_time(1:max_generic_timer)
    write(iunit_debug, 51) 'max', max_time(1:max_generic_timer)
    write(iunit_debug, 51) 'std', std_dev(1:max_generic_timer)

    write(iunit_debug, '(78a)') ('-', i = 1, 78)

  else

! Send time stats to master node:

    call mpi_send(time_stats, max_generic_timer, mpi_double_precision, 0, &
                  mytaskid, lib_mpi_comm, err_code_mpi)
  end if

  if (igb .eq. 0) then

60 format(i4, f8.1, 2f10.1, f7.1, f10.1)
61 format(a4, f8.1, 2f10.1, f7.1, f10.1)

    if (iamoeba .eq. 0) then    ! "Standard" PME
      if (master) then
    
        fmt1 = '(a4, a8, 2a10, a7, a10)'
    
        do i = max_generic_timer + 1, max_pme_timer
          avg_time(i) = time_stats(i)
        end do
    
        write(iunit_debug, '(/, a, /)') &
          'PME NonBond Parallel Profiling - NonSetup CPU Seconds:'
        write(iunit_debug, fmt1)' ', ' ', 'D', ' ', ' ', ' '
        write(iunit_debug, fmt1)' ', 'P', 'i', 'R', ' ', ' '
        write(iunit_debug, fmt1)' ', 'a', 'r', 'e', 'L', ' '
        write(iunit_debug, fmt1)' ', 'i', 'e', 'c', 'o', ' '
        write(iunit_debug, fmt1)' ', 'r', 'c', 'i', 'a', 'T'
        write(iunit_debug, fmt1)'T', 'l', 't', 'p', 'd', 'o'
        write(iunit_debug, fmt1)'a', 'i', 'F', 'F', 'B', 't'
        write(iunit_debug, fmt1)'s', 's', 'r', 'r', 'a', 'a'
        write(iunit_debug, fmt1)'k', 't', 'c', 'c', 'l', 'l'
        write(iunit_debug, '(49a)') ('-', i = 1, 49)
        write(iunit_debug, 60) mytaskid, &
                           time_stats(cit_setup_timer) + &
                           time_stats(build_list_timer), &
                           time_stats(dir_frc_sum_timer) + &
                           time_stats(adjust_masked_timer) + &
                           time_stats(pme_misc_timer), &
                           time_stats(bspline_timer) + &
                           time_stats(grid_charges_timer) + &
                           time_stats(scalar_sum_timer) + &
                           time_stats(grad_sum_timer) + &
                           time_stats(fft_timer), &
                           time_stats(atm_reassign_timer) + &
                           time_stats(img_reassign_timer) + &
                           time_stats(fft_slab_reassign_timer), &
                           time_stats(cit_setup_timer) + &
                           time_stats(build_list_timer) + &
                           time_stats(dir_frc_sum_timer) + &
                           time_stats(adjust_masked_timer) + &
                           time_stats(pme_misc_timer) + &
                           time_stats(bspline_timer) + &
                           time_stats(grid_charges_timer) + &
                           time_stats(scalar_sum_timer) + &
                           time_stats(grad_sum_timer) + &
                           time_stats(fft_timer) + &
                           time_stats(atm_reassign_timer) + &
                           time_stats(img_reassign_timer) + &
                           time_stats(fft_slab_reassign_timer)
    
        do i = 1, numtasks - 1
    
    ! Receive time stat from each node in turn:
    
          call mpi_recv(time_stats_mpi_buf(max_generic_timer+1), &
                        max_pme_timer - max_generic_timer, &
                        mpi_double_precision, i, i, mpi_comm_world, stat_mpi, &
                        err_code_mpi)
    
    ! Compute total time, write times to logfile:
    
          do j = max_generic_timer + 1, max_pme_timer
            avg_time(j) = avg_time(j) + time_stats_mpi_buf(j)
          end do
    
        write(iunit_debug, 60) i, &
                           time_stats_mpi_buf(cit_setup_timer) + &
                           time_stats_mpi_buf(build_list_timer), &
                           time_stats_mpi_buf(dir_frc_sum_timer) + &
                           time_stats_mpi_buf(adjust_masked_timer) + &
                           time_stats_mpi_buf(pme_misc_timer), &
                           time_stats_mpi_buf(bspline_timer) + &
                           time_stats_mpi_buf(grid_charges_timer) + &
                           time_stats_mpi_buf(scalar_sum_timer) + &
                           time_stats_mpi_buf(grad_sum_timer) + &
                           time_stats_mpi_buf(fft_timer), &
                           time_stats_mpi_buf(atm_reassign_timer) + &
                           time_stats_mpi_buf(img_reassign_timer) + &
                           time_stats_mpi_buf(fft_slab_reassign_timer), &
                           time_stats_mpi_buf(cit_setup_timer) + &
                           time_stats_mpi_buf(build_list_timer) + &
                           time_stats_mpi_buf(dir_frc_sum_timer) + &
                           time_stats_mpi_buf(adjust_masked_timer) + &
                           time_stats_mpi_buf(pme_misc_timer) + &
                           time_stats_mpi_buf(bspline_timer) + &
                           time_stats_mpi_buf(grid_charges_timer) + &
                           time_stats_mpi_buf(scalar_sum_timer) + &
                           time_stats_mpi_buf(grad_sum_timer) + &
                           time_stats_mpi_buf(fft_timer) + &
                           time_stats_mpi_buf(atm_reassign_timer) + &
                           time_stats_mpi_buf(img_reassign_timer) + &
                           time_stats_mpi_buf(fft_slab_reassign_timer)
    
        end do
    
        write(iunit_debug, '(49a)') ('-', i = 1, 49)
    
    ! Print statistics:
    
        do i = max_generic_timer + 1, max_pme_timer
          avg_time(i) = avg_time(i)/numtasks
        end do
    
        write(iunit_debug, 61) 'avg', &
                           avg_time(cit_setup_timer) + &
                           avg_time(build_list_timer), &
                           avg_time(dir_frc_sum_timer) + &
                           avg_time(adjust_masked_timer) + &
                           avg_time(pme_misc_timer), &
                           avg_time(bspline_timer) + &
                           avg_time(grid_charges_timer) + &
                           avg_time(scalar_sum_timer) + &
                           avg_time(grad_sum_timer) + &
                           avg_time(fft_timer), &
                           avg_time(atm_reassign_timer) + &
                           avg_time(img_reassign_timer) + &
                           avg_time(fft_slab_reassign_timer), &
                           avg_time(cit_setup_timer) + &
                           avg_time(build_list_timer) + &
                           avg_time(dir_frc_sum_timer) + &
                           avg_time(adjust_masked_timer) + &
                           avg_time(pme_misc_timer) + &
                           avg_time(bspline_timer) + &
                           avg_time(grid_charges_timer) + &
                           avg_time(scalar_sum_timer) + &
                           avg_time(grad_sum_timer) + &
                           avg_time(fft_timer) + &
                           avg_time(atm_reassign_timer) + &
                           avg_time(img_reassign_timer) + &
                           avg_time(fft_slab_reassign_timer)
    
        write(iunit_debug, '(49a)') ('-', i = 1, 49)
    
      else
    
    ! Send time stats to master node:
    
        call mpi_send(time_stats(max_generic_timer+1), &
                      max_pme_timer - max_generic_timer, mpi_double_precision, &
                      0, mytaskid, mpi_comm_world, err_code_mpi)
      end if

    else        ! Amoeba PME:

      if (master) then
    
        fmt1 = '(a4, a8, 2a10, a7, a10)'
    
        do i = max_generic_timer + 1, max_amoeba_timer
          avg_time(i) = time_stats(i)
        end do
    
        write(iunit_debug, '(/, a, /)') &
          'Amoeba NonBond Parallel Profiling - NonSetup CPU Seconds:'
        write(iunit_debug, fmt1)' ', ' ', 'D', ' ', ' ', ' '
        write(iunit_debug, fmt1)' ', 'P', 'i', 'R', ' ', ' '
        write(iunit_debug, fmt1)' ', 'a', 'r', 'e', 'L', ' '
        write(iunit_debug, fmt1)' ', 'i', 'e', 'c', 'o', ' '
        write(iunit_debug, fmt1)' ', 'r', 'c', 'i', 'a', 'T'
        write(iunit_debug, fmt1)'T', 'l', 't', 'p', 'd', 'o'
        write(iunit_debug, fmt1)'a', 'i', 'F', 'F', 'B', 't'
        write(iunit_debug, fmt1)'s', 's', 'r', 'r', 'a', 'a'
        write(iunit_debug, fmt1)'k', 't', 'c', 'c', 'l', 'l'
        write(iunit_debug, '(49a)') ('-', i = 1, 49)
        write(iunit_debug, 60) mytaskid, &
                           time_stats(cit_setup_timer) + &
                           time_stats(build_list_timer), &
                           time_stats(dir_frc_sum_timer) + &
                           time_stats(adjust_masked_timer) + &
                           time_stats(pme_misc_timer), &
                           time_stats(bspline_timer) + &
                           time_stats(grid_charges_timer) + &
                           time_stats(scalar_sum_timer) + &
                           time_stats(grad_sum_timer) + &
                           time_stats(fft_timer), &
                           time_stats(atm_reassign_timer) + &
                           time_stats(img_reassign_timer) + &
                           time_stats(fft_slab_reassign_timer), &
                           time_stats(cit_setup_timer) + &
                           time_stats(build_list_timer) + &
                           time_stats(dir_frc_sum_timer) + &
                           time_stats(adjust_masked_timer) + &
                           time_stats(pme_misc_timer) + &
                           time_stats(bspline_timer) + &
                           time_stats(grid_charges_timer) + &
                           time_stats(scalar_sum_timer) + &
                           time_stats(grad_sum_timer) + &
                           time_stats(fft_timer) + &
                           time_stats(atm_reassign_timer) + &
                           time_stats(img_reassign_timer) + &
                           time_stats(fft_slab_reassign_timer)
    
        do i = 1, numtasks - 1
    
    ! Receive time stat from each node in turn:
    
          call mpi_recv(time_stats_mpi_buf(max_generic_timer+1), &
                        max_amoeba_timer - max_generic_timer, &
                        mpi_double_precision, i, i, mpi_comm_world, stat_mpi, &
                        err_code_mpi)
    
    ! Compute total time, write times to logfile:
    
          do j = max_generic_timer + 1, max_amoeba_timer
            avg_time(j) = avg_time(j) + time_stats_mpi_buf(j)
          end do
    
        write(iunit_debug, 60) i, &
                           time_stats_mpi_buf(cit_setup_timer) + &
                           time_stats_mpi_buf(build_list_timer), &
                           time_stats_mpi_buf(dir_frc_sum_timer) + &
                           time_stats_mpi_buf(adjust_masked_timer) + &
                           time_stats_mpi_buf(pme_misc_timer), &
                           time_stats_mpi_buf(bspline_timer) + &
                           time_stats_mpi_buf(grid_charges_timer) + &
                           time_stats_mpi_buf(scalar_sum_timer) + &
                           time_stats_mpi_buf(grad_sum_timer) + &
                           time_stats_mpi_buf(fft_timer), &
                           time_stats_mpi_buf(atm_reassign_timer) + &
                           time_stats_mpi_buf(img_reassign_timer) + &
                           time_stats_mpi_buf(fft_slab_reassign_timer), &
                           time_stats_mpi_buf(cit_setup_timer) + &
                           time_stats_mpi_buf(build_list_timer) + &
                           time_stats_mpi_buf(dir_frc_sum_timer) + &
                           time_stats_mpi_buf(adjust_masked_timer) + &
                           time_stats_mpi_buf(pme_misc_timer) + &
                           time_stats_mpi_buf(bspline_timer) + &
                           time_stats_mpi_buf(grid_charges_timer) + &
                           time_stats_mpi_buf(scalar_sum_timer) + &
                           time_stats_mpi_buf(grad_sum_timer) + &
                           time_stats_mpi_buf(fft_timer) + &
                           time_stats_mpi_buf(atm_reassign_timer) + &
                           time_stats_mpi_buf(img_reassign_timer) + &
                           time_stats_mpi_buf(fft_slab_reassign_timer)
    
        end do
    
        write(iunit_debug, '(49a)') ('-', i = 1, 49)
    
    ! Print statistics:
    
        do i = max_generic_timer + 1, max_amoeba_timer
          avg_time(i) = avg_time(i)/numtasks
        end do
    
        write(iunit_debug, 61) 'avg', &
                           avg_time(cit_setup_timer) + &
                           avg_time(build_list_timer), &
                           avg_time(dir_frc_sum_timer) + &
                           avg_time(adjust_masked_timer) + &
                           avg_time(pme_misc_timer), &
                           avg_time(bspline_timer) + &
                           avg_time(grid_charges_timer) + &
                           avg_time(scalar_sum_timer) + &
                           avg_time(grad_sum_timer) + &
                           avg_time(fft_timer), &
                           avg_time(atm_reassign_timer) + &
                           avg_time(img_reassign_timer) + &
                           avg_time(fft_slab_reassign_timer), &
                           avg_time(cit_setup_timer) + &
                           avg_time(build_list_timer) + &
                           avg_time(dir_frc_sum_timer) + &
                           avg_time(adjust_masked_timer) + &
                           avg_time(pme_misc_timer) + &
                           avg_time(bspline_timer) + &
                           avg_time(grid_charges_timer) + &
                           avg_time(scalar_sum_timer) + &
                           avg_time(grad_sum_timer) + &
                           avg_time(fft_timer) + &
                           avg_time(atm_reassign_timer) + &
                           avg_time(img_reassign_timer) + &
                           avg_time(fft_slab_reassign_timer)
    
        write(iunit_debug, '(49a)') ('-', i = 1, 49)
    
      else
    
    ! Send time stats to master node:
    
        call mpi_send(time_stats(max_generic_timer+1), &
                      max_amoeba_timer - max_generic_timer, &
                      mpi_double_precision, &
                      0, mytaskid, mpi_comm_world, err_code_mpi)
      end if

    end if

  else if (igb .eq. 1 .or. igb .eq. 2 .or. igb .eq. 5 .or. igb .eq. 7 .or. igb .eq. 20) then

70 format(i4, 5f10.1)
71 format(a4, 5f10.1)

    if (master) then
  
      fmt1 = '(a4, 5a10)'
  
      do i = max_generic_timer + 1, max_gb_timer
        avg_time(i) = time_stats(i)
      end do
  
      write(iunit_debug, '(/, a, /)') &
        'GB NonBond Parallel Profiling - NonSetup CPU Seconds:'
      write(iunit_debug, fmt1)' ', ' ', ' ', 'O', 'D', ' '
      write(iunit_debug, fmt1)' ', ' ', ' ', 'f', 'i', ' '
      write(iunit_debug, fmt1)' ', 'R', ' ', 'f', 's', 'T'
      write(iunit_debug, fmt1)'T', 'a', 'D', 'D', 't', 'o'
      write(iunit_debug, fmt1)'a', 'd', 'i', 'i', 'r', 't'
      write(iunit_debug, fmt1)'s', 'i', 'a', 'a', 'i', 'a'
      write(iunit_debug, fmt1)'k', 'i', 'g', 'g', 'b', 'l'
      write(iunit_debug, '(54a)') ('-', i = 1, 54)
      write(iunit_debug, 70) mytaskid, &
                         time_stats(calc_gb_rad_timer), &
                         time_stats(calc_gb_diag_timer), &
                         time_stats(calc_gb_offdiag_timer), &
                         time_stats(dist_gb_rad_timer), &
                         time_stats(calc_gb_rad_timer) + &
                         time_stats(calc_gb_diag_timer) + &
                         time_stats(calc_gb_offdiag_timer) + &
                         time_stats(dist_gb_rad_timer)

      do i = 1, numtasks - 1
  
  ! Receive time stat from each node in turn:
  
        call mpi_recv(time_stats_mpi_buf(max_generic_timer+1), &
                      max_gb_timer - max_generic_timer, &
                      mpi_double_precision, i, i, lib_mpi_comm, stat_mpi, &
                      err_code_mpi)
  
  ! Compute total time, write times to iunit_debug:
  
        do j = max_generic_timer + 1, max_gb_timer
          avg_time(j) = avg_time(j) + time_stats_mpi_buf(j)
        end do
  
      write(iunit_debug, 70) i, &
                         time_stats_mpi_buf(calc_gb_rad_timer), &
                         time_stats_mpi_buf(calc_gb_diag_timer), &
                         time_stats_mpi_buf(calc_gb_offdiag_timer), &
                         time_stats_mpi_buf(dist_gb_rad_timer), &
                         time_stats_mpi_buf(calc_gb_rad_timer) + &
                         time_stats_mpi_buf(calc_gb_diag_timer) + &
                         time_stats_mpi_buf(calc_gb_offdiag_timer) + &
                         time_stats_mpi_buf(dist_gb_rad_timer)
      end do
  
      write(iunit_debug, '(54a)') ('-', i = 1, 54)
  
  ! Print statistics:
  
      do i = max_generic_timer + 1, max_gb_timer
        avg_time(i) = avg_time(i)/numtasks
      end do
  
      write(iunit_debug, 71) 'avg', &
                         avg_time(calc_gb_rad_timer), &
                         avg_time(calc_gb_diag_timer), &
                         avg_time(calc_gb_offdiag_timer), &
                         avg_time(dist_gb_rad_timer), &
                         avg_time(calc_gb_rad_timer) + &
                         avg_time(calc_gb_diag_timer) + &
                         avg_time(calc_gb_offdiag_timer) + &
                         avg_time(dist_gb_rad_timer)
  
      write(iunit_debug, '(54a)') ('-', i = 1, 54)
  
    else
  
  ! Send time stats to master node:
  
      call mpi_send(time_stats(max_generic_timer+1), &
                    max_gb_timer - max_generic_timer, mpi_double_precision, &
                    0, mytaskid, lib_mpi_comm, err_code_mpi)
    end if

  end if

  return

end subroutine profile_parallel_cpu

!*******************************************************************************
!
! Subroutine:  zero_time
!
! Description:  Zero generic profiling timer used as base for timing each
!               step.  This should be done on entering code in each time
!               step of the simulations.
!*******************************************************************************

subroutine zero_time()

  implicit none

  call second(generic_time1)

  return

end subroutine zero_time

!*******************************************************************************
!
! Subroutine:  update_time
!
! Description:  Update generic profiling timer.
!              
!*******************************************************************************

subroutine update_time(timer_id)

  implicit none

  integer       :: timer_id

  call second(generic_time2)
  time_stats(timer_id) = time_stats(timer_id) + generic_time2 - generic_time1
  generic_time1 = generic_time2

  return

end subroutine update_time

!*******************************************************************************
!
! Subroutine:  zero_pme_time
!
! Description:  Zero pme profiling timer used as base for timing each
!               step.  This should be done on entering pme code in each time
!               step of the simulations.
!*******************************************************************************

subroutine zero_pme_time()

  implicit none

  call second(pme_time1)

  return

end subroutine zero_pme_time

!*******************************************************************************
!
! Subroutine:  update_pme_time
!
! Description:  Update master profiling timer.
!              
!*******************************************************************************

subroutine update_pme_time(timer_id)

  implicit none

  integer       :: timer_id

  call second(pme_time2)
  time_stats(timer_id) = time_stats(timer_id) + pme_time2 - pme_time1
  pme_time1 = pme_time2

  return

end subroutine update_pme_time

!*******************************************************************************
!
! Subroutine:  zero_gb_time
!
! Description:  Zero gb profiling timer used as base for timing each
!               step.  This should be done on entering gb code in each time
!               step of the simulations.
!*******************************************************************************

subroutine zero_gb_time()

  implicit none

  call second(gb_time1)

  return

end subroutine zero_gb_time

!*******************************************************************************
!
! Subroutine:  update_gb_time
!
! Description:  Update master profiling timer.
!              
!*******************************************************************************

subroutine update_gb_time(timer_id)

  implicit none

  integer       :: timer_id

  call second(gb_time2)
  time_stats(timer_id) = time_stats(timer_id) + gb_time2 - gb_time1
  gb_time1 = gb_time2

  return

end subroutine update_gb_time

! A crude but effective way to measure times, intended for temporary performance
! testing in development.

!*******************************************************************************
!
! Subroutine:  init_test_timers
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine init_test_timers

  use parallel_dat_mod

  implicit none

  call_cnt(:) = 0
  io_bytes(:) = 0.d0
  elapsed_cpu(:) = 0.d0
  elapsed_wall_sec(:) = 0
  elapsed_wall_usec(:) = 0

  ! The way this code currently works, since the master does i/o, if a call
  ! is not assigned an id string in the master, it is basically not assigned
  ! a string that will be printed.  A feature...  You will see the slave side
  ! of bcasts, recip force calls, etc. not getting id'd for this reason.
  if (master) then
    tt_id_str(:) = 'Not called from master'
  else
    tt_id_str(:) = ''
  end if

  return

end subroutine init_test_timers

!*******************************************************************************
!
! Subroutine:  enable_test_timers
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine enable_test_timers

  implicit none

  test_timers_enabled = .true.

  return

end subroutine enable_test_timers

!*******************************************************************************
!
! Subroutine:  disable_test_timers
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine disable_test_timers

  implicit none

  test_timers_enabled = .false.

  return

end subroutine disable_test_timers

!*******************************************************************************
!
! Subroutine:  start_test_timer
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine start_test_timer(timer_id, id_str, byte_cnt)

  implicit none

  integer               :: timer_id
  character(len=*)      :: id_str
  integer               :: byte_cnt

  if (.not. test_timers_enabled) return

  if (timer_id .gt. time_test_max_cnt) then
    write(error_msg,*) 'TIME_TEST TIMER_ID VALUE TOO HIGH!'
    call mol_mech_error
  end if

  io_bytes(timer_id) = io_bytes(timer_id) + dble(byte_cnt)

  if (call_cnt(timer_id) .eq. 0) tt_id_str(timer_id) = id_str

  call_cnt(timer_id) = call_cnt(timer_id) + 1

  call second(start_cpu(timer_id))

  call get_wall_time(start_wall_sec(timer_id), start_wall_usec(timer_id))

  return

end subroutine start_test_timer

!*******************************************************************************
!
! Subroutine:  stop_test_timer
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine stop_test_timer(timer_id)

  implicit none

  integer       timer_id

  double precision      stop_cpu
  integer               stop_wall_sec, stop_wall_usec

  if (.not. test_timers_enabled) return

  call second(stop_cpu)

  call get_wall_time(stop_wall_sec, stop_wall_usec)

  elapsed_cpu(timer_id) = elapsed_cpu(timer_id) + &
                          (stop_cpu - start_cpu(timer_id))

  if (stop_wall_usec .lt. start_wall_usec(timer_id)) then
    stop_wall_usec = stop_wall_usec + 1000000
    stop_wall_sec = stop_wall_sec - 1
  end if

  elapsed_wall_sec(timer_id) = elapsed_wall_sec(timer_id) + &
                               (stop_wall_sec - start_wall_sec(timer_id))

  elapsed_wall_usec(timer_id) = elapsed_wall_usec(timer_id) + &
                               (stop_wall_usec - start_wall_usec(timer_id))

  if (elapsed_wall_usec(timer_id) .ge. 1000000) then
    elapsed_wall_usec(timer_id) = elapsed_wall_usec(timer_id) - 1000000
    elapsed_wall_sec(timer_id) = elapsed_wall_sec(timer_id) + 1
  end if

  return

end subroutine stop_test_timer

!*******************************************************************************
!
! Subroutine:  update_io_cnt
!
! Description: For use when you don't know i/o count until after call completes.
!              
!*******************************************************************************

subroutine update_io_cnt(timer_id, byte_cnt)

  implicit none

  integer       :: timer_id
  integer       :: byte_cnt

  if (.not. test_timers_enabled) return

  if (timer_id .gt. time_test_max_cnt) then
    write(error_msg,*) 'TIME_TEST TIMER_ID VALUE TOO HIGH!'
    call mol_mech_error
  end if

  io_bytes(timer_id) = io_bytes(timer_id) + dble(byte_cnt)

  return

end subroutine update_io_cnt

!*******************************************************************************
!
! Subroutine:  print_test_timers
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine print_test_timers_mpi

  use parallel_dat_mod

  implicit none

! Formal arguments:

! Local variables:

  integer                       :: alloc_failed
  integer                       :: i
  integer                       :: taskid
  integer                       :: err_code
  integer                       :: stat_mpi(mpi_status_size)

  ! These allocated arrays are only needed in the master:

  double precision, allocatable :: elapsed_cpu_slv(:)
  integer, allocatable          :: wall_msec_slv(:)
  double precision, allocatable :: io_bytes_slv(:)
  integer, allocatable          :: call_cnt_slv(:)

  ! These allocated arrays are needed in all processes:

  integer, allocatable          :: wall_msec(:)

  ! This auto storage is only used in master, but we let all processors push
  ! the stack...

  double precision              :: total_elapsed_cpu
  integer                       :: total_wall_msec
  double precision              :: total_io_bytes
  integer                       :: total_call_cnt

  double precision              :: elapsed_cpu_tot(time_test_max_cnt)
  integer                       :: wall_msec_tot(time_test_max_cnt)
  double precision              :: io_bytes_tot(time_test_max_cnt)
  integer                       :: call_cnt_tot(time_test_max_cnt)
  double precision              :: elapsed_cpu_min(time_test_max_cnt)
  integer                       :: wall_msec_min(time_test_max_cnt)
  double precision              :: io_bytes_min(time_test_max_cnt)
  integer                       :: call_cnt_min(time_test_max_cnt)
  double precision              :: elapsed_cpu_max(time_test_max_cnt)
  integer                       :: wall_msec_max(time_test_max_cnt)
  double precision              :: io_bytes_max(time_test_max_cnt)
  integer                       :: call_cnt_max(time_test_max_cnt)

  if (mytaskid .eq. 0) then

    ! Allocate storage for receiving stuff from slaves; we use allocation
    ! because some mpi implementations have fits over buffers on the stack.

	if( allocated(elapsed_cpu_slv)) deallocate(elapsed_cpu_slv)
	if( allocated(wall_msec_slv)) deallocate(wall_msec_slv)
	if( allocated(io_bytes_slv)) deallocate(io_bytes_slv)
	if( allocated(call_cnt_slv)) deallocate(call_cnt_slv)
	if( allocated(wall_msec)) deallocate(wall_msec)

    allocate(elapsed_cpu_slv(time_test_max_cnt), &
             wall_msec_slv(time_test_max_cnt), &
             io_bytes_slv(time_test_max_cnt), &
             call_cnt_slv(time_test_max_cnt), &
             wall_msec(time_test_max_cnt), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    ! Report Master's own times, and start accumulating totals:

    total_elapsed_cpu = 0.d0
    total_wall_msec = 0
    total_io_bytes = 0.d0
    total_call_cnt = 0

    do i = 1, time_test_max_cnt
      wall_msec(i) = elapsed_wall_sec(i) * 1000 + elapsed_wall_usec(i) / 1000
      total_elapsed_cpu = total_elapsed_cpu + elapsed_cpu(i)
      total_wall_msec = total_wall_msec + wall_msec(i)
      total_io_bytes = total_io_bytes + io_bytes(i)
      total_call_cnt = total_call_cnt + call_cnt(i)
      if (call_cnt(i) .gt. 0) then
        write(mdout,100) tt_id_str(i), mytaskid, i, elapsed_cpu(i), &
                         wall_msec(i), io_bytes(i)/1000.d0, call_cnt(i)
      end if

      elapsed_cpu_tot(i) = elapsed_cpu(i)
      wall_msec_tot(i) = wall_msec(i)
      io_bytes_tot(i) = io_bytes(i)
      call_cnt_tot(i) = call_cnt(i)
  
      elapsed_cpu_min(i) = elapsed_cpu(i)
      wall_msec_min(i) = wall_msec(i)
      io_bytes_min(i) = io_bytes(i)
      call_cnt_min(i) = call_cnt(i)

      elapsed_cpu_max(i) = elapsed_cpu(i)
      wall_msec_max(i) = wall_msec(i)
      io_bytes_max(i) = io_bytes(i)
      call_cnt_max(i) = call_cnt(i)

    end do

    ! Get and report Slave's times, accumulating totals into Master's arrays:

    do taskid = 1, numtasks - 1

      call mpi_recv(elapsed_cpu_slv, time_test_max_cnt, &
                    mpi_double_precision, taskid, taskid, lib_mpi_comm, &
                    stat_mpi, err_code)

      call mpi_recv(wall_msec_slv, time_test_max_cnt, &
                    mpi_integer, taskid, taskid, lib_mpi_comm, &
                    stat_mpi, err_code)

      call mpi_recv(io_bytes_slv, time_test_max_cnt, &
                    mpi_double_precision, taskid, taskid, lib_mpi_comm, &
                    stat_mpi, err_code)

      call mpi_recv(call_cnt_slv, time_test_max_cnt, &
                    mpi_integer, taskid, taskid, lib_mpi_comm, &
                    stat_mpi, err_code)

      do i = 1, time_test_max_cnt
        if (call_cnt_slv(i) .gt. 0) then
          write(mdout,100) tt_id_str(i), taskid, i, elapsed_cpu_slv(i), &
                           wall_msec_slv(i), io_bytes_slv(i)/1000.d0, &
                           call_cnt_slv(i)

          elapsed_cpu_tot(i) = elapsed_cpu_tot(i) + elapsed_cpu_slv(i)
          elapsed_cpu_min(i) = min(elapsed_cpu_min(i), elapsed_cpu_slv(i))
          elapsed_cpu_min(i) = min(elapsed_cpu_min(i), elapsed_cpu_slv(i))

          wall_msec_tot(i) = wall_msec_tot(i) + wall_msec_slv(i)
          wall_msec_min(i) = min(wall_msec_min(i), wall_msec_slv(i))
          wall_msec_max(i) = max(wall_msec_max(i), wall_msec_slv(i))

          io_bytes_tot(i) = io_bytes_tot(i) + io_bytes_slv(i)
          io_bytes_min(i) = min(io_bytes_min(i), io_bytes_slv(i))
          io_bytes_max(i) = max(io_bytes_max(i), io_bytes_slv(i))

          call_cnt_tot(i) = call_cnt_tot(i) + call_cnt_slv(i)
          call_cnt_min(i) = min(call_cnt_min(i), call_cnt_slv(i))
          call_cnt_max(i) = max(call_cnt_max(i), call_cnt_slv(i))

          total_elapsed_cpu = total_elapsed_cpu + elapsed_cpu_slv(i)
          total_wall_msec = total_wall_msec + wall_msec_slv(i)
          total_io_bytes = total_io_bytes + io_bytes_slv(i)
          total_call_cnt = total_call_cnt + call_cnt_slv(i)
        end if
      end do

    end do

100 format(/, '| DBG: ', a, ':', &
           /, '|   Task =', i3, ', Call ID =', i3, ', CPU sec = ', f10.4, &
           ', Wall msec = ', i8, /, '|   IO KB = ', f14.1, ', Call cnt = ', i7)
       
    ! Write out totals avgs, mins, maxs for all calls:

    do i = 1, time_test_max_cnt
      if (call_cnt_tot(i) .gt. 0) then
        write(mdout,110) i, tt_id_str(i), &
                         call_cnt_tot(i), elapsed_cpu_tot(i), &
                         wall_msec_tot(i), io_bytes_tot(i)/1000.d0
        write(mdout,120) call_cnt_tot(i) / numtasks, &
                         elapsed_cpu_tot(i) / dble(numtasks), &
                         wall_msec_tot(i) / numtasks, &
                         io_bytes_tot(i) / (1000.d0 * dble(numtasks))
        write(mdout,130) call_cnt_min(i), elapsed_cpu_min(i), &
                         wall_msec_min(i), io_bytes_min(i)/1000.d0
        write(mdout,140) call_cnt_max(i), elapsed_cpu_max(i), &
                         wall_msec_max(i), io_bytes_max(i)/1000.d0
        ! write % of total wall time:

        write(mdout,150) i, &
          dble(wall_msec_tot(i)) / dble(total_wall_msec) * 100.d0

        ! write % of total i/o:

        write(mdout,160) i, &
          io_bytes_tot(i) / total_io_bytes * 100.d0

      end if
    end do

110 format(/, &
           '| Timer Call ID ', i3, ', Ident: ', a, /, &
           '|', /, &
           '|               Call cnt      CPU sec   Wall msec         IO KB', &
           /,&
           '| Totals:       ', i8, 3x, f10.4, 3x, i9, 3x, f11.1)
120 format('| Per Task Avg: ', i8, 3x, f10.4, 3x, i9, 3x, f11.1)
130 format('| Per Task Min: ', i8, 3x, f10.4, 3x, i9, 3x, f11.1)
140 format('| Per Task Max: ', i8, 3x, f10.4, 3x, i9, 3x, f11.1)
150 format('|', /, &
           '| Percent of job wall msec for call ', i3, ' = ', f6.2)
160 format('| Percent of job IO KB for call ', i3, ' =     ', f6.2)

  else

    ! Allocate storage for calculating wall_msec:

	if( allocated(wall_msec)) deallocate(wall_msec)

    allocate(wall_msec(time_test_max_cnt), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    ! Send time stats to master node:

    call mpi_send(elapsed_cpu, time_test_max_cnt, mpi_double_precision, 0, &
                  mytaskid, lib_mpi_comm, err_code)

    do i = 1, time_test_max_cnt
      wall_msec(i) = elapsed_wall_sec(i) * 1000 + elapsed_wall_usec(i) / 1000
    end do

    call mpi_send(wall_msec, time_test_max_cnt, mpi_integer, 0, &
                  mytaskid, lib_mpi_comm, err_code)

    call mpi_send(io_bytes, time_test_max_cnt, mpi_double_precision, 0, &
                  mytaskid, lib_mpi_comm, err_code)

    call mpi_send(call_cnt, time_test_max_cnt, mpi_integer, 0, &
                  mytaskid, lib_mpi_comm, err_code)

  end if

  return

end subroutine print_test_timers_mpi

subroutine print_test_timers

  implicit none

! Local variables:

  integer       :: i
  integer       :: wall_msec

  do i = 1, time_test_max_cnt

    wall_msec = elapsed_wall_sec(i) * 1000 + elapsed_wall_usec(i) / 1000

    if (call_cnt(i) .gt. 0) then
      write(mdout,100) tt_id_str(i), 0, i, elapsed_cpu(i), wall_msec, &
                   io_bytes(i)/1000.d0, call_cnt(i)
    end if

  end do

100 format( &
    /, '| DBG: ', a, ':', &
    /, '|   Task =', i3, ', Call ID =', i3, ', CPU sec = ', f10.4, &
       ', Wall msec = ', i8, /, '|   IO KB = ', f14.1, ', Call cnt = ', i7)

  return

end subroutine print_test_timers

end module timers_mod
