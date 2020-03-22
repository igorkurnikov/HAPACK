#include "copyright.i"

!*******************************************************************************
!
! Module:  img_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module img_mod

  implicit none

  ! Global type declarations:

  type img_rec
    double precision    :: x
    double precision    :: y
    double precision    :: z
    double precision    :: qterm ! charge in amber pme, sq_polinv in amoeba - modified from charge in PMEMD 10
!    double precision    :: charge
  end type img_rec

  integer, save, allocatable                    :: gbl_atm_img_map(:)
  integer, save, allocatable                    :: gbl_img_atm_map(:)
  type(img_rec), save, allocatable              :: gbl_img(:)
  integer, save, allocatable                    :: gbl_img_iac(:)
  integer, save, allocatable                    :: gbl_excl_img_flags(:)

  ! my_img_* covers the range of images you "own"; ie. forces are accumulated
  ! in the local process for those images.  These are assigned WITHOUT wrapping.
  ! used_img_* covers the range of images you "use"; ie. you may need
  ! coordinate information and will report nonbonded force information to the
  ! owner for some images in this range. This range may wrap through natom to 1.

  integer, save         :: my_img_lo, my_img_hi    ! low and high indexes of atoms that belong to a processor 
  integer, save         :: used_img_lo, used_img_hi

  logical, save         :: used_img_range_wraps ! from lo through natom to hi

contains

!*******************************************************************************
!
! Subroutine:  alloc_img_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_img_mem(atm_cnt)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt

! Local variables:

  integer               :: alloc_failed

  if( numtasks .eq. 1 ) then
  ! In an MPI version, this stuff will be set elsewhere.
      my_img_lo = 1
      my_img_hi = atm_cnt
      used_img_lo = 1
      used_img_hi = atm_cnt
  endif

  if( allocated(gbl_atm_img_map)) deallocate(gbl_atm_img_map)
  if( allocated(gbl_img_atm_map)) deallocate(gbl_img_atm_map)
  if( allocated(gbl_img)) deallocate(gbl_img)
  if( allocated(gbl_img_iac)) deallocate(gbl_img_iac)
  if( allocated(gbl_excl_img_flags)) deallocate(gbl_excl_img_flags)

  allocate(gbl_atm_img_map(atm_cnt), &
           gbl_img_atm_map(atm_cnt), &
           gbl_img(atm_cnt), &
           gbl_img_iac(atm_cnt), &
           gbl_excl_img_flags(atm_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  gbl_img = img_rec(0.d0,0.d0,0.d0,0.d0)

  gbl_atm_img_map(:) = 0
  gbl_img_atm_map(:) = 0
  gbl_img_iac(:) = 0
  gbl_excl_img_flags(:) = 0

  return

end subroutine alloc_img_mem

!*******************************************************************************
!
! Subroutine:  bcast_img_dat
!
! Description: <TBS>
!              defined only in MPI             
!*******************************************************************************

subroutine bcast_img_dat(atm_cnt)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer       :: atm_cnt

  ! Nothing to broadcast.  We just allocate storage in the non-master nodes.

  if (.not. master) then
    call alloc_img_mem(atm_cnt)
  end if

  ! The allocated data is not initialized from the master node.

  return

end subroutine bcast_img_dat

!*******************************************************************************
!
! Subroutine:  fill_tranvec
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine fill_tranvec(tranvec)

  use pbc_mod

  implicit none

  double precision      :: tranvec(3, 18) ! this is (1:3,0:17) externally.
  integer               :: iv, i0, i1, i2, i3

! This works for both orthogonal and nonorthogonal unit cells.

  iv = 0

  do i3 = -1, 0
    do i2 = -1, 1
      do i1 = -1, 1
        iv = iv + 1
        do i0 = 1, 3
          tranvec(i0, iv) = i1 * ucell(i0, 1) + &
                            i2 * ucell(i0, 2) + &
                            i3 * ucell(i0, 3)
        end do
      end do
    end do
  end do

  return

end subroutine fill_tranvec


!*******************************************************************************
!
! Subroutine:  find_img_range
!
! Description:  Find the range of images this process uses.
!               Used only in MPI
!*******************************************************************************

subroutine find_img_range(img_cnt, img_atm_map)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: img_cnt
  integer               :: img_atm_map(img_cnt)

! Local variables:

  integer               :: img_i
  integer               :: run_lo, run_hi, run_cnt
  integer               :: lo_run_lo, lo_run_hi, lo_run_cnt
  integer               :: hi_run_lo, hi_run_hi, hi_run_cnt
  logical               :: in_run

  run_cnt = 0
  lo_run_cnt = 0

  in_run = .false.
  
! Find the longest unused run of images with id's below my_img_lo:

  do img_i = 1, my_img_lo

    if (img_atm_map(img_i) .lt. 0) then   ! img not used
      if (in_run) then
        run_hi = img_i
      else                
        run_lo = img_i
        run_hi = img_i
        in_run = .true.
      end if
    else                                  ! img used
      if (in_run) then
        run_cnt = run_hi - run_lo + 1
        if (run_cnt .gt. lo_run_cnt) then
          lo_run_cnt = run_cnt
          lo_run_lo = run_lo
          lo_run_hi = run_hi
        end if
        in_run = .false.
      end if
    end if

  end do

  if (in_run) then
    run_cnt = run_hi - run_lo + 1
    if (run_cnt .gt. lo_run_cnt) then
      lo_run_cnt = run_cnt
      lo_run_lo = run_lo
      lo_run_hi = run_hi
    end if
  end if
  
  run_cnt = 0
  hi_run_cnt = 0

  in_run = .false.
  
! Find the longest unused run of images with id's above my_img_hi:

  do img_i = my_img_hi + 1, img_cnt

    if (img_atm_map(img_i) .lt. 0) then   ! img not used
      if (in_run) then
        run_hi = img_i
      else                
        run_lo = img_i
        run_hi = img_i
        in_run = .true.
      end if
    else                                  ! img used
      if (in_run) then
        run_cnt = run_hi - run_lo + 1
        if (run_cnt .gt. hi_run_cnt) then
          hi_run_cnt = run_cnt
          hi_run_lo = run_lo
          hi_run_hi = run_hi
        end if
        in_run = .false.
      end if
    end if

  end do

  if (in_run) then
    run_cnt = run_hi - run_lo + 1
    if (run_cnt .gt. hi_run_cnt) then
      hi_run_cnt = run_cnt
      hi_run_lo = run_lo
      hi_run_hi = run_hi
    end if
  end if

! Now check and see if the two ranges are joined.  If so, we can combine them.
! The algorithm used in this subroutine could miss an opportunity to join
! ranges (due to finding a larger lo or hi subrange away from the boundary), but
! this is fairly unlikely, will not affect correctness of results, and will
! probably have a pretty small impact on performance (it will cause the used
! range to be slightly larger than it actually is)

  if (lo_run_cnt + hi_run_cnt .eq. 0) then      ! highly unlikely, but...

    used_img_lo = 1
    used_img_hi = img_cnt
    used_img_range_wraps = .false.

  else if (lo_run_cnt .gt. 0 .and. &
           hi_run_cnt .gt. 0 .and. &
           lo_run_lo .eq. 1 .and. &
           hi_run_hi .eq. img_cnt) then         ! unused ranges fuse...

    used_img_lo = lo_run_hi + 1
    used_img_hi = hi_run_lo - 1
    used_img_range_wraps = .false.

  else if (lo_run_cnt .ge. hi_run_cnt) then     ! biggest unused range lo

    used_img_lo = lo_run_hi + 1
    used_img_hi = lo_run_lo - 1

    if (used_img_hi .le. 0) then
      used_img_hi = img_cnt
      used_img_range_wraps = .false.
    else
      used_img_range_wraps = .true.
    end if

  else
   
    used_img_lo = hi_run_hi + 1
    used_img_hi = hi_run_lo - 1

    if (used_img_lo .gt. img_cnt) then
      used_img_lo = 1
      used_img_range_wraps = .false.
    else
      used_img_range_wraps = .true.
    end if

  end if

! BEGIN DBG
! write(0,*)'DBG: task, used_img_lo,hi=', mytaskid, used_img_lo, used_img_hi
! write(0,*)'DBG: task, used_img_range_wraps=', mytaskid, used_img_range_wraps
! END DBG

  return

end subroutine find_img_range

end module img_mod
