#include "copyright.i"
!*******************************************************************************
!
! Module: file_io_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module file_io_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  amopen
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine amopen(lun, fname, fstat, ffmt, facc)

  use parallel_dat_mod
  use file_io_dat_mod

  implicit none

! Formal arguments:

  integer               :: lun     ! logical unit number
  character(*)          :: fname   ! file name 
  character(1)          :: fstat   ! status code:
                                   ! "N", "O", or "U" = new, old, unk.
  character(1)          :: ffmt    ! format code: "U", "F" = unform., form.
  character(1)          :: facc    ! access code:
                                   ! "R", "W", "A" = read, read/write, append
! Local variables:

  integer               :: ios     ! i/o status variable
  character(11)         :: kform   ! form keyword
  character(7)          :: stat    ! status keyword
 
  if (fstat .eq. 'N') then
    stat = 'NEW'
  else if (fstat .eq. 'O') then
    stat = 'OLD'
  else if (fstat .eq. 'U') then
    stat = 'UNKNOWN'
  else
    write(error_msg, '(/,2x,a,i4)') 'amopen: bogus fstat, unit ', lun
    call mol_mech_error
  end if

  if (ffmt .eq. 'U') then
    kform = 'UNFORMATTED'
  else if (ffmt .eq. 'F') then
    kform = 'FORMATTED'
  else
    write(error_msg, '(/,2x,a,i4)') 'amopen: bogus ffmt, unit', lun
    call mol_mech_error
  end if

  open(lun, file=fname, status=stat, form=kform, iostat=ios)
 
  if (ios .ne. 0) then

    if (lun .eq. mdout) then
#ifndef DUMB
      write(0, '(/,2x,a,i4,a,a)') 'Unit ', lun, ' Error on OPEN: ', fname
#endif
      close(0)
    else
      write(mdout, '(/,2x,a,i4,a,a)') 'Unit ', lun, ' Error on OPEN: ', fname
      close(unit = mdout)
#ifndef DUMB
      write(0, '(/,2x,a,i4,a,a)') 'Unit ', lun, ' Error on OPEN: ', fname
#endif
      close(0)
    end if
    
    call mol_mech_error

  end if

  rewind(lun)

  return

end subroutine amopen

end module file_io_mod

