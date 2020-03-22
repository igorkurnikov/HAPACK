!*******************************************************************************
!
! Module: file_io_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module file_io_dat_mod

  implicit none

! File-related data.  No need to broadcast this stuff, since only the master
! does i/o:

  integer, parameter            :: max_fn_len = 80
  character(max_fn_len), save   :: mdout_name
  character(max_fn_len), save   :: mdcrd_name
  character(max_fn_len), save   :: mdvel_name
  character(max_fn_len), save   :: mden_name
  character(max_fn_len), save   :: restrt_name
  character(max_fn_len), save   :: prmtop_name
  character, save       :: owrite

  character(300), save   :: error_msg
  character(80), save    :: result_str

! Logical unit numbers for pmemd files.

  integer, parameter    :: mdin      =  35
  integer, parameter    :: mdout     =  36
  integer, parameter    :: prmtop    =  8
  integer, parameter    :: inpcrd    =  9
  integer, parameter    :: refc      = 10
  ! who was 11?
  integer, parameter    :: mdcrd     = 12
  integer, parameter    :: mdvel     = 13
  ! who was 14?
  integer, parameter    :: mden      = 15
  integer, parameter    :: restrt    = 16
  integer, parameter    :: restrt2   = 17
  
  ! Enumerations for crdfile_type:

  integer, parameter                    :: crdfile_type_undefined = 0
  integer, parameter                    :: crdfile_type_old = 1
  integer, parameter                    :: crdfile_type_new = 2

end module file_io_dat_mod

subroutine write_mdout(str)

use file_io_dat_mod
        
    implicit none
	
	character(*),intent(in) :: str
	
!	write(*,*) 'write_mdout(): mdout = ',mdout
!	write(*,*) str
	write(mdout,*)str

end subroutine write_mdout

subroutine format_int_fort(i,str_fmt)

use file_io_dat_mod
        
    implicit none
	
	integer     ,intent(in) :: i
	character(*),intent(in) :: str_fmt
	
!	write(*,*) 'format_int_fort'
!	write(*,*) str
	write(result_str,str_fmt) i

end subroutine format_int_fort

subroutine format_double_fort(x,str_fmt)

use file_io_dat_mod
        
    implicit none
	
	double precision    ,intent(in) :: x
	character(*),intent(in) :: str_fmt
	character(20)   ::str_fmt_2
	
!	str_fmt_2 = "(F8.3)"
	
!	write(*,*) 'format_double_fort'
!	write(*,*) str
	write(result_str,str_fmt) x

end subroutine format_double_fort


