#include "copyright.i"

!*******************************************************************************
!
! Module: amoeba_flags_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module amoeba_flags_mod

  implicit none

  integer, parameter    :: proceed = 3
  integer, parameter    :: proceed_induce = 5
  integer, parameter    :: proceed_postinduce = 9

 ! Values for bits of do_flags:

  integer, parameter    :: valid_bit = 0
  integer, parameter    :: user_bit = 1
  integer, parameter    :: user_induce_bit = 2
  integer, parameter    :: user_postinduce_bit = 3

contains

!*******************************************************************************!
! Subroutine:  set_user_bit
!
! Description: <TBS>
!
!*******************************************************************************
subroutine set_user_bit(do_this, do_flag)

  implicit none

  integer, intent(in)   :: do_this
  integer, intent(out)  :: do_flag

  if (do_this .eq. 1) then
    do_flag = ibset(do_flag, user_bit)
  else
    do_flag = ibclr(do_flag, user_bit)
  end if

  return

end subroutine set_user_bit

end module amoeba_flags_mod
