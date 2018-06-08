!*******************************************************************************
!
! Internal Subroutine:  pack_nb_list
!
! Description: See get_nb_list comments above.
!              
!*******************************************************************************

subroutine pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                        full_eval_cnt, img_j_full_eval, ipairs, num_packed)

  implicit none

! Formal arguments:

  integer               :: ee_eval_cnt
  integer               :: img_j_ee_eval(*)
  integer               :: full_eval_cnt
  integer               :: img_j_full_eval(*)
  integer               :: ipairs(*)
  integer               :: num_packed

! Local variables:

  integer               :: i

! Check for enough room on the pairs list.  We also allow space for the
! two pairlist counters at the front of the list, and for one blank integer at
! the end of the list because some short_ene implementations may prefetch one
! integer and could thus run off the end of the list.  If DIRFRC_COMTRANS is 
! defined we also allow space for one translation flag at the front of the list.

#ifdef DIRFRC_COMTRANS

  if (num_packed + ee_eval_cnt + full_eval_cnt + 4 .le. ipairs_size) then

    if (common_tran) then
      ipairs(num_packed + 1) = 1
    else
      ipairs(num_packed + 1) = 0
    end if

    ipairs(num_packed + 2) = ee_eval_cnt
    ipairs(num_packed + 3) = full_eval_cnt

    if (ee_eval_cnt .gt. 0) then
      i = num_packed + 4
      ipairs(i:i + ee_eval_cnt - 1) = img_j_ee_eval(1:ee_eval_cnt)
    end if

    if (full_eval_cnt .gt. 0) then
      i = num_packed + 4 + ee_eval_cnt
      ipairs(i:i + full_eval_cnt - 1) = img_j_full_eval(1:full_eval_cnt)
    end if

    num_packed = num_packed + ee_eval_cnt + full_eval_cnt + 3

  else

    ifail = 1

  end if

#else

  if (num_packed + ee_eval_cnt + full_eval_cnt + 3 .le. ipairs_size) then

    ipairs(num_packed + 1) = ee_eval_cnt
    ipairs(num_packed + 2) = full_eval_cnt

    if (ee_eval_cnt .gt. 0) then
      i = num_packed + 3
      ipairs(i:i + ee_eval_cnt - 1) = img_j_ee_eval(1:ee_eval_cnt)
    end if

    if (full_eval_cnt .gt. 0) then
      i = num_packed + 3 + ee_eval_cnt
      ipairs(i:i + full_eval_cnt - 1) = img_j_full_eval(1:full_eval_cnt)
    end if

    num_packed = num_packed + ee_eval_cnt + full_eval_cnt + 2

  else

    ifail = 1

  end if

#endif /* DIRFRC_COMTRANS */

  return

end subroutine pack_nb_list

!*******************************************************************************
!
! Internal Subroutine:  pack_nb_list_skip_belly_pairs
!
! Description: See get_nb_list comments above.
!              
!*******************************************************************************

subroutine pack_nb_list_skip_belly_pairs(ee_eval_cnt, img_j_ee_eval, &
                                         full_eval_cnt, img_j_full_eval, &
                                         ipairs, igroup, img_atm_map, &
                                         num_packed)

  implicit none

! Formal arguments:

  integer               :: ee_eval_cnt
  integer               :: img_j_ee_eval(*)
  integer               :: full_eval_cnt
  integer               :: img_j_full_eval(*)
  integer               :: ipairs(*)
  integer               :: igroup(*)
  integer               :: img_atm_map(*)
  integer               :: num_packed

! Local variables:

  integer               :: i
  integer               :: img_j
  integer               :: new_ee_eval_cnt
  integer               :: new_full_eval_cnt
  integer, parameter    :: mask27 = Z"07FFFFFF"


  if (igroup(atm_i) .eq. 0) then

    new_ee_eval_cnt = 0
    new_full_eval_cnt = 0

#ifdef DIRFRC_COMTRANS
    if (common_tran) then

      do i = 1, ee_eval_cnt
        if (igroup(img_atm_map(img_j_ee_eval(i))) .ne. 0) then
          new_ee_eval_cnt = new_ee_eval_cnt + 1
          img_j_ee_eval(new_ee_eval_cnt) = img_j_ee_eval(i)
        end if
      end do

      do i = 1, full_eval_cnt
        if (igroup(img_atm_map(img_j_full_eval(i))) .ne. 0) then
          new_full_eval_cnt = new_full_eval_cnt + 1
          img_j_full_eval(new_full_eval_cnt) = img_j_full_eval(i)
        end if
      end do

    else
#endif

      do i = 1, ee_eval_cnt
        img_j = iand(img_j_ee_eval(i), mask27)
        if (igroup(img_atm_map(img_j)) .ne. 0) then
          new_ee_eval_cnt = new_ee_eval_cnt + 1
          img_j_ee_eval(new_ee_eval_cnt) = img_j_ee_eval(i)
        end if
      end do

      do i = 1, full_eval_cnt
        img_j = iand(img_j_full_eval(i), mask27)
        if (igroup(img_atm_map(img_j)) .ne. 0) then
          new_full_eval_cnt = new_full_eval_cnt + 1
          img_j_full_eval(new_full_eval_cnt) = img_j_full_eval(i)
        end if
      end do

#ifdef DIRFRC_COMTRANS
    end if
#endif

    ee_eval_cnt = new_ee_eval_cnt
    full_eval_cnt = new_full_eval_cnt

  end if

  call pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                    full_eval_cnt, img_j_full_eval, &
                    gbl_ipairs, num_packed)

  return

end subroutine pack_nb_list_skip_belly_pairs
