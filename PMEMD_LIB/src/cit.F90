#include "copyright.i"

!*******************************************************************************
!
! Module:  cit_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module cit_mod

  implicit none

  type cit_tbl_rec
    integer             :: img_lo
    integer             :: img_hi
  end type cit_tbl_rec

  type atm_lst_rec
    integer             :: idx
    integer             :: nxt
  end type atm_lst_rec

  integer, save         :: cit_tbl_x_dim
  integer, save         :: cit_tbl_y_dim
  integer, save         :: cit_tbl_z_dim

contains

!*******************************************************************************
!
! Subroutine:  setup_crd_idx_tbl
! 
! Description: partition atoms into rectangular buckets within unit cell
!              atm_list - lists of atoms that belongs to different buckets
!*******************************************************************************

subroutine setup_crd_idx_tbl(atm_cnt, fraction, crd_idx_lst_tbl, atm_lst)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: fraction(3, atm_cnt)
  integer               :: crd_idx_lst_tbl(0 : cit_tbl_x_dim - 1, &
                                           0 : cit_tbl_y_dim - 1, &
                                           0 : cit_tbl_z_dim - 1)

  ! We 0-base the following array for efficiency, but don't use atm_lst(0) so
  ! we can use 0 as the distinguished non-entry (NULL).

  type(atm_lst_rec)     :: atm_lst(0:atm_cnt)

! Local variables:

  integer               :: atm_id
  integer               :: last_idx
  integer               :: nxt_atm_lst_idx
  integer               :: x_idx, y_idx, z_idx
  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z

  integer               :: tail_idx_tbl(0 : cit_tbl_x_dim - 1, &
                                        0 : cit_tbl_y_dim - 1, &
                                        0 : cit_tbl_z_dim - 1)

! Pre-initialize as needed:

  nxt_atm_lst_idx = 1

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

  crd_idx_lst_tbl(:,:,:) = 0                ! Marks empty entries 

! Load the atom id's:

  do atm_id = 1, atm_cnt

    x_idx = int(fraction(1, atm_id) * scale_fac_x)
    y_idx = int(fraction(2, atm_id) * scale_fac_y)
    z_idx = int(fraction(3, atm_id) * scale_fac_z)

    if (crd_idx_lst_tbl(x_idx, y_idx, z_idx) .eq. 0) then   ! New list.

      crd_idx_lst_tbl(x_idx, y_idx, z_idx) = nxt_atm_lst_idx

    else        ! List already started.  Follow the chain and add a node:

      last_idx = tail_idx_tbl(x_idx, y_idx, z_idx)
      atm_lst(last_idx)%nxt = nxt_atm_lst_idx

    end if

    ! Add the node:

    atm_lst(nxt_atm_lst_idx)%idx = atm_id
    atm_lst(nxt_atm_lst_idx)%nxt = 0
    tail_idx_tbl(x_idx, y_idx, z_idx) = nxt_atm_lst_idx
    nxt_atm_lst_idx = nxt_atm_lst_idx + 1

  end do

  return

end subroutine setup_crd_idx_tbl

!*******************************************************************************
!
! Subroutine:  setup_cit_tbl_dims
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine set_cit_tbl_dims(box, list_cutoff)

  implicit none

! Formal arguments:

  double precision, intent(in)  :: box(3)
  double precision, intent(in)  :: list_cutoff

  cit_tbl_x_dim = max(3, floor(box(1) / (list_cutoff * 0.5d0)))
  cit_tbl_y_dim = max(3, floor(box(2) / (list_cutoff * 0.5d0)))
  cit_tbl_z_dim = max(3, floor(box(3) / (list_cutoff * 0.5d0)))

  return

end subroutine set_cit_tbl_dims

!*******************************************************************************
!
! Subroutine:  setup_cit
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine setup_cit(atm_cnt, box, fraction, crd_idx_tbl, img, atm_img_map, &
                     img_atm_map, atm_iac, img_iac, charge)

  use img_mod
  use pbc_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: box(3)
  double precision      :: fraction(3, atm_cnt)
  type(cit_tbl_rec)     :: crd_idx_tbl(0 : cit_tbl_x_dim - 1, &
                                       0 : cit_tbl_y_dim - 1, &
                                       0 : cit_tbl_z_dim - 1)

  type(img_rec)         :: img(atm_cnt)
  integer               :: atm_img_map(atm_cnt)
  integer               :: img_atm_map(atm_cnt)
  integer               :: atm_iac(atm_cnt)
  integer               :: img_iac(atm_cnt)
  double precision      :: charge(atm_cnt)

! Local variables:

  integer               :: atm_id, img_id
  integer               :: i, j, k
  integer               :: is_orthog_stk

  integer               :: img_lo, img_hi  ! used only in MPI

  integer               :: nxt_idx
  double precision      :: box_stk(3)
  double precision      :: ucell_stk(3, 3)
  double precision      :: f1, f2, f3
  integer               :: crd_idx_lst_tbl(0 : cit_tbl_x_dim - 1, &
                                           0 : cit_tbl_y_dim - 1, &
                                           0 : cit_tbl_z_dim - 1)

  type(atm_lst_rec)     :: atm_lst(0:atm_cnt)

  call setup_crd_idx_tbl(atm_cnt, fraction, crd_idx_lst_tbl, atm_lst)

  is_orthog_stk = is_orthog

  if (is_orthog_stk .ne. 0) then
    box_stk(:) = box(:)
  else
    ucell_stk(:, :) = ucell(:, :)
  end if
  
  if( numtasks .gt. 1) then  ! MPI code

    img_hi = 0

    do k = 0, cit_tbl_z_dim - 1
      do j = 0, cit_tbl_y_dim - 1
        do i = 0, cit_tbl_x_dim - 1

          nxt_idx = crd_idx_lst_tbl(i, j, k)

          if (nxt_idx .ne. 0) then

            img_hi = img_hi + 1
            img_lo = img_hi

            do
          ! BEGIN DBG
          ! if (atm_lst(nxt_idx)%idx .lt. 1 .or. &
          !     atm_lst(nxt_idx)%idx .gt. natom) then
          !   write(0,*)'DBG: Found bad value ', atm_lst(nxt_idx)%idx, &
          !             'in atm_lst!!!'
          ! end if
          ! END DBG
              atm_img_map(atm_lst(nxt_idx)%idx) = img_hi
              nxt_idx = atm_lst(nxt_idx)%nxt
              if (nxt_idx .ne. 0) then
                img_hi = img_hi + 1
              else
                exit
              end if
            end do

            crd_idx_tbl(i, j, k)%img_lo = img_lo
            crd_idx_tbl(i, j, k)%img_hi = img_hi

            nxt_idx = crd_idx_lst_tbl(i, j, k)

          ! Images not yet claimed, as indicated by negative atm_id.

            do img_id = img_lo, img_hi
              img_atm_map(img_id) = - atm_lst(nxt_idx)%idx
              nxt_idx = atm_lst(nxt_idx)%nxt
            end do

          else ! (nxt_idx == 0)
            crd_idx_tbl(i, j, k)%img_lo = 0
            crd_idx_tbl(i, j, k)%img_hi = -1
          end if ! (ntx_idx != 0)

        end do  ! i = 0, cit_tbl_x_dim - 1
      end do  ! j = 0, cit_tbl_y_dim - 1
    end do   ! k = 0, cit_tbl_z_dim - 1
    
  ! Now assign the images you own to yourself.
  
    do img_id = my_img_lo, my_img_hi

      atm_id = - img_atm_map(img_id)
      img_atm_map(img_id) = atm_id

      img(img_id)%qterm = charge(atm_id)
      img_iac(img_id) = atm_iac(atm_id)

      if (is_orthog_stk .ne. 0) then

        img(img_id)%x = fraction(1, atm_id) * box_stk(1)
        img(img_id)%y = fraction(2, atm_id) * box_stk(2)
        img(img_id)%z = fraction(3, atm_id) * box_stk(3)

      else

        f1 = fraction(1, atm_id)
        f2 = fraction(2, atm_id)
        f3 = fraction(3, atm_id)

! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
! we can simplify the expression in this critical inner loop

!******************************************************************************
!     img(img_id)%y = f1 * ucell_stk(2, 1) + &
!                     f2 * ucell_stk(2, 2) + &
!                     f3 * ucell_stk(2, 3)
!     img(img_id)%z = f1 * ucell_stk(3, 1) + &
!                     f2 * ucell_stk(3, 2) + &
!                     f3 * ucell_stk(3, 3)
!******************************************************************************

        img(img_id)%x = f1 * ucell_stk(1, 1) + &
                        f2 * ucell_stk(1, 2) + &
                        f3 * ucell_stk(1, 3)
        img(img_id)%y = f2 * ucell_stk(2, 2) + &
                        f3 * ucell_stk(2, 3)
        img(img_id)%z = f3 * ucell_stk(3, 3)

      end if

    end do  ! (img_id)

  ! Under MPI we also use the img_iac value as a flag indicating that the
  ! other values have been set.  A zero value indicates "unset", and is
  ! purely used in list building, where we don't want to claim an image until
  ! we know we are using it, but also don't want to be reinitializing the 
  ! img and img_iac values needlessly.

    img_iac(1 : my_img_lo - 1) = 0
    img_iac(my_img_hi + 1 : atm_cnt) = 0

  else !  (numtasks == 1)

    img_id = 0

    do k = 0, cit_tbl_z_dim - 1
      do j = 0, cit_tbl_y_dim - 1
        do i = 0, cit_tbl_x_dim - 1

          nxt_idx = crd_idx_lst_tbl(i, j, k)

          if (nxt_idx .ne. 0) then

            img_id = img_id + 1
            crd_idx_tbl(i, j, k)%img_lo = img_id

            do

              atm_id = atm_lst(nxt_idx)%idx
              atm_img_map(atm_id) = img_id
              img_atm_map(img_id) = atm_id

              img(img_id)%qterm = charge(atm_id)
              img_iac(img_id) = atm_iac(atm_id)

              if (is_orthog_stk .ne. 0) then

                img(img_id)%x = fraction(1, atm_id) * box_stk(1)
                img(img_id)%y = fraction(2, atm_id) * box_stk(2)
                img(img_id)%z = fraction(3, atm_id) * box_stk(3)

              else

                f1 = fraction(1, atm_id)
                f2 = fraction(2, atm_id)
                f3 = fraction(3, atm_id)

! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
! we can simplify the expression in this critical inner loop

!******************************************************************************
!             img(img_id)%y = f1 * ucell_stk(2, 1) + &
!                             f2 * ucell_stk(2, 2) + &
!                             f3 * ucell_stk(2, 3)
!             img(img_id)%z = f1 * ucell_stk(3, 1) + &
!                             f2 * ucell_stk(3, 2) + &
!                             f3 * ucell_stk(3, 3)
!******************************************************************************

                img(img_id)%x = f1 * ucell_stk(1, 1) + &
                                f2 * ucell_stk(1, 2) + &
                                f3 * ucell_stk(1, 3)
                img(img_id)%y = f2 * ucell_stk(2, 2) + &
                                f3 * ucell_stk(2, 3)
                img(img_id)%z = f3 * ucell_stk(3, 3)

              end if

              nxt_idx = atm_lst(nxt_idx)%nxt

              if (nxt_idx .ne. 0) then
                img_id = img_id + 1
              else
                crd_idx_tbl(i, j, k)%img_hi = img_id
                exit
              end if
            end do
          else !  (nxt_idx == 0)
            crd_idx_tbl(i, j, k)%img_lo = 0
            crd_idx_tbl(i, j, k)%img_hi = -1
          end if
        end do  !  i
      end do  !  j 
    end do  !  k

  end if !  (numtasks > 1)

  return

end subroutine setup_cit

!*******************************************************************************
!
! Subroutine:  amoeba_setup_cit_master
!
! Description: <TBS>
!              
!*******************************************************************************
!
! BUGBUG - Temporary code for mapping all images in the master in amoeba.  This
!          is done to take care of the valence terms, as a temporary measure.

subroutine amoeba_setup_cit_master(atm_cnt, box, fraction, crd_idx_tbl, img, &
                                   atm_img_map, img_atm_map, atm_iac, img_iac, &
                                   qterm)

  use img_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: box(3)
  double precision      :: fraction(3, atm_cnt)
  type(cit_tbl_rec)     :: crd_idx_tbl(0 : cit_tbl_x_dim - 1, &
                                       0 : cit_tbl_y_dim - 1, &
                                       0 : cit_tbl_z_dim - 1)

  type(img_rec)         :: img(atm_cnt)
  integer               :: atm_img_map(atm_cnt)
  integer               :: img_atm_map(atm_cnt)
  integer               :: atm_iac(atm_cnt)
  integer               :: img_iac(atm_cnt)
  double precision      :: qterm(atm_cnt)

! Local variables:

  integer               :: atm_id, img_id
  integer               :: i, j, k
  integer               :: is_orthog_stk
  integer               :: nxt_idx
  double precision      :: box_stk(3)
  double precision      :: ucell_stk(3, 3)
  double precision      :: f1, f2, f3
  integer               :: crd_idx_lst_tbl(0 : cit_tbl_x_dim - 1, &
                                           0 : cit_tbl_y_dim - 1, &
                                           0 : cit_tbl_z_dim - 1)

  type(atm_lst_rec)     :: atm_lst(0:atm_cnt)
  
  call setup_crd_idx_tbl(atm_cnt, fraction, crd_idx_lst_tbl, atm_lst) ! partition atoms into rectagulare cells within periodical box

  is_orthog_stk = is_orthog

  if (is_orthog_stk .ne. 0) then
    box_stk(:) = box(:)
  else
    ucell_stk(:, :) = ucell(:, :)
  end if

  img_id = 0

  do k = 0, cit_tbl_z_dim - 1
    do j = 0, cit_tbl_y_dim - 1
      do i = 0, cit_tbl_x_dim - 1

        nxt_idx = crd_idx_lst_tbl(i, j, k)

        if (nxt_idx .ne. 0) then

          img_id = img_id + 1
          crd_idx_tbl(i, j, k)%img_lo = img_id

          do

            atm_id = atm_lst(nxt_idx)%idx
            atm_img_map(atm_id) = img_id
            img_atm_map(img_id) = atm_id

            img(img_id)%qterm = qterm(atm_id)
            img_iac(img_id) = atm_iac(atm_id)

            if (is_orthog_stk .ne. 0) then

              img(img_id)%x = fraction(1, atm_id) * box_stk(1)
              img(img_id)%y = fraction(2, atm_id) * box_stk(2)
              img(img_id)%z = fraction(3, atm_id) * box_stk(3)

            else

              f1 = fraction(1, atm_id)
              f2 = fraction(2, atm_id)
              f3 = fraction(3, atm_id)

! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
! we can simplify the expression in this critical inner loop

!******************************************************************************
!             img(img_id)%y = f1 * ucell_stk(2, 1) + &
!                             f2 * ucell_stk(2, 2) + &
!                             f3 * ucell_stk(2, 3)
!             img(img_id)%z = f1 * ucell_stk(3, 1) + &
!                             f2 * ucell_stk(3, 2) + &
!                             f3 * ucell_stk(3, 3)
!******************************************************************************

              img(img_id)%x = f1 * ucell_stk(1, 1) + &
                              f2 * ucell_stk(1, 2) + &
                              f3 * ucell_stk(1, 3)
              img(img_id)%y = f2 * ucell_stk(2, 2) + &
                              f3 * ucell_stk(2, 3)
              img(img_id)%z = f3 * ucell_stk(3, 3)

            end if

            nxt_idx = atm_lst(nxt_idx)%nxt

            if (nxt_idx .ne. 0) then
              img_id = img_id + 1
            else
              crd_idx_tbl(i, j, k)%img_hi = img_id
              exit
            end if

          end do

        else
          crd_idx_tbl(i, j, k)%img_lo = 0
          crd_idx_tbl(i, j, k)%img_hi = -1
        end if

      end do
    end do
  end do
  
  return

end subroutine amoeba_setup_cit_master

end module cit_mod
