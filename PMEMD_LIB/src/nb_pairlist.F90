#include "copyright.i"

!*******************************************************************************
!
! Module:  nb_pairlist_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module nb_pairlist_mod

  implicit none

! Data that is not broadcast:

  double precision, save, allocatable   :: gbl_saved_imgcrd(:,:)
  integer, allocatable, save            :: gbl_ipairs(:)
  integer, save                         :: ipairs_size

contains

!*******************************************************************************
!
! Subroutine:  alloc_nb_pairlist_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_nb_pairlist_mem(atm_cnt, cutlist)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  double precision, intent(in)  :: cutlist

! Local variables:

  integer                       :: alloc_failed
  integer, parameter            :: ipairs_size_min = 1000
  double precision, parameter   :: ipairs_size_coef = 0.225d0 
 
  if( allocated(gbl_saved_imgcrd))  deallocate(gbl_saved_imgcrd)

  allocate(gbl_saved_imgcrd(3, atm_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

! Allocate pairs list:
    
! The following ipairs_size calc is a crude heuristic assuming that the
! number of nonbonded pairs roughly scales with the cutoff volume and
! the number of atoms in the system.  If MPI is running, the total is
! divided up among the processes.  The 2 or 3 is for the counters and flags at
! the head of each atom sublist.

!write(*,*) "nb_pairlist_mod::nb_pairlist() ipairs_size_coef = ",ipairs_size_coef
!write(*,*) "nb_pairlist_mod::nb_pairlist() cutlist = ",cutlist
!write(*,*) "nb_pairlist_mod::nb_pairlist() atm_cnt = ",atm_cnt

#ifdef DIRFRC_COMTRANS
  ipairs_size = int((ipairs_size_coef * cutlist**3 + 3.d0) * dble(atm_cnt))
#else
  ipairs_size = int((ipairs_size_coef * cutlist**3 + 2.d0) * dble(atm_cnt))
#endif

!write(*,*) "nb_pairlist_mod::nb_pairlist() pt 2 ipairs_size = ",ipairs_size

  ipairs_size = ipairs_size / numtasks

  if (ipairs_size .lt. ipairs_size_min) ipairs_size = ipairs_size_min

 ! write(*,*) "nb_pairlist_mod::nb_pairlist() ipairs_size = ",ipairs_size

  if( allocated(gbl_ipairs)) deallocate(gbl_ipairs)  
  allocate(gbl_ipairs(ipairs_size), stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error

  return

end subroutine alloc_nb_pairlist_mem
  
!*******************************************************************************
!
! Subroutine:  bcast_nb_pairlist_dat
!
! Description: <TBS>
! Used only in MPI             
!*******************************************************************************

subroutine bcast_nb_pairlist_dat(atm_cnt, cutlist)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  double precision, intent(in)  :: cutlist
  
  if (.not. master) then
    call alloc_nb_pairlist_mem(atm_cnt, cutlist)
  end if

  ! The allocated data is not initialized from the master node.

  return

end subroutine bcast_nb_pairlist_dat

!*******************************************************************************
!
! Subroutine:   get_nb_list
!
! Description:  Main routine for list building.  Note it is structured for
!               parallelism.
!
!*******************************************************************************

subroutine get_nb_list(atm_cnt, box, flat_cit, img, img_atm_map, excl_img_pairlst, &
                       ico, fraction, tranvec, atm_maskdata, atm_mask, &
                       atm_img_map, excl_img_flags, &
                       iac, img_iac, charge, igroup, ntypes, ibelly, &
                       es_cutoff, vdw_cutoff, skinnb, verbose, ifail)

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: box(3)
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)

  type(img_rec)         :: img(atm_cnt)
  integer               :: img_atm_map(atm_cnt)
  integer               :: excl_img_pairlst(*)
  integer               :: ico(*)
  double precision      :: fraction(3, atm_cnt)
  double precision      :: tranvec(1:3, 0:17)
  type(maskdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
  integer               :: atm_img_map(atm_cnt)
  integer               :: excl_img_flags(atm_cnt)
  integer               :: iac(atm_cnt)
  integer               :: img_iac(atm_cnt)
  double precision      :: charge(atm_cnt)
  integer               :: igroup(atm_cnt)
  integer               :: ntypes
  integer               :: ibelly
  double precision      :: es_cutoff
  double precision      :: vdw_cutoff
  double precision      :: skinnb
  integer               :: verbose
  integer               :: ifail

! Local variables:

  double precision      :: x_box, y_box, z_box
  double precision      :: ucell_stk(3, 3)
  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
  double precision      :: cutlist_sq
  double precision      :: cut_x_frac, cut_y_frac, cut_z_frac
  double precision      :: x_i, y_i, z_i
  integer               :: atm_i, img_i
  integer               :: i
  integer               :: atm_mask_idx
  integer               :: excl_img_j_cnt
  integer               :: excl_pairlst_len
  integer               :: num_packed
  integer               :: x_i_idx, y_i_idx, z_i_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: i_bkt, saved_i_bkt_img_j_hi
#ifdef DIRFRC_COMTRANS
  logical               :: common_tran
#endif /* DIRFRC_COMTRANS */
  logical               :: dont_skip_belly_pairs
  integer               :: iaci, ntypes_stk
  integer               :: ee_eval_cnt
  integer               :: full_eval_cnt
  integer               :: is_orthog_stk


! Variables use in the included files.  We would put this stuff in 
! contained subroutines, but it really whacks performance for itanium/ifort.

  double precision      :: f1, f2, f3
  double precision      :: dx, dy, dz
  double precision      :: x_j, y_j, z_j
  integer               :: cur_bkt, yz_bkt, z_bkt
  integer               :: cur_tran, yz_tran, z_tran
  integer               :: atm_j, img_j
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  double precision      :: r_sq
  integer               :: enc_img
  double precision      :: i_tranvec(1:3, 0:17)

! Variables used for double cutoffs:

  double precision      :: es_cut_sq
  logical               :: cutoffs_equal

  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: x_trans(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: y_trans(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
  integer               :: z_trans(0 : cit_tbl_z_dim * 2 - 1)

! integer               :: total_pairs  ! DBG

  integer               :: img_j_ee_eval(atm_cnt)
  integer               :: img_j_full_eval(atm_cnt)

  num_packed = 0        ! For "regular" pairlist
  excl_pairlst_len = 0  ! For "excluded" pairlist

! total_pairs = 0       ! DBG


  dont_skip_belly_pairs = ibelly .eq. 0

  ntypes_stk = ntypes

  is_orthog_stk = is_orthog

  x_box = box(1)
  y_box = box(2)
  z_box = box(3)

  ucell_stk(:,:) = ucell(:,:)

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

  es_cut_sq = (es_cutoff + skinnb) * (es_cutoff + skinnb)

  cutoffs_equal = (es_cutoff .eq. vdw_cutoff)

  cutlist_sq = (vdw_cutoff + skinnb) * (vdw_cutoff + skinnb)

  ! In the following "cut" calcs, we add 0.0001d0 to thwart rounding errors
  ! in the gridding.

  cut_x_frac = (vdw_cutoff + skinnb) / x_box + 0.0001d0
  cut_y_frac = (vdw_cutoff + skinnb) / y_box + 0.0001d0
  cut_z_frac = (vdw_cutoff + skinnb) / z_box + 0.0001d0

  if (is_orthog_stk .eq. 0) then
    cut_x_frac = cut_factor(1) * cut_x_frac
    cut_y_frac = cut_factor(2) * cut_y_frac
    cut_z_frac = cut_factor(3) * cut_z_frac
  end if

  ! Set up the bucket and translation index mapping arrays.  We really only need
  ! to do this each time if the cit table dimensions are changing, but there
  ! may be advantages to having only local variables anyway.  I should check
  ! it out...

  do i = 0, cit_tbl_x_dim - 1
    x_bkts(i) = i
    x_bkts(i + cit_tbl_x_dim) = i
    x_bkts(i + cit_tbl_x_dim + cit_tbl_x_dim) = i
    x_trans(i) = 0
    x_trans(i + cit_tbl_x_dim) = 1
    x_trans(i + cit_tbl_x_dim + cit_tbl_x_dim) = 2
  end do

  do i = 0, cit_tbl_y_dim - 1
    y_bkts(i) = i * cit_tbl_x_dim
    y_bkts(i + cit_tbl_y_dim) = i * cit_tbl_x_dim
    y_bkts(i + cit_tbl_y_dim + cit_tbl_y_dim) = i * cit_tbl_x_dim
    y_trans(i) = 0
    y_trans(i + cit_tbl_y_dim) = 3
    y_trans(i + cit_tbl_y_dim + cit_tbl_y_dim) = 6
  end do

  do i = 0, cit_tbl_z_dim - 1
    z_bkts(i) = i * cit_tbl_x_dim * cit_tbl_y_dim
    z_bkts(i + cit_tbl_z_dim) = i * cit_tbl_x_dim * cit_tbl_y_dim
    z_trans(i) = 0
    z_trans(i + cit_tbl_z_dim) = 9
  end do

  ! Loop over the atoms you own, finding their nonbonded partners...
  ! img_i is the atom for which you are finding the partner atoms, img_j

  do img_i = my_img_lo, my_img_hi       ! main atom loop

    excl_pairlst_len = excl_pairlst_len + 1
    excl_img_j_cnt = 0

    ee_eval_cnt = 0
    full_eval_cnt = 0
    iaci = ntypes_stk * (img_iac(img_i) - 1)

    atm_i = img_atm_map(img_i)

    ! Mark excluded j images:

    atm_mask_idx = atm_maskdata(atm_i)%maskptr
    do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%nummask
      excl_img_flags(atm_img_map(atm_mask(i))) = 1
    end do

    ! These are imaged coordinates:

    x_i = img(img_i)%x
    y_i = img(img_i)%y
    z_i = img(img_i)%z

    ! These are indexes into the 3D cit:

    x_i_idx = int(fraction(1, atm_i) * scale_fac_x)
    y_i_idx = int(fraction(2, atm_i) * scale_fac_y)
    z_i_idx = int(fraction(3, atm_i) * scale_fac_z)

    ! This is the flat cit bucket index:

    i_bkt = x_i_idx + cit_tbl_x_dim * (y_i_idx + z_i_idx * cit_tbl_y_dim)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = int((1.d0 + fraction(1, atm_i) - cut_x_frac) * scale_fac_x)
    x_bkts_hi = int((1.d0 + fraction(1, atm_i) + cut_x_frac) * scale_fac_x)
    y_bkts_lo = int((1.d0 + fraction(2, atm_i) - cut_y_frac) * scale_fac_y)
    y_bkts_hi = int((1.d0 + fraction(2, atm_i) + cut_y_frac) * scale_fac_y)
    z_bkts_lo = int((1.d0 + fraction(3, atm_i) - cut_z_frac) * scale_fac_z)
    z_bkts_hi = z_i_idx + cit_tbl_z_dim
                                        
    saved_i_bkt_img_j_hi = flat_cit(i_bkt)%img_hi
    flat_cit(i_bkt)%img_hi = img_i - 1

#ifdef DIRFRC_COMTRANS
    common_tran = x_trans(x_bkts_lo) .eq. 1 .and. &
                  x_trans(x_bkts_hi) .eq. 1 .and. &
                  y_trans(y_bkts_lo) .eq. 3 .and. &
                  y_trans(y_bkts_hi) .eq. 3 .and. &
                  z_trans(z_bkts_lo) .eq. 9 .and. &
                  z_trans(z_bkts_hi) .eq. 9

    if (cutoffs_equal) then
      if (common_tran) then
#include "find_img_pairs_comtran_1cut.i"
      else
#include "find_img_pairs_nocomtran_1cut.i"
      end if
    else
      if (common_tran) then
#include "find_img_pairs_comtran_2cut.i"
      else
#include "find_img_pairs_nocomtran_2cut.i"
      end if
    end if

#else

    if (cutoffs_equal) then
#include "find_img_pairs_nocomtran_1cut.i"
    else
#include "find_img_pairs_nocomtran_2cut.i"
    end if

#endif /* DIRFRC_COMTRANS */

    flat_cit(i_bkt)%img_hi = saved_i_bkt_img_j_hi

!   total_pairs = total_pairs + img_j_cnt       ! DBG

    if (dont_skip_belly_pairs) then
      call pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                        full_eval_cnt, img_j_full_eval, &
                        gbl_ipairs, num_packed)
    else
      call pack_nb_list_skip_belly_pairs(ee_eval_cnt, img_j_ee_eval, &
                                         full_eval_cnt, img_j_full_eval, &
                                         gbl_ipairs, igroup, img_atm_map, &
                                         num_packed)
    end if

    ! Clear excluded j images flags:

    do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%nummask
      excl_img_flags(atm_img_map(atm_mask(i))) = 0
    end do

    if (ifail .eq. 1) return

    excl_img_pairlst(excl_pairlst_len) = excl_img_j_cnt
    excl_pairlst_len = excl_pairlst_len + excl_img_j_cnt

  end do ! main atom loop

!  write(mdout, *) 'DBG: total_pairs, listtot = ', & 
!               total_pairs, num_packed
!  write(mdout, *) 'DBG: avg packing = ', &
!              dble(total_pairs)/dble(num_packed)

!  if (verbose .gt. 0) then
!     if (master) write(mdout, *) 'listtot = ', num_packed
!  end if

  return

contains

#include "pack_nb_list_dflt.i"

end subroutine get_nb_list

!*******************************************************************************
!
! Subroutine:   get_excl_img_list
!
! Description:  Main routine for list building.  Note it is structured for
!               parallelism.
!
!*******************************************************************************

!#ifdef MPI
subroutine get_excl_img_list(atm_cnt, atm_img_map, img_atm_map, crd, box, &
                             excl_img_pairlst, atm_maskdata, atm_mask, &
                             img_iac, fraction, img, charge, iac)
!#else
!subroutine get_excl_img_list(atm_cnt, atm_img_map, img_atm_map, crd, &
!                             excl_img_pairlst, atm_maskdata, atm_mask)
!#endif

  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: atm_img_map(atm_cnt)
  integer               :: img_atm_map(atm_cnt)
  double precision      :: crd(3, *)
  double precision      :: box(3)
  integer               :: excl_img_pairlst(*)
  type(maskdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
!#ifdef MPI
  integer               :: img_iac(atm_cnt)
  double precision      :: fraction(3, atm_cnt)
  type(img_rec)         :: img(atm_cnt)
  double precision      :: charge(atm_cnt)
  integer               :: iac(atm_cnt)
!#endif /* MPI */

! Local variables:

  integer               :: atm_i, img_i
  integer               :: atm_j, img_j
  integer               :: atm_mask_off, atm_mask_idx
  integer               :: excl_img_j_cnt
  integer               :: excl_pairlst_len
!#ifdef MPI
  integer               :: is_orthog_stk
  double precision      :: x_box, y_box, z_box
  double precision      :: ucell_stk(3, 3)
  double precision      :: x_j, y_j, z_j
  double precision      :: f1, f2, f3
!#endif /* MPI */

  if( numtasks .gt. 1) then
    is_orthog_stk = is_orthog
    x_box = box(1)
    y_box = box(2)
    z_box = box(3)

    ucell_stk(:,:) = ucell(:,:)
  endif

  excl_pairlst_len = 0

  ! Loop over the atoms you own, finding excluded images you should process.

  do img_i = my_img_lo, my_img_hi       ! main atom loop

    excl_pairlst_len = excl_pairlst_len + 1
    excl_img_j_cnt = 0

    atm_i = img_atm_map(img_i)

    atm_mask_off = atm_maskdata(atm_i)%maskptr

    do atm_mask_idx = atm_mask_off + 1, &
                      atm_mask_off + atm_maskdata(atm_i)%nummask

      atm_j = atm_mask(atm_mask_idx)

      if( numtasks .gt. 1) then
        ! Under mpi, a process claims excluded pairs where the real z crd
        ! value for atm_j is less than that of atm_i.  If by some chance equal,
        ! the img id is used as a tie-breaker.

         if (crd(3, atm_j) .lt. crd(3, atm_i)) then
            img_j = atm_img_map(atm_j)
            excl_img_j_cnt = excl_img_j_cnt + 1
            excl_img_pairlst(excl_pairlst_len + excl_img_j_cnt) = img_j
            if (img_atm_map(img_j) .lt. 0) then
                img_atm_map(img_j) = - img_atm_map(img_j)
                if (img_iac(img_j) .eq. 0) then
                    if (is_orthog_stk .ne. 0) then
                        x_j = fraction(1, atm_j) * x_box
                        y_j = fraction(2, atm_j) * y_box
                        z_j = fraction(3, atm_j) * z_box
                    else
                        f1 = fraction(1, atm_j)
                        f2 = fraction(2, atm_j)
                        f3 = fraction(3, atm_j)
                        ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
                        ! we can simplify the expression in this critical inner loop
                        x_j = f1 * ucell_stk(1, 1) + f2 * ucell_stk(1, 2) + &
                        f3 * ucell_stk(1, 3)
                        y_j = f2 * ucell_stk(2, 2) + f3 * ucell_stk(2, 3)
                        z_j = f3 * ucell_stk(3, 3)
                    end if
                    img(img_j)%x = x_j
                    img(img_j)%y = y_j
                    img(img_j)%z = z_j
                    img(img_j)%qterm = charge(atm_j)
                    img_iac(img_j) = iac(atm_j)
                end if
            end if
         else if (crd(3, atm_j) .gt. crd(3, atm_i)) then
            ! The other guy gets this one...
         else
            img_j = atm_img_map(atm_j)
            if (img_j .lt. img_i) then
                excl_img_j_cnt = excl_img_j_cnt + 1
                excl_img_pairlst(excl_pairlst_len + excl_img_j_cnt) = img_j
                if (img_atm_map(img_j) .lt. 0) then
                    img_atm_map(img_j) = - img_atm_map(img_j)
                    if (img_iac(img_j) .eq. 0) then
                        if (is_orthog_stk .ne. 0) then
                            x_j = fraction(1, atm_j) * x_box
                            y_j = fraction(2, atm_j) * y_box
                            z_j = fraction(3, atm_j) * z_box
                        else
                            f1 = fraction(1, atm_j)
                            f2 = fraction(2, atm_j)
                            f3 = fraction(3, atm_j)
                            ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
                            ! we can simplify the expression in this critical inner loop
                            x_j = f1 * ucell_stk(1, 1) + f2 * ucell_stk(1, 2) + &
                            f3 * ucell_stk(1, 3)
                            y_j = f2 * ucell_stk(2, 2) + f3 * ucell_stk(2, 3)
                            z_j = f3 * ucell_stk(3, 3)
                        end if
                        img(img_j)%x = x_j
                        img(img_j)%y = y_j
                        img(img_j)%z = z_j
                        img(img_j)%qterm = charge(atm_j)
                        img_iac(img_j) = iac(atm_j)
                    end if
                end if
            end if
         end if
      else
         ! Single process just maps atms to imgs:
         if (atm_i .lt. atm_j) then
            img_j = atm_img_map(atm_j)
            excl_img_j_cnt = excl_img_j_cnt + 1
            excl_img_pairlst(excl_pairlst_len + excl_img_j_cnt) = img_j
         end if
      endif 

    end do

    excl_img_pairlst(excl_pairlst_len) = excl_img_j_cnt
    excl_pairlst_len = excl_pairlst_len + excl_img_j_cnt

  end do ! main atom loop

  return

end subroutine get_excl_img_list



!*******************************************************************************
!
! Subroutine:  check_my_atom_movement
!
! Description:  This routine checks if any atom on this task's atom list has
!               moved more than half the skin distance.  This is intended for
!               use only with parallel MD.
!
!*******************************************************************************

subroutine check_my_atom_movement(crd, saved_crd, box, saved_box, my_atm_lst, my_atm_cnt, skinnb, ntp, &
                                  new_list)

  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  double precision      :: box(3)
  double precision      :: saved_box(3)
  integer               :: my_atm_lst(*)
  integer               :: my_atm_cnt
  double precision      :: skinnb
  integer               :: ntp
  logical               :: new_list

! Local variables:

  integer               :: atm_lst_idx
  double precision      :: dx, dy, dz
  double precision      :: box_dx, box_dy, box_dz
  double precision      :: distance_limit
  integer               :: n

! If anybody moves more than half the nbskin added to cutoff in generating the
! pairlist, update the pairlist.

  distance_limit = 0.25d0 * skinnb * skinnb

  if (ntp .eq. 0) then  ! Constant volume

    do atm_lst_idx = 1, my_atm_cnt
      n = my_atm_lst(atm_lst_idx)

      dx = crd(1, n) - saved_crd(1, n)
      dy = crd(2, n) - saved_crd(2, n)
      dz = crd(3, n) - saved_crd(3, n)

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  else  ! Constant pressure scaling.

    box_dx = box(1) / saved_box(1)
    box_dy = box(2) / saved_box(2)
    box_dz = box(3) / saved_box(3)

    do atm_lst_idx = 1, my_atm_cnt
      n = my_atm_lst(atm_lst_idx)

      dx = crd(1, n) - saved_crd(1, n) * box_dx
      dy = crd(2, n) - saved_crd(2, n) * box_dy
      dz = crd(3, n) - saved_crd(3, n) * box_dz

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  end if

  new_list = .false.

  return

end subroutine check_my_atom_movement

!*******************************************************************************
!
! Subroutine:  save_used_atom_crds
!
! Description:  This is needed for the skin test for buffered pairlists, but
!               is also used in adjusting image coordinates between list
!               builds, hence the need to look at the send_atm_lst.
!               Used only in MPI
!*******************************************************************************

subroutine save_used_atom_crds(atm_cnt, crd, saved_crd, send_atm_lst, &
                               send_atm_cnts, off_tbl, my_atm_lst, my_atm_cnt)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: saved_crd(3, atm_cnt)
  integer               :: send_atm_lst(atm_cnt)
  integer               :: send_atm_cnts(0: numtasks - 1)
  integer               :: off_tbl(0:numtasks)
  integer               :: my_atm_lst(atm_cnt)
  integer               :: my_atm_cnt
  

! Local variables:

  integer               :: atm_lst_idx, atm_id
  integer               :: i
  integer               :: node, node_offset
  integer               :: send_cnt

    do atm_lst_idx = 1, my_atm_cnt
      atm_id = my_atm_lst(atm_lst_idx)
      saved_crd(1, atm_id) = crd(1, atm_id)
      saved_crd(2, atm_id) = crd(2, atm_id)
      saved_crd(3, atm_id) = crd(3, atm_id)
    end do

    do node = 0, numtasks - 1

      if (node .ne. mytaskid) then
        node_offset = off_tbl(node)
        send_cnt = send_atm_cnts(node)

        ! This also picks up stuff on the extra atoms list, which is probably
        ! not necessary but not harmful.

        do i = 1, send_cnt
          atm_id = send_atm_lst(node_offset + i)
          saved_crd(1, atm_id) = crd(1, atm_id)
          saved_crd(2, atm_id) = crd(2, atm_id)
          saved_crd(3, atm_id) = crd(3, atm_id)
        end do
      end if

    end do

  return

end subroutine save_used_atom_crds

!*******************************************************************************
!
! Subroutine:  check_all_atom_movement
!
! Description:  This routine checks if any atom has moved more than half the
!               skin distance.  This is intended for all minimizations and
!               for single processor MD.
!*******************************************************************************

subroutine check_all_atom_movement(atm_cnt, crd, saved_crd, box, saved_box, skinnb, ntp, new_list)

  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  double precision      :: box(3)
  double precision      :: saved_box(3)
  double precision      :: skinnb
  integer               :: ntp
  logical               :: new_list

! Local variables:

  double precision      :: dx, dy, dz
  double precision      :: box_dx, box_dy, box_dz
  double precision      :: distance_limit
  integer               :: n

! If anybody moves more than half the nbskin added to cutoff in generating the
! pairlist, update the pairlist.

  distance_limit = 0.25d0 * skinnb * skinnb

  if (ntp .eq. 0) then  ! Constant volume

    do n = 1, atm_cnt

      dx = crd(1, n) - saved_crd(1, n)
      dy = crd(2, n) - saved_crd(2, n)
      dz = crd(3, n) - saved_crd(3, n)

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  else  ! Constant pressure scaling.

    box_dx = box(1) / saved_box(1)
    box_dy = box(2) / saved_box(2)
    box_dz = box(3) / saved_box(3)

    do n = 1, atm_cnt

      dx = crd(1, n) - saved_crd(1, n) * box_dx
      dy = crd(2, n) - saved_crd(2, n) * box_dy
      dz = crd(3, n) - saved_crd(3, n) * box_dz

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  end if

  new_list = .false.
  return

end subroutine check_all_atom_movement

!*******************************************************************************
!
! Subroutine:  save_all_atom_crds
!
! Description:  This is needed for the skin test for buffered pairlists.
!*******************************************************************************

subroutine save_all_atom_crds(atm_cnt, crd, saved_crd)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: saved_crd(3, atm_cnt)

! Local variables:

  integer               :: i

  do i = 1, atm_cnt
    saved_crd(1, i) = crd(1, i)
    saved_crd(2, i) = crd(2, i)
    saved_crd(3, i) = crd(3, i)
  end do

  return

end subroutine save_all_atom_crds

!*******************************************************************************
!
! Subroutine:  save_imgcrds
!
! Description:  <TBS> 
!
!*******************************************************************************

subroutine save_imgcrds(img_cnt, img_lo, img_hi, img, saved_imgcrd)

  use img_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: img_cnt, img_lo, img_hi
  type(img_rec)         :: img(*)
  double precision      :: saved_imgcrd(3, *)

! Local variables:

  integer               :: i

  if( numtasks .gt. 1) then

  ! This routine saves some "unused" images which have uninitialized addresses.
  ! This is probably cheaper than checking the img_atm_map to insure that
  ! the image is in use, and is harmless.

    if (used_img_range_wraps) then
        do i = img_lo, img_cnt
            saved_imgcrd(1,i) = img(i)%x
            saved_imgcrd(2,i) = img(i)%y
            saved_imgcrd(3,i) = img(i)%z
        end do
        do i = 1, img_hi
            saved_imgcrd(1,i) = img(i)%x
            saved_imgcrd(2,i) = img(i)%y
            saved_imgcrd(3,i) = img(i)%z
        end do
    else
        do i = img_lo, img_hi
            saved_imgcrd(1,i) = img(i)%x
            saved_imgcrd(2,i) = img(i)%y
            saved_imgcrd(3,i) = img(i)%z
        end do
    end if
  else 
    do i = 1, img_cnt
        saved_imgcrd(1,i) = img(i)%x
        saved_imgcrd(2,i) = img(i)%y
        saved_imgcrd(3,i) = img(i)%z
    end do
  endif

  return

end subroutine save_imgcrds

!*******************************************************************************
!
! Subroutine:  adjust_imgcrds
!
! Description:  <TBS>
!              
!*******************************************************************************

subroutine adjust_imgcrds(img_cnt, img_lo, img_hi, img, img_atm_map, &
                          saved_imgcrd, crd, saved_crd, box, saved_box, ntp)

  use img_mod
  use pbc_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: img_cnt
  integer               :: img_lo, img_hi
  type(img_rec)         :: img(*)
  integer               :: img_atm_map(*)
  double precision      :: saved_imgcrd(3, *)
  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  double precision      :: box(3)
  double precision      :: saved_box(3)
  integer               :: ntp

! Local variables:

  integer               :: i, j
  double precision      :: box_del(3)

  if (ntp .eq. 0) then  ! Constant volume

    if( numtasks .gt. 1) then

        if (used_img_range_wraps) then

            do i = img_lo, img_cnt
                j = img_atm_map(i)
                if (j .gt. 0) then
                    img(i)%x = crd(1,j) + saved_imgcrd(1,i) - saved_crd(1,j)
                    img(i)%y = crd(2,j) + saved_imgcrd(2,i) - saved_crd(2,j)
                    img(i)%z = crd(3,j) + saved_imgcrd(3,i) - saved_crd(3,j)
                endif
            end do

            do i = 1, img_hi
                j = img_atm_map(i)
                if (j .gt. 0) then
                    img(i)%x = crd(1,j) + saved_imgcrd(1,i) - saved_crd(1,j)
                    img(i)%y = crd(2,j) + saved_imgcrd(2,i) - saved_crd(2,j)
                    img(i)%z = crd(3,j) + saved_imgcrd(3,i) - saved_crd(3,j)
                endif
            end do

        else

            do i = img_lo, img_hi
                j = img_atm_map(i)
                if (j .gt. 0) then
                    img(i)%x = crd(1,j) + saved_imgcrd(1,i) - saved_crd(1,j)
                    img(i)%y = crd(2,j) + saved_imgcrd(2,i) - saved_crd(2,j)
                    img(i)%z = crd(3,j) + saved_imgcrd(3,i) - saved_crd(3,j)
                endif
            end do
        end if

    else ! (numtasks == 1)

        do i = 1, img_cnt
            j = img_atm_map(i)
            img(i)%x = crd(1,j) + saved_imgcrd(1,i) - saved_crd(1,j)
            img(i)%y = crd(2,j) + saved_imgcrd(2,i) - saved_crd(2,j)
            img(i)%z = crd(3,j) + saved_imgcrd(3,i) - saved_crd(3,j)
        end do

    endif ! (numtasks > 1)

  else  ! (ntp != 0) Constant pressure scaling.

    box_del(:) = box(:) / saved_box(:)

    if( numtasks .gt. 1) then

        if (used_img_range_wraps) then

            do i = img_lo, img_cnt
                j = img_atm_map(i)
                if (j .gt. 0) then
                    img(i)%x = crd(1,j) + &
                        (saved_imgcrd(1,i) - saved_crd(1,j)) * box_del(1)
                    img(i)%y = crd(2,j) + &
                        (saved_imgcrd(2,i) - saved_crd(2,j)) * box_del(2)
                    img(i)%z = crd(3,j) + &
                        (saved_imgcrd(3,i) - saved_crd(3,j)) * box_del(3)
                endif
            end do

            do i = 1, img_hi
                j = img_atm_map(i)
                if (j .gt. 0) then
                    img(i)%x = crd(1,j) + &
                        (saved_imgcrd(1,i) - saved_crd(1,j)) * box_del(1)
                    img(i)%y = crd(2,j) + &
                        (saved_imgcrd(2,i) - saved_crd(2,j)) * box_del(2)
                    img(i)%z = crd(3,j) + &
                        (saved_imgcrd(3,i) - saved_crd(3,j)) * box_del(3)
                endif
            end do

        else !  ( !used_img_range_wraps )

            do i = img_lo, img_hi
                j = img_atm_map(i)
                if (j .gt. 0) then
                    img(i)%x = crd(1,j) + &
                        (saved_imgcrd(1,i) - saved_crd(1,j)) * box_del(1)
                    img(i)%y = crd(2,j) + &
                        (saved_imgcrd(2,i) - saved_crd(2,j)) * box_del(2)
                    img(i)%z = crd(3,j) + &
                        (saved_imgcrd(3,i) - saved_crd(3,j)) * box_del(3)
                endif
            end do

        end if !  ( used_img_range_wraps )

    else !  (numtasks == 1)

        do i = 1, img_cnt
            j = img_atm_map(i)
            img(i)%x = crd(1,j) + &
                   (saved_imgcrd(1,i) - saved_crd(1,j)) * box_del(1)
            img(i)%y = crd(2,j) + &
                   (saved_imgcrd(2,i) - saved_crd(2,j)) * box_del(2)
            img(i)%z = crd(3,j) + &
                   (saved_imgcrd(3,i) - saved_crd(3,j)) * box_del(3)
        end do

    end if ! ( numtasks > 1)

  end if !  (ntp == 0) 

  return

end subroutine adjust_imgcrds

end module nb_pairlist_mod
