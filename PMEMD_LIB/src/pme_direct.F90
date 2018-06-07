#include "copyright.i"

!*******************************************************************************
!
! Module:  pme_direct_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module pme_direct_mod

  implicit none

  ! The following storage is per-process common; ie., it SHOULD be
  ! broadcast from the master to the other processes!

  integer, parameter    :: pme_direct_int_cnt = 1

  integer                  mxeedtab

  common / pme_direct_int / mxeedtab

  save  :: / pme_direct_int /

  double precision, allocatable, save   :: gbl_eed_cub(:)

contains

!*******************************************************************************
!
! Subroutine:  init_pme_direct_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine init_pme_direct_dat()

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use prmtop_dat_mod
  use parallel_dat_mod

  implicit none

! Local variables:

  double precision      :: dxdr
  integer               :: i

! Setup dxdr map from r to x in table lookup of eed

  dxdr = ew_coeff
  
! For eed table, assume all nonbond distances are less than 1.5 * es_cutoff
! between nonbond updates; i.e. no excess motion; this is enforced by
! es_cutoff.

  mxeedtab = int(dxdr * eedtbdns * es_cutoff * 1.5d0)

  call alloc_pme_direct_mem()

  return

end subroutine init_pme_direct_dat

!*******************************************************************************
!
! Subroutine:  alloc_pme_direct_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_pme_direct_mem()

  implicit none

! Local variables:

  integer               :: alloc_failed

  if( allocated(gbl_eed_cub)) deallocate(gbl_eed_cub)

  allocate(gbl_eed_cub(4 * mxeedtab), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  gbl_eed_cub(:) = 0.d0

  return

end subroutine alloc_pme_direct_mem

!*******************************************************************************
!
! Subroutine:  bcast_pme_direct_dat
!
! Description: <TBS>
!              Used only in MPI              
!*******************************************************************************

subroutine bcast_pme_direct_dat

  use parallel_dat_mod

  implicit none

! Local variables:

  call mpi_bcast(mxeedtab, pme_direct_int_cnt, mpi_integer, 0, &
                 lib_mpi_comm, err_code_mpi)

  if (.not. master) then
    call alloc_pme_direct_mem()
  end if

  ! The allocated data is not initialized from the master node.

  return

end subroutine bcast_pme_direct_dat

!*******************************************************************************
!
! Subroutine:  pme_list
!
! Description:  Handles set-up and error checking for calling of get_nb_list
!               which creates the nonbond list.
!
!*******************************************************************************

subroutine pme_list(atm_cnt, crd, saved_crd, box, saved_box, igroup, excl_img_pairlst, atm_maskdata, atm_mask, &
                    atm_owner_map, my_atm_lst, my_atm_cnt, tranvec, imin_par  )

  use cit_mod
  use constraints_mod
  use pbc_mod
  use nb_pairlist_mod
  use pme_recip_mod
  use img_mod
  use loadbal_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use prmtop_dat_mod
  use timers_mod
  use amoeba_recip_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: saved_crd(3, atm_cnt)
  double precision      :: box(3)
  double precision      :: saved_box(3)
  integer               :: igroup(*)
  integer               :: imin_par
  integer               :: excl_img_pairlst(*)
  type(maskdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
  integer               :: atm_owner_map(atm_cnt)
  integer               :: my_atm_lst(atm_cnt) 
  integer               :: my_atm_cnt
  double precision      :: tranvec(*)

! Local variables:

  integer               :: ifail
  integer               :: alloc_failed
  logical               :: dont_skip_belly_pairs
  
  integer               :: i

  ! We 0-base the following array for efficiency, but don't use atm_lst(0) so
  ! we can use 0 as the distinguished non-entry (NULL).

  double precision      :: fraction(3, atm_cnt) ! in range 0.0 - +0.999...
  type(cit_tbl_rec)     :: crd_idx_tbl(0 : cit_tbl_x_dim - 1, &
                                       0 : cit_tbl_y_dim - 1, &
                                       0 : cit_tbl_z_dim - 1)

  ! Under amber 7, when both members of an atom pair are in the belly (ie.,
  ! fixed), we don't add them to the nonbonded pairlist.  This saves time but
  ! makes the energies look different (in a manner that actually does not
  ! affect the simulation). Amber 6 only ignores belly atoms for the "bonded"
  ! forces.  We replicate both behaviours here as appropriate.

  dont_skip_belly_pairs = ibelly .eq. 0
   
  saved_box(:) = box(:)     ! Needed when pressure scaling.
  
  if(numtasks .eq. 1) then
      call save_all_atom_crds(atm_cnt, crd, saved_crd)
  endif
 
  call get_fract_crds(atm_cnt, crd, fraction)

  if( numtasks .gt. 1 .and. iamoeba .eq. 1 .and. master ) then
     call amoeba_setup_cit_master(atm_cnt, box, fraction, crd_idx_tbl, gbl_img, &
                                   gbl_atm_img_map, gbl_img_atm_map, atm_iac, &
                                   gbl_img_iac, atm_qterm)
  else
    call setup_cit(atm_cnt, box, fraction, crd_idx_tbl, gbl_img, gbl_atm_img_map, &
                   gbl_img_atm_map, atm_iac, gbl_img_iac, atm_qterm)
  endif

  call update_pme_time(cit_setup_timer)

  do

    ifail = 0
    
    call get_nb_list(atm_cnt, box, crd_idx_tbl, gbl_img, &
                     gbl_img_atm_map, excl_img_pairlst, &
                     typ_ico, fraction, tranvec, atm_maskdata, atm_mask, &
                     gbl_atm_img_map, gbl_excl_img_flags, atm_iac, &
                     gbl_img_iac, atm_qterm, igroup, ntypes, ibelly, &
                     es_cutoff, vdw_cutoff, skinnb, verbose_pme, ifail)

    if (ifail .eq. 0) exit
        
    ! Deallocate the old array, grow it by 10%, reallocate,
    ! and go back up to the top to try again.

    deallocate(gbl_ipairs)

    ipairs_size = 11 * ipairs_size / 10       

	if( allocated(gbl_ipairs)) deallocate(gbl_ipairs)

    allocate(gbl_ipairs(ipairs_size), stat = alloc_failed)

    if (alloc_failed .ne. 0) then
      error_msg = 'FATAL dynamic memory allocation error in subroutine pme_list'
      call mol_mech_error
    end if

#ifdef PAIRLST_DBG
    write(mdout, '(/,a,5x,a,i8)') &
          '|', 'Nonbonded Pairs Reallocation:', ipairs_size
#endif

  end do
  
  call get_excl_img_list(atm_cnt, gbl_atm_img_map, gbl_img_atm_map, crd, box, &
                         excl_img_pairlst, atm_maskdata, atm_mask, &
                         gbl_img_iac, fraction, gbl_img, atm_qterm, atm_iac)

     
  if( numtasks .gt. 1) then
      if (is_orthog .ne. 0) then
          call claim_recip_imgs(atm_cnt, fraction, box, crd_idx_tbl, gbl_img, &
                                gbl_img_iac, gbl_img_atm_map)
      else
          call claim_recip_imgs_nonorthog(atm_cnt, fraction, crd_idx_tbl, gbl_img, &
                                gbl_img_iac, gbl_img_atm_map)
      end if
    
      ! Anything that claims images must be placed in front of this call. So
      ! far this includes:
      ! setup_cit()
      ! get_nb_list()
      ! claim_recip_imgs()

      call find_img_range(atm_cnt, gbl_img_atm_map)

      ! write(0,*)'DBG: task, img_lo, img_hi, range_wraps= ', &
      !           mytaskid, used_img_lo, used_img_hi, used_img_range_wraps


      ! BUGBUG - For amoeba, we temporarily don't use the used image-based atom
      !          sendlist.  The coordinate and force distribution scheme based on
      !          this is disabled until we reliably determine atom/image use.

      if (iamoeba .eq. 0) then      ! BUGBUG - Temporary conditional!
          call get_send_atm_lst(atm_cnt, gbl_img_atm_map, gbl_atm_img_map, &
                                atm_owner_map, gbl_atm_offsets, &
                                gbl_send_atm_lst, gbl_send_atm_cnts, &
                                used_img_lo, used_img_hi, used_img_range_wraps)
      end if        ! BUGBUG - Temporary conditional!

      ! BUGBUG - For amoeba, we temporarily save all atom coords because we don't
      !          yet know which ones we really use.

      if (iamoeba .eq. 0) then      ! BUGBUG - Temporary conditional!
          if (imin_par .eq. 0) then
              call save_used_atom_crds(atm_cnt, crd, saved_crd, &
                                       gbl_send_atm_lst, gbl_send_atm_cnts, &
                                       gbl_atm_offsets,my_atm_lst, my_atm_cnt)
          else
      ! Minimizations require all coords for skin check...
              call save_all_atom_crds(atm_cnt, crd, saved_crd)
          end if
      else
          call save_all_atom_crds(atm_cnt, crd, saved_crd)
      end if        ! BUGBUG - Temporary conditional!
      
  end if ! ( numtasks > 1)  

  call update_pme_time(build_list_timer)

  return

end subroutine pme_list

!*******************************************************************************
!
! Subroutine:  get_nb_energy
!
! Description:
!              
! The main routine for non bond energy (vdw and hbond) as well as direct part
! of ewald sum.  It is structured for parallelism.
!
!*******************************************************************************

subroutine get_nb_energy(img_frc, img, eed_cub, ipairs, tranvec, &
                         eed, evdw, ehb, eedvir, virial)

  use img_mod
  use timers_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use ene_frc_splines_mod

  implicit none

! Formal arguments:

  double precision              :: img_frc(3, *)
  type(img_rec), intent(in)     :: img(*)
  double precision, intent(in)  :: eed_cub(*)
#ifdef DIRFRC_NOVEC
  integer                       :: ipairs(*)
#else
  integer, intent(in)           :: ipairs(*)
#endif
  double precision, intent(in)  :: tranvec(1:3, 0:17)
  double precision              :: eed
  double precision              :: evdw
  double precision              :: ehb
  double precision              :: eedvir
  double precision              :: virial(3, 3)

! Local variables and parameters:

  double precision      del
  double precision      dxdr
  double precision      eedtbdns_stk
  double precision      eedvir_stk, eed_stk, evdw_stk, ehb_stk
  double precision      max_nb_cut2, es_cut2, es_cut
  double precision      x_i, y_i, z_i
  double precision      x_tran(1:3, 0:17)
  double precision      vxx, vxy, vxz, vyy, vyz, vzz
  integer               i
  integer               ipairs_idx
  integer               ntypes_stk
  integer               img_i
  integer               ee_eval_cnt
  integer               full_eval_cnt
#ifdef DIRFRC_COMTRANS
  integer               common_tran    ! flag - 1 if translation not needed
#endif /* DIRFRC_COMTRANS */
  logical               cutoffs_equal

  ntypes_stk = ntypes

  eedvir_stk = 0.d0
  eed_stk = 0.d0
  evdw_stk = 0.d0
  ehb_stk = 0.d0

  dxdr = ew_coeff
  eedtbdns_stk = eedtbdns
  del = 1.d0 / eedtbdns_stk
  max_nb_cut2 = vdw_cutoff * vdw_cutoff
  es_cut = es_cutoff
  es_cut2 = es_cut * es_cut
  cutoffs_equal = (vdw_cutoff .eq. es_cutoff)

  vxx = 0.d0
  vxy = 0.d0
  vxz = 0.d0
  vyy = 0.d0
  vyz = 0.d0
  vzz = 0.d0

  ipairs_idx = 1

  do img_i = my_img_lo, my_img_hi

#ifdef DIRFRC_COMTRANS
    ! Common translation (ie. no translation) flag is packed at
    ! the front of each sublist followed by the count(s) of sublist
    ! image pair entries.

    common_tran = ipairs(ipairs_idx)
    ipairs_idx = ipairs_idx + 1
#endif /* DIRFRC_COMTRANS */
    
    ! Electrostatic evaluation-only count followed by
    ! full evaluation count packed at the front of each pair sublist.
    
    ee_eval_cnt = ipairs(ipairs_idx)
!   write(0,*)'DBG: ee_eval_cnt =', ee_eval_cnt              
    full_eval_cnt = ipairs(ipairs_idx + 1)
!   write(0,*)'DBG: full_eval_cnt =', full_eval_cnt
    ipairs_idx = ipairs_idx + 2
    
    if (ee_eval_cnt + full_eval_cnt .gt. 0) then

      x_i = img(img_i)%x
      y_i = img(img_i)%y
      z_i = img(img_i)%z

#ifdef DIRFRC_COMTRANS
      if (common_tran .eq. 0) then
#endif /* DIRFRC_COMTRANS */
        ! We need all the translation vectors:
        do i = 0, 17
          x_tran(1, i) = tranvec(1, i) - x_i
          x_tran(2, i) = tranvec(2, i) - y_i
          x_tran(3, i) = tranvec(3, i) - z_i
        end do
#ifdef DIRFRC_COMTRANS
      else
        ! Just put the x,y,z values in the middle cell
        x_tran(1, 13) = - x_i
        x_tran(2, 13) = - y_i
        x_tran(3, 13) = - z_i
      end if
#endif /* DIRFRC_COMTRANS */

      if (cutoffs_equal) then
#ifdef DIRFRC_NOVEC
        call short_ene_novec(img_frc, img, efs_tbl, eed_cub, typ_ico, &
                             ipairs(ipairs_idx), gbl_img_iac, 
                             gbl_cn1, gbl_cn2, x_tran)
#else
        call short_ene_vec(img_frc, img, efs_tbl, eed_cub, typ_ico, &
                           ipairs(ipairs_idx), gbl_img_iac, &
                           gbl_cn1, gbl_cn2, x_tran)
#endif /* DIRFRC_NOVEC */
      else
#ifdef DIRFRC_NOVEC
        call short_ene_novec_2cut(img_frc, img, efs_tbl, eed_cub, typ_ico, &
                                  ipairs(ipairs_idx), gbl_img_iac, &
                                  gbl_cn1, gbl_cn2, x_tran)
#else
        call short_ene_vec_2cut(img_frc, img, efs_tbl, eed_cub, typ_ico, &
                                ipairs(ipairs_idx), gbl_img_iac, &
                                gbl_cn1, gbl_cn2, x_tran)
#endif /* DIRFRC_NOVEC */
      end if

      ipairs_idx = ipairs_idx + ee_eval_cnt + full_eval_cnt

    end if
  end do

  ! Save the energies:
  
  eedvir = eedvir_stk
  eed = eed_stk
  evdw = evdw_stk
  ehb = ehb_stk
  
  ! Save the virials.

  virial(1, 1) = vxx
  virial(1, 2) = vxy
  virial(2, 1) = vxy
  virial(1, 3) = vxz
  virial(3, 1) = vxz
  virial(2, 2) = vyy
  virial(2, 3) = vyz
  virial(3, 2) = vyz
  virial(3, 3) = vzz

  return

contains

#define BUILD_SHORT_ENE_VEC_2CUT
#include "short_ene_vec.i"
#undef BUILD_SHORT_ENE_VEC_2CUT
#include "short_ene_vec.i"
#define BUILD_SHORT_ENE_NOVEC_2CUT
#include "short_ene_novec.i"
#undef BUILD_SHORT_ENE_NOVEC_2CUT
#include "short_ene_novec.i"

end subroutine get_nb_energy

end module pme_direct_mod
