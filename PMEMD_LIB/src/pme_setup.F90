#include "copyright.i"

!*******************************************************************************
!
! Module: pme_setup_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module pme_setup_mod

  use gbl_datatypes_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  final_pme_setup
!
! Description:  Called after dynamic memory allocations to fill in some arrays
!               and perform other initial chores.  Called by ALL processes.
!*******************************************************************************

subroutine final_pme_setup(igroup, tranvec)

  use constraints_mod
  use img_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use alltasks_setup_mod
  use ene_frc_splines_mod
  use prmtop_dat_mod
  use timers_mod

  implicit none

! Formal arguments:

  integer  ::  igroup(*)
  double precision :: tranvec (*)

! Local variables:

  double precision              :: tim1, tim2
  double precision              :: max_erfc_relerr

  ! check_neutral() is currently executed by all processes.  This is not
  ! really necessary but we leave it here for now to keep output consistent
  ! with sander.

  call check_neutral

  call second(tim1)

  call fill_eed_table

  call fill_ef_spline_table

  if (iamoeba .eq. 0) call vdw_correct_setup

!  call load_list_mask(atm_numex, gbl_natex, natom, next + next, igroup)
  call load_list_mask(natom, atm_numex, next, gbl_natex, igroup)

  call fill_tranvec(tranvec)
  
  call chk_switch(max_erfc_relerr)
  call chk_ef_spline_tables(max_erfc_relerr)

  call second(tim2)

  return

end subroutine final_pme_setup

!*******************************************************************************
!
! Subroutine:  check_neutral
!
! Description:
!
! In general, the Ewald method is only truly applicable under conditions of
! charge neutrality.  When the system is not net neutral, the direct sum and
! reciprocal sums are not "beta" independent.  Regardless, the Ewald method can
! be applied with the ficticious assumption that there is a "uniform net
! neutralizing plasma".  
!
! This routine will remove any net charge resulting from conversion of the low
! precision parm topology charges in the case that the system is supposed to
! be net neutral.
!              
!*******************************************************************************

subroutine check_neutral

  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

  integer               i
  double precision      sum

  if (iamoeba .eq. 0) then

      sum = 0.d0

      do i = 1, natom
          sum = sum + atm_qterm(i)
      end do

!      if (master) write(*, '(/,5x,a,f12.8)') 'Sum of charges = ',  sum / 18.2223d0
    
      if (abs(sum/18.2223d0) .gt. 0.01) then
!         if (master) write(*, '(5x,a,/)') 'Assuming uniform neutralizing plasma'
      else
!         if (master) write(*, '(5x,a,/)') 'Forcing neutrality...'
         sum = sum / natom
         do i = 1, natom
            atm_qterm(i) = atm_qterm(i) - sum
         end do
      endif
  else
     ! BUGBUG - Need to do something for amoeba?
  endif

  return

end subroutine check_neutral

!*******************************************************************************
!
! Subroutine:  fill_eed_table
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine fill_eed_table

  use bspline_mod
  use pme_direct_mod
  use mdin_ewald_dat_mod

  implicit none

  double precision      del, x, switch, d_switch_dx

  integer               i, ibcbeg, ibcend

  double precision, dimension(mxeedtab) :: tau

  del = 1.d0 / eedtbdns

  x = 0.d0
  call get_ee_func(x, switch, d_switch_dx)
  gbl_eed_cub(2) = d_switch_dx
  x = (mxeedtab - 1) * del
  call get_ee_func(x, switch, d_switch_dx)
  gbl_eed_cub(2 + 4 * (mxeedtab - 1)) = d_switch_dx

  do i = 1, mxeedtab
    x = del * (i - 1)
    call get_ee_func(x, switch, d_switch_dx)
    tau(i) = x
    gbl_eed_cub(1 + 4 * (i - 1)) = switch
  end do

  ibcbeg = 1
  ibcend = 1
  call cubspl(tau, gbl_eed_cub, mxeedtab, ibcbeg, ibcend)

  return
  
end subroutine fill_eed_table

!*******************************************************************************
!
! Subroutine:  chk_switch
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine chk_switch(max_erfc_relerr)
!
! Compute estimated error for interpolation evaluation of forces and energy ???
!
  use pme_direct_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision, intent(out) :: max_erfc_relerr
  
! Local variables:

  double precision      switch
  double precision      aswitch
  double precision      d_switch_dx
  double precision      ad_switch_dx
  double precision      err
  double precision      max_err
  double precision      max_err_at
  double precision      derr
  double precision      max_derr
  double precision      max_derr_at
  double precision      df_err
  double precision      max_df_err
  double precision      max_df_err_at

  double precision      del
  double precision      dx
  double precision      dxdr
  double precision      half
  integer               i, j
  integer               ind
  integer               max, min
  double precision      small
  double precision      third
  double precision      x
  double precision      df_term
  double precision      adf_term

  half = 0.5d0
  third = 1.d0/3.d0
  dxdr = ew_coeff

!  if (master) then
!    write(mdout, *) &
!     '---------------------------------------------------'
!  end if

  small = 1.d-4
  max = dxdr * es_cutoff * eedtbdns
  min = dxdr * eedtbdns
  min = 1
  max_err = 0.d0
  max_derr = 0.d0
  max_df_err = 0.d0

  max_err_at = - 1.d0
  max_derr_at = - 1.d0
  max_df_err_at = - 1.d0

  del = 1.d0 / eedtbdns

!  if (master) then
!    write(mdout, 1000)
!    write(mdout, 1001) eedtbdns
!    write(mdout, 1002)
!  end if

!1000  format(1x, 'APPROXIMATING switch and d/dx switch', &
!                 ' using CUBIC SPLINE INTERPOLATION')

!1001  format(1x, 'using ', f8.1, ' points per unit in tabled values')

!1002  format(1x, 'TESTING RELATIVE ERROR over r ranging from ', '0.0 to cutoff')

  do i = min, max
    do j = 1, 10
      x = del * (i - 1) + j * 0.1d0 * del
      call get_ee_func(x, switch, d_switch_dx)

      ! cubic spline on switch is all we do:

      ind = eedtbdns * x

      dx = x - ind * del

      ind = 4 * ind

      aswitch = gbl_eed_cub(1 + ind) + dx * (gbl_eed_cub(2 + ind) + &
                dx * (gbl_eed_cub(3 + ind) + &
                dx * gbl_eed_cub(4 + ind) * third)*0.5d0)

      ad_switch_dx = gbl_eed_cub(2 + ind) + dx * (gbl_eed_cub(3 + ind) + &
                     dx * gbl_eed_cub(4 + ind) * half)


      err = abs(aswitch - switch) / (abs(switch) + small)

      ! The relative energy error is the same as err above.

      derr = abs(d_switch_dx - ad_switch_dx) /  (abs(d_switch_dx) + small)

      ! To get the relative force error, you have to consider how df is
      ! calculated using the two terms that have error (aswitch, ad_switch_dx).
      ! You can cancel terms across the difference, so there is no need to
      ! actually calculate df per se.  Note x is dxdr * delr.

      df_term = switch - d_switch_dx * x
      adf_term = aswitch - ad_switch_dx * x
      df_err = abs(df_term - adf_term) / (abs(df_term) + small)

      if (err .gt. max_err) then
        max_err = err
        max_err_at = x
      end if

      if (derr .gt. max_derr) then
        max_derr = derr
        max_derr_at = x
      end if

      if (df_err .gt. max_df_err) then
        max_df_err = df_err
        max_df_err_at = x
      end if

    end do
  end do
  
!  if (master) then
!    write(mdout, 1003)max_err, max_err_at
!    write(mdout, 1004)max_derr, max_derr_at
!    write(mdout, *) '---------------------------------------------------'
!  end if

  if (max_df_err .gt. max_err) then
    max_erfc_relerr = max_df_err
  else
    max_erfc_relerr = max_err
  end if

! write(mdout, *)'DBG: max df rel err =', max_df_err, ' at', max_df_err_at

!1003  format('| CHECK switch(x): max rel err = ', e12.4, '   at ', f10.6)
!1004  format('| CHECK d/dx switch(x): max rel err = ', e12.4, '   at ', f10.6)

  return

end subroutine chk_switch

!*******************************************************************************
!
! Subroutine:  get_ee_func
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine get_ee_func(x, switch, d_switch_dx)

  use gbl_constants_mod

  implicit none

  double precision      x, switch, d_switch_dx

! Get switch function multiplying the Coulomb interaction 1/r.
! r has been converted to x by x = dxdr*r for convenience (erfc for ewald):

  call derfcfun(x, switch)

  d_switch_dx = - 2.d0 * exp(-x * x)/sqrt(PI)

  return

end subroutine get_ee_func

!*******************************************************************************
!
! Subroutine:  vdw_correct_setup
!
! Description: <TBS>
!              
! Sets up the numbers of atoms in each vdw type.  Used for
! analytic pressure, energy, correction to vdw dispersion.
!
!*******************************************************************************

subroutine vdw_correct_setup

  use pme_force_mod
  use prmtop_dat_mod

  implicit none

  integer       n, j

  gbl_nvdwcls(:) = 0

  do n = 1, natom
    j = atm_iac(n)
    gbl_nvdwcls(j) = gbl_nvdwcls(j) + 1
  end do

  return

end subroutine vdw_correct_setup

!*******************************************************************************
!
! Subroutine:  load_list_mask
!
! Description:  This routine takes as input the numex() and natex() arrays from
!               the prmtop.  The prmtop format basically consists of a count of
!               excluded nonbond interactions for each atom (in numex()) and a
!               correlated concatenated "list of sublists" of the actual
!               exclusions in natex() (see the prmtop format doc to understand
!               the input - but the doc is not completely correct :-( ))).
!               The input format is flawed in that it is one-directional (only
!               listing exclusions for atoms with a higher atom number) and also
!               it is not possible to get data for an individual atom without
!               traversing the numex() array.  Here we basically double the list
!               (listing nonbonded interactions in both directions) and also
!               change the list data to allow finding the data for a given atom
!               in one probe (the maskdata_rec struct is used, and we make an
!               atm_nb_maskdata() array of this structure and an atm_nb_mask()
!               array that is equivalent to natex() but that accounts for
!               nonbonded interactions in both directions).
!              
!*******************************************************************************

subroutine load_list_mask(atm_cnt, numex, next, natex, igroup)

  use gbl_constants_mod
  use file_io_dat_mod
  use parallel_dat_mod
  use pme_force_mod, only : atm_nb_maskdata, atm_nb_mask
  use amoeba_adjust_mod, only : init_adjust_weight_idxs
  use mdin_ctrl_dat_mod, only : ibelly, iamoeba

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: numex(atm_cnt) ! Excluded atom count for each atom.
  integer               :: next           ! natex() size, from prmtop.
  integer               :: natex(next)    ! Excluded atom concatenated list.
  integer               :: igroup(*) 

! Local variables:

  integer               :: atm_i, atm_j
  integer               :: lst_idx, sublst_idx, num_sublst
  integer               :: mask_idx
  integer               :: offset
  integer               :: total_excl

! Double the mask to deal with our list generator

! Pass 1: get pointers, check size

  lst_idx = 0

  atm_nb_maskdata(:)%nummask = 0       ! array assignment

  do atm_i = 1, atm_cnt - 1         ! last atom never has any...
    num_sublst = numex(atm_i)
    do sublst_idx = 1, num_sublst
      atm_j = natex(lst_idx + sublst_idx)
      if (atm_j .gt. 0) then
        atm_nb_maskdata(atm_i)%nummask = atm_nb_maskdata(atm_i)%nummask + 1
        atm_nb_maskdata(atm_j)%nummask = atm_nb_maskdata(atm_j)%nummask + 1
      end if
    end do
    lst_idx = lst_idx + num_sublst
  end do

  total_excl = 0

  do atm_i = 1, atm_cnt
    total_excl = total_excl + atm_nb_maskdata(atm_i)%nummask
  end do

  if (total_excl .gt. 2 * next) then
    error_msg = 'Variable next has inconsistent value in prmtop!'
    call mol_mech_error
  end if

  offset = 0

  do atm_i = 1, atm_cnt
    atm_nb_maskdata(atm_i)%maskptr = offset
    offset = offset + atm_nb_maskdata(atm_i)%nummask
  end do

! Pass 2: fill mask array

  lst_idx = 0

  atm_nb_maskdata(:)%nummask = 0       ! array assignment
  
  if (ibelly .eq. 0) then
    do atm_i = 1, atm_cnt - 1
      num_sublst = numex(atm_i)
      do sublst_idx = 1, num_sublst
        atm_j = natex(lst_idx + sublst_idx)
        if (atm_j .gt. 0) then

          atm_nb_maskdata(atm_j)%nummask = atm_nb_maskdata(atm_j)%nummask + 1
          mask_idx = atm_nb_maskdata(atm_j)%maskptr + &
                     atm_nb_maskdata(atm_j)%nummask
          atm_nb_mask(mask_idx) = atm_i

          atm_nb_maskdata(atm_i)%nummask = atm_nb_maskdata(atm_i)%nummask + 1
          mask_idx = atm_nb_maskdata(atm_i)%maskptr + &
                     atm_nb_maskdata(atm_i)%nummask
          atm_nb_mask(mask_idx) = atm_j

        end if
      end do
      lst_idx = lst_idx + num_sublst
    end do
  else
    do atm_i = 1, atm_cnt - 1
      num_sublst = numex(atm_i)
      do sublst_idx = 1, num_sublst
        atm_j = natex(lst_idx + sublst_idx)
        if (atm_j .gt. 0) then
          if (igroup(atm_i) .ne. 0 .and. igroup(atm_j) .ne. 0) then

            atm_nb_maskdata(atm_j)%nummask = atm_nb_maskdata(atm_j)%nummask + 1
            mask_idx = atm_nb_maskdata(atm_j)%maskptr + &
                       atm_nb_maskdata(atm_j)%nummask
            atm_nb_mask(mask_idx) = atm_i

            atm_nb_maskdata(atm_i)%nummask = atm_nb_maskdata(atm_i)%nummask + 1
            mask_idx = atm_nb_maskdata(atm_i)%maskptr + &
            atm_nb_maskdata(atm_i)%nummask
            atm_nb_mask(mask_idx) = atm_j

          end if
        end if
      end do
      lst_idx = lst_idx + num_sublst
    end do
  end if

  if (iamoeba .ne. 0) call init_adjust_weight_idxs()
  
  return

end subroutine load_list_mask

!*******************************************************************************
!
! Subroutine:   open_pmemd_output_files
!
! Description:  Routine to open the dumping and restart files.
!*******************************************************************************

end module pme_setup_mod

subroutine start_mdout_log(box)
	
  use file_io_mod
  use prmtop_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use cit_mod
  use pbc_mod

  implicit none
  
  ! Formal arguments
  
  double precision      :: box(3)

  character(8)          :: date
  character(10)         :: time
  integer               :: itmp(3)
  double precision      :: rtmp(3)
  
  call amopen(mdout, mdout_name, 'U', 'F', 'W')
  
  call date_and_time(DATE=date, TIME=time) 

  write(mdout,'(12(a),/)') '| Run on ', date(5:6), '/', date(7:8), '/',  &
        date(1:4), ' at ', time(1:2), ':', time(3:4), ':', time(5:6)

  if (using_pme_potential) then

    itmp(1) = cit_tbl_x_dim
    itmp(2) = cit_tbl_y_dim
    itmp(3) = cit_tbl_z_dim

    write(mdout, '(a, 3i5)')'| Coordinate Index Table dimensions: ', &
                            itmp(1), itmp(2), itmp(3)

    rtmp(1) = box(1)/cit_tbl_x_dim
    rtmp(2) = box(2)/cit_tbl_y_dim
    rtmp(3) = box(3)/cit_tbl_z_dim

    write(mdout,'(a, 3f10.4, /)')'| Direct force subcell size = ', &
                                 rtmp(1), rtmp(2), rtmp(3)

  end if
  
  return

end subroutine start_mdout_log


