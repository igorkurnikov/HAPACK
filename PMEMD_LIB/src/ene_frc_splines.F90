#include "copyright.i"

!*******************************************************************************
!
! Module:  ene_frc_splines_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module ene_frc_splines_mod

  implicit none

  ! The following storage is per-process common; ie., it SHOULD be
  ! broadcast from the master to the other processes!

  integer, parameter    :: ene_frc_splines_int_cnt = 1

  ! The globally visible spline table will have two merged spline tables, each
  ! with four spline coefficients per tau value.  The spline table size here is
  ! essentially the dimension in tau for one table.

  integer                       efs_tbl_siz

  common / ene_frc_splines_int / efs_tbl_siz

  save  :: / ene_frc_splines_int /

  integer, parameter    :: ene_frc_splines_dbl_cnt = 1

  double precision              efs_tbl_dens

  common / ene_frc_splines_dbl / efs_tbl_dens

  save  :: / ene_frc_splines_dbl /

  double precision, parameter, private  :: small = 1.d-24

  double precision, parameter           :: r_zero_approx = 1.d-6

  double precision, save                :: lowest_efs_delr2

  double precision, allocatable, save   :: efs_tbl(:)

contains

!*******************************************************************************
!
! Subroutine:  init_ene_frc_splines_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine init_ene_frc_splines_dat()

  use mdin_ctrl_dat_mod

  implicit none

! Local variables:

  ! The table is extended by 1 to insure that table end errors don't effect
  ! the calculations.  All actual spline usage should be within
  ! es_cutoff * es_cutoff.

  efs_tbl_dens = 50.d0  ! optimal for 8.0 angstrom es_cutoff

  efs_tbl_siz = int(efs_tbl_dens) * &
                      (ceiling(es_cutoff * es_cutoff) + 1)

! write(6,*)'DBG: efs_tbl_dens = ', efs_tbl_dens
! write(6,*)'DBG: efs_tbl_siz = ', efs_tbl_siz

  call alloc_ene_frc_splines_mem()

  return

end subroutine init_ene_frc_splines_dat

!*******************************************************************************
!
! Subroutine:  alloc_ene_frc_splines_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_ene_frc_splines_mem()

  implicit none

! Local variables:

  integer               :: alloc_failed

  ! The efs_tbl contains two merged tables, each with 4 spline
  ! coefficients per tau value.

  if(allocated(efs_tbl)) deallocate(efs_tbl)
  allocate(efs_tbl(2 * 4 * efs_tbl_siz), stat = alloc_failed)
  if (alloc_failed .ne. 0) call setup_alloc_error
  efs_tbl(:) = 0.d0

  return

end subroutine alloc_ene_frc_splines_mem

!*******************************************************************************
!
! Subroutine:  bcast_ene_frc_splines_dat
!
! Description: <TBS>
!              Used only in MPI              
!*******************************************************************************

subroutine bcast_ene_frc_splines_dat

  use parallel_dat_mod

  implicit none

  call mpi_bcast(efs_tbl_siz, ene_frc_splines_int_cnt, mpi_integer, 0, &
                 lib_mpi_comm, err_code_mpi)

  call mpi_bcast(efs_tbl_dens, ene_frc_splines_dbl_cnt, mpi_double_precision, &
                 0, lib_mpi_comm, err_code_mpi)

  if (.not. master) then
    call alloc_ene_frc_splines_mem()
  end if

  ! The allocated data is not initialized from the master node.

  return

end subroutine bcast_ene_frc_splines_dat

!*******************************************************************************
!
! Subroutine:  fill_ef_spline_table
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine fill_ef_spline_table

  use bspline_mod
  use mdin_ewald_dat_mod

  implicit none

  double precision              :: beta
  double precision              :: del, u, ene, dedu, df, ddfdu
  integer                       :: i, j, k, ibcbeg, ibcend
  integer                       :: spline_tbl_siz
  double precision, parameter   :: half = 1.d0/2.d0
  double precision, parameter   :: sixth = 1.d0/6.d0

  double precision              :: tau(efs_tbl_siz)
  double precision              :: ene_spline_tbl(efs_tbl_siz * 4)
  double precision              :: frc_spline_tbl(efs_tbl_siz * 4)

  beta = ew_coeff
  del = 1.d0 / efs_tbl_dens

  ene_spline_tbl(:) = 0.d0
  frc_spline_tbl(:) = 0.d0

  call fill_ene_spline_table
  call fill_frc_spline_table

  ! Now merge the tables into the globally available table; also fix up the
  ! coefficients to avoid calcs in short_ene_efs

  do i = 1, efs_tbl_siz
    j = 8 * (i - 1)
    k = 4 * (i - 1)
    efs_tbl(1 + j) = ene_spline_tbl(1 + k)
    efs_tbl(2 + j) = ene_spline_tbl(2 + k)
    efs_tbl(3 + j) = ene_spline_tbl(3 + k) * half
    efs_tbl(4 + j) = ene_spline_tbl(4 + k) * sixth
    efs_tbl(5 + j) = frc_spline_tbl(1 + k)
    efs_tbl(6 + j) = frc_spline_tbl(2 + k)
    efs_tbl(7 + j) = frc_spline_tbl(3 + k) * half
    efs_tbl(8 + j) = frc_spline_tbl(4 + k) * sixth
  end do

  return

contains

!*******************************************************************************
!
! Internal Subroutine:  fill_ene_spline_table
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine fill_ene_spline_table

  implicit none

  ! We must estimate the energy when u ~= 0.  We do this at u = 1.d-12, or
  ! r = 1.d-6 angstrom which will produce a manageable but huge number.
  ! Basically, anything anywhere near this close together will cause the
  ! simulation to blow up, as we are basically dealing with a a function with r
  ! in the denominator, where r = sqrt(u).  Also, the energy errors
  ! very near 0 (say 0.01 angstrom**2) are actually larger if you make u
  ! smaller; this table is really not intended for use below 0.5 * 0.5 A**2.

  u = r_zero_approx * r_zero_approx
  call e_of_u(beta, u, ene)
  tau(1) = u
  ene_spline_tbl(1) = ene

  ! We don't set a low end derivative, due to the singularity at 0. This
  ! produces better results on the low end of the table.

  ibcbeg = 0    ! Don't use deriv at low end of table

  ! Set the high end derivative.

  u = (efs_tbl_siz - 1) * del
  call dedu_of_u(beta, u, dedu)
  ene_spline_tbl(2 + 4 * (efs_tbl_siz - 1)) = dedu

  ibcend = 1

  do i = 2, efs_tbl_siz
    u = del * (i - 1)
    call e_of_u(beta, u, ene)
    tau(i) = u
    ene_spline_tbl(1 + 4 * (i - 1)) = ene
  end do

  call cubspl(tau, ene_spline_tbl, efs_tbl_siz, ibcbeg, ibcend)

  return
  
end subroutine fill_ene_spline_table


!*******************************************************************************
!
! Internal Subroutine:  fill_frc_spline_table
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine fill_frc_spline_table

  implicit none

  ! We must estimate the force when u ~= 0.  We do this at u = 1.d-12, or
  ! r = 1.d-6 angstrom which will produce a manageable but huge number.
  ! Basically, anything anywhere near this close together will cause the
  ! simulation to blow up, as we are basically dealing with a a function with r
  ! in the denominator, where r = sqrt(u).

  u = r_zero_approx * r_zero_approx
  call df_of_u(beta, u, df)
  tau(1) = u
  frc_spline_tbl(1) = df

  ! We don't set a low end derivative, due to the singularity at 0. This
  ! produces better results on the low end of the table.

  ibcbeg = 0    ! Don't use deriv at low end of table

  ! Set the high end derivative.

  u = (efs_tbl_siz - 1) * del
  call ddfdu_of_u(beta, u, ddfdu)
  frc_spline_tbl(2 + 4 * (efs_tbl_siz - 1)) = ddfdu

  ibcend = 1

  do i = 2, efs_tbl_siz
    u = del * (i - 1)
    call df_of_u(beta, u, df)
    tau(i) = u
    frc_spline_tbl(1 + 4 * (i - 1)) = df
  end do

  call cubspl(tau, frc_spline_tbl, efs_tbl_siz, ibcbeg, ibcend)

  return
  
end subroutine fill_frc_spline_table

end subroutine fill_ef_spline_table

!*******************************************************************************
!
! Subroutine:  chk_ef_spline_tables
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine chk_ef_spline_tables(max_relerr)

  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision, intent(in) :: max_relerr

! Local variables:

  double precision      :: beta
  double precision      :: del
  double precision      :: ene_err_lim_exceeded_at
  double precision      :: frc_err_lim_exceeded_at
  
  if( .not. dirfrc_efs ) then 
	return
  end if

  beta = ew_coeff
  del = 1.d0 / efs_tbl_dens

!  if (master) then
!    write(mdout, 1000)
!    write(mdout, '(a)') &
!      '| APPROXIMATING direct energy using CUBIC SPLINE INTERPOLATION'
!    write(mdout, 1001) efs_tbl_dens
!  end if

  call chk_ene_spline(efs_tbl_siz, del, beta, efs_tbl, max_relerr, &
                      ene_err_lim_exceeded_at)

!  if (master) then
!    if (ene_err_lim_exceeded_at .ge. 0.d0) then
!      write(mdout, 1002) sqrt(ene_err_lim_exceeded_at)
!    else
!      write(mdout, 1003)
!    end if
!    write(mdout, '(a)') &
!      '| APPROXIMATING direct force using CUBIC SPLINE INTERPOLATION'
!    write(mdout, 1001) efs_tbl_dens
!  end if

  call chk_frc_spline(efs_tbl_siz, del, beta, efs_tbl, max_relerr, &
                      frc_err_lim_exceeded_at)

!  if (master) then
!    if (ene_err_lim_exceeded_at .ge. 0.d0) then
!      write(mdout, 1002) sqrt(frc_err_lim_exceeded_at)
!    else
!      write(mdout, 1003)
!    end if
!    write(mdout, 1000)
!  end if

!1000 format('|---------------------------------------------------')
!1001 format('|  with ', f6.1, ' points per unit in tabled values')
!1002 format('| Relative Error Limit not exceeded for r .gt. ', f6.2)
!1003 format('| Relative Error Limit not exceeded for r .gt. 0.d0') ! unlikely!

  ! store bottom usable end of ene-frc spline table

  if (ene_err_lim_exceeded_at .ge. 0.d0 .and. frc_err_lim_exceeded_at .ge. 0.d0) then
    lowest_efs_delr2 = max(ene_err_lim_exceeded_at, frc_err_lim_exceeded_at) + del
  else
    lowest_efs_delr2 = 0.d0
  end if
  
! if (master) write(6,*)'DBG: lowest_efs_delr2 = ', lowest_efs_delr2

  return

end subroutine chk_ef_spline_tables

!*******************************************************************************
!
! Subroutine:  chk_ene_spline
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine chk_ene_spline(tbl_siz, del, beta, ene_tbl, max_relerr, &
                          err_lim_exceeded_at)

  use mdin_ctrl_dat_mod

  implicit none

! Formal arguments:

  integer,          intent(in)  :: tbl_siz      ! logical table entries
  double precision, intent(in)  :: del          ! u increment per table entry
  double precision, intent(in)  :: beta         ! ewald coefficient
  double precision, intent(in)  :: ene_tbl(tbl_siz * 8)
  double precision, intent(in)  :: max_relerr
  double precision, intent(out) :: err_lim_exceeded_at

  ! err_lim_exceeded_at is the highest value of u with a relative error
  ! greater than max_relerr; if none, then it is -1.d0

! Local variables:

  double precision      :: u, du                ! u is synonym for delr2
                                                ! du applies to spline tbl
  double precision      :: ene, spline_ene
  double precision      :: err
  integer               :: i, j
  integer               :: ind
  integer               :: min_tau_idx, max_tau_idx

  min_tau_idx = 1 
  max_tau_idx = ceiling(es_cutoff * es_cutoff) * efs_tbl_dens

  ! Sanity check.
  if (max_tau_idx .gt. tbl_siz) max_tau_idx = tbl_siz

  err_lim_exceeded_at = - 1.d0

  do i = min_tau_idx, max_tau_idx
    ind = ishft((i - 1), 3)       ! 8 * (i - 1)
    do j = 1, 9
      u = dble(i - 1) * del + dble(j) * 0.1d0 * del
      du = dble(j) * 0.1d0 * del
      call e_of_u(beta, u, ene)

      spline_ene = (ene_tbl(1 + ind) + du * &
                   (ene_tbl(2 + ind) + du * &
                   (ene_tbl(3 + ind) + du * &
                    ene_tbl(4 + ind))))

      err = abs(spline_ene - ene) / (abs(ene) + small)
      if (err .gt. max_relerr) err_lim_exceeded_at = u
#ifdef EFS_DBG
      write(0,9001) err, u
!     write(0,9002) i,j,ind,du
      write(0,9003) ene, spline_ene, u
9001  format('| DBG_CHECK ene(r**2): rel err = ', e12.3, '   at ', f10.6)
9002  format('| i =', i6, ' j =', i6, ' ind =', i6, ' du =', f10.5)
9003  format('| ene =', f15.12, ' spline_ene =', f15.12, ' u =', f10.5)
#endif /* EFS_DBG */
    end do
  end do
  
  return

end subroutine chk_ene_spline

!*******************************************************************************
!
! Subroutine:  chk_frc_spline
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine chk_frc_spline(tbl_siz, del, beta, frc_tbl, max_relerr, &
                          err_lim_exceeded_at) 

  use mdin_ctrl_dat_mod

  implicit none

! Formal arguments:

  integer,          intent(in)  :: tbl_siz      ! logical table entries
  double precision, intent(in)  :: del          ! u increment per table entry
  double precision, intent(in)  :: beta         ! ewald coefficient
  double precision, intent(in)  :: frc_tbl(tbl_siz * 8)
  double precision, intent(in)  :: max_relerr
  double precision, intent(out) :: err_lim_exceeded_at

  ! err_lim_exceeded_at is the highest value of u with a relative error
  ! greater than max_relerr; if none, then it is -1.d0

! Local variables:

  double precision      :: u, du                ! u is synonym for delr2
                                                ! du applies to spline tbl
  double precision      :: df, spline_df
  double precision      :: err
  integer               :: i, j
  integer               :: ind
  integer               :: min_tau_idx, max_tau_idx

  min_tau_idx = 1
  max_tau_idx = ceiling(es_cutoff * es_cutoff) * efs_tbl_dens

  ! Sanity check.
  if (max_tau_idx .gt. tbl_siz) max_tau_idx = tbl_siz

  err_lim_exceeded_at = - 1.d0

  do i = min_tau_idx, max_tau_idx
    ind = ishft((i - 1), 3)       ! 8 * (i - 1)
    do j = 1, 9
      u = dble(i - 1) * del + dble(j) * 0.1d0 * del
      du = dble(j) * 0.1d0 * del
      call df_of_u(beta, u, df)

      spline_df = (frc_tbl(5 + ind) + du * &
                  (frc_tbl(6 + ind) + du * &
                  (frc_tbl(7 + ind) + du * &
                   frc_tbl(8 + ind))))

      err = abs(spline_df - df) / (abs(df) + small)
      if (err .gt. max_relerr) err_lim_exceeded_at = u

#ifdef EFS_DBG
      write(0,9001) err, u
!     write(0,9002) i,j,ind,du
      write(0,9003) df, spline_df, u
      9001  format('| DBG_CHECK df(r**2): rel err = ', e12.3, '   at ', f10.6)
      9002  format('| i =', i6, ' j =', i6, ' ind =', i6, ' du =', f10.5)
      9003  format('| df =', f15.12, ' spline_df =', f15.12, ' u =', f10.5)
#endif /* EFS_DBG */
    end do
  end do
  
  return

end subroutine chk_frc_spline

!*******************************************************************************
!
! Subroutine:  e_of_u
!
! Description:
!
! Result of e_of_u should be multiplied by charge product if you are
! actually calculating energies, but if you are filling a spline table
! you just use the result.
!              
!*******************************************************************************

subroutine e_of_u(beta, u, ene)

  implicit none

! Formal arguments:

  double precision, intent(in)  :: beta         ! ewald coefficient
  double precision, intent(in)  :: u            ! interatom distance squared
  double precision, intent(out) :: ene          ! calculated energy

! Local variables;

  double precision              :: switch       ! erfc(beta * r)
  double precision              :: r            ! interatom distance

  r = sqrt(u)
  call derfcfun(beta * r, switch)

  ene = switch / r

  return

end subroutine e_of_u


!*******************************************************************************
!
! Subroutine:  dedu_of_u
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine dedu_of_u(beta, u, dedu)

use gbl_constants_mod

  implicit none

! Formal arguments:

  double precision, intent(in)  :: beta         ! ewald coefficient
  double precision, intent(in)  :: u            ! interatom distance squared
  double precision, intent(out) :: dedu         ! calculated deriv of energy
                                                ! function

! Local variables;

  double precision              :: switch       ! erfc(beta * r)
  double precision              :: r            ! interatom distance

  r = sqrt(u)
  call derfcfun(beta * r, switch)

  dedu = -beta * exp(-(beta * beta * u)) / (sqrt(PI) * u) - &
          switch / (2.d0 * u * r)

  return

end subroutine dedu_of_u

!*******************************************************************************
!
! Subroutine:  df_of_u
!
! Description:
!              
! Result of df_of_u should be multiplied by charge product if you are
! actually calculating energies, but if you are filling a spline table
! you just use the result.
!*******************************************************************************

subroutine df_of_u(beta, u, df)

use gbl_constants_mod

  implicit none

! Formal arguments:

  double precision, intent(in)  :: beta         ! ewald coefficient
  double precision, intent(in)  :: u            ! interatom distance squared
  double precision, intent(out) :: df           ! calculated force

! Local variables;

  double precision              :: r            ! interatom distance
  double precision              :: switch       ! erfc(beta * r)

  r = sqrt(u)
  call derfcfun(beta * r, switch)

  df =  2.d0 * exp(-(beta * beta * u)) * beta / (sqrt(PI) * u) + switch / (u * r)

  return

end subroutine df_of_u

!*******************************************************************************
!
! Subroutine:  ddfdu_of_u
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine ddfdu_of_u(beta, u, ddfdu)

use gbl_constants_mod

  implicit none

! Formal arguments:

  double precision, intent(in)  :: beta         ! ewald coefficient
  double precision, intent(in)  :: u            ! interatom distance squared
  double precision, intent(out) :: ddfdu        ! calculated deriv of force
                                                ! function

! Local variables;

  double precision              :: r            ! interatom distance
  double precision              :: switch       ! erfc(beta * r)
  double precision              :: cterm        ! common exponential term

  r = sqrt(u)
  call derfcfun(beta * r, switch)

  cterm =  -exp(-(beta * beta * u)) * beta / (sqrt(PI) * u)

  ddfdu = cterm * (2.d0 * beta * beta + 3.d0 / u) - &
          3.d0 * switch / (2.d0 * u * u * r)

  return

end subroutine ddfdu_of_u

end module ene_frc_splines_mod



  