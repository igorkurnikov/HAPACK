!*******************************************************************************
!
! Module: mdin_ewald_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module mdin_ewald_dat_mod

use file_io_dat_mod

  implicit none

  integer, save                 ::  nfft1, nfft2, nfft3, bspl_order, netfrc, &
                                    vdwmeth, verbose_pme

  double precision, save        :: dsum_tol, ew_coeff, skinnb, eedtbdns, &
                                   fft_grids_per_ang
  

contains

!*******************************************************************************
!
! Subroutine:  compute_nfft
!
! Description:  This routine uses check_prime_factors to get the smallest
!               integer that will give at least the indicated grids per angstrom
!               given the specified box length (one axis, of course), with the
!               restriction that the integer have only prime factors of
!               2, 3, 5, and for fftw only, 7 (basically, check_prime_factors
!               is provided by the fft implementation being used, and that
!               determines the acceptable prime factors).
!*******************************************************************************

subroutine compute_nfft(box_len, grids_per_ang, nfft)

  use fft1d_mod

  implicit none

! Formal arguments:

  double precision, intent(in)  :: box_len
  double precision, intent(in)  :: grids_per_ang
  integer, intent(out)          :: nfft

! Local variables:

  integer               :: candidate_grid_cnt
  integer               :: ret_code

  candidate_grid_cnt = ceiling(box_len * grids_per_ang)

  do
    call check_prime_factors(candidate_grid_cnt, ret_code)
    if (ret_code .eq. 1) then
      nfft = candidate_grid_cnt
      exit
    end if
    candidate_grid_cnt = candidate_grid_cnt + 1
  end do

  return

end subroutine compute_nfft

!*******************************************************************************
!
! Subroutine:  compute_even_nfft
!
! Description:  This routine uses check_prime_factors to get the smallest
!               integer that will give at least the indicated grids per angstrom
!               given the specified box length (one axis, of course), with the
!               restriction that the integer have only prime factors of
!               2, 3, 5, and for fftw only, 7 (basically, check_prime_factors
!               is provided by the fft implementation being used, and that
!               determines the acceptable prime factors).  This routine also
!               REQUIRES that one of the prime factors is two.
!
!*******************************************************************************

subroutine compute_even_nfft(box_len, grids_per_ang, nfft)

  use fft1d_mod

  implicit none

! Formal arguments:

  double precision, intent(in)  :: box_len
  double precision, intent(in)  :: grids_per_ang
  integer, intent(out)          :: nfft

! Local variables:

  integer               :: candidate_grid_cnt
  integer               :: ret_code

  candidate_grid_cnt = ceiling(box_len * grids_per_ang)
  if (mod(candidate_grid_cnt, 2) .ne. 0) &
    candidate_grid_cnt = candidate_grid_cnt + 1
  do
    call check_prime_factors(candidate_grid_cnt, ret_code)
    if (ret_code .eq. 1) then
      nfft = candidate_grid_cnt
      exit
    end if
    candidate_grid_cnt = candidate_grid_cnt + 2
  end do

  return

end subroutine compute_even_nfft

!*******************************************************************************
!
! Subroutine:  compute_ew_coeff
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine compute_ew_coeff(ew_coeff_par, cutoff,dsum_tol_par)

  implicit none

! Formal arguments:

  double precision      :: ew_coeff_par
  double precision      :: cutoff
  double precision      :: dsum_tol_par

! Local variables:

  integer               :: i, n

  double precision      :: erfc_val
  double precision      :: term, x, xlo, xhi, y

! First get direct sum tolerance. How big must ew_coeff be to get
! terms outside the cutoff below tol?

  x = 0.5d0
  i = 0

  do
    x = 2.d0 * x
    i = i + 1
    y = x * cutoff
    call derfcfun(y, erfc_val)
    term = erfc_val/cutoff
    if (term .lt. dsum_tol_par) exit
  end do

! binary search tolerance is 2 to the -50th

  n = i + 50
  xlo = 0.d0
  xhi = x

  do i = 1, n
    x = (xlo + xhi)/2
    y = x * cutoff
    call derfcfun(y, erfc_val)
    term = erfc_val/cutoff
    if (term .ge. dsum_tol_par) then
      xlo = x
    else 
      xhi = x
    end if
  end do

  ew_coeff_par = x

  return

end subroutine compute_ew_coeff

subroutine set_add_pmemd_int_pars(vals)
! Set values of variables not in mdin_ctrl_int, mdin_ctrl_dbl and parmtop_int  
	
	use mdin_ctrl_dat_mod
	use parallel_dat_mod
	
	implicit none
	integer, intent(in)    :: vals(*)

	using_pme_potential   = vals(1)
	using_gb_potential    = vals(2) 
	dbg_atom_redistribution = vals(3)
	loadbal_verbose       = vals(4)  
	master                = vals(5)
	
end subroutine set_add_pmemd_int_pars

subroutine set_pme_int_pars_to_fort(vals)
! Set values of integer variables of mdin_ewald_dat_mod 
	
	use mdin_ctrl_dat_mod
	use parallel_dat_mod
	
	implicit none
	integer, intent(in)    :: vals(*)

	nfft1   = vals(1);
	nfft2   = vals(2);
	nfft3   = vals(3);
	bspl_order  = vals(4);
	netfrc      = vals(5);
	vdwmeth     = vals(6);
	verbose_pme = vals(7);
	
end subroutine set_pme_int_pars_to_fort

subroutine set_pme_dbl_pars_to_fort(vals)
! Set values of double variables of mdin_ewald_dat_mod 
	
	use mdin_ctrl_dat_mod
	use parallel_dat_mod
	
	implicit none
	double precision, intent(in)    :: vals(*)

	skinnb   = vals(1)
	dsum_tol = vals(2);
	ew_coeff = vals(3);
	eedtbdns = vals(4);
	fft_grids_per_ang = vals(5);
	
end subroutine set_pme_dbl_pars_to_fort

subroutine set_netfrc(netfrc_par)

    implicit none
	integer, intent(in)    :: netfrc_par

    netfrc = netfrc_par
    
end subroutine set_netfrc

end module mdin_ewald_dat_mod

