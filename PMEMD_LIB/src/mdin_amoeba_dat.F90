#include "copyright.i"

!*******************************************************************************!
! Module: mdin_amoeba_dat_mod
!
! Description: <TBS>
!
!*******************************************************************************

module mdin_amoeba_dat_mod

  implicit none

  integer, save          ::  do_amoeba_valence, do_amoeba_nonbond, do_bond, &
                             do_ureyb, do_reg_angle, do_trig_angle, do_opbend, &
                             do_torsion, do_pi_torsion, do_strbend, &
                             do_torsion_torsion, do_str_torsion, do_recip, &
                             do_adjust, do_direct, do_self, do_vdw, &
                             do_induced, do_vdw_taper, do_vdw_longrange, &
                             beeman_integrator, dipole_scf_iter_max, &
                             amoeba_verbose

  double precision,save  ::     compress, dipole_scf_tol, ee_dsum_cut, &
                                ee_damped_cut, sor_coefficient, &
                                thole_expon_coeff, vdw_taper

contains

!*******************************************************************************!
! Subroutine:  init_mdin_amoeba_dat
!
! Description: <TBS>
!
!*******************************************************************************

subroutine init_mdin_amoeba_dat

  use file_io_dat_mod
  use file_io_mod

  implicit none

  call am_val_set_user_bit(do_bond, do_ureyb, do_reg_angle, do_trig_angle, &
                           do_opbend, do_torsion, do_str_torsion, &
                           do_pi_torsion, do_strbend, &
                           do_torsion_torsion)

  call am_nonbond_set_user_bit(do_recip, do_adjust, do_direct, do_self, &
                               do_vdw, do_induced)

  return
  
end subroutine init_mdin_amoeba_dat


end module mdin_amoeba_dat_mod
