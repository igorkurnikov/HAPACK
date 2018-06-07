#include "copyright.i"

!*******************************************************************************
!
! Module: gb_force_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module gb_force_mod

  use gbl_datatypes_mod

  implicit none

  ! Potential energies, with breakdown, from GB.  This is intended to be the
  ! external interface to potential energy data produced by this module, in
  ! particular by subroutine gb_force().

  type gb_pot_ene_rec
    sequence
    double precision    :: total
    double precision    :: vdw_tot
    double precision    :: elec_tot
    double precision    :: gb
    double precision    :: bond
    double precision    :: angle
    double precision    :: dihedral
    double precision    :: vdw_14
    double precision    :: elec_14
    double precision    :: restraint
  end type gb_pot_ene_rec

  integer, parameter    :: gb_pot_ene_rec_size = 10

  type(gb_pot_ene_rec), parameter      :: null_gb_pot_ene_rec = &
    gb_pot_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)

contains

!*******************************************************************************
!
! Subroutine:  gb_force
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine gb_force(atm_cnt, crd, frc, mass, pot_ene, irespa, &
                    atm_jrc, atm_xc, atm_weight, igroup, belly_atm_cnt, natc, &
                    atm_owner_map, my_atm_lst, my_atm_cnt )

  use angles_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use dist_constr_mod
  use gb_ene_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use parallel_mod
  use prmtop_dat_mod
  use runfiles_mod
  use timers_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3, atm_cnt)
  double precision              :: frc(3, atm_cnt)
  double precision              :: mass(atm_cnt)
  type(gb_pot_ene_rec)          :: pot_ene
  integer, intent(in)           :: irespa
  integer                       :: atm_jrc(*)
  double precision              :: atm_xc(3, *)
  double precision              :: atm_weight(*)
  integer                       :: igroup(*)
  integer                       :: belly_atm_cnt
  integer                       :: natc
  integer                       :: atm_owner_map(atm_cnt) 
  integer                       :: my_atm_lst(atm_cnt)
  integer                       :: my_atm_cnt
  
  
! Local variables:

  integer                       :: i, j

! Zero energies that are stack or call parameters:

  pot_ene = null_gb_pot_ene_rec

! If no force calcs are to be done, clear the frc array and bag out now.

  if (ntf .eq. 8) then
    frc(:,:) = 0.d0
    return
  end if

  call zero_time()
  call zero_gb_time()

! Calculate the non-bonded contributions:

  frc(:,:) = 0.d0

  call gb_ene(crd, frc, atm_gb_radii, atm_gb_fs, atm_qterm, atm_iac, &
              typ_ico, atm_numex, gbl_natex, atm_cnt, belly_atm_cnt, &
              pot_ene%gb, pot_ene%elec_tot, pot_ene%vdw_tot, irespa)

  call update_time(nonbond_time)
  
! Calculate the other contributions:

  call gb_bonded_force(crd, frc, pot_ene, atm_jrc, atm_xc, atm_weight, natc, atm_owner_map )

  ! Sum up total potential energy for this task:

  pot_ene%total = pot_ene%vdw_tot + &
                  pot_ene%elec_tot + &
                  pot_ene%gb + &
                  pot_ene%bond + &
                  pot_ene%angle + &
                  pot_ene%dihedral + &
                  pot_ene%vdw_14 + &
                  pot_ene%elec_14 + &
                  pot_ene%restraint

  if( numtasks .gt. 1) then
! Distribute forces to atom owners; this is a reduce sum operation:
    
    call gb_frcs_distrib(atm_cnt, frc, dbl_mpi_recv_buf)

! Distribute energies to all processes.

    call gb_distribute_enes(pot_ene)

  endif

  return

end subroutine gb_force

!*******************************************************************************
!
! Subroutine:  gb_bonded_force
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine gb_bonded_force(crd, frc, pot_ene, atm_jrc, atm_xc, atm_weight, natc, atm_owner_map )

  use angles_mod
  use bonds_mod
  use constraints_mod
  use dist_constr_mod
  use dihedrals_mod
  use dynamics_dat_mod
  use mdin_ctrl_dat_mod
  use timers_mod

  implicit none

! Formal arguments:

  double precision              :: crd(3, *)
  double precision              :: frc(3, *)
  type(gb_pot_ene_rec)          :: pot_ene
  integer                       :: atm_jrc(*)
  double precision              :: atm_xc(3,*)
  double precision              :: atm_weight(*)
  integer                       :: natc
  integer                       :: atm_owner_map(*) 

! Local variables:

  ! These energy variables are temporaries, for summing. DON'T use otherwise!

  double precision              :: bond_ene
  double precision              :: angle_ene
  double precision              :: dihedral_ene
  double precision              :: vdw_14_ene
  double precision              :: elec_14_ene
  double precision              :: nb_constr_ene

  double precision              :: molvir(3, 3)         ! dummy argument
  double precision              :: e14_vir(3, 3)        ! dummy argument

! Bond energy contribution:

  if (ntf .le. 1) then
    if (cit_nbonh .gt. 0) then
      call get_bond_energy(cit_nbonh, cit_h_bond, crd, frc, bond_ene, atm_owner_map)
      pot_ene%bond = bond_ene
    end if
  end if

  if (ntf .le. 2) then
    if (cit_nbona .gt. 0) then
      call get_bond_energy(cit_nbona, cit_a_bond, crd, frc, bond_ene, atm_owner_map)
      pot_ene%bond = pot_ene%bond + bond_ene
    end if
  end if
  
  call update_time(bond_time)

! Angle energy contribution:

  if (ntf .le. 3) then
    if (cit_ntheth .gt. 0) then
      call get_angle_energy(cit_ntheth, cit_angle, crd, frc, angle_ene, atm_owner_map)
      pot_ene%angle = angle_ene
    end if
  end if

  if (ntf .le. 4) then
    if (cit_ntheta .gt. 0) then
      call get_angle_energy(cit_ntheta, cit_angle(cit_ntheth+1), &
                            crd, frc, angle_ene, atm_owner_map)
      pot_ene%angle = pot_ene%angle + angle_ene
    end if
  end if

  call update_time(angle_time)

! Dihedral energy contribution:

  if (ntf .le. 5) then
    if (cit_nphih .gt. 0) then
      call get_dihed_energy(atm_owner_map, cit_nphih, cit_dihed, crd, frc, &
                            dihedral_ene, vdw_14_ene, elec_14_ene, &
                            atm_rel_crd, molvir, e14_vir )
                                                       
      pot_ene%dihedral = dihedral_ene
      pot_ene%vdw_14 = vdw_14_ene
      pot_ene%elec_14 = elec_14_ene
    end if
  end if

  if (ntf .le. 6) then
    if (cit_nphia .gt. 0) then
      call get_dihed_energy(atm_owner_map, cit_nphia, cit_dihed(cit_nphih + 1), crd, frc, &
                            dihedral_ene, vdw_14_ene, elec_14_ene, &
                            atm_rel_crd, molvir, e14_vir )
                                                     
      pot_ene%dihedral = pot_ene%dihedral + dihedral_ene
      pot_ene%vdw_14 = pot_ene%vdw_14 + vdw_14_ene
      pot_ene%elec_14 = pot_ene%elec_14 + elec_14_ene
    end if
  end if

  call update_time(dihedral_time)

! Calculate the position restraints energy:

  if (natc .gt. 0) then
     call get_crd_constraint_energy(natc, pot_ene%restraint, atm_jrc, &
                                    crd, frc, atm_xc, atm_weight, atm_owner_map )
  endif
  
! Calculate non-bonded constraints energy

  if( cit_num_dist_constr .gt. 0 ) then
	 
	 call get_dist_constr_energy(cit_num_dist_constr, cit_dist_constr, crd, frc, nb_constr_ene, atm_owner_map)
	 pot_ene%restraint = pot_ene%restraint + nb_constr_ene
  
  endif

  return

end subroutine gb_bonded_force

!*******************************************************************************
!
! Subroutine:  gb_distribute_enes
!
! Description: <TBS>
!              Used only in MPI
!*******************************************************************************

subroutine gb_distribute_enes(pot_ene)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  type(gb_pot_ene_rec) :: pot_ene

! Local variables:

  type(gb_pot_ene_rec), save    :: dat_in, dat_out

  dat_in = pot_ene

  call mpi_allreduce(dat_in%total, dat_out%total, &
                     gb_pot_ene_rec_size, mpi_double_precision, &
                     mpi_sum, lib_mpi_comm, err_code_mpi)

  pot_ene = dat_out

  return

end subroutine gb_distribute_enes

end module gb_force_mod

subroutine gb_force_proxy(atm_cnt,crd,frc,mass,si,ncalls, &
                          atm_jrc, atm_xc, atm_weight, igroup, belly_atm_cnt, natc, &
                          atm_owner_map, my_atm_lst, my_atm_cnt )

   use gb_force_mod
   use state_info_mod

   integer           ::   atm_cnt
   double precision  ::   crd(3,atm_cnt)
   double precision  ::   frc(3,atm_cnt)
   double precision  ::   mass(atm_cnt)
   double precision  ::   si(si_cnt)
   integer           ::   ncalls
   integer           ::   atm_jrc(*)
   double precision  ::   atm_xc(3,*)
   double precision  ::   atm_weight(*) 
   integer           ::   igroup(*)
   integer           ::   belly_atm_cnt
   integer           ::   natc
   integer           ::   atm_owner_map(atm_cnt)
   integer           ::   my_atm_lst(*)
   integer           ::   my_atm_cnt
   
   
   type(gb_pot_ene_rec)  :: gb_pot_ene
  
   integer atm_cnt_new

   si(:) = 0.0d0
   frc(:,:) = 0.0d0
   
   atm_cnt_new = atm_cnt
   
!   write(*,*) " gb_force_proxy() pt 1 "

    call gb_force(atm_cnt_new, crd, frc, mass, gb_pot_ene, ncalls, &
                  atm_jrc, atm_xc, atm_weight, igroup, belly_atm_cnt, natc, &
                  atm_owner_map, my_atm_lst, my_atm_cnt)

    si(si_pot_ene) = gb_pot_ene%total
    si(si_vdw_ene) = gb_pot_ene%vdw_tot
    si(si_elect_ene) = gb_pot_ene%elec_tot
    si(si_hbond_ene) = gb_pot_ene%gb                ! temporary hack
    si(si_bond_ene) = gb_pot_ene%bond
    si(si_angle_ene) = gb_pot_ene%angle
    si(si_dihedral_ene) = gb_pot_ene%dihedral
    si(si_vdw_14_ene) = gb_pot_ene%vdw_14
    si(si_elect_14_ene) = gb_pot_ene%elec_14
    si(si_restraint_ene) = gb_pot_ene%restraint
    si(si_pme_err_est) = 0.d0
    
!    write(*,*) " gb_force_proxy() pt end "

end subroutine gb_force_proxy

