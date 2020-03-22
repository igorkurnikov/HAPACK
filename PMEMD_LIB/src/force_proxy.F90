!
! Proxy function to call force subroutines
!
subroutine pme_force_proxy(atm_cnt,crd,saved_crd,box,saved_box,vel,frc,mass,new_list, &
           atm_jrc, atm_xc, atm_weight, igroup, natc, si, virial, ekcmt, pme_err_est, &
           atm_owner_map, my_atm_lst, my_atm_cnt, tranvec, imin_par )

   use pme_force_mod
   use amoeba_force_mod
   use state_info_mod
   use parallel_dat_mod
   use img_mod
   use mdin_ctrl_dat_mod, only : iamoeba
   
   implicit none

   ! Formal Arguments:

   integer           ::   atm_cnt
   double precision  ::   crd(3,atm_cnt)
   double precision  ::   saved_crd(3,atm_cnt)
   double precision  ::   box(3)
   double precision  ::   saved_box(3)
   double precision  ::   vel(3,atm_cnt)
   double precision  ::   frc(3,atm_cnt)
   double precision  ::   mass(atm_cnt)
   integer           ::   new_list
   integer           ::   atm_jrc(*)
   double precision  ::   atm_xc(3,*)
   double precision  ::   atm_weight(*) 
   integer           ::   igroup(*)
   integer           ::   natc
   double precision  ::   si(si_cnt)
   double precision    :: virial(3)            ! Only used for MD
   double precision    :: ekcmt(3)             ! Only used for MD
   double precision    :: pme_err_est          ! Only used for MD
   integer             :: atm_owner_map(atm_cnt)
   integer             :: my_atm_lst(*)
   integer             :: my_atm_cnt
   double precision    :: tranvec(*)
   integer             :: imin_par
   
   type(pme_pot_ene_rec)  :: pme_pot_ene
   type(amba_pot_ene_rec) :: amba_pot_ene 
	
	logical new_list_loc;
	integer i
	
	if(new_list .eq. 0) then
		new_list_loc = .false.
	else
		new_list_loc = .true.
	endif 

    si(:) = 0.0d0
    frc(:,:) = 0.0d0
    
    virial(:) = 0.0d0
    ekcmt(:)  = 0.0d0
    pme_err_est = 0.0d0
                                              
   if(iamoeba .eq. 1) then
        if( numtasks .gt. 1) then
            call amoeba_force_mpi(atm_cnt, crd, saved_crd, box, saved_box, frc, mass, gbl_img_atm_map, gbl_atm_img_map, &
                                  atm_owner_map, my_atm_lst, my_atm_cnt, new_list_loc, amba_pot_ene, si(si_diprms), & 
                                  si(si_dipiter), virial, &
                                  atm_jrc, atm_xc, atm_weight, igroup, natc, tranvec, imin_par)
        
        else
            call amoeba_force_uni(atm_cnt, crd, saved_crd, box, saved_box, frc, mass, gbl_img_atm_map, gbl_atm_img_map, &
                                  atm_owner_map, my_atm_lst, my_atm_cnt, new_list_loc , amba_pot_ene, si(si_diprms), & 
                                  si(si_dipiter), virial, &
                                  atm_jrc, atm_xc, atm_weight, igroup, natc, tranvec, imin_par)
        endif
    else
        call pme_force(atm_cnt, crd, saved_crd, box, saved_box, vel, frc, mass, gbl_img_atm_map, gbl_atm_img_map, &
                       atm_owner_map, my_atm_lst, my_atm_cnt, & 
                       new_list_loc, pme_pot_ene, virial, ekcmt,pme_err_est, atm_jrc, atm_xc, atm_weight, igroup, & 
                       natc, tranvec, imin_par)
    endif

   ! Store energy terms in state info array for printout.

    if(iamoeba .eq. 1) then
        si(si_pot_ene)   = amba_pot_ene%total
        si(si_vdw_ene)   = amba_pot_ene%vdw
        si(si_elect_ene) = amba_pot_ene%elec
        si(si_hbond_ene) = amba_pot_ene%hbond
        si(si_bond_ene)  = amba_pot_ene%bond
        si(si_angle_ene)     = amba_pot_ene%angle
        si(si_dihedral_ene)  = amba_pot_ene%dihedral
        si(si_vdw_14_ene)    = amba_pot_ene%vdw_14
        si(si_elect_14_ene)  = amba_pot_ene%elec_14
        si(si_restraint_ene) = amba_pot_ene%restraint  
        si(si_polar)         = amba_pot_ene%polar  
    else
        si(si_pot_ene) = pme_pot_ene%total
        si(si_vdw_ene) = pme_pot_ene%vdw_tot
        si(si_elect_ene) = pme_pot_ene%elec_tot
        si(si_hbond_ene) = pme_pot_ene%hbond
        si(si_bond_ene) = pme_pot_ene%bond
        si(si_angle_ene) = pme_pot_ene%angle
        si(si_dihedral_ene) = pme_pot_ene%dihedral
        si(si_vdw_14_ene) = pme_pot_ene%vdw_14
        si(si_elect_14_ene) = pme_pot_ene%elec_14
        si(si_restraint_ene) = pme_pot_ene%restraint
    endif

end subroutine pme_force_proxy
