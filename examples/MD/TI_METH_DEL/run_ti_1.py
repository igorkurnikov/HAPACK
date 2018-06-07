#
# This example script to run thermodynamic Integration 
#

pmset_1 = HaMolSet()
pmset_1.SetName("methane_del_wat_1_ti")
pmset_2 = HaMolSet()
pmset_1.LoadHarlemFile("methane_ch0_wat_1.hlm")
pmset_2.LoadHarlemFile("methane_dum_wat_1.hlm")
p_mm_mod_1 = pmset_1.GetMolMechMod(1)
p_mm_mod_1.Initialize()
p_mm_mod_1.LoadAmberRestartFile("methane_del_wat_lmb_1_t500.rst")
p_mm_mod_2 = pmset_2.GetMolMechMod(1)
p_mm_mod_2.Initialize()
p_mm_mod_1.period_bcond = CONST_VOL
p_mm_mod_1.pressure_reg_method = NO_CRD_SCALING
p_mm_mod_1.nonb_cutoff_dist = 6.0/BOHR_TO_ANG
p_mm_mod_2.period_bcond = CONST_VOL
p_mm_mod_2.pressure_reg_method = NO_CRD_SCALING
p_mm_mod_2.nonb_cutoff_dist = 6.0/BOHR_TO_ANG
p_mm_mod_1.klambda_ti = 6.0
p_mm_mod_1.lambda_ti  = 0.04691
#p_mm_mod_1.lambda_ti  = 0.23076
#p_mm_mod_1.lambda_ti  = 0.5
#p_mm_mod_1.lambda_ti  = 0.76923
#p_mm_mod_1.lambda_ti  = 0.95308
p_mm_mod_1.ti_sync_freq = 10
p_mm_mod_1.length_md_run  =  10000000
p_mm_mod_1.nbstp_wrt_rstrt =  5000
p_mm_mod_1.nbstp_wrt_coord =  1000
p_mm_mod_1.nbstp_wrt_ener =  1000
p_mm_mod_1.print_freq          =  100
p_mm_mod_1.RunTI(p_mm_mod_2)


