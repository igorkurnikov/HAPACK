#
# This example script to run thermodynamic Integration 
#

pmset_1 = HaMolSet()
pmset_1.SetName("his_wat_ti")
pmset_2 = HaMolSet()
pmset_1.LoadHarlemFile("his_neutr_wat_1.hlm")
pmset_2.LoadHarlemFile("his_prot_wat_1.hlm")
p_mm_mod_1 = pmset_1.GetMolMechMod(1)
p_mm_mod_1.Initialize()
#p_mm_mod_1.LoadAmberRestartFile("his_neutr_wat_ti_lmb_0_3.rst")
p_mm_mod_1.shake_constr = H_ATOM_SHAKE
p_mm_mod_1.md_time_step = 0.002
p_mm_mod_1.temp_control_method = CONST_TEMP_LANGEVIN
p_mm_mod_2 = pmset_2.GetMolMechMod(1)
p_mm_mod_2.Initialize()
p_mm_mod_1.length_md_run  =  5000
p_mm_mod_1.wrt_rstrt_freq =  500
p_mm_mod_1.wrt_ener_freq =  100
p_mm_mod_1.wrt_log_freq   =  100
p_mm_mod_1.p_ti_mod.num_lmb = 3
p_mm_mod_1.p_ti_mod.cur_idx_lmb = 0
p_mm_mod_1.RunTI(p_mm_mod_2)
delta_g = p_mm_mod_1.p_ti_mod.CalcDeltaG()
p_mm_mod_1.p_ti_mod.ReduceDvDlData(1000)


