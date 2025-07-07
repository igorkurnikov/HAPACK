#
# This example script to run thermodynamic Integration 
#
pmset = GetCurMolSet()
pmset.SetName("methane_wat_1_ti")
p_mm_mod = pmset.GetMolMechMod(1)
p_ti_mod = p_mm_mod.p_ti_mod
p_ti_mod.num_lmb = 5
p_ti_mod.num_equilib_points = 2000
delta_g = p_ti_mod.CalcDeltaG()
p_ti_mod.ReduceDvDlData(10000)
#print "Delta G = ",delta_g


