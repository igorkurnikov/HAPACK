#
# This example script to run thermodynamic Integration 
#

mset_1 = HaMolSet()
mset_1.SetName("methane_wat_1_ti")
mset_2 = HaMolSet()
mset_1.LoadHarlemFile("methane_wat_1.hlm")
mset_2.LoadHarlemFile("methane_wat_1_ch1.hlm")
mm_mod_1 = mset_1.GetMolMechMod(1)
mm_mod_1.Initialize()
mm_mod_2 = mset_2.GetMolMechMod(1)
mm_mod_2.Initialize()
model_2 = mm_mod_2.GetMolMechModel()
ti_mod = mm_mod_1.p_ti_mod
mm_mod_1.SetNumMDSteps(10000)
mm_mod_1.SetWrtRstrtFreq(500)
mm_mod_1.SetWrtEnerFreq(100)
mm_mod_1.SetWrtLogFreq(100)
ti_mod.SetNumLambda(3)
mm_mod_1.RunTI(model_2)
delta_g = ti_mod.CalcDeltaG()
#print "Delta G = ",delta_g

