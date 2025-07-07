#
# This example script to run thermodynamic Integration 
#

mset_1 = HaMolSet()
mset_1.SetName("mg_wat_1_ti")
mset_2 = HaMolSet()
mset_1.LoadHarlemFile("mg_neutral_wat_1.hlm")
mset_2.LoadHarlemFile("mg_wat_1.hlm")
mm_mod_1 = mset_1.GetMolMechMod(1)
mm_mod_1.InitMolMechModel(AMOEBA)
mm_mod_2 = mset_2.GetMolMechMod(1)
mm_mod_2.InitMolMechModel(AMOEBA)
model_2 = mm_mod_2.GetMolMechModel()
ti_mod = mm_mod_1.GetTISimMod()
#mm_mod_1.SetPerBoundaryCondType(CONST_VOL)
mm_mod_1.SetNumMDSteps(500000)
mm_mod_1.SetWrtRstrtFreq(500)
mm_mod_1.SetWrtEnerFreq(10)
mm_mod_1.SetWrtLogFreq(10)
ti_mod.SetNumLambda(5)
ti_mod.SetCurIdxLambda(0)
ti_mod.SetMaxLambdaIdx(5)
mm_mod_1.RunTI(model_2)
ti_mod.SetNumEqPoints(0)
delta_g = ti_mod.CalcDeltaG()
#print "Delta G = ",delta_g

