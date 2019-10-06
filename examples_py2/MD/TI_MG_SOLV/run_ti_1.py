#
# This example script to run thermodynamic Integration 
#

mset_1 = HaMolSet()
mset_2 = HaMolSet()
mset_1.LoadHarlemFile("mg_neutral_wat_1.hlm")
mset_2.LoadHarlemFile("mg_wat_1.hlm")
mset_1.SetName("mg_wat_1_ti")
mm_mod_1 = mset_1.GetMolMechMod(1)
mm_mod_1.InitMolMechModel(AMBER_10)
mm_mod_2 = mset_2.GetMolMechMod(1)
mm_mod_2.InitMolMechModel(AMBER_10)
model_2 = mm_mod_2.GetMolMechModel()
ti_mod = mm_mod_1.GetTISimMod()

#mm_mod_1.SetPerBoundaryCondType(CONST_VOL)
mm_mod_1.SetNumMDSteps(10000)
mm_mod_1.SetWrtRstrtFreq(500)
mm_mod_1.SetWrtEnerFreq(10)
mm_mod_1.SetWrtLogFreq(10)
#parameters of TI job (indexes need to check)
ti_mod.SetNumLambda(5)
ti_mod.SetCurIdxLambda(0)
ti_mod.SetMaxLambdaIdx(5)
# run TI
# this will save all dV/dlambda to file 
# name_lmb_CurrentLambda_TotalLambda_dvdl.dat
#
mm_mod_1.RunTI(model_2)
# stat TI calculation from step  # parameter
ti_mod.SetNumEqPoints(0)
delta_g = ti_mod.CalcDeltaG()
#print "Delta G = ", delta_g

