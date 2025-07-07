#
# This example script to run thermodynamic Integration 
# HIS(D) ->HIS(PROT)  capped by ACE and NME, ACE NME charge distributions from AMBER-03
#
mset_1 = HaMolSet()
mset_2 = HaMolSet()
mset_1.LoadHarlemFile("his_cap_wat_n_2.hlm")
mset_2.LoadHarlemFile("his_cap_wat_p_2.hlm")
mset_1.SetName("his_w")
mm_mod_1 = mset_1.GetMolMechMod(1)
mm_mod_1.InitMolMechModel(AMBER_10)
mm_mod_2 = mset_2.GetMolMechMod(1)
mm_mod_2.InitMolMechModel(AMBER_10)
model_2 = mm_mod_2.GetMolMechModel()
ti_mod = mm_mod_1.GetTISimMod()
#p_mm_mod_1.LoadAmberRestartFile("dtox_h172_lmb_05_t278.rst")
mm_mod_1.SetShakeConstr(H_ATOM_SHAKE)
mm_mod_1.SetMDTimeStep(0.002)
mm_mod_1.SetNumMDSteps(10000)
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

