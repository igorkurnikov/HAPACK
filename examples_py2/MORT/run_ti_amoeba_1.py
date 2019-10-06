#
# This example script to run thermodynamic Integration 
#

mset_1 = HaMolSet()
mset_1.SetName("agl_solv_3_ti")
mset_2 = HaMolSet()
mset_1.LoadHarlemFile("agl_solv_3.hlm")
mset_2.LoadHarlemFile("agl_solv_3.hlm")
mm_mod_1 = mset_1.GetMolMechMod(1)
mm_mod_1.InitMolMechModel(AMOEBA)
ti_simm = mm_mod_1.GetTISimulator()
#mm_mod_1.LoadAmberRestartFile("agl.rst")
#mm_mod_1.shake_constr = H_ATOM_SHAKE
mm_mod_1.SetMDTimeStep(0.001)
#mm_mod_1.SetTempCtrlMethod(CONST_TEMP_LANGEVIN)
mm_mod_2 = mset_2.GetMolMechMod(1)
mm_mod_2.InitMolMechModel(AMOEBA)
mm_model_2 = mm_mod_2.GetMolMechModel()
mm_mod_1.SetNumMDSteps(1000)
mm_mod_1.SetWrtRstrtFreq(100)
mm_mod_1.SetWrtMDTrajFreq(100)
mm_mod_1.SetWrtLogFreq(10)
ti_simm.SetNumLambda(3)
ti_simm.SetCurIdxLambda(0)
mm_mod_1.RunTI(mm_model_2)
ti_simm.SetNumEqPoints(0)
delta_g = ti_simm.CalcDeltaG()
ti_simm.ReduceDvDlData(500)


