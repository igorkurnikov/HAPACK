mset = GetCurMolSet()
mm_mod = mset.GetMolMechMod(1)
mm_model = mm_mod.GetMolMechModel()
#mm_model.InitModel(AMBER_10)
mm_model.InitModel(AMOEBA)
mm_model.SetTholeExponCoef(0.0)
mm_model.SetRestrainedAtoms("MG_AT")
mm_model.SetRestrRefCrdFromStr("9.1 9.1 9.1")
mm_mod.SetMMRunType(MD_RUN)
mm_mod.SetNumMDSteps(10000000)
mm_mod.SetWrtRstrtFreq(10000)
mm_mod.SetWrtEnerFreq(0)
mm_mod.SetWrtMDTrajFreq(10)
mm_mod.SetWrtLogFreq(1000)
mm_mod.Run()



