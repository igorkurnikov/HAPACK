mset = GetCurMolSet()
mm_mod = mset.GetMolMechMod(1)
mm_model = mm_mod.GetMolMechModel()
mm_model.InitModel(AMOEBA)
mm_model.SetRestrainedAtoms("MG_AT")
mm_model.SetRestrRefCrdFromXYZFile("mg_pos_constr.xyz")
mm_model.SetDistConstrFromFile("constr_rmo8_1.9_k100.dat")
mm_mod.SetMMRunType(MD_RUN)
#mm_model.SetNBCutDist(6.0)
mm_mod.SetNumMDSteps(100)
mm_mod.SetWrtRstrtFreq(100)
mm_mod.SetWrtMDTrajFreq(100)
mm_mod.SetWrtLogFreq(10)
mm_mod.SetWrtConstrFreq(1)
mm_mod.Run()



