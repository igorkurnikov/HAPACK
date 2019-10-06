mset = GetCurMolSet()
print(mset.GetNAtoms())
mm_mod = mset.GetMolMechMod(1)
mm_model = mm_mod.GetMolMechModel()
mm_mod.Initialize()
mm_mod.SetMMRunType(MD_RUN)
mm_mod.SetRunInternal() 
mm_model.SetNBCutDist(6.0)
mm_mod.SetNumMDSteps(1000)
mm_mod.SetWrtRstrtFreq(100)
mm_mod.SetWrtMDTrajFreq(100)
mm_mod.SetWrtLogFreq(10)
mm_mod.Run()



