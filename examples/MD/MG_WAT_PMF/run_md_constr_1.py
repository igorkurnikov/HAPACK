mset = GetCurMolSet()
dist = 2.5
prefix = mset.GetName()
prefix += "_" + ("%.2f" % dist)
mm_mod = mset.GetMolMechMod(1)
mm_model = mm_mod.GetMolMechModel()
mm_model.InitModel(AMBER_10)
#mm_model.InitModel(AMOEBA)
#mm_model.SetTholeExponCoef(0.095)
mm_model.SetRestrainedAtoms("MG_AT")
mm_model.SetRestrRefCrdFromStr("9.1 9.1 9.1")
mm_mod.SetMMRunType(MD_RUN)
mm_mod.SetNumMDSteps(10000)
mm_mod.SetWrtRstrtFreq(10000)
mm_mod.SetWrtEnerFreq(0)
mm_mod.SetWrtMDTrajFreq(1000)
mm_mod.SetWrtLogFreq(1000)
mm_mod.SetWrtConstrFreq(10)
for dist in (1.9,2.5):
  mm_model.SetHarmConstraint("$MG$MG1.MG","$SOLVENT$HOH8.O",dist,100.0)
  mm_mod.Run()



