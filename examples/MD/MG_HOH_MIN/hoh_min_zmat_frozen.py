mset = GetCurMolSet()
zmat = mset.GetZMat()
mm_mod = mset.GetMolMechMod(1)
mm_mod.SetZMatMin()
mm_mod.SetMaxNumMinimSteps(500)
mm_model = mm_mod.GetMolMechModel()
mm_model.InitModel(AMBER_10)
mm_model.SetRestrainedAtoms("O_AT")
mm_model.SetRestrRefCrdFromStr("0.0 0.0 0.0")
elc_a = zmat.GetCrdByTag("aHOH")
#ang = 1.824217
ang = 1.7
print ang * RAD_TO_DEG
elc_a.SetValue(ang)
#elc.SetFrozen()
mm_mod.RunMinEne()
ene = mm_mod.GetUnConstrEne()
print " ang = ", ang, " ene = ", ene

