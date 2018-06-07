mset = HaMolSet()
mset.LoadHarlemFile("MG_HOH1_1.hlm")
zmat = mset.GetZMat()
mm_mod = mset.GetMolMechMod(1)
mm_mod.SetZMatMin()
mm_mod.SetMaxNumMinimSteps(1000)
mm_model = mm_mod.GetMolMechModel()
mm_model.InitModel(AMOEBA)
mm_model.SetTholeExponCoef(0.095)  # for magnesium in AMOEBA
mm_model.SetRestrainedAtoms("MG_AT")
mm_model.SetRestrRefCrdFromStr("0.0 0.0 0.0")
dist_arr = []
ene_arr = []
dist = 1.5
first = 1
while( dist < 5.1):
  elc = zmat.GetCrdByTag("rMO1")
  elc.SetValue(dist)
  elc.SetFrozen()
  first = 0
  mm_mod.RunMinEne()
  ene = mm_mod.GetUnConstrEne()
  print " dist = ", dist, " ene = ", ene
  dist_arr.append(dist)
  ene_arr.append(ene)
  dist = dist + 0.05
n = len(dist_arr)
for i in range(n):
  print dist_arr[i]," ",ene_arr[i]
  

