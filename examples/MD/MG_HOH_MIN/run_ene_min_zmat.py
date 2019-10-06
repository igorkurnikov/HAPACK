mset = HaMolSet()
mset.LoadHarlemFile("MG_HOH1_1.hlm")
#zmat = mset.GetZMat()
#zmat.InitAllXYZ()
mm_mod = mset.GetMolMechMod(1)
mm_mod.SetZMatMin()
mm_mod.SetMaxNumMinimSteps(5000)
mm_model = mm_mod.GetMolMechModel()
mm_model.InitModel(AMBER_10)
mm_model.SetRestrainedAtoms("MG_AT")
mm_model.SetRestrRefCrdFromStr("0.0 0.0 0.0")
dist_arr = []
ene_arr = []
dist = 1.3
first = 1
while( dist < 1.31):
  mm_model.SetHarmConstraint("MG1.MG","HOH2.O",dist,100.0)
  if(not first): mm_model.UpdateConstraints()
  first = 0
  mm_mod.RunMinEne()
  ene = mm_mod.GetUnConstrEne()
  print(" dist = ", dist, " ene = ", ene)
  dist_arr.append(dist)
  ene_arr.append(ene)
  dist = dist + 0.05
n = len(dist_arr)
for i in range(n):
  print(dist_arr[i]," ",ene_arr[i])
  

