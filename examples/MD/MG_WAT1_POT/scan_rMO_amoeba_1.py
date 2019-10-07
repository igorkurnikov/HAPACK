mset = GetCurMolSet()
zm = mset.GetZMat()
mm_mod = mset.GetMolMechMod(1)
mm_model = mm_mod.GetMolMechModel()
#mm_model.SetFFTGridsPerAng(0.025)
mm_model.SetUseMortLib(0)
mm_model.SetPMECoef( 0.000015 )
mm_model.InitModel(AMOEBA)
#mm_model.SetCalcRecipSpaceInter(0)
#mm_model.SetCalcSelfInter(0)
mm_model.SetTholeExponCoef(0.0)
#mm_model.SetTholeExponCoef(0.39)
#mm_model.SetTholeExponCoef(10.0)
#mm_model.SetTholeExponCoef(0.095)
mm_model.SetDipoleScfTol(0.000001)
b1 = zm.GetCrdByTag("rMO")
b1.SetFrozen()
dist = 1.5
dist_arr = []
ene_arr = []
while( dist < 5.01 ):
  b1.SetValue(dist)
  zm.SetAtomCrd()
  mm_mod.CalcEnergy()
  ene = mm_mod.GetEne()
  print("dist = ", dist, " ene = ",ene)
  dist_arr.append(dist)
  ene_arr.append(ene)
  dist = dist + 0.05
n = len(dist_arr)
for i in range(n):
  print(dist_arr[i]," ",ene_arr[i])

