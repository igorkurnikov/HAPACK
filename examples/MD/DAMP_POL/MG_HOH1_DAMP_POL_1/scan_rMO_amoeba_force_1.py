mset = GetCurMolSet()
zm = mset.GetZMat()
aptr1 = mset.GetAtomByRef("HOH2.O")
mm_mod = mset.GetMolMechMod(1)
mm_model = mm_mod.GetMolMechModel()
amber_driver = mm_mod.p_amber_driver
mm_model.InitModel(AMOEBA)
#mm_model.SetTholeExponCoef(0.39)
mm_model.SetUseMortLib(0)
mm_model.SetPMECoef( 0.000015 )
mm_model.SetTholeExponCoef(0.0)
mm_model.SetDipoleScfTol(0.0001)
#mm_model.SetDipoleScfTol(0.01)
b1 = zm.GetCrdByTag("rMO")
b1.SetFrozen()
dist = 1.5
delt = 0.001
print("  dist     ene       f_z_dir      f_z_ene ")
while( dist < 5.01 ):
  b1.SetValue(dist)
  zm.SetAtomCrd()
  z0 = aptr1.GetZ()
  mm_mod.CalcEnergy()
  ene0= mm_mod.GetEne()
  f_z_dir = amber_driver.atm_frc.GetVal_idx0(3*1 + 2)
  aptr1.SetZ(z0 - 0.5*delt)
  mm_mod.CalcEnergy()
  ene1 = mm_mod.GetEne()
  aptr1.SetZ(z0 + 0.5*delt) 
  mm_mod.CalcEnergy()
  ene2 = mm_mod.GetEne()
  f_z_ene = -(ene2 - ene1)/delt
  print("%6.3f  %12.6f %9.3f  %9.3f " % (dist, ene0, f_z_dir, f_z_ene)) 
  dist = dist + 0.05
print("end")

