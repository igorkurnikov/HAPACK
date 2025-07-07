mset = GetCurMolSet()
mm_mod = mset.GetMolMechMod(1)
mm_model = mm_mod.GetMolMechModel()
mm_model.InitModel(AMOEBA)
mm_mod.CalcEnergy()
dip = mm_model.GetTotIndDipole()
dip.Scale(EL_ANG_TO_DEBYE)
dip1 = mm_model.GetTotIndDipole1()
dip1.Scale(EL_ANG_TO_DEBYE)
dip2 = mm_model.GetTotIndDipole2()
dip2.Scale(EL_ANG_TO_DEBYE)
print " tot ind dipole = %10.6f %10.6f %10.6f "  % (dip.GetX(), dip.GetY(), dip.GetZ())
print " tot ind dipole1 = %10.6f %10.6f %10.6f " % (dip1.GetX(), dip1.GetY(), dip1.GetZ())
print " tot ind dipole2 = %10.6f %10.6f %10.6f " % (dip2.GetX(), dip2.GetY(), dip2.GetZ())




