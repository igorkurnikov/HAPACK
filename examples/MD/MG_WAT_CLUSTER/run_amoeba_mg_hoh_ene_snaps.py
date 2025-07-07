mset = GetCurMolSet()
mm_mod = mset.GetMolMechMod(1)
mm_model = mm_mod.GetMolMechModel()
#mm_model.InitModel(AMBER_10)
mm_model.InitModel(AMOEBA)
#mm_model.SetTholeExponCoef(0.0)
mm_model.SetTholeExponCoef(0.39)
outf = open("mg_hoh7_snaps_amoeba_ene.dat","w")
for snap in mset.GetCrdSnapshots():
  sname = snap.GetName()
  snum   = int(sname[4:])
  mset.SetCrdFromSnapshot(snap)
  mset.per_bc.SetBox(100.0,100.0,100.0)
  mm_mod.CalcEnergy()
  print snum, " ", mm_mod.GetEne()
  print >> outf, snum, " ", mm_mod.GetEne()
outf.close()
print "End"



