mset = GetCurMolSet()
qc_mod = mset.GetQCMod(1)
zm = qc_mod.GetZMat()
b1 = zm.GetCrdByTag("roh1")
b1.SetFrozen()
qc_mod.InitBasis("3-21G")
qc_mod.SetEneMinCalc()
dist = 0.7
dist_arr = []
ene_arr = []
mset.DeleteCrdSnapshots()
snap_prefix = "SNAP_HF_321G_RMO_"
while( dist < 1.5 ):
  b1.SetValue(dist)
  prefix = "h2o_" + ( '%.2f' %  dist )
  snap_name = snap_prefix  + ("%0.2f" % dist)
  qc_mod.SetPrefix(prefix)
  qc_mod.Run()
  mset.AddCrdSnapshot(snap_name)
  ene = qc_mod.GetHFEne()
  print("ene = ",qc_mod.GetEne())
  print("dist = ", dist, " HF ene = ",ene)
  dist_arr.append(dist)
  ene_arr.append(ene)
  dist = dist + 0.05
n = len(dist_arr)
for i in range(n):
  print(dist_arr[i]," ",ene_arr[i])

