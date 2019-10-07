pmset = GetCurMolSet()
qcmod = pmset.GetQCMod(1)
cvar.HaQCMod_int_engine = INT_ENGINE_GAUSS
ovlp_mat = HaMat_double()
r1 = HaOperR()
rmats = HaMat_doubleArr()
r1.EvalGauBasisSet(qcmod.AtBasis,rmats)
HaQCMod_CalcOvlpMat(qcmod.AtBasis,qcmod.AtBasis,ovlp_mat)
nr = ovlp_mat.num_rows()
nc = ovlp_mat.num_cols()
for i in range(nr):
  for  j in range(nc):
     print("%8.5f" % (ovlp_mat.GetVal(i+1,j+1)), end=' ')
  print(" ")

