pmset = GetCurMolSet()
qcmod = pmset.GetQCMod(1)
#HaQCMod_SetQCIntEngine(INT_ENGINE_IPACK)
ovlp_mat = HaMat_double()
HaBasisSet_CalcOvlpMat(qcmod.AtBasis,qcmod.AtBasis,ovlp_mat)
nr = ovlp_mat.num_rows()
nc = ovlp_mat.num_cols()
print nr,nc
for i in range(nr):
  for  j in range(nc):
     print "%8.5f" % (ovlp_mat.GetVal(i+1,j+1)),
  print " "

