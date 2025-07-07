mset1 = HaMolSet()
mset2 = HaMolSet()
mset1.LoadHarlemFile("HEM_FRAG1.hlm")
mset2.LoadHarlemFile("HEM_FRAG1_2.hlm")
crd1 = mset1.GetCrdArray()
rot_mat = HaMat_double()
trans_vec = HaVec_double()
eps = GetSuperimposeMat( crd1, mset2, rot_mat, trans_vec )
print " eps = ", eps
print " Rot Mat : "
for i in range(3):
  print " %9.4f %9.4f %9.4f " % ( rot_mat.GetVal_idx0(0,i),   rot_mat.GetVal_idx0(1,i),  rot_mat.GetVal_idx0(2,i) )
print " "
qc_mod1 = mset1.GetQCMod(1)
qc_mod2 = mset2.GetQCMod(1)
qc_mod1.InitBasis("3-21G_M1")
qc_mod2.InitBasis("3-21G_M1")
qc_mod1.LoadDataFromFChk("HEM_FRAG1.fchk")
qc_mod2.LoadDataFromFChk("HEM_FRAG1_2.fchk")
tr_mat = HaMat_double() 
qc_mod1.AtBasis.GetTransfMat( tr_mat, rot_mat )
mo_cf_tr = HaMat_double()
matmult( mo_cf_tr, tr_mat, qc_mod1.MO_coef )
nb = qc_mod1.MO_coef.num_rows()
imo = 85
for i in range(nb):
  print "%4d %9.5f %9.5f " % (i, mo_cf_tr.GetVal_idx0(i,imo-1), qc_mod2.MO_coef.GetVal_idx0(i,imo-1))
print "end"
  
  
  