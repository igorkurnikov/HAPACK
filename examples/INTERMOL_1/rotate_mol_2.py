pmset = GetCurMolSet()
pmol = pmset.GetMoleculeNum(0)
aptr = pmset.GetAtomByRef("$GLY2$GLY2.C")
phi   = doublep()
cost = doublep()
psi    = doublep()
trans  = Vec3D()
pmol.GetPosEulerTrans(phi,cost,psi,trans)
phi_0   = -0.01
cost_0 = 0.99
psi_0   = 2.50201
rmat    = HaMat_double()
rmat.newsize(3,3)
Rot3D_EulerToRotMat(phi_0, cost_0 , psi_0, rmat )
delt = 0.005
rmat_delt    = HaMat_double()
rmat_delt.newsize(3,3)
Rot3D_EulerToRotMat(0.0, 0.999, 0.0, rmat_delt)
print rmat_delt.GetVal(1,1), rmat_delt.GetVal(1,2), rmat_delt.GetVal(1,3)
print rmat_delt.GetVal(2,1), rmat_delt.GetVal(2,2), rmat_delt.GetVal(2,3)
print rmat_delt.GetVal(3,1), rmat_delt.GetVal(3,2), rmat_delt.GetVal(3,3)
tmp = HaMat_double()
Rot3D_RotMatToEuler(rmat,phi,cost,psi)
for i in range(80):
  pmol.SetPosition(rmat,trans)
  print i*delt, aptr.GetX(), aptr.GetY(), aptr.GetZ()
#  print i*delt, phi.value(), cost.value(), psi.value()
  matmult(tmp,rmat,rmat_delt)
  matmult_T1(rmat,rmat_delt,tmp)
  Rot3D_RotMatToEuler(rmat,phi,cost,psi)
print "End of Script"
