pmset = GetCurMolSet()
pmol = pmset.GetMoleculeNum(0)
aptr = pmset.GetAtomByRef("$GLY2$GLY2.C")
print "Old Atom Crd:", aptr.GetX(), aptr.GetY(), aptr.GetZ()
phi   = doublep()
cost = doublep()
psi    = doublep()
trans  = Vec3D()
pmol.GetPosEulerTrans(phi,cost,psi,trans)
phi_0   = -1.0
cost_0 = 0.99
psi_0   = 2.5
phi.assign(phi_0)
cost.assign(cost_0)
psi.assign(psi_0)
pmol.SetPosEulerTrans(phi.value(),cost.value(),psi.value(),trans)
rmat_0 = HaMat_double()
rmat_0.newsize(3,3)
rmat_1 = HaMat_double()
rmat_1.newsize(3,3)
rmat     = HaMat_double()
dmat    = HaMat_double()
tmp      = HaMat_double() 
Rot3D_EulerToRotMat(phi_0, cost_0, psi_0, rmat_0)
Rot3D_EulerToRotMat(phi_0, cost_0+0.001, psi_0, rmat_1)
matmult_T2(dmat,rmat_1,rmat_0)
rmat = rmat_0
print 0.0, aptr.GetX(), aptr.GetY(), aptr.GetZ(),phi.value(),cost.value(),psi.value()
for i in range(50):
  matmult(tmp,dmat,rmat)
  rmat = tmp
  pmol.SetPosition(rmat,trans)
  Rot3D_RotMatToEuler(rmat,phi,cost,psi)
  print (i+1)*0.001, aptr.GetX(), aptr.GetY(), aptr.GetZ(), phi.value(), cost.value(), psi.value()
print "End of Script"
