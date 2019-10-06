pmset = GetCurMolSet()
pmol = pmset.GetMoleculeNum(0)
aptr = pmset.GetAtomByRef("$GLY2$GLY2.C")
print "Old Atom Crd:", aptr.GetX(), aptr.GetY(), aptr.GetZ()
phi   = doublep()
cost = doublep()
psi    = doublep()
trans  = Vec3D()
pmol.GetPosEulerTrans(phi,cost,psi,trans)
phi_0   = -0.5
cost_0 = 0.99
psi_0   = 2.0
phi.assign(phi_0)
cost.assign(cost_0)
psi.assign(psi_0)
pmol.SetPosEulerTrans(phi.value(),cost.value(),psi.value(),trans)
print "phi=", phi.value()
print "cost=", cost.value()
print "psi=", psi.value()
print trans.GetX(), trans.GetY(),trans.GetZ()
print "Atom Crd:", aptr.GetX(), aptr.GetY(), aptr.GetZ()
for i in range(20):
  delt = 0.001*i
  phi.assign(phi_0)
  cost.assign(cost_0) 
  psi.assign(psi_0)
  Rot3D_IncrEulerAng(phi,cost,psi,-0.0,delt, 0.0)
  pmol.SetPosEulerTrans(phi.value(),cost.value(),psi.value(),trans)
  print delt, aptr.GetX(), aptr.GetY(), aptr.GetZ(), phi.value(), cost.value(), psi.value()
print "End of Script"
