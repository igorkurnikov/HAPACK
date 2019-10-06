pmset = GetCurMolSet()
pmol = pmset.GetMoleculeNum(0)
p_phi = doublep()
p_cost = doublep()
p_psi   = doublep()
trans = Vec3D()
pmol.GetPosEulerTrans(p_phi, p_cost, p_psi, trans)
print("phi=", p_phi.value())
print("cost=", p_cost.value())
print("psi=", p_psi.value())
print("trans=", trans.GetX(), trans.GetY(), trans.GetZ())
p_phi.assign(1.01)
p_cost.assign(0.2)
p_psi.assign(1.03)
pmol.SetPosEulerTrans(p_phi.value(), p_cost.value(), p_psi.value(), trans)
pmset.RefreshAllViews(RFApply)



 
