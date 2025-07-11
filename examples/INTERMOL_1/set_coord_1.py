pmset = GetCurMolSet()
mol1 = pmset.GetMoleculeNum(0)
mol2 =  pmset.GetMoleculeNum(1)
trans1 = Vec3D()
trans2 = Vec3D()
trans1.SetX(-42.980601746)
trans1.SetY(14.605993063)
trans1.SetZ(34.303955895)
trans2.SetX(1.220677444)
trans2.SetY(54.426379915)
trans2.SetZ(42.317604249)
mol1.SetPosEulerTrans( 1.036376496,0.891240042,3.099665229, trans1 )
mol2.SetPosEulerTrans( 3.062027926,0.692409655,-3.133207412, trans2 )
pmset.RefreshAllViews(RFApply)
