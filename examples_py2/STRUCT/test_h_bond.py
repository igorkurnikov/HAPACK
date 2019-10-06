pmset = GetCurMolSet()
at1 = pmset.GetAtomByRef(ASP52:2.N)
at2 = pmset.GetAtomByRef(ASN44:2.O)
ires = pmset.AreHBonded(at1,at2)
