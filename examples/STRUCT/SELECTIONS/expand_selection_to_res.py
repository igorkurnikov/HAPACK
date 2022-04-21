mset = GetCurMolSet()
ritr = ResidueIteratorMolSet(mset)
for r in ritr:
  if( r.HasSelectedAtoms() ):
    print(r.GetRef())
    r.SelectAtomsAll()

 
 