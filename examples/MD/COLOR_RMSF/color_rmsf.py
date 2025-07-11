mset = GetCurMolSet()
grp = mset.GetAtomGroupByID("ALL_CA")
print(grp.GetNAtoms())
aitr = AtomIteratorAtomGroup(grp)
finp = open("rmsf_atom_avg_MD.out","r")
aptr = aitr.GetFirstAtom()
while( aptr):
  str = finp.readline()
  str = str.strip()
  words = str.split()
  res = aptr.GetHostRes()
  aitr_r = AtomIteratorAtomGroup(res)
  aptr_r = aitr_r.GetFirstAtom()
  while(aptr_r):
    aptr_r.tempf =  float(words[1])
    aptr_r = aitr_r.GetNextAtom()
  print(aptr.GetRef(),words[0], aptr.tempf)
  aptr = aitr.GetNextAtom()
finp.close()
