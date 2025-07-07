#
#  Script to superimpose two molecules aligning their sequences 
#  if they are in random order
#
outf=open("eps.out","w")
nm1 = []
nm2 = []
nm1.append("pept_1.pdb")
nm1.append("tyr_1.pdb")
nm2.append("pept_2.pdb")
nm2.append("tyr_2.pdb")
nn = len(nm1)
for i in range(nn):
  pmset = HaMolSet()
  pmset.FetchFile(FormatPDB,nm1[i])
  pmset.FetchFile(FormatPDB,nm2[i])
  pmol1 = pmset.GetMoleculeNum(0)
  pmol2 = pmset.GetMoleculeNum(1)
  print pmol1, pmol2
  aitr_m1 = AtomIteratorMolecule(pmol1)
  atl1 = AtomList()
  aptr = aitr_m1.GetFirstAtom()
  while( aptr != None):
    atl1.InsertAtom(aptr)
    aptr = aitr_m1.GetNextAtom()
  eps = pmset.AlignOverlapMol(atl1,pmol2)
  print " eps = ",eps
  print >> outf, nm1[i], nm2[i], eps
outf.close()