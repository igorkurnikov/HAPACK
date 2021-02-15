#
#  Script to superimpose two molecules aligning their sequences 
#  if they are in random order
#
from molset import *
outf=open("eps.out","w")
nm1 = []
nm2 = []
nm1.append("pept_1.pdb")
nm1.append("tyr_1.pdb")
nm2.append("pept_2.pdb")
nm2.append("tyr_2.pdb")
nn = len(nm1)
for i in range(nn):
  pmset = MolSet()
  pmset.FetchFile(FormatPDB,nm1[i])
  pmset.FetchFile(FormatPDB,nm2[i])
  pmol1 = pmset.GetMolByIdx(0)
  pmol2 = pmset.GetMolByIdx(1)
  aitr_m1 = AtomIteratorMolecule(pmol1)
  g1 = AtomGroup()
  aptr = aitr_m1.GetFirstAtom()
  while( aptr != None):
    g1.InsertAtom(aptr)
    aptr = aitr_m1.GetNextAtom()
  eps = pmset.AlignOverlapMol(g1,pmol2)
  print(" eps = ",eps)
  print(nm1[i], nm2[i], eps, file=outf)
outf.close()