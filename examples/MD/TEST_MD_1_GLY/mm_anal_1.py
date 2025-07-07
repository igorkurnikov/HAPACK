pmset = GetCurMolSet()
mm_mod = pmset.GetMolMechMod(1)
nn =  mm_mod.Dihedrals.size()
for i in range(nn):
  print mm_mod.Dihedrals[i].calc_14



