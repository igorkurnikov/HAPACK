#-----------------------------------------------------------------------------
# gen_restr_CA_9ang.py
#-------------------------------------------------------------------------------
# Generate Arbalest Resraints File (Type 1) for protein CA atoms for 
# residues with closet distance to the ligand > 9 Ang
#  
# usage: harlem_nogui TYK2_EJM31.hin --script gen_restr_CA_9ang.py  
#
fn_rules = "atom_restr_ca_9ang_rules.xml"
fn_list  = "atom_restr_ca_9ang_list.xml"
mset = GetCurMolSet()
pApp.RasMolCmd("select *.ca")
atgrp_ca = mset.SetAtomGroupFromSelection("CA_9ANG")
pApp.RasMolCmd("select within(9.0,LIG)")
del_atoms= []
aitr = AtomIteratorAtomGroup(atgrp_ca)
aptr = aitr.GetFirstAtom()
while(aptr):
  res = aptr.GetHostRes()
  if(res.HasSelectedAtoms()): del_atoms.append(aptr)
  aptr = aitr.GetNextAtom()
for aptr in del_atoms:
  atgrp_ca.DeleteAtom(aptr)

mm_mod = mset.GetMolMechMod(1)
mm_model = mm_mod.GetMolMechModel()
mm_model.SetRestrainedAtoms("CA_9ANG")
print("Saving Arbalest restraints to: ", fn_rules, " and ",fn_list)
mm_model.SaveAtomRestrArbalestIndForm(fn_rules,fn_list)
print("Done!")
quit()
