#------------------------------------------------------------------------------
# gen_restr_heavy_atoms.py
#-------------------------------------------------------------------------------
# Generate Arbalest Resraints File (Type 1) for all heavy atoms of the protein  
#
# usage: harlem_nogui TYK2_EJM31.hin -script gen_restr_heavy_atoms.py  
#
fn_rules = "atom_restr_heavy_rules.xml"
fn_list  = "atom_restr_heavy_list.xml"
mset = GetCurMolSet()
aitr = AtomIteratorMolSet(mset)
aptr = aitr.GetFirstAtom()
while(aptr):
  if(aptr.GetElemNo() != 1): aptr.Select() 
  aptr = aitr.GetNextAtom()
atgrp = mset.SetAtomGroupFromSelection("HEAVY_ATOMS")
mm_mod = mset.GetMolMechMod(1)
mm_model = mm_mod.GetMolMechModel()
mm_model.SetRestrainedAtoms("HEAVY_ATOMS")
print("Saving Arbalest restraints to: ", fn_rules, " and ",fn_list)
mm_model.SaveAtomRestrArbalestIndForm(fn_rules,fn_list)
print("Done!")
quit()
