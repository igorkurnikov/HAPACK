# Setup ARROW FF v.5.14 atom types on the protein   
# and save to HIN file
#
# example usage: harlem_nogui THROMBIN.pdb -script typify_prot_arrow.py 
#
# OUTPUT: protein_arrow,hin
#
mset = GetCurMolSet()
mol_editor = mset.GetMolEditor(1)
mol_editor.RenameAtomsToAmber(mset)
mol_editor.FixBondsUsingTempl(mset)
mm_mod = mset.GetMolMechMod(1)
mm_mod.InitMolMechModel(ARROW_5_14_CT)
print("Saving ARROW typified file protein_arrow.hin")
mset.SaveHINFile("protein_arrow.hin")
print("Done")
quit()


