from molset import *

mset = MolSet()
mset.LoadPDBFile("6PTI.pdb")   # Load Protein atomic coordinates

medit = MolEditor()
medit.DeleteExtraAtoms(mset)   # Delete aotoms not found in residue templates 
medit.AddMissingAtoms(mset)    # Add Hydrogens and other missing atoms from residue templates
medit.FixBondsUsingTempl(mset) # Fix bonds that may incorrectly be build when loading PDB
medit.Solvate(mset)            # Solvate molecule in water box 

# Energy minimization

mm_mod = mset.GetMolMechMod(True)
mm_mod.InitMolMechModel(AMBER_99_SB) # Set Protein Force field to ff98sb  
mm_mod.SetMMExternalProg(PMEMD_18)   # Set external MM simulation program to AMBER 18 pmemd binary
mm_mod.SetMMRunType(MIN_RUN)        
mm_mod.SaveAllInpFiles()
mm_mod.Run()                         # Run Energy minimization
mm_mod.UpdateMolInfo()               # Load Molecular coordinates from the AMBER restart file    

# Molecular Dynamics

mm_mod.SetMMRunType(MD_RUN)
mm_mod.SetNumMDSteps(1000)
mm_mod.SetWrtRstrtFreq(100)
mm_mod.SetWrtMDTrajFreq(100)
mm_mod.SetWrtLogFreq(10)
mm_mod.SaveAllInpFiles()
mm_mod.Run()                         # Run Molecular Dynamics 

# Compute RMSD of the trajectory 

