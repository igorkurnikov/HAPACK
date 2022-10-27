import os
from molset import *

cdir = os.getcwd()

mset = MolSet()
mset.LoadPDBFile("6PTI.pdb")
medit = MolEditor()

medit.AddMissingAtoms(mset)  # Add Hydrogens and other missing atoms from residue templates
medit.Solvate(mset)          # Solvate molecule in water box 

# Energy minimization

mm_mod = mset.GetMolMechMod(True)
mm_mod.InitMolMechModel(AMBER_99_SB)
mm_mod.SetMMExternalProg(PMEMD_18)
mm_mod.SetMMRunType(MIN_RUN)
mm_mod.p_amber_driver.SaveAllInpFiles()
mm_mod.Run()
mm_mod.p_amber_driver.LoadAmberRestartFile("MOLSET.rst")

# Molecular Dynamics

mm_mod.SetMMRunType(MD_RUN)
#mm_mod.SetRunInternal(0)
mm_model = mm_mod.GetMolMechModel()
mm_model.SetNBCutDist(9.0)
mm_mod.SetNumMDSteps(1000)
mm_mod.SetWrtRstrtFreq(100)
mm_mod.SetWrtMDTrajFreq(100)
mm_mod.SetWrtLogFreq(10)
mm_mod.p_amber_driver.SaveAllInpFiles()
mm_mod.Run()

# Compute RMSD of the trajectory 

