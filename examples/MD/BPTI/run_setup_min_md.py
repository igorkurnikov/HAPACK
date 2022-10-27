from molset import *

mset = MolSet()
mset.LoadPDBFile("6PTI.pdb")   # Load Protein atomic coordinates

medit = MolEditor()
medit.DeleteExtraAtoms(mset)   # Delete atoms not found in residue templates 
medit.AddMissingAtoms(mset)    # Add Hydrogens and other missing atoms from residue templates
medit.FixBondsUsingTempl(mset) # Fix bonds that may incorrectly be build when loading PDB
medit.Solvate(mset)            # Solvate molecule in water box 

# Compute pKa using Continuum Electrostatics model and set atomic charges acccording to pH:

prot_mod = mset.GetProtonRedoxMod(True)
pApp = GetHarlemApp()
pApp.RasMolCmd("select protein")
prot_mod.CalcPKaForSelection(True)
prot_mod.SetPH(6.0)
prot_mod.SetChargesForCurrentPH()

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

# Compute RMSD of CA atoms along the MD trajectory 

trj_anal_mod = mm_mod.GetTrajAnalMod()       # Get MD trajectory analysis module
rmsd_ag = trj_anal_mod.GetRMSDAgent(True);   # Get an Agent  to compute RMSD along MD trajectory
pApp.RasMolCmd("define CA_ATOMS *.CA")       # define Atom Group with C-alpha atoms of the protein
rmsd_ag.SetAtomsFit("CA_ATOMS");
rmsd_ag.SetAtomsRMSD("CA_ATOMS")
trj_anal_mod.AnalyzeTrajectory(True)         # Analyse MD trajectory computing RMSD of CA_ATOMS 




