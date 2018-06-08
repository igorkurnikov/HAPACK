#
# This script allows you to restrict the trajectory view to only the first 
#  molecule in the harlem molecule list. This is useful for example when you only 
#  want to map the trajectory of a molecule that has been solvated in a box of water.
#
# The trajectory coordinates are saved to the file "red_traj.mdcrd" in your
#  file directory.
#
mset = HaMolSet()
mset.LoadHarlemFile("methane_wat_1.hlm")
mol = mset.GetMoleculeNum(0) #  Selects the first molecule in the system
traj_anal = mset.GetTrajAnalMod(1)
traj_anal.SetAmberMDCrdTraj("methane_wat_1.mdcrd")
traj_anal.SetPtBegin(10)
traj_anal.SetPtEnd(50)
traj_anal.SetPtStep(5)
traj_anal.SetAlignToFstPt(0)
traj_anal.SetWrapCrd(1)
traj_anal.ReduceAmberMDTraj("red_traj.mdcrd",mol)


