#
# This script to analyze MD trajectory
#
mset = HaMolSet()
mset.LoadHarlemFile("mg_wat_1.hlm")
traj_anal = mset.GetTrajAnalMod(1)
traj_anal.SetAmberMDCrdTraj("mg_wat_1.mdcrd")
traj_anal.SetPtStep(2)
traj_anal.traj_script = "save_crd_1.py"
traj_anal.AnalyzeTrajectory()


