#
# This script to compute G(r)  between Mg2+ and H2O along MD trajectory
#
mset = HaMolSet()
mset.LoadHarlemFile("mg_wat_1.hlm")
traj_anal = mset.GetTrajAnalMod(1)
traj_anal.SetAmberMDCrdTraj("mg_wat_1.mdcrd")
at_corr_ag =  traj_anal.GetAtomCorrAgent(1)
at_corr_ag.SetAtGroup1ByExpr("*.MG")
at_corr_ag.SetAtGroup2ByExpr("*.O")
at_corr_ag.SetDistRange(1.5,6.0,181)
traj_anal.SetWrapCrd(1)
traj_anal.SetPtEnd(10000)
traj_anal.AnalyzeTrajectory()



