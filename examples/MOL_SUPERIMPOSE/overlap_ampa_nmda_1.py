#
#  Script to superimpose NMDA model and AMPA 
#  over CA atoms of M1 and M3 helicies of 4 subunits
#
mset = GetCurMolSet()
grp_ampa = mset.GetAtomGroupByID("M1_M3_CA_AMPA")
grp_nmda = mset.GetAtomGroupByID("M1_M3_CA_NMDA")
mset.OverlapMol(grp_ampa,grp_nmda)
mset.RefreshAllViews(RFApply)