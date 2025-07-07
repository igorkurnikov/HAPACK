#
# This script reduce and fit MD trajectory for atom group 
#
if(script_status == SCRIPT_START): 
# Code to initialize the script on the trajectory

  iunit = 120
  mtraj_fname = "GLUR2_NOWAT_t19.1_89.1_fit_prot.mdcrd"
#  grp_name_fit = "CA_ATOMS"
  grp_name_fit = "PROT"
  grp_name_tr  = "PROT"
 
  mset = GetCurMolSet()
  grp_fit = mset.GetAtomGroupByID(grp_name_fit)
  grp_tr = mset.GetAtomGroupByID(grp_name_tr)
  crd_init = HaVec_double()
  grp_fit.SaveCrdToArray(crd_init)
  
  rot_mat = HaMat_double()
  trans_vec = HaVec_double()
  
  MMDriverAmber.OpenAmberMDTrajFortran(iunit,mtraj_fname)

elif(script_status == SCRIPT_STOP):
#  Code to end the script after reading of MD trajectory is finished.
  
  print "\nEnd Trajectory Conversion\n"
  MMDriverAmber.CloseAmberMDTrajFortran(iunit);

else: 
#  Actual code that saves transformed trajectory snapshorts
  eps =  PointContainer_GetSuperimposeMat( crd_init, grp_fit, rot_mat, trans_vec);
  if( eps > 0.0):
    print ipt_curr, "eps = ", eps
    grp_tr.Transform(rot_mat, trans_vec)
  else:
    print "Error to fit coordinates"
  MMDriverAmber.WriteCrdToAmberMDTrajFortran( iunit, grp_tr );
  print "Snapshot saved to ", mtraj_fname
	  