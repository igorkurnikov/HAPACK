#
# This script reduce and fit MD trajectory for atom group 
#
if(script_status == SCRIPT_START): 
# Code to initialize the script on the trajectory

  iunit = 120
  mtraj_fname = "GLUR2_NOWAT_t49.1_59.1_nofit.mdcrd"
  mset = GetCurMolSet()
  grp_tr = mset.GetAtomGroupByID("PROT")
  
  MMDriverAmber.OpenAmberMDTrajFortran(iunit,mtraj_fname)

elif(script_status == SCRIPT_STOP):
#  Code to end the script after reading of MD trajectory is finished.
  
  print("\nEnd Trajectory Conversion\n")
  MMDriverAmber.CloseAmberMDTrajFortran(iunit);

else: 
#  Actual code that saves transformed trajectory snapshorts
  MMDriverAmber.WriteCrdToAmberMDTrajFortran( iunit, grp_tr );
  print("Snapshot saved to ", mtraj_fname)
	  