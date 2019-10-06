#
# This script center and warp molecules along MD trajectory 
#

if(script_status == SCRIPT_START): 
# Code to initialize the script on the trajectory

  iunit = 120
  mtraj_fname = "traj_mod.mdcrd"
#  grp_name = "CA_ATOMS"
  grp_name = "GLF_CONTACT_GRP"
  
  mset = GetCurMolSet()
  grp = mset.GetAtomGroupByID(grp_name)
  cnt = grp.GetAverageCoord()
  cnt.SetX(35.0)
  cnt.SetY(35.0)
  cnt.SetZ(35.0)
  i = 0 # Sets counter "i" beginning at 0.
  
  MMDriverAmber.OpenAmberMDTrajFortran(iunit,mtraj_fname)

elif(script_status == SCRIPT_STOP):
#  Code to end the script after reading of MD trajectory is finished.
  
  print "\nEnd Trajectory Conversion\n"
  MMDriverAmber.CloseAmberMDTrajFortran(iunit);

else: 
#  Actual code that saves transformed trajectory snapshorts

  if(i % 1 == 0):  # The value used with the modulus for "i" determines how often to save snapshots
    mset.WrapAndCenter(grp_name,cnt)
    MMDriverAmber.WriteCrdToAmberMDTrajFortran( iunit, mset );
    print "Snapshot saved to ", mtraj_fname
  i= i+1
	  