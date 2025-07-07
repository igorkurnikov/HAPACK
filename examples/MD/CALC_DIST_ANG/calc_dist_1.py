#
#  This script computes atom atom distances 
#  along the MD trajectory 
#  Output is in dist1.dat file
#
if(script_status == SCRIPT_START): 
# Init the script
  
  pmset = GetCurMolSet()
  outf = open("dist1.dat","w")
  aptr1 = pmset.GetAtomByRef("TYR3.CE1")
  aptr2 = pmset.GetAtomByRef("TYR3.CZ")

elif(script_status == SCRIPT_STOP):
#  Code to end the script after MD trajectory is finished.
  
  outf.close()
  print "\nTrajectory Playback Completed\n"

else: 
#  Computations on snapshots:

    dist = Vec3D_CalcDistance(aptr1,aptr2) # Calculates atom-atom distance (Ang)
    print >>outf, "%5d %9.4f" % (idx_curr_pt, dist)


