#
#  This script computes valence angle values 
#  along the MD trajectory 
#  Outputs is in ang1.dat file
#
if(script_status == SCRIPT_START): 
# Init the script
  
  pmset = GetCurMolSet()
  outf = open("ang1.dat","w")
  aptr1 = pmset.GetAtomByRef("TYR3.CE1")
  aptr2 = pmset.GetAtomByRef("TYR3.CZ")
  aptr3 = pmset.GetAtomByRef("TYR3.OH")
  i = 0 # initialize snapshot counter

elif(script_status == SCRIPT_STOP):
#  Code to end the script after MD trajectory is finished.
  
  outf.close()
  print("\nTrajectory Playback Completed\n")

else: 
#  Computations on snapshots:

    i = i+1 # Increments label by a value for each file name.
    ang = Vec3D_CalcAngle(aptr1,aptr2,aptr3)*RAD_TO_DEG #  Calculates the angle
    print("%5d %9.4f" % (i, ang), file=outf)
    


