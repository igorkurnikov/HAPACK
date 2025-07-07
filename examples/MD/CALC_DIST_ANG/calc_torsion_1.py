#
#  This script computes torsional angle values 
#  along the MD trajectory 
#  Outputs is in tor1.dat file
#
if(script_status == SCRIPT_START): 
# Init the script
  
  pmset = GetCurMolSet()
  outf = open("tor1.dat","w")
  aptr1 = pmset.GetAtomByRef("TYR3.CE1")
  aptr2 = pmset.GetAtomByRef("TYR3.CZ")
  aptr3 = pmset.GetAtomByRef("TYR3.OH")
  aptr4= pmset.GetAtomByRef("TYR3.HH")
  i = 0 # initialize snapshot counter

elif(script_status == SCRIPT_STOP):
#  Code to end the script after MD trajectory is finished.
  
  outf.close()
  print "\nTrajectory Playback Completed\n"

else: 
#  Computations on snapshots:

    i = i+1 
    tor = Vec3D_CalcTorsion(aptr1,aptr2,aptr3,aptr4)*RAD_TO_DEG #  Calculates the torsional angle
    print >>outf, "%5d %9.4f" % (i, tor)


