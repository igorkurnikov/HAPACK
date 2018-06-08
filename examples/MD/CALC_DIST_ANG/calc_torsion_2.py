#
# This script to compute dihedral angles along MD trajectory 
#
if(script_status == SCRIPT_START):
  pmset = GetCurMolSet()
  tors = []
  t1 = ["TYR3.CE1","TYR3.CZ","TYR3.OH","TYR3.HH"]
  tors.append(t1)
  outf = open("tors.dat","w")
  i = 0

elif(script_status == SCRIPT_STOP):
  outf.close()
  print "End dihedrals calculations"

else:
  i = i + 1
  print >> outf,"%8d" % i,
  for t in tors:
    aptr1 = pmset.GetAtomByRef(t[0])
    aptr2 = pmset.GetAtomByRef(t[1])
    aptr3 = pmset.GetAtomByRef(t[2])
    aptr4 = pmset.GetAtomByRef(t[3])
    tv = Vec3D_CalcTorsion(aptr1,aptr2,aptr3,aptr4)*RAD_TO_DEG
    print >>outf, " %12.4f " % tv, 
  print >> outf, "  "

