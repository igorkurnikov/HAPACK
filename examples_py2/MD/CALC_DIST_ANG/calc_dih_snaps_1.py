#
# This script calculates dihedrals for crd snapshots 
#
mset = GetCurMolSet()
tors = []
t1 = ["$CYT$TYR67:A.CA","$CYT$TYR67:A.CB","$CYT$TYR67:A.CG","$CYT$TYR67:A.CD2"]
tors.append(t1)
t2 = ["$CYT$RUI119:B.C03","$CYT$RUI119:B.C05","$CYT$RUI119:B.N46","$CYT$RUI119:B.C47"]
tors.append(t2)
t3 = ["$CYT$CYS66:A.SG","$CYT$RUI119:B.C49","$CYT$RUI119:B.C47","$CYT$RUI119:B.N46"]
tors.append(t3)
outf = open("tors.dat","w")
for snap in mset.GetCrdSnapshots():
  mset.SetCrdFromSnapshot(snap)
  sname = snap.GetName()
  snum   = int(sname[4:])
  if( snum % 10 !=  0):
    continue
  print "%8d" % snum,
  print >> outf,"%8d" % snum,
  for t in tors:
    aptr1 = mset.GetAtomByRef(t[0])
    aptr2 = mset.GetAtomByRef(t[1])
    aptr3 = mset.GetAtomByRef(t[2])
    aptr4 = mset.GetAtomByRef(t[3])
    tv = Vec3D_CalcTorsion(aptr1,aptr2,aptr3,aptr4)*RAD_TO_DEG
    print " %12.4f " % tv, 
    print >>outf, " %12.4f " % tv, 
  print "  "
  print >> outf, "  "
outf.close()

