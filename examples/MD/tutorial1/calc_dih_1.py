#
#	This script is used to compute the values of a dihedral angle along the MD 
#  trajectory after defining the atoms involved in the angle.
#  It outputs these angle values in tor1.dat.
#

pmset = GetCurMolSet()
aptr1 = pmset.GetAtomByRef("TYR3.CE1")
aptr2 = pmset.GetAtomByRef("TYR3.CZ")
aptr3 = pmset.GetAtomByRef("TYR3.OH")
aptr4= pmset.GetAtomByRef("TYR3.HH")
#  In the GetAtomByRef command you select the atoms involved in the dihedral angle

outf = open("tor1.dat","a")
tor = Vec3D_CalcTorsion(aptr1,aptr2,aptr3,aptr4) #  Calculates the angle
print(tor, file=outf)
outf.close()


