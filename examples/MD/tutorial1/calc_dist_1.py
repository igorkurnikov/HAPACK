#
#	This script is used to calculate the distance between two atoms along the MD
#  trajectory after defining the two atoms.
#  It outputs the values in the file dist1.dat.
#

pmset = GetCurMolSet()
aptr1 = pmset.GetAtomByRef("TYR3.O")
aptr2 = pmset.GetAtomByRef("LYS2.NZ")
#  In GetAtomBy Ref you define the two atoms you want to measure the distance between

outf = open("dist1.dat","a")
dist = Vec3D_CalcDistance(aptr1,aptr2) #Calculates the distance
print(dist, file=outf)
outf.close()

