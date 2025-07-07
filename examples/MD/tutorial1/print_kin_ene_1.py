#
#	This script computes the kinetic energy of the first molecule of the system 
#  along steps of the MD trajectory.
#  It outputs these values in the file kin_ene_1.dat.
#

pmset = GetCurMolSet()
mm_mod = pmset.GetMolMechMod(1) #  Selects the first molecule of the system
mm_info = mm_mod.p_mm_info
outf = open("kin_ene_1.dat","a")
print >>outf,  mm_info.kin_ene #  Prints the kinetic energy of the molecule
outf.close()

