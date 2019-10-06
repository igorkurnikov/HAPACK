#	
#	This script is used to compute the total energy of the first molecule in Harlem
#  molecule list along the MD trajectory. It outputs these values in the file
#  ene_1.dat
#

pmset = GetCurMolSet()
mm_mod = pmset.GetMolMechMod(1) # Selects the first molecule in the system
mm_info = mm_mod.p_mm_info
outf = open("ene_1.dat","a")
print(mm_info.tot_energy, file=outf)  #  Prints the current energy of a step
outf.close()
