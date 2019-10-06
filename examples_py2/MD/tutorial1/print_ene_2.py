#
#  	This script prints the total energy, potential energy, kinetic energy, and temperature
#  along the steps of the MD trajectory.  
#  It reads this information from  the ***.mden file, and outputs to the file ene_2.dat.
#  
#  You can find an example of formatting output with the PYTHON %  operator to create lists. 
#  The increment of a variable 't' tracks MD time.
#

if(script_status == SCRIPT_START): #  Script initialization
  pmset = GetCurMolSet()
  mm_mod = pmset.GetMolMechMod(1) #  Selects first molecule of the system
  mm_info = mm_mod.p_mm_info
  outf = open("ene_2.dat","w")
  t = 0.0 # Time variable. Starts at 0.
  
  print >>outf, " time   tot_ene    pot_ene    kin_ene  temperature"
  #  Prints out labels for the lists.

elif(script_status == SCRIPT_STOP): 
#  When MD trajectory analysis is complete the output file is closed.
  print "\nTotal Energy, Potential Energy, Kinetic Energy,"
  print "and Temperature Data saved to ene_2.dat"
  outf.close()
  

else:
  print >>outf, "%6.3f %9.3f %9.3f %9.3f %9.3f" %  (t, mm_info.tot_energy, mm_info.pot_ene, mm_info.kin_ene, mm_info.temp)
  #  Prints out the value while placing spaces in between each value to create neat lists
  
  t = t + 0.01 #  Increments time at each step


