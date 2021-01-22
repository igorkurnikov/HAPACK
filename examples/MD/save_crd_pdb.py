#
#   This script is used to save Harlem file format snapshots (.hlm) along the MD 
#  trajectory. These files are saved in the same directory as your molecule files. 
#

if(script_status == SCRIPT_START): 
# Code to initialize the script on the trajectory
  
  pmset = GetCurMolSet()
  i = 0 # Sets a label "i" beginning at 0.
  

elif(script_status == SCRIPT_STOP):
#  Code to end the script after MD trajectory is finished.
  
  print("\nTrajectory Playback Completed\n")
  print("PDB snapshots saved to file directory")


else: 
#  Actual code that saves the PDB files.

    i = i+1 # Increments label by a value for each file name.
    
    if(i % 10 == 0):
      # The value used with the modulus for "i" determines how often to take a 
      #  snapshot. Here it is set to take a snapshot every 100 MD steps.
      #  Change this value as desired.
          
      fname = "snap" + str(10000+i) + ".pdb" #  Sets the file name of the snapshot
      print("Snapshot saved to ", fname)
      pmset.SavePDBFile(fname) # Writes the PDB file