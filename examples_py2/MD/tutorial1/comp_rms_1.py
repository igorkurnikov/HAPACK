#
#	This script computes the RMS deviation of the first molecule of the system
#  along the steps of the MD trajectory compared to its initial coordinates. 
#  It saves these deviation values to the file rms_1.dat.
#  

if(script_status == SCRIPT_START):
#	Initial part of the script - opens file rms_1.dat and sets up lists of atoms 
#  with initial coordinates and current coordinates of the first molecule.

  print "Script start"
  pmset = GetCurMolSet()
  pmol = pmset.GetFirstMolecule()
  atl1 = AtomGroup()
  atl2 = AtomGroup()
  rot_mat = HaMat_double()
  rot_mat.newsize(3,3)
  trans_vec = HaVec_double()
  trans_vec.newsize(3)
  p_eps = doublep()
  aitr = AtomIteratorMolecule(pmol)
  aptr = aitr.GetFirstAtom()
  arr_at1 = []  
  
  while aptr != None:            # A cycle on atoms of the first molecule
    new_at = HaAtom()
    new_at.SetX( aptr.GetX())
    new_at.SetY( aptr.GetY())
    new_at.SetZ( aptr.GetZ())
    atl1.InsertAtom(new_at)	 # Sets up a spot for a new value	
    arr_at1.append(new_at)    # This line is needed so that PYTHON 
					 # will not delete "new_at" automatically
    atl2.InsertAtom(aptr)
    aptr = aitr.GetNextAtom()
  outf = open("rms_1.dat","w")
  t = 0.0
  print >>outf, " time   rms"


elif(script_status == SCRIPT_STOP):
# When MD trajectory analysis is complete this command closes the output file
  
  print "\nRMS Deviation Data saved to rms_1.dat"
  outf.close()
  
  
else:
#  	Main part of the script: 
#  This is where RMS deviation values in the points of MD trajectory are computed

  PointCollection_GetSuperimposeMat( atl1, atl2, rot_mat, trans_vec, p_eps)
  rms = p_eps.value()
  print >>outf, "%6.3f %9.3f " %  (t, rms)
  t = t + 0.01        
#  The value used with t is the time step used between MD snapshots. This value
#  should be changed depending on the time step.
