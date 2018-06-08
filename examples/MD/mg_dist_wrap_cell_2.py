#
#	This script is to save coordinated oxygen atoms around MG ion  
#
if(script_status == SCRIPT_START): 
# Code to initialize the script on the trajectory
  mset = GetCurMolSet()
  editor = mset.GetMolEditor(1)
  editor.WrapToUnitCell(mset,mset.per_bc)
  i = 0 # Sets a label "i" beginning at 0.
  dist_cut = 3.0
  flig = open("lig_status.dat","w")
  atoms_ox = []
  status = []
  status_prev = []
  status_new = []
  for aptr in mset:
    if( aptr.GetName() == "MG"):
      at_ion = aptr
      print aptr.GetRef()
    if( aptr.GetName() == "O"):
      atoms_ox.append(aptr)
  idx = -1
  for atox in atoms_ox:
    idx = idx+1
    dist = Vec3D_CalcDistance(at_ion,atox)
    if( dist < dist_cut ): 
      status.append(idx)
      status_prev.append(idx)

elif(script_status == SCRIPT_STOP):
#  Code to end the script after MD trajectory is finished.
  
  print " N water = ",len(atoms_ox)
  print " N_lig_last =",len(status)
  print " N_lig_prev =",len(status_prev)
  flig.close()
  print "\nTrajectory Playback Completed\n"
  
else: 
  i = i+1 # Increments label by a value for each file name.
  editor.WrapToUnitCell(mset,mset.per_bc)
    
  status_new = []
  if(i % 10 == 0):
    idx = -1
    for atox in atoms_ox:
      idx = idx + 1
      dist = Vec3D_CalcDistance(at_ion,atox)
      if (dist < dist_cut):
        status_new.append(idx)
    diff_prev_last = 0
    nlig_prev = len(status_prev)
    if( len(status) != nlig_prev ):
      diff_prev_last = 1
    else:
      for m in range(0,nlig_prev):
        if( status_prev[m] != status[m] ):
          diff_prev_last = 1
    diff_prev_new  = 0  
    if( len(status_new) != nlig_prev ):
      diff_prev_new = 1
    else:
      for m in range(0,nlig_prev):
        if( status_prev[m] != status_new[m] ):
          diff_prev_new = 1
    switch_lig = 0
    if( diff_prev_last == 1 ):
      switch_lig = 1
    if( diff_prev_new == 0 ):
      switch_lig = 0
    
    print >> flig, "%3d   " % switch_lig,
    for idx in  status:
      print >> flig, "%5d " % idx,
    print >> flig, "   "
    status_prev = status
    status      = status_new

