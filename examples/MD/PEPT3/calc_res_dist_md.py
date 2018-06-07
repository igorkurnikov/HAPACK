#
# Script to compute minimal distance between residue sidechain
# along MD trajectory
#
###########################################
if(script_status == SCRIPT_START): 
  mset = GetCurMolSet()
  res1 = mset.GetResByRef("LYS1")
  res2 = mset.GetResByRef("TYR2")
  bb_at_names = ["N","CA","C","O"]
  fout = open("res_dist.dat","w")
  idx_t = 0

elif(script_status == SCRIPT_STOP):
  fout.close()
  print "\nTrajectory Playback Completed\n"
  
else: 
  idx_t = idx_t + 1
  dist_min = 9999.9
  for aptr1 in res1:
    if( aptr1.GetName() in bb_at_names):
      continue
    for aptr2 in res2:
      if( aptr2.GetName() in bb_at_names):
        continue
      dist = Vec3D_CalcDistance(aptr1,aptr2)
      if( dist < dist_min):
        dist_min = dist
  print idx_t,dist_min
  print >> fout,"%6d %9.4f " % (idx_t,dist_min)




