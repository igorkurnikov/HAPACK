#
# Calculate h-bonds with loose criteria (6.0 Ang between donor/acceptor)
# and compute averages along MD trajectory
#
import math
from math import *
import molsetc
from molsetc import *
pmset = GetCurMolSet()
script_status = cvar.script_status
# to execute start part 
#script_status = SCRIPT_CONTINUE
#script_status = SCRIPT_START
#script_status = SCRIPT_STOP
if(script_status == SCRIPT_START):
  outf = open("hydroph_pdz1_NDSLL_44ns.dat","w")
#tanya test___________________________  
  #outf2 = open("list_hbonds.dat","w") 
# the end of test_______________________  
  atmitr=AtomIteratorMolSet(pmset)
  atm_tr = atmitr.GetFirstAtom()
  atnm_map = {}  # create map of atom references to atom sequence numbers
  i = 0
  while(atm_tr != None):
    i= i+1
    ref = atm_tr.GetRef()
    atnm_map[ref] = i 
    atm_tr = atmitr.GetNextAtom()
  atl = AtomList()
  d_list = []  #  donor atoms of possible hydrophobic bonds
  a_list = []  #  acceptor atoms of possible hydrophobic bonds
  i = 0
  a_arr = [] #array for acceptors involved to the hydrofobic thether
  res_arr = []
  atm_tr = atmitr.GetFirstAtom()
  while(atm_tr != None):  
#  at_num = atm_tr.GetSerNo()
    i = i + 1
    print i
#    if(i > 5000): break
    at_name = atm_tr.GetName()
#    at_don  = atm_tr.IsHBDonor()
    at_don = atm_tr.GetElemNo()
    res = atm_tr.GetHostRes()
#    print "B"
    if(atm_tr.GetElemNo() == 1 or atm_tr.GetElemNo() == 7 or atm_tr.GetElemNo() == 8):
#     if(at_name =="CA" or at_name =="C"):
       atm_tr = atmitr.GetNextAtom()
       continue
    if(at_name =="CA" or at_name =="C"): 
       atm_tr = atmitr.GetNextAtom()
       continue     
    don_ref = atm_tr.GetRef()
    print don_ref
    a_arr.append(atm_tr)
    res_arr.append(res) 
    atm_tr = atmitr.GetNextAtom()
  na = len(a_arr)
  print "na =", na
  atm_tr = atmitr.GetFirstAtom()
  for i in range(na): 
       for j in range(i+1,na):
          res1 = HaAtom_GetHostRes(a_arr[i])
	  res2 = HaAtom_GetHostRes(a_arr[j])
	  res1_a = MatPoint_GetRef(a_arr[i]) 
	  res2_a = MatPoint_GetRef(a_arr[j])
	  dist = Vec3D_CalcDistance(a_arr[i],a_arr[j], ANGSTROM_U) 
#	  if( res_arr[i] == res_arr[j] and dist >6.0):	  #Igor'
	  if( res1== res2):
	     continue
	  if (dist>6.0):
	      continue
#	  if( res1 == res2): 
#	    continue
#	  if (dist>6.0):
#             continue 
#	  dist = Vec3D_CalcDistance(a_arr[i],a_arr[j], ANGSTROM_U)   
          d_list.append(a_arr[i])         # add donor atom of the possible H-bond  
          a_list.append(a_arr[j]) # add acceptor atom of the possible H-bond 
#          print "TEST"	 
  nn = len(d_list)   # number of possible hydrophobic bonds
  print "N_HB_possible = ",nn
  dist1_arr = []      # accumulate sum of distances between donor and acceptor 
  dist1_2_arr = []   # accumulate sum of squares of distances between donor and acceptor 
  on_arr = []      # number of MD snapshots this h-bond is ON for normal geometric criteria
  n_md_pt = 0    # number of MD-points 
  for i in range(nn):       # initialize arrays for average distances and angles
    dist1_arr.append(0.0)
    dist1_2_arr.append(0.0)
    on_arr.append(0)
elif(script_status == SCRIPT_STOP):
  nn = len(d_list)
  for i in range(nn):
    n_on = on_arr[i]
    if( n_on > 0 and n_md_pt > 0):
      dist1_av = dist1_arr[i]/n_md_pt
      dist1_rms =  sqrt(dist1_2_arr[i]/n_md_pt - dist1_av*dist1_av)
      dref = d_list[i].GetRef()
      aref = a_list[i].GetRef()
      p_on = float(n_on)/float(n_md_pt)  # compute fraction of MD snapshots H-bond is on, float convert integer to real
      idx_d  = atnm_map[dref]
      idx_a  = atnm_map[aref]
      print >> outf, "%20.20s  %20.20s" %  (dref, aref),  
      print >> outf, "%5d %5d" %  (idx_d, idx_a), 
      print >> outf, " %5d %6.3f %8.3f %8.3f " %  (n_on, p_on,dist1_av, dist1_rms)
  outf.close()
else:
  n_md_pt = n_md_pt + 1
  nn = len(d_list)  # number of possible H-bonds
  print "point ", n_md_pt
  for i in range(nn):   # cycle over possible H-bonds
    dist1 = Vec3D_CalcDistance(a_list[i],d_list[i])*BOHR_TO_ANG
    dist1_arr[i] = dist1_arr[i] + dist1
    dist1_2_arr[i] = dist1_2_arr[i] + dist1*dist1
    if( dist1 < 4.0):
      on_arr[i] = on_arr[i] + 1   
  print "end_script"
  
