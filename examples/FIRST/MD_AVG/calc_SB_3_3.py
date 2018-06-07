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
  outf = open("SB_pdz1_NDSLL_44ns.dat","w")
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
  d_list = []  #  donor atoms of possible H-bonds
  h_list = []  #  hydrogen atoms of possible H-bonds
  a_list = []  #  acceptor atoms of possible H-bonds
  i = 0
  atm_tr = atmitr.GetFirstAtom()
  while(atm_tr != None):  
#  at_num = atm_tr.GetSerNo()
    i = i + 1
    print i
#    if(i > 100): break
    at_name = atm_tr.GetName()
    at_don  = atm_tr.IsHBDonor()
    res = atm_tr.GetHostRes()
    res_name_d = res.GetName()
    if(not at_don or res_name_d == "HOH" or atm_tr.GetElemNo() == 1):
      atm_tr = atmitr.GetNextAtom()
      continue
    if (res_name_d !="LYS" and res_name_d !="ARG" and res_name_d !="HIS"):
      atm_tr = atmitr.GetNextAtom()
      continue
    if(at_name =="NE" or at_name =="NH1" or at_name =="NH2" or at_name =="NZ" or at_name =="ND1"):
 #     atm_tr = atmitr.GetNextAtom()
  #    continue    
#      
      don_ref = atm_tr.GetRef()
      atm_tr.GetBondedAtoms(atl) # Get List of atoms bonded to atm_tr
      h_arr = []                          # array for hydrogens connected to donor
      aptr = atl.GetFirstAtom()
      while( aptr != None):
         if(aptr.GetElemNo() == 1):
          h_arr.append(aptr)          # add hydrogen to array 
         aptr = atl.GetNextAtom()
      atmitr_acc=AtomIteratorMolSet(pmset)
      atm_tr_acceptor = AtomIteratorMolSet_GetFirstAtom(atmitr_acc)
      while(atm_tr_acceptor !=None):
       res = HaAtom_GetHostRes(atm_tr_acceptor)
       res_name = HaResidue_GetName(res)   
       atom_acceptor = HaAtom_GetName(atm_tr_acceptor)       
       if(not HaAtom_IsHBAcceptor(atm_tr_acceptor) or res_name == "HOH"):
         atm_tr_acceptor = AtomIteratorMolSet_GetNextAtom(atmitr_acc)
	 continue
       if (atom_acceptor == "OE1" or atom_acceptor == "OE2" or atom_acceptor == "OD1" or atom_acceptor == "OD2"):
        dist = Vec3D_CalcDistance(atm_tr_acceptor,atm_tr, ANGSTROM_U)
        if( dist > 6.0):
         atm_tr_acceptor = AtomIteratorMolSet_GetNextAtom(atmitr_acc)
	 continue
        nh = len(h_arr)        # number of of hydrogens attached to the donor 
        for ih in range(nh):  
          atm_tr_h = h_arr[ih]        # pointer to hydrogen in the array of hydrogens
	  href = atm_tr_h.GetRef()  # reference(name) of the hydrogen
	  acc_ref = MatPoint_GetRef(atm_tr_acceptor) 
#	  print don_ref, href, acc_ref, " -", atnm_map[don_ref],atnm_map[href],atnm_map[acc_ref]
#   Add  possible hydrogen bond (with  donor - acceptor distance < 6.0)
          d_list.append(atm_tr)         # add donor atom of the possible H-bond  
	  h_list.append(atm_tr_h)      # add hydrogen atom of the possible H-bond  
	  a_list.append(atm_tr_acceptor) # add acceptor atom of the possible H-bond  
       atm_tr_acceptor = AtomIteratorMolSet_GetNextAtom(atmitr_acc)  # end of cycle on acceptor atoms
    atm_tr = atmitr.GetNextAtom()  # end of cycle on donor atoms 
  nn = len(d_list)   # number of possible H-bonds
  print "N_HB_possible = ",nn
  dist1_arr = []      # accumulate sum of distances between donor and acceptor 
  dist1_2_arr = []   # accumulate sum of squares of distances between donor and acceptor 
  dist2_arr = []
  dist2_2_arr = []
  ang_arr = []
  ang_2_arr = []
  on_arr = []      # number of MD snapshots this h-bond is ON for normal geometric criteria
  n_md_pt = 0    # number of MD-points 
  for i in range(nn):       # initialize arrays for average distances and angles
    dist1_arr.append(0.0)
    dist1_2_arr.append(0.0)
    dist2_arr.append(0.0)
    dist2_2_arr.append(0.0)
    ang_arr.append(0.0)
    ang_2_arr.append(0.0)
    on_arr.append(0)
 #      n = 0
elif(script_status == SCRIPT_STOP):
  nn = len(d_list)
  for i in range(nn):
    n_on = on_arr[i]
    if( n_on > 0 and n_md_pt > 0):
      dist1_av = dist1_arr[i]/n_md_pt
      dist1_rms =  sqrt(dist1_2_arr[i]/n_md_pt - dist1_av*dist1_av)
      dist2_av = dist2_arr[i]/n_md_pt
      dist2_rms = sqrt(dist2_2_arr[i]/n_md_pt - dist2_av*dist2_av)
      ang_av = (ang_arr[i]/n_md_pt)
      ang_rms = sqrt(ang_2_arr[i]/n_md_pt - ang_av*ang_av)
      ang_av = ang_av * Rad2Deg
      ang_rms = ang_rms * Rad2Deg
      aptr_d = d_list[i] 
      aptr_h = h_list[i]
      aptr_a = a_list[i]
      dref = aptr_d.GetRef()
      href = aptr_h.GetRef()
      aref = MatPoint_GetRef(aptr_a)
      p_on = float(n_on)/float(n_md_pt)  # compute fraction of MD snapshots H-bond is on, float convert integer to real
      idx_d  = atnm_map[dref]
      idx_h  = atnm_map[href]
      idx_a  = atnm_map[aref]
      print >> outf, "%20.20s  %20.20s  %20.20s" %  (dref, href, aref), 
      print >> outf, "%5d %5d %5d" %  (idx_d,idx_h,idx_a), 
      print >> outf, " %5d %6.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f " %  (n_on, p_on,dist1_av, dist2_av,ang_av,dist1_rms, dist2_rms,ang_rms)
  outf.close()
#  outf2.close()
else:
  n_md_pt = n_md_pt + 1
  nn = len(d_list)  # number of possible H-bonds
  print "point ", n_md_pt
  for i in range(nn):   # cycle over possible H-bonds
    aptr_d = d_list[i]   # donor atom of a possible H-bond
    aptr_h = h_list[i]   # hydrogen atom of a possible H-bond
    aptr_a = a_list[i]   # donor atom of a possible H-bond
    dist1 = Vec3D_CalcDistance(aptr_d,aptr_a)*BOHR_TO_ANG
    dist2 = Vec3D_CalcDistance(aptr_h,aptr_a)*BOHR_TO_ANG
    ang   = Vec3D_CalcAngle (aptr_d, aptr_h, aptr_a)
    dist1_arr[i] = dist1_arr[i] + dist1
    dist1_2_arr[i] = dist1_2_arr[i] + dist1*dist1
    dist2_arr[i] = dist2_arr[i] + dist2
    dist2_2_arr[i] = dist2_2_arr[i] + dist2*dist2
    ang_arr[i] = ang_arr[i] + ang
    ang_2_arr[i] = ang_2_arr[i] + ang*ang
    ref_d = aptr_d.GetRef()
    ref_h = aptr_h.GetRef()
    ref_a = MatPoint_GetRef(aptr_a)
#    print ref_d, "-",ref_h, "-",ref_a, "   ", 
#   print "%8.3f %8.3f %8.3f " % (dist1,dist2,ang)
    if( dist1 < 4.6 and dist2 < 3.6 and ang > 1.396 and ang < 3.14):
      on_arr[i] = on_arr[i] + 1   
    #print dref, href,aref  
  print "end_script"
  
