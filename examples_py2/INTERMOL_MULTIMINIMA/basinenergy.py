#--------------------------------------------------------------------------------------------------------------------------
# This script generates the rmsd vs energy, density around single minima, states during the simulation for each molecule
# and states vs snapshoots files. They can be used for analysis of the multiple minima jumping problem.
# The toy model has 6 symmetric minimum conformations.
#
# The required input files are:
# all_mc_pts.dat          accepted conformations of the MC trajectory
# empir_param.dat         contains the van der Waals radius
# bundle50A.hlm           the toy model topology
# min(i).pdb              for i= 1,2,...6 contains the minimum energy conformations of the toy model

# The output files are:
# rmsd_state(i).out for i=1,2,...6        the rmsd vs energy file for each i posible conformation
# density_state(i).dat for i=1,2,...6     the density of points around the i minima (from the rmsd vs energy files)
# eff_ene_file.dat                        the energy profile of the MC trajectory
# minimum(i).hlm for i=1,2,...6           the minimum energy configurations reached in each i conformation explored
#--------------------------------------------------------------------------------------------------------------------------

import os
nm1 = []

# name of the trajectory file, format (snaphoot,x,y,z,q1,q2,q3,q4)
trj_file = "all_mc_pts.dat"
#name of file with rmsd values
out_file1=open("rmsd_state1.out","w")
out_file2=open("rmsd_state2.out","w")
out_file3=open("rmsd_state3.out","w")
out_file4=open("rmsd_state4.out","w")
out_file5=open("rmsd_state5.out","w")
out_file6=open("rmsd_state6.out","w")
out_zone1=open("density_state1.dat","w")
out_zone2=open("density_state2.dat","w")
out_zone3=open("density_state3.dat","w")
out_zone4=open("density_state4.dat","w")
out_zone5=open("density_state5.dat","w")
out_zone6=open("density_state6.dat","w")
out_file7=open("accepted_states.dat","w")
out_file8=open("snap_states.dat","w")

ene_aux = open("energy_snap.dat","w")

# function that counts the numbers of conformations around the minima
# In this case counts points in the rectangle areas defined by the local
# minima value and -2 Kcal/mol by 0 - 1.5, 1.5- 3.5, 3.5 -5.5,..., 9.5 - 11.5 Angstroms

def tabledat(ene1,eps1,zone,min_ene1, min_rmsd1):
  #print ene1, eps1
  if ene1 < -2.0 and eps1 < 1.5:
      zone[0] +=1
      #print zone[0]
  if ene1 < -2.0 and 1.5 <= eps1 <3.5:
      zone[1] +=1
      #print zone[1]
  if ene1 < -2.0 and 3.5 <= eps1 < 5.5:
      zone[2] +=1
      #print zone[2]
  if ene1 < -2.0 and 5.5 <= eps1 < 7.5:
      zone[3] +=1
      #print zone[3]
  if ene1 < -2.0 and 7.5 <= eps1 < 9.5:
      zone[4] +=1
      #print zone[4]
  if ene1 < -2.0 and 9.5 <= eps1 < 11.5:
      zone[5] +=1
      #print zone[5]
  if ene1 < min_ene1:
      min_ene1 = ene1
      min_rmsd1 = eps1
  return zone, min_ene1, min_rmsd1

# what HLM file to use for fitting?
nm1.append("bundle50A.hlm") # HLM file 

# load HLM file on pmset1 object
ifile = open(trj_file, "r")
pmset1 = HaMolSet()
pmsetx = HaMolSet()
pmset1.FetchFile(FormatHarlem,nm1[0])
pmsetx.FetchFile(FormatHarlem,nm1[0])
# change the identification number of residues
n=0
aitr_m1 = AtomIteratorMolSet(pmsetx)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
    at_name = aptr.GetName()
    if(at_name =="N"):	
	    n= n+1
    pres = aptr.GetHostRes()
    pres.serno = n    
    aptr = aitr_m1.GetNextAtom()

#save the structure in PDB format
pmsetx.SavePDBFile("temp.pdb")

#creates two HaMolset objects from the PDB file
pmset = HaMolSet()
pmset2 =HaMolSet()
pmset.FetchFile(FormatPDB,"temp.pdb")   # initial snapshoot for the RMS calculation
pmset2.FetchFile(FormatPDB,"temp.pdb")  # (current snapshoot)for the RMS calculation

pmseta1 = HaMolSet()
pmseta2 = HaMolSet()
pmseta3 = HaMolSet()
pmseta4 = HaMolSet()
pmseta5 = HaMolSet()
pmseta6 = HaMolSet()
pmseta1.FetchFile(FormatPDB,"min1.pdb")
pmseta2.FetchFile(FormatPDB,"min2.pdb")
pmseta3.FetchFile(FormatPDB,"min3.pdb")
pmseta4.FetchFile(FormatPDB,"min4.pdb")
pmseta5.FetchFile(FormatPDB,"min5.pdb")
pmseta6.FetchFile(FormatPDB,"min6.pdb")

# change the identification number of residues
n=0
aitr_m1 = AtomIteratorMolSet(pmseta1)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
    at_name = aptr.GetName()
    if(at_name =="N"):	
	    n= n+1
    pres = aptr.GetHostRes()
    pres.serno = n    
    aptr = aitr_m1.GetNextAtom()

n=0
aitr_m1 = AtomIteratorMolSet(pmseta2)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
    at_name = aptr.GetName()
    if(at_name =="N"):	
	    n= n+1
    pres = aptr.GetHostRes()
    pres.serno = n    
    aptr = aitr_m1.GetNextAtom()


n=0
aitr_m1 = AtomIteratorMolSet(pmseta3)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
    at_name = aptr.GetName()
    if(at_name =="N"):	
	    n= n+1
    pres = aptr.GetHostRes()
    pres.serno = n    
    aptr = aitr_m1.GetNextAtom()


n=0
aitr_m1 = AtomIteratorMolSet(pmseta4)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
    at_name = aptr.GetName()
    if(at_name =="N"):	
	    n= n+1
    pres = aptr.GetHostRes()
    pres.serno = n    
    aptr = aitr_m1.GetNextAtom()


n=0
aitr_m1 = AtomIteratorMolSet(pmseta5)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
    at_name = aptr.GetName()
    if(at_name =="N"):	
	    n= n+1
    pres = aptr.GetHostRes()
    pres.serno = n    
    aptr = aitr_m1.GetNextAtom()


n=0
aitr_m1 = AtomIteratorMolSet(pmseta6)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
    at_name = aptr.GetName()
    if(at_name =="N"):	
	    n= n+1
    pres = aptr.GetHostRes()
    pres.serno = n    
    aptr = aitr_m1.GetNextAtom()

#gets the number of molecules in pmset1
nmol = pmset1.GetNMol()
print "nmol" , nmol

#initializes the Empirical Potential method 
empmod = pmset1.GetEmpiricalMod(1)
empmod.Initialize()

#creates auxiliary arrays
atl1 = AtomGroup()
atl2 = AtomGroup()
atl3 = AtomGroup()
atl4 = AtomGroup()
atl5 = AtomGroup()
atl6 = AtomGroup()
prue = AtomGroup()
atm_ca_arr_set1 = []
atm_ca_arr_set2 = []
new_at = HaAtom()

#creates the first atom group with the template structure atl1
#the template atom group in this case will contains CA atoms
#the identity of atoms can be modified ie: ("CA" and "N")

aitr_m1 = AtomIteratorMolSet(pmseta1)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
    at_name = aptr.GetName()
    if(at_name =="CA") :
        atl1.InsertAtom(aptr)
    aptr = aitr_m1.GetNextAtom()

aitr_m1 = AtomIteratorMolSet(pmseta2)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
    at_name = aptr.GetName()
    if(at_name =="CA") :
        atl2.InsertAtom(aptr)
    aptr = aitr_m1.GetNextAtom()

aitr_m1 = AtomIteratorMolSet(pmseta3)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
    at_name = aptr.GetName()
    if(at_name =="CA") :
        atl3.InsertAtom(aptr)
    aptr = aitr_m1.GetNextAtom()

aitr_m1 = AtomIteratorMolSet(pmseta4)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
    at_name = aptr.GetName()
    if(at_name =="CA") :
        atl4.InsertAtom(aptr)
    aptr = aitr_m1.GetNextAtom()

aitr_m1 = AtomIteratorMolSet(pmseta5)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
    at_name = aptr.GetName()
    if(at_name =="CA") :
        atl5.InsertAtom(aptr)
    aptr = aitr_m1.GetNextAtom()

aitr_m1 = AtomIteratorMolSet(pmseta6)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
    at_name = aptr.GetName()
    if(at_name =="CA") :
        atl6.InsertAtom(aptr)
    aptr = aitr_m1.GetNextAtom()

#creates and array of CA atoms from pmset2 object(atm_ca_arr_set2)
aitr_m1 = AtomIteratorMolSet(pmset2)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
	at_name = aptr.GetName()
 	if(at_name =="CA") :
	      atm_ca_arr_set2.append( aptr)
	aptr = aitr_m1.GetNextAtom()

#creates an array of CA atoms from pmset1 object(atm_ca_arr_set1)
aitr_m1 = AtomIteratorMolSet(pmset1)
aptr = aitr_m1.GetFirstAtom()
while( aptr != None):
	at_name = aptr.GetName()
	if(at_name =="CA") :
		atm_ca_arr_set1.append( aptr)
	aptr = aitr_m1.GetNextAtom()

#creates the translation and quaternion vectors
#creates qang,qx,qy,qz components
trans = Vec3D()
quat = Quaternion()
qang = doublep()
qx = doublep()
qy = doublep()
qz = doublep()

#start to read the trajectory file
snap = 0 #count the number of snapshoots
str = ifile.readline() #read each line of the trajectory file 
state_vec=[0,0,0,0]
state_vec1=[0,1,2,3]
state_vec2=[0,1,3,2]
state_vec3=[0,2,1,3]
state_vec4=[0,2,3,1]
state_vec5=[0,3,2,1]
state_vec6=[0,3,1,2]
zone1 = [0,0,0,0,0,0]
zone2 = [0,0,0,0,0,0]
zone3 = [0,0,0,0,0,0]
zone4 = [0,0,0,0,0,0]
zone5 = [0,0,0,0,0,0]
zone6 = [0,0,0,0,0,0]
min_ene = [10000000,10000000,10000000,10000000,10000000,10000000]
min_rmsd = [0,0,0,0,0,0]
rate=0
rate_aux =0
while (1):
         snap = snap +1
         coord_arr=[]
         imol=0
         # Set positions of rigid bodies in pmset1 objects
         while ( imol < nmol):   
            nums = str.split() #read the coordinate values for each rigid body
            qang.assign(float(nums[4]))
            qx.assign(float(nums[5]))
            qy.assign(float(nums[6]))
            qz.assign(float(nums[7]))
            trans.SetX( float(nums[1]) )
            trans.SetY( float(nums[2]))
            trans.SetZ( float(nums[3]))
            state_vec[imol]=(int(nums[8]))
            #get the molecule imol from pmset1 object
            pmol1 = pmset1.GetMoleculeNum(imol)
            quat.SetQuaternion(qang,qx,qy,qz)
            #translate and rotate the molecule imol
            pmol1.SetQuaternionTrans(quat, trans)
            str = ifile.readline()
            imol = imol+1 #while loop
         print >> out_file7, snap, state_vec[0], state_vec[1], state_vec[2], state_vec[3] 
         if(snap > 0 and snap%1 == 0):

          #creates atl2 Atomgroup that will contains the CA atoms of the current snapshoot
          atl_curr = AtomGroup()

          #gets the x,y,z coordinates for each CA atom in pmset1 (already translated and rotated) 
          for i in range (len(atm_ca_arr_set1) ):
                new_at = atm_ca_arr_set1[i]
                coord_arr.append(new_at.GetX() )
                coord_arr.append(new_at.GetY() )
                coord_arr.append(new_at.GetZ() )

          #update the x,y,z coordinates for each CA atom in pmset2 (pmset2 now translated and rotated)
          #add all CA atoms from pmset2 to atl2 group
          ind =0
          for i in range (len(atm_ca_arr_set2) ):
                new_at = atm_ca_arr_set2[i]
                new_at.SetX(coord_arr[ind])
                ind = ind+1
                new_at.SetY(coord_arr[ind])
                ind = ind+1
                new_at.SetZ(coord_arr[ind])
                ind = ind+1
                atl_curr.InsertAtom(new_at)

          #calculates the RMSD and energy value
          enep = empmod.ScoreEnergy()
          
          if (state_vec1==state_vec):  #123
              eps = pmseta1.OverlapMol(atl1, atl_curr)
              #print "state_vec1\n"
              print >> out_file8, snap, 1
              print >> out_file1, snap, eps, enep
              temp = min_ene[0]
              zone1,min_ene[0],min_rmsd[0] = tabledat(enep,eps,zone1,min_ene[0],min_rmsd[0])
              if (temp !=min_ene[0]):
                  pmset1.SaveHarlemFile("minimum1.hlm")
              rate_aux =rate_aux+1
              
          if (state_vec2==state_vec):  #132
              eps = pmseta2.OverlapMol(atl2, atl_curr)
              #print "state_vec2\n"
              print >> out_file8, snap, 2
              print >> out_file2, snap, eps, enep
              temp = min_ene[1]
              zone2,min_ene[1],min_rmsd[1] = tabledat(enep,eps,zone2,min_ene[1],min_rmsd[1])
              if (temp !=min_ene[1]):
                  pmset1.SaveHarlemFile("minimum2.hlm")
              rate_aux =rate_aux+1
              
          if (state_vec3==state_vec):  #213
              eps = pmseta3.OverlapMol(atl3, atl_curr)
              #print "state_vec3\n"
              print >> out_file8, snap, 3
              print >> out_file3, snap, eps, enep
              temp = min_ene[2]
              zone3,min_ene[2],min_rmsd[2] = tabledat(enep,eps,zone3,min_ene[2],min_rmsd[2])
              if (temp !=min_ene[2]):
                  pmset1.SaveHarlemFile("minimum3.hlm")
              rate_aux =rate_aux+1
              
          if (state_vec4==state_vec):  #231
              eps = pmseta4.OverlapMol(atl4, atl_curr)
              #print "state_vec4\n"
              print >> out_file8, snap, 4
              print >> out_file4, snap, eps, enep
              temp = min_ene[3]
              zone4,min_ene[3],min_rmsd[3] = tabledat(enep,eps,zone4,min_ene[3],min_rmsd[3])
              if (temp !=min_ene[3]):
                  pmset1.SaveHarlemFile("minimum4.hlm")
              rate_aux =rate_aux+1
              
          if (state_vec5==state_vec):  #321
              eps = pmseta5.OverlapMol(atl5, atl_curr)
              #print "state_vec5\n",
              print >> out_file8, snap, 5
              print >> out_file5, snap, eps, enep
              temp = min_ene[4]
              zone5,min_ene[4],min_rmsd[4] = tabledat(enep,eps,zone5,min_ene[4],min_rmsd[4])
              if (temp !=min_ene[4]):
                  pmset1.SaveHarlemFile("minimum5.hlm")
              rate_aux =rate_aux+1
              
          if (state_vec6==state_vec):  #312
              eps = pmseta6.OverlapMol(atl6, atl_curr)
              #print "state_vec6\n"
              print >> out_file8, snap, 6
              print >> out_file6, snap, eps, enep
              temp = min_ene[5]
              zone6,min_ene[5],min_rmsd[5] = tabledat(enep,eps,zone6,min_ene[5],min_rmsd[5])
              if (temp !=min_ene[5]):
                  pmset1.SaveHarlemFile("minimum6.hlm")
              rate_aux =rate_aux+1
          print >> ene_aux, snap, enep    
          rate=rate+1
         if str == '':
             print "rate = ", rate
             print "rate_aux =", rate_aux
             print >> out_zone1, zone1[0], zone1[1], zone1[2], zone1[3], zone1[4], zone1[5]
             print >> out_zone2, zone2[0], zone2[1], zone2[2], zone2[3], zone2[4], zone2[5]
             print >> out_zone3, zone3[0], zone3[1], zone3[2], zone3[3], zone3[4], zone3[5]
             print >> out_zone4, zone4[0], zone4[1], zone4[2], zone4[3], zone4[4], zone4[5]
             print >> out_zone5, zone5[0], zone5[1], zone5[2], zone5[3], zone5[4], zone5[5]
             print >> out_zone6, zone6[0], zone6[1], zone6[2], zone6[3], zone6[4], zone6[5]
             print "snap = " , snap
             break

out_file1.close()
out_file2.close()
out_file3.close()
out_file4.close()
out_file5.close()
out_file6.close()
out_file7.close()
out_file8.close()
out_zone1.close()
out_zone2.close()
out_zone3.close()
out_zone4.close()
out_zone5.close()
out_zone6.close()
ifile.close()
ene_aux.close()
print "Work completed"


