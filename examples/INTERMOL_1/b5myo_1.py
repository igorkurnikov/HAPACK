import molset
from molset import *

def calc_charge(pmset):
# Calculate the total charge of the molecular set
# as a sum of partial atomic charges of the atoms
   aitr = AtomIteratorMolSet(pmset)
   aptr = aitr.GetFirstAtom()

   ch = 0.0
   while(aptr != None):
       ch = ch + aptr.GetCharge()
       aptr = aitr.GetNextAtom()
   print "MolSet charge = ",ch
   return ch

def read_float_col(fname,ncol):
   ifile = open(fname,"r")
   array = []
   while (1):
     str = ifile.readline()
     if(str == ""):
         break
     nums = str.split()
     array.append( float( nums[ncol] ))
   ifile.close()
   return array
      
def set_exp_pk(pmset):
# Set pKa of the protonatable groups in myoglobin and cyt b5
# pKa values are for Sperm Whale but close to horse heart
   rptr = pmset.GetResByRef("$MYO_1$GLY1")
   altst = rptr.GetAltChemState(0)
   if(altst != None):
      altst.pk = 7.47
   rptr = pmset.GetResByRef("$MYO_1$HIS24")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS36")
   altst = rptr.GetAltChemState(0)
   altst.pk = 7.67
   rptr = pmset.GetResByRef("$MYO_1$HIS48")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.42
   rptr = pmset.GetResByRef("$MYO_1$HIS64")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS81")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.65
   rptr = pmset.GetResByRef("$MYO_1$HIS82")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS93")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS97")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.3
   rptr = pmset.GetResByRef("$MYO_1$HIS113")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.51
   rptr = pmset.GetResByRef("$MYO_1$HIS116")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.66
   rptr = pmset.GetResByRef("$MYO_1$HIS119")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.38
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.5
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(1)
   altst.pk = 4.5
   rptr = pmset.GetResByRef("$CYTB5_1$HIS15:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 8.47
   rptr = pmset.GetResByRef("$CYTB5_1$HIS26:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.92
   rptr = pmset.GetResByRef("$CYTB5_1$HIS39:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$CYTB5_1$HIS63:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$CYTB5_1$HIS80:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.3
   rptr = pmset.GetResByRef("$CYTB5_1$HEM201:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.0
   rptr = pmset.GetResByRef("$CYTB5_1$HEM201:A")
   altst = rptr.GetAltChemState(1)
   altst.pk = 5.9
   return 1

def set_exp_pk_shift_prop(pmset):
# Set pKa of the protonatable groups in myoglobin and cyt b5
# pKa values are for Sperm Whale but close to horse heart
# shift pKa of one of cytb5 propionate to 7.0
   rptr = pmset.GetResByRef("$MYO_1$GLY1")
   altst = rptr.GetAltChemState(0)
   if(altst != None):
      altst.pk = 7.47
   rptr = pmset.GetResByRef("$MYO_1$HIS24")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS36")
   altst = rptr.GetAltChemState(0)
   altst.pk = 7.67
   rptr = pmset.GetResByRef("$MYO_1$HIS48")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.42
   rptr = pmset.GetResByRef("$MYO_1$HIS64")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS81")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.65
   rptr = pmset.GetResByRef("$MYO_1$HIS82")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS93")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS97")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.3
   rptr = pmset.GetResByRef("$MYO_1$HIS113")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.51
   rptr = pmset.GetResByRef("$MYO_1$HIS116")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.66
   rptr = pmset.GetResByRef("$MYO_1$HIS119")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.38
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.5
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(1)
   altst.pk = 4.5
   rptr = pmset.GetResByRef("$CYTB5_1$HIS15:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 8.47
   rptr = pmset.GetResByRef("$CYTB5_1$HIS26:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.92
   rptr = pmset.GetResByRef("$CYTB5_1$HIS39:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$CYTB5_1$HIS63:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$CYTB5_1$HIS80:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.3
   rptr = pmset.GetResByRef("$CYTB5_1$HEM201:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.0
   rptr = pmset.GetResByRef("$CYTB5_1$HEM201:A")
   altst = rptr.GetAltChemState(1)
   altst.pk = 5.9
   return 1

def set_myo_pk_nat_1(pmset):
# Set pKa  of the protonatable groups in myoglobin  
# as computed for native myoglobin ion_str = 0.0
# grid = 121  
# 
   rptr = pmset.GetResByRef("$MYO_1$HIS24")
   altst = rptr.GetAltChemState(0)
   altst.pk = 3.51
   rptr = pmset.GetResByRef("$MYO_1$HIS36")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.33
   rptr = pmset.GetResByRef("$MYO_1$HIS48")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.19
   rptr = pmset.GetResByRef("$MYO_1$HIS64")
   altst = rptr.GetAltChemState(0)
   altst.pk = 1.59
   rptr = pmset.GetResByRef("$MYO_1$HIS81")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.33
   rptr = pmset.GetResByRef("$MYO_1$HIS82")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.09
   rptr = pmset.GetResByRef("$MYO_1$HIS93")
   altst = rptr.GetAltChemState(0)
   altst.pk = 2.06
   rptr = pmset.GetResByRef("$MYO_1$HIS97")
   altst = rptr.GetAltChemState(0)
   altst.pk = 7.01
   rptr = pmset.GetResByRef("$MYO_1$HIS113")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.07
   rptr = pmset.GetResByRef("$MYO_1$HIS116")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.27
   rptr = pmset.GetResByRef("$MYO_1$HIS119")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.16
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.55
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(1)
   altst.pk = 3.67
   return 1

def set_myo_pk_dme_1(pmset):
# Set pKa  of the protonatable groups in myoglobin  
# as computed for DME myoglobin ion_str = 0.0
# grid = 121  
# 
   rptr = pmset.GetResByRef("$MYO_1$HIS24")
   altst = rptr.GetAltChemState(0)
   altst.pk = 3.45
   rptr = pmset.GetResByRef("$MYO_1$HIS36")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.29
   rptr = pmset.GetResByRef("$MYO_1$HIS48")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.99
   rptr = pmset.GetResByRef("$MYO_1$HIS64")
   altst = rptr.GetAltChemState(0)
   altst.pk = -0.04
   rptr = pmset.GetResByRef("$MYO_1$HIS81")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.28
   rptr = pmset.GetResByRef("$MYO_1$HIS82")
   altst = rptr.GetAltChemState(0)
   altst.pk = 3.99
   rptr = pmset.GetResByRef("$MYO_1$HIS93")
   altst = rptr.GetAltChemState(0)
   altst.pk = 1.26
   rptr = pmset.GetResByRef("$MYO_1$HIS97")
   altst = rptr.GetAltChemState(0)
   altst.pk = 1.69
   rptr = pmset.GetResByRef("$MYO_1$HIS113")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.05
   rptr = pmset.GetResByRef("$MYO_1$HIS116")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.25
   rptr = pmset.GetResByRef("$MYO_1$HIS119")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.1
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(0)
   altst.pk = 30.0
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(1)
   altst.pk = 30.0
   return 1

def set_myo_pk_s92d_1(pmset):
# Set pKa  of the protonatable groups in myoglobin  
# as computed for myoglobin S92D ion_str = 0.0
# grid = 121  
# 
   rptr = pmset.GetResByRef("$MYO_1$HIS24")
   altst = rptr.GetAltChemState(0)
   altst.pk = 3.53
   rptr = pmset.GetResByRef("$MYO_1$HIS36")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.34
   rptr = pmset.GetResByRef("$MYO_1$HIS48")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.2
   rptr = pmset.GetResByRef("$MYO_1$HIS64")
   altst = rptr.GetAltChemState(0)
   altst.pk = 1.42
   rptr = pmset.GetResByRef("$MYO_1$HIS81")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.36
   rptr = pmset.GetResByRef("$MYO_1$HIS82")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.18
   rptr = pmset.GetResByRef("$MYO_1$HIS93")
   altst = rptr.GetAltChemState(0)
   altst.pk = 2.0
   rptr = pmset.GetResByRef("$MYO_1$HIS97")
   altst = rptr.GetAltChemState(0)
   altst.pk = 10.69
   rptr = pmset.GetResByRef("$MYO_1$HIS113")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.07
   rptr = pmset.GetResByRef("$MYO_1$HIS116")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.27
   rptr = pmset.GetResByRef("$MYO_1$HIS119")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.27
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.02
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(1)
   altst.pk = 3.57
   return 1

def set_myo_pk_s92d_2(pmset):
# Set pKa  of the protonatable groups in myoglobin  
# computed by shift of exp pKa 
# 
   rptr = pmset.GetResByRef("$MYO_1$HIS24")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS36")
   altst = rptr.GetAltChemState(0)
   altst.pk = 7.76
   rptr = pmset.GetResByRef("$MYO_1$HIS48")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.51
   rptr = pmset.GetResByRef("$MYO_1$HIS64")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS81")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.89
   rptr = pmset.GetResByRef("$MYO_1$HIS82")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS93")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS97")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.6
   rptr = pmset.GetResByRef("$MYO_1$HIS113")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.76
   rptr = pmset.GetResByRef("$MYO_1$HIS116")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.72
   rptr = pmset.GetResByRef("$MYO_1$HIS119")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.62
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.6
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(1)
   altst.pk = 4.6 
   return 1

def set_myo_pk_v67r_1(pmset):
# Set pKa  of the protonatable groups in myoglobin  
# as computed for myoglobin V67R ion_str = 0.0
# grid = 121  
# 
   rptr = pmset.GetResByRef("$MYO_1$HIS24")
   altst = rptr.GetAltChemState(0)
   altst.pk = 3.4
   rptr = pmset.GetResByRef("$MYO_1$HIS36")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.3
   rptr = pmset.GetResByRef("$MYO_1$HIS48")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.12
   rptr = pmset.GetResByRef("$MYO_1$HIS64")
   altst = rptr.GetAltChemState(0)
   altst.pk = 0.08
   rptr = pmset.GetResByRef("$MYO_1$HIS81")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.3
   rptr = pmset.GetResByRef("$MYO_1$HIS82")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.04
   rptr = pmset.GetResByRef("$MYO_1$HIS93")
   altst = rptr.GetAltChemState(0)
   altst.pk = 2.51
   rptr = pmset.GetResByRef("$MYO_1$HIS97")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.53
   rptr = pmset.GetResByRef("$MYO_1$HIS113")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.04
   rptr = pmset.GetResByRef("$MYO_1$HIS116")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.25
   rptr = pmset.GetResByRef("$MYO_1$HIS119")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.09
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.06
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(1)
   altst.pk = 3.37
   return 1
   
def set_myo_pk_v67r_dme_1(pmset):
# Set pKa  of the protonatable groups in myoglobin  
# as computed for myoglobin V67R-DME ion_str = 0.0
# grid = 121  
# 
   rptr = pmset.GetResByRef("$MYO_1$HIS24")
   altst = rptr.GetAltChemState(0)
   altst.pk = 3.32
   rptr = pmset.GetResByRef("$MYO_1$HIS36")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.25
   rptr = pmset.GetResByRef("$MYO_1$HIS48")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.89
   rptr = pmset.GetResByRef("$MYO_1$HIS64")
   altst = rptr.GetAltChemState(0)
   altst.pk = -1.9
   rptr = pmset.GetResByRef("$MYO_1$HIS81")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.24
   rptr = pmset.GetResByRef("$MYO_1$HIS82")
   altst = rptr.GetAltChemState(0)
   altst.pk = 3.9
   rptr = pmset.GetResByRef("$MYO_1$HIS93")
   altst = rptr.GetAltChemState(0)
   altst.pk = 0.552
   rptr = pmset.GetResByRef("$MYO_1$HIS97")
   altst = rptr.GetAltChemState(0)
   altst.pk = 1.24
   rptr = pmset.GetResByRef("$MYO_1$HIS113")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.02
   rptr = pmset.GetResByRef("$MYO_1$HIS116")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.23
   rptr = pmset.GetResByRef("$MYO_1$HIS119")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.02
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(0)
   altst.pk = 30.0
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(1)
   altst.pk = 30.0
   return 1

def set_myo_pk_v67r_dme_2(pmset):
# Set pKa  of the protonatable groups in myoglobin  
# computed by shift of exp pKa 
# 
   rptr = pmset.GetResByRef("$MYO_1$HIS24")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS36")
   altst = rptr.GetAltChemState(0)
   altst.pk = 7.67
   rptr = pmset.GetResByRef("$MYO_1$HIS48")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.20
   rptr = pmset.GetResByRef("$MYO_1$HIS64")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS81")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.77
   rptr = pmset.GetResByRef("$MYO_1$HIS82")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS93")
   altst = rptr.GetAltChemState(0)
   altst.pk = -10.0
   rptr = pmset.GetResByRef("$MYO_1$HIS97")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.6
   rptr = pmset.GetResByRef("$MYO_1$HIS113")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.71
   rptr = pmset.GetResByRef("$MYO_1$HIS116")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.68
   rptr = pmset.GetResByRef("$MYO_1$HIS119")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.37
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.6
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(1)
   altst.pk = 4.6 
   return 1

def set_myo_pk_s92d_dme_1(pmset):
# Set pKa  of the protonatable groups in myoglobin  
# as computed for myoglobin V67R-DME ion_str = 0.0
# grid = 121  
# 
   rptr = pmset.GetResByRef("$MYO_1$HIS24")
   altst = rptr.GetAltChemState(0)
   altst.pk = 3.5
   rptr = pmset.GetResByRef("$MYO_1$HIS36")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.33
   rptr = pmset.GetResByRef("$MYO_1$HIS48")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.04
   rptr = pmset.GetResByRef("$MYO_1$HIS64")
   altst = rptr.GetAltChemState(0)
   altst.pk = 0.45
   rptr = pmset.GetResByRef("$MYO_1$HIS81")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.32
   rptr = pmset.GetResByRef("$MYO_1$HIS82")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.13
   rptr = pmset.GetResByRef("$MYO_1$HIS93")
   altst = rptr.GetAltChemState(0)
   altst.pk = 3.37
   rptr = pmset.GetResByRef("$MYO_1$HIS97")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.4
   rptr = pmset.GetResByRef("$MYO_1$HIS113")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.06
   rptr = pmset.GetResByRef("$MYO_1$HIS116")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.27
   rptr = pmset.GetResByRef("$MYO_1$HIS119")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.14
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(0)
   altst.pk = 30.0
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(1)
   altst.pk = 30.0
   return 1

def set_calc_part_pk_1(pmset):
# Set pK of the protonatable groups in myoglobin and cyt b5
# derived from calculations on isolated compounds with grid 65 and ionic strength 0.02
   rptr = pmset.GetResByRef("$MYO_1$HIS24")
   altst = rptr.GetAltChemState(0)
   altst.pk = 2.80
   rptr = pmset.GetResByRef("$MYO_1$HIS36")
   altst = rptr.GetAltChemState(0)
   altst.pk = 7.03
   rptr = pmset.GetResByRef("$MYO_1$HIS48")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.14
   rptr = pmset.GetResByRef("$MYO_1$HIS64")
   altst = rptr.GetAltChemState(0)
   altst.pk = -5.27
   rptr = pmset.GetResByRef("$MYO_1$HIS81")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.58
   rptr = pmset.GetResByRef("$MYO_1$HIS82")
   altst = rptr.GetAltChemState(0)
   altst.pk = 3.36
   rptr = pmset.GetResByRef("$MYO_1$HIS93")
   altst = rptr.GetAltChemState(0)
   altst.pk = -16.29
   rptr = pmset.GetResByRef("$MYO_1$HIS97")
   altst = rptr.GetAltChemState(0)
   altst.pk = 3.36
   rptr = pmset.GetResByRef("$MYO_1$HIS113")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.86
   rptr = pmset.GetResByRef("$MYO_1$HIS116")
   altst = rptr.GetAltChemState(0)
   altst.pk = 2.93
   rptr = pmset.GetResByRef("$MYO_1$HIS119")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.79
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.46
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(1)
   altst.pk = 4.09
   rptr = pmset.GetResByRef("$CYTB5_1$HIS15:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.14
   rptr = pmset.GetResByRef("$CYTB5_1$HIS26:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.84
   rptr = pmset.GetResByRef("$CYTB5_1$HIS39:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = -18.9
   rptr = pmset.GetResByRef("$CYTB5_1$HIS63:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = -13.23
   rptr = pmset.GetResByRef("$CYTB5_1$HIS80:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.88
   rptr = pmset.GetResByRef("$CYTB5_1$HEM201:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 8.11
   rptr = pmset.GetResByRef("$CYTB5_1$HEM201:A")
   altst = rptr.GetAltChemState(1)
   altst.pk = 5.18
   return 1

def set_calc_complex_pk_1(pmset):
# Set pK of the protonatable groups in myoglobin and cyt b5
# derived from calculations on tb5myo_6.hlm complex
#  with grid 65 and ionic strength 0.02
#
   rptr = pmset.GetResByRef("$MYO_1$HIS24")
   altst = rptr.GetAltChemState(0)
   altst.pk = 1.97
   rptr = pmset.GetResByRef("$MYO_1$HIS36")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.58
   rptr = pmset.GetResByRef("$MYO_1$HIS48")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.76
   rptr = pmset.GetResByRef("$MYO_1$HIS64")
   altst = rptr.GetAltChemState(0)
   altst.pk = -5.44
   rptr = pmset.GetResByRef("$MYO_1$HIS81")
   altst = rptr.GetAltChemState(0)
   altst.pk = 7.09
   rptr = pmset.GetResByRef("$MYO_1$HIS82")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.14
   rptr = pmset.GetResByRef("$MYO_1$HIS93")
   altst = rptr.GetAltChemState(0)
   altst.pk = -16.16
   rptr = pmset.GetResByRef("$MYO_1$HIS97")
   altst = rptr.GetAltChemState(0)
   altst.pk = 1.77
   rptr = pmset.GetResByRef("$MYO_1$HIS113")
   altst = rptr.GetAltChemState(0)
   altst.pk = 5.81
   rptr = pmset.GetResByRef("$MYO_1$HIS116")
   altst = rptr.GetAltChemState(0)
   altst.pk = 3.02
   rptr = pmset.GetResByRef("$MYO_1$HIS119")
   altst = rptr.GetAltChemState(0)
   altst.pk = 2.67
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(0)
   altst.pk = 7.00
   rptr = pmset.GetResByRef("$MYO_1$HEM154")
   altst = rptr.GetAltChemState(1)
   altst.pk = 4.92
   rptr = pmset.GetResByRef("$CYTB5_1$HIS15:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.13
   rptr = pmset.GetResByRef("$CYTB5_1$HIS26:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 4.96
   rptr = pmset.GetResByRef("$CYTB5_1$HIS39:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = -18.35
   rptr = pmset.GetResByRef("$CYTB5_1$HIS63:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = -15.64
   rptr = pmset.GetResByRef("$CYTB5_1$HIS80:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.25
   rptr = pmset.GetResByRef("$CYTB5_1$HEM201:A")
   altst = rptr.GetAltChemState(0)
   altst.pk = 6.02
   rptr = pmset.GetResByRef("$CYTB5_1$HEM201:A")
   altst = rptr.GetAltChemState(1)
   altst.pk = 8.76
   return 1

def get_traj_point(traj_file,mol):
#  Read next point of MC trajectory and set molecular coordinates
#
    str = traj_file.readline()
    crds = str.split()
    x = float(crds[1])
    y = float(crds[2])
    z = float(crds[3])

    phi = float(crds[4])
    cost = float(crds[5])
    psi = float(crds[6])

    mol.SetPosEulerTrans(phi,cost,psi,x,y,z)
    return 1


def get_mutants_el_ene(pmset):
#
#
    elmod = pmset.GetElectrostMod(1)
    immod = pmset.GetInterMolMod(1)
    enes = [0.0,0.0,0.0,0.0,0.0,0.0]
    atoms = []
    init_ch = []
    
    for i in range(9):
	 atoms.append(None)
	 init_ch.append(0.0)

    atoms[0] = pmset.GetAtomByRef("$MYO_1$SER92.CB")
    atoms[1] = pmset.GetAtomByRef("$MYO_1$VAL67.CG1")
    atoms[2] = pmset.GetAtomByRef("$MYO_1$VAL67.CG2")
    atoms[3] = pmset.GetAtomByRef("$MYO_1$HEM154.CGD")
    atoms[4] = pmset.GetAtomByRef("$MYO_1$HEM154.O1D")
    atoms[5] = pmset.GetAtomByRef("$MYO_1$HEM154.O2D")
    atoms[6] = pmset.GetAtomByRef("$MYO_1$HEM154.CGA")
    atoms[7] = pmset.GetAtomByRef("$MYO_1$HEM154.O1A")
    atoms[8] = pmset.GetAtomByRef("$MYO_1$HEM154.O2A")

    for i in range(0,9):
	init_ch[i] = atoms[i].GetCharge()

    pel_ene = ptrcreate("double")
    immod.el_pot_field.clear()
    immod.CalcElStaticInter(pel_ene)
    enes[0] = ptrvalue(pel_ene)

# DME:
    for i in range(3,9):
        atoms[i].SetCharge(0.0)
    immod.el_pot_field.clear()
    immod.CalcElStaticInter(pel_ene)
    enes[1] = ptrvalue(pel_ene)
    
# S92D   
    for i in range(3,9):
        atoms[i].SetCharge(init_ch[i])
    atoms[0].SetCharge(-0.7883)

    immod.el_pot_field.clear()
    immod.CalcElStaticInter(pel_ene)
    enes[2] = ptrvalue(pel_ene)

# S92D- DME 
    for i in range(3,9):
        atoms[i].SetCharge(0.0)
    immod.el_pot_field.clear()
    immod.CalcElStaticInter(pel_ene)
    enes[3] = ptrvalue(pel_ene)

# V67R
    for i in range(0,9):
        atoms[i].SetCharge(init_ch[i])

    atoms[1].SetCharge(0.1808)
    atoms[2].SetCharge(0.1808)
    immod.el_pot_field.clear()
    immod.CalcElStaticInter(pel_ene)
    enes[4] = ptrvalue(pel_ene)

# V67R-DME
    for i in range(3,9):
        atoms[i].SetCharge(0.0)

    immod.el_pot_field.clear()
    immod.CalcElStaticInter(pel_ene)
    enes[5] = ptrvalue(pel_ene)
# Restore charges and return:

    for i in range(0,9):
        atoms[i].SetCharge(init_ch[i])
    return enes

def get_mutants_inter_ene_traj(pmset):
#
#
    elmod = pmset.GetElectrostMod(1)
    immod = pmset.GetInterMolMod(1)
    immod.calc_et_rate = 0
    immod.play_back_flag = 1
    atoms = []
    init_ch = []
    enes_all = []
    
    for i in range(9):
	 atoms.append(None)
	 init_ch.append(0.0)

    atoms[0] = pmset.GetAtomByRef("$MYO_1$SER92.CB")
    atoms[1] = pmset.GetAtomByRef("$MYO_1$VAL67.CG1")
    atoms[2] = pmset.GetAtomByRef("$MYO_1$VAL67.CG2")
    atoms[3] = pmset.GetAtomByRef("$MYO_1$HEM154.CGD")
    atoms[4] = pmset.GetAtomByRef("$MYO_1$HEM154.O1D")
    atoms[5] = pmset.GetAtomByRef("$MYO_1$HEM154.O2D")
    atoms[6] = pmset.GetAtomByRef("$MYO_1$HEM154.CGA")
    atoms[7] = pmset.GetAtomByRef("$MYO_1$HEM154.O1A")
    atoms[8] = pmset.GetAtomByRef("$MYO_1$HEM154.O2A")

    for i in range(0,9):
	init_ch[i] = atoms[i].GetCharge()

    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)
    
# DME:
    for i in range(3,9):
        atoms[i].SetCharge(0.0)
    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)
     
# S92D   
    for i in range(3,9):
        atoms[i].SetCharge(init_ch[i])
    atoms[0].SetCharge(-0.7883)

    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)   

# S92D- DME 
    for i in range(3,9):
        atoms[i].SetCharge(0.0)
    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)       

# V67R
    for i in range(0,9):
        atoms[i].SetCharge(init_ch[i])

    atoms[1].SetCharge(0.1808)
    atoms[2].SetCharge(0.1808)
    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)       
   
# V67R-DME
    for i in range(3,9):
        atoms[i].SetCharge(0.0)

    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)       
    
# Restore charges and return:

    for i in range(0,9):
        atoms[i].SetCharge(init_ch[i])

    return enes_all

def get_mutants_inter_ene_traj_2(pmset):
#
#  Compute energies on mutations using pKa of His computed for
#  for myoh_3.hlm geometry grid=121, ion str = 0 
#
    elmod = pmset.GetElectrostMod(1)
    immod = pmset.GetInterMolMod(1)
    immod.calc_et_rate = 0
    immod.play_back_flag = 1
    atoms = []
    init_ch = []
    enes_all = []
        
    set_myo_pk_native_1(pmset)
    pmset.SetChargesForPH(7.0)
    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)
    
# DME:
    set_myo_pk_dme_1(pmset)
    pmset.SetChargesForPH(7.0)
    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)
     
# S92D   
    set_myo_pk_s92d_1(pmset)
    pmset.SetChargesForPH(7.0)
    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)   

# S92D- DME 
    set_myo_pk_s92d_dme_1(pmset)
    pmset.SetChargesForPH(7.0)
    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)       

# V67R
    set_myo_pk_v67r_1(pmset)
    pmset.SetChargesForPH(7.0)
    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)       
   
# V67R-DME
    set_myo_pk_v67r_dme_1(pmset)
    pmset.SetChargesForPH(7.0)
    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)       
    
# Restore charges and return:

    return enes_all



def get_mutants_ie_calc_pk(pmset):
#
#
    elmod = pmset.GetElectrostMod(1)
    immod = pmset.GetInterMolMod(1)
    immod.calc_et_rate = 0
    immod.play_back_flag = 1
    atoms = []
    init_ch = []
    enes_all = []
    
    for i in range(9):
	 atoms.append(None)
	 init_ch.append(0.0)

    atoms[0] = pmset.GetAtomByRef("$MYO_1$SER92.CB")
    atoms[1] = pmset.GetAtomByRef("$MYO_1$VAL67.CG1")
    atoms[2] = pmset.GetAtomByRef("$MYO_1$VAL67.CG2")
    atoms[3] = pmset.GetAtomByRef("$MYO_1$HEM154.CGD")
    atoms[4] = pmset.GetAtomByRef("$MYO_1$HEM154.O1D")
    atoms[5] = pmset.GetAtomByRef("$MYO_1$HEM154.O2D")
    atoms[6] = pmset.GetAtomByRef("$MYO_1$HEM154.CGA")
    atoms[7] = pmset.GetAtomByRef("$MYO_1$HEM154.O1A")
    atoms[8] = pmset.GetAtomByRef("$MYO_1$HEM154.O2A")

    rptr = pmset.GetResByRef("$MYO_1$HEM154")
    altmyo_1 = rptr.GetAltChemState(0)
    altmyo_2 = rptr.GetAltChemState(1)


    for i in range(0,9):
	init_ch[i] = atoms[i].GetCharge()

    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)
    
# DME:
    immod.el_pot_field.clear()
    altmyo_1.inactive_flag = 1
    altmyo_1.pk = 30.0
    altmyo_2.inactive_flag = 1
    altmyo_2.pk = 30.0
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)
    altmyo_1.inactive_flag = 0
    altmyo_2.inactive_flag = 0
    
# S92D   
    for i in range(3,9):
        atoms[i].SetCharge(init_ch[i])
    atoms[0].SetCharge(-0.7883)
    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)   

# S92D- DME 
    immod.el_pot_field.clear()
    altmyo_1.inactive_flag = 1
    altmyo_1.pk = 30.0
    altmyo_2.inactive_flag = 1
    altmyo_2.pk = 30.0
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)       
    altmyo_1.inactive_flag = 0
    altmyo_2.inactive_flag = 0

# V67R
    for i in range(0,9):
        atoms[i].SetCharge(init_ch[i])

    atoms[1].SetCharge(0.1808)
    atoms[2].SetCharge(0.1808)
    immod.el_pot_field.clear()
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)       
   
# V67R-DME
    immod.el_pot_field.clear()
    altmyo_1.inactive_flag = 1
    altmyo_1.pk = 30.0
    altmyo_2.inactive_flag = 1
    altmyo_2.pk = 30.0
    mol2 = pmset.GetMoleculeNum(1)
    mol2.SetPosEulerTrans(-2.275495064,0.671218618,3.007353402,16.141789846,38.292624483, 18.491243967)
    immod.RunMCDock()
    enev = read_float_col("eff_ene_file.dat",1)
    enes_all.append(enev)       
    altmyo_1.inactive_flag = 0
    altmyo_2.inactive_flag = 0

# Restore charges and return:

    for i in range(0,9):
        atoms[i].SetCharge(init_ch[i])
	
    return enes_all

  

