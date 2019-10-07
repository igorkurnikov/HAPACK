if(script_status == SCRIPT_START):
  pmset = GetCurMolSet()
  mm_mod = pmset.GetMolMechMod(1)
  pmset_f = HaMolSet()
  pmset_f.LoadHarlemFile("azu83_120P30.hlm")
  molm_f   = pmset_f.GetMolByName("AZU_83_1_ORIG")
  molext_f = pmset_f.GetMolByName("EXTERNAL_CHARGES")
# 
  atmap1 = AtomMapping(pmset,molm_f)
  atmap1.Map2to1ByAtomRef()
  atmap1.PrintInfo(1)
#  
  atmap2= AtomMapping(molm_f,molext_f)
  atmap2.Map2to1ByAtomDistance()
  atmap2.PrintInfo(1)
  
  atmap1.SyncAtomCrd2From1()
  atmap2.SyncAtomCrd2From1()
#
  ns = 0
elif(script_status == SCRIPT_STOP):
#
  print("Stop MD Analysis")
else:
  ns = ns + 1
  str_ns = str(ns)
  if( ns < 10): 
    str_ns = "0" + str_ns
  if( ns < 100): 
    str_ns = "0" + str_ns
  fname = "snap" + str_ns + ".hlm"
  if( ns % 100 == 0):
    atmap1.SyncAtomCrd2From1()
    atmap2.SyncAtomCrd2From1()
    pmset_f.SaveHarlemFile(fname)

 