#script_status = SCRIPT_START
if(script_status == SCRIPT_START):
  pmset_md = GetCurMolSet()
  pmset_pka = HaMolSet()
  pmset_pka.LoadHarlemFile("dtox_tdom_hid_prot_his.hlm")
  atmap = AtomMapping(pmset_md,pmset_pka)
  atmap.Map2to1ByAtomRef()
  atmap.PrintInfo(1)
  atmap.SyncAtomCrd2From1()
  elmod = pmset_pka.GetElectrostMod(1)
  prot_mod = pmset_pka.GetProtonRedoxMod(1)
  elmod.nx = 121
  elmod.ny = 121
  elmod.nz = 121
  elmod.epsi = 6.0
  elmod.epsout = 80.0
  elmod.rionst = 0.15
  prot_mod.SetStdPKa_G1()
  fout = open("pka_md.out","w")
  prot_mod.multi_site_pop_method.SetWithValue(MC_MULTI_SITE_CALC)
  prot_mod.n_mc_cyc = 100000
  prot_mod.SetPH(6.5)
#  prot_mod.CalcPKaForSelection()
  nalt = prot_mod.alt_chem_states.size()
  print nalt
  ipt = 0
#  for i in range(nalt):
#    print i,prot_mod.alt_chem_states[i].pk
elif(script_status == SCRIPT_STOP):
  fout.close()
  print "end of md script"
else:
  atmap.SyncAtomCrd2From1()
  prot_mod.CalcPKaForSelection()
  print >>fout, "%5d " % ipt,
  for ialt in range(nalt):
    print >>fout, "%7.3f " % prot_mod.alt_chem_states[ialt].pk,
  print >> fout, "  "
  ipt = ipt+1

  


 