# Calculation of pKa  and EP
pmset = GetCurMolSet()
pmset.set_std_redox_pot = 1
pmset.SetStdResPK_G1()
pmset.SetChargesForPH(7.0)
electr_mod = pmset.GetElectrostMod(1)
electr_mod.rionst = 0.02
electr_mod.epsi = 4.0
electr_mod.nx = 65
electr_mod.ny = 65
electr_mod.nz = 65
pmset.model_titration = 0
pmset.ph = 7.0
pmset.e0 = -0.1
pmset.save_alt_st_inter = 1
pmset.read_alt_st_inter = 0
pmset.save_only_redox_titr = 1
pmset.n_mc_cyc = 1000000
pmset.CalcPKaForSelection()
ritr = ResidueIteratorMolSet(pmset)
rptr = ritr.GetFirstRes()
while (rptr != None):
   nst = rptr.GetNumAltStates()
   if(nst > 0):
     for i in range(nst):
       altst = rptr.GetAltChemState(i)
       print((rptr.GetRef()).c_str(),altst.pk) 
   rptr = ritr.GetNextRes()




