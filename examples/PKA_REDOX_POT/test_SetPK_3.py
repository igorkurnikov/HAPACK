# Calculation of pK shifts for all resiudes in the molecular set
# 
pmset = GetCurMolSet()
#pmset.SetStdResPK()
pmset.SetChargesForPH(7.0)
ritr = ResidueIteratorMolSet(pmset)
electr_mod = pmset.GetElectrostMod(1)
electr_mod.rionst = 0.02
rptr = ritr.GetFirstRes()
while (rptr != None):
   nalt = rptr.GetNumAltStates()
   print(rptr.GetRef().c_str())
   electr_mod.CalcResPK(rptr,None)
   rptr = ritr.GetNextRes()
rptr = ritr.GetFirstRes()
while (rptr != None):
   if(rptr.GetName() == "HIS" or rptr.GetName() == "HEM"):
      print(rptr.PrintPK())
   rptr = ritr.GetNextRes()





