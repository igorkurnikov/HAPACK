# Calculation of pKa shifts for myoglobin and cytb5 in cytb5-myo complex
# 
pmset = GetCurMolSet()
pmset.SetStdResPK()
pmset.SetChargesForPH(7.0)
electr_mod = pmset.GetElectrostMod(1)
for res_ref in("$MYO_1$HIS24","$MYO_1$HIS36","$MYO_1$HIS48",\
"$MYO_1$HIS64","$MYO_1$HIS81","$MYO_1$HIS82","$MYO_1$HIS97", \
"$MYO_1$HIS113","$MYO_1$HIS116","$MYO_1$HIS119","$MYO_1$HEM154",  \
"$CYTB5_1$HIS15:A","$CYTB5_1$HIS26:A","$CYTB5_1$HIS39:A","$CYTB5_1$HIS63:A","$CYTB5_1$HIS80:A",\
"$CYTB5_1$HEM201:A"):
   rptr = pmset.GetResByRef(res_ref)
   electr_mod.CalcResPK(rptr)
ritr = ResidueIteratorMolSet(pmset)
rptr = ritr.GetFirstRes()
while (rptr != None):
   if(rptr.GetName() == "HIS" or rptr.GetName() == "HEM" ):
      print rptr.PrintPK()
   rptr = ritr.GetNextRes()
