# Calculation of pK shifts in cytb5 / myoglobin complex 
#
def calc_charge(pmset):
   aitr = AtomIteratorMolSet(pmset)
   aptr = aitr.GetFirstAtom()
   ch = 0.0
   while(aptr != None):
       ch = ch + aptr.GetCharge()
       aptr = aitr.GetNextAtom()
   print("MolSet charge = ",ch)
   return 1

pmset = GetCurMolSet()
#calc_charge(pmset)
#pmset.SetStdResPK()
#pmset.SetChargesForPH(7.0)
#calc_charge(pmset)
electr_mod = pmset.GetElectrostMod(1)
rptr = pmset.GetResByRef("$CYTB5_1$HEM201:A")
rptr.SetStdPK()
altst= rptr.GetAltResState(0)
print("std_pk_0 = ",altst.std_pk)
altst= rptr.GetAltResState(1)
print("std_pk_1 = ",altst.std_pk)
altst= rptr.GetAltResState(0)
print("std_pk_2 = ",altst.std_pk)
rptr.PrintPK()
#rptr = pmset.GetResByRef("$MYO_1$HIS24")
#altst= rptr.GetAltResState(0)
#print "His 24 std_pk = ",altst.std_pk

#electr_mod.CalcResPK(rptr,None)
#pmset.SetChargesForPH(7.0)
#for res_ref in("$CYTB5_1$HIS15","$CYTB5_1$HIS26", "$CYTB5_1$HIS39",\
#                     "$CYTB5_1$HIS63","$CYTB5_1$HIS80","$CYTB5_1$HEM201"):
#   rptr = pmset.GetResByRef(res_ref)
#   electr_mod.CalcResPK(rptr)
ritr = ResidueIteratorMolSet(pmset)
rptr = ritr.GetFirstRes()
while (rptr != None):
   if(rptr.GetName() == "HIS" or rptr.GetName() == "HEM"):
     rptr.PrintPK()
   rptr = ritr.GetNextRes()






