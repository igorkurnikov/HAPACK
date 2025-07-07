# Test setup of PH dependent charges
# load cytb5h_2.hlm 
pmset = GetCurMolSet()
pmset.SetStdResPK()
pmset.SetChargesForPH(5.0)
ritr = ResidueIteratorMolSet(pmset)
rptr = ritr.GetFirstRes()
while (rptr != None):
   rptr.PrintPK()
   rptr = ritr.GetNextRes()
rptr = pmset.GetResByRef("$CYTB5_1$HIS15")
aitr = AtomIteratorAtomSet(rptr)
print "charges HIS:"
aptr = aitr.GetFirstAtom()
while (aptr != None):
    print aptr.GetName(), " ",aptr.GetCharge()
    aptr = aitr.GetNextAtom()
rptr.InterpolResParams("HIS","HIS#PROTA",0.5)
print "charges intermediate between HIS and HIS#PROTA"
aptr = aitr.GetFirstAtom()
while (aptr != None):
    print aptr.GetName(), " ",aptr.GetCharge()
    aptr = aitr.GetNextAtom()




