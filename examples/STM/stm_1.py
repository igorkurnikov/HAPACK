file_res = open("res.out","w") 
for i in range(100,105):
   pmset = MolSet()
   fname = "pu-10.0-" + str(i) + "-l.mol" 
   pmset.FetchFile(FormatMDL,0,fname)
   print(i," Na = ",pmset.GetNAtoms())
   qcmod = pmset.GetQCMod(1)
   qcmod.InitBasis("minb6g")
   print("NB = ",qcmod.GetNBfunc())
   stmmod = pmset.GetSTMMod(1)
   stmmod.num_mol_ao = 56
   stmmod.CalcTMatr1()
   print(fname,stmmod.tmat_el, file=file_res) 
   pmset.DeleteAll()
file_res.close()
