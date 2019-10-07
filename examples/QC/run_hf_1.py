mset = GetCurMolSet()
qc_mod = mset.GetQCMod(1)
qc_mod.InitBasis("3-21G")
qc_mod.Run()
print("ene = ",qc_mod.GetEne())
print("HF ene = ",qc_mod.GetHFEne())

