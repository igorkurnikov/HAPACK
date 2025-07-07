pmset = GetCurMolSet()
dna_mod = pmset.GetNuclAcidMod(1)
dna_mod.SetSeq("GC")
dna_mod.SetFFtype("AMBER94")
#dna_mod.nsym_unit.SetVal(1,10)
#dna_mod.max_iter = 2000
#dna_mod.ene_per_unit_flag = 1
#dna_mod.sup_helix_rad = 60.0
dna_mod.BuildNuclAcid()
dna_mod.MinEne()