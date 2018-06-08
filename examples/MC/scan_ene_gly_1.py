pmset = GetCurMolSet()
mol_editor = pmset.GetMolEditor(1)
mm_mod = pmset.GetMolMechMod(1)
mm_info  = mm_mod.p_mm_info

at1 = pmset.GetAtomByRef("$PEPT$GLY1.N")
at2 = pmset.GetAtomByRef("$PEPT$GLY1.CA")
at3 = pmset.GetAtomByRef("$PEPT$GLY1.C")
at4 = pmset.GetAtomByRef("$PEPT$GLY2.N")
fout = open("ene.out","w")
for i in range(0,36):
  ang = 10.0*i
  mol_editor.SetTorsion(at1,at2,at3,at4,ang*DEG_TO_RAD)
  mm_mod.CalcEnergy()
  print ang, mm_info.pot_ene
  print >> fout, ang, mm_info.pot_ene
  pmset.RefreshAllViews(RFApply)
fout.close()

