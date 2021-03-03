import math
import random
pmset = GetCurMolSet()
mol_editor = pmset.GetMolEditor()
mm_mod = pmset.GetMolMechMod()
mm_info  = mm_mod.p_mm_info

at1 = pmset.GetAtomByRef("$PEPT$GLY1.N")
at2 = pmset.GetAtomByRef("$PEPT$GLY1.CA")
at3 = pmset.GetAtomByRef("$PEPT$GLY1.C")
at4 = pmset.GetAtomByRef("$PEPT$GLY2.N")
fout = open("ene.out","w")   # output file for MC trajectory and energies
kt = 1    # kT in kcal/mol
nstep = 1000000    # maximum number of MC steps 
ang = 180.0      # inital value of the torsion angle       (degrees)
dang = 10.0        # maximum step along the angle coordinate (degrees)
npr_freq = 100    # step frequency to print info to console 
ene = 100000.0

for i in range(0,nstep):
  ang_old = ang
  ang = ang + dang*(random.random() -0.5 )
  if( ang > 360.0):
    ang = ang - 360.0
  if( ang < 0.0):
    ang = ang + 360.0
  mol_editor.SetTorsion(at1,at2,at3,at4,ang*DEG_TO_RAD)
  mm_mod.CalcEnergy()
  is_acc = 0
  if( mm_info.pot_ene < ene ):
    is_acc = 1
  else:
    delt_ene = mm_info.pot_ene - ene
    dex = math.exp(- delt_ene/kt)
    rnd = random.random()
    if( dex > rnd):
      is_acc = 1
  if( is_acc == 1):
    ene =  mm_info.pot_ene
    pmset.RefreshAllViews(RFApply)
  else:
    ang = ang_old
  print(i,ang,ene, file=fout)
  if( i % npr_freq == 0):
    print(i, "  ", ang, "  ",mm_info.pot_ene)
fout.close()

