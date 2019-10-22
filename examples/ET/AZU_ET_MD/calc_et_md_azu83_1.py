if(script_status == SCRIPT_START):
  pmset = GetCurMolSet()
  mm_mod = pmset.GetMolMechMod(1)
  pmset_f = MolSet()
  pmset_f.LoadHarlemFile("azu83_120P30.hlm")
  aitr_f = AtomIteratorMolSet(pmset_f)
  outf = open("hda_1.dat","w")
  ns = 0
  print(" time   hda", file=outf)
  qcmod = pmset_f.GetQCMod(1)
  qcmod.ndo_method = ZINDO_1
  qcmod.SetCharge(2)
  qcmod.max_it_avg = 20
  etmod = pmset_f.GetETCouplMod(1)
elif(script_status == SCRIPT_STOP):
#  After MD trajectory analysis - close output file
  outf.close()
else:
  aptr_f = aitr_f.GetFirstAtom() 
  while( aptr_f != None):   # set coordinates for not terminating atoms
    res = aptr_f.GetHostRes()
    res_name = res.GetName()
    if(res_name != "TRM"):
      atref = aptr_f.GetRef()
      aptr = pmset.GetAtomByRef(atref)
      if(aptr != None):
        aptr_f.SetX( aptr.GetX())
	aptr_f.SetY( aptr.GetY())
	aptr_f.SetZ( aptr.GetZ())
    aptr_f = aitr_f.GetNextAtom()
  aptr_f = aitr_f.GetFirstAtom()
  while( aptr_f != None):  # set coordinates for  terminating atoms
    res = aptr_f.GetHostRes()
    res_name = res.GetName()
    if(res_name == "TRM"):
      bnd_at_f1 = AtomGroup()
      aptr_f.GetBondedAtoms(bnd_at_f1)
      aitr_bf1 = AtomIteratorAtomGroup( bnd_at_f1 )
      aptr_f2 = aitr_bf1.GetFirstAtom()
      bnd_at_f2 = AtomGroup()
      aptr_f2.GetBondedAtoms(bnd_at_f2)
      aitr_bf2 = AtomIteratorAtomGroup( bnd_at_f1 )
      aptr_b = aitr_bf2.GetFirstAtom()
      vec = Vec3D()
      vec.SetX(0.0)
      vec.SetY(0.0)
      vec.SetZ(0.0)
      while( aptr_b != None):
        if( aptr_b.GetX() != aptr_f.GetX()):
	  vd = Vec3D()
	  vd.SetX( aptr_b.GetX() - aptr_f2.GetX())
	  vd.SetY( aptr_b.GetY() - aptr_f2.GetY())
	  vd.SetZ( aptr_b.GetZ() - aptr_f2.GetZ())
	  vd.normalize()
	  vec.SetX( vec.GetX() + vd.GetX())
	  vec.SetY( vec.GetY()  + vd.GetY())
	  vec.SetZ( vec.GetZ() + vd.GetZ())
        aptr_b = aitr_bf2.GetNextAtom()
      vec.normalize()
      vec.ScaleCoord(-1.08/BOHR_TO_ANG)
      aptr_f.SetX( aptr_f2.GetX() + vec.GetX())
      aptr_f.SetY( aptr_f2.GetY()  + vec.GetY())
      aptr_f.SetZ( aptr_f2.GetZ() + vec.GetZ())
    aptr_f = aitr_f.GetNextAtom()
  ns = ns + 1
  str_ns = str(ns)
  if( ns < 10): 
    str_ns = "0" + str_ns
  if( ns < 100): 
    str_ns = "0" + str_ns
  fname = "snap" + str_ns + ".hlm"
  pmset_f.SaveHarlemFile(fname)
  qcmod.RunCNDOThread()
  etmod.CopyEigVecsFromMO()
  idx_mo = HaVec_int()
  nmo = 8
  idx_mo.newsize(nmo)
  for i in range(nmo):
    idx_mo.SetVal(i+1,ihomo - nmo + i + 1 )
  etmod.AddRedoxOrbFromEigVec(idx_mo)
  etmod.SetTunEne(-0.2)
  etmod.CalcHDAfromGF()
  nr = etmod.da_coupl_val.num_rows()
  nc = etmod.da_coupl_val.num_cols()
  print(" HDA from GF: ")
  for ir in range(nr):
    print(ns, end=' ', file=outf)
    for ic in range(nc):
         print("%12.6e" % etmod.da_coupl_val.GetVal(ir+1,ic+1), end=' ', file=outf)
         print(etmod.da_coupl_val.GetVal(ir+1,ic+1), end=' ')
  print(" ", file=outf)

