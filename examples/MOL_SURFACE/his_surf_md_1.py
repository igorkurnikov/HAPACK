if(script_status == SCRIPT_START): 
  pmset = GetCurMolSet()
  pApp.RasMolCmd("select not hoh & not na")
  his_at_names = ["CG","ND1","HD1","CE1","HE1","NE2","CD2","HD2"]
  his_area = dict()
  fout = open("his_sa.dat","w")
  i = 0
elif(script_status == SCRIPT_STOP):
  fout.close()
  print "end md analysis"
else:
  #pmset.CalcMolSurface(HaSurface.SEXCL_SURF)
  pmset.CalcMolSurface(HaSurface.SACCESS_SURF)
  his_area[23] = 0.0
  his_area[51] = 0.0
  his_area[57] = 0.0
  his_area[122] = 0.0
  his_area[123] = 0.0
  his_area[172] = 0.0
  for aptr in pmset:
    pres = aptr.GetHostRes()
    nm_at = aptr.GetName() 
    ires = pres.GetSerNo()
    nm_res = pres.GetName()
    if(nm_res == "HIS" and nm_at in his_at_names):
      his_area[ires] += aptr.solv_access_area
  i = i + 1
  print >>fout, i, his_area[23],his_area[51],his_area[57],
  print >>fout, his_area[122],his_area[123],his_area[172]
  print "end pt"
