pmset = GetCurMolSet()
pApp.RasMolCmd("select not hoh & not na")
#pmset.CalcMolSurface(HaSurface.SEXCL_SURF)
pmset.CalcMolSurface(HaSurface.SACCESS_SURF)
his_area = dict()
his_area[23] = 0.0
his_area[51] = 0.0
his_area[57] = 0.0
his_area[122] = 0.0
his_area[123] = 0.0
his_area[172] = 0.0
his_at_names = ["CG","ND1","HD1","CE1","HE1","NE2","CD2","HD2"]
for aptr in pmset:
  pres = aptr.GetHostRes()
  nm_at = aptr.GetName() 
  ires = pres.GetSerNo()
  nm_res = pres.GetName()
  if(nm_res == "HIS" and nm_at in his_at_names):
    his_area[ires] += aptr.solv_access_area
for ires in his_area.keys():
  print "his ",ires,"  SA = ",his_area[ires]
print ""
print "end script"
      
