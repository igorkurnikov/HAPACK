mset = GetCurMolSet()
i = 0
fout = open("atom_id.dat","w")
for aptr in (mset):
  i = i+1
  print >>fout, i, "  ",aptr.GetRef()
fout.close()
print "end of script"

