mset = GetCurMolSet()
i = 0
fout = open("atom_id.dat","w")
for aptr in (mset):
  i = i+1
  print(i, "  ",aptr.GetRef(), file=fout)
fout.close()
print("end of script")

