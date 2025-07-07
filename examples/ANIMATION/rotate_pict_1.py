pmset = GetCurMolSet()
pview = pmset.GetActiveMolView()
fname = ""
for i in range(1,19):
  print i
  fname = "test_" + str(1000+i) + ".pic"
  pApp.RasMolCmd("rotate Y 20")
  pview.WritePICTFile(fname)
print "end"
