pmset = GetCurMolSet()
pview = pmset.GetActiveMolView()
fname = ""
for i in range(1,19):
  print i
  fname = "test_" + str(1000+i) + ".gif"
  pApp.RasMolCmd("rotate Y 20")
  pview.WriteGIFFile(fname)
print "end"
