import time
pmset = GetCurMolSet()
aitr = AtomIteratorMolSet(pmset)
aptr = aitr.GetFirstAtom()
pview = pmset.GetActiveMolView()
fname = ""
for i in range(1,19):
  print i, aptr.x,aptr.y,aptr.z
  fname = "test_" + str(1000+i) + ".jpeg"
  pApp.RasMolCmd("rotate Y 20")
  pview.ApplyTransform()
  pview.DrawFrame()
  pview.WriteJPEGFile(fname)
print "end"
