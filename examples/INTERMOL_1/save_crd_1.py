i= i+1
if(i % 10 == 0):
   fname = "myoh_3_box_" + str(i) + ".hlm"
   print(" Save to ", fname)
   pmset = GetCurMolSet() 
   pmset.SaveHarlemFile(fname)
