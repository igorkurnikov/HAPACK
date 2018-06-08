import os
for i in range(1,126):
  prefix1 = "E66C_Y67F_frag2_"
  prefix2 = "E66C_Y67F_frag2_6311G_N0083_"
  fn1 = prefix1 + str(i) + "0.chk"
  fn2 = prefix2 + str(i) + "0.chk"
  cmd = "mv " + fn1 + " " + fn2
  print cmd
  os.system(cmd)