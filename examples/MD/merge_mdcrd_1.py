fout = open("traj_full.mdcrd","w")
for i in range(1,13):
  sn = str(i)
  if( i < 10):
    sn = "0" + sn
  fname = "mg_wat_1_ch_1.6_" + sn + ".mdcrd"
  print(fname)
  finp = open(fname,"r")
  il = 0
  for line in finp:
    il = il + 1
    if( il == 1 and i != 1 ):  
      continue
    fout.write( line )
  finp.close()
fout.close()

  
