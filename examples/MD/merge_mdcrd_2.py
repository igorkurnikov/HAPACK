import sys
nf = len( sys.argv ) - 1
fout = open("traj_full.mdcrd","w")
for i in range(nf):
  fname = sys.argv[i+1]
  print fname
  finp = open(fname,"r")
  il = 0
  for line in finp:
    il = il + 1
    if( il == 1 and i != 1 ):  
      continue
    fout.write( line )
  finp.close()
fout.close()

  
