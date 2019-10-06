finp = open("dudl_vals.dat","r")
fout = open("dudl_red.dat","w")
n = 0
sum = 0.0
nred = 100000
for line in finp:
  vals = line.split()
  if( len(vals) != 2):
    continue
  idx = int(vals[0])
  n = n+1
  sum = sum + float(vals[1])
  if( n % nred == 0):
    sum = sum/nred
    print("%d %16.9f " % (n,sum), file=fout)
    sum = 0.0
finp.close()
fout.close()

  