finp = open("dudl_vals_lmb0769_t0_100.dat","r")
fout = open("dudl_red.dat","w")
n = 0
sum = 0.0
nred = 10000
for line in finp:
  vals = line.split()
  if( len(vals) != 2):
    continue
  idx = int(vals[0])
  n = n+1
  sum = sum + float(vals[1])
  if( n % nred == 0):
    sum = sum/nred
    print >>fout, "%d %16.9f " % (n,sum)
    sum = 0.0
finp.close()
fout.close()

  