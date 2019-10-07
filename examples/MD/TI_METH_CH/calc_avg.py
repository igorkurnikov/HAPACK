finp = open("dudl_vals.dat","r")
n = 0
sum = 0.0
for line in finp:
  vals = line.split()
  if( len(vals) != 2):
    continue
  idx = int(vals[0])
  if(idx < 2000):
    continue
  if(idx > 500000):
    continue
  n = n+1
  sum = sum + float(vals[1])
finp.close()
print("n =",n)
print("average =",sum/n)

  