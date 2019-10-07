import sys
if( len(sys.argv) != 6 ):
  print("usage reduce_mdcrd_simple  in.mdcrd  out.mdcrd  nat has_Pbox reduce_factor ") 
fname_inp = sys.argv[1] 
fname_out = sys.argv[2]
nat = int( sys.argv[3] )
has_pbox = int( sys.argv[4] )
reduce_factor = int( sys.argv[5] )
print("fname_inp = ",fname_inp)
print("fname_out = ",fname_out)
print("nat =", nat)
print("has_pbox =", has_pbox)
print("reduce_factor =", reduce_factor)

finp = open(fname_inp,"r")
fout = open(fname_out,"w")
ncrd = nat * 3
nwrite = ncrd/10
if( ncrd % 10 != 0):
  nwrite = nwrite + 1
if( has_pbox ):
  nwrite = nwrite + 1
nskip = nwrite * ( reduce_factor - 1)
print("nwrite = ", nwrite)
print("nskip = ",nskip)
itot = 0
iwrite = 0
iskip = 0
for line in finp:
    itot = itot + 1
    if( itot == 1  ):
        fout.write( line )
        continue
    if( iwrite < nwrite ):
        fout.write( line )
        iwrite = iwrite + 1
        continue
    if( iskip < nskip ):
        iskip = iskip + 1
        continue
    fout.write( line )
    iwrite = 1
    iskip = 0
finp.close()
fout.close()

  
