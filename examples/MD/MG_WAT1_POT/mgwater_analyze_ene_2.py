from subprocess import *
def calc_ene():
  ene = 0.0
  p = Popen(["/cygdrive/c/tinker/bin/analyze","mgwater1","E"],stdout=PIPE)
  for line in p.stdout: 
    if( line.find("Total Potential Energy") != -1 ):
      ene = float( line[35:50] )
  return ene

def set_dist( dist ):
  #dist = 1.88013
  #dist = 2.0
  xmg = 1.933960 - dist 
  fxyz = open("mgwater1.xyz","w")
  print("     4  Water + Magnesium Ion ", file=fxyz)
  print("     1  O      1.933960    0.000000    0.000000     36     2     3  ", file=fxyz)
  print("     2  H      2.518390    0.770996    0.000000    37     1        ", file=fxyz)
  print("     3  H      2.518390   -0.770996   0.000000    37     1        ", file=fxyz)
  print("%s %12.6f %s " % ("    4  Mg+", xmg ,"0.000000    0.000000    11              "), file=fxyz)
  fxyz.close()

dist = 1.5
while (dist < 5.01 ):
  set_dist(dist)
  ene = calc_ene()
  print(dist, ene)
  dist = dist + 0.05
