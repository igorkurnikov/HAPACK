fname_log = "pept_1_box.out"
fname_ene = "ene.out"
par_names = ["NSTEP","TIME(PS)","Etot","EPtot","VOLUME"]
info = {}

def read_next_pt(flog):
  in_pt = 0
  flag_avg = 0
  flag_rms = 0
  info.clear()
  for line in flog:
    if(len(line) == 0): 
      return -1
    if(line.find("A V E R A G E S   O V E R")  != -1):
      flag_avg = 1
    if(line.find("R M S  F L U C T U A T I O N S")  != -1):
      flag_rms = 1
    if(line.find("NSTEP =") != -1):
      in_pt = 1	
    if(in_pt == 1): 
      if(line.find("---------") > 0 ):
        if( flag_avg == 1 or flag_rms == 1):
	  return 2
        return 1
      words = line.split()
      n = len(words)
      i = 0
      while( i < n):
        if( words[i] in par_names and (i+2) < n ):
          info[words[i]] = words[i+2]
          i = i + 2
        i = i + 1   
  return -1
  
flog  = open(fname_log,"r")
fene = open(fname_ene,"w")
ires = 1
nnm = len(par_names)
for i in range(nnm):
  name = par_names[i]
  print name, "   ",
  print >> fene, name, "   ",
print " "
print >> fene," "
      
while( ires > 0):
  ires = read_next_pt(flog)
  if( ires > 0 ):
    if( ires == 2):
      continue
    for i in range(nnm):
      name = par_names[i]
      val = "99999"
      if(info.has_key(name)):
        val = info[name]
      print val, "  ",
      print >> fene, val, "   ",
    print " "
    print >> fene, "   "
flog.close()
fene.close()
