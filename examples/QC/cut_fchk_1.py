import sys
if( len(sys.argv) != 3 ):
  print "usage python " + sys.argv[0] + " file.fchk " + " short_file.fchk"
  sys.exit()
fname_in = sys.argv[1]
fname_out = sys.argv[2]
finp = open(fname_in,"r")
fout = open(fname_out,"w")
sections_yes = ["Number of atoms","Info1-9","Charge","Multiplicity", \
            "Number of electrons","Number of basis functions", \
                "Atomic numbers","Number of independent functions"]
sections_no = ["Total SCF Density"]
nread = 2
write_line = True
nl = 0
for line in finp:
  nl = nl + 1
  if( nl > 20000000 ):
    continue
  line.rstrip("\n")
  if( nread > 0 ):
    if( write_line ): print >> fout, line,
    nread = nread - 1
    continue
  print line,
  name_s = line[:43].strip()
  type = line[43]
  ne_s = line[47:49] 
  if( ne_s == "N="):
    snum = line[51:61].strip()
    num = int(snum)
    if( type == "I"):
      nread = num/6
      if( num % 6 != 0):
        nread = nread + 1
    if( type == "R"):
      nread = num/5
      if( num % 5 != 0):
        nread = nread + 1
  print "nread = ",nread
  if( name_s not in sections_no ):
    write_line = True
  else:
    write_line = False
  if( write_line): print >> fout, line,

  
    
