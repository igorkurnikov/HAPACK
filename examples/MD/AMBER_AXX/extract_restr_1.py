#
# Script to exctract time dependence of values of specific distance restraints
# from AMBER restraints file
#
#########################################
# input params:
fname_in  = "1-out-NMR_win1.txt"
fname_out = "restr_len.dat"
n_restr_tot = 58          # total number of restraints
idx_restr_needed = [0,1]  # indexes of restraints needed (0-based)
screen_print = 1
#
n_restr_needed = len(idx_restr_needed)
n_lines_pt = n_restr_tot / 9 
if( n_lines_pt % 9 != 0 ):
  n_lines_pt = n_lines_pt + 1
try:
  idx_pt = 0
  fin  = open(fname_in,"r")
  fout = open(fname_out,"w")
  while(1):
    idx_pt = idx_pt + 1
    line = fin.readline()
    restr_val = []
    for i in range(n_lines_pt):
      line = fin.readline()
      tokens = line.split()
      for tok in tokens:
        restr_val.append(tok)
    if( len(restr_val) < 1 ):
      break
    if( screen_print):
      print "pt ", idx_pt," Num restr read ",len(restr_val)
      print " needed restr values ",
      for idx_c in idx_restr_needed:
        print restr_val[idx_c],"  ",
      print " "
      print >> fout,idx_pt,"  ",
    for idx_c in idx_restr_needed:
      print >> fout, restr_val[idx_c],"   ",
    print >> fout, "  "
finally:
  fin.close()
  fout.close()
print "End of script"

    
  
