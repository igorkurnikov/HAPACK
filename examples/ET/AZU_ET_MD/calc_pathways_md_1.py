if(script_status == SCRIPT_START):
  pmset = GetCurMolSet()
  etmod   = pmset.GetETCouplMod(1)
  etmod.pathways_calc_type = BEST_PATH
  outf = open("pathways_md_1.dat","w")
  ns = 0
  print(" time   pathways_hda", file=outf)
elif(script_status == SCRIPT_STOP):
#  After MD trajectory analysis - close output file
  outf.close()
else:
  if( ns % 10 == 0):
    etmod.InitiatePathwaysGraph()
    etmod.path_coupl_calc()
    print(ns, etmod.best_path_coupl, file=outf)
  ns = ns+1


