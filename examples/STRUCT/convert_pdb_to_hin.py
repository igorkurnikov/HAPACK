#!/usr/bin/env python3
import sys
if(len(sys.argv) < 3):
  print("usage: " + sys.argv[0] + "  inp.pdb  out.hin ")
  exit(1)
from molset import *
fn_inp = sys.argv[1]
fn_out = sys.argv[2]
mset = HaMolSet()
mset.LoadPDBFile(fn_inp)
mset.SaveHINFile(fn_out)
print("Done!")
