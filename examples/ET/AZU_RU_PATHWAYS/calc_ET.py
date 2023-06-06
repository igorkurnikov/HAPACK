#
# Compute Electron transfer coupling between chosen donor and acceptor atoms 
#
import sys
import os
import molset

if( len(sys.argv) < 4): 
    print(f"usage: {sys.argv[0]} protein.pdb  donor_expr  acceptor_expr")
    exit(0)
    
fname_pdb      = sys.argv[1]
donor_expr    = sys.argv[2]
acceptor_expr = sys.argv[3]

mset = molset.MolSet()
mset.LoadPDBFile(fname_pdb)
etmod = mset.GetETCouplMod(True)

mset.SelectAtomsExpr( donor_expr )
donor_grp = mset.SetAtomGroupFromSelection("DONOR")

if( donor_grp.GetNAtoms() == 0 ):
    print(f"Can not find donor atoms for expression {donor_expr} ")
    exit(1)
    
mset.SelectAtomsExpr( acceptor_expr )
acceptor_grp = mset.SetAtomGroupFromSelection("ACCEPTOR")

if( acceptor_grp.GetNAtoms() == 0 ):
    print(f"Can not find donor atoms for expression {acceptor_expr} ")
    exit(1)
    
    
etmod.path_coupl_calc()



