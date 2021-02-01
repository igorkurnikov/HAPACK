import os
import shutil
import tempfile
NUMPY_IMPORTED = 0
MDTRAJ_IMPORTED = 0
try:
    import numpy as np
    NUMPY_IMPORTED = 1
except:
    pass
#    print("numpy is not available - some functionaliy will be turned off ")

try:
    import mdtraj as md
    MDTRAJ_IMPORTED = 1
except:
    pass
#    print("mdtraj is not available - some functionaliy will be turned off ")

import molset

def MolSet_to_mdtraj_top( mset : molset.MolSet ):
    if( MDTRAJ_IMPORTED == 0 ): return None
    temp_dir = tempfile.mkdtemp(prefix='molset_temp_')
    print("temp_dir =",temp_dir)
    temp_fname = os.path.join(temp_dir,"temp.pdb")
    ires = mset.SavePDBFile(temp_fname)
    topology = md.load_pdb(temp_fname).topology
    if(os.path.exists(temp_dir)): shutil.rmtree(temp_dir)
    return topology

def MolSet_crd_from_frame( mset : molset.MolSet, t ):
    if( MDTRAJ_IMPORTED == 0 ): return None
    print(t.time)
    for i,at in enumerate(mset):
        at.SetX(float(t.xyz[0,i,0])*10.0)
        at.SetY(float(t.xyz[0,i,1])*10.0)
        at.SetZ(float(t.xyz[0,i,2])*10.0)



