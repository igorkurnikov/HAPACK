import os
import shutil
import tempfile
import numpy as np
import molset
import mdtraj as md

def MolSet_to_mdtraj_top( mset : molset.MolSet ):
   temp_dir = tempfile.mkdtemp(prefix='molset_temp_')
   print("temp_dir =",temp_dir)
   temp_fname = os.path.join(temp_dir,"temp.pdb")
   ires = mset.SavePDBFile(temp_fname)
   topology = md.load_pdb(temp_fname).topology
   if(os.path.exists(temp_dir)): shutil.rmtree(temp_dir)
   return topology

def MolSet_crd_from_frame( mset : molset.MolSet, t : md.Trajectory ):
    print(t.time)
    for i,at in enumerate(mset):
        print("i =",i)
        at.SetX(float(t.xyz[0,i,0]))
        at.SetY(float(t.xyz[0,i,1]))
        at.SetZ(float(t.xyz[0,i,2]))
    aitr = molset.AtomIteratorMolSet(mset)
    at = aitr.GetFirstAtom()
    print(at.GetX(), at.GetY(), at.GetZ())


