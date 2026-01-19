"""  Utilities to interact with mdtraj package """

import os
import shutil
import tempfile
import re

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
    """ Convert MolSet to mdtraj.Topology """
    if( MDTRAJ_IMPORTED == 0 ): return None
    temp_dir = tempfile.mkdtemp(prefix='molset_temp_')
    print("temp_dir =",temp_dir)
    temp_fname = os.path.join(temp_dir,"temp.pdb")
    ires = mset.SavePDBFile(temp_fname)
    topology = md.load_pdb(temp_fname).topology
    if(os.path.exists(temp_dir)): shutil.rmtree(temp_dir)
    return topology

def MolSet_crd_from_frame( mset : molset.MolSet, t : md.Trajectory ):
    """ set coordinates of Molset from mdtraj trajectory frame """
    if( MDTRAJ_IMPORTED == 0 ): return None
    print(t.time)
    for i,at in enumerate(mset):
        at.SetX(float(t.xyz[0,i,0])*10.0)
        at.SetY(float(t.xyz[0,i,1])*10.0)
        at.SetZ(float(t.xyz[0,i,2])*10.0)
        
def get_at_idx(top: md.Topology, resi: int, at_name: str):
    """ 
    get atom index (0-based) in mdtraj Topology 
    by residue index (0-based) and atom name 
    """
    if( MDTRAJ_IMPORTED == 0 ): return -1
    at_arr = [at for at in top.residue(resi-1).atoms_by_name(at_name)]
    if( len(at_arr) == 0 ): return -1
    return at_arr[0].index

def get_atidx_mdtraj(top : md.Topology, atid : str):
    at_idx = 0
    for i,atom in enumerate(top.atoms):
        if str(atom) == atid: 
            at_idx = i
            break
    return at_idx

def get_atidx_mdtraj_m(top : md.Topology, atid : str):
    # This does not work when there is a missing residue
    # as resid assumes order number of the residue
    #
    tokens = atid.split('.')
    resid_f = tokens[0]
    atname = tokens[1]
    match = re.search(r'\d+', resid_f)
    resid = int(match.group()) - 1
    at_idx = top.select(f"resid {resid} and name {atname}")[0]
    return at_idx

def get_atom_pair_idx_mdtraj(top : md.Topology, atpair_str: str):
    tokens = atpair_str.strip().split(" ")
    at_idx_0 = get_atidx_mdtraj(top,tokens[0])
    at_idx_1 = get_atidx_mdtraj(top,tokens[1])
    return([at_idx_0,at_idx_1])
    
def get_dist_at_pairs( trj, atpair_str_list ):
    atom_pairs = []
    for s in atpair_str_list:
        at_pair = get_atom_pair_idx_mdtraj(trj.top,s )
        atom_pairs.append(at_pair)
    distances = md.compute_distances(trj,atom_pairs)*10.0
    return distances

def make_pair_idx( top: md.Topology,  atid_arr : list[str] ):
    """ 
    make atom pair index array ( 0-based )  
    from string array of atom ids in RASMOL format: Trp29.CG 
    for mdtraj.Topology  
    
    return: numpy array atom pair indexes (dtype=int) shape (-1,2)
    """    
    pair_idx = np.empty(0,dtype=int)
    col_names = []
    for i in range( len(atid_arr)//2 ):
        atid_1 = atid_arr[2*i]
        atid_2 = atid_arr[2*i+1]
        tokens1 = atid_1.split(".")
        nr1 = int( ''.join(re.findall(r'\d+',tokens1[0])) ) 
        atn1 = tokens1[1]
        tokens2 = atid_2.split(".")
        nr2 = int( ''.join(re.findall(r'\d+',tokens2[0])) )
        atn2 = tokens2[1]
        pair_idx = np.append( pair_idx,get_at_idx(top,nr1,atn1) ) 
        pair_idx = np.append( pair_idx,get_at_idx(top,nr2,atn2) )  
        col_names.append(atid_1 + " - " + atid_2)
    pair_idx = np.reshape(pair_idx,(-1,2))
    return(pair_idx,col_names)

def make_dihedral_idx( top: md.Topology, atid_arr : list[str] ):
    """ 
    make atom quadruplets index array (0-based) 
    from string array of atom id in RASMOL format: Trp29.CG 
    for mdtraj.Topology  

    return: numpy array atom pair indexes (dtype=int) shape (-1,4)
    """   
    dihedral_idx = np.empty(0,dtype=int)
    col_names = []
    for i in range( len(atid_arr)//4 ):
        atid_1 = atid_arr[4*i]
        atid_2 = atid_arr[4*i+1]
        atid_3 = atid_arr[4*i+2]
        atid_4 = atid_arr[4*i+3]
        tokens1 = atid_1.split(".")
        nr1 = int( ''.join(re.findall(r'\d+',tokens1[0])) ) 
        atn1 = tokens1[1]
        tokens2 = atid_2.split(".")
        nr2 = int( ''.join(re.findall(r'\d+',tokens2[0])) )
        atn2 = tokens2[1]
        tokens3 = atid_3.split(".")
        nr3 = int( ''.join(re.findall(r'\d+',tokens3[0])) )
        atn3 = tokens3[1]
        tokens4 = atid_4.split(".")
        nr4 = int( ''.join(re.findall(r'\d+',tokens4[0])) )
        atn4 = tokens4[1]
        dihedral_idx = np.append( dihedral_idx,get_at_idx(top,nr1,atn1) ) 
        dihedral_idx = np.append( dihedral_idx,get_at_idx(top,nr2,atn2) )  
        dihedral_idx = np.append( dihedral_idx,get_at_idx(top,nr3,atn3) ) 
        dihedral_idx = np.append( dihedral_idx,get_at_idx(top,nr4,atn4) ) 
        col_names.append(atid_1 + " - " + atid_2 + " - " + atid_3 + " - " + atid_4)
    dihedral_idx = np.reshape(dihedral_idx,(-1,4))
    return(dihedral_idx,col_names)






