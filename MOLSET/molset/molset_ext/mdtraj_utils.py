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

def load_trajectory(fname, top):
    """
    Load an MD trajectory via mdtraj, tolerating AMBER NetCDF ('coordinates'
    variable) as well as OpenMM/OpenFE NetCDF ('positions' variable).

    Returns an mdtraj.Trajectory. Raises on unrecoverable failure.
    """
    if MDTRAJ_IMPORTED == 0:
        raise RuntimeError("mdtraj is not available")

    # First attempt: mdtraj's native auto-dispatch.
    try:
        return md.load(fname, top=top)
    except (KeyError, OSError, ValueError, TypeError) as exc:
        # Only retry for NetCDF-like files; for other failures, rethrow.
        ext = os.path.splitext(fname)[1].lower()
        if ext not in ('.nc', '.ncdf', '.netcdf'):
            raise

        # OpenMM NetCDF fallback: variable 'positions' in nanometers.
        try:
            import netCDF4
        except ImportError:
            raise RuntimeError(
                "Failed to load NetCDF via mdtraj AMBER path (KeyError 'coordinates') "
                "and the 'netCDF4' Python package is not installed for OpenMM fallback. "
                "Install with: pip install netCDF4"
            ) from exc

        d = netCDF4.Dataset(fname, 'r')
        try:
            if 'positions' in d.variables:
                pos_var = d.variables['positions']
                fill_value = getattr(pos_var, '_FillValue', None)
                print("NetCDF 'positions' shape:", pos_var.shape, "dims:", pos_var.dimensions,
                      "fill_value:", fill_value)
                xyz = np.array(pos_var[:], dtype=np.float32)
                # OpenMM multi-state / alchemical reporters produce 4D:
                # (n_iterations, n_replicas, n_atoms, 3). Collapse to 3D by
                # picking replica 0.
                if xyz.ndim == 4:
                    replica_index = 0
                    print(f"4D positions detected; selecting replica index {replica_index} "
                          f"(of {xyz.shape[1]} replicas)")
                    xyz = xyz[:, replica_index, :, :]
                elif xyz.ndim != 3:
                    raise RuntimeError(
                        f"Unexpected 'positions' ndim={xyz.ndim}, shape={xyz.shape}"
                    )

                # Filter out uninitialized frames (NetCDF checkpoint mode stores
                # positions only at some iterations; the rest are fill value
                # ~9.97e36). Detect frames where all coords are at or near fill.
                fill_threshold = 1.0e30  # anything above this is fill
                per_frame_max = np.abs(xyz).max(axis=(1, 2))
                valid_mask = per_frame_max < fill_threshold
                n_valid = int(valid_mask.sum())
                n_total = xyz.shape[0]
                if n_valid == 0:
                    raise RuntimeError(
                        f"All {n_total} frames in the NetCDF trajectory are unwritten "
                        f"fill-value frames. This often happens with OpenFE checkpoint-mode "
                        f"files: the 'positions' variable is allocated for all iterations "
                        f"but populated only at checkpoint iterations. Try a trajectory "
                        f"saved in full (non-checkpoint) mode."
                    )
                if n_valid < n_total:
                    print(f"Filtering fill-value frames: keeping {n_valid} of {n_total} "
                          f"actual populated frames")
                    xyz = xyz[valid_mask]
                # Stash the valid-frames mask for time/cell filtering below.
                _valid_mask = valid_mask
            elif 'coordinates' in d.variables:
                # AMBER-style under an exotic schema that mdtraj rejected.
                xyz = np.array(d.variables['coordinates'][:], dtype=np.float32)
                xyz = xyz / 10.0  # AMBER is Angstrom, mdtraj expects nm
                _valid_mask = None
            else:
                raise RuntimeError(
                    "NetCDF file %s has neither 'positions' nor 'coordinates' "
                    "variable. Found: %s" % (fname, list(d.variables.keys()))
                )

            # time (collapse multi-state dim if present, then apply valid mask)
            if 'time' in d.variables:
                time = np.array(d.variables['time'][:], dtype=np.float32)
                if time.ndim == 2:
                    time = time[:, 0]
                if _valid_mask is not None and len(time) == len(_valid_mask):
                    time = time[_valid_mask]
            else:
                time = np.arange(len(xyz), dtype=np.float32)

            # box (optional) — strip state dim if present, apply valid mask
            cell_lengths = None
            cell_angles = None
            if 'cell_lengths' in d.variables:
                cell_lengths = np.array(d.variables['cell_lengths'][:], dtype=np.float32)
                if cell_lengths.ndim == 3:
                    cell_lengths = cell_lengths[:, 0, :]
                if _valid_mask is not None and len(cell_lengths) == len(_valid_mask):
                    cell_lengths = cell_lengths[_valid_mask]
                prog = getattr(d, 'program', '')
                if 'amber' in prog.lower():
                    cell_lengths = cell_lengths / 10.0
            if 'cell_angles' in d.variables:
                cell_angles = np.array(d.variables['cell_angles'][:], dtype=np.float32)
                if cell_angles.ndim == 3:
                    cell_angles = cell_angles[:, 0, :]
                if _valid_mask is not None and len(cell_angles) == len(_valid_mask):
                    cell_angles = cell_angles[_valid_mask]
        finally:
            d.close()

        trj = md.Trajectory(xyz, top, time=time)
        if cell_lengths is not None and cell_angles is not None:
            trj.unitcell_lengths = cell_lengths
            trj.unitcell_angles = cell_angles
        return trj

def MolSet_crd_from_frame( mset : molset.MolSet, t : md.Trajectory ):
    """ set coordinates of Molset from mdtraj trajectory frame """
    if( MDTRAJ_IMPORTED == 0 ): return None
    print("MolSet_crd_from_frame: frame time =", t.time, "xyz shape =", t.xyz.shape)
    # Diagnostic sample of actual data
    print(f"  first 3 atom coords (nm): {t.xyz[0, :3, :]}")
    print(f"  xyz min/max/mean: {t.xyz.min():.4f} / {t.xyz.max():.4f} / {t.xyz.mean():.4f}")
    n_mset = sum(1 for _ in mset)
    n_traj = t.xyz.shape[1]
    print(f"  MolSet atoms = {n_mset}, trajectory atoms = {n_traj}")
    if n_traj < n_mset:
        print(f"  WARNING: trajectory has fewer atoms than MolSet ({n_traj} < {n_mset}); "
              f"leaving trailing atoms unchanged")
    elif n_traj > n_mset:
        print(f"  INFO: trajectory has extra atoms ({n_traj} > {n_mset}); "
              f"using first {n_mset}")
    n_use = min(n_mset, n_traj)
    for i, at in enumerate(mset):
        if i >= n_use:
            break
        at.SetX(float(t.xyz[0, i, 0]) * 10.0)
        at.SetY(float(t.xyz[0, i, 1]) * 10.0)
        at.SetZ(float(t.xyz[0, i, 2]) * 10.0)
        
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






