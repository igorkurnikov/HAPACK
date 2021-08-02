RDKIT_IMPORTED = 0
try:
    from rdkit.Chem import rdchem
    from rdkit.Chem import rdmolfiles
    import rdkit.Geometry 
    RDKIT_IMPORTED = 1
except:
    pass
#    print("rdkit is not available - some functionaliy will be turned off ")

import molset

def MolSet_to_rdkit_Mol( mset: molset.MolSet ):
  if( RDKIT_IMPORTED == 0 ): return None
  #str_pdb = mset.SavePDBToString()
  #m = rdmolfiles.MolFromPDBBlock(str_pdb, removeHs=False )
  m = rdchem.Mol()
  m_rw = rdchem.RWMol(m)
  conf = rdchem.Conformer()
  at_idx_map = mset.GetAtomSeqNumMap()
  for i,at in enumerate(mset):
#    at_idx_map[at] = i
    at_rd = rdchem.Atom(at.GetStdSymbol())
    m_rw.AddAtom(at_rd)
    conf.SetAtomPosition(i,rdkit.Geometry.Point3D(at.GetX(),at.GetY(),at.GetZ()) )
  bitr = mset.GetBondIterator()
  for bnd in bitr:
    at1 = bnd.GetFirstAtom()
    at2 = bnd.GetSecondAtom()
    idx1 = at_idx_map[at1]
    idx2 = at_idx_map[at2]
    bnd_type = rdchem.BondType.UNSPECIFIED
    if( bnd.IsSingle() ): bnd_type = rdchem.BondType.SINGLE
    if( bnd.IsDouble() ): bnd_type = rdchem.BondType.DOUBLE
    if( bnd.IsTriple() ): bnd_type = rdchem.BondType.TRIPLE
    if( bnd.IsAromatic() ): bnd_type = rdchem.BondType.AROMATIC
    if( bnd.IsVirtual() ): continue
    m_rw.AddBond(idx1,idx2, order = bnd_type)
  m_rw.AddConformer(conf)
  return m_rw


