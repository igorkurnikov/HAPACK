#
# Classes to interact between Scripps py_babel and HARLEM
#
#import PyBabel
#from PyBabel import *
#import molset
#from molset import *
#import sys
#sys.path.append("c:\\scripps_pylib")
#sys.path.append("c:\\python20")
#sys.path.append("c:\\scripps_pylib\\PyBabel")

import PyBabel.atomTypes
from PyBabel.atomTypes import AtomHybridization

import sys

class scr_atom:    
    def __init__(self):
        self.coords = [0.0,0.0,0.0]
        self.bonds  = []
        self.element = "X"
        self.ref =  ""
	self.name = "atom"
	self.number = 0

class scr_bond:
    def __init__(self):
        self.atom1 = scr_atom()
        self.atom2 = scr_atom()

def set_scr_atoms(pmset):
    scr_atoms = []
    scr_bonds = []
    aitr = AtomIteratorMolSet(pmset)
    aptr = aitr.GetFirstAtom()
    idx = 0
    at_idx_map = {}
    while (aptr != None):
        satom1 = scr_atom()
        scr_atoms.append(satom1)
        scr_atoms[idx].coords = [aptr.GetX(),aptr.GetY(),aptr.GetZ()]
	elemno = aptr.GetElemNo()
	scr_atoms[idx].element = (aptr.GetStdSymbol()).c_str()
	atref = aptr.GetRef()
	scr_atoms[idx].name = atref
	scr_atoms[idx].ref = atref
	scr_atoms[idx].number = idx+1
	at_idx_map[atref] = idx
	idx = idx + 1
	aptr = aitr.GetNextAtom()
	
    bitr = BondIteratorMolSet(pmset)
    bptr = bitr.GetFirstBond()
    while(bptr != None):
       at1= bptr.srcatom
       atref1 = at1.GetRef()
       at2= bptr.dstatom
       atref2 = at2.GetRef()
#       print "bond: ",atref1,"-",atref2
       idx1=  at_idx_map[atref1]
       idx2 = at_idx_map[atref2]
#       print "bond: ",idx1,"-",idx2
       sat1 = scr_atoms[idx1]
       sat2 = scr_atoms[idx2]
       sbond = scr_bond()
       sbond.atom1 = sat1
       sbond.atom2 = sat2
       scr_bonds.append(sbond)
       scr_atoms[idx1].bonds.append(sbond)
       scr_atoms[idx2].bonds.append(sbond)
       bptr = bitr.GetNextBond()
    res = []
    res.append(scr_atoms);
    res.append(scr_bonds);    
    return (scr_atoms,scr_bonds)
	
    