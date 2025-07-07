#sys.path.append("c:\\scripps_pylib")
from PyBabel.atomTypes import *
pmset = GetCurMolSet()
(atoms,bonds) = set_scr_atoms(pmset)
print len(atoms[2].bonds)
print atoms[1].ref
babel = AtomHybridization()
babel.assignHybridization(atoms)
from PyBabel.gasteiger import Gasteiger
Gast = Gasteiger()
Gast.compute(atoms) 
for i in  range(len(atoms)):
    at = atoms[i]
    print i,"  ",at.ref,"  ",at.gast_charge,"  ",at.babel_type 
