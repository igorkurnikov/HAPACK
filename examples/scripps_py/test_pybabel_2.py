sys.path.append("c:\\scripps_pylib")
from PyBabel.atomTypes import *
pmset = GetCurMolSet()
(atoms,bonds) = set_scr_atoms(pmset)

babel = AtomHybridization()
babel.assignHybridization(atoms)

from PyBabel.cycle import RingFinder
r = RingFinder()
r.findRings(atoms,bonds)
r.printRings()
print("num rings = ",r.ringCount)


