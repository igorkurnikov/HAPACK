# Symmetry transformation of the protein crystal
# using coordinates of the main periodic unit 
# and symmetry operations' rotational matricies and translational vectors
#
# use: Load File pdz1.pdb  that contains  coordinates
# periodic unit of the pdb file 1gq4.pdb (pdz1 domain)
#
# the script will load several times unit cell from pdz1.pdb 
# and make required symmetry transformations generating several symmetry 
# equivalent molecule positions
# translate along symmetry axes 
# change applied transformations at the end of the script 
#
import math
pmset = GetCurMolSet()
#  Set rotation matrix and transpose vector 
#  for symmetry operation 1
tmat  = []
dsp    = []
#unit cell translations
trans  = []
# unit cell pars:
a = 50.373/BOHR_TO_ANG
b = 50.373/BOHR_TO_ANG
c = 66.631/BOHR_TO_ANG  
alpha = 90.00*DEG_TO_RAD
beta  = 90.00*DEG_TO_RAD
gamma = 120.00*DEG_TO_RAD

for i in range(6):
  tmat.append(HaMat_double(3,3))
  dsp.append(HaVec_double(3,0.0))
  trans.append(HaVec_double(3,0.0))
  
#  Set rotation matrix and transpose vector 
#  for symmetry operation 1-5
tmat[0].SetFromStr( "-0.500000  0.866025  0.000000  \
			    -0.866025 -0.500000  0.000000  \
				0.000000  0.000000  1.000000")
dsp[0].SetVal_idx0(2,22.21033/BOHR_TO_ANG)
#
tmat[1].SetFromStr( "-0.500000 -0.866025  0.000000  \
			     0.866025 -0.500000  0.000000  \
			     0.000000  0.000000  1.000000")
dsp[1].SetVal_idx0(2,44.42067/BOHR_TO_ANG)
#
tmat[2].SetFromStr( "-0.500000 0.866025  0.000000  \
			     0.866025  0.500000  0.000000  \
			     0.000000  0.000000  -1.000000")
#  
tmat[3].SetFromStr( "1.000000     0.00000    0.000000  \
			    0.000000   -1.00000    0.000000  \
			    0.000000    0.00000   -1.000000")
dsp[3].SetVal_idx0(2,44.42067/BOHR_TO_ANG) 
#  
tmat[4].SetFromStr( "-0.500000 -0.866025  0.000000  \
			    -0.866025   0.500000  0.000000  \
			     0.000000  0.000000   -1.000000")
dsp[4].SetVal_idx0(2,22.21033/BOHR_TO_ANG)
#
tmat[5].SetFromStr( "1.000000     0.00000    0.000000  \
			    0.000000     1.00000    0.000000  \
			    0.000000     0.00000    1.000000")
#
umat = tmat[5]
#
#unit cell translations
trans[0].SetVal_idx0(0,a) 
trans[1].SetVal_idx0(0,-a) 
trans[2].SetVal_idx0(0,b*math.cos(gamma))
trans[2].SetVal_idx0(1,b*math.sin(gamma))
trans[3].SetVal_idx0(0,-b*math.cos(gamma))
trans[3].SetVal_idx0(1,-b*math.sin(gamma))
trans[4].SetVal_idx0(2,c) 
trans[5].SetVal_idx0(2,-c) 
#
for i in range(6):
  for j in range(4):
    pmset.LoadPDBFile("pdz1.pdb")
    nn = pmset.GetNMol() - 1
    pmol = pmset.GetMoleculeNum(nn)
    pmol.Transform(tmat[i],dsp[i])
    pmol.Transform(umat,trans[j])

pview = pmset.GetActiveMolView()
pview.InitialTransform()
pview.DefaultRepresentation()
pmset.RefreshAllViews( RFApply | RFColour | RFInitial )


