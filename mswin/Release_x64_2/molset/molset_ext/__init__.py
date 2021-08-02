from .traj_utils import *
from .rdkit_utils import *
#import molset23

class MolSet(molset.molsetc.MolSet):

    def test_print(self):
        print("test_print")

    def to_mdtraj_top(self):
        return MolSet_to_mdtraj_top( self )

    def set_crd_from_mdtraj_frame(self, frame ):
        return MolSet_crd_from_frame( self , frame )



