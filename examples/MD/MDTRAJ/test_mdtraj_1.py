import os
import mdtraj as md
from molset import MolSet

mset = MolSet()
print(os.getcwd())
path_script = os.path.dirname(os.path.realpath(__file__))
os.chdir(path_script)
print(os.getcwd())
fpath = os.path.join("..","mg_wat_1.hlm")
mset.LoadHarlemFile(fpath)
topology = mset.to_mdtraj_top()
print(topology)
fpath_traj = os.path.join("..","mg_wat_1.mdcrd")
for frame in md.iterload(fpath_traj, top=topology, chunk = 1, stride = 50):
    print(frame)
    mset.set_crd_from_mdtraj_frame(frame)


print("Done")


