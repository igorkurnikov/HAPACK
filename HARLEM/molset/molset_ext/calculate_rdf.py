import mdtraj as md
import molset
import traj_utils

mset = molset.MolSet()
mset.LoadNRGFile("config.nrg")
top_m = traj_utils.MolSet_to_mdtraj_top(mset)

t = md.load_xyz("n256_npt_equil_100ps.pos_0.xyz", top=top_m)
pairs_o = top_m.select_pairs("name O1","name O1")
(r,g_r) = md.compute_rdf(t,pairs=pairs_o)
print(g_r)
