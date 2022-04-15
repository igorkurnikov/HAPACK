import mdtraj as md
import molset
import traj_utils

mset = molset.MolSet()
mset.LoadNRGFile("config.nrg")
top_m = traj_utils.MolSet_to_mdtraj_top(mset)

t = md.load_xyz("n256_npt_equil_100ps.pos_0.xyz", top=top_m)
t.save("n256_npt_equil_100ps_R4_0.trr")

t = md.load_xyz("n256_npt_equil_100ps.pos_1.xyz", top=top_m)
t.save("n256_npt_equil_100ps_R4_1.trr")

t = md.load_xyz("n256_npt_equil_100ps.pos_2.xyz", top=top_m)
t.save("n256_npt_equil_100ps_R4_2.trr")

t = md.load_xyz("n256_npt_equil_100ps.pos_3.xyz", top=top_m)
t.save("n256_npt_equil_100ps_R4_3.trr")



