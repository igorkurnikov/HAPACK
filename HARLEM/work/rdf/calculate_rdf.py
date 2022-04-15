import mdtraj as md
import molset
import traj_utils

mset = molset.MolSet()
mset.LoadNRGFile("config.nrg")
top_m = traj_utils.MolSet_to_mdtraj_top(mset)

t = md.load("R4_b25_all.trr", top=top_m)
pairs_o = top_m.select_pairs("name O","name O")

(r,gr) = md.compute_rdf(t,pairs_o)
n = len(r)
fout = open("gr_o_r4_mdtraj.dat","w")
for i in range(len(r)):
    print(f" {r[i]:6.3f} {gr[i]:9.5f} ")
    fout.write(f" {r[i]:6.3f} {gr[i]:9.5f} \n")
fout.close()

t = md.load("md.trr", top=top_m)
(r,gr) = md.compute_rdf(t,pairs_o)
n = len(r)
fout = open("gr_o_md_mdtraj.dat","w")
for i in range(len(r)):
    print(f" {r[i]:6.3f} {gr[i]:9.5f} ")
    fout.write(f" {r[i]:6.3f} {gr[i]:9.5f} \n")
fout.close()



