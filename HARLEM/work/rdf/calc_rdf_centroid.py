import mdtraj as md
import molset
import numpy as np

mset = molset.MolSet()
mset.LoadNRGFile("config.nrg")
top_m = mset.to_mdtraj_top()
# top_m = traj_util.MolSet_to_mdtraj_top(mset)

ta = []
ta.append( md.load("n256_npt_equil_100ps_R4_0.trr", top=top_m) )
ta.append( md.load("n256_npt_equil_100ps_R4_1.trr", top=top_m) )
ta.append( md.load("n256_npt_equil_100ps_R4_2.trr", top=top_m) )
ta.append( md.load("n256_npt_equil_100ps_R4_3.trr", top=top_m) )

t = ta[0]
t.xyz = np.add( t.xyz, ta[1].xyz )
t.xyz = np.add( t.xyz, ta[2].xyz )
t.xyz = np.add( t.xyz, ta[3].xyz )

np.multiply(t.xyz,0.25)

pairs_o = top_m.select_pairs("name O","name O")
(r,gr) = md.compute_rdf(t,pairs_o)

n = len(r)
fout = open("gr_o_r4_centr_mdtraj.dat","w")
for i in range(len(r)):
    print(f" {r[i]:6.3f} {gr[i]:9.5f} ")
    fout.write(f" {r[i]:6.3f} {gr[i]:9.5f} \n")
fout.close()

