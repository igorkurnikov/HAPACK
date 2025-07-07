This directory should contain next files:

runMC_Empirical.py
basinenergy.py
empir_param.dat
equilibrated.hlm
min1.pdb
min2.pdb
min3.pdb
min4.pdb
min5.pdb
min6.pdb
bundle50A.hlm


1. First, load equilibrated.hlm into HARLEM and run runMC_Empirical.py script.
2. Next, load and run basinenergy.py
3. It is going to generate next files:

rmsd_state(i).out for i=1,2,...6
density_state(i).dat for i=1,2,...6
eff_ene_file.dat
minimum(i).hlm for i=1,2,...6
all_mc_pts.dat  accepted trial conformations
