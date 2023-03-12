import parmed as pmd
amber = pmd.load_file('NRAS_15H_4.top', 'NRAS_15H_4.inp')
#save a GROMACS topology and GRO file
amber.save('NRAS_15H_4_gmx.top')
amber.save('NRAS_15H_4.gro')

