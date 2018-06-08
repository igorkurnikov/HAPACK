
source /etc/profile

export PATH=${HOME}/HARLEM/bin:${PATH}

cd $PBS_O_WORKDIR

harlem_nogui_mpi -n 8 mg_wat_1.hlm -script run_md_1.py 

