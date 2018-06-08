import os
os.system("del restart_MC.dat")
pmset = GetCurMolSet()
##empmod = pmset.GetEmpiricalMod(1)
intmod = pmset.GetInterMolMod(1)
##empmod.Initialize()
intmod.Initialize()
intmod.tr_ratio= 0.95
intmod.ang_ratio = 0.01 #1.0
intmod.empirical_flag=1
intmod.num_mc_steps =25
##intmod.RunMCQuantSampling()

## Makes non-uniform stepsize
intmod.nreplicas =1
intmod.rem_steps = 4000
intmod.rem_flag =1
##intmod.MC_temp = 250
intmod.temperature_max =300.0
##intmod.xy_mc_flag = 1 # Run RunQuasiREM only in xy plane
intmod.RunQuasiREM()
print "###END####"
