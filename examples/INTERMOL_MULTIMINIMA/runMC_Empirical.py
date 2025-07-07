pmset = GetCurMolSet()
empmod = pmset.GetEmpiricalMod(1)
intmod = pmset.GetInterMolMod(1)
empmod.Initialize()
intmod.Initialize()
intmod.tr_ratio =3.0
intmod.ang_ratio = 0.064
intmod.empirical_flag=1
intmod.num_mc_steps = 1000
intmod.RunMCQuantSampling()
