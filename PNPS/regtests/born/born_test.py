import pnps
import os
import math


def theoretic_born_enery(diel_from, diel_to, q, R):
    BOHR_TO_ANG = 0.529177249  # Bohrs to Angstroms
    HARTREE_TO_KT = 1059.642  # HARTREE TO kT at 298K

    de = -0.5 * HARTREE_TO_KT * (1 / diel_from - 1 / diel_to) * q * q / (R / BOHR_TO_ANG)
    return de


def solve_p(diel_ion, diel_sol, x=0.0, y=0.0, z=0.0, q=1.0, R=2.0, GridSize=64, BC="CoulBC", PBC=None):
    if PBC is None:
        PBC = [False, False, False]
    # Create Continuum representation
    contworld = pnps.ContWorld(
        GridSize=[GridSize, GridSize, GridSize],
        GridScale=4.0,
        Qions=[1, -1],
        PBC=PBC
    )

    # setup builder
    builder = pnps.GetWorldBuilder(
        DielConstBulk=diel_sol,
        Cbulk=0.0,
        DiffusionModel="Plane",
        BoundaryCondition=BC,
        MakeDielectricMap=True,
        MakeChargeMap=True,
        MakeDiffusionMap=False,
        MakeConcentrationMap=False,
        MakeSoftRepultionMap=False
    )
    # add atoms
    builder.addAtoms(
        DielConst=diel_ion,
        AtomsPQR=[[x, y, z, q, R]]
    )
    # Do building
    builder.BuildContWorld(contworld)

    # Write created maps to file system
    # contworld.WriteNodeIndexing("S1_System.systop")
    # Clean-up
    del builder

    pnps.SolveP(contworld,
                MaxIterations=-2,
                Convergence=0.0,
                Relaxation=1.0)

    # contworld.WritePotential("S2_Pot.bin")
    E = contworld.SystemEnergy

    del contworld
    return E


def test_born_solvation_energy():
    diel_ion = 2.0
    diel_sol = 80.0
    x = 0.0
    y = 0.0
    z = 0.0
    q = 1.0
    R = 2.0
    grid_scale = 4  # 2 grid/A
    n = 3

    # dEsim = solve_p(diel_ion,diel_sol,x, y, z, q, R)-solve_p(diel_ion,diel_ion,x, y, z, q, R)

    dEsim = 0.0
    dEsimSQ = 0.0
    h = 1.0 / grid_scale
    tot_num = 0
    for ix in range(n):
        for iy in range(ix + 1):
            for iz in range(iy + 1):
                x = ix * h / (n + 1.0)
                y = iy * h / (n + 1.0)
                z = iz * h / (n + 1.0)
                dE = solve_p(diel_ion, diel_sol, x, y, z, q, R, 64) - solve_p(diel_ion, diel_ion, x, y, z, q, R, 64)
                dEsim += dE
                dEsimSQ += dE * dE
                tot_num += 1

    dEsim = dEsim / tot_num
    dEsimSQ = dEsimSQ / tot_num
    sd = dEsimSQ - dEsim * dEsim
    sd = math.sqrt(abs(sd))
    dEtheor = theoretic_born_enery(diel_ion, diel_sol, q, R)

    print("Simulated solvation energy:       %.3f +- %.3f kT" % (dEsim, sd))
    print("Theoretical solvation energy:     %.3f" % dEtheor)
    abs_diff = dEsim - dEtheor
    print("Absolute difference (sim-theory): %.3f" % abs_diff)
    rel_diff = abs_diff / dEtheor
    print("Relative difference (sim-theory): %.3f" % rel_diff)

    ref_dEsim = -69.12774138940426
    ref_sd = 0.14783895743446462
    ref_abs_diff = -0.7877441990622032
    ref_rel_diff = 0.011526839792927715

    def abs_rel_diff(val, ref):
        return abs((val - ref) / ref)

    print("Difference with Reference:")
    print("\tSimulated solvation energy:       %.3e" % (dEsim-ref_dEsim))
    print("\tSimulated solvation energy, StDev:       %.3e" % (sd - ref_sd))
    print("Relative difference with Reference:")
    print("\tSimulated solvation energy:       %.3e" % ((dEsim - ref_dEsim)/ref_dEsim))
    print("\tSimulated solvation energy, StDev:       %.3e" % ((sd - ref_sd)/ref_sd))


    print(dEsim)
    print(sd)
    print(abs_diff)
    print(rel_diff)
    print(abs(sd-ref_sd)/ref_dEsim)
    print(abs(abs_diff - ref_abs_diff) / ref_dEsim)

    assert abs_rel_diff(dEsim, ref_dEsim) < 1.0e-6
    assert abs((sd-ref_sd)/ref_dEsim) < 1.0e-6
    assert abs(abs_diff - ref_abs_diff) / ref_dEsim < 1.0e-6


def test_born_solvation_energy_oddgrid():
    diel_ion = 2.0
    diel_sol = 80.0
    x = 0.0
    y = 0.0
    z = 0.0
    q = 1.0
    R = 2.0
    grid_scale = 4  # 2 grid/A
    n = 3

    # dEsim = solve_p(diel_ion,diel_sol,x, y, z, q, R)-solve_p(diel_ion,diel_ion,x, y, z, q, R)

    dEsim = 0.0
    dEsimSQ = 0.0
    h = 1.0 / grid_scale
    tot_num = 0
    for ix in range(n):
        for iy in range(ix + 1):
            for iz in range(iy + 1):
                x = ix * h / (n + 1.0)
                y = iy * h / (n + 1.0)
                z = iz * h / (n + 1.0)
                dE = solve_p(diel_ion, diel_sol, x, y, z, q, R, 65) - solve_p(diel_ion, diel_ion, x, y, z, q, R, 65)
                dEsim += dE
                dEsimSQ += dE * dE
                tot_num += 1

    dEsim = dEsim / tot_num
    dEsimSQ = dEsimSQ / tot_num
    sd = dEsimSQ - dEsim * dEsim
    sd = math.sqrt(abs(sd))
    dEtheor = theoretic_born_enery(diel_ion, diel_sol, q, R)

    print("Simulated solvation energy:       %.3f +- %.3f kT" % (dEsim, sd))
    print("Theoretical solvation energy:     %.3f" % dEtheor)
    abs_diff = dEsim - dEtheor
    print("Absolute difference (sim-theory): %.3f" % abs_diff)
    rel_diff = abs_diff / dEtheor
    print("Relative difference (sim-theory): %.3f" % rel_diff)

    ref_dEsim = -69.12786293519918
    ref_sd = 0.1477743546953387
    ref_abs_diff = -0.787865744857
    ref_rel_diff = 0.0115286183384

    def abs_rel_diff(val, ref):
        return abs((val - ref) / ref)

    print("Difference with Reference:")
    print("\tSimulated solvation energy:       %.3e" % (dEsim-ref_dEsim))
    print("\tSimulated solvation energy, StDev:       %.3e" % (sd - ref_sd))
    print("Relative difference with Reference:")
    print("\tSimulated solvation energy:       %.3e" % ((dEsim - ref_dEsim)/ref_dEsim))
    print("\tSimulated solvation energy, StDev:       %.3e" % ((sd - ref_sd)/ref_sd))


    print(dEsim)
    print(sd)
    print(abs_diff)
    print(rel_diff)
    print(abs(sd-ref_sd)/ref_dEsim)
    print(abs(abs_diff - ref_abs_diff) / ref_dEsim)

    assert abs_rel_diff(dEsim, ref_dEsim) < 1.0e-6
    assert abs((sd-ref_sd)/ref_dEsim) < 1.0e-6
    assert abs(abs_diff - ref_abs_diff) / ref_dEsim < 1.0e-6

def test_born_solvation_energy_pbc():
    """
    Check PBC by symmetry
    """
    pnpsapp = pnps.PNPSApp.GetPNPSApp()
    pnpsapp.SetNumOfThreads(1)

    diel_ion = 2.0
    diel_sol = 80.0
    x = 0.0
    y = 0.0
    z = 0.0
    q = 1.0
    R = 2.0
    grid_scale = 4  # grid/A

    PBS_List = [
        [True, False, False],
        [False, True, False],
        [False, False, True],
        [True, True, False],
        [False, True, True],
        [True, False, True],
    ]

    E_List = [PBC+[solve_p(diel_ion, diel_sol, x, y, z, q, R, 64, "ZeroBC", PBC[:3])] for PBC in PBS_List]

    print("%6s %6s %6s %20s" % ("PBC_X", "PBC_Y", "PBC_Z", "E"))
    for x, y, z, E in E_List:
        print("%6s %6s %6s %20.7f" % (str(x), str(y), str(z), E))

    assert abs(E_List[0][3] - E_List[1][3]) < 1e-3
    assert abs(E_List[0][3] - E_List[2][3]) < 1e-3
    assert abs(E_List[3][3] - E_List[4][3]) < 1e-3
    assert abs(E_List[3][3] - E_List[5][3]) < 2e-3


if __name__ == "__main__":
    test_born_solvation_energy()
    test_born_solvation_energy_oddgrid()
    test_born_solvation_energy_pbc()
