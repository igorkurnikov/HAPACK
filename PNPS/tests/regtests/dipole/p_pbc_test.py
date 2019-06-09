import pnps
import math


def solve_p(PBC, Relaxation=1.0, MaxIterations=-2, atom1=None, atom2=None, DielBulk=2):
    """
    ZeroBC is set intentionally to boost PBC effects
    Also bulk dielectric constant default set to 2
    """
    if atom1 is None:
        atom1 = [1, 2, 0, 1, 2]
    if atom2 is None:
        atom2 = [-0.5, -3, -1, -1, 3]
    # Create Continuum representation
    print(PBC)
    contworld = pnps.ContWorld(
        GridSize=[64, 64, 64],
        GridScale=2.0,
        PBC=PBC,
        Qions=[1, -1]
    )

    # setup builder
    builder = pnps.GetWorldBuilder(
        DielConstBulk=1,
        Cbulk=0.0,
        DiffusionModel="Plane",
        BoundaryCondition="ZeroBC",
        MakeDielectricMap=True,
        MakeChargeMap=True,
        MakeDiffusionMap=False,
        MakeConcentrationMap=False,
        MakeSoftRepultionMap=False
    )
    # add atoms
    builder.addAtoms(
        DielConst=1,
        AtomsPQR=[
            atom1,
            atom2
        ]
    )
    # Do building
    builder.BuildContWorld(contworld)

    # Write created maps to file system
    # contworld.WriteNodeIndexing("S1_System_%s_%s_%s_systop" % tuple(str(pbc) for pbc in PBC))
    # Clean-up
    del builder

    pnps.SolveP(contworld,
                MaxIterations=MaxIterations,
                Convergence=0.0,
                Relaxation=Relaxation)

    # contworld.WritePotential("Pot_%s_%s_%s.bin" % tuple(str(pbc) for pbc in PBC))
    E = contworld.SystemEnergy

    del contworld
    return E


def test_on_dipole():
    pnpsapp = pnps.PNPSApp.GetPNPSApp()
    pnpsapp.SetNumOfThreads(1)

    PBS_List = [
        [True, False, False, 1.0, -2, 3456.4548547],
        [False, True, False, 1.0, -2, 3456.4073694],
        [False, False, True, 1.0, -2, 3456.4533899],
        [True, True, False, 1.0, -2, 3456.5276086],
        [False, True, True, 1.0, -2, 3456.5283411],
        [True, False, True, 1.0, -2, 3456.5964563],
        [True, True, True, 1.9, 600, 3456.6812952]
    ]

    E_List = [PBC+[solve_p(PBC[:3], PBC[3], PBC[4])] for PBC in PBS_List]

    print("%6s %6s %6s %20s %20s" % ("PBC_X", "PBC_Y", "PBC_Z", "E", "dE"))
    for x, y, z, _, _, Eref, E in E_List:
        print("%6s %6s %6s %20.7f %20.7f" % (str(x), str(y), str(z), E, E - Eref))
    for x, y, z, _, _, Eref, E in E_List:
        assert abs(E - Eref) < 0.002


def test_on_dipole_symmetry():
    """
    Set bulk diel constant to 1 same as atoms to avoid inconsistencies with probe rolling

    The test is done for two atoms along one of the axis (three cases x y and z)
    The tests are ordered is such way that energy should be same for i and i+6 and i+12
    within each group 1=2 and 4=5
    """
    pnpsapp = pnps.PNPSApp.GetPNPSApp()
    pnpsapp.SetNumOfThreads(1)

    PBS_List = [
        [True, False, False, -1.0, -2, [2, 0, 0, 1, 2], [-3, 0, 0, -1, 3]],
        [False, True, False, -1.0, -2, [2, 0, 0, 1, 2], [-3, 0, 0, -1, 3]],
        [False, False, True, -1.0, -2, [2, 0, 0, 1, 2], [-3, 0, 0, -1, 3]],
        [False, True, True, -1.0, -2, [2, 0, 0, 1, 2], [-3, 0, 0, -1, 3]],
        [True, True, False, -1.0, -2, [2, 0, 0, 1, 2], [-3, 0, 0, -1, 3]],
        [True, False, True, -1.0, -2, [2, 0, 0, 1, 2], [-3, 0, 0, -1, 3]],

        [False, True, False, -1.0, -2, [0, 2, 0, 1, 2], [0, -3, 0, -1, 3]],
        [True, False, False, -1.0, -2, [0, 2, 0, 1, 2], [0, -3, 0, -1, 3]],
        [False, False, True, -1.0, -2, [0, 2, 0, 1, 2], [0, -3, 0, -1, 3]],
        [True, False, True, -1.0, -2, [0, 2, 0, 1, 2], [0, -3, 0, -1, 3]],
        [True, True, False, -1.0, -2, [0, 2, 0, 1, 2], [0, -3, 0, -1, 3]],
        [False, True, True, -1.0, -2, [0, 2, 0, 1, 2], [0, -3, 0, -1, 3]],

        [False, False, True, -1.0, -2, [0, 0, 2, 1, 2], [0, 0, -3, -1, 3]],
        [True, False, False, -1.0, -2, [0, 0, 2, 1, 2], [0, 0, -3, -1, 3]],
        [False, True, False, -1.0, -2, [0, 0, 2, 1, 2], [0, 0, -3, -1, 3]],
        [True, True, False, -1.0, -2, [0, 0, 2, 1, 2], [0, 0, -3, -1, 3]],
        [False, True, True, -1.0, -2, [0, 0, 2, 1, 2], [0, 0, -3, -1, 3]],
        [True, False, True, -1.0, -2, [0, 0, 2, 1, 2], [0, 0, -3, -1, 3]],
    ]

    results = [[solve_p(PBC[:3], PBC[3], PBC[4], PBC[5], PBC[6], 1)] + PBC[:3] for PBC in PBS_List]

    print("%6s %6s %6s %20s %20s" % ("PBC_X", "PBC_Y", "PBC_Z", "E", "dE"))
    for E, x, y, z in results:
        print("%6s %6s %6s %20.3f %20.7f" % (str(x), str(y), str(z), E, 0))
    for i in range(6):
        assert abs(results[i][0] - results[i + 6][0]) < 0.002
        assert abs(results[i][0] - results[i + 12][0]) < 0.002
    for i in [0, 6, 12]:
        assert abs(results[i+1][0] - results[i + 2][0]) < 0.002
        assert abs(results[i+4][0] - results[i + 5][0]) < 0.002


if __name__ == "__main__":
    test_on_dipole()
    test_on_dipole_symmetry()
