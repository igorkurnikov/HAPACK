import pnps
import math


def solve_p(PBC, Relaxation=1.0, MaxIterations=-2):
    """
    ZeroBC is set intentionally to boost PBC effects
    """
    # Create Continuum representation
    print PBC
    contworld = pnps.ContWorld(
        GridSize=[65, 65, 65],
        GridScale=2.0,
        PBC=PBC,
        Qions=[1, -1]
    )

    # setup builder
    builder = pnps.GetWorldBuilder(
        DielConstBulk=80,
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
        DielConst=2,
        AtomsPQR=[
            [1, 2, 0, 1, 2],
            [-0.5, -3, -1, -1, 3]
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
        [True, False, False, 1.0, -2, 1666.7860186],
        [False, True, False, 1.0, -2, 1666.7852862],
        [False, False, True, 1.0, -2, 1666.7855303 ],
        [True, True, False, 1.0, -2, 1681.4860314],
        [False, True, True, 1.0, -2, 1663.5953446],
        [True, True, True, 1.9, 600, 1666.7893145]
    ]

    E_List = [PBC+[solve_p(PBC[:3], PBC[3], PBC[4])] for PBC in PBS_List]

    print "%6s %6s %6s %20s %20s" % ("PBC_X", "PBC_Y", "PBC_Z", "E", "dE")
    for x, y, z, _, _, Eref, E in E_List:
        print "%6s %6s %6s %20.7f %20.7f" % (str(x), str(y), str(z), E, E - Eref)
    for x, y, z, _, _, Eref, E in E_List:
        assert abs(E - Eref) < 0.001


if __name__ == "__main__":
    test_on_dipole()
