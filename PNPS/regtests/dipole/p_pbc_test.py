import pnps
import math


def solve_p(PBC, Relaxation=1.0, MaxIterations=-2, DielConstBulk=80):
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
        DielConstBulk=DielConstBulk,
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
        [[True, False, False], 1.00, -2, 1666.7859575, 1777.6373415],
        [[False, True, False], 1.00, -2, 1666.7852251, 1777.6024294],
        [[False, False, True], 1.00, -2, 1666.7855913, 1777.6392336],
        [[True, True, False], 1.00, -2, 1681.4857872, 1741.9863391],
        [[True, False, True], 1.00, -2, 1666.2342607, 1772.9649780],
        [[False, True, True], 1.00, -2, 1663.5952836, 1782.8773930],
        [[True, True, True], 1.90, 600, 1666.7890093, 1777.8196535]
    ]
    E_List = []
    for PBC, Relaxation, MaxIterations, refEsol, refEvac in PBS_List:
        Esol = solve_p(PBC, Relaxation, MaxIterations, 80.0)
        Evac = solve_p(PBC, Relaxation, MaxIterations, 1.0)
        E_List.append(
            (PBC, Relaxation, MaxIterations, Esol, Evac, refEsol, refEvac)
        )

    print "%6s %6s %6s %14s %10s %14s %10s %10s %10s" % ("PBC_X", "PBC_Y", "PBC_Z",
                                                         "Esol", "dEsolRef", "Evac", "dEvacRef", "dE", "ddEref")
    for PBC, Relaxation, MaxIterations, Esol, Evac, refEsol, refEvac in E_List:
        dE = Esol - Evac
        dEref = refEsol - refEvac
        print "%6s %6s %6s %14.7f %10.7f %14.7f %10.7f %10.7f %10.7f" % (str(PBC[0]), str(PBC[1]), str(PBC[2]),
                                                                         Esol, Esol - refEsol,
                                                                         Evac, Evac - refEvac,
                                                                         dE, dE - dEref)
    ReferenceGeneration = False
    if ReferenceGeneration:
        # Reference generation
        for PBC, Relaxation, MaxIterations, Esol, Evac, refEsol, refEvac in E_List:
            print "[[%6s, %6s, %6s], %4.2f, %4d, %14.7f, %14.7f]," % (str(PBC[0]), str(PBC[1]), str(PBC[2]),
                                                                    Relaxation, MaxIterations,
                                                                    Esol, Evac)
    else:
        for PBC, Relaxation, MaxIterations, Esol, Evac, refEsol, refEvac in E_List:
            dE = Esol - Evac
            dEref = refEsol - refEvac

            # there are only 2 charges and floats has about 7 significant digits
            # so because the energy is on order of 1666.7890093 it is
            # reasonable to set tolerance to 0.005

            assert abs(Esol - refEsol) < 0.005
            assert abs(Evac - refEvac) < 0.005
            assert abs(dE - dEref) < 0.005


if __name__ == "__main__":
    test_on_dipole()
