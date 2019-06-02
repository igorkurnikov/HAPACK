
def abs_rel_diff(val, ref):
    return abs((val - ref) / ref)


def run_ahl_poisson(m_datadir, NumOfThreads=None):
    import pnps
    import os
    import time

    pnpsapp = pnps.PNPSApp.GetPNPSApp()
    if NumOfThreads is not None:
        pnpsapp.SetNumOfThreads(NumOfThreads)
    else:
        NumOfThreads = pnpsapp.GetNumOfThreads()

    print("Number of threads to use: %d" % pnpsapp.GetNumOfThreads())

    start = time.time()
    # Create Continuum representation
    contworld = pnps.ContWorld(
        GridSize=[171, 171, 171],
        GridScale=1.0,
        Qions=[1, -1]
    )

    # setup builder
    builder = pnps.GetWorldBuilder(
        DielConstBulk=80.0,
        DiffCoefBulk=[19.57, 20.32],
        Cbulk=0.0,
        RprobeDiel=1.4,
        RprobeDiff=0.5,
        Rions=[2.0, 2.0],
        DiffusionModel="Plane",
        BoundaryCondition="ZeroBC",
        MakeDielectricMap=True,
        MakeChargeMap=True,
        MakeDiffusionMap=False,
        MakeConcentrationMap=True,
        MakeSoftRepultionMap=False
    )
    # add atoms
    builder.addAtoms(
        DielConst=4.0,
        PQR=os.path.join(str(m_datadir), "7ahlc_std_amber_pdb2pqr.pqr")
    )
    # add membrane with hole, note that second value in R is large
    builder.addTube(
        DielConst=2.0,
        x=0.0, y=0.0,
        z=[-43.0, 15.0],
        R=[15.0, 250.5]
    )

    # Do building
    builder.BuildContWorld(contworld)

    # Write created maps to file system
    # contworld.WriteNodeIndexing("S1_System.systop")
    # Clean-up
    del builder

    startP = time.time()
    pnps.SolveP(contworld,
                MaxIterations=-2,
                Convergence=0.0,
                Relaxation=1.0)
    endP = time.time()
    # contworld.WritePotential("S2_Pot.bin")
    E = contworld.SystemEnergy

    del contworld

    end = time.time()
    return NumOfThreads, E, end - start, endP - startP


Eref = 1.63854188269094e+05


def test_ahl_poisson(shared_datadir):

    nt, E, t, tP = run_ahl_poisson(shared_datadir)

    dE = abs_rel_diff(E, Eref)

    print("Absolute relative difference with reference is %s" % dE)

    assert dE < 1.0e-7


def test_ahl_poisson_parallel(shared_datadir, NumOfThreads=None):
    if NumOfThreads is None:
        NumOfThreads=(1, 2, 4)

    results = []

    for nt in NumOfThreads:
        results.append(run_ahl_poisson(shared_datadir, nt))

    #Eref = results[0][1]
    print("%2s %20s %10s %10s %10s" % ("NT", "E", "dE", "Time", "Poisson Time"))
    for nt, E, t, tP in results:
        print("%2d %20.12e %10.3e %10.3f %10.3f" % (nt, E, abs(E - Eref), t, tP))
        dE = abs_rel_diff(E, Eref)
        assert dE < 1.0e-7


if __name__ == "__main__":
    # This is for stand-alone/debug execution
    import inspect
    import os
    import sys
    cur_dir = os.path.dirname(os.path.abspath(inspect.getframeinfo(inspect.currentframe()).filename))

    if len(sys.argv) == 1:
        test_ahl_poisson(os.path.join(cur_dir, "data"))

    if len(sys.argv) == 1:
        NumOfThreads = (1, 2, 4, 8)
    else:
        NumOfThreads = (int(nt) for nt in sys.argv[1:])

    test_ahl_poisson_parallel(os.path.join(cur_dir, "data"), NumOfThreads)
