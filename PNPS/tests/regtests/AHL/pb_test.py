
def test_ahl_poisson(shared_datadir):
    import pnps
    import os

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
        Cbulk=1.0,
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
        PQR=os.path.join(str(shared_datadir), "7ahlc_std_amber_pdb2pqr.pqr")
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

    pnps.SolvePB(contworld,
                 MaxIterationsLPB=-4,
                 MaxIterationsNPB=-4,
                 Convergence=0.0,
                 Relaxation=1.0)

    # contworld.WritePotential("S2_Pot.bin")
    E = contworld.SystemEnergy

    del contworld

    Eref = 1.63786767198194e+05

    def abs_rel_diff(val, ref):
        return abs((val - ref) / ref)

    dE = E - Eref
    dErel = abs_rel_diff(E, Eref)

    print("Difference with reference is %s" % dE)
    print("Absolute relative difference with reference is %s" % dErel)

    assert dErel < 3.0e-7


if __name__ == "__main__":
    # This is for stand-alone/debug execution
    import inspect
    import os
    cur_dir = os.path.dirname(os.path.abspath(inspect.getframeinfo(inspect.currentframe()).filename))
    test_ahl_poisson(os.path.join(cur_dir, "data"))
