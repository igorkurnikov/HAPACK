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
    PQR=os.path.join("..", "7ahlc_std_amber_pdb2pqr.pqr")
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
contworld.WriteNodeIndexing("S1_System.systop")
# Clean-up
del builder

pnps.SolveP(contworld,
        MaxIterations=-2,
        Convergence=0.0,
        Relaxation=1.0)

contworld.WritePotential("S2_Pot.bin")

del contworld
