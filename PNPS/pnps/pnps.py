from pnpsll import pnpsll

from pnpsll.pnpsll import ContWorld, GetWorldBuilder, PNPSApp


def SolveP(contworld,
           MaxIterations=-1,
           Convergence=0.0,
           Relaxation=1.0,
           ConvergenceCheck=20,
           Solver=0,
           Verbose=True,
           QmobMod=2,
           MinIterations=0,
           ConvFacMaxHistory=1
           ):
    """
    Solve
    Input Parameters:
        MaxIterations=int, default=-1
            Maximal number of iterations, negative number means automaticly determined

    Returned value:
        Exit status
    """
    p = pnpsll.PoissonSolver()
    params = {
        'MaxIterations': MaxIterations,
        'Convergence': Convergence,
        'Relaxation': Relaxation,
        'ConvergenceCheck': ConvergenceCheck,
        'Solver': Solver,  # ,SolverStr
        'Verbose': Verbose,
        'QmobMod': QmobMod,
        'MinIterations': MinIterations,
        'ConvFacMaxHistory': ConvFacMaxHistory
    }
    p.LoadParamFromPyDict(params)

    p.SetContWorld(contworld)
    p.InitSolver()
    p.Solve()
    del p


def SolvePB(contworld,
            MaxIterationsLPB=-1,
            MaxIterationsNPB=-1,
            Convergence=0.0,
            Relaxation=1.0,
            ConvergenceCheck=20,
            Solver=0,
            Verbose=True
            ):
    """
    Solve
    Input Parameters:
        MaxIterationsLPB=int, default=-1
            Maximal number of iterations, negative number means automaticly determined

    Returned value:
        Exit status
    """
    pb = pnpsll.PoissonBoltzmannSolver()
    pb.MaxIterationsLPB = MaxIterationsLPB
    pb.MaxIterationsNPB = MaxIterationsNPB
    pb.Convergence = Convergence
    pb.Relaxation = Relaxation
    pb.ConvergenceCheck = ConvergenceCheck
    pb.solver = Solver
    if Verbose:
        pb.verbose = True
    else:
        pb.verbose = False

    pb.ShowParameters()
    pb.SetRelaxation(Relaxation)
    if pb.dielectricZSSUM is not None:
        pb.ShowProperties()

    pb.SetContWorld(contworld)
    pb.InitSolver()
    pb.Solve()
    del pb
