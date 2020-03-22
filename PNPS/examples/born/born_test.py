import pnps
import os
import math

def theoretic_born_enery(diel_from, diel_to, q, R):
    BOHR_TO_ANG       = 0.529177249  # Bohrs to Angstroms
    HARTREE_TO_KT     = 1059.642  # HARTREE TO kT at 298K
    
    de=-0.5*HARTREE_TO_KT*(1/diel_from-1/diel_to)*q*q/(R/BOHR_TO_ANG)
    return de


def solve_p(diel_ion, diel_sol, x=0.0, y=0.0, z=0.0, q=1.0, R=2.0):
    # Create Continuum representation
    contworld = pnps.ContWorld(
        GridSize=[65, 65, 65],
        GridScale=4.0,
        Qions=[1, -1]
    )
    
    # setup builder
    builder = pnps.GetWorldBuilder(
        DielConstBulk=diel_sol,
        Cbulk=0.0,
        DiffusionModel="Plane",
        BoundaryCondition="CoulBC",
        MakeDielectricMap=True,
        MakeChargeMap=True,
        MakeDiffusionMap=False,
        MakeConcentrationMap=False,
        MakeSoftRepultionMap=False
    )
    # add atoms
    builder.addAtoms(
        DielConst=diel_ion,
        AtomsPQR=[[x,y,z,q,R]]
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
    E=contworld.SystemEnergy
    
    del contworld
    return E

def test_born_solvation_energy():
    diel_ion=2.0
    diel_sol=80.0
    x=0.0
    y=0.0
    z=0.0
    q=1.0
    R=2.0
    grid_scale=4  # 2 grid/A
    n=3

    #dEsim = solve_p(diel_ion,diel_sol,x, y, z, q, R)-solve_p(diel_ion,diel_ion,x, y, z, q, R)

    dEsim=0.0
    dEsimSQ=0.0
    h=1.0/grid_scale
    tot_num=0
    for ix in range(n):
        for iy in range(ix+1):
            for iz in range(iy+1):
            x=ix*h/(n+1.0)
            y=iy*h/(n+1.0)
            z=iz*h/(n+1.0)
            dE = solve_p(diel_ion,diel_sol,x, y, z, q, R)-solve_p(diel_ion,diel_ion,x, y, z, q, R)
            dEsim += dE
            dEsimSQ += dE*dE
            tot_num+=1

    dEsim = dEsim/tot_num
    dEsimSQ = dEsimSQ/tot_num
    sd = dEsimSQ-dEsim*dEsim
    sd = math.sqrt(abs(sd))
    dEtheor = theoretic_born_enery(diel_ion, diel_sol, q, R)

    print("Simulated solvation energy:       %.3f +- %.3f kT" % (dEsim,sd))
    print("Theoretical solvation energy:     %.3f" % dEtheor)
    abs_diff=dEsim-dEtheor
    print("Absolute difference (sim-theory): %.3f" % abs_diff)
    rel_diff=abs_diff/dEtheor
    print("Relative difference (sim-theory): %.3f" % rel_diff)

    ref_dEsim = -69.12786293519918
    ref_sd = 0.1477743546953387
    ref_abs_diff = -0.787865744857
    ref_rel_diff = 0.0115286183384

    def abs_rel_diff(ref,val):
        return abs((val - ref)/ref)

    assert abs_rel_diff(ref_dEsim, dEsim) < 1.0e-4
    assert abs_rel_diff(ref_sd, sd) < 1.0e-4
    assert abs_rel_diff(ref_abs_diff, abs_diff) < 1.0e-4

if __name__ == "__main__":
    test_born_solvation_energy()
