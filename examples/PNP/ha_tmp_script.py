from time import *

#common parameters
mset=GetCurMolSet()

GridSize=65
GridScale=1.0
epsi=2.0;
epsout=80.0
epsout2=2.0
rionst=0.0
exrad=2.0
radprb=1.4

print "PNPMod"
StartPNP=time()

pnpmod=mset.GetPNPMod(1)

#Set Size of the Grid
intArray_setitem(pnpmod.InitWorldD.GridSize,0,GridSize)
intArray_setitem(pnpmod.InitWorldD.GridSize,1,GridSize)
intArray_setitem(pnpmod.InitWorldD.GridSize,2,GridSize)
#Set GridScale
pnpmod.InitWorldD.GridScale=GridScale
#Set Periodic Boundary Condition
boolArray_setitem(pnpmod.InitWorldD.PBC,0,0)
boolArray_setitem(pnpmod.InitWorldD.PBC,1,0)
boolArray_setitem(pnpmod.InitWorldD.PBC,2,0)
#where to put molecule PNPMod.DontMove/PNPMod.MoveToCenter/PNPMod.UseOffset
pnpmod.SetMoleculePosition(PNPMod.MoveToCenter)
#set boundary contition
pnpmod.BoundaryCond=1
#set dielectrical constant 0 is for Protein, 1 is for solvent
pnpmod.DielDiffTablesD.Epsilon[0]=epsi
pnpmod.DielDiffTablesD.Epsilon[1]=epsout
#set Ionic Strenrgh
pnpmod.BulkD.Concentration[0]=rionst
#Add Membrane  PNPMod.MemT_None/PNPMod.MemT_MembraneZ/PNPMod.MemT_Tube
#pnpmod.SetMembraneType(PNPMod.MemT_MembraneZ,4.0,-4.0,-2.5)
#pnpmod.SetMembraneType(PNPMod.MemT_Tube,4.0,-4.0,-2.5,1.0,0.0,2.0)
pnpmod.RunPoissonBolzmann()
EsolPNP=pnpmod.GeometryWorld.SystemEnergy

pnpmod.DielDiffTablesD.Epsilon[1]=epsout2
pnpmod.BulkD.Concentration[0]=0.0

pnpmod.RunPoissonBolzmann()
EvacPNP=pnpmod.GeometryWorld.SystemEnergy

dEPNP=EsolPNP-EvacPNP
EndPNP=time()
OperTimePNP=EndPNP-StartPNP

print "ElectostMod"
StartElMod=time()
elmod=mset.GetElectrostMod(1)
elmod.nx=GridSize
elmod.ny=GridSize
elmod.nz=GridSize
elmod.perfil=pnpmod.ConvertGridSizeGridScaleToPerFil(GridSize,GridScale)
elmod.rionst=rionst
elmod.exrad=exrad
elmod.radprb=radprb
elmod.epsi=epsi
elmod.epsout=epsout
elmod.boundary=3
elmod.nlit=-1

elmod.run(RUN_FOREGROUND)
EsolElMod=elmod.tot_ene

elmod.epsout=epsout2
elmod.rionst=0.0

elmod.run(RUN_FOREGROUND)
EvacElMod=elmod.tot_ene

dEElMod=EsolElMod-EvacElMod
EndElMod=time()
OperTimeElMod=EndElMod-StartElMod

print "APBS"
StartAPBS=time()
apbsmod=mset.GetAPBSMod(1)
apbsmod.nx=GridSize
apbsmod.ny=GridSize
apbsmod.nz=GridSize
apbsmod.GridScale=GridScale
apbsmod.rionst=rionst
apbsmod.exrad=exrad
apbsmod.radprb=radprb
apbsmod.epsi=epsi
apbsmod.epsout=epsout
apbsmod.chgm="spl0"
apbsmod.bcfl="mdh"
apbsmod.nlev=0
apbsmod.Run()

EsolAPBS=apbsmod.E

apbsmod.epsout=epsout2
apbsmod.rionst=0.0
apbsmod.Run()

EvacAPBS=apbsmod.E

dEAPBS=EsolAPBS-EvacAPBS
EndAPBS=time()
OperTimeAPBS=EndAPBS-StartAPBS

print "Solver       Esol \tEvac \tdE \tTime"
#print "Solver       Esol_w_Salt \tEsol \tdE"
print "PNPMOD       ",EsolPNP,"\t",EvacPNP,"\t",dEPNP,"\t",OperTimePNP
print "ElectrostMod ",EsolElMod,"\t",EvacElMod,"\t",dEElMod,"\t",OperTimeElMod
print "APBS         ",EsolAPBS,"\t",EvacAPBS,"\t",dEAPBS,"\t",OperTimeAPBS
