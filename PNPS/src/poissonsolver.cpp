//
// C++ Implementation: poissonsolver
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#include "poissonsolver.h"
#include "tinyxml.h"
#include "pnpdebug.h"
#include "contworld.h"
#include "math.h"
#include "pnpconstants.h"
#include "mapio.h"

#include "pnputil.h"
#define FMUL(X,Y) (X*Y)
#define FADD(X,Y) (X+Y)
#define FMAF(X,Y,Z) (X*Y+Z)

#include <stdlib.h>

#include <assert.h>

#include <algorithm>


PoissonSolver::PoissonSolver()
 : GenericSolver()
{
	InitZero();
}


PoissonSolver::~PoissonSolver()
{
	Clear();
}
int PoissonSolver::InitZero()
{
	PnpsObject::SetName("PoissonSolver");
	SolverStr.push_back("Auto");
	SolverStr.push_back("NodeIndexBased");
	SolverStr.push_back("ArrayDirect");
	SolverStr.push_back("PNPC");
	
	World=NULL;
	
	ConvergenceCheck=20;
	QmobMod=0;
	MaxIterations=300;
	Convergence=0.0;
	Relaxation=1.9;
	solver=0;
	verbose=true;
	QmobMod=2;
	MinIterations=0;

	 //Initial values of variables:
	totalChange=0;
	relativeChange=0;
	totalEnergy=0;
	totalEnergyInd=0;
	verbose=true;
	
	dielectricZSSUM=NULL;
	dielectricZDBSUM=NULL;
	
	dielectricXS = NULL;
	dielectricYS = NULL;
	dielectricZS = NULL;
	dielectricZSSUM = NULL;
	dielectricXmS = NULL;
	dielectricYmS = NULL;
	dielectricZmS = NULL;
	QstS = NULL;
	dielectricXDB = NULL;
	dielectricYDB = NULL;
	dielectricZDB = NULL;
	dielectricZDBSUM = NULL;
	dielectricXmDB = NULL;
	dielectricYmDB = NULL;
	dielectricZmDB = NULL;
	dielectricCh = NULL;
	
	IndexNoSingular=NULL;
	IndexDielBoarder=NULL;
	IndexCharge=NULL;
	IndexSingular=NULL;
	//ChargeSum=NULL;
	Qst=NULL;
	PhiCharge=NULL;
	PhiSingular=NULL;
	QmobNum[0]=0;
	QmobNum[1]=0;
	QmobNum[2]=0;
	IndexQmob=NULL;
	dielectricChMob=NULL;
	Qmob=NULL;
	
	QmobDielBoarderNum[0]=0;
	QmobDielBoarderNum[1]=0;
	QmobDielBoarderNum[2]=0;
	IndexQmobDielBoarder=NULL;
	QmobDielBoarder=NULL;
	dielectricXQmobDB=NULL;
	dielectricYQmobDB=NULL;
	dielectricZQmobDB=NULL;
	dielectricXmQmobDB=NULL;
	dielectricYmQmobDB=NULL;
	dielectricZmQmobDB=NULL;
	dielectricZQmobDBSUM=NULL;
	
	QmobDielBoarderQstNum[0]=0;
	QmobDielBoarderQstNum[1]=0;
	QmobDielBoarderQstNum[2]=0;
	IndexQmobDielBoarderQst=NULL;
	QmobDielBoarderQst=NULL;
	QstQmobDielBoarderQst=NULL;
	dielectricXQmobDBQst=NULL;
	dielectricYQmobDBQst=NULL;
	dielectricZQmobDBQst=NULL;
	dielectricXmQmobDBQst=NULL;
	dielectricYmQmobDBQst=NULL;
	dielectricZmQmobDBQst=NULL;
	dielectricZQmobDBSUMQst=NULL;
	
	//QmobFlag=NULL;
	
	ChargeSum=NULL;
	
	ChargesEnergy=NULL;
	potential=NULL;
	
	WayToCalcSystemEnergy=0;
	ConvFacMaxHistory=1;

	CalcVolume=NULL;
    PotBW=nullptr;
	return EXIT_SUCCESS;
}
int PoissonSolver::Clear()
{
//	DeleteCArray(IonsQ);
//	DeleteCVecArray(C,NIonsTypes);
//	DeleteObjByPnt(NIndexing);
	DeleteCArray(IndexQmob);
	DeleteCArray(Qmob);
	//DeleteCArray(QmobFlag);
	DeleteCArray(dielectricChMob);
	DeleteCArray(dielectricXS);
	DeleteCArray(dielectricYS);
	DeleteCArray(dielectricZS);
	DeleteCArray(dielectricZSSUM);
	DeleteCArray(dielectricXmS);
	DeleteCArray(dielectricYmS);
	DeleteCArray(dielectricZmS);
	DeleteCArray(QstS);
	DeleteCArray(dielectricXDB);
	DeleteCArray(dielectricYDB);
	DeleteCArray(dielectricZDB);
	DeleteCArray(dielectricZDBSUM);
	DeleteCArray(dielectricXmDB);
	DeleteCArray(dielectricYmDB);
	DeleteCArray(dielectricZmDB);
	DeleteCArray(dielectricCh);
	DeleteCArray(IndexNoSingular);
	DeleteCArray(IndexDielBoarder);
	DeleteCArray(IndexCharge);
	DeleteCArray(IndexSingular);
	DeleteCArray(Qst);
	
	DeleteCArray(IndexQmobDielBoarder);
	DeleteCArray(QmobDielBoarder);
	DeleteCArray(dielectricXQmobDB);
	DeleteCArray(dielectricYQmobDB);
	DeleteCArray(dielectricZQmobDB);
	DeleteCArray(dielectricXmQmobDB);
	DeleteCArray(dielectricYmQmobDB);
	DeleteCArray(dielectricZmQmobDB);
	DeleteCArray(dielectricZQmobDBSUM);
	
	DeleteCArray(IndexQmobDielBoarderQst);
	DeleteCArray(QmobDielBoarderQst);
	DeleteCArray(QstQmobDielBoarderQst);
	DeleteCArray(dielectricXQmobDBQst);
	DeleteCArray(dielectricYQmobDBQst);
	DeleteCArray(dielectricZQmobDBQst);
	DeleteCArray(dielectricXmQmobDBQst);
	DeleteCArray(dielectricYmQmobDBQst);
	DeleteCArray(dielectricZmQmobDBQst);
	DeleteCArray(dielectricZQmobDBSUMQst);
	
	DeleteCArray(ChargesEnergy);
	DeleteCArray(PhiCharge);
	DeleteCArray(PhiSingular);
	
	DeleteCArray(CalcVolume);
#	ifdef PNPDOUBLE
	DeleteCArray(potential);
#	endif
    DeleteObjByPnt(PotBW);
	return EXIT_SUCCESS;
}
int PoissonSolver::LoadParamFromPyDict(PyObject *dict)
{
	Clear();
	int i,gridPoint;
	//Read Primary Parameters
	haPyDict_GetItemAsInt(dict,"MaxIterations",&MaxIterations);
	haPyDict_GetItemAsFloat(dict,"Relaxation",&Relaxation);
	haPyDict_GetItemAsFloat(dict,"Convergence",&Convergence);
	haPyDict_GetItemAsFloat(dict,"Tolerance",&Convergence);
	haPyDict_GetItemAsInt(dict,"ConvergenceCheck",&ConvergenceCheck);
	haPyDict_GetItemAsInt(dict,"Solver",&solver);
	haPyDict_GetItemAsBool(dict,"Verbose",&verbose);

	haPyDict_GetItemAsInt(dict,"QmobMod",&QmobMod);
	haPyDict_GetItemAsInt(dict,"MinIterations",&MinIterations);
	haPyDict_GetItemAsInt(dict,"ConvFacMaxHistory",&ConvFacMaxHistory);

	ShowParameters();
	
	SetRelaxation(Relaxation);
	
	if(dielectricZSSUM!=NULL)ShowProperties();
	
	return EXIT_SUCCESS;
}
int PoissonSolver::SetRelaxation(float _Relaxation)
{
	int i,gridPoint;
	Relaxation=_Relaxation;
	//Scale Parameters to Internal Units
	//Calculate Derivative Parameters
	if(dielectricZSSUM!=NULL)for(i=0;i<SingularNum[2];i++){
		gridPoint=IndexSingular[i];
		dielectricZSSUM[i]=dielectricXS[i]+dielectricYS[i]+dielectricZS[i]+dielectricXmS[i]+dielectricYmS[i]+dielectricZmS[i];
		dielectricZSSUM[i]=Relaxation/dielectricZSSUM[i];
	}
	if(dielectricZDBSUM!=NULL)for(i=0;i<DielBoarderNum[2];i++){
		gridPoint=IndexDielBoarder[i];
		dielectricZDBSUM[i]=dielectricXDB[i]+dielectricYDB[i]+dielectricZDB[i]+dielectricXmDB[i]+dielectricYmDB[i]+dielectricZmDB[i];
		dielectricZDBSUM[i]=Relaxation/dielectricZDBSUM[i];
	}
	if(dielectricZQmobDBSUM!=NULL)for(i=0;i<QmobDielBoarderNum[2];i++){
		gridPoint=IndexQmobDielBoarder[i];
		dielectricZQmobDBSUM[i]=dielectricXQmobDB[i]+dielectricYQmobDB[i]+dielectricZQmobDB[i]+dielectricXmQmobDB[i]+dielectricYmQmobDB[i]+dielectricZmQmobDB[i];
		dielectricZQmobDBSUM[i]=Relaxation/dielectricZQmobDBSUM[i];
	}
	om2 = Relaxation;
	om1 = 1.0-om2;
	om2d6 = om2/6.0;
	
	if(dielectricZSSUM!=NULL)ShowProperties();
	return EXIT_SUCCESS;
}
int PoissonSolver::SetContWorld(ContWorld* _world)
{
	int i;
	World=_world;
	GridScale = World->GridScale;
	GS_X = World->GridSize[0];
	GS_Y = World->GridSize[1];
	GS_Z = World->GridSize[2];
	GS_XY = GS_X*GS_Y;
	GS_XYZ = GS_XY*GS_Z;
	if(World->Potential == NULL){
		if(!(World->Potential = new float[GS_XYZ])){
			fprintf(stderr,"ERROR 104: No memory available\n");
			exit(104);
		}
		for(i=0;i<GS_XYZ;i++)World->Potential[i]=0.0;
	}
#	ifndef PNPDOUBLE
	potential = World->Potential;
#	else
	
#	endif

    if (PotBW == nullptr) {
        PotBW = new FieldBW();
    }
    if (!PotBW->SameSize(World->GridSize)) {
        PotBW->Init(World->GridSize);
    }

	return EXIT_SUCCESS;
}
int PoissonSolver::ShowParameters()
{
	DbgPrint2("CPoisson::ShowParameters\n");
	
	pnpPrintGroup0("\nParameters of Poisson solver set up\n");
	pnpPrintGroup0("    MaxIterations:.................... %d\n", MaxIterations);
	pnpPrintGroup0("    Convergence:...................... %.8g kT\n", Convergence);
	pnpPrintGroup0("    Relaxation:....................... %.5g\n", Relaxation);
	
	return EXIT_SUCCESS;
}
int PoissonSolver::ShowProperties()
{
	DbgPrint2("CPoisson::ShowProperties\n");
	return EXIT_SUCCESS;
}
int PoissonSolver::InitSolver()
{
  //Solver{Auto=0,NodeIndexBased=1,ArrayDirect=2,PNPC=3};
  //Solver solver;
	int status;
	bool bGuessNumberOfIteration=false;
	if(Relaxation<0.0||MaxIterations<0)
	{
		bGuessNumberOfIteration=true;
		Relaxation=1.0;
	}
	if(solver==Auto)
	{
		if(World->NIndexing==NULL)status=InitSolverAD();
		else status=InitSolverNIB();
	}
	else if(solver==NodeIndexBased)
		status=InitSolverNIB();
	else if(solver==ArrayDirect)
		status=InitSolverAD();
	else if(solver==PNPC)
		status=InitSolverW();
	
	int itmp;
	if(bGuessNumberOfIteration)
		itmp=GuessNumberOfIteration();
	if(MaxIterations<0)
	{
		MaxIterations=itmp*abs(MaxIterations);
		if(MaxIterations<60)
			MaxIterations=60;
	}
	return status;
}
int PoissonSolver::InitSolverW()
{
	DbgPrint0("CPoisson::InitSolverW()");
	long count = 0;
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	long gridSizeXY;
	long gridPoint;
	int i,j,k;
	float * dielectric[3];
	float * staticCharge;
	long jgrid,kgrid;

	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	
	poissonSolverData = new PoissonSolverDataW();
	assert(poissonSolverData!=NULL);
	
	poissonSolverData->maxIterations = MaxIterations;
	poissonSolverData->convergence = Convergence;
	poissonSolverData->relaxation = Relaxation;
	
	World->CheckArrays("P",true);
	
	if(World->C!=NULL)
		poissonSolverData->hasDynamicCharges = 1;
	else
		poissonSolverData->hasDynamicCharges = 0;
	
	
	for(i=0;i<3;i++)
	{
		dielectric[i] = World->Epsilon[i];
	}
	staticCharge = World->Qstat;

	for(k=1;k<gridSizeZ-1;k++) {
		kgrid = k*gridSizeXY;
		for(j=1;j<gridSizeY-1;j++) {
			jgrid = kgrid+j*gridSizeX;
			for(i=1;i<gridSizeX-1;i++) {
				gridPoint = jgrid+i;
				if(dielectric[0][gridPoint]!=dielectric[0][gridPoint-1] || dielectric[0][gridPoint]!=dielectric[1][gridPoint] || dielectric[0][gridPoint]!=dielectric[1][gridPoint-gridSizeX] || dielectric[0][gridPoint]!=dielectric[2][gridPoint] || dielectric[0][gridPoint]!=dielectric[2][gridPoint-gridSizeXY] || (poissonSolverData->hasDynamicCharges==0 && staticCharge[gridPoint]!=0))
					count++;
			}
		}
	}

	poissonSolverData->borderPoints = new long[count+1];
	assert(poissonSolverData->borderPoints!=NULL);
	poissonSolverData->om2InverseDielectricSum = new float [count];
	assert(poissonSolverData->om2InverseDielectricSum!=NULL);
	if(!poissonSolverData->hasDynamicCharges) {
		poissonSolverData->typeOfBorderPoint = new short[count];
		assert(poissonSolverData->typeOfBorderPoint!=NULL);
	}
	else
		poissonSolverData->typeOfBorderPoint = NULL;

	count = 0;

	for(k=1;k<gridSizeZ-1;k++) {
		kgrid = k*gridSizeXY;
		for(j=1;j<gridSizeY-1;j++) {
			jgrid = kgrid+j*gridSizeX;
			for(i=1;i<gridSizeX-1;i++) {
				gridPoint = jgrid+i;
				if(dielectric[0][gridPoint]!=dielectric[0][gridPoint-1] || dielectric[0][gridPoint]!=dielectric[1][gridPoint] || dielectric[0][gridPoint]!=dielectric[1][gridPoint-gridSizeX] || dielectric[0][gridPoint]!=dielectric[2][gridPoint] || dielectric[0][gridPoint]!=dielectric[2][gridPoint-gridSizeXY]) {
					if(poissonSolverData->hasDynamicCharges==0) {
						if(staticCharge[gridPoint]!=0)
							poissonSolverData->typeOfBorderPoint[count] = POISSON_SOLVER_CHARGED_UNCOMMON_DIELECTRIC_POINT;
						else 
							poissonSolverData->typeOfBorderPoint[count] = POISSON_SOLVER_UNCOMMON_DIELECTRIC_POINT;
					}
					poissonSolverData->om2InverseDielectricSum[count] = Relaxation/(dielectric[0][gridPoint]+dielectric[0][gridPoint-1]+dielectric[1][gridPoint]+dielectric[1][gridPoint-gridSizeX]+dielectric[2][gridPoint]+dielectric[2][gridPoint-gridSizeXY]);
					poissonSolverData->borderPoints[count] = gridPoint;
					count++;
				}
				else if(poissonSolverData->hasDynamicCharges==0 && staticCharge[gridPoint]!=0) {
					poissonSolverData->typeOfBorderPoint[count] = POISSON_SOLVER_CHARGED_POINT;
					poissonSolverData->borderPoints[count] = gridPoint;
					count++;
				}
			}
		}
	}
	poissonSolverData->borderPoints[count] = gridSizeXY*gridSizeZ; 
	return EXIT_SUCCESS;
}
int PoissonSolver::InitSolverNIB()
{
	DbgPrint2("CPoisson::InitSolverNIB\n");
	if(World==NULL)return NULL;
	
	//if(World->C!=NULL)SetQmobForPNP();
	//temp vars
	int IType;
	int i,j,k,kgrid,jgrid,BlackOrWhite;
	int GrdPnt;
	
	NodeIndex* NIndex=World->NIndexing->NIndex;
	//count b/w regions
	DielBoarderNum[0]=0;DielBoarderNum[1]=0;DielBoarderNum[2]=0;
	ChargeNum[0]=0;ChargeNum[1]=0;ChargeNum[2]=0;
	
	//correction for spesific type of calculation, P(no additional charge) P(NP)(dynamic charge)
	unsigned int specChargeMask=NodeIndexing::ChargeMask;
	//!@todo not forget to change when will do pnp
//   if(World->D!=NULL)//i.e. P(NP)
//   {
//     specChargeMask=ChargeMask|DiffMask;
//   }
	unsigned int ChargeDielBoarderMask=specChargeMask|NodeIndexing::DielBoarderMask;
	unsigned int BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;
	unsigned int DielBoarderMask=NodeIndexing::DielBoarderMask;
	DbgPrint0("DielBoarderMask=%X\n",DielBoarderMask);
	DbgPrint0("ChargeDielBoarderMask=%X\n",ChargeDielBoarderMask);
	DbgPrint0("BlackAndWhiteMask=%X\n",BlackAndWhiteMask);
	DbgPrint0("specChargeMask=%X\n",specChargeMask);
	
	int QmobTot=QmobNum[2]+QmobDielBoarderNum[2]+QmobDielBoarderQstNum[2];
	bool bQmobHere;
	bool bCalcVolume;
	assert(CalcVolume == nullptr);// currently not working
	bCalcVolume=true;
	
	for(k=1;k<GS_Z-1;k++)
	{
		kgrid = k*GS_XY;
		for(j=1;j<GS_Y-1;j++)
		{
			jgrid = kgrid+j*GS_X;
			for(i=1;i<GS_X-1;i++)
			{
				GrdPnt = jgrid+i;
				if(CalcVolume!= nullptr)
					bCalcVolume=CalcVolume[GrdPnt];
				if(bCalcVolume)
				{
					if(NIndex[GrdPnt]&specChargeMask)
					{
						ChargeNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
						ChargeNum[2]++;
					}
					if(NIndex[GrdPnt]&DielBoarderMask)
					{
						DielBoarderNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
						DielBoarderNum[2]++;
					}
				}
			}
		}
	}
	DbgPrint1("InitSolver:Total: MyrankWorld->MyRank=%d icharge=% d    \niDielBoarder=%d  QmobNum=%d\n",
						World->MyRank,ChargeNum[2],DielBoarderNum[2],QmobNum[2]);
	DbgPrint0("TotalNodes: %d\n",ChargeNum[2]+DielBoarderNum[2]+SingularNum[2]+NoSingularNum[2]+QmobNum[2]);
	DbgPrint1("InitSolver:Blach cell:: MyrankWorld->MyRank=%d icharge=% d    \niDielBoarder=%d QmobNum=%d\n",
						World->MyRank,ChargeNum[1],DielBoarderNum[1],QmobNum[2]);
  

	if(!(IndexDielBoarder = new int[DielBoarderNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(IndexCharge = new int[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(Qst = new float[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(dielectricCh = new float[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(PhiCharge = new float[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}

	//Allocate dielectric helping array
	dielectricXDB = new float[DielBoarderNum[2]];
	dielectricYDB = new float[DielBoarderNum[2]];
	dielectricZDB = new float[DielBoarderNum[2]];
	dielectricXmDB = new float[DielBoarderNum[2]];
	dielectricYmDB = new float[DielBoarderNum[2]];
	dielectricZmDB = new float[DielBoarderNum[2]];

	if (!dielectricXDB || !dielectricYDB || !dielectricZDB || !dielectricXmDB || !dielectricYmDB || !dielectricZmDB) {
		fprintf(stderr, "ERROR 204: No memory available\n");
		exit(204);
	}

	for(i=0;i<ChargeNum[2];i++)
		PhiCharge[i]=0.0f;
	
    //fill indexes of nodes types
	int iCharge1=0,iDielBoarder1=0;
	int iCharge2=ChargeNum[1],iDielBoarder2=DielBoarderNum[1];
	int QCount=0;
	float *Q=World->NIndexing->Q;
	float *Eps=World->NIndexing->Eps;
	double q=0.0f;
	for(i=0;i<World->NIndexing->QNum;i++)
	{
		q+=Q[i];
	}
	DbgPrint0("q=%f QCount=%d GridScale=%f\n",(float)q/4.0f/M_PI/World->GridScale,World->NIndexing->QNum,World->GridScale);
	q=0.0f;

    int bwStrideX = PotBW->StrideX;
    int bwStrideXY = PotBW->StrideXY;

	for(k=1;k<GS_Z-1;k++)
	{
		kgrid = k*GS_XY;
		for(j=1;j<GS_Y-1;j++)
		{
			jgrid = kgrid+j*GS_X;
			for(i=1;i<GS_X-1;i++)
			{
				GrdPnt = jgrid+i;
				if(CalcVolume!=NULL)
					bCalcVolume=CalcVolume[GrdPnt];
				if(bCalcVolume)
				{
					if(NIndex[GrdPnt]&specChargeMask)
					{
                        int iu = i/2 + j * bwStrideX + k * bwStrideXY;
						int iCharge;
						if (NIndex[GrdPnt] & BlackAndWhiteMask)
						{
							iCharge = iCharge1;
							iCharge1++;
						}
						else
						{
							iCharge = iCharge2;
							iCharge2++;
						}
						IndexCharge[iCharge] = iu;
						if (NIndex[GrdPnt] & DielBoarderMask){
							dielectricCh[iCharge] = \
								Eps[(NIndex[GrdPnt] & NodeIndexing::Epsilon0) >> NodeIndexing::Epsilon0Sft] + \
								Eps[(NIndex[GrdPnt] & NodeIndexing::Epsilon1) >> NodeIndexing::Epsilon1Sft] + \
								Eps[(NIndex[GrdPnt] & NodeIndexing::Epsilon2) >> NodeIndexing::Epsilon2Sft] + \
								Eps[(NIndex[GrdPnt - 1] & NodeIndexing::Epsilon0) >> NodeIndexing::Epsilon0Sft] + \
								Eps[(NIndex[GrdPnt - GS_X] & NodeIndexing::Epsilon1) >> NodeIndexing::Epsilon1Sft] + \
								Eps[(NIndex[GrdPnt - GS_XY] & NodeIndexing::Epsilon2) >> NodeIndexing::Epsilon2Sft];
							dielectricCh[iCharge] /= 6.0f;
						}
						else {
							dielectricCh[iCharge] = Eps[(NIndex[GrdPnt] & NodeIndexing::Epsilon0) >> NodeIndexing::Epsilon0Sft];
						}

						Qst[iCharge] = Q[QCount] / dielectricCh[iCharge];
						q += Q[QCount];

					}
					if(NIndex[GrdPnt]&DielBoarderMask)
					{
						int iDielBoarder;
                        int iu = i / 2 + j * bwStrideX + k * bwStrideXY;
						if(NIndex[GrdPnt]&BlackAndWhiteMask){
							iDielBoarder = iDielBoarder1;
							iDielBoarder1++;
						}
						else{
							iDielBoarder = iDielBoarder2;
							iDielBoarder2++;
						}
						IndexDielBoarder[iDielBoarder]=iu;

						dielectricXDB[iDielBoarder] = Eps[(NIndex[GrdPnt] & NodeIndexing::Epsilon0) >> NodeIndexing::Epsilon0Sft];
						dielectricYDB[iDielBoarder] = Eps[(NIndex[GrdPnt] & NodeIndexing::Epsilon1) >> NodeIndexing::Epsilon1Sft];
						dielectricZDB[iDielBoarder] = Eps[(NIndex[GrdPnt] & NodeIndexing::Epsilon2) >> NodeIndexing::Epsilon2Sft];
						dielectricXmDB[iDielBoarder] = Eps[(NIndex[GrdPnt - 1] & NodeIndexing::Epsilon0) >> NodeIndexing::Epsilon0Sft];
						dielectricYmDB[iDielBoarder] = Eps[(NIndex[GrdPnt - GS_X] & NodeIndexing::Epsilon1) >> NodeIndexing::Epsilon1Sft];
						dielectricZmDB[iDielBoarder] = Eps[(NIndex[GrdPnt - GS_XY] & NodeIndexing::Epsilon2) >> NodeIndexing::Epsilon2Sft];

						double diel_sum = dielectricXDB[iDielBoarder] + dielectricYDB[iDielBoarder] + dielectricZDB[iDielBoarder] + \
							dielectricXmDB[iDielBoarder] + dielectricYmDB[iDielBoarder] + dielectricZmDB[iDielBoarder];

						diel_sum /= 6.0;
						
						dielectricXDB[iDielBoarder] = dielectricXDB[iDielBoarder] / diel_sum - 1.0f;
						dielectricYDB[iDielBoarder] = dielectricYDB[iDielBoarder] / diel_sum - 1.0f;
						dielectricZDB[iDielBoarder] = dielectricZDB[iDielBoarder] / diel_sum - 1.0f;
						dielectricXmDB[iDielBoarder] = dielectricXmDB[iDielBoarder] / diel_sum - 1.0f;
						dielectricYmDB[iDielBoarder] = dielectricYmDB[iDielBoarder] / diel_sum - 1.0f;
						dielectricZmDB[iDielBoarder] = dielectricZmDB[iDielBoarder] / diel_sum - 1.0f;

					}
						
				}
				if(NIndex[GrdPnt]&specChargeMask)
					QCount++;
			}
		}
	}
	DbgPrint0("q=%f QCount=%d",(float)q/4/M_PI/World->GridScale,QCount);
	return EXIT_SUCCESS;
}
int PoissonSolver::InitSolverAD()
{
      /// @todo implement me
  //std::istrstream iStr(ProcedureCommand);
  //char tmp[MAP_IO_STRING_LENGTH];
  
  //Read parameters
	DbgPrint2("Renew Singularities List\n");

	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	int gridSizeXY,gridSizeXYZ;
	int gridPoint;
	int i,j,k;
	float * dielectric[3];
	float * staticCharge;
	int jgrid,kgrid;
	int iCharge,iDielBoarder,iSingular,iNoSingular;
	int BlackOrWhite;
	float *positiveCharge;
	float *negativeCharge;
  
	int *typeOfBorderPoint=NULL;
	
	World->CheckArrays("P",true);
	
	if(World->C==NULL)
	{
		positiveCharge = NULL;
		negativeCharge = NULL;
	}
	else
	{
		positiveCharge = World->C[0];
		negativeCharge = World->C[1];
	}
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeZ*gridSizeXY;

	if(dielectricXS!=NULL){delete [] dielectricXS;dielectricXS=NULL;}
	if(dielectricYS!=NULL){delete [] dielectricYS;dielectricYS=NULL;}
	if(dielectricZS!=NULL){delete [] dielectricZS;dielectricZS=NULL;}
	if(dielectricZSSUM!=NULL){delete [] dielectricZSSUM;dielectricZSSUM=NULL;}
	if(dielectricXmS!=NULL){delete [] dielectricXmS;dielectricXmS=NULL;}
	if(dielectricYmS!=NULL){delete [] dielectricYmS;dielectricYmS=NULL;}
	if(dielectricZmS!=NULL){delete [] dielectricZmS;dielectricZmS=NULL;}
	if(dielectricXDB!=NULL){delete [] dielectricXDB;dielectricXDB=NULL;}
	if(dielectricYDB!=NULL){delete [] dielectricYDB;dielectricYDB=NULL;}
	if(dielectricZDB!=NULL){delete [] dielectricZDB;dielectricZDB=NULL;}
	if(dielectricZDBSUM!=NULL){delete [] dielectricZDBSUM;dielectricZDBSUM=NULL;}
	if(dielectricXmDB!=NULL){delete [] dielectricXmDB;dielectricXmDB=NULL;}
	if(dielectricYmDB!=NULL){delete [] dielectricYmDB;dielectricYmDB=NULL;}
	if(dielectricZmDB!=NULL){delete [] dielectricZmDB;dielectricZmDB=NULL;}
	if(dielectricCh!=NULL){delete [] dielectricCh;dielectricCh=NULL;}
	if(IndexNoSingular!=NULL){delete [] IndexNoSingular;IndexNoSingular=NULL;}
	if(IndexDielBoarder!=NULL){delete [] IndexDielBoarder;IndexDielBoarder=NULL;}
	if(IndexCharge!=NULL){delete [] IndexCharge;IndexCharge=NULL;}
	if(IndexSingular!=NULL){delete [] IndexSingular;IndexSingular=NULL;}
	if(ChargeSum!=NULL){delete [] ChargeSum;ChargeSum=NULL;}
  

  //Allocate Potential
	if(World->Potential == NULL){
		if(!(World->Potential = new float[gridSizeXYZ])){
			fprintf(stderr,"ERROR 104: No memory available\n");
			exit(104);
		}
		for(j=0;j<GS_XYZ;j++)World->Potential[j]=0.0;
	}
  

	for(i=0;i<3;i++)
		dielectric[i] = World->Epsilon[i];
  
	if(!(typeOfBorderPoint = new int[GS_XYZ])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}


	if(!(ChargeSum = new float[GS_XYZ])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	} 
	staticCharge=ChargeSum;
  
	if(QmobMod==0||World->D==NULL)
	{
		for(i=0;i<GS_XYZ;i++)
			if(World->Qstat[i]!=0)staticCharge[i]=1.0;
	}
	else {
		for(i=0;i<GS_XYZ;i++)
			if(World->Qstat[i]!=0||World->D[0][i]>0||World->D[1][i]>0)staticCharge[i]=1.0;
	}
	for(k=0;k<GS_XYZ;k++){
    //dielectricSum[k] = 6*dielectric[0][k];
		typeOfBorderPoint[k] = Boarder;
	}

  
	NoSingularNum[0]=0;NoSingularNum[1]=0;NoSingularNum[2]=0;
	SingularNum[0]=0;SingularNum[1]=0;SingularNum[2]=0;
	DielBoarderNum[0]=0;DielBoarderNum[1]=0;DielBoarderNum[2]=0;
	ChargeNum[0]=0;ChargeNum[1]=0;ChargeNum[2]=0;
  
  //int *Temp=new int[gridSizeXYZ];
	for(k=1;k<gridSizeZ-1;k++) {
		kgrid = k*gridSizeXY;
		for(j=1;j<gridSizeY-1;j++) {
			jgrid = kgrid+j*gridSizeX;
			for(i=1;i<gridSizeX-1;i++) {
				gridPoint = jgrid+i;
				BlackOrWhite=k+j+i+World->startBlackAndWhite;
				typeOfBorderPoint[gridPoint] = NoSingular;
        
        //if(BlackOrWhite%2==0)Temp[gridPoint]=0;
        //else Temp[gridPoint]=gridPoint;
				if(dielectric[0][gridPoint]!=dielectric[0][gridPoint-1] ||
							 dielectric[0][gridPoint]!=dielectric[1][gridPoint]||
							 dielectric[0][gridPoint]!=dielectric[1][gridPoint-gridSizeX] ||
							 dielectric[0][gridPoint]!=dielectric[2][gridPoint] ||
							 dielectric[0][gridPoint]!=dielectric[2][gridPoint-gridSizeXY])
				{
          //dielectricSum[gridPoint] = dielectric[0][gridPoint] + dielectric[0][gridPoint-1] + dielectric[1][gridPoint] + dielectric[1][gridPoint-gridSizeX] + dielectric[2][gridPoint] + dielectric[2][gridPoint-gridSizeXY];
					if(staticCharge[gridPoint]!=0){
						typeOfBorderPoint[gridPoint] = ChargeAndDielBoarder;
						if(BlackOrWhite%2==0)SingularNum[1]++;
						SingularNum[2]++;
					}
					else{
						typeOfBorderPoint[gridPoint] = DielBoarder;
						if(BlackOrWhite%2==0)DielBoarderNum[1]++;
						DielBoarderNum[2]++;
					}
				}
				else if(staticCharge[gridPoint]!=0) {
					typeOfBorderPoint[gridPoint] = Charge;
					if(BlackOrWhite%2==0)ChargeNum[1]++;
					ChargeNum[2]++;
				}
				if(typeOfBorderPoint[gridPoint] == NoSingular){
					if(BlackOrWhite%2==0)NoSingularNum[1]++;
					NoSingularNum[2]++;
				}
			}
		}
	}
  //if(World->MyRank==0)World->writeMapFloat("tmpde0.gz",World->Epsilon[1],gridSizeXYZ,1);
  
	DbgPrint1("Total: MyrankWorld->MyRank=%d icharge=% d    \niDielBoarder=%d iSingular=%d iNoSingular=%d\n",
						World->MyRank,ChargeNum[2],DielBoarderNum[2],SingularNum[2],NoSingularNum[2]);
	DbgPrint1("Blach cell:: MyrankWorld->MyRank=%d icharge=% d    \niDielBoarder=%d iSingular=%d iNoSingular=%d\n",
						World->MyRank,ChargeNum[1],DielBoarderNum[1],SingularNum[1],NoSingularNum[1]);
    
	if(!(IndexNoSingular = new int[NoSingularNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(IndexDielBoarder = new int[DielBoarderNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(IndexCharge = new int[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(IndexSingular = new int[SingularNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(PhiSingular = new float[SingularNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(PhiCharge = new float[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	for(i=0;i<ChargeNum[2];i++)
		PhiCharge[i]=0.0f;
	for(i=0;i<SingularNum[2];i++)
		PhiSingular[i]=0.0f;
	
	iNoSingular=-1;
	iSingular=-1; 
	iCharge=-1;
	iDielBoarder=-1; 
    //for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint=gridPoint+2){
	for(k=1;k<gridSizeZ-1;k++) {
		kgrid = k*gridSizeXY;
		for(j=1;j<gridSizeY-1;j++) {
			jgrid = kgrid+j*gridSizeX;
			for(i=1;i<gridSizeX-1;i++) {
				gridPoint = jgrid+i;
				BlackOrWhite=k+j+i+World->startBlackAndWhite;
				if(BlackOrWhite%2==0){
					switch(typeOfBorderPoint[gridPoint]){
						case NoSingular:
							iNoSingular++;
							IndexNoSingular[iNoSingular]=gridPoint;
							break;
						case Charge:
							iCharge++;
							IndexCharge[iCharge]=gridPoint;
							break;
						case DielBoarder:
							iDielBoarder++;
							IndexDielBoarder[iDielBoarder]=gridPoint;
							break;
						case ChargeAndDielBoarder:
							iSingular++;
							IndexSingular[iSingular]=gridPoint;
							break;
						case Boarder:
							break;   
					}
				}
			}
		}
	}
    //for(gridPoint=1;gridPoint<gridSizeXYZ;gridPoint=gridPoint+2){
	for(k=1;k<gridSizeZ-1;k++) {
		kgrid = k*gridSizeXY;
		for(j=1;j<gridSizeY-1;j++) {
			jgrid = kgrid+j*gridSizeX;
			for(i=1;i<gridSizeX-1;i++) {
				gridPoint = jgrid+i;
				BlackOrWhite=k+j+i+World->startBlackAndWhite;
				if(BlackOrWhite%2==1){
					switch(typeOfBorderPoint[gridPoint]){
						case NoSingular:
							iNoSingular++;
							IndexNoSingular[iNoSingular]=gridPoint;
							break;
						case Charge:
							iCharge++;
							IndexCharge[iCharge]=gridPoint;
							break;
						case DielBoarder:
							iDielBoarder++;
							IndexDielBoarder[iDielBoarder]=gridPoint;
							break;
						case ChargeAndDielBoarder:
							iSingular++;
							IndexSingular[iSingular]=gridPoint;
							break;
						case Boarder:
							break;   
					}
				}
			}
		}
	}
  //Free memory. SAVE it for children!
	if(typeOfBorderPoint!=NULL)
	{
		delete [] typeOfBorderPoint;
		typeOfBorderPoint=NULL;
	}
  
	//
  
	dielectricXDB = new float[DielBoarderNum[2]];
	dielectricYDB = new float[DielBoarderNum[2]];
	dielectricZDB = new float[DielBoarderNum[2]];
	dielectricXmDB = new float[DielBoarderNum[2]];
	dielectricYmDB = new float[DielBoarderNum[2]];
	dielectricZmDB = new float[DielBoarderNum[2]];
	dielectricZDBSUM = new float[DielBoarderNum[2]];
	dielectricXS = new float[SingularNum[2]];
	dielectricYS = new float[SingularNum[2]];
	dielectricZS = new float[SingularNum[2]];
	dielectricXmS = new float[SingularNum[2]];
	dielectricYmS = new float[SingularNum[2]];
	dielectricZmS = new float[SingularNum[2]];
	dielectricZSSUM = new float[SingularNum[2]];
  
	if(!(dielectricCh = new float[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
  
	if(!dielectricXDB||!dielectricYDB||!dielectricZDB||!dielectricZDBSUM||!dielectricXS||!dielectricYS||!dielectricZS||!dielectricZSSUM){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!dielectricXmDB||!dielectricYmDB||!dielectricZmDB||!dielectricXmS||!dielectricYmS||!dielectricZmS){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
  
	for(i=0;i<SingularNum[2];i++){
		gridPoint=IndexSingular[i];
		dielectricXS[i]=dielectric[0][gridPoint];
		dielectricYS[i]=dielectric[1][gridPoint];
		dielectricZS[i]=dielectric[2][gridPoint];
		dielectricXmS[i]=dielectric[0][gridPoint-1];
		dielectricYmS[i]=dielectric[1][gridPoint-gridSizeX];
		dielectricZmS[i]=dielectric[2][gridPoint-gridSizeXY];
		dielectricZSSUM[i]=dielectricXS[i]+dielectricYS[i]+dielectricZS[i]+dielectricXmS[i]+dielectricYmS[i]+dielectricZmS[i];
		dielectricZSSUM[i]=Relaxation/dielectricZSSUM[i];
	}
	for(i=0;i<ChargeNum[2];i++){
		gridPoint=IndexCharge[i];
    //denominator[gridPoint]=temp1;
		dielectricCh[i]=dielectric[0][gridPoint];
    //staticCharge[gridPoint]/=dielectricCh[i];
	}
	for(i=0;i<DielBoarderNum[2];i++){
		gridPoint=IndexDielBoarder[i];
		dielectricXDB[i]=dielectric[0][gridPoint];
		dielectricYDB[i]=dielectric[1][gridPoint];
		dielectricZDB[i]=dielectric[2][gridPoint];
		dielectricXmDB[i]=dielectric[0][gridPoint-1];
		dielectricYmDB[i]=dielectric[1][gridPoint-gridSizeX];
		dielectricZmDB[i]=dielectric[2][gridPoint-gridSizeXY];
		dielectricZDBSUM[i]=dielectricXDB[i]+dielectricYDB[i]+dielectricZDB[i]+dielectricXmDB[i]+dielectricYmDB[i]+dielectricZmDB[i];
		dielectricZDBSUM[i]=Relaxation/dielectricZDBSUM[i];
	}
	return EXIT_SUCCESS;
}
int PoissonSolver::SetQmobFromConcentration()
{
	int i, IType;
	int GrdPnt;
	float **C=World->C;
	double q;
	
	for(i=QmobNum[0];i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		Qmob[i]=(C[0][GrdPnt]-C[1][GrdPnt])/dielectricChMob[i];
	}
	for(i=QmobDielBoarderNum[0];i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		QmobDielBoarder[i]=(C[0][GrdPnt]-C[1][GrdPnt]);
	}
	for(i=QmobDielBoarderQstNum[0];i<QmobDielBoarderQstNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarderQst[i];
		QmobDielBoarderQst[i]=(C[0][GrdPnt]-C[1][GrdPnt]);
	}
	q=0.0;
	for(i=0;i<QmobNum[2];i++)
	{
		q+=Qmob[i]*dielectricChMob[i];
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		q+=QmobDielBoarder[i];
	}
	for(i=0;i<QmobDielBoarderQstNum[2];i++)
	{
		q+=QmobDielBoarderQst[i];
	}
	DbgPrint0("Mobile Charges induced: q=%f QmobNum=%d\n",(float)q/4/M_PI/World->GridScale);
	return EXIT_SUCCESS;
}
int PoissonSolver::SetQmobFromConcentrationDouble()
{
	int i, IType;
	int GrdPnt;
	double **C=World->CDouble;
	double q;
	
	for(i=QmobNum[0];i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		q=(C[0][GrdPnt]-C[1][GrdPnt])/dielectricChMob[i];
		Qmob[i]=q;
	}
	for(i=QmobDielBoarderNum[0];i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		q=(C[0][GrdPnt]-C[1][GrdPnt]);
		QmobDielBoarder[i]=q;
	}
	for(i=QmobDielBoarderQstNum[0];i<QmobDielBoarderQstNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarderQst[i];
		q=(C[0][GrdPnt]-C[1][GrdPnt]);
		QmobDielBoarderQst[i]=q;
	}
	q=0.0;
	for(i=0;i<QmobNum[2];i++)
	{
		q+=Qmob[i]*dielectricChMob[i];
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		q+=QmobDielBoarder[i];
	}
	for(i=0;i<QmobDielBoarderQstNum[2];i++)
	{
		q+=QmobDielBoarderQst[i];
	}
	DbgPrint0("Mobile Charges induced: q=%f QmobNum=%d\n",(float)q/4/M_PI/World->GridScale);
	return EXIT_SUCCESS;
}
int PoissonSolver::SetQmobForPNP()
{
	int i,j,k, i1, IType;
	int ix, iy, iz;
	int GrdPnt;
	int jgrid,kgrid;
	int GS_X = World->GridSize[0];
	int GS_Y = World->GridSize[1];
	int GS_Z = World->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int itmp1,itmp2;
	int iQmob,iQmob2;
	int iQmobDielBoarder,iQmobDielBoarder2;
	int iQmobDielBoarderQst,iQmobDielBoarderQst2;
	
	NodeIndexing* NIndexing=World->NIndexing;
	NodeIndex* NIndex=NIndexing->NIndex;
	
	NodeIndex DiffBoarderMask=NodeIndexing::DiffIon0BoarderMask;
	NodeIndex BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;
	NodeIndex DielBoarderMask=NodeIndexing::DielBoarderMask;
	
	NodeIndex specChargeMask=NodeIndexing::ChargeMask;
	NodeIndex ChargeDielBoarderMask=specChargeMask|NodeIndexing::DielBoarderMask;
	
	PNP_EXIT_FAIL_NULL(World->C,"Concentration Maps do not exist\n");
	for(IType=0;IType<World->NIonsTypes;IType++)
	{
		PNP_EXIT_FAIL_NULL(World->C[IType],"One of Concentration Maps does not exist\n");
	}
	PNP_EXIT_FAIL_NULL(World->NIndexing,"World->NIndexing does not exist\n");
	
	if(World->D==NULL)
	{
		for(IType=0;IType<World->NIonsTypes;IType++)
		{
			for(i=0;i<GS_XYZ;i++)
			{
					if(NIndexing->GetDiffFloat(IType,i)==0.0)World->C[IType][i]=0.0f;
			}
		}
	}
	else
	{
		for(IType=0;IType<World->NIonsTypes;IType++)
		{
			for(i=0;i<GS_XYZ;i++)
			{
				if(World->D[IType][i]==0.0)World->C[IType][i]=0.0f;
			}
		}
	}
	QmobNum[0]=0;
	QmobNum[1]=0;
	QmobNum[2]=0;
	QmobDielBoarderNum[0]=0;
	QmobDielBoarderNum[1]=0;
	QmobDielBoarderNum[2]=0;
	int countChargeDielBoarder=0;
	
	bool bCalcVolume=true;
	bool bQmobHere;
	//!@todo dising for 2 ions type
	for(k=1;k<GS_Z-1;k++)
	{
		kgrid = k*GS_XY;
		for(j=1;j<GS_Y-1;j++)
		{
			jgrid = kgrid+j*GS_X;
			for(i=1;i<GS_X-1;i++)
			{
				GrdPnt = jgrid+i;
				if(CalcVolume!=NULL)
					bCalcVolume=CalcVolume[GrdPnt];
				
				bQmobHere=false;
				for(IType=0;IType<World->NIonsTypes;IType++)
					if(World->C[IType][GrdPnt]>0.0)bQmobHere=true;
				
				if(bCalcVolume&&bQmobHere)
				{
					if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
					{
						QmobDielBoarderQstNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
						QmobDielBoarderQstNum[2]++;
					}
					else if(NIndex[GrdPnt]&DielBoarderMask)
					{
						QmobDielBoarderNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
						QmobDielBoarderNum[2]++;
					}
					else
					{
						QmobNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
						QmobNum[2]++;
					}
					
				}
			}
		}
	}
	if(countChargeDielBoarder>0)
		pnpError("countChargeDielBoarder=%d>0 this situation is not imlemented\n",countChargeDielBoarder);
	if(QmobNum[2]>0)
	{
		if(IndexQmob==NULL)IndexQmob=new int[QmobNum[2]];
		if(Qmob==NULL)Qmob=new float[QmobNum[2]];
		if(dielectricChMob==NULL)dielectricChMob=new float[QmobNum[2]];
	}
	if(QmobDielBoarderNum[2]>0)
	{
		if(IndexQmobDielBoarder==NULL)IndexQmobDielBoarder=new int[QmobDielBoarderNum[2]];
		if(QmobDielBoarder==NULL)QmobDielBoarder=new float[QmobDielBoarderNum[2]];
		dielectricXQmobDB = new float[QmobDielBoarderNum[2]];
		dielectricYQmobDB = new float[QmobDielBoarderNum[2]];
		dielectricZQmobDB = new float[QmobDielBoarderNum[2]];
		dielectricXmQmobDB = new float[QmobDielBoarderNum[2]];
		dielectricYmQmobDB = new float[QmobDielBoarderNum[2]];
		dielectricZmQmobDB = new float[QmobDielBoarderNum[2]];
		dielectricZQmobDBSUM = new float[QmobDielBoarderNum[2]];
	}
	if(QmobDielBoarderQstNum[2]>0)
	{
		if(IndexQmobDielBoarderQst==NULL)IndexQmobDielBoarderQst=new int[QmobDielBoarderQstNum[2]];
		if(QmobDielBoarderQst==NULL)QmobDielBoarderQst=new float[QmobDielBoarderQstNum[2]];
		if(QstQmobDielBoarderQst==NULL)QstQmobDielBoarderQst=new float[QmobDielBoarderQstNum[2]];
		dielectricXQmobDBQst = new float[QmobDielBoarderQstNum[2]];
		dielectricYQmobDBQst = new float[QmobDielBoarderQstNum[2]];
		dielectricZQmobDBQst = new float[QmobDielBoarderQstNum[2]];
		dielectricXmQmobDBQst = new float[QmobDielBoarderQstNum[2]];
		dielectricYmQmobDBQst = new float[QmobDielBoarderQstNum[2]];
		dielectricZmQmobDBQst = new float[QmobDielBoarderQstNum[2]];
		dielectricZQmobDBSUMQst = new float[QmobDielBoarderQstNum[2]];
	}
	
	fprintf(stdout,"		QmobNum:............... [%d,%d,%d]\n", QmobNum[0], QmobNum[1], QmobNum[2]);
	fprintf(stdout,"		QmobDielBoarderNum:.... [%d,%d,%d]\n", QmobDielBoarderNum[0], QmobDielBoarderNum[1], QmobDielBoarderNum[2]);
	fprintf(stdout,"		QmobDielBoarderNumQst:. [%d,%d,%d]\n", QmobDielBoarderQstNum[0], QmobDielBoarderQstNum[1], QmobDielBoarderQstNum[2]);
	
	float *Eps=World->NIndexing->Eps;
	iQmob=0;
	iQmob2=QmobNum[1];
	iQmobDielBoarder=0;
	iQmobDielBoarder2=QmobDielBoarderNum[1];
	iQmobDielBoarderQst=0;
	iQmobDielBoarderQst2=QmobDielBoarderQstNum[1];
	int QCount=0;
	for(k=1;k<GS_Z-1;k++) 
	{
		kgrid = k*GS_XY;
		for(j=1;j<GS_Y-1;j++) 
		{
			jgrid = kgrid+j*GS_X;
			for(i=1;i<GS_X-1;i++) 
			{
				GrdPnt = jgrid+i;
				if(NIndex[GrdPnt]&specChargeMask)
					QCount++;
				
				if(CalcVolume!=NULL)
					bCalcVolume=CalcVolume[GrdPnt];
				
				bQmobHere=false;
				for(IType=0;IType<World->NIonsTypes;IType++)
					if(World->C[IType][GrdPnt]>0.0)bQmobHere=true;
				
				if(bCalcVolume&&bQmobHere)
				{
					if(NIndex[GrdPnt]&BlackAndWhiteMask)
					{
						if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
						{
							IndexQmobDielBoarderQst[iQmobDielBoarderQst]=GrdPnt;
							dielectricXQmobDBQst[iQmobDielBoarderQst]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYQmobDBQst[iQmobDielBoarderQst]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZQmobDBQst[iQmobDielBoarderQst]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricXmQmobDBQst[iQmobDielBoarderQst]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYmQmobDBQst[iQmobDielBoarderQst]=Eps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZmQmobDBQst[iQmobDielBoarderQst]=Eps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricZQmobDBSUMQst[iQmobDielBoarderQst]=dielectricXQmobDBQst[iQmobDielBoarderQst]+dielectricYQmobDBQst[iQmobDielBoarderQst]+dielectricZQmobDBQst[iQmobDielBoarderQst]+dielectricXmQmobDBQst[iQmobDielBoarderQst]+dielectricYmQmobDBQst[iQmobDielBoarderQst]+dielectricZmQmobDBQst[iQmobDielBoarderQst];
							dielectricZQmobDBSUMQst[iQmobDielBoarderQst]=Relaxation/dielectricZQmobDBSUMQst[iQmobDielBoarderQst];
							QstQmobDielBoarderQst[iQmobDielBoarderQst]=NIndexing->Q[QCount];
							iQmobDielBoarderQst++;
						}
						else if(NIndex[GrdPnt]&DielBoarderMask)
						{
							IndexQmobDielBoarder[iQmobDielBoarder]=GrdPnt;
							dielectricXQmobDB[iQmobDielBoarder]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYQmobDB[iQmobDielBoarder]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZQmobDB[iQmobDielBoarder]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricXmQmobDB[iQmobDielBoarder]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYmQmobDB[iQmobDielBoarder]=Eps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZmQmobDB[iQmobDielBoarder]=Eps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricZQmobDBSUM[iQmobDielBoarder]=dielectricXQmobDB[iQmobDielBoarder]+dielectricYQmobDB[iQmobDielBoarder]+dielectricZQmobDB[iQmobDielBoarder]+dielectricXmQmobDB[iQmobDielBoarder]+dielectricYmQmobDB[iQmobDielBoarder]+dielectricZmQmobDB[iQmobDielBoarder];
							dielectricZQmobDBSUM[iQmobDielBoarder]=Relaxation/dielectricZQmobDBSUM[iQmobDielBoarder];
							iQmobDielBoarder++;
						}
						else
						{
							IndexQmob[iQmob]=GrdPnt;
							dielectricChMob[iQmob]=World->NIndexing->GetDielFloat(0,GrdPnt);
							iQmob++;
						}
					}
					else
					{
						if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
						{
							IndexQmobDielBoarderQst[iQmobDielBoarderQst2]=GrdPnt;
							dielectricXQmobDBQst[iQmobDielBoarderQst2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYQmobDBQst[iQmobDielBoarderQst2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZQmobDBQst[iQmobDielBoarderQst2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricXmQmobDBQst[iQmobDielBoarderQst2]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYmQmobDBQst[iQmobDielBoarderQst2]=Eps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZmQmobDBQst[iQmobDielBoarderQst2]=Eps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricZQmobDBSUMQst[iQmobDielBoarderQst2]=dielectricXQmobDBQst[iQmobDielBoarderQst2]+dielectricYQmobDBQst[iQmobDielBoarderQst2]+dielectricZQmobDBQst[iQmobDielBoarderQst2]+dielectricXmQmobDBQst[iQmobDielBoarderQst2]+dielectricYmQmobDBQst[iQmobDielBoarderQst2]+dielectricZmQmobDBQst[iQmobDielBoarderQst2];
							dielectricZQmobDBSUMQst[iQmobDielBoarderQst2]=Relaxation/dielectricZQmobDBSUMQst[iQmobDielBoarderQst2];
							QstQmobDielBoarderQst[iQmobDielBoarderQst2]=NIndexing->Q[QCount];
							iQmobDielBoarderQst2++;
						}
						else if(NIndex[GrdPnt]&DielBoarderMask)
						{
							IndexQmobDielBoarder[iQmobDielBoarder2]=GrdPnt;
							dielectricXQmobDB[iQmobDielBoarder2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYQmobDB[iQmobDielBoarder2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZQmobDB[iQmobDielBoarder2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricXmQmobDB[iQmobDielBoarder2]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							dielectricYmQmobDB[iQmobDielBoarder2]=Eps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
							dielectricZmQmobDB[iQmobDielBoarder2]=Eps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
							dielectricZQmobDBSUM[iQmobDielBoarder2]=dielectricXQmobDB[iQmobDielBoarder2]+dielectricYQmobDB[iQmobDielBoarder2]+dielectricZQmobDB[iQmobDielBoarder2]+dielectricXmQmobDB[iQmobDielBoarder2]+dielectricYmQmobDB[iQmobDielBoarder2]+dielectricZmQmobDB[iQmobDielBoarder2];
							dielectricZQmobDBSUM[iQmobDielBoarder2]=Relaxation/dielectricZQmobDBSUM[iQmobDielBoarder2];
							iQmobDielBoarder2++;
						}
						else
						{
							IndexQmob[iQmob2]=GrdPnt;
							dielectricChMob[iQmob2]=World->NIndexing->GetDielFloat(0,GrdPnt);
							iQmob2++;
						}
					}
				}
			}
		}
	}
	SetQmobFromConcentration();
	return EXIT_SUCCESS;
}


int PoissonSolver::Solve()
{
	if(solver==Auto)
	{
		if(World->NIndexing==NULL)return PoissonSolverAD();
		else return PoissonSolverNIB();
	}
	else if(solver==NodeIndexBased)
		return PoissonSolverNIB();
	else if(solver==ArrayDirect)
		return PoissonSolverAD();
	else if(solver==PNPC)
		return PoissonSolverW();
	return EXIT_FAILURE;
}
int PoissonSolver::PoissonSolverAD()
{
	float gridScale;
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	int iteration;
	int i,j,k;
	int gridPoint;
	float om1,om2,om2d6;
	float * potential;
	float * staticCharge;
	int gridSizeXY;
	int gridSizeXYZ;
	float temp1,temp2,temp3,temp4,temp5,temp6,temp7;
	float dynamicChargeFactor,IonStrengthFactor;
	float fpoh;
	double totalEnergyOld=totalEnergy;
	bool *PeriodicBoundaryCondition;
	float *positiveCharge;
	float *negativeCharge;
	
	float * dielectric[3];
	
	
	for(i=0;i<3;i++)
		dielectric[i] = World->Epsilon[i];
	
	if(!(World->Potential)||!(World->Qstat)){
		fprintf(stderr,"ERROR 110: Arrays NOT yet initialize\n");
		exit(105);
	}

	//cout<<"TRACE: CPoisson::poissonSolver() START\n";
	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;
	
	PeriodicBoundaryCondition=World->PBC;
	
	potential = World->Potential;
	if(World->C==NULL)
	{
		positiveCharge = NULL;
		negativeCharge = NULL;
	}
	else
	{
		positiveCharge = World->C[0];
		negativeCharge = World->C[1];
	}

	om2 = Relaxation;
	om1 = 1.0-om2;
	om2d6 = om2/6.0;
	
	//poissonBoltzmannSolver: assuming cation and anion concentration profiles are the same
	//poissonBoltzmannSolver: using cation concentration profile to calculate Debye lengths
	fpoh = 4.0*M_PI*gridScale;
	
	/*convert from charge density to M then convert from M to k'2
	* (or inverse debyle length squ					ared).
	*/
	dynamicChargeFactor = (float)1.0/(COANGS*4.0*M_PI*DFACT*DFACT);
	IonStrengthFactor = 1.0/(4.0*M_PI*gridScale);
	
	staticCharge=ChargeSum;
	
	//For speading calc denom first
	if(positiveCharge==0)
		for(k=0;k<gridSizeXYZ;k++){
		staticCharge[k]=World->Qstat[k];
		}
		else
			for(k=0;k<gridSizeXYZ;k++){
			staticCharge[k]=World->Qstat[k]+positiveCharge[k]-negativeCharge[k];
			}
			for(i=0;i<ChargeNum[2];i++){
				gridPoint=IndexCharge[i];
				staticCharge[gridPoint]/=dielectricCh[i];
			}
	
			for(iteration=1;iteration<=MaxIterations;iteration++) {
				for(j=0;j<=1;j++){
					for(i=SingularNum[j];i<SingularNum[j+1];i++){
						gridPoint=IndexSingular[i];
						temp1 = potential[gridPoint+1]*dielectricXS[i];
						temp2 = potential[gridPoint-1]*dielectricXmS[i];
						temp3 = potential[gridPoint+gridSizeX]*dielectricYS[i];
						temp4 = potential[gridPoint-gridSizeX]*dielectricYmS[i];
						temp5 = potential[gridPoint+gridSizeXY]*dielectricZS[i];
						temp6 = potential[gridPoint-gridSizeXY]*dielectricZmS[i];				
						temp1 = temp1+temp2;
						temp2 = temp3+temp4;
						temp3 = temp5+temp6;
						temp1 = temp1+temp2;
						temp2 = temp3+staticCharge[gridPoint];								
						temp1 = dielectricZSSUM[i]*(temp1+temp2);
						potential[gridPoint] = potential[gridPoint]*om1+temp1;
					}
					for(i=ChargeNum[j];i<ChargeNum[j+1];i++){
						gridPoint=IndexCharge[i];
						temp1 = potential[gridPoint+1]+potential[gridPoint-1];
						temp2 = potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX];
						temp3 = potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY];
						temp4 = potential[gridPoint]*om1;
						temp5 = om2d6*(temp1+temp2+temp3+staticCharge[gridPoint]);
				//temp6 = denominator[gridPoint]*(temp3+staticCharge[gridPoint]);
						potential[gridPoint] = temp4+temp5;
					}
					for(i=DielBoarderNum[j];i<DielBoarderNum[j+1];i++){
						gridPoint=IndexDielBoarder[i];
						temp1 = potential[gridPoint+1]*dielectricXDB[i];
						temp2 = potential[gridPoint-1]*dielectricXmDB[i];
						temp3 = potential[gridPoint+gridSizeX]*dielectricYDB[i];
						temp4 = potential[gridPoint-gridSizeX]*dielectricYmDB[i];
						temp5 = potential[gridPoint+gridSizeXY]*dielectricZDB[i];
						temp6 = potential[gridPoint-gridSizeXY]*dielectricZmDB[i];
						temp7 = potential[gridPoint]*om1;
						potential[gridPoint] = temp7+dielectricZDBSUM[i]*(temp1+temp2+temp3+temp4+temp5+temp6);
					}
					for(i=NoSingularNum[j];i<NoSingularNum[j+1];i++){
						gridPoint=IndexNoSingular[i];
						potential[gridPoint] = potential[gridPoint]*om1 + om2d6 * (potential[gridPoint+1] + potential[gridPoint-1] + potential[gridPoint+gridSizeX] + potential[gridPoint-gridSizeX] + potential[gridPoint+gridSizeXY] + potential[gridPoint-gridSizeXY]);
					}
					World->BorderExchange(potential);
				}
		
		//World->BorderExchange(potential);
		
				if((verbose&&(iteration%ConvergenceCheck==0))||iteration==MaxIterations)
				{
					temp1=totalEnergy;
					totalEnergy=CalculateEnergyPAD(fpoh,potential,staticCharge,dielectricCh, IndexCharge, IndexSingular, ChargeNum[2], SingularNum[2]);
					totalChange = fabs(totalEnergy-temp1);
					relativeChange=totalChange/totalEnergy;
					ConvFac=totalChange;
					
					if(verbose)
						pnpPrintGroup0("<PoissonIterations Nit=\"%8d\" E=\"%20.16le\" dE=\"%17.8lg\" rel.E=\"%17.8le\" ConvFac=\"%17.8le\"/>\n", iteration, totalEnergy, totalChange, relativeChange,ConvFac);

// #ifdef WIN32
//			 if(totalEnergy>1E13)
// #else
//			 if(totalEnergy>1E13||totalEnergy==NAN)
// #endif
					if(totalEnergy>1E13)
					{
						fprintf(stdout,"totalEnergy>1E13 || totalEnergy==NAN />\n", totalEnergy);
						return EXIT_FAILURE;
					}
					
					if(ConvFac<Convergence&&iteration>MinIterations)iteration = MaxIterations+1;
				}
			}
	//Calculate last energy
			totalEnergy=CalculateEnergyPAD(fpoh,potential,staticCharge,dielectricCh, IndexCharge, IndexSingular,ChargeNum[2], SingularNum[2]);
	
			if(verbose)
				pnpPrintGroup0("<PoissonFinal E=\"%.10e\" Eerr=\"%.10e\" Niter=\"%d\"/>\n", totalEnergy, totalChange, iteration-1);
			return EXIT_SUCCESS;
}
float PoissonSolver::CalculateEnergyPAD(float fpoh,float *Potential,float *StaticCharge,float *Epsilon,int *IndexCharge, int *IndexSingular,int ChargeNum,int SingularNum)
{
	int i,gridPoint;
	double Energy = 0.0;

	for(i=0;i<ChargeNum;i++){
		gridPoint=IndexCharge[i];
		Energy+= Potential[gridPoint]*StaticCharge[gridPoint]*Epsilon[i];
	}
	for(i=0;i<SingularNum;i++){
		gridPoint=IndexSingular[i];
		Energy+= Potential[gridPoint]*StaticCharge[gridPoint];
	}
#ifdef MPI_PARALLEL
	if(World->NProcs!=1)
	{
		MPI::Status	status;
		double EnergyProc;
		if(World->MyRank==0)
		{
			for(i=1;i<World->NProcs;i++)
			{
				pnpsapp->MyComGroup.Recv(&EnergyProc, 1, MPI::DOUBLE, i, SEND_ENERGY, status);
				Energy+=EnergyProc;
			}
		}
		else
		{
			pnpsapp->MyComGroup.Send(&Energy, 1, MPI::DOUBLE, 0, SEND_ENERGY);
		}
	}
#endif
	//Energy=Energy/(fpoh*2.0);
	World->SystemEnergy=Energy/(fpoh*2.0);
	return Energy/(fpoh*2.0);
}
int PoissonSolver::PoissonSolverW()
{
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	float gridScale;
	int iteration;
	int i,j,k;
	long gridPoint;
	float om1,om2,om2six;
	float * potential;
	float * potentialToo;
	float * dielectric[3];
	long gridSizeXY;
	long gridSizeXYZ;
	double totalChange;
	double totalEnergy;
	float change;
	long nextBorderPoint;
	long currentBorderCount;
	long * borderPoints;
	short * typeOfBorderPoint;
	float * om2InverseDielectricSum;
	float * positiveCharge;
	float * negativeCharge;
	float * chargeSum;
	long kgrid;
	long gridPointTemp;
	float temp1,temp2,temp3,temp4,temp5,temp6,temp7;
	int maxIterations;
	float convergence;
	long gridPointMax;
	float fpoh;

//	 assert(output!=NULL);
//	 assert(output->dielectricMap[0]!=NULL);
//	 assert(output->dielectricMap[1]!=NULL);
//	 assert(output->dielectricMap[2]!=NULL);
//	 assert(output->potentialMap!=NULL);

	assert(poissonSolverData!=NULL);
	assert(poissonSolverData->borderPoints!=NULL);

	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;

	borderPoints = poissonSolverData->borderPoints;
	typeOfBorderPoint = poissonSolverData->typeOfBorderPoint;
	om2InverseDielectricSum = poissonSolverData->om2InverseDielectricSum;

	maxIterations = poissonSolverData->maxIterations;
	convergence = poissonSolverData->convergence;
	
	potential = World->Potential;
	potentialToo = potential;
	
	chargeSum = World->Qstat;
	if(World->C!=NULL)
	{
		positiveCharge = World->C[0];
		negativeCharge = World->C[1];
	}
	else
	{
		positiveCharge = NULL;
		negativeCharge = NULL;
		poissonSolverData->hasDynamicCharges=false;
	}
	
	for(i=0;i<3;i++)
		dielectric[i] = World->Epsilon[i];

	om2 = poissonSolverData->relaxation;
	om2six = om2/6;
	om1 = 1-om2;

	fpoh = 4*M_PI*gridScale;

	if(maxIterations>=POISSON_SOLVER_CONVERGENCE_CHECK) {
		fprintf(stdout,"poissonSolver: energy change (kT)		total energy (kT)		iteration\n");
		fprintf(stdout,"poissonSolver:-----------------------------------------------------\n");
	}

	if(poissonSolverData->hasDynamicCharges) {
		nextBorderPoint = borderPoints[0];
		currentBorderCount = 0;
		
		for(k=1;k<gridSizeZ-1;k++) {
			kgrid = k*gridSizeXY;
			for(j=1;j<gridSizeY-1;j++) {
				gridPoint = kgrid+j*gridSizeX;
				gridPointMax = gridPoint+gridSizeX-1;
				if(gridPointMax<nextBorderPoint) {
					for(gridPoint++;gridPoint<gridPointMax;gridPoint++) {
						chargeSum[gridPoint]+=positiveCharge[gridPoint]-negativeCharge[gridPoint];
						chargeSum[gridPoint]*=om2six/dielectric[0][gridPoint];
					}
				}
				else {
					for(gridPoint++;gridPoint<gridPointMax;gridPoint++) {
						if(gridPoint<nextBorderPoint) {
							chargeSum[gridPoint]+=positiveCharge[gridPoint]-negativeCharge[gridPoint];
							chargeSum[gridPoint]*=om2six/dielectric[0][gridPoint];
						}
						else {
							chargeSum[gridPoint]+=positiveCharge[gridPoint]-negativeCharge[gridPoint];
							currentBorderCount++;
							nextBorderPoint = borderPoints[currentBorderCount];
						}
					}
				}
			}
		}

		for(iteration=1;iteration<=MaxIterations;iteration++) {
			nextBorderPoint = borderPoints[0];
			currentBorderCount = 0;

			if(iteration%ConvergenceCheck==0) {
				totalChange = 0;
				totalEnergy = 0;
				for(k=1;k<gridSizeZ-1;k++) {
					kgrid = k*gridSizeXY;
					for(j=1;j<gridSizeY-1;j++) {
						gridPoint = kgrid+j*gridSizeX;
						gridPointMax = gridPoint+gridSizeX-1;
						while(nextBorderPoint<gridPointMax) {
							for(gridPoint++;gridPoint<nextBorderPoint;gridPoint++) {
								temp1 = potential[gridPoint+1]+potential[gridPoint-1];
								temp2 = potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX];
								temp3 = potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY];
								temp4 = potential[gridPoint]*om2;
								temp5 = om2six*(temp1+temp2);
								temp6 = om2six*temp3-temp4;
								change = temp5+temp6+chargeSum[gridPoint];
								potentialToo[gridPoint] = potential[gridPoint]+change;
								totalChange+=change*chargeSum[gridPoint]*dielectric[0][gridPoint]/om2six;
								totalEnergy+=potential[gridPoint]*chargeSum[gridPoint]*dielectric[0][gridPoint]/om2six;
							}
							temp1 = potential[gridPoint+1]*dielectric[0][gridPoint];
							temp2 = potential[gridPoint-1]*dielectric[0][gridPoint-1];
							temp3 = potential[gridPoint+gridSizeX]*dielectric[1][gridPoint];
							temp4 = potential[gridPoint-gridSizeX]*dielectric[1][gridPoint-gridSizeX];
							temp5 = potential[gridPoint+gridSizeXY]*dielectric[2][gridPoint];
							temp6 = potential[gridPoint-gridSizeXY]*dielectric[2][gridPoint-gridSizeXY];
							temp7 = potential[gridPoint]*om2;
							change = om2InverseDielectricSum[currentBorderCount]*(temp1+temp2+temp3+temp4+temp5+temp6+chargeSum[gridPoint])-temp7;
							potentialToo[gridPoint] = potential[gridPoint]+change;
							totalChange+=change*chargeSum[gridPoint];
							totalEnergy+=potential[gridPoint]*chargeSum[gridPoint];
							currentBorderCount++;
							nextBorderPoint = borderPoints[currentBorderCount];
						}
						for(gridPoint++;gridPoint<gridPointMax;gridPoint++) {
							temp1 = potential[gridPoint+1]+potential[gridPoint-1];
							temp2 = potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX];
							temp3 = potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY];
							temp4 = potential[gridPoint]*om2;
							temp5 = om2six*(temp1+temp2);
							temp6 = om2six*temp3-temp4;
							change = temp5+temp6+chargeSum[gridPoint];
							potentialToo[gridPoint] = potential[gridPoint]+change;
							totalChange+=change*chargeSum[gridPoint]*dielectric[0][gridPoint]/om2six;
							totalEnergy+=potential[gridPoint]*chargeSum[gridPoint]*dielectric[0][gridPoint]/om2six;
						}
					}
				}
	
				totalEnergy = totalEnergy/(fpoh*2);
				totalChange = fabs(totalChange/(fpoh*2));
				fprintf(stdout,"poissonSolver: %18e		%17e	 %10i\n",totalChange,totalEnergy,iteration);
	
				if(totalChange<Convergence)
					iteration = MaxIterations+1;
			}
			else {
				for(k=1;k<gridSizeZ-1;k++) {
					kgrid = k*gridSizeXY;
					for(j=1;j<gridSizeY-1;j++) {
						gridPoint = kgrid+j*gridSizeX;
						gridPointMax = gridPoint+gridSizeX-1;
						while(nextBorderPoint<gridPointMax) {
							gridPointTemp = gridPoint;
							for(gridPoint++;gridPoint<nextBorderPoint;gridPoint+=2) {
								potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY])+chargeSum[gridPoint];
							}
							for(gridPoint=gridPointTemp+2;gridPoint<nextBorderPoint;gridPoint+=2) {
								potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY])+chargeSum[gridPoint];
							}
							gridPoint = nextBorderPoint;
							temp1 = potential[gridPoint+1]*dielectric[0][gridPoint];
							temp2 = potential[gridPoint-1]*dielectric[0][gridPoint-1];
							temp3 = potential[gridPoint+gridSizeX]*dielectric[1][gridPoint];
							temp4 = potential[gridPoint-gridSizeX]*dielectric[1][gridPoint-gridSizeX];
							temp5 = potential[gridPoint+gridSizeXY]*dielectric[2][gridPoint];
							temp6 = potential[gridPoint-gridSizeXY]*dielectric[2][gridPoint-gridSizeXY];
							temp7 = potential[gridPoint]*om1;
							potentialToo[gridPoint] = temp7+om2InverseDielectricSum[currentBorderCount]*(temp1+temp2+temp3+temp4+temp5+temp6+chargeSum[gridPoint]);
							currentBorderCount++;
							nextBorderPoint = borderPoints[currentBorderCount];
						}
						gridPointTemp = gridPoint;
						for(gridPoint++;gridPoint<gridPointMax;gridPoint+=2) {
							potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY])+chargeSum[gridPoint];
						}
						for(gridPoint=gridPointTemp+2;gridPoint<gridPointMax;gridPoint+=2) {
							potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY])+chargeSum[gridPoint];
						}
					}
				}
			}
		}
		
		nextBorderPoint = borderPoints[0];
		currentBorderCount = 0;
		
		for(k=1;k<gridSizeZ-1;k++) {
			kgrid = k*gridSizeXY;
			for(j=1;j<gridSizeY-1;j++) {
				gridPoint = kgrid+j*gridSizeX;
				gridPointMax = gridPoint+gridSizeX-1;
				if(gridPointMax<nextBorderPoint) {
					for(gridPoint++;gridPoint<gridPointMax;gridPoint++) {
						chargeSum[gridPoint]/=om2six/dielectric[0][gridPoint];
						chargeSum[gridPoint]-=positiveCharge[gridPoint]-negativeCharge[gridPoint];
					}
				}
				else {
					for(gridPoint++;gridPoint<gridPointMax;gridPoint++) {
						if(gridPoint<nextBorderPoint) {
							chargeSum[gridPoint]/=om2six/dielectric[0][gridPoint];
							chargeSum[gridPoint]-=positiveCharge[gridPoint]-negativeCharge[gridPoint];
						}
						else {
							chargeSum[gridPoint]-=positiveCharge[gridPoint]-negativeCharge[gridPoint];
							currentBorderCount++;
							nextBorderPoint = borderPoints[currentBorderCount];
						}
					}
				}
			}
		}
	}
	else {
		currentBorderCount = 0;
		
		for(gridPoint=borderPoints[currentBorderCount];gridPoint<gridSizeXYZ;gridPoint=borderPoints[++currentBorderCount]) {
			if(typeOfBorderPoint[currentBorderCount]==POISSON_SOLVER_CHARGED_POINT) {
				chargeSum[gridPoint]/=dielectric[0][gridPoint];
			}
		}
			
		for(iteration=1;iteration<=maxIterations;iteration++) {
			nextBorderPoint = borderPoints[0];
			currentBorderCount = 0;

			if(iteration%ConvergenceCheck==0) {
				totalChange = 0;
				totalEnergy = 0;
				for(k=1;k<gridSizeZ-1;k++) {
					kgrid = k*gridSizeXY;
					for(j=1;j<gridSizeY-1;j++) {
						gridPoint = kgrid+j*gridSizeX;
						gridPointMax = gridPoint+gridSizeX-1;
						while(nextBorderPoint<gridPointMax) {
							gridPointTemp = gridPoint;
							for(gridPoint++;gridPoint<nextBorderPoint;gridPoint+=2) {
								potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
							}
							for(gridPoint=gridPointTemp+2;gridPoint<nextBorderPoint;gridPoint+=2) {
								potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
							}
							gridPoint = nextBorderPoint;
							if(typeOfBorderPoint[currentBorderCount]==POISSON_SOLVER_CHARGED_POINT) {
								temp1 = potential[gridPoint+1]+potential[gridPoint-1];
								temp2 = potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX];
								temp3 = potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY];
								temp4 = potential[gridPoint]*om2;
								temp5 = om2six*(temp1+temp2);
								temp6 = om2six*(temp3+chargeSum[gridPoint]);
								change = temp5+temp6-temp4;
								potentialToo[gridPoint] = potential[gridPoint]+change;
								totalChange+=change*chargeSum[gridPoint]*dielectric[0][gridPoint];
								totalEnergy+=potential[gridPoint]*chargeSum[gridPoint]*dielectric[0][gridPoint];
							}
							else if(typeOfBorderPoint[currentBorderCount]==POISSON_SOLVER_UNCOMMON_DIELECTRIC_POINT) {
								temp1 = potential[gridPoint+1]*dielectric[0][gridPoint];
								temp2 = potential[gridPoint-1]*dielectric[0][gridPoint-1];
								temp3 = potential[gridPoint+gridSizeX]*dielectric[1][gridPoint];
								temp4 = potential[gridPoint-gridSizeX]*dielectric[1][gridPoint-gridSizeX];
								temp5 = potential[gridPoint+gridSizeXY]*dielectric[2][gridPoint];
								temp6 = potential[gridPoint-gridSizeXY]*dielectric[2][gridPoint-gridSizeXY];
								temp7 = potential[gridPoint]*om1;
								potentialToo[gridPoint] = temp7+om2InverseDielectricSum[currentBorderCount]*(temp1+temp2+temp3+temp4+temp5+temp6);
							}
							else {
								temp1 = potential[gridPoint+1]*dielectric[0][gridPoint];
								temp2 = potential[gridPoint-1]*dielectric[0][gridPoint-1];
								temp3 = potential[gridPoint+gridSizeX]*dielectric[1][gridPoint];
								temp4 = potential[gridPoint-gridSizeX]*dielectric[1][gridPoint-gridSizeX];
								temp5 = potential[gridPoint+gridSizeXY]*dielectric[2][gridPoint];
								temp6 = potential[gridPoint-gridSizeXY]*dielectric[2][gridPoint-gridSizeXY];
								change = om2InverseDielectricSum[currentBorderCount]*(temp1+temp2+temp3+temp4+temp5+temp6+chargeSum[gridPoint])-potential[gridPoint]*om2;
								potentialToo[gridPoint] = potential[gridPoint]+change;
								totalChange+=change*chargeSum[gridPoint];
								totalEnergy+=potential[gridPoint]*chargeSum[gridPoint];
							}
							currentBorderCount++;
							nextBorderPoint = borderPoints[currentBorderCount];
						}
						gridPointTemp = gridPoint;
						for(gridPoint=gridPointTemp+1;gridPoint<gridPointMax;gridPoint+=2) {
							potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
						}
						for(gridPoint=gridPointTemp+2;gridPoint<gridPointMax;gridPoint+=2) {
							potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
						}
					}
				}
	
				totalEnergy = totalEnergy/(fpoh*2);
				totalChange = fabs(totalChange/(fpoh*2));
				relativeChange = totalChange/totalEnergy;
				
				if(verbose)
					pnpPrintGroup0("<PoissonIterations Nit=\"%8d\" E=\"%20.16le\" dE=\"%17.8lg\" rel.E=\"%17.8le\"/>\n", iteration, totalEnergy, totalChange, relativeChange);
				if(totalChange<convergence)
					iteration = maxIterations+1;
			}
			else {
				for(k=1;k<gridSizeZ-1;k++) {
					kgrid = k*gridSizeXY;
					for(j=1;j<gridSizeY-1;j++) {
						gridPoint = kgrid+j*gridSizeX;
						gridPointMax = gridPoint+gridSizeX-1;
						while(nextBorderPoint<gridPointMax) {
							gridPointTemp = gridPoint;
							for(gridPoint++;gridPoint<nextBorderPoint;gridPoint+=2) {
								potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
							}
							for(gridPoint=gridPointTemp+2;gridPoint<nextBorderPoint;gridPoint+=2) {
								potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
							}
							gridPoint = nextBorderPoint;
							if(typeOfBorderPoint[currentBorderCount]==POISSON_SOLVER_CHARGED_POINT) {
								temp1 = potential[gridPoint+1]+potential[gridPoint-1];
								temp2 = potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX];
								temp3 = potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY];
								temp4 = potential[gridPoint]*om1;
								temp5 = om2six*(temp1+temp2);
								temp6 = om2six*(temp3+chargeSum[gridPoint]);
								potentialToo[gridPoint] = temp4+temp5+temp6;
							}
							else if(typeOfBorderPoint[currentBorderCount]==POISSON_SOLVER_UNCOMMON_DIELECTRIC_POINT) {
								temp1 = potential[gridPoint+1]*dielectric[0][gridPoint];
								temp2 = potential[gridPoint-1]*dielectric[0][gridPoint-1];
								temp3 = potential[gridPoint+gridSizeX]*dielectric[1][gridPoint];
								temp4 = potential[gridPoint-gridSizeX]*dielectric[1][gridPoint-gridSizeX];
								temp5 = potential[gridPoint+gridSizeXY]*dielectric[2][gridPoint];
								temp6 = potential[gridPoint-gridSizeXY]*dielectric[2][gridPoint-gridSizeXY];
								temp7 = potential[gridPoint]*om1;
								potentialToo[gridPoint] = temp7+om2InverseDielectricSum[currentBorderCount]*(temp1+temp2+temp3+temp4+temp5+temp6);
							}
							else {
								temp1 = potential[gridPoint+1]*dielectric[0][gridPoint];
								temp2 = potential[gridPoint-1]*dielectric[0][gridPoint-1];
								temp3 = potential[gridPoint+gridSizeX]*dielectric[1][gridPoint];
								temp4 = potential[gridPoint-gridSizeX]*dielectric[1][gridPoint-gridSizeX];
								temp5 = potential[gridPoint+gridSizeXY]*dielectric[2][gridPoint];
								temp6 = potential[gridPoint-gridSizeXY]*dielectric[2][gridPoint-gridSizeXY];
								potentialToo[gridPoint] = potential[gridPoint]*om1+om2InverseDielectricSum[currentBorderCount]*(temp1+temp2+temp3+temp4+temp5+temp6+chargeSum[gridPoint]);
							}
							currentBorderCount++;
							nextBorderPoint = borderPoints[currentBorderCount];
						}
						gridPointTemp = gridPoint;
						for(gridPoint=gridPointTemp+1;gridPoint<gridPointMax;gridPoint+=2) {
							potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
						}
						for(gridPoint=gridPointTemp+2;gridPoint<gridPointMax;gridPoint+=2) {
							potentialToo[gridPoint] = potential[gridPoint]*om1+om2six*(potential[gridPoint+1]+potential[gridPoint-1]+potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX]+potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY]);
						}
					}
				}
			}
		}

		currentBorderCount = 0;
		
		for(gridPoint=borderPoints[currentBorderCount];gridPoint<gridSizeXYZ;gridPoint=borderPoints[currentBorderCount++]) {
			if(typeOfBorderPoint[currentBorderCount]==POISSON_SOLVER_CHARGED_POINT) {
				chargeSum[gridPoint]*=dielectric[0][gridPoint];
			}
		}
	}


	return EXIT_SUCCESS;
}
int PoissonSolver::PoissonSolverNIB(bool ckenergy)
{
	//local vars
#ifndef PNPDOUBLE
	float gridScale = World->GridScale;
#else
	double gridScale = World->GridScale;
#endif
	int GS_X = World->GridSize[0];
	int GS_Y = World->GridSize[1];
	int GS_Z = World->GridSize[2];
	int GS_XY = GS_X*GS_Y;
	int GS_XYZ = GS_XY*GS_Z;
	bool *PeriodicBoundaryCondition = World->PBC;
#ifndef PNPDOUBLE
	//vars
	float fpoh = 4.0*M_PI*gridScale;
	float dynamicChargeFactor = 1.0/(COANGS*4.0*M_PI*DFACT*DFACT);
	float IonStrengthFactor = 1.0/(4.0*M_PI*gridScale);
	float om2 = Relaxation;
	float om1 = 1.0-om2;
	float om2d6 = om2/6.0;
	double totalEnergyOld=totalEnergy;
	//arrays
#else
	//vars
	double fpoh = 4.0*M_PI*gridScale;
	double dynamicChargeFactor = 1.0/(COANGS*4.0*M_PI*DFACT*DFACT);
	double IonStrengthFactor = 1.0/(4.0*M_PI*gridScale);
	double om2 = Relaxation;
	double om1 = 1.0-om2;
	double om2d6 = om2/6.0;
	
	double totalEnergyOld=totalEnergy;
	//arrays
	if(potential == NULL)
	{
		potential = new double[GS_XYZ];
		for(i=0;i<GS_XYZ;i++)
			potential[i]= (double)World->Potential[i];
	}
#endif

	
	//Check out is everything 
	if(!(World->Potential)){
		fprintf(stderr,"ERROR 110: Arrays NOT yet initialize\n");
		exit(105);
	}
	int countConvFacHistory=0;
	double *ConvFacHistory=NULL;
	if(ConvFacMaxHistory>1)
	{
		ConvFacHistory=new double[ConvFacMaxHistory];
		for(int i=0;i<ConvFacMaxHistory;i++)
			ConvFacHistory[i]=1e10;
	}

    PotBW->SetFromField(potential);

    int bwStrideX = PotBW->StrideX;
    int bwStrideXY = PotBW->StrideXY;
    int H_X = (GS_X+1)/2;
    int GS_X_odd = GS_X % 2;


    int *DielBoarderHmX=new int[DielBoarderNum[2]];

    for (int j = 0; j <= 1; j++) {
        for (int i = DielBoarderNum[j]; i < DielBoarderNum[j + 1]; i++) {
            int iu = IndexDielBoarder[i];
            int not_j = !j;

            int iz = iu / bwStrideXY;
            int iy = iu % bwStrideXY / bwStrideX;
            //int ix = GrdPnt % GS_X;
            //assert(GrdPnt == ix + iy * GS_X + iz * GS_XY);

            DielBoarderHmX[i] = (not_j + iy + iz) % 2 - 1;
        }
    }

	//if 0 continue if > 0 it is good if <0 it is bad
	int ReturnStatus = 0;
	int TotalIterations;
    const int GS_Step=4;

	#pragma omp parallel
	{
		int iteration;
		int i,j, k;
		int GrdPnt;

#ifndef PNPDOUBLE
		float tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
#else
		double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
#endif
		// check that black and white nodes matches
		int FirstNode = 0;
		GrdPnt = 1+1*GS_X+1*GS_XY;
		assert(World->NIndexing->NIndex[GrdPnt] & NodeIndexing::BlackAndWhiteMask == 1);
		//Iteration itself
		for (iteration = 1; iteration <= MaxIterations; iteration++)
		{
            float *PotU;// potential to update
            float *PotR;// potential to read

            int HmX, HpX;

            //calculation over black and white nodes
            for (j = 0; j <= 1; j++) {
                int not_j = !j;
                if (j == 0) {
                    // for calculation over black nodes
                    PotU = PotBW->B;
                    PotR = PotBW->W;
                }
                else {
                    // for calculation over white nodes
                    PotU = PotBW->W;
                    PotR = PotBW->B;
                    
                }

                    
#pragma omp for schedule(dynamic)
                for (int iz_big = 1; iz_big<GS_Z - 1; iz_big +=GS_Step) {
                    int iz_fin = (std::min)(GS_Z - 1, iz_big + GS_Step);
                    for (int iy_big = 1; iy_big<GS_Y - 1; iy_big+=GS_Step) {
                        int iy_fin=(std::min)(GS_Y - 1, iy_big + GS_Step);
                        
                        for (int iz = iz_big; iz<iz_fin; iz++) {
                            for (int iy = iy_big; iy<iy_fin; iy++) {
                                int iu0 = (j + iy + iz) % 2;
                                int iuS = iu0 + iy * bwStrideX + iz * bwStrideXY;
                                int iuF = iuS + H_X - 2 - (iu0 & GS_X_odd);
                            
                                HmX = (not_j + iy + iz) % 2 - 1;
                                HpX = HmX+1;

                                #pragma ivdep
                                for (int iu = iuS; iu<= iuF; ++iu) {
                                    PotU[iu] = PotU[iu] * om1 + om2d6 * (PotR[iu + HmX] + PotR[iu + HpX] + \
                                        PotR[iu - bwStrideX] + PotR[iu + bwStrideX] + PotR[iu - bwStrideXY] + PotR[iu + bwStrideXY]);
                                }
                            }
                        }
                    }
                }
                

#pragma omp for
                for (i = ChargeNum[j]; i < ChargeNum[j + 1]; i++) {
                    int iu = IndexCharge[i];
                    PotU[iu] += om2d6 * Qst[i];
                }
                

#pragma omp for
                for (i = DielBoarderNum[j]; i < DielBoarderNum[j + 1]; i++) {
                    int iu = IndexDielBoarder[i];
                    PotU[iu] += om2d6 * (\
                        PotR[iu + DielBoarderHmX[i]] * dielectricXmDB[i] + \
                        PotR[iu + DielBoarderHmX[i]+1] * dielectricXDB[i] + \
                        PotR[iu - bwStrideX] * dielectricYmDB[i] + \
                        PotR[iu + bwStrideX] * dielectricYDB[i] + \
                        PotR[iu - bwStrideXY] * dielectricZmDB[i] + \
                        PotR[iu + bwStrideXY] * dielectricZDB[i]);
                }
#pragma omp for
                for (i = QmobNum[j]; i < QmobNum[j + 1]; i++) {
                    int iu = IndexQmob[i];
                    PotU[iu] += om2d6 * Qmob[i];
                }

                PotBW->BorderExchange(World->PBC);
            }
            


			//checking and printing energy

			#pragma omp master
			{
				TotalIterations = iteration;
				if (ckenergy)
				{
					if ((verbose && (iteration%ConvergenceCheck == 0)) || iteration == MaxIterations)
					{
                        //PotBW->SetField(potential);
						CalcSystemEnergy(iteration);
						relativeChange = totalChange / totalEnergy;
						if (verbose)
						{
							if (iteration / ConvergenceCheck <= 1)
							{
								pnpPrintGroup0("P     =========================================================================\n");
								pnpPrintGroup0("P      %9s %22s %12s %12s %12s\n", "Iteration", "Energy,kT", "dE", "rel.E", "ConvFac");
								pnpPrintGroup0("P     -------------------------------------------------------------------------\n");
							}
							pnpPrintGroup0("P      %9d %22.14e %12.4e %12.4e %12.4e\n", iteration, totalEnergy, totalChange, relativeChange, ConvFac);
						}
						if (totalEnergy > 1E13)
						{
							pnpPrintGroup0("P     -------------------------------------------------------------------------\n");
							pnpPrintGroup0("P      ERROR: Poisson Solver has diverged, try smaller relaxation\n");
							pnpPrintGroup0("P     =========================================================================\n");
							ReturnStatus = -1;
						}
						if (ConvFacMaxHistory > 1)
						{
							ConvFacHistory[countConvFacHistory] = ConvFac;
							countConvFacHistory++;
							if (countConvFacHistory >= ConvFacMaxHistory)
								countConvFacHistory = 0;

							bool convereged = true;
							for (i = 0; i < ConvFacMaxHistory; i++)
								convereged = convereged && (ConvFacHistory[i] <= Convergence);

							if (convereged) {
								ReturnStatus = 1;
							}
						}
						else if (ConvFac<Convergence&&iteration>MinIterations) {
							ReturnStatus = 1;
						};
					}
				}
			}
            #pragma omp barrier
			if (ReturnStatus!=0) {
				break;
			}
		}
	}
    PotBW->SetField(potential);
	if(verbose&&ckenergy&&ReturnStatus>=0)
	{
		pnpPrintGroup0("P     -------------------------------------------------------------------------\n");
		pnpPrintGroup0("P      Results: E=%.14e Niter=%d\n", totalEnergy, TotalIterations);
		pnpPrintGroup0("P     =========================================================================\n");
		//pnpPrintGroup0("<PoissonFinal E=\"%.10e\" Eerr=\"%.10e\" Niter=\"%d\"/>\n", totalEnergy, totalChange, iteration-1);
	}
#ifdef PNPDOUBLE
	for(i=0;i<GS_XYZ;i++)
		World->Potential[i] = (float)potential[i];
	delete [] potential;
#endif
	
	DeleteCArray(ConvFacHistory);
	
	if((Convergence>0.0 && ConvFac>Convergence)|| ReturnStatus<0)
		return EXIT_FAILURE;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergy(int iteration)
{
// 	CalcSystemEnergyFloat(iteration);
	if(WayToCalcSystemEnergy==0)
	{
		//if(PMFCalculation0!=NULL)
		//	CalcSystemEnergy4PMF(iteration);
		//else
			CalcSystemEnergyStdDevPhi(iteration);
	}
	else
	{
		if(WayToCalcSystemEnergy==1)
			CalcSystemEnergyStdDevPhi(iteration);
	}
	//CalcSystemEnergyDouble(iteration);
		//CalcSystemEnergyMaxPhiChange(iteration);
// 	CalcSystemEnergyLongDouble(iteration);
// 	CalcSystemEnergy0(iteration);
// 	CalcSystemEnergy1(iteration);
 	//CalcSystemEnergy3(iteration);
	//if(iteration==MaxIterations)
	//	CalcSystemEnergyAnalizer(iteration);
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergyDouble(int iteration)
{
	int i,GrdPnt;
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double OldTotalEnergy=totalEnergy;
	float fpoh = 4.0*M_PI*World->GridScale;
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		
		tmp=double(potential[GrdPnt])*double(Qst[i])*double(dielectricCh[i]);
		EnergyCharge+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=double(QstS[i])*double(potential[GrdPnt]);
		EnergySingular+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		tmp=double(potential[GrdPnt])*double(Qmob[i])*double(dielectricChMob[i]);
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		tmp=double(potential[GrdPnt])*double(QmobDielBoarder[i]);
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	DbgPrint0("CSE_Double  :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n"
			,(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));
#ifdef MPI_PARALLEL
	int dest;
	double locEnergyCharge;
	double locEnergySingular;
	double locEnergyQmob;
	double locSumAbs,locSumSQ;
	if(World->MyRank==0)
	{
		for(dest=1;dest<World->NProcs;dest++)
		{
			pnpsapp->MyComGroup.Recv(&locEnergyCharge, 1, MPI::DOUBLE, dest, 1);
			pnpsapp->MyComGroup.Recv(&locEnergySingular, 1, MPI::DOUBLE, dest, 2);
			pnpsapp->MyComGroup.Recv(&locEnergyQmob, 1, MPI::DOUBLE, dest, 3);
			
			
			EnergyCharge+=locEnergyCharge;
			EnergySingular+=locEnergySingular;
			EnergyQmob+=locEnergyQmob;
			pnpsapp->MyComGroup.Recv(&locSumAbs, 1, MPI::DOUBLE, dest, 4);
			pnpsapp->MyComGroup.Recv(&locSumSQ, 1, MPI::DOUBLE, dest, 5);
			SumSQ+=locSumSQ;
			SumAbs+=locSumAbs;
		}
		DbgPrint0("CSE_DoubleF :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n"
			,(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));
	}
	else
	{
		pnpsapp->MyComGroup.Send(&EnergyCharge, 1, MPI::DOUBLE, 0, 1);
		pnpsapp->MyComGroup.Send(&EnergySingular, 1, MPI::DOUBLE, 0, 2);
		pnpsapp->MyComGroup.Send(&EnergyQmob, 1, MPI::DOUBLE, 0, 3);
		pnpsapp->MyComGroup.Send(&SumSQ, 1, MPI::DOUBLE, 0, 4);
		pnpsapp->MyComGroup.Send(&SumAbs, 1, MPI::DOUBLE, 0, 5);
	}
	
	pnpsapp->MyComGroup.Bcast(&EnergyCharge, 1, MPI::DOUBLE, 0);
	pnpsapp->MyComGroup.Bcast(&EnergySingular, 1, MPI::DOUBLE, 0);
	pnpsapp->MyComGroup.Bcast(&EnergyQmob, 1, MPI::DOUBLE, 0);
	pnpsapp->MyComGroup.Bcast(&SumSQ, 1, MPI::DOUBLE, 0);
	pnpsapp->MyComGroup.Bcast(&SumAbs, 1, MPI::DOUBLE, 0);
	
#endif

	totalEnergy=(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergyLongDouble(int iteration)
{
	int i,GrdPnt;
	static long double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	long double Dev=0.0;
	long double tmp,SumSQ=0.0,SumAbs=0.0;
	long double EnergyCharge = 0.0;
	long double EnergySingular = 0.0;
	long double EnergyQmob = 0.0;
	long double OldTotalEnergy=totalEnergy;
	long double  fpoh = 4.0*M_PI*World->GridScale;
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		tmp=potential[GrdPnt];
		tmp=Qst[i];
		tmp=dielectricCh[i];
		
		tmp=potential[GrdPnt]*Qst[i]*dielectricCh[i];
		EnergyCharge+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=QstS[i]*potential[GrdPnt];
		EnergySingular+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		tmp=potential[GrdPnt]*Qmob[i]*dielectricChMob[i];
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		tmp=potential[GrdPnt]*QmobDielBoarder[i];
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	DbgPrint0("CSE_LDouble :E=%.16Le Eq=%.10Le Esing=%.10Le EQmob=%.10lle\n",(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));

	totalEnergy=(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergyFloat(int iteration)
{
	int i,GrdPnt;
	static float oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	float Dev=0.0;
	float tmp,SumSQ=0.0,SumAbs=0.0;
	float EnergyCharge = 0.0;
	float EnergySingular = 0.0;
	float EnergyQmob = 0.0;
	float OldTotalEnergy=totalEnergy;
	float  fpoh = 4.0*M_PI*World->GridScale;
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		tmp=potential[GrdPnt];
		tmp=Qst[i];
		tmp=dielectricCh[i];
		
		tmp=potential[GrdPnt]*Qst[i]*dielectricCh[i];
		EnergyCharge+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=QstS[i]*potential[GrdPnt];
		EnergySingular+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		tmp=potential[GrdPnt]*Qmob[i]*dielectricChMob[i];
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		tmp=potential[GrdPnt]*QmobDielBoarder[i];
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	
	DbgPrint0("CSE_Float   :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n"
			,float((EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0)),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));

	totalEnergy=(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}

int PoissonSolver::CalcSystemEnergyMaxPhiChange(int iteration)
{
	int i,GrdPnt;
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double OldTotalEnergy=totalEnergy;
	double  fpoh = 4.0*M_PI*World->GridScale;
	float dphi=0.0,maxdphi=0.0;
	int GrdPntMaxDPhi=-1;
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		
		tmp=double(potential[GrdPnt])*double(Qst[i])*double(dielectricCh[i]);
		EnergyCharge+=tmp;
		dphi=fabs(potential[GrdPnt]-PhiCharge[i]);
		PhiCharge[i]=potential[GrdPnt];
		if(dphi>maxdphi)
		{
			maxdphi=dphi;
			GrdPntMaxDPhi=GrdPnt;
		}
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=double(QstS[i])*double(potential[GrdPnt]);
		EnergySingular+=tmp;
		dphi=fabs(potential[GrdPnt]-PhiSingular[i]);
		PhiSingular[i]=potential[GrdPnt];
		if(dphi>maxdphi)
		{
			maxdphi=dphi;
			GrdPntMaxDPhi=GrdPnt;
		}
	}

	//DbgPrint0("CSE_0       :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n",(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));
	DbgPrint0("CSE_MaxPhiChange %.5e at %d\n",maxdphi,GrdPntMaxDPhi);
	totalEnergy=(EnergyCharge+EnergySingular)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac=maxdphi;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergyStdDevPhi(int iteration)
{
	int i;
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double OldTotalEnergy=totalEnergy;
	double  fpoh = 4.0*M_PI*World->GridScale;
	float dphi=0.0,maxdphi=0.0;
	int GrdPntMaxDPhi=-1;
	

    float *Pot;

    for (int j = 0; j <= 1; j++) {
        if (j == 0) {
            Pot = PotBW->B;
        }
        else {
            Pot = PotBW->W;
        }
	    for(i=ChargeNum[j];i<ChargeNum[j+1];i++)
	    {
            int ip = IndexCharge[i];
		
		    tmp=double(Pot[ip])*double(Qst[i])*double(dielectricCh[i]);
		    EnergyCharge+=tmp;
		    dphi= Pot[ip] -PhiCharge[i];
		    SumSQ+=dphi*dphi;
		    PhiCharge[i]= Pot[ip];
	    }
    }

	//DbgPrint0("CSE_0       :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n",(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));
	//DbgPrint0("CSE_MaxPhiChange %.5e at %d\n",maxdphi,GrdPntMaxDPhi);
	totalEnergy=EnergyCharge/(fpoh*2.0);
	int ChargedNodes=ChargeNum[2];
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac=sqrt(SumSQ/double(ChargedNodes));
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergy0(int iteration)
{
	int i,GrdPnt;
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double OldTotalEnergy=totalEnergy;
	double  fpoh = 4.0*M_PI*World->GridScale;
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		
		tmp=double(potential[GrdPnt])*double(Qst[i])*double(dielectricCh[i]);
		EnergyCharge+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=double(QstS[i])*double(potential[GrdPnt]);
		EnergySingular+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		tmp=double(potential[GrdPnt])*double(Qmob[i])*double(dielectricChMob[i]);
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		tmp=double(potential[GrdPnt])*double(QmobDielBoarder[i]);
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	DbgPrint0("CSE_0       :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n",(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));

	totalEnergy=(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergy3(int iteration)
{
	int i,GrdPnt;
	int AvrOver=10;
	static int Counter=0;
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,tmp2,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double OldTotalEnergy=totalEnergy;
	double  fpoh = 4.0*M_PI*World->GridScale;
	
	
	if(ChargesEnergy==NULL)
	{
		ChargesEnergy=new double[AvrOver];
		for(i=0;i<AvrOver;i++)ChargesEnergy[i]=-1e300;
		Counter=0;
	}
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		
		tmp=double(potential[GrdPnt])*double(Qst[i])*double(dielectricCh[i]);
		EnergyCharge+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=double(QstS[i])*double(potential[GrdPnt]);
		EnergySingular+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		tmp=double(potential[GrdPnt])*double(Qmob[i])*double(dielectricChMob[i]);
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		tmp=double(potential[GrdPnt])*double(QmobDielBoarder[i]);
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	
	

	totalEnergy=(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	ChargesEnergy[Counter]=totalEnergy;
	Counter++;
	if(Counter>=AvrOver)Counter=0;
	
	tmp=0.0;
	tmp2=0.0;
	for(i=0;i<AvrOver;i++)
	{
		tmp+=ChargesEnergy[i];
		tmp2+=ChargesEnergy[i]*ChargesEnergy[i];
	}
	tmp/=AvrOver;
	tmp2/=AvrOver;
	tmp2 = sqrt(fabs(tmp2-tmp*tmp));
	
	DbgPrint0("CSE_3       :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n",(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));
	DbgPrint0("CSE_3       :E=%.16e Eq=%.10e\n",tmp,tmp2);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	
	
	return EXIT_SUCCESS;
}
int compare_doubles( const void* a, const void* b ) 
{
	double* arg1 = (double*) a;
	double* arg2 = (double*) b;
	if( *arg1 < *arg2 ) return -1;
	else if( *arg1 == *arg2 ) return 0;
	else return 1;
}
int PoissonSolver::CalcSystemEnergy1(int iteration)
{
	int i,GrdPnt;
	int TotCharge=ChargeNum[2]+SingularNum[2]+QmobNum[2]+QmobDielBoarderNum[2],TotChargeCount;
	
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double OldTotalEnergy=totalEnergy;
	double  fpoh = 4.0*M_PI*World->GridScale;
	
	if(ChargesEnergy==NULL)
	{
		ChargesEnergy=new double[TotCharge];
	}
	TotChargeCount=0;
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		
		tmp=double(potential[GrdPnt])*double(Qst[i])*double(dielectricCh[i]);
		
		ChargesEnergy[TotChargeCount]=tmp;
		TotChargeCount++;
		
		EnergyCharge+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
		
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=double(QstS[i])*double(potential[GrdPnt]);
		ChargesEnergy[TotChargeCount]=tmp;
		TotChargeCount++;
		
		EnergySingular+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		tmp=double(potential[GrdPnt])*double(Qmob[i])*double(dielectricChMob[i]);
		ChargesEnergy[TotChargeCount]=tmp;
		TotChargeCount++;
		
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		tmp=double(potential[GrdPnt])*double(QmobDielBoarder[i]);
		ChargesEnergy[TotChargeCount]=tmp;
		TotChargeCount++;
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	//sort ChargesEnergy
	qsort( ChargesEnergy, TotCharge, sizeof(double), compare_doubles );
	tmp=0.0;
	for(i=0;i<TotCharge;i++)
	{
		tmp+=ChargesEnergy[i];
	}
	DbgPrint0("CSE_1       :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n",tmp/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));

	totalEnergy=(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::CalcSystemEnergyAnalizer(int iteration)
{
	int i,GrdPnt;
	static long double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	long double Dev=0.0;
	long double tmp,SumSQ=0.0,SumAbs=0.0;
	long double EnergyCharge = 0.0;
	long double EnergySingular = 0.0;
	long double EnergyQmob = 0.0;
	long double OldTotalEnergy=totalEnergy;
	long double ftmp, tmpAvr=0.0,tmpMin=0.0,tmpMax=0.0,tmpSD=0.0;
	long double  fpoh = 4.0*M_PI*World->GridScale;
	
	if(ChargeNum[2]>0)
	{
		tmp=fabs(potential[IndexCharge[0]]*Qst[0]*dielectricCh[0]);
		tmpMin=tmp;
		tmpMax=tmp;
	}
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		
		tmp=potential[GrdPnt]*Qst[i]*dielectricCh[i];
		EnergyCharge+=tmp;
		
		ftmp=fabs(tmp);
		tmpAvr+=ftmp;
		tmpSD+=ftmp*ftmp;
		
		if(ftmp<tmpMin)tmpMin=ftmp;
		if(ftmp>tmpMax)tmpMax=ftmp;
		
		//DbgPrint0("%.16Le %.16Le\n",tmp,ftmp);
		
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=QstS[i]*potential[GrdPnt];
		EnergySingular+=tmp;
		
		ftmp=fabs(tmp);
		tmpAvr+=ftmp;
		tmpSD+=ftmp*ftmp;
		
		if(ftmp<tmpMin)tmpMin=ftmp;
		if(ftmp>tmpMax)tmpMax=ftmp;
		
		//DbgPrint0("%.16Le %.16Le\n",tmp,ftmp);
		
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobNum[2];i++)
	{
		GrdPnt=IndexQmob[i];
		tmp=potential[GrdPnt]*Qmob[i]*dielectricChMob[i];
		
		ftmp=fabs(tmp);
		tmpAvr+=ftmp;
		tmpSD+=ftmp*ftmp;
		
		if(ftmp<tmpMin)tmpMin=ftmp;
		if(ftmp>tmpMax)tmpMax=ftmp;
		
		//DbgPrint0("%.16Le %.16Le\n",tmp,ftmp);
		
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	for(i=0;i<QmobDielBoarderNum[2];i++)
	{
		GrdPnt=IndexQmobDielBoarder[i];
		tmp=potential[GrdPnt]*QmobDielBoarder[i];
		
		ftmp=fabs(tmp);
		tmpAvr+=ftmp;
		tmpSD+=ftmp*ftmp;
		
		if(ftmp<tmpMin)tmpMin=ftmp;
		if(ftmp>tmpMax)tmpMax=ftmp;
		
		//DbgPrint0("%.16Le %.16Le\n",tmp,ftmp);
		
		EnergyQmob+=tmp;
		SumSQ+=tmp*tmp;
		SumAbs+=fabs(tmp);
	}
	tmpAvr=tmpAvr/(ChargeNum[2]+SingularNum[2]+QmobNum[2]+QmobDielBoarderNum[2]);
	tmpSD=tmpSD/(ChargeNum[2]+SingularNum[2]+QmobNum[2]+QmobDielBoarderNum[2]);
	tmpSD=sqrt(tmpSD-tmpAvr*tmpAvr);
	tmpAvr/=fpoh;
	tmpSD/=fpoh;
	tmpMin/=fpoh;
	tmpMax/=fpoh;
	DbgPrint0("fpoh %.7Le\n",fpoh);
	DbgPrint0("CSE_Analize :Avr(|Ei|)=%.5Le SD(|Ei|)=%.5Le Min(|Ei|)=%.10Le Max(|Ei|)=%.10Le\n",tmpAvr,tmpSD,tmpMin,tmpMax);
	DbgPrint0("CSE_LDouble :E=%.16Le Eq=%.10Le Esing=%.10Le EQmob=%.10lle\n",(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));

	totalEnergy=(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac = sqrt(fabs(oldSumSQ-SumSQ))/(fpoh*2.0);
	Dev=oldSumAbs-SumAbs;
	//if(oldSumSQ==0.0)fprintf(stdout,"PT Nit E dE RMSD dEabs\n");
	//if(verbose)fprintf(stdout,"PT %8d %17.12lg %17.12lg %17.12lg %17.12lg\n", iteration, totalEnergy, totalEnergy-oldTotEn,totalChange, Dev);
	
	oldSumSQ=SumSQ;
	oldSumAbs=SumAbs;
	oldTotEn=totalEnergy;
	World->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
int PoissonSolver::GuessNumberOfIteration()
{
	int GrdPnt;
	int i,j,k;
	int ix,iy,iz;
	float *sn[3];
	double tmp;
	for(i=0;i<3;i++)
	{
		sn[i]=new float [World->GridSize[i]];
		if(World->PBC[i])
		{
			for(j=0;j<World->GridSize[i];j++)
			{
				sn[i][j]=1.0/sqrt(float(World->GridSize[i]));
			}
		}
		else
		{
			sn[i][0]=0.0;
			sn[i][World->GridSize[i]-1]=0.0;
			for(j=1;j<World->GridSize[i]-1;j++)
			{
				tmp=M_PI*float(j)/float(World->GridSize[i]-1);
				sn[i][j]=sqrt(2.0)*sin(tmp)/sqrt(float(World->GridSize[i]-1));
			}
		}
	}
#ifndef PNPDOUBLE
	float *tmpPot=new float[GS_XYZ];
#else
	double *tmpPot=new double[GS_XYZ];
#endif
	
	//float *tmpPot2=new float[GS_XYZ];
	
	for(ix=0;ix<World->GridSize[0];ix++)
		for(iy=0;iy<World->GridSize[1];iy++)
			for(iz=0;iz<World->GridSize[2];iz++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		tmpPot[GrdPnt]=sn[0][ix]*sn[1][iy]*sn[2][iz];
		//tmpPot2[GrdPnt]=tmpPot[GrdPnt];
	}
#ifndef PNPDOUBLE
	float *PotPointer=potential;
#else
	double *PotPointer=potential;
#endif
	
	potential=tmpPot;
	int solverLoc;
	if(solver==Auto)
	{
		if(World->NIndexing==NULL)solverLoc=ArrayDirect;
		else solverLoc=NodeIndexBased;
	}
	else
	{
		solverLoc=solver;
	}
	
	if(solverLoc==NodeIndexBased)
	{
		int oldQmobNum[3];
		int oldChargeNum[3];
		int oldConvFacMaxHistory = ConvFacMaxHistory;
		int oldMaxIterations = MaxIterations;

		for (int ii = 0; ii < 3; ii++) {
			oldQmobNum[ii] = QmobNum[ii];
			oldChargeNum[ii] = ChargeNum[ii];
			QmobNum[ii]=0;
			ChargeNum[ii]=0;
		}
		ConvFacMaxHistory = 1;
		MaxIterations = 1;

		PoissonSolverNIB(false);

		for (int ii = 0; ii < 3; ii++) {
			QmobNum[ii] = oldQmobNum[ii];
			ChargeNum[ii] = oldChargeNum[ii];
		}
		ConvFacMaxHistory = oldConvFacMaxHistory;
		MaxIterations = oldMaxIterations;
	}
	else if(solverLoc==ArrayDirect)
	{
		pnpError("Cann't estemate number of interaction for ArrayDirect solver\n");
		return 0;
	}
	else
	{
		pnpError("Unknown type of solver\n");
		return 0;
	}
	
	NodeIndex* NIndex=World->NIndexing->NIndex;
	unsigned int BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;
	tmp=0.0;
	for(ix=0;ix<World->GridSize[0];ix++)
		for(iy=0;iy<World->GridSize[1];iy++)
			for(iz=0;iz<World->GridSize[2];iz++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		if(!(NIndex[GrdPnt]&BlackAndWhiteMask))
		{
		tmp+=potential[GrdPnt]*sn[0][ix]*sn[1][iy]*sn[2][iz];
		}
	}
	float spec=2.0*tmp;
    int iter;
    if(spec>=1.0){
        pnpPrint("Error: Found Gauss-Seidel Spectral Radius is bigger then 1: %f\n", spec);
        pnpPrint("       Please select relaxation parameter and number of iteration manually\n");
        exit(1);
    }
    else {
	    pnpPrint("Gauss-Seidel Spectral Radius: %f\n",spec);
	    iter=int(7.8/log(1.0 + sqrt(1.0-spec)));
	    pnpPrint("Estimated Number of Iterations to Convergence: %d\n",iter);
	    DeleteCVecArray(sn,3);
	    delete [] tmpPot;
	    potential=PotPointer;
	    Relaxation=2.0/(1.0 + sqrt(1.0 - spec));
	    pnpPrint("Estimated Relaxation Coefficient: %f\n",Relaxation);
	    SetRelaxation(Relaxation);
    }
	return iter;
}
