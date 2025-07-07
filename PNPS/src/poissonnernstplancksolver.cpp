//
// C++ Implementation: poissonnernstplancksolver
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

#include "poissonnernstplancksolver.h"
#include "nernstplanksolver.h"
#include "poissonsolver.h"
#include "tinyxml.h"
#include "pnpdebug.h"
#include "contworld.h"
#include "math.h"
#include "pnpconstants.h"
#include "mapio.h"
#include "pnputil.h"

#include <iostream>
#include <Python.h>

using namespace std;

PoissonNernstPlanckSolver::PoissonNernstPlanckSolver()
{
	PMFWeightModeStr.push_back("None");
	PMFWeightModeStr.push_back("Linear");
	InitZero();
}


PoissonNernstPlanckSolver::~PoissonNernstPlanckSolver()
{
	
	Clear();
}
int PoissonNernstPlanckSolver::InitZero()
{
	PnpsObject::SetName("PoissonNernstPlanckSolver");
	World=NULL;
	Poisson=NULL;
	NernstPlank=NULL;
	ItotErr=5550.0;
	IposErr=5550.0;
	InegErr=5550.0;
	positiveCurrentProfile=NULL;
	negativeCurrentProfile=NULL;
	currentDimension=2;
	ConvergenceCheck=10;PMFWeightMode=0;
	bLimitCurrentCalc=false;
	LimitCurrentCalcZ[0]=0;
	LimitCurrentCalcZ[1]=0;
	SaveMemory=false;
	bDouble=false;
	verbose=true;
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::Clear()
{
	DeleteCArray(positiveCurrentProfile);
	DeleteCArray(negativeCurrentProfile);
	DeleteObjByPnt(Poisson);
	DeleteObjByPnt(NernstPlank);
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::LoadParamFromPyDict(PyObject *dict)
{
	Clear();
	int i,gridPoint;
	//Read Primary Parameters
	// =haPyDict_GetItemValueAsBool(dict,"",);
	// =haPyDict_GetItemValueAsInt(dict,"",);
	// =haPyDict_GetItemValueAsFloat(dict,"",);

	MaxIterations=haPyDict_GetItemValueAsInt(dict,"MaxIterations",MaxIterations);
	verbose=haPyDict_GetItemValueAsBool(dict,"Verbose",verbose);
	ConvergenceCheck=haPyDict_GetItemValueAsInt(dict,"ConvergenceCheck",ConvergenceCheck);
	PMFWeightMode=haPyDict_GetItemValueAsInt(dict,"PMFWeightMode",PMFWeightMode);
	SaveMemory=haPyDict_GetItemValueAsBool(dict,"SaveMemory",SaveMemory);

	currentDimension=2;
	//Scale Parameters to Internal Units
	//Calculate Derivative Parameters
	if(Poisson==NULL)
	{
		Poisson = new PoissonSolver();
		Poisson->LoadParamFromPyDict(PyDict_GetItemString(dict,"P_Param"));
	}
	else Poisson->LoadParamFromPyDict(PyDict_GetItemString(dict,"P_Param"));
  
	if(NernstPlank==NULL)
	{
		NernstPlank = new NernstPlankSolver();
		NernstPlank->LoadParamFromPyDict(PyDict_GetItemString(dict,"NP_Param"));
	}
	else NernstPlank->LoadParamFromPyDict(PyDict_GetItemString(dict,"NP_Param"));
	
	PyObject *p=NULL;
	if((p=PyDict_GetItemString(dict,"LimitCurrentCalcZ"))!=NULL)
	{
		bLimitCurrentCalc=true;
		haPy_SetCFloatArrFromListOfFloat(p,LimitCurrentCalcZ);
	}
	else
	{
		bLimitCurrentCalc=false;
		LimitCurrentCalcZ[0]=0;
		LimitCurrentCalcZ[1]=0;
	}
	bDouble=haPyDict_GetItemValueAsBool(dict,"bDouble",bDouble);
	
	Poisson->verbose=verbose;
	NernstPlank->verbose=verbose;
	
	
	ShowParameters();
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::SetContWorld(ContWorld* _world)
{
	int i;
	World=_world;
	if(Poisson!=NULL)
	{
		Poisson->SetContWorld(World);
	}
	if(NernstPlank!=NULL)
	{
		NernstPlank->SetContWorld(World);
	}
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::InitSolver()
{
	pnpPrint("\nPoisson-Nernst-Plank Solver\n");
	int iteration;
	float Tolerance;
	int i,j,k,i1,i2,gridPoint;
	float gridScale;
	int * gridSize;
  
	float *PotentialTMP, *Potential;
	float om1,om2;
	int iCurrentDimension;
	float currentSum,currentCoordinate;
	float currentDev,current;
	float currentDevPos,currentPos;
	float currentDevNeg,currentNeg;
	float maxCurrent,maxCurrentCoordinate;
	
	iCurrentDimension = currentDimension;

	gridScale = World->GridScale;
	gridSize = World->GridSize;
	Potential = World->Potential;
	
	int IDem;
	if( World->MyRank==0) IDem=World->GridSizeGlobal[iCurrentDimension]-1;
	else IDem=gridSize[iCurrentDimension]-1;
	
	if(!(positiveCurrentProfile = new double[IDem])){
		pnpError("ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(negativeCurrentProfile = new double[IDem])){
		pnpError("ERROR 204: No memory available\n");
		exit(204);
	}
	if(World->PMF!=NULL)pnpPrint("Calculate with PMF correction\n");
	
	if(World->D!=NULL&&World->NIndexing!=NULL)
		World->NIndexing->SetIonAccess(World->D);
	
	if(World->C!=NULL&&World->NIndexing!=NULL)Poisson->SetQmobForPNP();
	
	NernstPlank->InitSolver();
	Poisson->InitSolver();
	//if(SaveMemory)
	//	DeleteObjByPnt(World->NIndexing);
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::Solve()
{
	if(bDouble)
		SolveDouble();
	else
		SolveSingle();
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::SolveSingle()
{
	pnpPrint("\nPoisson-Nernst-Plank Solver\n");
	int iteration;
	float Tolerance;
	int i,j,k,i1,i2,gridPoint;
	float gridScale;
	int * gridSize;
	float *PotentialTMP, *Potential;
  //assert(pnpSolverData!=NULL);
#ifdef MPI_PARALLEL
	MPI::Status  status;
#endif
	
	gridScale = World->GridScale;
	gridSize = World->GridSize;
  //PotentialTMP = World->PotentialTMP;
	Potential = World->Potential;
	
	
	for(iteration=1;iteration<=MaxIterations;iteration++) 
	{
		int PoissonStatus,NernstPlankStatus;

		assert(Poisson);
		PoissonStatus=Poisson->Solve();
		if(PoissonStatus!=EXIT_SUCCESS)
		{
			pnpPrintGroup0("PNPSR -------------------------------------------------------------------------\n");
			pnpPrintGroup0("PNPSR  ERROR: Poisson solver has failed.\n");
			pnpPrintGroup0("PNPSR =========================================================================\n");
			return EXIT_FAILURE;
		}
		
		if(PMFWeightMode==1)
		{
			NernstPlank->PMFWeight=double(iteration)/double(MaxIterations);
			pnpPrint0("PMFWeight=%f\n",NernstPlank->PMFWeight);
		}
		NernstPlankStatus=NernstPlank->Solve();
    //NernstPlankStatus=NernstPlank->NernstPlanckSolverN2tmp();
		if(Poisson->solver!=Poisson->ArrayDirect)
			Poisson->SetQmobFromConcentration();
		if(NernstPlankStatus!=EXIT_SUCCESS)
		{
			pnpPrintGroup0("PNPSR -------------------------------------------------------------------------\n");
			pnpPrintGroup0("PNPSR  ERROR: Nernst-Planck solver has failed.\n");
			pnpPrintGroup0("PNPSR =========================================================================\n");
			return EXIT_FAILURE;
		}


		if(World->MyRank==0)
		{
			if(iteration==1)
			{
				pnpPrintGroup0("PNPSR =========================================================================\n");
				pnpPrintGroup0("PNPSR  %5s %16s %22s %12s %12s\n","Iter.", "MaxConcChange, M", "Energy,kT","dE","rel.E");
				pnpPrintGroup0("PNPSR -------------------------------------------------------------------------\n");
			}
			if(iteration!=1 && iteration%ConvergenceCheck==1)
			{
				pnpPrintGroup0("PNPSR -------------------------------------------------------------------------\n");
				pnpPrintGroup0("PNPSR  %5s %16s %22s %12s %12s\n","Iter.", "MaxConcChange, M", "Energy,kT","dE","rel.E");
				pnpPrintGroup0("PNPSR -------------------------------------------------------------------------\n");
			}
			pnpPrintGroup0("PNPSR  %5d %16.4e %22.14e %12.4e %12.4e\n", iteration, NernstPlank->MaxChange, Poisson->totalEnergy, Poisson->totalChange, Poisson->relativeChange);
			//pnpPrint("<PNPiteration Nit=\"%d\" MaxChange=\"%16.8g\" E=\"%16.8lg\" dE=\"%16.8lg\" rel_dE=\"%16.8lg\"/>\n", iteration, NernstPlank->MaxChange, Poisson->totalEnergy, Poisson->totalChange, Poisson->relativeChange);
		}
		if(iteration%ConvergenceCheck==0)
		{
			CombineAndPrintCurrents(iteration);
		}
	}
	if(iteration%ConvergenceCheck!=0)
		CombineAndPrintCurrents(iteration);

	pnpPrintGroup0("PNPSR -------------------------------------------------------------------------\n");
	pnpPrintGroup0("PNPSR  Poisson-Nernst-Planck solver reached max iterations.\n");
	pnpPrintGroup0("PNPSR =========================================================================\n");

	return EXIT_SUCCESS;
}
VectorField3D* PoissonNernstPlanckSolver::CalcCartI(int ion)
{
	pnpPrint("PoissonNernstPlanckSolver::CalcCartI ion=%d\n",ion);
	int * GS;
	int GS_X,GS_Y,GS_Z;
	int GS_XY;
	int GS_XYZ;
	int gridPoint;
	float gridScale;
	float fpoh;
	float flcon;
	float ftot;
	
	int i,j,k;
	int gridp,gridn;
	float pot;
	int gridInc;
	float *C[2];
	float * potential;
	float *UTMP=NULL;
	float *TMP=NULL;
	float **D=World->D;
	int kgrid,jgrid;
	
	NodeIndex* NIndex=World->NIndexing->NIndex;
	NodeIndexing* NIndexing=World->NIndexing;
	NodeIndex DiffBoarderMask=NodeIndexing::DiffIon0BoarderMask;
	NodeIndex BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;
	
	GS = World->GridSize;
	GS_X = GS[0];
	GS_Y = GS[1];
	GS_Z = GS[2];
	GS_XY = GS[0]*GS[1];
	GS_XYZ = GS[2]*GS_XY;
	gridScale = World->GridScale;
	
	C[0] = World->C[0];
	C[1] = World->C[1];
	potential = World->Potential;
	
	if(NernstPlank->UTMPSingle==NULL)
		NernstPlank->UTMPSingle=new float[GS_XYZ];
	if(NernstPlank->TMPSingle==NULL)
		NernstPlank->TMPSingle=new float[GS_XYZ];
	
	UTMP=NernstPlank->UTMPSingle;
	TMP=NernstPlank->TMPSingle;
	
	VectorField3D* J=new VectorField3D(GS,gridScale,3);
	J->FillValue(0.0);
	float **V=J->V;
	
	if(World->PMF!=NULL)
		for(i=0;i<GS_XYZ;i++)
			UTMP[i]=0.5*(World->IonsQ[ion]*potential[i]+World->PMF[ion][i]);
	else
		for(i=0;i<GS_XYZ;i++)
			UTMP[i]=0.5*World->IonsQ[ion]*potential[i];
	
	fpoh = 4*M_PI*gridScale;
	flcon = 1.6*gridScale*gridScale*1000;
	ftot = flcon/fpoh;
	
	float d,dc,sc;
	double dtmp;
	//Calculate horizontal flew
	int ax;
	for(ax=0;ax<3;ax++)
	{
		switch(ax)
		{
			case 0:
				gridInc=1;
				break;
			case 1:
				gridInc=GS_X;
				break;
			case 2:
				gridInc=GS_XY;
				break;
		}
		ftot = -1.0*World->IonsQ[ion]*flcon/fpoh;
		if(World->D!=NULL)
		{
			for(k=0;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=0;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=0;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+gridInc;
						if(NIndexing->GetDiffFloat(ion,gridPoint)>0.0)
						{
							pot = UTMP[gridp]-UTMP[gridPoint];
							d = NIndexing->GetDiffFloat(ion,gridp)>0.0?					0.5*(D[ion][gridPoint]+D[ion][gridp]):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							V[ax][gridPoint]=ftot*d*(dc+pot*sc);
						}
					}
				}
			}
		}
		else
		{
			for(k=1;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				//dtmp=0.0;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+gridInc;
						if(NIndexing->GetDiffFloat(ion,gridPoint)>0.0)
						{
							pot = UTMP[gridp]-UTMP[gridPoint];
							d = NIndexing->GetDiffFloat(ion,gridp)>0.0?					0.5*(NIndexing->GetDiffFloat(ion,gridPoint)+NIndexing->GetDiffFloat(ion,gridp)):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							V[ax][gridPoint]=ftot*d*(dc+pot*sc);
						}
					}
				}
			}
		}
	}
	return J;
}
VectorField3D* PoissonNernstPlanckSolver::CalcIinout(VectorField3D* CartI)
{
	int *GS = World->GridSize;
	float gridScale = World->GridScale;
	int GS_X = GS[0];
	int GS_Y = GS[1];
	int GS_Z = GS[2];
	int GS_XY = GS[0]*GS[1];
	int GS_XYZ = GS[2]*GS_XY;
	int ix, iy, iz, GrdPnt;
	VectorField3D* J=new VectorField3D(GS,gridScale,3);
	J->FillValue(0.0);
	for(iz=1;iz<GS_Z-1;iz++)
		for(iy=1;iy<GS_Y-1;iy++)
			for(ix=1;ix<GS_X-1;ix++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		J->V[0][GrdPnt]=0.5*(CartI->V[0][GrdPnt-1]+CartI->V[0][GrdPnt]);
		J->V[1][GrdPnt]=0.5*(CartI->V[1][GrdPnt-GS_X]+CartI->V[1][GrdPnt]);
		J->V[2][GrdPnt]=0.5*(CartI->V[2][GrdPnt-GS_XY]+CartI->V[2][GrdPnt]);
	}
	return J;
}
VectorField3D* PoissonNernstPlanckSolver::CalcAvrI(VectorField3D* I,int iavr)
{
	int *GS = World->GridSize;
	float gridScale = World->GridScale;
	int GS_X = GS[0];
	int GS_Y = GS[1];
	int GS_Z = GS[2];
	int GS_XY = GS[0]*GS[1];
	int GS_XYZ = GS[2]*GS_XY;
	int ix, iy, iz, GrdPnt;
	VectorField3D* J=new VectorField3D(GS,gridScale,3);
	J->FillValue(0.0);
	for(iz=1;iz<GS_Z-1;iz++)
		for(iy=1;iy<GS_Y-1;iy++)
			for(ix=1;ix<GS_X-1;ix++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		J->V[0][GrdPnt]=0.5*(I->V[0][GrdPnt-1]+I->V[0][GrdPnt]);
		J->V[1][GrdPnt]=0.5*(I->V[1][GrdPnt-GS_X]+I->V[1][GrdPnt]);
		J->V[2][GrdPnt]=0.5*(I->V[2][GrdPnt-GS_XY]+I->V[2][GrdPnt]);
	}
	return J;
}
int PoissonNernstPlanckSolver::SolveDouble()
{
	pnpPrint("\nPoisson-Nernst-Plank Solver Double\n");
	int iteration;
	float Tolerance;
	int i,j,k,i1,i2,gridPoint,GrdPnt,ion;
	float gridScale;
	int * gridSize;
	float *PotentialTMP, *Potential;
  //assert(pnpSolverData!=NULL);
#ifdef MPI_PARALLEL
	MPI::Status  status;
#endif
	
	gridScale = World->GridScale;
	gridSize = World->GridSize;
  //PotentialTMP = World->PotentialTMP;
	Potential = World->Potential;
	if(World->PotentialDouble==NULL)
	{
		World->PotentialDouble=new double[World->GS_XYZ];
		for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
			World->PotentialDouble[GrdPnt]=World->Potential[GrdPnt];
	}
	if(World->CDouble==NULL)
	{
		World->CDouble=new double*[World->NIonsTypes];
		for(ion=0;ion<World->NIonsTypes;ion++)
		{
			World->CDouble[ion]=new double[World->GS_XYZ];
			for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
				World->CDouble[ion][GrdPnt]=World->C[ion][GrdPnt];
		}
	}
	for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
		World->PotentialDouble[GrdPnt]=World->Potential[GrdPnt];
	for(ion=0;ion<World->NIonsTypes;ion++)
	{
		for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
			World->CDouble[ion][GrdPnt]=World->C[ion][GrdPnt];
	}
	
	for(iteration=1;iteration<=MaxIterations;iteration++) 
	{
		int PoissonStatus,NernstPlankStatus;

		assert(Poisson);
		PoissonStatus=Poisson->Solve();
		if(PoissonStatus!=EXIT_SUCCESS)
		{
			return EXIT_FAILURE;
		}
		for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
			World->PotentialDouble[GrdPnt]=World->Potential[GrdPnt];
		if(PMFWeightMode==1)
		{
			NernstPlank->PMFWeight=double(iteration)/double(MaxIterations);
			pnpPrint0("PMFWeight=%f\n",NernstPlank->PMFWeight);
		}
		NernstPlankStatus=NernstPlank->NernstPlanckSolverDDouble();
    //NernstPlankStatus=NernstPlank->NernstPlanckSolverN2tmp();
		if(Poisson->solver!=Poisson->ArrayDirect)
			Poisson->SetQmobFromConcentrationDouble();
		if(NernstPlankStatus!=EXIT_SUCCESS)
		{
			return EXIT_FAILURE;
		}


		if(World->MyRank==0)
			pnpPrint("<PNPiteration Nit=\"%d\" MaxChange=\"%16.8g\" E=\"%16.8lg\" dE=\"%16.8lg\" rel_dE=\"%16.8lg\"/>\n", iteration, NernstPlank->MaxChange, Poisson->totalEnergy, Poisson->totalChange, Poisson->relativeChange);

		if(iteration%ConvergenceCheck==0)
		{
			for(ion=0;ion<World->NIonsTypes;ion++)
			{
				for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
					World->C[ion][GrdPnt]=World->CDouble[ion][GrdPnt];
			}
			CombineAndPrintCurrentsDouble(iteration);
		}
	}
	for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
		World->Potential[GrdPnt]=World->PotentialDouble[GrdPnt];
	for(ion=0;ion<World->NIonsTypes;ion++)
	{
		for(GrdPnt=0;GrdPnt<World->GS_XYZ;GrdPnt++)
			World->C[ion][GrdPnt]=World->CDouble[ion][GrdPnt];
	}
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::CombineAndPrintCurrents(int iteration)
{
	int i,j,k,i1,i2,gridPoint;
#ifdef MPI_PARALLEL
	MPI::Status  status;
#endif
	int iCurrentDimension;
	iCurrentDimension = currentDimension;
	
	int iLimitCurrentCalcZ[2];
	int diLimitCurrentCalcZ;
	if(bLimitCurrentCalc)
	{
		for(i=0;i<2;i++)
		{
			iLimitCurrentCalcZ[i]=int(World->ConvFloatToGlobIntUnitsZ(LimitCurrentCalcZ[i])+0.5);
			if(iLimitCurrentCalcZ[i]<1)
				iLimitCurrentCalcZ[i]=1;
			if(iLimitCurrentCalcZ[i]>World->GridSizeGlobal[2]-2)
				iLimitCurrentCalcZ[i]=World->GridSizeGlobal[2]-2;
		}
		diLimitCurrentCalcZ=iLimitCurrentCalcZ[1]-iLimitCurrentCalcZ[0]+1;
	}
	
	NernstPlank->nernstPlanckSolverCurrentProfile(iCurrentDimension,positiveCurrentProfile,negativeCurrentProfile);
	float Phi0=0.0,Phi1=0.0;
#ifdef MPI_PARALLEL
	pnpsapp->MyComGroup.Barrier();
	if(World->MyRank==0)
	{
		for(i=1;i<World->NProcs;i++)
		{
			World->GetBorder(&i1,&i2,i);
			pnpsapp->MyComGroup.Recv(positiveCurrentProfile+i1, i2-i1, MPI::DOUBLE, i, 0,status);
			pnpsapp->MyComGroup.Recv(negativeCurrentProfile+i1, i2-i1, MPI::DOUBLE, i, 1,status);
		}
	}
	else
	{
		pnpsapp->MyComGroup.Send(positiveCurrentProfile, World->GridSize[iCurrentDimension]-1, MPI::DOUBLE, 0, 0);
		pnpsapp->MyComGroup.Send(negativeCurrentProfile, World->GridSize[iCurrentDimension]-1, MPI::DOUBLE, 0, 1);
	}
	Phi0=World->Potential[0];
	Phi1=World->Potential[World->GS_XYZ-1];
			
	pnpsapp->MyComGroup.Bcast(&Phi0, 1, MPI::DOUBLE, 0);
	pnpsapp->MyComGroup.Bcast(&Phi1, 1, MPI::DOUBLE, World->NProcs-1);
#else
	DbgPrint0("Phi[0]=%g Phi[%d]=%g MyRank=%d\n",World->Potential[0],World->GS_XYZ-1,World->Potential[World->GS_XYZ-1],World->MyRank);
	Phi0=World->Potential[0];
	Phi1=World->Potential[World->GS_XYZ-1];
#endif
			
	if(World->MyRank==0)
	{
		double currentSum,currentCoordinate;
		double maxCurrent = 0.0;
		double maxCurrentCoordinate = 0.0;
		double currentDev = 0.0;
		double current = 0.0;
		double currentDevPos = 0.0;
		double currentPos = 0.0;
		double currentDevNeg = 0.0;
		double currentNeg = 0.0;

		pnpPrintGroup0("PNPSR-Current =================================================================\n");
		pnpPrintGroup0("PNPSR-Current  iteration %5d\n",iteration);
		pnpPrintGroup0("PNPSR-Current -----------------------------------------------------------------\n");
		pnpPrintGroup0("PNPSR-Current  %8s %16s %16s %16s\n","z, Ang.", "I0, pA", "I1, pA", "Itot, pA");
		pnpPrintGroup0("PNPSR-Current -----------------------------------------------------------------\n");

		//pnpPrint("<PNPCurrent iter=\"%5d\" Format=\"%%6.3f %%16.6f %%16.6f %%16.6f\" Note=\"z Ipos Ineg Itot\">\n",iteration);
		for(k=0;k<World->GridSizeGlobal[iCurrentDimension]-1;k++)
		{
			currentSum = positiveCurrentProfile[k]+negativeCurrentProfile[k];
			pnpPrintGroup0("PNPSR-Current  %8.3f %16.6f %16.6f %16.6f\n", (k-World->GridSizeOriginal[2]/2+0.5)/World->GridScale, positiveCurrentProfile[k], negativeCurrentProfile[k], currentSum);
			//pnpPrint("%6.3f %16.6f %16.6f %16.6f\n",(k-World->GridSizeOriginal[2]/2+0.5)/World->GridScale, positiveCurrentProfile[k], negativeCurrentProfile[k], currentSum);
			currentDev += currentSum*currentSum;
			current += currentSum;
			currentDevPos += positiveCurrentProfile[k]*positiveCurrentProfile[k];
			currentPos += positiveCurrentProfile[k];
			currentDevNeg += negativeCurrentProfile[k]*negativeCurrentProfile[k];
			currentNeg += negativeCurrentProfile[k];
			if(fabs(maxCurrent)<fabs(currentSum)) {
				maxCurrent = currentSum;
				maxCurrentCoordinate = k;
			}
		}
		pnpPrintGroup0("PNPSR-Current -----------------------------------------------------------------\n");

		//pnpPrint("</PNPCurrent>\n");
		currentDev /= World->GridSizeGlobal[iCurrentDimension]-1;
		current /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentDevPos /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentPos /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentDevNeg /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentNeg /= World->GridSizeGlobal[iCurrentDimension]-1;
		ItotErr=currentDev-current*current;
		IposErr=currentDevPos-currentPos*currentPos;
		InegErr=currentDevNeg-currentNeg*currentNeg;
		DbgPrint0("ItotErr=%g IposErr=%g InegErr=%g\n",ItotErr, IposErr, InegErr);
		ItotErr=sqrt(ItotErr);
		IposErr=sqrt(IposErr);
		InegErr=sqrt(InegErr);
		pnpPrintGroup0("PNPSR-Current  maxCurrent = %f pA at z = %f\n",maxCurrent,maxCurrentCoordinate);
		pnpPrintGroup0("PNPSR-Current -----------------------------------------------------------------\n");
		pnpPrintGroup0("PNPSR-Current  Current ( Total ) = %14.7f +- %14.7f pA\n",current, ItotErr);
		pnpPrintGroup0("PNPSR-Current  Current ( Ion0  ) = %14.7f +- %14.7f pA\n",currentPos, IposErr);
		pnpPrintGroup0("PNPSR-Current  Current ( Ion1  ) = %14.7f +- %14.7f pA\n",currentNeg, InegErr);
		pnpPrintGroup0("PNPSR-Current  Applied potential = %14.7f mV\n",(Phi0-Phi1)*CONFAC);
		pnpPrintGroup0("PNPSR-Current -----------------------------------------------------------------\n");
		Itot=current;
		Ipos=currentPos;
		Ineg=currentNeg;
		if(bLimitCurrentCalc)
		{
			currentDev = 0.0;
			current = 0.0;
			currentDevPos = 0.0;
			currentPos = 0.0;
			currentDevNeg = 0.0;
			currentNeg = 0.0;
			for(k=iLimitCurrentCalcZ[0];k<=iLimitCurrentCalcZ[1];k++)
			{
				currentSum = positiveCurrentProfile[k]+negativeCurrentProfile[k];
				currentDev += currentSum*currentSum;
				current += currentSum;
				currentDevPos += positiveCurrentProfile[k]*positiveCurrentProfile[k];
				currentPos += positiveCurrentProfile[k];
				currentDevNeg += negativeCurrentProfile[k]*negativeCurrentProfile[k];
				currentNeg += negativeCurrentProfile[k];
			}
			currentDev /= diLimitCurrentCalcZ;
			current /= diLimitCurrentCalcZ;
			currentDevPos /= diLimitCurrentCalcZ;
			currentPos /= diLimitCurrentCalcZ;
			currentDevNeg /= diLimitCurrentCalcZ;
			currentNeg /= diLimitCurrentCalcZ;
			ItotErr=sqrt(currentDev-current*current);
			IposErr=sqrt(currentDevPos-currentPos*currentPos);
			InegErr=sqrt(currentDevNeg-currentNeg*currentNeg);
			pnpPrintGroup0("PNPSR-Current  Current along Z region: Z = [ %g, %g ]\n", LimitCurrentCalcZ[0], LimitCurrentCalcZ[1]);
			pnpPrintGroup0("PNPSR-Current  CurrentZlimit(Total) = %14.7f +- %14.7f pA\n",current, ItotErr);
			pnpPrintGroup0("PNPSR-Current  CurrentZlimit(Ion0 ) = %14.7f +- %14.7f pA\n",currentPos, IposErr);
			pnpPrintGroup0("PNPSR-Current  CurrentZlimit(Ion1 ) = %14.7f +- %14.7f pA\n",currentNeg, InegErr);
			pnpPrintGroup0("PNPSR-Current  Applied potential    = %14.7f mV\n",(Phi0-Phi1)*CONFAC);
			//pnpPrintGroup0("PNPSR-Current -----------------------------------------------------------------\n");
			//pnpPrint0("<PNPiterCurSummary Zlimit=\"%g %g\" Current=\"%.7e +- %.7e pA\", CurrentPos=\"%.7e +- %.7e pA\" CurrentNeg=\"%.7e +- %.7e pA\" dU=\"%.7e mV\"/>\n",LimitCurrentCalcZ[0],LimitCurrentCalcZ[1],current, ItotErr,currentPos, IposErr,currentNeg, InegErr,(Phi0-Phi1)*CONFAC);
			Itot=current;
			Ipos=currentPos;
			Ineg=currentNeg;
		}
		pnpPrintGroup0("PNPSR-Current =================================================================\n");
	}
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::CombineAndPrintCurrentsDouble(int iteration)
{
	int i,j,k,i1,i2,gridPoint;
#ifdef MPI_PARALLEL
	MPI::Status  status;
#endif
	int iCurrentDimension;
	iCurrentDimension = currentDimension;
	
	int iLimitCurrentCalcZ[2];
	int diLimitCurrentCalcZ;
	if(bLimitCurrentCalc)
	{
		for(i=0;i<2;i++)
		{
			iLimitCurrentCalcZ[i]=int(World->ConvFloatToGlobIntUnitsZ(LimitCurrentCalcZ[i])+0.5);
			if(iLimitCurrentCalcZ[i]<1)
				iLimitCurrentCalcZ[i]=1;
			if(iLimitCurrentCalcZ[i]>World->GridSizeGlobal[2]-2)
				iLimitCurrentCalcZ[i]=World->GridSizeGlobal[2]-2;
		}
		diLimitCurrentCalcZ=iLimitCurrentCalcZ[1]-iLimitCurrentCalcZ[0]+1;
	}
	
	NernstPlank->nernstPlanckSolverCurrentProfile(iCurrentDimension,positiveCurrentProfile,negativeCurrentProfile);
	float Phi0=0.0,Phi1=0.0;
#ifdef MPI_PARALLEL
	pnpsapp->MyComGroup.Barrier();
	if(World->MyRank==0)
	{
		for(i=1;i<World->NProcs;i++)
		{
			World->GetBorder(&i1,&i2,i);
			pnpsapp->MyComGroup.Recv(positiveCurrentProfile+i1, i2-i1, MPI::DOUBLE, i, 0,status);
			pnpsapp->MyComGroup.Recv(negativeCurrentProfile+i1, i2-i1, MPI::DOUBLE, i, 1,status);
		}
	}
	else
	{
		pnpsapp->MyComGroup.Send(positiveCurrentProfile, World->GridSize[iCurrentDimension]-1, MPI::DOUBLE, 0, 0);
		pnpsapp->MyComGroup.Send(negativeCurrentProfile, World->GridSize[iCurrentDimension]-1, MPI::DOUBLE, 0, 1);
	}
	Phi0=World->Potential[0];
	Phi1=World->Potential[World->GS_XYZ-1];
			
	pnpsapp->MyComGroup.Bcast(&Phi0, 1, MPI::DOUBLE, 0);
	pnpsapp->MyComGroup.Bcast(&Phi1, 1, MPI::DOUBLE, World->NProcs-1);
#else
	DbgPrint0("Phi[0]=%g Phi[%d]=%g MyRank=%d\n",World->Potential[0],World->GS_XYZ-1,World->Potential[World->GS_XYZ-1],World->MyRank);
	Phi0=World->Potential[0];
	Phi1=World->Potential[World->GS_XYZ-1];
#endif
			
	if(World->MyRank==0)
	{
		double currentSum,currentCoordinate;
		double maxCurrent = 0.0;
		double maxCurrentCoordinate = 0.0;
		double currentDev = 0.0;
		double current = 0.0;
		double currentDevPos = 0.0;
		double currentPos = 0.0;
		double currentDevNeg = 0.0;
		double currentNeg = 0.0;
		pnpPrint("<PNPCurrent iter=\"%5d\" Format=\"%%6.3f %%16.6f %%16.6f %%16.6f\" Note=\"z Ipos Ineg Itot\">\n",iteration);
		for(k=0;k<World->GridSizeGlobal[iCurrentDimension]-1;k++)
		{
			currentSum = positiveCurrentProfile[k]+negativeCurrentProfile[k];
			pnpPrint("%6.3f %16.6f %16.6f %16.6f\n",(k-World->GridSizeOriginal[2]/2+0.5)/World->GridScale, positiveCurrentProfile[k], negativeCurrentProfile[k], currentSum);
			currentDev += currentSum*currentSum;
			current += currentSum;
			currentDevPos += positiveCurrentProfile[k]*positiveCurrentProfile[k];
			currentPos += positiveCurrentProfile[k];
			currentDevNeg += negativeCurrentProfile[k]*negativeCurrentProfile[k];
			currentNeg += negativeCurrentProfile[k];
			if(fabs(maxCurrent)<fabs(currentSum)) {
				maxCurrent = currentSum;
				maxCurrentCoordinate = k;
			}
		}
		pnpPrint("</PNPCurrent>\n");
		currentDev /= World->GridSizeGlobal[iCurrentDimension]-1;
		current /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentDevPos /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentPos /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentDevNeg /= World->GridSizeGlobal[iCurrentDimension]-1;
		currentNeg /= World->GridSizeGlobal[iCurrentDimension]-1;
		ItotErr=currentDev-current*current;
		IposErr=currentDevPos-currentPos*currentPos;
		InegErr=currentDevNeg-currentNeg*currentNeg;
		DbgPrint0("ItotErr=%g IposErr=%g InegErr=%g\n",ItotErr, IposErr, InegErr);
		ItotErr=sqrt(ItotErr);
		IposErr=sqrt(IposErr);
		InegErr=sqrt(InegErr);
		pnpPrint0("<pnpSolver: maxCurrent = %f pA at %c = %f\n",maxCurrent,currentDimension,maxCurrentCoordinate);
		pnpPrint0("<PNPiterCurSummary Current=\"%.7e +- %.7e pA\", CurrentPos=\"%.7e +- %.7e pA\" CurrentNeg=\"%.7e +- %.7e pA\" dU=\"%.7e mV\"/>\n",current, ItotErr,currentPos, IposErr,currentNeg, InegErr,(Phi0-Phi1)*CONFAC);
		Itot=current;
		Ipos=currentPos;
		Ineg=currentNeg;
		if(bLimitCurrentCalc)
		{
			currentDev = 0.0;
			current = 0.0;
			currentDevPos = 0.0;
			currentPos = 0.0;
			currentDevNeg = 0.0;
			currentNeg = 0.0;
			for(k=iLimitCurrentCalcZ[0];k<=iLimitCurrentCalcZ[1];k++)
			{
				currentSum = positiveCurrentProfile[k]+negativeCurrentProfile[k];
				currentDev += currentSum*currentSum;
				current += currentSum;
				currentDevPos += positiveCurrentProfile[k]*positiveCurrentProfile[k];
				currentPos += positiveCurrentProfile[k];
				currentDevNeg += negativeCurrentProfile[k]*negativeCurrentProfile[k];
				currentNeg += negativeCurrentProfile[k];
			}
			currentDev /= diLimitCurrentCalcZ;
			current /= diLimitCurrentCalcZ;
			currentDevPos /= diLimitCurrentCalcZ;
			currentPos /= diLimitCurrentCalcZ;
			currentDevNeg /= diLimitCurrentCalcZ;
			currentNeg /= diLimitCurrentCalcZ;
			ItotErr=sqrt(currentDev-current*current);
			IposErr=sqrt(currentDevPos-currentPos*currentPos);
			InegErr=sqrt(currentDevNeg-currentNeg*currentNeg);
			pnpPrint0("<PNPiterCurSummary Zlimit=\"%g %g\" Current=\"%.7e +- %.7e pA\", CurrentPos=\"%.7e +- %.7e pA\" CurrentNeg=\"%.7e +- %.7e pA\" dU=\"%.7e mV\"/>\n",LimitCurrentCalcZ[0],LimitCurrentCalcZ[1],current, ItotErr,currentPos, IposErr,currentNeg, InegErr,(Phi0-Phi1)*CONFAC);
			Itot=current;
			Ipos=currentPos;
			Ineg=currentNeg;
		}
	}
	return EXIT_SUCCESS;
}
int PoissonNernstPlanckSolver::ShowParameters()
{
	cout << "\nParameters of Poisson-Nernst-Plank solver\n";
	cout << "    MaxIterations:.................... "<< MaxIterations <<"\n";
	return EXIT_SUCCESS;
}
