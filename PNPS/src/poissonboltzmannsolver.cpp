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

#include "poissonboltzmannsolver.h"
#include "tinyxml.h"
#include "pnpdebug.h"
#include "contworld.h"
#include "math.h"
#include "pnpconstants.h"
#include "mapio.h"

#include <stdlib.h>

PoissonBoltzmannSolver::PoissonBoltzmannSolver()
 : GenericSolver()
{
	InitZero();
}


PoissonBoltzmannSolver::~PoissonBoltzmannSolver()
{
	Clear();
}
int PoissonBoltzmannSolver::InitZero()
{
	PnpsObject::SetName("PoissonBoltzmannSolver");
	SolverStr.push_back("Auto");
	SolverStr.push_back("NodeIndexBased");
	SolverStr.push_back("ArrayDirect");
	
	m_ContWorld=NULL;
	
	ConvergenceCheck=20;
	MaxIterationsLPB=300;
	MaxIterationsNPB=300;
	Convergence=0.0;
	Relaxation=1.9;
	
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
	Qst=NULL;
	PhiCharge=NULL;
	PhiSingular=NULL;
	PBZoneNum[0]=0;
	PBZoneNum[1]=0;
	PBZoneNum[2]=0;
	IndexPBZone=NULL;
	om2d6LPB=0.0;
	potential=NULL;

	IonicStrength=0.0;
	return EXIT_SUCCESS;
}
int PoissonBoltzmannSolver::Clear()
{
//	DeleteCArray(IonsQ);
//	DeleteCVecArray(C,NIonsTypes);
//	DeleteObjByPnt(NIndexing);
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
	DeleteCArray(IndexPBZone);
	DeleteCArray(PhiCharge);
	DeleteCArray(PhiSingular);
	return EXIT_SUCCESS;
}
int PoissonBoltzmannSolver::SetRelaxation(float _Relaxation)
{
	int i,gridPoint;
	Relaxation=_Relaxation;
	//Scale Parameters to Internal Units
	//Calculate Derivative Parameters
	if(dielectricZSSUM!=NULL)
	{
		for(i=0;i<SingularNum[2];i++)
		{
			gridPoint=IndexSingular[i];
			dielectricZSSUM[i]=dielectricXS[i]+dielectricYS[i]+dielectricZS[i]+dielectricXmS[i]+dielectricYmS[i]+dielectricZmS[i];
			dielectricZSSUM[i]=Relaxation/dielectricZSSUM[i];
		}
	}
	if(dielectricZDBSUM!=NULL)
	{
		for(i=0;i<DielBoarderNum[2];i++)
		{
			gridPoint=IndexDielBoarder[i];
			dielectricZDBSUM[i]=dielectricXDB[i]+dielectricYDB[i]+dielectricZDB[i]+dielectricXmDB[i]+dielectricYmDB[i]+dielectricZmDB[i];
			dielectricZDBSUM[i]=Relaxation/dielectricZDBSUM[i];
		}
	}
	om2 = Relaxation;
	om1 = 1.0-om2;
	om2d6 = om2/6.0;
	
	//PB
	if(m_ContWorld!=NULL)
	{
		float A = 2.0 * IonicStrength/m_ContWorld->NIndexing->Eps[ (m_ContWorld->NIndexing->NIndex[0]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft ];
		om2d6LPB=om2/(6.0+A);
	}
	return EXIT_SUCCESS;
}
int PoissonBoltzmannSolver::SetContWorld(ContWorld* _world)
{
	int i;
	m_ContWorld=_world;
	GridScale = m_ContWorld->GridScale;
	GS_X = m_ContWorld->GridSize[0];
	GS_Y = m_ContWorld->GridSize[1];
	GS_Z = m_ContWorld->GridSize[2];
	GS_XY = GS_X*GS_Y;
	GS_XYZ = GS_XY*GS_Z;
	if(m_ContWorld->Potential == NULL){
		if(!(m_ContWorld->Potential = new float[GS_XYZ])){
			fprintf(stderr,"ERROR 104: No memory available\n");
			exit(104);
		}
		for(i=0;i<GS_XYZ;i++)m_ContWorld->Potential[i]=0.0;
	}
	potential = m_ContWorld->Potential;
	float I=0.0;
	if(m_ContWorld->NIonsTypes==0)
	{
		I=0.0;
	}
	else
	{
		for(i=1;i<GS_XYZ;i++)
		{
			if(m_ContWorld->NIndexing->GetConcFloat(0,i)!=0.0){
				I=m_ContWorld->NIndexing->GetConcFloat(0,i);
				break;
			}
		}
	}
	IonicStrength=I;
	/*
	for(i=0;i<NodeIndexMaxValues;i++)
		if(m_ContWorld->NIndexing->C[i]>0.0)
	{
		I=m_ContWorld->NIndexing->C[i];
		break;
	}*/
	
	float A = 2.0 * IonicStrength/m_ContWorld->NIndexing->Eps[ (m_ContWorld->NIndexing->NIndex[0]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft ];
	om2d6LPB=om2/(6.0+A);
	return EXIT_SUCCESS;
}
int PoissonBoltzmannSolver::ShowParameters()
{
	DbgPrint2("PoissonBoltzmannSolver::ShowParameters\n");
	
	pnpPrintGroup0("\nParameters of Poisson solver set up\n");
	pnpPrintGroup0("    MaxIterationsLPB:................. %d\n", MaxIterationsLPB);
	pnpPrintGroup0("    MaxIterationsNPB:................. %d\n", MaxIterationsNPB);
	pnpPrintGroup0("    Convergence:...................... %.8g kT\n", Convergence);
	pnpPrintGroup0("    Relaxation:....................... %.5g\n", Relaxation);
	
	return EXIT_SUCCESS;
}
int PoissonBoltzmannSolver::ShowProperties()
{
	DbgPrint2("PoissonBoltzmannSolver::ShowProperties\n");
	return EXIT_SUCCESS;
}
int PoissonBoltzmannSolver::InitSolver()
{
	int status;
	bool bGuessNumberOfIteration=false;
	if(Relaxation<0.0||MaxIterationsLPB<0)
	{
		bGuessNumberOfIteration=true;
		Relaxation=1.0;
	}
  //Solver{Auto=0,NodeIndexBased=1,ArrayDirect=2,PNPC=3};
  //Solver solver;
	//if(solver==Auto)
	//{
	//	if(m_ContWorld->NIndexing==NULL)return InitSolverAD();
	//	else return InitSolverNIB();
	//}
	//else if(solver==NodeIndexBased)
	status=InitSolverNIB();
	//else if(solver==ArrayDirect)
	//	InitSolverAD();
	int itmp;
	if(bGuessNumberOfIteration)
		itmp=GuessNumberOfIteration();
	if(MaxIterationsLPB<0)
	{
		MaxIterationsLPB=itmp*abs(MaxIterationsLPB);
		if(MaxIterationsLPB<60)
			MaxIterationsLPB=60;
	}
	if(MaxIterationsNPB<0)
	{
		MaxIterationsNPB=itmp*abs(MaxIterationsNPB);
		if(MaxIterationsNPB<60)
			MaxIterationsNPB=60;
	}
	return status;
}
int PoissonBoltzmannSolver::InitSolverNIB()
{
	DbgPrint2("PoissonBoltzmannSolver::InitSolverNIB\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld,"m_ContWorld is not initialized\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld->Potential,"Potential is not initialize at m_ContWorld\n");
	
// 	if(m_ContWorld->C==NULL)
// 	{
// 		pnpWarning("Explicit concentrations of mobile ions are not present will try to get from NIndexing\n");
// 		PNP_EXIT_FAIL_NULL(m_ContWorld->NIndexing,"NIndexing is not initialize at m_ContWorld and thus cannot get concentrations for mobile ions\n");
// 		PNP_EXIT_ON_FAIL(m_ContWorld->SetInitConcentrationFromNIndexing());
// 		pnpPrint("Concentrations of mobile ions initiated from NIndexing and value of C[0]=%f\n",m_ContWorld->C[0][0]);
// 	}
	//temp vars
	int i,j,k,kgrid,jgrid,BlackOrWhite;
	int GrdPnt;
	NodeIndexing* NIndexing=m_ContWorld->NIndexing;
	NodeIndex* NIndex=m_ContWorld->NIndexing->NIndex;
  
//	Clear();
	
	//count b/w regions
	NoSingularNum[0]=0;NoSingularNum[1]=0;NoSingularNum[2]=0;
	SingularNum[0]=0;SingularNum[1]=0;SingularNum[2]=0;
	DielBoarderNum[0]=0;DielBoarderNum[1]=0;DielBoarderNum[2]=0;
	ChargeNum[0]=0;ChargeNum[1]=0;ChargeNum[2]=0;
	PBZoneNum[0]=0;PBZoneNum[1]=0;PBZoneNum[2]=0;
	
	//correction for spesific type of calculation, P(no additional charge) P(NP)(dynamic charge)
	unsigned int specChargeMask=NodeIndexing::ChargeMask;
	//!@todo not forget to change when will do pnp
//   if(m_ContWorld->D!=NULL)//i.e. P(NP)
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
	

	int CountQmobOnDielBoarder=0;
	int CountQmobOnQst=0;
	int CountQmobOnDielBoarderAndQst=0;
	for(k=1;k<GS_Z-1;k++)
	{
		kgrid = k*GS_XY;
		for(j=1;j<GS_Y-1;j++)
		{
			jgrid = kgrid+j*GS_X;
			for(i=1;i<GS_X-1;i++)
			{
				GrdPnt = jgrid+i;
				if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
				{
					
					SingularNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
					SingularNum[2]++;
				}
				else if(NIndex[GrdPnt]&specChargeMask)
				{
					ChargeNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
					ChargeNum[2]++;
				}
				else if(NIndex[GrdPnt]&DielBoarderMask)
				{
					//DbgPrint0("NIndex[%d]=%X\n",GrdPnt,NIndex[GrdPnt]);
					unsigned int t1=NIndex[GrdPnt];
					DielBoarderNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
					DielBoarderNum[2]++;
				}
				//PB
				else if(NIndexing->GetConcFloat(0,GrdPnt)>0.0)
				{
					PBZoneNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
					PBZoneNum[2]++;
				}
				else
				{
					NoSingularNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
					NoSingularNum[2]++;
				}
			}
		}
	}
	DbgPrint0("CountQmobOnDielBoarder=%d CountQmobOnQst=%d CountQmobOnDielBoarderAndQst=%d\n",CountQmobOnDielBoarder, CountQmobOnQst, CountQmobOnDielBoarderAndQst);
	if(CountQmobOnQst>0)
	{
		fprintf(stderr,"ERROR\nERROR situation then Qmob On Qst is not implemented yet\n");
	}
		
	DbgPrint1("InitSolver:Total: MyrankWorld->MyRank=%d icharge=% d    \niDielBoarder=%d iSingular=%d iNoSingular=%d PBZoneNum=%d\n",
						m_ContWorld->MyRank,ChargeNum[2],DielBoarderNum[2],SingularNum[2],NoSingularNum[2],PBZoneNum[2]);
	DbgPrint0("TotalNodes: %d\n",ChargeNum[2]+DielBoarderNum[2]+SingularNum[2]+NoSingularNum[2]+PBZoneNum[2]);
	DbgPrint1("InitSolver:Blach cell:: MyrankWorld->MyRank=%d icharge=% d    \niDielBoarder=%d iSingular=%d iNoSingular=%d PBZoneNum=%d\n",
						m_ContWorld->MyRank,ChargeNum[1],DielBoarderNum[1],SingularNum[1],NoSingularNum[1],PBZoneNum[2]);
  
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
	if(!(Qst = new float[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(dielectricCh = new float[ChargeNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(IndexSingular = new int[SingularNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!(QstS = new float[SingularNum[2]])){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	//PB
	if(!(IndexPBZone = new int[PBZoneNum[2]])){
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
	
  //fill indexes of nodes types
	int iCharge=0,iDielBoarder=0,iSingular=0,iNoSingular=0;
	int iCharge2=ChargeNum[1],iDielBoarder2=DielBoarderNum[1];
	int iSingular2=SingularNum[1],iNoSingular2=NoSingularNum[1];
	//PB
	int iPBZone=0,iPBZone2=PBZoneNum[1];
	
	int QCount=0;
	float *Q=m_ContWorld->NIndexing->Q;
	float *Eps=m_ContWorld->NIndexing->Eps;
	double q=0.0f;
	for(i=0;i<m_ContWorld->NIndexing->QNum;i++)
	{
		q+=Q[i];
	}
	DbgPrint0("q=%f QCount=%d GridScale=%f\n",(float)q/4/M_PI/m_ContWorld->GridScale,m_ContWorld->NIndexing->QNum,m_ContWorld->GridScale);
	q=0.0f;
	
	for(k=1;k<GS_Z-1;k++)
	{
		kgrid = k*GS_XY;
		for(j=1;j<GS_Y-1;j++)
		{
			jgrid = kgrid+j*GS_X;
			for(i=1;i<GS_X-1;i++)
			{
				GrdPnt = jgrid+i;
				
				if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
				{
					if(NIndex[GrdPnt]&BlackAndWhiteMask)
					{
						IndexSingular[iSingular]=GrdPnt;
						QstS[iSingular]=Q[QCount];
						q+=Q[QCount];
						QCount++;
						iSingular++;
					}
					else
					{
						IndexSingular[iSingular2]=GrdPnt;
						QstS[iSingular2]=Q[QCount];
						q+=Q[QCount];
						QCount++;
						iSingular2++;
					}
				}
				else if(NIndex[GrdPnt]&specChargeMask)
				{
					if(NIndex[GrdPnt]&BlackAndWhiteMask)
					{
						IndexCharge[iCharge]=GrdPnt;
						dielectricCh[iCharge]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
						Qst[iCharge]=Q[QCount]/dielectricCh[iCharge];
						q+=Q[QCount];
						QCount++;
						iCharge++;
					}
					else
					{
						IndexCharge[iCharge2]=GrdPnt;
						dielectricCh[iCharge2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
						Qst[iCharge2]=Q[QCount]/dielectricCh[iCharge2];
						q+=Q[QCount];
						QCount++;
						iCharge2++;
					}
				}
				else if(NIndex[GrdPnt]&DielBoarderMask)
				{
					if(NIndex[GrdPnt]&BlackAndWhiteMask)
					{
						IndexDielBoarder[iDielBoarder]=GrdPnt;
						iDielBoarder++;
					}
					else
					{
						IndexDielBoarder[iDielBoarder2]=GrdPnt;
						iDielBoarder2++;
					}
				}
				//PB
				else if(NIndexing->GetConcFloat(0,GrdPnt)>0.0)
				{
					if(NIndex[GrdPnt]&BlackAndWhiteMask)
					{
						IndexPBZone[iPBZone]=GrdPnt;
						iPBZone++;
					}
					else
					{
						IndexPBZone[iPBZone2]=GrdPnt;
						iPBZone2++;
					}
				}
				else
				{
					if(NIndex[GrdPnt]&BlackAndWhiteMask)
					{
						IndexNoSingular[iNoSingular]=GrdPnt;
						iNoSingular++;
					}
					else
					{
						IndexNoSingular[iNoSingular2]=GrdPnt;
						iNoSingular2++;
					}
				}
				//if(NIndex[GrdPnt]&NodeIndexing::ChargeMask)QCount++;
			}
		}
	}
	DbgPrint0("q=%f QCount=%d",(float)q/4/M_PI/m_ContWorld->GridScale,QCount);
  //Allocate dielectric helping array
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
  

  
	if(!dielectricXDB||!dielectricYDB||!dielectricZDB||!dielectricZDBSUM||!dielectricXS||!dielectricYS||!dielectricZS||!dielectricZSSUM){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
	if(!dielectricXmDB||!dielectricYmDB||!dielectricZmDB||!dielectricXmS||!dielectricYmS||!dielectricZmS){
		fprintf(stderr,"ERROR 204: No memory available\n");
		exit(204);
	}
  
  //float *Epsilon[3];
  //for(i=0;i<3;i++)Epsilon[i] = m_ContWorld->Epsilon[i];
  
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		dielectricXS[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		dielectricYS[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		dielectricZS[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		dielectricXmS[i]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		dielectricYmS[i]=Eps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		dielectricZmS[i]=Eps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		dielectricZSSUM[i]=dielectricXS[i]+dielectricYS[i]+dielectricZS[i]+dielectricXmS[i]+dielectricYmS[i]+dielectricZmS[i];
		dielectricZSSUM[i]=Relaxation/dielectricZSSUM[i];
	}
	for(i=0;i<DielBoarderNum[2];i++)
	{
		GrdPnt=IndexDielBoarder[i];
		dielectricXDB[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		dielectricYDB[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		dielectricZDB[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		dielectricXmDB[i]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		dielectricYmDB[i]=Eps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		dielectricZmDB[i]=Eps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		dielectricZDBSUM[i]=dielectricXDB[i]+dielectricYDB[i]+dielectricZDB[i]+dielectricXmDB[i]+dielectricYmDB[i]+dielectricZmDB[i];
		dielectricZDBSUM[i]=Relaxation/dielectricZDBSUM[i];
	}
	//PB
	float A = 2.0 * IonicStrength/m_ContWorld->NIndexing->Eps[ (m_ContWorld->NIndexing->NIndex[0]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft ];
	om2d6LPB=om2/(6.0+A);
	return EXIT_SUCCESS;
}

// int PoissonBoltzmannSolver::InitSolverAD()
// {
//       /// @todo implement me
//   //std::istrstream iStr(ProcedureCommand);
//   //char tmp[MAP_IO_STRING_LENGTH];
//   
//   //Read parameters
// 	DbgPrint2("Renew Singularities List\n");
// 
// 	int gridSizeX;
// 	int gridSizeY;
// 	int gridSizeZ;
// 	int gridSizeXY,gridSizeXYZ;
// 	int gridPoint;
// 	int i,j,k;
// 	float * dielectric[3];
// 	float * staticCharge;
// 	int jgrid,kgrid;
// 	int iCharge,iDielBoarder,iSingular,iNoSingular;
// 	int BlackOrWhite;
// 	float *positiveCharge;
// 	float *negativeCharge;
//   
// 	int *typeOfBorderPoint=NULL;
// 	
// 	m_ContWorld->CheckArrays("P",true);
// 	
// 	if(m_ContWorld->C==NULL)
// 	{
// 		positiveCharge = NULL;
// 		negativeCharge = NULL;
// 	}
// 	else
// 	{
// 		positiveCharge = m_ContWorld->C[0];
// 		negativeCharge = m_ContWorld->C[1];
// 	}
// 	gridSizeX = m_ContWorld->GridSize[0];
// 	gridSizeY = m_ContWorld->GridSize[1];
// 	gridSizeZ = m_ContWorld->GridSize[2];
// 	gridSizeXY = gridSizeX*gridSizeY;
// 	gridSizeXYZ = gridSizeZ*gridSizeXY;
// 
// 	if(dielectricXS!=NULL){delete [] dielectricXS;dielectricXS=NULL;}
// 	if(dielectricYS!=NULL){delete [] dielectricYS;dielectricYS=NULL;}
// 	if(dielectricZS!=NULL){delete [] dielectricZS;dielectricZS=NULL;}
// 	if(dielectricZSSUM!=NULL){delete [] dielectricZSSUM;dielectricZSSUM=NULL;}
// 	if(dielectricXmS!=NULL){delete [] dielectricXmS;dielectricXmS=NULL;}
// 	if(dielectricYmS!=NULL){delete [] dielectricYmS;dielectricYmS=NULL;}
// 	if(dielectricZmS!=NULL){delete [] dielectricZmS;dielectricZmS=NULL;}
// 	if(dielectricXDB!=NULL){delete [] dielectricXDB;dielectricXDB=NULL;}
// 	if(dielectricYDB!=NULL){delete [] dielectricYDB;dielectricYDB=NULL;}
// 	if(dielectricZDB!=NULL){delete [] dielectricZDB;dielectricZDB=NULL;}
// 	if(dielectricZDBSUM!=NULL){delete [] dielectricZDBSUM;dielectricZDBSUM=NULL;}
// 	if(dielectricXmDB!=NULL){delete [] dielectricXmDB;dielectricXmDB=NULL;}
// 	if(dielectricYmDB!=NULL){delete [] dielectricYmDB;dielectricYmDB=NULL;}
// 	if(dielectricZmDB!=NULL){delete [] dielectricZmDB;dielectricZmDB=NULL;}
// 	if(dielectricCh!=NULL){delete [] dielectricCh;dielectricCh=NULL;}
// 	if(IndexNoSingular!=NULL){delete [] IndexNoSingular;IndexNoSingular=NULL;}
// 	if(IndexDielBoarder!=NULL){delete [] IndexDielBoarder;IndexDielBoarder=NULL;}
// 	if(IndexCharge!=NULL){delete [] IndexCharge;IndexCharge=NULL;}
// 	if(IndexSingular!=NULL){delete [] IndexSingular;IndexSingular=NULL;}
// 	if(ChargeSum!=NULL){delete [] ChargeSum;ChargeSum=NULL;}
//   
// 
//   //Allocate Potential
// 	if(m_ContWorld->Potential == NULL){
// 		if(!(m_ContWorld->Potential = new float[gridSizeXYZ])){
// 			fprintf(stderr,"ERROR 104: No memory available\n");
// 			exit(104);
// 		}
// 		for(j=0;j<GS_XYZ;j++)m_ContWorld->Potential[j]=0.0;
// 	}
//   
// 
// 	for(i=0;i<3;i++)
// 		dielectric[i] = m_ContWorld->Epsilon[i];
//   
// 	if(!(typeOfBorderPoint = new int[GS_XYZ])){
// 		fprintf(stderr,"ERROR 204: No memory available\n");
// 		exit(204);
// 	}
// 
// 
// 	if(!(ChargeSum = new float[GS_XYZ])){
// 		fprintf(stderr,"ERROR 204: No memory available\n");
// 		exit(204);
// 	} 
// 	staticCharge=ChargeSum;
//   
// 	if(QmobMod==0||m_ContWorld->D==NULL)
// 	{
// 		for(i=0;i<GS_XYZ;i++)
// 			if(m_ContWorld->Qstat[i]!=0)staticCharge[i]=1.0;
// 	}
// 	else {
// 		for(i=0;i<GS_XYZ;i++)
// 			if(m_ContWorld->Qstat[i]!=0||m_ContWorld->D[0][i]>0||m_ContWorld->D[1][i]>0)staticCharge[i]=1.0;
// 	}
// 	for(k=0;k<GS_XYZ;k++){
//     //dielectricSum[k] = 6*dielectric[0][k];
// 		typeOfBorderPoint[k] = Boarder;
// 	}
// 
//   
// 	NoSingularNum[0]=0;NoSingularNum[1]=0;NoSingularNum[2]=0;
// 	SingularNum[0]=0;SingularNum[1]=0;SingularNum[2]=0;
// 	DielBoarderNum[0]=0;DielBoarderNum[1]=0;DielBoarderNum[2]=0;
// 	ChargeNum[0]=0;ChargeNum[1]=0;ChargeNum[2]=0;
//   
//   //int *Temp=new int[gridSizeXYZ];
// 	for(k=1;k<gridSizeZ-1;k++) {
// 		kgrid = k*gridSizeXY;
// 		for(j=1;j<gridSizeY-1;j++) {
// 			jgrid = kgrid+j*gridSizeX;
// 			for(i=1;i<gridSizeX-1;i++) {
// 				gridPoint = jgrid+i;
// 				BlackOrWhite=k+j+i+m_ContWorld->startBlackAndWhite;
// 				typeOfBorderPoint[gridPoint] = NoSingular;
//         
//         //if(BlackOrWhite%2==0)Temp[gridPoint]=0;
//         //else Temp[gridPoint]=gridPoint;
// 				if(dielectric[0][gridPoint]!=dielectric[0][gridPoint-1] ||
// 							 dielectric[0][gridPoint]!=dielectric[1][gridPoint]||
// 							 dielectric[0][gridPoint]!=dielectric[1][gridPoint-gridSizeX] ||
// 							 dielectric[0][gridPoint]!=dielectric[2][gridPoint] ||
// 							 dielectric[0][gridPoint]!=dielectric[2][gridPoint-gridSizeXY])
// 				{
//           //dielectricSum[gridPoint] = dielectric[0][gridPoint] + dielectric[0][gridPoint-1] + dielectric[1][gridPoint] + dielectric[1][gridPoint-gridSizeX] + dielectric[2][gridPoint] + dielectric[2][gridPoint-gridSizeXY];
// 					if(staticCharge[gridPoint]!=0){
// 						typeOfBorderPoint[gridPoint] = ChargeAndDielBoarder;
// 						if(BlackOrWhite%2==0)SingularNum[1]++;
// 						SingularNum[2]++;
// 					}
// 					else{
// 						typeOfBorderPoint[gridPoint] = DielBoarder;
// 						if(BlackOrWhite%2==0)DielBoarderNum[1]++;
// 						DielBoarderNum[2]++;
// 					}
// 				}
// 				else if(staticCharge[gridPoint]!=0) {
// 					typeOfBorderPoint[gridPoint] = Charge;
// 					if(BlackOrWhite%2==0)ChargeNum[1]++;
// 					ChargeNum[2]++;
// 				}
// 				if(typeOfBorderPoint[gridPoint] == NoSingular){
// 					if(BlackOrWhite%2==0)NoSingularNum[1]++;
// 					NoSingularNum[2]++;
// 				}
// 			}
// 		}
// 	}
//   //if(m_ContWorld->MyRank==0)m_ContWorld->writeMapFloat("tmpde0.gz",m_ContWorld->Epsilon[1],gridSizeXYZ,1);
//   
// 	DbgPrint1("Total: MyrankWorld->MyRank=%d icharge=% d    \niDielBoarder=%d iSingular=%d iNoSingular=%d\n",
// 						m_ContWorld->MyRank,ChargeNum[2],DielBoarderNum[2],SingularNum[2],NoSingularNum[2]);
// 	DbgPrint1("Blach cell:: MyrankWorld->MyRank=%d icharge=% d    \niDielBoarder=%d iSingular=%d iNoSingular=%d\n",
// 						m_ContWorld->MyRank,ChargeNum[1],DielBoarderNum[1],SingularNum[1],NoSingularNum[1]);
//     
// 	if(!(IndexNoSingular = new int[NoSingularNum[2]])){
// 		fprintf(stderr,"ERROR 204: No memory available\n");
// 		exit(204);
// 	}
// 	if(!(IndexDielBoarder = new int[DielBoarderNum[2]])){
// 		fprintf(stderr,"ERROR 204: No memory available\n");
// 		exit(204);
// 	}
// 	if(!(IndexCharge = new int[ChargeNum[2]])){
// 		fprintf(stderr,"ERROR 204: No memory available\n");
// 		exit(204);
// 	}
// 	if(!(IndexSingular = new int[SingularNum[2]])){
// 		fprintf(stderr,"ERROR 204: No memory available\n");
// 		exit(204);
// 	}
// 	if(!(PhiSingular = new float[SingularNum[2]])){
// 		fprintf(stderr,"ERROR 204: No memory available\n");
// 		exit(204);
// 	}
// 	if(!(PhiCharge = new float[ChargeNum[2]])){
// 		fprintf(stderr,"ERROR 204: No memory available\n");
// 		exit(204);
// 	}
// 	for(i=0;i<ChargeNum[2];i++)
// 		PhiCharge[i]=0.0f;
// 	for(i=0;i<SingularNum[2];i++)
// 		PhiSingular[i]=0.0f;
// 	
// 	iNoSingular=-1;
// 	iSingular=-1; 
// 	iCharge=-1;
// 	iDielBoarder=-1; 
//     //for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint=gridPoint+2){
// 	for(k=1;k<gridSizeZ-1;k++) {
// 		kgrid = k*gridSizeXY;
// 		for(j=1;j<gridSizeY-1;j++) {
// 			jgrid = kgrid+j*gridSizeX;
// 			for(i=1;i<gridSizeX-1;i++) {
// 				gridPoint = jgrid+i;
// 				BlackOrWhite=k+j+i+m_ContWorld->startBlackAndWhite;
// 				if(BlackOrWhite%2==0){
// 					switch(typeOfBorderPoint[gridPoint]){
// 						case NoSingular:
// 							iNoSingular++;
// 							IndexNoSingular[iNoSingular]=gridPoint;
// 							break;
// 						case Charge:
// 							iCharge++;
// 							IndexCharge[iCharge]=gridPoint;
// 							break;
// 						case DielBoarder:
// 							iDielBoarder++;
// 							IndexDielBoarder[iDielBoarder]=gridPoint;
// 							break;
// 						case ChargeAndDielBoarder:
// 							iSingular++;
// 							IndexSingular[iSingular]=gridPoint;
// 							break;
// 						case Boarder:
// 							break;   
// 					}
// 				}
// 			}
// 		}
// 	}
//     //for(gridPoint=1;gridPoint<gridSizeXYZ;gridPoint=gridPoint+2){
// 	for(k=1;k<gridSizeZ-1;k++) {
// 		kgrid = k*gridSizeXY;
// 		for(j=1;j<gridSizeY-1;j++) {
// 			jgrid = kgrid+j*gridSizeX;
// 			for(i=1;i<gridSizeX-1;i++) {
// 				gridPoint = jgrid+i;
// 				BlackOrWhite=k+j+i+m_ContWorld->startBlackAndWhite;
// 				if(BlackOrWhite%2==1){
// 					switch(typeOfBorderPoint[gridPoint]){
// 						case NoSingular:
// 							iNoSingular++;
// 							IndexNoSingular[iNoSingular]=gridPoint;
// 							break;
// 						case Charge:
// 							iCharge++;
// 							IndexCharge[iCharge]=gridPoint;
// 							break;
// 						case DielBoarder:
// 							iDielBoarder++;
// 							IndexDielBoarder[iDielBoarder]=gridPoint;
// 							break;
// 						case ChargeAndDielBoarder:
// 							iSingular++;
// 							IndexSingular[iSingular]=gridPoint;
// 							break;
// 						case Boarder:
// 							break;   
// 					}
// 				}
// 			}
// 		}
// 	}
//   //Free memory. SAVE it for children!
// 	if(typeOfBorderPoint!=NULL)
// 	{
// 		delete [] typeOfBorderPoint;
// 		typeOfBorderPoint=NULL;
// 	}
//   
// 	//
//   
// 	dielectricXDB = new float[DielBoarderNum[2]];
// 	dielectricYDB = new float[DielBoarderNum[2]];
// 	dielectricZDB = new float[DielBoarderNum[2]];
// 	dielectricXmDB = new float[DielBoarderNum[2]];
// 	dielectricYmDB = new float[DielBoarderNum[2]];
// 	dielectricZmDB = new float[DielBoarderNum[2]];
// 	dielectricZDBSUM = new float[DielBoarderNum[2]];
// 	dielectricXS = new float[SingularNum[2]];
// 	dielectricYS = new float[SingularNum[2]];
// 	dielectricZS = new float[SingularNum[2]];
// 	dielectricXmS = new float[SingularNum[2]];
// 	dielectricYmS = new float[SingularNum[2]];
// 	dielectricZmS = new float[SingularNum[2]];
// 	dielectricZSSUM = new float[SingularNum[2]];
//   
// 	if(!(dielectricCh = new float[ChargeNum[2]])){
// 		fprintf(stderr,"ERROR 204: No memory available\n");
// 		exit(204);
// 	}
//   
// 	if(!dielectricXDB||!dielectricYDB||!dielectricZDB||!dielectricZDBSUM||!dielectricXS||!dielectricYS||!dielectricZS||!dielectricZSSUM){
// 		fprintf(stderr,"ERROR 204: No memory available\n");
// 		exit(204);
// 	}
// 	if(!dielectricXmDB||!dielectricYmDB||!dielectricZmDB||!dielectricXmS||!dielectricYmS||!dielectricZmS){
// 		fprintf(stderr,"ERROR 204: No memory available\n");
// 		exit(204);
// 	}
//   
// 	for(i=0;i<SingularNum[2];i++){
// 		gridPoint=IndexSingular[i];
// 		dielectricXS[i]=dielectric[0][gridPoint];
// 		dielectricYS[i]=dielectric[1][gridPoint];
// 		dielectricZS[i]=dielectric[2][gridPoint];
// 		dielectricXmS[i]=dielectric[0][gridPoint-1];
// 		dielectricYmS[i]=dielectric[1][gridPoint-gridSizeX];
// 		dielectricZmS[i]=dielectric[2][gridPoint-gridSizeXY];
// 		dielectricZSSUM[i]=dielectricXS[i]+dielectricYS[i]+dielectricZS[i]+dielectricXmS[i]+dielectricYmS[i]+dielectricZmS[i];
// 		dielectricZSSUM[i]=Relaxation/dielectricZSSUM[i];
// 	}
// 	for(i=0;i<ChargeNum[2];i++){
// 		gridPoint=IndexCharge[i];
//     //denominator[gridPoint]=temp1;
// 		dielectricCh[i]=dielectric[0][gridPoint];
//     //staticCharge[gridPoint]/=dielectricCh[i];
// 	}
// 	for(i=0;i<DielBoarderNum[2];i++){
// 		gridPoint=IndexDielBoarder[i];
// 		dielectricXDB[i]=dielectric[0][gridPoint];
// 		dielectricYDB[i]=dielectric[1][gridPoint];
// 		dielectricZDB[i]=dielectric[2][gridPoint];
// 		dielectricXmDB[i]=dielectric[0][gridPoint-1];
// 		dielectricYmDB[i]=dielectric[1][gridPoint-gridSizeX];
// 		dielectricZmDB[i]=dielectric[2][gridPoint-gridSizeXY];
// 		dielectricZDBSUM[i]=dielectricXDB[i]+dielectricYDB[i]+dielectricZDB[i]+dielectricXmDB[i]+dielectricYmDB[i]+dielectricZmDB[i];
// 		dielectricZDBSUM[i]=Relaxation/dielectricZDBSUM[i];
// 	}
// //   if(MemoryCarefullUsage)
// //   {
// //     delete [] m_ContWorld->Epsilon[0];
// //     m_ContWorld->Epsilon[0]=NULL;
// //     delete [] m_ContWorld->Epsilon[1];
// //     m_ContWorld->Epsilon[1]=NULL;
// //     delete [] m_ContWorld->Epsilon[2];
// //     m_ContWorld->Epsilon[2]=NULL;
// //   }
// 	return EXIT_SUCCESS;
// }

int PoissonBoltzmannSolver::Solve()
{
	
	/*if(solver==Auto)
	{
		if(m_ContWorld->NIndexing==NULL)return PoissonBoltzmannSolverAD();
		else return PoissonBoltzmannSolverNIB();
	}
	else if(solver==NodeIndexBased)
		return PoissonBoltzmannSolverNIB();
	else if(solver==ArrayDirect)
		return PoissonBoltzmannSolverAD();*/
	
	PNP_EXIT_ON_FAIL_MES(LinearPBSolverNIB(MaxIterationsLPB),"Something was wrong during LinearPBSolverNIB\n");
	PNP_EXIT_ON_FAIL_MES(NonlinearPBSolverNIB(MaxIterationsNPB),"Something was wrong during NonlinearPBSolverNIB\n");
	return EXIT_SUCCESS;
}
// int PoissonBoltzmannSolver::PoissonBoltzmannSolverAD()
// {
// 	float gridScale;
// 	int gridSizeX;
// 	int gridSizeY;
// 	int gridSizeZ;
// 	int iteration;
// 	int i,j,k;
// 	int gridPoint;
// 	float om1,om2,om2d6;
// 	float * potential;
// 	float * staticCharge;
// 	int gridSizeXY;
// 	int gridSizeXYZ;
// 	float temp1,temp2,temp3,temp4,temp5,temp6,temp7;
// 	float dynamicChargeFactor,IonStrengthFactor;
// 	float fpoh;
// 	double totalEnergyOld=totalEnergy;
// 	bool *PeriodicBoundaryCondition;
// 	float *positiveCharge;
// 	float *negativeCharge;
// 	
// 	float * dielectric[3];
// 	
// 	
// 	for(i=0;i<3;i++)
// 		dielectric[i] = m_ContWorld->Epsilon[i];
// 	
// 	if(!(m_ContWorld->Potential)||!(m_ContWorld->Qstat)){
// 		fprintf(stderr,"ERROR 110: Arrays NOT yet initialize\n");
// 		exit(105);
// 	}
// 
// 	//cout<<"TRACE: CPoisson::poissonSolver() START\n";
// 	gridScale = m_ContWorld->GridScale;
// 	gridSizeX = m_ContWorld->GridSize[0];
// 	gridSizeY = m_ContWorld->GridSize[1];
// 	gridSizeZ = m_ContWorld->GridSize[2];
// 	gridSizeXY = gridSizeX*gridSizeY;
// 	gridSizeXYZ = gridSizeXY*gridSizeZ;
// 	
// 	PeriodicBoundaryCondition=m_ContWorld->PBC;
// 	
// 	potential = m_ContWorld->Potential;
// 	if(m_ContWorld->C==NULL)
// 	{
// 		positiveCharge = NULL;
// 		negativeCharge = NULL;
// 	}
// 	else
// 	{
// 		positiveCharge = m_ContWorld->C[0];
// 		negativeCharge = m_ContWorld->C[1];
// 	}
// 
// 	om2 = Relaxation;
// 	om1 = 1.0-om2;
// 	om2d6 = om2/6.0;
// 	
// 	//poissonBoltzmannSolver: assuming cation and anion concentration profiles are the same
// 	//poissonBoltzmannSolver: using cation concentration profile to calculate Debye lengths
// 	fpoh = 4.0*M_PI*gridScale;
// 	
// 	/*convert from charge density to M then convert from M to k'2
// 	* (or inverse debyle length squ					ared).
// 	*/
// 	dynamicChargeFactor = (float)1.0/(COANGS*4.0*M_PI*DFACT*DFACT);
// 	IonStrengthFactor = 1.0/(4.0*M_PI*gridScale);
// 	
// 	staticCharge=ChargeSum;
// 	
// 	//For speading calc denom first
// 	if(positiveCharge==0)
// 		for(k=0;k<gridSizeXYZ;k++){
// 		staticCharge[k]=m_ContWorld->Qstat[k];
// 		}
// 		else
// 			for(k=0;k<gridSizeXYZ;k++){
// 			staticCharge[k]=m_ContWorld->Qstat[k]+positiveCharge[k]-negativeCharge[k];
// 			}
// 			for(i=0;i<ChargeNum[2];i++){
// 				gridPoint=IndexCharge[i];
// 				staticCharge[gridPoint]/=dielectricCh[i];
// 			}
// 	
// 			for(iteration=1;iteration<=MaxIterations;iteration++) {
// 				for(j=0;j<=1;j++){
// 					for(i=SingularNum[j];i<SingularNum[j+1];i++){
// 						gridPoint=IndexSingular[i];
// 						temp1 = potential[gridPoint+1]*dielectricXS[i];
// 						temp2 = potential[gridPoint-1]*dielectricXmS[i];
// 						temp3 = potential[gridPoint+gridSizeX]*dielectricYS[i];
// 						temp4 = potential[gridPoint-gridSizeX]*dielectricYmS[i];
// 						temp5 = potential[gridPoint+gridSizeXY]*dielectricZS[i];
// 						temp6 = potential[gridPoint-gridSizeXY]*dielectricZmS[i];				
// 						temp1 = temp1+temp2;
// 						temp2 = temp3+temp4;
// 						temp3 = temp5+temp6;
// 						temp1 = temp1+temp2;
// 						temp2 = temp3+staticCharge[gridPoint];								
// 						temp1 = dielectricZSSUM[i]*(temp1+temp2);
// 						potential[gridPoint] = potential[gridPoint]*om1+temp1;
// 					}
// 					for(i=ChargeNum[j];i<ChargeNum[j+1];i++){
// 						gridPoint=IndexCharge[i];
// 						temp1 = potential[gridPoint+1]+potential[gridPoint-1];
// 						temp2 = potential[gridPoint+gridSizeX]+potential[gridPoint-gridSizeX];
// 						temp3 = potential[gridPoint+gridSizeXY]+potential[gridPoint-gridSizeXY];
// 						temp4 = potential[gridPoint]*om1;
// 						temp5 = om2d6*(temp1+temp2+temp3+staticCharge[gridPoint]);
// 				//temp6 = denominator[gridPoint]*(temp3+staticCharge[gridPoint]);
// 						potential[gridPoint] = temp4+temp5;
// 					}
// 					for(i=DielBoarderNum[j];i<DielBoarderNum[j+1];i++){
// 						gridPoint=IndexDielBoarder[i];
// 						temp1 = potential[gridPoint+1]*dielectricXDB[i];
// 						temp2 = potential[gridPoint-1]*dielectricXmDB[i];
// 						temp3 = potential[gridPoint+gridSizeX]*dielectricYDB[i];
// 						temp4 = potential[gridPoint-gridSizeX]*dielectricYmDB[i];
// 						temp5 = potential[gridPoint+gridSizeXY]*dielectricZDB[i];
// 						temp6 = potential[gridPoint-gridSizeXY]*dielectricZmDB[i];
// 						temp7 = potential[gridPoint]*om1;
// 						potential[gridPoint] = temp7+dielectricZDBSUM[i]*(temp1+temp2+temp3+temp4+temp5+temp6);
// 					}
// 					for(i=NoSingularNum[j];i<NoSingularNum[j+1];i++){
// 						gridPoint=IndexNoSingular[i];
// 						potential[gridPoint] = potential[gridPoint]*om1 + om2d6 * (potential[gridPoint+1] + potential[gridPoint-1] + potential[gridPoint+gridSizeX] + potential[gridPoint-gridSizeX] + potential[gridPoint+gridSizeXY] + potential[gridPoint-gridSizeXY]);
// 					}
// 					m_ContWorld->BorderExchange(potential);
// 				}
// 		
// 		//m_ContWorld->BorderExchange(potential);
// 		
// 				if(verbose||iteration==MaxIterations)if(iteration%ConvergenceCheck==0)
// 				{
// 					temp1=totalEnergy;
// 					totalEnergy=CalculateEnergyPAD(fpoh,potential,staticCharge,dielectricCh, IndexCharge, IndexSingular, ChargeNum[2], SingularNum[2]);
// 					totalChange = fabs(totalEnergy-temp1);
// 					relativeChange=totalChange/totalEnergy;
// 					ConvFac=totalChange;
// 					
// 					if(verbose)
// 						pnpPrintGroup0("<PoissonIterations Nit=\"%8d\" E=\"%20.16le\" dE=\"%17.8lg\" rel.E=\"%17.8le\" ConvFac=\"%17.8le\"/>\n", iteration, totalEnergy, totalChange, relativeChange,ConvFac);
// 
// // #ifdef WIN32
// //			 if(totalEnergy>1E13)
// // #else
// //			 if(totalEnergy>1E13||totalEnergy==NAN)
// // #endif
// 					if(totalEnergy>1E13)
// 					{
// 						fprintf(stdout,"totalEnergy>1E13 || totalEnergy==NAN />\n", totalEnergy);
// 						return EXIT_FAILURE;
// 					}
// 					
// 					if(ConvFac<Convergence&&iteration>MinIterations)iteration = MaxIterations+1;
// 				}
// 			}
// 	//Calculate last energy
// 			totalEnergy=CalculateEnergyPAD(fpoh,potential,staticCharge,dielectricCh, IndexCharge, IndexSingular,ChargeNum[2], SingularNum[2]);
// 	
// 			if(verbose)
// 				pnpPrintGroup0("<PoissonFinal E=\"%.10e\" Eerr=\"%.10e\" Niter=\"%d\"/>\n", totalEnergy, totalChange, iteration-1);
// 			return EXIT_SUCCESS;
// }

int PoissonBoltzmannSolver::LinearPBSolverNIB(int Niter)
{
	DbgPrint2("PoissonBoltzmannSolver::LinearPBSolver\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld,"m_ContWorld is not initialized\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld->Potential,"Potential is not initialize at m_ContWorld\n");
	
  //tmp var
	int iteration;
	int i,j;
	int GrdPnt;
	float tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7;
	
	//local vars
	float gridScale = m_ContWorld->GridScale;
	int GS_X = m_ContWorld->GridSize[0];
	int GS_Y = m_ContWorld->GridSize[1];
	int GS_Z = m_ContWorld->GridSize[2];
	int GS_XY = GS_X*GS_Y;
	int GS_XYZ = GS_XY*GS_Z;
	bool *PeriodicBoundaryCondition = m_ContWorld->PBC;
	//vars
	float fpoh = 4.0*M_PI*gridScale;
	float dynamicChargeFactor = 1.0/(COANGS*4.0*M_PI*DFACT*DFACT);
	float IonStrengthFactor = 1.0/(4.0*M_PI*gridScale);
	float om2 = Relaxation;
	float om1 = 1.0-om2;
	float om2d6 = om2/6.0;
	double totalEnergyOld=totalEnergy;
	
	
  //Iteration itself
	for(iteration=1;iteration<=Niter;iteration++)
	{
    //calculation over black and white nodes
		for(j=0;j<=1;j++)
		{
			for(i=SingularNum[j];i<SingularNum[j+1];i++)
			{
				GrdPnt=IndexSingular[i];
				tmp1 = potential[GrdPnt+1]*dielectricXS[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmS[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYS[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmS[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZS[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmS[i];				
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+tmp4;
				tmp3 = tmp5+tmp6;
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+QstS[i];								
				tmp1 = dielectricZSSUM[i]*(tmp1+tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1+tmp1;
			}
			for(i=ChargeNum[j];i<ChargeNum[j+1];i++){
				GrdPnt=IndexCharge[i];
				tmp1 = potential[GrdPnt+1]+potential[GrdPnt-1];
				tmp2 = potential[GrdPnt+GS_X]+potential[GrdPnt-GS_X];
				tmp3 = potential[GrdPnt+GS_XY]+potential[GrdPnt-GS_XY];
				tmp4 = potential[GrdPnt]*om1;
				tmp5 = om2d6*(tmp1+tmp2+tmp3+Qst[i]);
				//tmp6 = denominator[GrdPnt]*(tmp3+staticCharge[GrdPnt]);
				potential[GrdPnt] = tmp4+tmp5;
			}
			for(i=DielBoarderNum[j];i<DielBoarderNum[j+1];i++)
			{
				GrdPnt=IndexDielBoarder[i];
				tmp1 = potential[GrdPnt+1]*dielectricXDB[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmDB[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYDB[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmDB[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZDB[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmDB[i];
				tmp7 = potential[GrdPnt]*om1;
				potential[GrdPnt] = tmp7+dielectricZDBSUM[i]*(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6);
				//potential[GrdPnt]+=100.0;
			}
			for(i=NoSingularNum[j];i<NoSingularNum[j+1];i++)
			{
				GrdPnt=IndexNoSingular[i];
				potential[GrdPnt] = potential[GrdPnt]*om1 + om2d6 * (potential[GrdPnt+1] + potential[GrdPnt-1] + potential[GrdPnt+GS_X] + potential[GrdPnt-GS_X] + potential[GrdPnt+GS_XY] + potential[GrdPnt-GS_XY]);
				//potential[GrdPnt]+=10.0;
			}
      //PB
			for(i=PBZoneNum[j];i<PBZoneNum[j+1];i++)
			{
				GrdPnt=IndexPBZone[i];
				potential[GrdPnt] = potential[GrdPnt]*om1 + om2d6LPB * (potential[GrdPnt+1] + potential[GrdPnt-1] + potential[GrdPnt+GS_X] + potential[GrdPnt-GS_X] + potential[GrdPnt+GS_XY] + potential[GrdPnt-GS_XY]);
			}
			m_ContWorld->BorderExchange(potential);
		}
    //checking and printing energy
		if((verbose&&(iteration%ConvergenceCheck==0))||iteration==Niter)
		{
			CalcSystemEnergy(iteration);
			relativeChange=totalChange/totalEnergy;
			
			//if(verbose)
			//	pnpPrintGroup0("<PoissonBoltzmannIterations Nit=\"%8d\" E=\"%20.16e\" dE=\"%.4e\" rel.E=\"%.4e\" ConvFac=\"%.4e\"/>\n", iteration, totalEnergy, totalChange, relativeChange,ConvFac);
			if(iteration/ConvergenceCheck<=1)
				{
					pnpPrintGroup0("LPB   =========================================================================\n");
					pnpPrintGroup0("LPB    %9s %22s %12s %12s %12s\n","Iteration", "Energy,kT","dE","rel.E","ConvFac");
					pnpPrintGroup0("LPB   -------------------------------------------------------------------------\n");
				}
				pnpPrintGroup0("LPB    %9d %22.14e %12.4e %12.4e %12.4e\n", iteration, totalEnergy, totalChange, relativeChange,ConvFac);
			if(totalEnergy>1E13)
			{
				//fprintf(stdout,"totalEnergy>1E13 || totalEnergy==NAN />\n");
				pnpPrintGroup0("LPB   -------------------------------------------------------------------------\n");
				pnpPrintGroup0("LPB    ERROR: Linear Poisson-Boltzmann Solver has diverged, try smaller relaxation\n");
				pnpPrintGroup0("LPB   =========================================================================\n");
				
				return EXIT_FAILURE;
			}
			if(ConvFac<Convergence)iteration = Niter+1;
		}
	}
	if(verbose)
	{
		//pnpPrintGroup0("<PoissonBoltzmannFinal E=\"%.10e\" Eerr=\"%.10e\" Niter=\"%d\"/>\n", totalEnergy, totalChange, iteration-1);
		pnpPrintGroup0("LPB   -------------------------------------------------------------------------\n");
		pnpPrintGroup0("LPB    Linear Poisson-Boltzmann Solver\n");
		pnpPrintGroup0("LPB    Results: E=%.14e Niter=%d\n", totalEnergy, iteration-1);
		pnpPrintGroup0("LPB   =========================================================================\n");
	}
	if(MaxIterationsNPB==0&&Convergence>0.0&&ConvFac>Convergence)
	{
		pnpError("Convergence>0.0&&ConvFac>Convergence and no NPB steps");
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
int PoissonBoltzmannSolver::NonlinearPBSolverNIB(int Niter)
{
	DbgPrint2("PoissonBoltzmannSolver::NonlinearPBSolver\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld,"m_ContWorld is not initialized\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld->Potential,"Potential is not initialize at m_ContWorld\n");
	
  //tmp var
	int iteration;
	int i,j;
	int GrdPnt;
	float tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7;
	
	//local vars
	float gridScale = m_ContWorld->GridScale;
	int GS_X = m_ContWorld->GridSize[0];
	int GS_Y = m_ContWorld->GridSize[1];
	int GS_Z = m_ContWorld->GridSize[2];
	int GS_XY = GS_X*GS_Y;
	int GS_XYZ = GS_XY*GS_Z;
	bool *PeriodicBoundaryCondition = m_ContWorld->PBC;
	//vars
	float fpoh = 4.0*M_PI*gridScale;
	float dynamicChargeFactor = 1.0/(COANGS*4.0*M_PI*DFACT*DFACT);
	float IonStrengthFactor = 1.0/(4.0*M_PI*gridScale);
	float om2 = Relaxation;
	float om1 = 1.0-om2;
	float om2d6 = om2/6.0;
	double totalEnergyOld=totalEnergy;
	
	//PB
	float OneSixth=1.0/6.0,OneTwentyeth=1.0/20.0;
	float A=2.0*IonicStrength/m_ContWorld->NIndexing->Eps[(m_ContWorld->NIndexing->NIndex[0]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
	float om2d6PB=om2/(6.0+A);
	float coef=fpoh*COANGS/(m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
	DbgPrint2("I=%f[int]=%f[M] A=%f Epsilon=%f[int]=%f\n",IonicStrength,IonicStrength/coef,A,m_ContWorld->NIndexing->Eps[(m_ContWorld->NIndexing->NIndex[0]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft],m_ContWorld->NIndexing->Eps[(m_ContWorld->NIndexing->NIndex[0]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft]*EPKT);
	
  //Iteration itself
	for(iteration=1;iteration<=Niter;iteration++)
	{
    //calculation over black and white nodes
		for(j=0;j<=1;j++)
		{
			for(i=SingularNum[j];i<SingularNum[j+1];i++)
			{
				GrdPnt=IndexSingular[i];
				tmp1 = potential[GrdPnt+1]*dielectricXS[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmS[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYS[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmS[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZS[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmS[i];				
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+tmp4;
				tmp3 = tmp5+tmp6;
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+QstS[i];								
				tmp1 = dielectricZSSUM[i]*(tmp1+tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1+tmp1;
			}
			for(i=ChargeNum[j];i<ChargeNum[j+1];i++){
				GrdPnt=IndexCharge[i];
				tmp1 = potential[GrdPnt+1]+potential[GrdPnt-1];
				tmp2 = potential[GrdPnt+GS_X]+potential[GrdPnt-GS_X];
				tmp3 = potential[GrdPnt+GS_XY]+potential[GrdPnt-GS_XY];
				tmp4 = potential[GrdPnt]*om1;
				tmp5 = om2d6*(tmp1+tmp2+tmp3+Qst[i]);
				//tmp6 = denominator[GrdPnt]*(tmp3+staticCharge[GrdPnt]);
				potential[GrdPnt] = tmp4+tmp5;
			}
			for(i=DielBoarderNum[j];i<DielBoarderNum[j+1];i++)
			{
				GrdPnt=IndexDielBoarder[i];
				tmp1 = potential[GrdPnt+1]*dielectricXDB[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmDB[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYDB[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmDB[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZDB[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmDB[i];
				tmp7 = potential[GrdPnt]*om1;
				potential[GrdPnt] = tmp7+dielectricZDBSUM[i]*(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6);
				//potential[GrdPnt]+=100.0;
			}
			for(i=NoSingularNum[j];i<NoSingularNum[j+1];i++)
			{
				GrdPnt=IndexNoSingular[i];
				potential[GrdPnt] = potential[GrdPnt]*om1 + om2d6 * (potential[GrdPnt+1] + potential[GrdPnt-1] + potential[GrdPnt+GS_X] + potential[GrdPnt-GS_X] + potential[GrdPnt+GS_XY] + potential[GrdPnt-GS_XY]);
				//potential[GrdPnt]+=10.0;
			}
      //PB
			for(i=PBZoneNum[j];i<PBZoneNum[j+1];i++)
			{
				GrdPnt=IndexPBZone[i];
				tmp1=potential[GrdPnt]*potential[GrdPnt];
				tmp2=1.0+OneSixth*(tmp1+OneTwentyeth*tmp1*tmp1);
				om2d6PB=om2/(6.0+A*tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1 + om2d6PB * (potential[GrdPnt+1] + potential[GrdPnt-1] + potential[GrdPnt+GS_X] + potential[GrdPnt-GS_X] + potential[GrdPnt+GS_XY] + potential[GrdPnt-GS_XY]);
			}
			m_ContWorld->BorderExchange(potential);
		}
    //checking and printing energy
		if((verbose&&(iteration%ConvergenceCheck==0))||iteration==Niter)
		{
			CalcSystemEnergy(iteration);
			relativeChange=totalChange/totalEnergy;
			
			if(iteration/ConvergenceCheck<=1)
				{
					pnpPrintGroup0("NLPB  =========================================================================\n");
					pnpPrintGroup0("NLPB   %9s %22s %12s %12s %12s\n","Iteration", "Energy,kT","dE","rel.E","ConvFac");
					pnpPrintGroup0("NLPB  -------------------------------------------------------------------------\n");
				}
				pnpPrintGroup0("NLPB   %9d %22.14e %12.4e %12.4e %12.4e\n", iteration, totalEnergy, totalChange, relativeChange,ConvFac);
			if(totalEnergy>1E13)
			{
				//fprintf(stdout,"totalEnergy>1E13 || totalEnergy==NAN />\n");
				pnpPrintGroup0("NLPB  -------------------------------------------------------------------------\n");
				pnpPrintGroup0("NLPB   ERROR: Linear Poisson-Boltzmann Solver has diverged, try smaller relaxation\n");
				pnpPrintGroup0("NLPB  =========================================================================\n");
				
				return EXIT_FAILURE;
			}
			if(ConvFac<Convergence)iteration = Niter+1;
		}
	}
	if(verbose)
	{
		//pnpPrintGroup0("<PoissonBoltzmannFinal E=\"%.10e\" Eerr=\"%.10e\" Niter=\"%d\"/>\n", totalEnergy, totalChange, iteration-1);
		pnpPrintGroup0("NLPB  -------------------------------------------------------------------------\n");
		pnpPrintGroup0("NLPB   Non-Linear Poisson-Boltzmann Solver\n");
		pnpPrintGroup0("NLPB   Results: E=%.14e Niter=%d\n", totalEnergy, iteration-1);
		pnpPrintGroup0("NLPB  =========================================================================\n");
	}
	if(Convergence>0.0&&ConvFac>Convergence)
		return EXIT_FAILURE;
  
	return EXIT_SUCCESS;
}
int PoissonBoltzmannSolver::CalcSystemEnergy(int iteration)
{
// 	CalcSystemEnergyFloat(iteration);
	
	//if(PMFCalculation0!=NULL)
	//	CalcSystemEnergy4PMF(iteration);
	//else
		//CalcSystemEnergyDouble(iteration);
		CalcSystemEnergyStdDevPhi(iteration);
		//CalcSystemEnergyMaxPhiChange(iteration);
// 	CalcSystemEnergyLongDouble(iteration);
// 	CalcSystemEnergy0(iteration);
// 	CalcSystemEnergy1(iteration);
 	//CalcSystemEnergy3(iteration);
	//if(iteration==MaxIterations)
	//	CalcSystemEnergyAnalizer(iteration);
	return EXIT_SUCCESS;
}
int PoissonBoltzmannSolver::CalcSystemEnergyStdDevPhi(int iteration)
{
	int i,GrdPnt;
	static double oldSumSQ=0.0,oldSumAbs=0.0,oldTotEn=0.0;
	double Dev=0.0;
	double tmp,SumSQ=0.0,SumAbs=0.0;
	double EnergyCharge = 0.0;
	double EnergySingular = 0.0;
	double EnergyQmob = 0.0;
	double OldTotalEnergy=totalEnergy;
	double  fpoh = 4.0*M_PI*m_ContWorld->GridScale;
	float dphi=0.0,maxdphi=0.0;
	int GrdPntMaxDPhi=-1;
	
	for(i=0;i<ChargeNum[2];i++)
	{
		GrdPnt=IndexCharge[i];
		
		tmp=double(potential[GrdPnt])*double(Qst[i])*double(dielectricCh[i]);
		EnergyCharge+=tmp;
		dphi=potential[GrdPnt]-PhiCharge[i];
		SumSQ+=dphi*dphi;
		PhiCharge[i]=potential[GrdPnt];
	}
	for(i=0;i<SingularNum[2];i++)
	{
		GrdPnt=IndexSingular[i];
		tmp=double(QstS[i])*double(potential[GrdPnt]);
		EnergySingular+=tmp;
		dphi=potential[GrdPnt]-PhiSingular[i];
		SumSQ+=dphi*dphi;
		PhiSingular[i]=potential[GrdPnt];
	}

	//DbgPrint0("CSE_0       :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n",(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));
	//DbgPrint0("CSE_MaxPhiChange %.5e at %d\n",maxdphi,GrdPntMaxDPhi);
	totalEnergy=(EnergyCharge+EnergySingular)/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac=sqrt(SumSQ/double(ChargeNum[2]+SingularNum[2]));
	m_ContWorld->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
/*float PoissonBoltzmannSolver::CalculateEnergyPAD(float fpoh,float *Potential,float *StaticCharge,float *Epsilon,int *IndexCharge, int *IndexSingular,int ChargeNum,int SingularNum)
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
	if(m_ContWorld->NProcs!=1)
	{
		MPI::Status	status;
		double EnergyProc;
		if(m_ContWorld->MyRank==0)
		{
			for(i=1;i<m_ContWorld->NProcs;i++)
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
	m_ContWorld->SystemEnergy=Energy/(fpoh*2.0);
	return Energy/(fpoh*2.0);
}*/
int PoissonBoltzmannSolver::GuessNumberOfIteration()
{
	int GrdPnt;
	int i,j,k;
	int ix,iy,iz;
	float *sn[3];
	double tmp;
	for(i=0;i<3;i++)
	{
		sn[i]=new float [m_ContWorld->GridSize[i]];
		if(m_ContWorld->PBC[i])
		{
			for(j=0;j<m_ContWorld->GridSize[i];j++)
			{
				sn[i][j]=1.0/sqrt(float(m_ContWorld->GridSize[i]));
			}
		}
		else
		{
			sn[i][0]=0.0;
			sn[i][m_ContWorld->GridSize[i]-1]=0.0;
			for(j=1;j<m_ContWorld->GridSize[i]-1;j++)
			{
				tmp=M_PI*float(j)/float(m_ContWorld->GridSize[i]-1);
				sn[i][j]=sqrt(2.0)*sin(tmp)/sqrt(float(m_ContWorld->GridSize[i]-1));
			}
		}
	}
	float *tmpPot=new float[GS_XYZ];
	//float *tmpPot2=new float[GS_XYZ];
	
	for(ix=0;ix<m_ContWorld->GridSize[0];ix++)
		for(iy=0;iy<m_ContWorld->GridSize[1];iy++)
			for(iz=0;iz<m_ContWorld->GridSize[2];iz++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		tmpPot[GrdPnt]=sn[0][ix]*sn[1][iy]*sn[2][iz];
		//tmpPot2[GrdPnt]=tmpPot[GrdPnt];
	}
	float *PotPointer=potential;
	potential=tmpPot;
	
		float tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7;
		float om2 = 1.0;
		float om1 = 0.0;
		float om2d6 = 1.0/6.0;
		for(j=0;j<=1;j++){
			for(i=SingularNum[j];i<SingularNum[j+1];i++){
				GrdPnt=IndexSingular[i];
				tmp1 = potential[GrdPnt+1]*dielectricXS[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmS[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYS[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmS[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZS[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmS[i];				
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+tmp4;
				tmp3 = tmp5+tmp6;
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3;
				tmp1 = dielectricZSSUM[i]*(tmp1+tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1+tmp1;
			}
			for(i=ChargeNum[j];i<ChargeNum[j+1];i++){
				GrdPnt=IndexCharge[i];
				tmp1 = potential[GrdPnt+1]+potential[GrdPnt-1];
				tmp2 = potential[GrdPnt+GS_X]+potential[GrdPnt-GS_X];
				tmp3 = potential[GrdPnt+GS_XY]+potential[GrdPnt-GS_XY];
				tmp4 = potential[GrdPnt]*om1;
				tmp5 = om2d6*(tmp1+tmp2+tmp3);
				potential[GrdPnt] = tmp4+tmp5;
			}
			for(i=DielBoarderNum[j];i<DielBoarderNum[j+1];i++){
				GrdPnt=IndexDielBoarder[i];
				tmp1 = potential[GrdPnt+1]*dielectricXDB[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmDB[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYDB[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmDB[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZDB[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmDB[i];
				tmp7 = potential[GrdPnt]*om1;
				potential[GrdPnt] = tmp7+dielectricZDBSUM[i]*(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6);
			}
			for(i=NoSingularNum[j];i<NoSingularNum[j+1];i++){
				GrdPnt=IndexNoSingular[i];
				potential[GrdPnt] = potential[GrdPnt]*om1 + om2d6 * (potential[GrdPnt+1] + potential[GrdPnt-1] + potential[GrdPnt+GS_X] + potential[GrdPnt-GS_X] + potential[GrdPnt+GS_XY] + potential[GrdPnt-GS_XY]);
			}
			//PB
			for(i=PBZoneNum[j];i<PBZoneNum[j+1];i++)
			{
				GrdPnt=IndexPBZone[i];
				potential[GrdPnt] = potential[GrdPnt]*om1 + om2d6LPB * (potential[GrdPnt+1] + potential[GrdPnt-1] + potential[GrdPnt+GS_X] + potential[GrdPnt-GS_X] + potential[GrdPnt+GS_XY] + potential[GrdPnt-GS_XY]);
			}
			m_ContWorld->BorderExchange(potential);
		}
	
	NodeIndex* NIndex=m_ContWorld->NIndexing->NIndex;
	unsigned int BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;
	tmp=0.0;
	for(ix=0;ix<m_ContWorld->GridSize[0];ix++)
		for(iy=0;iy<m_ContWorld->GridSize[1];iy++)
			for(iz=0;iz<m_ContWorld->GridSize[2];iz++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		if(!(NIndex[GrdPnt]&BlackAndWhiteMask))
		{
			tmp+=potential[GrdPnt]*sn[0][ix]*sn[1][iy]*sn[2][iz];
		}
	}
	float spec=2.0*tmp;
	pnpPrint("Gauss-Seidel Spectral Radius: %f\n",spec);
	int iter=int(7.8/log(1.0 + sqrt(1.0-spec)));
	pnpPrint("Estimated Number of Iterations to Convergence: %d\n",iter);
	DeleteCVecArray(sn,3);
	delete [] tmpPot;
	potential=PotPointer;
	Relaxation=2.0/(1.0 + sqrt(1.0 - spec));
	pnpPrint("Estimated Relaxation Coefficient: %f\n",Relaxation);
	SetRelaxation(Relaxation);
	return iter;
}
