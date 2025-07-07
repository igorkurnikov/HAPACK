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

#include "pbwithljsolver.h"
#include "tinyxml.h"
#include "pnpdebug.h"
#include "contworld.h"
#include "math.h"
#include "pnpconstants.h"
#include "mapio.h"
#include "pnputil.h"

#include <Python.h>
#include <stdlib.h>
#include <float.h>

PBwithLJSolver::PBwithLJSolver()
 : GenericSolver()
{
	InitZero();
}


PBwithLJSolver::~PBwithLJSolver()
{
	Clear();
}
int PBwithLJSolver::InitZero()
{
	PnpsObject::SetName("PBwithLJSolver");
	SolverStr.push_back("Auto");
	SolverStr.push_back("NodeIndexBased");
	SolverStr.push_back("ArrayDirect");
	
	m_ContWorld=NULL;
	
	ConvergenceCheck=20;
	MaxIterations=300;
	Convergence=0.0;
	Relaxation=1.9;
	
	 //Initial values of variables:
	totalChange=0;
	relativeChange=0;
	totalEnergy=0;
	totalEnergyInd=0;
	solver=0;
	verbose=true;
	bPBExpAll=false;
	bAnalyseExplosion=false;
	
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
	//QPBZone=NULL;
	
	PBLJZoneNum[0]=0;
	PBLJZoneNum[1]=0;
	PBLJZoneNum[2]=0;
	IndexPBLJZone=NULL;
	
	PBLJQstZoneNum[0]=0;
	PBLJQstZoneNum[1]=0;
	PBLJQstZoneNum[2]=0;
	IndexPBLJQstZone=NULL;
	QstPBLJ=NULL;
	dielectricPBLJQst=NULL;
	
	PBLJDBZoneNum[0]=0;
	PBLJDBZoneNum[1]=0;
	PBLJDBZoneNum[2]=0;
	IndexPBLJDBZone=NULL;
	dielectricXPBLJDB=NULL;
	dielectricYPBLJDB=NULL;
	dielectricZPBLJDB=NULL;
	dielectricZPBLJDBSUM=NULL;
	dielectricXmPBLJDB=NULL;
	dielectricYmPBLJDB=NULL;
	dielectricZmPBLJDB=NULL;
	
	PBLJDBQstZoneNum[0]=0;
	PBLJDBQstZoneNum[1]=0;
	PBLJDBQstZoneNum[2]=0;
	IndexPBLJDBQstZone=NULL;
	QstPBLJDB=NULL;
	dielectricXPBLJDBQst=NULL;
	dielectricYPBLJDBQst=NULL;
	dielectricZPBLJDBQst=NULL;
	dielectricZPBLJDBQstSUM=NULL;
	dielectricXmPBLJDBQst=NULL;
	dielectricYmPBLJDBQst=NULL;
	dielectricZmPBLJDBQst=NULL;
	int i;
	for(i=0;i<PBwithLJConvFacMaxHistory;i++)
		ConvFacHistory[i]=1e10;
	potential=NULL;
	return EXIT_SUCCESS;
}
int PBwithLJSolver::Clear()
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
	DeleteCArray(IndexPBLJZone);
	DeleteCArray(IndexPBLJDBZone);
	DeleteCArray(dielectricXPBLJDB);
	DeleteCArray(dielectricYPBLJDB);
	DeleteCArray(dielectricZPBLJDB);
	DeleteCArray(dielectricZPBLJDBSUM);
	DeleteCArray(dielectricXmPBLJDB);
	DeleteCArray(dielectricYmPBLJDB);
	DeleteCArray(dielectricZmPBLJDB);
	DeleteCArray(IndexPBLJQstZone);
	DeleteCArray(QstPBLJ);
	DeleteCArray(dielectricPBLJQst);
	DeleteCArray(IndexPBLJDBQstZone);
	DeleteCArray(QstPBLJDB);
	DeleteCArray(dielectricXPBLJDBQst);
	DeleteCArray(dielectricYPBLJDBQst);
	DeleteCArray(dielectricZPBLJDBQst);
	DeleteCArray(dielectricZPBLJDBQstSUM);
	DeleteCArray(dielectricXmPBLJDBQst);
	DeleteCArray(dielectricYmPBLJDBQst);
	DeleteCArray(dielectricZmPBLJDBQst);
	DeleteCArray(PhiCharge);
	DeleteCArray(PhiSingular);
	
	
	return EXIT_SUCCESS;
}
int PBwithLJSolver::LoadParamFromPyDict(PyObject *dict)
{
	PyObject *p;

	haPyDict_GetItemAsInt(dict,"MaxIterations",&MaxIterations);
	haPyDict_GetItemAsFloat(dict,"Relaxation",&Relaxation);
	haPyDict_GetItemAsFloat(dict,"Convergence",&Convergence);
	haPyDict_GetItemAsFloat(dict,"Tolerance",&Convergence);
	haPyDict_GetItemAsInt(dict,"ConvergenceCheck",&ConvergenceCheck);
	haPyDict_GetItemAsInt(dict,"solver",&solver);
	haPyDict_GetItemAsBool(dict,"Verbose",&verbose);
	haPyDict_GetItemAsBool(dict,"bPBExpAll",&bPBExpAll);
	haPyDict_GetItemAsBool(dict,"bAnalyseExplosion",&bAnalyseExplosion);

	DbgPrint2("MaxIterations=%d  Relaxation=%f Convergence=%f ConvergenceCheck=%d\n",MaxIterations,Relaxation,Convergence,ConvergenceCheck);
	DbgPrint2("solver=%d  Verbose=%d bPBExpAll=%d bAnalyseExplosion=%d\n",solver,(int)verbose,(int)bPBExpAll,(int)bAnalyseExplosion);

	ShowParameters();
	SetRelaxation(Relaxation);
	if(dielectricZSSUM!=NULL)ShowProperties();
	return EXIT_SUCCESS;
}
int PBwithLJSolver::SetRelaxation(float _Relaxation)
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
	if(dielectricZPBLJDBSUM!=NULL)
	{
		for(i=0;i<PBLJDBZoneNum[2];i++)
		{
			gridPoint=IndexPBLJDBZone[i];
			dielectricZPBLJDBSUM[i]=dielectricXPBLJDB[i]+dielectricYPBLJDB[i]+dielectricZPBLJDB[i]+dielectricXmPBLJDB[i]+dielectricYmPBLJDB[i]+dielectricZmPBLJDB[i];
			dielectricZPBLJDBSUM[i]=Relaxation/dielectricZPBLJDBSUM[i];
		}
	}
	if(dielectricZPBLJDBQstSUM!=NULL)
	{
		for(i=0;i<PBLJDBQstZoneNum[2];i++)
		{
			gridPoint=IndexPBLJDBQstZone[i];
			dielectricZPBLJDBQstSUM[i]=dielectricXPBLJDBQst[i]+dielectricYPBLJDBQst[i]+dielectricZPBLJDBQst[i]+dielectricXmPBLJDBQst[i]+dielectricYmPBLJDBQst[i]+dielectricZmPBLJDBQst[i];
			dielectricZPBLJDBQstSUM[i]=Relaxation/dielectricZPBLJDBQstSUM[i];
		}
	}
	om2 = Relaxation;
	om1 = 1.0-om2;
	om2d6 = om2/6.0;
	
	return EXIT_SUCCESS;
}
int PBwithLJSolver::SetContWorld(ContWorld* _world)
{
	int i,j;
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
	
	if(m_ContWorld->PMF==NULL)
	{
		pnpWarning("PMF is not present, will assume that you run PB with diffrent ion sizes\n");
		m_ContWorld->PMF=new float*[m_ContWorld->NIonsTypes];
		for(j=0;j<m_ContWorld->NIonsTypes;j++)
		{
			m_ContWorld->PMF[j]=new float[GS_XYZ];
			for(i=0;i<GS_XYZ;i++)
				m_ContWorld->PMF[j][i]=0.0;
		}
	}
	return EXIT_SUCCESS;
}
int PBwithLJSolver::ShowParameters()
{
	DbgPrint2("PBwithLJSolver::ShowParameters\n");
	
	pnpPrintGroup0("\nParameters of PBSR solver set up\n");
	pnpPrintGroup0("    MaxIterations:.................... %d\n", MaxIterations);
	pnpPrintGroup0("    Convergence:...................... %.3e\n", Convergence);
	pnpPrintGroup0("    Relaxation:....................... %.5g\n", Relaxation);
	
	return EXIT_SUCCESS;
}
int PBwithLJSolver::ShowProperties()
{
	DbgPrint2("PBwithLJSolver::ShowProperties\n");
	return EXIT_SUCCESS;
}
int PBwithLJSolver::InitSolver()
{
	int status;
	bool bGuessNumberOfIteration=false;
	status=InitSolverNIB();
	return status;
}
int PBwithLJSolver::InitSolverNIB()
{
	DbgPrint2("PBwithLJSolver::InitSolverNIB\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld,"m_ContWorld is not initialized\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld->Potential,"Potential is not initialize at m_ContWorld\n");
	
	//temp vars
	int i,j,k,kgrid,jgrid,BlackOrWhite;
	int GrdPnt;
	NodeIndexing* NIndexing=m_ContWorld->NIndexing;
	NodeIndex* NIndex=m_ContWorld->NIndexing->NIndex;
	
	
	
	float **Vlj=m_ContWorld->PMF;
	
	for(i=0;i<PBwithLJConvFacMaxHistory;i++)
		ConvFacHistory[i]=1e10;
//	Clear();
	
	//count b/w regions
	NoSingularNum[0]=0;NoSingularNum[1]=0;NoSingularNum[2]=0;
	SingularNum[0]=0;SingularNum[1]=0;SingularNum[2]=0;
	DielBoarderNum[0]=0;DielBoarderNum[1]=0;DielBoarderNum[2]=0;
	ChargeNum[0]=0;ChargeNum[1]=0;ChargeNum[2]=0;
	PBZoneNum[0]=0;PBZoneNum[1]=0;PBZoneNum[2]=0;
	PBLJZoneNum[0]=0;PBLJZoneNum[1]=0;PBLJZoneNum[2]=0;
	PBLJQstZoneNum[0]=0;PBLJQstZoneNum[1]=0;PBLJQstZoneNum[2]=0;
	PBLJDBZoneNum[0]=0;PBLJDBZoneNum[1]=0;PBLJDBZoneNum[2]=0;
	PBLJDBQstZoneNum[0]=0;PBLJDBQstZoneNum[1]=0;PBLJDBQstZoneNum[2]=0;
//	PBDBNum[0]=0;PBDBNum[1]=0;PBDBNum[2]=0;
	
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
	DbgPrint0("BnW(NIndex[0])=%d\n",int(NIndex[0]&BlackAndWhiteMask));
#ifdef MPI_PARALLEL
	pnpsapp->MyComGroup.Barrier();
#endif
	for(k=1;k<GS_Z-1;k++)
	{
		kgrid = k*GS_XY;
		for(j=1;j<GS_Y-1;j++)
		{
			jgrid = kgrid+j*GS_X;
			for(i=1;i<GS_X-1;i++)
			{
				GrdPnt = jgrid+i;
				if(NIndexing->GetConcFloat(0,GrdPnt)>0.0 ||  NIndexing->GetConcFloat(1,GrdPnt)>0.0)
				{
					if((Vlj[0][GrdPnt]!=0.0 && Vlj[1][GrdPnt]!=0.0)||(NIndexing->GetConcFloat(0,GrdPnt)==0.0 ||  NIndexing->GetConcFloat(1,GrdPnt)==0.0)||bPBExpAll)
					{
						if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
						{
							PBLJDBQstZoneNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							PBLJDBQstZoneNum[2]++;
						}
						else if(NIndex[GrdPnt]&specChargeMask)
						{
							PBLJQstZoneNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							PBLJQstZoneNum[2]++;
						}
						else if(NIndex[GrdPnt]&DielBoarderMask)
						{
							PBLJDBZoneNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							PBLJDBZoneNum[2]++;
						}
						else
						{
							PBLJZoneNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							PBLJZoneNum[2]++;
						}
					}
					else
					{
						if(NIndexing->GetConcFloat(0,GrdPnt)==0.0 ||  NIndexing->GetConcFloat(1,GrdPnt)==0.0)
						{
							//pnpError("node %d: only one of the ions is present on the node(nonlinear statdard PB region), current version of PBLJ can not solve this case\n",GrdPnt);
							
						}
						else if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
						{
							//pnpError("node %d: Static charge, dielectric boarder and mobile ions are on the same node(nonlinear statdard PB region), current version of PBLJ can not solve this case\n",GrdPnt);
							PBLJDBQstZoneNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							PBLJDBQstZoneNum[2]++;
						}
						else if(NIndex[GrdPnt]&specChargeMask)
						{
							//pnpError("node %d: Static charge and mobile ions are on the same node(nonlinear statdard PB region), current version of PBLJ can not solve this case\n",GrdPnt);
							PBLJQstZoneNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							PBLJQstZoneNum[2]++;
						}
						else if(NIndex[GrdPnt]&DielBoarderMask)
						{
							//pnpError("node %d: dielectric boarder and mobile ions are on the same node(nonlinear statdard PB region), current version of PBLJ can not solve this case\n",GrdPnt);
							PBLJDBZoneNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							PBLJDBZoneNum[2]++;
						}
						else
						{
							PBZoneNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
							PBZoneNum[2]++;
						}
					}
				}
				else
				{
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
						DielBoarderNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
						DielBoarderNum[2]++;
					}
					else
					{
						NoSingularNum[1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
						NoSingularNum[2]++;
					}
				}
			}
		}
	}
	pnpPrint0("Initialization of PBLJ Solver:\n");
	pnpPrint0("	Qst+Diel.Boarder:........... %d\n",SingularNum[2]);
	pnpPrint0("	Qst:........................ %d\n",ChargeNum[2]);
	pnpPrint0("	Diel.Boarder:............... %d\n",DielBoarderNum[2]);
	pnpPrint0("	No Qst + No Diel.Boarder:... %d\n",NoSingularNum[2]);
	pnpPrint0("	Qmob+Qst+Diel.Boarder(PBLJ): %d\n",PBLJDBQstZoneNum[2]);
	pnpPrint0("	Qmob+Qst(PBLJ):............. %d\n",PBLJQstZoneNum[2]);
	pnpPrint0("	Qmob+Diel.Boarder(PBLJ):.... %d\n",PBLJDBZoneNum[2]);
	pnpPrint0("	Qmob(PBLJ):................. %d\n",PBLJZoneNum[2]);
	pnpPrint0("	Qmob(PB):................... %d\n",PBZoneNum[2]);
	i=SingularNum[2] + ChargeNum[2] + DielBoarderNum[2] + NoSingularNum[2] + PBLJDBZoneNum[2] + PBLJZoneNum[2] + PBZoneNum[2]+PBLJDBQstZoneNum[2]+PBLJQstZoneNum[2];
	j=(GS_X-2)*(GS_Y-2)*(GS_Z-2);
	pnpPrint0("Totally %d nodes for PBLJ, total calculatable nodes: %d\n", i,j);
	if(i!=j)
		pnpError("There are some nodes which is not in the calculations, it is ok if only part of the system is caclualted, in other cases probably such situation is not yet implemented\n");
	pnpPrint0("\n");
	pnpPrint0("Distrebution among black and white nodes\n");
	pnpPrint0("	Qst+Diel.Boarder:........... %15d %15d\n",SingularNum[1],SingularNum[2]-SingularNum[1]);
	pnpPrint0("	Qst:........................ %15d %15d\n",ChargeNum[1],ChargeNum[2]-ChargeNum[1]);
	pnpPrint0("	Diel.Boarder:............... %15d %15d\n",DielBoarderNum[1],DielBoarderNum[2]-DielBoarderNum[1]);
	pnpPrint0("	No Qst + No Diel.Boarder:... %15d %15d\n",NoSingularNum[1],NoSingularNum[2]-NoSingularNum[1]);
	pnpPrint0("	Qmob+Qst+Diel.Boarder(PBLJ): %15d %15d\n",PBLJDBQstZoneNum[1],PBLJDBQstZoneNum[2]-PBLJDBQstZoneNum[1]);
	pnpPrint0("	Qmob+Qst(PBLJ):............. %15d %15d\n",PBLJQstZoneNum[1],PBLJQstZoneNum[2]-PBLJQstZoneNum[1]);
	pnpPrint0("	Qmob+Diel.Boarder(PBLJ):.... %15d %15d\n",PBLJDBZoneNum[1],PBLJDBZoneNum[2]-PBLJDBZoneNum[1]);
	pnpPrint0("	Qmob(PBLJ):................. %15d %15d\n",PBLJZoneNum[1],PBLJZoneNum[2]-PBLJZoneNum[1]);
	pnpPrint0("	Qmob(PB):................... %15d %15d\n",PBZoneNum[1],PBZoneNum[2]-PBZoneNum[1]);
	
	PNP_EXIT_FAIL_NULL((IndexNoSingular = new int[NoSingularNum[2]]), "No memory available\n");
	
	PNP_EXIT_FAIL_NULL((IndexDielBoarder = new int[DielBoarderNum[2]]), "No memory available\n");
	
	PNP_EXIT_FAIL_NULL((IndexCharge = new int[ChargeNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((Qst = new float[ChargeNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricCh = new float[ChargeNum[2]]), "No memory available\n");
	
	PNP_EXIT_FAIL_NULL((IndexSingular = new int[SingularNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((QstS = new float[SingularNum[2]]), "No memory available\n");
	
	
	PNP_EXIT_FAIL_NULL((IndexPBZone = new int[PBZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((IndexPBLJZone = new int[PBLJZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((IndexPBLJQstZone = new int[PBLJQstZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((QstPBLJ = new float[PBLJQstZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricPBLJQst = new float[PBLJQstZoneNum[2]]), "No memory available\n");
	
	PNP_EXIT_FAIL_NULL((IndexPBLJDBQstZone = new int[PBLJDBQstZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((QstPBLJDB = new float[PBLJDBQstZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((IndexPBLJDBZone = new int[PBLJDBZoneNum[2]]), "No memory available\n");
	
	PNP_EXIT_FAIL_NULL((PhiSingular = new float[SingularNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((PhiCharge = new float[ChargeNum[2]]), "No memory available\n");
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
	int iPBLJZone=0,iPBLJZone2=PBLJZoneNum[1];
	int iPBLJQstZone=0,iPBLJQstZone2=PBLJQstZoneNum[1];
	int iPBLJDBZone=0,iPBLJDBZone2=PBLJDBZoneNum[1];
	int iPBLJDBQstZone=0,iPBLJDBQstZone2=PBLJDBQstZoneNum[1];
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
				if(NIndexing->GetConcFloat(0,GrdPnt)>0.0 ||  NIndexing->GetConcFloat(1,GrdPnt)>0.0)
				{
					if((Vlj[0][GrdPnt]!=0.0 && Vlj[1][GrdPnt]!=0.0)||(NIndexing->GetConcFloat(0,GrdPnt)==0.0 ||  NIndexing->GetConcFloat(1,GrdPnt)==0.0)||bPBExpAll)
					{
						if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
						{
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexPBLJDBQstZone[iPBLJDBQstZone]=GrdPnt;
								QstPBLJDB[iPBLJDBQstZone]=Q[QCount];
								iPBLJDBQstZone++;
							}
							else
							{
								IndexPBLJDBQstZone[iPBLJDBQstZone2]=GrdPnt;
								QstPBLJDB[iPBLJDBQstZone2]=Q[QCount];
								iPBLJDBQstZone2++;
							}
						}
						else if(NIndex[GrdPnt]&specChargeMask)
						{
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexPBLJQstZone[iPBLJQstZone]=GrdPnt;
								dielectricPBLJQst[iPBLJQstZone]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
								QstPBLJ[iPBLJQstZone]=Q[QCount]/dielectricPBLJQst[iPBLJQstZone];
								iPBLJQstZone++;
							}
							else
							{
								IndexPBLJQstZone[iPBLJQstZone2]=GrdPnt;
								dielectricPBLJQst[iPBLJQstZone2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
								QstPBLJ[iPBLJQstZone2]=Q[QCount]/dielectricPBLJQst[iPBLJQstZone2];
								iPBLJQstZone2++;
							}
						}
						else if(NIndex[GrdPnt]&DielBoarderMask)
						{
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexPBLJDBZone[iPBLJDBZone]=GrdPnt;
								iPBLJDBZone++;
							}
							else
							{
								IndexPBLJDBZone[iPBLJDBZone2]=GrdPnt;
								iPBLJDBZone2++;
							}
						}
						else
						{
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexPBLJZone[iPBLJZone]=GrdPnt;
								iPBLJZone++;
							}
							else
							{
								IndexPBLJZone[iPBLJZone2]=GrdPnt;
								iPBLJZone2++;
							}
						}
					}
					else
					{
						if(NIndexing->GetConcFloat(0,GrdPnt)==0.0 ||  NIndexing->GetConcFloat(1,GrdPnt)==0.0)
						{
							//pnpError("node %d: only one of the ions is present on the node(nonlinear statdard PB region), current version of PBLJ can not solve this case\n",GrdPnt);
							
						}
						else if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
						{
							//pnpError("node %d: Static charge, dielectric boarder and mobile ions are on the same node(nonlinear statdard PB region), current version of PBLJ can not solve this case\n",GrdPnt);
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexPBLJDBQstZone[iPBLJDBQstZone]=GrdPnt;
								QstPBLJDB[iPBLJDBQstZone]=Q[QCount];
								iPBLJDBQstZone++;
							}
							else
							{
								IndexPBLJDBQstZone[iPBLJDBQstZone2]=GrdPnt;
								QstPBLJDB[iPBLJDBQstZone2]=Q[QCount];
								iPBLJDBQstZone2++;
							}
						}
						else if(NIndex[GrdPnt]&specChargeMask)
						{
							//pnpError("node %d: Static charge and mobile ions are on the same node(nonlinear statdard PB region), current version of PBLJ can not solve this case\n",GrdPnt);
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexPBLJQstZone[iPBLJQstZone]=GrdPnt;
								dielectricPBLJQst[iPBLJQstZone]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
								QstPBLJ[iPBLJQstZone]=Q[QCount]/dielectricPBLJQst[iPBLJQstZone];
								iPBLJQstZone++;
							}
							else
							{
								IndexPBLJQstZone[iPBLJQstZone2]=GrdPnt;
								dielectricPBLJQst[iPBLJQstZone2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
								QstPBLJ[iPBLJQstZone2]=Q[QCount]/dielectricPBLJQst[iPBLJQstZone2];
								iPBLJQstZone2++;
							}
						}
						else if(NIndex[GrdPnt]&DielBoarderMask)
						{
							//pnpError("node %d: dielectric boarder and mobile ions are on the same node(nonlinear statdard PB region), current version of PBLJ can not solve this case\n",GrdPnt);
							if(NIndex[GrdPnt]&BlackAndWhiteMask)
							{
								IndexPBLJDBZone[iPBLJDBZone]=GrdPnt;
								iPBLJDBZone++;
							}
							else
							{
								IndexPBLJDBZone[iPBLJDBZone2]=GrdPnt;
								iPBLJDBZone2++;
							}
						}
						else
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
					}
				}
				else
				{
					if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
					{
						if(NIndex[GrdPnt]&BlackAndWhiteMask)
						{
							IndexSingular[iSingular]=GrdPnt;
							QstS[iSingular]=Q[QCount];
							q+=Q[QCount];
							iSingular++;
						}
						else
						{
							IndexSingular[iSingular2]=GrdPnt;
							QstS[iSingular2]=Q[QCount];
							q+=Q[QCount];
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
							iCharge++;
						}
						else
						{
							IndexCharge[iCharge2]=GrdPnt;
							dielectricCh[iCharge2]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
							Qst[iCharge2]=Q[QCount]/dielectricCh[iCharge2];
							q+=Q[QCount];
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
				}
				if(NIndex[GrdPnt]&specChargeMask)
					QCount++;
			}
		}
	}
	DbgPrint0("q=%f QCount=%d",(float)q/4/M_PI/m_ContWorld->GridScale,QCount);
  //Allocate dielectric helping array
	PNP_EXIT_FAIL_NULL((dielectricXDB = new float[DielBoarderNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricYDB = new float[DielBoarderNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricZDB = new float[DielBoarderNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricXmDB = new float[DielBoarderNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricYmDB = new float[DielBoarderNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricZmDB = new float[DielBoarderNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricZDBSUM = new float[DielBoarderNum[2]]), "No memory available\n");
	
	PNP_EXIT_FAIL_NULL((dielectricXS = new float[SingularNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricYS = new float[SingularNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricZS = new float[SingularNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricXmS = new float[SingularNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricYmS = new float[SingularNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricZmS = new float[SingularNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricZSSUM = new float[SingularNum[2]]), "No memory available\n");
	
	PNP_EXIT_FAIL_NULL((dielectricXPBLJDB = new float[PBLJDBZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricYPBLJDB = new float[PBLJDBZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricZPBLJDB = new float[PBLJDBZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricXmPBLJDB = new float[PBLJDBZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricYmPBLJDB = new float[PBLJDBZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricZmPBLJDB = new float[PBLJDBZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricZPBLJDBSUM = new float[PBLJDBZoneNum[2]]), "No memory available\n");
	
	PNP_EXIT_FAIL_NULL((dielectricXPBLJDBQst = new float[PBLJDBQstZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricYPBLJDBQst = new float[PBLJDBQstZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricZPBLJDBQst = new float[PBLJDBQstZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricXmPBLJDBQst = new float[PBLJDBQstZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricYmPBLJDBQst = new float[PBLJDBQstZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricZmPBLJDBQst = new float[PBLJDBQstZoneNum[2]]), "No memory available\n");
	PNP_EXIT_FAIL_NULL((dielectricZPBLJDBQstSUM = new float[PBLJDBQstZoneNum[2]]), "No memory available\n");
	
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
	for(i=0;i<PBLJDBZoneNum[2];i++)
	{
		GrdPnt=IndexPBLJDBZone[i];
		dielectricXPBLJDB[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		dielectricYPBLJDB[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		dielectricZPBLJDB[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		dielectricXmPBLJDB[i]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		dielectricYmPBLJDB[i]=Eps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		dielectricZmPBLJDB[i]=Eps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		dielectricZPBLJDBSUM[i]=dielectricXPBLJDB[i]+dielectricYPBLJDB[i]+dielectricZPBLJDB[i]+dielectricXmPBLJDB[i]+dielectricYmPBLJDB[i]+dielectricZmPBLJDB[i];
		dielectricZPBLJDBSUM[i]=Relaxation/dielectricZPBLJDBSUM[i];
	}
	for(i=0;i<PBLJDBQstZoneNum[2];i++)
	{
		GrdPnt=IndexPBLJDBQstZone[i];
		dielectricXPBLJDBQst[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		dielectricYPBLJDBQst[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		dielectricZPBLJDBQst[i]=Eps[(NIndex[GrdPnt]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		dielectricXmPBLJDBQst[i]=Eps[(NIndex[GrdPnt-1]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft];
		dielectricYmPBLJDBQst[i]=Eps[(NIndex[GrdPnt-GS_X]&NodeIndexing::Epsilon1)>>NodeIndexing::Epsilon1Sft];
		dielectricZmPBLJDBQst[i]=Eps[(NIndex[GrdPnt-GS_XY]&NodeIndexing::Epsilon2)>>NodeIndexing::Epsilon2Sft];
		dielectricZPBLJDBQstSUM[i]=dielectricXPBLJDBQst[i]+dielectricYPBLJDBQst[i]+dielectricZPBLJDBQst[i]+dielectricXmPBLJDBQst[i]+dielectricYmPBLJDBQst[i]+dielectricZmPBLJDBQst[i];
		dielectricZPBLJDBQstSUM[i]=Relaxation/dielectricZPBLJDBQstSUM[i];
	}
	return EXIT_SUCCESS;
}
int PBwithLJSolver::Solve()
{
	if(bAnalyseExplosion)
		return AnalyseExplosion();
	DbgPrint2("PBwithLJSolver::Solve\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld,"m_ContWorld is not initialized\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld->Potential,"Potential is not initialize at m_ContWorld\n");
	
  //tmp var
	int iteration;
	int i,j;
	int GrdPnt;
	float tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7;
	float **Vlj=m_ContWorld->PMF;
	float *VljK=m_ContWorld->PMF[0];
	float *VljCl=m_ContWorld->PMF[1];
	
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
	float I;
	for(i=0;i<NodeIndexMaxValues;i++)
		if(m_ContWorld->NIndexing->C[i]>0.0)
	{
		I=m_ContWorld->NIndexing->C[i];
		break;
	}
	
	float OneSixth=1.0/6.0,OneTwentyeth=1.0/20.0;
	int iEpsOut= (m_ContWorld->NIndexing->NIndex[0]&NodeIndexing::Epsilon0) >> NodeIndexing::Epsilon0Sft;
#ifdef MPI_PARALLEL
	//will broadcast A 
	//becouse Epsilon at node 0 is equal to Eps out only at first process
	pnpsapp->MyComGroup.Bcast(&iEpsOut, 1, MPI::INT, 0);
	//pnpsapp->MyComGroup.Bcast(&A, 1, MPI::FLOAT, 0);
#endif
	float EpsOut=m_ContWorld->NIndexing->Eps[iEpsOut];
	float A=I/EpsOut;
	float A2=2.0*A;
	float om2d6PB=om2/(6.0+A2);
	float coef=fpoh*COANGS/(m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
	DbgPrint2("I=%f[int]=%f[M] A=%f Epsilon=%f[int]=%f\n",I,I/coef,A,EpsOut,EpsOut*EPKT);
	int countConvFacHistory=0;
  //Iteration itself
	for(iteration=1;iteration<=MaxIterations;iteration++)
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
				om2d6PB=om2/(6.0+A2*tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1 + om2d6PB * (potential[GrdPnt+1] + potential[GrdPnt-1] + potential[GrdPnt+GS_X] + potential[GrdPnt-GS_X] + potential[GrdPnt+GS_XY] + potential[GrdPnt-GS_XY]);
			}
			for(i=PBLJDBZoneNum[j];i<PBLJDBZoneNum[j+1];i++)
			{
				GrdPnt=IndexPBLJDBZone[i];
				tmp1 = potential[GrdPnt+1]*dielectricXPBLJDB[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmPBLJDB[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYPBLJDB[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmPBLJDB[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZPBLJDB[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmPBLJDB[i];				
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+tmp4;
				tmp3 = tmp5+tmp6;
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3-I*(exp(potential[GrdPnt]-VljCl[GrdPnt])-exp(-potential[GrdPnt]-VljK[GrdPnt]));
				tmp1 = dielectricZPBLJDBSUM[i]*(tmp1+tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1+tmp1;
			}
			for(i=PBLJDBQstZoneNum[j];i<PBLJDBQstZoneNum[j+1];i++)
			{
				GrdPnt=IndexPBLJDBQstZone[i];
				tmp1 = potential[GrdPnt+1]*dielectricXPBLJDBQst[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmPBLJDBQst[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYPBLJDBQst[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmPBLJDBQst[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZPBLJDBQst[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmPBLJDBQst[i];				
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+tmp4;
				tmp3 = tmp5+tmp6;
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+QstPBLJDB[i]-I*(exp(potential[GrdPnt]-VljCl[GrdPnt])-exp(-potential[GrdPnt]-VljK[GrdPnt]));
				tmp1 = dielectricZPBLJDBQstSUM[i]*(tmp1+tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1+tmp1;
			}
			for(i=PBLJQstZoneNum[j];i<PBLJQstZoneNum[j+1];i++)
			{
				GrdPnt=IndexPBLJQstZone[i];
				tmp1 = potential[GrdPnt+1]+potential[GrdPnt-1];
				tmp2 = potential[GrdPnt+GS_X]+potential[GrdPnt-GS_X];
				tmp3 = potential[GrdPnt+GS_XY]+potential[GrdPnt-GS_XY];
				tmp4 = potential[GrdPnt]*om1;
				tmp5 = om2d6*(tmp1+tmp2+tmp3+QstPBLJ[i]-A*(exp(potential[GrdPnt]-VljCl[GrdPnt])-exp(-potential[GrdPnt]-VljK[GrdPnt])));
				potential[GrdPnt] = tmp4+tmp5;
			}
			for(i=PBLJZoneNum[j];i<PBLJZoneNum[j+1];i++)
			{
				GrdPnt=IndexPBLJZone[i];
				tmp1 = potential[GrdPnt+1]+potential[GrdPnt-1];
				tmp2 = potential[GrdPnt+GS_X]+potential[GrdPnt-GS_X];
				tmp3 = potential[GrdPnt+GS_XY]+potential[GrdPnt-GS_XY];
				tmp4 = potential[GrdPnt]*om1;
				tmp5 = om2d6*(tmp1+tmp2+tmp3-A*(exp(potential[GrdPnt]-VljCl[GrdPnt])-exp(-potential[GrdPnt]-VljK[GrdPnt])));
				potential[GrdPnt] = tmp4+tmp5;
			}
			m_ContWorld->BorderExchange(potential);
		}
    //checking and printing energy
		if((verbose&&(iteration%ConvergenceCheck==0))||iteration==MaxIterations)
		{
			CalcSystemEnergy(iteration);
			relativeChange=totalChange/totalEnergy;
			
			if(verbose)
			{
				if(iteration/ConvergenceCheck<=1)
				{
					pnpPrintGroup0("PBSR  =========================================================================\n");
					pnpPrintGroup0("PBSR   %9s %22s %12s %12s %12s\n","Iteration", "Energy,kT","dE","rel.E","ConvFac");
					pnpPrintGroup0("PBSR  -------------------------------------------------------------------------\n");
				}
				pnpPrintGroup0("PBSR   %9d %22.14e %12.4e %12.4e %12.4e\n", iteration, totalEnergy, totalChange, relativeChange,ConvFac);
				//pnpPrintGroup0("<PBwithLJIterations Nit=\"%9d\" E=\"%20.16e\" dE=\"%.4e\" rel.E=\"%.4e\" ConvFac=\"%.4e\"/>\n", iteration, totalEnergy, totalChange, relativeChange,ConvFac);
			}
#if defined(_MSC_VER)
			if(fabs(totalEnergy) > 1E13 )
#else
			if(totalEnergy>1E13 || (!isfinite(totalEnergy)))
#endif
			{
				pnpPrintGroup0("PBSR  -------------------------------------------------------------------------\n");
				pnpPrintGroup0("PBSR   ERROR: PBSR has diverged, try smaller relaxation\n");
				pnpPrintGroup0("PBSR  =========================================================================\n");
				return EXIT_FAILURE;
			}
			ConvFacHistory[countConvFacHistory]=ConvFac;
			countConvFacHistory++;
			if(countConvFacHistory>=PBwithLJConvFacMaxHistory)
				countConvFacHistory=0;
			bool convereged=true;
			for(i=0;i<PBwithLJConvFacMaxHistory;i++)
				convereged=convereged&&(ConvFacHistory[i]<=Convergence);
			if(convereged)break;//iteration = MaxIterations+1;
/*			if(!pnpsapp->HaveTimeToRun())
			{
				pnpPrint("Run out of time\n");
				break;
			}*/
		}
	}

	pnpPrintGroup0("PBSR  -------------------------------------------------------------------------\n");
	if(verbose)
		pnpPrintGroup0("PBSR   Results: E=%.14e Niter=%d\n", totalEnergy, iteration-1);
	pnpPrintGroup0("PBSR  =========================================================================\n");
	if(Convergence>0.0&&ConvFac>Convergence)
		return EXIT_FAILURE;
  
	return EXIT_SUCCESS;
}
int pnpsIsFinite(float v)
{
	#if defined(_MSC_VER)
	if(_finite(v))
	#else
	if(isfinite(v))
	#endif
	{
		if(v<1.0e8&&v>-1.0e8)
			return 1;
		else
			return 0;
	}
	else
		return 0;
}
int PBwithLJSolver::AnalyseExplosion()
{
	DbgPrint2("PBwithLJSolver::AnalyseExplosion\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld,"m_ContWorld is not initialized\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld->Potential,"Potential is not initialize at m_ContWorld\n");
	
  //tmp var
	int iteration;
	int i,j;
	int GrdPnt;
	float tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7;
	float **Vlj=m_ContWorld->PMF;
	float *VljK=m_ContWorld->PMF[0];
	float *VljCl=m_ContWorld->PMF[1];
	
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
	
	//vf
	VectorField3D VF(GS_X,GS_Y,GS_Z, gridScale, 1);
	VF.FillValue(0.0);
	float *err=VF.V[0];
	//PB
	float I;
	for(i=0;i<NodeIndexMaxValues;i++)
		if(m_ContWorld->NIndexing->C[i]>0.0)
	{
		I=m_ContWorld->NIndexing->C[i];
		break;
	}
	
	float OneSixth=1.0/6.0,OneTwentyeth=1.0/20.0;
	int iEpsOut= (m_ContWorld->NIndexing->NIndex[0]&NodeIndexing::Epsilon0) >> NodeIndexing::Epsilon0Sft;
#ifdef MPI_PARALLEL
	//will broadcast A 
	//becouse Epsilon at node 0 is equal to Eps out only at first process
	pnpsapp->MyComGroup.Bcast(&iEpsOut, 1, MPI::INT, 0);
	//pnpsapp->MyComGroup.Bcast(&A, 1, MPI::FLOAT, 0);
#endif
	float EpsOut=m_ContWorld->NIndexing->Eps[iEpsOut];
	float A=I/EpsOut;
	float A2=2.0*A;
	float om2d6PB=om2/(6.0+A2);
	float coef=fpoh*COANGS/(m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
	DbgPrint2("I=%f[int]=%f[M] A=%f Epsilon=%f[int]=%f\n",I,I/coef,A,EpsOut,EpsOut*EPKT);
	int countConvFacHistory=0;
  //Iteration itself
	for(iteration=1;iteration<=MaxIterations;iteration++)
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
				if(!pnpsIsFinite(potential[GrdPnt]))
				{
					int ix=GrdPnt%GS_X;
					int iy=(GrdPnt%GS_XY)/GS_X;
					int iz=GrdPnt/GS_XY;
					pnpError("Pot at GrdPnt(SingularNum[%d]=%d) = %d [%d %d %d] is not finite (%f)\n",j,i,GrdPnt,ix,iy,iz,potential[GrdPnt]);
					err[GrdPnt]=1.0;
				}
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
				if(!pnpsIsFinite(potential[GrdPnt]))
				{
					int ix=GrdPnt%GS_X;
					int iy=(GrdPnt%GS_XY)/GS_X;
					int iz=GrdPnt/GS_XY;
					pnpError("Pot at GrdPnt(ChargeNum[%d]=%d) = %d [%d %d %d] is not finite (%f)\n",j,i,GrdPnt,ix,iy,iz,potential[GrdPnt]);
					err[GrdPnt]=1.0;
				}
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
				if(!pnpsIsFinite(potential[GrdPnt]))
				{
					int ix=GrdPnt%GS_X;
					int iy=(GrdPnt%GS_XY)/GS_X;
					int iz=GrdPnt/GS_XY;
					pnpError("Pot at GrdPnt(DielBoarderNum[%d]=%d) = %d [%d %d %d] is not finite (%f)\n",j,i,GrdPnt,ix,iy,iz,potential[GrdPnt]);
					err[GrdPnt]=1.0;
				}
			}
			for(i=NoSingularNum[j];i<NoSingularNum[j+1];i++)
			{
				GrdPnt=IndexNoSingular[i];
				potential[GrdPnt] = potential[GrdPnt]*om1 + om2d6 * (potential[GrdPnt+1] + potential[GrdPnt-1] + potential[GrdPnt+GS_X] + potential[GrdPnt-GS_X] + potential[GrdPnt+GS_XY] + potential[GrdPnt-GS_XY]);
				//potential[GrdPnt]+=10.0;
				if(!pnpsIsFinite(potential[GrdPnt]))
				{
					int ix=GrdPnt%GS_X;
					int iy=(GrdPnt%GS_XY)/GS_X;
					int iz=GrdPnt/GS_XY;
					pnpError("Pot at GrdPnt(NoSingularNum[%d]=%d) = %d [%d %d %d] is not finite (%f)\n",j,i,GrdPnt,ix,iy,iz,potential[GrdPnt]);
					err[GrdPnt]=1.0;
				}
			}
      //PB
			for(i=PBZoneNum[j];i<PBZoneNum[j+1];i++)
			{
				GrdPnt=IndexPBZone[i];
				tmp1=potential[GrdPnt]*potential[GrdPnt];
				tmp2=1.0+OneSixth*(tmp1+OneTwentyeth*tmp1*tmp1);
				om2d6PB=om2/(6.0+A2*tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1 + om2d6PB * (potential[GrdPnt+1] + potential[GrdPnt-1] + potential[GrdPnt+GS_X] + potential[GrdPnt-GS_X] + potential[GrdPnt+GS_XY] + potential[GrdPnt-GS_XY]);
				if(!pnpsIsFinite(potential[GrdPnt]))
				{
					int ix=GrdPnt%GS_X;
					int iy=(GrdPnt%GS_XY)/GS_X;
					int iz=GrdPnt/GS_XY;
					pnpError("Pot at GrdPnt(PBZoneNum[%d]=%d) = %d [%d %d %d] is not finite (%f)\n",j,i,GrdPnt,ix,iy,iz,potential[GrdPnt]);
					err[GrdPnt]=1.0;
				}
			}
			for(i=PBLJDBZoneNum[j];i<PBLJDBZoneNum[j+1];i++)
			{
				GrdPnt=IndexPBLJDBZone[i];
				tmp1 = potential[GrdPnt+1]*dielectricXPBLJDB[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmPBLJDB[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYPBLJDB[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmPBLJDB[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZPBLJDB[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmPBLJDB[i];				
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+tmp4;
				tmp3 = tmp5+tmp6;
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3-I*(exp(potential[GrdPnt]-VljCl[GrdPnt])-exp(-potential[GrdPnt]-VljK[GrdPnt]));
				tmp1 = dielectricZPBLJDBSUM[i]*(tmp1+tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1+tmp1;
				if(!pnpsIsFinite(potential[GrdPnt]))
				{
					int ix=GrdPnt%GS_X;
					int iy=(GrdPnt%GS_XY)/GS_X;
					int iz=GrdPnt/GS_XY;
					pnpError("Pot at GrdPnt(PBLJDBZoneNum[%d]=%d) = %d [%d %d %d] is not finite (%f)\n",j,i,GrdPnt,ix,iy,iz,potential[GrdPnt]);
					err[GrdPnt]=1.0;
				}
			}
			for(i=PBLJDBQstZoneNum[j];i<PBLJDBQstZoneNum[j+1];i++)
			{
				GrdPnt=IndexPBLJDBQstZone[i];
				tmp1 = potential[GrdPnt+1]*dielectricXPBLJDBQst[i];
				tmp2 = potential[GrdPnt-1]*dielectricXmPBLJDBQst[i];
				tmp3 = potential[GrdPnt+GS_X]*dielectricYPBLJDBQst[i];
				tmp4 = potential[GrdPnt-GS_X]*dielectricYmPBLJDBQst[i];
				tmp5 = potential[GrdPnt+GS_XY]*dielectricZPBLJDBQst[i];
				tmp6 = potential[GrdPnt-GS_XY]*dielectricZmPBLJDBQst[i];				
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+tmp4;
				tmp3 = tmp5+tmp6;
				tmp1 = tmp1+tmp2;
				tmp2 = tmp3+QstPBLJDB[i]-I*(exp(potential[GrdPnt]-VljCl[GrdPnt])-exp(-potential[GrdPnt]-VljK[GrdPnt]));
				tmp1 = dielectricZPBLJDBQstSUM[i]*(tmp1+tmp2);
				potential[GrdPnt] = potential[GrdPnt]*om1+tmp1;
				if(!pnpsIsFinite(potential[GrdPnt]))
				{
					int ix=GrdPnt%GS_X;
					int iy=(GrdPnt%GS_XY)/GS_X;
					int iz=GrdPnt/GS_XY;
					pnpError("Pot at GrdPnt(PBLJDBQstZoneNum[%d]=%d) = %d [%d %d %d] is not finite (%f)\n",j,i,GrdPnt,ix,iy,iz,potential[GrdPnt]);
					err[GrdPnt]=1.0;
				}
			}
			for(i=PBLJQstZoneNum[j];i<PBLJQstZoneNum[j+1];i++)
			{
				GrdPnt=IndexPBLJQstZone[i];
				tmp1 = potential[GrdPnt+1]+potential[GrdPnt-1];
				tmp2 = potential[GrdPnt+GS_X]+potential[GrdPnt-GS_X];
				tmp3 = potential[GrdPnt+GS_XY]+potential[GrdPnt-GS_XY];
				tmp4 = potential[GrdPnt]*om1;
				tmp5 = om2d6*(tmp1+tmp2+tmp3+QstPBLJ[i]-A*(exp(potential[GrdPnt]-VljCl[GrdPnt])-exp(-potential[GrdPnt]-VljK[GrdPnt])));
				potential[GrdPnt] = tmp4+tmp5;
				if(!pnpsIsFinite(potential[GrdPnt]))
				{
					int ix=GrdPnt%GS_X;
					int iy=(GrdPnt%GS_XY)/GS_X;
					int iz=GrdPnt/GS_XY;
					pnpError("Pot at GrdPnt(PBLJQstZoneNum[%d]=%d) = %d [%d %d %d] is not finite (%f)\n",j,i,GrdPnt,ix,iy,iz,potential[GrdPnt]);
					err[GrdPnt]=1.0;
				}
			}
			for(i=PBLJZoneNum[j];i<PBLJZoneNum[j+1];i++)
			{
				GrdPnt=IndexPBLJZone[i];
				float PotOld=potential[GrdPnt];
				tmp1 = potential[GrdPnt+1]+potential[GrdPnt-1];
				tmp2 = potential[GrdPnt+GS_X]+potential[GrdPnt-GS_X];
				tmp3 = potential[GrdPnt+GS_XY]+potential[GrdPnt-GS_XY];
				tmp4 = potential[GrdPnt]*om1;
				tmp6 = exp(potential[GrdPnt]-VljCl[GrdPnt])-exp(-potential[GrdPnt]-VljK[GrdPnt]);
				tmp5 = om2d6*(tmp1+tmp2+tmp3-A*tmp6);
				potential[GrdPnt] = tmp4+tmp5;
				if(!pnpsIsFinite(potential[GrdPnt]))
				{
					int ix=GrdPnt%GS_X;
					int iy=(GrdPnt%GS_XY)/GS_X;
					int iz=GrdPnt/GS_XY;
					pnpError("Pot at GrdPnt(PBLJZoneNum[%d]=%d) = %d [%d %d %d] is not finite (%f)\n",j,i,GrdPnt,ix,iy,iz,potential[GrdPnt]);
					pnpError("\texp=%f pot=%f VljK=%f VljCl=%f\n",tmp6,PotOld,VljK[GrdPnt],VljCl[GrdPnt]);
					pnpError("\texp(-pot-VljK)=%f exp(pot-VljCl)=%f\n",exp(-PotOld-VljK[GrdPnt]),exp(PotOld-VljCl[GrdPnt]));
					err[GrdPnt]=1.0;
				}
			}
			m_ContWorld->BorderExchange(potential);
		}
    //checking and printing energy
		if((verbose&&(iteration%ConvergenceCheck==0))||iteration==MaxIterations)
		{
			CalcSystemEnergy(iteration);
			relativeChange=totalChange/totalEnergy;
			
			if(verbose)
				pnpPrintGroup0("<PBwithLJIterations Nit=\"%8d\" E=\"%20.16e\" dE=\"%.4e\" rel.E=\"%.4e\" ConvFac=\"%.4e\"/>\n", iteration, totalEnergy, totalChange, relativeChange,ConvFac);
			if(totalEnergy>1E13)
			{
				fprintf(stdout,"totalEnergy>1E13 || totalEnergy==NAN />\n");
				return EXIT_FAILURE;
			}
			ConvFacHistory[countConvFacHistory]=ConvFac;
			countConvFacHistory++;
			if(countConvFacHistory>=PBwithLJConvFacMaxHistory)
				countConvFacHistory=0;
			bool convereged=true;
			for(i=0;i<PBwithLJConvFacMaxHistory;i++)
				convereged=convereged&&(ConvFacHistory[i]<=Convergence);
			if(convereged)break;//iteration = MaxIterations+1;
			/*if(!pnpsapp->HaveTimeToRun())
			{
				pnpPrint("Run out of time\n");
				break;
			}*/
		}
	}
	VF.WriteToFile("ErrPBwithLJSolver.bin");
	if(verbose)
		pnpPrintGroup0("<PBwithLJSolverFinal E=\"%.10e\" Eerr=\"%.10e\" Niter=\"%d\"/>\n", totalEnergy, totalChange, iteration-1);
	
	if(Convergence>0.0&&ConvFac>Convergence)
		return EXIT_FAILURE;
  
	return EXIT_SUCCESS;
}
int PBwithLJSolver::CalcSystemEnergy(int iteration)
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
int PBwithLJSolver::CalcSystemEnergyStdDevPhi(int iteration)
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
	for(i=0;i<PBLJQstZoneNum[2];i++)
	{
		GrdPnt=IndexPBLJQstZone[i];
		tmp=double(potential[GrdPnt])*double(QstPBLJ[i])*double(dielectricPBLJQst[i]);
		EnergySingular+=tmp;
	}
	for(i=0;i<PBLJDBQstZoneNum[2];i++)
	{
		GrdPnt=IndexPBLJDBQstZone[i];
		tmp=double(QstPBLJDB[i])*double(potential[GrdPnt]);
		EnergySingular+=tmp;
	}
	//DbgPrint0("CSE_0       :E=%.16e Eq=%.10e Esing=%.10e EQmob=%.10e\n",(EnergyCharge+EnergySingular+EnergyQmob)/(fpoh*2.0),(EnergyCharge)/(fpoh*2.0),(EnergySingular)/(fpoh*2.0),(EnergyQmob)/(fpoh*2.0));
	//DbgPrint0("CSE_MaxPhiChange %.5e at %d\n",maxdphi,GrdPntMaxDPhi);
	int ChargedNodes=ChargeNum[2]+SingularNum[2]+PBLJQstZoneNum[2]+PBLJDBQstZoneNum[2];
	//totalEnergy=(EnergyCharge+EnergySingular)/(fpoh*2.0);
	totalEnergy=EnergyCharge+EnergySingular;
	//DbgPrint0("totalEnergy=%g %g %g\n",totalEnergy,EnergyCharge,EnergySingular);
#ifdef MPI_PARALLEL
	int dest;
	pnpsapp->MyComGroup.Barrier();
	
	if(m_ContWorld->MyRank==0)
	{
		//DbgPrint0("totalEnergy=%g\n",totalEnergy);
		for(dest=1;dest<m_ContWorld->NProcs;dest++)
		{
			pnpsapp->MyComGroup.Recv(&tmp, 1, MPI::DOUBLE, dest, 0);
			//DbgPrint0("totalEnergy+=%g\n",tmp);
			totalEnergy+=tmp;
			pnpsapp->MyComGroup.Recv(&tmp, 1, MPI::DOUBLE, dest, 0);
			SumSQ+=tmp;
			pnpsapp->MyComGroup.Recv(&i, 1, MPI::INT, dest, 0);
			ChargedNodes+=i;
		}
	}
	else
	{
		pnpsapp->MyComGroup.Send(&totalEnergy, 1, MPI::DOUBLE, 0, 0);
		pnpsapp->MyComGroup.Send(&SumSQ, 1, MPI::DOUBLE, 0, 0);
		pnpsapp->MyComGroup.Send(&ChargedNodes, 1, MPI::INT, 0, 0);
	}
#endif
	totalEnergy=totalEnergy/(fpoh*2.0);
	totalChange=fabs(totalEnergy-OldTotalEnergy);
	ConvFac=sqrt(SumSQ/double(ChargedNodes));
#ifdef MPI_PARALLEL
	//int dest;
	pnpsapp->MyComGroup.Barrier();
			
	if(m_ContWorld->MyRank==0)
	{
		for(dest=1;dest<m_ContWorld->NProcs;dest++)
		{
			pnpsapp->MyComGroup.Send(&totalEnergy, 1, MPI::DOUBLE, dest, 0);
			pnpsapp->MyComGroup.Send(&SumSQ, 1, MPI::DOUBLE, dest, 0);
			pnpsapp->MyComGroup.Send(&ChargedNodes, 1, MPI::INT, dest, 0);
			pnpsapp->MyComGroup.Send(&ConvFac, 1, MPI::DOUBLE, dest, 0);
		}
	}
	else
	{
		pnpsapp->MyComGroup.Recv(&totalEnergy, 1, MPI::DOUBLE, 0, 0);
		pnpsapp->MyComGroup.Recv(&SumSQ, 1, MPI::DOUBLE, 0, 0);
		pnpsapp->MyComGroup.Recv(&ChargedNodes, 1, MPI::INT, 0, 0);
		pnpsapp->MyComGroup.Recv(&ConvFac, 1, MPI::DOUBLE, 0, 0);
	}
#endif
	m_ContWorld->SystemEnergy=totalEnergy;
	return EXIT_SUCCESS;
}
