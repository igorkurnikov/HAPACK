//
// C++ Implementation: pnpinterfaces
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifdef MPI_PARALLEL
#  include "mpi.h"
#endif

#include "pnputil.h"
#include "contworld.h"
#include "pnpdebug.h"
#include "math.h"
#include "pnpconstants.h"
#include "tinyxml.h"
#include "mapio.h"
#include "pbwithljsolver.h"
#include "nernstplanksolver.h"
#include <string.h>
#include <Python.h>

PNPUtil::PNPUtil()
{
}

PNPUtil::~PNPUtil()
{
}
int PNPUtil::RemoveLargePMFfromPNP(ContWorld *world,float LargePMF)
{
	pnpPrint0("PNPUtil::RemoveLargePMFfromPNP\n");
	int GS_X=world->GridSize[0];
	int GS_Y=world->GridSize[1];
	int GS_Z=world->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int i,j,k,ion;
	int GridPoint;

	int count=0;
	for(ion=0;ion<world->NIonsTypes;ion++)
	{
		//float* D=NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion0);
		for(i=0;i<GS_X;i++)
			for(j=0;j<GS_Y;j++)
				for(k=0;k<GS_Z;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(world->NIndexing->GetDiffFloat(ion,GridPoint)!=0.0
				&&world->PMF[ion][GridPoint]>LargePMF)
			{
				world->NIndexing->SetDiffToZero(ion,GridPoint);
				count++;
			}
		}
		pnpPrint0("%d points for ion %d was removed due to PMF>=%g\n",count,ion,LargePMF);
	}
	
	world->NIndexing->CalcDiffBoarder();
	return EXIT_SUCCESS;
}
int PNPUtil::RemoveSmallCfromPNP(ContWorld *world,float SmallC)
{
	pnpPrint0("PNPUtil::RemoveSmallCfromPNP\n");
	int GS_X=world->GridSize[0];
	int GS_Y=world->GridSize[1];
	int GS_Z=world->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int i,j,k,ion;
	int GridPoint;
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	
	int count=0;
	for(ion=0;ion<world->NIonsTypes;ion++)
	{
		//float* D=NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion0);
		for(i=0;i<GS_X;i++)
			for(j=0;j<GS_Y;j++)
				for(k=0;k<GS_Z;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(world->NIndexing->GetDiffFloat(ion,GridPoint)!=0.0
						&&world->C[ion][GridPoint]<=SmallC)
			{
				world->NIndexing->SetDiffToZero(ion,GridPoint);
				count++;
			}
		}
	}
	world->NIndexing->CalcDiffBoarder();
	pnpPrint0("%d points was removed due to C<=%g M (%g [int])\n",count,SmallC/coef,SmallC);
	return EXIT_SUCCESS;
}
int PNPUtil::RemoveSmallCfromPNPNearDiffBoarder(ContWorld *world,float SmallC)
{
	pnpPrint0("<RemoveSmallCfromPNPNearDiffBoarder>\n");
	int GS_X=world->GridSize[0];
	int GS_Y=world->GridSize[1];
	int GS_Z=world->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int i,j,k,ion;
	int GrdPnt;
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	
	NodeIndexing *NIndexing=world->NIndexing;
	
	
	for(ion=0;ion<world->NIonsTypes;ion++)
	{
		int count=0;
		int countTot=0;
		int cycle=0;
		do
		{
			count=0;
			for(i=0;i<GS_X;i++)
				for(j=0;j<GS_Y;j++)
					for(k=0;k<GS_Z;k++)
			{
				GrdPnt=i+j*GS_X+k*GS_XY;
				if(NIndexing->GetDiffFloat(ion,GrdPnt)!=0.0&& (NIndexing->GetDiffFloat(ion,GrdPnt+1)==0.0 || NIndexing->GetDiffFloat(ion,GrdPnt-1)==0.0 || NIndexing->GetDiffFloat(ion,GrdPnt+GS_X)==0.0 || NIndexing->GetDiffFloat(ion,GrdPnt-GS_X)==0.0 || NIndexing->GetDiffFloat(ion,GrdPnt+GS_XY)==0.0 || NIndexing->GetDiffFloat(ion,GrdPnt-GS_XY)==0.0) &&world->C[ion][GrdPnt]<=SmallC)
				{
					world->NIndexing->SetDiffToZero(ion,GrdPnt);
					count++;
					countTot++;
				}
			}
			cycle++;
		}
		while(count>0);
		pnpPrint0("for ion %d: %d points was removed due to C<=%g M (%g [int]) in %d cycles\n",ion,countTot,SmallC/coef,SmallC,cycle);
	}
	world->NIndexing->CalcDiffBoarder();
	pnpPrint0("</RemoveSmallCfromPNPNearDiffBoarder>\n");
	return EXIT_SUCCESS;
}
int PNPUtil::SetInternalCtoZero(ContWorld *world)
{
	pnpPrint0("PNPUtil::SetInternalCtoZero\n");
	int GS_X=world->GridSize[0];
	int GS_Y=world->GridSize[1];
	int GS_Z=world->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int i,j,k,ion;
	int GridPoint;

	int count=0;
	
	if(world->C==NULL)
		world->SetInitConcentrationFromNIndexing();
	for(ion=0;ion<world->NIonsTypes;ion++)
	{
		//float* D=NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion0);
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			world->C[ion][GridPoint]=0.0;
		}
	}
	return EXIT_SUCCESS;
}
int PNPUtil::RemoveLargedPMFfromPNP(ContWorld *m_ContWorld,float dPMF)
{
	pnpPrint0("PNPUtil::RemoveLargedPMFfromPNP\n");
	int i,j,k,ion;
	int GridPoint;
	double dPMFx;
	float **V=m_ContWorld->PMF;
	int GS_X=m_ContWorld->GridSize[0];
	int GS_Y=m_ContWorld->GridSize[1];
	int GS_Z=m_ContWorld->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int cycles=5000;
	for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
	{
		
		int iCycles;
		for(iCycles=0;iCycles<cycles;cycles++)
		{
			int Count=0;
			for(i=1;i<GS_X-1;i++)
				for(j=1;j<GS_Y-1;j++)
					for(k=1;k<GS_Z-1;k++)
			{
				GridPoint=i+j*GS_X+k*GS_XY;
				if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint)!=0.0)
				{
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-1)==0.0)
					{//boarder at -x
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint+1]);
						if(dPMFx>=dPMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+1)==0.0)
					{//boarder at +x
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint-1]);
						if(dPMFx>=dPMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_X)==0.0)
					{//boarder at -y
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint+GS_X]);
						if(dPMFx>=dPMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_X)==0.0)
					{//boarder at +y
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint-GS_X]);
						if(dPMFx>=dPMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_XY)==0.0)
					{//boarder at -z
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint+GS_XY]);
						if(dPMFx>=dPMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_XY)==0.0)
					{//boarder at +z
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint-GS_XY]);
						if(dPMFx>=dPMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							Count++;
						}
					}
				}
			}
			pnpPrint0("Removed %d nodes, for ion %d\n",Count,ion);
			if(Count==0)break;
		}
	}
	m_ContWorld->NIndexing->CalcDiffBoarder();
	return EXIT_SUCCESS;
}
bool PNPUtil::ShouldRemoveLargedPMFandLargePMFfromPNP(ContWorld *m_ContWorld,float dPMF,float PMF,int ion, int GrdPnt)
{
	float **C=m_ContWorld->C;
	float **V=m_ContWorld->PMF;
	int GS_X=m_ContWorld->GridSize[0];
	int GS_Y=m_ContWorld->GridSize[1];
	int GS_Z=m_ContWorld->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	NodeIndexing *NIndexing=m_ContWorld->NIndexing;
	
	float dix[6],dixt;
	dix[0] = NIndexing->GetDiffFloat(ion,GrdPnt+1)>0.0?
			0.5*(NIndexing->GetDiffFloat(ion,GrdPnt)+NIndexing->GetDiffFloat(ion,GrdPnt+1)):0.0;
	dix[1] = NIndexing->GetDiffFloat(ion,GrdPnt-1)>0.0?
			0.5*(NIndexing->GetDiffFloat(ion,GrdPnt)+NIndexing->GetDiffFloat(ion,GrdPnt-1)):0.0;
	dix[2] = NIndexing->GetDiffFloat(ion,GrdPnt+GS_X)>0.0?
			0.5*(NIndexing->GetDiffFloat(ion,GrdPnt)+NIndexing->GetDiffFloat(ion,GrdPnt+GS_X)):0.0;
	dix[3] = NIndexing->GetDiffFloat(ion,GrdPnt-GS_X)>0.0?
			0.5*(NIndexing->GetDiffFloat(ion,GrdPnt)+NIndexing->GetDiffFloat(ion,GrdPnt-GS_X)):0.0;
	dix[4] = NIndexing->GetDiffFloat(ion,GrdPnt+GS_XY)>0.0?
			0.5*(NIndexing->GetDiffFloat(ion,GrdPnt)+NIndexing->GetDiffFloat(ion,GrdPnt+GS_XY)):0.0;
	dix[5] = NIndexing->GetDiffFloat(ion,GrdPnt-GS_XY)>0.0?
			0.5*(NIndexing->GetDiffFloat(ion,GrdPnt)+NIndexing->GetDiffFloat(ion,GrdPnt-GS_XY)):0.0;
	dix[0] *= C[ion][GrdPnt+1];
	dix[1] *= C[ion][GrdPnt-1];
	dix[2] *= C[ion][GrdPnt+GS_X];
	dix[3] *= C[ion][GrdPnt-GS_X];
	dix[4] *= C[ion][GrdPnt+GS_XY];
	dix[5] *= C[ion][GrdPnt-GS_XY];
	dixt = dix[0]+dix[1]+dix[2]+dix[3]+dix[4]+dix[5];
	
	float PsiSum=dix[0]*V[ion][GrdPnt+1] + dix[1]*V[ion][GrdPnt-1] + dix[2]*V[ion][GrdPnt+GS_X] + dix[3]*V[ion][GrdPnt-GS_X] + dix[4]*V[ion][GrdPnt+GS_XY] + dix[5]*V[ion][GrdPnt-GS_XY];
	float Criteria=2.0*dixt+PsiSum-dixt*V[ion][GrdPnt];
	if(Criteria<=0.0)
		return true;
	else
		return false;
}
int PNPUtil::RemoveLargedPMFandLargePMFfromPNP(ContWorld *m_ContWorld,float dPMF,float PMF)
{
	pnpPrint0("PNPUtil::RemoveLargedPMFandLargePMFfromPNP\n");
	pnpPrint0("will remove points where dPMF>=%g and PMF>=%g\n",dPMF,PMF);
	int i,j,k,ion;
	int GridPoint;
	double dPMFx;
	float **V=m_ContWorld->PMF;
	int GS_X=m_ContWorld->GridSize[0];
	int GS_Y=m_ContWorld->GridSize[1];
	int GS_Z=m_ContWorld->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int cycles=5000;
	int TotalPointsRemoved=0;
	int *RemoveNode=new int[GS_XYZ];
	for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
	{
		
		int iCycles;
		for(iCycles=0;iCycles<cycles;cycles++)
		{
			for(i=0;i<GS_XYZ;i++)
				RemoveNode[i]=0;
			int Count=0;
			for(i=1;i<GS_X-1;i++)
				for(j=1;j<GS_Y-1;j++)
					for(k=1;k<GS_Z-1;k++)
			{
				GridPoint=i+j*GS_X+k*GS_XY;
				if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint)!=0.0&& ( m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-1)==0.0 || m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+1)==0.0 || m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_X)==0.0 || m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_X)==0.0 || m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_XY)==0.0 || m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_XY)==0.0))
				{
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-1)!=0.0)
					{//boarder at -x
						if(ShouldRemoveLargedPMFandLargePMFfromPNP(m_ContWorld, dPMF, PMF,ion,GridPoint-1))
						{
							RemoveNode[GridPoint-1]++;
							//m_ContWorld->NIndexing->SetDiffToZero(ion,);
							//Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+1)!=0.0)
					{//boarder at +x
						if(ShouldRemoveLargedPMFandLargePMFfromPNP(m_ContWorld, dPMF, PMF,ion,GridPoint+1))
						{
							RemoveNode[GridPoint+1]++;
							//m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint+1);
							//Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_X)!=0.0)
					{//boarder at -y
						if(ShouldRemoveLargedPMFandLargePMFfromPNP(m_ContWorld, dPMF, PMF,ion,GridPoint-GS_X))
						{
							RemoveNode[GridPoint-GS_X]++;
							//m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint-GS_X);
							//Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_X)!=0.0)
					{//boarder at +y
						if(ShouldRemoveLargedPMFandLargePMFfromPNP(m_ContWorld, dPMF, PMF,ion,GridPoint+GS_X))
						{
							RemoveNode[GridPoint+GS_X]++;
							//m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint+GS_X);
							//Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_XY)!=0.0)
					{//boarder at -z
						if(ShouldRemoveLargedPMFandLargePMFfromPNP(m_ContWorld, dPMF, PMF,ion,GridPoint-GS_XY))
						{
							RemoveNode[GridPoint-GS_XY]++;
							//m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint-GS_XY);
							//Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_XY)!=0.0)
					{//boarder at +z
						if(ShouldRemoveLargedPMFandLargePMFfromPNP(m_ContWorld, dPMF, PMF,ion,GridPoint+GS_XY))
						{
							RemoveNode[GridPoint+GS_XY]++;
							//m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint+GS_XY);
							//Count++;
						}
					}
				}
			}
			for(i=0;i<GS_XYZ;i++)
			{
				if(RemoveNode[i]>0)
				{
					m_ContWorld->NIndexing->SetDiffToZero(ion,i);
					Count++;
				}
			}
			pnpPrint0("Removed %d nodes, for ion %d\n",Count,ion);
			TotalPointsRemoved+=Count;
			if(Count==0)break;
		}
		pnpPrint0("Totally removed %d nodes, for ion %d\n",TotalPointsRemoved,ion);
	}
	delete [] RemoveNode;
	m_ContWorld->NIndexing->CalcDiffBoarder();
	return EXIT_SUCCESS;
}
int PNPUtil::RemoveLargedPMFandLargePMFfromPNPOld(ContWorld *m_ContWorld,float dPMF,float PMF)
{
	pnpPrint0("PNPUtil::RemoveLargedPMFandLargePMFfromPNP\n");
	pnpPrint0("will remove points where dPMF>=%g and PMF>=%g\n",dPMF,PMF);
	int i,j,k,ion;
	int GridPoint;
	double dPMFx;
	float **V=m_ContWorld->PMF;
	int GS_X=m_ContWorld->GridSize[0];
	int GS_Y=m_ContWorld->GridSize[1];
	int GS_Z=m_ContWorld->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int cycles=5000;
	int TotalPointsRemoved=0;
	int Count;
	for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
	{
		
		int iCycles;
		for(iCycles=0;iCycles<cycles;cycles++)
		{
			Count=0;
			for(i=1;i<GS_X-1;i++)
				for(j=1;j<GS_Y-1;j++)
					for(k=1;k<GS_Z-1;k++)
			{
				GridPoint=i+j*GS_X+k*GS_XY;
				if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint)!=0.0)
				{
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-1)==0.0)
					{//boarder at -x
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint+1]);
						if(dPMFx>=dPMF&&V[ion][GridPoint]>=PMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+1)==0.0)
					{//boarder at +x
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint-1]);
						if(dPMFx>=dPMF&&V[ion][GridPoint]>=PMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_X)==0.0)
					{//boarder at -y
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint+GS_X]);
						if(dPMFx>=dPMF&&V[ion][GridPoint]>=PMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_X)==0.0)
					{//boarder at +y
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint-GS_X]);
						if(dPMFx>=dPMF&&V[ion][GridPoint]>=PMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_XY)==0.0)
					{//boarder at -z
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint+GS_XY]);
						if(dPMFx>=dPMF&&V[ion][GridPoint]>=PMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_XY)==0.0)
					{//boarder at +z
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint-GS_XY]);
						if(dPMFx>=dPMF&&V[ion][GridPoint]>=PMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							Count++;
						}
					}
				}
			}
			pnpPrint0("Removed %d nodes, for ion %d\n",Count,ion);
			TotalPointsRemoved+=Count;
			if(Count==0)break;
		}
		pnpPrint0("Totally removed %d nodes, for ion %d\n",TotalPointsRemoved,ion);
	}
	
	m_ContWorld->NIndexing->CalcDiffBoarder();
	return EXIT_SUCCESS;
}

int PNPUtil::RemoveCavitiesAtDiffusionMap(ContWorld *m_ContWorld)
{
	int ion;
	int GSX=m_ContWorld->GridSize[0];
	int GSY=m_ContWorld->GridSize[1];
	int GSZ=m_ContWorld->GridSize[2];
	int i,j,k,i0;
	int GSXY=GSX*GSY;
	int GSXYZ=GSXY*GSZ;
	int GrdPnt;
	int gridpx,gridmx;
	int gridpy,gridmy;
	int gridpz,gridmz;
	

	for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
	{
		int count=1;
		int count0=0;
		int count1=1;
		int Cycles=0;
		int CyclesTot=0;
		int countNoCavities=0;
		//pnpPrint0("IonField[%d]=%d\n",ion,(unsigned int)m_ContWorld->NIndexing->IonField[ion]);
		float *Field=m_ContWorld->NIndexing->GetCMap(NodeIndexing::DiffConst,(NodeIndexing::NodeIndexDescriptor)m_ContWorld->NIndexing->IonField[ion]);
		
		for(i=1;i<GSX-1;i++)
			for(j=1;j<GSY-1;j++)
				for(k=1;k<GSZ-1;k++)
		{
			GrdPnt=i+j*GSX+k*GSXY;
			Field[GrdPnt]=-Field[GrdPnt];
		}
		int *list=new int[GSXYZ];
		for(i=0;i<GSXYZ;i++)list[i]=-1;
		list[0]=1+GSX+GSXY;
		Field[1+GSX+GSXY]=-Field[1+GSX+GSXY];
		do
		{
			for(i0=count0;i0<count1;i0++)
			{
				GrdPnt=list[i0];
				k=GrdPnt/GSXY;
				j=(GrdPnt%GSXY)/GSX;
				i=GrdPnt%GSX;
				
				gridpx = GrdPnt+1;
				gridmx = GrdPnt-1;
				gridpy = GrdPnt+GSX;
				gridmy = GrdPnt-GSX;
				gridpz = GrdPnt+GSXY;
				gridmz = GrdPnt-GSXY;
				if(Field[gridpx]<0.0&&i<GSX-1)
				{
					Field[gridpx]=-Field[gridpx];
					list[count]=gridpx;
					count++;
				}
				if(Field[gridmx]<0.0&&i>1)
				{
					Field[gridmx]=-Field[gridmx];
					list[count]=gridmx;
					count++;
				}
				if(Field[gridpy]<0.0&&j<GSY-1)
				{
					Field[gridpy]=-Field[gridpy];
					list[count]=gridpy;
					count++;
				}
				if(Field[gridmy]<0.0&&j>1)
				{
					Field[gridmy]=-Field[gridmy];
					list[count]=gridmy;
					count++;
				}
				if(Field[gridpz]<0.0&&k<GSZ-1)
				{
					Field[gridpz]=-Field[gridpz];
					list[count]=gridpz;
					count++;
				}
				if(Field[gridmz]<0.0&&k>1)
				{
					Field[gridmz]=-Field[gridmz];
					list[count]=gridmz;
					count++;
				}
			}
			count0=count1;
			count1=count;
			Cycles++;
		}
		while(count1>count0);
		for(i=0;i<GSXYZ;i++)
		{
			if(Field[i]<0.0)m_ContWorld->NIndexing->SetDiffToZero(ion,i);
		}
		delete [] list;
		delete [] Field;
		DbgPrint0("Made %d Cycles for cavity removing, %d nodes accesible\n",Cycles,count);
		
	}
	m_ContWorld->NIndexing->CalcDiffBoarder();
	return EXIT_SUCCESS;
}
int PNPUtil::ConvertPBLJresultsToDynamicCharge(ContWorld *m_ContWorld)
{
	int i;	
	float ch1=0.0,ch2=0.0;
	float fpoh= 4*M_PI*m_ContWorld->GridScale;
	float IStrength;
	float temp=fpoh*COANGS/(m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
	fprintf(stdout,"\nConvert IonStrength(from PBLJ) To Dynamic Charge\n");
	
	if(m_ContWorld->Potential==NULL)
	{
		pnpPrint("no potential will set to zero\n");
		m_ContWorld->Potential = new float[m_ContWorld->GS_XYZ];
		for(i=0;i<m_ContWorld->GS_XYZ;i++)
		{
			m_ContWorld->Potential[i]=0.0;
		}
	}
	if(m_ContWorld->C==NULL)
	{
		m_ContWorld->C=new float*[2];
		m_ContWorld->C[0]=NULL;
		m_ContWorld->C[1]=NULL;
	}
	if(m_ContWorld->C[0]==NULL)
	{
		if(!(m_ContWorld->C[0] = new float[m_ContWorld->GS_XYZ])){
			fprintf(stderr,"ERROR 104: No memory available\n");
			exit(104);
		}
	}
	if(m_ContWorld->C[1]==NULL)
	{
		if(!(m_ContWorld->C[1] = new float[m_ContWorld->GS_XYZ])){
			fprintf(stderr,"ERROR 104: No memory available\n");
			exit(104);
		}
	}
	for(i=0;i<NodeIndexMaxValues;i++)
	{
		if(m_ContWorld->NIndexing->C[i]!=0.0)
		{
			IStrength=m_ContWorld->NIndexing->C[i];
			break;
		}
	}
	pnpPrint("Will use bulk concentration: %f\n",IStrength/temp);
	
	for(i=0;i<m_ContWorld->GS_XYZ;i++)
	{
		if(m_ContWorld->NIndexing->GetConcFloat(0,i)>0.0)
		{
			m_ContWorld->C[0][i]=IStrength*exp(-m_ContWorld->Potential[i]-m_ContWorld->PMF[0][i]);		
			ch1+=m_ContWorld->C[0][i];
		}
		else
		{
			m_ContWorld->C[0][i]=0.0;
		}
		if(m_ContWorld->NIndexing->GetConcFloat(1,i)>0.0)
		{	
			m_ContWorld->C[1][i]=IStrength*exp(m_ContWorld->Potential[i]-m_ContWorld->PMF[1][i]);
			ch2+=m_ContWorld->C[1][i];
		}
		else
		{
			m_ContWorld->C[1][i]=0.0;
		}
	}
	fprintf(stdout,"Charge Is %g	+ (- %g )= %g\n",ch1/fpoh,ch2/fpoh,(ch1-ch2)/fpoh);
	return EXIT_SUCCESS;
}
#ifdef HARLEM_MOD
HaVec_float PNPUtil::GetRlim(ContWorld *world,const TiXmlElement *RlimElt)
{
	int GS_X=world->GridSize[0];
	int GS_Y=world->GridSize[1];
	int GS_Z=world->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int i,j,k,ion;
	int GridPoint;
	
	std::vector<float> zRlim;
	std::vector<float> Rlim;
	
	RlimElt->GetTwoVectorFloatElement("Rlim",&zRlim,&Rlim);
	
	float* AbsZ=new float[GS_Z];
	int Zhalf=GS_Z/2;
	for(k=0;k<GS_Z;k++)
	{
		AbsZ[k]=float(k-Zhalf)/world->GridScale;
		//pnpPrint("%d %f\n",k,AbsZ[k]);
	}
	int kStart,kEnd;
	int CurVal;
	float fac;
	
	HaVec_float RlimF;
	RlimF.newsize(GS_Z);
	
	k=0;
	while(AbsZ[k]<=zRlim[0])
	{
		RlimF[k]=Rlim[0];
		k++;
	}
	kStart=k;
	k=GS_Z-1;
	while(AbsZ[k]>=zRlim[zRlim.size()-1])
	{
		RlimF[k]=Rlim[zRlim.size()-1];
		k--;
	}
	kEnd=k;
	CurVal=0;
	for(k=kStart;k<=kEnd;k++)
	{
		if(AbsZ[k]>=zRlim[CurVal+1])CurVal++;
		fac=(Rlim[CurVal+1]-Rlim[CurVal])/(zRlim[CurVal+1]-zRlim[CurVal]);
		RlimF[k]=Rlim[CurVal]+(AbsZ[k]-zRlim[CurVal])*fac;
	}
	
	return RlimF;
}
#endif
int PNPUtil::ScaleDiffusionInTheChannel(ContWorld *world,float x0,float y0,std::vector<float>* zRlim,std::vector<float>* Rlim,std::vector<float>* zDiffScale,std::vector<float>* DiffScale)
{
	int GS_X=world->GridSize[0];
	int GS_Y=world->GridSize[1];
	int GS_Z=world->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int i,j,k,ion;
	int GridPoint;
	
	pnpPrint("i zDiffScale DiffScale\n");
	for(i=0;i<zDiffScale->size();i++)
	{
		pnpPrint("%d %f %f\n",i,(*zDiffScale)[i],(*DiffScale)[i]);
	}
	pnpPrint("i zRlim Rlim\n");
	for(i=0;i<zRlim->size();i++)
	{
		pnpPrint("%d %f %f\n",i,(*zRlim)[i],(*Rlim)[i]);
	}
	
	float* AbsZ=new float[GS_Z];
	int Zhalf=GS_Z/2;
	for(k=0;k<GS_Z;k++)
	{
		AbsZ[k]=float(k-Zhalf)/world->GridScale;
		//pnpPrint("%d %f\n",k,AbsZ[k]);
	}
	int kStart,kEnd;
	int CurVal;
	float fac;
	
	float* RlimF=new float[GS_Z];
	k=0;
	while(AbsZ[k]<=(*zRlim)[0])
	{
		RlimF[k]=(*Rlim)[0];
		k++;
	}
	kStart=k;
	k=GS_Z-1;
	while(AbsZ[k]>=(*zRlim)[zRlim->size()-1])
	{
		RlimF[k]=(*Rlim)[zRlim->size()-1];
		k--;
	}
	kEnd=k;
	CurVal=0;
	for(k=kStart;k<=kEnd;k++)
	{
		if(AbsZ[k]>=(*zRlim)[CurVal+1])CurVal++;
		fac=((*Rlim)[CurVal+1]-(*Rlim)[CurVal])/((*zRlim)[CurVal+1]-(*zRlim)[CurVal]);
		RlimF[k]=(*Rlim)[CurVal]+(AbsZ[k]-(*zRlim)[CurVal])*fac;
	}
	
	float* DiffScaleF=new float[GS_Z];
	k=0;
	while(AbsZ[k]<=(*zDiffScale)[0])
	{
		DiffScaleF[k]=(*DiffScale)[0];
		k++;
	}
	kStart=k;
	k=GS_Z-1;
	while(AbsZ[k]>=(*zDiffScale)[zDiffScale->size()-1])
	{
		DiffScaleF[k]=(*DiffScale)[zDiffScale->size()-1];
		k--;
	}
	kEnd=k;
	CurVal=0;
	for(k=kStart;k<=kEnd;k++)
	{
		if(AbsZ[k]>=(*zDiffScale)[CurVal+1])CurVal++;
		fac=((*DiffScale)[CurVal+1]-(*DiffScale)[CurVal])/((*zDiffScale)[CurVal+1]-(*zDiffScale)[CurVal]);
		DiffScaleF[k]=(*DiffScale)[CurVal]+(AbsZ[k]-(*zDiffScale)[CurVal])*fac;
	}
	
	pnpPrint("z[grid] Rlim DiffScale\n");
	for(k=0;k<GS_Z;k++)
	{
		pnpPrint("%d %f %f %f\n",k,AbsZ[k],RlimF[k],DiffScaleF[k]);
	}
	
	int status=ScaleDiffusionInTheChannel(world, x0, y0, RlimF,DiffScaleF);
	
	delete [] DiffScaleF;
	delete [] RlimF;
	delete [] AbsZ;
	return status;
}
int PNPUtil::SetCtoZeroWhereDZero(ContWorld *m_ContWorld)
{
	pnpPrint0("<PNPUtil::SetCtoZeroWhereDZero>\n");
	
	int GS_X=m_ContWorld->GridSize[0];
	int GS_Y=m_ContWorld->GridSize[1];
	int GS_Z=m_ContWorld->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int i,j,k,ion;
	int GridPoint;
	for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
	{
		for(i=0;i<GS_X;i++)
			for(j=0;j<GS_Y;j++)
				for(k=0;k<GS_Z;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint)==0.0)
			{
				m_ContWorld->C[ion][GridPoint]=0.0;
			}
		}
	}
	m_ContWorld->NIndexing->CalcDiffBoarder();
	pnpPrint0("</PNPUtil::SetCtoZeroWhereDZero>\n");
	return EXIT_SUCCESS;
}
int PNPUtil::SetDZeroWithIndex(ContWorld* _ContWorld)
{
	int i, j;
	pnpPrint0("<PNPUtil::SetDZeroWithIndex>\n");
	for (i = 0; i<_ContWorld->NIonsTypes; i++)
	{
		int count = 0;
		for (j = 0; j<_ContWorld->GS_XYZ; j++)
		{
			if (_ContWorld->NIndexing->GetDiffFloat(i, j) == 0.0 && _ContWorld->D[i][j] != 0.0)
			{
				_ContWorld->D[i][j] = 0.0;
				count++;
			}
		}
		pnpPrint0("%d point removed for ions %d\n", count, i);
	}
	pnpPrint0("</PNPUtil::SetDZeroWithIndex>\n");
	return EXIT_SUCCESS;
}
int PNPUtil::SetDzeroAtEps(ContWorld *m_ContWorld,int iEps)
{
	pnpPrint0("<PNPUtil::SetDzeroAtEps iEps=\"%d\">\n",iEps);
	
	int GS_X=m_ContWorld->GridSize[0];
	int GS_Y=m_ContWorld->GridSize[1];
	int GS_Z=m_ContWorld->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int i,j,k,ion;
	int GridPoint;
	
	NodeIndexing* NI=m_ContWorld->NIndexing;
	float GridScale=m_ContWorld->GridScale;
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS/(GridScale*GridScale*GridScale);
	
	
	iEps=iEps-1;
	pnpPrint0("Eps[iEps]=%g\n",NI->Eps[iEps]*EPKT);
	for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
	{
		int count=0;
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(NI->GetDiffFloat(ion,GridPoint)!=0.0)
			{
				if(NI->GetDiel(0,GridPoint)==iEps 
							|| NI->GetDiel(0,GridPoint-1)==iEps
							|| NI->GetDiel(1,GridPoint)==iEps
							|| NI->GetDiel(1,GridPoint-GS_X)==iEps
							|| NI->GetDiel(2,GridPoint)==iEps
							|| NI->GetDiel(2,GridPoint-GS_XY)==iEps)
				{
					NI->SetDiffToZero(ion,GridPoint);
					count++;
				}
			}
		}
		pnpPrint0("for ion %d removed %d nodes\n",ion,count);
	}
	m_ContWorld->NIndexing->CalcDiffBoarder();
	pnpPrint0("</PNPUtil::SetDzeroAtEps>\n");
	return EXIT_SUCCESS;
}
int PNPUtil::ScaleDiffusionInTheChannel(ContWorld *world,float x0,float y0,float* Rlim,float *DiffScale)
{
	pnpPrint0("PNPUtil::RemoveLargePMFfromPNP\n");
	//Rlim=Rlim*world->GridScale;
	x0=world->ConvFloatToGlobIntUnitsX(x0);
	y0=world->ConvFloatToGlobIntUnitsX(y0);
	
	int GS_X=world->GridSize[0];
	int GS_Y=world->GridSize[1];
	int GS_Z=world->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int i,j,k,ion;
	int GridPoint;
	if(world->D==NULL)
	{
		pnpPrint0("There is no world->D will create initial diffusion from world->NIndexing\n");
		PNP_EXIT_FAIL_NULL(world->NIndexing,"There is no world->NIndexing\n");
		world->D=new float*[world->NIonsTypes];
		if(world->NIonsTypes>=1)
			world->D[0]=world->NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion0);
		if(world->NIonsTypes>=2)
			world->D[1]=world->NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion1);
		if(world->NIonsTypes>=3)
			world->D[2]=world->NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion2);
		if(world->NIonsTypes>=4)
			world->D[3]=world->NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion3);
	}
	float RSQ,R;
	float RlimSQ;
	for(ion=0;ion<world->NIonsTypes;ion++)
	{
		//float* D=NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion0);
		for(i=0;i<GS_X;i++)
			for(j=0;j<GS_Y;j++)
				for(k=0;k<GS_Z;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			RSQ=(i-x0)*(i-x0)+(j-y0)*(j-y0);
			RlimSQ=Rlim[k]*Rlim[k]*world->GridScale*world->GridScale;
			if(RSQ<=RlimSQ&&world->D[ion][GridPoint]>0.0)
			{
				world->D[ion][GridPoint]*=DiffScale[k];
			}
		}
	}
	world->NIndexing->CalcDiffBoarder();
	
	return EXIT_SUCCESS;
}
int PNPUtil::RemoveQfromNI(ContWorld *m_ContWorld)
{
	pnpPrint0("<NPMaskBuilder::RemoveQfromNI>\n");
	m_ContWorld->NIndexing->SetChargeMapToZero();
	pnpPrint0("</NPMaskBuilder::RemoveQfromNI>\n");
	return EXIT_SUCCESS;
}
int PNPUtil::ShouldRemoveNodeBasedOnNPCriteria(ContWorld *m_ContWorld,int ion, int GrdPnt,float Relaxation,float MaxdC, bool RemoveNegC)
{
	float **C=m_ContWorld->C;
	float *Potential=m_ContWorld->Potential;
	float **VLJ=m_ContWorld->PMF;
	float *IonsQ=m_ContWorld->IonsQ;
	
	int GS_X=m_ContWorld->GridSize[0];
	int GS_Y=m_ContWorld->GridSize[1];
	int GS_Z=m_ContWorld->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	NodeIndexing *NIndexing=m_ContWorld->NIndexing;
	
	float dix[6],dixt;
	dix[0] = NIndexing->GetDiffFloat(ion,GrdPnt+1)>0.0?
			0.5*(NIndexing->GetDiffFloat(ion,GrdPnt)+NIndexing->GetDiffFloat(ion,GrdPnt+1)):0.0;
	dix[1] = NIndexing->GetDiffFloat(ion,GrdPnt-1)>0.0?
			0.5*(NIndexing->GetDiffFloat(ion,GrdPnt)+NIndexing->GetDiffFloat(ion,GrdPnt-1)):0.0;
	dix[2] = NIndexing->GetDiffFloat(ion,GrdPnt+GS_X)>0.0?
			0.5*(NIndexing->GetDiffFloat(ion,GrdPnt)+NIndexing->GetDiffFloat(ion,GrdPnt+GS_X)):0.0;
	dix[3] = NIndexing->GetDiffFloat(ion,GrdPnt-GS_X)>0.0?
			0.5*(NIndexing->GetDiffFloat(ion,GrdPnt)+NIndexing->GetDiffFloat(ion,GrdPnt-GS_X)):0.0;
	dix[4] = NIndexing->GetDiffFloat(ion,GrdPnt+GS_XY)>0.0?
			0.5*(NIndexing->GetDiffFloat(ion,GrdPnt)+NIndexing->GetDiffFloat(ion,GrdPnt+GS_XY)):0.0;
	dix[5] = NIndexing->GetDiffFloat(ion,GrdPnt-GS_XY)>0.0?
			0.5*(NIndexing->GetDiffFloat(ion,GrdPnt)+NIndexing->GetDiffFloat(ion,GrdPnt-GS_XY)):0.0;
	dix[0] *= C[ion][GrdPnt+1];
	dix[1] *= C[ion][GrdPnt-1];
	dix[2] *= C[ion][GrdPnt+GS_X];
	dix[3] *= C[ion][GrdPnt-GS_X];
	dix[4] *= C[ion][GrdPnt+GS_XY];
	dix[5] *= C[ion][GrdPnt-GS_XY];
	dixt = dix[0]+dix[1]+dix[2]+dix[3]+dix[4]+dix[5];
	
	float psi[7];
	psi[0] = 0.5*(IonsQ[ion]*Potential[GrdPnt]+VLJ[ion][GrdPnt]);
	psi[1] = 0.5*(IonsQ[ion]*Potential[GrdPnt+1]+VLJ[ion][GrdPnt+1]);
	psi[2] = 0.5*(IonsQ[ion]*Potential[GrdPnt-1]+VLJ[ion][GrdPnt-1]);
	psi[3] = 0.5*(IonsQ[ion]*Potential[GrdPnt+GS_X]+VLJ[ion][GrdPnt+GS_X]);
	psi[4] = 0.5*(IonsQ[ion]*Potential[GrdPnt-GS_X]+VLJ[ion][GrdPnt-GS_X]);
	psi[5] = 0.5*(IonsQ[ion]*Potential[GrdPnt+GS_XY]+VLJ[ion][GrdPnt+GS_XY]);
	psi[6] = 0.5*(IonsQ[ion]*Potential[GrdPnt-GS_XY]+VLJ[ion][GrdPnt-GS_XY]);
	
	
	float fpoh = 4*M_PI*m_ContWorld->GridScale;
	float conv = (m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale)/COANGS;
	float coef=fpoh*COANGS/(m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
	float coef2=COANGS/(m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
	float cloc[7];
	cloc[0] =C[ion][GrdPnt]/fpoh;
	cloc[1] =C[ion][GrdPnt+1]/fpoh;
	cloc[2] =C[ion][GrdPnt-1]/fpoh;
	cloc[3] =C[ion][GrdPnt+GS_X]/fpoh;
	cloc[4] =C[ion][GrdPnt-GS_X]/fpoh;
	cloc[5] =C[ion][GrdPnt+GS_XY]/fpoh;
	cloc[6] =C[ion][GrdPnt-GS_XY]/fpoh;
	
	float dphi[6];
	dphi[0] = dix[0]*psi[1];
	dphi[1] = dix[1]*psi[2];
	dphi[2] = dix[2]*psi[3];
	dphi[3] = dix[3]*psi[4];
	dphi[4] = dix[4]*psi[5];
	dphi[5] = dix[5]*psi[6];
	
	float phit = dphi[0] + dphi[1] + dphi[2] + dphi[3] + dphi[4] + dphi[5] - dixt*psi[0];
	float TMP = 1.0/(dixt-phit);
	
	float rodlt1 = cloc[1]*dix[0] + cloc[2]*dix[1] + cloc[3]*dix[2] + cloc[4]*dix[3] + cloc[5]*dix[4] + cloc[6]*dix[5];
	float rophi1 = cloc[1]*dphi[0] + cloc[2]*dphi[1] + cloc[3]*dphi[2] + cloc[4]*dphi[3] + cloc[5]*dphi[4] + cloc[6]*dphi[5];
	float term1 = rodlt1*(1.0-psi[0])+rophi1;
	
	float change = Relaxation*(term1*TMP-cloc[0]);
	float newC=cloc[0]+change;
	
	//float PsiSum=dix[0]*V[ion][GrdPnt+1] + dix[1]*V[ion][GrdPnt-1] + dix[2]*V[ion][GrdPnt+GS_X] + dix[3]*V[ion][GrdPnt-GS_X] + dix[4]*V[ion][GrdPnt+GS_XY] + dix[5]*V[ion][GrdPnt-GS_XY];
	//float Criteria=2.0*dixt+PsiSum-dixt*V[ion][GrdPnt];
	//if(newC<=0.0||fabs(change)/coef2>0.0010)
	if(MaxdC>0.0)
	{
		if(newC<=0.0 && RemoveNegC)
			return 1;
		if(fabs(change)>MaxdC/fpoh)
			return 2;
		else
			return 0;
	}
	else
	{
		if(newC<=0.0)
			return 1;
		else
			return 0;
	}
}
int PNPUtil::RemoveNodesFromPNPBasedOnNPCriteria(ContWorld *m_ContWorld,float Relaxation, int MaxCycles,float MaxdC, bool RemoveNegC)
{
	pnpPrint0("<PNPUtil::RemoveNodesFromPNPBasedOnNPCriteria>\n");
	float dPMF=2, PMF=10;
	//pnpPrint0("will remove points where dPMF>=%g and PMF>=%g\n",dPMF,PMF);
	
	int i,j,k,ion;
	int GridPoint;
	double dPMFx;
	float **V=m_ContWorld->PMF;
	int GS_X=m_ContWorld->GridSize[0];
	int GS_Y=m_ContWorld->GridSize[1];
	int GS_Z=m_ContWorld->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int TotalPointsRemoved=0;
	
	int *RemoveNode=new int[GS_XYZ];
	for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
	{
		int TotalPointsRemovedForThisIon=0;
		int iCycles;
		for(iCycles=0;iCycles<MaxCycles;iCycles++)
		{
			for(i=0;i<GS_XYZ;i++)
				RemoveNode[i]=0;
			
			for(i=1;i<GS_X-1;i++)
				for(j=1;j<GS_Y-1;j++)
					for(k=1;k<GS_Z-1;k++)
			{
				GridPoint=i+j*GS_X+k*GS_XY;
				if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint)!=0.0&& ( m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-1)==0.0 || m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+1)==0.0 || m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_X)==0.0 || m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_X)==0.0 || m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_XY)==0.0 || m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_XY)==0.0))
				{
					RemoveNode[GridPoint]=ShouldRemoveNodeBasedOnNPCriteria(m_ContWorld, ion,GridPoint, Relaxation,MaxdC,RemoveNegC);
					/*if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-1)!=0.0)
					{//boarder at -x
						if(ShouldRemoveNodeBasedOnNPCriteria(m_ContWorld, ion,GridPoint-1, Relaxation,MaxdC))
						{
							RemoveNode[GridPoint-1]++;
							//m_ContWorld->NIndexing->SetDiffToZero(ion,);
							//Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+1)!=0.0)
					{//boarder at +x
						if(ShouldRemoveNodeBasedOnNPCriteria(m_ContWorld, ion,GridPoint+1, Relaxation,MaxdC))
						{
							RemoveNode[GridPoint+1]++;
							//m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint+1);
							//Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_X)!=0.0)
					{//boarder at -y
						if(ShouldRemoveNodeBasedOnNPCriteria(m_ContWorld, ion,GridPoint-GS_X, Relaxation,MaxdC))
						{
							RemoveNode[GridPoint-GS_X]++;
							//m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint-GS_X);
							//Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_X)!=0.0)
					{//boarder at +y
						if(ShouldRemoveNodeBasedOnNPCriteria(m_ContWorld, ion,GridPoint+GS_X, Relaxation,MaxdC))
						{
							RemoveNode[GridPoint+GS_X]++;
							//m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint+GS_X);
							//Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_XY)!=0.0)
					{//boarder at -z
						if(ShouldRemoveNodeBasedOnNPCriteria(m_ContWorld, ion,GridPoint-GS_XY, Relaxation,MaxdC))
						{
							RemoveNode[GridPoint-GS_XY]++;
							//m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint-GS_XY);
							//Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_XY)!=0.0)
					{//boarder at +z
						if(ShouldRemoveNodeBasedOnNPCriteria(m_ContWorld, ion,GridPoint+GS_XY, Relaxation,MaxdC))
						{
							RemoveNode[GridPoint+GS_XY]++;
							//m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint+GS_XY);
							//Count++;
						}
					}*/
				}
			}
			int Count;
			int CountNegC=0,CountBigdC=0;
			for(i=0;i<GS_XYZ;i++)
			{
				if(RemoveNode[i]==1)
				{
					m_ContWorld->NIndexing->SetDiffToZero(ion,i);
					CountNegC++;
				}
				else if(RemoveNode[i]==2)
				{
					m_ContWorld->NIndexing->SetDiffToZero(ion,i);
					CountBigdC++;
				}
			}
			Count=CountBigdC+CountNegC;
			pnpPrint0("Removed %d nodes(NegC) and %d (BigdC), for ion %d\n",CountNegC,CountBigdC,ion);
			TotalPointsRemovedForThisIon+=Count;
			TotalPointsRemoved+=Count;
			if(Count==0)break;
		}
		pnpPrint0("Totally removed %d nodes, for ion %d\n",TotalPointsRemovedForThisIon,ion);
	}
	delete [] RemoveNode;
	m_ContWorld->NIndexing->CalcDiffBoarder();
	pnpPrint0("</PNPUtil::RemoveNodesFromPNPBasedOnNPCriteria>\n");
	return TotalPointsRemoved;
}
/*bool PNPUtil::ProcessPNPUtilCmds(ContWorld *world,const TiXmlElement *Elt)
{
	if(m_ContWorld!=NULL)
	{
		Status=m_ContWorld->RemoveDiffusionPointsAtNegativeC();
		if(Status==EXIT_FAILURE)return EXIT_FAILURE;
	}
	else
	{
		pnpError("m_ContWorld is not initialized\n");
		return false;
	}
}*/
/////////////////////////////////////////////////////////////////////
NPMaskBuilder::NPMaskBuilder(ContWorld *world)
{
	m_ContWorld=world;
	NPMask=NULL;
}
NPMaskBuilder::~NPMaskBuilder()
{
	DeleteObjByPnt(NPMask);
}
int NPMaskBuilder::InitNPMask()
{
	NPMask=new VectorField3D(m_ContWorld->GridSize,m_ContWorld->GridScale,m_ContWorld->NIonsTypes);
	NPMask->FillValue(0.0);
	return EXIT_SUCCESS;
}
int NPMaskBuilder::SetToNIDiffusion()
{
	pnpPrint0("<NPMaskBuilder::SetToNIDiffusion>\n");
	int GS_X=m_ContWorld->GridSize[0];
	int GS_Y=m_ContWorld->GridSize[1];
	int GS_Z=m_ContWorld->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	NodeIndexing *NIndexing=m_ContWorld->NIndexing;
	
	int i,j,k;
	int ion, GrdPnt;
	
	for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
	{
		for(i=0;i<GS_X;i++)
			for(j=0;j<GS_Y;j++)
				for(k=0;k<GS_Z;k++)
		{
			GrdPnt=i+j*GS_X+k*GS_XY;
			if(NIndexing->GetDiffFloat(ion,GrdPnt)!=0.0)
			{
				NPMask->V[ion][GrdPnt]=1.0;
			}
		}
	}
	pnpPrint0("</NPMaskBuilder::SetToNIDiffusion>\n");
	return EXIT_SUCCESS;
}
int NPMaskBuilder::RemoveTubeRegion(int ion,float X,float Y,float Z0, float Z1,float R0,float R1)
{
	pnpPrint0("<NPMaskBuilder::RemoveTubeRegion>\n");
	pnpPrint0("Input:\n");
	pnpPrint0("\tIon= %d\n",ion);
	pnpPrint0("\tX  = %g [int] = %g [A]\n", X, m_ContWorld->ConvGlobIntToGlobExtUnitsX(X));
	pnpPrint0("\tY  = %g [int] = %g [A]\n", Y, m_ContWorld->ConvGlobIntToGlobExtUnitsY(Y));
	pnpPrint0("\tZ0 = %g [int] = %g [A]\n", Z0, m_ContWorld->ConvGlobIntToGlobExtUnitsZ(Z0));
	pnpPrint0("\tZ1 = %g [int] = %g [A]\n", Z1, m_ContWorld->ConvGlobIntToGlobExtUnitsZ(Z1));
	pnpPrint0("\tR0 = %g [int] = %g [A]\n", R0, R0/m_ContWorld->GridScale);
	pnpPrint0("\tR1 = %g [int] = %g [A]\n", R1, R1/m_ContWorld->GridScale);
	
	float R0SQ=R0*R0;
	float R1SQ=R1*R1;
	
	int GS_X=m_ContWorld->GridSizeGlobal[0];
	int GS_Y=m_ContWorld->GridSizeGlobal[1];
	int GS_Z=m_ContWorld->GridSizeGlobal[2];
	int GS_XY=GS_X*GS_Y;
	
	int start[3],end[3];
	
	start[0]=(int)(X-R1-1.0);
	end[0]=(int)(X+R1+1.0);
	start[1]=(int)(Y-R1-1.0);
	end[1]=(int)(Y+R1+1.0);
	start[2]=(int)(Z0-1.0);
	end[2]=(int)(Z1+1.0);
	
	if(start[0]<0)start[0]=0;
	if(end[0]>GS_X-1)end[0]=GS_X-1;
	if(start[1]<0)start[1]=0;
	if(end[1]>GS_Y-1)end[1]=GS_Y-1;
	if(start[2]<0)start[2]=0;
	if(end[2]>GS_Z-1)end[2]=GS_Z-1;
	
	int GrdPnt;
	int ix,iy,iz;
	float RSQ;
	int count=0;
	
	for(ix=start[0];ix<=end[0];ix++)
		for(iy=start[1];iy<=end[1];iy++)
			for(iz=start[2];iz<=end[2];iz++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		RSQ=(ix-X)*(ix-X)+(iy-Y)*(iy-Y);
		
		if(R0SQ<=RSQ&&RSQ<=R1SQ&&Z0<=iz&&iz<=Z1)
		{
			if(NPMask->V[ion][GrdPnt]>0.0)
			{
				NPMask->V[ion][GrdPnt]=0.0;
				count++;
			}
		}
	}
	pnpPrint0("Output:\n");
	pnpPrint0("\t%d nodes was removed\n", count);
	pnpPrint0("</NPMaskBuilder::RemoveTubeRegion>\n");
	return EXIT_SUCCESS;
}
int NPMaskBuilder::RemoveHighC(float HighC)
{
	pnpPrint0("<NPMaskBuilder::RemoveHighC>\n");
	int GS_X=m_ContWorld->GridSizeGlobal[0];
	int GS_Y=m_ContWorld->GridSizeGlobal[1];
	int GS_Z=m_ContWorld->GridSizeGlobal[2];
	int GS_XY=GS_X*GS_Y;
	int GrdPnt;
	int ion;
	int ix,iy,iz;
	float RSQ;
	int count=0;
	
	float fpoh= 4*M_PI*m_ContWorld->GridScale;
	float coef=fpoh*COANGS/(m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
	
	pnpPrint0("Input:\n");
	pnpPrint0("\tHighC= %f\n",HighC/coef);
	
	for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
	{
		for(ix=2;ix<=GS_X-3;ix++)
			for(iy=2;iy<=GS_Y-3;iy++)
				for(iz=2;iz<=GS_Z-3;iz++)
		{
			GrdPnt=ix+iy*GS_X+iz*GS_XY;
			if(m_ContWorld->C[ion][GrdPnt]>=HighC&&NPMask->V[ion][GrdPnt]>0.0)
			{
				NPMask->V[ion][GrdPnt]=0.0;
				NPMask->V[ion][GrdPnt+1]=0.0;
				NPMask->V[ion][GrdPnt-1]=0.0;
				NPMask->V[ion][GrdPnt+GS_X]=0.0;
				NPMask->V[ion][GrdPnt-GS_X]=0.0;
				NPMask->V[ion][GrdPnt+GS_XY]=0.0;
				NPMask->V[ion][GrdPnt-GS_XY]=0.0;
				
				count++;
			}
		}
	}
	pnpPrint0("Output:\n");
	pnpPrint0("\t%d points removed\n",count);
	pnpPrint0("</NPMaskBuilder::RemoveHighC>\n");
	return EXIT_SUCCESS;
}
int NPMaskBuilder::ReadNPMask(const char * filename)
{
	pnpPrint0("<NPMaskBuilder::ReadNPMask>\n");
	NPMask->ReadFromFile(filename,1.0);
	pnpPrint0("</NPMaskBuilder::ReadNPMask>\n");
	return EXIT_SUCCESS;
}
int NPMaskBuilder::WriteNPMask(const char * filename)
{
	pnpPrint0("<NPMaskBuilder::WriteNPMask>\n");
	NPMask->WriteToFile(filename,1.0);
	pnpPrint0("</NPMaskBuilder::WriteNPMask>\n");
	return EXIT_SUCCESS;
}
int NPMaskBuilder::CmdNPMaskBuilder(const TiXmlElement *Elt)
{
	PNP_EXIT_FAIL_NULL(m_ContWorld,"NPMaskBuilder::CmdNPMaskBuilder m_ContWorld is not initialized\n");
	pnpPrint0("<NPMaskBuilder::CmdNPMaskBuilder>\n");
	InitNPMask();
	
	const TiXmlElement *CldElt;
	char ctmp[256];
	CldElt=Elt->FirstChildElement();
	do
	{
		if(strcmp("SetToNIDiffusion",CldElt->Value())==0)
			SetToNIDiffusion();
		else if(strcmp("ReadNPMask",CldElt->Value())==0)
		{
			ReadNPMask(CldElt->Attribute("filename"));
		}
		else if(strcmp("WriteNPMask",CldElt->Value())==0)
		{
			WriteNPMask(CldElt->Attribute("filename"));
		}
		else if(strcmp("RemoveTubeRegion",CldElt->Value())==0)
		{
			int ion;
			float XY[2], Z[2], R[2];
			CldElt->GetIntAttribute("Ion",&ion);
			CldElt->GetArrOfFloatAttribute("XY",XY,2);
			CldElt->GetArrOfFloatAttribute("Z",Z,2);
			CldElt->GetArrOfFloatAttribute("R",R,2);
			
			if(ion<0||ion>10||R[0]<0.0||R[1]<0.0)
			{
				pnpError("RemoveTubeRegion something wrong with parameters\n");
				pnpError("Ion=%d XY=[%g %g][A A] Z=[%g %g][A A] R=[%g %g][A A]\n", ion, XY[0],XY[1], Z[0], Z[1], R[0], R[1]);
			}
			else
			{
				XY[0]=m_ContWorld->ConvFloatToGlobIntUnitsX(XY[0]);
				XY[1]=m_ContWorld->ConvFloatToGlobIntUnitsY(XY[1]);
				Z[0]=m_ContWorld->ConvFloatToGlobIntUnitsZ(Z[0]);
				Z[1]=m_ContWorld->ConvFloatToGlobIntUnitsZ(Z[1]);
				R[0]*=m_ContWorld->GridScale;
				R[1]*=m_ContWorld->GridScale;
				
				RemoveTubeRegion(ion, XY[0],XY[1], Z[0], Z[1], R[0], R[1]);
			}
		}
		else if(strcmp("RemoveHighC",CldElt->Value())==0)
		{
			float HighC;
			float fpoh= 4*M_PI*m_ContWorld->GridScale;
			float coef=fpoh*COANGS/(m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
			CldElt->GetFloatAttribute("HighC",&HighC);
			RemoveHighC(coef*HighC);
		}
		CldElt=CldElt->NextSiblingElement();
	}
	while(CldElt!=NULL);
	pnpPrint0("</NPMaskBuilder::CmdNPMaskBuilder>\n");
	return EXIT_SUCCESS;
}
/////////////////////////////////////////////////////////////////////
IAVCalc::IAVCalc(ContWorld *world)
{
	m_ContWorld=world;
	PBLJ=NULL;
	NP=NULL;
	//NPMask=NULL;
}
IAVCalc::~IAVCalc()
{
	DeleteObjByPnt(PBLJ);
	DeleteObjByPnt(NP);
}
int pnpsIsBadNum(float v)
{
	#if defined(_MSC_VER)
	if(_finite(v))
	#else
	if(isfinite(v))
	#endif
	{
		if(v<1.0e9&&v>-1.0e9)
			return 0;
		else
			return 1;
	}
	else
		return 1;
}
int IAVCalc::CmdIAVCalcDict(PyObject *dict)
{
	PNP_EXIT_FAIL_NULL(m_ContWorld,"IAVCalc::CmdIAVCalc m_ContWorld is not initialized\n");
	pnpPrint0("<IAVCalc::CmdIAVCalc>\n");
	
	
	//const TiXmlElement *CldElt;
	//char ctmp[256];
	
	//const TiXmlElement *PBLJElt=NULL;
	//const TiXmlElement *NPElt=NULL;
	
	
	int i;
	int ion,GrdPnt;
	int GS_X=m_ContWorld->GridSizeGlobal[0];
	int GS_Y=m_ContWorld->GridSizeGlobal[1];
	int GS_Z=m_ContWorld->GridSizeGlobal[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_X*GS_Y*GS_Z;
	
	
	//Read Parameters for IAVCalc
	
	
	VectorField3D *VCpres=new VectorField3D(m_ContWorld->GridSize,m_ContWorld->GridScale,m_ContWorld->NIonsTypes);
	float **Cpres=VCpres->V;
	
	float MaxdC=haPyDict_GetItemValueAsFloat(dict,"MaxdC",10000.0);
	
	float fpoh= 4*M_PI*m_ContWorld->GridScale;
	float coef=fpoh*COANGS/(m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
	pnpPrint0("MaxdC=%g M\n");
	MaxdC*=coef;
	pnpPrint0("MaxdC=%g [int]\n");
	int MaxCycle=haPyDict_GetItemValueAsInt(dict,"MaxCycles",1000);

	float Relaxation;

	bool RunPBLJwNPCriteriaCycle=haPyDict_GetItemValueAsBool(dict,"RunPBLJwNPCriteriaCycle",false);
	Relaxation=haPyDict_GetItemValueAsFloat(dict,"RelaxationNPCriteria",1.0);
	int Cycles4onePBLJ=haPyDict_GetItemValueAsInt(dict,"Cycles4onePBLJ",1);
	bool RemoveNegC=haPyDict_GetItemValueAsBool(dict,"RemoveNegC",true);
	
	if(RunPBLJwNPCriteriaCycle)
	{
		PyObject *PBLJdict=PyDict_GetItemString(dict,"PBSR_Param");
		PNP_EXIT_FAIL_NULL(PBLJdict,"IAVCalc::CmdIAVCalc can not find PBwithLJSolver description in input file\n");
		int NodesRemoved;
		int Cycles=0;
		do
		{
			
			//PBLJ
			PBLJ=new PBwithLJSolver();
			PBLJ->LoadParamFromPyDict(PBLJdict);
			PBLJ->SetContWorld(m_ContWorld);
			PBLJ->InitSolver();
			PBLJ->Solve();
			DeleteObjByPnt(PBLJ);

			PNPUtil::ConvertPBLJresultsToDynamicCharge(m_ContWorld);
			
			//Remove Bad nodes
			NodesRemoved=0;
			NodesRemoved+=PNPUtil::RemoveNodesFromPNPBasedOnNPCriteria(m_ContWorld, Relaxation,Cycles4onePBLJ,MaxdC,RemoveNegC);
			m_ContWorld->NIndexing->RemoveBadDiffusionPoints();
			PNPUtil::RemoveCavitiesAtDiffusionMap(m_ContWorld);
			PNPUtil::SetCtoZeroWhereDZero(m_ContWorld);
			Cycles++;
			pnpPrint0("RunPBLJwNPCriteriaCycle:: Cycle=%d\n",Cycles);
			
			pnpPrint0("SaveTMPfiles:\n");
			//check values first
			int status=0;
			for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
			{
				pnpsIsBadNum(m_ContWorld->Potential[GrdPnt]);
				pnpsIsBadNum(m_ContWorld->C[0][GrdPnt]);
				pnpsIsBadNum(m_ContWorld->C[1][GrdPnt]);
			}
			
			if(status==0)
			{
				pnpPrint0("status=%d\n",status);
				m_ContWorld->WriteNodeIndexing("tmpIAVCalc_NI.gz");
				m_ContWorld->WritePotential("tmpIAVCalc_p.bin");
				m_ContWorld->WriteDiffusion("tmpIAVCalc_D.bin");
				m_ContWorld->WriteDynamicCharge("tmpIAVCalc_C.bin");
			}
			else
			{
				pnpPrint0("status=%d\n",status);
				pnpPrint0("Potential has bad numbers wouldn't save tmp files\n");
			}
		}
		while(NodesRemoved>0&&Cycles<MaxCycle);
		if(NodesRemoved>0&&Cycles<MaxCycle)
			pnpError("made 100 Cycles and nodes still removed\n");
	}
	
	bool RunDeepPBLJwNPCycle=haPyDict_GetItemValueAsBool(dict,"RunDeepPBLJwNPCycle",false);
	
	int RunSuccessfullCycles=haPyDict_GetItemValueAsInt(dict,"RunSuccessfullCycles",40);
	
	
	float dC;
	if(RunDeepPBLJwNPCycle)
	{
		pnpPrint0("<RunPBLJwNPCriteriaCycle>\n");
		PyObject *PBLJdict=PyDict_GetItemString(dict,"PBSR_Param");
		PNP_EXIT_FAIL_NULL(PBLJdict,"IAVCalc::CmdIAVCalc can not find PBwithLJSolver description in input file\n");
		PyObject *NPdict=PyDict_GetItemString(dict,"NP_Param");
		PNP_EXIT_FAIL_NULL(NPdict,"IAVCalc::CmdIAVCalc can not find NernstPlankSolver description in input file\n");
		//NP
		for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
			for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
				Cpres[ion][GrdPnt]=m_ContWorld->C[ion][GrdPnt];
		
		
		
		int NodesToRemove0=0,NodesToRemove1=0;
		int NegNodes0=0,NegNodes1=0;
		int SuccessfullCycles=0;
		
		int Cycle=0;
		do
		{
			if(NodesToRemove0+NodesToRemove1>0)
			{
				m_ContWorld->NIndexing->CalcDiffBoarder();
				m_ContWorld->NIndexing->RemoveBadDiffusionPoints();
				PNPUtil::RemoveCavitiesAtDiffusionMap(m_ContWorld);
				PNPUtil::SetDZeroWithIndex(m_ContWorld);
				PNPUtil::SetCtoZeroWhereDZero(m_ContWorld);
				
				PBLJ=new PBwithLJSolver();
				PBLJ->LoadParamFromPyDict(PBLJdict);
				PBLJ->SetContWorld(m_ContWorld);
				PBLJ->InitSolver();
				PBLJ->Solve();
				DeleteObjByPnt(PBLJ);

				PNPUtil::ConvertPBLJresultsToDynamicCharge(m_ContWorld);
				
				//Save new good C
				for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
					for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
						Cpres[ion][GrdPnt]=m_ContWorld->C[ion][GrdPnt];
				
				DeleteObjByPnt(NP);
				
			}
			if(NP==NULL)
			{
				NP=new NernstPlankSolver();
				NP->LoadParamFromPyDict(NPdict);
				NP->SetContWorld(m_ContWorld);
				NP->InitSolver();
			}
			
			//MakeHalfStep
			NodesToRemove0=0;
			NegNodes0=0;
			for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
			{
				NP->NernstPlanckSolverDpre(ion);
				NP->NernstPlanckSolverDiterHalfStep(ion,true,0);
				NP->NernstPlanckSolverDpost(ion);
				for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
				{
					dC=fabs(m_ContWorld->C[ion][GrdPnt]-Cpres[ion][GrdPnt]);
					if(m_ContWorld->C[ion][GrdPnt]<0.0||dC>MaxdC)
					{
						NegNodes0++;
					}
				}
				for(i=NP->SingularNum[ion][0];i<NP->SingularNum[ion][1];i++)
				{
					GrdPnt=NP->IndexSingular[ion][i];
					dC=fabs(m_ContWorld->C[ion][GrdPnt]-Cpres[ion][GrdPnt]);
					if(m_ContWorld->C[ion][GrdPnt]<0.0||dC>MaxdC)
					{
						if(NP->dix[ion][0][i]==0.0 || NP->dix[ion][1][i]==0.0 || NP->dix[ion][2][i]==0.0 || NP->dix[ion][3][i]==0.0 || NP->dix[ion][4][i]==0.0 || NP->dix[ion][5][i]==0.0)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GrdPnt);
							NodesToRemove0++;
						}
					}
				}
			}
			//if(NegNodes0>0)//restore old good C
			if(NodesToRemove0>0)
				for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
					for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
						m_ContWorld->C[ion][GrdPnt]=Cpres[ion][GrdPnt];
			
			NegNodes1=0;
			NodesToRemove1=0;
			for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
			{
				NP->NernstPlanckSolverDpre(ion);
				NP->NernstPlanckSolverDiterHalfStep(ion,true,1);
				NP->NernstPlanckSolverDpost(ion);
				for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
				{
					dC=fabs(m_ContWorld->C[ion][GrdPnt]-Cpres[ion][GrdPnt]);
					if(m_ContWorld->C[ion][GrdPnt]<0.0||dC>MaxdC)
					{
						NegNodes1++;
					}
				}
				for(i=NP->SingularNum[ion][1];i<NP->SingularNum[ion][2];i++)
				{
					GrdPnt=NP->IndexSingular[ion][i];
					dC=fabs(m_ContWorld->C[ion][GrdPnt]-Cpres[ion][GrdPnt]);
					if(m_ContWorld->C[ion][GrdPnt]<0.0||dC>MaxdC)
					{
						if(NP->dix[ion][0][i]==0.0 || NP->dix[ion][1][i]==0.0 || NP->dix[ion][2][i]==0.0 || NP->dix[ion][3][i]==0.0 || NP->dix[ion][4][i]==0.0 || NP->dix[ion][5][i]==0.0)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GrdPnt);
							NodesToRemove1++;
						}
					}
				}
			}
			//if(NegNodes0+NegNodes1>0)//restore old good C
			if(NodesToRemove0+NodesToRemove1>0)
				for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
					for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
						m_ContWorld->C[ion][GrdPnt]=Cpres[ion][GrdPnt];
			
			Cycle++;
			pnpPrint0("RunPBLJwNPCriteriaCycle:: Cycle=%d\n",Cycle);
			pnpPrint0("RunPBLJwNPCriteriaCycle:: NegNodes0+NegNodes1=%d\n",NegNodes0+NegNodes1);
			pnpPrint0("RunPBLJwNPCriteriaCycle:: NodesToRemove0+NodesToRemove1=%d\n", NodesToRemove0+NodesToRemove1);
			if(NodesToRemove0+NodesToRemove1==0)
			{
				SuccessfullCycles++;
				pnpPrint0("RunPBLJwNPCriteriaCycle:: Complite SuccessfullCycle number: %d\n",SuccessfullCycles);
				if(NegNodes0+NegNodes1>0)
					pnpPrint0("But still have negative C nodes %d\n",NegNodes0+NegNodes1);
			}
			else
			{
				SuccessfullCycles=0;
				pnpPrint0("RunPBLJwNPCriteriaCycle:: Will removed %d boarder points, total negative C points %d\n",NodesToRemove0+NodesToRemove1,NegNodes0+NegNodes1);
			}
			
		}
		while(SuccessfullCycles<RunSuccessfullCycles&&Cycle<MaxCycle);
		
		DeleteObjByPnt(NP);
		
		pnpPrint0("</RunPBLJwNPCriteriaCycle>\n");
		
	}
	/*CldElt=Elt->FirstChildElement();
	do
	{
		//if(strcmp("SetToNIDiffusion",CldElt->Value())==0)
		//	SetToNIDiffusion();
		//else if(strcmp("ReadNPMask",CldElt->Value())==0)
		//{
		//	ReadNPMask(CldElt->Attribute("filename"));
		//}
		CldElt=CldElt->NextSiblingElement();
	}
	while(CldElt!=NULL);*/
	DeleteObjByPnt(VCpres);
	pnpPrint0("</IAVCalc::CmdIAVCalc>\n");
	return EXIT_SUCCESS;
}
/////////////////////////////////////////////////////////////////////
int haPyDict_GetItemAsBool(PyObject *dict, const char *key, bool *v)
{
	PyObject *p=PyDict_GetItemString(dict,key);
	if(p!=NULL)
	{
		*v=(bool)PyInt_AsLong(p);
	}
	return EXIT_SUCCESS;
}
int haPyDict_GetItemAsInt(PyObject *dict, const char *key, int *v)
{
	PyObject *p=PyDict_GetItemString(dict,key);
	if(p!=NULL)
	{
		*v=(int)PyInt_AsLong(p);
	}
	return EXIT_SUCCESS;
}
int haPyDict_GetItemAsFloat(PyObject *dict, const char *key, float *v)
{
	PyObject *p=PyDict_GetItemString(dict,key);
	if(p!=NULL)
	{
		*v=(float)PyFloat_AsDouble(p);
	}
	return EXIT_SUCCESS;
}
int haPyDict_GetItemAsString(PyObject *dict, const char *key, std::string *v)
{
	PyObject *p=PyDict_GetItemString(dict,key);
	if(p!=NULL)
	{
		*v=PyString_AsString(p);
	}
	return EXIT_SUCCESS;
}

bool haPyDict_GetItemValueAsBool(PyObject *dict, const char *key, bool vdefault)
{
	bool v=vdefault;
	PyObject *p=PyDict_GetItemString(dict,key);
	if(p!=NULL)
	{
		v=(bool)PyInt_AsLong(p);
	}
	return v;
}
int haPyDict_GetItemValueAsInt(PyObject *dict, const char *key, int vdefault)
{
	int v=vdefault;
	PyObject *p=PyDict_GetItemString(dict,key);
	if(p!=NULL)
	{
		v=(int)PyInt_AsLong(p);
	}
	return v;
}
float haPyDict_GetItemValueAsFloat(PyObject *dict, const char *key, float vdefault)
{
	float v=vdefault;
	PyObject *p=PyDict_GetItemString(dict,key);
	if(p!=NULL)
	{
		v=(float)PyFloat_AsDouble(p);
	}
	return v;
}
char* haPyDict_GetItemValueAsString(PyObject *dict, const char *key, const char *vdefault)
{
	char *v;
	PyObject *p=PyDict_GetItemString(dict,key);
	if(p!=NULL)
	{
		v=PyString_AsString(p);
	}
	else
	{
		int i=strlen(vdefault);
		v=new char[i+1];
		strcpy(v,vdefault);
	}
	return v;
}
int haPy_SetCFloatArrFromListOfFloat(PyObject *vlist, float *carr)
{
	int n=PyList_Size(vlist);
	int i;
	for(i=0;i<n;i++)
	{
		carr[i]=(float)PyFloat_AsDouble(PyList_GetItem(vlist,i));
	}
	return EXIT_SUCCESS;
}