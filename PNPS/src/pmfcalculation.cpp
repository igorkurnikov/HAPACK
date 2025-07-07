//
// C++ Implementation: pmfcalculation
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
#include "mpi.h"
#endif

#include "pmfcalculation.h"
#include "contworld.h"
#include "tinyxml.h"
#include "pnpdebug.h"
#include "buildworldni.h"
#include "poissonsolver.h"
#include "pmfcalculation.h"
#include "math.h"
#include "pnpconstants.h"
#include "pnps.h"
#include "mapio.h"
#include <vector>


PMFCalculation * PMFCalculation0;

PMFCalculation::PMFCalculation()
{
	InitZero();
	PMFCalculation0=this;
}


PMFCalculation::~PMFCalculation()
{
	Clear();
}
int PMFCalculation::InitZero()
{
	m_ContWorld=NULL;
	m_ContWorldRough=NULL;
	//m_PoissonSolver=NULL;
	//m_BuildWorldNI=NULL;
	NIons=0;
	m_Surf=NULL;
	m_iVtmp=NULL;
	m_Qstat=NULL;
	
	Npnts=NULL;
	IonsForPMF=NULL;
	PMFPnts=NULL;
	PMF=NULL;
	PMFCalc=NULL;
	PMFxyz[0]=NULL;
	PMFxyz[1]=NULL;
	PMFxyz[2]=NULL;
	PMFPntsStart=NULL;
	
	GproteinLoad=false;
	GproteinQnum=0;
	GproteinQst=NULL;
	GproteinPhi=NULL;
	GproteinQndx=NULL;
	GproteinRoughLoad=false;
	GproteinRoughQnum=0;
	GproteinRoughQst=NULL;
	GproteinRoughPhi=NULL;
	GproteinRoughQndx=NULL;
	GionLoad=false;
	GionQnum=NULL;
	GionQst=NULL;
	GionPhi=NULL;
	GionQndx=NULL;
	grdpntQst=NULL;
	Gsubtract=NULL;
	
	WriteMapsForEachPosition=false;
	return EXIT_SUCCESS;
}
int PMFCalculation::Clear()
{
	DeleteObjByPnt(m_ContWorld);
	DeleteObjByPnt(m_ContWorldRough);
	DeleteCArray(m_Qstat);
	
	
	DeleteCArray(IonsForPMF);
	DeleteCVecArray(PMFPnts,NIons);
	DeleteCVecArray(PMF,NIons);
	DeleteCVecArray(PMFCalc,NIons);
	DeleteCArray(Npnts);

	DeleteCVecArray(m_Surf,3);
	DeleteCVecArray(m_iVtmp,3);
	
	DeleteCVecArray(PMFxyz[0],NIons);
	DeleteCVecArray(PMFxyz[1],NIons);
	DeleteCVecArray(PMFxyz[2],NIons);
	
	DeleteCArray(grdpntQst);
	DeleteCArray(Gsubtract);
	return EXIT_SUCCESS;
}
#if 0
int PMFCalculation::Run(const TiXmlElement *RootElt)
{
	const TiXmlElement *CldElt;
	
	CldElt=RootElt->FirstChildElement("IonsForPMF");
	PNP_EXIT_FAIL_NULL(CldElt,"ERROR: There is no <IonsForPMF> in <PMFCalculation>\n");
	PNP_EXIT_ON_FAIL(LoadIonsForPMF(CldElt));
	
	DbgPrint0("Load ContWorld parameters\n");
	CldElt=RootElt->FirstChildElement("ContWorld");
	PNP_EXIT_FAIL_NULL(CldElt,"There is no <ContWorld> in <PMFCalculation>\n");
	PNP_EXIT_FAIL_NULL((m_ContWorld=SingleWorldCalcController::CmdContWorld(CldElt)),
		"Can't initialize m_ContWorld\n");
	
	DbgPrint0("LoadGprotein if needed\n");
	if(RootElt->Attribute("LoadGprotein")!=NULL)
	{
		PNP_EXIT_ON_FAIL(LoadGprotein(RootElt->Attribute("LoadGprotein")));
	}
	if(RootElt->GetBoolAttribute("WriteMapsForEachPosition",&WriteMapsForEachPosition)==EXIT_FAILURE)
		WriteMapsForEachPosition=false;
	
	DbgPrint0("Load ions positions for calculation\n");
	CldElt=RootElt->FirstChildElement("PointsForPMF");
	PNP_EXIT_FAIL_NULL(CldElt,"ERROR: There is no <PointsForPMF> in <PMFCalculation>\n");
	if(pnpsapp->AmIBigBrother())
		PNP_EXIT_ON_FAIL(LoadPointsForPMF(CldElt));
	DbgPrint0("Check if something was already calculated\n");
	if(RootElt->Attribute("TempResultsOut")!=NULL)
	{
		sprintf(TempResultsOut,"%s\0",RootElt->Attribute("TempResultsOut"));
		//Check if something was already calculated
		if(pnpsapp->GetMyAbsRank()==pnpsapp->GetMaster())
		{
			bool OnePntPerNode=false;
			int count;
			
			RootElt->GetBoolAttribute("OnePntPerNode",&OnePntPerNode);
			if(OnePntPerNode)
				count=LoadCalculatedPointsOnePntPerNode(TempResultsOut);
			else
				count=LoadCalculatedPoints(TempResultsOut);
			if(count>0)
				RearrangePoints();
		}
	}
	else
	{
		sprintf(TempResultsOut,"TempResultsOut.dat\0");
	}
#ifdef MPI_PARALLEL
	char Temp[PNP_MAP_IO_STRING_LENGTH];
	sprintf(Temp,"%s\0",TempResultsOut);
	pnpsapp->AddMyGroupNumberToFileName(TempResultsOut,Temp);
#endif
	
	DbgPrint0("Distrebute points over processes\n");
	ReDistrebutePntsOverProcesses();
	
	if(RootElt->Attribute("Type")!=NULL)
	{
		if(strcmp("Simple",RootElt->Attribute("Type"))==0)
			return RunSimple(RootElt);
		else if(strcmp("Caching",RootElt->Attribute("Type"))==0)
			return RunCashing(RootElt);
		else if(strcmp("SimpleFocusing",RootElt->Attribute("Type"))==0)
			return RunSimpleFocusing(RootElt);
		else if(strcmp("NoIonSizeNoProteinCharge",RootElt->Attribute("Type"))==0)
			return RunNoIonSizeNoProteinCharge(RootElt);
		else
			return RunSimple(RootElt);
	}
	else
	{
		return RunSimple(RootElt);
	}
}
#endif
int PMFCalculation::RunNoIonSizeNoProteinCharge(const TiXmlElement *RootElt)
{
	pnpPrint("<PMFCalculation::RunNoIonSizeNoProteinCharge>\n");
	const TiXmlElement *EltContWorld=RootElt->FirstChildElement("ContWorld");
	const TiXmlElement *EltBuildWorld=RootElt->FirstChildElement("BuildWorldNI");
	const TiXmlElement *EltPoisson=RootElt->FirstChildElement("PoissonSolver");
	
	PNP_EXIT_FAIL_NULL(EltBuildWorld,"ERROR: There is no <BuildWorldNI> in <PMFCalculation>\n");
	PNP_EXIT_FAIL_NULL(EltPoisson,"ERROR: There is no <PoissonSolver> in <PMFCalculation>\n");
	
	int i,j,k;
	int iElt;
	
	BuildWorldNI *m_BuildWorldNI=CmdBuildWorld(EltBuildWorld,m_ContWorld);
	//remove Charge from protein
	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
	{
		if(m_BuildWorldNI->GOElms[iElt]->GetStdStrName()=="AtomsParameters")
		{
			GOAtoms* goatoms=(GOAtoms*)m_BuildWorldNI->GOElms[iElt];
			for(i=0;i<goatoms->q.size();i++)
				goatoms->q[i]=0.0;
		}
	}
	//Testing on first ion energy
	//m_BuildWorldNI->GOElms.resize(m_BuildWorldNI->GOElms.size()+1);
	//m_BuildWorldNI->GOElms[m_BuildWorldNI->GOElms.size()-1]=IonsForPMF[0];
	//SetIonCoor(0,0);
	//m_BuildWorldNI->BuildContWorld(m_ContWorld);
	//PoissonSolver* m_PoissonSolver=SingleWorldCalcController::CmdPoissonSolver(EltPoisson,m_ContWorld);
	
	//Create initial m_ContWorld
	m_BuildWorldNI->BuildContWorld(m_ContWorld);
	
	//Change Parameters for internal 
	m_BuildWorldNI->ChangeUnits(m_ContWorld,true);
	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
	{
		m_BuildWorldNI->GOElms[iElt]->ChangeUnits(m_ContWorld,true);
	}
	
	GOAtoms* goIon;
	
	
	//do Calculation
	int start,end,npnts,pnts4proc;
	for(i=0;i<NIons;i++)
	{
		
		goIon=IonsForPMF[i];
		
		m_BuildWorldNI->GOElms[m_BuildWorldNI->GOElms.size()-1]=goIon;
		
		npnts=Npnts[i]-PMFPntsStart[i];
		start=PMFPntsStart[i]+npnts*pnpsapp->GetMyGroupNumber()/pnpsapp->GetNumberOfGroups();
		end=PMFPntsStart[i]+npnts*(pnpsapp->GetMyGroupNumber()+1)/pnpsapp->GetNumberOfGroups();
		pnts4proc=end-start;
		if(pnpsapp->GetMyGroupNumber()==pnpsapp->GetNumberOfGroups()-1)end=PMFPntsStart[i]+npnts;
		DbgPrint0("SimplePMFCalculation: start=%d end=%d npnts=%d pnts4proc=%d\n", start, end, npnts, pnts4proc);
		
		DefClock0;
		
		for(j=start;j<end;j++)
		{
			pnpPrint0("PMFCalc[%d][0] %d %d\n",i,PMFCalc[i][j],PMFCalc[i][start]);
			if(PMFCalc[i][j]==false)
			{
				SetIonCoor(i,j);
				DbgPrint0("Start PMPPnts[%d][%d]=%d, r=[%f,%f,%f][Ext.coor.] R=%f A\n",i,j,PMFPnts[i][j],goIon->r[0][0],goIon->r[1][0],goIon->r[2][0],goIon->R[0]);
				goIon->ChangeUnits(m_ContWorld,true);
				DbgPrint0("      PMPPnts[%d][%d]=%d, r=[%f,%f,%f][int.coor.] R=%f grids\n",i,j,PMFPnts[i][j],goIon->r[0][0],goIon->r[1][0],goIon->r[2][0],goIon->R[0]);
				StartClock0;
				
				m_BuildWorldNI->BuildContWorldCharges(m_ContWorld,m_ContWorld->NIndexing);
				
				PoissonSolver* m_PoissonSolver;
				m_PoissonSolver=new PoissonSolver();
				m_PoissonSolver->WayToCalcSystemEnergy=1;
				m_PoissonSolver->LoadXML(EltPoisson);
				m_PoissonSolver->SetContWorld(m_ContWorld);
				m_PoissonSolver->InitSolver();
				m_PoissonSolver->Solve();
				
				PMFCalc[i][j]=true;
				PMF[i][j]=m_ContWorld->SystemEnergy;
				
				PrintSIPResultSingle(i,j,PMF[i][j],m_PoissonSolver->totalChange,m_PoissonSolver->ConvFac);
				
				DeleteObjByPnt(m_PoissonSolver);
				goIon->ChangeUnits(m_ContWorld,false);
				StopClockWMes0("one point calculation");
			}
		}
	}
	m_BuildWorldNI->GOElms.pop_back();
	
	//Change Parameters to external
	m_BuildWorldNI->ChangeUnits(m_ContWorld,false);
	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
	{
		m_BuildWorldNI->GOElms[iElt]->ChangeUnits(m_ContWorld,false);
	}
	fprintf(stdout,"RunNoIonSizeNoProteinCharge %d/%d DONE\n",pnpsapp->GetMyGroupNumber(),pnpsapp->GetNumberOfGroups());
	pnpPrint("</PMFCalculation::RunNoIonSizeNoProteinCharge>\n");
	return EXIT_SUCCESS;
}
int PMFCalculation::RunSimple(const TiXmlElement *RootElt)
{
	/*const TiXmlElement *EltContWorld=RootElt->FirstChildElement("ContWorld");
	const TiXmlElement *EltBuildWorld=RootElt->FirstChildElement("BuildWorldNI");
	const TiXmlElement *EltPoisson=RootElt->FirstChildElement("PoissonSolver");
	
	PNP_EXIT_FAIL_NULL(EltBuildWorld,"ERROR: There is no <BuildWorldNI> in <PMFCalculation>\n");
	PNP_EXIT_FAIL_NULL(EltPoisson,"ERROR: There is no <PoissonSolver> in <PMFCalculation>\n");
	
	m_BuildWorldNI=CmdBuildWorld(EltBuildWorld,m_ContWorld);
	
	//do Calculation
	int i,j,k;
	//Change Parameters for internal 
	int iElt;
	m_BuildWorldNI->ChangeUnits(m_ContWorld,true);
	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
	{
		m_BuildWorldNI->GOElms[iElt]->ChangeUnits(m_ContWorld,true);
	}
	
	GOAtoms* goIon;
	m_BuildWorldNI->GOElms.resize(m_BuildWorldNI->GOElms.size()+1);
	
	
	int start,end,npnts,pnts4proc;
	for(i=0;i<NIons;i++)
	{
		
		goIon=IonsForPMF[i];
		goIon->ChangeUnits(m_ContWorld,true);
		m_BuildWorldNI->GOElms[m_BuildWorldNI->GOElms.size()-1]=goIon;
		
		npnts=Npnts[i]-PMFPntsStart[i];
		start=PMFPntsStart[i]+npnts*pnpsapp->GetMyGroupNumber()/pnpsapp->GetNumberOfGroups();
		end=PMFPntsStart[i]+npnts*(pnpsapp->GetMyGroupNumber()+1)/pnpsapp->GetNumberOfGroups();
		pnts4proc=end-start;
		if(pnpsapp->GetMyGroupNumber()==pnpsapp->GetNumberOfGroups()-1)end=PMFPntsStart[i]+npnts;
		DbgPrint0("SimplePMFCalculation: start=%d end=%d npnts=%d pnts4proc=%d\n", start, end, npnts, pnts4proc);
		
		DefClock0;
		
		for(j=start;j<end;j++)
		{
			pnpPrint0("PMFCalc[%d][0] %d %d\n",i,PMFCalc[i][j],PMFCalc[i][start]);
			if(PMFCalc[i][j]==false)
			{
				SetIonCoor(i,j);
				fprintf(stdout,"Start PMPPnts[%d][%d]=%d, r=[%f,%f,%f][int.coor.]\n",i,j,PMFPnts[i][j],goIon->r[0][0],goIon->r[1][0],goIon->r[2][0]);
				
				StartClock0;
				
				BuildContWorld(m_ContWorld,m_BuildWorldNI);
				
				//init grdpntQst and Gsubtract
				DeleteCArray(grdpntQst);
				DeleteCArray(Gsubtract);
				grdpntQst=m_ContWorld->NIndexing->GetChargeIndex();
				if(GproteinLoad||GionLoad)
				{
					if(Gsubtract==NULL)
					{
						Gsubtract=new double[m_ContWorld->NIndexing->QNum];
						for(k=0;k<m_ContWorld->NIndexing->QNum;k++)
							Gsubtract[k]=0.0;
					}
		
					if(GproteinLoad)
					{
						if(GproteinQnum>m_ContWorld->NIndexing->QNum-1)
						{
							pnpWarning("Number of partial charges in Gprotein is bigger then partial charges in protein-ion - 1\n");
						}
						double qphi;
						int iGprotein;
						int count=0;
						int j0;
						for(j0=0;j0<GproteinQnum;j0++)
						{
							while(GproteinQndx[j0]!=grdpntQst[count])
							{
								count++;
							}
							qphi=0.0;
							iGprotein=-1;
							qphi=GproteinQst[j0]*GproteinPhi[j0];
							Gsubtract[count]=-qphi;
						}
					}
				}
				m_PoissonSolver=SingleWorldCalcController::CmdPoissonSolver(EltPoisson,m_ContWorld);
				
				PMFCalc[i][j]=true;
				PMF[i][j]=m_ContWorld->SystemEnergy;
				
				PrintSIPResultSingle(i,j,PMF[i][j],m_PoissonSolver->totalChange,m_PoissonSolver->ConvFac);
				
				DeleteObjByPnt(m_PoissonSolver);
				StopClockWMes0("one point calculation");
			}
		}
		goIon->ChangeUnits(m_ContWorld,false);
	}
	m_BuildWorldNI->GOElms.pop_back();
	
	//Change Parameters to external
	m_BuildWorldNI->ChangeUnits(m_ContWorld,false);
	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
	{
		m_BuildWorldNI->GOElms[iElt]->ChangeUnits(m_ContWorld,false);
	}
	fprintf(stdout,"SimplePMFCalculation %d/%d DONE\n",pnpsapp->GetMyGroupNumber(),pnpsapp->GetNumberOfGroups());*/
	return EXIT_SUCCESS;
}
int PMFCalculation::RunCashing(const TiXmlElement *RootElt)
{
	return EXIT_SUCCESS;
}
int PMFCalculation::RunSimpleFocusing(const TiXmlElement *RootElt)
{
	const TiXmlElement *EltContWorld=RootElt->FirstChildElement("ContWorld");
	const TiXmlElement *EltContWorldRough=RootElt->FirstChildElement("ContWorldRough");
	const TiXmlElement *EltBuildWorld=RootElt->FirstChildElement("BuildWorldNI");
	const TiXmlElement *EltPoissonRough=RootElt->FirstChildElement("PoissonSolverRough");
	const TiXmlElement *EltPoisson=RootElt->FirstChildElement("PoissonSolver");
	
	PNP_EXIT_FAIL_NULL(EltContWorld,"ERROR: There is no <ContWorld> in <PMFCalculation>\n");
	PNP_EXIT_FAIL_NULL(EltBuildWorld,"ERROR: There is no <BuildWorldNI> in <PMFCalculation>\n");
	PNP_EXIT_FAIL_NULL(EltPoisson,"ERROR: There is no <PoissonSolver> in <PMFCalculation>\n");
	
	PNP_EXIT_FAIL_NULL(EltContWorldRough,"ERROR: There is no <ContWorldRough> in <PMFCalculation>\n");
	PNP_EXIT_FAIL_NULL(EltPoissonRough,"ERROR: There is no <PoissonSolverRough> in <PMFCalculation>\n");
	
	m_ContWorldRough=SingleWorldCalcController::CmdContWorld(EltContWorldRough);
	DbgPrint0("LoadGproteinRough if needed\n");
	if(RootElt->Attribute("LoadGproteinRough")!=NULL)
	{
		PNP_EXIT_ON_FAIL(LoadGproteinRough(RootElt->Attribute("LoadGproteinRough")));
	}
	
	//output for rough results
	if(RootElt->Attribute("TempResultsOutRough")!=NULL)
		sprintf(TempResultsOutRough,"%s\0",RootElt->Attribute("TempResultsOutRough"));
	else
		sprintf(TempResultsOut,"TempResultsOutRough.dat\0");
#ifdef MPI_PARALLEL
	char Temp[PNP_MAP_IO_STRING_LENGTH];
	sprintf(Temp,"%s\0",TempResultsOutRough);
	pnpsapp->AddMyGroupNumberToFileName(TempResultsOutRough,Temp);
#endif
	
	BuildWorldNI* m_BuildWorld=CmdBuildWorld(EltBuildWorld,m_ContWorld);
	BuildWorldNI* m_BuildWorldRough=CmdBuildWorld(EltBuildWorld,m_ContWorldRough);
	PoissonSolver* m_PoissonSolver=NULL;
	
	VectorField3D * PotRough = new VectorField3D(m_ContWorldRough->GridSize, m_ContWorldRough->GridScale, 1);
	m_ContWorldRough->Potential=PotRough->V[0];
	
	VectorField3D * Pot = new VectorField3D(m_ContWorld->GridSize, m_ContWorld->GridScale, 1);
	m_ContWorld->Potential=Pot->V[0];
	//do Calculation
	int i,j,k;
	//Change Parameters for internal 
	int iElt;
	m_BuildWorld->ChangeUnits(m_ContWorld,true);
	m_BuildWorldRough->ChangeUnits(m_ContWorldRough,true);
	for(iElt=0;iElt<m_BuildWorld->GOElms.size();iElt++)
	{
		m_BuildWorld->GOElms[iElt]->ChangeUnits(m_ContWorld,true);
		m_BuildWorldRough->GOElms[iElt]->ChangeUnits(m_ContWorldRough,true);
	}
	
	GOAtoms* goIon=new GOAtoms();
	m_BuildWorld->GOElms.resize(m_BuildWorld->GOElms.size()+1);
	m_BuildWorldRough->GOElms.resize(m_BuildWorldRough->GOElms.size()+1);
	
	char filenameMapsOut[PNP_MAP_IO_STRING_LENGTH];
	int start,end,npnts,pnts4proc;
	for(i=0;i<NIons;i++)
	{
		m_BuildWorld->GOElms[m_BuildWorld->GOElms.size()-1]=goIon;
		m_BuildWorldRough->GOElms[m_BuildWorldRough->GOElms.size()-1]=goIon;
		
		npnts=Npnts[i]-PMFPntsStart[i];
		start=PMFPntsStart[i]+npnts*pnpsapp->GetMyGroupNumber()/pnpsapp->GetNumberOfGroups();
		end=PMFPntsStart[i]+npnts*(pnpsapp->GetMyGroupNumber()+1)/pnpsapp->GetNumberOfGroups();
		pnts4proc=end-start;
		if(pnpsapp->GetMyGroupNumber()==pnpsapp->GetNumberOfGroups()-1)
			end=PMFPntsStart[i]+npnts;
		DbgPrint0("SimplePMFCalculation: start=%d end=%d npnts=%d pnts4proc=%d\n", start, end, npnts, pnts4proc);
		
		DefClock0;
		
		for(j=start;j<end;j++)
		{
			pnpPrint0("PMFCalc[%d][0] %d %d\n",i,PMFCalc[i][j],PMFCalc[i][start]);
			if(PMFCalc[i][j]==false)
			{
				goIon->Copy(IonsForPMF[i]);
				SetIonCoor(goIon,i,j);
				DbgPrint0("Start PMPPnts[%d][%d]=%d, r=[%f,%f,%f][Ext.coor.] R=%f A\n",i,j,PMFPnts[i][j],goIon->r[0][0],goIon->r[1][0],goIon->r[2][0],goIon->R[0]);
				
				goIon->ChangeUnits(m_ContWorldRough,true);
				DbgPrint0("      PMPPnts[%d][%d]=%d, r=[%f,%f,%f][int.coor.] R=%f grids\n",i,j,PMFPnts[i][j],goIon->r[0][0],goIon->r[1][0],goIon->r[2][0],goIon->R[0]);
				
				StartClock0;
				
				BuildContWorld(m_ContWorldRough,m_BuildWorldRough);
				if(WriteMapsForEachPosition)
				{
					pnpsapp->AddProcessNumberToFileName(filenameMapsOut,"SIP_Rough_NI.gz","_i",i,NIons);
					pnpsapp->AddProcessNumberToFileName(filenameMapsOut,filenameMapsOut,"_p",j,Npnts[i]);
					m_ContWorldRough->WriteNodeIndexing(filenameMapsOut,1);
				}
				//init grdpntQst and Gsubtract
				DeleteCArray(grdpntQst);
				DeleteCArray(Gsubtract);
				grdpntQst=m_ContWorldRough->NIndexing->GetChargeIndex();
				if(GproteinRoughLoad)
				{
					Gsubtract=new double[m_ContWorldRough->NIndexing->QNum];
					for(k=0;k<m_ContWorldRough->NIndexing->QNum;k++)
						Gsubtract[k]=0.0;
					
					if(GproteinRoughQnum>m_ContWorldRough->NIndexing->QNum-1)
					{
						pnpWarning("Number of partial charges in Gprotein is bigger then partial charges in protein-ion - 1\n");
					}
					double qphi;
					int iGprotein;
					int count=0;
					int j0;
					for(j0=0;j0<GproteinRoughQnum;j0++)
					{
						while(GproteinRoughQndx[j0]!=grdpntQst[count])
						{
							count++;
						}
						qphi=0.0;
						iGprotein=-1;
						qphi=GproteinRoughQst[j0]*GproteinRoughPhi[j0];
						Gsubtract[count]=-qphi;
					}
				}
				m_PoissonSolver=SingleWorldCalcController::CmdPoissonSolver(EltPoissonRough,m_ContWorldRough);
				DeleteCArray(grdpntQst);
				DeleteCArray(Gsubtract);
				
				if(WriteMapsForEachPosition)
				{
					pnpsapp->AddProcessNumberToFileName(filenameMapsOut,"SIP_Rough_P.bin","_i",i,NIons);
					pnpsapp->AddProcessNumberToFileName(filenameMapsOut,filenameMapsOut,"_p",j,Npnts[i]);
					m_ContWorldRough->WritePotential(filenameMapsOut,1);
				}
				
				PrintSIPResultSingleRough(i,j,m_ContWorldRough->SystemEnergy,m_PoissonSolver->totalChange,m_PoissonSolver->ConvFac);
				
				DeleteObjByPnt(m_PoissonSolver);
				if(j==start)Pot->InterpolateFromExt(PotRough);
				else Pot->InterpolateBoarderFromExt(PotRough);
				//PMFCalculation0=this;
				goIon->Copy(IonsForPMF[i]);
				SetIonCoor(goIon,i,j);
				DbgPrint0("Start PMPPnts[%d][%d]=%d, r=[%f,%f,%f][Ext.coor.] R=%f A\n",i,j,PMFPnts[i][j],goIon->r[0][0],goIon->r[1][0],goIon->r[2][0],goIon->R[0]);
				
				goIon->ChangeUnits(m_ContWorld,true);
				DbgPrint0("      PMPPnts[%d][%d]=%d, r=[%f,%f,%f][int.coor.] R=%f grids\n",i,j,PMFPnts[i][j],goIon->r[0][0],goIon->r[1][0],goIon->r[2][0],goIon->R[0]);
				BuildContWorld(m_ContWorld,m_BuildWorld);
				if(WriteMapsForEachPosition)
				{
					pnpsapp->AddProcessNumberToFileName(filenameMapsOut,"SIP_NI.gz","_i",i,NIons);
					pnpsapp->AddProcessNumberToFileName(filenameMapsOut,filenameMapsOut,"_p",j,Npnts[i]);
					m_ContWorld->WriteNodeIndexing(filenameMapsOut,1);
				}
				//init grdpntQst and Gsubtract
				DeleteCArray(grdpntQst);
				DeleteCArray(Gsubtract);
				grdpntQst=m_ContWorld->NIndexing->GetChargeIndex();
				if(GproteinLoad)
				{
					Gsubtract=new double[m_ContWorld->NIndexing->QNum];
					for(k=0;k<m_ContWorld->NIndexing->QNum;k++)
						Gsubtract[k]=0.0;
					
					if(GproteinQnum>m_ContWorld->NIndexing->QNum-1)
					{
						pnpWarning("Number of partial charges in Gprotein is bigger then partial charges in protein-ion - 1\n");
					}
					double qphi;
					int iGprotein;
					int count=0;
					int j0;
					for(j0=0;j0<GproteinQnum;j0++)
					{
						while(GproteinQndx[j0]!=grdpntQst[count])
						{
							count++;
							if(count>m_ContWorld->NIndexing->QNum-1)
							{
								pnpError("Grid position of charges from loaded Gprotein is not coinside with protein charges in protein-ion system, check grid sizes\n");
								return EXIT_FAILURE;
							}
						}
						qphi=0.0;
						iGprotein=-1;
						qphi=GproteinQst[j0]*GproteinPhi[j0];
						Gsubtract[count]=-qphi;
					}
				}
				m_PoissonSolver=SingleWorldCalcController::CmdPoissonSolver(EltPoisson,m_ContWorld);
				DeleteCArray(grdpntQst);
				DeleteCArray(Gsubtract);
				
				PMFCalc[i][j]=true;
				PMF[i][j]=m_ContWorld->SystemEnergy;
				
				PrintSIPResultSingle(i,j,PMF[i][j],m_PoissonSolver->totalChange,m_PoissonSolver->ConvFac);
				
				DeleteObjByPnt(m_PoissonSolver);
				StopClockWMes0("one point calculation");
				
				if(WriteMapsForEachPosition)
				{
					pnpsapp->AddProcessNumberToFileName(filenameMapsOut,"SIP_P.bin","_i",i,NIons);
					pnpsapp->AddProcessNumberToFileName(filenameMapsOut,filenameMapsOut,"_p",j,Npnts[i]);
					m_ContWorld->WritePotential(filenameMapsOut,1);
				}
				goIon->ChangeUnits(m_ContWorld,false);
			}
		}
	}
	m_BuildWorld->GOElms.pop_back();
	m_BuildWorldRough->GOElms.pop_back();
	//Change Parameters to external
	m_BuildWorld->ChangeUnits(m_ContWorld,false);
	m_BuildWorldRough->ChangeUnits(m_ContWorldRough,false);
	for(iElt=0;iElt<m_BuildWorld->GOElms.size();iElt++)
	{
		m_BuildWorld->GOElms[iElt]->ChangeUnits(m_ContWorld,false);
		m_BuildWorldRough->GOElms[iElt]->ChangeUnits(m_ContWorldRough,false);
	}
	
	fprintf(stdout,"SimplePMFCalculation %d/%d DONE\n",pnpsapp->GetMyGroupNumber(),pnpsapp->GetNumberOfGroups());
	m_ContWorldRough->Potential=NULL;
	delete PotRough;
	m_ContWorld->Potential=NULL;
	delete Pot;
	
	delete goIon;
	delete m_ContWorldRough;
	delete m_BuildWorld;
	delete m_BuildWorldRough;
	return EXIT_SUCCESS;
}
int PMFCalculation::RunOld(const TiXmlElement *RootElt)
{
	/*DefClock0;
	int i,j,k;
	const TiXmlElement *CldElt;
	const TiXmlElement *PoissonElt;
	const TiXmlElement *BuildWorldElt;
	
	if(RootElt==NULL)return EXIT_FAILURE;
	
	//Load All Parameters
	CldElt=RootElt->FirstChildElement("IonsForPMF");
	PNP_EXIT_FAIL_NULL(CldElt,"ERROR: There is no <IonsForPMF> in <PMFCalculation>\n");
	LoadIonsForPMF(CldElt);
	
	CldElt=RootElt->FirstChildElement("ContWorld");
	PNP_EXIT_FAIL_NULL(CldElt,"ERROR: There is no <ContWorld> in <PMFCalculation>\n");
	m_ContWorld=SingleWorldCalcController::CmdContWorld(CldElt);
	
	BuildWorldElt=RootElt->FirstChildElement("BuildWorldNI");
	PNP_EXIT_FAIL_NULL(BuildWorldElt,"ERROR: There is no <BuildWorldNI> in <PMFCalculation>\n");
	m_BuildWorldNI=CmdBuildWorld(BuildWorldElt,m_ContWorld);
	
	PoissonElt=RootElt->FirstChildElement("PoissonSolver");
	PNP_EXIT_FAIL_NULL(CldElt,"ERROR: There is no <PoissonSolver> in <PMFCalculation>\n");
	
	CldElt=RootElt->FirstChildElement("PointsForPMF");
	PNP_EXIT_FAIL_NULL(CldElt,"ERROR: There is no <PointsForPMF> in <PMFCalculation>\n");
	
	if(pnpsapp->GetMyAbsRank()==pnpsapp->GetMaster())
		LoadPointsForPMF(CldElt);
	
	//Some more Parameters
	int Type=0;
	if(RootElt->Attribute("Type")!=NULL)
	{
		if(strcmp("Simple",RootElt->Attribute("Type"))==0)
			Type=0;
		else if(strcmp("Caching",RootElt->Attribute("Type"))==0)
			Type=1;
		else if(strcmp("SimpleFocusing",RootElt->Attribute("Type"))==0)
			Type=2;
	}
	
	DbgPrint0("LoadGprotein if needed\n");
	if(RootElt->Attribute("LoadGprotein")!=NULL)
	{
		LoadGprotein(RootElt->Attribute("LoadGprotein"));
	}
	
	
	DbgPrint0("Check if something was already calculated\n");
	if(RootElt->Attribute("TempResultsOut")!=NULL)
	{
		sprintf(TempResultsOut,"%s\0",RootElt->Attribute("TempResultsOut"));
		//Check if something was already calculated
		if(pnpsapp->GetMyAbsRank()==pnpsapp->GetMaster())
		{
			bool OnePntPerNode=false;
			int count;
			
			RootElt->GetBoolAttribute("OnePntPerNode",&OnePntPerNode);
			if(OnePntPerNode)
				count=LoadCalculatedPointsOnePntPerNode(TempResultsOut);
			else
				count=0;//LoadCalculatedPoints(TempResultsOut);
			
			if(count>0)
				RearrangePoints();
		}
	}
	else
	{
		sprintf(TempResultsOut,"TempResultsOut.dat\0");
	}
	DbgPrint0("Distrebute points over processes\n");
	ReDistrebutePntsOverProcesses();
#ifdef MPI_PARALLEL
	char Temp[PNP_MAP_IO_STRING_LENGTH];
	sprintf(Temp,"%s\0",TempResultsOut);
	pnpsapp->AddMyGroupNumberToFileName(TempResultsOut,Temp);
#endif
	//do Calculation
#ifdef MPI_PARALLEL
	
	switch(Type)
	{
		case 0:
			SimplePMFCalculationMPI0(RootElt,PoissonElt);
			break;
		case 1:
			//CachingPMFCalculationMPI1(RootElt,PoissonElt);
			break;
	};
#else
	switch(Type)
	{
		case 0:
			SimplePMFCalculation(PoissonElt);
			break;
		case 1:
			//CachingPMFCalculation(PoissonElt);
			break;
	};
#endif
	pnpPrint0("PMFCalculation::Run DONE\n");*/
	return EXIT_SUCCESS;
}
int PMFCalculation::ReDistrebutePntsOverProcesses()
{
	int i,j;
#ifdef MPI_PARALLEL
	if(pnpsapp->GetMyAbsRank()==pnpsapp->GetMaster())
	{
		for(i=0;i<NIons;i++)
		{
			pnpsapp->MyComSuperGroup.Bcast(&(Npnts[i]), 1, MPI::INT, pnpsapp->GetMaster());
			
			pnpsapp->MyComSuperGroup.Bcast(PMFPnts[i], Npnts[i], MPI::INT, pnpsapp->GetMaster());
			for(j=0;j<3;j++)
				pnpsapp->MyComSuperGroup.Bcast(PMFxyz[j][i], Npnts[i], MPI::FLOAT, pnpsapp->GetMaster());
			pnpsapp->MyComSuperGroup.Bcast(PMF[i], Npnts[i], MPI::DOUBLE, pnpsapp->GetMaster());
			
			pnpsapp->MyComSuperGroup.Bcast(PMFPntsStart,NIons,MPI::INT, pnpsapp->GetMaster());
			
			int *PMFCalcBOOL=new int[Npnts[i]];
			for(j=0;j<Npnts[i];j++)
				PMFCalcBOOL[j]=PMFCalc[i][j];
			pnpsapp->MyComSuperGroup.Bcast(PMFCalcBOOL, Npnts[i], MPI::INT, pnpsapp->GetMaster());
			delete [] PMFCalcBOOL;
		}
	}
	else
	{
		for(i=0;i<NIons;i++)
		{
			pnpsapp->MyComSuperGroup.Bcast(&(Npnts[i]), 1, MPI::INT, pnpsapp->GetMaster());
			
			PMFPnts[i]=new int[Npnts[i]];
			for(j=0;j<3;j++)
				PMFxyz[j][i]=new float[Npnts[i]];
			PMF[i]=new double[Npnts[i]];
			PMFCalc[i]=new bool[Npnts[i]];
			
			pnpsapp->MyComSuperGroup.Bcast(PMFPnts[i], Npnts[i], MPI::INT, pnpsapp->GetMaster());
			for(j=0;j<3;j++)
				pnpsapp->MyComSuperGroup.Bcast(PMFxyz[j][i], Npnts[i], MPI::FLOAT, pnpsapp->GetMaster());
			pnpsapp->MyComSuperGroup.Bcast(PMF[i], Npnts[i], MPI::DOUBLE, pnpsapp->GetMaster());
			
			pnpsapp->MyComSuperGroup.Bcast(PMFPntsStart,NIons,MPI::INT, pnpsapp->GetMaster());
			
			int *PMFCalcBOOL=new int[Npnts[i]];
			pnpsapp->MyComSuperGroup.Bcast(PMFCalcBOOL, Npnts[i], MPI::INT, pnpsapp->GetMaster());
			for(j=0;j<Npnts[i];j++)
				PMFCalc[i][j]=PMFCalcBOOL[j];
			delete [] PMFCalcBOOL;
		}
	}
#endif
	return EXIT_SUCCESS;
}
int PMFCalculation::GetPntPos(int ion,int node)
{
	int iPos=-1;
	int i;
	
	for(i=0;i<Npnts[ion];i++)
	{
		if(PMFPnts[ion][i]==node)return i;
	}
	
	return iPos;
}
int PMFCalculation::LoadCalculatedPointsOnePntPerNode(const char *FileName)
{
	DbgPrint0("PMFCalculation::LoadCalculatedPointsOnePntPerNode\n");
	int count=0;
	
	pnpPrint0("Will read calculated points from %s\n",FileName);
	FILE *out=fopen(FileName,"rt");
	
	int irc[3];
	float frc[3];
	float rabs[3];
	int ion;
	int iPos;
	int ClstNode;
	double E;
	double dE;
	double Conv;
	char str[PNP_MAP_IO_STRING_LENGTH];
	
	if(out==NULL)
	{
		pnpWarning("Can't read file:%s\n",FileName);
		return 0;
	}
	
	while(fgets(str,PNP_MAP_IO_STRING_LENGTH,out)!=NULL)
	{
		if(sscanf(str,"%d %d %d %g %g %g %le %le %le", &ion, &iPos, &ClstNode, &rabs[0], &rabs[1], &rabs[2], &E, &dE, &Conv)==9)
		{
			iPos=GetPntPos(ion, ClstNode);
			if(iPos>=0)
			{
				if(PMFCalc[ion][iPos]==false)
				{
					PMF[ion][iPos]=E;
					PMFCalc[ion][iPos]=true;
					count++;
				}
			}
		}
	}
	fclose(out);
	
	pnpPrint0("%d positions was already calculated\n",count);
	return count;
}
int PMFCalculation::RearrangePoints()
{
	int ion;
	int count;
	int j,k;
	//rearrange poinst for calculation, calculated first, needed to calculate last
	for(ion=0;ion<NIons;ion++)
	{
		count=0;
		for(j=0;j<Npnts[ion];j++)
		{
			if(PMFCalc[ion][j]==true)
				count++;
		}
		if(count==Npnts[ion])
		{
			pnpPrint0("Position of ion %d was complitly calculated\n",ion);
			PMFPntsStart[ion]=Npnts[ion];
		}
		else
		{
			pnpPrint0("For ion %d %d position of %d was calculated\n",ion,count,Npnts[ion]);
			PMFPntsStart[ion]=count;
			
			for(j=0;j<Npnts[ion];j++)
			{
				if(PMFCalc[ion][j]==false)
				{
					k=j+1;
					while(PMFCalc[ion][k]==false&&k<Npnts[ion])
						k++;
					if(k>=Npnts[ion])
						break;
					//move k to j and j to k
					//pnpPrint0("will move %d to %d %d %d\n",k,
					bool btmp;
					float ftmp;
					double dtmp;
					int itmp;
					//flag index
					btmp=PMFCalc[ion][j];
					PMFCalc[ion][j]=PMFCalc[ion][k];
					PMFCalc[ion][k]=btmp;
					//PMF
					dtmp=PMF[ion][j];
					PMF[ion][j]=PMF[ion][k];
					PMF[ion][k]=dtmp;
					//grid index
					itmp=PMFPnts[ion][j];
					PMFPnts[ion][j]=PMFPnts[ion][k];
					PMFPnts[ion][k]=itmp;
					//xyz coor
					ftmp=PMFxyz[0][ion][j];
					PMFxyz[0][ion][j]=PMFxyz[0][ion][k];
					PMFxyz[0][ion][k]=ftmp;
					
					ftmp=PMFxyz[1][ion][j];
					PMFxyz[1][ion][j]=PMFxyz[1][ion][k];
					PMFxyz[1][ion][k]=ftmp;
					
					ftmp=PMFxyz[2][ion][j];
					PMFxyz[2][ion][j]=PMFxyz[2][ion][k];
					PMFxyz[2][ion][k]=ftmp;
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
int PMFCalculation::LoadCalculatedPoints(const char *FileName)
{
	DbgPrint2("PMFCalculation::LoadCalculatedPoints\n");
	int i,j,k;
	int ion,pnt,PosNum;
	double dPMF;
	int count=0;
	char str[100];
	pnpPrint0("Will read calculated points from %s\n",FileName);
	FILE *out=fopen(FileName,"rt");
	
	if(out==NULL)
	{
		pnpError("Can't read file:%s\n",FileName);
		return 0;
	}
	while(fgets(str,100,out)!=NULL)
	{
		if(sscanf(str,"%d %d %le",&ion,&pnt,&dPMF)==3)
		{
			PosNum=-1;
			for(i=0;i<Npnts[ion];i++)
			{
				if(PMFPnts[ion][i]==pnt)
				{
					PosNum=i;
					break;
				}
			}
			if(PosNum>=0)
			{
				PMF[ion][PosNum]=dPMF;
				PMFCalc[ion][PosNum]=true;
				count++;
			}
		}
	}
	fclose(out);
		
	pnpPrint0("%d positions was already calculated\n",count);
	//rearrange poinst for calculation, calculated first, needed to calculate last
	for(ion=0;ion<NIons;ion++)
	{
		count=0;
		for(j=0;j<Npnts[ion];j++)
		{
			if(PMFCalc[ion][j]==true)
				count++;
		}
		if(count==Npnts[ion])
		{
			pnpPrint0("Position of ion %d was complitly calculated\n",ion);
			PMFPntsStart[ion]=Npnts[ion];
		}
		else
		{
			pnpPrint0("For ion %d %d position of %d was calculated\n",ion,count,Npnts[ion]);
			PMFPntsStart[ion]=count;
			
			for(j=0;j<Npnts[ion];j++)
			{
				if(PMFCalc[ion][j]==false)
				{
					k=j+1;
					while(PMFCalc[ion][k]==false&&k<Npnts[ion])
						k++;
					if(k>=Npnts[ion])
						break;
					//move k to j and j to k
					//pnpPrint0("will move %d to %d %d %d\n",k,
					bool btmp;
					float ftmp;
					double dtmp;
					int itmp;
					//flag index
					btmp=PMFCalc[ion][j];
					PMFCalc[ion][j]=PMFCalc[ion][k];
					PMFCalc[ion][k]=btmp;
					//PMF
					dtmp=PMF[ion][j];
					PMF[ion][j]=PMF[ion][k];
					PMF[ion][k]=dtmp;
					//grid index
					itmp=PMFPnts[ion][j];
					PMFPnts[ion][j]=PMFPnts[ion][k];
					PMFPnts[ion][k]=itmp;
					//xyz coor
					ftmp=PMFxyz[0][ion][j];
					PMFxyz[0][ion][j]=PMFxyz[0][ion][k];
					PMFxyz[0][ion][k]=ftmp;
					
					ftmp=PMFxyz[1][ion][j];
					PMFxyz[1][ion][j]=PMFxyz[1][ion][k];
					PMFxyz[1][ion][k]=ftmp;
					
					ftmp=PMFxyz[2][ion][j];
					PMFxyz[2][ion][j]=PMFxyz[2][ion][k];
					PMFxyz[2][ion][k]=ftmp;
				}
			}
		}
	}
	
	return EXIT_SUCCESS;
}
int PMFCalculation::SimplePMFCalculation(const TiXmlElement *PoissonElt)
{
	/*int i,j,k;
	
	fprintf(stdout,"SimplePMFCalculation\n");
	//Change Parameters for internal 
	int iElt;
	m_BuildWorldNI->ChangeUnits(m_ContWorld,true);
	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
	{
		m_BuildWorldNI->GOElms[iElt]->ChangeUnits(m_ContWorld,true);
	}
	
	GOAtoms* goIon;
	m_BuildWorldNI->GOElms.resize(m_BuildWorldNI->GOElms.size()+1);
	for(i=0;i<NIons;i++)
	{
		
		goIon=IonsForPMF[i];
		goIon->ChangeUnits(m_ContWorld,true);
		m_BuildWorldNI->GOElms[m_BuildWorldNI->GOElms.size()-1]=goIon;
		for(j=0;j<Npnts[i];j++)
		{
			if(PMFCalc[i][j]==false)
			{
				DefClock0;
				StartClock0;
				SetIonCoor(i,j);
				fprintf(stdout,"Start PMPPnts[%d][%d]=%d, r=[%f,%f,%f][int.coor.]\n",i,j,PMFPnts[i][j],goIon->r[0][0],goIon->r[1][0],goIon->r[2][0]);
				BuildContWorld(m_ContWorld,m_BuildWorldNI);
				
				//init grdpntQst and Gsubtract
				DeleteCArray(grdpntQst);
				DeleteCArray(Gsubtract);
				grdpntQst=m_ContWorld->NIndexing->GetChargeIndex();
				if(GproteinLoad||GionLoad)
				{
					if(Gsubtract==NULL)
					{
						Gsubtract=new double[m_ContWorld->NIndexing->QNum];
						for(k=0;k<m_ContWorld->NIndexing->QNum;k++)
							Gsubtract[k]=0.0;
					}
		
					if(GproteinLoad)
					{
						if(GproteinQnum>m_ContWorld->NIndexing->QNum-1)
						{
							pnpWarning("Number of partial charges in Gprotein is bigger then partial charges in protein-ion - 1\n");
						}
						double qphi;
						int iGprotein;
						int count=0;
						int j0;
						for(j0=0;j0<GproteinQnum;j0++)
						{
							while(GproteinQndx[j0]!=grdpntQst[count])
							{
								count++;
							}
							qphi=0.0;
							iGprotein=-1;
							qphi=GproteinQst[j0]*GproteinPhi[j0];
							Gsubtract[count]=-qphi;
						}
					}
				}
				
				m_PoissonSolver=SingleWorldCalcController::CmdPoissonSolver(PoissonElt,m_ContWorld);
				
				PMFCalc[i][j]=true;
				PMF[i][j]=m_ContWorld->SystemEnergy;
				
				PrintSIPResultSingle(i,j,PMF[i][j],m_PoissonSolver->totalChange,m_PoissonSolver->ConvFac);
				
				DeleteObjByPnt(m_PoissonSolver);
				StopClockWMes0("one point calculation");
			}
		}
		goIon->ChangeUnits(m_ContWorld,false);
	}
	m_BuildWorldNI->GOElms.pop_back();
	//Change Parameters to external
	m_BuildWorldNI->ChangeUnits(m_ContWorld,false);
	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
	{
		m_BuildWorldNI->GOElms[iElt]->ChangeUnits(m_ContWorld,false);
	}*/
	return EXIT_SUCCESS;
}
int PMFCalculation::SimplePMFCalculationMPI0(const TiXmlElement *PMFElt,const TiXmlElement *PoissonElt)
{
#ifdef MPI_PARALLEL
	/*int i,j,k;	
	
	//fprintf(stdout,"SimplePMFCalculation %d/%d %s\n",pnpsapp->MyRank,pnpsapp->NProcs,TempResultsOut);
	//Change Parameters for internal 
	int iElt;
	m_BuildWorldNI->ChangeUnits(m_ContWorld,true);
	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
	{
		m_BuildWorldNI->GOElms[iElt]->ChangeUnits(m_ContWorld,true);
	}
	
	GOAtoms* goIon;
	m_BuildWorldNI->GOElms.resize(m_BuildWorldNI->GOElms.size()+1);
	
	
	int start,end,npnts,pnts4proc;
	for(i=0;i<NIons;i++)
	{
		
		goIon=IonsForPMF[i];
		goIon->ChangeUnits(m_ContWorld,true);
		m_BuildWorldNI->GOElms[m_BuildWorldNI->GOElms.size()-1]=goIon;
		
		npnts=Npnts[i]-PMFPntsStart[i];
		start=PMFPntsStart[i]+npnts*pnpsapp->GetMyGroupNumber()/pnpsapp->GetNumberOfGroups();
		end=PMFPntsStart[i]+npnts*(pnpsapp->GetMyGroupNumber()+1)/pnpsapp->GetNumberOfGroups();
		pnts4proc=end-start;
		if(pnpsapp->GetMyGroupNumber()==pnpsapp->GetNumberOfGroups()-1)end=PMFPntsStart[i]+npnts;
		DbgPrint0("SimplePMFCalculation: start=%d end=%d npnts=%d pnts4proc=%d\n", start, end, npnts, pnts4proc);
		
		DefClock0;
		
		for(j=start;j<end;j++)
		{
			pnpPrint0("PMFCalc[%d][0] %d %d\n",i,PMFCalc[i][j],PMFCalc[i][start]);
			if(PMFCalc[i][j]==false)
			{
				SetIonCoor(i,j);
				fprintf(stdout,"Start PMPPnts[%d][%d]=%d, r=[%f,%f,%f][int.coor.]\n",i,j,PMFPnts[i][j],goIon->r[0][0],goIon->r[1][0],goIon->r[2][0]);
				
				StartClock0;
				
				BuildContWorld(m_ContWorld,m_BuildWorldNI);
				
				//init grdpntQst and Gsubtract
				DeleteCArray(grdpntQst);
				DeleteCArray(Gsubtract);
				grdpntQst=m_ContWorld->NIndexing->GetChargeIndex();
				if(GproteinLoad||GionLoad)
				{
					if(Gsubtract==NULL)
					{
						Gsubtract=new double[m_ContWorld->NIndexing->QNum];
						for(k=0;k<m_ContWorld->NIndexing->QNum;k++)
							Gsubtract[k]=0.0;
					}
		
					if(GproteinLoad)
					{
						if(GproteinQnum>m_ContWorld->NIndexing->QNum-1)
						{
							pnpWarning("Number of partial charges in Gprotein is bigger then partial charges in protein-ion - 1\n");
						}
						double qphi;
						int iGprotein;
						int count=0;
						int j0;
						for(j0=0;j0<GproteinQnum;j0++)
						{
							while(GproteinQndx[j0]!=grdpntQst[count])
							{
								count++;
							}
							qphi=0.0;
							iGprotein=-1;
							qphi=GproteinQst[j0]*GproteinPhi[j0];
							Gsubtract[count]=-qphi;
						}
					}
				}
				m_PoissonSolver=SingleWorldCalcController::CmdPoissonSolver(PoissonElt,m_ContWorld);
				
				PMFCalc[i][j]=true;
				PMF[i][j]=m_ContWorld->SystemEnergy;
				
				PrintSIPResultSingle(i,j,PMF[i][j],m_PoissonSolver->totalChange,m_PoissonSolver->ConvFac);
				
				DeleteObjByPnt(m_PoissonSolver);
				StopClockWMes0("one point calculation");
			}
		}
		goIon->ChangeUnits(m_ContWorld,false);
	}
	m_BuildWorldNI->GOElms.pop_back();
	
	//Change Parameters to external
	m_BuildWorldNI->ChangeUnits(m_ContWorld,false);
	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
	{
		m_BuildWorldNI->GOElms[iElt]->ChangeUnits(m_ContWorld,false);
	}
	fprintf(stdout,"SimplePMFCalculation %d/%d DONE\n",pnpsapp->GetMyGroupNumber(),pnpsapp->GetNumberOfGroups());*/
#endif
	return EXIT_SUCCESS;
}
// int PMFCalculation::CachingPMFCalculation(const TiXmlElement *PoissonElt)
// {
// 	int i,j;
// 	int iElt;
// 	
// 	//Change Parameters for internal 
// 	m_BuildWorldNI->ChangeUnits(m_ContWorld,true);
// 	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
// 	{
// 		m_BuildWorldNI->GOElms[iElt]->ChangeUnits(m_ContWorld,true);
// 	}
// 	
// 	fprintf(stdout,"CachingPMFCalculation\n");
// 	GOAtoms* goIon;
// 	PreBuildContWorld(m_ContWorld,m_BuildWorldNI);
// 	for(i=0;i<NIons;i++)
// 	{
// 		pnpPrint0("CachingPMFCalculation for ion %d\n",i);
// 		goIon=IonsForPMF[i];
// 		goIon->ChangeUnits(m_ContWorld,true);
// 		for(j=0;j<Npnts[i];j++)
// 		{
// 			if(PMFCalc[i][j]==false)
// 			{
// 				SetIonCoor(i,j);
// 				
// 				
// 				DbgPrint0("Start PMPPnts[%d][%d]=%d, r=[%f,%f,%f][int.coor.] q = %f [int.charge]\n"
// 						,i,j,PMFPnts[i][j],goIon->r[0][0],goIon->r[1][0],goIon->r[2][0],goIon->q[0]);
// 				FinalizePreBuildContWorld(m_ContWorld,m_BuildWorldNI,goIon);
// 				m_PoissonSolver=SingleWorldCalcController::CmdPoissonSolver(PoissonElt,m_ContWorld);
// 				
// 				PMFCalc[i][j]=true;
// 				PMF[i][j]=m_ContWorld->SystemEnergy;
// 				
// 				PrintSIPResultSingle(i,j,PMF[i][j],m_PoissonSolver->totalChange,m_PoissonSolver->ConvFac);
// 				
// 				DeleteObjByPnt(m_PoissonSolver);
// 			}
// 		}
// 		goIon->ChangeUnits(m_ContWorld,false);
// 	}
// 	
// 	//Change Parameters to external
// 	m_BuildWorldNI->ChangeUnits(m_ContWorld,false);
// 	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
// 	{
// 		m_BuildWorldNI->GOElms[iElt]->ChangeUnits(m_ContWorld,false);
// 	}
// 	
// 	return EXIT_SUCCESS;
// }
void PMFCalculation::PrintSIPResultSingle(int ion, int iPos,double E,double dE,double Conv)
{
	float rabs[3]={PMFxyz[0][ion][iPos],PMFxyz[1][ion][iPos],PMFxyz[2][ion][iPos]};
	int ClstNode=PMFPnts[ion][iPos];

	FILE* PMFTMPOUT=fopen(TempResultsOut,"at");
	
	fprintf(PMFTMPOUT,"%d\t%d\t%d\t%.5g\t%.5g\t%.5g\t%.16e\t%.4e\t%.4e\n"
			,ion,iPos,ClstNode,rabs[0],rabs[1],rabs[2],E,dE,Conv);
	fclose(PMFTMPOUT);
	
	pnpPrint0("<pmfresults iIon=\"%d\" iPosition=\"%d\" ClstNode=\"%d\" r=\"%.6g, %.6g, %.6g\"  PMF=\" %.16e\" dE=\"%.4e\" ConvFac=\"%.4e\"/>\n"
			,ion,iPos,ClstNode,rabs[0],rabs[1],rabs[2],E,dE,Conv);
}
void PMFCalculation::PrintSIPResultSingleRough(int ion, int iPos,double E,double dE,double Conv)
{
	float rabs[3]={PMFxyz[0][ion][iPos],PMFxyz[1][ion][iPos],PMFxyz[2][ion][iPos]};
	int ClstNode=m_ContWorldRough->ConvGlobalExternalXYZToGrdPnt(PMFxyz[0][ion][iPos], PMFxyz[1][ion][iPos], PMFxyz[2][ion][iPos]);
	

	FILE* PMFTMPOUT=fopen(TempResultsOutRough,"at");
	
	fprintf(PMFTMPOUT,"%d\t%d\t%d\t%.5g\t%.5g\t%.5g\t%.16e\t%.4e\t%.4e\n"
			,ion,iPos,ClstNode,rabs[0],rabs[1],rabs[2],E,dE,Conv);
	fclose(PMFTMPOUT);
	
	pnpPrint0("<pmfresults iIon=\"%d\" iPosition=\"%d\" ClstNode=\"%d\" r=\"%.6g, %.6g, %.6g\"  PMF=\" %.16e\" dE=\"%.4e\" ConvFac=\"%.4e\"/>\n"
			,ion,iPos,ClstNode,rabs[0],rabs[1],rabs[2],E,dE,Conv);
}
// int PMFCalculation::CachingPMFCalculationMPI0(const TiXmlElement *PMFElt,const TiXmlElement *PoissonElt)
// {
// #ifdef MPI_PARALLEL
// 	int i,j,np;
// 	
// 	int iElt;
// 	//Change Parameters for internal 
// 	m_BuildWorldNI->ChangeUnits(m_ContWorld,true);
// 	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
// 	{
// 		m_BuildWorldNI->GOElms[iElt]->ChangeUnits(m_ContWorld,true);
// 	}
// 	
// 	DbgPrint0("CachingPMFCalculation\n");
// 	GOAtoms* goIon;
// 	PreBuildContWorld(m_ContWorld,m_BuildWorldNI);
// 
// 	MPI::COMM_WORLD.Barrier();
// 	for(np=0;np<pnpsapp->NProcs;np++)
// 	{
// 		if(pnpsapp->MyRank==np)
// 			DbgPrint0("CachingPMFCalculationMPI0:Ready to start PMF Calculation\n");
// 		MPI::COMM_WORLD.Barrier();
// 	}
// 	int start,end,npnts,pnts4proc;
// 	for(i=0;i<NIons;i++)
// 	{
// 		goIon=IonsForPMF[i];
// 		goIon->ChangeUnits(m_ContWorld,true);
// 		
// 		npnts=Npnts[i];
// 		pnts4proc=npnts/pnpsapp->NProcs;
// 		start=pnts4proc*pnpsapp->MyRank;
// 		end=pnts4proc*(pnpsapp->MyRank+1);
// 		if(pnpsapp->MyRank==pnpsapp->NProcs-1)end=npnts;
// 		DbgPrint0("start=%d end=%d npnts=%d pnts4proc=%d\n",start,end,npnts,pnts4proc);
// 		
// 		for(j=start;j<end;j++)
// 		{
// 			if(PMFCalc[i][j]==false)
// 			{
// 				SetIonCoor(i,j);
// 				fprintf(stdout,"Start PMPPnts[%d][%d]=%d, r=[%f,%f,%f][int.coor.]\n"
// 						,i,j,PMFPnts[i][j],goIon->r[0][0],goIon->r[1][0],goIon->r[2][0]);
// 				FinalizePreBuildContWorld(m_ContWorld,m_BuildWorldNI,goIon);
// 				m_PoissonSolver=PNPSapp::CmdPoissonSolver(PoissonElt,m_ContWorld);
// 				
// 				PMFCalc[i][j]=true;
// 				PMF[i][j]=m_ContWorld->SystemEnergy;
// 				
// 				PrintSIPResultSingle(i,j,PMF[i][j],m_PoissonSolver->totalChange,m_PoissonSolver->ConvFac);
// 				
// 				DeleteObjByPnt(m_PoissonSolver);
// 			}
// 		}
// 		goIon->ChangeUnits(m_ContWorld,false);
// 	}
// 	m_BuildWorldNI->ChangeUnits(m_ContWorld,false);
// 	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
// 	{
// 		m_BuildWorldNI->GOElms[iElt]->ChangeUnits(m_ContWorld,false);
// 	}
// 	MPI::COMM_WORLD.Barrier();
// 	for(np=0;np<pnpsapp->NProcs;np++)
// 	{
// 		if(pnpsapp->MyRank==np)
// 			DbgPrint0("DONE WITH PMF CALCULATION\n");
// 		MPI::COMM_WORLD.Barrier();
// 	}
// #endif
// 	return EXIT_SUCCESS;
// }
// int PMFCalculation::CachingPMFCalculationMPI1(const TiXmlElement *PMFElt,const TiXmlElement *PoissonElt)
// {
// #ifdef MPI_PARALLEL
// 	int i,j,np;
// 	
// 	int iElt;
// 	//Change Parameters for internal 
// 	m_BuildWorldNI->ChangeUnits(m_ContWorld,true);
// 	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
// 	{
// 		m_BuildWorldNI->GOElms[iElt]->ChangeUnits(m_ContWorld,true);
// 	}
// 	
// 	DbgPrint0("CachingPMFCalculation\n");
// 	GOAtoms* goIon;
// 	
// 	PreBuildContWorld1(m_ContWorld,m_BuildWorldNI);
// 
// 	MPI::COMM_WORLD.Barrier();
// 	for(np=0;np<pnpsapp->NProcs;np++)
// 	{
// 		if(pnpsapp->MyRank==np)
// 			DbgPrint0("CachingPMFCalculationMPI0:Ready to start PMF Calculation\n");
// 		MPI::COMM_WORLD.Barrier();
// 	}
// 	
// 	int start,end,npnts,pnts4proc;
// 	for(i=0;i<NIons;i++)
// 	{
// 		goIon=IonsForPMF[i];
// 		goIon->ChangeUnits(m_ContWorld,true);
// 		
// 		npnts=Npnts[i];
// 		pnts4proc=npnts/pnpsapp->NProcs;
// 		start=pnts4proc*pnpsapp->MyRank;
// 		end=pnts4proc*(pnpsapp->MyRank+1);
// 		if(pnpsapp->MyRank==pnpsapp->NProcs-1)end=npnts;
// 		DbgPrint0("start=%d end=%d npnts=%d pnts4proc=%d\n",start,end,npnts,pnts4proc);
// 		
// 		DefClock0;
// 		
// 		for(j=start;j<end;j++)
// 		{
// 			if(PMFCalc[i][j]==false)
// 			{
// 				StartClock0;
// 				SetIonCoor(i,j);
// 				fprintf(stdout,"Start PMPPnts[%d][%d]=%d, r=[%f,%f,%f][int.coor.]\n"
// 						,i,j,PMFPnts[i][j],goIon->r[0][0],goIon->r[1][0],goIon->r[2][0]);
// 				FinalizePreBuildContWorld1(m_ContWorld,m_BuildWorldNI,goIon);
// 				StopClockWMes0("BuildContWorld");
// 				StartClock0;
// 				
// 				m_PoissonSolver=new PoissonSolver();
// 				m_PoissonSolver->LoadXML(PoissonElt);
// 				m_PoissonSolver->SetContWorld(m_ContWorld);
// 				m_PoissonSolver->InitSolver();
// 				if(m_PoissonSolver->Solve()==EXIT_FAILURE)
// 				{
// 					int ig;
// 					for(ig=0;ig<m_ContWorld->GS_XYZ;ig++)
// 					{
// 						m_ContWorld->Potential[ig]=0.0;
// 					}
// 					if(m_PoissonSolver->Solve()==EXIT_FAILURE)
// 					{
// 						int ig;
// 						for(ig=0;ig<m_ContWorld->GS_XYZ;ig++)
// 						{
// 							m_ContWorld->Potential[ig]=0.0;
// 						}
// 					
// 					}
// 				}
// 				
// 				StopClockWMes0("PoissonSolver");
// 				
// 				PMFCalc[i][j]=true;
// 				PMF[i][j]=m_ContWorld->SystemEnergy;
// 				
// 				PrintSIPResultSingle(i,j,PMF[i][j],m_PoissonSolver->totalChange,m_PoissonSolver->ConvFac);
// 				
// 				DeleteObjByPnt(m_PoissonSolver);
// 			}
// 		}
// 		goIon->ChangeUnits(m_ContWorld,false);
// 	}
// 	m_BuildWorldNI->ChangeUnits(m_ContWorld,false);
// 	for(iElt=0;iElt<m_BuildWorldNI->GOElms.size();iElt++)
// 	{
// 		m_BuildWorldNI->GOElms[iElt]->ChangeUnits(m_ContWorld,false);
// 	}
// 	MPI::COMM_WORLD.Barrier();
// 	for(np=0;np<pnpsapp->NProcs;np++)
// 	{
// 		if(pnpsapp->MyRank==np)
// 			DbgPrint0("DONE WITH PMF CALCULATION\n");
// 		MPI::COMM_WORLD.Barrier();
// 	}
// #endif
// 	return EXIT_SUCCESS;
// }
// int PMFCalculation::CachingPMFCalculationMPI(const TiXmlElement *PMFElt,const TiXmlElement *PoissonElt)
// {
// #ifdef MPI_PARALLEL
// // 	MPI::Status	status;
// // 	int i,j;
// // 	int ion,PosNum,SlaveProc;
// // 	int ion2,PosNum2;
// // 	int ToDo;
// // 	double mPMF;
// // 	bool HaveSomethingToDo=true;
// // 	fprintf(stdout,"CachingPMFCalculation\n");
// // 	DbgPrint0("CachingPMFCalculation\n");
// // 	if(pnpsapp->MyRank==0)
// // 	{
// // 		ion=0;
// // 		PosNum=0;
// // 		while(PMFCalc[ion][PosNum]!=false)
// // 		{
// // 			PosNum++;
// // 			if(PosNum>=Npnts[ion])
// // 			{
// // 				ion++;
// // 				PosNum=0;
// // 				if(ion>=NIons)
// // 				{
// // 					HaveSomethingToDo=false;
// // 					break;
// // 				}
// // 			}
// // 		}
// // 		
// // 		//Cycle for tasks sending
// // 		while(HaveSomethingToDo)
// // 		{
// // 			if(pnpsapp->NProcs==1)
// // 			{
// // 				fprintf(stdout,"ERROR: Only 1 proccess. Noprocesses for calculation\n");
// // 				HaveSomethingToDo=false;
// // 				break;
// // 			}
// // 			//What Process and what to do
// // 			DbgPrint0("Ready To Start sending points\n");
// // 			MPI::COMM_WORLD.Recv(&SlaveProc, 1, MPI::INT, MPI::ANY_SOURCE, MES_WHO);
// // 			DbgPrint0("Ready To Start sending points: Proc:%d\n",SlaveProc);
// // 			ToDo=MES_HAVETOCALC;
// // 			MPI::COMM_WORLD.Send(&ToDo, 1, MPI::INT, SlaveProc, MES_WHATTODO);
// // 			
// // 			MPI::COMM_WORLD.Recv(&ToDo, 1, MPI::INT, SlaveProc, MES_WHATTODO);
// // 			DbgPrint0("Proc:%d want to do %d\n",SlaveProc,ToDo);
// // 			
// // 			if(ToDo==MES_WANT_CALC)
// // 			{
// // 				DbgPrint0("Looking for Job to Calculate for Proc=%d\n",SlaveProc,ToDo);
// // 				while(PMFCalc[ion][PosNum]!=false)
// // 				{
// // 					PosNum++;
// // 					if(PosNum>=Npnts[ion])
// // 					{
// // 						ion++;
// // 						PosNum=0;
// // 						if(ion>=NIons)
// // 						{
// // 							HaveSomethingToDo=false;
// // 							break;
// // 						}
// // 					}
// // 				}
// // 				//FindWhatToCalculate
// // 				DbgPrint0("Will assign ion=%d point[%d]=%d for Proc=%d\n",ion,PosNum,PMFPnts[ion][PosNum],SlaveProc);
// // 				MPI::COMM_WORLD.Send(&ion, 1, MPI::INT, 0, MES_ION);
// // 				MPI::COMM_WORLD.Send(&PosNum, 1, MPI::INT, 0, MES_PNT);
// // 				DbgPrint0("Have sent task to Proc=%d\n",SlaveProc);
// // 			}
// // 			else if(ToDo==MES_HAVE_RESULTS)
// // 			{
// // 				MPI::COMM_WORLD.Recv(&ion2, 1, MPI::INT, SlaveProc, MES_ION);
// // 				MPI::COMM_WORLD.Recv(&PosNum2, 1, MPI::INT, SlaveProc, MES_PNT);
// // 				MPI::COMM_WORLD.Recv(&mPMF, 1, MPI::DOUBLE, SlaveProc, MES_PMF);
// // 				
// // 				PMFCalc[ion2][PosNum2]=true;
// // 				PMF[ion2][PosNum2]=mPMF;
// // 				
// // 				PrintSIPResultSingle(i,j,PMF[i][j],m_PoissonSolver->totalChange,m_PoissonSolver->ConvFac);
// // 				FILE* PMFTMPOUT=fopen(TempResultsOut,"at");
// // 				if(PMFTMPOUT!=NULL)
// // 				{
// // 					fprintf(PMFTMPOUT,"%3d %10d %20.16le\n"
// // 						,ion2,PMFPnts[ion2][PosNum2],PMF[ion2][PosNum2]);
// // 					fclose(PMFTMPOUT);
// // 				}
// // 				
// // 				fprintf(stdout,"<PMFTMPOUT Ion=\"%d\" Pnt=\"%d\" PMF=\" %.16le\"/>\n"
// // 						,ion2,PMFPnts[ion2][PosNum2],PMF[ion2][PosNum2]);
// // 			}
// // 		}
// // 		int HaveSentToProc=0;
// // 		while(HaveSentToProc<pnpsapp->NProcs-1)
// // 		{
// // 			MPI::COMM_WORLD.Recv(&SlaveProc, 1, MPI::INT, MPI::ANY_SOURCE, MES_WHO);
// // 			ToDo=MES_DONE;
// // 			MPI::COMM_WORLD.Send(&ToDo, 1, MPI::INT, SlaveProc, MES_WHATTODO);
// // 			HaveSentToProc++;
// // 		}
// // 		//SendToAllThatDone
// // 	}
// // 	else
// // 	{
// // 		//preproceccing
// // 		GOAtoms* goIon;
// // 		DbgPrint0("PreBuildContWorld\n");
// // 		PreBuildContWorld(m_ContWorld,m_BuildWorldNI);
// // 		for(i=0;i<NIons;i++)
// // 			IonsForPMF[i]->ChangeUnits(m_ContWorld,true);
// // 		DbgPrint0("PreBuildContWorld done\n");
// // 		while(HaveSomethingToDo)
// // 		{
// // 			//Get What To Calc
// // 			//connect to main proccess
// // 			DbgPrint0("Asking for point to calculate\n");
// // 			MPI::COMM_WORLD.Send(&(pnpsapp->MyRank), 1, MPI::INT, 0, MES_WHO);
// // 			MPI::COMM_WORLD.Recv(&ToDo, 1, MPI::INT, 0, MES_WHATTODO);
// // 			DbgPrint0("Recv status %d from main proccess\n",ToDo);
// // 			if(ToDo==MES_DONE)
// // 			{
// // 				HaveSomethingToDo=false;
// // 				break;
// // 			}
// // 			ToDo=MES_WANT_CALC;
// // 			MPI::COMM_WORLD.Send(&ToDo, 1, MPI::INT, 0, MES_WHATTODO);
// // 			//now get parameters ion PosNum
// // 			MPI::COMM_WORLD.Recv(&ion, 1, MPI::INT, 0, MES_ION);
// // 			MPI::COMM_WORLD.Recv(&PosNum, 1, MPI::INT, 0, MES_PNT);
// // 			DbgPrint0("Going to Calculate Ion=%d at point %d\n",ion,PosNum);
// // 			//Calculation
// // 			goIon=IonsForPMF[ion];
// // 			SetIonCoor(ion,PosNum);
// // 			FinalizePreBuildContWorld(m_ContWorld,m_BuildWorldNI,goIon);
// // 			m_PoissonSolver=PNPSapp::CmdPoissonSolver(PoissonElt,m_ContWorld);
// // 			DeleteObjByPnt(m_PoissonSolver);
// // 			PMFCalc[i][j]=true;
// // 			PMF[i][j]=m_ContWorld->SystemEnergy;
// // 			mPMF=m_ContWorld->SystemEnergy;
// // 			fprintf(stdout,"<PMFTMPOUT Ion=\"%d\" Pnt=\"%d\" PMF=\" %.16le\"/>\n",i,PMFPnts[i][j],PMF[i][j]);
// // 			//Send Results Back
// // 			//connect to main proccess
// // 			MPI::COMM_WORLD.Send(&(pnpsapp->MyRank), 1, MPI::INT, 0, MES_WHO);
// // 			MPI::COMM_WORLD.Recv(&ToDo, 1, MPI::INT, 0, MES_WHATTODO);
// // 			if(ToDo==MES_DONE)
// // 			{
// // 				HaveSomethingToDo=false;
// // 				break;
// // 			}
// // 			ToDo=MES_WANT_CALC;
// // 			MPI::COMM_WORLD.Send(&ToDo, 1, MPI::INT, 0, MES_WHATTODO);
// // 			//now get parameters
// // 			MPI::COMM_WORLD.Send(&ion, 1, MPI::INT, 0, MES_ION);
// // 			MPI::COMM_WORLD.Send(&PosNum, 1, MPI::INT, 0, MES_PNT);
// // 			MPI::COMM_WORLD.Send(&mPMF, 1, MPI::DOUBLE, 0, MES_PMF);
// // 		}
// // 		//postproceccing
// // 		
// // 		for(i=0;i<NIons;i++)
// // 			IonsForPMF[i]->ChangeUnits(m_ContWorld,false);
// // 	}
// #endif
// 	return EXIT_SUCCESS;
// }
int PMFCalculation::SetIonCoor(GOAtoms* goIon, int IonType,int iPos)
{
	goIon->r[0][0]=PMFxyz[0][IonType][iPos];
	goIon->r[1][0]=PMFxyz[1][IonType][iPos];
	goIon->r[2][0]=PMFxyz[2][IonType][iPos];
	return EXIT_SUCCESS;
}
int PMFCalculation::SetIonCoor(int IonType,int iPos)
{
	GOAtoms* goIon=IonsForPMF[IonType];
	goIon->r[0][0]=PMFxyz[0][IonType][iPos];
	goIon->r[1][0]=PMFxyz[1][IonType][iPos];
	goIon->r[2][0]=PMFxyz[2][IonType][iPos];
	
	/*
	GOAtoms* goIon=IonsForPMF[IonType];
	int CurrentIonPos=PMFPnts[IonType][iPos];
	int rsent[3],i,gridPoint;
	int gsX=m_ContWorld->GridSize[0];
	int gsXY=m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1];
	rsent[0]=CurrentIonPos%gsX;
	rsent[2]=CurrentIonPos/gsXY;
	rsent[1]=(CurrentIonPos-rsent[2]*gsXY)/gsX;
		
	gridPoint=rsent[0]+rsent[1]*gsX+rsent[2]*gsXY;
	if(gridPoint!=CurrentIonPos)fprintf(stdout,"ERROR for some reason gridPoint!=r[CurrentIon]\n");
		
	goIon->r[0][0]=rsent[0];
	goIon->r[1][0]=rsent[1];
	goIon->r[2][0]=rsent[2];*/
	return EXIT_SUCCESS;
}

BuildWorldNI* PMFCalculation::CmdBuildWorld(const TiXmlElement *Elt,ContWorld* _ContWorld)
{
	BuildWorldNI* _BuildWorldNI;
	if(_ContWorld!=NULL)
	{
		_BuildWorldNI=new BuildWorldNI();
		_BuildWorldNI->LoadXML(Elt);
	}
	else
	{
		fprintf(stderr,"ERROR: There is no grid set yet\n");
	}
	return _BuildWorldNI;
}
int PMFCalculation::LoadGion(const char * filename)
{
// 	int i;
// 	float fpoh= 4*M_PI*m_ContWorld->GridScale;
// 	
// 	pnpPrint0("Load Qst and Phi for Gion from %s\n",filename);
// 	
// 	FILE *in=fopen(filename,"rb");
// 	PNP_EXIT_FAIL_NULL(in,"Can't read file\n");
// 	
// 	fread(&GproteinQnum,sizeof(int),1,in);
// 	pnpPrint0("\t%d points to read\n",GproteinQnum);
// 	
// 	GionQst = new float[GproteinQnum];
// 	GionPhi = new float[GproteinQnum];
// 	GionQndx= new int[GproteinQnum];
// 	
// 	fread(GionQndx,sizeof(int),GproteinQnum,in);
// 	fread(GionQst,sizeof(float),GproteinQnum,in);
// 	fread(GionPhi,sizeof(float),GproteinQnum,in);
// 	fclose(in);
// 	
// 	for(i=0;i<GproteinQnum;i++)
// 	{
// 		GionQst[i]*=fpoh;
// 	}
// 	GproteinLoad=true;
	return EXIT_SUCCESS;
}
int PMFCalculation::LoadGprotein(const char * filename)
{
	PNP_EXIT_FAIL_NULL(m_ContWorld, "m_ContWorld is not initialize yet\n");
	
	int i;
	float fpoh= 4*M_PI*m_ContWorld->GridScale;
	
	pnpPrint0("Load Qst and Phi for Gprotein from %s\n",filename);
	
	if(pnpsapp->GetMyAbsRank()==pnpsapp->GetMaster())
	{
		FILE *in=fopen(filename,"rb");
		PNP_EXIT_FAIL_NULL(in,"Can't read file\n");
		
		fread(&GproteinQnum,sizeof(int),1,in);
		pnpPrint0("\t%d points to read\n",GproteinQnum);
		
		GproteinQst = new float[GproteinQnum];
		GproteinPhi = new float[GproteinQnum];
		GproteinQndx= new int[GproteinQnum];
		
		fread(GproteinQndx,sizeof(int),GproteinQnum,in);
		fread(GproteinQst,sizeof(float),GproteinQnum,in);
		fread(GproteinPhi,sizeof(float),GproteinQnum,in);
		fclose(in);
		
		for(i=0;i<GproteinQnum;i++)
		{
			GproteinQst[i]*=fpoh;
		}
		GproteinLoad=true;
#ifdef MPI_PARALLEL
		pnpsapp->MyComSuperGroup.Bcast(&GproteinQnum, 1, MPI::INT, pnpsapp->GetMaster());
		
		pnpsapp->MyComSuperGroup.Bcast(GproteinQst, GproteinQnum, MPI::FLOAT, pnpsapp->GetMaster());
		pnpsapp->MyComSuperGroup.Bcast(GproteinPhi, GproteinQnum, MPI::FLOAT, pnpsapp->GetMaster());
		pnpsapp->MyComSuperGroup.Bcast(GproteinQndx, GproteinQnum, MPI::INT, pnpsapp->GetMaster());
#endif
	}
	else
	{
#ifdef MPI_PARALLEL
		pnpsapp->MyComSuperGroup.Bcast(&GproteinQnum, 1, MPI::INT, pnpsapp->GetMaster());
		pnpPrint0("\t%d points to read\n",GproteinQnum);
		
		GproteinQst = new float[GproteinQnum];
		GproteinPhi = new float[GproteinQnum];
		GproteinQndx= new int[GproteinQnum];
		
		pnpsapp->MyComSuperGroup.Bcast(GproteinQst, GproteinQnum, MPI::FLOAT, pnpsapp->GetMaster());
		pnpsapp->MyComSuperGroup.Bcast(GproteinPhi, GproteinQnum, MPI::FLOAT, pnpsapp->GetMaster());
		pnpsapp->MyComSuperGroup.Bcast(GproteinQndx, GproteinQnum, MPI::INT, pnpsapp->GetMaster());
#endif
		GproteinLoad=true;
	}
	DbgPrint2("PMFCalculation::LoadGprotein DONE\n");
	return EXIT_SUCCESS;
}
int PMFCalculation::LoadGproteinRough(const char * filename)
{
	PNP_EXIT_FAIL_NULL(m_ContWorldRough, "m_ContWorldRough is not initialize yet\n");
	
	int i;
	float fpoh= 4*M_PI*m_ContWorldRough->GridScale;
	
	pnpPrint0("Load Qst and Phi for Gprotein from %s\n",filename);
	
	if(pnpsapp->GetMyAbsRank()==pnpsapp->GetMaster())
	{
		FILE *in=fopen(filename,"rb");
		PNP_EXIT_FAIL_NULL(in,"Can't read file\n");
		
		fread(&GproteinRoughQnum,sizeof(int),1,in);
		pnpPrint0("\t%d points to read\n",GproteinRoughQnum);
		
		GproteinRoughQst = new float[GproteinRoughQnum];
		GproteinRoughPhi = new float[GproteinRoughQnum];
		GproteinRoughQndx= new int[GproteinRoughQnum];
		
		fread(GproteinRoughQndx,sizeof(int),GproteinRoughQnum,in);
		fread(GproteinRoughQst,sizeof(float),GproteinRoughQnum,in);
		fread(GproteinRoughPhi,sizeof(float),GproteinRoughQnum,in);
		fclose(in);
		
		for(i=0;i<GproteinRoughQnum;i++)
		{
			GproteinRoughQst[i]*=fpoh;
		}
		GproteinRoughLoad=true;
#ifdef MPI_PARALLEL
		pnpsapp->MyComSuperGroup.Bcast(&GproteinRoughQnum, 1, MPI::INT, pnpsapp->GetMaster());
		
		pnpsapp->MyComSuperGroup.Bcast(GproteinRoughQst, GproteinRoughQnum, MPI::FLOAT, pnpsapp->GetMaster());
		pnpsapp->MyComSuperGroup.Bcast(GproteinRoughPhi, GproteinRoughQnum, MPI::FLOAT, pnpsapp->GetMaster());
		pnpsapp->MyComSuperGroup.Bcast(GproteinRoughQndx, GproteinRoughQnum, MPI::INT, pnpsapp->GetMaster());
#endif
	}
	else
	{
#ifdef MPI_PARALLEL
		pnpsapp->MyComSuperGroup.Bcast(&GproteinRoughQnum, 1, MPI::INT, pnpsapp->GetMaster());
		pnpPrint0("\t%d points to read\n",GproteinRoughQnum);
		
		GproteinRoughQst = new float[GproteinQnum];
		GproteinRoughPhi = new float[GproteinQnum];
		GproteinRoughQndx= new int[GproteinQnum];
		
		pnpsapp->MyComSuperGroup.Bcast(GproteinRoughQst, GproteinRoughQnum, MPI::FLOAT, pnpsapp->GetMaster());
		pnpsapp->MyComSuperGroup.Bcast(GproteinRoughPhi, GproteinRoughQnum, MPI::FLOAT, pnpsapp->GetMaster());
		pnpsapp->MyComSuperGroup.Bcast(GproteinRoughQndx, GproteinRoughQnum, MPI::INT, pnpsapp->GetMaster());
#endif
		GproteinRoughLoad=true;
	}
	DbgPrint2("PMFCalculation::LoadGprotein DONE\n");
	return EXIT_SUCCESS;
}
int PMFCalculation::LoadIonsForPMF(const TiXmlElement *Elt)
{
	const TiXmlElement *CldElt;
	int i;
	NIons=0;
	CldElt=Elt->FirstChildElement("AtomsParameters");
	while(CldElt!=NULL)
	{
		NIons++;
		CldElt=CldElt->NextSiblingElement("AtomsParameters");
	}
	
	
	pnpPrint0("NIons=%d\n",NIons);
	IonsForPMF=new GOAtoms*[NIons];
	
	i=0;
	CldElt=Elt->FirstChildElement("AtomsParameters");
	if(CldElt!=NULL)
	{
		do
		{
			GOAtoms* Atoms=new GOAtoms();
			Atoms->LoadXML(CldElt);
			IonsForPMF[i]=Atoms;
			CldElt=CldElt->NextSiblingElement("AtomsParameters");
			i++;
		}
		while(CldElt!=NULL);
	}
	Npnts=new int[NIons];
	PMFPnts=new int*[NIons];
	PMF=new double*[NIons];
	PMFCalc=new bool*[NIons];
	PMFxyz[0]=new float*[NIons];
	PMFxyz[1]=new float*[NIons];
	PMFxyz[2]=new float*[NIons];
	PMFPntsStart=new int[NIons];
	for(i=0;i<NIons;i++)
	{
		Npnts[i]=NULL;
		PMFPnts[i]=NULL;
		PMF[i]=NULL;
		PMFCalc[i]=NULL;
		PMFxyz[0][i]=NULL;
		PMFxyz[1][i]=NULL;
		PMFxyz[2][i]=NULL;
		PMFPntsStart[i]=0;
		
	}
	return EXIT_SUCCESS;
}
int PMFCalculation::ReadPointsFromPDB(int ion,const char *filename)
{
	int k;
	std::vector<float> r[3];
	FILE* fpdb = fopen(filename,"rt");
	PNP_EXIT_FAIL_NULL1(fpdb,"Can't read %s\n",filename);
	char Record[1024];
	while( fgets(Record,1024,fpdb) )
	{
		if(strncmp("ATOM",Record,4)==0)
		{
			//float tmpx,tmpy,tmpz,tmpq,tmpr;
			//sscanf(Record+27,"%f %f %f %f %f",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			double tmpx,tmpy,tmpz,tmpq,tmpr;
			sscanf(Record+27,"%lf %lf %lf",&tmpx,&tmpy,&tmpz);
				//pnpPrint0("\t%f %f %f %f %f\n",tmpx,tmpy,tmpz,tmpq,tmpr);
			r[0].push_back(float(tmpx));
			r[1].push_back(float(tmpy));
			r[2].push_back(float(tmpz));
		}
	}
	fclose(fpdb);
	pnpPrint0("%d atoms was read from %s for ion %d\n",r[0].size(),filename,ion);
	Npnts[ion]=r[0].size();
	PMFxyz[0][ion]=new float[Npnts[ion]];
	PMFxyz[1][ion]=new float[Npnts[ion]];
	PMFxyz[2][ion]=new float[Npnts[ion]];
	for(k=0;k<Npnts[ion];k++)
	{
		PMFxyz[0][ion][k]=r[0][k];
		PMFxyz[1][ion][k]=r[1][k];
		PMFxyz[2][ion][k]=r[2][k];
	}
	r[0].clear();
	r[1].clear();
	r[2].clear();
	
	//convert coordinates to internal, and also calculate grid point(closest)
	PMFPnts[ion]=new int[Npnts[ion]];
	for(k=0;k<Npnts[ion];k++)
	{
		PMFPnts[ion][k]=-1;
	}
	
	PNP_EXIT_FAIL_NULL(m_ContWorld, "m_ContWorld is not initialize yet\n");
	int j;
	int irc[3];
	float frc[3];
	for(j=0;j<3;j++)
	{
		irc[j]=m_ContWorld->GridSizeGlobal[j]/2;
		frc[j]=irc[j];
	}
	float interR[3];
	for(k=0;k<Npnts[ion];k++)
	{
		PMFPnts[ion][k]=m_ContWorld->ConvGlobalExternalXYZToGrdPnt(PMFxyz[0][ion][k], PMFxyz[1][ion][k], PMFxyz[2][ion][k]);
	}
	
	return EXIT_SUCCESS;
}
int PMFCalculation::LoadPointsForPMF(const TiXmlElement *Elt)
{
	int i,j,k;
	int NionsAtPnts;
// 	DbgPrint0("PMFPntsType=%s\n",Elt->Attribute("PntsType"));
// 	//Define which type of points we have
// 	if(Elt->Attribute("PntsType")!=NULL)
// 	{
// 		if(strcmp("grd",Elt->Attribute("PntsType"))==0)
// 			PMFPntsType=PMFPntsTypeGrid;
// 		else if(strcmp("xyz",Elt->Attribute("PntsType"))==0)
// 			PMFPntsType=PMFPntsTypeXYZ;
// 		else
// 			PMFPntsType=PMFPntsTypeGrid;
// 	}
// 	else
// 	{
		PMFPntsType=PMFPntsTypePDB;
// 	}
	//Read points
	/*if(PMFPntsType==PMFPntsTypeGrid)
	{
		NionsAtPnts=Elt->GetNumOfElement("Points");
		const TiXmlElement *PointsElt=Elt->FirstChildElement("Points");
		for(i=0;i<NionsAtPnts;i++)
		{
			int NpntsAtPnts;
			int *ipnt=NULL;
			PointsElt->GetArrOfIntFromText(&ipnt,&NpntsAtPnts);
			Npnts[i]=NpntsAtPnts;
			PMFPnts[i]=ipnt;
			//get coordinates
			int gsX=m_ContWorld->GridSize[0];
			int gsXY=m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1];
			
			for(k=0;k<Npnts[i];k++)
			{
				PMFxyz[0][i][k]=PMFPnts[i][k]%gsX;
				PMFxyz[2][i][k]=PMFPnts[i][k]/gsXY;
				PMFxyz[1][i][k]=(PMFPnts[i][k]-PMFxyz[2][i][k]*gsXY)/gsX;
			}
			PointsElt=PointsElt->NextSiblingElement("Points");
		}
	}
	if(PMFPntsType==PMFPntsTypeXYZ)
	{
		int irc[3];
		float frc[3];
		for(j=0;j<3;j++)
		{
			irc[j]=m_ContWorld->GridSizeGlobal[j]/2;
			frc[j]=irc[j];
		}
		NionsAtPnts=Elt->GetNumOfElement("Points");
		const TiXmlElement *PointsElt=Elt->FirstChildElement("Points");
		for(i=0;i<NionsAtPnts;i++)
		{
			int NpntsAtPnts;
			int *ipnt=NULL;
			PointsElt->GetArrOfFloatCoordFromText(&PMFxyz[0][i],&PMFxyz[1][i],&PMFxyz[2][i],&NpntsAtPnts);
			for(k=0;k<NpntsAtPnts;k++)
			{
				pnpPrint0("%d %f %f %f\n",i,PMFxyz[0][i][k],PMFxyz[1][i][k],PMFxyz[2][i][k]);
			}
			Npnts[i]=NpntsAtPnts;
			
			//convert coordinates to internal, and also calculate grid point(closest)
			PMFPnts[i]=new int[Npnts[i]];
			for(k=0;k<Npnts[i];k++)
			{
				for(j=0;j<3;j++)
					PMFxyz[j][i][k]=PMFxyz[j][i][k]*m_ContWorld->GridScale+frc[j];
				PMFPnts[i][k]=int((PMFxyz[0][i][k]+0.5))+int((PMFxyz[1][i][k]+0.5))*m_ContWorld->GS_X+int((PMFxyz[2][i][k]+0.5))*m_ContWorld->GS_XY;
			}
			
			PointsElt=PointsElt->NextSiblingElement("Points");
		}
		
		
	}*/
	if(PMFPntsType==PMFPntsTypePDB)
	{
		NionsAtPnts=0;
		
		const char *filename;
		filename=Elt->Attribute("PDB");
		if(filename!=NULL)
		{
			NionsAtPnts=1;
			ReadPointsFromPDB(0,filename);
		}
		else
		{
			filename=Elt->Attribute("PDB0");
			if(filename!=NULL)
			{
				NionsAtPnts=1;
				ReadPointsFromPDB(0,filename);
			}
			filename=Elt->Attribute("PDB1");
			if(filename!=NULL)
			{
				NionsAtPnts++;
				ReadPointsFromPDB(1,filename);
			}
		}
	}
	//If both ions have the same set of points
	if(NionsAtPnts<NIons)
	{
		pnpWarning("NionsAtPnts<NIons will use same set of points for all ions\n");
		for(i=NionsAtPnts;i<NIons;i++)
		{
			Npnts[i]=Npnts[0];
			PMFPnts[i]=new int[Npnts[i]];
			for(j=0;j<3;j++)
				PMFxyz[j][i]=new float[Npnts[i]];
			
			for(j=0;j<Npnts[i];j++)
				PMFPnts[i][j]=PMFPnts[0][j];
			for(k=0;k<Npnts[i];k++)
				for(j=0;j<3;j++)
					PMFxyz[j][i][k]=PMFxyz[j][0][k];
		}
	}
	//Set initial values for results
	for(i=0;i<NIons;i++)
	{
		PMF[i]=new double[Npnts[i]];
		PMFCalc[i]=new bool[Npnts[i]];
	}
	for(i=0;i<NIons;i++)
	{
		for(j=0;j<Npnts[i];j++)
		{
			PMFCalc[i][j]=false;
			PMF[i][j]=0.0;
		}
	}
	return EXIT_SUCCESS;
}
int PMFCalculation::BuildContWorld(ContWorld* world,BuildWorldNI* _BuildWorldNI)
{
	DbgPrint1("BuildWorldNI::BuildContWorld()\n");
	int i,j,k;
	float Displacement[3]={0.0,0.0,0.0};
	float* Surf=new float[world->GridSizeXYZGlobal*4];
	int* iVtmp=new int[world->GridSizeXYZGlobal];
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	int iElt;
	//Change Parameters for internal 
// 	_BuildWorldNI->ChangeUnits(world,true);
// 	fprintf(stdout,"IonsForPMF[0]->r=[%f,%f,%f]\n",IonsForPMF[0]->r[0][0],IonsForPMF[0]->r[1][0],IonsForPMF[0]->r[2][0]);
// 	GOAtoms* atms=(GOAtoms*)_BuildWorldNI->GOElms[0];
// 	fprintf(stdout,"atms[0]->r=[%f,%f,%f]\n",atms->r[0][0],atms->r[1][0],atms->r[2][0]);
// 	for(iElt=0;iElt<_BuildWorldNI->GOElms.size();iElt++)
// 	{
// 		_BuildWorldNI->GOElms[iElt]->ChangeUnits(world,true);
// 	}
//	fprintf(stdout,"IonsForPMF[0]->r=[%f,%f,%f]\n",IonsForPMF[0]->r[0][0],IonsForPMF[0]->r[1][0],IonsForPMF[0]->r[2][0]);
//	fprintf(stdout,"atms[0]->r=[%f,%f,%f]\n",atms->r[0][0],atms->r[1][0],atms->r[2][0]);
	NodeIndexing *NIndexing;
	DeleteObjByPnt(world->NIndexing);
	world->NIndexing=new NodeIndexing();
	NIndexing=world->NIndexing;
	
	unsigned int GS[3]={world->GridSizeGlobal[0],world->GridSizeGlobal[1],world->GridSizeGlobal[2]};
	NIndexing->SetNNodes(GS,world->GridScale);
	//Set Values Table
	
	for(i=0;i<_BuildWorldNI->DielNum;i++)NIndexing->Eps[i]=_BuildWorldNI->DielConst[i];
	
	
	for(i=0;i<_BuildWorldNI->DiffusionNum;i++)
	{
		NIndexing->D[i]=_BuildWorldNI->DiffusionConst[i];
		if(_BuildWorldNI->DiffusionConst[i]==0.0)NIndexing->C[i]=0.0;
		else NIndexing->C[i]=_BuildWorldNI->BulkParam->C[0];
	}
	//
	
	for(i=0;i<_BuildWorldNI->DielNum;i++)DbgPrint1("DielConst[%d]=%f\n",i,_BuildWorldNI->DielConst[i]);
	for(i=0;i<3;i++)
	{
		DbgPrint1("Build Dielectric Map %d\n",i);
		for(j=0;j<world->GridSizeXYZGlobal;j++)
			iVtmp[j]=_BuildWorldNI->BulkParam->Epsilon;

		for(j=0;j<4*world->GridSizeXYZGlobal;j++)Surf[j]=-BIGDISTANSE;
		for(j=3;j<4*world->GridSizeXYZGlobal;j=j+4)Surf[j]=BIGDISTANSE;
		
		Displacement[0]=0.0;
		Displacement[1]=0.0;
		Displacement[2]=0.0;
		Displacement[i]=-0.5;
		
		for(iElt=0;iElt<_BuildWorldNI->GOElms.size();iElt++)
		{
			if(_BuildWorldNI->GOElms[iElt]->Epsilon!=_BuildWorldNI->BulkParam->Epsilon)
				_BuildWorldNI->GOElms[iElt]->BuildPreMaps(_BuildWorldNI,world,iVtmp,_BuildWorldNI->GOElms[iElt]->Epsilon,_BuildWorldNI->BulkParam->Epsilon,Displacement,0.0,_BuildWorldNI->Rwat,Surf,0);
		}
		
		_BuildWorldNI->FinalazeSEV(world,iVtmp,_BuildWorldNI->BulkParam->Epsilon,_BuildWorldNI->Rwat,Surf);
		//Move indexes to 0 based
		for(j=0;j<world->GridSizeXYZGlobal;j++)
			iVtmp[j]--;
		switch(i)
		{
			case 0:
				NIndexing->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon0,NodeIndexing::Epsilon0Sft);
				break;
			case 1:
				NIndexing->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon1,NodeIndexing::Epsilon1Sft);
				break;
			case 2:
				NIndexing->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon2,NodeIndexing::Epsilon2Sft);
				break;
		}
		
		
	}
	
	NIndexing->CalcDielBoarder();
	fprintf(stdout,"iVtmp %p %d %d\n",iVtmp,iVtmp[0],iVtmp[100]);
	delete [] iVtmp;
	delete [] Surf;
	
	//charges
	float* Qstat = new float[world->GridSizeXYZGlobal];
	for(i=0;i<world->GridSizeXYZGlobal;i++)Qstat[i]=0.0;

	for(iElt=0;iElt<_BuildWorldNI->GOElms.size();iElt++)
	{
		_BuildWorldNI->GOElms[iElt]->BuildPreMapsCharges(_BuildWorldNI,world,Qstat);
	}
	double Q=0.0;
	for(i=0;i<world->GridSizeXYZGlobal;i++)Q+=Qstat[i];
	NIndexing->SetChargeMapFromArray(Qstat);
	delete [] Qstat;
	DbgPrint1("Charge of Map: %.6g [e]\n",Q/(4.0*M_PI*world->GridScale));
	NIndexing->PBC[0]=world->PBC[0];
	NIndexing->PBC[1]=world->PBC[1];
	NIndexing->PBC[2]=world->PBC[2];
	
// 	//Change Parameters back to extenal then needed
// 	_BuildWorldNI->ChangeUnits(world,false);
// 	for(iElt=0;iElt<_BuildWorldNI->GOElms.size();iElt++)
// 	{
// 		_BuildWorldNI->GOElms[iElt]->ChangeUnits(world,false);
// 	}
	return EXIT_SUCCESS;
}

int PMFCalculation::PreBuildContWorld(ContWorld* world,BuildWorldNI* _BuildWorldNI)
{
	DbgPrint1("BuildWorldNI::BuildContWorld()\n");
	int i,j,k;
	float Displacement[3]={0.0,0.0,0.0};
	
	m_Surf=new float*[3];
	m_iVtmp=new int*[3];
	for(i=0;i<3;i++)
	{
		m_Surf[i]=new float[world->GridSizeXYZGlobal*4];
		m_iVtmp[i]=new int[world->GridSizeXYZGlobal];
	}
	m_Qstat = new float[world->GridSizeXYZGlobal];
	
	float* Surf=NULL;
	int* iVtmp=NULL;
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	int iElt;
	

	
	for(i=0;i<_BuildWorldNI->DielNum;i++)
		DbgPrint1("DielConst[%d]=%f\n",i,_BuildWorldNI->DielConst[i]);
	for(i=0;i<3;i++)
	{
		Surf=m_Surf[i];
		iVtmp=m_iVtmp[i];
		DbgPrint1("Build Dielectric Map %d\n",i);
		for(j=0;j<world->GridSizeXYZGlobal;j++)
			iVtmp[j]=_BuildWorldNI->BulkParam->Epsilon;

		for(j=0;j<4*world->GridSizeXYZGlobal;j++)Surf[j]=-BIGDISTANSE;
		for(j=3;j<4*world->GridSizeXYZGlobal;j=j+4)Surf[j]=BIGDISTANSE;
		
		Displacement[0]=0.0;
		Displacement[1]=0.0;
		Displacement[2]=0.0;
		Displacement[i]=-0.5;
		
		for(iElt=0;iElt<_BuildWorldNI->GOElms.size();iElt++)
		{
			if(_BuildWorldNI->GOElms[iElt]->Epsilon!=_BuildWorldNI->BulkParam->Epsilon)
				_BuildWorldNI->GOElms[iElt]->BuildPreMaps(_BuildWorldNI,world,iVtmp,_BuildWorldNI->GOElms[iElt]->Epsilon,_BuildWorldNI->BulkParam->Epsilon,Displacement,0.0,_BuildWorldNI->Rwat,Surf,0);
		}
	}
	//charges
	
	
	float* Qstat = m_Qstat;
	for(i=0;i<world->GridSizeXYZGlobal;i++)Qstat[i]=0.0;

	for(iElt=0;iElt<_BuildWorldNI->GOElms.size();iElt++)
	{
		_BuildWorldNI->GOElms[iElt]->BuildPreMapsCharges(_BuildWorldNI,world,Qstat);
	}
	double Q=0.0;
	for(i=0;i<world->GridSizeXYZGlobal;i++)Q+=Qstat[i];
	DbgPrint1("Charge of Map: %.6g [e]\n",Q/(4.0*M_PI*world->GridScale));
	
	//Change Parameters back to extenal then needed
	
	DbgPrint1("PreBuildContWorld done\n");
	return EXIT_SUCCESS;
}
int PMFCalculation::PreBuildContWorld1(ContWorld* world,BuildWorldNI* _BuildWorldNI)
{
	DbgPrint1("BuildWorldNI::BuildContWorld()\n");
	int i,j,k;
	float Displacement[3]={0.0,0.0,0.0};
	
	float* Surf=new float[world->GridSizeXYZGlobal*4];
	int* iVtmp=new int[world->GridSizeXYZGlobal];
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	int iElt;
	
	for(i=0;i<_BuildWorldNI->DielNum;i++)
		DbgPrint1("DielConst[%d]=%f\n",i,_BuildWorldNI->DielConst[i]);
	for(i=0;i<3;i++)
	{
		DbgPrint1("Build Dielectric Map %d\n",i);
		for(j=0;j<world->GridSizeXYZGlobal;j++)
			iVtmp[j]=_BuildWorldNI->BulkParam->Epsilon;

		for(j=0;j<4*world->GridSizeXYZGlobal;j++)Surf[j]=-BIGDISTANSE;
		for(j=3;j<4*world->GridSizeXYZGlobal;j=j+4)Surf[j]=BIGDISTANSE;
		
		Displacement[0]=0.0;
		Displacement[1]=0.0;
		Displacement[2]=0.0;
		Displacement[i]=-0.5;
		
		for(iElt=0;iElt<_BuildWorldNI->GOElms.size();iElt++)
		{
			if(_BuildWorldNI->GOElms[iElt]->Epsilon!=_BuildWorldNI->BulkParam->Epsilon)
				_BuildWorldNI->GOElms[iElt]->BuildPreMaps(_BuildWorldNI,world,iVtmp,_BuildWorldNI->GOElms[iElt]->Epsilon,_BuildWorldNI->BulkParam->Epsilon,Displacement,0.0,_BuildWorldNI->Rwat,Surf,0);
		}
		//Save Surf and iVtmp
		IdSurf[i]=pnpsapp->DropArray((Bytef*)Surf,4*world->GridSizeXYZGlobal*sizeof(float));
		IdiVtmp[i]=pnpsapp->DropArray((Bytef*)iVtmp,world->GridSizeXYZGlobal*sizeof(int));
	}
	delete [] Surf;
	delete [] iVtmp;
	//charges	
	float* Qstat = new float[world->GridSizeXYZGlobal];
	for(i=0;i<world->GridSizeXYZGlobal;i++)Qstat[i]=0.0;

	for(iElt=0;iElt<_BuildWorldNI->GOElms.size();iElt++)
	{
		_BuildWorldNI->GOElms[iElt]->BuildPreMapsCharges(_BuildWorldNI,world,Qstat);
	}
	double Q=0.0;
	for(i=0;i<world->GridSizeXYZGlobal;i++)Q+=Qstat[i];
	DbgPrint1("Charge of Map: %.6g [e]\n",Q/(4.0*M_PI*world->GridScale));
	
	
	//Save Qstat
	IdQstat=pnpsapp->DropArray((Bytef*)Qstat,world->GridSizeXYZGlobal*sizeof(float));
	delete [] Qstat;	
	//Change Parameters back to extenal then needed
	
	DbgPrint1("PreBuildContWorld done\n");
	return EXIT_SUCCESS;
}
int PMFCalculation::FinalizePreBuildContWorld(ContWorld* world,BuildWorldNI* _BuildWorldNI,GOAtoms* goIon)
{
	DbgPrint1("BuildWorldNI::FinalizePreBuildContWorld()\n");
	int i,j,k;
	float Displacement[3]={0.0,0.0,0.0};
	float* Surf=new float[world->GridSizeXYZGlobal*4];
	int* iVtmp=new int[world->GridSizeXYZGlobal];
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	int iElt;

	NodeIndexing *NIndexing;
	DeleteObjByPnt(world->NIndexing);
	world->NIndexing=new NodeIndexing();
	NIndexing=world->NIndexing;
	
	unsigned int GS[3]={world->GridSizeGlobal[0],world->GridSizeGlobal[1],world->GridSizeGlobal[2]};
	NIndexing->SetNNodes(GS,world->GridScale);
	//Set Values Table
	
	for(i=0;i<_BuildWorldNI->DielNum;i++)NIndexing->Eps[i]=_BuildWorldNI->DielConst[i];
	
	
	for(i=0;i<_BuildWorldNI->DiffusionNum;i++)
	{
		NIndexing->D[i]=_BuildWorldNI->DiffusionConst[i];
		if(_BuildWorldNI->DiffusionConst[i]==0.0)NIndexing->C[i]=0.0;
		else NIndexing->C[i]=_BuildWorldNI->BulkParam->C[0];
	}
	//
	
	for(i=0;i<_BuildWorldNI->DielNum;i++)
		DbgPrint1("DielConst[%d]=%f\n",i,_BuildWorldNI->DielConst[i]);
	for(i=0;i<3;i++)
	{
		DbgPrint1("Build Dielectric Map %d\n",i);
		for(j=0;j<world->GridSizeXYZGlobal;j++)
			iVtmp[j]=m_iVtmp[i][j];

		for(j=0;j<4*world->GridSizeXYZGlobal;j++)
			Surf[j]=m_Surf[i][j];
		
		Displacement[0]=0.0;
		Displacement[1]=0.0;
		Displacement[2]=0.0;
		Displacement[i]=-0.5;
		
		goIon->BuildPreMaps(_BuildWorldNI,world,iVtmp,goIon->Epsilon,_BuildWorldNI->BulkParam->Epsilon,Displacement,0.0,_BuildWorldNI->Rwat,Surf,0);
		
		_BuildWorldNI->FinalazeSEV(world,iVtmp,_BuildWorldNI->BulkParam->Epsilon,_BuildWorldNI->Rwat,Surf);
		//Move indexes to 0 based
		for(j=0;j<world->GridSizeXYZGlobal;j++)
			iVtmp[j]--;
		switch(i)
		{
			case 0:
				NIndexing->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon0,NodeIndexing::Epsilon0Sft);
				break;
			case 1:
				NIndexing->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon1,NodeIndexing::Epsilon1Sft);
				break;
			case 2:
				NIndexing->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon2,NodeIndexing::Epsilon2Sft);
				break;
		}
		
		
	}
	
	NIndexing->CalcDielBoarder();
	fprintf(stdout,"iVtmp %p %d %d\n",iVtmp,iVtmp[0],iVtmp[100]);
	delete [] iVtmp;
	delete [] Surf;
	
	//charges
	float* Qstat = new float[world->GridSizeXYZGlobal];
	for(i=0;i<world->GridSizeXYZGlobal;i++)Qstat[i]=m_Qstat[i];

	goIon->BuildPreMapsCharges(_BuildWorldNI,world,Qstat);
	double Q=0.0;
	for(i=0;i<world->GridSizeXYZGlobal;i++)Q+=Qstat[i];
	NIndexing->SetChargeMapFromArray(Qstat);
	delete [] Qstat;
	DbgPrint1("Charge of Map: %.6g [e]\n",Q/(4.0*M_PI*world->GridScale));
	NIndexing->PBC[0]=world->PBC[0];
	NIndexing->PBC[1]=world->PBC[1];
	NIndexing->PBC[2]=world->PBC[2];
	
	//init grdpntQst and Gsubtract
	DeleteCArray(grdpntQst);
	DeleteCArray(Gsubtract);
	grdpntQst=world->NIndexing->GetChargeIndex();
	if(GproteinLoad||GionLoad)
	{
		if(Gsubtract==NULL)
		{
			Gsubtract=new double[world->NIndexing->QNum];
			for(k=0;k<world->NIndexing->QNum;k++)
				Gsubtract[k]=0.0;
		}
		
		if(GproteinLoad)
		{
			if(GproteinQnum>world->NIndexing->QNum-1)
			{
				pnpWarning("Number of partial charges in Gprotein is bigger then partial charges in protein-ion - 1\n");
			}
			double qphi;
			int iGprotein;
			int count=0;
			for(j=0;j<GproteinQnum;j++)
			{
				while(GproteinQndx[j]!=grdpntQst[count])
				{
					count++;
				}
				qphi=0.0;
				iGprotein=-1;
				qphi=GproteinQst[j]*GproteinPhi[j];
				Gsubtract[count]=-qphi;
			}
		}
	}
	
	return EXIT_SUCCESS;
}
int PMFCalculation::FinalizePreBuildContWorld1(ContWorld* world,BuildWorldNI* _BuildWorldNI,GOAtoms* goIon)
{
	DbgPrint1("BuildWorldNI::FinalizePreBuildContWorld1()\n");
	int i,j,k;
	float Displacement[3]={0.0,0.0,0.0};
	float* Surf=new float[world->GridSizeXYZGlobal*4];
	int* iVtmp=new int[world->GridSizeXYZGlobal];
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	int iElt;

	NodeIndexing *NIndexing;
	DeleteObjByPnt(world->NIndexing);
	world->NIndexing=new NodeIndexing();
	NIndexing=world->NIndexing;
	
	unsigned int GS[3]={world->GridSizeGlobal[0],world->GridSizeGlobal[1],world->GridSizeGlobal[2]};
	NIndexing->SetNNodes(GS,world->GridScale);
	//Set Values Table
	
	for(i=0;i<_BuildWorldNI->DielNum;i++)NIndexing->Eps[i]=_BuildWorldNI->DielConst[i];
	
	
	for(i=0;i<_BuildWorldNI->DiffusionNum;i++)
	{
		NIndexing->D[i]=_BuildWorldNI->DiffusionConst[i];
		if(_BuildWorldNI->DiffusionConst[i]==0.0)NIndexing->C[i]=0.0;
		else NIndexing->C[i]=_BuildWorldNI->BulkParam->C[0];
	}
	//
	
	for(i=0;i<_BuildWorldNI->DielNum;i++)
		DbgPrint1("DielConst[%d]=%f\n",i,_BuildWorldNI->DielConst[i]);
	for(i=0;i<3;i++)
	{
		DbgPrint1("Build Dielectric Map %d\n",i);
		//for(j=0;j<world->GridSizeXYZGlobal;j++)
		//	iVtmp[j]=m_iVtmp[i][j];
		pnpsapp->PullArray((Bytef*)iVtmp,IdiVtmp[i]);
		pnpsapp->PullArray((Bytef*)Surf,IdSurf[i]);
		//for(j=0;j<4*world->GridSizeXYZGlobal;j++)
		//	Surf[j]=m_Surf[i][j];
		
		Displacement[0]=0.0;
		Displacement[1]=0.0;
		Displacement[2]=0.0;
		Displacement[i]=-0.5;
		
		goIon->BuildPreMaps(_BuildWorldNI,world,iVtmp,goIon->Epsilon,_BuildWorldNI->BulkParam->Epsilon,Displacement,0.0,_BuildWorldNI->Rwat,Surf,0);
		
		_BuildWorldNI->FinalazeSEV(world,iVtmp,_BuildWorldNI->BulkParam->Epsilon,_BuildWorldNI->Rwat,Surf);
		//Move indexes to 0 based
		for(j=0;j<world->GridSizeXYZGlobal;j++)
			iVtmp[j]--;
		switch(i)
		{
			case 0:
				NIndexing->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon0,NodeIndexing::Epsilon0Sft);
				break;
			case 1:
				NIndexing->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon1,NodeIndexing::Epsilon1Sft);
				break;
			case 2:
				NIndexing->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon2,NodeIndexing::Epsilon2Sft);
				break;
		}
		
		
	}
	
	NIndexing->CalcDielBoarder();
	fprintf(stdout,"iVtmp %p %d %d\n",iVtmp,iVtmp[0],iVtmp[100]);
	delete [] iVtmp;
	delete [] Surf;
	
	//charges
	float* Qstat = new float[world->GridSizeXYZGlobal];
	pnpsapp->PullArray((Bytef*)Qstat,IdQstat);

	goIon->BuildPreMapsCharges(_BuildWorldNI,world,Qstat);
	double Q=0.0;
	for(i=0;i<world->GridSizeXYZGlobal;i++)Q+=Qstat[i];
	NIndexing->SetChargeMapFromArray(Qstat);
	delete [] Qstat;
	DbgPrint1("Charge of Map: %.6g [e]\n",Q/(4.0*M_PI*world->GridScale));
	

	
	
	NIndexing->PBC[0]=world->PBC[0];
	NIndexing->PBC[1]=world->PBC[1];
	NIndexing->PBC[2]=world->PBC[2];
	
	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
PMFProcessing::PMFProcessing()
{
	
}
PMFProcessing::~PMFProcessing()
{
}
int PMFProcessing::InitZero()
{
	Nions=0;
	GridSize[0]=1;
	GridSize[1]=1;
	GridSize[2]=1;
	GS_X=1;GS_Y=1;GS_Z=1;GS_XY=1;GS_XYZ=1;
	GridScale=1.0;
		
	PMF=NULL;
	Eion=NULL;
	Eprot=0.0;
	m_ContWorld=NULL;
	EavrIon=NULL;
	Mask=NULL;
	MaskIon=NULL;
	return EXIT_SUCCESS;
}
int PMFProcessing::Clear()
{
	DeleteCVecArray(PMF,Nions);
	DeleteCVecArray(Mask,Nions);
	DeleteCVecArray(Eion,Nions);
	DeleteCVecArray(MaskIon,Nions);
	DeleteCArray(EavrIon);
	return EXIT_SUCCESS;
}
int PMFProcessing::SetContWorld(ContWorld* _ContWorld)
{
	PNP_EXIT_FAIL_NULL(_ContWorld,"PMFProcessing::SetContWorld: ContWorld==NULL\n");
	
	m_ContWorld=_ContWorld;
	GridSize[0]=m_ContWorld->GridSize[0];
	GridSize[1]=m_ContWorld->GridSize[1];
	GridSize[2]=m_ContWorld->GridSize[2];
	GS_X=GridSize[0];
	GS_Y=GridSize[1];
	GS_Z=GridSize[2];
	GS_XY=GS_X*GS_Y;
	GS_XYZ=GS_XY*GS_Z;
	GridScale=m_ContWorld->GridScale;
	Nions=m_ContWorld->NIonsTypes;
	return EXIT_SUCCESS;
}
int PMFProcessing::CmdPMFProcess(const TiXmlElement *Elt)
{
	int i,j;
	char FileIon[128]="\0";
	char FileProtIon[128]="\0";
	char FilePMFout[128]="\0";
	bool IonsHaveSameSize;
	float ConvFac;
	bool UseConvFac;
	bool PrintDisconvPoint;
	double dRstep;
	int Rstep;
	int StreachIter;
	int SmoothIter;
	bool ConvertPoissonPotToPMF;
	//Read Parameters
	Elt->GetCStrAttribute("FileProtIon",FileProtIon);
	Elt->GetCStrAttribute("FileIon",FileIon);
	if(Elt->GetIntAttribute("Nions",&Nions)!=EXIT_SUCCESS)
		Nions=m_ContWorld->NIonsTypes;
	Elt->GetDoubleAttribute("Eprot",&Eprot);
	if(Elt->GetBoolAttribute("IonsHaveSameSize",&IonsHaveSameSize)!=EXIT_SUCCESS)IonsHaveSameSize=false;
	//if(Elt->GetBoolAttribute("PrintDisconvPoint",&PrintDisconvPoint)!=EXIT_SUCCESS)PrintDisconvPoint=false;
	if(Elt->GetFloatAttribute("ConvFac",&ConvFac)==EXIT_SUCCESS)
		UseConvFac=true;
	else
	{
		UseConvFac=false;
		ConvFac=0.0;
	}
	if(Elt->GetDoubleAttribute("Rstep",&dRstep)==EXIT_SUCCESS)
	{
#if CC == pgCC
		Rstep=int(dRstep*GridScale+0.5);
#else
		Rstep=lround(dRstep*GridScale);
#endif

	}
	else
	{
		Rstep=0;
	}
	bool LinearInterpolation;
	if(Elt->GetBoolAttribute("LinearInterpolation",&LinearInterpolation)!=EXIT_SUCCESS)
		LinearInterpolation=true;
	
	if(Elt->GetIntAttribute("StreachIter",&StreachIter)!=EXIT_SUCCESS)
		StreachIter=0;
	if(Elt->GetIntAttribute("SmoothIter",&SmoothIter)!=EXIT_SUCCESS)
		SmoothIter=0;
	
	if(Elt->GetBoolAttribute("ConvertPoissonPotToPMF",&ConvertPoissonPotToPMF)!=EXIT_SUCCESS)ConvertPoissonPotToPMF=false;
	
	//Npnts=new int[Nions];
	PMF=new double*[Nions];
	Eion=new double*[Nions];
	EavrIon=new double[Nions];
	Mask=new int*[Nions];
	MaskIon=new int*[Nions];
	
	for(i=0;i<Nions;i++)
	{
		EavrIon[i]=0.0;
		PMF[i]=new double[GS_XYZ];
		Mask[i]=new int[GS_XYZ];
		Eion[i]=new double[GS_XYZ];
		MaskIon[i]=new int[GS_XYZ];
		for(j=0;j<GS_XYZ;j++)
		{
			PMF[i][j]=0.0;
			Mask[i][j]=0;
			Eion[i][j]=0.0;
			MaskIon[i][j]=0;
		}
	}
	pnpPrint0("Will read Eprot-ion\n");
	ReadEnergyFrom5ColumnFile(FileProtIon,PMF,Mask,ConvFac);
	pnpPrint0("Will read Eion\n");
	ReadEnergyFrom5ColumnFile(FileIon,Eion,MaskIon,ConvFac);
	
	CalculateAvrEion(IonsHaveSameSize);
	CombineEprotionEprotEion(IonsHaveSameSize);
	
	for(i=0;i<Nions;i++)
	{
		SetNbhoodToCalcNode(PMF[i],Mask[i],Rstep);
		if(LinearInterpolation)
			LineInterpolation0(PMF[i],Mask[i],4);
		
		double PMFatDzero=0.0;
		if(Elt->GetDoubleAttribute("PMFatDzero",&PMFatDzero)==EXIT_SUCCESS)
		{
			SetPMFatDzero(PMF,Mask,i,PMFatDzero,MaskPMFInterpol,0);
		}
		double RemovePMFmoreThen=0.0;
		if(Elt->GetDoubleAttribute("RemoveLargePMF",&RemovePMFmoreThen)==EXIT_SUCCESS)
		{
			RemoveLargePMF(PMF,Mask,i,RemovePMFmoreThen);
		}
	}
	
	if(ConvertPoissonPotToPMF)
		CmdConvertPoissonPotToPMF(PMF,Mask,MaskPMFInterpol,true);
	
	
	
	for(i=0;i<Nions;i++)
	{
		AverageThrLaplas1(PMF[i],Mask[i],StreachIter,2);
		
		AverageThrLaplas(PMF[i],Mask[i],SmoothIter,-1);
	}
	//CmdPMFProcess_Interpolate0(PMF,Mask);
	//CmdPMFProcess_Interpolate1(PMF,Mask);
	
	ShowProperties();
	
	//Write PMF
	m_ContWorld->PMF=new float*[Nions];
	for(i=0;i<Nions;i++)
	{
		m_ContWorld->PMF[i]=new float[GS_XYZ];
		for(j=0;j<m_ContWorld->GS_XYZ;j++)
			m_ContWorld->PMF[i][j] = PMF[i][j];
	}
	
	VectorIntField3D vInt(GridSize, float(GridScale), Nions,Mask);
	vInt.WriteToFile("Mask.gz");
	return EXIT_SUCCESS;
}
int PMFProcessing::CopyPMFtoContWorld()
{
	int i,j;
	m_ContWorld->PMF=new float*[Nions];
	for(i=0;i<Nions;i++)
	{
		m_ContWorld->PMF[i]=new float[GS_XYZ];
		for(j=0;j<m_ContWorld->GS_XYZ;j++)
			m_ContWorld->PMF[i][j] = PMF[i][j];
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::CmdConvertPoissonPotToPMF(double** V,int **_Mask,int MaskValueToSet,bool onlyWhereDnonZero)
{
	int i,j;
	pnpPrint0("PMFProcessing::CmdConvertPoissonPotToPMF\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld->Potential,"Potential is not initiate\n");
		
	for(i=0;i<Nions;i++)
	{
		DbgPrint0("IonsQ[%d]=%f\n",i,m_ContWorld->IonsQ[i]);
		if(onlyWhereDnonZero)
		{
			for(j=0;j<m_ContWorld->GS_XYZ;j++)
			{
				if(_Mask[i][j]==0&&m_ContWorld->NIndexing->GetDiffFloat(i,j)!=0.0)
				{
					V[i][j]=m_ContWorld->Potential[j]*m_ContWorld->IonsQ[i];
					_Mask[i][j]=MaskValueToSet;
				}
			}
		}
		else
		{
			for(j=0;j<m_ContWorld->GS_XYZ;j++)
			{
				if(_Mask[i][j]==0)
				{
					V[i][j]=m_ContWorld->Potential[j]*m_ContWorld->IonsQ[i];
					_Mask[i][j]=MaskValueToSet;
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::SetEionPolynom4(double** dPMF,int** Mask,double a4,double a3,double a2,double a1,double a0,double Rlim)
{
	pnpPrint0("PMFProcessing::SetEionPolynom4\n");
	pnpPrint0("\ta4=%e\n",a4);
	pnpPrint0("\ta3=%e\n",a3);
	pnpPrint0("\ta2=%e\n",a2);
	pnpPrint0("\ta1=%e\n",a1);
	pnpPrint0("\ta0=%e\n",a0);
	pnpPrint0("\tRlim=%e\n",Rlim);
	int i,j,k,pnt,ion;
	double x,y,z,r;
	
	int irc[3];
	float frc[3];
	for(j=0;j<3;j++)
	{
		irc[j]=m_ContWorld->GridSizeGlobal[j]/2;
		frc[j]=irc[j];
	}
	for(ion=0;ion<Nions;ion++)
		for(i=0;i<GS_X;i++)
			for(j=0;j<GS_Y;j++)
				for(k=0;k<GS_Z;k++)
	{
		pnt=i+j*GS_X+k*GS_XY;
		x=i-irc[0];
		y=j-irc[1];
		z=k-irc[2];
		r=sqrt(x*x+y*y+z*z);
		if(r<=Rlim)
		{
			dPMF[ion][pnt]=a4*r*r*r*r+a3*r*r*r+a2*r*r+a1*r+a0;
			Mask[ion][pnt]=MaskPMFCalc;
		}
	}
	
	return EXIT_SUCCESS;
}
int PMFProcessing::SetEionPolynom6(double** dPMF,int** Mask,double a6,double a5,double a4,double a3,double a2,double a1,double a0,double Rlim)
{
	pnpPrint0("PMFProcessing::SetEionPolynom6\n");
	pnpPrint0("\ta6=%e\n",a6);
	pnpPrint0("\ta5=%e\n",a5);
	pnpPrint0("\ta4=%e\n",a4);
	pnpPrint0("\ta3=%e\n",a3);
	pnpPrint0("\ta2=%e\n",a2);
	pnpPrint0("\ta1=%e\n",a1);
	pnpPrint0("\ta0=%e\n",a0);
	pnpPrint0("\tRlim=%e\n",Rlim);
	int i,j,k,pnt,ion;
	double x,y,z,r;
	
	int irc[3];
	float frc[3];
	for(j=0;j<3;j++)
	{
		irc[j]=m_ContWorld->GridSizeGlobal[j]/2;
		frc[j]=irc[j];
	}
	for(ion=0;ion<Nions;ion++)
		for(i=0;i<GS_X;i++)
			for(j=0;j<GS_Y;j++)
				for(k=0;k<GS_Z;k++)
	{
		pnt=i+j*GS_X+k*GS_XY;
		x=i-irc[0];
		y=j-irc[1];
		z=k-irc[2];
		r=sqrt(x*x+y*y+z*z);
		if(r<=Rlim)
		{
			dPMF[ion][pnt]=a6*r*r*r*r*r*r+a5*r*r*r*r*r+a4*r*r*r*r+a3*r*r*r+a2*r*r+a1*r+a0;
			Mask[ion][pnt]=MaskPMFCalc;
		}
	}
	
	return EXIT_SUCCESS;
}
int PMFProcessing::SetEionPolynom8(double** dPMF,int** Mask,double a8,double a7,double a6,double a5,double a4,double a3,double a2,double a1,double a0,double Rlim)
{
	pnpPrint0("PMFProcessing::SetEionPolynom8\n");
	pnpPrint0("\ta8=%e\n",a8);
	pnpPrint0("\ta7=%e\n",a7);
	pnpPrint0("\ta6=%e\n",a6);
	pnpPrint0("\ta5=%e\n",a5);
	pnpPrint0("\ta4=%e\n",a4);
	pnpPrint0("\ta3=%e\n",a3);
	pnpPrint0("\ta2=%e\n",a2);
	pnpPrint0("\ta1=%e\n",a1);
	pnpPrint0("\ta0=%e\n",a0);
	pnpPrint0("\tRlim=%e\n",Rlim);
	int i,j,k,pnt,ion;
	double x,y,z,r;
	
	int irc[3];
	float frc[3];
	for(j=0;j<3;j++)
	{
		irc[j]=m_ContWorld->GridSizeGlobal[j]/2;
		frc[j]=irc[j];
	}
	for(ion=0;ion<Nions;ion++)
		for(i=0;i<GS_X;i++)
			for(j=0;j<GS_Y;j++)
				for(k=0;k<GS_Z;k++)
	{
		pnt=i+j*GS_X+k*GS_XY;
		x=i-irc[0];
		y=j-irc[1];
		z=k-irc[2];
		r=sqrt(x*x+y*y+z*z);
		if(r<=Rlim)
		{
			dPMF[ion][pnt]=a8*r*r*r*r*r*r*r*r+a7*r*r*r*r*r*r+a6*r*r*r*r*r*r+a5*r*r*r*r*r+a4*r*r*r*r+a3*r*r*r+a2*r*r+a1*r+a0;
			Mask[ion][pnt]=MaskPMFCalc;
		}
	}
	
	return EXIT_SUCCESS;
}
int PMFProcessing::ReadEnergyFrom9ColumnFileCopyOctants(const char * FileName,double** dPMF,int** Mask,float ConvFac)
{
	//FileProtIon
	pnpPrint0("FileName %s dPMF = %p Mask = %p\n", FileName,dPMF,Mask);
	int i,j;
	int ion,pnt;
	double dtmp;
	float drel,dconvfac;
	char str[PNP_MAP_IO_STRING_LENGTH];
	
	int irc[3];
	float frc[3];
	for(j=0;j<3;j++)
	{
		irc[j]=m_ContWorld->GridSizeGlobal[j]/2;
		frc[j]=irc[j];
	}
	
	FILE* in=NULL;
	in=fopen(FileName,"r");
	if(in!=NULL)
	{
		int* count=new int[Nions];
		int* countHiConv=new int[Nions];
		double* Eavr=new double[Nions];
		double* Esd=new double[Nions];
		double* Emin=new double[Nions];
		double* Emax=new double[Nions];
		
		pnpPrint0("Open %s for reading energies\n",FileName);
		for(i=0;i<Nions;i++)
		{
			count[i]=0;
			countHiConv[i]=0;
			Eavr[i]=0.0;
			Esd[i]=0.0;
		}
		while(fgets(str,PNP_MAP_IO_STRING_LENGTH,in)!=NULL)
		{
			int iPos;
			float rabs[3];
			float rabs0[3];
			if(sscanf(str,"%d %d %d %g %g %g %le %e %e", &ion, &iPos, &pnt, &rabs[0], &rabs[1], &rabs[2], &dtmp, &drel, &dconvfac)==9)
			{
				if(dconvfac<=ConvFac)
				{
					pnt=m_ContWorld->ConvGlobalExternalXYZToGrdPnt(rabs[0], rabs[1], rabs[2]);
					dPMF[ion][pnt]=dtmp;
					Mask[ion][pnt]=MaskPMFCalc;
					
					if(count[ion]==0)
					{
						Emin[ion]=dtmp;
						Emax[ion]=dtmp;
					}
					
					count[ion]++;
					Eavr[ion]+=dtmp;
					Esd[ion]+=dtmp*dtmp;
					if(dtmp>Emax[ion])Emax[ion]=dtmp;
					if(dtmp<Emin[ion])Emin[ion]=dtmp;
					
					//other octants:
					//2
					pnt=m_ContWorld->ConvGlobalExternalXYZToGrdPnt(-rabs[0], rabs[1], rabs[2]);
					dPMF[ion][pnt]=dtmp;
					Mask[ion][pnt]=MaskPMFCalc;
					//3
					pnt=m_ContWorld->ConvGlobalExternalXYZToGrdPnt(rabs[0], -rabs[1], rabs[2]);
					dPMF[ion][pnt]=dtmp;
					Mask[ion][pnt]=MaskPMFCalc;
					//4
					pnt=m_ContWorld->ConvGlobalExternalXYZToGrdPnt(-rabs[0], -rabs[1], rabs[2]);
					dPMF[ion][pnt]=dtmp;
					Mask[ion][pnt]=MaskPMFCalc;
					//1
					pnt=m_ContWorld->ConvGlobalExternalXYZToGrdPnt(rabs[0], rabs[1], -rabs[2]);
					dPMF[ion][pnt]=dtmp;
					Mask[ion][pnt]=MaskPMFCalc;
					//2
					pnt=m_ContWorld->ConvGlobalExternalXYZToGrdPnt(-rabs[0], rabs[1], -rabs[2]);
					dPMF[ion][pnt]=dtmp;
					Mask[ion][pnt]=MaskPMFCalc;
					//3
					pnt=m_ContWorld->ConvGlobalExternalXYZToGrdPnt(rabs[0], -rabs[1], -rabs[2]);
					dPMF[ion][pnt]=dtmp;
					Mask[ion][pnt]=MaskPMFCalc;
					//4
					pnt=m_ContWorld->ConvGlobalExternalXYZToGrdPnt(-rabs[0], -rabs[1], -rabs[2]);
					dPMF[ion][pnt]=dtmp;
					Mask[ion][pnt]=MaskPMFCalc;
					
				}
				else
				{
					countHiConv[ion]++;
//					if(PrintDisconvPoint)pnpPrint0("ion %d pnt %d\n",ion,pnt);
				}
			}
			else
			{
				pnpPrint0("WARNING:Unknown line: %s :WARNING\n",str);
			}
		}
		fclose(in);
		j=0;
		for(ion=0;ion<Nions;ion++)
		{
			j+=count[ion];
			Eavr[ion]=Eavr[ion]/count[ion];
			Esd[ion]=sqrt((Esd[ion]/count[ion])-Eavr[ion]*Eavr[ion]);
			pnpPrint0("Energies for Ion [%d]\n",ion);
			pnpPrint0("\t%d points was added.\n",count[ion]);
			//if(IonsHaveSameSize&&countHiConv[ion]>0)
			pnpPrint0("\t%d was not added due to high convergens factor(> %g)\n",countHiConv[ion],ConvFac);
			pnpPrint0("\tEavr=%.12g Esd=%.12g Emin=%.12g Emax=%.12g\n"
					, Eavr[ion], Esd[ion], Emin[ion], Emax[ion]);
			pnpPrint0("\tEavr-Eprot=%.12g Emin-Eprot=%.12g Emax-Eprot=%.12g\n"
					, Eavr[ion]-Eprot, Emin[ion]-Eprot, Emax[ion]-Eprot);
		}
		pnpPrint0("\t%d total reads\n",j);
		DeleteCArray(count);
		DeleteCArray(countHiConv);
		DeleteCArray(Eavr);
		DeleteCArray(Esd);
		DeleteCArray(Emin);
		DeleteCArray(Emax);
	}
	else
	{
		pnpPrint0("WARNING: Can't open %s to read\n",FileName);
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::ReadEnergyFrom9ColumnFile(const char * FileName,double** dPMF,int** Mask,float ConvFac)
{
	//FileProtIon
	pnpPrint0("FileName %s dPMF = %p Mask = %p\n", FileName,dPMF,Mask);
	int i,j;
	int ion,pnt;
	double dtmp;
	float drel,dconvfac;
	char str[PNP_MAP_IO_STRING_LENGTH];
	
	FILE* in=NULL;
	in=fopen(FileName,"r");
	if(in!=NULL)
	{
		int* count=new int[Nions];
		int* countHiConv=new int[Nions];
		double* Eavr=new double[Nions];
		double* Esd=new double[Nions];
		double* Emin=new double[Nions];
		double* Emax=new double[Nions];
		
		pnpPrint0("Open %s for reading energies\n",FileName);
		for(i=0;i<Nions;i++)
		{
			count[i]=0;
			countHiConv[i]=0;
			Eavr[i]=0.0;
			Esd[i]=0.0;
		}
		while(fgets(str,PNP_MAP_IO_STRING_LENGTH,in)!=NULL)
		{
			int iPos;
			float rabs[3];
			if(sscanf(str,"%d %d %d %g %g %g %le %e %e", &ion, &iPos, &pnt, &rabs[0], &rabs[1], &rabs[2], &dtmp, &drel, &dconvfac)==9)
			{
				if(dconvfac<=ConvFac)
				{
					pnt=m_ContWorld->ConvGlobalExternalXYZToGrdPnt(rabs[0], rabs[1], rabs[2]);
					dPMF[ion][pnt]=dtmp;
					Mask[ion][pnt]=MaskPMFCalc;
					
					if(count[ion]==0)
					{
						Emin[ion]=dtmp;
						Emax[ion]=dtmp;
					}
					
					count[ion]++;
					Eavr[ion]+=dtmp;
					Esd[ion]+=dtmp*dtmp;
					if(dtmp>Emax[ion])Emax[ion]=dtmp;
					if(dtmp<Emin[ion])Emin[ion]=dtmp;
				}
				else
				{
					countHiConv[ion]++;
//					if(PrintDisconvPoint)pnpPrint0("ion %d pnt %d\n",ion,pnt);
				}
			}
			else
			{
				pnpPrint0("WARNING:Unknown line: %s :WARNING\n",str);
			}
		}
		fclose(in);
		j=0;
		for(ion=0;ion<Nions;ion++)
		{
			j+=count[ion];
			Eavr[ion]=Eavr[ion]/count[ion];
			Esd[ion]=sqrt((Esd[ion]/count[ion])-Eavr[ion]*Eavr[ion]);
			pnpPrint0("Energies for Ion [%d]\n",ion);
			pnpPrint0("\t%d points was added.\n",count[ion]);
			//if(IonsHaveSameSize&&countHiConv[ion]>0)
			pnpPrint0("\t%d was not added due to high convergens factor(> %g)\n",countHiConv[ion],ConvFac);
			pnpPrint0("\tEavr=%.12g Esd=%.12g Emin=%.12g Emax=%.12g\n"
					, Eavr[ion], Esd[ion], Emin[ion], Emax[ion]);
			pnpPrint0("\tEavr-Eprot=%.12g Emin-Eprot=%.12g Emax-Eprot=%.12g\n"
					, Eavr[ion]-Eprot, Emin[ion]-Eprot, Emax[ion]-Eprot);
		}
		pnpPrint0("\t%d total reads\n",j);
		DeleteCArray(count);
		DeleteCArray(countHiConv);
		DeleteCArray(Eavr);
		DeleteCArray(Esd);
		DeleteCArray(Emin);
		DeleteCArray(Emax);
	}
	else
	{
		pnpPrint0("WARNING: Can't open %s to read\n",FileName);
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::ReadEnergyFrom9ColumnFileSubEion(const char * FileName,double** dPMF,int** Mask,float ConvFac,double** Eion)
{
	//FileProtIon
	pnpPrint0("FileName %s dPMF = %p Mask = %p\n", FileName,dPMF,Mask);
	int i,j;
	int ion,pnt;
	double dtmp;
	float drel,dconvfac;
	char str[PNP_MAP_IO_STRING_LENGTH];
	
	FILE* in=NULL;
	in=fopen(FileName,"r");
	if(in!=NULL)
	{
		int* count=new int[Nions];
		int* countHiConv=new int[Nions];
		double* Eavr=new double[Nions];
		double* Esd=new double[Nions];
		double* Emin=new double[Nions];
		double* Emax=new double[Nions];
		
		pnpPrint0("Open %s for reading energies\n",FileName);
		for(i=0;i<Nions;i++)
		{
			count[i]=0;
			countHiConv[i]=0;
			Eavr[i]=0.0;
			Esd[i]=0.0;
		}
		while(fgets(str,PNP_MAP_IO_STRING_LENGTH,in)!=NULL)
		{
			int iPos;
			float rabs[3];
			if(sscanf(str,"%d %d %d %g %g %g %le %e %e", &ion, &iPos, &pnt, &rabs[0], &rabs[1], &rabs[2], &dtmp, &drel, &dconvfac)==9)
			{
				if(dconvfac<=ConvFac)
				{
					pnt=m_ContWorld->ConvGlobalExternalXYZToGrdPnt(rabs[0], rabs[1], rabs[2]);
					dPMF[ion][pnt]=dtmp-Eion[ion][pnt];
					Mask[ion][pnt]=MaskPMFCalc;
					
					if(count[ion]==0)
					{
						Emin[ion]=dtmp;
						Emax[ion]=dtmp;
					}
					
					count[ion]++;
					Eavr[ion]+=dtmp;
					Esd[ion]+=dtmp*dtmp;
					if(dtmp>Emax[ion])Emax[ion]=dtmp;
					if(dtmp<Emin[ion])Emin[ion]=dtmp;
				}
				else
				{
					countHiConv[ion]++;
//					if(PrintDisconvPoint)pnpPrint0("ion %d pnt %d\n",ion,pnt);
				}
			}
			else
			{
				pnpPrint0("WARNING:Unknown line: %s :WARNING\n",str);
			}
		}
		fclose(in);
		j=0;
		for(ion=0;ion<Nions;ion++)
		{
			j+=count[ion];
			Eavr[ion]=Eavr[ion]/count[ion];
			Esd[ion]=sqrt((Esd[ion]/count[ion])-Eavr[ion]*Eavr[ion]);
			pnpPrint0("Energies for Ion [%d]\n",ion);
			pnpPrint0("\t%d points was added.\n",count[ion]);
			//if(IonsHaveSameSize&&countHiConv[ion]>0)
			pnpPrint0("\t%d was not added due to high convergens factor(> %g)\n",countHiConv[ion],ConvFac);
			pnpPrint0("\tEavr=%.12g Esd=%.12g Emin=%.12g Emax=%.12g\n"
					, Eavr[ion], Esd[ion], Emin[ion], Emax[ion]);
			pnpPrint0("\tEavr-Eprot=%.12g Emin-Eprot=%.12g Emax-Eprot=%.12g\n"
					, Eavr[ion]-Eprot, Emin[ion]-Eprot, Emax[ion]-Eprot);
		}
		pnpPrint0("\t%d total reads\n",j);
		DeleteCArray(count);
		DeleteCArray(countHiConv);
		DeleteCArray(Eavr);
		DeleteCArray(Esd);
		DeleteCArray(Emin);
		DeleteCArray(Emax);
	}
	else
	{
		pnpPrint0("WARNING: Can't open %s to read\n",FileName);
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::ReadEnergyFrom5ColumnFile(const char * FileName,double** dPMF,int** Mask,float ConvFac)
{
	//FileProtIon
	pnpPrint0("FileName %s dPMF = %p Mask = %p\n", FileName,dPMF,Mask);
	int i,j;
	int ion,pnt;
	double dtmp;
	float drel,dconvfac;
	char str[PNP_MAP_IO_STRING_LENGTH];
	
	FILE* in=NULL;
	in=fopen(FileName,"r");
	if(in!=NULL)
	{
		int* count=new int[Nions];
		int* countHiConv=new int[Nions];
		double* Eavr=new double[Nions];
		double* Esd=new double[Nions];
		double* Emin=new double[Nions];
		double* Emax=new double[Nions];
		
		pnpPrint0("Open %s for reading energies\n",FileName);
		for(i=0;i<Nions;i++)
		{
			count[i]=0;
			countHiConv[i]=0;
			Eavr[i]=0.0;
			Esd[i]=0.0;
		}
		while(fgets(str,PNP_MAP_IO_STRING_LENGTH,in)!=NULL)
		{
			if(sscanf(str,"%d %d %le %g %g",&ion,&pnt,&dtmp,&drel,&dconvfac)==5)
			{
				if(dconvfac<=ConvFac)
				{
					dPMF[ion][pnt]=dtmp;
					Mask[ion][pnt]=MaskPMFCalc;
					
					if(count[ion]==0)
					{
						Emin[ion]=dtmp;
						Emax[ion]=dtmp;
					}
					
					count[ion]++;
					Eavr[ion]+=dtmp;
					Esd[ion]+=dtmp*dtmp;
					if(dtmp>Emax[ion])Emax[ion]=dtmp;
					if(dtmp<Emin[ion])Emin[ion]=dtmp;
				}
				else
				{
					countHiConv[ion]++;
//					if(PrintDisconvPoint)pnpPrint0("ion %d pnt %d\n",ion,pnt);
				}
			}
			else
			{
				pnpPrint0("WARNING:Unknown line: %s :WARNING\n",str);
			}
		}
		fclose(in);
		j=0;
		for(ion=0;ion<Nions;ion++)
		{
			j+=count[ion];
			Eavr[ion]=Eavr[ion]/count[ion];
			Esd[ion]=sqrt((Esd[ion]/count[ion])-Eavr[ion]*Eavr[ion]);
			pnpPrint0("Energies for Ion [%d]\n",ion);
			pnpPrint0("\t%d points was added.\n",count[ion]);
			//if(IonsHaveSameSize&&countHiConv[ion]>0)
				pnpPrint0("\t%d was not added due to high convergens factor(> %g)\n",countHiConv[ion],ConvFac);
			pnpPrint0("\tEavr=%.12g Esd=%.12g Emin=%.12g Emax=%.12g\n"
					, Eavr[ion], Esd[ion], Emin[ion], Emax[ion]);
			pnpPrint0("\tEavr-Eprot=%.12g Emin-Eprot=%.12g Emax-Eprot=%.12g\n"
					, Eavr[ion]-Eprot, Emin[ion]-Eprot, Emax[ion]-Eprot);
		}
		pnpPrint0("\t%d total reads\n",j);
		DeleteCArray(count);
		DeleteCArray(countHiConv);
		DeleteCArray(Eavr);
		DeleteCArray(Esd);
		DeleteCArray(Emin);
		DeleteCArray(Emax);
	}
	else
	{
		pnpPrint0("WARNING: Can't open %s to read\n",FileName);
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
/*
{
	in=fopen(FileIon,"r");
	if(in!=NULL)
	{
		int* count=new int[Nions];
		int* countHiConv=new int[Nions];
		double* Eavr=new double[Nions];
		double* Esd=new double[Nions];
		double* Emin=new double[Nions];
		double* Emax=new double[Nions];
		
		pnpPrint0("Open %s for ion energy in solvent\n",FileIon);
		for(i=0;i<Nions;i++)
		{
			count[i]=0;
			countHiConv[i]=0;
			Eavr[i]=0.0;
			Esd[i]=0.0;
		}
		while(fgets(str,100,in)!=NULL)
		{
			if(sscanf(str,"%d %d %le %g %g",&ion,&pnt,&dtmp,&drel,&dconvfac)==5)
			{
				if((IonsHaveSameSize==true&&dconvfac<=ConvFac)||(IonsHaveSameSize==false))
				{
					dIon[ion][pnt]=dtmp;
					MaskIon[ion][pnt]=1;
					
					if(count[ion]==0)
					{
						Emin[ion]=dtmp;
						Emax[ion]=dtmp;
					}
					
					count[ion]++;
					Eavr[ion]+=dtmp;
					Esd[ion]+=dtmp*dtmp;
					if(dtmp>Emax[ion])Emax[ion]=dtmp;
					if(dtmp<Emin[ion])Emin[ion]=dtmp;
				}
			}
			if(IonsHaveSameSize==true&&dconvfac>=ConvFac)
			{
				countHiConv[ion]++;
			}
		}
		fclose(in);
		for(ion=0;ion<Nions;ion++)
		{
			Eavr[ion]=Eavr[ion]/count[ion];
			Esd[ion]=sqrt((Esd[ion]/count[ion])-Eavr[ion]*Eavr[ion]);
			
			Eion[ion]=Eavr[ion];
			
			pnpPrint0("Ion[%d]\n",ion);
			pnpPrint0("\t%d points was added\n",count[ion]);
			//if(IonsHaveSameSize&&countHiConv[ion]>0)
				pnpPrint0("\t%d was not added due to high convergens factor(> %g)\n",countHiConv[ion],ConvFac);
			pnpPrint0("\tEavr=%.12lg Esd=%.12lg Emin=%.12lg Emax=%.12lg\n"
				, Eavr[ion], Esd[ion], Emin[ion], Emax[ion]);
		}
		DeleteCArray(count);
		DeleteCArray(countHiConv);
		DeleteCArray(Eavr);
		DeleteCArray(Esd);
		DeleteCArray(Emin);
		DeleteCArray(Emax);
	}
	else
	{
		pnpPrint0("WARNING: Can't open %s (FileIon)\n",FileIon);
	}
	pnpPrint0("Eprot=%lf\n",Eprot);
	for(ion=0;ion<Nions;ion++)
	{
		pnpPrint0("Eion[%d]=%lf\n",ion,Eion[ion]);
	}*/
	//
int PMFProcessing::CalculateAvrEion(bool IonsHaveSameSize)
{
	int ion,j,_Nions;
	
	
	if(IonsHaveSameSize)_Nions=1;
	else _Nions=Nions;
	//Combine dPMF and dIon
	for(ion=0;ion<_Nions;ion++)
	{
		int count=0;
		double tmp=0.0;
		for(j=0;j<m_ContWorld->GS_XYZ;j++)
		{
			if(MaskIon[ion][j]==1)
			{
				tmp+=Eion[ion][j];
				count++;
			}
		}
		EavrIon[ion]=tmp/double(count);
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::CombineEprotionEprotEion(bool IonsHaveSameSize)
{
	int ion,j;
	//Combine dPMF and dIon
	for(ion=0;ion<Nions;ion++)
	{
		int count0=0,count1=0;
		for(j=0;j<m_ContWorld->GS_XYZ;j++)
		{
			if(Mask[ion][j]==1)
			{
				if(IonsHaveSameSize)
				{
					if(MaskIon[0][j]==MaskPMFCalc)
					{
						PMF[ion][j]=PMF[ion][j]-Eion[0][j]-Eprot;
						count1++;
					}
					else
					{
						PMF[ion][j]=PMF[ion][j]-EavrIon[0]-Eprot;
						count0++;
					}
				}
				else
				{
					if(MaskIon[ion][j]==MaskPMFCalc)
					{
						PMF[ion][j]=PMF[ion][j]-Eion[ion][j]-Eprot;
						count1++;
					}
					else
					{
						PMF[ion][j]=PMF[ion][j]-EavrIon[ion]-Eprot;
						count0++;
					}
				}
			}
		}
		pnpPrint0("Substruct energy of ion[%d]. Using one energy %d, using calculated at grid %d \n",ion,count0,count1);
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::CmdPMFProcess_Interpolate1(double** V,int **Mask)
{
	//
	int i,j,ix,iy,iz;
	for(i=0;i<Nions;i++)
		for(ix=1;ix<GS_X-1;ix++)
			for(iy=1;iy<GS_Y-1;iy++)
				for(iz=1;iz<GS_Z-1;iz++)
	{
		j=ix+iy*GS_X+iz*GS_XY;
		if(Mask[i][j]==0)
		{
			if( (Mask[i][j-1]!=0)&&
							(Mask[i][j+1]!=0)&&
							(Mask[i][j-m_ContWorld->GS_X]!=0)&&
							(Mask[i][j+m_ContWorld->GS_X]!=0)&&
							(Mask[i][j-m_ContWorld->GS_XY]!=0)&&
							(Mask[i][j+m_ContWorld->GS_XY]!=0))
			{
				PMF[i][j] = (PMF[i][j+1]+PMF[i][j-1]+PMF[i][j+GS_X]+PMF[i][j-GS_X]+PMF[i][j+GS_XY]+PMF[i][j-GS_XY])/6.0;
				Mask[i][j]=3;
			}
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::ShowProperties()
{
	int i, j;
	//Show properties

	for(i=0;i<Nions;i++)
	{
		int count0=0;
		int count1=0;
		double Eavr0=0.0,Esd0=0.0,Emin0,Emax0;
		double Eavr1=0.0,Esd1=0.0,Emin1,Emax1;
		for(j=0;j<m_ContWorld->GS_XYZ;j++)
		{
			if(Mask[i][j]==1)
			{
				if(count1==0)
				{
					Emin1=PMF[i][j];
					Emax1=PMF[i][j];
				}
				
				count1++;
				
				Eavr1+=PMF[i][j];
				Esd1+=PMF[i][j]*PMF[i][j];
				if(PMF[i][j]>Emax1)Emax1=PMF[i][j];
				if(PMF[i][j]<Emin1)Emin1=PMF[i][j];
			}
			
			{
				if(count0==0)
				{
					Emin0=PMF[i][j];
					Emax0=PMF[i][j];
				}
				
				count0++;
				
				Eavr0+=PMF[i][j];
				Esd0+=PMF[i][j]*PMF[i][j];
				if(PMF[i][j]>Emax0)Emax0=PMF[i][j];
				if(PMF[i][j]<Emin0)Emin0=PMF[i][j];
			}
		}
		Eavr0=Eavr0/count0;
		Esd0=sqrt((Esd0/count0)-Eavr0*Eavr0);
		
		Eavr1=Eavr1/count1;
		Esd1=sqrt((Esd1/count1)-Eavr1*Eavr1);
		
		pnpPrint0("Properties of finile PMFs ion=%d\n",i);
		pnpPrint0("\tCalculated PMF\n");
		
		pnpPrint0("\t\t%d points from calculations was added\n", count1);
		pnpPrint0("\t\tEavr1=%.12lg Esd1=%.12lg Emin1=%.12lg Emax1=%.12lg\n", 
			Eavr1, Esd1, Emin1, Emax1);
		
		pnpPrint0("\tOver all grids of the system\n");
		pnpPrint0("\t\t%d total grids in system\n", count0);
		pnpPrint0("\t\tEavr0=%.12lg Esd0=%.12lg Emin0=%.12lg Emax0=%.12lg\n", 
			Eavr0, Esd0, Emin0, Emax0);
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::CmdPMFProcessOld(const TiXmlElement *Elt)
{
	DefClock0;
	fprintf(stdout,"<ResultsPMFProcess>\n");
	StartClock0;
	int i,j;
	char FileIon[128]="\0";
	char FileProtIon[128]="\0";
	char FilePMFout[128]="\0";
	double **dPMF;
	int Nions,*Npnts;
	int **Mask;
	int PMFpnts;
	double Eprot;
	Elt->GetCStrAttribute("FileProtIon",FileProtIon);
	Elt->GetCStrAttribute("FileIon",FileIon);
	Elt->GetCStrAttribute("PMFout",FilePMFout);
	Elt->GetIntAttribute("Nions",&Nions);
	Elt->GetDoubleAttribute("Eprot",&Eprot);
	Npnts=new int[Nions];
	dPMF=new double*[Nions];
	Mask=new int*[Nions];
	for(i=0;i<Nions;i++)
	{
		dPMF[i]=new double[m_ContWorld->GS_XYZ];
		Mask[i]=new int[m_ContWorld->GS_XYZ];
		for(j=0;j<m_ContWorld->GS_XYZ;j++)
		{
			dPMF[i][j]=0.0;
			Mask[i][j]=0;
			
		}
	}
	
	int ion,pnt;
	double dtmp;
	char str[256];
	
	FILE* in=NULL;
	in=fopen(FileProtIon,"r");
	if(in!=NULL)
	{
		while(fgets(str,100,in)!=NULL)
		{
			if(sscanf(str,"%d %d %le",&ion,&pnt,&dtmp)==3)
			{
				dPMF[ion][pnt]=dtmp-Eprot;
				Mask[ion][pnt]=1;
			}
		}
		fclose(in);
	}
	else
		pnpPrint0("WARNING: Can't open %s (FileProtIon)\n",FileProtIon);
	
	in=fopen(FileIon,"r");
	if(in!=NULL)
	{
		while(fgets(str,100,in)!=NULL)
		{
			if(sscanf(str,"%d %d %le",&ion,&pnt,&dtmp)==3)
			{
				dPMF[ion][pnt]=dPMF[ion][pnt]-dtmp;
			}
		}
		fclose(in);
	}
	else
		pnpPrint0("WARNING: Can't open %s (FileIon)\n",FileIon);
	
	VectorField3D Vec(m_ContWorld->GridSize, m_ContWorld->GridScale, Nions);
	for(i=0;i<Nions;i++)
	{
		for(j=0;j<m_ContWorld->GS_XYZ;j++)
			Vec.V[i][j] = dPMF[i][j];
	}
	
	int ix,iy,iz;
	for(i=0;i<Nions;i++)
		for(ix=1;ix<m_ContWorld->GS_X-1;ix++)
			for(iy=1;iy<m_ContWorld->GS_Y-1;iy++)
				for(iz=1;iz<m_ContWorld->GS_Z-1;iz++)
	{
		j=ix+iy*m_ContWorld->GS_X+iz*m_ContWorld->GS_XY;
		if(Mask[i][j]==0&&m_ContWorld->NIndexing->GetDiffFloat(i,j)!=0.0)
		{
			//DbgPrint0("%d %d %d %d %d %d\n",Mask[i][j-1],Mask[i][j+1],Mask[i][j-m_ContWorld->GS_X],Mask[i][j+m_ContWorld->GS_X],Mask[i][j-m_ContWorld->GS_XY],Mask[i][j+m_ContWorld->GS_XY]);
			
			if( (Mask[i][j-1]!=0)&&
									(Mask[i][j+1]!=0)&&
									(Mask[i][j-m_ContWorld->GS_X]!=0)&&
									(Mask[i][j+m_ContWorld->GS_X]!=0)&&
									(Mask[i][j-m_ContWorld->GS_XY]!=0)&&
									(Mask[i][j+m_ContWorld->GS_XY]!=0))
			{
				Vec.V[i][j] = (Vec.V[i][j+1]+Vec.V[i][j-1]+Vec.V[i][j+m_ContWorld->GS_X]+Vec.V[i][j-m_ContWorld->GS_X]+Vec.V[i][j+m_ContWorld->GS_XY]+Vec.V[i][j-m_ContWorld->GS_XY])/6.0;
				Mask[i][j]=3;
			}
		}
	}
	//Potential
	if(m_ContWorld->Potential!=NULL)
	{
		int iq[2]={+1.0,-1.0};
		for(i=0;i<Nions;i++)
		{
			for(j=0;j<m_ContWorld->GS_XYZ;j++)
			{
				if(Mask[i][j]==0)
				{
					Vec.V[i][j] = iq[i]*m_ContWorld->Potential[j];
					Mask[i][j]=10;
				}
			}
		}
	}
	//Write PMF
	Vec.WriteToFile(FilePMFout);
	DeleteCArray(Npnts);
	DeleteCVecArray(dPMF,Nions);
	DeleteCVecArray(Mask,Nions);
	StopClockWMes0("PMFProcess");
	fprintf(stdout,"</ResultsPMFPMFProcess>\n");
	
	return EXIT_SUCCESS;
}
bool PMFProcessing::CmdPMFProcess_CheckPnt(int pnt)
{
	if(pnt>=0&&pnt<m_ContWorld->GS_XYZ)
		return true;
	else
		return false;
}
int PMFProcessing::CmdPMFProcess_Interpolate0(double** V,int **_Mask)
{
	int i,j;
	int shft;
	for(i=0;i<Nions;i++)
	{
		for(j=0;j<m_ContWorld->GS_XYZ;j++)
		{
			if(Mask[i][j]==1)
			{
				CmdPMFProcess_Interpolate0_1(-1,i,j,V,_Mask);
				CmdPMFProcess_Interpolate0_1(1,i,j,V,_Mask);
				CmdPMFProcess_Interpolate0_1(-m_ContWorld->GS_X,i,j,V,_Mask);
				CmdPMFProcess_Interpolate0_1(m_ContWorld->GS_X,i,j,V,_Mask);
				CmdPMFProcess_Interpolate0_1(-m_ContWorld->GS_XY,i,j,V,_Mask);
				CmdPMFProcess_Interpolate0_1(m_ContWorld->GS_XY,i,j,V,_Mask);
				
				CmdPMFProcess_Interpolate0_1(-1-m_ContWorld->GS_X-m_ContWorld->GS_XY,i,j,V,_Mask);
				CmdPMFProcess_Interpolate0_1(-1+m_ContWorld->GS_X-m_ContWorld->GS_XY,i,j,V,_Mask);
				CmdPMFProcess_Interpolate0_1(1-m_ContWorld->GS_X-m_ContWorld->GS_XY,i,j,V,_Mask);
				CmdPMFProcess_Interpolate0_1(1+m_ContWorld->GS_X-m_ContWorld->GS_XY,i,j,V,_Mask);
				CmdPMFProcess_Interpolate0_1(-1-m_ContWorld->GS_X+m_ContWorld->GS_XY,i,j,V,_Mask);
				CmdPMFProcess_Interpolate0_1(-1+m_ContWorld->GS_X+m_ContWorld->GS_XY,i,j,V,_Mask);
				CmdPMFProcess_Interpolate0_1(1-m_ContWorld->GS_X+m_ContWorld->GS_XY,i,j,V,_Mask);
				CmdPMFProcess_Interpolate0_1(1+m_ContWorld->GS_X+m_ContWorld->GS_XY,i,j,V,_Mask);
			}
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::CmdPMFProcess_Interpolate0_1(int shft,int i,int j,double** V,int **_Mask)
{
	if(CmdPMFProcess_CheckPnt(j+shft))
		if(m_ContWorld->NIndexing->GetDiffFloat(i,j+shft)!=0.0)
	{
		if(V[i][j+shft] !=0.0)
			V[i][j+shft]=(V[i][j]+V[i][j+shft])*0.5;
		else
			V[i][j+shft]=V[i][j];
		_Mask[i][j+shft]=2;
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::SetNbhoodToCalcNodeAllIons(double** _V,int **_Mask,int Rstep)
{
		int i;
		for(i=0;i<Nions;i++)
		{
			SetNbhoodToCalcNode(_V[i],_Mask[i],Rstep);
		}
	return EXIT_SUCCESS;
}
int PMFProcessing::SetNbhoodToCalcNode(double* V,int *_Mask,int Rstep)
{
	int i,j,k,i1,j1,k1,i2,j2,k2,i3,j3,k3;
	int GridPoint,GridPoint1,GridPoint2,GridPointEnd;
	int count;
	float x,x1,x2;
	for(i=Rstep;i<GS_X-Rstep;i++)
		for(j=Rstep;j<GS_Y-Rstep;j++)
			for(k=Rstep;k<GS_Z-Rstep;k++)
	{
		GridPoint=i+j*GS_X+k*GS_XY;
		if(_Mask[GridPoint]==1)
		{
			for(i1=i-Rstep;i1<=i+Rstep;i1++)
				for(j1=j-Rstep;j1<=j+Rstep;j1++)
					for(k1=k-Rstep;k1<=k+Rstep;k1++)
			{
				GridPoint1=i1+j1*GS_X+k1*GS_XY;
				if(_Mask[GridPoint1]!=MaskPMFCalc)
				{
					if(_Mask[GridPoint]==MaskPMFInterpol)
					{
						V[GridPoint1]=0.5*(V[GridPoint1]-V[GridPoint]);
						_Mask[GridPoint1]=MaskPMFInterpol;
					}
					else
					{
						V[GridPoint1]=V[GridPoint];
						_Mask[GridPoint1]=MaskPMFInterpol;
					}
				}
			}
		}
	}
	
	return EXIT_SUCCESS;
}
int PMFProcessing::SetGrid(int gsX,int gsY,int gsZ,float gridscale,int m_Nion)
{
	GS_X=gsX;
	GS_Y=gsY;
	GS_Z=gsZ;
	GridScale=gridscale;
	Nions=m_Nion;
	GridSize[0]=gsX;
	GridSize[1]=gsY;
	GridSize[2]=gsZ;
	
	GS_XY=GS_X*GS_Y;
	GS_XYZ=GS_X*GS_Y*GS_Z;
	
	return EXIT_SUCCESS;
}
int PMFProcessing::SetNbhoodToCalcNode_VF(VectorField3D *VF_PMF,VectorIntField3D *VF_Mask,int ion,int Rstep)
{
	float *V=VF_PMF->V[ion];
	int *_Mask=VF_Mask->V[ion];
	int i,j,k,i1,j1,k1,i2,j2,k2,i3,j3,k3;
	int GridPoint,GridPoint1,GridPoint2,GridPointEnd;
	int count;
	float x,x1,x2;
	for(i=Rstep;i<GS_X-Rstep;i++)
		for(j=Rstep;j<GS_Y-Rstep;j++)
			for(k=Rstep;k<GS_Z-Rstep;k++)
	{
		GridPoint=i+j*GS_X+k*GS_XY;
		if(_Mask[GridPoint]==1)
		{
			for(i1=i-Rstep;i1<=i+Rstep;i1++)
				for(j1=j-Rstep;j1<=j+Rstep;j1++)
					for(k1=k-Rstep;k1<=k+Rstep;k1++)
			{
				GridPoint1=i1+j1*GS_X+k1*GS_XY;
				if(_Mask[GridPoint1]!=MaskPMFCalc)
				{
					if(_Mask[GridPoint]==MaskPMFInterpol)
					{
						V[GridPoint1]=0.5*(V[GridPoint1]-V[GridPoint]);
						_Mask[GridPoint1]=MaskPMFInterpol;
					}
					else
					{
						V[GridPoint1]=V[GridPoint];
						_Mask[GridPoint1]=MaskPMFInterpol;
					}
				}
			}
		}
	}
	
	return EXIT_SUCCESS;
}
int PMFProcessing::SetMaskToWherePMF_VF(VectorField3D *VF_PMF,float val,VectorIntField3D *VF_Mask,int MaskToSet)
{
	int GridPoint,ion;
	int i,j,k;
	int GS_XYZ=VF_PMF->GridSize[0]*VF_PMF->GridSize[1]*VF_PMF->GridSize[2];
	for(ion=0;ion<VF_Mask->Nelem;ion++)
	{
		for(GridPoint=0;GridPoint<GS_XYZ;GridPoint++)
		{
			if(VF_PMF->V[ion][GridPoint]==val)
			{
				VF_Mask->V[ion][GridPoint]=MaskToSet;
			}
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::SetPMFatInd(double** V,int **_Mask,int ion,double PMFVal, int MaskWhereSet)
{
	DbgPrint0("SetPMFatInd\n");

	int i,j,k;
	int GridPoint;
	for(i=0;i<GS_X;i++)
		for(j=0;j<GS_Y;j++)
			for(k=0;k<GS_Z;k++)
	{
		GridPoint=i+j*GS_X+k*GS_XY;
		if(_Mask[ion][GridPoint]==MaskWhereSet)
		{
			V[ion][GridPoint]=PMFVal;
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::SetPMFatDzero(double** V,int **_Mask,int ion,double PMFatDzero,int MaskValueToSet,int MaskWhereSet)
{
	DbgPrint0("PMFatDzero = %g \n",PMFatDzero);
	if(m_ContWorld==NULL)
	{
		pnpPrint0("Error: m_ContWorld=NULL, initialize ContWorld first\n");
		return EXIT_FAILURE;
	}
	if(m_ContWorld->NIndexing==NULL)
	{
		pnpPrint0("Error: m_ContWorld->NIndexing=NULL, initialize ContWorld complitly\n");
		return EXIT_FAILURE;
	}
	int i,j,k;
	int GridPoint;
	for(i=0;i<GS_X;i++)
		for(j=0;j<GS_Y;j++)
			for(k=0;k<GS_Z;k++)
	{
		GridPoint=i+j*GS_X+k*GS_XY;
		if(_Mask[ion][GridPoint]==MaskWhereSet&&m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint)==0.0)
		{
			V[ion][GridPoint]=PMFatDzero;
			_Mask[ion][GridPoint]=MaskValueToSet;
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::SetPMFatDzero_VF(VectorField3D *VF_PMF,float val,VectorIntField3D *VF_Mask,int ion,float PMFatDzero,int MaskValueToSet,int MaskWhereSet)
{
	float **V=VF_PMF->V;
	int **_Mask=VF_Mask->V;
	DbgPrint0("PMFatDzero = %g \n",PMFatDzero);
	if(m_ContWorld==NULL)
	{
		pnpPrint0("Error: m_ContWorld=NULL, initialize ContWorld first\n");
		return EXIT_FAILURE;
	}
	if(m_ContWorld->NIndexing==NULL)
	{
		pnpPrint0("Error: m_ContWorld->NIndexing=NULL, initialize ContWorld complitly\n");
		return EXIT_FAILURE;
	}
	int i,j,k;
	int GridPoint;
	for(i=0;i<GS_X;i++)
		for(j=0;j<GS_Y;j++)
			for(k=0;k<GS_Z;k++)
	{
		GridPoint=i+j*GS_X+k*GS_XY;
		if(_Mask[ion][GridPoint]==MaskWhereSet&&m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint)==0.0)
		{
			V[ion][GridPoint]=PMFatDzero;
			_Mask[ion][GridPoint]=MaskValueToSet;
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::SetPMFatDnonzero(double** V,int **_Mask,int ion,double PMFatDzero,int MaskValueToSet,int MaskWhereSet)
{
	DbgPrint0("PMFatDzero = %g \n",PMFatDzero);
	if(m_ContWorld==NULL)
	{
		pnpPrint0("Error: m_ContWorld=NULL, initialize ContWorld first\n");
		return EXIT_FAILURE;
	}
	if(m_ContWorld->NIndexing==NULL)
	{
		pnpPrint0("Error: m_ContWorld->NIndexing=NULL, initialize ContWorld complitly\n");
		return EXIT_FAILURE;
	}
	int i,j,k;
	int GridPoint;
	for(i=0;i<GS_X;i++)
		for(j=0;j<GS_Y;j++)
			for(k=0;k<GS_Z;k++)
	{
		GridPoint=i+j*GS_X+k*GS_XY;
		if(_Mask[ion][GridPoint]==MaskWhereSet&&m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint)>0.0)
		{
			V[ion][GridPoint]=PMFatDzero;
			_Mask[ion][GridPoint]=MaskValueToSet;
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::SetSIPwhereDnotZero(double** V,int **_Mask,int ion,double SIP,float Z0,float Z1,int MaskValueToSet,int MaskWhereSet)
{
	DbgPrint0("SetSIPwhereDnotZero %g \n",SIP);
	if(m_ContWorld==NULL)
	{
		pnpPrint0("Error: m_ContWorld=NULL, initialize ContWorld first\n");
		return EXIT_FAILURE;
	}
	if(m_ContWorld->NIndexing==NULL)
	{
		pnpPrint0("Error: m_ContWorld->NIndexing=NULL, initialize ContWorld complitly\n");
		return EXIT_FAILURE;
	}
	int i,j,k;
	int GridPoint;
	i=GS_Z/2;
	Z0+=i;
	Z1+=i;
	DbgPrint0("SetSIPwhereDnotZero Z %g %g \n",Z0,Z1);
	for(i=0;i<GS_X;i++)
		for(j=0;j<GS_Y;j++)
			for(k=0;k<GS_Z;k++)
	{
		if(k>=Z0&&k<=Z1)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(_Mask[ion][GridPoint]==MaskWhereSet&&m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint)>0.0)
			//if(_Mask[ion][GridPoint]==MaskWhereSet)
			{
				V[ion][GridPoint]=SIP;
				_Mask[ion][GridPoint]=MaskValueToSet;
			}
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::SetBorder(int **_Mask,int Mask0,int MaskWith,int MaskValToSet)
{
	int ion;
	int i,j,k,i1,i2,j2,k2,i3,j3,k3;
	int GridPoint,GridPoint1,GridPoint2,GridPointEnd;
	int count,ix,ix1,ix2;
	float x,x1,x2;
	
	int *tmpMap=new int[GS_XYZ];
	
	
	for(ion=0;ion<Nions;ion++)
	{
		for(GridPoint=0;GridPoint<GS_XYZ;GridPoint++)
		tmpMap[GridPoint]=0;
		//set border & remove 2
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(_Mask[ion][GridPoint]==Mask0)
			{
				if(_Mask[ion][GridPoint-1]==MaskWith||_Mask[ion][GridPoint+1]==MaskWith||
							_Mask[ion][GridPoint-GS_X]==MaskWith||_Mask[ion][GridPoint+GS_X]==MaskWith||
							_Mask[ion][GridPoint-GS_XY]==MaskWith||_Mask[ion][GridPoint+GS_XY]==MaskWith)
				{
					tmpMap[GridPoint]=MaskValToSet;
				}
			}
		}
		for(GridPoint=0;GridPoint<GS_XYZ;GridPoint++)
		{
			if(tmpMap[GridPoint]==MaskValToSet)
			{
				_Mask[ion][GridPoint]=MaskValToSet;
			}
		}
	}
	delete [] tmpMap;
	return EXIT_SUCCESS;
}
int PMFProcessing::LineInterpolation0AllIons(double** _V,int **_Mask,int Rstep)
{
	int i;
	for(i=0;i<Nions;i++)
	{
		LineInterpolation0(_V[i],_Mask[i],Rstep);
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::LineInterpolation0(double* V,int *_Mask,int Rstep)
{
	int i,j,k,i1,i2,j2,k2,i3,j3,k3;
	int GridPoint,GridPoint1,GridPoint2,GridPointEnd;
	int count,ix,ix1,ix2;
	float x,x1,x2;
	
	//set border & remove 2
	for(i=1;i<GS_X-1;i++)
		for(j=1;j<GS_Y-1;j++)
			for(k=1;k<GS_Z-1;k++)
	{
		GridPoint=i+j*GS_X+k*GS_XY;
		if(_Mask[GridPoint]>0)
		{
			if(_Mask[GridPoint-1]==0||_Mask[GridPoint+1]==0||
						_Mask[GridPoint-GS_X]==0||_Mask[GridPoint+GS_X]==0||
						_Mask[GridPoint-GS_XY]==0||_Mask[GridPoint+GS_XY]==0)
			{
				_Mask[GridPoint]=MaskRegionBoarder;
			}
		}
	}
	
	for(i=1;i<GS_X-1;i++)
		for(j=1;j<GS_Y-1;j++)
			for(k=1;k<GS_Z-1;k++)
	{
		GridPoint=i+j*GS_X+k*GS_XY;
		if(_Mask[GridPoint]==MaskPMFInterpol)
		{
			//V[GridPoint]=0.0;
			//_Mask[GridPoint]=0;
		}
		else
		{
			_Mask[GridPoint]=-_Mask[GridPoint];
		}
	}
	
	for(i1=0;i1<1;i1++)
	{
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint1=i+j*GS_X+k*GS_XY;
			GridPointEnd=(j+1)*GS_X+k*GS_XY;
			count=0;
			if(_Mask[GridPoint1]<0)
			{
				GridPoint2=GridPoint1+1;
				while(_Mask[GridPoint2]==MaskPMFInterpol&&GridPoint2<GridPointEnd)
				{
					GridPoint2++;
				}
				if(GridPoint2<GridPointEnd&&GridPoint2-GridPoint1<=Rstep)
				{
					x1=(float)GridPoint1;
					x2=(float)GridPoint2;
					for(GridPoint=GridPoint1+1;GridPoint<GridPoint2;GridPoint++)
					{
						x=(float)GridPoint;
						V[GridPoint]=V[GridPoint1]+(V[GridPoint2]-V[GridPoint1])*(x-x1)/(x2-x1);
						//_Mask[GridPoint]=MaskPMFInterpol;
					}
				}
			}
		}
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint1=i+j*GS_X+k*GS_XY;
			GridPointEnd=i+(GS_Y-1)*GS_X+k*GS_XY;
			count=0;
			ix1=j;
			if(_Mask[GridPoint1]<0)
			{
				ix2=ix1+1;
				GridPoint2=i+ix2*GS_X+k*GS_XY;
				
				while(_Mask[GridPoint2]==MaskPMFInterpol&&GridPoint2<GridPointEnd)
				{
					ix2++;
					GridPoint2=i+ix2*GS_X+k*GS_XY;
				}
				if(GridPoint2<GridPointEnd&&ix2-ix1<=Rstep)
				{
					x1=(float)ix1;
					x2=(float)ix2;
					for(ix=x1;ix<x2;ix++)
					{
						x=(float)ix;
						GridPoint=i+ix*GS_X+k*GS_XY;
						V[GridPoint]=V[GridPoint1]+(V[GridPoint2]-V[GridPoint1])*(x-x1)/(x2-x1);
						//_Mask[GridPoint]=MaskPMFInterpol;
					}
				}
			}
		}
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint1=i+j*GS_X+k*GS_XY;
			GridPointEnd=i+j*GS_X+(GS_Z-1)*GS_XY;
			count=0;
			ix1=k;
			if(_Mask[GridPoint1]<0)
			{
				ix2=ix1+1;
				GridPoint2=i+j*GS_X+ix2*GS_XY;
				
				while(_Mask[GridPoint2]==MaskPMFInterpol&&GridPoint2<GridPointEnd)
				{
					ix2++;
					GridPoint2=i+j*GS_X+ix2*GS_XY;
				}
				if(GridPoint2<GridPointEnd&&ix2-ix1<=Rstep)
				{
					x1=(float)ix1;
					x2=(float)ix2;
					for(ix=x1;ix<x2;ix++)
					{
						x=(float)ix;
						GridPoint=i+j*GS_X+ix*GS_XY;
						V[GridPoint]=V[GridPoint1]+(V[GridPoint2]-V[GridPoint1])*(x-x1)/(x2-x1);
						//_Mask[GridPoint]=MaskPMFInterpol;
					}
				}
			}
		}
	}
	for(i=1;i<GS_X-1;i++)
		for(j=1;j<GS_Y-1;j++)
			for(k=1;k<GS_Z-1;k++)
	{
		GridPoint=i+j*GS_X+k*GS_XY;
		if(_Mask[GridPoint]<0)
		{
			_Mask[GridPoint]=-_Mask[GridPoint];
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::LineInterpolation1AllIons(double** _V,int **_Mask,int Rstep)
{
	int i;
	for(i=0;i<Nions;i++)
	{
		LineInterpolation1(_V,_Mask,i,Rstep);
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::LineInterpolation1(double** V,int **_Mask,int ion,int Rstep)
{
	pnpPrint("PMFProcessing::LineInterpolation1 ion: %d  Rstep: %d\n",ion,Rstep);
	PNP_EXIT_FAIL_NULL(V,"V is not initialize\n");
	PNP_EXIT_FAIL_NULL(_Mask,"_Mask is not initialize\n");
	
	int i,j,k,i1,i2,j2,k2,i3,j3,k3;
	int GridPoint,GridPoint1,GridPoint2,GridPointEnd;
	int count,ix,ix1,ix2;
	float x,x1,x2;
	int CountTotPntInter=0;
	
	for(i1=0;i1<1;i1++)
	{
		//X
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint1=i+j*GS_X+k*GS_XY;
			count=0;
			if(_Mask[ion][GridPoint1]==MaskPMFCalc||_Mask[ion][GridPoint1]==MaskPMFInterpol)
			{
				i2=i+1;j2=j;k2=k;
				GridPoint2=i2+j2*GS_X+k2*GS_XY;
				count=1;
				while(_Mask[ion][GridPoint2]!=MaskPMFCalc&&_Mask[ion][GridPoint2]!=MaskPMFInterpol&&count<=Rstep)
				{
					i2++;
					GridPoint2=i2+j2*GS_X+k2*GS_XY;
					count++;
				}
				i3=i2;
				if((_Mask[ion][GridPoint2]==MaskPMFCalc||_Mask[ion][GridPoint2]==MaskPMFInterpol)&&count<=Rstep)
				{
					x1=(float)i;
					x2=(float)i3;
					for(i2=i+1;i2<i3;i2++)
					{
						GridPoint=i2+j2*GS_X+k2*GS_XY;
						x=(float)i2;
						V[ion][GridPoint]=V[ion][GridPoint1]+(V[ion][GridPoint2]-V[ion][GridPoint1])*(x-x1)/(x2-x1);
						_Mask[ion][GridPoint]=MaskPMFInterpolTemp;
						CountTotPntInter++;
					}
				}
			}
		}
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(_Mask[ion][GridPoint]==MaskPMFInterpolTemp)
				_Mask[ion][GridPoint]=MaskPMFInterpol;
		}
		//Y
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint1=i+j*GS_X+k*GS_XY;
			count=0;
			if(_Mask[ion][GridPoint1]==MaskPMFCalc||_Mask[ion][GridPoint1]==MaskPMFInterpol)
			{
				i2=i;j2=j+1;k2=k;
				GridPoint2=i2+j2*GS_X+k2*GS_XY;
				count=1;
				while(_Mask[ion][GridPoint2]!=MaskPMFCalc&&_Mask[ion][GridPoint2]!=MaskPMFInterpol&&count<=Rstep)
				{
					j2++;
					GridPoint2=i2+j2*GS_X+k2*GS_XY;
					count++;
				}
				j3=j2;
				if((_Mask[ion][GridPoint2]==MaskPMFCalc||_Mask[ion][GridPoint2]==MaskPMFInterpol)&&count<=Rstep)
				{
					x1=(float)j;
					x2=(float)j3;
					for(j2=j+1;j2<j3;j2++)
					{
						GridPoint=i2+j2*GS_X+k2*GS_XY;
						x=(float)j2;
						V[ion][GridPoint]=V[ion][GridPoint1]+(V[ion][GridPoint2]-V[ion][GridPoint1])*(x-x1)/(x2-x1);
						_Mask[ion][GridPoint]=MaskPMFInterpolTemp;
						CountTotPntInter++;
					}
				}
			}
		}
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(_Mask[ion][GridPoint]==MaskPMFInterpolTemp)
				_Mask[ion][GridPoint]=MaskPMFInterpol;
		}
		//Y
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint1=i+j*GS_X+k*GS_XY;
			count=0;
			if(_Mask[ion][GridPoint1]==MaskPMFCalc||_Mask[ion][GridPoint1]==MaskPMFInterpol)
			{
				i2=i;j2=j;k2=k+1;
				GridPoint2=i2+j2*GS_X+k2*GS_XY;
				count=1;
				while(_Mask[ion][GridPoint2]!=MaskPMFCalc&&_Mask[ion][GridPoint2]!=MaskPMFInterpol&&count<=Rstep)
				{
					k2++;
					GridPoint2=i2+j2*GS_X+k2*GS_XY;
					count++;
				}
				k3=k2;
				if((_Mask[ion][GridPoint2]==MaskPMFCalc||_Mask[ion][GridPoint2]==MaskPMFInterpol)&&count<=Rstep)
				{
					x1=(float)k;
					x2=(float)k3;
					for(k2=k+1;k2<k3;k2++)
					{
						GridPoint=i2+j2*GS_X+k2*GS_XY;
						x=(float)k2;
						V[ion][GridPoint]=V[ion][GridPoint1]+(V[ion][GridPoint2]-V[ion][GridPoint1])*(x-x1)/(x2-x1);
						_Mask[ion][GridPoint]=MaskPMFInterpolTemp;
						CountTotPntInter++;
					}
				}
			}
			
		}
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(_Mask[ion][GridPoint]==MaskPMFInterpolTemp)
				_Mask[ion][GridPoint]=MaskPMFInterpol;
		}
	}
	pnpPrint("Total Points interpolated: %d\n",CountTotPntInter);
	return CountTotPntInter;
}
int PMFProcessing::LineInterpolation1_VF(VectorField3D *VF_PMF,VectorIntField3D *VF_Mask,int ion,int Rstep)
{
	float **V=VF_PMF->V;
	int **_Mask=VF_Mask->V;
	pnpPrint("PMFProcessing::LineInterpolation1 ion: %d  Rstep: %d\n",ion,Rstep);
	PNP_EXIT_FAIL_NULL(V,"V is not initialize\n");
	PNP_EXIT_FAIL_NULL(_Mask,"_Mask is not initialize\n");
	
	int i,j,k,i1,i2,j2,k2,i3,j3,k3;
	int GridPoint,GridPoint1,GridPoint2,GridPointEnd;
	int count,ix,ix1,ix2;
	float x,x1,x2;
	int CountTotPntInter=0;
	
	for(i1=0;i1<1;i1++)
	{
		//X
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint1=i+j*GS_X+k*GS_XY;
			count=0;
			if(_Mask[ion][GridPoint1]==MaskPMFCalc||_Mask[ion][GridPoint1]==MaskPMFInterpol)
			{
				i2=i+1;j2=j;k2=k;
				GridPoint2=i2+j2*GS_X+k2*GS_XY;
				count=1;
				while(_Mask[ion][GridPoint2]!=MaskPMFCalc&&_Mask[ion][GridPoint2]!=MaskPMFInterpol&&count<=Rstep)
				{
					i2++;
					GridPoint2=i2+j2*GS_X+k2*GS_XY;
					count++;
				}
				i3=i2;
				if((_Mask[ion][GridPoint2]==MaskPMFCalc||_Mask[ion][GridPoint2]==MaskPMFInterpol)&&count<=Rstep)
				{
					x1=(float)i;
					x2=(float)i3;
					for(i2=i+1;i2<i3;i2++)
					{
						GridPoint=i2+j2*GS_X+k2*GS_XY;
						x=(float)i2;
						V[ion][GridPoint]=V[ion][GridPoint1]+(V[ion][GridPoint2]-V[ion][GridPoint1])*(x-x1)/(x2-x1);
						_Mask[ion][GridPoint]=MaskPMFInterpolTemp;
						CountTotPntInter++;
					}
				}
			}
		}
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(_Mask[ion][GridPoint]==MaskPMFInterpolTemp)
				_Mask[ion][GridPoint]=MaskPMFInterpol;
		}
		//Y
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint1=i+j*GS_X+k*GS_XY;
			count=0;
			if(_Mask[ion][GridPoint1]==MaskPMFCalc||_Mask[ion][GridPoint1]==MaskPMFInterpol)
			{
				i2=i;j2=j+1;k2=k;
				GridPoint2=i2+j2*GS_X+k2*GS_XY;
				count=1;
				while(_Mask[ion][GridPoint2]!=MaskPMFCalc&&_Mask[ion][GridPoint2]!=MaskPMFInterpol&&count<=Rstep)
				{
					j2++;
					GridPoint2=i2+j2*GS_X+k2*GS_XY;
					count++;
				}
				j3=j2;
				if((_Mask[ion][GridPoint2]==MaskPMFCalc||_Mask[ion][GridPoint2]==MaskPMFInterpol)&&count<=Rstep)
				{
					x1=(float)j;
					x2=(float)j3;
					for(j2=j+1;j2<j3;j2++)
					{
						GridPoint=i2+j2*GS_X+k2*GS_XY;
						x=(float)j2;
						V[ion][GridPoint]=V[ion][GridPoint1]+(V[ion][GridPoint2]-V[ion][GridPoint1])*(x-x1)/(x2-x1);
						_Mask[ion][GridPoint]=MaskPMFInterpolTemp;
						CountTotPntInter++;
					}
				}
			}
		}
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(_Mask[ion][GridPoint]==MaskPMFInterpolTemp)
				_Mask[ion][GridPoint]=MaskPMFInterpol;
		}
		//Y
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint1=i+j*GS_X+k*GS_XY;
			count=0;
			if(_Mask[ion][GridPoint1]==MaskPMFCalc||_Mask[ion][GridPoint1]==MaskPMFInterpol)
			{
				i2=i;j2=j;k2=k+1;
				GridPoint2=i2+j2*GS_X+k2*GS_XY;
				count=1;
				while(_Mask[ion][GridPoint2]!=MaskPMFCalc&&_Mask[ion][GridPoint2]!=MaskPMFInterpol&&count<=Rstep)
				{
					k2++;
					GridPoint2=i2+j2*GS_X+k2*GS_XY;
					count++;
				}
				k3=k2;
				if((_Mask[ion][GridPoint2]==MaskPMFCalc||_Mask[ion][GridPoint2]==MaskPMFInterpol)&&count<=Rstep)
				{
					x1=(float)k;
					x2=(float)k3;
					for(k2=k+1;k2<k3;k2++)
					{
						GridPoint=i2+j2*GS_X+k2*GS_XY;
						x=(float)k2;
						V[ion][GridPoint]=V[ion][GridPoint1]+(V[ion][GridPoint2]-V[ion][GridPoint1])*(x-x1)/(x2-x1);
						_Mask[ion][GridPoint]=MaskPMFInterpolTemp;
						CountTotPntInter++;
					}
				}
			}
			
		}
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(_Mask[ion][GridPoint]==MaskPMFInterpolTemp)
				_Mask[ion][GridPoint]=MaskPMFInterpol;
		}
	}
	pnpPrint("Total Points interpolated: %d\n",CountTotPntInter);
	return CountTotPntInter;
}
int PMFProcessing::LineInterpolationWithDiff(double** V,int **_Mask,int ion,int Rstep)
{
	pnpPrint("PMFProcessing::LineInterpolationWithDiff ion: %d  Rstep: %d\n",ion,Rstep);
	PNP_EXIT_FAIL_NULL(V,"V is not initialize\n");
	PNP_EXIT_FAIL_NULL(_Mask,"V is not initialize\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld,"m_ContWorld is not initialize\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld->NIndexing,"m_ContWorld->NIndexing is not initialize\n");
	int i,j,k,i1,i2,j2,k2,i3,j3,k3;
	int GridPoint,GridPoint1,GridPoint2,GridPointEnd;
	int count,ix,ix1,ix2;
	float x,x1,x2;
	int CountTotPntInter=0;
	
	for(i1=0;i1<1;i1++)
	{
		//X
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint1=i+j*GS_X+k*GS_XY;
			count=0;
			if((_Mask[ion][GridPoint1]==MaskPMFCalc||_Mask[ion][GridPoint1]==MaskPMFInterpol)&&m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint1)>0.0)
			{
				i2=i+1;j2=j;k2=k;
				GridPoint2=i2+j2*GS_X+k2*GS_XY;
				count=1;
				while(_Mask[ion][GridPoint2]!=MaskPMFCalc&&_Mask[ion][GridPoint2]!=MaskPMFInterpol&&count<=Rstep&&m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint2)>0.0)
				{
					i2++;
					GridPoint2=i2+j2*GS_X+k2*GS_XY;
					count++;
				}
				i3=i2;
				if((_Mask[ion][GridPoint2]==MaskPMFCalc||_Mask[ion][GridPoint2]==MaskPMFInterpol)&&count<=Rstep&&m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint2)>0.0)
				{
					x1=(float)i;
					x2=(float)i3;
					for(i2=i+1;i2<i3;i2++)
					{
						GridPoint=i2+j2*GS_X+k2*GS_XY;
						x=(float)i2;
						V[ion][GridPoint]=V[ion][GridPoint1]+(V[ion][GridPoint2]-V[ion][GridPoint1])*(x-x1)/(x2-x1);
						_Mask[ion][GridPoint]=MaskPMFInterpolTemp;
						CountTotPntInter++;
					}
				}
			}
		}
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(_Mask[ion][GridPoint]==MaskPMFInterpolTemp)
				_Mask[ion][GridPoint]=MaskPMFInterpol;
		}
		//Y
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint1=i+j*GS_X+k*GS_XY;
			count=0;
			if((_Mask[ion][GridPoint1]==MaskPMFCalc||_Mask[ion][GridPoint1]==MaskPMFInterpol)&&m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint1)>0.0)
			{
				i2=i;j2=j+1;k2=k;
				GridPoint2=i2+j2*GS_X+k2*GS_XY;
				count=1;
				while(_Mask[ion][GridPoint2]!=MaskPMFCalc&&_Mask[ion][GridPoint2]!=MaskPMFInterpol&&count<=Rstep&&m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint2)>0.0)
				{
					j2++;
					GridPoint2=i2+j2*GS_X+k2*GS_XY;
					count++;
				}
				j3=j2;
				if((_Mask[ion][GridPoint2]==MaskPMFCalc||_Mask[ion][GridPoint2]==MaskPMFInterpol)&&count<=Rstep&&m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint2)>0.0)
				{
					x1=(float)j;
					x2=(float)j3;
					for(j2=j+1;j2<j3;j2++)
					{
						GridPoint=i2+j2*GS_X+k2*GS_XY;
						x=(float)j2;
						V[ion][GridPoint]=V[ion][GridPoint1]+(V[ion][GridPoint2]-V[ion][GridPoint1])*(x-x1)/(x2-x1);
						_Mask[ion][GridPoint]=MaskPMFInterpolTemp;
						CountTotPntInter++;
					}
				}
			}
		}
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(_Mask[ion][GridPoint]==MaskPMFInterpolTemp)
				_Mask[ion][GridPoint]=MaskPMFInterpol;
		}
		//Y
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint1=i+j*GS_X+k*GS_XY;
			count=0;
			if((_Mask[ion][GridPoint1]==MaskPMFCalc||_Mask[ion][GridPoint1]==MaskPMFInterpol)&&m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint1)>0.0)
			{
				i2=i;j2=j;k2=k+1;
				GridPoint2=i2+j2*GS_X+k2*GS_XY;
				count=1;
				while(_Mask[ion][GridPoint2]!=MaskPMFCalc&&_Mask[ion][GridPoint2]!=MaskPMFInterpol&&count<=Rstep&&m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint2)>0.0)
				{
					k2++;
					GridPoint2=i2+j2*GS_X+k2*GS_XY;
					count++;
				}
				k3=k2;
				if((_Mask[ion][GridPoint2]==MaskPMFCalc||_Mask[ion][GridPoint2]==MaskPMFInterpol)&&count<=Rstep&&m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint2)>0.0)
				{
					x1=(float)k;
					x2=(float)k3;
					for(k2=k+1;k2<k3;k2++)
					{
						GridPoint=i2+j2*GS_X+k2*GS_XY;
						x=(float)k2;
						V[ion][GridPoint]=V[ion][GridPoint1]+(V[ion][GridPoint2]-V[ion][GridPoint1])*(x-x1)/(x2-x1);
						_Mask[ion][GridPoint]=MaskPMFInterpolTemp;
						CountTotPntInter++;
					}
				}
			}
			
		}
		for(i=1;i<GS_X-1;i++)
			for(j=1;j<GS_Y-1;j++)
				for(k=1;k<GS_Z-1;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(_Mask[ion][GridPoint]==MaskPMFInterpolTemp)
				_Mask[ion][GridPoint]=MaskPMFInterpol;
		}
	}
	pnpPrint("Total Points interpolated: %d\n",CountTotPntInter);
	return CountTotPntInter;
}
int PMFProcessing::LineInterpolation(double* V,int *_Mask)
{
	int i,j,k,i1,i2,j2,k2,i3,j3,k3;
	int GridPoint,GridPoint1,GridPoint2,GridPointEnd;
	int count;
	float x,x1,x2;
	for(i=1;i<GS_X-1;i++)
		for(j=1;j<GS_Y-1;j++)
			for(k=1;k<GS_Z-1;k++)
	{	
		count=0;
		if(_Mask[GridPoint1]==1)
		{
			GridPoint2=GridPoint1+1;
			while(_Mask[GridPoint2]!=1&&GridPoint2<GridPointEnd)
			{
				GridPoint2++;
			}
			if(GridPoint2<GridPointEnd)
			{
				x1=(float)GridPoint1;
				x2=(float)GridPoint2;
				for(GridPoint=GridPoint1+1;GridPoint<GridPoint2;GridPoint++)
				{
					x=(float)GridPoint;
					V[GridPoint]=V[GridPoint1]+(V[GridPoint]-V[GridPoint1])*(x-x1)/(x2-x1);
				}
			}
		}
	}/*
	for(i=1;i<GS_X-1;i++)
		for(j=1;j<GS_Y-1;j++)
			for(k=1;k<GS_Z-1;k++)
	{	for(i=1;i<GS_X-1;i++)
		for(j=1;j<GS_Y-1;j++)
			for(k=1;k<GS_Z-1;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
		GridPoint=i+j*GS_X+k*GS_XY;
		
		count=0;
		if(D[ion][GridPoint]!=0&&V[GridPoint]!=0)
		{
			i2=GridPoint+GridSize[0];
			j3=j;
			GridPointEnd=i+(GridSize[1]-1)*GridSize[0]+k*GridSize[0]*GridSize[1];
			while(D[ion][i2]!=0&&V[i2]==0&&i2<GridPointEnd)
			{
				i2+=GridSize[0];
				j3++;
			}
			if(i2<GridPointEnd)
			{
				if(D[ion][i2]==0)
				{
					for(i3=GridPoint+GridSize[0];i3<i2;i3=i3+GridSize[0])
					{
						V[i3]=V[GridPoint];
					}
				}
				else if(V[i2]!=0)
				{
					z1=(float)j;
					z2=(float)j3;
					k3=j+1;
					for(i3=GridPoint+GridSize[0];i3<i2;i3=i3+GridSize[0])
					{
						x=(float)k3;
						V[i3]=V[GridPoint]+(V[i2]-V[GridPoint])*(x-z1)/(z2-z1);
						k3++;
					}
				}
			}
		}
	}
	for(i=1;i<GridSize[0]-1;i++)for(j=1;j<GridSize[1]-1;j++)for(k=1;k<GridSize[2]-1;k++)
	{
		GridPoint=i+j*GridSize[0]+k*GridSize[0]*GridSize[1];
		
		count=0;
		if(D[ion][GridPoint]!=0&&V[GridPoint]!=0)
		{
			i2=GridPoint+GridSize[0]*GridSize[1];
			j3=k;
			GridPointEnd=i+j*GridSize[0]+(GridSize[2]-1)*GridSize[0]*GridSize[1];
			while(D[ion][i2]!=0&&V[i2]==0&&i2<GridPointEnd)
			{
				i2+=GridSize[0]*GridSize[1];
				j3++;
			}
			if(i2<GridPointEnd)
			{
				if(D[ion][i2]==0)
				{
					for(i3=GridPoint+GridSize[0]*GridSize[1];i3<i2;i3=i3+GridSize[0]*GridSize[1])
					{
						V[i3]=V[GridPoint];
					}
				}
				else if(V[i2]!=0)
				{
					z1=(float)k;
					z2=(float)j3;
					k3=k+1;
					for(i3=GridPoint+GridSize[0]*GridSize[1];i3<i2;i3=i3+GridSize[0]*GridSize[1])
					{
						x=(float)k3;
						V[i3]=V[GridPoint]+(V[i2]-V[GridPoint])*(float)(x-z1)/(float)(z2-z1);
						k3++;
					}
				}
			}
		}
	}*/
	return EXIT_SUCCESS;
}
int PMFProcessing::AverageThrLaplasAllIons(double** V,int **_Mask,int iter,int MaskNotToDo)
{
	int i;
	for(i=0;i<Nions;i++)
	{
		AverageThrLaplas(V[i],_Mask[i],iter,MaskNotToDo);
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::AverageThrLaplas(double* V,int *_Mask,int iter,int MaskNotToDo)
{
  //!Not very right now. Solve FDM laplasian iter times for RF field
	DbgPrint2("RFAverage. iter:%d\n",iter);
	int i,j,k,i1,GridPoint,BnW;
	
	for(i1=0;i1<iter;i1++)
	{
		for(BnW=0;BnW<2;BnW++)
			for(GridPoint=BnW;GridPoint<GS_XYZ;GridPoint=GridPoint+2)
		{
			i=GridPoint%GS_X;
			k=GridPoint/GS_XY;
			j=(GridPoint%GS_XY)/GS_X;
			
			if(i>0&&j>0&&k>0&&i<GS_X-1&&j<GS_Y-1&&k<GS_Z-1)
				if(_Mask[GridPoint]!=MaskNotToDo)
			{
				V[GridPoint]=(V[GridPoint]*6.0+V[GridPoint+1]+V[GridPoint-1]+V[GridPoint-GS_X]+V[GridPoint+GS_X]+V[GridPoint+GS_XY]+V[GridPoint-GS_XY])/12.0;
			}
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::ConvertIndex(int **_Mask,int MaskFrom,int MaskTo)
{
	int i,GridPoint;
	for(i=0;i<Nions;i++)
	for(GridPoint=0;GridPoint<GS_XYZ;GridPoint++)
	{
		if(_Mask[i][GridPoint]==MaskFrom)
			_Mask[i][GridPoint]=MaskTo;
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::ConvertIndexWithinIntZ(int **_Mask,int MaskFrom,int MaskTo,int iZ0,int iZ1)
{
	int ion,i,j,k,GridPoint;
	for(ion=0;ion<Nions;ion++)
		for(i=0;i<GS_X;i++)
			for(j=0;j<GS_Y;j++)
				for(k=iZ0;k<=iZ1;k++)
	{
		GridPoint=i+j*GS_X+k*GS_XY;
		if(_Mask[i][GridPoint]==MaskFrom)
			_Mask[i][GridPoint]=MaskTo;
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::AnalizePMF(double** V,int **_Mask)
{
	return EXIT_SUCCESS;
}
int PMFProcessing::RemoveBadDiffusionPoints()
{
	m_ContWorld->NIndexing->RemoveBadDiffusionPoints();
	return EXIT_SUCCESS;
}
int PMFProcessing::RemoveDiffusionPointsAtHighC(float HighC)
{
	int IType,i,k,j,GrdPnt;
	PNP_EXIT_FAIL_NULL(m_ContWorld,"m_ContWorld is not initialized\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld->C,"concentration is not initialize at m_ContWorld\n");
	for(IType=0;IType<m_ContWorld->NIonsTypes;IType++)
		PNP_EXIT_FAIL_NULL(m_ContWorld->C[IType],"concentration is not initialize at m_ContWorld\n");
	
	float fpoh= 4*M_PI*m_ContWorld->GridScale;
	float coef=fpoh*COANGS / (m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
	HighC=HighC*coef;
	
	int DiffZeroInd=m_ContWorld->NIndexing->GetDiffZeroInd();
	int count=0;
	for(IType=0;IType<m_ContWorld->NIonsTypes;IType++)
	{
		for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
		{
			if(m_ContWorld->NIndexing->GetDiffFloat(IType,GrdPnt)!=0.0)
			{
				if(m_ContWorld->C[IType][GrdPnt]>HighC)
				{
					m_ContWorld->NIndexing->SetIonField(IType,GrdPnt,DiffZeroInd);
					count++;
				}
			}
		}
	}
	if(count>0)m_ContWorld->NIndexing->CalcDiffBoarder();
	pnpPrint("%d points was remoove becouse of its concentration was > %g\n",count,HighC/coef);
	return EXIT_SUCCESS;
}
int PMFProcessing::CorrectFlexDiffusionWithNIndexDiff()
{
	int IType,i,k,j,GrdPnt;
	PNP_EXIT_FAIL_NULL(m_ContWorld,"m_ContWorld is not initialized\n");
	PNP_EXIT_FAIL_NULL(m_ContWorld->D,"flexible diffusion is not initialize at m_ContWorld\n");
	for(IType=0;IType<m_ContWorld->NIonsTypes;IType++)
		PNP_EXIT_FAIL_NULL(m_ContWorld->D[IType],"flexible diffusion is not initialize at m_ContWorld\n");
	for(IType=0;IType<m_ContWorld->NIonsTypes;IType++)
	{
		for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
		{
			if(m_ContWorld->NIndexing->GetDiffFloat(IType,GrdPnt)==0.0)
			{
				m_ContWorld->D[IType][GrdPnt]=0.0;
			}
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::RemovePMFwithHighdPMF(double** V,int **_Mask,double dPMF,int MaskToSet,int cycles)
{
	pnpPrint0("PNPUtil::RemoveLargePMFfromPNP\n");
	int i,j,k,ion;
	int GridPoint;
	double dPMFx;
	
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
							_Mask[ion][GridPoint]=MaskToSet;
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+1)==0.0)
					{//boarder at +x
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint-1]);
						if(dPMFx>=dPMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							_Mask[ion][GridPoint]=MaskToSet;
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_X)==0.0)
					{//boarder at -y
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint+GS_X]);
						if(dPMFx>=dPMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							_Mask[ion][GridPoint]=MaskToSet;
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_X)==0.0)
					{//boarder at +y
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint-GS_X]);
						if(dPMFx>=dPMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							_Mask[ion][GridPoint]=MaskToSet;
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_XY)==0.0)
					{//boarder at -z
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint+GS_XY]);
						if(dPMFx>=dPMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							_Mask[ion][GridPoint]=MaskToSet;
							Count++;
						}
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_XY)==0.0)
					{//boarder at +z
						dPMFx=fabs(V[ion][GridPoint]-V[ion][GridPoint-GS_XY]);
						if(dPMFx>=dPMF)
						{
							m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
							_Mask[ion][GridPoint]=MaskToSet;
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
int PMFProcessing::RemoveLargePMFfromPNP(double** V,int **_Mask,double LargePMF,int MaskToSet)
{
	pnpPrint0("PNPUtil::RemoveLargePMFfromPNP\n");
	int i,j,k,ion;
	int GridPoint;
	
	for(ion=0;ion<m_ContWorld->NIonsTypes;ion++)
	{
		//float* D=NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion0);
		for(i=0;i<GS_X;i++)
			for(j=0;j<GS_Y;j++)
				for(k=0;k<GS_Z;k++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint)!=0.0
				&&V[ion][GridPoint]>LargePMF)
			{
				m_ContWorld->NIndexing->SetDiffToZero(ion,GridPoint);
				_Mask[ion][GridPoint]=MaskToSet;
			}
		}
	}
	m_ContWorld->NIndexing->CalcDiffBoarder();
	return EXIT_SUCCESS;
}
int PMFProcessing::RemoveLargePMF(double** V,int **_Mask,int ion,double PMFVal)
{
	int i,j,k;
	int GridPoint;
	for(i=0;i<GS_X;i++)
		for(j=0;j<GS_Y;j++)
			for(k=0;k<GS_Z;k++)
	{
		GridPoint=i+j*GS_X+k*GS_XY;
		if(V[ion][GridPoint]>PMFVal)
		{
			V[ion][GridPoint]=PMFVal;
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::AverageThrLaplasDAllIons(double** V,int **_Mask,int iter,int MaskAmong)
{
	int ion;
	for(ion=0;ion<Nions;ion++)
	{
		int i,j,k,i1,GridPoint,BnW;
		float count;
	
		for(i1=0;i1<iter;i1++)
		{
			for(BnW=0;BnW<2;BnW++)
				for(GridPoint=BnW;GridPoint<GS_XYZ;GridPoint=GridPoint+2)
			{
				i=GridPoint%GS_X;
				k=GridPoint/GS_XY;
				j=(GridPoint%GS_XY)/GS_X;
			
			//GridPoint=i+j*GS_X+k*GS_XY;
				count=6.0;
			
				if(_Mask[ion][GridPoint]==MaskAmong)
					if(i>0&&j>0&&k>0&&i<GS_X-1&&j<GS_Y-1&&k<GS_Z-1)
				{
					V[ion][GridPoint]*=6.0;
					//if(_Mask[GridPoint+1]!=0)
					//{
					V[ion][GridPoint]+=V[ion][GridPoint+1];
						count++;
					//}
					//if(_Mask[GridPoint-1]!=0)
					//{
						V[ion][GridPoint]+=V[ion][GridPoint-1];
						count++;
					//}
					//if(_Mask[GridPoint+GS_X]!=0)
					//{
						V[ion][GridPoint]+=V[ion][GridPoint+GS_X];
						count++;
					//}
					//if(_Mask[GridPoint-GS_X]!=0)
					//{
						V[ion][GridPoint]+=V[ion][GridPoint-GS_X];
						count++;
					//}
					//if(_Mask[GridPoint+GS_XY]!=0)
					//{
						V[ion][GridPoint]+=V[ion][GridPoint+GS_XY];
						count++;
					//}
					//if(_Mask[GridPoint-GS_XY]!=0)
					//{
						V[ion][GridPoint]+=V[ion][GridPoint-GS_XY];
						count++;
					//}
						V[ion][GridPoint]/=count;
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::AverageThrLaplasDnoZeroAllIons(double** V,int **_Mask,int iter,int MaskAmong)
{
	int ion;
	for(ion=0;ion<Nions;ion++)
	{
		int i,j,k,i1,GridPoint,BnW;
		float count;
	
		for(i1=0;i1<iter;i1++)
		{
			for(BnW=0;BnW<2;BnW++)
				for(GridPoint=BnW;GridPoint<GS_XYZ;GridPoint=GridPoint+2)
			{
				i=GridPoint%GS_X;
				k=GridPoint/GS_XY;
				j=(GridPoint%GS_XY)/GS_X;
			
			//GridPoint=i+j*GS_X+k*GS_XY;
				count=6.0;
			
				if(_Mask[ion][GridPoint]==MaskAmong)
					if(i>0&&j>0&&k>0&&i<GS_X-1&&j<GS_Y-1&&k<GS_Z-1)
				{
					V[ion][GridPoint]*=6.0;
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+1)!=0.0)
					{
						V[ion][GridPoint]+=V[ion][GridPoint+1];
						count++;
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-1)!=0.0)
					{
						V[ion][GridPoint]+=V[ion][GridPoint-1];
						count++;
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_X)!=0.0)
					{
						V[ion][GridPoint]+=V[ion][GridPoint+GS_X];
						count++;
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_X)!=0.0)
					{
						V[ion][GridPoint]+=V[ion][GridPoint-GS_X];
						count++;
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+GS_XY)!=0.0)
					{
						V[ion][GridPoint]+=V[ion][GridPoint+GS_XY];
						count++;
					}
					if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-GS_XY)!=0.0)
					{
						V[ion][GridPoint]+=V[ion][GridPoint-GS_XY];
						count++;
					}
					V[ion][GridPoint]/=count;
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::AverageThrLaplas1AllIons(double** V,int **_Mask,int iter,int MaskAmong)
{
	int i;
	for(i=0;i<Nions;i++)
	{
		AverageThrLaplas1(V[i],_Mask[i],iter,MaskAmong);
	}
	return EXIT_SUCCESS;
}
int PMFProcessing::AverageThrLaplas1(double* V,int *_Mask,int iter,int MaskAmong)
{
  //!Not very right now. Solve FDM laplasian iter times for RF field
	DbgPrint2("AverageThrLaplas1. iter:%d\n",iter);
	int i,j,k,i1,GridPoint,BnW;
	float count;
	
	for(i1=0;i1<iter;i1++)
	{
		for(BnW=0;BnW<2;BnW++)
			for(GridPoint=BnW;GridPoint<GS_XYZ;GridPoint=GridPoint+2)
		{
			i=GridPoint%GS_X;
			k=GridPoint/GS_XY;
			j=(GridPoint%GS_XY)/GS_X;
			
			//GridPoint=i+j*GS_X+k*GS_XY;
			count=6.0;
			
			if(_Mask[GridPoint]==MaskAmong)
				if(i>0&&j>0&&k>0&&i<GS_X-1&&j<GS_Y-1&&k<GS_Z-1)
			{
				V[GridPoint]*=6.0;
				if(_Mask[GridPoint+1]!=0)
				{
					V[GridPoint]+=V[GridPoint+1];
					count++;
				}
				if(_Mask[GridPoint-1]!=0)
				{
					V[GridPoint]+=V[GridPoint-1];
					count++;
				}
				if(_Mask[GridPoint+GS_X]!=0)
				{
					V[GridPoint]+=V[GridPoint+GS_X];
					count++;
				}
				if(_Mask[GridPoint-GS_X]!=0)
				{
					V[GridPoint]+=V[GridPoint-GS_X];
					count++;
				}
				if(_Mask[GridPoint+GS_XY]!=0)
				{
					V[GridPoint]+=V[GridPoint+GS_XY];
					count++;
				}
				if(_Mask[GridPoint-GS_XY]!=0)
				{
					V[GridPoint]+=V[GridPoint-GS_XY];
					count++;
				}
				V[GridPoint]/=count;
			}
		}
	}
	return EXIT_SUCCESS;
}

SIPPointsSearch::SIPPointsSearch()
{
	InitZero();
}
SIPPointsSearch::~SIPPointsSearch()
{
	Clear();
}
int SIPPointsSearch::InitZero()
{
	return EXIT_SUCCESS;
}
int SIPPointsSearch::Clear()
{
	return EXIT_SUCCESS;
}
int SIPPointsSearch::SetContWorld(ContWorld* _ContWorld)
{
	m_ContWorld=_ContWorld;
	return EXIT_SUCCESS;
}
int SIPPointsSearch::CmdSIPPointsSearch(const TiXmlElement *Elt)
{
	
	if(strcmp("ChannelPoints",Elt->Attribute("Type"))==0)
	{
		float x,y,z[2];
		float rstep;
		if(Elt->GetFloatAttribute("x",&x)==EXIT_FAILURE)
			x=0.0f;
		if(Elt->GetFloatAttribute("y",&y)==EXIT_FAILURE)
			y=0.0f;
		if(Elt->GetArrOfFloatAttribute("z",z,2)==EXIT_FAILURE)
		{
			z[0]=-100.0f;
			z[1]=100.0f;
		}
		if(Elt->GetFloatAttribute("rstep",&rstep)==EXIT_FAILURE)
			rstep=0.0f;
		if(Elt->Attribute("FileName")!=NULL)
			ChannelPoints(0,x,y,z[0],z[1],rstep,Elt->Attribute("FileName"));
		if(Elt->Attribute("FileName0")!=NULL)
			ChannelPoints(0,x,y,z[0],z[1],rstep,Elt->Attribute("FileName0"));
		if(Elt->Attribute("FileName1")!=NULL)
			ChannelPoints(1,x,y,z[0],z[1],rstep,Elt->Attribute("FileName1"));
	}
	else if(strcmp("AlongZ",Elt->Attribute("Type"))==0)
	{
		float x,y,z[2];
		float rstep;
		if(Elt->GetFloatAttribute("x",&x)==EXIT_FAILURE)
			x=0.0f;
		if(Elt->GetFloatAttribute("y",&y)==EXIT_FAILURE)
			y=0.0f;
		if(Elt->GetArrOfFloatAttribute("z",z,2)==EXIT_FAILURE)
		{
			z[0]=-100.0f;
			z[1]=100.0f;
		}
		if(Elt->GetFloatAttribute("rstep",&rstep)==EXIT_FAILURE)
			rstep=0.0f;
		AlongZ(x,y,z[0],z[1],rstep,Elt->Attribute("FileName"));
	}
	else if(strcmp("DZeroBoarder",Elt->Attribute("Type"))==0)
	{
		float x,y,z[2];
		float r;
		int ion;
		if(Elt->GetFloatAttribute("x",&x)==EXIT_FAILURE)
			x=0.0f;
		if(Elt->GetFloatAttribute("y",&y)==EXIT_FAILURE)
			y=0.0f;
		if(Elt->GetArrOfFloatAttribute("z",z,2)==EXIT_FAILURE)
		{
			z[0]=-100.0f;
			z[1]=100.0f;
		}
		if(Elt->GetFloatAttribute("r",&r)==EXIT_FAILURE)
			r=m_ContWorld->GridSize[0]/m_ContWorld->GridScale/2.0;
		if(Elt->GetIntAttribute("ion",&ion)==EXIT_FAILURE)
			ion=0;
		DZeroBoarder(x,y,z[0],z[1],r,ion,Elt->Attribute("FileName"));
	}
	else if(strcmp("Cylinder",Elt->Attribute("Type"))==0)
	{
		float x,y,z[2];
		float r,rstep;
		int ion;
		if(Elt->GetFloatAttribute("x",&x)==EXIT_FAILURE)
			x=0.0f;
		if(Elt->GetFloatAttribute("y",&y)==EXIT_FAILURE)
			y=0.0f;
		if(Elt->GetArrOfFloatAttribute("z",z,2)==EXIT_FAILURE)
		{
			z[0]=-100.0f;
			z[1]=100.0f;
		}
		if(Elt->GetFloatAttribute("r",&r)==EXIT_FAILURE)
			r=m_ContWorld->GridSize[0]/m_ContWorld->GridScale/2.0;
		if(Elt->GetIntAttribute("ion",&ion)==EXIT_FAILURE)
			ion=0;
		if(Elt->GetFloatAttribute("rstep",&rstep)==EXIT_FAILURE)
			rstep=m_ContWorld->GridSize[0]/m_ContWorld->GridScale/2.0;
		Cylinder(x,y,z[0],z[1],r,rstep,ion,Elt->Attribute("FileName"));
	}
	else if(strcmp("Box",Elt->Attribute("Type"))==0)
	{
		float x[2],y[2],z[2];
		float rstep;
		int ion;
		if(Elt->GetArrOfFloatAttribute("x",x,2)==EXIT_FAILURE)
		{
			x[0]=-m_ContWorld->GridSize[0]/m_ContWorld->GridScale/2.0;
			x[1]=-x[0];
		}
		if(Elt->GetArrOfFloatAttribute("y",y,2)==EXIT_FAILURE)
		{
			y[0]=-m_ContWorld->GridSize[1]/m_ContWorld->GridScale/2.0;
			y[1]=-y[0];
		}
		if(Elt->GetArrOfFloatAttribute("z",z,2)==EXIT_FAILURE)
		{
			z[0]=-m_ContWorld->GridSize[2]/m_ContWorld->GridScale/2.0;
			z[1]=-z[1];
		}
		if(Elt->GetFloatAttribute("rstep",&rstep)==EXIT_FAILURE)
			rstep=m_ContWorld->GridSize[0]/m_ContWorld->GridScale/2.0;
		if(Elt->GetIntAttribute("ion",&ion)==EXIT_FAILURE)
			ion=0;
		Box(x[0],y[0],z[0],x[1],y[1],z[1],rstep,ion,Elt->Attribute("FileName"));
	}
	else if(strcmp("CombineFiles",Elt->Attribute("Type"))==0)
	{
		float dr;
		if(Elt->GetFloatAttribute("dr",&dr)==EXIT_FAILURE)
			dr=0.1/m_ContWorld->GridScale;
		
		CombineFiles(Elt->Attribute("FileName0"),Elt->Attribute("FileName1"),Elt->Attribute("FileNameOut"),dr);
	}
	return EXIT_SUCCESS;
}
bool r0withinr1(float x0, float y0, float z0,float x1, float y1, float z1,float dr)
{
	float t_rSQ=(x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1);
	if(t_rSQ<=dr*dr)return true;
	return false;
}
int SIPPointsSearch::CombineFiles(const char *FileName0, const char *FileName1,const char *FileNameOut,float dr)
{
	pnpPrint0("DZeroBoarder: FileName0=%s FileName1=%s FileNameOut=%s dr=%g\n",FileName0,FileName1,FileNameOut,dr);
	int count;
	std::vector<float> r0[3];
	FILE* fpdb = fopen(FileName0,"rt");
	PNP_EXIT_FAIL_NULL1(fpdb,"Can't read %s\n",FileName0);
	char Record[1024];
	while( fgets(Record,1024,fpdb) )
	{
		if(strncmp("ATOM",Record,4)==0)
		{
			//float tmpx,tmpy,tmpz,tmpq,tmpr;
			//sscanf(Record+27,"%f %f %f %f %f",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			double tmpx,tmpy,tmpz,tmpq,tmpr;
			sscanf(Record+27,"%lf %lf %lf",&tmpx,&tmpy,&tmpz);
			//pnpPrint0("\t%f %f %f %f %f\n",tmpx,tmpy,tmpz,tmpq,tmpr);
			r0[0].push_back(float(tmpx));
			r0[1].push_back(float(tmpy));
			r0[2].push_back(float(tmpz));
		}
	}
	fclose(fpdb);
	pnpPrint0("%d atoms was read from %s\n",r0[0].size(),FileName0);
	fpdb = fopen(FileName1,"rt");
	PNP_EXIT_FAIL_NULL1(fpdb,"Can't read %s\n",FileName1);
	count=0;
	while( fgets(Record,1024,fpdb) )
	{
		if(strncmp("ATOM",Record,4)==0)
		{
			//float tmpx,tmpy,tmpz,tmpq,tmpr;
			//sscanf(Record+27,"%f %f %f %f %f",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			double tmpx,tmpy,tmpz,tmpq,tmpr;
			sscanf(Record+27,"%lf %lf %lf",&tmpx,&tmpy,&tmpz);
			//pnpPrint0("\t%f %f %f %f %f\n",tmpx,tmpy,tmpz,tmpq,tmpr);
			r0[0].push_back(float(tmpx));
			r0[1].push_back(float(tmpy));
			r0[2].push_back(float(tmpz));
			count++;
		}
	}
	fclose(fpdb);
	pnpPrint0("%d atoms was read from %s\n",count,FileName1);
	
	pnpPrint0("%d atoms before sorting\n",r0[0].size());
	//sorting etc,
	int i,j;
	
	
	//writing out
	int count1=1;
	int TotalPnts=0;
	FILE *out=fopen(FileNameOut,"w");
	
	for(i=0;i<r0[0].size();i++)
	{
		count=0;
		for(j=0;j<i;j++)
		{
			count+=r0withinr1(r0[0][i],r0[1][i],r0[2][i],r0[0][j],r0[1][j],r0[2][j],dr);
		}
		if(count==0)
		{
			fprintf(out,"ATOM  %5d %4s %c%3s%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f              \n",count1,"K  ",' ',"K  ",' ',count1,' ',r0[0][i],r0[1][i],r0[2][i],1.0,1.0);
			count1++;
			if(count1>1000)count1=1;
			TotalPnts++;
		}
	}
	fclose(out);
	pnpPrint0("Total points record: %d\n",TotalPnts);
	return EXIT_SUCCESS;
}
int SIPPointsSearch::Box(float fx0,float fy0,float fz0, float fx1,float fy1,float fz1, float frstep,int ion,const char *FileName)
{
	pnpPrint0("Box: x=[%f, %f] y=[%f, %f] z=[%f, %f] r=%f ion=%d FileName=%s\n",fx0,fx1,fy0,fy1,fz0,fz1,frstep,ion,FileName);
	
	int BoxSize[3];
	int i,j,k;
	
	for(i=0;i<3;i++)
	{
		BoxSize[i]=(int)(m_ContWorld->GridSize[i]/2);
	}
	
	fx0=fx0*m_ContWorld->GridScale+BoxSize[0];
	fx1=fx1*m_ContWorld->GridScale+BoxSize[0];
	fy0=fy0*m_ContWorld->GridScale+BoxSize[1];
	fy1=fy1*m_ContWorld->GridScale+BoxSize[1];
	fz0=fz0*m_ContWorld->GridScale+BoxSize[2];
	fz1=fz1*m_ContWorld->GridScale+BoxSize[2];
	frstep=frstep*m_ContWorld->GridScale;
	
	int x0=int(fx0+0.5);
	int x1=int(fx1+0.5);
	int y0=int(fy0+0.5);
	int y1=int(fy1+0.5);
	int z0=int(fz0+0.5);
	int z1=int(fz1+0.5);
	int rstep=int(frstep+0.5);
	
	
	int rstart[3];
	int rend[3];
	rstart[0]=x0;
	rend[0]=x1;
	rstart[1]=y0;
	rend[1]=y1;
	rstart[2]=z0;
	rend[2]=z1;
	
	for(i=0;i<3;i++)
	{
		if(rstart[i]<1)rstart[i]=1;
		if(rend[i]>m_ContWorld->GridSize[i]-2)rend[i]=m_ContWorld->GridSize[i]-2;
	}
	
	DbgPrint2("\trstart=[%d %d %d]\n",rstart[0],rstart[1],rstart[2]);
	DbgPrint2("\trend=[%d %d %d]\n",rend[0],rend[1],rend[2]);
	
	
	int GridPoint;
	double ftmpx,ftmpy,ftmpz;
	int itmpx,itmpy,itmpz;
	int gsX=m_ContWorld->GridSize[0];
	int gsXY=m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1];
	
	FILE *out=fopen(FileName,"w");
	
	int count1=1;
	int TotalPnts=0;
	for(k=rstart[2];k<=rend[2];k=k+rstep)
		for(j=rstart[1];j<=rend[1];j=j+rstep)
			for(i=rstart[0];i<=rend[0];i=i+rstep)
	{
		GridPoint=i+j*gsX+k*gsXY;
		if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint)>0.0)
		{
				itmpx=i-BoxSize[0];
				itmpy=j-BoxSize[1];
				itmpz=k-BoxSize[2];
				ftmpx=double(itmpx)/m_ContWorld->GridScale;
				ftmpy=double(itmpy)/m_ContWorld->GridScale;
				ftmpz=double(itmpz)/m_ContWorld->GridScale;
				fprintf(out,"ATOM  %5d %4s %c%3s%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f              \n",count1,"K  ",' ',"K  ",' ',count1,' ',ftmpx,ftmpy,ftmpz,1.0,1.0);
				count1++;
				if(count1>1000)count1=1;
				TotalPnts++;
		}
	}
	pnpPrint0("Total points found: %d\n",TotalPnts);
	fclose(out);
	return EXIT_SUCCESS;
}
int SIPPointsSearch::Cylinder(float fx,float fy,float fz1,float fz2,float fr,float frstep,int ion,const char *FileName)
{
	pnpPrint0("Cylinder: x=%f y=%f z=[%f, %f] r=%f rstep=%f ion=%d FileName=%s\n",fx,fy,fz1,fz2,fr,frstep,ion,FileName);
	
	int BoxSize[3];
	int i,j,k;
	
	for(i=0;i<3;i++)
	{
		BoxSize[i]=(int)(m_ContWorld->GridSize[i]/2);
	}
	
	fx=fx*m_ContWorld->GridScale+BoxSize[0];
	fy=fy*m_ContWorld->GridScale+BoxSize[1];
	fz1=fz1*m_ContWorld->GridScale+BoxSize[2];
	fz2=fz2*m_ContWorld->GridScale+BoxSize[2];
	fr=fr*m_ContWorld->GridScale;
	frstep=frstep*m_ContWorld->GridScale;
	
	int x=int(fx+0.5);
	int y=int(fy+0.5);
	int z1=int(fz1+0.5);
	int z2=int(fz2+0.5);
	int r=int(fr+0.5);
	int rstep=int(frstep+0.5);
	
	float rSQ=fr*fr;
	float tmp_rSQ;
			
	int rstart[3];
	int rend[3];
	rstart[0]=x-r;
	rend[0]=x+r;
	rstart[1]=y-r;
	rend[1]=y+r;
	rstart[2]=z1;
	rend[2]=z2;
	
	for(i=0;i<3;i++)
	{
		if(rstart[i]<1)rstart[i]=1;
		if(rend[i]>m_ContWorld->GridSize[i]-2)rend[i]=m_ContWorld->GridSize[i]-2;
	}
	
	pnpPrint0("\tx=%f y=%f z=[%f, %f] r=%f\n",fx,fy,fz1,fz2,fr);
	pnpPrint0("\tx=%d y=%d z=[%d, %d] r=%d\n",x,y,z1,z2,r);
	pnpPrint0("\trstart=[%d %d %d]\n",rstart[0],rstart[1],rstart[2]);
	pnpPrint0("\trend=[%d %d %d]\n",rend[0],rend[1],rend[2]);
	
	
	int GridPoint;
	double ftmpx,ftmpy,ftmpz;
	int itmpx,itmpy,itmpz;
	int gsX=m_ContWorld->GridSize[0];
	int gsXY=m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1];
	
	FILE *out=fopen(FileName,"w");
	
	int count1=1;
	int TotalPnts=0;
	for(k=rstart[2];k<=rend[2];k++)
		for(j=rstart[1];j<=rend[1];j++)
			for(i=rstart[0];i<=rend[0];i++)
	{
		if((i-x)%rstep==0&&(j-y)%rstep==0&&(k-z1)%rstep==0)
		{
			tmp_rSQ=(i-fx)*(i-fx)+(j-fy)*(j-fy);
			if(tmp_rSQ<=rSQ)
			{
				GridPoint=i+j*gsX+k*gsXY;
				if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint)>0.0)
				{
					itmpx=i-BoxSize[0];
					itmpy=j-BoxSize[1];
					itmpz=k-BoxSize[2];
					ftmpx=double(itmpx)/m_ContWorld->GridScale;
					ftmpy=double(itmpy)/m_ContWorld->GridScale;
					ftmpz=double(itmpz)/m_ContWorld->GridScale;
					fprintf(out,"ATOM  %5d %4s %c%3s%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f              \n",count1,"K  ",' ',"K  ",' ',count1,' ',ftmpx,ftmpy,ftmpz,1.0,1.0);
					count1++;
					if(count1>1000)count1=1;
					TotalPnts++;
				}
			}
		}
	}
	pnpPrint0("Total points found: %d\n",TotalPnts);
	fclose(out);
	return EXIT_SUCCESS;
}

int SIPPointsSearch::DZeroBoarder(float fx,float fy,float fz1,float fz2,float fr,int ion,const char *FileName)
{
	pnpPrint0("DZeroBoarder: x=%f y=%f z=[%f, %f] r=%f FileName=%s\n",fx,fy,fz1,fz2,fr,FileName);
	
	int BoxSize[3];
	int i,j,k;
	
	for(i=0;i<3;i++)
	{
		BoxSize[i]=(int)(m_ContWorld->GridSize[i]/2);
	}
	
	fx=fx*m_ContWorld->GridScale+BoxSize[0];
	fy=fy*m_ContWorld->GridScale+BoxSize[1];
	fz1=fz1*m_ContWorld->GridScale+BoxSize[2];
	fz2=fz2*m_ContWorld->GridScale+BoxSize[2];
	fr=fr*m_ContWorld->GridScale;
	
	int x=int(fx+0.5);
	int y=int(fy+0.5);
	int z1=int(fz1+0.5);
	int z2=int(fz2+0.5);
	int r=int(fr+0.5);
	
	float rSQ=fr*fr;
	float tmp_rSQ;
			
	int rstart[3];
	int rend[3];
	rstart[0]=x-r;
	rend[0]=x+r;
	rstart[1]=y-r;
	rend[1]=y+r;
	rstart[2]=z1;
	rend[2]=z2;
	
	for(i=0;i<3;i++)
	{
		if(rstart[i]<1)rstart[i]=1;
		if(rend[i]>m_ContWorld->GridSize[i]-2)rend[i]=m_ContWorld->GridSize[i]-2;
	}
	
	pnpPrint0("\tx=%f y=%f z=[%f, %f] r=%f\n",fx,fy,fz1,fz2,fr);
	pnpPrint0("\tx=%d y=%d z=[%d, %d] r=%d\n",x,y,z1,z2,r);
	pnpPrint0("\trstart=[%d %d %d]\n",rstart[0],rstart[1],rstart[2]);
	pnpPrint0("\trend=[%d %d %d]\n",rend[0],rend[1],rend[2]);
	
	
	int GridPoint;
	double ftmpx,ftmpy,ftmpz;
	int itmpx,itmpy,itmpz;
	int gsX=m_ContWorld->GridSize[0];
	int gsXY=m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1];
	
	FILE *out=fopen(FileName,"w");
	
	int count1=1;
	int TotalPnts=0;
	for(k=rstart[2];k<=rend[2];k++)
		for(j=rstart[1];j<=rend[1];j++)
			for(i=rstart[0];i<=rend[0];i++)
	{
		tmp_rSQ=(i-fx)*(i-fx)+(j-fy)*(j-fy);
		if(tmp_rSQ<=rSQ)
		{
			GridPoint=i+j*gsX+k*gsXY;
			if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint)>0.0)
			{
				
				if(m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+1)==0.0
							 ||m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-1)==0.0
							 ||m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+gsX)==0.0
							 ||m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-gsX)==0.0
							 ||m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint+gsXY)==0.0
							 ||m_ContWorld->NIndexing->GetDiffFloat(ion,GridPoint-gsXY)==0.0
					)
				{
					itmpx=i-BoxSize[0];
					itmpy=j-BoxSize[1];
					itmpz=k-BoxSize[2];
					ftmpx=double(itmpx)/m_ContWorld->GridScale;
					ftmpy=double(itmpy)/m_ContWorld->GridScale;
					ftmpz=double(itmpz)/m_ContWorld->GridScale;
					fprintf(out,"ATOM  %5d %4s %c%3s%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f              \n",count1,"K  ",' ',"K  ",' ',count1,' ',ftmpx,ftmpy,ftmpz,1.0,1.0);
					count1++;
					if(count1>1000)count1=1;
					TotalPnts++;
				}
			}
		}
	}
	pnpPrint0("Total points found: %d\n",TotalPnts);
	fclose(out);
	return EXIT_SUCCESS;
}
int SIPPointsSearch::ChannelPoints(int ion, float fx,float fy,float fz1,float fz2,float frstep,const char *FileName)
{
	pnpPrint0("ChannelPoints: x=%f y=%f z=[%f, %f] rstep=%f FileName=%s\n",fx,fy,fz1,fz2,frstep,FileName);
	vector<int> VecPoints;
	int BoxSize[3];
	int i,j,k,GridPoint,GridPoint2D,isurface,zadd;
	for(i=0;i<3;i++)
	{
		BoxSize[i]=(int)(m_ContWorld->GridSize[i]/2);
	}
	int x=(int)(fx*m_ContWorld->GridScale)+BoxSize[0];
	int y=(int)(fy*m_ContWorld->GridScale)+BoxSize[1];
	int z1=(int)(fz1*m_ContWorld->GridScale)+BoxSize[2];
	int z2=(int)(fz2*m_ContWorld->GridScale)+BoxSize[2];
	int rstep=(int)(frstep*m_ContWorld->GridScale);
  
	
	float surface,ftmp1;
	int *Map3D, *list;
	int count0,count1,count;
	NodeIndex* NIndex=m_ContWorld->NIndexing->NIndex;
	float *Dni=m_ContWorld->NIndexing->D;
	NodeIndexing* NI=m_ContWorld->NIndexing;
	
	Map3D = new int[m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1]*m_ContWorld->GridSize[2]];
	list = new int[m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1]*m_ContWorld->GridSize[2]];
	
	DbgPrint0("RFWriteRNDChannelPoints \n\tx = %d y = %d z = [%d,%d] step =%d\n",x,y,z1,z2,rstep);
	DbgPrint0("Surface Along The Z Axis in Chanel\n");
	//DbgPrint0("GridZ\tGridSurface\tZ, A\tS, A^2\n"); 
  //Build Map3D along the pore, where if D[i]==0 then Map3d[i]=0 else Map3d[i]=1
	for(k=z1;k<z2;k++)
	{
		count0=0;
		count1=1;
		i=x;j=y;
		GridPoint2D=i+j*m_ContWorld->GridSize[0];
		zadd=k*m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1];
		list[count0]=GridPoint2D;
		;
		Map3D[GridPoint2D+zadd]=1;
		do
		{
			j=list[count0]/m_ContWorld->GridSize[0];
			i=list[count0]%m_ContWorld->GridSize[0];
			if(i>0)
				if((NI->GetDiffFloat(0,list[count0]-1+zadd)>0.0)
					&&(NI->GetDiffFloat(1,list[count0]-1+zadd)>0.0)
					&&(Map3D[list[count0]-1+zadd]<=0))
			{
				Map3D[list[count0]-1+zadd]=1;
				list[count1]=list[count0]-1;
				count1++;
			}
			if(i<m_ContWorld->GridSize[0]-1)
				if((NI->GetDiffFloat(ion,list[count0]-1+zadd)>0.0)
								&&(Map3D[list[count0]+1+zadd]<=0))
			{
				Map3D[list[count0]+1+zadd]=1;
				list[count1]=list[count0]+1;
				count1++;
			}
			if(j>0)
				if((NI->GetDiffFloat(ion,list[count0]-m_ContWorld->GridSize[0]+zadd)>0.0)
								&&(Map3D[list[count0]-m_ContWorld->GridSize[0]+zadd]<=0))
			{
				Map3D[list[count0]-m_ContWorld->GridSize[0]+zadd]=1;
				list[count1]=list[count0]-m_ContWorld->GridSize[0];
				count1++;
			}
			if(j<m_ContWorld->GridSize[1]-1)
				if((NI->GetDiffFloat(ion,list[count0]+m_ContWorld->GridSize[0]+zadd)>0.0)
								&&(Map3D[list[count0]+m_ContWorld->GridSize[0]+zadd]<=0))
			{
				Map3D[list[count0]+m_ContWorld->GridSize[0]+zadd]=1;
				list[count1]=list[count0]+m_ContWorld->GridSize[0];
				count1++;
			}
			count0++;
		}
		while(count0!=count1);
	}
	count0=0;
	count1=0;
	for(k=0;k<m_ContWorld->GridSize[2];k++)
		for(j=0;j<m_ContWorld->GridSize[1];j++)
			for(i=0;i<m_ContWorld->GridSize[0];i++)
	{
		GridPoint=i+j*m_ContWorld->GridSize[0]+k*m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1];
		if(Map3D[GridPoint]>0)
		{
			if(((i-x)%rstep==0)&&((j-y)%rstep==0)&&((k-z1)%rstep==0))
			{
				list[count0]=GridPoint;
				count1++;
			}
			count0++;
		}
	}
	
	FILE *out=fopen(FileName,"w");
	count0=0;
	count1=1;
	int TotalPnts=0;
	double ftmpx,ftmpy,ftmpz;
	int itmpx,itmpy,itmpz;
	int gsX=m_ContWorld->GridSize[0];
	int gsXY=m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1];
	for(k=0;k<m_ContWorld->GridSize[2];k++)
		for(j=0;j<m_ContWorld->GridSize[1];j++)
			for(i=0;i<m_ContWorld->GridSize[0];i++)
	{
		GridPoint=i+j*m_ContWorld->GridSize[0]+k*m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1];
		if(Map3D[GridPoint]>0)
		{
			if(((i-x)%rstep==0)&&((j-y)%rstep==0)&&((k-z1)%rstep==0))
			{
				list[count0]=GridPoint;
				//fprintf(out,"%d ",GridPoint);
				//itmpx=GridPoint%gsX;
				//itmpz=GridPoint/gsXY;
				//itmpy=(GridPoint-itmpz*gsXY)/gsX;
				itmpx=i-BoxSize[0];
				itmpy=j-BoxSize[1];
				itmpz=k-BoxSize[2];
				ftmpx=double(itmpx)/m_ContWorld->GridScale;
				ftmpy=double(itmpy)/m_ContWorld->GridScale;
				ftmpz=double(itmpz)/m_ContWorld->GridScale;
				fprintf(out,"ATOM  %5d %4s %c%3s%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f              \n",count1,"K  ",' ',"K  ",' ',count1,' ',ftmpx,ftmpy,ftmpz,1.0,1.0);
				//"ATOM      1  N   ALA A   1     -12.005 -12.238  25.958  1.00  0.00"
				VecPoints.push_back(GridPoint);
				count1++;
				if(count1>1000)count1=1;
				TotalPnts++;
			}
			count0++;
		}
	}
	fprintf(out,"TER\nEND\n");
	fclose(out);
	
	pnpPrint0("TotalPnts %d\n",TotalPnts);
	delete [] Map3D;
	delete [] list;
	return EXIT_SUCCESS;
}
int SIPPointsSearch::AlongZ(float fx,float fy,float fz1,float fz2,float frstep,const char *FileName)
{
	FILE *out=fopen(FileName,"w");
	int i,N;
	int count1=1;
	double ftmpx=fx,ftmpy=fy,ftmpz=fz1,dtmpz;
	dtmpz=fz2-fz1;
	N=dtmpz/frstep;
	for(i=0;i<=N;i++)
	{
		ftmpz=fz1+frstep*i;
		fprintf(out,"ATOM  %5d %4s %c%3s%c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f              \n",count1,"K  ",' ',"K  ",' ',count1,' ',ftmpx,ftmpy,ftmpz,1.0,1.0);
	}
	fprintf(out,"TER\nEND\n");
	fclose(out);
	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
PMFProcTools::PMFProcTools()
{
	
}
PMFProcTools::~PMFProcTools()
{
}
int PMFProcTools::SetMaskWherePMFNotEqVal(VectorIntField3D* Mask,int SetMask,VectorField3D* PMF,float NotEqVal)
{
	pnpPrint("SetMaskWherePMFNotEqVal\n");
	int iElm;
	unsigned int GrdPnt;
	unsigned int GS_X = PMF->GridSize[0];
	unsigned int GS_Y = PMF->GridSize[1];
	unsigned int GS_Z = PMF->GridSize[2];
	unsigned int GS_XY=GS_X*GS_Y;
	unsigned int GS_XYZ=GS_XY*GS_Z;
	for(iElm=0;iElm<PMF->Nelem;iElm++)
	{
		for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
		{
			if(PMF->V[iElm][GrdPnt]!=NotEqVal)
				Mask->V[iElm][GrdPnt]=SetMask;
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcTools::AddToMaskWherePMFNotEqVal(VectorIntField3D* Mask,int AddMask,VectorField3D* PMF,float NotEqVal)
{
	pnpPrint("SetMaskWherePMFNotEqVal\n");
	int iElm;
	unsigned int GrdPnt;
	unsigned int GS_X = PMF->GridSize[0];
	unsigned int GS_Y = PMF->GridSize[1];
	unsigned int GS_Z = PMF->GridSize[2];
	unsigned int GS_XY=GS_X*GS_Y;
	unsigned int GS_XYZ=GS_XY*GS_Z;
	for(iElm=0;iElm<PMF->Nelem;iElm++)
	{
		for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
		{
			if(PMF->V[iElm][GrdPnt]!=NotEqVal)
				Mask->V[iElm][GrdPnt]+=AddMask;
		}
	}
	return EXIT_SUCCESS;
}
int PMFProcTools::SmoothLevelOffThrLaplas(VectorField3D* PMF,VectorIntField3D* Mask,int iters,int iMaskWhereSmooth)
{
	pnpPrint("SmoothLevelOffThrLaplas iters:%d\n",iters);
	int ix,iy,iz,i1,BnW;
	float count;
	
	int iElm;
	unsigned int GrdPnt;
	unsigned int GS_X = PMF->GridSize[0];
	unsigned int GS_Y = PMF->GridSize[1];
	unsigned int GS_Z = PMF->GridSize[2];
	unsigned int GS_XY=GS_X*GS_Y;
	unsigned int GS_XYZ=GS_XY*GS_Z;
	
	for(iElm=0;iElm<PMF->Nelem;iElm++)
	{
		float *V=PMF->V[iElm];
		for(i1=0;i1<iters;i1++)
		{
			for(BnW=0;BnW<2;BnW++)
			{
				for(ix=1;ix<GS_X-1;ix++)
					for(iy=1;iy<GS_Y-1;iy++)
						for(iz=1;iz<GS_Z-1;iz++)
				{
					GrdPnt=ix+iy*GS_X+iz*GS_XY;
					if((ix+iy+iz)%2==BnW)
					{
						if(Mask->V[iElm][GrdPnt]==iMaskWhereSmooth)
						{
							V[GrdPnt]=(V[GrdPnt]*6.0+V[GrdPnt+1]+V[GrdPnt-1]+V[GrdPnt-GS_X]+V[GrdPnt+GS_X]+V[GrdPnt-GS_XY]+V[GrdPnt+GS_XY])/12.0;
						}
					}
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////


