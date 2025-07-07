//
// C++ Implementation: nernstplanksolver
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

#include "nernstplanksolver.h"
#include "tinyxml.h"
#include "pnpdebug.h"
#include "contworld.h"
#include "math.h"
#include "pnpconstants.h"
#include "mapio.h"
#include "pnputil.h"

#include <Python.h>

using namespace std;


NernstPlankSolver::NernstPlankSolver()
 : PnpsObject(), GenericSolver()
{
	InitZero();
}


NernstPlankSolver::~NernstPlankSolver()
{
	Clear();
}
int NernstPlankSolver::InitZero()
{
	PnpsObject::SetName("NernstPlankSolver");
	SolverStr.push_back("Auto");
	SolverStr.push_back("NodeIndexBased");
	SolverStr.push_back("ArrayDirect");
	SolverStr.push_back("PNPC");
	
	World=NULL;
	ConvergenceCheck=20;
	//CurrentCalc=2000;
	PMFWeight=-1.0;
	solver=0;
	writeStatQmob=false;
	bConcCorrection=false;
	PMFWeight=-1.0;
	int i,j;
  //Initial values of variables:
	verbose=true;
	UTMPSingle = NULL;
	TMPSingle  = NULL;
	UTMPDouble = NULL;
	TMPDouble  = NULL;
		
	QmobMask=NULL;
	writeStatQmob=false;
	solver=Auto;
	for(i=0;i<2;i++)
	{
		dixt[i] = NULL;
		for(j=0;j<6;j++)
		{
			dix[i][j] = NULL;
		}
	}
	for(j=0;j<6;j++)
	{
		dphi[j] = NULL;
	}
	bConcCorrection=false;
	//bDouble=false;
	CalcVolume=NULL;
	return EXIT_SUCCESS;
}
int NernstPlankSolver::Clear()
{
	int i,j,IType;
	
	DeleteCArray(CalcVolume);
	
	if(World!=NULL)
	{
		DeleteCVecArray(QmobMask,World->NIonsTypes);
		DeleteCVecArray(IndexNoSingular,World->NIonsTypes);
		DeleteCVecArray(IndexSingular,World->NIonsTypes);
		DeleteCArray(UTMPDouble);
		DeleteCArray(TMPDouble);
		DeleteCArray(UTMPSingle);
		DeleteCArray(TMPSingle);
		DeleteCVecArray(dixt,World->NIonsTypes);
		
		for(IType=0;IType<2;IType++)
		{
			for(j=0;j<6;j++)if(dix[IType][j]!=0)delete [] dix[IType][j];
			
		}
		for(j=0;j<6;j++)if(dphi[j]!=0)delete [] dphi[j];
	}
	return EXIT_SUCCESS;
}
int NernstPlankSolver::SaveXML(TiXmlElement* Elt, HaContext* p_ctxt )
{
	return EXIT_SUCCESS;
}
int NernstPlankSolver::LoadXML(const TiXmlElement* Elt, HaContext* p_ctxt )
{
	if(Elt==NULL)
	{
		fprintf(stderr,"ERROR: Can't find <%s>\n",PnpsObject::GetCStrName());
		return EXIT_FAILURE;
	}
	if(strcmp(PnpsObject::GetCStrName(),Elt->Value()))
	{
		fprintf(stderr,"ERROR: Wrong XML Element %s, expecting %s\n",Elt->Value(),PnpsObject::GetCStrName());
		return EXIT_FAILURE;
	}
	Clear();
	int i,gridPoint;
	//Read Primary Parameters
	Elt->GetIntAttribute("MaxIterations",&MaxIterations);
	Elt->GetFloatAttribute("Convergence",&Convergence);
	Elt->GetFloatAttribute("Relaxation",&Relaxation);  
	
	if(Elt->GetIntAttribute("ConvergenceCheck",&ConvergenceCheck)!=EXIT_SUCCESS)ConvergenceCheck=20;
	//if(Elt->GetIntAttribute("CurrentCalc",&CurrentCalc)!=EXIT_SUCCESS)CurrentCalc=2000;
	
	if(Elt->GetStdStrIndex("Solver",&solver,SolverStr)!=EXIT_SUCCESS)solver=0;
	
	if(Elt->GetBoolAttribute("writeStatQmob",&writeStatQmob)!=EXIT_SUCCESS)
		writeStatQmob=false;
	if(Elt->GetStdStrAttribute("QmobMaskFile",&QmobMaskFile)!=EXIT_SUCCESS)
	{
		QmobMaskFile="";
		writeStatQmob=false;
	}
	if(Elt->GetBoolAttribute("bConcCorrection",&bConcCorrection)!=EXIT_SUCCESS)
		bConcCorrection=false;
	
	
	if(Elt->GetBoolAttribute("Verbose",&verbose)!=EXIT_SUCCESS)verbose=true;
	
	if(Elt->GetFloatAttribute("PMFWeight",&PMFWeight)!=EXIT_SUCCESS)PMFWeight=-1.0;
	else if(PMFWeight>1.0)PMFWeight=-1.0;
	
	//if(Elt->GetBoolAttribute("bDouble",&bDouble)!=EXIT_SUCCESS)
	//	bDouble=false;
	
	ShowParameters();
	return EXIT_SUCCESS;
}
int NernstPlankSolver::LoadParamFromPyDict(PyObject *dict)
{
	if(dict==NULL)
	{
		fprintf(stderr,"ERROR: Can't find <%s>\n",PnpsObject::GetCStrName());
		return EXIT_FAILURE;
	}
	Clear();
	int i,gridPoint;
	//Read Primary Parameters
	haPyDict_GetItemAsInt(dict,"MaxIterations",&MaxIterations);
	haPyDict_GetItemAsFloat(dict,"Relaxation",&Relaxation);
	haPyDict_GetItemAsFloat(dict,"Convergence",&Convergence);
	haPyDict_GetItemAsFloat(dict,"Tolerance",&Convergence);
	haPyDict_GetItemAsInt(dict,"ConvergenceCheck",&ConvergenceCheck);
	haPyDict_GetItemAsInt(dict,"solver",&solver);
	haPyDict_GetItemAsBool(dict,"Verbose",&verbose);

	haPyDict_GetItemAsBool(dict,"writeStatQmob",&writeStatQmob);
	haPyDict_GetItemAsBool(dict,"bConcCorrection",&bConcCorrection);
	haPyDict_GetItemAsString(dict,"QmobMaskFile",&QmobMaskFile);
	haPyDict_GetItemAsFloat(dict,"PMFWeight",&PMFWeight);

	
	
	//if(Elt->GetBoolAttribute("bDouble",&bDouble)!=EXIT_SUCCESS)
	//	bDouble=false;
	
	ShowParameters();
	return EXIT_SUCCESS;
}
int NernstPlankSolver::SetContWorld(ContWorld* _world)
{
	int i,j;
	World=_world;
	if(World->C==NULL)
	{
		World->SetInitConcentrationFromNIndexing();
	}
	else if(World->C[0]==NULL)
	{
		World->SetInitConcentrationFromNIndexing();
	}
	if(this->QmobMaskFile!="")
	{
		pnpPrint("Loading QmobMaskFile from %s\n",QmobMaskFile.c_str());
			
		if(QmobMask!=NULL)delete [] QmobMask;
			
		VectorField3D* VField3DContr=new VectorField3D();
#ifndef MPI_PARALLEL
		VField3DContr->ReadFromFileAddPBC(QmobMaskFile.c_str(), 1.0,World->PBC[0],World->PBC[1],World->PBC[2]);
#else
		if(pnpsapp->GetNumberOfGroups()==pnpsapp->GetTotalProc())
			VField3DContr->ReadFromFileAddPBC(QmobMaskFile.c_str(), 1.0,World->PBC[0],World->PBC[1],World->PBC[2]);
		//if(pnpsapp->GetNumProcsInMyGroup()==pnpsapp->GetTotalProc())
		else
			VField3DContr->ReadFromFileAddPBCandSplitSDMPI4FDaZ(QmobMaskFile.c_str(), 1.0,World->PBC,World->MyGlobalZ0,World->MyGlobalZ1);
#endif
		QmobMask=new int*[World->NIonsTypes];
		for(j=0;j<World->NIonsTypes;j++)
		{
			QmobMask[j]=new int[World->GS_XYZ];
			for(i=0;i<World->GS_XYZ;i++)
			{
				if(VField3DContr->V[j][i]>0.0)
					QmobMask[j][i]=1;
				else
					QmobMask[j][i]=0;
			}
		}
		
		
		delete VField3DContr;
		
		if(writeStatQmob)
		{
			VectorField3D* VField3DContr2=new VectorField3D(World->GridSize,World->GridScale,World->NIonsTypes);
			for(j=0;j<VField3DContr2->GetNelem();j++)
			{
				for(i=0;i<World->GS_XYZ;i++)
				{
					if(World->NIndexing->GetDiffFloat(j,i)>0.0&&QmobMask[j][i]==1)
						VField3DContr2->V[j][i]=1.0;
					else
						VField3DContr2->V[j][i]=0.0;
				}
			}
			VField3DContr2->WriteToFile("StatQmob.bin",1.0,1);
			delete VField3DContr2;
		}
	}
	return EXIT_SUCCESS;
}
int NernstPlankSolver::ShowParameters()
{
	fprintf(stdout,"\nParameters of Nernst-Plank solver\n");
	fprintf(stdout,"    MaxIterations:.................... %d\n", MaxIterations);
	fprintf(stdout,"    Convergence:...................... %g\n", Convergence);
	fprintf(stdout,"    Relaxation:....................... %g\n", Relaxation);
	return EXIT_SUCCESS;
}
int NernstPlankSolver::InitSolver()
{
//	 int i;
//	 for(i=0;i<World->NIonsTypes;i++){
//		 if(UTMP[i] == NULL)if(!(UTMP[i] = new float[World->GridSizeXYZ])){
//			 cerr << "ERROR 104: No memory available\n";
//			 exit(104);
//		 }
//		 if(TMP[i] == NULL)if(!(TMP[i] = new float[World->GridSizeXYZ])){
//			 cerr << "ERROR 104: No memory available\n";
//			 exit(104);
//		 }
//	 }
	
	if(solver==Auto)
	{
		if(World->NIndexing==NULL)return InitSolverW();
		else if(World->D!=NULL)return InitSolverDflex();
		//else if(World->D!=NULL)return InitSolverW();
		else return InitSolverD();
	}
	else
	{
		if(solver==NodeIndexBased)
			return InitSolverD();
		else if(solver==ArrayDirect)
			return InitSolverAD();
		else if(solver==PNPC)
			return InitSolverW();
		else return InitSolverD();
	}
	
}
int NernstPlankSolver::InitSolverW()
{
	DbgPrint2("NernstPlankSolver::InitSolverW()\n");
	//int nernstPlanckSolverSetup(NernstPlanckSolverData ** nernstPlanckSolverData, int maxIterations, float convergence, float relaxation, int expInterpolation, Output * output) {
	int i,j,k;
	long gridPoint;
	long jgrid,kgrid;
	float * diffusion;
	int	gridSizeX;
	int	gridSizeY;
	int	gridSizeZ;
	long gridSizeXY;
	long count;
	long borderCount;
	int expInterpolation=0;
	/*assert(output!=NULL);
	assert(output->diffusionMap!=NULL);
	assert(output->staticChargeMap!=NULL);
	assert(output->positiveDynamicChargeMap!=NULL);
	assert(output->negativeDynamicChargeMap!=NULL);*/
	if(World->D==NULL)
	{
		World->D=new float*[2];
		World->D[0]=World->NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion0);
		World->D[1]=NULL;
	}
	//WriteMapGZ("diff.gz",World->D[0],World->GridSize,1.0);
	diffusion = World->D[0];
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;

	//nernstPlanckSolverData = (NernstPlanckSolverDataW*)malloc(sizeof(NernstPlanckSolverDataW));
	nernstPlanckSolverData = new NernstPlanckSolverDataW;
	assert(nernstPlanckSolverData!=NULL);

	nernstPlanckSolverData->maxIterations = MaxIterations;
	nernstPlanckSolverData->convergence = Convergence;
	nernstPlanckSolverData->relaxation = Relaxation;
	nernstPlanckSolverData->expInterpolation = expInterpolation;

	for(k=1;k<gridSizeZ-1;k++) {
		kgrid = k*gridSizeXY;
		for(j=1;j<gridSizeY-1;j++) {
			jgrid = kgrid+j*gridSizeX;
			for(i=1;i<gridSizeX-1;i++) {
				gridPoint = jgrid+i;
				if(diffusion[gridPoint]!=0 && diffusion[gridPoint+1]==0 && diffusion[gridPoint-1]==0 && diffusion[gridPoint+gridSizeX]==0 && diffusion[gridPoint-gridSizeX]==0 && diffusion[gridPoint+gridSizeXY]==0 && diffusion[gridPoint-gridSizeXY]==0)
				{
					diffusion[gridPoint]=0;
					DbgPrint2("Bed Diffusional point at [%d] have removed it.",gridPoint);
				}
			}
		}
	}

	count = 0;
	
	for(k=1;k<gridSizeZ-1;k++) {
		kgrid = k*gridSizeXY;
		for(j=1;j<gridSizeY-1;j++) {
			jgrid = kgrid+j*gridSizeX;
			for(i=1;i<gridSizeX-1;i++) {
				gridPoint = jgrid+i;
				if(diffusion[gridPoint-1]>0 && diffusion[gridPoint]<=0)
					count++;
			}
			if(diffusion[jgrid+gridSizeX-2]>0)
				count++;
		}
	}
	DbgPrint0("nernstPlanckSolver count: %d\n",count);
	nernstPlanckSolverData->start = new long[count];
	assert(nernstPlanckSolverData->start!=NULL);
	nernstPlanckSolverData->end = new long[count+1];
	assert(nernstPlanckSolverData->end!=NULL);

	count=0;
	for(k=1;k<gridSizeZ-1;k++) {
		kgrid = k*gridSizeXY;
		for(j=1;j<gridSizeY-1;j++) {
			jgrid = kgrid+j*gridSizeX;
			if(diffusion[jgrid]>0)
				nernstPlanckSolverData->start[count]=jgrid+1;
			for(i=1;i<gridSizeX-1;i++) {
				gridPoint = jgrid+i;
				if(diffusion[gridPoint-1]<=0 && diffusion[gridPoint]>0)
					nernstPlanckSolverData->start[count] = gridPoint;
				if(diffusion[gridPoint-1]>0 && diffusion[gridPoint]<=0) {
					nernstPlanckSolverData->end[count]=gridPoint;
					count++;
				}
			}
			if(diffusion[jgrid+gridSizeX-2]>0) {
				nernstPlanckSolverData->end[count] = jgrid+gridSizeX-1;
				count++;
			}
		}
	}
	nernstPlanckSolverData->end[count] = gridSizeZ*gridSizeXY;

	count = 0;
	borderCount=0;
	while(nernstPlanckSolverData->end[count]<gridSizeZ*gridSizeXY) {
		for(gridPoint=nernstPlanckSolverData->start[count];gridPoint<nernstPlanckSolverData->end[count];gridPoint++) {
			if(diffusion[gridPoint]!=diffusion[gridPoint+1] || diffusion[gridPoint]!=diffusion[gridPoint-1] || diffusion[gridPoint]!=diffusion[gridPoint+gridSizeX] || diffusion[gridPoint]!=diffusion[gridPoint-gridSizeX] || diffusion[gridPoint]!=diffusion[gridPoint+gridSizeXY] || diffusion[gridPoint]!=diffusion[gridPoint-gridSizeXY])
				borderCount++;
		}
		count++;
	}
	DbgPrint0("nernstPlanckSolver borderCount: %d\n",borderCount);
	
	nernstPlanckSolverData->borderPoints = new long[borderCount+1];
	assert(nernstPlanckSolverData->borderPoints!=NULL);
	nernstPlanckSolverData->dixt = new float[borderCount];
	assert(nernstPlanckSolverData->dixt!=NULL);
	for(i=0;i<6;i++) {
		nernstPlanckSolverData->dix[i] = new float[borderCount];
		assert(nernstPlanckSolverData->dix[i]!=NULL);
	}

	count = 0;
	borderCount=0;
	while(nernstPlanckSolverData->end[count]<gridSizeZ*gridSizeXY) {
		for(gridPoint=nernstPlanckSolverData->start[count];gridPoint<nernstPlanckSolverData->end[count];gridPoint++) {
			if(diffusion[gridPoint]!=diffusion[gridPoint+1] || diffusion[gridPoint]!=diffusion[gridPoint-1] || diffusion[gridPoint]!=diffusion[gridPoint+gridSizeX] || diffusion[gridPoint]!=diffusion[gridPoint-gridSizeX] || diffusion[gridPoint]!=diffusion[gridPoint+gridSizeXY] || diffusion[gridPoint]!=diffusion[gridPoint-gridSizeXY]) {
				nernstPlanckSolverData->borderPoints[borderCount] = gridPoint;
				nernstPlanckSolverData->dix[0][borderCount] = diffusion[gridPoint+1]>0?0.5*(diffusion[gridPoint]+diffusion[gridPoint+1]):0;
				nernstPlanckSolverData->dix[1][borderCount] = diffusion[gridPoint-1]>0?0.5*(diffusion[gridPoint]+diffusion[gridPoint-1]):0;
				nernstPlanckSolverData->dix[2][borderCount] = diffusion[gridPoint+gridSizeX]>0?0.5*(diffusion[gridPoint]+diffusion[gridPoint+gridSizeX]):0;
				nernstPlanckSolverData->dix[3][borderCount] = diffusion[gridPoint-gridSizeX]>0?0.5*(diffusion[gridPoint]+diffusion[gridPoint-gridSizeX]):0;
				nernstPlanckSolverData->dix[4][borderCount] = diffusion[gridPoint+gridSizeXY]>0?0.5*(diffusion[gridPoint]+diffusion[gridPoint+gridSizeXY]):0;
				nernstPlanckSolverData->dix[5][borderCount] = diffusion[gridPoint-gridSizeXY]>0?0.5*(diffusion[gridPoint]+diffusion[gridPoint-gridSizeXY]):0;
				nernstPlanckSolverData->dixt[borderCount] = nernstPlanckSolverData->dix[0][borderCount]+nernstPlanckSolverData->dix[1][borderCount]+nernstPlanckSolverData->dix[2][borderCount]+nernstPlanckSolverData->dix[3][borderCount]+nernstPlanckSolverData->dix[4][borderCount]+nernstPlanckSolverData->dix[5][borderCount];
				borderCount++;
			}
		}
		count++;
	}
	nernstPlanckSolverData->borderPoints[borderCount] = gridSizeZ*gridSizeXY;
	return EXIT_SUCCESS;
}
int NernstPlankSolver::InitSolverD()
{
	DbgPrint2("NernstPlankSolver::InitSolverD()\n");
	int i,j,k, i1, IType;
	int ix, iy, iz;
	int GrdPnt;
	int jgrid,kgrid;
	//float * diffusion;
	//float **D;
	//int NIonsTypes;
	int GS_X = World->GridSize[0];
	int GS_Y = World->GridSize[1];
	int GS_Z = World->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	//int *tmp;
	int itmp1,itmp2;
	//int BlackOrWhite;
	
	NodeIndexing* NIndexing=World->NIndexing;
	NodeIndex* NIndex=NIndexing->NIndex;
	
	NodeIndex DiffBoarderMask=NodeIndexing::DiffIon0BoarderMask;
	NodeIndex BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;
	
//	 UTMP[0]=0;
//	 UTMP[1]=0;
//	 TMP[0]=0;
//	 TMP[1]=0;
	//for(i=0;i<2;i++)for(j=0;j<6;j++)dphi[i][j]=NULL;
	MaxChange=0.0;
	
	
	if(!World->C[0]){
		cerr << "ERROR 311: C[0] Map Not init\n";
		exit(104);
	}
	if(!World->C[1]){
		cerr << "ERROR 311: C[1] Map Not init\n";
		exit(104);
	}
	if(!World->NIndexing){
		cerr << "ERROR 311: World->NIndexing Map Not init\n";
		exit(104);
	}
	if(!World->Potential){
		cerr << "WARNING 311: Potential Map Not init, will init and set to zero values\n";
		World->Potential=new float[GS_XYZ];
		for(i=0;i<GS_XYZ;i++)World->Potential[i]=0.0;
	}
	
	int DiffZeroInd=NIndexing->GetDiffZeroInd();

	NIndexing->RemoveBadDiffusionPoints();
//#endif
	for(IType=0;IType<World->NIonsTypes;IType++)
		for(i=0;i<GS_XYZ;i++)
			if(NIndexing->GetDiffFloat(IType,i)==0.0)World->C[IType][i]=0.0f;
	
	
	for(IType=0;IType<World->NIonsTypes;IType++)
	{
		DiffBoarderMask=NIndexing->DiffBoarderMask[IType];
		//Calculate number of nodes with singularity and not and exclude nodes with D=0
		
		NoSingularNum[IType][0] = 0;
		NoSingularNum[IType][1] = 0;
		NoSingularNum[IType][2] = 0;
		SingularNum[IType][0] = 0;
		SingularNum[IType][1] = 0;
		SingularNum[IType][2] = 0;
		if(QmobMask==NULL)
		{
			for(k=1;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						GrdPnt = jgrid+i;
						if(NIndexing->GetDiffFloat(IType,GrdPnt)>0.0)
						{
							if(NIndex[GrdPnt]&DiffBoarderMask)
							{
								SingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
								SingularNum[IType][2]++;
							}
							else{
								NoSingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
								NoSingularNum[IType][2]++;
							}
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
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						GrdPnt = jgrid+i;
						if(QmobMask[IType][GrdPnt]>0)
							if(NIndexing->GetDiffFloat(IType,GrdPnt)>0.0)
						{
							if(NIndex[GrdPnt]&DiffBoarderMask)
							{
								SingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
								SingularNum[IType][2]++;
							}
							else{
								NoSingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
								NoSingularNum[IType][2]++;
							}
						}
					}
				}
			}
		}
			
		fprintf(stdout,"Nernst-Plank Solver:\n");
		fprintf(stdout,"		IonsQ:........................ %g\n", World->IonsQ[IType]);
		fprintf(stdout,"		SingularNum:...................... %d\n", SingularNum[IType][2]);
		fprintf(stdout,"		NoSingularNum:.................... %d\n", NoSingularNum[IType][2]);
		fprintf(stdout,"		IndexSingularParts:............... [%d,%d,%d]\n",
						SingularNum[IType][0],SingularNum[IType][1],SingularNum[IType][2]);
		fprintf(stdout,"		IndexNoSingularParts:............. [%d,%d,%d]\n",NoSingularNum[IType][0],NoSingularNum[IType][1],NoSingularNum[IType][2]);
		
		if(!(IndexSingular[IType] = new int[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		if(!(IndexNoSingular[IType] = new int[NoSingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		//cout<<"NP3sleep 10 sec\n";
		//sleep(10);
		for(i=0;i<6;i++)if(!(dix[IType][i] = new float[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		if(!(dixt[IType] = new float[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		int iSingular=0,iNoSingular=0;
		int iSingular2=SingularNum[IType][1],iNoSingular2=NoSingularNum[IType][1];
		if(QmobMask==NULL)
		{
			for(k=1;k<GS_Z-1;k++) 
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++) 
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++) 
					{
						GrdPnt = jgrid+i;
						if(NIndexing->GetDiffFloat(IType,GrdPnt)>0.0)
						{
							if(NIndex[GrdPnt]&DiffBoarderMask)
							{
								if(NIndex[GrdPnt]&BlackAndWhiteMask)
								{
									IndexSingular[IType][iSingular]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iSingular++;
								}
								else
								{
									IndexSingular[IType][iSingular2]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular2, GrdPnt);
									iSingular2++;
								}
							}
							else
							{
								if(NIndex[GrdPnt]&BlackAndWhiteMask)
								{
									IndexNoSingular[IType][iNoSingular]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular++;
								}
								else
								{
									IndexNoSingular[IType][iNoSingular2]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular2++;
								}
							}
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
				for(j=1;j<GS_Y-1;j++) 
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++) 
					{
						GrdPnt = jgrid+i;
						if(NIndexing->GetDiffFloat(IType,GrdPnt)>0.0)
						{
							
							if(QmobMask[IType][GrdPnt]>0)
								if(NIndex[GrdPnt]&DiffBoarderMask)
							{
								if(NIndex[GrdPnt]&BlackAndWhiteMask)
								{
									IndexSingular[IType][iSingular]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iSingular++;
								}
								else
								{
									IndexSingular[IType][iSingular2]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular2, GrdPnt);
									iSingular2++;
								}
							}
							else
							{
								if(NIndex[GrdPnt]&BlackAndWhiteMask)
								{
									IndexNoSingular[IType][iNoSingular]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular++;
								}
								else
								{
									IndexNoSingular[IType][iNoSingular2]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular2++;
								}
							}
						}
					}
				}
			}
		}
		
		for(i=0;i<SingularNum[IType][2];i++)
		{
			GrdPnt=IndexSingular[IType][i];
			dix[IType][0][i] = NIndexing->GetDiffFloat(IType,GrdPnt+1)>0.0?
					0.5*(NIndexing->GetDiffFloat(IType,GrdPnt)+NIndexing->GetDiffFloat(IType,GrdPnt+1)):0.0;
			dix[IType][1][i] = NIndexing->GetDiffFloat(IType,GrdPnt-1)>0.0?
					0.5*(NIndexing->GetDiffFloat(IType,GrdPnt)+NIndexing->GetDiffFloat(IType,GrdPnt-1)):0.0;
			dix[IType][2][i] = NIndexing->GetDiffFloat(IType,GrdPnt+GS_X)>0.0?
					0.5*(NIndexing->GetDiffFloat(IType,GrdPnt)+NIndexing->GetDiffFloat(IType,GrdPnt+GS_X)):0.0;
			dix[IType][3][i] = NIndexing->GetDiffFloat(IType,GrdPnt-GS_X)>0.0?
					0.5*(NIndexing->GetDiffFloat(IType,GrdPnt)+NIndexing->GetDiffFloat(IType,GrdPnt-GS_X)):0.0;
			dix[IType][4][i] = NIndexing->GetDiffFloat(IType,GrdPnt+GS_XY)>0.0?
					0.5*(NIndexing->GetDiffFloat(IType,GrdPnt)+NIndexing->GetDiffFloat(IType,GrdPnt+GS_XY)):0.0;
			dix[IType][5][i] = NIndexing->GetDiffFloat(IType,GrdPnt-GS_XY)>0.0?
					0.5*(NIndexing->GetDiffFloat(IType,GrdPnt)+NIndexing->GetDiffFloat(IType,GrdPnt-GS_XY)):0.0;
			dixt[IType][i] = dix[IType][0][i]+dix[IType][1][i]+dix[IType][2][i]+dix[IType][3][i]+dix[IType][4][i]+dix[IType][5][i];
		}
	}
	
	if(World->PMF!=NULL)fprintf(stdout,"Calculate NP with PMF\n");
	/*for(i=0;i<World->GridSizeXYZ;i++)tmp[i]=0;
	for(i=0;i<SingularNum;i++)tmp[IndexSingular[i]]=1;
	for(i=0;i<NoSingularNum;i++)tmp[IndexNoSingular[i]]=2;
	World->writeMapInt("tmp.gz",tmp,World->GridSizeXYZ,1);
	delete [] tmp;*/
	
	//SAVE memory for children! Delete excess of D[0]
	//if(World->D[0]!=0)delete [] World->D[0];
	//cout<<"sleep 10 sec\n";
	//sleep(10);
	return 1;
}
int NernstPlankSolver::InitSolverDflex()
{
	DbgPrint2("NernstPlankSolver::InitSolverDflex()\n");
	int i,j,k, i1, IType;
	int ix, iy, iz;
	int GrdPnt;
	int jgrid,kgrid;
	//float * diffusion;
	//float **D;
	//int NIonsTypes;
	int GS_X = World->GridSize[0];
	int GS_Y = World->GridSize[1];
	int GS_Z = World->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	//int *tmp;
	int itmp1,itmp2;
	float **D=World->D;
	//int BlackOrWhite;
	
	//If volume there run calculation is not set, make a fake one
	
	
	NodeIndexing* NIndexing=World->NIndexing;
	NodeIndex* NIndex=NIndexing->NIndex;
	
	NodeIndex DiffBoarderMask=NodeIndexing::DiffIon0BoarderMask;
	NodeIndex BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;

//	 UTMP[0]=0;
//	 UTMP[1]=0;
//	 TMP[0]=0;
//	 TMP[1]=0;
	//for(i=0;i<2;i++)for(j=0;j<6;j++)dphi[i][j]=NULL;
	MaxChange=0.0;
	
	
	if(!World->C[0]){
		cerr << "ERROR 311: C[0] Map Not init\n";
		exit(104);
	}
	if(!World->C[1]){
		cerr << "ERROR 311: C[1] Map Not init\n";
		exit(104);
	}
	if(!World->NIndexing){
		cerr << "ERROR 311: World->NIndexing Map Not init\n";
		exit(104);
	}
	if(!World->Potential){
		cerr << "WARNING 311: Potential Map Not init, will init and set to zero values\n";
		World->Potential=new float[GS_XYZ];
		for(i=0;i<GS_XYZ;i++)World->Potential[i]=0.0;
	}
	

	/*if(!World->PotentialTMP){
	if(!(World->PotentialTMP = new float[World->GridSizeXYZ])){
	cerr << "ERROR 104: No memory available\n";
	exit(104);
}
}*/
//cout<<"NP2sleep 10 sec\n";
	//sleep(10);
	
	//D=World->D;
	for(IType=0;IType<World->NIonsTypes;IType++)
	{
		for(k=1;k<GS_Z-1;k++)
		{
			kgrid = k*GS_XY;
			for(j=1;j<GS_Y-1;j++)
			{
				jgrid = kgrid+j*GS_X;
				for(i=1;i<GS_X-1;i++)
				{
					GrdPnt = jgrid+i;
					if(D[IType][GrdPnt]!=0.0 && D[IType][GrdPnt+1]==0.0 && D[IType][GrdPnt-1]==0.0 && D[IType][GrdPnt+GS_X]==0.0 && D[IType][GrdPnt-GS_X]==0.0 && D[IType][GrdPnt+GS_XY]==0.0 &&D[IType][GrdPnt-GS_XY]==0.0)
					{
						D[IType][GrdPnt]=0.0;
						DbgPrint2("Bad Diffusional point at [%d] have removed it.",GrdPnt);
					}
				}
			}
		}
	}

	for(IType=0;IType<World->NIonsTypes;IType++)
		for(i=0;i<GS_XYZ;i++)
			if(World->D[IType][i]==0.0)World->C[IType][i]=0.0f;
	
	bool bCalcVolume;
	bCalcVolume=true;
	
	for(IType=0;IType<World->NIonsTypes;IType++)
	{
		//DiffBoarderMask=NIndexing->DiffBoarderMask[IType];
		//Calculate number of nodes with singularity and not and exclude nodes with D=0
		
		NoSingularNum[IType][0] = 0;
		NoSingularNum[IType][1] = 0;
		NoSingularNum[IType][2] = 0;
		SingularNum[IType][0] = 0;
		SingularNum[IType][1] = 0;
		SingularNum[IType][2] = 0;
		if(QmobMask==NULL)
		{
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
							if(D[IType][GrdPnt]>0.0){
								if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
								{
									SingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
									SingularNum[IType][2]++;
								}
								else{
									NoSingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
									NoSingularNum[IType][2]++;
								}
							}
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
							if(QmobMask[IType][GrdPnt]>0&&D[IType][GrdPnt]>0.0)
							{
								if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
								{
									SingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
									SingularNum[IType][2]++;
								}
								else{
									NoSingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
									NoSingularNum[IType][2]++;
								}
							}
						}
					}
				}
			}
		}
			
		fprintf(stdout,"Nernst-Plank Solver:\n");
		fprintf(stdout,"		IonsQ:........................ %g\n", World->IonsQ[IType]);
		fprintf(stdout,"		SingularNum:...................... %d\n", SingularNum[IType][2]);
		fprintf(stdout,"		NoSingularNum:.................... %d\n", NoSingularNum[IType][2]);
		fprintf(stdout,"		IndexSingularParts:............... [%d,%d,%d]\n",
						SingularNum[IType][0],SingularNum[IType][1],SingularNum[IType][2]);
		fprintf(stdout,"		IndexNoSingularParts:............. [%d,%d,%d]\n",NoSingularNum[IType][0],NoSingularNum[IType][1],NoSingularNum[IType][2]);
		
		if(!(IndexSingular[IType] = new int[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		if(!(IndexNoSingular[IType] = new int[NoSingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		//cout<<"NP3sleep 10 sec\n";
		//sleep(10);
		for(i=0;i<6;i++)if(!(dix[IType][i] = new float[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		if(!(dixt[IType] = new float[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		int iSingular=0,iNoSingular=0;
		int iSingular2=SingularNum[IType][1],iNoSingular2=NoSingularNum[IType][1];
		if(QmobMask==NULL)
		{
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
							if(D[IType][GrdPnt]>0.0)
							{
								if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
								{
									if(NIndex[GrdPnt]&BlackAndWhiteMask)
									{
										IndexSingular[IType][iSingular]=GrdPnt;
										//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
										iSingular++;
									}
									else
									{
										IndexSingular[IType][iSingular2]=GrdPnt;
										//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular2, GrdPnt);
										iSingular2++;
									}
								}
								else
								{
									if(NIndex[GrdPnt]&BlackAndWhiteMask)
									{
										IndexNoSingular[IType][iNoSingular]=GrdPnt;
										//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
										iNoSingular++;
									}
									else
									{
										IndexNoSingular[IType][iNoSingular2]=GrdPnt;
										//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
										iNoSingular2++;
									}
								}
							}
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
							if(D[IType][GrdPnt]>0.0&&QmobMask[IType][GrdPnt]>0)
							{
								if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
													D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
								{
									if(NIndex[GrdPnt]&BlackAndWhiteMask)
									{
										IndexSingular[IType][iSingular]=GrdPnt;
										//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
										iSingular++;
									}
									else
									{
										IndexSingular[IType][iSingular2]=GrdPnt;
										//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular2, GrdPnt);
										iSingular2++;
									}
								}
								else
								{
									if(NIndex[GrdPnt]&BlackAndWhiteMask)
									{
										IndexNoSingular[IType][iNoSingular]=GrdPnt;
										//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
										iNoSingular++;
									}
									else
									{
										IndexNoSingular[IType][iNoSingular2]=GrdPnt;
										//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
										iNoSingular2++;
									}
								}
							}
						}
					}
				}
			}
		}
		
		for(i=0;i<SingularNum[IType][2];i++)
		{
			GrdPnt=IndexSingular[IType][i];
			dix[IType][0][i] = NIndexing->GetDiffFloat(IType,GrdPnt+1)>0.0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt+1]):0.0;
			dix[IType][1][i] = NIndexing->GetDiffFloat(IType,GrdPnt-1)>0.0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt-1]):0.0;
			dix[IType][2][i] = NIndexing->GetDiffFloat(IType,GrdPnt+GS_X)>0.0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt+GS_X]):0.0;
			dix[IType][3][i] = NIndexing->GetDiffFloat(IType,GrdPnt-GS_X)>0.0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt-GS_X]):0.0;
			dix[IType][4][i] = NIndexing->GetDiffFloat(IType,GrdPnt+GS_XY)>0.0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt+GS_XY]):0.0;
			dix[IType][5][i] = NIndexing->GetDiffFloat(IType,GrdPnt-GS_XY)>0.0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt-GS_XY]):0.0;
			dixt[IType][i] = dix[IType][0][i]+dix[IType][1][i]+dix[IType][2][i]+dix[IType][3][i]+dix[IType][4][i]+dix[IType][5][i];
		}
	}
	
	if(World->PMF!=NULL)fprintf(stdout,"Calculate NP with PMF\n");
	/*for(i=0;i<World->GridSizeXYZ;i++)tmp[i]=0;
	for(i=0;i<SingularNum;i++)tmp[IndexSingular[i]]=1;
	for(i=0;i<NoSingularNum;i++)tmp[IndexNoSingular[i]]=2;
	World->writeMapInt("tmp.gz",tmp,World->GridSizeXYZ,1);
	delete [] tmp;*/
	
	//SAVE memory for children! Delete excess of D[0]
	//if(World->D[0]!=0)delete [] World->D[0];
	//cout<<"sleep 10 sec\n";
	//sleep(10);
	return 1;
}
/////////////////////////////////////////////////////////////////////////////
int NernstPlankSolver::InitSolverAD()
{
	DbgPrint2("NernstPlankSolver::InitSolverAD()\n");
	int i,j,k, i1, IType;
	int ix, iy, iz;
	int GrdPnt;
	int jgrid,kgrid;
	//float * diffusion;
	//float **D;
	//int NIonsTypes;
	int GS_X = World->GridSize[0];
	int GS_Y = World->GridSize[1];
	int GS_Z = World->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	//int *tmp;
	int itmp1,itmp2;
	
	int BlackOrWhite;
	
	//for(i=0;i<2;i++)for(j=0;j<6;j++)dphi[i][j]=NULL;
	MaxChange=0.0;
	
	if(World->D==NULL)
	{
		if(World->NIndexing!=NULL)
		{
			pnpPrint0("WARNING: Will convert Diffusion from NIndexing\n");
			World->D=new float*[World->NIonsTypes];
			if(World->NIonsTypes>=1)
				World->D[0]=World->NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion0);
			if(World->NIonsTypes>=2)
				World->D[1]=World->NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion1);
			if(World->NIonsTypes>=3)
				World->D[2]=World->NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion2);
			if(World->NIonsTypes>=4)
				World->D[3]=World->NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion3);
		}
		else
		{
			fprintf(stderr, "ERR:		Diffusion Map NOT yet initialize\n");
		}
	}
	float **D=World->D;
	float **C=World->C;
	if(!World->C[0]){
		cerr << "ERROR 311: C[0] Map Not init\n";
		exit(104);
	}
	if(!World->C[1]){
		cerr << "ERROR 311: C[1] Map Not init\n";
		exit(104);
	}
	if(!World->Potential){
		cerr << "WARNING 311: Potential Map Not init, will init and set to zero values\n";
		World->Potential=new float[GS_XYZ];
		for(i=0;i<GS_XYZ;i++)World->Potential[i]=0.0;
	}
	
	for(IType=0;IType<World->NIonsTypes;IType++)
	{
		for(k=1;k<GS_Z-1;k++)
		{
			kgrid = k*GS_XY;
			for(j=1;j<GS_Y-1;j++)
			{
				jgrid = kgrid+j*GS_X;
				for(i=1;i<GS_X-1;i++)
				{
					GrdPnt = jgrid+i;
					if(D[IType][GrdPnt]!=0.0 && D[IType][GrdPnt+1]==0.0 && D[IType][GrdPnt-1]==0.0 && D[IType][GrdPnt+GS_X]==0.0 && D[IType][GrdPnt-GS_X]==0.0 && D[IType][GrdPnt+GS_XY]==0.0 &&D[IType][GrdPnt-GS_XY]==0.0)
					{
						D[IType][GrdPnt]=0.0;
						C[IType][GrdPnt]=0.0;
						DbgPrint2("Bed Diffusional point at [%d] have removed it.",GrdPnt);
					}
				}
			}
		}
	}
	for(IType=0;IType<World->NIonsTypes;IType++)
		for(i=0;i<GS_XYZ;i++)
			if(D[IType][GrdPnt]==0.0)C[IType][i]=0.0f;
	
	
	for(IType=0;IType<World->NIonsTypes;IType++)
	{
		//DiffBoarderMask=NIndexing->DiffBoarderMask[IType];
		//Calculate number of nodes with singularity and not and exclude nodes with D=0
		
		NoSingularNum[IType][0] = 0;
		NoSingularNum[IType][1] = 0;
		NoSingularNum[IType][2] = 0;
		SingularNum[IType][0] = 0;
		SingularNum[IType][1] = 0;
		SingularNum[IType][2] = 0;
		if(QmobMask==NULL)
		{
			for(k=1;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						GrdPnt = jgrid+i;
						
						if(D[IType][GrdPnt]>0.0){
							BlackOrWhite=(k+j+i+World->startBlackAndWhite)%2;
							if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
							{
								SingularNum[IType][1]+=BlackOrWhite;
								SingularNum[IType][2]++;
							}
							else{
								NoSingularNum[IType][1]+=BlackOrWhite;
								NoSingularNum[IType][2]++;
							}
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
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						GrdPnt = jgrid+i;
						if(QmobMask[IType][GrdPnt]>0)
							if(D[IType][GrdPnt]>0.0)
						{
							BlackOrWhite=(k+j+i+World->startBlackAndWhite)%2;
							if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
							{
								SingularNum[IType][1]+=BlackOrWhite;
								SingularNum[IType][2]++;
							}
							else{
								NoSingularNum[IType][1]+=BlackOrWhite;
								NoSingularNum[IType][2]++;
							}
						}
					}
				}
			}
		}
			
		fprintf(stdout,"Nernst-Plank Solver:\n");
		fprintf(stdout,"		IonsQ:........................ %g\n", World->IonsQ[IType]);
		fprintf(stdout,"		SingularNum:...................... %d\n", SingularNum[IType][2]);
		fprintf(stdout,"		NoSingularNum:.................... %d\n", NoSingularNum[IType][2]);
		fprintf(stdout,"		IndexSingularParts:............... [%d,%d,%d]\n",
						SingularNum[IType][0],SingularNum[IType][1],SingularNum[IType][2]);
		fprintf(stdout,"		IndexNoSingularParts:............. [%d,%d,%d]\n",NoSingularNum[IType][0],NoSingularNum[IType][1],NoSingularNum[IType][2]);
		
		if(!(IndexSingular[IType] = new int[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		if(!(IndexNoSingular[IType] = new int[NoSingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		//cout<<"NP3sleep 10 sec\n";
		//sleep(10);
		for(i=0;i<6;i++)if(!(dix[IType][i] = new float[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		if(!(dixt[IType] = new float[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		int iSingular=0,iNoSingular=0;
		int iSingular2=SingularNum[IType][1],iNoSingular2=NoSingularNum[IType][1];
		if(QmobMask==NULL)
		{
			for(k=1;k<GS_Z-1;k++) 
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++) 
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++) 
					{
						GrdPnt = jgrid+i;
						if(D[IType][GrdPnt]>0.0)
						{
							BlackOrWhite=(k+j+i+World->startBlackAndWhite)%2;
							if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
							{
								if(BlackOrWhite)
								{
									IndexSingular[IType][iSingular]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iSingular++;
								}
								else
								{
									IndexSingular[IType][iSingular2]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular2, GrdPnt);
									iSingular2++;
								}
							}
							else
							{
								if(BlackOrWhite)
								{
									IndexNoSingular[IType][iNoSingular]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular++;
								}
								else
								{
									IndexNoSingular[IType][iNoSingular2]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular2++;
								}
							}
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
				for(j=1;j<GS_Y-1;j++) 
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++) 
					{
						GrdPnt = jgrid+i;
						if(D[IType][GrdPnt]>0.0)
						{
							BlackOrWhite=(k+j+i+World->startBlackAndWhite)%2;
							if(QmobMask[IType][GrdPnt]>0)
								if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
							{
								if(BlackOrWhite)
								{
									IndexSingular[IType][iSingular]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iSingular++;
								}
								else
								{
									IndexSingular[IType][iSingular2]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular2, GrdPnt);
									iSingular2++;
								}
							}
							else
							{
								if(BlackOrWhite)
								{
									IndexNoSingular[IType][iNoSingular]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular++;
								}
								else
								{
									IndexNoSingular[IType][iNoSingular2]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular2++;
								}
							}
						}
					}
				}
			}
		}
		
		for(i=0;i<SingularNum[IType][2];i++)
		{
			GrdPnt=IndexSingular[IType][i];
			dix[IType][0][i] = D[IType][GrdPnt+1]>0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt+1]):0;
			dix[IType][1][i] = D[IType][GrdPnt-1]>0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt-1]):0;
			dix[IType][2][i] = D[IType][GrdPnt+GS_X]>0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt+GS_X]):0;
			dix[IType][3][i] = D[IType][GrdPnt-GS_X]>0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt-GS_X]):0;
			dix[IType][4][i] = D[IType][GrdPnt+GS_XY]>0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt+GS_XY]):0;
			dix[IType][5][i] = D[IType][GrdPnt-GS_XY]>0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt-GS_XY]):0;
			dixt[IType][i] = dix[IType][0][i]+dix[IType][1][i]+dix[IType][2][i]+dix[IType][3][i]+dix[IType][4][i]+dix[IType][5][i];
		}
	}
	
	if(World->PMF!=NULL)fprintf(stdout,"Calculate NP with PMF\n");
	return 1;
}
int NernstPlankSolver::Solve()
{
	if(solver==Auto)
	{
		if(World->NIndexing==NULL)return NernstPlanckSolverW();
		else if(World->D!=NULL)return NernstPlanckSolverD();
		else return NernstPlanckSolverD();
	}
	else
	{
		if(solver==NodeIndexBased)
			return NernstPlanckSolverD();
		else if(solver==ArrayDirect)
			return NernstPlanckSolverD();
		else if(solver==PNPC)
			return NernstPlanckSolverW();
		else return NernstPlanckSolverD();
	}
}
int NernstPlankSolver::NernstPlanckSolverW()
{
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	float gridScale;
	int iteration;
	int i;
	long gridPoint;
	float om1,om2;
	float * potential;
	float * diffusion;
	long gridSizeXY;
	long gridSizeXYZ;
	float MaxChange;
	float change;
	int maxIterations;
	float convergence;
	float fpoh;
	float * positiveCharge;
	float * negativeCharge;
	float rodlt1,rodlt2;
	float phi0,phi1,phi2,phi3,phi4,phi5,phi6;
	float rophi1,rophi2;
	float term1,term2;
	float phit;
	float phit1,phit2;
	float denom1,denom2;
	long count;
	long * start;
	long * end;
	float * dixt;
	float * dix[6];
	long * borderPoints;
	long borderCount;
	long nextBorderPoint;
	long gridpx,gridnx;
	long gridpy,gridny;
	long gridpz,gridnz;
	float conv;
	long * delphiPoints;
	float * delphiCoeff[12];
	double delphiExp;
	long delphiCount;
	long nextDelphiPoint;
	long nextSpecialPoint;
	float maxDelphi;
	int expInterpolation;
	/*double ghalf;*/

	/*assert(output!=NULL);
	assert(output->diffusionMap!=NULL);
	assert(output->potentialMap!=NULL);
	assert(output->positiveDynamicChargeMap!=NULL);
	assert(output->negativeDynamicChargeMap!=NULL);*/

	assert(nernstPlanckSolverData!=NULL);

	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;

	maxIterations = nernstPlanckSolverData->maxIterations;
	convergence = nernstPlanckSolverData->convergence;
	expInterpolation =	nernstPlanckSolverData->expInterpolation;
	start = nernstPlanckSolverData->start;
	end = nernstPlanckSolverData->end;
	borderPoints = nernstPlanckSolverData->borderPoints;
	dixt = nernstPlanckSolverData->dixt;
	for(i=0;i<6;i++)
		dix[i] = nernstPlanckSolverData->dix[i];
	
	diffusion = World->D[0];
	potential = World->Potential;
	positiveCharge = World->C[0];
	negativeCharge = World->C[1];
	
	om1 = nernstPlanckSolverData->relaxation;
	om2 = 1-om1;

	fpoh = 4*M_PI*gridScale;
	conv = (gridScale*gridScale*gridScale)/COANGS;

	for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++) {
		potential[gridPoint]*=0.5;
		positiveCharge[gridPoint]/=fpoh;
		negativeCharge[gridPoint]/=fpoh;
	}

	count = 0;
	delphiCount = 0;
	maxDelphi=0;
	while(expInterpolation && end[count]<gridSizeXYZ) {
		for(gridPoint=start[count];gridPoint<end[count];gridPoint++) {
			gridpx = gridPoint+1;
			gridnx = gridPoint-1;
			gridpy = gridPoint+gridSizeX;
			gridny = gridPoint-gridSizeX;
			gridpz = gridPoint+gridSizeXY;
			gridnz = gridPoint-gridSizeXY;
			maxDelphi = maxDelphi>fabs(potential[gridPoint]-potential[gridpx])?maxDelphi:fabs(potential[gridPoint]-potential[gridpx]);
			maxDelphi = maxDelphi>fabs(potential[gridPoint]-potential[gridnx])?maxDelphi:fabs(potential[gridPoint]-potential[gridnx]);
			maxDelphi = maxDelphi>fabs(potential[gridPoint]-potential[gridpy])?maxDelphi:fabs(potential[gridPoint]-potential[gridpy]);
			maxDelphi = maxDelphi>fabs(potential[gridPoint]-potential[gridny])?maxDelphi:fabs(potential[gridPoint]-potential[gridny]);
			maxDelphi = maxDelphi>fabs(potential[gridPoint]-potential[gridpz])?maxDelphi:fabs(potential[gridPoint]-potential[gridpz]);
			maxDelphi = maxDelphi>fabs(potential[gridPoint]-potential[gridnz])?maxDelphi:fabs(potential[gridPoint]-potential[gridnz]);
			if(fabs(potential[gridPoint]-potential[gridpx])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) { 
				delphiCount++;
			}
			else if(fabs(potential[gridPoint]-potential[gridnx])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
				delphiCount++;
			}
			else if(fabs(potential[gridPoint]-potential[gridpy])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
				delphiCount++;
			}
			else if(fabs(potential[gridPoint]-potential[gridny])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
				delphiCount++;
			}
			else if(fabs(potential[gridPoint]-potential[gridpz])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
				delphiCount++;
			}
			else if(fabs(potential[gridPoint]-potential[gridnz])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
				delphiCount++;
			}
		}
		count++;
	}

	if(expInterpolation) {
		fprintf(stdout,"nernstPlanckSolver: there are %li delphi points\n",delphiCount);
		fprintf(stdout,"nernstPlanckSolver: maxDelphi = %f\n",maxDelphi);
	}

	delphiPoints = new long[delphiCount+1];
	assert(delphiPoints!=NULL);
	for(i=0;i<12 && expInterpolation;i++) {
		delphiCoeff[i] = new float[delphiCount];
		assert(delphiCoeff[i]!=NULL);
	}
	
	count = 0;
	delphiCount = 0;
	while(expInterpolation && end[count]<gridSizeXYZ) {
		for(gridPoint=start[count];gridPoint<end[count];gridPoint++) {
			gridpx = gridPoint+1;
			gridnx = gridPoint-1;
			gridpy = gridPoint+gridSizeX;
			gridny = gridPoint-gridSizeX;
			gridpz = gridPoint+gridSizeXY;
			gridnz = gridPoint-gridSizeXY;
			if(fabs(potential[gridPoint]-potential[gridpx])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE)
				delphiPoints[delphiCount++] = gridPoint;
			else if(fabs(potential[gridPoint]-potential[gridnx])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE)
				delphiPoints[delphiCount++] = gridPoint;
			else if(fabs(potential[gridPoint]-potential[gridpy])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE)
				delphiPoints[delphiCount++] = gridPoint;
			else if(fabs(potential[gridPoint]-potential[gridny])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE)
				delphiPoints[delphiCount++] = gridPoint;
			else if(fabs(potential[gridPoint]-potential[gridpz])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE)
				delphiPoints[delphiCount++] = gridPoint;
			else if(fabs(potential[gridPoint]-potential[gridnz])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE)
				delphiPoints[delphiCount++] = gridPoint;
			if(delphiCount>0 && delphiPoints[delphiCount-1]==gridPoint) {
				delphiCount--;
				for(i=0;i<12;i++)
					delphiCoeff[i][delphiCount] = 1;
				if(fabs(potential[gridPoint]-potential[gridpx])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
					delphiExp = exp(potential[gridPoint]-potential[gridpx]);
					/*ghalf = 2*(1-delphiExp)/(1-delphiExp*delphiExp);
					delphiCoeff[0][delphiCount] = ghalf;
					delphiCoeff[6][delphiCount] = 2-ghalf;*/
					delphiCoeff[0][delphiCount] = 1/(0.5*(1+delphiExp));
					delphiCoeff[6][delphiCount] = 1/(0.5*(1+1/delphiExp));
					/*delphiCoeff[12][delphiCount] = 2*(potential[gridpx]-potential[gridPoint])/(1/delphiExp-delphiExp);*/
				}
				if(fabs(potential[gridPoint]-potential[gridnx])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
					delphiExp = exp(potential[gridPoint]-potential[gridnx]);
					/*ghalf = 2*(1-delphiExp)/(1-delphiExp*delphiExp);
					delphiCoeff[1][delphiCount] = ghalf;
					delphiCoeff[7][delphiCount] = 2-ghalf;*/
					delphiCoeff[1][delphiCount] = 1/(0.5*(1+delphiExp));
					delphiCoeff[7][delphiCount] = 1/(0.5*(1+1/delphiExp));
					/*delphiCoeff[13][delphiCount] = 2*(potential[gridnx]-potential[gridPoint])/(1/delphiExp-delphiExp);*/
				}
				if(fabs(potential[gridPoint]-potential[gridpy])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
					delphiExp = exp(potential[gridPoint]-potential[gridpy]);
					/*ghalf = 2*(1-delphiExp)/(1-delphiExp*delphiExp);
					delphiCoeff[2][delphiCount] = ghalf;
					delphiCoeff[8][delphiCount] = 2-ghalf;*/
					delphiCoeff[2][delphiCount] = 1/(0.5*(1+delphiExp));
					delphiCoeff[8][delphiCount] = 1/(0.5*(1+1/delphiExp));
					/*delphiCoeff[14][delphiCount] = 2*(potential[gridpy]-potential[gridPoint])/(1/delphiExp-delphiExp);*/
				}
				if(fabs(potential[gridPoint]-potential[gridny])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
					delphiExp = exp(potential[gridPoint]-potential[gridny]);
					/*ghalf = 2*(1-delphiExp)/(1-delphiExp*delphiExp);
					delphiCoeff[3][delphiCount] = ghalf;
					delphiCoeff[9][delphiCount] = 2-ghalf;*/
					delphiCoeff[3][delphiCount] = 1/(0.5*(1+delphiExp));
					delphiCoeff[9][delphiCount] = 1/(0.5*(1+1/delphiExp));
					/*delphiCoeff[15][delphiCount] = 2*(potential[gridny]-potential[gridPoint])/(1/delphiExp-delphiExp);*/
				}
				if(fabs(potential[gridPoint]-potential[gridpz])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
					delphiExp = exp(potential[gridPoint]-potential[gridpz]);
					/*ghalf = 2*(1-delphiExp)/(1-delphiExp*delphiExp);
					delphiCoeff[4][delphiCount] = ghalf;
					delphiCoeff[10][delphiCount] = 2-ghalf;*/
					delphiCoeff[4][delphiCount] = 1/(0.5*(1+delphiExp));
					delphiCoeff[10][delphiCount] = 1/(0.5*(1+1/delphiExp));
					/*delphiCoeff[16][delphiCount] = 2*(potential[gridpz]-potential[gridPoint])/(1/delphiExp-delphiExp);*/
				}
				if(fabs(potential[gridPoint]-potential[gridnz])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
					delphiExp = exp(potential[gridPoint]-potential[gridnz]);
					/*ghalf = 2*(1-delphiExp)/(1-delphiExp*delphiExp);
					delphiCoeff[5][delphiCount] = ghalf;
					delphiCoeff[11][delphiCount] = 2-ghalf;*/
					delphiCoeff[5][delphiCount] = 1/(0.5*(1+delphiExp));
					delphiCoeff[11][delphiCount] = 1/(0.5*(1+1/delphiExp));
					/*delphiCoeff[17][delphiCount] = 2*(potential[gridnz]-potential[gridPoint])/(1/delphiExp-delphiExp);*/
				}
				delphiCount++;
			}
		}
		count++;
	}
	delphiPoints[delphiCount] = gridSizeXYZ;

	for(iteration=1;iteration<=maxIterations;iteration++) {
		count=0;
		borderCount = 0;
		delphiCount = 0;
		nextBorderPoint = borderPoints[borderCount];
		nextDelphiPoint = delphiPoints[delphiCount];
		nextSpecialPoint = nextBorderPoint<nextDelphiPoint?nextBorderPoint:nextDelphiPoint;

		if(iteration%NERNST_PLANCK_SOLVER_CONVERGENCE_CHECK==0) {
			MaxChange = 0;
			while(end[count]<gridSizeXYZ) {
				if(end[count]<nextSpecialPoint) {
					for(gridPoint=start[count];gridPoint<end[count];gridPoint++) {
						gridpx = gridPoint+1;
						gridnx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridny = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridnz = gridPoint-gridSizeXY;
						rodlt1 = positiveCharge[gridpx]+positiveCharge[gridnx]+positiveCharge[gridpy]+positiveCharge[gridny]+positiveCharge[gridpz]+positiveCharge[gridnz];
						rodlt2 = negativeCharge[gridpx]+negativeCharge[gridnx]+negativeCharge[gridpy]+negativeCharge[gridny]+negativeCharge[gridpz]+negativeCharge[gridnz];
						phit = potential[gridpx]+potential[gridnx]+potential[gridpy]+potential[gridny]+potential[gridpz]+potential[gridnz]-6*potential[gridPoint];
						change = om1*((rodlt1*(1-potential[gridPoint])+positiveCharge[gridpx]*potential[gridpx]+positiveCharge[gridnx]*potential[gridnx]+positiveCharge[gridpy]*potential[gridpy]+positiveCharge[gridny]*potential[gridny]+positiveCharge[gridpz]*potential[gridpz]+positiveCharge[gridnz]*potential[gridnz])/(6-phit)-positiveCharge[gridPoint]);
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						positiveCharge[gridPoint]+=change;
						change = om1*((rodlt2*(1+potential[gridPoint])-negativeCharge[gridpx]*potential[gridpx]-negativeCharge[gridnx]*potential[gridnx]-negativeCharge[gridpy]*potential[gridpy]-negativeCharge[gridny]*potential[gridny]-negativeCharge[gridpz]*potential[gridpz]-negativeCharge[gridnz]*potential[gridnz])/(6+phit)-negativeCharge[gridPoint]);
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						negativeCharge[gridPoint]+=change;
					}
				}
				else {
					for(gridPoint=start[count];gridPoint<end[count];gridPoint++) {
						gridpx = gridPoint+1;
						gridnx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridny = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridnz = gridPoint-gridSizeXY;
						if(gridPoint<nextSpecialPoint) {
							rodlt1 = positiveCharge[gridpx]+positiveCharge[gridnx]+positiveCharge[gridpy]+positiveCharge[gridny]+positiveCharge[gridpz]+positiveCharge[gridnz];
							rodlt2 = negativeCharge[gridpx]+negativeCharge[gridnx]+negativeCharge[gridpy]+negativeCharge[gridny]+negativeCharge[gridpz]+negativeCharge[gridnz];
							rophi1 = positiveCharge[gridpx]*potential[gridpx]+positiveCharge[gridnx]*potential[gridnx]+positiveCharge[gridpy]*potential[gridpy]+positiveCharge[gridny]*potential[gridny]+positiveCharge[gridpz]*potential[gridpz]+positiveCharge[gridnz]*potential[gridnz];
							rophi2 = -negativeCharge[gridpx]*potential[gridpx]-negativeCharge[gridnx]*potential[gridnx]-negativeCharge[gridpy]*potential[gridpy]-negativeCharge[gridny]*potential[gridny]-negativeCharge[gridpz]*potential[gridpz]-negativeCharge[gridnz]*potential[gridnz];
							term1 = rodlt1*(1-potential[gridPoint])+rophi1;
							term2 = rodlt2*(1+potential[gridPoint])+rophi2;
							phit = potential[gridpx]+potential[gridnx]+potential[gridpy]+potential[gridny]+potential[gridpz]+potential[gridnz]-6*potential[gridPoint];
							denom1 = 6-phit;
							denom2 = 6+phit;
						}
						else if(nextBorderPoint==nextDelphiPoint) {
							rodlt1 = positiveCharge[gridpx]*dix[0][borderCount]+positiveCharge[gridnx]*dix[1][borderCount]+positiveCharge[gridpy]*dix[2][borderCount]+positiveCharge[gridny]*dix[3][borderCount]+positiveCharge[gridpz]*dix[4][borderCount]+positiveCharge[gridnz]*dix[5][borderCount];
							rodlt2 = negativeCharge[gridpx]*dix[0][borderCount]+negativeCharge[gridnx]*dix[1][borderCount]+negativeCharge[gridpy]*dix[2][borderCount]+negativeCharge[gridny]*dix[3][borderCount]+negativeCharge[gridpz]*dix[4][borderCount]+negativeCharge[gridnz]*dix[5][borderCount];
							phi0 = potential[gridPoint];
							phi1 = dix[0][borderCount]*(potential[gridpx]-phi0);
							phi2 = dix[1][borderCount]*(potential[gridnx]-phi0);
							phi3 = dix[2][borderCount]*(potential[gridpy]-phi0);
							phi4 = dix[3][borderCount]*(potential[gridny]-phi0);
							phi5 = dix[4][borderCount]*(potential[gridpz]-phi0);
							phi6 = dix[5][borderCount]*(potential[gridnz]-phi0);
							rophi1 = positiveCharge[gridpx]*phi1*delphiCoeff[0][delphiCount]+positiveCharge[gridnx]*phi2*delphiCoeff[1][delphiCount]+positiveCharge[gridpy]*phi3*delphiCoeff[2][delphiCount]+positiveCharge[gridny]*phi4*delphiCoeff[3][delphiCount]+positiveCharge[gridpz]*phi5*delphiCoeff[4][delphiCount]+positiveCharge[gridnz]*phi6*delphiCoeff[5][delphiCount];
							rophi2 = -negativeCharge[gridpx]*phi1*delphiCoeff[6][delphiCount]-negativeCharge[gridnx]*phi2*delphiCoeff[7][delphiCount]-negativeCharge[gridpy]*phi3*delphiCoeff[8][delphiCount]-negativeCharge[gridny]*phi4*delphiCoeff[9][delphiCount]-negativeCharge[gridpz]*phi5*delphiCoeff[10][delphiCount]-negativeCharge[gridnz]*phi6*delphiCoeff[11][delphiCount];
							term1 = rodlt1+rophi1;
							term2 = rodlt2+rophi2;
							phit1 = phi1*delphiCoeff[6][delphiCount]+phi2*delphiCoeff[7][delphiCount]+phi3*delphiCoeff[8][delphiCount]+phi4*delphiCoeff[9][delphiCount]+phi5*delphiCoeff[10][delphiCount]+phi6*delphiCoeff[11][delphiCount];
							phit2 = phi1*delphiCoeff[0][delphiCount]+phi2*delphiCoeff[1][delphiCount]+phi3*delphiCoeff[2][delphiCount]+phi4*delphiCoeff[3][delphiCount]+phi5*delphiCoeff[4][delphiCount]+phi6*delphiCoeff[5][delphiCount];
							denom1 = dixt[borderCount]-phit1;
							denom2 = dixt[borderCount]+phit2;
							borderCount++;
							delphiCount++;
							nextBorderPoint = borderPoints[borderCount];
							nextDelphiPoint = delphiPoints[delphiCount];
							nextSpecialPoint = nextBorderPoint<nextDelphiPoint?nextBorderPoint:nextDelphiPoint;
						}
						else if(gridPoint==nextBorderPoint) {
							rodlt1 = positiveCharge[gridpx]*dix[0][borderCount]+positiveCharge[gridnx]*dix[1][borderCount]+positiveCharge[gridpy]*dix[2][borderCount]+positiveCharge[gridny]*dix[3][borderCount]+positiveCharge[gridpz]*dix[4][borderCount]+positiveCharge[gridnz]*dix[5][borderCount];
							rodlt2 = negativeCharge[gridpx]*dix[0][borderCount]+negativeCharge[gridnx]*dix[1][borderCount]+negativeCharge[gridpy]*dix[2][borderCount]+negativeCharge[gridny]*dix[3][borderCount]+negativeCharge[gridpz]*dix[4][borderCount]+negativeCharge[gridnz]*dix[5][borderCount];
							phi0 = potential[gridPoint];
							phi1 = dix[0][borderCount]*potential[gridpx];
							phi2 = dix[1][borderCount]*potential[gridnx];
							phi3 = dix[2][borderCount]*potential[gridpy];
							phi4 = dix[3][borderCount]*potential[gridny];
							phi5 = dix[4][borderCount]*potential[gridpz];
							phi6 = dix[5][borderCount]*potential[gridnz];
							rophi1 = positiveCharge[gridpx]*phi1+positiveCharge[gridnx]*phi2+positiveCharge[gridpy]*phi3+positiveCharge[gridny]*phi4+positiveCharge[gridpz]*phi5+positiveCharge[gridnz]*phi6;
							rophi2 = -negativeCharge[gridpx]*phi1-negativeCharge[gridnx]*phi2-negativeCharge[gridpy]*phi3-negativeCharge[gridny]*phi4-negativeCharge[gridpz]*phi5-negativeCharge[gridnz]*phi6;
							term1 = rodlt1*(1-phi0)+rophi1;
							term2 = rodlt2*(1+phi0)+rophi2;
							phit = phi1+phi2+phi3+phi4+phi5+phi6-dixt[borderCount]*phi0;
							denom1 = dixt[borderCount]-phit;
							denom2 = dixt[borderCount]+phit;
							borderCount++;
							nextBorderPoint = borderPoints[borderCount];
							nextSpecialPoint = nextBorderPoint<nextDelphiPoint?nextBorderPoint:nextDelphiPoint;
						}
						else {
							rodlt1 = positiveCharge[gridpx]+positiveCharge[gridnx]+positiveCharge[gridpy]+positiveCharge[gridny]+positiveCharge[gridpz]+positiveCharge[gridnz];
							rodlt2 = negativeCharge[gridpx]+negativeCharge[gridnx]+negativeCharge[gridpy]+negativeCharge[gridny]+negativeCharge[gridpz]+negativeCharge[gridnz];
							phi0 = potential[gridPoint];
							phi1 = potential[gridpx]-phi0;
							phi2 = potential[gridnx]-phi0;
							phi3 = potential[gridpy]-phi0;
							phi4 = potential[gridny]-phi0;
							phi5 = potential[gridpz]-phi0;
							phi6 = potential[gridnz]-phi0;
							rophi1 = positiveCharge[gridpx]*phi1*delphiCoeff[0][delphiCount]+positiveCharge[gridnx]*phi2*delphiCoeff[1][delphiCount]+positiveCharge[gridpy]*phi3*delphiCoeff[2][delphiCount]+positiveCharge[gridny]*phi4*delphiCoeff[3][delphiCount]+positiveCharge[gridpz]*phi5*delphiCoeff[4][delphiCount]+positiveCharge[gridnz]*phi6*delphiCoeff[5][delphiCount];
							rophi2 = -negativeCharge[gridpx]*phi1*delphiCoeff[6][delphiCount]-negativeCharge[gridnx]*phi2*delphiCoeff[7][delphiCount]-negativeCharge[gridpy]*phi3*delphiCoeff[8][delphiCount]-negativeCharge[gridny]*phi4*delphiCoeff[9][delphiCount]-negativeCharge[gridpz]*phi5*delphiCoeff[10][delphiCount]-negativeCharge[gridnz]*phi6*delphiCoeff[11][delphiCount];
							term1 = rodlt1+rophi1;
							term2 = rodlt2+rophi2;
							phit1 = phi1*delphiCoeff[6][delphiCount]+phi2*delphiCoeff[7][delphiCount]+phi3*delphiCoeff[8][delphiCount]+phi4*delphiCoeff[9][delphiCount]+phi5*delphiCoeff[10][delphiCount]+phi6*delphiCoeff[11][delphiCount];
							phit2 = phi1*delphiCoeff[0][delphiCount]+phi2*delphiCoeff[1][delphiCount]+phi3*delphiCoeff[2][delphiCount]+phi4*delphiCoeff[3][delphiCount]+phi5*delphiCoeff[4][delphiCount]+phi6*delphiCoeff[5][delphiCount];
							denom1 = 6-phit1;
							denom2 = 6+phit2;
							delphiCount++;
							nextDelphiPoint = delphiPoints[delphiCount];
							nextSpecialPoint = nextBorderPoint<nextDelphiPoint?nextBorderPoint:nextDelphiPoint;
						}
						change = om1*(term1/denom1-positiveCharge[gridPoint]);
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						positiveCharge[gridPoint]+=change;
						change = om1*(term2/denom2-negativeCharge[gridPoint]);
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						negativeCharge[gridPoint]+=change;
					}
				}
				count++;
			}
			MaxChange*=conv;
			fprintf(stdout,"nernstPlanckSolver: MaxChange = %e M at iteration %i\n",MaxChange,iteration);

			if(MaxChange<convergence)
				iteration = maxIterations+1;
		}
		else {
			while(end[count]<gridSizeXYZ) {
				if(end[count]<nextSpecialPoint) {
					for(gridPoint=start[count];gridPoint<end[count];gridPoint++) {
						gridpx = gridPoint+1;
						gridnx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridny = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridnz = gridPoint-gridSizeXY;
						rodlt1 = positiveCharge[gridpx]+positiveCharge[gridnx]+positiveCharge[gridpy]+positiveCharge[gridny]+positiveCharge[gridpz]+positiveCharge[gridnz];
						rodlt2 = negativeCharge[gridpx]+negativeCharge[gridnx]+negativeCharge[gridpy]+negativeCharge[gridny]+negativeCharge[gridpz]+negativeCharge[gridnz];
						phit = potential[gridpx]+potential[gridnx]+potential[gridpy]+potential[gridny]+potential[gridpz]+potential[gridnz]-6*potential[gridPoint];
						positiveCharge[gridPoint] = om2*positiveCharge[gridPoint]+om1*(rodlt1*(1-potential[gridPoint])+positiveCharge[gridpx]*potential[gridpx]+positiveCharge[gridnx]*potential[gridnx]+positiveCharge[gridpy]*potential[gridpy]+positiveCharge[gridny]*potential[gridny]+positiveCharge[gridpz]*potential[gridpz]+positiveCharge[gridnz]*potential[gridnz])/(6-phit);
						negativeCharge[gridPoint] = om2*negativeCharge[gridPoint]+om1*(rodlt2*(1+potential[gridPoint])-negativeCharge[gridpx]*potential[gridpx]-negativeCharge[gridnx]*potential[gridnx]-negativeCharge[gridpy]*potential[gridpy]-negativeCharge[gridny]*potential[gridny]-negativeCharge[gridpz]*potential[gridpz]-negativeCharge[gridnz]*potential[gridnz])/(6+phit);
					}
				}
				else {
					for(gridPoint=start[count];gridPoint<end[count];gridPoint++) {
						gridpx = gridPoint+1;
						gridnx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridny = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridnz = gridPoint-gridSizeXY;
						if(gridPoint<nextSpecialPoint) {
							rodlt1 = positiveCharge[gridpx]+positiveCharge[gridnx]+positiveCharge[gridpy]+positiveCharge[gridny]+positiveCharge[gridpz]+positiveCharge[gridnz];
							rodlt2 = negativeCharge[gridpx]+negativeCharge[gridnx]+negativeCharge[gridpy]+negativeCharge[gridny]+negativeCharge[gridpz]+negativeCharge[gridnz];
							rophi1 = positiveCharge[gridpx]*potential[gridpx]+positiveCharge[gridnx]*potential[gridnx]+positiveCharge[gridpy]*potential[gridpy]+positiveCharge[gridny]*potential[gridny]+positiveCharge[gridpz]*potential[gridpz]+positiveCharge[gridnz]*potential[gridnz];
							rophi2 = -negativeCharge[gridpx]*potential[gridpx]-negativeCharge[gridnx]*potential[gridnx]-negativeCharge[gridpy]*potential[gridpy]-negativeCharge[gridny]*potential[gridny]-negativeCharge[gridpz]*potential[gridpz]-negativeCharge[gridnz]*potential[gridnz];
							term1 = rodlt1*(1-potential[gridPoint])+rophi1;
							term2 = rodlt2*(1+potential[gridPoint])+rophi2;
							phit = potential[gridpx]+potential[gridnx]+potential[gridpy]+potential[gridny]+potential[gridpz]+potential[gridnz]-6*potential[gridPoint];
							denom1 = 6-phit;
							denom2 = 6+phit;
						}
						else if(nextBorderPoint==nextDelphiPoint) {
							rodlt1 = positiveCharge[gridpx]*dix[0][borderCount]+positiveCharge[gridnx]*dix[1][borderCount]+positiveCharge[gridpy]*dix[2][borderCount]+positiveCharge[gridny]*dix[3][borderCount]+positiveCharge[gridpz]*dix[4][borderCount]+positiveCharge[gridnz]*dix[5][borderCount];
							rodlt2 = negativeCharge[gridpx]*dix[0][borderCount]+negativeCharge[gridnx]*dix[1][borderCount]+negativeCharge[gridpy]*dix[2][borderCount]+negativeCharge[gridny]*dix[3][borderCount]+negativeCharge[gridpz]*dix[4][borderCount]+negativeCharge[gridnz]*dix[5][borderCount];
							phi0 = potential[gridPoint];
							phi1 = dix[0][borderCount]*(potential[gridpx]-phi0);
							phi2 = dix[1][borderCount]*(potential[gridnx]-phi0);
							phi3 = dix[2][borderCount]*(potential[gridpy]-phi0);
							phi4 = dix[3][borderCount]*(potential[gridny]-phi0);
							phi5 = dix[4][borderCount]*(potential[gridpz]-phi0);
							phi6 = dix[5][borderCount]*(potential[gridnz]-phi0);
							rophi1 = positiveCharge[gridpx]*phi1*delphiCoeff[0][delphiCount]+positiveCharge[gridnx]*phi2*delphiCoeff[1][delphiCount]+positiveCharge[gridpy]*phi3*delphiCoeff[2][delphiCount]+positiveCharge[gridny]*phi4*delphiCoeff[3][delphiCount]+positiveCharge[gridpz]*phi5*delphiCoeff[4][delphiCount]+positiveCharge[gridnz]*phi6*delphiCoeff[5][delphiCount];
							rophi2 = -negativeCharge[gridpx]*phi1*delphiCoeff[6][delphiCount]-negativeCharge[gridnx]*phi2*delphiCoeff[7][delphiCount]-negativeCharge[gridpy]*phi3*delphiCoeff[8][delphiCount]-negativeCharge[gridny]*phi4*delphiCoeff[9][delphiCount]-negativeCharge[gridpz]*phi5*delphiCoeff[10][delphiCount]-negativeCharge[gridnz]*phi6*delphiCoeff[11][delphiCount];
							term1 = rodlt1+rophi1;
							term2 = rodlt2+rophi2;
							phit1 = phi1*delphiCoeff[6][delphiCount]+phi2*delphiCoeff[7][delphiCount]+phi3*delphiCoeff[8][delphiCount]+phi4*delphiCoeff[9][delphiCount]+phi5*delphiCoeff[10][delphiCount]+phi6*delphiCoeff[11][delphiCount];
							phit2 = phi1*delphiCoeff[0][delphiCount]+phi2*delphiCoeff[1][delphiCount]+phi3*delphiCoeff[2][delphiCount]+phi4*delphiCoeff[3][delphiCount]+phi5*delphiCoeff[4][delphiCount]+phi6*delphiCoeff[5][delphiCount];
							denom1 = dixt[borderCount]-phit1;
							denom2 = dixt[borderCount]+phit2;
							borderCount++;
							delphiCount++;
							nextBorderPoint = borderPoints[borderCount];
							nextDelphiPoint = delphiPoints[delphiCount];
							nextSpecialPoint = nextBorderPoint<nextDelphiPoint?nextBorderPoint:nextDelphiPoint;
						}
						else if(gridPoint==nextBorderPoint) {
							rodlt1 = positiveCharge[gridpx]*dix[0][borderCount]+positiveCharge[gridnx]*dix[1][borderCount]+positiveCharge[gridpy]*dix[2][borderCount]+positiveCharge[gridny]*dix[3][borderCount]+positiveCharge[gridpz]*dix[4][borderCount]+positiveCharge[gridnz]*dix[5][borderCount];
							rodlt2 = negativeCharge[gridpx]*dix[0][borderCount]+negativeCharge[gridnx]*dix[1][borderCount]+negativeCharge[gridpy]*dix[2][borderCount]+negativeCharge[gridny]*dix[3][borderCount]+negativeCharge[gridpz]*dix[4][borderCount]+negativeCharge[gridnz]*dix[5][borderCount];
							phi0 = potential[gridPoint];
							phi1 = dix[0][borderCount]*potential[gridpx];
							phi2 = dix[1][borderCount]*potential[gridnx];
							phi3 = dix[2][borderCount]*potential[gridpy];
							phi4 = dix[3][borderCount]*potential[gridny];
							phi5 = dix[4][borderCount]*potential[gridpz];
							phi6 = dix[5][borderCount]*potential[gridnz];
							rophi1 = positiveCharge[gridpx]*phi1+positiveCharge[gridnx]*phi2+positiveCharge[gridpy]*phi3+positiveCharge[gridny]*phi4+positiveCharge[gridpz]*phi5+positiveCharge[gridnz]*phi6;
							rophi2 = -negativeCharge[gridpx]*phi1-negativeCharge[gridnx]*phi2-negativeCharge[gridpy]*phi3-negativeCharge[gridny]*phi4-negativeCharge[gridpz]*phi5-negativeCharge[gridnz]*phi6;
							term1 = rodlt1*(1-phi0)+rophi1;
							term2 = rodlt2*(1+phi0)+rophi2;
							phit = phi1+phi2+phi3+phi4+phi5+phi6-dixt[borderCount]*phi0;
							denom1 = dixt[borderCount]-phit;
							denom2 = dixt[borderCount]+phit;
							borderCount++;
							nextBorderPoint = borderPoints[borderCount];
							nextSpecialPoint = nextBorderPoint<nextDelphiPoint?nextBorderPoint:nextDelphiPoint;
						}
						else {
							rodlt1 = positiveCharge[gridpx]+positiveCharge[gridnx]+positiveCharge[gridpy]+positiveCharge[gridny]+positiveCharge[gridpz]+positiveCharge[gridnz];
							rodlt2 = negativeCharge[gridpx]+negativeCharge[gridnx]+negativeCharge[gridpy]+negativeCharge[gridny]+negativeCharge[gridpz]+negativeCharge[gridnz];
							phi0 = potential[gridPoint];
							phi1 = potential[gridpx]-phi0;
							phi2 = potential[gridnx]-phi0;
							phi3 = potential[gridpy]-phi0;
							phi4 = potential[gridny]-phi0;
							phi5 = potential[gridpz]-phi0;
							phi6 = potential[gridnz]-phi0;
							rophi1 = positiveCharge[gridpx]*phi1*delphiCoeff[0][delphiCount]+positiveCharge[gridnx]*phi2*delphiCoeff[1][delphiCount]+positiveCharge[gridpy]*phi3*delphiCoeff[2][delphiCount]+positiveCharge[gridny]*phi4*delphiCoeff[3][delphiCount]+positiveCharge[gridpz]*phi5*delphiCoeff[4][delphiCount]+positiveCharge[gridnz]*phi6*delphiCoeff[5][delphiCount];
							rophi2 = -negativeCharge[gridpx]*phi1*delphiCoeff[6][delphiCount]-negativeCharge[gridnx]*phi2*delphiCoeff[7][delphiCount]-negativeCharge[gridpy]*phi3*delphiCoeff[8][delphiCount]-negativeCharge[gridny]*phi4*delphiCoeff[9][delphiCount]-negativeCharge[gridpz]*phi5*delphiCoeff[10][delphiCount]-negativeCharge[gridnz]*phi6*delphiCoeff[11][delphiCount];
							term1 = rodlt1+rophi1;
							term2 = rodlt2+rophi2;
							phit1 = phi1*delphiCoeff[6][delphiCount]+phi2*delphiCoeff[7][delphiCount]+phi3*delphiCoeff[8][delphiCount]+phi4*delphiCoeff[9][delphiCount]+phi5*delphiCoeff[10][delphiCount]+phi6*delphiCoeff[11][delphiCount];
							phit2 = phi1*delphiCoeff[0][delphiCount]+phi2*delphiCoeff[1][delphiCount]+phi3*delphiCoeff[2][delphiCount]+phi4*delphiCoeff[3][delphiCount]+phi5*delphiCoeff[4][delphiCount]+phi6*delphiCoeff[5][delphiCount];
							denom1 = 6-phit1;
							denom2 = 6+phit2;
							delphiCount++;
							nextDelphiPoint = delphiPoints[delphiCount];
							nextSpecialPoint = nextBorderPoint<nextDelphiPoint?nextBorderPoint:nextDelphiPoint;
						}
						positiveCharge[gridPoint] = om2*positiveCharge[gridPoint]+om1*term1/denom1;
						negativeCharge[gridPoint] = om2*negativeCharge[gridPoint]+om1*term2/denom2;
					}
				}
				count++;
			}
		}
	}

	free(delphiPoints);
	for(i=0;i<12 && expInterpolation;i++)
		free(delphiCoeff[i]);

	for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++) {
		potential[gridPoint]*=2;
		positiveCharge[gridPoint]*=fpoh;
		negativeCharge[gridPoint]*=fpoh;
	}

	return EXIT_SUCCESS;
}
int NernstPlankSolver::CheckSystem()
{
	int i,j,k;
	int IonType;
	int gridPoint;
	
	VectorIntField3D *iV=new VectorIntField3D(World->GridSize,World->GridScale,World->NIonsTypes);
	iV->FillValue(0);
	for(IonType=0;IonType<World->NIonsTypes;IonType++)
	{
		for(i=0;i<NoSingularNum[IonType][2];i++)
		{
			gridPoint=IndexNoSingular[IonType][i];
			iV->V[IonType][gridPoint]=1;
		}
		for(i=0;i<SingularNum[IonType][2];i++)
		{
			gridPoint=IndexSingular[IonType][i];
			iV->V[IonType][gridPoint]=2;
		}
	}
	
	int kgrid,jgrid,GrdPnt;
	int GS_X=World->GS_X;
	int GS_Y=World->GS_Y;
	int GS_Z=World->GS_Z;
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_X*GS_Y*GS_Z;
	
	for(IonType=0;IonType<World->NIonsTypes;IonType++)
	{
		for(k=1;k<GS_Z-1;k++)
		{
			for(j=1;j<GS_Y-1;j++)
			{
				for(i=1;i<GS_X-1;i++)
				{
					GrdPnt = i+j*GS_X+k*GS_XY;
					
					if(World->C[IonType][GrdPnt]==0.0&&iV->V[IonType][GrdPnt]!=0)
					{
						pnpPrint0("World->C[%d][%d]==0.0&&iV->V[%d][%d]!=0\n",IonType,GrdPnt,IonType,GrdPnt);
					}
					if(World->NIndexing->GetDiffFloat(IonType,GrdPnt)==0.0&&iV->V[IonType][GrdPnt]!=0)
					{
						pnpPrint0("NIndexing->GetDiffFloat(%d,%d)==0.0&&iV->V[%d][%d]!=0\n",IonType,GrdPnt,IonType,GrdPnt);
					}
				}
			}
		}
	}
	if(World->PMF!=NULL)
	{
		
		
		pnpPrint0("PMF present, will check values:\n");
		for(IonType=0;IonType<World->NIonsTypes;IonType++)
		{
			float dPMFx,dPMFy,dPMFz;
			int countPMFx=0,countPMFy=0,countPMFz=0;
			for(k=1;k<GS_Z;k++)
			{
				for(j=1;j<GS_Y;j++)
				{
					for(i=1;i<GS_X;i++)
					{
						GrdPnt = i+j*GS_X+k*GS_XY;
						dPMFx=World->PMF[IonType][GrdPnt]-World->PMF[IonType][GrdPnt-1];
						dPMFy=World->PMF[IonType][GrdPnt]-World->PMF[IonType][GrdPnt-GS_X];
						dPMFz=World->PMF[IonType][GrdPnt]-World->PMF[IonType][GrdPnt-GS_XY];
						
						if((dPMFx*dPMFx)>4.0)countPMFx++;
						if((dPMFy*dPMFy)>4.0)countPMFy++;
						if((dPMFz*dPMFz)>4.0)countPMFz++;
					}
				}
			}
			pnpPrint0("\tIon %d: countPMFx=%d,countPMFy=%d,countPMFz=%d Throught All Space\n"
					,IonType,countPMFx,countPMFy,countPMFz);
		}
		for(IonType=0;IonType<World->NIonsTypes;IonType++)
		{
			float dPMFx,dPMFy,dPMFz;
			int countPMFx=0,countPMFy=0,countPMFz=0;
			for(k=1;k<GS_Z;k++)
			{
				for(j=1;j<GS_Y;j++)
				{
					for(i=1;i<GS_X;i++)
					{
						GrdPnt = i+j*GS_X+k*GS_XY;
						if(World->NIndexing->GetDiffFloat(IonType,GrdPnt)!=0.0)
						{
							dPMFx=World->PMF[IonType][GrdPnt]-World->PMF[IonType][GrdPnt-1];
							dPMFy=World->PMF[IonType][GrdPnt]-World->PMF[IonType][GrdPnt-GS_X];
							dPMFz=World->PMF[IonType][GrdPnt]-World->PMF[IonType][GrdPnt-GS_XY];
							
							if((dPMFx*dPMFx)>4.0&&World->NIndexing->GetDiffFloat(IonType,GrdPnt-1)!=0.0)
								countPMFx++;
							if((dPMFy*dPMFy)>4.0&&World->NIndexing->GetDiffFloat(IonType,GrdPnt-GS_X)!=0.0)
								countPMFy++;
							if((dPMFz*dPMFz)>4.0&&World->NIndexing->GetDiffFloat(IonType,GrdPnt-GS_XY)!=0.0)
								countPMFz++;
						}
					}
				}
			}
			pnpPrint0("\tIon %d: countPMFx=%d,countPMFy=%d,countPMFz=%d Throught D!=0\n"
					,IonType,countPMFx,countPMFy,countPMFz);
		}
		pnpPrint0("PMF checked\n");
	}
	delete iV;
	return EXIT_SUCCESS;
}
int NernstPlankSolver::NernstPlanckSolverD()
{
	int status;
	PNP_EXIT_FAIL_NULL(World,"World is not initialize\n");
	int IonType;
	float MaxChangeTMP=0.0;
	for(IonType=0;IonType<World->NIonsTypes;IonType++)
	{
		NernstPlanckSolverDpre(IonType);
		status=NernstPlanckSolverDiter(IonType,true);
		NernstPlanckSolverDpost(IonType);
		if(MaxChange>MaxChangeTMP)MaxChangeTMP=MaxChange;
	}
	MaxChange=MaxChangeTMP;
	return status;
}
int NernstPlankSolver::NernstPlanckSolverDpre(int IonType)
{
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	float gridScale;
	int iteration;
	int i,j,k,i1;
	int gridPoint;
	float om1,om2;
	float * potential;
	//float * diffusion;
	int gridSizeXY;
	int gridSizeXYZ;
	float change;
	float fpoh;
	//float * positiveCharge;
	//float * negativeCharge;
	float rodlt1,rodlt2;
	float phi0,phi1,phi2,phi3,phi4,phi5,phi6;
	float rophi1,rophi2;
	float term1,term2;
	float phit;
	float phit1,phit2;
	//float denom1,denom2;
	int gridpx,gridmx;
	int gridpy,gridmy;
	int gridpz,gridmz;
	float conv;
	bool *PeriodicBoundaryCondition;
	int itmp1,itmp2;
	float *C[2];
	//int IonType;
	
	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;	
	//diffusion = World->D[0];
	potential = World->Potential;
	//positiveCharge = World->C[0];
	//negativeCharge = World->C[1];
	C[0] = World->C[0];
	C[1] = World->C[1];
	PeriodicBoundaryCondition=World->PBC;
	
	om1 = Relaxation;
	om2 = 1-om1;

	fpoh = 4*M_PI*World->GridScale;
	conv = (World->GridScale*World->GridScale*World->GridScale)/COANGS;

	if(UTMPSingle==NULL)
		UTMPSingle = new float[gridSizeXYZ];
	if(TMPSingle == NULL)
	{
		TMPSingle = new float[gridSizeXYZ];
		for(j=0;j<gridSizeXYZ;j++)TMPSingle[j]=0.0f;
	}
	float *UTMP=UTMPSingle;
	float *TMP=TMPSingle;
	
	int MaxSingularNum=0,IonType2;
	for(IonType2=0;IonType2<World->NIonsTypes;IonType2++)
	{
		if(MaxSingularNum<SingularNum[IonType2][2])
			MaxSingularNum=SingularNum[IonType2][2];
	}
	for(i=0;i<6;i++)
		if(dphi[i] == NULL)
			if(!(dphi[i] = new float[MaxSingularNum]))
	{
		cerr << "ERROR 204: No memory available\n";
		exit(204);
	}
	//Set up the Potential energy
	if(World->PMF!=NULL)
	{
		if(potential!=NULL)
		{
			if(PMFWeight>=0.0)
			{
				//for(i=0;i<World->NIonsTypes;i++)
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[j]=0.5*(World->IonsQ[IonType]*potential[j]+World->PMF[IonType][j]*PMFWeight);
			}
			else
			{
				//for(i=0;i<World->NIonsTypes;i++)
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[j]=0.5*(World->IonsQ[IonType]*potential[j]+World->PMF[IonType][j]);
			}
		}
		else
		{
			if(PMFWeight>=0.0)
			{
				//for(i=0;i<World->NIonsTypes;i++)
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[j]=0.5*(World->PMF[IonType][j])*PMFWeight;
			}
			else
			{
				//for(i=0;i<World->NIonsTypes;i++)
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[j]=0.5*(World->PMF[IonType][j]);
			}
		}
	}
	else
		for(i=0;i<World->NIonsTypes;i++)
			for(j=0;j<gridSizeXYZ;j++)UTMP[j]=0.5*(World->IonsQ[IonType]*potential[j]);

	//for(IonType=0;IonType<World->NIonsTypes;IonType++)
	//{
		for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
			C[IonType][gridPoint]/=fpoh;

		for(i1=0;i1<NoSingularNum[IonType][2];i1++)
		{
			gridPoint=IndexNoSingular[IonType][i1];
			gridpx = gridPoint+1;
			gridmx = gridPoint-1;
			gridpy = gridPoint+gridSizeX;
			gridmy = gridPoint-gridSizeX;
			gridpz = gridPoint+gridSizeXY;
			gridmz = gridPoint-gridSizeXY;
	
			phit = UTMP[gridpx] + UTMP[gridmx] + UTMP[gridpy] + UTMP[gridmy] + UTMP[gridpz] + UTMP[gridmz] - 6.0*UTMP[gridPoint];
			TMP[gridPoint]=1.0/(6.0-phit);
		}
		for(i=0;i<SingularNum[IonType][2];i++)
		{
			gridPoint=IndexSingular[IonType][i];
			gridpx = gridPoint+1;
			gridmx = gridPoint-1;
			gridpy = gridPoint+gridSizeX;
			gridmy = gridPoint-gridSizeX;
			gridpz = gridPoint+gridSizeXY;
			gridmz = gridPoint-gridSizeXY;
		
			dphi[0][i] = dix[IonType][0][i]*UTMP[gridpx];
			dphi[1][i] = dix[IonType][1][i]*UTMP[gridmx];
			dphi[2][i] = dix[IonType][2][i]*UTMP[gridpy];
			dphi[3][i] = dix[IonType][3][i]*UTMP[gridmy];
			dphi[4][i] = dix[IonType][4][i]*UTMP[gridpz];
			dphi[5][i] = dix[IonType][5][i]*UTMP[gridmz];
			phit = dphi[0][i] + dphi[1][i] + dphi[2][i] + dphi[3][i] + dphi[4][i] + dphi[5][i] - dixt[IonType][i]*UTMP[gridPoint];
			TMP[gridPoint] = 1.0/(dixt[IonType][i]-phit);
		}
	//}
	return EXIT_SUCCESS;
}
int NernstPlankSolver::NernstPlanckSolverDiter(int IonType,bool calcchange)
{
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	float gridScale;
	int iteration;
	int i,j,k,i1;
	int gridPoint;
	float om1,om2;
	float * potential;
	//float * diffusion;
	int gridSizeXY;
	int gridSizeXYZ;
	float change;
	float fpoh;
	//float * positiveCharge;
	//float * negativeCharge;
	float rodlt1,rodlt2;
	float phi0,phi1,phi2,phi3,phi4,phi5,phi6;
	float rophi1,rophi2;
	float term1,term2;
	float phit;
	float phit1,phit2;
	//float denom1,denom2;
	int gridpx,gridmx;
	int gridpy,gridmy;
	int gridpz,gridmz;
	float conv;
	bool *PeriodicBoundaryCondition;
	int itmp1,itmp2;
	float *C[2];
	//int IonType;
	
	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;	
	//diffusion = World->D[0];
	potential = World->Potential;
	//positiveCharge = World->C[0];
	//negativeCharge = World->C[1];
	C[0] = World->C[0];
	C[1] = World->C[1];
	float *UTMP=UTMPSingle;
	float *TMP=TMPSingle;
	PeriodicBoundaryCondition=World->PBC;
	
	om1 = Relaxation;
	om2 = 1-om1;

	fpoh = 4*M_PI*World->GridScale;
	conv = (World->GridScale*World->GridScale*World->GridScale)/COANGS;
	
	for(iteration=1;iteration<=MaxIterations;iteration++)
	{
		if(calcchange&&((verbose&&(iteration%ConvergenceCheck==0))||(iteration==MaxIterations)))
		{
			MaxChange = 0;
			//for(IonType=0;IonType<World->NIonsTypes;IonType++)
			//{
				int countNegative=0;
				for(j=0;j<2;j++)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
					
						rodlt1 = C[IonType][gridpx] + C[IonType][gridmx] + C[IonType][gridpy] + C[IonType][gridmy] + C[IonType][gridpz] + C[IonType][gridmz];
						change = ( rodlt1 * ( 1.0 - UTMP[gridPoint] ) + C[IonType][gridpx]*UTMP[gridpx] + C[IonType][gridmx]*UTMP[gridmx] + C[IonType][gridpy]*UTMP[gridpy] + C[IonType][gridmy]*UTMP[gridmy] + C[IonType][gridpz]*UTMP[gridpz] + C[IonType][gridmz]*UTMP[gridmz]) * TMP[gridPoint];
						change -= C[IonType][gridPoint];
						change *= om1;
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						C[IonType][gridPoint]+=change;
						
//						 if(C[IonType][gridPoint]<0.0)
//						 {
//							 fprintf(stderr,"C[%d][%d]=%g<0.0\n",IonType,gridPoint,[IonType][gridPoint]);
//							 return EXIT_FAILURE;
//						 }
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
						
						rodlt1 = C[IonType][gridpx]*dix[IonType][0][i] + C[IonType][gridmx]*dix[IonType][1][i] + C[IonType][gridpy]*dix[IonType][2][i] + C[IonType][gridmy]*dix[IonType][3][i] + C[IonType][gridpz]*dix[IonType][4][i] + C[IonType][gridmz]*dix[IonType][5][i];
						rophi1 = C[IonType][gridpx]*dphi[0][i] + C[IonType][gridmx]*dphi[1][i] + C[IonType][gridpy]*dphi[2][i] + C[IonType][gridmy]*dphi[3][i] + C[IonType][gridpz]*dphi[4][i] + C[IonType][gridmz]*dphi[5][i];
						term1 = rodlt1*(1.0-UTMP[gridPoint])+rophi1;
						change = om1*(term1*TMP[gridPoint]-C[IonType][gridPoint]);
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						C[IonType][gridPoint]+=change;
					}
					if(bConcCorrection)
					{
						for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
						{
							gridPoint=IndexNoSingular[IonType][i];
							if(C[IonType][gridPoint]<0.0)
							{
								C[IonType][gridPoint]=0.0;
								if(countNegative<=10)
									PntsWithNegativeConc[countNegative]=gridPoint;
								countNegative++;
							}
						}
						for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++)
						{
							gridPoint=IndexSingular[IonType][i];
							if(C[IonType][gridPoint]<0.0)
							{
								C[IonType][gridPoint]=0.0;
								if(countNegative<10)
									PntsWithNegativeConc[countNegative]=gridPoint;
								countNegative++;
							}
						}
					}
					World->BorderExchange(C[IonType]);
				}
				if(bConcCorrection&&countNegative>0)
				{
					pnpWarning("IonType %d have %d points with value<0.0\n",IonType,countNegative);
					if(countNegative<10)
					{
						pnpWarning("\tpnts = ");
						for(i=0;i<countNegative;i++)
							pnpPrint(" %d",PntsWithNegativeConc[i]);
						pnpPrint("\n");
					}
				}
			//}
			MaxChange*=conv;
#ifdef MPI_PARALLEL
			int dest;
			pnpsapp->MyComGroup.Barrier();
			
			if(World->MyRank==0)
			{
				for(dest=1;dest<World->NProcs;dest++)
				{
					pnpsapp->MyComGroup.Recv(&change, 1, MPI::FLOAT, dest, 0);
					MaxChange = MaxChange>change?MaxChange:change;
					MaxChange = -MaxChange<change?MaxChange:-change;
				}
			}
			else
			{
				pnpsapp->MyComGroup.Send(&MaxChange, 1, MPI::FLOAT, 0, 0);
			}
#endif
			if(verbose)if(World->MyRank==0)
			{
				fprintf(stdout, "<NernstPlankIteration Nit=\"%6d\"	dC=\"%16.8g\"/>\n", iteration, MaxChange);
			}
#if defined(_MSC_VER) || defined(__DECCXX) || (__GNUC__ < 3)
			if(MaxChange>1E13)
#else
				if(MaxChange>1E13||MaxChange==NAN)
#endif
			{
				pnpError("MaxChange=%f will finish\n",MaxChange);
				return EXIT_FAILURE;
			}
			if(MaxChange<Convergence)
				iteration = MaxIterations+1;
		}
		else
		{
			//for(IonType=0;IonType<World->NIonsTypes;IonType++)
			//{
				for(j=0;j<2;j++)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
						rodlt1 = C[IonType][gridpx] + C[IonType][gridmx] + C[IonType][gridpy] + C[IonType][gridmy] + C[IonType][gridpz] + C[IonType][gridmz];
						change = ( rodlt1 * ( 1.0 - UTMP[gridPoint] ) + C[IonType][gridpx]*UTMP[gridpx] + C[IonType][gridmx]*UTMP[gridmx] + C[IonType][gridpy]*UTMP[gridpy] + C[IonType][gridmy]*UTMP[gridmy] + C[IonType][gridpz]*UTMP[gridpz] + C[IonType][gridmz]*UTMP[gridmz]) * TMP[gridPoint];
						change -= C[IonType][gridPoint];
						change *= om1;
						C[IonType][gridPoint]+=change;
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++){
						gridPoint=IndexSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
	
						rodlt1 = C[IonType][gridpx]*dix[IonType][0][i] + C[IonType][gridmx]*dix[IonType][1][i] + C[IonType][gridpy]*dix[IonType][2][i] + C[IonType][gridmy]*dix[IonType][3][i] + C[IonType][gridpz]*dix[IonType][4][i] + C[IonType][gridmz]*dix[IonType][5][i];
						rophi1 = C[IonType][gridpx]*dphi[0][i] + C[IonType][gridmx]*dphi[1][i] + C[IonType][gridpy]*dphi[2][i] + C[IonType][gridmy]*dphi[3][i] + C[IonType][gridpz]*dphi[4][i] + C[IonType][gridmz]*dphi[5][i];
						term1 = rodlt1*(1.0-UTMP[gridPoint])+rophi1;
						change = om1*(term1*TMP[gridPoint]-C[IonType][gridPoint]);
						C[IonType][gridPoint]+=change;
					}
					if(bConcCorrection)
					{
						for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
						{
							gridPoint=IndexNoSingular[IonType][i];
							if(C[IonType][gridPoint]<0.0)
							{
								C[IonType][gridPoint]=0.0;
							}
						}
						for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++){
							gridPoint=IndexSingular[IonType][i];
							if(C[IonType][gridPoint]<0.0)
							{
								C[IonType][gridPoint]=0.0;
							}
						}
					}
					World->BorderExchange(C[IonType]);
				}
			//}
		}
	}
	return EXIT_SUCCESS;
}
int NernstPlankSolver::NernstPlanckSolverDiterHalfStep(int IonType,bool calcchange,int BlackOrWhite)
{
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	float gridScale;
	int iteration;
	int i,j,k,i1;
	int gridPoint;
	float om1,om2;
	float * potential;
	//float * diffusion;
	int gridSizeXY;
	int gridSizeXYZ;
	float change;
	float fpoh;
	//float * positiveCharge;
	//float * negativeCharge;
	float rodlt1,rodlt2;
	float phi0,phi1,phi2,phi3,phi4,phi5,phi6;
	float rophi1,rophi2;
	float term1,term2;
	float phit;
	float phit1,phit2;
	//float denom1,denom2;
	int gridpx,gridmx;
	int gridpy,gridmy;
	int gridpz,gridmz;
	float conv;
	bool *PeriodicBoundaryCondition;
	int itmp1,itmp2;
	float *C[2];
	//int IonType;
	
	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;	
	//diffusion = World->D[0];
	potential = World->Potential;
	//positiveCharge = World->C[0];
	//negativeCharge = World->C[1];
	C[0] = World->C[0];
	C[1] = World->C[1];
	float *UTMP=UTMPSingle;
	float *TMP=TMPSingle;
	PeriodicBoundaryCondition=World->PBC;
	
	om1 = Relaxation;
	om2 = 1-om1;

	fpoh = 4*M_PI*World->GridScale;
	conv = (World->GridScale*World->GridScale*World->GridScale)/COANGS;
	
	for(iteration=1;iteration<=MaxIterations;iteration++)
	{
			MaxChange = 0;
			//for(IonType=0;IonType<World->NIonsTypes;IonType++)
			//{
			int countNegative=0;
			j=BlackOrWhite;
				for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
				{
					gridPoint=IndexNoSingular[IonType][i];
					gridpx = gridPoint+1;
					gridmx = gridPoint-1;
					gridpy = gridPoint+gridSizeX;
					gridmy = gridPoint-gridSizeX;
					gridpz = gridPoint+gridSizeXY;
					gridmz = gridPoint-gridSizeXY;
					
					rodlt1 = C[IonType][gridpx] + C[IonType][gridmx] + C[IonType][gridpy] + C[IonType][gridmy] + C[IonType][gridpz] + C[IonType][gridmz];
					change = ( rodlt1 * ( 1.0 - UTMP[gridPoint] ) + C[IonType][gridpx]*UTMP[gridpx] + C[IonType][gridmx]*UTMP[gridmx] + C[IonType][gridpy]*UTMP[gridpy] + C[IonType][gridmy]*UTMP[gridmy] + C[IonType][gridpz]*UTMP[gridpz] + C[IonType][gridmz]*UTMP[gridmz]) * TMP[gridPoint];
					change -= C[IonType][gridPoint];
					change *= om1;
					MaxChange = MaxChange>change?MaxChange:change;
					MaxChange = -MaxChange<change?MaxChange:-change;
					C[IonType][gridPoint]+=change;
						
//						 if(C[IonType][gridPoint]<0.0)
//						 {
//							 fprintf(stderr,"C[%d][%d]=%g<0.0\n",IonType,gridPoint,[IonType][gridPoint]);
//							 return EXIT_FAILURE;
//						 }
				}
				for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++)
				{
					gridPoint=IndexSingular[IonType][i];
					gridpx = gridPoint+1;
					gridmx = gridPoint-1;
					gridpy = gridPoint+gridSizeX;
					gridmy = gridPoint-gridSizeX;
					gridpz = gridPoint+gridSizeXY;
					gridmz = gridPoint-gridSizeXY;
						
					rodlt1 = C[IonType][gridpx]*dix[IonType][0][i] + C[IonType][gridmx]*dix[IonType][1][i] + C[IonType][gridpy]*dix[IonType][2][i] + C[IonType][gridmy]*dix[IonType][3][i] + C[IonType][gridpz]*dix[IonType][4][i] + C[IonType][gridmz]*dix[IonType][5][i];
					rophi1 = C[IonType][gridpx]*dphi[0][i] + C[IonType][gridmx]*dphi[1][i] + C[IonType][gridpy]*dphi[2][i] + C[IonType][gridmy]*dphi[3][i] + C[IonType][gridpz]*dphi[4][i] + C[IonType][gridmz]*dphi[5][i];
					term1 = rodlt1*(1.0-UTMP[gridPoint])+rophi1;
					change = om1*(term1*TMP[gridPoint]-C[IonType][gridPoint]);
					MaxChange = MaxChange>change?MaxChange:change;
					MaxChange = -MaxChange<change?MaxChange:-change;
					C[IonType][gridPoint]+=change;
				}
				if(bConcCorrection)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						if(C[IonType][gridPoint]<0.0)
						{
							C[IonType][gridPoint]=0.0;
							if(countNegative<=10)
								PntsWithNegativeConc[countNegative]=gridPoint;
							countNegative++;
						}
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexSingular[IonType][i];
						if(C[IonType][gridPoint]<0.0)
						{
							C[IonType][gridPoint]=0.0;
							if(countNegative<10)
								PntsWithNegativeConc[countNegative]=gridPoint;
							countNegative++;
						}
					}
				}
				World->BorderExchange(C[IonType]);
			
			if(bConcCorrection&&countNegative>0)
			{
				pnpWarning("IonType %d have %d points with value<0.0\n",IonType,countNegative);
				if(countNegative<10)
				{
					pnpWarning("\tpnts = ");
					for(i=0;i<countNegative;i++)
						pnpPrint(" %d",PntsWithNegativeConc[i]);
					pnpPrint("\n");
				}
			}
			//}
			MaxChange*=conv;
#ifdef MPI_PARALLEL
			int dest;
			pnpsapp->MyComGroup.Barrier();
			
			if(World->MyRank==0)
			{
				for(dest=1;dest<World->NProcs;dest++)
				{
					pnpsapp->MyComGroup.Recv(&change, 1, MPI::FLOAT, dest, 0);
					MaxChange = MaxChange>change?MaxChange:change;
					MaxChange = -MaxChange<change?MaxChange:-change;
				}
			}
			else
			{
				pnpsapp->MyComGroup.Send(&MaxChange, 1, MPI::FLOAT, 0, 0);
			}
#endif
			if(verbose)if(World->MyRank==0)
			{
				fprintf(stdout, "<NernstPlankIteration Nit=\"%6d\"	dC=\"%16.8g\"/>\n", iteration, MaxChange);
			}
#if defined(_MSC_VER) || defined(__DECCXX) || (__GNUC__ < 3)
			if(MaxChange>1E13)
#else
				if(MaxChange>1E13||MaxChange==NAN)
#endif
			{
				pnpError("MaxChange=%f will finish\n",MaxChange);
				return EXIT_FAILURE;
			}
			if(MaxChange<Convergence)
				iteration = MaxIterations+1;
	}
	return EXIT_SUCCESS;
}
int NernstPlankSolver::NernstPlanckSolverDpost(int IonType)
{
	int i1,gridPoint;
	int gridSizeXYZ=World->GS_XYZ;
	float fpoh = 4*M_PI*World->GridScale;
	
	//for(i1=0;i1<World->NIonsTypes;i1++)
		for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
			World->C[IonType][gridPoint]*=fpoh;
	return EXIT_SUCCESS;
}

int NernstPlankSolver::NernstPlanckSolverDDouble()
{
	int status;
	PNP_EXIT_FAIL_NULL(World,"World is not initialize\n");
	int IonType;
	float MaxChangeTMP=0.0;
	for(IonType=0;IonType<World->NIonsTypes;IonType++)
	{
		NernstPlanckSolverDpreDouble(IonType);
		status=NernstPlanckSolverDiterDouble(IonType,true);
		NernstPlanckSolverDpostDouble(IonType);
		if(MaxChange>MaxChangeTMP)MaxChangeTMP=MaxChange;
	}
	MaxChange=MaxChangeTMP;
	return status;
}
int NernstPlankSolver::NernstPlanckSolverDpreDouble(int IonType)
{
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	double gridScale;
	int iteration;
	int i,j,k,i1;
	int gridPoint;
	double om1,om2;
	double * potential;
	//float * diffusion;
	int gridSizeXY;
	int gridSizeXYZ;
	double change;
	double fpoh;
	//float * positiveCharge;
	//float * negativeCharge;
	double rodlt1,rodlt2;
	double phi0,phi1,phi2,phi3,phi4,phi5,phi6;
	double rophi1,rophi2;
	double term1,term2;
	double phit;
	double phit1,phit2;
	//float denom1,denom2;
	int gridpx,gridmx;
	int gridpy,gridmy;
	int gridpz,gridmz;
	double conv;
	bool *PeriodicBoundaryCondition;
	int itmp1,itmp2;
	double *C[2];
	//int IonType;
	
	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;	
	//diffusion = World->D[0];
	potential = World->PotentialDouble;
	//positiveCharge = World->C[0];
	//negativeCharge = World->C[1];
	C[0] = World->CDouble[0];
	C[1] = World->CDouble[1];
	PeriodicBoundaryCondition=World->PBC;
	
	om1 = Relaxation;
	om2 = 1-om1;

	fpoh = 4*M_PI*World->GridScale;
	conv = (World->GridScale*World->GridScale*World->GridScale)/COANGS;

	if(UTMPDouble==NULL)
		UTMPDouble = new double[gridSizeXYZ];
	if(TMPDouble == NULL)
	{
		TMPDouble = new double[gridSizeXYZ];
		for(j=0;j<gridSizeXYZ;j++)TMPDouble[j]=0.0f;
	}
	double *UTMP=UTMPDouble;
	double *TMP=TMPDouble;
	
	int MaxSingularNum=0,IonType2;
	for(IonType2=0;IonType2<World->NIonsTypes;IonType2++)
	{
		if(MaxSingularNum<SingularNum[IonType2][2])
			MaxSingularNum=SingularNum[IonType2][2];
	}
	for(i=0;i<6;i++)
		if(dphi[i] == NULL)
			if(!(dphi[i] = new float[MaxSingularNum]))
	{
		cerr << "ERROR 204: No memory available\n";
		exit(204);
	}
	//Set up the Potential energy
	if(World->PMF!=NULL)
	{
		if(potential!=NULL)
		{
			if(PMFWeight>=0.0)
			{
				//for(i=0;i<World->NIonsTypes;i++)
				for(j=0;j<gridSizeXYZ;j++)
					UTMP[j]=0.5*(World->IonsQ[IonType]*potential[j]+World->PMF[IonType][j]*PMFWeight);
			}
			else
			{
				//for(i=0;i<World->NIonsTypes;i++)
				for(j=0;j<gridSizeXYZ;j++)
					UTMP[j]=0.5*(World->IonsQ[IonType]*potential[j]+World->PMF[IonType][j]);
			}
		}
		else
		{
			if(PMFWeight>=0.0)
			{
				//for(i=0;i<World->NIonsTypes;i++)
				for(j=0;j<gridSizeXYZ;j++)
					UTMP[j]=0.5*(World->PMF[IonType][j])*PMFWeight;
			}
			else
			{
				//for(i=0;i<World->NIonsTypes;i++)
				for(j=0;j<gridSizeXYZ;j++)
					UTMP[j]=0.5*(World->PMF[IonType][j]);
			}
		}
	}
	else
		for(i=0;i<World->NIonsTypes;i++)
			for(j=0;j<gridSizeXYZ;j++)UTMP[j]=0.5*(World->IonsQ[IonType]*potential[j]);

	//for(IonType=0;IonType<World->NIonsTypes;IonType++)
	//{
	for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
		C[IonType][gridPoint]/=fpoh;

	for(i1=0;i1<NoSingularNum[IonType][2];i1++)
	{
		gridPoint=IndexNoSingular[IonType][i1];
		gridpx = gridPoint+1;
		gridmx = gridPoint-1;
		gridpy = gridPoint+gridSizeX;
		gridmy = gridPoint-gridSizeX;
		gridpz = gridPoint+gridSizeXY;
		gridmz = gridPoint-gridSizeXY;
	
		phit = UTMP[gridpx] + UTMP[gridmx] + UTMP[gridpy] + UTMP[gridmy] + UTMP[gridpz] + UTMP[gridmz] - 6.0*UTMP[gridPoint];
		TMP[gridPoint]=1.0/(6.0-phit);
	}
	for(i=0;i<SingularNum[IonType][2];i++)
	{
		gridPoint=IndexSingular[IonType][i];
		gridpx = gridPoint+1;
		gridmx = gridPoint-1;
		gridpy = gridPoint+gridSizeX;
		gridmy = gridPoint-gridSizeX;
		gridpz = gridPoint+gridSizeXY;
		gridmz = gridPoint-gridSizeXY;
		
		dphi[0][i] = dix[IonType][0][i]*UTMP[gridpx];
		dphi[1][i] = dix[IonType][1][i]*UTMP[gridmx];
		dphi[2][i] = dix[IonType][2][i]*UTMP[gridpy];
		dphi[3][i] = dix[IonType][3][i]*UTMP[gridmy];
		dphi[4][i] = dix[IonType][4][i]*UTMP[gridpz];
		dphi[5][i] = dix[IonType][5][i]*UTMP[gridmz];
		phit = dphi[0][i] + dphi[1][i] + dphi[2][i] + dphi[3][i] + dphi[4][i] + dphi[5][i] - dixt[IonType][i]*UTMP[gridPoint];
		TMP[gridPoint] = 1.0/(dixt[IonType][i]-phit);
	}
	//}
	return EXIT_SUCCESS;
}
int NernstPlankSolver::NernstPlanckSolverDiterDouble(int IonType,bool calcchange)
{
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	double gridScale;
	int iteration;
	int i,j,k,i1;
	int gridPoint;
	double om1,om2;
	double * potential;
	//float * diffusion;
	int gridSizeXY;
	int gridSizeXYZ;
	double change;
	double fpoh;
	//float * positiveCharge;
	//float * negativeCharge;
	double rodlt1,rodlt2;
	double phi0,phi1,phi2,phi3,phi4,phi5,phi6;
	double rophi1,rophi2;
	double term1,term2;
	double phit;
	double phit1,phit2;
	//float denom1,denom2;
	int gridpx,gridmx;
	int gridpy,gridmy;
	int gridpz,gridmz;
	double conv;
	bool *PeriodicBoundaryCondition;
	int itmp1,itmp2;
	double *C[2];
	//int IonType;
	
	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;	
	//diffusion = World->D[0];
	//potential = World->PotentialDouble;
	//positiveCharge = World->C[0];
	//negativeCharge = World->C[1];
	C[0] = World->CDouble[0];
	C[1] = World->CDouble[1];
	double *UTMP=UTMPDouble;
	double *TMP=TMPDouble;
	PeriodicBoundaryCondition=World->PBC;
	
	om1 = Relaxation;
	om2 = 1-om1;

	fpoh = 4*M_PI*World->GridScale;
	conv = (World->GridScale*World->GridScale*World->GridScale)/COANGS;
	
	for(iteration=1;iteration<=MaxIterations;iteration++)
	{
		if(calcchange&&((verbose&&(iteration%ConvergenceCheck==0))||(iteration==MaxIterations)))
		{
			MaxChange = 0;
			//for(IonType=0;IonType<World->NIonsTypes;IonType++)
			//{
			int countNegative=0;
			for(j=0;j<2;j++)
			{
				for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
				{
					gridPoint=IndexNoSingular[IonType][i];
					gridpx = gridPoint+1;
					gridmx = gridPoint-1;
					gridpy = gridPoint+gridSizeX;
					gridmy = gridPoint-gridSizeX;
					gridpz = gridPoint+gridSizeXY;
					gridmz = gridPoint-gridSizeXY;
					
					rodlt1 = C[IonType][gridpx] + C[IonType][gridmx] + C[IonType][gridpy] + C[IonType][gridmy] + C[IonType][gridpz] + C[IonType][gridmz];
					change = ( rodlt1 * ( 1.0 - UTMP[gridPoint] ) + C[IonType][gridpx]*UTMP[gridpx] + C[IonType][gridmx]*UTMP[gridmx] + C[IonType][gridpy]*UTMP[gridpy] + C[IonType][gridmy]*UTMP[gridmy] + C[IonType][gridpz]*UTMP[gridpz] + C[IonType][gridmz]*UTMP[gridmz]) * TMP[gridPoint];
					change -= C[IonType][gridPoint];
					change *= om1;
					MaxChange = MaxChange>change?MaxChange:change;
					MaxChange = -MaxChange<change?MaxChange:-change;
					C[IonType][gridPoint]+=change;
						
//						 if(C[IonType][gridPoint]<0.0)
//						 {
//							 fprintf(stderr,"C[%d][%d]=%g<0.0\n",IonType,gridPoint,[IonType][gridPoint]);
//							 return EXIT_FAILURE;
//						 }
				}
				for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++)
				{
					gridPoint=IndexSingular[IonType][i];
					gridpx = gridPoint+1;
					gridmx = gridPoint-1;
					gridpy = gridPoint+gridSizeX;
					gridmy = gridPoint-gridSizeX;
					gridpz = gridPoint+gridSizeXY;
					gridmz = gridPoint-gridSizeXY;
						
					rodlt1 = C[IonType][gridpx]*dix[IonType][0][i] + C[IonType][gridmx]*dix[IonType][1][i] + C[IonType][gridpy]*dix[IonType][2][i] + C[IonType][gridmy]*dix[IonType][3][i] + C[IonType][gridpz]*dix[IonType][4][i] + C[IonType][gridmz]*dix[IonType][5][i];
					rophi1 = C[IonType][gridpx]*dphi[0][i] + C[IonType][gridmx]*dphi[1][i] + C[IonType][gridpy]*dphi[2][i] + C[IonType][gridmy]*dphi[3][i] + C[IonType][gridpz]*dphi[4][i] + C[IonType][gridmz]*dphi[5][i];
					term1 = rodlt1*(1-UTMP[gridPoint])+rophi1;
					change = om1*(term1*TMP[gridPoint]-C[IonType][gridPoint]);
					MaxChange = MaxChange>change?MaxChange:change;
					MaxChange = -MaxChange<change?MaxChange:-change;
					C[IonType][gridPoint]+=change;
				}
				if(bConcCorrection)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						if(C[IonType][gridPoint]<0.0)
						{
							C[IonType][gridPoint]=0.0;
							if(countNegative<=10)
								PntsWithNegativeConc[countNegative]=gridPoint;
							countNegative++;
						}
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexSingular[IonType][i];
						if(C[IonType][gridPoint]<0.0)
						{
							C[IonType][gridPoint]=0.0;
							if(countNegative<=10)
								PntsWithNegativeConc[countNegative]=gridPoint;
							countNegative++;
						}
					}
				}
				World->BorderExchangeDouble(C[IonType]);
			}
			if(bConcCorrection&&countNegative>0)
			{
				pnpWarning("IonType %d have %d points with value<0.0\n",IonType,countNegative);
				if(countNegative<=10)
				{
					pnpWarning("\tpnts = ");
					for(i=0;i<countNegative;i++)
						pnpPrint(" %d",PntsWithNegativeConc[i]);
					pnpPrint("\n");
				}
			}
			//}
			MaxChange*=conv;
#ifdef MPI_PARALLEL
			int dest;
			pnpsapp->MyComGroup.Barrier();
			
			if(World->MyRank==0)
			{
				for(dest=1;dest<World->NProcs;dest++)
				{
					pnpsapp->MyComGroup.Recv(&change, 1, MPI::FLOAT, dest, 0);
					MaxChange = MaxChange>change?MaxChange:change;
					MaxChange = -MaxChange<change?MaxChange:-change;
				}
			}
			else
			{
				pnpsapp->MyComGroup.Send(&MaxChange, 1, MPI::FLOAT, 0, 0);
			}
#endif
			if(verbose)if(World->MyRank==0)
			{
				fprintf(stdout, "<NernstPlankIteration Nit=\"%6d\"	dC=\"%16.8g\"/>\n", iteration, MaxChange);
			}
#if defined(_MSC_VER) || defined(__DECCXX) || (__GNUC__ < 3)
			if(MaxChange>1E13)
#else
				if(MaxChange>1E13||MaxChange==NAN)
#endif
			{
				pnpError("MaxChange=%f will finish\n",MaxChange);
				return EXIT_FAILURE;
			}
			if(MaxChange<Convergence)
				iteration = MaxIterations+1;
		}
		else
		{
			//for(IonType=0;IonType<World->NIonsTypes;IonType++)
			//{
			for(j=0;j<2;j++)
			{
				for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
				{
					gridPoint=IndexNoSingular[IonType][i];
					gridpx = gridPoint+1;
					gridmx = gridPoint-1;
					gridpy = gridPoint+gridSizeX;
					gridmy = gridPoint-gridSizeX;
					gridpz = gridPoint+gridSizeXY;
					gridmz = gridPoint-gridSizeXY;
					rodlt1 = C[IonType][gridpx] + C[IonType][gridmx] + C[IonType][gridpy] + C[IonType][gridmy] + C[IonType][gridpz] + C[IonType][gridmz];
					change = ( rodlt1 * ( 1.0 - UTMP[gridPoint] ) + C[IonType][gridpx]*UTMP[gridpx] + C[IonType][gridmx]*UTMP[gridmx] + C[IonType][gridpy]*UTMP[gridpy] + C[IonType][gridmy]*UTMP[gridmy] + C[IonType][gridpz]*UTMP[gridpz] + C[IonType][gridmz]*UTMP[gridmz]) * TMP[gridPoint];
					change -= C[IonType][gridPoint];
					change *= om1;
					C[IonType][gridPoint]+=change;
				}
				for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++){
					gridPoint=IndexSingular[IonType][i];
					gridpx = gridPoint+1;
					gridmx = gridPoint-1;
					gridpy = gridPoint+gridSizeX;
					gridmy = gridPoint-gridSizeX;
					gridpz = gridPoint+gridSizeXY;
					gridmz = gridPoint-gridSizeXY;
	
					rodlt1 = C[IonType][gridpx]*dix[IonType][0][i] + C[IonType][gridmx]*dix[IonType][1][i] + C[IonType][gridpy]*dix[IonType][2][i] + C[IonType][gridmy]*dix[IonType][3][i] + C[IonType][gridpz]*dix[IonType][4][i] + C[IonType][gridmz]*dix[IonType][5][i];
					rophi1 = C[IonType][gridpx]*dphi[0][i] + C[IonType][gridmx]*dphi[1][i] + C[IonType][gridpy]*dphi[2][i] + C[IonType][gridmy]*dphi[3][i] + C[IonType][gridpz]*dphi[4][i] + C[IonType][gridmz]*dphi[5][i];
					term1 = rodlt1*(1-UTMP[gridPoint])+rophi1;
					change = om1*(term1*TMP[gridPoint]-C[IonType][gridPoint]);
					C[IonType][gridPoint]+=change;
				}
				if(bConcCorrection)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						if(C[IonType][gridPoint]<0.0)
						{
							C[IonType][gridPoint]=0.0;
						}
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++){
						gridPoint=IndexSingular[IonType][i];
						if(C[IonType][gridPoint]<0.0)
						{
							C[IonType][gridPoint]=0.0;
						}
					}
				}
				World->BorderExchangeDouble(C[IonType]);
			}
			//}
		}
	}
	return EXIT_SUCCESS;
}
int NernstPlankSolver::NernstPlanckSolverDpostDouble(int IonType)
{
	int i1,gridPoint;
	int gridSizeXYZ=World->GS_XYZ;
	float fpoh = 4*M_PI*World->GridScale;
	
	//for(i1=0;i1<World->NIonsTypes;i1++)
	for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
		World->CDouble[IonType][gridPoint]*=fpoh;
	return EXIT_SUCCESS;
}

int NernstPlankSolver::NernstPlanckSolverDOld()
{
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	float gridScale;
	int iteration;
	int i,j,k,i1;
	int gridPoint;
	float om1,om2;
	float * potential;
	//float * diffusion;
	int gridSizeXY;
	int gridSizeXYZ;
	float change;
	float fpoh;
	//float * positiveCharge;
	//float * negativeCharge;
	float rodlt1,rodlt2;
	float phi0,phi1,phi2,phi3,phi4,phi5,phi6;
	float rophi1,rophi2;
	float term1,term2;
	float phit;
	float phit1,phit2;
	//float denom1,denom2;
	int gridpx,gridmx;
	int gridpy,gridmy;
	int gridpz,gridmz;
	float conv;
	bool *PeriodicBoundaryCondition;
	int itmp1,itmp2;
	float *C[2];
	int IonType;
	
	//CheckSystem();

	if(!World){
		cerr << "ERROR 311: No World\n";
		exit(104);
	}
	
	if(verbose&&World->MyRank==0)
	{
		//fprintf(stdout,"\nNernst-Plank Solver\n");
	}
	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;	
	//diffusion = World->D[0];
	potential = World->Potential;
	//positiveCharge = World->C[0];
	//negativeCharge = World->C[1];
	C[0] = World->C[0];
	C[1] = World->C[1];
	PeriodicBoundaryCondition=World->PBC;
	
	om1 = Relaxation;
	om2 = 1-om1;

	fpoh = 4*M_PI*World->GridScale;
	conv = (World->GridScale*World->GridScale*World->GridScale)/COANGS;

	
	if(UTMPSingle==NULL)
		UTMPSingle = new float[gridSizeXYZ];
	if(TMPSingle == NULL)
	{
		TMPSingle = new float[gridSizeXYZ];
		for(j=0;j<gridSizeXYZ;j++)TMPSingle[j]=0.0f;
	}
	float *UTMP=UTMPSingle;
	float *TMP=TMPSingle;
	
	for(IonType=0;IonType<World->NIonsTypes;IonType++)
	{
		int MaxSingularNum=0,IonType2;
		for(IonType2=0;IonType2<World->NIonsTypes;IonType2++)
		{
			if(MaxSingularNum<SingularNum[IonType2][2])
				MaxSingularNum=SingularNum[IonType2][2];
		}
		for(i=0;i<6;i++)
			if(dphi[i] == NULL)
				if(!(dphi[i] = new float[MaxSingularNum]))
		{
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		//Set up the Potential energy
		if(World->PMF!=NULL)
		{
			if(potential!=NULL)
			{
				if(PMFWeight>=0.0)
				{
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[j]=0.5*(World->IonsQ[i]*potential[j]+World->PMF[i][j]*PMFWeight);
				}
				else
				{
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[j]=0.5*(World->IonsQ[i]*potential[j]+World->PMF[i][j]);
				}
			}
			else
			{
				if(PMFWeight>=0.0)
				{
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[j]=0.5*(World->PMF[i][j])*PMFWeight;
				}
				else
				{
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[j]=0.5*(World->PMF[i][j]);
				}
			}
		}
		else
			for(j=0;j<gridSizeXYZ;j++)UTMP[j]=0.5*(World->IonsQ[i]*potential[j]);

	
		for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)C[IonType][gridPoint]/=fpoh;

		for(i1=0;i1<NoSingularNum[IonType][2];i1++)
		{
			gridPoint=IndexNoSingular[IonType][i1];
			gridpx = gridPoint+1;
			gridmx = gridPoint-1;
			gridpy = gridPoint+gridSizeX;
			gridmy = gridPoint-gridSizeX;
			gridpz = gridPoint+gridSizeXY;
			gridmz = gridPoint-gridSizeXY;
	
			phit = UTMP[gridpx] + UTMP[gridmx] + UTMP[gridpy] + UTMP[gridmy] + UTMP[gridpz] + UTMP[gridmz] - 6.0*UTMP[gridPoint];
			TMP[gridPoint]=1.0/(6.0-phit);
		}
		for(i=0;i<SingularNum[IonType][2];i++)
		{
			gridPoint=IndexSingular[IonType][i];
			gridpx = gridPoint+1;
			gridmx = gridPoint-1;
			gridpy = gridPoint+gridSizeX;
			gridmy = gridPoint-gridSizeX;
			gridpz = gridPoint+gridSizeXY;
			gridmz = gridPoint-gridSizeXY;
		
			dphi[0][i] = dix[IonType][0][i]*UTMP[gridpx];
			dphi[1][i] = dix[IonType][1][i]*UTMP[gridmx];
			dphi[2][i] = dix[IonType][2][i]*UTMP[gridpy];
			dphi[3][i] = dix[IonType][3][i]*UTMP[gridmy];
			dphi[4][i] = dix[IonType][4][i]*UTMP[gridpz];
			dphi[5][i] = dix[IonType][5][i]*UTMP[gridmz];
			phit = dphi[0][i] + dphi[1][i] + dphi[2][i] + dphi[3][i] + dphi[4][i] + dphi[5][i] - dixt[IonType][i]*UTMP[gridPoint];
			TMP[gridPoint] = 1.0/(dixt[IonType][i]-phit);
		}
		
		for(iteration=1;iteration<=MaxIterations;iteration++)
		{
			if((verbose&&(iteration%ConvergenceCheck==0))||(iteration==MaxIterations))
			{
				MaxChange = 0;
				for(j=0;j<2;j++)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
					
						rodlt1 = C[IonType][gridpx] + C[IonType][gridmx] + C[IonType][gridpy] + C[IonType][gridmy] + C[IonType][gridpz] + C[IonType][gridmz];
						change = ( rodlt1 * ( 1.0 - UTMP[gridPoint] ) + C[IonType][gridpx]*UTMP[gridpx] + C[IonType][gridmx]*UTMP[gridmx] + C[IonType][gridpy]*UTMP[gridpy] + C[IonType][gridmy]*UTMP[gridmy] + C[IonType][gridpz]*UTMP[gridpz] + C[IonType][gridmz]*UTMP[gridmz]) * TMP[gridPoint];
						change -= C[IonType][gridPoint];
						change *= om1;
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						C[IonType][gridPoint]+=change;
						
//						 if(C[IonType][gridPoint]<0.0)
//						 {
//							 fprintf(stderr,"C[%d][%d]=%g<0.0\n",IonType,gridPoint,[IonType][gridPoint]);
//							 return EXIT_FAILURE;
//						 }
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
						
						rodlt1 = C[IonType][gridpx]*dix[IonType][0][i] + C[IonType][gridmx]*dix[IonType][1][i] + C[IonType][gridpy]*dix[IonType][2][i] + C[IonType][gridmy]*dix[IonType][3][i] + C[IonType][gridpz]*dix[IonType][4][i] + C[IonType][gridmz]*dix[IonType][5][i];
						rophi1 = C[IonType][gridpx]*dphi[0][i] + C[IonType][gridmx]*dphi[1][i] + C[IonType][gridpy]*dphi[2][i] + C[IonType][gridmy]*dphi[3][i] + C[IonType][gridpz]*dphi[4][i] + C[IonType][gridmz]*dphi[5][i];
						term1 = rodlt1*(1-UTMP[gridPoint])+rophi1;
						change = om1*(term1*TMP[gridPoint]-C[IonType][gridPoint]);
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						C[IonType][gridPoint]+=change;
					}
					World->BorderExchange(C[IonType]);
				}
				MaxChange*=conv;
#ifdef MPI_PARALLEL
				int dest;
				pnpsapp->MyComGroup.Barrier();
				
				if(World->MyRank==0)
				{
					for(dest=1;dest<World->NProcs;dest++)
					{
						pnpsapp->MyComGroup.Recv(&change, 1, MPI::FLOAT, dest, 0);
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
					}
				}
				else
				{
					pnpsapp->MyComGroup.Send(&MaxChange, 1, MPI::FLOAT, 0, 0);
				}
#endif
				if(verbose)if(World->MyRank==0)
				{
					fprintf(stdout, "<NernstPlankIteration Nit=\"%6d\"	dC=\"%16.8g\" ion=\"%d\"/>\n", iteration, MaxChange,IonType);
				}
#if defined(_MSC_VER) || defined(__DECCXX) || (__GNUC__ < 3)
				if(MaxChange>1E13)
#else
				if(MaxChange>1E13||MaxChange==NAN)
#endif
				{
					for(i1=0;i1<World->NIonsTypes;i1++)for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
							C[i1][gridPoint]*=fpoh;
					return EXIT_FAILURE;
				}
				if(MaxChange<Convergence)
					iteration = MaxIterations+1;
			}
			else
			{
				for(j=0;j<2;j++)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
						rodlt1 = C[IonType][gridpx] + C[IonType][gridmx] + C[IonType][gridpy] + C[IonType][gridmy] + C[IonType][gridpz] + C[IonType][gridmz];
						change = ( rodlt1 * ( 1.0 - UTMP[gridPoint] ) + C[IonType][gridpx]*UTMP[gridpx] + C[IonType][gridmx]*UTMP[gridmx] + C[IonType][gridpy]*UTMP[gridpy] + C[IonType][gridmy]*UTMP[gridmy] + C[IonType][gridpz]*UTMP[gridpz] + C[IonType][gridmz]*UTMP[gridmz]) * TMP[gridPoint];
						change -= C[IonType][gridPoint];
						change *= om1;
						C[IonType][gridPoint]+=change;
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++){
						gridPoint=IndexSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
	
						rodlt1 = C[IonType][gridpx]*dix[IonType][0][i] + C[IonType][gridmx]*dix[IonType][1][i] + C[IonType][gridpy]*dix[IonType][2][i] + C[IonType][gridmy]*dix[IonType][3][i] + C[IonType][gridpz]*dix[IonType][4][i] + C[IonType][gridmz]*dix[IonType][5][i];
						rophi1 = C[IonType][gridpx]*dphi[0][i] + C[IonType][gridmx]*dphi[1][i] + C[IonType][gridpy]*dphi[2][i] + C[IonType][gridmy]*dphi[3][i] + C[IonType][gridpz]*dphi[4][i] + C[IonType][gridmz]*dphi[5][i];
						term1 = rodlt1*(1-UTMP[gridPoint])+rophi1;
						change = om1*(term1*TMP[gridPoint]-C[IonType][gridPoint]);
						C[IonType][gridPoint]+=change;
					}
					World->BorderExchange(C[IonType]);
				}
			}
		}
		for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
			C[IonType][gridPoint]*=fpoh;
	}
	return EXIT_SUCCESS;
}
int NernstPlankSolver::NernstPlanckSolverN2tmp()
{
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	float gridScale;
	int iteration;
	int i,j,k,i1;
	int gridPoint;
	float om1,om2;
	float * potential;
	//float * diffusion;
	int gridSizeXY;
	int gridSizeXYZ;
	float change;
	float fpoh;
	//float * positiveCharge;
	//float * negativeCharge;
	float rodlt1,rodlt2;
	float phi0,phi1,phi2,phi3,phi4,phi5,phi6;
	float rophi1,rophi2;
	float term1,term2;
	float phit;
	float phit1,phit2;
	//float denom1,denom2;
	int gridpx,gridmx;
	int gridpy,gridmy;
	int gridpz,gridmz;
	float conv;
	bool *PeriodicBoundaryCondition;
	int itmp1,itmp2;
	float *C[2];
	int IonType;
	float _MaxChange=0.0;

	
	if(!World){
		cerr << "ERROR 311: No World\n";
		exit(104);
	}
	
	if(verbose&&World->MyRank==0)
	{
		//fprintf(stdout,"\nNernst-Plank Solver\n");
	}
	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;	
	//diffusion = World->D[0];
	potential = World->Potential;
	//positiveCharge = World->C[0];
	//negativeCharge = World->C[1];
	C[0] = World->C[0];
	C[1] = World->C[1];
	PeriodicBoundaryCondition=World->PBC;
	
	om1 = Relaxation;
	om2 = 1-om1;

	fpoh = 4*M_PI*World->GridScale;
	conv = (World->GridScale*World->GridScale*World->GridScale)/COANGS;

	
	if(UTMPSingle==NULL)
		UTMPSingle = new float[gridSizeXYZ];
	if(TMPSingle == NULL)
	{
		TMPSingle = new float[gridSizeXYZ];
		for(j=0;j<gridSizeXYZ;j++)TMPSingle[j]=0.0f;
	}
	float *UTMP=UTMPSingle;
	float *TMP=TMPSingle;
	
	float *utmp=UTMP;
	float *tmp=TMP;
	float *Ci;
	float *dphi0,*dphi1,*dphi2,*dphi3,*dphi4,*dphi5;
	
	//Set up the Potential energy
	for(IonType=0;IonType<World->NIonsTypes;IonType++)
	{
		int MaxSingularNum=0,IonType2;
		for(IonType2=0;IonType2<World->NIonsTypes;IonType2++)
		{
			if(MaxSingularNum<SingularNum[IonType2][2])
				MaxSingularNum=SingularNum[IonType2][2];
		}
		for(i=0;i<6;i++)
			if(dphi[i] == NULL)
				if(!(dphi[i] = new float[MaxSingularNum]))
		{
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		dphi0=dphi[0];
		dphi1=dphi[1];
		dphi2=dphi[2];
		dphi3=dphi[3];
		dphi4=dphi[4];
		dphi5=dphi[5];
		Ci=C[IonType];
		
		if(World->PMF!=NULL)
			for(j=0;j<gridSizeXYZ;j++)utmp[j]=0.5*(World->IonsQ[IonType]*potential[j]+World->PMF[IonType][j]);
		else
			for(j=0;j<gridSizeXYZ;j++)utmp[j]=0.5*(World->IonsQ[IonType]*potential[j]);

	
		for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)C[IonType][gridPoint]/=fpoh;

		for(i1=0;i1<NoSingularNum[IonType][2];i1++)
		{
			gridPoint=IndexNoSingular[IonType][i1];
			gridpx = gridPoint+1;
			gridmx = gridPoint-1;
			gridpy = gridPoint+gridSizeX;
			gridmy = gridPoint-gridSizeX;
			gridpz = gridPoint+gridSizeXY;
			gridmz = gridPoint-gridSizeXY;
	
			phit = utmp[gridpx] + utmp[gridmx] + utmp[gridpy] + utmp[gridmy] + utmp[gridpz] + utmp[gridmz] - 6.0*utmp[gridPoint];
			tmp[gridPoint]=1.0/(6.0-phit);
		}
		for(i=0;i<SingularNum[IonType][2];i++)
		{
			gridPoint=IndexSingular[IonType][i];
			gridpx = gridPoint+1;
			gridmx = gridPoint-1;
			gridpy = gridPoint+gridSizeX;
			gridmy = gridPoint-gridSizeX;
			gridpz = gridPoint+gridSizeXY;
			gridmz = gridPoint-gridSizeXY;
		
			dphi0[i] = dix[IonType][0][i]*utmp[gridpx];
			dphi1[i] = dix[IonType][1][i]*utmp[gridmx];
			dphi2[i] = dix[IonType][2][i]*utmp[gridpy];
			dphi3[i] = dix[IonType][3][i]*utmp[gridmy];
			dphi4[i] = dix[IonType][4][i]*utmp[gridpz];
			dphi5[i] = dix[IonType][5][i]*utmp[gridmz];
			phit = dphi0[i] + dphi1[i] + dphi2[i] + dphi3[i] + dphi4[i] + dphi5[i] - dixt[IonType][i]*utmp[gridPoint];
			tmp[gridPoint] = 1.0/(dixt[IonType][i]-phit);
		}
		for(iteration=1;iteration<=MaxIterations;iteration++)
		{
			if((verbose&&(iteration%ConvergenceCheck==0))||(iteration==MaxIterations))
			{
				MaxChange = 0;
				for(j=0;j<2;j++)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
					
						rodlt1 = Ci[gridpx] + Ci[gridmx] + Ci[gridpy] + Ci[gridmy] + Ci[gridpz] + Ci[gridmz];
						change = ( rodlt1 * ( 1.0 - utmp[gridPoint] ) + Ci[gridpx]*utmp[gridpx] + Ci[gridmx]*utmp[gridmx] + Ci[gridpy]*utmp[gridpy] + Ci[gridmy]*utmp[gridmy] + Ci[gridpz]*utmp[gridpz] + Ci[gridmz]*utmp[gridmz]) * tmp[gridPoint];
						change -= Ci[gridPoint];
						change *= om1;
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						Ci[gridPoint]+=change;
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
						
						rodlt1 = Ci[gridpx]*dix[IonType][0][i] + Ci[gridmx]*dix[IonType][1][i] + Ci[gridpy]*dix[IonType][2][i] + Ci[gridmy]*dix[IonType][3][i] + Ci[gridpz]*dix[IonType][4][i] + Ci[gridmz]*dix[IonType][5][i];
						rophi1 = Ci[gridpx]*dphi0[i] + Ci[gridmx]*dphi1[i] + Ci[gridpy]*dphi2[i] + Ci[gridmy]*dphi3[i] + Ci[gridpz]*dphi4[i] + Ci[gridmz]*dphi5[i];
						term1 = rodlt1*(1-utmp[gridPoint])+rophi1;
						change = om1*(term1*tmp[gridPoint]-Ci[gridPoint]);
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						Ci[gridPoint]+=change;
					}
					World->BorderExchange(Ci);
				}
				MaxChange*=conv;
				if(verbose||iteration==MaxIterations)
				{
					if(World->NProcs!=1){
					//int tag=SEND_ENERGY;
					/*MPI::COMM_WORLD.Barrier();
						MPI::Status	status;
						float MaxChangeProc;
						if(World->MyRank==0){
						for(i=1;i<World->NProcs;i++){
						MPI::COMM_WORLD.Recv(&MaxChangeProc, 1, MPI::FLOAT, i, SEND_MAXCHANGE, status);
						MaxChange = MaxChange>MaxChangeProc?MaxChange:MaxChangeProc;
					}
					}
						else{
						MPI::COMM_WORLD.Send(&MaxChange, 1, MPI::FLOAT, 0, SEND_MAXCHANGE);
					}
						MPI::COMM_WORLD.Barrier();*/
					}
					if(verbose)if(World->MyRank==0){
						fprintf(stdout, "<NernstPlankIteration Nit=\"%6d\"	dC[%d]=\"%16.8g\"/>\n", iteration, IonType,MaxChange);
		//MaxChange=0;
					}
#if defined(_MSC_VER) || defined(__DECCXX) || (__GNUC__ < 3)
					if(MaxChange>1E13)
#else
						if(MaxChange>1E13||MaxChange==NAN)
#endif
					{
						for(i1=0;i1<World->NIonsTypes;i1++)
							for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
								C[i1][gridPoint]*=fpoh;
						return EXIT_FAILURE;
					}
				}
	
				if(MaxChange<Convergence)
					iteration = MaxIterations+1;
			}
			else
			{
				for(j=0;j<2;j++)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
						rodlt1 = Ci[gridpx] + Ci[gridmx] + Ci[gridpy] + Ci[gridmy] + Ci[gridpz] + Ci[gridmz];
						change = ( rodlt1 * ( 1.0 - utmp[gridPoint] ) + Ci[gridpx]*utmp[gridpx] + Ci[gridmx]*utmp[gridmx] + Ci[gridpy]*utmp[gridpy] + Ci[gridmy]*utmp[gridmy] + Ci[gridpz]*utmp[gridpz] + Ci[gridmz]*utmp[gridmz]) * tmp[gridPoint];
						change -= Ci[gridPoint];
						change *= om1;
						Ci[gridPoint]+=change;
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++){
						gridPoint=IndexSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
	
						rodlt1 = Ci[gridpx]*dix[IonType][0][i] + Ci[gridmx]*dix[IonType][1][i] + Ci[gridpy]*dix[IonType][2][i] + Ci[gridmy]*dix[IonType][3][i] + Ci[gridpz]*dix[IonType][4][i] + Ci[gridmz]*dix[IonType][5][i];
						rophi1 = Ci[gridpx]*dphi0[i] + Ci[gridmx]*dphi1[i] + Ci[gridpy]*dphi2[i] + Ci[gridmy]*dphi3[i] + Ci[gridpz]*dphi4[i] + Ci[gridmz]*dphi5[i];
						term1 = rodlt1*(1-utmp[gridPoint])+rophi1;
						change = om1*(term1*tmp[gridPoint]-Ci[gridPoint]);
						Ci[gridPoint]+=change;
					}
				}
				World->BorderExchange(Ci);
			}
		}
		if(MaxChange>_MaxChange)_MaxChange=MaxChange;
	}
	MaxChange=_MaxChange;
	for(i1=0;i1<World->NIonsTypes;i1++)
		for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
			C[i1][gridPoint]*=fpoh;

	return EXIT_SUCCESS;
}
VectorField3D* NernstPlankSolver::CalcJ()
{
	int * GS;
	int GS_X;
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
	float *UTMP=UTMPSingle;
	float *TMP=TMPSingle;
	
	NodeIndex* NIndex=World->NIndexing->NIndex;
	NodeIndexing* NIndexing=World->NIndexing;
	NodeIndex DiffBoarderMask=NodeIndexing::DiffIon0BoarderMask;
	NodeIndex BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;
	
	GS = World->GridSize;
	GS_X = GS[0];
	GS_XY = GS[0]*GS[1];
	GS_XYZ = GS[2]*GS_XY;
	gridScale = World->GridScale;
	
	C[0] = World->C[0];
	C[1] = World->C[1];
	potential = World->Potential;
	
	VectorField3D* J=new VectorField3D(GS,gridScale,8);
	float **V=J->V;
	for(j=0;j<World->NIonsTypes;j++)
	{
		if(World->PMF!=NULL)
			for(i=0;i<GS_XYZ;i++)
				UTMP[i]=0.5*(World->IonsQ[j]*potential[i]+World->PMF[j][i]);
		else
			for(i=0;i<GS_XYZ;i++)
				UTMP[i]=0.5*World->IonsQ[j]*potential[i];
		
		fpoh = 4*M_PI*gridScale;
		flcon = 1.6*gridScale*gridScale*1000;
		ftot = flcon/fpoh;
		
		float d,dc,sumc;
		
		for(k=0;k<3;k++)
		{
			switch(k)
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
			ftot = -1.0*World->IonsQ[j]*flcon/fpoh;
			for(i=0;i<NoSingularNum[j][2];i++)
			{
				gridPoint=IndexNoSingular[j][i];
				gridp = gridPoint+gridInc;
				pot = UTMP[gridp]-UTMP[gridPoint];
				d=NIndexing->GetDiffFloat(j,gridPoint);
				dc=C[j][gridp]-C[j][gridPoint];
				sumc=(C[j][gridp]+C[j][gridPoint]);
				V[k+j*4][gridPoint]=ftot*d*(dc+pot*sumc);
			
			}
			for(i=0;i<SingularNum[j][2];i++)
			{
				gridPoint=IndexSingular[j][i];
				gridp = gridPoint+gridInc;
				pot = UTMP[gridp]-UTMP[gridPoint];
				d=0.5*(NIndexing->GetDiffFloat(j,gridPoint)+NIndexing->GetDiffFloat(j,gridp));
				dc=C[j][gridp]-C[j][gridPoint];
				sumc=(C[j][gridp]+C[j][gridPoint]);
				V[k+j*4][gridPoint]=ftot*d*(dc+pot*sumc);
			//V[k+j*4][gridPoint]=ftot*dix[j][2*k][i]*(C[j][gridp]-C[j][gridPoint]+pot*(C[j][gridp]+C[j][gridPoint]));
			}
		}
		for(gridPoint=0;gridPoint<GS_XYZ;gridPoint++)
		{
			V[3+j*4][gridPoint]=sqrt(V[j*4][gridPoint]*V[j*4][gridPoint]+V[1+j*4][gridPoint]*V[1+j*4][gridPoint]+V[2+j*4][gridPoint]*V[2+j*4][gridPoint]);
		}
	}
	return J;
}

void NernstPlankSolver::nernstPlanckSolverCurrentProfile(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile)
{
	if(solver==Auto)
	{
		if(World->NIndexing==NULL)
			nernstPlanckSolverCurrentProfileW(dim, positiveCurrentProfile, negativeCurrentProfile);
		else if(World->D!=NULL)
			nernstPlanckSolverCurrentProfileDflex(dim, positiveCurrentProfile, negativeCurrentProfile);
		else nernstPlanckSolverCurrentProfileN(dim, positiveCurrentProfile, negativeCurrentProfile);
		
	}
	else
	{
		if(solver==NodeIndexBased)
			nernstPlanckSolverCurrentProfileN(dim, positiveCurrentProfile, negativeCurrentProfile);
		else if(solver==ArrayDirect)
			nernstPlanckSolverCurrentProfileAD(dim, positiveCurrentProfile, negativeCurrentProfile);
		else if(solver==PNPC)
			nernstPlanckSolverCurrentProfileW(dim, positiveCurrentProfile, negativeCurrentProfile);
		else
			nernstPlanckSolverCurrentProfileN(dim, positiveCurrentProfile, negativeCurrentProfile);
	}
	
}
void NernstPlankSolver::nernstPlanckSolverCurrentProfileW(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile)
{
	int * gridSize;
	long gridSizeXY;
	long gridSizeXYZ;
	long gridPoint;
	float gridScale;
	float fpoh;
	float flcon;
	float * potential;
	float * positive;
	float * negative;
	long * start;
	long * end;
	long * borderPoints;
	float * diffusion;
	float * dix[6];
	int i;
	long gridp,gridn;
	float potp,potn;
	int grid[3];
	long gridInc;
	long count;
	long borderCount;
	long nextBorderPoint;
	float delphiExp;
	int expInterpolation;
	
	gridSize = World->GridSize;
	gridSizeXY = gridSize[0]*gridSize[1];
	gridSizeXYZ = gridSize[2]*gridSizeXY;
	gridScale = World->GridScale;
	positive = World->C[0];
	negative = World->C[1];
	diffusion = World->D[0];
	potential = World->Potential;
	
	fpoh = 4*M_PI*gridScale;
	flcon = 1.6*gridScale*gridScale*1000;

	start = nernstPlanckSolverData->start;
	end = nernstPlanckSolverData->end;
	expInterpolation = nernstPlanckSolverData->expInterpolation;
	borderPoints = nernstPlanckSolverData->borderPoints;
	for(i=0;i<6;i++)
		dix[i] = nernstPlanckSolverData->dix[i];

	if(dim==0)
		gridInc = 1;
	else if(dim==1)
		gridInc = gridSize[0];
	else
		gridInc = gridSizeXY;
	
	for(i=0;i<gridSize[dim]-1;i++) {
		positiveCurrentProfile[i] = 0;
		negativeCurrentProfile[i] = 0;
	}

	for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++) {
		potential[gridPoint]*=0.5;
		positive[gridPoint]/=fpoh;
		negative[gridPoint]/=fpoh;
	}

	count=0;
	borderCount=0;
	nextBorderPoint = borderPoints[borderCount];

	while(end[count]<gridSizeXYZ) {
		for(gridPoint=start[count];gridPoint<end[count];gridPoint++) {
			grid[2] = gridPoint/gridSizeXY;
			grid[1] = (gridPoint%gridSizeXY)/gridSize[0];
			grid[0] = (gridPoint%gridSizeXY)%gridSize[0];
			gridp = gridPoint+gridInc;
			gridn = gridPoint-gridInc;
			potp = potential[gridp]-potential[gridPoint];
			potn = potential[gridPoint]-potential[gridn];
			if(gridPoint<nextBorderPoint) {
				if(grid[dim]==gridSize[dim]-2) {
					if(expInterpolation && fabs(potential[gridp]-potential[gridPoint])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
						delphiExp = exp(potp);
						positiveCurrentProfile[grid[dim]]+=diffusion[gridPoint]*(2*potp/(delphiExp-1/delphiExp))*(positive[gridp]-positive[gridPoint]+potp*(positive[gridp]/(0.5*(1+1/delphiExp))+positive[gridPoint]/(0.5*(1+delphiExp))));
						negativeCurrentProfile[grid[dim]]+=diffusion[gridPoint]*(2*potp/(delphiExp-1/delphiExp))*(negative[gridp]-negative[gridPoint]-potp*(negative[gridp]/(0.5*(1+delphiExp))+negative[gridPoint]/(0.5*(1+1/delphiExp))));
					}
					else {
						positiveCurrentProfile[grid[dim]]+=diffusion[gridPoint]*(positive[gridp]-positive[gridPoint]+potp*(positive[gridp]+positive[gridPoint]));
						negativeCurrentProfile[grid[dim]]+=diffusion[gridPoint]*(negative[gridp]-negative[gridPoint]-potp*(negative[gridp]+negative[gridPoint]));
					}
				}
				if(expInterpolation && fabs(potential[gridPoint]-potential[gridn])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
					delphiExp = exp(potn);
					positiveCurrentProfile[grid[dim]-1]+=diffusion[gridPoint]*(2*potn/(delphiExp-1/delphiExp))*(positive[gridPoint]-positive[gridn]+potn*(positive[gridPoint]/(0.5*(1+1/delphiExp))+positive[gridn]/(0.5*(1+delphiExp))));
					negativeCurrentProfile[grid[dim]-1]+=diffusion[gridPoint]*(2*potn/(delphiExp-1/delphiExp))*(negative[gridPoint]-negative[gridn]-potn*(negative[gridPoint]/(0.5*(1+delphiExp))+negative[gridn]/(0.5*(1+1/delphiExp))));
				}
				else {
					positiveCurrentProfile[grid[dim]-1]+=diffusion[gridPoint]*(positive[gridPoint]-positive[gridn]+potn*(positive[gridPoint]+positive[gridn]));
					negativeCurrentProfile[grid[dim]-1]+=diffusion[gridPoint]*(negative[gridPoint]-negative[gridn]-potn*(negative[gridPoint]+negative[gridn]));
				}
			}
			else {
				if(grid[dim]==gridSize[dim]-2) {
					if(expInterpolation && fabs(potential[gridp]-potential[gridPoint])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
						delphiExp = exp(potp);
						positiveCurrentProfile[grid[dim]]+=dix[2*dim][borderCount]*(2*potp/(delphiExp-1/delphiExp))*(positive[gridp]-positive[gridPoint]+potp*(positive[gridp]/(0.5*(1+1/delphiExp))+positive[gridPoint]/(0.5*(1+delphiExp))));
						negativeCurrentProfile[grid[dim]]+=dix[2*dim][borderCount]*(2*potp/(delphiExp-1/delphiExp))*(negative[gridp]-negative[gridPoint]-potp*(negative[gridp]/(0.5*(1+delphiExp))+negative[gridPoint]/(0.5*(1+1/delphiExp))));
					}
					else {
						positiveCurrentProfile[grid[dim]]+=dix[2*dim][borderCount]*(positive[gridp]-positive[gridPoint]+potp*(positive[gridp]+positive[gridPoint]));
						negativeCurrentProfile[grid[dim]]+=dix[2*dim][borderCount]*(negative[gridp]-negative[gridPoint]-potp*(negative[gridp]+negative[gridPoint]));
					}
				}
				if(expInterpolation && fabs(potential[gridPoint]-potential[gridn])>NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE) {
					delphiExp = exp(potn);
					positiveCurrentProfile[grid[dim]-1]+=dix[2*dim+1][borderCount]*(2*potn/(delphiExp-1/delphiExp))*(positive[gridPoint]-positive[gridn]+potn*(positive[gridPoint]/(0.5*(1+1/delphiExp))+positive[gridn]/(0.5*(1+delphiExp))));
					negativeCurrentProfile[grid[dim]-1]+=dix[2*dim+1][borderCount]*(2*potn/(delphiExp-1/delphiExp))*(negative[gridPoint]-negative[gridn]-potn*(negative[gridPoint]/(0.5*(1+delphiExp))+negative[gridn]/(0.5*(1+1/delphiExp))));
				}
				else {
					positiveCurrentProfile[grid[dim]-1]+=dix[2*dim+1][borderCount]*(positive[gridPoint]-positive[gridn]+potn*(positive[gridPoint]+positive[gridn]));
					negativeCurrentProfile[grid[dim]-1]+=dix[2*dim+1][borderCount]*(negative[gridPoint]-negative[gridn]-potn*(negative[gridPoint]+negative[gridn]));
				}
				borderCount++;
				nextBorderPoint = borderPoints[borderCount];
			}
		}
		count++;
	}

	for(i=0;i<gridSize[dim]-1;i++) {
		positiveCurrentProfile[i]*=-flcon;
		negativeCurrentProfile[i]*=flcon;
	}

	for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++) {
		potential[gridPoint]*=2;
		positive[gridPoint]*=fpoh;
		negative[gridPoint]*=fpoh;
	}
}
void NernstPlankSolver::nernstPlanckSolverCurrentProfileN(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile)
{
	int * GS;
	int GS_X,GS_Y,GS_Z;
	int GS_XY;
	int GS_XYZ;
	int gridPoint;
	float gridScale;
	float fpoh;
	float flcon;
	
	int i,j,k,ion;
	int kgrid,jgrid;
	int GrdPnt;
	int gridp,gridn;
	float pot;
	int grid[3];
	int gridInc;
	float *C[2];
	double *I[2];
	float * potential;
	
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
	I[0] = positiveCurrentProfile;
	I[1] = negativeCurrentProfile;
	potential = World->Potential;
	
	
	fpoh = 4*M_PI*gridScale;
	flcon = 1.6*gridScale*gridScale*1000;
	flcon/=fpoh;

	gridInc = GS_XY;
	
	for(i=0;i<World->NIonsTypes;i++)
		for(j=0;j<GS_Z-1;j++)
			I[i][j]=0.0;

	if(UTMPSingle == NULL)if(!(UTMPSingle = new float[GS_XYZ])){
		cerr << "ERROR 104: No memory available\n";
		exit(104);
	}
	
	float d,dc,sc,I1,I2;
	float *utmp=UTMPSingle;
	for(ion=0;ion<World->NIonsTypes;ion++)
	{	
		if(World->PMF!=NULL)
		{
			if(potential!=NULL)
			{
				if(PMFWeight>=0.0)
				{
					for(j=0;j<GS_XYZ;j++)
						utmp[j]=0.5*(World->IonsQ[ion]*potential[j]+World->PMF[ion][j]*PMFWeight);
				}
				else
				{
					for(j=0;j<GS_XYZ;j++)
						utmp[j]=0.5*(World->IonsQ[ion]*potential[j]+World->PMF[ion][j]);
				}
			}
			else
			{
				if(PMFWeight>=0.0)
				{
					for(j=0;j<GS_XYZ;j++)
						utmp[j]=0.5*(World->PMF[ion][j])*PMFWeight;
				}
				else
				{
					for(j=0;j<GS_XYZ;j++)
						utmp[j]=0.5*(World->PMF[ion][j]);
				}
			}
		}
		else
			for(j=0;j<GS_XYZ;j++)
				utmp[j]=0.5*(World->IonsQ[ion]*potential[j]);
		if(QmobMask==NULL)
		{
			for(k=0;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+GS_XY;
						if(NIndexing->GetDiffFloat(ion,gridPoint)>0.0)
						{
							pot = utmp[gridp]-utmp[gridPoint];
							d = NIndexing->GetDiffFloat(ion,gridPoint+GS_XY)>0.0?					0.5*(NIndexing->GetDiffFloat(ion,gridPoint)+NIndexing->GetDiffFloat(ion,gridPoint+GS_XY)):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							I[ion][k]+=d*(dc+pot*sc);
						}
					}
				}
			}
		}
		else
		{
			for(k=0;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+GS_XY;
						if(NIndexing->GetDiffFloat(ion,gridPoint)>0.0&&QmobMask[ion][gridPoint]>0)
						{
							pot = utmp[gridp]-utmp[gridPoint];
							d = NIndexing->GetDiffFloat(ion,gridPoint+GS_XY)>0.0?					0.5*(NIndexing->GetDiffFloat(ion,gridPoint)+NIndexing->GetDiffFloat(ion,gridPoint+GS_XY)):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							I[ion][k]+=d*(dc+pot*sc);
						}
					}
				}
			}
		}
	}
	for(i=0;i<GS[dim]-1;i++)
	{
		positiveCurrentProfile[i]*=-flcon;
		negativeCurrentProfile[i]*=flcon;
	}
}
void NernstPlankSolver::nernstPlanckSolverCurrentProfileAD(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile)
{
	int * GS;
	int GS_X,GS_Y,GS_Z;
	int GS_XY;
	int GS_XYZ;
	int gridPoint;
	float gridScale;
	float fpoh;
	float flcon;
	
	int i,j,k,ion;
	int kgrid,jgrid;
	int GrdPnt;
	int gridp,gridn;
	float pot;
	int grid[3];
	int gridInc;
	float *C[2];
	double *I[2];
	float * potential;
	float **D=World->D;
	
	GS = World->GridSize;
	GS_X = GS[0];
	GS_Y = GS[1];
	GS_Z = GS[2];
	GS_XY = GS[0]*GS[1];
	GS_XYZ = GS[2]*GS_XY;
	gridScale = World->GridScale;
	
	C[0] = World->C[0];
	C[1] = World->C[1];
	I[0] = positiveCurrentProfile;
	I[1] = negativeCurrentProfile;
	potential = World->Potential;


	fpoh = 4*M_PI*gridScale;
	flcon = 1.6*gridScale*gridScale*1000;
	flcon/=fpoh;

	gridInc = GS_XY;
	
	for(i=0;i<World->NIonsTypes;i++)
		for(j=0;j<GS_Z-1;j++)
			I[i][j]=0.0;

	if(UTMPSingle == NULL)if(!(UTMPSingle = new float[GS_XYZ])){
		cerr << "ERROR 104: No memory available\n";
		exit(104);
	}
	
	float d,dc,sc,I1,I2;
	float *utmp=UTMPSingle;
	for(ion=0;ion<World->NIonsTypes;ion++)
	{	
		if(World->PMF!=NULL)
		{
			if(potential!=NULL)
			{
				if(PMFWeight>=0.0)
				{
					for(j=0;j<GS_XYZ;j++)
						utmp[j]=0.5*(World->IonsQ[ion]*potential[j]+World->PMF[ion][j]*PMFWeight);
				}
				else
				{
					for(j=0;j<GS_XYZ;j++)
						utmp[j]=0.5*(World->IonsQ[ion]*potential[j]+World->PMF[ion][j]);
				}
			}
			else
			{
				if(PMFWeight>=0.0)
				{
					for(j=0;j<GS_XYZ;j++)
						utmp[j]=0.5*(World->PMF[ion][j])*PMFWeight;
				}
				else
				{
					for(j=0;j<GS_XYZ;j++)
						utmp[j]=0.5*(World->PMF[ion][j]);
				}
			}
		}
		else
			for(j=0;j<GS_XYZ;j++)
				utmp[j]=0.5*(World->IonsQ[ion]*potential[j]);
		
		if(QmobMask==NULL)
		{
			for(k=0;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+GS_XY;
						if(D[ion][gridPoint]>0.0f)
						{
							pot = utmp[gridp]-utmp[gridPoint];
							d = D[ion][gridPoint+GS_XY]>0.0f?					0.5*(D[ion][gridPoint]+D[ion][gridPoint+GS_XY]):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							I[ion][k]+=d*(dc+pot*sc);
						}
					}
				}
			}
		}
		else
		{
			for(k=0;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+GS_XY;
						if(D[ion][gridPoint]>0.0f&&QmobMask[ion][gridPoint]>0)
						{
							pot = utmp[gridp]-utmp[gridPoint];
							d = D[ion][gridPoint+GS_XY]>0.0f?					0.5*(D[ion][gridPoint]+D[ion][gridPoint+GS_XY]):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							I[ion][k]+=d*(dc+pot*sc);
						}
					}
				}
			}
		}
	}
	for(i=0;i<GS[dim]-1;i++)
	{
		positiveCurrentProfile[i]*=-flcon;
		negativeCurrentProfile[i]*=flcon;
	}
}
void NernstPlankSolver::nernstPlanckSolverCurrentProfileDflex(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile)
{
	int * GS;
	int GS_X,GS_Y,GS_Z;
	int GS_XY;
	int GS_XYZ;
	int gridPoint;
	float gridScale;
	float fpoh;
	float flcon;
	
	int i,j,k,ion;
	int kgrid,jgrid;
	int GrdPnt;
	int gridp,gridn;
	float pot;
	int grid[3];
	int gridInc;
	float *C[2];
	double *I[2];
	float * potential;
	float **D=World->D;
	
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
	I[0] = positiveCurrentProfile;
	I[1] = negativeCurrentProfile;
	potential = World->Potential;


	
	
	fpoh = 4*M_PI*gridScale;
	flcon = 1.6*gridScale*gridScale*1000;
	flcon/=fpoh;

	gridInc = GS_XY;
	
	for(i=0;i<World->NIonsTypes;i++)
		for(j=0;j<GS_Z-1;j++)
			I[i][j]=0.0;

	if(UTMPSingle == NULL)if(!(UTMPSingle = new float[GS_XYZ])){
		cerr << "ERROR 104: No memory available\n";
		exit(104);
	}
	
	float d,dc,sc,I1,I2;
	float *utmp=UTMPSingle;
	for(ion=0;ion<World->NIonsTypes;ion++)
	{	
		if(World->PMF!=NULL)
		{
			if(potential!=NULL)
			{
				if(PMFWeight>=0.0)
				{
					for(j=0;j<GS_XYZ;j++)
						utmp[j]=0.5*(World->IonsQ[ion]*potential[j]+World->PMF[ion][j]*PMFWeight);
				}
				else
				{
					for(j=0;j<GS_XYZ;j++)
						utmp[j]=0.5*(World->IonsQ[ion]*potential[j]+World->PMF[ion][j]);
				}
			}
			else
			{
				if(PMFWeight>=0.0)
				{
					for(j=0;j<GS_XYZ;j++)
						utmp[j]=0.5*(World->PMF[ion][j])*PMFWeight;
				}
				else
				{
					for(j=0;j<GS_XYZ;j++)
						utmp[j]=0.5*(World->PMF[ion][j]);
				}
			}
		}
		else
			for(j=0;j<GS_XYZ;j++)
				utmp[j]=0.5*(World->IonsQ[ion]*potential[j]);
		if(QmobMask==NULL)
		{
			for(k=0;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+GS_XY;
						if(NIndexing->GetDiffFloat(ion,gridPoint)>0.0)
						{
							pot = utmp[gridp]-utmp[gridPoint];
							d = NIndexing->GetDiffFloat(ion,gridPoint+GS_XY)>0.0?					0.5*(D[ion][gridPoint]+D[ion][gridPoint+GS_XY]):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							I[ion][k]+=d*(dc+pot*sc);
						}
					}
				}
			}
		}
		else
		{
			for(k=0;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+GS_XY;
						if(NIndexing->GetDiffFloat(ion,gridPoint)>0.0&&QmobMask[ion][gridPoint]>0)
						{
							pot = utmp[gridp]-utmp[gridPoint];
							d = NIndexing->GetDiffFloat(ion,gridPoint+GS_XY)>0.0?					0.5*(D[ion][gridPoint]+D[ion][gridPoint+GS_XY]):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							I[ion][k]+=d*(dc+pot*sc);
						}
					}
				}
			}
		}
	}
	for(i=0;i<GS[dim]-1;i++)
	{
		positiveCurrentProfile[i]*=-flcon;
		negativeCurrentProfile[i]*=flcon;
	}
}

//======================================DOUBLE
#ifdef PNPDOUBLE
NernstPlankSolverDouble::NernstPlankSolverDouble()
	: HaObject(), GenericSolver()
{
	InitZero();
}


NernstPlankSolverDouble::~NernstPlankSolverDouble()
{
	Clear();
}
int NernstPlankSolverDouble::InitZero()
{
	HaObject::SetName("NernstPlankSolverDouble");
	SolverStr.push_back("Auto");
	SolverStr.push_back("NodeIndexBased");
	SolverStr.push_back("ArrayDirect");
	SolverStr.push_back("PNPC");
	
	World=NULL;
	ConvergenceCheck=20;
	PMFWeight=-1.0;
	int i,j;
  //Initial values of variables:
	verbose=true;
	for(i=0;i<2;i++)
	{
		UTMP[i] = NULL;
		TMP[i]  = NULL;
	}
	QmobMask=NULL;
	solver=Auto;
	for(i=0;i<2;i++)
	{
		dixt[i] = NULL;
		for(j=0;j<6;j++)
		{
			dix[i][j] = NULL;
			dphi[i][j] = NULL;
		}
	}
	return EXIT_SUCCESS;
}
int NernstPlankSolverDouble::Clear()
{
	int i,j,IType;
	DeleteCArray(QmobMask);
	if(World!=NULL)
	{
		DeleteCVecArray(IndexNoSingular,World->NIonsTypes);
		DeleteCVecArray(IndexSingular,World->NIonsTypes);
		DeleteCArray(UTMP,World->NIonsTypes);
		DeleteCVecArray(TMP,World->NIonsTypes);
		DeleteCVecArray(dixt,World->NIonsTypes);
		for(IType=0;IType<2;IType++)
		{
			for(j=0;j<6;j++)if(dix[IType][j]!=0)delete [] dix[IType][j];
			for(j=0;j<6;j++)if(dphi[IType][j]!=0)delete [] dphi[IType][j];
		}
	}
}
int NernstPlankSolverDouble::SaveXML(TiXmlElement* Elt, HaContext* p_ctxt )
{
	return EXIT_SUCCESS;
}
int NernstPlankSolverDouble::LoadXML(const TiXmlElement* Elt, HaContext* p_ctxt )
{
	if(Elt==NULL)
	{
		fprintf(stderr,"ERROR: Can't find <%s>\n",HaObject::GetCStrName());
		return EXIT_FAILURE;
	}
	if(strcmp(HaObject::GetCStrName(),Elt->Value()))
	{
		fprintf(stderr,"ERROR: Wrong XML Element %s, expecting %s\n",Elt->Value(),HaObject::GetCStrName());
		return EXIT_FAILURE;
	}
	Clear();
	int i,gridPoint;
	//Read Primary Parameters
	Elt->GetIntAttribute("MaxIterations",&MaxIterations);
	Elt->GetDoubleAttribute("Convergence",&Convergence);
	Elt->GetDoubleAttribute("Relaxation",&Relaxation);  
	
	if(Elt->GetIntAttribute("ConvergenceCheck",&ConvergenceCheck)!=EXIT_SUCCESS)ConvergenceCheck=20;
	if(Elt->GetStdStrIndex("Solver",&solver,SolverStr)!=EXIT_SUCCESS)solver=0;
	if(Elt->GetStdStrAttribute("QmobMaskFile",&QmobMaskFile)!=EXIT_SUCCESS)QmobMaskFile="";
	if(Elt->GetBoolAttribute("Verbose",&verbose)!=EXIT_SUCCESS)verbose=true;
	
	if(Elt->GetDoubleAttribute("PMFWeight",&PMFWeight)!=EXIT_SUCCESS)PMFWeight=-1.0;
	else if(PMFWeight>1.0)PMFWeight=-1.0;
	
	ShowParameters();
	return EXIT_SUCCESS;
}
int NernstPlankSolverDouble::SetContWorld(ContWorld* _world)
{
	int i;
	World=_world;
	if(World->C==NULL)
	{
		World->SetInitConcentrationFromNIndexing();
	}
	else if(World->C[0]==NULL)
	{
		World->SetInitConcentrationFromNIndexing();
	}
	return EXIT_SUCCESS;
}
int NernstPlankSolverDouble::ShowParameters()
{
	fprintf(stdout,"\nParameters of Nernst-Plank solver\n");
	fprintf(stdout,"    MaxIterations:.................... %d\n", MaxIterations);
	fprintf(stdout,"    Convergence:...................... %g\n", Convergence);
	fprintf(stdout,"    Relaxation:....................... %g\n", Relaxation);
	return EXIT_SUCCESS;
}
int NernstPlankSolverDouble::InitSolver()
{
//	 int i;
//	 for(i=0;i<World->NIonsTypes;i++){
//		 if(UTMP[i] == NULL)if(!(UTMP[i] = new double[World->GridSizeXYZ])){
//			 cerr << "ERROR 104: No memory available\n";
//			 exit(104);
//		 }
//		 if(TMP[i] == NULL)if(!(TMP[i] = new double[World->GridSizeXYZ])){
//			 cerr << "ERROR 104: No memory available\n";
//			 exit(104);
//		 }
//	 }
	
	if(solver==Auto)
	{
		if(World->NIndexing==NULL)return InitSolverAD();
		else if(World->D!=NULL)return InitSolverDflex();
		else return InitSolverD();
	}
	else
	{
		if(solver==NodeIndexBased)
			return InitSolverD();
		else if(solver==ArrayDirect)
			return InitSolverAD();
		else return InitSolverD();
	}
	
}
int NernstPlankSolverDouble::InitSolverD()
{
	int i,j,k, i1, IType;
	int ix, iy, iz;
	int GrdPnt;
	int jgrid,kgrid;
	//double * diffusion;
	//double **D;
	//int NIonsTypes;
	int GS_X = World->GridSize[0];
	int GS_Y = World->GridSize[1];
	int GS_Z = World->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	//int *tmp;
	int itmp1,itmp2;
	//int BlackOrWhite;
	
	NodeIndexing* NIndexing=World->NIndexing;
	NodeIndex* NIndex=NIndexing->NIndex;
	
	NodeIndex DiffBoarderMask=NodeIndexing::DiffIon0BoarderMask;
	NodeIndex BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;
	
//	 UTMP=0;
//	 UTMP[1]=0;
//	 TMP=0;
//	 TMP[1]=0;
	for(i=0;i<2;i++)for(j=0;j<6;j++)dphi[i][j]=NULL;
	MaxChange=0.0;
	
	if(!World->CDouble)
	{
		PNP_EXIT_FAIL_NULL(World->C,"ERROR 311: C Map Not init\n");
		PNP_EXIT_FAIL_NULL(World->C[0],"ERROR 311: C[0] Map Not init\n");
		PNP_EXIT_FAIL_NULL(World->C[1],"ERROR 311: C[1] Map Not init\n");
		
		World->CDouble=new double*[World->NIonsTypes];
		for(j=0;j<World->NIonsTypes;j++)
		{
			World->CDouble[j]=new double[GS_XYZ];
			for(i=0;i<GS_XYZ;i++)
				World->CDouble[j][i]=World->C[j][i];
		}
		
	}
	PNP_EXIT_FAIL_NULL(World->NIndexing,"ERROR 311: World->NIndexing Map Not init\n");
	
	if(!World->PotentialDouble){
		cerr << "WARNING 311: Potential Map Not init, will init and set to zero values\n";
		World->PotentialDouble=new double[GS_XYZ];
		for(i=0;i<GS_XYZ;i++)World->PotentialDouble[i]=0.0;
	}
	

	/*if(!World->PotentialTMP){
	if(!(World->PotentialTMP = new double[World->GridSizeXYZ])){
	cerr << "ERROR 104: No memory available\n";
	exit(104);
}
}*/
//cout<<"NP2sleep 10 sec\n";
	//sleep(10);
	
	//D=World->D;
	int DiffZeroInd=NIndexing->GetDiffZeroInd();
	
	for(IType=0;IType<World->NIonsTypes;IType++)
	{
		for(k=1;k<GS_Z-1;k++)
		{
			kgrid = k*GS_XY;
			for(j=1;j<GS_Y-1;j++)
			{
				jgrid = kgrid+j*GS_X;
				for(i=1;i<GS_X-1;i++)
				{
					GrdPnt = jgrid+i;
					if(NIndexing->GetDiffDouble(IType,GrdPnt)!=0.0 && NIndexing->GetDiffDouble(IType,GrdPnt+1)==0.0 && NIndexing->GetDiffDouble(IType,GrdPnt-1)==0.0 && NIndexing->GetDiffDouble(IType,GrdPnt+GS_X)==0.0 && NIndexing->GetDiffDouble(IType,GrdPnt-GS_X)==0.0 && NIndexing->GetDiffDouble(IType,GrdPnt+GS_XY)==0.0 &&NIndexing->GetDiffDouble(IType,GrdPnt-GS_XY)==0.0)
					{
						NIndexing->SetIonField(IType,GrdPnt,DiffZeroInd);
						DbgPrint2("Bed Diffusional point at [%d] have removed it.",GrdPnt);
					}
				}
			}
		}
	}
	for(IType=0;IType<World->NIonsTypes;IType++)
		for(i=0;i<GS_XYZ;i++)
			if(NIndexing->GetDiffDouble(IType,i)==0.0)World->C[IType][i]=0.0f;
	
	
	for(IType=0;IType<World->NIonsTypes;IType++)
	{
		DiffBoarderMask=NIndexing->DiffBoarderMask[IType];
		//Calculate number of nodes with singularity and not and exclude nodes with D=0
		
		NoSingularNum[IType][0] = 0;
		NoSingularNum[IType][1] = 0;
		NoSingularNum[IType][2] = 0;
		SingularNum[IType][0] = 0;
		SingularNum[IType][1] = 0;
		SingularNum[IType][2] = 0;
		if(QmobMask==NULL)
		{
			for(k=1;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						GrdPnt = jgrid+i;
						if(NIndexing->GetDiffDouble(IType,GrdPnt)>0.0)
						{
							if(NIndex[GrdPnt]&DiffBoarderMask)
							{
								SingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
								SingularNum[IType][2]++;
							}
							else{
								NoSingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
								NoSingularNum[IType][2]++;
							}
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
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						GrdPnt = jgrid+i;
						if(QmobMask[GrdPnt]>0)
							if(NIndexing->GetDiffDouble(IType,GrdPnt)>0.0)
						{
							if(NIndex[GrdPnt]&DiffBoarderMask)
							{
								SingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
								SingularNum[IType][2]++;
							}
							else{
								NoSingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
								NoSingularNum[IType][2]++;
							}
						}
					}
				}
			}
		}
			
		fprintf(stdout,"Nernst-Plank Solver:\n");
		fprintf(stdout,"		IonsQ:........................ %g\n", World->IonsQ[IType]);
		fprintf(stdout,"		SingularNum:...................... %d\n", SingularNum[IType][2]);
		fprintf(stdout,"		NoSingularNum:.................... %d\n", NoSingularNum[IType][2]);
		fprintf(stdout,"		IndexSingularParts:............... [%d,%d,%d]\n",
						SingularNum[IType][0],SingularNum[IType][1],SingularNum[IType][2]);
		fprintf(stdout,"		IndexNoSingularParts:............. [%d,%d,%d]\n",NoSingularNum[IType][0],NoSingularNum[IType][1],NoSingularNum[IType][2]);
		
		if(!(IndexSingular[IType] = new int[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		if(!(IndexNoSingular[IType] = new int[NoSingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		//cout<<"NP3sleep 10 sec\n";
		//sleep(10);
		for(i=0;i<6;i++)if(!(dix[IType][i] = new double[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		if(!(dixt[IType] = new double[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		int iSingular=0,iNoSingular=0;
		int iSingular2=SingularNum[IType][1],iNoSingular2=NoSingularNum[IType][1];
		if(QmobMask==NULL)
		{
			for(k=1;k<GS_Z-1;k++) 
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++) 
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++) 
					{
						GrdPnt = jgrid+i;
						if(NIndexing->GetDiffDouble(IType,GrdPnt)>0.0)
						{
							if(NIndex[GrdPnt]&DiffBoarderMask)
							{
								if(NIndex[GrdPnt]&BlackAndWhiteMask)
								{
									IndexSingular[IType][iSingular]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iSingular++;
								}
								else
								{
									IndexSingular[IType][iSingular2]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular2, GrdPnt);
									iSingular2++;
								}
							}
							else
							{
								if(NIndex[GrdPnt]&BlackAndWhiteMask)
								{
									IndexNoSingular[IType][iNoSingular]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular++;
								}
								else
								{
									IndexNoSingular[IType][iNoSingular2]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular2++;
								}
							}
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
				for(j=1;j<GS_Y-1;j++) 
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++) 
					{
						GrdPnt = jgrid+i;
						if(NIndexing->GetDiffDouble(IType,GrdPnt)>0.0)
						{
							
							if(QmobMask[GrdPnt]>0)
								if(NIndex[GrdPnt]&DiffBoarderMask)
							{
								if(NIndex[GrdPnt]&BlackAndWhiteMask)
								{
									IndexSingular[IType][iSingular]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iSingular++;
								}
								else
								{
									IndexSingular[IType][iSingular2]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular2, GrdPnt);
									iSingular2++;
								}
							}
							else
							{
								if(NIndex[GrdPnt]&BlackAndWhiteMask)
								{
									IndexNoSingular[IType][iNoSingular]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular++;
								}
								else
								{
									IndexNoSingular[IType][iNoSingular2]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular2++;
								}
							}
						}
					}
				}
			}
		}
		
		for(i=0;i<SingularNum[IType][2];i++)
		{
			GrdPnt=IndexSingular[IType][i];
			dix[IType][0][i] = NIndexing->GetDiffDouble(IType,GrdPnt+1)>0.0?
					0.5*(NIndexing->GetDiffDouble(IType,GrdPnt)+NIndexing->GetDiffDouble(IType,GrdPnt+1)):0.0;
			dix[IType][1][i] = NIndexing->GetDiffDouble(IType,GrdPnt-1)>0.0?
					0.5*(NIndexing->GetDiffDouble(IType,GrdPnt)+NIndexing->GetDiffDouble(IType,GrdPnt-1)):0.0;
			dix[IType][2][i] = NIndexing->GetDiffDouble(IType,GrdPnt+GS_X)>0.0?
					0.5*(NIndexing->GetDiffDouble(IType,GrdPnt)+NIndexing->GetDiffDouble(IType,GrdPnt+GS_X)):0.0;
			dix[IType][3][i] = NIndexing->GetDiffDouble(IType,GrdPnt-GS_X)>0.0?
					0.5*(NIndexing->GetDiffDouble(IType,GrdPnt)+NIndexing->GetDiffDouble(IType,GrdPnt-GS_X)):0.0;
			dix[IType][4][i] = NIndexing->GetDiffDouble(IType,GrdPnt+GS_XY)>0.0?
					0.5*(NIndexing->GetDiffDouble(IType,GrdPnt)+NIndexing->GetDiffDouble(IType,GrdPnt+GS_XY)):0.0;
			dix[IType][5][i] = NIndexing->GetDiffDouble(IType,GrdPnt-GS_XY)>0.0?
					0.5*(NIndexing->GetDiffDouble(IType,GrdPnt)+NIndexing->GetDiffDouble(IType,GrdPnt-GS_XY)):0.0;
			dixt[IType][i] = dix[IType][0][i]+dix[IType][1][i]+dix[IType][2][i]+dix[IType][3][i]+dix[IType][4][i]+dix[IType][5][i];
		}
	}
	
	if(World->PMF!=NULL)fprintf(stdout,"Calculate NP with PMF\n");
	/*for(i=0;i<World->GridSizeXYZ;i++)tmp[i]=0;
	for(i=0;i<SingularNum;i++)tmp[IndexSingular[i]]=1;
	for(i=0;i<NoSingularNum;i++)tmp[IndexNoSingular[i]]=2;
	World->writeMapInt("tmp.gz",tmp,World->GridSizeXYZ,1);
	delete [] tmp;*/
	
	//SAVE memory for children! Delete excess of D[0]
	//if(World->D[0]!=0)delete [] World->D[0];
	//cout<<"sleep 10 sec\n";
	//sleep(10);
	return 1;
}
int NernstPlankSolverDouble::InitSolverDflex()
{
	int i,j,k, i1, IType;
	int ix, iy, iz;
	int GrdPnt;
	int jgrid,kgrid;
	//double * diffusion;
	//double **D;
	//int NIonsTypes;
	int GS_X = World->GridSize[0];
	int GS_Y = World->GridSize[1];
	int GS_Z = World->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	//int *tmp;
	int itmp1,itmp2;
	float **D=World->D;
	//int BlackOrWhite;
	
	NodeIndexing* NIndexing=World->NIndexing;
	NodeIndex* NIndex=NIndexing->NIndex;
	
	NodeIndex DiffBoarderMask=NodeIndexing::DiffIon0BoarderMask;
	NodeIndex BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;

//	 UTMP=0;
//	 UTMP[1]=0;
//	 TMP=0;
//	 TMP[1]=0;
	for(i=0;i<2;i++)for(j=0;j<6;j++)dphi[i][j]=NULL;
	MaxChange=0.0;
	
	
	if(!World->CDouble)
	{
		PNP_EXIT_FAIL_NULL(World->C,"ERROR 311: C Map Not init\n");
		PNP_EXIT_FAIL_NULL(World->C[0],"ERROR 311: C[0] Map Not init\n");
		PNP_EXIT_FAIL_NULL(World->C[1],"ERROR 311: C[1] Map Not init\n");
		
		World->CDouble=new double*[World->NIonsTypes];
		for(j=0;j<World->NIonsTypes;j++)
		{
			World->CDouble[j]=new double[GS_XYZ];
			for(i=0;i<GS_XYZ;i++)
				World->CDouble[j][i]=World->C[j][i];
		}
		
	}
	PNP_EXIT_FAIL_NULL(World->NIndexing,"ERROR 311: World->NIndexing Map Not init\n");
	
	if(!World->PotentialDouble){
		cerr << "WARNING 311: Potential Map Not init, will init and set to zero values\n";
		World->PotentialDouble=new double[GS_XYZ];
		for(i=0;i<GS_XYZ;i++)World->PotentialDouble[i]=0.0;
	}
	

	/*if(!World->PotentialTMP){
	if(!(World->PotentialTMP = new double[World->GridSizeXYZ])){
	cerr << "ERROR 104: No memory available\n";
	exit(104);
}
}*/
//cout<<"NP2sleep 10 sec\n";
	//sleep(10);
	
	//D=World->D;
	for(IType=0;IType<World->NIonsTypes;IType++)
	{
		for(k=1;k<GS_Z-1;k++)
		{
			kgrid = k*GS_XY;
			for(j=1;j<GS_Y-1;j++)
			{
				jgrid = kgrid+j*GS_X;
				for(i=1;i<GS_X-1;i++)
				{
					GrdPnt = jgrid+i;
					if(D[IType][GrdPnt]!=0.0 && D[IType][GrdPnt+1]==0.0 && D[IType][GrdPnt-1]==0.0 && D[IType][GrdPnt+GS_X]==0.0 && D[IType][GrdPnt-GS_X]==0.0 && D[IType][GrdPnt+GS_XY]==0.0 &&D[IType][GrdPnt-GS_XY]==0.0)
					{
						D[IType][GrdPnt]=0.0;
						DbgPrint2("Bed Diffusional point at [%d] have removed it.",GrdPnt);
					}
				}
			}
		}
	}

	for(IType=0;IType<World->NIonsTypes;IType++)
		for(i=0;i<GS_XYZ;i++)
			if(World->D[IType][i]==0.0)World->C[IType][i]=0.0f;
	
	for(IType=0;IType<World->NIonsTypes;IType++)
	{
		//DiffBoarderMask=NIndexing->DiffBoarderMask[IType];
		//Calculate number of nodes with singularity and not and exclude nodes with D=0
		
		NoSingularNum[IType][0] = 0;
		NoSingularNum[IType][1] = 0;
		NoSingularNum[IType][2] = 0;
		SingularNum[IType][0] = 0;
		SingularNum[IType][1] = 0;
		SingularNum[IType][2] = 0;
		if(QmobMask==NULL)
		{
			for(k=1;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						GrdPnt = jgrid+i;
						if(NIndexing->GetDiffDouble(IType,GrdPnt)>0.0){
							if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
							{
								SingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
								SingularNum[IType][2]++;
							}
							else{
								NoSingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
								NoSingularNum[IType][2]++;
							}
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
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						GrdPnt = jgrid+i;
						if(QmobMask[GrdPnt]>0)
							if(NIndexing->GetDiffDouble(IType,GrdPnt)>0.0)
						{
							if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
							{
								SingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
								SingularNum[IType][2]++;
							}
							else{
								NoSingularNum[IType][1]+=NIndex[GrdPnt]&BlackAndWhiteMask;
								NoSingularNum[IType][2]++;
							}
						}
					}
				}
			}
		}
			
		fprintf(stdout,"Nernst-Plank Solver:\n");
		fprintf(stdout,"		IonsQ:........................ %g\n", World->IonsQ[IType]);
		fprintf(stdout,"		SingularNum:...................... %d\n", SingularNum[IType][2]);
		fprintf(stdout,"		NoSingularNum:.................... %d\n", NoSingularNum[IType][2]);
		fprintf(stdout,"		IndexSingularParts:............... [%d,%d,%d]\n",
						SingularNum[IType][0],SingularNum[IType][1],SingularNum[IType][2]);
		fprintf(stdout,"		IndexNoSingularParts:............. [%d,%d,%d]\n",NoSingularNum[IType][0],NoSingularNum[IType][1],NoSingularNum[IType][2]);
		
		if(!(IndexSingular[IType] = new int[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		if(!(IndexNoSingular[IType] = new int[NoSingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		//cout<<"NP3sleep 10 sec\n";
		//sleep(10);
		for(i=0;i<6;i++)if(!(dix[IType][i] = new double[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		if(!(dixt[IType] = new double[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		int iSingular=0,iNoSingular=0;
		int iSingular2=SingularNum[IType][1],iNoSingular2=NoSingularNum[IType][1];
		if(QmobMask==NULL)
		{
			for(k=1;k<GS_Z-1;k++) 
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++) 
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++) 
					{
						GrdPnt = jgrid+i;
						if(NIndexing->GetDiffDouble(IType,GrdPnt)>0.0)
						{
							if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
							{
								if(NIndex[GrdPnt]&BlackAndWhiteMask)
								{
									IndexSingular[IType][iSingular]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iSingular++;
								}
								else
								{
									IndexSingular[IType][iSingular2]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular2, GrdPnt);
									iSingular2++;
								}
							}
							else
							{
								if(NIndex[GrdPnt]&BlackAndWhiteMask)
								{
									IndexNoSingular[IType][iNoSingular]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular++;
								}
								else
								{
									IndexNoSingular[IType][iNoSingular2]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular2++;
								}
							}
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
				for(j=1;j<GS_Y-1;j++) 
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++) 
					{
						GrdPnt = jgrid+i;
						if(NIndexing->GetDiffDouble(IType,GrdPnt)>0.0)
						{
							
							if(QmobMask[GrdPnt]>0)
								if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
							{
								if(NIndex[GrdPnt]&BlackAndWhiteMask)
								{
									IndexSingular[IType][iSingular]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iSingular++;
								}
								else
								{
									IndexSingular[IType][iSingular2]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular2, GrdPnt);
									iSingular2++;
								}
							}
							else
							{
								if(NIndex[GrdPnt]&BlackAndWhiteMask)
								{
									IndexNoSingular[IType][iNoSingular]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular++;
								}
								else
								{
									IndexNoSingular[IType][iNoSingular2]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular2++;
								}
							}
						}
					}
				}
			}
		}
		
		for(i=0;i<SingularNum[IType][2];i++)
		{
			GrdPnt=IndexSingular[IType][i];
			dix[IType][0][i] = NIndexing->GetDiffDouble(IType,GrdPnt+1)>0.0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt+1]):0.0;
			dix[IType][1][i] = NIndexing->GetDiffDouble(IType,GrdPnt-1)>0.0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt-1]):0.0;
			dix[IType][2][i] = NIndexing->GetDiffDouble(IType,GrdPnt+GS_X)>0.0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt+GS_X]):0.0;
			dix[IType][3][i] = NIndexing->GetDiffDouble(IType,GrdPnt-GS_X)>0.0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt-GS_X]):0.0;
			dix[IType][4][i] = NIndexing->GetDiffDouble(IType,GrdPnt+GS_XY)>0.0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt+GS_XY]):0.0;
			dix[IType][5][i] = NIndexing->GetDiffDouble(IType,GrdPnt-GS_XY)>0.0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt-GS_XY]):0.0;
			dixt[IType][i] = dix[IType][0][i]+dix[IType][1][i]+dix[IType][2][i]+dix[IType][3][i]+dix[IType][4][i]+dix[IType][5][i];
		}
	}
	
	if(World->PMF!=NULL)fprintf(stdout,"Calculate NP with PMF\n");
	/*for(i=0;i<World->GridSizeXYZ;i++)tmp[i]=0;
	for(i=0;i<SingularNum;i++)tmp[IndexSingular[i]]=1;
	for(i=0;i<NoSingularNum;i++)tmp[IndexNoSingular[i]]=2;
	World->writeMapInt("tmp.gz",tmp,World->GridSizeXYZ,1);
	delete [] tmp;*/
	
	//SAVE memory for children! Delete excess of D[0]
	//if(World->D[0]!=0)delete [] World->D[0];
	//cout<<"sleep 10 sec\n";
	//sleep(10);
	return 1;
}
/////////////////////////////////////////////////////////////////////////////
int NernstPlankSolverDouble::InitSolverAD()
{
	int i,j,k, i1, IType;
	int ix, iy, iz;
	int GrdPnt;
	int jgrid,kgrid;
	//double * diffusion;
	//double **D;
	//int NIonsTypes;
	int GS_X = World->GridSize[0];
	int GS_Y = World->GridSize[1];
	int GS_Z = World->GridSize[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	//int *tmp;
	int itmp1,itmp2;
	
	int BlackOrWhite;
	
	for(i=0;i<2;i++)for(j=0;j<6;j++)dphi[i][j]=NULL;
	MaxChange=0.0;
	
	if(World->D==NULL)
	{
		if(World->NIndexing!=NULL)
		{
			pnpPrint0("WARNING: Will convert Diffusion from NIndexing\n");
			World->D=new float*[World->NIonsTypes];
			if(World->NIonsTypes>=1)
				World->D[0]=World->NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion0);
			if(World->NIonsTypes>=2)
				World->D[1]=World->NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion1);
			if(World->NIonsTypes>=3)
				World->D[2]=World->NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion2);
			if(World->NIonsTypes>=4)
				World->D[3]=World->NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion3);
		}
		else
		{
			fprintf(stderr, "ERR:		Diffusion Map NOT yet initialize\n");
		}
	}
	float **D=World->D;
	double **C=World->CDouble;
	if(!World->CDouble[0]){
		cerr << "ERROR 311: C[0] Map Not init\n";
		exit(104);
	}
	if(!World->CDouble[1]){
		cerr << "ERROR 311: C[1] Map Not init\n";
		exit(104);
	}
	if(!World->PotentialDouble){
		cerr << "WARNING 311: Potential Map Not init, will init and set to zero values\n";
		World->PotentialDouble=new double[GS_XYZ];
		for(i=0;i<GS_XYZ;i++)World->PotentialDouble[i]=0.0;
	}
	
	for(IType=0;IType<World->NIonsTypes;IType++)
	{
		for(k=1;k<GS_Z-1;k++)
		{
			kgrid = k*GS_XY;
			for(j=1;j<GS_Y-1;j++)
			{
				jgrid = kgrid+j*GS_X;
				for(i=1;i<GS_X-1;i++)
				{
					GrdPnt = jgrid+i;
					if(D[IType][GrdPnt]!=0.0 && D[IType][GrdPnt+1]==0.0 && D[IType][GrdPnt-1]==0.0 && D[IType][GrdPnt+GS_X]==0.0 && D[IType][GrdPnt-GS_X]==0.0 && D[IType][GrdPnt+GS_XY]==0.0 &&D[IType][GrdPnt-GS_XY]==0.0)
					{
						D[IType][GrdPnt]=0.0;
						C[IType][GrdPnt]=0.0;
						DbgPrint2("Bed Diffusional point at [%d] have removed it.",GrdPnt);
					}
				}
			}
		}
	}
	for(IType=0;IType<World->NIonsTypes;IType++)
		for(i=0;i<GS_XYZ;i++)
			if(D[IType][GrdPnt]==0.0)C[IType][i]=0.0f;
	
	
	for(IType=0;IType<World->NIonsTypes;IType++)
	{
		//DiffBoarderMask=NIndexing->DiffBoarderMask[IType];
		//Calculate number of nodes with singularity and not and exclude nodes with D=0
		
		NoSingularNum[IType][0] = 0;
		NoSingularNum[IType][1] = 0;
		NoSingularNum[IType][2] = 0;
		SingularNum[IType][0] = 0;
		SingularNum[IType][1] = 0;
		SingularNum[IType][2] = 0;
		if(QmobMask==NULL)
		{
			for(k=1;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						GrdPnt = jgrid+i;
						
						if(D[IType][GrdPnt]>0.0){
							BlackOrWhite=(k+j+i+World->startBlackAndWhite)%2;
							if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
							{
								SingularNum[IType][1]+=BlackOrWhite;
								SingularNum[IType][2]++;
							}
							else{
								NoSingularNum[IType][1]+=BlackOrWhite;
								NoSingularNum[IType][2]++;
							}
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
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						GrdPnt = jgrid+i;
						if(QmobMask[GrdPnt]>0)
							if(D[IType][GrdPnt]>0.0)
						{
							BlackOrWhite=(k+j+i+World->startBlackAndWhite)%2;
							if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
							{
								SingularNum[IType][1]+=BlackOrWhite;
								SingularNum[IType][2]++;
							}
							else{
								NoSingularNum[IType][1]+=BlackOrWhite;
								NoSingularNum[IType][2]++;
							}
						}
					}
				}
			}
		}
			
		fprintf(stdout,"Nernst-Plank Solver:\n");
		fprintf(stdout,"		IonsQ:........................ %g\n", World->IonsQ[IType]);
		fprintf(stdout,"		SingularNum:...................... %d\n", SingularNum[IType][2]);
		fprintf(stdout,"		NoSingularNum:.................... %d\n", NoSingularNum[IType][2]);
		fprintf(stdout,"		IndexSingularParts:............... [%d,%d,%d]\n",
						SingularNum[IType][0],SingularNum[IType][1],SingularNum[IType][2]);
		fprintf(stdout,"		IndexNoSingularParts:............. [%d,%d,%d]\n",NoSingularNum[IType][0],NoSingularNum[IType][1],NoSingularNum[IType][2]);
		
		if(!(IndexSingular[IType] = new int[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		if(!(IndexNoSingular[IType] = new int[NoSingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		//cout<<"NP3sleep 10 sec\n";
		//sleep(10);
		for(i=0;i<6;i++)if(!(dix[IType][i] = new double[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		if(!(dixt[IType] = new double[SingularNum[IType][2]])){
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
		int iSingular=0,iNoSingular=0;
		int iSingular2=SingularNum[IType][1],iNoSingular2=NoSingularNum[IType][1];
		if(QmobMask==NULL)
		{
			for(k=1;k<GS_Z-1;k++) 
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++) 
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++) 
					{
						GrdPnt = jgrid+i;
						if(D[IType][GrdPnt]>0.0)
						{
							BlackOrWhite=(k+j+i+World->startBlackAndWhite)%2;
							if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
												D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
							{
								if(BlackOrWhite)
								{
									IndexSingular[IType][iSingular]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iSingular++;
								}
								else
								{
									IndexSingular[IType][iSingular2]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular2, GrdPnt);
									iSingular2++;
								}
							}
							else
							{
								if(BlackOrWhite)
								{
									IndexNoSingular[IType][iNoSingular]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular++;
								}
								else
								{
									IndexNoSingular[IType][iNoSingular2]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular2++;
								}
							}
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
				for(j=1;j<GS_Y-1;j++) 
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++) 
					{
						GrdPnt = jgrid+i;
						if(D[IType][GrdPnt]>0.0)
						{
							BlackOrWhite=(k+j+i+World->startBlackAndWhite)%2;
							if(QmobMask[GrdPnt]>0)
								if(D[IType][GrdPnt]!=D[IType][GrdPnt+1] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt-1] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt+GS_X] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt-GS_X] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt+GS_XY] ||
													 D[IType][GrdPnt]!=D[IType][GrdPnt-GS_XY])
							{
								if(BlackOrWhite)
								{
									IndexSingular[IType][iSingular]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iSingular++;
								}
								else
								{
									IndexSingular[IType][iSingular2]=GrdPnt;
									//fprintf(stdout,"IndexSingular[%d][%d]=%d\n",IType,iSingular2, GrdPnt);
									iSingular2++;
								}
							}
							else
							{
								if(BlackOrWhite)
								{
									IndexNoSingular[IType][iNoSingular]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular++;
								}
								else
								{
									IndexNoSingular[IType][iNoSingular2]=GrdPnt;
									//fprintf(stdout,"IndexNoSingular[%d][%d]=%d\n",IType,iSingular, GrdPnt);
									iNoSingular2++;
								}
							}
						}
					}
				}
			}
		}
		
		for(i=0;i<SingularNum[IType][2];i++)
		{
			GrdPnt=IndexSingular[IType][i];
			dix[IType][0][i] = D[IType][GrdPnt+1]>0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt+1]):0;
			dix[IType][1][i] = D[IType][GrdPnt-1]>0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt-1]):0;
			dix[IType][2][i] = D[IType][GrdPnt+GS_X]>0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt+GS_X]):0;
			dix[IType][3][i] = D[IType][GrdPnt-GS_X]>0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt-GS_X]):0;
			dix[IType][4][i] = D[IType][GrdPnt+GS_XY]>0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt+GS_XY]):0;
			dix[IType][5][i] = D[IType][GrdPnt-GS_XY]>0?
					0.5*(D[IType][GrdPnt]+D[IType][GrdPnt-GS_XY]):0;
			dixt[IType][i] = dix[IType][0][i]+dix[IType][1][i]+dix[IType][2][i]+dix[IType][3][i]+dix[IType][4][i]+dix[IType][5][i];
		}
	}
	
	if(World->PMF!=NULL)fprintf(stdout,"Calculate NP with PMF\n");
	return 1;
}
int NernstPlankSolverDouble::Solve()
{
	if(solver==Auto)
	{
		if(World->D!=NULL)return NernstPlanckSolverD();
		else return NernstPlanckSolverD();
	}
	else
	{
		if(solver==NodeIndexBased)
			return NernstPlanckSolverD();
		else if(solver==ArrayDirect)
			return NernstPlanckSolverD();
		else return NernstPlanckSolverD();
	}
}
int NernstPlankSolverDouble::CheckSystem()
{
	int i,j,k;
	int IonType;
	int gridPoint;
	
	VectorIntField3D *iV=new VectorIntField3D(World->GridSize,World->GridScale,World->NIonsTypes);
	iV->FillValue(0);
	for(IonType=0;IonType<World->NIonsTypes;IonType++)
	{
		for(i=0;i<NoSingularNum[IonType][2];i++)
		{
			gridPoint=IndexNoSingular[IonType][i];
			iV->V[IonType][gridPoint]=1;
		}
		for(i=0;i<SingularNum[IonType][2];i++)
		{
			gridPoint=IndexSingular[IonType][i];
			iV->V[IonType][gridPoint]=2;
		}
	}
	
	int kgrid,jgrid,GrdPnt;
	int GS_X=World->GS_X;
	int GS_Y=World->GS_Y;
	int GS_Z=World->GS_Z;
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_X*GS_Y*GS_Z;
	
	for(IonType=0;IonType<World->NIonsTypes;IonType++)
	{
		for(k=1;k<GS_Z-1;k++)
		{
			for(j=1;j<GS_Y-1;j++)
			{
				for(i=1;i<GS_X-1;i++)
				{
					GrdPnt = i+j*GS_X+k*GS_XY;
					
					if(World->C[IonType][GrdPnt]==0.0&&iV->V[IonType][GrdPnt]!=0)
					{
						pnpPrint0("World->C[%d][%d]==0.0&&iV->V[%d][%d]!=0\n",IonType,GrdPnt,IonType,GrdPnt);
					}
					if(World->NIndexing->GetDiffDouble(IonType,GrdPnt)==0.0&&iV->V[IonType][GrdPnt]!=0)
					{
						pnpPrint0("NIndexing->GetDiffDouble(%d,%d)==0.0&&iV->V[%d][%d]!=0\n",IonType,GrdPnt,IonType,GrdPnt);
					}
				}
			}
		}
	}
	if(World->PMF!=NULL)
	{
		
		
		pnpPrint0("PMF present, will check values:\n");
		for(IonType=0;IonType<World->NIonsTypes;IonType++)
		{
			double dPMFx,dPMFy,dPMFz;
			int countPMFx=0,countPMFy=0,countPMFz=0;
			for(k=1;k<GS_Z;k++)
			{
				for(j=1;j<GS_Y;j++)
				{
					for(i=1;i<GS_X;i++)
					{
						GrdPnt = i+j*GS_X+k*GS_XY;
						dPMFx=World->PMF[IonType][GrdPnt]-World->PMF[IonType][GrdPnt-1];
						dPMFy=World->PMF[IonType][GrdPnt]-World->PMF[IonType][GrdPnt-GS_X];
						dPMFz=World->PMF[IonType][GrdPnt]-World->PMF[IonType][GrdPnt-GS_XY];
						
						if((dPMFx*dPMFx)>4.0)countPMFx++;
						if((dPMFy*dPMFy)>4.0)countPMFy++;
						if((dPMFz*dPMFz)>4.0)countPMFz++;
					}
				}
			}
			pnpPrint0("\tIon %d: countPMFx=%d,countPMFy=%d,countPMFz=%d Throught All Space\n"
					,IonType,countPMFx,countPMFy,countPMFz);
		}
		for(IonType=0;IonType<World->NIonsTypes;IonType++)
		{
			double dPMFx,dPMFy,dPMFz;
			int countPMFx=0,countPMFy=0,countPMFz=0;
			for(k=1;k<GS_Z;k++)
			{
				for(j=1;j<GS_Y;j++)
				{
					for(i=1;i<GS_X;i++)
					{
						GrdPnt = i+j*GS_X+k*GS_XY;
						if(World->NIndexing->GetDiffDouble(IonType,GrdPnt)!=0.0)
						{
							dPMFx=World->PMF[IonType][GrdPnt]-World->PMF[IonType][GrdPnt-1];
							dPMFy=World->PMF[IonType][GrdPnt]-World->PMF[IonType][GrdPnt-GS_X];
							dPMFz=World->PMF[IonType][GrdPnt]-World->PMF[IonType][GrdPnt-GS_XY];
							
							if((dPMFx*dPMFx)>4.0&&World->NIndexing->GetDiffDouble(IonType,GrdPnt-1)!=0.0)
								countPMFx++;
							if((dPMFy*dPMFy)>4.0&&World->NIndexing->GetDiffDouble(IonType,GrdPnt-GS_X)!=0.0)
								countPMFy++;
							if((dPMFz*dPMFz)>4.0&&World->NIndexing->GetDiffDouble(IonType,GrdPnt-GS_XY)!=0.0)
								countPMFz++;
						}
					}
				}
			}
			pnpPrint0("\tIon %d: countPMFx=%d,countPMFy=%d,countPMFz=%d Throught D!=0\n"
					,IonType,countPMFx,countPMFy,countPMFz);
		}
		pnpPrint0("PMF checked\n");
	}
	delete iV;
	return EXIT_SUCCESS;
}
int NernstPlankSolverDouble::NernstPlanckSolverD()
{
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	double gridScale;
	int iteration;
	int i,j,k,i1;
	int gridPoint;
	double om1,om2;
	double * potential;
	//double * diffusion;
	int gridSizeXY;
	int gridSizeXYZ;
	double change;
	double fpoh;
	//double * positiveCharge;
	//double * negativeCharge;
	double rodlt1,rodlt2;
	double phi0,phi1,phi2,phi3,phi4,phi5,phi6;
	double rophi1,rophi2;
	double term1,term2;
	double phit;
	double phit1,phit2;
	//double denom1,denom2;
	int gridpx,gridmx;
	int gridpy,gridmy;
	int gridpz,gridmz;
	double conv;
	bool *PeriodicBoundaryCondition;
	int itmp1,itmp2;
	double *C[2];
	int IonType;
	
	//CheckSystem();

	if(!World){
		cerr << "ERROR 311: No World\n";
		exit(104);
	}
	
	if(verbose&&World->MyRank==0)
	{
		//fprintf(stdout,"\nNernst-Plank Solver\n");
	}
	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;	
	//diffusion = World->D[0];
	potential = World->PotentialDouble;
	//positiveCharge = World->C[0];
	//negativeCharge = World->C[1];
	C[0] = World->CDouble[0];
	C[1] = World->CDouble[1];
	PeriodicBoundaryCondition=World->PBC;
	
	om1 = Relaxation;
	om2 = 1-om1;

	fpoh = 4*M_PI*World->GridScale;
	conv = (World->GridScale*World->GridScale*World->GridScale)/COANGS;

	for(i=0;i<World->NIonsTypes;i++){
		if(UTMP[i] == NULL)if(!(UTMP[i] = new double[gridSizeXYZ])){
			cerr << "ERROR 104: No memory available\n";
			exit(104);
		}
		if(TMP[i] == NULL)
		{
			if(!(TMP[i] = new double[gridSizeXYZ])){
				cerr << "ERROR 104: No memory available\n";
				exit(104);
			}
			for(j=0;j<gridSizeXYZ;j++)TMP[i][j]=0.0f;
		}
	}
	for(j=0;j<World->NIonsTypes;j++)
	{
		for(i=0;i<6;i++)
			if(dphi[j][i] == NULL)
				if(!(dphi[j][i] = new double[SingularNum[j][2]]))
		{
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
	}
	//Set up the Potential energy
	if(World->PMF!=NULL)
	{
		if(potential!=NULL)
		{
			if(PMFWeight>=0.0)
			{
				for(i=0;i<World->NIonsTypes;i++)
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[i][j]=0.5*(World->IonsQ[i]*potential[j]+World->PMF[i][j]*PMFWeight);
			}
			else
			{
				for(i=0;i<World->NIonsTypes;i++)
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[i][j]=0.5*(World->IonsQ[i]*potential[j]+World->PMF[i][j]);
			}
		}
		else
		{
			if(PMFWeight>=0.0)
			{
				for(i=0;i<World->NIonsTypes;i++)
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[i][j]=0.5*(World->PMF[i][j])*PMFWeight;
			}
			else
			{
				for(i=0;i<World->NIonsTypes;i++)
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[i][j]=0.5*(World->PMF[i][j]);
			}
		}
	}
	else
		for(i=0;i<World->NIonsTypes;i++)
			for(j=0;j<gridSizeXYZ;j++)UTMP[i][j]=0.5*(World->IonsQ[i]*potential[j]);

	for(IonType=0;IonType<World->NIonsTypes;IonType++)
	{
		for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)C[IonType][gridPoint]/=fpoh;

		for(i1=0;i1<NoSingularNum[IonType][2];i1++)
		{
			gridPoint=IndexNoSingular[IonType][i1];
			gridpx = gridPoint+1;
			gridmx = gridPoint-1;
			gridpy = gridPoint+gridSizeX;
			gridmy = gridPoint-gridSizeX;
			gridpz = gridPoint+gridSizeXY;
			gridmz = gridPoint-gridSizeXY;
	
			phit = UTMP[IonType][gridpx] + UTMP[IonType][gridmx] + UTMP[IonType][gridpy] + UTMP[IonType][gridmy] + UTMP[IonType][gridpz] + UTMP[IonType][gridmz] - 6.0*UTMP[IonType][gridPoint];
			TMP[IonType][gridPoint]=1.0/(6.0-phit);
		}
		for(i=0;i<SingularNum[IonType][2];i++)
		{
			gridPoint=IndexSingular[IonType][i];
			gridpx = gridPoint+1;
			gridmx = gridPoint-1;
			gridpy = gridPoint+gridSizeX;
			gridmy = gridPoint-gridSizeX;
			gridpz = gridPoint+gridSizeXY;
			gridmz = gridPoint-gridSizeXY;
		
			dphi[IonType][0][i] = dix[IonType][0][i]*UTMP[IonType][gridpx];
			dphi[IonType][1][i] = dix[IonType][1][i]*UTMP[IonType][gridmx];
			dphi[IonType][2][i] = dix[IonType][2][i]*UTMP[IonType][gridpy];
			dphi[IonType][3][i] = dix[IonType][3][i]*UTMP[IonType][gridmy];
			dphi[IonType][4][i] = dix[IonType][4][i]*UTMP[IonType][gridpz];
			dphi[IonType][5][i] = dix[IonType][5][i]*UTMP[IonType][gridmz];
			phit = dphi[IonType][0][i] + dphi[IonType][1][i] + dphi[IonType][2][i] + dphi[IonType][3][i] + dphi[IonType][4][i] + dphi[IonType][5][i] - dixt[IonType][i]*UTMP[IonType][gridPoint];
			TMP[IonType][gridPoint] = 1.0/(dixt[IonType][i]-phit);
		}
	}
//	 VectorField3D VF1(World->GridSize,World->GridScale,2,UTMP);
//	 VF1.WriteToFile("UTMP.bin");
//	 VectorField3D VF2(World->GridSize,World->GridScale,2,TMP);
//	 VF2.WriteToFile("TMP.bin");
//	 VectorField3D VF3(World->GridSize,World->GridScale,2,World->PMF);
//	 VF3.WriteToFile("RF.bin");
//	 VectorField3D VI(World->GridSize,World->GridScale,2);
//	 for(IonType=0;IonType<World->NIonsTypes;IonType++)
//	 {
//		 for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
//			 VI.V[IonType][gridPoint]=0;
	// 
//		 for(i1=0;i1<NoSingularNum[IonType][2];i1++)
//		 {
//			 gridPoint=IndexNoSingular[IonType][i1];
//			 VI.V[IonType][gridPoint]=1;
//		 }
//		 for(i=0;i<SingularNum[IonType][2];i++)
//		 {
//			 gridPoint=IndexSingular[IonType][i];
//			 VI.V[IonType][gridPoint]=2;
//		 }
//	 }
//	 VI.WriteToFile("ind.bin");
	for(iteration=1;iteration<=MaxIterations;iteration++)
	{
		if((verbose&&(iteration%ConvergenceCheck==0))||(iteration==MaxIterations))
		{
			MaxChange = 0;
			for(j=0;j<2;j++)
			{
				for(IonType=0;IonType<World->NIonsTypes;IonType++)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
					
						rodlt1 = C[IonType][gridpx] + C[IonType][gridmx] + C[IonType][gridpy] + C[IonType][gridmy] + C[IonType][gridpz] + C[IonType][gridmz];
						change = ( rodlt1 * ( 1.0 - UTMP[IonType][gridPoint] ) + C[IonType][gridpx]*UTMP[IonType][gridpx] + C[IonType][gridmx]*UTMP[IonType][gridmx] + C[IonType][gridpy]*UTMP[IonType][gridpy] + C[IonType][gridmy]*UTMP[IonType][gridmy] + C[IonType][gridpz]*UTMP[IonType][gridpz] + C[IonType][gridmz]*UTMP[IonType][gridmz]) * TMP[IonType][gridPoint];
						change -= C[IonType][gridPoint];
						change *= om1;
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						C[IonType][gridPoint]+=change;
						
//						 if(C[IonType][gridPoint]<0.0)
//						 {
//							 fprintf(stderr,"C[%d][%d]=%g<0.0\n",IonType,gridPoint,[IonType][gridPoint]);
//							 return EXIT_FAILURE;
//						 }
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
						
						rodlt1 = C[IonType][gridpx]*dix[IonType][0][i] + C[IonType][gridmx]*dix[IonType][1][i] + C[IonType][gridpy]*dix[IonType][2][i] + C[IonType][gridmy]*dix[IonType][3][i] + C[IonType][gridpz]*dix[IonType][4][i] + C[IonType][gridmz]*dix[IonType][5][i];
						rophi1 = C[IonType][gridpx]*dphi[IonType][0][i] + C[IonType][gridmx]*dphi[IonType][1][i] + C[IonType][gridpy]*dphi[IonType][2][i] + C[IonType][gridmy]*dphi[IonType][3][i] + C[IonType][gridpz]*dphi[IonType][4][i] + C[IonType][gridmz]*dphi[IonType][5][i];
						term1 = rodlt1*(1-UTMP[IonType][gridPoint])+rophi1;
						change = om1*(term1*TMP[IonType][gridPoint]-C[IonType][gridPoint]);
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						C[IonType][gridPoint]+=change;
					}
				}
				for(i1=0;i1<World->NIonsTypes;i1++)
					World->BorderExchangeDouble(C[i1]);
			}
			MaxChange*=conv;
#ifdef MPI_PARALLEL
			int dest;
			pnpsapp->MyComGroup.Barrier();
			
			if(World->MyRank==0)
			{
				for(dest=1;dest<World->NProcs;dest++)
				{
					pnpsapp->MyComGroup.Recv(&change, 1, MPI::DOUBLE, dest, 0);
					MaxChange = MaxChange>change?MaxChange:change;
					MaxChange = -MaxChange<change?MaxChange:-change;
				}
			}
			else
			{
				pnpsapp->MyComGroup.Send(&MaxChange, 1, MPI::DOUBLE, 0, 0);
			}
#endif
			if(verbose)if(World->MyRank==0)
			{
				fprintf(stdout, "<NernstPlankIteration Nit=\"%6d\"	dC=\"%16.8g\"/>\n", iteration, MaxChange);
			}
#if defined(_MSC_VER) || defined(__DECCXX) || (__GNUC__ < 3)
			if(MaxChange>1E13)
#else
				if(MaxChange>1E13||MaxChange==NAN)
#endif
			{
				for(i1=0;i1<World->NIonsTypes;i1++)for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
						C[i1][gridPoint]*=fpoh;
				return EXIT_FAILURE;
			}
			if(MaxChange<Convergence)
				iteration = MaxIterations+1;
		}
		else
		{
			for(j=0;j<2;j++)
			{
				for(IonType=0;IonType<World->NIonsTypes;IonType++)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
						rodlt1 = C[IonType][gridpx] + C[IonType][gridmx] + C[IonType][gridpy] + C[IonType][gridmy] + C[IonType][gridpz] + C[IonType][gridmz];
						change = ( rodlt1 * ( 1.0 - UTMP[IonType][gridPoint] ) + C[IonType][gridpx]*UTMP[IonType][gridpx] + C[IonType][gridmx]*UTMP[IonType][gridmx] + C[IonType][gridpy]*UTMP[IonType][gridpy] + C[IonType][gridmy]*UTMP[IonType][gridmy] + C[IonType][gridpz]*UTMP[IonType][gridpz] + C[IonType][gridmz]*UTMP[IonType][gridmz]) * TMP[IonType][gridPoint];
						change -= C[IonType][gridPoint];
						change *= om1;
						C[IonType][gridPoint]+=change;
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++){
						gridPoint=IndexSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
	
						rodlt1 = C[IonType][gridpx]*dix[IonType][0][i] + C[IonType][gridmx]*dix[IonType][1][i] + C[IonType][gridpy]*dix[IonType][2][i] + C[IonType][gridmy]*dix[IonType][3][i] + C[IonType][gridpz]*dix[IonType][4][i] + C[IonType][gridmz]*dix[IonType][5][i];
						rophi1 = C[IonType][gridpx]*dphi[IonType][0][i] + C[IonType][gridmx]*dphi[IonType][1][i] + C[IonType][gridpy]*dphi[IonType][2][i] + C[IonType][gridmy]*dphi[IonType][3][i] + C[IonType][gridpz]*dphi[IonType][4][i] + C[IonType][gridmz]*dphi[IonType][5][i];
						term1 = rodlt1*(1-UTMP[IonType][gridPoint])+rophi1;
						change = om1*(term1*TMP[IonType][gridPoint]-C[IonType][gridPoint]);
						C[IonType][gridPoint]+=change;
					}
				}
				for(i1=0;i1<World->NIonsTypes;i1++)
					World->BorderExchangeDouble(C[i1]);
			}
		}
	}

	for(i1=0;i1<World->NIonsTypes;i1++)
		for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
			C[i1][gridPoint]*=fpoh;

	return EXIT_SUCCESS;
}
int NernstPlankSolverDouble::NernstPlanckSolverDOld()
{
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	double gridScale;
	int iteration;
	int i,j,k,i1;
	int gridPoint;
	double om1,om2;
	double * potential;
	//double * diffusion;
	int gridSizeXY;
	int gridSizeXYZ;
	double change;
	double fpoh;
	//double * positiveCharge;
	//double * negativeCharge;
	double rodlt1,rodlt2;
	double phi0,phi1,phi2,phi3,phi4,phi5,phi6;
	double rophi1,rophi2;
	double term1,term2;
	double phit;
	double phit1,phit2;
	//double denom1,denom2;
	int gridpx,gridmx;
	int gridpy,gridmy;
	int gridpz,gridmz;
	double conv;
	bool *PeriodicBoundaryCondition;
	int itmp1,itmp2;
	double *C[2];
	int IonType;
	
	//CheckSystem();

	if(!World){
		cerr << "ERROR 311: No World\n";
		exit(104);
	}
	
	if(verbose&&World->MyRank==0)
	{
		//fprintf(stdout,"\nNernst-Plank Solver\n");
	}
	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;	
	//diffusion = World->D[0];
	potential = World->PotentialDouble;
	//positiveCharge = World->C[0];
	//negativeCharge = World->C[1];
	C[0] = World->CDouble[0];
	C[1] = World->CDouble[1];
	PeriodicBoundaryCondition=World->PBC;
	
	om1 = Relaxation;
	om2 = 1-om1;

	fpoh = 4*M_PI*World->GridScale;
	conv = (World->GridScale*World->GridScale*World->GridScale)/COANGS;

	for(i=0;i<World->NIonsTypes;i++){
		if(UTMP[i] == NULL)if(!(UTMP[i] = new double[gridSizeXYZ])){
			cerr << "ERROR 104: No memory available\n";
			exit(104);
		}
		if(TMP[i] == NULL)
		{
			if(!(TMP[i] = new double[gridSizeXYZ])){
				cerr << "ERROR 104: No memory available\n";
				exit(104);
			}
			for(j=0;j<gridSizeXYZ;j++)TMP[i][j]=0.0f;
		}
	}
	for(j=0;j<World->NIonsTypes;j++)
	{
		for(i=0;i<6;i++)
			if(dphi[j][i] == NULL)
				if(!(dphi[j][i] = new double[SingularNum[j][2]]))
		{
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
	}
	//Set up the Potential energy
	if(World->PMF!=NULL)
	{
		if(potential!=NULL)
		{
			if(PMFWeight>=0.0)
			{
				for(i=0;i<World->NIonsTypes;i++)
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[i][j]=0.5*(World->IonsQ[i]*potential[j]+World->PMF[i][j]*PMFWeight);
			}
			else
			{
				for(i=0;i<World->NIonsTypes;i++)
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[i][j]=0.5*(World->IonsQ[i]*potential[j]+World->PMF[i][j]);
			}
		}
		else
		{
			if(PMFWeight>=0.0)
			{
				for(i=0;i<World->NIonsTypes;i++)
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[i][j]=0.5*(World->PMF[i][j])*PMFWeight;
			}
			else
			{
				for(i=0;i<World->NIonsTypes;i++)
					for(j=0;j<gridSizeXYZ;j++)
						UTMP[i][j]=0.5*(World->PMF[i][j]);
			}
		}
	}
	else
		for(i=0;i<World->NIonsTypes;i++)
			for(j=0;j<gridSizeXYZ;j++)UTMP[i][j]=0.5*(World->IonsQ[i]*potential[j]);

	for(IonType=0;IonType<World->NIonsTypes;IonType++)
	{
		for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)C[IonType][gridPoint]/=fpoh;

		for(i1=0;i1<NoSingularNum[IonType][2];i1++)
		{
			gridPoint=IndexNoSingular[IonType][i1];
			gridpx = gridPoint+1;
			gridmx = gridPoint-1;
			gridpy = gridPoint+gridSizeX;
			gridmy = gridPoint-gridSizeX;
			gridpz = gridPoint+gridSizeXY;
			gridmz = gridPoint-gridSizeXY;
	
			phit = UTMP[IonType][gridpx] + UTMP[IonType][gridmx] + UTMP[IonType][gridpy] + UTMP[IonType][gridmy] + UTMP[IonType][gridpz] + UTMP[IonType][gridmz] - 6.0*UTMP[IonType][gridPoint];
			TMP[IonType][gridPoint]=1.0/(6.0-phit);
		}
		for(i=0;i<SingularNum[IonType][2];i++)
		{
			gridPoint=IndexSingular[IonType][i];
			gridpx = gridPoint+1;
			gridmx = gridPoint-1;
			gridpy = gridPoint+gridSizeX;
			gridmy = gridPoint-gridSizeX;
			gridpz = gridPoint+gridSizeXY;
			gridmz = gridPoint-gridSizeXY;
		
			dphi[IonType][0][i] = dix[IonType][0][i]*UTMP[IonType][gridpx];
			dphi[IonType][1][i] = dix[IonType][1][i]*UTMP[IonType][gridmx];
			dphi[IonType][2][i] = dix[IonType][2][i]*UTMP[IonType][gridpy];
			dphi[IonType][3][i] = dix[IonType][3][i]*UTMP[IonType][gridmy];
			dphi[IonType][4][i] = dix[IonType][4][i]*UTMP[IonType][gridpz];
			dphi[IonType][5][i] = dix[IonType][5][i]*UTMP[IonType][gridmz];
			phit = dphi[IonType][0][i] + dphi[IonType][1][i] + dphi[IonType][2][i] + dphi[IonType][3][i] + dphi[IonType][4][i] + dphi[IonType][5][i] - dixt[IonType][i]*UTMP[IonType][gridPoint];
			TMP[IonType][gridPoint] = 1.0/(dixt[IonType][i]-phit);
		}
	}
//	 VectorField3D VF1(World->GridSize,World->GridScale,2,UTMP);
//	 VF1.WriteToFile("UTMP.bin");
//	 VectorField3D VF2(World->GridSize,World->GridScale,2,TMP);
//	 VF2.WriteToFile("TMP.bin");
//	 VectorField3D VF3(World->GridSize,World->GridScale,2,World->PMF);
//	 VF3.WriteToFile("RF.bin");
//	 VectorField3D VI(World->GridSize,World->GridScale,2);
//	 for(IonType=0;IonType<World->NIonsTypes;IonType++)
//	 {
//		 for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
//			 VI.V[IonType][gridPoint]=0;
	// 
//		 for(i1=0;i1<NoSingularNum[IonType][2];i1++)
//		 {
//			 gridPoint=IndexNoSingular[IonType][i1];
//			 VI.V[IonType][gridPoint]=1;
//		 }
//		 for(i=0;i<SingularNum[IonType][2];i++)
//		 {
//			 gridPoint=IndexSingular[IonType][i];
//			 VI.V[IonType][gridPoint]=2;
//		 }
//	 }
//	 VI.WriteToFile("ind.bin");
	for(iteration=1;iteration<=MaxIterations;iteration++)
	{
		if((verbose&&(iteration%ConvergenceCheck==0))||(iteration==MaxIterations))
		{
			MaxChange = 0;
			for(j=0;j<2;j++)
			{
				for(IonType=0;IonType<World->NIonsTypes;IonType++)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
					
						rodlt1 = C[IonType][gridpx] + C[IonType][gridmx] + C[IonType][gridpy] + C[IonType][gridmy] + C[IonType][gridpz] + C[IonType][gridmz];
						change = ( rodlt1 * ( 1.0 - UTMP[IonType][gridPoint] ) + C[IonType][gridpx]*UTMP[IonType][gridpx] + C[IonType][gridmx]*UTMP[IonType][gridmx] + C[IonType][gridpy]*UTMP[IonType][gridpy] + C[IonType][gridmy]*UTMP[IonType][gridmy] + C[IonType][gridpz]*UTMP[IonType][gridpz] + C[IonType][gridmz]*UTMP[IonType][gridmz]) * TMP[IonType][gridPoint];
						change -= C[IonType][gridPoint];
						change *= om1;
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						C[IonType][gridPoint]+=change;
						
//						 if(C[IonType][gridPoint]<0.0)
//						 {
//							 fprintf(stderr,"C[%d][%d]=%g<0.0\n",IonType,gridPoint,[IonType][gridPoint]);
//							 return EXIT_FAILURE;
//						 }
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
						
						rodlt1 = C[IonType][gridpx]*dix[IonType][0][i] + C[IonType][gridmx]*dix[IonType][1][i] + C[IonType][gridpy]*dix[IonType][2][i] + C[IonType][gridmy]*dix[IonType][3][i] + C[IonType][gridpz]*dix[IonType][4][i] + C[IonType][gridmz]*dix[IonType][5][i];
						rophi1 = C[IonType][gridpx]*dphi[IonType][0][i] + C[IonType][gridmx]*dphi[IonType][1][i] + C[IonType][gridpy]*dphi[IonType][2][i] + C[IonType][gridmy]*dphi[IonType][3][i] + C[IonType][gridpz]*dphi[IonType][4][i] + C[IonType][gridmz]*dphi[IonType][5][i];
						term1 = rodlt1*(1-UTMP[IonType][gridPoint])+rophi1;
						change = om1*(term1*TMP[IonType][gridPoint]-C[IonType][gridPoint]);
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						C[IonType][gridPoint]+=change;
					}
				}
				for(i1=0;i1<World->NIonsTypes;i1++)
					World->BorderExchangeDouble(C[i1]);
			}
			MaxChange*=conv;
#ifdef MPI_PARALLEL
			int dest;
			pnpsapp->MyComGroup.Barrier();
			
			if(World->MyRank==0)
			{
				for(dest=1;dest<World->NProcs;dest++)
				{
					pnpsapp->MyComGroup.Recv(&change, 1, MPI::FLOAT, dest, 0);
					MaxChange = MaxChange>change?MaxChange:change;
					MaxChange = -MaxChange<change?MaxChange:-change;
				}
			}
			else
			{
				pnpsapp->MyComGroup.Send(&MaxChange, 1, MPI::FLOAT, 0, 0);
			}
#endif
			if(verbose)if(World->MyRank==0)
			{
				fprintf(stdout, "<NernstPlankIteration Nit=\"%6d\"	dC=\"%16.8g\"/>\n", iteration, MaxChange);
			}
#if defined(_MSC_VER) || defined(__DECCXX) || (__GNUC__ < 3)
			if(MaxChange>1E13)
#else
				if(MaxChange>1E13||MaxChange==NAN)
#endif
			{
				for(i1=0;i1<World->NIonsTypes;i1++)for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
						C[i1][gridPoint]*=fpoh;
				return EXIT_FAILURE;
			}
			if(MaxChange<Convergence)
				iteration = MaxIterations+1;
		}
		else
		{
			for(j=0;j<2;j++)
			{
				for(IonType=0;IonType<World->NIonsTypes;IonType++)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
						rodlt1 = C[IonType][gridpx] + C[IonType][gridmx] + C[IonType][gridpy] + C[IonType][gridmy] + C[IonType][gridpz] + C[IonType][gridmz];
						change = ( rodlt1 * ( 1.0 - UTMP[IonType][gridPoint] ) + C[IonType][gridpx]*UTMP[IonType][gridpx] + C[IonType][gridmx]*UTMP[IonType][gridmx] + C[IonType][gridpy]*UTMP[IonType][gridpy] + C[IonType][gridmy]*UTMP[IonType][gridmy] + C[IonType][gridpz]*UTMP[IonType][gridpz] + C[IonType][gridmz]*UTMP[IonType][gridmz]) * TMP[IonType][gridPoint];
						change -= C[IonType][gridPoint];
						change *= om1;
						C[IonType][gridPoint]+=change;
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++){
						gridPoint=IndexSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
	
						rodlt1 = C[IonType][gridpx]*dix[IonType][0][i] + C[IonType][gridmx]*dix[IonType][1][i] + C[IonType][gridpy]*dix[IonType][2][i] + C[IonType][gridmy]*dix[IonType][3][i] + C[IonType][gridpz]*dix[IonType][4][i] + C[IonType][gridmz]*dix[IonType][5][i];
						rophi1 = C[IonType][gridpx]*dphi[IonType][0][i] + C[IonType][gridmx]*dphi[IonType][1][i] + C[IonType][gridpy]*dphi[IonType][2][i] + C[IonType][gridmy]*dphi[IonType][3][i] + C[IonType][gridpz]*dphi[IonType][4][i] + C[IonType][gridmz]*dphi[IonType][5][i];
						term1 = rodlt1*(1-UTMP[IonType][gridPoint])+rophi1;
						change = om1*(term1*TMP[IonType][gridPoint]-C[IonType][gridPoint]);
						C[IonType][gridPoint]+=change;
					}
				}
				for(i1=0;i1<World->NIonsTypes;i1++)
					World->BorderExchangeDouble(C[i1]);
			}
		}
	}

	for(i1=0;i1<World->NIonsTypes;i1++)
		for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
			C[i1][gridPoint]*=fpoh;

	return EXIT_SUCCESS;
}
int NernstPlankSolverDouble::NernstPlanckSolverN2tmp()
{
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	double gridScale;
	int iteration;
	int i,j,k,i1;
	int gridPoint;
	double om1,om2;
	double * potential;
	//double * diffusion;
	int gridSizeXY;
	int gridSizeXYZ;
	double change;
	double fpoh;
	//double * positiveCharge;
	//double * negativeCharge;
	double rodlt1,rodlt2;
	double phi0,phi1,phi2,phi3,phi4,phi5,phi6;
	double rophi1,rophi2;
	double term1,term2;
	double phit;
	double phit1,phit2;
	//double denom1,denom2;
	int gridpx,gridmx;
	int gridpy,gridmy;
	int gridpz,gridmz;
	double conv;
	bool *PeriodicBoundaryCondition;
	int itmp1,itmp2;
	double *C[2];
	int IonType;
	double _MaxChange=0.0;

	
	if(!World){
		cerr << "ERROR 311: No World\n";
		exit(104);
	}
	
	if(verbose&&World->MyRank==0)
	{
		//fprintf(stdout,"\nNernst-Plank Solver\n");
	}
	gridScale = World->GridScale;
	gridSizeX = World->GridSize[0];
	gridSizeY = World->GridSize[1];
	gridSizeZ = World->GridSize[2];
	gridSizeXY = gridSizeX*gridSizeY;
	gridSizeXYZ = gridSizeXY*gridSizeZ;	
	//diffusion = World->D[0];
	potential = World->PotentialDouble;
	//positiveCharge = World->C[0];
	//negativeCharge = World->C[1];
	C[0] = World->CDouble[0];
	C[1] = World->CDouble[1];
	PeriodicBoundaryCondition=World->PBC;
	
	om1 = Relaxation;
	om2 = 1-om1;

	fpoh = 4*M_PI*World->GridScale;
	conv = (World->GridScale*World->GridScale*World->GridScale)/COANGS;

	
	if(UTMP == NULL)if(!(UTMP = new double[gridSizeXYZ])){
		cerr << "ERROR 104: No memory available\n";
		exit(104);
	}
	if(TMP == NULL)if(!(TMP = new double[gridSizeXYZ])){
		cerr << "ERROR 104: No memory available\n";
		exit(104);
	}
	double *utmp=UTMP;
	double *tmp=TMP;
	double *Ci;
	double *dphi0,*dphi1,*dphi2,*dphi3,*dphi4,*dphi5;
	for(j=0;j<World->NIonsTypes;j++)
	{
		for(i=0;i<6;i++)if(dphi[j][i] == NULL)if(!(dphi[j][i] = new double[SingularNum[j][2]]))
		{
			cerr << "ERROR 204: No memory available\n";
			exit(204);
		}
	}
	//Set up the Potential energy
	for(IonType=0;IonType<World->NIonsTypes;IonType++)
	{
		dphi0=dphi[IonType][0];
		dphi1=dphi[IonType][1];
		dphi2=dphi[IonType][2];
		dphi3=dphi[IonType][3];
		dphi4=dphi[IonType][4];
		dphi5=dphi[IonType][5];
		Ci=C[IonType];
		
		if(World->PMF!=NULL)
			for(j=0;j<gridSizeXYZ;j++)utmp[j]=0.5*(World->IonsQ[IonType]*potential[j]+World->PMF[IonType][j]);
		else
			for(j=0;j<gridSizeXYZ;j++)utmp[j]=0.5*(World->IonsQ[IonType]*potential[j]);

	
		for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)C[IonType][gridPoint]/=fpoh;

		for(i1=0;i1<NoSingularNum[IonType][2];i1++)
		{
			gridPoint=IndexNoSingular[IonType][i1];
			gridpx = gridPoint+1;
			gridmx = gridPoint-1;
			gridpy = gridPoint+gridSizeX;
			gridmy = gridPoint-gridSizeX;
			gridpz = gridPoint+gridSizeXY;
			gridmz = gridPoint-gridSizeXY;
	
			phit = utmp[gridpx] + utmp[gridmx] + utmp[gridpy] + utmp[gridmy] + utmp[gridpz] + utmp[gridmz] - 6.0*utmp[gridPoint];
			tmp[gridPoint]=1.0/(6.0-phit);
		}
		for(i=0;i<SingularNum[IonType][2];i++)
		{
			gridPoint=IndexSingular[IonType][i];
			gridpx = gridPoint+1;
			gridmx = gridPoint-1;
			gridpy = gridPoint+gridSizeX;
			gridmy = gridPoint-gridSizeX;
			gridpz = gridPoint+gridSizeXY;
			gridmz = gridPoint-gridSizeXY;
		
			dphi0[i] = dix[IonType][0][i]*utmp[gridpx];
			dphi1[i] = dix[IonType][1][i]*utmp[gridmx];
			dphi2[i] = dix[IonType][2][i]*utmp[gridpy];
			dphi3[i] = dix[IonType][3][i]*utmp[gridmy];
			dphi4[i] = dix[IonType][4][i]*utmp[gridpz];
			dphi5[i] = dix[IonType][5][i]*utmp[gridmz];
			phit = dphi0[i] + dphi1[i] + dphi2[i] + dphi3[i] + dphi4[i] + dphi5[i] - dixt[IonType][i]*utmp[gridPoint];
			tmp[gridPoint] = 1.0/(dixt[IonType][i]-phit);
		}
		for(iteration=1;iteration<=MaxIterations;iteration++)
		{
			if((verbose&&(iteration%ConvergenceCheck==0))||(iteration==MaxIterations))
			{
				MaxChange = 0;
				for(j=0;j<2;j++)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
					
						rodlt1 = Ci[gridpx] + Ci[gridmx] + Ci[gridpy] + Ci[gridmy] + Ci[gridpz] + Ci[gridmz];
						change = ( rodlt1 * ( 1.0 - utmp[gridPoint] ) + Ci[gridpx]*utmp[gridpx] + Ci[gridmx]*utmp[gridmx] + Ci[gridpy]*utmp[gridpy] + Ci[gridmy]*utmp[gridmy] + Ci[gridpz]*utmp[gridpz] + Ci[gridmz]*utmp[gridmz]) * tmp[gridPoint];
						change -= Ci[gridPoint];
						change *= om1;
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						Ci[gridPoint]+=change;
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
						
						rodlt1 = Ci[gridpx]*dix[IonType][0][i] + Ci[gridmx]*dix[IonType][1][i] + Ci[gridpy]*dix[IonType][2][i] + Ci[gridmy]*dix[IonType][3][i] + Ci[gridpz]*dix[IonType][4][i] + Ci[gridmz]*dix[IonType][5][i];
						rophi1 = Ci[gridpx]*dphi0[i] + Ci[gridmx]*dphi1[i] + Ci[gridpy]*dphi2[i] + Ci[gridmy]*dphi3[i] + Ci[gridpz]*dphi4[i] + Ci[gridmz]*dphi5[i];
						term1 = rodlt1*(1-utmp[gridPoint])+rophi1;
						change = om1*(term1*tmp[gridPoint]-Ci[gridPoint]);
						MaxChange = MaxChange>change?MaxChange:change;
						MaxChange = -MaxChange<change?MaxChange:-change;
						Ci[gridPoint]+=change;
					}
					World->BorderExchangeDouble(Ci);
				}
				MaxChange*=conv;
				if(verbose||iteration==MaxIterations)
				{
					if(World->NProcs!=1){
					//int tag=SEND_ENERGY;
					/*MPI::COMM_WORLD.Barrier();
						MPI::Status	status;
						double MaxChangeProc;
						if(World->MyRank==0){
						for(i=1;i<World->NProcs;i++){
						MPI::COMM_WORLD.Recv(&MaxChangeProc, 1, MPI::FLOAT, i, SEND_MAXCHANGE, status);
						MaxChange = MaxChange>MaxChangeProc?MaxChange:MaxChangeProc;
					}
					}
						else{
						MPI::COMM_WORLD.Send(&MaxChange, 1, MPI::FLOAT, 0, SEND_MAXCHANGE);
					}
						MPI::COMM_WORLD.Barrier();*/
					}
					if(verbose)if(World->MyRank==0){
						fprintf(stdout, "<NernstPlankIteration Nit=\"%6d\"	dC[%d]=\"%16.8g\"/>\n", iteration, IonType,MaxChange);
		//MaxChange=0;
					}
#if defined(_MSC_VER) || defined(__DECCXX) || (__GNUC__ < 3)
					if(MaxChange>1E13)
#else
						if(MaxChange>1E13||MaxChange==NAN)
#endif
					{
						for(i1=0;i1<World->NIonsTypes;i1++)
							for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
								C[i1][gridPoint]*=fpoh;
						return EXIT_FAILURE;
					}
				}
	
				if(MaxChange<Convergence)
					iteration = MaxIterations+1;
			}
			else
			{
				for(j=0;j<2;j++)
				{
					for(i=NoSingularNum[IonType][j];i<NoSingularNum[IonType][j+1];i++)
					{
						gridPoint=IndexNoSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
						rodlt1 = Ci[gridpx] + Ci[gridmx] + Ci[gridpy] + Ci[gridmy] + Ci[gridpz] + Ci[gridmz];
						change = ( rodlt1 * ( 1.0 - utmp[gridPoint] ) + Ci[gridpx]*utmp[gridpx] + Ci[gridmx]*utmp[gridmx] + Ci[gridpy]*utmp[gridpy] + Ci[gridmy]*utmp[gridmy] + Ci[gridpz]*utmp[gridpz] + Ci[gridmz]*utmp[gridmz]) * tmp[gridPoint];
						change -= Ci[gridPoint];
						change *= om1;
						Ci[gridPoint]+=change;
					}
					for(i=SingularNum[IonType][j];i<SingularNum[IonType][j+1];i++){
						gridPoint=IndexSingular[IonType][i];
						gridpx = gridPoint+1;
						gridmx = gridPoint-1;
						gridpy = gridPoint+gridSizeX;
						gridmy = gridPoint-gridSizeX;
						gridpz = gridPoint+gridSizeXY;
						gridmz = gridPoint-gridSizeXY;
	
						rodlt1 = Ci[gridpx]*dix[IonType][0][i] + Ci[gridmx]*dix[IonType][1][i] + Ci[gridpy]*dix[IonType][2][i] + Ci[gridmy]*dix[IonType][3][i] + Ci[gridpz]*dix[IonType][4][i] + Ci[gridmz]*dix[IonType][5][i];
						rophi1 = Ci[gridpx]*dphi0[i] + Ci[gridmx]*dphi1[i] + Ci[gridpy]*dphi2[i] + Ci[gridmy]*dphi3[i] + Ci[gridpz]*dphi4[i] + Ci[gridmz]*dphi5[i];
						term1 = rodlt1*(1-utmp[gridPoint])+rophi1;
						change = om1*(term1*tmp[gridPoint]-Ci[gridPoint]);
						Ci[gridPoint]+=change;
					}
				}
				World->BorderExchangeDouble(Ci);
			}
		}
		if(MaxChange>_MaxChange)_MaxChange=MaxChange;
	}
	MaxChange=_MaxChange;
	for(i1=0;i1<World->NIonsTypes;i1++)
		for(gridPoint=0;gridPoint<gridSizeXYZ;gridPoint++)
			C[i1][gridPoint]*=fpoh;

	return EXIT_SUCCESS;
}
HaPyDoubleVectorField3D* NernstPlankSolverDouble::CalcJ()
{
	int * GS;
	int GS_X;
	int GS_XY;
	int GS_XYZ;
	int gridPoint;
	double gridScale;
	double fpoh;
	double flcon;
	double ftot;
	
	int i,j,k;
	int gridp,gridn;
	double pot;
	int gridInc;
	double *C[2];
	double * potential;
	
	
	NodeIndex* NIndex=World->NIndexing->NIndex;
	NodeIndexing* NIndexing=World->NIndexing;
	NodeIndex DiffBoarderMask=NodeIndexing::DiffIon0BoarderMask;
	NodeIndex BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;
	
	GS = World->GridSize;
	GS_X = GS[0];
	GS_XY = GS[0]*GS[1];
	GS_XYZ = GS[2]*GS_XY;
	gridScale = World->GridScale;
	
	C[0] = World->CDouble[0];
	C[1] = World->CDouble[1];
	potential = World->PotentialDouble;

	HaPyDoubleVectorField3D* J=new HaPyDoubleVectorField3D(GS,gridScale,8);
	double **V=J->V;
	if(World->PMF!=NULL)
		for(i=0;i<World->NIonsTypes;i++){
		for(j=0;j<GS_XYZ;j++)UTMP[i][j]=0.5*(World->IonsQ[i]*potential[j]+World->PMF[i][j]);
		}
		else
		{
			for(i=0;i<World->NIonsTypes;i++)
			{
				for(j=0;j<GS_XYZ;j++)UTMP[i][j]=0.5*World->IonsQ[i]*potential[j];
			}
		}
	
		fpoh = 4*M_PI*gridScale;
		flcon = 1.6*gridScale*gridScale*1000;
		ftot = flcon/fpoh;

//	 for(i=0;i<World->NIonsTypes;i++)
//		 for(gridPoint=0;gridPoint<GS_XYZ;gridPoint++)
//				 C[i][gridPoint]/=fpoh;
	//cout<<dim<<"fdfg\n";
		double d,dc,sumc;
		for(j=0;j<World->NIonsTypes;j++)
		{
			for(k=0;k<3;k++)
			{
				switch(k)
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
				ftot = -1.0*World->IonsQ[j]*flcon/fpoh;
				for(i=0;i<NoSingularNum[j][2];i++)
				{
					gridPoint=IndexNoSingular[j][i];
					gridp = gridPoint+gridInc;
					pot = UTMP[j][gridp]-UTMP[j][gridPoint];
					d=NIndexing->GetDiffDouble(j,gridPoint);
					dc=C[j][gridp]-C[j][gridPoint];
					sumc=(C[j][gridp]+C[j][gridPoint]);
					V[k+j*4][gridPoint]=ftot*d*(dc+pot*sumc);
				
				}
				for(i=0;i<SingularNum[j][2];i++)
				{
					gridPoint=IndexSingular[j][i];
					gridp = gridPoint+gridInc;
					pot = UTMP[j][gridp]-UTMP[j][gridPoint];
					d=0.5*(NIndexing->GetDiffDouble(j,gridPoint)+NIndexing->GetDiffDouble(j,gridp));
					dc=C[j][gridp]-C[j][gridPoint];
					sumc=(C[j][gridp]+C[j][gridPoint]);
					V[k+j*4][gridPoint]=ftot*d*(dc+pot*sumc);
				//V[k+j*4][gridPoint]=ftot*dix[j][2*k][i]*(C[j][gridp]-C[j][gridPoint]+pot*(C[j][gridp]+C[j][gridPoint]));
				}
			}
		}
		for(j=0;j<World->NIonsTypes;j++)
			for(gridPoint=0;gridPoint<GS_XYZ;gridPoint++)
		{
			V[3+j*4][gridPoint]=sqrt(V[j*4][gridPoint]*V[j*4][gridPoint]+V[1+j*4][gridPoint]*V[1+j*4][gridPoint]+V[2+j*4][gridPoint]*V[2+j*4][gridPoint]);
		}
		return J;
}
void NernstPlankSolverDouble::nernstPlanckSolverCurrentProfile(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile)
{
	if(solver==Auto)
	{
		if(World->NIndexing==NULL)
			nernstPlanckSolverCurrentProfileAD(dim, positiveCurrentProfile, negativeCurrentProfile);
		else if(World->D!=NULL)
			nernstPlanckSolverCurrentProfileDflex(dim, positiveCurrentProfile, negativeCurrentProfile);
		else nernstPlanckSolverCurrentProfileN(dim, positiveCurrentProfile, negativeCurrentProfile);
		
	}
	else
	{
		if(solver==NodeIndexBased)
			nernstPlanckSolverCurrentProfileN(dim, positiveCurrentProfile, negativeCurrentProfile);
		else if(solver==ArrayDirect)
			nernstPlanckSolverCurrentProfileAD(dim, positiveCurrentProfile, negativeCurrentProfile);
		else
			nernstPlanckSolverCurrentProfileN(dim, positiveCurrentProfile, negativeCurrentProfile);
	}
	
}
void NernstPlankSolverDouble::nernstPlanckSolverCurrentProfileN(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile)
{
	int * GS;
	int GS_X,GS_Y,GS_Z;
	int GS_XY;
	int GS_XYZ;
	int gridPoint;
	double gridScale;
	double fpoh;
	double flcon;
	
	int i,j,k,ion;
	int kgrid,jgrid;
	int GrdPnt;
	int gridp,gridn;
	double pot;
	int grid[3];
	int gridInc;
	double *C[2];
	double *I[2];
	double * potential;
	
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
	I[0] = positiveCurrentProfile;
	I[1] = negativeCurrentProfile;
	potential = World->Potential;


	
	
	fpoh = 4*M_PI*gridScale;
	flcon = 1.6*gridScale*gridScale*1000;
	flcon/=fpoh;

	gridInc = GS_XY;
	
	for(i=0;i<World->NIonsTypes;i++)
		for(j=0;j<GS_Z-1;j++)
			I[i][j]=0.0;

	if(UTMP == NULL)if(!(UTMP = new double[GS_XYZ])){
		cerr << "ERROR 104: No memory available\n";
		exit(104);
	}
	
	double d,dc,sc,I1,I2;
	double *utmp=UTMP;
	for(ion=0;ion<World->NIonsTypes;ion++)
	{	
		if(World->PMF!=NULL)
		{
			for(j=0;j<GS_XYZ;j++)utmp[j]=0.5*(World->IonsQ[ion]*potential[j]+World->PMF[ion][j]);
		}
		else
		{
			for(j=0;j<GS_XYZ;j++)utmp[j]=0.5*World->IonsQ[ion]*potential[j];
		}
		if(QmobMask==NULL)
		{
			for(k=0;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+GS_XY;
						if(NIndexing->GetDiffDouble(ion,gridPoint)>0.0)
						{
							pot = utmp[gridp]-utmp[gridPoint];
							d = NIndexing->GetDiffDouble(ion,gridPoint+GS_XY)>0.0?					0.5*(NIndexing->GetDiffDouble(ion,gridPoint)+NIndexing->GetDiffDouble(ion,gridPoint+GS_XY)):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							I[ion][k]+=d*(dc+pot*sc);
						}
					}
				}
			}
		}
		else
		{
			for(k=0;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+GS_XY;
						if(NIndexing->GetDiffDouble(ion,gridPoint)>0.0&&QmobMask[gridPoint]>0)
						{
							pot = utmp[gridp]-utmp[gridPoint];
							d = NIndexing->GetDiffDouble(ion,gridPoint+GS_XY)>0.0?					0.5*(NIndexing->GetDiffDouble(ion,gridPoint)+NIndexing->GetDiffDouble(ion,gridPoint+GS_XY)):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							I[ion][k]+=d*(dc+pot*sc);
						}
					}
				}
			}
		}
	}
	for(i=0;i<GS[dim]-1;i++)
	{
		positiveCurrentProfile[i]*=-flcon;
		negativeCurrentProfile[i]*=flcon;
	}
}
void NernstPlankSolverDouble::nernstPlanckSolverCurrentProfileAD(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile)
{
	int * GS;
	int GS_X,GS_Y,GS_Z;
	int GS_XY;
	int GS_XYZ;
	int gridPoint;
	double gridScale;
	double fpoh;
	double flcon;
	
	int i,j,k,ion;
	int kgrid,jgrid;
	int GrdPnt;
	int gridp,gridn;
	double pot;
	int grid[3];
	int gridInc;
	double *C[2];
	double *I[2];
	double * potential;
	double **D=World->D;
	
	GS = World->GridSize;
	GS_X = GS[0];
	GS_Y = GS[1];
	GS_Z = GS[2];
	GS_XY = GS[0]*GS[1];
	GS_XYZ = GS[2]*GS_XY;
	gridScale = World->GridScale;
	
	C[0] = World->C[0];
	C[1] = World->C[1];
	I[0] = positiveCurrentProfile;
	I[1] = negativeCurrentProfile;
	potential = World->Potential;


	fpoh = 4*M_PI*gridScale;
	flcon = 1.6*gridScale*gridScale*1000;
	flcon/=fpoh;

	gridInc = GS_XY;
	
	for(i=0;i<World->NIonsTypes;i++)
		for(j=0;j<GS_Z-1;j++)
			I[i][j]=0.0;

	if(UTMP == NULL)if(!(UTMP = new double[GS_XYZ])){
		cerr << "ERROR 104: No memory available\n";
		exit(104);
	}
	
	double d,dc,sc,I1,I2;
	double *utmp=UTMP;
	for(ion=0;ion<World->NIonsTypes;ion++)
	{	
		if(World->PMF!=NULL)
		{
			for(j=0;j<GS_XYZ;j++)utmp[j]=0.5*(World->IonsQ[ion]*potential[j]+World->PMF[ion][j]);
		}
		else
		{
			for(j=0;j<GS_XYZ;j++)utmp[j]=0.5*World->IonsQ[ion]*potential[j];
		}
		if(QmobMask==NULL)
		{
			for(k=0;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+GS_XY;
						if(D[ion][gridPoint]>0.0f)
						{
							pot = utmp[gridp]-utmp[gridPoint];
							d = D[ion][gridPoint+GS_XY]>0.0f?					0.5*(D[ion][gridPoint]+D[ion][gridPoint+GS_XY]):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							I[ion][k]+=d*(dc+pot*sc);
						}
					}
				}
			}
		}
		else
		{
			for(k=0;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+GS_XY;
						if(D[ion][gridPoint]>0.0f&&QmobMask[gridPoint]>0)
						{
							pot = utmp[gridp]-utmp[gridPoint];
							d = D[ion][gridPoint+GS_XY]>0.0f?					0.5*(D[ion][gridPoint]+D[ion][gridPoint+GS_XY]):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							I[ion][k]+=d*(dc+pot*sc);
						}
					}
				}
			}
		}
	}
	for(i=0;i<GS[dim]-1;i++)
	{
		positiveCurrentProfile[i]*=-flcon;
		negativeCurrentProfile[i]*=flcon;
	}
}
void NernstPlankSolverDouble::nernstPlanckSolverCurrentProfileDflex(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile)
{
	int * GS;
	int GS_X,GS_Y,GS_Z;
	int GS_XY;
	int GS_XYZ;
	int gridPoint;
	double gridScale;
	double fpoh;
	double flcon;
	
	int i,j,k,ion;
	int kgrid,jgrid;
	int GrdPnt;
	int gridp,gridn;
	double pot;
	int grid[3];
	int gridInc;
	double *C[2];
	double *I[2];
	double * potential;
	double **D=World->D;
	
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
	I[0] = positiveCurrentProfile;
	I[1] = negativeCurrentProfile;
	potential = World->Potential;


	
	
	fpoh = 4*M_PI*gridScale;
	flcon = 1.6*gridScale*gridScale*1000;
	flcon/=fpoh;

	gridInc = GS_XY;
	
	for(i=0;i<World->NIonsTypes;i++)
		for(j=0;j<GS_Z-1;j++)
			I[i][j]=0.0;

	if(UTMP == NULL)if(!(UTMP = new double[GS_XYZ])){
		cerr << "ERROR 104: No memory available\n";
		exit(104);
	}
	
	double d,dc,sc,I1,I2;
	double *utmp=UTMP;
	for(ion=0;ion<World->NIonsTypes;ion++)
	{	
		if(World->PMF!=NULL)
		{
			for(j=0;j<GS_XYZ;j++)utmp[j]=0.5*(World->IonsQ[ion]*potential[j]+World->PMF[ion][j]);
		}
		else
		{
			for(j=0;j<GS_XYZ;j++)utmp[j]=0.5*World->IonsQ[ion]*potential[j];
		}
		if(QmobMask==NULL)
		{
			for(k=0;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+GS_XY;
						if(NIndexing->GetDiffDouble(ion,gridPoint)>0.0)
						{
							pot = utmp[gridp]-utmp[gridPoint];
							d = NIndexing->GetDiffDouble(ion,gridPoint+GS_XY)>0.0?					0.5*(D[ion][gridPoint]+D[ion][gridPoint+GS_XY]):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							I[ion][k]+=d*(dc+pot*sc);
						}
					}
				}
			}
		}
		else
		{
			for(k=0;k<GS_Z-1;k++)
			{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
					jgrid = kgrid+j*GS_X;
					for(i=1;i<GS_X-1;i++)
					{
						gridPoint = jgrid+i;
						gridp = gridPoint+GS_XY;
						if(NIndexing->GetDiffDouble(ion,gridPoint)>0.0&&QmobMask[gridPoint]>0)
						{
							pot = utmp[gridp]-utmp[gridPoint];
							d = NIndexing->GetDiffDouble(ion,gridPoint+GS_XY)>0.0?					0.5*(D[ion][gridPoint]+D[ion][gridPoint+GS_XY]):0;
							dc = C[ion][gridp]-C[ion][gridPoint];
							sc = C[ion][gridp]+C[ion][gridPoint];
							I[ion][k]+=d*(dc+pot*sc);
						}
					}
				}
			}
		}
	}
	for(i=0;i<GS[dim]-1;i++)
	{
		positiveCurrentProfile[i]*=-flcon;
		negativeCurrentProfile[i]*=flcon;
	}
}
#endif


