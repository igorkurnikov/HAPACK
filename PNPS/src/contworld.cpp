//
// C++ Implementation: contworld
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

#include "contworld.h"
#include "pnpdebug.h"
#include "pnpconstants.h"
#include "mapio.h"
#include <string>
#include <sstream>
#include <tinyxml.h>
#include <zlib.h>
#include <math.h>


#ifdef HARLEM_MOD
	#include "hasurface.h"
	#include "halinalg.h"
#endif
#include <stdexcept>

using namespace std;
///////////////////////////////////////////////////////////////////////////////
MapsIOData::MapsIOData()
{
	ModeStr.push_back("r");
	ModeStr.push_back("w");
	ModeStr.push_back("a");
	
	AddGroupStringSuffix=false;
	AddGroupNumberSuffix=false;
}
MapsIOData::~MapsIOData()
{
}
///////////////////////////////////////////////////////////////////////////////
NodeIndexing::NodeIndexing()
{
	NIndex=NULL;
	Q=NULL;
	NIonsTypes=2;
	
	//setup fields
	int i;
	
	for(i=0;i<MaxIonTypes;i++)
	{
		
	}
	EpsilonField[0]=Epsilon0;
	EpsilonField[1]=Epsilon1;
	EpsilonField[2]=Epsilon2;
	EpsilonFieldSft[0]=Epsilon0Sft;
	EpsilonFieldSft[1]=Epsilon1Sft;
	EpsilonFieldSft[2]=Epsilon2Sft;
	IonField[0]=Ion0;
	IonField[1]=Ion1;
	IonField[2]=Ion2;
	IonField[3]=Ion3;
	IonFieldSft[0]=Ion0Sft;
	IonFieldSft[1]=Ion1Sft;
	IonFieldSft[2]=Ion2Sft;
	IonFieldSft[3]=Ion3Sft;
	DiffBoarderMask[0]=DiffIon0BoarderMask;
	DiffBoarderMask[1]=DiffIon1BoarderMask;
	DiffBoarderMask[2]=DiffIon2BoarderMask;
	DiffBoarderMask[3]=DiffIon3BoarderMask;
	DiffBoarderMaskSft[0]=DiffIon0BoarderMaskSft;
	DiffBoarderMaskSft[1]=DiffIon1BoarderMaskSft;
	DiffBoarderMaskSft[2]=DiffIon2BoarderMaskSft;
	DiffBoarderMaskSft[3]=DiffIon3BoarderMaskSft;
	
	InitZero();
}
NodeIndexing::~NodeIndexing()
{
	if(NIndex!=NULL)
	{
		delete [] NIndex;
		NIndex=NULL;
	}
	if(Q!=NULL)
	{
		delete [] Q;
		Q=NULL;
	}
}
int NodeIndexing::InitZero()
{
	unsigned int i;
	for(i=0;i<NodeIndexMaxValues;i++)
	{
		Eps[i]=0.0f;
		C[i]=0.0f;
		D[i]=0.0f;
	}
	PBC[0]=false;
	PBC[1]=false;
	PBC[2]=false;
	GridScale=0.0f;
	GridSize[0]=0;
	GridSize[1]=0;
	GridSize[2]=0;
	if(NIndex!=NULL)
	{
		delete [] NIndex;
		NIndex=NULL;
	}
	if(Q!=NULL)
	{
		delete [] Q;
		Q=NULL;
	}
	//QFormat=QFormatStraightSequence;
	//ChargeNum[0]=0;
	//ChargeNum[1]=0;
	//ChargeNum[2]=0;
	QNum=0;
	return EXIT_SUCCESS;
}
int NodeIndexing::SetNoPBC()
{
	//check if already NoPBC
	if(PBC[0]==false&&PBC[1]==false&&PBC[2]==false)return EXIT_SUCCESS;
	
	int i,j,k,gridPoint1,gridPoint2;
	int GridSizeOld[3]={GridSize[0],GridSize[1],GridSize[2]};
	int GridSizeNew[3]={GridSize[0],GridSize[1],GridSize[2]};
	if(PBC[0])GridSizeNew[0]-=2;
	if(PBC[1])GridSizeNew[1]-=2;
	if(PBC[2])GridSizeNew[2]-=2;
	int GridSizeNewXYZ=GridSizeNew[0]*GridSizeNew[1]*GridSizeNew[2];
	//Move data
	int itmp1=GridSizeNew[0];
	int itmp2=itmp1*GridSizeNew[1];
	int itmp3=GridSizeOld[0];
	int itmp4=itmp3*GridSizeOld[1];
	int itmp6,itmp5;
	
	
	//Nindex
	NodeIndex* NIndexNew=new NodeIndex[GridSizeNewXYZ];
	for(k=0;k<GridSizeNew[2];k++)
		for(j=0;j<GridSizeNew[1];j++)
			for(i=0;i<GridSizeNew[0];i++)
	{
		gridPoint1 = i + itmp1 * j + itmp2 *k;
		gridPoint2 = i + itmp3 * j + itmp4 *k;
		if(PBC[0])gridPoint2+=1;
		if(PBC[1])gridPoint2+=itmp3;
		if(PBC[2])gridPoint2+=itmp4;
		NIndexNew[gridPoint1]=NIndex[gridPoint2];
	}
	delete [] NIndex;
	NIndex=NIndexNew;
	//Q
	//!@todo here i belive that there is no charge on boarde, thus do not do anything
	//BlackNwhite
	//should work
	//set new values
	PBC[0]=false;
	PBC[1]=false;
	PBC[2]=false;
	GridSize[0]=GridSizeNew[0];
	GridSize[1]=GridSizeNew[1];
	GridSize[2]=GridSizeNew[2];
	
	return EXIT_SUCCESS;
}
int NodeIndexing::SetPBC(bool pbcX,bool pbcY,bool pbcZ)
{
	if(PBC[0]==pbcX&&PBC[1]==pbcY&&PBC[2]==pbcZ)return EXIT_SUCCESS;
	if(PBC[0]!=false||PBC[1]!=false||PBC[2]!=false)
	{
		SetNoPBC();
	}
	int GridSizeOld[3]={GridSize[0],GridSize[1],GridSize[2]};
	int GridSizeNew[3]={GridSize[0],GridSize[1],GridSize[2]};
	PBC[0]=pbcX;
	PBC[1]=pbcY;
	PBC[2]=pbcZ;
	if(PBC[0])GridSizeNew[0]+=2;
	if(PBC[1])GridSizeNew[1]+=2;
	if(PBC[2])GridSizeNew[2]+=2;
	
	//Move data
	int i,j,k;
	int gridPointOld,gridPointNew;
	int gridPoint1,gridPoint2;
	int GS_X_Old=GridSizeOld[0];
	int GS_XY_Old=GridSizeOld[0]*GridSizeOld[1];
	int GS_XYZ_Old=GridSizeOld[0]*GridSizeOld[1]*GridSizeOld[2];
	int GS_X_New=GridSizeNew[0];
	int GS_Y_New=GridSizeNew[1];
	int GS_Z_New=GridSizeNew[2];
	int GS_XY_New=GridSizeNew[0]*GridSizeNew[1];
	int GS_XYZ_New=GridSizeNew[0]*GridSizeNew[1]*GridSizeNew[2];
	int BnW=(NIndex[0]&BlackAndWhiteMask)>>BlackAndWhiteMaskSft;
	
		//set new values
	GridSize[0]=GridSizeNew[0];
	GridSize[1]=GridSizeNew[1];
	GridSize[2]=GridSizeNew[2];
	//Nindex
	NodeIndex* NIndexNew=new NodeIndex[GS_XYZ_New];
	DbgPrint0("NIndexNew %p %p\n",NIndexNew,NIndex);
	for(i=0;i<GS_XYZ_Old;i++)
		NIndexNew[i]=NIndex[i];
	
	delete [] NIndex;
	NIndex=NIndexNew;
	DbgPrint0("NIndexNew %p %p\n",NIndexNew,NIndex);
	ConvertToPBC(NIndex,GridSize[0]-2*pbcX,GridSize[1]-2*pbcY,GridSize[2]-2*pbcZ, pbcX, pbcY, pbcZ);
	
	/*
	//copy body
	for(k=GridSizeOld[2]-1;k>-1;k--)
		for(j=GridSizeOld[1]-1;j>-1;j--)
			for(i=GridSizeOld[0]-1;i>-1;i--)
	{
		gridPointOld = i + GS_X_Old * j + GS_XY_Old *k;
		gridPointNew = i + GS_X_New * j + GS_XY_New *k;
		if(PBC[0])gridPointNew+=1;
		if(PBC[1])gridPointNew+=GS_X_New;
		if(PBC[2])gridPointNew+=GS_XY_New;
		NIndexNew[gridPointNew]=NIndex[gridPointOld];
	}
	//copy sides;
	if(PBC[0])
		for(k=0;k<GridSizeNew[2];k++)
			for(j=0;j<GridSizeNew[1];j++)
	{
		gridPoint1=GS_X_New * j + GS_XY_New * k;
		gridPoint2=gridPoint1 + GS_X_New - 1;
		
		NIndexNew[gridPoint1]=NIndexNew[gridPoint2-1];
		NIndexNew[gridPoint2]=NIndexNew[gridPoint1+1];
	}
	if(PBC[1])
		for(k=0;k<GridSizeNew[2];k++)
			for(i=0;i<GridSizeNew[0];i++)
	{
		gridPoint1=i + GS_XY_New * k;
		gridPoint2=gridPoint1 + GS_X_New*(GS_Y_New - 1);
		
		NIndexNew[gridPoint1]=NIndexNew[gridPoint2-GS_X_New];
		NIndexNew[gridPoint2]=NIndexNew[gridPoint1+GS_X_New];
	}
	if(PBC[2])
		for(i=0;i<GridSizeNew[0];i++)
			for(j=0;j<GridSizeNew[1];j++)
	{
		gridPoint1=i + GS_X_New * j;
		gridPoint2=gridPoint1 + GS_XY_New*(GS_Z_New - 1);
		
		NIndexNew[gridPoint1]=NIndexNew[gridPoint2-GS_XY_New];
		NIndexNew[gridPoint2]=NIndexNew[gridPoint1+GS_XY_New];
	}
	*/
	

	
	CalcDielBoarder();
	CalcDiffBoarder();
	//BlackNwhite
	if(PBC[0])BnW++;
	if(PBC[1])BnW++;
	if(PBC[2])BnW++;
	BnW=BnW%2;
	SetBlackAndWhiteNodes(BnW);
	
	//Q
	//!@todo here i belive that there is no charge on boarde, thus do not do anything
	

	return EXIT_SUCCESS;
}
float NodeIndexing::GetDielConstInUse(int i)
{
	return Eps[i]*EPKT;
}
void NodeIndexing::SetDielConstInUse(int i, float Val)
{
	Eps[i]=Val/EPKT;
}
int NodeIndexing::SetNNodes(unsigned int *gridsize,float gridscale)
{
	unsigned int i;
	
	if(NIndex!=NULL)
	{
		delete NIndex;
		NIndex=NULL;
	}
	for(i=0;i<3;i++)GridSize[i]=gridsize[i];
	GridScale=gridscale;
	unsigned int NNodes=GridSize[0]*GridSize[1]*GridSize[2];
	NIndex=new NodeIndex[NNodes];
	
	for(i=0;i<NNodes;i++)
	{
		NIndex[i]=0;
	}
	SetBlackAndWhiteNodes();
	return EXIT_SUCCESS;
}
int NodeIndexing::SetIndexFieldFromIntArray(int *arr,NodeIndexDescriptor mask,NodeIndexDescriptor sft)
{
	if(arr==NULL)return EXIT_FAILURE;
	
	unsigned int i;
	NodeIndex Cur;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	unsigned int ClearMask=~mask;
	//clean field
	for(i=0;i<GS_XYZ;i++) NIndex[i]=NIndex[i]&ClearMask;
	
	for(i=0;i<GS_XYZ;i++)
	{
		//!@todo somehow set NIndex[i]&~mask+arr[i]<<sft
		Cur=arr[i]<<sft;
		NIndex[i]+=Cur;
	}
	
	return EXIT_SUCCESS;
}
int* NodeIndexing::GetIntArrayFromIndexField(NodeIndexDescriptor FieldType, NodeIndex Mask)
{
	if(FieldType<DielConst&&FieldType>=Charge)
	{
		fprintf(stderr,"Cannot operate with such FieldType\n");
		return NULL;
	}
	int i;
	int GridSizeXYZ= GridSize[0]*GridSize[1]*GridSize[2];
	NodeIndexDescriptor Sft=GetShtFromMask(Mask);

	int *Map=new int[GridSizeXYZ];
	for(i=0;i<GridSizeXYZ;i++)
	{
		Map[i]=(NIndex[i]&Mask)>>Sft;
	}
	return Map;
}
int NodeIndexing::SetIonAccess(float **DiffusionsMaps)
{
	if(DiffusionsMaps==NULL)return EXIT_FAILURE;
	
	int ion;
	unsigned int i;
	NodeIndex Cur;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	NodeIndex ClearMask,mask,sft;
	
	for(i=0;i<NodeIndexMaxValues;i++)
		D[i]=0.0;
	
	for(ion=0;ion<NIonsTypes;ion++)
	{
		if(ion==0)
		{
			mask=Ion0;
			sft=Ion0Sft;
		}
		if(ion==1)
		{
			mask=Ion1;
			sft=Ion1Sft;
		}
		if(ion==2)
		{
			mask=Ion2;
			sft=Ion2Sft;
		}
		if(ion==3)
		{
			mask=Ion3;
			sft=Ion3Sft;
		}
		ClearMask=~mask;
		//clean field
		for(i=0;i<GS_XYZ;i++) NIndex[i]=NIndex[i]&ClearMask;
		
		//find nonzero value for diffusion
		float vDiffusion=0.0;
		i=0;
		while(vDiffusion==0.0&&i<GS_XYZ)
		{
			vDiffusion=DiffusionsMaps[ion][i];
			i++;
		}
		D[ion+1]=vDiffusion;
		
		//fill ion fild with non zero values
		for(i=0;i<GS_XYZ;i++)
		{
			//!@todo somehow set NIndex[i]&~mask+arr[i]<<sft
			if(DiffusionsMaps[ion][i]>0.0)
			{
				Cur=(ion+1)<<sft;
				NIndex[i]+=Cur;
			}
		}
	}
	CalcDiffBoarder();
	return EXIT_SUCCESS;
}
NodeIndexing::NodeIndexDescriptor NodeIndexing::GetShtFromMask(NodeIndex Mask) const
{
	switch(Mask)
	{
		case Epsilon0:
			return Epsilon0Sft;
		case Epsilon1:
			return Epsilon1Sft;
		case Epsilon2:
			return Epsilon2Sft;
		case Ion0:
			return Ion0Sft;
		case Ion1:
			return Ion1Sft;
		case Ion2:
			return Ion2Sft;
		case Ion3:
			return Ion3Sft;
		case BlackAndWhiteMask:
			return BlackAndWhiteMaskSft;
		case ChargeMask:
			return ChargeMaskSft;
		case DielBoarderMask:
			return DielBoarderMaskSft;
		case DiffIon0BoarderMask:
			return DiffIon1BoarderMaskSft;
		case DiffIon1BoarderMask:
			return DiffIon2BoarderMaskSft;
		case DiffIon2BoarderMask:
			return DiffIon0BoarderMaskSft;
		case DiffIon3BoarderMask:
			return DiffIon3BoarderMaskSft;
			
		default:
			fprintf(stderr,"Cannot find shift for such mask, set to 0\n");
			return (NodeIndexDescriptor)0;
	}
}
float* NodeIndexing::GetCMap(NodeIndexDescriptor FieldType, NodeIndex Mask)
{
	if(FieldType<DielConst&&FieldType>Charge)
	{
		fprintf(stderr,"Cannot operate with such FieldType\n");
		return NULL;
	}
	int i;
	int GridSizeXYZ= GridSize[0]*GridSize[1]*GridSize[2];
	
	float *Vmap;
	NodeIndexDescriptor Sft=GetShtFromMask(Mask);
	if(FieldType==DielConst)Vmap=Eps;
	if(FieldType==DiffConst)Vmap=D;
	if(FieldType==Conc)Vmap=C;
	if(FieldType!=Charge)
	{
		float *Map=new float[GridSizeXYZ];
		for(i=0;i<GridSizeXYZ;i++)
		{
			Map[i]=Vmap[(NIndex[i]&Mask)>>Sft];
		}
		return Map;
	}
	else
	{
		return GetChargeArray();
	}
}
#ifdef HARLEM_MOD
HaField3D* NodeIndexing::GetHaField3D(NodeIndexDescriptor FieldType, NodeIndexDescriptor Mask)
{
	if(FieldType<DielConst&&FieldType>Charge)
	{
		fprintf(stderr,"Cannot operate with such FieldType\n");
		return NULL;
	}
	int i;
	int GridSizeXYZ= GridSize[0]*GridSize[1]*GridSize[2];
	//float fpoh= 4*M_PI*GridScale;
	//float coef=fpoh*COANGS/(GridScale*GridScale*GridScale);
	
	float Vmap[NodeIndexMaxValues];
	NodeIndexDescriptor Sft=GetShtFromMask(Mask);
	if(FieldType==DielConst)
	{
		for(i=0;i<NodeIndexMaxValues;i++)
		{
			Vmap[i]=EPKT*Eps[i];
		}
	}
	if(FieldType==DiffConst)
	{
		for(i=0;i<NodeIndexMaxValues;i++)
		{
			Vmap[i]=D[i];
		}
	}
	if(FieldType==Conc)
	{
		for(i=0;i<NodeIndexMaxValues;i++)
		{
			Vmap[i]=C[i];
		}
	}
	
	HaField3D *Field3D=new HaField3D();
	Field3D->SetDimensions(GridSize[0],GridSize[1],GridSize[2]);
	Field3D->SetCenterAsZero(GridScale);
	switch(Mask)
	{
		case Epsilon0:
			Field3D->SetName("Epsilon0");
			Field3D->ShiftGridCorners((double)0.5/GridScale,0,0);
			break;
		case Epsilon1:
			Field3D->SetName("Epsilon1");
			Field3D->ShiftGridCorners(0,(double)0.5/GridScale,0);
			break;
		case Epsilon2:
			Field3D->SetName("Epsilon2");
			Field3D->ShiftGridCorners(0,0,(double)0.5/GridScale);
			break;
		case Ion0:
			if(FieldType==DiffConst)Field3D->SetName("DiffConstIon0");
			else Field3D->SetName("ConcIon0");
			break;
		case Ion1:
			if(FieldType==DiffConst)Field3D->SetName("DiffConstIon1");
			else Field3D->SetName("ConcIon1");
			break;
		case Ion2:
			if(FieldType==DiffConst)Field3D->SetName("DiffConstIon2");
			else Field3D->SetName("ConcIon2");
			break;
		case Ion3:
			if(FieldType==DiffConst)Field3D->SetName("DiffConstIon3");
			else Field3D->SetName("ConcIon3");
			break;
		case ChargeMask:
			Field3D->SetName("ChargeStatic");
			break;
	}
	float *Map=Field3D->GetValPtr(0,0,0);
	if(FieldType!=Charge)
	{
		
		for(i=0;i<GridSizeXYZ;i++)
		{
			Map[i]=Vmap[(NIndex[i]&Mask)>>Sft];
		}
		return Field3D;
	}
	else
	{
		unsigned int count=0;
		float fpoh = 4*M_PI*GridScale;
		for(i=0;i<GridSizeXYZ;i++)
		{
			if(NIndex[i]&ChargeMask)
			{
				Map[i]=Q[count]/fpoh;
				count++;
			}
			else
			{
				Map[i]=0.0f;
			}
		}
		return Field3D;
	}
}
#endif

int NodeIndexing::ReadFromFile4SDMPI(const char *filename,ContWorld *World)
{
	NodeIndexing* NIndexing4Read=new NodeIndexing();
	NIndexing4Read->ReadFromFileAddPBC(filename,PBC[0],PBC[1],PBC[2]);
	this->GetMyPart4MPI(World,NIndexing4Read);
	delete NIndexing4Read;

	return EXIT_SUCCESS;
}
int NodeIndexing::ReadFromFile(const char *filename)
{
	gzFile file;
	char str[PNP_MAP_IO_STRING_LENGTH];
	int i;
	int GS[3];
	
	file = gzopen(filename,"rb");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	//Header
	if(Z_NULL==gzgets(file,str,PNP_MAP_IO_STRING_LENGTH)) {
		fprintf(stderr,"readMap: problem reading from file: %s\n", filename);
		return EXIT_FAILURE;
	}
	//string StrHeader(string);
	//istringstream ins(str);
	// @todo check ReadFromFile
	TiXmlElement *header=new TiXmlElement("NodeIndexing");
	header->Parse(str, NULL, TIXML_ENCODING_LEGACY);
	//ins>>*header;
	
	header->GetArrOfIntAttribute("GridSize",GS,3);
	
	GridSize[0]=GS[0];
	GridSize[1]=GS[1];
	GridSize[2]=GS[2];
	
	header->GetFloatAttribute("GridScale",&GridScale);

	if(header->GetIntAttribute("NIonsTypes",&NIonsTypes)==EXIT_FAILURE)
	{
		NIonsTypes=2;
	}
	IonsQ=new float[NIonsTypes];
	if(header->GetArrOfFloatAttribute("IonsQ",IonsQ,NIonsTypes)==EXIT_FAILURE)
	{
		NIonsTypes=2;
		delete [] IonsQ;
		IonsQ=new float[NIonsTypes];
		IonsQ[0]=1.0;
		IonsQ[1]=-1.0;
	}

	float Vtmp[NodeIndexMaxValues];
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS/(GridScale*GridScale*GridScale);
	
	
	header->GetArrOfFloatAttribute("Eps",Vtmp,NodeIndexMaxValues);
	for(i=0;i<NodeIndexMaxValues;i++)Eps[i]=Vtmp[i]/EPKT;
	
	header->GetArrOfFloatAttribute("D",D,NodeIndexMaxValues);
	
	header->GetArrOfFloatAttribute("C",Vtmp,NodeIndexMaxValues);
	for(i=0;i<NodeIndexMaxValues;i++)C[i]=Vtmp[i]*coef;
	
	//HaXML::SetAtribute(header,"TableType","TwoColumns");
	//int QNum;
	header->GetIntAttribute("QNum",&QNum);
	int BnW;
	header->GetIntAttribute("BnW",&BnW);
	delete header;
	
	unsigned int NNodes=GridSize[0]*GridSize[1]*GridSize[2];
	if(NIndex!=NULL)delete [] NIndex;
	NIndex=new NodeIndex[NNodes];
	ReadIndexGZfromTwoColumns(file,NIndex,NNodes);
	SetBlackAndWhiteNodes(BnW);
	Q=new float[QNum];
	double q=0.0;
	ReadMapGZOneColumns(file,Q,QNum);
	for(i=0;i<QNum;i++)
	{
		
		Q[i]*=fpoh;
		q+=Q[i];
		
	}
	DbgPrint0("q=%f QNum=%d\n",(float)q/fpoh,QNum);
	gzclose(file);
	
	DbgPrint0("Recalc. Boarders\n");
	CalcDielBoarder();
	CalcDiffBoarder();
	
	return EXIT_SUCCESS;
}
int NodeIndexing::ConvertToPBC(NodeIndex *V,int GS_X,int GS_Y,int GS_Z,bool pbcX,bool pbcY,bool pbcZ)
{
	if(pbcX==false&&pbcY==false&&pbcZ==false)
		return EXIT_SUCCESS;
	int GS[3]={GS_X, GS_Y, GS_Z};
	int GSwPBC[3]={GS_X+pbcX*2, GS_Y+pbcY*2, GS_Z+pbcZ*2};
	int i,j,k,gridPoint1,gridPoint2;
	int itmp1=GS[0];
	int itmp2=itmp1*GS[1];
	int itmp3=GSwPBC[0];
	int itmp4=itmp3*GSwPBC[1];
	int itmp6,itmp5;
	int GS_XYZ=GS_X*GS_Y*GS_Z;
	int GSwPBC_XYZ=GSwPBC[0]*GSwPBC[1]*GSwPBC[2];
	//If periodic Boundary condition then rearrange map
	if(GSwPBC_XYZ!=GS_XYZ)
	{
		//Move data
		for(k=GS[2]-1;k>-1;k--)
			for(j=GS[1]-1;j>-1;j--)
				for(i=GS[0]-1;i>-1;i--)
		{
			gridPoint1 = i + itmp1 * j + itmp2 *k;
			gridPoint2 = i + itmp3 * j + itmp4 *k;
			if(pbcX)gridPoint2+=1;
			if(pbcY)gridPoint2+=itmp3;
			if(pbcZ)gridPoint2+=itmp4;
			V[gridPoint2]=V[gridPoint1];
		}
		//copy sides;
		if(pbcX)
			for(k=0;k<GSwPBC[2];k++)
				for(j=0;j<GSwPBC[1];j++)
		{
			itmp1=itmp3 * j + itmp4 * k;
			itmp2=itmp1 + itmp3 - 2;
			V[itmp1]=V[itmp2];
			V[itmp2+1]=V[itmp1+1];
		}
		if(pbcY)
			for(k=0;k<GSwPBC[2];k++)
				for(i=0;i<GSwPBC[0];i++)
		{
			itmp1=i + itmp4 *k;
			itmp2=itmp1+itmp4 - itmp3;
			V[itmp1]=V[itmp2 - itmp3];
			V[itmp2]=V[itmp1 + itmp3];
		}
		if(pbcZ)
		{
			itmp1=itmp4*(GSwPBC[2]-2);
			itmp2=itmp1+itmp4;
			for(i=0;i<itmp4;i++)
			{
				V[i]=V[i+itmp1];
				V[i+itmp2]=V[i+itmp4];
			}
		}
	}
	return EXIT_SUCCESS;
}
int NodeIndexing::ReadFromFileAddPBC(const char *filename,bool pbcX,bool pbcY,bool pbcZ)
{
	this->ReadFromFile(filename);
	this->SetPBC(pbcX, pbcY, pbcZ);
	return EXIT_SUCCESS;
}

int NodeIndexing::WriteToFile(const char *filename)
{
	DbgPrint0("NodeIndexing::WriteIndexFile(%s)\n",filename);
	bool _PBC[3]={PBC[0],PBC[1],PBC[2]};
	SetNoPBC();
	//Prepare Values
	unsigned int i;
	float Vtmp[NodeIndexMaxValues];
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS/(GridScale*GridScale*GridScale);
	unsigned int GSXYZ=GridSize[0]*GridSize[1]*GridSize[2];
	int BnW;
	//Build XML Header
	TiXmlElement* header=new TiXmlElement("NodeIndexing");
	//int SetArrOfIntAttribute(const char* name,int* vals,int n);
	int GS[3]={GridSize[0],GridSize[1],GridSize[2]};
	header->SetArrOfIntAttribute("GridSize",GS,3);
	header->SetFloatAttribute("GridScale",GridScale);

	header->SetIntAttribute("NIonsTypes",NIonsTypes);
	header->SetArrOfFloatAttribute("IonsQ",IonsQ,NIonsTypes);
	
	for(i=0;i<NodeIndexMaxValues;i++)Vtmp[i]=EPKT*Eps[i];
	header->SetArrOfFloatAttribute("Eps",Vtmp,NodeIndexMaxValues);
	header->SetArrOfFloatAttribute("D",D,NodeIndexMaxValues);
	for(i=0;i<NodeIndexMaxValues;i++)Vtmp[i]=C[i]/coef;
	header->SetArrOfFloatAttribute("C",Vtmp,NodeIndexMaxValues);
	//int QNum=ChargeNum[2];
	header->SetIntAttribute("QNum",QNum);
	if(GSXYZ>0)BnW=(NIndex[0]&BlackAndWhiteMask)>>BlackAndWhiteMaskSft;
	else BnW=0;
	header->SetIntAttribute("BnW",BnW);
	header->SetAttribute("Comment","Epsilon,Diel,Flags in TwoColumns; Charges in 1");
	
	gzFile file;
	file = gzopen(filename,"wb");
	if(file==NULL) {
		fprintf(stderr,"ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	string StrHeader;
	StrHeader<<*header;
	gzprintf(file,"%s\n",StrHeader.c_str());
	//Set BlackAndWhiteMask flag to 0 for bigger compression
	int notBnW=~BlackAndWhiteMask;
	for(i=0;i<GSXYZ;i++)
	{
		NIndex[i]=NIndex[i]&notBnW;
	}
	WriteIndexGZinTwoColumns(file,NIndex,GSXYZ);
	//restore BlackAndWhiteMask flag
	SetBlackAndWhiteNodes(BnW);
	for(i=0;i<QNum;i++)
	{
		Q[i]/=fpoh;
	}
	WriteMapGZOneColumns(file,Q,QNum);
	for(i=0;i<QNum;i++)
	{
		Q[i]*=fpoh;
	}
	gzclose(file);
	delete header;
	SetPBC(_PBC[0],_PBC[1],_PBC[2]);
	return EXIT_SUCCESS;
}
int NodeIndexing::GetCentralPartFromNIExtMyPart4MPI(NodeIndexing* NIExt)
{
	int i;
	this->NIonsTypes=NIExt->NIonsTypes;
	for(i=0;i<NodeIndexMaxValues;i++)
	{
		this->Eps[i]=NIExt->Eps[i];
		this->D[i]=NIExt->D[i];
		this->C[i]=NIExt->C[i];
	}
	
	int CentralNode[3];
	int CentralNodeExt[3];
	int ixExt,iyExt,izExt,GrdPntExt;
	
	for(i=0;i<3;i++)
	{
		CentralNode[i]=GridSize[i]/2;
		CentralNodeExt[i]=NIExt->GridSize[i]/2;
	}

	int ix,iy,iz,GrdPnt;
	for(ix=0;ix<GridSize[0];ix++)
		for(iy=0;iy<GridSize[1];iy++)
			for(iz=0;iz<GridSize[2];iz++)
	{
		GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
		ixExt=ix-CentralNode[0]+CentralNodeExt[0];
		iyExt=iy-CentralNode[1]+CentralNodeExt[1];
		izExt=iz-CentralNode[2]+CentralNodeExt[2];
		GrdPntExt=ixExt+iyExt*NIExt->GridSize[0]+izExt*NIExt->GridSize[0]*NIExt->GridSize[1];
		this->NIndex[GrdPnt]=NIExt->NIndex[GrdPntExt];
	}
	
	//startBlackAndWhite=this->NIndex[0]&NodeIndexing::BlackAndWhiteMask;
	
	//Transmit Q
	QNum=0;
	for(GrdPnt=0;GrdPnt<GridSize[0]*GridSize[1]*GridSize[2];GrdPnt++)
	{
		if(NIndex[GrdPnt]&ChargeMask)
		{
			QNum++;
		}
	}
	DeleteCArray(Q);
	if(QNum>0)
		this->Q=new float[QNum];
	int QExtCount=0,Qcount=0;
	for(ixExt=0;ixExt<NIExt->GridSize[0];ixExt++)
		for(iyExt=0;iyExt<NIExt->GridSize[1];iyExt++)
			for(izExt=0;izExt<NIExt->GridSize[2];izExt++)
	{
		ix=ixExt+CentralNode[0]-CentralNodeExt[0];
		iy=iyExt+CentralNode[1]-CentralNodeExt[1];
		iz=izExt+CentralNode[2]-CentralNodeExt[2];
		GrdPntExt=ixExt+iyExt*NIExt->GridSize[0]+izExt*NIExt->GridSize[0]*NIExt->GridSize[1];
		GrdPnt=ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1];
		
		
		if(NIExt->NIndex[GrdPntExt]&ChargeMask)
		{
			if(ix>=0&&ix<GridSize[0]&&iy>=0&&iy<GridSize[1]&&iz>=0&&iz<GridSize[2])
			{
				Q[Qcount]=NIExt->Q[QExtCount];
				Qcount++;
			}
			QExtCount++;
		}
	}
	this->CheckBoarder4Q();
	return EXIT_SUCCESS;
}
int NodeIndexing::GetMyPart4MPI(ContWorld *World,NodeIndexing* NIndexGlobal)
{
#ifdef MPI_PARALLEL
	int i,j,k;
	int MyRank = pnpsapp->MyComGroup.Get_rank();
	int NProcs = pnpsapp->MyComGroup.Get_size();
	pnpsapp->MyComGroup.Barrier();
	if(MyRank==0)
	{
		//NIndexGlobal->WriteToFile("Global.gz");
		unsigned int gridsize[3]={World->GridSize[0],World->GridSize[1],World->GridSize[2]};
		this->SetNNodes(gridsize,World->GridScale);
					
		pnpsapp->MyComGroup.Bcast(&(NIndexGlobal->NIonsTypes), 1, MPI::INT, 0);
		pnpsapp->MyComGroup.Bcast(NIndexGlobal->Eps, NodeIndexMaxValues, MPI::FLOAT, 0);
		pnpsapp->MyComGroup.Bcast(NIndexGlobal->D, NodeIndexMaxValues, MPI::FLOAT, 0);
		pnpsapp->MyComGroup.Bcast(NIndexGlobal->C, NodeIndexMaxValues, MPI::FLOAT, 0);
					
					
		int dest=0;
		unsigned int specChargeMask=NodeIndexing::ChargeMask;
		unsigned int ChargeDielBoarderMask=specChargeMask|NodeIndexing::DielBoarderMask;
		unsigned int BlackAndWhiteMask=NodeIndexing::BlackAndWhiteMask;
		unsigned int DielBoarderMask=NodeIndexing::DielBoarderMask;
		int G_XY=NIndexGlobal->GridSize[0]*NIndexGlobal->GridSize[1];
		
		int sz1,sz2,StartArr,EndArr,SizeArr;
		World->GetBorder(&sz1,&sz2,dest);
		StartArr=sz1*GridSize[0]*GridSize[1];
		EndArr=(sz2+1)*GridSize[0]*GridSize[1];
		SizeArr=EndArr-StartArr;
		DbgPrint0("GetMyPart4MPI: SizeArr[%d]=%d z0=%d z1=%d BnW=%d\n",dest,SizeArr,sz1,sz2,int(NIndexGlobal->NIndex[sz1*G_XY]&BlackAndWhiteMask));
		for(dest=1;dest<NProcs;dest++)
		{
			//NIndex
			
			World->GetBorder(&sz1,&sz2,dest);
			StartArr=sz1*GridSize[0]*GridSize[1];
			EndArr=(sz2+1)*GridSize[0]*GridSize[1];
			SizeArr=EndArr-StartArr;
			DbgPrint0("GetMyPart4MPI: SizeArr[%d]=%d z0=%d z1=%d BnW=%d\n",dest,SizeArr,sz1,sz2,int(NIndexGlobal->NIndex[sz1*G_XY]&BlackAndWhiteMask));
			pnpsapp->MyComGroup.Send((int*)(NIndexGlobal->NIndex)+StartArr, SizeArr, MPI::INT, dest, 0);
		}
		//Q
		int *StartQArr=new int[NProcs];
		int *EndQArr=new int[NProcs];
		int *SizeQArr=new int[NProcs];
		//!@todo now we exclude only XY charges but it shouldnt metter if initial NIndex was ok
		int *StartLocSys=new int[NProcs];
		int *EndLocSys=new int[NProcs];
		DbgPrint0("StartLocSys=%p\n",StartLocSys);
		for(dest=0;dest<NProcs;dest++)
		{
			World->GetBorder(&i,&j,dest);
			StartLocSys[dest]=(i+1)*G_XY;
			EndLocSys[dest]=j*G_XY-1;
			StartQArr[dest]=0;
			EndQArr[dest]=0;
			SizeQArr[dest]=0;
			//DbgPrint0("Border[%d]=%d %d, %d %d, %d\n",dest,i,j,StartLocSys[dest],EndLocSys[dest],G_XY);
		}
		int totGS_XYZ = NIndexGlobal->GridSize[0] * NIndexGlobal->GridSize[1] * NIndexGlobal->GridSize[2];
		int count;
					
		
		
		for(dest=0;dest<NProcs;dest++)
		{
			count=0;
			for(i=StartLocSys[dest];i<=EndLocSys[dest];i++)
			{
				if(NIndexGlobal->NIndex[i]&NIndexGlobal->ChargeMask)
				{
					count++;
				}
			}
			SizeQArr[dest]=count;
		}
		dest=0;
		count=0;
		for(i=StartLocSys[0];i<=EndLocSys[NProcs-1];i++)
		{
			if(i==EndLocSys[dest])
			{
				EndQArr[dest]=count-1;
				dest++;
			}
			if(NIndexGlobal->NIndex[i]&NIndexGlobal->ChargeMask)
			{
				count++;
			}
		}
		StartQArr[0]=0;
		for(dest=1;dest<NProcs;dest++)
		{
			if(SizeQArr[dest]>0)
				StartQArr[dest]=EndQArr[dest-1]+1;
			else
				StartQArr[dest]=EndQArr[dest];
		}
		/*for(dest=0;dest<NProcs;dest++)
		{
			count=0;
			//!@todo optimize this cycle
			for(i=StartLocSys[0];i<=EndLocSys[dest];i++)
			{
				if(i==StartLocSys[dest])
					StartQArr[dest]=count;
				if(i==EndLocSys[dest])
					EndQArr[dest]=count;
				if(NIndexGlobal->NIndex[i]&NIndexGlobal->ChargeMask)
				{
					count++;
				}
			}
		}*/
		for(dest=0;dest<NProcs;dest++)
		{
			World->GetBorder(&i,&j,dest);
			DbgPrint0("StartQArr[%d]=%d,%d ,%d LocSys[%d]=%d,%d, z=[%d,%d]\n",dest,StartQArr[dest],EndQArr[dest],SizeQArr[dest],dest,StartLocSys[dest],EndLocSys[dest],i,j);
		}
		pnpsapp->MyComGroup.Barrier();
		for(dest=1;dest<NProcs;dest++)
		{
			pnpsapp->MyComGroup.Send(&SizeQArr[dest], 1, MPI::INT, dest, 0);
			DbgPrint0("SizeArr[%d]=%d\n",dest,SizeQArr[dest]);
			if(SizeQArr[dest]>0)pnpsapp->MyComGroup.Send((NIndexGlobal->Q)+StartQArr[dest], SizeQArr[dest], MPI::FLOAT, dest, 0);
		}
					
		//now about me
		this->NIonsTypes=NIndexGlobal->NIonsTypes;
		for(i=0;i<NodeIndexMaxValues;i++)
		{
			this->Eps[i]=NIndexGlobal->Eps[i];
			this->D[i]=NIndexGlobal->D[i];
			this->C[i]=NIndexGlobal->C[i];
		}
		{
			//int sz1,sz2,StartArr,EndArr,SizeArr;
			World->GetBorder(&sz1,&sz2,0);
			StartArr=sz1*GridSize[0]*GridSize[1];
			EndArr=(sz2+1)*GridSize[0]*GridSize[1];
			SizeArr=EndArr-StartArr;
			DbgPrint0("StartArr %d %d\n",StartArr,EndArr);
			for(i=StartArr;i<EndArr;i++)
			{
				this->NIndex[i]=NIndexGlobal->NIndex[i];
			}
			this->QNum=SizeQArr[0];
			DbgPrint0("QNum %d %d\n",this->QNum,EndQArr[0]);
			if(this->QNum>0)
			{
				this->Q=new float[SizeArr];
				for(i=StartQArr[0];i<=EndQArr[0];i++)
				{
					this->Q[i]=NIndexGlobal->Q[i];
				}
			}
			else
			{
				this->Q=NULL;
			}
		}
					
		World->startBlackAndWhite=this->NIndex[0]&NodeIndexing::BlackAndWhiteMask;
		DbgPrint0("StartLocSys=%p\n",StartLocSys);
		delete [] StartLocSys;
		delete [] EndLocSys;
		delete [] StartQArr;
		delete [] EndQArr;
		delete [] SizeQArr;
		
		this->CheckBoarder4Q();
		
	}
	else
	{
		int iz1,iz2;
		World->GetBorder(&iz1,&iz2,MyRank);
		unsigned int gridsize[3]={World->GridSize[0],World->GridSize[1],World->GridSize[2]};
		this->SetNNodes(gridsize,World->GridScale);
				
		pnpsapp->MyComGroup.Bcast(&(this->NIonsTypes), 1, MPI::INT, 0);
		pnpsapp->MyComGroup.Bcast(this->Eps, NodeIndexMaxValues, MPI::FLOAT, 0);
		pnpsapp->MyComGroup.Bcast(this->D, NodeIndexMaxValues, MPI::FLOAT, 0);
		pnpsapp->MyComGroup.Bcast(this->C, NodeIndexMaxValues, MPI::FLOAT, 0);
				
		int SizeArr;
		SizeArr=(iz2-iz1+1)*GridSize[0]*GridSize[1];
		DbgPrint0("SizeArr=%d\n",SizeArr);
		pnpsapp->MyComGroup.Recv((int*)(this->NIndex), SizeArr, MPI::INT, 0, 0);
		World->startBlackAndWhite=this->NIndex[0]&NodeIndexing::BlackAndWhiteMask;
				
		pnpsapp->MyComGroup.Barrier();
		pnpsapp->MyComGroup.Recv(&SizeArr, 1, MPI::INT, 0, 0);
		this->QNum=SizeArr;
		DbgPrint0("QNum=%d\n",SizeArr);
		if(SizeArr>0)
		{
			this->Q=new float[SizeArr];
			pnpsapp->MyComGroup.Recv((this->Q), SizeArr, MPI::FLOAT, 0, 0);
		}
		this->CheckBoarder4Q();
	}
	DbgPrint0("GetMyPart4MPI: BnW=%d\n",int(this->NIndex[0]&NodeIndexing::BlackAndWhiteMask));
	pnpsapp->MyComGroup.Barrier();
#endif
	return EXIT_SUCCESS;
}
// int NodeIndexing::CollectGlobalSysFromMPIParts(CWorld *World,NodeIndexing* NIndexGlobal)
// {
// 	return EXIT_SUCCESS;
// }
int NodeIndexing::SetBlackAndWhiteNodes(int FirstNode)
{
	if(NIndex==NULL)
	{
		fprintf(stderr,"ERROR 102: int NodeIndexing::SetBlackAndWhiteNodes(int FirstNode): NIndex not init\n");
		return EXIT_FAILURE;
	}
	//temp vars
	int BlackOrWhite;
	unsigned int i,j,k,kgrid,jgrid;
	unsigned int GrdPnt;
	unsigned int GS_X = GridSize[0];
	unsigned int GS_Y = GridSize[1];
	unsigned int GS_Z = GridSize[2];
	unsigned int GS_XY=GS_X*GS_Y;
	unsigned int GS_XYZ=GS_XY*GS_Z;
	int notBnW=~BlackAndWhiteMask;
	//Black or white
	for(k=0;k<GS_Z;k++)
	{
		kgrid = k*GS_XY;
		for(j=0;j<GS_Y;j++)
		{
			jgrid = kgrid+j*GS_X;
			for(i=0;i<GS_X;i++)
			{
				GrdPnt = jgrid+i;
				BlackOrWhite=k+j+i+FirstNode;
				//!@todo it defenitly can be more effective
				NIndex[GrdPnt] = (NIndex[GrdPnt]&notBnW)|((BlackOrWhite%2)<<BlackAndWhiteMaskSft);
			}
		}
	}
	return EXIT_SUCCESS;
}
int NodeIndexing::CalcDielBoarder()
{
	DbgPrint1("NodeIndexing::CalcDielBoarder()\n");
	if(NIndex==NULL)
	{
		fprintf(stderr,"ERROR 102: int NodeIndexing::SetBlackAndWhiteNodes(int FirstNode): NIndex not init\n");
		return EXIT_FAILURE;
	}
	//temp vars
	int BlackOrWhite;
	unsigned int i,j,k,kgrid,jgrid;
	unsigned int GrdPnt;
	unsigned int GS_X = GridSize[0];
	unsigned int GS_Y = GridSize[1];
	unsigned int GS_Z = GridSize[2];
	unsigned int GS_XY=GS_X*GS_Y;
	unsigned int GS_XYZ=GS_XY*GS_Z;
	//clean dielboarder
	unsigned int UnDielBoarder=~DielBoarderMask;
	int bcount=0;
	
	for(i=0;i<GS_XYZ;i++) NIndex[i]=NIndex[i]&UnDielBoarder;
	for(k=1;k<GS_Z-1;k++)
	{
		kgrid = k*GS_XY;
		for(j=1;j<GS_Y-1;j++)
		{
			jgrid = kgrid+j*GS_X;
			for(i=1;i<GS_X-1;i++)
			{
				
				GrdPnt = jgrid+i;
//				int sum= ((NIndex[GrdPnt]&Epsilon0)>>Epsilon0Sft) + ((NIndex[GrdPnt]&Epsilon1)>>Epsilon1Sft) + ((NIndex[GrdPnt]&Epsilon2)>>Epsilon2Sft) + ((NIndex[GrdPnt-1]&Epsilon0)>>Epsilon0Sft) +  ((NIndex[GrdPnt-GS_X]&Epsilon1)>>Epsilon1Sft) + ((NIndex[GrdPnt-GS_XY]&Epsilon2)>>Epsilon2Sft);
				if(((NIndex[GrdPnt]&Epsilon0)>>Epsilon0Sft)!=((NIndex[GrdPnt-1]&Epsilon0)>>Epsilon0Sft) ||
					 ((NIndex[GrdPnt]&Epsilon0)>>Epsilon0Sft)!=((NIndex[GrdPnt]&Epsilon1)>>Epsilon1Sft)||
					 ((NIndex[GrdPnt]&Epsilon0)>>Epsilon0Sft)!=((NIndex[GrdPnt-GS_X]&Epsilon1)>>Epsilon1Sft) ||
					 ((NIndex[GrdPnt]&Epsilon0)>>Epsilon0Sft)!=((NIndex[GrdPnt]&Epsilon2)>>Epsilon2Sft) ||
					 ((NIndex[GrdPnt]&Epsilon0)>>Epsilon0Sft)!=((NIndex[GrdPnt-GS_XY]&Epsilon2)>>Epsilon2Sft))
				//if(sum!=6*((NIndex[GrdPnt]&Epsilon0)>>Epsilon0Sft))
				{
					NIndex[GrdPnt] = NIndex[GrdPnt]|DielBoarderMask;
					bcount++;
				}
			}
		}
	}
	DbgPrint0("CalcDielBoarder %d\n",bcount);
	return EXIT_SUCCESS;
}
int NodeIndexing::CalcDiffBoarder()
{
	DbgPrint1("NodeIndexing::CalcDiffBoarder()\n");
	if(NIndex==NULL)
	{
		fprintf(stderr,"ERROR 102: int NodeIndexing::CalcDiffBoarder(int FirstNode): NIndex not init\n");
		return EXIT_FAILURE;
	}
	//temp vars
	int BlackOrWhite;
	unsigned int i,j,k,kgrid,jgrid,ion;
	unsigned int GrdPnt;
	unsigned int GS_X = GridSize[0];
	unsigned int GS_Y = GridSize[1];
	unsigned int GS_Z = GridSize[2];
	unsigned int GS_XY=GS_X*GS_Y;
	unsigned int GS_XYZ=GS_XY*GS_Z;
	//clean dielboarder
	unsigned int UnDiffBoarder;
	unsigned int DiffBoarderMask;

	for(ion=0;ion<NIonsTypes;ion++)
	{
		switch(ion)
		{
			case 0:
				DiffBoarderMask=DiffIon0BoarderMask;
				break;
			case 1:
				DiffBoarderMask=DiffIon1BoarderMask;
				break;
			case 2:
				DiffBoarderMask=DiffIon2BoarderMask;
				break;
			case 3:
				DiffBoarderMask=DiffIon3BoarderMask;
				break;
		}
		UnDiffBoarder=~DiffBoarderMask;
		//unset diffboarder
		for(i=0;i<GS_XYZ;i++) NIndex[i]=NIndex[i]&UnDiffBoarder;
		
		//check for ziro's set all to C[0]/D[0]
// 		for(i=0;i<GS_XYZ;i++)
// 		{
// 			if(GetDiffFloat(ion,i)==0.0)SetIonField(ion,i,0);
// 		}
		//remove bad Diffuson nodes
		for(k=1;k<GS_Z-1;k++)
		{
				kgrid = k*GS_XY;
				for(j=1;j<GS_Y-1;j++)
				{
						jgrid = kgrid+j*GS_X;
						for(i=1;i<GS_X-1;i++)
						{
								GrdPnt = jgrid+i;
								if(GetIonField(ion,GrdPnt)!=0 && GetIonField(ion,GrdPnt+1)==0 && GetIonField(ion,GrdPnt-1)==0 && GetIonField(ion,GrdPnt+GS_X)==0 && GetIonField(ion,GrdPnt-GS_X)==0 && GetIonField(ion,GrdPnt+GS_XY)==0 &&GetIonField(ion,GrdPnt-GS_XY)==0)
												SetIonField(ion,GrdPnt,0);
						}
				}
		}
		for(k=1;k<GS_Z-1;k++)
		{
			kgrid = k*GS_XY;
			for(j=1;j<GS_Y-1;j++)
			{
				jgrid = kgrid+j*GS_X;
				for(i=1;i<GS_X-1;i++)
				{
					GrdPnt = jgrid+i;
					if(GetIonField(ion,GrdPnt)!=GetIonField(ion,GrdPnt+1) || GetIonField(ion,GrdPnt)!=GetIonField(ion,GrdPnt-1) || GetIonField(ion,GrdPnt)!=GetIonField(ion,GrdPnt+GS_X) || GetIonField(ion,GrdPnt)!=GetIonField(ion,GrdPnt-GS_X) || GetIonField(ion,GrdPnt)!=GetIonField(ion,GrdPnt+GS_XY) || GetIonField(ion,GrdPnt)!=GetIonField(ion,GrdPnt-GS_XY))
					{
							NIndex[GrdPnt] = NIndex[GrdPnt]|DiffBoarderMask;
					}
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
int NodeIndexing::SetSCharDielMapForGAPS(signed char **CDiel,float *gapsEps)
{
	int i,j;
	for(i=0;i<4;i++)
	{
		DbgPrint0("CDiel[%d]=%p\n",i,CDiel[i]);
	}
	for(i=0;i<NodeIndexMaxValues;i++)
	{
		DbgPrint0("Eps[%d]=\t%f=\t%f\n",i,Eps[i],gapsEps[i]);
	}
	unsigned int GrdPnt;
	unsigned int GS_X = GridSize[0];
	unsigned int GS_Y = GridSize[1];
	unsigned int GS_Z = GridSize[2];
	unsigned int GS_XY=GS_X*GS_Y;
	unsigned int GS_XYZ=GS_XY*GS_Z;
	
	int iEpsConv[NodeIndexMaxValues];
	int iGapsZero;
	for(j=0;j<NodeIndexMaxValues;j++)
	{
		if(gapsEps[j]==0.0)
		{
			iGapsZero=j;
			break;
		}
	}
	for(i=0;i<NodeIndexMaxValues;i++)
	{
		if(Eps[i]!=0.0)
		{
			for(j=0;j<NodeIndexMaxValues;j++)
			{
				if(fabs(Eps[i]-gapsEps[j])<0.0000001)
				{
					iEpsConv[i]=j;
					break;
				}
			}
		}
		else
		{
			iEpsConv[i]=iGapsZero;
		}
	}
	for(i=0;i<NodeIndexMaxValues;i++)
	{
		DbgPrint0("iEpsConv[%d]=%d\n",i,iEpsConv[i]);
	}
	
	CDiel[0]=new signed char[GS_XYZ];
	CDiel[1]=new signed char[GS_XYZ];
	CDiel[2]=new signed char[GS_XYZ];
	CDiel[3]=new signed char[GS_XYZ];
	
	for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
	{
		CDiel[0][GrdPnt]=iEpsConv[GetDiel(0,GrdPnt)];
		CDiel[1][GrdPnt]=iEpsConv[GetDiel(1,GrdPnt)];
		CDiel[2][GrdPnt]=iEpsConv[GetDiel(2,GrdPnt)];
		CDiel[3][GrdPnt]=0;
	}
	for(i=0;i<4;i++)
	{
		DbgPrint0("CDiel[%d]=%p\n",i,CDiel[i]);
	}
	return EXIT_SUCCESS;
}
int NodeIndexing::CheckNodeIndex()
{
	unsigned int i,j,k,kgrid,jgrid,ion;
	unsigned int GrdPnt;
	unsigned int GS_X = GridSize[0];
	unsigned int GS_Y = GridSize[1];
	unsigned int GS_Z = GridSize[2];
	unsigned int GS_XY=GS_X*GS_Y;
	unsigned int GS_XYZ=GS_XY*GS_Z;
	
	unsigned int specChargeMask=ChargeMask;
	unsigned int ChargeDielBoarderMask=specChargeMask|DielBoarderMask;
	unsigned int BlackAndWhiteMask=BlackAndWhiteMask;
	unsigned int DielBoarderMask=DielBoarderMask;
	
	int CntQst=0;
	int CntQmob=0;
	int CntQstDielBoarder=0;
	int CntQmobDielBoarder=0;
	int CntRest=0;
	int GS[3]={GridSize[0],GridSize[1],GridSize[2]};
	VectorIntField3D* Vi=new VectorIntField3D(GS,GridScale,1);
	for(k=1;k<GS_Z-1;k++)
		for(j=1;j<GS_Y-1;j++)
			for(i=1;i<GS_X-1;i++)
	{
		GrdPnt = i+j*GS_X+k*GS_XY;
		if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
		{
			Vi->V[0][GrdPnt]=1;
		}
		else
		{
			Vi->V[0][GrdPnt]=0;
		}
	}
	Vi->WriteToFile("qstboarder.gz");
	delete Vi;
	/*for(k=1;k<GS_Z-1;k++)
		for(j=1;j<GS_Y-1;j++)
			for(i=1;i<GS_X-1;i++)
	{
		GrdPnt = i+j*GS_X+k*GS_Y;
		if((NIndex[GrdPnt]&ChargeDielBoarderMask)==ChargeDielBoarderMask)
		{
			
		}
		else if(NIndex[GrdPnt]&specChargeMask)
		{
			
		}
		else if(NIndex[GrdPnt]&DielBoarderMask)
		{
			 
		}
		else
		{
			CntRest++;
		}
	}*/
	return EXIT_SUCCESS;
}
int NodeIndexing::RemoveBadDiffusionPoints()
{
	int i,j,k;
	int IType;
	int GrdPnt;
	int GS_X=GridSize[0];
	int GS_XY=GridSize[0]*GridSize[1];
	int count,countstot=0;
	int DiffZeroInd=GetDiffZeroInd();
	int NodesWithZeroD;
	int countZeroDiffAround[4]={0,0,0,0};
	int countOneNonZeroDiffAround[4]={0,0,0,0};
	int countCycle[4]={0,0,0,0};
	
	pnpPrint("<NodeIndexing::RemoveBadDiffusionPoints>\n");
	for(IType=0;IType<NIonsTypes;IType++)
	{
		countCycle[IType]=0;
		do
		{
			count=0;
			for(k=1;k<GridSize[2]-1;k++)
			{
				for(j=1;j<GridSize[1]-1;j++)
				{
					for(i=1;i<GridSize[0]-1;i++)
					{
						GrdPnt = i+j*GS_X+k*GS_XY;
						if(GetDiffFloat(IType,GrdPnt)!=0.0 && GetDiffFloat(IType,GrdPnt+1)==0.0 && GetDiffFloat(IType,GrdPnt-1)==0.0 && GetDiffFloat(IType,GrdPnt+GS_X)==0.0 && GetDiffFloat(IType,GrdPnt-GS_X)==0.0 && GetDiffFloat(IType,GrdPnt+GS_XY)==0.0 &&GetDiffFloat(IType,GrdPnt-GS_XY)==0.0)
						{
							SetIonField(IType,GrdPnt,DiffZeroInd);
							//DbgPrint2("Bad Diffusional point at %d [%d %d %d] for ion %d have removed it. Because all nodes around have zero diffusion\n",GrdPnt,i,j,k,IType);
							count++;
							countstot++;
							countZeroDiffAround[IType]++;
						}
						if(GetDiffFloat(IType,GrdPnt)!=0.0)
						{
							NodesWithZeroD=0;
							if(GetDiffFloat(IType,GrdPnt+1)==0.0)NodesWithZeroD++;
							if(GetDiffFloat(IType,GrdPnt-1)==0.0)NodesWithZeroD++;
							if(GetDiffFloat(IType,GrdPnt+GS_X)==0.0)NodesWithZeroD++;
							if(GetDiffFloat(IType,GrdPnt-GS_X)==0.0)NodesWithZeroD++;
							if(GetDiffFloat(IType,GrdPnt+GS_XY)==0.0)NodesWithZeroD++;
							if(GetDiffFloat(IType,GrdPnt-GS_XY)==0.0)NodesWithZeroD++;
							if(NodesWithZeroD==5)
							{
								SetIonField(IType,GrdPnt,DiffZeroInd);
								//DbgPrint2("Bad Diffusional point at %d [%d %d %d] for ion %d have removed it. Because only one node around hove not zero diffusion\n",GrdPnt,i,j,k,IType);
								count++;
								countstot++;
								countOneNonZeroDiffAround[IType]++;
							}
						}
					}
				}
			}
			countCycle[IType]++;
		}
		while(count>0);
	}
	if(countstot>0)
	{
		CalcDiffBoarder();
		
		pnpPrint("Totally %d points was removed\n",countstot);
		for(IType=0;IType<NIonsTypes;IType++)
		{
			pnpPrint("\tFor ion %d removed in %d cycles\n",IType,countCycle[IType]);
			pnpPrint("\t                    %d points where all points around have zero diffusion\n",countZeroDiffAround[IType]);
			pnpPrint("\t                    %d points where only one points has non zero diffusion\n",countOneNonZeroDiffAround[IType]);
		}
	}
	pnpPrint("</NodeIndexing::RemoveBadDiffusionPoints>\n");
	return EXIT_SUCCESS;
}
int NodeIndexing::SetChargeMapFromArrayNoQonBoarder(float *q)
{
	//temp vars
	unsigned int i,j,k,grd;
	unsigned int GS_X=GridSize[0];
	unsigned int GS_XY=GridSize[0]*GridSize[1];
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	
	unsigned int ClearMask=~ChargeMask;
	for(i=0;i<GS_XYZ;i++)NIndex[i] = NIndex[i]&ClearMask;
	
	QNum=0;
	for(k=1;k<GridSize[2]-1;k++)
		for(j=1;j<GridSize[1]-1;j++)
			for(i=1;i<GridSize[0]-1;i++)
	{
		grd=i+j*GS_X+k*GS_XY;
		if(q[grd]!=0.0)
		{
			QNum++;
		}
	}
	
	DbgPrint0("NodeIndexing::QNum=%d\n",QNum);
	if(Q!=NULL)delete [] Q;
	Q=new float[QNum];
	
	QNum=0;
	for(k=1;k<GridSize[2]-1;k++)
		for(j=1;j<GridSize[1]-1;j++)
			for(i=1;i<GridSize[0]-1;i++)
	{
		grd=i+j*GS_X+k*GS_XY;
		if(q[grd]!=0.0)
		{
			Q[QNum]=q[grd];
			NIndex[grd] = NIndex[grd]|ChargeMask;
			QNum++;
		}
	}
	//ChargeNum[2]=QNum;
	//QFormat=QFormatStraightSequence;
	return EXIT_SUCCESS;
}
int NodeIndexing::SetChargeMapToZero()
{
	int i;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	QNum=0;
	unsigned int ClearMask=~ChargeMask;
	for(i=0;i<GS_XYZ;i++)
	{
		NIndex[i] = NIndex[i]&ClearMask;
	}
	//ChargeNum[2]=QNum;
	return EXIT_SUCCESS;
}
int NodeIndexing::SetChargeMapFromArray(float *q)
{
	//temp vars
	unsigned int i;
	unsigned int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	QNum=0;
	for(i=0;i<GS_XYZ;i++)
	{
		if(q[i]!=0.0)QNum++;
	}
	DbgPrint0("NodeIndexing::QNum=%d\n",QNum);
	if(Q!=NULL)delete [] Q;
	Q=new float[QNum];
	QNum=0;
	unsigned int ClearMask=~ChargeMask;
	for(i=0;i<GS_XYZ;i++)
	{
		NIndex[i] = NIndex[i]&ClearMask;
		if(q[i]!=0.0)
		{
			Q[QNum]=q[i];
			NIndex[i] = NIndex[i]|ChargeMask;
			QNum++;
		}
	}
	//ChargeNum[2]=QNum;
	//QFormat=QFormatStraightSequence;
	return EXIT_SUCCESS;
}
int* NodeIndexing::GetChargeIndex()
{
	int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	int *qndx=new int[QNum];
	int i,count=0;
	
	for(i=0;i<GS_XYZ;i++)
	{
		if(NIndex[i]&ChargeMask)
		{
			qndx[count]=i;
			count++;
		}
	}
	return qndx;
}
float* NodeIndexing::GetChargeArray()
{
	int GS_XYZ=GridSize[0]*GridSize[1]*GridSize[2];
	float *q=new float[GS_XYZ];
	unsigned int i,count=0;
	
	//if(QFormat==QFormatStraightSequence)
	//{
		for(i=0;i<GS_XYZ;i++)
		{
			if(NIndex[i]&ChargeMask)
			{
				q[i]=Q[count];
				count++;
			}
			else
			{
				q[i]=0.0f;
			}
		}
//	 }
//	 else
//	 {
//	 }
	//for(i=0;i<ChargeNum[3];i++)q[IndexCharge[i]]=Q[i];
	return q;
}
/*int NodeIndexing::ConverQtoQFormatBlackNWhiteSequences()
{
	if(QFormat==QFormatBlackNWhiteSequences)return EXIT_SUCCESS;
	return EXIT_SUCCESS;
}
int NodeIndexing::ConverQtoQFormatStraightSequence()
{
	if(QFormat==QFormatStraightSequence)return EXIT_SUCCESS;
	return EXIT_SUCCESS;
}*/
int NodeIndexing::InsertSphereInDielMap(float *r, float R, int Value)
{
	if(NIndex==NULL)
	{
		fprintf(stderr,"ERROR 102: NodeIndexing::InsertSphereInDielMap: NIndex not init\n");
		return EXIT_FAILURE;
	}
	int i,j;
	//int iR=(int)roundf(R+1.0);
	int iR=(int)(R+1.5);
	int ix,iy,iz,ClosestPoint[3];
	float RSQ=R*R,distanseSQ;
	
	for(j=0;j<3;j++)
	{
		r[j]-=0.5;
		for(i=0;i<3;i++)ClosestPoint[i]=(int)(r[i]+0.5);
		for(ix=ClosestPoint[0]-iR;ix<=ClosestPoint[0]+iR;ix++)
			for(iy=ClosestPoint[1]-iR;iy<=ClosestPoint[1]+iR;iy++)
				for(iz=ClosestPoint[2]-iR;iz<=ClosestPoint[2]+iR;iz++)
		{
			distanseSQ=(ix-r[0])*(ix-r[0])+(iy-r[1])*(iy-r[1])+(iz-r[2])*(iz-r[2]);
			if(distanseSQ<=RSQ){
				SetDiel(j,ix+iy*GridSize[0]+iz*GridSize[0]*GridSize[1],Value);
			}
		}
		r[j]+=0.5;
	}
	return EXIT_SUCCESS;
}
int NodeIndexing::CheckBoarder4Q()
{
	int i,j,k;
	int itmp1,itmp2;
	int gridSizeXY=GridSize[0]*GridSize[1];
	int tag;
	int count=0;
	unsigned int ClearMask=~ChargeMask;
	//YZ face
	for(k=0;k<GridSize[2];k++)
		for(j=0;j<GridSize[1];j++)
	{
		itmp1=GridSize[0] * j + gridSizeXY * k;
		itmp2=itmp1 + GridSize[0] - 1;
		if(NIndex[itmp1]&ChargeMask)
		{
			NIndex[itmp1]=NIndex[itmp1]&ClearMask;
			count++;
		}
		if(NIndex[itmp2]&ChargeMask)
		{
			NIndex[itmp2]=NIndex[itmp2]&ClearMask;
			count++;
		}
	}
	DbgPrint0("Have found %d on the boardes (it is ok if it is MPI version)\n",count);
	//XZ face
	for(k=0;k<GridSize[2];k++)
		for(i=0;i<GridSize[0];i++)
	{
		itmp1=i + gridSizeXY *k;
		itmp2=i + GridSize[0]*(GridSize[1]-1)+ gridSizeXY *k;
		if(NIndex[itmp1]&ChargeMask)
		{
			NIndex[itmp1]=NIndex[itmp1]&ClearMask;
			count++;
		}
		if(NIndex[itmp2]&ChargeMask)
		{
			NIndex[itmp2]=NIndex[itmp2]&ClearMask;
			count++;
		}
	}
	DbgPrint0("Have found %d on the boardes (it is ok if it is MPI version)\n",count);
	//XY face
	itmp2=gridSizeXY*(GridSize[2]-1);
	for(i=0;i<gridSizeXY;i++)
	{
		if(NIndex[i]&ChargeMask)
		{
			NIndex[i]=NIndex[i]&ClearMask;
			count++;
		}
		if(NIndex[i+itmp2]&ChargeMask)
		{
			NIndex[i+itmp2]=NIndex[i+itmp2]&ClearMask;
			count++;
		}
	}
	DbgPrint0("Have found %d on the boardes (it is ok if it is MPI version)\n",count);
	return EXIT_SUCCESS;
}
void NodeIndexing::SetDiffToZero(int ion,int node)
{
	int iDzero=GetDiffZeroInd();
	NodeIndex mask=IonField[ion];
	NodeIndex sft=IonFieldSft[ion];
	unsigned int ClearMask=~mask;
	NIndex[node]=NIndex[node]&ClearMask;
	
	NodeIndex Cur=iDzero<<sft;
	NIndex[node]+=Cur;
}
///////////////////////////////////////////////////////////////////////////////
ContWorld::ContWorld()
{
	InitZero();
}
ContWorld::ContWorld(int _GridSizeX,int _GridSizeY,int _GridSizeZ,float _GridScale, bool PBCX,bool PBCY,bool PBCZ,float q1,float q2)
{
	InitZero();
	SetContWorldTwoMobIons(_GridSizeX, _GridSizeY, _GridSizeZ, _GridScale, PBCX,PBCY,PBCZ,q1, q2);
}
ContWorld::ContWorld(int* m_GridSize,float m_GridScale, bool *m_PBC,int m_NIonsTypes,float *m_Qions)
{
	InitZero();
	SetContWorld(m_GridSize, m_GridScale, m_PBC, m_NIonsTypes, m_Qions);
}
ContWorld::~ContWorld()
{
	Clear();
}
int ContWorld::InitZero()
{
	PnpsObject::SetName("ContWorld");
	BoundaryStr.push_back("Zero");
//	BoundaryStr.push_back("Coul");
	MyRank = 0;
	NProcs = 1;
	Epsilon[0]=NULL;Epsilon[1]=NULL;Epsilon[2]=NULL;
	C=NULL;
	D=NULL;
	PMF=NULL;
	Qstat=NULL;
	Potential=NULL;
	IonsQ=NULL;
	NIonsTypes=0;
	//Defult Parameters:
	MPISpreadWorld=0;
	startBlackAndWhite=0;
	NIndexing=NULL;
	PBC[0]=PBC[1]=PBC[2]=false;
	GridSize[0]=GridSize[1]=GridSize[2]=0;
	GridSizeGlobal[0]=GridSizeGlobal[1]=GridSizeGlobal[2]=0;
	GridSizeOriginal[0]=GridSizeOriginal[1]=GridSizeOriginal[2]=0;
	
	MyGlobalZ0=-1;
	MyGlobalZ1=-1;
	CDouble=NULL;
	PotentialDouble=NULL;
	return EXIT_SUCCESS;
}
int ContWorld::Clear()
{
	DeleteCArray(IonsQ);
	DeleteCArray(Epsilon[0]);
	DeleteCArray(Epsilon[1]);
	DeleteCArray(Epsilon[2]);
	DeleteCArray(Qstat);
	DeleteCArray(Potential);
	
	DeleteCVecArray(C,NIonsTypes);
	DeleteCVecArray(D,NIonsTypes);
	DeleteCVecArray(PMF,NIonsTypes);
	
	DeleteCArray(PotentialDouble);
	DeleteCVecArray(CDouble,NIonsTypes);
	
	DeleteObjByPnt(NIndexing);
	return EXIT_SUCCESS;
}
int ContWorld::SetContWorld(int* m_GridSize,float m_GridScale, bool *m_PBC,int m_NIonsTypes,float *m_Qions)
{
	int i;
	GridSizeOriginal[0]=m_GridSize[0];
	GridSizeOriginal[1]=m_GridSize[1];
	GridSizeOriginal[2]=m_GridSize[2];
	GridScale=m_GridScale;
	PBC[0]=m_PBC[0];
	PBC[1]=m_PBC[1];
	PBC[2]=m_PBC[2];

	NIonsTypes=m_NIonsTypes;
	IonsQ=new float[m_NIonsTypes];
	for(i=0;i<NIonsTypes;i++)
		IonsQ[i]=m_Qions[i];
	return CalculateDerivativeParameters();
}
int ContWorld::SetContWorldNoMobIons(int _GridSizeX,int _GridSizeY,int _GridSizeZ,float _GridScale, bool PBCX,bool PBCY,bool PBCZ)
{
	GridSizeOriginal[0]=_GridSizeX;
	GridSizeOriginal[1]=_GridSizeY;
	GridSizeOriginal[2]=_GridSizeZ;
	GridScale=_GridScale;
	PBC[0]=PBCX;
	PBC[1]=PBCY;
	PBC[2]=PBCZ;
	NIonsTypes=0;
	IonsQ=NULL;
	return CalculateDerivativeParameters();
}
int ContWorld::SetContWorldTwoMobIons(int _GridSizeX,int _GridSizeY,int _GridSizeZ,float _GridScale, bool PBCX,bool PBCY,bool PBCZ,float q1,float q2)
{
	GridSizeOriginal[0]=_GridSizeX;
	GridSizeOriginal[1]=_GridSizeY;
	GridSizeOriginal[2]=_GridSizeZ;
	GridScale=_GridScale;
	PBC[0]=PBCX;
	PBC[1]=PBCY;
	PBC[2]=PBCZ;
	NIonsTypes=2;
	IonsQ=new float[2];
	IonsQ[0]=q1;
	IonsQ[1]=q2;
	return CalculateDerivativeParameters();
}
int ContWorld::Print()
{
	int i;
	char bbool[2][10]={"False\0","True\0"};
	pnpPrintGroup0("CW    =ContWorld===============================================================\n");
	pnpPrintGroup0("CW    Original Grid Size:.................... [%d, %d, %d]\n",GridSizeOriginal[0],GridSizeOriginal[1],GridSizeOriginal[2]);
	pnpPrintGroup0("CW    Global Grid Size (with PBC padding):... [%d, %d, %d]\n",GridSizeGlobal[0],GridSizeGlobal[1],GridSizeGlobal[2]);
	//pnpPrintGroup0("CW    Local Grid Size:....................... [%d, %d, %d]\n",GridSize[0],GridSize[1],GridSize[2]);
	pnpPrintGroup0("CW    Grid Scale:............................ %.3f grids/A\n",GridScale);
	pnpPrintGroup0("CW    Periodic Boundary Conditions........... [%s, %s, %s]\n",bbool[PBC[0]],bbool[PBC[1]],bbool[PBC[2]]);
	pnpPrintGroup0("CW    Number of Mobile Ions Types:........... %d\n",NIonsTypes);

	pnpPrintGroup0("CW    Charge of Mobile Ions:................. [");
	for(i=0;i<NIonsTypes;i++){
		pnpPrintGroup0("%.1f",IonsQ[i]);
		if(i<NIonsTypes-1) pnpPrintGroup0(", ");
		else pnpPrintGroup0("]\n");
	}
		

	pnpPrintGroup0("CW    =========================================================================\n");
	return EXIT_SUCCESS;
}
int ContWorld::CalculateDerivativeParameters()
{
	//PeriodicBoundaryCondition
	if(PBC[0])
	{
		GridSizeGlobal[0] = GridSizeOriginal[0]+2;
		//fprintf(stdout, "		Increase gridsizeX at 2 due for periodical boundary contition\n");
	}
	else GridSizeGlobal[0] = GridSizeOriginal[0];
	if(PBC[1])
	{
		GridSizeGlobal[1] = GridSizeOriginal[1]+2;
		//fprintf(stdout,"		Increase gridsizeY at 2 due for periodical boundary contition\n");
	}
	else GridSizeGlobal[1] = GridSizeOriginal[1];
	if(PBC[2])
	{
		GridSizeGlobal[2] = GridSizeOriginal[2]+2;
		//fprintf(stdout,"		Increase gridsizeZ at 2 due for periodical boundary contition\n");
	}
	else GridSizeGlobal[2] = GridSizeOriginal[2];
	
	/*GridSizeGlobal[0] = GridSizeOriginal[0];
	GridSizeGlobal[1] = GridSizeOriginal[1];
	GridSizeGlobal[2] = GridSizeOriginal[2];*/
	//
/*#ifdef MPI_PARALLEL
#ifdef MDMPI
	GetBorder(&MyGlobalZ0,&MyGlobalZ1,MyRank);
	GridSize[0]=GridSizeGlobal[0];
	GridSize[1]=GridSizeGlobal[1];
	GridSize[2]=GridSizeGlobal[2];
#else
	GetBorder(&MyGlobalZ0,&MyGlobalZ1,MyRank);
	GridSize[2]=MyGlobalZ0-MyGlobalZ1+1;
#endif
#else
	GetBorder(&MyGlobalZ0,&MyGlobalZ1,MyRank);
	GridSize[0]=GridSizeGlobal[0];
	GridSize[1]=GridSizeGlobal[1];
	GridSize[2]=GridSizeGlobal[2];
#endif*/
	GetBorder(&MyGlobalZ0,&MyGlobalZ1,MyRank);
	GridSize[0]=GridSizeGlobal[0];
	GridSize[1]=GridSizeGlobal[1];
	GridSize[2]=MyGlobalZ1-MyGlobalZ0+1;
	//
	GridSizeXYZGlobal = GridSizeGlobal[0] * GridSizeGlobal[1] * GridSizeGlobal[2];
	GridSizeXYZOriginal = GridSizeOriginal[0] * GridSizeOriginal[1] * GridSizeOriginal[2];
	GS_XYZ = GridSize[0] * GridSize[1] * GridSize[2];
	GS_XY=GridSize[0] * GridSize[1];
	GS_X=GridSize[0];GS_Y=GridSize[1];GS_Z=GridSize[2];
	return EXIT_SUCCESS;
}
int ContWorld::ReadMaps(MapsIOData* Dt)
{
	int i,j;
	char ctmp2[PNP_MAP_IO_STRING_LENGTH];
	int GS[3];
	const char* ctmp;
	DefClock0;
	
	
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS/(GridScale*GridScale*GridScale);
//#ifndef MPI_PARALLEL
	if(Dt->NodeIndexingFile!="")
	{
		ReadNodeIndexing(Dt->NodeIndexingFile.c_str());
	}
	if(Dt->PotentialMapFile!="")
	{
		ReadPotential(Dt->PotentialMapFile.c_str());
	}
	if(Dt->DynamicChargeMapFile!="")
	{
		ReadDynamicCharge(Dt->DynamicChargeMapFile.c_str());
	}
	if(Dt->StaticChargeMapFile!="")
	{
		pnpError("Not Implemented (StaticChargeMap Reading)");
	}
	if(Dt->DielectricConstantMap!="")
	{
		pnpError("Not Implemented (DielectricConstantMap Reading)");
	}
	if(Dt->DiffusionMapFile!="")
	{
		ReadDiffusion(Dt->DiffusionMapFile.c_str());
	}
	if(Dt->PMFMapFile!="")
	{
		ReadPMF(Dt->PMFMapFile.c_str(),Dt->PMFMapFile2.c_str());
	}
	return EXIT_SUCCESS;
}
int ContWorld::ReadNodeIndexing(const char * filename)
{
#ifndef MPI_PARALLEL
	DeleteObjByPnt(NIndexing);
	pnpPrint("\tRead Node Indexing from:.................. %s\n",filename);
	DefClock0;
	StartClock0;
	NodeIndexing *NIndexingTMP=new NodeIndexing();
	NIndexingTMP->ReadFromFileAddPBC(filename,0,0,0);

	//i.e. if contworld not initiated then initiate from NodeIndexing
	if(GridSize[0]==0)
	{
		if(NIndexingTMP->NIonsTypes==2)
		{
			SetContWorldTwoMobIons(NIndexingTMP->GridSize[0], NIndexingTMP->GridSize[1], NIndexingTMP->GridSize[2], NIndexingTMP->GridScale, PBC[0],PBC[1],PBC[2], NIndexingTMP->IonsQ[0], NIndexingTMP->IonsQ[1]);
		}
		else
		{
			pnpError("So far can handle only two mobile ions\n");
		}
	}

	if(fabs(NIndexingTMP->GridScale-GridScale)/GridScale>0.0001)
	{
		pnpError("Grid scale of NodeIndexing is not coinside with grid scale of ContWorld\n");
		pnpError("NodeIndexing->GridScale=%g ContWorld->GridScale=%g\n",NIndexingTMP->GridScale,GridScale);
	}
	if(NIndexingTMP->GridScale==GridScale &&
			NIndexingTMP->GridSize[0]==GridSize[0] &&
			NIndexingTMP->GridSize[1]==GridSize[1] &&
			NIndexingTMP->GridSize[2]==GridSize[2])
	{
		NIndexing=NIndexingTMP;
	}
	else
	{
		pnpWarning("Grid of NodeIndexing is not coinside with ContWorld\n");
		pnpWarning("will assume you want to calculate only part of the system\n");
		pnpWarning("this=[%d %d %d] NodeIndexing=[%d %d %d]\n",GridSize[0],GridSize[1],GridSize[2],\
			NIndexingTMP->GridSize[0],NIndexingTMP->GridSize[1],NIndexingTMP->GridSize[2]);
		
		unsigned int gridsize[3]={GridSize[0],GridSize[1],GridSize[2]};
		NIndexing=new NodeIndexing();
		NIndexing->SetNNodes(gridsize,GridScale);
		NIndexing->GetCentralPartFromNIExtMyPart4MPI(NIndexingTMP);
		NIndexing->SetBlackAndWhiteNodes(startBlackAndWhite);
		delete NIndexingTMP;
	}
	StopClockWMes0("Read Node Indexing");
#else
	DeleteObjByPnt(NIndexing);
	if(pnpsapp->AmIGroupLeader())
		pnpPrint("\tRead Node Indexing from:.................. %s\n",filename);
	DefClock0;
	StartClock0;
	NodeIndexing *NIndexingTMP=NULL;
	NodeIndexing *NIndexingTMP2=NULL;
	if(pnpsapp->AmIGroupLeader())
	{
		NIndexingTMP=new NodeIndexing();
		NIndexingTMP->ReadFromFileAddPBC(filename,PBC[0],PBC[1],PBC[2]);
		if(fabs(NIndexingTMP->GridScale-GridScale)/GridScale>0.0001)
		{
			pnpError("Grid scale of NodeIndexing is not coinside with grid scale of ContWorld\n");
			pnpError("NodeIndexing->GridScale=%g ContWorld->GridScale=%g\n",NIndexingTMP->GridScale,GridScale);
			return EXIT_FAILURE;
		}
		if(NIndexingTMP->GridScale==GridScale &&
				NIndexingTMP->GridSize[0]==GridSizeGlobal[0] &&
				NIndexingTMP->GridSize[1]==GridSizeGlobal[1] &&
				NIndexingTMP->GridSize[2]==GridSizeGlobal[2])
		{
			
		}
		else
		{
			pnpWarning("Grid of NodeIndexing is not coinside with ContWorld\n");
			pnpWarning("will assume you want to calculate only part of the system\n");
			
			unsigned int gridsize[3]={GridSizeGlobal[0],GridSizeGlobal[1],GridSizeGlobal[2]};
			NIndexingTMP2=new NodeIndexing();
			NIndexingTMP2->SetNNodes(gridsize,GridScale);
			NIndexingTMP2->GetCentralPartFromNIExtMyPart4MPI(NIndexingTMP);
			NIndexingTMP2->SetBlackAndWhiteNodes(startBlackAndWhite);
			delete NIndexingTMP;
			NIndexingTMP=NIndexingTMP2;
		}
	}
	if(pnpsapp->GetNumProcsInMyGroup()==1)
		NIndexing=NIndexingTMP;
	else
	{
		NIndexing=new NodeIndexing();
		NIndexing->GetMyPart4MPI(this,NIndexingTMP);
		DeleteObjByPnt(NIndexingTMP);
	}
	StopClockWMes0("Read Node Indexing");
	
#endif
	return EXIT_SUCCESS;
}
int ContWorld::WriteContTop(const char * filename,int opt)
{
	return WriteNodeIndexing(filename,opt);
}
int ContWorld::WriteNodeIndexing(const char * filename,int opt)
{
	if(opt==1)//Original grid
	{
#ifndef MPI_PARALLEL
		pnpPrint0("\tWrite Node Indexing to:.................. %s\n",filename);
		NIndexing->SetNoPBC();
		NIndexing->WriteToFile(filename);
		NIndexing->SetPBC(PBC[0],PBC[1],PBC[2]);
#else
		if(pnpsapp->GetNumProcsInMyGroup()==1)
		{
			char FileNameLoc[PNP_MAP_IO_STRING_LENGTH];
			//pnpsapp->AddMyGroupNumberToFileName(FileNameLoc,filename);
			pnpPrint0("\tWrite Node Indexing of proc %d to:....... %s\n", pnpsapp->GetMyAbsRank(), FileNameLoc);
			NIndexing->SetNoPBC();
			NIndexing->WriteToFile(FileNameLoc);
			NIndexing->SetPBC(PBC[0],PBC[1],PBC[2]);
		}
#endif
		return EXIT_SUCCESS;
	}
	else if (opt==2)//Local grid
	{
#ifndef MPI_PARALLEL
		pnpPrint0("\tWrite Node Indexing to:.................. %s\n",filename);
		NIndexing->WriteToFile(filename);
#else
		char FileNameLoc[PNP_MAP_IO_STRING_LENGTH];
		pnpsapp->AddMyAbsProcNumberToFileName(FileNameLoc,filename);
		pnpPrint0("\tWrite Node Indexing of proc %d to:....... %s\n", pnpsapp->GetMyAbsRank(), FileNameLoc);
		NIndexing->WriteToFile(FileNameLoc);
#endif
		return EXIT_SUCCESS;
	}
	else
	{
		pnpError("Option for writing NodeIndexing is not known (%d)\n",opt);
		return EXIT_FAILURE;
	}
}
int ContWorld::ReadDiffusion(const char * filename)
{
	pnpPrint0("\tRead Diffusion from:......... %s\n", filename);
	int i,j;
	DeleteCVecArray(D,NIonsTypes);
	
	VectorField3D* VF3D=new VectorField3D();
#ifndef MPI_PARALLEL
	VF3D->ReadFromFileAddPBC(filename, 1.0,PBC[0],PBC[1],PBC[2]);
	VF3D=CheckGridOfVectorField3D(VF3D,true,"Diffusion",filename);
#else
	if(pnpsapp->GetNumberOfGroups()==pnpsapp->GetTotalProc())
	{
		VF3D->ReadFromFileAddPBC(filename, 1.0,PBC[0],PBC[1],PBC[2]);
		VF3D=CheckGridOfVectorField3D(VF3D,true,"Diffusion",filename);
	}
	else
	{
		VectorField3D* VF3D2=NULL;
		if(pnpsapp->AmIGroupLeader())
		{
			VF3D2=new VectorField3D();
			VF3D2->ReadFromFileAddPBC(filename, 1.0,PBC[0],PBC[1],PBC[2]);
			VF3D2=CheckGridOfVectorField3D(VF3D2,true,"Diffusion",filename);
		}
		VF3D->SplitExtVF3DSDMPI4FDaZ(VF3D2,MyGlobalZ0,MyGlobalZ1);
		delete VF3D2;
	}
#endif
	if(VF3D->GetNelem()==NIonsTypes)
	{
		D=VF3D->V;
		VF3D->amode=VectorField3D::EXTERNAL_ALLOC;
		delete VF3D;
	}
	else if(VF3D->GetNelem()==1)
	{
		pnpWarning("only one component at Diffusion will assume is the same for all ions\n");
		D=new float*[NIonsTypes];
		for(j=0;j<NIonsTypes;j++)
		{
			D[j]=new float[GS_XYZ];
			for(i=0;i<GS_XYZ;i++)
				D[j][i]=VF3D->V[0][i];
		}
		delete VF3D;
	}
	else
	{
		pnpWarning("Something wrong with Diffusions\n");
	}
	return EXIT_SUCCESS;
}
int ContWorld::WriteDiffusion(const char * filename,int opt)
{
	int i;
	const char * filename2=filename;
	PNP_EXIT_FAIL_NULL(D,"Diffusion Map NOT yet initialize\n");
	for(i=0;i<NIonsTypes;i++)
		PNP_EXIT_FAIL_NULL1(D[i],"Diffusion Map for ion %d NOT yet initialize\n",i);
	
	if(opt==1)//Original grid
	{
		VectorField3D *VField3DContr=new VectorField3D(GridSize,GridScale,NIonsTypes,D);
#ifndef MPI_PARALLEL
		pnpPrint0("\tWrite Diffusion Map to:.. %s\n",filename);
		VField3DContr->WriteToFileRemovePBC(filename,1.0,1,PBC[0],PBC[1],PBC[2]);
#else
		if(pnpsapp->GetNumProcsInMyGroup()==1)
		{
			char FileNameLoc[PNP_MAP_IO_STRING_LENGTH];
			//pnpsapp->AddMyGroupNumberToFileName(FileNameLoc,filename);
			pnpPrint0("\tWrite Diffusion Map of proc %d to: %s\n", pnpsapp->GetMyAbsRank(), FileNameLoc);
			VField3DContr->WriteToFileRemovePBC(FileNameLoc,1.0,1,PBC[0],PBC[1],PBC[2]);
		}
		else
		{
			pnpPrint0("\tWrite Diffusion Map to:.. %s\n",filename);
			VField3DContr->WriteToFileCombineSDMPIDistrRemovePBC(filename,1.0,1,GridSizeGlobal,PBC,MyGlobalZ0,MyGlobalZ1);
		}
#endif
		delete VField3DContr;
		return EXIT_SUCCESS;
	}
	else if (opt==2)//Local grid
	{
		VectorField3D *VField3DContr=new VectorField3D(GridSize,GridScale,NIonsTypes,PMF);
#ifndef MPI_PARALLEL
		pnpPrint0("\tWrite Diffusion Map to:.. %s\n",filename);
		VField3DContr->WriteToFile(filename,1.0,1);
#else
		char FileNameLoc[PNP_MAP_IO_STRING_LENGTH];
		pnpsapp->AddMyAbsProcNumberToFileName(FileNameLoc,filename);
		pnpPrint0("\tWrite Diffusion Map of proc %d to: %s\n", pnpsapp->GetMyAbsRank(), FileNameLoc);
		VField3DContr->WriteToFile(FileNameLoc,1.0,1);
#endif
		delete VField3DContr;
		return EXIT_SUCCESS;
	}
	else
	{
		pnpError("Option for writing Diffusion is not known (%d)\n",opt);
		return EXIT_FAILURE;
	}
}
VectorField3D* ContWorld::CheckGridOfVectorField3D(VectorField3D* VF3D,bool deleteOldVF3DifNotSame,const char * VF3DDescription,const char * ReadFromFile)
{
	int i,j,k;
	bool ChangedVF3D=false;
	VectorField3D* VF3D2;
	if(fabs(VF3D->GridScale-GridScale)>GridScale*0.001)
	{
		pnpWarning("Grid Size and/or Grid Scale of %s from %s is not coinside with ContWorld\n",VF3DDescription,ReadFromFile);
		pnpWarning("Grid Size: %s=[%d, %d, %d] [g, g, g], ContWorld(Global)=[%d, %d, %d] [g, g, g]\n",ReadFromFile,VF3D->GridSize[0],VF3D->GridSize[1],VF3D->GridSize[2],GridSizeGlobal[0],GridSizeGlobal[1],GridSizeGlobal[2]);
		pnpWarning("Grid Scale: %s=%g g/A, ContWorld=%g g/A\n",ReadFromFile,VF3D->GridScale,GridScale);
		
		
		if(VF3D->GridSize[0]/VF3D->GridScale>=GridSize[0]/GridScale&&
				 VF3D->GridSize[1]/VF3D->GridScale>=GridSize[1]/GridScale&&
				 VF3D->GridSize[2]/VF3D->GridScale>=GridSize[2]/GridScale)
		{
			pnpWarning("Will assume you are doing focusing\n");
			VF3D2=new VectorField3D(GridSizeGlobal,GridScale,VF3D->Nelem);
			VF3D2->InterpolateFromExt(VF3D);
			ChangedVF3D=true;
		}
		else
		{
			pnpWarning("Will assume you are reading internal part of system calculated with bigger gridscale, work on single CPU\n");
			VF3D2=new VectorField3D(GridSizeGlobal,GridScale,VF3D->Nelem);
			VF3D2->InterpolateFromExt(VF3D);
			ChangedVF3D=true;
		}
	}
	else if(VF3D->GridSize[0]!=GridSizeGlobal[0] || VF3D->GridSize[1]!=GridSizeGlobal[1] || VF3D->GridSize[2]!=GridSizeGlobal[2])
	{
		pnpWarning("Grid Size of %s from %s is not coinside with ContWorld\n",VF3DDescription,ReadFromFile);
		pnpWarning("Grid Size: %s=[%d, %d, %d] [g, g, g], ContWorld(Global)=[%d, %d, %d] [g, g, g]\n",ReadFromFile,VF3D->GridSize[0],VF3D->GridSize[1],VF3D->GridSize[2],GridSizeGlobal[0],GridSizeGlobal[1],GridSizeGlobal[2]);
		pnpWarning("Will assume you are calculating only part of the system\n");
		VF3D2=new VectorField3D(GridSizeGlobal,GridScale,VF3D->Nelem);
		VF3D2->InterpolateFromExt(VF3D);
		ChangedVF3D=true;
	}
	if(ChangedVF3D)
	{
		if(deleteOldVF3DifNotSame)
			delete VF3D;
		return VF3D2;
	}
	else
	{
		return VF3D;
	}
}
int ContWorld::ReadPotentialChargeZRange(const char * filename,float Z0,float Z1)
{
	VectorField3D* VF3D=ReadVectorField3D(filename,1.0);
	int i,j;
	int count=0;
	Z0=ConvFloatToGlobIntUnitsZ(Z0);
	Z1=ConvFloatToGlobIntUnitsZ(Z1);
	int ix,iy,iz,GrdPnt;
	int iZ0,iZ1;
	iZ0=roundf(Z0);
	iZ1=roundf(Z1);
	if(iZ0<0)iZ0=0;
	if(iZ1>GS_Z-1)iZ1=GS_Z-1;
	if((iZ0>iZ1)||(iZ0>GS_Z-1)||(iZ1<0))
	{
		pnpWarning("z-range out of simulation box for ReadDynamicChargeZRange\n");
		delete VF3D;
		return EXIT_SUCCESS;
	}
		
	ix=GS_X-1;
	for(iz=iZ0;iz<=iZ1;iz++)
		for(iy=0;iy<GS_Y;iy++)
			for(ix=0;ix<GS_X;ix++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		Potential[GrdPnt]=VF3D->V[0][GrdPnt];
	}
	return EXIT_SUCCESS;
}
int ContWorld::ReadPotential(const char * filename)
{
	pnpPrint0("\tRead potential from:......... %s\n", filename);
	
	DeleteCArray(Potential);
	int i,j;
	
	VectorField3D* VF3D=new VectorField3D();
#ifndef MPI_PARALLEL
	VF3D->ReadFromFileAddPBC(filename, 1.0,PBC[0],PBC[1],PBC[2]);
	VF3D=CheckGridOfVectorField3D(VF3D,true,"potential",filename);
#else
	if(pnpsapp->GetNumberOfGroups()==pnpsapp->GetTotalProc())
	{
		VF3D->ReadFromFileAddPBC(filename, 1.0,PBC[0],PBC[1],PBC[2]);
		VF3D=CheckGridOfVectorField3D(VF3D,true,"potential",filename);
	}
	else
	{
		VectorField3D* VF3D2=NULL;
		if(pnpsapp->AmIGroupLeader())
		{
			VF3D2=new VectorField3D();
			VF3D2->ReadFromFileAddPBC(filename, 1.0,PBC[0],PBC[1],PBC[2]);
			VF3D2=CheckGridOfVectorField3D(VF3D2,true,"potential",filename);
		}
		VF3D->SplitExtVF3DSDMPI4FDaZ(VF3D2,MyGlobalZ0,MyGlobalZ1);
		delete VF3D2;
	}
#endif
	Potential=VF3D->V[0];
	VF3D->amode=VectorField3D::EXTERNAL_ALLOC;
	delete [] VF3D->V;
	delete VF3D;
	
	return EXIT_SUCCESS;
}
int ContWorld::WritePotential(const char * filename,int opt)
{
	PNP_EXIT_FAIL_NULL(Potential,"Potential Map NOT yet initialize\n");

	float *tmp[1];
	tmp[0]=Potential;
	
	
	if(opt==1)//Original grid
	{
		VectorField3D *VField3DContr=new VectorField3D(GridSize,GridScale,1,tmp);
#ifndef MPI_PARALLEL
		pnpPrint0("\tWrite Potential Map to:................ %s\n",filename);
		VField3DContr->WriteToFileRemovePBC(filename,1.0,1,PBC[0],PBC[1],PBC[2]);
#else
		if(pnpsapp->GetNumProcsInMyGroup()==1)
		{
			char FileNameLoc[PNP_MAP_IO_STRING_LENGTH];
			//pnpsapp->AddMyGroupNumberToFileName(FileNameLoc,filename);
			pnpPrint0("\tWrite Potential Map of proc %d to:..... %s\n", pnpsapp->GetMyAbsRank(), filename);
			VField3DContr->WriteToFileRemovePBC(filename,1.0,1,PBC[0],PBC[1],PBC[2]);
		}
		else
		{
			pnpPrintMain0("\tWrite Potential Map to:................ %s\n",filename);
			VField3DContr->WriteToFileCombineSDMPIDistrRemovePBC(filename,1.0,1,GridSizeGlobal,PBC,MyGlobalZ0,MyGlobalZ1);
		}
#endif
		delete VField3DContr;
		return EXIT_SUCCESS;
	}
	else if (opt==2)//Local grid
	{
		VectorField3D *VField3DContr=new VectorField3D(GridSize,GridScale,1,tmp);
#ifndef MPI_PARALLEL
		pnpPrint0("\tWrite Potential Map to:................ %s\n",filename);
		VField3DContr->WriteToFile(filename,1.0,1);
#else
		//char FileNameLoc[PNP_MAP_IO_STRING_LENGTH];
		//pnpsapp->AddMyAbsProcNumberToFileName(FileNameLoc,filename);
		pnpPrint0("\tWrite Potential Map of proc %d to:..... %s\n", pnpsapp->GetMyAbsRank(), filename);
		VField3DContr->WriteToFile(filename,1.0,1);
#endif
		delete VField3DContr;
		return EXIT_SUCCESS;
	}
	else
	{
		pnpError("Option for writing Potential Map is not known (%d)\n",opt);
		return EXIT_FAILURE;
	}
	return EXIT_FAILURE;
}
int ContWorld::ReadPMF(const char * filename,const char * filename2)
{
	pnpPrint0("\tRead PMF from:......... %s\n", filename);
	
	int i,j;
	DeleteCVecArray(PMF,NIonsTypes);
	
	VectorField3D* VF3D=new VectorField3D();
#ifndef MPI_PARALLEL
	VF3D->ReadFromFileAddPBC(filename, 1.0,PBC[0],PBC[1],PBC[2]);
	VF3D=CheckGridOfVectorField3D(VF3D,true,"PMF",filename);
	
	if(strcmp(filename2,"")!=0)
	{
		pnpPrint0("\twill add PMF from:......... %s\n", filename2);
		VectorField3D* VF3D2=new VectorField3D();
		VF3D2->ReadFromFileAddPBC(filename2, 1.0,PBC[0],PBC[1],PBC[2]);
		VF3D2=CheckGridOfVectorField3D(VF3D2,true,"PMF",filename2);
		for(j=0;j<VF3D->Nelem;j++)
		{
			for(i=0;i<GridSizeXYZGlobal;i++)
				VF3D->V[j][i]+=VF3D2->V[j][i];
		}
		delete VF3D2;
	}
#else
	if(pnpsapp->GetNumberOfGroups()==pnpsapp->GetTotalProc())
	{
		VF3D->ReadFromFileAddPBC(filename, 1.0,PBC[0],PBC[1],PBC[2]);
		VF3D=CheckGridOfVectorField3D(VF3D,true,"PMF",filename);
	}
	else
	{
		VectorField3D* VF3D2=NULL;
		if(pnpsapp->AmIGroupLeader())
		{
			VF3D2=new VectorField3D();
			VF3D2->ReadFromFileAddPBC(filename, 1.0,PBC[0],PBC[1],PBC[2]);
			VF3D2=CheckGridOfVectorField3D(VF3D2,true,"PMF",filename);
		}
		VF3D->SplitExtVF3DSDMPI4FDaZ(VF3D2,MyGlobalZ0,MyGlobalZ1);
		delete VF3D2;
	}
#endif
	if(VF3D->GetNelem()==NIonsTypes)
	{
		PMF=VF3D->V;
		VF3D->amode=VectorField3D::EXTERNAL_ALLOC;
		delete VF3D;
	}
	else if(VF3D->GetNelem()==1)
	{
		pnpWarning("only one component at PMF will assume it is electrostatic potential\n");
		PMF=new float*[NIonsTypes];
		for(j=0;j<NIonsTypes;j++)
		{
			PMF[j]=new float[GS_XYZ];
			for(i=0;i<GS_XYZ;i++)
				PMF[j][i]=VF3D->V[0][i]*IonsQ[j];
		}
		delete VF3D;
	}
	else
	{
		pnpWarning("Something wrong with PMF\n");
	}
	
	/*int sz1,sz2,GrdPnt;
	GetBorder(&sz1,&sz2,pnpsapp->GetMyRankInGroup());
	for(i=0;i<GS_Z;i++)
	{
		GrdPnt=int(GS_X/2)+int(GS_Y/2)*GS_X+i*GS_XY;
		pnpPrint("%d ",i+sz1);
		for(j=0;j<NIonsTypes;j++)
			pnpPrint("%g ",PMF[j][GrdPnt]);
		pnpPrint("\n");
	}*/
	/*char filenameTMP_NI[1024]="TMP.bin";
	pnpsapp->AddMyAbsProcNumberToFileName(filenameTMP_NI,filenameTMP_NI);
	WritePMF(filenameTMP_NI,2);*/
	return EXIT_SUCCESS;
}
int ContWorld::WritePMF(const char * filename,int opt)
{
	int i;
	const char * filename2=filename;
	PNP_EXIT_FAIL_NULL(PMF,"Potential of Mean Field Map NOT yet initialize\n");
	for(i=0;i<NIonsTypes;i++)
		PNP_EXIT_FAIL_NULL1(PMF[i],"Potential of Mean Field Map for ion %d NOT yet initialize\n",i);
	
	if(opt==1)//Original grid
	{
		VectorField3D *VField3DContr=new VectorField3D(GridSize,GridScale,NIonsTypes,PMF);
#ifndef MPI_PARALLEL
		pnpPrint0("\tWrite Potential of Mean Field Map to:.. %s\n",filename);
		VField3DContr->WriteToFileRemovePBC(filename,1.0,1,PBC[0],PBC[1],PBC[2]);
#else
		if(pnpsapp->GetNumProcsInMyGroup()==1)
		{
			char FileNameLoc[PNP_MAP_IO_STRING_LENGTH];
			//pnpsapp->AddMyGroupNumberToFileName(FileNameLoc,filename);
			pnpPrint0("\tWrite Potential of Mean Field Map of proc %d to: %s\n", pnpsapp->GetMyAbsRank(), FileNameLoc);
			VField3DContr->WriteToFileRemovePBC(FileNameLoc,1.0,1,PBC[0],PBC[1],PBC[2]);
		}
		else
		{
			pnpPrint0("\tWrite Potential of Mean Field Map to:.. %s\n",filename);
			VField3DContr->WriteToFileCombineSDMPIDistrRemovePBC(filename,1.0,1,GridSizeGlobal,PBC,MyGlobalZ0,MyGlobalZ1);
		}
#endif
		delete VField3DContr;
		return EXIT_SUCCESS;
	}
	else if (opt==2)//Local grid
	{
		VectorField3D *VField3DContr=new VectorField3D(GridSize,GridScale,NIonsTypes,PMF);
#ifndef MPI_PARALLEL
		pnpPrint0("\tWrite Potential of Mean Field Map to:.. %s\n",filename);
		VField3DContr->WriteToFile(filename,1.0,1);
#else
		char FileNameLoc[PNP_MAP_IO_STRING_LENGTH];
		pnpsapp->AddMyAbsProcNumberToFileName(FileNameLoc,filename);
		pnpPrint0("\tWrite Potential of Mean Field Map of proc %d to: %s\n", pnpsapp->GetMyAbsRank(), FileNameLoc);
		VField3DContr->WriteToFile(FileNameLoc,1.0,1);
#endif
		delete VField3DContr;
		return EXIT_SUCCESS;
	}
	else
	{
		pnpError("Option for writing NodeIndexing is not known (%d)\n",opt);
		return EXIT_FAILURE;
	}
	return EXIT_FAILURE;
}
VectorField3D* ContWorld::ReadVectorField3D(const char * filename,float coef)
{
	pnpPrint0("\tRead VectorField3D from:......... %s\n", filename);
	int i,j;
	
	VectorField3D* VF3D=new VectorField3D();
#ifndef MPI_PARALLEL
	VF3D->ReadFromFileAddPBC(filename, coef,PBC[0],PBC[1],PBC[2]);
	VF3D=CheckGridOfVectorField3D(VF3D,true,"ions' concentrations",filename);
#else
	if(pnpsapp->GetNumberOfGroups()==pnpsapp->GetTotalProc())
	{
		VF3D->ReadFromFileAddPBC(filename, coef,PBC[0],PBC[1],PBC[2]);
		VF3D=CheckGridOfVectorField3D(VF3D,true,"ions' concentrations",filename);
	}
	else
	{
		VectorField3D* VF3D2=NULL;
		if(pnpsapp->AmIGroupLeader())
		{
			VF3D2=new VectorField3D();
			VF3D2->ReadFromFileAddPBC(filename, coef,PBC[0],PBC[1],PBC[2]);
			VF3D2=CheckGridOfVectorField3D(VF3D2,true,"ions' concentrations",filename);
		}
		VF3D->SplitExtVF3DSDMPI4FDaZ(VF3D2,MyGlobalZ0,MyGlobalZ1);
		delete VF3D2;
	}
#endif
	return VF3D;
}
int ContWorld::ReadDynamicChargeZRange(const char * filename,float Z0,float Z1)
{
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS/(GridScale*GridScale*GridScale);
	VectorField3D* VF3D=ReadVectorField3D(filename,coef);
	int i,j;
	int count=0;
	if(Z0>Z1)
	{
		float tmp=Z0;
		Z0=Z1;
		Z1=tmp;
	}
	
	DbgPrint0("ReadDynamicChargeZRange Z0 = %g A, Z1 = %g A\n",Z0,Z1);
	Z0=ConvFloatToGlobIntUnitsZ(Z0);
	Z1=ConvFloatToGlobIntUnitsZ(Z1);
	
	if(NIndexing!=NULL)
	{
		for(j=0;j<NIonsTypes;j++)
		{
			for(i=0;i<GS_XYZ;i++)
			{
				if(NIndexing->GetDiffFloat(j,i)==0.0&&C[j][i]>0.0)
				{
					count++;
					VF3D->V[j][i]=0.0;
				}
			}
		}
		if(count>0)
			pnpWarning("%d positions was found where C!=0 and D==0, have set them to 0\n",count);
	}
	int ix,iy,iz,GrdPnt,GrdPntGlob;
	int iGlobZ0,iGlobZ1;
	iGlobZ0=roundf(Z0);
	iGlobZ1=roundf(Z1);
	DbgPrint0("ReadDynamicChargeZRange Z0 = %d [int Glob], Z1 = %d [int Glob]\n",iGlobZ0,iGlobZ1);
	
	if(iGlobZ0>GridSizeGlobal[2]-1||iGlobZ1<0)
	{
		pnpWarning("z-range out of simulation box no data will be readed\n");
		delete VF3D;
		return EXIT_SUCCESS;
	}
	
	if(iGlobZ0<0)iGlobZ0=0;
	if(iGlobZ1>GridSizeGlobal[2]-1)iGlobZ1=GridSizeGlobal[2]-1;
	
	DbgPrint0("ReadDynamicChargeZRange Z0 = %d [int Glob], Z1 = %d [int Glob]\n",iGlobZ0,iGlobZ1);
	
	int MyZ1, MyZ2;
	GetBorder(&MyZ1,&MyZ2,0);
	int iZ0,iZ1;
	iZ0=iGlobZ0-MyZ1;
	iZ1=iGlobZ1-MyZ1;
	if(iZ0>GridSize[2]-1||iZ1<0)
	{
		pnpPrint0("z-range out of this portion of simulation box\n");
		delete VF3D;
		return EXIT_SUCCESS;
	}
	DbgPrint0("ReadDynamicChargeZRange Z0 = %d [int loc], Z1 = %d [int loc]\n",iZ0,iZ1);
	if(iZ0<0)iZ0=0;
	if(iZ1>GS_Z-1)iZ1=GS_Z-1;
	DbgPrint0("ReadDynamicChargeZRange Z0 = %d [int loc], Z1 = %d [int loc]\n",iZ0,iZ1);
	
	
	
	for(j=0;j<NIonsTypes;j++)
		for(iz=iZ0;iz<=iZ1;iz++)
			for(iy=0;iy<GS_Y;iy++)
				for(ix=0;ix<GS_X;ix++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		C[j][GrdPnt]=VF3D->V[j][GrdPnt];
	}
	delete VF3D;
	return EXIT_SUCCESS;
}
int ContWorld::ReadDynamicCharge(const char * filename)
{
	pnpPrint0("\tRead Dynamic Charge from:......... %s\n", filename);
	int i,j;
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS/(GridScale*GridScale*GridScale);
	DeleteCVecArray(C,NIonsTypes);
	
	VectorField3D* VF3D=new VectorField3D();
#ifndef MPI_PARALLEL
	VF3D->ReadFromFileAddPBC(filename, coef,PBC[0],PBC[1],PBC[2]);
	VF3D=CheckGridOfVectorField3D(VF3D,true,"ions' concentrations",filename);
#else
	if(pnpsapp->GetNumberOfGroups()==pnpsapp->GetTotalProc())
	{
		VF3D->ReadFromFileAddPBC(filename, coef,PBC[0],PBC[1],PBC[2]);
		VF3D=CheckGridOfVectorField3D(VF3D,true,"ions' concentrations",filename);
	}
	else
	{
		VectorField3D* VF3D2=NULL;
		if(pnpsapp->AmIGroupLeader())
		{
			VF3D2=new VectorField3D();
			VF3D2->ReadFromFileAddPBC(filename, coef,PBC[0],PBC[1],PBC[2]);
			VF3D2=CheckGridOfVectorField3D(VF3D2,true,"ions' concentrations",filename);
		}
		VF3D->SplitExtVF3DSDMPI4FDaZ(VF3D2,MyGlobalZ0,MyGlobalZ1);
		delete VF3D2;
	}
#endif
	if(VF3D->GetNelem()==NIonsTypes)
	{
		C=VF3D->V;
		VF3D->amode=VectorField3D::EXTERNAL_ALLOC;
		delete VF3D;
	}
	else if(VF3D->GetNelem()==1)
	{
		pnpWarning("only one component for concentrations, will copy it for others ions types\n");
		C=new float*[NIonsTypes];
		for(j=0;j<NIonsTypes;j++)
		{
			C[j]=new float[GS_XYZ];
			for(i=0;i<GS_XYZ;i++)
				C[j][i]=VF3D->V[0][i];
		}
		delete VF3D;
	}
	//check if C==0 where D==0
	int count=0;
	DbgPrint0("GS_XYZ=%d\n",GS_XYZ);
	if(NIndexing!=NULL)
	{
		for(j=0;j<NIonsTypes;j++)
		{
			for(i=0;i<GS_XYZ;i++)
			{
				if(NIndexing->GetDiffFloat(j,i)==0.0&&C[j][i]>0.0)
				{
					count++;
					C[j][i]=0.0;
				}
			}
		}
		if(count>0)
			pnpWarning("%d positions was found where C!=0 and D==0, have set them to 0\n",count);
	}
	
	return EXIT_SUCCESS;
}
int ContWorld::WriteDynamicCharge(const char * filename,int opt)
{
	int i;
	const char * filename2=filename;
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS/(GridScale*GridScale*GridScale);
	if(C==NULL)
	{
		pnpWarning("Dynamic Charge Map NOT yet initialize will set initial guess\n");
		SetInitConcentrationFromNIndexing();
	}
	
	PNP_EXIT_FAIL_NULL(C,"Dynamic Charge Map NOT yet initialize\n");
	for(i=0;i<NIonsTypes;i++)
		PNP_EXIT_FAIL_NULL1(C[i],"Dynamic Charge Map for ion %d NOT yet initialize\n",i);
	
	if(opt==1)//Original grid
	{
		VectorField3D *VField3DContr=new VectorField3D(GridSize,GridScale,NIonsTypes,C);
#ifndef MPI_PARALLEL
		pnpPrint0("\tWrite Dynamic Charge Map to:........... %s\n",filename);
		VField3DContr->WriteToFileRemovePBC(filename,1.0/coef,1,PBC[0],PBC[1],PBC[2]);
#else
		if(pnpsapp->GetNumProcsInMyGroup()==1)
		{
			//char FileNameLoc[PNP_MAP_IO_STRING_LENGTH];
			//pnpsapp->AddMyGroupNumberToFileName(FileNameLoc,filename);
			pnpPrint0("\tWrite Dynamic Charge Map of proc %d to: %s\n", pnpsapp->GetMyAbsRank(), filename);
			VField3DContr->WriteToFileRemovePBC(filename,1.0/coef,1,PBC[0],PBC[1],PBC[2]);
		}
		else
		{
			pnpPrint0("\tWrite Dynamic Charge Map to:........... %s\n",filename);
			VField3DContr->WriteToFileCombineSDMPIDistrRemovePBC(filename,1.0/coef,1,GridSizeGlobal,PBC,MyGlobalZ0,MyGlobalZ1);
		}
#endif
		delete VField3DContr;
		return EXIT_SUCCESS;
	}
	else if (opt==2)//Local grid
	{
		VectorField3D *VField3DContr=new VectorField3D(GridSize,GridScale,NIonsTypes,C);
#ifndef MPI_PARALLEL
		pnpPrint0("\tWrite Dynamic Charge Map to:........... %s\n",filename);
		VField3DContr->WriteToFile(filename,1.0/coef,1);
#else
		//char FileNameLoc[PNP_MAP_IO_STRING_LENGTH];
		//pnpsapp->AddMyAbsProcNumberToFileName(FileNameLoc,filename);
		pnpPrint0("\tWrite Dynamic Charge Map of proc %d to: %s\n", pnpsapp->GetMyAbsRank(), filename);
		VField3DContr->WriteToFile(filename,1.0/coef,1);
#endif
		delete VField3DContr;
		return EXIT_SUCCESS;
	}
	else
	{
		pnpError("Option for writing NodeIndexing is not known (%d)\n",opt);
		return EXIT_FAILURE;
	}
	return EXIT_FAILURE;
}
int ContWorld::WriteMaps(MapsIOData* Dt)
{
	pnpPrint0("Going to write maps\n");
	//!@todo converting PBC implimenting MPI in VectorField3D
	const char* ctmp;
	VectorField3D* VField3DContr=NULL;
	float fpoh = 4*M_PI*GridScale;
	float temp = fpoh*COANGS/(GridScale*GridScale*GridScale);
	
	DbgPrint2("DielectricConstantMap : \"%s\"\n",Dt->DielectricConstantMap.c_str());
	DbgPrint2("DiffusionMapFile : \"%s\"\n",Dt->DiffusionMapFile.c_str());
	DbgPrint2("StaticChargeMapFile : \"%s\"\n",Dt->StaticChargeMapFile.c_str());
	DbgPrint2("PMFMapFile : \"%s\"\n",Dt->PMFMapFile.c_str());
	DbgPrint2("NodeIndexingFile : \"%s\"\n",Dt->NodeIndexingFile.c_str());
	DbgPrint2("PotentialMapFile : \"%s\"\n",Dt->PotentialMapFile.c_str());
	DbgPrint2("DynamicChargeMapFile : \"%s\"\n",Dt->DynamicChargeMapFile.c_str());
	
	if(Dt->NodeIndexingFile!="")
	{
		WriteNodeIndexing(Dt->NodeIndexingFile.c_str());
	}
	if(Dt->PotentialMapFile!="")
	{
		WritePotential(Dt->PotentialMapFile.c_str());
	}
	if(Dt->DynamicChargeMapFile!="")
	{
		WriteDynamicCharge(Dt->DynamicChargeMapFile.c_str());
	}
	if(Dt->StaticChargeMapFile!="")
	{
		if(Qstat)
		{
			ctmp=Dt->StaticChargeMapFile.c_str();
			fprintf(stdout,"		Write Static Charge Map to:................ %s\n",ctmp);
			float *tmp[1];
			tmp[0]=Qstat;
			VField3DContr=new VectorField3D(GridSize,GridScale,1,tmp);
			VField3DContr->WriteToFileRemovePBC(ctmp,1.0/fpoh,2,PBC[0],PBC[1],PBC[2]);
			delete VField3DContr;VField3DContr=NULL;
		}
		else{
			fprintf(stderr, "		StaticCharge Map NOT yet initialize\n");
		}
	}
	if(Dt->DielectricConstantMap!="")
	{
		if(Epsilon[0]&&Epsilon[1]&&Epsilon[2]) 
		{
			ctmp=Dt->DielectricConstantMap.c_str();
			fprintf(stdout,"		Write DielectricConstant Map to:........... %s\n",ctmp);
			VField3DContr=new VectorField3D(GridSize,GridScale,3,Epsilon);
			VField3DContr->WriteToFileRemovePBC(ctmp,EPKT,2,PBC[0],PBC[1],PBC[2]);
			delete VField3DContr;VField3DContr=NULL;
		}
		else fprintf(stderr, "ERR:		DielectricConstant NOT yet initialize\n");
	}
	if(Dt->DiffusionMapFile!="")
	{
		if(D!=NULL)
		{
			WriteDiffusion(Dt->DiffusionMapFile.c_str());
		}
		else 
		{
			fprintf(stderr, "ERR:		Diffusion Map NOT yet initialize\n");
			if(NIndexing!=NULL)
			{
				pnpPrint0("Will convert from NIndexing\n");
				D=new float*[NIonsTypes];
				if(NIonsTypes>=1)
					D[0]=NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion0);
				if(NIonsTypes>=2)
					D[1]=NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion1);
				if(NIonsTypes>=3)
					D[2]=NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion2);
				if(NIonsTypes>=4)
					D[3]=NIndexing->GetCMap(NodeIndexing::DiffConst,NodeIndexing::Ion3);
				
				WriteDiffusion(Dt->DiffusionMapFile.c_str());
			}
		}
	}
	if(Dt->PMFMapFile!="")
	{
		WritePMF(Dt->PMFMapFile.c_str());
	}
	pnpPrint0("Writting done\n");
	return EXIT_SUCCESS;
}
int ContWorld::SaveQstPhi(const char * filename)
{
	PNP_EXIT_FAIL_NULL(NIndexing,"NodeIndexing not initialize\n");
	PNP_EXIT_FAIL_NULL(Potential,"Potential not initialize\n");
	
	pnpPrint0("Save Qst and Phi to %s\n",filename);
	
	int i;
	float fpoh= 4*M_PI*GridScale;
	
	float *lQst = new float[NIndexing->QNum];
	float *lPhi = new float[NIndexing->QNum];
	int *qndx=NIndexing->GetChargeIndex();
	for(i=0;i<NIndexing->QNum;i++)
	{
		lQst[i]=NIndexing->Q[i]/fpoh;
	}
	for(i=0;i<NIndexing->QNum;i++)
	{
		lPhi[i]=Potential[qndx[i]];
	}
	//lPhi
	FILE *out=fopen(filename,"wb");
	fwrite(&(NIndexing->QNum),sizeof(int),1,out);
	fwrite(qndx,sizeof(int),NIndexing->QNum,out);
	fwrite(lQst,sizeof(float),NIndexing->QNum,out);
	fwrite(lPhi,sizeof(float),NIndexing->QNum,out);
	fclose(out);
	DeleteCArray(lQst);
	DeleteCArray(lPhi);
	DeleteCArray(qndx);
	return EXIT_SUCCESS;
}
int ContWorld::BorderExchange(float *ArrayToExchange)
{
	int i,j,k;
	int itmp1,itmp2;
	int gridSizeXY=GridSize[0]*GridSize[1];
#ifdef MPI_PARALLEL
	MPI::Status	status;
#endif
	int tag;
	if(PBC[0])
		for(k=0;k<GridSize[2];k++)
			for(j=0;j<GridSize[1];j++)
	{
		itmp1=GridSize[0] * j + gridSizeXY * k;
		itmp2=itmp1 + GridSize[0] - 2;
		ArrayToExchange[itmp1]=ArrayToExchange[itmp2];
		ArrayToExchange[itmp2+1]=ArrayToExchange[itmp1+1];
	}
	if(PBC[1])
		for(k=0;k<GridSize[2];k++)
			for(i=0;i<GridSize[0];i++)
	{
		itmp1=i + gridSizeXY *k;
		itmp2=itmp1+gridSizeXY - GridSize[0];
		ArrayToExchange[itmp1]=ArrayToExchange[itmp2 - GridSize[0]];
		ArrayToExchange[itmp2]=ArrayToExchange[itmp1 + GridSize[0]];
	}
	if(NProcs==1)
	{
		if(PBC[2]){
			itmp1=gridSizeXY*(GridSize[2]-2);
			itmp2=itmp1+gridSizeXY;
			for(i=0;i<gridSizeXY;i++){
				ArrayToExchange[i]=ArrayToExchange[i+itmp1];
				ArrayToExchange[i+itmp2]=ArrayToExchange[i+gridSizeXY];
			}
		}
	}
	else
	{
#ifdef MPI_PARALLEL
		if(MyRank%2==0)
		{
			if(MyRank<NProcs-1)
			{
				tag=DO_ONE;
				pnpsapp->MyComGroup.Send(ArrayToExchange+(GridSize[2]-2)*gridSizeXY, gridSizeXY, MPI::FLOAT, MyRank+1, tag);
				tag=DO_TWO;
				pnpsapp->MyComGroup.Recv(ArrayToExchange+(GridSize[2]-1)*gridSizeXY, gridSizeXY, MPI::FLOAT, MyRank+1, tag, status);
			}
			if(MyRank!=0)
			{
				tag=DO_THREE;
				pnpsapp->MyComGroup.Recv(ArrayToExchange, gridSizeXY, MPI::FLOAT, MyRank-1, tag, status);
				tag=DO_FOUR;
				pnpsapp->MyComGroup.Send(ArrayToExchange+gridSizeXY, gridSizeXY, MPI::FLOAT, MyRank-1, tag);
			}
		}
		else
		{
			tag=DO_ONE;
			pnpsapp->MyComGroup.Recv(ArrayToExchange, gridSizeXY, MPI::FLOAT, MyRank-1, tag, status);
			tag=DO_TWO;
			pnpsapp->MyComGroup.Send(ArrayToExchange+gridSizeXY, gridSizeXY, MPI::FLOAT, MyRank-1, tag);
			if(MyRank<NProcs-1)
			{
				tag=DO_THREE;
				pnpsapp->MyComGroup.Send(ArrayToExchange+(GridSize[2]-2)*gridSizeXY, gridSizeXY, MPI::FLOAT, MyRank+1, tag);
				tag=DO_FOUR;
				pnpsapp->MyComGroup.Recv(ArrayToExchange+(GridSize[2]-1)*gridSizeXY, gridSizeXY, MPI::FLOAT, MyRank+1, tag, status);
			}
		}
			//exchange by PBC over z axis
		if(PBC[2])
		{
			if(MyRank==0)
			{
				tag=DO_FIVE;
				pnpsapp->MyComGroup.Send(ArrayToExchange+gridSizeXY, gridSizeXY, MPI::FLOAT, NProcs-1, tag);
				tag=DO_SIX;
				pnpsapp->MyComGroup.Recv(ArrayToExchange, gridSizeXY, MPI::FLOAT, NProcs-1, tag, status);
			}
			else if(MyRank==NProcs-1)
			{
				tag=DO_FIVE;
				pnpsapp->MyComGroup.Recv(ArrayToExchange+(GridSize[2]-1)*gridSizeXY, gridSizeXY, MPI::FLOAT, 0, tag, status);
				tag=DO_SIX;
				pnpsapp->MyComGroup.Send(ArrayToExchange+(GridSize[2]-2)*gridSizeXY, gridSizeXY, MPI::FLOAT, 0, tag);
			}
		}
#endif
	}
	return EXIT_SUCCESS;
}
int ContWorld::BorderExchangeInt(int *ArrayToExchange)
{
	int i,j,k;
	int itmp1,itmp2;
	int gridSizeXY=GridSize[0]*GridSize[1];
#ifdef MPI_PARALLEL
	MPI::Status	status;
#endif
	int tag;
	if(PBC[0])
		for(k=0;k<GridSize[2];k++)
			for(j=0;j<GridSize[1];j++)
	{
		itmp1=GridSize[0] * j + gridSizeXY * k;
		itmp2=itmp1 + GridSize[0] - 2;
		ArrayToExchange[itmp1]=ArrayToExchange[itmp2];
		ArrayToExchange[itmp2+1]=ArrayToExchange[itmp1+1];
	}
	if(PBC[1])
		for(k=0;k<GridSize[2];k++)
			for(i=0;i<GridSize[0];i++)
	{
		itmp1=i + gridSizeXY *k;
		itmp2=itmp1+gridSizeXY - GridSize[0];
		ArrayToExchange[itmp1]=ArrayToExchange[itmp2 - GridSize[0]];
		ArrayToExchange[itmp2]=ArrayToExchange[itmp1 + GridSize[0]];
	}
	if(NProcs==1)
	{
		if(PBC[2]){
			itmp1=gridSizeXY*(GridSize[2]-2);
			itmp2=itmp1+gridSizeXY;
			for(i=0;i<gridSizeXY;i++){
				ArrayToExchange[i]=ArrayToExchange[i+itmp1];
				ArrayToExchange[i+itmp2]=ArrayToExchange[i+gridSizeXY];
			}
		}
	}
	return EXIT_SUCCESS;
}
int ContWorld::BorderExchangeDouble(double *ArrayToExchange)
{
	int i,j,k;
	int itmp1,itmp2;
	int gridSizeXY=GridSize[0]*GridSize[1];
	int tag;
	if(PBC[0])
		for(k=0;k<GridSize[2];k++)
			for(j=0;j<GridSize[1];j++)
	{
		itmp1=GridSize[0] * j + gridSizeXY * k;
		itmp2=itmp1 + GridSize[0] - 2;
		ArrayToExchange[itmp1]=ArrayToExchange[itmp2];
		ArrayToExchange[itmp2+1]=ArrayToExchange[itmp1+1];
	}
	if(PBC[1])
		for(k=0;k<GridSize[2];k++)
			for(i=0;i<GridSize[0];i++)
	{
		itmp1=i + gridSizeXY *k;
		itmp2=itmp1+gridSizeXY - GridSize[0];
		ArrayToExchange[itmp1]=ArrayToExchange[itmp2 - GridSize[0]];
		ArrayToExchange[itmp2]=ArrayToExchange[itmp1 + GridSize[0]];
	}
	if(NProcs==1)
	{
		if(PBC[2]){
			itmp1=gridSizeXY*(GridSize[2]-2);
			itmp2=itmp1+gridSizeXY;
			for(i=0;i<gridSizeXY;i++){
				ArrayToExchange[i]=ArrayToExchange[i+itmp1];
				ArrayToExchange[i+itmp2]=ArrayToExchange[i+gridSizeXY];
			}
		}
	}
	else{
//to be implemented
	}
	return EXIT_SUCCESS;
}
int ContWorld::SetInitConcentrationFromNIndexing()
{
	int i;
	if(C==NULL)
	{
		C=new float*[NIonsTypes];
		for(i=0;i<NIonsTypes;i++)
		{
			C[i]=NULL;
		}
	}
	for(i=0;i<NIonsTypes;i++)
	{
		if(C[i]==NULL)C[i]=NIndexing->GetCMap(NodeIndexing::Conc,NIndexing->IonField[i]);
	}
	
	return EXIT_SUCCESS;
}
int ContWorld::convertRealToOriginalGrid(float *mapIn, float *mapOut){
	/*!
	Convert Arrays from calculative box (GridSize[]) to original box (GridSizeOriginal[]). The difference appears because of periodic boundary condition. mapIn and mapOut can be the same array, after rearragement just mapOut have sence.
	 */
	int i,j,k,gridPoint1,gridPoint2;
	int itmp1=GridSizeOriginal[0];
	int itmp2=itmp1*GridSizeOriginal[1];
	int itmp3=GridSize[0];
	int itmp4=itmp3*GridSize[1];
	int itmp6,itmp5;
	
	//If periodic Boundary condition then rearrange map
	if(GS_XYZ!=GridSizeXYZOriginal){		
		//Move data
		for(k=0;k<GridSizeOriginal[2];k++)
			for(j=0;j<GridSizeOriginal[1];j++)
				for(i=0;i<GridSizeOriginal[0];i++)
		{
			gridPoint1 = i + itmp1 * j + itmp2 *k;
			gridPoint2 = i + itmp3 * j + itmp4 *k;
			if(PBC[0])gridPoint2+=1;
			if(PBC[1])gridPoint2+=itmp3;
			if(PBC[2])gridPoint2+=itmp4;
			mapOut[gridPoint1]=mapIn[gridPoint2];
		}
	}
	else if(mapIn!=mapOut){
		for(i=0;i<GridSizeXYZOriginal;i++)mapOut[i]=mapIn[i];
	}
	return EXIT_SUCCESS;
}
int ContWorld::convertOriginalToRealGrid(float *mapIn, float *mapOut)
{
	/*!
	Convert Arrays from original box (GridSizeOriginal[]) to calculative box (GridSize[]). The difference appears because of periodic boundary condition. mapIn and mapOut can be the same array, after rearragement just mapOut have sence.
	 */
	int i,j,k,gridPoint1,gridPoint2;
	int itmp1=GridSizeOriginal[0];
	int itmp2=itmp1*GridSizeOriginal[1];
	int itmp3=GridSizeGlobal[0];
	int itmp4=itmp3*GridSizeGlobal[1];
	int itmp6,itmp5;
	//If periodic Boundary condition then rearrange map
	if(GridSizeXYZGlobal!=GridSizeXYZOriginal){
		//Move data
		for(k=GridSizeOriginal[2]-1;k>-1;k--)
			for(j=GridSizeOriginal[1]-1;j>-1;j--)
				for(i=GridSizeOriginal[0]-1;i>-1;i--)
		{
			gridPoint1 = i + itmp1 * j + itmp2 *k;
			gridPoint2 = i + itmp3 * j + itmp4 *k;
			if(PBC[0])gridPoint2+=1;
			if(PBC[1])gridPoint2+=itmp3;
			if(PBC[2])gridPoint2+=itmp4;
			mapOut[gridPoint2]=mapIn[gridPoint1];
		}
		//copy sides;
		if(PBC[0])
			for(k=0;k<GridSizeGlobal[2];k++)
				for(j=0;j<GridSizeGlobal[1];j++)
		{
			itmp1=itmp3 * j + itmp4 * k;
			itmp2=itmp1 + itmp3 - 2;
			mapOut[itmp1]=mapOut[itmp2];
			mapOut[itmp2+1]=mapOut[itmp1+1];
		}
		if(PBC[1])
			for(k=0;k<GridSizeGlobal[2];k++)
				for(i=0;i<GridSizeGlobal[0];i++)
		{
			itmp1=i + itmp4 *k;
			itmp2=itmp1+itmp4 - itmp3;
			mapOut[itmp1]=mapOut[itmp2 - itmp3];
			mapOut[itmp2]=mapOut[itmp1 + itmp3];
		}
		if(PBC[2]){
			itmp1=itmp4*(GridSizeGlobal[2]-2);
			itmp2=itmp1+itmp4;
			for(i=0;i<itmp4;i++){
				mapOut[i]=mapOut[i+itmp1];
				mapOut[i+itmp2]=mapOut[i+itmp4];
			}
		}
	}
	else if(mapIn!=mapOut){
		for(i=0;i<GridSizeXYZOriginal;i++)mapOut[i]=mapIn[i];
	}
	return EXIT_SUCCESS;
}
int ContWorld::convertOriginalToRealGrid(int *mapIn, int *mapOut)
{
	/*!
	Convert Arrays from original box (GridSizeOriginal[]) to calculative box (GridSize[]). The difference appears because of periodic boundary condition. mapIn and mapOut can be the same array, after rearragement just mapOut have sence.
	 */
	int i,j,k,gridPoint1,gridPoint2;
	int itmp1=GridSizeOriginal[0];
	int itmp2=itmp1*GridSizeOriginal[1];
	int itmp3=GridSizeGlobal[0];
	int itmp4=itmp3*GridSizeGlobal[1];
	int itmp6,itmp5;
	//If periodic Boundary condition then rearrange map
	if(GridSizeXYZGlobal!=GridSizeXYZOriginal){
		//Move data
		for(k=GridSizeOriginal[2]-1;k>-1;k--)
			for(j=GridSizeOriginal[1]-1;j>-1;j--)
				for(i=GridSizeOriginal[0]-1;i>-1;i--)
		{
			gridPoint1 = i + itmp1 * j + itmp2 *k;
			gridPoint2 = i + itmp3 * j + itmp4 *k;
			if(PBC[0])gridPoint2+=1;
			if(PBC[1])gridPoint2+=itmp3;
			if(PBC[2])gridPoint2+=itmp4;
			mapOut[gridPoint2]=mapIn[gridPoint1];
		}
		//copy sides;
		if(PBC[0])
			for(k=0;k<GridSizeGlobal[2];k++)
				for(j=0;j<GridSizeGlobal[1];j++)
		{
			itmp1=itmp3 * j + itmp4 * k;
			itmp2=itmp1 + itmp3 - 2;
			mapOut[itmp1]=mapOut[itmp2];
			mapOut[itmp2+1]=mapOut[itmp1+1];
		}
		if(PBC[1])
			for(k=0;k<GridSizeGlobal[2];k++)
				for(i=0;i<GridSizeGlobal[0];i++)
		{
			itmp1=i + itmp4 *k;
			itmp2=itmp1+itmp4 - itmp3;
			mapOut[itmp1]=mapOut[itmp2 - itmp3];
			mapOut[itmp2]=mapOut[itmp1 + itmp3];
		}
		if(PBC[2]){
			itmp1=itmp4*(GridSizeGlobal[2]-2);
			itmp2=itmp1+itmp4;
			for(i=0;i<itmp4;i++){
				mapOut[i]=mapOut[i+itmp1];
				mapOut[i+itmp2]=mapOut[i+itmp4];
			}
		}
	}
	else if(mapIn!=mapOut){
		for(i=0;i<GridSizeXYZOriginal;i++)mapOut[i]=mapIn[i];
	}
	return EXIT_SUCCESS;
}
int ContWorld::GetBorder(int *z1,int *z2,int ProcRank)
{
	int i,i1,i2;
	
	*z1=0;
	*z2=GridSizeGlobal[2]-1;
#ifdef MPI_PARALLEL
	if(MPISpreadWorld==0)
	{
		i=0;
		i1=i*(GridSizeGlobal[2]/NProcs);
		i2=(i+1)*(GridSizeGlobal[2]/NProcs)-1;
		if(i!=0)i1--;
		if(i!=NProcs-1)i2++;
		if((i2-i1+1)%2==0)i2++;
		if(i==NProcs-1)i2=GridSizeGlobal[2]-1;
		
		for(i=1;i<=ProcRank;i++)
		{
			i1=i2-1;		
			i2=(i+1)*(GridSizeGlobal[2]/NProcs)-1;
			if(i!=NProcs-1)i2++;
			if((i2-i1+1)%2==0)i2++;
			if(i==NProcs-1)i2=GridSizeGlobal[2]-1;
		}
		*z1=i1;
		*z2=i2;
	}
	else
	{
		*z1=MPISpreadWorld[NProcs-1][ProcRank*2+1];
		*z2=MPISpreadWorld[NProcs-1][ProcRank*2+2];
	}
#endif
	return EXIT_SUCCESS;
}
int ContWorld::CheckArrays(const char* ArSet,bool BuildIfNI)
{
	int i;
	if(BuildIfNI)
	{
		for(i=0;i<3;i++)
		{
			if(Epsilon[i]==NULL)
			{
				Epsilon[i]=NIndexing->GetCMap(NodeIndexing::DielConst,NIndexing->EpsilonField[i]);
			}
		}
		if(Qstat==NULL)
		{
			Qstat=NIndexing->GetChargeArray();
		}
		if(Potential == NULL){
			if(!(Potential = new float[GS_XYZ])){
				fprintf(stderr,"ERROR 104: No memory available\n");
				exit(104);
			}
			for(i=0;i<GS_XYZ;i++)Potential[i]=0.0;
		}
	}
	return EXIT_SUCCESS;
}
int ContWorld::CheckSystem()
{
	int i,j,k;
	int IonType,GrdPnt;
	pnpPrint0("Will check the system.\n");
	if(PMF!=NULL)
	{
		pnpPrint0("* PMF present, will check values:\n");
		for(IonType=0;IonType<NIonsTypes;IonType++)
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
						dPMFx=PMF[IonType][GrdPnt]-PMF[IonType][GrdPnt-1];
						dPMFy=PMF[IonType][GrdPnt]-PMF[IonType][GrdPnt-GS_X];
						dPMFz=PMF[IonType][GrdPnt]-PMF[IonType][GrdPnt-GS_XY];
						
						if((dPMFx*dPMFx)>4.0)countPMFx++;
						if((dPMFy*dPMFy)>4.0)countPMFy++;
						if((dPMFz*dPMFz)>4.0)countPMFz++;
					}
				}
			}
			pnpPrint0("\tIon %d: dPMF>2kT countPMFx=%d,countPMFy=%d,countPMFz=%d Throught All Space\n"
					,IonType,countPMFx,countPMFy,countPMFz);
		}
		{
			//Chech dPMF>2kT At Plases where D!=0
			int countBiggerThen2kT=0;
			VectorIntField3D *iV=new VectorIntField3D(GridSize,GridScale,NIonsTypes);
			iV->FillValue(0);
			for(IonType=0;IonType<NIonsTypes;IonType++)
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
							if(NIndexing->GetDiffFloat(IonType,GrdPnt)!=0.0)
							{
								dPMFx=PMF[IonType][GrdPnt]-PMF[IonType][GrdPnt-1];
								dPMFy=PMF[IonType][GrdPnt]-PMF[IonType][GrdPnt-GS_X];
								dPMFz=PMF[IonType][GrdPnt]-PMF[IonType][GrdPnt-GS_XY];
								
								if((dPMFx*dPMFx)>4.0&&NIndexing->GetDiffFloat(IonType,GrdPnt-1)!=0.0)
								{
									countPMFx++;
									iV->V[IonType][GrdPnt]+=1;
								}
								if((dPMFy*dPMFy)>4.0&&NIndexing->GetDiffFloat(IonType,GrdPnt-GS_X)!=0.0)
								{
									countPMFy++;
									iV->V[IonType][GrdPnt]+=10;
								}
								if((dPMFz*dPMFz)>4.0&&NIndexing->GetDiffFloat(IonType,GrdPnt-GS_XY)!=0.0)
								{
									countPMFz++;
									iV->V[IonType][GrdPnt]+=100;
								}
							}
						}
					}
				}
				pnpPrint0("\tIon %d: dPMF>2kT countPMFx=%d,countPMFy=%d,countPMFz=%d Throught D!=0\n"
						,IonType,countPMFx,countPMFy,countPMFz);
				countBiggerThen2kT+=countPMFx+countPMFy+countPMFz;
			}
			if(countBiggerThen2kT>0)
			{
				pnpPrint0("\tSince you have places there |dPMF_r|>2 kT will write this place to PNP_SYSCHECK_PMF_bigger_2kT_Mask.g\n");
				iV->WriteToFile("PNP_SYSCHECK_PMF_bigger_2kT_Mask.gz",2);
			}
			delete iV;
		}
		pnpPrint0("= PMF checked\n");
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
	return EXIT_SUCCESS;
}
float ContWorld::AdditionalPotantialAlongZ(float Z0,float Z1,float Phi0,float Phi1,int InternalZ)
{
	int GlobalZ;
	int irc[3],z0,z1,j;
	float DPhi=Phi1-Phi0;
	
	GlobalZ=MyGlobalZ0+InternalZ;
	for(j=0;j<3;j++)irc[j]=GridSizeGlobal[j]/2;
	Z0*=GridScale;
	Z1*=GridScale;
	Z0+=(float)((int)GridSizeGlobal[2]/2);
	Z1+=(float)((int)GridSizeGlobal[2]/2);
	z0=(int)(Z0);
	z1=(int)(Z1+1.0);
	if(GlobalZ<=z0)return Phi0;
	else if(GlobalZ>=z1)return Phi1;
	else
	{
		return Phi0+DPhi*(float)(GlobalZ-z0)/(float)(z1-z0);
	}
}
int ContWorld::AddPotential(float Z0,float Z1,float Phi0,float Phi1)
{
	DbgPrint0("AddPotential Z=[%g,%g] Phi=[%g,%g]\n",Z0,Z1,Phi0,Phi1);
	int i,j,k,GridPoint;
	float Phi;
	//from mV to kT
	Phi0/=CONFAC;
	Phi1/=CONFAC;
  
	if(Potential==NULL)
	{
		Potential=new float[GS_XYZ];
		for(i=0;i<GS_XYZ;i++)
		{
			Potential[i]=0.0;
		}
	}
  
	for(k=0;k<GridSize[2];k++)
	{
		Phi=AdditionalPotantialAlongZ(Z0, Z1, Phi0, Phi1,k);
		for(i=0;i<GridSize[0];i++)
			for(j=0;j<GridSize[1];j++)
		{
			GridPoint=i+j*GS_X+k*GS_XY;
			Potential[GridPoint]+=Phi;
		}
	}
	return EXIT_SUCCESS;
}
int ContWorld::ApplyAsymConcByScaling(float z,float* Cz0,float* Cz1)
{
	int ion,i,j,k,GridPoint;
	pnpPrint0("ApplyAsymConc Z=%f Cz0=[",z);
	for(i=0;i<NIonsTypes;i++)
		pnpPrint0(" %f",Cz0[i]);
	pnpPrint0("] Cz1=[");
	for(i=0;i<NIonsTypes;i++)
		pnpPrint0(" %f",Cz1[i]);
	pnpPrint0("]\n");
	
	float fZ=ConvFloatToLocIntUnitsZ(z);
	float *fCz0=new float[NIonsTypes];
	float *fCz1=new float[NIonsTypes];
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS/(GridScale*GridScale*GridScale);
	for(i=0;i<NIonsTypes;i++)
	{
		fCz0[i]=coef*Cz0[i];
		fCz1[i]=coef*Cz1[i];
	}
	DbgPrint0("Internal Units:ApplyAsymConc Z=%f Cz0=[",fZ);
	for(i=0;i<NIonsTypes;i++)
		DbgPrint0(" %f ",fCz0[i]);
	DbgPrint0("] Cz1=[");
	for(i=0;i<NIonsTypes;i++)
		DbgPrint0(" %f ",fCz1[i]);
	DbgPrint0("]\n");
	
	if(C==NULL)
	{
		SetInitConcentrationFromNIndexing();
	}
	float Cold=C[0][0];
	if(D==NULL)
	{
		for(ion=0;ion<NIonsTypes;ion++)
		{
			for(k=0;k<GridSize[2];k++)
			{
				float tmpC;
				if(k<=fZ)
					tmpC=fCz0[ion];
				else
					tmpC=fCz1[ion];
				
				for(i=0;i<GridSize[0];i++)
					for(j=0;j<GridSize[1];j++)
				{
					GridPoint=i+j*GS_X+k*GS_XY;
					if(NIndexing->GetDiffFloat(ion,GridPoint)>0.0)
						C[ion][GridPoint]*=tmpC/Cold;
					else
						C[ion][GridPoint]=0.0;
				}
			}
		}
	}
	else
	{
		for(ion=0;ion<NIonsTypes;ion++)
		{
			for(k=0;k<GridSize[2];k++)
			{
				float tmpC;
				if(k<=fZ)
					tmpC=fCz0[ion];
				else
					tmpC=fCz1[ion];
				
				for(i=0;i<GridSize[0];i++)
					for(j=0;j<GridSize[1];j++)
				{
					GridPoint=i+j*GS_X+k*GS_XY;
					if(D[ion][GridPoint]>0.0)
						C[ion][GridPoint]*=tmpC/Cold;
					else
						C[ion][GridPoint]=0.0;
				}
			}
		}
	}	
	delete [] fCz0;
	delete [] fCz1;
	
	return EXIT_SUCCESS;
}
int ContWorld::ApplyAsymConc(float z,float* Cz0,float* Cz1)
{
	int ion,i,j,k,GridPoint;
	pnpPrint0("ApplyAsymConc Z=%f Cz0=[",z);
	for(i=0;i<NIonsTypes;i++)
		pnpPrint0(" %f",Cz0[i]);
	pnpPrint0("] Cz1=[");
	for(i=0;i<NIonsTypes;i++)
		pnpPrint0(" %f",Cz1[i]);
	pnpPrint0("]\n");
	
	float fZ=ConvFloatToLocIntUnitsZ(z);
	float *fCz0=new float[NIonsTypes];
	float *fCz1=new float[NIonsTypes];
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS/(GridScale*GridScale*GridScale);
	for(i=0;i<NIonsTypes;i++)
	{
		fCz0[i]=coef*Cz0[i];
		fCz1[i]=coef*Cz1[i];
	}
	DbgPrint0("Internal Units:ApplyAsymConc Z=%f Cz0=[",fZ);
	for(i=0;i<NIonsTypes;i++)
		DbgPrint0(" %f ",fCz0[i]);
	DbgPrint0("] Cz1=[");
	for(i=0;i<NIonsTypes;i++)
		DbgPrint0(" %f ",fCz1[i]);
	DbgPrint0("]\n");
	
	if(C==NULL)
	{
		C=new float*[NIonsTypes];
		for(i=0;i<NIonsTypes;i++)
			C[i]=NULL;
	}
	for(i=0;i<NIonsTypes;i++)
	{
		if(C[i]==NULL)
		{
			C[i]=new float[GS_XYZ];
			for(j=0;j<GS_XYZ;j++)
				C[i][j]=0.0;
		}
	}
	if(D==NULL)
	{
		for(ion=0;ion<NIonsTypes;ion++)
		{
			for(k=0;k<GridSize[2];k++)
			{
				float tmpC;
				if(k<=fZ)
					tmpC=fCz0[ion];
				else
					tmpC=fCz1[ion];
				
				for(i=0;i<GridSize[0];i++)
					for(j=0;j<GridSize[1];j++)
				{
					GridPoint=i+j*GS_X+k*GS_XY;
					if(NIndexing->GetDiffFloat(ion,GridPoint)>0.0)
						C[ion][GridPoint]=tmpC;
					else
						C[ion][GridPoint]=0.0;
				}
			}
		}
	}
	else
	{
		for(ion=0;ion<NIonsTypes;ion++)
		{
			for(k=0;k<GridSize[2];k++)
			{
				float tmpC;
				if(k<=fZ)
					tmpC=fCz0[ion];
				else
					tmpC=fCz1[ion];
				
				for(i=0;i<GridSize[0];i++)
					for(j=0;j<GridSize[1];j++)
				{
					GridPoint=i+j*GS_X+k*GS_XY;
					if(D[ion][GridPoint]>0.0)
						C[ion][GridPoint]=tmpC;
					else
						C[ion][GridPoint]=0.0;
				}
			}
		}
	}	
	delete [] fCz0;
	delete [] fCz1;
	
	return EXIT_SUCCESS;
}
int ContWorld::ConvertIonStrengthToDynamicCharge()
{
	int i;	
	float ch1=0.0,ch2=0.0;
	float fpoh= 4*M_PI*GridScale;
	float IStrength;
	float temp=fpoh*COANGS/(GridScale*GridScale*GridScale);
	fprintf(stdout,"\nConvert IonStrength To Dynamic Charge\n");
	if(C==NULL)
	{
		C=new float*[2];
		C[0]=NULL;
		C[1]=NULL;
	}
	if(C[0]==NULL)
	{
		if(!(C[0] = new float[GS_XYZ])){
			fprintf(stderr,"ERROR 104: No memory available\n");
			exit(104);
		}
	}
	if(C[1]==NULL)
	{
		if(!(C[1] = new float[GS_XYZ])){
			fprintf(stderr,"ERROR 104: No memory available\n");
			exit(104);
		}
	}
	for(i=0;i<NodeIndexMaxValues;i++)
	{
		if(NIndexing->C[i]!=0.0)
		{
			IStrength=NIndexing->C[i];
			break;
		}
	}
	pnpPrint("Will use bulk concentration: %f\n",IStrength/temp);
	
	for(i=0;i<GS_XYZ;i++)
	{
		if(NIndexing->GetConcFloat(0,i)>0.0)
		{
			C[0][i]=IStrength*exp(-Potential[i]);		
			C[1][i]=IStrength*exp(Potential[i]);
			ch1+=C[0][i];
			ch2+=C[1][i];
		}
		else
		{
			C[0][i]=0.0;
			C[1][i]=0.0;
		}
	}
	fprintf(stdout,"Charge Is %g	+ (- %g )= %g\n",ch1/fpoh,ch2/fpoh,(ch1-ch2)/fpoh);
	return EXIT_SUCCESS;
}
int ContWorld::RemoveDiffusionPoints(int IType,int* pnts,int n)
{
	int i,k,j,GrdPnt;
	pnpPrint("<ContWorld::RemoveDiffusionPoints ion=\"%d\" pnts=\"",IType);
	for(i=0;i<n;i++)
	{
		pnpPrint("%d",pnts[i]);
		if(i<n-1)pnpPrint(" ");
	}
	pnpPrint("\">\n");
	pnpPrint("will set to zero diffusion for ion %d at points:\n",IType);
	for(i=0;i<n;i++)
	{
		pnpPrint("\tD(pnt=%d)=%f\n",pnts[i],NIndexing->GetDiffFloat(IType,pnts[i]));
	}
	PNP_EXIT_FAIL_NULL(NIndexing,"NIndexing is not initialized\n");
	
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS / (GridScale*GridScale*GridScale);
	
	
	int DiffZeroInd=NIndexing->GetDiffZeroInd();
	int count=0;
	for(i=0;i<n;i++)
	{
		GrdPnt=pnts[i];
		if(NIndexing->GetDiffFloat(IType,GrdPnt)!=0.0)
		{
			NIndexing->SetIonField(IType,GrdPnt,DiffZeroInd);
			count++;
		}
		else
		{
			pnpPrint("Diffusion is already zero for ion=%d pnt=%d D=%f\n",IType,GrdPnt,NIndexing->GetDiffFloat(IType,GrdPnt));
		}
	}
	
	pnpPrint("%d points was removed for ion %d\n",count,IType);
	
	if(count>0)
	{
		NIndexing->RemoveBadDiffusionPoints();
	}
	
	
	
	pnpPrint("New values are for diffusion for ion %d:\n",IType);
	for(i=0;i<n;i++)
	{
		pnpPrint("\tD(pnt=%d)=%f\n",pnts[i],NIndexing->GetDiffFloat(IType,pnts[i]));
	}
	pnpPrint("</ContWorld::RemoveDiffusionPoints>\n");
	return EXIT_SUCCESS;
}
int ContWorld::RemoveDiffusionPointsAtNegativeC()
{
	int IType,i,k,j,GrdPnt;
	PNP_EXIT_FAIL_NULL(NIndexing,"NIndexing is not initialized\n");
	PNP_EXIT_FAIL_NULL(C,"concentration is not initialize at m_ContWorld\n");
	for(IType=0;IType<NIonsTypes;IType++)
		PNP_EXIT_FAIL_NULL(C[IType],"concentration is not initialize at m_ContWorld\n");
	
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS / (GridScale*GridScale*GridScale);
	
	
	int DiffZeroInd=NIndexing->GetDiffZeroInd();
	int count=0;
	for(IType=0;IType<NIonsTypes;IType++)
	{
		for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
		{
			if(NIndexing->GetDiffFloat(IType,GrdPnt)!=0.0)
			{
				if(C[IType][GrdPnt]<0.0)
				{
					NIndexing->SetIonField(IType,GrdPnt,DiffZeroInd);
					count++;
				}
			}
		}
	}
	if(count>0)
	{
		NIndexing->CalcDiffBoarder();
		NIndexing->RemoveBadDiffusionPoints();
	}
	
	pnpPrint("%d points was remove because of its concentration was < 0.0\n",count);
	return EXIT_SUCCESS;
}
bool ContWorld::FindImplicitMembranePosition(double *z0,double *z1)
{
	bool FloatEpsilon=false;
	bool IntEpsilon=false;
	if(this->Epsilon[0]!=NULL && this->Epsilon[1]!=NULL && this->Epsilon[2]!=NULL)
		FloatEpsilon=true;
	if(this->NIndexing!=NULL)
		IntEpsilon=true;

	if(FloatEpsilon)
	{
		pnpError("FindImplicitMembranePosition is not implemented for variable epsilon(ContWorld->Epsilon[])\n");
		return false;
	}
	else if(IntEpsilon)
	{
		int iz;
		int GS_X=this->NIndexing->GridSize[0];
		int GS_Y=this->NIndexing->GridSize[0];
		int GS_Z=this->NIndexing->GridSize[0];
		int GS_XY=GS_X*GS_Y;
		int GS_XYZ=GS_X*GS_Y*GS_Z;

		iz=0;
		int iEps0=this->NIndexing->GetDiel(0,GS_XY*iz);
		iz=GS_Z-1;
		int iEps1=this->NIndexing->GetDiel(0,GS_XY*iz);

		if(iEps0!=iEps1)
		{
			pnpError("FindImplicitMembranePosition: dielectric constant is difrent in cis/trans compartments\n");
			return false;
		}

		//find bottom z
		float fz0;
		for(iz=0;iz<2*GS_Z;iz++)
		{
			//iz%2==0 epsilon which shifted in y
			//iz%2==1 epsilon which shifted in z
			int iEps=this->NIndexing->GetDiel(1+iz%2,GS_XY*(iz/2));
			if(iEps!=iEps0)
				fz0=iz*0.5;
		}
		//find top z
		float fz1;
		for(iz=2*GS_Z-1;iz>=0;iz--)
		{
			//iz%2==0 epsilon which shifted in y
			//iz%2==1 epsilon which shifted in z
			int iEps=this->NIndexing->GetDiel(1+iz%2,GS_XY*(iz/2));
			if(iEps!=iEps1)
				fz1=iz*0.5;
		}
		//printf("FindImplicitMembranePosition %g %g\n",z0,z1);

		*z0=ConvGlobIntToGlobExtUnitsZ(fz0);
		*z1=ConvGlobIntToGlobExtUnitsZ(fz1);

		printf("Found implicit membrane at z %g %g\n",*z0,*z1);
		return true;
	}

	return false;
}
ContWorld* ReadContWorldFromNodeIndexing(const char * node_indexing_filename)
{
	ContWorld* cw=new ContWorld();
	cw->ReadNodeIndexing(node_indexing_filename);
	return cw;
}

#ifdef HARLEM_MOD
HaVec_float* ContWorld::GetNumberOfIonsAlongChannel(int ion, float Xcenter, float Ycenter,HaVec_float* RlimAlZ)
{
	PNP_EXIT_NULL_ON_NULL(C,"Concentration is not initialized\n");
	PNP_EXIT_NULL_ON_NULL(C[ion],"Concentration is not initialized\n");
	int ix,iy,iz,GrdPnt;
	int Rcenter[3]={GS_X/2,GS_Y/2,GS_Z/2};
	
	float xC=Xcenter*GridScale+Rcenter[0];
	float yC=Ycenter*GridScale+Rcenter[1];
	float r;
	HaVec_float* IonsAlongZ=new HaVec_float(GS_Z,0.0f);
	for(iz=0;iz<GS_Z;iz++)
		for(ix=0;ix<GS_X;ix++)
			for(iy=0;iy<GS_Y;iy++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		r=(float(ix)-xC)*(float(ix)-xC)+(float(iy)-yC)*(float(iy)-yC);
		if(r<=(RlimAlZ->GetVal_idx0(iz))*(RlimAlZ->GetVal_idx0(iz))*GridScale*GridScale)
			IonsAlongZ->SetVal_idx0(iz,IonsAlongZ->GetVal_idx0(iz)+C[ion][GrdPnt]);
	}
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS/(GridScale*GridScale*GridScale);
	for(iz=0;iz<GS_Z;iz++)
	{
		IonsAlongZ->SetVal_idx0(iz,COANGS*IonsAlongZ->GetVal_idx0(iz)/(coef*GridScale*GridScale*GridScale));
	}
	return IonsAlongZ;
}
HaVec_float* ContWorld::GetRadiusAlongZ(int ion, float Xcenter, float Ycenter, HaVec_float* RlimAlZ)
{
	int ix,iy,iz,GrdPnt;
	int Rcenter[3]={GS_X/2,GS_Y/2,GS_Z/2};
	
	float xC=Xcenter*GridScale+Rcenter[0];
	float yC=Ycenter*GridScale+Rcenter[1];
	float r;
	HaVec_float* RAlongZ=new HaVec_float(GS_Z,0.0f);
	if(D!=NULL)
	{
		for(iz=0;iz<GS_Z;iz++)
			for(ix=0;ix<GS_X;ix++)
				for(iy=0;iy<GS_Y;iy++)
		{
			GrdPnt=ix+iy*GS_X+iz*GS_XY;
			r=(float(ix)-xC)*(float(ix)-xC)+(float(iy)-yC)*(float(iy)-yC);
			if(D[ion][GrdPnt]>0.0)
				if(r<=(RlimAlZ->GetVal_idx0(iz))*(RlimAlZ->GetVal_idx0(iz))*GridScale*GridScale)
					RAlongZ->SetVal_idx0(iz,RAlongZ->GetVal_idx0(iz)+1.0);
		}
	}
	else if(NIndexing!=NULL)
	{
		for(iz=0;iz<GS_Z;iz++)
			for(ix=0;ix<GS_X;ix++)
				for(iy=0;iy<GS_Y;iy++)
		{
			GrdPnt=ix+iy*GS_X+iz*GS_XY;
			r=(float(ix)-xC)*(float(ix)-xC)+(float(iy)-yC)*(float(iy)-yC);
			if(NIndexing->GetDiffFloat(ion,GrdPnt)>0.0)
				if(r<=(RlimAlZ->GetVal_idx0(iz))*(RlimAlZ->GetVal_idx0(iz))*GridScale*GridScale)
					RAlongZ->SetVal_idx0(iz,RAlongZ->GetVal_idx0(iz)+1.0);
		}
	}
	else
	{
		pnpError("Diffusion is not initialized\n");
		return NULL;
	}
	for(iz=0;iz<GS_Z;iz++)
	{
		RAlongZ->SetVal_idx0(iz,sqrt(RAlongZ->GetVal_idx0(iz)/M_PI)/GridScale);
	}
	return RAlongZ;
}
int ContWorld::CorrectDiffusion(int ion, float Xcenter, float Ycenter, float Z0, float Z1, HaVec_float* RlimAlZ,HaVec_float* DscaleAlZ)
{
	PNP_EXIT_NULL_ON_NULL(D,"Diffusion is not initialized\n");
	PNP_EXIT_NULL_ON_NULL(D[ion],"Diffusion is not initialized\n");
	int ix,iy,iz,GrdPnt;
	int Rcenter[3]={GS_X/2,GS_Y/2,GS_Z/2};
	
	float xC=Xcenter*GridScale+Rcenter[0];
	float yC=Ycenter*GridScale+Rcenter[1];
	float r;
	int iz0=Z0*GridScale+Rcenter[2];
	int iz1=Z1*GridScale+Rcenter[2];
	if(iz0<0)iz0=0;
	if(iz1>GS_Z-1)iz1=GS_Z-1;
	HaVec_float* IonsAlongZ=new HaVec_float(GS_Z,0.0f);
	for(iz=iz0;iz<=iz1;iz++)
		for(ix=0;ix<GS_X;ix++)
			for(iy=0;iy<GS_Y;iy++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		r=(float(ix)-xC)*(float(ix)-xC)+(float(iy)-yC)*(float(iy)-yC);
		if(r<=(RlimAlZ->GetVal_idx0(iz))*(RlimAlZ->GetVal_idx0(iz))*GridScale*GridScale)
			if(D[ion][GrdPnt]>0.0)
				D[ion][GrdPnt]=D[ion][GrdPnt]*DscaleAlZ->GetVal_idx0(iz);
	}
	return EXIT_SUCCESS;
}
HaVec_float* ContWorld::GetDavrAlongZ(int ion, float Xcenter, float Ycenter, HaVec_float* RlimAlZ)
{
	int ix,iy,iz,GrdPnt;
	int Rcenter[3]={GS_X/2,GS_Y/2,GS_Z/2};
	
	float xC=Xcenter*GridScale+Rcenter[0];
	float yC=Ycenter*GridScale+Rcenter[1];
	float r;
	int count;
	double sum;
	HaVec_float* AlongZ=new HaVec_float(GS_Z,0.0f);
	if(D!=NULL)
	{
		for(iz=0;iz<GS_Z;iz++)
		{
			count=0;
			sum=0.0;
			for(ix=0;ix<GS_X;ix++)
				for(iy=0;iy<GS_Y;iy++)
			{
				GrdPnt=ix+iy*GS_X+iz*GS_XY;
				r=(float(ix)-xC)*(float(ix)-xC)+(float(iy)-yC)*(float(iy)-yC);
				if(D[ion][GrdPnt]>0.0)
				{
					if(r<=(RlimAlZ->GetVal_idx0(iz))*(RlimAlZ->GetVal_idx0(iz))*GridScale*GridScale)
					{
						sum+=D[ion][GrdPnt];
						count++;
					}
				}
			}
			if(count>0)
				AlongZ->SetVal_idx0(iz,sum/double(count));
			else
				AlongZ->SetVal_idx0(iz,0.0);
		}
	}
	else if(NIndexing!=NULL)
	{
		for(iz=0;iz<GS_Z;iz++)
		{
			count=0;
			sum=0.0;
			for(ix=0;ix<GS_X;ix++)
				for(iy=0;iy<GS_Y;iy++)
			{
				GrdPnt=ix+iy*GS_X+iz*GS_XY;
				r=(float(ix)-xC)*(float(ix)-xC)+(float(iy)-yC)*(float(iy)-yC);
				if(NIndexing->GetDiffFloat(ion,GrdPnt)>0.0)
				{
					if(r<=(RlimAlZ->GetVal_idx0(iz))*(RlimAlZ->GetVal_idx0(iz))*GridScale*GridScale)
					{
						sum+=NIndexing->GetDiffFloat(ion,GrdPnt);
						count++;
					}
				}
			}
			if(count>0)
				AlongZ->SetVal_idx0(iz,sum/double(count));
			else
				AlongZ->SetVal_idx0(iz,0.0);
		}
	}
	else
	{
		pnpError("Diffusion is not initialized\n");
		return NULL;
	}
	//for(iz=0;iz<GS_Z;iz++)
	//{
	//	RAlongZ->SetVal_idx0(iz,sqrt(RAlongZ->GetVal_idx0(iz)/M_PI)/GridScale);
	//}
	return AlongZ;
}
double ContWorld::GetNumberOfIonsInChannel(int ion, float Xcenter, float Ycenter, float Z0, float Z1, HaVec_float* RlimAlZ)
{
	PNP_EXIT_NULL_ON_NULL(C,"Concentration is not initialized\n");
	PNP_EXIT_NULL_ON_NULL(C[ion],"Concentration is not initialized\n");
	int ix,iy,iz,GrdPnt;
	int Rcenter[3]={GS_X/2,GS_Y/2,GS_Z/2};
	
	float xC=Xcenter*GridScale+Rcenter[0];
	float yC=Ycenter*GridScale+Rcenter[1];
	float r;
	HaVec_float* IonsAlongZ=new HaVec_float(GS_Z,0.0f);
	for(iz=0;iz<GS_Z;iz++)
		for(ix=0;ix<GS_X;ix++)
			for(iy=0;iy<GS_Y;iy++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		r=(float(ix)-xC)*(float(ix)-xC)+(float(iy)-yC)*(float(iy)-yC);
		if(r<=(RlimAlZ->GetVal_idx0(iz))*(RlimAlZ->GetVal_idx0(iz))*GridScale*GridScale)
			IonsAlongZ->SetVal_idx0(iz,IonsAlongZ->GetVal_idx0(iz)+C[ion][GrdPnt]);
	}
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS/(GridScale*GridScale*GridScale);
	for(iz=0;iz<GS_Z;iz++)
	{
		IonsAlongZ->SetVal_idx0(iz,COANGS*IonsAlongZ->GetVal_idx0(iz)/(coef*GridScale*GridScale*GridScale));
	}
	double Nions=0;
	int iz0=Z0*GridScale+Rcenter[2];
	int iz1=Z1*GridScale+Rcenter[2];
	if(iz0<0)iz0=0;
	if(iz1>GS_Z-1)iz1=GS_Z-1;
	for(iz=iz0;iz<=iz1;iz++)
	{
		Nions+=IonsAlongZ->GetVal_idx0(iz);
	}
	delete IonsAlongZ;
	return Nions;
}
HaVec_float* ContWorld::GetPotValusAlongZ(int ix,int iy)
{
	PNP_EXIT_NULL_ON_NULL(Potential,"Concentration is not initialized\n");
	int iz,GrdPnt;
	HaVec_float* VAlongZ=new HaVec_float(GS_Z,0.0f);
	for(iz=0;iz<GS_Z;iz++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		VAlongZ->SetVal_idx0(iz,Potential[GrdPnt]);
	}
	return VAlongZ;
}
HaVec_float* ContWorld::GetConcValusAlongZ(int ion,int ix,int iy)
{
	PNP_EXIT_NULL_ON_NULL(C,"Concentration is not initialized\n");
	PNP_EXIT_NULL_ON_NULL(C[ion],"Concentration is not initialized\n");
	int iz,GrdPnt;
	float fpoh= 4*M_PI*GridScale;
	float coef=fpoh*COANGS/(GridScale*GridScale*GridScale);
	HaVec_float* VAlongZ=new HaVec_float(GS_Z,0.0f);
	for(iz=0;iz<GS_Z;iz++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		VAlongZ->SetVal_idx0(iz,C[ion][GrdPnt]/coef);
	}
	return VAlongZ;
}
HaVec_float* ContWorld::GetPMFValusAlongZ(int ion,int ix,int iy)
{
	PNP_EXIT_NULL_ON_NULL(PMF,"PMF is not initialized\n");
	PNP_EXIT_NULL_ON_NULL(PMF[ion],"PMF is not initialized\n");
	int iz,GrdPnt;
	HaVec_float* VAlongZ=new HaVec_float(GS_Z,0.0f);
	for(iz=0;iz<GS_Z;iz++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		VAlongZ->SetVal_idx0(iz,PMF[ion][GrdPnt]);
	}
	return VAlongZ;
}
int ContWorld::SetPotentialToZeroWhereDiffZero(int ion)
{
	PNP_EXIT_NULL_ON_NULL(Potential,"Potential is not initialized\n");
	int ix,iy,iz,GrdPnt;
	int Rcenter[3]={GS_X/2,GS_Y/2,GS_Z/2};
	
	if(D!=NULL)
	{
		for(iz=0;iz<GS_Z;iz++)
			for(ix=0;ix<GS_X;ix++)
				for(iy=0;iy<GS_Y;iy++)
		{
			GrdPnt=ix+iy*GS_X+iz*GS_XY;
			if(D[ion][GrdPnt]==0.0)
			{
				Potential[GrdPnt]=0.0;
			}
		}
	}
	else if(NIndexing!=NULL)
	{
		for(iz=0;iz<GS_Z;iz++)
			for(ix=0;ix<GS_X;ix++)
				for(iy=0;iy<GS_Y;iy++)
		{
			GrdPnt=ix+iy*GS_X+iz*GS_XY;
			if(NIndexing->GetDiffFloat(ion,GrdPnt)==0.0)
			{
				Potential[GrdPnt]=0.0;
			}
		}
	}
	else
	{
		pnpError("Diffusion is not initialized\n");
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
HaVec_float* ContWorld::GetAvrPotAlongPore(int ion, float Xcenter, float Ycenter, HaVec_float* RlimAlZ)
{
	PNP_EXIT_NULL_ON_NULL(Potential,"Potential is not initialized\n");
	int ix,iy,iz,GrdPnt;
	int Rcenter[3]={GS_X/2,GS_Y/2,GS_Z/2};
	
	float xC=Xcenter*GridScale+Rcenter[0];
	float yC=Ycenter*GridScale+Rcenter[1];
	float r;
	int count;
	double sum;
	HaVec_float* AlongZ=new HaVec_float(GS_Z,0.0f);
	if(D!=NULL)
	{
		for(iz=0;iz<GS_Z;iz++)
		{
			count=0;
			sum=0.0;
			for(ix=0;ix<GS_X;ix++)
				for(iy=0;iy<GS_Y;iy++)
			{
				GrdPnt=ix+iy*GS_X+iz*GS_XY;
				r=(float(ix)-xC)*(float(ix)-xC)+(float(iy)-yC)*(float(iy)-yC);
				if(D[ion][GrdPnt]>0.0)
				{
					if(r<=(RlimAlZ->GetVal_idx0(iz))*(RlimAlZ->GetVal_idx0(iz))*GridScale*GridScale)
					{
						sum+=Potential[GrdPnt];
						count++;
					}
				}
			}
			if(count>0)
				AlongZ->SetVal_idx0(iz,sum/double(count));
			else
				AlongZ->SetVal_idx0(iz,0.0);
		}
	}
	else if(NIndexing!=NULL)
	{
		for(iz=0;iz<GS_Z;iz++)
		{
			count=0;
			sum=0.0;
			for(ix=0;ix<GS_X;ix++)
				for(iy=0;iy<GS_Y;iy++)
			{
				GrdPnt=ix+iy*GS_X+iz*GS_XY;
				r=(float(ix)-xC)*(float(ix)-xC)+(float(iy)-yC)*(float(iy)-yC);
				if(NIndexing->GetDiffFloat(ion,GrdPnt)>0.0)
				{
					if(r<=(RlimAlZ->GetVal_idx0(iz))*(RlimAlZ->GetVal_idx0(iz))*GridScale*GridScale)
					{
						sum+=Potential[GrdPnt];
						count++;
					}
				}
			}
			if(count>0)
				AlongZ->SetVal_idx0(iz,sum/double(count));
			else
				AlongZ->SetVal_idx0(iz,0.0);
		}
	}
	else
	{
		pnpError("Diffusion is not initialized\n");
		return NULL;
	}
	//for(iz=0;iz<GS_Z;iz++)
	//{
	//	RAlongZ->SetVal_idx0(iz,sqrt(RAlongZ->GetVal_idx0(iz)/M_PI)/GridScale);
	//}
	return AlongZ;
}
#endif
