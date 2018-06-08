//
// C++ Implementation: buildworldni
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
#  include "mpi.h"
#endif

#ifdef HARLEM_MOD
    #include "hamolset.h"
#endif

#include "buildworldni.h"

#include "pnpdebug.h"
#include "contworld.h"
#include "math.h"
#include "pnpconstants.h"
#include "mapio.h"
#include <vector>
#include <string>
#include "pnpstructs.h"

#include <iostream>
using namespace std;

///////////////////////////////////////////////////////////////////////////////
GenericGeometricalObject::GenericGeometricalObject()
{
	InitZero();
}
GenericGeometricalObject::GenericGeometricalObject(const GenericGeometricalObject *GGO)
{
	InitZero();
	
	GenericGeometricalObject::Copy(GGO);
}
GenericGeometricalObject::~GenericGeometricalObject()
{
	Clear();
}
int GenericGeometricalObject::Print(BuildWorldNI *buildWorld)
{
	int i;
	pnpPrintGroup0("BW      %s\n",name.c_str());

	if(buildWorld==NULL)
	{
		pnpPrintGroup0("BW        Dielectric Constant index:......... %d\n",Epsilon);
		pnpPrintGroup0("BW        Mobile Ions Diffusion index:....... [");
		for(i=0;i<NIonsTypes;i++)
		{
			pnpPrintGroup0("%d",IonsD[i]);
			if(i<NIonsTypes-1) pnpPrintGroup0(", ");
			else pnpPrintGroup0("]\n");
		}
	}
	else
	{
		pnpPrintGroup0("BW        Dielectric Constant:............... %.3f\n",buildWorld->DielConst[Epsilon-1]);
		pnpPrintGroup0("BW        Mobile Ions Diffusion Coefficient:. [");
		for(i=0;i<NIonsTypes;i++)
		{
			pnpPrintGroup0("%.3f",buildWorld->DiffusionConst[IonsD[i]-1]);
			if(i<NIonsTypes-1) pnpPrintGroup0(", ");
			else pnpPrintGroup0("] 1E-6 cm^2/s\n");
		}
	}
	pnpPrintGroup0("BW        Mobile Ions Concentration.......... [");
	for(i=0;i<NIonsTypes;i++)
	{
		pnpPrintGroup0("%.3f M",C[i]);
		if(i<NIonsTypes-1) pnpPrintGroup0(", ");
		else pnpPrintGroup0("]\n");
	}
	return EXIT_SUCCESS;
}
int GenericGeometricalObject::Copy(const GenericGeometricalObject *GGO)
{
	NIonsTypes=GGO->NIonsTypes;
	Epsilon=GGO->Epsilon;
	old_world=GGO->old_world;
	InternalUnit=GGO->InternalUnit;
	
	DeleteCArray(IonsD);
	DeleteCArray(C);
	IonsD=new int[NIonsTypes];
	C=new float[NIonsTypes];
	int i;
	for(i=0;i<NIonsTypes;i++)
	{
		IonsD[i]=GGO->IonsD[i];
		C[i]=GGO->C[i];
	}
	return EXIT_SUCCESS;
}
int GenericGeometricalObject::InitZero()
{
	PnpsObject::SetName("GenericGeometricalObject");
	NIonsTypes=0;
	Epsilon=0;
	IonsD=NULL;
	C=NULL;
	old_world=NULL;
	InternalUnit=false;
	ID_GridBegin=0;
	ID_GridEnd=0;
	return EXIT_SUCCESS;
}
int GenericGeometricalObject::SetGridIDs(int startIDfrom)
{
	ID_GridBegin=startIDfrom;
	ID_GridEnd=startIDfrom;
	return ID_GridEnd+1;
}
int GenericGeometricalObject::Clear()
{
	DeleteCArray(IonsD);
	DeleteCArray(C);
	return EXIT_SUCCESS;
}
int GenericGeometricalObject::SetBoundaryCondition(BuildWorldNI *Builder,ContWorld* world)
{
	return EXIT_SUCCESS;
}
int GenericGeometricalObject::setPar(int iEps, int iD0, int iD1,float C0,float C1)
{
	Clear();
	Epsilon=iEps+1;
	NIonsTypes=2;
	DeleteCArray(IonsD);
	IonsD=new int[NIonsTypes];
	IonsD[0]=iD0+1;
	IonsD[1]=iD1+1;
	
	DeleteCArray(C);
	C=new float[NIonsTypes];
	C[0]=C0;
	if(C1<0.0)
		C[1]=C0;
	else
		C[1]=C1;
	
	InternalUnit=false;
	old_world=NULL;
	return PNPS_EXIT_SUCCESS;
}
int GenericGeometricalObject::setParam(int iEps,int m_NIonsTypes, int *iD,float *m_C)
{
	Clear();
	int i;
	Epsilon=iEps+1;
	NIonsTypes=m_NIonsTypes;
	DeleteCArray(IonsD);
	IonsD=new int[NIonsTypes];
	for(i=0;i<NIonsTypes;i++)
		IonsD[i]=iD[i]+1;
	
	DeleteCArray(C);
	C=new float[NIonsTypes];
	for(i=0;i<NIonsTypes;i++)
		C[i]=m_C[i];
	
	InternalUnit=false;
	old_world=NULL;
	return PNPS_EXIT_SUCCESS;
}

int GenericGeometricalObject::ChangeUnits(ContWorld* world,bool ToInternal)
{
	if(world==NULL)
	{
		fprintf(stderr,"ERROR: ChangeUnits: world==NULL\n");
		return EXIT_FAILURE;
	}	
	if(old_world!=NULL&&old_world!=world&&ToInternal==false)
	{
		fprintf(stderr,"ERROR: ChangeUnits: Can't return to External Parameters using diffrent ContWorld\n");
		return EXIT_FAILURE;
	}
	//check if it is already right units
	if(InternalUnit==ToInternal&&old_world!=NULL)return EXIT_SUCCESS;
	
	//convert units
	int i;
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	
	if(ToInternal)
	{
		for(i=0;i<NIonsTypes;i++)
		{
			C[i]=C[i]*coef;
		}
	}
	else
	{
		for(i=0;i<NIonsTypes;i++)
		{
			C[i]=C[i]/coef;
		}
	}
	old_world=world;
	InternalUnit=ToInternal;
	return EXIT_SUCCESS;
}
int GenericGeometricalObject::BuildPreMaps(BuildWorldNI *Builder,ContWorld* world,int *Field,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth,float *Surf,int ParmMask)
{
	return EXIT_SUCCESS;
}
int GenericGeometricalObject::BuildPreMapsMemLim(PartialWorldForSAS *pw,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth)
{
	return EXIT_SUCCESS;
}
int GenericGeometricalObject::BuildLJRepultionMap(BuildWorldNI *Builder,ContWorld* world,int Ion,float *V,float LimitVlj)
{
	return EXIT_SUCCESS;
}
int GenericGeometricalObject::BuildDistMaps(BuildWorldNI *Builder,ContWorld* world,int *Field,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth,float *Surf)
{
	return EXIT_SUCCESS;
}
int GenericGeometricalObject::BuildDistMapsFromAtomsCenter(BuildWorldNI *Builder,ContWorld* world,float *Dist,float Rmax,float DoNotConsiderAtomsSmallerThen)
{
	return EXIT_SUCCESS;
}
int GenericGeometricalObject::BuildPreMapsCharges(BuildWorldNI *Builder,ContWorld* world,float *Field)
{
	return EXIT_SUCCESS;
}
int GenericGeometricalObject::RotateGGO(double *n, double cosa, double sina)
{
	pnpWarning("Cann't rotate this GenericGeometricalObject\n");
	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
GOAtoms::GOAtoms()
{
	InitZero();
}
GOAtoms::GOAtoms(const GOAtoms * goAtoms)
	:GenericGeometricalObject()
{
	InitZero();
	
	GOAtoms::Copy(goAtoms);
}
GOAtoms::~GOAtoms()
{
	Clear();
}
int GOAtoms::Copy(const GOAtoms *goAtoms)
{
	GenericGeometricalObject::Copy(goAtoms);
	ChargeDist=goAtoms->ChargeDist;
	ChargeDistN=goAtoms->ChargeDistN;
	
	int i;
	for(i=0;i<3;i++)
		Offset[i]=goAtoms->Offset[i];
	
	NAtoms=goAtoms->NAtoms;
	
	r[0].resize(NAtoms);
	r[1].resize(NAtoms);
	r[2].resize(NAtoms);
	R.resize(NAtoms);
	q.resize(NAtoms);
	
	for(i=0;i<NAtoms;i++)
	{
		r[0][i]=goAtoms->r[0][i];
		r[1][i]=goAtoms->r[1][i];
		r[2][i]=goAtoms->r[2][i];
		R[i]=goAtoms->R[i];
		q[i]=goAtoms->q[i];
	}
	MakePreRoll=goAtoms->MakePreRoll;
	return EXIT_SUCCESS;
}
int GOAtoms::InitZero()
{
	GenericGeometricalObject::InitZero();
	PnpsObject::SetName("AtomsParameters");
	ChargeDistStrs.push_back("ClosestNode");
	ChargeDistStrs.push_back("Linear1");
	ChargeDistStrs.push_back("Linear2");
	ChargeDistStrs.push_back("Linear3");
	
	ChargeDist=-1;
	ChargeDistN=-1;
	MakePreRoll=false;
	IAVMethod=0;

	r_x=&r[0];
	r_y=&r[1];
	r_z=&r[2];
	return EXIT_SUCCESS;
}
int GOAtoms::Clear()
{
	GenericGeometricalObject::Clear();
	return EXIT_SUCCESS;
}
int GOAtoms::Print(BuildWorldNI *buildWorld)
{
	GenericGeometricalObject::Print(buildWorld);
	pnpPrintGroup0("BW        Charge Distrebution Method:........ %s\n",ChargeDistStrs[ChargeDist].c_str());
	pnpPrintGroup0("BW        Atoms Number:...................... %d\n",NAtoms);
	return EXIT_SUCCESS;
}
int GOAtoms::DefineChargeDist(const char *str)
{
	if(strcmp("1-Node",str)==0)
		return ClosestNode;
	if(strcmp("ClosestNode",str)==0)
		return ClosestNode;
	if(strcmp("Linear1",str)==0)
	{
		ChargeDistN=1;
		return Linear;
	}
	if(strcmp("8-Nods",str)==0)
	{
		ChargeDistN=1;
		return Linear;
	}
	if(strcmp("Linear2",str)==0)
	{
		ChargeDistN=2;
		return Linear;
	}
	if(strcmp("Linear3",str)==0)
	{
		ChargeDistN=3;
		return Linear;
	}
	if(strcmp("Linear",str)==0)
	{
		return Linear;
	}
	if(strcmp("Cone",str)==0)
		return Cone;
	pnpPrint0("WARNING:was not able to find charge distrebution way will use Linear1\n");
	
	ChargeDistN=1;
	return Linear;
}
int GOAtoms::setPar(int iEps, int iD0, int iD1,float C0,float C1)
{
	Clear();
	Epsilon=iEps+1;
	NIonsTypes=2;
	DeleteCArray(IonsD);
	IonsD=new int[NIonsTypes];
	IonsD[0]=iD0+1;
	IonsD[1]=iD1+1;
	
	DeleteCArray(C);
	C=new float[NIonsTypes];
	C[0]=C0;
	if(C1<0.0)
		C[1]=C0;
	else
		C[1]=C1;
	
	InternalUnit=false;
	old_world=NULL;
	return PNPS_EXIT_SUCCESS;
}
int GOAtoms::SavePQR(const char *filename)
{
	return EXIT_SUCCESS;
}
int GOAtoms::LoadDelphiStyle(const char *PDB,const char *CRG,const char *SIZ,bool center)
{
	pnpError("NOT WORKING YET\n");
	int i,j;
	char Record[1024];
	char AtmName[6],ResName[5],OldResName[5];
	double tmpx,tmpy,tmpz,tmpq=0.0,tmpr=2.0;
	std::vector< std::string > AtmNameTMP;
	std::vector< double > ValTMP;
	bool Comments;
	
	
	
	//read Radius DB
	FILE* fsiz = fopen(SIZ,"rt");
	PNP_EXIT_FAIL_NULL1(fsiz,"Can't read %s\n",SIZ);
	
	std::vector< std::string > ResNameRadDB;
	std::vector< std::vector<std::string> > AtmNameRadDB;
	std::vector< std::vector<double> > RadDB;
	
	Comments=true;
	do
	{
		fgets(Record,1024,fsiz);
		if(Record[0]=='!')
			pnpPrint("Comments in %s: %s",SIZ,Record);
		else
			Comments=false;
	}
	while(Comments);
	
	pnpPrint("Table Type: %s",Record);
	if(strncmp(Record,"atom__res_radius_",17)!=0)
	{
		pnpError("unknown table format %d",strncmp(Record,"atom__res_radius_",17));
		return EXIT_FAILURE;
	}
	
	OldResName[0]='O';OldResName[1]='L';OldResName[2]='D';OldResName[3]='\0';
	while( fgets(Record,1024,fsiz) )
	{
		strncpy(AtmName,Record,6);AtmName[5]='\0';
		strncpy(ResName,Record+6,4);ResName[3]='\0';
		sscanf(Record+10,"%lf",&tmpr);
		//pnpPrint("rec|%.4s|%.6s|%f|\n",ResName,AtmName,tmpr);
		//pnpPrint("cur res %s %s %d\n",ResName,OldResName,strcmp(OldResName,ResName));
		if(strncmp(OldResName,ResName,4)!=0)
		{
			
			if(ValTMP.size()>0)
			{
				ResNameRadDB.push_back(OldResName);
				AtmNameRadDB.push_back(AtmNameTMP);
				RadDB.push_back(ValTMP);
				
				//pnpPrint("|%.4s|\n",OldResName);
				//for(j=0;j<ValTMP.size();j++)
				//	pnpPrint("\t|%.6s|%f|\n",AtmNameTMP[j].c_str(),ValTMP[j]);
				
				AtmNameTMP.resize(0);
				ValTMP.resize(0);
			}
			strncpy(OldResName,ResName,4);
		}
		
		AtmNameTMP.push_back(AtmName);
		ValTMP.push_back(tmpr);
	}
	//push back last entry
	ResNameRadDB.push_back(OldResName);
	AtmNameRadDB.push_back(AtmNameTMP);
	RadDB.push_back(ValTMP);
	AtmNameTMP.resize(0);
	ValTMP.resize(0);
	fclose(fsiz);
	/*for(i=0;i<ResNameRadDB.size();i++)
	{
		pnpPrint("|%.4s|\n",ResNameRadDB[i].c_str());
		for(j=0;j<AtmNameRadDB[i].size();j++)
			pnpPrint("\t|%.6s|%f|\n",AtmNameRadDB[i][j].c_str(),RadDB[i][j]);
	}*/
	
	//read charge DB
	FILE* fcrg = fopen(CRG,"rt");
	PNP_EXIT_FAIL_NULL1(fcrg,"Can't read %s\n",CRG);
	
	std::vector< std::string > ResNameCrgDB;
	std::vector< std::vector<std::string> > AtmNameCrgDB;
	std::vector< std::vector<double> > CrgDB;
	
	Comments=true;
	do
	{
		fgets(Record,1024,fsiz);
		if(Record[0]=='!')
			pnpPrint("Comments in %s: %s",SIZ,Record);
		else
			Comments=false;
	}
	while(Comments);
	
	pnpPrint("Table Type: %s",Record);
	if(strncmp(Record,"atom__resnumbc_charge_",17)!=0)
	{
		pnpError("unknown table format %d",strncmp(Record,"atom__resnumbc_charge_",22));
		return EXIT_FAILURE;
	}
	
	OldResName[0]='O';OldResName[1]='L';OldResName[2]='D';OldResName[3]='\0';
	while( fgets(Record,1024,fcrg) )
	{
		strncpy(AtmName,Record,6);AtmName[5]='\0';
		strncpy(ResName,Record+6,4);ResName[3]='\0';
		sscanf(Record+15,"%lf",&tmpq);
		//pnpPrint("rec|%.4s|%.6s|%f|\n",ResName,AtmName,tmpr);
		//pnpPrint("cur res %s %s %d\n",ResName,OldResName,strcmp(OldResName,ResName));
		if(strncmp(OldResName,ResName,4)!=0)
		{
			
			if(ValTMP.size()>0)
			{
				ResNameCrgDB.push_back(OldResName);
				AtmNameCrgDB.push_back(AtmNameTMP);
				CrgDB.push_back(ValTMP);
				
				//pnpPrint("|%.4s|\n",OldResName);
				//for(j=0;j<ValTMP.size();j++)
				//	pnpPrint("\t|%.6s|%f|\n",AtmNameTMP[j].c_str(),ValTMP[j]);
				
				AtmNameTMP.resize(0);
				ValTMP.resize(0);
			}
			strncpy(OldResName,ResName,4);
		}
		
		AtmNameTMP.push_back(AtmName);
		ValTMP.push_back(tmpq);
	}
	//push back last entry
	ResNameCrgDB.push_back(OldResName);
	AtmNameCrgDB.push_back(AtmNameTMP);
	CrgDB.push_back(ValTMP);
	AtmNameTMP.resize(0);
	ValTMP.resize(0);
	fclose(fcrg);
	/*for(i=0;i<ResNameCrgDB.size();i++)
	{
		pnpPrint("|%.4s|\n",ResNameCrgDB[i].c_str());
		for(j=0;j<AtmNameCrgDB[i].size();j++)
			pnpPrint("\t|%.6s|%f|\n",AtmNameCrgDB[i][j].c_str(),CrgDB[i][j]);
	}*/
	
	
	FILE* fpdb = fopen(PDB,"rt");

	PNP_EXIT_FAIL_NULL1(fpdb,"Can't read %s\n",PDB);

	
	while( fgets(Record,1024,fpdb) )
	{
		if(strncmp("ATOM",Record,4)==0||strncmp("HETATM",Record,6)==0)
		{
			strncpy(AtmName,Record+13,5);AtmName[4]='\0';
			strncpy(ResName,Record+17,3);ResName[3]='\0';
			//float tmpx,tmpy,tmpz,tmpq,tmpr;
			//sscanf(Record+27,"%f %f %f %f %f",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			sscanf(Record+27,"%lf %lf %lf",&tmpx,&tmpy,&tmpz);
			//pnpPrint0("\t%f %f %f %f %f\n",tmpx,tmpy,tmpz,tmpq,tmpr);
			r[0].push_back(float(tmpx));
			r[1].push_back(float(tmpy));
			r[2].push_back(float(tmpz));
			
			//find R in DB
			tmpr=0.0;
			int iRes=0;
			while(strncmp(ResNameRadDB[iRes].c_str(),ResName,3)!=0)
			{
				iRes++;
				if(ResNameRadDB.size()<=iRes)
				{
					iRes=-1;
					break;
				}
			}
			//if not in DB try to find with wild card
			if(iRes<0)
			{
				iRes=0;
				while(strncmp(ResNameRadDB[iRes].c_str(),"   ",3)!=0)
				{
					iRes++;
					if(ResNameRadDB.size()<=iRes)
					{
						iRes=-1;
						break;
					}
				}
			}
			int iAtm=0;
			if(iRes>=0)
			{
				while(strncmp(AtmNameRadDB[iRes][iAtm].c_str(),AtmName,4)!=0)
				{
					iAtm++;
					if(AtmNameRadDB[iRes].size()<=iAtm)
					{
						iAtm=-1;
						break;
					}
				}
				if(iAtm>=0)
					tmpr=RadDB[iRes][iAtm];
				else
				{
					
					//try to use wildcard
					for(i=3;i>0;i--)
					{
						iAtm=0;
						while(strncmp(AtmNameRadDB[iRes][iAtm].c_str(),AtmName,i)!=0)
						{
							iAtm++;
							if(AtmNameRadDB[iRes].size()<=iAtm)
							{
								iAtm=-1;
								break;
							}
						}
						if(iAtm>=0)
						{
							tmpr=RadDB[iRes][iAtm];
							break;
						}
					}
					if(iAtm>=0)
						pnpPrint("R:Was using a wildcard for atomtype %s in residue %s = %.4s %.5s %f\n",AtmName,ResName,ResNameRadDB[iRes].c_str(),AtmNameRadDB[iRes][iAtm].c_str(),tmpr);
					else
						pnpPrint("R:Can not find atomtype %s in residue %s in db, using value %f\n",AtmName,ResName,tmpr);
				}
			}
			else
			{
				pnpPrint("R:Can not find residue %s in db, using value %f\n",ResName,tmpr);
			}
			//find Q in DB
			tmpq=0.0;
			iRes=0;
			while(strncmp(ResNameCrgDB[iRes].c_str(),ResName,3)!=0)
			{
				iRes++;
				if(ResNameCrgDB.size()<=iRes)
				{
					iRes=-1;
					break;
				}
			}
			//if not in DB try to find with wild card
			/*if(iRes<0)
			{
				iRes=0;
				while(strncmp(ResNameCrgDB[iRes].c_str(),"   ",3)!=0)
				{
					iRes++;
					if(ResNameCrgDB.size()<=iRes)
					{
						iRes=-1;
						break;
					}
				}
			}*/
			iAtm=0;
			if(iRes>=0)
			{
				while(strncmp(AtmNameCrgDB[iRes][iAtm].c_str(),AtmName,4)!=0)
				{
					iAtm++;
					if(AtmNameCrgDB[iRes].size()<=iAtm)
					{
						iAtm=-1;
						break;
					}
				}
				if(iAtm>=0)
					tmpq=CrgDB[iRes][iAtm];
				else
				{
					
					//try to use wildcard
					/*for(i=3;i>0;i--)
					{
						iAtm=0;
						while(strncmp(AtmNameCrgDB[iRes][iAtm].c_str(),AtmName,i)!=0)
						{
							iAtm++;
							if(AtmNameCrgDB[iRes].size()<=iAtm)
							{
								iAtm=-1;
								break;
							}
						}
						if(iAtm>=0)
						{
							tmpq=CrgDB[iRes][iAtm];
							break;
						}
					}
					if(iAtm>=0)
						pnpPrint("Was using a wildcard for atomtype %s in residue %s = %.4s %.5s %f\n",AtmName,ResName,ResNameCrgDB[iRes].c_str(),AtmNameCrgDB[iRes][iAtm].c_str(),tmpq);
					else*/
						pnpPrint("Q:Can not find atomtype %s in residue %s in db, using value %f\n",AtmName,ResName,tmpq);
				}
			}
			else
			{
				pnpPrint("Q:Can not find residue %s in db, using value %f\n",ResName,tmpq);
			}
			q.push_back(float(tmpq));
			R.push_back(float(tmpr));
			pnpPrint("PRQ: %.27s %f %f\n",Record,tmpr,tmpq);
		}
		
	}
	fclose(fpdb);
	pnpPrint0("%d atoms was read from %s\n",r[0].size(),PDB);
	return EXIT_SUCCESS;
}
int GOAtoms::LoadPQR(const char *filename)
{
	FILE* fpqr = fopen(filename,"rt");

	PNP_EXIT_FAIL_NULL1(fpqr,"Can't read %s\n",filename);

	char Record[1024];
	while( fgets(Record,1024,fpqr) )
	{
		if(strncmp("ATOM",Record,4)==0||strncmp("HETATM",Record,6)==0)
		{
			//float tmpx,tmpy,tmpz,tmpq,tmpr;
			//sscanf(Record+27,"%f %f %f %f %f",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			double tmpx,tmpy,tmpz,tmpq,tmpr;
			sscanf(Record+27,"%lf %lf %lf %lf %lf",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			//pnpPrint0("\t%f %f %f %f %f\n",tmpx,tmpy,tmpz,tmpq,tmpr);
			r[0].push_back(float(tmpx));
			r[1].push_back(float(tmpy));
			r[2].push_back(float(tmpz));
			q.push_back(float(tmpq));
			R.push_back(float(tmpr));
		}
		
	}
	fclose(fpqr);
	pnpPrint0("%d atoms was read from %s\n",r[0].size(),filename);
	NAtoms=R.size();
	return EXIT_SUCCESS;
}
int GOAtoms::LoadAtomicParam(const char *filename)
{
	FILE* fpqr = fopen(filename,"rt");

	PNP_EXIT_FAIL_NULL1(fpqr,"Can't read %s\n",filename);
	IER.resize(NIonsTypes);
	SRA.resize(NIonsTypes);
	SRN.resize(NIonsTypes);
	char Record[1024];
	while( fgets(Record,1024,fpqr) )
	{
		if(!(Record[0]=='@' || Record[0]=='#'))
		{
			//float tmpx,tmpy,tmpz,tmpq,tmpr;
			//sscanf(Record+27,"%f %f %f %f %f",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			int tmpResNum, tmpAtmNum;
			char tmpResName[32],tmpAtmName[32];
			double tmpx,tmpy,tmpz,tmpq,tmpr;
			//ResNum  ResName    AtmNum  AtmName               x               y               z               Q           Rdiel 
			//Rier1           Rier2           SR-A1           SR-N1           SR-A2           SR-N2
			double tmpRier1,tmpRier2;
			double tmpSR_A1, tmpSR_N1, tmpSR_A2, tmpSR_N2;
			sscanf(Record,"%d %s %d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&tmpResNum, tmpResName, &tmpAtmNum, tmpAtmName, &tmpx,&tmpy,&tmpz,&tmpq,&tmpr,&tmpRier1,&tmpRier2,&tmpSR_A1, &tmpSR_N1, &tmpSR_A2, &tmpSR_N2);
			//pnpPrint0("\t%f %f %f %f %f\n",tmpx,tmpy,tmpz,tmpq,tmpr);
			r[0].push_back(float(tmpx));
			r[1].push_back(float(tmpy));
			r[2].push_back(float(tmpz));
			q.push_back(float(tmpq));
			R.push_back(float(tmpr));
			IER[0].push_back(float(tmpRier1));
			IER[1].push_back(float(tmpRier2));
			SRA[0].push_back(float(tmpSR_A1));
			SRN[0].push_back(float(tmpSR_N1));
			SRA[1].push_back(float(tmpSR_A2));
			SRN[1].push_back(float(tmpSR_N2));
		}
		
	}
	fclose(fpqr);
	pnpPrint0("%d atoms was read from %s\n",r[0].size(),filename);
	NAtoms=R.size();
	IAVMethod=3;
	return EXIT_SUCCESS;
}
int GOAtoms::SetCenterToOrigin()
{
	float xmin=r[0][0];
	float xmax=r[0][0];
	float ymin=r[1][0];
	float ymax=r[1][0];
	float zmin=r[2][0];
	float zmax=r[2][0];
	
	int i;
	for(i=1;i<R.size();i++)
	{
		if(xmin>r[0][i])xmin=r[0][i];
		if(xmax<r[0][i])xmax=r[0][i];
		if(ymin>r[1][i])ymin=r[1][i];
		if(ymax<r[1][i])ymax=r[1][i];
		if(zmin>r[2][i])zmin=r[2][i];
		if(zmax<r[2][i])zmax=r[2][i];
	}
	pnpPrint("MinMax Coordinated [%.3f,%.3f,%.3f]-[%.3f,%.3f,%.3f] to 0\n",xmin,ymin,zmin,xmax,ymax,zmax);
	float xc=0.5*(xmax+xmin);
	float yc=0.5*(ymax+ymin);
	float zc=0.5*(zmax+zmin);
	pnpPrint("Moving center[%.3f,%.3f,%.3f] to 0\n",xc,yc,zc);
	for(i=0;i<R.size();i++)
	{
		r[0][i]-=xc;
		r[1][i]-=yc;
		r[2][i]-=zc;
	}
	return 1;
}
int GOAtoms::SetCenterToOriginWithRad()
{
	float xmin=r[0][0]-R[0];
	float xmax=r[0][0]+R[0];
	float ymin=r[1][0]-R[0];
	float ymax=r[1][0]+R[0];
	float zmin=r[2][0]-R[0];
	float zmax=r[2][0]+R[0];
	
	int i;
	for(i=1;i<R.size();i++)
	{
		if(xmin>r[0][i]-R[i])xmin=r[0][i]-R[i];
		if(xmax<r[0][i]+R[i])xmax=r[0][i]+R[i];
		if(ymin>r[1][i]-R[i])ymin=r[1][i]-R[i];
		if(ymax<r[1][i]+R[i])ymax=r[1][i]+R[i];
		if(zmin>r[2][i]-R[i])zmin=r[2][i]-R[i];
		if(zmax<r[2][i]+R[i])zmax=r[2][i]+R[i];
	}
	pnpPrint("MinMax Coordinated [%.3f,%.3f,%.3f]-[%.3f,%.3f,%.3f]\n",xmin,ymin,zmin,xmax,ymax,zmax);
	pnpPrint("Range [%.3f,%.3f,%.3f]\n",xmax-xmin,ymax-ymin,zmax-zmin);
	float xc=0.5*(xmax+xmin);
	float yc=0.5*(ymax+ymin);
	float zc=0.5*(zmax+zmin);
	pnpPrint("Moving center[%.3f,%.3f,%.3f] to 0\n",xc,yc,zc);
	for(i=0;i<R.size();i++)
	{
		r[0][i]-=xc;
		r[1][i]-=yc;
		r[2][i]-=zc;
	}
	return 1;
}
int GOAtoms::LoadPRQ(const char *filename)
{
	FILE* fpqr = fopen(filename,"rt");

	PNP_EXIT_FAIL_NULL1(fpqr,"Can't read %s\n",filename);

	char Record[1024];
	while( fgets(Record,1024,fpqr) )
	{
		if(strncmp("ATOM",Record,4)==0||strncmp("HETATM",Record,6)==0)
		{
			//float tmpx,tmpy,tmpz,tmpq,tmpr;
			//sscanf(Record+27,"%f %f %f %f %f",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			double tmpx,tmpy,tmpz,tmpq,tmpr;
			sscanf(Record+27,"%lf %lf %lf %lf %lf",&tmpx,&tmpy,&tmpz,&tmpr,&tmpq);
			//pnpPrint0("\t%f %f %f %f %f\n",tmpx,tmpy,tmpz,tmpq,tmpr);
			r[0].push_back(float(tmpx));
			r[1].push_back(float(tmpy));
			r[2].push_back(float(tmpz));
			q.push_back(float(tmpq));
			R.push_back(float(tmpr));
		}
		
	}
	fclose(fpqr);
	pnpPrint0("%d atoms was read from %s\n",r[0].size(),filename);
	return EXIT_SUCCESS;
}
int GOAtoms::LoadPQRonlyQR(const char *filename)
{
	FILE* fpqr = fopen(filename,"rt");

	PNP_EXIT_FAIL_NULL1(fpqr,"Can't read %s\n",filename);

	char Record[1024];
	int i=0;
	while( fgets(Record,1024,fpqr) )
	{
		if(strncmp("ATOM",Record,4)==0||strncmp("HETATM",Record,6)==0)
		{
			//float tmpx,tmpy,tmpz,tmpq,tmpr;
			//sscanf(Record+27,"%f %f %f %f %f",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			double tmpx,tmpy,tmpz,tmpq,tmpr;
			sscanf(Record+27,"%lf %lf %lf %lf %lf",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			q[i]=float(tmpq);
			R[i]=float(tmpr);
			i++;
			if(i>NAtoms)
			{
				pnpError("number atoms in %s more then already initiated will ignore the rest\n",filename);
				break;
			}
		}
		
	}
	if(i!=NAtoms)
	{
		pnpError("number atoms in %s (%d) is not coinside with already initiated(%d)\n",filename,i,NAtoms);
	}
	fclose(fpqr);
	pnpPrint0("qr for %d atoms was read from %s\n",r[0].size(),filename);
	return EXIT_SUCCESS;
}
int GOAtoms::LoadPRE(const char *filename)
{
	FILE* fpqr = fopen(filename,"rt");

	PNP_EXIT_FAIL_NULL1(fpqr,"Can't read %s\n",filename);

	char Record[1024];
	while( fgets(Record,1024,fpqr) )
	{
		if(strncmp("ATOM",Record,4)==0||strncmp("HETATM",Record,6)==0)
		{
			//float tmpx,tmpy,tmpz,tmpq,tmpr;
			//sscanf(Record+27,"%f %f %f %f %f",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			double tmpx,tmpy,tmpz,tmpq,tmpr;
			sscanf(Record+27,"%lf %lf %lf %lf %lf",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			//pnpPrint0("\t%f %f %f %f %f\n",tmpx,tmpy,tmpz,tmpq,tmpr);
			//r[0].push_back(float(tmpx));
			//r[1].push_back(float(tmpy));
			//r[2].push_back(float(tmpz));
			HalfSigma.push_back(float(tmpq));
			FourEpsilon.push_back(float(tmpr));
		}
		
	}
	fclose(fpqr);
	pnpPrint0("%d atoms was read from %s\n",r[0].size(),filename);
	return EXIT_SUCCESS;
}
int GOAtoms::LoadPAB(const char *filename)
{
	FILE* fpqr = fopen(filename,"rt");

	PNP_EXIT_FAIL_NULL1(fpqr,"Can't read %s\n",filename);

	char Record[1024];
	LJA.resize(NIonsTypes);
	LJB.resize(NIonsTypes);
	while( fgets(Record,1024,fpqr) )
	{
		if(strncmp("ATOM",Record,4)==0||strncmp("HETATM",Record,6)==0)
		{
			//float tmpx,tmpy,tmpz,tmpq,tmpr;
			//sscanf(Record+27,"%f %f %f %f %f",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			double tmpA1,tmpB1,tmpA2,tmpB2;
			sscanf(Record+27,"%le %le %le %le",&tmpA1,&tmpB1,&tmpA2,&tmpB2);
			LJA[0].push_back(float(tmpA1));
			LJB[0].push_back(float(tmpB1));
			LJA[1].push_back(float(tmpA2));
			LJB[1].push_back(float(tmpB2));
		}
		
	}
	fclose(fpqr);
	pnpPrint0("LoadPAB: %d atoms was read from %s\n",LJA[0].size(),filename);
	IAVMethod=2;
	return EXIT_SUCCESS;
}
int GOAtoms::LoadPAN(const char *filename)
{
	FILE* fpqr = fopen(filename,"rt");

	PNP_EXIT_FAIL_NULL1(fpqr,"Can't read %s\n",filename);

	char Record[1024];
	SRA.resize(NIonsTypes);
	SRN.resize(NIonsTypes);
	while( fgets(Record,1024,fpqr) )
	{
		if(strncmp("ATOM",Record,4)==0||strncmp("HETATM",Record,6)==0)
		{
			//float tmpx,tmpy,tmpz,tmpq,tmpr;
			//sscanf(Record+27,"%f %f %f %f %f",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			double tmpA1,tmpB1,tmpA2,tmpB2;
			sscanf(Record+27,"%le %le %le %le",&tmpA1,&tmpB1,&tmpA2,&tmpB2);
			SRA[0].push_back(float(tmpA1));
			SRN[0].push_back(float(tmpB1));
			SRA[1].push_back(float(tmpA2));
			SRN[1].push_back(float(tmpB2));
		}
		
	}
	fclose(fpqr);
	pnpPrint0("LoadPAN: %d atoms was read from %s\n",SRA[0].size(),filename);
	IAVMethod=3;
	return EXIT_SUCCESS;
}
int GOAtoms::LoadIER(const char *filename)
{
	IER.resize(2);
	
	FILE* fpqr = fopen(filename,"rt");

	PNP_EXIT_FAIL_NULL1(fpqr,"Can't read %s\n",filename);

	char Record[1024];
	while( fgets(Record,1024,fpqr) )
	{
		if(strncmp("ATOM",Record,4)==0||strncmp("HETATM",Record,6)==0)
		{
			//float tmpx,tmpy,tmpz,tmpq,tmpr;
			//sscanf(Record+27,"%f %f %f %f %f",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
			double tmp1,tmp2;
			sscanf(Record+27,"%lf %lf",&tmp1,&tmp2);
			//pnpPrint0("\t%f %f %f %f %f\n",tmpx,tmpy,tmpz,tmpq,tmpr);
			//r[0].push_back(float(tmpx));
			//r[1].push_back(float(tmpy));
			//r[2].push_back(float(tmpz));
			IER[0].push_back(float(tmp1));
			IER[1].push_back(float(tmp2));
		}
		
	}
	fclose(fpqr);
	pnpPrint0("LoadIER: %d atoms was read from %s\n",r[0].size(),filename);
	IAVMethod=1;
	return EXIT_SUCCESS;
}
#ifdef HARLEM_MOD
int GOAtoms::LoadHaMolSet(HaMolSet *molset)
{
	HaAtom* aptr;
	AtomIteratorMolSet aitr(molset);

	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if(molset->p_save_opt_default->save_selected && !aptr->Selected())
			continue;
				
		r[0].push_back(float(aptr->GetX()));
		r[1].push_back(float(aptr->GetY()));
		r[2].push_back(float(aptr->GetZ()));
		q.push_back(float(aptr->GetCharge()));
		R.push_back(float(aptr->radius));
	}
	NAtoms=R.size();
	return EXIT_SUCCESS;
}
#endif
int GOAtoms::ChangeUnits(ContWorld* world,bool ToInternal)
{
	if(world==NULL)
	{
		fprintf(stderr,"ERROR: ChangeUnits: world==NULL\n");
		return EXIT_FAILURE;
	}	
	if(old_world!=NULL&&old_world!=world&&ToInternal==false)
	{
		fprintf(stderr,"ERROR: ChangeUnits: Can't return to External Parameters using diffrent ContWorld\n");
		return EXIT_FAILURE;
	}
	//check if it is already right units(is made in pre)
	if(InternalUnit==ToInternal&&old_world!=NULL)return EXIT_SUCCESS;
	
	//Do GGO stuff
	GenericGeometricalObject::ChangeUnits(world,ToInternal);
	
	//convert units
	int i,j;
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	
	double GS2=world->GridScale*world->GridScale;
	double GS6=GS2*GS2*GS2;
	double GS12=GS6*GS6;
	
	int irc[3];
	float frc[3];
	for(j=0;j<3;j++)
	{
		irc[j]=world->GridSizeGlobal[j]/2;
		frc[j]=irc[j];
	}
	if(ToInternal)
	{
		//printf("ChangeUnits: %f %f %f %f\n",r[0][0],r[1][0],r[2][0],R[0]);
		for(j=0;j<3;j++)
			for(i=0;i<NAtoms;i++)
				r[j][i]=r[j][i]*world->GridScale+frc[j];
		
		for(i=0;i<NAtoms;i++)
		{
			R[i]*=world->GridScale;
			q[i]*=fpoh;
		}
		for(i=0;i<HalfSigma.size();i++)
		{
			HalfSigma[i]*=world->GridScale;
		}
		for(j=0;j<LJA.size();j++)
			for(i=0;i<LJA[j].size();i++)
		{
			LJA[j][i]=LJA[j][i]*GS12;
		}
		for(j=0;j<LJB.size();j++)
			for(i=0;i<LJB[j].size();i++)
		{
			LJB[j][i]=LJB[j][i]*GS6;
		}
		for(i=0;i<IER.size();i++)
			for(j=0;j<IER[i].size();j++)
				IER[i][j]*=world->GridScale;
		for(i=0;i<SRA.size();i++)
			for(j=0;j<SRA[i].size();j++)
				SRA[i][j]*=world->GridScale;
		//printf("ChangeUnits: %f %f %f %f\n",r[0][0],r[1][0],r[2][0],R[0]);
	}
	else
	{
		for(j=0;j<3;j++)
			for(i=0;i<NAtoms;i++)
				r[j][i]=(r[j][i]-frc[j])/world->GridScale;
		
		for(i=0;i<NAtoms;i++)
		{
			R[i]/=world->GridScale;
			q[i]/=fpoh;
		}
		for(i=0;i<HalfSigma.size();i++)
		{
			HalfSigma[i]/=world->GridScale;
		}
		for(j=0;j<LJA.size();j++)
			for(i=0;i<LJA[j].size();i++)
		{
			LJA[j][i]=LJA[j][i]/GS12;
		}
		for(j=0;j<LJB.size();j++)
			for(i=0;i<LJB[j].size();i++)
		{
			LJB[j][i]=LJB[j][i]/GS6;
		}
		for(i=0;i<IER.size();i++)
			for(j=0;j<IER[i].size();j++)
				IER[i][j]/=world->GridScale;
		for(i=0;i<SRA.size();i++)
			for(j=0;j<SRA[i].size();j++)
				SRA[i][j]/=world->GridScale;
	}
	return EXIT_SUCCESS;
}
int GOAtoms::BuildDistMapsFromAtomsCenter(BuildWorldNI *Builder,ContWorld* world,float *Dist,float Rmax,float DoNotConsiderAtomsSmallerThen)
{
	DbgPrint2("GOAtoms::BuildDistMapsFromAtomsCenter\n");
	int i,j,k,gridpoint,rint[3],gridPoint2;
	int GridSizeXY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
	float RSQ,RsmSQ,RtmpSQ,Ri,Rsm,Rtmp;
	float RleftSQ,RrightSQ;
	float rf[3];
	int iR,iRsm;
	int ix,iy,iz;
	float ftmp;
	float *vtmp;
	
	int *GridSize=world->GridSizeGlobal;
	int GSX=GridSize[0];
	int GSXY=GridSize[0]*GridSize[1];
	int GSXYZ=GridSize[0]*GridSize[1]*GridSize[2];
	
	int start[3];
	int end[3];
	
	for(i=0;i<NAtoms;i++)
	{
		if(R[i]<DoNotConsiderAtomsSmallerThen)
			continue;
		Ri=Rmax;
		iR=(int)(Ri+0.5);
		RSQ=Ri*Ri;
		for(k=0;k<3;k++){
			rf[k]=r[k][i];
			rint[k]=(int)(rf[k]+0.5);
			start[k]=rint[k]-iR;
			end[k]=rint[k]+iR;
			if(start[k]<0)start[k]=0;
			if(end[k]>GridSize[k]-1)start[k]=GridSize[k]-1;
		}
			
		for(ix=start[0];ix<=end[0];ix++)
			for(iy=start[1];iy<=end[1];iy++)
				for(iz=start[2];iz<=end[2];iz++)
		{
			gridpoint=ix+iy*GSX+iz*GSXY;
			RtmpSQ=(rf[0]-ix)*(rf[0]-ix)+(rf[1]-iy)*(rf[1]-iy)+(rf[2]-iz)*(rf[2]-iz);
			RtmpSQ=sqrt(RtmpSQ);
			if(RtmpSQ<Rmax)
				if(RtmpSQ<Dist[gridpoint])//RtmpSQ<=RSQ
			{
				Dist[gridpoint]=RtmpSQ;
			}
		}
	}
	
	return EXIT_SUCCESS;
}
int GOAtoms::BuildPreMaps(BuildWorldNI *Builder,ContWorld* world,int *Field,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth,float *Surf,int ParmMask)
{
	DbgPrint2("GOAtoms::BuildPreMaps(iValue=%d iBulkValue=%d)\n"
			,iValue,iBulkValue);
// 	DbgPrint2("\t Rion=%f[grids] Rsmoth=%f[grids] \nDispl=[%f,%f,%f][grids,grids,grids]\n", Rion, Rsmoth, Displ[0], Displ[1], Displ[2]);
// 	DbgPrint2("\t r[0]=[%f,%f,%f] R=%f q=%f Natom=%d\n"
// 			,r[0][0],r[1][0],r[2][0],R[0],q[0],NAtoms);
// 	DbgPrint2("\t r[0]=[%f,%f,%f] R=%f q=%f Natom=%d\n"
// 			,r[0][NAtoms-1],r[1][NAtoms-1],r[2][NAtoms-1],R[NAtoms-1],q[NAtoms-1],NAtoms);
	int i,j,k,gridpoint,rint[3],gridPoint2;
	int GridSizeXY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
	float RSQ,RsmSQ,RtmpSQ,Ri,Rsm,Rtmp;
	float RleftSQ,RrightSQ;
	float rf[3];
	int iR,iRsm;
	int ix,iy,iz;
	float ftmp;
	float *vtmp;
	
	int *GridSize=world->GridSizeGlobal;
	int GSX=GridSize[0];
	int GSXY=GridSize[0]*GridSize[1];
	int GSXYZ=GridSize[0]*GridSize[1]*GridSize[2];
	
// 	float Rmax=0.0;
// 	for(i=0;i<NAtoms;i++)
// 	{
// 		if(Rmax<R[i])Rmax=R[i];
// 	}
// 	Rmax+=Rion;
// 	float RmaxWpr=Rmax+Rsmoth;
// 	float iRmaxWpr=int(float(RmaxWpr+1.5));
// 	
	int start[3];
	int end[3];
	
	bool bIER=false;
	int Ion=0;
	if(ParmMask>=1&&IER.size()>0)
		if(IER[0].size()>0)
	{
		bIER=true;
		Ion=ParmMask-1;
	}
	DbgPrint0("bIER=%d",bIER);
	if(Rsmoth==0.0f)
	{
		DbgPrint0("Rsmoth==0.0f");
		for(i=0;i<NAtoms;i++)
		{
			if(bIER)
				Ri=IER[Ion][i];
			else
				Ri=R[i]+Rion;
			
			iR=(int)(Ri+0.5);
			RSQ=Ri*Ri;
			for(k=0;k<3;k++){
				rf[k]=r[k][i]+Displ[k];
				rint[k]=(int)(rf[k]+0.5);
				start[k]=rint[k]-iR;
				end[k]=rint[k]+iR;
				if(start[k]<0)start[k]=0;
				if(end[k]>GridSize[k]-1)end[k]=GridSize[k]-1;
			}
			
			for(ix=start[0];ix<=end[0];ix++)
				for(iy=start[1];iy<=end[1];iy++)
					for(iz=start[2];iz<=end[2];iz++)
			{
				gridpoint=ix+iy*GSX+iz*GSXY;
				RtmpSQ=(rf[0]-ix)*(rf[0]-ix)+(rf[1]-iy)*(rf[1]-iy)+(rf[2]-iz)*(rf[2]-iz);
				vtmp=Surf+4*gridpoint;
				if(RtmpSQ<=RSQ)//RtmpSQ<=RSQ
				{
					Field[gridpoint]=iValue;
					vtmp[0]=-BIGDISTANSE;
				}
			}
		}
	}
	else
	{
		for(i=0;i<NAtoms;i++)
		{
			
			if(bIER)
				Ri=IER[Ion][i];
			else
				Ri=R[i]+Rion;
			if(Ri==0.0){continue;}
			
			
			Rsm=Ri+Rsmoth;
			iR=(int)(Ri+0.5);
			iRsm=(int)(Rsm+1.0);
			RleftSQ=Rsm-0.707106781f;
			RrightSQ=Rsm+0.707106781f;
			RleftSQ*=RleftSQ;
			RrightSQ*=RrightSQ;
			RSQ=Ri*Ri;
			RsmSQ=Rsm*Rsm;
			
			for(k=0;k<3;k++)
			{
				rf[k]=r[k][i]+Displ[k];
				rint[k]=(int)(rf[k]+0.5);
				start[k]=rint[k]-iRsm;
				end[k]=rint[k]+iRsm;
				if(start[k]<0)start[k]=0;
				if(end[k]>GridSize[k]-1)end[k]=GridSize[k]-1;
			}
			for(ix=start[0];ix<=end[0];ix++)
				for(iy=start[1];iy<=end[1];iy++)
					for(iz=start[2];iz<=end[2];iz++)
			{
				gridpoint=ix+iy*GSX+iz*GSXY;
				RtmpSQ=(rf[0]-ix)*(rf[0]-ix)+(rf[1]-iy)*(rf[1]-iy)+(rf[2]-iz)*(rf[2]-iz);
				vtmp=Surf+4*gridpoint;
				if(RtmpSQ<=RSQ)//RtmpSQ<=RSQ
				{
					Field[gridpoint]=iValue;
					vtmp[0]=-BIGDISTANSE;
					vtmp[3]=0.0;
				}
				else//RSQ<RtmpSQ
				{
					Rtmp=sqrt(RtmpSQ);
					if(RtmpSQ<=RsmSQ)//RSQ<RtmpSQ<RleftSQ
					{
						if(Field[gridpoint]==iBulkValue)
						{
							vtmp[0]=-BIGDISTANSE;
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp-Ri;
						}
						else if(Field[gridpoint]<0&&Rtmp-Ri<vtmp[3])
						{
							vtmp[0]=-BIGDISTANSE;
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp-Ri;
						}
					}
					if(RleftSQ<=RtmpSQ&&RtmpSQ<=RrightSQ)//RleftSQ<=RtmpSQ<=RrightSQ
					{
						if(vtmp[0]>-100.0f)//e.i. is intersection
						{
							if(Rtmp-Ri<vtmp[3])
							{
								ftmp=Rsm/Rtmp;
								vtmp[0]=rf[0]+ftmp*(ix-rf[0]);
								vtmp[1]=rf[1]+ftmp*(iy-rf[1]);
								vtmp[2]=rf[2]+ftmp*(iz-rf[2]);
								vtmp[3]=Rtmp-Ri;
							}
						}
						else if(Field[gridpoint]==iBulkValue)
						{
							ftmp=Rsm/Rtmp;
							vtmp[0]=rf[0]+ftmp*(ix-rf[0]);
							vtmp[1]=rf[1]+ftmp*(iy-rf[1]);
							vtmp[2]=rf[2]+ftmp*(iz-rf[2]);
							vtmp[3]=Rtmp-Ri;
						}
					}
				}
			}
		}
		if(MakePreRoll)
		{
			
			DbgPrint0("MakePreRoll\n");
			int lim=(int)(Rsmoth+0.5);
			int itmp=GSXYZ*4;
			int iSolAcFlag=1000000000;
			RsmSQ=Rsmoth*Rsmoth;
			for(i=0;i<itmp;i=i+4)
			{
				if(Surf[i]>-100.0)
				{
					vtmp=Surf+i;
					for(k=0;k<3;k++)
					{
						rint[k]=(int)(vtmp[k]+0.5);
						start[k]=rint[k]-lim;
						end[k]=rint[k]+lim;
						if(start[k]<0)start[k]=0;
						if(end[k]>GridSize[k]-1)end[k]=GridSize[k]-1;
					}
					for(ix=start[0];ix<=end[0];ix++)
						for(iy=start[1];iy<=end[1];iy++)
							for(iz=start[2];iz<=end[2];iz++)
					{
						gridpoint=ix+iy*GSX+iz*GSXY;
						RtmpSQ=(vtmp[0]-ix)*(vtmp[0]-ix)+(vtmp[1]-iy)*(vtmp[1]-iy)+(vtmp[2]-iz)*(vtmp[2]-iz);
						if(RtmpSQ<RsmSQ&&Field[gridpoint]<0)
						{
							Field[gridpoint]=iSolAcFlag;
						}
					}
				}
			}
			for(i=0;i<GSXYZ;i++)
			{
				if(Field[i]<0)
				{
					Field[i]=-Field[i];
				}
				else if(Field[i]==iSolAcFlag)
				{
					Field[i]=-iValue;
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
int GOAtoms::BuildPreMapsMemLim(PartialWorldForSAS *pw,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth)
{
	DbgPrint2("GOAtoms::BuildPreMaps(iValue=%d iBulkValue=%d)\n"
			,iValue,iBulkValue);
// 	DbgPrint2("\t Rion=%f[grids] Rsmoth=%f[grids] \nDispl=[%f,%f,%f][grids,grids,grids]\n", Rion, Rsmoth, Displ[0], Displ[1], Displ[2]);
// 	DbgPrint2("\t r[0]=[%f,%f,%f] R=%f q=%f Natom=%d\n"
// 			,r[0][0],r[1][0],r[2][0],R[0],q[0],NAtoms);
// 	DbgPrint2("\t r[0]=[%f,%f,%f] R=%f q=%f Natom=%d\n"
// 			,r[0][NAtoms-1],r[1][NAtoms-1],r[2][NAtoms-1],R[NAtoms-1],q[NAtoms-1],NAtoms);
	int i,j,k,gridpoint,rint[3],gridPoint2;
	//int GridSizeXY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
	float RSQ,RsmSQ,RtmpSQ,Ri,Rsm,Rtmp;
	float RleftSQ,RrightSQ;
	float rf[3];
	int iR,iRsm;
	int ix,iy,iz;
	float ftmp;
	float *vtmp;
	
	
	float *Surf=pw->Surf;
	
	int GSX=pw->GS_X;
	int GSY=pw->GS_Y;
	int GSZ=pw->locGS_Z;
	int GSXY=pw->GS_X*pw->GS_Y;
	int GSXYZ=pw->GS_X*pw->GS_Y*pw->locGS_Z;
	
	int *Field=pw->Field+GSXY*pw->locR0_Z;
	int GridSize[3]={GSX,GSY,GSZ};
	
// 	float Rmax=0.0;
// 	for(i=0;i<NAtoms;i++)
// 	{
// 		if(Rmax<R[i])Rmax=R[i];
// 	}
// 	Rmax+=Rion;
// 	float RmaxWpr=Rmax+Rsmoth;
// 	float iRmaxWpr=int(float(RmaxWpr+1.5));
	// 	
	int start[3];
	int end[3];
	
	bool bIER=false;
	int Ion=0;
// 	if(ParmMask>=1&&IER.size()>0)
// 		if(IER[0].size()>0)
// 	{
// 		bIER=true;
// 		Ion=ParmMask-1;
// 	}
	DbgPrint0("bIER=%d",bIER);
	if(Rsmoth==0.0f)
	{
		DbgPrint0("Rsmoth==0.0f");
		for(i=0;i<NAtoms;i++)
		{
			if(bIER)
				Ri=IER[Ion][i];
			else
				Ri=R[i]+Rion;
			
			iR=(int)(Ri+0.5);
			RSQ=Ri*Ri;
			for(k=0;k<3;k++){
				rf[k]=r[k][i]+Displ[k];
				rint[k]=(int)(rf[k]+0.5);
				start[k]=rint[k]-iR;
				end[k]=rint[k]+iR;
				if(start[k]<0)start[k]=0;
				if(end[k]>GridSize[k]-1)end[k]=GridSize[k]-1;
			}
			
			for(ix=start[0];ix<=end[0];ix++)
				for(iy=start[1];iy<=end[1];iy++)
					for(iz=start[2];iz<=end[2];iz++)
			{
				gridpoint=ix+iy*GSX+iz*GSXY;
				RtmpSQ=(rf[0]-ix)*(rf[0]-ix)+(rf[1]-iy)*(rf[1]-iy)+(rf[2]-iz)*(rf[2]-iz);
				vtmp=Surf+4*gridpoint;
				if(RtmpSQ<=RSQ)//RtmpSQ<=RSQ
				{
					Field[gridpoint]=iValue;
					vtmp[0]=-BIGDISTANSE;
				}
			}
		}
	}
	else
	{
		for(i=0;i<NAtoms;i++)
		{
			
			if(bIER)
				Ri=IER[Ion][i];
			else
				Ri=R[i]+Rion;
			if(Ri==0.0){continue;}
			
			
			Rsm=Ri+Rsmoth;
			iR=(int)(Ri+0.5);
			iRsm=(int)(Rsm+1.0);
			RleftSQ=Rsm-0.707106781f;
			RrightSQ=Rsm+0.707106781f;
			RleftSQ*=RleftSQ;
			RrightSQ*=RrightSQ;
			RSQ=Ri*Ri;
			RsmSQ=Rsm*Rsm;
			
			for(k=0;k<3;k++)
			{
				rf[k]=r[k][i]+Displ[k];
				rint[k]=(int)(rf[k]+0.5);
				start[k]=rint[k]-iRsm;
				end[k]=rint[k]+iRsm;
			}
			bool AtomReachable=true;
			for(k=0;k<3;k++)
			{
				if(end[k]<0)AtomReachable=false;
				if(start[k]>GridSize[k]-1)AtomReachable=false;
			}
			if(!AtomReachable)
			{
			//DbgPrint0("Atom Not Here\n");
				continue;
			}
		
			for(k=0;k<3;k++)
			{
				if(start[k]<0)start[k]=0;
				if(end[k]>GridSize[k]-1)end[k]=GridSize[k]-1;
			}
			
			for(ix=start[0];ix<=end[0];ix++)
				for(iy=start[1];iy<=end[1];iy++)
					for(iz=start[2];iz<=end[2];iz++)
			{
				gridpoint=ix+iy*GSX+iz*GSXY;
				RtmpSQ=(rf[0]-ix)*(rf[0]-ix)+(rf[1]-iy)*(rf[1]-iy)+(rf[2]-iz)*(rf[2]-iz);
				vtmp=Surf+4*gridpoint;
				if(RtmpSQ<=RSQ)//RtmpSQ<=RSQ
				{
					Field[gridpoint]=iValue;
					vtmp[0]=-BIGDISTANSE;
					vtmp[3]=0.0;
				}
				else//RSQ<RtmpSQ
				{
					Rtmp=sqrt(RtmpSQ);
					if(RtmpSQ<=RsmSQ)//RSQ<RtmpSQ<RleftSQ
					{
						if(Field[gridpoint]==iBulkValue)
						{
							vtmp[0]=-BIGDISTANSE;
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp-Ri;
						}
						else if(Field[gridpoint]<0&&Rtmp-Ri<vtmp[3])
						{
							vtmp[0]=-BIGDISTANSE;
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp-Ri;
						}
					}
					if(RleftSQ<=RtmpSQ&&RtmpSQ<=RrightSQ)//RleftSQ<=RtmpSQ<=RrightSQ
					{
						if(vtmp[0]>-100.0f)//e.i. is intersection
						{
							if(Rtmp-Ri<vtmp[3])
							{
								ftmp=Rsm/Rtmp;
								vtmp[0]=rf[0]+ftmp*(ix-rf[0]);
								vtmp[1]=rf[1]+ftmp*(iy-rf[1]);
								vtmp[2]=rf[2]+ftmp*(iz-rf[2]);
								vtmp[3]=Rtmp-Ri;
							}
						}
						else if(Field[gridpoint]==iBulkValue)
						{
							ftmp=Rsm/Rtmp;
							vtmp[0]=rf[0]+ftmp*(ix-rf[0]);
							vtmp[1]=rf[1]+ftmp*(iy-rf[1]);
							vtmp[2]=rf[2]+ftmp*(iz-rf[2]);
							vtmp[3]=Rtmp-Ri;
						}
					}
				}
			}
		}
		/*
		if(MakePreRoll)
		{
			
			DbgPrint0("MakePreRoll\n");
			int lim=(int)(Rsmoth+0.5);
			int itmp=GSXYZ*4;
			int iSolAcFlag=1000000000;
			RsmSQ=Rsmoth*Rsmoth;
			for(i=0;i<itmp;i=i+4)
			{
				if(Surf[i]>-100.0)
				{
					vtmp=Surf+i;
					for(k=0;k<3;k++)
					{
						rint[k]=(int)(vtmp[k]+0.5);
						start[k]=rint[k]-lim;
						end[k]=rint[k]+lim;
						if(start[k]<0)start[k]=0;
						if(end[k]>GridSize[k]-1)end[k]=GridSize[k]-1;
					}
					for(ix=start[0];ix<=end[0];ix++)
						for(iy=start[1];iy<=end[1];iy++)
							for(iz=start[2];iz<=end[2];iz++)
					{
						gridpoint=ix+iy*GSX+iz*GSXY;
						RtmpSQ=(vtmp[0]-ix)*(vtmp[0]-ix)+(vtmp[1]-iy)*(vtmp[1]-iy)+(vtmp[2]-iz)*(vtmp[2]-iz);
						if(RtmpSQ<RsmSQ&&Field[gridpoint]<0)
						{
							Field[gridpoint]=iSolAcFlag;
						}
					}
				}
			}
			for(i=0;i<GSXYZ;i++)
			{
				if(Field[i]<0)
				{
					Field[i]=-Field[i];
				}
				else if(Field[i]==iSolAcFlag)
				{
					Field[i]=-iValue;
				}
			}
		}
		*/
	}
	return EXIT_SUCCESS;
}
/*{
	DbgPrint2("GOAtoms::BuildPreMapsMemLim(iValue=%d iBulkValue=%d)\n"
			,iValue,iBulkValue);
	int i,j,k,gridpoint,gridPoint2;
	float RSQ,RsmSQ,RtmpSQ,Ri,Rsm,Rtmp;
	float RleftSQ,RrightSQ;
	float rf[3];
	int iR,iRsm;
	int ix,iy,iz;
	float ftmp;
	float *vtmp;
	float *Surf=pw->Surf;
	
	
	
	int GridSize[3]={pw->GS_X,pw->GS_Y,pw->locGS_Z};
	int GS_X=pw->GS_X;
	int GS_XY=pw->GS_X*pw->GS_Y;
	int GS_XYZ=pw->GS_X*pw->GS_Y*pw->locGS_Z;
	
	int *Field=pw->Field+GS_XY*pw->locR0_Z;
	//int FieldShift=GS_XY*pw->locR0_Z;
	
	int start[3];
	int end[3];
	int rint[3];
	
	for(i=0;i<NAtoms;i++)
	{
		
		Ri=R[i]+Rion;
		if(Ri==0.0){continue;}
			
			
		Rsm=Ri+Rsmoth;
		iR=(int)(Ri+0.5);
		iRsm=(int)(Rsm+1.0);
		RleftSQ=Rsm-0.707106781f;
		RrightSQ=Rsm+0.707106781f;
		RleftSQ*=RleftSQ;
		RrightSQ*=RrightSQ;
		RSQ=Ri*Ri;
		RsmSQ=Rsm*Rsm;
			
		for(k=0;k<3;k++)
		{
			rf[k]=r[k][i]+Displ[k];
			rint[k]=(int)(rf[k]+0.5);
			start[k]=rint[k]-iRsm;
			end[k]=rint[k]+iRsm;
		}
		//DbgPrint0("rf[%d]=%f %f %f\n",i,rf[0],rf[1],rf[2]);
		//DbgPrint0("start[%d]=%d %d %d\n",i,start[0],start[1],start[2]);
		//DbgPrint0("end[%d]=%d %d %d\n",i,end[0],end[1],end[2]);
		bool AtomReachable=true;
		for(k=0;k<3;k++)
		{
			if(end[k]<0)AtomReachable=false;
			if(start[k]>GridSize[k]-1)AtomReachable=false;
		}
		if(!AtomReachable)
		{
			//DbgPrint0("Atom Not Here\n");
			continue;
		}
		
		for(k=0;k<3;k++)
		{
			if(start[k]<0)start[k]=0;
			if(end[k]>GridSize[k]-1)end[k]=GridSize[k]-1;
		}
		//DbgPrint0("start[%d]=%d %d %d\n",i,start[0],start[1],start[2]);
		//DbgPrint0("end[%d]=%d %d %d\n",i,end[0],end[1],end[2]);
		for(ix=start[0];ix<=end[0];ix++)
			for(iy=start[1];iy<=end[1];iy++)
				for(iz=start[2];iz<=end[2];iz++)
		{
			gridpoint=ix+iy*GS_X+iz*GS_XY;
			RtmpSQ=(rf[0]-ix)*(rf[0]-ix)+(rf[1]-iy)*(rf[1]-iy)+(rf[2]-iz)*(rf[2]-iz);
			vtmp=Surf+4*gridpoint;
			if(RtmpSQ<=RSQ)//RtmpSQ<=RSQ
			{
				Field[gridpoint]=iValue;
				vtmp[0]=-BIGDISTANSE;
				vtmp[3]=0.0;
			}
			else//RSQ<RtmpSQ
			{
				Rtmp=sqrt(RtmpSQ);
				if(RtmpSQ<=RsmSQ)//RSQ<RtmpSQ<RleftSQ
				{
					if(Field[gridpoint]==iBulkValue)
					{
						vtmp[0]=-BIGDISTANSE;
						Field[gridpoint]=-iValue;
						vtmp[3]=Rtmp-Ri;
					}
					else if(Field[gridpoint]<0&&Rtmp-Ri<vtmp[3])
					{
						vtmp[0]=-BIGDISTANSE;
						Field[gridpoint]=-iValue;
						vtmp[3]=Rtmp-Ri;
					}
				}
				if(RleftSQ<=RtmpSQ&&RtmpSQ<=RrightSQ)//RleftSQ<=RtmpSQ<=RrightSQ
				{
					if(vtmp[0]>-100.0f)//e.i. is intersection
					{
						if(Rtmp-Ri<vtmp[3])
						{
							ftmp=Rsm/Rtmp;
							vtmp[0]=rf[0]+ftmp*(ix-rf[0]);
							vtmp[1]=rf[1]+ftmp*(iy-rf[1]);
							vtmp[2]=rf[2]+ftmp*(iz-rf[2]);
							vtmp[3]=Rtmp-Ri;
						}
					}
					else if(Field[gridpoint]==iBulkValue)
					{
						ftmp=Rsm/Rtmp;
						vtmp[0]=rf[0]+ftmp*(ix-rf[0]);
						vtmp[1]=rf[1]+ftmp*(iy-rf[1]);
						vtmp[2]=rf[2]+ftmp*(iz-rf[2]);
						vtmp[3]=Rtmp-Ri;
					}
				}
			}
		}
	}
	return EXIT_SUCCESS;
}*/
int GOAtoms::BuildLJRepultionMap(BuildWorldNI *Builder,ContWorld* world,int Ion,float *V,float LimitVlj)
{
	PNP_EXIT_FAIL_NULL(Builder,"BuildWorldNI==NULL\n");
	PNP_EXIT_FAIL_NULL(world,"ContWorld==NULL\n");
	
	DbgPrint2("GOAtoms::BuildLJRepultionMap\n");
	int i,j,k,gridpoint,rint[3],gridPoint2;
	int GridSizeXY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
	float RSQ,RsmSQ,RtmpSQ,Ri,Rsm,Rtmp;
	float RleftSQ,RrightSQ;
	float rf[3];
	int iR,iRsm;
	int ix,iy,iz;
	float ftmp;
	float *vtmp;
	
	int *GridSize=world->GridSizeGlobal;
	double GS2=world->GridScale*world->GridScale;
	double GS6=GS2*GS2*GS2;
	double GS12=GS6*GS6;
	int GSX=GridSize[0];
	int GSXY=GridSize[0]*GridSize[1];
	int GSXYZ=GridSize[0]*GridSize[1]*GridSize[2];
	int start[3];
	
	int end[3];
	
	if(IAVMethod==2)
	{
		for(i=0;i<NAtoms;i++)
		{
			double Sigma;
			double LimitVljSigma;
			double A=LJA[Ion][i];
			double B=LJB[Ion][i];
			if(A==0.0||B==0.0)
				Sigma=0.0;
			else
				Sigma=pow(A/B,1.0/6.0);
			LimitVljSigma=LimitVlj*Sigma;
			
			//DbgPrint0("atm %d A [%e %e] B[%e %e] sigma [%g %g]\n",i,A,A/GS12,B,B/GS6,Sigma,Sigma/world->GridScale);
			
			float LimitVljSigmaSQ=LimitVljSigma*LimitVljSigma;
			
			double R6,R12;
			
			if(A==0.0||B==0.0)
				continue;
			Ri=LimitVljSigma;
			iR=(int)(Ri+0.5);
			RSQ=Ri*Ri;
			for(k=0;k<3;k++){
				rf[k]=r[k][i];
				rint[k]=(int)(rf[k]+0.5);
				start[k]=rint[k]-iR;
				end[k]=rint[k]+iR;
				if(start[k]<0)start[k]=0;
				if(end[k]>GridSize[k]-1)start[k]=GridSize[k]-1;
			}
				
			for(ix=start[0];ix<=end[0];ix++)
				for(iy=start[1];iy<=end[1];iy++)
					for(iz=start[2];iz<=end[2];iz++)
			{
				gridpoint=ix+iy*GSX+iz*GSXY;
				RtmpSQ=(rf[0]-ix)*(rf[0]-ix)+(rf[1]-iy)*(rf[1]-iy)+(rf[2]-iz)*(rf[2]-iz);
				if(RtmpSQ<LimitVljSigmaSQ)
				{
					R6=RtmpSQ*RtmpSQ*RtmpSQ;
					R12=R6*R6;
					
					V[gridpoint]+=(A/R12)-(B/R6);
#if defined(_MSC_VER) || defined(__DECCXX) || (__GNUC__ < 3)
					if(V[gridpoint]>1.0e2||RtmpSQ<1.0)
#else
					if(V[gridpoint]>1.0e2||RtmpSQ<1.0||V[gridpoint]==NAN||V[gridpoint]==INFINITY)
#endif
						V[gridpoint]=1.0e2;
				}
			}
		}
	}
	else if(IAVMethod==3)//i.e. AN
	{
		for(i=0;i<NAtoms;i++)
		{
			double Sigma;
			double LimitVljSigma;
			double A=SRA[Ion][i];
			double N=SRN[Ion][i];
			if(A==0.0||N==0.0)
				Sigma=0.0;
			else
				Sigma=A;
			LimitVljSigma=Sigma;
			
			//DbgPrint0("atm %d A [%e %e] N[%e] sigma [%g %g]\n",i,A,A/world->GridScale,N,Sigma,Sigma/world->GridScale);
			
			float LimitVljSigmaSQ=LimitVljSigma*LimitVljSigma;
			
			double R6,R12;
			
			if(A==0.0||N==0.0)
				continue;
			Ri=LimitVljSigma;
			iR=(int)(Ri+0.5);
			RSQ=Ri*Ri;
			for(k=0;k<3;k++){
				rf[k]=r[k][i];
				rint[k]=(int)(rf[k]+0.5);
				start[k]=rint[k]-iR;
				end[k]=rint[k]+iR;
				if(start[k]<0)start[k]=0;
				if(end[k]>GridSize[k]-1)start[k]=GridSize[k]-1;
			}
				
			for(ix=start[0];ix<=end[0];ix++)
				for(iy=start[1];iy<=end[1];iy++)
					for(iz=start[2];iz<=end[2];iz++)
			{
				gridpoint=ix+iy*GSX+iz*GSXY;
				RtmpSQ=(rf[0]-ix)*(rf[0]-ix)+(rf[1]-iy)*(rf[1]-iy)+(rf[2]-iz)*(rf[2]-iz);
				if(RtmpSQ<LimitVljSigmaSQ)
				{
					Rtmp=sqrt(RtmpSQ);
					
					V[gridpoint]+=pow(A/Rtmp,N);
#if defined(_MSC_VER) || defined(__DECCXX) || (__GNUC__ < 3)
					if(V[gridpoint]>1.0e2||RtmpSQ<1.0)
#else
						if(V[gridpoint]>1.0e2||RtmpSQ<1.0||V[gridpoint]==NAN||V[gridpoint]==INFINITY)
#endif
							V[gridpoint]=1.0e2;
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
#if defined(WITH_CUDA)
extern "C" int GOAtoms_SetCoulombicBC(int *GridSize,int *BldBCatPlane,float4 *rq,int Nq,float eps,float *Potential);
#endif
int GOAtoms::SetBoundaryConditionOnCuda(BuildWorldNI *Builder,ContWorld* world)
{
	PNP_EXIT_FAIL_NULL(Builder,"BuildWorldNI==NULL\n");
	PNP_EXIT_FAIL_NULL(world,"ContWorld==NULL\n");
	
	if(Builder->BoundaryCondition==BuildWorldNI::CoulBC)
	{
		int i,j,k,natom,iq;
		int GrPnt;
		int GS_X=world->GridSize[0];
		int GS_Y=world->GridSize[0];
		int GS_Z=world->GridSize[0];
		int GS_XY=world->GridSize[0]*world->GridSize[1];
		float r1;
		float eps;
		float fpoh= 4*M_PI*world->GridScale;
		float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
		
		pnpPrint("Setting up CoulombicBC\n");
		
		int iStart,jStart,kStart;
		int iEnd,jEnd,kEnd;
		float kDH=0.0;
		if(Builder->MakeConcentrationMap==true)//i.e for PB or PNP
		{
			float I;
			for(i=0;i<NodeIndexMaxValues;i++)
				if(world->NIndexing->C[i]>0.0)
			{
				I=world->NIndexing->C[i];
				break;
			}
			float IM=I/coef;
			float A = 2.0 * I/world->NIndexing->Eps[ (world->NIndexing->NIndex[0]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft ];
			kDH=sqrt(A);
		}
		if(kDH==0.0)
		{
#if defined(WITH_CUDA)
			//Here will do with cuda
			float4 *rq=new float4[NAtoms];
			int Nq=0;
			for(j=0;j<NAtoms;j++)
			{
				if(q[j]!=0.0f)
				{
					rq[Nq].x=r[0][j];
					rq[Nq].y=r[1][j];
					rq[Nq].z=r[2][j];
					rq[Nq].w=q[j];
					Nq++;
				}
			}
			
			
			/*
			if(Builder->BldBCatPlane[2])
			{
			for(i=0;i<GS_X;i++)
			for(j=0;j<GS_Y;j++)//XY
			{
			k=0;
			GrPnt=i+j*GS_X+k*GS_XY;
			eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
					
			for(iq=0;iq < Nq;iq++)
			{
			r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
			world->Potential[GrPnt]+=rq[iq].w/(eps*r1);
		}
					
			k=GS_Z-1;
			GrPnt=i+j*GS_X+k*GS_XY;
			eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
					
			for(iq=0;iq < Nq;iq++)
			{
			r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
			world->Potential[GrPnt]+=rq[iq].w/(eps*r1);
		}
		}
		}
			if(Builder->BldBCatPlane[0])
			{
			jStart=0;
			kStart=0;
			jEnd=GS_Y;
			kEnd=GS_Z;
			if(Builder->BldBCatPlane[2])
			{
			kStart=1;
			kEnd=GS_Z-1;
		}
			for(j=jStart;j<jEnd;j++)//YZ
			for(k=kStart;k<kEnd;k++)
			{
			i=0;
			GrPnt=i+j*GS_X+k*GS_XY;
			eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
					
			for(iq=0;iq < Nq;iq++)
			{
			r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
			world->Potential[GrPnt]+=rq[iq].w/(eps*r1);
		}
					
			i=GS_X-1;
			GrPnt=i+j*GS_X+k*GS_XY;
			eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
					
			for(iq=0;iq < Nq;iq++)
			{
			r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
			world->Potential[GrPnt]+=rq[iq].w/(eps*r1);
		}
		}
		}
			if(Builder->BldBCatPlane[1])
			{
			iStart=0;
			kStart=0;
			iEnd=GS_X;
			kEnd=GS_Z;
			if(Builder->BldBCatPlane[2])
			{
			kStart=1;
			kEnd=GS_Z-1;
		}
			if(Builder->BldBCatPlane[0])
			{
			iStart=1;
			iEnd=GS_X-1;
		}
			for(i=iStart;i<iEnd;i++)//ZX
			for(k=kStart;k<kEnd;k++)
			{
			j=0;
			GrPnt=i+j*GS_X+k*GS_XY;
			eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
					
			for(iq=0;iq < Nq;iq++)
			{
			r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
			world->Potential[GrPnt]+=rq[iq].w/(eps*r1);
		}
					
			j=GS_Y-1;
			GrPnt=i+j*GS_X+k*GS_XY;
			eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
					
			for(iq=0;iq < Nq;iq++)
			{
			r1=sqrt((i-rq[iq].x)*(i-rq[iq].x) + (j-rq[iq].y)*(j-rq[iq].y) + (k-rq[iq].z)*(k-rq[iq].z));
			world->Potential[GrPnt]+=rq[iq].w/(eps*r1);
		}
		}
		}
*/
			float t_eps=world->NIndexing->GetDielFloat(0,0)*4.0*M_PI;
			int t_BldBCatPlane[3]={Builder->BldBCatPlane[0],Builder->BldBCatPlane[1],Builder->BldBCatPlane[2]};
 
			GOAtoms_SetCoulombicBC(world->GridSize,t_BldBCatPlane, rq, Nq,t_eps,world->Potential);
 
 			DeleteCArray(rq);
#endif
		}
		else
		{
			if(Builder->BldBCatPlane[2])
			{
				for(i=0;i<GS_X;i++)
					for(j=0;j<GS_Y;j++)//XY
				{
					k=0;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
						
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]*exp(-kDH*r1)/(eps*r1);
					}
						
					k=GS_Z-1;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
						
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]*exp(-kDH*r1)/(eps*r1);
					}
				}
			}
			if(Builder->BldBCatPlane[0])
			{
				jStart=0;
				kStart=0;
				jEnd=GS_Y;
				kEnd=GS_Z;
				if(Builder->BldBCatPlane[2])
				{
					kStart=1;
					kEnd=GS_Z-1;
				}
				for(j=jStart;j<jEnd;j++)//YZ
					for(k=kStart;k<kEnd;k++)
				{
					i=0;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
						
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]*exp(-kDH*r1)/(eps*r1);
					}
						
					i=GS_X-1;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
						
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]*exp(-kDH*r1)/(eps*r1);
					}
				}
			}
			if(Builder->BldBCatPlane[1])
			{
				iStart=0;
				kStart=0;
				iEnd=GS_X;
				kEnd=GS_Z;
				if(Builder->BldBCatPlane[2])
				{
					kStart=1;
					kEnd=GS_Z-1;
				}
				if(Builder->BldBCatPlane[0])
				{
					iStart=1;
					iEnd=GS_X-1;
				}
				for(i=iStart;i<iEnd;i++)//ZX
					for(k=kStart;k<kEnd;k++)
				{
					j=0;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
						
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]*exp(-kDH*r1)/(eps*r1);
					}
						
					j=GS_Y-1;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
						
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]*exp(-kDH*r1)/(eps*r1);
					}
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
int GOAtoms::SetBoundaryCondition(BuildWorldNI *Builder,ContWorld* world)
{
	PNP_EXIT_FAIL_NULL(Builder,"BuildWorldNI==NULL\n");
	PNP_EXIT_FAIL_NULL(world,"ContWorld==NULL\n");
	
	if(Builder->BoundaryCondition==BuildWorldNI::CoulBC)
	{
		int i,j,k,natom;
		int GrPnt;
		int GS_X=world->GridSize[0];
		int GS_Y=world->GridSize[1];
		int GS_Z=world->GridSize[2];
		int GS_XY=world->GridSize[0]*world->GridSize[1];
		float r1;
		float eps;
		float fpoh= 4*M_PI*world->GridScale;
		float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
		
		pnpPrint("Setting up CoulombicBC\n");
		
		int iStart,jStart,kStart;
		int iEnd,jEnd,kEnd;
		float kDH=0.0;
		if(Builder->MakeConcentrationMap==true)//i.e for PB or PNP
		{
			float I;
			I=Builder->BulkParam->C[0]*coef;
			/*for(i=0;i<NodeIndexMaxValues;i++)
				if(world->NIndexing->C[i]>0.0)
			{
				I=world->NIndexing->C[i];
				break;
			}*/
			float IM=I/coef;
			float A = 2.0 * I/world->NIndexing->Eps[ (world->NIndexing->NIndex[0]&NodeIndexing::Epsilon0)>>NodeIndexing::Epsilon0Sft ];
			kDH=sqrt(A);
		}
		if(kDH==0.0)
		{
			
			if(Builder->BldBCatPlane[2])
			{
				for(i=1;i<GS_X;i++)
					for(j=1;j<GS_Y;j++)//XY
				{
					k=0;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
					
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]/(eps*r1);
					}
					
					k=GS_Z-1;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
					
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]/(eps*r1);
					}
				}
			}
			if(Builder->BldBCatPlane[0])
			{
				jStart=0;
				kStart=0;
				jEnd=GS_Y;
				kEnd=GS_Z;
				if(Builder->BldBCatPlane[2])
				{
					kStart=1;
					kEnd=GS_Z-1;
				}
				for(j=jStart;j<jEnd;j++)//YZ
					for(k=kStart;k<kEnd;k++)
				{
					i=0;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
					
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]/(eps*r1);
					}
					
					i=GS_X-1;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
					
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]/(eps*r1);
					}
				}
			}
			if(Builder->BldBCatPlane[1])
			{
				iStart=0;
				kStart=0;
				iEnd=GS_X;
				kEnd=GS_Z;
				if(Builder->BldBCatPlane[2])
				{
					kStart=1;
					kEnd=GS_Z-1;
				}
				if(Builder->BldBCatPlane[0])
				{
					iStart=1;
					iEnd=GS_X-1;
				}
				for(i=iStart;i<iEnd;i++)//ZX
					for(k=kStart;k<kEnd;k++)
				{
					j=0;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
					
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]/(eps*r1);
					}
					
					j=GS_Y-1;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
					
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]/(eps*r1);
					}
				}
			}
		}
		else
		{
			if(Builder->BldBCatPlane[2])
			{
				for(i=0;i<GS_X;i++)
					for(j=0;j<GS_Y;j++)//XY
				{
					k=0;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
						
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]*exp(-kDH*r1)/(eps*r1);
					}
						
					k=GS_Z-1;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
						
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]*exp(-kDH*r1)/(eps*r1);
					}
				}
			}
			if(Builder->BldBCatPlane[0])
			{
				jStart=0;
				kStart=0;
				jEnd=GS_Y;
				kEnd=GS_Z;
				if(Builder->BldBCatPlane[2])
				{
					kStart=1;
					kEnd=GS_Z-1;
				}
				for(j=jStart;j<jEnd;j++)//YZ
					for(k=kStart;k<kEnd;k++)
				{
					i=0;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
						
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]*exp(-kDH*r1)/(eps*r1);
					}
						
					i=GS_X-1;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
						
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]*exp(-kDH*r1)/(eps*r1);
					}
				}
			}
			if(Builder->BldBCatPlane[1])
			{
				iStart=0;
				kStart=0;
				iEnd=GS_X;
				kEnd=GS_Z;
				if(Builder->BldBCatPlane[2])
				{
					kStart=1;
					kEnd=GS_Z-1;
				}
				if(Builder->BldBCatPlane[0])
				{
					iStart=1;
					iEnd=GS_X-1;
				}
				for(i=iStart;i<iEnd;i++)//ZX
					for(k=kStart;k<kEnd;k++)
				{
					j=0;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
						
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]*exp(-kDH*r1)/(eps*r1);
					}
						
					j=GS_Y-1;
					GrPnt=i+j*GS_X+k*GS_XY;
					eps=world->NIndexing->GetDielFloat(0,GrPnt)*4*M_PI;
						
					for(natom=0;natom < NAtoms;natom++)
					{
						r1=sqrt((i-r[0][natom])*(i-r[0][natom]) + (j-r[1][natom])*(j-r[1][natom]) + (k-r[2][natom])*(k-r[2][natom]));
						world->Potential[GrPnt]+=q[natom]*exp(-kDH*r1)/(eps*r1);
					}
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
int GOAtoms::BuildPreMapsCharges(BuildWorldNI *Builder,ContWorld* world,float *Field)
{
	if(ChargeDist==ClosestNode)
	{
		BuildPreMapsChargesClosestNode(Builder,world,Field);
	}
	else if(ChargeDist==Linear)
	{
		if(ChargeDistN==1)
			BuildPreMapsChargesLinear1(Builder,world,Field);
		else
			BuildPreMapsChargesLinearN(Builder,world,Field,ChargeDistN);
	}
	else if(ChargeDist==Cone)
	{
		BuildPreMapsChargesCone(Builder,world,Field,ChargeDistN);
	}
	return EXIT_SUCCESS;
}
int GOAtoms::BuildPreMapsChargesClosestNode(BuildWorldNI *Builder,ContWorld* world,float *Field)
{
	int i,j,GridPoint;
	int GSX=world->GridSizeGlobal[0];
	int GSY=world->GridSizeGlobal[1];
	int GSZ=world->GridSizeGlobal[2];
	int GSXY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
	int GSXYZ = world->GridSizeGlobal[0] * world->GridSizeGlobal[1] * world->GridSizeGlobal[2];
	int ir[3];
	
	InfoPrint("Charge Distrebution Method is Single Grid\n");
	for(i=0;i<NAtoms;i++)
	{
		for(j=0;j<3;j++)
			ir[j]=(int)roundf(r[j][i]);
		GridPoint=ir[0]+ir[1]*GSX+ir[2]*GSXY;
		if(ir[0]<0||ir[1]<0||ir[2]<0||ir[0]>GSX-1||ir[1]>GSY-1||ir[2]>GSZ-1)
		{
			pnpWarning("Atom %d is not fitted in box\n",i);
		}
		else
		{
			Field[GridPoint]+=q[i];
		}
	}
	
	return EXIT_SUCCESS;
}
int GOAtoms::BuildPreMapsChargesLinear1(BuildWorldNI *Builder,ContWorld* world,float *Field)
{
	int i,j,GridPoint,gr[3];
	int GSX=world->GridSizeGlobal[0];
	int GSY=world->GridSizeGlobal[1];
	int GSZ=world->GridSizeGlobal[2];
	int GSXY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
	int GSXYZ=world->GridSizeGlobal[0]*world->GridSizeGlobal[1]*world->GridSizeGlobal[2];
	int ir[3];
	float Q=0;
	float f,F;
	float rgrid[3];

	InfoPrint("Charge Distrebution Method is Quick And Simple 8 Grids (Klapper, I., Hagstrom, R., Fine, R., Sharp, K., Honig, B. (1986). Focusing of electric fields in the active site of Cu-Zn Superoxide Dismutase: Effects of ionic strength and amino-acid modification. Proteins 1, p 47.)\n");
	int countA=0,countC=0;
	for(i=0;i<NAtoms;i++)
	{
		for(j=0;j<3;j++)ir[j]=(int)r[j][i];
		GridPoint=ir[0]+ir[1]*GSX+ir[2]*GSXY;
		
		if(ir[0]<0||ir[1]<0||ir[2]<0||ir[0]>GSX-1||ir[1]>GSY-1||ir[2]>GSZ-1)
		{
			if(countA<200)
				pnpWarning("Atom %d is not fitted in box\n",i);
			countA++;
		}
		else
		{
			F=0.0;
			for(gr[0]=ir[0];gr[0]<ir[0]+2;gr[0]++)
				for(gr[1]=ir[1];gr[1]<ir[1]+2;gr[1]++)
					for(gr[2]=ir[2];gr[2]<ir[2]+2;gr[2]++)
			{
				GridPoint=gr[0]+gr[1]*GSX+gr[2]*GSXY;
				if(gr[0]>0&&gr[1]>0&&gr[2]>0&&gr[0]<GSX-1&&gr[1]<GSY-1&&gr[2]<GSZ-1)
				{
					for(j=0;j<3;j++)
					{
						rgrid[j]=float(gr[j])-r[j][i];
						rgrid[j]=fabs(rgrid[j]);
					}
					f=(1.0-rgrid[0])*(1.0-rgrid[1])*(1.0-rgrid[2]);
					Field[GridPoint]+=q[i]*f;
					F+=f;
					Q+=q[i]*f;
				}
				else
				{
					countC++;
				}
			}
		}
	}
	if(countA>0)
		pnpWarning("Totally %d atoms was outside of the box\n",countA);
	if(countC>0)
		pnpWarning("%d partial charge points was out of the box\n",countC);
	return EXIT_SUCCESS;
}
int GOAtoms::BuildPreMapsChargesLinearN(BuildWorldNI *Builder,ContWorld* world,float *Field, int N)
{
	int i,j,GridPoint,gr[3];
	int GSX=world->GridSizeGlobal[0];
	int GSY=world->GridSizeGlobal[1];
	int GSZ=world->GridSizeGlobal[2];
	int GSXY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
	int GSXYZ=world->GridSizeGlobal[0]*world->GridSizeGlobal[1]*world->GridSizeGlobal[2];
	int ir[3];
	float Q=0;
	float f,F;
	float rgrid[3];
	float fN=float(N);
	float fM=1.0/(fN*fN*fN*fN*fN*fN);

	InfoPrint("Trilinear interpolation with interger radius %d, i.e. total nodes in each demention is %d or %d\n",N, 2*N-1, 2*N);
	
	for(i=0;i<NAtoms;i++)
	{
		for(j=0;j<3;j++)ir[j]=(int)r[j][i];
		GridPoint=ir[0]+ir[1]*GSX+ir[2]*GSXY;
		
		if(ir[0]<0||ir[1]<0||ir[2]<0||ir[0]>GSX-1||ir[1]>GSY-1||ir[2]>GSZ-1)
		{
			pnpWarning("Atom %d is not fitted in box\n",i);
		}
		else
		{
			F=0.0;
			for(gr[0]=ir[0]-N+1;gr[0]<ir[0]+N+1;gr[0]++)
				for(gr[1]=ir[1]-N+1;gr[1]<ir[1]+N+1;gr[1]++)
					for(gr[2]=ir[2]-N+1;gr[2]<ir[2]+N+1;gr[2]++)
			{
				GridPoint=gr[0]+gr[1]*GSX+gr[2]*GSXY;
				for(j=0;j<3;j++)
				{
					rgrid[j]=(float)gr[j]-r[j][i];
					rgrid[j]=fabs(rgrid[j]);
				}
				f=(fN-rgrid[0])*(fN-rgrid[1])*(fN-rgrid[2])*fM;
				Field[GridPoint]+=q[i]*f;
				F+=f;
				Q+=q[i]*f;
			}
		}
	}
	return EXIT_SUCCESS;
}
int GOAtoms::BuildPreMapsChargesCone(BuildWorldNI *Builder,ContWorld* world,float *Field, int N)
{
	int i,j,GridPoint,gr[3],grLoc[3],GridPointLoc;
	int GSX=world->GridSizeGlobal[0];
	int GSY=world->GridSizeGlobal[1];
	int GSZ=world->GridSizeGlobal[2];
	int GSXY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
	int GSXYZ=world->GridSizeGlobal[0]*world->GridSizeGlobal[1]*world->GridSizeGlobal[2];
	int ir[3];
	float Q=0;
	float f,F;
	float rgrid[3];
	float fN=float(N);
	float fM=1.0/(fN*fN*fN*fN*fN*fN);

	float fR,fRSQ;
	InfoPrint("cone interpolation with interger radius %d, i.e. total nodes in each demention is %d or %d\n",N, 2*N-1, 2*N);
	
	float *ChDist=new float[2*N*2*N*2*N];
	
	for(i=0;i<NAtoms;i++)
	{
		for(j=0;j<2*N*2*N*2*N;j++)ChDist[j]=0.0;
		
		for(j=0;j<3;j++)ir[j]=(int)r[j][i];
		GridPoint=ir[0]+ir[1]*GSX+ir[2]*GSXY;
		
		if(ir[0]<0||ir[1]<0||ir[2]<0||ir[0]>GSX-1||ir[1]>GSY-1||ir[2]>GSZ-1)
		{
			pnpWarning("Atom %d is not fitted in box\n",i);
		}
		else
		{
			F=0.0;
			for(grLoc[0]=0;grLoc[0]<2*N;grLoc[0]++)
				for(grLoc[1]=0;grLoc[1]<2*N;grLoc[1]++)
					for(grLoc[2]=0;grLoc[2]<2*N;grLoc[2]++)
			{
				GridPointLoc=grLoc[0]+grLoc[1]*2*N+gr[2]*4*N*N;
				
				gr[0]=ir[0]-N+1+grLoc[0];
				gr[1]=ir[1]-N+1+grLoc[1];
				gr[2]=ir[2]-N+1+grLoc[2];				
				GridPoint=gr[0]+gr[1]*GSX+gr[2]*GSXY;
				
				fRSQ=(float(gr[0])-r[0][i])*(float(gr[0])-r[0][i])+(float(gr[1])-r[1][i])*(float(gr[1])-r[1][i])+(float(gr[2])-r[2][i])*(float(gr[2])-r[2][i]);
				fR=sqrt(fRSQ);
				f=fN-fR;
				if(f<0.0)f=0.0;
				
				ChDist[GridPointLoc]=f;
				F+=f;
				
			}
			for(grLoc[0]=0;grLoc[0]<2*N;grLoc[0]++)
				for(grLoc[1]=0;grLoc[1]<2*N;grLoc[1]++)
					for(grLoc[2]=0;grLoc[2]<2*N;grLoc[2]++)
			{
				GridPointLoc=grLoc[0]+grLoc[1]*2*N+gr[2]*4*N*N;
				
				gr[0]=ir[0]-N+1+grLoc[0];
				gr[1]=ir[1]-N+1+grLoc[1];
				gr[2]=ir[2]-N+1+grLoc[2];				
				GridPoint=gr[0]+gr[1]*GSX+gr[2]*GSXY;
				
				//ChDist[GridPointLoc]=f;
				Field[GridPoint]+=q[i]*ChDist[GridPointLoc]/F;
				//Q+=q[i]*ChDist[GridPointLoc]/F;
				
			}
		}
	}
	delete [] ChDist;
	return EXIT_SUCCESS;
}
int GOAtoms::RotateGGO(double *n, double cosa, double sina)
{
	int i,j;
	double rtemp[3];
	
	for(i=0;i<NAtoms;i++)
	{
		for(j=0;j<3;j++)
			rtemp[j]=r[j][i];
		RotateVecDouble(rtemp,n,cosa,sina);
		for(j=0;j<3;j++)
			r[j][i]=rtemp[j];
	}
	return EXIT_SUCCESS;
}
int GOAtoms::SetGridIDs(int startIDfrom)
{
	ID_GridBegin=startIDfrom;
	ID_GridEnd=startIDfrom+NAtoms-1;
	return ID_GridEnd+1;
}
///////////////////////////////////////////////////////////////////////////////
GOTube::GOTube()
{
	InitZero();
}
GOTube::~GOTube()
{
	Clear();
}
int GOTube::InitZero()
{
	GenericGeometricalObject::InitZero();
	PnpsObject::SetName("TubeParameters");
	HalfSigma=0.0;
	FourEpsilon=0.0;
	return EXIT_SUCCESS;
}
int GOTube::Clear()
{
	GenericGeometricalObject::Clear();
	return EXIT_SUCCESS;
}
int GOTube::Print(BuildWorldNI *buildWorld)
{
	GenericGeometricalObject::Print(buildWorld);
	pnpPrintGroup0("BW        Centeral Axis x,y Coordinates:..... [%.3f A, %.3f A]\n",XY[0],XY[1]);
	pnpPrintGroup0("BW        Tube Z Boarders:................... [%.3f A, %.3f A]\n",Z[0],Z[1]);
	pnpPrintGroup0("BW        Tube Radii:........................ [%.3f A, %.3f A]\n",R[0],R[1]);
	return EXIT_SUCCESS;
}
int GOTube::setTubeParam(float X,float Y,float Z0,float Z1,float R0,float R1)
{
	XY[0]=X;
	XY[1]=Y;
	Z[0]=Z0;
	Z[1]=Z1;
	R[0]=R0;
	R[1]=R1;
	return EXIT_SUCCESS;
}
int GOTube::ChangeUnits(ContWorld* world,bool ToInternal)
{
	if(world==NULL)
	{
		fprintf(stderr,"ERROR: ChangeUnits: world==NULL\n");
		return EXIT_FAILURE;
	}	
	if(old_world!=NULL&&old_world!=world&&ToInternal==false)
	{
		fprintf(stderr,"ERROR: ChangeUnits: Can't return to External Parameters using diffrent ContWorld\n");
		return EXIT_FAILURE;
	}
	//check if it is already right units(is made in pre)
	if(InternalUnit==ToInternal&&old_world!=NULL)return EXIT_SUCCESS;
	
	//Do GGO stuff
	GenericGeometricalObject::ChangeUnits(world,ToInternal);

	//convert units
	int i,j;
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	
	int irc[3];
	float frc[3];
	for(j=0;j<3;j++)
	{
		irc[j]=world->GridSizeGlobal[j]/2;
		frc[j]=irc[j];
	}
	if(ToInternal)
	{
		XY[0]=XY[0]*world->GridScale+frc[0];
		XY[1]=XY[1]*world->GridScale+frc[1];
		Z[0]=Z[0]*world->GridScale+frc[2];
		Z[1]=Z[1]*world->GridScale+frc[2];
		R[0]*=world->GridScale;
		R[1]*=world->GridScale;
	}
	else
	{
		XY[0]=(XY[0]-frc[0])/world->GridScale;
		XY[1]=(XY[1]-frc[1])/world->GridScale;
		Z[0]=(Z[0]-frc[2])/world->GridScale;
		Z[1]=(Z[1]-frc[2])/world->GridScale;
		R[0]/=world->GridScale;
		R[1]/=world->GridScale;
	}
	return EXIT_SUCCESS;
}
int GOTube::BuildDistMapsFromAtomsCenter(BuildWorldNI *Builder,ContWorld* world,float *Dist,float Rmax,float DoNotConsiderAtomsSmallerThen)
{
	DbgPrint2("GOTube::BuildDistMapsFromAtomsCenter\n");
	
	float Ratom=1.6*world->GridScale;
	float X=XY[0];
	float Y=XY[1];
	float Z1=Z[0]+Ratom;
	float Z2=Z[1]-Ratom;
	float R1=R[0]+Ratom;
	float R2=R[1]-Ratom;
	
	int i,j,k,gridpoint;
	int ix,iy,iz;
	int GSX=world->GridSizeGlobal[0];
	int GSY=world->GridSizeGlobal[1];
	int GSZ=world->GridSizeGlobal[2];
	int GSXY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
	int GSXYZ=world->GridSizeGlobal[0]*world->GridSizeGlobal[1]*world->GridSizeGlobal[2];
	float RSQ,Rtmp,fz;
	float RSQ1=R1*R1;
	float RSQ2=R2*R2;
	
	int start[3],end[3];
  //GODonut
	/*GODonut *goDonut=new GODonut(X,Y,(float)Z1,R1,0.0,1.0);
	goDonut->Epsilon=Epsilon;
	goDonut->NIonsTypes=NIonsTypes;
	goDonut->IonsD=new int[NIonsTypes];
	goDonut->C=new float[NIonsTypes];
	for(i=0;i<NIonsTypes;i++)
	{
		goDonut->IonsD[i]=IonsD[i];
		goDonut->C[i]=C[i];
	}
	goDonut->X=X;
	goDonut->Y=Y;
	goDonut->Z=Z1;
	goDonut->R=R1;
	goDonut->r=0.0f;
	goDonut->BuildDistMapsFromAtomsCenter(Builder, world,Dist,Rmax, DoNotConsiderAtomsSmallerThen);
	goDonut->Z=Z2;
	goDonut->BuildDistMapsFromAtomsCenter(Builder, world,Dist,Rmax, DoNotConsiderAtomsSmallerThen); 
	goDonut->R=R2;
	goDonut->Z=Z1;
	goDonut->BuildDistMapsFromAtomsCenter(Builder, world,Dist,Rmax, DoNotConsiderAtomsSmallerThen); 
	goDonut->Z=Z2;
	goDonut->BuildDistMapsFromAtomsCenter(Builder, world,Dist,Rmax, DoNotConsiderAtomsSmallerThen);
  */
  
	start[0]=(int)(X-R2-Rmax-1);
	end[0]=(int)(X+R2+Rmax+1);
	start[1]=(int)(Y-R2-Rmax-1);
	end[1]=(int)(Y+R2+Rmax+1);
	start[2]=(int)(Z1-Rmax-1);
	end[2]=(int)(Z2+Rmax+1);
	if(start[0]<0)start[0]=0;
	if(end[0]>GSX-1)end[0]=GSX-1;
	if(start[1]<0)start[1]=0;
	if(end[1]>GSY-1)end[1]=GSY-1;
	if(start[2]<0)start[2]=0;
	if(end[2]>GSZ-1)end[2]=GSZ-1;
  
 
  //@to do if hole small
  
	//cout<<"start ["<<start[0]<<","<<start[1]<<","<<start[2]<<"]\n";
	//cout<<"end ["<<end[0]<<","<<end[1]<<","<<end[2]<<"]\n";
	for(ix=start[0];ix<=end[0];ix++)
		for(iy=start[1];iy<=end[1];iy++)
			for(iz=start[2];iz<=end[2];iz++)
	{
		gridpoint=ix+iy*GSX+iz*GSXY;
		RSQ=(ix-X)*(ix-X)+(iy-Y)*(iy-Y);
		Rtmp=sqrt(RSQ);
		if(Z1<=iz&&iz<=Z2)
		{
			if(RSQ1<=RSQ&&RSQ<=RSQ2)
			{
				Dist[gridpoint]=0.0;
			}
			/*else if (RSQ<RSQ1&&R1-Rtmp<Dist[gridpoint])
			{
				Dist[gridpoint]=R1-Rtmp;
			}
			else if (Rtmp-R2<Dist[gridpoint])//i.e RSQ>RSQ2
			{
				Dist[gridpoint]=Rtmp-R2;
			}*/
		}
		else if (iz<Z1)
		{
			if(RSQ1<=RSQ&&RSQ<=RSQ2&&Z1-iz<Dist[gridpoint])
			{
				Dist[gridpoint]=Z1-iz;
			}
			/*else if (RSQ<RSQ1)
			{
				Dist[gridpoint]=R1-sqrt(RSQ);
			}
			else //i.e RSQ>RSQ2
			{
				Dist[gridpoint]=sqrt(RSQ)-R2;
			}*/
		}
		else //i.e. iz>Z2
		{
			if(RSQ1<=RSQ&&RSQ<=RSQ2&&iz-Z2<Dist[gridpoint])
			{
				Dist[gridpoint]=iz-Z2;
			}
			/*else if (RSQ<RSQ1)
			{
				Dist[gridpoint]=R1-sqrt(RSQ);
			}
			else //i.e RSQ>RSQ2
			{
				Dist[gridpoint]=sqrt(RSQ)-R2;
			}*/
		}
	}
	//delete goDonut;
	return EXIT_SUCCESS;
}
int GOTube::BuildLJRepultionMap(BuildWorldNI *Builder,ContWorld* world,int Ion,float *V,float LimitVlj)
{
	DbgPrint2("GOTube::BuildDistMapsFromAtomsCenter\n");
	DbgPrint2("impliment me\n");
	/*
	float HSigma=IonsHalfSigma+HalfSigma;
	float LimitVljSigma=LimitVlj*HSigma;
	float FEpsilon=sqrt(IonsFourEpsilon*FourEpsilon);
	float HSigmaSQ=HSigma*HSigma;
	float HSigma6=HSigmaSQ*HSigmaSQ*HSigmaSQ;
	float HSigma12=HSigma6*HSigma6;
	float RtmpSQ,R6,R12;
	
	float Ratom=HalfSigma;
	//float Ratom=0.0;
	float X=XY[0];
	float Y=XY[1];
	float Z1=Z[0]+Ratom;
	float Z2=Z[1]-Ratom;
	float R1=R[0]+Ratom;
	float R2=R[1]-Ratom;
	
	int i,j,k,gridpoint;
	int ix,iy,iz;
	int GSX=world->GridSizeGlobal[0];
	int GSY=world->GridSizeGlobal[1];
	int GSZ=world->GridSizeGlobal[2];
	int GSXY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
	int GSXYZ=world->GridSizeGlobal[0]*world->GridSizeGlobal[1]*world->GridSizeGlobal[2];
	float RSQ,Rtmp,fz;
	float RSQ1=R1*R1;
	float RSQ2=R2*R2;
	
	int start[3],end[3];

  
	start[0]=(int)(X-R2-LimitVljSigma-1);
	end[0]=(int)(X+R2+LimitVljSigma+1);
	start[1]=(int)(Y-R2-LimitVljSigma-1);
	end[1]=(int)(Y+R2+LimitVljSigma+1);
	start[2]=(int)(Z1-LimitVljSigma-1);
	end[2]=(int)(Z2+LimitVljSigma+1);
	if(start[0]<0)start[0]=0;
	if(end[0]>GSX-1)end[0]=GSX-1;
	if(start[1]<0)start[1]=0;
	if(end[1]>GSY-1)end[1]=GSY-1;
	if(start[2]<0)start[2]=0;
	if(end[2]>GSZ-1)end[2]=GSZ-1;
  
 
  //@to do if hole small
  
	cout<<"start ["<<start[0]<<","<<start[1]<<","<<start[2]<<"]\n";
	cout<<"end ["<<end[0]<<","<<end[1]<<","<<end[2]<<"]\n";
	for(ix=start[0];ix<=end[0];ix++)
		for(iy=start[1];iy<=end[1];iy++)
			for(iz=start[2];iz<=end[2];iz++)
	{
		gridpoint=ix+iy*GSX+iz*GSXY;
		RSQ=(ix-X)*(ix-X)+(iy-Y)*(iy-Y);
		Rtmp=sqrt(RSQ);
		if(Z1<=iz&&iz<=Z2)
		{
			if(RSQ1<=RSQ&&RSQ<=RSQ2)
			{
				V[gridpoint]+=1.0e6;
			}
		}
		else if (iz<Z1)
		{
			if(RSQ1<=RSQ&&RSQ<=RSQ2)
			{
				Rtmp=Z1-iz;
				if(Rtmp<LimitVljSigma)
				{
					RtmpSQ=Rtmp*Rtmp;
					R6=RtmpSQ*RtmpSQ*RtmpSQ;
					R12=R6*R6;
					V[gridpoint]+=FEpsilon*((HSigma12/R12)-(HSigma6/R6));
				}
			}
		}
		else //i.e. iz>Z2
		{
			if(RSQ1<=RSQ&&RSQ<=RSQ2)
			{
				Rtmp=iz-Z2;
				if(Rtmp<LimitVljSigma)
				{
					RtmpSQ=Rtmp*Rtmp;
					R6=RtmpSQ*RtmpSQ*RtmpSQ;
					R12=R6*R6;
					V[gridpoint]+=FEpsilon*((HSigma12/R12)-(HSigma6/R6));

				}
			}
		}
#		if defined(_MSC_VER) || defined(__DECCXX) || (__GNUC__ < 3)
		if(V[gridpoint]>1.0e8||RtmpSQ<1.0)
#		else
		if(V[gridpoint]>1.0e8||RtmpSQ<1.0||V[gridpoint]==NAN||V[gridpoint]==INFINITY)
#		endif
			V[gridpoint]=1.0e8;
	}
	//delete goDonut;*/
	return EXIT_SUCCESS;
}
int GOTube::BuildPreMaps(BuildWorldNI *Builder,ContWorld* world,int *Field,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth,float *Surf,int ParmMask)
{
	DbgPrint2("GOTube::BuildPreMaps(iValue=%d iBulkValue=%d)\n"
			,iValue,iBulkValue);
	float X=XY[0];
	float Y=XY[1];
	float Z1=Z[0];
	float Z2=Z[1];
	float R1=R[0];
	float R2=R[1];
	X+=Displ[0];
	Y+=Displ[1];
	Z1+=Displ[2];
	Z2+=Displ[2];
	
	int i,j,k,gridpoint;
	int ix,iy,iz;
	int GSX=world->GridSizeGlobal[0];
	int GSY=world->GridSizeGlobal[1];
	int GSZ=world->GridSizeGlobal[2];
	int GSXY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
	int GSXYZ=world->GridSizeGlobal[0]*world->GridSizeGlobal[1]*world->GridSizeGlobal[2];
	float RSQ,Rtmp,fz;
	float Ri1=R1-Rion;
	float Ri2=R2+Rion;
	float RSQ1=Ri1*Ri1;
	float RSQ2=Ri2*Ri2;
	float Rsm1=Ri1-Rsmoth;
	float Rsm2=Ri2+Rsmoth;
	float RsmSQ1=Rsm1*Rsm1;
	float RsmSQ2=Rsm2*Rsm2;
	float RsmSQ1s=(Rsm1-0.707106781f)*(Rsm1-0.707106781f);
	float RsmSQ1b=(Rsm1+0.707106781f)*(Rsm1+0.707106781f);
	float RsmSQ2s=(Rsm2-0.707106781f)*(Rsm2-0.707106781f);
	float RsmSQ2b=(Rsm2+0.707106781f)*(Rsm2+0.707106781f);
	float RoSQ1=R1*R1;
	float RoSQ2=R2*R2;
  
	float fz1=Z1-Rion;
	float fz2=Z2+Rion;
	float fz1up=fz1-Rsmoth+0.5;
	float fz1dn=fz1-Rsmoth-0.5;
	float fz2up=fz2+Rsmoth+0.5;
	float fz2dn=fz2+Rsmoth-0.5;
	float fzsm1=fz1-Rsmoth;
	float fzsm2=fz2+Rsmoth;
	int z1=(int)(fz1+0.5);
	int z2=(int)(fz2+0.5);
	int zsm1=(int)(fzsm1+0.5);
	int zsm2=(int)(fzsm2+0.5);
	int start[3],end[3];
	float DNull[3]={0.0,0.0,0.0},ftmp;
	float *vtmp;
  
	DbgPrint0("GOTube:: Rion=%f Rsmoth=%f\n",Rion,Rsmoth);
	DbgPrint0("GOTube:: R1=%f R2=%f\n",R1,R2);
	DbgPrint0("GOTube:: Ri1=%f Ri2=%f\n",Ri1,Ri2);
	DbgPrint0("GOTube:: Rsm1=%f Rsm2=%f\n",Rsm1,Rsm2);
  
	DbgPrint0("GOTube:: fz1=%f fz2=%f\n",fz1,fz2);
	DbgPrint0("GOTube:: z1=%d z2=%d\n",z1,z2);
	DbgPrint0("GOTube:: fzsm1=%f fzsm2=%f\n",fzsm1,fzsm2);
	DbgPrint0("GOTube:: zsm1=%d zsm2=%d\n",zsm1,zsm2);
  //GODonut
	GODonut *goDonut=new GODonut(X,Y,(float)z1,R1-Rion,0.0,1.0);
	goDonut->Epsilon=Epsilon;
	goDonut->NIonsTypes=NIonsTypes;
	goDonut->IonsD=new int[NIonsTypes];
	goDonut->C=new float[NIonsTypes];
	for(i=0;i<NIonsTypes;i++)
	{
		goDonut->IonsD[i]=IonsD[i];
		goDonut->C[i]=C[i];
	}
	goDonut->X=X;
	goDonut->Y=Y;
	goDonut->Z=Z1;
	goDonut->R=R1;
	goDonut->r=0.0f;
	goDonut->BuildPreMaps(Builder, world,Field, iValue, iBulkValue,Displ,Rion, Rsmoth,Surf,0);
	goDonut->Z=Z2;
	goDonut->BuildPreMaps(Builder, world,Field, iValue, iBulkValue,Displ,Rion, Rsmoth,Surf,0);
	goDonut->R=R2;
	goDonut->Z=Z1;
	goDonut->BuildPreMaps(Builder, world,Field, iValue, iBulkValue,Displ,Rion, Rsmoth,Surf,0);
	goDonut->Z=Z2;
	goDonut->BuildPreMaps(Builder, world,Field, iValue, iBulkValue,Displ,Rion, Rsmoth,Surf,0);
  
  
	start[0]=(int)(X-Rsm2-1);
	end[0]=(int)(X+Rsm2+1);
	start[1]=(int)(Y-Rsm2-1);
	end[1]=(int)(Y+Rsm2+1);
	start[2]=zsm1-1;
	end[2]=zsm2+1;
	if(start[0]<0)start[0]=0;
	if(end[0]>GSX-1)end[0]=GSX-1;
	if(start[1]<0)start[1]=0;
	if(end[1]>GSY-1)end[1]=GSY-1;
	if(start[2]<0)start[2]=0;
	if(end[2]>GSZ-1)end[2]=GSZ-1;
  
 
  //@to do if hole small
  
	//cout<<"start ["<<start[0]<<","<<start[1]<<","<<start[2]<<"]\n";
	//cout<<"end ["<<end[0]<<","<<end[1]<<","<<end[2]<<"]\n";
	for(ix=start[0];ix<=end[0];ix++)
		for(iy=start[1];iy<=end[1];iy++)
			for(iz=start[2];iz<=end[2];iz++)
	{
		gridpoint=ix+iy*GSX+iz*GSXY;
		RSQ=(ix-X)*(ix-X)+(iy-Y)*(iy-Y);
		fz=iz;
		if(RsmSQ1s<=RSQ&&RSQ<=RsmSQ2b&&fz1dn<=fz&&fz<=fz2up)
		{
			vtmp=Surf+4*gridpoint;
			if(z1<=iz&&iz<=z2)
			{
        //DbgPrint0("GOTube:: fz=%f\n",fz);
				if(RSQ1<=RSQ&&RSQ<=RSQ2)
				{
					if(!((iz>Z2&&RSQ<RoSQ1)||(iz<Z1&&RSQ<RoSQ1)||(iz>Z2&&RSQ>RoSQ2)||(iz<Z1&&RSQ>RoSQ2)))
						if(Field[gridpoint]==iBulkValue||Field[gridpoint]<0)
					{
						Field[gridpoint]=iValue;
						vtmp[0]=-BIGDISTANSE;
						vtmp[3]=0.0;
					}
				}
				else 
				{
					if(vtmp[0]>-100.0f)//e.i. is possible intersection
					{
						if(RsmSQ1s<=RSQ&&RSQ<=RsmSQ1b)//e.i. is intersection
						{
							ftmp=sqrt(RSQ);
							Rtmp=Ri1-ftmp;
							if(Rtmp<vtmp[3])
							{
								Field[gridpoint]=-iValue;
								ftmp=Rsm1/ftmp;
								vtmp[0]=X+ftmp*(ix-X);
								vtmp[1]=Y+ftmp*(iy-Y);
								vtmp[2]=iz;
								vtmp[3]=Rtmp;
							}
						}
						else if(RsmSQ2s<=RSQ&&RSQ<=RsmSQ2b)//e.i. is intersection
						{
							ftmp=sqrt(RSQ);
							Rtmp=ftmp-Ri2;
							if(Rtmp<vtmp[3])
							{
								Field[gridpoint]=-iValue;
								ftmp=Rsm2/ftmp;
								vtmp[0]=X+ftmp*(ix-X);
								vtmp[1]=Y+ftmp*(iy-Y);
								vtmp[2]=iz;
								vtmp[3]=Rtmp;
							}
						}
						else//e.i. is no intersection
						{
              //filling Field
							if(RsmSQ1<=RSQ&&RSQ<RSQ1)
							{
								if(Field[gridpoint]==iBulkValue)Field[gridpoint]=-iValue;
								else if(Field[gridpoint]<0)
								{
									Rtmp=Ri1-sqrt(RSQ);
									if(Rtmp<vtmp[3])
									{
										vtmp[3]=Rtmp;
										Field[gridpoint]=-iValue;
									}
								}
							}
							else if(RSQ2<RSQ&&RSQ<=RsmSQ2)
							{
								if(Field[gridpoint]==iBulkValue)Field[gridpoint]=-iValue;
								else if(Field[gridpoint]<0)
								{
									Rtmp=sqrt(RSQ)-Ri2;
									if(Rtmp<vtmp[3])
									{
										vtmp[3]=Rtmp;
										Field[gridpoint]=-iValue;
									}
								}
							}
							vtmp[0]=-BIGDISTANSE;
						}
					}
					else //no intersection still
					{
            //filling Making Boarder
						if(Field[gridpoint]==iBulkValue)
						{
							if(RsmSQ1s<=RSQ&&RSQ<=RsmSQ1b)
							{
								ftmp=sqrt(RSQ);
								Rtmp=Ri1-ftmp;
                
								ftmp=Rsm1/ftmp;
								vtmp[0]=X+ftmp*(ix-X);
								vtmp[1]=Y+ftmp*(iy-Y);
								vtmp[2]=iz;
								vtmp[3]=Rtmp;
							}
							else if(RSQ<=RsmSQ2b&&RSQ>=RsmSQ2s)
							{
								ftmp=sqrt(RSQ);
								Rtmp=ftmp-Ri2;
								ftmp=Rsm2/ftmp;
								vtmp[0]=X+ftmp*(ix-X);
								vtmp[1]=Y+ftmp*(iy-Y);
								vtmp[2]=iz;
								vtmp[3]=Rtmp;
							}
						}
            //filling Field
						if(RsmSQ1<=RSQ&&RSQ<RSQ1)
						{
							if(Field[gridpoint]==iBulkValue)Field[gridpoint]=-iValue;
							else if(Field[gridpoint]<0)
							{
								Rtmp=Ri1-sqrt(RSQ);
								if(Rtmp<vtmp[3])
								{
									vtmp[3]=Rtmp;
									Field[gridpoint]=-iValue;
								}
							}
						}
						else if(RSQ2<RSQ&&RSQ<=RsmSQ2)
						{
							if(Field[gridpoint]==iBulkValue)Field[gridpoint]=-iValue;
							else if(Field[gridpoint]<0)
							{
								Rtmp=sqrt(RSQ)-Ri2;
								if(Rtmp<vtmp[3])
								{
									vtmp[3]=Rtmp;
									Field[gridpoint]=-iValue;
								}
							}
						}
					}
				}
			}//end of if(iz>=z1&&iz<=z2)
			else
			{
        //top and bottom
				if(RSQ>=RSQ1&&RSQ<=RSQ2)
				{
					if(vtmp[0]>-100.0f)//e.i. is possible intersection
					{
						if(fz<=fz1up&&fz>=fz1dn)//e.i. is intersection
						{
							Rtmp=fz1-fz;
							if(Rtmp<vtmp[3])
							{
								Field[gridpoint]=-iValue;
								vtmp[0]=ix;
								vtmp[1]=iy;
								vtmp[2]=fz1-Rsmoth;
								vtmp[3]=Rtmp;
							}
						}
						else if(fz<=fz2up&&fz>=fz2dn)//e.i. is intersection
						{
							Rtmp=fz-fz2;
							if(Rtmp<vtmp[3])
							{
								Field[gridpoint]=-iValue;
								vtmp[0]=ix;
								vtmp[1]=iy;
								vtmp[2]=fz2+Rsmoth;
								vtmp[3]=Rtmp;
							}
						}
						else//e.i. is no intersection
						{
            //filling Field
							if(fz<fz1&&fz>=fzsm1)
							{
								Rtmp=fz1-fz;
								if(Field[gridpoint]==iBulkValue)
								{
									Field[gridpoint]=-iValue;
									vtmp[3]=Rtmp;
								}
								else if(Rtmp<vtmp[3])
								{
									Field[gridpoint]=-iValue;
									vtmp[3]=Rtmp;
								}
							}
							else if(fz>fz2&&fz<=fzsm2)
							{
								Rtmp=fz-fz2;
								if(Field[gridpoint]==iBulkValue)
								{
									Field[gridpoint]=-iValue;
									vtmp[3]=Rtmp;
								}
								else if(Rtmp<vtmp[3])
								{
									Field[gridpoint]=-iValue;
									vtmp[3]=Rtmp;
								}
							}
							vtmp[0]=-BIGDISTANSE;
						}
					}
					else
					{
          //filling Making Boarder
						if(Field[gridpoint]==iBulkValue)
						{
							if(fz<=fz2up&&fz>=fz2dn)
							{
								Rtmp=fz-fz2;
								vtmp[0]=ix;
								vtmp[1]=iy;
              //vtmp[2]=fz2+Rsmoth;
								vtmp[2]=iz;
								vtmp[3]=Rtmp;
							}
							else if(fz<=fz1up&&fz>=fz1dn)
							{
								Rtmp=fz1-fz;
								vtmp[0]=ix;
								vtmp[1]=iy;
              //vtmp[2]=iz;
								vtmp[2]=fz1-Rsmoth;
								vtmp[3]=Rtmp;
							}
						}
          //filling Field
						if(fz<fz1&&fz>=fzsm1)
						{
							Rtmp=fz1-fz;
							if(Field[gridpoint]==iBulkValue)
							{
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
							else if(Rtmp<vtmp[3]&&Field[gridpoint]<0)
							{
              
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
						}
						else if(fz>fz2&&fz<=fzsm2)
						{
							Rtmp=fz-fz2;
							if(Field[gridpoint]==iBulkValue)
							{
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
							else if(Rtmp<vtmp[3]&&Field[gridpoint]<0)
							{
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
						}
					}
				}
			}
		}//end of if(RsmSQ1s<=RSQ&&RSQ<=RsmSQ2b&&fz1dn<=fz&&fz<=fz2up)
	}
	X-=Displ[0];
	Y-=Displ[1];
	Z1-=Displ[2];
	Z2-=Displ[2];
	delete goDonut;
	return EXIT_SUCCESS;
}
int GOTube::BuildPreMapsMemLim(PartialWorldForSAS *pw,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth)
{
	DbgPrint2("GOTube::BuildPreMaps(iValue=%d iBulkValue=%d)\n"
			,iValue,iBulkValue);
	float X=XY[0];
	float Y=XY[1];
	float Z1=Z[0];
	float Z2=Z[1];
	float R1=R[0];
	float R2=R[1];
	X+=Displ[0];
	Y+=Displ[1];
	Z1+=Displ[2];
	Z2+=Displ[2];
	
	int i,j,k,gridpoint;
	int ix,iy,iz;
	float *Surf=pw->Surf;
	
	int GSX=pw->GS_X;
	int GSY=pw->GS_Y;
	int GSZ=pw->locGS_Z;
	int GSXY=pw->GS_X*pw->GS_Y;
	int GSXYZ=pw->GS_X*pw->GS_Y*pw->locGS_Z;
	
	int *Field=pw->Field+GSXY*pw->locR0_Z;
	
	float RSQ,Rtmp,fz;
	float Ri1=R1-Rion;
	float Ri2=R2+Rion;
	float RSQ1=Ri1*Ri1;
	float RSQ2=Ri2*Ri2;
	float Rsm1=Ri1-Rsmoth;
	float Rsm2=Ri2+Rsmoth;
	float RsmSQ1=Rsm1*Rsm1;
	float RsmSQ2=Rsm2*Rsm2;
	float RsmSQ1s=(Rsm1-0.707106781f)*(Rsm1-0.707106781f);
	float RsmSQ1b=(Rsm1+0.707106781f)*(Rsm1+0.707106781f);
	float RsmSQ2s=(Rsm2-0.707106781f)*(Rsm2-0.707106781f);
	float RsmSQ2b=(Rsm2+0.707106781f)*(Rsm2+0.707106781f);
	float RoSQ1=R1*R1;
	float RoSQ2=R2*R2;
  
	float fz1=Z1-Rion;
	float fz2=Z2+Rion;
	float fz1up=fz1-Rsmoth+0.5;
	float fz1dn=fz1-Rsmoth-0.5;
	float fz2up=fz2+Rsmoth+0.5;
	float fz2dn=fz2+Rsmoth-0.5;
	float fzsm1=fz1-Rsmoth;
	float fzsm2=fz2+Rsmoth;
	int z1=(int)(fz1+0.5);
	int z2=(int)(fz2+0.5);
	int zsm1=(int)(fzsm1+0.5);
	int zsm2=(int)(fzsm2+0.5);
	int start[3],end[3];
	float DNull[3]={0.0,0.0,0.0},ftmp;
	float *vtmp;
  
	cout<<"Rsmoth "<<Rsmoth<<"\n";
	cout<<"zsm "<<zsm1<<" "<<zsm2<<"\n";
	cout<<"z "<<z1<<" "<<z2<<"\n";
	cout<<"Rsm "<<RsmSQ1<<" "<<RsmSQ2<<"\n";
	cout<<"RSQ "<<RSQ1<<" "<<RSQ2<<"\n";
  
	DbgPrint0("GOTube:: Rion=%f Rsmoth=%f\n",Rion,Rsmoth);
	DbgPrint0("GOTube:: R1=%f R2=%f\n",R1,R2);
	DbgPrint0("GOTube:: Ri1=%f Ri2=%f\n",Ri1,Ri2);
	DbgPrint0("GOTube:: Rsm1=%f Rsm2=%f\n",Rsm1,Rsm2);
  
	DbgPrint0("GOTube:: fz1=%f fz2=%f\n",fz1,fz2);
	DbgPrint0("GOTube:: z1=%d z2=%d\n",z1,z2);
	DbgPrint0("GOTube:: fzsm1=%f fzsm2=%f\n",fzsm1,fzsm2);
	DbgPrint0("GOTube:: zsm1=%d zsm2=%d\n",zsm1,zsm2);
  //GODonut
	GODonut *goDonut=new GODonut(X,Y,(float)z1,R1-Rion,0.0,1.0);
	goDonut->Epsilon=Epsilon;
	goDonut->NIonsTypes=NIonsTypes;
	goDonut->IonsD=new int[NIonsTypes];
	goDonut->C=new float[NIonsTypes];
	for(i=0;i<NIonsTypes;i++)
	{
		goDonut->IonsD[i]=IonsD[i];
		goDonut->C[i]=C[i];
	}
	goDonut->X=X;
	goDonut->Y=Y;
	goDonut->Z=Z1;
	goDonut->R=R1;
	goDonut->r=0.0f;
	goDonut->BuildPreMapsMemLim(pw,iValue,iBulkValue,Displ, Rion, Rsmoth);
	goDonut->Z=Z2;
	goDonut->BuildPreMapsMemLim(pw,iValue,iBulkValue,Displ, Rion, Rsmoth);
	goDonut->R=R2;
	goDonut->Z=Z1;
	goDonut->BuildPreMapsMemLim(pw,iValue,iBulkValue,Displ, Rion, Rsmoth);
	goDonut->Z=Z2;
	goDonut->BuildPreMapsMemLim(pw,iValue,iBulkValue,Displ, Rion, Rsmoth);
  
  
	start[0]=(int)(X-Rsm2-1);
	end[0]=(int)(X+Rsm2+1);
	start[1]=(int)(Y-Rsm2-1);
	end[1]=(int)(Y+Rsm2+1);
	start[2]=zsm1-1;
	end[2]=zsm2+1;
	if(start[0]<0)start[0]=0;
	if(end[0]>GSX-1)end[0]=GSX-1;
	if(start[1]<0)start[1]=0;
	if(end[1]>GSY-1)end[1]=GSY-1;
	if(start[2]<0)start[2]=0;
	if(end[2]>GSZ-1)end[2]=GSZ-1;
  
 
  //@to do if hole small
  
	//cout<<"start ["<<start[0]<<","<<start[1]<<","<<start[2]<<"]\n";
	//cout<<"end ["<<end[0]<<","<<end[1]<<","<<end[2]<<"]\n";
	for(ix=start[0];ix<=end[0];ix++)
		for(iy=start[1];iy<=end[1];iy++)
			for(iz=start[2];iz<=end[2];iz++)
	{
		gridpoint=ix+iy*GSX+iz*GSXY;
		RSQ=(ix-X)*(ix-X)+(iy-Y)*(iy-Y);
		fz=iz;
		if(RsmSQ1s<=RSQ&&RSQ<=RsmSQ2b&&fz1dn<=fz&&fz<=fz2up)
		{
			vtmp=Surf+4*gridpoint;
			if(z1<=iz&&iz<=z2)
			{
        //DbgPrint0("GOTube:: fz=%f\n",fz);
				if(RSQ1<=RSQ&&RSQ<=RSQ2)
				{
					if(!((iz>Z2&&RSQ<RoSQ1)||(iz<Z1&&RSQ<RoSQ1)||(iz>Z2&&RSQ>RoSQ2)||(iz<Z1&&RSQ>RoSQ2)))
						if(Field[gridpoint]==iBulkValue||Field[gridpoint]<0)
					{
						Field[gridpoint]=iValue;
						vtmp[0]=-BIGDISTANSE;
						vtmp[3]=0.0;
					}
				}
				else 
				{
					if(vtmp[0]>-100.0f)//e.i. is possible intersection
					{
						if(RsmSQ1s<=RSQ&&RSQ<=RsmSQ1b)//e.i. is intersection
						{
							ftmp=sqrt(RSQ);
							Rtmp=Ri1-ftmp;
							if(Rtmp<vtmp[3])
							{
								Field[gridpoint]=-iValue;
								ftmp=Rsm1/ftmp;
								vtmp[0]=X+ftmp*(ix-X);
								vtmp[1]=Y+ftmp*(iy-Y);
								vtmp[2]=iz;
								vtmp[3]=Rtmp;
							}
						}
						else if(RsmSQ2s<=RSQ&&RSQ<=RsmSQ2b)//e.i. is intersection
						{
							ftmp=sqrt(RSQ);
							Rtmp=ftmp-Ri2;
							if(Rtmp<vtmp[3])
							{
								Field[gridpoint]=-iValue;
								ftmp=Rsm2/ftmp;
								vtmp[0]=X+ftmp*(ix-X);
								vtmp[1]=Y+ftmp*(iy-Y);
								vtmp[2]=iz;
								vtmp[3]=Rtmp;
							}
						}
						else//e.i. is no intersection
						{
              //filling Field
							if(RsmSQ1<=RSQ&&RSQ<RSQ1)
							{
								if(Field[gridpoint]==iBulkValue)Field[gridpoint]=-iValue;
								else if(Field[gridpoint]<0)
								{
									Rtmp=Ri1-sqrt(RSQ);
									if(Rtmp<vtmp[3])
									{
										vtmp[3]=Rtmp;
										Field[gridpoint]=-iValue;
									}
								}
							}
							else if(RSQ2<RSQ&&RSQ<=RsmSQ2)
							{
								if(Field[gridpoint]==iBulkValue)Field[gridpoint]=-iValue;
								else if(Field[gridpoint]<0)
								{
									Rtmp=sqrt(RSQ)-Ri2;
									if(Rtmp<vtmp[3])
									{
										vtmp[3]=Rtmp;
										Field[gridpoint]=-iValue;
									}
								}
							}
							vtmp[0]=-BIGDISTANSE;
						}
					}
					else //no intersection still
					{
            //filling Making Boarder
						if(Field[gridpoint]==iBulkValue)
						{
							if(RsmSQ1s<=RSQ&&RSQ<=RsmSQ1b)
							{
								ftmp=sqrt(RSQ);
								Rtmp=Ri1-ftmp;
                
								ftmp=Rsm1/ftmp;
								vtmp[0]=X+ftmp*(ix-X);
								vtmp[1]=Y+ftmp*(iy-Y);
								vtmp[2]=iz;
								vtmp[3]=Rtmp;
							}
							else if(RSQ<=RsmSQ2b&&RSQ>=RsmSQ2s)
							{
								ftmp=sqrt(RSQ);
								Rtmp=ftmp-Ri2;
								ftmp=Rsm2/ftmp;
								vtmp[0]=X+ftmp*(ix-X);
								vtmp[1]=Y+ftmp*(iy-Y);
								vtmp[2]=iz;
								vtmp[3]=Rtmp;
							}
						}
            //filling Field
						if(RsmSQ1<=RSQ&&RSQ<RSQ1)
						{
							if(Field[gridpoint]==iBulkValue)Field[gridpoint]=-iValue;
							else if(Field[gridpoint]<0)
							{
								Rtmp=Ri1-sqrt(RSQ);
								if(Rtmp<vtmp[3])
								{
									vtmp[3]=Rtmp;
									Field[gridpoint]=-iValue;
								}
							}
						}
						else if(RSQ2<RSQ&&RSQ<=RsmSQ2)
						{
							if(Field[gridpoint]==iBulkValue)Field[gridpoint]=-iValue;
							else if(Field[gridpoint]<0)
							{
								Rtmp=sqrt(RSQ)-Ri2;
								if(Rtmp<vtmp[3])
								{
									vtmp[3]=Rtmp;
									Field[gridpoint]=-iValue;
								}
							}
						}
					}
				}
			}//end of if(iz>=z1&&iz<=z2)
			else
			{
        //top and bottom
				if(RSQ>=RSQ1&&RSQ<=RSQ2)
				{
					if(vtmp[0]>-100.0f)//e.i. is possible intersection
					{
						if(fz<=fz1up&&fz>=fz1dn)//e.i. is intersection
						{
							Rtmp=fz1-fz;
							if(Rtmp<vtmp[3])
							{
								Field[gridpoint]=-iValue;
								vtmp[0]=ix;
								vtmp[1]=iy;
								vtmp[2]=fz1-Rsmoth;
								vtmp[3]=Rtmp;
							}
						}
						else if(fz<=fz2up&&fz>=fz2dn)//e.i. is intersection
						{
							Rtmp=fz-fz2;
							if(Rtmp<vtmp[3])
							{
								Field[gridpoint]=-iValue;
								vtmp[0]=ix;
								vtmp[1]=iy;
								vtmp[2]=fz2+Rsmoth;
								vtmp[3]=Rtmp;
							}
						}
						else//e.i. is no intersection
						{
            //filling Field
							if(fz<fz1&&fz>=fzsm1)
							{
								Rtmp=fz1-fz;
								if(Field[gridpoint]==iBulkValue)
								{
									Field[gridpoint]=-iValue;
									vtmp[3]=Rtmp;
								}
								else if(Rtmp<vtmp[3])
								{
									Field[gridpoint]=-iValue;
									vtmp[3]=Rtmp;
								}
							}
							else if(fz>fz2&&fz<=fzsm2)
							{
								Rtmp=fz-fz2;
								if(Field[gridpoint]==iBulkValue)
								{
									Field[gridpoint]=-iValue;
									vtmp[3]=Rtmp;
								}
								else if(Rtmp<vtmp[3])
								{
									Field[gridpoint]=-iValue;
									vtmp[3]=Rtmp;
								}
							}
							vtmp[0]=-BIGDISTANSE;
						}
					}
					else
					{
          //filling Making Boarder
						if(Field[gridpoint]==iBulkValue)
						{
							if(fz<=fz2up&&fz>=fz2dn)
							{
								Rtmp=fz-fz2;
								vtmp[0]=ix;
								vtmp[1]=iy;
              //vtmp[2]=fz2+Rsmoth;
								vtmp[2]=iz;
								vtmp[3]=Rtmp;
							}
							else if(fz<=fz1up&&fz>=fz1dn)
							{
								Rtmp=fz1-fz;
								vtmp[0]=ix;
								vtmp[1]=iy;
              //vtmp[2]=iz;
								vtmp[2]=fz1-Rsmoth;
								vtmp[3]=Rtmp;
							}
						}
          //filling Field
						if(fz<fz1&&fz>=fzsm1)
						{
							Rtmp=fz1-fz;
							if(Field[gridpoint]==iBulkValue)
							{
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
							else if(Rtmp<vtmp[3]&&Field[gridpoint]<0)
							{
              
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
						}
						else if(fz>fz2&&fz<=fzsm2)
						{
							Rtmp=fz-fz2;
							if(Field[gridpoint]==iBulkValue)
							{
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
							else if(Rtmp<vtmp[3]&&Field[gridpoint]<0)
							{
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
						}
					}
				}
			}
		}//end of if(RsmSQ1s<=RSQ&&RSQ<=RsmSQ2b&&fz1dn<=fz&&fz<=fz2up)
	}
	X-=Displ[0];
	Y-=Displ[1];
	Z1-=Displ[2];
	Z2-=Displ[2];
	delete goDonut;
	return EXIT_SUCCESS;
}
int GOTube::RotateGGO(double *n, double cosa, double sina)
{
	int i,j;
	double rtemp[3];
	pnpWarning("Can rotate GOTube only around Z\n");
	
	rtemp[0]=XY[0];
	rtemp[1]=XY[1];
	rtemp[2]=0.0;
	RotateVecDouble(rtemp,n,cosa,sina);
	XY[0]=rtemp[0];
	XY[1]=rtemp[1];
	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
GOMembraneZ::GOMembraneZ()
{
	InitZero();
}
GOMembraneZ::~GOMembraneZ()
{
	Clear();
}
int GOMembraneZ::InitZero()
{
	GenericGeometricalObject::InitZero();
	PnpsObject::SetName("MembraneZParameters");
	return EXIT_SUCCESS;
}
int GOMembraneZ::Clear()
{
	GenericGeometricalObject::Clear();
	return EXIT_SUCCESS;
}
int GOMembraneZ::Print(BuildWorldNI *buildWorld)
{
	GenericGeometricalObject::Print(buildWorld);
	pnpPrintGroup0("BW        Z Boarders:........................ [%.3f A, %.3f A]\n",Z[0],Z[1]);
	return EXIT_SUCCESS;
}
int GOMembraneZ::ChangeUnits(ContWorld* world,bool ToInternal)
{
	if(world==NULL)
	{
		fprintf(stderr,"ERROR: ChangeUnits: world==NULL\n");
		return EXIT_FAILURE;
	}	
	if(old_world!=NULL&&old_world!=world&&ToInternal==false)
	{
		fprintf(stderr,"ERROR: ChangeUnits: Can't return to External Parameters using diffrent ContWorld\n");
		return EXIT_FAILURE;
	}
	//check if it is already right units(is made in pre)
	if(InternalUnit==ToInternal&&old_world!=NULL)return EXIT_SUCCESS;
	
	//Do GGO stuff
	GenericGeometricalObject::ChangeUnits(world,ToInternal);

	//convert units
	int i,j;
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	
	int irc[3];
	float frc[3];
	for(j=0;j<3;j++)
	{
		irc[j]=world->GridSizeGlobal[j]/2;
		frc[j]=irc[j];
	}
	if(ToInternal)
	{
		Z[0]=Z[0]*world->GridScale+frc[2];
		Z[1]=Z[1]*world->GridScale+frc[2];
	}
	else
	{
		Z[0]=(Z[0]-frc[2])/world->GridScale;
		Z[1]=(Z[1]-frc[2])/world->GridScale;
	}
	return EXIT_SUCCESS;
}
int GOMembraneZ::BuildPreMaps(BuildWorldNI *Builder,ContWorld* world,int *Field,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth,float *Surf,int ParmMask)
{
	DbgPrint2("GOMembraneZ::BuildPreMaps(iValue=%d iBulkValue=%d)\n"
			,iValue,iBulkValue);
	float Z1=Z[0];
	float Z2=Z[1];
	Z1+=Displ[2];
	Z2+=Displ[2];
	
	int i,j,k,gridpoint;
	int ix,iy,iz;
	int GSX=world->GridSizeGlobal[0];
	int GSY=world->GridSizeGlobal[1];
	int GSZ=world->GridSizeGlobal[2];
	int GSXY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
	int GSXYZ=world->GridSizeGlobal[0]*world->GridSizeGlobal[1]*world->GridSizeGlobal[2];
	float Rtmp,fz;
	float fz1=Z1-Rion;
	float fz2=Z2+Rion;
	float fz1up=fz1-Rsmoth+0.5;
	float fz1dn=fz1-Rsmoth-0.5;
	float fz2up=fz2+Rsmoth+0.5;
	float fz2dn=fz2+Rsmoth-0.5;
	float fzsm1=fz1-Rsmoth;
	float fzsm2=fz2+Rsmoth;
	int z1=(int)(fz1+0.5);
	int z2=(int)(fz2+0.5);
	int zsm1=(int)(fzsm1+0.5);
	int zsm2=(int)(fzsm2+0.5);
	int start[3],end[3];
	float DNull[3]={0.0,0.0,0.0},ftmp;
	float *vtmp;

	DbgPrint0("GOMembraneZ:: fz1=%f fz2=%f\n",fz1,fz2);
	DbgPrint0("GOMembraneZ:: z1=%d z2=%d\n",z1,z2);
	DbgPrint0("GOMembraneZ:: fzsm1=%f fzsm2=%f\n",fzsm1,fzsm2);
	DbgPrint0("GOMembraneZ:: zsm1=%d zsm2=%d\n",zsm1,zsm2);
	
	start[0]=0;
	end[0]=GSX-1;
	start[1]=0;
	end[1]=GSY-1;
	start[2]=zsm1-1;
	end[2]=zsm2+1;
	if(start[2]<0)start[2]=0;
	if(end[2]>GSZ-1)end[2]=GSZ-1;
  
 
  //@to do if hole small
  
	//cout<<"start ["<<start[0]<<","<<start[1]<<","<<start[2]<<"]\n";
	//cout<<"end ["<<end[0]<<","<<end[1]<<","<<end[2]<<"]\n";
	for(ix=start[0];ix<=end[0];ix++)
		for(iy=start[1];iy<=end[1];iy++)
			for(iz=start[2];iz<=end[2];iz++)
	{
		gridpoint=ix+iy*GSX+iz*GSXY;
		
		fz=iz;
		if(fz1dn<=fz&&fz<=fz2up)
		{
			vtmp=Surf+4*gridpoint;
			if(z1<=iz&&iz<=z2)
			{
				if(Field[gridpoint]==iBulkValue||Field[gridpoint]<0)
				{
					Field[gridpoint]=iValue;
					vtmp[0]=-BIGDISTANSE;
					vtmp[3]=0.0;
				}
			}//end of if(iz>=z1&&iz<=z2)
			else
			{
				//top and bottom
				if(vtmp[0]>-100.0f)//e.i. is possible intersection
				{
					if(fz<=fz1up&&fz>=fz1dn)//e.i. is intersection
					{
						Rtmp=fz1-fz;
						if(Rtmp<vtmp[3])
						{
							Field[gridpoint]=-iValue;
							vtmp[0]=ix;
							vtmp[1]=iy;
							vtmp[2]=fz1-Rsmoth;
							vtmp[3]=Rtmp;
						}
					}
					else if(fz<=fz2up&&fz>=fz2dn)//e.i. is intersection
					{
						Rtmp=fz-fz2;
						if(Rtmp<vtmp[3])
						{
							Field[gridpoint]=-iValue;
							vtmp[0]=ix;
							vtmp[1]=iy;
							vtmp[2]=fz2+Rsmoth;
							vtmp[3]=Rtmp;
						}
					}
					else//e.i. is no intersection
					{
					//filling Field
						if(fz<fz1&&fz>=fzsm1)
						{
							Rtmp=fz1-fz;
							if(Field[gridpoint]==iBulkValue)
							{
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
							else if(Rtmp<vtmp[3])
							{
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
						}
						else if(fz>fz2&&fz<=fzsm2)
						{
							Rtmp=fz-fz2;
							if(Field[gridpoint]==iBulkValue)
							{
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
							else if(Rtmp<vtmp[3])
							{
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
						}
						vtmp[0]=-BIGDISTANSE;
					}
				}
				else
				{
					//filling Making Boarder
					if(Field[gridpoint]==iBulkValue)
					{
						if(fz<=fz2up&&fz>=fz2dn)
						{
							Rtmp=fz1-fz;
							vtmp[0]=ix;
							vtmp[1]=iy;
						//vtmp[2]=fz2+Rsmoth;
							vtmp[2]=iz;
							vtmp[3]=Rtmp;
						}
						else if(fz<=fz1up&&fz>=fz1dn)
						{
							Rtmp=iz-fz2;
							vtmp[0]=ix;
							vtmp[1]=iy;
						//vtmp[2]=iz;
							vtmp[2]=fz1-Rsmoth;
							vtmp[3]=Rtmp;
						}
					}
					//filling Field
					if(fz<fz1&&fz>=fzsm1)
					{
						Rtmp=fz1-fz;
						if(Field[gridpoint]==iBulkValue)
						{
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp;
						}
						else if(Rtmp<vtmp[3]&&Field[gridpoint]<0)
						{
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp;
						}
					}
					else if(fz>fz2&&fz<=fzsm2)
					{
						Rtmp=iz-fz2;
						if(Field[gridpoint]==iBulkValue)
						{
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp;
						}
						else if(Rtmp<vtmp[3]&&Field[gridpoint]<0)
						{
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp;
						}
					}
				}
			}
		}//end of if(RsmSQ1s<=RSQ&&RSQ<=RsmSQ2b&&fz1dn<=fz&&fz<=fz2up)
	}
	return EXIT_SUCCESS;
}
int GOMembraneZ::BuildPreMapsMemLim(PartialWorldForSAS *pw,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth)
{
	DbgPrint2("GOMembraneZ::BuildPreMaps(iValue=%d iBulkValue=%d)\n"
			,iValue,iBulkValue);
	
	float Z1=Z[0];
	float Z2=Z[1];
	Z1+=Displ[2];
	Z2+=Displ[2];
	
	int i,j,k,gridpoint;
	int ix,iy,iz;
	
	float *Surf=pw->Surf;
	
	int GSX=pw->GS_X;
	int GSY=pw->GS_Y;
	int GSZ=pw->locGS_Z;
	int GSXY=pw->GS_X*pw->GS_Y;
	int GSXYZ=pw->GS_X*pw->GS_Y*pw->locGS_Z;
	
	int *Field=pw->Field+GSXY*pw->locR0_Z;
	
	float Rtmp,fz;
	float fz1=Z1-Rion;
	float fz2=Z2+Rion;
	float fz1up=fz1-Rsmoth+0.5;
	float fz1dn=fz1-Rsmoth-0.5;
	float fz2up=fz2+Rsmoth+0.5;
	float fz2dn=fz2+Rsmoth-0.5;
	float fzsm1=fz1-Rsmoth;
	float fzsm2=fz2+Rsmoth;
	int z1=(int)(fz1+0.5);
	int z2=(int)(fz2+0.5);
	int zsm1=(int)(fzsm1+0.5);
	int zsm2=(int)(fzsm2+0.5);
	int start[3],end[3];
	float DNull[3]={0.0,0.0,0.0},ftmp;
	float *vtmp;

	DbgPrint0("GOMembraneZ:: fz1=%f fz2=%f\n",fz1,fz2);
	DbgPrint0("GOMembraneZ:: z1=%d z2=%d\n",z1,z2);
	DbgPrint0("GOMembraneZ:: fzsm1=%f fzsm2=%f\n",fzsm1,fzsm2);
	DbgPrint0("GOMembraneZ:: zsm1=%d zsm2=%d\n",zsm1,zsm2);
	
	start[0]=0;
	end[0]=GSX-1;
	start[1]=0;
	end[1]=GSY-1;
	start[2]=zsm1-1;
	end[2]=zsm2+1;
	if(start[2]<0)start[2]=0;
	if(end[2]>GSZ-1)end[2]=GSZ-1;
  
 
  //@to do if hole small
  
	//cout<<"start ["<<start[0]<<","<<start[1]<<","<<start[2]<<"]\n";
	//cout<<"end ["<<end[0]<<","<<end[1]<<","<<end[2]<<"]\n";
	for(ix=start[0];ix<=end[0];ix++)
		for(iy=start[1];iy<=end[1];iy++)
			for(iz=start[2];iz<=end[2];iz++)
	{
		gridpoint=ix+iy*GSX+iz*GSXY;
		
		fz=iz;
		if(fz1dn<=fz&&fz<=fz2up)
		{
			vtmp=Surf+4*gridpoint;
			if(z1<=iz&&iz<=z2)
			{
				if(Field[gridpoint]==iBulkValue||Field[gridpoint]<0)
				{
					Field[gridpoint]=iValue;
					vtmp[0]=-BIGDISTANSE;
					vtmp[3]=0.0;
				}
			}//end of if(iz>=z1&&iz<=z2)
			else
			{
				//top and bottom
				if(vtmp[0]>-100.0f)//e.i. is possible intersection
				{
					if(fz<=fz1up&&fz>=fz1dn)//e.i. is intersection
					{
						Rtmp=fz1-fz;
						if(Rtmp<vtmp[3])
						{
							Field[gridpoint]=-iValue;
							vtmp[0]=ix;
							vtmp[1]=iy;
							vtmp[2]=fz1-Rsmoth;
							vtmp[3]=Rtmp;
						}
					}
					else if(fz<=fz2up&&fz>=fz2dn)//e.i. is intersection
					{
						Rtmp=fz-fz2;
						if(Rtmp<vtmp[3])
						{
							Field[gridpoint]=-iValue;
							vtmp[0]=ix;
							vtmp[1]=iy;
							vtmp[2]=fz2+Rsmoth;
							vtmp[3]=Rtmp;
						}
					}
					else//e.i. is no intersection
					{
					//filling Field
						if(fz<fz1&&fz>=fzsm1)
						{
							Rtmp=fz1-fz;
							if(Field[gridpoint]==iBulkValue)
							{
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
							else if(Rtmp<vtmp[3])
							{
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
						}
						else if(fz>fz2&&fz<=fzsm2)
						{
							Rtmp=fz-fz2;
							if(Field[gridpoint]==iBulkValue)
							{
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
							else if(Rtmp<vtmp[3])
							{
								Field[gridpoint]=-iValue;
								vtmp[3]=Rtmp;
							}
						}
						vtmp[0]=-BIGDISTANSE;
					}
				}
				else
				{
					//filling Making Boarder
					if(Field[gridpoint]==iBulkValue)
					{
						if(fz<=fz2up&&fz>=fz2dn)
						{
							Rtmp=fz1-fz;
							vtmp[0]=ix;
							vtmp[1]=iy;
						//vtmp[2]=fz2+Rsmoth;
							vtmp[2]=iz;
							vtmp[3]=Rtmp;
						}
						else if(fz<=fz1up&&fz>=fz1dn)
						{
							Rtmp=iz-fz2;
							vtmp[0]=ix;
							vtmp[1]=iy;
						//vtmp[2]=iz;
							vtmp[2]=fz1-Rsmoth;
							vtmp[3]=Rtmp;
						}
					}
					//filling Field
					if(fz<fz1&&fz>=fzsm1)
					{
						Rtmp=fz1-fz;
						if(Field[gridpoint]==iBulkValue)
						{
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp;
						}
						else if(Rtmp<vtmp[3]&&Field[gridpoint]<0)
						{
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp;
						}
					}
					else if(fz>fz2&&fz<=fzsm2)
					{
						Rtmp=iz-fz2;
						if(Field[gridpoint]==iBulkValue)
						{
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp;
						}
						else if(Rtmp<vtmp[3]&&Field[gridpoint]<0)
						{
							Field[gridpoint]=-iValue;
							vtmp[3]=Rtmp;
						}
					}
				}
			}
		}//end of if(RsmSQ1s<=RSQ&&RSQ<=RsmSQ2b&&fz1dn<=fz&&fz<=fz2up)
	}
	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
GODonut::GODonut()
{
	InitZero();
}
GODonut::GODonut(float x, float y,float z,float _R,float _r,float Scale)
{
	InitZero();
	R=_R*Scale;
	r=_r*Scale;
	Z=z*Scale;
	X=x*Scale;
	Y=y*Scale;
}
GODonut::~GODonut()
{
	Clear();
}
int GODonut::InitZero()
{
	GenericGeometricalObject::InitZero();
	PnpsObject::SetName("GODonut");
	return EXIT_SUCCESS;
}
int GODonut::Clear()
{
	GenericGeometricalObject::Clear();
	return EXIT_SUCCESS;
}
int GODonut::ChangeUnits(ContWorld* world,bool ToInternal)
{
	if(world==NULL)
	{
		fprintf(stderr,"ERROR: ChangeUnits: world==NULL\n");
		return EXIT_FAILURE;
	}	
	if(old_world!=NULL&&old_world!=world&&ToInternal==false)
	{
		fprintf(stderr,"ERROR: ChangeUnits: Can't return to External Parameters using diffrent ContWorld\n");
		return EXIT_FAILURE;
	}
	//check if it is already right units(is made in pre)
	if(InternalUnit==ToInternal&&old_world!=NULL)return EXIT_SUCCESS;
	
	//Do GGO stuff
	GenericGeometricalObject::ChangeUnits(world,ToInternal);
	
	//convert units
	int i,j;
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	
	int irc[3];
	float frc[3];
	for(j=0;j<3;j++)
	{
		irc[j]=world->GridSizeGlobal[j]/2;
		frc[j]=irc[j];
	}
	if(ToInternal)
	{
		X=X*world->GridScale+frc[0];
		Y=Y*world->GridScale+frc[1];
		Z=Z*world->GridScale+frc[2];
		R*=world->GridScale;
		r*=world->GridScale;
	}
	else
	{
		X=(X-frc[0])/world->GridScale;
		Y=(Y-frc[1])/world->GridScale;
		Z=(Z-frc[2])/world->GridScale;
		R/=world->GridScale;
		r/=world->GridScale;
	}
	return EXIT_SUCCESS;
}
int GODonut::BuildPreMaps(BuildWorldNI *Builder,ContWorld* world,int *Field,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth,float *Surf,int ParmMask)
{
	int ix,iy,iz,gridpoint;
	int GSX=world->GridSizeGlobal[0];
	int GSY=world->GridSizeGlobal[1];
	int GSZ=world->GridSizeGlobal[2];
	int GSXY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
	int GSXYZ=world->GridSizeGlobal[0]*world->GridSizeGlobal[1]*world->GridSizeGlobal[2];
	int start[3],end[3];
	float RSQ,rSQ,x,y,z,Rtmp,r1SQ,rsm1SQ,ftmp;
	float RSQ1,RSQ2;
	float *vtmp;
	float Rxy[3],RxyDSQ,RxyD,Rsur[3];
	X+=Displ[0];
	Y+=Displ[1];
	Z+=Displ[2];
	
	float r1=r+Rion;
	float rsm1=r+Rion+Rsmoth;
	
	float RSQ1b=R+rsm1+0.707106781f;
	RSQ1b*=RSQ1b;
	float RSQ1s=R-rsm1-0.707106781f;
	RSQ1s*=RSQ1s;
	
	RSQ1=R-rsm1;
	RSQ1*=RSQ1;
	RSQ2=R+rsm1;
	RSQ2*=RSQ2;
	
	r1SQ=r1*r1;
	rsm1SQ=rsm1*rsm1;
	
	float rsm1SQs=(rsm1-0.707106781f)*(rsm1-0.707106781f);
	float rsm1SQb=(rsm1+0.707106781f)*(rsm1+0.707106781f);
	
	start[0]=(int)(X-(R+r+Rion+Rsmoth)+0.5)-1;
	end[0]=(int)(X+(R+r+Rion+Rsmoth)+0.5)+1;
	start[1]=(int)(Y-(R+r+Rion+Rsmoth)+0.5)-1;
	end[1]=(int)(Y+(R+r+Rion+Rsmoth)+0.5)+1;
	start[2]=(int)(Z-(r+Rion+Rsmoth)+0.5)-1;
	end[2]=(int)(Z+(r+Rion+Rsmoth)+0.5)+1;
	
	if(start[0]<0)start[0]=0;
	if(end[0]>GSX-1)end[0]=GSX-1;
	if(start[1]<0)start[1]=0;
	if(end[1]>GSY-1)end[1]=GSY-1;
	if(start[2]<0)start[2]=0;
	if(end[2]>GSZ-1)end[2]=GSZ-1;
	
	for(ix=start[0];ix<=end[0];ix++)
		for(iy=start[1];iy<=end[1];iy++)
			for(iz=start[2];iz<=end[2];iz++)
	{
		gridpoint=ix+iy*GSX+iz*GSXY;
		x=ix-X;
		y=iy-Y;
		z=iz-Z;
		RSQ=x*x+y*y+z*z;
		RxyDSQ=x*x+y*y;
		RxyD=sqrt(RxyDSQ);
		vtmp=Surf+4*gridpoint;
		if(RSQ1s<=RxyDSQ&&RxyDSQ<=RSQ1b)
		{
			rSQ=R-RxyD;
			rSQ=rSQ*rSQ+z*z;
			
			if(rSQ<=r1SQ){
				if(Field[gridpoint]==iBulkValue||Field[gridpoint]<0)
				{
					Field[gridpoint]=iValue;
					vtmp[0]=-BIGDISTANSE;
				}
			}
			else if(rSQ<=rsm1SQs){
				Rtmp=sqrt(rSQ)-r-Rion;
				if(Field[gridpoint]==iBulkValue)
				{
					Field[gridpoint]=-iValue;
					vtmp[0]=-BIGDISTANSE;
					vtmp[3]=Rtmp;
				}
				else if(Field[gridpoint]<0)
				{
					if(Rtmp<vtmp[0])
					{
						Field[gridpoint]=-iValue;
						vtmp[0]=-BIGDISTANSE;
						vtmp[3]=Rtmp;
					}
				}
			}
			else if(rSQ<=rsm1SQb)//e.i. is possible boarder
			{
				if(vtmp[0]>-100.0f)//e.i. is possible intersection
				{
					Rtmp=sqrt(rSQ)-r-Rion;
					if(Rtmp<vtmp[0])
					{
						Field[gridpoint]=-iValue;
						Rxy[0]=X+x*R/RxyD;
						Rxy[1]=Y+y*R/RxyD;
						Rxy[2]=Z;
						Rsur[0]=ix-Rxy[0];
						Rsur[1]=iy-Rxy[1];
						Rsur[2]=iz-Rxy[2];
						ftmp=sqrt(rSQ);
						Rtmp=ftmp-r-Rion;
						ftmp=rsm1/ftmp;
						vtmp[0]=Rxy[0]+ftmp*Rsur[0];
						vtmp[1]=Rxy[1]+ftmp*Rsur[1];
						vtmp[2]=Rxy[2]+ftmp*Rsur[2];
						vtmp[3]=Rtmp;
					}
				}
				else//e.i. simple boarder
				{
					if(Field[gridpoint]==iBulkValue)
					{
						Field[gridpoint]=-iValue;
						Rxy[0]=X+x*R/RxyD;
						Rxy[1]=Y+y*R/RxyD;
						Rxy[2]=Z;
						Rsur[0]=ix-Rxy[0];
						Rsur[1]=iy-Rxy[1];
						Rsur[2]=iz-Rxy[2];
						ftmp=sqrt(rSQ);
						Rtmp=ftmp-r-Rion;
						ftmp=rsm1/ftmp;
						vtmp[0]=Rxy[0]+ftmp*Rsur[0];
						vtmp[1]=Rxy[1]+ftmp*Rsur[1];
						vtmp[2]=Rxy[2]+ftmp*Rsur[2];
						vtmp[3]=Rtmp;
					}
					else if(Field[gridpoint]<0)
					{
						Rtmp=sqrt(rSQ)-r-Rion;
						if(Rtmp<vtmp[0])
						{
							vtmp[0]=-BIGDISTANSE;
							vtmp[3]=Rtmp;
						}
					}
				}
			}
		}
	}
	X-=Displ[0];
	Y-=Displ[1];
	Z-=Displ[2];
	return EXIT_SUCCESS;
}
int GODonut::BuildPreMapsMemLim(PartialWorldForSAS *pw,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth)
{
	int ix,iy,iz,gridpoint;
	float *Surf=pw->Surf;
	
	int GSX=pw->GS_X;
	int GSY=pw->GS_Y;
	int GSZ=pw->locGS_Z;
	int GSXY=pw->GS_X*pw->GS_Y;
	int GSXYZ=pw->GS_X*pw->GS_Y*pw->locGS_Z;
	
	int *Field=pw->Field+GSXY*pw->locR0_Z;
	int start[3],end[3];
	float RSQ,rSQ,x,y,z,Rtmp,r1SQ,rsm1SQ,ftmp;
	float RSQ1,RSQ2;
	float *vtmp;
	float Rxy[3],RxyDSQ,RxyD,Rsur[3];
	X+=Displ[0];
	Y+=Displ[1];
	Z+=Displ[2];
	
	float r1=r+Rion;
	float rsm1=r+Rion+Rsmoth;
	
	float RSQ1b=R+rsm1+0.707106781f;
	RSQ1b*=RSQ1b;
	float RSQ1s=R-rsm1-0.707106781f;
	RSQ1s*=RSQ1s;
	
	RSQ1=R-rsm1;
	RSQ1*=RSQ1;
	RSQ2=R+rsm1;
	RSQ2*=RSQ2;
	
	r1SQ=r1*r1;
	rsm1SQ=rsm1*rsm1;
	
	float rsm1SQs=(rsm1-0.707106781f)*(rsm1-0.707106781f);
	float rsm1SQb=(rsm1+0.707106781f)*(rsm1+0.707106781f);
	
	start[0]=(int)(X-(R+r+Rion+Rsmoth)+0.5)-1;
	end[0]=(int)(X+(R+r+Rion+Rsmoth)+0.5)+1;
	start[1]=(int)(Y-(R+r+Rion+Rsmoth)+0.5)-1;
	end[1]=(int)(Y+(R+r+Rion+Rsmoth)+0.5)+1;
	start[2]=(int)(Z-(r+Rion+Rsmoth)+0.5)-1;
	end[2]=(int)(Z+(r+Rion+Rsmoth)+0.5)+1;
	
	if(start[0]<0)start[0]=0;
	if(end[0]>GSX-1)end[0]=GSX-1;
	if(start[1]<0)start[1]=0;
	if(end[1]>GSY-1)end[1]=GSY-1;
	if(start[2]<0)start[2]=0;
	if(end[2]>GSZ-1)end[2]=GSZ-1;
	
	for(ix=start[0];ix<=end[0];ix++)
		for(iy=start[1];iy<=end[1];iy++)
			for(iz=start[2];iz<=end[2];iz++)
	{
		gridpoint=ix+iy*GSX+iz*GSXY;
		x=ix-X;
		y=iy-Y;
		z=iz-Z;
		RSQ=x*x+y*y+z*z;
		RxyDSQ=x*x+y*y;
		RxyD=sqrt(RxyDSQ);
		vtmp=Surf+4*gridpoint;
		if(RSQ1s<=RxyDSQ&&RxyDSQ<=RSQ1b)
		{
			rSQ=R-RxyD;
			rSQ=rSQ*rSQ+z*z;
			
			if(rSQ<=r1SQ){
				if(Field[gridpoint]==iBulkValue||Field[gridpoint]<0)
				{
					Field[gridpoint]=iValue;
					vtmp[0]=-BIGDISTANSE;
				}
			}
			else if(rSQ<=rsm1SQs){
				Rtmp=sqrt(rSQ)-r-Rion;
				if(Field[gridpoint]==iBulkValue)
				{
					Field[gridpoint]=-iValue;
					vtmp[0]=-BIGDISTANSE;
					vtmp[3]=Rtmp;
				}
				else if(Field[gridpoint]<0)
				{
					if(Rtmp<vtmp[0])
					{
						Field[gridpoint]=-iValue;
						vtmp[0]=-BIGDISTANSE;
						vtmp[3]=Rtmp;
					}
				}
			}
			else if(rSQ<=rsm1SQb)//e.i. is possible boarder
			{
				if(vtmp[0]>-100.0f)//e.i. is possible intersection
				{
					Rtmp=sqrt(rSQ)-r-Rion;
					if(Rtmp<vtmp[0])
					{
						Field[gridpoint]=-iValue;
						Rxy[0]=X+x*R/RxyD;
						Rxy[1]=Y+y*R/RxyD;
						Rxy[2]=Z;
						Rsur[0]=ix-Rxy[0];
						Rsur[1]=iy-Rxy[1];
						Rsur[2]=iz-Rxy[2];
						ftmp=sqrt(rSQ);
						Rtmp=ftmp-r-Rion;
						ftmp=rsm1/ftmp;
						vtmp[0]=Rxy[0]+ftmp*Rsur[0];
						vtmp[1]=Rxy[1]+ftmp*Rsur[1];
						vtmp[2]=Rxy[2]+ftmp*Rsur[2];
						vtmp[3]=Rtmp;
					}
				}
				else//e.i. simple boarder
				{
					if(Field[gridpoint]==iBulkValue)
					{
						Field[gridpoint]=-iValue;
						Rxy[0]=X+x*R/RxyD;
						Rxy[1]=Y+y*R/RxyD;
						Rxy[2]=Z;
						Rsur[0]=ix-Rxy[0];
						Rsur[1]=iy-Rxy[1];
						Rsur[2]=iz-Rxy[2];
						ftmp=sqrt(rSQ);
						Rtmp=ftmp-r-Rion;
						ftmp=rsm1/ftmp;
						vtmp[0]=Rxy[0]+ftmp*Rsur[0];
						vtmp[1]=Rxy[1]+ftmp*Rsur[1];
						vtmp[2]=Rxy[2]+ftmp*Rsur[2];
						vtmp[3]=Rtmp;
					}
					else if(Field[gridpoint]<0)
					{
						Rtmp=sqrt(rSQ)-r-Rion;
						if(Rtmp<vtmp[0])
						{
							vtmp[0]=-BIGDISTANSE;
							vtmp[3]=Rtmp;
						}
					}
				}
			}
		}
	}
	X-=Displ[0];
	Y-=Displ[1];
	Z-=Displ[2];
	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
BuildWorldNI::BuildWorldNI()
{
	InitZero();
}

BuildWorldNI::~BuildWorldNI()
{
	Clear();
}

int BuildWorldNI::InitZero()
{
	PnpsObject::SetName("BuildWorldNI");
	BoundaryStr.push_back("Zero");
	BoundaryStr.push_back("Coul");
	BoundaryCondition=ZeroBC;
	BldBCatPlane[0]=true;
	BldBCatPlane[1]=true;
	BldBCatPlane[2]=true;
	
	StdStrDiffusionModes.push_back("Plain");
	StdStrDiffusionModes.push_back("Exp");
	
	
	MakeDielectricMap=true;
	MakeDiffusionMap=true;
	MakeConcentrationMap=true;
	MakeChargeMap=true;
	MakeLJRepultion=false;
	AddElPotToPMF=false;
	
	DiffusionMode=Plain;
	Kexpdiff=0.62;
	DiffRmaxAffect=10.0;
			
	Rwat=1.4f;
	Rsmooth=0.0f;
	NIonsTypes=0;
	IonsR=NULL;
	IonsHalfSigma=NULL;
	IonsFourEpsilon=NULL;
	IonsLJA=NULL;
	IonsLJB=NULL;
	BulkParam=NULL;
	DielConst=NULL;
	DiffusionConst=NULL;
	
	RemovingCavities=true;
	RemovingCavitiesOnDielectricMap=false;
	
	LimitVlj=1.0;
	
	SaveDistMap="";
	InternalUnit=false;
	old_world=NULL;
	//GOElm.
	BuildUsingGPU=false;
	MemoryLimitOnOneArray=-1;
	return EXIT_SUCCESS;
}
int BuildWorldNI::Clear()
{	
	DeleteCArray(IonsR);
	DeleteCArray(IonsHalfSigma);
	DeleteCArray(IonsFourEpsilon);
	DeleteCArray(IonsLJA);
	DeleteCArray(IonsLJB);
	DeleteCArray(DielConst);
	DeleteCArray(DiffusionConst);
	DeleteObjByPnt(BulkParam);
	while(GOElms.size()>0)
	{
		DeleteObjByPnt(GOElms[GOElms.size()-1]);
		GOElms.pop_back();
	}
	return EXIT_SUCCESS;
}
int BuildWorldNI::Print()
{
	int i;
	pnpPrintGroup0("BW    =BuildWorldNI============================================================\n");
	//pnpPrintGroup0("BW    Number of Mobile Ions Types:........... %d\n",NIonsTypes);
	
	pnpPrintGroup0("BW    Radii of Mobile Ions:.................. [");
	for(i=0;i<NIonsTypes;i++)
	{
		pnpPrintGroup0("%.3f A",IonsR[i]);
		if(i<NIonsTypes-1) pnpPrintGroup0(", ");
		else pnpPrintGroup0("]\n");
	}
	pnpPrintGroup0("BW    Water Probe Radius for Diel. Map:...... %.3f A\n",Rwat);
	pnpPrintGroup0("BW    Smoothing Radius for Diffusion Maps:... %.3f A\n",Rsmooth);
	
	pnpPrintGroup0("BW    Boundary Condition:.................... %s\n",BoundaryStr[BoundaryCondition].c_str());
	pnpPrintGroup0("BW    Geometry Objects:\n");

	BulkParam->Print(this);
	int iElt;
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		GOElms[iElt]->Print(this);
	}
	pnpPrintGroup0("BW    =========================================================================\n");
	return EXIT_SUCCESS;
}
int BuildWorldNI::addGGO(GenericGeometricalObject* ggo)
{
	GOElms.push_back(ggo);
	return PNPS_EXIT_SUCCESS;
}
//Setter for Python Interface
int BuildWorldNI::setEpsilonValues(float e0, float e1, float e2, float e3, \
	float e4, float e5, float e6, float e7, \
	float e8, float e9, float e10, float e11, \
	float e12, float e13, float e14)
{
	DielNum=15;
	DeleteCArray(DielConst);
	DielConst=new float[DielNum];
	
	DielConst[0]=e0;
	DielConst[1]=e1;
	DielConst[2]=e2;
	DielConst[3]=e3;
	DielConst[4]=e4;
	DielConst[5]=e5;
	DielConst[6]=e6;
	DielConst[7]=e7;
	DielConst[8]=e8;
	DielConst[9]=e9;
	DielConst[10]=e10;
	DielConst[11]=e11;
	DielConst[12]=e12;
	DielConst[13]=e13;
	DielConst[14]=e14;
	
	int i;
	
	for(i=DielNum-1;i>=0;i--)
	{
		if(DielConst[i]!=0.0)break;
		DielNum--;
	}
	
	/*for(i=0;i<DielNum;i++)
	{
		pnpPrint("%f\n",DielConst[i]);
	}*/
	return PNPS_EXIT_SUCCESS;
}
int BuildWorldNI::setDiffusionValues(float D0, float D1, float D2, float D3, \
		float D4, float D5, float D6)
{
	DiffusionNum=7;
	DeleteCArray(DiffusionConst);
	DiffusionConst=new float[DiffusionNum];
	
	DiffusionConst[0]=D0;
	DiffusionConst[1]=D1;
	DiffusionConst[2]=D2;
	DiffusionConst[3]=D3;
	DiffusionConst[4]=D4;
	DiffusionConst[5]=D5;
	DiffusionConst[6]=D6;
	
	int i;
	
	for(i=DiffusionNum-1;i>=0;i--)
	{
		if(DiffusionConst[i]!=0.0)
		{
			if(DiffusionNum!=7)DiffusionNum++;
			break;
		}
		DiffusionNum--;
	}
	for(i=0;i<DiffusionNum;i++)
	{
		if(DiffusionConst[i]==0.0)
		{
			iDzero=i;
			pnpPrint("iDzero=%d %f\n",iDzero,DiffusionConst[iDzero]);
			break;
		}
	}
	for(i=0;i<DiffusionNum;i++)
	{
		pnpPrint("%f\n",DiffusionConst[i]);
	}
	return PNPS_EXIT_SUCCESS;
}
int BuildWorldNI::getEpsilonValueIndex(float epsilon)
{
	int i;
	if(DielConst==NULL)
	{
		DielNum=15;
		DielConst=new float[DielNum];
		for(i=0;i<DielNum;i++)
			DielConst[i]=0.0;
	}
	int iDiel=-1;
	for(i=0;i<DielNum;i++)
	{
		if(epsilon==DielConst[i])
			return i;
	}
	for(i=0;i<DielNum;i++)
	{
		if(DielConst[i]==0.0)
		{
			DielConst[i]=epsilon;
			return i;
		}
	}
	return -1;
}
int BuildWorldNI::getDiffusionValueIndex(float D0)
{
	int i;
	if(DiffusionConst==NULL)
	{
		DiffusionNum=7;
		DiffusionConst=new float[DiffusionNum];
		for(i=0;i<DiffusionNum;i++)
			DiffusionConst[i]=0.0;
		iDzero=0;
	}
	for(i=0;i<DiffusionNum;i++)
	{
		if(D0==DiffusionConst[i])
			return i;
	}
	for(i=1;i<DiffusionNum;i++)
	{
		if(DiffusionConst[i]==0.0)
		{
			DiffusionConst[i]=D0;
			return i;
		}
	}
	return -1;
}
int BuildWorldNI::setRwat(float m_Rwat)
{
	Rwat=m_Rwat;
	return PNPS_EXIT_SUCCESS;
}
int BuildWorldNI::setRsmooth(float m_Rsmooth)
{
	Rsmooth=m_Rsmooth;
	return PNPS_EXIT_SUCCESS;
}
int BuildWorldNI::setIonsR(float R0, float R1)
{
	NIonsTypes=2;
	DeleteCArray(IonsR);
	IonsR=new float[NIonsTypes];
	IonsR[0]=R0;
	IonsR[1]=R1;
	return PNPS_EXIT_SUCCESS;
}
int BuildWorldNI::setIonsRadii(float *m_R)
{
	int i;
	DeleteCArray(IonsR);
	IonsR=new float[NIonsTypes];
	for(i=0;i<NIonsTypes;i++)
		IonsR[i]=m_R[i];
	return PNPS_EXIT_SUCCESS;
}
int BuildWorldNI::setDiffusionMode(int DiffMode)
{
	if(!(DiffMode==Plain || DiffMode==Exp))
	{
		pnpError("Cannot recognize diffusion type\n");
		return PNPS_EXIT_FAILURE;
	}
	DiffusionMode=DiffMode;
	return PNPS_EXIT_SUCCESS;
}
int BuildWorldNI::setExpDifPar(float m_Kexpdiff, float m_Rexpdiff, float m_DiffRmaxAffect)
{
	Kexpdiff=m_Kexpdiff;
	Rexpdiff=m_Rexpdiff;
	DiffRmaxAffect=m_DiffRmaxAffect;
	return PNPS_EXIT_SUCCESS;
}
int BuildWorldNI::setMakeDielectricMap(int bVal)
{
	MakeDielectricMap=bVal;
	return PNPS_EXIT_SUCCESS;
}
int BuildWorldNI::setMakeChargeMap(int bVal)
{
	MakeChargeMap=bVal;
	return PNPS_EXIT_SUCCESS;
}
int BuildWorldNI::setMakeDiffusionMap(int bVal)
{
	MakeDiffusionMap=bVal;
	return PNPS_EXIT_SUCCESS;
}
int BuildWorldNI::setMakeConcentrationMap(int bVal)
{
	MakeConcentrationMap=bVal;
	return PNPS_EXIT_SUCCESS;
}
int BuildWorldNI::setMakeSoftRepultionMap(int bVal)
{
	MakeLJRepultion=bVal;
	//IAVMethod=3;
	return PNPS_EXIT_SUCCESS;
}
int BuildWorldNI::setBulkPar(int iEps, int iD0, int iD1,float C0,float C1)
{
	BulkParam=new GenericGeometricalObject();
	BulkParam->SetName("BulkParameters");
	BulkParam->setPar(iEps, iD0, iD1, C0, C1);
	return PNPS_EXIT_SUCCESS;
}
int BuildWorldNI::setBulkParam(int iEps, int* iD, float* C)
{
	BulkParam=new GenericGeometricalObject();
	BulkParam->SetName("BulkParameters");
	BulkParam->setParam(iEps, NIonsTypes, iD, C);
	return PNPS_EXIT_SUCCESS;
}

/*<AtomsParameters EpsilonInd="2"
		IonsDInd="3 3"
		IonsCVal="0.0 0.0"
		ChargeDist="Linear1"
		Offset="0.0 0.0 0.0"
		FileName="7ahlc_std_amber_pdb2pqr.pqr"
		FilePAN="7ahl.pan"
		MakePreRoll="true"
			/>*/
/*int BuildWorldNI::setEpsilonValues(PyObject *EpsilonValues)
{
	//is it a list
	if(!PyList_Check(EpsilonValues))
	{
		pnpError("BuildWorldNI::setEpsilonValues : EpsilonValues is not a list\n");
		return PNPS_EXIT_FAILURE;
	}
	int Nval=PyList_Size(EpsilonValues);
	//check number of epsilons
	if(Nval>NodeIndexMaxValues-1)
	{
		pnpError("BuildWorldNI::setEpsilonValues : len(EpsilonValues)=%d>%d . Can not handle it\n",Nval,NodeIndexMaxValues-1);
		return PNPS_EXIT_FAILURE;
	}
	DeleteCArray(DielConst);
	DielConst=new float[Nval+1];
	int i;
	for(i=0;i<Nval;i++)
	{
		PyObject* pyeps=PyList_GetItem(EpsilonValues, i);
		float v=0.0;
		if(PyFloat_Check(pyeps))
		{
			v=(float)PyFloat_AsDouble(pyeps);
		}
		else if(PyInt_Check(pyeps))
		{
			v=(float)PyInt_AsLong(pyeps);
		}
		else if(PyLong_Check(pyeps))
		{
			v=(float)PyLong_AsLong(pyeps);
		}
		else
		{
			pnpError("BuildWorldNI::setEpsilonValues : EpsilonValues[%d] is not a real number ( ",i);
			PyObject_Print(pyeps,stderr,Py_PRINT_RAW);
			fprintf(stderr," )\n");
			return PNPS_EXIT_FAILURE;
		}
		pnpPrint("%f\n",v);
		DielConst[i]=v;
	}
	//BuildWorldNI *bld=new BuildWorldNI();
	//PyListObject
	pnpPrint("Hello");
	return PNPS_EXIT_SUCCESS;
}*/
int BuildWorldNI::ChangeUnits(ContWorld* world,bool ToInternal)
{
	if(old_world!=NULL&&old_world!=world&&ToInternal==false)
	{
		fprintf(stderr,"ERROR: ChangeUnits: Can't return to External Parameters using diffrent ContWorld\n");
		return EXIT_FAILURE;
	}
	//check if it is already right units
	if(InternalUnit==ToInternal&&old_world!=NULL)return EXIT_SUCCESS;
	//convert units
	int i,j;
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	double GS2=world->GridScale*world->GridScale;
	double GS6=GS2*GS2*GS2;
	double GS12=GS6*GS6;
	if(ToInternal)
	{
		Rwat*=world->GridScale;
		Rsmooth*=world->GridScale;
		DiffRmaxAffect*=world->GridScale;
		for(i=0;i<DielNum;i++)
		{
			DielConst[i]/=EPKT;
		}
		for(i=0;i<NIonsTypes;i++)
		{
			IonsR[i]*=world->GridScale;
		}
		if(IonsHalfSigma!=NULL)
			for(i=0;i<NIonsTypes;i++)
		{
			IonsHalfSigma[i]*=world->GridScale;
		}
		if(IonsLJA!=NULL)
			for(i=0;i<NIonsTypes;i++)
		{
			IonsLJA[i]=float(double(IonsLJA[i])*GS12);
		}
		if(IonsLJB!=NULL)
			for(i=0;i<NIonsTypes;i++)
		{
			IonsLJB[i]=float(double(IonsLJB[i])*GS6);
		}
	}
	else
	{
		Rwat/=world->GridScale;
		Rsmooth/=world->GridScale;
		DiffRmaxAffect/=world->GridScale;
		for(i=0;i<DielNum;i++)
		{
			DielConst[i]*=EPKT;
		}
		for(i=0;i<NIonsTypes;i++)
		{
			IonsR[i]/=world->GridScale;
		}
		if(IonsHalfSigma!=NULL)
			for(i=0;i<NIonsTypes;i++)
		{
			IonsHalfSigma[i]/=world->GridScale;
		}
		if(IonsLJA!=NULL)
			for(i=0;i<NIonsTypes;i++)
		{
			IonsLJA[i]=float(double(IonsLJA[i])/GS12);
		}
		if(IonsLJB!=NULL)
			for(i=0;i<NIonsTypes;i++)
		{
			IonsLJB[i]=float(double(IonsLJB[i])/GS6);
		}
	}
	old_world=world;
	InternalUnit=ToInternal;
	return EXIT_SUCCESS;
}
int BuildWorldNI::BuildDiffExp(ContWorld* world)
{
	DbgPrint1("CGeometryWorld::BuildDiffExp()\n");
	int i,j,k;
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	float Dbulk;
	float BigDist=DiffRmaxAffect;
	float BigDistAngstrem=BigDist/world->GridScale;
	world->D=new float*[world->NIonsTypes];
	//world->C=new float*[world->NIonsTypes];
	
	//float r0=2.2;
	//float a=0.62;
	float rH=1.11*world->GridScale;
	
	int iElt;
	VectorField3D *dist=NULL;
	//if(SaveDistMap!="")
	//{
		dist=new VectorField3D(world->GridSizeGlobal, world->GridScale, world->NIonsTypes);
	//}
	
	for(i=0;i<world->NIonsTypes;i++)
	{
		DbgPrint1("Build Diffusion Map %d\n",i);
		
		
		for(j=0;j<world->GridSizeXYZGlobal;j++)dist->V[i][j]=BIGDISTANSE;
		
		for(iElt=0;iElt<GOElms.size();iElt++)
		{
			GOElms[iElt]->BuildDistMapsFromAtomsCenter(this,world,dist->V[i],BigDist, rH);
		}
		for(j=0;j<world->GridSizeXYZGlobal;j++)
			dist->V[i][j]/=world->GridScale;
		//replace all 0.0 values with 0 index
		world->D[i]=new float[world->GridSizeXYZGlobal];
		
		for(j=0;j<world->GridSizeXYZGlobal;j++)
		{
			Dbulk=world->NIndexing->GetDiffFloat(i,j);
			if(Dbulk>0.0f)
			{
				if(dist->V[i][j]<BigDistAngstrem)
				{
					if(dist->V[i][j]<Rexpdiff)
						world->D[i][j]=0.0;
					else
						world->D[i][j]=Dbulk*(1.0-exp(-Kexpdiff*(dist->V[i][j]-Rexpdiff)));
				}
				else
					world->D[i][j]=Dbulk;
			}
			else
			{
				world->D[i][j]=0.0;
			}
		}
	}
	if(SaveDistMap!="")
	{
		dist->WriteToFile(SaveDistMap.c_str(),1.0);
	}
	delete dist;
	return EXIT_SUCCESS;
}
int BuildWorldNI::BuildDiffExpOld(float K0, ContWorld* world)
{
	DbgPrint1("CGeometryWorld::BuildDiffExp()\n");
	int i,j,k;
	float Displacement[3]={0.0,0.0,0.0};
	float* Surf=new float[world->GridSizeXYZGlobal*4];
	int* iVtmp=new int[world->GridSizeXYZGlobal];
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	float Dbulk;
	float BigDist=7.0*world->GridScale;
	world->D=new float*[world->NIonsTypes];
	world->C=new float*[world->NIonsTypes];
	
	float r0=2.2*world->GridScale;
	float rH=1.11*world->GridScale;
	
	int iElt;
	VectorField3D *dist=NULL;
	if(SaveDistMap!="")
	{
		 dist=new VectorField3D(world->GridSizeGlobal, world->GridScale, world->NIonsTypes);
	}
	
	for(i=0;i<world->NIonsTypes;i++)
	{
		DbgPrint1("Build Diffusion Map %d\n",i);
		
		for(j=0;j<world->GridSizeXYZGlobal;j++)iVtmp[j]=BulkParam->IonsD[i];
		
		for(j=0;j<4*world->GridSizeXYZGlobal;j++)Surf[j]=-BIGDISTANSE;
		for(j=3;j<4*world->GridSizeXYZGlobal;j=j+4)Surf[j]=BIGDISTANSE;
		
		
		Displacement[0]=0.0;
		Displacement[1]=0.0;
		Displacement[2]=0.0;
		
		for(iElt=0;iElt<GOElms.size();iElt++)
		{
			if(GOElms[iElt]->IonsD[i]!=BulkParam->IonsD[i])
			{
				if(GOElms[iElt]->GetStdStrName()=="AtomsParameters")
				{
					GOAtoms* goa=(GOAtoms*)GOElms[iElt];
					bool makepreroll=goa->MakePreRoll;
					goa->MakePreRoll=false;
					GOElms[iElt]->BuildPreMaps(this,world,iVtmp,GOElms[iElt]->IonsD[i],BulkParam->IonsD[i],Displacement,0.0,DiffRmaxAffect,Surf,0);
					goa->MakePreRoll=makepreroll;
					
				}
				else
					GOElms[iElt]->BuildPreMaps(this,world,iVtmp,GOElms[iElt]->IonsD[i],BulkParam->IonsD[i],Displacement,Rwat,DiffRmaxAffect,Surf,0);
			}
			//GGO->BuildDistMaps(iVtmp,GGO->IonDiffusion[i],BulkGGO->IonDiffusion[i]+1,Displacement,GGO->IonRadius[i],BigDist,Surf);
		}
		if(dist!=NULL)
		{
			for(j=0;j<world->GridSizeXYZGlobal;j++)
				dist->V[i][j]=Surf[3+4*j];
		}
		for(j=0;j<world->GridSizeXYZGlobal;j++)
		{
			if(iVtmp[j]<0.0f)
				iVtmp[j]=BulkParam->IonsD[i];
		}
		//FinalazeSEVDiff(iVtmp,BulkGGO->IonDiffusion[i],RDsmooth,Surf);
			
		//replace all 0.0 values with 0 index
		world->D[i]=new float[world->GridSizeXYZGlobal];
		
		for(j=0;j<world->GridSizeXYZGlobal;j++)
		{
			Dbulk=world->NIndexing->GetDiffFloat(i,j);
			if(Dbulk>0.0f)
			{
				if(Surf[3+4*j]<BigDist)
					world->D[i][j]=Dbulk*(1.0-exp(-K0*Surf[3+4*j]/world->GridScale));
				else
					world->D[i][j]=Dbulk;
			}
			else
			{
				world->D[i][j]=0.0;
			}
			//if(iVtmp[j]>0)C[i][j]=BulkGGO->C[i]*coef;
		}
	}
	if(SaveDistMap!="")
	{
		dist->WriteToFile(SaveDistMap.c_str(),1.0);
	}
	delete [] iVtmp;
	delete [] Surf;
	
	return EXIT_SUCCESS;
}
#if defined(WITH_CUDA)
extern "C" int FinalazeSEVOnCUDA(int *GridSize,int *Field,int iBulkValue,float Rsmooth,float *Surf);
extern "C" int FinalazeSEVOnCUDA2(int *GridSize,int *Field,int iBulkValue,float Rsmooth,float3 *surf_points, int Nsurf_points);
extern "C" int BuildAtomsDielPreMapsOnCUDA(GOAtomsStruct* atms,float *Displ);
#endif
int BuildWorldNI::BuildDielMapsOnCuda(ContWorld* world,NodeIndexing *NIndexingNew)
{
#if defined(WITH_CUDA)
	DbgPrint1("BuildWorldNI::BuildDielMapsOnCuda\n");
	int i,j,k;
	float Displacement[3]={0.0,0.0,0.0};
	float* Surf=NULL;
	int* iVtmp=NULL;
	int iElt;
	
	if(Surf==NULL)
		Surf=new float[world->GridSizeXYZGlobal*4];
	if(iVtmp==NULL)
		iVtmp=new int[world->GridSizeXYZGlobal];
			
	for(i=0;i<DielNum;i++)DbgPrint1("DielConst[%d]=%f\n",i,DielConst[i]);
	for(i=0;i<3;i++)
	{
		DbgPrint1("Build Dielectric Map %d\n",i);
		for(j=0;j<world->GridSizeXYZGlobal;j++)
			iVtmp[j]=BulkParam->Epsilon;
				
		
		for(j=0;j<4*world->GridSizeXYZGlobal;j++)Surf[j]=-BIGDISTANSE;
		for(j=3;j<4*world->GridSizeXYZGlobal;j=j+4)Surf[j]=BIGDISTANSE;
				
		Displacement[0]=0.0;
		Displacement[1]=0.0;
		Displacement[2]=0.0;
		Displacement[i]=-0.5;
				
		DefClock0;
		StartClock0;
		for(iElt=0;iElt<GOElms.size();iElt++)
		{
			if(GOElms[iElt]->Epsilon!=BulkParam->Epsilon)
			{
				if(GOElms[iElt]->GetStdStrName()=="AtomsParameters")
				{
					GOAtoms *goatoms=(GOAtoms*)GOElms[iElt];
					int Natms=goatoms->NAtoms;
					
					GOAtomsStruct* cuda_atms=GOAtomsStruct_Create(world->GS_X,world->GS_Y,world->GS_Z,Natms,Rwat);
					
					cuda_atms->Surf=Surf;
					cuda_atms->iVtmp=iVtmp;
					cuda_atms->iDiel=goatoms->Epsilon;
					cuda_atms->iDielBulk=BulkParam->Epsilon;
					/*int count=0;
					for(j=0;j<Natms;j++)
					{
						if(goatoms->R[j]>0.0f)
						{
							cuda_atms->r[0][count]=goatoms->r[0][j];
							cuda_atms->r[1][count]=goatoms->r[1][j];
							cuda_atms->r[2][count]=goatoms->r[2][j];
							cuda_atms->R[count]=goatoms->R[j];
							count++;
						}
					}
					goatoms->NAtoms=count;
					Natms=goatoms->NAtoms;*/
					for(j=0;j<Natms;j++)
					{
						cuda_atms->r[0][j]=goatoms->r[0][j];
						cuda_atms->r[1][j]=goatoms->r[1][j];
						cuda_atms->r[2][j]=goatoms->r[2][j];
						cuda_atms->R[j]=goatoms->R[j];
					}
#if defined(WITH_CUDA)
					BuildAtomsDielPreMapsOnCUDA(cuda_atms,Displacement);
#endif
					//GOElms[iElt]->BuildPreMaps(this,world,iVtmp,GOElms[iElt]->Epsilon,BulkParam->Epsilon,Displacement,0.0,Rwat,Surf,0);
					cuda_atms=GOAtomsStruct_Delete(cuda_atms);
					
				}
				else
				{
					pnpPrint("Cannot do this element, can do only atoms\n");
				}
			}
		}
		StopClockWMes0("BuildPreMaps");
		StartClock0;
		
#if defined(WITH_CUDA)
		FinalazeSEVOnCUDA(world->GridSize,iVtmp,BulkParam->Epsilon,Rwat,Surf);
#endif
		//FinalazeSEV(world,iVtmp,BulkParam->Epsilon,Rwat,Surf);
		StopClockWMes0("FinalazeSEV");
				//Move indexes to 0 based
		for(j=0;j<world->GridSizeXYZGlobal;j++)
			iVtmp[j]--;
		switch(i)
		{
			case 0:
				NIndexingNew->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon0,NodeIndexing::Epsilon0Sft);
				break;
			case 1:
				NIndexingNew->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon1,NodeIndexing::Epsilon1Sft);
				break;
			case 2:
				NIndexingNew->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon2,NodeIndexing::Epsilon2Sft);
				break;
		}
				
				
	}
	NIndexingNew->CalcDielBoarder();
	DeleteCArray(iVtmp);
	DeleteCArray(Surf);
#endif
	return EXIT_SUCCESS;
}
int BuildWorldNI::BuildContWorld(ContWorld* world)
{
	DbgPrint1("BuildWorldNI::BuildContWorld()\n");
	DbgPrint1("MakeDiffusionMap=%d\n",(int)MakeDiffusionMap);
	int i,j,k;
	float Displacement[3]={0.0,0.0,0.0};
	float* Surf=NULL;
	int* iVtmp=NULL;
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	int iElt;
	
	//Change Parameters for internal 
	this->ChangeUnits(world,true);
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		GOElms[iElt]->ChangeUnits(world,true);
	}
	
	NodeIndexing *NIndexingOld,*NIndexingNew;
	
	NIndexingOld=world->NIndexing;
	world->NIndexing=new NodeIndexing();
	NIndexingNew=world->NIndexing;
	
	unsigned int GS[3]={world->GridSizeGlobal[0],world->GridSizeGlobal[1],world->GridSizeGlobal[2]};
	NIndexingNew->SetNNodes(GS,world->GridScale);
	NIndexingNew->NIonsTypes=world->NIonsTypes;
	NIndexingNew->IonsQ=new float[NIndexingNew->NIonsTypes];
	for(i=0;i<NIndexingNew->NIonsTypes;i++)
		NIndexingNew->IonsQ[i]=world->IonsQ[i];

	//Set Values Table
	
	for(i=0;i<DielNum;i++)NIndexingNew->Eps[i]=DielConst[i];
	
	
	for(i=0;i<DiffusionNum;i++)
	{
		NIndexingNew->D[i]=DiffusionConst[i];
		if(DiffusionConst[i]==0.0)NIndexingNew->C[i]=0.0;
		else NIndexingNew->C[i]=BulkParam->C[0]*coef;
	}
	//
	if(MakeDielectricMap)
	{
		for(i=0;i<DielNum;i++)
		{
			DbgPrint1("DielConst[%d]=%f\n",i,DielConst[i]);
		}
		
		if(iVtmp==NULL)
			iVtmp=new int[world->GridSizeXYZGlobal];
		
		//Check for memory limited calculations
		bool DoMemoryLimited=false;
		
		int NsubWorlds;//Number of sub worlds
		PartialWorldForSAS *SubWs=NULL;
		
		if(Surf!=NULL)
			DeleteCArray(Surf);
		
		if(MemoryLimitOnOneArray>0)
		{
			int locGS_X;
			int locGS_Y;
			int locGS_Z;
			
			locGS_X=world->GridSizeGlobal[0];
			locGS_Y=world->GridSizeGlobal[1];
			locGS_Z=world->GridSizeGlobal[2];
			
			int locGS_XY=locGS_X*locGS_Y;
			
			if( (locGS_XY*locGS_Z*4) <= (MemoryLimitOnOneArray*1024*1024) )
			{
				pnpPrint("No need to save memory, one temp. array will take %f MB\n",float(locGS_XY*locGS_Z*4)/1024.0/1024.0);
			}
			else
			{
				DoMemoryLimited=true;
				locGS_Z=MemoryLimitOnOneArray*1024*256/locGS_XY;
				//locGS_Z=world->GridSizeGlobal[2]/2;
				//locGS_Z=world->GridSizeGlobal[2];
				NsubWorlds=world->GridSizeGlobal[2]/locGS_Z;
				if(world->GridSizeGlobal[2]%locGS_Z!=0)NsubWorlds++;
				
				SubWs=new PartialWorldForSAS[NsubWorlds];
				Surf=new float[locGS_XY*locGS_Z*4];
				
				int isubw;
				for(isubw=0;isubw<NsubWorlds;isubw++)
				{
					SubWs[isubw].GS_X=world->GridSizeGlobal[0];
					SubWs[isubw].GS_Y=world->GridSizeGlobal[1];
					SubWs[isubw].GS_Z=world->GridSizeGlobal[2];
					SubWs[isubw].GS_XY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
					
					SubWs[isubw].locGS_Z=locGS_Z;
					if(isubw==NsubWorlds-1)
						SubWs[isubw].locGS_Z=world->GridSizeGlobal[2]-locGS_Z*(NsubWorlds-1);
					
					SubWs[isubw].locR0_Z=isubw*locGS_Z;
					SubWs[isubw].Field=iVtmp;
					SubWs[isubw].Surf=Surf;
					SubWs[isubw].Nsurf_points=0;
					SubWs[isubw].surf_points=NULL;
				}
				DbgPrint0("Split World into %d parts:\n",NsubWorlds);
				for(isubw=0;isubw<NsubWorlds;isubw++)
				{
					DbgPrint0("locGS=[%4d %4d %4d] ",SubWs[isubw].GS_X,SubWs[isubw].GS_Y,SubWs[isubw].locGS_Z);
					DbgPrint0("locR0Z=%4d\n",SubWs[isubw].locR0_Z);
				}
			}
			
		}
		
		if(!DoMemoryLimited)
		{
			Surf=new float[world->GridSizeXYZGlobal*4];
		}
		
		if(DoMemoryLimited)
		{
			for(i=0;i<3;i++)
			{
				DbgPrint1("Build Dielectric Map %d\n",i);
				for(j=0;j<world->GridSizeXYZGlobal;j++)
					iVtmp[j]=BulkParam->Epsilon;
				
				Displacement[0]=0.0;
				Displacement[1]=0.0;
				Displacement[2]=0.0;
				Displacement[i]=-0.5;
				
				DefClock0;
				StartClock0;
				
				int isubw;
				//Do preroll if needed
				for(iElt=0;iElt<GOElms.size();iElt++)
				{
					if(GOElms[iElt]->Epsilon!=BulkParam->Epsilon)
					{
						if(GOElms[iElt]->GetStdStrName()=="AtomsParameters")
						{
							GOAtoms *goatoms=(GOAtoms*)GOElms[iElt];
							if(goatoms->MakePreRoll)
							{
								DbgPrint0("Doing preroll\n");
								for(isubw=0;isubw<NsubWorlds;isubw++)
								{
									Displacement[0]=0.0;
									Displacement[1]=0.0;
									Displacement[2]=0.0;
									Displacement[i]=-0.5;
									Displacement[2]-=float(SubWs[isubw].locR0_Z);
									SubWs[isubw].Nsurf_points=0;
									DeleteCArray(SubWs[isubw].surf_points);
					
									int locXYZ=SubWs[isubw].GS_XY*SubWs[isubw].locGS_Z;
									for(j=0;j<4*locXYZ;j++)Surf[j]=-BIGDISTANSE;
									for(j=3;j<4*locXYZ;j=j+4)Surf[j]=BIGDISTANSE;
									
									goatoms->BuildPreMapsMemLim(&SubWs[isubw],GOElms[iElt]->Epsilon,BulkParam->Epsilon,Displacement,0.0,Rwat);
									//calc surf points
									for(j=0;j<locXYZ*4;j=j+4)
										if(Surf[j]>-100.0)
											SubWs[isubw].Nsurf_points++;
									
									DbgPrint0("Nsurf_points[%d]=%d\n",isubw,SubWs[isubw].Nsurf_points);
									
									if(SubWs[isubw].Nsurf_points>0)
									{
										SubWs[isubw].surf_points=new float3[SubWs[isubw].Nsurf_points];
										
										k=0;
										for(j=0;j<locXYZ*4;j=j+4)
											if(Surf[j]>-100.0)
										{
											SubWs[isubw].surf_points[k].x=Surf[j];
											SubWs[isubw].surf_points[k].y=Surf[j+1];
											SubWs[isubw].surf_points[k].z=Surf[j+2]+float(SubWs[isubw].locR0_Z);
											k++;
										}
									}
								}
								//
								if(BuildUsingGPU)
								{
#if defined(WITH_CUDA)
									for(isubw=0;isubw<NsubWorlds;isubw++)
									{
										if(SubWs[isubw].Nsurf_points>0)
											FinalazeSEVOnCUDA2(world->GridSize,iVtmp,BulkParam->Epsilon,Rwat,SubWs[isubw].surf_points,SubWs[isubw].Nsurf_points);
										SubWs[isubw].Nsurf_points=0;
										DeleteCArray(SubWs[isubw].surf_points);
									}
#else
									pnpWarning("Compiled without CUDA, will use CPU\n");
									for(isubw=0;isubw<NsubWorlds;isubw++)
									{
										if(SubWs[isubw].Nsurf_points>0)
											FinalazeSEV2(world,iVtmp,BulkParam->Epsilon,Rwat,SubWs[isubw].surf_points,SubWs[isubw].Nsurf_points);
										SubWs[isubw].Nsurf_points=0;
										DeleteCArray(SubWs[isubw].surf_points);
									}
#endif
								}
								else
								{
									for(isubw=0;isubw<NsubWorlds;isubw++)
									{
										if(SubWs[isubw].Nsurf_points>0)
											FinalazeSEV2(world,iVtmp,BulkParam->Epsilon,Rwat,SubWs[isubw].surf_points,SubWs[isubw].Nsurf_points);
										SubWs[isubw].Nsurf_points=0;
										DeleteCArray(SubWs[isubw].surf_points);
									}
								}
								for(j=0;j<world->GridSizeXYZGlobal;j++)
								{
									if(iVtmp[j]<0)
									{
										iVtmp[j]=-iVtmp[j];
									}
								}
// 								int *V[1]={iVtmp};
// 								VectorIntField3D VI(world->GridSize,world->GridScale,1,V);
// 								VI.WriteToFile("iVtmp.gz");
								DbgPrint0("Done prerolling\n");
							}
						}
					}
				}
				for(isubw=0;isubw<NsubWorlds;isubw++)
				{
					Displacement[0]=0.0;
					Displacement[1]=0.0;
					Displacement[2]=0.0;
					Displacement[i]=-0.5;
					Displacement[2]-=float(SubWs[isubw].locR0_Z);
					SubWs[isubw].Nsurf_points=0;
					DeleteCArray(SubWs[isubw].surf_points);
					
					int locXYZ=SubWs[isubw].GS_XY*SubWs[isubw].locGS_Z;
					for(j=0;j<4*locXYZ;j++)Surf[j]=-BIGDISTANSE;
					for(j=3;j<4*locXYZ;j=j+4)Surf[j]=BIGDISTANSE;
					
					for(iElt=0;iElt<GOElms.size();iElt++)
					{
						if(GOElms[iElt]->Epsilon!=BulkParam->Epsilon)
						{
							GOElms[iElt]->BuildPreMapsMemLim(&SubWs[isubw],GOElms[iElt]->Epsilon,BulkParam->Epsilon,Displacement,0.0,Rwat);
						}
							
					}
					//calc surf points
					for(j=0;j<locXYZ*4;j=j+4)
						if(Surf[j]>-100.0)
							SubWs[isubw].Nsurf_points++;
					
					DbgPrint0("Nsurf_points[%d]=%d\n",isubw,SubWs[isubw].Nsurf_points);
					
					if(SubWs[isubw].Nsurf_points>0)
					{
						SubWs[isubw].surf_points=new float3[SubWs[isubw].Nsurf_points];
						
						k=0;
						for(j=0;j<locXYZ*4;j=j+4)
							if(Surf[j]>-100.0)
						{
							SubWs[isubw].surf_points[k].x=Surf[j];
							SubWs[isubw].surf_points[k].y=Surf[j+1];
							SubWs[isubw].surf_points[k].z=Surf[j+2]+float(SubWs[isubw].locR0_Z);
							k++;
						}
					}
// 					for(k=0;k<SubWs[isubw].Nsurf_points;k++)
// 					{
// 						DbgPrint0("sp[%d]=%f %f %f\n",k,SubWs[isubw].surf_points[k].x,SubWs[isubw].surf_points[k].y,SubWs[isubw].surf_points[k].z);
// 					}
				}
				StopClockWMes0("BuildPreMaps");
				StartClock0;
				
				if(BuildUsingGPU)
				{
#if defined(WITH_CUDA)
					for(isubw=0;isubw<NsubWorlds;isubw++)
					{
						if(SubWs[isubw].Nsurf_points>0)
							FinalazeSEVOnCUDA2(world->GridSize,iVtmp,BulkParam->Epsilon,Rwat,SubWs[isubw].surf_points,SubWs[isubw].Nsurf_points);
						SubWs[isubw].Nsurf_points=0;
						DeleteCArray(SubWs[isubw].surf_points);
					}
#else
					pnpWarning("Compiled without CUDA, will use CPU\n");
					for(isubw=0;isubw<NsubWorlds;isubw++)
					{
						if(SubWs[isubw].Nsurf_points>0)
							FinalazeSEV2(world,iVtmp,BulkParam->Epsilon,Rwat,SubWs[isubw].surf_points,SubWs[isubw].Nsurf_points);
						SubWs[isubw].Nsurf_points=0;
						DeleteCArray(SubWs[isubw].surf_points);
					}
#endif
				}
				else
				{
					for(isubw=0;isubw<NsubWorlds;isubw++)
					{
						if(SubWs[isubw].Nsurf_points>0)
							FinalazeSEV2(world,iVtmp,BulkParam->Epsilon,Rwat,SubWs[isubw].surf_points,SubWs[isubw].Nsurf_points);
						SubWs[isubw].Nsurf_points=0;
						DeleteCArray(SubWs[isubw].surf_points);
					}
				}
				for(j=0;j<world->GridSizeXYZGlobal;j++)
				{
					if(iVtmp[j]<0)
					{
						iVtmp[j]=-iVtmp[j];
					}
				}
				//RemovingCavitiesOnDielectricMap
				if(RemovingCavitiesOnDielectricMap)
				{
					DbgPrint2("RemovingCavitiesOnDielectricMap\n");
					RemovingCavitiesAtValues(world->GridSize[0],world->GridSize[1],world->GridSize[2],iVtmp, RemCavOnDielWhere2Look ,RemCavOnDielFillWith);
				}
				StopClockWMes0("FinalazeSEV");
				//Move indexes to 0 based
				for(j=0;j<world->GridSizeXYZGlobal;j++)
					iVtmp[j]--;
				switch(i)
				{
					case 0:
						NIndexingNew->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon0,NodeIndexing::Epsilon0Sft);
						break;
					case 1:
						NIndexingNew->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon1,NodeIndexing::Epsilon1Sft);
						break;
					case 2:
						NIndexingNew->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon2,NodeIndexing::Epsilon2Sft);
						break;
				}
			}
			NIndexingNew->CalcDielBoarder();
			DeleteCArray(Surf);
		}
		else //i.e not Memory Limited
		{
			for(i=0;i<3;i++)
			{
				DbgPrint1("Build Dielectric Map %d\n",i);
				for(j=0;j<world->GridSizeXYZGlobal;j++)
					iVtmp[j]=BulkParam->Epsilon;
				
				
				for(j=0;j<4*world->GridSizeXYZGlobal;j++)Surf[j]=-BIGDISTANSE;
				for(j=3;j<4*world->GridSizeXYZGlobal;j=j+4)Surf[j]=BIGDISTANSE;
				
				Displacement[0]=0.0;
				Displacement[1]=0.0;
				Displacement[2]=0.0;
				Displacement[i]=-0.5;
				
				DefClock0;
				StartClock0;
				for(iElt=0;iElt<GOElms.size();iElt++)
				{
					if(GOElms[iElt]->Epsilon!=BulkParam->Epsilon)
						GOElms[iElt]->BuildPreMaps(this,world,iVtmp,GOElms[iElt]->Epsilon,BulkParam->Epsilon,Displacement,0.0,Rwat,Surf,0);
				}
				StopClockWMes0("BuildPreMaps");
				//StartClock0;
				
				if(BuildUsingGPU)
				{
#if defined(WITH_CUDA)
					FinalazeSEVOnCUDA(world->GridSize,iVtmp,BulkParam->Epsilon,Rwat,Surf);
#else
					pnpWarning("Compiled without CUDA, will use CPU\n");
					FinalazeSEV(world,iVtmp,BulkParam->Epsilon,Rwat,Surf);
#endif
				}
				else
					FinalazeSEV(world,iVtmp,BulkParam->Epsilon,Rwat,Surf);
				
				StopClockWMes0("FinalazeSEV");
				//Move indexes to 0 based
				for(j=0;j<world->GridSizeXYZGlobal;j++)
					iVtmp[j]--;
				switch(i)
				{
					case 0:
						NIndexingNew->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon0,NodeIndexing::Epsilon0Sft);
						break;
					case 1:
						NIndexingNew->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon1,NodeIndexing::Epsilon1Sft);
						break;
					case 2:
						NIndexingNew->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Epsilon2,NodeIndexing::Epsilon2Sft);
						break;
				}
			}
		}
		NIndexingNew->CalcDielBoarder();
		//StopClockWMes0("Make Diel Maps");
		DeleteCArray(SubWs);
	}
	int locNIonsTypes=NIonsTypes;
	if(MakeDiffusionMap==false&&MakeConcentrationMap==true)locNIonsTypes=1;//i.e maps for PB
	
	DbgPrint1("MakeDiffusionMap=%d locNIonsTypes=%d\n", (int)MakeDiffusionMap, locNIonsTypes);
	if(MakeDiffusionMap==true||MakeConcentrationMap==true)
		for(i=0;i<locNIonsTypes;i++)
	{
		if(!BuildUsingGPU)
		{
			if(Surf==NULL)
				Surf=new float[world->GridSizeXYZGlobal*4];
			if(iVtmp==NULL)
				iVtmp=new int[world->GridSizeXYZGlobal];
		}
		DbgPrint1("Build Diffusion Map %d\n",i);
		
		for(j=0;j<world->GridSizeXYZGlobal;j++)iVtmp[j]=BulkParam->IonsD[i];
		
		for(j=0;j<4*world->GridSizeXYZGlobal;j++)Surf[j]=-BIGDISTANSE;
		for(j=3;j<4*world->GridSizeXYZGlobal;j=j+4)Surf[j]=BIGDISTANSE;
		
		
		Displacement[0]=0.0;
		Displacement[1]=0.0;
		Displacement[2]=0.0;
		for(iElt=0;iElt<GOElms.size();iElt++)
		{
			
			if(GOElms[iElt]->IonsD[i]!=BulkParam->IonsD[i])
			{
				//DbgPrint1("GOElms[iElt]->IonsD[i]=%d, BulkParam->IonsD[i]=%d\n",GOElms[iElt]->IonsD[i],BulkParam->IonsD[i]);
				GOElms[iElt]->BuildPreMaps(this,world,iVtmp,GOElms[iElt]->IonsD[i],BulkParam->IonsD[i],Displacement,IonsR[i],Rsmooth,Surf,i+1);
			}
		}
		
		FinalazeSEVDiff(world,iVtmp,BulkParam->IonsD[i],Rsmooth,Surf);
		//FinalazeSEV(world,iVtmp,BulkParam->IonsD[i],Rwat,Surf);
		//Move indexes to 0 based
		for(j=0;j<world->GridSizeXYZGlobal;j++)
			iVtmp[j]--;
		
		switch(i)
		{
			case 0:
				NIndexingNew->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Ion0,NodeIndexing::Ion0Sft);
				break;
			case 1:
				NIndexingNew->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Ion1,NodeIndexing::Ion1Sft);
				break;
			case 2:
				NIndexingNew->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Ion2,NodeIndexing::Ion2Sft);
				break;
			case 3:
				NIndexingNew->SetIndexFieldFromIntArray(iVtmp,NodeIndexing::Ion3,NodeIndexing::Ion3Sft);
				break;
		}
	}
	if(MakeDiffusionMap==true&&MakeConcentrationMap==true)//e.i. for PNP
	{
		NIndexingNew->CalcDiffBoarder();
		//SetInitConcentrationFromNIndexing();
	}
	
	DeleteCArray(iVtmp);
	DeleteCArray(Surf);
	
	if(DiffusionMode==Exp)//Exp
		BuildDiffExp(world);
	//charges
	if(MakeChargeMap)
	{
		BuildContWorldCharges(world,NIndexingNew);
	}
	NIndexingNew->PBC[0]=world->PBC[0];
	NIndexingNew->PBC[1]=world->PBC[1];
	NIndexingNew->PBC[2]=world->PBC[2];
	//MakeLJRepultion
	if(MakeLJRepultion)
	{
		DbgPrint0("MakeLJRepultion\n");
		world->PMF=new float*[world->NIonsTypes];
		for(i=0;i<world->NIonsTypes;i++)
		{
			world->PMF[i]=new float[world->GridSizeXYZGlobal];
			for(j=0;j<world->GridSizeXYZGlobal;j++)
				world->PMF[i][j]=0.0;
			for(iElt=0;iElt<GOElms.size();iElt++)
			{
				GOElms[iElt]->BuildLJRepultionMap(this,world,i,world->PMF[i],LimitVlj);
			}
		}
	}
	if(AddElPotToPMF)
	{
		if(world->Potential!=NULL)
		{
		for(i=0;i<world->NIonsTypes;i++)
			for(j=0;j<world->GridSizeXYZGlobal;j++)
				world->PMF[i][j]+=world->Potential[j]*world->IonsQ[i];
		}
		else
		{
			pnpError("Can not add electrostatic potential because it is not initialized\n");
		}
	}
	//BoundaryCondition
	if(BoundaryCondition==CoulBC)
	{
		DefClock0;
		StartClock0;
		pnpPrint("Calculate Coulombic Boundary Conditions\n");
		pnpPrint("\tBldBCatPlane [%d %d %d]\n",BldBCatPlane[0],BldBCatPlane[1],BldBCatPlane[2]);
		if(world->Potential==NULL)
			world->Potential=new float[world->GridSizeXYZGlobal];
		for(i=0;i<world->GridSizeXYZGlobal;i++)
			world->Potential[i]=0.0;
		
		for(iElt=0;iElt<GOElms.size();iElt++)
		{
			if(GOElms[iElt]->GetStdStrName()=="AtomsParameters")
			{
				GOAtoms *goatoms=(GOAtoms*)GOElms[iElt];
				if(BuildUsingGPU)
				{
#if defined(WITH_CUDA)
					goatoms->SetBoundaryConditionOnCuda(this,world);
#else
					pnpWarning("Compiled without CUDA, will use CPU\n");
					goatoms->SetBoundaryCondition(this,world);
#endif
				}
				else
					goatoms->SetBoundaryCondition(this,world);
			}
			else
			{
				GOElms[iElt]->SetBoundaryCondition(this,world);
			}
		}
		StopClockWMes0("Calculation of Coulombic BC");
	}
	else
	{
		pnpPrint("Boundary Conditions: Zero\n");
	}
	//Change Parameters back to extenal then needed
	this->ChangeUnits(world,false);
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		GOElms[iElt]->ChangeUnits(world,false);
	}
	return EXIT_SUCCESS;
}
int BuildWorldNI::BuildContWorldCharges(ContWorld* world,NodeIndexing* NIndexingNew)
{
	int i,j,k,iElt;
	float* Qstat = new float[world->GridSizeXYZGlobal];
	for(i=0;i<world->GridSizeXYZGlobal;i++)Qstat[i]=0.0;
	
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		GOElms[iElt]->BuildPreMapsCharges(this,world,Qstat);
	}
	double Q=0.0;
	for(i=0;i<world->GridSizeXYZGlobal;i++)Q+=Qstat[i];
	DbgPrint1("Charge of Map: %.6g [e]\n",Q/(4.0*M_PI*world->GridScale));
	Q=0.0;
	for(i=1;i<world->GridSize[0]-1;i++)
		for(j=1;j<world->GridSize[1]-1;j++)
			for(k=1;k<world->GridSize[2]-1;k++)
	{
		Q+=Qstat[i+j*world->GridSize[0]+k*world->GridSize[0]*world->GridSize[1]];
	}
	DbgPrint1("Charge of Map: %.6g [e]\n",Q/(4.0*M_PI*world->GridScale));
	NIndexingNew->SetChargeMapFromArrayNoQonBoarder(Qstat);
	
	delete [] Qstat;
	
	return EXIT_SUCCESS;
}
int BuildWorldNI::FinalazeSEV(ContWorld* world,int *Field,int iBulkValue,float Rsmooth,float *Surf)
{
	DbgPrint2("BuildSES iBulkValue=%d\n",iBulkValue);
	int i,j,k,ix,iy,iz,gridpoint,rint[3];
	float RsmSQ,RtmpSQ;
	int *GridSize=world->GridSizeGlobal;
	int GS_X = world->GridSizeGlobal[0];
	int GS_Y = world->GridSizeGlobal[1];
	int GS_Z = world->GridSizeGlobal[2];
	int GSX=GridSize[0];
	int GSXY=GridSize[0]*GridSize[1];
	int GSXYZ=GridSize[0]*GridSize[1]*GridSize[2];
	int start[3];
	int end[3];
	int lim=(int)(Rsmooth+0.5);
	int itmp=GSXYZ*4;
	float *vtmp;
	
	RsmSQ=Rsmooth*Rsmooth;
	for(i=0;i<itmp;i=i+4)
	{
		if(Surf[i]>-100.0)
		{
			//DbgPrint2("BuildSES %d\n",i/4);
			vtmp=Surf+i;
			for(k=0;k<3;k++)
			{
				rint[k]=(int)(vtmp[k]+0.5);
				start[k]=rint[k]-lim;
				end[k]=rint[k]+lim;
				if(start[k]<0)start[k]=0;
				if(end[k]>GridSize[k]-1)end[k]=GridSize[k]-1;
			}
			for(ix=start[0];ix<=end[0];ix++)
				for(iy=start[1];iy<=end[1];iy++)
					for(iz=start[2];iz<=end[2];iz++)
			{
				gridpoint=ix+iy*GSX+iz*GSXY;
				RtmpSQ=(vtmp[0]-ix)*(vtmp[0]-ix)+(vtmp[1]-iy)*(vtmp[1]-iy)+(vtmp[2]-iz)*(vtmp[2]-iz);
				if(RtmpSQ<RsmSQ&&Field[gridpoint]<0)
				{
					Field[gridpoint]=iBulkValue;
				}
			}
		}
	}
	for(i=0;i<GSXYZ;i++)
	{
		if(Field[i]<0)
		{
			Field[i]=-Field[i];
		}
	}
	//RemovingCavitiesOnDielectricMap
	if(RemovingCavitiesOnDielectricMap)
	{
		DbgPrint2("RemovingCavitiesOnDielectricMap\n");
		RemovingCavitiesAtValues(GS_X,GS_Y,GS_Z,Field, RemCavOnDielWhere2Look ,RemCavOnDielFillWith);
	}
	return EXIT_SUCCESS;
}
int BuildWorldNI::FinalazeSEV2(ContWorld* world,int *Field,int iBulkValue,float Rsmooth,float3 *surf_points, int Nsurf_points)
{
	DbgPrint2("BuildSES iBulkValue=%d\n",iBulkValue);
	int i;
	int4 GS={world->GridSize[0],world->GridSize[1],world->GridSize[2],world->GridSize[0]*world->GridSize[1]};
	int GSXYZ=world->GridSize[0]*world->GridSize[1]*world->GridSize[2];
	int ix,iy,iz,gridpoint;
	float RsmSQ,RtmpSQ;
	
	int3 rint;
	int3 start;
	int3 end;
	
	RsmSQ=Rsmooth*Rsmooth;
	for(i=0;i<Nsurf_points;i++)
	{
		float RsmSQ=Rsmooth*Rsmooth;
		int lim=(int)(Rsmooth+0.5);
		
		rint.x=(int)(surf_points[i].x+0.5);
		rint.y=(int)(surf_points[i].y+0.5);
		rint.z=(int)(surf_points[i].z+0.5);
		
		start.x=rint.x-lim;
		start.y=rint.y-lim;
		start.z=rint.z-lim;
		
		end.x=rint.x+lim;
		end.y=rint.y+lim;
		end.z=rint.z+lim;
		
		if(start.x<0)start.x=0;
		if(start.y<0)start.y=0;
		if(start.z<0)start.z=0;
		
		if(end.x>GS.x-1)end.x=GS.x-1;
		if(end.y>GS.y-1)end.x=GS.y-1;
		if(end.z>GS.z-1)end.x=GS.z-1;
		
		for(ix=start.x;ix<=end.x;ix++)
			for(iy=start.y;iy<=end.y;iy++)
				for(iz=start.z;iz<=end.z;iz++)
		{
			gridpoint=ix+iy*GS.x+iz*GS.w;
			RtmpSQ=(surf_points[i].x-ix)*(surf_points[i].x-ix)+(surf_points[i].y-iy)*(surf_points[i].y-iy)+(surf_points[i].z-iz)*(surf_points[i].z-iz);
			if(RtmpSQ<RsmSQ&&Field[gridpoint]<0)
			{
				Field[gridpoint]=iBulkValue;
			}
		}
	}
	return EXIT_SUCCESS;
}
int BuildWorldNI::FinalazeSEVDiff(ContWorld* world,int *Field,int iBulkValue,float Rsmooth,float *Surf)
{
	DbgPrint2("BuildSES iBulkValue=%d\n",iBulkValue);
	int i,j,k,ix,iy,iz,gridpoint,rint[3];
	float RsmSQ,RtmpSQ;
	int GSX=world->GridSizeGlobal[0];
	int GSXY=world->GridSizeGlobal[0]*world->GridSizeGlobal[1];
	int GSXYZ=world->GridSizeGlobal[0]*world->GridSizeGlobal[1]*world->GridSizeGlobal[2];
	int start[3];
	int end[3];
	int lim=(int)(Rsmooth+0.5);
	int itmp=GSXYZ*4;
	float *vtmp;
	
	int GS_X = world->GridSizeGlobal[0];
	int GS_Y = world->GridSizeGlobal[1];
	int GS_Z = world->GridSizeGlobal[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_XY*GS_Z;
	int GrdPnt;
	int jgrid,kgrid;
	
	RsmSQ=Rsmooth*Rsmooth;
	for(i=0;i<itmp;i=i+4)
	{
		if(Surf[i]>-100.0)
		{
			//DbgPrint2("BuildSES %d\n",i/4);
			vtmp=Surf+i;
			for(k=0;k<3;k++)
			{
				rint[k]=(int)(vtmp[k]+0.5);
				start[k]=rint[k]-lim;
				end[k]=rint[k]+lim;
				if(start[k]<0)start[k]=0;
				if(end[k]>world->GridSizeGlobal[k]-1)end[k]=world->GridSizeGlobal[k]-1;
			}
			for(ix=start[0];ix<=end[0];ix++)
				for(iy=start[1];iy<=end[1];iy++)
					for(iz=start[2];iz<=end[2];iz++)
			{
				gridpoint=ix+iy*GSX+iz*GSXY;
				RtmpSQ=(vtmp[0]-ix)*(vtmp[0]-ix)+(vtmp[1]-iy)*(vtmp[1]-iy)+(vtmp[2]-iz)*(vtmp[2]-iz);
				if(RtmpSQ<RsmSQ&&Field[gridpoint]<0)
					Field[gridpoint]=iBulkValue;
			}
		}
	}
	for(i=0;i<GSXYZ;i++)
	{
		if(Field[i]<0)Field[i]=-Field[i];
	}
	//Set All D=0.0 to one index iD=0
	//replace all 0.0 values with iDzero index
	int iDzero1based;
	for(k=1;k<=DiffusionNum;k++)
	{
		if(DiffusionConst[k-1]==0.0)
		{
			iDzero1based=k;
			break;
		}
	}
	for(k=iDzero+1;k<=DiffusionNum;k++)
	{
		if(DiffusionConst[k-1]==0.0)
			for(j=0;j<world->GridSizeXYZGlobal;j++)
		{
			if(Field[j]==k)Field[j]=iDzero1based;
		}
	}
	//Removing Cavities,Box bourder has no Cavities
	if(RemovingCavities)
		RemovingCavitiesList(GS_X,GS_Y,GS_Z,Field, iDzero1based);
	//removing bad nodes
	int count=0;
	int Cycles=0;
	do
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
					if(Field[GrdPnt]>0)
					{
						if(Field[GrdPnt+1]==0 && Field[GrdPnt-1]==0 && Field[GrdPnt+GS_X]==0 && Field[GrdPnt-GS_X]==0 && Field[GrdPnt+GS_XY]==0 &&Field[GrdPnt-GS_XY]==0)
						{
							Field[GrdPnt]=iDzero1based;
							count++;
						}
					}
				}
			}
		}
		Cycles++;
		DbgPrint0("%d Cycle. Found %d	1-node cavity nodes, have removed it\n",Cycles,count);
	}
	while(count>0);
	DbgPrint0("Made %d Cycles for 1-node cavity removing\n",Cycles);
	return EXIT_SUCCESS;
}
int BuildWorldNI::RemovingCavitiesList(int GSX,int GSY,int GSZ,int *Field, int iDzero)
{
	int i,j,k,i0;
	int GSXY=GSX*GSY;
	int GSXYZ=GSXY*GSZ;
	int GrdPnt;
	int gridpx,gridmx;
	int gridpy,gridmy;
	int gridpz,gridmz;
	int count=0;
	int count0=0;
	int count1=8;
	int Cycles=0;
	int CyclesTot=0;
	int countNoCavities=0;

	for(i=1;i<GSX-1;i++)
		for(j=1;j<GSY-1;j++)
			for(k=1;k<GSZ-1;k++)
	{
		GrdPnt=i+j*GSX+k*GSXY;
		if(Field[GrdPnt]!=iDzero)Field[GrdPnt]=-Field[GrdPnt];
	}
	int *list=new int[GSXYZ];
	for(i=0;i<GSXYZ;i++)list[i]=-1;


	for(i=1;i<GSX-1;i++)
		for(j=1;j<GSY-1;j++)
			for(k=1;k<GSZ-1;k++)
	{
		k=1;
		GrdPnt=i+j*GSX+k*GSXY;
		if(Field[GrdPnt]!=iDzero)
		{
			list[count]=GrdPnt;
			Field[GrdPnt]=-Field[GrdPnt];
			count++;
		}
		k=GSZ-2;
		GrdPnt=i+j*GSX+k*GSXY;
		if(Field[GrdPnt]!=iDzero)
		{
			list[count]=GrdPnt;
			Field[GrdPnt]=-Field[GrdPnt];
			count++;
		}
	}
	count1=count;
	
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
			if(Field[gridpx]<0&&i<GSX-1)
			{
				Field[gridpx]=-Field[gridpx];
				list[count]=gridpx;
				count++;
			}
			if(Field[gridmx]<0&&i>1)
			{
				Field[gridmx]=-Field[gridmx];
				list[count]=gridmx;
				count++;
			}
			if(Field[gridpy]<0&&j<GSY-1)
			{
				Field[gridpy]=-Field[gridpy];
				list[count]=gridpy;
				count++;
			}
			if(Field[gridmy]<0&&j>1)
			{
				Field[gridmy]=-Field[gridmy];
				list[count]=gridmy;
				count++;
			}
			if(Field[gridpz]<0&&k<GSZ-1)
			{
				Field[gridpz]=-Field[gridpz];
				list[count]=gridpz;
				count++;
			}
			if(Field[gridmz]<0&&k>1)
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
	delete [] list;
	DbgPrint0("Made %d Cycles for cavity removing, %d nodes accesible\n",Cycles,count);
	for(i=0;i<GSXYZ;i++)
	{
		if(Field[i]<0)Field[i]=iDzero;
	}
	return EXIT_SUCCESS;
}
int BuildWorldNI::RemovingCavitiesAtValues(int GSX,int GSY,int GSZ,int *Field, int iWhereToLook, int iFillWith)
{
	int i,j,k,i0;
	int GSXY=GSX*GSY;
	int GSXYZ=GSXY*GSZ;
	int GrdPnt;
	int gridpx,gridmx;
	int gridpy,gridmy;
	int gridpz,gridmz;
	int count=1;
	int count0=0;
	int count1=1;
	int Cycles=0;
	int CyclesTot=0;
	int countNoCavities=0;

	for(i=1;i<GSX-1;i++)
		for(j=1;j<GSY-1;j++)
			for(k=1;k<GSZ-1;k++)
	{
		GrdPnt=i+j*GSX+k*GSXY;
		if(Field[GrdPnt]==iWhereToLook)Field[GrdPnt]=-Field[GrdPnt];
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
			if(Field[gridpx]<0&&i<GSX-1)
			{
				Field[gridpx]=-Field[gridpx];
				list[count]=gridpx;
				count++;
			}
			if(Field[gridmx]<0&&i>1)
			{
				Field[gridmx]=-Field[gridmx];
				list[count]=gridmx;
				count++;
			}
			if(Field[gridpy]<0&&j<GSY-1)
			{
				Field[gridpy]=-Field[gridpy];
				list[count]=gridpy;
				count++;
			}
			if(Field[gridmy]<0&&j>1)
			{
				Field[gridmy]=-Field[gridmy];
				list[count]=gridmy;
				count++;
			}
			if(Field[gridpz]<0&&k<GSZ-1)
			{
				Field[gridpz]=-Field[gridpz];
				list[count]=gridpz;
				count++;
			}
			if(Field[gridmz]<0&&k>1)
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
	delete [] list;
	DbgPrint0("Made %d Cycles for cavity removing, %d nodes accesible\n",Cycles,count);
	for(i=0;i<GSXYZ;i++)
	{
		if(Field[i]<0)Field[i]=iFillWith;
	}
	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
BuildWorldEu::BuildWorldEu()
{
	InitZero();
}


BuildWorldEu::~BuildWorldEu()
{
	Clear();
}

int BuildWorldEu::InitZero()
{
	BuildWorldNI::InitZero();
	PnpsObject::SetName("BuildWorldEu");
	TotalGridIds=-1;
	DielIdx=NULL;
	GGOIdx=NULL;
	DoNotChangeEpsFlag=1000000000;
	return EXIT_SUCCESS;
}
int BuildWorldEu::Clear()
{	
	BuildWorldNI::Clear();
	return EXIT_SUCCESS;
}
int BuildWorldEu::BuildContWorld(ContWorld* world)
{
	DbgPrint1("BuildWorldEu::BuildContWorld()\n");
	int i,j,k;
	int* iEpsilon[3];
	iEpsilon[0]=new int[world->GridSizeXYZGlobal];
	iEpsilon[1]=new int[world->GridSizeXYZGlobal];
	iEpsilon[2]=new int[world->GridSizeXYZGlobal];
	//int* iC=new int[world->GridSizeXYZGlobal];
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	int iElt;
	
	//Change Parameters for internal 
	this->ChangeUnits(world,true);
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		GOElms[iElt]->ChangeUnits(world,true);
	}
	
	
	NodeIndexing *NIndexingOld,*NIndexingNew;
	
	NIndexingOld=world->NIndexing;
	world->NIndexing=new NodeIndexing();
	NIndexingNew=world->NIndexing;
	
	unsigned int GS[3]={world->GridSizeGlobal[0],world->GridSizeGlobal[1],world->GridSizeGlobal[2]};
	NIndexingNew->SetNNodes(GS,world->GridScale);
	//Set Values Table
	
	for(i=0;i<DielNum;i++)NIndexingNew->Eps[i]=DielConst[i];
	
	
	for(i=0;i<DiffusionNum;i++)
	{
		NIndexingNew->D[i]=DiffusionConst[i];
		if(DiffusionConst[i]==0.0)NIndexingNew->C[i]=0.0;
		else NIndexingNew->C[i]=BulkParam->C[0]*coef;
	}
	
	for(i=0;i<DielNum;i++)DbgPrint1("DielConst[%d]=%f\n",i,DielConst[i]);
	DbgPrint1("Build Dielectric Map %d\n",i);
	
	//Check grid ID
	TotalGridIds=0;
	TotalGridIds=BulkParam->SetGridIDs(TotalGridIds);
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		TotalGridIds=GOElms[iElt]->SetGridIDs(TotalGridIds);
	}
	DielIdx=new int[TotalGridIds];
	GGOIdx=new GenericGeometricalObject*[TotalGridIds];
	DielIdx[0]=BulkParam->Epsilon;
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		for(i=GOElms[iElt]->ID_GridBegin;i<=GOElms[iElt]->ID_GridEnd;i++)
		{
			DielIdx[i]=GOElms[iElt]->Epsilon;
			GGOIdx[i]=GOElms[iElt];
		}
	}
	for(j=0;j<world->GridSizeXYZGlobal;j++)
	{
		iEpsilon[0][j]=BulkParam->ID_GridBegin;
		iEpsilon[1][j]=BulkParam->ID_GridBegin;
		iEpsilon[2][j]=BulkParam->ID_GridBegin;
	}
	DefClock0;
	StartClock0;
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		if(GOElms[iElt]->Epsilon!=BulkParam->Epsilon)
		{
			if(GOElms[iElt]->GetStdStrName()=="AtomsParameters")
				BuildPreMapsAtoms((GOAtoms*)GOElms[iElt],world,iEpsilon,Rwat);
			if(GOElms[iElt]->GetStdStrName()=="TubeParameters")
				BuildPreMapsTube((GOTube*)GOElms[iElt],world,iEpsilon,Rwat);
		}
	}
	StopClockWMes0("BuildPreMaps");
	StartClock0;
		
	FinalazeSEV(world,iEpsilon);
	StopClockWMes0("FinalazeSEV");
	//Move indexes to 0 based
	for(i=0;i<3;i++)
		for(j=0;j<world->GridSizeXYZGlobal;j++)
	{
		//if(iEpsilon[i][j]==0)iEpsilon[i][j]=BulkParam->Epsilon;
		//else iEpsilon[i][j]=goatoms->Epsilon;
		iEpsilon[i][j]--;
	}
	NIndexingNew->SetIndexFieldFromIntArray(iEpsilon[0],NodeIndexing::Epsilon0,NodeIndexing::Epsilon0Sft);
	NIndexingNew->SetIndexFieldFromIntArray(iEpsilon[1],NodeIndexing::Epsilon1,NodeIndexing::Epsilon1Sft);
	NIndexingNew->SetIndexFieldFromIntArray(iEpsilon[2],NodeIndexing::Epsilon2,NodeIndexing::Epsilon2Sft);
	NIndexingNew->CalcDielBoarder();
	
	delete [] iEpsilon[0];
	delete [] iEpsilon[1];
	delete [] iEpsilon[2];
	DeleteCArray(DielIdx);
	//delete [] iC;
	//charges
	if(MakeChargeMap)
	{
		float* Qstat = new float[world->GridSizeXYZGlobal];
		for(i=0;i<world->GridSizeXYZGlobal;i++)Qstat[i]=0.0;
	
		for(iElt=0;iElt<GOElms.size();iElt++)
		{
			GOElms[iElt]->BuildPreMapsCharges(this,world,Qstat);
		}
		double Q=0.0;
		for(i=0;i<world->GridSizeXYZGlobal;i++)Q+=Qstat[i];
		NIndexingNew->SetChargeMapFromArray(Qstat);
		delete [] Qstat;
		DbgPrint1("Charge of Map: %.6g [e]\n",Q/(4.0*M_PI*world->GridScale));
	}
	NIndexingNew->PBC[0]=world->PBC[0];
	NIndexingNew->PBC[1]=world->PBC[1];
	NIndexingNew->PBC[2]=world->PBC[2];
	
	
	if(BoundaryCondition==CoulBC)
	{
		pnpPrint("Calculate Coulombic Boundary Conditions\n");
		pnpPrint("\tBldBCatPlane [%d %d %d]\n",BldBCatPlane[0],BldBCatPlane[1],BldBCatPlane[2]);
		if(world->Potential==NULL)
			world->Potential=new float[world->GridSizeXYZGlobal];
		for(i=0;i<world->GridSizeXYZGlobal;i++)
			world->Potential[i]=0.0;
		
		for(iElt=0;iElt<GOElms.size();iElt++)
		{
			GOElms[iElt]->SetBoundaryCondition(this,world);
		}
	}
	else
	{
		pnpPrint("Boundary Conditions: Zero\n");
	}
	//Change Parameters back to extenal then needed
	this->ChangeUnits(world,false);
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		GOElms[iElt]->ChangeUnits(world,false);
	}
	return EXIT_SUCCESS;
}
int BuildWorldEu::BuildPreMapsAtoms(GOAtoms *goatoms,ContWorld* world,int **iEpsilon,float Rpr)
{
	DbgPrint0("BuildWorldEu::BuildPreMapsAtoms\n");
	int GrdPnt2;
	int i,j,k;
	int ix,iy,iz;
	int ix2,iy2,iz2;
	float Rion=0.0;
	int GS_X=world->GS_X;
	int GS_Y=world->GS_Y;
	int GS_Z=world->GS_Z;
	int GS_XY=world->GS_XY;
	int ID_Start=goatoms->ID_GridBegin;
	float d1,d2;
	
	int iatom;
	float rcx,rcy,rcz;
	float fx2,fy2,fz2;
	for(iatom=0;iatom<goatoms->NAtoms;iatom++)
	{
		float R=goatoms->R[iatom]+Rion;
		
		if(R==0.0)continue;
		
		rcx=goatoms->r[0][iatom];
		rcy=goatoms->r[1][iatom];
		rcz=goatoms->r[2][iatom];
		ix=roundf(rcx);
		iy=roundf(rcy);
		iz=roundf(rcz);
		
		float RpRpr=R+Rpr;
		float RpRprSQ=RpRpr*RpRpr;
		float RpRprEpsSQmQ=RpRprSQ-0.25;
		float RSQ=R*R;
		float RepsSQmQ=RSQ-0.25;
		int iRmax=int(float(RpRpr+1.5));
		
		int startX=ix-iRmax;
		int startY=iy-iRmax;
		int startZ=iz-iRmax;
		int endX=ix+iRmax;
		int endY=iy+iRmax;
		int endZ=iz+iRmax;
		if(startX<0)startX=0;
		if(startY<0)startY=0;
		if(startZ<0)startZ=0;
		if(endX>=GS_X)endX=GS_X-1;
		if(endY>=GS_Y)endY=GS_Y-1;
		if(endZ>=GS_Z)endZ=GS_Z-1;
		
		for(ix2=startX;ix2<=endX;ix2++)
			for(iy2=startY;iy2<=endY;iy2++)
				for(iz2=startZ;iz2<=endZ;iz2++)
		{
			fx2=(float)ix2;
			fy2=(float)iy2;
			fz2=(float)iz2;
			float drx=fx2-rcx;
			float dry=fy2-rcy;
			float drz=fz2-rcz;
			float drSQ=drx*drx+dry*dry+drz*drz;
			float drSQx=drSQ+drx;
			float drSQy=drSQ+dry;
			float drSQz=drSQ+drz;
			GrdPnt2=ix2+iy2*GS_X+iz2*GS_XY;
			
			if(drSQx<RepsSQmQ)
			{
				iEpsilon[0][GrdPnt2]=iatom+ID_Start+DoNotChangeEpsFlag;
			}
			else if(drSQx<RpRprEpsSQmQ)
			{
				if(iEpsilon[0][GrdPnt2]==BulkParam->ID_GridBegin)
				{
					iEpsilon[0][GrdPnt2]=iatom+ID_Start;
				}
				else if(iEpsilon[0][GrdPnt2]<DoNotChangeEpsFlag)
				{
					d1=GGOIdx[iatom+ID_Start]->GetDistToSurf(fx2,fy2,fz2,Rion,Rpr,iatom+ID_Start);
					d2=GGOIdx[iEpsilon[0][GrdPnt2]]->GetDistToSurf(fx2,fy2,fz2,Rion,Rpr,iEpsilon[0][GrdPnt2]);
					if(d1<d2)
						iEpsilon[0][GrdPnt2]=iatom+ID_Start;
				}
			}
			if(drSQy<RepsSQmQ)
			{
				iEpsilon[1][GrdPnt2]=iatom+ID_Start+DoNotChangeEpsFlag;
			}
			else if(drSQy<RpRprEpsSQmQ)
			{
				if(iEpsilon[1][GrdPnt2]==BulkParam->ID_GridBegin)
				{
					iEpsilon[1][GrdPnt2]=iatom+ID_Start;
				}
				else if(iEpsilon[1][GrdPnt2]<DoNotChangeEpsFlag)
				{
					d1=GGOIdx[iatom+ID_Start]->GetDistToSurf(fx2,fy2,fz2,Rion,Rpr,iatom+ID_Start);
					d2=GGOIdx[iEpsilon[1][GrdPnt2]]->GetDistToSurf(fx2,fy2,fz2,Rion,Rpr,iEpsilon[1][GrdPnt2]);
					if(d1<d2)
						iEpsilon[1][GrdPnt2]=iatom+ID_Start;
				}
			}
			if(drSQz<RepsSQmQ)
			{
				iEpsilon[2][GrdPnt2]=iatom+ID_Start+DoNotChangeEpsFlag;
			}
			else if(drSQz<RpRprEpsSQmQ)
			{
				if(iEpsilon[2][GrdPnt2]==BulkParam->ID_GridBegin)
				{
					iEpsilon[2][GrdPnt2]=iatom+ID_Start;
				}
				else if(iEpsilon[2][GrdPnt2]<DoNotChangeEpsFlag)
				{
					d1=GGOIdx[iatom+ID_Start]->GetDistToSurf(fx2,fy2,fz2,Rion,Rpr,iatom+ID_Start);
					d2=GGOIdx[iEpsilon[2][GrdPnt2]]->GetDistToSurf(fx2,fy2,fz2,Rion,Rpr,iEpsilon[2][GrdPnt2]);
					if(d1<d2)
						iEpsilon[2][GrdPnt2]=iatom+ID_Start;
				}
			}
		}
	}
	
	return EXIT_SUCCESS;
}
int BuildWorldEu::BuildPreMapsTube(GOTube *gt,ContWorld* world,int **iEpsilon,float Rpr)
{
	DbgPrint0("BuildWorldEu::BuildPreMapsAtoms\n");
	int GrdPnt;
	int i,j,k;
	int ix,iy,iz;
	float radmax2=0.0;
	float Rion=0.0;
	int GS_X=world->GS_X;
	int GS_Y=world->GS_Y;
	int GS_Z=world->GS_Z;
	int GS_XY=world->GS_XY;
	int GS_XYZ=world->GS_XYZ;
	int ID_Start=gt->ID_GridBegin;
	
	gt->tX=gt->XY[0];
	gt->tY=gt->XY[1];
	gt->tZ1=gt->Z[0]-Rion;
	gt->tZ2=gt->Z[1]+Rion;
	gt->tR1=gt->R[0]-Rion;
	gt->tR2=gt->R[1]+Rion;
	
	gt->tZ1mRpr=gt->tZ1-Rpr;
	gt->tZ2pRpr=gt->tZ2+Rpr;
	
	float fx,fy,fz;
	
	for(ix=0;ix<GS_X;ix++)
		for(iy=0;iy<GS_Y;iy++)
			for(iz=0;iz<GS_Z;iz++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		fx=(float)ix;
		fy=(float)iy;
		fz=(float)iz;
		if((gt->tZ1mRpr<fz)&&(fz<gt->tZ1))
		{
			iEpsilon[0][GrdPnt]=ID_Start;
			iEpsilon[1][GrdPnt]=ID_Start;
			iEpsilon[2][GrdPnt]=ID_Start;
		}
		else if((gt->tZ1<=fz)&&(fz<=gt->tZ2))
		{
			iEpsilon[0][GrdPnt]=ID_Start+DoNotChangeEpsFlag;
			iEpsilon[1][GrdPnt]=ID_Start+DoNotChangeEpsFlag;
			iEpsilon[2][GrdPnt]=ID_Start+DoNotChangeEpsFlag;
		}
		else if((gt->tZ2<fz)&&(fz<gt->tZ2pRpr))
		{
			iEpsilon[0][GrdPnt]=ID_Start;
			iEpsilon[1][GrdPnt]=ID_Start;
			iEpsilon[2][GrdPnt]=ID_Start;
		}
	}
	return EXIT_SUCCESS;
}
int BuildWorldEu::FinalazeSEV(ContWorld* world,int **iEpsilon)
{
	DbgPrint0("BuildWorldEu::FinalazeSEV\n");
	int GrdPnt,GrdPnt2;
	int i,j,k;
	int ix,iy,iz;
	int ix2,iy2,iz2;
	int GS_X=world->GS_X;
	int GS_Y=world->GS_Y;
	int GS_Z=world->GS_Z;
	int GS_XY=world->GS_XY;
	int GS_XYZ=world->GS_XYZ;
	
	
	float Rpr=Rwat;
	float RprSQ=Rpr*Rpr;
	float RepsSQmQ=RprSQ-0.25;
	int iRmax=int(float(Rwat+1.5));
	
	
	int GridIDs[6];
	float dr[3];
	float fx,fy,fz;
	float fx2,fy2,fz2;
	float d;
	int count0=0;
	int iEps;
	for(iEps=0;iEps<3;iEps++)
	{
		for(ix=1;ix<GS_X-1;ix++)
			for(iy=1;iy<GS_Y-1;iy++)
				for(iz=1;iz<GS_Z-1;iz++)
		{
			GrdPnt=ix+iy*GS_X+iz*GS_XY;
			fx=(float)ix;
			fy=(float)iy;
			fz=(float)iz;
			GridIDs[0]=abs(iEpsilon[iEps][GrdPnt])%DoNotChangeEpsFlag;
			GridIDs[1]=abs(iEpsilon[iEps][GrdPnt-1])%DoNotChangeEpsFlag;
			GridIDs[2]=abs(iEpsilon[iEps][GrdPnt-GS_X])%DoNotChangeEpsFlag;
			GridIDs[3]=abs(iEpsilon[iEps][GrdPnt-GS_XY])%DoNotChangeEpsFlag;
			GridIDs[4]=abs(iEpsilon[iEps][GrdPnt+1])%DoNotChangeEpsFlag;
			GridIDs[5]=abs(iEpsilon[iEps][GrdPnt+GS_X])%DoNotChangeEpsFlag;
			GridIDs[6]=abs(iEpsilon[iEps][GrdPnt+GS_XY])%DoNotChangeEpsFlag;
			int sum=DielIdx[GridIDs[0]] + DielIdx[GridIDs[1]] + DielIdx[GridIDs[2]] +  DielIdx[GridIDs[3]]+DielIdx[GridIDs[4]] + DielIdx[GridIDs[5]] + DielIdx[GridIDs[6]];
			//if((sum!=4*DielIdx[GridIDs[0]]))
			if(GridIDs[0]!=BulkParam->ID_GridBegin&&sum!=7*DielIdx[GridIDs[0]])
			{
				
				count0++;
				int GridID=GridIDs[0];
				/*float dmax=1000000.0;
				for(i=0;i<4;i++)
				{
					if(GridIDs[i]!=BulkParam->ID_GridBegin)
					{
						d=GGOIdx[GridIDs[i]]->GetDistToSurf(fx,fy,fz,0.0,Rwat,GridIDs[i]);
						if(d<dmax)
						{
							GridID=GridIDs[i];
							dmax=d;
						}
					}
				}*/
				
				GGOIdx[GridID]->GetClosestSurfPoint(&fx,&fy,&fz,0.0,Rwat,GridID);
				
				int startX=ix-iRmax;
				int startY=iy-iRmax;
				int startZ=iz-iRmax;
				int endX=ix+iRmax;
				int endY=iy+iRmax;
				int endZ=iz+iRmax;
				if(startX<0)startX=0;
				if(startY<0)startY=0;
				if(startZ<0)startZ=0;
				if(endX>=GS_X)endX=GS_X-1;
				if(endY>=GS_Y)endY=GS_Y-1;
				if(endZ>=GS_Z)endZ=GS_Z-1;
				
				for(ix2=startX;ix2<=endX;ix2++)
					for(iy2=startY;iy2<=endY;iy2++)
						for(iz2=startZ;iz2<=endZ;iz2++)
				{
					fx2=(float)ix2;
					fy2=(float)iy2;
					fz2=(float)iz2;
					dr[0]=fx2-fx;
					dr[1]=fy2-fy;
					dr[2]=fz2-fz;
					float drSQ=dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
					GrdPnt2=ix2+iy2*GS_X+iz2*GS_XY;
					if((drSQ+dr[iEps]<RepsSQmQ)&&(iEpsilon[iEps][GrdPnt2]<DoNotChangeEpsFlag))
						iEpsilon[iEps][GrdPnt2]=-abs(iEpsilon[iEps][GrdPnt2]);
				}
			}
		}
	}
	DbgPrint0("BigBoarders %d\n",count0);
	for(i=0;i<3;i++)
		for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
	{
		if(iEpsilon[i][GrdPnt]<0)iEpsilon[i][GrdPnt]=BulkParam->ID_GridBegin;
		iEpsilon[i][GrdPnt]=DielIdx[iEpsilon[i][GrdPnt]%DoNotChangeEpsFlag];
	}
	return EXIT_SUCCESS;
}
int BuildWorldEu::FinalazeSEV2(ContWorld* world,int **iEpsilon)
{
	DbgPrint0("BuildWorldEu::FinalazeSEV\n");
	int GrdPnt,GrdPnt2;
	int i,j,k;
	int ix,iy,iz;
	int ix2,iy2,iz2;
	int GS_X=world->GS_X;
	int GS_Y=world->GS_Y;
	int GS_Z=world->GS_Z;
	int GS_XY=world->GS_XY;
	int GS_XYZ=world->GS_XYZ;
	
	
	float Rpr=Rwat;
	float RprSQ=Rpr*Rpr;
	float RepsSQmQ=RprSQ-0.25;
	int iRmax=int(float(Rwat+1.5));
	
	
	int GridIDs[6];
	float rcx,rcy,rcz;
	float fx,fy,fz;
	float fx2,fy2,fz2;
	float d;
	int count0=0;
	int iatom;
	for(ix=1;ix<GS_X-1;ix++)
		for(iy=1;iy<GS_Y-1;iy++)
			for(iz=1;iz<GS_Z-1;iz++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		fx=(float)ix;
		fy=(float)iy;
		fz=(float)iz;
		GridIDs[0]=abs(iEpsilon[0][GrdPnt])%DoNotChangeEpsFlag;
		GridIDs[1]=abs(iEpsilon[1][GrdPnt])%DoNotChangeEpsFlag;
		GridIDs[2]=abs(iEpsilon[2][GrdPnt])%DoNotChangeEpsFlag;
		GridIDs[3]=abs(iEpsilon[0][GrdPnt-1])%DoNotChangeEpsFlag;
		GridIDs[4]=abs(iEpsilon[1][GrdPnt-GS_X])%DoNotChangeEpsFlag;
		GridIDs[5]=abs(iEpsilon[2][GrdPnt-GS_XY])%DoNotChangeEpsFlag;
		
		int sum=DielIdx[GridIDs[0]] + DielIdx[GridIDs[1]] + DielIdx[GridIDs[2]] +  DielIdx[GridIDs[3]] + DielIdx[GridIDs[4]] + DielIdx[GridIDs[5]];
		if(sum!=6*DielIdx[GridIDs[0]])
		{
			count0++;
			int GridID=GridIDs[0];
			float dmax=1000000.0;
			for(i=0;i<6;i++)
			{
				if(GridIDs[i]!=BulkParam->ID_GridBegin)
				{
					d=GGOIdx[GridIDs[i]]->GetDistToSurf(fx,fy,fz,0.0,Rwat,GridIDs[i]);
					if(d<dmax)
					{
						GridID=GridIDs[i];
						dmax=d;
					}
				}
			}
			GGOIdx[GridID]->GetClosestSurfPoint(&fx,&fy,&fz,0.0,Rwat,GridID);
			
			int startX=ix-iRmax;
			int startY=iy-iRmax;
			int startZ=iz-iRmax;
			int endX=ix+iRmax;
			int endY=iy+iRmax;
			int endZ=iz+iRmax;
			if(startX<0)startX=0;
			if(startY<0)startY=0;
			if(startZ<0)startZ=0;
			if(endX>=GS_X)endX=GS_X-1;
			if(endY>=GS_Y)endY=GS_Y-1;
			if(endZ>=GS_Z)endZ=GS_Z-1;
			
			for(ix2=startX;ix2<=endX;ix2++)
				for(iy2=startY;iy2<=endY;iy2++)
					for(iz2=startZ;iz2<=endZ;iz2++)
			{
				fx2=(float)ix2;
				fy2=(float)iy2;
				fz2=(float)iz2;
				float drx=fx2-fx;
				float dry=fy2-fy;
				float drz=fz2-fz;
				float drSQ=drx*drx+dry*dry+drz*drz;
				GrdPnt2=ix2+iy2*GS_X+iz2*GS_XY;
				if((drSQ+drx<RepsSQmQ)&&(iEpsilon[0][GrdPnt2]<DoNotChangeEpsFlag))
					iEpsilon[0][GrdPnt2]=-abs(iEpsilon[0][GrdPnt2]);
				if((drSQ+dry<RepsSQmQ)&&(iEpsilon[0][GrdPnt2]<DoNotChangeEpsFlag))
					iEpsilon[1][GrdPnt2]=-abs(iEpsilon[1][GrdPnt2]);
				if((drSQ+drz<RepsSQmQ)&&(iEpsilon[0][GrdPnt2]<DoNotChangeEpsFlag))
					iEpsilon[2][GrdPnt2]=-abs(iEpsilon[2][GrdPnt2]);
			}
		}
	}
	DbgPrint0("BigBoarders %d\n",count0);
	for(i=0;i<3;i++)
		for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
	{
		if(iEpsilon[i][GrdPnt]<0)iEpsilon[i][GrdPnt]=BulkParam->ID_GridBegin;
		iEpsilon[i][GrdPnt]=DielIdx[iEpsilon[i][GrdPnt]%DoNotChangeEpsFlag];
	}
	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
BuildWorldCmp::BuildWorldCmp()
{
	InitZero();
}


BuildWorldCmp::~BuildWorldCmp()
{
	Clear();
}

int BuildWorldCmp::InitZero()
{
	BuildWorldNI::InitZero();
	PnpsObject::SetName("BuildWorldEu");
	TotalGridIds=-1;
	DielIdx=NULL;
	GGOIdx=NULL;
	return EXIT_SUCCESS;
}
int BuildWorldCmp::Clear()
{	
	BuildWorldNI::Clear();
	return EXIT_SUCCESS;
}
int BuildWorldCmp::BuildContWorld(ContWorld* world)
{
	DbgPrint1("BuildWorldEu::BuildContWorld()\n");
	int i,j,k;
	int* iEpsilon[3];
	iEpsilon[0]=new int[world->GridSizeXYZGlobal];
	iEpsilon[1]=new int[world->GridSizeXYZGlobal];
	iEpsilon[2]=new int[world->GridSizeXYZGlobal];
	//int* iC=new int[world->GridSizeXYZGlobal];
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	int iElt;
	
	//Change Parameters for internal 
	this->ChangeUnits(world,true);
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		GOElms[iElt]->ChangeUnits(world,true);
	}
	
	
	NodeIndexing *NIndexingOld,*NIndexingNew;
	
	NIndexingOld=world->NIndexing;
	world->NIndexing=new NodeIndexing();
	NIndexingNew=world->NIndexing;
	
	unsigned int GS[3]={world->GridSizeGlobal[0],world->GridSizeGlobal[1],world->GridSizeGlobal[2]};
	NIndexingNew->SetNNodes(GS,world->GridScale);
	//Set Values Table
	
	for(i=0;i<DielNum;i++)NIndexingNew->Eps[i]=DielConst[i];
	
	
	for(i=0;i<DiffusionNum;i++)
	{
		NIndexingNew->D[i]=DiffusionConst[i];
		if(DiffusionConst[i]==0.0)NIndexingNew->C[i]=0.0;
		else NIndexingNew->C[i]=BulkParam->C[0]*coef;
	}
	
	for(i=0;i<DielNum;i++)DbgPrint1("DielConst[%d]=%f\n",i,DielConst[i]);
	DbgPrint1("Build Dielectric Map %d\n",i);
	
	//Check grid ID
	TotalGridIds=0;
	TotalGridIds=BulkParam->SetGridIDs(TotalGridIds);
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		TotalGridIds=GOElms[iElt]->SetGridIDs(TotalGridIds);
	}
	DielIdx=new int[TotalGridIds];
	GGOIdx=new GenericGeometricalObject*[TotalGridIds];
	DielIdx[0]=BulkParam->Epsilon;
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		for(i=GOElms[iElt]->ID_GridBegin;i<=GOElms[iElt]->ID_GridEnd;i++)
		{
			DielIdx[i]=GOElms[iElt]->Epsilon;
			GGOIdx[i]=GOElms[iElt];
		}
	}
	for(j=0;j<world->GridSizeXYZGlobal;j++)
	{
		iEpsilon[0][j]=BulkParam->ID_GridBegin;
		iEpsilon[1][j]=BulkParam->ID_GridBegin;
		iEpsilon[2][j]=BulkParam->ID_GridBegin;
	}
	DefClock0;
	StartClock0;
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		if(GOElms[iElt]->Epsilon!=BulkParam->Epsilon)
		{
			if(GOElms[iElt]->GetStdStrName()=="AtomsParameters")
				BuildPreMapsAtoms((GOAtoms*)GOElms[iElt],world,iEpsilon,Rwat);
		}
	}
	StopClockWMes0("BuildPreMaps");
	StartClock0;
		
	FinalazeSEV(world,iEpsilon);
	StopClockWMes0("FinalazeSEV");
	//Move indexes to 0 based
	for(i=0;i<3;i++)
		for(j=0;j<world->GridSizeXYZGlobal;j++)
	{
		//if(iEpsilon[i][j]==0)iEpsilon[i][j]=BulkParam->Epsilon;
		//else iEpsilon[i][j]=goatoms->Epsilon;
		iEpsilon[i][j]--;
	}
	NIndexingNew->SetIndexFieldFromIntArray(iEpsilon[0],NodeIndexing::Epsilon0,NodeIndexing::Epsilon0Sft);
	NIndexingNew->SetIndexFieldFromIntArray(iEpsilon[1],NodeIndexing::Epsilon1,NodeIndexing::Epsilon1Sft);
	NIndexingNew->SetIndexFieldFromIntArray(iEpsilon[2],NodeIndexing::Epsilon2,NodeIndexing::Epsilon2Sft);
	NIndexingNew->CalcDielBoarder();
	
	delete [] iEpsilon[0];
	delete [] iEpsilon[1];
	delete [] iEpsilon[2];
	DeleteCArray(DielIdx);
	//delete [] iC;
	//charges
	if(MakeChargeMap)
	{
		float* Qstat = new float[world->GridSizeXYZGlobal];
		for(i=0;i<world->GridSizeXYZGlobal;i++)Qstat[i]=0.0;
	
		for(iElt=0;iElt<GOElms.size();iElt++)
		{
			GOElms[iElt]->BuildPreMapsCharges(this,world,Qstat);
		}
		double Q=0.0;
		for(i=0;i<world->GridSizeXYZGlobal;i++)Q+=Qstat[i];
		NIndexingNew->SetChargeMapFromArray(Qstat);
		delete [] Qstat;
		DbgPrint1("Charge of Map: %.6g [e]\n",Q/(4.0*M_PI*world->GridScale));
	}
	NIndexingNew->PBC[0]=world->PBC[0];
	NIndexingNew->PBC[1]=world->PBC[1];
	NIndexingNew->PBC[2]=world->PBC[2];
	
	
	if(BoundaryCondition==CoulBC)
	{
		pnpPrint("Calculate Coulombic Boundary Conditions\n");
		pnpPrint("\tBldBCatPlane [%d %d %d]\n",BldBCatPlane[0],BldBCatPlane[1],BldBCatPlane[2]);
		if(world->Potential==NULL)
			world->Potential=new float[world->GridSizeXYZGlobal];
		for(i=0;i<world->GridSizeXYZGlobal;i++)
			world->Potential[i]=0.0;
		
		for(iElt=0;iElt<GOElms.size();iElt++)
		{
			GOElms[iElt]->SetBoundaryCondition(this,world);
		}
	}
	else
	{
		pnpPrint("Boundary Conditions: Zero\n");
	}
	//Change Parameters back to extenal then needed
	this->ChangeUnits(world,false);
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		GOElms[iElt]->ChangeUnits(world,false);
	}
	return EXIT_SUCCESS;
}
int BuildWorldCmp::BuildPreMapsAtoms(GOAtoms *goatoms,ContWorld* world,int **iEpsilon,float Rpr)
{
	DbgPrint0("BuildWorldCmp::BuildPreMapsAtoms\n");
	int GrdPnt2;
	int i,j,k;
	int ix,iy,iz;
	int ix2,iy2,iz2;
	float Rion=0.0;
	int GS_X=world->GS_X;
	int GS_Y=world->GS_Y;
	int GS_Z=world->GS_Z;
	int GS_XY=world->GS_XY;
	int ID_Start=goatoms->ID_GridBegin;
	float d1,d2;
	
	int iatom;
	float rcx,rcy,rcz;
	float fx2,fy2,fz2;
	for(iatom=0;iatom<goatoms->NAtoms;iatom++)
	{
		float R=goatoms->R[iatom]+Rion;
		
		if(R==0.0)continue;
		
		rcx=goatoms->r[0][iatom];
		rcy=goatoms->r[1][iatom];
		rcz=goatoms->r[2][iatom];
		ix=roundf(rcx);
		iy=roundf(rcy);
		iz=roundf(rcz);
		
		float RpRpr=R+Rpr;
		float RpRprSQ=RpRpr*RpRpr;
		float RpRprEpsSQmQ=RpRprSQ-0.25;
		float RSQ=R*R;
		float RepsSQmQ=RSQ-0.25;
		int iRmax=int(float(RpRpr+1.5));
		
		int startX=ix-iRmax;
		int startY=iy-iRmax;
		int startZ=iz-iRmax;
		int endX=ix+iRmax;
		int endY=iy+iRmax;
		int endZ=iz+iRmax;
		if(startX<0)startX=0;
		if(startY<0)startY=0;
		if(startZ<0)startZ=0;
		if(endX>=GS_X)endX=GS_X-1;
		if(endY>=GS_Y)endY=GS_Y-1;
		if(endZ>=GS_Z)endZ=GS_Z-1;
		
		for(ix2=startX;ix2<=endX;ix2++)
			for(iy2=startY;iy2<=endY;iy2++)
				for(iz2=startZ;iz2<=endZ;iz2++)
		{
			fx2=(float)ix2;
			fy2=(float)iy2;
			fz2=(float)iz2;
			float drx=fx2-rcx;
			float dry=fy2-rcy;
			float drz=fz2-rcz;
			float drSQ=drx*drx+dry*dry+drz*drz;
			float drSQx=drSQ+drx;
			float drSQy=drSQ+dry;
			float drSQz=drSQ+drz;
			GrdPnt2=ix2+iy2*GS_X+iz2*GS_XY;
			
			if(drSQx<RpRprEpsSQmQ)
			{
				iEpsilon[0][GrdPnt2]=iatom+ID_Start;
			}
			if(drSQy<RpRprEpsSQmQ)
			{
				iEpsilon[1][GrdPnt2]=iatom+ID_Start;
			}
			if(drSQz<RpRprEpsSQmQ)
			{
				iEpsilon[2][GrdPnt2]=iatom+ID_Start;
			}
		}
	}
	
	return EXIT_SUCCESS;
}
int BuildWorldCmp::FinalazeSEV(ContWorld* world,int **iEpsilon)
{
	DbgPrint0("BuildWorldEu::FinalazeSEV\n");
	int GrdPnt,GrdPnt2;
	int i,j,k;
	int ix,iy,iz;
	int ix2,iy2,iz2;
	int GS_X=world->GS_X;
	int GS_Y=world->GS_Y;
	int GS_Z=world->GS_Z;
	int GS_XY=world->GS_XY;
	int GS_XYZ=world->GS_XYZ;
	
	
	float Rpr=Rwat;
	float RprSQ=Rpr*Rpr;
	float RepsSQmQ=RprSQ-0.25;
	int iRmax=int(float(Rwat+1.5));
	
	
	int GridIDs[6];
	float rcx,rcy,rcz;
	float fx,fy,fz;
	float fx2,fy2,fz2;
	float d;
	int count0=0;
	int iatom;
	for(ix=1;ix<GS_X-1;ix++)
		for(iy=1;iy<GS_Y-1;iy++)
			for(iz=1;iz<GS_Z-1;iz++)
	{
		GrdPnt=ix+iy*GS_X+iz*GS_XY;
		fx=(float)ix;
		fy=(float)iy;
		fz=(float)iz;
		GridIDs[0]=abs(iEpsilon[0][GrdPnt]);
		GridIDs[1]=abs(iEpsilon[1][GrdPnt]);
		GridIDs[2]=abs(iEpsilon[2][GrdPnt]);
		GridIDs[3]=abs(iEpsilon[0][GrdPnt-1]);
		GridIDs[4]=abs(iEpsilon[1][GrdPnt-GS_X]);
		GridIDs[5]=abs(iEpsilon[2][GrdPnt-GS_XY]);
		
		int sum=DielIdx[GridIDs[0]] + DielIdx[GridIDs[1]] + DielIdx[GridIDs[2]] +  DielIdx[GridIDs[3]] + DielIdx[GridIDs[4]] + DielIdx[GridIDs[5]];
		if(sum!=6*DielIdx[GridIDs[0]])
		{
			count0++;
			int GridID=GridIDs[0];
			float dmax=1000000.0;
			for(i=0;i<6;i++)
			{
				if(GridIDs[i]!=BulkParam->ID_GridBegin)
				{
					d=GGOIdx[GridIDs[i]]->GetDistSQToUnitCenter(fx,fy,fz,0.0,Rwat,GridIDs[i]);
					if(d<dmax)
					{
						GridID=GridIDs[i];
						dmax=d;
					}
				}
			}
			GGOIdx[GridID]->GetClosestSurfPoint(&fx,&fy,&fz,0.0,Rwat,GridID);
			
			int startX=ix-iRmax;
			int startY=iy-iRmax;
			int startZ=iz-iRmax;
			int endX=ix+iRmax;
			int endY=iy+iRmax;
			int endZ=iz+iRmax;
			if(startX<0)startX=0;
			if(startY<0)startY=0;
			if(startZ<0)startZ=0;
			if(endX>=GS_X)endX=GS_X-1;
			if(endY>=GS_Y)endY=GS_Y-1;
			if(endZ>=GS_Z)endZ=GS_Z-1;
			
			for(ix2=startX;ix2<=endX;ix2++)
				for(iy2=startY;iy2<=endY;iy2++)
					for(iz2=startZ;iz2<=endZ;iz2++)
			{
				fx2=(float)ix2;
				fy2=(float)iy2;
				fz2=(float)iz2;
				float drx=fx2-fx;
				float dry=fy2-fy;
				float drz=fz2-fz;
				float drSQ=drx*drx+dry*dry+drz*drz;
				GrdPnt2=ix2+iy2*GS_X+iz2*GS_XY;
				if(drSQ+drx<RepsSQmQ)
					iEpsilon[0][GrdPnt2]=-abs(iEpsilon[0][GrdPnt2]);
				if(drSQ+dry<RepsSQmQ)
					iEpsilon[1][GrdPnt2]=-abs(iEpsilon[1][GrdPnt2]);
				if(drSQ+drz<RepsSQmQ)
					iEpsilon[2][GrdPnt2]=-abs(iEpsilon[2][GrdPnt2]);
			}
		}
	}
	DbgPrint0("BigBoarders %d\n",count0);
	for(i=0;i<3;i++)
		for(GrdPnt=0;GrdPnt<GS_XYZ;GrdPnt++)
	{
		if(iEpsilon[i][GrdPnt]<0)iEpsilon[i][GrdPnt]=BulkParam->ID_GridBegin;
		iEpsilon[i][GrdPnt]=DielIdx[iEpsilon[i][GrdPnt]];
	}
	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
BuildWorldScaled::BuildWorldScaled()
{
	InitZero();
}
BuildWorldScaled::~BuildWorldScaled()
{
	Clear();
}
int BuildWorldScaled::InitZero()
{
	PnpsObject::SetName("BuildWorldScaled");
	
	BoundaryCondition=0;
	Rwat=1.4f;
	BulkParam=NULL;
	DielNum=0;
	DielConst=NULL;
	return EXIT_SUCCESS;
}
int BuildWorldScaled::Clear()
{
	return EXIT_SUCCESS;
}
int BuildWorldScaled::BuildContWorld(ContWorld* world)
{
	DbgPrint1("BuildWorldScaled::BuildContWorld()\n");
	int i,j,k;
	float Displacement[3]={0.0,0.0,0.0};
	//float* Surf=new float[world->GridSizeXYZGlobal*4];
	//int* iVtmp=new int[world->GridSizeXYZGlobal];
	float fpoh= 4*M_PI*world->GridScale;
	float coef=fpoh*COANGS/(world->GridScale*world->GridScale*world->GridScale);
	int iElt;
	
	//Change Parameters for internal 
// 	this->ChangeUnits(world,true);
// 	for(iElt=0;iElt<GOElms.size();iElt++)
// 	{
// 		GOElms[iElt]->ChangeUnits(world,true);
// 	}
	world->NIndexing=new NodeIndexing();
	NodeIndexing *NIndexingNew=world->NIndexing;
	
	unsigned int GS[3]={world->GridSizeGlobal[0],world->GridSizeGlobal[1],world->GridSizeGlobal[2]};
	NIndexingNew->SetNNodes(GS,world->GridScale);
	//Set Values Table
	
	for(i=0;i<DielNum;i++)NIndexingNew->Eps[i]=DielConst[i];
	
	BuildContWorldDielMaps(world,NIndexingNew);
	BuildContWorldCharges(world,NIndexingNew);
	
	return EXIT_SUCCESS;
}
int BuildWorldScaled::BuildContWorldCharges(ContWorld* world,NodeIndexing* NIndexingNew)
{
	int i,iElt;
	float* Qstat = new float[world->GridSizeXYZGlobal];
	for(i=0;i<world->GridSizeXYZGlobal;i++)Qstat[i]=0.0;
	
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		GOElms[iElt]->ChangeUnits(world,true);
		GOElms[iElt]->BuildPreMapsCharges(NULL,world,Qstat);
		GOElms[iElt]->ChangeUnits(world,false);
	}
	double Q=0.0;
	for(i=0;i<world->GridSizeXYZGlobal;i++)Q+=Qstat[i];
	NIndexingNew->SetChargeMapFromArray(Qstat);
	delete [] Qstat;
	DbgPrint1("Charge of Map: %.6g [e]\n",Q/(4.0*M_PI*world->GridScale));
	return EXIT_SUCCESS;
}
/* Print n as a binary number */
void printbitssimple(int n) {
	int i;
	for(i = 31; i>=0; i--) {
		if((1<<i)&n) printf("1");
		else printf("0");
		if((i)%8==0)printf("|");
	}
}

int BuildWorldScaled::BuildContWorldDielMaps(ContWorld* world,NodeIndexing* NIndexingNew)
{
	NodeIndex* NIndex=NIndexingNew->NIndex;
	NodeIndex ClearMask=~(NodeIndexing::Epsilon0|NodeIndexing::Epsilon1|NodeIndexing::Epsilon2);
	printf("NodeIndexing::Epsilon0=");printbitssimple(NodeIndexing::Epsilon0);printf("\n");
	printf("NodeIndexing::Epsilon1=");printbitssimple(NodeIndexing::Epsilon1);printf("\n");
	printf("NodeIndexing::Epsilon2=");printbitssimple(NodeIndexing::Epsilon2);printf("\n");
	printf("ClearMask             =");printbitssimple(ClearMask);printf("\n");
	NodeIndex EpsBulk=0;
	int Val=BulkParam->Epsilon;
	EpsBulk+=Val<<NodeIndexing::Epsilon0Sft;
	EpsBulk+=Val<<NodeIndexing::Epsilon1Sft;
	EpsBulk+=Val<<NodeIndexing::Epsilon2Sft;
	printf("EpsBulk               =");printbitssimple(EpsBulk);printf("\n");
	
	int j;
	for(j=0;j<world->GridSizeXYZGlobal;j++)
		NIndex[j]=(NIndex[j]&ClearMask)+EpsBulk;
	
	
	/*int i,iElt;
	float* Qstat = new float[world->GridSizeXYZGlobal];
	for(i=0;i<world->GridSizeXYZGlobal;i++)Qstat[i]=0.0;
	
	for(iElt=0;iElt<GOElms.size();iElt++)
	{
		GOElms[iElt]->ChangeUnits(world,true);
		GOElms[iElt]->BuildPreMapsCharges(NULL,world,Qstat);
		GOElms[iElt]->ChangeUnits(world,false);
	}
	double Q=0.0;
	for(i=0;i<world->GridSizeXYZGlobal;i++)Q+=Qstat[i];
	NIndexingNew->SetChargeMapFromArray(Qstat);
	delete [] Qstat;
	DbgPrint1("Charge of Map: %.6g [e]\n",Q/(4.0*M_PI*world->GridScale));*/
	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
DielDiffMapsPatcher::DielDiffMapsPatcher(ContWorld* world)
{
	m_ContWorld=world;
	eps[0]=NULL;
	eps[1]=NULL;
	eps[2]=NULL;

	refDiff[0]=NULL;
	refDiff[1]=NULL;
	refDiff[2]=NULL;
	refDiff[3]=NULL;
}
DielDiffMapsPatcher::~DielDiffMapsPatcher()
{
	DeleteCArray(eps[0]);
	DeleteCArray(eps[1]);
	DeleteCArray(eps[2]);

	DeleteCArray(refDiff[0]);
	DeleteCArray(refDiff[1]);
	DeleteCArray(refDiff[2]);
	DeleteCArray(refDiff[3]);
}
void DielDiffMapsPatcher::GetIntDielMaps()
{
	DeleteCArray(eps[0]);
	DeleteCArray(eps[1]);
	DeleteCArray(eps[2]);

	NodeIndexing *NI=m_ContWorld->NIndexing;
	
	eps[0]=NI->GetIntArrayFromIndexField(NI->DielConst,NI->Epsilon0);
	eps[1]=NI->GetIntArrayFromIndexField(NI->DielConst,NI->Epsilon1);
	eps[2]=NI->GetIntArrayFromIndexField(NI->DielConst,NI->Epsilon2);
}
void DielDiffMapsPatcher::PushIntDielMaps()
{
	if(m_ContWorld==NULL || eps[0]==NULL || eps[1]==NULL || eps[2]==NULL)
	{
		pnpError("Cann't push some values are not initiated\n");
		return;
	}
	NodeIndexing *NI=m_ContWorld->NIndexing;
	
	NI->SetIndexFieldFromIntArray(eps[0],NodeIndexing::Epsilon0,NodeIndexing::Epsilon0Sft);
	NI->SetIndexFieldFromIntArray(eps[1],NodeIndexing::Epsilon1,NodeIndexing::Epsilon1Sft);
	NI->SetIndexFieldFromIntArray(eps[2],NodeIndexing::Epsilon2,NodeIndexing::Epsilon2Sft);

	DeleteCArray(eps[0]);
	DeleteCArray(eps[1]);
	DeleteCArray(eps[2]);
}
void DielDiffMapsPatcher::PatchDielMaps(int epsToOver,int epsNew,float x, float y, float z1, float z2,float m_R)
{
	if(m_ContWorld==NULL || eps[0]==NULL || eps[1]==NULL || eps[2]==NULL)
	{
		pnpError("Cann't patch some values are not initiated\n");
		return;
	}
	NodeIndexing *NI=m_ContWorld->NIndexing;
	pnpPrint("===============================================================================\n");
    pnpPrint("Processing Element : Diel.Const %.3f -> %.3f\n",EPKT*NI->Eps[epsToOver],EPKT*NI->Eps[epsNew]);
    pnpPrint("\tx,y = %.3f\t%.3f\tz1,z2 = %.3f\t%.3f\t R = %.3f\n",x, y, z1, z2, m_R);

	int i,j,k;
	int GS_X=m_ContWorld->GridSize[0];
	int GS_XY=m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1];
	int GS_XYZ=m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1]*m_ContWorld->GridSize[2];

	//convert units
	
	int irc[3];
	float frc[3];
	for(j=0;j<3;j++)
	{
		irc[j]=m_ContWorld->GridSizeGlobal[j]/2;
		frc[j]=irc[j];
	}
	//ToInternal
	x=x*m_ContWorld->GridScale+frc[0];
	y=y*m_ContWorld->GridScale+frc[1];
	z1=z1*m_ContWorld->GridScale+frc[2];
	z2=z2*m_ContWorld->GridScale+frc[2];
	m_R=m_R*m_ContWorld->GridScale;

	for(i=0;i<3;i++)
	{
		int ix,iy,iz,GrdPnt;
		float X,Y,Z1,Z2;
		int iR;
		float R,RSQ,RtmpSQ;
		float Displ[3]={0.0,0.0,0.0};
		int start[3];
		int end[3];
		Displ[i]=-0.5;

		iR=(int)(m_R+0.5);
		RSQ=m_R*m_R;

		X=x+Displ[0];
		Y=y+Displ[1];
		Z1=z1+Displ[2];
		Z2=z2+Displ[2];

		start[0]=(int)(X+0.5-iR);
		end[0]=(int)(X+0.5+iR);
		start[1]=(int)(Y+0.5-iR);
		end[1]=(int)(Y+0.5+iR);
		start[2]=(int)(Z1-1.5);
		end[2]=(int)(Z2+1.5);

		for(k=0;k<3;k++)
		{
			if(start[k]<0)start[k]=0;
			if(end[k]>m_ContWorld->GridSize[k]-1)end[k]=m_ContWorld->GridSize[k]-1;
		}
		
		for(ix=start[0];ix<=end[0];ix++)
			for(iy=start[1];iy<=end[1];iy++)
				for(iz=start[2];iz<=end[2];iz++)
		{
			GrdPnt=ix+iy*GS_X+iz*GS_XY;
			//eps[i][GrdPnt]=epsNew;
			//printf("eps[%d]=%d\n",i,GrdPnt);
			if(eps[i][GrdPnt]==epsToOver)
			{
				RtmpSQ=(x-ix)*(x-ix)+(y-iy)*(y-iy);
				if(RtmpSQ<=RSQ)//RtmpSQ<=RSQ
				{
					eps[i][GrdPnt]=epsNew;
				}
			}
		}

	}
}
void DielDiffMapsPatcher::InitNewDiff()
{
	DeleteCArray(refDiff[0]);
	DeleteCArray(refDiff[1]);
	DeleteCArray(refDiff[2]);
	DeleteCArray(refDiff[3]);

	NodeIndexing *NI=m_ContWorld->NIndexing;
	int i,GrdPnr;
	int initDiffFromIAV=0;
	if(m_ContWorld->D==NULL)initDiffFromIAV=1;
	else
	{
		if(m_ContWorld->D[0]==NULL)initDiffFromIAV=1;
		if(m_ContWorld->D[1]==NULL)initDiffFromIAV=1;
	}

	if(initDiffFromIAV)
	{
		m_ContWorld->D=new float*[m_ContWorld->NIonsTypes];
		for(i=0;i<m_ContWorld->NIonsTypes;i++)
		{
			refDiff[i]=NI->GetCMap(NodeIndexing::DiffConst, NI->IonField[i]);

			m_ContWorld->D[i]=new float[m_ContWorld->GridSizeXYZGlobal];
			for(GrdPnr=0;GrdPnr<m_ContWorld->GridSizeXYZGlobal;GrdPnr++)
				m_ContWorld->D[i][GrdPnr]=refDiff[i][GrdPnr];
		}
	}
	else
	{
		for(i=0;i<m_ContWorld->NIonsTypes;i++)
		{
			refDiff[i]=m_ContWorld->D[i];
			m_ContWorld->D[i]=new float[m_ContWorld->GridSizeXYZGlobal];
			for(GrdPnr=0;GrdPnr<m_ContWorld->GridSizeXYZGlobal;GrdPnr++)
				m_ContWorld->D[i][GrdPnr]=refDiff[i][GrdPnr];
		}
	}
}

void DielDiffMapsPatcher::DelRefDiff()
{
	DeleteCArray(refDiff[0]);
	DeleteCArray(refDiff[1]);
	DeleteCArray(refDiff[2]);
	DeleteCArray(refDiff[3]);
}
void DielDiffMapsPatcher::PatchDiffMaps(float Dscale1K, float Dscale2K, float Dscale1Cl, float Dscale2Cl,float x, float y, float z1, float z2,float m_R)
{
	if(m_ContWorld==NULL || refDiff[0]==NULL || refDiff[1]==NULL)
	{
		pnpError("Cann't patch some values are not initiated\n");
		return;
	}
	NodeIndexing *NI=m_ContWorld->NIndexing;
	pnpPrint("===============================================================================\n");
    pnpPrint("Processing Element :\n");
	pnpPrint("\tScaling Diffusion: %.3f -> %.3f\n",Dscale1K,Dscale2K);
	pnpPrint("\tScaling Diffusion: %.3f -> %.3f\n",Dscale1Cl,Dscale2Cl);
    pnpPrint("\tScaling at: x,y = %.3f\t%.3f\tz1,z2 = %.3f\t%.3f\t R = %.3f\n",x, y, z1, z2, m_R);

	int i,j,k;
	int GS_X=m_ContWorld->GridSize[0];
	int GS_XY=m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1];
	int GS_XYZ=m_ContWorld->GridSize[0]*m_ContWorld->GridSize[1]*m_ContWorld->GridSize[2];

	//convert units
	
	int irc[3];
	float frc[3];
	for(j=0;j<3;j++)
	{
		irc[j]=m_ContWorld->GridSizeGlobal[j]/2;
		frc[j]=irc[j];
	}
	//ToInternal
	x=x*m_ContWorld->GridScale+frc[0];
	y=y*m_ContWorld->GridScale+frc[1];
	z1=z1*m_ContWorld->GridScale+frc[2];
	z2=z2*m_ContWorld->GridScale+frc[2];
	m_R=m_R*m_ContWorld->GridScale;


	int ix,iy,iz,GrdPnt;
	int iR;
	float R,RSQ,RtmpSQ;
	float DscaleK, DscaleCl;
	int start[3];
	int end[3];

	iR=(int)(m_R+0.5);
	RSQ=m_R*m_R;

	start[0]=(int)(x+0.5-iR);
	end[0]=(int)(x+0.5+iR);
	start[1]=(int)(y+0.5-iR);
	end[1]=(int)(y+0.5+iR);
	start[2]=(int)(z1-1.5);
	end[2]=(int)(z2+1.5);

	for(k=0;k<3;k++)
	{
		if(start[k]<0)start[k]=0;
		if(end[k]>m_ContWorld->GridSize[k]-1)end[k]=m_ContWorld->GridSize[k]-1;
	}
		
	for(ix=start[0];ix<=end[0];ix++)
	{
		for(iy=start[1];iy<=end[1];iy++)
		{
			
			RtmpSQ=(x-ix)*(x-ix)+(y-iy)*(y-iy);
			if(RtmpSQ<=RSQ)
			{
				for(iz=start[2];iz<=end[2];iz++)
				{
					if(iz>=z1 && iz <=z2)
					{
						GrdPnt=ix+iy*GS_X+iz*GS_XY;

						if(Dscale1K==Dscale2K)
							DscaleK=Dscale1K;
						else
						{
							DscaleK=Dscale1K+(iz-z1)*(Dscale2K-Dscale1K)/(z2-z1);
						}
						if(Dscale1K==Dscale2K)
							DscaleCl=Dscale1Cl;
						else
						{
							DscaleCl=Dscale1Cl+(iz-z1)*(Dscale2Cl-Dscale1Cl)/(z2-z1);
						}
						m_ContWorld->D[0][GrdPnt]=DscaleK*refDiff[0][GrdPnt];
						m_ContWorld->D[1][GrdPnt]=DscaleCl*refDiff[1][GrdPnt];
					}
				}
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////////