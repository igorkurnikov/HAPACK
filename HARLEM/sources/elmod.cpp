//
// C++ Implementation: elmod
//
// Description: 
//
//
// Author: mikola <mikola@linux>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#define HARLEM_MPI 1
#include <mpi.h>

#include <regex>

#include <boost/filesystem/path.hpp>
#include <boost/algorithm/string.hpp>


#include "hamolset.h"

#include "elmod.h"
#include "protonredox.h"
#include "pnpdebug.h"
#include "pnpsapp.h"
#include "contworld.h"
#include "buildworldni.h"
#include "poissonsolver.h"
#include "poissonboltzmannsolver.h"
#include "pbwithljsolver.h"
#include "mapio.h"
#include "pnputil.h"
#include "pnpconstants.h"
#include "haconsts.h"

#ifdef WITH_CUDA
extern "C" {
#include "cps.h"
#include "cps_prep.h"
#include "cps_util.h"
}
#endif


#if defined(_MSC_VER)
#include <process.h>
#endif

#include "command.h"
#include "stdlib.h"
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <float.h>
#include <math.h>
#include "electrostmod.h"
#include "hamolecule.h"

#include "etcoupl.h"
#include "hasurface.h"
#include "hamolview.h"
#include "harlemapp.h"
#include "randomip.h"
#include "tinyxml.h"
//#include "haconsts.h"

#include <sstream>

#include "mpi.h"
///////////////////////////////////////////////////////////////////////////////
HaAtomParmEntry::HaAtomParmEntry(const char* cName,double _dbl0,double _dbl1)
{
	Name=cName;
	dbl0=_dbl0;
	dbl1=_dbl1;
}
HaAtomParmEntry::~HaAtomParmEntry()
{
}
void HaAtomParmEntry::SetValues(double _dbl0,double _dbl1)
{
	dbl0=_dbl0;
	dbl1=_dbl1;
}
///////////////////////////////////////////////////////////////////////////////
HaResParmEntry::HaResParmEntry(const char* cName, const char* cMod)
{
	Name=cName;
	Mod=cMod;
}
HaResParmEntry::~HaResParmEntry()
{
}
HaAtomParmEntry* HaResParmEntry::SetAtom(const char* cName,double _dbl0,double _dbl1)
{
	HaAtomParmEntry* atm=GetAtom(cName);
	if(atm==NULL)
	{
		atm=new HaAtomParmEntry(cName, _dbl0, _dbl1);
		Atms.push_back(atm);
	}
	else
	{
		atm->SetValues(_dbl0, _dbl1);
	}
	return atm;
}
HaAtomParmEntry* HaResParmEntry::GetAtom(const char* cName)
{
	int i;
	std::string CompName=cName;
	for(i=0;i<Atms.size();i++)
	{
		if(Atms[i]->Name==CompName)
			return Atms[i];
	}
	return NULL;
}
int HaResParmEntry::GetAtomListNum(const char* cName)
{
	int i;
	std::string CompName=cName;
	for(i=0;i<Atms.size();i++)
	{
		if(Atms[i]->Name==CompName)
			return i;
	}
	return -1;
}
void HaResParmEntry::DelAtom(const char* cName)
{
	int i=GetAtomListNum(cName);
	if(i<0)
	{
		PrintLog("Error: There is no atom with name %s in the residue\n",cName);
		return;
	}
	delete Atms[i];
	Atms.erase(Atms.begin()+i);
}
void HaResParmEntry::CopyAtomsFromRes(HaResParmEntry* CpRes)
{
	if(CpRes==NULL)return;
	int i;
	for(i=0;i<CpRes->Atms.size();i++)
	{
		SetAtom(CpRes->Atms[i]->Name.c_str(),CpRes->Atms[i]->dbl0,CpRes->Atms[i]->dbl1);
	}
}
void HaResParmEntry::PrintInfo()
{
	PrintLog("  Name:%s\t Mod:%s\n",Name.c_str(),Mod.c_str());
	int i;
	for(i=0;i<Atms.size();i++)
	{
		PrintLog("    %6s%12.5f%12.5f\n",Atms[i]->Name.c_str(),Atms[i]->dbl0,Atms[i]->dbl1);
	}
}
void HaResParmEntry::PrintSumDbl0()
{
	int i;
	double dtmp=0.0;
	for(i=0;i<Atms.size();i++)
	{
		dtmp+=Atms[i]->dbl0;
	}
	PrintLog("  SumDbl0[%s#%s]=%15.5f\n",Name.c_str(),Mod.c_str(),dtmp);
}
///////////////////////////////////////////////////////////////////////////////
HaAtomsParmFF::HaAtomsParmFF(const char* cName, const char* cRef, const char* cNotes)
{
	Name=cName;
	Ref=cRef;
	Notes=cNotes;
}
HaAtomsParmFF::~HaAtomsParmFF()
{
	
}
HaResParmEntry* HaAtomsParmFF::NewRes(const char* cName, const char* cMod)
{
	HaResParmEntry* res=new HaResParmEntry(cName,cMod);
	Ress.push_back(res);
	return res;
}
HaResParmEntry* HaAtomsParmFF::GetRes(const char* cName, const char* cMod)
{
	int i;
	std::string CompName=cName;
	std::string CompMod=cMod;
	for(i=0;i<Ress.size();i++)
	{
		if(Ress[i]->Name==CompName&&Ress[i]->Mod==CompMod)
			return Ress[i];
	}
	return NULL;
}
HaResParmEntry* HaAtomsParmFF::SetRes(const char* cName, const char* cMod)
{
	HaResParmEntry* res=GetRes(cName,cMod);
	if(res==NULL)
	{
		res=NewRes(cName,cMod);
	}
	return res;
}
HaResParmEntry* HaAtomsParmFF::CopyRes(const char* cName, const char* cMod, const char* cNameFrom, const char* cModFrom)
{
	HaResParmEntry* res=SetRes(cName,cMod);
	HaResParmEntry* resFrom=GetRes(cNameFrom,cModFrom);
	if(resFrom==NULL)
	{
		PrintLog("Error: Can not find Residue %s#%s in FF:%s (to copy from)\n", cNameFrom, cModFrom, Name.c_str());
	}
	else
	{
		int i;
		for(i=0;i<resFrom->Atms.size();i++)
		{
			res->SetAtom(resFrom->Atms[i]->Name.c_str(),resFrom->Atms[i]->dbl0,resFrom->Atms[i]->dbl1);
		}
	}
	return res;
}
void HaAtomsParmFF::PrintFFInfo()
{
	PrintLog("Name:%s\n",Name.c_str());
	PrintLog("Reference:%s\n",Ref.c_str());
	PrintLog("Notes:%s\n",Notes.c_str());
}
void HaAtomsParmFF::PrintEntries()
{
	int i;
	for(i=0;i<Ress.size();i++)
	{
		Ress[i]->PrintInfo();
	}
	
}
void HaAtomsParmFF::PrintSumDbl0()
{
	int i;
	for(i=0;i<Ress.size();i++)
	{
		Ress[i]->PrintSumDbl0();
	}
}
///////////////////////////////////////////////////////////////////////////////
HaAtomsParmDB::HaAtomsParmDB(const char* cName)
{
	Name=cName;
}
HaAtomsParmDB::~HaAtomsParmDB()
{
	int i;
	for(i=0;i<FFs.size();i++)
	{
		delete FFs[i];
	}
	FFs.resize(0);
}
HaAtomsParmFF* HaAtomsParmDB::NewFF(const char* Name, const char* Ref, const char* Notes)
{
	HaAtomsParmFF* ff=new HaAtomsParmFF(Name,Ref,Notes);
	FFs.push_back(ff);
	return ff;
}
HaAtomsParmFF* HaAtomsParmDB::GetFF(const char* cName)
{
	int i;
	std::string CompName=cName;
	
	for(i=0;i<FFs.size();i++)
	{
		//PrintLog("%s =? %s\n",FFs[i]->Name.c_str(),CompName.c_str());
		if(FFs[i]->Name==CompName)
			return FFs[i];
	}
	return NULL;
}
int HaAtomsParmDB::NumFF()
{
	return FFs.size();
}
HaAtomsParmFF* HaAtomsParmDB::GetFFbyNum(int i)
{
	if(i<FFs.size())
	{
		return FFs[i];
	}
	else
	{
		return NULL;
	}
}
void HaAtomsParmDB::CompareAtomNamesInFFs(const char* cName,const char* cNameRef)
{
	HaAtomsParmFF* ff=GetFF(cName);
	HaAtomsParmFF* ffref=GetFF(cNameRef);
	PrintLog("Comparing %s with %s\n",cName,cNameRef);
	if(ff==NULL)
	{
		PrintLog("Error: Can not find FF with name %s\n",cName);
		return;
	}
	if(ffref==NULL)
	{
		PrintLog("Error: Can not find FF with name %s\n",cNameRef);
		return;
	}
	int ir;
	for(ir=0;ir<ff->Ress.size();ir++)
	{
		HaResParmEntry* res=ff->Ress[ir];
		HaResParmEntry* rref=ffref->GetRes(res->Name.c_str(),res->Mod.c_str());
		if(ffref==NULL)
		{
			PrintLog("Error: Can not find residue: %s:%s#%s\n",cNameRef,res->Name.c_str(),res->Mod.c_str());
		}
		else
		{
			int ia;
			if(res->Atms.size()!=rref->Atms.size())
			{
				PrintLog("Error: Number of atoms in res:%s#%s between two FF is different\n",res->Name.c_str(),res->Mod.c_str());
			}
			for(ia=0;ia<res->Atms.size();ia++)
			{
				HaAtomParmEntry* a=res->Atms[ia];
				HaAtomParmEntry* aref=rref->GetAtom(a->Name.c_str());
				if(aref==NULL)
					PrintLog("Error: There is no atom:%s in ref. FF:%s\n",a->Name.c_str(),cNameRef);
			}
		}
	}
}
void HaAtomsParmDB::PrintSideToSideFFs(const char* cName,const char* cNameRef)
{
	HaAtomsParmFF* ff=GetFF(cName);
	HaAtomsParmFF* ffref=GetFF(cNameRef);
	PrintLog("Printing side to side params for FF:%s and same res/atoms from FF:%s\n",cName,cNameRef);
	if(ff==NULL)
	{
		PrintLog("Error: Can not find FF with name %s\n",cName);
		return;
	}
	if(ffref==NULL)
	{
		PrintLog("Error: Can not find FF with name %s\n",cNameRef);
		return;
	}
	int ir;
	for(ir=0;ir<ff->Ress.size();ir++)
	{
		HaResParmEntry* res=ff->Ress[ir];
		HaResParmEntry* rref=ffref->GetRes(res->Name.c_str(),res->Mod.c_str());
		if(ffref==NULL)
		{
			PrintLog("Error: Can not find residue: %s:%s#%s\n",cNameRef,res->Name.c_str(),res->Mod.c_str());
		}
		else
		{
			int ia;
			PrintLog("  Name:%s\t Mod:%s\n",res->Name.c_str(),res->Mod.c_str());
			if(res->Atms.size()!=rref->Atms.size())
			{
				PrintLog("Error: Number of atoms in res:%s#%s between two FF is different\n",res->Name.c_str(),res->Mod.c_str());
			}
			for(ia=0;ia<res->Atms.size();ia++)
			{
				HaAtomParmEntry* a=res->Atms[ia];
				HaAtomParmEntry* aref=rref->GetAtom(a->Name.c_str());
				if(aref==NULL)
					PrintLog("Error: There is no atom:%s in ref. FF:%s\n",a->Name.c_str(),cNameRef);
				else
					PrintLog("    %6s |%12.5f%12.5f |%12.5f%12.5f\n",a->Name.c_str(),a->dbl0,a->dbl1,aref->dbl0,aref->dbl1);
			}
		}
	}
}
void HaAtomsParmDB::SetAtomsParam(MolSet* pmset, const char *FFName)
{
	HaAtomsParmFF* ff=GetFF(FFName);
	if(ff==NULL)
	{
		PrintLog("Error: Can not find FF with name %s\n",FFName);
		return;
	}
	PrintLog("\nGoing to use following FF:\n");
	ff->PrintFFInfo();
	
	HaAtom* aptr;
	HaResidue* rptr;
	ResidueIteratorMolSet ritr(pmset);
	for(rptr = ritr.GetFirstRes(); rptr; rptr = ritr.GetNextRes() )
	{
		if(rptr->HasSelectedAtoms())
		{
			HaResParmEntry* res=ff->GetRes(rptr->GetName(),rptr->GetNameModifier());
			if(res==NULL)
			{
				PrintLog("Error: Can not find Residue %s#%s in FF:%s\n", rptr->GetName(), rptr->GetNameModifier(), FFName);
			}
			else
			{
				AtomIteratorResidue aitr_res(rptr);
				for(aptr = aitr_res.GetFirstAtom(); aptr; aptr = aitr_res.GetNextAtom())
				{
					if(pmset->save_opt_default.save_selected && !aptr->Selected())
						continue;
					HaAtomParmEntry* atm=res->GetAtom(aptr->GetName());
					if(atm==NULL)
					{
						PrintLog("Error: Can not find atoms:%s in Residue %s#%s of FF:%s\n", aptr->GetName(), rptr->GetName(), rptr->GetNameModifier(), FFName);
					}
					else
					{
						aptr->SetCharge(atm->dbl0);
						aptr->radius=atm->dbl1;
					}
					
				}
			}
		}
	}
	PrintLog("\n");
}
void HaAtomsParmDB::PrintFFsInfo()
{
	int i;
	PrintLog("DB:%s\n",Name.c_str());
	for(i=0;i<FFs.size();i++)
	{
		PrintLog("\n");
		FFs[i]->PrintFFInfo();
	}
	PrintLog("\n");
}
///////////////////////////////////////////////////////////////////////////////
HaAtomsParmDB* QRDB=NULL;
HaAtomsParmDB* GetQRDB()
{
	if(QRDB==NULL)
	{
		PrintLog("QRDB is not initiated will start new DB\n");
		QRDB=new HaAtomsParmDB("QRDB");
	}
	//PrintLog("QRDB %p\n",QRDB);
	return QRDB;
}
void DeleteQRDB()
{
	if(QRDB!=NULL)
	{
		delete QRDB;
		QRDB=NULL;
	}
}
#ifdef ELMOD_COMPILE
ElMod::ElMod(MolSet* new_phost_mset)
	:ElectrostMod(new_phost_mset,COMP_MOD_EL)
{
	m_ContWorld=NULL;
	MoleculePosition=MoveToCenter;
	MembraneType=MemT_None;
	epsmem=2.0;
	
	BoxSetMethod=SetBoxPerFil;
	GridScale=0.0;
	OffsetX=0.0;
	OffsetY=0.0;
	OffsetZ=0.0;
	
	memZ0=-10.0;
	memZ1=10.0;
	memR0=10.0;
	memR1=200.0;
	memX=0.0;
	memY=0.0;
	
	BldType=BldNI;
	
	PNPSApp::InitPNPSApp();
	RoughWorld=NULL;
	Relaxation=-1.0;
	Convergence=0.0;
	
	RemovingCavitiesOnDielectricMap=false;
	RemCavOnDielFillWith=2;
	RemCavOnDielWhere2Look=1;
	
	PQRFile="";
	PREFile="";
	bPBSR=false;
	bConvIndConc=false;
	IonsHalfSigma.newsize(2);
	IonsHalfSigma.SetVal_idx0(0,(float)2.658);//K+
	IonsHalfSigma.SetVal_idx0(1,(float)1.868);//Cl-
	IonsFourEpsilon.newsize(2);
	IonsFourEpsilon.SetVal_idx0(0,(float)0.0005536);//K+
	IonsFourEpsilon.SetVal_idx0(1,(float)0.00467522);//Cl-
	
	t_x=NULL;
	t_y=NULL;
	t_z=NULL;
	
	GPUAccelerated=false;
	
	ChargeDist=GOAtoms::Linear;
	ChargeDistPar=1;
	MemoryLimitOnOneArray=-1;
	
	bInsertProbe=false;
	probe_x=0.0;
	probe_y=0.0;
	probe_z=0.0;
	probe_q=0.0;
	probe_R=0.0;
	probe_eps=2.0;
	Phi=0.0;
	Theta=0.0;
	bCalMolSet=true;
}


ElMod::~ElMod()
{
}

void ElMod::clear()
{
	DeleteCoorTmpArr();
	ElectrostMod::clear();
}
ElMod* ElMod::GetRoughWorld()
{
	if(RoughWorld==NULL)
	{
		RoughWorld=new ElMod(phost_mset);
	}
	return RoughWorld;
}
void ElMod::InsertProbe(double m_x,double m_y,double m_z,double m_q,double m_R,double m_eps)
{
	bInsertProbe=true;
	probe_x=m_x;
	probe_y=m_y;
	probe_z=m_z;
	probe_q=m_q;
	probe_R=m_R;
	probe_eps=m_eps;
}
void ElMod::UnsetProbe()
{
	bInsertProbe=false;
}
void ElMod::SetTranslateAndRotate(double m_x,double m_y,double m_z,double m_Phi,double m_Theta)
{
	OffsetX=m_x;
	OffsetY=m_y;
	OffsetZ=m_z;
	Phi=m_Phi;
	Theta=m_Theta;
}
bool ElMod::run(RunMode rmode)
{
	printf("ElMod::run\n");
	DefClock0;
	StartClock0;
	PNP_EXIT_FAIL_NULL(phost_mset,"ElMod does not assign to any MolSet\n");
	
	if(InitContWorld()==EXIT_FAILURE)
	{
		printf("ElMod::Error Can not initialize ContWorld\n");
		return false;
	}
	
	if(boundary==2)//Focusing
	{
		PNP_EXIT_FAIL_NULL(RoughWorld,"System for focusing is not set\n");
		RoughWorld->run(rmode);
		ContWorld* m_ContWorldRough=RoughWorld->m_ContWorld;
		float *tmp[1],*tmp2[1];
		tmp[0]=m_ContWorldRough->Potential;
		tmp2[0]=m_ContWorld->Potential;
		
		VectorField3D * Pot = new VectorField3D(m_ContWorldRough->GridSize, m_ContWorldRough->GridScale, 1, tmp);
		VectorField3D * PotFine = new VectorField3D(m_ContWorld->GridSize, m_ContWorld->GridScale, 1,tmp2);
		PotFine->InterpolateFromExt(Pot);
		delete PotFine;
		delete Pot;
		//m_ContWorldRough->Potential=NULL;
		//delete m_ContWorldRough;
		//RoughWorld->m_ContWorld=NULL;
		
		
		
	//PotFine->FillValue(0.0);
		
		//Pot->amode=VectorField3D::INTERNAL_ALLOC;
		
		//delete [] tmp[0];
		
	}
	
	if(BuildWorld()==EXIT_FAILURE)
	{
		printf("ElMod::Error Can not Build ContWorld\n");
		return false;
	}
	
	if(rionst==0.0)
	{
		if(SolvePoisson()==EXIT_FAILURE)
		{
			printf("ElMod::Error Can not solve Poisson Equation\n");
			return false;
		}
	}
	else
	{
		if(bPBSR)
		{
			if(SolvePoissonBoltzmannSoftRepultion()==EXIT_FAILURE)
			{
				printf("ElMod::Error Can not solve Poisson-Boltzmann Equation\n");
				return false;
			}
		}
		else
		{
			if(SolvePoissonBoltzmann()==EXIT_FAILURE)
			{
				printf("ElMod::Error Can not solve Poisson-Boltzmann Equation\n");
				return false;
			}
		}
	}
	StopClockWMes0("running electrostatic calculation");
	//return ElectrostMod::run(rmode);
	return true;
}
#ifdef WITH_CUDA
bool ElMod::runCPS(RunMode rmode)
{
	#if defined(WITH_CPS)
	if(GridScale==0.0)
	{
		pnpError("Cannot calculate GridScale yet\n");
	}
	
	DefClock0;
	StartClock0;
	//CPS_init();
	//Prepere atoms
	CPS_Atoms* a=new CPS_Atoms;
	a->Nq=0;
	a->Nsp=0;
	a->iEps=2;
	HaAtom* aptr;
	AtomIteratorMolSet aitr(phost_mset);
	
	//count atoms with Q and R
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if(phost_mset->m_save_selected && !aptr->Selected())
			continue;
		if(aptr->radius>0.0)
			a->Nsp++;
		if(aptr->GetCharge()!=0.0)
			a->Nq++;
	}
	a->R=new float4[a->Nsp];
	a->Q=new float4[a->Nq];
	int countR=0;
	int countQ=0;
	if(MoleculePosition==MoveToCenter)
	{
		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if(phost_mset->m_save_selected && !aptr->Selected())
				continue;
			if(aptr->radius>0.0)
			{
				a->R[countR].x=float(aptr->GetX()*BOHR_TO_ANG-gr_cent[0]);
				a->R[countR].y=float(aptr->GetY()*BOHR_TO_ANG-gr_cent[1]);
				a->R[countR].z=float(aptr->GetZ()*BOHR_TO_ANG-gr_cent[2]);
				a->R[countR].w=float(aptr->radius * BOHR_TO_ANG);
				countR++;
			}
			if(aptr->GetCharge()!=0.0)
			{
				a->Q[countQ].x=float(aptr->GetX()*BOHR_TO_ANG-gr_cent[0]);
				a->Q[countQ].y=float(aptr->GetY()*BOHR_TO_ANG-gr_cent[1]);
				a->Q[countQ].z=float(aptr->GetZ()*BOHR_TO_ANG-gr_cent[2]);
				a->Q[countQ].w=float(aptr->GetCharge());
				countQ++;
			}
		}
	}
	else if(MoleculePosition==DontMove)
	{
		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if(phost_mset->m_save_selected && !aptr->Selected())
				continue;
			if(aptr->radius>0.0)
			{
				a->R[countR].x=float(aptr->GetX()*BOHR_TO_ANG);
				a->R[countR].y=float(aptr->GetY()*BOHR_TO_ANG);
				a->R[countR].z=float(aptr->GetZ()*BOHR_TO_ANG);
				a->R[countR].w=float(aptr->radius * BOHR_TO_ANG);
				countR++;
			}
			if(aptr->GetCharge()!=0.0)
			{
				a->Q[countQ].x=float(aptr->GetX()*BOHR_TO_ANG);
				a->Q[countQ].y=float(aptr->GetY()*BOHR_TO_ANG);
				a->Q[countQ].z=float(aptr->GetZ()*BOHR_TO_ANG);
				a->Q[countQ].w=float(aptr->GetCharge());
				countQ++;
			}
		}
	}
	else if(MoleculePosition==UseOffset)
	{
		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if(phost_mset->m_save_selected && !aptr->Selected())
				continue;
			if(aptr->radius>0.0)
			{
				a->R[countR].x=float(aptr->GetX()*BOHR_TO_ANG-OffsetX);
				a->R[countR].y=float(aptr->GetY()*BOHR_TO_ANG-OffsetY);
				a->R[countR].z=float(aptr->GetZ()*BOHR_TO_ANG-OffsetZ);
				a->R[countR].w=float(aptr->radius * BOHR_TO_ANG);
				countR++;
			}
			if(aptr->GetCharge()!=0.0)
			{
				a->Q[countQ].x=float(aptr->GetX()*BOHR_TO_ANG-OffsetX);
				a->Q[countQ].y=float(aptr->GetY()*BOHR_TO_ANG-OffsetY);
				a->Q[countQ].z=float(aptr->GetZ()*BOHR_TO_ANG-OffsetZ);
				a->Q[countQ].w=float(aptr->GetCharge());
				countQ++;
			}
		}
	}
	else if(MoleculePosition==TranslateAndRotate)
	{
		Vec3D vPhi(0.0,0.0,1.0);
		Vec3D vTheta(0.0,1.0,0.0);
		double cosPhi=cos(Phi);
		double sinPhi=sin(Phi);
		double cosTheta=cos(-Theta);
		double sinTheta=sin(-Theta);
			
		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if(phost_mset->m_save_selected && !aptr->Selected())
				continue;
			Vec3D m_r((aptr->GetX()+OffsetX)*BOHR_TO_ANG, (aptr->GetY()+OffsetY)*BOHR_TO_ANG, (aptr->GetZ()+OffsetZ)*BOHR_TO_ANG);
				
			m_r.Rotate(vTheta, cosTheta, sinTheta);
			m_r.Rotate(vPhi, cosPhi, sinPhi);
			
			
			if(aptr->radius>0.0)
			{
				a->R[countR].x=float(m_r.GetX());
				a->R[countR].y=float(m_r.GetY());
				a->R[countR].z=float(m_r.GetZ());
				a->R[countR].w=float(aptr->radius * BOHR_TO_ANG);
				countR++;
			}
			if(aptr->GetCharge()!=0.0)
			{
				a->Q[countQ].x=float(m_r.GetX());
				a->Q[countQ].y=float(m_r.GetY());
				a->Q[countQ].z=float(m_r.GetZ());
				a->Q[countQ].w=float(aptr->GetCharge());
				countQ++;
			}
		}
	}
	int i;
	if(fabs(epsi-epsout)<0.001)
		a->iEps=1;
	else
		a->iEps=2;
	a->Walloc=0;
	StopClockWMes0("atoms prep");
	StartClock0;
	//Run Solver
	//Set CPSworld
	CPS_Solve_Input inp;
	inp.GS_X=nx;
	inp.GS_Y=ny;
	inp.GS_Z=nz;
	inp.GridScale=GridScale;
	inp.iEpsBulk=1;
	
	for(i=0;i<CPS_MAX_EPS_NUM;i++)
		inp.Eps[i]=0.0;
	inp.Eps[1]=epsout;
	inp.Eps[2]=epsi;
	inp.a=a;
	//Bondaries
	if(boundary==0)
		inp.BC=0;
	else if(boundary==3)
	{
		inp.BC=1;
		inp.BC_X=1;
		inp.BC_Y=1;
		inp.BC_Z=1;
	}
	else
		ErrorInMod("ElMod::runCPS", "Can not handle this type of bondaries\n");
	//
	inp.Rprobe= radprb;// * BOHR_TO_ANG;
	//solver params
	inp.MaxIter=nlit;
	inp.Rel=1.9;
	inp.Tol=0.0;
	//What PU to use
	if(GPUAccelerated)
		inp.UsePU=1;
	else
		inp.UsePU=0;
	//CPS_Solve
	CPS_Solver_Results* res=CPS_Solve(&inp);
	StopClockWMes0("CPS_Solve");
	
	//save potential
	float **fPot=new float*[1];
	int GS[3]={nx,ny,nz};
	fPot[0]=CPS_World_GetRegPot_Host(res->w);
	VectorField3D vPot(GS, inp.GridScale, 1,fPot);
	vPot.amode=VectorField3D::INTERNAL_ALLOC;
	vPot.WriteToFile("Pot3.bin");
	VectorField3D vPotRef("Pot2.bin");
	pnpPrint("RMSD for everything\n");
	vPot.RMSD(&vPotRef);
	pnpPrint("RMSD for internal\n");
	vPot.RMSDInternal(&vPotRef);
	//delete stuf
	CPS_Solver_Results_Free(res);
	CPS_Atoms_Free_Host(a);
	#endif
	return true;
}
#endif
int ElMod::InitContWorld()
{
	bool PBCX=(bool)iper[0];
	bool PBCY=(bool)iper[1];
	bool PBCZ=(bool)iper[2];
	
	//
	int na = 0;
	HaAtom* aptr;
	AtomIteratorMolSet aitr(phost_mset);

	bool calc_range = true;

	if( (this->max_coord.GetX() - this->min_coord.GetX()) > 0.001 &&
				(this->max_coord.GetY() - this->min_coord.GetY()) > 0.001 &&
				(this->max_coord.GetZ() - this->min_coord.GetZ()) > 0.001 ) 
	{
		calc_range = false;
	}

	double xmin1, ymin1,zmin1,xmax1,ymax1,zmax1;

	if(calc_range)
	{
		xmin1 = 10000.0;
		ymin1 = 10000.0;
		zmin1 = 10000.0;
		xmax1 = -10000.0;
		ymax1 = -10000.0;
		zmax1 = -10000.0;
	}
	else
	{
		xmin1 = this->min_coord.GetX_Ang();
		ymin1 = this->min_coord.GetY_Ang();
		zmin1 = this->min_coord.GetZ_Ang();

		xmax1 = this->max_coord.GetX_Ang();
		ymax1 = this->max_coord.GetY_Ang();
		zmax1 = this->max_coord.GetZ_Ang();
	}

	double xn2[3],rad33,chrgv4;
	
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if(phost_mset->p_save_opt_default->save_selected && !aptr->Selected())
			continue;
		xn2[0] = aptr->GetX();
		xn2[1] = aptr->GetY();
		xn2[2] = aptr->GetZ();
		rad33    = aptr->radius;
		chrgv4  = aptr->GetCharge();
		
		if( calc_range)
		{
			double tmp;
			tmp = xn2[0] - rad33; if(tmp < xmin1) xmin1 = tmp;
			tmp = xn2[1] - rad33; if(tmp < ymin1) ymin1 = tmp;
			tmp = xn2[2] - rad33; if(tmp < zmin1) zmin1 = tmp;
			tmp = xn2[0] + rad33; if(tmp > xmax1) xmax1 = tmp;
			tmp = xn2[1] + rad33; if(tmp > ymax1) ymax1 = tmp;
			tmp = xn2[2] + rad33; if(tmp > zmax1) zmax1 = tmp;
		}
		na++;
	}
	// Add boundary atoms
	//double gr_cent[3];
	gr_cent[0] = (xmin1 + xmax1)/2.0;
	gr_cent[1] = (ymin1 + ymax1)/2.0;
	gr_cent[2] = (zmin1 + zmax1)/2.0;

	double range = (xmax1 - xmin1);
	if( (ymax1 - ymin1) > range ) range = (ymax1 - ymin1);
	if( (zmax1 - zmin1) > range ) range = (zmax1 - zmin1);

	int ngrid = this->nx;

	range = range*100.0/this->perfil;

	float scale;

	if( this->perfil > 0.000001)  
	{
		scale = (ngrid - 1)/range;
	}
	else
	{
		scale = 0.5;
	}
	
	if(BoxSetMethod==SetBoxGridScale)
		scale=GridScale;
	else// i.e BoxSetMethod==SetBoxPerFil
		GridScale=scale;
	
	if(MoleculePosition==MoveToCenter)
	{
		double gxmin = (gr_cent[0] - range*0.5);
		double gymin = (gr_cent[1] - range*0.5);
		double gzmin = (gr_cent[2] - range*0.5);
		double gxmax = (gr_cent[0] + range*0.5);
		double gymax = (gr_cent[1] + range*0.5);
		double gzmax = (gr_cent[2] + range*0.5);
			
		this->el_pot_map.SetDimensions(nx,ny,nz);
		this->el_pot_map.SetGridCornersCoord( gxmin, gymin, gzmin,gxmax, gymax, gzmax); 
	}
	else//e.i. DontMove or UseOffset
	{
		double xrange=double(nx)/GridScale;
		double yrange=double(ny)/GridScale;
		double zrange=double(nz)/GridScale;
		double gxmin = - xrange*0.5;
		double gymin = - yrange*0.5;
		double gzmin = - zrange*0.5;
		double gxmax = xrange*0.5;
		double gymax = yrange*0.5;
		double gzmax = zrange*0.5;
			
		this->el_pot_map.SetDimensions(nx,ny,nz);
		this->el_pot_map.SetGridCornersCoord( gxmin, gymin, gzmin, gxmax, gymax, gzmax); 
	}
	if(m_ContWorld!=NULL)
	{
		m_ContWorld->Potential=NULL;
		//if(bConvIndConc)
		//	m_ContWorld->C=NULL;
		delete m_ContWorld;
	}
	m_ContWorld=new ContWorld();
	PNP_EXIT_FAIL_NULL(m_ContWorld,"Cannot initialize ContWorld\n");
	
	if(rionst==0.0)
	{
		PNP_EXIT_ON_FAIL_MES(m_ContWorld->SetContWorldNoMobIons(nx, ny, nz, scale, PBCX,PBCY, PBCZ),"Cannot Set Parameters for ContWorld");
	}
	else
	{
		PNP_EXIT_ON_FAIL_MES(m_ContWorld->SetContWorldTwoMobIons(nx, ny, nz, scale, PBCX,PBCY, PBCZ,1.0,-1.0),"Cannot Set Parameters for ContWorld");
	}
	
	m_ContWorld->Potential=this->el_pot_map.GetFieldPtr();
	pnpPrint("GridSize [%d %d %d] GridScale %f\n",m_ContWorld->GridSize[0], m_ContWorld->GridSize[1], m_ContWorld->GridSize[2],m_ContWorld->GridScale);
	
	//m_ContWorld->
	return EXIT_SUCCESS;
}
void ElMod::SetMemTube(double m_x, double m_y, double m_z0, double m_z1, double m_R0, double m_R1, double m_eps)
{
	MembraneType=MemT_Tube;
	memX=m_x;
	memY=m_y;
	memZ0=m_z0;
	memZ1=m_z1;
	memR0=m_R0;
	memR1=m_R1;
	epsmem=m_eps;
}
void ElMod::SetMemNone()
{
	MembraneType=MemT_None;
}
int ElMod::BuildWorld()
{
	int i;
	BuildWorldNI* m_BuildWorldNI;
	if(BldType==BldNI)
	{
		m_BuildWorldNI=new BuildWorldNI();
		m_BuildWorldNI->BuildUsingGPU=GPUAccelerated;
		m_BuildWorldNI->MemoryLimitOnOneArray=MemoryLimitOnOneArray;
	}
	else if(BldType==BldEu)
		m_BuildWorldNI=new BuildWorldEu();
	else if(BldType==BldCmp)
		m_BuildWorldNI=new BuildWorldCmp();
	
	m_BuildWorldNI->Rwat=radprb;
	m_BuildWorldNI->Rsmooth=0.0f;
	
	if(bPBSR)
		m_BuildWorldNI->NIonsTypes=2;
	else
		m_BuildWorldNI->NIonsTypes=1;
	
	m_BuildWorldNI->IonsR=new float[m_BuildWorldNI->NIonsTypes];
	for(i=0;i<m_BuildWorldNI->NIonsTypes;i++)
		m_BuildWorldNI->IonsR[i]=exrad;
	
	m_BuildWorldNI->DielNum=4;
	m_BuildWorldNI->DielConst=new float[m_BuildWorldNI->DielNum];
	m_BuildWorldNI->DielConst[0]=epsout;
	m_BuildWorldNI->DielConst[1]=epsi;
	m_BuildWorldNI->DielConst[2]=epsmem;
	m_BuildWorldNI->DielConst[3]=probe_eps;
	
	//check which dielectric constants are the same
	int iEpsOut=1,iEpsIn=2,iEpsMem=3,iEpsProbe=4;
	iEpsOut=1;
	if(fabs(epsi-epsout)<1.e-3)
		iEpsIn=iEpsOut;
	if(MembraneType!=MemT_None)
	{
		if(fabs(epsi-epsmem)<1.e-3)
			iEpsMem=iEpsIn;
		if(fabs(epsout-epsmem)<1.e-3)
			iEpsMem=iEpsOut;
	}
	if(fabs(probe_eps-epsout)<1.e-3)
		iEpsProbe=iEpsOut;
	if(fabs(probe_eps-epsi)<1.e-3)
		iEpsProbe=iEpsIn;
	if(MembraneType!=MemT_None)
		if(fabs(probe_eps-epsmem)<1.e-3)
			iEpsProbe=iEpsMem;
	
	
	m_BuildWorldNI->DiffusionNum=3;
	m_BuildWorldNI->DiffusionConst=new float[m_BuildWorldNI->DiffusionNum];
	m_BuildWorldNI->DiffusionConst[0]=20.0;
	m_BuildWorldNI->DiffusionConst[1]=20.0;
	m_BuildWorldNI->DiffusionConst[2]=0.0;
	
	
	m_BuildWorldNI->DiffusionMode=BuildWorldNI::Plain;
	
	if(boundary==0)//0=Zero
		m_BuildWorldNI->BoundaryCondition=BuildWorldNI::ZeroBC;
	else if(boundary==1 || boundary==3)//1=Apprx Coul,3=Coul
		m_BuildWorldNI->BoundaryCondition=BuildWorldNI::CoulBC;
	else if(boundary==2)//Focusing, focusing is done in elmod
		m_BuildWorldNI->BoundaryCondition=BuildWorldNI::ZeroBC;
	else
	{
		pnpError("ElMod cannot use this boundary mode, will use Coul BC\n");
		m_BuildWorldNI->BoundaryCondition=BuildWorldNI::CoulBC;
	}
	
	m_BuildWorldNI->MakeDielectricMap=true;
	if(rionst==0.0)
	{
		m_BuildWorldNI->MakeDiffusionMap=false;
		m_BuildWorldNI->MakeConcentrationMap=false;
	}
	else
	{
		m_BuildWorldNI->MakeDiffusionMap=true;
		m_BuildWorldNI->MakeConcentrationMap=true;
	}
	if(bPBSR)
	{
		m_BuildWorldNI->MakeLJRepultion=true;
		m_BuildWorldNI->IonsHalfSigma=new float[2];
		m_BuildWorldNI->IonsFourEpsilon=new float[2];
		
		m_BuildWorldNI->IonsHalfSigma[0]=IonsHalfSigma.GetVal_idx0(0);
		m_BuildWorldNI->IonsHalfSigma[1]=IonsHalfSigma.GetVal_idx0(1);
		m_BuildWorldNI->IonsFourEpsilon[0]=IonsFourEpsilon.GetVal_idx0(0);
		m_BuildWorldNI->IonsFourEpsilon[1]=IonsFourEpsilon.GetVal_idx0(1);
	}
	m_BuildWorldNI->MakeChargeMap=true;
	m_BuildWorldNI->RemovingCavities=true;
	//m_BuildWorldNI->NoStaticCharges=false;
	
	m_BuildWorldNI->RemovingCavitiesOnDielectricMap=RemovingCavitiesOnDielectricMap;
	m_BuildWorldNI->RemCavOnDielFillWith=RemCavOnDielFillWith;
	m_BuildWorldNI->RemCavOnDielWhere2Look=RemCavOnDielWhere2Look;
	//Load Elements
	//Bulk Parameters
	m_BuildWorldNI->BulkParam=new GenericGeometricalObject();
	m_BuildWorldNI->BulkParam->SetName("BulkParameters");
	m_BuildWorldNI->BulkParam->Epsilon=iEpsOut;
	m_BuildWorldNI->BulkParam->IonsD=new int[m_BuildWorldNI->NIonsTypes];
	if(bPBSR)
	{
		m_BuildWorldNI->BulkParam->IonsD[0]=1;
		m_BuildWorldNI->BulkParam->IonsD[1]=2;
	}
	else
		m_BuildWorldNI->BulkParam->IonsD[0]=1;
	
	m_BuildWorldNI->BulkParam->C=new float[m_BuildWorldNI->NIonsTypes];
	for(i=0;i<m_BuildWorldNI->NIonsTypes;i++)
		m_BuildWorldNI->BulkParam->C[i]=rionst;
	m_BuildWorldNI->BulkParam->InternalUnit=false;
	m_BuildWorldNI->BulkParam->old_world=NULL;
	
	//Others Elements
	//GOAtom
	if(bCalMolSet)
	{
		GOAtoms* m_GOAtoms=new GOAtoms();
		m_GOAtoms->Epsilon=iEpsIn;
		m_GOAtoms->IonsD=new int[m_BuildWorldNI->NIonsTypes];
		for(i=0;i<m_BuildWorldNI->NIonsTypes;i++)
			m_GOAtoms->IonsD[i]=3;
		m_GOAtoms->C=new float[m_BuildWorldNI->NIonsTypes];
		for(i=0;i<m_BuildWorldNI->NIonsTypes;i++)
			m_GOAtoms->C[i]=0.0;
		m_GOAtoms->InternalUnit=false;
		m_GOAtoms->old_world=NULL;
		if(MembraneType>0)
			m_GOAtoms->MakePreRoll=true;
		
		
		m_GOAtoms->ChargeDist=ChargeDist;
		m_GOAtoms->ChargeDistN=ChargeDistPar;
		m_GOAtoms->Offset[0]=0.0;
		m_GOAtoms->Offset[1]=0.0;
		m_GOAtoms->Offset[2]=0.0;
		
		HaAtom* aptr;
		AtomIteratorMolSet aitr(phost_mset);
		
		if(MoleculePosition==MoveToCenter)
		{
			for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
			{
				if(phost_mset->p_save_opt_default->save_selected && !aptr->Selected())
					continue;
				
				m_GOAtoms->r[0].push_back(float(aptr->GetX() - gr_cent[0]));
				m_GOAtoms->r[1].push_back(float(aptr->GetY() - gr_cent[1]));
				m_GOAtoms->r[2].push_back(float(aptr->GetZ() - gr_cent[2]));
				m_GOAtoms->q.push_back(float(aptr->GetCharge()));
				m_GOAtoms->R.push_back(float(aptr->radius));
			}
			if(bInsertProbe&&iEpsProbe==iEpsIn)
			{
				m_GOAtoms->r[0].push_back(float(probe_x - gr_cent[0]));
				m_GOAtoms->r[1].push_back(float(probe_y - gr_cent[1]));
				m_GOAtoms->r[2].push_back(float(probe_z - gr_cent[2]));
				m_GOAtoms->q.push_back(float(probe_q));
				m_GOAtoms->R.push_back(float(probe_R));
			}
		}
		else if(MoleculePosition==DontMove)
		{
			for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
			{
				if(phost_mset->p_save_opt_default->save_selected && !aptr->Selected())
					continue;
				
				m_GOAtoms->r[0].push_back(float(aptr->GetX()));
				m_GOAtoms->r[1].push_back(float(aptr->GetY()));
				m_GOAtoms->r[2].push_back(float(aptr->GetZ()));
				m_GOAtoms->q.push_back(float(aptr->GetCharge()));
				m_GOAtoms->R.push_back(float(aptr->radius));
			}
			if(bInsertProbe&&iEpsProbe==iEpsIn)
			{
				m_GOAtoms->r[0].push_back(float(probe_x));
				m_GOAtoms->r[1].push_back(float(probe_y));
				m_GOAtoms->r[2].push_back(float(probe_z));
				m_GOAtoms->q.push_back(float(probe_q));
				m_GOAtoms->R.push_back(float(probe_R));
			}
		}
		else if(MoleculePosition==UseOffset)
		{
			for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
			{
				if(phost_mset->p_save_opt_default->save_selected && !aptr->Selected())
					continue;
				
				m_GOAtoms->r[0].push_back(float(aptr->GetX()-OffsetX));
				m_GOAtoms->r[1].push_back(float(aptr->GetY()-OffsetY));
				m_GOAtoms->r[2].push_back(float(aptr->GetZ()-OffsetZ));
				m_GOAtoms->q.push_back(float(aptr->GetCharge()));
				m_GOAtoms->R.push_back(float(aptr->radius));
			}
			if(bInsertProbe&&iEpsProbe==iEpsIn)
			{
				m_GOAtoms->r[0].push_back(float(probe_x-OffsetX));
				m_GOAtoms->r[1].push_back(float(probe_y-OffsetY));
				m_GOAtoms->r[2].push_back(float(probe_z-OffsetZ));
				m_GOAtoms->q.push_back(float(probe_q));
				m_GOAtoms->R.push_back(float(probe_R));
			}
		}
		else if(MoleculePosition==TranslateAndRotate)
		{
			Vec3D vPhi(0.0,0.0,1.0);
			Vec3D vTheta(0.0,1.0,0.0);
			double cosPhi=cos(Phi);
			double sinPhi=sin(Phi);
			double cosTheta=cos(-Theta);
			double sinTheta=sin(-Theta);
			
			for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
			{
				if(phost_mset->p_save_opt_default->save_selected && !aptr->Selected())
					continue;
				Vec3D m_r(aptr->GetX()+OffsetX,aptr->GetY()+OffsetY,aptr->GetZ()+OffsetZ);
				
				m_r.Rotate(vTheta, cosTheta, sinTheta);
				m_r.Rotate(vPhi, cosPhi, sinPhi);
				
				m_GOAtoms->r[0].push_back(float(m_r.GetX()));
				m_GOAtoms->r[1].push_back(float(m_r.GetY()));
				m_GOAtoms->r[2].push_back(float(m_r.GetZ()));
				
				m_GOAtoms->q.push_back(float(aptr->GetCharge()));
				m_GOAtoms->R.push_back(float(aptr->radius));
			}
			if(bInsertProbe&&iEpsProbe==iEpsIn)
			{
				Vec3D m_r(probe_x+OffsetX, probe_y+OffsetY, probe_z+OffsetZ );
				
				m_r.Rotate(vTheta, cosTheta, sinTheta);
				m_r.Rotate(vPhi, cosPhi, sinPhi);
				
				m_GOAtoms->r[0].push_back(float(m_r.GetX()));
				m_GOAtoms->r[1].push_back(float(m_r.GetY()));
				m_GOAtoms->r[2].push_back(float(m_r.GetZ()));
				m_GOAtoms->q.push_back(float(probe_q));
				m_GOAtoms->R.push_back(float(probe_R));
			}
		}
		m_GOAtoms->NAtoms=m_GOAtoms->R.size();
		
		if(PQRFile!="")
		{
			PrintLog("Load Q and R from file: %s\n",PQRFile.c_str());
			m_GOAtoms->LoadPQRonlyQR(PQRFile.c_str());
		}
		if(bPBSR)
		{
			if(PREFile!="")
			{
				PrintLog("Load params for SR from file: %s\n",PREFile.c_str());
				m_GOAtoms->LoadPRE(PREFile.c_str());
			}
			else
			{
				PrintLog("ERROR no file with params for SR is given\n");
				return EXIT_FAILURE;
			}
		}
		
		m_BuildWorldNI->GOElms.push_back(m_GOAtoms);
	}
	//insert probe ion
	if((bInsertProbe&&iEpsProbe!=iEpsIn)||!bCalMolSet)
	{
		GOAtoms* m_GOProbe=new GOAtoms();
		m_GOProbe->Epsilon=iEpsProbe;
		m_GOProbe->IonsD=new int[m_BuildWorldNI->NIonsTypes];
		for(i=0;i<m_BuildWorldNI->NIonsTypes;i++)
			m_GOProbe->IonsD[i]=3;
		m_GOProbe->C=new float[m_BuildWorldNI->NIonsTypes];
		for(i=0;i<m_BuildWorldNI->NIonsTypes;i++)
			m_GOProbe->C[i]=0.0;
		m_GOProbe->InternalUnit=false;
		m_GOProbe->old_world=NULL;
		//if(MembraneType>0)
		//	m_GOProbe->MakePreRoll=true;
	
	
		m_GOProbe->ChargeDist=ChargeDist;
		m_GOProbe->ChargeDistN=ChargeDistPar;
		m_GOProbe->Offset[0]=0.0;
		m_GOProbe->Offset[1]=0.0;
		m_GOProbe->Offset[2]=0.0;
			
		if(MoleculePosition==MoveToCenter)
		{
			m_GOProbe->r[0].push_back(float(probe_x-gr_cent[0]));
			m_GOProbe->r[1].push_back(float(probe_y-gr_cent[1]));
			m_GOProbe->r[2].push_back(float(probe_z-gr_cent[2]));
			m_GOProbe->q.push_back(float(probe_q));
			m_GOProbe->R.push_back(float(probe_R));
		}
		else if(MoleculePosition==DontMove)
		{
			m_GOProbe->r[0].push_back(float(probe_x));
			m_GOProbe->r[1].push_back(float(probe_y));
			m_GOProbe->r[2].push_back(float(probe_z));
			m_GOProbe->q.push_back(float(probe_q));
			m_GOProbe->R.push_back(float(probe_R));
		}
		else if(MoleculePosition==UseOffset)
		{
			m_GOProbe->r[0].push_back(float(probe_x-OffsetX));
			m_GOProbe->r[1].push_back(float(probe_y-OffsetY));
			m_GOProbe->r[2].push_back(float(probe_z-OffsetZ));
			m_GOProbe->q.push_back(float(probe_q));
			m_GOProbe->R.push_back(float(probe_R));
		}
		else if(MoleculePosition==TranslateAndRotate)
		{
			Vec3D vPhi(0.0,0.0,1.0);
			Vec3D vTheta(0.0,1.0,0.0);
			double cosPhi=cos(Phi);
			double sinPhi=sin(Phi);
			double cosTheta=cos(-Theta);
			double sinTheta=sin(-Theta);
		
			Vec3D m_r( probe_x+OffsetX, probe_y+OffsetY, probe_z+OffsetZ);
		
			m_r.Rotate(vTheta, cosTheta, sinTheta);
			m_r.Rotate(vPhi, cosPhi, sinPhi);
		
			m_GOProbe->r[0].push_back(float(m_r.GetX()));
			m_GOProbe->r[1].push_back(float(m_r.GetY()));
			m_GOProbe->r[2].push_back(float(m_r.GetZ()));
			m_GOProbe->q.push_back(float(probe_q));
			m_GOProbe->R.push_back(float(probe_R));
		}
		m_GOProbe->NAtoms=m_GOProbe->R.size();
		m_BuildWorldNI->GOElms.push_back(m_GOProbe);
		
	}
	//Membrane
	if(MembraneType==MemT_MembraneZ)
	{
		GOMembraneZ* membraneZ=new GOMembraneZ();
		membraneZ->Epsilon=iEpsMem;
		membraneZ->IonsD=new int[m_BuildWorldNI->NIonsTypes];
		for(i=0;i<m_BuildWorldNI->NIonsTypes;i++)
			membraneZ->IonsD[i]=3;
		membraneZ->C=new float[m_BuildWorldNI->NIonsTypes];
		for(i=0;i<m_BuildWorldNI->NIonsTypes;i++)
			membraneZ->C[i]=0.0;
		membraneZ->InternalUnit=false;
		membraneZ->old_world=NULL;
		
		membraneZ->Z[0]=memZ0;
		membraneZ->Z[1]=memZ1;
		
		m_BuildWorldNI->GOElms.push_back(membraneZ);
	}
	else if(MembraneType==MemT_Tube)
	{
		GOTube* tube=new GOTube();
		tube->Epsilon=iEpsMem;
		tube->IonsD=new int[m_BuildWorldNI->NIonsTypes];
		for(i=0;i<m_BuildWorldNI->NIonsTypes;i++)
			tube->IonsD[i]=3;
		tube->C=new float[m_BuildWorldNI->NIonsTypes];
		for(i=0;i<m_BuildWorldNI->NIonsTypes;i++)
			tube->C[i]=0.0;
		tube->InternalUnit=false;
		tube->old_world=NULL;
		
		if(MoleculePosition==TranslateAndRotate)
		{
			Vec3D vPhi(0.0,0.0,1.0);
			Vec3D vTheta(0.0,1.0,0.0);
			double cosPhi=cos(Phi);
			double sinPhi=sin(Phi);
			double cosTheta=cos(-Theta);
			double sinTheta=sin(-Theta);
			
			Vec3D m_r0(memX+OffsetX, memY+OffsetY, memZ0+OffsetZ );
			Vec3D m_r1(memX+OffsetX, memY+OffsetY, memZ1+OffsetZ );
			
			if(Theta!=0.0)
			{
				printf("ERROR: Theta!=0.0 can not do this for MemT_Tube\n");
			}
			m_r0.Rotate(vPhi, cosPhi, sinPhi);
			m_r1.Rotate(vPhi, cosPhi, sinPhi);
			
			tube->XY[0]=float(m_r0.GetX());
			tube->XY[1]=float(m_r0.GetY());
			tube->Z[0]=float(m_r0.GetZ());
			tube->Z[1]=float(m_r1.GetZ());
			tube->R[0]=memR0;
			tube->R[1]=memR1;
		}
		else
		{
			tube->Z[0]=memZ0;
			tube->Z[1]=memZ1;
			tube->R[0]=memR0;
			tube->R[1]=memR1;
			tube->XY[0]=memX;
			tube->XY[1]=memY;
		}
		
		m_BuildWorldNI->GOElms.push_back(tube);
	}
	
// 	double RotateAroundZby=0.0;
// 	if(RotateAroundZby!=0.0)
// 	{
// 		pnpPrint0("Will rotate all elements arount z by angle %g\n",RotateAroundZby);
// 		RotateAroundZby=RotateAroundZby*M_PI/180.0;
// 		double VecZ[3]={0.0,0.0,0.0};
// 		int iElt;
// 		m_GOAtoms->RotateGGO(VecZ,cos(RotateAroundZby),sin(RotateAroundZby));
// 	}
	
	m_BuildWorldNI->BuildContWorld(m_ContWorld);
	delete m_BuildWorldNI;
	
	//m_ContWorld->WriteDiffusion("TMP_Diff.bin");
	//m_ContWorld->WriteNodeIndexing("TMP_NI.gz");
	//m_ContWorld->WritePMF("TMP_PMF.bin");
	return EXIT_SUCCESS;
}
int ElMod::SolvePoisson()
{
	PoissonSolver* m_PoissonSolver=new PoissonSolver();
	
	m_PoissonSolver->MinIterations=0;
// 	if(nlit<=0)
// 	{
// 		m_PoissonSolver->MaxIterations=300;
// 		pnpError("ElMod Cannot estimate number of iteration needed, will try 300\n");
// 	}
// 	else
// 	{
		m_PoissonSolver->MaxIterations=nlit;
// 	}
	m_PoissonSolver->Convergence=Convergence;
	m_PoissonSolver->Relaxation=Relaxation;
	m_PoissonSolver->ConvergenceCheck=20;
	m_PoissonSolver->solver=0;
	m_PoissonSolver->verbose=true;
	m_PoissonSolver->QmobMod=2;	
	
	m_PoissonSolver->ShowParameters();
	//Scale Parameters to Internal Units
	//Calculate Derivative Parameters
	m_PoissonSolver->om2 = m_PoissonSolver->Relaxation;
	m_PoissonSolver->om1 = 1.0-m_PoissonSolver->om2;
	m_PoissonSolver->om2d6 = m_PoissonSolver->om2/6.0;
	
	
	m_PoissonSolver->SetContWorld(m_ContWorld);
	m_PoissonSolver->InitSolver();
	m_PoissonSolver->Solve();
	delete m_PoissonSolver;
	
	tot_ene=m_ContWorld->SystemEnergy;
	return EXIT_SUCCESS;
}
int ElMod::SolvePoissonBoltzmann()
{
	PoissonBoltzmannSolver* m_PoissonBoltzmannSolver=new PoissonBoltzmannSolver();
	
// 	if(nlit<=0)
// 	{
// 		m_PoissonBoltzmannSolver->MaxIterationsLPB=300;
// 		pnpError("ElMod Cannot estimate number of iteration needed, will try 300\n");
// 	}
// 	else
// 	{
		m_PoissonBoltzmannSolver->MaxIterationsLPB=nlit;
// 	}
	m_PoissonBoltzmannSolver->MaxIterationsNPB=nnit;
	m_PoissonBoltzmannSolver->Convergence=Convergence;
	m_PoissonBoltzmannSolver->Relaxation=Relaxation;
	m_PoissonBoltzmannSolver->ConvergenceCheck=20;
	m_PoissonBoltzmannSolver->verbose=true;
	
	m_PoissonBoltzmannSolver->ShowParameters();
	//Scale Parameters to Internal Units
	//Calculate Derivative Parameters
	m_PoissonBoltzmannSolver->om2 = m_PoissonBoltzmannSolver->Relaxation;
	m_PoissonBoltzmannSolver->om1 = 1.0-m_PoissonBoltzmannSolver->om2;
	m_PoissonBoltzmannSolver->om2d6 = m_PoissonBoltzmannSolver->om2/6.0;
	
	
	m_PoissonBoltzmannSolver->SetContWorld(m_ContWorld);
	m_PoissonBoltzmannSolver->InitSolver();
	m_PoissonBoltzmannSolver->Solve();
	
	if(bConvIndConc)
	{
		
		double gxmin = this->el_pot_map.GetXmin();
		double gymin = this->el_pot_map.GetYmin();
		double gzmin = this->el_pot_map.GetZmin();
		double gxmax = this->el_pot_map.GetXmax();
		double gymax = this->el_pot_map.GetYmax();
		double gzmax = this->el_pot_map.GetZmax();
		
		this->ConcMap0.SetDimensions(nx,ny,nz);
		this->ConcMap0.SetGridCornersCoord( gxmin, gymin, gzmin, gxmax, gymax, gzmax);
		this->ConcMap1.SetDimensions(nx,ny,nz);
		this->ConcMap1.SetGridCornersCoord( gxmin, gymin, gzmin, gxmax, gymax, gzmax);
		
		if(m_ContWorld->C!=NULL)
		{
			int i;
			for(i=0;i<2;i++)
			{
				if(m_ContWorld->C[i]!=NULL)
					delete [] m_ContWorld->C[i];
			}
			delete [] m_ContWorld->C;
		}
		m_ContWorld->C=new float*[2];
		m_ContWorld->C[0]=this->ConcMap0.GetFieldPtr();
		m_ContWorld->C[1]=this->ConcMap1.GetFieldPtr();
		PNPUtil::ConvertPBLJresultsToDynamicCharge(m_ContWorld);
		delete [] m_ContWorld->C;
		m_ContWorld->C=NULL;
		float fpoh= 4.0*M_PI*m_ContWorld->GridScale;
		float coef=fpoh*COANGS/(m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
		this->ConcMap0.MultiplyByValues(1.0/coef);
		this->ConcMap1.MultiplyByValues(1.0/coef);
	}
	
	delete m_PoissonBoltzmannSolver;
	
	tot_ene=m_ContWorld->SystemEnergy;
	return EXIT_SUCCESS;
}
int ElMod::SolvePoissonBoltzmannSoftRepultion()
{
	PBwithLJSolver* m_PBSR=new PBwithLJSolver();
	
	m_PBSR->MaxIterations=nlit;
	m_PBSR->Convergence=Convergence;
	m_PBSR->Relaxation=Relaxation;
	m_PBSR->ConvergenceCheck=20;
	m_PBSR->verbose=true;
	
	m_PBSR->ShowParameters();
	m_PBSR->SetRelaxation(m_PBSR->Relaxation);
	
	m_PBSR->SetContWorld(m_ContWorld);
	m_PBSR->InitSolver();
	m_PBSR->Solve();
	
	if(bConvIndConc)
	{
		
		double gxmin = this->el_pot_map.GetXmin();
		double gymin = this->el_pot_map.GetYmin();
		double gzmin = this->el_pot_map.GetZmin();
		double gxmax = this->el_pot_map.GetXmax();
		double gymax = this->el_pot_map.GetYmax();
		double gzmax = this->el_pot_map.GetZmax();
		
		this->ConcMap0.SetDimensions(nx,ny,nz);
		this->ConcMap0.SetGridCornersCoord( gxmin, gymin, gzmin, gxmax, gymax, gzmax);
		this->ConcMap1.SetDimensions(nx,ny,nz);
		this->ConcMap1.SetGridCornersCoord( gxmin, gymin, gzmin, gxmax, gymax, gzmax);
		
		if(m_ContWorld->C!=NULL)
		{
			int i;
			for(i=0;i<2;i++)
			{
				if(m_ContWorld->C[i]!=NULL)
					delete [] m_ContWorld->C[i];
			}
			delete [] m_ContWorld->C;
		}
		m_ContWorld->C=new float*[2];
		m_ContWorld->C[0]=this->ConcMap0.GetFieldPtr();
		m_ContWorld->C[1]=this->ConcMap1.GetFieldPtr();
		PNPUtil::ConvertPBLJresultsToDynamicCharge(m_ContWorld);
		delete [] m_ContWorld->C;
		m_ContWorld->C=NULL;
		float fpoh= 4.0*M_PI*m_ContWorld->GridScale;
		float coef=fpoh*COANGS/(m_ContWorld->GridScale*m_ContWorld->GridScale*m_ContWorld->GridScale);
		this->ConcMap0.MultiplyByValues(1.0/coef);
		this->ConcMap1.MultiplyByValues(1.0/coef);
	}
	
	delete m_PBSR;
	
	tot_ene=m_ContWorld->SystemEnergy;
	return EXIT_SUCCESS;
}
int ElMod::SetNinaBelowRouxRadii(bool only_selected_atoms)
{
	PNP_EXIT_FAIL_NULL(phost_mset,"ElMod does not assign to any MolSet\n");
	
	PrintLog("Setting the atomic radii as described at:\n");
	PrintLog("\tM. Nina, D. Beglov and B. Roux ,\n");
	PrintLog("\tAtomic radii for continuum electrostatics calculations\n");
	PrintLog("\tbased on molecular dynamics free energy simulations.\n");
	PrintLog("\tJ. Phys. Chem. B 101 (1997), pp. 5239–5248.\n");
	
	HaAtom* aptr;
	HaResidue* rptr;
	AtomIteratorMolSet aitr(phost_mset);

	for(aptr= aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if(only_selected_atoms && !aptr->Selected())
			continue;
		rptr=aptr->GetHostRes();
		std::string AtomName=aptr->GetName();
		std::string ResName=rptr->GetName();
		std::string ResModName=rptr->GetNameModifier();
		int elno = aptr->GetElemNo();
		//double R=-1.0;
		aptr->radius=-1.0;
		//Adopted from
		//  http://thallium.bsd.uchicago.edu/RouxLab/downloads/charmm/pbeq/radius.str
		// 	* Atomic radius derived from solvent electrostatic charge distribution
		// 	* Tested with free energy perturbation with explicit solvent
		// 	* Authors:  Mafalda Nina, Dmitrii Belogv, and Benoit Roux
		// 	* University of Montreal, June 1996.
		// 	* M. Nina and B. Roux. Atomic Radii for Continuum Electrostatics Calculations based on 
		// 	* Molecular Dynamics Free Energy Simulations. J. Phys. Chem. B 101: 5239-5248 (1997).
		// 	*
		// 	
		// 	prnlev 0
		// 	
		// 	! UPDATES:
		// 	! --------
		// 	!
		// 	! GLU and ASP modified December 1st, 1998 by Mafalda Nina
		// 	!
		// 	! Protonated histidine HSP has been included, January 1999 by Mafalda Nina
		// 	! dG_elec = -68.15 kcal/mol (MD/FES) or -68.10 kcal/mol (PBEQ)
		// 	!
		// 	! TEA and ions added by Benoit Roux, January 1999.
		// 	!
		// 	! sodium added by Benoit Roux, November 2000
		//Set to zero all H radii
		if(aptr->IsHydrogen())
		{
			aptr->radius=0.0;
			continue;
		}
		if(AtomName.StartsWith("H"))
		{
			PrintLog("WARNING: will assume that %d %s of Res %d %s is Hydrogen\n",aptr->GetSerNo(),aptr->GetName(),rptr->GetSerNo(),rptr->GetName());
			aptr->radius=0.0;
			continue;
		}

		//Patches CT3 N-Methylamide C-terminus
		//        ACE acetylated N-terminus (ACP for PRO)
		//scalar wmain set 2.06 sele (type CAY .or. type CAT) end
		//scalar wmain set 2.04 sele type CY end
		//scalar wmain set 1.52 sele type OY end
		//scalar wmain set 2.23 sele type NT end
		//scalar wmain set 1.40 sele type OT* end ! for COO- terminus
		if(ResModName=="NT")
		{
			if(AtomName=="N")
				aptr->radius = 2.23;
		}
		if(ResModName=="CT")
		{
			if((AtomName=="O")||(AtomName=="OXT"))
				aptr->radius = 1.4;
		}
		if(aptr->IsAminoBackbone())//Backbone
		{
			if(AtomName=="C")
				aptr->radius = 2.04;
			else if(AtomName=="O")
				aptr->radius = 1.52;
			else if(AtomName=="N")
				aptr->radius = 2.23;
			else if(AtomName=="CA")
			{
				if(ResName=="GLY")
					aptr->radius = 2.38;
				else
					aptr->radius = 2.86;
			}
		}
		else//i.e side chain and not hydrogen
		{
			//Carbons 
			if(AtomName=="CB")
				aptr->radius = 2.67;
			if(AtomName.StartsWith("CG"))//CG*
			{
				if((ResName=="ARG") || (ResName=="GLN") || (ResName=="ILE") || (ResName=="LYS") || (ResName=="MET") || (ResName=="PHE") || (ResName=="THR") || (ResName=="TRP") || (ResName=="VAL") || (ResName=="LEU") || (ResName=="TYR"))//probably LEU ant TYR also set here
				{
					aptr->radius = 2.46;
				}
				if(ResName=="HIS")
				{
					if(ResModName=="PROT")
						aptr->radius = 2.46;
					else if(ResModName=="EPSILON")
						aptr->radius = 2.46;
					else
						aptr->radius = 2.46;
				}
			}
			if((AtomName=="CG")&&(ResName=="GLU"))
				aptr->radius = 2.77;
			if(AtomName.StartsWith("CD"))//CD*
			{
				if((ResName=="ARG") || (ResName=="ILE") || (ResName=="LEU") || (ResName=="LYS"))
				{
					aptr->radius = 2.44;
				}
			}
			if((AtomName=="CD")||(AtomName=="CG"))
				if((ResName=="ASP") || (ResName=="GLU") || (ResName=="ASN") || (ResName=="GLN"))
					aptr->radius = 1.98;
			if(ResName=="PRO")
			{
				if((AtomName=="CB")||(AtomName=="CG")||(AtomName=="CD"))
					aptr->radius = 1.98;
			}
			if((ResName=="TYR")||(ResName=="PHE"))
			{
				if((AtomName=="CE1")||(AtomName=="CE2")||(AtomName=="CD1")||(AtomName=="CD2")||(AtomName=="CZ"))
					aptr->radius = 2.00;
			}
			if(ResName=="TRP")
			{
				if((AtomName=="CE2")||(AtomName=="CE3")||(AtomName=="CD1")||(AtomName=="CD2")||(AtomName=="CZ2")||(AtomName=="CZ3")||(AtomName=="CH2"))
					aptr->radius = 1.78;
			}
			if(ResName=="MET"&&(AtomName=="CE"))
				aptr->radius = 2.1;
			if((AtomName=="CZ") && (ResName=="ARG"))
				aptr->radius = 2.8;
			if((AtomName=="CE") && (ResName=="LYS"))
				aptr->radius = 2.8;
			if(ResName=="HIS"&&((AtomName=="CE1")||(AtomName=="CD2")))
			{
				if(ResModName=="PROT")
					aptr->radius = 1.98;
				else if(ResModName=="EPSILON")
					aptr->radius = 1.98;
				else
					aptr->radius = 1.98;
			}
			//Oxygens
			if((ResName=="GLU") || (ResName=="ASP"))
			{
				if(AtomName.StartsWith("OE")||AtomName.StartsWith("OD"))
					aptr->radius = 1.40;
			}
			if((ResName=="ASN") || (ResName=="GLN"))
			{
				if(AtomName.StartsWith("OE")||AtomName.StartsWith("OD"))
					aptr->radius = 1.42;
			}
			if((ResName=="SER") || (ResName=="THR"))
			{
				if(AtomName.StartsWith("OG"))
					aptr->radius = 1.64;
			}
			if(ResName=="TYR")
			{
				if(AtomName=="OH")
					aptr->radius = 1.85;
			}
			if((ResName=="TIP3")||(ResName=="WAT"))//for explicit water molecules
			{
				if(AtomName[0]=='O')
					aptr->radius = 2.2;
			}
			//Nitrogens
			if(ResName=="HIS"&&((AtomName=="NE2")||(AtomName=="ND1")))
			{
				if(ResModName=="PROT")
					aptr->radius = 2.3;
				else if(ResModName=="EPSILON")
					aptr->radius = 1.8;
				else
					aptr->radius = 1.8;
			}
			if((ResName=="ARG")&&((AtomName.StartsWith("NH")||(AtomName=="NE"))))
					aptr->radius = 2.13;
			if((ResName=="LYS") && (AtomName=="NZ"))
				aptr->radius = 2.13;
			if((ResName=="GLN") && (AtomName=="NE2"))
				aptr->radius = 2.15;
			if((ResName=="ASN") && (AtomName=="ND2"))
				aptr->radius = 2.15;
			if((ResName=="TRP")&&(AtomName=="NE1"))
				aptr->radius = 2.4;
			//Sulphur
			if((ResName=="MET")||(ResName=="CYS"))
			{
				if(AtomName.StartsWith("S"))
					aptr->radius = 2.0;
			}
			//Ions
			//scalar wmain set 2.035 select resname POT end
			//	!potassium ion K+
			//scalar wmain set 2.035 select resname CLA end
			//	!chloride ion Cl-
			//scalar wmain set 1.66 select resname SOD end
			//	!sodium ion Na+
			//else if(AtomName=="O")
			//	aptr->radius = 1.52;
		}
		if(aptr->radius <0.0)
		{
			PrintLog("WARNING: Can not assing Radius for atom: %d %s of Res %d %s will apply default\n",aptr->GetSerNo(),aptr->GetName(),rptr->GetSerNo(),rptr->GetName());
			//Set heavy atoms to average default values 
		if(AtomName.StartsWith("C"))
			aptr->radius = 2.3;
		else if(AtomName.StartsWith("O"))
			aptr->radius = 1.8;
		else if(AtomName.StartsWith("N"))
			aptr->radius = 2.3;
		else if(AtomName.StartsWith("S"))
			aptr->radius = 2.3;
		else 
			PrintLog("WARNING: \tCan not assing default Radius.\n");
		}
	}
	return EXIT_SUCCESS;
}
int ElMod::LoadConcMaps(const char *filename)
{
	VectorField3D* VField=new VectorField3D(filename);
	int i,j;
	for(i=0;i<2;i++)
	{
		HaField3D* conc;
		if(i==0)conc=&(ConcMap0);
		else conc=&(ConcMap1);
		
		std::string Name="";
		Name<<"Conc["<<i<<"]";
		conc->SetName( Name.ToStdString() );
		conc->SetDimensions(VField->GridSize[0],VField->GridSize[1],VField->GridSize[2]);
		conc->SetCenterAsZero(VField->GridScale);
		
		int GS_XYZ=VField->GridSize[0]*VField->GridSize[1]*VField->GridSize[2];
		float *HVec=NULL;
		float *VPar=VField->V[i];
		HVec=conc->GetFieldPtr();
		for(j=0;j<GS_XYZ;j++)HVec[j]=VPar[j];
	}
	delete VField;
	return EXIT_SUCCESS;
}
int ElMod::LoadPQR(const char *filename, bool LoadXYZ,bool LoadQ,bool LoadR)
{
	FILE* fpqr = fopen(filename,"rt");

	PNP_EXIT_FAIL_NULL1(fpqr,"Can't read %s\n",filename);

	char Record[1024];
	char AtmName[5];
	char ResName[5];
	
	HaChain  *chain;
	HaResidue  *group;
	HaAtom  *aptr;
	int count=0;
	
	MoleculesType::iterator mol_itr;
	for( mol_itr=phost_mset->HostMolecules.begin(); mol_itr != phost_mset->HostMolecules.end(); mol_itr++)
	{
		ChainIteratorMolecule ch_itr(*mol_itr);
		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			for(group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes())
			{
				AtomIteratorAtomGroup aitr_group(group);
				for(aptr= aitr_group.GetFirstAtom(); aptr; aptr = aitr_group.GetNextAtom())
				{
					if( !phost_mset->p_save_opt_default->save_selected || aptr->Selected())
					{
						std::string atname(aptr->GetName());
						if(atname.size() < 4)
							atname.insert(0," ");

						int k;
						for(k=0; k < 3; k++)
						{
							if(atname.size() < 4)
								atname+= " ";
						}
						
						if(atname.size() > 4)
							atname = atname.substr(0,3);
						
						std::string res_name(group->GetName());

						for(k=0; k < 3; k++)
						{
							if(res_name.size() < 3)
								res_name.insert(0," ");
						}
						
						if(res_name.size() > 3)
							res_name = res_name.substr(0,2);
						
						bool HaveARecord=false;
						do
						{
							fgets(Record,1024,fpqr);
							if(strncmp("ATOM",Record,4)==0||strncmp("HETATM",Record,6)==0)
								HaveARecord=true;
						}
						while(!HaveARecord);
						
						count++;
						
						double tmpx,tmpy,tmpz,tmpq,tmpr;
						sscanf(Record+27,"%lf %lf %lf %lf %lf",&tmpx,&tmpy,&tmpz,&tmpq,&tmpr);
						
						strncpy(AtmName, Record+12, 4);AtmName[4]='\0';
						strncpy(ResName, Record+17, 3);ResName[3]='\0';
						//printf("a=|%s|%s|\n",AtmName,ResName);
						
						if(atname==AtmName && res_name==ResName)
						{
							aptr->SetX(tmpx);
							aptr->SetY(tmpy);
							aptr->SetZ(tmpz);
							aptr->SetCharge(tmpq);
							aptr->radius=tmpr;
						}
						else
						{
							pnpError("Atom(s)/Residue(s) names from %f do not match molset(%d: %s<->%s; %s<->%s).\n", filename, count, atname.c_str(), AtmName, res_name.c_str(), ResName);
						}
					}
				}
			}
		}
	}
	fclose(fpqr);
	return EXIT_SUCCESS;
}
int ElMod::SaveCoorToTmpArr()
{
	int NAtoms=phost_mset->GetNAtoms();
	DeleteCoorTmpArr();
	t_x=new double[NAtoms];
	t_y=new double[NAtoms];
	t_z=new double[NAtoms];
	
	HaAtom* aptr;
	AtomIteratorMolSet aitr(phost_mset);
	
	int count=0;
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if(phost_mset->p_save_opt_default->save_selected && !aptr->Selected())
			continue;
			
		t_x[count]=aptr->GetX();
		t_y[count]=aptr->GetY();
		t_z[count]=aptr->GetZ();
		
		count++;
	}
	return TRUE;
}
int ElMod::LoadCoorFromTmpArr()
{
	
	HaAtom* aptr;
	AtomIteratorMolSet aitr(phost_mset);
	
	int count=0;
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if(phost_mset->p_save_opt_default->save_selected && !aptr->Selected())
			continue;
			
		aptr->SetX(t_x[count]);
		aptr->SetY(t_y[count]);
		aptr->SetZ(t_z[count]);
		
		count++;
	}
	return TRUE;
}
int ElMod::DeleteCoorTmpArr()
{
	DeleteCArray(t_x);
	DeleteCArray(t_y);
	DeleteCArray(t_z);
	return TRUE;
}
int ElMod::RotateMolSet(Vec3D& n,double a)
{
	HaAtom* aptr;
	AtomIteratorMolSet aitr(phost_mset);
	
	double cosa=cos(a);
	double sina=sin(a);
	
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if(phost_mset->p_save_opt_default->save_selected && !aptr->Selected())
			continue;
			
		aptr->Rotate(n, cosa, sina);
	}
	return TRUE;
}
void ElMod::TranslateAndRotateMolSet(double m_x,double m_y,double m_z,double m_Phi,double m_Theta)
{
	double m_OffsetX=m_x;
	double m_OffsetY=m_y;
	double m_OffsetZ=m_z;
	Vec3D vPhi(0.0,0.0,1.0);
	Vec3D vTheta(0.0,1.0,0.0);
	double cosPhi=cos(m_Phi);
	double sinPhi=sin(m_Phi);
	double cosTheta=cos(-m_Theta);
	double sinTheta=sin(-m_Theta);
	
	HaAtom* aptr;
	AtomIteratorMolSet aitr(phost_mset);
	
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if(phost_mset->p_save_opt_default->save_selected && !aptr->Selected())
			continue;
		Vec3D m_r(aptr->GetX()+m_OffsetX, aptr->GetY()+m_OffsetY, aptr->GetZ()+m_OffsetZ);
		
		m_r.Rotate(vTheta, cosTheta, sinTheta);
		m_r.Rotate(vPhi, cosPhi, sinPhi);
		
		aptr->SetX(m_r.GetX());
		aptr->SetY(m_r.GetY());
		aptr->SetZ(m_r.GetZ());
	}
}
int ElMod::GetMinMaxWithR()
{
	AtomIteratorMolSet aitr(phost_mset);
	HaAtom* aptr=aitr.GetFirstAtom();
	
	double xmin=aptr->GetX()-aptr->radius;
	double ymin=aptr->GetY()-aptr->radius;
	double zmin=aptr->GetZ()-aptr->radius;
	
	double xmax=aptr->GetX()+aptr->radius;
	double ymax=aptr->GetY()+aptr->radius;
	double zmax=aptr->GetZ()+aptr->radius;
	
	for(aptr=aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if(phost_mset->p_save_opt_default->save_selected && !aptr->Selected())
			continue;
		if(xmin>aptr->GetX()-aptr->radius)xmin=aptr->GetX()-aptr->radius;
		if(ymin>aptr->GetY()-aptr->radius)ymin=aptr->GetY()-aptr->radius;
		if(zmin>aptr->GetZ()-aptr->radius)zmin=aptr->GetZ()-aptr->radius;
	
		if(xmax<aptr->GetX()+aptr->radius)xmax=aptr->GetX()+aptr->radius;
		if(ymax<aptr->GetY()+aptr->radius)ymax=aptr->GetY()+aptr->radius;
		if(zmax<aptr->GetZ()+aptr->radius)zmax=aptr->GetZ()+aptr->radius;
	}
	printf("MinMaxWithR is [%.7f %.7f %.7f]-[%.7f %.7f %.7f]\n", xmin,ymin,zmin,xmax,ymax,zmax);
	//aptr->radius=tmpr;
	return TRUE;
}
int ElMod::CalcRAroundZero()
{
	AtomIteratorMolSet aitr(phost_mset);
	HaAtom* aptr=aitr.GetFirstAtom();
	
	double R=0;
	double x,y,z;
	double f1,f2,Rtmp;
	
	for(aptr=aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if(phost_mset->p_save_opt_default->save_selected && !aptr->Selected())
			continue;
		f1=sqrt(aptr->GetX()*aptr->GetX()+aptr->GetY()*aptr->GetY()+aptr->GetZ()*aptr->GetZ());
		f2=(f1+aptr->radius)/f1;
		x=aptr->GetX()*f2;
		y=aptr->GetY()*f2;
		z=aptr->GetZ()*f2;
		
		Rtmp=sqrt(x*x+y*y+z*z);
		if(Rtmp>R)R=Rtmp;
	}
	printf("RAroundZero is %.7f\n", R);
	return TRUE;
}
#endif
////////////////////////////////////////////////////////////////////////////
pKaCalcMod::pKaCalcMod(MolSet* new_phost_mset)
	:HaCompMod(COMP_MOD_PKA_CALC,new_phost_mset)
{
	pHmin=1.0;
	pHmax=16.0;
	pHstep=0.5;
	pKaCalcMethod=SCF_MULTI_SITE_CALC;
	MC_PKA_CALC_N_mc_cyc=10000;
	SCF_PKA_CALC_max_iter = 10000;
	SCF_PKA_CALC_pop_err_max = 0.0001;
	pKaCalcMethodStr.push_back("SCF_MULTI_SITE_CALC");
	pKaCalcMethodStr.push_back("MC_MULTI_SITE_CALC");
	SaveIntermediatePotNNI=false;
	SaveIntermediateResults=false;

	MolSet* pmset = GetMolSet();
	p_prot_rdx_mod = pmset->GetProtonRedoxMod(true);
}


pKaCalcMod::~pKaCalcMod()
{
}
int pKaCalcMod::PrintResWithAltProtState()
{
	PNP_EXIT_FAIL_NULL(phost_mset,"pKaCalcMod does not assign to any MolSet\n");
	pnpPrint("PrintResWithAltProtState\n");
	int old_save_selected=phost_mset->save_opt_default.save_selected;
	phost_mset->save_opt_default.save_selected = FALSE;
	
	HaResidue* pres;
	ResidueIteratorMolSet ritr(phost_mset);
	pnpPrint("ResNum ResName NumAltSt\n");
	for(pres = ritr.GetFirstRes(); pres; pres = ritr.GetNextRes())
	{
		std::string full_res_name = pres->GetFullName();
		int nst = p_prot_rdx_mod->GetNumResAltChemStates(pres);
		if(nst>0)
			pnpPrint("%d %s %d\n",pres->GetSerNo(),full_res_name.c_str(),nst);
	}
	
	phost_mset->save_opt_default.save_selected=old_save_selected;
	return EXIT_SUCCESS;
}
int pKaCalcMod::PrintPopulation()
{
	PNP_EXIT_FAIL_NULL(phost_mset,"pKaCalcMod does not assign to any MolSet\n");
	pnpPrint("PrintPopulation\n");
	
	
	int ipH, NpH=roundf((pHmax-pHmin)/pHstep)+1;
	double pH;
	
	pnpPrint("pH");
	int i;
	for( ipH = 0; ipH < NpH; ipH ++ )
	{
		pH=pHmin+ipH*pHstep;
		pnpPrint("\t%.4g",pH);
		
	}
	pnpPrint("\n");
	for(i=0; i < NumberOfAltStates; i++)
	{
		pnpPrint("%s",AltNames[i].c_str());
		for( ipH = 0; ipH < NpH; ipH ++ )
		{
			pnpPrint("\t%.10e",Pop.r0(i,ipH));
		}
		pnpPrint("\n");
	}
	return EXIT_SUCCESS;
}
int pKaCalcMod::PrintPop4Homooligamer(bool PrintOnlyIfAltPopMoreThenSmth,float Smth)
{
	PNP_EXIT_FAIL_NULL(phost_mset,"pKaCalcMod does not assign to any MolSet\n");
	pnpPrint("PrintPop4Homooligamer\n");
	int NChains=phost_mset->GetNChains();
	int AltResInChain=NumberOfAltStates/NChains;
	pnpPrint("NChains=%d AltResInChain=%d\n",NChains,AltResInChain);
	
	
	
	int ipH, NpH=roundf((pHmax-pHmin)/pHstep)+1;
	double pH;
	
	int i,iAltSt,iCh;
	AltChemState* alt_st;
	HaResidue* pres;
	for( ipH = 0; ipH < NpH; ipH ++ )
	{
		pH=pHmin+ipH*pHstep;
		pnpPrint("pH=%.4g\n",pH);
		
		pnpPrint("ResNum\tResName\tResMod\tAltSt");
		for(iCh=0; iCh < NChains; iCh++)
		{
			iAltSt=iCh*AltResInChain;
			alt_st=(AltChemState*) act_chem_states[iAltSt];
			pres=(HaResidue*)alt_st->GetHostAtomGroup();
			pnpPrint("\t%c",pres->GetHostChain()->ident);
		}
		pnpPrint("\tPopAvr\tPopSD\tResMore0_45\n");
		
		for(i=0; i < AltResInChain; i++)
		{
			alt_st=(AltChemState*) act_chem_states[i];
			pres=(HaResidue*)alt_st->GetHostAtomGroup();
			
			
			int ResMore0_45=0;
			double PopAvr=0.0,PopSD=0.0;
			double PopMin=Pop.r0(i,ipH);
			double PopMax=Pop.r0(i,ipH);
			for(iCh=0; iCh < NChains; iCh++)
			{
				iAltSt=i+iCh*AltResInChain;
				if(Pop.r0(iAltSt,ipH)>0.45)
					ResMore0_45++;
				PopAvr+=Pop.r0(iAltSt,ipH);
				PopSD+=Pop.r0(iAltSt,ipH)*Pop.r0(iAltSt,ipH);
				if(Pop.r0(iAltSt,ipH)>PopMax)PopMax=Pop.r0(iAltSt,ipH);
				if(Pop.r0(iAltSt,ipH)<PopMin)PopMin=Pop.r0(iAltSt,ipH);
			}
			PopAvr/=double(NChains);
			PopSD=sqrt((PopSD/double(NChains))-PopAvr*PopAvr);
			
			if((!PrintOnlyIfAltPopMoreThenSmth) || (PrintOnlyIfAltPopMoreThenSmth&&PopMax>=Smth))
			{
				pnpPrint("%d\t%s\t%s\t%s", pres->GetSerNo(),  pres->GetName(), pres->GetNameModifier(),alt_st->id.c_str());
				for(iCh=0; iCh < NChains; iCh++)
				{
					iAltSt=i+iCh*AltResInChain;
					pnpPrint("\t%.10e",Pop.r0(iAltSt,ipH));
				}
				pnpPrint("\t%.10e\t%.10e\t%d\n",PopAvr,PopSD,ResMore0_45);
			}
		}
// 		for(i=0; i < NumberOfAltStates; i++)
// 		{
// 			pnpPrint("%s",AltNames[i].c_str());
// 			for( ipH = 0; ipH < NpH; ipH ++ )
// 			{
// 				pnpPrint("\t%.10e",Pop.r0(i,ipH));
// 			}
// 			pnpPrint("\n");
// 			AltChemState* alt_st= pres->GetAltChemState(ist);
// 		}
	}
	return EXIT_SUCCESS;
}
int pKaCalcMod::PrintResults()
{
	PNP_EXIT_FAIL_NULL(phost_mset,"pKaCalcMod does not assign to any MolSet\n");
	pnpPrint("PrintResults\n");
	int old_save_selected=phost_mset->save_opt_default.save_selected;
	phost_mset->save_opt_default.save_selected=FALSE;
	
	
	pnpPrint("E1 = %.14e kT\n",E1);
	int iAltSt=0;
	HaResidue* pres;
	ResidueIteratorMolSet ritr(phost_mset);
	pnpPrint("iAltSt ResNum ResName NumAltSt pKaStd pKaIntr pKa pKaFromPop E2 E3 E4 ddG\n");
	for(pres = ritr.GetFirstRes(); pres; pres = ritr.GetNextRes())
	{
		std::string full_res_name = pres->GetFullName();
		int ist;
		int nst = p_prot_rdx_mod->GetNumResAltChemStates(pres);
		if(nst>0)
			for(ist = 0; ist < nst; ist++)
		{
			AltChemState* alt_st = p_prot_rdx_mod->GetResAltChemState(pres,ist);
			pnpPrint("%d %d %s %d %.3f %.3f %.3f %.3f %.14e %.14e %.14e %.14e %.14e\n",iAltSt,pres->GetSerNo(),full_res_name.c_str(),nst,pKaStd[iAltSt],pKaIntr[iAltSt],pKa[iAltSt],pKaFromPop[iAltSt], E2[iAltSt], E3[iAltSt], E4[iAltSt],ddG[iAltSt]);
			iAltSt++;
		}
			
			
	}
	int ist1, ist2;
	pnpPrint("Interaction Matrix:\n");
	for(ist1=0; ist1 < NumberOfAltStates; ist1++)
	{
		pnpPrint("\t%s",AltNames[ist1].c_str());
	}
	pnpPrint("\n");
	for(ist1 = 0 ; ist1 < NumberOfAltStates ; ist1++)
	{
		pnpPrint("%s",AltNames[ist1].c_str());
		for(ist2 = 0 ; ist2 < NumberOfAltStates ; ist2++)
		{
			pnpPrint("\t%.14e",inter_mat.r0(ist1,ist2));
		}
		pnpPrint("\n");
	}
	
	int ipH, NpH=roundf((pHmax-pHmin)/pHstep)+1;
	double pH;
	
	pnpPrint("pH");
	int i;
	for(i=0; i < NumberOfAltStates; i++)
	{
		pnpPrint("\t%s",AltNames[i].c_str());
	}
	pnpPrint("\n");
	for( ipH = 0; ipH < NpH; ipH ++ )
	{
		pH=pHmin+ipH*pHstep;
		pnpPrint("%.4g",pH);
		for(i=0; i < NumberOfAltStates; i++)
		{
			pnpPrint("\t%.10e",Pop.r0(i,ipH));
		}
		pnpPrint("\n");
	}
	
	phost_mset->save_opt_default.save_selected=old_save_selected;
	return EXIT_SUCCESS;
}
int pKaCalcMod::CalcpKaUsingElectrostMod()
{
	pnpPrint("Number OfAltStates : %d\n",NumberOfAltStates);
	
	int myrank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	
	MakeAltStList();
	RunCalcUsingElectrostMod();
	
	if(myrank==0)
	{
		CalcIntrpKa();
		CalcpKaWithInteraction();
	}
	return EXIT_SUCCESS;
}
int pKaCalcMod::CalcpKaUsingElMod()
{
	pnpPrint("Number OfAltStates : %d\n",NumberOfAltStates);
	
	int myrank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	
	MakeAltStList();
	RunCalcUsingElMod();
	
	if(myrank==0)
	{
		CalcIntrpKa();
		CalcpKaWithInteraction();
	}
	return EXIT_SUCCESS;
}
int pKaCalcMod::ReadCalculatedEnergies(const char *filename)
{
	FILE *in=fopen(filename,"r");
	//char *string=new char[;
	//printf ("Insert your full address: ");
	//gets (string);
	fclose(in);
	return EXIT_SUCCESS;
}
int pKaCalcMod::MakeAltStList()
{
	PNP_EXIT_FAIL_NULL(phost_mset,"pKaCalcMod does not assign to any MolSet\n");
	int old_save_selected=phost_mset->save_opt_default.save_selected;
	phost_mset->save_opt_default.save_selected = TRUE;
	
	AtomIteratorMolSet aitr(phost_mset);
	HaAtom* aptr;
	AltChemState* alt_st;
	HaResidue* rptr;
	ResidueIteratorMolSet ritr(phost_mset);
	
	for(aptr = aitr.GetFirstAtom();aptr; aptr = aitr.GetNextAtom())
	{
		if(aptr->Selected())
			sel_atoms.InsertAtom(aptr);
	}
	
	char CTMP[256];
	for(rptr = ritr.GetFirstRes(); rptr; rptr = ritr.GetNextRes() )
	{
		if(rptr->HasSelectedAtoms() && ( p_prot_rdx_mod->GetNumResAltChemStates(rptr) > 0) )
		{
			int ist;
			int nst = p_prot_rdx_mod->GetNumResAltChemStates(rptr);
			for(ist = 0; ist < nst; ist++)
			{
				AltChemState* alt_st = p_prot_rdx_mod->GetResAltChemState(rptr,ist);
				if(alt_st == NULL) continue;
				if(!alt_st->active_flag) continue;
				
				alt_st->SetAltCharges(0.0);
				act_chem_states.push_back(alt_st);
				if(rptr->GetHostChain()->ident==' ')
					sprintf(CTMP,"%d_%s\0",rptr->GetSerNo(),rptr->GetFullName().c_str());
				else
					sprintf(CTMP,"%d_%s_%c\0",rptr->GetSerNo(),rptr->GetFullName().c_str(),rptr->GetHostChain()->ident);
				string TNAME;
				TNAME+=CTMP;
				TNAME+="_AltState::";
				TNAME+=alt_st->id.c_str();
				/*if(alt_st->alt_state_type ==  AltChemState::PROTONATED)
					TNAME+="_AltChemState::PROTONATED";
				else if(alt_st->alt_state_type ==  AltChemState::UNPROTONATED)
					TNAME+="_AltChemState::UNPROTONATED";
				else if( alt_st->alt_state_type ==  AltChemState::OXIDIZED)
					TNAME+="_AltChemState::OXIDIZED";
				else if( alt_st->alt_state_type ==  AltChemState::REDUCED)
					TNAME+="_AltChemState::REDUCED";*/
				AltNames.push_back(TNAME);
			}
		}
	}
	NumberOfAltStates=act_chem_states.size();
	
	//Init Energy Vectors and interaction energies
	E1=0.0;
	E2.newsize(NumberOfAltStates);
	E2.set(0.0);
	E3.newsize(NumberOfAltStates);
	E3.set(0.0);
	E4.newsize(NumberOfAltStates);
	E4.set(0.0);
	ddG.newsize(NumberOfAltStates);
	ddG.set(0.0);
	pKa.newsize(NumberOfAltStates);
	pKa.set(0.0);
	pKaIntr.newsize(NumberOfAltStates);
	pKaIntr.set(0.0);
	dpKa.newsize(NumberOfAltStates);
	dpKa.set(0.0);
	pKa.newsize(NumberOfAltStates);
	pKa.set(0.0);
	pKaStd.newsize(NumberOfAltStates);
	pKaStd.set(0.0);
	pKaFromPop.newsize(NumberOfAltStates);
	pKaFromPop.set(0.0);
	inter_mat.newsize(NumberOfAltStates,NumberOfAltStates);
	inter_mat.set(0.0);
	
	E1done=0;
	E2done.newsize(NumberOfAltStates);
	E2done.set(0);
	E3done.newsize(NumberOfAltStates);
	E3done.set(0);
	E4done.newsize(NumberOfAltStates);
	E4done.set(0);
	inter_mat_done.newsize(NumberOfAltStates);
	inter_mat_done.set(0);
	
	int i;
	for(i=0; i < NumberOfAltStates; i++)
	{
		alt_st = (AltChemState*) act_chem_states[i];
		
		pKa[i]=alt_st->std_pk;
		pKaStd[i]=alt_st->std_pk;
		pKaIntr[i]=alt_st->std_pk;
	}
	
	phost_mset->save_opt_default.save_selected=old_save_selected;
	return EXIT_SUCCESS;
}
int pKaCalcMod::RunCalcUsingElectrostMod()
{
	PNP_EXIT_FAIL_NULL(phost_mset,"pKaCalcMod does not assign to any MolSet\n");
	int old_save_selected=phost_mset->save_opt_default.save_selected;
	phost_mset->save_opt_default.save_selected = TRUE;
	
	int myrank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	
	bool bres;
	double xmin, xmax, ymin, ymax, zmin, zmax;
	phost_mset->GetMinMaxCrd(xmin, ymin, zmin, xmax, ymax, zmax);
	
	ElectrostMod* m_el=phost_mset->GetElectrostMod(true);
	
	m_el->SetBoundaryAtoms(xmin, ymin, zmin, xmax, ymax, zmax);
	
	//To make sure that only right atoms was selected
	HaAtom* aptr;
	phost_mset->SelectAtoms(&sel_atoms);
	
	int iAltSt=0;
	AltChemState* alt_st;
	
	//E1 all res in their std protonation
	for(iAltSt=0; iAltSt < NumberOfAltStates; iAltSt++)
	{
		alt_st = (AltChemState*) act_chem_states[iAltSt];
		alt_st->SetAltCharges(0.0);
	}
	bres = m_el->run(RUN_FOREGROUND);
	E1 = m_el->tot_ene;
	
	//E2 Res iAltSt in alt state
	for(iAltSt=myrank;iAltSt<NumberOfAltStates;iAltSt=iAltSt+nprocs)
	{
		alt_st = (AltChemState*) act_chem_states[iAltSt];
		
		HaAtom* aptr;
		HaResidue* pres = (HaResidue*) alt_st->GetHostAtomGroup();
		
		map<void*,double,less<void*> > map_aptr_ch;
		AtomIteratorResidue aitr_res(pres);
	
		for(aptr = aitr_res.GetFirstAtom(); aptr; aptr = aitr_res.GetNextAtom())
		{
			map_aptr_ch[(void*)aptr] = aptr->GetCharge();
		}
		
		alt_st->SetAltCharges(1.0);
		
		bres = m_el->run(RUN_FOREGROUND);
		E2[iAltSt] = m_el->tot_ene;
		
		alt_st->SetAltCharges(0.0);
		
		for(aptr = aitr_res.GetFirstAtom(); aptr; aptr = aitr_res.GetNextAtom())
		{
			double ch_at = map_aptr_ch[(void*)aptr];
			aptr->SetCharge(ch_at);
		}
	}
	//E3 - reference system (solvent) unmodified chemical state  (in kT)
	//E4 - reference system (solvent) modified   chemical state  (in kT)
	for(iAltSt=myrank;iAltSt<NumberOfAltStates;iAltSt=iAltSt+nprocs)
	{
		alt_st = (AltChemState*) act_chem_states[iAltSt];
		
		HaAtom* aptr;
		HaResidue* pres = (HaResidue*) alt_st->GetHostAtomGroup();
		
		map<void*,double,less<void*> > map_aptr_ch;
		AtomIteratorResidue aitr_res(pres);
		
		phost_mset->UnSelectAtomsAll();
		for(aptr = aitr_res.GetFirstAtom(); aptr; aptr = aitr_res.GetNextAtom())
		{
			map_aptr_ch[(void*)aptr] = aptr->GetCharge();
			aptr->Select();
		}
		
		alt_st->SetAltCharges(0.0);
				
		bres = m_el->run(RUN_FOREGROUND);
		E3[iAltSt] = m_el->tot_ene;
				
		alt_st->SetAltCharges(1.0);
				
		bres = m_el->run(RUN_FOREGROUND);
		E4[iAltSt] = m_el->tot_ene;
				
		alt_st->SetAltCharges(0.0);
		
		for(aptr = aitr_res.GetFirstAtom(); aptr; aptr = aitr_res.GetNextAtom())
		{
			double ch_at = map_aptr_ch[(void*)aptr];
			aptr->SetCharge(ch_at);
		}
	}
	
	//interaction
	int na = sel_atoms.size();
	HaVec_double old_charges(na);
	
	phost_mset->SelectAtoms(&sel_atoms);
	AtomIteratorAtomGroup aitr_sel(&sel_atoms);
	int i = 0;
	for(aptr = aitr_sel.GetFirstAtom(); aptr; aptr = aitr_sel.GetNextAtom())
	{
		old_charges[i] = aptr->GetCharge();
		i++;
		aptr->SetCharge(0.0);
	}
	
	vector< AtomDoubleMap > sites; // changes of atomic charges during protonation/deprotonation or oxidation/reduction 
	double ch;

	for(i=0; i < NumberOfAltStates; i++)
	{
		alt_st = (AltChemState*) act_chem_states[i];
		HaResidue* rptr = (HaResidue*) alt_st->GetHostAtomGroup();
		
		sites.push_back( AtomDoubleMap("temp"));
		AtomDoubleMap& smap = sites.back();
		
		alt_st->SetAltCharges(-1.0);
		
		AtomIteratorResidue aitr(rptr);
		
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			ch = aptr->GetCharge();
			if(fabs(ch) > 0.00000001)
			{
				smap[aptr] = ch;
			}	
			aptr->SetCharge(0.0);
		}
	}

	int ist1, ist2;
	
	for(ist1 = myrank ; ist1 < NumberOfAltStates ; ist1=ist1+nprocs)
	{
		AtomDoubleMap& smap1 = sites[ist1];
		AtomDoubleMap::iterator mitr;
		for(mitr = smap1.begin(); mitr != smap1.end(); mitr++)
		{
			aptr = (HaAtom*)(*mitr).first;
			ch = (*mitr).second;
			aptr->SetCharge(ch);
		}
		
		  //if(pnp)bres = pnpmod->run(RUN_FOREGROUND);
		bres = m_el->run(RUN_FOREGROUND);
		  
		if(!bres)
		{
			PrintLog("Error in MolSet::CalcPKaForSelection()\n");
			PrintLog("Failed Continuum electrostatics calculations");
		}
		
		for(ist2 = 0; ist2 < NumberOfAltStates; ist2++)
		{ 
			if(ist2 == ist1) continue;
			AtomDoubleMap& smap2 = sites[ist2];
			
			double ene = 0.0;
			for(mitr = smap2.begin(); mitr != smap2.end(); mitr++)
			{
				aptr = (HaAtom*)(*mitr).first;
				ch   = (*mitr).second;
				double phi;
				phi = m_el->el_pot_map.GetInterpolValAtPoint(aptr->GetX(),aptr->GetY(),aptr->GetZ());
				ene += ch*phi;
			}
			inter_mat.r0(ist1,ist2) = ene;
		}
		for(mitr = smap1.begin(); mitr != smap1.end(); mitr++)
		{
			aptr = (HaAtom*)(*mitr).first;
			ch = (*mitr).second;
			aptr->SetCharge(0.0);
		}
	}
	
	i= 0;
	for(aptr = aitr_sel.GetFirstAtom(); aptr; aptr = aitr_sel.GetNextAtom())
	{
		aptr->SetCharge(old_charges[i]);
		i++;
	}
	
	//exenge energies in case of parrallel
	if(nprocs>1)
	{
#if defined(HARLEM_MPI)
		MPI_Barrier(MPI_COMM_WORLD);
		
		double *vE2=E2.begin();
		double *vE3=E3.begin();
		double *vE4=E4.begin();
		
		MPI_Bcast(&E1, 1, MPI_DOUBLE,0, MPI_COMM_WORLD);
		for(iAltSt=0;iAltSt<NumberOfAltStates;iAltSt++)
		{
			MPI_Bcast(&(vE2[iAltSt]), 1, MPI_DOUBLE,iAltSt%nprocs, MPI_COMM_WORLD);
			MPI_Bcast(&(vE3[iAltSt]), 1, MPI_DOUBLE,iAltSt%nprocs, MPI_COMM_WORLD);
			MPI_Bcast(&(vE4[iAltSt]), 1, MPI_DOUBLE,iAltSt%nprocs, MPI_COMM_WORLD);
		}
		double dtmp;
		for(ist1 = 0 ; ist1 < NumberOfAltStates ; ist1++)
		{
			for(ist2 = 0; ist2 < NumberOfAltStates; ist2++)
			{
				dtmp=inter_mat.r0(ist1,ist2);
				MPI_Bcast(&dtmp, 1, MPI_DOUBLE,iAltSt%nprocs, MPI_COMM_WORLD);
				inter_mat.r0(ist1,ist2)=dtmp;
			}
		}
#endif
	}
	
	phost_mset->SelectAtoms(&sel_atoms);
	phost_mset->save_opt_default.save_selected=old_save_selected;
	return EXIT_SUCCESS;
}
int pKaCalcMod::RunCalcUsingElMod()
{
#ifdef ELMOD_COMPILE
	PNP_EXIT_FAIL_NULL(phost_mset,"pKaCalcMod does not assign to any MolSet\n");
	int old_save_selected=phost_mset->p_save_opt_default->save_selected;
	phost_mset->p_save_opt_default->save_selected = TRUE;
	
	int myrank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	
	char filenameTMP[256],filenameTMPOUT[256];
	pnpsapp->AddProcessNumberToFileName(filenameTMPOUT,"TMPPKA_.dat","p_",myrank,nprocs);
	bool bres;
	double xmin, xmax, ymin, ymax, zmin, zmax;
	phost_mset->GetMinMaxCrd(xmin, ymin, zmin, xmax, ymax, zmax);
	
	ElMod* m_el=phost_mset->GetElMod(false);
	
	m_el->SetBoundaryAtoms(xmin, ymin, zmin, xmax, ymax, zmax);
	
	//To make sure that only right atoms was selected
	HaAtom* aptr;
	phost_mset->SelectAtoms(&sel_atoms);
	
	int iAltSt=0;
	AltChemState* alt_st;
	
	//E1 all res in their std protonation
	if(myrank==0&&E1done==0)
	{
		for(iAltSt=0; iAltSt < NumberOfAltStates; iAltSt++)
		{
			alt_st = (AltChemState*) act_chem_states[iAltSt];
			alt_st->SetAltCharges(0.0);
		}
		bres = m_el->run(RUN_FOREGROUND);
		E1 = m_el->tot_ene;
		if(SaveIntermediatePotNNI)
		{
			m_el->m_ContWorld->WriteNodeIndexing("NI_E1.gz");
			m_el->m_ContWorld->WritePotential("Pot_E1.bin");
			if(m_el->boundary==2)//Focusing
			{
				m_el->RoughWorld->m_ContWorld->WriteNodeIndexing("Rough_NI_E1.gz");
				m_el->RoughWorld->m_ContWorld->WritePotential("Rough_Pot_E1.bin");
			}
		}
		if(SaveIntermediateResults)
		{
			FILE *OUTLOC=fopen(filenameTMPOUT,"a");
			fprintf(OUTLOC,"E1 %.14e\n",E1);
			fclose(OUTLOC);
		}
		
	}
	//E2 Res iAltSt in alt state
	for(iAltSt=myrank;iAltSt<NumberOfAltStates;iAltSt=iAltSt+nprocs)
	{
		alt_st = (AltChemState*) act_chem_states[iAltSt];
		
		HaAtom* aptr;
		HaResidue* pres = (HaResidue*) alt_st->GetHostAtomGroup();
		
		map<void*,double,less<void*> > map_aptr_ch;
		AtomIteratorResidue aitr_res(pres);
	
		for(aptr = aitr_res.GetFirstAtom(); aptr; aptr = aitr_res.GetNextAtom())
		{
			map_aptr_ch[(void*)aptr] = aptr->GetCharge();
		}
		
		alt_st->SetAltCharges(1.0);
		if(E2done[iAltSt]==0)
		{
			bres = m_el->run(RUN_FOREGROUND);
			E2[iAltSt] = m_el->tot_ene;
			
			if(SaveIntermediatePotNNI)
			{
				pnpsapp->AddProcessNumberToFileName(filenameTMP,"NI_.gz","E2_",iAltSt,NumberOfAltStates);
				m_el->m_ContWorld->WriteNodeIndexing(filenameTMP);
				pnpsapp->AddProcessNumberToFileName(filenameTMP,"Pot_.bin","E2_",iAltSt,NumberOfAltStates);
				m_el->m_ContWorld->WritePotential(filenameTMP);
				if(m_el->boundary==2)//Focusing
				{
					pnpsapp->AddProcessNumberToFileName(filenameTMP,"Rough_NI_.gz","E2_",iAltSt,NumberOfAltStates);
					m_el->RoughWorld->m_ContWorld->WriteNodeIndexing(filenameTMP);
					pnpsapp->AddProcessNumberToFileName(filenameTMP,"Rough_Pot_.bin","E2_",iAltSt,NumberOfAltStates);
					m_el->RoughWorld->m_ContWorld->WritePotential(filenameTMP);
				}
			}
			if(SaveIntermediateResults)
			{
				FILE *OUTLOC=fopen(filenameTMPOUT,"a");
				fprintf(OUTLOC,"E2 %d %.14e\n",iAltSt,E2[iAltSt]);
				fclose(OUTLOC);
			}
		}
		alt_st->SetAltCharges(0.0);
		
		for(aptr = aitr_res.GetFirstAtom(); aptr; aptr = aitr_res.GetNextAtom())
		{
			double ch_at = map_aptr_ch[(void*)aptr];
			aptr->SetCharge(ch_at);
		}
	}
	//E3 - reference system (solvent) unmodified chemical state  (in kT)
	//E4 - reference system (solvent) modified   chemical state  (in kT)
	int MembraneUsed=m_el->MembraneType;
	m_el->MembraneType=ElMod::MemT_None;
	for(iAltSt=myrank;iAltSt<NumberOfAltStates;iAltSt=iAltSt+nprocs)
	{
		alt_st = (AltChemState*) act_chem_states[iAltSt];
		
		HaAtom* aptr;
		HaResidue* pres = (HaResidue*) alt_st->GetHostAtomGroup();
		
		map<void*,double,less<void*> > map_aptr_ch;
		AtomIteratorResidue aitr_res(pres);
		
		phost_mset->UnSelectAtomsAll();
		for(aptr = aitr_res.GetFirstAtom(); aptr; aptr = aitr_res.GetNextAtom())
		{
			map_aptr_ch[(void*)aptr] = aptr->GetCharge();
			aptr->Select();
		}
		
		alt_st->SetAltCharges(0.0);
		
		if(E3done[iAltSt]==0)
		{
			bres = m_el->run(RUN_FOREGROUND);
			E3[iAltSt] = m_el->tot_ene;
			
			if(SaveIntermediatePotNNI)
			{
				pnpsapp->AddProcessNumberToFileName(filenameTMP,"NI_.gz","E3_",iAltSt,NumberOfAltStates);
				m_el->m_ContWorld->WriteNodeIndexing(filenameTMP);
				pnpsapp->AddProcessNumberToFileName(filenameTMP,"Pot_.bin","E3_",iAltSt,NumberOfAltStates);
				m_el->m_ContWorld->WritePotential(filenameTMP);
				if(m_el->boundary==2)//Focusing
				{
					pnpsapp->AddProcessNumberToFileName(filenameTMP,"Rough_NI_.gz","E3_",iAltSt,NumberOfAltStates);
					m_el->RoughWorld->m_ContWorld->WriteNodeIndexing(filenameTMP);
					pnpsapp->AddProcessNumberToFileName(filenameTMP,"Rough_Pot_.bin","E3_",iAltSt,NumberOfAltStates);
					m_el->RoughWorld->m_ContWorld->WritePotential(filenameTMP);
				}
			}
			if(SaveIntermediateResults)
			{
				FILE *OUTLOC=fopen(filenameTMPOUT,"a");
				fprintf(OUTLOC,"E3 %d %.14e\n",iAltSt,E3[iAltSt]);
				fclose(OUTLOC);
			}
		}
		alt_st->SetAltCharges(1.0);
		
		if(E4done[iAltSt]==0)
		{
			bres = m_el->run(RUN_FOREGROUND);
			E4[iAltSt] = m_el->tot_ene;
			
			if(SaveIntermediatePotNNI)
			{
				pnpsapp->AddProcessNumberToFileName(filenameTMP,"NI_.gz","E4_",iAltSt,NumberOfAltStates);
				m_el->m_ContWorld->WriteNodeIndexing(filenameTMP);
				pnpsapp->AddProcessNumberToFileName(filenameTMP,"Pot_.bin","E4_",iAltSt,NumberOfAltStates);
				m_el->m_ContWorld->WritePotential(filenameTMP);
				if(m_el->boundary==2)//Focusing
				{
					pnpsapp->AddProcessNumberToFileName(filenameTMP,"Rough_NI_.gz","E4_",iAltSt,NumberOfAltStates);
					m_el->RoughWorld->m_ContWorld->WriteNodeIndexing(filenameTMP);
					pnpsapp->AddProcessNumberToFileName(filenameTMP,"Rough_Pot_.bin","E4_",iAltSt,NumberOfAltStates);
					m_el->RoughWorld->m_ContWorld->WritePotential(filenameTMP);
				}
			}
			if(SaveIntermediateResults)
			{
				FILE *OUTLOC=fopen(filenameTMPOUT,"a");
				fprintf(OUTLOC,"E4 %d %.14e\n",iAltSt,E4[iAltSt]);
				fclose(OUTLOC);
			}
		}
		alt_st->SetAltCharges(0.0);
		
		for(aptr = aitr_res.GetFirstAtom(); aptr; aptr = aitr_res.GetNextAtom())
		{
			double ch_at = map_aptr_ch[(void*)aptr];
			aptr->SetCharge(ch_at);
		}
	}
	m_el->MembraneType=MembraneUsed;
	//interaction
	int na = sel_atoms.size();
	HaVec_double old_charges(na);
	
	phost_mset->SelectAtoms(&sel_atoms);
	AtomIteratorAtomGroup aitr_sel(&sel_atoms);
	int i = 0;
	for(aptr = aitr_sel.GetFirstAtom(); aptr; aptr = aitr_sel.GetNextAtom())
	{
		old_charges[i] = aptr->GetCharge();
		i++;
		aptr->SetCharge(0.0);
	}
	
	vector< AtomDoubleMap > sites; // changes of atomic charges during protonation/deprotonation or oxidation/reduction 
	double ch;

	for(i=0; i < NumberOfAltStates; i++)
	{
		alt_st = (AltChemState*) act_chem_states[i];
		HaResidue* rptr = (HaResidue*) alt_st->GetHostAtomGroup();
		
		sites.push_back( AtomDoubleMap("temp"));
		AtomDoubleMap& smap = sites.back();
		
		alt_st->SetAltCharges(-1.0);
		
		AtomIteratorResidue aitr(rptr);
		
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			ch = aptr->GetCharge();
			if(fabs(ch) > 0.00000001)
			{
				smap[aptr] = ch;
			}	
			aptr->SetCharge(0.0);
		}
	}

	int ist1, ist2;
	
	for(ist1 = myrank ; ist1 < NumberOfAltStates ; ist1=ist1+nprocs)
	{
		AtomDoubleMap& smap1 = sites[ist1];
		AtomDoubleMap::iterator mitr;
		for(mitr = smap1.begin(); mitr != smap1.end(); mitr++)
		{
			aptr = (HaAtom*)(*mitr).first;
			ch = (*mitr).second;
			aptr->SetCharge(ch);
		}
		
		  //if(pnp)bres = pnpmod->run(RUN_FOREGROUND);
		
		if(inter_mat_done[ist1]==0)
		{
			bres = m_el->run(RUN_FOREGROUND);
			
			if(SaveIntermediatePotNNI)
			{
				pnpsapp->AddProcessNumberToFileName(filenameTMP,"NI_.gz","Iter_",iAltSt,NumberOfAltStates);
				m_el->m_ContWorld->WriteNodeIndexing(filenameTMP);
				pnpsapp->AddProcessNumberToFileName(filenameTMP,"Pot_.bin","Iter_",iAltSt,NumberOfAltStates);
				m_el->m_ContWorld->WritePotential(filenameTMP);
				if(m_el->boundary==2)//Focusing
				{
					pnpsapp->AddProcessNumberToFileName(filenameTMP,"Rough_NI_.gz","Iter_",iAltSt,NumberOfAltStates);
					m_el->RoughWorld->m_ContWorld->WriteNodeIndexing(filenameTMP);
					pnpsapp->AddProcessNumberToFileName(filenameTMP,"Rough_Pot_.bin","Iter_",iAltSt,NumberOfAltStates);
					m_el->RoughWorld->m_ContWorld->WritePotential(filenameTMP);
				}
			}
			
			if(!bres)
			{
				PrintLog("Error in MolSet::CalcPKaForSelection()\n");
				PrintLog("Failed Continuum electrostatics calculations");
			}
			
			for(ist2 = 0; ist2 < NumberOfAltStates; ist2++)
			{ 
				if(ist2 == ist1) continue;
				AtomDoubleMap& smap2 = sites[ist2];
				
				double ene = 0.0;
				for(mitr = smap2.begin(); mitr != smap2.end(); mitr++)
				{
					aptr = (HaAtom*)(*mitr).first;
					ch   = (*mitr).second;
					double phi;
					phi = m_el->el_pot_map.GetInterpolValAtPoint(aptr->GetX(),aptr->GetY(),aptr->GetZ());
					ene += ch*phi;
				}
				inter_mat.r0(ist1,ist2) = ene;
			}
			if(SaveIntermediateResults)
			{
				FILE *OUTLOC=fopen(filenameTMPOUT,"a");
				fprintf(OUTLOC,"In %d",ist1);
				for(ist2 = 0; ist2 < NumberOfAltStates; ist2++)
				{
					fprintf(OUTLOC," %.14e",inter_mat.r0(ist1,ist2));
				}
				fprintf(OUTLOC,"\n",ist1);
				fclose(OUTLOC);
			}
		}
		for(mitr = smap1.begin(); mitr != smap1.end(); mitr++)
		{
			aptr = (HaAtom*)(*mitr).first;
			ch = (*mitr).second;
			aptr->SetCharge(0.0);
		}
	}
	
	i= 0;
	for(aptr = aitr_sel.GetFirstAtom(); aptr; aptr = aitr_sel.GetNextAtom())
	{
		aptr->SetCharge(old_charges[i]);
		i++;
	}
	
	//exenge energies in case of parrallel
	if(nprocs>1)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		
		double *vE2=E2.begin();
		double *vE3=E3.begin();
		double *vE4=E4.begin();
		
		MPI_Bcast(&E1, 1, MPI_DOUBLE,0, MPI_COMM_WORLD);
		for(iAltSt=0;iAltSt<NumberOfAltStates;iAltSt++)
		{
			MPI_Bcast(&(vE2[iAltSt]), 1, MPI_DOUBLE,iAltSt%nprocs, MPI_COMM_WORLD);
			MPI_Bcast(&(vE3[iAltSt]), 1, MPI_DOUBLE,iAltSt%nprocs, MPI_COMM_WORLD);
			MPI_Bcast(&(vE4[iAltSt]), 1, MPI_DOUBLE,iAltSt%nprocs, MPI_COMM_WORLD);
		}
		double dtmp;
		for(ist1 = 0 ; ist1 < NumberOfAltStates ; ist1++)
		{
			for(ist2 = 0; ist2 < NumberOfAltStates; ist2++)
			{
				dtmp=inter_mat.r0(ist1,ist2);
				MPI_Bcast(&dtmp, 1, MPI_DOUBLE,ist1%nprocs, MPI_COMM_WORLD);
				inter_mat.r0(ist1,ist2)=dtmp;
			}
		}
	}
	
	phost_mset->SelectAtoms(&sel_atoms);
	phost_mset->p_save_opt_default->save_selected = old_save_selected;
#endif
	return EXIT_SUCCESS;
}
int pKaCalcMod::CalcIntrpKa()
{
	PNP_EXIT_FAIL_NULL(phost_mset,"pKaCalcMod does not assign to any MolSet\n");
	int old_save_selected=phost_mset->save_opt_default.save_selected;
	phost_mset->save_opt_default.save_selected = TRUE;
	
	
	int iAltSt=0;
	bool bres;
	phost_mset->SelectAtomsAll();
	HaResidue* pres;
	ResidueIteratorMolSet ritr(phost_mset);
	
	for(pres = ritr.GetFirstRes(); pres; pres = ritr.GetNextRes())
	{
		int ist;
		int nst = p_prot_rdx_mod->GetNumResAltChemStates(pres);
		if(nst>0)
			for(ist = 0; ist < nst; ist++)
		{
			AltChemState* alt_st = p_prot_rdx_mod->GetResAltChemState(pres,ist);
			ddG[iAltSt]=(E2[iAltSt]-E1)-(E4[iAltSt]-E3[iAltSt]);
			if(alt_st->alt_state_type == alt_st->alt_state_type.PROTONATED)
			{
				dpKa[iAltSt] = -ddG[iAltSt]/log(10.0);
				pKaIntr[iAltSt]=pKaStd[iAltSt]+dpKa[iAltSt];
				pKa[iAltSt]=pKaIntr[iAltSt];
			}
			else if(alt_st->alt_state_type == alt_st->alt_state_type.UNPROTONATED)
			{
				dpKa[iAltSt] = ddG[iAltSt]/log(10.0);
				pKaIntr[iAltSt]=pKaStd[iAltSt]+dpKa[iAltSt];
				pKa[iAltSt]=pKaIntr[iAltSt];
			}
			iAltSt++;
		}
	}
	
	
	phost_mset->save_opt_default.save_selected=old_save_selected;
	return EXIT_SUCCESS;
}
int pKaCalcMod::RoutineCalcAvgPopMC(HaMat_double& inter_mat, HaVec_double& alt_pop, HaVec_double& delt_e,int N_mc_cyc)
//!
//! interaction matrix - on the diagonal energy difference of alternative and unmodifed state when 
//! alt_pop[j] - population of aternative states  
//! delt_e[j]  - converged self-consistent difference between alternative and unmodifed state 
//!
{
	int nst = inter_mat.num_rows();

	delt_e.newsize(nst);
	delt_e = 0.0;
	alt_pop.newsize(nst); // population of alternative states
	alt_pop = 0.0;

	HaVec_double state_vec(nst);

	int i;
	for(i= 0; i < nst; i++)
	{
		if( inter_mat.r0(i,i) < 0.0)
		{
			state_vec[i] = 1.0;
		}
		else
		{
			state_vec[i] = 0.0;
		}
	}

	int icyc;
	int nstep = 0;

	Random rand_num_gen(5);
	
	for( icyc = 0; icyc < N_mc_cyc; icyc++)
	{
		if(icyc%10000==0)
			pnpPrint("RoutineCalcAvgPopMC iter:%d\n",icyc);
		for( i = 0; i < nst; i++)
		{
			int accept = TRUE;
			double delt_e = dot_double(state_vec.begin(),&inter_mat.r0(0,i), nst);
			if( state_vec[i] > 0.5 )
			{
				delt_e = -delt_e;
			}
			else
			{
				delt_e += inter_mat.r0(i,i);
			}
			
			if( delt_e > 0)
			{
				double dt = exp(-delt_e);
				double rn = rand_num_gen();
				if( dt < rn) accept = FALSE;
			}
			if( accept )
			{
				nstep++;
				if( state_vec[i] > 0.5) state_vec[i] = 0.0;
				else state_vec[i] = 1.0;
				int k;
				for(k = 0; k < nst; k++)
				{
					alt_pop[k] += state_vec[k];
				}
			}
		}
	}
	
	if(nstep > 0)
	{
		for(i = 0; i < nst; i++)
		{
			alt_pop[i] = alt_pop[i]/nstep;
			if( (1.0 - alt_pop[i]) < DBL_EPSILON  ) 
			{
				delt_e[i] = -100.0;
			}
			else if( alt_pop[i] < DBL_EPSILON ) 
			{
				delt_e[i] = +100.0;
			}
			else
			{
				delt_e[i] = log(1.0/alt_pop[i] - 1.0);
			}
		}
	}
	return EXIT_SUCCESS;
}
int pKaCalcMod::RoutineCalcAvgPopSCF(HaMat_double& inter_mat, HaVec_double& alt_pop, HaVec_double& delt_e,int max_iter,double pop_err_max)
{
	int nst = inter_mat.num_rows();

	delt_e.newsize(nst);
	delt_e = 0.0;
	alt_pop.newsize(nst); // population of alternative states
	alt_pop = 0.0;

	int iter = 0;
	int i,j;

	int converged = FALSE;

	for(iter = 0 ; iter < max_iter; iter++)
	{
//		PrintLog("Iter %4d :\n",iter); 
//		PrintLog("Popul    :\n"); 
		
//		PrintLog("pKa      :\n"); 
		double pop_err = 0.0;
		for(i=0; i < nst; i++)
		{
			delt_e[i] = inter_mat.r0(i,i);
			for(j = 0; j < nst; j++)
			{
				if(j == i) continue;
				delt_e[i] += inter_mat.r0(i,j) * alt_pop[j];
			}
		}

		for(i=0; i<nst;i++)
		{
			double pop_old = alt_pop[i];
			double pop_new = 1.0/(1.0 + exp(delt_e[i]));
//			PrintLog(" %d %8.3f",i,alt_pop[i]);
			if( iter > 0)
			{
				alt_pop[i] = pop_old*0.6 + pop_new*0.4;
			}
			else
			{
				alt_pop[i] = pop_new;
			}
			if( fabs(pop_new - pop_old) > pop_err) pop_err = fabs(pop_new - pop_old);
		}
//		PrintLog("\n");         

//		PrintLog(" Max population change: %9.5f \n",pop_err);

		if(pop_err < pop_err_max)
		{
//			PrintLog(" Population values have converged \n");
			converged = TRUE;
			break;
		}
	}
	
	if(!converged)
	{
		PrintLog(" Maximum Number of population convergence iteration exceeded");
	}

	return EXIT_SUCCESS;
}
int pKaCalcMod::CalcpKaWithInteraction()
{
	int i;
	AltChemState* alt_st;
	int ipH, NpH=roundf((pHmax-pHmin)/pHstep)+1;
	double pH;
	
	HaVec_double alt_pop(NumberOfAltStates), delt_e(NumberOfAltStates);
	Pop.newsize(NumberOfAltStates,NpH);
	Pop.set(0.0);
	
	PrintLog("CalcpKaWithInteraction\n");
	if(pKaCalcMethod == SCF_MULTI_SITE_CALC)
	{
		PrintLog("\tMethod: SCF_MULTI_SITE_CALC\n");
	}
	else
	{
		PrintLog("\tMethod: MC_MULTI_SITE_CALC\n");
	}
	for( ipH = 0; ipH < NpH; ipH ++ )
	{
		pH=pHmin+ipH*pHstep;
		pnpPrint("pH=%.3g\n",pH);
		for(i=0; i < NumberOfAltStates; i++)
		{
			double de;
			alt_st = (AltChemState*) act_chem_states[i];
			if(alt_st->alt_state_type ==  alt_st->alt_state_type.PROTONATED)
			{
				de = log(10.0)*(pH - pKaIntr[i]); // Diff of energies of alternative and reference states in kT
			}
			else if(alt_st->alt_state_type ==  alt_st->alt_state_type.UNPROTONATED)
			{
				de = log(10.0)*(pKaIntr[i] - pH);
			}
			inter_mat.r0(i,i) = de;
		}
		
		if(pKaCalcMethod == SCF_MULTI_SITE_CALC)
		{
			RoutineCalcAvgPopSCF(inter_mat, alt_pop, delt_e,SCF_PKA_CALC_max_iter,SCF_PKA_CALC_pop_err_max);
		}
		else
		{
			RoutineCalcAvgPopMC(inter_mat, alt_pop, delt_e,MC_PKA_CALC_N_mc_cyc);
		}
		
		//copy population to Pop
		for(i=0; i < NumberOfAltStates; i++)
		{
			Pop.r0(i,ipH)=alt_pop[i];
		}
		
		if(ipH==NpH/2)
		{
			PrintLog(" Energy shifts of alternative states energies at pH %f \n",pH);
			for(i=0; i < NumberOfAltStates; i++)
			{
				PrintLog(" %4d self_ene= %12.6f (kT)  tot_shift = %12.6f (kT) \n",
								i,inter_mat.r0(i,i), delt_e[i]);
				alt_st = (AltChemState*) act_chem_states[i];
				if(alt_st->alt_state_type ==  alt_st->alt_state_type.PROTONATED)
				{
					pKa[i] = pH - delt_e[i]/log(10.0);
				}
				else if(alt_st->alt_state_type ==  alt_st->alt_state_type.UNPROTONATED)
				{
					pKa[i] = pH + delt_e[i]/log(10.0); 
				}
			}
		}
	}
	//pKa based on population
	for(i=0; i < NumberOfAltStates; i++)
	{
		alt_st = (AltChemState*) act_chem_states[i];
		if(alt_st->alt_state_type ==  alt_st->alt_state_type.PROTONATED)
		{
			if(Pop.r0(i,0)<0.5)
				pKaFromPop[i]=-10.0;
			else if(Pop.r0(i,NpH-1)>0.5)
				pKaFromPop[i]=30.0;
			else
			{
				for( ipH = 0; ipH < NpH-1; ipH ++ )
				{
					if(Pop.r0(i,ipH)>0.5&&Pop.r0(i,ipH+1)<=0.5)
					{
						pH=pHmin+ipH*pHstep;
						pKaFromPop[i]=(0.5-Pop.r0(i,ipH+1))*pHstep/(Pop.r0(i,ipH)-Pop.r0(i,ipH+1))+pH;
					}
				}
			}
		}
		else if(alt_st->alt_state_type ==  alt_st->alt_state_type.UNPROTONATED)
		{
			if(Pop.r0(i,0)>0.5)
				pKaFromPop[i]=-10.0;
			else if(Pop.r0(i,NpH-1)<0.5)
				pKaFromPop[i]=30.0;
			else
			{
				for( ipH = 0; ipH < NpH-1; ipH ++ )
				{
					if(Pop.r0(i,ipH)<0.5&&Pop.r0(i,ipH+1)>=0.5)
					{
						pH=pHmin+ipH*pHstep;
						pKaFromPop[i]=(0.5-Pop.r0(i,ipH+1))*pHstep/(Pop.r0(i,ipH)-Pop.r0(i,ipH+1))+pH;
					}
				}
			}
		}
	}
	return EXIT_SUCCESS;
}
int pKaCalcMod::WritepKaCalcModToFile(const char *filename)
{
	TiXmlDocument Doc(filename);
	TiXmlElement Elt("pKaCalcMod");
	
	WritepKaCalcModToXmlElement(&Elt);
	
	Doc.InsertEndChild(Elt);
	Doc.SaveFile();
	return EXIT_SUCCESS;
}
int pKaCalcMod::ReadpKaCalcModFromFile(const char *filename)
{
	TiXmlDocument Doc(filename);
	Doc.LoadFile();
	
	TiXmlElement *Elt=Doc.RootElement();
	
	ReadpKaCalcModFromXmlElement(Elt);
	
	return EXIT_SUCCESS;
}
int pKaCalcMod::WritepKaCalcModToXmlElement(TiXmlElement *RootElt)
{
	RootElt->SetIntAttribute("NumberOfAltStates",NumberOfAltStates);
	
	RootElt->SetDoubleAttribute("pHmin",pHmin);
	RootElt->SetDoubleAttribute("pHmax",pHmax);
	RootElt->SetDoubleAttribute("pHstep",pHstep);
	
	RootElt->SetStdStrIndex("pKaCalcMethod",pKaCalcMethod,pKaCalcMethodStr);
	
	RootElt->SetIntAttribute("MC_PKA_CALC_N_mc_cyc",MC_PKA_CALC_N_mc_cyc);
	RootElt->SetIntAttribute("SCF_PKA_CALC_max_iter",SCF_PKA_CALC_max_iter);
	RootElt->SetDoubleAttribute("SCF_PKA_CALC_pop_err_max",SCF_PKA_CALC_pop_err_max);
	
	RootElt->SetDoubleAttribute("E1",E1);
	
	RootElt->SetVectorStringElement("AltNames",&AltNames);
	
	RootElt->SetHaVec_doubleElement("E2",&E2,"%.14e");
	RootElt->SetHaVec_doubleElement("E3",&E3,"%.14e");
	RootElt->SetHaVec_doubleElement("E4",&E4,"%.14e");
	RootElt->SetHaVec_doubleElement("ddG",&ddG,"%.14e");
	
	RootElt->SetHaVec_doubleElement("pKaStd",&pKaStd,"%.14e");
	RootElt->SetHaVec_doubleElement("pKaIntr",&pKaIntr,"%.14e");
	RootElt->SetHaVec_doubleElement("pKaFromPop",&pKaFromPop,"%.14e");
	RootElt->SetHaVec_doubleElement("pKa",&pKa,"%.14e");
	RootElt->SetHaVec_doubleElement("dpKa",&dpKa,"%.14e");
	
	RootElt->SetHaMat_doubleElement("inter_mat",&inter_mat);
	
	RootElt->SetHaMat_doubleElement("Pop",&Pop);
	
	return EXIT_SUCCESS;
}
int pKaCalcMod::ReadpKaCalcModFromXmlElement(const TiXmlElement *RootElt)
{
	RootElt->GetIntAttribute("NumberOfAltStates",&NumberOfAltStates);
	
	RootElt->GetDoubleAttribute("pHmin",&pHmin);
	RootElt->GetDoubleAttribute("pHmax",&pHmax);
	RootElt->GetDoubleAttribute("pHstep",&pHstep);
	
	RootElt->GetStdStrIndex("pKaCalcMethod",&pKaCalcMethod,pKaCalcMethodStr);
	
	RootElt->GetIntAttribute("MC_PKA_CALC_N_mc_cyc",&MC_PKA_CALC_N_mc_cyc);
	RootElt->GetIntAttribute("SCF_PKA_CALC_max_iter",&SCF_PKA_CALC_max_iter);
	RootElt->GetDoubleAttribute("SCF_PKA_CALC_pop_err_max",&SCF_PKA_CALC_pop_err_max);
	
	RootElt->GetDoubleAttribute("E1",&E1);
	
	RootElt->GetVectorStringElement("AltNames",&AltNames);
	
	RootElt->GetHaVec_doubleElement("E2",&E2);
	RootElt->GetHaVec_doubleElement("E3",&E3);
	RootElt->GetHaVec_doubleElement("E4",&E4);
	RootElt->GetHaVec_doubleElement("ddG",&ddG);
	
	RootElt->GetHaVec_doubleElement("pKaStd",&pKaStd);
	RootElt->GetHaVec_doubleElement("pKaIntr",&pKaIntr);
	RootElt->GetHaVec_doubleElement("pKaFromPop",&pKaFromPop);
	RootElt->GetHaVec_doubleElement("pKa",&pKa);
	RootElt->GetHaVec_doubleElement("dpKa",&dpKa);
	
	RootElt->GetHaMat_doubleElement("inter_mat",&inter_mat);
	
	RootElt->GetHaMat_doubleElement("Pop",&Pop);
	
	return EXIT_SUCCESS;
}
int pKaCalcMod::SetAltSt4ResInHomoolgmr(int ResNum)
{
	int i;
	
	pnpPrint("Electrostatic charge and dipole befor changing:\n");
	phost_mset->CalcDipole();
	for(i=0; i < NumberOfAltStates; i++)
	{
		AltChemState* alt_st=(AltChemState*) act_chem_states[i];
		HaResidue*pres=(HaResidue*)alt_st->GetHostAtomGroup();
		if(pres->GetSerNo()==ResNum)
		{
			pnpPrint("Set alternative state for %d\t%s\t%s\t%s\n", pres->GetSerNo(),  pres->GetName(), pres->GetNameModifier(),alt_st->id.c_str());
			alt_st->SetAltCharges(1.0);
		}
	}
	pnpPrint("Electrostatic charge and dipole after changing:\n");
	phost_mset->CalcDipole();
	return EXIT_SUCCESS;
}
////////////////////////////////////////////////////////////////////////////
PNPMod::PNPMod(MolSet* new_phost_mset)
	:HaCompMod(COMP_MOD_PNP,new_phost_mset)
{
#ifdef PNP_DEPRECATED
	PNPSApp::InitPNPSApp();
#endif
	SetNIonsTypes(2);
	SetIonName(0,"K");
	SetIonName(1,"Cl");
	
	LJA=NULL;
	LJB=NULL;
}

PNPMod::~PNPMod()
{
}

int PNPMod::SetNIonsTypes(int m_NIonsTypes)
{
	NIonsTypes=m_NIonsTypes;
	IonNames.resize(NIonsTypes);
	return EXIT_SUCCESS;
}
int PNPMod::SetIonName(int ion,const char * name)
{
	if(ion<NIonsTypes)
	{
		IonNames[ion]=name;
		return EXIT_SUCCESS;
	}
	else
	{
		PrintLog("Error in PNPMod::SetIonName() ion = % d > NIonsTypes = % d\n", ion, NIonsTypes);
		return EXIT_FAILURE;
	}
}
const char * PNPMod::GetIonName(int ion)
{
	if(ion<NIonsTypes)
	{
		return IonNames[ion].c_str();
	}
	else
	{
		PrintLog("Error in PNPMod::GetIonName() ion = % d > NIonsTypes = % d\n", ion, NIonsTypes);
		return NULL;
	}
}
int PNPMod::ReadAMBER94FF()
{
	boost::filesystem::path amber_94_ff_path(pApp->res_db_dir);
	amber_94_ff_path /= "amber_94_ff.dat";
	int status = ReadAMBERFF(amber_94_ff_path.string().c_str());
	
	//fill synonyms
	int i;
	int AtomTypeNumber;
	
	AtomTypeNumber=GetAtomTypeNumber(std::string("N"));
	char Nsyn[13][4]={"NA","N2","N*","NC","NB","N3","NP","NO","N4","N5","N6","N7","N8"};
	for(i=0; i < 13; i++) 
	{
		AtomTypesDB.push_back(std::string(Nsyn[i]));
		HalfSigmaDB.push_back(HalfSigmaDB[AtomTypeNumber]);
		FourEpsilonDB.push_back(FourEpsilonDB[AtomTypeNumber]);
	}
	
	AtomTypeNumber=GetAtomTypeNumber(std::string("C"));
	char Csyn[15][4]={"C*","CA","CB","CC","CN","CM","CK","CQ","CW","CV","CR","CA","CX","CY","CD"};
	for(i=0; i < 15; i++) 
	{
		AtomTypesDB.push_back(std::string(Csyn[i]));
		HalfSigmaDB.push_back(HalfSigmaDB[AtomTypeNumber]);
		FourEpsilonDB.push_back(FourEpsilonDB[AtomTypeNumber]);
	}
	
	
	for(i=0; i < AtomTypesDB.size(); i++) 
	{
		std::cout << AtomTypesDB[i] << " "  << HalfSigmaDB[i] << " "  << FourEpsilonDB[i] << " "  << endl;
	}
	return status;
}
int PNPMod::ReadAMBERFF(const char* filename)
{
	AtomTypesDB.clear();
	HalfSigmaDB.clear();
	FourEpsilonDB.clear();
	//Read Amber 94 FF
	std::string amber_94_ff_file(filename);
	
	char strline[1024]="000000000";
	
	FILE *fin=fopen(amber_94_ff_file.c_str(),"rt");
	if(fin==NULL)
	{
		std::string ErrorMsg="Can not open "+amber_94_ff_file;
		ErrorInMod("PNPMod::SavePREFile", ErrorMsg.c_str());
		return EXIT_FAILURE;
	}
	
	//find section with sigmas
	while ((!feof(fin))&&(!(strline[0]=='M'&&strline[1]=='O'&&strline[2]=='D'&&strline[3]=='4'))) 
	{
		fgets(strline,1023,fin);
	}
	if(!(strline[0]=='M'&&strline[1]=='O'&&strline[2]=='D'&&strline[3]=='4'))
	{
		std::string ErrorMsg="Can not find section with VdW params at "+amber_94_ff_file;
		ErrorInMod("PNPMod::SavePREFile", ErrorMsg.c_str());
		fclose(fin);
		return EXIT_FAILURE;
	}
	//read database
	char atomname[5];
	double r,eps;
	
	while ((!feof(fin))&&(!(strline[0]=='E'&&strline[1]=='N'&&strline[2]=='D'))) 
	{
		fgets(strline,1023,fin);
		if(sscanf(strline,"%s %lf %lf",atomname,&r,&eps)==3)
		{
			AtomTypesDB.push_back(std::string(atomname));
			HalfSigmaDB.push_back(r);
			FourEpsilonDB.push_back(eps/HaConsts::kT_to_kcal_mol);
		}
			//PrintLog(strline);
	}
	/*int i;
	for(i=0; i < AtomTypesDB.size(); i++) 
	{
		std::cout << AtomTypesDB[i] << " "  << HalfSigmaDB[i] << " "  << FourEpsilonDB[i] << " "  << endl;
	}*/
	return EXIT_SUCCESS;
}
int PNPMod::ReadOPLSFF()
{
	boost::filesystem::path ffoplsaanb_itp_path = pApp->res_db_dir;
	ffoplsaanb_itp_path /= "ffoplsaanb.itp";
	boost::filesystem::path ffoplsaa_rtp_path = pApp->res_db_dir;
	ffoplsaa_rtp_path /= "ffoplsaa.rtp";
	
	int status;
	status = ReadOPLSitp(ffoplsaanb_itp_path.string().c_str());
	status = ReadOPLSrtp(ffoplsaa_rtp_path.string().c_str());
	return EXIT_SUCCESS;
}
int PNPMod::ReadOPLSitp(const char *filename)
{
	FILE *fin=fopen(filename,"rt");
	if(fin==NULL)
	{
		std::string ErrorMsg="Can not open "+ std::string(filename);
		ErrorInMod("PNPMod::ReadOPLSitp", ErrorMsg.c_str());
		return EXIT_FAILURE;
	}
	int i;
	char strline[1024];
	//find [ atomtypes ]
	while (!feof(fin))
	{
		fgets(strline,1023,fin);
		if(strline[0]==';')continue;
		
		std::string wstr(strline);
		boost::algorithm::trim(wstr);
		if(wstr[0] == '[')
		{
			i=wstr.find(";");
			if(i != string::npos) wstr=wstr.substr(i);
                                boost::replace_all(wstr," ","");
                                boost::replace_all(wstr,"[","");
                                boost::replace_all(wstr,"]","");
			// wstr = std::regex_replace(wstr, std::regex(" "), "");
			// wstr = std::regex_replace(wstr, std::regex("["), "");
			// wstr = std::regex_replace(wstr, std::regex("]"), "");
			boost::algorithm::trim(wstr);
			if(boost::iequals(wstr,"atomtypes"))
			{
				PrintLog("found %s\n",wstr.c_str());
				break;
			}
			
		}
	}
	//read params
	while (!feof(fin))
	{
		fgets(strline,1023,fin);
		if(strline[0]==';')continue;
		if(strline[0]=='[')break;
		
		char atmtype[12];
		char name[5];
		int bond_type;
		double mass;
		double charge;
		char ptype;
		double sigma;
		double epsilon;
		if(sscanf(strline,"%s %s %d %lg %lg %c %lg %lg", atmtype, name, &bond_type, &mass, &charge, &ptype, &sigma, &epsilon)==8)
		{
			//PrintLog("found %s %s %d %g %g %c %g %g\n",atmtype, name, bond_type, mass, charge, ptype, sigma, epsilon);
			OPLSLJAtomTypes.push_back(atmtype);
			OPLSLJAtomName.push_back(name);
			OPLSLJSigma.push_back(sigma*10.0);
			OPLSLJEpsilon.push_back(epsilon*HaConsts::kJ_mol_to_kT);
		}
	}
	fclose(fin);
	return EXIT_SUCCESS;
}
int PNPMod::PrintOPLSLJSigmaEpsilon()
{
	int iat;
	PrintLog("AtomName AtomTypes Sigma[A] Epsilon[kT]\n");
	for(iat=0;iat<OPLSLJAtomName.size();iat++)
	{
		PrintLog("%10s %10s %14.10e %14.10e\n",OPLSLJAtomName[iat].c_str(),OPLSLJAtomTypes[iat].c_str(),OPLSLJSigma[iat],OPLSLJEpsilon[iat]);
	}
	return EXIT_SUCCESS;
}
int PNPMod::PrintOPLSLJAB()
{
	int iat;
	PrintLog("AtomName AtomTypes A[kT*A^12] B[kT*A^6]\n");
	for(iat=0;iat<OPLSLJAtomName.size();iat++)
	{
		double A,B;
		double s2=OPLSLJSigma[iat]*OPLSLJSigma[iat];
		double s6=s2*s2*s2;
		
		A = 4.0*OPLSLJEpsilon[iat]*s6*s6;
		B = 4.0*OPLSLJEpsilon[iat]*s6;
		PrintLog("%10s %10s %14.10e %14.10e\n",OPLSLJAtomName[iat].c_str(),OPLSLJAtomTypes[iat].c_str(),A,B);
	}
	return EXIT_SUCCESS;
}
int PNPMod::ReadOPLSrtp(const char *filename)
{
	FILE *fin=fopen(filename,"rt");
	if(fin==NULL)
	{
		std::string ErrorMsg="Can not open "+std::string(filename);
		ErrorInMod("PNPMod::ReadOPLSrtp", ErrorMsg.c_str());
		return EXIT_FAILURE;
	}
	int i;
	char strline[1024];
	fgets(strline,1023,fin);
	while (!feof(fin))
	{
		//find residue
		bool foundres=false;
		while ((!feof(fin))&&(!foundres))
		{
			std::string wstr(strline);
			boost::algorithm::trim(wstr);
			if(wstr[0] == '[')
			{
				i=wstr.find(";");
				if(i != std::string::npos) wstr=wstr.substr(i);
                                boost::replace_all(wstr," ","");
                                boost::replace_all(wstr,"[","");
                                boost::replace_all(wstr,"]","");
				// wstr = std::regex_replace(wstr, std::regex(" "), "");
				// wstr = std::regex_replace(wstr, std::regex("["), "");
				// wstr = std::regex_replace(wstr, std::regex("]"), "");

				if( wstr == boost::to_upper_copy(wstr) )
				{
					//here is resedue
					boost::algorithm::trim(wstr);
					PrintLog("found residue %s\n",wstr.c_str());
					foundres=true;
					OPLSResNames.push_back(wstr);
				}
			}
			fgets(strline,1023,fin);
		}
		//find atoms section
		bool foundatms=false;
		while ((!feof(fin))&&(!foundatms))
		{
			std::string wstr(strline);
			boost::algorithm::trim(wstr);
			if(wstr[0] == '[')
			{
				i=wstr.find(";");
				if(i != std::string::npos) wstr=wstr.substr(i);
                                boost::replace_all(wstr," ","");
                                boost::replace_all(wstr,"[","");
                                boost::replace_all(wstr,"]","");
				// wstr = std::regex_replace(wstr, std::regex(" "), "");
				// wstr = std::regex_replace(wstr, std::regex("["), "");
				// wstr = std::regex_replace(wstr, std::regex("]"), "");
				boost::algorithm::trim(wstr);
				if( wstr == "atoms" )
				{
					//here is atom
					//wresname=wstr;
					//PrintLog("found atoms %s\n",wstr.c_str());
					foundatms=true;
				}
			}
			fgets(strline,1023,fin);
		}
		//read params
		bool stillatoms=true;
		std::vector<std::string> locOPLSResAtomName;
		std::vector<std::string> locOPLSResAtomTypes;
		std::vector<double> locOPLSResAtomCharge;
		while ((!feof(fin))&&stillatoms)
		{
			std::string wstr(strline);
			boost::algorithm::trim(wstr);
			if(wstr[0] == '[')
			{
				stillatoms=false;
			}
			else if(! (wstr[0] == ';'))
			{
				char atmtype[12];
				char name[5];
				double charge;
				int funnynumber;
				if(sscanf(strline,"%s %s %lg %d", name, atmtype, &charge, &funnynumber)==4)
				{
					locOPLSResAtomName.push_back(name);
					locOPLSResAtomTypes.push_back(atmtype);
					locOPLSResAtomCharge.push_back(charge);
				}
				fgets(strline,1023,fin);
			}
			else
			{
				fgets(strline,1023,fin);
			}
			
		}
		OPLSResAtomName.push_back(locOPLSResAtomName);
		OPLSResAtomTypes.push_back(locOPLSResAtomTypes);
		OPLSResAtomCharge.push_back(locOPLSResAtomCharge);
	}
	fclose(fin);
	/*int ires,iat;
	for(ires=0;ires<OPLSResNames.size();ires++)
	{
		PrintLog("RES: %s\n",OPLSResNames[ires].c_str());
		for(iat=0;iat<OPLSResAtomName[ires].size();iat++)
		{
			PrintLog("ATM: %s %s %g\n",OPLSResAtomName[ires][iat].c_str(), OPLSResAtomTypes[ires][iat].c_str(), OPLSResAtomCharge[ires][iat]);
		}
	}*/
	
	return EXIT_SUCCESS;
}
int PNPMod::ReadIER(const char *filename,bool AddToDB)
{
	PrintLog("PNPMod::ReadIER from %s\n",filename);
	FILE *fin=fopen(filename,"rt");
	if(fin==NULL)
	{
		std::string ErrorMsg="Can not open "+std::string(filename);
		ErrorInMod("PNPMod::ReadOPLSrtp", ErrorMsg.c_str());
		return EXIT_FAILURE;
	}
	int i;
	char strline[1024];
	if(!AddToDB)
	{
		IERResNames.clear();
		for(i=0;i<IERAtomName.size();i++)
			IERAtomName[i].clear();
		IERAtomName.clear();
		for(i=0;i<IERRadiusK.size();i++)
			IERRadiusK[i].clear();
		IERRadiusK.clear();
		for(i=0;i<IERRadiusCl.size();i++)
			IERRadiusCl[i].clear();
		IERRadiusCl.clear();
	}
	
	fgets(strline,1023,fin);
	while (!feof(fin))
	{
		//find residue
		bool foundres=false;
		
		while ((!feof(fin))&&(!foundres))
		{
			std::string wstr(strline);
			boost::algorithm::trim(wstr);
			if(wstr[0] == '[')
			{
				i=wstr.find(";");
				if (i != string::npos) wstr = wstr.substr(i);
                                boost::replace_all(wstr," ","");
                                boost::replace_all(wstr,"[","");
                                boost::replace_all(wstr,"]","");
				// wstr = std::regex_replace(wstr, std::regex(" "), "");
				// wstr = std::regex_replace(wstr, std::regex("["), "");
				// wstr = std::regex_replace(wstr, std::regex("]"), "");
				if( wstr == boost::to_upper_copy(wstr) )
				{
					//here is residue
					boost::algorithm::trim(wstr);
					PrintLog("found residue %s\n",wstr.c_str());
					foundres=true;
					IERResNames.push_back(wstr);
				}
			}
			if(!foundres)
				fgets(strline,1023,fin);
		}
		//read params
		bool stillatoms=true;
		std::vector<std::string> locResAtomName;
		std::vector<double> locRadiusK;
		std::vector<double> locRadiusCl;
		
		fgets(strline,1023,fin);
		while ((!feof(fin))&&stillatoms)
		{
			std::string wstr(strline);
			boost::algorithm::trim(wstr);
			
			
			if(wstr[0] == '[')
			{
				stillatoms=false;
			}
			else if(!(wstr[0] == ';'))
			{
				char name[5];
				double rK,rCl;
				if(sscanf(strline,"%s %lg %lg", name, &rK, &rCl)==3)
				{
					locResAtomName.push_back(name);
					locRadiusK.push_back(rK);
					locRadiusCl.push_back(rCl);
				}
			}
			if(stillatoms)
				fgets(strline,1023,fin);
		}
		IERAtomName.push_back(locResAtomName);
		IERRadiusK.push_back(locRadiusK);
		IERRadiusCl.push_back(locRadiusCl);
	}
	fclose(fin);
	/*int ires,iat;
	for(ires=0;ires<IERResNames.size();ires++)
	{
		PrintLog("RES: %s\n",IERResNames[ires].c_str());
		for(iat=0;iat<IERAtomName[ires].size();iat++)
		{
				PrintLog("ATM: %s %g %g\n",IERAtomName[ires][iat].c_str(), IERRadiusK[ires][iat], IERRadiusK[ires][iat]);
		}
	}*/
	PrintLog("PNPMod::Finished reading from %s\n",filename);
	return EXIT_SUCCESS;
}
int PNPMod::SaveIER(const char* filename, bool OnlyHeavyAtoms)
{
	double x, y, z, rK, rCl;
	HaChain  *chain;
	HaResidue  *group;
	HaAtom  *aptr;
	int count;
	
	FILE* DataFile = fopen( filename, "w" );
	if( !DataFile )
	{
		PrintLog("\n");
		return( False );
	}
	
	count = 1;
	MoleculesType::iterator mol_itr;
	for( mol_itr=phost_mset->HostMolecules.begin(); mol_itr != phost_mset->HostMolecules.end(); mol_itr++)
	{
		ChainIteratorMolecule ch_itr(*mol_itr);
		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			for(group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes())
			{
				AtomIteratorAtomGroup aitr_group(group);
				for(aptr= aitr_group.GetFirstAtom(); aptr; aptr = aitr_group.GetNextAtom())
				{
					if( !phost_mset->save_opt_default.save_selected || aptr->Selected())
					{
//						if( prev && (chain->ident!=ch) )
//							fprintf( DataFile, "TER   %5d      %.3s %c%4d \n",
//							count++, prev->GetName(), ch, prev->serno);
						
						if( aptr->flag&HeteroFlag )
						{
							fputs("HETATM",DataFile);
						}
						else
							fputs("ATOM  ",DataFile);
						
						std::string atname(aptr->GetName());
						//std::string atname(aptr->GetFFSymbol());
						if(atname.size() < 4)
							atname.insert(0," ");

						int k;
						for(k=0; k < 3; k++)
						{
							if(atname.size() < 4)
								atname+= " ";
						}
						
						if(atname.size() > 4)
							atname = atname.substr(0,3);
						
						std::string res_name(group->GetName());

						for(k=0; k < 3; k++)
						{
							if(res_name.size() < 3)
								res_name.insert(0," ");
						}
						
						if(res_name.size() > 3)
							res_name = res_name.substr(0,2);
						
						
						fprintf( DataFile, "%5d %.4s %.3s %c%4d ",
										 count++, atname.c_str(), res_name.c_str(),
																					 chain->ident, group->serno );
						
						GetIER(aptr,&rK,&rCl,OnlyHeavyAtoms);
						fprintf(DataFile,"%9.3f%9.3f\n",rK,rCl);
						
					}
				}
			}
			fprintf( DataFile, "TER  \n");
		}
	}
	fputs("END   \n",DataFile);
	fclose( DataFile );
	return( True );
}
int PNPMod::GetResNumAtIERDB(std::string* ResName)
{
	int myres=-1;
	int ires;
	for(ires=0;ires<IERResNames.size();ires++)
	{
		if( *ResName == IERResNames[ires] )
		{
			myres=ires;
			break;
		}
	}
	return myres;
}
int PNPMod::GetAtmNumOfResAtIERDB(int myres,std::string* AtmName)
{
	int myat=-1;
	int iat;
	for(iat=0;iat<IERAtomName[myres].size();iat++)
	{
		if(*AtmName == IERAtomName[myres][iat] )
		{
			//PrintLog("RES: %s found\n",AtmName.c_str());
			//PrintLog("ATM: %s %s %g\n",OPLSResAtomName[myres][iat].c_str(), OPLSResAtomTypes[myres][iat].c_str(), OPLSResAtomCharge[myres][iat]);
			myat=iat;
			break;
		}
	}
	return myat;
}
int PNPMod::GetIER(HaAtom  *aptr, double *rK, double *rCl, bool OnlyHeavyAtoms)
{
	*rK=0.0;
	*rCl=0.0;
	
	HaResidue* res=aptr->GetHostRes();
	HaChain* chn=res->GetHostChain();
	int ResPos=0;//0-inside 1-NTerm -1-CTerm
	
	if(res->GetSerNo()==1)
		ResPos=1;
	if(res->GetSerNo()==chn->GetNRes())
		ResPos=-1;
	
	//Find Res
	std::string ResName(res->GetName());
	int myres=GetResNumAtIERDB(&ResName);
	if(myres==-1)
	{
		PrintLog("Error: resedue %s not found\n",ResName.c_str());
		return FALSE;
	}
	//Find Atom
	std::string AtmName=aptr->GetName();
	int myat=GetAtmNumOfResAtIERDB(myres,&AtmName);
	//Modify Terminal values
	if(ResPos==1)
	{
		//PrintLog("resedue %s SerNo=%d ResPos=%d\n",res->GetName(),res->GetSerNo(),ResPos);
		std::string NTermStr("NTERMALA");
		int NTerm=GetResNumAtIERDB(&NTermStr);
		if(NTerm==-1)
		{
			PrintLog("Error: resedue NTERMALA not found\n",ResName.c_str());
		}
		else
		{
			bool ReSetmyat=false;
			if(AtmName == "N" )ReSetmyat=true;
			if(AtmName == "CA")ReSetmyat=true;
			if(AtmName == "C" )ReSetmyat=true;
			if(AtmName == "O" )ReSetmyat=true;
			if(ReSetmyat)
			{
				myres=NTerm;
				myat=GetAtmNumOfResAtIERDB(myres,&AtmName);
			}
		}
	}
	else if(ResPos==-1)
	{
		//PrintLog("resedue %s SerNo=%d ResPos=%d\n",res->GetName(),res->GetSerNo(),ResPos);
		std::string CTermStr("CTERMALA");
		
		int CTerm=GetResNumAtIERDB(&CTermStr);
		if(CTerm==-1)
		{
			PrintLog("Error: resedue CTERMALA not found\n",ResName.c_str());
		}
		else
		{
			bool ReSetmyat=false;
			if(AtmName == "N" )ReSetmyat=true;
			if(AtmName == "CA")ReSetmyat=true;
			if(AtmName == "C" )ReSetmyat=true;
			if(AtmName == "O" )ReSetmyat=true;
			if(AtmName == "OXT")ReSetmyat=true;
			if(ReSetmyat)
			{
				myres=CTerm;
				myat=GetAtmNumOfResAtIERDB(myres,&AtmName);
			}
		}
	}
	if(myat==-1)
	{
		if(OnlyHeavyAtoms)
		{
			if(aptr->IsHydrogen())
			{
				return TRUE;
			}
			else
			{
				if(ResPos==0)
					PrintLog("Error: atom %s of resedue %s not found\n",AtmName.c_str(),ResName.c_str());
				else
                                        PrintLog("Error: atom %s of resedue %s not found, this is terminal residue so check db of terminals\n",AtmName.c_str(),ResName.c_str());
				return FALSE;
			}
		}
		else
		{
			if(ResPos==0)
				PrintLog("Error: atom %s of resedue %s not found\n",AtmName.c_str(),ResName.c_str());
			else
				PrintLog("Error: atom %s of resedue %s not found, this is terminal residue so check db of terminals\n",AtmName.c_str(),ResName.c_str());
			return FALSE;
		}
	}
	*rK=IERRadiusK[myres][myat];
	*rCl=IERRadiusCl[myres][myat];
	return TRUE;
}
int PNPMod::ReadPANDB(const char *filename,bool AddToDB)
{
	PrintLog("PNPMod::ReadPANDB from %s\n",filename);
	FILE *fin=fopen(filename,"rt");
	if(fin==NULL)
	{
		std::string ErrorMsg="Can not open "+std::string(filename);
		ErrorInMod("PNPMod::ReadPANDB", ErrorMsg.c_str());
		return EXIT_FAILURE;
	}
	int i;
	char strline[1024];
	if(!AddToDB)
	{
		SR_AN_ResNames.clear();
		for(i=0;i<SR_AN_AtomName.size();i++)
			SR_AN_AtomName[i].clear();
		SR_AN_AtomName.clear();
		for(i=0;i<SR_A_K.size();i++)
			SR_A_K[i].clear();
		SR_A_K.clear();
		for(i=0;i<SR_A_Cl.size();i++)
			SR_A_Cl[i].clear();
		SR_A_Cl.clear();
		for(i=0;i<SR_N_K.size();i++)
			SR_N_K[i].clear();
		SR_N_K.clear();
		for(i=0;i<SR_N_Cl.size();i++)
			SR_N_Cl[i].clear();
		SR_N_Cl.clear();
	}
	
	fgets(strline,1023,fin);
	while (!feof(fin))
	{
		//find residue
		bool foundres=false;
		
		while ((!feof(fin))&&(!foundres))
		{
			std::string wstr(strline);
			boost::algorithm::trim(wstr);
			if(wstr[0] == '[')
			{
				i = wstr.find(";");
				if (i != std::string::npos) wstr = wstr.substr(i);
                                boost::replace_all(wstr," ","");
                                boost::replace_all(wstr,"[","");
                                boost::replace_all(wstr,"]","");
				// wstr = std::regex_replace(wstr, std::regex(" "), "");
				// wstr = std::regex_replace(wstr, std::regex("["), "");
				// wstr = std::regex_replace(wstr, std::regex("]"), "");
				if(wstr == boost::to_upper_copy(wstr))
				{
					//here is resedue
					boost::algorithm::trim(wstr);
					PrintLog("found residue %s\n",wstr.c_str());
					foundres=true;
					SR_AN_ResNames.push_back(wstr);
				}
			}
			if(!foundres)
				fgets(strline,1023,fin);
		}
		//read params
		bool stillatoms=true;
		std::vector<std::string> locResAtomName;
		std::vector<double> locAK;
		std::vector<double> locACl;
		std::vector<double> locNK;
		std::vector<double> locNCl;
		
		fgets(strline,1023,fin);
		while ((!feof(fin))&&stillatoms)
		{
			std::string wstr(strline);
			boost::algorithm::trim(wstr);
			
			
			if(wstr[0] == '[')
			{
				stillatoms=false;
			}
			else if(!(wstr[0] == ';'))
			{
				char name[5];
				double AK,ACl,NK,NCl;
				if(sscanf(strline,"%s %lg %lg %lg %lg", name, &AK,&NK,&ACl,&NCl)==5)
				{
					locResAtomName.push_back(name);
					locAK.push_back(AK);
					locNK.push_back(NK);
					locACl.push_back(ACl);
					locNCl.push_back(NCl);
				}
			}
			if(stillatoms)
				fgets(strline,1023,fin);
		}
		SR_AN_AtomName.push_back(locResAtomName);
		SR_A_K.push_back(locAK);
		SR_A_Cl.push_back(locACl);
		SR_N_K.push_back(locNK);
		SR_N_Cl.push_back(locNCl);
	}
	fclose(fin);
	/*int ires,iat;
	for(ires=0;ires<IERResNames.size();ires++)
	{
	PrintLog("RES: %s\n",IERResNames[ires].c_str());
	for(iat=0;iat<IERAtomName[ires].size();iat++)
	{
	PrintLog("ATM: %s %g %g\n",IERAtomName[ires][iat].c_str(), IERRadiusK[ires][iat], IERRadiusK[ires][iat]);
}
}*/
	PrintLog("PNPMod::Finished reading from %s\n",filename);
	return EXIT_SUCCESS;
}
int PNPMod::SavePAN(const char* filename, bool OnlyHeavyAtoms)
{
	double x, y, z;
	double AK,ACl,NK,NCl;
	HaChain  *chain;
	HaResidue  *group;
	HaAtom  *aptr;
	int count;
	
	FILE* DataFile = fopen( filename, "w" );
	if( !DataFile )
	{
		PrintLog("\n");
		return( False );
	}
	
	count = 1;
	MoleculesType::iterator mol_itr;
	for( mol_itr=phost_mset->HostMolecules.begin(); mol_itr != phost_mset->HostMolecules.end(); mol_itr++)
	{
		ChainIteratorMolecule ch_itr(*mol_itr);
		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			for(group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes())
			{
				AtomIteratorAtomGroup aitr_group(group);
				for(aptr= aitr_group.GetFirstAtom(); aptr; aptr = aitr_group.GetNextAtom())
				{
					if( !phost_mset->save_opt_default.save_selected || aptr->Selected())
					{
//						if( prev && (chain->ident!=ch) )
//							fprintf( DataFile, "TER   %5d      %.3s %c%4d \n",
//							count++, prev->GetName(), ch, prev->serno);
						
						if( aptr->flag&HeteroFlag )
						{
							fputs("HETATM",DataFile);
						}
						else
							fputs("ATOM  ",DataFile);
						
						std::string atname(aptr->GetName());
						//std::string atname(aptr->GetFFSymbol());
						if(atname.size() < 4)
							atname.insert(0," ");

						int k;
						for(k=0; k < 3; k++)
						{
							if(atname.size() < 4)
								atname+= " ";
						}
						
						if(atname.size() > 4)
							atname = atname.substr(0,3);
						
						std::string res_name(group->GetName());

						for(k=0; k < 3; k++)
						{
							if(res_name.size() < 3)
								res_name.insert(0," ");
						}
						
						if(res_name.size() > 3)
							res_name = res_name.substr(0,2);
						
						
						fprintf( DataFile, "%5d %.4s %.3s %c%4d ",
										 count++, atname.c_str(), res_name.c_str(),
																					 chain->ident, group->serno );
						
						GetSR_AN(aptr,&AK,&NK,&ACl,&NCl,OnlyHeavyAtoms);
						fprintf(DataFile,"%9.3f%9.3f%9.3f%9.3f\n",AK,NK,ACl,NCl);
						
					}
				}
			}
			fprintf( DataFile, "TER  \n");
		}
	}
	fputs("END   \n",DataFile);
	fclose( DataFile );
	return( True );
}
int PNPMod::AssignPAN(bool OnlyHeavyAtoms)
{
	double x, y, z;
	double AK,ACl,NK,NCl;
	HaChain  *chain;
	HaResidue  *group;
	HaAtom  *aptr;
	int count;
	
	mSR_A_K.clear();
	mSR_A_Cl.clear();
	mSR_N_K.clear();
	mSR_N_Cl.clear();
	
	count = 1;
	MoleculesType::iterator mol_itr;
	for( mol_itr=phost_mset->HostMolecules.begin(); mol_itr != phost_mset->HostMolecules.end(); mol_itr++)
	{
		ChainIteratorMolecule ch_itr(*mol_itr);
		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			for(group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes())
			{
				AtomIteratorAtomGroup aitr_group(group);
				for(aptr= aitr_group.GetFirstAtom(); aptr; aptr = aitr_group.GetNextAtom())
				{
					if( !phost_mset->save_opt_default.save_selected || aptr->Selected())
					{
//						if( prev && (chain->ident!=ch) )
//							fprintf( DataFile, "TER   %5d      %.3s %c%4d \n",
//							count++, prev->GetName(), ch, prev->serno);
						
						
						std::string atname(aptr->GetName());
						//std::string atname(aptr->GetFFSymbol());
						if(atname.size() < 4)
							atname.insert(0," ");

						int k;
						for(k=0; k < 3; k++)
						{
							if(atname.size() < 4)
								atname+= " ";
						}
						
						if(atname.size() > 4)
							atname = atname.substr(0,3);
						
						std::string res_name(group->GetName());

						for(k=0; k < 3; k++)
						{
							if(res_name.size() < 3)
								res_name.insert(0," ");
						}
						
						if(res_name.size() > 3)
							res_name = res_name.substr(0,2);
						
						
						//fprintf( DataFile, "%5d %.4s %.3s %c%4d ",
						//				 count++, atname.c_str(), res_name.c_str(),
						//															 chain->ident, group->serno );
						
						GetSR_AN(aptr,&AK,&NK,&ACl,&NCl,OnlyHeavyAtoms);
						//fprintf(DataFile,"%9.3f%9.3f%9.3f%9.3f\n",AK,NK,ACl,NCl);
						mSR_A_K.push_back(AK);
						mSR_A_Cl.push_back(ACl);
						mSR_N_K.push_back(NK);
						mSR_N_Cl.push_back(NCl);
						
					}
				}
			}
		}
	}
	return( True );
}
int PNPMod::GetSR_AN(HaAtom  *aptr, double *AK, double *NK, double *ACl, double *NCl, bool OnlyHeavyAtoms)
{
	*AK=0.0;
	*NK=0.0;
	*ACl=0.0;
	*NCl=0.0;
	
	HaResidue* res=aptr->GetHostRes();
	HaChain* chn=res->GetHostChain();
	int ResPos=0;//0-inside 1-NTerm -1-CTerm
	
	std::string nameCT("CT");
	std::string nameNT("NT");
	std::string NameModifier=res->GetNameModifier();
	
	if( NameModifier == nameNT )
		ResPos=1;
	if( NameModifier == nameCT )
		ResPos=-1;
	
	//Find Res
	std::string ResName(res->GetName());
	int myres=GetResNumAtSR_AN_DB(&ResName);
	if(myres==-1)
	{
		PrintLog("Error: resedue %s not found\n",ResName.c_str());
		return FALSE;
	}
	//Find Atom
	std::string AtmName=aptr->GetName();
	int myat=GetAtmNumOfResAtSR_AN_DB(myres,&AtmName);
	//Modify Terminal values
	if(ResPos==1)
	{
		//PrintLog("resedue %s SerNo=%d ResPos=%d\n",res->GetName(),res->GetSerNo(),ResPos);
		std::string NTermStr("NTERMALA");
		int NTerm=GetResNumAtSR_AN_DB(&NTermStr);
		if(NTerm==-1)
		{
			PrintLog("Error: resedue NTERMALA not found\n",ResName.c_str());
		}
		else
		{
			bool ReSetmyat=false;
			if(AtmName == "N" )ReSetmyat=true;
			if(AtmName == "CA")ReSetmyat=true;
			if(AtmName == "C" )ReSetmyat=true;
			if(AtmName == "O" )ReSetmyat=true;
			if(ReSetmyat)
			{
				myres=NTerm;
				myat=GetAtmNumOfResAtSR_AN_DB(myres,&AtmName);
			}
		}
	}
	else if(ResPos==-1)
	{
		//PrintLog("resedue %s SerNo=%d ResPos=%d\n",res->GetName(),res->GetSerNo(),ResPos);
		std::string CTermStr("CTERMALA");
		
		int CTerm=GetResNumAtSR_AN_DB(&CTermStr);
		if(CTerm==-1)
		{
			PrintLog("Error: resedue CTERMALA not found\n",ResName.c_str());
		}
		else
		{
			bool ReSetmyat=false;
			if(AtmName == "N"  )ReSetmyat=true;
			if(AtmName == "CA" )ReSetmyat=true;
			if(AtmName == "C"  )ReSetmyat=true;
			if(AtmName == "O"  )ReSetmyat=true;
			if(AtmName == "OXT")ReSetmyat=true;
			if(ReSetmyat)
			{
				myres=CTerm;
				myat=GetAtmNumOfResAtSR_AN_DB(myres,&AtmName);
			}
		}
	}
	if(myat==-1)
	{
		if(OnlyHeavyAtoms)
		{
			if(aptr->IsHydrogen())
			{
				return TRUE;
			}
			else
			{
				if(ResPos==0)
					PrintLog("Error: atom %s of residue %s not found\n",AtmName.c_str(),ResName.c_str());
				else
					PrintLog("Error: atom %s of residue %s not found, this is terminal residue so check db of terminals\n",AtmName.c_str(),ResName.c_str());
				return FALSE;
			}
		}
		else
		{
			if(ResPos==0)
				PrintLog("Error: atom %s of residue %s not found\n",AtmName.c_str(),ResName.c_str());
			else
				PrintLog("Error: atom %s of residue %s not found, this is terminal residue so check db of terminals\n",AtmName.c_str(),ResName.c_str());
			return FALSE;
		}
	}
	*AK=SR_A_K[myres][myat];
	*NK=SR_N_K[myres][myat];
	*ACl=SR_A_Cl[myres][myat];
	*NCl=SR_N_Cl[myres][myat];
	return TRUE;
}
int PNPMod::GetResNumAtSR_AN_DB(std::string* ResName)
{
	int myres=-1;
	int ires;
	for(ires=0;ires<SR_AN_ResNames.size();ires++)
	{
		if(*ResName == SR_AN_ResNames[ires])
		{
			myres=ires;
			break;
		}
	}
	return myres;
}
int PNPMod::GetAtmNumOfResAtSR_AN_DB(int myres,std::string* AtmName)
{
	int myat=-1;
	int iat;
	for(iat=0;iat<SR_AN_AtomName[myres].size();iat++)
	{
		if( *AtmName == SR_AN_AtomName[myres][iat] )
		{
			//PrintLog("RES: %s found\n",AtmName.c_str());
			//PrintLog("ATM: %s %s %g\n",OPLSResAtomName[myres][iat].c_str(), OPLSResAtomTypes[myres][iat].c_str(), OPLSResAtomCharge[myres][iat]);
			myat=iat;
			break;
		}
	}
	return myat;
}

int PNPMod::GetAtomTypeNumber(std::string atomname)
{
	static int count=0;
	int i;
	
	for(i=0; i < AtomTypesDB.size(); i++) 
	{
		if(AtomTypesDB[i].compare(atomname)==0)
		{
			return i;
		}
	}
	std::cout << "Can not find atom with FFsymbol: " << atomname << endl;
	return -1;
}
double PNPMod::GetHalfSigma(std::string atomname)
{
	int AtomTypeNumber=GetAtomTypeNumber(atomname);
	if(AtomTypeNumber<0)
		return 0.0;
	else
		return HalfSigmaDB[AtomTypeNumber];
}
double PNPMod::GetFourEpsilon(std::string atomname)
{
	int AtomTypeNumber=GetAtomTypeNumber(atomname);
	if(AtomTypeNumber<0)
		return 0.0;
	else
		return FourEpsilonDB[AtomTypeNumber];
}
int PNPMod::SavePREFile(const char* filename)
{
	double x, y, z;
	HaChain  *chain;
	HaResidue  *group;
	HaAtom  *aptr;
	int count;
	
	FILE* DataFile = fopen( filename, "w" );
	if( !DataFile )
	{
		PrintLog("\n");
		return( False );
	}
	
	count = 1;
	MoleculesType::iterator mol_itr;
	for( mol_itr=phost_mset->HostMolecules.begin(); mol_itr != phost_mset->HostMolecules.end(); mol_itr++)
	{
		ChainIteratorMolecule ch_itr(*mol_itr);
		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			for(group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes())
			{
				AtomIteratorAtomGroup aitr_group(group);
				for(aptr= aitr_group.GetFirstAtom(); aptr; aptr = aitr_group.GetNextAtom())
				{
					if( !phost_mset->save_opt_default.save_selected || aptr->Selected())
					{
//						if( prev && (chain->ident!=ch) )
//							fprintf( DataFile, "TER   %5d      %.3s %c%4d \n",
//							count++, prev->GetName(), ch, prev->serno);
						
						if( aptr->flag&HeteroFlag )
						{
							fputs("HETATM",DataFile);
						}
						else
							fputs("ATOM  ",DataFile);
						
						std::string atname(aptr->GetName());
						//std::string atname(aptr->GetFFSymbol());
						if(atname.size() < 4)
							atname.insert(0," ");

						int k;
						for(k=0; k < 3; k++)
						{
							if(atname.size() < 4)
								atname+= " ";
						}
						
						if(atname.size() > 4)
							atname = atname.substr(0,3);
						
						std::string res_name(group->GetName());

						for(k=0; k < 3; k++)
						{
							if(res_name.size() < 3)
								res_name.insert(0," ");
						}
						
						if(res_name.size() > 3)
							res_name = res_name.substr(0,2);
						
						
						fprintf( DataFile, "%5d %.4s %.3s %c%4d    ",
										 count++, atname.c_str(), res_name.c_str(),
																					 chain->ident, group->serno );
						
						HaMolView* pView = phost_mset->GetActiveMolView();
						if(phost_mset->save_opt_default.save_transform && pView != NULL)
						{
							pView->GetTransfCoord(aptr->GetX_Ang(),aptr->GetY_Ang(),aptr->GetZ_Ang(),x,y,z);
						}
						else
						{
							x = aptr->GetX_Ang();
							y = aptr->GetY_Ang();
							z = aptr->GetZ_Ang();
						}
						
						fprintf(DataFile,"%8.3f%8.3f%8.3f",x,y,z);
						
						std::string ffname(aptr->GetFFSymbol());
						fprintf(DataFile,"%7.3f%6.3f\n",GetHalfSigma(ffname),GetFourEpsilon(ffname));
						
					}
				}
			}
			fprintf( DataFile, "TER  \n");
		}
	}
	fputs("END   \n",DataFile);
	fclose( DataFile );
	return( True );
}
int PNPMod::SetLJABfromAMBERFF()
{
	AtmsSet.clear();
	DeleteObjByPnt(LJA);
	DeleteObjByPnt(LJB);
	
	//Find epsstar and rstar for Ions
	double *epsstar=new double[NIonsTypes];
	double *rstar=new double[NIonsTypes];
	
	int i;
	for(i=0;i<NIonsTypes;i++)
	{
		int iDB=GetAtomTypeNumber(IonNames[i]);
		if(iDB>=0)
		{
			epsstar[i]=FourEpsilonDB[iDB];
			rstar[i]=HalfSigmaDB[iDB];
			PrintLog("epsstar[%d]=%g rstar[%d]=%g\n",i,epsstar[i],i,rstar[i]);
		}
	}
	
	//Set A and B of Protein with Ions
	//Make AtmsSet list with atoms to set
	HaChain  *chain;
	HaResidue  *group;
	HaAtom  *aptr;
	MoleculesType::iterator mol_itr;
	for( mol_itr=phost_mset->HostMolecules.begin(); mol_itr != phost_mset->HostMolecules.end(); mol_itr++)
	{
		ChainIteratorMolecule ch_itr(*mol_itr);
		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			for(group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes())
			{
				AtomIteratorAtomGroup aitr_group(group);
				for(aptr= aitr_group.GetFirstAtom(); aptr; aptr = aitr_group.GetNextAtom())
				{
					if( !phost_mset->save_opt_default.save_selected || aptr->Selected())
					{
						AtmsSet.push_back(aptr);
					}
				}
			}
		}
	}
	//
	LJA=new HaMat_double(NIonsTypes,AtmsSet.size());
	LJB=new HaMat_double(NIonsTypes,AtmsSet.size());
	LJA->set(0.0);
	LJB->set(0.0);
	//Calc A and B
	int ion,atm;
	for(atm=0;atm<AtmsSet.size();atm++)
	{
		std::string ffname(AtmsSet[atm]->GetFFSymbol());
		int iDB=GetAtomTypeNumber(ffname);
		if(iDB>=0)
		{
			double EpsStarAtm=FourEpsilonDB[iDB];
			double RStarAtm=HalfSigmaDB[iDB];
			for(ion=0;ion<NIonsTypes;ion++)
			{
				double eps=sqrt(epsstar[ion]*EpsStarAtm);
				double r=rstar[ion]+RStarAtm;
				double r3=r*r*r;
				double r6=r3*r3;
				double r12=r6*r6;
				LJA->r0(ion,atm)=eps*r12;
				LJB->r0(ion,atm)=2.0*eps*r6;
			}
		}
	}
	
	delete [] epsstar;
	delete [] rstar;
	return EXIT_SUCCESS;
}
int PNPMod::PrintLJAB()
{
	int ion,atm;
	HaChain  *chain;
	HaResidue  *group;
	HaAtom  *aptr;
	PrintLog("%5s%5s%5s%5s","Atm","AtmN","Res","ResN");
	for(ion=0;ion<NIonsTypes;ion++)
	{
		PrintLog("%14s%14s%10s","A","B","sigma");
	}
	PrintLog("\n");
	for(atm=0;atm<AtmsSet.size();atm++)
	{
		aptr=AtmsSet[atm];
		group=aptr->GetHostRes();
		chain=aptr->GetHostChain();
		PrintLog("%5s%5d%5s%5d",aptr->GetName(),atm+1
				,group->GetName(),group->serno);
		for(ion=0;ion<NIonsTypes;ion++)
		{
			double A=LJA->r0(ion,atm);
			double B=LJB->r0(ion,atm);
			double sigma;
			if(A==0.0)
				sigma=0.0;
			else
				sigma=pow(A/B,1.0/6.0);
			
			PrintLog("%14.6e%14.6e%10.6f",A,B,sigma);
		}
		PrintLog("\n");
	}
	return EXIT_SUCCESS;
}
int PNPMod::SavePABFile(const char* filename)
{
	double x, y, z,A,B;
	HaChain  *chain;
	HaResidue  *group;
	HaAtom  *aptr;
	int count=1;
	
	FILE* DataFile = fopen( filename, "w" );
	if( !DataFile )
	{
		PrintLog("\n");
		return( False );
	}
	int ion,atm;
	for(atm=0;atm<AtmsSet.size();atm++)
	{
		aptr=AtmsSet[atm];
		group=aptr->GetHostRes();
		chain=aptr->GetHostChain();
		
		if( aptr->flag&HeteroFlag )
		{
			fputs("HETATM",DataFile);
		}
		else
			fputs("ATOM  ",DataFile);
						
		std::string atname(aptr->GetName());
		
		if(atname.size() < 4)
			atname.insert(0," ");
		
		int k;
		for(k=0; k < 3; k++)
		{
			if(atname.size() < 4)
				atname+= " ";
		}
		
		if(atname.size() > 4)
			atname = atname.substr(0,3);
		
		std::string res_name(group->GetName());
		
		for(k=0; k < 3; k++)
		{
			if(res_name.size() < 3)
				res_name.insert(0," ");
		}
		
		if(res_name.size() > 3)
			res_name = res_name.substr(0,2);
		
		fprintf( DataFile, "%5d %.4s %.3s %c%4d    ",
			count++, atname.c_str(), res_name.c_str(),
			chain->ident, group->serno );
		
		for(ion=0;ion<NIonsTypes;ion++)
		{
			fprintf(DataFile,"%14.6e%14.6e",LJA->r0(ion,atm),LJB->r0(ion,atm));
		}
		fprintf(DataFile,"\n");
	}
	fputs("END   \n",DataFile);
	fclose( DataFile );
	return EXIT_SUCCESS;
}
int PNPMod::GetOPLSEpsilonSigma(const char *AtmNM, const char *ResNM, double *Eps, double *Sgm)
{
	*Eps=0.0;
	*Sgm=0.0;
	//Find Res
	std::string ResName(ResNM);
	
	int myres=-1;
	int ires;
	for(ires=0;ires<OPLSResNames.size();ires++)
	{
		if( ResName == OPLSResNames[ires] )
		{
			//PrintLog("RES: %s found\n",ResName.c_str());
			myres=ires;
			break;
		}
	}
	if(myres==-1)
	{
		PrintLog("Error: residue %s not found\n",ResName.c_str());
		return EXIT_FAILURE;
	}
	//Find Atom
	std::string AtmName(AtmNM);
	
	int myat=-1;
	int iat;
	for(iat=0;iat<OPLSResAtomName[myres].size();iat++)
	{
		if(AtmName == OPLSResAtomName[myres][iat] )
		{
			//PrintLog("RES: %s found\n",AtmName.c_str());
			//PrintLog("ATM: %s %s %g\n",OPLSResAtomName[myres][iat].c_str(), OPLSResAtomTypes[myres][iat].c_str(), OPLSResAtomCharge[myres][iat]);
			myat=iat;
			break;
		}
	}
	if(myat==-1)
	{
		PrintLog("Error: atom %s of resedue %s not found\n",AtmName.c_str(),ResName.c_str());
		return EXIT_FAILURE;
	}
	//Find AtomType in LJDB
	std::string AtmType=OPLSResAtomTypes[myres][myat];
	int myatLJ=-1;
	for(iat=0;iat<OPLSLJAtomTypes.size();iat++)
	{
		if(AtmType == OPLSLJAtomTypes[iat] )
		{
			//PrintLog("RES: %s found\n",AtmName.c_str());
			//PrintLog("ATM: %s %s %g\n",OPLSResAtomName[myres][iat].c_str(), OPLSResAtomTypes[myres][iat].c_str(), OPLSResAtomCharge[myres][iat]);
			myatLJ=iat;
			break;
		}
	}
	if(myatLJ==-1)
	{
		PrintLog("Error: type of %s for atom %s of resedue %s not found\n",AtmType.c_str(), AtmName.c_str(), ResName.c_str());
		return EXIT_FAILURE;
	}
	double s2=OPLSLJSigma[myatLJ]*OPLSLJSigma[myatLJ];
	double s6=s2*s2*s2;
	
	*Eps=OPLSLJEpsilon[myatLJ];
	*Sgm=OPLSLJSigma[myatLJ];
	return EXIT_SUCCESS;
}
int PNPMod::SetLJABfromOPLS()
{
	AtmsSet.clear();
	DeleteObjByPnt(LJA);
	DeleteObjByPnt(LJB);
	//Find epsstar and rstar for Ions
	double *epsstar=new double[NIonsTypes];
	double *rstar=new double[NIonsTypes];
	
	int i;
	for(i=0;i<NIonsTypes;i++)
	{
		double EpsStarAtm, RStarAtm;
		GetOPLSEpsilonSigma(IonNames[i].c_str(),IonNames[i].c_str(), &EpsStarAtm, &RStarAtm);
		
		epsstar[i]=EpsStarAtm;
		rstar[i]=RStarAtm;
		PrintLog("epsstar[%d]=%g rstar[%d]=%g\n",i,epsstar[i],i,rstar[i]);
	}
	//Set A and B of Protein with Ions
	//Make AtmsSet list with atoms to set
	HaChain  *chain;
	HaResidue  *group;
	HaAtom  *aptr;
	MoleculesType::iterator mol_itr;
	for( mol_itr=phost_mset->HostMolecules.begin(); mol_itr != phost_mset->HostMolecules.end(); mol_itr++)
	{
		ChainIteratorMolecule ch_itr(*mol_itr);
		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			for(group=ritr_ch.GetFirstRes(); group; group= ritr_ch.GetNextRes())
			{
				AtomIteratorAtomGroup aitr_group(group);
				for(aptr= aitr_group.GetFirstAtom(); aptr; aptr = aitr_group.GetNextAtom())
				{
					if( !phost_mset->save_opt_default.save_selected || aptr->Selected())
					{
						AtmsSet.push_back(aptr);
					}
				}
			}
		}
	}
	//
	LJA=new HaMat_double(NIonsTypes,AtmsSet.size());
	LJB=new HaMat_double(NIonsTypes,AtmsSet.size());
	LJA->set(0.0);
	LJB->set(0.0);
	//Calc A and B
	int ion,atm;
	for(atm=0;atm<AtmsSet.size();atm++)
	{
		std::string ffname(AtmsSet[atm]->GetFFSymbol());
		double EpsStarAtm, RStarAtm;
		GetOPLSEpsilonSigma(AtmsSet[atm]->GetName(), AtmsSet[atm]->GetHostRes()->GetName(), &EpsStarAtm, &RStarAtm);
		
		if(EpsStarAtm>0.0)
		{
			for(ion=0;ion<NIonsTypes;ion++)
			{
				double eps=sqrt(epsstar[ion]*EpsStarAtm);
				double r=0.5*(rstar[ion]+RStarAtm);
				double r3=r*r*r;
				double r6=r3*r3;
				double r12=r6*r6;
				LJA->r0(ion,atm)=4.0*eps*r12;
				LJB->r0(ion,atm)=4.0*eps*r6;
			}
		}
	}
	
	delete [] epsstar;
	delete [] rstar;
	return EXIT_SUCCESS;
}
int PNPMod::SavePREFreeFile(const char* filename)
{
	double x, y, z;
	HaChain  *chain;
	HaResidue  *group;
	HaAtom  *aptr;
	int count;
	
	FILE* DataFile = fopen( filename, "w" );
	if( !DataFile )
	{
		PrintLog("\n");
		return( False );
	}
	
	count = 1;
	MoleculesType::iterator mol_itr;
	for( mol_itr=phost_mset->HostMolecules.begin(); mol_itr != phost_mset->HostMolecules.end(); mol_itr++)
	{
		ChainIteratorMolecule ch_itr(*mol_itr);
		for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
		{
			ResidueIteratorChain ritr_ch(chain);
			for(group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes())
			{
				AtomIteratorAtomGroup aitr_group(group);
				for(aptr= aitr_group.GetFirstAtom(); aptr; aptr = aitr_group.GetNextAtom())
				{
					if( !phost_mset->save_opt_default.save_selected || aptr->Selected())
					{
//						if( prev && (chain->ident!=ch) )
//							fprintf( DataFile, "TER   %5d      %.3s %c%4d \n",
//							count++, prev->GetName(), ch, prev->serno);
						
						if( aptr->flag&HeteroFlag )
						{
							fputs("HETATM",DataFile);
						}
						else
							fputs("ATOM  ",DataFile);
						
						std::string atname(aptr->GetName());
						//std::string atname(aptr->GetFFSymbol());
						if(atname.size() < 4)
							atname.insert(0," ");

						int k;
						for(k=0; k < 3; k++)
						{
							if(atname.size() < 4)
								atname+= " ";
						}
						
						if(atname.size() > 4)
							atname = atname.substr(0,3);
						
						std::string res_name(group->GetName());

						for(k=0; k < 3; k++)
						{
							if(res_name.size() < 3)
								res_name.insert(0," ");
						}
						
						if(res_name.size() > 3)
							res_name = res_name.substr(0,2);
						
						
						fprintf( DataFile, "%5d %.4s %.3s %c%4d    ",
										 count++, atname.c_str(), res_name.c_str(),
																					 chain->ident, group->serno );
						
						HaMolView* pView = phost_mset->GetActiveMolView();
						if(phost_mset->save_opt_default.save_transform && pView != NULL)
						{
							pView->GetTransfCoord(aptr->GetX_Ang(), aptr->GetY_Ang(), aptr->GetZ_Ang(),x,y,z);
						}
						else
						{
							x = aptr->GetX_Ang();
							y = aptr->GetY_Ang();
							z = aptr->GetZ_Ang();
						}
						
						fprintf(DataFile,"%8.3f%8.3f%8.3f",x,y,z);
						
						std::string ffname(aptr->GetFFSymbol());
						fprintf(DataFile,"%10.6f%14.6e\n",GetHalfSigma(ffname),GetFourEpsilon(ffname));
						
					}
				}
			}
			fprintf( DataFile, "TER  \n");
		}
	}
	fputs("END   \n",DataFile);
	fclose( DataFile );
	return( True );
}
int PNPMod::RunPNPSFromString(const char* string)
{
	/*std::istringstream ins(string);
	TiXmlElement* Elt=new TiXmlElement("DontKnowName");
	ins>>*Elt;
	pnpsapp->RunFromTiXmlElement(Elt);*/
	return EXIT_FAILURE;
}
///////////////////////////////////////////////////////////////////////////////
#ifdef ELMOD_COMPILE
ElModRadDist::ElModRadDist(MolSet* new_phost_mset)
{
	pmset=new_phost_mset;
	elmod=pmset->GetElMod(FALSE);
	
	Rmin=0.0;
	Rmax=0.0;
	BinSize=0.0;
	nbins=0;
	gcont[0]=NULL;
	gcont[1]=NULL;
}
ElModRadDist::~ElModRadDist()
{
	int i;
	for(i=0;i<AtomsList.size();i++)
	{
		delete gcont[0][i];
		delete gcont[1][i];
	}
	delete [] gcont[0];
	delete [] gcont[1];
}
int ElModRadDist::CalcRadDistByGrid()
{
	if(pmset==NULL)
	{
		PrintLog("pmset is not given");
		return EXIT_FAILURE;
	}
	if(elmod==NULL)
	{
		PrintLog("elmod is not run yet");
		return EXIT_FAILURE;
	}
	
	gcont[0]=new HaVec_double*[AtomsList.size()];
	gcont[1]=new HaVec_double*[AtomsList.size()];
	
	int i,j,k;
	for(i=0;i<AtomsList.size();i++)
	{
		gcont[0][i]=new HaVec_double(nbins);
		gcont[0][i]->set(0.0);
		gcont[1][i]=new HaVec_double(nbins);
		gcont[1][i]->set(0.0);
	}
	
	HaAtom** HaAtomsList=new HaAtom*[AtomsList.size()];
	
	HaAtom* aptr;
	AtomIteratorMolSet aitr(pmset);
	for(i=0;i<AtomsList.size();i++)
	{
		j=1;
		HaAtomsList[i]=NULL;
		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if(j==AtomsList[i])
			{
				HaAtomsList[i]=aptr;
				break;
			}
			j++;
		}
	}
	
	float Cbulk=elmod->rionst;
	
	int GSX=elmod->ConcMap0.GetNx();
	int GSY=elmod->ConcMap0.GetNy();
	int GSZ=elmod->ConcMap0.GetNz();
	int GSXY=GSX*GSY;
	int GridSize[3]={GSX,GSY,GSZ};
	float GridScale=1.0/elmod->ConcMap0.GetXstep();
	
	float locRmax=Rmax*GridScale;
	float locRmin=Rmin*GridScale;
	float locBinSize=BinSize*GridScale;
	int iRmax=int(locRmax)+1;
	
	PrintLog("GS=[%d %d %d] GridScale=%g\n",GSX,GSY,GSZ,GridScale);
	
	int *icount=new int[nbins];
	
	for(i=0;i<AtomsList.size();i++)
	{
		aptr=HaAtomsList[i];
		
		float Ratom[3];
		Ratom[0]=aptr->GetX()*GridScale+(float)int(GSX/2);
		Ratom[1]=aptr->GetY()*GridScale+(float)int(GSY/2);
		Ratom[2]=aptr->GetZ()*GridScale+(float)int(GSZ/2);
		
		int start[3];
		int end[3];
		for(k=0;k<3;k++)
		{
			start[k]=Ratom[k]-iRmax;
			end[k]=Ratom[k]+iRmax;
			if(start[k]<0)start[k]=0;
			if(end[k]>GridSize[k]-1)end[k]=GridSize[k]-1;
		}
		for(k=0;k<2;k++)
		{
			float* C;
			if(k==0)C=elmod->ConcMap0.GetFieldPtr();
			else C=elmod->ConcMap1.GetFieldPtr();
			
			for(j=0;j<nbins;j++)
				icount[j]=0;
			
			double *gpnt=gcont[k][i]->begin();
			
			int GrdPnt;
			int ix,iy,iz;
			int ibin;
			float RtmpSQ,Rtmp;
			for(ix=start[0];ix<=end[0];ix++)
				for(iy=start[1];iy<=end[1];iy++)
					for(iz=start[2];iz<=end[2];iz++)
			{
				GrdPnt=ix+iy*GSX+iz*GSXY;
				Rtmp=(sqrt((Ratom[0]-ix)*(Ratom[0]-ix)+(Ratom[1]-iy)*(Ratom[1]-iy)+(Ratom[2]-iz)*(Ratom[2]-iz))-locRmin)/locBinSize;
				ibin=int(Rtmp+0.5);
				if(ibin>0 && ibin<nbins)
				{
					gpnt[ibin]+=C[GrdPnt]/Cbulk;
					icount[ibin]++;
				}
			}
			for(j=0;j<nbins;j++)
			{
				if(icount[j]>0)
					gpnt[j]/=icount[j];
			}
		}
		
	}
	//
	//
	
	delete [] HaAtomsList;
	return EXIT_SUCCESS;
}
int ElModRadDist::CalcWelRadDistByGrid()
{
/*  if(pmset==NULL)
  {
    PrintLog("pmset is not given");
    return EXIT_FAILURE;
  }
  if(elmod==NULL)
  {
    PrintLog("elmod is not run yet");
    return EXIT_FAILURE;
  }
  
  wel=new HaVec_double*[AtomsList.size()];
  
  int i,j,k;
  for(i=0;i<AtomsList.size();i++)
  {
    wel[i]=new HaVec_double(nbins);
    wel[i]->set(0.0);
  }
  
  HaAtom** HaAtomsList=new HaAtom*[AtomsList.size()];
  
  HaAtom* aptr;
  AtomIteratorMolSet aitr(pmset);
  for(i=0;i<AtomsList.size();i++)
  {
    j=1;
    HaAtomsList[i]=NULL;
    for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
    {
      if(j==AtomsList[i])
      {
        HaAtomsList[i]=aptr;
        break;
      }
      j++;
    }
  }
  
  float Cbulk=elmod->rionst;
  
  int GSX=elmod->ConcMap0.GetNx();
  int GSY=elmod->ConcMap0.GetNy();
  int GSZ=elmod->ConcMap0.GetNz();
  int GSXY=GSX*GSY;
  int GridSize[3]={GSX,GSY,GSZ};
  float GridScale=1.0/elmod->ConcMap0.GetXstep();
  
  float locRmax=Rmax*GridScale;
  float locRmin=Rmin*GridScale;
  float locBinSize=BinSize*GridScale;
  int iRmax=int(locRmax)+1;
  
  PrintLog("GS=[%d %d %d] GridScale=%g\n",GSX,GSY,GSZ,GridScale);
  
  int *icount=new int[nbins];
  
  for(i=0;i<AtomsList.size();i++)
  {
    aptr=HaAtomsList[i];
    
    float Ratom[3];
    Ratom[0]=aptr->GetX()*GridScale+(float)int(GSX/2);
    Ratom[1]=aptr->GetY()*GridScale+(float)int(GSY/2);
    Ratom[2]=aptr->GetZ()*GridScale+(float)int(GSZ/2);
    
    int start[3];
    int end[3];
    for(k=0;k<3;k++)
    {
      start[k]=Ratom[k]-iRmax;
      end[k]=Ratom[k]+iRmax;
      if(start[k]<0)start[k]=0;
      if(end[k]>GridSize[k]-1)end[k]=GridSize[k]-1;
    }
    for(k=0;k<2;k++)
    {
      float* C;
      if(k==0)C=elmod->ConcMap0.GetFieldPtr();
      else C=elmod->ConcMap1.GetFieldPtr();
      
      for(j=0;j<nbins;j++)
        icount[j]=0;
      
      double *gpnt=gcont[k][i]->begin();
      
      int GrdPnt;
      int ix,iy,iz;
      int ibin;
      float RtmpSQ,Rtmp;
      for(ix=start[0];ix<=end[0];ix++)
        for(iy=start[1];iy<=end[1];iy++)
          for(iz=start[2];iz<=end[2];iz++)
      {
        GrdPnt=ix+iy*GSX+iz*GSXY;
        Rtmp=(sqrt((Ratom[0]-ix)*(Ratom[0]-ix)+(Ratom[1]-iy)*(Ratom[1]-iy)+(Ratom[2]-iz)*(Ratom[2]-iz))-locRmin)/locBinSize;
        ibin=int(Rtmp+0.5);
        if(ibin>0 && ibin<nbins)
        {
          gpnt[ibin]+=C[GrdPnt]/Cbulk;
          icount[ibin]++;
        }
      }
      for(j=0;j<nbins;j++)
      {
        if(icount[j]>0)
          gpnt[j]/=icount[j];
      }
    }
    
  }
  //
  //
  
  delete [] HaAtomsList;*/
  return EXIT_SUCCESS;
}
int ElModRadDist::CalcRadDist()
{
/*	if(pmset==NULL)
	{
		PrintLog("pmset is not given");
		return EXIT_FAILURE;
	}
	if(elmod==NULL)
	{
		PrintLog("elmod is not run yet");
		return EXIT_FAILURE;
	}
	
	int i=1;
	HaAtom* aptr;
	AtomIteratorMolSet aitr(pmset);
	
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if(i==atom)break;
		i++;
	}
	float Ratom[3];
	Ratom[0]=aptr->GetX();
	Ratom[1]=aptr->GetY();
	Ratom[2]=aptr->GetZ();
	
	PrintLog("ElModRadDist::CalcRadDist\n");
	PrintLog(" Rmin=%g A Rmax=%g A BinSize=%g A heff=%g A\n",Rmin, Rmax, BinSize, heff);
	PrintLog(" Rmin=%g Rmax=%g BinSize=%g heff=%g\n",Rmin, Rmax, BinSize, heff);
	aptr->Print_info(cout,1);
	
	float Cbulk=elmod->rionst;
	
	for(i=0;i<2;i++)
	{
		HaField3D* conc;
		if(i==0)conc=&(elmod->ConcMap0);
		else conc=&(elmod->ConcMap1);
		
		HaVec_double m_gcont(nbins*2);
		HaVec_int m_count(nbins*2);
		
		float BoxMin[3],BoxMax[3];
		BoxMin[0] = conc->GetXmin();
		BoxMin[1] = conc->GetYmin();
		BoxMin[2] = conc->GetZmin();
		BoxMax[0] = conc->GetXmax();
		BoxMax[1] = conc->GetYmax();
		BoxMax[2] = conc->GetZmax();
		
		PrintLog("BoxMin=[%g %g %g]A3 BoxMax=[%g %g %g]A3 \n", BoxMin[0], BoxMin[1], BoxMin[2], BoxMax[0], BoxMax[1], BoxMax[2]);
		
		int k;
		float fx,fy,fz;
		float start[3];
		float end[3];
		float ax,ay,az;
		ax=Ratom[0];
		ay=Ratom[1];
		az=Ratom[2];
		
		PrintLog("atom=[%g %g %g]\n", Ratom[0], Ratom[1], Ratom[2]);
		float BinSizeRev=1.0/BinSize;
		
		int nsteps=1;
		double dsteps=1.0;
		for(k=0;k<3;k++)
		{
			start[k]=Ratom[k]-Rmax;
			end[k]=Ratom[k]+Rmax;
			if(start[k]<BoxMin[k])start[k]=BoxMin[k];
			if(end[k]>BoxMax[k])end[k]=BoxMax[k];
			nsteps*=int((end[k]-start[k])/heff);
			dsteps*=(end[k]-start[k])/heff;
		}
		PrintLog("Have to make %d cycles (dsteps=%g)\n",nsteps,dsteps);
		PrintLog("start=[%g %g %g]A3 end=[%g %g %g]A3 \n", start[0], start[1], start[2], end[0], end[1], end[2]);
		PrintLog("start=[%g %g %g] end=[%g %g %g]\n", start[0], start[1], start[2], end[0], end[1], end[2]);
		PrintLog("BinSizeRev=%g\n",BinSizeRev);
		
		m_gcont.set(0.0);
		m_count.set(0);
		
		int count=0;
		for(fx=start[0];fx<=end[0];fx+=heff)
			for(fy=start[1];fy<=end[1];fy+=heff)
				for(fz=start[2];fz<=end[2];fz+=heff)
		{
			float r;
			r=sqrt((ax-fx)*(ax-fx)+(ay-fy)*(ay-fy)+(az-fz)*(az-fz))*BinSizeRev;
			m_gcont[(int)r]+=conc->GetInterpolValAtPoint(fx,fy,fz)/Cbulk;
			m_count[(int)r]++;
			count++;
			if(count%500000==0)PrintLog("finished %d cycles\n",count);
		}
		for(k=0;k<nbins;k++)
		{
			if(m_count[k]>0)
				gcont[i][k]=m_gcont[k]/m_count[k];
			else
				gcont[i][k]=0.0;
		}
	}*/
	return EXIT_SUCCESS;
}
int ElModRadDist::Setheff(float _heff)
{
	heff=_heff;
	return EXIT_SUCCESS;
}
int ElModRadDist::SetMinMax(float min,float max,float binsize)
{
	Rmin=min;
	Rmax=max;
	BinSize=binsize;
	nbins=roundf((max-min)/binsize);
	nbins++;
	//Rleft.resize(nbins+1);
	//int i;
	//for(i=0;i<nbins+1;i++)
	//{
	//	Rleft[i]=Rmin+i*binsize;
	//}
	return EXIT_SUCCESS;
}
HaVec_double* ElModRadDist::Getgcont(int ion,int atm)
{
  return gcont[ion][atm];
}
int ElModRadDist::PrintRadDist(const char *filename)
{
	PrintLog("ElModRadDist::PrintRadDist to %s\n",filename);
	int i,j;
	float BinSizeH=BinSize/2.0;
	float r;
	FILE *out=fopen(filename,"w");
	
	
	fprintf(out,"r_A");
	HaAtom* aptr;
	AtomIteratorMolSet aitr(pmset);
	for(i=0;i<AtomsList.size();i++)
	{
		j=1;
		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if(j==AtomsList[i])
			{
				fprintf(out," g%s_%d-K w%s_%d-K g%s_%d-Cl w%s_%d-Cl", aptr->GetName(), AtomsList[i], aptr->GetName(), AtomsList[i], aptr->GetName(), AtomsList[i], aptr->GetName(), AtomsList[i]);
				break;
			}
			j++;
		}
	}
	fprintf(out,"\n");
	
	for(i=0;i<nbins;i++)
	{
		r= Rmin+BinSize*i;
		fprintf(out,"%.5g",r);
		for(j=0;j<AtomsList.size();j++)
		{
			float wK=11.0,wCl=11.0;
			if(gcont[0][j]->GetVal_idx0(i)>0.0)
				wK=-log(gcont[0][j]->GetVal_idx0(i));
			if(gcont[1][j]->GetVal_idx0(i)>0.0)
				wCl=-log(gcont[1][j]->GetVal_idx0(i));
			fprintf(out," %.4g %.4g %.4g %.4g", gcont[0][j]->GetVal_idx0(i), wK, gcont[1][j]->GetVal_idx0(i), wCl);
		}
		fprintf(out,"\n",r);
	}
	fclose(out);
	return EXIT_SUCCESS;
}
#endif
