//
// C++ Interface: buildworldni
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef BUILDWORLDNI_H
#define BUILDWORLDNI_H

#ifdef SWIG
#include <Python.h>
#endif

#include <vector>

#include "pnps.h"
#include "math.h"

#define BIGDISTANSE 10000
class ContWorld;
class BuildWorldNI;
class NodeIndexing;

#include "pnpstructs.h"

#ifdef HARLEM_MOD
class HaMolSet;
#endif

#include <string>

//!PartialWorldForSAS for memory limited search of SAS
typedef struct
{
	int GS_X;
	int GS_Y;
	int GS_Z;
	int GS_XY;
	int locGS_Z;
	int locR0_Z;
	
	float *Surf;
	int *Field;
	float3 *surf_points;
	int Nsurf_points;
	
}PartialWorldForSAS;

class GenericGeometricalObject:public PnpsObject
{
	public:
		GenericGeometricalObject();
		GenericGeometricalObject(const GenericGeometricalObject *GGO);
		~GenericGeometricalObject();
		
		friend class ElMod;
		
		virtual int Copy(const GenericGeometricalObject *GGO);
		virtual int InitZero();//!<initialize zero values
		virtual int Clear();//!< delete all internal arrays and objects
	public:
		int Epsilon;///< Index of Dielectric Constant
		int NIonsTypes;///< here is dim[IonDiffusion]
		int *IonsD;///< Diffusion Indexing
		float *C;///< Concentration of ions
	public:
		virtual int setPar(int iEps, int iD0, int iD1,float C0,float C1=-1.0);
		virtual int Print(BuildWorldNI *buildWorld=NULL);

#ifdef SWIG
%rename (cxxSetParam) setParam;
		%pythoncode %{
		def setParam(self,iEps, iD, C, NIonsTypes=None):
			"""
			set parameters
			Input Parameters:
				iEps = int
				iD = [float]*NIonsTypes
				C = float or [float]*NIonsTypes
			"""
			m_NIonsTypes=NIonsTypes
			if m_NIonsTypes!=None:
				self.NIonsTypes=m_NIonsTypes
			if m_NIonsTypes==None:
				if not (type(iD) is list or type(iD) is tuple):
					print "Error: iD should be list or tuple"
					return None
				m_NIonsTypes=len(iD)
			
			m_iD=new_intArray(self.NIonsTypes)
			if type(iD) is list or type(iD) is tuple:
				for i in range(self.NIonsTypes):intArray_setitem(m_iD,i,iD[i])
			else:
				for i in range(self.NIonsTypes):intArray_setitem(m_iD,i,iD)
			m_C = new_floatArray(self.NIonsTypes)
			if type(C) is list or type(C) is tuple:
				for i in range(self.NIonsTypes):floatArray_setitem(m_C,i,C[i])
			else:
				for i in range(self.NIonsTypes):floatArray_setitem(m_C,i,C)
			self.cxxSetParam(iEps,m_NIonsTypes,m_iD,m_C)
			delete_intArray(m_iD)
			delete_floatArray(m_C)
		%}
#endif
		virtual int setParam(int iEps,int m_NIonsTypes, int *iD,float *m_C);
		
		virtual int ChangeUnits(ContWorld* world,bool ToInternal);///<Change units for internal or external use, internal units are bound to grid parameters
		//!making premaps for dielectric constant and diffusion
		virtual int BuildPreMaps(BuildWorldNI *Builder,ContWorld* world,int *Field,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth,float *Surf,int ParmMask);
		virtual int BuildDistMapsFromAtomsCenter(BuildWorldNI *Builder,ContWorld* world,float *Dist,float Rmax,float DoNotConsiderAtomsSmallerThen);
		virtual int BuildDistMaps(BuildWorldNI *Builder,ContWorld* world,int *Field,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth,float *Surf);
		virtual int BuildPreMapsMemLim(PartialWorldForSAS *pw,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth);
		virtual int BuildPreMapsCharges(BuildWorldNI *Builder,ContWorld* world,float *Field);
		virtual int BuildLJRepultionMap(BuildWorldNI *Builder,ContWorld* world,int Ion,float *V,float LimitVlj);
		virtual int SetBoundaryCondition(BuildWorldNI *Builder,ContWorld* world);
		virtual int RotateGGO(double *n, double cosa, double sina);
		//!< rotate vector around unit vector n by angle with cosa and sina
		void RotateVecDouble(double *Rext,double *n, double cosa, double sina)
		{
			double Rout[3];
			Rout[0] = (n[0]*n[0] + (1-n[0]*n[0])*cosa)*Rext[0] +
					(n[0]*n[1]*(1-cosa)  -  n[2]*sina)*Rext[1] + 
					(n[0]*n[2]*(1-cosa)  +  n[1]*sina)*Rext[2];
			Rout[1] = (n[0]*n[1]*(1-cosa)  +  n[2]*sina)*Rext[0]+
					(n[1]*n[1] + (1-n[1]*n[1])*cosa)*Rext[1]+
					(n[1]*n[2]*(1-cosa)  -  n[0]*sina)*Rext[2];
			Rout[2] = (n[0]*n[2]*(1-cosa) - n[1]*sina)*Rext[0]+
					(n[1]*n[2]*(1-cosa) + n[0]*sina)*Rext[1]+
					(n[2]*n[2] + (1-n[2]*n[2])*cosa)*Rext[2];
			Rext[0]=Rout[0];
			Rext[1]=Rout[1];
			Rext[2]=Rout[2];
		}
		//!< rotate vector around unit vector n by angle with cosa and sina
		void RotateVecFloat(float *Rext,float *n, float cosa, float sina)
		{
			float Rout[3];
			Rout[0] = (n[0]*n[0] + (1-n[0]*n[0])*cosa)*Rext[0] +
					(n[0]*n[1]*(1-cosa)  -  n[2]*sina)*Rext[1] + 
					(n[0]*n[2]*(1-cosa)  +  n[1]*sina)*Rext[2];
			Rout[1] = (n[0]*n[1]*(1-cosa)  +  n[2]*sina)*Rext[0]+
					(n[1]*n[1] + (1-n[1]*n[1])*cosa)*Rext[1]+
					(n[1]*n[2]*(1-cosa)  -  n[0]*sina)*Rext[2];
			Rout[2] = (n[0]*n[2]*(1-cosa) - n[1]*sina)*Rext[0]+
					(n[1]*n[2]*(1-cosa) + n[0]*sina)*Rext[1]+
					(n[2]*n[2] + (1-n[2]*n[2])*cosa)*Rext[2];
			Rext[0]=Rout[0];
			Rext[1]=Rout[1];
			Rext[2]=Rout[2];
		}
		int ID_GridBegin;//!first ID of GGO for grid building
		int ID_GridEnd;//!last ID of GGO for grid building
		virtual int SetGridIDs(int startIDfrom);//!took grid ID, return next available number
		virtual float GetDistToSurf(float fx,float fy,float fz, float Ri, float Rpr,int GridID)
		{return 0.0;}
		virtual float GetDistSQToUnitCenter(float fx,float fy,float fz, float Ri, float Rpr,int GridID)
		{return 0.0;}
		virtual void GetClosestSurfPoint(float *fx,float *fy,float *fz, float Ri, float Rpr,int GridID)
		{}
	protected:
		ContWorld* old_world;
		bool InternalUnit;
};
class GOAtoms:public GenericGeometricalObject
{
	public:
		GOAtoms();
		GOAtoms(const GOAtoms * goAtoms);
		~GOAtoms();
		friend class ElMod;
		virtual int Copy(const GOAtoms *goAtoms);
		virtual int InitZero();//!<initialize zero values
		virtual int Clear();//!< delete all internal arrays and objects
		virtual int Print(BuildWorldNI *buildWorld=NULL);
		friend class PMFCalculation;
	public:
		enum{ClosestNode=0,Linear=1,Cone=2};
	public:
		int NAtoms;
		std::vector<float> *r_x,*r_y,*r_z;//!<pointers to r[] for python interfae
		std::vector<float> r[3];
		std::vector<float> R,q;
		std::vector<float> HalfSigma,FourEpsilon;
		std::vector< std::vector<double> > LJA,LJB;
		std::vector< std::vector<double> > SRA,SRN;
		std::vector< std::vector<float> > IER;
		
		int DefineChargeDist(const char *str);
		int ChargeDist;
		int ChargeDistN;
		float Offset[3];
		bool MakePreRoll;
	public:
		virtual int setPar(int iEps, int iD0, int iD1,float C0,float C1=-1.0);
		//!Save atoms in pqr format
		int SavePQR(const char *filename);
		//!Load atoms from pqr format
		int LoadPQR(const char *filename);
		int LoadPQRonlyQR(const char *filename);
		int LoadPRQ(const char *filename);
		//!Load atoms from DELPHI input file
		int LoadDelphiStyle(const char *PDB,const char *CRG,const char *SIZ,bool center);
		//!Load atoms from pre format (AMBER FF for VdW)
		int LoadPRE(const char *filename);
		//!Load param for atoms from pab format (A,B coef for LJ)
		int LoadPAB(const char *filename);
		//!Load param for ion exclusion
		int LoadIER(const char *filename);
		//!Load param for soft repultion
		int LoadPAN(const char *filename);
		int LoadAtomicParam(const char *filename);
#ifdef HARLEM_MOD
		int LoadHaMolSet(HaMolSet *molset);
#endif

		int IAVMethod;//!0-HWR0(FF+IonSize), 1-HWR, 2 - LJ, 3 - AN
		
		int SetCenterToOrigin();
		int SetCenterToOriginWithRad();
		
		virtual int ChangeUnits(ContWorld* world,bool ToInternal);///<Change units for internal or external use, internal units are bound to grid parameters
		
		virtual int BuildPreMaps(BuildWorldNI *Builder,ContWorld* world,int *Field,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth,float *Surf,int ParmMask);
		virtual int BuildPreMapsMemLim(PartialWorldForSAS *pw,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth);
		virtual int BuildDistMapsFromAtomsCenter(BuildWorldNI *Builder,ContWorld* world,float *Dist,float Rmax,float DoNotConsiderAtomsSmallerThen);
		virtual int BuildPreMapsCharges(BuildWorldNI *Builder,ContWorld* world,float *Field);
		virtual int BuildLJRepultionMap(BuildWorldNI *Builder,ContWorld* world,int Ion,float *V,float LimitVlj);
		virtual int SetBoundaryCondition(BuildWorldNI *Builder,ContWorld* world);
		virtual int SetBoundaryConditionOnCuda(BuildWorldNI *Builder,ContWorld* world);
		virtual int RotateGGO(double *n, double cosa, double sina);
		virtual int SetGridIDs(int startIDfrom);
		virtual float GetDistToSurf(float fx,float fy,float fz, float Ri, float Rpr, int GridID)
		{
			int iatom=GridID-ID_GridBegin;
			float d = (r[0][iatom]-fx)*(r[0][iatom]-fx) + (r[1][iatom]-fy)*(r[1][iatom]-fy) + (r[2][iatom]-fz)*(r[2][iatom]-fz);
			return sqrt(d)-Ri-R[iatom];
		}
		virtual float GetDistSQToUnitCenter(float fx,float fy,float fz, float Ri, float Rpr,int GridID)
		{
			int iatom=GridID-ID_GridBegin;
			return (r[0][iatom]-fx)*(r[0][iatom]-fx) + (r[1][iatom]-fy)*(r[1][iatom]-fy) + (r[2][iatom]-fz)*(r[2][iatom]-fz);
		}
		virtual void GetClosestSurfPoint(float *fx,float *fy,float *fz, float Ri, float Rpr,int GridID)
		{
			int iatom=GridID-ID_GridBegin;
			float d = (r[0][iatom]-*fx)*(r[0][iatom]-*fx) + (r[1][iatom]-*fy)*(r[1][iatom]-*fy) + (r[2][iatom]-*fz)*(r[2][iatom]-*fz);
			float f=(R[iatom]+Rpr+Ri)/sqrt(d);
			
			*fx=r[0][iatom] + f*(*fx-r[0][iatom]);
			*fy=r[1][iatom] + f*(*fy-r[1][iatom]);
			*fz=r[2][iatom] + f*(*fz-r[2][iatom]);
		}
	private:
		int BuildPreMapsChargesClosestNode(BuildWorldNI *Builder,ContWorld* world,float *Field);
		int BuildPreMapsChargesLinear1(BuildWorldNI *Builder,ContWorld* world,float *Field);
		int BuildPreMapsChargesLinearN(BuildWorldNI *Builder,ContWorld* world,float *Field,int N);
		int BuildPreMapsChargesCone(BuildWorldNI *Builder,ContWorld* world,float *Field,int N);
	public:
		std::vector<std::string> ChargeDistStrs;
};
class GOTube:public GenericGeometricalObject
{
	public:
		GOTube();
		~GOTube();
		friend class BuildWorldEu;
		virtual int InitZero();//!<initialize zero values
		virtual int Clear();//!< delete all internal arrays and objects
		virtual int Print(BuildWorldNI *buildWorld=NULL);
	public:
		float XY[2], Z[2], R[2];
		float HalfSigma,FourEpsilon;
	public:
		int setTubeParam(float X,float Y,float Z0,float Z1,float R0,float R1);
		
		virtual int ChangeUnits(ContWorld* world,bool ToInternal);///<Change units for internal or external use, internal units are bound to grid parameters
		virtual int BuildPreMaps(BuildWorldNI *Builder,ContWorld* world,int *Field,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth,float *Surf,int ParmMask);
		virtual int BuildPreMapsMemLim(PartialWorldForSAS *pw,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth);
		virtual int BuildDistMapsFromAtomsCenter(BuildWorldNI *Builder,ContWorld* world,float *Dist,float Rmax,float DoNotConsiderAtomsSmallerThen);
		virtual int BuildLJRepultionMap(BuildWorldNI *Builder,ContWorld* world,int Ion,float *V,float LimitVlj);
		
		virtual int RotateGGO(double *n, double cosa, double sina);
		virtual float GetDistToSurf(float fx,float fy,float fz, float Ri, float Rpr, int GridID)
		{
			float d;
			if((tZ1mRpr<fz)&&(fz<tZ1))
			{
				d=tZ1-fz;
			}
			else if((tZ1<=fz)&&(fz<=tZ2))
			{
				d=100000.0;
			}
			else if((tZ2<fz)&&(fz<tZ2pRpr))
			{
				d=fz-tZ2;
			}
			return d;
		}
		virtual void GetClosestSurfPoint(float *fx,float *fy,float *fz, float Ri, float Rpr,int GridID)
		{
			if((tZ1mRpr<*fz)&&(*fz<tZ1))
			{
				*fz=tZ1mRpr;
			}
			else if((tZ1<=*fz)&&(*fz<=tZ2))
			{
				*fz=-100000.0;
			}
			else if((tZ2<*fz)&&(*fz<tZ2pRpr))
			{
				*fz=tZ2pRpr;
			}
			//*fx=r[0][iatom] + f*(*fx-r[0][iatom]);
			//*fy=r[1][iatom] + f*(*fy-r[1][iatom]);
			//*fz=r[2][iatom] + f*(*fz-r[2][iatom]);
		}
	protected:
		float tX;
		float tY;
		float tZ1;
		float tZ2;
		float tR1;
		float tR2;
		float tZ1mRpr;
		float tZ2pRpr;
};
class GODonut:public GenericGeometricalObject
{
	public:
		GODonut();
		GODonut(float x, float y,float z,float R,float r,float Scale);
		~GODonut();
		virtual int InitZero();//!<initialize zero values
		virtual int Clear();//!< delete all internal arrays and objects
	public:
		float R,r;
		float X,Y,Z;
	public:		
		virtual int ChangeUnits(ContWorld* world,bool ToInternal);///<Change units for internal or external use, internal units are bound to grid parameters
		virtual int BuildPreMaps(BuildWorldNI *Builder,ContWorld* world,int *Field,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth,float *Surf,int ParmMask);
		virtual int BuildPreMapsMemLim(PartialWorldForSAS *pw,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth);
};
class GOMembraneZ:public GenericGeometricalObject
{
	public:
		GOMembraneZ();
		~GOMembraneZ();
		virtual int InitZero();//!<initialize zero values
		virtual int Clear();//!< delete all internal arrays and objects
		virtual int Print(BuildWorldNI *buildWorld=NULL);
	public:
		float Z[2];
	public:
		virtual int ChangeUnits(ContWorld* world,bool ToInternal);///<Change units for internal or external use, internal units are bound to grid parameters
		virtual int BuildPreMaps(BuildWorldNI *Builder,ContWorld* world,int *Field,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth,float *Surf,int ParmMask);
		virtual int BuildPreMapsMemLim(PartialWorldForSAS *pw,int iValue,int iBulkValue,float *Displ,float Rion,float Rsmoth);
};

//class PyObject;

/**
	Build Node Indexing of ContWorld
	@author Nikolay Simakov <nsimakov@andrew.cmu.edu>
*/
class BuildWorldNI:public PnpsObject
{
	public:
		BuildWorldNI();
		

		~BuildWorldNI();
		
		
		
		friend class PMFCalculation;//!<PMFCalculation use this class to build system
		friend class ElMod;
		friend class GOAtoms;
		virtual int InitZero();//!<initialize zero values
		virtual int Clear();//!< delete all internal arrays and objects
		virtual int Print();
	protected:
		
		GenericGeometricalObject* BulkParam;
		std::vector<GenericGeometricalObject*> GOElms;
	public:
		//Setters for Python Interfase
		//int setEpsilonValues(PyObject *EpsilonValues);
		int setEpsilonValues(float e0=0.0, float e1=0.0, float e2=0.0, float e3=0.0, \
			float e4=0.0, float e5=0.0, float e6=0.0, float e7=0.0, \
			float e8=0.0, float e9=0.0, float e10=0.0, float e11=0.0, \
			float e12=0.0, float e13=0.0, float e14=0.0);
		int setDiffusionValues(float D0=0.0, float D1=0.0, float D2=0.0, float D3=0.0, \
			float D4=0.0, float D5=0.0, float D6=0.0);
		//! get index of dielectric constant value, if this particular dielectric constant is not present add it and return the index 
		int getEpsilonValueIndex(float epsilon);
		//! get index of diffusion coeficient value, if this particular dielectric constant is not present add it and return the index 
		int getDiffusionValueIndex(float D0);
		int setIonsR(float R0, float R1);
		int setRwat(float m_Rwat);
		int setRsmooth(float m_Rsmooth);
		int setDiffusionMode(int DiffMode);
		int setExpDifPar(float m_Kexpdiff, float m_Rexpdiff, float m_DiffRmaxAffect);
		int setMakeDielectricMap(int bVal);
		int setMakeChargeMap(int bVal);
		int setMakeDiffusionMap(int bVal);
		int setMakeConcentrationMap(int bVal);
		int setMakeSoftRepultionMap(int bVal);
		int setBulkPar(int iEps, int iD0, int iD1,float C0,float C1=-1.0);
#ifdef SWIG
%rename (cxxSetBulkParam) setBulkParam;
%rename (cxxSetIonsRadii) setIonsRadii;
		%pythoncode %{
		def setBulkParam(self, iEps, iD, C):
			"""
			set parameters for bulk
			Input Parameters:
				iEps = int
				iD = [float]*NIonsTypes
				C = float or [float]*NIonsTypes
			"""
			if not (type(iD) is list or type(iD) is tuple):
				print "Error: iD should be list or tuple"
				return None
			if len(iD)!=self.NIonsTypes:
				print "Error: len(iD)!=self.NIonsTypes"
				return None
			m_iD=new_intArray(self.NIonsTypes)
			for i in range(self.NIonsTypes):intArray_setitem(m_iD,i,iD[i])
			m_C = new_floatArray(self.NIonsTypes)
			if type(C) is list or type(C) is tuple:
				for i in range(self.NIonsTypes):floatArray_setitem(m_C,i,C[i])
			else:
				for i in range(self.NIonsTypes):floatArray_setitem(m_C,i,C)
			self.cxxSetBulkParam(iEps,m_iD,m_C)
			delete_intArray(m_iD)
			delete_floatArray(m_C)
		def setIonsRadii(self, m_R):
			"""
			set parameters for bulk
			Input Parameters:
				m_R = [float]*NIonsTypes
			"""
			if len(m_R)!=self.NIonsTypes:
				print "Error: len(m_R)!=self.NIonsTypes"
				return None
			mm_R=new_floatArray(self.NIonsTypes)
			for i in range(self.NIonsTypes):floatArray_setitem(mm_R,i,m_R[i])
			self.cxxSetIonsRadii(iEps,mm_R)
			delete_floatArray(mm_R)
			%}
#endif
		int setBulkParam(int iEps, int* iD, float* C);
		int setIonsRadii(float *m_R);
		int addGGO(GenericGeometricalObject* ggo);
	public:
		
		float Rwat,Rsmooth;///<Solvent Radius, Smothing Radius
		int BoundaryCondition;///< Boundary Condition, is separeted from world one
		bool BldBCatPlane[3];
		int NIonsTypes; ///< Number of Ion's Types in the system now only 2 is working, should be the same as in ContWorld
		float *IonsR; ///< Charge of each ion type, now [1.0, -1.0] is tested
		float *IonsHalfSigma;///<as in AMBER
		float *IonsFourEpsilon;///<as in AMBER
		float *IonsLJA;///<A coef for LJ
		float *IonsLJB;///<B coef for LJ
		int DielNum;///<Number of Dielectric constant Values
		float *DielConst;///<Array of Dielectric constant Values
		int DiffusionNum;///<Number of Diffusion constant Values
		float *DiffusionConst;///<Array of Diffusion constant Values
		int DiffusionMode;///<Diffusion Modes Plain or Exp
		float Kexpdiff;///<DiffusionMode==Exp D=D0(1-exp(k*r-R))
		float Rexpdiff;
		float DiffRmaxAffect;
		int iDzero;///<for python interface, 0-based
		std::string SaveDistMap;
		
		bool MakeDielectricMap;
		bool MakeDiffusionMap;
		bool MakeConcentrationMap;
		bool MakeChargeMap;
		bool MakeLJRepultion;//!<Make maps with LJ potential, will store it in PMF
		bool AddElPotToPMF;
		
		bool RemovingCavities;
		bool RemovingCavitiesOnDielectricMap;
		int RemCavOnDielFillWith;
		int RemCavOnDielWhere2Look;
		
		float LimitVlj;//!<limit Vlj by use of LimitVlj*Sigma, if =1 then just repultion
		
		std::vector<std::string> BoundaryStr;
		enum {ZeroBC=0,CoulBC=1};
		
		bool BuildUsingGPU;
		int MemoryLimitOnOneArray;//!< Memory Limit on temporal arrays in MB (64 is a choice if <2GB totally)
	private:
		
	public:
		int ChangeUnits(ContWorld* world,bool ToInternal);///<Change units for internal or external use, internal units are bound to grid parameters

#ifdef SWIG
%rename (cxxBuildContWorld) BuildContWorld;
		%pythoncode %{
		def BuildContWorld(self,world,Verbose=True):
			"""Build continious world
			Input Parameters:
				world=ContWorld
					instance of ContWorld
				Verbose=bool, default=True
					if set will print BuildWorld configuration
			"""
			if Verbose:
				self.Print()
			self.cxxBuildContWorld(world)
		%}
#endif
		virtual int BuildContWorld(ContWorld* world);
		int BuildDielMapsOnCuda(ContWorld* world,NodeIndexing *NIndexingNew);
	protected:
		int BuildContWorldCharges(ContWorld* world,NodeIndexing* NIndexingNew);
	public:
		int BuildDiffExp(ContWorld* world);
		int BuildDiffExpOld(float K0, ContWorld* world);
		virtual int FinalazeSEV(ContWorld* world,int *Field,int iBulkValue,float Rsmooth,float *Surf);
		virtual int FinalazeSEV2(ContWorld* world,int *Field,int iBulkValue,float Rsmooth,float3 *surf_points, int Nsurf_points);
		int FinalazeSEVDiff(ContWorld* world,int *Field,int iBulkValue,float Rsmooth,float *Surf);
		int RemovingCavitiesList(int GSX,int GSY,int GSZ,int *Field,int iDzero);
		int RemovingCavitiesAtValues(int GSX,int GSY,int GSZ,int *Field, int iWhereToLook, int iFillWith);
	public:
		enum {Plain=0,Exp=1};
		std::vector<std::string> StdStrDiffusionModes;
	protected:
		bool InternalUnit;
		ContWorld* old_world;
#ifdef SWIG
	%pythoncode %{
		def setBulk(self,**kwargs):
			iEps = kwargs.get("iDielConst")
			if iEps==None:
				print "Error: iDielConst must be defined"
			
			iD = kwargs.get("iDiffCoef")
			if iD==None:
				print "Error: iDiffCoef must be defined"
			if len(iD)!=2:
				print "Error: iDiffCoef must have two values"
			
			C = kwargs.get("C",0.1)
			self.setBulkPar(iEps, iD[0], iD[1], C)
			
		def addAtoms(self,
				iDielConst=None,
				DielConst=2.0,
				MolSet=None,
				PQR=None,
			    PAN=None,
				AtomsPQR=None,
			    AtomicParam=None,
				**kwargs):
			"""
			add atoms to the system
			Input Parameters:
				DielConst=float, default=2.0
					dielectric constant
				MolSet=HaMolSet, default None
					Load coordinates, charges and radii from Harlem MolSet
					Note that if HaMolSet.p_save_opt_default.save_selected is set to 1
					only selected atoms will be loaded
				PQR=string, default None
					Load coordinates, charges and radii from PQR file
				PAN = string, default None
					Load coordinates, charges and radii and SR parameters from PAN file
				AtomicParam = string, default None
					Load coordinates, charges and radii and SR parameters from AtomicParam file
				AtomsPQR=list of list with x,y,z,q and r, default=None
					Load coordinates, charges and radii from list of values.
					For example:
						AtomsPQR=[[0.0,0.0,-1.0,-1.0,2.0],[0.0,0.0,1.0,1.0,2.0]]

			Returned value:
				GOAtoms - atoms container
			"""
			#print "BuildWorldNI::addAtoms()"
			
			if iDielConst==None:
				iDielConst=self.getEpsilonValueIndex(DielConst)
			if iDielConst==None:
				print "Error: iDielConst must be defined"
			
			a=GOAtoms()
			a.setParam(iDielConst, self.iDzero,0.0,NIonsTypes=self.NIonsTypes)
			
			a.ChargeDist=a.Linear;
			a.ChargeDistN=1;

			atomsLoaded=False
			
			if PQR!=None:
				a.LoadPQR(PQR)
				atomsLoaded=True
			if PAN is not None:
				a.LoadPAN(PAN)
				atomsLoaded=True
			if AtomicParam is not None:
				a.LoadAtomicParam(AtomicParam)
				atomsLoaded=True
			if MolSet!=None:
				a.LoadHaMolSet(MolSet)
				atomsLoaded=True
			if AtomsPQR!=None:
				for atom in AtomsPQR:
					x,y,z,q,r=atom
					a.r_x.push_back(x)
					a.r_y.push_back(y)
					a.r_z.push_back(z)
					a.q.push_back(q)
					a.R.push_back(r)
				a.NAtoms=a.R.size();
				atomsLoaded=True
			if atomsLoaded==False:
				print "Error: Can not load atoms parameters!"
			a.MakePreRoll=True
			#for key in kwargs:
			#	print "keyword arg: %s: %s" % (key, kwargs[key])
			self.addGGO(a)
			a.thisown = 0
			return a
		
		def addTube(self,
				iDielConst=None,
				DielConst=2.0,
				x=0.0,
				y=0.0,
				z=[-12.0,12.0],
				R=[4.0,1000.0]):
			"""
			add tube or membrane with cyllindric hole to the system
			Input Parameters:
				DielConst=float, default=2.0
					dielectric constant
				x=float, default=0.0
				y=float, default=0.0
					x,y - coordinate of central symmetry axle of tube
				z=[float,float], default=[-12.0,12.0]
					beginning and end of tube
				R=float or [float,float], default=4.0
					internal radius of cillindric hole in the membrane or internal and external radii of the tube
			Returned value:
				GOTube - tube container
			"""
			#print "BuildWorldNI::addTube()"
			
			if iDielConst==None:
				iDielConst=self.getEpsilonValueIndex(DielConst)
			if iDielConst==None:
				print "Error: iDielConst must be defined"
			
			t=GOTube()
			t.setParam(iDielConst, self.iDzero,0.0,NIonsTypes=self.NIonsTypes)
			
			if type(R) is list or type(R) is tuple:
				m_R=R
			else:
				m_R=[R,1000.0]

			t.setTubeParam(x,y,z[0],z[1],m_R[0],m_R[1]);
			self.addGGO(t)
			t.thisown = 0
			return t
		def addMembraneZ(self,
				iDielConst=None,
				DielConst=2.0,
				z=[-12.0,12.0]):
			"""
			add membrane slab
			Input Parameters:
				DielConst=float, default=2.0
					dielectric constant
				z=[float,float], default=[-12.0,12.0]
					beginning and end of tube
			Returned value:
				GOMembraneZ- tube container
			"""
			#print "BuildWorldNI::addTube()"
			
			if iDielConst==None:
				iDielConst=self.getEpsilonValueIndex(DielConst)
			if iDielConst==None:
				print "Error: iDielConst must be defined"
			
			mz=GOMembraneZ()
			mz.setParam(iDielConst, self.iDzero,0.0,NIonsTypes=self.NIonsTypes)
			
			if not (type(z) is list or type(z) is tuple):
				print "Error: z must be type of two values"

			floatArray_setitem(mz.Z,0,z[0])
			floatArray_setitem(mz.Z,1,z[1])

			self.addGGO(mz)
			mz.thisown = 0
			return mz
		%}
#endif
};
#ifdef SWIG
%pythoncode %{
def GetBuilderForPNPSR(**kwargs):
	"""
		depreciated
	GetBuilderForPNPSR - get system continuum representation builder

	Parameters
		MakeDielectricMap=Bool, default=True
			Create dielectric constant distribution map
		MakeChargeMap=Bool, default=True
			Create static charge distribution map
		MakeDiffusionMap=Bool, default=True
			Create diffusion coefficient distribution map
		MakeConcentrationMap=Bool, default=True
			Create mobile ions concentration distribution map
		MakeSoftRepultionMap=Bool, default=True
			Create soft repultion potential between static atoms and mobile ions
	DielConstValues = None
	DiffusionValues = None

	iDielConst=None #Example: 1 , i.e. DielConstValues[1]
	iDiffCoef=None #Example:[2,3]
	C=0.1,

	RprobeDiel = 1.4
	RprobeDiff = 0.5
	Rions      = [2.0,2.0]

	DiffusionModel = "Plane"  |"MD"
	Kexpdiff=0.62 #options for MD diffusion
	DiffRmaxAffect=10.0 #options for MD diffusion

	
	"""
	print kwargs
	DielConstValues = kwargs.get("DielConstValues",None)
	print DielConstValues
	if DielConstValues==None:
		print "Error: DielConstValues must be defined"
		return None
	DiffusionValues = kwargs.get("DiffusionValues",None)
	if DiffusionValues==None:
		print "Error: DiffusionValues must be defined"
		return None
	
	bld=BuildWorldNI()
	#Set epsilon values which will be used in ContWorld
	#In objects instead of real dielectric constant is used index to the element of this array
	a=[0.0]*15
	for i in xrange(0,len(DielConstValues)):
		a[i]=float(DielConstValues[i])
	bld.setEpsilonValues(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], a[13], a[14]);
	#Set diffusion values (Bulk + Zero) which will be used in ContWorld
	#In objects instead of real diffusion constant is used index to the element of this array
	a=[0.0]*7
	for i in xrange(0,len(DiffusionValues)):
		a[i]=float(DiffusionValues[i])
	bld.setDiffusionValues(a[0], a[1], a[2], a[3], a[4], a[5], a[6])


	
	#Rprobe for creation of dielectric maps
	RprobeDiel = kwargs.get("RprobeDiel",1.4)
	bld.setRwat(RprobeDiel)
	#Rsmoth for smothing IAV; since we use PNP-SR our initial guess will molecular surphace
	
	RprobeDiff = kwargs.get("RprobeDiff",0.5)
	bld.setRsmooth(RprobeDiff)
	#Radii of mobile ion, for creation of IAV; now we will use explicit IER
	Rions = kwargs.get("Rions",[2.0,2.0])
	bld.setIonsR(Rions[0],Rions[1]) # R0, R1
	
	#Type of distance depentend diffusion:
	# BuildWorldNI.Plain - Bulk diffusion throughout IAV
	# BuildWorldNI.Exp - Diffusion which exponentially decay upon aproatching protein atoms
	DiffusionModel = kwargs.get("DiffusionModel","Plane")
	if DiffusionModel=="Plane":
		bld.setDiffusionMode(BuildWorldNI.Plain)
	else:
		bld.setDiffusionMode(BuildWorldNI.Exp)
	# Parameters for exponential decay
	Kexpdiff = kwargs.get("Kexpdiff",0.62)
	Rexpdiff = kwargs.get("Rexpdiff",2.2)
	DiffRmaxAffect = getattr(kwargs,"DiffRmaxAffect",10.0)
	bld.setExpDifPar(Kexpdiff,Rexpdiff,DiffRmaxAffect)
	
	#following flags describe what maps to calculate
	MakeDielectricMap=kwargs.get("MakeDielectricMap",True)
	MakeChargeMap=kwargs.get("MakeChargeMap",True)
	MakeDiffusionMap=kwargs.get("MakeDiffusionMap",True)
	MakeConcentrationMap=kwargs.get("MakeConcentrationMap",True)
	MakeSoftRepultionMap=kwargs.get("MakeSoftRepultionMap",True)

	if MakeDielectricMap : MakeDielectricMap=1
	else : MakeDielectricMap=0

	if MakeChargeMap : MakeChargeMap=1
	else : MakeChargeMap=0

	if MakeDiffusionMap : MakeDiffusionMap=1
	else : MakeDiffusionMap=0

	if MakeConcentrationMap : MakeConcentrationMap=1
	else : MakeConcentrationMap=0

	if MakeSoftRepultionMap : MakeSoftRepultionMap=1
	else : MakeSoftRepultionMap=0

	bld.setMakeDielectricMap(MakeDielectricMap)
	bld.setMakeChargeMap(MakeChargeMap)
	bld.setMakeDiffusionMap(MakeDiffusionMap)
	bld.setMakeConcentrationMap(MakeConcentrationMap)
	bld.setMakeSoftRepultionMap(MakeSoftRepultionMap)

	iEps = kwargs.get("iDielConst",None)
	if iEps==None:
		print "Error: iDielConst must be defined"
			
	iD = kwargs.get("iDiffCoef",None)
	if iD==None:
		print "Error: iDiffCoef must be defined"
	if len(iD)!=2:
		print "Error: iDiffCoef must have two values"
			
	C = kwargs.get("C",0.1)
	bld.setBulkPar(iEps, iD[0], iD[1], C)
	
    #Why?
	bld.thisown = 0
	return bld
def GetWorldBuilder(
		DielConstBulk=80.0,
		DiffCoefBulk=[19.57,20.32],
		Cbulk=1.0,
		RprobeDiel = 1.4,
		RprobeDiff = 0.5,
		Rions      = [2.0,2.0],
		DiffusionModel = "Plane",
		KexpDiff=0.62,
		RexpDiff=2.2,
		DiffRmaxAffect=10.0,
		BoundaryCondition="ZeroBC",
		RemoveIonsCavities=True,
		MakeDielectricMap=True,
		MakeChargeMap=True,
		MakeDiffusionMap=True,
		MakeConcentrationMap=True,
		MakeSoftRepultionMap=False):
	"""
		 
	GetBuilderForPNPSR - get system continuum representation builder

	Input Parameters:
		Bulk properties setup:
			DielConstBulk=float, default=80.0,
				Dielectric constant of a bulk
			DiffCoefBulk=[float, float], dim=Number of mobile ions, default=[19.57,20.32], Units:[1E-6 cm^2/s]
				Diffusion coefficient of mobile ions in bulk. Units:[1E-6 cm^2/s].
				Default values corresponds for KCl
			Cbulk=float or [float, float], dim=Number of mobile ions,default=1.0, Units:[M]
				Concentration of mobile ions in bulk
			BoundaryCondition=string, "ZeroBC" or "CoulBC", default="ZeroBC"
				Boundary conditions
		Solvent and mobile ions radii:
			RprobeDiel=float, default=1.4, Units:[Angstrom]
				Solvent molecule radius for solvent accessible volume calculation,
				used in dielectic constant distribution calculation.
			RprobeDiff=float, default=0.5, Units:[Angstrom]
				Smoothing probe radius for smoothing ions accessible volume
			Rions=[float, float], dim=Number of mobile ions,default=[2.0,2.0], Units:[Angstrom]
				Mobile ions radii for ion accessible volume calculation
		Diffusion model and its parameters:
			DiffusionModel = "Plane" or "MD", default="Plane"
				Diffusion model
			Kexpdiff=float, default=0.62
				option for MD diffusion
			Rexpdiff=float, default=2.2, Units:[Angstrom]
				option for MD diffusion
			DiffRmaxAffect=float, default=10.0, Units:[Angstrom]
				option for MD diffusion
		Distribution map to building flags:
			MakeDielectricMap=Bool, default=True
				Create dielectric constant distribution map
			MakeChargeMap=Bool, default=True
				Create static charge distribution map
			MakeDiffusionMap=Bool, default=True
				Create diffusion coefficient distribution map
			MakeConcentrationMap=Bool, default=True
				Create mobile ions concentration distribution map
			MakeSoftRepultionMap=Bool, default=False
				Create soft repultion potential between static atoms and mobile ions
	    Other options:
	        RemoveIonsCavities=Bool, default=True remove mobile ions cavities
	Returned value:
			In case of successfull execution return BuildWorldNI instance, otherwise None

	"""
	bld=BuildWorldNI()

	#get number of mobile ions types
	if not (type(DiffCoefBulk) is list or type(DiffCoefBulk) is tuple):
		print "Error: DiffCoefBulk should be list or tuple!"
		return None

	NIonsTypes=len(DiffCoefBulk)
	if len(DiffCoefBulk)!=len(Rions):
		print "Error: Different number of mobile ions types: len(DiffCoefBulk)!=len(Rions)"
		return None
	if type(Cbulk) is list or type(Cbulk) is tuple:
		if len(DiffCoefBulk)!=len(Cbulk):
			print "Error: Different number of mobile ions types: len(DiffCoefBulk)!=len(Cbulk)"
			return None
	bld.NIonsTypes=NIonsTypes
	
	#Rprobe for creation of dielectric maps
	bld.setRwat(RprobeDiel)
	#Rsmoth for smothing IAV; since we use PNP-SR our initial guess will molecular surphace
	bld.setRsmooth(RprobeDiff)
	#Radii of mobile ion, for creation of IAV; now we will use explicit IER
	bld.setIonsR(Rions[0],Rions[1]) # R0, R1
	
	#Type of distance depentend diffusion:
	# BuildWorldNI.Plain - Bulk diffusion throughout IAV
	# BuildWorldNI.Exp - Diffusion which exponentially decay upon aproatching protein atoms
	if DiffusionModel=="Plane":
		bld.setDiffusionMode(BuildWorldNI.Plain)
	elif DiffusionModel=="MD":
		bld.setDiffusionMode(BuildWorldNI.Exp)
	else:
		print "Error: Unknown Diffusion Model!"
		return None
	# Parameters for exponential decay
	bld.setExpDifPar(KexpDiff,RexpDiff,DiffRmaxAffect)
	
	#following flags describe what maps to calculate
	if MakeDielectricMap : MakeDielectricMap=1
	else : MakeDielectricMap=0

	if MakeChargeMap : MakeChargeMap=1
	else : MakeChargeMap=0

	if MakeDiffusionMap : MakeDiffusionMap=1
	else : MakeDiffusionMap=0

	if MakeConcentrationMap : MakeConcentrationMap=1
	else : MakeConcentrationMap=0

	if MakeSoftRepultionMap : MakeSoftRepultionMap=1
	else : MakeSoftRepultionMap=0

	bld.setMakeDielectricMap(MakeDielectricMap)
	bld.setMakeChargeMap(MakeChargeMap)
	bld.setMakeDiffusionMap(MakeDiffusionMap)
	bld.setMakeConcentrationMap(MakeConcentrationMap)
	bld.setMakeSoftRepultionMap(MakeSoftRepultionMap)

	iEps = bld.getEpsilonValueIndex(DielConstBulk)
	
	iD = []
	for iDiff in DiffCoefBulk:
		iD.append(bld.getDiffusionValueIndex(iDiff))

	bld.setBulkParam(iEps, iD, Cbulk)

	if BoundaryCondition=="ZeroBC":
		bld.BoundaryCondition=0
	elif BoundaryCondition=="CoulBC":
		bld.BoundaryCondition=1
	else:
		bld.BoundaryCondition=0
		print "Error: Unknown Boundary Condition! Should be ZeroBC or CoulBC"
		return None
	if RemoveIonsCavities:
		bld.RemovingCavities=True
	else:
		bld.RemovingCavities=False
    #Why?
	bld.thisown = 0
	return bld
%}
#endif
class BuildWorldEu:public BuildWorldNI
{
	public:
		BuildWorldEu();

		~BuildWorldEu();
		
		friend class PMFCalculation;//!<PMFCalculation use this class to build system
		friend class ElMod;
		virtual int InitZero();//!<initialize zero values
		virtual int Clear();//!< delete all internal arrays and objects
		virtual int BuildContWorld(ContWorld* world);
		virtual int BuildPreMapsAtoms(GOAtoms *goatoms,ContWorld* world,int **iEpsilon,float Rpr);
		virtual int BuildPreMapsTube(GOTube *gotube,ContWorld* world,int **iEpsilon,float Rpr);
		virtual int FinalazeSEV(ContWorld* world,int **iEpsilon);
		virtual int FinalazeSEV2(ContWorld* world,int **iEpsilon);
	protected:
		int TotalGridIds;
		int * DielIdx;
		GenericGeometricalObject **GGOIdx;
		int DoNotChangeEpsFlag;
};
class BuildWorldCmp:public BuildWorldNI
{
	public:
		BuildWorldCmp();

		~BuildWorldCmp();
		
		friend class PMFCalculation;//!<PMFCalculation use this class to build system
		friend class ElMod;
		virtual int InitZero();//!<initialize zero values
		virtual int Clear();//!< delete all internal arrays and objects
		virtual int BuildContWorld(ContWorld* world);
		virtual int BuildPreMapsAtoms(GOAtoms *goatoms,ContWorld* world,int **iEpsilon,float Rpr);
		virtual int FinalazeSEV(ContWorld* world,int **iEpsilon);
	protected:
		int TotalGridIds;
		int * DielIdx;
		GenericGeometricalObject **GGOIdx;
};

/**
	Build Maps using scaled algorithm
	@author Nikolay Simakov <nsimakov@andrew.cmu.edu>
 */
class BuildWorldScaled:public PnpsObject
{
	public:
		BuildWorldScaled();
		~BuildWorldScaled();
		
		virtual int InitZero();//!<initialize zero values
		virtual int Clear();//!< delete all internal arrays and objects
		
		virtual int BuildContWorld(ContWorld* world);
		int BuildContWorldCharges(ContWorld* world,NodeIndexing* NIndexingNew);
		int BuildContWorldDielMaps(ContWorld* world,NodeIndexing* NIndexingNew);
		
		float Rwat;
		int DielNum;///<Number of Dielectric constant Values
		float *DielConst;///<Array of Dielectric constant Values
		int BoundaryCondition;
		
		GenericGeometricalObject* BulkParam;
		std::vector<GenericGeometricalObject*> GOElms;
};

/**
	Dielectric constant and diffusion coefficient maps patcher
 */
class DielDiffMapsPatcher:public PnpsObject
{
	public:
		DielDiffMapsPatcher(ContWorld* world);
		~DielDiffMapsPatcher();

		ContWorld* m_ContWorld;

		int *eps[3];
		//!Convert the world->NIndex Diel Maps to Int Maps
		void GetIntDielMaps();
		//!Push the changes Back
		void PushIntDielMaps();
		//!PatchDielMaps
		void PatchDielMaps(int epsToOver,int epsNew,float x, float y, float z1, float z2,float R);

		float *refDiff[4];
		//!Copy current diffusion to ref and initiate new diffusion maps from it
		void InitNewDiff();
		//!Delete refDiff
        void DelRefDiff();
		//!PatchDiffMaps
        void PatchDiffMaps(float Dscale1K, float Dscale2K, float Dscale1Cl, float Dscale2Cl,float x, float y, float z1, float z2,float R);
#ifdef SWIG
	%pythoncode %{
		def setDielConstInUse(self,DielConstInUse):
			for i in range(14):
				self.m_ContWorld.NIndexing.SetDielConstInUse(i,DielConstInUse[i])
		%}
#endif
};

#endif
