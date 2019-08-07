//
// C++ Interface: contworld
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef CONTWORLD_H
#define CONTWORLD_H

#include "pnps.h"
#include <vector>
#include <string>
#include <zlib.h>

// #define DIELECTRIC_CONSTANT_0 1000
// #define DIELECTRIC_CONSTANT_1 1001
// #define DIELECTRIC_CONSTANT_2 1002
// #define DIELECTRIC_CONSTANT 1010
// #define STATIC_CHARGE 1003
// #define DYNAMIC_CHARGE 1009
// #define POTENTIAL 1007
// #define DIFFUSION 1008
// #define RFMAP 1011
// #define TYPE_OF_ARRAY 4
// #define ALL_ARRAYS_WAS_SENDED 5
#define DO_ONE 2001
#define DO_TWO 2002
#define DO_THREE 2003
#define DO_FOUR 2004
#define DO_FIVE 2005
#define DO_SIX 2006
#define SEND_ENERGY 3006
// #define FILE_MAP_GZ 4001///< map gz pnp_c style
// #define FILE_MAP_MBN 4002///< map binary
// #define FILE_MAP_MB2 4003///< map binary 2 coloumn style
// #define FILE_MAP_IGB 4004///< map binary igb style for flounder

#ifdef HARLEM_MOD
class HaField3D;
//class HaVec_float;
#endif
class VectorField3D;
class HaVec_float;

class MapsIOData : public PnpsObject
{
	public:
		MapsIOData();
		~MapsIOData();
	public:
		std::string DielectricConstantMap;
		std::string DiffusionMapFile;
		std::string DynamicChargeMapFile;
		std::string StaticChargeMapFile;
		std::string PotentialMapFile;
		std::string PMFMapFile;
		std::string PMFMapFile2;
		std::string NodeIndexingFile;
		int Mode;
	public:
		enum {Read=0,Write=1,Auto=2};
		std::vector<std::string> ModeStr;
		bool AddGroupStringSuffix;
		bool AddGroupNumberSuffix;
};

typedef unsigned int NodeIndex;
const int NodeIndexMaxValues = 16;
const int MaxIonTypes = 4;

class ContWorld;
class VectorField3D;
/*! 
NodeIndexing is class for fast and light storing dielectrical,diffusional and concentration fields.
will use instead of float* becouse good for solver

Field masks for indexed maps, for 4 types of ions (up to 16 valuse)
maps<int,float> Epsilon (15 - variable dielectric const not impliment yet)
maps<int,float,float> Diffusion,Concentration (presumably 0 for 0, 15 for variable diffusion not impliment yet)
int=4 byte
|76543210|76543210|76543210|76543210|

|00000000|00000000|00001111|00000000|- Epsilon0, 00 00 0F 00
|00000000|00000000|11110000|00000000|- Epsilon1, 00 00 F0 00
|00000000|00001111|00000000|00000000|- Epsilon2, 00 0F 00 00
|00000000|01110000|00000000|00000000|- Ion0,		 00 70 00 00
|00000011|10000000|00000000|00000000|- Ion1,		 03 80 00 00
|00011100|00000000|00000000|00000000|- Ion2,		 1C 00 00 00
|11100000|00000000|00000000|00000000|- Ion3,		 E0 00 00 00
*/
class NodeIndexing
{
	public:
		NodeIndexing();
		~NodeIndexing();
		
		int InitZero();
	public:
		//!each field is of size 4bite
		enum NodeIndexDescriptor
		{
			Epsilon0=0x00000F00,
			Epsilon1=0x0000F000,
			Epsilon2=0x000F0000,
			Ion0=		0x00700000,
			Ion1=		0x03800000,
			Ion2=		0x1C000000,
			Ion3=		0xE0000000,
			Epsilon0Sft=8,
			Epsilon1Sft=12,
			Epsilon2Sft=16,
			Ion0Sft=		20,
			Ion1Sft=		23,
			Ion2Sft=		26,
			Ion3Sft=		29,
			Flexible=	 0xF,
			
			BlackAndWhiteMask=		0x00000001,//0000 0001
			ChargeMask=					 0x00000002,//0000 0010
			DielBoarderMask=			0x00000004,//0000 0100
			DiffIon0BoarderMask=	0x00000008,//0000 1000
			DiffIon1BoarderMask=	0x00000010,//0001 0000
			DiffIon2BoarderMask=	0x00000020,//0010 0000
			DiffIon3BoarderMask=	0x00000040,//0100 0000
			BlackAndWhiteMaskSft=0,
			ChargeMaskSft=			 1,
			DielBoarderMaskSft=	2,
			DiffIon0BoarderMaskSft=	3,
			DiffIon1BoarderMaskSft=	4,
			DiffIon2BoarderMaskSft=	5,
			DiffIon3BoarderMaskSft=	6,
			
			DielConst=1,
			DiffConst=2,
			Conc=3,
			Charge=4
		};
//		enum 
//		{
			//QFormatStraightSequence=11,
			//QFormatBlackNWhiteSequences=12,
//		};
		NodeIndex EpsilonField[MaxIonTypes];
		NodeIndex EpsilonFieldSft[MaxIonTypes];
		NodeIndex IonField[MaxIonTypes];
		NodeIndex IonFieldSft[MaxIonTypes];
		NodeIndex DiffBoarderMask[MaxIonTypes];
		NodeIndex DiffBoarderMaskSft[MaxIonTypes];
		
		
		unsigned int GridSize[3];
		float GridScale;

		float Eps[NodeIndexMaxValues];
		float D[NodeIndexMaxValues],C[NodeIndexMaxValues];
		int NIonsTypes;
		float *IonsQ; ///< Charge of each ion type, now [1.0, -1.0] is tested
		
		NodeIndex* NIndex;//!<properties of node epsilon, diffusion concentration(initial or for PB)
		//int QFormat;
		float *Q;//!<value of charge
		//int ChargeNum[3];//!<ChargeNum[0]=0,ChargeNum[1] - white block started,ChargeNum[2]white block finished. when QFormat==QFormatBlackNWhiteSequences ChargeNum[1]=0,ChargeNum[2]=size(Q)/size(float)
		//int ConverQtoQFormatBlackNWhiteSequences();
		//int ConverQtoQFormatStraightSequence();
		bool PBC[3];//!<Periodic Boundary Condition, GridSize[3] shows current state if there is PBC on axes X then GridSize[0] will encrease by 2, saving to file is always in original GridSize

		int SetNoPBC();//!<set Periodic Boundary Condition to [false, false, false]
		int SetPBC(bool pbcX,bool pbcY,bool pbcZ);
		
		//vars for internal usage
		int QNum;
	protected:
		NodeIndexDescriptor GetShtFromMask(NodeIndex Mask) const;
	public:
		int GetNx() const { return GridSize[0]; } //!< return number of grid points along X axes 
		int GetNy() const { return GridSize[1]; } //!< return number of grid points along Y axes 
		int GetNz() const { return GridSize[2]; } //!< return number of grid points along Z axes 
		//! Return Diel Cons of index i in external units (for python)
		float GetDielConstInUse(int i);
		//! Set Diel Cons of for i using external units as input (for python)
		void SetDielConstInUse(int i, float Val);
		int SetNNodes(unsigned int *gridsize,float gridscale);
		int SetBlackAndWhiteNodes(int FirstNode=0);
		int SetIndexFieldFromIntArray(int *arr,NodeIndexDescriptor mask,NodeIndexDescriptor sft);
		int SetIonAccess(float **DiffusionsMaps);
		float* GetCMap(NodeIndexDescriptor FieldType, NodeIndex mask);
		int* GetIntArrayFromIndexField(NodeIndexDescriptor FieldType, NodeIndex mask);
		float* GetField(NodeIndexDescriptor FieldType, NodeIndexDescriptor mask, float *Map = NULL);
		#ifdef HARLEM_MOD
		HaField3D* GetHaField3D(NodeIndexDescriptor FieldType, NodeIndexDescriptor mask);
		#endif
		int CalcDielBoarder();
		int CalcDiffBoarder();
		
		int InsertSphereInDielMap(float *r, float R, int Value);
		
		//since the border is not in the cycle we do not need charged points on the boarde, couse it will lead to index switches in Poisson Part 
		int CheckBoarder4Q();
		int CheckNodeIndex();
		int RemoveBadDiffusionPoints();
		
		int SetChargeMapFromArray(float *q);
		int SetChargeMapFromArrayNoQonBoarder(float *q);
		int SetChargeMapToZero();
		float* GetChargeArray();
		int* GetChargeIndex();
		int ReadFromFile(const char *filename);
		int ConvertToPBC(NodeIndex *V,int GS_X,int GS_Y,int GS_Z,bool pbcX,bool pbcY,bool pbcZ);
		int ReadFromFileAddPBC(const char *filename,bool pbcX,bool pbcY,bool pbcZ);
		int ReadFromFile4SDMPI(const char *filename,ContWorld *World);
		//!will write index table in semi XML way
		int WriteToFile(const char *filename);
		int GetMyPart4MPI(ContWorld *World,NodeIndexing* NIndexGlobal=NULL);//!< Get part of system for parallel processing NIndexGlobal shouldnt be NULL for master node
		int GetCentralPartFromNIExtMyPart4MPI(NodeIndexing* NIExt);
		//int CollectGlobalSysFromMPIParts(CWorld *World,NodeIndexing* NIndexGlobal=NULL);
		inline int GetDiel(int i,int node)
		{
			int Val;
			switch(i)
			{
				case 0:
					Val=(NIndex[node]&Epsilon0)>>Epsilon0Sft;
					return Val;
				case 1:
					Val=(NIndex[node]&Epsilon1)>>Epsilon1Sft;
					return Val;
				case 2:
					Val=(NIndex[node]&Epsilon2)>>Epsilon2Sft;
					return Val;
				default:
					Val=(NIndex[node]&Epsilon0)>>Epsilon0Sft;
					return Val;
			}
		}
		inline float GetDielFloat(int i,int node)
		{
			return Eps[GetDiel(i,node)];
		}
		int SetSCharDielMapForGAPS(signed char **CDiel,float *gapsEps);
		inline void SetDiel(int i,int node,int Val)
		{
			unsigned int ClearMask;
			NodeIndex Cur;
			
			switch(i)
			{
				case 0:
					ClearMask=~Epsilon0;
					NIndex[node]=NIndex[node]&ClearMask;
					Cur=Val<<Epsilon0Sft;
					NIndex[node]+=Cur;
					break;
				case 1:
					ClearMask=~Epsilon1;
					NIndex[node]=NIndex[node]&ClearMask;
					Cur=Val<<Epsilon1Sft;
					NIndex[node]+=Cur;
					break;
				case 2:
					ClearMask=~Epsilon2;
					NIndex[node]=NIndex[node]&ClearMask;
					Cur=Val<<Epsilon2Sft;
					NIndex[node]+=Cur;
					break;
			}
		}
		inline int GetIonField(int i,int node)
		{
			return (NIndex[node]&IonField[i])>>IonFieldSft[i];
		}
		inline float GetDiffFloat(int i,int node)
		{
			return D[GetIonField(i,node)];
		}
		inline double GetDiffDouble(int i,int node)
		{
			return double(D[GetIonField(i,node)]);
		}
		inline float GetConcFloat(int i,int node)
		{
			return C[GetIonField(i,node)];
		}
		inline void SetIonField(int i,int node,int Val)
		{
			unsigned int ClearMask=~IonField[i];
			NodeIndex Cur=Val<<IonFieldSft[i];
			NIndex[node]=NIndex[node]&ClearMask;
			NIndex[node]+=Cur;
		}
		inline int GetDiffZeroInd()
		{
			int i;
			for(i=0;i<NodeIndexMaxValues;i++)
			{
				if(D[i]==0)return i;
			}
			return -1;
		}
		void SetDiffToZero(int ion,int node);
};

#ifdef SWIG
%nodefaultctor ContWorld;
#endif
/**
	This class combain all properties of the calculation box and basic operation on it (writing, reding arrays and also Most(if not all) of the MPI routines are here)
	@author Nikolay Simakov <nsimakov@andrew.cmu.edu>
*/
#ifdef SWIG
%feature("docstring") ContWorld "Class for continuum representation of system"
#endif
class ContWorld:public PnpsObject
{
	public:
		ContWorld();
		//!depreciated Default constructor for system with two types of ions
		ContWorld(int _GridSizeX,int _GridSizeY,int _GridSizeZ,float _GridScale, bool PBCX,bool PBCY,bool PBCZ,float q1,float q2);
		//!Default constructor
		ContWorld(int* m_GridSize,float m_GridScale, bool *m_PBC,int m_NIonsTypes,float *m_Qions);
#ifdef SWIG
		%pythoncode %{
		def __init__(self,
				GridSize=[65,65,65],
				GridScale=2.0,
				PBC=[False,False,False],
				Qions=None,
				Verbose=True):
			"""ContWorld constructor
			Input Parameters:
				GridSize=[int,int,int], default=[65,65,65]
					System grid size
				GridScale=float, default=2.0
					Grid scale in grids/Angstroms
				PBC=[Bool, Bool, Bool], default=[False,False,False]
					Set periodicity in X,Y,Z dimentions
				Qions=None|[float, float], default=None
					charges of mobile ions
					Examples:
						None for no module charges
						[1.0,-1.0] for 2 types of mobile ions with unit charges +1 and -1
				Verbose=bool, default=True
					Print world parameters
			Returned value:
				ContWorld - instance of ContWorld
			
			Example:
				contworld=pnpmod.ContWorld(
					GridSize=[65,65,65],
					GridScale=1.0,
					PBC=[True,True,False],
					Qions=[1.0,-1.0]
					)
			"""
			#validate input
			#GridSize
			if not (type(GridSize) is list or type(GridSize) is tuple):
				raise TypeError("ContWorld.__init__: GridSize must be list of three integers")
			if len(GridSize)!=3:
				raise TypeError("ContWorld.__init__:GridSize must be list of three integers")
			for i in range(3):
				if not isinstance( GridSize[i], int ):
					raise TypeError("ContWorld.__init__: GridSize must be list of three integers")

			#GridScale
			if isinstance(GridScale, int):
				GridScale=float(GridScale)
			if not isinstance(GridScale,  float):
				raise TypeError("ContWorld.__init__: GridScale must be float")
			if GridScale<=0.0:
				raise ValueError("ContWorld.__init__: GridScale must be positive")

			#PBC
			if not (type(PBC) is list or type(PBC) is tuple):
				raise TypeError("ContWorld.__init__: PBC must be list of three bools (True|False)")
			if len(PBC)!=3:
				raise TypeError("ContWorld.__init__: PBC must be list of three bools (True|False)")
			for i in range(3):
				if isinstance( PBC[i], int):
					PBC[i]=bool(PBC[i])
				if not isinstance( PBC[i], bool):
					raise TypeError("ContWorld.__init__: PBC must be list of three bools (True|False)")

			#Qions
			if Qions!=None:
				if not (type(Qions) is list or type(Qions) is tuple):
					raise TypeError("ContWorld.__init__: Qions must be list of floats")
				for i in range(len(Qions)):
					if isinstance( Qions[i], int):
						Qions[i]=float(Qions[i])
					if not isinstance(Qions[i], float):
						raise TypeError("ContWorld.__init__: Qions must be list of floats")
				Qtot=0.0
				for i in range(len(Qions)):
					Qtot+=Qions[i]
				if Qtot>0.0001 or Qtot < -0.0001:
					raise ValueError("ContWorld.__init__: sum of Qions must be zero")
			#Verbose
			if isinstance(Verbose, int):
				Verbose=bool(Verbose)
			if not isinstance( Verbose, ( bool ) ):
				raise TypeError("ContWorld.__init__: Verbose must be bools (True|False)")


			#Convert python to c++/c
			m_GridSize=new_intArray(3)
			for i in range(3):intArray_setitem(m_GridSize,i,GridSize[i])
			m_PBC=new_boolArray(3)
			
			for i in range(3):
				boolArray_setitem(m_PBC,i,PBC[i])
				#add padding for periodicity
				#rem if PBC[i]:
				#	intArray_setitem(m_GridSize,i,intArray_getitem(m_GridSize,i)+2)
			m_NIonsTypes=0
			m_Qions=None
			if Qions!=None and len(Qions)>0:
				m_NIonsTypes=len(Qions)
				m_Qions=new_floatArray(3)
				for i in range(m_NIonsTypes):floatArray_setitem(m_Qions,i,Qions[i])
			#swig auto generated stuff
			this = _pnpsll.new_ContWorld(m_GridSize, GridScale, m_PBC, m_NIonsTypes, m_Qions)
			try: self.this.append(this)
			except: self.this = this
			#delete temporary arrays
			if m_NIonsTypes!=0:delete_floatArray(m_Qions)
			delete_boolArray(m_PBC)
			delete_intArray(m_GridSize)

			if Verbose:
				self.Print()
		%}
#endif
		~ContWorld();
		
		int InitZero();
		int Clear();
	//protected:
	public:
		//MPI variables:
		int MyRank;	///< Rank of process
		int NProcs;	///< Number of processes
		int MyGlobalZ0,MyGlobalZ1;///<Number of first and last Z-coordinate of this process
		int **MPISpreadWorld;///< the way of spleating the system over processes when executed in parrallel, do i use it
	public:
		int NIonsTypes; ///< Number of Ion's Types in the system now only 2 is working
		float *IonsQ; ///< Charge of each ion type, now [1.0, -1.0] is tested
		
		int GridSize[3];///< Local grid size, including boarders for periodic boundary condition and interprocesses layer
		int GridSizeGlobal[3];///< Total grid size over all nodes, including boarders for periodic boundary condition
		int GridSizeOriginal[3];///< Original grid size, without boarders for periodic boundary condition
		
		int GS_XYZ,GS_XY,GS_X,GS_Y,GS_Z;///< Local grid size, including boarders for periodic boundary condition and interprocesses layer
		int GridSizeXYZGlobal;///< Total grid size over all nodes, including boarders for periodic boundary condition
		int GridSizeXYZOriginal;///< Original grid size, without boarders for periodic boundary condition
		
		int startBlackAndWhite;///<Is the first node of grid is black or white (chess board)
		
		float GridScale;///< Scale of the grid [grids/A], GridScale=1/h, where h is distense between two nodes of the grid
		
		NodeIndexing* NIndexing;///< Indexed map for dielectric constand, diffusion, initial concentration and static charge, it is alternative to *Epsilon[3],*StaticCharge,**D
		
		float *Epsilon[3];///< Dielectric Constant. Elements of vector show the displacement by h/2 in x, y and z direction. Ext. units:[1]. Int. units:[e/kT]
		float *Qstat;///< Static charges. Ext. units:[e/grid]. Int. units:[4*M_PI*GridScale*e]
		float **C;///< Ions' concentrations. Ext. units:[M]. Int. units:[COANGS[e/A^3]*4*M_PI*GridScale[1/A]*/(GridScale[1/A]*GridScale[1/A]*GridScale[1/A]]. Internal units is directly used in Poisson Solver
		double **CDouble;
		float **D;///<Diffusion. Units:[1E-6 cm^2/s].
		float **PMF;///< Potential of Mean Force [kT]
		float *Potential;///< Electrostatical potential [kT/e]
		
		double *PotentialDouble;
		
		bool PBC[3];///< Periodic Boundary Condition on x,y,z direction
		int BoundaryCondition;///< Boundary Condition, is separeted from world building one
		

		
		
		double SystemEnergy;///< Total Energy of the system
	protected:
		std::vector<std::string> BoundaryStr;
	public:
		//!Set Grid Size with correction for PBC
		int SetContWorld(int* m_GridSize,float m_GridScale, bool *m_PBC,int m_NIonsTypes,float *m_Qions);
		//!depreciated Set Grid Size without correction for PBC
		int SetContWorldNoMobIons(int _GridSizeX,int _GridSizeY,int _GridSizeZ,float _GridScale, bool PBCX,bool PBCY,bool PBCZ);
		//depreciated
		int SetContWorldTwoMobIons(int _GridSizeX,int _GridSizeY,int _GridSizeZ,float _GridScale, bool PBCX,bool PBCY,bool PBCZ,float q1,float q2);

		//!Print ContWorld Summary
		int Print();
		

		int ReadMaps(MapsIOData* Dt);
		int WriteMaps(MapsIOData* Dt);
		
		//!Check the grid of VF3D with this world and will return the same object if
		//!gid coinside or will interpolate and return new object,
		VectorField3D* CheckGridOfVectorField3D(VectorField3D* VF3D,bool deleteOldVF3DifNotSame,const char * VF3DDescription,const char * ReadFromFile);
		
		int SaveQstPhi(const char * filename);
		//!Read Potential from file
		int ReadPotential(const char * filename);
		int ReadPotentialChargeZRange(const char * filename,float Z0,float Z1);
		//!Write Potential to file, opt: 1 Original grid, 2 Local as is
		int WritePotential(const char * filename,int opt=1);
		//!Read PMF from file
		int ReadPMF(const char * filename,const char * filename2);
		//!Write PMF to file, opt: 1 Original grid, 2 Local as is
		int WritePMF(const char * filename,int opt=1);
		//!Read NodeIndexing from file
		int ReadNodeIndexing(const char * filename);
		//!Write NodeIndexing to file, opt: 1 Original grid, 2 Local as is
		int WriteNodeIndexing(const char * filename,int opt=1);
		int WriteContTop(const char * filename,int opt=1);
		//!Read PMF from file
		int ReadDynamicCharge(const char * filename);
		int ReadDynamicChargeZRange(const char * filename,float Z0,float Z1);
		//!Write PMF to file, opt: 1 Original grid, 2 Local as is
		int WriteDynamicCharge(const char * filename,int opt=1);
		int ReadDiffusion(const char * filename);
		int WriteDiffusion(const char * filename,int opt=1);
		VectorField3D* ReadVectorField3D(const char * filename,float coef);
		
		int BorderExchange(float *ArrayToExchange);
		int BorderExchangeInt(int *ArrayToExchange);
		int BorderExchangeDouble(double *ArrayToExchange);
		//!Check Arrays for presens(Epsilon, Diffusion) if BuildIfNI and NodeIndex is present will create arrays from it
		int CheckArrays(const char* ArSet,bool BuildIfNI);
		//!Check System for diffrent parameters
		int CheckSystem();
	public:
		int SetInitConcentrationFromNIndexing();
		
		int convertRealToOriginalGrid(float *mapIn, float *mapOut);
		int convertOriginalToRealGrid(float *mapIn, float *mapOut);
		int convertOriginalToRealGrid(int *mapIn, int *mapOut);
		inline float ConvFloatToGlobIntUnitsX(float Xext)
		{
			return Xext*GridScale+(float)int(GridSizeGlobal[0]/2);
		}
		inline float ConvFloatToLocIntUnitsX(float Xext)
		{
			return Xext*GridScale+(float)int(GridSizeGlobal[0]/2);
		}
		inline float ConvFloatToGlobIntUnitsY(float Yext)
		{
			return Yext*GridScale+(float)int(GridSizeGlobal[1]/2);
		}
		inline float ConvFloatToLocIntUnitsY(float Yext)
		{
			return Yext*GridScale+(float)int(GridSizeGlobal[1]/2);
		}
		inline float ConvFloatToGlobIntUnitsZ(float Zext)
		{
			return Zext*GridScale+(float)int(GridSizeGlobal[2]/2);
		}
		inline float ConvFloatToLocIntUnitsZ(float Zext)
		{
			return Zext*GridScale+(float)int(GridSizeGlobal[2]/2)-MyGlobalZ0;
		}
		inline int ConvGlobalIntenalXYZToGrdPnt(float Xint,float Yint,float Zint)
		{
			return int((Xint+0.5))+int((Yint+0.5))*GS_X+int((Zint+0.5))*GS_XY;
		}
		inline int ConvGlobalExternalXYZToGrdPnt(float Xext,float Yext,float Zext)
		{
			return ConvGlobalIntenalXYZToGrdPnt(ConvFloatToGlobIntUnitsX(Xext), ConvFloatToGlobIntUnitsY(Yext), ConvFloatToGlobIntUnitsZ(Zext));
		}
		inline float ConvGlobIntToGlobExtUnitsX(float Xint)
		{
			return (Xint-(float)int(GridSizeGlobal[0]/2))/GridScale;
		}
		inline float ConvGlobIntToGlobExtUnitsY(float Yint)
		{
			return (Yint-(float)int(GridSizeGlobal[1]/2))/GridScale;
		}
		inline float ConvGlobIntToGlobExtUnitsZ(float Zint)
		{
			return (Zint-(float)int(GridSizeGlobal[2]/2))/GridScale;
		}
		int GetBorder(int *z1,int *z2,int ProcRank);
	protected:
		int CalculateDerivativeParameters();///< run after setting of main parameters
	public:
		int AddPotential(float Z0,float Z1,float Phi0,float Phi1);
#ifdef SWIG
%rename (cxxApplyAsymConc) ApplyAsymConc;
		%pythoncode %{
		def ApplyAsymConc(self,z,C0,C1):
			"""ApplyAsymConc - Set diffrent concentration of mobile ions below and above coordinate z
				Parameters:
					z=float
						coordinate for concentration switching, usually middle of membrane is a good choice
					C0=float|[float,float] - single float number or list of floats
						concentration of mobile ions below z, if single float value is given will use that concentration for all ions
						otherwise will use specific concentration for each ion
					C1=float|[float,float] - single float number or list of floats
						concentration of mobile ions above z, if single float value is given will use that concentration for all ions
						otherwise will use specific concentration for each ion
			"""
			#convert python to c++/c
			Cz0 = new_floatArray(self.NIonsTypes)
			if type(C0) is list or type(C0) is tuple:
				for i in range(self.NIonsTypes):floatArray_setitem(Cz0,i,C0[i])
			else:
				for i in range(self.NIonsTypes):floatArray_setitem(Cz0,i,C0)
			Cz1 = new_floatArray(self.NIonsTypes)
			if type(C1) is list or type(C1) is tuple:
				for i in range(self.NIonsTypes):floatArray_setitem(Cz1,i,C1[i])
			else:
				for i in range(self.NIonsTypes):floatArray_setitem(Cz1,i,C1)
			
			self.cxxApplyAsymConc(z,Cz0,Cz1)

			#clean-up
			delete_floatArray(Cz1)
			delete_floatArray(Cz0)
		%}
#endif
		int ApplyAsymConc(float z,float* Cz0,float* Cz1);
		int ApplyAsymConcByScaling(float z,float* Cz0,float* Cz1);
	private:
		float AdditionalPotantialAlongZ(float Z1,float Z2,float Phi0,float Phi1,int InternalZ);
	public:
		int ConvertIonStrengthToDynamicCharge();
		int RemoveDiffusionPointsAtNegativeC();
		int RemoveDiffusionPoints(int IType,int* pnts,int n);
		//#ifdef SWIG
		//%apply float *OUTPUT { float *FindImplicitMembranePosition_z0,float *FindImplicitMembranePosition_z1 };
		//bool FindImplicitMembranePosition(float *FindImplicitMembranePosition_z0,float *FindImplicitMembranePosition_z1);
		//#else
		bool FindImplicitMembranePosition(double *z0,double *z1);
		//#endif
		
#		ifdef HARLEM_MOD
		HaVec_float* GetNumberOfIonsAlongChannel(int ion, float Xcenter, float Ycenter, HaVec_float* RlimAlZ);
		double GetNumberOfIonsInChannel(int ion, float Xcenter, float Ycenter, float Z0, float Z1, HaVec_float* RlimAlZ);
		HaVec_float* GetRadiusAlongZ(int ion, float Xcenter, float Ycenter, HaVec_float* RlimAlZ);
		HaVec_float* GetDavrAlongZ(int ion, float Xcenter, float Ycenter, HaVec_float* RlimAlZ);
		HaVec_float* GetPotValusAlongZ(int ix,int iy);
		HaVec_float* GetConcValusAlongZ(int ion,int ix,int iy);
		HaVec_float* GetPMFValusAlongZ(int ion,int ix,int iy);
		int CorrectDiffusion(int ion, float Xcenter, float Ycenter, float Z0, float Z1, HaVec_float* RlimAlZ,HaVec_float* DscaleAlZ);
		int SetPotentialToZeroWhereDiffZero(int ion);
		HaVec_float* GetAvrPotAlongPore(int ion, float Xcenter, float Ycenter, HaVec_float* RlimAlZ);
#		endif

	#ifdef SWIG
	%pythoncode %{
		def AddPotentialAuto(self, PotDiff):
			z0=new_doublepf()
			z1=new_doublepf()
			if not self.FindImplicitMembranePosition(z0,z1):
				return None
			LimitCurrentCalcZ0=doublepf_value(z0)
			LimitCurrentCalcZ1=doublepf_value(z1)
			self.AddPotential(LimitCurrentCalcZ0,LimitCurrentCalcZ1,PotDiff,0.0)
			return [LimitCurrentCalcZ0,LimitCurrentCalcZ1]
	%}
	#endif
};
#ifdef SWIG
//%clearnodefaultctor;// ContWorld;
#endif

ContWorld* ReadContWorldFromNodeIndexing(const char * node_indexing_filename);

#ifdef SWIG
%pythoncode %{
def GetContWorld(**kwargs):
	"""depriciated"""
	GridSize = kwargs.get("GridSize")
	if GridSize==None:
		print("Error: GridSize must be defined")
		return None
	GridScale = kwargs.get("GridScale")
	if GridScale==None:
		print("Error: GridScale must be defined")
		return None
	
	PBC = kwargs.get("PBC",[0,0,0])
	Qions = kwargs.get("Qions",[1.0,-1.0])
	
	w=ContWorld(GridSize[0],GridSize[1],GridSize[2],
		GridScale,
		PBC[0],PBC[1],PBC[2],
		Qions[0],Qions[1])
	w.thisown = 0
	return w

def LoadContWorld(**kwargs):
	"""
	Options:
	SysTop - filename of continious topology
	Potential     - potential
	Concentration    - concentration
	Diffusion - Diffusion
	PMF     - filename for PMF of SR
	"""
	SysTop=kwargs.get("SysTop",None)
	

	
	if SysTop == None:
		print("Error: ContTop must be specified")
		return None
	w=ReadContWorldFromNodeIndexing(str(SysTop))

	Potential=kwargs.get("Potential",None)
	if Potential != None:
		w.ReadPotential(str(Potential))

	Concentration=kwargs.get("Concentration",None)
	if Concentration != None:
		w.ReadDynamicCharge(str(Concentration))

	Diffusion=kwargs.get("Diffusion",None)
	if Diffusion != None:
		w.ReadDiffusion(str(Diffusion))

	PMF=kwargs.get("PMF",None)
	if PMF != None:
		w.ReadPMF(str(PMF),"")

	return w
%}
#endif
#endif
