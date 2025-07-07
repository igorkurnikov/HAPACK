//
// C++ Interface: elmod
//
// Description: 
//
//
// Author: mikola <mikola@linux>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ELMOD_H
#define ELMOD_H

#include "electrostmod.h"
#include "haatgroup.h"
#include <wx/string.h>
#include <vector>
class ContWorld;
class TiXmlElement;
class HaAtom;
class ProtonRedoxMod;

class HaAtomParmEntry
{
	public:
		HaAtomParmEntry(const char* cName,double _dbl0,double _dbl1);
		~HaAtomParmEntry();
		
		void SetValues(double _dbl0,double _dbl1);
		
		double dbl0;
		double dbl1;
		
		std::string Name;
};
class HaResParmEntry
{
	public:
		HaResParmEntry(const char* cName, const char* cMod);
		~HaResParmEntry();
		
		//!if atom name already in the list then the values will be rewrited, otherwize new entry will append to the list
		HaAtomParmEntry* SetAtom(const char* cName,double _dbl0,double _dbl1);
		HaAtomParmEntry* GetAtom(const char* cName);
		int GetAtomListNum(const char* cName);
		void DelAtom(const char* cName);
		
		void CopyAtomsFromRes(HaResParmEntry* CpRes);
		
		void PrintInfo();
		void PrintSumDbl0();
		
		std::vector<HaAtomParmEntry*> Atms;
		
		std::string Name;
		std::string Mod;
};
class HaAtomsParmFF
{
	public:
		HaAtomsParmFF(const char* cName, const char* cRef, const char* cNotes);
		~HaAtomsParmFF();
		
		HaResParmEntry* NewRes(const char* cName, const char* cMod);
		//!if resedue present then return current else will allocate new
		HaResParmEntry* SetRes(const char* cName, const char* cMod);
		HaResParmEntry* GetRes(const char* cName, const char* cMod);
		HaResParmEntry* CopyRes(const char* cName, const char* cMod,const char* cNameFrom, const char* cModFrom);
		
		void PrintFFInfo();
		void PrintSumDbl0();
		void PrintEntries();
		
		std::vector<HaResParmEntry*> Ress;
		
		std::string Name;
		std::string Ref;
		std::string Notes;
};
class HaAtomsParmDB
{
	public:
		HaAtomsParmDB(const char* cName);
		~HaAtomsParmDB();
		
		HaAtomsParmFF* NewFF(const char* cName, const char* cRef, const char* cNotes);
		HaAtomsParmFF* GetFF(const char* cName);
		int NumFF();
		HaAtomsParmFF* GetFFbyNum(int i);

		void SetAtomsParam(HaMolSet* pmset, const char *FFName);
		
		//! compare atom names for residues in cName with that in cNameRef
		void CompareAtomNamesInFFs(const char* cName,const char* cNameRef);
		//! print params for res in FF cName and the same in cNameRef
		void PrintSideToSideFFs(const char* cName,const char* cNameRef);
		
		void PrintFFsInfo();
		
		std::vector<HaAtomsParmFF*> FFs;
		
		std::string Name;
};
HaAtomsParmDB* GetQRDB();
void DeleteQRDB();
/**
Electrostatic Module, replacement of ElectrostMod, mimic its behaviour

@author mikola
*/
class ElMod : public ElectrostMod
{
	public:
		ElMod(HaMolSet* new_phost_mset = NULL);
		
		~ElMod();
		
		friend class pKaCalcMod;
		
		virtual void clear();
		virtual bool run(RunMode rmode=RUN_FOREGROUND);
		#ifdef WITH_CUDA
		virtual bool runCPS(RunMode rmode=RUN_FOREGROUND);//!<Run CUDA Accelerated Poisson Solver
		#endif
		
		//Extra parameters
		int BoxSetMethod;
		enum {SetBoxPerFil=0,SetBoxGridScale=1};
		double GridScale;
		
		bool bCalMolSet;
		int MoleculePosition;
		enum 
		{
			MoveToCenter=0,
			DontMove=1,
			UseOffset=2,
			TranslateAndRotate=3
		};
		double OffsetX,OffsetY,OffsetZ;
		
		int MembraneType;
		double epsmem;
		double memZ0,memZ1;
		double memR0,memR1;
		double memX,memY;
		enum 
		{
			MemT_None=0,
			MemT_MembraneZ=1,
			MemT_Tube=2
		};
		int BldType;
		enum
		{
			BldNI=0,
			BldEu=1,
			BldCmp=2
		};
		void SetMemTube(double m_x, double m_y, double m_z0, double m_z1, double m_R0, double m_R1, double m_eps);
		void SetMemNone();
		float Relaxation;
		float Convergence;
		
		bool RemovingCavitiesOnDielectricMap;
		int RemCavOnDielFillWith;
		int RemCavOnDielWhere2Look;
		
		bool bInsertProbe;
		double probe_eps;
		double probe_x,probe_y,probe_z,probe_q,probe_R;
		double Phi,Theta;
		void InsertProbe(double m_x,double m_y,double m_z,double m_q,double m_R,double m_eps);
		void UnsetProbe();
		
		void SetTranslateAndRotate(double m_x,double m_y,double m_z,double m_Phi,double m_Theta);
	public:
		ContWorld* m_ContWorld;
		ElMod* GetRoughWorld();
	public:
		double gr_cent[3];
		int InitContWorld();
		int BuildWorld();
		int SolvePoisson();
		int SolvePoissonBoltzmann();
		int SolvePoissonBoltzmannSoftRepultion();
		//BuildWorldNI* m_BuildWorldNI;
		//PoissonSolver* m_PoissonSolver;
		//PoissonBoltzmannSolver* m_PoissonBoltzmannSolver;
		ElMod* RoughWorld;
	public:
		int SetNinaBelowRouxRadii(bool only_selected_atoms=false);
	public:
		bool bPBSR;
		bool bConvIndConc;
		HaVec_float IonsHalfSigma;
		HaVec_float IonsFourEpsilon;
		std::string PQRFile;//!<will load charges and radii from this file
		std::string PREFile;//!<will load params for soft repulsion from this file
		
		int LoadConcMaps(const char *filename);
		
		//!Load coordinates, charges and radii to host molset
		int LoadPQR(const char *filename, bool LoadXYZ=true,bool LoadQZ=true,bool LoadRZ=true);
		
		//!Save coordinates to temporary arrays, just check atom name and resname
		int SaveCoorToTmpArr();
		int LoadCoorFromTmpArr();
		int DeleteCoorTmpArr();
		
		int RotateMolSet(Vec3D& n,double a);
		void TranslateAndRotateMolSet(double m_x,double m_y,double m_z,double m_Phi,double m_Theta);
		
		//!calculate MinMax including R
		int GetMinMaxWithR();
		//!Calculate biggest radius which enclose the system (starts at origin)
		int CalcRAroundZero();
		
		double *t_x,*t_y,*t_z;
		
		bool GPUAccelerated;
		int ChargeDist;
		int ChargeDistPar;
		int MemoryLimitOnOneArray;
};

class pKaCalcMod : public  HaCompMod
{
	public:
		pKaCalcMod(HaMolSet* new_phost_mset = NULL);
		virtual ~pKaCalcMod();
		
		int PrintResWithAltProtState();
		int PrintResults();
		int PrintPopulation();
		int PrintPop4Homooligamer(bool PrintOnlyIfAltPopMoreThenSmth=false,float Smth=0.1);
		
		int CalcpKaUsingElectrostMod();
		int CalcpKaUsingElMod();
		
		int RunCalcUsingElectrostMod();
		int RunCalcUsingElMod();
		
		int MakeAltStList();
		//!ReadCalculatedEnergies from saving, good if run was aborted
		int ReadCalculatedEnergies(const char *filename);
		int CalcIntrpKa();
		int CalcpKaWithInteraction();
		
		int WritepKaCalcModToFile(const char *filename);
		int ReadpKaCalcModFromFile(const char *filename);
		int WritepKaCalcModToXmlElement(TiXmlElement *RootElt);
		int ReadpKaCalcModFromXmlElement(const TiXmlElement *RootElt);
		
		int SetAltSt4ResInHomoolgmr(int ResNum);
		
		int NumberOfAltStates;
		double E1;
		HaVec_double E2,E3,E4,ddG;
		HaVec_double pKa, pKaFromPop,pKaIntr,pKaStd,dpKa;
		HaMat_double inter_mat;
		
		//in case if run was aborted will reload values which was already done
		int E1done;
		HaVec_int E2done,E3done,E4done;
		HaVec_int inter_mat_done;
		
		HaMat_double Pop;
		
		double pHmin;
		double pHmax;
		double pHstep;
		int pKaCalcMethod;
		enum {SCF_MULTI_SITE_CALC=0,MC_MULTI_SITE_CALC=1};
		std::vector<std::string> pKaCalcMethodStr;
		int MC_PKA_CALC_N_mc_cyc;
		int SCF_PKA_CALC_max_iter;
		double SCF_PKA_CALC_pop_err_max;
		
		bool SaveIntermediatePotNNI;
		bool SaveIntermediateResults;
		
	protected:
		ProtonRedoxMod* p_prot_rdx_mod;

		AtomGroup sel_atoms;
		VecPtr act_chem_states;
		vector<string> AltNames;
		
		//! interaction matrix - on the diagonal energy difference of alternative and unmodifed state when 
		//! alt_pop[j] - population of aternative states  
		//! delt_e[j]  - converged self-consistent difference between alternative and unmodifed state 
		int RoutineCalcAvgPopMC(HaMat_double& inter_mat, HaVec_double& alt_pop, HaVec_double& delt_e,int N_mc_cyc);

		//! interaction matrix - on the diagonal energy difference of alternative and unmodifed state when 
		//! alt_pop[j] - population of aternative states  
		//! delt_e[j]  - converged self-consistent difference between alternative and unmodifed state 
		int RoutineCalcAvgPopSCF(HaMat_double& inter_mat, HaVec_double& alt_pop, HaVec_double& delt_e,int max_iter = 100,double pop_err_max = 0.0001);
};
class PNPMod : public  HaCompMod
{
	public:
		PNPMod(HaMolSet* new_phost_mset = NULL);
		~PNPMod();
	protected:
		int NIonsTypes;
		std::vector<std::string> IonNames;
		
		HaMat_double *LJA;
		HaMat_double *LJB;
		std::vector< HaAtom* > AtmsSet;
		
	public:
		int SetNIonsTypes(int m_NIonsTypes);
		int SetIonName(int ion,const char * name);
		const char * GetIonName(int ion);
		int PrintLJAB();
		int SavePABFile(const char* filename);
	public:
		//!Read Lennard-Jones parameters from AMBER94 FF
		int ReadAMBERFF(const char* filename);
		int ReadAMBER94FF();
		int GetAtomTypeNumber(std::string atomname);
		double GetHalfSigma(std::string atomname);
		double GetFourEpsilon(std::string atomname);
		//!Save pdb with Lennard-Jones EpsilonStr[kT] and RStar[A] from AMBER FF
		int SavePREFile(const char* filename);
		int SavePREFreeFile(const char* filename);
		
		int SetLJABfromAMBERFF();
	public:
		std::vector<std::string> AtomTypesDB;
		std::vector<double> HalfSigmaDB;//!<Half of Sigma DB[A]
		std::vector<double> FourEpsilonDB;//!<Four Epsilon DB[kT]
		//std::vector<std::string> AtomTypesDB;
		
	public:
		int SetLJABfromOPLS();
		int ReadOPLSFF();
		int PrintOPLSLJSigmaEpsilon();
		int PrintOPLSLJAB();//!<A = 4*epsilon*sigma^12 [kT A^12]; B = 4*epsilon*sigma^6 [kT A^6] 
	protected:
		int ReadOPLSitp(const char *filename);//!<read non-bonded params
		int ReadOPLSrtp(const char *filename);//!<read atoms in residues
		
		int GetOPLSEpsilonSigma(const char *AtmNM, const char *ResNM, double *Eps, double *Sgm);
		
		std::vector<wxString> OPLSLJAtomTypes;
		std::vector<wxString> OPLSLJAtomName;
		std::vector<double> OPLSLJSigma;//!<Sigma DB[A]
		std::vector<double> OPLSLJEpsilon;//!<Epsilon DB[kT]
		//res-atoms name
		std::vector<wxString> OPLSResNames;
		std::vector< std::vector<wxString> > OPLSResAtomName;
		std::vector< std::vector<wxString> > OPLSResAtomTypes;
		std::vector< std::vector<double> > OPLSResAtomCharge;
	
	public:
		int SaveIER(const char* filename, bool OnlyHeavyAtoms=true);
		int ReadIER(const char *filename,bool AddToDB);//!<read ion exclusion radii
		
		int GetIER(HaAtom  *aptr, double *rK, double *rCl, bool OnlyHeavyAtoms);
		int GetResNumAtIERDB(wxString* ResName);
		int GetAtmNumOfResAtIERDB(int myres,wxString* AtmName);
		
		std::vector<wxString> IERResNames;
		std::vector< std::vector<wxString> > IERAtomName;
		std::vector< std::vector<double> > IERRadiusK;
		std::vector< std::vector<double> > IERRadiusCl;
	public:
		int SavePAN(const char* filename, bool OnlyHeavyAtoms=true);
		int AssignPAN(bool OnlyHeavyAtoms=true);
		int ReadPANDB(const char *filename,bool AddToDB);//!<read ion exclusion radii
		
		int GetSR_AN(HaAtom  *aptr, double *AK, double *NK, double *ACl, double *NCl, bool OnlyHeavyAtoms);
		int GetResNumAtSR_AN_DB(wxString* ResName);
		int GetAtmNumOfResAtSR_AN_DB(int myres,wxString* AtmName);
		
		std::vector<wxString> SR_AN_ResNames;
		std::vector< std::vector<wxString> > SR_AN_AtomName;
		std::vector< std::vector<double> > SR_A_K;
		std::vector< std::vector<double> > SR_A_Cl;
		std::vector< std::vector<double> > SR_N_K;
		std::vector< std::vector<double> > SR_N_Cl;

		std::vector< float > mSR_A_K;
		std::vector< float > mSR_A_Cl;
		std::vector< float > mSR_N_K;
		std::vector< float > mSR_N_Cl;
	public:
		ContWorld* m_ContWorld;
		int RunPNPSFromString(const char* string);
};
class ElModRadDist
{
	public:
		ElModRadDist(HaMolSet* new_phost_mset);
		~ElModRadDist();
		HaMolSet *pmset;
		ElMod *elmod;
		int atom;
		
		int SetMinMax(float min,float max,float binsize);
		int Setheff(float _heff);
		int PrintRadDist(const char *filename);
  public:
		float Rmin;
		float Rmax;
		float BinSize;
		float heff;//!<effective step size for integration
		int nbins;
		//HaVec_double Rleft;
	public:
		HaVec_double** gcont[2];
    HaVec_double* Getgcont(int ion,int atm);
		HaVec_int AtomsList;
		
		int CalcRadDist();
		
		int CalcRadDistByGrid();
    
    HaVec_double** wel;
    int CalcWelRadDistByGrid();
};
#endif
