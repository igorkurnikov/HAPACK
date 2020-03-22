/*! \file haqchem.h

     Basic Classes to perform Quantum Chemical Calculations in HARLEM.

     \author Igor Kurnikov 
	 \date 1998-2002
*/

#ifndef HAQCHEM_H
#define HAQCHEM_H

class HaAtBasDB;
class HaPseudoPotDB;
class HaField3D;
class HaMatDB;

class InternalBasis;

#include "haatombasis.h"
#include "halocorb.h"
#include "haatgroup.h"
#include "command.h"
#include "hacompmod.h"
#include "qc_params.h"

typedef std::vector<GauAtomBasis> AtBasisType;
class GauFile;
class ZMatCrd;

namespace harlem
{
class RunOptions;
};

namespace harlem
{
	namespace qc
	{
		enum NDO_METHOD {CNDO_2 = 1, INDO_2, ZINDO_1, ZINDO_S, HUCKEL};
		enum SCRF_METHOD { SCRF_NO = 0, SCRF_PCM };
		enum WAVE_FUN_TYPE {HARTREE_FOCK, NDO, EXTENDED_HUCKEL, MP2, DFT_B3LYP, CCSD_T };
	};
};

//! \brief Computational module to perform Quantum Chemical calculations
//! \nosubgrouping
class HaQCMod: public HaCompMod
{
	friend class HaGaussMod;
	friend class QChemParDlgWX;
	friend class HaTests;

public:

	HaQCMod(MolSet* new_phost_mset);
	virtual ~HaQCMod();

	void SetStdParams();

	bool Print_info( ostream &sout, const int level) const;
	virtual int SaveXMLToStream(std::ostream& os, const harlem::SaveOptions* popt = NULL ) const;

//! \name Atom Center functions:  
//@{
	int  GetNumCnt() const;                           //!< Get Number of Centers (atoms and point charges)
	bool GetCntCharges(HaVec_double& charges) const;  //! Get a Vector of Center Charges
//@}

//! \name Molecule Electronic Wave function manipulation functions: 
//@{
	void SetCharge(int charge); //!< Set Total charge of the molecule 
    void SetMult(int mult);     //!< Set Multiplicity of the molecule 
    int GetCharge() const;        //!< return charge of the molecule
    int GetMult()   const;        //!< return the molecule multiplicity 

	void SetField(double el_x, double el_y, double el_z); //!< Set Electric field in the system in MM units electron/(Ang*Ang)
	void SetFieldAU(double el_x, double el_y, double el_z); //!< Set Electric field in the system in AU units (electron/(Ang*Ang))

    int GetNelectr() const;       //!< return the number of electrons in the molecule
	int GetNumAlphaEl(int active_bas=0) const; //!< Get number of Alpha Active Electrons in AtBasis or ActBas
	int GetNumBetaEl(int active_bas=0) const;  //!< Get number of Beta Active Electrons  in AtBasis or ActBas

	inline int GetNumOccMO() const { return ( this->GetNumAlphaEl()); } //!< Get the number of occupied MOs
	inline int GetNumVacMO() const { return ( this->GetNumMO() - this->GetNumAlphaEl()); } //!< Get the number of unoccupied MOs
	int GetNumMO()    const; //!< Get the number of MOs

	bool SetWaveFunType(const char* str_wf_type); //!< Set Wave Function type
//@}

//! \name Basis set manipulation functions:
//@{
    bool InitBasis(GauFile& gfile);           //!< load Basis Set from Gaussian rwf file ( not completed yet)
	bool InitBasis(const std::string& bname); //!< Initiate Basis set with a given name on all atoms 

    int FBasFunPos(const HaAtom* ref_aptr) const;      //!< get the position of the first basis function of the Atom 
 	const HaAtom* GetAtomOfAO(const int idx_AO) const; //!< Get the atom of the basis function

//	GauBasisSet* GetCurBasisSet() { return &AtBasis; } //!< Return pointer to a current Basis Set 
    int GetNBfunc() const;        //!< Return the number of basis functions 

	std::string GetBasName() const { return m_bas_name; }  //!< get the basis set name
	bool UsePseudoPot() const;    //!< Check if Pseudo Potential is Used
//@}

//! \name Local Active Orbital manipulation functions
//@{
	bool InitLocOrb(const char* setid ); //!< Initiate Local(Active) Orbitals with a method given by setid idenfificator  

	int GetNActiveOrb() const;           //!< Get the number of active orbitals

	bool IsLocOrbFullBasis();   //!< Check if Local Orbital basis coincide with Full Basis set

	int GetLocOrbIdxOfGrp(const std::string& gid , HaVec_int & ilgr ) const; //!< get the list of indexes of local orbitals of the group

	bool ExtractLocOrbSubMat(const std::string & gid1, const std::string & gid2,
		                 const HaMat_double & ActOrbMat, 
						 HaMat_double & ActOrbSubMat) const; 

	bool InsertLocOrbSubMat(const std::string & gid1, const std::string & gid2,
		                  HaMat_double & ActOrbMat, 
						  const HaMat_double & ActOrbSubMat) const; 

	const char* GetLocOrbSetID() const { return m_loc_orb_set_id.c_str(); }
//@}

//! \name Display related functions and data:
//@{
	static bool EvalLinCombOnGrid(const HaVec_double& orb_coef, ArrayOrb3D& bas_set, HaField3D& mo_grid); //!< Build a grid of MO values 
	VecPtr CreateOrbContour(const HaVec_double& orb_coef, ArrayOrb3D& bas_set, const double mo_isolvl=0.1, const int ngrid= 11);
	bool CreateMOcontour (const int imo, const double mo_isolvl=0.1, const int ngrid= 11);	
	
	int m_grid_size;
//@}

//! \name Wave function manipulation functions:
//@{
	bool BuildFockMatFromMOs( HaMat_double& fock_matrix, double cut_ene = -100000.0); //!< Build Fockian matrix from Molecular orbitals and their energies

	static QCIntEngineType GetQCIntEngine();       //!< Get Current Quantum Chemical Intgral Engine Type
	static void SetQCIntEngine(const QCIntEngineType& int_engine_new ); //!< Set Quantum Chemical Intgral Engine Type
	
	static QCIntEngineType int_engine; //!< integral engine IPACK, GAUSSIAN if available or others 	
//@}

//! \name Gaussian specific functions:
//@{
	bool InitBasOvlp();                   //!< Load AO basis Overlap Matrix
    bool Init1eDens(GauFile& gfile);      //!< Load 1e Density from Gaussian rwf file
	bool InitMOs(GauFile& gfile);         //!< Load MO coeffcients and energies
	static int LoadGauCom(const GauBasisSet& gbas);     //!< Load necessary Gaussian commons for Gaussian basis set gbas

	bool LoadDataFromFChk(const char* fname);  //!< Load data from Gaussian Formatted Checkpoint File

	bool load_mo_flag; //!< flag to load molecular orbitals from Gaussian Formatted Checkpoint File

	int ProjMatToActBas(HaMat_double & fmat, HaMat_double & fmat_lb); //!< transform matrix from the full basis to the active basis set
	int CalcEPfromMO(HaMat_double& gm,double ene); //!< Compute electron propagator matrix from MO

	HaMat_double& GetOvlpMat();  //!< Get Reference to the overlap matrix of the main basis recomputing if needed
//@}

//! \name Hartree Fock Calculations parameters
//@{
	int max_scf_iter;    //!< Maximum number of SCF iterations
	int max_it_avg;      //!< Number of SCF iteration with averaging of density matrix
	int max_it_noavg;    //!< Number of SCF iteration without averaging of density matrix
    int iuhf;            //!< UHF flag  
	double conv_dm;      //!< SCF convergence on density matrix
	double temp0_fermi;  //!< initial temperature to form density matrix from molecular orbitals 
	int iter_temp;       //!< number of step used to switch off effective temperature for density matrix formation

    int guess_only;        //!< compute only guess 
	int set_guess_from_mos; //!< Set init SCF guess from current MOs
//@}

#ifdef USE_IPACK
	int TestIPack1(); 
	int TestIPack2();
	int TestRandomGen();
#endif
	
	static int InitIPack(); //!< initialize timer and other objects needed for IPACK


//! \name GAUSSIAN LIB PARAMS
//@{
	static int max_gauss_mem; //!< Maximum memory to allocate for Gaussian functions	
	static void set_max_gauss_mem(int new_max_mem); //!< Set Maximum memory to allocate for Gaussian functions
//@}

//! \name  Job Control Parameters 
//@{
public:

	int Run(const harlem::RunOptions* popt = NULL ); //!< Run calculations 
    int StopCalc(); //!< Stop running caclulations

	int PrepGauss();   //!< Prepare execution of Gaussian subroutines
    int RunCNDO(); //!< Run CNDO/2 INDO/2 or ZINDO/S calculations as a thread if possible 
	int RunCNDOThread(); //!< Run CNDO/2 or INDO/2 or ZINDO/S calculations 
	int RunExtHuckel();  //!< Run Extended Huckel Calculations

    int stop_calc_flag; //!< flag to stop running calculations

	void SetSinglePtCalc();   //!< Set Single Point QChem calculations
	void SetEneMinCalc();     //!< Set Energy minimization calculations
	void SetTransStateCalc(); //!< Set Transition State Energy calculations

	bool IsSinglePtCalc() const; //!< Check if Single Point Calculations
	bool IsEneMinCalc() const;   //!< Check if Energy minimization Calculations
	bool IsTransStateCalc() const; //!< Check if Transition State Energy calculations

	int CalcEnergy(); //!< Calculate Energy (Perform Single Point Calculation) for the system using current QChem method
	int RunMinEne();  //!< Run Energy Minimization Calculations 

	bool IsGenBasisSet() const;  //!< Is Basis Set considered to be generic in calculations
	void SetBasisSetGen(); //!<  Set Basis set as generic (not-standard) in calculations

	bool IsUsingSymmetry() const; //!< use symmetry of the molecule in calculations if possible
	void SetUseSymmetry( int use_symmetry_par ); //!< Set to use( not to use) molcular symmetery in calculations if possible 

	void SetHF();    //!< Set Hartree-Fock calculations
	void SetMP2();   //!< Set MP2 calculations
	void SetDFT();   //!< Set DFT calculations with default hamiltonian (B3LYP)
	void SetB3LYP(); //!< Set DFT calculations with B3LYP functional
	void SetCCSD_T(); //!< Set Couple Cluster CCSD(T) calculations

	bool IsHF() const;    //!< Check if HARTREE-FOCK Calculations
	bool IsDFT() const;   //!< Check if DFT Calculations
	bool IsB3LYP() const; //!< Check if DFT calculations with B3LYP hamiltonian
	bool IsMP2() const;   //!< Check if MP2 Calculations
	bool IsCCSD_T() const; //!< Check if CCSD(T) calculations

	void SetCalcPolar( bool set_par = true ) { calc_polar = set_par; }
	bool ToCalcPolar() const { return calc_polar; } 

	bool IsUsingSCRF() const; //!< Use SCRF solvent model
	void SetSCRF(); //!< Set Default SCRF (self-consistent reaction field model of the solvent ) method
	void SetSCRFMethod( int set_method_par ); //!< Set SCRF (self-consistent reaction field model of the solvent ) method 

	harlem::qc::SCRF_METHOD scrf_method; //!< SCRF (self-consistent reaction field model of the solvent ) method  

	int SetExtCharge( HaAtom* aptr, double ch);       //!< Set external charge at position of atom aptr
	int SetExtCharge( const std::string& at_ref, double ch); //!< Set external charge at position of atom with reference at_ref
	void SetExtChCrdOffset( double crd_offset ); //! Set shift in coordinates of external charges relative atom positions to avoid infinite energies (default 0.01)

	harlem::qc::WAVE_FUN_TYPE wave_fun_type; //!< Type of electronic function calculations
	int ndo_method; //!< method for NDO calculations = 1 CNDO/2, =2 INDO_2, =3 ZINDO_1, =4 ZINDO_S, = 5 - HUCKEL

	ZMatCrd* GetZMat(); //!< Get Z-matrix associated with the module

	void SetPrefix( const std::string& prefix ); //!< Set prefix for files created during calculations
	std::string GetPrefix() const; //!< //!< Get prefix for files created during calculations

protected:

	std::string prefix;

	int FormDenMat(const HaMat_double& cmo, double* pa, int nel, double* pene_mo = NULL, double temp_fermi = 0.0); //!< Form density matrix in pa
	static double ZIndoGInt(double r, double k1, double k2);
    static double RRIntSTO(int it, int n1, double a1, int n2, double a2, int n3, 
				           double a3, int n4, double a4); //!< Calculate radial Coulomb integral for STOs
	int NDOExp(int ia , int l, int ovlp, int& n, HaVec_double& cf, HaVec_double& exp);
	int CNDOInteg(int iprint, int iatom,int natoms,
				   const HaVec_int& ian, const HaMat_double& c,
				   double* ss, double* ss_scaled, double* gss, double* gsd, double* gdd,
				   const HaVec_int& ilst_bf_at, const HaVec_int& ifst_bf_at);

    int SetNDOCore(int natoms, int nbzdo,
	               double* fmat,double* gss, double* gsd, double* gdd, int* ian,
                   int* ifst_bf_at, int* ilst_bf_at,
				   double* core_ch, int* iat_bf, double* fa); //!< Set NDO core hamiltonian
	
	int SetINDOAtCoulMat(int natoms, HaVec_double& cm, int ia, int elem, double* gss, double* gsd, double* gdd); //!< set matrix of INDO coulomb matrix elements 

public:

	int InitHuckParsStd();  //!< Initiate standard Extended Huckel parameters
	int InitHuckParsVela(); //!< Initiate Extended Huckel parameters
	int InitHuckHam(HaMat_double& hmat, HaMat_double& ss, ArrayOrb3D& bas); //!< Init Huckel Hamiltonian 
	double GetNDOValEl(int elem, double& ns_val, double& np_val, double& nd_val,double& nf_val); //!< Get the number of valence electron and electron dencities of S,P,D valence shells
    int FormFockNDO(int natoms, int* ian, const int* ifst_bf_at, const int* ilst_bf_at,
					double* da, double* db, double* fm, 
					const double* gss, const double* gsd,const double* gdd,
					vector<HaVec_double>& at_coul_int);
	double CalcNucRepEne(int natoms, const HaVec_int& ian, const HaMat_double& c, 
		                 const  HaVec_double& core_ch, const double* gss);
	
	GauBasisSet AtBasis;          //!< Current Gaussian Basis Set set in the module
	ArrayOrb3D* ActBas;           //!< Active Basis Set 

	int allocated_act_basis;   //!< flag to indicate that active basis is internally allocated and has to be deleted when the module is deleted

	std::string m_bas_name;
	std::string m_loc_orb_set_id;

	int charge;           //!< total charge
	int mult;             //!< multiplicity

	Vec3D GetDipole() const       { return dipole; } 
	Vec3D GetTotDipole() const    { return dipole; } 
	HaVec_double GetQpole() const { return qpole; }
	HaVec_double GetTotQpole() const { return qpole; } 
	HaVec_double GetPolarTensor() const { return polar; } 

	static HaMatDB* p_ndo_pars_db;

	void SetEne( double ene);         //!< Set Total Energy of the system
	void SetHFEne( double ene_hf );   //!< Set Hartree-Fock energy of the system
	void SetDFTEne( double ene_dft ); //!< Set DFT energy of the system
	
	double GetEne()    const; //!< Get total energy of the system 
	double GetHFEne()  const; //!< Get Hartree-Fock energy of the system
	double GetDFTEne() const; //!< Get DFT energy of the system
	
	HaMat_double ovlp_mat;
	HaMat_double huck_ham;
	HaMat_double MO_coef; 
	HaVec_double MOene;

protected:

	double ene;     //!< Total energy of the system 
	double ene_hf;  //!< Hartree-Fock Energy 
	double ene_dft; //!< DFT energy

	Vec3D el_field; //!< External Dipole Electric Field ( electron/(Ang*Ang)) 
	double distr_ext_chrg;     //!< Solvent charge distributed over points of dot surface 

	AtomDoubleMap ext_chrg; //!< map of extra external charges associated to atoms 
	double ext_ch_crd_offset; //!< shift in coordinates of external charges relative atom positions to avoid infinite energies ( default 0.01 )

	Vec3D dipole;        //!< total dipole     (electron * Ang )
	HaVec_double qpole;  //!< total quadrupole (electron * Ang * Ang )
	HaVec_double polar;  //!< polarizability tensor ( Ang**3) ( alpha_xx, alpha_yy, alpha_zz, alpha_xy, alpha_xz, alpha_yz )
	
	HaVec_double ene_history;
	
	enum {SINGLE_PT_CALC, ENERGY_MIN_CALC, TRANS_STATE_CALC } calc_type; //!< Type of quantum chemical calculations
	int use_symmetry; //!< parameter use symmetry in calculations  

	static harlem::RunOptions run_opt_default;

	bool calc_polar; //!< Calculate polarizabilties flag

};

extern "C"
{
extern void cbincf_(double* bin_cf, double* fact);
extern double ssz_(int* nn1, int* nn2, int* ll1, int* ll2, int* m,
				   double* amu, double* bmu, double* fact, double* bincf,int* lg);  
extern double gint1_(double* exp1,double* exp2, int* n1, int* n2, double* r, 
					 double* fact,double* bincf);
}

#ifdef HAQCHEM_CPP

bool   qc_db_init_flag=false;
HaPseudoPotDB pseudo_db;
HaAtBasDB bas_db;

#else

extern bool   qc_db_init_flag;
extern HaPseudoPotDB pseudo_db;
extern HaAtBasDB bas_db;


#endif // end if HAQCHEM_CPP


#endif // end if !defined(HAQCHEM_H) 
