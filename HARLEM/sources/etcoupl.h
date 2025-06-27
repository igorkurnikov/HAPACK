/*! \file etcoupl.h
 
    Classes to calculate Electron Transfer donor/acceptor coupling in HARLEM.
 
    \author Igor Kurnikov
    \date 1998-2002

*/


#ifndef ETCOUPL_H
#define ETCOUPL_H

#include "halocorb.h"
#include "hacompmod.h"
#include "haatgroup.h"

class HaQCMod;
class HaMolecule;
class MolSet;
class HaAtom;
class HaMat_double;
class LinCombOrb3D;


class ETNode;

class ETEdge
{
public:
	ETEdge();
	ETEdge(int new_inode1, int new_inode2);
	~ETEdge() { }

	bool operator==(const ETEdge& rhs) const;
	bool operator< (const ETEdge& rhs) const;

	double coupling; // value of the log coupling for the connection
//	double dist; // distance in ang

	int inode1;  // indexes of nodes in the vector of node ETCouplMod
	int inode2;

protected:
	void SetDefaultParam();
};


class ETPath
//!< class for coupling path between two nodes 
{
public:
	ETPath();
	virtual ~ETPath() { }

	bool clear();
	bool empty() { return trace.empty(); }

	double coupling; // coupiling value
	std::list<int> trace; // list of node indexes of the best path

};

class PathStep
{
public: 
	PathStep();
	PathStep(double new_coupling, int new_destination, int new_source);
	virtual ~PathStep();
	
	bool operator < (const PathStep & rhs) const;

	double coupling;     //!< coupling to donor in PATHWAYS calculations
	int destination;     //!< index of Destination   
	int source;          //!< and source nodes in the vector of nodes of ETCouplMod

};

const int BEST_PATH = 0, COUPL_MAP = 1; //!< for ETCouplMod::calc_type
const int HAM_S_DIP_TR = 0, HAM_TR = 1; //!< for ETCouplMod::ham_trunc_type
enum REDOX_ORB_TYPE {REDOX_ORB_DONOR = 0, REDOX_ORB_ACCEPTOR};

//! The module to perform computations of donor/acceptor ET coupling
class ETCouplMod: public HaCompMod
{
public:

	ETCouplMod(MolSet* new_phost_mset = NULL);
	~ETCouplMod();

	HaQCMod* GetQCMod() { return ptr_qc_mod; }  //!< Get Quantum Chemical module associated with this module

// Virtuals of HaCompMod:

	virtual void SetDebugLevel(int new_debug_level);
	virtual int OnDelAtoms(AtomContainer& del_atoms);  //!< Modify module content to react to deleted atoms 

	bool Clear(); //!< delete internal arrays and reset parameteres 
	              //!< to default values

//! \name PATHWAYS Calculations
//@{
	int pathways_calc_type;  //!< Type of PATHWAYS calculations BEST_PATH = 0, COUPL_MAP = 1;
	
	bool path_coupl_calc(); //!< Perform PATHWAYS calculations (best pathway or ET Coupling Map) 
	bool select_important(double thresh); //!< Find atoms that are on all pathways that are within certain threshold of the value
	                                      //!< of the strongest coupled path

	bool calc_intermol_path_coupl(); //!< Calculate PATHWAYS coupling between two molecules  
	
	bool InitiatePathwaysGraph();           //!< Initiate the atom connection graph for PATHWAYS calculations
	double CalcAtomContactCoupling(HaAtom* aptr1, HaAtom* aptr2);  //!< Set Pathways coupling (log) for an Atom-Atom Contact
	void ClearPathwaysGraph();              //!< Clear Graph for PATHWAYS Calculations
    int ColorMolSurfETCoupl();              //!< Color Molecular Surface by ET Coupling to Donor

	AtomGroup nodes; //!< Set of nodes of the graph in PATHWAYS calculations
	VecPtr    edges; //!< Vector of list of edges of the nodes

	bool pathways_graph_init_flag; //!< flag to indicate that the connectivity graph should be recomputed before PATHWAYS calculations are performed

	bool m_hbond_paths_flag; //!< Flag to use special functional form for H-Bond contact  

	std::vector<PathStep> best_path; //!< Vector of step that make best coupled path in PATHWAYS calculations
	std::vector<PathStep> coupl_map; //!< the map of best path coupling to the donor

    double pw_nb_decay;     //!< PATHWAYS non-bond exponential distance decay parameter
	double pw_nb_min_dist;  //!< PATHWAYS non-bond minimal distance
	double pw_hb_decay;     //!< PATHWAYS H-bond exponential distance decay parameter
	double pw_hb_min_dist;  //!< PATHWAYS H-bond minimal distance
	double pw_ln_cov_decay; //!< PATHWAYS ln of covalent bond decay parameter

	double best_path_coupl; //!< last calculated best path value
	int log_calc_result;    //!< Flag to log calculation results
    
    double pw_nb_decay_intermol; //!< non-bond exp decay param for intermolecular contacts, if ==0 then use pw_nb_decay 
	double nb_dist_limit;        //!< maximum distance at which non-bonded contact is computed

    int rebuild_mol_coupl_map;

	PtrDoubleMap mol1_coupl_map;
	PtrDoubleMap mol2_coupl_map;
//@}

//! \name Distance-based estimates of Donor/Acceptor interactions
//@{
	double calc_edge_dist(); //!< Calculate edge-to-edge distance between donor and acceptor groups 	
	bool DuttonModelCalc();  //!< Calculations of H_DA using Dutton Model
// Dutton Model Parameters:     
	double rho;      //!< Dimesionless average density of atoms in the donor/acceptor direction
	double beta;     //!< Average decay exponent of ET rate 
	double dim_less_coupling; //!< Dimesionless ET coupling
	double max_rate; //!< Maximal ET rate (activationless -  DeltaG = -lambda)
//@}

//! \name Quantum Chemical Calculation of donor/acceptor interactions
//@{
	bool CalcGFDonAccOrb();     //!< Calculate Green-Function Matrix elements from MOs
	bool CalcGFDonAccOrbHeff(); //!< Calculate Green-Function Matrix elements from Heff
    int  PrintOvlpElem();    //!< Print Overlap Matrix elements for selected orbitals
	int  PrintHeffElem();    //!< Print Heff Matrix elements for selected orbitals

    bool SetDAdipoleMat();                   //!< Set the matrix of interaction of the molecule with the donor/acceptor electric field
	double GetDAfield() const { return DA_field; } //!< Get the value of the donor/acceptor electric field
	bool SetDAfield(double field);           //!< Set the value of the donor/acceptor electric field
	void SetTunEne(double tun_ene_new) { tun_ene =  tun_ene_new; } //!< Set tunneling energy 
	double GetTunEne() const { return tun_ene; }                   //!< return current tunneling energy
	
	bool CalcHDAEneSplit(); //!< Find Donor/Acceptor couplings between donor and acceptor localized eigenvectors using minimal energy splitting algorithm
    int AddRedoxOrbFromEigVec(const HaVec_int& mo_idx); //!< Add truncated Eigen Vectors (specified indexes of MO) to the the list donor/acceptor orbitals
	int GetRedoxOrbsFromFrag( MolSet* pfrag ); //!< Get Redox Orbitals from the fragment donor/acceptor orbitals
	bool FindRedoxOrbsOvlpEigVecs( StrIntMap& lbl_idx_map, StrDoubleMap& ovlp_val_map,  REDOX_ORB_TYPE redox_orb_type = REDOX_ORB_DONOR ); //!< Find indexes of eigen vectors that donor/acceptor orbitals overlap most
	bool FindRedoxOrbSpaceOvlpEigVecs( StrIntMap& lbl_idx_map, StrDoubleMap& ovlp_val_map, HaVec_double& eigv_space_max_ovlp_val, REDOX_ORB_TYPE redox_orb_type = REDOX_ORB_DONOR ); //!< Find indexes of eigen vectors overlapping most with the space of donor/acceptor orbitals
    HaVec_int FindDonAccEigVecs(int idx_don, int idx_acc); //!< find indexes (0-based) of eigenvectors maximally overlaped with donor/acceptor orbitals idx_don and idx_acc 
	bool RotateRedoxOrb( const HaMat_double& rot_mat ); //!< Transform expansion coefficients of donor/acceptor orbitals for rotation defined by rot_mat  

	bool CreateEigVecContour(int idx, double flvl, int grid_size); //!< build an isosurface of the eigenvector
	bool PrintEigVecCoef(int idx);       //!< Print to STDOUT expansion coefficient of the eigenvector idx 
	bool ScanEigEneField(int first_eig_val, int last_eig_val, double ifield_val, double ffield_val, 
						 double step_val); //!< Scan eigenvalues versus donor/acceptor field
		
	bool CalcHDAfromGF();             //!< Compute donor/acceptor coupling using donor/acceptor Green-Function  
	bool CalcHDAPert(double& hda_coupl);  //!< Compute donor/acceptor coupling using perturbation formula

	int CopyEigVecsFromMO(); //!< Copy Eigenvectors and Eigenvalues from loaded MO coefficients and energies

	bool RecalcHeff(); //!< Recompute the effective 1e Hamiltonian matrix
	bool DiagHeff();   //!< diagonalize the effective 1e Hamiltonian matrix
	int ZeroLongInter(double cutoff); //!< Zero interactions in the hamiltonian beyond a cutoff distance
	int SaveHeffXml(FILE* file_out); //!< Save one-electron hamiltonian to a file
	int LoadFragmHeffXml(FILE* file_inp); //!< Load one-electron hamiltonian of the fragment 

	HaMat_double& GetActBasOvlpMat(); //!< Get overlap matrix of active basis set

	int ham_trunc_type; //!< hamiltonian truncation type (heff_mat,ssl,DA_dipole) or just heff_mat

    HaMat_double heff_mat; //!< Effective one-electron Hamiltonian of the system 
	HaMat_double ssl;      //!< Overlap matrix in local orbitals
	LinCombOrb3D eigv;     //!< Eigenstates of the Effective 1e hamiltonian
	HaVec_double enel;     //!< Eigenvalues of the Effective 1e hamiltonian

    HaMat_double heff_pert_mat; //!< matrix of estimated perturbation of heff matrix elements (Used in D&Q calculations)
    int use_pert_mat;           //!< flag to use perturbation matirx to costruct effective hamiltonian

	LinCombOrb3D donor_orbs; //!< donor orbitals
	LinCombOrb3D acc_orbs;   //!< acceptor orbitals

	HaVec_int ieig_don;    //!< Eigen vector indexes corresponding to the donor orbitals
	HaVec_int ieig_acc;    //!< Eigen vector indexes corresponding to the acceptor orbitals

	HaMat_double don_acc_gf; //!< donor/acceptor orbitals green-function matrix

	HaMat_double el_field_min;  //!< values electrical field minimizing energy spliting
	HaMat_double da_coupl_val;  //!< donor/acceptor coupling values between donor & acceptor orbitals 

	HaMat_double extern_field; //!< Matrix of the of the external field to be added to effective hamiltonian

	HaMat_double DA_dipole; //!< Matrix of Electrical Dipole in the direction from the donor
	                        //!< to acceptor (defined on active local orbitals)

	double DA_field;  //!< initial value of the electrical field in the donor-acceptor direction

	int set_dab_huck_inter;  //!< Flag for interactions between donor/acceptor
	                         //!< and bridge orbitals to be calculated using Extended Huckel 

	double tun_ene;          //!< tunneling energy for Green function and perturbation calculations

    HaQCMod*  ptr_qc_mod;

//@}

//! \name Divide-and-Conquer Construction of the model electronic hamiltonian
//@{
	bool PutSubMatToDB();           //!< Save group-group submatricies of the effective hamilonian to the database
	bool GetSubMatFromDB();         //!< Get group-group submatricies of the effective hamilonian from the database

	bool PrintProtectMat() const; //!< print the diagonal of the protection matrix

	HaMat_double protect_mat;  //!< Matrix of protection numbers characterizing how well particular block (group-group)
	                           //!< of the hamiltonian matrix is represented by patched hamiltonian obtained in Divide-and-Conquer
	                           //!< procedure

	std::string db_file_name;      //!< The name of the database file of the effective hamiltonian submatrices in D&Q calculations
//@}


};



#endif /* !ETCOUPL_H */
