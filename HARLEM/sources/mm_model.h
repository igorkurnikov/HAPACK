 /*!  \file mm_model.h

    Classes to setup Molecular Mechanics Model

    \author Igor Kurnikov 
    \date 2010-

*/
#ifndef MM_MODEL_H
#define MM_MODEL_H

#include "haatgroup.h"
#include "mm_params.h"
#include "mm_force_field.h"

#if defined(SWIG) 
%template(vector_MMDihedral)  vector<MMDihedral>; 
%template(vector_VdWContact)  vector<AtomContact>;
%template(vector_MMBond)      vector<MMBond>;	
#endif

class AmberMMModel;
namespace mort
{
	class molecule_t;
};

class MolMechModel
{
public:
	MolMechModel(HaMolMechMod* p_mm_mod_new);
	virtual ~MolMechModel();

	int InitModel(const ForceFieldType& = MMForceField::ff_type_default); //!< Initiate Molecular Mechanics Model with force field name ff_name_par
	int UpdateModel(); //!< Initialize MM Model with default FF if it is not initialized or set internal model arrays if needed after model modification  
	int Clear(); //!< Clear model
	int SetStdParams(); //!< Set Standard values for parameters

	int UpdateConstraints(); //!< Update internal arrays for constraints for internal MM simulations
	int UpdateConstraints_2(); //!< Second part of Update of internal arrays for constraints for internal MM simulations ( for MPI )
	int ClearMortModel(); //!< Clear MORT molecule(model) associated with Molecular Mechanics model
	
	void Bcast(MPI_Comm& comm); //!< Broadcast Class data over MPI Communicator 

	virtual int SaveXMLToStream(std::ostream& os, const harlem::SaveOptions* popt = NULL ) const; //!< Save model data to stream in XML format   
	
	MolSet* GetMolSet() { return pmset; }
	const MolSet* GetMolSet() const { return pmset; }

	ForceFieldType ff_type;  //!< Force Field Type (AMBER94, CHARMM22 etc..) 
	int to_init_mm_model;    //!< Flag to intialize MM Model 

	int IsAmoebaFF() const; //!< Is AMOEBA Force Field is set for the model

    enum PARAM_SET_METHOD { NOT_SET = 0, SET_DEFAULT = 1, SET_FF_FIELD = 2, SET_RES_TEMPL = 3, SET_SPEC = 4 }; //!< different parameter set methods

	int GetNA() const { return Atoms.size(); } //!< Get Number of Atoms

	vector<HaAtom*> Atoms;                        //!<  Atoms 
	set<MMBond, less<MMBond> > MBonds;            //!<  Valence Bonds
	set<MMValAngle, less<MMValAngle> > ValAngles; //!<  Valence angles
	vector<MMDihedral>   Dihedrals;               //!<  Dihedrals
	vector<MMDihedral>   ImprDihedrals;           //!<  Improper dihedral angles
		
	vector<AtomSet> excluded_atom_list;    //!< Excluded atom list: atoms for which non-bonded calculations are not computed 
	vector<AtomSet> nonbond_contact_list;  //!< Non-bonded contacts lists: atoms for which non-bonded interactions are computed

	AtomIntMap& GetAtIdxMap(int recalc = FALSE); //!< Get a map of HaAtom* to indexes in Atoms array ( 0-based ). Optionally recalculate. 
	AtomIntMap  at_idx_map; //!< Map of HaAtom* to indexes in Atoms array

	MMBond*     GetMMBond(HaAtom* pt1, HaAtom* pt2); //!< get bond between atoms
	MMValAngle* GetValAngle(HaAtom* pt1, HaAtom* pt2, HaAtom* pt3); //!< get valence angle between Atoms
	MMDihedral* GetDihedral(HaAtom* pt1, HaAtom* pt2, HaAtom* pt3,HaAtom* pt4); //!< get dihedral angle for atoms
	MMDihedral* GetImprDihedral(HaAtom* pt1, HaAtom* pt2, HaAtom* pt3, HaAtom* pt4); //!< get dihedral angle for atoms

	int SetMMBond(HaAtom* pt1,HaAtom* pt2,double r0,double fc,int set_type = NOT_SET); //!< Set valence bond 
	int SetValAngle(HaAtom* pt1,HaAtom* pt2,HaAtom* pt3,double a0,double fc,int set_type = NOT_SET); //!< Set valence angle
	MMDihedral* AddImprDihedral(HaAtom* pt1, HaAtom* pt2, HaAtom* pt3, HaAtom* pt4);   //!< Add improper angle

	//! map of Improper angles arranged by the residue of the 3rd atom
	multimap<unsigned long, unsigned long, less<unsigned long> > res_impr_dih_map; 

	void SetUseMortLib( int set_par ); //!< Set Parameters of the model Using MORT library 

	int build_nb_contact_list_flag;  //!< Flag to build non-bond contact list during initialization of the module
	int init_charges_flag;           //!< if true, attempt to set charges on the atoms during initialization as if they are zero
    int setup_params_from_mort_flag; //!< Setup MM model using Mort library functions

	int SetCoarseGrainedOPEPParams(); //!< Set MM parameters for coarse-grained model of AA, OPEP ver. 3.0 Maupetit et al 2007
	int SetCoarseGrainedDNAParams();  //!< Set MM parameters for coarse-grained model of DNA, De Pablo et al 
	int SetCoarseGrainedAAParams();   //!< Set MM parameters for coarse-grained model of AA, Bond and Sansom 2006 jose August, 2008

	int SetStdValParams();       //!< Set Standard parameters for valence bonds, angles and dihedrals
	int SetStdVdWParams();       //!< Set Standard Van-der-Waals parameters

	int  AddAtomsToExcludedAtomList(HaAtom* aptr1, HaAtom* aptr2, PtrIntMap& atoms_idx); //!< Add Atoms to excluded atom list taking into account atom order
	bool BuildExcludedAtomList();   //!< Build the excluded atom list 
	bool BuildNonBondContactList();  //!< Build the non-bonded contact list
        
	int BuildGrpGrpExcludedList(AtomContainer* group1, AtomContainer* group2); //!< Build a non-bonded excluded atom list for atoms that belong to two groups
	int BuildGrpGrpNonBondList(AtomContainer* group1, AtomContainer* group2);  //!< Build a non-bonded contact list for atoms that belong to two groups
	
	int GetResImprAngles(HaResidue* pres, list<MMDihedral*>& res_idih_list); //!< Get all improper angles with a third atom which belong to the 
	                                                                         //!< certain residue
	                                                                         //!< return the number of improper dihedrals found

	int Set14interDihFlags(); //!< Set flags for dihedrals not to calculate 1-4 interactions for atoms participating in bonds or valence angles or several dihedrals
	
	int SetBoundaryBox(double offset = 0.0); //!< Set periodic boundary box boundaries are placed offset bohrs from VdW surface of the molecular set

	bool CalcNonBondPt(const HaAtom* pt1, const HaAtom* pt2, double& vdw_at_ene, double& el_at_ene);  //!< Calc nonbond energy between two Atoms

//! \name Potential functions:
//@{
public:
	
	void SetOmitInteractionsParam(const OmitInteractionsParam& omit_interactions_new); //!< Set Option to omit certain atom interactions
	void SetMMElectrMethod( const MMElectrMethod& electr_method_new ); //!< Set method to treat electrostatic interactions

	void SetCalcDirectInter( bool set_par = true );     //!< Calculate direct interaction terms
	void SetCalcRecipSpaceInter( bool set_par = true ); //!< Calculate reciprocal space interaction terms  
	void SetCalcSelfInter( bool set_par = true );       //!< Calculate self interaction terms 
	void SetCalcAdjustInter( bool set_par = true );     //!< Calculate adjust interaction terms 
	void SetCalcVdWInter( bool set_par = true );        //!< Calculate Van-der-Waals interaction terms 
	void SetCalcInducedInter( bool set_par = true );    //!< Calculate Electrical Induced interaction terms 
	 
	double GetScale14Electr() const;       //!< Get scale coef for 1-4 coulomb inter.
	double GetScale14VdW() const;          //!< Get scale coef for 1-4 coulomb inter.
	double GetNBCutDist() const;           //!< Get cutoff distance for non-bond interactions
	double GetDielConst() const;           //!< Get dielectric constant
	double GetIonStrength() const;         //!< Get ionic strength
	
	int GetDipoleScfIterMax()  const; //!< Get Maximal number of SCF iterations of induced dipoles 
	double GetDipoleScfTol()   const; //!< Get Tolerance for convergence of induced dipoles
	double GetEEDsumCut()      const; //!< 
	double GetEEDampedCut()    const; //!< 
	double GetSorCoef()        const; //!< 
	double GetTholeExponCoef() const; //!< Get Thole's damping factor for interaction of induced dipoles
	double GetVdwTaper()       const; //!< 

	void SetScale14Electr( double scale_14_electr_new ); //!< Set scale coef (divide by) for 1-4 coulomb inter.
	void SetScale14VdW( double scale_14_vdw_new ); //!< Set Scale coef (divide by) for 1-4 vdw inter
	void SetScale14Inter( double scale_14); //!< Set Scale coef (divide by) for 1-4 Non-bonded interactions ( electr and VdW
	void SetNBCutDist( double nb_cut_dist_new );    //!< Set cutoff distance for non-bond interactions
	void SetDielConst( double diel_const_new );     //!< Set dielectric constant
	void SetIonStrength( double ion_strength_new ); //!< Set ionic strength  

	void SetDipoleScfIterMax( int dipole_scf_iter_max_new ); //!< Set Maximal number of SCF iterations of induced dipoles 
	void SetDipoleScfTol( double dipole_scf_tol_new );   //!< Set Tolerance for convergence of induced dipoles 
	void SetEEDsumCut( double ee_dsum_cut_new );                 
    void SetEEDampedCut( double ee_damped_cut_new ); 
	void SetSorCoef( double sor_coef_new );
    void SetTholeExponCoef( double thole_expon_coeff_new );  //!< Set Thole's damping factor for interaction of induced dipoles
	void SetVdwTaper( double vdw_taper_new );  //!<  

	OmitInteractionsParam omit_interactions;   //!< Omit certain atom interactions     to eliminate???

	MMElectrMethod electr_method; //!< method to treat electrostatic interactions

	double scale_14_electr; //!< Scale coef (divide by) for 1-4 coulomb inter.
	double scale_14_vdw;    //!< Scale coef (divide by) for 1-4 vdw inter. 
	double nb_cut_dist;     //!< cutoff distance for non-bond interactions
	                          
	double diel_const;        //!< Dielectric constant
    double ion_strength;      //!< Ionic Strength 

	int nonb_list_flag;        //!< Flag to generate non-bonded atom list, NTNB in AMBER
	int neutral_end_hydr_flag; //!< neutralize end hydrogen (needed for DNA) ICHDNA in AMBER

	enum {CALC_VDW_NO = 0,CALC_VDW_NORMAL,CALC_VDW_NO_ATTRACT} calc_vdw_flag; //!< Flag to compute Van-der-Waals interactions
	int calc_electr_flag;     //!< Flag to compute Electrostatic interactions

	int dipole_scf_iter_max; //!< Maximal number of SCF iterations of induced dipoles 
	
	double dipole_scf_tol;   //!< Tolerance for convergence of induced dipoles 
	double ee_dsum_cut;      
    double ee_damped_cut; 
	double sor_coefficient;
    double thole_expon_coeff;  //!< Thole's damping factor for interaction between polarizable atoms
	double vdw_taper;
    
	enum {NO_SOFT_REPULSION = 0, SOFT_REPULSION_NO_HBOND = 1, SOFT_REPULSION_ALL = 2} soft_repulsion; //!< Soft repulsions option, ISFTRP in AMBER
	double soft_repulsion_const; //!< Constant K for in K(r0^2- r^2)^2 in soft repulsion
	                             //!< term, RWELL in AMBER
//@}

//! \name parameters for Generalized Born Calculations
//@{
public:
	enum GenBornParamType  
	{
		BONDI_GB_PARAM = 0,      
		AMBER6_BONDI_GB_PARAM = 1, //!< Amber 6, JACS 122:2489 (2000)
		MOD_BONDI_GB_PARAM = 2,    //!< Biopolymers 56: 275 (2001)
		HUO_KOLLMAN_GB_PARAM = 3,
		JAYARAM_GB_PARAM     = 4,  //!< Jayaram et al. 'GB'
		MOD_JAYARAM_GB_PARAM = 5,  //!< Jayaram et al. 'MGB'
		HN_MOD_BONDI_GB_PARAM = 6
	} gb_param_type; //!<  type of Generalized Born calculations from enum GenBornParamType
 	static double GetGBAtomRad(HaAtom* aptr, int gb_param_type = MOD_BONDI_GB_PARAM);
	static double GetGBAtomScreening(HaAtom* aptr, int gb_param_type = MOD_BONDI_GB_PARAM);
//@}

//! \name Restraints and Belly dynamics control:
//@{
	int SetMovingAtoms( std::string mov_atom_array_name); //!< Freeze positions of all atoms except of the specified atom array (belly dynamics) 
	int SetMoveAll();                                     //!< Unfreeze all atoms

	AtomGroup* GetMovingAtoms(); //!< Get Atom Array of atoms allowed to move (=NULL if all atoms allowed to move)

	std::string moving_atoms;    //!< name of the Atom Array of atoms allowed to move in the dynamics =ALL_ATOMS for normal runs 

	int SetRestrainedAtoms( const char* restr_atom_group_name); //!< set name of atom group of harmonically restrained atoms 
	AtomGroup* GetRestrAtoms(); //!< Get Atom Group with harmonically restrained positions
	int SaveAtomRestrArbalestIndForm( std::string restr_desc_fname, std::string restr_list_fname); //!< Save 
	int SetRestrRefCrdFromXYZFile( const char* ref_crd_file_name_new ); //!< Set Reference coordinates for restrained atoms from TINKER XYZ file  
	int SetRestrRefCrdFromStr( const std::string& crd_str ); //!< Set reference coordinates for restrained atoms from TINKER XYZ file ( white space separated coordinates)
    int SetAtomRestrForceConst(double restr_const_new); //!< Set force constant for atom position restraints (kcal/mol/Ang^2)
	double GetAtomRestrForceConst(); //!< Get value of the force constant for atom position restraints (kcal/mol/Ang^2)

	enum RestrRefCrdType { RESTR_REFC_CURRENT_CRD = 0, RESTR_REFC_XYZ_CRD_FILE, RESTR_REFC_XYZ_CRD_STR } restr_ref_crd_type; //!< Restraints reference coordinates type  
	std::string restrained_atoms;    //!< name of the atom group of restrained atoms in Constrained Dynamics
	Vec3DValArray restr_ref_coords;  //!< Reference coordinates of restrained atoms 

	int SetDistConstrFromFile( const char* constr_file_name ); //!< Set Harmonic Atom Distance constraints from file 
	int GetNumHarmConstr() const; //!< Get the number of harmonic constraints

	std::vector<AtomContact>  DistConstraints;  //!<  Atom Distance Constraints as VdW energy terms, Harmonic terms etc 

protected:
	double atom_restr_const;        //!< force constant for atom position restraints (kcal/mol/Ang^2)
public:
	bool SetHBondRestraints(double force_const); //!< put harmonic constraints on the H-Bonds 
	int AddHarmConstraint(HaAtom* atom1,HaAtom* atom2, double eq_dist, double force_const); //!< Add Harmonic Constraint
	int SetHarmConstraint(HaAtom* aptr1,HaAtom* aptr2, double eq_dist, double force_const); //!< Set Harmonic Constraint (Update if already exists)
	int SetHarmConstraint(const std::string& at_ref1, const::std::string& at_ref2, double eq_dist, double force_const);

	bool ClearConstraints(); //!< Delete All Harmonic constraints


//@}

//! \name Water cap treatment:
//@{
	int water_cap_flag; //!< Flag to control Cap option, correspond to IVCAP in AMBER
	int cap_atom_num;   //!< the Cap atom pointer(serial num) MATCAP in AMBER
	double cap_fconst;  //!< the force constant for the CAP restraint pot. FCAP in AMBER
	                    //!< default (1.5 kcal/Ang*Ang ??)
//@}

//! \name Particle Mesh Ewald:
//@{
	void SetFFTGridsPerAng( double fft_grids_per_ang ); //!< Set Number of FFT grids per Ang in PME calculations

	int pmesh_ewald_flag;   //!< Turns on the Particle Mesh Ewald, IEWALD in AMBER
	
	int pme_grid_nx;   //!< size of charge grid in PME calculations, if 0 set from the periodical box size 
	int pme_grid_ny;
	int pme_grid_nz;   //!< NFFTX,NFFTY,NFFTZ in AMBER

    int pme_spline_order;      //!< the order of the B-spline interpolation in PME, SPLINE_ORDER in AMBER
    int pme_neutralize_system; //!< Force neutralization of the unit cell, ISCHARGED in AMBER (reverse 0 - 1)
    int pme_verbose;           //!< turn on voluminous output of info about PME run, VERBOSE in AMBER
	int exact_ewald_flag;      //!< Calculate Ewald calculation is run, EXACT_EWALD in AMBER
	int vdw_correction_flag;   //!< if =1 Apply analytical correction to VdW interaction energy beyond cutoff distance  
	double pme_dsum_tol;       //!< The width of the direct sum part of the Ewald sum, DSUM_TOL in AMBER
	double pme_ew_coeff;       //!< Ewald Coefficient
	double pme_eedtbdns;       //!< PME eedtbdns parameter  
	double skin_nb;            //!< size of intermediate region beyond non-bonded cutoff distance used to track non-bonded list rebuild
	double fft_grids_per_ang;  //!< number of PME FFT grid points per Ang 
	int subtract_avg_force_flag; //!< Subtract average force due to analytic forces in PME (1 (subtract), 0(not subtract), -1 (subtract for MD) not subtract for MIN

	void SetPMECoef( double pme_ew_coeff_new ); //!< Set PME coefficient (if = 0 it will be computed from pme_dsum_tol )

//@}

	AmberMMModel* p_amber_model; //!< MolMechModel in AMBER representation
	HaMolMechMod* p_mm_mod;
	MolSet*     pmset;         //!< MolSet for which molecular model is built - can contain extra atoms(force centers) or miss some atoms etc...
	mort::molecule_t* p_mort_model; //!< Molecular Mechanics model in the MORT library format
	
	int UpdateDataFromFort();

	int GetIndDipolesFromFort(); //!< Fill Array of induced dipole from fortran
	int GetMultipolesFromFort(); //!< Fill Array of multipoles from fortran

	void PrintIndDipoles(); //!<  Print the table of Induced Dipoles to console
	void PrintMultipoles(); //!<  Print the table of Multipoles to console
	void PrintTotMultipoles(const Vec3D* pt_orig = NULL ); //!< Print Multipoles of the whole molecular system to console

	Vec3D GetTotIndDipole(const Vec3D* pt_orig = NULL);  //!< Get total induced dipole moment (electron * Ang)
	Vec3D GetTotIndDipole1(const Vec3D* pt_orig = NULL); //!< Get total induced dipole moment from ind_dip_d (electron * Ang)
	Vec3D GetTotIndDipole2(const Vec3D* pt_orig = NULL); //!< Get total induced dipole moment from ind_dip_p (electron * Ang)
	Vec3D GetTotDipole(const Vec3D* pt_orig = NULL);    //!< Get total dipole moment         (electron * Ang)
	HaVec_double GetTotQpole(const Vec3D* pt_orig = NULL);     //!< Get total quadrupole moment     (electron * Ang * Ang)

	HaVec_double ind_dip_p; //!< Induced dipole values p? 
	HaVec_double ind_dip_d; //!< Induced dipole values d?
	HaVec_double global_multipole; //!< Multipole array

	friend class MolMechDlgWX;
	friend class MMDriverAmber;
	friend class AmberMMModel;
	friend class HaMolMechMod;
	friend class HaMolMembraneMod;
	friend class HaInterMolMod;
};

#endif // end if !defined(MM_MODEL_H) 
