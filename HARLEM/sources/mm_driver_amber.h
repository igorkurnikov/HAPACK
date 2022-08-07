/*!  \file mm_driver_amber.h

    Molecular Mechanics simulations using AMBER package  

    \author Igor Kurnikov 
    \date 2008-

*/

#if !defined MM_DRIVER_AMBER_H
#define MM_DRIVER_AMBER_H

#include "hamolmech.h"

class TimerAmber;
class MMBond;
class MMValAngle;
class MMDihedral;
class AmberMMModel;
class AtomFFParam;


class MMDriverAmber : public MMDriver
{
public:
	MMDriverAmber(HaMolMechMod* p_mm_mod_new);
	virtual ~MMDriverAmber();

	virtual std::string GetClassName() { return "MMDriverAmber"; }

	virtual int CalcEnergy() { return FALSE; } //!< Calculate energy of the system and save resluts to p_mm_info member of p_mm_mod

	TimerAmber*  p_tm;    //!< Timer object associated with the driver
	AmberMMModel* p_amber_model; //!< Molecular Mechanics model in AMBER representation
	
public:

	std::string sander_exe_fname;     //!< AMBER(sander) executable file name

	int amber_version;  //!< version of the AMBER program  
    int min_md_cpp_flag; //!< flag to do MIN and MD using C++ main functions 

	MPI_Comm   driver_mpi_comm;  //!< MPI Communicator associated with the driver
	MPI_Group  driver_mpi_group; //!< MPI Group        associated with the driver

	int SaveAmberRunFile();  //!< Save script to run AMBER on UNIX 
	int SaveAmberInpFile();  //!< Save AMBER Input Parameters File
	int SaveAmberTopFile();  //!< Save AMBER topology file (*.top) 
	int SaveAmberCrdFile();  //!< Save file with initial coordinates
	int SaveAmberRstFile(const char* fname); //!< Save current geometry in AMBER restart format
	int SaveAllInpFiles();  //!< Save all AMBER input files
	int DeleteOutputFiles(); //!< Delete AMBER Output Files 

	int RunAmberProg(int sync); //!< Run External AMBER (Sander or PMEMD) executable ( if sync == TRUE - no return before AMBER process exited)

	int LoadAmberRestartFile(const char* rst_file_name); //!< Load coordinates from the AMBER restart file
	int LoadAmberMDInfoFile();  //!< Load info from AMBER mdinfo file

// from parallel_dat_mod:
	int master;        //!< parameter indicating is this master node
	int numtasks;      //!< Number of MPI processes allocated for the problem ( = 1 - serial run ) 
	int mytaskid;      //!< MPI rank of the process ( = 0 for master)

	int my_atm_cnt;        //!< Number of atoms assigned to the node  MOVE TO AMBER MODEL 
	HaVec_int gbl_my_atm_lst;    //!< List of indexes atoms(1-based) assigned to the node
	HaVec_int gbl_atm_owner_map; //!< map of MPI task ids that atoms belong to  

	int SetMPICommAllProcs(); //!< Set MPI communicator of the driver to MPI_COMM_WORLD  
	
    void InitParallelDatMod(); //!< Init parallel_dat_mod parameters in PMEMD_LIB to be able call MPI functions

	void ResizeCrdVelFrcArrays(int natom); //!< Resize atm_crd, atm_vel, atm_last_vel and atm_frc arrays and initialize with 0.0

	void BcastCrd(MPI_Comm& comm); //!< Broadcast atm_crd Array
	void BcastVel(MPI_Comm& comm); //!< Broadcast atm_vel Array
	void BcastFrc(MPI_Comm& comm); //!< Broadcast atm_frc Array
	void BcastPBox(MPI_Comm& comm); //!< Broadcast Periodical Box

	void SetupMasterNode(); //!< Perform FORTRAN MM LIB Setup specific for Master Node
	void SetupShakePars();  //!< Setup Amber Shake Params
	void SetupPME();  //!< Setup Amber PME calculations
	void SetupGB();   //!< Setup Amber GB (Generalized Born) calculations
	void ValidateCntrParams(); //!< Validate compatibility of control parameters for MD/Min simulations
	void PrintMDCntrData();    //!< Print MD/Min control Data to Output File
	void InitPMEParams();      //!< Initialize parameters of PME model  
	void InitAddMDCtrlParams();      //!< Setup Additional control parameters dependent on chosen MM Simulation method
	void BCastAddMDCtrlParams(MPI_Comm& comm);  //!< BroadCast Additional MD Control Parameters
	void SetAddMDCtrParamsFortran();             //!< Set Additional MD control parameters to Fortran structures

	void AllGatherVec(HaVec_double& vec); //!< Gather Atom properties Vector (coord,vel,frc...) on MPI Nodes 
	
	void PrintLogMDOUT(const char* format, ... ); //!< Print to MDOUT
	void FlushMDOUT(); //!< Flush MDOUT

//    void SetEnergyParamsArrays(); //!< Set Structures Defining Energy Functional Model in PMEMD 
	static int stack_limit;
	void ResetStackLimits(); //!< Reset the stack limits if possible
	int InitSimulationsStep2();   //!< Init structures and output files for MD and Minimization simulations 
    int InitCtrlParams();   //!< Init control parameters	

	int SaveAmberTopToStream(std::ostream& os);  //!< Save AMBER TOPOLOGY to stream
	void SetAtomCrdToInternalArrays();   //!< Set Atom Coordinates  of MolMechnMod to atm_crd and pbc_box arrays of the MMDriverAmber module
	void GetAtomCrdFromInternalArrays(); //!< Get Atom Cooordinates of MolMechnMod from atm_crd and pbc_box arrays of the MMDriverAmber module

	void SaveModelToFortran();    //!< Set All Fortran Structures Describing the model
	void SetPBoxDataToFortran();    //!< Set Periodical Boundary Data to Fortran Structures
	void GetPBoxDataFromFortran();  //!< Get Periodical Boundary Data from Fortran Structures

	int CheckForStop(); //!< Check if Stop of calculations was requested
	int OpenOutputFiles();      //!< Open output files for simulations
	int CloseOutputFiles();     //!< Close output files for simulations

	void GBForce(int atm_cnt, HaVec_double& crd, HaVec_double& frc,
							      HaVec_double& si, int ncalls);     //!< Compute Force and energy in GB approximation
	void PMEForce(int atm_cnt, HaVec_double& crd, HaVec_double& vel, HaVec_double& frc,
							       int new_list, HaVec_double& si );  //!< Compute Force and energy in PME approximation

	void CalcForceAndEne(int new_list, int ncalls); //!< Calculate Potential Energy and Forces  

	void CalcKinEne(int only_cur_vel = FALSE ); //!< Calculate kinetic energies - eke (at current t), ekph (at t+1/2 dt), ekpbs( cross betwee atm_vel and atm_last_vel) 
                                                //!< if only_cur_vel = TRUE - do calculation only with atm_vel assuming the are current                                                                                        

	void IncrementVelAndCrd(int nstep, HaVec_double& crd, HaVec_double& vel, HaVec_double& frc, int& reset_velocities); //!< Increment Velocities and Coordinates According to current Forces and Velocities 
	
	void CalcPressure(HaVec_double& si); //!< Compute pressure using computed virial and kinetic energy 
	void ScaleCoordConstPress(HaVec_double& crd, HaVec_double& box, HaVec_double& si); //!< Scale Coordinates and periodic box at const pressure run, 
	                                                                                 // si - contains presssure components at [SI_PRESS_0],... 
   
	void PropagateVelHalfStepBack(HaVec_double& vel, HaVec_double& frc); //!< Propagate velocities a half step back using forces frc   

	void SaveVelToLastVel(int& all_vels_valid); //!< Save atm_vel to atm_last_vel

	void InitPBC();  //!< Initialize periodical box parameters
	void CalcNumDegFreedom(); //!< Compute the number of degrees of freedom for solvent and solute and total 
	void CalcCenMassVel(HaVec_double& crd, HaVec_double& vel, HaVec_double& xcm, HaVec_double& vcm, HaVec_double& ocm); //!< Calculate position and velocity of the center of mass
    void RemoveCOMVelocity(); //!< Subtract Center of Mass velocity from atomic velocities 
	void RemoveCOMVelAndResetCenter(int& all_crds_valid, int& all_vels_valid, Vec3D& sys_cnt); //!< Subtract Center of Mass velocity and reset center of the system
	void GetPosition(HaVec_double& crd, Vec3D& cnt ); //!< Compute center of the system cnt  
	void RePosition (HaVec_double& crd, Vec3D& cnt);  //!< Recenter the system at cnt
	void GrdMax(HaVec_double& frc,int& iatmax, double& fdmax); //!< get index coordinate with max value(fdmax) of frc
	void ZeroVelFrozenAtoms(); //!< Zero Velocities of frozen(belly) atoms
	void ZeroFrcFrozenAtoms(); //!< Zero Forces of frozen(belly) atoms

	void CollectCoords(int& new_list, int& nstep, int& collect_crds, int& all_crds_valid); //!< Collect coordinates from mpi nodes if needed
	void CollectVelocities(int& new_list, int& nstep, int& collect_crds, int& all_vels_valid); //!< Collect velocities from mpi nodes if needed 
	void ScaleVelConstTemp(); //!< Scale velocities to maintain constant temperature ref_temp	
	int CheckForNewNonBondList(); //!< Check if new non-bonded atom list is necessary to build on the next force calculation
	int CheckAllAtomMovement(int atm_cnt, HaVec_double& crd); //!< Check if any atom moved by more than 1/2 skin_nb compare to coordinates saved in 
	                                                              //!< gbl_atm_saved_crd and gbl_saved_box
	void RunMinMaster();      //!< Run Minimization on a master node 
	void RunMinSlave();       //!< Run minimization on a slave node  
	void RunMD();             //!< Run MD simulations 
	void CalcCurrEne();       //!< Calculate Energy and force for current coordinates

	void SetMMInfo(MMSysInfo& info, HaVec_double& si);  //!< Set MMSysInfo structure from sys_info[] array 
	void SetPrmTopIntFortran();   //!< Set Fortran Common /prm_top_int/
	void SetMDinCtrlDblFortran(); //!< Set Fortran Common /mdin_ctrl_dbl/ 
	void SetMDinCtrlIntFortran(); //!< Set Fortran Common /mdin_ctrl_int/
	void SetAddIntParsFortran();     //!< Set Additional Integer Fortran Parameters
	void SetPMEParsFortran();     //!< Set Parameters of PME model to Fortran 
	static void ModifyFormatVal(double& val,const std::string& format); //!< Modify double according to format to emulate saving/reading to top and crd files

	void PrintMDEneMDOUT(HaVec_double& si_vec,int nstep, double t); //!< Print Energy Components at the point of MD trajectory to MDOUT
	void PrintMinEneMDOUT(HaVec_double& si_vec,int nstep,double rms, double fdmax, int iatmax, std::string labmax);          //!< Print Energy Components at the point of Energy Minimization to MDOUT

	static void TestFFT1(); //!< Test fortran subroutines to perform FFT

	static int OpenAmberMDTrajFortran(int iunit, const std::string& fname, bool write_flag = true, bool formatted = true);
	static int WriteCrdToAmberMDTrajFortran( int iunit, PointContainer* pt_cont, bool save_box = true, bool formatted = true);
	static int CloseAmberMDTrajFortran(int iunit);

public:

	//! \name AMBER input and output files:
//@{
	std::string amber_run_file;  //!< Name of AMBER UNIX run file
	std::string amber_inp_file;  //!< Name of AMBER input file
	std::string amber_top_file;  //!< Name of AMBER topology file
	std::string amber_crd_file;  //!< Name of AMBER initial coordinates file

	std::string amber_out_file;        //!< Name of the Amber out file
	std::string amber_rst_file;        //!< Name of AMBER restart file

	std::string amber_trj_coord_file;   //!< Name of the AMBER MD coordinates trajectory file
	std::string amber_trj_vel_file;     //!< Name of the AMBER MD velocity trajectory file
	std::string amber_trj_ene_file;     //!< Name of the AMBER MD energy trajectory file
	std::string amber_constr_crd_file;  //!< Name of the AMBER coordinates constraints file

	std::string amber_log_file;        //!< Name of the Amber log file

	FILE* trj_coord_fp;   //!< File pointer for AMBER MD coordinates trajectory file
	FILE* trj_vel_fp;     //!< File pointer for AMBER MD velocities trajectory file
	FILE* trj_ene_fp;     //!< File pointer for AMBER MD velocities trajectory file
	
//	std::vector<fpos_t> coord_pt_pos; //!< Starting Positions of MD points in coordinate trajectory file 
//	std::vector<fpos_t> vel_pt_pos;   //!< Starting Positions of MD points in velocity trajectory file 
//	std::vector<fpos_t> ene_pt_pos;   //!< Starting Positions of MD points in energy trajectory file
	
//@}
	std::string title; //!< Title of AMBER job
	
// parameters from /mdin_ctrl_int/

	int ntave; 
	int nsnb; 
	int nrespa; //!< flag to indicate several steps without non-bonded force recalculation
	int ntp; 
    int ntc; 
	int jfastw; 
	int nrespai; 
	double ndfmin;   //!< Number of degrees of freedom to remove when computing temperature from kinetic energy 

	double t; //!< current time of MD trajectory
		
	int mdin;  //!< Fortran unit for mdin (= 35)
	int mdout; //!< Fortran unit for mdout (= 36)
	int mdcrd; //!< Fortran unit for MD trajectory
    int mdvel; //!< Fortran unit for velocities of MD trajectory

// Additional fortran parameters:
	
	int mdinfo_flush_interval; //!< Flush frequency of mdinfo in sec
	int mdout_flush_interval;  //!< Flush frequency of mdout in sec 
	int dbg_atom_redistribution;  //!< flag Don't do atom redivision at every list  build
                                      // (this is really only good for debugging the atom redivision code).
	int loadbal_verbose;       //!< Loadbalancing verbosity.  Default is minimal loadbalancing info in log file.  Values of
                                   //!< 1,2,3... increase amount of info dumped. This variable only has an effect for parallel runs.

	int next_mdout_flush_sec;  //!< time (seconds from reference) when to flush mdout next time 

    HaVec_double atm_crd;  //!< Current Atomic Coordinates in the module
	HaVec_double atm_frc;  //!< Atomic forces 
	HaVec_double atm_vel;       //!< Atom velocities
	HaVec_double atm_last_vel;  //!< Atom velocities on the previous step
	HaVec_double gbl_atm_saved_crd;  //!< Atomic Coordinates from last update of non-bonded list
	HaVec_double gbl_saved_box;  //!< Dimensions of periodic box since last update of non-bonded list

	HaVec_double sys_info; //!< Energy and other information of the system

	int collect_crds;      //!< indicates that atm_crd components needs to be collected from slave nodes     used only in MPI
	int all_crds_valid;    //!< indicates that all components atm_crd are valid (collected from slave nodes) Used only in MPI
	int all_vels_valid;    //!< indicates that all components atm_vel are valid (collected from slave nodes) Used only in MPI

// Periodical Box Parameters:

	HaVec_double pbc_box; //!< Dimensions of the unit cell
	double pbc_alpha;     //!< Angle between periodic unit direction vectors 1 and 2
	double pbc_beta;      //!< Angle between periodic unit direction vectors 1 and 3
	double pbc_gamma;     //!< Angle between periodic unit direction vectors 2 and 3
	int    is_orthog;     //!< Flag to show if periodical box is orthogonal
	int    is_octahedral;   //!< Flag to show if periodical box is octahedral
	HaVec_double recip;     //!< 
	HaVec_double ucell;     //!< 
	HaVec_double cut_factor;  
	HaVec_double reclng;
	double uc_volume;     //!< Volume of the unit cell 
	double uc_sphere;     //!< Largest sphere to fit in unit cell

	HaVec_double tranvec; //!< Translation vectors for unit cell


  int dirfrc_efs;        //!< flag corresponding to DIRFRC in PMEMD - need work to remove DIRFRC_EFS from PMEMD 
  int emulate_ext_amber; //!< Modify ff parameters and atom coordinates to emulate running external pmemd from save .crd and .top file

  enum STATE_INFO { SI_TOT_ENE = 0, 
		              SI_KIN_ENE = 1, 
					  SI_SOLUTE_KIN_ENE = 2,
		              SI_SOLVENT_KIN_ENE = 3,
                      SI_VOLUME = 4, 
					  SI_TOT_PRESS = 5,
                      SI_TOT_EKCMT = 6, // total kinetic energy of centers of masses of molecular fragments 
                      SI_TOT_VIRIAL = 7,
                      SI_POT_ENE    = 8,
                      SI_VDW_ENE    = 9,
                      SI_ELECT_ENE  = 10,
                      SI_HBOND_ENE = 11,
                      SI_BOND_ENE = 12,
                      SI_ANGLE_ENE = 13,
                      SI_DIHEDRAL_ENE = 14,
                      SI_VDW_14_ENE = 15,
                      SI_ELECT_14_ENE = 16,
                      SI_RESTRAINT_ENE = 17,
                      SI_DV_DLAMBDA = 18,
                      SI_DENSITY = 19,
					  SI_PME_ERR_EST = 20,
					  SI_PRESS_0     = 21,
	                  SI_PRESS_1     = 22,
	                  SI_PRESS_2     = 23,
					  SI_EKCMT_0     = 24,  // X - component of kinetic energy of centers of masses of molecular fragments 
	                  SI_EKCMT_1     = 25,  // Y - component of kinetic energy of centers of masses of molecular fragments 
	                  SI_EKCMT_2     = 26,  // Z - component of kinetic energy of centers of masses of molecular fragments 
	                  SI_VIR_0       = 27,
	                  SI_VIR_1       = 28,
	                  SI_VIR_2       = 29,
					  SI_KIN_ENE_PLUS_HALF_DT  = 30,
					  SI_KIN_ENE_MINUS_HALF_DT = 31,
					  SI_KIN_ENE_PBS           = 32,
					  SI_TEMP                  = 33,
					  SI_TEMP_SOLUTE           = 34,
					  SI_TEMP_SOLVENT          = 35,
					  SI_POLAR                 = 36,
					  SI_DIPITER               = 37,
					  SI_DIPRMS                = 38
	                  };

  static const int    SI_CNT = 39;

};


class AmberMMModel
//!< Molecular Mechanics Model in a representation matching AMBER program data structures 
{
public:
	AmberMMModel(MolMechModel* p_mm_model);
	virtual ~AmberMMModel();

	friend class MolMechModel;

	int  UpdateAmberData(); //!< Set AMBER Model data structures from Molecular Mechnics Model in HARLEM representation
	int  InitAmberModelAmoeba(); //!< Set Amoeba Model From Molecular Mechnics Model in HARLEM representation 
	
	int SetAtomNum(int natom_new); //!< Set number of atoms in the model and resize arrays correspondingly 
	
	void FindResMolPartition();  //!< Determine partition of atoms to residues (gbl_res_atoms etc) and molecules (atm_nsp)
    void CalcAddDihParams();     //!< Compute Additional parameters for vectorized Dihedral angles processing gbl_gamc, gbl_gams, gbl_ipn, gbl_fmn
	int SetAtomPosRestrData();  //!< Set data arrays for Harmonic Atomic Position Restraints
	
	void Bcast(MPI_Comm& comm); //!< Broadcast Class Data
	void BcastAtmMass(MPI_Comm& comm); //!< Broadcast atom masses Array
	
	void Clear();  //!< Clear Amber data structures
	void SetUpdateDataFlag(int to_update_flag_new = TRUE ); //!< Set flag to update AMBER data structures from MolMechModel before using them
private:
	int  to_update_amber_data; //!< Flag to update AMBER data structures from MolMechModel before using them

public:
	int ntf;   //!< Omit interaction parameter. To eliminate?
// Parameters from /parmtop_int/

	int natom;   //!< Number of Atoms
	int ntypes;  //!< Number of different VdW parameters 
	int nbonh;   //!< Number of bonds containing hydrogen
	int ntheth;  //!< Number of Valence angles containing hydrogen
	int nphih;   //!< Number of dihedral angles containing hydrogen ?
	int next;    //!< The length of atom non-bonded list
	int nres;    //!< Number of Amber residues (may be different from HARLEM a 1-atom res not supported) 
    int nbona;   //!< Number of bonds not containing hydrogen
	int ntheta;  //!< Number of Valence angles not containing hydrogen
	int nphia;  //!< Number of dihedral angles not containing hydrogen ?
	int numbnd; //!< Number of Valence Bond parameters
	int numang; //!< Number of Valence Angle parameters
	int nptra;  //!< Number of Dihedral Angle parameters 
	int nphb;   //!< Number of H-bonds
	int nspm;   //!< Total number of molecules (solute + solvent)
	int nttyp;   //!< ntypes * (ntypes + 1) / 2; 
    int bonda_idx; //!< first index of bonds not containg hydrogen i
	int anglea_idx; //!< first index of angles not containing hydrogen i
	int diheda_idx; //!< first index of dihedrals not containing hydrogen i
    int gbl_bond_allocsize;   //!< Total of number of bonds allocated 
	int gbl_angle_allocsize;  //!< Total of number of angles allocated
    int gbl_dihed_allocsize;  //!< Total of number of dihedrals allocated
	
	int ibelly;
	
	int igb; 
	int alpb; 
	int rbornstat; 
    int gbsa; 

	double dielc;
	double es_cutoff; 
	double vdw_cutoff; //!< Cutoff for VdW interactions
    double scnb; 
	double scee; 
	double intdiel; 
	double extdiel; 
	double saltcon; 
	double cut_inner; 
	double gb_cutoff; 
    double gb_alpha; 
	double gb_beta; 
	double gb_gamma; 
	double gb_fs_max; 
    double gb_kappa; 
	double gb_neckscale; 
	double arad;
    double bbox_xmin; 
	double bbox_ymin; 
	double bbox_zmin; 
    double bbox_xmax; 
	double bbox_ymax; 
	double bbox_zmax;

	double rgbmax; 
	double offset; 
	double surften; 

// AMOEBA Force field parameters
	int iamoeba; //!< Amoeba force field flag

	int do_amoeba_valence; 
	int do_amoeba_nonbond; 
	int do_bond; 
    int do_ureyb; 
	int do_reg_angle; 
	int do_trig_angle; 
	int do_opbend; 
    int do_torsion; 
	int do_pi_torsion; 
	int do_strbend; 
    int do_torsion_torsion; 
	int do_str_torsion; 
	int do_recip; 
    int do_adjust; 
	int do_direct; 
	int do_self; 
	int do_vdw; 
    int do_induced; 
	int do_vdw_taper; 
	int do_vdw_longrange; 
    int beeman_integrator; 
    int amoeba_verbose;
	 
// Constraints data from / constraints_dat_int / 
	int natc; 
	int belly_atm_cnt;

	int using_pme_potential;   //!< Use PME to calculate non-bonded interactions (calculatable)
	int using_gb_potential;    //!< Use Generalized Born model to computed non-bonded interactions (calculatable)

	int max_res_size; //!< Maximal residue size (not used?) 
	int n_solute_res; //!< Number of solute residues (not used?)
	int n_solute_mol; //!< Number of solute molecules (not used?)

// dynamics_dat parameters

	double tmass;         //!< Total mass of the atoms

	HaVec_int gbl_res_atms; //!< Array of indexes of first atoms of AMBER residues
	HaVec_int atm_iac;      //!< (Fortran based 1) indexes of atoms in the array of atom types
	HaVec_int typ_ico;      //!< Array (ntypes x ntypes) of VdW parameters indexes
	HaVec_int atm_nsp;      //!< the atom submolecule index array
	StrVec atm_igraph;          //!< Atom Names
	StrVec atm_isymbl;          //!< Atom force field symbols
	StrVec atm_itree;           //!< Atom tree symbols
	HaVec_double  atm_charge;   //!< The atom partial charge array for amber pme scaled by diel_const (not used in AMOEBA calculations).
	HaVec_double  atm_mass;     //!< Array of Atomic Masses
	HaVec_double  atm_mass_inv; //!< Arrays of Inverse Atomic Masses
	HaVec_int atm_numex;    //!< Array of total numbers of excluded atoms for atoms  
	HaVec_int gbl_natex;    //!< Excluded Atoms Arrays
	HaVec_double  gbl_cn1;      //!< local copy of array of VdW coefficents (at R^12)
    HaVec_double  gbl_cn2;      //!< local copy of array of VdW coefficents (at R^6 )
	HaVec_double  gbl_asol;     //!< local copy of array of H-bond coefficents (at R^12) 
	HaVec_double  gbl_bsol;     //!< local copy of array of H-bond coefficents (at R^10)
	HaVec_double  atm_gb_radii; //!< local copy of array of Generalized Born Atom radii
	HaVec_double  atm_gb_fs;    //!< local copy of array of Generalized Born Screening parameters

	HaVec_double gbl_rk;    //!< array of Valence bonds force constants 
	HaVec_double gbl_req;   //!< array of Valence bonds equilibrium distances
	HaVec_double gbl_tk;    //!< array of Valence Angles force constants 
	HaVec_double gbl_teq;   //!< array of Valence Angles equilibrium angles
	HaVec_double gbl_pk;    //!< array of dihedral angles force constants 
	HaVec_double gbl_pn;    //!< array of dihedral angles periodicity
	HaVec_double gbl_phase; //!< array of dihedral angles phases
	HaVec_double  gbl_gamc; //!< preprocessed array of dihedral parameters
	HaVec_double  gbl_gams; //!< preprocessed array of dihedral parameters
	HaVec_int gbl_ipn;  //!< preprocessed array of dihedral parameters
	HaVec_double  gbl_fmn;  //!< preprocessed array of dihedral parameters

	HaVec_int  gbl_bond;  //!< Array of Valence Bonds mimicking PMEMD structures
	HaVec_int  gbl_angle; //!< Array of Valence Angles mimicking PMEMD structures
	HaVec_int  gbl_dihed; //!< Array of Dihedral Angles mimicking PMEMD structures

	std::vector<HaVec_double> loc_bond_params; //!< bond parameters used in the module
	std::map<HaVec_double, int, less<HaVec_double> > bnd_par_idx_map; //!< map Bond Parameters to indexes in loc_bond_params

	std::vector<HaVec_double> loc_val_angle_params; //!< used valence Angle parameteres 
	std::map<HaVec_double,int, less<HaVec_double> > vang_par_idx_map; //!< map of indexes val angle parameteres in loc_val_angle_params
 	
	std::vector<HaVec_double> loc_dih_ang_par; //!< used dihedral Angle parameters 
	std::map<HaVec_double, int, less<HaVec_double> > dang_par_idx_map; //!< map of Dihedral Angle to indexes in loc_dih_ang_par

	std::vector<HaResidue*> amber_residues;      //!< List of pointers of HARLEM residues that correspond to AMBER residues 
	std::vector<int>        nat_amber_residues;  //!< Number of atoms in AMBER residues
	StrVec res_labels;  //!< AMBER residue labels

	vector<HaVec_double> loc_point_params;                   //!< Array of Locally used Point Parameters
	map<HaVec_double,int, less<HaVec_double> > ppar_idx_map; //!< Map between Point Params and their position in loc_point_params

	double num_deg;         //!< Total number of degrees of freedom
	double num_deg_solute;  //!< Number of solute degrees 
	double num_deg_solvent; //!< Number of solvent degrees of freedom

	HaVec_double  atm_xc;      //!< References Coordinates of Atoms in harmonic constrains
	HaVec_double  atm_weight;  //!< Weights of atomic constraints
	HaVec_int     atm_jrc;     //!< Array of indexes of constrained atoms (ntc total) 
 
	HaVec_int atm_igroup;   //!< Array of indexes of belly(frozen) atoms  need to write code to set values
	int SetMovingAtomsData(); //!< Set Moving Atoms Array (belly dynamics) from MolMechModel data

	int num_dist_constr;              //!< The number of Atom-Atom Distance constraints
	HaVec_int    dist_constr_idx;     //!< (3 X num_dist_constr) Distance constraints atom indexes and types 
	HaVec_double dist_constr_params;  //!< (2 X num_dist_constr) Distance constraints double parameters
	int SetDistConstrData();             //!< Set Distance Constraints arrays from MolMechModel data 

// AMOEBA parameters:
	HaVec_int atm_poltype;    //!< AMOEBA_ATOM_TYPE_INDEX
	HaVec_int atm_element;    //!< AMOEBA_ATOMIC_NUMBER  (Atom elements)
	StrVec  atm_class_idx;    //!< AMOEBA_ATOM_CLASS_INDEX  

	int n_bond_amoeba;         //!< Number of AMOEBA Valence bonds
	int n_bond_amoeba_params;  //!< Number of AMOEBA Valence bonds parameters 
	HaVec_int  gbl_bond_amoeba;     //!< Array of AMOEBA Valence Bonds ( 2 indexes of atoms and parm_idx ) x N Valence Bonds
	vector<HaVec_double>    bond_amoeba_params;  //!< Parameters of AMOEBA Valence Bonds 
	int bond_amoeba_ftab_degree;         //!< AMOEBA_REGULAR_BOND_FTAB_DEGREE
	HaVec_double bond_amoeba_ftab_coef;  //!< AMOEBA_REGULAR_BOND_FTAB_COEFFS

	int n_urey_bond;             //!< Number of Urey-Bradley bonds
	int n_urey_bond_params;      //!< Number of Urey-Bradley bond parameteres 
	HaVec_int   gbl_bond_urey;     //!< Array of Urey-Bradley bonds   ( 2 indexes of atoms and parm_idx ) x N Urey-Bradley bonds
	vector<HaVec_double> bond_urey_params;  //!< Parameters of Urey-Bradley bonds 
	int bond_urey_ftab_degree;              //!< AMOEBA_UREY_BRADLEY_BOND_FTAB_DEGREE
	HaVec_double bond_urey_ftab_coef;       //!< AMOEBA_UREY_BRADLEY_BOND_FTAB_COEFFS

	int n_angle_amoeba;             //!< The number of regular amoeba valence angles 
	int n_angle_amoeba_params;      //!< The number of regular amoeba valence angles parameters 
	HaVec_int  gbl_angle_amoeba_reg;            //!< Array of regular amoeba valence angles ( 3 indexes of atoms and parm_idx ) x N regular angles
	vector<HaVec_double>  angle_amoeba_params;  //!< Parameters of amoeba valence angles ( regular and trigonal )
	int angle_amoeba_ftab_degree;               //!< AMOEBA_REGULAR_ANGLE_FTAB_DEGREE and AMOEBA_TRIGONAL_ANGLE_FTAB_DEGREE
	HaVec_double angle_amoeba_ftab_coef;        //!< AMOEBA_REGULAR_ANGLE_FTAB_COEFFS and AMOEBA_TRIGONAL_ANGLE_FTAB_DEGREE

	int n_trig_angles;                  //!< The number of trigonal angles
	HaVec_int  gbl_angle_amoeba_trig;   //!< Array of amoeba trigonal angles  ( 4 indexes of atoms and parm_idx ) x N trigonal angles
	
	int n_opbend_angles;                //!< The number of opbend angles
	int n_opbend_angles_params;         //!< The number of opbend angles parameters
	HaVec_int  gbl_opbend_angle;        //!< Array of opbend angles  4 indexes of atoms and parm_idx ) x N opbend angles
	HaVec_double opbend_angle_params;   //!< Parameters of opbend angles

	int n_tors_amoeba;                         //!< The number of Torsion angles in AMOEBA
	int n_tors_amoeba_params;                  //!< The number of Torsion anlgle parameters in AMOEBA 
	HaVec_int  gbl_amoeba_tors_angle;          //!< Array of AMOEBA torsion angles  ( 4 indexes of atoms and parm_idx ) x N torsion angles
	vector<HaVec_double>  tors_amoeba_params;  //!< Parameters of AMOEBA torsion angles

	int n_pi_torsions;                     //!< The number of PI torsion angles
	int n_pi_torsions_params;              //!< The number of PI torsion angle parameters
	HaVec_int  gbl_pi_tors_angle;          //!< Array of PI torsion angles  ( 4 indexes of atoms and parm_idx ) x N PI torsion angles
	vector<HaVec_double>  pi_tors_params;  //!< Parameters of PI torsion angles

	int n_stretch_bend;                    //!< The number Stretch-bend terms
	int n_stretch_bend_params;             //!< The number of Stretch-bend parameters
	HaVec_int  gbl_str_bend_angle;         //!< Array of Stretch-bend angles  ( 3 indexes of atoms and parm_idx ) x N stretch/bend angles
	vector<HaVec_double>  str_bend_params; //!< Parameters of Stretch-bend angles

	int n_tors_tors;                        //!< The number of torsion-torsion terms
	int n_tors_tors_params;                 //!< The number of torsion-torsion parameters
	HaVec_int gbl_tors_tors;                //!< Array of Torsion-torsion terms	
	HaVec_int tors_tors_id_params;          //!< ids of torsion- torsion parameters 
	vector< vector<HaVec_double> >  tors_tors_params; //!< Torsion-torsion parameters
	                                                  // dimension of external vector is the number of torsion-torsion parameters
	                                                  // dimension of internal vector is 6 
	                                                  // 0-idx internal vector is values of torsion 1 in the tables
	                                                  // 1-idx internal vector is values of torsion 2 in the tables
	                                                  // 2-idx - AMOEBA_TORSION_TORSION_TORTOR_TABLE_0X_FUNC
                                                      // 3-idx - AMOEBA_TORSION_TORSION_TORTOR_TABLE_0X_DFUNC_DANGLE1
	                                                  // 4-idx - AMOEBA_TORSION_TORSION_TORTOR_TABLE_0X_DFUNC_DANGLE2
	                                                  // 5-idx - AMOEBA_TORSION_TORSION_TORTOR_TABLE_0X_D2FUNC_DANGLE1_DANGLE2

	HaVec_int atm_amoeba_vdw_type;     //!< Atom Amoeba VdW type
    HaVec_int atm_parent_id;           //!< IDs of parent atoms ( the same atom of non-hydrogen and atom connected to for hydrogen ) 
	HaVec_double atm_parent_weight;    //!< Weight of parent atom in AMOEBA VdW calculations

	double vdw_buffer_delta; //!< VDW_BUFFER_DELTA
	double vdw_buffer_gamma; //!< VDW_BUFFER_GAMMA
	int n_vdw_params;                //!< Number of VdW atom types
	HaVec_double amoeba_vdw_rstars;  //!< (n_vdw_params*n_vdw_params) R(*)(min energy distance) parameters of AMOEBA VdW interaction
	HaVec_double amoeba_vdw_depths;  //!< (n_vdw_params*n_vdw_params) E_min(min VdW energy of interaction of two atom types ) parameters of AMOEBA VdW interaction

	int num_local_multipoles;      //!< The number of local multipoles
	int num_chiral_frames;         //!< The number of chiral frames
	int num_reg_frames;            //!< The number of regular frames
	HaVec_double atm_multipoles;   //!< Atom multipoles array (10 X natom)
	HaVec_int atm_chiral_frames;   //!< AMOEBA_CHIRAL_FRAME_LIST  (3 X N_chiral)
	HaVec_int atm_reg_frames;      //!< AMOEBA_FRAME_DEF_LIST     (5 X N_regular)

	int AddAtomFrames(HaAtom* aptr, AtomFFParam* p_at_ff, StrAtomMap* p_templ_atname_to_res_map ); //!< Add frame vectors for the atom using Template Atom FF parameters 

	int num_adjust_list;             //!< Length of adjust list
	HaVec_int atm_adjust_list;       //!< AMOEBA_ADJUST_LIST 

	HaVec_double adjust_vdw_weights;    //!< AMOEBA_ADJUST_VDW_WEIGHTS_LIST
	HaVec_double adjust_mpole_weights;  //!< AMOEBA_ADJUST_MPOLE_WEIGHTS_LIST
	HaVec_double adjust_direct_weights; //!< AMOEBA_ADJUST_DIRECT_WEIGHTS_LIST
	HaVec_double adjust_polar_weights;  //!< AMOEBA_ADJUST_POLAR_WEIGHTS_LIST
	HaVec_double adjust_mutual_weights; //!< AMOEBA_ADJUST_MUTUAL_WEIGHTS_LIST

	HaVec_double atm_polar;        //!< Atomic polarizabilities (AMOEBA_POLARIZABILITY_LIST)
	HaVec_double atm_hpolar;       //!< Atomic hyper polarizabilities 
	HaVec_double atm_screen_polar; //!< Atomic polarizability screening coef ( sqrt(thole_coef) )
	HaVec_double damp_polar_strength; //!< if > 0.0, atom damp polarization of neighboring atoms 
	HaVec_double damp_polar_sensitivity; //!< if > 0.0, atom polarization is damped by the presence of other atoms
	HaVec_double damp_polar_rad;   //!< atom radii for damping atom polarization
	HaVec_double atm_qterm;      //!< sq_polinv array (atom reverse sqrt of polarizabilities?) for Amoeba.
	HaVec_int    is_polarizable; //!< Flag of indicators that atom is polarizable or not

	MMDriverAmber* p_amber_driver; //!< AMBER Mol Mech Driver class 
private:
	MolSet*     pmset;
	MolMechModel* p_mm_model;      //!< MolMechModel corresponding to the class
	HaMolMechMod* p_mm_mod;        //!< HaMolMech Module associted with the model
};


class TimerAmber
{
public:
	TimerAmber(MMDriverAmber* p_mm_driver_new);
	virtual ~TimerAmber();

	void InitTimers();     //!< Init timing variable, call at the beginning of simulations
	void EndSetupTimers(); //!< Save timings at the end of the setup
	void EndRunTimers();   //!< Save timings at the end of the run
	void PrintTimings();   //!< Print Simulation timings to STDOUT

	void ZeroTime(); //!< Reset timer  ??
	void UpdateTime(int itimer_type); //!< set timer of certain type

	int detailed_timing_flag; //!< flag to use detailed performance timing     

	double run_start_cputime;       
	double run_setup_end_cputime;
	double run_end_cputime;
	double run_start_walltime;
	double run_setup_end_walltime;
	double run_end_walltime;

	// from timer_mod should be synchronized

	static const int FCVE_DIST_TIME =  1;
	static const int NONBOND_TIME   =  2;
	static const int BOND_TIME      =  3;
	static const int ANGLE_TIME     =  4;
	static const int DIHEDRAL_TIME  =  5;
	static const int SHAKE_TIME     =  6;
	static const int RUNMD_TIME     =  7;
	static const int OTHER_TIME     =  8;
	static const int NONSETUP_TIME  =  9;
	static const int MAX_GENERIC_TIMER =  9;

	// PME timer constants; use these in update_pme_time().

	static const int CIT_SETUP_TIMER         = MAX_GENERIC_TIMER + 1;
	static const int BUILD_LIST_TIMER        = MAX_GENERIC_TIMER + 2;
	static const int BSPLINE_TIMER           = MAX_GENERIC_TIMER + 3;
	static const int GRID_CHARGES_TIMER      = MAX_GENERIC_TIMER + 4;
	static const int SCALAR_SUM_TIMER        = MAX_GENERIC_TIMER + 5;
	static const int GRAD_SUM_TIMER          = MAX_GENERIC_TIMER + 6;
	static const int FFT_TIMER               = MAX_GENERIC_TIMER + 7;
	static const int DIR_FRC_SUM_TIMER       = MAX_GENERIC_TIMER + 8;
	static const int ADJUST_MASKED_TIMER     = MAX_GENERIC_TIMER + 9;
	static const int PME_MISC_TIMER          = MAX_GENERIC_TIMER + 10;
	static const int ATM_REASSIGN_TIMER      = MAX_GENERIC_TIMER + 11;
	static const int IMG_REASSIGN_TIMER      = MAX_GENERIC_TIMER + 12;
	static const int FFT_SLAB_REASSIGN_TIMER = MAX_GENERIC_TIMER + 13;
	static const int MAX_PME_TIMER           = MAX_GENERIC_TIMER + 13;

	// GB timer constants; use these in update_gb_time().

	static const int CALC_GB_RAD_TIMER       = MAX_GENERIC_TIMER + 1;
	static const int CALC_GB_DIAG_TIMER      = MAX_GENERIC_TIMER + 2;
	static const int CALC_GB_OFFDIAG_TIMER   = MAX_GENERIC_TIMER + 3;
	static const int DIST_GB_RAD_TIMER       = MAX_GENERIC_TIMER + 4;
	static const int MAX_GB_TIMER            = MAX_GENERIC_TIMER + 4;

	// Maximum timer field must be kept in sync...

	static const int MAX_TIMER               = MAX_PME_TIMER;

	MMDriverAmber* p_mm_driver;
	AmberMMModel*  p_amber_model;

};


#endif // !define MM_DRIVER_AMBER_H

