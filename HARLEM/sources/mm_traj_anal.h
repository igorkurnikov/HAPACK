/*!  \file mm_traj_anal.h

    Classes to Analyze MD trajectories

    \author Igor Kurnikov 
    \date 2010-

*/
#ifndef MM_TRAJ_ANAL_H
#define MM_TRAJ_ANAL_H

#include "trajanal.h"

#ifndef _MSC_VER
class XTCTraj;
#endif
// In development just suggested interfaces:
class MDTrajectory
{
	public:
		MDTrajectory(MolSet* new_pmset);
		virtual ~MDTrajectory();
		
		enum {FormatUnk=0,AMBER_CRD=1,GMX_XTC=2};
		int format; //!< Format of MD trajectory
	
		int Open();
		int Close();
		int ReadNextFrame();
		//int ReadFirstTrajPoint()  { return TRUE; }
		//int ReadNextTrajPoint()   { return TRUE; }
		//int ReadTrajPointNum()    { return TRUE; }
		//int FindNumberTrajPoint() { return TRUE; }
		bool AMBER_CRD_withBoxSize;
		
		int RefreshAllViews();
	public:
		std::string CrdFileName;    //!< Name of the coordinates trajectory file
		std::string VelFileName;    //!< Name of the velocity trajectory file
		std::string EneFileName;    //!< Name of the energy trajectory file
		
		FILE* CrdFile;     //!< File pointer for AMBER MD coordinates trajectory file
		FILE* VelFile;     //!< File pointer for AMBER MD velocities trajectory file
		FILE* EneFile;     //!< File pointer for AMBER MD velocities trajectory file

	protected:
		MolSet* pmset;
		vector<HaAtom*>    Atoms;          //!<  Atoms and other force and mass centers

		HaVec_double cur_crd;
		HaVec_double cur_vel;
		HaVec_double cur_ene;
		
		double CurTime;
		double dt;
		double Time0;
		int frame;
	
		int OpenAMBER_CRD(const char* filename);
		
		int ReadNextFrameAMBER_CRD();
#ifndef _MSC_VER
		XTCTraj *xtctrj;
#endif
		int OpenXTC(const char* filename);
		int ReadNextFrameXTC();	
};

#ifndef _MSC_VER
class XDR;
class XTCTraj
{
	public:
		XTCTraj();
		~XTCTraj();
		
		std::string FileName;
		
		int Open();
		int Close();
		//int ReadFirstFrame();
		int ReadNextFrame();
		
		HaVec_float *r;
		int natoms;
		int step;
		float time;
		float box[3][3];
	protected:
		int XTC_MAGIC;
		int magic;
		FILE *fp;
		XDR  *xdr;
		
		int XTC_CHECK(char *str,bool bResult);
		int xdr_r2f_double(XDR *xdrs,double *r,bool bRead);
		int xdr_r2f_float(XDR *xdrs,float *r,bool bRead);
		//! Check magic number
		int check_xtc_magic();
		//! read header
		int xtc_header();
		int xtc_coord();
		int xdr3dfcoord(float *fp, int *size, float *precision);
		int sizeofint(const int size);
		int sizeofints( const int num_of_ints, unsigned int sizes[]);
		int receivebits(int buf[], int num_of_bits);
		void receiveints(int buf[], const int num_of_ints, int num_of_bits, unsigned int sizes[], int nums[]);
};
#endif

class RMSDAgent;
class AtomCorrAgent;
class MolMechModel;
class MMDriverAmber;
class AmberMMModel;

#ifndef PyObject_HEAD
  struct _object;
  typedef _object PyObject;
#endif

class MDTrajAnalMod
//!< Module to perform analysis of MD trajectory
{
public:
	MDTrajAnalMod(HaMolMechMod* p_mm_mod_new);
	virtual ~MDTrajAnalMod();

	enum MDAnalScriptStages { SCRIPT_CONTINUE = 0, SCRIPT_START = 1, SCRIPT_STOP = 2}; //!< Stages of MD Analysis Script execution
	enum ReadPBoxMode  { PBOX_DO_NOT_READ = 0,  PBOX_READ = 1,  PBOX_READ_IF_SET = 2 };  //!< Modes to read periodical boundary information from MD trajectory
	enum WritePBoxMode { PBOX_DO_NOT_WRITE = 0, PBOX_WRITE = 1, PBOX_WRITE_IF_SET = 2 }; //!< Modes to read periodical boundary information from MD trajectory

	int AnalyzeTrajectory(int sync = TRUE ); //!< playback and analyze MD trajectory ( if(!sync) - asynchronously in a separate thread) 
	int AnalyzeTrajectoryInternal();  //!< Playback and analyze MD trajectory

	int OpenAmberTrajFilesToRead();   //!< Open AMBER trajectory files (crd,vel,ene) to prepare reading trajectory
	int CloseAmberTrajFiles();    //!< Close AMBER trajectory files if open

	int BuildTrajIndex();  //!< Build index of MD trajectory points for fast access
	int ReadTrajPoint();  //!< Read next MD trajectory point 
	int LoadCurrPt();     //!< Load current MD trajectory point (ipt_curr)

	int ConvArbalestTrajToAmber( const std::string& md_traj_arbalest,  const std::string& md_traj_amber);
	int ReduceAmberMDTraj(const char* sh_traj_name, PointContainer* p_sub_group = NULL ); //!< Reduce trajectory to a subset of the molecular set 

	void SetAmberMDCrdTraj(const std::string& md_crd_fname ); //!< Set File name of AMBER MD trajectory (coordinates)
	void SetAmberMDVelTraj(const std::string& md_vel_fname ); //!< Set File name of AMBER MD trajectory (velocities)
	void SetAmberMDEneTraj(const std::string& md_ene_fname ); //!< Set File name of AMBER MD trajectory (energies)

	void SetAlignToFstPt( bool set_par = true );      //!< Set Option to align coordinates to those of the first trajectory point 
	void SetAlignToCurrentCrd( bool set_par = true ); //!< Set Option to align coordinates to current coordinates of the molecular set
	void SetReadPBox( int set_par = TRUE );     //!< Set Option to read periodical boundary box from MD trajectory 
	void SetWritePBox( int set_par = TRUE );    //!< Set Option to write periodical boundary box to reduced MD trajectory
	void SetWrapCrd( int set_par = TRUE ); //!< Set Option to wrap coordinates upon reading MD trajectory point 
	void SetPtBegin(int npt_step);         //!< Set trajectory point to start analysis (1-based, default = 1) 
	int GetPtBegin() const; //!< Get index of trajectory point to start analysis 
	void SetPtStep(int npt_step); //!< Set number of trajectory points to step between analyzed points (default = 1) 
	int GetPtStep() const;       //!< Get number of trajectory points to step between analyzed points
	void SetPtEnd(int npt_step); //!< Set index of last trajectory point to perform analysis ( 1-based )
	int GetPtEnd() const; //!< Get index of trajectory point to end analysis 
	int GetCurrPtIdx() const; //!< Get index of the current point

	double delay_time;   //!< delay time between steps in s 
	int ipt_curr;        //!< Index of the current MD point  ( 1-based )

	std::string traj_script;  //!< Python script to execute at points of trajectory
	std::vector<TrajAnalAgent*> agents; //!< Array of classes to analyze trajectories
	std::vector<PyObject*> py_agents;   //!< Array of python classes derived from TrajAnalAgent to analyze MD trajectory

	int AddAgent( TrajAnalAgent* p_agent);    //!< Add Agent for trajectory analysis
	int DeleteAgent( TrajAnalAgent* p_agent); //!< Delete Agent for trajectory analysis
	int AddPythonAgent( PyObject* p_obj); //!< Add Python Class derived from TrajAnalAgent

	void PrintAgents(); //!< Print to stdout Trajectory Analysis Agents associated with the module 

	TrajAnalAgent*  GetTrajAnalAgent(const char* agent_class_name, int create_flag); //!< Get Trajectory Analysis by class name 
	RMSDAgent*      GetRMSDAgent( int create_flag = FALSE);      //!< Get Atom Superimpose MD analysis agent 
	AtomCorrAgent*  GetAtomCorrAgent( int create_flag = FALSE);  //!< Get Atom Pairwise Correlation agent 
	
//	int DeleteTrajAnalAgent(TrajAnalAgent* p_ag); //!< Delete Trajectory Analysis Agent from the list

	friend class MolMechDlgWX;

protected:
	HaMolMechMod* p_mm_mod;
	MolMechModel* p_mm_model;
	MMDriverAmber* p_amber_driver;
	MolSet* pmset;

	int npt_begin;  //!< The number of trajectory points to skip in the beginning of the trajectory
	int npt_end;    //!< Index ot the last MD trajectory point to analyze    
	int npt_step;   //!< number of trajectory points to step between analyzed points (default = 1) 
	bool align_to_first_pt;    //!< Set to align trajectory to the coordinates of the first point while performing analysis 
	bool align_to_current_crd; //!< Set to align trajectory to the current coordinates of the molecular set while performing analysis 
	int read_pbox;     //!< Read periodical boundary info  0 - do not read, 1 - read, 2 - read if associated molecular set has a periodical box
	int write_pbox;    //!< Write periodical boundary info 0 - do not write, 1 - write , 2 - write if molecular set has a periodical box
	int wrap_crd;      //!< wrap coordinates into periodical box upon reading 
	
	std::vector<fpos_t> pt_pos; //!< Starting position of MD points in MD trajectory file 
};

namespace swig {
	const int SCRIPT_CONTINUE = MDTrajAnalMod::SCRIPT_CONTINUE;
	const int SCRIPT_START    = MDTrajAnalMod::SCRIPT_START;
	const int SCRIPT_STOP     = MDTrajAnalMod::SCRIPT_STOP;
}


//! Class to manage IO along MD trajectory
class MDTrajectoryIOAgent
{
public:

	MDTrajectoryIOAgent(MMDriverAmber* p_mm_driver_new);
	virtual ~MDTrajectoryIOAgent();

	virtual std::string GetClassName() const; //!< Get Class Name 

	virtual int Init();      //!< Initiate analysis of trajectory assume initial point info is supplied  
	virtual int AnalyzePt(int nstep, HaVec_double& si); //!< Compute properties for a trajectory point 
	virtual int Finalize(int nstep);                                         //!< Finalize trajectory analysis

	MMDriverAmber* p_mm_driver;
	AmberMMModel*  p_amber_model;
	HaMolMechMod*  p_mm_mod;
	MolMechModel*  p_mm_model;

	HaVec_double sit;      //!< Hold accumulated energies
	HaVec_double sit2;     //!< Hold accumulated energies squared
	HaVec_double sit_tmp;  //!< Hold accumulated energies for averages  every ntave steps   
	HaVec_double sit2_tmp; //!< Hold accumulated energies squared for averages every ntave steps 

	int nvalid;      //!< Number of steps where energies were evaluated
	int n_saved_crd; //!< Number of coordinates to save to trajectory and restart files

protected:
	
	ofstream constr_trj_fs;
	
};


class RMSDAgent : public TrajAnalAgent
//!< Agent to superimpose atoms on a reference coordinates along MD trajectory and compute relevant statistics
{
public:
	RMSDAgent();
	virtual ~RMSDAgent();

//! \name  Overide virtuals from TrajAnalAgent
//@{
	virtual std::string GetClassName() const { return "RMSDAgent"; }  //!< Get Class Name
	virtual int IsActive() const;           //!< Check if agent is activated
	virtual void SetActive(int active_flag); //!< Set active status (TRUE or FALSE)

	virtual int Init(TrajPointInfo* ppt_info = NULL);  //!< Initiate analysis of trajectory 
	virtual int AnalyzePt(TrajPointInfo* ppt_info = NULL); //!< Compute properties for the trajectory point 
	virtual int Finalize();                         //!< Finalize trajectory analysis
//@}
	int SetAtomsFit(const std::string& at_group_name); //!< Set group name for atoms to fit
	int SetAtomsFit(AtomContainer& active_atoms ); //!< Set atoms to compute RMSD, reference coordinates are set to current fit_atoms coordinates
	
	int SetAtomsRMSD(const std::string& at_group_name); //!< Set group name for atoms to compute RMSD 
	int SetAtomsRMSD(AtomContainer& active_atoms_new); //!< Set atoms to compute RMSD for (active_atoms_new), reference coordinates are set to current active_atoms coordinates

	int SetRefCrdFit(PointContainer& ref_coords, int ref_crd_type = REFC_SPECIAL );  //!< Set reference coordinates for fit
	int SetRefCrdFitFromXYZFile( const std::string& ref_crd_file_name );   //!< Set Reference Coordinates to fit From XYZ File
	int SetRefCrdFitFromAtomGroup( const std::string& at_grp_id, MolSet* pmset_ref); //!< Set Reference Coordinates to fit from Atom Group Name with ID at_grp_id of molecular set pmset_ref_new 

	int SetRefCrdRMSD(PointContainer& ref_coords, int ref_crd_type = REFC_SPECIAL );  //!< Set reference coordinates to compute RMSD
	int SetRefCrdRMSDFromXYZFile( const std::string& ref_crd_file_name );             //!< Set Reference Coordinates to compute RMSD From XYZ File
	int SetRefCrdRMSDFromAtomGroup( const std::string& at_grp_id, MolSet* pmset_ref); //!< Set Reference Coordinates to compute RMSD from Atom Group Name with ID at_grp_id of molecular set pmset_ref_new 

	int SetMolSet(MolSet* pmset); //!< Set Molecular Set for atoms to use in fitting coordinates and RMSD calculations 

	std::string fname_rmsd_out;      //!< File to output RMSD of active_atoms coordinates from reference coordinates 
	std::string fname_rmsd_atom_out; //!< File to output averaged indvidual atom RMSD of active_atoms coordinates from reference coordinates
	std::string fname_rmsf_atom_out; //!< File to output RMSF of active_atoms coordinates from averaged fit coordinates
	std::string avg_crd_file_name;    //!< The name of the output file with averaged coordinates ( xyz format ) 
	
	enum RefCrdType { REFC_CURRENT_CRD = 0, REFC_FIRST_PT, REFC_XYZ_CRD_FILE, REFC_ATOM_ARRAY_ID, REFC_SPECIAL } ref_crd_fit_type; //!< reference coordinates type for atoms to fit  
	RefCrdType ref_crd_rmsd_type;  //!< reference coordinates type for atoms to compute RMSD
	
	int calc_rmsd_per_atom_flag; //!< Flag to compute RMSD per atom 
	int calc_rmsf_per_atom_flag; //!< Flag to compute RMSF per atom 
	int calc_avg_crd_flag;       //!< Flag to calculate average coordinates 


protected:

	FILE* f_rmsd_out;
	int is_initiated_flag;  //!< Flag to indicate that agent was initiated 

	int npt; //!< number of points processed

	AtomGroup fit_atoms;  //!< atoms to fit 
	AtomGroup rmsd_atoms; //!< atoms to compute RMSD

	Vec3DValArray ref_coords_fit;   //!< reference coordinates to fit atoms
	Vec3DValArray ref_coords_rmsd;  //!< reference coordinates to compute RMSD

	MolSet* pmset;      //!< Molecular Set for atoms to fit and atoms to compute RMSD
	MolSet* pmset_ref;  //!< Molecular Set of reference atoms (optional) 

	int active_flag; //!< Active Status flag (TRUE or FALSE)

	HaVec_double rmsd_per_atom; //!< Accumulated rmsd per atom
	HaVec_double rmsf_per_atom; //!< RMSF per atom
	HaVec_double avg_coords;   //!< Accumulated averaged coordinates of active atoms
	HaVec_double sq_coords;    //!< Accumulated squares of coordinates of active atoms

};

class AtomCorrAgent : public TrajAnalAgent
//!< Class to compute pairwise corrleation functions between atoms along MD trajectory 
{
public:
	AtomCorrAgent( MolSet* pmset);
	virtual ~AtomCorrAgent();

//! \name  Overide virtuals from TrajAnalAgent
//@{
	virtual std::string GetClassName() const { return "AtomCorrAgent"; }  //!< Get Class Name
	virtual int IsActive() const;           //!< Check if agent is activated
	virtual void SetActive(int active_flag); //!< Set active status (TRUE or FALSE)

	virtual int Init(TrajPointInfo* ppt_info = NULL);  //!< Initiate analysis of trajectory assume initial point info is supplied  
	virtual int AnalyzePt(TrajPointInfo* ppt_info = NULL); //!< Compute properties for a trajectory point 
	virtual int Finalize();                         //!< Finalize trajectory analysis
	
	int SetAtGroup1ByExpr(const std::string& expr); //! Set First Group of atoms for correlation analysis
	int SetAtGroup2ByExpr(const std::string& expr); //! SetSecond group of atoms for correlation analysis

	void SetDistRange( double rmin, double rmax, int nr = 0 ); //!< Set minimal and maximum interatomic distances to computer g(r) optinally number of points on r ( 0 - choose automatically)

	HaVec_double GetRCut() const; //!< Get vector of cutoff radii used to compute interatomic distances population (nr dimension)

	HaVec_double GetGR()   const { return gr;   }     //!< Get computer G(r) function values. Dimension nr.
	HaVec_double GetNAvg() const { return navg; }     //!< Get averaged number of atoms within cutoff distances. Dimension nr.  
	HaVec_double GetDeltN2() const { return deltn2; } //!< Get fluctuations of number of atoms within cutoff distances. Dimension nr.

protected:

	int npt_proc; //!< number of points processed so far

	AtomGroup     at_grp_1;  //!< Atom Group 1 to use in correlation analysis   
	AtomGroup     at_grp_2;  //!< Atom Group 2 to use in correlation analysis

	double rmin;    //!< minimal interatomic distance
	double rmax;    //!< maximal interatomic distance
	int    nr;      //!< number of cutoff interatomic distances to compute statistics

	HaVec_double ra; //!< Array of interatom distances to compute statistics

	HaVec_int nr_acc;   //!< accumulated sum of number of atoms within radius in ra array
	HaVec_int n2r_acc;  //!< accumulated sum of squared number of atoms within radius value in ra array 

	HaVec_double gr;   //!< G(r) values at values (r[i] + r[i+1])/2   nr dimension. Last point iz zero.
	HaVec_double navg; //!< Computed average number of atoms within radius values in ra array
	HaVec_double deltn2; //!< Computed fluctuation of number of atoms within radius values in ra array

	std::string fname_out; //!< name of the output results file

	MolSet* pmset;

//@}

};











#endif // end of #define MM_TRAJ_ANAL_H