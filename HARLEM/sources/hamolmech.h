/*!  \file hamolmech.h

    Classes to perform Molecular Mechanics simulations

    \author Igor Kurnikov 
    \date 1999-2002

*/
#ifndef HAMOLMECH_H
#define HAMOLMECH_H

#include "hastl.h"
#include "hastring.h"
#include "command.h"
#include "hacompmod.h"
#include "halinalg.h"
#include "vec3d.h"
#include "haatom.h"
#include "mm_params.h"
#include "mm_force_field.h"

class AtomGroup;
class HaField3D;
class TrajAnalAgent;
class MDTrajectoryIOAgent;
class RMSDAgent;
class MMSysInfo;
class wxThread;     
class MMDriverAmber;
class MMDriverTinker;
class MMDriverGromacs;
class MolMechEvtHandler;
class MDSimMod;
class MinEneMod;
class TISimMod;
class MDTrajAnalMod;
class MolMechModel;
class MolMechDlgWX;

namespace harlem
{
	class RunOptions;
}

const int HA_MOL_MECH_EVENT = 3;

const int MOL_MECH_ID_TEST1 = 3000;
const int MOL_MECH_ID_TEST2 = 3001;
const int MM_DRIVER_AMBER_RUN_INTERNAL = 3002;
const int MM_SET_MPI_COMM_ALL_PROCS = 3003;
const int MM_MOD_SET_MPI_COMM_SPLIT_2 = 3004;
const int MM_INIT_SIMULATIONS_STEP_2 = 3005;
const int MM_MOD_INIT_MIXED_HAMILTONIAN = 3006;
const int MM_UPDATE_CONSTR_2 = 3007;

//! \brief Computational module to perform Molecular Mechanics computations
//!  \nosubgrouping
class HaMolMechMod : public HaCompMod
{
public:

//! \name MM model description and initialization  
//@{
	HaMolMechMod(MolSet* new_pmset);
	virtual ~HaMolMechMod();

	int SetStdParams();
	int Initialize();        //!< Init module - build MM model and prepare MM simulations

	int InitMolMechModel(const ForceFieldType& ff_type = MMForceField::ff_type_default);  //!< Build MM model with force field of a given type
	
	int InitMMSimulations();           //!< Init MM Simulations (Single or Mixed Hamiltonians dep on run_ti flag) 
	int InitSingleHamMMSimulations();  //!< Init MM Simulations with a single hamiltonian inside HARLEM
	int InitMixedHamSimulations(MolMechModel* p_mm_model_2); //!< Init Simulations with mixed hamiltonian (TI and other)
	int InitMixedHamSimulations_node(MolMechModel* p_mm_model_2); //!< Init Simulations with mixed hamiltonian (TI and other) executed on MPI node

	void BcastCtrlParams(MPI_Comm& comm); //!< Broadcast simulation control parameters

	int to_init_simulations;  //!< Flag to initiate module structures at next calls to module functions
	int to_stop_simulations;  //!< Flag to stop current simulations

	int SaveXMLToStream(std::ostream& os, const harlem::SaveOptions* popt = NULL ) const;  //!< Write Molecular Mechanics module data in XML format to stream
//@}

//! \name MM run control
//@{
public:	
	void SetMMRunType( const MMRunType& run_type_new); //!< Set Molecular Mechanics run type
	void SetRunInternal( int run_internal_flag_new = TRUE); //!< Set MM Run Internally or using external program 
	void SetMMExternalProg( const MMExternalProg& ext_mm_prog_new); //!< Set External MM Program to be used to run MM simulations	
	
	int RunCtrlThread(); //!< Run Thread Controlling MM Simulations
	int Run( const harlem::HashMap* popt = NULL ); //!< Run MM simulations if(sync) - synchroneously ( controlled by run_type and run_internal_flag )
	
	int RunMinEne( const harlem::HashMap* popt = NULL ); //!< Run Energy minimization
	int RunMD( const harlem::HashMap* popt = NULL ); //!< Run Molecular Dynamics simulations
	int RunTI(MolMechModel* p_mm_model_2);  //!< Run TI Calculations
	int RunExternal( const harlem::HashMap* popt = NULL ); //!< Run External MM Program ( if popt->ToRunSync() == true - no return before external program exit)
	int RunInternal(); //!< Run Internal MM Calculations of a given type

	int ControlCalc(); //!< Function to control MM Simulations (run by a control thread)

	void RunInternal_node(); //!< function to Run MM Simulations on MPI node

	int CalcEnergy();       //!< Calculate Current MM Energy of the system
	int CalcEnergySimple(); //!< Simple function to Calculate MM energy of the system ( not fully implemented)

	void PrintEneStr(MMSysInfo& info,std::string& str_out); //!< Print energy components to string
	void PrintEneStrAccurate(MMSysInfo& info,std::string& str_out); //!< Print accurate energy components to string
	void PrintLogEne(); //!< Print to log current MM energy and its components

	int StopCalc();      //!< Stop Molecular Mechanics Calculations or trajectory analysis
	int UpdateMolInfo(); //!< Update Current Molecular Coordinates and Energy Info from external simulation program or working internal data structures
	int UpdateMolView(); //!< Update Molecular View with current molecular coordinates and energy info 

	double GetEne() const; //!< Get Potential energy of the system 
	double GetTotEne() const; //!< Get Total Energy of the system (potential + kinetic)  
	double GetPotEne() const; //!< Get Potential Energy of the system 
	double GetConstrEne() const;  //!< Get Contraints energy
	double GetUnConstrEne() const;  //!< Get Potential Energy of the system - energy of constraints

	int LoadAmberRestartFile(const std::string& rst_file_name); //!< Load Restart File in AMBER format 	 

private:
	int run_internal_flag;      //!< if TRUE run MM calculation internally 
	MMRunType run_type;         //!< Molecular Mechanics run type 
	MMExternalProg ext_mm_prog; //!< xternal MM Program to be used to run MM simulations

	static harlem::RunOptions run_opt_default;

	int run_ti;                 //!< Flag to indicate TI run

	long ext_proc_id;  //!< process ID for external MM program 
	wxThread* ctrl_thread; //!< Control Thread for Molecular Mechanics simulations
	wxThread* run_thread;  //!< Thread for Running Molecular Mechanics simulations inside HARLEM
public:
	bool internal_mm_running;  //!< Flag to indicate that internal MM thread is running
	bool ctrl_thread_running;  //!< Flag to indicate that control MM thread is running

//@}

//! \name TI control functions
//@{
public:
	int CheckModelsForTI(MolMechModel* p_mm_model_1, MolMechModel* p_mm_model_2); //!< Check consistency of MM models for TI calculations
	static void CallMMFunctionOnSlaves(int id); //!< Call MM function on Slave Nodes with the given id using remote HA_MOL_MECH_EVENT event
	
	double lambda_ti;             //!< Factor to combine two hamiltonians in TI calculations 

	int SetMPICommSplit2();   //!< Set MPI communicators splitting MPI_COMM_WORLD into two equal sets of processors 
	MPI_Comm inter_model_comm;  //!< MPI communicator between corresponding nodes for two MM models in mixed hamiltonian 
	int inter_model_rank;       //!< rank of the processor in inter_model_comm
	MPI_Comm single_job_comm;   //!< communicator between all nodes of the MM job(with pure or mixed hamiltonian)
	int single_job_rank;        //!< rank of the processor in single_job_comm

private:
	MolMechModel* p_axx_mm_model; //!<  second MM Hamiltonian for mixed hamiltonian calculations
//@}
public:

//! \name MM Model functions
//@{
public:
	MolMechModel* GetMolMechModel(); //!< Get Molecular Mechanics Model associated with the Module
	const MolMechModel* GetMolMechModel() const; //!< Get Molecular Mechanics Model associated with the Module (const version) 

private:
 	MolMechModel* p_mm_model; //!< Molecular Mechanics Model
	virtual int OnDelAtoms(AtomContainer& del_atoms);  //!< Modify module content to react to deleted atoms (from HaCompMod)
	bool Print_info(ostream& sout, const int level);   //!< Print Module info
//@} 

//! \name Event Processing
//@{
	int ProcessEvent(int type, int id); //!< Process Event ( a-la wxEvent) 
//@}

//! \name Simulator and Driver modules performing different MM simulations 
//@{
public:
	MDSimMod*      GetMDSimMod();    //!< Get Module controling MD simulations
	MinEneMod*     GetMinEneMod();   //!< Get Module controling Energy Minimization simulations
	TISimMod*      GetTISimMod();    //!< Get Module controling TI simulations
	MDTrajAnalMod* GetTrajAnalMod(); //!< Get Module to analyze MD trajectories

	MMSysInfo* p_mm_info;  //!< Energy and other parameters of the MM system

	std::string mm_driver_name; //!< Name of the current MM driver

	MMDriverAmber*  p_amber_driver;  //!< Driver class for AMBER  calculations
	MMDriverTinker* p_tinker_driver; //!< Driver class for TINKER calculations
	MMDriverGromacs* p_gromacs_driver; //!< Driver class for GROMACS calculations

private:
	MDSimMod* p_md_mod;    //!< Module controling MD simulations
	MinEneMod* p_min_mod;  //!< Module controling Energy Minimization simulations
    TISimMod* p_ti_mod;    //!< Module controling TI simulations
//@}

//! \name  MD trajectory analysis functions
//@{
	MDTrajAnalMod* p_traj_anal_mod; //!< Module to analyze MD trajectories
	MDTrajectoryIOAgent* p_traj_io_agent;  //!< Agent to Save Info along MD trajectory

	int update_view_flag;     //!< Update Molecular View during MD analysis
	double update_view_interval; //!< minimal time period (in sec) to update molecular view during Molecular Mechanics run or MD trajectory analysis 
//@}
	
//! \name input coord and vel params:
//@{
	MMReadInitCrdType init_read_coord; //!< type of initial coords - corresponds to NTX in Amber 
	int restart_flag; //!< restart MM job - IREST in AMBER
//@}
//! \name output params:
//@{
public:
	void SetPrefix(const std::string& prefix = ""); //!< Set Prefix for output files  
	std::string GetPrefix() const; //!< Get Current Prefix for Output Files

	void SetRestrtFileFormat( const CrdFormatParam& format); //!< Set format of Restart file
	void SetMDtrajFileFormat( const CrdFormatParam& format); //!< Set format for MD trajectory file

	void SetWrtLogFreq(int wrt_freq);    //!< Set frequency (num MD steps) to write to log file 
	void SetWrtRstrtFreq(int wrt_freq);  //!< Set frequency to write restart file
	void SetWrtMDTrajFreq(int wrt_freq, int save_vel = FALSE );  //!< Set frequency to write MD trajectory ( coordinates, energy and velocities(optional))
	void SetWrtCoordFreq(int wrt_freq); //!< Set frequency to write MD coord file
	void SetWrtVelFreq(int wrt_freq);   //!< Set frequency to write MD velocities file 
	void SetWrtEnerFreq(int wrt_freq);  //!< Set frequency to write MD energy file
	void SetWrtConstrFreq(int wrt_freq); //!< Set frequency to write MD constraints info file  

private:
	std::string prefix; //!< Prefix for output files 

	CrdFormatParam write_coord_format; //!< format of final coord and vel output, NTXO in AMBER
	CrdFormatParam traj_wrt_format; //!< format for MD trajectory, IOUTFM in AMBER 

	std::string constr_trj_fname;       //!< Name of the file for trajectory of constraints values

	int wrt_log_freq;        //!< frequency (num MD steps) to write log file 
	int wrt_rstrt_freq;      //!< frequency to write restart file 
	int wrt_coord_freq;      //!< frequency to write MD coord file 	
	int wrt_vel_freq;        //!< frequency to write MD velocities file 
	int wrt_ener_freq;       //!< frequency to write MD energy file 
	int wrt_constr_freq;     //!< frequency to write MD constraints info file

	int limit_wrt_atoms;      //!< NTWPRT of AMBER, = 0 all atoms are written, < 0 only the solute
	                          //!< > 0 - first limit_wrt_atoms are written
	int wrap_coord;           //!< wrap coordinates to the main unit cell - corresponds to amber IWRAP
//@}

//! \name Energy minimization parameters:
//@{
public:

	void SetEneMinMethod( const EneMinMethod& method ); //!< Set type of energy minimization
	void SetMaxNumMinimSteps( int max_num_minim_steps_new);    //!< Set the maximal number of energy minimization steps
	void SetNumSteepDescentSteps( int nsteps );         //!< Set the number of steepest descent minimization steps 
	void SetInitMinStep( double init_min_step_new );    //!< Set initial size for minimization step 
	void SetGradCnvrgVal( double grad_cnvrg_val_new );   //!<  Set convergence criterium for energy gradient in energy minimization
	void SetZMatMin( bool set_par = true ); //!< Set energy minimization using Z-matrix
	bool IsZMatMin() const; //!< Check if energy minimization will Z-matrix
	
private:
	EneMinMethod min_type; //!< type of energy minimization 
	bool zmat_min; //!< flag to indicate energy minimization using Z-matrix  

    int max_num_minim_steps;         //!< Maximal number of energy minimization steps
	int num_steep_descent_steps; //!< The number of steepest descent steps 

	double init_min_step;   //!< Initial minimization step 
	double grad_cnvrg_val;  //!< Convergence criterium for energy gradient in energy minimization
//@}

//! \name MD parameters:
//@{
public:
	void SetNumMDSteps( int num_md_steps_new ); //!< Set the number of MD steps
	void SetRemoveInitRBMotion( int remove_init_motion = TRUE ); //!< Set flag to remove initial translational and rotational motion
	void SetRemoveRBMotionFreq( int freq ); //!< Set frequency to remove rigid body translational and rotational motion of the system

	void SetStartVelMethod( const StartVelMethod& start_vel_method_new); //!< Set method to generate start velocity
	void SetStartTime( double start_time_new );    //!< Set start time(ps) of MD trajectory
	void SetMDTimeStep( double md_time_step_new ); //!< Set MD time step(ps)
	void SetNBListUpdateFreq(int freq); //!< Set Frequency (md steps) to update non-bonded atoms list
	void SetPerBoundaryCondType( const PerBoundaryCondType& type); //!<  Set periodical boundary conditions type

private:
	int num_md_steps;               //!< The number of MD steps (The length of the MD run) (number of MD steps)
    int remove_init_rb_motion_flag; //!< Flag to remove initial rigid body translational and rotational motion flag
	int remove_rb_motion_freq;      //!< Frequency to remove rigid body translational and rotational motion of the system

	StartVelMethod start_vel_method; //!< Method to generate start velocity 
	double start_time;      //!< Start time(ps) of MD trajectory 
	double md_time_step;    //!< MD time step(ps) 
	int nb_list_update_freq; //!< Frequency to update non-bonded atoms list
	PerBoundaryCondType period_bcond;  //!< periodical boundary conditions type
//@}

//! \name Temperature regulation
//@{
public:
	void SetRefTemp( double ref_temp_new );   //!< Set Reference temperature to keep system at 
	void SetInitTemp( double init_temp_new ); //!< Set Initial temperature of the system
	void SetLangevinDumpConst( double langevin_dump_const_new ); //!< Set Langevin dynamics dumping constant (in ps^-1)
	void SetRandomSeed( int random_seed_new ); //!< Set a seed for random number generator
	void SetScaleInitVel( double scale_init_vel_new ); //!< Set a factor to scale initial velocities
	void SetTempCtrlMethod( const TempCtrlMethod& method); //!< Set temperature control method
	void SetTempDeviation( double temp_deviation_new ); //!< Set minimal temperature deviations trigering adjustment temperature algorithm
	void SetTempRelaxTimeSolute( double temp_relax_time_solute_new );   //!< Set Temperature relaxation time for solute
	void SetTempRelaxTimeSolvent( double temp_relax_time_solvent_new ); //!< Set Temperature relaxation time for solvent
	void SetVelLimit( double vel_limit_new ); //!< Set limiting value for atom velocities

private:
	double ref_temp;       //!< Reference temperature to keep system 
	double init_temp;      //!< Initial temperature of the system
	double langevin_dump_const; //!< Langevin dynamics dumping constant (in ps^-1) 
	int random_seed;       //!< seed for random number generator
	double scale_init_vel; //!< a factor to scale initial velocities
	TempCtrlMethod temp_control_method; //!< Switch for scaling temperature method (NTT in AMBER)
	int rand_vel_freq; //!< Frequency for velocity randomization 
	                       
	int last_solute_atom;           //!< last solute atom. ISOLVP in AMBER
	double temp_deviation;          //!< temperature deviations trigering adjustment temperature algorithm
	double temp_relax_time_solute;  //!< Temperature relaxation time for solute
	double temp_relax_time_solvent; //!< Temperature relaxation time for solvent
	double vel_limit;               //!< limiting value for atom velocities
//@}

//! \name Pressure regulation:
//@{
public:
	void SetPressureRegMethod( const PressureRegMethod& method); //!< Set pressure regulation method 
	void SetRefPressure( double ref_pressure_new );    //!< Set Reference Pressure 
	void SetCompressibility( double compressibility_new ); //!< Set Compressibility to use in pressure regulation algorithm 
	void SetPressRelaxTime( double press_relax_time_new ); //!< Set Pressure relaxation time (in ps)

private:
	PressureRegMethod pressure_reg_method; //!< Pressure regulation method 
	double ref_pressure;     //!< Reference pressure 
	double compressibility;  //!< compressibility in 1.0E-06/bar(water)
	double press_relax_time; //!< Pressure relaxation time (in ps)
//@}

//! \name SHAKE bond length constraints:
//@{
public:
	void SetShakeConstr( const MMShakeParam& shake_constr_new); //!< Set Shake bond constraint option
	void SetShakeTol( double shake_tol_new );  //!< Set geometrical tolerance for coordinates resetting in SHAKE algorithm

private:
	MMShakeParam shake_constr; //!< Shake bond constraint option
	double shake_tol; //!< relative geometrical tolerance for coordinate resetting in SHAKE algorithm
//@}

//! \name Special Water treatment:
//@{
private:
	int solute_solvent_image_flag; //!< flag does solute sees 
	                               //!< solvent images, IMGSLT in AMBER
                                   //!< if = 0, IMGSLT = 1
	int remove_solute_nb_cut_flag; //!< flag to remove non-bonded cutoff from the solute
	                               //!< IFTRES in AMBER, if =1, IFTRES = 0   

	enum{ NORMAL_FAST_WATER = 0,REDEF_WATER_NAMES = 1,ONLY_SHAKE_FAST_WATER = 2,
		  ONLY_SHAKE_FAST_WATER_REDEF_NAMES = 3, NO_FAST_WATER = 4} fast_water_method; //!< Fast Water definition flag for special treatment 
	                                                                                   //!< of TIP3P waters, JFASTW in AMBER
//@}

//! \name Interactions with graphical interface
//@{
public: 

#if !defined(HA_NOGUI)
	MolMechDlgWX* p_mm_dlg; 
#endif
	void OnChangePeriodicity(); //!< Function to execute on creating/deleting periodical box

//@}

	friend class MolMechModel;
	friend class MMDriverAmber;
	friend class MMDriverGromacs;
	friend class AmberMMModel;
	friend class MMRunInternalThread;
	friend class MMCtrlThread;
	friend class AnalyzeMDTrajThread;
	friend class MolMechDlgWX;
	friend class MDTrajectoryIOAgent;
	friend class MDTrajAnalMod;
	friend class HaInterMolMod;
	friend class MolMechEvtHandler;
	friend class HaMolMembraneMod;
	friend class PerBoundaryCondType;
	friend class TimerAmber;

public: // Tests:
	static void TestSaveAmoebaTopFile1(); //!< Test Saving Amoeba Top File ( Setup of MORT model and save using HARLEM functions) 
	static void TestSaveAmoebaTopFile2(); //!< Test Saving Amoeba Top File ( Setup using MORT model and save using MORT functions)

};

class MMSysInfo
{
public:
	MMSysInfo(HaMolMechMod* p_mm_mod_new);
	virtual ~MMSysInfo();

	void clear();

	int    nstep;  //!< Current num step of Minimization of MD trajectory
	double time;       //!< Current trajectory time
	double temp;       //!< Current simulation temperature
	double temp_solute; //!< Solute temperature
	double temp_solv;   //!< Solvent temperature
	double press;      //!< Current pressure in MD
	double pres_x;     //!< Pressure along X axis 
	double pres_y;     //!< Pressure along Y axis
	double pres_z;     //!< Pressure along Z axis
	double press_scale_solute;
	double press_scale_solvent;

	double tot_energy;    //!< Current computed MM energy (kcal/mol)

	double kin_ene;    //!< Current kinetic energy
	double pot_ene;    //!< Current Potential energy
	double bond_ene;   //!< Current Valence Bond energy
	double vang_ene;   //!< Valence Angle Energy
	double dihed_ene;  //!< Current Dihedral Angle energy
	
	double vdw_ene;       //!< Current VdW energy 
	double vdw_ene_14;    //!< Current VdW 1-4  energy 
	double vdw_ene_nb;    //!< Current VdW not 1-4  energy 

	double electr_ene;    //!< Current Electrostatic energy 
	double electr_ene_14; //!< Current Electrostatic 1-4 energy 
	double electr_ene_nb; //!< Current Electrostatic non 1-4 energy 
	double polar_ene;     //!< Polarization energy
	double polar_dip_iter; //!< Interaction energy of induced dipoles?          
	double polar_dip_rms;  //!< ??				   

	double gb_ene;        //!< Generalized Born Energy
	double hbond_ene;     //!< Energy of H-bonds

	double constraints_ene; //!< Energy of constraints

	double epol;           //!< Electronic polarization energy
	double e3body;         //!< 3-Body terms enregy

	double kin_ene_plus_half_dt;   //!< Kinetic energy corresponding to atomic velocities at point t + 1/2 dt
	double kin_ene_minus_half_dt;  //!< Kinetic energy corresponding to atomic velocities at point t - 1/2 dt
	double kin_ene_pbs;  //!< Kinetic energy corresponding to cross product of atomic velocities at point t - 1/2dt and t + 1/2dt

	double kin_ene_com;   //!< kinetic energy of centers of mass of molecular fragments ( for pressure calculations)
	double kin_ene_com_x; //!< X component of kinetic energy of centers of mass of molecular fragments  
	double kin_ene_com_y; //!< X component of kinetic energy of  centers of mass of molecular fragments 
	double kin_ene_com_z; //!< X component of kinetic energy of  centers of mass of molecular fragments 

	double kin_ene_solute;  //!< Kinetic energy of the solute
	double kin_ene_solvent; //!< Kinetic energy of the solvent 

	double virial_tot;  //!< Virial of the system
	double virial_x;    //!< Virial of the system along X axis
	double virial_y;    //!< Virial of the system along Y axis
	double virial_z;    //!< Virial of the system along Z axis

	double volume;         //!< Volume of the system
	double density;        //!< Average density in g/mL

	double dv_dlambda;  //!< Current dV/dl in TI calculations

	double av_perm_moment; //!< Permanent dipole moment ?
	double av_ind_moment;  //!< Induced dipole moment ?
	double av_tot_moment;  //!< Total dipole moment ?
	
	double pme_err_est;    //!< estimation of PME error 
	double rms_ene;        //!< RMS of energy during optimization
	double grad_ene_max;   //!< Maximum gradient of energy during optimization

private:
	HaMolMechMod* p_mm_mod;
};

class MDSimMod
//! Class to perform MD simulations
{
public:
	MDSimMod(HaMolMechMod* p_mm_mod_new);
	virtual ~MDSimMod();

protected:
	HaMolMechMod* p_mm_mod;
};

class MinEneMod
//! Class to perform Energy minimization simulations
{
public:
	MinEneMod(HaMolMechMod* p_mm_mod_new);
	virtual ~MinEneMod();

protected:
	HaMolMechMod* p_mm_mod;
};

class TISimMod
//!< Call to perform Thermodynamic Integration Calculations of free energy changes 
{
public:
	TISimMod(HaMolMechMod* p_mm_mod_new);
	virtual ~TISimMod();

	friend class MMDriverAmber;
	friend class HaMolMechMod;
private:
	int ti_sync_freq;           //!< Frequency to force synchronization of coordinates and velocities for 2 hamiltonians in TI calculations          
	double klambda_ti;          //!< power factor in mixing H(lmb) = H_1*(1-lmb)**k + H_0*(1 - (1-lmb)**k)  

	int num_lmb; //!< Number of lambdas to compute TI
	int cur_idx_lmb; //!< Current labmda index (from 0 to num_lmb-1)
	int max_idx;     //!< maximum lambda index (if < 0 all lambda indexes to num_lmb will be computed )  
	enum {START_CALC=0,CONTINUE_CALC} calc_mode; //!< mode of simulation start over or continue previous if avalable restart and dvdl_files
	int num_equilib_points; //!< Number of points to skip 

public:
	void SetNumLambda(int num_lmb_new); //!< Set Number of lambdas in TI calculations (and standard output file names)
	int  GetNumLambda();  //!< Return number of lambdas in TI caclulations
	void SetCurIdxLambda(int cur_idx_lmb_new); //!< Set Current labmda index (from 0 to num_lmb-1)
	void SetMaxLambdaIdx(int max_idx_new); //!< Set maximum lambda index (if < 0 all lambda indexes to num_lmb will be computed )  
	void SetTI_OutputFileNames(); //!< Set output file names corresponding to a given lambda index
	double GetLambdaByIdx(int idx); //!< Get value of lambda corresponding to current num_lmb and idx (0-bazed)
	double GetCurLambda(); //!< Get Current Lambda value
	double GetIntegWtByIdx(int idx); //!< Get value of weight in Numerical Integration corresponding to current num_lmb and idx (0-bazed)
	void SetNumEqPoints(int num_eq_pt_new); //!< Set the number of points to skip when computing <dV/dL>
	
	void CollectForceAndEneTI(HaVec_double& si); //!< Collect Forces and energies for TI calculations
	void SincCrdAndVelTI(); //!< Sync Coordinates and Velocities for 2 hamiltonians of TI calculations

	int ComputeDvDlAvg(); //!< Compute DV/DL values averages for different lambdas from accumulated values in dv/dl output files 
	double CalcDeltaG(int recalc_dvdl_avg = TRUE);  //!< Calculate Delta G using free computed 

	void ReduceDvDlData(int n_avg, const char* file_name = "dvdl_red.dat"); //!< reduce by averaging dvdl data to a file

	std::string GetFilePrefixIdx(int idx); //!< Get File Prefix for index 
	std::string GetCurFilePrefix(); //!< Get File Prefix corresponding to current cur_idx_lmb
	
	std::string file_prefix; //!< prefix for output files ( as prefix_lmb_2_5_dvdl.dat for dv/dl values files, prefix_lmb_2_5.mdcrd for MD trajectory etc....)
	FILE* file_dvdl;            //!< File to save DV/DL values
	HaMolMechMod* p_mm_mod;

	static HaVec_int allowed_num_lmb; //!< Allowed values of number of lambdas for integration 

	HaVec_double dvdl_avg; //!< Array of averages of dv/dl for lambda indexes from 1 to num_lmb
	double delta_g;        //!< Computed Delta G value
};

class MMDriver
{
public:
	virtual std::string GetClassName() = 0;

	virtual int CalcEnergy() = 0; //!< Calculate energy of the system and save results to p_mm_info member of p_mm_mod 
	
	HaMolMechMod* p_mm_mod;
	MolMechModel* p_mm_model;
	MolSet*     pmset;

	int to_save_input_files; //!< Flag to indicate that input files are needed to be saved before running external program
};

#endif // end if !defined(HAMOLMECH_H) 

