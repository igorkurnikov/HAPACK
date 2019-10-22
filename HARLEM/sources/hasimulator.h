/*! \file hasimulator.h
 
    General Monte Carlo Simulator Module  

    \author Igor Kurnikov
    \date 1999-2009

*/

#ifndef HASIMULATOR_H
#define HASIMULATOR_H

class TrajPointInfo;
class TrajAnalAgent;
class TraceMolAgent;
class TrajIOAgent;
class TraceMolAgent;
class UpdateMolViewNotifyAgent;
class HaEnergyFunc;
class Random;
class wxThread;

namespace harlem
{
	class Coord;
};

class MCSimulator
//!< Basic MC Simulator
{
public:
	MCSimulator();
	virtual ~MCSimulator();

	long num_mc_steps;        //!< the maximum number of MC iterations, default = 100
	int stop_calc_flag;       //!< it not zero stop docking calculations
	
	harlem::Coord* p_crd; //!< Current Coordinates of the simulated system 

//! \name Control of replay of MC or BD trajectory
//@{
	double delay_time;     //!< delay time between steps in s 
	long npt_begin;  //!< The number of trajectory points to skip in the beginning of the trajectory
	long npt_step;   //!< The number of trajectory points to step between points displayed and processed
	long npt_end;    //!< Index ot the last trajectory point to analyze  
	
	int dont_calc_ene_flag;   //!< if TRUE do not calculate intermolecular energy (just playback coordinates) 
//@}
    wxThread* sim_thread;

	virtual int SetEnergyFunc( HaEnergyFunc* p_ene_func_new); //!< Set Energy Functional
	virtual int InitEnergyFunc(); //!< Initialize if necessary Energy Functional associated with the Simulator 
	virtual double ComputeEnergy(harlem::Coord* pcrd); //!< Compute energy (in kcal/mol) for coordinates defined by pcrd
	virtual int SetInitPoint(harlem::Coord* pcrd_new = NULL);  //!< Set Initial coordinates of the system in simulations

	virtual void SetTemperature(double temp_new); //!< Set Simulation Temperature in K 
	virtual double GetTemperature() const; //!< Get Current Simulation Temperature in K 

	int RunMC();          //!< Run Monte Carlo simulations 
	int RunMCThread();    //!< Run Monte Carlo simulation in a separate thread 
	int PauseMC();        //!< Pause MC simulations 
	int ResumeMC();       //!< Resume MC simulations 
	int StopMC();         //!< Stop MC simulations
	
	virtual int IncrementCrd(harlem::Coord* pcrd) = 0; //!< Make a MC move changing coordinates of the system 
	virtual int SetCoord(harlem::Coord* pcrd) = 0;     //!< Set Coordinates of the system for pcrd

	int InitTrajAnalysis(TrajPointInfo* ppt_info);     //!< Initialize data structure for analysis of MC or BD trajectory
	int ComputePropTrajPoint(TrajPointInfo* ppt_info); //!< Compute properties at the trajectory point
    int FinalizeTrajAnalysis(); //!< Close data structures and files at the end of trajectory run or analysis
//	int SkipTrajPoints(TrajPointInfo* ppt_info,int nskip = 1); //!< Skip nskip trajectory points

	int AnalyzeTrajectory();   //!< Playback and compute properties along trajectory

	vector<TrajAnalAgent*> agents; //!< Array of classes to analyze trajectories

	TrajIOAgent* GetTrajectoryIOAgent(); //!< Set parameters to read/write energy parameters along the trajectory
	int DeleteTrajAnalAgent(TrajAnalAgent* p_ag); //!< Delete Trajectory Analysis Agent from the list

    TraceMolAgent* GetTrajectoryTraceAgent(int create_agent = FALSE);     //!< Get agent for trajectory trace parameters
	UpdateMolViewNotifyAgent* GetMolViewNotifyAgent(int create_agent = FALSE);  //!< Get agent to notify of molecular coordinate changes

	HaEnergyFunc* p_ene_func;

	virtual MolSet* GetMolSet(); //!< Get Associated Molecular Set
	MolSet* pmset;

protected: 
	double kT_MC;           //!< Temperature for MC calculations in (kcal/mol) 
	Random*        p_rand_num_gen;
	
};


#endif // end if !defined(HASIMULATOR_H)
