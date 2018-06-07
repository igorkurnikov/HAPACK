/*! \file trajectoryanalysis.h
 
    Classes for analysis of simulation trajectories

    \author Igor Kurnikov
    \date 2009

*/

#ifndef TRAJANAL_H
#define TRAJANAL_H

class MCSimulator;

namespace harlem
{
	class Coord;
};

class TrajPointInfo
//! Class to provide information about trajectory point to TrajAnalAgent
{
public:
    TrajPointInfo();
	virtual ~TrajPointInfo(); 

	virtual void Clear(); //!< Set content of the Point info to a default state  

	int do_skip;     //!< Flag to supress computations for a given point 
	int ipt;         //!< index of the trajectory point
	harlem::Coord* pcrd;   //!< Current Trajectory point 
	int is_accepted; //!< if point was accepted
	harlem::Coord* pcrd_rejected; //!< rejected point
	double tot_energy;      //!< Energy of current trajectory point
	double energy_rejected; //!< Energy of rejected point
};

class TrajAnalAgent
//! Abstract class for object to Compute properties along Simulation Trajectory 
{
public:
	virtual std::string GetClassName() const { return "TrajAnalAgent"; };  //!< Get Class Name 
	virtual int IsActive() const  { return TRUE; }                         //!< Check if module is activated
	virtual void SetActive(int active_flag) {}  //!< Set active status (TRUE or FALSE)

	virtual int Init( TrajPointInfo* ppt_info = NULL ) { return TRUE; } //!< Initiate analysis of trajectory
	virtual int AnalyzePt(TrajPointInfo* ppt_info = NULL ) = 0;        //!< Compute properties for a trajectory point 
	virtual int Finalize() { return TRUE;}                     //!< Finalize trajectory analysis
};

class TrajIOAgent : public TrajAnalAgent
{
public:
	TrajIOAgent(MCSimulator* p_sim_new);
	virtual ~TrajIOAgent();

	virtual std::string GetClassName() const { return "TrajIOAgent"; }
	virtual int IsActive() const;            //!< Check if module is activated
	virtual void SetActive(int active_flag);  //!< Set active status (TRUE or FALSE)

	virtual int Init(TrajPointInfo* ppt_info);
	virtual int AnalyzePt(TrajPointInfo* ppt_info);
	virtual int Finalize();

	int traj_io_mode;           //!<  from enum TRAJ_IO_MODES
	int output_rejected_points; //!< In MC computations output all generated points (including rejected) into file traj_all_pts_file_name

	std::string traj_file_name;          //!< Trajectory file name 
    std::string traj_all_pts_file_name;  //!< Trajectory all points file name 
	std::string traj_ene_file_name;      //!< Trajectory energy file name

	int save_image_seq_gif;  // flag to save snapshots of trajectory as GIF files
	int save_image_seq_pict; // flag to save snapshots of trajectory as PICT files

	int npt;            //!< Number of analyzed points
	double average_ene; //!< Average energy

	void SetReadCoord  (int set_on = TRUE); //!< Set Read Coordinates of the trajectory on or off  
	void SetWriteCoord (int set_on = TRUE); //!< Set Write Coordinates of the trajectory on or off  
	void SetReadEnergy (int set_on = TRUE); //!< Set Read Energy at the trajectory points on or off 
	void SetWriteEnergy(int set_on = TRUE); //!< Set Write Energy at the trajectory points on or off  

	int IsReadCoord()   const; //!< Check if Coordinates will be read at each call to AnalyzePt()
	int IsWriteCoord()  const; //!< Check if Coordinates will be written at each call to AnalyzePt()
    int IsReadEnergy()  const; //!< Check if Energy will be read at each call to AnalyzePt()
	int IsWriteEnergy() const; //!< Check if Energy will be written at each call to AnalyzePt()
	
	enum TRAJ_IO_MODES {ENERGY_WRITE = 0x00, ENERGY_READ = 0x01, COORD_WRITE = 0x02, COORD_READ = 0x04};

	MCSimulator* p_sim;
private:

	std::fstream ene_file;
	std::fstream traj_file;
	std::fstream all_points_file;

	int active_flag;
};

class TraceMolAgent : public TrajAnalAgent
{
public:
	TraceMolAgent(HaMolSet* pmset_new);
	virtual ~TraceMolAgent();

	virtual std::string GetClassName() const { return "TraceMolAgent"; }
	virtual int IsActive() const;            //!< Check if module is activated
	virtual void SetActive(int active_flag);  //!< Set active status (TRUE or FALSE)

	virtual int Init(TrajPointInfo* ppt_info);
	virtual int AnalyzePt(TrajPointInfo* ppt_info);
	virtual int Finalize();

	HaMolSet* pmset;
	AtomGroup traced_atoms;     //!< Atoms traced during Trajectory replay or simulations

	HaMolecule* trace_mol;
	HaResidue*  trace_res;
	HaChain*    trace_chain;
	int itr_res;

protected:
	int active_flag;
};

class UpdateMolViewNotifyAgent : public TrajAnalAgent
{
public:
	UpdateMolViewNotifyAgent(HaMolSet* pmset_new);
	virtual ~UpdateMolViewNotifyAgent();

	virtual std::string GetClassName() const { return "UpdateMolViewNotifyAgent"; }
	virtual int IsActive() const;            //!< Check if module is activated
	virtual void SetActive(int active_flag);  //!< Set active status (TRUE or FALSE)

	virtual int Init(TrajPointInfo* ppt_info);
	virtual int AnalyzePt(TrajPointInfo* ppt_info);
	virtual int Finalize();

	HaMolSet* pmset;
	int is_moved;
	double update_interval; //!< Minimal interval in seconds to notify molecular views to update molecular image
	unsigned long next_update_time;

protected:
	int active_flag;
};


#endif // !defined(TRAJANAL_H)


