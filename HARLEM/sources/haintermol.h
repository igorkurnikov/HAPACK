/*!  \file haintermol.h

    Classes to model intermolecular interactions in HARLEM  
 
    \author Igor Kurnikov 
    \date 1999-2009
*/

#ifndef HAINTERMOL_H
#define HAINTERMOL_H

#include "hastl.h"
#include "hacompmod.h"
#include "hasurface.h"
#include "hasimulator.h"

class HaField3D;

const int NO_ELECTR=0;                 //!< Do not compute electrostatic interactions
const int CONTINUUM_ELECTR = 1;        //!< Compute electrostatic interactions using PB calculations of the complex and its parts
const int COULOMB_ELECTR = 2;          //!< Compute electrostatic interactions using Coloumb law
const int CHARGES_IN_FIELD_ELECTR = 3; //!< Compute electrostatic interactions as the energy of the charges of one molecule in the field of the another

class InterMolMCSimulator;
class InterMolEnergyMinimizer;
class InterMolRepExchSimulator;
class ProtonRedoxMod;

namespace harlem
{
	class RigidBodyCoord;
	class Coord;
};

class Random;
class TraceMolAgent;
class TrajAnalAgent;
class TrajIOAgent;
class UpdateMolViewNotifyAgent;

//! Class to simulate docking between two molecules
class HaInterMolMod : public HaCompMod
{
public:
	HaInterMolMod(MolSet* new_phost_mset);
	virtual ~HaInterMolMod();

	void SetStdParams(); //!< Set Standard Module parameters

//! \name Setup of Interacting Atom Groups 
//@{
	int Initialize();           //!< Attempt to set Interacting groups From Molecules (if not set) and initialize Energy Functional   
	int ClearInternalStruct();  //!< Clear internal structures of the module needed for Intermolecular Energy calculations

	int SetInteractGroupsFromMolecules();        //!< Set Molecules in the Molecular Set as Interacting groups (should be > 1 molecules) 
	int ClearInteractGroups();                   //!< Clear Atom Interaction groups  
	int AddInteractGroup(AtomContainer* p_atgrp); //!< Add Interaction Group

	int module_to_init_flag;      //!< Flag to indicate that module component need to be initialized
//@}
	int SetCoord(harlem::Coord* pcrd);                 //!< Set Coordinates of the system from HaCoord object
	int SetRigidBodyCoord(harlem::RigidBodyCoord* pcrd); //!< Set Rigid Body Coordinates of the system from RigidBodyCoord object

	double CalculateMMEnergy(); //!< Calculate Intermolecular Energy using Molecular Mechanics model
	double CalcElStaticInter(); //!< Calculate Intermolecular Electrostatic Interactions between molecules using Continuum dielectric approach  	                          
    double CalcContElectrEne(std::vector<AtomContainer*> inter_groups); //!< Calculate electrostatic energy beween two groups of atoms	
   
	bool CalcEffInterEne();    //!< Compute effective interaction energy including ET terms
	
	HaMat_double Hessian(int energy_type, VecPtr pMol); // Hessian matrix
	HaVec_double Jacobian(int energy_type, VecPtr pMol); // Jacobian mantrix
	int NormalModes(int energy_type, VecPtr ptmol); // Calculation of normal modes of the system
	
	int to_build_nb_contact_list; //!< flag to build atom nonbond interaction list at the call to functions that use it
	int to_build_intermol_excl_atom_list; //!< flag to build excluded atom list for intermolecular energy calculations

	double cur_intermol_ene;  //!< current energy or effective energy
	double electr_inter_ene;  //!< electrostatic interaction energy
	double vdw_inter_ene;     //!< van-der-Waals interaction energy
	double add_eff_ene;       //!< effective addition to the intermolecular energy from ET rate 
	
    int InitMolecularFields(); //!< Initiate electrostatic and vdW repulsion field around the molecule 
	double CalcChargesInFieldEne(); //!< Calculate intermolecular energy using charges(mol 2) in the field(mol 1) approximation 
	void SetElectrModel(int new_elecr_model_idx); //!< Set electrostatic energy computation method

	int compute_pk;    //!< flag to compute dynamical changes of pKs of ionizable groups in binding
    int electr_model;  //!< method to compute electrostatic energy (NO_ELECTR=0, CONTINUUM_ELECTR = 1, COULOMB_ELECTR = 2, CHARGES_IN_FIELD_ELECTR = 3)
	int calc_et_rate;  //!< Flag to add to the energy the effective addition from ET rate, default = 0 
	int empirical_flag;  //!< run simulations with empirical potential

public:
	InterMolMCSimulator*      p_mc_sim;        //!< MC Simulator
	InterMolEnergyMinimizer*  p_ene_minimizer; //!< Energy minimizer
	InterMolRepExchSimulator* p_rex_sim;       //!< Replica Exchange Simulator
	ProtonRedoxMod*           p_prot_rdx_mod;  //!< Protonation and Redox equilibrium module 

	std::vector<AtomContainer*>  interact_groups; //!< Interacting Atom Groups partitioning Molecular Se 

	HaField3D el_pot_field;   //!< Electrostatic field around the immobile molecule
	HaField3D vdw_pot_field;  //!< Effective VdW 'repulsion' field around the immobile molecule  
    
	HaMat_double normalmode_vec;
	HaVec_double normalmode_val;

private:
};



class InterMolMCSimulator : public MCSimulator
{
public:
	InterMolMCSimulator( HaInterMolMod* p_im_mod_new);
	virtual ~InterMolMCSimulator(); 

	virtual void SetStdParams(); //!< Set Default Simulation Parameters

	bool freeze_first_mol;  //!< flag to freeze first molecule (interacting group) 
	
	double ang_ratio;         //!< scale factor to generate random rotation,  default = 0.02
	double tr_ratio;          //!< scale factor to generate random translation, default = 0.3

	int equil_conf_vol_vdw;   //!< if TRUE run equilibration trajectory in a confined volume with volume constraint
	int amber_flag;      //!< run energy minimization of the system on each step of MC trajectory
	int rex_flag;        //!< turn on replica exchange
	int xy_mc_flag;      //!< flag to run RunMCEmpiricalXY()

    int mc_steps_betw_loc_min;  //!< Number of MC steps between local minimization runs 
	
	double x_orig;
	double y_orig;
	double z_orig;

	HaInterMolMod* GetInterMolMod(); //!< Return Intermolecular Module Associated with the simulator

	virtual int InitEnergyFunc(); //!< Initialize if necessary Energy Functional associated with the Simulator 
	virtual double ComputeEnergy(harlem::Coord* pcrd); //!< Compute energy (in kcal/mol) for coordinates defined by pcrd
	virtual int SetInitPoint(harlem::Coord* pcrd_new = NULL);  //!< Set Initial coordinates of the system in simulations

	virtual int IncrementCrd(harlem::Coord* pcrd); //!< Make a MC move changing coordinates of the system 
	virtual int SetCoord(harlem::Coord* pcrd);     //!< Set Coordinates of the system for pcrd

	int RunMCEmpirical() ;     //!< run MC using empirical scoring function
	int RunMCQuantSampling();  //!< run MC in XYZ using N cylinders and quantitative sampling analysis
	int RunMCEmpiricalXY();    //!< run MC in XY plane using smooth cylinders
    int RunMCEmpiricalNMA() ;  //!<  run MC using empirical scoring function and Normal Modes
	int RunQuasiREM();         //!< Run "quasi" replica exchange on single processor

	int SetDiscretizedMoves(); //!< Set Discretized Moves
	int IsDiscretizedMoves();  //!< Check if MC moves will be discretized 

private: 
	HaInterMolMod* p_im_mod;
	
};

class InterMolEnergyMinimizer
{
public:
	InterMolEnergyMinimizer(HaInterMolMod* p_inter_mol_new) { p_inter_mol = p_inter_mol_new;}
	virtual ~InterMolEnergyMinimizer() {}

	virtual void SetStdParams(); //!< Set Default Simulation Parameters

	int MinimizeEnergy(int energy_type, VecPtr ptmol); //!< Minimize energy 

	int LineSearch(int energy_type, VecPtr ptmol, HaVec_double g, HaVec_double p, double *f, double stpmax);
	int SteepestDescentMinimizer(int nsteps);  //!< Steepest Descent Minimizer
	double GoldenSectionSearch(HaVec_double& xold, Vec3DValArray& force_array, Vec3DValArray& torque_array, double& tol); //!< Line Search
	int StepAlongGradient(HaVec_double& xold, double alpha, Vec3DValArray& force_array, Vec3DValArray& torque_array); //!< Step along gradient for all rigid bodies 

private: 
	HaInterMolMod* p_inter_mol;
};

class InterMolRepExchSimulator
{
public:
	InterMolRepExchSimulator(HaInterMolMod* p_inter_mol_new) { p_inter_mol = p_inter_mol_new;}
	virtual ~InterMolRepExchSimulator() {}

	virtual void SetStdParams(); //!< Set Default Simulation Parameters
	
//	double CalcEnergy( RigidBodyCoord& crd); //!< Calculate Energy for given coordinates

	int nreplicas ; //!< number of replicas
	int rem_steps ; //!< number of REM steps
	double temperature_max ; //!< REM high temperature
	int vary_temperature_flag; //!< turn on temperature variations

	int ireplica ;  //!< current replica
	HaMat_double position_mat; //!< Matrix with coordinates of molecules (used in replica exchange)
	HaVec_double energy_arr;   //!< Array of energies (used in replica exchange)
	HaVec_int exchange_arr;    //!< Array of replica numbers
	std::string MC_traj_file_replica_basename ; //!< File with MC trajectory for replica exchange
	std::string MC_energy_file_replica_basename ; //!< File with MC energies for replica exchange
	std::string MC_rst_file_basename; //!< File with MC coordinates to restart replica exchange
	int n_playback_replica;
	HaVec_double p_acc_ratio; //!< point acceptance ratio

private: 
	HaInterMolMod* p_inter_mol;
};




#endif // end if !defined(HAINTERMOL_H) 
