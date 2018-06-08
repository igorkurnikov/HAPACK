/*!  \file haempirical.h

    Classes to calculate scores  
 
    \author Igor Kurnikov and Kirill Speransky
    \date 2005-2006
*/

#ifndef HASCOREFUNCTION_H
#define HASCOREFUNCTION_H

class HaMolecule;
class wxThread;

#include "hastl.h"
#include "hastring.h"
#include "command.h"
#include "hacompmod.h"
#include "halinalg.h"
#include "vec3d.h"
#include "haatom.h"
#include "wx/tokenzr.h"
#include "hasurface.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "hamolmech.h" // jose

class HaEmpiricalMod : public HaCompMod
{
public:
	HaEmpiricalMod(HaMolSet* new_phost_mset);
	~HaEmpiricalMod(); 
	int SetStdParams();
	int Initialize();
	double ScoreEnergy();
	double GeometryScoreEnergy(); // Total score energy from geometry components of function (Repulsion, Box, Angle)
	double PenaltyConstraints();
 	double PenaltyPackDistance();
	double PenaltyPackAngle();
	int CalcPackAngleForceTorque(Vec3DValArray& torque_array); //Caculate Torque due to angel b/n two long axix of the helices
	double PenaltyContact();
	double PenaltyBured();
	double PenaltyPairwise();
	double PenaltyVDW();
	double PenaltyVDW_Bured();
	double PenaltyDensity();
	double PenaltySymmetry();
	double PenaltySolventAccessible();
	double PenaltyHelicePack();
//	double PenaltyLoop();
	int LoadEmpConstrains();
	int LoadSolventAccessibleAtoms();
	int EstablishChains();
	double SoftSqrWellPotential(double& current_value, double& average, double& stdev, double& weight);
	double SqrPotential(double& current_value, double& average, double& stdev, double& weight);
	Vec3DValArray CenterOfMass();
//	Vec3D CenterOfMassHelix(HaMolecule* pMol);
	//HaVec_double FindAxes(HaMolecule* pMol);
	Vec3DValArray FindAxes();
	int LoadEmpParam();
	int CbettaSetUp();
	Vec3D FindCentralAxis();
	int Neighborhood();
    double CheckNeighbor(int i, int j);
	void ResidueTypeList();
	int LineSegments(); // Computes matrix of line segments assosiated with helix
	double  HarmonicEnergy(); // Simple harmonik energy between centars of mass of molecules
	double ToyEnergy(); // Simple toy model 4 points of contacts rigid bodies 
	double MinEnergy(); // Global mimum energy potential
	double PenaltyCentralAttract(); //Keep system near [0, 0, 0] point with (r0- r^2)^2 potential
	int CalcForceCentralAttract(Vec3DValArray& force_cntl_array); // Evaluate the force from (r0- r^2)^2 potential 
	void CalculateCoarseGrainedBackbone(); // Calculate position of Coarse Grained Backbone (4 Calpha combined)
	int CalcRepulForceTorque(Vec3DValArray& force_array, Vec3DValArray& torque_array); // Evaluate the force and torque on each rigid body ((sigma/r)^4 potential)
	double PenaltyRepulsion(); // Evaluate the repulsive energy for each rigid body ( (sigma/r)^4 potential)
	double GetRepulEnergy(HaAtom* aptr1, HaAtom* aptr2); // Evaluate (sigma/r)^4 potential
	double GetRepulDerivative(HaAtom* aptr1, HaAtom* aptr2); // Evaluate Derivative of sigma/r^4 potential
	int GetMaxDimension(); 
	int GetGeomCenter();
	int GetGeomCenterToy(); // function for cylinder model
	int QuantSampling();
	double LJ_ene();
	double LJState_ene();
	double LennardJonesEnergy();
	double BuriedEnergy();
	int InitCylinders();  // functions for cylinder model
	double LJ_eneCylinder(); // functions for cylinder model
	double BuriedEnergyCyl(); // functions for cylinder model
	double PenaltyRepulsionCyl(); // functions for cylinder model
//Help functions form geometry library
	/*
double Segments_Dist_3D ( Vec3D p1, Vec3D p2, Vec3D p3, Vec3D p4);
void Segment_Point_Near_3D ( Vec3D p1, Vec3D p2, Vec3D p,
  Vec3D pn, double *dist, double *t );
double Segment_Point_Dist_3D ( Vec3D p1, Vec3D p2, Vec3D p );
bool Dvec_Eq ( int n, Vec3D a1, Vec3D a2 );
void Dvec_Copy ( int n, Vec3D a1, Vec3D a2 );
double D_Max ( double x, double y );
double D_Min ( double x, double y );
bool Minquad ( double x1, double y1, double x2, double y2, double x3, double y3, 
  double *xmin, double *ymin );
int Parabola_Ex ( double x1, double y1, double x2, double y2, double x3, 
  double y3, double *x, double *y );
  */
double Segments_Dist_3D ( double p1[3], double p2[3], double p3[3],   double p4[3] );
void Segment_Point_Near_3D ( double p1[3], double p2[3], double p[3],  double pn[3], double *dist, double *t);
double Segment_Point_Dist_3D ( double p1[3], double  p2[3], double  p[3]);
bool Dvec_Eq ( int n, double  a1[], double  a2[]  );
void Dvec_Copy ( int n, double  a1[], double  a2[] );
double D_Max ( double x, double y );
double D_Min ( double x, double y );
bool Minquad ( double x1, double y1, double x2, double y2, double x3, double y3, 
  double *xmin, double *ymin );
int Parabola_Ex ( double x1, double y1, double x2, double y2, double x3, 
  double y3, double *x, double *y );

//
	int module_to_init_flag;
    double sigma_constr; // deviation of distance in constrains distances, A
	double pack_dist_com ;// mean distance (center of mass) in packing distance
	double pack_dist_axis ;// mean distance (between axises) in packing distance
	double sigma_pack_dist_com ; // standard deviation of distance (center of mass) in packing distance
	double sigma_pack_dist_axis;// standard deviation of distance (between axises) in packing distance
	double pack_angle;// mean of packing angle
	double sigma_pack_angle;// standard deviation of packing angle
	double weight_constraints; // weight of distance constrains
    double weight_pack_distance; // weight of packing distance
	double weight_pack_angle; // weight of packing angle
	double	dist_contct_com ; // mean distance (center of mass) in helix contact
	double sigma_dist_contct_com ; // standard deviation (center of mass) of helix conact distance

	int num_contact_com ; // number of closest helixes
	double weight_num_contact;
	double weight_vdw; // weight of vdw repulsion
	double weight_bured ; // weight of bured penalty
	double pack_dens; // packing density
	double sigma_pack_dens; // standard deviation of packing density
    double weight_pack_dens;// weight of  packing density
	double sigma_sym ; // deviation of symmetry
	double weight_sym ; // weight of symmetry
	double face_up_bound  ;// solvent accesible face angle
	double weight_sa ;// weight of solvent accesible face
	double dist_neighborhood; // distance for checking neighborhood molecules
	double lenght_factor; // define how one side is longer than the other (in symmetry)
	HaVec_int chain_arr; // array of chain indexes
	int nchain ;  // Number of chains in molecule
	HaVec_int state_old; // states along the trajectory 
	HaVec_int state_new; // new states 
	HaVec_int exch12;    // exchanges along the trajectory
	int curr_state;

protected:
	int npt_dist; // Number of distance constraints 
	int n_sa_atoms ; // Number of solvent accessible atoms 
	HaVec_double emp_dist; // array of distances recorded from distance constrains file
	VecPtr atm_dist; // array of atom pointers for distance constrains
	VecPtr atm2_dist; // array of atom pointers for distance constrains
	VecPtr atm_solacces; // array of atom pointers of solvent accesible residues
	VecPtr atm_sc_array; // array of sidechain atom pointers
	VecPtr atm_ca_array; // array of Calpha atom pointers
	Vec3DValArray atm_cgbb_array; // array of Coarse Grained Backbone (4 Calpha combined) pointers
	HaVec_double rad_cgbb_array; // array of Coarse Grained Backbone (4 Calpha combined) radia
    Vec3DValArray center_arr; // array of helix centers of mass
	Vec3DValArray center_arr_t; // array of cylinders center of mass
//	Vec3D center_helix; // helix center of the mass 
	Vec3D geom_center; // Geometrical center of hamolset
	Vec3D geom_center_t; // Geometrical center of hamolset
    HaVec_double axis_vec; //  Longest axis of molecule
	Vec3DValArray axis_arr;//  Array of helix longest axes
    HaVec_double la_value; //Array of lipid accessibility values (LA scale)
    HaVec_double la_weight_value; //Array of lipid accessibility weights (LA scale)
	wxArrayString residue_arr; //Array of residue names (LA scale)
	HaVec_double sc_vdwradius; // UNRES SC radius
	wxArrayString residue_unres_arr ; //// UNRES residues
	wxArrayString pairwise_name_arr ; //// Pairwise potential, residue raw names
	HaField3D pairwise_energy_arr; //// Pairwise potential, residue energy values
	HaField3D pairwise_energy_arr_sa; //// Pairwise potential, residue energy values for solvent accesible residues
	int com_flag; // flag to perform calculation of function
	Vec3D central_axis ; // central axis of the molset
	double weight_bringtocenter; // weight of bring to center function
	double pack_dens_calculated; // calculated packing density value
	HaMat_double neighbor_mat;  // matrix of 	neighbor  molecules;
	HaVec_int mol_res_correspond; // corresprondence between atom number and resudue type for pairwise scale
	HaVec_int mol_res_correspond_la; // corresprondence between atom number and resudue type for accesibility scale
	HaVec_int marker_res_sa; // corresprondence between atom number and resudue type for accesibility scale
	HaVec_int first_res_mol; // gives the first atom of each chain 
	VecPtr segment_vec; //!< contains the ends of helix segments
	double atomic_volume;
	HaVec_int topology_arr; // aray of molset topology (in or out)
	HaVec_float radii;
	double maximum_dimension; // set max dimension of the simulation box
};

// HaMolMembraneMod Class
// knowledge based potential and future coarse grained models developed for membrane potential 
// ScoreEnergy() function used for Intermolecular interactions, jose
//class HaMolMechMod;  //for composition instead of inheritance
class HaMolMembraneMod : public HaCompMod
{
private:
	HaMolMechMod* MolMechModule;
//protected:
//	HaMolMembraneMod* r;
public:
	HaMolMembraneMod(HaMolSet* new_phost_mset = NULL);
	~HaMolMembraneMod();

	int SetStdParams();
    int module_to_init_flag; //!< Flag to initiate module sturctures at the next calls to module functions
	int module_to_init_HaMolMechMod; //!< Flag to initiate HaMolMechMod structure
	double nonbond_DFIRE_cutoff_dist; //!< Cutoff distance for DFIRE_SCM pairwise potential  @ jose October 22, 2008
	bool pairwiseDfire_flag; //!< Flag of calculations of DFIRE_SCM Pairwise interactions  @ jose October 22, 2008
	int pairwiseDfire_core;//!< Flag of DFIRE_SCM core calculations
	int pairwiseDfire_sa;//!< Flag of DFIRE_SCM solvent accesible calculations
	bool build_nb_coarsegrained_contact_list; //!< Flag to build nonbonded contact list of coarse grained force centroids @ jose November 10, 2008
	int display_results_flag; //!< Flag to show calculation results

	int SetCoarseGrainedDFireCoreParams(); //!< Set CG parameters of AA, DFIRE_SCM CORE @ jose October 22, 2008

	int Initialize();//!< Init module - build MM model and force field 
	vector<HaAtom*> AtomsCentroids; //!<  Atoms and other force and mass centers
	vector<HaResidue*> Residues; //!<  contains residues   jose October 22, 2008
	vector<HaAtom*> LipidInterfaceAtoms; //!<  Atoms located on the interface 
	vector<HaAtom*> CentreAtoms; //!<  Centre Atoms of each molecule
	int ClearMembraneModel(); //!< Clear Membrane Model

	bool BuildNonBondSCContactList(); //!< Build the non-bonded SC contact list for coarse grained calculations jose
	vector<AtomSet> nonbond_SC_contact_list; //!< Non-bonded SC contact list: SC representation atoms of different molecules for which non-bonded interactions are computed according to SCM @ jose October 22, 2008
	bool BuildNonBondCAContactList(); //!< Build the non-bonded Calpha contact list for coarse grained calculations 
	vector<AtomSet> nonbond_CA_contact_list; //!< Non-bonded Calpha contact list: SC representation atoms of different molecules for which non-bonded interactions are computed according to Calpha

	bool BuildClashAtomList(); //!< Build a clash atom list from CA and SC atom force centroids @ jose November 11, 2008
	vector<AtomSet> nonbond_atom_clash_list; //!< Non-bonded atom contact list: atom representation of CA & SC from different molecules for which repulsion forces are computed @ jose November 11, 2008

	int LoadDFireCoreParams();   //!< Initialize DFIRE_SCM CORE parameters @ jose October 22, 2008	
	StrDoubleMap la_value; //!< Map of Adamian Empirical Lipid potential @ jose October 22, 2008
	StrDoubleMap la_weight_value; //!< Map of Adamian Empirical Lipid potential weight @ jose October 22, 2008
	StrDoubleMap sc_vdwradius; //!< Map of SideChain vdW radius @ jose October 22, 2008
    StrVecMap pairwise_energy_vec; //!< Map of residue-residue interactions energy vector to residue-residue key for solvent accesible residues @ jose November 6, 2008
	StrVecMap pairwise_energy_vec_sa; //!< Map of residue-residue interactions energy vector to residue-residue key @ jose November 6, 2008		

	int ScoreEnergy(); //!< Calculates Total Energy of the system
	double pairwise_ene_cg; //!< Energy of pairwise sidechain coarse grainded interactions 
	double vdw_at_repul; //!< Energy of soft alpha carbons and SC atom force centroids van der Waals repulsion
	double vdw_ene_repul; //!< Total energy of soft alpha carbons and SC atom force centroids van der Waals repulsion 
	double lipid_polar_ene; //!< Energy of simple model of heads and core lipid interaction
	double constraint_ene_mol; //!< Energy constraint term for each molecule
	double constraint_ene; //!< Energy term which restraint the molecules to a sphere of radius MAX_radius
	double MAX_radius; //!< Radius of the sphere in which molecules are constrained
	double tot_energy; //!< Total Energy Kcal/mol
	bool CalcVdwRep( HaAtom* pt1,  HaAtom* pt2, double& vdw_at_ene); //!< Energy of soft force centroids van der Waals repulsion
	double PairwiseEnergy( HaAtom* pt1, HaAtom* pt2);//!< Get the pairwise energy between two force centroids
	bool vdw_x4_flag; //!< Flag for repulsion lennard jones 4 functional 
	bool vdw_x12_flag; //!< Flag for repulsion lennard jones 12 functional
	bool vdw_x4_f_flag; //!< Flag for repulsion lennard jones 4 functional soft
	bool vdw_HardSphere; //!< Flag for hard sphere Shakhnovich 2000
	int LoopClosure(); //!< Loop closure constraint
	int EntropySCM(); //<! Side chain entropies

	int OPEP(); //<! Evaluates OPEP Force Field
	int SetCoarseGrainedOPEPParams(); //!< Set MM parameters for coarse-grained model of AA, OPEP ver. 3.0 Maupetit et al 2007
	double SoftSqrWellPotentiala(double& a , double& b, double& c, double& d); //<!
	int CbettaSetUp();
	int LineSegments();
	Vec3DValArray FindAxes();
	double angle_stat;
	double stdev;
	double angle_pack;
	Vec3DValArray axis_arr;//  Array of helix longest axes
	VecPtr segment_vec; //!< contains the ends of helix segments
	VecPtr atm_sc_array; // array of sidechain atom pointers
	VecPtr atm_ca_array; // array of Calpha atom pointers
	wxArrayString residue_unres_arr ; //// UNRES residues
	double anglevar;

protected:

};
#endif // end if !defined(HASCOREFUNCTION_H) 
