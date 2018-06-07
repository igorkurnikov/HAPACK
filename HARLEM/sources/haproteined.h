/*! \file haproteined.h

    Classes to perform Essential Dynamics and clustering analysis

    \author Igor Kurnikov  
    \date 2011-

*/
#ifndef HA_PROTEIN_ED
#define HA_PROTEIN_ED

#include "hatypes.h"

class AtomGroup;
class MDTrajectory;

class CollectCrdAnalType : public HaEnum1
//! Extended Enum class for Collective Coordinates Analysis Type 
{
public:
	CollectCrdAnalType();
	virtual ~CollectCrdAnalType();

	enum Value { CCRD_SPATIAL_PLATO = 0, CCRD_SPATIAL_RV = 1, 
		         CCRD_TEMPORAL_PCA = 2, CCRD_TEMPORAL_ICA = 3, CCRD_TEMPORAL_SFA = 4, CCRD_TEMPORAL_ISFA = 5 };

	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	CollectCrdAnalType(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
};

namespace swig {
	const CollectCrdAnalType CCRD_SPATIAL_PLATO   = CollectCrdAnalType::CCRD_SPATIAL_PLATO;
	const CollectCrdAnalType CCRD_SPATIAL_RV      = CollectCrdAnalType::CCRD_SPATIAL_RV;
	const CollectCrdAnalType CCRD_TEMPORAL_PCA    = CollectCrdAnalType::CCRD_TEMPORAL_PCA;
	const CollectCrdAnalType CCRD_TEMPORAL_ICA    = CollectCrdAnalType::CCRD_TEMPORAL_ICA;
	const CollectCrdAnalType CCRD_TEMPORAL_SFA    = CollectCrdAnalType::CCRD_TEMPORAL_SFA;
	const CollectCrdAnalType CCRD_TEMPORAL_ISFA   = CollectCrdAnalType::CCRD_TEMPORAL_ISFA;
}

class CollectCrdAnalMod: public HaCompMod
//!  Class to perform Collective Coordinates Analysis of MD trajectory
{
public:
	CollectCrdAnalMod(HaMolSet* new_phost_mset);
	virtual ~CollectCrdAnalMod();
	
	void SetActiveAtomGroup(const std::string& atgrp_name); //!< Set the name of the active atom group
	AtomGroup*   GetActiveAtomGroup();     //!< Get Active Atom Group
	std::string  GetActiveAtomGroupName(); //!< Get Active Atom Group Name  

	int ConvertMDTrajToPlato(); //!< Convert MD trajectory to PLATO format
	int LoadPlatoOutputFile();  //!< Load info from PLATO output file
	int SavePlatoInputFiles();  //!< Save Input files for PLATO program
	int RunPlato(); //!< Run PLATO program
	int CalcTimeProj();     //!< Compute projections of Eigen vectors along trajectory 
	int SaveTimeProjFile(); //!< Save Time Projections to File
	int ShiftAlongEigenVec( int idx_vec, double shift_val ); //!< Shift structure along Eigen Vector idx_vec

	CollectCrdAnalType anal_type; //!< Type of Collective Coordinates Analysis

	int num_clusters;        //!< The number of atom clusters in spatial analysis 
	int num_eig_vec;         //!< The number of computed eigen vectors in temporal analysis or eigen vectors of similarity matrix
	int num_time_proj;       //!< The number of computed time projections in Temporal analysis 
	int sim_matrix_flag;     //!< Flag to compute Similarity matrix
 
	std::string md_traj_fname_plato; //!< The name of MD trajectory file in PLATO format
	std::string active_atgrp;        //!< The name of the atom group to analyze
	std::string plato_input_fname;   //!< PLATO input  file name
	std::string plato_output_fname;  //!< PLATO output file name
	std::string plato_run_fname;     //!< PLATO run file name
	std::string plato_root_dir;      //!< Root directory of PLATO installation
	std::string time_proj_fname;     //!< The name of the file to store time projections of eigen vectors along MD trajectory

	MDTrajectory* p_md_traj;

	enum { LOAD_WVECTORS, LOAD_MIXING_WVECTORS, LOAD_CORRELATIONS} load_vector_type;

	std::vector<double>       eigen_vals;
	std::vector<HaVec_double> eigen_vecs;
	std::vector<HaVec_double> time_projections;
};

#endif


