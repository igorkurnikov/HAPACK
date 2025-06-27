/*! \file protonredox.h

    Classes to study protonation and redox equilibrium

    \author Igor Kurnikov  
    \date 2010-

*/
#if !defined(PROTONREDOX_H)
#define PROTONREDOX_H

#include "hatypes.h"
#include "hacompmod.h" 
#include "haatgroup.h"


class AltChemStateType : public HaEnum1
//! Enum class for Type of Alternative Chemical State 
{
public:
	AltChemStateType();
	virtual ~AltChemStateType();

	enum Value { PROTONATED = 0, UNPROTONATED = 1, REDUCED = 2, OXIDIZED = 3 };

	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	AltChemStateType(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
};

namespace swig {
	const int PROTONATED    = AltChemStateType::PROTONATED;
	const int UNPROTONATED  = AltChemStateType::UNPROTONATED;
	const int REDUCED       = AltChemStateType::REDUCED;
	const int OXIDIZED      = AltChemStateType::OXIDIZED;
}

class AltChemState
//! Class to define an alternative protonation (or redox) state of the AtomGroup
{
public:
	AltChemState();
    AltChemState(AtomGroup* new_host_atom_set);
	AltChemState(const AltChemState& ref_state);
	virtual ~AltChemState();

	AtomGroup* GetHostAtomGroup() { return host_atom_group; }
	void SetStdParam();       //!< Set Standard parameters

	int SetAltCharges(double weight); //!< Set Charges on the AtomGroup with a weight of the alternative chemical state 
	int SetAltChForAtom(const char* at_name, double new_ch); //!< Set Alternative charge for atom with name

	std::string   id;            //!< Identificator of the Alternative Chemical State
	PtrDoubleMap chmap;       //!< map of changes of atom charges during the transition to the Alternative Chemical State
	std::string mod_atom_name;   //!< Atom name to be protonated
	AltChemStateType  alt_state_type;  //!< type of alternative chemical state 
    double   pk;              //!< pKa value for a transition 
	double   std_pk;          //!< pKa at standard conditions
	int      active_flag;     //!< activity flag for the transition 
	
	AtomGroup* host_atom_group;  //!< AtomGroup the alternative chemical state is defined for
};

#if defined(SWIG) 
%template(vector_AltChemState)  std::vector<AltChemState*>; 
#endif

class MultiSitePopulationMethod : public HaEnum1
//! Enum class for Method to compute population of multiple titration sites 
{
public:
	MultiSitePopulationMethod();
	virtual ~MultiSitePopulationMethod(); 

	enum Value { SCF_MULTI_SITE_CALC = 0, MC_MULTI_SITE_CALC = 1, PFUNC_MULTI_SITE_CALC = 2 };

	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	MultiSitePopulationMethod(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
};

namespace swig {
	const int SCF_MULTI_SITE_CALC    = MultiSitePopulationMethod::SCF_MULTI_SITE_CALC;
	const int MC_MULTI_SITE_CALC     = MultiSitePopulationMethod::MC_MULTI_SITE_CALC;
	const int PFUNC_MULTI_SITE_CALC  = MultiSitePopulationMethod::PFUNC_MULTI_SITE_CALC;
}


class ProtonRedoxMod: public HaCompMod
//!  Class to compute properties of chemical tranformation in the system such as protonation and redox changes
{
public:
	ProtonRedoxMod(MolSet* new_phost_mset);
	virtual ~ProtonRedoxMod();
 
	std::vector<AltChemState*> alt_chem_states; //!< Alternative chemical states of the system 

	void ClearResAltChemStates(HaResidue* pres); //!< Clear Alternative Chemical States for Residue
	void ClearAltChemStates(); //!< Clear all Alternative chemical states of the system

	int GetNumResAltChemStates(HaResidue* pres);   //!< Get the number of alternative chemical states for residue
	AltChemState* GetResAltChemState(HaResidue* pres, int alt_state_idx); //!< Get alternate chemical state of residue by index
	AltChemState* GetResAltChemStateByAtName(HaResidue* pres, const char* at_name); //!< Get alternate chemical state of residue by protonated atom name

	static int heme_model; //!< Model for heme redox center =0 - extra charge distributed, =1 located on Fe
    void PrintResPKa(HaResidue* pres); //!< Print current pKa values for alternative protonation states
	
	AltChemState* AddAltChemState(AtomGroup* pgrp); //!< Add Empty Alternative Chemical State for Residue pres

	int SetStdResPKa(HaResidue* pres, int set_redox_pot); //!< Set standard pKa and E_1/2 values for the residue
	int SetStdResPKa_G1(HaResidue* pres, int set_redox_pot);  //!< Set standard pKa and E_1/2 values for most protonation active residues (HIS,HEM,terminal groups) 

	void SetStdPKa(); //!< Set Standard pKa values for protonation and deprotonation of residues 
	void SetStdPKa_G1(); //!< Set Standard pKa values for most activily protonatable residues(HIS,HEM, terminal groups, etc) 
	void CalcPKaForSelection(); //!< Calc PKa for residues in the selected part using interaction matricies
	void CalcPKaForSelection(bool pnp);

	bool SetStdPKforAtName(HaResidue* pres, const char* at_name, double std_pk_new); //!< Set standard pk value for the protonation involving atom at_name
	
	static int CalcAvgPopMC (HaMat_double& inter_mat,  HaVec_double& avg_st_pop, HaVec_double& alt_st_ene,int N_mc_cyc=10000);  //!< Calculate average alt state population using Monte-Carlo simulations
	static int CalcAvgPopSCF(HaMat_double& inter_mat,  HaVec_double& avg_st_pop, HaVec_double& alt_st_ene); //!< Calculate self-consistenly average population of the alt states
	static int CalcAvgPopPFunc(HaMat_double& inter_mat,  HaVec_double& avg_st_pop, HaVec_double& alt_st_ene); //!< Direct calculation of average populations of alternative states using partition function

	static int TestCalcPopFun(); //!< Test functions for calculations of average state polulations

	void SetAltStatesActive(int set_flag); //!< Set Alternative Residue states active (or inactive) for selected residues
    
	bool SetResChargesForPH( HaResidue* pres, double pH_val); //!< Set Residue Atomic Charges corresponding to a given pH value and current pKa of its alternative chemical state
	bool SetChargesForPH( double pH_val); //!< Set Atomic charges for a certain pH
	bool SetChargesForCurrentPH();        //!< Set Atomic charged corresponding to a current pH
	void SetPH(double new_ph); //!< Set current pH value in the system
	double GetPH() const;      //!< Get current pH value in the system

	int set_std_redox_pot;    //!< flag to set standard redox-potentials for redox groups
	int save_alt_st_inter;    //!< flag to save alternative states interaction matrix to "ALT_ST_INTER.DAT"
	int read_alt_st_inter;    //!< flag to read alternative states interaction matrix from "ALT_ST_INTER.DAT"
    int save_titration_data;      //!< model redox and pH titration
	int save_only_redox_titr; //!< save only redox potential titration
	
    int n_mc_cyc;      //!< number of Monte-Carlo cycles when computing pKa and redox-potentials
	MultiSitePopulationMethod multi_site_pop_method; //!< method to compute population of multiple protonation and redox sites 
	
	double ph_min;  //!< minimal value of pH in titration simulation
	double ph_max;  //!< maximal value of pH in titration simulation
	double ph_step; //!< step on pH in titration simulation

	double e0_min;  //!< minimal value of electrode potential E0 in titration simulation
	double e0_max;  //!< maximal value of E0 in titration simulation
	double e0_step; //!< step on E0 in titration simulation
   
	double e0; //!< Current external electrode potential (V)
	double ph; //!< Current pH value in the system
    
protected:
	std::multimap<AtomGroup*, AltChemState*> res_altst_map;
	typedef std::multimap<AtomGroup*, AltChemState*> altst_map_type;
};


#endif // !PROTONREDOX_H
