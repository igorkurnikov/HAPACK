/*!  \file mm_elements.h

    Molecular Mechanics Model elements 

    \author Igor Kurnikov 
    \date 2010-

*/
#ifndef MM_ELEMENTS_H
#define MM_ELEMENTS_H

#include "haatgroup.h"
#include "mm_params.h"


class AtomFFParam
//! Atom Force Field parameters 
{
public:
	AtomFFParam();
	//AtomFFParam(HaAtom* aptr_ref_new);
	virtual ~AtomFFParam(); 

	void Clear(); //!< Clear content
	int LoadXml(const TiXmlElement* xml_element, std::string at_name, int option=0 ); //!< Load Atom FF Template from XML element

	int HasDipole() const;   //!< Check if atom has electic dipole parameter
	int HasQPole()  const;   //!< Check if atom has electric quadrupole parameter
	int HasPolar()  const;   //!< Check if atom has polarizability parameters 
	bool HasHPolar() const;   //!< Check if atom has hyperpolarizability parameters 
	int HasFrameAtomNames() const; //!< Check of atom has defined atoms names for local reference frame
	int IsFrameSet()  const; //!< Check if atotms defining local refernce frame has been set
	bool AlterPolarizability() const; //!< Indicate if atom is affecting polarizabilty of surrounding atoms
	bool IsPolarPositionDep() const; //!< Indicate if atom polarizability if affected by the presence of other atoms
	bool HasScreenPolar() const;   //!<   Indicate if polarizabilty screening coefficent ( sqrt(thole_coef) ) is set
	
	double GetCharge() const;        //!< Get Electric Charge
	HaVec_double GetDipole() const;  //!< Get Electric Dipole
	HaVec_double GetQPole()  const;  //!< Get Electric Quadrupole
	int IsChiralFrame() const;       //!< Check if local frame is chiral
	int IsBisectFrame() const;       //!< Check if local frame vector Z is defined as a bisector 
	bool SetFrameFromAtomNames(HaAtom* aptr_ref);  //!< Set Pointers to atoms defining reference frame using stored array of atom names

public:
	double charge;          //!< electrical charge of the atom
	HaVec_double dipole;    //!< electrical dipole of the atom
	HaVec_double qpole;     //!< electrical quadrupole of the atom

	double polar;             //!< polarizability of the atom
	double hpolar;            //!< hyperpolarizability of the atom 
	double damp_polar_strength;    //!< Strength of damping polarizability of surrounding atoms ( if > 0 polarizability of neigbouring atoms will be reduced)
	double damp_polar_sensitivity; //!< Sensitivity to damping atom polarizability ( if > 0 polarizability of atom will be affected by the presence of other atoms)
	double damp_polar_rad;         //!< Atom radius for polarizabiltiy damping interactions
	double screen_polar;           //!< Polarity screenning coeffient ( if > 0 than sqrt(thole_coef) or not set )

	double mass;  //!< Atom Mass 
	double rad_vdw;  //!< VdW radius of the atom
	double ene_vdw;  //!< VdW energy of the atom

	//std::string at_name;    //!< Atom name in the residue or residue force field template 
	std::string ff_symbol;  //!< Force field symbol of the atom
	std::string ff_polar_symbol;  //!<  Force field symbol of the atom for electrostatic and polarization parameters
	
	StrVec    frame_atom_names;  //!< Array of atom names defining local frame of the atom 
	AtomGroup frame_atoms;       //!< Array of atom pointer defing local frame of the atom
	int bisect_flag;             //!< Flag to specify that Z-axis of the local frame of the atom is defined as a bisector 

	// HaAtom* aptr_ref;
};

class ResFFTemplate
//! Class for Residue Force Field Template
{
public:
	ResFFTemplate(HaResidue* p_res_templ_new);
	virtual ~ResFFTemplate();

	void Clear(); //!< Clear All content

	bool LoadXml(const TiXmlElement* xml_element, int option=0 ); //!< Load Residue FF Template from XML element

	shared_ptr<AtomFFParam> GetAtomFFParam(const std::string& at_name); //!< Get Atom Force Field Parameters by Atom Name 
	bool SetAtomFFParam(const std::string& at_name, shared_ptr<AtomFFParam> p_at_ff_param); //!< Set Atom Force Field parameters

	void SetResFFVersion(const std::string& res_ff_version_new); //< Set Version string of the Residue force field template 
	string GetFullName(); //!< Get Full Name of the residue template

	HaResidue* GetResTemplate(); //!< Get Residue structure template associated with the Force Field residue template
	
	vector<StrVec> bonds;
	vector<StrVec> angles;
	vector<StrVec> dihedrals;
	vector<StrVec> improper_dihedrals;

private:
	string res_name;       //!< Full Residue Name in Residue Database
	string res_ff_version; //!< Identifier for the version of the force field for the residue 

	HaResidue* p_res_templ;
	map<string,shared_ptr<AtomFFParam>> at_name_ff_param_map; //!< Map of atom names to atomic force field parameters
};

class MMBond
//! Valence bond in the Molecular Mechanics model
{
public:
	MMBond();
	MMBond(HaAtom* new_pt1, HaAtom* new_pt2);
	MMBond(const MMBond& bond_ref);
	~MMBond();

	bool operator==(const MMBond& rhs) const;
	bool operator< (const MMBond& rhs) const;

	double r0; // Equilibrium distance (in Ang)
	double fc; // force constant

	int set_type; 

public:
	HaAtom* pt1;
	HaAtom* pt2;

};

class MMValAngle
//! Valence Angle in Molecular Mechanics Module
{
public:
	MMValAngle();
	MMValAngle(HaAtom* new_pt1, HaAtom* new_pt2, 
			HaAtom* new_pt3);
	virtual ~MMValAngle();

public:
	HaAtom* pt1;
	HaAtom* pt2;
	HaAtom* pt3;

	double a0; //!< Equilibrium angle (in degrees)
	double fc; //!< force constant

	int set_type; 

	bool operator == (const MMValAngle& rhs) const;
	bool operator <  (const MMValAngle& rhs) const;
};

class MMDihedral
//! Dihedral Angle in Molecular Mechanics Module
{
public:
	MMDihedral();
	MMDihedral(HaAtom* new_pt1, HaAtom* new_pt2, 
			HaAtom* new_pt3, HaAtom* new_pt4, bool improper_flag = false);
	virtual ~MMDihedral();

public:
	bool improper; //!< flag to indicate improper torsional angle
	bool calc_14;  //!< flag to whether to calculate 1-4 interactions for end atoms of the torsion

	HaAtom* pt1;
	HaAtom* pt2;
	HaAtom* pt3;
	HaAtom* pt4;

	int set_type; 

// Torsional potential is expressed as: 	
// E_tors = ( PK/IDIVF) * (1 + cos(PN*phi - PHASE))

	int GetNTerms() const { return pk.size(); } //! Get the number of terms in a dihedral angle 
	void ClearParams() { pn.clear(); phase.clear(); pk.clear(); idivf.clear(); }  //!< clear FF parameters
	
	vector<double> pn;    //!< Periodicity of the dihedral of the given type 
	vector<double> phase; //!< phase of the dihedral of the given type
	vector<double> pk;    //!< force constant (potential depth) 
	vector<double> idivf; //!< constant to multiply pk to get a potential depth   

	int AddTerm( double pn_new, double phase_new, double pk_new, double idivf_new = 1.0 );

	bool operator == (const MMDihedral& rhs) const;
	bool operator <  (const MMDihedral& rhs) const;

};

typedef vector<MMDihedral> MMDihedralArray;

class AtomContact
//! Atom-Atom contact(constraint) harmonic, VdW (6-12) or (10-12) or Coulomb (q_i*q_j/r_ij) in the Molecular Mechanics model
{
public:
	AtomContact();
	AtomContact(HaAtom* new_pt1, HaAtom* new_pt2, const AtomContactType& cnt_type = AtomContactType::HARMONIC_CNT );
	AtomContact(const AtomContact& at_cnt_ref);
	~AtomContact();

	bool operator==(const AtomContact& rhs) const;
	bool operator< (const AtomContact& rhs) const;

	HaVec_double cf; //! Coef params of Atom-Atom contact - the meaning depends on cnt_type

//  cnt_type == HARMONIC CNT
//  cf[0] - eq_dist (X0) (Ang) 
//  cf[1] - force_const (K) (kcal/mol/(Ang*Ang)) 
//
//  cnt_type == VDW_CNT_6_12 or cnt_type == VDW_CNT_6_12_NO_REP
//  cf[0] - coef at R^-12  ( kcal/mol * Ang^6)
//  cf[1] - coef at R^-6   ( kcal/mol * Ang^12)
//
//  cnt_type == VDW_CNT_10_12:
//  cf[0] - coef at R^-12  ( kcal/mol * Ang^6)
//  cf[1] - coef at R^-10  ( kcal/mol * Ang^12) 
	
	AtomContactType cnt_type; //!< Atom-Atom contact (distance constraint) type  

	int SetParamsEneR( double ene_min, double rmin); //!< Set parameters with minimal energy value (absolute value in kcal/mol) and optimal distance (Ang)  
	double GetRMin() const;
	double GetEneMin() const;
	double GetHarmForceConst() const;

	bool IsHarmonic() const { return (cnt_type == cnt_type.HARMONIC_CNT); }
	
	int set_type; 

public:
	HaAtom* pt1;
	HaAtom* pt2;
};


#endif // end if !defined(MM_ELEMENTS_H ) 