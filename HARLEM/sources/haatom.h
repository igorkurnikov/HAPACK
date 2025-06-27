/*! \file haatom.h

    Classes to define Atom object in HARLEM.
    
	\author Igor Kurnikov
    \date 1997-2002

*/

#ifndef HAATOM_H
#define HAATOM_H

#include "haconst.h"  

#include "hastl.h"
#include "hastring.h"

#include "vec3d.h"

class HaResidue;
class HaChain;
class HaMolecule;
class MolSet;
class ChemGroup;
class AtomGroup;
class HaBond;
class AtomFFParam;

template class std::vector<HaBond*>;

const int DUMMY_ELEM = 99; //!< element number used for dummy charges
                                                                              
/* Atom Flags */
const int SphereFlag    =  0x02;     //!< Sphere representation
const int HeteroFlag    =  0x04;     //!< HETATM record
const int HydrogenFlag  =  0x08;     //!< Hydrogen atom 
const int NormAtomFlag  =  0x10;     //!< Flag for a atom if it is not Hydrogen and not HeteroAtom
const int NonBondFlag   =  0x20;
const int BreakFlag     =  0x40;     //!< Break in backbone
const int DisplayedFlag =  0x80;     //
const int HBDonorFlag   =  0x100;    //!< Indicate Hydrogen Donor status
const int HBAcceptorFlag=  0x200;    //!< indicate hydrogen Acceptor status

enum  HYBRIDIZATION { NO_HYBRID, SP_HYBRID, SP2_HYBRID, SP3_HYBRID };

class TiXmlElement;

//! \brief Class to define Atom object
//! \nosubgrouping
class HaAtom : public Vec3D
{
public:
//! \name Constructors, Destructors, Set Parameters 
//@{
  HaAtom();
  HaAtom(const HaAtom& atom_ref);
  virtual ~HaAtom();

  bool SetParamFrom(const HaAtom& atom_ref);
//@}

  int GetElemNo() const;                 //!< Get atom element
  void SetElemNo(const int new_elem_no); //!< Set element of the atom
 
  void SetName(const std::string& atname );       //!< Set atom name
  void SetNameFast(const std::string& atname );   //!< Set atom name is sure it is valid, upper case and no blanks
  const char* GetName() const;           //!< Get the name of the atom

  std::string GetStdSymbol() const;  //!< Get Std Symbol of the atom 

//  bool operator==(const HaAtom& rhs) const;
//  bool operator< (const HaAtom& rhs) const;

//! \name Property Flag Manipulation Utilities
//@{
  virtual void Select();                     //!< Select atom
  virtual void UnSelect();                   //!< Unselect atom
  virtual int Selected() const;              //!< Check if atom is selected
  void SelectBondedHydrogens();              //!< Select Hydrogens bonded to the atom

  bool IsDrawSphere() const;
  void SetDrawSphere(bool set_mode);

  void SetDisplayed( bool set_mode); //!< Set Displayed in Wireframe model
  bool IsDisplayed() const;          //!< (cross if no bonds attached)

  bool IsDummy() const;   //!< Check if atom is dummy( elemno = DUMMY_ELEM )
  void SetDummy();        //!< Set Atom to be dummy ( elemno = DUMMY_ELEM )
  
  bool IsProxy() const;   //!< (for use in residue templates) Atom is used replaced by an atom of the another residue connected to a given residue
  void SetProxy(bool proxy_flag_new = true);  //!< Set proxy_flag

  std::string GetReplacedAtName() const; //!< Get ID string describing atom replaced by this proxy atom ( if IsProxy() == true )
  void SetReplacedAtName( const std::string& repl_atom_name ); //!< Set ID string of the atom replaced by the given proxy atom  

  virtual bool IsHydrogen() const; //!< Check if atom is Hydrogen

  std::string GetDescription() const; //!< Get Atom Description    ( used in residue templates )
  void SetDescription( const std::string& desc ); //!< Set Atom description ( used in residue templates )

//@}
  
//! \name H-Bond donor/acceptor status manupulation
//@{
  bool IsHBDonor() const;              //!< Check if the atom if Hydrogen bond donor
  bool IsHBAcceptor() const;           //!< Check if the atom if Hydrogen bond donor
  void SetHBDonor( bool set_mode);     //!< Set hydrogen bond donor flag for the atom
  void SetHBAcceptor( bool set_mode);  //!< Set hydrogen bond acceptor flag for the atom
  int SetHBStatus(const char* hb_status_str); //!< Set donor/acceptor 
  int GetHBStatusTextStr(std::string& hb_status_str);
//@}

//! \name Molecule Structure Related functions:
//@{
  HaResidue* GetHostRes();             //!< Get the residue Atom belongs to
  const HaResidue* GetHostRes() const;
  HaChain*       GetHostChain();
  const HaChain* GetHostChain() const;
  ChemGroup* GetHostChemGroup();
  HaMolecule* GetHostMol();
  const HaMolecule* GetHostMol() const;
  MolSet* GetHostMolSet();
  const MolSet* GetHostMolSet() const;

  void SetHostRes(HaResidue* new_phost_res);

  using BondIterator       = std::vector<std::shared_ptr<HaBond>>::iterator;
  using BondIterator_const = std::vector<std::shared_ptr<HaBond>>::const_iterator;

  BondIterator Bonds_begin() { return bonds.begin(); }
  BondIterator Bonds_end()   { return bonds.end();   }

  BondIterator_const Bonds_begin() const { return bonds.begin(); }
  BondIterator_const Bonds_end()   const { return bonds.end();   }
 	
  std::vector<std::shared_ptr<HaBond>>& GetBonds() { return bonds;} //!< Get Covalent bonds of the atom

  int GetBondedAtoms(AtomGroup& bonded_atoms); //!< get atoms bonded to the given atom
  AtomGroup GetBondedAtoms();                  //!< get atoms bonded to the given atom
  int GetHBondAcc(AtomGroup& hbonded_acc_atoms); //!< get hydrogen bond acceptor atoms H-bonded to the given H-Bond donor atom
  bool IsBonded(HaAtom& at2);
  
  void RemoveBond(HaBond* pb); //!< Remove Bond from atom bond array 

  static bool CreateBond(HaAtom* aptr1, HaAtom* aptr2); //!< Create bond between 2 atoms
  static bool DeleteBond(HaAtom* aptr1, HaAtom* aptr2); //!< Delete bond between 2 atoms
  
  static HaAtom* AddAtomFromTempl( HaAtom* aptr2, HaAtom* aptr3, HaAtom* aptr4, 
		                           const HaAtom* aptr_templ, const HaAtom* aptr_templ_2, const HaAtom* aptr_templ_3, 
							       const HaAtom* aptr_templ_4); //!< add atom corresponding to atom template (aptr_temp), neigbour atoms and their templates

  static bool SetCoordSubstH(const HaAtom* aptr1, const HaAtom* aptr2, HaAtom* haptr); //!< generate coordinates of a hydrogen connected to aptr1 in the direction of aptr2

  static int GetReachableAtoms(AtomGroup& block_atoms, HaAtom* aptr2, AtomGroup& reach_atoms, int& loop); //!< Get Atoms Accessible moving along bond graph from a given atom
//@}

  static int RegisterAtName( const std::string& at_name); //!< return Reference number of the atom name 
  static void FillStdAtomTypes();     //!< Fill list of names of atoms with standard names
  static StrVec ElemDesc;   //!< list of names of atoms
  static StrIntMap at_name_refno_map;  //!< map of atom names to ref_no 

  bool IsSameName(const HaAtom* aptr_ref) const; //!< Check if atom has the same name
  bool IsMatch(const HaAtom* atempl) const; //!< Check if Atom match template atom 
  
//! \name Output and Related functions:
//@{
  static int AtTypeFromLbl(const std::string & Label); //!< Determine atom element number from the atom label

  bool Print_info(std::ostream &sout, const int level) const;

  enum AtomRefType{ ATOMREF_FULL = 0,  //!< Full atom reference include molecule name  
	                ATOMREF_STD,       //!< Include molecule name only if there are more than 1 molecule in the molecular set
					ATOMREF_NO_MOL,    //!< Do not include molecule name in the atom reference
					ATOMREF_NO_RES, 
					ATOMREF_ELEM_NAME,
					ATOMREF_ELEM_NO}; //!< Atom Reference Type to use in GetRef and FillRef

  std::string GetRef(int ref_type = ATOMREF_FULL ) const;              //!< Get string with an Atom Reference ID
  virtual bool FillRef(char* buf, int ref_type = ATOMREF_FULL ) const; //!< Fill Atom Reference

  int GetSerNo() const; //!< Get Serial number of the atom in the Molecular Set
//@}

//! \name Atom identification functions:
//@{
  static int GetElemNoFromName(const std::string& at_name, const HaResidue* pres = NULL); //!< Get Element Number from element Symbol
  static int GetElemNoFromChar(char ch_fst);        //!< Get Element Number from the first char of Atom Name
  static std::string GetStdSymbolElem(int elem);    //!< Get Standard symbol of the element 
//@}
 
//! \name Force Field related functions 
//@{
  static double ElemVDWRadius(int elem, bool united_atom_flag = false); //!< Get Standard Van der Waals radius of the Element  
  static double StdElemMass(int elem);               //!< Get Atomic mass of the element from the periodical table
  double GetStdMass() const;                         //!< Get Atomic mass of the element of the atom from the periodical table         

  static double ElemDuttonRadius(int elem); //!< Get Atom Radii used in Dutton ET model

  const std::string GetFFSymbol() const;             //!< Get Force Field Symbol of the Atom
  void SetFFSymbol(const std::string& new_ff_symbol); //!< Set Force Field Symbol of the Atom
  
  double GetCharge() const;                  //!< Get Atom Charge (a.u.)
  void SetCharge(double new_charge);         //!< Set Atom Charge (a.u.)
  
  double GetMass() const;                    //!< Get atom mass 
  void SetMass(double new_mass);             //!< Set atom mass 

  double GetVdWRad() const;                  //!< Get atom VdW radius ( Ang )
  void SetVdWRad(double new_rad_vdw);        //!< Set atom VdW radius 

  double GetVdWEne() const;                  //!< Get VdW Energy parameter of the Atom ( kcal/mol ) ( attraction energy of two atoms at 2*r_vdw distance )
  void SetVdWEne(double new_ene_vdw);        //!< Set VdW Energy parameter of the Atom

  std::shared_ptr<AtomFFParam> ps_ff_par;

  // double mass;                 //!< Mass of the Atom in Atomic units
  // double vdw_rad;   //!< atom VdW radius (in Bohr) to compute Van-der-Waals interactions 
  // double ew;        //!< mimimum energy of vdW interaction of the two identical atoms kcal/mol

  std::vector<std::string> comments; //!< Additional comments read from the HIN file
//@}

  TiXmlElement* AddXml(TiXmlElement* parent_element, const char* name = "", int option=0) const; //!< Add a description of the atom to the XML element
  virtual int LoadXml(const TiXmlElement* xml_element, int option=0); //!< Load atom description from the XML element

//!\name Atom Classification:
//@{
  bool IsAlphaCarbon()     const ;
  bool IsSugarPhosphate()  const ;
  bool IsAminoBackbone()   const ;
  bool IsShapelyBackbone() const ;
  bool IsNucleicBackbone() const ;
  bool IsShapelySpecial()  const ;
  bool IsCysteineSulfur() const ;
//@}

//! \name ATOM structure parameters
//@{
  short  x, y, z;              //!< Image Coordinates
  double tempf;                //!< Temperature Factor used also to pass misc atom properties 
  short  col;                  //!< Atom Colour
  std::string label;           //!< Atom Label format
  unsigned char elemno;        //!< Atomic Number
  int refno;                   //!< ElemDesc index number
  unsigned int flag;           //!< Database flags
  double radius;               //!< Atom Radius in Ang to use in Electrostatics calculations
  double image_radius;         //!< Atom Radius in Ang for plotting on the screen
  short  irad;                 //!< Atom Image Radius in pixels 
  double solv_access_area;     //!< Solvent accessible area of the atom (in Ang*Ang)
  
  HYBRIDIZATION hybrid;      //!< Hybridization of the atom

  HYBRIDIZATION GetHybrid() const { return hybrid; } //!< 
  int GetHybridTextStr(std::string& hybrid_str) const; 
  int SetHybrid( const std::string& hybrid_str );              //!< Set Atom Hybridization

  static double StdBondLen(const HaAtom* aptr1, const HaAtom* aptr2); //!< Get Standard Bond length (in Ang) for two atoms 

  int GetNBonds() const;      //!< Get number of bonds of the atom 

//@}

 friend class MolSet;

protected:
  HaResidue* phost_res;   //!< the residue atom belongs to

  std::vector<std::shared_ptr<HaBond>> bonds;  //!< Array of bonds the atom  

  bool proxy_flag;                //!< flag to indicate that atom is used in residue template in place of an atom of an another residue bonded to 
  std::string replaced_atom_name; //!< if (replace_atom_flag) indicate atom id this atom is replacing in the molecular structure
  std::string desc;               //!< description of the atom ( in the residue template )

  void Clear();
};

typedef std::vector<HaAtom*> HaAtomVector;
  
class AtomDoubleMap : public std::map<HaAtom*,double>
{
public:
    
	AtomDoubleMap(const char* new_name = "");
	virtual ~AtomDoubleMap();

	double GetValue(HaAtom* aptr) const;
	int SetValue(HaAtom* aptr, double new_val);
	const char* GetName() const { return name.c_str(); }
	void SetName(const char* new_name) { name = new_name; }

protected:
	std::string name;
};

class AtomAtomMap : public std::map<HaAtom*,HaAtom*>
{
public:
	AtomAtomMap();
	virtual ~AtomAtomMap();

	HaAtom* GetValue(HaAtom* aptr);
	void SetValue(HaAtom* aptr, HaAtom* val);
};


/* CPK Colour Indices
 *  0 Light Grey    1 Sky Blue      2 Red           3 Yellow
 *  4 White         5 Pink          6 Golden Rod    7 Blue
 *  8 Orange        9 Dark Grey    10 Brown        11 Purple
 * 12 Deep Pink    13 Green        14 Fire Brick   15 Mid Green
 */

const int MAXELEMNO = 104;
typedef struct {
           char   symbol[2];
           double covalrad;  // Covalent radius if Bohrs
           double vdwrad;    // Van-der-Waals Radius in Bohrs
		   double mass;      // Atomic Mass from Periodical Table
           int cpkcol;
           char *name;
        } ElemStruct;


#if !defined(HAATOM_CPP)

extern ElemStruct ElementArr[MAXELEMNO];

#endif


#endif /* !HAATOM_H */

