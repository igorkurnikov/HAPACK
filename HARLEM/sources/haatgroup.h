/*!  \file haatgroup.h

   Classes to define atomic groups, residues and chains of residues
   
   \author Igor Kurnikov  
   \date 1998-2002

*/
#ifndef HAATGROUP_H
#define HAATGROUP_H


#include "haio.h"
#include "tinyxml.h"
#include "rapidxml.hpp"
#include "haconst.h"
#include "hacoord.h"
#include "haatom.h"
#include "habond.h"
#include "halinalg.h"


class AtomExpr;
typedef std::set<HaAtom*, less<HaAtom*> > AtomSet;
typedef std::map<HaAtom*, int, less<HaAtom*> > AtomIntMap;
typedef std::map<const HaAtom*, int, less<const HaAtom*> > CAtomIntMap;
// typedef std::map<HaAtom*, HaAtom*, less<HaAtom*> > AtomAtomMap;
typedef std::multimap<HaAtom*, HaAtom*, less<HaAtom*> > AtomAtomMultiMap;
typedef std::map<std::string, HaAtom*, less<std::string> > StrAtomMap;

//#if !defined(RAPIDXML_HPP_INCLUDED)
//	namespace rapidxml { template<class Ch = char> class xml_node; }
//#endif


class AtomIterator : public PointIterator
//! Abstract class for an iterator on collections of atoms 
{
public:

	AtomIterator() {}
	virtual ~AtomIterator() {}

	virtual HaAtom* GetFirstAtom() = 0; //!< Get First Atom in the collection 
	virtual HaAtom* GetNextAtom() = 0;  //!< Get Next  Atom in the collection

	Vec3D* GetFirstPt()  { return (Vec3D*) this->GetFirstAtom(); }
	Vec3D* GetNextPt()   { return (Vec3D*) this->GetNextAtom(); }

};

class AtomIterator_const : public PointIterator_const
//! Abstract class for a constant iterator on collections of atoms 
{
public:

	AtomIterator_const() {}
	virtual ~AtomIterator_const() {}

	virtual const HaAtom* GetFirstAtom() = 0; //!< Get First Atom in the collection 
	virtual const HaAtom* GetNextAtom() = 0;  //!< Get Next Atom in the collection

	const Vec3D* GetFirstPt()  { return (const Vec3D*) this->GetFirstAtom(); }
	const Vec3D* GetNextPt()   { return (const Vec3D*) this->GetNextAtom(); }
	
};

class AtomIteratorGen;

class AtomLoadOptions : public harlem::HashMap
//! Class to define options for loading atoms from stream or file
{
public:
	AtomLoadOptions();
	AtomLoadOptions( const AtomLoadOptions& ref);
	virtual ~AtomLoadOptions();

	virtual void Copy( const harlem::HashMap &ref); 
	virtual harlem::HashMap* clone() const; //!< Create a separate copy of the object and return its pointer 

	void SetStdOptions();

	const std::string& GetDefaultMolName() const { return mol_name_default; } 
	void SetDefaultMolName( const std::string& mol_name) { mol_name_default = mol_name; }

	int ConvertResNames();  //!<  Residue Names conversion option: = 0 - NO COnversion, = 1 Standard conversion: ( HIP -> HIS#PROT,  HIE -> HIS#EPSILON etc ), 
	void SetConvertResNames( int convert_res_names_opt ); //!< Set Residue Conversion option: = 0 - NO COnversion, = 1 Standard conversion: ( HIP -> HIS#PROT,  HIE -> HIS#EPSILON etc )

	bool ToCalcBonds() const { return calc_bonds; } 
	void SetCalcBonds( bool set_par = true ) { calc_bonds = set_par; }

	bool UniqueAtNames() const  { return unique_atom_names; }
	void SetUniqueAtNames ( bool set_par = true ) { unique_atom_names = set_par; }

	friend class ChooseMolFileDlg;
	 
protected:
	std::string mol_name_default;
	bool calc_bonds;  //!< Flag to compute bonds for loaded atoms 
	bool unique_atom_names;  //!< Flag to set unique names for loaded atoms 

};

class AtomSaveOptions : public harlem::SaveOptions
//! Class to define options for saving atoms to stream or file
{
public:
	AtomSaveOptions();
	AtomSaveOptions( const AtomSaveOptions& ref);
	virtual ~AtomSaveOptions();

	virtual void Copy( const harlem::HashMap& ref );
	virtual harlem::HashMap* clone() const; //!< Create a separate copy of the object and return its pointer 

	void SetStdOptions();

	int save_selected;  //!< Flag to save only selected atoms 
	int save_connect;   //!< Flag to save atom covalent bonds connection info 
	int save_transform; //!< Flag to save atom using transformation matrix of the current Molecular View
	int save_atom_ref;  //!< Flag to save atom reference at the end of atom line 
	int save_amber_pdb; //!< Flag to save PDB files with residue names and atom names matching AMBER database
	int save_sep_wat_mol; //!< Flag to save water as separate molecules 

	HaAtom::AtomRefType at_ref_type; //!< Type of the atom reference to save at the end of the atom line

protected:
};

class  AtomContainer: public PointContainer
//! Abstract class for a collection of Atoms
{
public:
    virtual AtomIterator* GetAtomIteratorPtr() = 0; //!< get atom iterator for a given atom collection
	virtual int GetNAtoms() const = 0;   
	virtual int IsMember(const HaAtom* aptr) const = 0;  //!< check if atom belongs to the collection

    AtomIteratorGen __iter__(); //!< Get Atom Iterator ( for Python compatibility )
	AtomIteratorGen GetAtomIterator(); //!< Get General Atom Iterator 

// Virtuals from PointContainer: 
	virtual int IsAtomCollection() const { return TRUE;} //!< check if Point Collection consist of atoms

// Manipulate position and orientation of Atom Collection in 3D space 

	bool GetStdRotMat(HaMat_double& rot_mat); //!< Get Standard Rotation matrix for a given Atom Collection's orientation
	bool GetStdMomInertRotMat(HaMat_double& rot_mat); //!< Get Standard Rotation matrix for a given Atom Collection's orientation based on Principal Moments of Inertia:
	int GetStdPosition(HaMat_double& rot_std, Vec3D& trans_std); //!< Determine std rot matrix and translation for a given molecule
    int SetPosition(const HaMat_double& rot_new, const Vec3D& trans_new);    //!< Set molecular orientations corresponding to a given std rotation and translation
	int	GetStdPositionMomInertia(HaMat_double& rot_std, Vec3D& trans_std);
	int SetPositionMomInertia(const HaMat_double& rot_new, const Vec3D& trans_new);

	int RotateAtoms( const HaMat_double& rot_mat, const Vec3D& cnt); //!< Rotate Atoms around point cnt by rotation matrix rot 
	int TranslateAtoms(const Vec3D& tr_vec );                        //!< Translate Atoms by vector tr_vec

    int SetPosEulerTrans( double phi, double cos_theta, double psi, const Vec3D& trans);  //!< Set the molecule position using Euler angles and translation vectors
	void GetPosEulerTrans( double& phi, double& cos_theta, double& psi, Vec3D& trans); //!< get std translational and rotational (Euler angles) coordinates for a given molecule position
	int SetQuaternionTrans(const Quaternion& q, const Vec3D& trans); //!< Set molecular orientations corresponding to a given quaternion and translation
	void GetQuaternionTrans(Quaternion& q,Vec3D& trans); //!< get std translational and rotational (quaternion) coordinates for a given molecule position

	int SetIntCoordFromStr(const char* int_crd_str); //!< Set Internal coordinates reading them from the string in seq x,y,z, phi,cost,psi

	int SaveXYZFile( const char* fout_name,  const AtomSaveOptions* p_opt = NULL ); //!< Save AtomContainer to XYZ file
	virtual int SaveXYZStream(std::ostream& sout, const AtomSaveOptions* p_opt = NULL ); //!< Save AtomContainer to output stream

	int SaveGROFile(const char* fout_name, const AtomSaveOptions* p_opt = NULL); //!< Save AtomContainer to GROMACS GRO file
	virtual int SaveGROStream(std::ostream& sout, const AtomSaveOptions* p_opt = NULL); //!< Save AtomContainer to GROMACS GRO file stream
};


class AtomIteratorGen
//! General Atom Iterator for a collection of Atoms
{
public: 
	AtomIteratorGen( AtomContainer* pat_cont_new);
	AtomIteratorGen( const AtomIteratorGen& ref);
	virtual ~AtomIteratorGen();

	HaAtom* GetFirstAtom() { return pat_itr->GetFirstAtom(); }
	HaAtom* GetNextAtom()  { return pat_itr->GetNextAtom();  }
	int     GetNAtoms()    { return pat_cont->GetNAtoms();   }

#if defined(SWIG) 
%exception {
try {
	$action
	} catch(std::out_of_range) {
		PyErr_SetString(PyExc_StopIteration,"End of Atoms in the Atom Collection");
		return NULL;
	}
}
#endif
	HaAtom* next(); //!<  Return next atom in the sequence (first on the first call) throw std::out_of_range() if no more atoms ( Python compatibility )
	HaAtom* __next__();

#if defined(SWIG)
%exception;
#endif
	AtomIteratorGen __iter__() const; //!< Get a copy of the iterator ( Python compatibility )

private:
	AtomIterator*   pat_itr;
	AtomContainer*  pat_cont;

	int first_called;
};

class AtomIteratorAtomGroup : public AtomIterator
//! Atom iterator class to browse atoms of the Atom Group
{
public:
	AtomIteratorAtomGroup(AtomGroup* new_p_at_group);
	virtual ~AtomIteratorAtomGroup();
	
	HaAtom* GetFirstAtom(); //!< Return the first atom of the Atom Group (=NULL if no atoms) 
	HaAtom* GetNextAtom();  //!< Return Next atom in the sequence (=NULL if no more atoms)
	
protected:
	vector<HaAtom*>::iterator aitrm; 
	AtomGroup* p_at_group;
};

class AtomIteratorAtomGroup_const : public AtomIterator_const
//! Atom iterator class to browse atoms of the Atom Group
{
public:
	AtomIteratorAtomGroup_const(const AtomGroup* new_p_at_group);
	virtual ~AtomIteratorAtomGroup_const();
	
	const HaAtom* GetFirstAtom(); //!< Return the first atom of the Atom Group (=NULL if no atoms) 
	const HaAtom* GetNextAtom();  //!< Return Next atom in the sequence (=NULL if no more atoms)
	
protected:
	vector<HaAtom*>::const_iterator aitrm; 
	
	const AtomGroup* p_at_group;
};


class AtomGroup : public vector<HaAtom*>, public AtomContainer
//! Class to define a group of atoms
{	
public:
	AtomGroup() {}
	AtomGroup(const AtomGroup& ref_atset);
	AtomGroup(AtomExpr* expr, MolSet* pmset);
	virtual ~AtomGroup();

// Overidables of AtomContainer:

	virtual AtomIterator* GetAtomIteratorPtr();  //!< Return pointer to the corresponing Atom Iterator
	int GetNAtoms() const;                   //!< Return the number of atoms in the list
	int IsMember(const HaAtom* aptr) const;  //!< Check if the atom is a member of the group

// Overidables of PointContainer:

    virtual PointIterator*       GetPointIteratorPtr();
	virtual PointIterator_const* GetPointIteratorPtr() const; 
	int GetNumPt() const  { return size(); }
	
	HaAtom* GetAtomByName(const std::string& at_name); //!< Get Atom of the Group by name 
	HaAtom* GetAtomByIdx( size_t idx ) { return this->at(idx); } //!< Get Atom of the Group by index

	bool InsertAtom(HaAtom* aptr); //!< Insert Atom into Atom Group
	bool DeleteAtom(HaAtom* aptr); //!< Delete Atom from the Atom Group
	int DeleteAtoms(const PtrSet& ptr_set); //!< Delete Atoms from the Group that belong to ptr_set; return the number of atoms deleted
	int DelSelAtoms(); //!< Delete Selected Atoms from the Atom Group, return number of atoms deleted

// Manipulation of Atom Group using Atom expressions:

	void SetFromExpr(AtomExpr* expr, MolSet* pmset );       //!< Build Atom Group from expression 
	void SetFromExprStr(const char* expr_str, MolSet* pmset ); //!< Build Atom Group using expression string
	void AddFromExpr(AtomExpr* expr, MolSet* pmset );          //!< Add to the Atom Group atoms that satisfy expression 
	void AddFromExprStr(const char* expr_str, MolSet* pmset ); //!< Add to the Atom Group atoms that satisfy expression given by the string
	void DeleteAtomsExpr(AtomExpr* expr, MolSet* pmset );     //!< Delete Atoms from the Atom Group atoms that satisfy expression 
	void DeleteAtomsExprStr(const char* expr_str, MolSet* pmset ); //!< Delete Atoms from the Atom Group using expression string
	void KeepOnlyAtomsExpr(AtomExpr* expr, MolSet* pmset );     //!< Keep only Atoms from the Atom Group that satisfy expression 
	void KeepOnlyAtomsExprStr(const char* expr_str, MolSet* pmset ); //!< Keep only Atoms from the Atom Group that satisfy expression given by string

    const char* GetID() const;           //!< Return Group ID
	void SetID(const std::string& new_id);  //!< Set ID of the Atom Set

	int HasSelectedAtoms(); //!< Check if some of the atoms of the group are selected
	void SelectAtomsAll();   //!< Select All Atoms in the residue

protected:
	std::string id;	
};

class MolSet;

// Beware AddAtom is defined in winbase.h (WIN32)!!

class ChemGroup : public AtomGroup
//!  Class to define Chemical(Functional) Atom Group object
{
public:
  ChemGroup();
  ChemGroup(MolSet* new_phost_mset, const char* new_id= "" );
  virtual ~ChemGroup();
          
  bool FillRef(char* buf) const; //!< Fill Text Reference of the group 

  double GetProtect() const;    //!< return screening factor 
  bool SetProtect(const double new_protect);
  
  friend class MolSet;
 
  static std::string GetIDFromRef(const std::string& buf);

  bool Print_info(ostream &sout, const int level) const;

protected:

   MolSet* phost_mset;

   double protect; // screening factor of the group in the fragment
                  // from 0 to 1;  1.0 - non-perturbed group  
};

class HaChain;

// enum ResCodes {R_ALA=0, R_GLY, R_LEU, R_SER, R_VAL, R_THR, 
//                         R_LYS, R_ASP, R_ILE, R_ASN, R_GLU, 
//					       R_PRO, R_ARG, R_PHE, R_GLN, R_TYR,
//	    		            R_HIS, R_CYS, R_MET, R_TRP};

//! Atom Types to add
enum ADD_ATOM_TYPE 
{                    
    ADD_ALL_ATOMS = 0x0,
    ADD_HYDROGENS = 0x1,
    ADD_POLAR_HYDROGENS = 0x2,
	ADD_HEAVY_ATOMS = 0x4 
};

class ForceFieldType;
typedef AtomIteratorAtomGroup AtomIteratorResidue;
typedef std::list<AtomGroup> AtomGroupList;

//!  Class to define Residue in a polymer or biopolymer chain   
class HaResidue : public AtomGroup
{
public:
	HaResidue();
//	HaResidue(const HaResidue& ref_res);
	HaResidue(HaChain* new_phost_ch);
	virtual ~HaResidue();
	
	bool SetParamFrom(const HaResidue& res_ref); //!< Copy residue parameters from another residue
	
	int GetNAtomsNonProxy() const; //!< Get the number of Non-proxy atoms in the residue

	HaAtom* AddNewAtom(); //!< Add new atom to the residue

	HaAtom* GetAtomByName(const std::string& atname);      //!< Get Atom of the residue by its name
	const HaAtom* GetAtomByName(const std::string& atname) const;   //!< Get Atom of the residue by its name (const version)
	 
	HaChain* GetHostChain() { return phost_ch; } //!< Get the chain of the residue belongs to 
	bool SetHostChain(HaChain* new_phost) { phost_ch=new_phost; return true;} //!< Set the chain of the residue belongs to 
	
	HaMolecule* GetHostMol();  //!< Get the molecule the residue belongs to
	const HaMolecule* GetHostMol() const;  //!< Get const pointer to the molecule the residue belongs to

	MolSet* GetHostMolSet();                //!< Get the molecule set the residue belongs to
	const MolSet* GetHostMolSet() const;  //!< Get const pointer to the molecular set the residue belongs to

	HaResidue* GetNextResInChain();        //!< Get Next Residue in Chain
	HaResidue* GetNextResInChain() const;  //!< Get Next Residue in Chain
	HaResidue* GetPrevResInChain();        //!< Get Previous Residue in Chain
	HaResidue* GetPrevResInChain() const;  //!< Get Previous Residue in Chain

	bool IsBonded(HaResidue* res2) ; //!< Check if this residue bonded to another residue 
	int HasBackBHBond(HaResidue* res2); //!< Check if amino acid residue has a backbone H-bond to another aminoacid 

	const char* GetName() const;  //!< Get the basic name of the residue
	void SetName(const std::string& res_name, int convert_res_names = 0); //!< Set Residue Name 
	void SetNameModifier( const std::string& new_name_mod) { NameModifier = new_name_mod; } 
    
	const char* GetNameModifier() { return NameModifier.c_str(); }  //!< Get Residue Name Modifier
	std::string GetFullName() const; //!< Get Name of the residue with the modifier 

//! Extract Short Residue name from the full residue name (with modifier)  
	static std::string GetResNameFromFullName(const char* res_full_name);

	std::string GetRef() const;  //!< Get the text reference of the residue
	virtual bool FillRef(char* buf,int mode = 0) const; //!< Write the text reference of the residue to the buffer

	AtomIntMap  GetAtomSeqNumMap(); //!< Get the map of atoms to sequence atom numbers in the molecule
	CAtomIntMap GetAtomSeqNumMap() const; //<! Get the map of atoms to sequence atom numbers in the molecule - const version

	bool SetUniqueAtomNames();                  //!< Set unique names to atoms of the residue 
	std::string GetUniqueAtomName(int elem_no); //!< Get unique atom name for element elem_no 
	bool SplitResidue();

    int SetStdCharges();     //!< Set atomic charges corresponding to the template of the residue
	int InterpolResParams(const char* res_name_1, const char* res_name_2, double weight_1); //!< Make the residue parameters intermediate between two templates 

	HaResidue* GetTemplate(); //!< Get template corresponding to a current residue full name
	int CheckStruct(); //!< Check if the structure of the residue correspond to a template
	int CheckStructMortLib(const ForceFieldType& ff_type); //!< Check the structure of the residue against residue database in Mort Library Format
	
	int AddMissingAtoms(ADD_ATOM_TYPE atom_type ); //!< Add Atoms to the residue based on the residue template in DB
	int AddWaterHydrogens();                      //!< Add Hydrogens to a water molecule residue

protected:
	int AddMissingAtoms_2(HaResidue* prtempl, ADD_ATOM_TYPE atom_type);

public: 

	int GetSerNo() const { return serno; }  //!< Get a serial number if the residue in the chain

	static const char* GetResNameInTable(const int j); //!< Get residue name in the table of residue names ResNames by index
// Structure definition functions:

	bool IsAmino()       const;  //!< Check if the residue is an aminoacid 
	bool IsAminoNucleo() const;  //!< Check if the residue is an aminoacid or nucleic acid base
	bool IsNucleo()      const;
	bool IsProtein()     const; 
	bool IsDNA()         const;
	bool IsSolvent()     const;
	bool IsWater()       const;
	bool IsIon()         const;
	bool IsPyrimidine()  const;
	bool IsPurine()      const;
	bool IsRNA()         const;
	bool IsProline()     const;
	bool IsHistidine()   const;
	bool IsCysteine()    const;
	bool IsAdenine()     const;
	bool IsCytosine()    const;
	bool IsGuanine()     const;
	bool IsThymine()     const;
	bool IsCoenzyme()    const;
	bool IsTerm()        const;

	int CalcStdCrdSys(int fit_std_geom = FALSE); //!< Calculate standard coordinate system

// Geometry Calculations:

	static double CalcPhiAngle(const HaResidue* prev, const HaResidue* curr);
	static double CalcPsiAngle(const HaResidue* curr, const HaResidue* next );
	
	HaChain*  phost_ch;          //!< chain the residue belongs to 
	short serno;                 //!< Residue serial number   
	double width;                //!< Ribbon Width          
	short col1;                  //!< Ribbon Colour #1      
	short col2;                  //!< Ribbon Colour #2      
	char  insert;                //!< PDB insertion code    
	int   refno;                 //!< Residue Name index in ResNames table  
	unsigned char  struc;        //!< Secondary Structure   
	unsigned char  flag;         //!< Database flags        

	std::string NameModifier;    //!< Modifier of the residue name

//    vector<AltChemState> alt_res_states; //!<  Alternative Residue States (protonated, deprotonated etc) 
    Vec3DValArray std_crd_sys;       //!< Standard coordinate system 

	static int RegisterResName( const std::string& res_name); //!< return Reference number of the residue name
	static void InitStdResNames(); //!< Initialize standard residue names
	static void InitResSynonym(); //!< Initialize standard residue name synonyms

	static std::vector<std::string> ResNames;               //!< table of Residue Names 
	static StrIntMap res_name_refno_map;       //!< map of residue names to ref_no 
	static StrStrMap ResSynonym_to_std;        //!< the table of Residue Name Synonyms from AMBER to STD names( with modifies) 
	static StrStrMap ResSynonym_std_to_AMBER;  //!< the table of Residue Name Synonyms from STD names( with modifies) to AMBER

protected:
	void Clear();	
};

class AtomIteratorChain : public AtomIterator
//! Atom iterator class to browse atoms of a chain (HaChain class)
{
public:
	AtomIteratorChain();
	AtomIteratorChain(HaChain* new_chain);
	virtual ~AtomIteratorChain();

	int SetForChain(HaChain* new_chain);
	
	HaAtom* GetFirstAtom(); //!< Return the first atom of the chain (=NULL if no atoms) 
	HaAtom* GetNextAtom();  //!< Return Next atom in the chain (=NULL if no more atoms)
	
protected:
	std::vector<HaResidue*>::iterator ritrm;
	vector<HaAtom*>::iterator aitr;
	
	HaChain* chain;
};

class AtomIteratorChain_const : public AtomIterator_const
//! Atom iterator class to browse atoms of a chain (HaChain class)
{
public:
	AtomIteratorChain_const();
	AtomIteratorChain_const(const HaChain* new_chain);
	virtual ~AtomIteratorChain_const();

	int SetForChain(const HaChain* new_chain);
	
	const HaAtom* GetFirstAtom(); //!< Return the first atom of the chain (=NULL if no atoms) 
	const HaAtom* GetNextAtom();  //!< Return Next atom in the chain (=NULL if no more atoms)
	
protected:
	std::vector<HaResidue*>::const_iterator ritrm;
	vector<HaAtom*>::const_iterator aitr;               
	
	const HaChain* chain;
};

class ResidueIteratorChain
//! Residue iterator class to browse residues of the chain
{
public:
	ResidueIteratorChain(HaChain* new_chain);
	ResidueIteratorChain(const ResidueIteratorChain& ritr_ref);
	virtual ~ResidueIteratorChain();
	
	HaResidue* GetFirstRes(); //!< Return the first atom of the Chain (=NULL if no residues) 
	HaResidue* GetNextRes();  //!< Return Next atom in the chain (=NULL if no more residues)
	
protected:
    std::vector<HaResidue*>::iterator res_itr;
	HaChain* chain;
};

class HaChain : public AtomContainer
//!< Class to define chain of residues 
{
public:
  HaChain();
  HaChain(HaMolecule* new_phost_mol, const char new_ident=' ');

  virtual ~HaChain();

  bool SetParamFrom(const HaChain& chain_ref); //!< copy parameters from reference a reference chain 
  
  HaResidue* AddResidue(int res_ser_no);  //!< Add residue to the chain with the serial number
  int GetUniqResSerNo(int term_res_flag=0) const; //!< Get residue id number not yet in the chain

  bool SetUniqueResNo();  //!< Set unique id(serial) numbers to the residues in the chain

  HaResidue* GetFirstRes(); //!< Get First Residue of the chain 
  HaResidue* GetResBySerNo(const int res_ser_no); //!< Get the residue by its serial number

  int GetNRes() const { return res_arr.size(); }  //!< Get the number of residues in the chain

// Overidable of AtomContainer:

	AtomIterator* GetAtomIteratorPtr() { return new AtomIteratorChain(this); }
	AtomIterator_const* GetAtomIteratorPtr() const { return new AtomIteratorChain_const(this); }
	int GetNAtoms() const; //!< Return the number of atoms in the chain
	int IsMember(const HaAtom* aptr) const;       //!< Check if the atom is a member of the chain

// Overidable of PointContainer:

    PointIterator* GetPointIteratorPtr();
	PointIterator_const* GetPointIteratorPtr() const;
	int GetNumPt() const;

// RASMOL Chain structure parameters
   
  char ident;                      //!< Chain identifier             
 
  vector<HaResidue*> res_arr;   //!< Array of Residues in the chain
  multimap<int, HaResidue*, less<int> >  res_map; //!< Map of serial numbers to Residues 
 
  HaMolecule* GetHostMol() { return phost_mol; }  //!< Get the molecule chain belongs to
  const HaMolecule* GetHostMol() const { return (const HaMolecule*) phost_mol; } //!< Get the molecule chain belongs to (const version)
  void SetMolHost(HaMolecule* new_phost_mol) { phost_mol=new_phost_mol; }

protected:
  
  void SetDefaultParam();

  HaMolecule* phost_mol; //!< The Molecule chain belongs to
};

class PeriodicUnitInfo
//! class to describe periodicity of the system
{
public:

	PeriodicUnitInfo();
	PeriodicUnitInfo(const PeriodicUnitInfo& ref);
	virtual ~PeriodicUnitInfo();

	int SetBox(double a, double b, double c, 
		       double alpha_n = 90.0*DEG_TO_RAD, double beta_n = 90.0*DEG_TO_RAD, double gamma_n = 90.0*DEG_TO_RAD); //!< Set parameters of periodic boundary

	int IsSet() const;          //!< Check if Periodicity Information is set
	int IsOrthogonal() const;   //!< Check if the periodical unit cell is orthogonal
	bool IsOctahedron() const;   //!< Check if the periodical unit cell is octahedral
	int IsValid() const;        //!< Check if Periodical Boundary box params are valid
	void Set( bool to_set = true );     //!< Set/Unset Periodical Boundary Conditions
	void SetOctahedron( bool to_set = true ); //!< Set/Unset Octahedron Boundary conditions  
	int SetStdBox(AtomContainer* at_coll); //!< Compute Periodical Box for AtomContainer 

	std::string spacegroup;
	int orthogonal_flag;
	bool octahedron_flag;
	
	Vec3DValArray ucell;         //!< Unit cell vectors
	Vec3DValArray recip_ucell;   //!< Reciprocal cell unit vectors

	double GetA() const { return pbc_a; } //!< Get dimension of the unit cell in the a(x) direction  (Ang)
	double GetB() const { return pbc_b; } //!< Get dimension of the unit cell in the b(y) direction  (Ang)
	double GetC() const { return pbc_c; } //!< Get dimension of the unit cell in the c(z) direction  (Ang)
	double GetAlpha() const { return alpha; } //!< Get Angle between a and b  (radians)
	double GetBeta()  const { return beta; } //!<  Get Angle between a and c  (radians)
	double GetGamma() const { return gamma; } //!< Get Angle between b and c  (radians)

protected:
	
	double pbc_a;  //!< Dimension of the unit cell in the a(x) direction  (Ang)
	double pbc_b;  //!< Dimension of the unit cell in the b(y) direction  (Ang)
	double pbc_c;  //!< Dimension of the unit cell in the c(z) direction  (Ang)

	double alpha;  //!< Angle between a and b  (radians)
	double beta;   //!< Angle between a and c  (radians)
	double gamma;  //!< Angle between b and c  (radians)

public:
	double ucell_vol;  //!< Unit cell volume (Ang**3)
	double ucell_sph;  //!< radius of the largest sphere to fit the unit cell

	double cut_factor[3];
    double reclng[3];
};


class CrdSnapshot
//!< Snapshot of atom coordinates of Atom Container
{
public:
	CrdSnapshot(AtomContainer* p_at_coll_new, const std::string& name_new = "SNAPSHOT");
	virtual ~CrdSnapshot();

	void Clear(); //!< clear content of the snapshot

	int SaveCrd(const HaVec_double& crd);   //!< Save coordinate array to the snapshot
	int SavePBox(const HaVec_double& pbox); //!< Save Periodical Pox info to the snapshot
	int SaveCurrentAtomCrd();   //!< Save current atom coordinates (including Periodical box is available) of the host Atom Container to the snapshot 
	int SetAtomCrd();           //!< Set coordinates of host atom collection from coordinates values of the snapshot
	int IsValid() const;  //!< Check if snapshot coordinates are valid( dimension corresponds to the size of the AtomContainer p_at_cont)
	int HasPBox() const;  //!< Check if snapshot has periodical box info 
	std::string GetName() const;        //!< Get Name (ID) of the snapshot
	void SetName(const std::string& name_new); //!< Set Name (ID) of the snapshot
	std::string GetDesc() const;        //!< Get Descrition of the snapshot
	void SetDesc(const std::string& desc_new); //!< Set Description of the snapshot

	HaVec_double GetCrd()  const { return crd;  } //!< Get stored atom coordinate values
	HaVec_double GetPBox() const { return pbox; } //!< Get stored periodical box dimensions

	int LoadXMLNode( rapidxml::xml_node<>* node, const harlem::HashMap* popt ); //!< Load Snapshot from XML node

	int OnDelAtoms( AtomContainer& at_del ); //!< Modify Coordinate snapshot when atoms were deleted

	AtomContainer* p_at_cont; 

protected:
	std::string name;  //!< Snapshot Name (ID)
	std::string desc;  //!< Snapshot description
	HaVec_double crd;  //!< Atom coordinates of the container
	HaVec_double pbox; //!< periodical box info, if size=3 only box side sizes if size = 6 - directional angles (in radians)
};

class CrdSnapshotIterator
//!< Iterator on CrdSnapshot objects associated with Molecular Set
{
public:
	CrdSnapshotIterator( MolSet* pmset );
	CrdSnapshotIterator(const CrdSnapshotIterator& ref);
	virtual ~CrdSnapshotIterator();

	CrdSnapshotIterator& operator=( const CrdSnapshotIterator& ref);
#if defined(SWIG) 
%exception {
try {
	$action
	} catch(std::out_of_range) {
		PyErr_SetString(PyExc_StopIteration,"End of Coordinate Snapshots");
		return NULL;
	}
}
#endif
	CrdSnapshot* next(); //!<  Return next coordinate snapshot  (throw std::out_of_tange if no more atoms)  for Python compatibility
#if defined(SWIG)
%exception;
#endif
	CrdSnapshotIterator __iter__() const; //!< Get a copy of the iterator ( for Python compatibility )

//	CrdSnapshot* GetFirst(); //!< Get First Atom Coordinate Array associated with Molecular Set
//	CrdSnapshot* GetNext();  //!< Get Next  Atom Coordinate Array associated with Molecular Set

protected:
	bool first_call;
	std::vector<CrdSnapshot*>::iterator itr_curr;
	std::vector<CrdSnapshot*>::iterator itr_end;
};

class MutationMap
//!< Mutaton Map between atoms of two molecules
{
public:

	MutationMap(MolSet* pmset);
	virtual ~MutationMap();

	void Clear();
	int LoadArbalestMutMap( std::string fname ) noexcept; //!< Load Mutation Map from Arbalest XML file

	std::map<HaAtom*, HaAtom*> atom_atom_map;

	HaMolecule* pmol1;
	HaMolecule* pmol2;

	MolSet* pmset;
};


#endif /* !HAATGROUP_H */

