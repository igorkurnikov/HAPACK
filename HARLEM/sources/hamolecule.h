/*! \file hamolecule.h

    Classes to define a molecule object in HARLEM.

    \author Igor Kurnikov 
 
    \date 1997-2008

*/
#ifndef HAMOLECULE_H
#define HAMOLECULE_H

#include "hastring.h"
#include "hastl.h"
#include "halinalg.h"
#include "haatom.h"
#include "habond.h"
#include "haatgroup.h"
#include "object3d.h"
#include "command.h"

class HaAtom;
class HaBond;
class HaQCMod;
class HaMolView;
class GauFile;

const int SourceNone = 0;
const int SourcePDB  = 1;
const int SourceCalc = 2;

				  
const int FeatHelix = 1;
const int FeatSheet = 2;
const int FeatTurn  = 3;


class SecStructElement
{
public:
	SecStructElement();
	virtual ~SecStructElement();

	int init, term;
	char chain;
	char type;

	bool operator == (const SecStructElement& rhs);
	bool operator < (const SecStructElement& rhs);

protected:
};

class ResidueIteratorMolecule;

//!  Class to define a single Molecule object
class HaMolecule: public Object3D, public AtomContainer
{
public:
	HaMolecule(MolSet* new_phost_mset, std::string new_name="MOL");
	HaMolecule(const HaMolecule& Mol_ref);
    virtual ~HaMolecule();

	MolSet* GetHostMolSet() { return phost_mset; }
	const MolSet* GetHostMolSet() const { return phost_mset; }

	bool AddMolCopy(HaMolecule& Mol_ref, bool create_new_chain = true, AtomAtomMap* ptr_atom_map = NULL);

	static int AttachFragment( HaAtom* catom_host, HaAtom* catom_frag );
	int CombineMolecules(HaMolecule* frag_mol, HaAtom* catom_host, HaAtom* catom_frag );


protected:
	MolSet* phost_mset;

	std::string mol_name;
	int serno; 

public:

	std::string GetName() const; //!< Get Molecule Name
	int GetSerNo() const;    //!< Get Molecule Serial Number in the Molset
	std::string GetRef() const;  //!<  get a text reference for a molecule 
	bool FillRef(char* buf,int mode = 0) const; //!< Fill string with a molecule text reference

//! virtual function from Object3D

	virtual bool SetObjName(const char* new_name);
	
	virtual int RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt ); //!< Rotate object around the center
	virtual int RotateObjFromWorld( const HaMat_double& rot_mat, const Vec3D& cnt ); // !< Rotate object with respect to world

	virtual int Translate( const Vec3D& tr_vec ); //!< Translate the molecule by (dx,dy,dz)

//! Display Related functions:

	int SetAtomScreenCoord(HaMolView* pview);
	
	bool SetUniqueAtomNames();

//! Atoms Related functions:

	bool InitAtoms(GauFile& gfile);      //!< load atom coordinates from Gaussian rwf file	                                      

	HaAtom* GetAtomByRef(const char* at_ref);             //!< Get an atom of the molecule by its text reference

// virtual overidable finctions of AtomContainer

    virtual int GetNAtoms() const;                   //!< return the number of atoms in molecule
	virtual AtomIterator* GetAtomIteratorPtr();         //!< Create AtomIterator for the molecule
    virtual int IsMember(const HaAtom* aptr) const;  //!< Check if the atom belongs to the Molecule

// Overidables of PointContainer:

    virtual PointIterator*       GetPointIteratorPtr();
	virtual PointIterator_const* GetPointIteratorPtr() const; 
	int GetNumPt() const  { return GetNAtoms(); }

// Bonds Related functions:

	int GetNBonds()  const;  //!< Get Number of Covalent bonds  
	int GetNHBonds() const;  //!< Get Number of Hydrogen bonds  
	int GetNSSBonds() const;  //!< Get Number of SS bonds  

public:

	int SetTermResNames(); //!< Set proper terminal residue name modfiers (N & C terminal for proteins) (3' and 5' for Nucleic Acids)
	int SetHISNames(); //!< Set His Modifier Name based upon the presence of HD1 and HE2 hydrogens
    int SetCysBridgeNames(); //!< Set proper residue name for cys involved in disulfide bridges 

public:   	
    
	bool Print_info(ostream &sout, const int level) const; //!< Print information about the molecule

// Residues Related Functions:

	HaResidue* AddChainAndResidue(); //!< Create first chain and first 'RES' residue
	int GetNRes() const;     //!< get the number of residues in the molecule

	HaResidue* GetResByRef(const std::string& res_str); //!< Get a residue of the molecule by its text reference

// Chain Related functions	
	int GetNChains() { return Chains.size(); } //!< Get the number of chains in the molecule
	HaChain* AddChain(char ident); //!< Add Blank Chain to the list with identificator ident

	char GetChainIdentMax(); //!< Get Maximal char identificator of the molecule
	HaChain* GetFirstChain(); //!< get a pointer to the first Chain of Molecule

	HaChain* GetChain(const char chain_id ); //!< Get a chain by id 

	bool FixChainsIdent(); //!< set unique chain ident on the molecule

// Secondary Structure manipulation functions:
public:
	void DescribeMolecule();

	void Renumber(int start );
	AtomIntMap GetAtomSeqNumMap();  //!< Get the map of atoms to sequence atom numbers in the molecule
	CAtomIntMap GetAtomSeqNumMap() const; //!< Get the map of atoms to sequence atom numbers in the molecule - Const version
	void DescribeSequence();

	static int SeqFormat; //!< Format Sequence in DesribeSequence() function
	
	SecStructElement* AddFeature();       //!< Add Secondary Structure element of a given type
	void DeleteFeatures(const int itype); //!< Delete Secondary Structure elements of the given type
	int GetNumFeatures(const int itype) const;
	void UpdateFeature(SecStructElement* ptr, int mask);
	void ProcessFeatures();
	
	bool IsSecStructFound() { return sec_struct_found; }
	bool sec_struct_found;  //!<  flag to indicate that secondary structure is found 

	int structsource;  //!< indicate whether secondary structure were loaded from PDB or computed
	
	list<SecStructElement> Features;
	
// Molecule data:
private:
	list<HaChain> Chains;      //!< List of Chains	

public:
	
	std::string classification;
	std::string identcode;

public:

    friend class ChainIteratorMolecule;
	friend class AtomIteratorMolecule;
	friend class AtomIteratorMolecule_const;
	friend class ResidueIteratorMolecule;
	friend class ChainIteratorMolSet;
	friend class AtomIteratorMolSet;
	friend class ResidueIteratorMolSet;
	friend class ResidueIteratorMolSet_const;
	friend class MolSet;

	//	typedef AtomIteratorMolecule AtomIterator; //!<  Atom iterator type for the molecule
	typedef ResidueIteratorMolecule ResidueIterator; //!<  Residue iterator type for the molecule

};


class AtomIteratorMolecule: public AtomIterator
//! Atom iterator class to browse atoms of the molecule
{
public:

	AtomIteratorMolecule(HaMolecule* new_pMol);
	virtual ~AtomIteratorMolecule();
	
	HaAtom* GetFirstAtom(); //!< Return the first atom of the Molecule (=NULL if no atoms) 
	HaAtom* GetNextAtom();  //!< Return Next atom in the molecule (=NULL if no more atoms)

protected:

	AtomIteratorChain aitr;
	list<HaChain>::iterator ch_itr;

	HaMolecule* pMol;
};

class AtomIteratorMolecule_const: public AtomIterator_const
//! Atom iterator class to browse atoms of the molecule
{
public:

	AtomIteratorMolecule_const(const HaMolecule* new_pMol);
	virtual ~AtomIteratorMolecule_const();

	const HaAtom* GetFirstAtom(); //!< Return the first atom of the Molecule (=NULL if no atoms) 
	const HaAtom* GetNextAtom();  //!< Return Next atom in the molecule (=NULL if no more atoms)
	
protected:

	AtomIteratorChain_const aitr;
	list<HaChain>::const_iterator ch_itr;
	
	const HaMolecule* pMol;
};

class ResidueIteratorMolecule
//! Residue iterator class to browse residues of the molecule
{
public:
	ResidueIteratorMolecule(HaMolecule* new_pmol);
    ResidueIteratorMolecule(const ResidueIteratorMolecule& ritr_ref);
	virtual ~ResidueIteratorMolecule();
	
	HaResidue* GetFirstRes(); //!< Return the first residue of the molecule (=NULL if no residues) 
	HaResidue* GetNextRes();  //!< Return next residue of the molecule (=NULL if no more residues)
	
protected:
	vector<HaResidue*>::iterator res_itr;
	list<HaChain>::iterator ch_itr;
	
	HaMolecule* pmol;
};

class ChainIteratorMolecule
//! Chain iterator class to browse chains of the molecule
{
public:
	ChainIteratorMolecule(HaMolecule* new_pmol);
    ChainIteratorMolecule(const ChainIteratorMolecule& chitr_ref);
	virtual ~ChainIteratorMolecule();
	
	HaChain* GetFirstChain(); //!< Return the first Chain of the Molecule (=NULL if no chains) 
	HaChain* GetNextChain();  //!< Return next chain in the molecule (=NULL if no more chains)
	
protected:
	list<HaChain>::iterator ch_itr;
	
	HaMolecule* pmol;
};

/*=================*/
/*  Database Flags */
/*=================*/

const int SelectFlag   =  0x01;
const int AllAtomFlag  =  0x1c;
const int HelixFlag    =  0x03;


/* Group Flags */
const int CystineFlag    = 0x01;     //!< Disulphide bonded cysteine /
const int StrandFlag     = 0x02;     //!< Strands representation      
const int DashStrandFlag = 0x04;     //!< Dash Strands representation 
const int RibbonFlag     = 0x08;     //!< Solid Ribbon representation 
const int TraceFlag      = 0x10;     //!< Smooth trace representation 
const int CartoonFlag    = 0x20;     //!< Richardson protein cartoon  
const int DotsFlag       = 0x40;     //!< Dotted trace representation 

/* Structure Flags */
const int Helix3Flag  =    0x01;     //!< 3,10-Helix structure       
const int Helix4Flag  =    0x02;     //!< Alpha Helix structure      
const int Helix5Flag  =    0x03;     //!< 5-Helix structure          
const int SheetFlag   =    0x04;     //!< Beta Sheet structure       
const int TurnFlag    =    0x08;     //!< Turn Secondary structure   

#endif /* !HAMOLECULE_H */



