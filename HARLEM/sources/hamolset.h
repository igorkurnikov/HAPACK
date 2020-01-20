 /*!  \file hamolset.h
  
     Class to define set of molecules in HARLEM and accompaniing iterator classes
 
    \author Igor Kurnikov  
    \date 1999-2002
*/
#if !defined(HAMOLSET_H)
#define HAMOLSET_H

#include "hastl.h"
#include "rapidxml.hpp"
#include "command.h"
#include "haatgroup.h"
#include "object3d.h"

class HaMolecule;
class HaAtom;
class HaBond;
class ZMatCrd;
class HaMat_double;
class HaMolView;
class HaQCMod;
class ChemGroup;
class HaCompMod;
class ETCouplMod;
class HaDaltonMod;
class HaGaussMod;
class ElectrostMod;
class PNPMod;
class ElMod;
class pKaCalcMod;
class HaInterMolMod;
class HaMolMechMod;
class HaScatterMod;
class Object3D;
class HaDisplayedSurface;
class StmMod;
class NuclAcidMod;
class HaZindoMod;
class ProtonRedoxMod;
class HaEmpiricalMod;
class MDTrajAnalMod;
class HaMolMembraneMod; //jose
class APBSMod;
class MolEditor;
class HaFlexMod;
class CollectCrdAnalMod;

class AtomIteratorMolSet;
class AtomIteratorMolSet_const;
class BondIteratorMolSet;
class ResidueIteratorMolSet;
class ResidueIteratorMolSet_const;
class ChainIteratorMolSet;
class MolSetEvtHandler; 
class CrdSnapshotIterator;

class FragmentCreatePars;

class MolViewWX;

class ForceFieldType;
namespace mort { class molecule_t; }

//#if !defined(RAPIDXML_HPP_INCLUDED)
//	namespace rapidxml { template<class Ch = char> class xml_node; }
//#endif

//! \brief Class to describe a collection of molecules - Central Class of HARLEM
//! \nosubgrouping
class MolSet : public AtomContainer
{
	friend class MolSetParDlg;
	
public:
//! \name Constructors and Destructors:
//@{
	MolSet();
	virtual ~MolSet();
//@}
//! \name Debug Control 
//@{
	enum DEBUG_FLAG_PARAMS { file_reading_debug = 0x0001 };
	int debug_flag;
//! Print information about the molecule set
	bool Print_info(ostream &sout, const int level); //!< Print molecular set info
//@}

//! \name Axxiliary classes:
//@{
//	typedef AtomIteratorMolSet       AtomIterator;
//  typedef AtomIteratorMolSet_const AtomIterator_const;
	typedef std::vector<HaBond*>    BondIterator;
	typedef ResidueIteratorMolSet   ResidueIterator;
//@}
//! \name RASMOL-type Text Command Processing:
//@{
	int ExecuteCommand(CmdParser& cmd_pr);
	int ExecuteShowCommand(CmdParser& cmd_pr);
//@}
//! \name wxWidgets command processor:
//@{
	MolSetEvtHandler* p_evt_h; //!< Handler for wxWidgets events in MolSet  
//@}

//! \name Superimpose molecules functions:
//@{
	double OverlapMol(AtomGroup& fmolatset, AtomGroup& smolatset); //!< Superimpose two molecules containing atom sets fmolatset and smolatset
    double AlignOverlapMol(AtomGroup& fmolatset, HaMolecule* pMol2, 
		PtrPtrMap* fit=NULL, HaVec_double* p_trans=NULL, HaMat_double* p_rot=NULL); //!< Align Sequences and Overlap two molecules containing atom sets fmolatset and smolatset, return RMS of atoms
//@}
//! \name Input/Output geometry from/to files
//@{
	static AtomLoadOptions* p_load_opt_default; //!< Atom Load Options structure
	static AtomSaveOptions* p_save_opt_default; //!< Atom Save Options structure

	int FetchFile(int format, const char* file_name ); //!< Load file in a given format (FormatPDB,FormatHarlem etc..)
	int LoadHarlemFile (const char* fname, const AtomLoadOptions* popt = NULL);    //!< Load Molecules from file in HARLEM (*.HLM) (OLD or new XML) format
	int LoadAmberPrepFile(const char* fname);        //!< Load Molecule in AMBER PREP (*.IN) format
	int LoadAmberTopFile(const char* fname);         //!< Load Molecule in AMBER TOP (*.TOP;) format and matching *.crd or *.rst format
	int LoadRWFMolecule (const char* fname);         //!< Load structure in binary Gaussian checkpoint (*.rwf, *.chk) format
	int LoadPDBFile(const char* fname, int flag=0);  //!< Load molecule in PDB format (flag = 1 - NMR extension)
	int LoadMol2File(const char* fname);             //!< Load molecule in TRIPOS *.mol format 
	int LoadMDLFile(const char* fname);              //!< Load molecule in TRIPOS *.mdl format 
	int LoadXYZFile(const char* fname, const AtomLoadOptions* popt = NULL);   //!< Load molecule in TINKER XYZ format ( idx, at_nm,x,y,z, ff_idx, at_bond_1, at_bond_2,...)
	int LoadHINFile(const char* fname, const AtomLoadOptions* popt = NULL);   //!< Load molecule in Arbalest HIN format 

	int LoadXYZStream( std::istream& is, const AtomLoadOptions* popt = NULL); //!< Load coordinates from stream
	int LoadHINStream( std::istream& is, const AtomLoadOptions* popt = NULL); //!< Load Molecules from the stream in Arbalest HIN format 
	int LoadXMLStream (std::istream& is, const AtomLoadOptions* popt = NULL); //!< Load from file in HARLEM XML format
	int LoadXMLNode( rapidxml::xml_node<>* node, const AtomLoadOptions* popt = NULL ); //!< Load MolSet data from RAPID XML node
	
	int LoadOldHarlemFile(FILE* fp ); //!< Load file in OLD Harlem format

	int SetCoordFromFile(const char* fname, int iform = FormatGUESS); //!< Set Molecular Set atom coordinates from File (in pdb, xyz or other format) 
	int SetCrdFromArray( const HaVec_double& crd_arr ); //!< Set Molecular Set atom coordinates for coordinate array

	int SavePDBFile(const char* filename );      //!< Save molecules into a file in PDB format
	int SaveHarlemFile(const char* filename, const AtomSaveOptions* popt = NULL );  //!< Save molecules into a file in current HARLEM format (XML) (*.hmlx)
	int SaveOldHarlemFile(const char* filename );   //!< Save molecules into a file in OLD HARLEM format (*.hml)
	int SaveXYZRadFile(const char* filename );   //!< Save file with lines (x,y,z, atom_radius) for MSMS input
	int SaveHINFile(const char* filename );      //!< Save molecules into a file in Arbalest HIN format
	
	int SavePQRFile(const char* filename, bool SaveChainLetter=true);//!< Save molecule into PQR format file
	int SavePQRFreeFile(const char* filename);
	
	TiXmlElement* AddXml(TiXmlElement* parent_element,const char* name = "", int option=0) const; //!< Add Minimal Molecular Set Descripion to XML element
	int SaveXML(FILE* file_out, int option=0) const;              //!< Save Molecular Set Description to a file
	int SaveOldHarlemStream(std::ostream& os); //!< Save MolSet in Old Harlem Format to stream
	int SaveHINToStream(std::ostream& os) const; //!< Save MolSet in Abalest HIN Format to stream
	int SaveXMLToStream(std::ostream& os, const AtomSaveOptions* popt = NULL ) const;  //!< Write MolSet data in XML format to stream
	  
	ZMatCrd* GetZMat( const harlem::HashMap* popt = NULL); //!< Get Z-Matrix for the molecular set   

protected:
	void ProcessPDBAtom(const std::string& line, int heta, IntPtrMap& id_at_map, HaMolecule* pMol, HaChain* &pch_cur, HaResidue* &pres_cur);

//@}

//! \name Molecular Structure Modifications
//@{
public:

	static MolSet* CurMolSet;                            //!< Currently active Molecular Set
	std::vector<HaMolecule*> HostMolecules;                   //!< Molecules in the molecular set  
	std::multimap<int, HaMolecule*, less<int> > serno_mol_map; //!< map of serial numbers of molecules to molecular pointers
	std::multimap<std::string, HaMolecule*, less<std::string> > name_mol_map; //!< map of names of molecules to molecular pointers
	 
	HaMolecule* AddNewMolecule( int mol_ser_no = -1 );  //!< function to add a new molecule to the set 
	HaMolecule* GetFirstMolecule();                  //!< Get First Molecule in the set
	HaMolecule* GetMolByIdx(int imol);               //!< Get Molecule  by index (0-based)
	HaMolecule* GetMolByIdx0(int imol);              //!< Get Molecule by index (0-based)
	HaMolecule* GetMolByIdx1(int imol);              //!< Get Molecule by index (1-based)

	HaMolecule* GetMolByName(const char* mol_name);  //!< Get first molecule with a given name
	const HaMolecule* GetMolByName(const char* mol_name) const;  //!< Get first molecule with a given name
	HaMolecule* GetMolByRef(const char* mol_ref);  //!< Get molecule by reference ( NAME + IDX(1-based) or IDX )
	const HaMolecule* GetMolByRef(const char* mol_ref) const;  //!< Get molecule by reference ( NAME + IDX(1-based) or IDX )

	void DeleteAll();                    //!< Delete all the contents of the molecular set 
	bool DeleteMol(HaMolecule* pMol);    //!< Delete molecule
	bool DeleteAtomWithRef(const char* atref); //!< Delete Atom identified by reference
	bool DeleteAtom(HaAtom* aptr);       //!< Delete Atom
	bool DeleteAtoms(AtomContainer& atset);    //!< Delete Atoms

	void OnAtomSeqChange(); //!<  Invalidate atom indexes when atoms are added, deleted or rearranged
	void OnChangePeriodicity(); //!< Notify Views and modules that periodical box was create or deleted

	int RenumberSelectedRes(int start_num = 1); //!< Renumber Selected Residues (should be in one chain) 
//@}
//! \name Iterations over atoms: implementations of AtomContainer and PointContainer functions:
//@{
	virtual PointIterator*       GetPointIteratorPtr();
	virtual PointIterator_const* GetPointIteratorPtr() const;
	virtual int GetNumPt() const; 
	virtual AtomIterator*       GetAtomIteratorPtr(); 
	virtual AtomIterator_const* GetAtomIteratorPtr() const; 
	virtual int IsMember(const HaAtom* aptr) const;
//@}
//! \name Molecular Structure Info
//@{
	int GetNMol() const;                 //!< Get the number of Molecules in the molecular set
	int GetNRes() const;                 //!< Get the number of Residues in the molecular set
	int GetNChains() const;              //!< Get the number of chains in the molecular set
	virtual int GetNAtoms() const;               //!< Get the total number of atoms in the molecule set 
	int GetNDumAtoms() const;            //!< Get the number of dummy atoms (elemno = DUMMY_ELEM )
	int GetNBonds() const;               //!< Get the total number of bonds in the molecule set
	int GetNHBonds() const;              //!< Get the number of H-bonds in the Molecular Set
	int GetNSSBonds() const;             //!< Get the number of SS-bonds in the Molecular Set
	int GetNBackbBonds() const;          //!< Get the number of backbone bonds in the Molecular Set
	
	double FindClosestContact(HaAtom* atc1,HaAtom* atc2); //!< Find closest Atom contact in the molecular Set

	HaResidue* GetResByRef(const char* res_ref);  //!< get a residue by reference
	HaAtom* GetAtomByRef(const std::string& at_ref);     //!< get a single atom by reference
	bool GetAtomsByRef(const char* at_ref, AtomGroup& at_set); //!< get a group of atoms by reference
	
	HaAtom* GetAtomBySeqNum(int seq_num); //!< Get Atom By Sequence Number (0-based) 
	int GetSeqNumForAtom( HaAtom* aptr);  //!< Get Sequence Number (0-based) for Atom

	void InitAtomIdx();  //!< Init Atom Array and Atom Sequence Number Map

	AtomIntMap  GetAtomSeqNumMap();        //!< Get a map of atom pointers to their current sequence numbers 
	CAtomIntMap GetAtomSeqNumMap() const;  //!< Get a map of atom pointers to their current sequence numbers 
	
//@}
//! \name Selection functions:
//@{
    void SelectAtomsAll();    //!< Select all atoms in the set
	void SelectAtoms(AtomContainer* atom_coll); //!< Select Atoms of the Atom collection and Unselect other atoms 
    void SelectAtomsMask( int mask );       //!< Select Atoms with flags satisfying mask
	void SelectAtomsExprObj( AtomExpr* expr ); //!< Select Atoms satisfying expression
	void SelectAtomsExpr( const char* expr_str); //!< Select Atoms using string with RASMOL-like selection expression ("select" keyword is added in the function)
	void SelectAtomsInBoundaryBox();   //!< Select Atoms within Boundary Box
	void UnSelectAtomsAll();  //!< Unselest all atoms in the set
	void RevertAtomSelection();  //!< Select All Unselected atoms and UnSelect all Selected Atoms
	void DisplaySelectCount();         //!< Print the number of selected atoms
//@}
//! \name Covalent Bonds, H-Bonds, SS-bonds, Backbone Manipulation:
//@{
	int AreHBonded(HaAtom* src, HaAtom* dst) const; //!< check if atoms are H-Bonded (src - H-donor and dst H-acceptor 
    HaHBond* AddHBond(HaAtom* src, HaAtom* dst);  //!< Create a hydrogen bond 

    void CreateHydrogenBond(HaAtom* src, HaAtom* dst, int energy, int offset);

	HaBond*  AddBond( HaAtom* src, HaAtom* dst ); //!< Create a covalent bond
    bool DeleteBond(HaAtom* src, HaAtom* dst);    //!< Delete Covalent bond between atoms

	void ClearBackbone();  //!< Clear Backbone bonds

	set<HaHBond, less<HaHBond> > HBonds;        //!< Hydrogen Bonds of the molecular set
	std::vector<HaBond*>  Bonds;    //!< Valence bonds of the molecular set
	std::vector<HaBond*>  BackboneBonds; //!< Backbone Bonds in the molecular set 

	bool SSBonds_found;  //!< flag to show if SS-bonds have been found 	
	bool HBonds_found;   //!< flag to show if H-bonds  have been found
	bool to_find_backb;  //!< Flag to find Backbone bonds before using them

protected:
	ZMatCrd*   p_zmat;        //!< ZMatrix of the molecular set
	MolEditor* p_mol_editor;  //!< Module to edit molecular structure

    AtomIntMap at_seq_num_map; //< atom pointer to index (0-based) map - for fast finding atom indexes in the Molecular Set
	AtomGroup  at_array;       //< atom array for access to

public:
//@}

//! \name Periodical Boundary Conditions:
//@{
	PeriodicUnitInfo* per_bc; //!< Periodic Unit Information for the Molecular Set

	int WrapToUnitCell(); //!< Wrap coodinates of the molecular set to Periodical Unit Cell
	int WrapAndCenter( const std::string& grp_name, const Vec3D& cnt_crd ); ////!< Wrap Coorinates with AtomGroup grp_name centered at coordinates cnt_crd
//@}
//! \name Names of Molecules and Molecular Set
//@{
	void SetName(const char* new_name); //!< Set the name of the molecular set
	const char* GetName() const;              //!< Get the name of the molecular set
	std::string GetUniqueMolName(const std::string& suggest_name); //!< Modify a suggested name of the molecule to avoid a conflict with the existing names of the molecules in the set

	std::string name_mset; //!< the name of the molecular set
//@}

//! \name Get Iterators on Molecular Structure Elements:
//@{
	ResidueIteratorMolSet GetResidueIterator(); // Get Residue Iterator 
//@}

//! \name Chemical(Functional) Groups:
//@{
   	int GetNChemGroups() const;    //!< Get total number of chemical groups
	ChemGroup* AddBlankChemGroup(const std::string& gid = "" );   //!< Add an empty chemical groups
	bool DeleteChemGroup(const std::string& gid );                //!< Delete Chemical group
	bool DeleteChemGroupPtr( ChemGroup* grp_ptr );   //!< Delete Chemical Group if pointer is known
	bool SetChemGrpSelected(const std::string& gid);            //!< Set Selected Atoms to a chemical group
	ChemGroup& GetChemGroupByIdx(int index) ;                  //!< Get Chemical Group by index
	ChemGroup* GetChemGroupByID(const std::string & gid);         //!< Get Chemical Group by id
//	const ChemGroup& GetChemGroupByIdx(int index) const ;      //!< Get Chemical Group by id
	ChemGroup* GetChemGroupByAtom(const HaAtom* aptr);         //!< Get Chemical Group containg atom
	
	bool SetStdChemGroups();  //!< Set Standard Chemical Groups for the molecular set
	void RenumberGrp();       //!< Reset Group number to consequitive starting from 10 

	bool CheckUniqChemGrpID(const std::string& gid);
	std::string GetUniqChemGrpID(int buf_reg_flag);

	list<ChemGroup> ChemGroups; //!< Chemical Functional Groups
	VecPtr chemg_idx;
	typedef list<ChemGroup> ChemGroupsType;

//@}

//! \name Named Atom Groups:
//@{
	AtomGroup* AddAtomGroup( const char* id = ""); //!< Add new Atom Group with id
    AtomGroup* GetAtomGroupByID( const char* id); //!< Get Atom Group By id

	AtomGroup* SetAtomGroupFromSelection( const char* id); //!< Set Atom Group with id by selected atoms of the molset (create new group if needed)

	list<AtomGroup> NamedAtomGroups; //!< Named Atom Groups
	typedef list<AtomGroup> NamedAtomGroupsType;

	bool DeleteAtomGroup(const char* id );         //!< Delete Atom Group with id
	bool DeleteAtomGroupPtr( AtomGroup* atgrp_ptr );     //!< Delete Atom Group by pointer
    int CreateAxxMol(const char* mol_name, const char* id); //!< create axxiliary molecule from the group of atoms to set external charges or force centers
//@}

//! \name Atom Coordinate Snapshots
//@{
	CrdSnapshotIterator GetCrdSnapshots(); //!< Get iterator over coordinate snapshots of the molecular set
	vector<CrdSnapshot*> crd_snapshots; //!< Coordinate Snapshots
	void DeleteCrdSnapshots(); //!< Delete All Coordinate Snapshots
	CrdSnapshot* AddCrdSnapshot(const std::string& snap_name_new = ""); //!< add atom coordinate snapshot for given Molecular Set 
	CrdSnapshot* AddCrdSnapshotForGroup(const std::string& grp_name, const std::string& snap_name_new = ""); //!< add atom coordinate snapshot for atom group of the given Molecular Set 
	CrdSnapshot* GetCrdSnapshotByName(const char* snp_name, bool create = false ); //!< Get Atom Coordinate Array By Name
	int DeleteCrdSnapshot( CrdSnapshot* psnap ); //!< Delete Snapshot 

	int SetCrdFromSnapshot( CrdSnapshot* psnap ); //!< Set Coordinates of the Molecular Set from Snapshot
	int SetCrdFromSnapshot( const std::string& snap_name ); //!< Set Coordinates of the Molecular Set from Snapshot with name snap name

	int SaveCrdSnapshots(const std::string& fname, const harlem::HashMap* popt = NULL ) const; //!< Save Crd snapshots to file in XML format
	int SaveCrdSnapshots(std::ostream& os, const harlem::HashMap* popt = NULL ) const;   //!< Save Crd snapshots to stream in XML format
	int LoadCrdSnapshots(const std::string& fname, const harlem::HashMap* popt = NULL ); //!< Load Crd snapshots from file in XML format

//@}  

//! \name Secondary Structure Elements:
//@{
	int DescribeSecStruct(); //!< Print the description of secondary structure elements
	int PrintHBonds();       //!< Print list of H-Bonds in the system
//@}

//! \name molecular Surfaces & Volume manipulations
//@{
	HaDisplayedSurface* GetMolSurface(int create_flag = FALSE);  //!< Get Molecular Surface (Object3D with name "VdW_surface" )
	HaDisplayedSurface* CalcMolSurface(int surf_type = 1); //!< Calculate Molecular Surface
	HaDisplayedSurface* CalcMolSurfDens(); //!< Calculate Molecular Surface computing isosurface of density
	int CalcSolventAccessArea(); //!< Calculate solvent access area of selected atoms using GeoBall library
	int CreateExcludedVolumeMol(); //< Create EXCLUDED_VOLUME molecule with excluded volume Grid Points for selected atoms
	int SaveCrdExclVolArb(); //!< Save EXCLUDED_VOLUME molecule coordinates In Arbalest Exclude Volume Format
//@}

public:
//! \name Atomic Parameters (charges. force-field symbols etc):
//@{
	double  CalculatePotential( double x, double y, double z ); //!< calculate electrostatic potential 
	                                                            //!< in a point from atomic point charges and eps=10.0
	bool CalcDipole();    //!< Calculate and print charge and dipole of selected atoms of mol set 
	bool SetVdwRadii();   //!< Set atomic radii to Standard VdW  Radii
	bool SetParseRadii(); //!< Set atomic radii to PARSE radii set
	bool SetHPPRadii();   //!< Set atomic radii used on H++ server
	
    vector<AtomDoubleMap> ChargeMaps; //!< The vector of atom-charge maps that 

	AtomDoubleMap* GetChargeMapByName(const char* map_name); //!< Get a atom-charge map by the name
	AtomDoubleMap* CreateChargeMap(const char* map_name);    //!< Create an atom-charge map with a given name 
	int SetChargeMapByCurrentCharges(const char* map_name);        //!< Set a atom-charge name by the current atomic chagres of the molecular set 
	int SetChargesFromChargeMap(AtomDoubleMap* charge_map);  //!< Set atomic charges from the atom-charge map 
//@}
//! \name Fragmentation:
//@{	
	MolSet* parent_mset;     //!< The pointer to the parent Molecular Set ( != NULL if molset is a fragment)
	MolSet* CreateFragmentFromSelection(const char* frag_name, FragmentCreatePars* params = NULL);
	
	vector<MolSet*> Fragments;   //!< Vector of fragments
	PtrPtrMap frag_atom_maps;      //!< map of fragments to mappings of atoms of fragments to atoms of the molecular set
	
	int AssociateFragment(MolSet* frag);     //!< Add an existing molecular set as a fragment
	int ReleaseFragment(MolSet* frag);  //!< Detach the fragment from the molecular set without deleting it 
	int DeleteFragment(MolSet* frag);   //!< Delete fragment 
	int ReleaseAllFragments();            //!< Detach all fragment from the molecular set without deleting them
	int DeleteAllFragments();             //!< Delete all fragments of the molecular set
	int BuildFragmentAtomMap(MolSet* frag, AtomAtomMap& frag_atom_map); //<! build the map of correspondence between the atoms of the fragment and those of the parent molecule
    int SelectAtomsMatchingFragment(MolSet* frag); //!< Select Atoms matching atoms of the fragment
	int IsFragment(const MolSet* pmset); //!< Check if the Molecular Set is the fragment of the given Molecular Set
	int FragmentIdx(const MolSet* pmset); //!< Return index of the fragment in Fragments & frag_atom_map arrays, (-1) if pmset is not a fragment of this Molecular Set
	int SyncFragmentCoord(MolSet* frag);  //!< Synchronize Cooordinates of the fragment with current coordinates of the Molecular Set 
//@}
//! \name Computational Modules:
//@{
	HaCompMod*       GetCompModule( int mtype, bool create_module = false); //!< Get a computational module of a given type associated with the molecular set 
	const HaCompMod* GetCompModule( int mtype ) const; //!< Get a computational module of a given type associated with the molecular set (const version ) 

	HaQCMod*       GetQCMod( const bool create_module = false);       //!< Get Quantum Chemical module associated with the molecular set 
	const HaQCMod* GetQCMod() const;   //!< Get Quantum Chemical module associated with the molecular set (const version)   
	ETCouplMod*       GetETCouplMod( const bool create_module = false);  //!< Get Electron Transfer computational module associated with the molecular set 	
	HaGaussMod*    GetGaussMod( const bool create_module = false);    //!< Get Gaussian qchem package interaction module associated with the molecular set 
	HaDaltonMod*      GetDaltonMod( const bool create_module = false);   //!< Get Dalton qchem package  interaction module associated with the molecular set 
	ElectrostMod*   GetElectrostMod( const bool create_module = false);   //!< Get Continuum Electrostatic computational module associated with the molecular set 
#ifdef ELMOD_COMPILE
	ElMod*   GetElMod( const bool create_module = false);   //!< Get Continuum Electrostatic computational module associated with the molecular set 
#endif
	pKaCalcMod* GetpKaCalcMod( const bool create_module = false);
	PNPMod*   GetPNPMod( const bool create_module = false);   //!< Get Continuum Electrostatic computational module associated with the molecular set (PNP)
	APBSMod*       GetAPBSMod( const bool create_module = false);
	HaInterMolMod* GetInterMolMod( const bool create_module = false); //!< Get Intermolecular interaction module associated with the molecular set 
	HaMolMechMod*  GetMolMechMod( const bool create_module = false);  //!< Get Molecular Mechanics module associated with the molecular set 
	const HaMolMechMod* GetMolMechMod() const;  //!< Get Molecular Mechanics module associated with the molecular set ( const version ) 
	MDTrajAnalMod* GetTrajAnalMod( bool create_module = false); //!< Get MD trajectory analysis module
	HaScatterMod*  GetScatterMod( const bool create_module = false);  //!< Get the Electron Scattering module associated with the molecular set 
	StmMod*          GetSTMMod( const bool create_module = false);            //!< Get the Electron Scattering module associated with the molecular set 
	NuclAcidMod*     GetNuclAcidMod( const bool create_module = false); //!< Get Nucleic Acid modeling module 
	HaZindoMod*      GetZindoMod( const bool create_module = false);    //!< Get ZINDO interface module associated with the molecular set 
	ProtonRedoxMod*  GetProtonRedoxMod( const bool create_module = false);    //!<  Get Chemical Transformations Module
    HaEmpiricalMod* GetEmpiricalMod( const bool create_module = false);    //!<  Get  Empirical calculations module
    HaMolMembraneMod* GetMolMembraneMod(const bool create_module = false); //!< Get Membrane model calculations module jose
	HaFlexMod*        GetFlexMod(const bool create_module = false);        //!< Get Molecular Flexibility Module
	MolEditor*        GetMolEditor(const bool create_module = false);      //!< Get Molecular Editor Module
	CollectCrdAnalMod*   GetCollectCrdAnalMod(const bool create_module = false); //!< Get Cluster Analysis Module

	vector<HaCompMod*> CompModules; //!< List of Computational modulles associated with the molecular set
//@}
public:
//! \name Molecular Display functions:
//@{
	virtual int AnnounceGeomChange(); //!< Call to update views to reflect changes in the molecular geometry
	void RefreshAllViews(long lHint = 0L); //!< Refresh Molecular Views associated with molecular set
	HaMolView* GetActiveMolView();   //!< Get Active HaMolView object associated with molecular set 
	const HaMolView* GetActiveMolView() const;   //!< Get Active HaMolView object associated with molecular set (const version)

	HaMolView* mset_pview;  //!< A pointer to the Molecular View object associated with the molecular set 
	MolViewWX* canvas_wx;   //!< A pointer to the wxWindows Molecular Canvas View associated with the molecular set
	bool AddObject3D(Object3D* new_view_object);   //!< Add a new 3D object to the list
	bool DeleteObject3D(Object3D* pobj);           //!< Delete 3D object from the list using the pointer 
	bool DeleteObject3D(const std::string& obj_name);  //!< Delete 3D object from the list using its name

	list<Object3D*> ViewObjects; //!< the list of 3D objects

	void ClearPickedAtoms(); //!< Clear Picked Atoms Set
	AtomGroup picked_atoms;

	StrVec info_str; //!< vector of string to display on the screen
//@}

//! \name GLEAP and MORT library interactions
//@{
	int SetMortMol(mort::molecule_t& mort_mol, const ForceFieldType& ff_type ); //!< Set Content of mort::molecule_t using info of the Molecular Set 
	int SavePDBMort(const char* fname); //!< Save PDB file using MORT Library
	int SaveSDFMort(const char* fname); //!< Save SDF file using MORT Library
//@}

};

typedef vector<HaMolecule*> MoleculesType;

class AtomIteratorMolecule;
class AtomIteratorMolecule_const;


class AtomIteratorMolSet: public AtomIterator
//! Atom iterator class to browse atoms of the molecular set
{
public:
	AtomIteratorMolSet(MolSet* new_pmset);
	AtomIteratorMolSet(const AtomIteratorMolSet& ref);
	virtual ~AtomIteratorMolSet();
	
// Implementation of AtomIterator functions 

	virtual PointIterator* clone() const; //!<  Create a copy of the iterator with the same state (from PointIterator )
	virtual HaAtom* GetFirstAtom(); //!< Return the first atom of the Molecular Set (=NULL if no atoms) 
	virtual HaAtom* GetNextAtom();  //!< Return Next atom in the sequence (=NULL if no more atoms)
	
	HaAtom* operator *() { return *aitr_res; }

#if defined(SWIG) 
%exception {
try {
	$action
	} catch(std::out_of_range) {
		PyErr_SetString(PyExc_StopIteration,"End of Atoms in the MolSet");
		return NULL;
	}
}
#endif
	HaAtom* next(); //!<  Return next atom in the sequence (first on the first call). Throw std::out_of_range() if no more atoms (Python compatibility)
#if defined(SWIG)
%exception;
#endif
	AtomIteratorMolSet __iter__() const; //!< Get a copy of the iterator ( for Python compatibility )

protected:
	std::vector<HaAtom*>::iterator aitr_res;
	std::vector<HaAtom*>::iterator aitr_res_end;
	ResidueIteratorMolSet* pritr;
	bool first_called; //!< flag to indicate that iterator was already called at least one time
};


class AtomIteratorMolSet_const : public AtomIterator_const
//! Const Atom iterator class to browse atoms of the molecular set
{
public:
	AtomIteratorMolSet_const(const MolSet* new_pmset);
	AtomIteratorMolSet_const(const AtomIteratorMolSet_const& ref);
	virtual ~AtomIteratorMolSet_const();

	virtual PointIterator_const* clone() const; //!<  Create a copy of the iterator with the same state (from PointIterator_const )
	virtual const HaAtom* GetFirstAtom(); //!< Return the first atom of the Molecular Set (=NULL if no atoms) 
	virtual const HaAtom* GetNextAtom();  //!< Return Next atom in the sequence (=NULL if no more atoms)
	
protected:
    vector<HaAtom*>::const_iterator aitr_res;
	vector<HaAtom*>::const_iterator aitr_res_end;
	ResidueIteratorMolSet_const* pritr;
};


class ResidueIteratorMolSet
//! Residue iterator class to browse residues of the molecular set
{
public:
	ResidueIteratorMolSet(MolSet* new_pmset);
	ResidueIteratorMolSet(const ResidueIteratorMolSet& ritr_ref);
	virtual ~ResidueIteratorMolSet();
	
//	ResidueIteratorMolSet* clone(); //!< Create a copy of the iterator with the same state
	HaResidue* GetFirstRes(); //!< Return the first residue of the molecular Set (=NULL if no atoms) 
	HaResidue* GetNextRes();  //!< Return next residue in the molecular set (=NULL if no more atoms)
	HaResidue* GetCurrRes();  //!< Get Residue Pointer corresponding to Current Iterator State

#if defined(SWIG) 
%exception{
	try {
		$action
		} catch (std::out_of_range) {
	  PyErr_SetString(PyExc_StopIteration,"End of Atoms in the Atom Collection");
	  return NULL;
      }
}
#endif
	HaResidue* next(); //!<  Return next residue in the sequence (first on the first call) throw std::out_of_range() if there are no more residues ( Python compatibility )
	HaResidue* __next__();

#if defined(SWIG)
%exception;
#endif
	ResidueIteratorMolSet __iter__() const; //!< Get a copy of the iterator ( Python compatibility )

protected:
	vector<HaResidue*>::iterator res_itr;
	list<HaChain>::iterator ch_itr;
	MoleculesType::iterator mol_itr;

	MoleculesType::iterator mol_itr_begin;
	MoleculesType::iterator mol_itr_end;

	int first_called;
};

class ResidueIteratorMolSet_const
//! Const Residue iterator class to browse residues of the molecular set
{
public:
	ResidueIteratorMolSet_const(const MolSet* new_pmset);
	virtual ~ResidueIteratorMolSet_const();
	
	const HaResidue* GetFirstRes(); //!< Return the first residue of the Molecular Set (=NULL if no atoms) 
	const HaResidue* GetNextRes();  //!< Return the next residue in the sequence (=NULL if no more atoms)
	
protected:
	vector<HaResidue*>::const_iterator res_itr;
	list<HaChain>::const_iterator ch_itr;
	MoleculesType::const_iterator mol_itr;

	MoleculesType::const_iterator mol_itr_begin;
	MoleculesType::const_iterator mol_itr_end;
};

class ChainIteratorMolSet
//! Chain iterator class to browse chains of the molecular set
{
public:
	ChainIteratorMolSet(MolSet* new_pmset);
	virtual ~ChainIteratorMolSet();
	
	HaChain* GetFirstChain(); //!< Return the first chain of the molecular Set (=NULL if no atoms) 
	HaChain* GetNextChain();  //!< Return next chain in the molecular set (=NULL if no more atoms)
	
protected:
	list<HaChain>::iterator ch_itr;
	MoleculesType::iterator mol_itr;
	
	MolSet* pmset;
};

class ChemGroupIterator
//! Iterator class to browse Chemical Groups of the molecular set
{
public:
	ChemGroupIterator(MolSet* new_pmset);
	virtual ~ChemGroupIterator();
	
	ChemGroup* GetFirst(); //!< Return the first Chemical Group of the Molecular Set (=NULL if no Chemical Groups) 
	ChemGroup* GetNext();  //!< Return next Chemical Group in the Molecular Set (=NULL if no more Chemical Groups)

protected:
    list<ChemGroup>::iterator CurGroupItr;	
	MolSet* pmset;
};

class AtomGroupIteratorMolSet
//! Iterator class to browse Named Atom Groups of the molecular set
{
public:
	AtomGroupIteratorMolSet(MolSet* new_pmset);
	virtual ~AtomGroupIteratorMolSet();
	
	AtomGroup* GetFirst(); //!< Return the first named Atom Group of the Molecular Set (=NULL if no named Atom Groups) 
	AtomGroup* GetNext();  //!< Return next named Atom Group (=NULL if no more Atom Groups)
	
protected:
	list<AtomGroup>::iterator CurListItr;	
	MolSet* pmset;
};

class AtomGroupIteratorMolSet_const
//! Iterator class to browse Named Atom Groups of the molecular set (const version)
 {
public:
	AtomGroupIteratorMolSet_const(const MolSet* new_pmset);
	virtual ~AtomGroupIteratorMolSet_const();
	
	const AtomGroup* GetFirst(); //!< Return the first named Atom Group of the Molecular Set (=NULL if no named Atom Groups) 
	const AtomGroup* GetNext() ; //!< Return next named Atom Group (=NULL if no more Atom Groups)
	
protected:
	list<AtomGroup>::const_iterator CurListItr;	
	const MolSet* pmset;
};


class BondIteratorMolSet
//! Bond iterator class to browse bonds of the molecular set
{
public:
	BondIteratorMolSet(MolSet* new_pmset);
	virtual ~BondIteratorMolSet();
	
	HaBond* GetFirstBond(); //!< Return the first bond of the Molecular Set (=NULL if no bonds) 
	HaBond* GetNextBond();  //!< Return Next bond in the sequence (=NULL if no more bonds)
	
protected:
	vector<HaBond*>::iterator bitrm;
	MoleculesType::iterator mol_itr;
	
	MolSet* pmset;
};

class HBondIteratorMolSet
//! Bond iterator class to browse Hydrogen Bonds of the molecular set
{
public:
	HBondIteratorMolSet(MolSet* new_pmset);
	virtual ~HBondIteratorMolSet();
	
	HaHBond* GetFirstBond(); //!< Return the first Hydrogen bond or SS bond of the Molecular Set (=NULL if no H-bonds(SS bonds)) 
	HaHBond* GetNextBond();  //!< Return next Hydrogen bond or SS bond of the Molecular (=NULL if no more H-bonds(SS bonds))
	
protected:
	set< HaHBond,less<HaHBond> >::iterator bitrm;
	MoleculesType::iterator mol_itr;
		
	MolSet* pmset;
};


class FragmentCreatePars
//! axxiliary class to specify parameters for fragment creation
{
public:
	FragmentCreatePars();
	int add_hydr;
};

class PyAccMolSetProp
//! class for python accellerated access to molset properties
{
public:
	PyAccMolSetProp(MolSet* new_pmset);
	~PyAccMolSetProp();

	MolSet* pmset;

	std::vector<int>* GetAtomsSerNoAsVec();
	std::vector<int>* GetResidueSerNoAsVec();
	std::vector<double>* GetAtomsChargeAsVec();
	std::vector<double>* GetAtomsRadiusAsVec();
	std::vector<std::string>* GetAtomsNameAsVec();
	std::vector<std::string>* GetResidueNameAsVec();
	std::vector<double>* GetAtomsCoorXAsVec();
	std::vector<double>* GetAtomsCoorYAsVec();
	std::vector<double>* GetAtomsCoorZAsVec();
	std::vector<double>* GetAtomsIonExcludedRadiusAsVec(double Rion);

	void WriteAtomParamFileForPNP(const char *filename,\
	std::vector<int>* ResidueSerNo, std::vector<std::string>* ResidueName, std::vector<int>* AtomsSerNo, std::vector<std::string>* AtomsName, \
	std::vector<double>* AtomsCoorX, std::vector<double>* AtomsCoorY, std::vector<double>* AtomsCoorZ,\
	std::vector<double>* AtomsCharge, std::vector<double>* AtomsRadius,  std::vector<double>* AtomsIER1, std::vector<double>* AtomsIER2,\
	std::vector<float>* SR_A_K, std::vector<float>* SR_N_K, std::vector<float>* SR_A_Cl, std::vector<float>* SR_N_Cl);
};

extern "C" {
#if defined(HAMOLSET_CPP)
	MolSet* GetCurMolSet() { return MolSet::CurMolSet; }
	void SetCurMolSet(MolSet* pmset);
#else
	extern MolSet* GetCurMolSet();
	extern void SetCurMolSet(MolSet* pmset);
#endif
}

#endif // end !defined(HAMOLSET_H)


