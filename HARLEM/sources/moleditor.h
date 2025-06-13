/*! \file moleditor.h

     Classes to edit molecular structures 

    \author Igor Kurnikov
    \date 2008-
*/
#ifndef MOLEDITOR_H_HARLEM
#define MOLEDITOR_H_HARLEM

class MolEditor
//!< Class to edit molecules  
{
public:
	MolEditor();
	virtual ~MolEditor();

	void Init(); //!< Init with standard parameters

	vector<HaAtom> FindMissingHydrogens( HaAtom* aptr );

	int DeleteExtraAtoms(MolSet* pmset);   //!< Delete atoms not found in residues templates
	int AddMissingAtoms(MolSet* pmset);    //!< Add Missing atoms in the structure using residue templates
	int AddHydrogens(MolSet* pmset);       //!< Add hydrogens
	int AddPolarHydrogens(MolSet* pmset);  //!< Add Polar hydrogens involved in hydrogen-bonding
	int AddHydrogensHybrid(MolSet* pmset); //!< Add Hydrogens using atom hybridization information
	int SetAtomElemFromTempl(MolSet* pmset){return 0;} //!< Set Atom Elements from Residue Templates
	int SetHBondDonAccStatus(AtomContainer* p_at_coll); //!< Set Hydrogen Bond donor Acceptor status for atoms in a collection
	int SetStdAtomicParams(MolSet* pmset, int at_params_type); //!< Set atomic parameters according to flag at_params_type 
	int SetFormalAtChrgFromTempl(MolSet* pmset); //!< Set Formal Atomic Charges From Residue Templates
	static int ClearAtomFFParams(MolSet* pmset); //!< Clear Atomic Force Field Parameters 
	int FixBondsUsingTempl(MolSet* pmset); //!< Fix Bonds of the selected residues according to residue templates
	int OrderAtomsInRes(MolSet* pmset);    //!< order atoms in residues according to residue templates
	int RenameAtomsFlexToAmber(HaResidue* pres); //!< Rename Atoms of Nucleic Acid Residue from FLEX to AMBER notation 
	static int FixFlexDNA(HaMolecule* pMol); //!< Fix Structure of DNA molecule build in FLEX module 
	int RenameAtomsToAmber(MolSet* pmset); //!< Rename Atoms of Residues according to AMBER databases   
	int RenameAtomsToGromacs(MolSet* pmset); //!< Rename Atoms of Residues according to GROMACS databases  
	int ConvertWaterArrowVB(MolSet* pmset); //!< Convert Water Residues to Separate Molecules, make H-H virtual bond and Rename Water atoms
	int ConvertWaterFastAmber(MolSet* pmset); //!< Coonvert Water Residues to Amber Fast Water ( make bonds between H atoms , rename water residues to WAT )
	
	int CreateCovBonds(AtomContainer* at_col); //!< Compute Covalent Bonds based on atom-atom distances and atomic VdW radii  
	static int BondIfClose(HaAtom* sptr, HaAtom* dptr); //!< Bond two atoms if they they close enough (relative to their VdW radii)

	int SetBondDist(HaAtom* aptr1, HaAtom* aptr2, double new_dist);             //!< Set bond distance between atoms to a new value 
    int SetAngle(HaAtom* aptr1, HaAtom* aptr2, HaAtom* aptr3, double ang_new);  //!< Set Valence angle (in rad)
    int SetTorsion(HaAtom* aptr1, HaAtom* aptr2, HaAtom* aptr3, HaAtom* aptr4, double tors_new); //!< Set Dihedral angle (in rad)

	int FindHBondsAtomCollection( AtomContainer* p_at_coll, vector<HaHBond>& hbonds); //!< Compute H-Bonds for the collection of atoms
    bool CalcHBonds(MolSet* pmset, bool recalc=false);  //!< Compute Hydrogen Bonds for the molecular Set
	int IsValidHBond(HaHBond* p_hb); //!< Check if hydrogen bond satisfy current criteria 
	void CalcHydrogenBonds(MolSet* pmset);            //!< Compute Hydrogen Bonds for the molecular Set - Old,not used, to incorporate
	void FindDisulphideBridges(MolSet* pmset);             //!< Find Disulphide bridges
	int CalcNucleicHBonds(HaChain* chn1 );   //!< Find backbone hydrogen bonds for a nucleic acid chain
	int CalcProteinHBonds(HaChain* chn1 );   //!< Find backbone hydrogen bonds for a protein chain
	int FindBackbone (MolSet* pmset);  //!< Compute Backbone Bonds for the molecular Set
	int UpdateBackBone(MolSet* pmset); //!< Recompute Backbone Bonds for the molecular Set if to_find_backb flag is set

    int DetermineSecStructure(HaMolecule* pmol, int flag); //!< Compute secondary structure of the molecule
	int FindAlphaHelix( HaMolecule* pmol, int pitch, int flag );
	int FindTurnStructure( HaMolecule* pmol);
	int FindBetaTurns( HaMolecule* pmol);
    int FindBetaSheets( HaMolecule* pmol);

	void SetAlphaHelix(MolSet* pmset); //!< Set selected part of the protein to an ideal alpha helix conformation

	double max_bond_length;      //!< Maximal covalent bond length (in Bohrs) - used in covalent bond finding 
	double max_hbond;            //!< max D/A dist to search for H-Bonds (Bohrs)
	double max_da_dist_no_acc;   //!< max D/A dist for H -Bond including N,O acceptor atoms (Bohrs)
	double max_da_dist_s_acc;    //!< max D/A dist for H -bond including S atoms  (Bohrs)
	double max_ha_dist_no_acc;   //!< max H/A dist for H -Bond including N,O acceptor atoms (Bohrs)
	double max_ha_dist_s_acc;    //!< max H/A dist for H -Bond including S acceptor atoms  (Bohrs)
	double max_hda_angle;        //!< max H - donor - Acceptor angle in radians (Bohrs)
	bool m_calc_s_hbonds_flag;   //!< flag to compute H-Bonds involving sulfur

	int CenterAtOrigin(AtomContainer* at_cont); //!< Center atom container at the coordinate origin
	int CenterAtOriginWithRad(AtomContainer* at_cont); //!< Center atom container at the coordinate origin taking into account atom VdW radii 
	
	int Solvate(MolSet* pmset); //!< Solvate Molecular Set 
	int CenterSoluteInSolvent(MolSet* pmset); //!< Center solute in the center of periodic box
	int CenterMolInPBox(MolSet* pmset);  //!< Center Molecule in the center of periodic box
	void AddIons(MolSet* pmset, int n_na, int n_cl); //!< Add Na+ and CL- atoms 
	int ReplicatePeriodBox(MolSet* pmset, int nx, int ny, int nz); //!< Replicate Molecular Set Along X,Y,Z axes using box periodical boundary info
    int WrapToUnitCell(AtomContainer* at_cont, PeriodicUnitInfo* per_info); //!< Wrap Atoms in the periodical system to the Unit Cell  
    int DeleteOverlapMols(MolSet* pmset, AtomGroup& at_coll); //!< Remove all molecules (connected groups of atoms) that overlap at_coll 
	int SplitToMolecules(AtomContainer* p_at_coll, vector<AtomGroup>& mols); //!< Split Atom Collection into groups of bonded atoms (molecules) 

	

	int MergeMolecules(HaMolecule* pMol1, HaMolecule* pMol2); //!< Merge molecule 2 into molecule 1

	HaMolecule* CreateTransAlk(MolSet* pmset, const int nunit, const std::string& name="ALK");   //!< Create trans-alkane with the length of nunit   
	HaMolecule* CreateSurf(MolSet* pmset, const int num_layers, const std::string& name="GOLD"); //!< Create Several layers of the metal surface
	bool Create2DMolArray(MolSet* pmset, HaMolecule* pMol_ref, 
						   const double deltx, const double delty, const int nx, const int ny, 
						   const double alpha, const double tilt); //!< Create two-dimensional layer of molecules
	
	bool AddElectrSurf(MolSet* pmset, int add_surf_below_flag, int add_surf_top_flag, int add_atom_top_flag,
		               int add_atom_below_flag); //!< Add molecular surface or a single donor/acceptor atom  to the structure

	std::string solv_name;       //!< Solvent to solvate molecular set
	double solv_buffer_dist;     //!< thickness of the solvent buffer around the solute 
	double min_solute_solv_dist; //!< minimal distance of added solvent molecules to the solute
	int num_na_add;  //!< Number of sodium ions to add
	int num_cl_add;  //!< Number of chloride ions to add
};

//! constants to specify different types of atom parameters to set  
enum AtomParams {
   BACKBONE_CHRG = 0x0001,             //!< Set backbone charges
   PROT_CHARGED_GROUPS_CHRG = 0x0002,  //!< Set charges only on the charged residues
   ZERO_CHRG =     0x0004,            //!< Set all charges zero
   AMBER_ALL_ATOM_CHRGS = 0x0008,     //!< Set All Atom Amber Force field charges
   AMBER_ALL_ATOM_FF_SYMBOLS = 0x0010, //!< Set All Atom Amber Force field atom symbols
   AMBER_ALL_ATOM_MASSES = 0x0020,     //!< Set All Atom Amber Force fieldAtomic Masses 
   ATOM_ELEM_FROM_TEMPL = 0x0040,      //!< Set Atom Elements From Residue templates
   ATOM_MASSES_ELEMENT = 0x0080,       //!< Set Atomic Masses corresponding to elements
   ATOM_HBOND_DA_STATUS = 0x0100        //!< Set H-Bond Donor/Acceptor status for atoms
}; 


#endif  /* !DEFINED MOLEDITOR_H_HARLEM */

