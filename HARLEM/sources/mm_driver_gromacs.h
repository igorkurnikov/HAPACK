/*!  \file mm_driver_gromacs.h

    Molecular Mechanics simulations using GROMACS package  

    \author Igor Kurnikov 
    \date 2018-

*/ 

#if !defined MM_DRIVER_GROMACS_H
#define MM_DRIVER_GROMACS_H

#include "hamolmech.h"

class MMDriverGromacs : public MMDriver
{
public:
	MMDriverGromacs(HaMolMechMod* p_mm_mod_new);
	virtual ~MMDriverGromacs();

	virtual std::string GetClassName() { return "MMDriverGromacs"; }

	virtual int CalcEnergy() { return FALSE;} //!< Calculate energy of the system and save results to p_mm_info member of p_mm_mod
	virtual int SaveAllInpFiles(); //!< Save input files for Gromacs

	void PartitionAtomsToMolecules(); //!< Partition Atoms to GROMACS molecules 

	void SetFileNamesWithPrefix(std::string prefix); //!< Set GROMACS input and output file names with prefix
	
	bool SaveMdpFile();        //!< Save MM run parameters in GROMACS MDP format to File
	int  SaveGromacsTopFile();  //!< Save Molecular System Topology in GROMACS format to File
	bool SaveInitCrdFiles();    //!< Save Initial Coordinates Files for the Run - possibly from previous MD or MIN runs 
	bool SaveRunFile();        //!< Save GROMACS Run script
	
	int SaveMdpToStream( std::ostream& os );    //!< Save MM run parameters in GROMACS MDP format to std::stream
	int SaveGromacsTopToStream( std::ostream& os );    //!< Save Molecular Sysytem Topology in GROMACS format to std::stream

	int  SaveAtomTypesToStream(std::ostream& os);  //!< Save [ atomtypes ] section of GROMACS topology
	bool SaveAtomsToStream    (std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map); //!< Save [ atoms ] section of GROMACS topology 
	bool SaveBondsToStream    (std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map); //!< Save [ bonds ] section of GROMACS topology
	bool Save14PairsToStream    (std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map); //!< Save [ pairs ] section of GROMACS topology  ( pairs of atoms with 1-4 intramolecular interactions )
	bool SaveAnglesToStream   (std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map); //!< Save [ angles ] section of GROMACS topology
	bool SaveDihedralsToStream(std::ostream& os, AtomGroup& group, AtomIntMap& at_idx_map); //!< Save [ dihedrals ] (proper & improper) section of GROMACS topology

	std::string inp_fname;
	std::string top_fname;
	std::string init_crd_fname;
	std::string restr_crd_fname;
	std::string trj_fname;
	std::string ene_fname;
	std::string run_fname;
	std::string tpr_fname;

	static std::set<std::string> std_gmx_mols;

protected:

	std::vector<AtomGroup> gmx_mol_partition;
	MolSet*     pmset;
	MolMechModel* p_mm_model;      //!< MolMechModel corresponding to the class
	HaMolMechMod* p_mm_mod;        //!< HaMolMech Module associated with the model
};


#endif // !defined MM_DRIVER_GROMACS_H