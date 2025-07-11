/*!  \file mm_driver_arbalest.h

    Molecular Mechanics simulations using ARBALEST package  

    \author Igor Kurnikov 
    \date 2025-

*/ 

#if !defined MM_DRIVER_ARBALEST_H
#define MM_DRIVER_ARBALEST_H

#include "hamolmech.h"

class MMDriverArbalest : public MMDriver
{
public:
	MMDriverArbalest(MolSet* pmset);
	MMDriverArbalest(HaMolMechMod* p_mm_mod_new);
	virtual ~MMDriverArbalest();

	virtual std::string GetClassName() { return "MMDriverArbalest"; }

	virtual int CalcEnergy() { return FALSE;} //!< Calculate energy of the system and save results to p_mm_info member of p_mm_mod
	virtual int SaveAllInpFiles(); //!< Save input files for Arbalest

	bool InitForceField(std::string ff_name); //!< Initialize Force Field on the system with a name string

	void SetFileNamesWithPrefix(std::string prefix); //!< Set ARBALEST input and output file names with prefix
	
	bool SaveConfigFile();      //!< Save MM run parameters in ARBALEST CONFIG format to File
	bool SaveRunFiles();        //!< Save ARBALEST Run script
	
	int SaveConfigToStream( std::ostream& os );    //!< Save MM run parameters in ARBALEST format to std::stream

	std::vector<std::string> GetPosRestraintsDescAndList(std::set<HaAtom*>& atoms_saved); //!< Get Positional Restraints Description and List in Arbalest Format
	bool SavePosRestraintsStream(std::ostream& os_desc, std::ostream& os_list, std::set<HaAtom*>& saved_atoms); //!< Write Positional Restraints Description and List in Arbalest Format to Streams

	bool SaveMolDefToStream(std::ostream& os , HaMolecule* pmol, std::set<std::string>& mol_defined, std::string pos_restr_desc); //!< Save Description of a molecule ( specified by 0-based index) to a stream in Arbalest Config format
	bool SaveStdMolDefToStream(std::ostream& os, std::string mol_name); //!< Save a description of a standard molecule ( with HIN files saved in Input/HIN )

	std::string arbalest_exe; // !< ARBALEST executable path
	std::string config_fname;
	std::string output_dir_name;
	std::string trj_fname;
	std::string ene_fname;
	std::string run_fname;

	const std::set<std::string> std_mol_names = { "HOH","NA","CL" };

protected:


};


#endif // !defined MM_DRIVER_GROMACS_H
