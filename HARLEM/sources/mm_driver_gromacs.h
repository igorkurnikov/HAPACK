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
	
	int SaveTopFile();                      //!< Save Molecular Sysytem Topology in GROMACS format to File
	int SaveTopToStream(std::ostream& os);  //!< Save Molecular Sysytem Topology in GROMACS format to Stream

	int SaveAtomTypesToStream(std::ostream& os); //!< Save [ atomtypes ] section of GROMACS topology
	int SaveAtomsToStream(std::ostream& os);     //!< Save [ atoms ] section of GROMACS topology 
	int SaveBondsToStream(std::ostream& os);     //!< Save [ bonds ] section of GROMACS topology
	int SavePairsToStream(std::ostream& os);     //!< Save [ pairs ] section of GROMACS topology  ( pairs of atoms with 1-4 intramolecular interactions )
	int SaveAnglesToStream(std::ostream& os);    //!< Save [ angles ] section of GROMACS topology
	int SaveDihedralsToStream(std::ostream& os); //!< Save [ dihedrals ] (proper & improper) section of GROMACS topology

protected:

	MolSet*     pmset;
	MolMechModel* p_mm_model;      //!< MolMechModel corresponding to the class
	HaMolMechMod* p_mm_mod;        //!< HaMolMech Module associated with the model
};


#endif // !defined MM_DRIVER_GROMACS_H