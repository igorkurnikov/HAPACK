/*!  \file mm_driver_tinker.h

    Molecular Mechanics simulations using TINKER package  

    \author Igor Kurnikov 
    \date 2010-

*/ 

#if !defined MM_DRIVER_TINKER_H
#define MM_DRIVER_TINKER_H

#include "hamolmech.h"

class MMDriverTinker : public MMDriver
{
public:
	MMDriverTinker(HaMolMechMod* p_mm_mod_new);
	virtual ~MMDriverTinker();

	virtual std::string GetClassName() { return "MMDriverTinker"; }

	virtual int CalcEnergy() { return FALSE;} //!< Calculate energy of the system and save results to p_mm_info member of p_mm_mod
	virtual int SaveAllInpFiles(); //!< Save input files for Tinker 

	int InitFFNumMap(); //!< Initiate Tinker FF maps
	int LoadFFNumMapFromFile( const char* fname); //!< Load map of atom names to Tinker ff numbers from XML file 

	StrIntMap atn_ff_num_map; //!< Map of full atom names ( as AMBER94:ALA#NT.CA ) to TINKER ff number
	int GetAtFFNum(const char* force_field_name, const char* res_templ_name, const char* at_name); //!< Get TINKER atom force field number 

protected:

	int to_init_ff_map; //!< Flag to indicate that atom name -> Tinker atom ff numbers should be initiated 

};


#endif // !defined MM_DRIVER_TINKER_H