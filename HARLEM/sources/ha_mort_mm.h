/*! \file ha_mort_mm.h

   Classes and Function to define interface bewtween Molecular Mechanics parts MORT library and HARLEM.

  \author Igor Kurnikov 
  \date 2010-

*/
#if !defined(HA_MORT_MM_H)
#define HA_MORT_MM_H

namespace mort
{
	class database_t;
	class molecule_t;
}

class MortForceField
{
public:
	MortForceField();
	virtual ~MortForceField(); 

	mort::database_t* p_mdb;    //!< Database of residues in MORT format
	mort::molecule_t* p_atomff; //!< Parameters of force field
	mort::molecule_t* p_poleff; //!< Polarization related parameters of the force field

};

#endif //!defined(HA_MORT_MM_H)