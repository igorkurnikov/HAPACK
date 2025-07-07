/*! \file hazindo.h

    Interface to ZINDO from HARLEM.

    \author Igor Kurnikov  
    \date 1998-2004

*/
#ifndef HAZINDO_H
#define HAZINDO_H

#include "hacompmod.h" 

class HaZindoMod: public HaCompMod
//!  Class to control semiempirical quantum chemical calculations with ZINDO program   
{
public:
	HaZindoMod(HaMolSet* new_phost_mset);
	virtual ~HaZindoMod();

	int SaveParamFile(const char* fname); //!< Save Zindo Input Parameter File

protected:

	HaQCMod* p_qc_mod;     //!< The pointer to Quantum Chemical module associated with the Gaussian module 

};


#endif // end of ifndef HAZINDO_H

