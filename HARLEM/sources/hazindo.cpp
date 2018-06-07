/*!
    Classes to provide interface to ZINDO from HARLEM.

    Implementation

    \author Igor Kurnikov 
    \date 2004-
*/

#define HAZINDO_CPP

#include <mpi.h>

#ifdef _MSC_VER
#include <process.h>
#endif

#include "hamolset.h"
#include "hacompmod.h"
#include "haqchem.h"
#include "hazindo.h"

HaZindoMod::HaZindoMod(HaMolSet* new_phost_mset):
HaCompMod(COMP_MOD_ZINDO, new_phost_mset)
{
	if(new_phost_mset != NULL) p_qc_mod = new_phost_mset->GetQCMod(true); 
}

HaZindoMod::~HaZindoMod()
{

}

int 
HaZindoMod::SaveParamFile(const char* fname)
{
   FILE* fpar = fopen(fname,"w");
   if(fpar == NULL) return FALSE;  
   
   fclose(fpar);
   return TRUE;
}