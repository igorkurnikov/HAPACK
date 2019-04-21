/*!  hacompmod.cpp
    
    Base class of the computational modules in HARLEM 

    \author Igor Kurnikov
	\date 1999-2002
*/


#include "stdlib.h"
#if defined(_WIN32)
#include <process.h>
#endif

#include "hampi.h"
#include "hacompmod.h"
#include "electrostmod.h"
#include "etcoupl.h"
#include "haqchem.h"
#include "hagaussian.h"
#include "hadalton.h"
#include "haintermol.h"
#include "hamolmech.h"
#include "hascattermod.h"
#include "stmmod.h"
#include "nuclacidmod.h"
#include "hazindo.h"
#include "protonredox.h"
#include "haempirical.h"
#include "apbsmod.h"
#include "elmod.h"
#include "haflexmod.h"
#include "haproteined.h"

HaCompMod::HaCompMod( const int new_mtype, HaMolSet* new_phost_mset):
mtype(new_mtype)
{
	phost_mset = new_phost_mset;
	debug_level = 1;
}



HaCompMod::~HaCompMod()
{
	
}


HaCompMod* HaCompMod::CreateCompMod( const int mtype, HaMolSet* new_phost_mset )
{
	HaCompMod* pmod=NULL;

	if( mtype == COMP_MOD_ELECTROST )
	{
		pmod= new ElectrostMod(new_phost_mset);
	}
	else if ( mtype == COMP_MOD_EL )
	{
#ifdef ELMOD_COMPILE
		pmod = new ElMod(new_phost_mset);
#endif
	}
	else if ( mtype == COMP_MOD_PKA_CALC)
	{
		pmod = new pKaCalcMod(new_phost_mset);
	}
	else if ( mtype == COMP_MOD_PNP )
	{
		pmod = new PNPMod(new_phost_mset);
  }
	else if ( mtype == COMP_MOD_ET_COUPL )
	{
		pmod = new ETCouplMod(new_phost_mset);
	}
	else if ( mtype == COMP_MOD_QCHEM )
	{
		pmod = new HaQCMod(new_phost_mset);
	}
	else if ( mtype == COMP_MOD_GAUSSIAN )
	{
		pmod = new HaGaussMod(new_phost_mset);
	}
	else if ( mtype == COMP_MOD_DALTON )
	{
		pmod = new HaDaltonMod(new_phost_mset);
	}
	else if ( mtype == COMP_MOD_INTERMOL )
	{
		pmod = new HaInterMolMod(new_phost_mset);
	}
	else if ( mtype == COMP_MOD_MOLMECH )
	{
		pmod = new HaMolMechMod(new_phost_mset);
	}
	else if ( mtype == COMP_MOD_SCATTER )
	{
		pmod = new HaScatterMod(new_phost_mset);
	}
	else if ( mtype == COMP_MOD_STM )
	{
		pmod = new StmMod(new_phost_mset);
	}
	else if ( mtype == COMP_MOD_NUCL_ACID )
	{
		pmod = new NuclAcidMod(new_phost_mset);
	}
	else if ( mtype == COMP_MOD_ZINDO )
	{
		pmod = new HaZindoMod(new_phost_mset);
	}
	else if ( mtype == COMP_MOD_PROTON_REDOX )
	{
		pmod = new ProtonRedoxMod(new_phost_mset);
	}
	else if ( mtype == COMP_MOD_EMPIRICAL)
	{
		pmod = new HaEmpiricalMod(new_phost_mset);
	}
	else if (mtype == COMP_MOD_MEMBRANE) // jose
	{
		pmod = new HaMolMembraneMod(new_phost_mset);
	}
	else if ( mtype == COMP_MOD_APBS)
	{
		pmod = new APBSMod(new_phost_mset);
	}
	else if( mtype == COMP_MOD_FLEX )
	{
		pmod = new HaFlexMod(new_phost_mset);
	}
	else if( mtype == COMP_MOD_CLUSTER_ANAL )
	{
		pmod = new CollectCrdAnalMod(new_phost_mset);
	}
	else
	{
		pmod= NULL;
	}
	
	return pmod;
}

void HaCompMod::SetDebugLevel(int new_debug_level)
{
	debug_level = new_debug_level;
}

int HaCompMod::SaveXMLToStream(std::ostream& os, const harlem::SaveOptions* popt ) const
{
	return FALSE;
}

int HaCompMod::OnDelAtoms(AtomContainer& del_atoms)
{
	return TRUE;
}
