/*! \file hacompmod.h
 
     base class for computational modules in HARLEM  

    \author Igor Kurnikov
    \date 1999-2007

*/

#ifndef HACOMPMOD_H
#define HACOMPMOD_H

#include "haconst.h"
#include "hastl.h"
#include "hastring.h"
#include "hatypes.h"

class HaMolSet;
class AtomContainer;

const int COMP_MOD_ELECTROST   = 0x0001; //!< mtype for the continuum electrostatics module
const int COMP_MOD_ET_COUPL   = 0x0002; //!< mtype for the electron transfer module
const int COMP_MOD_QCHEM      = 0x0003; //!< mtype for the quantum chemical module
const int COMP_MOD_GAUSSIAN   = 0x0004; //!< mtype for the Gaussian interaction module
const int COMP_MOD_DALTON     = 0x0005; //!< mtype for the Dalton interaction module
const int COMP_MOD_INTERMOL   = 0x0006; //!< mtype for the intermolecular interaction module
const int COMP_MOD_MOLMECH    = 0x0007; //!< mtype for the molecular mechanics module
const int COMP_MOD_SCATTER    = 0x0008; //!< mtype for the scattering computations module
const int COMP_MOD_STM        = 0x0009; //!< mtype for the STM simulations module
const int COMP_MOD_NUCL_ACID  = 0x0010; //!< mtype for the DNA modeling module
const int COMP_MOD_ZINDO      = 0x0011; //!< mtype for ZINDO interface module
const int COMP_MOD_PROTON_REDOX = 0x0012; //!< mtype for Chemical Transformations module
const int COMP_MOD_EMPIRICAL  = 0x0013;  //!< mtype for Empirical calculation module
const int COMP_MOD_PNP  = 0x0014;        //!< mtype for the continuum electrostatics module(PNP)
const int COMP_MOD_APBS  = 0x0015;        //!< mtype for APBS
const int COMP_MOD_EL  = 0x0016;          //!< mtype for new continuum electrostatics module
const int COMP_MOD_PKA_CALC  = 0x0017;    //!< mtype for pKa Calculation module
const int COMP_MOD_MEMBRANE  = 0x0018;    //!< mtype for Membrane model energy calculations jose
const int COMP_MOD_FLEX      = 0x0019;    //!< mtype for Molecular Flexibility Module
const int COMP_MOD_CLUSTER_ANAL = 0x0020; //!< mtype for Clustering Analysis Module

namespace harlem 
{ 	
	class SaveOptions;

	class RunOptions : public harlem::HashMap
	{
	public:
		RunOptions() { SetStdOptions(); }
		RunOptions( const RunOptions& ref ) { SetStdOptions(); Copy(ref); }
		virtual ~RunOptions() {}

		void SetStdOptions() { SetRunSync(true); SetSaveInpFile(true); SetLoadOutput(true); }

		virtual void Copy( const harlem::HashMap& ref ) 
		{  
			harlem::HashMap::Copy(ref);
			const RunOptions* pref = dynamic_cast<const RunOptions*>(&ref);
			if( pref != NULL )
			{
				SetRunSync( pref->ToRunSync() ); 
				SetSaveInpFile( pref->ToSaveInpFile() ); 
				SetLoadOutput( pref->ToLoadOutput() );
			}
		}

		harlem::HashMap* clone() const { return new RunOptions(*this); }  

		void SetRunSync( bool set_flag = true ) { run_sync_flag = set_flag; } 
		bool ToRunSync() const { return run_sync_flag; } 

		void SetSaveInpFile( bool set_flag = true ) { save_inp_file_flag = set_flag; }
		bool ToSaveInpFile() const { return save_inp_file_flag; }

		void SetLoadOutput( bool set_flag = true ) { load_output_flag = set_flag; }
		bool ToLoadOutput() const { return load_output_flag; }

	protected:
		bool run_sync_flag; //!< Run Program or function synchroneously 
		bool save_inp_file_flag; //!< Save Input Files for the program before launching the program
		bool load_output_flag;   //!< Load result information from output files
	};
}

 
class HaCompMod
//! Parent class for computational modules in HARLEM
{
public:
	HaCompMod(const int new_mtype, HaMolSet* new_phost_mset = NULL); //!< Constructor - supply type in new_mtype
	virtual ~HaCompMod();

	static HaCompMod* CreateCompMod( const int mtype, HaMolSet* new_phost_mset = NULL); //!< Create computational module of the given type

	bool SetMolHost(HaMolSet* new_phost_mset) { phost_mset = new_phost_mset; return true; }
	HaMolSet* GetMolSet() { return phost_mset; }             //!< Get Molecular Set associated with the module
	const HaMolSet* GetMolSet() const { return phost_mset; } //!< Get Molecular Set associated with the module

	int GetType() { return mtype; } //!< return type of the computational module

	virtual void SetDebugLevel(int new_debug_level);
	virtual int SaveXMLToStream(std::ostream& os, const harlem::SaveOptions* popt = NULL ) const; 

	virtual int OnDelAtoms(AtomContainer& del_atoms); //!< Modify module content to react to deleted atoms 
	
	int debug_level;
	
protected:
	HaMolSet* phost_mset; //!< pointer to the Molecular Set associated with the module
	const int mtype;      //!< the type of the computational module 
	
};

#endif // end if !defined(HACOMPMOD_H)
