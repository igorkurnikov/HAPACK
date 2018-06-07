//  hadalton.h
//
//  Classes 
//  to provide interface to DALTON from HARLEM.
//
//  definitions
//
//  Igor Kurnikov , University of Pittsburgh 
//
//  Created:    July 8 1998
//
#ifndef HADALTON_H
#define HADALTON_H


class HaQCMod;

#include "hastring.h"
#include "hacompmod.h"
#include "hagaussian.h" // for enum RunMode

#if 0 // comment DaltonDXFile for now, don't know how to calculate 
      // record length necessary for using opendx,readx,writedx functions

class DaltonDXFile
//  Class to define Hydrogen Bond object 
{
public:

	int  open(const char* fname, const char* mode);         
	int  read();          
	int  write();
	int  close();

protected:

	int lugdi;
	int lrgdi;
	int nbgdi;

};

#endif

class HaDaltonMod :  public HaCompMod
//!  Class to control Quantum Chemical calculations using Dalton program   
{
public:

	HaDaltonMod(HaMolSet* new_phost_mset);
	virtual ~HaDaltonMod();

	void SetStdFileNames();
	void SetStdJobFlags();

	bool SetInpFilePrefix(const char* new_prefix) { return false; }
	bool Set_mol_file_name(const char* mf_name);
	bool Save_mol_file();

	bool Set_param_file_name(const char* pf_name);
	bool Save_param_file();

	bool run(const RunMode rmode);
	
	void Set_rotangle_calc();
	void Set_london_magdip_calc();

	bool test_1();

protected:
    
	std::string mol_file_name;
	std::string param_file_name;
	HaQCMod* p_qc_mod;

// Flags for a job types:
	bool ABALNR_skip;
	
// Fill different sections of input files:

	void Fill_section_general(ofstream& pfile);
	void Fill_section_wave_fun(ofstream& pfile);
	void Fill_subsection_HF_param(ofstream& pfile);
	void Fill_section_properties(ofstream& pfile);
	void Fill_subsection_ABALNR(ofstream& pfile);

};



#endif /* !HADALTON_H */
