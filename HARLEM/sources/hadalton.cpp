//  hadalton.cpp
//
//  Classes 
//  to provide interface to DALTON from HARLEM.
//
//  Implementation
//
//  Igor Kurnikov , University of Pittsburgh 
//
//  Revisions: May 28 1998
//

#define HADALTON_CPP

#include <mpi.h>

#include <assert.h>

#if defined(_MSC_VER)
#include <process.h>
#include <errno.h>
#endif

#include "f2c.h"
#include "haio.h"
#include "hamolecule.h"
#include "hamolset.h"
#include "haqchem.h"
#include "hadalton.h"
#include "tokens.h"

extern "C" {

extern void opendx_(int *ludx, char *name__, int *nelem, 
	int *nrec, char *status, int *lrdx, int *nbdx, 
	logical *olddx, ftnlen name_len, ftnlen status_len);

extern logical finddx_(int *lu, int *lrdx, int *i__, int *len, 
	                   int *ivec);

extern void readdx_(int *lu, int *lrdx, int *i__, 
	int *len, int *ivec);

extern void writdx_(int *lu, int *lrdx, int *i__, 
	int *len, int *ivec);

}


HaDaltonMod::HaDaltonMod(MolSet* new_phost_mset): 
HaCompMod(COMP_MOD_DALTON, new_phost_mset)
{
	if(new_phost_mset != NULL) 
		p_qc_mod = new_phost_mset->GetQCMod(true);
	SetStdFileNames();
	SetStdJobFlags();
}

HaDaltonMod::~HaDaltonMod()
{

}


void
HaDaltonMod::SetStdFileNames()
{
	Set_mol_file_name("MOLECULE.INP");
	Set_param_file_name("DALTON.INP");
}

void
HaDaltonMod::SetStdJobFlags()
{
	ABALNR_skip=false;
}


bool
HaDaltonMod::run(const RunMode rmode)
{
	int result,mode;
	Save_mol_file();
	Save_param_file();
#if defined(_MSC_VER) 
	if(rmode == RUN_BACKGROUND)
		mode= _P_NOWAIT;
	else if(rmode == RUN_FOREGROUND )
		mode= _P_WAIT;

	result=_spawnlp(mode,"daltonx",NULL);
	if(result == -1)
	{
		std::cout << " Error in submitting DALTON job " << std::endl;
		if(errno == E2BIG)
		{
			std::cout << " Argument list exceeds 1024 bytes " << std::endl;
		}
		else if( errno == EINVAL)
		{
			std::cout << " mode argument is invalid " << std::endl;
		}
		else if( errno == ENOENT)
		{
			std::cout << " File or path is not found " << std::endl;
		}
		else if( errno == ENOEXEC)
		{
			std::cout << " Specified file is not executable or has invalid executable-file format " << std::endl;
		}
		else if( errno == ENOMEM )
		{
			std::cout << " Not enough memory is available to execute new process " << std::endl;
		}
	}		
#else
	system("dalton.x");
#endif

	return true;
}


bool 
HaDaltonMod::Set_mol_file_name(const char* mf_name)
{
	mol_file_name = mf_name;
	return true;
}

bool HaDaltonMod::Set_param_file_name(const char* pf_name)
{
	param_file_name = pf_name;
	return true;
}


bool HaDaltonMod::Save_mol_file()
{
	char buf[120];
	std::ofstream mfile(mol_file_name.c_str());
	mfile << "ATOMBASIS" << std::endl;
	mfile << " Harlem generated Molecule file " << std::endl;
	mfile << " Test Molecule " << std::endl;

	char cbas_smb=' ';
	int na=      p_qc_mod->GetNumCnt();
	int charge=  p_qc_mod->GetCharge();

	MolSet* phmol_set= p_qc_mod->GetMolSet();
	assert(phmol_set != NULL);
	int nsym_gen= 0;
	char kasym[3][3];
	int i,j;
	for(i=0; i <= 3; i++)
	{
		for(j =0; j <= 3; j++)
			kasym[i][j]=' ';
	}
	int iang=1;
	double thr=0.0;
	
	sprintf(buf," %4d%3d%2d         %1d%10.2f",na,charge,nsym_gen,iang,thr);

	buf[0]=cbas_smb;

	HaAtom* aptr;

	mfile << buf << std::endl;
	AtomIteratorMolSet aitr(phmol_set);
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		double at_chg=(double)aptr->GetElemNo();
		int at_num=1;
		std::string bname="3-21G";
		sprintf(buf,"      %4.0f%5d      %s",at_chg,at_num,bname.c_str()); 
		mfile << buf << std::endl;
		std::string at_lbl=aptr->GetStdSymbol();
		double x = aptr->GetX_Ang();
		double y = aptr->GetY_Ang();
		double z = aptr->GetZ_Ang();
		sprintf(buf,"%s    %20.15f %20.15f %20.15f ",at_lbl.c_str(), x,y,z);
		mfile << buf << std::endl;
	}

	mfile << std::endl << std::endl;

	mfile.close();
	return true;
}



bool HaDaltonMod::Save_param_file()
{
	std::ofstream pfile(param_file_name.c_str());
	Fill_section_general(pfile);
	Fill_section_wave_fun(pfile);
	Fill_section_properties(pfile);
	pfile << "*END OF INPUT" << std::endl;
	pfile.close();
	return true;
}

void HaDaltonMod::Set_rotangle_calc()
{
	
}

void HaDaltonMod::Set_london_magdip_calc()
{
	Set_rotangle_calc();
	ABALNR_skip=true;
}

void HaDaltonMod::Fill_section_general(std::ofstream& pfile)
{
	pfile << "**GENERAL" << std::endl;
	pfile << ".RUN PROPERTIES" << std::endl;
	pfile << ".DIRECT" << std::endl;
}

void
HaDaltonMod::Fill_section_wave_fun(std::ofstream& pfile)
{
	pfile << "**WAVE FUNCTIONS" << std::endl;
	pfile << ".HF" << std::endl;
	Fill_subsection_HF_param(pfile);
}

void
HaDaltonMod::Fill_subsection_HF_param(std::ofstream& pfile)
{
	pfile << "*HF INPUT" << std::endl;
	pfile << ".THRESH" << std::endl;
	pfile << "1.0D-8" << std::endl;
}


void
HaDaltonMod::Fill_section_properties(std::ofstream& pfile)
{
	pfile << "**PROPERTIES" << std::endl;
	pfile << ".NOCMC" << std::endl;
	pfile << ".VROA"  << std::endl;
	pfile << ".ABALNR" << std::endl;
	Fill_subsection_ABALNR(pfile);
}

void
HaDaltonMod::Fill_subsection_ABALNR(std::ofstream& pfile)
{
	pfile << "*ABALNR" << std::endl;
	if(ABALNR_skip) pfile << ".SKIP"   << std::endl;
	pfile << ".ROTANG" << std::endl;
	pfile << ".THRLNR" << std::endl;
	pfile << "5.0D-8"  << std::endl;
    pfile << ".FREQUE" << std::endl;
    pfile << " 1 "     << std::endl;
    pfile << " 0.0773179772 " << std::endl;
}

bool
HaDaltonMod::test_1()
{
	return true;
}

#if 0

int 
DaltonDXfile::open(const char* fname, const char* mode)
{
    opendx_(const char* fname, const char* mode);
}

int 
DaltonDXfile::read()
{
	readdx_();
}

int 
DaltonDXFile::write()
{
	writedx_();
}

#endif



