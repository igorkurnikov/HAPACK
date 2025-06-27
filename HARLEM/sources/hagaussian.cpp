/*!
    Classes to provide interface to Gaussian from HARLEM.

    Implementation

    \author Igor Kurnikov 
    \date 1998-
*/

#define HAGAUSSIAN_CPP

#include <mpi.h>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include <wx/filename.h>
#include <wx/process.h>

#include <assert.h>
#include "f2c.h"
#include "haio.h"
#include "halinalg.h"
#include "haatgroup.h"
#include "hamolecule.h"
#include "hamolset.h"
#include "hacoord.h"
#include "haintcrd.h"
#include "haqchem.h"
#include "hagaussian.h"
#include "hapseudopot.h"
#include "tokens.h"
#include "hasurface.h"
#include "gaufile.h"
#include "math.h"


harlem::RunOptions HaGaussMod::run_opt_default;

HaGaussMod::HaGaussMod(MolSet* pmset_new): HaCompMod(COMP_MOD_GAUSSIAN, pmset_new )
{
	if( pmset_new != NULL) p_qc_mod = pmset_new->GetQCMod(true); 
	SetStdFileNames();
	SetStdJobFlags();
}

HaGaussMod::~HaGaussMod()
{

}

void HaGaussMod::SetStdFileNames()
{
	gaussian_version = "09";
	gaussian_exe     = (std::string)"g" + gaussian_version;
}

void HaGaussMod::SetStdJobFlags()
{
	pseudo_pot_flag = false;

	SetNumSharedMemCores(2);         
	SetNumProc(1);
	SetMaxMem(1000);

	SetLoadNonOptGeom( true );
	SetLoadGeomZMatOrient( true );
	SetLoadGeomStdOrient ( false );
	SetReadInitGeomChkFile( false );
	SetReadHFGuessChkFile( false );
	SetNoStdOrient( false );
	SetSaveBasisSetGen( false );

}


int HaGaussMod::SetFilePrefix(const char* new_prefix)
{
	inp_file_prefix = new_prefix;
	return TRUE;
}

std::string HaGaussMod::GetFilePrefix() const 
{
	std::string prefix_loc = inp_file_prefix;
	if( prefix_loc.empty() ) prefix_loc = p_qc_mod->GetPrefix();
	return prefix_loc;
}

std::string HaGaussMod::GetInpFileName() const
{
	return (GetFilePrefix() + ".gjf");
}


std::string HaGaussMod::GetCHKFileName() const
{
	return (GetFilePrefix() + ".chk");
}

std::string HaGaussMod::GetFCHKFileName() const
{
	return (GetFilePrefix() + ".fchk");
}

std::string HaGaussMod::GetRWFFileName() const
{
	return (GetFilePrefix() + ".rwf");
}

std::string HaGaussMod::GetOutFileName() const
{
	return (GetFilePrefix() + ".out");
}

bool HaGaussMod::SaveInpFile()
{
    MolSet* pmset = this->GetMolSet();
	std::string str;

	std::ofstream inp_file_stream( GetInpFileName().c_str() );
	if(inp_file_stream.fail())
	{
		PrintLog("Error in HaGaussMod::SaveInpFile() \n");
		PrintLog("Fail to to open file %s \n",GetInpFileName().c_str());
		return false;
	}
	try
	{
		FillSectionProcCommands(inp_file_stream);
		FillSectionJob(inp_file_stream);
		FillSectionCoord(inp_file_stream);
		FillSectionExtCharges(inp_file_stream);
		FillSectionBasis(inp_file_stream);
		FillSectionInitMO(inp_file_stream);
	}
	catch( std::exception& ex)
	{
		PrintLog("Error in HaGaussMod::SaveInpFile() \n");
		PrintLog("%s \n",ex.what());
		return false;
	}

	return true;
}

class GaussianProcess : public wxProcess
{
public:
	GaussianProcess() {  p_qc_mod = NULL; }

	HaQCMod* p_qc_mod;
	
	virtual void OnTerminate(int pid, int status)
	{
		if(p_qc_mod) p_qc_mod->StopCalc();
		PrintLog("GAUSSIAN Process Has Stopped \n");
	}
};

int HaGaussMod::Run( const harlem::RunOptions* popt_par )
{
	const harlem::RunOptions* popt = popt_par;
	if( popt == NULL ) popt = &run_opt_default; 

	int result;
	if( popt->ToSaveInpFile() ) SaveInpFile();

	wxString cmd_line = gaussian_exe.c_str();
	cmd_line += " ";
	cmd_line += GetInpFileName();

	PrintLog(" HaGaussMod::Run()  cmd_line: ");
	PrintLog(" %s \n", cmd_line.ToStdString().c_str() );

	GaussianProcess* p_gauss_proc = new GaussianProcess();
	p_gauss_proc->p_qc_mod = this->p_qc_mod;

	int res;
	if( popt->ToRunSync() )
	{
		res = wxExecute(cmd_line,wxEXEC_SYNC,p_gauss_proc);
	}
	else
	{
		res = wxExecute(cmd_line,wxEXEC_ASYNC,p_gauss_proc);
	}

	if( popt->ToLoadOutput() ) LoadOutput();

	return TRUE;
}

int HaGaussMod::LoadOutput()
{
	LoadOutFile( GetOutFileName() );
	return TRUE;	
}

int HaGaussMod::LoadOutFile(const std::string& fname)
{
	std::ifstream is(fname.c_str());
	if(is.fail()) 
	{
		PrintLog(" Error to open Gaussian out file %s \n", fname.c_str() );
		return FALSE;
	}
	return LoadOutFromStream( is );
}

static int skip_lines( std::istream& is, std::string& line, int n )
{
	int i;
	for( i = 0; i < n; i++)
	{
		std::getline(is,line);
		if( is.eof() ) return FALSE;
	}
	return TRUE;
}
	
int HaGaussMod::LoadOutFromStream( std::istream& is )
{
	std::string line;
	try
	{	
		enum { SEARCH_MODE = 0, READ_OPT_GEOM = 1, READ_SUMMARY = 2 } read_mode0;
		read_mode0 = SEARCH_MODE;

		MolSet* pmset = this->GetMolSet();
		int na = pmset->GetNAtoms();

		if( p_qc_mod->IsEneMinCalc() ) p_qc_mod->ene_history.clear(); 

		HaVec_double crd_tmp;
		std::vector<std::string> str_vec;
		int ires;
		int i;

		bool opt_geom = true;

		for(;;)
		{	
			std::getline(is,line);
			if( is.eof() ) break;

			if( read_mode0 == SEARCH_MODE )
			{
				if( line.find("Optimization completed.") != std::string::npos ) 
				{
					std::getline(is,line);
					if( line.find( "Stationary point found") == std::string::npos )
					{
						PrintLog("Warning in HaGaussMod::LoadOutFromStream() \n");
						PrintLog("Optimization is complete but stationary point is not found \n");
						continue;
					}
					opt_geom = true;
					read_mode0 = READ_OPT_GEOM;
					continue;
				}
				if( load_non_opt_geom && line.find("Optimization stopped.") != std::string::npos ) 
				{
					read_mode0 = READ_OPT_GEOM;
					opt_geom = false;
					continue;
				}
				if( line.find("SCF Done:") != std::string::npos ) 
				{
					size_t pi = line.find("=");
					size_t pe = line.find("A.U.");
					std::string ene_str = line.substr(pi+1, pe - pi - 2);
					double ene = 0.0; 
					if( harlem::IsFloat(ene_str) ) ene = harlem::ToDouble(ene_str);

					PrintLog(" SCF ene = %12.6f\n",ene);
					if( p_qc_mod->IsEneMinCalc() ) p_qc_mod->ene_history.push_back(ene);

					if( p_qc_mod->IsHF() || p_qc_mod->IsMP2() )  p_qc_mod->SetHFEne(ene);
					if( p_qc_mod->IsDFT() ) p_qc_mod->SetDFTEne(ene);

					if( p_qc_mod->IsHF() || p_qc_mod->IsDFT() ) p_qc_mod->SetEne(ene);
				}
				if( line.find("l9999") != std::string::npos )
				{
					std::string summary;
					bool start_summary = false;
					for(;;)
					{
						std::getline(is,line);
						if( is.eof() ) break;
						boost::trim(line);
						bool find_vert = (line.find("|") != std::string::npos );
						if( !start_summary && !find_vert ) continue;
						if( line.empty() ) 
						{
							if( start_summary ) break;
							continue;
						}
						
						summary += line;
						if( line[line.size() - 1] == '@' ) break; 
					}
					LoadOutSummary(summary);
				}

				if( line.find("Dipole moment") != std::string::npos ) 
//				if( line.size() > 14 && line.substr(14) == " Dipole moment" ) 
				{
					std::getline(is,line);
					boost::trim(line);
					boost::split(str_vec,line,boost::is_any_of(" "),boost::token_compress_on);

					if( str_vec.size() < 6 || !harlem::IsFloat(str_vec[1]) ||
						!harlem::IsFloat(str_vec[3]) || !harlem::IsFloat(str_vec[5]) ) 
					{
						PrintLog(" Error reading Dipole moment from Gaussian Output file \n" );
						continue;
					}
					p_qc_mod->dipole[0] = atof( str_vec[1].c_str() )/EL_ANG_TO_DEBYE;
					p_qc_mod->dipole[1] = atof( str_vec[3].c_str() )/EL_ANG_TO_DEBYE;
					p_qc_mod->dipole[2] = atof( str_vec[5].c_str() )/EL_ANG_TO_DEBYE;
					continue;
				}
				if( line.find("Quadrupole moment") != std::string::npos ) 
//				else if( line.size() > 18 && line.substr(18) == " Quadrupole moment") 
				{
					std::getline(is,line);
					boost::trim(line);
					boost::split(str_vec,line,boost::is_any_of(" "),boost::token_compress_on);

					if( str_vec.size() < 6 || !harlem::IsFloat(str_vec[1]) ||
						!harlem::IsFloat(str_vec[3]) || !harlem::IsFloat(str_vec[5]) ) 
					{
						PrintLog(" Error reading Quadrupole moment from Gaussian Output file \n" );
						continue;
					}
					p_qc_mod->qpole[0] = atof( str_vec[1].c_str() )/EL_ANG_TO_DEBYE;
					p_qc_mod->qpole[1] = atof( str_vec[3].c_str() )/EL_ANG_TO_DEBYE;
					p_qc_mod->qpole[2] = atof( str_vec[5].c_str() )/EL_ANG_TO_DEBYE;
					
					std::getline(is,line);
					boost::trim(line);
					boost::split(str_vec,line,boost::is_any_of(" "),boost::token_compress_on);

					if( str_vec.size() < 6 || !harlem::IsFloat(str_vec[1]) ||
						!harlem::IsFloat(str_vec[3]) || !harlem::IsFloat(str_vec[5]) ) 
					{
						PrintLog(" Error reading Quadrupole moment from Gaussian Output file \n" );
						continue;
					}
					p_qc_mod->qpole[3] = atof( str_vec[1].c_str() )/EL_ANG_TO_DEBYE;
					p_qc_mod->qpole[4] = atof( str_vec[3].c_str() )/EL_ANG_TO_DEBYE;
					p_qc_mod->qpole[5] = atof( str_vec[5].c_str() )/EL_ANG_TO_DEBYE;

					continue;
				}
				if( p_qc_mod->IsMP2() )
				{
					if( line.find("EUMP2 =") != std::string::npos )
					{
						size_t pi = line.find("EUMP2 =");
						pi += 8;
						std::string ene_str = line.substr(pi);
						boost::trim(ene_str);
						size_t pos_d = ene_str.find('D');
						if( pos_d != std::string::npos ) ene_str[pos_d] = 'E'; 

						double ene_mp2 = harlem::ToDouble(ene_str);

						PrintLog( " mp2 ene = %12.6f \n",ene_mp2);

						p_qc_mod->SetEne(ene_mp2);
					}
				}
			} // end if(read_mode == SEARCH_MODE)
			else if( read_mode0 == READ_OPT_GEOM ) 
			{
				if( load_geom_zmat_orient && line.find("Z-Matrix orientation:") != std::string::npos )
				{
					ires = skip_lines(is, line, 4 );
					if( !ires ) throw std::runtime_error("Failed to read Optimized geometry in Z-matrix Orientation (skipping lines)");

					ZMatCrd* pzm = p_qc_mod->GetZMat();
					if( pzm == NULL || pzm->IsEmpty() ) throw std::runtime_error("Failed to read Optimized geometry in Z-matrix Orientation (No Zmatrix present)");

					int nz = pzm->GetNZ();
					for( i = 0; i < nz; i++ )
					{
						std::getline(is,line);
						if( is.eof() ) throw std::runtime_error("Failed to read optimized geometry in Z-matrix Orientation (eof) ");
						boost::trim(line);
						boost::split(str_vec,line,boost::is_any_of(" ,\r\n"),boost::token_compress_on);
						if( str_vec.size() != 6 ) throw std::runtime_error("Failed to read optimized geometry in Z-matrix Orientation ( #tokens != 6)");
			
						double x = atof( str_vec[3].c_str() ); 
						double y = atof( str_vec[4].c_str() ); 
						double z = atof( str_vec[5].c_str() ); 
						pzm->SetXYZCrd(i,x,y,z);
					}
					pzm->SetFromAtomCrd();

					std::string snap_name = "OPT_GEOM_QCHEM";
					if( !opt_geom )  snap_name = "NOPT_GEOM_QCHEM";
					CrdSnapshot* p_snap = pmset->GetCrdSnapshotByName( snap_name.c_str() ,true);
					if (p_snap == NULL ) throw std::runtime_error("Error to add optimized coordinate snapshot");
		
					ires = p_snap->SaveCurrentAtomCrd();
					if( !ires )  throw std::runtime_error("Error to save optimized coordinates to snapshot"); 
				
					if( !load_geom_std_orient ) 
					{
						read_mode0 = SEARCH_MODE;
					}
				}
				if( load_geom_std_orient && line.find("Standard orientation:") != std::string::npos )
				{
					ires = skip_lines(is, line, 4 );
					if( !ires ) throw std::runtime_error("Failed to read Optimized geometry (skipping lines)");

					crd_tmp.clear();
					crd_tmp.reserve(3*na);
					for( i = 0; i < na; i++ )
					{
						std::getline(is,line);
						if( is.eof() ) throw std::runtime_error("Failed to read optimized geometry in Standard Orientation (eof) ");
						boost::trim(line);
						boost::split(str_vec,line,boost::is_any_of(" ,\r\n"),boost::token_compress_on);
						if( str_vec.size() != 6 ) throw std::runtime_error("Failed to read optimized geometry in Standard Orientation ( #tokens != 6)");
						double fval;
						fval = atof( str_vec[3].c_str() ); crd_tmp.push_back( fval );
						fval = atof( str_vec[4].c_str() ); crd_tmp.push_back( fval );
						fval = atof( str_vec[5].c_str() ); crd_tmp.push_back( fval );
					}
					std::string snap_name = "OPT_GEOM_QCHEM_STD_ORT";
					if( !opt_geom )  snap_name = "NOPT_GEOM_QCHEM_STD_ORT";
					CrdSnapshot* p_snap = pmset->GetCrdSnapshotByName( snap_name.c_str(),true);
					if (p_snap == NULL ) throw std::runtime_error("Error to add optimized coordinate snapshot");
		
					ires = p_snap->SaveCrd(crd_tmp);
					if( !ires )  throw std::runtime_error("Error to save optimized coordinates to snapshot"); 
				
					p_snap->SetAtomCrd();
					read_mode0 = SEARCH_MODE;
				}
			}
		} // end for(;;) line reading cycle 
	}
	catch( const std::exception& ex)
	{
		PrintLog("Error in HaGaussMod::LoadOutFromStream() \n");
		PrintLog("On line: %s \n", line.c_str());
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return TRUE;
}

int HaGaussMod::LoadOutSummary( std::string summary_str )
{
//	PrintLog(" HaGaussMod::LoadOutSummary() pt 1 \n");
	boost::trim(summary_str);
	std::vector<std::string> str_vec;
	boost::split(str_vec,summary_str,boost::is_any_of("|"),boost::token_compress_on);

	std::vector<std::string> str_vec2;

	int i = 0;
	for( i = 0; i < str_vec.size(); i++ )
	{
		std::string str = str_vec[i];
		if( boost::istarts_with(str,"Polar=") )
		{
			PrintLog("Polar String: %s \n", str.c_str());
			std::string text = str.substr(6);
			boost::split(str_vec2,text,boost::is_any_of(","),boost::token_compress_on);
			if( str_vec2.size() != 6 )
			{
				PrintLog("Error in HaGaussMod::LoadOutSummary() \n");
				PrintLog("Can't read polarization info %s\n",str.c_str());
			}
			else
			{
				p_qc_mod->polar.resize(6);
				double cf = BOHR_TO_ANG*BOHR_TO_ANG*BOHR_TO_ANG;
				p_qc_mod->polar[0] = atof(str_vec2[0].c_str())*cf;
				p_qc_mod->polar[1] = atof(str_vec2[2].c_str())*cf;
				p_qc_mod->polar[2] = atof(str_vec2[5].c_str())*cf;
				p_qc_mod->polar[3] = atof(str_vec2[1].c_str())*cf;
				p_qc_mod->polar[4] = atof(str_vec2[3].c_str())*cf;
				p_qc_mod->polar[5] = atof(str_vec2[4].c_str())*cf;
			}
		}
		//if( boost::istarts_with(str,"Dipole=") )
		//{
		//	PrintLog("Dipole String: %s \n", str.c_str());
		//	boost::split(str_vec2,str.substr(7),boost::is_any_of(","));
		//	if( str_vec2.size() != 3 )
		//	{
		//		PrintLog("Error in HaGaussMod::LoadOutSummary() \n");
		//		PrintLog("Can't read dipole moment info %s\n",str.c_str());
		//	}
		//	else
		//	{
		//		double cf = BOHR_TO_ANG;
		//		p_qc_mod->dipole[0] = atof(str_vec2[0].c_str())*cf;
		//		p_qc_mod->dipole[1] = atof(str_vec2[1].c_str())*cf;
		//		p_qc_mod->dipole[2] = atof(str_vec2[2].c_str())*cf;
		//	}
		//}
	}
	return TRUE;
}

int HaGaussMod::RunFormChk(const char* fname_chk, const char* fname_fchk )
{
	GaussianProcess* p_gauss_proc = new GaussianProcess();
	p_gauss_proc->p_qc_mod = this->p_qc_mod;

	wxString cmd_line = "formchk ";
	cmd_line += fname_chk;
	cmd_line += " ";
	cmd_line += fname_fchk;

	PrintLog(" HaGaussMod::RunFormChk()  cmd_line: ");
	PrintLog(" %s \n", cmd_line.ToStdString().c_str() );

	int res = wxExecute(cmd_line,wxEXEC_SYNC,p_gauss_proc);

	PrintLog(" Finished converting .chk to .fchk file \n"); 
	return TRUE;
}

void HaGaussMod::FillSectionProcCommands(std::ostream& os) const
{
	if( n_sh_cores > 1 ) 
	{
		os << "%NProcShared=" << n_sh_cores << std::endl;
	}
	os << "%mem=" << max_mem << "MB" << std::endl;
	os << "%chk=" << GetCHKFileName() << std::endl;
}

void HaGaussMod::FillSectionJob(std::ostream& os) const
{
	if(p_qc_mod == NULL ) throw std::runtime_error("QChem Module is not set up");
	
	const MolSet* pmset = this->GetMolSet();	
	os << "#p ";

	std::string bas_name = p_qc_mod->GetBasName();
	int i;

	switch(p_qc_mod->wave_fun_type)
	{
		case harlem::qc::HARTREE_FOCK:
			os << "hf/";
			if( save_basis_set_gen ) 
			{
				os << "gen ";
			}
			else
			{
				os << bas_name << "  ";
			}
			os << "SCF=TIGHT ";
			break;
		case harlem::qc::NDO:
			if( p_qc_mod->ndo_method == harlem::qc::CNDO_2)
			{
				os << "CNDO ";
			}
			else if( p_qc_mod->ndo_method == harlem::qc::INDO_2)
			{
                os << "INDO ";
			}
			else
			{
                os << "INDO ";
			}
			break;
		case harlem::qc::MP2:
			os << " MP2/" ;
			if( save_basis_set_gen ) 
			{
				os << "gen ";
			}
			else
			{
				os << bas_name << "  ";
			}
			
			break;
		case harlem::qc::DFT_B3LYP:
			os << "B3LYP/";
			if( save_basis_set_gen ) 
			{
				os << "gen ";
			}
			else
			{
				os << bas_name << "  ";
			}
			break;
		case harlem::qc::CCSD_T:
			os << "CCSD(T)/";
			if( save_basis_set_gen ) 
			{
				os << "gen ";
			}
			else
			{
				os << bas_name << "  ";
			}
			break;
		default: 
			break;
	}

	if( read_init_geom_chk_file )
	{
		os << "geom=checkpoint ";
	}

	if( read_hf_guess_chk_file )
	{
		os << "guess=read ";
	}	

	if( p_qc_mod->IsEneMinCalc() )
	{
		ZMatCrd* pzm = p_qc_mod->GetZMat();
		if( pzm && !pzm->IsEmpty() ) 
		{
			os << " opt=z-mat ";
		}
		else
		{
			os << " opt ";
		}
	}

	if( p_qc_mod->ToCalcPolar() )
	{
		os << "polar "; 
	}

	for( i = 0; i < 3; i++ )
	{
		if( fabs( p_qc_mod->el_field[i] ) > 0.0000001 )
		{
			int ifld = p_qc_mod->el_field[i]*10000*BOHR_TO_ANG*BOHR_TO_ANG; 
			os << "Field=";
			if( i == 0 ) os << "X";
			if( i == 1 ) os << "Y";
			if( i == 2 ) os << "Z";
			
			if( ifld < 0 )
			{
				os << "-";
				ifld = -ifld;
			}
			else
			{
				os << "+";
			}
			os << ifld << " ";
		}
	}

	if( no_std_orient )
	{
		os << "NOSYMM ";
	}

	if( p_qc_mod->set_guess_from_mos && p_qc_mod->MO_coef.num_cols() > 0 )
	{
		os << "guess=cards ";
	}
    
	const HaMolecule* pmol_ext_ch = pmset->GetMolByName("EXTERNAL_CHARGES");
	if( pmol_ext_ch != NULL || !p_qc_mod->ext_chrg.empty() )
	{
		os << " CHARGE ";
	}
	
	if( !p_qc_mod->IsUsingSymmetry() ) os << "NOSYMM ";
	
	if(p_qc_mod->UsePseudoPot()) os << " pseudo=cards ";

	if( !add_kw_str.empty() ) 
	{
		os << " " << add_kw_str;
	}
	
	os << std::endl;
	os << "  " << std::endl;
	os << " This is Gaussian file created by HARLEM " << std::endl;
	os << "  " << std::endl;
}


void HaGaussMod::FillSectionCoord(std::ostream& os) const
{
	char buf[120];
	MolSet* phmol_set= p_qc_mod->GetMolSet();

	sprintf(buf,"  %d %d ",p_qc_mod->GetCharge(),p_qc_mod->GetMult());
	os << buf << std::endl;

	if( read_init_geom_chk_file )
	{
		os << "  " << std::endl;
		return;
	}

	ZMatCrd* p_zmat = p_qc_mod->GetZMat();

	if( p_zmat->IsEmpty() )
	{
		HaAtom* aptr;
		AtomIteratorMolSet aitr(phmol_set);
		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			HaMolecule* pmol = aptr->GetHostMol();
			if( strcmp(pmol->GetObjName(),"EXTERNAL_CHARGES") == 0) continue;

			os << " " << aptr->GetStdSymbol().c_str();
			sprintf(buf," %20.15f  %20.15f  %20.15f ",aptr->GetX_Ang(),aptr->GetY_Ang(),aptr->GetZ_Ang());
			os << buf << std::endl;  
		}
	}
	else
	{
		os << p_zmat->SaveToString();
	}
	os << "  " << std::endl;
}

void HaGaussMod::FillSectionBasis(std::ostream& os) const
{
	if(p_qc_mod->wave_fun_type == harlem::qc::NDO) return;
	if( !p_qc_mod->IsGenBasisSet() && !save_basis_set_gen ) return;
	
	std::set<GauAtomBasis> bs_set;
	std::set<const HaPseudoPot*> ppot_set;
	AtBasisType::iterator bitr;

	for(bitr= p_qc_mod->AtBasis.at_bas_vec.begin(); bitr != p_qc_mod->AtBasis.at_bas_vec.end() ; bitr++)
	{
		if(bs_set.find(*bitr) == bs_set.end()) // atom basis is not encountered yet
		{			
			bs_set.insert(*bitr);
			(*bitr).SaveGaussianInp(os);
			os << "****" << std::endl;
		}
		if((*bitr).IsSetPseudoPot())
		{
			const HaPseudoPot* ppot= (*bitr).GetPseudoPot();
			if(ppot_set.find(ppot) == ppot_set.end())
			{
				ppot_set.insert(ppot);
			}
		}
	}
	if(!ppot_set.empty())
	{
		os <<  "    " << std::endl;

		std::set<const HaPseudoPot*>::iterator pitr;
		for(pitr=ppot_set.begin(); pitr != ppot_set.end(); pitr++)
		{
			(*pitr)->SaveGaussInp(os);
		}
	}
	os <<  "    " << std::endl;
}

void HaGaussMod::FillSectionExtCharges(std::ostream& os) const
{
	char buf[120];
	MolSet* pmset = p_qc_mod->GetMolSet();
	const HaMolecule* pmol_ext_ch = pmset->GetMolByName("EXTERNAL_CHARGES");

	if( pmol_ext_ch == NULL && p_qc_mod->ext_chrg.empty() ) return;

	if( fabs(p_qc_mod->distr_ext_chrg) > 0.0001 ) // distributed charge as a dot surface
	{
		DotStruct* dot_surf = NULL;
		std::list<Object3D*>::iterator oitr;
		for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end();)
		{
			if( (*oitr)->GetObjType() == OBJ3D_DOT_SURFACE )
			{
				dot_surf = (DotStruct*) (*oitr);
				break;
			}
			oitr++;
		}
		if( dot_surf == NULL ) return; 
			
		int save_dot_ch = TRUE;
		int np = 0;
			
		np = dot_surf->GetCount();
		if(np == 0) return;
			
		double ch = p_qc_mod->distr_ext_chrg/np;
				
		std::vector<HaDot>::iterator ditr;
		for(ditr= dot_surf->dots.begin(); ditr != dot_surf->dots.end(); ditr++)
		{
			HaDot& cur_dot = *ditr;
			sprintf(buf,"%16.9f %16.9f %16.9f %16.9f ",cur_dot.GetX(),cur_dot.GetY(),cur_dot.GetZ(),ch);
			os << buf << std::endl;
		}
	}
	else if( pmol_ext_ch != NULL )
	{
        AtomIteratorMolecule_const aitr(pmol_ext_ch);
		const HaAtom* aptr;
		for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom() )
		{
			double ch = aptr->GetCharge();
			sprintf(buf,"%16.9f %16.9f %16.9f %16.9f ",aptr->GetX() + 0.01,aptr->GetY() + 0.01,aptr->GetZ() + 0.01, ch);
			os << buf << std::endl;
		}
	}
	AtomDoubleMap::iterator ch_itr = p_qc_mod->ext_chrg.begin();

	for( ; ch_itr != p_qc_mod->ext_chrg.end(); ch_itr++ )
	{
		HaAtom* aptr = (*ch_itr).first;
		double ch    = (*ch_itr).second;

		double crd_shift = p_qc_mod->ext_ch_crd_offset;
		sprintf(buf,"%16.9f %16.9f %16.9f %16.9f ",aptr->GetX() + crd_shift,aptr->GetY() + crd_shift,aptr->GetZ() + crd_shift,ch);
		os << buf << std::endl;
	}
	os << " " << std::endl;
}

void HaGaussMod::FillSectionInitMO(std::ostream& os) const
{
	if( !p_qc_mod->set_guess_from_mos || p_qc_mod->MO_coef.num_cols() == 0 ) return;
	char buf[120];

	os << "(5E16.8)" << std::endl;
	int nc = p_qc_mod->MO_coef.num_cols();
	int nb = p_qc_mod->MO_coef.num_rows();
	if( p_qc_mod->AtBasis.GetNBfunc() != nb ) 
		throw std::runtime_error(" Number of Basis Set functions do not equal to the number of rows in MO_coef matrix used as an initial guess");  
	int i;
	for(i = 1; i <= nc; i++)
	{
		sprintf(buf,"%5d",i);
		os << buf << std::endl;
		HaVec_double cf;
		cf.set_ext_alloc( &p_qc_mod->MO_coef.r1(1,i),nb ); 
		write_double_array_chuncks(os,cf,5,"%16.8e");
	}
	os << "  " << std::endl;
}

void HaGaussMod::PrintCurBcommon()
{
#if defined(GAUSSVER)

	char buf[120];
	printf("\n Common block B: \n");
	printf(" Description of the basis set \n\n");
	printf(" the number of primitive shells : nshell = %d \n",b_.nshell); 
	printf(" the highest angular q.n. in b.s: maxtyp = %d \n",b_.maxtyp);
	printf("\n Characteristics of the shells: \n");
	printf("\n   N    AON      AOS SHELLA SHELLT SHELLC SHELLN   X         Y         Z \n");  
	
	int i;
	for( i=0; i < b_.nshell ; i++)
	{
		printf( "%4d %8x  %4d  %4d  %4d  %4d  %4d %10.5f %10.5f %10.5f \n",
			i+1, 
			b_.aon[i],
			b_.aos[i],
			b_.shella[i],
			b_.shellt[i],
			b_.shellc[i],
			b_.shelln[i],
			b_.x[i],
			b_.y[i],
			b_.z[i]);
	}
	cout << endl;	
	
	cout << " Gaussian exponents and expansion coef. of S and P orb.: " << endl;	
	
	for( i=0; i < b_.nshell ; i++)
	{
		cout << " Shell # " << (i+1) << endl;
		cout << " Gaussian exponents and expansion coef. of S and P orb.: " << endl;
		int ng= b_.shelln[i];
		int ib_sp=b_.shella[i]-1;
		for(int j=0; j < ng; j++)
		{
			sprintf(buf," %16.9f  %16.9f  %16.9f ", b_.exx[ib_sp],b_.c1[ib_sp],b_.c2[ib_sp]);
			cout << buf << endl;
			ib_sp++;
		}
		cout << endl;
	}	
#endif
}

//
// For Symmetrical matricies in Gaussian (as Density matrix for example):
//
// double & 
// GauMat1e::operator()(const int i, const int j)
// 1 element access
// {
//	if(smode == SMODE_SYMM)
//		// get index of low triangle stored matrix
//		return( (i >= j) ? *(Data.begin() + i*(i-1)/2 + j-1) 
//	                	 : *(Data.begin() + j*(j-1)/2 + i-1)  );
//	return( Data(i,j));
// }

//GauMat1e & 
//GauMat1e::operator=(const HaMat_double& A)
//{
//	nbasis=A.num_cols();
//	
//	int k=0;
//	if(smode==SMODE_SYMM)
//	{
//		for(int i=1; i <= nbasis; i++)
//		{
//			for(int j=1; j <= i; j++)
//			{
//				*(Data.begin() + k)=A(i,j);
//				k++;
//			}
//		}
//	}
//	else
//		Data=A;
//	return *this;
//}


void HaGaussMod::SetNumSharedMemCores(int n_sh_cores_new)
{
	n_sh_cores = n_sh_cores_new; 
}

void HaGaussMod::SetNumProc(int n_proc_new)
{
	n_proc = n_proc_new;
}

void HaGaussMod::SetMaxMem(int max_mem_new )
{
	max_mem = max_mem_new;
}

int  HaGaussMod::GetNumSharedMemCores() const
{
	return n_sh_cores;
}

int  HaGaussMod::GetNumProc() const
{
	return n_proc;
}

int  HaGaussMod::GetMaxMem()  const
{
	return max_mem;
}	

void  HaGaussMod::SetLoadNonOptGeom( bool set_par )
{
	load_non_opt_geom = set_par;
}

void HaGaussMod::SetLoadGeomZMatOrient( bool set_par )
{
	load_geom_zmat_orient = set_par;
}
	
void HaGaussMod::SetLoadGeomStdOrient ( bool set_par )
{
	load_geom_std_orient = set_par;
}

void HaGaussMod::SetReadInitGeomChkFile( bool set_par )
{
	read_init_geom_chk_file = set_par;
}

void HaGaussMod::SetReadHFGuessChkFile( bool set_par )
{
	read_hf_guess_chk_file = set_par;
}

void HaGaussMod::SetNoStdOrient( bool set_par )
{
	no_std_orient = set_par;
}

void HaGaussMod::SetSaveBasisSetGen( bool set_par  )
{
	save_basis_set_gen = set_par;
}

void HaGaussMod::SetAddKWStr( const std::string& add_kw_str_par )
{
	add_kw_str = add_kw_str_par;
}
