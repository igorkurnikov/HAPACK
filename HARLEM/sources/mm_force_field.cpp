/*! \file mm_force_field.cpp

    Classes to describe Molecular Mechanics force field in HARLEM
 
    \author Igor Kurnikov 
    \date 2010-
*/

#define MM_FORCE_FIELD_CPP

#include "mpi.h"

#include <math.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <wx/filename.h>

#include "tinyxml.h"

#include "harlemapp.h"
#include "haresdb.h"
#include "hamolmech.h"
#include "mm_force_field.h"
#include "mm_elements.h"

#include "ha_mort_mm.h"

ForceFieldType MMForceField::ff_type_default = ForceFieldType::AMBER_10;
std::vector<MMForceField*> MMForceField::ff_arr;

StrVec  MMForceField::resff_files_add;        
StrVec  MMForceField::tinker_param_files_add; 
StrVec  MMForceField::amber_param_files_add;

MMForceField::MMForceField()
{
	p_mort_ff = NULL;
	Clear();
}
	
MMForceField::~MMForceField()
{
	Clear();
}

int MMForceField::SetFFType( const ForceFieldType& ff_type_par)
{
	ff_type = ff_type_par;
	Clear();
	SetDefaultVdW();
	SetDefaultParamFiles();
	return TRUE;
}

ForceFieldType MMForceField::GetFFType() const
{
	return ff_type;
}

void MMForceField::Clear()
{
	symb_mass_map.clear();
	symb_ppar_map.clear();
	symb_hpar_map.clear();
	bond_param_map.clear();    
	vang_param_map.clear();
	dih_param_map.clear();
	impdih_param_map.clear(); 

	if( p_mort_ff != NULL) delete p_mort_ff;
	if( res_ff_templates.size() > 0 )
	{
		int n = res_ff_templates.size();
		for( int i = 0; i < n; i++ )
		{
			delete res_ff_templates[i];
		}
		res_ff_templates.clear();
		res_name_ff_templ_map.clear();
	}
}

void MMForceField::SetDefaultVdW()
{
	HaVec_double ppar(2);
	ppar[0] = 1.908;
	ppar[1] = 0.086;
	symb_ppar_map["XXST"] = ppar; //!< Set Default VdW parameters

	ppar[0] = 0.0;
	ppar[1] = 0.0;
	symb_ppar_map["DC"] = ppar; //!< Set Default VdW parameters for Dummy Atoms
	symb_ppar_map["DN"] = ppar; //!< Set Default VdW parameters for Dummy Atoms
	symb_ppar_map["DO"] = ppar; //!< Set Default VdW parameters for Dummy Atoms
	symb_ppar_map["DH"] = ppar; //!< Set Default VdW parameters for Dummy Atoms
}

int MMForceField::IsMortFFInitiated()
{
	if( p_mort_ff == NULL) return FALSE;
	return TRUE;
}
 
int MMForceField::Init()
{
	int ires; 
	std::string file_name;

	using namespace boost::filesystem;
 
//	cout << "MMForceField::Init() pt 1 \n" << "Current PATH:" << current_path() << std::endl;

	directory_iterator ditr( current_path() );
	for( ; ditr != directory_iterator(); ditr++)
	{
		file_name = ditr->path().filename().string();
		// PrintLog(" file_name = %s \n",file_name.c_str() );
		if( boost::starts_with(file_name,"frcmod") )
		{
			amber_param_files.push_back(file_name);
			PrintLog(" Add AMBER parameter file %s in the working directory \n", file_name.c_str());
		}
		else if( boost::starts_with(file_name,"resff") ) 
		{
			resff_files.push_back( file_name );
			PrintLog(" Add Residue Force Field file %s in the working directory \n", file_name.c_str());
		}
		else if( boost::starts_with(file_name,"tinker_param") ) 
		{
			tinker_param_files.push_back( file_name );
			PrintLog(" Add Tinker Force Field Parameter file %s in the working directory \n", file_name.c_str());
		}
	}

	try
	{
		int i;
		for(i = 0; i < amber_param_files.size(); i++)
		{
			file_name = amber_param_files[i];
			PrintLog(" Force field parameter file (AMBER FORMAT) = %s \n",file_name.c_str());

			bool file_exists = ::wxFileExists(file_name.c_str());
			if( !file_exists ) throw std::runtime_error(" Parameter file " + file_name + " doesn't exist ");
				
			ires = LoadAmberParamFile(file_name);
			if( !ires ) throw std::runtime_error(" Error loading Amber Parameter file " + file_name );
		}

		for(i = 0; i < tinker_param_files.size(); i++)
		{
			file_name = tinker_param_files[i];
			PrintLog(" Force field parameter file (TINKER FORMAT) = %s \n",file_name.c_str());

			bool file_exists = ::wxFileExists(file_name.c_str());
			if( !file_exists ) throw std::runtime_error(" Parameter file " + file_name + " doesn't exist ");
				
			ires = LoadTinkerParamFile(file_name);
			if( !ires ) throw std::runtime_error(" Error loading Amber Parameter file " + file_name );
		}

		for(i = 0; i < resff_files.size(); i++)
		{
			file_name = resff_files[i];
			PrintLog(" Residue Force field file = %s \n",file_name.c_str());

			bool file_exists = ::wxFileExists(file_name.c_str());
			if( !file_exists ) throw std::runtime_error(" Residue FF file " + file_name + " doesn't exist ");
				
			ires = LoadResFFTemplateXMLFile( file_name.c_str() );
			if( !ires ) throw std::runtime_error(" Error loading Residue FF file " + file_name );
		}

//		ires = InitMortFF();
//		if( !ires ) throw std::runtime_error( (std::string)" Error to initialize MORT structures for force field " + ff_type.label() );
	}
	catch( std::exception& ex )
	{
		PrintLog("Error in MMForceField::Init() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}

	return TRUE;
}

MMForceField* MMForceField::GetMMForceField(const ForceFieldType& ff_type_par, int create )
{
	ForceFieldType ff_type_loc = ff_type_par;
	if( ff_type_loc == ForceFieldType::UNKNOWN_FF ) ff_type_loc = MMForceField::ff_type_default;
	
	int i;
	for(i = 0; i < ff_arr.size(); i++)
	{
		if(ff_arr[i]->GetFFType() == ff_type_loc) return ff_arr[i];
	}

	if(!create) return NULL;

	MMForceField* p_ff = new MMForceField();
	p_ff->SetFFType(ff_type_loc);

	int ires = p_ff->Init();
	if( !ires )
	{
		delete p_ff;
		p_ff = NULL;
	}
	if( p_ff ) ff_arr.push_back(p_ff);
	return p_ff;
}

double MMForceField::FindAtomMassFromSymbol( const std::string& ats )
{
	if( symb_mass_map.count(ats) == 0 ) return -1.0;
	return symb_mass_map[ats];
}

HaVec_double MMForceField::FindPointParamFromSymbol(const char* ff_symb)
{
	std::string at_symbol = ff_symb;

	HaVec_double ppar;
	map<std::string,HaVec_double,less<std::string> >::iterator pitr;

	pitr = symb_ppar_map.find(at_symbol);
	if(pitr != symb_ppar_map.end() )
	{
		ppar = (*pitr).second;
	}

	return ppar;	
}

static int StrMatch(const std::string& str, const std::string& templ_str)
{
	if( templ_str == "X" ) 
		return 1;
	if( str == 	templ_str) 
		return 2;
	return 0;
}

HaVec_double MMForceField::FindBondParamFromSymbol(const char* as1,const char* as2)
{
	std::string at_symbol_1 = as1;
	std::string at_symbol_2 = as2;
	
	if(at_symbol_2 < at_symbol_1) at_symbol_2.swap(at_symbol_1);

	boost::trim(at_symbol_1); boost::to_upper(at_symbol_1);
	boost::trim(at_symbol_2); boost::to_upper(at_symbol_2);

	HaVec_double best_match;

	std::string label;
	map<std::string,HaVec_double,less<std::string> >::iterator mitr;

	int ifound = FALSE;
    label = at_symbol_1 + "-" + at_symbol_2;
	mitr = bond_param_map.find(label);
	if(mitr != bond_param_map.end())
	{
		ifound = TRUE;
		best_match = (*mitr).second;
	}

	std::string ats_1;
	std::string ats_2;

	if(!ifound) 
	{
		ats_1 = "X";
		ats_2 = at_symbol_2;
		if( ats_2 < ats_1 )
		{
			ats_1 = at_symbol_2;
			ats_2 = "X";
		}
		label = ats_1 + "-" + ats_2;
		mitr = bond_param_map.find(label);
	    if(mitr != bond_param_map.end())
		{
			ifound = TRUE;
	       	best_match = (*mitr).second;
		}
	}

	if(!ifound) 
	{
		ats_1 = "X";
		ats_2 = at_symbol_1;
		if( ats_2 < ats_1 )
		{
			ats_1 = at_symbol_1;
			ats_2 = "X";
		}
		label = ats_1 + "-" + ats_2;
		mitr = bond_param_map.find(label);
	    if(mitr != bond_param_map.end())
		{
			ifound = TRUE;
	       	best_match = (*mitr).second;
		}
	}

	return best_match;
}

HaVec_double MMForceField::FindHBondParamFromSymbol(const char* as1,const char* as2)
{
	std::string at_symbol_1 = as1;
	std::string at_symbol_2 = as2;
	if(at_symbol_2 < at_symbol_1) at_symbol_2.swap(at_symbol_1);

	boost::trim(at_symbol_1); boost::to_upper(at_symbol_1);
	boost::trim(at_symbol_2); boost::to_upper(at_symbol_2);

	HaVec_double best_match;

	std::string label;
	map<std::string,HaVec_double,less<std::string> >::iterator mitr;

	int ifound = FALSE;
    label = at_symbol_1 + "-" + at_symbol_2;
	mitr = symb_hpar_map.find(label);
	if(mitr != symb_hpar_map.end())
	{
		ifound = TRUE;
		best_match = (*mitr).second;
	}

return best_match;
}

HaVec_double MMForceField::FindValAngleParamFromSymbol(const char* ats1,const char* ats2,const char* ats3)
{
	std::string at_symbol_1 = ats1;
	std::string at_symbol_2 = ats2;
	std::string at_symbol_3 = ats3;

	if(at_symbol_3 < at_symbol_1) at_symbol_3.swap(at_symbol_1);

	boost::trim(at_symbol_1); boost::to_upper(at_symbol_1);
	boost::trim(at_symbol_2); boost::to_upper(at_symbol_2);
	boost::trim(at_symbol_3); boost::to_upper(at_symbol_3);
	
	HaVec_double best_match;

	std::string label;
	map<std::string,HaVec_double,less<std::string> >::iterator mitr;
	
	int ifound = FALSE;

	std::string ats_1 = at_symbol_1;
	std::string ats_2 = at_symbol_2;
	std::string ats_3 = at_symbol_3;

	label = ats_1 + "-" + ats_2 + "-" + ats_3;
	mitr = vang_param_map.find(label);
	if(mitr != vang_param_map.end())
	{
		best_match = (*mitr).second;
		ifound = TRUE;
	}

	if(!ifound) 
	{
		ats_1 = "X";
		ats_3 = at_symbol_3;
		if( ats_3 < ats_1 )
		{
			ats_1 = at_symbol_3;
			ats_3 = "X";
		}
		label = ats_1 + "-" + ats_2 + "-" + ats_3;
	    mitr = vang_param_map.find(label);
	    if(mitr != vang_param_map.end())
		{
		    best_match = (*mitr).second;
			ifound = TRUE;
		}
	}

	if(!ifound) 
	{
		ats_1 = "X";
		ats_3 = at_symbol_1;
		if( ats_3 < ats_1 )
		{
			ats_1 = at_symbol_1;
			ats_3 = "X";
		}
		label = ats_1 + "-" + ats_2 + "-" + ats_3;
	    mitr = vang_param_map.find(label);
	    if(mitr != vang_param_map.end())
		    best_match = (*mitr).second;
	}

	return best_match;

}

HaVec_double MMForceField::FindDihedralParamFromSymbol(const char* as1,const char* as2,const char* as3,const char* as4, bool improper_flag)
{
	std::string at_symbol_1 = as1;
	std::string at_symbol_2 = as2;
	std::string at_symbol_3 = as3;
	std::string at_symbol_4 = as4;
	
	if( (at_symbol_4 < at_symbol_1)  || 
		(( at_symbol_4 == at_symbol_1) &&  (at_symbol_2 > at_symbol_3))   )
	{
		at_symbol_4.swap(at_symbol_1);
		at_symbol_3.swap(at_symbol_2);
	}

	boost::trim(at_symbol_1); boost::to_upper(at_symbol_1);
	boost::trim(at_symbol_2); boost::to_upper(at_symbol_2);
	boost::trim(at_symbol_3); boost::to_upper(at_symbol_3);
	boost::trim(at_symbol_4); boost::to_upper(at_symbol_4);
	
	std::string ats_1;
	std::string ats_2;
	std::string ats_3;
	std::string ats_4;

	HaVec_double best_match;

	map<std::string,HaVec_double,less<std::string> >  *ptr_dih_map;

	if(improper_flag)
		ptr_dih_map = &impdih_param_map;
	else
		ptr_dih_map = &dih_param_map;

	std::string label;
	map<std::string,HaVec_double,less<std::string> >::iterator mitr;

	int ifound = FALSE;
	label = at_symbol_1 + "-" + at_symbol_2 + "-" + at_symbol_3 + "-" + at_symbol_4;
	mitr = ptr_dih_map->find(label);
	if(mitr != ptr_dih_map->end())
	{
		ifound = TRUE;
		best_match = (*mitr).second;
	}

	if(!ifound) 
	{
		ats_1 = "X";
		ats_2 = at_symbol_2;
		ats_3 = at_symbol_3;
		ats_4 = at_symbol_4;

		if( (ats_4 < ats_1) || ((ats_4 == ats_1) && (at_symbol_3 < at_symbol_2)) )
		{
			ats_1 = at_symbol_4;
			ats_2 = at_symbol_3;
			ats_3 = at_symbol_2;
			ats_4 = "X";
		}

		label = ats_1 + "-" + ats_2 + "-" + ats_3 + "-" + ats_4;
	    mitr = ptr_dih_map->find(label);
	    if(mitr != ptr_dih_map->end())
		{
			ifound = TRUE;
		    best_match = (*mitr).second;
		}
	}

	if(!ifound) 
	{
		ats_1 = at_symbol_1;
		ats_2 = at_symbol_2;
		ats_3 = at_symbol_3;
		ats_4 = "X";

		if( (ats_4 < ats_1) || ((ats_4 == ats_1) && (at_symbol_3 < at_symbol_2)) )
		{
			ats_1 = "X";
			ats_2 = at_symbol_3;
			ats_3 = at_symbol_2;
			ats_4 = at_symbol_1;
		}

		label = ats_1 + "-" + ats_2 + "-" + ats_3 + "-" + ats_4;
	    mitr = ptr_dih_map->find(label);
	    if(mitr != ptr_dih_map->end())
		{
			ifound = TRUE;
		    best_match = (*mitr).second;
		}
	}

	if(!ifound) 
	{
		ats_1 = "X";
		ats_2 = at_symbol_2;
		ats_3 = at_symbol_3;
		ats_4 = "X";

		if( (at_symbol_3 < at_symbol_2) )
		{
			ats_2 = at_symbol_3;
			ats_3 = at_symbol_2;
		}

		label = ats_1 + "-" + ats_2 + "-" + ats_3 + "-" + ats_4;
	    mitr = ptr_dih_map->find(label);
	    if(mitr != ptr_dih_map->end())
		{
			ifound = TRUE;
		    best_match = (*mitr).second;
		}
	}

	if(!ifound) 
	{
		ats_1 = "X";
		ats_2 = "X";
		ats_3 = at_symbol_3;
		ats_4 = at_symbol_4;

		if( (ats_4 < ats_1) || ((ats_4 == ats_1) && (ats_3 < ats_2)) )
		{
			ats_1 = at_symbol_4;
			ats_2 = at_symbol_3;
			ats_3 = "X";
			ats_4 = "X";
		}

		label = ats_1 + "-" + ats_2 + "-" + ats_3 + "-" + ats_4;
	    mitr = ptr_dih_map->find(label);
	    if(mitr != ptr_dih_map->end())
		{
			ifound = TRUE;
		    best_match = (*mitr).second;
		}
	}

	if(!ifound) 
	{
		std::string ats_1 = "X";
		std::string ats_2 = "X";
		std::string ats_3 = at_symbol_2;
		std::string ats_4 = at_symbol_1;

		if( (ats_4 < ats_1) || ((ats_4 == ats_1) && (ats_3 < ats_2)) )
		{
			ats_1 = at_symbol_1;
			ats_2 = at_symbol_2;
			ats_3 = "X";
			ats_4 = "X";
		}

		label = ats_1 + "-" + ats_2 + "-" + ats_3 + "-" + ats_4;
	    mitr = ptr_dih_map->find(label);
	    if(mitr != ptr_dih_map->end())
		{
			ifound = TRUE;
		    best_match = (*mitr).second;
		}
	}

	return best_match;
}

const int ASCII_SPACE = 32;
const int ASCII_MINUS = 45;

static void get_str_from_stream(istream& is,std::string& str)
{
	int ich;
	str.clear();
	int in_str = FALSE;
	for(;;)
	{
		ich = is.get();
		if(is.fail()) return;
		if(in_str)
		{
			if(isspace(ich) || ich == ASCII_MINUS) return;
			str += (char) ich;
		}
		else
		{
			if(isspace(ich) || ich == ASCII_MINUS) continue;
			str += (char) ich;
			in_str = TRUE;
		}
	}
}

int MMForceField::LoadAmberParamFile(const std::string& ff_param_fname )
{
	char buf[256]; 
	std::string str_inp;
	StrVec tokens;

	try
	{
		std::ifstream is_f(ff_param_fname.c_str());
		if(is_f.fail()) throw std::runtime_error(" Can't open file " + ff_param_fname );
		
		enum READ_MODE { READ_ATOM_MASSES, READ_BOND_PARAMS,  
			READ_ANGLE_PARAMS, READ_DIHEDRAL_PARAMS, 
			READ_IMPR_DIHEDRAL_PARAMS, READ_WATER_PARAMS,
			READ_VDW_ATOM_SYNONYMS,
			READ_VDW_PARAMS,READ_H_BOND_PARAMS} read_mode;

		std::getline(is_f,str_inp);
		if( is_f.fail() ) throw std::runtime_error("Error reading first line in file " + ff_param_fname );
	
		read_mode = READ_ATOM_MASSES;

		multimap<std::string, std::string, less<std::string> >vdw_synonyms;
		int i;

		for(;;)
		{
			std::getline(is_f,str_inp);
			boost::trim(str_inp);

			if(is_f.eof() ) return TRUE;
			if(is_f.fail()) throw std::runtime_error("Error Reading line from file" + ff_param_fname );
			
			if( str_inp == "MASS")             { read_mode = READ_ATOM_MASSES; continue; }
			if( str_inp == "BOND")             { read_mode = READ_BOND_PARAMS; continue; }
			if( str_inp.substr(0,4) == "ANGL") { read_mode = READ_ANGLE_PARAMS; continue; }
			if( str_inp.substr(0,4) == "DIHE") { read_mode = READ_DIHEDRAL_PARAMS; continue; }
			if( str_inp.substr(0,4) == "IMPR") { read_mode = READ_IMPR_DIHEDRAL_PARAMS; continue; }
			if( str_inp.substr(0,4) == "NONB") { read_mode = READ_VDW_PARAMS; continue; }

			std::istringstream is(str_inp);

			if( read_mode == READ_ATOM_MASSES)
			{
				if( str_inp.empty() )
				{
					read_mode = READ_BOND_PARAMS;
					continue;
				}

				std::string ats;
				get_str_from_stream(is,ats); if( is.fail()) throw std::runtime_error("Error Reading Atom Mass Params");

				boost::to_upper(ats);
				
				double amass;

				is >> amass;
				if(is.fail()) throw std::runtime_error("Error Reading Atom Mass Params");

				symb_mass_map[ats] = amass;
			}
			else if( read_mode == READ_BOND_PARAMS)
			{
				if( str_inp.empty() )
				{
					read_mode = READ_ANGLE_PARAMS;
					continue;
				}

				if( str_inp.find("-") == std::string::npos) continue; // skip line   C   H   HO  N   NA  NB  NC  N2  NT  N2  N3  N*  O   OH  OS  P   O2  what is it?

				std::string at1s,at2s;

				get_str_from_stream(is,at1s); if( is.fail()) throw std::runtime_error("Error Reading Bond Params");
				get_str_from_stream(is,at2s); if( is.fail()) throw std::runtime_error("Error Reading Bond Params");

				boost::to_upper(at1s);
				boost::to_upper(at2s);

				if(at2s < at1s) at1s.swap(at2s);
				
				double fconst;
				double eq_dist;

				is >> fconst;
				is >> eq_dist;
				if(is.fail()) throw std::runtime_error("Error Reading Bond Params");

				HaVec_double bpar(2);
				bpar[0] = eq_dist;
				bpar[1] = fconst;
				std::string label = at1s + "-" + at2s;

				bond_param_map[label] = bpar;
			}
			else if( read_mode == READ_ANGLE_PARAMS)
			{
				if( str_inp.empty() )
				{
					read_mode = READ_DIHEDRAL_PARAMS;
					continue;
				}

				std::string at1s,at2s,at3s;

				get_str_from_stream(is,at1s); if( is.fail()) throw std::runtime_error("Error Reading Valence Angle Params");
				get_str_from_stream(is,at2s); if( is.fail()) throw std::runtime_error("Error Reading Valence Angle Params");
				get_str_from_stream(is,at3s); if( is.fail()) throw std::runtime_error("Error Reading Valence Angle Params");

				boost::to_upper(at1s);
				boost::to_upper(at2s);
				boost::to_upper(at3s);

				if(at3s < at1s) at1s.swap(at3s);
			
				double fconst;
				double eq_ang;

				is >> fconst;
				is >> eq_ang;
				if(is.fail()) throw std::runtime_error("Error Reading Valence Angle Params");
				
				HaVec_double vpar(2);
				vpar[0] = eq_ang;
				vpar[1] = fconst;
				std::string label = at1s + "-" + at2s + "-" + at3s;

				vang_param_map[label] = vpar;
			}
			else if( read_mode == READ_DIHEDRAL_PARAMS)
			{
				if( str_inp.empty() )
				{
					read_mode = READ_IMPR_DIHEDRAL_PARAMS;
					continue;
				}

				std::string at1s,at2s,at3s,at4s;

				get_str_from_stream(is,at1s); if( is.fail()) throw std::runtime_error("Error Reading Dihedral Angle Params");
				get_str_from_stream(is,at2s); if( is.fail()) throw std::runtime_error("Error Reading Dihedral Angle Params");
				get_str_from_stream(is,at3s); if( is.fail()) throw std::runtime_error("Error Reading Dihedral Angle Params");
				get_str_from_stream(is,at4s); if( is.fail()) throw std::runtime_error("Error Reading Dihedral Angle Params");

				boost::to_upper(at1s);
				boost::to_upper(at2s);
				boost::to_upper(at3s);
				boost::to_upper(at4s);

				if( (at1s > at4s) || ( (at1s == at4s) &&  (at2s > at3s) ) )
				{
					at1s.swap(at4s);
					at2s.swap(at3s);
				}

				int idiv_fc;
				double fconst;
				double ph;
				double nperiod;

				is >> idiv_fc;
				is >> fconst;
				is >> ph;
				is >> nperiod;

				if(is.fail()) throw std::runtime_error("Error Reading Dihedral Angle Params");
				
				HaVec_double tpl(100);

				int nt = 1;
				int idx = 0;
				tpl[idx] = fabs(nperiod);idx++;
				tpl[idx] = ph;           idx++;
				tpl[idx] = fconst;       idx++;
				tpl[idx] = idiv_fc;      idx++;

//				std::string at1s_2,at2s_2,at3s_2,at4s_2;

				while( nperiod  < 0 ) 
				{
					std::getline(is_f,str_inp);
					std::istringstream is2(str_inp.substr(11));
					is2 >> idiv_fc;
					is2 >> fconst;
					is2 >> ph;
					is2 >> nperiod;

					if(is2.fail()) throw std::runtime_error("Error Reading Dihedral Angle Params");

					nt++;
					tpl[idx] = fabs(nperiod);idx++;
					tpl[idx] = ph;           idx++;
					tpl[idx] = fconst;       idx++;
					tpl[idx] = idiv_fc;      idx++;
				}

				HaVec_double tpar(idx);
				for(int i = 0; i< idx; i++)
				{
					tpar[i] = tpl[i];
				}

				std::string label = at1s + "-" + at2s + "-" + at3s + "-" + at4s;

				//if( label == "H-N-C-O")
				//{
				//	PrintLog(" Dihedral Params for H-N-C-O: \n");
				//	int i;
				//	for( i = 0; i < tpar.size(); i++)
				//	{
				//		PrintLog(" %5i  %12.6f \n",i,tpar[i]);
				//	}
				//}

				dih_param_map[label] = tpar;
			}
			else if( read_mode == READ_IMPR_DIHEDRAL_PARAMS)
			{
				if( str_inp.empty() )
				{
					read_mode = READ_WATER_PARAMS;
					continue;
				}

				std::string at1s,at2s,at3s,at4s;

				get_str_from_stream(is,at1s); if( is.fail()) throw std::runtime_error("Error Reading Improper Dihedral Angle Params");
				get_str_from_stream(is,at2s); if( is.fail()) throw std::runtime_error("Error Reading Improper Dihedral Angle Params");
				get_str_from_stream(is,at3s); if( is.fail()) throw std::runtime_error("Error Reading Improper Dihedral Angle Params");
				get_str_from_stream(is,at4s); if( is.fail()) throw std::runtime_error("Error Reading Improper Dihedral Angle Params");

				boost::to_upper(at1s);
				boost::to_upper(at2s);
				boost::to_upper(at3s);
				boost::to_upper(at4s);

				if( (at1s > at4s) || ( (at1s == at4s) &&  (at2s > at3s) ) )
				{
					at1s.swap(at4s);
					at2s.swap(at3s);
				}

				double fconst;
				double ph;
				double nperiod;

				is >> fconst;
				is >> ph;
				is >> nperiod;

				if(is.fail()) throw std::runtime_error("Error Reading Improper Dihedral Angle Params");
				
				HaVec_double tpar(4);

				tpar[0] = fabs(nperiod);
				tpar[1] = ph;
				tpar[2] = fconst;
				tpar[3] = 1.0;

				std::string label = at1s + "-" + at2s + "-" + at3s + "-" + at4s;

				impdih_param_map[label] = tpar;
			}
			else if( read_mode == READ_WATER_PARAMS )
			{
				if( str_inp.empty() )
				{
					read_mode = READ_VDW_ATOM_SYNONYMS;
					continue;
				}

			}
			else if( read_mode == READ_VDW_ATOM_SYNONYMS)
			{   
				std::string ats1;
				is >> ats1;

				if( is.fail() ) 
				{
					read_mode = READ_VDW_PARAMS;
					continue;
				}
				for(;;) 
				{
					std::string ats2;
					is >> ats2;
					if( is.fail() ) break;
					multimap<std::string,std::string, less<std::string> >::value_type  str_pair(ats1,ats2);
					vdw_synonyms.insert(str_pair);
				}
			}
			else if( read_mode == READ_VDW_PARAMS)
			{
				//        if( isspace(str[0]) )
				//			{
				//				read_mode = READ_H_BOND_PARAMS;
				//				PrintLog("Reading H Bond Parameters \n");
				//				continue;
				//			}

				if( str_inp.substr(0,4) == "MOD4") { continue; }

				std::string ats1;
				is >> ats1;

				if( is.fail() ) break;

				double rad;
				double ene;
				is >> rad;
				is >> ene;
				if( is.fail() ) throw std::runtime_error("Error Reading VdW parameters");
				
				HaVec_double ppar(2);
				ppar[0] = rad;
				ppar[1] = ene;

//				PrintLog(" Read VDW params for %s \n", ats1.c_str());

				symb_ppar_map[ats1] = ppar;
				//PrintLog("%s %f\n",ats1.c_str(),rad);

				multimap<std::string, std::string, less<std::string> >::iterator syn_itr1, syn_itr2;
				syn_itr1 = vdw_synonyms.lower_bound(ats1);
				syn_itr2 = vdw_synonyms.upper_bound(ats1);
				if(syn_itr1 == vdw_synonyms.end())
					continue;
				for(;syn_itr1 != syn_itr2;syn_itr1++)
				{
					std::string& ats2 = (*syn_itr1).second;
					symb_ppar_map[ats2] = ppar;
				}
			}
			else if( read_mode == READ_H_BOND_PARAMS)  // not accessible currently
			{ 
				std::string ats4;
				is >> ats4;
				std::string ats5;
				is >> ats5;
				if( is.fail() )
				{
					PrintLog("Finished Reading Parameter Set \n");
					break;
				}

				double acoef;
				double bcoef;
				double n;
				is >> acoef;
				is >> bcoef;
				is >> n;

				if( is.fail() ) throw std::runtime_error("Error Reading H bond parameters");
				
				HaVec_double hpar(3);
				hpar[0] = acoef;
				hpar[1] = bcoef;
				hpar[2] = n; 
				std::string label = ats4 + "-" + ats5;
				symb_hpar_map[label] = hpar;
			}
		} // end for(;;) cycle on reading lines in the file
	} 
	catch( std::exception& e)
	{
		PrintLog("Exception in MMForceField::LoadAmberParamFile() %s \n", e.what());
		PrintLog("Error reading line: %s \n",str_inp.c_str() );
		return FALSE;
	}

	return TRUE;
}

int MMForceField::LoadTinkerParamFile(const std::string& ff_param_fname )
{
	char buf[256]; 
	std::string str_inp;
	StrVec tokens;
	std::string token;

	try
	{
		std::ifstream is_f(ff_param_fname.c_str());
		if(is_f.fail()) throw std::runtime_error(" Can't open file " + ff_param_fname );

		for(;;)
		{
			std::getline(is_f,str_inp);
			
			if(is_f.eof() ) return TRUE;
			if(is_f.fail()) throw std::runtime_error("Error Reading line from file" + ff_param_fname );

			boost::trim(str_inp);

			std::istringstream is(str_inp);
			// boost::split(tokens,str_inp,boost::is_any_of(" \""),boost::token_compress_on);
			is >> token;

			if( is.eof()  ) continue;
			if( is.fail() ) continue;

			if( token == "atom" )
			{


			}
			else if( token == "vdw" )
			{
				 
			}


		}

	}
	catch( const std::exception& e)
	{
		PrintLog("Exception in MMForceField::LoadTinkerParamFile() %s \n", e.what());
		PrintLog("Error reading line: %s \n",str_inp.c_str() );
		return FALSE;
	}

	return TRUE;
}

StrVec MMForceField::GetAmberParamFiles() const
{
	return amber_param_files;
}

void MMForceField::SetDefaultParamFiles()
{
	std::string res_dir = pApp->res_db_dir;
	std::string amber_parm_dir = res_dir + "amber_parm" + boost::filesystem::path::preferred_separator;
	
	amber_param_files.clear();
	resff_files.clear();

	if( ff_type == ForceFieldType::AMBER_94 )
	{
//		file_name += "amber_94_ff.dat";
		amber_param_files.push_back( amber_parm_dir + "parm94.dat"  );
		amber_param_files.push_back( amber_parm_dir + "frcmod.add1" );
		amber_param_files.push_back( amber_parm_dir + "frcmod.add2" );

		resff_files.push_back( res_dir + "resff_amber94_amino.xml" );
		resff_files.push_back( res_dir + "resff_amber_solvents.xml");
	}
	else if( ff_type == ForceFieldType::AMBER_99_SB )
	{
		amber_param_files.push_back( amber_parm_dir + "parm99.dat"    );
		amber_param_files.push_back( amber_parm_dir + "frcmod.ff99SB" );
		
		resff_files.push_back( res_dir + "resff_amber94_amino.xml"  );
		resff_files.push_back( res_dir + "resff_amber_solvents.xml" );
	}
	else if( ff_type == ForceFieldType::AMBER_99_BSC0 )
	{
		amber_param_files.push_back( amber_parm_dir + "parm99.dat"     );
		amber_param_files.push_back( amber_parm_dir + "frcmod.ff99SB"  );
		amber_param_files.push_back( amber_parm_dir + "frcmod.parmbsc0");

		resff_files.push_back( res_dir +  "resff_amber94_amino.xml");
		resff_files.push_back( res_dir +  "resff_amber_solvents.xml");
	}
	else if( ff_type == ForceFieldType::AMBER_03 )
	{
		amber_param_files.push_back(amber_parm_dir + "parm99.dat");
		amber_param_files.push_back(amber_parm_dir + "frcmod.ff03");

//		resff_files.push_back( res_dir +  "resff_amber03_amino.xml");
	} 
	else if( ff_type == ForceFieldType::AMBER_10 )
	{
		amber_param_files.push_back(amber_parm_dir + "parm10.dat");
		amber_param_files.push_back(amber_parm_dir + "frcmod.add1");
		amber_param_files.push_back(amber_parm_dir + "frcmod.add2");
		amber_param_files.push_back(amber_parm_dir + "frcmod.ionsjc_tip3p");
		amber_param_files.push_back(amber_parm_dir + "gaff.dat");
		
		resff_files.push_back( res_dir +  "resff_amber10_amino.xml");
		resff_files.push_back( res_dir +  "resff_amber10_nucleo.xml");
		resff_files.push_back( res_dir +  "resff_amber_solvents.xml");
	}
	else if( ff_type == ForceFieldType::AMOEBA )
	{	
		resff_files.push_back( res_dir +  "resff_amoeba_amino.xml");
		resff_files.push_back( res_dir +  "resff_amoeba_solvents.xml");

		tinker_param_files.push_back(amber_parm_dir + "amoebapro.prm");
	}
	else if( ff_type == ForceFieldType::ARROW_5_14_CT )
	{	
		resff_files.push_back( res_dir +  "resff_arrow_5.14_ct_amino.xml");
		resff_files.push_back(res_dir + "resff_arrow_solvents.xml");
	}
	else if (ff_type == ForceFieldType::ARROW_2_0)
	{
		resff_files.push_back(res_dir + "resff_arrow_2.0_amino.xml");
		resff_files.push_back(res_dir + "resff_arrow_solvents.xml");
	}

}

std::string MMForceField::GetAmberResName( const std::string& res_name_full, const ForceFieldType& ff_type )
{
	std::string res_name;
	std::string name_mod;
	std::string amber_res_name;

	size_t n = res_name_full.find("#");
	if( n == std::string::npos ) 
	{
		res_name = res_name_full;
		name_mod = "";
	}
	else if( n == (res_name_full.size() - 1) )
	{
		res_name = res_name_full.substr(0,(res_name_full.size() - 1) );
		name_mod = "";
	}
	else
	{
		res_name = res_name_full.substr(0,n);
		name_mod = res_name_full.substr(n+1);
	}

	amber_res_name = res_name;
	if( res_name == "HOH") amber_res_name = "WAT";
	if( res_name == "HIS") 
	{
		amber_res_name = "HID";
		if( name_mod == "EPSILON")          return "HIE";
		else if( name_mod == "EPSILON_NT" ) return "NHIE";
		else if( name_mod == "EPSILON_CT" ) return "CHIE";
		else if( name_mod == "PROT" )       return "HIP";
		else if( name_mod == "PROT_NT" )    return "NHIP";
		else if( name_mod == "PROT_CT" )    return "CHIP";
		else if( name_mod == "PROT_TO_HID") return "HIK";
	}
	else if( res_name == "CYS") 
	{
		if( name_mod == "UNPROT") return "CYX";
		else if( name_mod == "UNPROT_NT") return "NCYX";
		else if( name_mod == "UNPROT_CT") return "CCYX";
	}
	else if( res_name == "ASP") 
	{
		if( name_mod == "PROT") return "ASH";
	}
	else if( res_name == "GLU") 
	{
		if( name_mod == "PROT") return "GLH";
	}
	else if( res_name == "NA") 
	{
		if( name_mod == "NEUTRAL") return "NAN";
	}
	else if( res_name == "KA") 
	{
		if( name_mod == "NEUTRAL") return "KAN";
	} 
	else if( res_name == "MG") 
	{
		if( name_mod == "NEUTRAL") return "MGN";
	}
	else if( res_name == "CA") 
	{
		if( name_mod == "NEUTRAL") return "CAN";
	}

	if( !name_mod.empty() )
	{
		if( name_mod == "NT") 
		{
			std::string ns("N");
			amber_res_name = ns + amber_res_name;
		}
		else if( name_mod == "CT") 
		{
			std::string cs("C");
			amber_res_name = cs + amber_res_name;
		}
		else if( name_mod == "D") 
		{
			std::string ds("D");
			amber_res_name = ds + amber_res_name;
		}
		else if( name_mod == "D3") 
		{
			std::string ds("D");
			amber_res_name = ds + amber_res_name + "3";
		}
		else if( name_mod == "D5") 
		{
			std::string ds("D");
			amber_res_name = ds + amber_res_name + "5";
		}
		else if( name_mod == "DN") 
		{
			std::string ds("D");
			amber_res_name = ds + amber_res_name + "N";
		}
		else if( name_mod == "R3") 
		{
			amber_res_name = amber_res_name + "3";
		}
		else if( name_mod == "R5") 
		{
			amber_res_name = amber_res_name + "5";
		}
		else if( name_mod == "RN") 
		{
			amber_res_name = amber_res_name + "N";
		}
		else
		{
			amber_res_name += "_";
			amber_res_name += name_mod;
		}
	}
	
	return amber_res_name;
}

std::string MMForceField::GetAmberAtName( const std::string& at_name_par, const std::string& res_name_full, const ForceFieldType& ff_type )
{
	std::string atname = boost::to_upper_copy( at_name_par );
	boost::trim( atname );

	std::string at_name_amber = atname;
	
	std::string res_name;
	std::string name_mod;

	size_t n = res_name_full.find("#");
	if( n == std::string::npos ) 
	{
		res_name = res_name_full;
		name_mod = "";
	}
	else if( n == (res_name_full.size() - 1) )
	{
		res_name = res_name_full.substr(0,(res_name_full.size() - 1) );
		name_mod = "";
	}
	else
	{
		res_name = res_name_full.substr(0,n);
		name_mod = res_name_full.substr(n+1);
	}

	if( res_name == "A" || res_name == "G" || res_name == "C" 
		|| res_name == "T" || res_name == "U" )
	{
		if( atname == "O1P") at_name_amber = "OP1"; 
		else if( atname == "O2P" ) at_name_amber = "OP2";
		else if( atname == "O5X" ) at_name_amber = "O5'";
		else if( atname == "C5X" ) at_name_amber = "C5'";
		else if( atname == "H5X1" ) at_name_amber = "H5'";
		else if( atname == "H5X2" ) at_name_amber = "H5''";
		else if( atname == "C4X" ) at_name_amber = "C4'";
		else if( atname == "H4X" ) at_name_amber = "H4'";
		else if( atname == "O4X" ) at_name_amber = "O4'";
		else if( atname == "C1X" ) at_name_amber = "C1'";
		else if( atname == "H1X" ) at_name_amber = "H1'";
		else if( atname == "C3X" ) at_name_amber = "C3'";
		else if( atname == "H3X" ) at_name_amber = "H3'";
		else if( atname == "C2X" ) at_name_amber = "C2'";
		else if( atname == "H2X1" ) at_name_amber = "H2'";
		else if( atname == "H2X2" ) at_name_amber = "H2''";
		else if( atname == "O2X" ) at_name_amber = "O2'";
		else if( atname == "HOX2" ) at_name_amber = "HO2'";
		else if( atname == "O3X" ) at_name_amber = "O3'";
		else if( atname == "H5T" ) at_name_amber = "HO5'";
		else if( atname == "H3T" ) at_name_amber = "HO3'";
	}

//	if( res_name_full == "MET#NT")
//	{
//		if(at_name == "HG2") at_name_mort = "HG1";
//		if(at_name == "HG3") at_name_mort = "HG2";
//		if(at_name == "HB2") at_name_mort = "HB1";
//		if(at_name == "HB3") at_name_mort = "HB2";
//	}
	
	return at_name_amber;
}

std::string MMForceField::GetAtNameFromAmber( const std::string& atname_amber_par, const std::string& res_name_full )
{
	std::string atname_amber = boost::to_upper_copy( atname_amber_par );
	boost::trim( atname_amber );

	std::string atname = atname_amber;
	
	std::string res_name;
	std::string name_mod;

	size_t n = res_name_full.find("#");
	if( n == std::string::npos ) 
	{
		res_name = res_name_full;
		name_mod = "";
	}
	else if( n == (res_name_full.size() - 1) )
	{
		res_name = res_name_full.substr(0,(res_name_full.size() - 1) );
		name_mod = "";
	}
	else
	{
		res_name = res_name_full.substr(0,n);
		name_mod = res_name_full.substr(n+1);
	}

	if( res_name == "A" || res_name == "G" || res_name == "C" || res_name == "T" || res_name == "U" )
	{
		if( atname == "OP1") atname = "O1P";
		else if( atname == "OP2") atname = "O2P";
		else if( atname == "O5'") atname = "O5X";
		else if( atname == "C5'") atname = "C5X";
		else if( atname == "H5'")  atname = "H5X1";
		else if( atname == "H5''") atname = "H5X2";
		else if( atname == "C4'") atname = "C4X";
		else if( atname == "H4'") atname = "H4X";
		else if( atname == "O4'") atname = "O4X";
		else if( atname == "C1'") atname = "C1X";
		else if( atname == "H1'") atname = "H1X";
		else if( atname == "C3'") atname = "C3X";
		else if( atname == "H3'") atname = "H3X";
		else if( atname == "C2'") atname = "C2X";
		else if( atname == "H2'") atname = "H2X1";
		else if( atname == "H2''") atname = "H2X2";
		else if( atname == "O2'") atname = "O2X";
		else if( atname == "HO2'") atname = "HOX2";
		else if( atname == "O3'") atname = "O3X";
		else if( atname == "HO5'") atname = "H5T";
		else if( atname == "HO3'") atname = "H3T";
	}

	return atname;
}

int MMForceField::LoadResFFTemplateXMLFile(const char* fname)
{
	char buf[256];
//	PrintLog("MMForceField::LoadResFFTemplateXMLFile() pt 1 \n");
	try
	{
		TiXmlDocument doc;
		bool bres = doc.LoadFile(fname);

		if(!bres) throw std::runtime_error("Invalid XML file");

		const TiXmlElement* root_element;
		const TiXmlElement* data_element;

		root_element = doc.FirstChildElement();
		if( root_element == NULL ) throw std::runtime_error("No ROOT Element \n");
		std::string root_elem_type = root_element->Value();
		if( root_elem_type != "HARLEM_DATA" ) throw std::runtime_error("ROOT Element type is not HARLEM_DATA \n");

		HaResDB* p_res_db = HaResDB::GetDefaultResDB();

		data_element = root_element->FirstChildElement();
		while( data_element )
		{
			std::string data_elem_type = data_element->Value();
			if( data_elem_type != "resff" ) throw std::runtime_error("Data Element type is not resff \n");
			if( !data_element->CStrAttribute("rname") ) throw std::runtime_error("Empty Residue name of the template \n"); 
			std::string res_name = data_element->CStrAttribute("rname");
//			PrintLog(" MMForceField::LoadResFFTemplateXMLFile() pt 4  res_name = %s\n",res_name.c_str()); 
			HaResidue* p_res_templ = p_res_db->GetTemplateForResidue(res_name.c_str());
			if( p_res_templ == NULL ) throw std::runtime_error("No residue template with name " + res_name );
			
			ResFFTemplate* p_res_ff_templ = new ResFFTemplate(p_res_templ);

//			if( res_name_ff_templ_map.count(res_name) > 0 )
   //  		{
			//	p_res_ff_templ =  res_name_ff_templ_map[res_name];
			//	PrintLog(" Update residue force field template %s  from file %s \n", res_name.c_str(), fname );
			//}
			//else
			//{
			//	p_res_ff_templ = new ResFFTemplate(p_res_templ);
			//}

			if( !p_res_ff_templ->LoadXml( data_element ) )
			{
				delete p_res_ff_templ;
				throw std::runtime_error( "Error loading residue template " + res_name );
			}
			else
			{
				if( res_name_ff_templ_map.count(res_name) > 0 )
				{
					delete res_name_ff_templ_map[res_name];
					PrintLog(" Updating residue force field template %s from file %s \n", res_name.c_str(), fname );
				}
				res_name_ff_templ_map[res_name] = p_res_ff_templ;
			}
			data_element = data_element->NextSiblingElement();
		}
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in MMForceField::LoadResFFTemplateXMLFile() \n to process file %s \n", fname);
		PrintLog(" Error: %s \n",ex.what());
		return FALSE;
	}
//	PrintLog(" MMForceField::LoadResFFTemplateXMLFile() pt end \n");
	return TRUE;
}

ResFFTemplate* MMForceField::GetResidueTemplate(const std::string& full_res_name)
{
	std::map<std::string,ResFFTemplate*>::iterator mitr;
	mitr = res_name_ff_templ_map.find(full_res_name);
	if( mitr == res_name_ff_templ_map.end() ) return NULL;
	return (*mitr).second;
}