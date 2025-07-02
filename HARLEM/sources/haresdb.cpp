/*! \file haresdb.cpp

    Classes to manipulated Residues Database

    \author Igor Kurnikov 

    \date 2000-2002

*/

#define HARESDB_CPP

#include "tinyxml.h"

#include <filesystem>
#include <boost/algorithm/string.hpp>

#include "command.h"
#include "harlemapp.h"
#include "haresdb.h"
#include "hamolecule.h"

HaResDB* HaResDB::res_db = NULL;

HaResDB::HaResDB()
{
    SetName("RESDB");
	SetStdParams();
	Init();
}

HaResDB::~HaResDB()
{

}

void HaResDB::SetStdParams()
{
	
}

struct sort_fname_func
{
	bool operator() (std::string fn1,std::string fn2) 
	{	 
		std::string ext1 = harlem::GetExtFromFileName(fn1);
		std::string ext2 = harlem::GetExtFromFileName(fn2);

		if( ext1 != ext2 ) return (ext1 < ext2);
		return ( fn1 < fn2);
	}
} sort_fn;

int HaResDB::Init()
{
	namespace fs = std::filesystem;   

	HarlemApp* pApp_loc = GetHarlemApp();

	fs::directory_iterator ditr_main(pApp_loc->res_db_dir);

    PrintLog(" Initialize Residue Database \n");
	for (; ditr_main != fs::directory_iterator(); ditr_main++)
	{
		std::string file_name = ditr_main->path().filename().string();
		//		PrintLog(" file_name = %s \n",file_name.c_str() );
		if (boost::starts_with(file_name, "res_"))
		{
			//std::string full_name = pApp->res_db_dir + file_name;
			//res_db_files.push_back(full_name);
			res_db_files.push_back(file_name);
			if( ha_debug_level > 1) PrintLog(" Add residue template file %s  in the Residue Database directory \n", file_name);
		}
	}

	std::sort(res_db_files.begin(),res_db_files.end(),sort_fn);

	fs::path cur_path = fs::current_path();

	try 
	{
		fs::directory_iterator ditr(cur_path);

		for (; ditr != fs::directory_iterator(); ditr++)
		{
			std::string file_name = ditr->path().filename().string();
			//		PrintLog(" file_name = %s \n",file_name );
			if (boost::starts_with(file_name, "res_"))
			{
				std::string full_name = "./" + file_name;
				res_db_files.push_back(full_name);
				PrintLog(" Add residue template file %s  in the working directory \n", file_name);
			}
		}
	}
	catch (fs::filesystem_error& ex)
	{
		PrintLog("%s \n", ex.what());
	}

	for(std::string db_file_name : res_db_files)
	{
		fs::path db_file_path(pApp->res_db_dir);
		if( !boost::starts_with(db_file_name, ".") ) 
		{
			db_file_path /= db_file_name;
			db_file_name = db_file_path.string();
		}
		std::string ext_str = harlem::GetExtFromFileName( db_file_name );
		boost::to_lower(ext_str);
		
		if(!fs::exists(db_file_name) )
		{
			PrintLog("HaResDB::Init(): Can't find file %s \n", db_file_name);
			continue;
		}
		if( ext_str == "hlm" )
		{
			LoadHarlemFile(db_file_name);
		}
		else if (ext_str == "mol2")
		{
			LoadMol2File(db_file_name);
		}
	}

	for (HaMolecule* pmol : HostMolecules)  /// set map of residue templates 
	{
		std::string mol_name = pmol->GetName();

		if (HaResidue::ResSynonym_to_std.count(mol_name) > 0) mol_name = HaResidue::ResSynonym_to_std[mol_name];

		if (res_name_templ_map.count(mol_name) > 0)
		{
			PrintLog("HaResDB::Init():  Residue Template %s Is already in DB - overiding \n", mol_name);
		}
		res_name_templ_map[mol_name] = pmol;
	}

	// Modify templates with Instructions in XML files - Probably is too complicated - remove

	for (std::string db_file_name : res_db_files)
	{
		fs::path db_file_path(pApp->res_db_dir);
		if (!boost::starts_with(db_file_name, "."))
		{
			db_file_path /= db_file_name;
			db_file_name = db_file_path.string();
		}
		std::string ext_str = harlem::GetExtFromFileName(db_file_name);
		boost::to_lower(ext_str);

		if (ext_str != "xml") continue;

		if (!fs::exists(db_file_name))
		{
			PrintLog("HaResDB::Init(): Can't find file %s \n", db_file_name);
			continue;
		}
		LoadXMLFile(db_file_name);
	}

	return TRUE;
}

HaResDB* HaResDB::GetDefaultResDB()
{
	MolSet* pmset_save = GetCurMolSet();
	if( res_db == nullptr) res_db = new HaResDB;
	if( pmset_save) CurMolSet = pmset_save;
	return res_db;
}

HaMolecule* HaResDB::GetMolTemplForRes(const std::string& res_templ_name)
{
	std::string std_res_name = res_templ_name;
	if (HaResidue::ResSynonym_to_std.count(std_res_name) > 0) std_res_name = HaResidue::ResSynonym_to_std[std_res_name];
	
	if (res_name_templ_map.count(std_res_name) == 0) return nullptr;

	HaMolecule* p_templ = res_name_templ_map[std_res_name];
	return p_templ;
}

HaResidue* HaResDB::GetTemplateForResidue(const std::string& res_full_name)
{
	std::string std_res_name = res_full_name;
	if (HaResidue::ResSynonym_to_std.count(std_res_name) > 0) std_res_name = HaResidue::ResSynonym_to_std[std_res_name];

	std::string res_sh_name = HaResidue::GetResNameFromFullName(std_res_name);

	HaMolecule* templ_mol = GetMolTemplForRes(std_res_name);
	if( templ_mol == NULL)
	{
		PrintLog("Error in HaResDB::GetTemplateForResidue() \n Can't find template for residue %s in DB \n", res_full_name);
		return nullptr;
	}

	HaResidue* prtempl;
	ChainIteratorMolecule chitr_templ(templ_mol);
	HaChain* chain = chitr_templ.GetFirstChain();

// find in the residue template residue with the name of the current residue 
	ResidueIteratorChain ritr_ch(chain);
	for(prtempl = ritr_ch.GetFirstRes(); prtempl; prtempl= ritr_ch.GetNextRes())
	{
		if(!strcmp(prtempl->GetName(), res_sh_name.c_str() ) )
			break;
	}
	
	if (prtempl == NULL)
	{
		PrintLog(" HaResDB::GetTemplateForResidue \n Can't find residue named %s in the residue template %s",
			       res_sh_name, std_res_name);
		return NULL;
	}
	
	return prtempl;
}


HaAtom* HaResDB::GetTemplateForAtom(HaAtom* aptr)
{
	HaResidue* pres = aptr->GetHostRes();
	std::string full_res_name = pres->GetFullName();
	HaResidue* prtempl = GetTemplateForResidue(full_res_name);
	if(prtempl == NULL) return NULL;

	HaAtom* atempl = prtempl->GetAtomByName(aptr->GetName());

	HaMolecule* mol_templ = prtempl->GetHostMol();
	if( atempl == NULL )
	{
		PrintLog(" Can't find Atom named %s in the residue template %s \n",
			    aptr->GetName(), mol_templ->GetObjName() );
		return NULL;
	}

	return atempl;
}

int HaResDB::GetTemplResAtomMaps( HaResidue* pres, AtomAtomMap& res_to_templ_map, AtomAtomMap& templ_to_res_map)
{
	res_to_templ_map.clear();
	templ_to_res_map.clear();

	std::string res_fname = pres->GetFullName();
	HaMolecule* res_templ = GetMolTemplForRes(res_fname.c_str());
	
	if( res_templ == NULL ) return FALSE;

	HaResidue* prtempl;
	ChainIteratorMolecule chitr_templ(res_templ);
	HaChain* chain = chitr_templ.GetFirstChain();

// find residue in the residue template with the name of the current residue 
	ResidueIteratorChain ritr_ch(chain);
	for(prtempl = ritr_ch.GetFirstRes(); prtempl; prtempl= ritr_ch.GetNextRes())
	{
		if(!strcmp(prtempl->GetName(), pres->GetName()) )
			break;
	}

	if (prtempl == NULL)
	{
		PrintLog(" Error in HaResDB::GetTemplResAtomMaps() \n Can't find residue named %s in the residue template %s",
			    pres->GetName(), res_templ->GetObjName() );
		return FALSE;
	}
    
	HaAtom* atempl;
	HaAtom* aptr;

	AtomIteratorAtomGroup aitr_prtempl(prtempl);
	for(atempl = aitr_prtempl.GetFirstAtom(); atempl; atempl = aitr_prtempl.GetNextAtom())
	{
		aptr = pres->GetAtomByName(atempl->GetName());
		if( aptr == NULL )
		{
			if( atempl->IsProxy() ) continue;
			char buf[256];
			pres->FillRef(buf);
			PrintLog(" Error in HaResDB::GetTemplResAtomMaps() \n");
			PrintLog("Can't find atom named %s in the residue %s found in the template %s",
					atempl->GetName(), buf, res_templ->GetObjName() );
			templ_to_res_map[atempl] = NULL;
			continue;
		}
		templ_to_res_map[atempl] = aptr;
		res_to_templ_map[aptr] = atempl;
	}

	// map proxy atoms of the residue template:
	for(atempl = aitr_prtempl.GetFirstAtom(); atempl; atempl = aitr_prtempl.GetNextAtom())
	{
		if( !atempl->IsProxy() ) continue;
		std::string atempl_name = atempl->GetName();
//		PrintLog(" HaResDB::GetTemplResAtomMaps() pt 1  proxy atom: %s \n", atempl_name.c_str());
		if( templ_to_res_map.count( atempl ) > 0 ) continue; // atom is already mapped
		try
		{
			AtomGroup bonded_atoms;
			atempl->GetBondedAtoms( bonded_atoms );
			
			if(  bonded_atoms.size() == 0 ) throw std::runtime_error( " Proxy Atom does not have bonded atoms " );
			HaAtom* aptr_tb = bonded_atoms[0];
			std::string atb_name = aptr_tb->GetName();
//			PrintLog(" HaResDB::GetTemplResAtomMaps() pt 2  atom bonded to proxy atom: %s \n", atb_name.c_str());

			if( templ_to_res_map.count( aptr_tb ) == 0 ) throw std::runtime_error(" atom " + atb_name + " bonded to proxy atom is not mapped ");
			HaAtom* aptr_b = templ_to_res_map[aptr_tb];
			AtomGroup bonded_atoms_2;
			aptr_b->GetBondedAtoms( bonded_atoms_2 );
			AtomIteratorAtomGroup aitr_2( &bonded_atoms_2 );
			for( aptr = aitr_2.GetFirstAtom(); aptr; aptr = aitr_2.GetNextAtom() )
			{
				if( aptr->GetHostRes() == pres ) continue;
				std::string repl_at_name = atempl->GetReplacedAtName();
//				PrintLog(" HaResDB::GetTemplResAtomMaps() pt 3  expected replaced atom: %s \n", atb_name.c_str());
				std::string at_name = aptr->GetName();
//				PrintLog(" HaResDB::GetTemplResAtomMaps() pt 4  actual replaced atom: %s \n",   at_name.c_str());
				if( !repl_at_name.empty() && (repl_at_name != at_name) ) continue;
				break;
			}
			if( aptr == NULL ) 
			{
				PrintLog(" Warning: in HaResDB::GetTemplResAtomMaps() \n");
				PrintLog(" proxy_atom %s in residue template residue %s is not mapped  \n", atempl_name.c_str(), res_fname.c_str() );
				continue;
			}
		}
		catch( std::exception& ex )
		{
			PrintLog(" Error in HaResDB::GetTemplResAtomMaps() \n");
			PrintLog(" Setting atom map for residue %s  proxy_atom %s \n",res_fname.c_str(), atempl_name.c_str() );
			PrintLog(" Setting atom map for residue %s \n",res_fname.c_str() );
			PrintLog("%s\n",ex.what());
			continue;
		}
		templ_to_res_map[atempl] = aptr;
		res_to_templ_map[aptr] = atempl;
	}

	AtomIteratorAtomGroup aitr_pres(pres);
	for(aptr = aitr_pres.GetFirstAtom(); aptr; aptr = aitr_pres.GetNextAtom())
	{
		if( res_to_templ_map.find(aptr) != res_to_templ_map.end() )
		{
			atempl = prtempl->GetAtomByName(aptr->GetName() );
			if( atempl == NULL )
			{
				PrintLog(" Error in HaResDB::GetTemplResAtomMaps() \n");
				PrintLog("Can't find atom named %s in the template %s found in the residue %s",
						atempl->GetName(), res_templ->GetObjName(), pres->GetName() );
				res_to_templ_map[aptr] = NULL;
				continue;
			}
		}
	}

// find corresponding atoms to the terminating atoms of the residue template

	AtomIteratorMolecule aitr_templ(res_templ);
	for(atempl = aitr_templ.GetFirstAtom(); atempl; atempl = aitr_templ.GetNextAtom() )
	{
		if( prtempl->HasAtom(atempl) )
			continue;
		
		AtomGroup bonded_atoms;
		atempl->GetBondedAtoms(bonded_atoms);
		bool found_match = false;
        
		HaAtom* atempl2;
		AtomIteratorAtomGroup aitr_bonded_atoms(&bonded_atoms);
		for(atempl2 = aitr_bonded_atoms.GetFirstAtom(); atempl2; atempl2 = aitr_bonded_atoms.GetNextAtom() ) 
		{
			if( !prtempl->HasAtom(atempl2) ) continue;
			aptr = templ_to_res_map[atempl2];
			if( aptr == NULL ) continue;

			AtomGroup bonded_atoms_2;
			aptr->GetBondedAtoms(bonded_atoms_2);
            HaAtom* neib_atom;
			AtomIteratorAtomGroup aitr_bonded_atoms_2(&bonded_atoms_2);
			for( neib_atom = aitr_bonded_atoms_2.GetFirstAtom(); neib_atom; neib_atom = aitr_bonded_atoms_2.GetNextAtom() )
			{
				if( pres->HasAtom( neib_atom ) )
					continue;
				templ_to_res_map[atempl] = neib_atom;
				res_to_templ_map[neib_atom] = atempl;
				found_match = true;
				break;
			}
			if( found_match) break;
		}
		if(!found_match) templ_to_res_map[atempl] = NULL;
	}
	return TRUE;
}

int HaResDB::GetTemplResAtNameAtomMap( HaResidue* pres, StrAtomMap& templ_atname_to_res_map)
{
	templ_atname_to_res_map.clear();
	AtomAtomMap res_to_templ_map;
	AtomAtomMap templ_to_res_map;

	int ires = GetTemplResAtomMaps( pres, res_to_templ_map, templ_to_res_map);
	if( !ires ) return FALSE;

	AtomAtomMap::iterator mitr;
	for( mitr = templ_to_res_map.begin(); mitr != templ_to_res_map.end(); mitr++ )
	{
		HaAtom* atempl = (*mitr).first;
		HaAtom* aptr   = (*mitr).second;
		std::string atname_t = atempl->GetName();
		templ_atname_to_res_map[atname_t] = aptr;
	}
	return TRUE;
}

int HaResDB::LoadXMLFile( const std::string& fname )
{
	// PrintLog(" HaResDB::LoadXMLFile() fname = %s \n", fname);
	try
	{
		TiXmlDocument doc;
		bool bres = doc.LoadFile(fname); 

		if(!bres) throw std::runtime_error( " Invalid XML file "); 

		const TiXmlElement* root_element = doc.FirstChildElement();
		if( root_element == NULL ) throw std::runtime_error( "No ROOT Element ");
		std::string root_elem_tag = root_element->Value();
		if( root_elem_tag != "HARLEM_DATA" ) throw std::runtime_error("ROOT Element tag is not HARLEM_DATA ");

		const TiXmlElement* res_templ_element = root_element->FirstChildElement();
		while( res_templ_element )
		{
			std::string elem_tag = res_templ_element->Value();
			if( elem_tag == "res_templ" ) 
			{
				if( !res_templ_element->CStrAttribute("name") ) throw std::runtime_error("Empty name of the residue template " ); 
				std::string res_name = res_templ_element->CStrAttribute("name");
				HaResidue* p_res_templ = NULL;
				std::string action = "new";
				if( res_templ_element->CStrAttribute("action") ) action = res_templ_element->CStrAttribute("action");

				if( action == "modify" )
				{
					p_res_templ = this->GetTemplateForResidue( res_name );
					if( p_res_templ == NULL ) throw std::runtime_error( "No residue template with name " + res_name );
				}
				if( p_res_templ == NULL) throw std::runtime_error( "No residue template set to add " );
			
				const TiXmlElement* data_element = res_templ_element->FirstChildElement();
				while( data_element )
				{
					std::string data_elem_tag = data_element->Value();
					if( data_elem_tag == "atom")
					{
						std::string at_name = "";
						HaAtom* aptr = NULL;
						if( data_element->CStrAttribute("name") ) at_name = data_element->CStrAttribute("name");
						aptr = p_res_templ->GetAtomByName( at_name.c_str() );
						if( !aptr ) aptr = p_res_templ->AddNewAtom();
						aptr->LoadXml( data_element );
					}
					data_element = data_element->NextSiblingElement();
				}
			}
			res_templ_element = res_templ_element->NextSiblingElement();
		}
	}
	catch(std::exception& ex) 
	{
		PrintLog(" Error in HaResDB::LoadXMLFile() \n");
		PrintLog(" Loading Residue template file %s \n", fname.c_str() );
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return TRUE;
}
