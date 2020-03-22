/*!  \file mm_driver_tinker.cpp

    Molecular Mechanics simulations using TINKER package  

    \author Igor Kurnikov 
    \date 2010-

*/
#include <mpi.h>

#include <tinyxml.h>

#include <boost/algorithm/string.hpp>

#include "harlemapp.h"
#include "hamolmech.h"
#include "mm_driver_tinker.h"

MMDriverTinker::MMDriverTinker(HaMolMechMod* p_mm_mod_new)
{
	p_mm_mod = p_mm_mod_new;
	pmset = p_mm_mod->GetMolSet();

	to_save_input_files = TRUE;
}

MMDriverTinker::~MMDriverTinker()
{

}

int MMDriverTinker::SaveAllInpFiles()
{
	
	to_save_input_files = FALSE;
	return TRUE;
}

int MMDriverTinker::InitFFNumMap()
{
    DIR *dir = opendir(pApp->res_db_dir.c_str());    
    if(dir)
    {
        struct dirent *entry;
        while((entry = readdir(dir)) != 0)
        {
//            PrintLog( " res file %s \n",entry->d_name);
			std::string prm_file_name = entry->d_name;
			if( boost::starts_with(prm_file_name,"prm_") && boost::ends_with(prm_file_name,".xml") )
			{
				LoadFFNumMapFromFile(prm_file_name.c_str());
			}
        }
        closedir(dir);
		to_init_ff_map = FALSE;
    }
    else
    {
        ErrorInMod("HaResDB::Init()"," Can't find residue template directory "); 
		return FALSE;
    }
	return TRUE;
}

int MMDriverTinker::LoadFFNumMapFromFile( const char* fname)
{
	TiXmlDocument doc;
	bool bres = doc.LoadFile(fname);

	if(!bres)
	{
		PrintLog("Error in MMDriverTinker::LoadFFNumMapFromFile() \n");
		PrintLog("Loading file %s \n",fname);
		return FALSE;
	}

	TiXmlElement* root_elem = doc.FirstChildElement();
	if( root_elem == NULL )
	{
		PrintLog("Error in MMDriverTinker::LoadFFNumMapFromFile() \n");
		PrintLog("No ROOT Element in XML file %s \n",fname);
		return FALSE;
	}	

	std::string root_name = root_elem->Value();
	boost::to_upper(root_name);

	if( root_name != "FF_ATOM_SYMBOL_MAP")
	{
		PrintLog("Error in MMDriverTinker::LoadFFNumMapFromFile() \n");
		PrintLog("ROOT Element in XML file %s is not FF_ATOM_SYMBOL_MAP \n",fname);
		return FALSE;
	}
	std::string ff_name = root_elem->Attribute("ff_name");
	if( ff_name.empty() )
	{
		PrintLog("Error in MMDriverTinker::LoadFFNumMapFromFile() \n");
		PrintLog("No FF_NAME in ROOT Element in XML file %s \n",fname);
		return FALSE;
	}
	boost::to_upper(ff_name);

	TiXmlElement* p_data_elem;
	for(p_data_elem = root_elem->FirstChildElement("DATA"); p_data_elem; p_data_elem = p_data_elem->NextSiblingElement("DATA"))
	{
		std::string ff_smbl_type = p_data_elem->Attribute("ff_smbl");
		if(ff_smbl_type != "tinker_at_type") continue;
		std::vector<std::string> res_templ_names;
		TiXmlAttribute* p_attr;
		for( p_attr = p_data_elem->FirstAttribute(); p_attr; p_attr = p_attr->Next())
		{
			std::string name  = p_attr->Name();
			std::string value = p_attr->Value();
			boost::to_lower(name);
			boost::to_upper(value);
			if(name != "res") continue;
			res_templ_names.push_back(value);
		}
		int nres = res_templ_names.size();
		int ires;
		std::string str_map = p_data_elem->GetText();
		std::istringstream is(str_map);
		while( is.good() )
		{
			std::string at_name;
			int at_ff_num;
			is >> at_name;
			is >> at_ff_num;
			if ( is.bad() ) break;
			
			boost::to_upper(at_name);

			for(ires = 0; ires < nres; ires++)
			{
				std::string at_str_full = ff_name + ":";
				at_str_full += res_templ_names[ires];
				at_str_full += ".";
				at_str_full += at_name;
				atn_ff_num_map[at_str_full] = at_ff_num;
			}
		}
	}
	return TRUE;
}

int MMDriverTinker::GetAtFFNum(const char* force_field_name, const char* res_templ_name, const char* at_name)
{
	std::string ff_name_loc = force_field_name;
	std::string res_templ_name_loc = res_templ_name;
	std::string at_name_loc = at_name;

	boost::trim(ff_name_loc);
	boost::trim(res_templ_name_loc);
	boost::trim(at_name_loc);

	std::string at_str_full = ff_name_loc + ":";
	at_str_full += res_templ_name_loc + ".";
	at_str_full += at_name_loc;

	boost::to_upper(at_str_full);
	
	StrIntMap::iterator itr = atn_ff_num_map.find(at_str_full);
	if( itr == atn_ff_num_map.end()) return 0;

	int ff_num = (*itr).second;
	return ff_num;
}


//	std::auto_ptr<XercesDOMParser> p_parser(new XercesDOMParser);
//	try 
//	{  
//		p_parser->parse(fname);
//		DOMDocument* xmlDoc = m_ConfigFileParser->getDocument(); 
//		DOMElement* p_elem_root = xmlDoc->getDocumentElement();  
//		if( !p_elem_root ) 
//		{
//			PrintLog(" Error in MMDriverTinker::LoadFFNumMapFromFile() \n");
//			PrintLog(" Empty root element \n");
//			return FALSE;
//		}
//	}
//	catch( xercesc::XMLException& e ) 
//	{		
//		char* msg = XMLString::transcode( e.getMessage() ); 
//		PrintLog(" Error in MMDriverTinker::LoadFFNumMapFromFile() \n");
//		PrintLog(" Parsing file %s \n",fname);
//		PrintLog(" Error message: %s \n,msg);
//		XMLString::release( &msg );
//		return FALSE;
//	}  



//int MMDriverTinker::LoadFFNumMapFromNode( DOMElement* p_node )
//{
// Commented Xerces-c XML processing for now
//
//	if( p_node == NULL ) return FALSE;
//	try
//	{
//		DOMNodeList* p_child_list = p_node->getChildNodes();  
//		const XMLSize_t nn = p_child_list->getLength();
//		int i;
//		for(i = 0; i < nn; i++)
//		{
//			DOMNode* p_node = p_child_list->item(i);
//			if( p_node->getNodeType() &&  // true is not NULL  
//				p_node->getNodeType() == DOMNode::ELEMENT_NODE ) 
//			{
//				p_elem = dynamic_cast< xercesc::DOMElement* >( p_node );
//				if( p_elem->getTagName()
//
//			}
//
//		}
//	}
//	catch ( xercesc::XMLException& e ) 
//	{
//		char* msg = XMLString::transcode( e.getMessage() ); 
//		PrintLog(" Error in MMDriverTinker::LoadFFNumMapFromXMLElem() \n");
//		PrintLog(" Error message: %s \n,msg);
//		XMLString::release( &msg );
//		return FALSE;
//	}
//
//}