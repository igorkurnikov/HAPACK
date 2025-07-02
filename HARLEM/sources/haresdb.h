/*! \file haresdb.h

    Classes to manipulate Residue Database in HARLEM

    \author  Igor Kurnikov 
    \date 2000-2012

*/
#if !defined(HARESDB_H)
#define HARESDB_H

#include "hamolset.h"

class HaResDB : public MolSet
//! Class to define Residue Database 
{
public:
	HaResDB();
	virtual ~HaResDB();

protected:
	void SetStdParams();

public:
	int Init();          //!< Init residue database
	HaMolecule* GetMolTemplForRes(const std::string& res_templ_name);     //!< Get a pointer for a molecule with the residue template 
	HaResidue*  GetTemplateForResidue(const std::string& res_full_name);  //!< Get a template residue by full name
	HaAtom* GetTemplateForAtom(HaAtom* aptr);  //!< Get A template for an atom

//! 
	int GetTemplResAtomMaps( HaResidue* pres, AtomAtomMap& res_to_templ_map, AtomAtomMap& templ_to_res_map); //!< Get maps between Atom pointers of the residue and its template
	int GetTemplResAtNameAtomMap( HaResidue* pres, StrAtomMap& templ_atname_to_res_map);  //!< Get a map between residue template atom names and residue atom pointers 
	int LoadXMLFile( const std::string& fname ); //!< Load File with residue templates 

	static HaResDB* GetDefaultResDB(); //!< Get Default Residue DB 

public:
	std::vector<std::string> res_db_files; 
	std::map<std::string, HaMolecule*> res_name_templ_map; // map of full residue names to residue templates

protected:
	static HaResDB* res_db;
};


#endif // !defined(HARESDB_H)
