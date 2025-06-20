/*!  \file mm_force_field.h

    Classes to describe Molecular Mechanics force field in HARLEM

    \author Igor Kurnikov 
    \date 2010-

*/
#ifndef MM_FORCE_FIELD_H
#define MM_FORCE_FIELD_H

#include "mm_params.h"

class ResFFTemplate;
class MortForceField;
namespace mort
{
	class molecule_t;
}

class MMForceField
{
public:
	MMForceField();
	virtual ~MMForceField();

	int Init(); //!< Initialize force field

	int SetFFType( const ForceFieldType& ff_type);  //!< Set Force Field Type
	ForceFieldType GetFFType() const;                //!< Get Force Field Type

	void Clear(); //< Clear Force Field Data
	void SetDefaultVdW(); //!< Set Default VdW Parameters 
	void SetDefaultParamFiles(); //!< Set default FF Parameter file names  
	
	StrDoubleMap symb_mass_map;   //!< Map ff atom symbols to atomic mass
	map<string,vector<double>>  symb_ppar_map;     //!< Map strings to VdW atom parameters
	map<string, vector<double>>  symb_hpar_map;     //!< Map strings to H bond atom parameters
	map<string, vector<double>>  bond_param_map;    //!< Map strings to Valence Bonds parameters
	map<string, vector<double>>  vang_param_map;    //!< Map strings to Valence Angles parameters  
	map<string, vector<double>>  dih_param_map;     //!< Map strings to Dihedral Angles parameters
	map<string, vector<double>>  impdih_param_map;  //!< Map strings to Improper Dihedral Angles Parameters

	double FindAtomMassFromSymbol( const std::string& ats );    //!< Find atom mass for atom ff symbol ( return -1.0 if not defined )
	vector<double>  FindPointParamFromSymbol( const string& ats1);  //!< Find Force Field atom parameters for FFSymbol
	vector<double>  FindBondParamFromSymbol ( const string& ats1, const string& ats2);
	vector<double>  FindHBondParamFromSymbol (const string& ats1, const string& ats2);
	vector<double>  FindValAngleParamFromSymbol(const string& ats1, const string& ats2, const string& ats3);
	vector<double>  FindDihedralParamFromSymbol(const string& ats1, const string& ats2, const string& ats3, const string& ats4, bool improper_flag = false);
 
	int LoadAmberParamFile(const string& ff_param_fname );  //!< Load Force-Field parameters from file in AMBER format
	int LoadTinkerParamFile(const string& ff_param_fname );  //!< Load Force-Field parameters from file in TINKER format
	StrVec GetAmberParamFiles() const; //!< Get file names of Force Field parameter files in AMBER format

	static std::vector< MMForceField* > ff_arr;  //!< Array of initialized Force fields
	static MMForceField* GetMMForceField(const ForceFieldType& ff_type, int create = FALSE ); //!< Get and initialize MM Force Field by name, create empty if not foud and create = TRUE  
	static ForceFieldType ff_type_default;    //!< Default Force Field Type
	static std::string GetAmberResName( const string& full_res_name, const ForceFieldType& ff_type = ForceFieldType::UNKNOWN_FF ); //!< Get Mort Library residue name for particular full residue name ( with modifier )
	static std::string GetAmberAtName(  const string& at_name, const string& full_res_name, const ForceFieldType& ff_type = ForceFieldType::UNKNOWN_FF ); //!< Get AMBER Library atom name for particular full residue name ( with modifier )
	static std::string GetAtNameFromAmber( const std::string& at_name_amber, const std::string& full_res_name ); //!< Get AMBER library atom name from HARLEM atom name 

	int IsMortFFInitiated();   //!< Check if MORT Data structures of the force field are initiated 
	int InitMortFF();          //!< Init MORT database of residues and Force Field parameters  
	MortForceField* p_mort_ff; //!< MORT library objects that define force field 

	int SaveResFFTemplatesFromMort(const char* fname, MolSet* pmset_res); //!< Save force field of residues of pmset_res (set using Mort library) into XML file 
	int LoadResFFTemplateXMLFile(const char* fname); //!< Load XML file of Residue Force Field Templates 

	static void switch_amoeba_type( mort::molecule_t& mol, const mort::molecule_t& poleff );

	ResFFTemplate* GetResidueTemplate(const string& full_res_name); //!< Get Residue Force Field Template by a full Residue name

	static StrVec resff_files_add;        //!< Additional Residue force field templates files
	static StrVec tinker_param_files_add; //!< Additional force field parameters files in TINKER format
	static StrVec amber_param_files_add;  //!< Additional force field parameters files in AMBER format

private:
	ForceFieldType ff_type;  //!< Force Field Type (AMBER94, CHARMM22 etc..) 

	StrVec amber_param_files;   //!< Names of force field parameter files in AMBER format
	StrVec tinker_param_files;  //!< Names of force field parameter files in TINKER format
	StrVec resff_files;  //!< Names of files of Force field residue templates

	std::vector<ResFFTemplate*>           res_ff_templates;        //!< Residue Force Field Templates
	std::map<std::string,ResFFTemplate*>  res_name_ff_templ_map;   //!< Map of Residue Names to Residue Force Field Templates
};

#endif // end if !defined(MM_FORCE_FIELD_H) 