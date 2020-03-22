#ifndef HAFLEXMOD_H
#define HAFLEXMOD_H

#include "hastl.h"
#include "hastring.h"
#include "command.h"
#include "hacompmod.h"
#include "halinalg.h"
#include "vec3d.h"
#include "haatom.h"
#include "haatgroup.h"
#include "hahbhp.h"

#include <string>
using namespace std;

class MolEditor; 
class HaMolView;
class HBondAvg;
class HaHydrophobicTether;

class HaFlexMod: public HaCompMod 
{
public:
	HaFlexMod(MolSet* new_pmset=NULL);
	~HaFlexMod();

	void Init();

	void ClearHB();   //!< Clear Hydrogen bond array
	void ClearHPT();  //!< Cleat Hydrophobic tethers array

public:
    MolEditor* p_mol_editor;      //!< Pointer to Molecular Editor object (to find hydrogen bonds)
	double hydrophob_dist_cutoff; //!< Distance cutoff for hydrophobic Tethers (Ang)
	double hpt_duty_cycle_cutoff;  //!< Hydrophobic Tethers cutoff based on duty cycle value ( stay on percent )
	double hb_duty_cycle_cutoff;  //!< H-bond cutoff based on duty cycle value  ( stay on percent )
	double hb_energy_cutoff;      //!< H-bond energy cutoff ( kcal/mol) 

	std::string active_at_array_id; //!< ID of the active atom array  

	AtomIntMap atom_ref_map;

	vector<HBondAvg*> HBArray;
	vector<HaHydrophobicTether*> HPTArray;

	vector<HBondAvg*> selected_hb;              //!< selected hydrogen bonds
	vector<HaHydrophobicTether*> selected_hpt;  //!< selected Hydrophobic tethers

	std::string struct_file_name;
	std::string first_data_file_name;

	std::string md_traj_file_name;

	int SaveHBondFile();
	int SaveHydrophobTethersFile();
	int SaveStructFile();
	int SaveFirstInpFiles();
	int ReadFirstDataFile();
	int ComputeBondsDutyCycleAlongMD();
	int DeleteSelectedHBonds();
	int DeleteSelectedHPTethers();

	void AddAllAtomGroup();
	void createAtomRefMap();

	int Display(HaMolView* pView);

public:
	int FindHydrogenBonds();
	int FindHydrophobicContacts();
};

#endif


