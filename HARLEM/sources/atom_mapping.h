#ifndef ATOM_MAPPING_H
#define ATOM_MAPPING_H

#include "vec3d.h"
#include "haintcrd.h"

class HaAtom;
class AtomContainer;

class AtomMapping
//!< Class to set mapping rules of atoms of one atom collection to another
{
public:
	AtomMapping();
	AtomMapping(AtomContainer* p_ac_1, AtomContainer* p_ac_2);
	virtual ~AtomMapping();

	void PrintInfo(int detailed = FALSE); //!< Print description of the atom mappping

	AtomContainer* p_ac_1; //!< Atom Collection 1
	AtomContainer* p_ac_2; //!< Atom Collection 2

	std::map< HaAtom*, HaAtom*, less<HaAtom*> > atmap_2to1; //!< Map of atoms of Atom Collection 2 to those of Atom Collection 1 
	std::map< HaAtom*, HaAtom*, less<HaAtom*> > atmap_1to2; //!< Map of atoms of Atom Collection 1 to those of Atom Collection 2 

	std::vector<SingleAtomCrdRule*> SyncRules2from1; //!< Rules to synchronize coordinates of atom collection 2 from coordinates of the atom collection 1  
	std::vector<SingleAtomCrdRule*> SyncRules1from2; //!< Rules to synchronize coordinates of atom collection 1 from coordinates of the atom collection 2  

	void ClearSyncRules1from2(); //!< Clear Synchronization rules for atoms of Atom Collection 1 from those of Atom Collection 2
	void ClearSyncRules2from1(); //!< Clear Synchronization rules for atoms of Atom Collection 2 from those of Atom Collection 1
	
	int Map2to1ByAtomDistance(); //!< Map atoms of Atom Collection 2 to atoms of Atom Collection 1 based on atom-atom distance, bonding etc 
	int Map2to1ByAtomRef(); //!< Map atoms of Atom Collection 2 to atoms of Atom Collection 1 based on atom ids (names, residues names and numbers, molecule names)
	int AssociateAtomsRefAtoms(HaAtom* aptr_mng, HaAtom* aref_1, HaAtom* aref_2, HaAtom* aref_3, 
		 double dist = -1.0, double vang = -1.0, double torsion = -500.0, int priority = -1 ); //!< Map atom relative other 3 atoms

	int SyncAtomCrd1From2(); //!< Sync coordinates of atoms of Atom Collection 1 from those of Atom Collection 2
	int SyncAtomCrd2From1(); //!< Sync coordinates of atoms of Atom Collection 2 from those of Atom Collection 1

	static int BuildSyncRulesForMissingAtoms(AtomContainer& all_atoms, AtomContainer& known_atoms, 
		                                     std::vector<Pt3CrdRule*>& rules); //!< build rules for coordinates of atoms in Atom Collection (all_atoms) based on coordinates of atoms assumed known (atoms_known)
};

#endif //!< #ifndef ATOM_MAPPING_H