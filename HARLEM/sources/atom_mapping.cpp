/*!  \file atom_mapping.cpp

   Classes to define mapping of Atom Collection to each other and synchronization of their coordinates
   
   \author Igor Kurnikov  
   \date 2010-

*/

#include <memory>

#include "haio.h"

#include "rapidxml.hpp"

#include "haatom.h"
#include "haatgroup.h"
#include "atom_mapping.h"

AtomMapping::AtomMapping()
{
	p_ac_1   = NULL;
	p_ac_2   = NULL;
}

AtomMapping::AtomMapping(AtomContainer* p_ac_1_new, AtomContainer* p_ac_2_new)
{
    p_ac_1 = p_ac_1_new;
	p_ac_2 = p_ac_2_new;
}

AtomMapping::~AtomMapping()
{
	ClearSyncRules1from2();
	ClearSyncRules2from1();
}

void AtomMapping::ClearSyncRules1from2()
{
	int nr = SyncRules1from2.size();
	int i;
	for(i = 0; i < nr; i++)
	{
		delete SyncRules1from2[i];
	}
	SyncRules1from2.clear();
}

void AtomMapping::ClearSyncRules2from1()
{
	int nr = SyncRules2from1.size();
	int i;
	for(i = 0; i < nr; i++)
	{
		delete SyncRules2from1[i];
	}
	SyncRules2from1.clear();
}

int AtomMapping::SetAtom3PtSyncRule(HaAtom* aptr_mng, HaAtom* aref_1, HaAtom* aref_2, HaAtom* aref_3,
										double dist, double vang, double torsion, int priority )
{
	int idx_ac = 0;
	if( p_ac_1->IsMember(aptr_mng) ) idx_ac = 1;
	if( p_ac_2->IsMember(aptr_mng) ) idx_ac = 2;

	if(idx_ac == 0)
	{
		PrintLog("Error in AtomMapping::AssociateAtomsRefAtoms \n");
		PrintLog("atom does not belong to Atom Collection 1 and 2\n");
		return FALSE;
	}

	Pt3CrdRule* p_rule = new Pt3CrdRule( aptr_mng, aref_1, aref_2, aref_3);
	if( dist > 0.0 ) p_rule->SetBondLen(dist);
	if( vang > 0.0 ) p_rule->SetValAng(vang);
	if( torsion > 0.0) p_rule->SetDihAng( torsion );
	if( priority > 0 ) p_rule->SetPriority( priority );

	if( idx_ac == 1 )
	{
		SyncRules1from2.push_back(p_rule);
	}
	else
	{
		SyncRules2from1.push_back(p_rule);
	}
	return TRUE;
}

int AtomMapping::Map2to1ByAtomDistance()
{
	ClearSyncRules2from1();
	std::auto_ptr<AtomIterator> p_aitr_2(p_ac_2->GetAtomIteratorPtr());

	double xmin,xmax,ymin,ymax,zmin,zmax;

	p_ac_2->GetMinMaxCrd(xmin,ymin,zmin,xmax,ymax,zmax);

	BoxPartition part;

	part.SetDimensions(20,20,20);
	part.SetBoundaries(xmin - 5.0,ymin - 5.0,zmin - 5.0,xmax + 5.0,ymax + 5.0,zmax + 5.0);
	part.SetRegionRad(5.0);

	part.DistributePointsToCells( *p_ac_1 );

	std::set<HaAtom* ,less<HaAtom*> > atoms_mapped; 
	atmap_2to1.clear();
	HaAtom* aptr_2;
	
// Find Atoms having identical positions in Atom Colections 1 and 2
	
	for( aptr_2 = p_aitr_2->GetFirstAtom(); aptr_2; aptr_2 = p_aitr_2->GetNextAtom())
	{
		int idx_cell = part.GetPointCellIdx(aptr_2);
		
		VecPtr& pts_cell = part[idx_cell];
        int np = pts_cell.size();
		int i;
		for( i=0; i < np; i++)
		{
			HaAtom* aptr_1 = (HaAtom*) pts_cell[i];
			double dist2 = Vec3D::CalcDistanceSq(aptr_1,aptr_2);
			if( dist2 < 0.0001)
			{
				atmap_2to1[aptr_2] = aptr_1;
				atoms_mapped.insert(aptr_2);
				break;
			}
		}
	}

// Find atoms in Atom Collection 2 that terminate bonds found in Atom Collection 1
		
	for( aptr_2 = p_aitr_2->GetFirstAtom(); aptr_2; aptr_2 = p_aitr_2->GetNextAtom() )
	{
		if( atoms_mapped.count(aptr_2) > 0 ) continue;
		AtomGroup bonded_atoms;
		aptr_2->GetBondedAtoms(bonded_atoms);
		
		AtomIteratorAtomGroup aitr_b(&bonded_atoms);
		
		int found_mapping = FALSE;
		HaAtom* aptr_b;
		for(aptr_b = aitr_b.GetFirstAtom(); aptr_b; aptr_b = aitr_b.GetNextAtom())
		{
			if( atmap_2to1.count(aptr_b) == 0) continue;
			HaAtom* matched_atom = atmap_2to1[aptr_b]; // atom matching atom bonded to aptr_2
			AtomGroup bonded_atoms_2;
			matched_atom->GetBondedAtoms(bonded_atoms_2);
			AtomIteratorAtomGroup aitr_g(&bonded_atoms_2);
			HaAtom* aptr_g;
			
			for(aptr_g = aitr_g.GetFirstAtom(); aptr_g; aptr_g = aitr_g.GetNextAtom())
			{
				double val_ang = Vec3D::CalcAngle(aptr_2, matched_atom, aptr_g);
				if( val_ang < 0.05)
				{
					Pt2CrdRule* p_sync_rule = new Pt2CrdRule(aptr_2,matched_atom, aptr_g);
					p_sync_rule->SetBondLen(Vec3D::CalcDistance(aptr_2,matched_atom));
					p_sync_rule->SetPriority(1);
					SyncRules2from1.push_back(p_sync_rule);
					atoms_mapped.insert(aptr_2);
					found_mapping = TRUE;
					break;
				}
			}
			if( found_mapping) break;
		}
	}

	if( atoms_mapped.size() == p_ac_2->GetNAtoms() )
	{
		PrintLog("All Atoms of the child were matched to atoms of the parent \n");
		PrintLog("As coinciding atoms and atoms in the specified direction \n");
		return TRUE;
	}

// Find coordinate synchronization rules for atoms based coordinates of already mapped atoms
	AtomGroup mapped_atoms_arr;
	std::set<HaAtom*, less<HaAtom*> >::iterator aitr_m;
	for(aitr_m = atoms_mapped.begin(); aitr_m != atoms_mapped.end(); aitr_m++)
	{
		mapped_atoms_arr.push_back(*aitr_m);
	}

	std::vector<Pt3CrdRule*> add_rules;
	int all_atoms_matched = BuildSyncRulesForMissingAtoms(*p_ac_2, mapped_atoms_arr, add_rules);

	int nr = add_rules.size();
	int ir;
	for(ir = 0; ir < nr; ir++)
	{
		Pt3CrdRule* p_rule = add_rules[ir];
		HaAtom* aptr_mng = p_rule->GetManagedAtom();
		atoms_mapped.insert(aptr_mng);
		int ip = p_rule->GetPriority();
		p_rule->SetPriority(ip+2);
		SyncRules2from1.push_back(p_rule);
	}
	
	if( atoms_mapped.size() == p_ac_2->GetNAtoms() )
	{
		PrintLog("All Atoms of Atom Collection 2 have been matched to atoms of Atom Collection 1 \n");
	}
	else
	{
		PrintLog("Not All Atoms of the Atom Collection 2 were matched to atoms of Atom Collection 1 \n");
		return FALSE;
	}

	return TRUE;	
}

int AtomMapping::Map2to1ByAtomRef()
{
	ClearSyncRules2from1();
    std::auto_ptr<AtomIterator> p_aitr_1(p_ac_1->GetAtomIteratorPtr());
	std::auto_ptr<AtomIterator> p_aitr_2(p_ac_2->GetAtomIteratorPtr());

	AtomGroup mapped_atoms;
	map<std::string, HaAtom*, less<std::string> > ref_atom_map_1;
	HaAtom* aptr;
	for( aptr = p_aitr_1->GetFirstAtom(); aptr; aptr = p_aitr_1->GetNextAtom() )
	{
		std::string at_ref = aptr->GetRef();
		if( ref_atom_map_1.count(at_ref) > 0 ) 
		{
			PrintLog("Error in AtomMapping::Map2to1ByAtomRef() \n");
			PrintLog("Atom references in Atom Collection 1 are not unique \n");
			return FALSE;
		}
		ref_atom_map_1[at_ref] = aptr;
	}

	for( aptr = p_aitr_2->GetFirstAtom(); aptr; aptr = p_aitr_2->GetNextAtom() )
	{
		std::string at_ref = aptr->GetRef();
		if( ref_atom_map_1.count(at_ref) > 0 ) 
		{
			atmap_2to1[aptr] = ref_atom_map_1[at_ref];
			mapped_atoms.push_back(aptr);
		}
	}

	if( mapped_atoms.GetNAtoms() == p_ac_2->GetNAtoms() )
	{
		PrintLog("All Atoms of Atom Collection 2 have been matched to atoms of Atom Collection 1 with the same references \n");
		return TRUE;
	}

	std::vector<Pt3CrdRule*> add_rules;
	int all_atoms_matched = BuildSyncRulesForMissingAtoms(*p_ac_2, mapped_atoms, add_rules);

	int nr = add_rules.size();
	int ir;
	for(ir = 0; ir < nr; ir++)
	{
		Pt3CrdRule* p_rule = add_rules[ir];
		HaAtom* aptr_mng = p_rule->GetManagedAtom();
		mapped_atoms.push_back(aptr_mng);
		int ip = p_rule->GetPriority();
		p_rule->SetPriority(ip+2);
		SyncRules2from1.push_back(p_rule);
	}
	
	if( mapped_atoms.size() == p_ac_2->GetNAtoms() )
	{
		PrintLog("All Atoms of Atom Collection 2 have been matched to atoms of Atom Collection 1 \n");
	}
	else
	{
		PrintLog("Not All Atoms of the Atom Collection 2 were matched to atoms of Atom Collection 1 \n");
		return FALSE;
	}
	return TRUE;
}

int AtomMapping::BuildSyncRulesForMissingAtoms(AtomContainer& all_atoms, AtomContainer& known_atoms, 
		                                       std::vector<Pt3CrdRule*>& rules)
{
	rules.clear();

	std::set<HaAtom*> mapped_atoms, mapped_atoms_orig;
	std::auto_ptr<AtomIterator> p_aitr_known( known_atoms.GetAtomIteratorPtr() );
	HaAtom* aptr;
	for( aptr = p_aitr_known->GetFirstAtom(); aptr; aptr = p_aitr_known->GetNextAtom() )
	{
		mapped_atoms.insert(aptr);
	}
	
	mapped_atoms_orig = mapped_atoms;

	std::map<HaAtom*,AtomGroup, less<HaAtom*> > all_nb;      // all neigbors(bonded atoms) of the atom
	std::map<HaAtom*,AtomGroup, less<HaAtom*> > mapped_nb;   // mapped neigbors of the atom
	std::map<HaAtom*,AtomGroup, less<HaAtom*> > unmapped_nb; // unmapped neigbors of the atom
	
	std::auto_ptr<AtomIterator> p_aitr_all( all_atoms.GetAtomIteratorPtr() );
	
	for( aptr = p_aitr_all->GetFirstAtom(); aptr; aptr = p_aitr_all->GetNextAtom() )
	{
		AtomGroup empty_grp;
		all_nb[aptr]    = empty_grp;
		mapped_nb[aptr]   = empty_grp;
		unmapped_nb[aptr] = empty_grp;
		AtomGroup bonded_atoms;
		aptr->GetBondedAtoms(bonded_atoms);
		AtomIteratorAtomGroup aitr_b(&bonded_atoms);
		HaAtom* aptr_b;
		for( aptr_b = aitr_b.GetFirstAtom(); aptr_b; aptr_b = aitr_b.GetNextAtom())
		{
			all_nb[aptr].push_back(aptr_b);
			if( mapped_atoms.count(aptr_b) > 0 )
			{
				mapped_nb[aptr].push_back(aptr_b);
			}
			else
			{
				unmapped_nb[aptr].push_back(aptr_b);
			}
		}
	}
	for(;;)
	{
		int found_cnt_t1 = FALSE; // Process type 1 centers having 2 mapped neigbors and 1 or more unmapped atoms
		for( aptr = p_aitr_all->GetFirstAtom(); aptr; aptr = p_aitr_all->GetNextAtom() )
		{
			if( mapped_atoms.count(aptr) == 0 ) continue;
			if( unmapped_nb[aptr].size() == 0 ) continue;
			if( mapped_nb[aptr].size() < 2 ) continue;
			
			found_cnt_t1 = TRUE;
			int n = unmapped_nb[aptr].size();
			int i;
			for( i = 0; i < n; i++)
			{
				HaAtom* p_mng_atom = unmapped_nb[aptr][i];
				HaAtom* p_ref_1    = aptr;
				HaAtom* p_ref_2    = mapped_nb[aptr][0];
				HaAtom* p_ref_3    = mapped_nb[aptr][1];
				Pt3CrdRule* p_rule = new Pt3CrdRule(p_mng_atom,p_ref_1,p_ref_2,p_ref_3);
				rules.push_back(p_rule);
				mapped_atoms.insert(p_mng_atom);
				mapped_nb[aptr].push_back(p_mng_atom);
			}
			unmapped_nb[aptr].clear();
		}
		if( found_cnt_t1 ) continue;

		int found_cnt_t2 = FALSE; // Now process centers of type 2 having an unmapped neighbor, one mapped neighbor that in turn has at least 2 mapped neghbors
		for( aptr = p_aitr_all->GetFirstAtom(); aptr; aptr = p_aitr_all->GetNextAtom() )
		{
			if( mapped_atoms.count(aptr) == 0 ) continue;
			if( unmapped_nb[aptr].size() == 0 ) continue;
			if( mapped_nb[aptr].size() == 0 ) continue;
			
			HaAtom* aptr1 = mapped_nb[aptr][0];
			if( mapped_nb[aptr1].size() < 2) continue;
			found_cnt_t2 = TRUE;
			HaAtom* aptr2 = mapped_nb[aptr1][0];
			if( aptr2 == aptr ) aptr2 = mapped_nb[aptr1][1];

			int n = unmapped_nb[aptr].size();
			int i;
			for( i = 0; i < n; i++)
			{
				HaAtom* p_mng_atom = unmapped_nb[aptr][i];
				HaAtom* p_ref_1    = aptr;
				HaAtom* p_ref_2    = aptr1;
				HaAtom* p_ref_3    = aptr2;
				Pt3CrdRule* p_rule = new Pt3CrdRule(p_mng_atom,p_ref_1,p_ref_2,p_ref_3);
				rules.push_back(p_rule);
				mapped_atoms.insert(p_mng_atom);
				mapped_nb[aptr].push_back(p_mng_atom);
			}
			unmapped_nb[aptr].clear();
		}
		if( found_cnt_t2 ) continue;

		int found_cnt_t3 = FALSE; // Now process centers of type 3 having unmapped neighbor, one mapped neighbor that doesn't have any mapped neghbor
		for( aptr = p_aitr_all->GetFirstAtom(); aptr; aptr = p_aitr_all->GetNextAtom() )
		{
			if( mapped_atoms.count(aptr) == 0 ) continue;
			if( unmapped_nb[aptr].size() == 0 ) continue;
			if( mapped_nb[aptr].size()   == 0 ) continue;
			found_cnt_t3 = TRUE;
				
			HaAtom* aptr1 = mapped_nb[aptr][0];
			int n = unmapped_nb[aptr].size();
			int i;
			for( i = 0; i < n; i++)
			{
				HaAtom* p_mng_atom = unmapped_nb[aptr][i];
				HaAtom* p_ref_1    = aptr;
				HaAtom* p_ref_2    = aptr1;
				HaAtom* p_ref_3;
				if( i == 0 )
				{
					p_ref_3 = NULL;
				}
				else
				{
					p_ref_3 = unmapped_nb[aptr][0];
				}
				Pt3CrdRule* p_rule = new Pt3CrdRule(p_mng_atom,p_ref_1,p_ref_2,p_ref_3);
				rules.push_back(p_rule);
				mapped_atoms.insert(p_mng_atom);
				mapped_nb[aptr].push_back(p_mng_atom);
			}
			unmapped_nb[aptr].clear();
		}
		if( found_cnt_t3 ) continue;
		
		int found_cnt_t4 = FALSE; // Now process centers of type 4 having only unmapped neighbors 
		for( aptr = p_aitr_all->GetFirstAtom(); aptr; aptr = p_aitr_all->GetNextAtom() )
		{
			if( mapped_atoms.count(aptr) == 0 ) continue;
			if( unmapped_nb[aptr].size() == 0 ) continue;
			
			found_cnt_t4 = TRUE;
				
			int n = unmapped_nb[aptr].size();
			int i;
			for( i = 0; i < n; i++)
			{
				HaAtom* p_mng_atom = unmapped_nb[aptr][i];
				HaAtom* p_ref_1    = aptr;
				HaAtom* p_ref_2;
				HaAtom* p_ref_3;
				if( i == 0 ) p_ref_2 = NULL;
				else p_ref_2 = unmapped_nb[aptr][0];
				
				if( i < 2 ) p_ref_3 = NULL;
				else p_ref_3 = unmapped_nb[aptr][1];

				Pt3CrdRule* p_rule = new Pt3CrdRule(p_mng_atom,p_ref_1,p_ref_2,p_ref_3);
				rules.push_back(p_rule);
				mapped_atoms.insert(p_mng_atom);
				mapped_nb[aptr].push_back(p_mng_atom);
			}
			unmapped_nb[aptr].clear();
		}
		if( found_cnt_t4 ) continue;
		break;
	}

	int all_atoms_mapped = TRUE;
	for( aptr = p_aitr_known->GetFirstAtom(); aptr; aptr = p_aitr_known->GetNextAtom() )
	{
		if( mapped_atoms.count(aptr) == 0 ) all_atoms_mapped = FALSE;
	}

	if( !all_atoms_mapped )return FALSE;
	
	return TRUE;
	
	//int nseq = 0;
	//int cur_priority = 1;
	//for(;;)
	//{
	//	int found_unmapped_atom = FALSE;
	//	HaAtom* aptr_ch;
	//	for( aptr_ch = p_aitr_all->GetFirstAtom(); aptr_ch; aptr_ch = p_aitr_all->GetNextAtom() )
	//	{
	//		if( mapped_atoms.count(aptr_ch) == 0 ) continue;
	//		AtomGroup bonded_atoms;
	//		aptr_ch->GetBondedAtoms(bonded_atoms);
	//		AtomIteratorAtomGroup aitr_g(&bonded_atoms);
	//		HaAtom* aptr_g;
	//		
	//		AtomGroup bonded_atoms_mapped;
	//	    AtomGroup bonded_atoms_unmapped;

	//		for( aptr_g = aitr_g.GetFirstAtom(); aptr_g; aptr_g = aitr_g.GetNextAtom())
	//		{
	//			if( mapping_info.count(aptr_ch) > 0) 
	//			{
	//				bonded_atoms_mapped.InsertAtom(aptr_ch);
	//			}
	//			else
	//			{
	//				bonded_atoms_unmapped.InsertAtom(aptr_ch);
	//			}
	//		}

	//		if( bonded_atoms_unmapped.size() == 0 ) continue;

	//		HaAtom* aref_2 = NULL;
	//		HaAtom* aref_3 = NULL;

	//		AtomIteratorAtomGroup aitr_bm (&bonded_atoms_mapped);

	//		HaAtom* aptr_bm;
	//		aptr_bm = aitr_bm.GetFirstAtom();

	//		if( aptr_bm != NULL )  
	//		{
	//			aref_2  = aptr_bm;
	//			aptr_bm = aitr_bm.GetNextAtom();
	//			if( aptr_bm != NULL )
	//			{
	//				aref_3 = aptr_bm;
	//			}
	//		}

	//		if( aref_3 == NULL )
	//		{
	//			if( aref_2 == NULL) continue;
	//		
	//			AtomGroup bonded_atoms_2;
	//			aref_2->GetBondedAtoms(bonded_atoms_2);
	//			bonded_atoms_2.DeleteAtom(aptr_ch);

	//			AtomIteratorAtomGroup aitr_g2(&bonded_atoms_2);
	//			
	//			HaAtom* aptr_2;
	//			for( aptr_2 = aitr_g2.GetFirstAtom(); aptr_2; aptr_2 = aitr_g2.GetNextAtom())
	//			{
	//				if( mapping_info.count(aptr_2) > 0) 
	//				{
	//					aref_3 = aptr_2;
	//					break;
	//				}
	//			}
	//		}	

	//		if(aref_3 == NULL) continue;

	//		AtomIteratorAtomGroup aitr_um(&bonded_atoms_unmapped);
	//		HaAtom* aptr_um;

	//		int first_atom = TRUE;
	//		for( aptr_um = aitr_um.GetFirstAtom(); aptr_um; aptr_um = aitr_um.GetNextAtom())
	//		{
 //               AssociateAtomsRefAtoms(aptr_um, aptr_ch, aref_2, aref_3, nseq);
	//			found_unmapped_atom = TRUE;
	//			nseq++;
	//			if( first_atom )
	//			{
	//				first_atom = FALSE;
	//				aref_3 = aptr_um;
	//			}
	//		}
	//	}	
	//	if(!found_unmapped_atom) break;
	//}

	//all_atoms_matched = TRUE;

	//for( aptr_ch = aitr_ch.GetFirstAtom(); aptr_ch; aptr_ch = aitr_ch.GetNextAtom())
	//{
	//	if( mapping_info.count(aptr_ch) == 0 ) 
	//	{
	//		std::string atom_ref = aptr_ch->GetRef();
	//		PrintLog("Warning in AtomMapping::MapChildAtoms() \n");
	//		PrintLog("No matching was foung for atom %s ", atom_ref.c_str());
	//		all_atoms_matched = FALSE;
	//	}
	//}
	//
	//if( all_atoms_matched )
	//{
	//	PrintLog("All Atoms of the child were matched to atoms of the parent \n");
	//}
	//else
	//{
	//	PrintLog("Not All Atoms of the child were matched to atoms of the parent \n");
	//	return FALSE;
	//}
}

bool compare_sync_rules(CrdAssignRule* p_rule_1,CrdAssignRule* p_rule_2) { return ( p_rule_1->GetPriority() < p_rule_2->GetPriority()); }

int AtomMapping::SyncAtomCrd1From2()
{
	std::map< HaAtom*, HaAtom*, less<HaAtom*> >::iterator aitr_m;
	for( aitr_m = atmap_1to2.begin() ; aitr_m != atmap_1to2.end(); aitr_m++)
	{
		HaAtom* aptr_mng = (*aitr_m).first;
		HaAtom* aptr_ref = (*aitr_m).second;
		aptr_mng->SetCoordFrom(*aptr_ref);
	}

//	std::sort(SyncRules1from2.begin(), SyncRules1from2.end(), compare_sync_rules); 

	int nr = SyncRules1from2.size();
	int ir;
	for(ir = 0; ir < nr; ir++ )
	{
		SyncRules1from2[ir]->SetManagedAtomCrd();
	}
	
	return TRUE;
}

int AtomMapping::SyncAtomCrd2From1()
{
	std::map< HaAtom*, HaAtom*, less<HaAtom*> >::iterator aitr_m;
	for( aitr_m = atmap_2to1.begin() ; aitr_m != atmap_2to1.end(); aitr_m++)
	{
		HaAtom* aptr_mng = (*aitr_m).first;
		HaAtom* aptr_ref = (*aitr_m).second;
		aptr_mng->SetCoordFrom(*aptr_ref);
	}

//	std::sort(SyncRules2from1.begin(), SyncRules2from1.end(), compare_sync_rules); 

	int nr = SyncRules2from1.size();
	int ir;
	for(ir = 0; ir < nr; ir++ )
	{
		SyncRules2from1[ir]->SetManagedAtomCrd();
	}
	return TRUE;
}

void AtomMapping::PrintInfo(int detailed)
{
	int na1 = p_ac_1->GetNAtoms();
	int na2 = p_ac_2->GetNAtoms();
	int n_map_dir_2 = atmap_2to1.size();

	set<HaAtom*, less<HaAtom*> > unmapped_atoms;
	std::auto_ptr<AtomIterator> p_aitr_2( p_ac_2->GetAtomIteratorPtr() );

	HaAtom* aptr2;
	for( aptr2 = (*p_aitr_2).GetFirstAtom(); aptr2; aptr2 = (*p_aitr_2).GetNextAtom() )
	{
		unmapped_atoms.insert(aptr2);
	}
	
	std::map< HaAtom*, HaAtom*, less<HaAtom*> >::iterator mitr;
	for( mitr = atmap_2to1.begin(); mitr != atmap_2to1.end(); mitr++)
	{
		aptr2 = (*mitr).first;
		unmapped_atoms.erase(aptr2);
	}

	std::vector<SingleAtomCrdRule*>::iterator ritr;
	for( ritr = SyncRules2from1.begin(); ritr != SyncRules2from1.end(); ritr++)
	{
		aptr2 = (*ritr)->GetManagedAtom();
		unmapped_atoms.erase(aptr2);	
	}

	PrintLog(" Tot Number of Atoms in Collections 1 and 2 are : %d  %d \n", na1, na2 );
	PrintLog(" Number of directly mapped atoms of Collection 2: %d \n", n_map_dir_2);
	PrintLog(" Number of indirectly mapped atoms of Collection 2: %d \n", SyncRules2from1.size() );
	PrintLog(" Number of unmapped atoms of Collection 2: %d \n", unmapped_atoms.size() );	

	if( detailed )
	{
		PrintLog("\nUnmapped atoms: |n");
		set<HaAtom*, less<HaAtom*> >::iterator aitr_u;
		for( aitr_u = unmapped_atoms.begin(); aitr_u != unmapped_atoms.end(); aitr_u++)
		{
			PrintLog( "%s \n",(*aitr_u)->GetRef().c_str() );
		}
	}
}

