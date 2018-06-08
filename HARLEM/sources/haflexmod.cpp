#include <mpi.h>
#include "haflexmod.h"
#include <iostream>
#include <cstdio>

#include <wx/string.h>
#include <wx/filename.h>
#include <wx/process.h>

#include <assert.h>
#include <float.h>
#include <math.h>
#include "haio.h"
#include "tokens.h"
#include "harlemapp.h"
#include "haatgroup.h"
#include "hamolmech.h"
#include "mm_traj_anal.h"
#include "haatom.h"
#include "haatgroup.h"
#include "hamolset.h"
#include "hamolecule.h"
#include "moleditor.h"
#include "haresdb.h"
#include "hahbhp.h"
#include "canvas3d.h"
#include "hamolview.h"


HaFlexMod::HaFlexMod(HaMolSet* new_phost_mset):
HaCompMod(COMP_MOD_FLEX,new_phost_mset)
{
	p_mol_editor = new MolEditor();
	Init();
}

HaFlexMod::~HaFlexMod()
{
	if(p_mol_editor) delete p_mol_editor;

	ClearHB();
	ClearHPT();
}

void HaFlexMod::ClearHB()
{
	int nhb = HBArray.size();
	int i;
	for( i = 0; i < nhb; i++)
	{
		delete HBArray[i];
	}
	HBArray.clear();
}

void HaFlexMod::ClearHPT()
{
	int nph = HPTArray.size();
	int i;
	for( i = 0; i < nph; i++)
	{
		delete HPTArray[i];
	}
	HPTArray.clear();
}

void HaFlexMod::Init()
{
	struct_file_name.clear();
	first_data_file_name.clear();

	hydrophob_dist_cutoff = 3.9;
	hb_duty_cycle_cutoff  = 0.0;
	hpt_duty_cycle_cutoff  = 0.0;
	hb_energy_cutoff      = 0.0;
}

void HaFlexMod::AddAllAtomGroup()
{
	HaMolSet* pmset = GetMolSet();
	AtomGroup* atgrp = pmset->GetAtomGroupByID("ALLATOMLIST");
	if( atgrp == NULL) 
	{
		atgrp = pmset->AddAtomGroup("ALLATOMLIST");
	}
	atgrp->clear();
	AtomIteratorMolSet aitr(pmset);
	HaAtom* aptr;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		atgrp->InsertAtom(aptr);
	}	
}


void HaFlexMod::createAtomRefMap()
{
	atom_ref_map.clear();
	AtomIteratorMolSet aitr( this->GetMolSet() );

	HaAtom* aptr; 
	int i = 1;
	
	for(aptr = aitr.GetFirstAtom() ; aptr ; aptr = aitr.GetNextAtom())
	{
		atom_ref_map[aptr] = i;
		i++;
	}
	
	PrintLog("%d entries recorded in map\n", atom_ref_map.size());
	
}

int HaFlexMod::FindHydrogenBonds()
{
	HaMolSet* pmset = this->GetMolSet();
	
	AtomContainer* p_active_atoms = pmset;
	AtomGroup* p_act_at_array = pmset->GetAtomGroupByID(active_at_array_id.c_str());
	if( p_act_at_array != NULL) p_active_atoms = p_act_at_array;

	vector<HaHBond> hbond_arr_axx;
	p_mol_editor->FindHBondsAtomCollection(p_active_atoms,hbond_arr_axx );
	int n_hb = hbond_arr_axx.size();
	this->ClearHB();
	int i;
	for(i = 0; i < n_hb;i++)
	{
		HaHBond& hb_ref = hbond_arr_axx[i];
		HBondAvg* p_hb =  new HBondAvg( hb_ref.GetDonorAtom(), hb_ref.GetAcceptorAtom(), hb_ref.GetHAtom() );
		HBArray.push_back(p_hb);
	}
	
	return TRUE;
}

int HaFlexMod::FindHydrophobicContacts()
{
	vector<HaAtom*> donorAtomList;
	vector<HaAtom*> acceptorAtomList;
	vector<HaAtom*> donorAcceptorList;
	vector<double> donorAcceptorDistanceList;
	
	ClearHPT();
	selected_hpt.clear();

	PrintLog(" Hydrophobic distance cutoff %8.3f \n", hydrophob_dist_cutoff);
	
	AtomIteratorMolSet aitr_mset(phost_mset);
	HaAtom* iter;
	AtomGroup bonded_atoms;

	int nat = phost_mset->GetNAtoms();
	vector<HaAtom*>    np_atoms;
	np_atoms.reserve(nat);
	HaAtom* aptr;

	for(aptr = aitr_mset.GetFirstAtom(); aptr; aptr = aitr_mset.GetNextAtom())
	{
		int elem = aptr->GetElemNo();

		if( elem != 6 && elem != 16) continue;
		HaResidue* pres = aptr->GetHostRes();
		std::string res_name = pres->GetName();
		if( elem == 16 )
		{
			if( res_name == "MET" )
			{
				np_atoms.push_back(aptr);
			}
			continue;
		}

		std::string at_name = aptr->GetName();
		if( pres->IsAmino() && at_name == "CA")
		{
			np_atoms.push_back(aptr);
			continue;
		}
		
		aptr->GetBondedAtoms(bonded_atoms);
		AtomIteratorAtomGroup aitr_grp(&bonded_atoms);
		HaAtom* aptr2;
		for(aptr2 = aitr_grp.GetFirstAtom(); aptr2; aptr2 = aitr_grp.GetNextAtom())
		{
			if(aptr2->GetElemNo() != 6 && aptr2->GetElemNo() != 16 && aptr2->GetElemNo() != 1) break;
		}
		if( aptr2 == NULL) np_atoms.push_back(aptr);
	}		

	int nat_np = np_atoms.size();
	int i,j;

	double ph_dist_2 = hydrophob_dist_cutoff*hydrophob_dist_cutoff;

	for(i = 0; i < nat_np; i++)
	{
		HaAtom* aptr1 = np_atoms[i];
		for( j = 0; j < i; j++ )
		{
			HaAtom* aptr2 = np_atoms[j];
			if( aptr1->GetHostRes() == aptr2->GetHostRes() ) continue;
			Vec3D diff_v = (*aptr1) - (*aptr2);
			double dist2 = diff_v.length2();
			if(dist2 > ph_dist_2) continue;

			donorAtomList.push_back(aptr1);
			acceptorAtomList.push_back(aptr2);
			donorAcceptorDistanceList.push_back(sqrt(dist2));	
		}
	}

	for(int m = 0; m < donorAtomList.size(); m++) 
	{
		// PrintLog("DONOR      : %s\n", donorAtomList[m]->GetRef());
		// PrintLog("ACCEPTOR   : %s\n", acceptorAtomList[m]->GetRef());
		// HaHydrophobicTether(int f, int s, float minDis, float maxDis, float avgDis, float DC)
		HaHydrophobicTether* p_hpt = new HaHydrophobicTether(donorAtomList[m], acceptorAtomList[m]) ;
		p_hpt->SetAvgDistance(donorAcceptorDistanceList[m]);
		HPTArray.push_back(p_hpt);
	}
	return TRUE;	
}

int HaFlexMod::SaveHBondFile()
{
	if(atom_ref_map.size() == 0) createAtomRefMap();

	FILE* fp = fopen("hbonds.in","w");
	if(fp == NULL) return FALSE;

	int i;
	int n = HBArray.size();

	for(i=0; i < n; i++)
	{
		int id=-1;
		int ia=-1;
		HaAtom* aptr_d = HBArray[i]->GetDonorAtom();
		HaAtom* aptr_a = HBArray[i]->GetAcceptorAtom();
		if( hb_duty_cycle_cutoff > 0.01 && HBArray[i]->GetDutyCycle() < hb_duty_cycle_cutoff) continue;
		if( atom_ref_map.count(aptr_d) ) id = atom_ref_map[aptr_d]; 
		if( atom_ref_map.count(aptr_a) ) ia = atom_ref_map[aptr_a];
		fprintf(fp,"%d %d\n",id,ia);
	}
	fclose(fp);

	return TRUE;
}

	
int HaFlexMod::SaveHydrophobTethersFile()
{
	if(atom_ref_map.size() == 0) createAtomRefMap();

	FILE* fp = fopen("hphobes.in","w");
	if(fp == NULL) return FALSE;

	int i;
	int n = HPTArray.size();

	for(i=0; i < n; i++)
	{
		int i1=-1;
		int i2=-1;
		HaAtom* aptr_1 = HPTArray[i]->GetFirstAtom();
		HaAtom* aptr_2 = HPTArray[i]->GetSecondAtom();
		if( hpt_duty_cycle_cutoff > 0.01 && HPTArray[i]->GetDutyCycle() < hpt_duty_cycle_cutoff) continue;
		if( atom_ref_map.count(aptr_1) ) i1 = atom_ref_map[aptr_1]; 
		if( atom_ref_map.count(aptr_2) ) i2 = atom_ref_map[aptr_2];
		fprintf(fp,"%d %d\n",i1,i2);
	}
	fclose(fp);

	return TRUE;
}

int HaFlexMod::SaveStructFile()
{
	AddAllAtomGroup();
	if(struct_file_name.size() == 0)
	{
		struct_file_name = phost_mset->GetName();
		struct_file_name += "_fst";
		struct_file_name += ".hlm";
	}

	phost_mset->SaveOldHarlemFile(struct_file_name.c_str());
	
	return TRUE;
}

int HaFlexMod::SaveFirstInpFiles()
{
	SaveHBondFile();
	SaveHydrophobTethersFile();
	SaveStructFile();

	return TRUE;
}

int HaFlexMod::ReadFirstDataFile()
{
	char buf[256];
	PrintLog(" HaFlexMod::ReadFirstDataFile() 1 \n");
	if(first_data_file_name.size() == 0)
	{
		first_data_file_name = phost_mset->GetName();
		first_data_file_name += "_fst_data";
		first_data_file_name += ".txt";
	}

	list<AtomGroup>::iterator gitr = phost_mset->NamedAtomGroups.begin();
	
	for(; gitr != phost_mset->NamedAtomGroups.end(); )
	{
		std::string grp_name = (*gitr).GetID();
		if( grp_name.find("RIGIDCLUST") != std::string::npos )  
		{
			gitr = phost_mset->NamedAtomGroups.erase(gitr);
		}
		else
		{
			gitr++;
		}
	}

	ifstream fin(first_data_file_name.c_str());
	
	fin.getline (buf,256); 
	fin.getline (buf,256); 
	fin.getline (buf,256); 
	
	if( fin.eof() ) return FALSE;

	vector<AtomGroup> grp_vec;
	vector<HaAtom*> atom_vec;
	atom_vec.reserve(phost_mset->GetNAtoms());

	AtomIteratorMolSet aitr(phost_mset);
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr != NULL; aptr = aitr.GetNextAtom())
	{
		atom_vec.push_back(aptr);
	}
	
	for(;;)
	{	
		fin.getline(buf,256);
		if( fin.eof() ) break;
		int fst_num = 0;
		int orig_num = 0;
		int rigid_num = 0;
		int stressed_num = 0;
		int coll_mode = 0;
		
		istrstream line_s(buf);
		line_s >> fst_num; 
		line_s >> orig_num;
		line_s >> rigid_num;
		line_s >> stressed_num;
		line_s >> coll_mode;

		if( orig_num == 0 ) break; 
		
		if( grp_vec.size() < rigid_num ) grp_vec.resize(rigid_num);
		if( orig_num > atom_vec.size()) break;
		aptr = atom_vec[ orig_num - 1];
		grp_vec[rigid_num - 1].InsertAtom(aptr);
	}
	
	int ng = grp_vec.size();
	
	int i;
	for( i= 0; i < ng; i++)
	{
		std::string grp_name = "RIGIDCLST";
		sprintf(buf,"%d",i+1);
		grp_name += buf;
		grp_vec[i].SetID(grp_name.c_str());
		phost_mset->NamedAtomGroups.push_back(grp_vec[i]);
	}

	return TRUE;
}

int HaFlexMod::DeleteSelectedHBonds()
{
	int nsel = selected_hb.size();
	if( nsel == 0) return TRUE;
    
	int i;

	vector<HBondAvg*>::iterator bitr = HBArray.begin();
	for(; bitr != HBArray.end(); )
	{
		int found = FALSE;
		HBondAvg* p_hb = *bitr;
		for( i = 0; i < nsel; i++)
		{
			if( selected_hb[i] == p_hb ) found = TRUE;
		}
		if(found)
		{
			bitr = HBArray.erase(bitr);
		}
		else
		{
			bitr++;
		}
	}
	selected_hb.clear();
	return TRUE;
}

int HaFlexMod::DeleteSelectedHPTethers()
{
	int nsel = selected_hpt.size();
	if( nsel == 0) return TRUE;
    
	int i;

	vector<HaHydrophobicTether*>::iterator bitr = HPTArray.begin();
	for(; bitr != HPTArray.end(); )
	{
		int found = FALSE;
		HaHydrophobicTether* p_hpt = *bitr;
		for( i = 0; i < nsel; i++)
		{
			if( selected_hpt[i] == p_hpt ) found = TRUE;
		}
		if(found)
		{
			bitr = HPTArray.erase(bitr);
		}
		else
		{
			bitr++;
		}
	}
	selected_hpt.clear();
	return TRUE;
}

int HaFlexMod::ComputeBondsDutyCycleAlongMD()
{
//	PrintLog(" HaFlexMod::ComputeBondsDutyCycleAlongMD() pt 1 \n");
	HaMolSet* pmset = this->GetMolSet();
	int ires;
	MDTrajectory trj(pmset);
	trj.CrdFileName = md_traj_file_name.c_str();
	ires = trj.Open();
	if(!ires) 
	{
		PrintLog("Error to open Md Trajectory file %s ",trj.CrdFileName.c_str() );
		return FALSE;
	}

	double ph_dist_cutoff_save = hydrophob_dist_cutoff; 
	double hb_dist_cutoff_save = p_mol_editor->max_da_dist_no_acc;
	double hb_ang_cutoff_save  = p_mol_editor->max_hda_angle;

	hydrophob_dist_cutoff += 2.0;
	p_mol_editor->max_da_dist_no_acc += 2.0;
	p_mol_editor->max_hda_angle += 0.4;

// Find hydrogen bonds and hydrophobic tethers with looser criteria   

	FindHydrogenBonds();
	FindHydrophobicContacts();

// Restore hydrogen bond and hydrophobic tethers criteria 

	hydrophob_dist_cutoff = ph_dist_cutoff_save;
	p_mol_editor->max_da_dist_no_acc  = hb_dist_cutoff_save;
	p_mol_editor->max_hda_angle = hb_ang_cutoff_save;

	int n_hb = HBArray.size();
	int n_ph = HPTArray.size();

	vector<int> count_hb(n_hb,0);
	vector<int> count_ph(n_ph,0);

	int n_pt = 0;
	int i;

	for(;;)
	{
		int ires = trj.ReadNextFrame();
		if(!ires) break;

		n_pt++;
		PrintLog("MD point %d \n");

		for(i = 0; i < n_hb; i++)
		{
			if( p_mol_editor->IsValidHBond(HBArray[i]) )
				count_hb[i]++;
		}

		for(i = 0; i < n_ph; i++)
		{
			HaHydrophobicTether* pph = HPTArray[i];
			HaAtom* aptr_1 = pph->GetFirstAtom();
			HaAtom* aptr_2 = pph->GetSecondAtom();
			
			double dist = Vec3D::CalcDistance(aptr_1,aptr_2);

			if( dist < hydrophob_dist_cutoff)
			{
				count_ph[i]++;
			}
		}
	}

	PrintLog("Number of MD points read %d \n", n_pt);
	
	if(n_pt == 0) return FALSE;

	vector<HBondAvg*>::iterator hb_itr;
	i = 0;
	for(hb_itr = HBArray.begin(); hb_itr != HBArray.end(); )
	{
		HBondAvg* p_hb = *hb_itr;
		int ic = count_hb[i];
		if( ic > 0)
		{
			double duty_cycle = (double)ic/(double)(n_pt);
			p_hb->SetDutyCycle(duty_cycle);
			hb_itr++;
		}
		else
		{
			delete p_hb;
			hb_itr = HBArray.erase(hb_itr);
		}
		i++;
	}

	vector<HaHydrophobicTether*>::iterator hpt_itr;
	i = 0;
	for(hpt_itr = HPTArray.begin(); hpt_itr != HPTArray.end(); )
	{
		HaHydrophobicTether* p_hpt = *hpt_itr;
		int ic = count_ph[i];
		if( ic > 0)
		{
			double duty_cycle = (double)ic/(double)(n_pt);
			p_hpt->SetDutyCycle(duty_cycle);
			hpt_itr++;
		}
		else
		{
			delete p_hpt;
			hpt_itr = HPTArray.erase(hpt_itr);
		}
		i++;
	}

	return TRUE;
}

int HaFlexMod::Display(HaMolView* pView)
{
	int nhb = HBArray.size();
	int nph = HPTArray.size();

	HaColor yellow_col(255,255,0); // yellow
    HaColor green_col(0,255,0);   // green

	int irad = (int) (pView->Scale* 0.2); // radius of the cylinder 

	int i;
	int j;
	HaAtom* aptr1;
	HaAtom* aptr2;

	

	for(i = 0; i < nhb; i++)
	{
		aptr1 = HBArray[i]->GetDonorAtom();
		aptr2 = HBArray[i]->GetAcceptorAtom();
		if( aptr1->Selected() && aptr2->Selected())
		{
			pView->pCanv->ClipCylinder(aptr1->x, aptr1->y, aptr1->z, aptr2->x, aptr2->y, aptr2->z,
						               green_col.cidx,green_col.cidx,irad);
		}
	}
		
	for(i = 0; i < nph; i++)
	{
		aptr1 = HPTArray[i]->GetFirstAtom();
		aptr2 = HPTArray[i]->GetSecondAtom();
		if( aptr1->Selected() && aptr2->Selected())
		{
			pView->pCanv->ClipCylinder(aptr1->x, aptr1->y, aptr1->z, aptr2->x, aptr2->y, aptr2->z,
					                   yellow_col.cidx,yellow_col.cidx,irad);
		}
	}
	return TRUE;
}
