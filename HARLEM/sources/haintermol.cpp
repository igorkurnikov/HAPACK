/*  \file haintermol.cpp

    Classes to model intermolecular interactions in HARLEM  
 
    \author Igor Kurnikov 
    \date 1999-2002
*/

#define HAINTERMOL_CPP

#include <mpi.h>

#include <chrono>
#include <thread>

#include <float.h>
#include <math.h>
#include "hampi.h"

#include "rapidxml.hpp"

#include "harlemapp.h"
#include "tokens.h"
#include "hacoord.h"
#include "haatgroup.h"
#include "rigidbodycoord.h"
#include "haenefunc.h"
#include "trajanal.h"
#include "hasimulator.h"
#include "haintermol.h"
#include "hamolset.h"
#include "hamolmech.h"
#include "mm_driver_amber.h"
#include "hamolecule.h"
#include "electrostmod.h"
#include "protonredox.h"
#include "hamolview.h"
#include "etcoupl.h"
#include "haempirical.h"
#include "protonredox.h"

#include <stdlib.h>
#include <time.h>
#include "mpi.h"
  
#include "randomip.h"

using namespace harlem;

static float FMAX( float a, float b)
{
 if (a > b) return a;
 else return b;		 
}



HaInterMolMod::HaInterMolMod(MolSet* new_phost_mset):
HaCompMod(COMP_MOD_INTERMOL,new_phost_mset)
{
	p_mc_sim        = new InterMolMCSimulator(this);     
	p_ene_minimizer = new InterMolEnergyMinimizer(this);
	p_rex_sim       = new InterMolRepExchSimulator(this);
	
	MolSet* pmset = GetMolSet();
	p_prot_rdx_mod  = pmset->GetProtonRedoxMod(true);

	SetStdParams();
}

HaInterMolMod::~HaInterMolMod()
{
	delete(p_mc_sim);
	delete(p_ene_minimizer);
	delete(p_rex_sim);
}

void HaInterMolMod::SetStdParams()
{
	cur_intermol_ene = 0.0;
	electr_inter_ene = 0.0;
	vdw_inter_ene    = 0.0;

	calc_et_rate = FALSE;
	compute_pk = 0;
	empirical_flag = FALSE;

	to_build_nb_contact_list = TRUE;
	to_build_intermol_excl_atom_list = TRUE;

	electr_model = CHARGES_IN_FIELD_ELECTR;
}	

int HaInterMolMod::Initialize()
{
	int ires;
	ClearInternalStruct();
	if( interact_groups.size() < 2 )
	{
		ires = SetInteractGroupsFromMolecules();
		if( !ires)
		{
			PrintLog("Error in HaInterMolMod::Initialize() \n");
			PrintLog("Can not set Interacting Groups from Molecules \n");
			return FALSE;
		}
	}
		
	MolSet* pmset = GetMolSet();
	HaMolMechMod* pmm_mod = pmset->GetMolMechMod(true);

	if(p_mc_sim->amber_flag)
	{
		pmm_mod->p_amber_driver->SaveAllInpFiles();
	}

	if( electr_model == COULOMB_ELECTR ) 
	{
		pmm_mod->p_mm_model->calc_electr_flag = TRUE;
	}
	else if( electr_model == CHARGES_IN_FIELD_ELECTR ) 
	{
		pmm_mod->p_mm_model->calc_electr_flag = FALSE;
	}
	else 
	{
		pmm_mod->p_mm_model->calc_electr_flag = FALSE; // do not calculate electrostatic interaction 
	                                       // in the MM module
	}
	
	if(calc_et_rate)
	{
		ETCouplMod* ptr_et_coupl_mod= pmset->GetETCouplMod();
		if( ptr_et_coupl_mod == NULL )
		{
			PrintLog("Error In HaInterMolMod::Initialize()\n");
			PrintLog("ET Coupling module is not set \n");
			return FALSE;
		}
		
		AtomGroup* pdon = pmset->GetAtomGroupByID("DONOR");
		AtomGroup* pacc = pmset->GetAtomGroupByID("ACCEPTOR");

		if( !pdon || !pacc || pdon->empty() || pacc->empty() )
		{
			PrintLog("Error in HaInterMolMod::Initialize()");
			PrintLog(" ET Coupling module is not set \n");
			return FALSE;				
		}
	}
	
	return TRUE; 
}

int HaInterMolMod::ClearInternalStruct()
{
	el_pot_field.clear();
	vdw_pot_field.clear();
	return TRUE;
}

int
HaInterMolMod::SetInteractGroupsFromMolecules()
{
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	int imol;
	
	if(pmset->GetNMol() == 1)
	{
		PrintLog("Error in HaInterMolMod::Initialize() \n");
		PrintLog("The number of molecules in the Molecular Set is equal 1 \n");
		return FALSE;
	}

	ClearInteractGroups();

	interact_groups.resize(nmol);
	for( imol= 0; imol < nmol; imol++) 
	{
		HaMolecule* pmol = pmset->HostMolecules[imol];
		interact_groups[imol] = pmol;
	}
	
	return TRUE;
}

int 
HaInterMolMod::ClearInteractGroups()
{
	interact_groups.clear();
	to_build_nb_contact_list = TRUE;
	to_build_intermol_excl_atom_list = TRUE;

	return TRUE;
}

int 
HaInterMolMod::AddInteractGroup(AtomContainer* p_at_cont)
{
	if( p_at_cont == NULL ) 
	{
		PrintLog("Error in HaInterMolMod::AddInteractGroup() \n");
		PrintLog("p_at_cont == NULL \n");
		return FALSE;
	}
	
	MolSet* pmset = this->GetMolSet();
	AtomIteratorGen aitr(p_at_cont);
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if( aptr->GetHostMolSet() != pmset )
		{
			PrintLog("Error in HaInterMolMod::AddInteractGroup() \n");
			PrintLog("An atom of the atom container doesn't belong to the Molecular Set associated with Intermolecular Module \n");
			return FALSE;
		}
	}
	interact_groups.push_back(p_at_cont);
	to_build_nb_contact_list = TRUE;
	to_build_intermol_excl_atom_list = TRUE;

	return TRUE;
}

int HaInterMolMod::SetCoord(harlem::Coord* pcrd)
{
	if( pcrd->GetClassName() == "RigidBodyCoord" )
	{
		this->SetRigidBodyCoord( (RigidBodyCoord*) pcrd);
	}
	return TRUE;
}

int HaInterMolMod::SetRigidBodyCoord(RigidBodyCoord* pcrd)
{
	int ng = interact_groups.size();
	if( ng != pcrd->GetNumObj() )
	{
		throw "Error in HaInterMolMod::SetRigidBodyCoord() \n The number of Interacting groups does not equal to the number of objects in pcrd \n" ;
	}

	int i;
	for(i = 0; i < ng; i++)
	{
		Vec3D trans;
		double phi       = pcrd->GetPhi(i);
		double cos_theta = pcrd->GetCosTheta(i);
		double psi       = pcrd->GetPsi(i);
		trans[0]         = pcrd->GetTransX(i);
		trans[1]         = pcrd->GetTransY(i);
		trans[2]         = pcrd->GetTransZ(i);

		interact_groups[i]->SetPosEulerTrans(phi, cos_theta, psi, trans);
	}
	return TRUE;
}

double HaInterMolMod::CalculateMMEnergy()
{
	MolSet* pmset = GetMolSet();		
	HaMolMechMod* pmm_mod = pmset->GetMolMechMod(true);
	if(interact_groups.size() < 2)
	{
		int ires = SetInteractGroupsFromMolecules();
		if( !ires) 
		{
			PrintLog("Error in HaInterMolMod::CalculateMMEnergy() \n");
			PrintLog("Interacting Atom Groups are not set \n");
			return FALSE;
		}
	}
	if(pmm_mod->p_mm_model->to_init_mm_model )
	{
		pmm_mod->InitMolMechModel();
		to_build_nb_contact_list = TRUE;
		to_build_intermol_excl_atom_list = TRUE;
	}
	if(to_build_intermol_excl_atom_list) 
	{
		pmm_mod->p_mm_model->BuildGrpGrpExcludedList(interact_groups[0],interact_groups[1]);
		to_build_intermol_excl_atom_list = FALSE;
	}
	if(to_build_nb_contact_list)
	{
		pmm_mod->p_mm_model->BuildGrpGrpNonBondList(interact_groups[0],interact_groups[1]);
		to_build_nb_contact_list = FALSE;
	}

	pmm_mod->CalcEnergy();
	
    electr_inter_ene = pmm_mod->p_mm_info->electr_ene_nb;
	vdw_inter_ene    = pmm_mod->p_mm_info->vdw_ene_nb;
	cur_intermol_ene = electr_inter_ene + vdw_inter_ene;
	return cur_intermol_ene;
}

double
HaInterMolMod::CalcElStaticInter()
{	
	if(interact_groups.size() < 2)
	{
		int ires = SetInteractGroupsFromMolecules();
		if( !ires) 
		{ 
			PrintLog("Error in HaInterMolMod::CalcElStaticInter() \n");
			PrintLog("There are less than two interacting groups in the Intermolecular Module \n");
			return 0.0;
		}
	}

	if(interact_groups.size() > 2)
	{
		PrintLog("Error in HaInterMolMod::CalcElStaticInter() \n");
		PrintLog("There are more than two interacting groups in the Intermolecular Module \n");
	    return 0.0;
	}

	double el_inter_ene = 0.0; 

	MolSet* pmset = GetMolSet();		
	ElectrostMod* el_mod = pmset->GetElectrostMod(true);

	if(electr_model == CHARGES_IN_FIELD_ELECTR)
	{
		CalcChargesInFieldEne();
		el_inter_ene = electr_inter_ene;
	}
	else if(electr_model == COULOMB_ELECTR)
	{
		CalculateMMEnergy();
		el_inter_ene = electr_inter_ene;

	}
	else if(electr_model == CONTINUUM_ELECTR)
	{
		el_inter_ene = CalcContElectrEne(interact_groups);
	}
	return el_inter_ene;
}

double HaInterMolMod::CalcContElectrEne(std::vector<AtomContainer*> inter_groups)
{
	if(interact_groups.size() < 2)
	{
		int ires = SetInteractGroupsFromMolecules();
		if( !ires) 
		{ 
			PrintLog("Error in HaInterMolMod::CalcContElectrEne() \n");
			PrintLog("There are less than two interacting groups in the Intermolecular Module \n");
			return 0.0;
		}
	}
	double inter_energy = 0.0;

	unsigned int ngrp = inter_groups.size();
	if(ngrp < 2)
	{
		ErrorInMod("HaInterMolMod::CalcContElectrEne()",
			"The number of interacting groups is smaller than 2");
		return 0.0;
	}

	MolSet* pmset = GetMolSet();	
	ElectrostMod* el_mod = pmset->GetElectrostMod(true);
	ProtonRedoxMod* p_prot_redox_mod = pmset->GetProtonRedoxMod(true);
	if(el_mod == NULL) return 0.0;
	
	double xmin= 100000.0, xmax=-100000.0;
	double ymin= 100000.0, ymax=-100000.0;
	double zmin= 100000.0, zmax=-100000.0;

	unsigned int igrp;
	HaAtom* aptr;

// Find minimum and maximum coordinates

	for( igrp = 0; igrp < ngrp; igrp++)
	{
		AtomIteratorGen aitr(inter_groups[igrp]);
		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			double x, y, z;

			x = aptr->GetX();
			y = aptr->GetY();
			z = aptr->GetZ();

			if( x < xmin) xmin = x;
			if( x > xmax) xmax = x;
			if( y < ymin) ymin = y;
			if( y > ymax) ymax = y;
			if( z < zmin) zmin = z;
			if( z > zmax) zmax = z;
		}
	}
 
	el_mod->SetBoundaryAtoms(xmin, ymin, zmin, xmax, ymax, zmax);

	AtomIteratorMolSet aitr(phost_mset);
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		aptr->UnSelect();
	}

	AtomGroup active_atoms;

	for( igrp = 0; igrp < ngrp; igrp++)
	{
		AtomIteratorGen aitr(inter_groups[igrp]);
		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			aptr->Select();
			if(this->compute_pk)
			{
				active_atoms.InsertAtom(aptr);
			}
		}
	}

	std::vector<HaResidue*> res_list_tot;
	int nres = pmset->GetNRes();
	res_list_tot.reserve(nres);

	std::map<void*,double> pk_complex;
	std::map<void*,double> pk_part;

	HaResidue* pres;

    double cur_ph = p_prot_redox_mod->GetPH();

	if(this->compute_pk)
	{
		AtomIteratorAtomGroup as_itr(&active_atoms);
		HaResidue* old_pres = NULL;
		for(aptr = as_itr.GetFirstAtom(); aptr; aptr = as_itr.GetNextAtom())
		{
			pres = aptr->GetHostRes();
			if(pres == old_pres) continue;
			old_pres = pres;
			res_list_tot.push_back(pres);
		}
		
		int i; 
		int nrt = res_list_tot.size();

		pmset->SelectAtoms(&active_atoms);
		p_prot_redox_mod->CalcPKaForSelection();
		p_prot_redox_mod->SetChargesForCurrentPH();

		for(i=0; i < nrt; i++)
		{
		    int nalt = p_prot_rdx_mod->GetNumResAltChemStates(res_list_tot[i]);
			if(nalt == 0) continue;
			int j;
			for(j = 0; j < nalt; j++)
			{
				AltChemState* alt_st =p_prot_rdx_mod->GetResAltChemState(res_list_tot[i],j);
				pk_complex[(void*)alt_st] = alt_st->pk;
			}
		}
		pmset->SelectAtoms(&active_atoms);	
	}

	bool bres;

	pmset->save_opt_default.save_selected = TRUE; // only interacting molecules are saved

	double ene_cmplx;
	HaVec_double ene_part(ngrp);
	HaVec_double ene_part_2(ngrp);
	
	bres = el_mod->run(RUN_FOREGROUND);
	ene_cmplx = el_mod->tot_ene;

	if(!bres)
	{
		ErrorInMod("HaIntermolMod::CalcContElectrEne()",
				    "Error reading energy from DELPHI output");
		return false;
	}

	for( igrp = 0; igrp < ngrp; igrp++)
	{		
		pmset->SelectAtoms(inter_groups[igrp]);
			
		bres = el_mod->run(RUN_FOREGROUND);
		ene_part[igrp] = el_mod->tot_ene;
		
		if(!bres)
		{
			ErrorInMod("ElectrostMod::CalcETReorgEne()",
				"Error reading energy of the part from DELPHI output");
			return false;
		}	
	}

	if(this->compute_pk) // Compute pKs for isolated compounds 
	{
		for( igrp = 0; igrp < ngrp; igrp++)
		{
			std::vector<HaResidue*> res_list_part;
			res_list_part.reserve(nres);

			AtomIteratorGen as_itr(inter_groups[igrp]);

			HaResidue* old_res_1 = NULL;
			for(aptr = as_itr.GetFirstAtom(); aptr; aptr = as_itr.GetNextAtom())
			{
				HaResidue* pres = aptr->GetHostRes();
				if(pres == old_res_1) continue;
				old_res_1 = pres;	
				res_list_part.push_back(pres);
			}
			
			int i; 
			int nrp = res_list_part.size(); // the number of residues in the interacting part 
			double cur_ph = p_prot_redox_mod->GetPH();

			pmset->SelectAtoms(inter_groups[igrp]);
			p_prot_redox_mod->CalcPKaForSelection();
			
			for(i=0; i < nrp; i++)
			{
				int nalt = p_prot_rdx_mod->GetNumResAltChemStates(res_list_part[i]);
				if(nalt == 0) continue;
				int j;
				for(j = 0; j < nalt; j++)
				{
					AltChemState* alt_st = p_prot_rdx_mod->GetResAltChemState( res_list_part[i],j);
					pk_part[alt_st] = alt_st->pk;
				}
			}
		}
//		pmset->SetChargesForCurrentPH();
//		pmset->SelectAtoms(&active_atoms);

//		pmset->p_save_opt_default->save_selected = TRUE; // only interacting molecules are saved

//		bres = el_mod->run(RUN_FOREGROUND);
//		ene_cmplx_2 = el_mod->tot_ene;
//				
//		for( igrp = 0; igrp < ngrp; igrp++)
//		{		
//			pmset->SelectAtoms(inter_groups[igrp]);
//			
//			bres = el_mod->run(RUN_FOREGROUND);
//			ene_part_2[igrp] = el_mod->tot_ene;
//		}
//		pmset->SelectAtoms(&active_atoms);

//		double inter_ene_2 = ene_cmplx_2;
//		for( igrp = 0; igrp < ngrp; igrp++)
//		{
//			inter_ene_2 -= ene_part_2[igrp];  
//		}
//		PrintLog(" Electr inter ene for charge distr of isol molecules: %9.3f kT \n",inter_ene_2);
	}		

	phost_mset->save_opt_default.save_selected = FALSE;

	double activ_ch_ene = 0.0; // acivation energy assciated with the changes of the charges of the isolated compounds
		                       // to that of the complex

	if(this->compute_pk) // Compute activation energy associated with pK shifts
	{
		std::map<void*,double>::iterator sitr;
   
		FILE* fapp = fopen("pks_inter_mol.dat","a");

		for( sitr = pk_complex.begin(); sitr != pk_complex.end(); sitr++)
		{
			AltChemState* alt_res_st = (AltChemState*)(*sitr).first;
		    
			if(alt_res_st->active_flag) continue;
//			double ene_unpr_cplx;
//			double ene_pr_cplx;
//			double ene_unpr_part;
//			double ene_pr_part;

			double pk_p = pk_part[(void*)alt_res_st];         // pKs of individual molecules
			double pk_c = pk_complex[(void*)alt_res_st];      // pKs in the complex molecules

			if(fapp != NULL)
			{
				fprintf(fapp,"%8.3f %8.3f ; ",pk_p,pk_c);
			}

			HaResidue* pres= (HaResidue*) alt_res_st->GetHostAtomGroup();
			std::string label = pres->GetRef();
			PrintLog( " Residue : %s    pk individual: %9.3f   pk complex: %9.3f \n", label.c_str(),pk_p,pk_c);

//			ene_unpr_part = 0.0;
                                                     // Changes of the energies of the protonated and
                                                     // unportonated forms of the molecule 
			double ddG_p = log(10.0)*(cur_ph - pk_p);  // - for the isolated compounds (in kT)
			double ddG_c = log(10.0)*(cur_ph - pk_c);  // - for the complex (in kT)

// 1.0/(1.0 + exp(ddG_p)) = exp(-ddG_p)/(1.0 + exp(-ddG_p)) - population of unprotonated form 

			double d_ch = 1.0/(1.0 + exp(ddG_c)) - 1.0/(1.0 + exp(ddG_p)) ;

			double d_ene = + ddG_p * d_ch; // activation energy associated with the changes of the charges
			                                           // of the isolated compounds to the charged of the complex
			activ_ch_ene += d_ene;
		}
		if(fapp != NULL)
		{
			fprintf(fapp," \n");
			fclose(fapp);
		}
		
		PrintLog(" activation energy to change charges of the isolated compounds = %12.6f \n", activ_ch_ene);
	}
	
	
	inter_energy = ene_cmplx;
	if(debug_level >= 5)
		PrintLog( " Energy of the complex is %12.4f kT \n", ene_cmplx);
	for( igrp = 0; igrp < ngrp; igrp++)
	{
		inter_energy -= ene_part[igrp];  
		if(debug_level >= 5)
			PrintLog( " Energy of the molecule # %3d is %12.4f  kT\n ",(igrp+1), ene_part[igrp]);
	}
	
	inter_energy += activ_ch_ene;
	
	if(debug_level >= 5)
		PrintLog( " Interaction energy is %12.4f  kT\n", inter_energy);
	this->electr_inter_ene = inter_energy * KT_300_KCAL;
	
	return inter_energy;
	
}

void intermol_doc_run(InterMolMCSimulator* ptr_im_mc_sim)
{
	if ( ptr_im_mc_sim == NULL) return;
	HaInterMolMod* ptr_im_mod = ptr_im_mc_sim->GetInterMolMod();

	int rex_flag = ptr_im_mc_sim->rex_flag;
	int empirical_flag = ptr_im_mod->empirical_flag;
	if (!rex_flag && !empirical_flag)
	{
		ptr_im_mc_sim->RunMC();
	}
	else if (!rex_flag && empirical_flag)
	{
		//ptr_im_mod -> RunMCEmpirical();
		ptr_im_mc_sim->RunMCQuantSampling();
	}
	else if (rex_flag && empirical_flag)
	{
		ptr_im_mc_sim->RunQuasiREM();
	}
}


int InterMolMCSimulator::InitEnergyFunc()
{
	if(p_im_mod->module_to_init_flag)
	{
		int ires = p_im_mod->Initialize();
		if( !ires) return FALSE;
	}
	return TRUE;
}

double InterMolMCSimulator::ComputeEnergy(harlem::Coord* pcrd)
{
	p_im_mod->SetCoord(p_crd);
	p_im_mod->CalcEffInterEne();
	return p_im_mod->cur_intermol_ene;
}

int InterMolMCSimulator::SetInitPoint(harlem::Coord* p_crd_new)
{
	RigidBodyCoord* p_rb_crd;
	int ng = p_im_mod->interact_groups.size();
	if( p_crd == NULL )
	{
		p_crd = new RigidBodyCoord();
	}
	p_rb_crd = (RigidBodyCoord*)p_crd;
	if(p_crd_new == NULL)
	{
		p_rb_crd->SetNumObj(ng);
		p_rb_crd->SetFromCurrAtomCrd(p_im_mod->interact_groups);
		return TRUE;
	}
	if( p_crd_new->GetClassName() == "RigidBodyCoord" || p_crd_new->GetClassName() == "RigidBodyCoordDiscretized")
	{
		RigidBodyCoord* p_rb_crd_new = (RigidBodyCoord*) p_crd_new;
		int num_obj = p_rb_crd_new->GetNumObj();
		if( num_obj == ng )
		{
			PrintLog("Error in InterMolMCSimulator::SetInitPoint() \n");
			PrintLog("The number objects in pcrd_new %d is not equal \n to the number of interacting groups %d \n", num_obj,ng);
			return FALSE;
		}
		p_crd->SetFrom(p_crd_new);
	}
	if( p_crd_new->GetClassName() == "RigidBodyCoordDiscretized" )
	{
		this->SetDiscretizedMoves();
	}
	return TRUE;
}

int InterMolMCSimulator::SetDiscretizedMoves()
{
	if( IsDiscretizedMoves() ) return TRUE;
	RigidBodyCoordDiscretized* p_d_crd = new RigidBodyCoordDiscretized();
	if( p_crd != NULL ) 
	{
		p_d_crd->SetFrom(p_crd);
		delete p_crd;
		p_crd = p_d_crd;
		return TRUE;
	}
	p_d_crd->SetFromCurrAtomCrd(p_im_mod->interact_groups);
	return TRUE;
}

int InterMolMCSimulator::IsDiscretizedMoves()
{
	if( p_crd == NULL) return FALSE;
	if( p_crd->GetClassName() == "RigidBodyCoordDiscretized") return TRUE;
	return FALSE;
}



void 
HaInterMolMod::SetElectrModel(int new_electr_model_idx)
{
	if ( new_electr_model_idx == CONTINUUM_ELECTR )
		electr_model = CONTINUUM_ELECTR;
	else if ( new_electr_model_idx == COULOMB_ELECTR)
		electr_model = COULOMB_ELECTR;
	else if ( new_electr_model_idx == CHARGES_IN_FIELD_ELECTR)
		electr_model = CHARGES_IN_FIELD_ELECTR;
	else
		electr_model = NO_ELECTR;
}

int HaInterMolMod::InitMolecularFields()
{
	if(interact_groups.size() < 2)
	{
		PrintLog("Error In HaInterMolMod::InitMolecularFields() \n");
		PrintLog(" The Number of interaction Atom Groups is less than two \n");
		return FALSE;
	}
		
	MolSet* pmset     = GetMolSet(); 
	ElectrostMod* el_mod = pmset->GetElectrostMod(1);
	bool bres;

	AtomContainer* pmol1 = interact_groups[0];
	AtomContainer* pmol2 = interact_groups[1];

	if(pmol1 == NULL) return FALSE;

	double xmin,ymin,zmin,xmax,ymax,zmax;

	pmset->GetMinMaxCrd(xmin,ymin,zmin,xmax,ymax,zmax);
	el_mod->SetBoundaryAtoms(xmin, ymin, zmin, xmax, ymax, zmax);

	pmset->save_opt_default.save_selected = TRUE;

	HaAtom* aptr;
	AtomIteratorMolSet aitr_mset(pmset);
	for(aptr = aitr_mset.GetFirstAtom(); aptr; aptr = aitr_mset.GetNextAtom())
	{
		aptr->UnSelect();
	}
	AtomIteratorGen aitr_mol(pmol1);
	for(aptr = aitr_mol.GetFirstAtom(); aptr; aptr = aitr_mol.GetNextAtom())
	{
		aptr->Select();
	}

	bres = el_mod->run(RUN_FOREGROUND);
	if(!bres) return FALSE;
	
	el_pot_field.copy_from(el_mod->el_pot_map);

	pmol1->GetMinMaxCrd(xmin, ymin, zmin, xmax, ymax, zmax);
	
	xmin -= 10.0;
    ymin -= 10.0;
	zmin -= 10.0;
	xmax += 10.0;
	ymax += 10.0;
	zmax += 10.0;


	pmol2->GetMinMaxCrd(xmin, ymin, zmin, xmax, ymax, zmax);

 	xmin -= 20.0;
    ymin -= 20.0;
	zmin -= 20.0;
	xmax += 20.0;
	ymax += 20.0;
	zmax += 20.0;


    vdw_pot_field.clear();
	vdw_pot_field.SetDimensions(181,181,181);
	vdw_pot_field.SetGridCornersCoord(xmin, ymin, zmin, xmax, ymax, zmax);
    vdw_pot_field.FillZeros();

	for(aptr = aitr_mol.GetFirstAtom(); aptr; aptr = aitr_mol.GetNextAtom())
	{
		aptr->Select();
		double rad = el_mod->radprb + aptr->radius;
		double rad2 = rad*rad;

		int mx = (int)(rad/vdw_pot_field.GetXstep()) + 2;
		int my = (int)(rad/vdw_pot_field.GetYstep()) + 2;
		int mz = (int)(rad/vdw_pot_field.GetZstep()) + 2;
		
		int ix,iy,iz;

		vdw_pot_field.GetClosestGridPoint( aptr->GetX(), aptr->GetY(),aptr->GetZ(), ix, iy, iz);

		int mxmin = ((ix - mx) < 0)? 0 : (ix - mx);
		int mymin = ((iy - my) < 0)? 0 : (iy - my);
		int mzmin = ((iz - mz) < 0)? 0 : (iz - mz);

		int mxmax = ((ix + mx) > vdw_pot_field.GetNx()-1) ? vdw_pot_field.GetNx()-1 : (ix + mx);
		int mymax = ((iy + my) > vdw_pot_field.GetNy()-1) ? vdw_pot_field.GetNy()-1 : (iy + my);
		int mzmax = ((iz + mz) > vdw_pot_field.GetNz()-1) ? vdw_pot_field.GetNz()-1 : (iz + mz);

		for( iz= mzmin; iz <= mzmax; iz++)
		{
			for( iy= mymin; iy <= mymax; iy++)
			{
				for( ix= mxmin; ix <= mxmax; ix++)
				{
					float x,y,z;
                    vdw_pot_field.GetXYZ(x,y,z,ix,iy,iz);
					double d2 =  (x - aptr->GetX())*(x - aptr->GetX());
					d2 += (y - aptr->GetY())*(y - aptr->GetY());
					d2 += (z - aptr->GetZ())*(z - aptr->GetZ());
					if(d2 < rad2) vdw_pot_field.SetValue(ix,iy,iz,1.0);
				}
			}
		}
	}
	return TRUE;

}

double
HaInterMolMod::CalcChargesInFieldEne()
{
	if(interact_groups.size() < 2)
	{
		int ires = SetInteractGroupsFromMolecules();
		if( !ires) 
		{ 
			PrintLog("Error in HaInterMolMod::CalcContElectrEne() \n");
			PrintLog("There are less than two interacting groups in the Intermolecular Module \n");
			return 0.0;
		}
	}

	AtomContainer* pmol1 = interact_groups[0];
	AtomContainer* pmol2 = interact_groups[1];

	if(el_pot_field.GetNx() == 0) InitMolecularFields();

	AtomIteratorGen aitr(pmol2);

	HaAtom* aptr;
	
	double inter_ene = 0.0;

	vdw_inter_ene = 0.0;
	electr_inter_ene = 0.0;

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		float x = aptr->GetX(); 
		float y = aptr->GetY();
		float z = aptr->GetZ();

		double vdw_pot = vdw_pot_field.GetInterpolValAtPoint( x, y, z);
		double el_pot  = el_pot_field.GetInterpolValAtPoint ( x, y, z);

		if(vdw_pot > 0.5) 
		{
			vdw_inter_ene = 30000.0*vdw_pot;
			printf("LOOP vdw_inter_ene %6.2f \n", vdw_inter_ene);
		}

        electr_inter_ene += aptr->GetCharge() * el_pot;
	}
    
    electr_inter_ene *= KT_300_KCAL;
	printf("vdw_inter_ene %6.2f, electr_inter_ene %6.2f \n", vdw_inter_ene, electr_inter_ene);
	inter_ene = electr_inter_ene + vdw_inter_ene ;
	return inter_ene;
}





int InterMolMCSimulator::IncrementCrd(Coord* pcrd) 
{
	std::string crd_class_name = pcrd->GetClassName();
	if(crd_class_name == "RigidBodyCoordDiscretized" )
	{
		RigidBodyCoordDiscretized* pcrd_d = (RigidBodyCoordDiscretized*) pcrd;
		int n = pcrd_d->GetNumCrd();
		int i;
		for(i = 0; i < n; i++)
		{
			if( pcrd_d->IsCrdFrozen(i)) continue;
			double rn = p_rand_num_gen->GetRndNum();
			int idx_cur =  pcrd_d->crd_v_int[i];
			if( rn < 0.3333333 )
			{
				if( idx_cur > 0 )
				{
					pcrd_d->crd_v_int[i] = idx_cur - 1;
				}
				else
				{
					pcrd_d->crd_v_int[i] = idx_cur + 1;
				}
			}
			if( rn > 0.6666666 )
			{
				if( idx_cur < (pcrd_d->dim_crd[i] - 1) )
				{
					pcrd_d->crd_v_int[i] = idx_cur + 1;
				}
				else
				{
					pcrd_d->crd_v_int[i] = idx_cur - 1;
				}
			}
		}
		pcrd_d->ConvertDiscrCrdToFloat();
	}
	else
	{
		if( crd_class_name.find("RigidBodyCoord") == std::string::npos )
		{
			throw "Error in InterMolMCSimulator::IncrementCrd() \n crd_class_name != RigidBodyCoord \n";
		}
		RigidBodyCoord* pcrd_rb = (RigidBodyCoord*) pcrd;
		int nmol = pcrd_rb->GetNumObj();
		int im;
		for( im = 0; im < nmol; im++)
		{
			double phi       = pcrd_rb->GetPhi(im);      
		  	double cos_theta = pcrd_rb->GetCosTheta(im);
			double psi       = pcrd_rb->GetPsi(im);
			
			double incr = PI * ang_ratio *( p_rand_num_gen->GetRndNum() - 0.5 );

			Rot3D::IncrEulerAng(phi,cos_theta,psi, incr, 0,0); 
			incr = PI * ang_ratio *( p_rand_num_gen->GetRndNum() - 0.5 ); 
			Rot3D::IncrEulerAng(phi,cos_theta,psi, 0, incr,0); 
			incr = PI * ang_ratio *( p_rand_num_gen->GetRndNum() - 0.5 ); 
			Rot3D::IncrEulerAng(phi,cos_theta,psi, 0, 0, incr); 
			
			pcrd_rb->SetPhi(im,phi);
			pcrd_rb->SetCosTheta(im,cos_theta);
			pcrd_rb->SetPsi(im,psi);

			double x = pcrd_rb->GetTransX(im);
			double y = pcrd_rb->GetTransY(im);
			double z = pcrd_rb->GetTransZ(im);

			x += tr_ratio * (p_rand_num_gen->GetRndNum() - 0.5); 
			y += tr_ratio * (p_rand_num_gen->GetRndNum() - 0.5); 
			z += tr_ratio * (p_rand_num_gen->GetRndNum() - 0.5); 

			pcrd_rb->SetTransX(im,x);
			pcrd_rb->SetTransY(im,y);
			pcrd_rb->SetTransZ(im,z);
		}
	}
	return TRUE;
}

int InterMolMCSimulator::SetCoord(Coord* pcrd)
{
	return p_im_mod->SetCoord(pcrd);
}

HaInterMolMod* InterMolMCSimulator::GetInterMolMod() 
{ 
	return p_im_mod; 
}



int InterMolMCSimulator::RunMCEmpirical() 
{	
	return TRUE;
}


/*

int 
InterMolMCSimulator::RunMCEmpirical() 
{	
	if(p_im_mod->interact_groups.size() < 2)
	{
		PrintLog(" Error in HaInterMolMod::RunMCEmpirical() \n");
		PrintLog(" Less than 2 Interacting Atom Groups are set \n");
		return FALSE;
	}
	
	if(!rex_flag) srand( time(NULL) );

	double step_quat = 1*ang_ratio;
	Random rand_num_gen( rand() );
	RandomGauss rand_num_gen_q( rand() );
	stop_calc_flag = 0;
	int t =0;
	
	//if( interact_groups.size() > 1) 
	//	freeze_first_mol = true;
	
	MolSet* pmset = p_im_mod->GetMolSet();
	HaMolMechMod* pmm_mod = pmset->GetMolMechMod(true);
	HaEmpiricalMod* emp_mod = pmset->GetEmpiricalMod(true);
    

	if(emp_mod->module_to_init_flag)
	{
		emp_mod->Initialize();
	}
	
	double intermol_ene_accepted =0;
	
	if(amber_flag)
	{
		if(pmm_mod->module_to_init_flag)
		{	
			//pmm_mod->build_nb_contact_list_flag = FALSE;
			pmm_mod->Initialize();
			//pmm_mod->max_num_minim_steps=5;
			pmm_mod->SaveAllInpFiles();
		}
	}
	
	int nmol = p_im_mod->interact_groups.size();

	vector<Quaternion> vec_quat_old(nmol);
	vector<Quaternion> vec_quat_cur(nmol);

	Quaternion q_incr;
	Quaternion q_temp;

	double incr;
	double qang, qx, qy,qz, r1, r2, r3, r4, s2;

	clock_t update_time = clock();
	double ene_prev;
	double kt_kcal = MC_temp* BOLTZ * AVOG_NUM/(1000.0*CAL_TO_JOULES); 
	//PrintLog("MC_temp %4.1f \n", MC_temp);  // Remove
	int itr_res = 0;
	int imol;
	
	
	char buf[256];
	
    int iframe = 0;

	int ind = 0;
	
	fstream eff_ene_file;
	
	if(p_im_mod->dont_calc_ene_flag)
	{
		sprintf(buf,"%s%i%s", p_im_mod->p_rex_sim->MC_energy_file_replica_basename.c_str(), p_im_mod->p_rex_sim->n_playback_replica,".dat");
		eff_ene_file.open(buf, ios::in);
	}
	else
	{
		eff_ene_file.open("eff_ene_file.dat", ios::out);
	}

	fstream traj_file;
	fstream all_points_file;

	sprintf(buf,"%s%i%s", p_im_mod->p_rex_sim->MC_traj_file_replica_basename.c_str(), p_im_mod->p_rex_sim->ireplica,"REG.dat");
	traj_file.open(buf, ios::out|ios::app  );

	if(output_rejected_points)
	{
		all_points_file.open("all_mc_pts.dat",ios::out);
	}
	
	if( traj_file.fail())
	{
		ErrorInMod("HaInterMolMod::RunMCEmpirical()",
			"Can't open file to write MC trajectory ");
		return FALSE;
	}
	if( freeze_first_mol)
	{
		AtomContainer* pMol1 = p_im_mod->interact_groups[0];
		Vec3D trans1;
		
		pMol1->GetQuaternionTrans(vec_quat_cur[0], trans1);
		vec_quat_cur[0].GetQuaternion(qang, qx, qy,qz);
		
		sprintf(buf,"%6d  %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f",
			-1, trans1[0], trans1[1], trans1[2], qang, qx, qy,qz);
		
		traj_file << buf << endl;
	}

    HaMat_double cur_trans(3,nmol);
	Vec3D trans_v;
	
    HaMat_double trans_old(3,nmol);	

	double average_eff_ene = 0.0; 

	for( imol= 0; imol < nmol; imol++)
	{
		AtomContainer* pMol = p_im_mod->interact_groups[imol];
		pMol->GetQuaternionTrans(vec_quat_cur[imol], trans_v);
		
		cur_trans(1,imol+1) = trans_v[0];
		cur_trans(2,imol+1) = trans_v[1];
		cur_trans(3,imol+1) = trans_v[2];
	}

	int np_accept = 0; 
	int np_reject = 0;
	
	VecPtr pt_mol;
	pt_mol.resize(nmol);
	int j;
	for(j =0; j < nmol ; j++ )
	{
		AtomContainer* pMol = p_im_mod->interact_groups[j];
		pt_mol[j] = pMol;
	}

	int imode= 1;
	
	// Previous energy
	if( rex_flag){
		ene_prev = p_im_mod->p_rex_sim->energy_arr[ p_im_mod->p_rex_sim->exchange_arr[p_im_mod->p_rex_sim->ireplica] ];
	}
	else ene_prev = 10e9;
	for( int istep = 0; istep < num_mc_steps ; istep++)       //   MC steps
	{
		bool step_move_flag= false;

		if (istep ==0 || (istep+1)%100 ==0)	PrintLog("---MC step %d \n", istep);
		if(!rex_flag && istep%100 == 0) 
		{
			emp_mod->Neighborhood();
		}
		
		if( stop_calc_flag != 0 )
		{
			stop_calc_flag = 0;
			break;
		}
		
		int moved_flag = FALSE;

		for( imol = 0; imol < nmol; imol++)
		{
			if( freeze_first_mol && (imol == 0) )
					continue;

			AtomContainer* pMol = p_im_mod->interact_groups[imol];

					
			vec_quat_old[imol] = vec_quat_cur[imol];
					
			trans_old(1,imol+1) = cur_trans(1,imol+1);
			trans_old(2,imol+1) = cur_trans(2,imol+1);
			trans_old(3,imol+1) = cur_trans(3,imol+1);
			q_incr = vec_quat_cur[imol];

			for (;;)
			{
				r1 = rand_num_gen_q() * step_quat;
				r2 = rand_num_gen_q() * step_quat;
				r3 = rand_num_gen_q() * step_quat;
				s2=r1*r1+r2*r2+r3*r3;
				if (s2<1.0) break;
			}
			r4=sqrt(1-s2);
			q_temp.SetQuaternion(r4,r1,r2,r3);
			q_incr.operator *=(q_temp);

			vec_quat_cur[imol] = q_incr;
	 
			cur_trans(1,imol+1) += tr_ratio * (rand_num_gen() - 0.5); rand_num_gen();
			cur_trans(2,imol+1) += tr_ratio * (rand_num_gen() - 0.5); rand_num_gen();
			cur_trans(3,imol+1) += tr_ratio * (rand_num_gen() - 0.5); rand_num_gen();
				
								
			int ic = 0;
			if(!p_im_mod->discr_moves_flag || ic == 0)
			{
				Vec3D trans_l;
				trans_l[0] = cur_trans(1,imol+1);
				trans_l[1] = cur_trans(2,imol+1);
				trans_l[2] = cur_trans(3,imol+1);
					
				pMol->SetQuaternionTrans(vec_quat_cur[imol], trans_l);	
				moved_flag = TRUE;

				//SteepestDescentMinimizer(50*nmol);	
				if (!equil_conf_vol_vdw)
					p_im_mod->cur_intermol_ene = emp_mod->ScoreEnergy();
				else
					p_im_mod->cur_intermol_ene = emp_mod->GeometryScoreEnergy();
				if((istep+1)%100 == 0) PrintLog("Score Energy %10.3f imol %d\n", p_im_mod->cur_intermol_ene, imol);
			}
				
			bool is_accept = false;	

		    if( p_im_mod->cur_intermol_ene < ene_prev ) 
			{
				is_accept = true;
//				PrintLog("cur_intermol_ene %8.4f\n" , cur_intermol_ene);
			}
			else
			{
				double rand_dbl = rand_num_gen(); // rand()/RAND_MAX;
				double dtest = exp( (ene_prev - p_im_mod->cur_intermol_ene)/kt_kcal );

				if( dtest > rand_dbl )
				{
					is_accept = true;
				}
				else
				{
				//	PrintLog(" Point is rejected \n");
					np_reject++;
							
					ind =0;						
					if(moved_flag)
					{
						Vec3D trans_o;
						trans_o[0] = trans_old(1,imol+1);
						trans_o[1] = trans_old(2,imol+1);
						trans_o[2] = trans_old(3,imol+1);
							
						pMol->SetQuaternionTrans(vec_quat_old[imol], trans_o);
					}

					vec_quat_cur[imol] = vec_quat_old[imol];
					cur_trans(1,imol+1) = trans_old(1,imol+1);
					cur_trans(2,imol+1) = trans_old(2,imol+1);
					cur_trans(3,imol+1) = trans_old(3,imol+1);
							
				}
			}

			if( is_accept )
			{
				pmset->info_str.clear();
					
				ene_prev = p_im_mod->cur_intermol_ene;
				intermol_ene_accepted = p_im_mod->cur_intermol_ene;
				np_accept++;
					
				average_eff_ene+= p_im_mod->cur_intermol_ene; 
				eff_ene_file << np_accept << "  " << p_im_mod->cur_intermol_ene ;
				eff_ene_file << endl;

				if(output_rejected_points)
				{
					for (int ii = 0; ii < nmol; ii++)
					{
						vec_quat_cur[ii].GetQuaternion(qang, qx, qy, qz);
						sprintf(buf,"%6d  %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f",
						np_accept, cur_trans(1,ii+1), cur_trans(2,ii+1), cur_trans(3,ii+1), qang, qx, qy, qz);
						all_points_file << buf << endl; //added by Jose
					}
				}
			}
			if (is_accept) step_move_flag= true;
		} // End of for(imol)

		if( clock() > update_time || istep == (num_mc_steps - 1))
		{
			if(!moved_flag)
			{
				for( imol = 0; imol < nmol; imol++)
				{
					if( freeze_first_mol && (imol == 0) )
						continue;
					AtomContainer* pMol = p_im_mod->interact_groups[imol];
					Vec3D trans_l;
					trans_l[0] = cur_trans(1,imol+1);
					trans_l[1] = cur_trans(2,imol+1);
					trans_l[2] = cur_trans(3,imol+1);
					pMol->SetQuaternionTrans(vec_quat_cur[imol], trans_l);
					moved_flag = TRUE;
				}
			}				
			if (istep%10 ==0) pmset->AnnounceGeomChange();
			update_time = clock() + CLOCKS_PER_SEC*1 ;
			//PrintLog("update_time= %d", update_time);
		}
		if (intermol_ene_accepted != 0.0)
		{
			pmset->info_str.clear();
			sprintf(buf,"points accepted %d Total effective energy is %14.6f kcal/mol Replica %d",np_accept,intermol_ene_accepted, p_im_mod->p_rex_sim->ireplica);
			pmset->info_str.push_back(buf);
		}
		else
		{
			pmset->info_str.clear();
			sprintf(buf,"No points accepted, Replica %d", p_im_mod->p_rex_sim->ireplica);
			pmset->info_str.push_back(buf);
		}
	} // end MC step
//	sprintf(buf," Total effective energy is %14.6f kcal/mol ",intermol_ene_accepted);
//	PrintMessage(buf);
//	pmset->info_str.push_back(buf);
	
	t=0; 
	int jmol;

	for (jmol =0 ; jmol < nmol; jmol++)
	{
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, cur_trans(1,jmol+1));
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, cur_trans(2,jmol+1));
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, cur_trans(3,jmol+1));
		vec_quat_cur[jmol].GetQuaternion(qang, qx, qy,qz);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qang);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qx);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qy);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qz);
	}

	int n_visit = pt_visit.size();
	//PrintLog("Number of points visited %d \n", n_visit); 
	
	if( np_accept > 0) 
	{
		double acc_ratio = (double)np_accept/( (double) np_accept+ (double) np_reject);
		//PrintLog(" point acceptance ratio: %7.3f \n", acc_ratio );  
		//PrintLog(" average effective energy: %12.6e \n", average_eff_ene/np_accept );
		p_acc_ratio[ireplica] = acc_ratio ;
		if(average_eff_ene != 0) energy_arr[ireplica] = average_eff_ene/np_accept;
	}

	return True;
	
}

*/

int InterMolMCSimulator::RunMCQuantSampling()   
{	
	return TRUE;
}


/*

int InterMolMCSimulator::RunMCQuantSampling()   
   // Added by jose
{	// This is a copy of RunMCEmpirical()

	if(module_to_init_flag)
	{
		Initialize();
	}
	
	if(!rex_flag) srand( time(NULL) );

	double step_quat = 1*ang_ratio;
	Random rand_num_gen( rand() );
	RandomGauss rand_num_gen_q( rand() );
	stop_calc_flag = 0;
	int t =0;
	
	bool freeze_first_mol = false;
	//if( interact_groups.size() > 1) 
	//	freeze_first_mol = true;
	
	MolSet* pmset = GetMolSet();
//	HaMolMechMod* pmm_mod = pmset->GetMolMechMod(true); jose
	//HaEmpiricalMod* emp_mod = pmset->GetEmpiricalMod(true);
	HaMolMembraneMod* mem_mod = pmset->GetMolMembraneMod(true); //jose
	if (mem_mod->module_to_init_flag) //jose
	{
		mem_mod->Initialize();
//		mem_mod->SetCoarseGrainedDFireCoreParams();
		mem_mod->display_results_flag = FALSE;
	}
	
	double intermol_ene_accepted =0;
		
	int nmol = interact_groups.size();

	vector<Quaternion> vec_quat_old(nmol);
	vector<Quaternion> vec_quat_cur(nmol);

	Quaternion q_incr;
	Quaternion q_temp;

	double incr;
	double qang, qx, qy,qz, r1, r2, r3, r4;

	typedef std::chrono::system_clock Time;
	auto update_time = Time::now();

	double ene_prev;
	double kt_kcal = MC_temp* BOLTZ * AVOG_NUM/(1000.0*CAL_TO_JOULES); 
	//PrintLog("MC_temp %4.1f \n", MC_temp);  // Remove
	int itr_res = 0;
	int imol;
	
	
	char buf[256], buf1[256],ene_str[256];;
	
    int iframe = 0;

	int ind = 0;
	
	fstream eff_ene_file;
	
	if(dont_calc_ene_flag)
	{
		sprintf(buf,"%s%i%s", MC_energy_file_replica_basename.c_str(), n_playback_replica,".dat");
		eff_ene_file.open(buf, ios::in);
	}
	else
	{
		eff_ene_file.open("eff_ene_file.dat", ios::out);
	}
	
	fstream traj_file;
	fstream all_points_file;
	Vec3D trans1;
	double qang_r, qx_r, qy_r, qz_r, state_r;
	if (play_back_flag)
	{
		sprintf(buf,"%s%i%s", MC_traj_file_replica_basename.c_str(), ireplica,"REG.dat");
	//	sprintf(buf,"all_mc_pts.dat");
	// trajectory file
		traj_file.open(buf,ios::in  );
		if( traj_file.fail())
		{
			ErrorInMod("HaInterMolMod::RunMCQuantSampling()",
				       "Can't open file with MC trajectory ");
			return FALSE;
		}

// Read coordinates of the first molecule
		
		istrstream is(buf);
		int junk;

		if(npt_begin > 0)
		{
			for( int i = 0; i < npt_begin; i++)
				traj_file.getline(buf,255);
			if( dont_calc_ene_flag )
			{
				for(int i = 0; i < npt_begin; i++)
					eff_ene_file.getline(buf,255);
			}
		}
		if(traj_file.fail())
		{
			stop_calc_flag = TRUE;
			return FALSE;	
		}	
	}
	else
	{
		sprintf(buf,"%s%i%s", MC_traj_file_replica_basename.c_str(), ireplica,"REG.dat");
		traj_file.open(buf, ios::out|ios::app  );

		if(output_rejected_points)
		{
			all_points_file.open("all_mc_pts.dat",ios::out);
		}
		
		if( traj_file.fail())
		{
			ErrorInMod("HaInterMolMod::RunMCEmpirical()",
				"Can't open file to write MC trajectory ");
			return FALSE;
		}
		if( freeze_first_mol)
		{
			AtomContainer* pMol1 = interact_groups[0];
			Vec3D trans1;
			
			pMol1->GetQuaternionTrans(vec_quat_cur[0], trans1);
			vec_quat_cur[0].GetQuaternion(qang, qx, qy,qz);
			
			sprintf(buf,"%6d  %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f",
				-1, trans1[0], trans1[1], trans1[2], qang, qx, qy,qz);
			
			traj_file << buf << endl;
		}
	}
	
    HaMat_double cur_trans(3,nmol);
	Vec3D trans_v;
	
    HaMat_double trans_old(3,nmol);	

	for( imol= 0; imol < nmol; imol++)
	{
		AtomContainer* pMol = interact_groups[imol];
		pMol->GetQuaternionTrans(vec_quat_cur[imol], trans_v);
		
		cur_trans(1,imol+1) = trans_v[0];
		cur_trans(2,imol+1) = trans_v[1];
		cur_trans(3,imol+1) = trans_v[2];
	}

	int np_accept = 0; 
	int np_reject = 0;
	
	double average_eff_ene = 0.0;
	
	VecPtr pt_mol;
	pt_mol.resize(nmol);
	int j;
	for(j =0; j < nmol ; j++ )
	{
		AtomContainer* pMol = interact_groups[j];
		pt_mol[j] = pMol;
	}

	int imode= 1;
	
	// Previous energy
	if( rex_flag){
		ene_prev = energy_arr[ exchange_arr[ireplica] ];
	}
	else ene_prev = 10e9;
//	PrintLog("First ene_prev %2.3f \n",ene_prev);
	int count_step=0;
	bool flag_step = false;
	int count_step_i = 0;
	for( int istep = 0; istep < num_mc_steps ; istep++)       //   MC steps
	{
		bool step_move_flag= false;

		if (istep ==0 || (istep+1)%100 ==0)	PrintLog("---MC step %d \n", istep);
		if(!rex_flag && istep%100 == 0) 
		{
			//emp_mod ->Neighborhood(); remove for not cylinder
		}
		
		if( stop_calc_flag != 0 )
		{
			stop_calc_flag = 0;
			break;
		}
		
		int moved_flag = FALSE;

		for( imol = 0; imol < nmol; imol++)
		{
//sprintf(buf,"%2.3f %2.3f", ang_ratio, tr_ratio);
//check1 << buf << endl;
			if( freeze_first_mol && (imol == 0) )
					continue;

				AtomContainer* pMol = interact_groups[imol];
				if(!play_back_flag )
				{
					
					vec_quat_old[imol] = vec_quat_cur[imol];
					
					trans_old(1,imol+1) = cur_trans(1,imol+1);
					trans_old(2,imol+1) = cur_trans(2,imol+1);
					trans_old(3,imol+1) = cur_trans(3,imol+1);
					q_incr = vec_quat_cur[imol];

					if ((np_accept) % 10 ==0 && !flag_step)
					{
						flag_step = true;
						count_step_i = 0;
					}
					if (flag_step)
					{
						incr = PI*ang_ratio *( 2.0* rand_num_gen() - 1.0 );rand_num_gen();
						incr *=100;
						q_temp.QuaternionFromAxis(incr, 0, 1.0, 0);
						q_incr.operator *=(q_temp);
						count_step_i +=1;
						if ((count_step_i %10) ==0)
						{	count_step_i = 0;
							flag_step = false;
						}
					}
					else 
					{
					incr = PI*ang_ratio *( 2.0* rand_num_gen() - 1.0 );rand_num_gen();
					//incr = PI*ang_ratio *( rand_num_gen_q() );rand_num_gen_q();
					q_temp.QuaternionFromAxis(incr, 1.0, 0, 0);
					q_incr.operator *=(q_temp);
					incr = PI*ang_ratio *( 2.0* rand_num_gen() - 1.0 );rand_num_gen();
//					incr = PI*ang_ratio *( rand_num_gen_q() );rand_num_gen_q();
					q_temp.QuaternionFromAxis(incr, 0, 1.0, 0);
					q_incr.operator *=(q_temp); 
					incr = PI*ang_ratio *( 2.0* rand_num_gen() - 1.0 );rand_num_gen();
//					incr = PI*ang_ratio *( rand_num_gen_q() );rand_num_gen_q();
					q_temp.QuaternionFromAxis(incr,  0, 0, 1.0);
					q_incr.operator *=(q_temp);
					}


					vec_quat_cur[imol] = q_incr;
	 
						cur_trans(1,imol+1) += tr_ratio * (rand_num_gen() - 0.5); rand_num_gen();
						cur_trans(2,imol+1) += tr_ratio * (rand_num_gen() - 0.5); rand_num_gen();
						cur_trans(3,imol+1) += tr_ratio * (rand_num_gen() - 0.5); rand_num_gen();
				}
// playback
				if(play_back_flag)
				{
				// Read MC Trajectory point:
				//
					traj_file.getline(buf,255);
					if( dont_calc_ene_flag ) eff_ene_file.getline(ene_str,255);
					
					if(npt_step > 1)
					{
						for(int  i = 0; i < (npt_step-1); i++)
						{
							traj_file.getline(buf1,255);
							if( dont_calc_ene_flag )
							{
								eff_ene_file.getline(buf1,255);
							}
						}
					}
					if(traj_file.fail())
					{
						stop_calc_flag = TRUE;
						break;
					}
					istrstream is(buf);
					
					int junk;

					AtomContainer* pMol1 = interact_groups[imol];
					is >> junk;
					is >> trans1[0] >> trans1[1] >> trans1[2];
					is >> qang_r >> qx_r >> qy_r >> qz_r >> state_r;

					cur_trans(1,imol+1) = trans1[0];
					cur_trans(2,imol+1) = trans1[1];
					cur_trans(3,imol+1) = trans1[2];

					q_temp.SetQuaternion(qang_r, qx_r, qy_r, qz_r);

					vec_quat_cur[imol] = q_temp;

					if(is.fail())
					{
						stop_calc_flag = TRUE;
						break;
					}
				}
//
								
				if(!dont_calc_ene_flag)
				{
					int ic = 0;
					if(!discr_moves_flag || ic == 0)
					{
						Vec3D trans_l;
						trans_l[0] = cur_trans(1,imol+1);
						trans_l[1] = cur_trans(2,imol+1);
						trans_l[2] = cur_trans(3,imol+1);
					
						pMol->SetQuaternionTrans(vec_quat_cur[imol], trans_l);	
						moved_flag = TRUE;
						count_step++;
						//if (flag_step)
						//	mem_mod->build_nb_coarsegrained_contact_list = false;
						//else
						//	mem_mod->build_nb_coarsegrained_contact_list = true;
						
//						SteepestDescentMinimizer(50*nmol);	
						if (!equil_conf_vol_vdw)
						{
						//	cur_intermol_ene = emp_mod->ScoreEnergy();
						//	emp_mod->QuantSampling();
						//	cur_intermol_ene = emp_mod->LennardJonesEnergy();
//							cur_intermol_ene = emp_mod->ScoreEnergy();
							mem_mod->ScoreEnergy();
							cur_intermol_ene = mem_mod->tot_energy;
						}
						else
							//cur_intermol_ene = emp_mod->GeometryScoreEnergy();
						{
							mem_mod->pairwiseDfire_flag = false;
							mem_mod->ScoreEnergy();
							cur_intermol_ene = mem_mod->tot_energy;
						}
//						cur_intermol_ene = emp_mod->MinEnergy();
//						cur_intermol_ene = emp_mod->HarmonicEnergy();
//						cur_intermol_ene = emp_mod->ToyEnergy();
//							sprintf(buf," %2.3f ", cur_intermol_ene);
//				check4 << buf << endl;
						if((istep+1)%100 == 0) PrintLog("Score Energy %10.3f imol %d\n", cur_intermol_ene, imol);
					}
					
				}
				
				bool is_accept = false;	
				
				if(!play_back_flag)
				{
					if( cur_intermol_ene < ene_prev ) 
					{
						is_accept = true;
//						PrintLog("cur_intermol_ene %8.4f\n" , cur_intermol_ene);
					}
					else
					{
						double rand_dbl = rand_num_gen(); // rand()/RAND_MAX;
						double dtest = exp( (ene_prev - cur_intermol_ene)/kt_kcal );
						//double dtest = exp( (ene_prev - cur_intermol_ene)/MC_temp ) ;
//						PrintLog("ene_prev %8.4f, cur_intermol_ene %8.4f, rand_dbl %8.4f, DTEST= %6.4f \n" , ene_prev , cur_intermol_ene , rand_dbl ,dtest);  // Remove
						if( dtest > rand_dbl )
						{
							is_accept = true;
						}
						else
						{
						//	PrintLog(" Point is rejected \n");
							np_reject++;
							
							ind =0;						
							if(moved_flag)
							{
								Vec3D trans_o;
								trans_o[0] = trans_old(1,imol+1);
								trans_o[1] = trans_old(2,imol+1);
								trans_o[2] = trans_old(3,imol+1);
							
								pMol->SetQuaternionTrans(vec_quat_old[imol], trans_o);
							}

							
							vec_quat_cur[imol] = vec_quat_old[imol];
							cur_trans(1,imol+1) = trans_old(1,imol+1);
							cur_trans(2,imol+1) = trans_old(2,imol+1);
							cur_trans(3,imol+1) = trans_old(3,imol+1);
							
						}
					}
				}

				if (is_accept)
				{
					ene_prev = cur_intermol_ene;
					average_eff_ene += cur_intermol_ene;
				}
				//else
				//{
				//	average_eff_ene += ene_prev;
				//}


				if( (is_accept || play_back_flag))
				{
					pmset->info_str.clear();
					
					ene_prev = cur_intermol_ene;
					intermol_ene_accepted = cur_intermol_ene;
					np_accept++;

					
					if( !dont_calc_ene_flag)
					{
						//average_eff_ene+= cur_intermol_ene; 
						eff_ene_file << count_step << "  " << cur_intermol_ene ;
						eff_ene_file << endl;
					}
			}
		} // End of for(imol)

		if(delay_time > 0)
		{
			std::this_thread::sleep_for(std::chrono::milliseconds(delay_time*1000));
		}

		//if( Time::no > update_time || istep == (num_mc_steps - 1))
		if(1)
		{
			if(!moved_flag)
			{
				for( imol = 0; imol < nmol; imol++)
				{
					if( freeze_first_mol && (imol == 0) )
						continue;
					AtomContainer* pMol = interact_groups[imol];
					Vec3D trans_l;
					trans_l[0] = cur_trans(1,imol+1);
					trans_l[1] = cur_trans(2,imol+1);
					trans_l[2] = cur_trans(3,imol+1);
					pMol->SetQuaternionTrans(vec_quat_cur[imol], trans_l);
					moved_flag = TRUE;
				}
			}
			if(dont_calc_ene_flag)
			{
				pmset->RefreshAllViews(RFApply | RFRefresh);
				PrintLog("RefreshAllViews() \n");
				if(save_image_seq_gif || save_image_seq_pict)
				{
					iframe++;
					HaMolView* pview = pmset->GetActiveMolView();
					if(save_image_seq_gif)
					{
						sprintf(buf,"traj%6d.gif",(100000 + iframe));
						pview->WriteGIFFile(buf);
					}
					if(save_image_seq_pict)
					{
						sprintf(buf,"traj%6d.pic",(100000 + iframe));
						pview->WritePICTFile(buf);
					}
				}
			}
			//if (istep%10 ==0) 
			else 
			{
				pmset->AnnounceGeomChange();
				update_time = Time::now();
				update_time += std::chrono::milliseconds( 1000 ); 
			}
		}
		if (intermol_ene_accepted != 0.0)
		{
			pmset->info_str.clear();
			sprintf(buf,"points accepted %d Total effective energy is %14.6f kcal/mol Replica %d",np_accept,intermol_ene_accepted, ireplica);
			pmset->info_str.push_back(buf);
		}
		else
		{
			pmset->info_str.clear();
			sprintf(buf,"No points accepted, Replica %d", ireplica);
			pmset->info_str.push_back(buf);
		}
	} // end MC step
//	sprintf(buf," Total effective energy is %14.6f kcal/mol ",intermol_ene_accepted);
//	PrintMessage(buf);
//	pmset->info_str.push_back(buf);
	
	t=0; 
	int jmol;

	for (jmol =0 ; jmol < nmol; jmol++)
	{
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, cur_trans(1,jmol+1));
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, cur_trans(2,jmol+1));
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, cur_trans(3,jmol+1));
		vec_quat_cur[jmol].GetQuaternion(qang, qx, qy,qz);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qang);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qx);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qy);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qz);
	}

	int n_visit = pt_visit.size();
	//PrintLog("Number of points visited %d \n", n_visit); 
	
	if( np_accept > 0) 
	{
		double acc_ratio = (double)np_accept/( (double) np_accept+ (double) np_reject);
		//PrintLog(" point acceptance ratio: %7.3f \n", acc_ratio );  
		//PrintLog(" average effective energy: %12.6e \n", average_eff_ene/np_accept );
		p_acc_ratio[ireplica] = acc_ratio ;
		if(average_eff_ene != 0) energy_arr[ireplica] = average_eff_ene/np_accept; //jose
		//if (average_eff_ene!=0.0) energy_arr[ireplica] = average_eff_ene/double(count_step);
		//energy_arr[ireplica] = ene_prev;
	}

	return True;	
}

*/

int InterMolMCSimulator::RunMCEmpiricalXY() 
{
	return TRUE;
}

/*

int InterMolMCSimulator::RunMCEmpiricalXY() 
{	

	if(module_to_init_flag)
	{
		Initialize();
	}
	
	if(!rex_flag) srand( time(NULL) );

	Random rand_num_gen( rand() );
	stop_calc_flag = 0;
	int t =0;
	
	bool freeze_first_mol = false;
	//		if( interact_groups.size() > 1) 
	//freeze_first_mol = true;
	
	MolSet* pmset = GetMolSet();
	HaMolMechMod* pmm_mod = pmset->GetMolMechMod(true);
	HaEmpiricalMod* emp_mod = pmset->GetEmpiricalMod(true);
	if(emp_mod -> module_to_init_flag)
	{
		emp_mod -> Initialize();
	}
	
	double intermol_ene_accepted =0;
	
	if(amber_flag)
	{
		if(pmm_mod->module_to_init_flag)
		{	
			//pmm_mod->build_nb_contact_list_flag = FALSE;
			pmm_mod->Initialize();
			//pmm_mod->max_num_minim_steps=5;
			pmm_mod->SaveAllInpFiles();
		}
	}
	
	int nmol = interact_groups.size();


	vector<Quaternion> vec_quat_old(nmol);
	vector<Quaternion> vec_quat_cur(nmol);

	Quaternion q_incr;
	Quaternion q_temp;
	double incr;
	double qang, qx, qy,qz;

	clock_t update_time = clock();
	double ene_prev;
	double kt_kcal = MC_temp* BOLTZ * AVOG_NUM/(1000.0*CAL_TO_JOULES); 
	//PrintLog("MC_temp %4.1f \n", MC_temp);  // Remove
	int itr_res = 0;
	int imol;
	
	
	char buf[256];
	
    int iframe = 0;

	int ind = 0;
	
	fstream eff_ene_file;
	
	if(dont_calc_ene_flag)
	{
		sprintf(buf,"%s%i%s", MC_energy_file_replica_basename.c_str(), n_playback_replica,".dat");
		eff_ene_file.open(buf, ios::in);
	}
	else
	{
		eff_ene_file.open("eff_ene_file.dat", ios::out);
	}
	
	fstream traj_file;
	fstream all_points_file;
	

		sprintf(buf,"%s%i%s", MC_traj_file_replica_basename.c_str(), ireplica,"REG.dat");
		traj_file.open(buf, ios::out|ios::app  );

		if(output_rejected_points)
		{
			all_points_file.open("all_mc_pts.dat",ios::out);
		}
		
		if( traj_file.fail())
		{
			ErrorInMod("HaInterMolMod::RunMCEmpirical()",
				"Can't open file to write MC trajectory ");
			return FALSE;
		}
		if( freeze_first_mol)
		{
			AtomContainer* pMol1 = interact_groups[0];
			Vec3D trans1;
			
			pMol1->GetQuaternionTrans(vec_quat_cur[0], trans1);
			vec_quat_cur[0].GetQuaternion(qang, qx, qy,qz);
			
			sprintf(buf,"%6d  %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f",
				-1, trans1[0], trans1[1], trans1[2], qang, qx, qy,qz);
			
			traj_file << buf << endl;
		}

		
	
    HaMat_double cur_trans(3,nmol);
	Vec3D trans_v;
	
    HaMat_double trans_old(3,nmol);	



	for( imol= 0; imol < nmol; imol++)
	{
		AtomContainer* pMol = interact_groups[imol];
		pMol->GetQuaternionTrans(vec_quat_cur[imol], trans_v);
		
		cur_trans(1,imol+1) = trans_v[0];
		cur_trans(2,imol+1) = trans_v[1];
		cur_trans(3,imol+1) = trans_v[2];
	}

	int np_accept = 0; 
	int np_reject = 0;
	
	double average_eff_ene = 0.0;
	
	VecPtr pt_mol;
	pt_mol.resize(nmol);
	int j;
	for(j =0; j < nmol ; j++ )
	{
		AtomContainer* pMol = interact_groups[j];
		pt_mol[j] = pMol;
	}

	int imode= 1;
	
	// Previous energy
	if( rex_flag){
		ene_prev = energy_arr[ exchange_arr[ireplica] ];
	}
	else ene_prev = 10e9;
	
	//freeze_first_mol = true; //REMOVE
	for( int istep = 0; istep < num_mc_steps ; istep++)       //   MC steps
	{
		bool step_move_flag= false;

		if ((istep+1)%100 ==0)	PrintLog("---MC step %d \n", istep);
		if(!rex_flag && istep%100 == 0) 
		{
			emp_mod ->Neighborhood();
		}
		
		if( stop_calc_flag != 0 )
		{
			stop_calc_flag = 0;
			break;
		}
		
		int moved_flag = FALSE;
		
		for( imol = 0; imol < nmol; imol++)
		{
			if( freeze_first_mol && (imol == 0) )
				continue;

				AtomContainer* pMol = interact_groups[imol];

					
				vec_quat_old[imol] = vec_quat_cur[imol];					
				trans_old(1,imol+1) = cur_trans(1,imol+1);
				trans_old(2,imol+1) = cur_trans(2,imol+1);
				trans_old(3,imol+1) = cur_trans(3,imol+1);
//					q_incr = vec_quat_cur[imol];
//					incr = PI*ang_ratio *( 2*rand_num_gen() - 1.0 );rand_num_gen();
//					q_temp.QuaternionFromAxis(incr, 1.0, 0, 0);
//					q_incr.operator *=(q_temp);
//					incr = PI*ang_ratio *( 2*rand_num_gen() - 1.0 );rand_num_gen();
//					q_temp.QuaternionFromAxis(incr, 0, 1.0, 0);
//					q_incr.operator *=(q_temp);
//					incr = PI*ang_ratio *( 2*rand_num_gen() - 1.0 );rand_num_gen();
//					q_temp.QuaternionFromAxis(incr,  0, 0, 1.0);
//					q_incr.operator *=(q_temp);

//					vec_quat_cur[imol] = q_incr;



	 
				cur_trans(1,imol+1) += tr_ratio * (rand_num_gen() - 0.5); rand_num_gen();
				cur_trans(2,imol+1) += tr_ratio * (rand_num_gen() - 0.5); rand_num_gen();
//				cur_trans(3,imol+1) += tr_ratio * (rand_num_gen() - 0.5); rand_num_gen();

				
								
				int ic = 0;
				if(!discr_moves_flag || ic == 0)
				{
					Vec3D trans_l;
					trans_l[0] = 0;
					trans_l[1] = 0;
					trans_l[2] = 0;
					trans_l[0] = cur_trans(1,imol+1);
					trans_l[1] = cur_trans(2,imol+1);
					trans_l[2] = cur_trans(3,imol+1);
	
					pMol->SetQuaternionTrans(vec_quat_cur[imol], trans_l);	
					moved_flag = TRUE;

					cur_intermol_ene = emp_mod->GeometryScoreEnergy();
					if((istep+1)%100 == 0) PrintLog("Score Energy %10.3f \n", cur_intermol_ene);
				}
					
				
				bool is_accept = false;	
				
				if( cur_intermol_ene < ene_prev ) 
				{
					is_accept = true;
				}
				else
				{
					double rand_dbl = rand_num_gen(); // rand()/RAND_MAX;
					double dtest = exp( (ene_prev - cur_intermol_ene)/kt_kcal );
					if( dtest > rand_dbl )
					{
						is_accept = true;
					}
					else
					{
					//	PrintLog(" Point is rejected \n");
						np_reject++;
							
						ind =0;						
						if(moved_flag)
						{
							Vec3D trans_o;
							trans_o[0] = 0;
							trans_o[1] = 0;
							trans_o[2] = 0;
							trans_o[0] = trans_old(1,imol+1);
							trans_o[1] = trans_old(2,imol+1);
							trans_o[2] = trans_old(3,imol+1);
							
							pMol->SetQuaternionTrans(vec_quat_old[imol], trans_o);
						}

						vec_quat_cur[imol] = vec_quat_old[imol];
						cur_trans(1,imol+1) = trans_old(1,imol+1);
						cur_trans(2,imol+1) = trans_old(2,imol+1);
						cur_trans(3,imol+1) = trans_old(3,imol+1);	
					}
				}


				if( is_accept)
				{
					pmset->info_str.clear();
					
					ene_prev = cur_intermol_ene;
					intermol_ene_accepted = cur_intermol_ene;
					np_accept++;

					
					if( !dont_calc_ene_flag)
					{
						average_eff_ene+= cur_intermol_ene; 
						eff_ene_file << np_accept << "  " << cur_intermol_ene ;
						eff_ene_file << endl;
					}

					
					if(!play_back_flag)
					{
						if(output_rejected_points)
						{   
							for (int ii=0; ii < nmol; ii++)
							{
								vec_quat_cur[ii].GetQuaternion(qang, qx, qy, qz);
								sprintf(buf,"%6d  %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f",
								np_accept, cur_trans(1,ii+1), cur_trans(2,ii+1), cur_trans(3,ii+1), qang, qx, qy, qz);
								all_points_file << buf << " accepted " << endl;
							}
						}
					}
			}
			if (is_accept) step_move_flag= true;
		} // End of for(imol)

		double acc_ratio;

		if(!play_back_flag && step_move_flag)
		{
			for(imol  = 0; imol < nmol; imol++)
			{
				vec_quat_cur[imol].GetQuaternion(qang, qx, qy, qz);
				sprintf(buf,"%6d  %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f",
					imol+1, cur_trans(1,imol+1), cur_trans(2,imol+1), cur_trans(3,imol+1), qang, qx, qy, qz);
				traj_file << buf << endl;	
			}
		}


		if( clock() > update_time || istep == (num_mc_steps - 1))
		{
			if(!moved_flag)
			{
				for( imol = 0; imol < nmol; imol++)
				{
					if( freeze_first_mol && (imol == 0) )
						continue;
					AtomContainer* pMol = interact_groups[imol];
					Vec3D trans_l;
					trans_l[0] = cur_trans(1,imol+1);
					trans_l[1] = cur_trans(2,imol+1);
					trans_l[2] = cur_trans(3,imol+1);
					pMol->SetQuaternionTrans(vec_quat_cur[imol], trans_l);
					moved_flag = TRUE;
				}
			}				
			if (istep%10 ==0) pmset->AnnounceGeomChange();
			update_time = clock() + CLOCKS_PER_SEC*1 ;
		}
		if (intermol_ene_accepted != 0.0)
		{
			pmset->info_str.clear();
			sprintf(buf,"points accepted %d Total effective energy is %14.6f kcal/mol Replica %d",np_accept,intermol_ene_accepted, ireplica);
			pmset->info_str.push_back(buf);
		}
		else
		{
			pmset->info_str.clear();
			sprintf(buf,"No points accepted, Replica %d", ireplica);
			pmset->info_str.push_back(buf);
		}
	}

	t=0; 
	int jmol;

	for (jmol =0 ; jmol < nmol; jmol++)
	{
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, cur_trans(1,jmol+1));
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, cur_trans(2,jmol+1));
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, cur_trans(3,jmol+1));
		vec_quat_cur[jmol].GetQuaternion(qang, qx, qy,qz);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qang);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qx);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qy);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qz);
	}
	int n_visit = pt_visit.size();
	
	if( np_accept > 0) 
	{
		double acc_ratio = (double)np_accept/( (double) np_accept+ (double) np_reject);
		p_acc_ratio[ireplica] = acc_ratio ;
		if(average_eff_ene != 0) energy_arr[ireplica] = average_eff_ene/np_accept;
	}

	return True;
	
}

*/

int InterMolMCSimulator::RunMCEmpiricalNMA() 
{
	return TRUE;
}

/*

int InterMolMCSimulator::RunMCEmpiricalNMA() 
{	
	if(module_to_init_flag)
	{
		Initialize();
	}
	
	if(!rex_flag) srand( time(NULL) );
	
	Random rand_num_gen( rand() );
	stop_calc_flag = 0;
	int t =0;
	
	bool freeze_first_mol = false;
	//		if( interact_groups.size() > 1) 
	//freeze_first_mol = true;
	
	MolSet* pmset = GetMolSet();
	HaMolMechMod* pmm_mod = pmset->GetMolMechMod(true);
	HaEmpiricalMod* emp_mod = pmset->GetEmpiricalMod(true);
	if(emp_mod -> module_to_init_flag)
	{
		emp_mod -> Initialize();
	}
	
	double intermol_ene_accepted =0;
	
	if(amber_flag)
	{
		if(pmm_mod->module_to_init_flag)
		{	
			//pmm_mod->build_nb_contact_list_flag = FALSE;
			pmm_mod->Initialize();
			//pmm_mod->max_num_minim_steps=5;
			pmm_mod->SaveAllInpFiles();
		}
	}
	
	int nmol = interact_groups.size();
	
	
	vector<Quaternion> vec_quat_old(nmol);
	vector<Quaternion> vec_quat_cur(nmol);
	
	Quaternion q_incr;
	Quaternion q_temp;
	double incr;
	double qang, qx, qy,qz;
	
	clock_t update_time = clock();
	double ene_prev;
	double kt_kcal = MC_temp* BOLTZ * AVOG_NUM/(1000.0*CAL_TO_JOULES); 
	//PrintLog("MC_temp %4.1f \n", MC_temp);  // Remove
	int itr_res = 0;
	int imol;
	
	
	char buf[256];
	
	int iframe = 0;
	
	int ind = 0;
	
	fstream eff_ene_file;
	
	if(dont_calc_ene_flag)
	{
		sprintf(buf,"%s%i%s", MC_energy_file_replica_basename.c_str(), n_playback_replica,".dat");
		eff_ene_file.open(buf, ios::in);
	}
	else
	{
		eff_ene_file.open("eff_ene_file.dat", ios::out);
	}
	
	fstream traj_file;
	fstream all_points_file;
	
	sprintf(buf,"%s%i%s", MC_traj_file_replica_basename.c_str(), ireplica,"NMA.dat");
	traj_file.open(buf, ios::out|ios::app  );
	
	if(output_rejected_points)
	{
		all_points_file.open("all_mc_pts.dat",ios::out);
	}
	
	if( traj_file.fail())
	{
		ErrorInMod("HaInterMolMod::RunMCEmpiricalNMA()",
			"Can't open file to write MC trajectory ");
		return FALSE;
	}
	if( freeze_first_mol)
	{
		AtomContainer* pMol1 = interact_groups[0];
		Vec3D trans1;
		
		pMol1->GetQuaternionTrans(vec_quat_cur[0], trans1);
		vec_quat_cur[0].GetQuaternion(qang, qx, qy,qz);
		
		sprintf(buf,"%6d  %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f",
			-1, trans1[0], trans1[1], trans1[2], qang, qx, qy,qz);
		
		traj_file << buf << endl;
	}
	
	 
	 
    HaMat_double cur_trans(3,nmol);
	Vec3D trans_v;
	
    HaMat_double trans_old(3,nmol);	
	
	
	
	for( imol= 0; imol < nmol; imol++)
	{
		AtomContainer* pMol = interact_groups[imol];
		pMol->GetQuaternionTrans(vec_quat_cur[imol], trans_v);
		
		cur_trans(1,imol+1) = trans_v[0];
		cur_trans(2,imol+1) = trans_v[1];
		cur_trans(3,imol+1) = trans_v[2];
		
	}
	
	int np_accept = 0; 
	int np_reject = 0;
	
	double average_eff_ene = 0.0;
	
	HaVec_double normal_step(nmol*7);
	//	int nmode = nmol*7*0.3;
	HaVec_int  mode_index(nmol*7);
	int nmode = 0;
	double eigval_sum = 0;
	
	double coef;
	
	VecPtr pt_mol;
	pt_mol.resize(nmol);
	int j;
	for(j =0; j < nmol ; j++ )
	{
		AtomContainer* pMol = interact_groups[j];
		pt_mol[j] = pMol;
	}
	
	
	// Previous energy
	if( rex_flag)
	{
		ene_prev = energy_arr[ exchange_arr[ireplica] ];
	}
	else ene_prev = 10e9;
	//	PrintLog("First ene_prev %2.3f \n",ene_prev);
	
	int i;
	if(rex_flag) 
	{
		ind =1;
		for (i= 1; i <= nmol*7-1; i++)
		{
			mode_index(i)=0;
			eigval_sum += normalmode_val(i);
		}
		eigval_sum = eigval_sum/(nmol*7-1);
		
		for (i= 1; i <= nmol*7 -1; i++)
		{
			if (normalmode_val(i) < eigval_sum)
			{
				mode_index(ind) = i;
				ind++;
				nmode++;
			}
		}
		nmode =FMAX((float) nmode, (float) nmol*3);
		coef = (double) 1/(nmode*RAND_MAX);
	}
	
	for( int istep = 0; istep < num_mc_steps ; istep++)       //   MC steps
	{
		bool step_move_flag= false;
		if(!rex_flag && istep%10 == 0) 
		{
			emp_mod ->Neighborhood();
			//			if ( (cur_intermol_ene - min_ene) > 10.0)
			//ENERGY			NormalModes(2, pt_mol);
			NormalModes(3, pt_mol);
			ind =1;
			eigval_sum = 0;
			nmode =0;
			for (i= 1; i <= nmol*7 -1; i++)
			{
				mode_index(i)=0;
				eigval_sum += normalmode_val(i);
			}
			eigval_sum = eigval_sum/(nmol*7-1);
			
			for (i= 1; i <= nmol*7-1; i++)
			{
				if (normalmode_val(i) < eigval_sum)
				{
					mode_index(ind) = i;
					ind++;
					nmode++;
				}
			}
			nmode =FMAX((float) nmode, (float) nmol*3);
			coef = (double) 1/nmode;
			for (i= 1; i <= nmol*7; i++)
			{
				PrintLog("mode_index(%d) = %d , normalmode_val(%d) = %f\n", i, mode_index(i), i, normalmode_val(i));
			}
			PrintLog("nmode %d , eigval_sum %2.1f \n", nmode, eigval_sum);
			
		}
		
		int moved_flag = FALSE;
		
		//		int nmode = nmol*6;
		for (i= 1; i <= nmol*7; i++)
		{
			normal_step(i) = 0.0;
		}
		
		double total_ang_ratio = 0;
		double	total_tr_ratio = 0;
		int imode;
		//		imode = istep +1;
		if (istep%100 ==0) PrintLog("-()-()- NORMAL MODES step %d\n", istep);
		nmode = 14; //REMOVE
		for( imode = 1; imode  <= nmode; imode++)
		{
			ind =0;
			double rand_dbl = rand_num_gen() * coef;
			total_ang_ratio += ang_ratio*rand_dbl;
			total_tr_ratio += tr_ratio*rand_dbl;
			if(rand_dbl >0.5) rand_dbl = - rand_dbl;
			//for (i= 1; i <= nmol; i++)
			for (i= 1; i <= 1; i++)
			{
				normal_step(ind+1) += normalmode_vec(ind+1,imode); //* tr_ratio*rand_dbl;
				ind++;
				normal_step(ind+1) += normalmode_vec(ind+1,imode); //* tr_ratio*rand_dbl;
				ind++;
				normal_step(ind+1) += normalmode_vec(ind+1,imode); //* tr_ratio*rand_dbl;
				ind++;
				normal_step(ind+1) += normalmode_vec(ind+1,imode); // *tr_ratio*rand_dbl;
				ind++;
				normal_step(ind+1) += normalmode_vec(ind+1,imode); // *tr_ratio*rand_dbl;
				ind++;
				normal_step(ind+1) += normalmode_vec(ind+1,imode); // *tr_ratio*rand_dbl;
				ind++;
				normal_step(ind+1) += normalmode_vec(ind+1,imode); // *tr_ratio*rand_dbl;
				ind++;
				//PrintLog("tr_ratio*rand_dbl = %f  %f %f\n", tr_ratio*rand_dbl, tr_ratio, rand_dbl);
			}
		}
			for (i= 1; i <= nmol*7; i++)
			{
				PrintLog("normal_step(%d) = %f\n", i, normal_step(i));
			}

		ind=0;
		
		if(!play_back_flag )
		{
			for( imol = 0; imol < nmol; imol++) // Old position
			{
				vec_quat_old[imol] = vec_quat_cur[imol];
				
				trans_old(1,imol+1) = cur_trans(1,imol+1);
				trans_old(2,imol+1) = cur_trans(2,imol+1);
				trans_old(3,imol+1) = cur_trans(3,imol+1);
			}	
			ind=0;
			for(  imol= 0; imol < nmol; imol++)
			{
				q_incr = vec_quat_cur[imol];
				qang = normal_step(ind+1);
				ind++;
				qx = normal_step(ind+1);
				ind++;
				qy = normal_step(ind+1);
				ind++;
				qz = normal_step(ind+1);
				ind++;
				
				q_temp.SetQuaternion(qang, qx, qy, qz);
				q_temp.Normalize();
//				q_incr.operator *=(q_temp);
				vec_quat_cur[imol] = q_incr;
				
				cur_trans(1,imol+1) = cur_trans(1,imol+1) + normal_step(ind+1);
				ind++;
				cur_trans(2,imol+1) = cur_trans(2,imol+1) + normal_step(ind+1);
				ind++;
//				cur_trans(3,imol+1) = cur_trans(3,imol+1) + normal_step(ind+1);
				ind++;
			}
			
			
			
			if(!dont_calc_ene_flag)
			{
				int ic = 0;
				if(!discr_moves_flag || ic == 0)
				{
					Vec3D trans_l;
					for(  imol= 0; imol < nmol; imol++)
					{
						AtomContainer* pMol = interact_groups[imol];
						trans_l[0] = cur_trans(1,imol+1);
						trans_l[1] = cur_trans(2,imol+1);
						trans_l[2] = cur_trans(3,imol+1);
						
						pMol->SetQuaternionTrans(vec_quat_cur[imol], trans_l);	
					}
					moved_flag = TRUE;
					
					//						SteepestDescentMinimizer(50*nmol);	
					// cur_intermol_ene = emp_mod->ScoreEnergy();
					cur_intermol_ene = emp_mod-> HarmonicEnergy();

					if((istep+1)%1 == 0) PrintLog("Score Energy %10.3f imol %d\n", cur_intermol_ene, imol);
				}
				
			}
			
			bool is_accept = false;	
			
			if(!play_back_flag)
			{
				if( cur_intermol_ene < ene_prev ) 
				{
					is_accept = true;
				}
				else
				{
					double rand_dbl = rand_num_gen(); // rand()/RAND_MAX;
					double dtest = exp( (ene_prev - cur_intermol_ene)/kt_kcal );
					//double dtest = exp( (ene_prev - cur_intermol_ene)/MC_temp ) ;
					//	PrintLog("ene_prev %8.4f, cur_intermol_ene %8.4f, rand_dbl %8.4f, DTEST= %6.4f \n" , ene_prev , cur_intermol_ene , rand_dbl ,dtest);  // Remove
					if( dtest > rand_dbl )
					{
						is_accept = true;
					}
					else
					{
						//	PrintLog(" Point is rejected \n");
						np_reject++;
						
						if(output_rejected_points)
						{
							vec_quat_cur[imol].GetQuaternion(qang, qx, qy, qz);
							sprintf(buf,"%6d  %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f",
								np_accept, cur_trans(1,imol+1), cur_trans(2,imol+1), cur_trans(3,imol+1), qang, qx, qy, qz);
							all_points_file << buf << " rejected " << endl;
						}
						
						ind =0;						
						if(moved_flag)
						{
							for(  imol= 0; imol < nmol; imol++)
							{
								
								AtomContainer* pMol = interact_groups[imol];
								Vec3D trans_o;
								trans_o[0] = trans_old(1,imol+1);
								trans_o[1] = trans_old(2,imol+1);
								trans_o[2] = trans_old(3,imol+1);
								
								pMol->SetQuaternionTrans(vec_quat_old[imol], trans_o);
							}
						}
						
						for(  imol= 0; imol < nmol; imol++)
						{
							
							vec_quat_cur[imol] = vec_quat_old[imol];
							cur_trans(1,imol+1) = trans_old(1,imol+1);
							cur_trans(2,imol+1) = trans_old(2,imol+1);
							cur_trans(3,imol+1) = trans_old(3,imol+1);
						}
					}
				}
			}
			

			if( (is_accept || play_back_flag))
			{
				pmset->info_str.clear();
				
				ene_prev = cur_intermol_ene;
				intermol_ene_accepted = cur_intermol_ene;
				np_accept++;
				
				
				if( !dont_calc_ene_flag)
				{
					average_eff_ene+= cur_intermol_ene; 
					eff_ene_file << np_accept << "  " << cur_intermol_ene ;
					eff_ene_file << endl;
				}
				
				
				if(!play_back_flag)
				{
					if(output_rejected_points)
					{
						all_points_file << buf << " accepted " << endl;
					}
				}
			}
			if (is_accept) step_move_flag= true;
			
		}	
			
		if(!play_back_flag && step_move_flag)
		{
			for(imol  = 0; imol < nmol; imol++)
			{
				vec_quat_cur[imol].GetQuaternion(qang, qx, qy, qz);
				sprintf(buf,"%6d  %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f",
					imol+1, cur_trans(1,imol+1), cur_trans(2,imol+1), cur_trans(3,imol+1), qang, qx, qy, qz);
				traj_file << buf << endl;	
			}
		}
		
		
		if( clock() > update_time || istep == (num_mc_steps - 1))
		{
			if(!moved_flag)
			{
				for( imol = 0; imol < nmol; imol++)
				{
					if( freeze_first_mol && (imol == 0) )
						continue;
					AtomContainer* pMol = interact_groups[imol];
					Vec3D trans_l;
					trans_l[0] = cur_trans(1,imol+1);
					trans_l[1] = cur_trans(2,imol+1);
					trans_l[2] = cur_trans(3,imol+1);
					pMol->SetQuaternionTrans(vec_quat_cur[imol], trans_l);
					moved_flag = TRUE;
				}
			}				
			if (istep%10 ==0) pmset->AnnounceGeomChange();
			update_time = clock() + CLOCKS_PER_SEC*1 ;
		}
		if (intermol_ene_accepted != 0.0)
		{
			pmset->info_str.clear();
			sprintf(buf,"points accepted %d Total effective energy is %14.6f kcal/mol Replica %d",np_accept,intermol_ene_accepted, ireplica);
			pmset->info_str.push_back(buf);
		}
		else
		{
			pmset->info_str.clear();
			sprintf(buf,"No points accepted, Replica %d", ireplica);
			pmset->info_str.push_back(buf);
		}
	} // Enr of istep loop
	
	//	sprintf(buf," Total effective energy is %14.6f kcal/mol ",intermol_ene_accepted);
	//	PrintMessage(buf);
	//	pmset->info_str.push_back(buf);
	
	t=0; 
	int jmol;
	
	for (jmol =0 ; jmol < nmol; jmol++)
	{
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, cur_trans(1,jmol+1));
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, cur_trans(2,jmol+1));
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, cur_trans(3,jmol+1));
		vec_quat_cur[jmol].GetQuaternion(qang, qx, qy,qz);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qang);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qx);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qy);
		position_mat.SetVal_idx0(exchange_arr[ireplica],t++, qz);
	}
	
	int n_visit = pt_visit.size();
	//PrintLog("Number of points visited %d \n", n_visit); 
	
	if( np_accept > 0) 
	{
		double acc_ratio = (double)np_accept/( (double) np_accept+ (double) np_reject);
		//PrintLog(" point acceptance ratio: %7.3f \n", acc_ratio );  
		//PrintLog(" average effective energy: %12.6e \n", average_eff_ene/np_accept );
		p_acc_ratio[ireplica] = acc_ratio ;
		if(average_eff_ene != 0) energy_arr[ireplica] = average_eff_ene/np_accept;
	}
	
	return True;
	
}

*/

bool HaInterMolMod::CalcEffInterEne()
{
	cur_intermol_ene = 0.0;
	double el_inter_ene = 0.0;

	MolSet* pmset = GetMolSet();

	if(p_mc_sim->amber_flag )
	{
		HaMolMechMod* pmm_mod = pmset->GetMolMechMod();	
		if( pmm_mod == NULL) 
		{
			PrintLog(" Error in HaInterMolMod::CalcEffInterEne() \n");
			PrintLog(" Molecular Mechanics Module is not initialized \n");
			return false;
		}

		pmm_mod->SetRunType("MIN_RUN");
		pmm_mod->p_amber_driver->SaveAmberInpFile();
		pmm_mod->p_amber_driver->SaveAmberCrdFile();
		pmm_mod->Run();
		std::this_thread::sleep_for(std::chrono::milliseconds(2000));
		pmm_mod->SetRunType("MD_RUN");
		pmm_mod->p_amber_driver->SaveAmberInpFile();
		pmm_mod->p_amber_driver->SaveAmberCrdFile();
		pmm_mod->Run();
		std::this_thread::sleep_for(std::chrono::milliseconds(200));
	}

	if( electr_model == CHARGES_IN_FIELD_ELECTR )
	{
		cur_intermol_ene = CalcChargesInFieldEne();
	}
	else
	{
		CalculateMMEnergy();
		
		if( electr_model == CONTINUUM_ELECTR )
			el_inter_ene = CalcElStaticInter();
		
		cur_intermol_ene+= el_inter_ene * KT_300_KCAL;
	}

	add_eff_ene = 0.0;
	if(calc_et_rate)
	{
		ETCouplMod* ptr_et_coupl_mod= pmset->GetETCouplMod();

		int old_log_et = ptr_et_coupl_mod->log_calc_result;
		ptr_et_coupl_mod->log_calc_result = 0;
		ptr_et_coupl_mod->pathways_graph_init_flag = false;
		if( debug_level > 5)
			PrintLog(" pw_nb_decay = %12.6f",ptr_et_coupl_mod->pw_nb_decay);
		//				bool res = ptr_et_coupl_mod->path_coupl_calc();
		bool res = ptr_et_coupl_mod->calc_intermol_path_coupl();
		ptr_et_coupl_mod->log_calc_result = old_log_et;
		
		if( !res || (ptr_et_coupl_mod->best_path_coupl <  DBL_EPSILON) ) // don't allow zero coupling
			add_eff_ene = 10.0e5;
		else
		{
			if(debug_level > 5)
				PrintLog(" Pathways coupling value = %12.6e \n ", ptr_et_coupl_mod->best_path_coupl);
			add_eff_ene = -2*log(ptr_et_coupl_mod->best_path_coupl)*KT_300_KCAL;
		}
		cur_intermol_ene += add_eff_ene;
	}
	
	return true;
}

int 
InterMolMCSimulator::RunQuasiREM() 
{
	return TRUE;
}

/*

int 
InterMolMCSimulator::RunQuasiREM() 
{
	srand(time(NULL));
	MolSet* pmset = GetMolSet();
	int temp = 0 ; 
	double betta = 0.0;
	double temp_ene =0.0;
	int nmol = pmset->GetNMol();
	HaVec_double temperature_arr ;
	int j = 0 ;
	HaMat_double stats ;
	char buf[256];
	fstream replica_file;
	replica_file.open("replica_exchange.dat",ios::out );
	fstream traj_file;
	fstream ene_file;
	fstream rst_file;
	double ene_low = 10e9;
	double rem_acc_ratio = 1;
	double rem_acc_ratio1 = 1;
	double MC_temp_old = MC_temp;
	double temperature_max_old = temperature_max;
	time_t t1 = time(NULL);
	clock_t tt1 = clock();
    double  rand_dbl;
	double rand_expo;
	int n;
	int imol;
	double acc_ratio_coef =0;

	//HaEmpiricalMod* emp_mod = pmset->GetEmpiricalMod(true); jose
	HaMolMembraneMod* mem_mod = pmset->GetMolMembraneMod(true);
	
	for( ireplica= 0; ireplica < nreplicas; ireplica++) 
	{
		sprintf(buf,"%s%i%s", MC_traj_file_replica_basename.c_str(), ireplica,".dat");
		traj_file.open(buf,ios::out );
		sprintf(buf,"%s%i%s", MC_energy_file_replica_basename.c_str(), ireplica,".dat");
		ene_file.open(buf,ios::out );
		traj_file.close();
		ene_file.close();
	}
	
	if(module_to_init_flag)
	{
		Initialize();
	}

	temperature_arr.newsize(nreplicas);
	stats.newsize(nreplicas,rem_steps);
	position_mat.newsize(nreplicas,nmol*7);
	energy_arr.newsize(nreplicas);
	
	if (nreplicas == 1)
	{
		betta = log(temperature_max/MC_temp)/nreplicas;
		for ( n=0; n < nreplicas; n++)
		{
			temperature_arr[n]= MC_temp*exp(betta*n) ;
			//		temperature_arr[n]= MC_temp;
			sprintf(buf,"Replica %2d Temperature %2.1f", n, temperature_arr[n]);
			replica_file << buf << endl;
			//	PrintLog("temperature_arr[n] %3.3f, %3.3f\n", temperature_arr[n], betta);
		}
	}
	else
	{
		betta = log(temperature_max/MC_temp)/(nreplicas-1);
		for (n=0; n < nreplicas; n++)
		{
			temperature_arr[n]= MC_temp*exp(betta*n) ;
			//		temperature_arr[n]= MC_temp;
			sprintf(buf,"Replica %2d Temperature %2.1f", n, temperature_arr[n]);
			replica_file << buf << endl;
			//	PrintLog("temperature_arr[n] %3.3f, %3.3f\n", temperature_arr[n], betta);
		}
		
	}
	Vec3D trans_vec; // Array of translational coordinates of molecules
	vector<Quaternion> vec_quat_cur(nmol); // Arrays of rotational coordinates of molecules
	double qang, qx, qy, qz;
	


	int i;
	for (ireplica= 0; ireplica < nreplicas; ireplica++)
	{
		j=0;
		for(imol=0 ; imol < nmol; imol++)
		{
			HaMolecule* pMol = pmset->GetMolByIdx(imol);
			pMol->GetQuaternionTrans(vec_quat_cur[imol],trans_vec);
			position_mat.SetVal_idx0(ireplica,j++, trans_vec[0]);
			position_mat.SetVal_idx0(ireplica,j++, trans_vec[1]);
			position_mat.SetVal_idx0(ireplica,j++,trans_vec[2]);
			vec_quat_cur[imol].GetQuaternion(qang, qx, qy, qz);
			position_mat.SetVal_idx0(ireplica,j++,qang);
			position_mat.SetVal_idx0(ireplica,j++,qx);
			position_mat.SetVal_idx0(ireplica,j++,qy);
			position_mat.SetVal_idx0(ireplica,j++,qz);
		} 
		for(i= 0; i < rem_steps ; i++) 
		{
			stats.SetVal_idx0(ireplica,i,0);
		}
	}
	
	for( i= 0; i < nreplicas ; i++) 
	{
		exchange_arr[i] = i;
		stats[i] = 0;
		p_acc_ratio[i] = 0.5;
		energy_arr[i] = 10e9;
	}
	
	// Read coordinates from restart file (if it exists)
	wxString str;
	char* cline = NULL;
	sprintf(buf,"%s%s", MC_rst_file_basename.c_str(),".dat");
	FILE* finfo= fopen(buf,"r");
	if(finfo == NULL)
	{
		PrintLog("RESTART_MC.DAT does not exist\n");
	}
	else
	{
		cline = fgets(buf,255,finfo); 
		str = buf;  
		wxArrayString sub_str;
		wxStringTokenizer tkz(str," ");
		while ( tkz.HasMoreTokens())
		{
			wxString token = tkz.GetNextToken();
			sub_str.Add(token);
		}

		long nrepl, nmlc;
		bool res =sub_str[0].ToLong(&nrepl);
		res =sub_str[1].ToLong(&nmlc);
		for (ireplica= 0; ireplica < nrepl; ireplica++)
		{
			j=0;
			for( int imol= 0; imol < nmlc; imol++)
			{
				if (imol > nmol-1) continue;

				cline = fgets(buf,255,finfo); 
				
				str = buf;  
				wxArrayString sub_str;
				wxStringTokenizer tkz(str," ");
				while ( tkz.HasMoreTokens())
				{
					wxString token = tkz.GetNextToken();
					sub_str.Add(token);
				}
				
				bool res1 = sub_str[0].ToDouble(&trans_vec[0]) ;
				bool res2 = sub_str[1].ToDouble(&trans_vec[1]) ;
				bool res3 = sub_str[2].ToDouble(&trans_vec[2]) ;
				bool res4 = sub_str[3].ToDouble(&qang) ;
				bool res5 = sub_str[4].ToDouble(&qx) ;
				bool res6 = sub_str[5].ToDouble(&qy) ;
				bool res7 = sub_str[6].ToDouble(&qz) ;
				position_mat.SetVal_idx0(ireplica,j++, trans_vec[0]);
				position_mat.SetVal_idx0(ireplica,j++, trans_vec[1]);
				position_mat.SetVal_idx0(ireplica,j++, trans_vec[2]);
				position_mat.SetVal_idx0(ireplica,j++, qang);
				position_mat.SetVal_idx0(ireplica,j++, qx);
				position_mat.SetVal_idx0(ireplica,j++, qy);
				position_mat.SetVal_idx0(ireplica,j++, qz);
			}
		}
		fclose(finfo);
	}

	VecPtr pt_mol;
	pt_mol.resize(nmol);
	int i_step = 0;
	int n_steps = FALSE;
	int temp_replica = 0;
	for(j =0; j < nmol ; j++ )
	{
		AtomContainer* pMol = interact_groups[j];
		pt_mol[j] = pMol;
	}
	// Run loop
	int jstep;
	tt1 = clock();	
	for (jstep =0 ; jstep < rem_steps; jstep++)
	{
		//PrintLog("jstep = %d\n", jstep);
		int result;
		if (nreplicas > 1 && jstep >0 )
		{
			if(vary_temperature_flag)
			{
				int rem_acc_ratio_flag = FALSE;
				if(jstep > 0.1*rem_steps)
				{
					for ( ireplica =0 ; ireplica < nreplicas-1; ireplica++)
					{
						int nsteps = (int)(jstep - 0.1*rem_steps);
						for(i = nsteps; i <= jstep; i++)
						{
							rem_acc_ratio += stats.GetVal_idx0(ireplica+1,i);
						}
						rem_acc_ratio = rem_acc_ratio/(0.1*rem_steps);
						if(rem_acc_ratio < 0.1 )
						{
							rem_acc_ratio_flag = TRUE;
						}
					}
				}
				if(rem_acc_ratio_flag)
				{
					temperature_max = temperature_max - (temperature_max_old -MC_temp_old)*0.1;
					//					MC_temp = temperature_arr[0] + (temperature_max_old -MC_temp_old)*0.1;
					if(temperature_max/MC_temp_old < 1.0 )
					{
						//						MC_temp = temperature_arr[0];
						//						temperature_max = temperature_arr[0];
						temperature_max = MC_temp_old;
					}
					betta = log(temperature_max/MC_temp_old)/(nreplicas-1);
					for (n=0; n < nreplicas; n++)
					{
						temperature_arr[n]= MC_temp_old*exp(betta*n) ;
						//		temperature_arr[n]= MC_temp;
						//	sprintf(buf,"# %2d Temperature %2.1f", n, temperature_arr[n]);
						//		PrintLog(" jstep %d, %d stats[nreplicas-1]/jstep %3.2f \n", jstep, stats[nreplicas-1], (double)stats[nreplicas-1]/jstep);
						//	replica_file << buf << endl;
						//	PrintLog("temperature_arr[n] %3.3f, %3.3f\n", temperature_arr[n], betta);
					}
				}
				else
				{
					temperature_max = temperature_max_old;
					MC_temp = MC_temp_old;
					betta = log(temperature_max/MC_temp)/(nreplicas-1);
					for (n=0; n < nreplicas; n++)
					{
						temperature_arr[n]= MC_temp*exp(betta*n) ;
						//		temperature_arr[n]= MC_temp;
						//		sprintf(buf,"Replica %2d Temperature %2.1f", n, temperature_arr[n]);
						//		replica_file << buf << endl;
						//	PrintLog("temperature_arr[n] %3.3f, %3.3f\n", temperature_arr[n], betta);
					}
				}
			}
		}
//		clock_t tt2 = clock();
//		sprintf(buf,"Time tt2 %2.2f", (double)(tt2-tt1)/CLOCKS_PER_SEC);
//		replica_file << buf << endl;
		
		for ( ireplica =0 ; ireplica < nreplicas; ireplica++)
		{
//			clock_t tt3 = clock();
//			fstream component_file;
//			component_file.open("energy_components.dat",ios::out | ios::app);
//			sprintf(buf,"replica %d step %d", ireplica,jstep);
//			component_file << buf << endl;
//			component_file.close();
			
			if (jstep% 100 ==0 )PrintLog("\n##### REM step %i #### Replica %i (%i) ##################\n\n", jstep, ireplica,exchange_arr[ireplica]);
			j =0;
			MC_temp = temperature_arr[ireplica];
			for( imol= 0; imol < nmol; imol++)
			{
				HaMolecule* pMol = pmset->GetMolByIdx(imol);
				trans_vec[0] = position_mat.GetVal_idx0(exchange_arr[ireplica], j++);
				trans_vec[1] = position_mat.GetVal_idx0(exchange_arr[ireplica], j++);
				trans_vec[2] = position_mat.GetVal_idx0(exchange_arr[ireplica], j++);
				qang = position_mat.GetVal_idx0(exchange_arr[ireplica], j++);
				qx  = position_mat.GetVal_idx0(exchange_arr[ireplica], j++);
				qy  = position_mat.GetVal_idx0(exchange_arr[ireplica], j++);
				qz  = position_mat.GetVal_idx0(exchange_arr[ireplica], j++);
				vec_quat_cur[imol].SetQuaternion(qang,qx,qy,qz);
				result = pMol-> SetQuaternionTrans(vec_quat_cur[imol],trans_vec);
				//PrintLog("SETqang1 %3.2f, %3.2f,  %3.2f  %3.2f\n ", qang, qx, qy,qz);
				//PrintLog("1 %3.1f,2 %3.1f, 3 %3.1f \n ",position_mat.GetVal_idx0(exchange_arr[ireplica], j-3), position_mat.GetVal_idx0(exchange_arr[ireplica], j-2), position_mat.GetVal_idx0(exchange_arr[ireplica], j-1));
				//				PrintLog("rot_v[0] %3.1f, rot_v[1] %3.1f, rot_v[2] %3.1f, j %i \n", c_phi(imol+1), c_cos_theta(imol+1),c_psi(imol+1),  j );  	
			}
			//if ( (energy_arr[exchange_arr[ireplica]] - ene_low) > 100.0)
			//{
				//Minimizer(2, pt_mol);
			//}

//			if (jstep>0 && jstep%10 == 0 && ireplica ==  nreplicas-1)
//			{
//				NormalModes(2, pt_mol);
//			}

			//if (emp_mod -> module_to_init_flag) emp_mod -> Initialize(); jose
			//emp_mod -> Neighborhood();
			//if (mem_mod->module_to_init_flag) mem_mod -> Initialize(); jose
			//emp_mod -> Neighborhood();
			
//			clock_t tt4 = clock();
//			sprintf(buf,"Time tt3 %2.2f replica %d", (double)(tt4-tt3)/CLOCKS_PER_SEC, ireplica );
//			replica_file << buf << endl;

			
			rand_dbl = (double) rand()/RAND_MAX;
			if(rand_dbl <1.0) rand_expo = -log(1.0 - rand_dbl) ;
			else rand_expo = 5.0;
			//			tr_ratio = rand_expo*rand_expo*p_acc_ratio[exchange_arr[ireplica]];
			if(p_acc_ratio[exchange_arr[ireplica]] > 0.5 )
			{
				acc_ratio_coef = acc_ratio_coef + 0.1;
			}
			else acc_ratio_coef = acc_ratio_coef - 0.1 ;

			if (acc_ratio_coef > 100.0) // 100.0 jose
			{
				acc_ratio_coef = 0.0;
				printf("acc_ratio_coef MORE  %2.1f\n", acc_ratio_coef, p_acc_ratio[exchange_arr[ireplica]]);
			}
			if (acc_ratio_coef < 0.0)
			{
				acc_ratio_coef = 0.0;
				printf("acc_ratio_coef LESS  %2.1f %2.1f %d\n", acc_ratio_coef, p_acc_ratio[exchange_arr[ireplica]], exchange_arr[ireplica]);
			}

			tr_ratio = rand_expo*rand_expo + acc_ratio_coef;

			rand_dbl = (double) rand()/RAND_MAX;
			if(rand_dbl <1.0) rand_expo = -log(1.0 - rand_dbl) ;
			else rand_expo = 5.0;

			ang_ratio = (rand_expo*rand_expo*p_acc_ratio[exchange_arr[ireplica]]+ acc_ratio_coef*0.01)*0.1; 
			//		sprintf(buf,"%s%i%s", "MC_replica", ireplica,jstep,".hlm");
			//		pmset->SaveHarlemFile(buf);
			if(empirical_flag)
			{
				
					if (xy_mc_flag) RunMCEmpiricalXY();
					else RunMCQuantSampling(); // This should be a function with number of MC steps as input parameters
			}
			else
			{
				RunMCDock();
			}


			if (energy_arr[ireplica] < ene_low)
			{
				sprintf(buf,"%s%i_step%i%s", "MC_replica", ireplica,jstep,".hlm");
				//pmset->SaveHarlemFile(buf);
				ene_low = energy_arr[ireplica] ;
			}
//			clock_t tt5 = clock();
			
			//time_t t3=time(NULL);
			//printf("%d seconds elapsed between  t3-t1\n", t3-t1);
			j =0;
			for( imol= 0; imol < nmol; imol++)
			{
				trans_vec[0] = position_mat.GetVal_idx0(exchange_arr[ireplica], j++);
				trans_vec[1] = position_mat.GetVal_idx0(exchange_arr[ireplica], j++);
				trans_vec[2] = position_mat.GetVal_idx0(exchange_arr[ireplica], j++);
				qang = position_mat.GetVal_idx0(exchange_arr[ireplica], j++);
				qx  = position_mat.GetVal_idx0(exchange_arr[ireplica], j++);
				qy  = position_mat.GetVal_idx0(exchange_arr[ireplica], j++);
				qz  = position_mat.GetVal_idx0(exchange_arr[ireplica], j++);

				sprintf(buf,"%s%i%s", MC_traj_file_replica_basename.c_str(), ireplica,".dat");
				traj_file.open(buf,ios::out | ios::app );
				sprintf(buf,"%6d  %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f", nmol*jstep+imol,
					trans_vec[0], trans_vec[1], trans_vec[2], qang, qx, qy, qz);
				traj_file << buf << endl;
				traj_file.close();	
			}
			
			sprintf(buf,"%s%i%s", MC_energy_file_replica_basename.c_str(), ireplica,".dat");
			ene_file.open(buf,ios::out | ios::app);
			sprintf(buf,"%2d %5.3f",jstep, energy_arr[ireplica] );
			ene_file << buf << endl;
			ene_file.close();
			sprintf(buf,"3.1f %3.2f %d %3.2f %3.1f %3.2f %3.2f", temperature_arr[ireplica],energy_arr[ireplica], exchange_arr[ireplica],p_acc_ratio[exchange_arr[ireplica]],acc_ratio_coef, tr_ratio, ang_ratio);
			replica_file << buf << endl;
		} // END OF for(ireplica)
		
		for( ireplica= 0; ireplica < nreplicas -1; ireplica++) 
		{
			double deltaE = energy_arr[ireplica+1] - energy_arr[ireplica];
			double deltaBetta = (1.0/temperature_arr[ireplica+1])-(1.0/temperature_arr[ireplica]);
			rand_dbl = (double) rand()/RAND_MAX;
			double wFactor = - deltaBetta * deltaE ; 
			double dtest = exp(-wFactor);
			//				PrintLog("wFactor %3.2f rand_Number %3.2f dtest %3.2f &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n", wFactor,rand_dbl,dtest);
			if(wFactor < 0.0 || rand_dbl <  dtest )
			{
				temp= exchange_arr[ireplica];
				exchange_arr[ireplica] = exchange_arr[ireplica+1];
				exchange_arr[ireplica+1] = temp;
				stats.SetVal_idx0(ireplica+1,jstep,1.0);
			}
			
			//			PrintLog("ireplica %d, stats[ireplica+1] %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5\n", ireplica, stats[ireplica+1]);
		}

		//		sprintf(buf,"Time tt5 %2.2f", (double)(tt7-tt6)/CLOCKS_PER_SEC);
		
	 	//		for( ireplica= 0; ireplica < nreplicas -1; ireplica++) 
		//		{
		//			temp_ene = energy_arr[exchange_arr[ireplica]];
		//			energy_arr[exchange_arr[ireplica]]= energy_arr[exchange_arr[ireplica+1]];
		//			energy_arr[exchange_arr[ireplica+1]]= temp_ene;
		//		}
//		if (jstep == (rem_steps-1) || jstep%10 == 0)
		if ( jstep == (rem_steps-1) )
		{
			sprintf(buf,"%s.dat", MC_rst_file_basename.c_str());
			rst_file.open(buf,ios::out);
			sprintf(buf,"%d %d", nreplicas, nmol );
			rst_file << buf << endl; 
			for( ireplica= 0; ireplica < nreplicas; ireplica++) 
			{
				j=0;
				for( imol= 0; imol < nmol; imol++)
				{
					for( i= 0; i < 7; i++)
					{
						sprintf(buf,"%14.9f", position_mat.GetVal_idx0(exchange_arr[ireplica], j++) );
						rst_file << buf;
					}
					rst_file << endl;
				}
			}
			rst_file.close();
		}
		if (!result)
		{
			stop_calc_flag =1;
			break;
		}
		if ( jstep%10 == 0) pmset->AnnounceGeomChange();	
} // END OF for(jstep)

for( ireplica= 0; ireplica < nreplicas-1; ireplica++) 
{
	rem_acc_ratio = 0.0;
	for(i = 0; i < rem_steps; i++)
	{
		rem_acc_ratio += stats.GetVal_idx0(ireplica+1,i);
	}
	sprintf(buf,"Replicas %d - %d #ofExchanges %3.0f", ireplica, ireplica+1, rem_acc_ratio);
	replica_file << buf << endl;
}
sprintf(buf,"Number of REM steps %3d", rem_steps);
replica_file << buf << endl;
time_t t4=time(NULL);
sprintf(buf,"Time elapsed %d min. (%d sec.) \n", (t4-t1)/60, (t4-t1));
replica_file << buf << endl;
clock_t tt2 = clock();
	printf("TIME OF THE RUN: %2.2f  sec.\n", (double)(tt2-tt1)/CLOCKS_PER_SEC);

// Check random numbers
//	sprintf(buf,"step LOG(1-rand()   rand_num_gen() LOG(1-rand_num_gen())");
//	replica_file << buf << endl;
//	n=0;
//	for( j= 0; j < 10; j++) 
//	{
//	srand(time(NULL));
//	Random rand_num_gen(time(NULL));
//	std::this_thread::sleep_for(std::chrono::milliseconds(2000));
//
//	for( i= 0; i < 100; i++) 
//	{
//		rand_dbl = (double) rand()/RAND_MAX;
//		rand_expo = -log(1.0 - rand_dbl) ;
//		sprintf(buf,"%d %2.3f %2.3f %2.3f", ++n, rand_expo, rand_num_gen(),-log(1-rand_num_gen()));
//		replica_file << buf << endl;
//	}
//	}
replica_file.close();

return TRUE;
}

*/


int
HaInterMolMod::NormalModes(int energy_type, VecPtr ptmol)
{
	std::fstream eigen;
	eigen.open("eigenvector.dat", std::ios::out );

	if(module_to_init_flag)
	{
		Initialize();
	}
	MolSet* pmset = GetMolSet();

//  HaEmpiricalMod* emp_mod = pmset->GetEmpiricalMod(true);
	//	if (emp_mod ->inertia_axes.GetVal_idx0(0,0) == NULL )	emp_mod -> FindAxes();

	double incr = 0;
	int nmol = ptmol.size();
	int n_size = nmol*6;
	Vec3D trans_v;

	normalmode_vec.newsize(n_size,n_size);
	normalmode_val.newsize(n_size);
    int ind = 0;
	
	ind=0;
	HaVec_double phi_o(nmol), cos_theta_o(nmol), psi_o(nmol); 
	
	HaVec_double coord_old;
	coord_old.newsize(n_size);
	ind=0;
	
	HaMat_double hessian(n_size,n_size);
	PrintLog("Calculating Hessian ... ... ... \n");
	hessian = Hessian(energy_type, ptmol);
	
	HaMat_double::mat_sdiag(hessian, normalmode_vec, normalmode_val);
/*	PrintLog("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW\n");

	for( i = 1; i <= n_size; i++)
	{
		PrintLog("eig_val(%d)=%3.2f \n",i,normalmode_val(i));
		sprintf(buf,"eig_val(%d)=%3.2f ",i,normalmode_val(i));
		eigen << buf << endl;
	}
	PrintLog("Mode1  Mode2    Mode3    Mode4    Mode5    Mode6\n" );
	for( i = 1; i <= n_size; i++)
	{
		for( j = 1; j <= n_size; j++)
		{
			PrintLog("%4.6f ", normalmode_vec(i,j) );
			sprintf(buf,"%4.6f ", normalmode_vec(i,j) );
			eigen << buf ;
		}
		PrintLog("\n");
		sprintf(buf,"\n");
		eigen << buf ;
	}
	
	
	HaMat_double trans_o(3,nmol);	
	for(  i= 1; i <= n_size; i++)
	{
		ind=0;
		for(  imol= 0; imol < nmol; imol++)
		{
			pMol = (HaMolecule*) ptmol[imol];
			pMol->GetPosEulerTrans( cur_phi, cur_cos_theta, cur_psi,trans_v);
			
			phi_o(imol+1) = cur_phi;
			psi_o(imol+1) = cur_psi;
			cos_theta_o(imol+1) = cur_cos_theta;
			
			trans_o(1,imol+1) = trans_v[0];
			trans_o(2,imol+1) = trans_v[1];
			trans_o(3,imol+1) = trans_v[2];
			
			cur_phi = cur_phi + normalmode_vec(ind+1,i);
			//			PrintLog("ind%d cc %2.3f cur_cos_theta %2.3f \n", ind, cc(ind+1,i),cur_cos_theta);
			ind++;
			
			double cos_theta_new; 
			cos_theta_new = cur_cos_theta + normalmode_vec(ind+1,i);
			if(cos_theta_new > 1.0)
			{
				cos_theta_new = cur_cos_theta - normalmode_vec(ind+1,i);
			}
			else if(cos_theta_new < -1.0)
			{
				cos_theta_new = cur_cos_theta - normalmode_vec(ind+1,i);
			}
			ind++;
			
			cur_psi = cur_psi + normalmode_vec(ind+1,i);
			ind++;
			
			trans_v[0] = trans_v[0] + normalmode_vec(ind+1,i)*15;
			//			PrintLog("trans%d cc %2.3f%\n", ind, cc(ind+1,i));
			ind++;
			trans_v[1] = trans_v[1] + normalmode_vec(ind+1,i)*15;
			//			PrintLog("trans%d cc %2.3f%\n", ind, cc(ind+1,i));
			ind++;
			trans_v[2] = trans_v[2] + normalmode_vec(ind+1,i)*15;
			//			PrintLog("trans%d cc %2.3f%\n", ind,cc(ind+1,i));
			ind++;
			//			PrintLog("angles %2.3f %2.3f %2.3f\n", cur_phi, cur_cos_theta, cur_psi);
			pMol->SetPosEulerTrans(cur_phi, cos_theta_new, cur_psi,trans_v);
		}
		sprintf(buf,"Mode%i%s", i,".pdb");
		pmset->SavePDBFile(buf);
		for(  imol= 0; imol < nmol; imol++)
		{
			pMol = (HaMolecule*) ptmol[imol];
			trans_v[0]= trans_o(1,imol+1) ;
			trans_v[1]= trans_o(2,imol+1);
			trans_v[2]= trans_o(3,imol+1);
			pMol->SetPosEulerTrans(phi_o(imol+1), cos_theta_o(imol+1), psi_o(imol+1),trans_v);
		}
		pmset->RefreshAllViews(RFApply | RFRefresh);
	}
*/
	eigen.close();
	return TRUE;
}


HaMat_double
HaInterMolMod::Hessian(int energy_type, VecPtr ptmol)
{
	std::fstream eigen;
	eigen.open("eigenvector.dat", std::ios::out );
	std::fstream check;
	check.open("check_theta.dat", std::ios::out| std::ios::app);


	MolSet* pmset = GetMolSet();
	HaMolecule* pMol;
	HaEmpiricalMod* emp_mod = pmset->GetEmpiricalMod(true);
	if(module_to_init_flag)
	{
		Initialize();
	}
	if(emp_mod -> module_to_init_flag)
	{
		emp_mod ->Initialize();
	}

	typedef double (HaEmpiricalMod::*EnergyFuncPtn)();
	EnergyFuncPtn func_ptn; 
	if (energy_type == 1)
	{
		
		func_ptn = &HaEmpiricalMod::ScoreEnergy;
	}
	else if (energy_type == 2)
	{
		func_ptn = &HaEmpiricalMod::GeometryScoreEnergy;
	}
	else
	{
		func_ptn = &HaEmpiricalMod::HarmonicEnergy;
	}

	double var_prev;
//	var_prev = emp_mod->ScoreEnergy();
	var_prev = (emp_mod ->* func_ptn)();

	HaVec_double d_ene1, d_ene2;
	
	double incr = 0;
	int nmol = ptmol.size();
	
	Vec3D trans_v;

	int n_size = nmol*7;
	HaMat_double hessian;
	hessian.newsize(n_size,n_size);
	int ind = 0;
	
	int i,j;
	int ind1;
	ind=0;
	HaVec_double d_var; // Incriment in rotation and translation, First 4 for rotation
	d_var.newsize(n_size);
	int imol;
	for(imol= 0; imol < nmol; imol++)
	{
		d_var[ind++]= 0.001;
		d_var[ind++]= 0.001;
		d_var[ind++]= 0.001;
		d_var[ind++]= 0.001;
		d_var[ind++]= 0.001;
		d_var[ind++]= 0.001;
		d_var[ind++]= 0.001;
	}
	
	HaVec_double coord_old;
	coord_old.newsize(n_size);
	ind=0;
	int cos_flag = 1;

	int temp1 = 0;

	Quaternion quat;
	double qang, qx, qy,qz;
	ind=0;
	for(imol= 0; imol < nmol; imol++)
	{
		pMol = (HaMolecule*) ptmol[imol];
		pMol->GetQuaternionTrans(quat, trans_v);
		quat.GetQuaternion(qang, qx, qy,qz);
		coord_old[ind++]= qang;
		coord_old[ind++]= qx;
		coord_old[ind++]= qy;
		coord_old[ind++]= qz;
		coord_old[ind++]= trans_v[0];
		coord_old[ind++]= trans_v[1];
		coord_old[ind++]= trans_v[2];
	}
	
	//Diagonal elements
	HaVec_double coord_plus;
	coord_plus.newsize(n_size);
	HaVec_double coord_minus;
	coord_minus.newsize(n_size);
	HaVec_double coord_2plus;
	coord_2plus.newsize(n_size);
	HaVec_double coord_2minus;
	coord_2minus.newsize(n_size);
	//Off-diagonal elements;
	HaVec_double coord_plus_plus;
	coord_plus_plus.newsize(n_size);
	HaVec_double coord_minus_plus;
	coord_minus_plus.newsize(n_size);
	HaVec_double coord_plus_minus;
	coord_plus_minus.newsize(n_size);
	HaVec_double coord_minus_minus;
	coord_minus_minus.newsize(n_size);
	
	
	double var_plus, var_minus;
	double var_2plus, var_2minus;
	double var_plus_plus, var_minus_plus;
	double var_plus_minus, var_minus_minus;
	
	
	for( i= 0; i < n_size; i++)
	{
		for( j= 0; j < n_size; j++)
		{
			if (j == i) coord_2minus[j] = coord_old[j] - d_var[j]- d_var[j];
			else coord_2minus[j] = coord_old[j];
		}
		
		ind1= 0;
		for(imol= 0; imol < nmol; imol++)
		{
			pMol = (HaMolecule*) ptmol[imol];
			qang= coord_2minus[ind1++];
			qx= coord_2minus[ind1++];
			qy= coord_2minus[ind1++];
			qz= coord_2minus[ind1++];
			quat.SetQuaternion(qang, qx, qy,qz);
			quat.Normalize();
			trans_v[0] = coord_2minus[ind1++];
			trans_v[1] = coord_2minus[ind1++];
			trans_v[2] = coord_2minus[ind1++];
			pMol->SetQuaternionTrans(quat,trans_v);
		}
		var_2minus = (emp_mod ->* func_ptn)();
		
		ind1= 0;
		for(imol= 0; imol < nmol; imol++)
		{
			
			pMol = (HaMolecule*) ptmol[imol];
			qang = coord_old[ind1++];
			qx = coord_old[ind1++];
			qy = coord_old[ind1++];
			qz = coord_old[ind1++];
			quat.SetQuaternion(qang, qx, qy,qz);
			trans_v[0] = coord_old[ind1++];
			trans_v[1] = coord_old[ind1++];
			trans_v[2] = coord_old[ind1++];
			pMol->SetQuaternionTrans(quat,trans_v);
		}
		
		for( j= 0; j < n_size; j++)
		{
			if (j == i) coord_2plus[j] = coord_old[j] + d_var[j]+ d_var[j];
			else coord_2plus[j] = coord_old[j];
		}
		ind1= 0;
		for(imol= 0; imol < nmol; imol++)
		{
			pMol = (HaMolecule*) ptmol[imol];
			qang= coord_2plus[ind1++];
			qx= coord_2plus[ind1++];
			qy= coord_2plus[ind1++];
			qz= coord_2plus[ind1++];
			quat.SetQuaternion(qang, qx, qy,qz);
			quat.Normalize();
			trans_v[0] = coord_2plus[ind1++];
			trans_v[1] = coord_2plus[ind1++];
			trans_v[2] = coord_2plus[ind1++];
			pMol->SetQuaternionTrans(quat,trans_v);
		}
		var_2plus = (emp_mod ->* func_ptn)();
		
		ind1= 0;
		for(imol= 0; imol < nmol; imol++)
		{
			
			pMol = (HaMolecule*) ptmol[imol];
			qang = coord_old[ind1++];
			qx = coord_old[ind1++];
			qy = coord_old[ind1++];
			qz = coord_old[ind1++];
			quat.SetQuaternion(qang, qx, qy,qz);
			trans_v[0] = coord_old[ind1++];
			trans_v[1] = coord_old[ind1++];
			trans_v[2] = coord_old[ind1++];
			pMol->SetQuaternionTrans(quat,trans_v);
		}
		
		for( j= 0; j < n_size; j++)
		{
			if (j == i) 
			{
				coord_plus[j] = coord_old[j] + d_var[j];
			}
			else coord_plus[j] = coord_old[j];
		}
		
		ind1= 0;
		for(imol= 0; imol < nmol; imol++)
		{
			pMol = (HaMolecule*) ptmol[imol];
			qang= coord_plus[ind1++];
			qx= coord_plus[ind1++];
			qy= coord_plus[ind1++];
			qz= coord_plus[ind1++];
			quat.SetQuaternion(qang, qx, qy,qz);
			quat.Normalize();
			trans_v[0] = coord_plus[ind1++];
			trans_v[1] = coord_plus[ind1++];
			trans_v[2] = coord_plus[ind1++];
			pMol->SetQuaternionTrans(quat,trans_v);
		}
		var_plus = (emp_mod ->* func_ptn)();
		
		ind1= 0;
		for(imol= 0; imol < nmol; imol++)
		{
			
			pMol = (HaMolecule*) ptmol[imol];
			qang = coord_old[ind1++];
			qx = coord_old[ind1++];
			qy = coord_old[ind1++];
			qz = coord_old[ind1++];
			quat.SetQuaternion(qang, qx, qy,qz);
			trans_v[0] = coord_old[ind1++];
			trans_v[1] = coord_old[ind1++];
			trans_v[2] = coord_old[ind1++];
			pMol->SetQuaternionTrans(quat,trans_v);
		}
		
		for( j= 0; j < n_size; j++)
		{
			if (j == i) coord_minus[j] = coord_old[j] - d_var[j];
			else coord_minus[j] = coord_old[j];
		}
		
		ind1= 0;
		for(imol= 0; imol < nmol; imol++)
		{
			pMol = (HaMolecule*) ptmol[imol];
			qang= coord_minus[ind1++];
			qx= coord_minus[ind1++];
			qy= coord_minus[ind1++];
			qz= coord_minus[ind1++];
			quat.SetQuaternion(qang, qx, qy,qz);
			quat.Normalize();
			trans_v[0] = coord_minus[ind1++];
			trans_v[1] = coord_minus[ind1++];
			trans_v[2] = coord_minus[ind1++];
			pMol->SetQuaternionTrans(quat,trans_v);
		}
		var_minus = (emp_mod ->* func_ptn)();
		
		ind1= 0;
		for(imol= 0; imol < nmol; imol++)
		{
			
			pMol = (HaMolecule*) ptmol[imol];
			qang = coord_old[ind1++];
			qx = coord_old[ind1++];
			qy = coord_old[ind1++];
			qz = coord_old[ind1++];
			quat.SetQuaternion(qang, qx, qy,qz);
			trans_v[0] = coord_old[ind1++];
			trans_v[1] = coord_old[ind1++];
			trans_v[2] = coord_old[ind1++];
			pMol->SetQuaternionTrans(quat,trans_v);
		}
		
		double numer = - var_2plus +16*var_plus - 30*var_prev +16*var_minus - var_2minus;
		double denum = 12 *  d_var[i] *d_var[i];
		hessian.SetVal_idx0(i,i,numer/denum);
		//		PrintLog("numer/denum %f \n",numer/denum);
		int k;
		for( k= 0; k < n_size; k++)
		{
			if(k>i)
			{
				for( j= 0; j < n_size; j++)
				{
					coord_plus_plus[j] = coord_old[j];
					if (j == i) coord_plus_plus[j] = coord_old[j] + d_var[j];
					if (j == k) coord_plus_plus[j] = coord_old[j] + d_var[j];
				}
				ind1= 0;
				for(imol= 0; imol < nmol; imol++)
				{
					pMol = (HaMolecule*) ptmol[imol];
					qang= coord_plus_plus[ind1++];
					qx= coord_plus_plus[ind1++];
					qy= coord_plus_plus[ind1++];
					qz= coord_plus_plus[ind1++];
					quat.SetQuaternion(qang, qx, qy,qz);
					quat.Normalize();
					trans_v[0] = coord_plus_plus[ind1++];
					trans_v[1] = coord_plus_plus[ind1++];
					trans_v[2] = coord_plus_plus[ind1++];
					pMol->SetQuaternionTrans(quat,trans_v);
				}
				var_plus_plus = (emp_mod ->* func_ptn)();
				
				ind1= 0;
				for(imol= 0; imol < nmol; imol++)
				{
					
					pMol = (HaMolecule*) ptmol[imol];
					qang = coord_old[ind1++];
					qx = coord_old[ind1++];
					qy = coord_old[ind1++];
					qz = coord_old[ind1++];
					quat.SetQuaternion(qang, qx, qy,qz);
					trans_v[0] = coord_old[ind1++];
					trans_v[1] = coord_old[ind1++];
					trans_v[2] = coord_old[ind1++];
					pMol->SetQuaternionTrans(quat,trans_v);
				}
				
				for( j= 0; j < n_size; j++)
				{
					coord_plus_minus[j] = coord_old[j];
					if (j == i) coord_plus_minus[j] = coord_old[j] + d_var[j];
					if (j == k) coord_plus_minus[j] = coord_old[j] - d_var[j];
				}
				ind1= 0;
				for(imol= 0; imol < nmol; imol++)
				{
					pMol = (HaMolecule*) ptmol[imol];
					qang= coord_plus_minus[ind1++];
					qx= coord_plus_minus[ind1++];
					qy= coord_plus_minus[ind1++];
					qz= coord_plus_minus[ind1++];
					quat.SetQuaternion(qang, qx, qy,qz);
					quat.Normalize();
					trans_v[0] = coord_plus_minus[ind1++];
					trans_v[1] = coord_plus_minus[ind1++];
					trans_v[2] = coord_plus_minus[ind1++];
					pMol->SetQuaternionTrans(quat,trans_v);
				}
				var_plus_minus = (emp_mod ->* func_ptn)();
				
				ind1= 0;
				for(imol= 0; imol < nmol; imol++)
				{
					
					pMol = (HaMolecule*) ptmol[imol];
					qang = coord_old[ind1++];
					qx = coord_old[ind1++];
					qy = coord_old[ind1++];
					qz = coord_old[ind1++];
					quat.SetQuaternion(qang, qx, qy,qz);
					trans_v[0] = coord_old[ind1++];
					trans_v[1] = coord_old[ind1++];
					trans_v[2] = coord_old[ind1++];
					pMol->SetQuaternionTrans(quat,trans_v);
				}
				
				for( j= 0; j < n_size; j++)
				{
					coord_minus_minus[j] = coord_old[j];
					if (j == i) coord_minus_minus[j] = coord_old[j] - d_var[j];
					if (j == k) coord_minus_minus[j] = coord_old[j] - d_var[j];
					
				}
				ind1= 0;
				for(imol= 0; imol < nmol; imol++)
				{
					pMol = (HaMolecule*) ptmol[imol];
					qang= coord_minus_minus[ind1++];
					qx= coord_minus_minus[ind1++];
					qy= coord_minus_minus[ind1++];
					qz= coord_minus_minus[ind1++];
					quat.SetQuaternion(qang, qx, qy,qz);
					quat.Normalize();
					trans_v[0] = coord_minus_minus[ind1++];
					trans_v[1] = coord_minus_minus[ind1++];
					trans_v[2] = coord_minus_minus[ind1++];
					pMol->SetQuaternionTrans(quat,trans_v);
				}
				var_minus_minus = (emp_mod ->* func_ptn)();
				
				ind1= 0;
				for(imol= 0; imol < nmol; imol++)
				{
					
					pMol = (HaMolecule*) ptmol[imol];
					qang = coord_old[ind1++];
					qx = coord_old[ind1++];
					qy = coord_old[ind1++];
					qz = coord_old[ind1++];
					quat.SetQuaternion(qang, qx, qy,qz);
					trans_v[0] = coord_old[ind1++];
					trans_v[1] = coord_old[ind1++];
					trans_v[2] = coord_old[ind1++];
					pMol->SetQuaternionTrans(quat,trans_v);
				}
				
				for( j= 0; j < n_size; j++)
				{
					coord_minus_plus[j] = coord_old[j];
					if (j == i) coord_minus_plus[j] = coord_old[j] - d_var[j];
					if (j == k) coord_minus_plus[j] = coord_old[j] + d_var[j];
				}
				ind1= 0;
				for(imol= 0; imol < nmol; imol++)
				{
					pMol = (HaMolecule*) ptmol[imol];
					qang= coord_minus_plus[ind1++];
					qx= coord_minus_plus[ind1++];
					qy= coord_minus_plus[ind1++];
					qz= coord_minus_plus[ind1++];
					quat.SetQuaternion(qang, qx, qy,qz);
					quat.Normalize();
					trans_v[0] = coord_minus_plus[ind1++];
					trans_v[1] = coord_minus_plus[ind1++];
					trans_v[2] = coord_minus_plus[ind1++];
					pMol->SetQuaternionTrans(quat,trans_v);
				}
				var_minus_plus = (emp_mod ->* func_ptn)();
				
				ind1= 0;
				for(imol= 0; imol < nmol; imol++)
				{
					
					pMol = (HaMolecule*) ptmol[imol];
					qang = coord_old[ind1++];
					qx = coord_old[ind1++];
					qy = coord_old[ind1++];
					qz = coord_old[ind1++];
					quat.SetQuaternion(qang, qx, qy,qz);
					trans_v[0] = coord_old[ind1++];
					trans_v[1] = coord_old[ind1++];
					trans_v[2] = coord_old[ind1++];
					pMol->SetQuaternionTrans(quat,trans_v);
				}
				double numer = var_plus_plus - var_plus_minus - var_minus_plus + var_minus_minus;
				double denum = 4 *  d_var[i] *d_var[k];
				hessian.SetVal_idx0(i,k,numer/denum);
				//				PrintLog("numer/denumIK i%d k%d %f \n", i, k, numer/denum);
			}
		}
	}
	for( i= 0; i < n_size; i++)
	{
		for( j= 0; j< n_size; j++)
		{
			double temp = hessian.GetVal_idx0(i,j);
			hessian.SetVal_idx0(j,i,temp);
		}
	}
	
//	for( i= 0; i < nmol*6; i++)
//	{
//		for( j= 0; j< nmol*6; j++)
//		{
//			PrintLog("%2.2f ",hessian.GetVal_idx0(i,j) );
//		}
//		PrintLog("\n");
//	}
	eigen.close();
	check.close();
	return hessian;
}

HaVec_double
HaInterMolMod::Jacobian(int energy_type, VecPtr ptmol)
{
	if(module_to_init_flag)
	{
		Initialize();
	}
	MolSet* pmset = GetMolSet();
	HaMolecule* pMol;
	HaEmpiricalMod* emp_mod = pmset->GetEmpiricalMod(true);
	if(emp_mod -> module_to_init_flag)
	{
		emp_mod ->Initialize();
	}

	typedef double (HaEmpiricalMod::*EnergyFuncPtn)();
	EnergyFuncPtn func_ptn; 
	if (energy_type == 1)
	{
		
		func_ptn = &HaEmpiricalMod::ScoreEnergy;
	}
	else if (energy_type == 2)
	{
		func_ptn = &HaEmpiricalMod::GeometryScoreEnergy;
	}
	else
	{
		func_ptn = &HaEmpiricalMod::HarmonicEnergy;
	}

	double incr = 0;
	int nmol = ptmol.size();
	
	Vec3D trans_v;
	
	int n_size = nmol*7;
	HaVec_double jacobian;
	jacobian.newsize(n_size);
    int ind = 0;
	int i,j;
	int ind1;
	ind=0;
	HaVec_double d_var; // Incriment in rotation and translation, First 4 for rotation
	d_var.newsize(n_size);
	int imol;
	for(imol= 0; imol < nmol; imol++)
	{
		d_var[ind++]= 0.001;
		d_var[ind++]= 0.001;
		d_var[ind++]= 0.001;
		d_var[ind++]= 0.001;
		d_var[ind++]= 0.001;
		d_var[ind++]= 0.001;
		d_var[ind++]= 0.001;
	}

	HaVec_double phi_o(nmol), cos_theta_o(nmol), psi_o(nmol); 
	
	HaVec_double coord_old;
	coord_old.newsize(n_size);
	Quaternion quat;
	double qang, qx, qy,qz;
	ind=0;
	for(imol= 0; imol < nmol; imol++)
	{
		pMol = (HaMolecule*) ptmol[imol];
		//pMol->GetPosEulerTrans( cur_phi, cur_cos_theta, cur_psi,trans_v);
		pMol->GetQuaternionTrans(quat, trans_v);
		quat.GetQuaternion(qang, qx, qy,qz);
		coord_old[ind++]= qang;
		coord_old[ind++]= qx;
		coord_old[ind++]= qy;
		coord_old[ind++]= qz;
		coord_old[ind++]= trans_v[0];
		coord_old[ind++]= trans_v[1];
		coord_old[ind++]= trans_v[2];
	}
	
	
	//Diagonal elements
	HaVec_double coord_plus;
	coord_plus.newsize(n_size);
	HaVec_double coord_minus;
	coord_minus.newsize(n_size);

	double var_plus, var_minus;
	
	for( i= 0; i < n_size; i++)
	{
		for( j= 0; j < n_size; j++)
		{
			if (j == i) 
			{
				coord_plus[j] = coord_old[j] + d_var[j];
			}
			else coord_plus[j] = coord_old[j];
		}
		ind1= 0;
		for(imol= 0; imol < nmol; imol++)
		{
			pMol = (HaMolecule*) ptmol[imol];
			qang= coord_plus[ind1++];
			qx= coord_plus[ind1++];
			qy= coord_plus[ind1++];
			qz= coord_plus[ind1++];
			quat.SetQuaternion(qang, qx, qy,qz);
			quat.Normalize();
			trans_v[0] = coord_plus[ind1++];
			trans_v[1] = coord_plus[ind1++];
			trans_v[2] = coord_plus[ind1++];
			pMol->SetQuaternionTrans(quat,trans_v);
		}
		var_plus = (emp_mod ->* func_ptn)();
		ind1= 0;
		for(imol= 0; imol < nmol; imol++)
		{
			
			pMol = (HaMolecule*) ptmol[imol];
			qang = coord_old[ind1++];
			qx = coord_old[ind1++];
			qy = coord_old[ind1++];
			qz = coord_old[ind1++];
			quat.SetQuaternion(qang, qx, qy,qz);
			trans_v[0] = coord_old[ind1++];
			trans_v[1] = coord_old[ind1++];
			trans_v[2] = coord_old[ind1++];
			pMol->SetQuaternionTrans(quat,trans_v);
		}
		
		for( j= 0; j < n_size; j++)
		{
			if (j == i) coord_minus[j] = coord_old[j] - d_var[j];
			else coord_minus[j] = coord_old[j];
		}
		ind1= 0;
		for(imol= 0; imol < nmol; imol++)
		{
			pMol = (HaMolecule*) ptmol[imol];
			qang= coord_minus[ind1++];
			qx= coord_minus[ind1++];
			qy= coord_minus[ind1++];
			qz= coord_minus[ind1++];
			quat.SetQuaternion(qang, qx, qy,qz);
			quat.Normalize();
			trans_v[0] = coord_minus[ind1++];
			trans_v[1] = coord_minus[ind1++];
			trans_v[2] = coord_minus[ind1++];
			pMol->SetQuaternionTrans(quat,trans_v);
		}
		var_minus = (emp_mod ->* func_ptn)();
		
		ind1= 0;
		for(imol= 0; imol < nmol; imol++)
		{
			
			pMol = (HaMolecule*) ptmol[imol];
			qang = coord_old[ind1++];
			qx = coord_old[ind1++];
			qy = coord_old[ind1++];
			qz = coord_old[ind1++];
			quat.SetQuaternion(qang, qx, qy,qz);
			trans_v[0] = coord_old[ind1++];
			trans_v[1] = coord_old[ind1++];
			trans_v[2] = coord_old[ind1++];
			pMol->SetQuaternionTrans(quat,trans_v);
		}
		
		double numer1 = var_plus - var_minus;
		double denum1 = 2 * d_var[i];
		jacobian[i] = numer1/denum1;
	}

	return jacobian;
}

int InterMolEnergyMinimizer::MinimizeEnergy(int energy_type, VecPtr ptmol)
{
	std::fstream enefile;
	enefile.open("minimizedenergy.dat", std::ios::out| std::ios::app );
	
	std::fstream stepfile;
	stepfile.open("step_vector.dat", std::ios::out| std::ios::app );
	
	MolSet* pmset = p_inter_mol->GetMolSet();
	HaEmpiricalMod* emp_mod = pmset->GetEmpiricalMod(true);
	typedef double (HaEmpiricalMod::*EnergyFuncPtn)();
	EnergyFuncPtn func_ptn; 
	if (energy_type == 1)
	{
		
		func_ptn = &HaEmpiricalMod::ScoreEnergy;
	}
	else if (energy_type == 2)
	{
		func_ptn = &HaEmpiricalMod::GeometryScoreEnergy;
	}
		else
	{
		func_ptn = &HaEmpiricalMod::HarmonicEnergy;
	}

	//
	//
	//	double hhh = fp();
	//	PrintLog("POINTER ENERGY %2.2f \n", hhh);
	//	HaInterMolMod* mm = pmset->GetInterMolMod(true);
	//	double nnn = (emp_mod ->* FuncPtn)();
	//	cout<<"Wasup "<< nnn <<"\n";
	//	cout<<"Wasup "<<foo2(&fsq,11)<<"\n";
	//

	int nmol = ptmol.size();
	int i,j;
	double n_size = nmol*7;	
	HaMat_double hessian_minus_one(n_size,n_size);
	HaMat_double hessian(n_size,n_size);
	HaVec_double jacobian(n_size);
	NumVector<double> step_vector(n_size);
	double stepsize_tr = 1;
	double stepsize_rot = 1;

	Vec3D trans_v;
	char buf[256];
	double ene_prev, ene_cur;
	
	int opt_flag = FALSE;
	int istep;
	int result;
	double tollerance = 0.01;
	double min_ene_threshold = 0.5;
	int minimizer_steps = 100;

	ene_cur  = (emp_mod ->* func_ptn)();
	ene_prev = 10e16;

	if (ene_cur > min_ene_threshold) 
	{
		PrintLog("Calculating Hessian first ... ... ...\n");
		hessian = p_inter_mol->Hessian(energy_type, ptmol);
		PrintLog("Calculating Hessian first ... ... ... done \n");
		jacobian = p_inter_mol->Jacobian(energy_type, ptmol);
	}
	
	double f = 0.0 ;
	double delta_ene = 0;
	ene_cur  = (emp_mod ->* func_ptn)();
	PrintLog("Initial Energy is %14.6f kcal/mol \n", ene_cur);
	
	for (istep = 0; istep < minimizer_steps; istep++)
//	for (; ;)
	{
		if (istep >0) ene_cur = f;
		if (ene_cur <= min_ene_threshold)
		{
			pmset->info_str.pop_back();
			sprintf(buf,"Step %d Minimized energy is %14.6f kcal/mol ", istep, ene_cur);
			PrintMessage(buf);
			pmset->info_str.push_back(buf);
			return TRUE;
		}
		else
		{
			sprintf(buf,"cur= %4.6f ene_prev= %4.6f", ene_cur, ene_prev);
			enefile << buf << std::endl ;
			if (istep >0)
			{
				delta_ene = fabs(ene_cur - ene_prev);
				if (result == 0) // Sufficient function decrease.
				{
//					if ( delta_ene < tollerance && opt_flag) break;
					if ( delta_ene < tollerance)
					{
						PrintLog("Calculating Hessian, Sufficient function decrease ...\n");
						hessian = p_inter_mol->Hessian(energy_type, ptmol);
						opt_flag = TRUE;
					}
					jacobian = p_inter_mol->Jacobian(energy_type, ptmol);
				}
				else if (result == 1) // Convergence (in search of ZERO energy this result requeres checking if it is really ZERO).
				{
					PrintLog("result %d delta_ene %2.2f tollerance %2.2f\n",result, delta_ene, tollerance);
					if ( delta_ene < tollerance) break;
					PrintLog("Calculating Hessian, next step, close to Convergence ...\n");
					hessian = p_inter_mol->Hessian(energy_type, ptmol);
					PrintLog("Calculating Jacobian, next step, close to Convergence ...\n");
					jacobian = p_inter_mol->Jacobian(energy_type, ptmol);
					opt_flag = TRUE;

				}
			}	
			ene_prev = ene_cur;
			hessian_minus_one= hessian;
			double eigenval_min = 10e9;
			
			HaMat_double cc(n_size,n_size);
			HaVec_double eigval_hessian(n_size);
			HaMat_double::mat_sdiag(hessian, cc, eigval_hessian);
			
			for( i = 1; i <= n_size; i++)
			{
				if (eigval_hessian(i) < eigenval_min)  eigenval_min = eigval_hessian(i);
			}
			if (eigenval_min <0 ) 
			{
				for( i= 0; i < n_size; i++)
				{
					double temp_eig = hessian_minus_one.GetVal_idx0(i,i) + fabs(eigenval_min) + 1.0;
					hessian_minus_one.SetVal_idx0(i,i,temp_eig);
				}
				
			}
			//		HaVec_double eigval_hessian_minus_one(nmol*6);
			//		HaMat_double cc(nmol*6,nmol*6);
			//		HaMat_double::mat_sdiag(hessian_minus_one, cc, eigval_hessian_minus_one);
			//		for( i = 1; i <= nmol*6; i++)
			//		{
			//			PrintLog("eig_val^1(%d)=%3.2f \n",i,eigval_hessian_minus_one(i));
			//		}
			
			for( i= 0; i < n_size; i++)
			{
				for( j= 0; j< n_size; j++)
				{
					sprintf(buf,"%2.3f ", hessian.GetVal_idx0(i,j) );
					stepfile << buf ;
				}
				stepfile << std::endl ;
			}
			stepfile << "H^-1" << std::endl ;
			
			HaMat_double::mat_inverse(hessian_minus_one);
			for( i= 0; i < n_size; i++)
			{
				for( j= 0; j< n_size; j++)
				{
					sprintf(buf,"%2.3f ", hessian_minus_one.GetVal_idx0(i,j) );
					stepfile << buf ;
				}
				stepfile << std::endl ;
			}
			
			double len_vec=0;
			step_vector = matmult(hessian_minus_one, jacobian);
			for( i= 1; i <= n_size; i++)
			{
				step_vector(i)= - step_vector(i);
				//			sprintf(buf,"1step_vector(%d)= %4.6f,sqrt %4.6f", i, step_vector(i),sqrt(len_vec));
				//			stepfile << buf << endl ;
				len_vec += step_vector(i)*step_vector(i);
			}
			len_vec = 1/sqrt(len_vec);
			for( i= 1; i <= n_size; i++)
			{
				step_vector(i) *= len_vec;
				//			PrintLog("step_vector(%d) = %2.3f \n", i, step_vector(i) );
				sprintf(buf,"step_vector(%d)= %4.6f ,sqrt %4.6f", i, step_vector(i),sqrt(len_vec));
				stepfile << buf << std::endl ;
			}
			stepfile << std::endl ;
			double stpmax = 1.0;
			HaVec_double g;
			g.newsize(n_size);
			HaVec_double p;
			p.newsize(n_size);
			for( i= 0; i < n_size; i++)
			{
				g[i] = jacobian[i];
				p[i] = step_vector(i+1);
			}
			std::fstream file;
				file.open("gradient.dat", std::ios::out| std::ios::app );
				for( i= 0; i < n_size; i++)
				{
					sprintf(buf,"My i(%d) g = %2.3f p= %2.3f",i, g[i],p[i]);
					file << buf << std::endl ;
				}
				file.close();
			
			//		int res =0;
			result = LineSearch(energy_type, ptmol, g, p, &f, stpmax);
			
//			PrintLog("^^^ res = %d Energy after LineSearch %2.4f\n", result, f);
//			sprintf(buf,"Minimized%i%s", k,".pdb");
//			pmset->SavePDBFile(buf);
//		k++;
		}
		if (istep%50 == 0)
		{
			pmset->info_str.clear();
			sprintf(buf,"Step %d Minimized energy is %14.6f kcal/mol ", istep, ene_cur);
			PrintMessage(buf);
			pmset->info_str.push_back(buf);
		}
		
	}
	pmset->RefreshAllViews(RFApply | RFRefresh);
	enefile.close();
	stepfile.close();
	return TRUE;
}


int
InterMolEnergyMinimizer::LineSearch(int energy_type, VecPtr ptmol, HaVec_double g, HaVec_double p, double *f_return, double stpmax)
{
//	fstream enefile;
//	enefile.open("gradient.dat", ios::out|ios::app );
//	for( i= 0; i < 24; i++)
//	{
//		sprintf(buf,"LS i(%d) g = %2.3f p= %2.3f",i, g[i],p[i]);
//		enefile << buf << endl ;
//	}
//	enefile.close();
	// Given an n-dimensional point xold[1..n], the value of the function and gradient there, fold
	// and g[1..n], and a direction p[1..n], finds a new point x[1..n] along the direction p from
	// xold where the function func has decreased sufficiently. The new function value is returned
	// in f. stpmax is an input quantity that limits the length of the steps so that you do not try to
	// evaluate the function in regions where it is undefined or subject to overflow. p is usually the
	// Newton direction. The output quantity check is false (0) on a normal exit. It is true (1) when
	// x is too close to xold. In a minimization algorithm, this usually signals convergence and can
	// be ignored. However, in a zero-finding algorithm the calling program should check whether the
	// convergence is spurious. Some difficult problems may require double precision in this routine.
	double alf = 1.0e-4; // Ensures sufficient decrease in function value.
	double tolx = 1.0e-7; // Convergence criterion on .x.
	double f;
	MolSet* pmset = p_inter_mol->GetMolSet();
	HaEmpiricalMod* emp_mod = pmset->GetEmpiricalMod(true);
	HaMolecule* pMol; 
	int nmol = ptmol.size();
	int ind;
	double a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,test,tmplam;
	double fold;
	Vec3D trans_v;
	int check=0;
	int imol;
	int n_size = nmol*7;
	HaVec_double xold;
	xold.newsize(n_size);

	int i;
	Quaternion quat;
	double qang, qx, qy,qz;


	typedef double (HaEmpiricalMod::*EnergyFuncPtn)();
	EnergyFuncPtn func_ptn; 
	if (energy_type == 1)
	{
		
		func_ptn = &HaEmpiricalMod::ScoreEnergy;
	}
	else if (energy_type == 2)
	{
		func_ptn = &HaEmpiricalMod::GeometryScoreEnergy;
	}
	else
	{
		func_ptn = &HaEmpiricalMod::HarmonicEnergy;
	}

	for (sum=0.0, i=1; i<=n_size ;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax) for (i=1;i<=n_size;i++) p[i] *= stpmax/sum;  // Scale if attempted step is too big.
	for (slope=0.0,i=1;i<=n_size;i++) slope += g[i]*p[i];
	if (slope >= 0.0) 
	{
		PrintLog("Roundoff problem in LineSearch\n");
	}
	test=0.0;								// Compute Lambda min.
	ind =0;
	for(  imol= 0; imol < nmol; imol++)
	{
		pMol = (HaMolecule*) ptmol[imol];
		pMol->GetQuaternionTrans(quat, trans_v);
		quat.GetQuaternion(qang, qx, qy,qz);
		xold[ind] = qang;
		ind++;
		xold[ind] = qx; 
		ind++;
		xold[ind] = qy;
		ind++;
		xold[ind] = qz;
		ind++;
		xold[ind] = trans_v[0];
		ind++;
		xold[ind] = trans_v[1];
		ind++;
		xold[ind] = trans_v[2];
		ind++;
	}
	fold = (emp_mod ->* func_ptn)();
//	PrintLog("f OLD Energy %2.2f \n", fold);

	for (i=1; i<=n_size; i++) 
	{
		temp=fabs(p[i])/FMAX( (float)fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}

	alamin=tolx/test;
	alam= 1.0;					// Always try full Newton step first.
	double check_alam = alam;
	int k = 0;

	for (;;)
//	for(k = 0; k<1; k++)
	{						// Start of iteration loop.
//		PrintLog("Minimizer step size %2.6f \n", alam );
		ind = 0;
		for(  imol= 0; imol < nmol; imol++)
		{
			pMol = (HaMolecule*) ptmol[imol];
			qang = xold[ind] + p[ind]*alam;
			ind++;
			qx = xold[ind] + p[ind]*alam;
			ind++;
			qy = xold[ind] + p[ind]*alam;
			ind++;
			qz = xold[ind] + p[ind]*alam;
			ind++;
			quat.SetQuaternion(qang, qx, qy,qz);
			quat.Normalize();
			trans_v[0] = xold[ind] + p[ind]*alam;
			ind++;
			trans_v[1] = xold[ind] + p[ind]*alam;
			ind++;
			trans_v[2] = xold[ind] + p[ind]*alam;
			ind++;
			pMol->SetQuaternionTrans(quat,trans_v);
		}

		f= (emp_mod ->* func_ptn)(); // *f=(*func)(x); //
//		PrintLog("f Energy %2.2f \n", f);
		if (alam < alamin) 
		{  // Convergence on .x. For zero finding,
			// the calling program should
			// verify the convergence.
			check=1;
			*f_return = f;
			PrintLog("Converged \n");
			return check;
		} 
		else if (f <= fold+alf*alam*slope) // Sufficient function decrease.
		{
//			PrintLog("f %2.2f  fold %2.2f , step_size = %2.2f\n", f, fold, alam);
			*f_return = f;
			return check;
		}
		else 
		{									//	Backtrack.
			if (alam == check_alam) tmplam = -slope/(2.0*(f-fold-slope)); // First time.
			else 
			{									// Subsequent backtracks.
				rhs1 = f-fold-alam*slope;
				rhs2=f2-fold-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else
				{
					disc=b*b-3.0*a*slope;
					if (disc < 0.0) tmplam=0.5*alam;
					else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
					else tmplam=-slope/(b+sqrt(disc));
				}
				if (tmplam > 0.5*alam) tmplam=0.5*alam; 
			}
		}
		alam2=alam;
		f2 = f;
		alam=FMAX(tmplam,0.1*alam); 
//		k++;
	}
}

int
InterMolEnergyMinimizer::SteepestDescentMinimizer(int nsteps)
{
	MolSet* pmset = p_inter_mol->GetMolSet();
	HaEmpiricalMod* emp_mod = pmset->GetEmpiricalMod(true);
	HaMolecule* pMol; 
	int nmol = pmset->GetNMol();
	Vec3DValArray torque_array(nmol);
	Vec3DValArray force_array(nmol);
	Vec3DValArray force_cntrl_array(nmol);
	Vec3DValArray force_total_array(nmol);
	Vec3DValArray torque_ang_array(nmol);
	Vec3DValArray torque_total_array(nmol);
	Vec3D trans_v;
	HaMat_double rmat(3,3);
	int imol;
	int n_size = nmol*7;
	HaVec_double xold;
	xold.newsize(n_size);
	
	int i, ind, istep;
	Quaternion quat;
	Quaternion quat_old;
	
	
	double qang, qx, qy,qz;
	Vec3D quat_vec;
	Vec3D trans_vec;
	
	
	typedef double (HaEmpiricalMod::*EnergyFuncPtn)();
	EnergyFuncPtn func_ptn; 
	func_ptn = &HaEmpiricalMod::GeometryScoreEnergy;
	
	double alf = 1.0;
	double ene_new, ene_old;
	
	emp_mod ->CalcRepulForceTorque(force_array, torque_array);
	emp_mod ->CalcForceCentralAttract(force_cntrl_array);
	emp_mod ->CalcPackAngleForceTorque(torque_ang_array);

	for (imol=0; imol< nmol; imol++)
	{
		force_total_array[imol] = force_array[imol] + force_cntrl_array[imol]; 
		torque_total_array[imol] = torque_array[imol] + torque_ang_array[imol];
	}

//	trans_vec =force_array[0]; 
//	quat_vec = torque_array[0];
//	PrintLog("trans_vec %2.3f %2.3f %2.3f\n", trans_vec.GetX(), trans_vec[1],trans_vec[2]);
//	PrintLog("q111uat_vec %2.3f %2.3f %2.3f\n", quat_vec[0], quat_vec[1],quat_vec[2]);

	double tol = 1.0;
	double force_tol = 100.0*n_size;
	double total_force ;
	ene_new = (emp_mod ->* func_ptn)();

	for(istep=0; istep < nsteps; istep++)
	{
		total_force = 0;
		for (imol=0; imol< nmol; imol++)
		{
			trans_vec = force_total_array[imol]; 
			quat_vec = torque_total_array[imol];
			//quat_vec = torque_ang_array[imol];
//            total_force += trans_vec[0] * trans_vec[0] + trans_vec[1] * trans_vec[1] + trans_vec[2] * trans_vec[2];
		    total_force += quat_vec[0] * quat_vec[0] + quat_vec[1] * quat_vec[1] + quat_vec[2] * quat_vec[2];
//			PrintLog("trans_vec %2.3f %2.3f %2.3f\n", trans_vec.GetX(), trans_vec[1],trans_vec[2]);
//			PrintLog("quat_vec %2.3f %2.3f %2.3f\n", quat_vec[0], quat_vec[1],quat_vec[2]);
		}
		PrintLog("total_force <= force_tol %2.3f  %2.3f,  ALFA= %2.3f, step %d\n", total_force, force_tol, alf, istep);

		if (total_force <= force_tol) break;
		ind = 0;
		for(  imol= 0; imol < nmol; imol++)
		{
			pMol = pmset->GetMolByIdx(imol);
			pMol->GetQuaternionTrans(quat,trans_v);
			quat.GetQuaternion(qang, qx, qy,qz);
			xold[ind] = qang;
			ind++;
			xold[ind] = qx; 
			ind++;
			xold[ind] = qy;
			ind++;
			xold[ind] = qz;
			ind++;
			xold[ind] = trans_v[0];
			ind++;
			xold[ind] = trans_v[1];
			ind++;
			xold[ind] = trans_v[2];
			ind++;
		}


		alf = GoldenSectionSearch(xold, force_total_array, torque_total_array, tol);


		StepAlongGradient(xold, alf, force_total_array, torque_total_array);


		ene_old = ene_new;
		emp_mod ->CalcRepulForceTorque(force_array, torque_array);
		emp_mod ->CalcForceCentralAttract(force_cntrl_array);
		emp_mod ->CalcPackAngleForceTorque(torque_ang_array);
		for (imol=0; imol< nmol; imol++)
		{
			force_total_array[imol] = force_array[imol] + force_cntrl_array[imol]; 
			torque_total_array[imol] = torque_array[imol] + torque_ang_array[imol];
		}


		ene_new = (emp_mod ->* func_ptn)();
		tol = fabs( (ene_new - ene_old) / ene_new) * 10.0;
//		PrintLog("istep %d Tolerance %2.1f ene_new %2.1f - ene_old %2.1f\n", istep, tol, ene_new,  ene_old);

	}
	return TRUE;
}
	

double
InterMolEnergyMinimizer::GoldenSectionSearch(HaVec_double& xold, Vec3DValArray& force_array, Vec3DValArray& torque_array, double& tol)
{
	MolSet* pmset = p_inter_mol->GetMolSet();
	HaEmpiricalMod* emp_mod = pmset->GetEmpiricalMod(true);
	double c1, a, b;

	double alpha, alpha2;
	alpha =1.0;
	double f1, f2;

	HaMolecule* pMol; 
	int nmol = pmset->GetNMol();
	Vec3D trans_v;
	double qang, qx, qy, qz;
	int imol, ind;
	Quaternion quat;
	ind= 0;
	typedef double (HaEmpiricalMod::*EnergyFuncPtn)();
	EnergyFuncPtn func_ptn; 
	func_ptn = &HaEmpiricalMod::GeometryScoreEnergy;

	c1 = 1.0;
	a = -c1 - 2.0;
	b= c1 + 1.5;
	alpha = a + (1.0 - c1) * (b - a);
	StepAlongGradient(xold, alpha, force_array, torque_array);
	f1 = (emp_mod ->* func_ptn)();
	ind =0;
	for(  imol= 0; imol < nmol; imol++)
	{
		pMol = pmset->GetMolByIdx(imol);
		qang = xold[ind];
		ind++;
		qx = xold[ind]; 
		ind++;
		qy = xold[ind];
		ind++;
		qz = xold[ind];
		ind++;
		trans_v[0] = xold[ind];
		ind++;
		trans_v[1] = xold[ind];
		ind++;
		trans_v[2] = xold[ind];
		ind++;
		quat.SetQuaternion(qang, qx, qy,qz);
		pMol->SetQuaternionTrans(quat,trans_v);
	}
				
	alpha2 = a + (c1 * (b - a));
	StepAlongGradient(xold, alpha2, force_array, torque_array);
	f2 = (emp_mod ->* func_ptn)();
	ind =0;
	for(  imol= 0; imol < nmol; imol++)
	{
		pMol = pmset->GetMolByIdx(imol);
		qang = xold[ind];
		ind++;
		qx = xold[ind]; 
		ind++;
		qy = xold[ind];
		ind++;
		qz = xold[ind];
		ind++;
		trans_v[0] = xold[ind];
		ind++;
		trans_v[1] = xold[ind];
		ind++;
		trans_v[2] = xold[ind];
		ind++;
		quat.SetQuaternion(qang, qx, qy,qz);
		pMol->SetQuaternionTrans(quat,trans_v);
	}

	int i;
	for(i =0; i< 50; i++) 
	{
		if (fabs(f1 - f2) <= tol) break;
//		PrintLog("IF fabs(f1 - f2)= %2.3f <= tol= %2.3f \n", fabs(f1 - f2), tol);
		if (f1 > f2)
		{
            a = alpha;
			alpha = alpha2;
			f1 = f2;
			alpha2 = a + (c1 * (b - a));
			StepAlongGradient(xold, alpha2, force_array, torque_array);
			f2 = (emp_mod ->* func_ptn)();
			ind =0;
			for(  imol= 0; imol < nmol; imol++)
			{
				pMol = pmset->GetMolByIdx(imol);
				qang = xold[ind];
				ind++;
				qx = xold[ind]; 
				ind++;
				qy = xold[ind];
				ind++;
				qz = xold[ind];
				ind++;
				trans_v[0] = xold[ind];
				ind++;
				trans_v[1] = xold[ind];
				ind++;
				trans_v[2] = xold[ind];
				ind++;
				quat.SetQuaternion(qang, qx, qy,qz);
				pMol->SetQuaternionTrans(quat,trans_v);
			}
		}
		else
		{
            b = alpha2;
            alpha2 = alpha;
            f2 = f1;
            alpha = a + (1 - c1) * (b - a);
			StepAlongGradient(xold, alpha, force_array, torque_array);
			f1 =(emp_mod ->* func_ptn)();
			ind =0;
			for(  imol= 0; imol < nmol; imol++)
			{
				pMol = pmset->GetMolByIdx(imol);
				qang = xold[ind];
				ind++;
				qx = xold[ind]; 
				ind++;
				qy = xold[ind];
				ind++;
				qz = xold[ind];
				ind++;
				trans_v[0] = xold[ind];
				ind++;
				trans_v[1] = xold[ind];
				ind++;
				trans_v[2] = xold[ind];
				ind++;
				quat.SetQuaternion(qang, qx, qy,qz);
				pMol->SetQuaternionTrans(quat,trans_v);
			}
		}
	}


	return alpha;
}

int
InterMolEnergyMinimizer::StepAlongGradient(HaVec_double& xold, double alpha, Vec3DValArray& force_array, Vec3DValArray& torque_array)
{
	int move_flag =0;
	MolSet* pmset = p_inter_mol->GetMolSet();
	HaMolecule* pMol; 
	int nmol = pmset->GetNMol();
	Vec3D trans_v;
	Vec3D trans_vec;
	Vec3D quat_vec;
	double qang, qx, qy, qz;
	double qang_o, qx_o, qy_o, qz_o;
	HaMat_double rmat_old(3,3);

	int imol, ind;
	Quaternion quat;
	Quaternion quat_old;

	double len;
	ind =0;
	for(  imol= 0; imol < nmol; imol++)
	{
		quat_vec = torque_array[imol];
		trans_vec = force_array[imol];
		pMol = pmset->GetMolByIdx(imol);
		len = quat_vec.length();
		qang = len*alpha;
		len = 1.0/len;
		qang_o = xold[ind];
		ind++;
		qx = -quat_vec[0]*len;
		qx_o =  xold[ind];
		ind++;
		qy = -quat_vec[1]*len;
		qy_o = xold[ind];
		ind++;
		qz = -quat_vec[2]*len;
		qz_o = xold[ind];
		ind++;
		quat.SetQuaternion(qang, qx, qy, qz);

		quat_old.SetQuaternion(qang_o, qx_o, qy_o, qz_o);
		quat.Normalize();
		quat.operator *=(quat_old);
		len = 1.0/ trans_vec.length();
		trans_v[0] = xold[ind] - trans_vec[0]*len*alpha;
		ind++;
		trans_v[1] = xold[ind] - trans_vec[1]*len*alpha;
		ind++;
		trans_v[2] = xold[ind] - trans_vec[2]*len*alpha;
		ind++;
		pMol->SetQuaternionTrans(quat,trans_v);
//		pMol->SetQuaternionTrans(quat_old,trans_v);
	}
	return move_flag =1;
}

InterMolMCSimulator::InterMolMCSimulator( HaInterMolMod* p_im_mod_new) 
{ 
	p_im_mod = p_im_mod_new; 
	pmset = p_im_mod_new->GetMolSet();
	p_rand_num_gen = new Random(5);
	p_crd = NULL;

	TrajIOAgent* p_ene_ag = GetTrajectoryIOAgent();
	
	SetStdParams();
}

 InterMolMCSimulator::~InterMolMCSimulator()
 {
	if(p_rand_num_gen) delete p_rand_num_gen;
	if(p_crd != NULL) delete p_crd;
	
	int na = agents.size();
	int i;
	for(i=0; i < na; i++)
	{
		TrajAnalAgent* p_ag = agents[i];
		if( p_ag != NULL) delete p_ag;
	}
 }

void
InterMolMCSimulator::SetStdParams()
{
	freeze_first_mol = true; 
	
	ang_ratio = 0.02;
	tr_ratio  = 0.3;

	equil_conf_vol_vdw = FALSE;
	
	amber_flag = FALSE;
	mc_steps_betw_loc_min = 50;	
	rex_flag = FALSE;
	xy_mc_flag  = FALSE;
}

void
InterMolEnergyMinimizer::SetStdParams()
{
	
}

void
InterMolRepExchSimulator::SetStdParams()
{
	nreplicas = 2;
	ireplica = 0;
	exchange_arr[0] = 0;
	position_mat.newsize(nreplicas,10*7); // Move initialization somewhere else
	energy_arr.newsize(nreplicas);
	exchange_arr.newsize(nreplicas);
	p_acc_ratio.newsize(nreplicas); 
	rem_steps = 10;
	temperature_max = 300.0;
	MC_traj_file_replica_basename = "MC_traj_file_temperature";
	MC_energy_file_replica_basename = "MC_ene_temperature";
	MC_rst_file_basename = "restart_MC";
	n_playback_replica = 0;
	vary_temperature_flag = TRUE;
}

