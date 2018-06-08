/*!
     Classes to study protonation and redox equilibrium

    \author Igor Kurnikov 
    \date 2010-
*/

#define PROTONREDOX_CPP

#include "mpi.h"

#include <bitset>

#include <float.h>

#include <boost/algorithm/string.hpp>

#include "tinyxml.h"
#include "haxml.h"
#include "haobject.h"
#include "randomip.h"
#include "hamatdb.h"
#include "haatgroup.h"
#include "hamolecule.h"
#include "hamolset.h"
#include "moleditor.h"
#include "elmod.h"
#include "protonredox.h"

IntStrMap SetAltChemStateTypeLbls()
{
	IntStrMap lmap;
	lmap[AltChemStateType::PROTONATED]   = "PROTONATED";
	lmap[AltChemStateType::UNPROTONATED] = "UNPROTONATED";
	lmap[AltChemStateType::REDUCED]      = "REDUCED";
	lmap[AltChemStateType::OXIDIZED]     = "OXIDIZED";
	return lmap;
}

IntStrMap AltChemStateType::labels = SetAltChemStateTypeLbls();  

AltChemStateType::AltChemStateType()
{
	v_= PROTONATED;	
}
	
AltChemStateType::~AltChemStateType()
{

}

int AltChemStateType::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int AltChemStateType::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}


int ProtonRedoxMod::heme_model = 0;

AltChemState::AltChemState()
{
	SetStdParam();
}

AltChemState::AltChemState(AtomGroup* new_host_atom_group)
{
	SetStdParam();
	host_atom_group = new_host_atom_group;
}

AltChemState::AltChemState(const AltChemState& ref_state)
{
	id = ref_state.id;
	chmap = ref_state.chmap;
    mod_atom_name = ref_state.mod_atom_name;
	alt_state_type = ref_state.alt_state_type;
	pk = ref_state.pk;
	std_pk = ref_state.std_pk;
	active_flag = ref_state.active_flag;
	host_atom_group = ref_state.host_atom_group;
}

AltChemState::~AltChemState()
{

}

void AltChemState::SetStdParam()
{
	id = "GEN_ALT_STATE";
	mod_atom_name = "";
	alt_state_type = alt_state_type.PROTONATED;
	pk = 0.0;
	chmap.clear();
	std_pk = 0.0;
    active_flag = TRUE;
	host_atom_group = NULL;	
}

int AltChemState::SetAltChForAtom(const char* at_name, double new_ch)
{
    AtomIteratorAtomGroup aitr(host_atom_group);
	HaAtom* aptr = host_atom_group->GetAtomByName(at_name);
	if( aptr != NULL)
	{
       chmap[aptr] = new_ch;
	   return TRUE;
	}
    return FALSE;
}

int AltChemState::SetAltCharges(double weight)
//! if weight = -1.0 set charges equal to the difference between standard
//! and alternative charges
//! 
{
	PtrDoubleMap::iterator itr;

	for( itr = chmap.begin(); itr != chmap.end(); itr++)
	{
		double   delt_ch  = itr.GetVal();
		HaAtom* aptr = (HaAtom*) itr.GetKey();
		HaResidue* pres = aptr->GetHostRes();
		
        HaResidue* res_templ = pres->GetTemplate();
        if(res_templ == NULL)
		{
            PrintLog(" Error in AltChemState::SetAltCharges() \n");
			PrintLog(" Residue %s does not have a template \n", (pres->GetRef()).c_str());
			continue;
		}

		HaAtom* aptr_templ = res_templ->GetAtomByName(aptr->GetName());
		if(aptr_templ == NULL) 
		{
            PrintLog(" Error in AltChemState::SetAltCharges() \n");
			PrintLog(" Atom %s does not have a conterpart in residue template \n", (aptr->GetRef()).c_str());
			continue;	
		}
		
		double ch_std = aptr_templ->GetCharge();
		double ch_alt = ch_std + delt_ch;
		double ch_new;
		if( fabs(weight + 1.0) < 0.000001)
		{
			ch_new = delt_ch;
		}
		else
		{
			ch_new = ch_std*(1.0 - weight) + ch_alt*weight;
		}
		aptr->SetCharge( ch_new );
	}

	return TRUE;
}

IntStrMap SetMultiSitePopulationMethodLbls()
{
	IntStrMap lmap;
	lmap[MultiSitePopulationMethod::SCF_MULTI_SITE_CALC] = "SCF Multi Site Calc";
	lmap[MultiSitePopulationMethod::MC_MULTI_SITE_CALC]  = "Monte Carlo Multi Site Calc";
	return lmap;
}

IntStrMap MultiSitePopulationMethod::labels = SetMultiSitePopulationMethodLbls();  

MultiSitePopulationMethod::MultiSitePopulationMethod()
{
	v_= SCF_MULTI_SITE_CALC;	
}
	
MultiSitePopulationMethod::~MultiSitePopulationMethod()
{

}

int MultiSitePopulationMethod::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int MultiSitePopulationMethod::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}


ProtonRedoxMod::ProtonRedoxMod(HaMolSet* new_phost_mset):
HaCompMod(COMP_MOD_PROTON_REDOX, new_phost_mset)
{
	set_std_redox_pot = FALSE;
	save_alt_st_inter = TRUE;
	read_alt_st_inter = FALSE;
	save_titration_data = FALSE;
	save_only_redox_titr = FALSE;

    n_mc_cyc = 10000;
	multi_site_pop_method = multi_site_pop_method.SCF_MULTI_SITE_CALC;

	ph_min = 4.0;
    ph_max = 10.0;
	ph_step = 0.5;

	e0_min  = -0.3;
	e0_max  = +0.3;
	e0_step = +0.05;

	ph = 7.0;
	e0 = 0.0;
}

ProtonRedoxMod::~ProtonRedoxMod()
{

}

void ProtonRedoxMod::ClearResAltChemStates(HaResidue* pres)
{
	vector<AltChemState*>::iterator sitr = alt_chem_states.begin();
	for(; sitr != alt_chem_states.end(); )
	{
		AltChemState* p_alt_st = *sitr;
		if( p_alt_st->GetHostAtomGroup() == (AtomGroup*) pres )
		{
			sitr = alt_chem_states.erase(sitr);
		}
		else
		{
			sitr++;
		}
	}
}
	
void ProtonRedoxMod::ClearAltChemStates()
{
	int n = alt_chem_states.size();
	int i;
	for(i = 0; i < n; i++)
	{
		delete alt_chem_states[i];
	}
	alt_chem_states.clear();
}

int ProtonRedoxMod::GetNumResAltChemStates(HaResidue* pres)
{
	return res_altst_map.count((AtomGroup*)pres);
}

AltChemState* ProtonRedoxMod::GetResAltChemState(HaResidue* pres, int alt_state_idx)
{
	int ns = GetNumResAltChemStates(pres);
	if( alt_state_idx < 0 || alt_state_idx >= ns ) return NULL; 

	pair<altst_map_type::iterator,altst_map_type::iterator> eq_range_pair;
	eq_range_pair = res_altst_map.equal_range(pres);
	altst_map_type::iterator sitr;

	int idx = 0;
	for( sitr = eq_range_pair.first; sitr != eq_range_pair.second; sitr++)
	{
		AltChemState* p_alt_chem_st = (*sitr).second;
		if( idx == alt_state_idx)  return p_alt_chem_st;
		idx++;
	}
	return NULL;
}

AltChemState* ProtonRedoxMod::GetResAltChemStateByAtName(HaResidue* pres, const char* at_name)
{
	int ns = GetNumResAltChemStates(pres);
	if( ns == 0) return NULL;

	pair<altst_map_type::iterator,altst_map_type::iterator> eq_range_pair;
	eq_range_pair = res_altst_map.equal_range(pres);
	altst_map_type::iterator sitr;

	for( sitr = eq_range_pair.first; sitr != eq_range_pair.second; sitr++)
	{
		AltChemState* p_alt_chem_st = (*sitr).second;
		if(  strcmp(p_alt_chem_st->mod_atom_name.c_str(),at_name) == 0 ) return p_alt_chem_st;
	}
	return NULL;
}

void ProtonRedoxMod::PrintResPKa(HaResidue* pres)
{
	int ist;
	int nst = GetNumResAltChemStates(pres);
	std::string res_ref = pres->GetRef();
	if(nst == 0) 
	{
		PrintLog("Residue %s do not have alternative protonation states \n",res_ref.c_str());
		return;
	}
	std::string res_full_name = pres->GetFullName();

	PrintLog("pKas for %s \n",res_ref.c_str());

	for(ist = 0; ist < nst; ist++)
	{
		const AltChemState* alt_res_st = GetResAltChemState(pres,ist);
		if(alt_res_st == NULL) continue;
		PrintLog(" ( %s -> %s ) pK = %12.6f \n", res_full_name.c_str(), alt_res_st->id.c_str(),alt_res_st->pk);
	}
}

AltChemState* ProtonRedoxMod::AddAltChemState(AtomGroup* pgrp)
{
	if( pgrp == NULL) return NULL;
	AltChemState* p_alt_chem_st = new AltChemState(pgrp);
	alt_chem_states.push_back( p_alt_chem_st );
	altst_map_type::value_type pair(pgrp,p_alt_chem_st);
	res_altst_map.insert(pair);
	return p_alt_chem_st;
}


int ProtonRedoxMod::SetStdResPKa_G1(HaResidue* pres, int set_redox_pot)
{
	std::string res_full_name = pres->GetFullName();
	std::string res_sh_name = pres->GetName();
	std::string name_mod    = pres->GetNameModifier();

	ClearResAltChemStates(pres);

	if(res_sh_name == "HIS")
	{
		AltChemState* p_alt_chem_st = AddAltChemState(pres);
		if( name_mod == "PROT_TO_HID")
		{
			p_alt_chem_st->id       = "PROT_TO_HID->PROT";
			p_alt_chem_st->SetAltChForAtom("N",0.0678);
			p_alt_chem_st->SetAltChForAtom("H",0.0028);
			p_alt_chem_st->SetAltChForAtom("CA",-0.1542);
			p_alt_chem_st->SetAltChForAtom("HA",0.0331);
			p_alt_chem_st->SetAltChForAtom("CB",0.0048);
			p_alt_chem_st->SetAltChForAtom("HB2",0.0408);
			p_alt_chem_st->SetAltChForAtom("HB3",0.0408);
			p_alt_chem_st->SetAltChForAtom("CG",0.0254);
			p_alt_chem_st->SetAltChForAtom("ND1",0.2298);
			p_alt_chem_st->SetAltChForAtom("HD1",0.0217);
			p_alt_chem_st->SetAltChForAtom("CE1",-0.2227);
			p_alt_chem_st->SetAltChForAtom("HE1",0.1289);
			p_alt_chem_st->SetAltChForAtom("NE2",0.4009);
			p_alt_chem_st->SetAltChForAtom("HE2",0.3911);
			p_alt_chem_st->SetAltChForAtom("CD2",-0.2433);
			p_alt_chem_st->SetAltChForAtom("HD2",0.117);
			p_alt_chem_st->SetAltChForAtom("C",0.1368);
			p_alt_chem_st->SetAltChForAtom("O",-0.0215);
			p_alt_chem_st->mod_atom_name = "NE2";
			p_alt_chem_st->std_pk = 6.3;
		}
		else if( name_mod == "PROT_TO_HIE")
		{
			p_alt_chem_st->id       = "PROT_TO_HID->PROT";
			p_alt_chem_st->SetAltChForAtom("N",0.0678);
			p_alt_chem_st->SetAltChForAtom("H",0.0028);
			p_alt_chem_st->SetAltChForAtom("CA",-0.0773);
			p_alt_chem_st->SetAltChForAtom("HA",-0.0148);
			p_alt_chem_st->SetAltChForAtom("CB",-0.0340);
			p_alt_chem_st->SetAltChForAtom("HB2",0.0443);
			p_alt_chem_st->SetAltChForAtom("HB3",0.0443);
			p_alt_chem_st->SetAltChForAtom("CG",-0.1880);
			p_alt_chem_st->SetAltChForAtom("ND1",0.3919);
			p_alt_chem_st->SetAltChForAtom("HD1",0.3866);
			p_alt_chem_st->SetAltChForAtom("CE1",-0.1805);
			p_alt_chem_st->SetAltChForAtom("HE1",0.1246);
			p_alt_chem_st->SetAltChForAtom("NE2",0.1077);
			p_alt_chem_st->SetAltChForAtom("HE2",0.0572);
			p_alt_chem_st->SetAltChForAtom("CD2",0.1066);
			p_alt_chem_st->SetAltChForAtom("HD2",0.0455);
			p_alt_chem_st->SetAltChForAtom("C",0.1368);
			p_alt_chem_st->SetAltChForAtom("O",-0.0215);
			p_alt_chem_st->mod_atom_name = "ND1";
			p_alt_chem_st->std_pk = 6.5;
		}
		else if( name_mod == "EPSILON")
		{
			p_alt_chem_st->id       = "HIS->ND1_CH_PLUS_1";
			p_alt_chem_st->SetAltChForAtom("ND1",1.0);
			p_alt_chem_st->mod_atom_name = "ND1";
			p_alt_chem_st->std_pk = 6.5;
		}
		else
		{
			p_alt_chem_st->id       = "HIS->NE2_CH_PLUS_1";
			p_alt_chem_st->SetAltChForAtom("NE2",1.0);
			p_alt_chem_st->mod_atom_name = "NE2";
			p_alt_chem_st->std_pk = 6.3;
		} 
        p_alt_chem_st->alt_state_type =  p_alt_chem_st->alt_state_type.PROTONATED;
		p_alt_chem_st->pk = p_alt_chem_st->std_pk;
        
	}
	
	if(res_sh_name == "HEM")
	{
		{
			AltChemState* p_alt_chem_st = AddAltChemState(pres);
			p_alt_chem_st->id       = "HEM_PROT_1";
			p_alt_chem_st->SetAltChForAtom("CGA",0.2);	
			p_alt_chem_st->SetAltChForAtom("O1A",0.4);	
			p_alt_chem_st->SetAltChForAtom("O2A",0.4);	
			p_alt_chem_st->alt_state_type =  p_alt_chem_st->alt_state_type.PROTONATED;
			p_alt_chem_st->std_pk = 4.5;
			p_alt_chem_st->pk = p_alt_chem_st->std_pk;
			p_alt_chem_st->mod_atom_name = "O1A";
		}
		{
			AltChemState* p_alt_chem_st = AddAltChemState(pres);
			p_alt_chem_st->id       = "HEM_PROT_2";
			p_alt_chem_st->SetAltChForAtom("CGD",0.2);	
			p_alt_chem_st->SetAltChForAtom("O1D",0.4);	
			p_alt_chem_st->SetAltChForAtom("O2D",0.4);	
			p_alt_chem_st->alt_state_type =  p_alt_chem_st->alt_state_type.PROTONATED;
			p_alt_chem_st->std_pk = 4.5;
			p_alt_chem_st->pk = p_alt_chem_st->std_pk;
			p_alt_chem_st->mod_atom_name = "O1D";
		}
		if(set_redox_pot)
		{
			char buf[128];
			sprintf(buf,"%d",pres->GetSerNo());
			std::string serno_str = buf;
			boost::trim(serno_str);
			std::string grp_name = "HEM" + serno_str;
			std::string ch_str; 
			HaChain* chain = pres->GetHostChain();
			ch_str += chain->ident;

			if( ch_str != " ") grp_name += ch_str;
	        grp_name += "_HIS";
			HaMolecule* pmol = pres->GetHostMol();
			grp_name = ((std::string)pmol->GetObjName()) + "_" + grp_name; 
            HaMolSet* pmset = pmol->GetHostMolSet();
			AtomGroup* atgrp = pmset->GetAtomGroupByID( grp_name.c_str() );
			HaAtom* aptr;

            if( atgrp == NULL)
			{
				atgrp = pmset->AddAtomGroup(grp_name.c_str());
				AtomIteratorAtomGroup aitr(pres);
                for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
				{
                   atgrp->InsertAtom(aptr);
				}

			    HaAtom* p_fe = pres->GetAtomByName("FE");
			    AtomGroup bonded_atoms;
			    p_fe->GetBondedAtoms(bonded_atoms);

			    AtomIteratorAtomGroup aitr_b(&bonded_atoms);
			    HaAtom* aptr_b = aitr_b.GetFirstAtom();
			    for( aptr_b = aitr_b.GetFirstAtom(); aptr_b; aptr_b = aitr_b.GetNextAtom())
				{
                   HaResidue* pres_b = aptr_b->GetHostRes();
				   if( ((std::string) pres_b->GetName()) == (std::string) "HIS")
				   { 
					  HaAtom* aptr;
					  aptr = pres_b->GetAtomByName("NE2"); if(aptr)atgrp->InsertAtom(aptr);
					  aptr = pres_b->GetAtomByName("CG");  if(aptr)atgrp->InsertAtom(aptr);
					  aptr = pres_b->GetAtomByName("ND1"); if(aptr)atgrp->InsertAtom(aptr);
					  aptr = pres_b->GetAtomByName("CE1"); if(aptr)atgrp->InsertAtom(aptr);
					  aptr = pres_b->GetAtomByName("CD2"); if(aptr)atgrp->InsertAtom(aptr);
					  aptr = pres_b->GetAtomByName("CB");  if(aptr)atgrp->InsertAtom(aptr);
					  aptr = pres_b->GetAtomByName("HD2"); if(aptr)atgrp->InsertAtom(aptr);
					  aptr = pres_b->GetAtomByName("HD1"); if(aptr)atgrp->InsertAtom(aptr);
					  aptr = pres_b->GetAtomByName("HE1"); if(aptr)atgrp->InsertAtom(aptr);
					  aptr = pres_b->GetAtomByName("HB2"); if(aptr)atgrp->InsertAtom(aptr);
					  aptr = pres_b->GetAtomByName("HB3"); if(aptr)atgrp->InsertAtom(aptr);
				   }
				}
			}
			
			AltChemState* p_alt_chem_st = AddAltChemState(atgrp);
			p_alt_chem_st->id = "HEM_RED";
			if( heme_model == 1)
			{
				p_alt_chem_st->SetAltChForAtom("FE",-1.0);	
			}
			else
			{
			  p_alt_chem_st->SetAltChForAtom("FE",-0.24);	
			  p_alt_chem_st->SetAltChForAtom("NA",-0.04);	
			  p_alt_chem_st->SetAltChForAtom("NB",-0.04);	
			  p_alt_chem_st->SetAltChForAtom("NC",-0.04);	
			  p_alt_chem_st->SetAltChForAtom("ND",-0.04);	
			  p_alt_chem_st->SetAltChForAtom("C1A",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("C2A",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("C3A",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("C4A",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("CHA",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("C1B",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("C2B",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("C3B",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("C4B",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("CHB",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("C1C",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("C2C",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("C3C",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("C4C",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("CHC",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("C1D",-0.03);	
			  p_alt_chem_st->   SetAltChForAtom("C2D",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("C3D",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("C4D",-0.03);	
			  p_alt_chem_st->SetAltChForAtom("CHD",-0.03);
			}
			p_alt_chem_st->alt_state_type =  p_alt_chem_st->alt_state_type.REDUCED;
			p_alt_chem_st->std_pk = -0.2;
			p_alt_chem_st->pk = p_alt_chem_st->std_pk;
			p_alt_chem_st->mod_atom_name = "FE";
		}
	}

	if(name_mod == "NT")
	{
		AltChemState* p_alt_chem_st = AddAltChemState(pres);
		p_alt_chem_st->id       = "NEND_UNPROT";
		if(res_sh_name == "PRO")
		{
			p_alt_chem_st->SetAltChForAtom("N", -0.4);		
			p_alt_chem_st->SetAltChForAtom("H2",-0.3);	
			p_alt_chem_st->SetAltChForAtom("H3",-0.3);	
		}
		else
		{
			p_alt_chem_st->SetAltChForAtom("N", -0.25);	
			p_alt_chem_st->SetAltChForAtom("H1",-0.25);	
			p_alt_chem_st->SetAltChForAtom("H2",-0.25);	
			p_alt_chem_st->SetAltChForAtom("H3",-0.25);	
		}
        p_alt_chem_st->alt_state_type =  p_alt_chem_st->alt_state_type.UNPROTONATED;
        p_alt_chem_st->std_pk = 8.0;
		p_alt_chem_st->pk = p_alt_chem_st->std_pk;
        p_alt_chem_st->mod_atom_name = "N";
	}
	return TRUE;
}

int ProtonRedoxMod:: SetStdResPKa(HaResidue* pres, int set_redox_pot)
{
	std::string res_full_name = pres->GetFullName();
	std::string res_sh_name = pres->GetName();
	std::string name_mod    = pres->GetNameModifier();

	ClearResAltChemStates(pres);

	SetStdResPKa_G1(pres,set_redox_pot);

	if(res_sh_name == "ASP")
	{
		AltChemState* p_alt_chem_st = AddAltChemState(pres);
		p_alt_chem_st->id = "ASP_PROT";
		p_alt_chem_st->SetAltChForAtom("OD1",0.5);	
		p_alt_chem_st->SetAltChForAtom("OD2",0.5);	
        p_alt_chem_st->alt_state_type =  p_alt_chem_st->alt_state_type.PROTONATED;
        p_alt_chem_st->std_pk = 3.86;
		p_alt_chem_st->pk = p_alt_chem_st->std_pk;
        p_alt_chem_st->mod_atom_name = "OD1";
	}

	if(res_sh_name == "GLU")
	{
		AltChemState* p_alt_chem_st = AddAltChemState(pres);
		p_alt_chem_st->id = "GLU_PROT";
		p_alt_chem_st->SetAltChForAtom("OE1",0.5);	
		p_alt_chem_st->SetAltChForAtom("OE2",0.5);	
        p_alt_chem_st->alt_state_type =  p_alt_chem_st->alt_state_type.PROTONATED;
        p_alt_chem_st->std_pk = 4.07;
		p_alt_chem_st->pk = p_alt_chem_st->std_pk;
        p_alt_chem_st->mod_atom_name = "OE1";
	}

	if(res_sh_name == "LYS")
	{
		AltChemState* p_alt_chem_st = AddAltChemState(pres);
		p_alt_chem_st->id       = "LYS_UNPROT";
		p_alt_chem_st->SetAltChForAtom("NZ", -0.25);	
		p_alt_chem_st->SetAltChForAtom("HZ1",-0.25);	
		p_alt_chem_st->SetAltChForAtom("HZ2",-0.25);	
		p_alt_chem_st->SetAltChForAtom("HZ3",-0.25);	
        p_alt_chem_st->alt_state_type =  p_alt_chem_st->alt_state_type.UNPROTONATED;
        p_alt_chem_st->std_pk = 10.53;
		p_alt_chem_st->pk = p_alt_chem_st->std_pk;
        p_alt_chem_st->mod_atom_name = "NZ";
	}

	if(res_sh_name == "ARG")
	{
		AltChemState* p_alt_chem_st = AddAltChemState(pres);
		p_alt_chem_st->id       = "ARG_UNPROT";
 		p_alt_chem_st->SetAltChForAtom("NH1", -0.2);	
		p_alt_chem_st->SetAltChForAtom("HH11",-0.15);
		p_alt_chem_st->SetAltChForAtom("HH12",-0.15);
 		p_alt_chem_st->SetAltChForAtom("NH2", -0.2);	
		p_alt_chem_st->SetAltChForAtom("HH21",-0.15);
		p_alt_chem_st->SetAltChForAtom("HH22",-0.15);
        p_alt_chem_st->alt_state_type =  p_alt_chem_st->alt_state_type.UNPROTONATED;
        p_alt_chem_st->std_pk = 12.1;
		p_alt_chem_st->pk = p_alt_chem_st->std_pk;
        p_alt_chem_st->mod_atom_name = "NH1";
	}

	if(name_mod == "CT")
	{
		AltChemState* p_alt_chem_st = AddAltChemState(pres);
		p_alt_chem_st->id       = "CEND_PROT";
		p_alt_chem_st->SetAltChForAtom("O",0.5);	
		p_alt_chem_st->SetAltChForAtom("OXT",0.5);	
        p_alt_chem_st->alt_state_type =  p_alt_chem_st->alt_state_type.PROTONATED;
        p_alt_chem_st->std_pk = 3.1;
		p_alt_chem_st->pk = p_alt_chem_st->std_pk;
        p_alt_chem_st->mod_atom_name = "OXT";
	}
	return TRUE;
}

void ProtonRedoxMod::SetStdPKa()
{
	HaMolSet* pmset = GetMolSet();
	HaResidue* rptr;
    ResidueIteratorMolSet ritr(pmset);
    for(rptr = ritr.GetFirstRes(); rptr; rptr = ritr.GetNextRes() )
	{
		SetStdResPKa(rptr,set_std_redox_pot);
	}
}

void ProtonRedoxMod::SetStdPKa_G1()
{
	HaMolSet* pmset = GetMolSet();
	HaResidue* rptr;
    ResidueIteratorMolSet ritr(pmset);
    for(rptr = ritr.GetFirstRes(); rptr; rptr = ritr.GetNextRes() )
	{
		SetStdResPKa_G1(rptr,set_std_redox_pot);
	}
}

bool ProtonRedoxMod::SetStdPKforAtName(HaResidue* pres, const char* at_name,double std_pk_new)
{
	AltChemState* alt_res_st = GetResAltChemStateByAtName(pres,at_name);
	if( alt_res_st == NULL)
		return false;
	alt_res_st->std_pk = std_pk_new;
	return true;
}

void ProtonRedoxMod::SetAltStatesActive(int set_flag)
{
	HaMolSet* pmset = GetMolSet();
	HaResidue* rptr;
    ResidueIteratorMolSet ritr(pmset);
    for(rptr = ritr.GetFirstRes(); rptr; rptr = ritr.GetNextRes() )
	{
		if(rptr->HasSelectedAtoms() && (GetNumResAltChemStates(rptr) > 0) )
		{
			int nalt = GetNumResAltChemStates(rptr);
			int ialt;
			for(ialt = 0; ialt < nalt;ialt++)
			{
				AltChemState* alt_st = GetResAltChemState(rptr,ialt);
				if(set_flag)
				{
					alt_st->active_flag = TRUE;
				}
				else
				{
					alt_st->active_flag = FALSE;
				}
			}
		}
	}
}

bool ProtonRedoxMod::SetResChargesForPH( HaResidue* pres, double pH_val)
{
	int ist;
    int nst = GetNumResAltChemStates(pres);
	if(nst == 0) return false;

	pres->SetStdCharges();

	std::string res_fname = pres->GetFullName();

	for(ist=0; ist < nst; ist++)
	{
	    AltChemState* alt_res_st = GetResAltChemState(pres,ist);
		if(alt_res_st == NULL) continue;

		double delt_pk = alt_res_st->pk - pH_val;

		double ratio = exp(log(10.)*delt_pk); // ratio is large for small pH - favors protonated state
        
		double weight_prot = ratio/(1.0 + ratio); // weight of the protonated state

		if(alt_res_st->alt_state_type == alt_res_st->alt_state_type.UNPROTONATED)
		{
			weight_prot = 1.0/(1.0 + ratio);
		}

		alt_res_st->SetAltCharges(weight_prot);
	}
    return true;
}

void ProtonRedoxMod::CalcPKaForSelection()
{
	CalcPKaForSelection(false);
}

void ProtonRedoxMod::CalcPKaForSelection(bool pnp)
{
	int i;
	HaMolSet* pmset = GetMolSet();
	AtomGroup sel_atoms;
	AtomIteratorMolSet aitr(pmset);
	HaAtom* aptr;
	//bool pnp=true;
	
	for(aptr = aitr.GetFirstAtom();aptr; aptr = aitr.GetNextAtom())
	{
		if(aptr->Selected())
			sel_atoms.InsertAtom(aptr);
	}
	
	int na = sel_atoms.size();

	ElectrostMod* elmod;
	elmod = pmset->GetElectrostMod(1);
	
	HaResidue* rptr; 

    ResidueIteratorMolSet ritr(pmset);
    ResidueIteratorMolSet ritr2(pmset);
	VecPtr act_chem_states;

   //>mikola's adding
	TiXmlDocument Doc("OutCalcPKs.xml");
  
	TiXmlDeclaration Decl("1.0","UTF-8","yes");
	Doc.InsertEndChild(Decl);
  
	TiXmlElement RootElt("HaMolSet::CalcPKs");
	
	vector<string> AltNames;
//>mikola's adding
//
// Find active alternative chemical states (titratable or redox-active residues)
// and set standard charge distributions of these residues
//

	AltChemState* alt_st;

    for(rptr = ritr.GetFirstRes(); rptr; rptr = ritr.GetNextRes() )
	{
		if(rptr->HasSelectedAtoms() && (GetNumResAltChemStates(rptr) > 0) )
		{
			int ist;
			int nst = GetNumResAltChemStates(rptr);
			for(ist = 0; ist < nst; ist++)
			{
				AltChemState* alt_st = GetResAltChemState(rptr,ist);
				if(alt_st == NULL) continue;
				if(!alt_st->active_flag) continue;

				alt_st->SetAltCharges(0.0);
				act_chem_states.push_back(alt_st);
				//>mikola's adding
				std::string TNAME=rptr->GetFullName().c_str();
				TNAME+="_AltState::";
				TNAME+=alt_st->id.c_str();
				AltNames.push_back(TNAME);
				//<mikola's adding
			}	
		}
	}

	int nst = act_chem_states.size();

 //>mikola's adding
	HaVec_double vTemp(nst);
	char XmlAtrName[32];
 
	for(i=0; i < nst; i++)
	{
		alt_st = (AltChemState*) act_chem_states[i];
		if(alt_st->alt_state_type ==  alt_st->alt_state_type.PROTONATED)
			AltNames[i]+="_AltChemState::PROTONATED";
		else if(alt_st->alt_state_type ==  alt_st->alt_state_type.UNPROTONATED)
			AltNames[i]+="_AltChemState::UNPROTONATED";
		else if( alt_st->alt_state_type ==  alt_st->alt_state_type.OXIDIZED)
			AltNames[i]+="_AltChemState::OXIDIZED";
		else if( alt_st->alt_state_type ==  alt_st->alt_state_type.REDUCED)
			AltNames[i]+="_AltChemState::REDUCED";

		vTemp[i]=alt_st->std_pk;
	}
	HaXML::SetElement(&RootElt,"ReseduesWithAltStates",&AltNames);
	HaXML::SetElement(&RootElt,"pKstd",&vTemp);
 //<mikola's adding
	HaMat_double inter_mat(nst,nst,0.0); //!< Interaction matrix between alternative sites
	PrintLog("CalcPKasForSelection0\n");

	if(!read_alt_st_inter)
	{

		for(i=0; i < nst; i++)
		{
			alt_st = (AltChemState*) act_chem_states[i];
		  //if(pnp)pnpmod->CalcAltStatePK(alt_st,&sel_atoms);
			elmod->CalcAltStatePK(alt_st,&sel_atoms);
			pmset->SelectAtoms(&sel_atoms);
		}

 		HaVec_double old_charges(na); //!< An Array of current atomic charges  
		for(i = 0; i < na; i++)
		{
			aptr = sel_atoms[i];
			old_charges[i] = aptr->GetCharge();
			aptr->SetCharge(0.0);
		}

		double xmin, xmax, ymin, ymax, zmin, zmax;
		pmset->GetMinMaxCrd(xmin, ymin, zmin, xmax, ymax, zmax);

		elmod->SetBoundaryAtoms(xmin, ymin, zmin, xmax, ymax, zmax);

		vector< AtomDoubleMap > sites; // changes of atomic charges during protonation/deprotonation or oxidation/reduction 
		double ch;

		for(i=0; i < nst; i++) // fill sites - array of maps of atoms to  changes of atomic charges during protonation/deprotonation or oxidation/reduction 
		{
			alt_st = (AltChemState*) act_chem_states[i];
			rptr = (HaResidue*) alt_st->GetHostAtomGroup();

			sites.push_back( AtomDoubleMap("temp"));
			AtomDoubleMap& smap = sites.back();

			alt_st->SetAltCharges(-1.0);

			AtomIteratorResidue aitr_res(rptr);

			for(aptr = aitr_res.GetFirstAtom(); aptr; aptr = aitr_res.GetNextAtom())
			{
				ch = aptr->GetCharge();
				if(fabs(ch) > 0.00000001)
				{
					smap[aptr] = ch;
				}	
				aptr->SetCharge(0.0);
			}	
		}

		int ist1, ist2;

		for(ist1 = 0 ; ist1 < nst ; ist1++)
		{
			AtomDoubleMap& smap1 = sites[ist1];
			AtomDoubleMap::iterator mitr;
			for(mitr = smap1.begin(); mitr != smap1.end(); mitr++)
			{
				aptr = (HaAtom*)(*mitr).first;
				ch = (*mitr).second;
				aptr->SetCharge(ch);
			}
			bool bres;

			//if(pnp)bres = pnpmod->run(RUN_FOREGROUND);
			bres = elmod->run(RUN_FOREGROUND);

			if(!bres)
			{
				PrintLog("Error in HaMolSet::CalcPKaForSelection() \n");
				PrintLog("Failed Continuum electrostatics calculations \n");
			}

			for(ist2 = 0; ist2 < nst; ist2++)
			{ 
				if(ist2 == ist1) continue;
				AtomDoubleMap& smap2 = sites[ist2];

				double ene = 0.0;
				for(mitr = smap2.begin(); mitr != smap2.end(); mitr++)
				{
					aptr = (HaAtom*)(*mitr).first;
					ch   = (*mitr).second;
					double phi;
					phi = elmod->el_pot_map.GetInterpolValAtPoint(aptr->GetX(),aptr->GetY(),aptr->GetZ());
					ene += ch*phi;
				}
				inter_mat.r0(ist1,ist2) = ene;
			}

			for(mitr = smap1.begin(); mitr != smap1.end(); mitr++)
			{
				aptr = (HaAtom*)(*mitr).first;
				ch = (*mitr).second;
				aptr->SetCharge(0.0);
			}
		}
		PrintLog("Interaction Matrix:\n");
		//for(ist1 = 0 ; ist1 < nst ; ist1++) jose
		//{
		//	for(ist2 = 0 ; ist2 < nst ; ist2++)
		//	{
		//		PrintLog("%.14e",inter_mat.r0(ist1,ist2));
		//		if(ist2!=nst-1)
		//			PrintLog("\t");
		//		else
		//			PrintLog("\n");
		//	}
		//} jose comment
	//<< jose added
		fstream matrix_file1;
		fstream matrix_file2;
		matrix_file1.open("Matrix_int223.dat",ios::out | ios::app);
		matrix_file2.open("Matrix_int257.dat",ios::out | ios::app);
		for(ist1 = 0 ; ist1 < nst ; ist1++) 
		{   
			for(ist2 = 0 ; ist2 < nst ; ist2++)
			{
				PrintLog("%.14e",inter_mat.r0(ist1,ist2));
				//if (ist1==0) //223
				//{
					matrix_file1 << inter_mat.r0(ist1, ist2);
					if(ist2!=nst-1){
						PrintLog("\t");
						matrix_file1 << "\t";
					}
					else {
						PrintLog("\n");
						matrix_file1 << "\n";
					}
				//}
				//else { 
				//	matrix_file2 << inter_mat.r0(ist1, ist2);
				//	if(ist2!=nst-1){
				//		PrintLog("\t");
				//		matrix_file2 << "\t";
				//	}
				//	else {
				//		PrintLog("\n");
				//		matrix_file2 << "\n";
				//	}
				//}
			}
		} 
		matrix_file1.close();
		matrix_file2.close();

	//>> jose added

		for(i=0; i < na; i++)
		{
			sel_atoms[i]->SetCharge(old_charges[i]);
		}
	} // end if(!read_alt_st_inter)
	else
	{	
		HaMatDB inter_mat_db;
		int ires = inter_mat_db.open("ALT_ST_INTER_MAT.DAT","r");
		if(ires) 
		{
			inter_mat_db.GetMat("ALT_ST_INTER_MAT_1", inter_mat);
			if(inter_mat.num_rows() != nst) 
			{
				PrintLog("Error Reading interaction matrix of alternative states \n");
				PrintLog("Dimensions of the matrix %d do not equal to the number active states %d \n",
					inter_mat.num_rows(), nst);
				return;
			}
		}
		else
		{
			PrintLog(" Error to open file with matrix of alternative chemical states interactions \n");
			return;
		}

		for( i = 0; i < nst; i++)
		{
			alt_st = (AltChemState*) act_chem_states[i];
			alt_st->pk = inter_mat.r0(i,i);
		}
		inter_mat_db.close();
	}

	//return; //<< jose commented

	//>mikola's adding
	HaXML::SetElement(&RootElt,"Inter",&inter_mat);
	for(i=0; i < nst; i++)
	{
		alt_st = (AltChemState*) act_chem_states[i];
		vTemp[i]=alt_st->pk;
	}
	HaXML::SetElement(&RootElt,"pKnoIter",&vTemp);
	//<mikola's adding
	if(save_alt_st_inter)
	{
		HaMatDB inter_mat_db;
		int ires = inter_mat_db.open("ALT_ST_INTER_MAT.DAT","w");
		if(ires) 
		{
			for( i = 0; i < nst; i++)
			{
				alt_st = (AltChemState*) act_chem_states[i];
				inter_mat.r0(i,i) = alt_st->pk;
			}
			inter_mat_db.PutMat("ALT_ST_INTER_MAT_1", inter_mat);
			inter_mat_db.close();
			//<<jose added
		fstream matrix_file3;
		matrix_file3.open("Matrix_intMCK.dat",ios::out | ios::app);
		for(int ist1 = 0 ; ist1 < nst ; ist1++) 
		{   
			for(int ist2 = 0 ; ist2 < nst ; ist2++)
			{
				PrintLog("%.14e",inter_mat.r0(ist1,ist2));
				matrix_file3 << inter_mat.r0(ist1, ist2);
				if(ist2!=nst-1){
					PrintLog("\t");
					matrix_file3 << "\t";
				}
				else {
					PrintLog("\n");
					matrix_file3 << "\n";
				}
			}
		} 
		matrix_file3.close();

	//>> jose added
		}
	}

	FILE* file_pop; 

	if( save_titration_data)
	{
		file_pop = fopen("ALT_ST_POP.DAT","w");
		if(file_pop == NULL) return;
	}
	else
	{
		ph_min = ph;
		ph_max = ph + 0.000001;
		ph_step = 1.0;
		e0_min = e0;
		e0_max = e0+0.000001;
		e0_step = 1.0;
	}

	HaVec_double alt_pop(nst), delt_e(nst);

	double ph_1, e0_1;
	double ph_last, e0_last;

	for( ph_1 = ph_min; ph_1 < ph_max; ph_1 += ph_step )
	{
		for( e0_1 = e0_min; e0_1 < e0_max; e0_1 += e0_step )
		{
			ph_last = ph_1;
			e0_last = e0_1;
			for(i=0; i < nst; i++)
			{
				double de;
				alt_st = (AltChemState*) act_chem_states[i];
				if(alt_st->alt_state_type ==  alt_st->alt_state_type.PROTONATED)
				{
					de = log(10.0)*(ph_1 - alt_st->pk); // Diff of energies of alternative and reference states in kT
				}
				else if(alt_st->alt_state_type ==  alt_st->alt_state_type.UNPROTONATED)
				{
					de = log(10.0)*(alt_st->pk - ph_1);
				}
				else if( alt_st->alt_state_type ==  alt_st->alt_state_type.OXIDIZED)
				{
					de = EV_TO_KT*(alt_st->pk - e0_1);
				}
				else if( alt_st->alt_state_type ==  alt_st->alt_state_type.REDUCED)
				{
					de = EV_TO_KT*(e0_1 - alt_st->pk);
				}
				inter_mat.r0(i,i) = de;
			}

			if(multi_site_pop_method == multi_site_pop_method.SCF_MULTI_SITE_CALC)
			{
				CalcAvgPopSCF(inter_mat, alt_pop, delt_e);
			}
			else if( multi_site_pop_method == multi_site_pop_method.MC_MULTI_SITE_CALC || inter_mat.num_rows() >= 32 )
			{
				CalcAvgPopMC(inter_mat, alt_pop, delt_e,n_mc_cyc);
			}
			else 
			{
				CalcAvgPopPFunc(inter_mat, alt_pop, delt_e);
			}
			//>mikola's adding
			sprintf(XmlAtrName,"Inter_pH%.3lg",ph_1);
			HaXML::SetElement(&RootElt,XmlAtrName,&inter_mat);
			sprintf(XmlAtrName,"alt_pop_pH%.3lg",ph_1);
			HaXML::SetElement(&RootElt,XmlAtrName,&alt_pop);
			sprintf(XmlAtrName,"delt_e_pH%.3lg",ph_1);
			HaXML::SetElement(&RootElt,XmlAtrName,&delt_e);
			/*printf("pH:%g Matrix:\n ",ph_1);
			for(i=0; i < nst; i++)
			{
			for(j=0; j < nst; j++)
			{
			printf("%8.6le ",inter_mat.r0(i,j));
			}
			printf("\n");
			}*/
			//<mikola's adding

			if( save_titration_data)
			{
				fprintf(file_pop," %9.4f %9.4f  ",ph_1,e0_1);
				for(i = 0; i < nst; i++)
				{
					alt_st = (AltChemState*) act_chem_states[i];
					if( save_only_redox_titr && 
						(alt_st->alt_state_type ==  alt_st->alt_state_type.PROTONATED || 
						alt_st->alt_state_type ==  alt_st->alt_state_type.UNPROTONATED) ) continue; 
					fprintf(file_pop," %9.4f ",alt_pop[i]);
				}
				fprintf(file_pop,"\n");
			}
		} // end cycle on e0
	} // end cycle on ph

	PrintLog(" Energy shifts of alternative states energies \n");
	for(i=0; i < nst; i++)
	{
		PrintLog(" %4d self_ene= %12.6f (kT)  tot_shift = %12.6f (kT) \n",
			i,inter_mat.r0(i,i), delt_e[i]);
		alt_st = (AltChemState*) act_chem_states[i];
		if(alt_st->alt_state_type ==  alt_st->alt_state_type.PROTONATED)
		{
			alt_st->pk = ph_last - delt_e[i]/log(10.0);
		}
		else if(alt_st->alt_state_type ==  alt_st->alt_state_type.UNPROTONATED)
		{
			alt_st->pk = ph_last + delt_e[i]/log(10.0); 
		}
		else if( alt_st->alt_state_type ==  alt_st->alt_state_type.OXIDIZED)
		{
			alt_st->pk = e0_last + delt_e[i]/EV_TO_KT;
		}
		else if( alt_st->alt_state_type ==  alt_st->alt_state_type.REDUCED)
		{
			alt_st->pk = e0_last - delt_e[i]/EV_TO_KT;
		}
	}
	//>mikola's adding
	for(i=0; i < nst; i++)
	{
		alt_st = (AltChemState*) act_chem_states[i];
		vTemp[i]=alt_st->pk;
	}
	HaXML::SetElement(&RootElt,"pKwithIter",&vTemp);
	//<mikola's adding
	if( save_titration_data)
	{
		fclose(file_pop);
	}
	//>mikola's adding
	Doc.InsertEndChild(RootElt);
	Doc.SaveFile();
	//<mikola's adding
}

int ProtonRedoxMod::CalcAvgPopMC(HaMat_double& inter_mat, HaVec_double& avg_st_pop, HaVec_double& alt_st_ene,int N_mc_cyc)
//!
//! interaction matrix - hamiltonian for interaction between alternative states  
//! avg_st_pop[j]     -  average population of aternative states  
//! alt_st_ene[j]     -  effective energy difference corresponding to difference in average population of alternative and unmodifed states 
//!
{
	int nst = inter_mat.num_rows();

	if( N_mc_cyc <= 0) 
	{
		PrintLog(" Error in ProtonRedoxMod::CalcAvgPopMC()  N_mc_cyc = %d   <= 0 \n", N_mc_cyc);
		return FALSE;
	}

	alt_st_ene.newsize(nst);
	alt_st_ene = 0.0;
	avg_st_pop.newsize(nst); // Average population of alternative states
	avg_st_pop = 0.0;

	HaVec_double state_vec(nst);

	int i;
	for(i= 0; i < nst; i++)
	{
		if( inter_mat.r0(i,i) < 0.0)
		{
            state_vec[i] = 1.0;
		}
		else
		{
            state_vec[i] = 0.0;
		}
	}

    int icyc;
	int nstep = 0;

	Random rand_num_gen(5);
	
	for( icyc = 0; icyc < N_mc_cyc; icyc++)
	{
		for( i = 0; i < nst; i++) // Try to flip state one at a time as an elemental monte carlo move
		{
		   int accept = TRUE;
           double delt_e = dot_double(state_vec.begin(),&inter_mat.r0(0,i), nst);
		   if( state_vec[i] > 0.5 )
		   {
               delt_e = -delt_e;
		   }
		   else
		   {
               delt_e += inter_mat.r0(i,i);
		   }
              
		   if( delt_e > 0)
		   {
			   double dt = exp(-delt_e);
               double rn = rand_num_gen();
			   if( dt < rn) accept = FALSE;
		   }
		   if( accept )  state_vec[i] = (state_vec[i] > 0.5) ? 0.0 : 1.0; // flip state 
              
		   nstep++;
		   int k;
		   for(k = 0; k < nst; k++)
		   {
			   avg_st_pop[k] += state_vec[k];
		   }
		}
	}
    
	for(i = 0; i < nst; i++)
	{
		avg_st_pop[i] = avg_st_pop[i]/nstep;
		if( (1.0 - avg_st_pop[i]) < DBL_EPSILON  ) 
		{
			alt_st_ene[i] = -100.0;
		}
		else if( avg_st_pop[i] < DBL_EPSILON ) 
		{
			alt_st_ene[i] = +100.0;
		}
		else
		{
			alt_st_ene[i] = log(1.0/avg_st_pop[i] - 1.0);
		}
	}
	return TRUE;
}

int ProtonRedoxMod::CalcAvgPopPFunc(HaMat_double& inter_mat,  HaVec_double& avg_st_pop, HaVec_double& alt_st_ene)
{
	int nst = inter_mat.num_rows();

	if( nst == 0 ) 
	{
		PrintLog(" Error in ProtonRedoxMod::CalcAvgPopDir() \n");
		PrintLog(" Interaction matrix has zero dimensions \n ");
		return FALSE;
	}

	if( nst >= 32) 
	{
		PrintLog(" Error in ProtonRedoxMod::CalcAvgPopDir() \n");
		PrintLog(" Interaction matrix dimensions %d  are too large to direct population calculations max size = %d \n",
			       nst, 32 );
		return FALSE;
	}

	alt_st_ene.newsize(nst);
	alt_st_ene = 0.0;
	avg_st_pop.newsize(nst); // Average population of alternative states
	avg_st_pop = 0.0;

	HaVec_double state_vec(nst);

	int i;
	int j,k;
	unsigned long last = 1;
	
	for( i = 0; i < nst; i++)
	{
		last *= 2;
	}

	double part_fun = 0.0;

	for( i = 0; i < last; i++)
	{
		std::bitset<32> st_bit_vec(i);
		for(j = 0; j < nst; j++)
		{
			state_vec[j] = st_bit_vec[j] ? 1.0 : 0.0 ; 
		}
		double ene = 0.0;
		
		for( j = 0; j < nst; j++)
		{
			for( k = 0; k <= j; k++ )
			{
				ene += inter_mat.r0(k,j) *  state_vec[k] * state_vec[j]; 
			}
		}
		double frac = exp(-ene);
		part_fun += frac;
		for( j = 0; j < nst; j++)
		{
			avg_st_pop[j] += state_vec[j]* frac;
		}
	}

	for( j = 0; j < nst; j++ )
	{
		avg_st_pop[j] = avg_st_pop[j]/part_fun;
	}

	for(j = 0; j < nst; j++)
	{
		if( (1.0 - avg_st_pop[j]) < DBL_EPSILON  ) 
		{
			alt_st_ene[j] = -100.0;
		}
		else if( avg_st_pop[j] < DBL_EPSILON ) 
		{
			alt_st_ene[j] = +100.0;
		}
		else
		{
			alt_st_ene[j] = log(1.0/avg_st_pop[j] - 1.0);
		}
	}
	return TRUE;
}


int ProtonRedoxMod::CalcAvgPopSCF(HaMat_double& inter_mat, HaVec_double& alt_pop, HaVec_double& delt_e)
//!
//! interaction matrix - on the diagonal energy difference of alternative and unmodifed state when 
//! alt_pop[j] - population of aternative states  
//! delt_e[j]  - converged self-consistent difference between alternative and unmodifed state 
//!
{
	int nst = inter_mat.num_rows();

	delt_e.newsize(nst);
	delt_e = 0.0;
	alt_pop.newsize(nst); // population of alternative states
	alt_pop = 0.0;

	int iter = 0;
	int max_iter = 100;
	double pop_err_max = 0.0001;
	int i,j;

	int converged = FALSE;

	for(iter = 0 ; iter < max_iter; iter++)
	{
//		PrintLog("Iter %4d :\n",iter); 
//		PrintLog("Popul    :\n"); 
		
//		PrintLog("pKa      :\n"); 
		double pop_err = 0.0;
		for(i=0; i < nst; i++)
		{
			delt_e[i] = inter_mat.r0(i,i);
			for(j = 0; j < nst; j++)
			{
				if(j == i) continue;
				delt_e[i] += inter_mat.r0(i,j) * alt_pop[j];
			}
		}

		for(i=0; i<nst;i++)
		{
			double pop_old = alt_pop[i];
			double pop_new = 1.0/(1.0 + exp(delt_e[i]));
//			PrintLog(" %d %8.3f",i,alt_pop[i]);
			if( iter > 0)
			{
                alt_pop[i] = pop_old*0.6 + pop_new*0.4;
			}
			else
			{
                alt_pop[i] = pop_new;
			}
			if( fabs(pop_new - pop_old) > pop_err) pop_err = fabs(pop_new - pop_old);
		}
//		PrintLog("\n");         

//		PrintLog(" Max population change: %9.5f \n",pop_err);

		if(pop_err < pop_err_max)
		{
//			PrintLog(" Population values have converged \n");
			converged = TRUE;
			break;
		}
	}
	
	if(!converged)
	{
		PrintLog(" Maximum Number of population convergence iteration exceeded");
	}

	return TRUE;
}

int ProtonRedoxMod::TestCalcPopFun()
{
	int nst = 5;
	HaMat_double inter_mat(nst,nst);

	int i,j;
	for( i = 0; i < nst; i++)
	{
		inter_mat.r0(i,i) = 0.5*i;
	}
	
	for( i = 0; i < nst; i++ )
	{
		for( j = 0; j < i; j++)
		{
			inter_mat.r0(i,j) = 0.1 + 0.05 * j;
			inter_mat.r0(j,i) = 0.1 + 0.05 * j;
		}
	}

	HaVec_double avg_st_pop(nst); 
	HaVec_double alt_st_ene(nst);

	CalcAvgPopPFunc(inter_mat, avg_st_pop, alt_st_ene);
	
	PrintLog(" CalcAvgPopDir() results: \n");
	for( i = 0; i < nst; i++)
	{
		PrintLog(" %d  avg_pop=%16.9f  ene=%16.9f \n",  i, avg_st_pop[i],alt_st_ene[i]); 
	}
	PrintLog("\n");

	int n_mc = 10000000;
	CalcAvgPopMC(inter_mat, avg_st_pop, alt_st_ene,n_mc); 

	PrintLog(" CalcAvgPopMC() results: \n");
	for( i = 0; i < nst; i++)
	{
		PrintLog(" %d  avg_pop=%16.9f  ene=%16.9f \n",  i, avg_st_pop[i],alt_st_ene[i]); 
	}
	PrintLog("\n");

	CalcAvgPopSCF(inter_mat, avg_st_pop, alt_st_ene);

	PrintLog(" CalcAvgPopSCF() results: \n");
	for( i = 0; i < nst; i++)
	{
		PrintLog(" %d  avg_pop=%16.9f  ene=%16.9f \n",  i, avg_st_pop[i],alt_st_ene[i]); 
	}
	PrintLog("\n");

	return TRUE;
}

bool ProtonRedoxMod::SetChargesForPH(double pH_val)
{
	HaMolSet* pmset = GetMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);
	p_mol_editor->SetStdAtomicParams(pmset,AMBER_ALL_ATOM_CHRGS);
	HaResidue* rptr;
    ResidueIteratorMolSet ritr(pmset);
    for(rptr = ritr.GetFirstRes(); rptr; rptr = ritr.GetNextRes() )
	{
		SetResChargesForPH(rptr,pH_val);
	}
	return true;
}

bool ProtonRedoxMod::SetChargesForCurrentPH()
{
	return SetChargesForPH(GetPH());
}

void ProtonRedoxMod::SetPH(double new_ph)
{
	ph = new_ph;

}

double ProtonRedoxMod::GetPH() const
{
	return ph;
}
