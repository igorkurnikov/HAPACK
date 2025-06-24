/*! \file mm_model.cpp

    Classes to setup Molecular Mechanics Model
 
    \author Igor Kurnikov 
    \date 2010-
*/

#define MM_MODEL_CPP

#include <float.h>
#include <math.h>
#include "mpi.h"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include "hawx_add.h"
#include "FCMangle.h"

#include "harlemapp.h"
#include "hampi.h"
#include "hamolset.h"
#include "moleditor.h"
#include "hamolmech.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "mm_force_field.h"
#include "mm_driver_amber.h"
#include "haresdb.h"

extern "C" 
{
	extern void FC_FUNC_MODULE(amoeba_induced_mod,get_ind_dip)(int* natom, double* ind_dip_d_out, double* ind_dip_p_out );
	extern void FC_FUNC_MODULE(amoeba_multipoles_mod,get_global_multipole)( int* natom, double* global_multipole );
}

MolMechModel::MolMechModel(HaMolMechMod* p_mm_mod_new)
{
	p_mm_mod = p_mm_mod_new;
	pmset = p_mm_mod->GetMolSet();
	p_amber_model = new AmberMMModel(this);
	p_mort_model = NULL;

	electr_method.SetMMModel(this);

	to_init_mm_model = TRUE;
	SetStdParams();
}
	
MolMechModel::~MolMechModel()
{

}

int MolMechModel::SetStdParams()
{
	build_nb_contact_list_flag = TRUE;
	init_charges_flag = TRUE;
	SetUseMortLib(FALSE);

	SetDielConst(1.0);
	SetIonStrength(0.0);
	gb_param_type  = MOD_BONDI_GB_PARAM;
	SetScale14Electr(1.2);
	SetScale14VdW(2.0);

	SetNBCutDist(9.0);

	double tmp;
	if( pmset->per_bc->IsSet() )
	{
		tmp = nb_cut_dist; 
		p_mm_mod->period_bcond = p_mm_mod->period_bcond.CONST_PRES;
		if( nb_cut_dist > 0.35*fabs(pmset->per_bc->GetA()) ) tmp = 0.35*fabs(pmset->per_bc->GetA());
		if( nb_cut_dist > 0.35*fabs(pmset->per_bc->GetB()) ) tmp = 0.35*fabs(pmset->per_bc->GetB());
		if( nb_cut_dist > 0.35*fabs(pmset->per_bc->GetC()) ) tmp = 0.35*fabs(pmset->per_bc->GetC());

		MMDriverAmber::ModifyFormatVal(tmp,FLOAT_F8_3);
		SetNBCutDist(tmp);
	}

	nonb_list_flag = TRUE;
	calc_vdw_flag = CALC_VDW_NORMAL;
	calc_electr_flag = True;
	neutral_end_hydr_flag = 0;
	soft_repulsion = NO_SOFT_REPULSION;
	soft_repulsion_const = 0.0;

	subtract_avg_force_flag = -1;

	moving_atoms = "ALL_ATOMS";
    restrained_atoms = "NO_RESTRAINTS";
	if (pmset->GetAtomGroupByID("RESTRAINED_ATOMS")) SetRestrainedAtoms("RESTRAINED_ATOMS");


	restr_ref_crd_type = RESTR_REFC_CURRENT_CRD;

	atom_restr_const = 1.0;

// Water cap treatment:

	water_cap_flag = False;
	cap_atom_num = 0;
	cap_fconst = 1.5;

// Particle Mesh Ewald:
   
   pmesh_ewald_flag = 0; 
   pme_grid_nx = pme_grid_ny = pme_grid_nz = 0; //!< Set from box size  
   vdw_correction_flag = 1;
   pme_spline_order = 4;
   pme_neutralize_system = 1;
   pme_verbose = 0;
   exact_ewald_flag = 1;
   pme_dsum_tol = 1.0e-5; 
   pme_ew_coeff = 0.0;
   skin_nb      = 2.0;
   fft_grids_per_ang = 1.0;
   pme_eedtbdns = 5000.0;

   dipole_scf_iter_max = 50;

//   dipole_scf_tol = 0.01;
   dipole_scf_tol = 0.001;
   ee_dsum_cut = 7.0;
   ee_damped_cut = 4.5;
   sor_coefficient = 0.75;
   thole_expon_coeff = 0.0; //!< thole_ij = thole_i/sqrt(0.39) * thole_j/sqrt(0.39)  ???
// thole_expon_coeff = 0.39;
// thole_expon_coeff = 0.095;
   vdw_taper = 0.9;

   subtract_avg_force_flag = TRUE;
   ff_type = ForceFieldType::UNKNOWN_FF;

   return TRUE;
}

void MolMechModel::Bcast(MPI_Comm& comm)
{
	int ires;
//	PrintLog(" MolMechModel::Bcast() pt 1 \n");
	ires = MPI_Bcast(&omit_interactions,1,MPI_INT,0,comm); // AMBER ntf
	ires = MPI_Bcast(&nb_cut_dist,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&subtract_avg_force_flag,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&diel_const,1,MPI_DOUBLE,0,comm);      // AMBER dielc
	ires = MPI_Bcast(&scale_14_vdw,1,MPI_DOUBLE,0,comm);    // AMBER scnb
	ires = MPI_Bcast(&scale_14_electr,1,MPI_DOUBLE,0,comm); // AMBER scee
	ires = MPI_Bcast(&nb_cut_dist,1,MPI_DOUBLE,0,comm); // AMBER gb_cutoff  ??
	ires = MPI_Bcast(&electr_method,1,MPI_INT,0,comm);        // AMBER igb

	ires = MPI_Bcast(&ff_type.value(),1,MPI_INT,0,comm);  // Force Field type

	p_amber_model->Bcast(comm);
//	PrintLog(" MolMechModel::Bcast() pt end \n");
}

int MolMechModel::SaveXMLToStream(std::ostream& os, const harlem::SaveOptions* popt ) const
{
	const MolSet* pmset_loc = GetMolSet();
	os << "<mm_model mset=\"" << pmset_loc->GetName() << "\">" << std::endl;
		
	char buf[256];

	if( !ImprDihedrals.empty())
	{
		for(shared_ptr<MMDihedral> ps_dih : ImprDihedrals )
		{
			const MMDihedral& impr_dihedral = *ps_dih;
			if( impr_dihedral.pt1 == NULL || impr_dihedral.pt2 == NULL || 
				impr_dihedral.pt3 == NULL || impr_dihedral.pt4 == NULL )
			{
				continue;
			}
			os << "<impr>";

			const HaAtom* aptr = (const HaAtom*) impr_dihedral.pt1;
			aptr->FillRef(buf);
			os << buf << "  ";

			aptr = (HaAtom*) impr_dihedral.pt2;
			aptr->FillRef(buf);
			os << buf << "  ";

			aptr = (HaAtom*) impr_dihedral.pt3;
			aptr->FillRef(buf);
			os << buf << "  ";

			aptr = (HaAtom*) impr_dihedral.pt4;
			aptr->FillRef(buf);
			os << buf << "  " << std::endl;

			os << "</impr>" << std::endl;
		}
	}

	int nv = 0;
	set<MMBond, less<MMBond> >::const_iterator mbitr = MBonds.begin();

	for(; mbitr != MBonds.end(); mbitr++)
	{
		const MMBond& bnd = (const MMBond&) *mbitr;
		if( bnd.set_type == MolMechModel::SET_SPEC)
		{
			os << "<bond set=\"spec\" ";

			sprintf(buf,"r0=\"%12.6f\"  fc=\"%12.6f\" ", bnd.r0, bnd.fc);
			os << buf << ">" << std::endl;

			const HaAtom* aptr = (const HaAtom*) bnd.pt1;
			aptr->FillRef(buf);
			os << buf << "  ";

			aptr = (HaAtom*) bnd.pt2;
			aptr->FillRef(buf);

			os << buf;

			os << "</bond>" << std::endl;
		}
	}

	set<MMValAngle, less<MMValAngle> >::const_iterator vaitr = ValAngles.begin();

	for(; vaitr != ValAngles.end(); vaitr++)
	{
		const MMValAngle& vang = (const MMValAngle&) *vaitr;
		if( vang.set_type == MolMechModel::SET_SPEC)
		{
			os << "<angle set=\"spec\" ";

			sprintf(buf,"a0=\"%12.6f\"  fc=\"%12.6f\" ",vang.a0,vang.fc);
			os << buf << std::endl;
			os << ">";

			const HaAtom* aptr = (const HaAtom*) vang.pt1;
			aptr->FillRef(buf);
			os << buf << "  ";

			aptr = (HaAtom*) vang.pt2;
			aptr->FillRef(buf);
			os << buf << "  ";

			aptr = (HaAtom*) vang.pt3;
			aptr->FillRef(buf);
			os << buf << "  ";
			
			os << "</angle>" << std::endl;
		}  
	}
	os << "</mm_model>" << std::endl;
	return TRUE;
}


int MolMechModel::IsAmoebaFF() const
{
	return (ff_type == ForceFieldType::AMOEBA);
}


template<typename T>
using FourPtrs = std::array<T*, 4>;

// Helper: return a canonical copy (the smaller of a and reverse(a))
template<typename T>
FourPtrs<T> canonical(const FourPtrs<T>& a) {
	FourPtrs<T> rev = a;
	std::reverse(rev.begin(), rev.end());
	// choose lexicographically smaller
	if (std::lexicographical_compare(rev.begin(), rev.end(),
		a.begin(), a.end()))
		return rev;
	else
		return a;
}

// Hash functor that hashes the canonical form
template<typename T>
struct FourPtrsHash {
	std::size_t operator()(FourPtrs<T> const& a_) const noexcept {
		FourPtrs<T> a = canonical(a_);
		std::size_t seed = 0;
		std::hash<T*> hasher;
		for (T* p : a) {
			// boost‐style mix
			seed ^= hasher(p)
				+ 0x9e3779b9
				+ (seed << 6)
				+ (seed >> 2);
		}
		return seed;
	}
};

// Equality allows either direct or reverse match
template<typename T>
struct FourPtrsEq {
	bool operator()(FourPtrs<T> const& x,
		FourPtrs<T> const& y) const noexcept
	{
		// exact match?
		if (x == y) return true;
		// reverse match?
		for (int i = 0; i < 4; ++i)
			if (x[i] != y[3 - i]) return false;
		return true;
	}
};


using FourAtoms = FourPtrs<HaAtom*>;

int MolMechModel::InitModel(const ForceFieldType& ff_type_par )
{
	int ires;
	Clear();

	MMForceField* p_ff = MMForceField::GetMMForceField(ff_type_par,TRUE);
	if( p_ff == NULL)
	{
		PrintLog("Error in MolMechModel::InitModel() \n");
		PrintLog("Can not initiate Force Field %s \n",ff_type_par.label());
		return FALSE;
	}
	ff_type = p_ff->GetFFType();

	for(HaAtom* aptr : *pmset)
	{
		Atoms.push_back(aptr);
	}

	if( ff_type == ForceFieldType::AMOEBA ) // TMP IGOR should consider to turn it off
	{
		ClearMortModel();
		try
		{
			ires = pmset->SetMortMol(*p_mort_model, ff_type );
			if(!ires) throw std::runtime_error(" Error in SetMortMol() \n");
		
			ires = p_amber_model->InitAmberModelAmoeba();
			if(!ires) throw std::runtime_error(" Error in InitAmberModelAmoeba() \n");
			ires = p_amber_model->UpdateAmberData();
			if(!ires) throw std::runtime_error(" Error in UpdateAmberData() \n"); 
			
			to_init_mm_model = FALSE;
			if( pme_spline_order < 5 ) pme_spline_order = 5;
		}
		catch( std::exception& ex) 
		{	
			PrintLog(" Error in MolMechModel::InitModel() \n");
			PrintLog(" Error to set AMOEBA model for the molecular set \n");
			ff_type = ForceFieldType::UNKNOWN_FF;
			to_init_mm_model = TRUE;
			return FALSE;
		}
		return FALSE;
	}

// Map atoms to their conterparts in the residue templates

	AtomAtomMap pt_to_templ_map;        // Map Atoms to Atom Templates 
	AtomAtomMap pt_to_mut_templ_map;    // Map Atoms to Mutated Atom Templates 

	HaResDB* p_res_db = HaResDB::GetDefaultResDB(); 

	ResidueIteratorMolSet ritr(pmset);
	for (HaResidue* pres : ritr)
	{
		std::string res_fname = pres->GetFullName();
		// if (res_fname == "TRM") continue; //  Termination groups are treated separately
		HaResidue* p_res_templ = p_res_db->GetTemplateForResidue(res_fname);
		if (p_res_templ == NULL)
		{
			PrintLog(" Can't find residue template %s for residue %s\n",res_fname, pres->GetRef());
			PrintLog(" No atom FF parameters will be set for atoms of the residue \n");
			continue;
		}

		ResFFTemplate* p_res_ff_templ = p_ff->GetResidueTemplate(res_fname);
		if (p_res_ff_templ == NULL)
		{
			PrintLog(" Force field %s  doesn't have specific parameters for residue name %s\n", p_ff->GetFFType().label(), res_fname);
			PrintLog(" FF parameters will be set from residue template: %s for residue: %s \n", pres->GetFullName(), pres->GetRef());
		}

		std::string res_mut_fname;
		HaResidue* p_res_mut_templ = nullptr;
		ResFFTemplate* p_res_mut_ff_templ = nullptr;

		if (pres->IsAlchemicalTransformationSet())
		{
			p_mm_mod->run_ti = TRUE; 
			res_mut_fname = pres->p_res_transform->res_name_b;
			p_res_mut_templ = p_res_db->GetTemplateForResidue(res_mut_fname);
			if (!p_res_mut_templ)
			{
				PrintLog(" Can't find mutated residue template %s for residue %s\n", res_mut_fname, pres->GetRef());
				PrintLog(" the mutated state of the residue will not be set\n");
			}

			p_res_mut_ff_templ = p_ff->GetResidueTemplate(res_mut_fname);
			if (p_res_mut_templ && !p_res_mut_ff_templ)
			{
				PrintLog(" Force field %s doesn't have specific parameters for residue name %s\n", p_ff->GetFFType().label(), res_mut_fname);
				PrintLog(" FF parameters for mutated state of residue: %s will be set from residue template: %s \n", pres->GetRef(), res_mut_fname );
			}
		}

		for (HaAtom* aptr : *pres)  // Set Atom FF parameters from FF Residue templates
		{
			std::string at_name = aptr->GetName();

			HaAtom* atempl = p_res_templ->GetAtomByName(at_name);

			pt_to_templ_map[aptr] = atempl;

			double charge = atempl ? atempl->GetCharge() : 0.0;
			double mass = atempl ? atempl->GetMass() : 4.0;
			std::string ff_symb = atempl ? atempl->GetFFSymbol() : "DU";

			if (atempl && p_res_ff_templ != NULL)
			{
				shared_ptr<AtomFFParam> p_at_ff = p_res_ff_templ->GetAtomFFParam(aptr->GetName());
				if (p_at_ff != NULL)
				{
					charge = p_at_ff->GetCharge();
					ff_symb = p_at_ff->ff_symbol;
				}
				else
				{
					PrintLog(" Can't find parameters for atom %s in residue template %s ", aptr->GetRef(), p_res_ff_templ->GetFullName());
				}
			}
			aptr->SetFFSymbol(ff_symb);
			aptr->SetCharge(charge);
			aptr->SetMass(mass);
		}

		for (HaAtom* aptr : *pres)  // Set Atom FF parameters from FF Residue templates for mutated state
		{
			HaAtom* atempl_mut = nullptr;
			std::string at_name_mut;
			if (pres->IsAlchemicalTransformationSet() && p_res_mut_templ)
			{
				at_name_mut = pres->p_res_transform->at_names_b[aptr];
				atempl_mut = p_res_mut_templ->GetAtomByName(at_name_mut);
				if (atempl_mut) pt_to_mut_templ_map[aptr] = atempl_mut;
			}

			if (pres->IsAlchemicalTransformationSet() && atempl_mut)
			{
				shared_ptr<AtomFFParam> p_at_mut_params = make_shared<AtomFFParam>();
				p_at_mut_params->charge = atempl_mut->GetCharge();
				p_at_mut_params->ff_symbol = atempl_mut->GetFFSymbol();
				p_at_mut_params->mass = atempl_mut->GetMass();

				if (p_res_mut_ff_templ)
				{
					shared_ptr<AtomFFParam> p_at_mut_ff_templ = p_res_mut_ff_templ->GetAtomFFParam(at_name_mut);
					if (p_at_mut_ff_templ) *p_at_mut_params = *p_at_mut_ff_templ;
				}

				pres->p_res_transform->at_ff_params[aptr] = p_at_mut_params;
			}
		}

		for (HaAtom* aptr : *pres)
		{
			if (pt_to_templ_map.count(aptr) > 0) continue; // Set FF parameters for atoms not found in residue templates
			std::string at_name = aptr->GetName();
			if (boost::algorithm::starts_with(at_name, "HTM")) // special rules for Terminating hydrogens
			{
				if (aptr->IsHydrogen())
				{
					AtomGroup bonded_atoms = aptr->GetBondedAtoms();
					if (bonded_atoms.size() != 1)
					{
						PrintLog("TRM hydrogen atom %s has humber of bonds not equal 1 - FF Params are not set \n", aptr->GetRef());
						continue;
					}

					HaAtom* aptr_host = bonded_atoms[0];
					AtomGroup bonded_atoms_host = aptr_host->GetBondedAtoms();
					int num_hydrogens_host = 0;
					for (HaAtom* aptr_b : bonded_atoms_host)
					{
						if (aptr_b->IsHydrogen()) num_hydrogens_host++;
					}

					std::string at_host_ff_s = aptr_host->GetFFSymbol();

					if (bonded_atoms.size() != 4)
					{
						PrintLog("Atom %s bonded to Hydrogen TRM atom %s has %d bonds (should be 4 -  SP4 carbon) - FF Params are not set",
							aptr_host->GetRef(), aptr->GetRef(), bonded_atoms.size());
						continue;
					}

					if (aptr_host->GetElemNo() != 6)
					{
						PrintLog("Atom %s bonded to TRM atom %s is not carbon \n - FF Params are not set ",
							aptr_host->GetRef(), aptr->GetRef());
						continue;
					}

					if (ff_type == ForceFieldType::AMBER_94 || ff_type == ForceFieldType::AMBER_99_SB || ff_type == ForceFieldType::AMBER_99_BSC0 ||
						ff_type == ForceFieldType::AMBER_03 || ff_type == ForceFieldType::AMBER_10)
					{
						if (at_host_ff_s == "CT")
						{
							std::string h_at_s = "H1";
							for (HaAtom* aptr_b : bonded_atoms_host)
							{
								if (aptr->IsHydrogen() && aptr->GetFFSymbol() == "HC") h_at_s = "HC"; // terminal Methyl hydrogens should be all H1 or all HC
							}
							aptr->SetFFSymbol(h_at_s);
							aptr->SetCharge(0.0);
						}
					}
					else if (ff_type == ForceFieldType::ARROW_5_14_CT || ff_type == ForceFieldType::ARROW_2_0)
					{

						if (num_hydrogens_host == 2)
						{
							aptr_host->SetFFSymbol("C2");
							for (HaAtom* aptr_b : bonded_atoms)
								if (aptr_b->IsHydrogen()) aptr_b->SetFFSymbol("HC2");
						}

						if (num_hydrogens_host == 3)
						{
							aptr_host->SetFFSymbol("C3");
							for (HaAtom* aptr_b : bonded_atoms)
								if (aptr_b->IsHydrogen()) aptr_b->SetFFSymbol("HC3");
						}
					}
				}
				else
				{
					PrintLog("Terminating Atom %s is not Hydrogen - FF Params are not set\n", aptr->GetRef());
					continue;
				}
			}
		}

		for (HaAtom* aptr : *pres)
		{
			std::string at_name = aptr->GetName();
			if (boost::algorithm::starts_with(at_name, "DA")) // special rules for State A Dummy atoms in Mutating residues
			{
				double mass = 4.0;
				if (pres->p_res_transform->at_ff_params.count(aptr) > 0)
				{
					shared_ptr<AtomFFParam> pat_par_mut = pres->p_res_transform->at_ff_params[aptr];
					if (pat_par_mut)
					{
						if (pat_par_mut->mass > 0.5) mass = pat_par_mut->mass;
					}
				}
				aptr->SetFFSymbol("DU");
				aptr->SetCharge(0);
				aptr->SetMass(mass);
			}
		} // cycle on Atoms -DA dummy rules 

		if (p_res_ff_templ)  // Add improper Angles 
		{
			int n_impr = p_res_ff_templ->improper_dihedrals.size();
			int i;
			for (i = 0; i < n_impr; i++)
			{
				std::string at_ref = p_res_ff_templ->improper_dihedrals[i][0];
				HaAtom* aptr1 = pres->GetAtomByName(at_ref);
				if (aptr1 == NULL && ha_debug_level > 1) PrintLog(" Can't map atom %s \n", at_ref);

				at_ref = p_res_ff_templ->improper_dihedrals[i][1];
				HaAtom* aptr2 = pres->GetAtomByName(at_ref);
				if (aptr2 == NULL && ha_debug_level > 1) PrintLog(" Can't map atom %s \n", at_ref);;

				at_ref = p_res_ff_templ->improper_dihedrals[i][2];
				HaAtom* aptr3 = pres->GetAtomByName(at_ref);
				if (aptr3 == NULL && ha_debug_level > 1) PrintLog(" Can't map atom %s\n", at_ref);

				at_ref = p_res_ff_templ->improper_dihedrals[i][3];
				HaAtom* aptr4 = pres->GetAtomByName(at_ref);
				if (aptr4 == NULL && ha_debug_level > 1) PrintLog(" Can't map atom %s\n", at_ref);

				if (aptr1 && aptr2 && aptr3 && aptr4)
				{
					AddImprDihedral(aptr1, aptr2, aptr3, aptr4, false);
				}
				else
				{
					if( ha_debug_level > 1 ) PrintLog(" Warning: failed to set an improper angle from the residue FF template %s \n", res_fname);
				}
			}
		}

		if (pres->IsAlchemicalTransformationSet() && p_res_mut_ff_templ)
		{
			const int n_impr = p_res_mut_ff_templ->improper_dihedrals.size();
			int i;
			for (i = 0; i < n_impr; i++)
			{
				if (p_res_ff_templ->improper_dihedrals[i].size() != 4)
				{
					PrintLog(" Improper angle templates have number of atoms not equal 4 \n");
					continue;
				}
				std::string at_ref = p_res_mut_ff_templ->improper_dihedrals[i][0];
				HaAtom* aptr1 = pres->GetAtomByName(at_ref);
				if (aptr1 == NULL) PrintLog(" Can't map atom %s \n", at_ref.c_str());

				at_ref = p_res_mut_ff_templ->improper_dihedrals[i][1];
				HaAtom* aptr2 = pres->GetAtomByName(at_ref);
				if (aptr2 == NULL) PrintLog(" Can't map atom %s \n", at_ref.c_str());;

				at_ref = p_res_mut_ff_templ->improper_dihedrals[i][2];
				HaAtom* aptr3 = pres->GetAtomByName(at_ref);
				if (aptr3 == NULL) PrintLog(" Can't map atom %s\n", at_ref.c_str());

				at_ref = p_res_mut_ff_templ->improper_dihedrals[i][3];
				HaAtom* aptr4 = pres->GetAtomByName(at_ref);
				if (aptr4 == NULL) PrintLog(" Can't map atom %s\n", at_ref.c_str());

				if (aptr1 && aptr2 && aptr3 && aptr4)
				{
					AddImprDihedral(aptr1, aptr2, aptr3, aptr4, true);
				}
				else
				{
					PrintLog(" Warning: failed to set and improper angle from mutated residue FF template %s \n", res_mut_fname.c_str());
				}
			}
		}
	} // Cycle on residues


	for (HaAtom* aptr : Atoms)
	{
		std::string ff_s = aptr->GetFFSymbol();
		if (ff_s.empty())
		{
			PrintLog("No FF symbol was set for Atom %s \n", aptr->GetRef().c_str());
		}
	}

	if (ff_type == ForceFieldType::ARROW_5_14_CT || ff_type == ForceFieldType::ARROW_2_0)
	{
		to_init_mm_model = FALSE;
		return TRUE;
	}

	map<HaAtom*, set<HaAtom*>> conn_graph;      // Atom connection graph for the main state
	map<HaAtom*, set<HaAtom*>> conn_graph_mut;  // Atom connection graph for the mutated state

	BondIteratorMolSet bitr(pmset);
	for (HaBond* bptr = bitr.GetFirstBond(); bptr; bptr = bitr.GetNextBond())   // Set MM Bonds
	{
		std::string ff_s1 = bptr->srcatom->GetFFSymbol();
		std::string ff_s2 = bptr->dstatom->GetFFSymbol();
		boost::to_upper(ff_s1);
		boost::to_upper(ff_s2);

		if (ff_s1 == "DH" || ff_s2 == "DH") continue; // ignore dummy hydrogens
		if (ff_s1 == "DC" || ff_s2 == "DC") continue; // ignore dummy carbons
		if (ff_s1 == "DN" || ff_s2 == "DN") continue; // ignore dummy nitrogens
		if (ff_s1 == "DO" || ff_s2 == "DO") continue; // ignore dummy oxygens


		if (bptr->srcatom == bptr->dstatom)
		{
			PrintLog("A bond has the same source and destination atoms %s - %s ", bptr->srcatom->GetRef().c_str(), bptr->dstatom->GetRef().c_str());
			continue;
		}
		if (conn_graph[bptr->srcatom].count(bptr->dstatom) > 0)
		{
			PrintLog("MolMechModel::InitModel() Atoms %s %s are already connected \n", bptr->srcatom->GetRef().c_str(), bptr->dstatom->GetRef().c_str());
		}

		SetMMBond(bptr->srcatom, bptr->dstatom, 0.0, 0.0, NOT_SET);

		conn_graph[bptr->srcatom].insert(bptr->dstatom);
		conn_graph[bptr->dstatom].insert(bptr->srcatom);
	}

	for (HaResidue* pres : ritr)   // Set MM Bonds for the Mutated state
	{
		for (HaAtom* aptr : *pres)
		{
			if (pres->IsAlchemicalTransformationSet())
			{
				for (HaBond& bnd : pres->p_res_transform->bonds_b)
				{
					SetMMBond(bnd.srcatom, bnd.dstatom, 0.0, 0.0, NOT_SET, true);

					conn_graph_mut[bnd.srcatom].insert(bnd.dstatom);
					conn_graph_mut[bnd.dstatom].insert(bnd.srcatom);
				}
			}
		}
	}

	// FAST WATER:
	if (ff_type == ForceFieldType::AMOEBA || ff_type == ForceFieldType::AMBER_10 || ff_type == ForceFieldType::AMBER_03 || ff_type == ForceFieldType::AMBER_99_BSC0 ||
		ff_type == ForceFieldType::AMBER_99_SB)
	{
		for (HaResidue* pres : ritr )
		{
			if (pres->IsWater())
			{
				HaAtom* ph1 = nullptr;
				HaAtom* ph2 = nullptr;
				for (HaAtom* aptr : *pres )
				{
					if (aptr->IsHydrogen())
					{
						if (ph1 == nullptr)
						{
							ph1 = aptr;
							continue;
						}
						if (ph2 == nullptr)
						{
							ph2 = aptr;
							continue;
						}
					}
				}
				if (ph1 != nullptr & ph2 != nullptr)
				{
					if (conn_graph[ph1].count(ph2) > 0) continue; // H1 - H2 bond in water in already present  
					SetMMBond(ph1, ph2, 0.0, 0.0, NOT_SET);
					conn_graph[ph1].insert(ph2);
					conn_graph[ph2].insert(ph1);
				}
			}
		}
	}

	for(HaAtom* pt: Atoms) // Fill Arrays of Valence Angles and Dihedrals
	{
		set<HaAtom*>::iterator itr1;
		for( itr1 = conn_graph[pt].begin(); itr1 != conn_graph[pt].end(); itr1++) 
		{
			HaAtom* pt1 = *itr1;
			set<HaAtom*>::iterator itr2;
			for( itr2 = conn_graph[pt].begin() ; itr2 != itr1; itr2++)
			{
				HaAtom* pt2 = *itr2;
				SetValAngle(pt1,pt,pt2,0.0,0.0,NOT_SET);
				if( pt2 > pt )
				{
					for(HaAtom* pt3 : conn_graph[pt2] )
					{
						if( pt3 == pt || pt3 == pt1 ) continue;
						Dihedrals.push_back(std::make_shared<MMDihedral>(pt1,pt,pt2,pt3));
					}
				}
				if(pt1 > pt)
				{
					for(HaAtom* pt3 : conn_graph[pt1])
					{
						if( pt3 == pt || pt3 == pt2 ) continue;
						Dihedrals.push_back(std::make_shared<MMDihedral>(pt3,pt1,pt,pt2));
					}
				}
			}
		}
	}

	for (HaAtom* pt : Atoms) // Fill Arrays of Valence Angles and Dihedrals for the mutated state
	{
		set<HaAtom*>::iterator itr1;
		for (itr1 = conn_graph_mut[pt].begin(); itr1 != conn_graph_mut[pt].end(); itr1++)
		{
			HaAtom* pt1 = *itr1;
			set<HaAtom*>::iterator itr2;
			for (itr2 = conn_graph_mut[pt].begin(); itr2 != itr1; itr2++)
			{
				HaAtom* pt2 = *itr2;
				SetValAngle(pt1, pt, pt2, 0.0, 0.0, NOT_SET, true);
				if (pt2 > pt)
				{
					for (HaAtom* pt3 : conn_graph_mut[pt2])
					{
						if (pt3 == pt || pt3 == pt1) continue;
						Dihedrals_mut.push_back(std::make_shared<MMDihedral>(pt1, pt, pt2, pt3));
					}
				}
				if (pt1 > pt)
				{
					for (HaAtom* pt3 : conn_graph_mut[pt1])
					{
						if (pt3 == pt || pt3 == pt2) continue;
						Dihedrals_mut.push_back(std::make_shared<MMDihedral>(pt3, pt1, pt, pt2));
					}
				}
			}
		}
	}

	Set14interDihFlags();
	BuildExcludedAtomList();
//	if(build_nb_contact_list_flag)
//	{
//		PrintLog("Build Non-bonded contacts list \n");
//		BuildNonBondContactList();
//	}

	SetStdVdWParams();
	SetStdValParams();

	to_init_mm_model = FALSE;

	p_amber_model->UpdateAmberData();

	return True;
}

int MolMechModel::UpdateModel()
{
	int ires = TRUE;
	if( to_init_mm_model ) ires = InitModel();
	if( !ires ) return FALSE;
	if( p_amber_model->to_update_amber_data ) ires = p_amber_model->UpdateAmberData();
	return ires;
}



int MolMechModel::Clear()
{
	Atoms.clear();
	MBonds.clear();
//	ValAngles.clear();
	Dihedrals.clear();
	
//	ImprDihedrals.clear();

	DistConstraints.clear();

	excluded_atom_list.clear();
	nonbond_contact_list.clear();
	
	to_init_mm_model     = TRUE;
	p_amber_model->Clear();

	ClearMortModel();

	return TRUE;
}

AtomIntMap& MolMechModel::GetAtIdxMap(int recalc)
{
	if( recalc )
	{
		at_idx_map.clear();
		int na = Atoms.size();
		int i = 0;
		for( i = 0; i < na; i++)
		{
			HaAtom* aptr = Atoms[i];
			at_idx_map[aptr] = i;
		}
	}
	return at_idx_map;
}

AtomFFParam* MolMechModel::GetAtomFFParamMut(HaAtom* aptr)
{
	if (!aptr) return nullptr;
	HaResidue* pres = aptr->GetHostRes();
	
	if (!pres->IsAlchemicalTransformationSet()) return nullptr;
	
	if (pres->p_res_transform->at_ff_params.count(aptr) > 0)
		return pres->p_res_transform->at_ff_params[aptr].get();
	
	return nullptr;
}


MMBond* MolMechModel::GetMMBond(HaAtom* pt1, HaAtom* pt2, bool mutated_state)
{
	if( pt1 == NULL || pt2 == NULL) return NULL;

    MMBond bnd;
	
	if( pt1 < pt2)
	{
		bnd.pt1 = pt1;
		bnd.pt2 = pt2;
	}
	else
	{
		bnd.pt1 = pt1;
		bnd.pt2 = pt2;
	}

	set<MMBond>::iterator itr;
	
	if (mutated_state)
	{
		itr = MBonds_mut.find(bnd);
		if(itr == MBonds_mut.end()) return nullptr;
	}
	else
	{
		itr = MBonds.find(bnd);
		if(itr == MBonds.end()) return nullptr;
	}
	return (MMBond*)&(*itr);
}

MMValAngle* MolMechModel::GetValAngle(HaAtom* pt1, HaAtom* pt2, HaAtom* pt3, bool mutated_state)
{
	if( pt1 == NULL || pt2 == NULL || pt3 == NULL) return NULL;

	MMValAngle va;

	if( pt1 < pt3)
	{
		va.pt1 = pt1;
		va.pt2 = pt2;
		va.pt3 = pt3;
	}
	else
	{
		va.pt1 = pt3;
		va.pt2 = pt2;
		va.pt3 = pt1;
	}

	set<MMValAngle>::iterator itr;
	if (mutated_state)
	{
		itr = ValAngles_mut.find(va);
		if (itr == ValAngles_mut.end()) return nullptr;
	}
	else
	{
		itr = ValAngles.find(va);
		if (itr == ValAngles.end()) return nullptr;
	}
	return (MMValAngle*)&(*itr);
}

std::shared_ptr<MMDihedral> MolMechModel::GetDihedral(HaAtom* pt1, HaAtom* pt2, HaAtom* pt3, HaAtom* pt4, bool mutated_state )
{
	if( pt1 == nullptr || pt2 == nullptr || pt3 == nullptr || pt4 == nullptr) return nullptr;

	auto itr = Dihedrals.begin();
	auto itr_end = Dihedrals.end();
	if (mutated_state)
	{
		itr = Dihedrals_mut.begin();
		itr_end = Dihedrals_mut.end();
	}

	for(; itr != itr_end; itr++)
	{
		MMDihedral& dih = *(*itr);
		
		if( (dih.pt1 == pt1 && dih.pt2 == pt2 && dih.pt3 == pt3 && dih.pt4 == pt4) || 
			(dih.pt1 == pt4 && dih.pt3 == pt2 && dih.pt3 == pt2 && dih.pt4 == pt1) )
		{
			return *itr;
		}
	}
	return nullptr;
}

std::shared_ptr<MMDihedral> MolMechModel::GetImprDihedral(HaAtom* pt1, HaAtom* pt2, HaAtom* pt3, HaAtom* pt4, bool mutate_state ) 
{
	if( pt1 == nullptr || pt2 == nullptr || pt3 == nullptr || pt4 == nullptr) return nullptr;

	auto itr = ImprDihedrals.begin();
	auto itr_end = ImprDihedrals.end();
	if (mutate_state)
	{
		itr = ImprDihedrals_mut.begin();
		itr_end = ImprDihedrals_mut.end();
	}

	for( ; itr != itr_end; itr++)
	{
		MMDihedral& dih = *(*itr);
		
		if( dih.pt1 == pt1 && dih.pt2 == pt2 && dih.pt3 == pt3 && dih.pt4 == pt4 )
		{
			return *itr;
		}
	}
	return nullptr;
}

std::shared_ptr<MMDihedral> MolMechModel::AddImprDihedral(HaAtom* pt1, HaAtom* pt2, HaAtom* pt3, HaAtom* pt4, bool mutated_state )
{
	if( pt1 == nullptr || pt2 == nullptr || pt3 == nullptr || pt4 == nullptr)
	{
		ErrorInMod("MolMechModel::AddImprDihedral()",
			       " One of the pointers to HaAtom is NULL");
		return nullptr;
	}

	std::shared_ptr<MMDihedral> ps_dih;
	if (mutated_state)
	{
		ImprDihedrals_mut.push_back(std::make_shared<MMDihedral>(pt1, pt2, pt3, pt4, true));
		ps_dih = ImprDihedrals_mut.back();
	}
	else
	{
		ImprDihedrals.push_back(std::make_shared<MMDihedral>(pt1, pt2, pt3, pt4, true));
		ps_dih = ImprDihedrals.back();
	}
	
	return( ps_dih );
}

int MolMechModel::SetMMBond(HaAtom* pt1, HaAtom* pt2, double r0, double fc, int set_type, bool mutated_state )
{
	if( pt1 == NULL || pt2 == NULL  )
	{
		PrintLog("Error in MolMechModel::SetMMBond() \n");
		PrintLog(" One of the atom pointers is NULL \n");
		return FALSE;
	}
	
	HaAtom* pt1_loc = pt1;
	HaAtom* pt2_loc = pt2;

	if(pt2 < pt1)
	{
		pt1_loc  = pt2;
		pt2_loc  = pt1;
	}

	MMBond* pbnd = GetMMBond(pt1_loc,pt2_loc);
	
	if(pbnd == NULL)
	{
		MMBond& bnd = (MMBond&) *((MBonds.insert(MMBond(pt1_loc,pt2_loc))).first);
		pbnd = &bnd;
	}
	
	if( pbnd->set_type <= set_type)
	{
		pbnd->r0 = r0;
		pbnd->fc = fc;
		pbnd->set_type = set_type;
	}

	p_amber_model->SetUpdateDataFlag();

	return TRUE;
}

int MolMechModel::SetValAngle(HaAtom* pt1, HaAtom* pt2, HaAtom* pt3, double a0, double fc, int set_type, bool mutated_state)
{
	if( pt1 == NULL || pt2 == NULL || pt3 == NULL )
	{
		PrintLog("Error in HaMolMechMod::SetValAngle() \n");
		PrintLog(" One of atom pointers is NULL \n");
		return FALSE;
	}

	if( pt1 == pt2 || pt1 == pt3 || pt2 == pt3)
	{
		PrintLog(" Error in: HaMolMechMod::SetValAngle(): invalid valence angle \n");
		PrintLog(" pt1 == pt2 || pt1 == pt3 || pt2 == pt3 \n");
		return FALSE;
	}
	
	HaAtom* pt1_loc = pt1;
	HaAtom* pt2_loc = pt2;
	HaAtom* pt3_loc = pt3;

	if(pt3 < pt1)
	{
		pt1_loc  = pt3;
		pt3_loc = pt1;
	}

	MMValAngle* pang = GetValAngle(pt1_loc,pt2_loc,pt3_loc);
	
	if(pang == NULL)
	{
		MMValAngle& va = (MMValAngle&) *((ValAngles.insert(MMValAngle(pt1_loc,pt2_loc,pt3_loc))).first);
		pang = &va;
	}
	
	if( pang->set_type <= set_type)
	{
		pang->a0 = a0;
		pang->fc = fc;
		pang->set_type = set_type;
	}

	p_amber_model->SetUpdateDataFlag();

	return TRUE;
}

int MolMechModel::Set14interDihFlags()
{
	AtomAtomMultiMap forbid_14_map;
	AtomAtomMultiMap forbid_14_map_mut;
	
	for(const MMBond& bnd : MBonds)
	{
		HaAtom* pt1 = bnd.pt1;
		HaAtom* pt2 = bnd.pt2;
		AtomAtomMultiMap::value_type p(pt1,pt2);
		forbid_14_map.insert(p);
		AtomAtomMultiMap::value_type p2(pt2,pt1);
		forbid_14_map.insert(p2);
	}

	for (const MMBond& bnd : MBonds_mut)
	{
		HaAtom* pt1 = bnd.pt1;
		HaAtom* pt2 = bnd.pt2;
		AtomAtomMultiMap::value_type p(pt1, pt2);
		forbid_14_map_mut.insert(p);
		AtomAtomMultiMap::value_type p2(pt2, pt1);
		forbid_14_map_mut.insert(p2);
	}

	for( const MMValAngle& va : ValAngles)
	{
		HaAtom* pt1 = va.pt1;
		HaAtom* pt3 = va.pt3;
		AtomAtomMultiMap::value_type p(pt1,pt3);
		forbid_14_map.insert(p);
		AtomAtomMultiMap::value_type p2(pt3,pt1);
		forbid_14_map.insert(p2);
	}

	for (const MMValAngle& va : ValAngles_mut)
	{
		HaAtom* pt1 = va.pt1;
		HaAtom* pt3 = va.pt3;
		AtomAtomMultiMap::value_type p(pt1, pt3);
		forbid_14_map_mut.insert(p);
		AtomAtomMultiMap::value_type p2(pt3, pt1);
		forbid_14_map_mut.insert(p2);
	}

	for( shared_ptr<MMDihedral> sp_dih : Dihedrals)
	{
		HaAtom* pt1 = sp_dih->pt1;
		HaAtom* pt4 = sp_dih->pt4;

		AtomAtomMultiMap::iterator ditr1, ditr2;
		ditr1 = forbid_14_map.lower_bound(pt1);
		ditr2 = forbid_14_map.upper_bound(pt1);

		bool do_14 = true;
		for(; ditr1 != ditr2; ditr1++)
		{
			if( (*ditr1).second == pt4)
			{
				do_14 = false;
				break;
			}
		}
		if( do_14 ) // only one dihedral with given 1-4 atom pair will have calc_14 flag set to true
		{
			sp_dih->calc_14 = true;
			AtomAtomMultiMap::value_type p(pt1,pt4);
			forbid_14_map.insert(p);
			AtomAtomMultiMap::value_type p2(pt4,pt1);
			forbid_14_map.insert(p2);
		}
		else
		{
			sp_dih->calc_14 = false;
		}
	}

	for (shared_ptr<MMDihedral> sp_dih : Dihedrals_mut)
	{
		HaAtom* pt1 = sp_dih->pt1;
		HaAtom* pt4 = sp_dih->pt4;

		AtomAtomMultiMap::iterator ditr1, ditr2;
		ditr1 = forbid_14_map_mut.lower_bound(pt1);
		ditr2 = forbid_14_map_mut.upper_bound(pt1);

		bool do_14 = true;
		for (; ditr1 != ditr2; ditr1++)
		{
			if ((*ditr1).second == pt4)
			{
				do_14 = false;
				break;
			}
		}
		if (do_14) // only one dihedral with given 1-4 atom pair will have calc_14 flag set to true
		{
			sp_dih->calc_14 = true;
			AtomAtomMultiMap::value_type p(pt1, pt4);
			forbid_14_map_mut.insert(p);
			AtomAtomMultiMap::value_type p2(pt4, pt1);
			forbid_14_map_mut.insert(p2);
		}
		else
		{
			sp_dih->calc_14 = false;
		}
	}
	return True;
}

int MolMechModel::SetCoarseGrainedOPEPParams()  // jose 11/04/2008 under construction
{
	UpdateModel();

	vector<HaAtom*>::iterator aitr;
	for (aitr = Atoms.begin(); aitr != Atoms.end(); aitr++)
	{
		HaAtom* aptr = (*aitr);
		HaResidue* rptr = aptr->GetHostRes();
		//std::string res_name = rptr->GetName();
		std::string res_mod = rptr->GetNameModifier();

		if (res_mod != "CG" ) continue; // under decision

		std::string at_name = aptr->GetName();
		if (at_name == "X")
		{
			aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
			aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
		}
	}
	
	PrintLog("Set OPEP Bonds params \n");

	set<MMBond, less<MMBond> >::iterator bitr;
	for( bitr = MBonds.begin(); bitr != MBonds.end(); bitr++)
	{
        MMBond* bptr = (MMBond*)&(*bitr);
		HaAtom* aptr1 = bptr->pt1; 
		HaAtom* aptr2 = bptr->pt2;

		HaResidue* pres1 = aptr1->GetHostRes();
		HaResidue* pres2 = aptr2->GetHostRes();

		//std::string res_name_1 = pres1->GetFullName();
		//std::string res_name_2 = pres2->GetFullName();

		std::string res_mod1 = pres1->GetNameModifier();
		std::string res_mod2 = pres2->GetNameModifier();

		if (res_mod1 != "CG" || res_mod2 != "CG")	continue;

		std::string at_name1 = aptr1->GetName();
		std::string at_name2 = aptr2->GetName();

		if (( at_name1 == "X" && at_name2 == "CA") || ( at_name1 == "CA" && at_name2 == "X"))
		{
			bptr->set_type = SET_FF_FIELD;
			bptr->fc = 400.0;	//< Kcal/mol*A^2 Derreumaux JCP, 1999 
			bptr->r0 = 0.0;		//< future modification from OPEP parameter for each residue
		}
	}

	PrintLog("Set OPEP Valence params \n");

	set<MMValAngle, less<MMValAngle> >::iterator vaitr;
	for( vaitr = ValAngles.begin(); vaitr != ValAngles.end(); vaitr++)
	{
        MMValAngle* pang = (MMValAngle*) &(*vaitr);
		HaAtom* aptr1 = pang->pt1; 
		HaAtom* aptr2 = pang->pt2;
		HaAtom* aptr3 = pang->pt3;

		HaResidue* pres1 = aptr1->GetHostRes();
		HaResidue* pres2 = aptr2->GetHostRes();
		HaResidue* pres3 = aptr3->GetHostRes();

		//std::string res_name_1 = pres1->GetFullName();
		//std::string res_name_2 = pres2->GetFullName();
		//std::string res_name_3 = pres3->GetFullName();

		std::string res_mod1 = pres1->GetNameModifier();
		std::string res_mod2 = pres2->GetNameModifier();
		std::string res_mod3 = pres3->GetNameModifier();

		if (res_mod1 != "CG" || res_mod2 != "CG" || res_mod3 != "CG")	continue;

		std::string at_name1 = aptr1->GetName();
		std::string at_name2 = aptr2->GetName();
		std::string at_name3 = aptr3->GetName();

		if ((at_name1 == "X" && at_name2 == "CA" && at_name3 == "N") || (at_name1 == "N" && at_name2 == "CA" && at_name3 == "X"))
		{
			pang->set_type = SET_FF_FIELD;
			pang->fc = 12.0;	//< Kcal/mol*A^2 Derreumaux JCP, 1999
			pang->a0 = 0.0;		//< Degrees, OPEP parameter for each residue
		}
	}

	PrintLog("Set OPEP Dihedral Angles params \n");

	for( shared_ptr<MMDihedral> daitr : Dihedrals)
	{
		HaAtom* aptr1 = (*daitr).pt1; 
		HaAtom* aptr2 = (*daitr).pt2;
		HaAtom* aptr3 = (*daitr).pt3;
		HaAtom* aptr4 = (*daitr).pt4;

		HaResidue* pres1 = aptr1->GetHostRes();
		HaResidue* pres2 = aptr2->GetHostRes();
		HaResidue* pres3 = aptr3->GetHostRes();
		HaResidue* pres4 = aptr4->GetHostRes();

		//std::string res_name_1 = pres1->GetFullName();
		//std::string res_name_2 = pres2->GetFullName();
		//std::string res_name_3 = pres3->GetFullName();
		//std::string res_name_4 = pres4->GetFullName();

		std::string res_mod1 = pres1->GetNameModifier();
		std::string res_mod2 = pres2->GetNameModifier();
		std::string res_mod3 = pres3->GetNameModifier();
		std::string res_mod4 = pres4->GetNameModifier();

		if (res_mod1 != "CG" || res_mod2 != "CG" || res_mod3 != "CG" || res_mod4 != "CG")	continue;

		std::string at_name1 = aptr1->GetName();
		std::string at_name2 = aptr2->GetName();
		std::string at_name3 = aptr3->GetName();
		std::string at_name4 = aptr4->GetName();

		if ((at_name1 == "X" && at_name2 == "N" && at_name3 == "CA" && at_name4 == "C") ||
			(at_name1 == "C" && at_name2 == "CA" && at_name3 == "N" && at_name4 == "X"))
		{
			(*daitr).ClearParams();
			(*daitr).set_type = SET_FF_FIELD;
			(*daitr).pn.push_back(1.0);
			(*daitr).pk.push_back(3);		// Kcal/mol*rad^2 OPEP parameter Derreumaux JCP, 1999
			(*daitr).idivf.push_back(1.0);
			(*daitr).phase.push_back(0.0);	// OPEP geometric parameter parameter for each residue
		}
	}

	PrintLog("Set OPEP nonbonded params \n");

	// future coding

	PrintLog("Set OPEP hydrogen-bonded params \n");

	p_amber_model->SetUpdateDataFlag();

	return TRUE;
}

int MolMechModel::SetCoarseGrainedDNAParams()
{
	UpdateModel();

	double eps = 0.26; // epsilon energy parameter of the De Pablo force field
	double dcut = 6.86; // d_cut - distance to start rupulsive VdW interations in De pablo force field

	// Set MolMechMod params corresponding to De Pablo Coarse Grained DNA model  

	SetIonStrength(0.15);
	SetDielConst(80.0);
	electr_method = electr_method.SCREENED_COULOMB;
	SetNBCutDist(30.0);
	SetScale14VdW(1000000.0);

    vector<HaAtom*>::iterator aitr;
	for( aitr = Atoms.begin(); aitr != Atoms.end(); aitr++)
	{
		HaAtom* aptr    = (*aitr);
		HaResidue* pres = aptr->GetHostRes();
		std::string res_name = pres->GetFullName();

		if( res_name != "A#CG" && res_name != "T#CG" 
			&& res_name != "G#CG" &&  res_name != "C#CG")
		{
			continue; // residue of the atom is not course grained nucleotide
		}

		std::string at_name = aptr->GetName();
 		if( at_name == "SG")
		{
			aptr->SetVdWRad(0.5*dcut);
			aptr->SetVdWEne(eps);
			aptr->SetFFSymbol("SG");
		}
 		if( at_name == "PH")
		{
			aptr->SetVdWRad(0.5*dcut);
			aptr->SetVdWEne(eps);
			aptr->SetFFSymbol("PH");
		}
		if( at_name == "AB" )
		{
			aptr->SetVdWRad(0.5 * dcut);
			aptr->SetVdWEne(eps);
			aptr->SetFFSymbol("AB");
		}
		if( at_name == "TB" )
		{
			aptr->SetVdWRad(0.5 * dcut);
			aptr->SetVdWEne(eps);
			aptr->SetFFSymbol("TB");
		}
		if( at_name == "GB" )
		{
			aptr->SetVdWRad(0.5 * dcut);
			aptr->SetVdWEne(eps);
			aptr->SetFFSymbol("GB");
		}
		if( at_name == "CB" )
		{
			aptr->SetVdWRad(0.5 * dcut);
			aptr->SetVdWEne(eps);
			aptr->SetFFSymbol("CB");
		}
	}

	PrintLog(" Set DNA Bonds params \n");

	set<MMBond, less<MMBond> >::iterator bitr;
	for( bitr = MBonds.begin(); bitr != MBonds.end(); bitr++)
	{
        MMBond* bptr = (MMBond*)&(*bitr);
		HaAtom* aptr1 = bptr->pt1; 
		HaAtom* aptr2 = bptr->pt2;

		HaResidue* pres1 = aptr1->GetHostRes();
		HaResidue* pres2 = aptr2->GetHostRes();

		std::string res_name_1 = pres1->GetFullName();
		std::string res_name_2 = pres2->GetFullName();

		if( res_name_1 != "A#CG" && res_name_1 != "T#CG" 
			&& res_name_1 != "G#CG" &&  res_name_1 != "C#CG")
		{
			continue; // residue of the first atom is not course grained nucleotide
		}

		if( res_name_2 != "A#CG" && res_name_2 != "T#CG" 
			&& res_name_2 != "G#CG" &&  res_name_2 != "C#CG")
		{
			continue; // residue of the second atom is not course grained nucleotide
		}

		std::string at_name_1 = aptr1->GetName();
		std::string at_name_2 = aptr2->GetName();

		if( (at_name_1 == "SG" && at_name_2 == "PH") || 
			(at_name_1 == "PH" && at_name_2 == "SG") )
		{
                        
			bptr->set_type = SET_FF_FIELD;
			bptr->fc = 100.0*eps;
			if( pres1 == pres2 )
			{
				bptr->r0 = 3.899;
			}
			else
			{
				bptr->r0 = 3.559;
			}
		}

		if( (at_name_1 == "SG" && at_name_2 == "AB") ||
			(at_name_1 == "AB" && at_name_2 == "SG"))
		{
			bptr->set_type = SET_FF_FIELD;
			bptr->fc = 100.0*eps;
			bptr->r0 = 6.43;
		}

		if( (at_name_1 == "SG" && at_name_2 == "TB") ||
			(at_name_1 == "TB" && at_name_2 == "SG"))
		{
			bptr->set_type = SET_FF_FIELD;
			bptr->fc = 100.0*eps;
			bptr->r0 = 4.88;
		}

		if( (at_name_1 == "SG" && at_name_2 == "CB") ||
			(at_name_1 == "CB" && at_name_2 == "SG"))
		{
			bptr->set_type = SET_FF_FIELD;
			bptr->fc = 100.0*eps;
		        bptr->r0 = 4.921;
		}

		if( (at_name_1 == "SG" && at_name_2 == "GB") ||
			(at_name_1 == "GB" && at_name_2 == "SG"))
		{
			bptr->set_type = SET_FF_FIELD;
			bptr->fc = 100.0*eps;
			bptr->r0 = 6.392;
		}
	}

	PrintLog(" Set DNA Valence Angles params \n");

	set<MMValAngle, less<MMValAngle> >::iterator vaitr;
	for( vaitr = ValAngles.begin(); vaitr != ValAngles.end(); vaitr++)
	{
        MMValAngle* pang = (MMValAngle*) &(*vaitr);
		HaAtom* aptr1 = pang->pt1; 
		HaAtom* aptr2 = pang->pt2;
		HaAtom* aptr3 = pang->pt3;

		HaResidue* pres1 = aptr1->GetHostRes();
		HaResidue* pres2 = aptr2->GetHostRes();
		HaResidue* pres3 = aptr3->GetHostRes();

		std::string res_name_1 = pres1->GetFullName();
		std::string res_name_2 = pres2->GetFullName();
		std::string res_name_3 = pres3->GetFullName();


		if( res_name_1 != "A#CG" && res_name_1 != "T#CG" 
			&& res_name_1 != "G#CG" &&  res_name_1 != "C#CG")
		{
			continue; // residue of the first atom is not course grained nucleotide
		}

		if( res_name_2 != "A#CG" && res_name_2 != "T#CG" 
			&& res_name_2 != "G#CG" &&  res_name_2 != "C#CG")
		{
			continue; // residue of the second atom is not course grained nucleotide
		}

		if( res_name_3 != "A#CG" && res_name_3 != "T#CG" 
			&& res_name_3 != "G#CG" &&  res_name_3 != "C#CG")
		{
			continue; // residue of the third atom is not course grained nucleotide
		}

		std::string at_name_1 = aptr1->GetName();
		std::string at_name_2 = aptr2->GetName();
		std::string at_name_3 = aptr3->GetName();

		if( at_name_1 == "SG" && at_name_2 == "PH" && at_name_3 == "SG" )
		{
			pang->set_type = SET_FF_FIELD;
			pang->fc = 400.0*eps;
			pang->a0 = 94.49;
		}

		if( at_name_1 == "PH" && at_name_2 == "SG" && at_name_3 == "PH" )
		{
			pang->set_type = SET_FF_FIELD;
			pang->fc = 400.0*eps;
			pang->a0 = 120.15;
		}

		if( (at_name_1 == "PH" && at_name_2 == "SG" && at_name_3 == "AB") ||
			(at_name_1 == "AB" && at_name_2 == "SG" && at_name_3 == "PH") )
		{
			pang->set_type = SET_FF_FIELD;
			pang->fc = 400.0*eps;
			if( pres1 == pres3 )
			{
				pang->a0 = 113.13;
			}
			else
			{
				pang->a0 = 108.38;
			}
		}

		if( (at_name_1 == "PH" && at_name_2 == "SG" && at_name_3 == "TB") ||
			(at_name_1 == "TB" && at_name_2 == "SG" && at_name_3 == "PH") )
		{
			pang->set_type = SET_FF_FIELD;
			pang->fc = 400.0*eps;
			if( pres1 == pres3 )
			{
				pang->a0 = 102.79;
			}
			else
			{
				pang->a0 = 112.72;
			}
		}

		if( (at_name_1 == "PH" && at_name_2 == "SG" && at_name_3 == "CB") ||
			(at_name_1 == "CB" && at_name_2 == "SG" && at_name_3 == "PH") )
		{
			pang->set_type = SET_FF_FIELD;
			pang->fc = 400.0*eps;
			if( pres1 == pres3 )
			{
				pang->a0 = 103.49;
			}
			else
			{
				pang->a0 = 112.39;
			}
		}

		if( (at_name_1 == "PH" && at_name_2 == "SG" && at_name_3 == "GB") ||
			(at_name_1 == "GB" && at_name_2 == "SG" && at_name_3 == "PH") )
		{
			pang->set_type = SET_FF_FIELD;
			pang->fc = 400.0*eps;
			if( pres1 == pres3 )
			{
				pang->a0 = 113.52;
			}
			else
			{
				pang->a0 = 108.12;
			}
		}
	}

	PrintLog(" Set DNA Dihedral Angles params \n");

	for( shared_ptr<MMDihedral> daitr : Dihedrals )
	{
		HaAtom* aptr1 = (*daitr).pt1; 
		HaAtom* aptr2 = (*daitr).pt2;
		HaAtom* aptr3 = (*daitr).pt3;
		HaAtom* aptr4 = (*daitr).pt4;

		HaResidue* pres1 = aptr1->GetHostRes();
		HaResidue* pres2 = aptr2->GetHostRes();
		HaResidue* pres3 = aptr3->GetHostRes();
		HaResidue* pres4 = aptr4->GetHostRes();

		std::string res_name_1 = pres1->GetFullName();
		std::string res_name_2 = pres2->GetFullName();
		std::string res_name_3 = pres3->GetFullName();
		std::string res_name_4 = pres4->GetFullName();


		if( res_name_1 != "A#CG" && res_name_1 != "T#CG" 
			&& res_name_1 != "G#CG" &&  res_name_1 != "C#CG")
		{
			continue; // residue of the first atom is not a course grained nucleotide
		}

		if( res_name_2 != "A#CG" && res_name_2 != "T#CG" 
			&& res_name_2 != "G#CG" &&  res_name_2 != "C#CG")
		{
			continue; // residue of the second atom is not a course grained nucleotide
		}

		if( res_name_3 != "A#CG" && res_name_3 != "T#CG" 
			&& res_name_3 != "G#CG" &&  res_name_3 != "C#CG")
		{
			continue; // residue of the third atom is not a course grained nucleotide
		}

		if( res_name_4 != "A#CG" && res_name_4 != "T#CG" 
			&& res_name_4 != "G#CG" &&  res_name_4 != "C#CG")
		{
			continue; // residue of the fourth atom is not a course grained nucleotide
		}


		std::string at_name_1 = aptr1->GetName();
		std::string at_name_2 = aptr2->GetName();
		std::string at_name_3 = aptr3->GetName();
		std::string at_name_4 = aptr4->GetName();

		double eps_tors = -4.0*eps;

		if( (at_name_1 == "PH" && at_name_2 == "SG" && at_name_3 == "PH" && at_name_4 == "SG")  ||
			(at_name_1 == "SG" && at_name_2 == "PH" && at_name_3 == "SG" && at_name_4 == "PH") )
		{
			(*daitr).ClearParams();
			(*daitr).set_type = SET_FF_FIELD;
			(*daitr).pn.push_back(1.0);
			(*daitr).pk.push_back(eps_tors);
			(*daitr).idivf.push_back(1.0);
			if( pres2 == pres3 )
			{
				(*daitr).phase.push_back(-179.17);
			}
			else
			{
				(*daitr).phase.push_back(-154.80);
			}
		}

		if( (at_name_1 == "AB" && at_name_2 == "SG" && at_name_3 == "PH" && at_name_4 == "SG")  ||
			(at_name_1 == "SG" && at_name_2 == "PH" && at_name_3 == "SG" && at_name_4 == "AB") )
		{
			(*daitr).ClearParams();
			(*daitr).set_type = SET_FF_FIELD;
			(*daitr).pn.push_back(1.0);
			(*daitr).pk.push_back(eps_tors);
			(*daitr).idivf.push_back(1.0);
			if( pres2 == pres3 )
			{
				(*daitr).phase.push_back(50.69);
			}
			else
			{
				(*daitr).phase.push_back(-22.60);
			}
		}

		if( (at_name_1 == "TB" && at_name_2 == "SG" && at_name_3 == "PH" && at_name_4 == "SG")  ||
			(at_name_1 == "SG" && at_name_2 == "PH" && at_name_3 == "SG" && at_name_4 == "TB") )
		{
			(*daitr).ClearParams();
			(*daitr).set_type = SET_FF_FIELD;
			(*daitr).pn.push_back(1.0);
			(*daitr).pk.push_back(eps_tors);
			(*daitr).idivf.push_back(1.0);
			if( pres2 == pres3 )
			{
				(*daitr).phase.push_back(54.69);
			}
			else
			{
				(*daitr).phase.push_back(-33.42);
			}
		}

		if( (at_name_1 == "CB" && at_name_2 == "SG" && at_name_3 == "PH" && at_name_4 == "SG")  ||
			(at_name_1 == "SG" && at_name_2 == "PH" && at_name_3 == "SG" && at_name_4 == "CB") )
		{
			(*daitr).ClearParams();
			(*daitr).set_type = SET_FF_FIELD;
			(*daitr).pn.push_back(1.0);
			(*daitr).pk.push_back(eps_tors);
			(*daitr).idivf.push_back(1.0);
			if( pres2 == pres3 )
			{
				(*daitr).phase.push_back(54.50);
			}
			else
			{
				(*daitr).phase.push_back(-32.72);
			}
		}

		if( (at_name_1 == "GB" && at_name_2 == "SG" && at_name_3 == "PH" && at_name_4 == "SG")  ||
			(at_name_1 == "SG" && at_name_2 == "PH" && at_name_3 == "SG" && at_name_4 == "GB") )
		{
			(*daitr).ClearParams();
			(*daitr).set_type = SET_FF_FIELD;
			(*daitr).pn.push_back(1.0);
			(*daitr).pk.push_back(eps_tors);
			(*daitr).idivf.push_back(1.0);
			if( pres2 == pres3 )
			{
				(*daitr).phase.push_back(50.69);
			}
			else
			{
				(*daitr).phase.push_back(-22.30);
			}
		}
	}
  
    set<VecPtr> go_contacts;

	PrintLog("Set GO-type VdW contacts \n");

	HaAtom* aptr_sg_p   = NULL; //!< Sugar Pseudo Atom of the previous nucleotide 
	HaAtom* aptr_bs_p   = NULL; //!< Base Pseudo Atom of the previous nucleotide 
	HaAtom* aptr_bs_p2  = NULL; //!< Base Pseudo Atom of the previous-previous nucleotide 

	VecPtr cnt(2);

	MolSet::ResidueIterator ritr(pmset);
	HaResidue* pres;

	for( pres = ritr.GetFirstRes(); pres != NULL; pres = ritr.GetNextRes() )
	{
		std::string res_name = pres->GetFullName();
		if( res_name != "A#CG" && res_name != "T#CG" 
			&& res_name != "G#CG" &&  res_name != "C#CG")
		{
			aptr_sg_p = NULL;
			aptr_bs_p = NULL;
			aptr_bs_p2 = NULL;
			continue;
		}
		
		HaAtom* aptr_sg = pres->GetAtomByName("SG");
		HaAtom* aptr_ph = pres->GetAtomByName("PH");
		
		HaAtom* aptr_bs = NULL;
		if( res_name == "A#CG" ) aptr_bs = pres->GetAtomByName("AB");
		if( res_name == "T#CG" ) aptr_bs = pres->GetAtomByName("TB");
		if( res_name == "G#CG" ) aptr_bs = pres->GetAtomByName("GB");
		if( res_name == "C#CG" ) aptr_bs = pres->GetAtomByName("CB");

		if( aptr_bs == NULL || aptr_sg == NULL)
		{
			PrintLog("No Base or Sugar Pseudo Atoms in the nucleotide \n");
			aptr_sg_p = NULL;
			aptr_bs_p = NULL;
			aptr_bs_p2 = NULL;
			continue;
		}
	
		if( aptr_bs_p != NULL && aptr_sg_p != NULL && 
			aptr_ph != NULL && aptr_sg_p->IsBonded(*aptr_ph) )
		{
			
			AtomContact vdw_cnt(aptr_bs, aptr_bs_p, AtomContactType::VDW_CNT_6_12 );
			double dist = Vec3D::CalcDistance(aptr_bs, aptr_bs_p,ANGSTROM_U);
			vdw_cnt.SetParamsEneR(eps,dist);
			vdw_cnt.set_type = SET_FF_FIELD;
			DistConstraints.push_back(vdw_cnt);
	
			cnt[0] = aptr_bs;
			cnt[1] = aptr_bs_p;	    
			go_contacts.insert(cnt);

			std::string at_str_1 = aptr_bs->GetRef();
			std::string at_str_2 = aptr_bs_p->GetRef();
		
			PrintLog(" Add explicit VdW term between atoms %s and %s \n",
				at_str_1.c_str(), at_str_2.c_str());
		
		
			if( aptr_bs_p2 != NULL)
			{
				AtomContact vdw_cnt(aptr_bs, aptr_bs_p2, AtomContactType::VDW_CNT_6_12);
				double dist = Vec3D::CalcDistance(aptr_bs, aptr_bs_p2,ANGSTROM_U);
				vdw_cnt.SetParamsEneR(eps,dist);
				vdw_cnt.set_type = SET_FF_FIELD;
				DistConstraints.push_back(vdw_cnt);
		    			
				cnt[0] = aptr_bs;
				cnt[1] = aptr_bs_p2;
			    go_contacts.insert(cnt);

				std::string at_str_1 = aptr_bs->GetRef();
				std::string at_str_2 = aptr_bs_p2->GetRef();
			
				PrintLog(" Add explicit VdW term between atoms %s and %s \n",
					at_str_1.c_str(), at_str_2.c_str());
			}
			aptr_bs_p2 = aptr_bs_p;
		}
		else
		{
			aptr_bs_p2 = NULL;
		}
	
		aptr_bs_p = aptr_bs;
		aptr_sg_p = aptr_sg;
	}

// Hydrogen Bond Parameters for base pairing:
 
	PrintLog(" Add explicit H-Bonds terms to the force field \n");

	int i,j;
	int nat = Atoms.size();

	double cut2 = 9.0;
	cut2 = cut2*cut2;

	for( i=0; i < (nat-1); i++)
	{
		HaAtom* aptr1 = Atoms[i];
		HaResidue* pres1 = aptr1->GetHostRes();
		std::string res_name_1 = pres1->GetFullName();

		if( res_name_1 != "A#CG" && res_name_1 != "T#CG" 
			&& res_name_1 != "G#CG" &&  res_name_1 != "C#CG")
		{
			continue;  
		}
		
		std::string at_name_1 = aptr1->GetName();

		for( j= i+1; j < nat; j++)
		{
			HaAtom* aptr2 = Atoms[j];
			HaResidue* pres2 = aptr2->GetHostRes();
			std::string res_name_2 = pres2->GetFullName();

			if( res_name_2 != "A#CG" && res_name_2 != "T#CG" 
				&& res_name_2 != "G#CG" &&  res_name_2 != "C#CG")
			{
				continue;  
			}

			cnt[0] = aptr1;
			cnt[1] = aptr2;
			if(go_contacts.count(cnt) > 0) continue;
			cnt[0] = aptr2;
			cnt[1] = aptr1;
			if(go_contacts.count(cnt) > 0) continue;

			double dist = Vec3D::CalcDistanceSq(aptr1,aptr2);

			if( dist < cut2)
			{
				std::string at_name_1 = aptr1->GetName();
				std::string at_name_2 = aptr2->GetName();
                
				if( at_name_1 == "SG" || at_name_1 == "PH") continue;
				if( at_name_2 == "SG" || at_name_2 == "PH") continue;

                std::string hb_lbl;

				if( (at_name_1 == "AB" && at_name_2 == "TB") || 
                    (at_name_1 == "TB" && at_name_2 == "AB") )
				{
					hb_lbl = "AT";
				}

				if( (at_name_1 == "GB" && at_name_2 == "CB") || 
                    (at_name_1 == "CB" && at_name_2 == "GB") )
				{
					hb_lbl = "GC";
				}

				if( hb_lbl == "AT" || hb_lbl == "GC")
				{
					std::string at_str_1 = aptr1->GetRef();
					std::string at_str_2 = aptr2->GetRef();
					
					PrintLog(" Add MM hydrogen bond term between atoms %s and %s \n",
						       at_str_1.c_str(), at_str_2.c_str());

					AtomContact hbnd(aptr1,aptr2,AtomContactType::VDW_CNT_10_12);

					if( hb_lbl == "AT" )
					{
						double dist_at = 2.9002;
						hbnd.SetParamsEneR((2.0/3.0)*4.0*eps, dist_at);
					}

     				if( hb_lbl == "GC" )
					{
						double dist_gc = 2.8694;
						hbnd.SetParamsEneR(4.0*eps, dist_gc);
					}

					DistConstraints.push_back(hbnd);
				}
				else
				{
					AtomContact vdw_cnt(aptr1, aptr2, AtomContactType::VDW_CNT_6_12_NO_REP);
					vdw_cnt.SetParamsEneR(eps, 1.0);
					vdw_cnt.set_type = SET_FF_FIELD;
					DistConstraints.push_back(vdw_cnt);

					std::string at_str_1 = aptr1->GetRef();
					std::string at_str_2 = aptr2->GetRef();
		
					PrintLog(" Add explicit VdW term between atoms %s and %s \n",
							   at_str_1.c_str(), at_str_2.c_str());
				}
			}
		}
	}

	p_amber_model->SetUpdateDataFlag();
	BuildExcludedAtomList();

	return TRUE;	
}

int MolMechModel::SetCoarseGrainedAAParams()
{
	UpdateModel(); // load all parameters found on res_*.hlm files including charge, mass, ff-symbols, MBonds is created with r0,fc=0, NOT_SET, 
	double rcut = 12.0;  // cutoff radius in Angstroms
	double sigma = 4.7; // Lennard Jones effective size in Angstroms 
	return TRUE;
};

int MolMechModel::SetMoveAll()
{
	moving_atoms = "ALL_ATOMS";
    if( pmset->per_bc->IsSet() )
	{
		if( p_mm_mod->period_bcond == p_mm_mod->period_bcond.NO_PERIODICITY) 
		{
			p_mm_mod->period_bcond = p_mm_mod->period_bcond.CONST_PRES;
		}
	}
	return TRUE;
}

AtomGroup* MolMechModel::GetMovingAtoms()
{
	if(moving_atoms.length() == 0 || moving_atoms == "ALL_ATOMS")
	{
		return NULL;
	}
	AtomGroup* mv_atoms = pmset->GetAtomGroupByID(moving_atoms.c_str());
	return mv_atoms;
}

AtomGroup* MolMechModel::GetRestrAtoms()
{
	if(restrained_atoms.length() == 0 || restrained_atoms == "NO_RESTRAINTS")
	{
		return NULL;
	}	
	AtomGroup* restr_atoms = pmset->GetAtomGroupByID(restrained_atoms.c_str());
	return restr_atoms;
}

int MolMechModel::SaveAtomRestrArbalestIndForm(std::string restr_desc_fname, std::string restr_list_fname)
{
	char buf[128];
	std::string str_tmp;
	AtomGroup* restr_atoms = GetRestrAtoms();
	if (!restr_atoms || restr_atoms->GetNAtoms() == 0)
	{
		PrintLog("No Restrained Atoms defined");
		return FALSE;
	}

	int num_restr_atoms = restr_atoms->GetNAtoms();
	if (restr_ref_coords.size() != num_restr_atoms)
	{
		PrintLog("Error in MAmberMMModel::SaveAtomRestrArbalestIndForm() \n");
		PrintLog("Dimension of Reference Coordinates array is not equal \n");
		PrintLog("To the size of restrained atoms group \n");
		return FALSE;
	}

	double fc_restr = GetAtomRestrForceConst();
	sprintf(buf, "%12.6f", fc_restr);
	std::string fc_str = buf;

	std::ofstream os_desc;
	std::ofstream os_list;

	os_desc.open(restr_desc_fname.c_str());
	os_list.open(restr_list_fname.c_str());

	os_desc << "  <RestrainingRules> " << std::endl;
	os_list << "  <Restraints> " << std::endl;

	AtomIteratorAtomGroup aitr(restr_atoms);
	HaAtom* aptr;
	bool use_chain = false;
	aptr = aitr.GetFirstAtom();
	if (aptr)
	{
		HaChain* pchain_fst = aptr->GetHostChain();
		for (; aptr; aptr = aitr.GetNextAtom())
		{
			HaChain*   pchain = aptr->GetHostChain();
			if (pchain != pchain_fst)
			{
				use_chain = true;
				break;
			}
		}
	}

	MolSet* pmset = this->GetMolSet();
	ResidueIteratorMolSet ritr(pmset);
	std::map<HaAtom*, int> atom_res_seq_map; // map atoms to residue sequential number 
	int n_res_seq = 0;
	for (HaResidue* rptr : ritr)
	{
		n_res_seq++;
		for (HaAtom* aptr : *rptr)
			atom_res_seq_map[aptr] = n_res_seq;
	}
	

	int idx_r = -1;
	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
	{
		idx_r++;
		Vec3D crd_r = restr_ref_coords[idx_r];
		HaResidue* pres = aptr->GetHostRes();
		HaChain*   pchain = aptr->GetHostChain();
		HaMolecule* pmol = aptr->GetHostMol();

		std::string at_name = aptr->GetName();
		std::string res_name = pres->GetName(); 
		std::string name_mod = pres->GetNameModifier();
		std::string res_name_save = res_name;

		if( res_name == "HIS")
		{
			res_name_save = "HID";
			if( name_mod == "EPSILON") res_name_save = "HIE";
			if( name_mod == "PROT")    res_name_save = "HIP";
		}

		if (res_name == "CYS")
		{
			if (name_mod == "UNPROT")    res_name_save = "CYX";
		}

		//std::string res_n_str = harlem::ToString( pres->GetSerNo() );
		std::string res_n_str = std::to_string(atom_res_seq_map[aptr]); // use sequence number of the residue - it looks like Arbalest uses this

		std::string atom_id_arb;
		std::string atom_restr_id; 

		if( res_name.empty() )
		{
			atom_id_arb = at_name;
			atom_restr_id = "PositionRestraint_" + at_name; 
		}
		else
		{
			atom_id_arb = res_name_save + "(" + res_n_str + ")";
			if (use_chain)
			{
				atom_id_arb.push_back('\'');
				atom_id_arb.push_back(pchain->ident);
				atom_id_arb.push_back('\'');
			}
			atom_id_arb +=  ":" + at_name;
			atom_restr_id = "PositionRestraint_" + res_name_save + "_" + res_n_str;
			if (use_chain) atom_restr_id += std::string(1,pchain->ident);
			atom_restr_id += "_" + at_name;
			PrintLog("use_chain = %d atom_restr_id = %s \n", use_chain, atom_restr_id.c_str());
		}
		if( !os_desc.fail() )
		{
			os_desc << "      <RestraintDefinition IsRelativeCoords=\"false\" Title=\"" << atom_restr_id  << "\">" << std::endl;
			os_desc << "           <LinkToAtom_a>" << atom_id_arb << "</LinkToAtom_a>" << std::endl;
			sprintf(buf,"%8.3f",crd_r.GetX_Ang()); 
			os_desc << "            <Anchor_C X=\"" << boost::trim_copy( std::string(buf) ) << "\"";
			sprintf(buf,"%8.3f",crd_r.GetY_Ang());
			os_desc << " Y=\"" <<  boost::trim_copy( std::string(buf)) << "\"";
			sprintf(buf,"%8.3f",crd_r.GetZ_Ang());
			os_desc << " Z=\"" <<  boost::trim_copy( std::string(buf)) << "\"" << " />" << std::endl;

			os_desc <<  "           <ForceConstant_Ca>" <<  boost::trim_copy(fc_str) << "</ForceConstant_Ca> " << std::endl;
			os_desc <<  "           <EquilibriumLength_Ca>0</EquilibriumLength_Ca> " << std::endl;
			os_desc << "      </RestraintDefinition>" << std::endl;
		}
		if( !os_list.fail() )
		{
			os_list << "       <Restraint MolGroup=\"\" Title=\"" << atom_restr_id  << "\"  />" << std::endl;
		}
	}
	os_desc << "  </RestrainingRules> " << std::endl;
	os_list << "  </Restraints> " << std::endl;
	return TRUE;
}


int MolMechModel::SetMovingAtoms(std::string moving_atoms_name_new)
{
//	PrintLog(" MolMechModel::SetMovingAtoms() grp = %s\n",moving_atoms_name_new.c_str());
	AtomGroup* mv_atlist;
	std::string moving_atoms_old = moving_atoms;
	moving_atoms = "ALL_ATOMS";
	if( moving_atoms_name_new == "ALL_ATOMS" || moving_atoms_name_new.size() == 0) return TRUE;
	
	mv_atlist = pmset->GetAtomGroupByID(moving_atoms_name_new.c_str());
	if(mv_atlist != NULL)
	{
		moving_atoms = moving_atoms_name_new;
//		p_mm_mod->period_bcond = p_mm_mod->period_bcond.NO_PERIODICITY; // ??? is BELLY incompatible with periodical boundary conditions?
	}
	else
	{
		PrintLog("Error in MolMechModel::SetMovingAtoms() \n No atom array %s \n",moving_atoms_name_new.c_str()); 
	}

	if( moving_atoms != moving_atoms_old ) p_amber_model->SetUpdateDataFlag(TRUE);

	return TRUE;
}

bool MolMechModel::SetRestrainedAtoms(const std::string restr_atom_group_name)
{
	AtomGroup* p_atgrp_restr = NULL;
	std::string restr_atgrp_name = restr_atom_group_name;
	boost::trim(restr_atgrp_name);

	if( restr_atgrp_name == this->restrained_atoms ) return true;

	if( restr_atgrp_name == "NO_RESTRAINTS" ||  restr_atgrp_name.empty() )
	{
		restrained_atoms = "NO_RESTRAINTS";
	}
	else
	{
		p_atgrp_restr = pmset->GetAtomGroupByID(restr_atgrp_name.c_str());
	}

	if(p_atgrp_restr != NULL)
	{
		restrained_atoms = restr_atgrp_name;
		restr_ref_crd_type = RESTR_REFC_CURRENT_CRD;
		
		int na = p_atgrp_restr->GetNAtoms();
		restr_ref_coords.resize(na);

		int i;
		for(i = 0; i < na; i++)
		{
			restr_ref_coords[i].SetX( (*p_atgrp_restr)[i]->GetX() );
			restr_ref_coords[i].SetY( (*p_atgrp_restr)[i]->GetY() );
			restr_ref_coords[i].SetZ( (*p_atgrp_restr)[i]->GetZ() );
		}
	}
	else
	{
		PrintLog("Error in MolMechModel::SetRestrainedAtoms() \n Not found Atom Group: %s \n",restr_atgrp_name.c_str()); 
		restrained_atoms = "NO_RESTRAINTS";
	}
	
	if(p_amber_model->natom > 0 ) p_amber_model->SetAtomPosRestrData();
	return true;
}

int MolMechModel::SetRestrRefCrdFromXYZFile( const char* ref_crd_file_name_new )
{
	AtomGroup* p_restr_atoms = this->GetRestrAtoms();
	if( p_restr_atoms == NULL )
	{
		PrintLog(" Error in MolMechModel::SetRestrRefCrdFromXYZFile() \n");
		PrintLog(" Restrained Atoms Group is not set \n");
		return FALSE;
	}

	ifstream refc_fs( ref_crd_file_name_new );
	if( !refc_fs.good() )
	{	
		PrintLog(" Error in MolMechModel::SetRestrRefCrdFromXYZFile() \n");
		PrintLog(" Error to open file %s \n", ref_crd_file_name_new);
		return FALSE;
	}

	char buf[256];
	int na;
	refc_fs.getline(buf,256);
	istrstream ss(buf);
	ss >> na;
	if(!ss)
	{
		PrintLog(" Error in MolMechModel::SetRestrRefCrdFromXYZFile() \n");
		PrintLog(" Error reating atom number from file %s \n", ref_crd_file_name_new);
		return FALSE;
	}

	if( na != p_restr_atoms->GetNAtoms() )
	{
		PrintLog(" Error in MolMechModel::SetRestrRefCrdFromXYZFile() \n");
		PrintLog(" Number of atoms = %d in file %s \n", na, ref_crd_file_name_new);
		PrintLog(" Is different from the number of Restrained Atoms Group %d \n",p_restr_atoms->GetNAtoms()); 
		return FALSE;
	}

	restr_ref_coords.resize(na);
	int i;
	for(i = 0; i < na; i++)
	{
		refc_fs.getline(buf,256);
		istrstream line_s(buf);
		std::string id;
		int idx;
		double x,y,z;
		line_s >> idx;
		line_s >> id;
		line_s >> x;
		line_s >> y;
		line_s >> z;
		if( !line_s ) 
		{
			PrintLog(" Error in MolMechModel::SetRestrRefCrdFromXYZFile() \n");
			PrintLog(" Processing line %s \n", buf);
			restr_ref_coords.clear();
			return FALSE;
		}
		restr_ref_coords[i].SetX_Ang(x);
		restr_ref_coords[i].SetY_Ang(y);
		restr_ref_coords[i].SetZ_Ang(z);
	}
	restr_ref_crd_type = RESTR_REFC_XYZ_CRD_FILE;
	if(p_amber_model->natom > 0) p_amber_model->SetAtomPosRestrData();
	return TRUE;
}

int MolMechModel::SetRestrRefCrdFromStr( const std::string& crd_str )
{
	try
	{
		AtomGroup* p_restr_atoms = this->GetRestrAtoms();
		if( p_restr_atoms == NULL ) throw std::runtime_error(" Restrained Atoms Group is not set ");

		std::string crd_str_local = boost::trim_copy(crd_str);
		std::vector<std::string> str_vec;
		boost::split(str_vec,crd_str_local,boost::is_any_of(" "),boost::token_compress_on);

		int na = p_restr_atoms->GetNAtoms();
		int ncrd = str_vec.size();
		if( ncrd != 3*na ) throw std::runtime_error(" number of restr coords " + harlem::ToString(ncrd) + " nat = " + harlem::ToString(na));

		restr_ref_coords.resize(na);
		for( int i = 0; i < na; i++)
		{
			restr_ref_coords[i].SetX_Ang( atof(str_vec[3*i].c_str()) );
			restr_ref_coords[i].SetY_Ang( atof(str_vec[3*i+1].c_str()) );
			restr_ref_coords[i].SetZ_Ang( atof(str_vec[3*i+2].c_str()) );
		}
		restr_ref_crd_type = RESTR_REFC_XYZ_CRD_FILE;
		if(p_amber_model->natom > 0) p_amber_model->SetAtomPosRestrData();
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in MolMechModel::SetRestrRefCrdFromStr() \n");
		PrintLog(" %s\n",ex.what());
		return FALSE;
	}
	return TRUE;
	


}

int MolMechModel::SetAtomRestrForceConst(double restr_const_new)
{
	atom_restr_const = restr_const_new;
	if( p_amber_model->atm_weight.size() > 0 ) p_amber_model->atm_weight = restr_const_new;
	return TRUE;
}

double MolMechModel::GetAtomRestrForceConst()
{
	return atom_restr_const;
}

int MolMechModel::GetNumHarmConstr() const
{
	int i;
	int n = 0;
	for(i = 0; i < DistConstraints.size(); i++)
	{
		if( DistConstraints[i].cnt_type != AtomContactType::HARMONIC_CNT ) continue;
		n++;
	}
	return n;
}

int MolMechModel::SetDistConstrFromFile( const char* constr_file_name )
{
	std::ifstream ifs(constr_file_name);
	if( ifs.fail() )
	{
		PrintLog(" Error in MolMechModel::SetAtomDistConstrFromFile() \n");
		PrintLog(" Failed to open file %s \n", constr_file_name);
		return FALSE;
	}
	char buf[256];
	DistConstraints.clear();

	for(;;)
	{
		ifs.getline(buf,256);
		if( ifs.fail() ) return TRUE;
		if( buf[0] == '#' ||  buf[0] == '!' ) continue;

		std::string constr_type;
		std::string atid_1;
		std::string atid_2;

		std::istrstream iss(buf);

		iss >> constr_type;
		iss >> atid_1;
		iss >> atid_2;

		if( iss.fail() ) 
		{
			PrintLog(" Invalid line in Constraints File: %s \n",buf);
			continue;
		}
				
		int iat1, iat2;
		HaAtom* aptr1;
		HaAtom* aptr2;

		if( boost::algorithm::all( atid_1, boost::algorithm::is_digit() ) )
		{
			iat1 = atoi( atid_1.c_str() );
			if( iat1 < 0 || iat1 >= Atoms.size() ) 
			{
				PrintLog( " Invalid atom number in line %s \n", buf);
				continue;
			}
			aptr1 = Atoms[iat1];
		}
		else
		{
			aptr1 = pmset->GetAtomByRef(atid_1.c_str());
			if( aptr1 == NULL ) 
			{
				PrintLog( " Invalid atom reference in line %s \n", buf);
				continue;
			}
		}
		
		if( boost::algorithm::all( atid_2, boost::algorithm::is_digit() ) )
		{
			iat2 = atoi( atid_2.c_str() );
			if( iat2 < 0 || iat2 >= Atoms.size() ) 
			{
				PrintLog( " Invalid atom number in line %s \n", buf);
				continue;
			}
			aptr2 = Atoms[iat2];
		}
		else
		{
			aptr2 = pmset->GetAtomByRef(atid_2.c_str());
			if( aptr2 == NULL ) 
			{
				PrintLog( " Invalid atom reference in line %s \n", buf);
				continue;
			}
		}

		if( aptr1 == aptr2 ) 
		{
			PrintLog(" The same atoms in line %s \n", buf );
			continue;
		}

		AtomContact cnt(aptr1,aptr2);
		int ires = cnt.cnt_type.SetWithLabel(constr_type.c_str());
		if( !ires )
		{
			PrintLog(" Invalid constraint type  in line: %s \n",buf  );
			continue;
		}

		PrintLog(" MolMechModel::SetDistConstrFromFile pt 6 \n");
		
		double dist, fc;
		
		iss >> dist;

		if( iss.fail() )
		{
			dist = Vec3D::CalcDistance(aptr1,aptr2);
		}
		
		iss >> fc;

		if( cnt.cnt_type == AtomContactType::HARMONIC_CNT )
		{
			if( iss.fail() )
			{
				fc = 10.0;
			}
			AddHarmConstraint(aptr1,aptr2,dist,fc);
		}
	}
	return TRUE;
}

std::string MolMechModel::GetFFSymbolFromMolStruct(HaAtom* aptr)
{
	std::string atn = aptr->GetName();
	HaResidue* pres = aptr->GetHostRes();

	if (pres->IsAlchemicalTransformationSet())
	{
		AtomFFParam* at_ff_par_mut = pres->p_res_transform->GetAtomFFParamMut(aptr);
		if (at_ff_par_mut) return at_ff_par_mut->ff_symbol;
	}

	std::string ff_s = "DU";
	AtomGroup bonded_atoms = aptr->GetBondedAtoms();
	if (bonded_atoms.size() == 1) // so far only set rules for atoms with  one neighbor
	{
		HaAtom* aptr_host = bonded_atoms[0]; // Heavy Atom  a given Terminal atom is Attached
		std::string host_ff_s = aptr_host->GetFFSymbol();
		AtomGroup bonded_atoms_host = aptr_host->GetBondedAtoms();  // Atoms bonded to the Heavy Atom host
		int nh_host = 0;
		for (HaAtom* aptr_b : bonded_atoms_host)
			if (aptr_b->IsHydrogen()) nh_host++;

		if (host_ff_s == "CT") ff_s = "H1";
		for (HaAtom* aptr_b : bonded_atoms_host)
		{
			if (aptr_b->GetFFSymbol() == "HC") ff_s = "HC";
		}
	}
	return ff_s;  
}

int MolMechModel::SetStdValParams()
{
	MMForceField* p_ff = MMForceField::GetMMForceField(ff_type);
	if( p_ff == NULL)
	{
		PrintLog("Error in MolMechModel::SetStdValParams() \n");
		PrintLog("Can not initiate Force Field %s \n", ff_type.label() );
		return FALSE;
	}

	set<MMBond>::iterator mbitr = MBonds.begin();
	for(; mbitr != MBonds.end(); mbitr++)
	{
		MMBond& bnd = (MMBond&)(*mbitr);
		HaAtom* pt1 = bnd.pt1;
		HaAtom* pt2 = bnd.pt2;

		if( bnd.set_type < SET_FF_FIELD)
		{
			std::string at1_ff_s = pt1->GetFFSymbol();
			std::string at2_ff_s = pt2->GetFFSymbol();

			if (at1_ff_s == "DU") at1_ff_s = GetFFSymbolFromMolStruct(pt1);
			if (at2_ff_s == "DU") at2_ff_s = GetFFSymbolFromMolStruct(pt2);

			vector<double> bpar = p_ff->FindBondParamFromSymbol(at1_ff_s, at2_ff_s);
			if( bpar.size() > 1 )
			{
				bnd.r0 = bpar[0];
				bnd.fc = bpar[1];
				bnd.set_type = SET_FF_FIELD;
				continue;
			}
			
			PrintLog("Can't find Force field parameters for a bond between atoms \n %s %s \n with FF symbols %s %s \n",
				pt1->GetRef(),pt2->GetRef(), pt1->GetFFSymbol(), pt2->GetFFSymbol() );
			PrintLog("Set Default Parameters For the bond based on Current Geometry \n");
		
			double mbdist = Vec3D::CalcDistance(pt1,pt2,ANGSTROM_U);

			bnd.r0 = mbdist;
			bnd.fc = 100.0;
			bnd.set_type = SET_DEFAULT;
		}
	}

	mbitr = MBonds_mut.begin();
	for (; mbitr != MBonds_mut.end(); mbitr++)
	{
		MMBond& bnd = (MMBond&)(*mbitr);
		HaAtom* pt1 = bnd.pt1;
		HaAtom* pt2 = bnd.pt2;

		if (bnd.set_type < SET_FF_FIELD)
		{
			std::string at1_ff_s = pt1->GetFFSymbol();
			std::string at2_ff_s = pt2->GetFFSymbol();

			if(GetAtomFFParamMut(pt1) != nullptr) at1_ff_s = GetAtomFFParamMut(pt1)->ff_symbol;
			if(GetAtomFFParamMut(pt2) != nullptr) at2_ff_s = GetAtomFFParamMut(pt2)->ff_symbol;

			if (at1_ff_s == "DU") at2_ff_s = GetFFSymbolFromMolStruct(pt1);
			if (at2_ff_s == "DU") at2_ff_s = GetFFSymbolFromMolStruct(pt2);

			vector<double> bpar = p_ff->FindBondParamFromSymbol(at1_ff_s, at2_ff_s);
			if (bpar.size() > 1)
			{
				bnd.r0 = bpar[0];
				bnd.fc = bpar[1];
				bnd.set_type = SET_FF_FIELD;
				continue;
			}

			PrintLog("Can't find Force field parameters for a bond between mutated atoms \n %s %s \n with FF symbols %s %s \n",
				pt1->GetRef(), pt2->GetRef(), at1_ff_s, at2_ff_s);
			PrintLog("Set Standard Parameters For the bond \n");

			double mbdist = Vec3D::CalcDistance(pt1, pt2, ANGSTROM_U);

			bnd.r0 = mbdist;
			bnd.fc = 100.0;
			bnd.set_type = SET_DEFAULT;
		}
	}

	for(auto vitr = ValAngles.begin(); vitr != ValAngles.end(); vitr++)
	{
		MMValAngle& va = (MMValAngle&)(*vitr);
		HaAtom* pt1 = va.pt1;
		HaAtom* pt2 = va.pt2;
		HaAtom* pt3 = va.pt3;
		
		if( va.set_type < SET_FF_FIELD)
		{
			std::string at1_ff_s = pt1->GetFFSymbol();
			std::string at2_ff_s = pt2->GetFFSymbol();
			std::string at3_ff_s = pt3->GetFFSymbol();

			if (at1_ff_s == "DU") at1_ff_s = GetFFSymbolFromMolStruct(pt1);
			if (at2_ff_s == "DU") at2_ff_s = GetFFSymbolFromMolStruct(pt2);
			if (at3_ff_s == "DU") at3_ff_s = GetFFSymbolFromMolStruct(pt3);

			vector<double> vpar(2);
			vpar = p_ff->FindValAngleParamFromSymbol(at1_ff_s, at2_ff_s, at3_ff_s);
			if( vpar.size() > 1)
			{
				va.a0 = vpar[0];
				va.fc = vpar[1];
				va.set_type = SET_FF_FIELD;
				continue;
			}
			
			PrintLog("Can't find Force field parameters for a val angle between atoms \n %s %s %s \n with FF symbols %s %s %s \n",
					pt1->GetRef(), pt2->GetRef(), pt3->GetRef(),
				    at1_ff_s, at2_ff_s, at3_ff_s);
			PrintLog("Set Standard Parameters For the Valence Angle \n");
			
			double vangle = Vec3D::CalcAngle(pt1,pt2,pt3)/DEG_TO_RAD;
			
			va.a0 = vangle;
			va.fc = 63.0;
			va.set_type = SET_DEFAULT;		
		}
	}

	for (auto vitr = ValAngles_mut.begin(); vitr != ValAngles_mut.end(); vitr++)
	{
		MMValAngle& va = (MMValAngle&)(*vitr);
		HaAtom* pt1 = va.pt1;
		HaAtom* pt2 = va.pt2;
		HaAtom* pt3 = va.pt3;

		if (va.set_type < SET_FF_FIELD)
		{
			std::string at1_ff_s = pt1->GetFFSymbol();
			std::string at2_ff_s = pt2->GetFFSymbol();
			std::string at3_ff_s = pt3->GetFFSymbol();

			if (GetAtomFFParamMut(pt1) != nullptr) at1_ff_s = GetAtomFFParamMut(pt1)->ff_symbol;
			if (GetAtomFFParamMut(pt2) != nullptr) at2_ff_s = GetAtomFFParamMut(pt2)->ff_symbol;
			if (GetAtomFFParamMut(pt3) != nullptr) at3_ff_s = GetAtomFFParamMut(pt3)->ff_symbol;

			if (at1_ff_s == "DU") at1_ff_s = GetFFSymbolFromMolStruct(pt1);
			if (at2_ff_s == "DU") at2_ff_s = GetFFSymbolFromMolStruct(pt2);
			if (at3_ff_s == "DU") at3_ff_s = GetFFSymbolFromMolStruct(pt3);

			HaVec_double vpar(2);
			vpar = p_ff->FindValAngleParamFromSymbol(at1_ff_s, at2_ff_s, at3_ff_s);
			if (vpar.size() > 1)
			{
				va.a0 = vpar[0];
				va.fc = vpar[1];
				va.set_type = SET_FF_FIELD;
				continue;
			}

			PrintLog("Can't find Force field parameters (mutated state ) for a val angle between atoms \n %s %s %s \n with FF symbols %s %s %s \n",
				pt1->GetRef(), pt2->GetRef(), pt3->GetRef(),
				at1_ff_s , at2_ff_s, at3_ff_s );
			PrintLog("Set Standard Parameters For the Valence Angle \n");

			double vangle = Vec3D::CalcAngle(pt1, pt2, pt3) / DEG_TO_RAD;

			va.a0 = vangle;
			va.fc = 63.0;
			va.set_type = SET_DEFAULT;
		}
	}
	
	for( shared_ptr<MMDihedral> lditr : Dihedrals )
	{
		MMDihedral& dih = *lditr;

		HaAtom* pt1 = dih.pt1;
		HaAtom* pt2 = dih.pt2;
		HaAtom* pt3 = dih.pt3;
		HaAtom* pt4 = dih.pt4;
			
		if( dih.set_type <= SET_FF_FIELD)
		{
			dih.ClearParams();

			std::string at1_ff_s = pt1->GetFFSymbol();
			std::string at2_ff_s = pt2->GetFFSymbol();
			std::string at3_ff_s = pt3->GetFFSymbol();
			std::string at4_ff_s = pt4->GetFFSymbol();

			if (at1_ff_s == "DU") at1_ff_s = GetFFSymbolFromMolStruct(pt1);
			if (at2_ff_s == "DU") at2_ff_s = GetFFSymbolFromMolStruct(pt2);
			if (at3_ff_s == "DU") at3_ff_s = GetFFSymbolFromMolStruct(pt3);
			if (at4_ff_s == "DU") at4_ff_s = GetFFSymbolFromMolStruct(pt4);

			vector<double> dpar = p_ff->FindDihedralParamFromSymbol(at1_ff_s, at2_ff_s, at3_ff_s, at4_ff_s);
			if(dpar.size() > 1 )
			{
				int nt = dpar.size()/4;
				int it;
				for( it = 0; it < nt; it++)
				{
					dih.AddTerm(dpar[it*4],dpar[it*4+1],dpar[it*4+2],dpar[it*4+3]);
 				}
				dih.set_type = SET_FF_FIELD;
				continue;
			}

//	IGOR TMP		PrintLog("Can't find Force field parameters for a Dihedral angle between atoms \n  %s %s %s %s\n with FF symbols %s %s %s %s\n",
//					  buf1, buf2, buf3, buf4,
//					  pt1->FFSymbol,pt2->FFSymbol,pt3->FFSymbol,pt4->FFSymbol );
//			PrintLog("Set Standard Parameters For the Dihedral Angle \n");

			double dih_angle = Vec3D::CalcDihedral(pt1,pt2,pt3,pt4);	
			dih.AddTerm(4.0,dih_angle,0.0);
			dih.set_type = SET_DEFAULT;
		}
	}

	for (shared_ptr<MMDihedral> lditr : Dihedrals_mut)
	{
		MMDihedral& dih = *lditr;

		HaAtom* pt1 = dih.pt1;
		HaAtom* pt2 = dih.pt2;
		HaAtom* pt3 = dih.pt3;
		HaAtom* pt4 = dih.pt4;

		if (dih.set_type < SET_FF_FIELD)
		{
			dih.ClearParams();

			std::string at1_ff_s = pt1->GetFFSymbol();
			std::string at2_ff_s = pt2->GetFFSymbol();
			std::string at3_ff_s = pt3->GetFFSymbol();
			std::string at4_ff_s = pt4->GetFFSymbol();

			if (GetAtomFFParamMut(pt1) != nullptr) at1_ff_s = GetAtomFFParamMut(pt1)->ff_symbol;
			if (GetAtomFFParamMut(pt2) != nullptr) at2_ff_s = GetAtomFFParamMut(pt2)->ff_symbol;
			if (GetAtomFFParamMut(pt3) != nullptr) at3_ff_s = GetAtomFFParamMut(pt3)->ff_symbol;
			if (GetAtomFFParamMut(pt4) != nullptr) at4_ff_s = GetAtomFFParamMut(pt4)->ff_symbol;

			if (at1_ff_s == "DU") at1_ff_s = GetFFSymbolFromMolStruct(pt1);
			if (at2_ff_s == "DU") at2_ff_s = GetFFSymbolFromMolStruct(pt2);
			if (at3_ff_s == "DU") at3_ff_s = GetFFSymbolFromMolStruct(pt3);
			if (at4_ff_s == "DU") at4_ff_s = GetFFSymbolFromMolStruct(pt4);

			vector<double> dpar = p_ff->FindDihedralParamFromSymbol(at1_ff_s, at2_ff_s, at3_ff_s, at4_ff_s);
			if (dpar.size() > 1)
			{
				int nt = dpar.size() / 4;
				int it;
				for (it = 0; it < nt; it++)
				{
					dih.AddTerm(dpar[it * 4], dpar[it * 4 + 1], dpar[it * 4 + 2], dpar[it * 4 + 3]);
				}
				dih.set_type = SET_FF_FIELD;
				continue;
			}

			//	IGOR TMP		PrintLog("Can't find Force field parameters for a Dihedral angle between atoms \n  %s %s %s %s\n with FF symbols %s %s %s %s\n",
			//					  buf1, buf2, buf3, buf4,
			//					  pt1->FFSymbol.c_str(),pt2->FFSymbol.c_str(),pt3->FFSymbol.c_str(),pt4->FFSymbol.c_str() );
			//			PrintLog("Set Standard Parameters For the Dihedral Angle \n");

			double dih_angle = Vec3D::CalcDihedral(pt1, pt2, pt3, pt4);
			dih.AddTerm(4.0, dih_angle, 0.0);
			dih.set_type = SET_DEFAULT;
		}
	}

	for( shared_ptr<MMDihedral> lditr : ImprDihedrals)
	{
		MMDihedral& dih = *lditr;
		
		HaAtom* pt1 = (*lditr).pt1;
		HaAtom* pt2 = (*lditr).pt2;
		HaAtom* pt3 = (*lditr).pt3;
		HaAtom* pt4 = (*lditr).pt4;	
		
		if( dih.set_type < SET_FF_FIELD)
		{
			dih.ClearParams();

			std::string at1_ff_s = pt1->GetFFSymbol();
			std::string at2_ff_s = pt2->GetFFSymbol();
			std::string at3_ff_s = pt3->GetFFSymbol();
			std::string at4_ff_s = pt4->GetFFSymbol();

			if (at1_ff_s == "DU") at1_ff_s = GetFFSymbolFromMolStruct(pt1);
			if (at2_ff_s == "DU") at2_ff_s = GetFFSymbolFromMolStruct(pt2);
			if (at3_ff_s == "DU") at3_ff_s = GetFFSymbolFromMolStruct(pt3);
			if (at4_ff_s == "DU") at4_ff_s = GetFFSymbolFromMolStruct(pt4);

			vector<double> dpar = p_ff->FindDihedralParamFromSymbol(at1_ff_s, at2_ff_s, at3_ff_s, at4_ff_s, true);
			if(dpar.size() > 1  )
			{
				dih.AddTerm(dpar[0],dpar[1],dpar[2],dpar[3]);
				dih.set_type = SET_FF_FIELD;
				continue;
			}
				
			PrintLog("Can't find Force field parameters for a Improper Dihedral angle between atoms:\n  %s %s %s %s\n with FF symbols %s %s %s %s\n",
				pt1->GetRef(), pt2->GetRef(), pt3->GetRef(), pt4->GetRef(),
				pt1->GetFFSymbol(),pt2->GetFFSymbol(),pt3->GetFFSymbol(),pt4->GetFFSymbol() );
		
			PrintLog("Set Standard Parameters For the Improper Dihedral Angle \n");
		
			double dih_angle = Vec3D::CalcDihedral(pt1,pt2,pt3,pt4);
			dih.AddTerm(2.0, 180.0, 20.0); 	
			dih.set_type = SET_DEFAULT;
		}
	}

	for (shared_ptr<MMDihedral> lditr : ImprDihedrals_mut)
	{
		MMDihedral& dih = *lditr;

		HaAtom* pt1 = (*lditr).pt1;
		HaAtom* pt2 = (*lditr).pt2;
		HaAtom* pt3 = (*lditr).pt3;
		HaAtom* pt4 = (*lditr).pt4;

		if (dih.set_type < SET_FF_FIELD)
		{
			dih.ClearParams();
			std::string at1_ff_s = pt1->GetFFSymbol();
			std::string at2_ff_s = pt2->GetFFSymbol();
			std::string at3_ff_s = pt3->GetFFSymbol();
			std::string at4_ff_s = pt4->GetFFSymbol();

			if (GetAtomFFParamMut(pt1) != nullptr) at1_ff_s = GetAtomFFParamMut(pt1)->ff_symbol;
			if (GetAtomFFParamMut(pt2) != nullptr) at2_ff_s = GetAtomFFParamMut(pt2)->ff_symbol;
			if (GetAtomFFParamMut(pt3) != nullptr) at3_ff_s = GetAtomFFParamMut(pt3)->ff_symbol;
			if (GetAtomFFParamMut(pt4) != nullptr) at4_ff_s = GetAtomFFParamMut(pt4)->ff_symbol;

			if (at1_ff_s == "DU") at1_ff_s = GetFFSymbolFromMolStruct(pt1);
			if (at2_ff_s == "DU") at2_ff_s = GetFFSymbolFromMolStruct(pt2);
			if (at3_ff_s == "DU") at3_ff_s = GetFFSymbolFromMolStruct(pt3);
			if (at4_ff_s == "DU") at4_ff_s = GetFFSymbolFromMolStruct(pt4);

			vector<double> dpar = p_ff->FindDihedralParamFromSymbol(at1_ff_s, at2_ff_s, at3_ff_s, at4_ff_s, true);
			if (dpar.size() > 1)
			{
				dih.AddTerm(dpar[0], dpar[1], dpar[2], dpar[3]);
				dih.set_type = SET_FF_FIELD;
				continue;
			}

			PrintLog("Can't find Force field parameters for a Improper Dihedral angle between atoms:\n  %s %s %s %s\n with FF symbols %s %s %s %s\n",
				pt1->GetRef(), pt2->GetRef(), pt3->GetRef(), pt4->GetRef(),
				pt1->GetFFSymbol(), pt2->GetFFSymbol(), pt3->GetFFSymbol(), pt4->GetFFSymbol());

			PrintLog("Set Standard Parameters For the Improper Dihedral Angle \n");

			double dih_angle = Vec3D::CalcDihedral(pt1, pt2, pt3, pt4);
			dih.AddTerm(2.0, 180.0, 20.0);
			dih.set_type = SET_DEFAULT;
		}
	}

	return TRUE;
}

int MolMechModel::SetStdVdWParams()
{
	MMForceField* p_ff = MMForceField::GetMMForceField(ff_type);
	if( p_ff == NULL)
	{
		PrintLog("Error in MolMechModel::SetStdVdWParams() \n");
		PrintLog("Can not initiate Force Field %s \n",ff_type.label() );
		return FALSE;
	}

	vector<HaAtom*>::iterator pitr;

	for(pitr = Atoms.begin(); pitr != Atoms.end(); pitr++)
	{
		HaAtom* aptr = *pitr;
		
		vector<double> ppar = p_ff->FindPointParamFromSymbol(aptr->GetFFSymbol());

		if(ppar.size() > 1)
		{
			aptr->SetVdWRad(ppar[0]);
			aptr->SetVdWEne(ppar[1]);
			continue;
		}

		if (aptr->GetFFSymbol() == "DU")  // Set Zero VDW parameters for DUMMY atoms 
		{
			aptr->SetVdWRad(1.0);
			aptr->SetVdWEne(0.0);
			continue;
		}

		PrintLog("Can't Find VdW parameters for Atom %s FF symbol = \'%s\' \n",aptr->GetRef(),aptr->GetFFSymbol());
        PrintLog("Will Set Standard VdW parameters \n");

		std::string elem_symbol = aptr->GetStdSymbol(); 
//		aptr->SetFFSymbol(elem_symbol);
		ppar = p_ff->FindPointParamFromSymbol(elem_symbol);
		if(ppar.size() > 1)
		{
			aptr->SetVdWRad(ppar[0]);
			aptr->SetVdWEne(ppar[1]);
		}
		else
		{
//			aptr->SetFFSymbol("XXST");
			aptr->SetVdWRad(1.908);
			aptr->SetVdWEne(0.086);
		}
	}
	return TRUE;
}	

int MolMechModel::AddAtomsToExcludedAtomList(HaAtom* aptr1, HaAtom* aptr2, PtrIntMap& atoms_idx, bool mutated_state )
{
	if( aptr1== NULL || aptr2 == NULL) return FALSE;

	if( atoms_idx.count(aptr1) == 0 ) return FALSE;
	if( atoms_idx.count(aptr2) == 0 ) return FALSE;

	int ipt1 = atoms_idx[aptr1];
	int ipt2 = atoms_idx[aptr2];
			
	if (mutated_state)
	{
		if (ipt2 > ipt1) excluded_atom_list_mut[ipt1].insert(aptr2);
		if (ipt1 > ipt2) excluded_atom_list_mut[ipt2].insert(aptr1);
	}
	else
	{
		if (ipt2 > ipt1) excluded_atom_list[ipt1].insert(aptr2);
		if (ipt1 > ipt2) excluded_atom_list[ipt2].insert(aptr1);
	}

	return TRUE;
}


bool MolMechModel::BuildExcludedAtomList()
{
	excluded_atom_list.clear();
	excluded_atom_list_mut.clear();
	
	int nn = Atoms.size();
	set<HaAtom*> tmp_set;

	excluded_atom_list.resize(nn,tmp_set);
	excluded_atom_list_mut.resize(nn, tmp_set);

	PtrIntMap pt_indx_map;
	
	int i;
	for( i = 0; i < nn; i++)
	{
		pt_indx_map[ Atoms[i] ] = i;
	}

// include into the excluded list atoms separated by one, two or
// three bonds
// 
	for(const MMBond& bond : MBonds )
	{
		AddAtomsToExcludedAtomList(bond.pt1, bond.pt2, pt_indx_map);
	}
	for (const MMBond& bond : MBonds_mut)
	{
		AddAtomsToExcludedAtomList(bond.pt1, bond.pt2, pt_indx_map, true);
	}

	for(const MMValAngle& vang : ValAngles)
	{
		AddAtomsToExcludedAtomList(vang.pt1, vang.pt2,pt_indx_map);
		AddAtomsToExcludedAtomList(vang.pt1, vang.pt3,pt_indx_map);
		AddAtomsToExcludedAtomList(vang.pt2, vang.pt3,pt_indx_map);
	}

	for (const MMValAngle& vang : ValAngles_mut)
	{
		AddAtomsToExcludedAtomList(vang.pt1, vang.pt2, pt_indx_map, true);
		AddAtomsToExcludedAtomList(vang.pt1, vang.pt3, pt_indx_map, true);
		AddAtomsToExcludedAtomList(vang.pt2, vang.pt3, pt_indx_map, true);
	}
		
	for( shared_ptr<MMDihedral> sp_dih: Dihedrals)
	{	
		AddAtomsToExcludedAtomList(sp_dih->pt1, sp_dih->pt2,pt_indx_map);
		AddAtomsToExcludedAtomList(sp_dih->pt1, sp_dih->pt3,pt_indx_map);
		AddAtomsToExcludedAtomList(sp_dih->pt1, sp_dih->pt4,pt_indx_map);
		AddAtomsToExcludedAtomList(sp_dih->pt2, sp_dih->pt3,pt_indx_map);
		AddAtomsToExcludedAtomList(sp_dih->pt2, sp_dih->pt4,pt_indx_map);
		AddAtomsToExcludedAtomList(sp_dih->pt3, sp_dih->pt4,pt_indx_map);
	}

	for (shared_ptr<MMDihedral> sp_dih : Dihedrals_mut)
	{
		AddAtomsToExcludedAtomList(sp_dih->pt1, sp_dih->pt2, pt_indx_map, true);
		AddAtomsToExcludedAtomList(sp_dih->pt1, sp_dih->pt3, pt_indx_map, true);
		AddAtomsToExcludedAtomList(sp_dih->pt1, sp_dih->pt4, pt_indx_map, true);
		AddAtomsToExcludedAtomList(sp_dih->pt2, sp_dih->pt3, pt_indx_map, true);
		AddAtomsToExcludedAtomList(sp_dih->pt2, sp_dih->pt4, pt_indx_map, true);
		AddAtomsToExcludedAtomList(sp_dih->pt3, sp_dih->pt4, pt_indx_map, true);
	}

	vector<AtomContact>::iterator vdw_itr;

	for( vdw_itr = DistConstraints.begin(); vdw_itr != DistConstraints.end(); vdw_itr++ )
	{
		AtomContact& vdw_cnt = *vdw_itr;
		if( vdw_cnt.cnt_type == AtomContactType::VDW_CNT_6_12 || vdw_cnt.cnt_type == AtomContactType::VDW_CNT_6_12_NO_REP 
			|| vdw_cnt.cnt_type == AtomContactType::VDW_CNT_10_12 )
		{
			AddAtomsToExcludedAtomList(vdw_cnt.pt1, vdw_cnt.pt2, pt_indx_map);
		}
	}

	return true;
}

bool MolMechModel::BuildNonBondContactList()
{
	nonbond_contact_list.clear();
	nonbond_contact_list_mut.clear();

	int i,j;
	int nn = Atoms.size();

	if(excluded_atom_list.size() != nn)
	{
		ErrorInMod("MolMechModel::BuildNonBondContactList()",
			"The size of excluded atom list is not equal to the number of Atoms");
		return false;
	}

	if (excluded_atom_list_mut.size() != nn)
	{
		ErrorInMod("MolMechModel::BuildNonBondContactList()",
			"The size of mutated state excluded atom list is not equal to the number of Atoms");
		return false;
	}

	double cut2= nb_cut_dist * nb_cut_dist;
	nonbond_contact_list.resize(nn);
	nonbond_contact_list_mut.resize(nn);

	for(i = 0 ; i < nn-1; i++)
	{
		HaAtom* pt = Atoms[i];
		double x1 = pt->GetX();
		double y1 = pt->GetY();
		double z1 = pt->GetZ();

        set<HaAtom*>& pt_nonb_contact_list = nonbond_contact_list[i];
		set<HaAtom*>& pt_nonb_contact_list_mut = nonbond_contact_list_mut[i];
		for(j = i+1; j < nn; j++)
		{
			HaAtom* pt2 = Atoms[j];
			double xx = pt2->GetX() - x1;
			double r2 = xx*xx;
			if( r2 > cut2)
				continue;
			double yy = pt2->GetY() - y1;
			r2 += yy*yy;
			if( r2 > cut2)
				continue;
			double zz = pt2->GetZ() - z1;
			r2 += zz*zz;
			if( r2 > cut2)
				continue;

			pt_nonb_contact_list.insert(pt2);
			pt_nonb_contact_list_mut.insert(pt2);
		}

		for(HaAtom* aptr: excluded_atom_list[i])
		{
			pt_nonb_contact_list.erase(aptr);
		}
		for (HaAtom* aptr : excluded_atom_list_mut[i])
		{
			pt_nonb_contact_list_mut.erase(aptr);
		}
	}
	return true;
}

bool MolMechModel::has_mut_atoms_in_group(const AtomGroup& group)
{
	bool has_mut_atoms = false;
	for (HaAtom* aptr : group)
	{
		if (GetAtomFFParamMut(aptr) != nullptr ) {
			has_mut_atoms = true;
			break;
		}
	}
	return has_mut_atoms; 
}

int MolMechModel::BuildGrpGrpExcludedList(AtomContainer* group1, AtomContainer* group2)
{
	int nn = Atoms.size();
	excluded_atom_list.clear();
	excluded_atom_list.resize(nn);

	set<HaAtom*,less<HaAtom*> >all_atom_set;
	AtomIntMap& pt_indx_map = GetAtIdxMap(TRUE); 
	
	int i;
	HaAtom* aptr;
	for(i = 0; i < nn; i++)
	{
		all_atom_set.insert( Atoms[i]);
	}
	
	AtomIteratorGen aitr1(group1);
	AtomIteratorGen aitr2(group2);
	set<HaAtom*,less<HaAtom*> >group1_atom_set;
	set<HaAtom*,less<HaAtom*> >group2_atom_set;

	for(aptr= aitr1.GetFirstAtom(); aptr != NULL; aptr= aitr1.GetNextAtom())
	{
		if(all_atom_set.count(aptr) == 0 )
		{
			PrintLog("Warning in MolMechModel::BuildGrpGrpExcludedList() \n");
			PrintLog("One of the atoms in the atom colection does not belong Atoms of HaMolMechMod \n");
			continue;
		}
		group1_atom_set.insert(aptr);
	}
	for(aptr= aitr2.GetFirstAtom(); aptr != NULL; aptr= aitr2.GetNextAtom())
	{
		if(all_atom_set.count(aptr) == 0 )
		{
			PrintLog("Warning in MolMechModel::BuildGrpGrpExcludedList() \n");
			PrintLog("One of the atoms in the atom colection does not belong Atoms of HaMolMechMod \n");
			continue;
		}
		group2_atom_set.insert(aptr);
	}
	
	set<HaAtom*,less<HaAtom*> >::iterator aitr_set_1;
	set<HaAtom*,less<HaAtom*> >::iterator aitr_set_2;

	for(aitr_set_1 = group1_atom_set.begin(); aitr_set_1 != group1_atom_set.end(); aitr_set_1++)
	{
		HaAtom* aptr1 = *aitr_set_1;
		int idx = pt_indx_map[aptr1];
		for(aitr_set_2 = group1_atom_set.begin(); aitr_set_2 != group1_atom_set.end(); aitr_set_2++)
		{
			HaAtom* aptr2 = *aitr_set_2;
			if( aptr1 == aptr2) continue;
			excluded_atom_list[idx].insert(aptr2);
		}
	}
	for(aitr_set_1 = group2_atom_set.begin(); aitr_set_1 != group2_atom_set.end(); aitr_set_1++)
	{
		HaAtom* aptr1 = *aitr_set_1;
		int idx = pt_indx_map[aptr1];
		for( aitr_set_2 = group2_atom_set.begin(); aitr_set_2 != group2_atom_set.end(); aitr_set_2++)
		{
			HaAtom* aptr2 = *aitr_set_2;
			if( aptr1 == aptr2) continue;
			excluded_atom_list[idx].insert(aptr2);
		}
	}
	return TRUE;
}

int MolMechModel::BuildGrpGrpNonBondList(AtomContainer* group1, AtomContainer* group2)
{
	nonbond_contact_list.clear();
	int nn = Atoms.size();
	nonbond_contact_list.resize(nn);
	
	set<HaAtom*, less<HaAtom*> > Points_set;

	int i;
	for(i = 0; i < nn; i++)
	{
		Points_set.insert(Atoms[i]);
	}

// Map HaAtom* to their indexes in Atoms array
	AtomIntMap& pt_indx_map  = GetAtIdxMap(TRUE); 
	
	for(i= 0; i < nn; i++)
	{
		nonbond_contact_list[i].clear();
	}	

	HaAtom* aptr1;
	HaAtom* aptr2;

	AtomIteratorGen aitr1(group1);
	AtomIteratorGen aitr2(group2);

	for(aptr1= aitr1.GetFirstAtom(); aptr1 != NULL; aptr1= aitr1.GetNextAtom())
	{
		if(aptr1->GetHostMolSet() != pmset )
		{
			PrintLog("Error in MolMechModel::BuildGrpGrpNonBondList() \n");
			PrintLog("One of the atoms in the atom coleection doesnt't belong to the Molecular Set Associated with the module \n");
			return FALSE;
		}
	}
	for(aptr2= aitr2.GetFirstAtom(); aptr2 != NULL; aptr2= aitr2.GetNextAtom())
	{
		if(aptr2->GetHostMolSet() != pmset )
		{
			PrintLog("Error in MolMechModel::BuildGrpGrpNonBondList() \n");
			PrintLog("One of the atoms in the atom coleection doesnt't belong to the Molecular Set Associated with the module \n");
			return FALSE;
		}
	}

	for(aptr1= aitr1.GetFirstAtom(); aptr1 != NULL; aptr1= aitr1.GetNextAtom())
	{
		set<HaAtom*, less<HaAtom*> >::iterator pitr;
			pitr = Points_set.find(aptr1);
		if(pitr == Points_set.end())
		{
			continue;	
		}
		HaAtom* pt1 = *pitr;
		
		for(aptr2 = aitr2.GetFirstAtom(); aptr2 != NULL; aptr2 = aitr2.GetNextAtom())
		{
			if(aptr1 == aptr2)
			{
				ErrorInMod("MolMechModel::BuildGrpGrpNonBondList()",
					" Groups have common atoms ");
				return false;
			}
			pitr = Points_set.find(aptr2);
			if(pitr == Points_set.end())
			{
				continue;	
			}
			HaAtom* pt2 = *pitr;
			
//			if( Vec3D::CalcDistance(pt1, pt2) < 10.0 ) // TEMPORAL! Hard coded nonb-list 
//			{
				if(pt1 < pt2) 
				{
					int ipt1 = pt_indx_map[pt1];
					nonbond_contact_list[ipt1].insert(pt2);
				}
				else if(pt1 > pt2)
				{
					int ipt2 = pt_indx_map[pt2];
					nonbond_contact_list[ipt2].insert(pt1);				
				}
//			}
		}
	}
	return TRUE;
}

bool MolMechModel::SetHBondRestraints(double force_const)
{	
	if(!pmset->HBonds_found) 
	{
		MolEditor* p_mol_editor = pmset->GetMolEditor(true);
		p_mol_editor->CalcHBonds(pmset);
	}
		  
	set<HaHBond, less<HaHBond> >::iterator hbnd_itr; //to scan hbonds list
	for(hbnd_itr = pmset->HBonds.begin(); hbnd_itr !=pmset->HBonds.end(); hbnd_itr++)
	{
		HaAtom* pt1; HaAtom* pt2; //the two atom pointers
		pt1 = (*hbnd_itr).dst;
		pt2 = (*hbnd_itr).src;
		if(pt1->Selected() && pt2->Selected()) //check if they are in selected set
		{
			double mbdist = Vec3D::CalcDistance(pt1,pt2,ANGSTROM_U);
			this->AddHarmConstraint(pt1,pt2,mbdist,force_const);
		}
	}

	return true;
}

int MolMechModel::AddHarmConstraint(HaAtom* aptr1,HaAtom* aptr2, double eq_dist, double force_const)
{
	if(aptr1 == NULL || aptr2 == NULL || aptr1 == aptr2 )
	{
		ErrorInMod("HaMolMechMod::AddHarmConstraint()"," Atoms of the constraint are not defined or are the same");
		return FALSE;
	}

	AtomContact cnt(aptr1,aptr2,AtomContactType::HARMONIC_CNT);
	cnt.SetParamsEneR(force_const,eq_dist);
	DistConstraints.push_back(cnt);

	p_amber_model->SetUpdateDataFlag();

	return TRUE;
}

int MolMechModel::SetHarmConstraint(HaAtom* aptr1,HaAtom* aptr2, double eq_dist, double force_const)
{
	if(aptr1 == NULL || aptr2 == NULL || aptr1 == aptr2 )
	{
		ErrorInMod("HaMolMechMod::SetHarmConstraint()"," Atoms of the constraint are not defined or are the same");
		return FALSE;
	}
	AtomContact* p_c = NULL;
	std::vector<AtomContact>::iterator citr;

	for( citr = DistConstraints.begin(); citr != DistConstraints.end(); citr++ )
	{
		AtomContact& cnt = *citr;
		if( !cnt.IsHarmonic() ) continue;
		if( cnt.pt1 != aptr1 )
		{
			if( cnt.pt1 != aptr2 ) continue; 
			if( cnt.pt2 != aptr1 ) continue;
		}
		else
		{
			if( cnt.pt2 != aptr2 ) continue;
		}
		p_c = &cnt;
		break;
	}
	if( p_c != NULL )
	{
		p_c->SetParamsEneR(force_const,eq_dist);
		p_amber_model->SetUpdateDataFlag();
		return TRUE;
	}
	return AddHarmConstraint(aptr1, aptr2, eq_dist, force_const);
}

int MolMechModel::SetHarmConstraint(const std::string& at_ref1, const::std::string& at_ref2, double eq_dist, double force_const)
{
	MolSet* pmset = this->GetMolSet();
	
	HaAtom* aptr1 = pmset->GetAtomByRef(at_ref1.c_str());
	HaAtom* aptr2 = pmset->GetAtomByRef(at_ref2.c_str());
	
	if( aptr1 == NULL ) 
	{
		PrintLog("Error in MolMechModel::SetHarmConstraint() \n");
		PrintLog("Atom Not found: %s \n",at_ref1.c_str());
		return FALSE;
	}

	if( aptr2 == NULL ) 
	{
		PrintLog("Error in MolMechModel::SetHarmConstraint() \n");
		PrintLog("Atom Not found: %s \n",at_ref2.c_str());
		return FALSE;
	}

	return SetHarmConstraint( aptr1, aptr2, eq_dist, force_const );
}

bool MolMechModel::ClearConstraints()
{
	DistConstraints.clear();
	p_amber_model->SetDistConstrData();

	return true;
}

int MolMechModel::UpdateConstraints()
{
//	PrintLog(" MolMechModel::UpdateConstraints() pt 1 \n");
	p_amber_model->SetDistConstrData();
//	if( p_mm_mod->p_amber_driver->numtasks > 1 ) HaMolMechMod::CallMMFunctionOnSlaves(MM_UPDATE_CONSTR_2);
	if( pApp->mpi_driver->nprocs > 1 ) HaMolMechMod::CallMMFunctionOnSlaves(MM_UPDATE_CONSTR_2);
	UpdateConstraints_2();
	return TRUE;
}

void MolMechModel::SetFFTGridsPerAng( double fft_grids_per_ang_par )
{
	fft_grids_per_ang = fft_grids_per_ang_par;
}


bool MolMechModel::CalcNonBondPt(const HaAtom* pt1, const HaAtom* pt2, 
			                     double& vdw_at_ene, double& el_at_ene)
{
	double r2 = 0.0;
	double r6, rr;
		
	double tmp=pt1->GetX()- pt2->GetX();
	r2 += tmp*tmp;
	tmp=pt1->GetY()- pt2->GetY();
	r2 += tmp*tmp;
	tmp=pt1->GetZ()- pt2->GetZ();
	r2 += tmp*tmp;
	
	if( calc_electr_flag )
		rr = sqrt(r2);
	
	r6= r2*r2*r2;
	
	if( calc_vdw_flag) 
	{
		double sigma6 = pt1->GetVdWRad() + pt2->GetVdWRad();
		sigma6*= (sigma6*sigma6);
		sigma6*= sigma6; 
		tmp = sigma6/r6;
		
		double estar = sqrt(pt1->GetVdWEne() * pt2->GetVdWEne());
		
		// E_vdW= sqrt(ew_i* ew_j)( (sigma_i+sigma_j)^12/r^12 - 
		//                         	 2.0* (sigma_i+sigma_j)^6/r^6 )
		
		vdw_at_ene = tmp*tmp;
		
		if(calc_vdw_flag == CALC_VDW_NORMAL)
		{
			vdw_at_ene = estar*( tmp*tmp - 2.0* tmp);
		}
		else if(calc_vdw_flag == CALC_VDW_NO_ATTRACT )
		{
			vdw_at_ene = estar*tmp*tmp;
		}
	}
	else
		vdw_at_ene = 0.0;
	
	if( calc_electr_flag )
	{
		double ch1 = pt1->GetCharge();
		double ch2 = pt2->GetCharge();
		
		if( ( fabs(ch1) > 1e-6 ) && ( fabs(ch2) > 1e-6 ) )
		{
			el_at_ene = HARTREE_TO_KCAL*ch1*ch2/(diel_const*rr);
			
//			PrintLog(" contr to el.ene ch1= %9.3f, ch2= %9.3f, dist= %9.3f, at_ene = %9.3f \n",
//				ch1, ch2, rr, el_at_ene );
		}
		else
			el_at_ene = 0.0;
	}	
	else
		el_at_ene = 0.0;

	return true;
}




int MolMechModel::SetBoundaryBox(double offset)
{
	double min_size_x =10.0e9 ;
	double min_size_y =10.0e9;
	double min_size_z =10.0e9;
	double max_size_x =-10.0e9;
	double max_size_y =-10.0e9;
	double max_size_z =-10.0e9;
	double box_size_x =0 ;
	double box_size_y =0 ;
	double box_size_z =0 ;
	double rad_vdw =0 ;
	
	if( !pmset->per_bc->IsSet() )
	{
		HaAtom* aptr;
		AtomIteratorMolSet aitr(pmset);
		for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			if(aptr ->GetX() < min_size_x)
			{
				min_size_x =aptr ->GetX();
				rad_vdw = 	aptr ->radius;
			}
			if(aptr ->GetY() < min_size_y) min_size_y =aptr ->GetY() ;
			if(aptr ->GetZ() < min_size_z) min_size_z =aptr ->GetZ() ;
			if(aptr ->GetX() > max_size_x) max_size_x =aptr ->GetX() ;
			if(aptr ->GetY() > max_size_y) max_size_y =aptr ->GetY() ;
			if(aptr ->GetZ() > max_size_z) max_size_z =aptr ->GetZ() ;
		}

		double px = max_size_x - min_size_x + rad_vdw + 2*offset;
		double py = max_size_y - min_size_y + rad_vdw + 2*offset;
		double pz = max_size_z - min_size_z + rad_vdw + 2*offset;

		pmset->per_bc->SetBox(px,py,pz);
	}
	return TRUE;
}

void MolMechModel::SetOmitInteractionsParam(const OmitInteractionsParam& omit_interactions_new)
{
	omit_interactions = omit_interactions_new;
}

void MolMechModel::SetMMElectrMethod( const MMElectrMethod& electr_method_new  )
{
	electr_method = electr_method_new;
}

void  MolMechModel::SetCalcDirectInter( bool set_par )
{
	p_amber_model->do_direct =  ( set_par ? 1 : 0 );
}

void MolMechModel::SetCalcRecipSpaceInter( bool set_par )
{
	p_amber_model->do_recip = ( set_par ? 1 : 0 );
}

void MolMechModel::SetCalcSelfInter( bool set_par )
{
	p_amber_model->do_self = ( set_par ? 1 : 0 );
}

void MolMechModel::SetCalcAdjustInter( bool set_par )
{
	p_amber_model->do_adjust = ( set_par ? 1 : 0 );
}

void MolMechModel::SetCalcVdWInter( bool set_par )
{
	p_amber_model->do_vdw = ( set_par ? 1 : 0 );
}

void MolMechModel::SetCalcInducedInter( bool set_par )
{
	p_amber_model->do_induced = ( set_par ? 1 : 0 );
}

double MolMechModel::GetScale14Electr() const
{
	return scale_14_electr;
}

double MolMechModel::GetScale14VdW() const
{
	return scale_14_vdw;
}

double MolMechModel::GetNBCutDist() const
{	
	return nb_cut_dist;
}

double MolMechModel::GetDielConst() const
{
	return diel_const;
}

double MolMechModel::GetIonStrength() const
{
	return ion_strength;
}

int MolMechModel::GetDipoleScfIterMax() const
{
	return dipole_scf_iter_max;
}

double MolMechModel::GetDipoleScfTol() const
{
	return dipole_scf_tol;
}
	
double MolMechModel::GetEEDsumCut() const
{
	return ee_dsum_cut;
}

double MolMechModel::GetEEDampedCut() const
{
	return ee_damped_cut;
}

double MolMechModel::GetSorCoef() const
{
	return sor_coefficient;
}
	
double MolMechModel::GetTholeExponCoef() const
{
	return thole_expon_coeff;
}
	
double MolMechModel::GetVdwTaper() const
{
	return vdw_taper;
}

void MolMechModel::SetScale14Electr( double scale_14_electr_new )
{
	scale_14_electr = scale_14_electr_new;
	p_amber_model->scee  = scale_14_electr;

}

void MolMechModel::SetScale14VdW( double scale_14_vdw_new )
{
	scale_14_vdw = scale_14_vdw_new;
	p_amber_model->scnb  = scale_14_vdw;
}

void MolMechModel::SetScale14Inter( double scale_14)
{
	SetScale14Electr(scale_14);
	SetScale14VdW(scale_14);
}
	
void MolMechModel::SetNBCutDist( double nb_cut_dist_new )
{
	nb_cut_dist = nb_cut_dist_new;
	p_amber_model->gb_cutoff = nb_cut_dist; 
	p_amber_model->es_cutoff  = nb_cut_dist;   
	p_amber_model->vdw_cutoff = nb_cut_dist; 
}

void MolMechModel::SetDielConst( double diel_const_new )
{
	diel_const = diel_const_new;
	p_amber_model->dielc = diel_const;
}

void MolMechModel::SetIonStrength( double ion_strength_new )
{
	ion_strength = ion_strength_new;
}

void MolMechModel::SetDipoleScfIterMax( int dipole_scf_iter_max_new )
{
	dipole_scf_iter_max = dipole_scf_iter_max_new;
}

void MolMechModel::SetDipoleScfTol( double dipole_scf_tol_new ) 
{
	dipole_scf_tol = dipole_scf_tol_new;
}

void MolMechModel::SetEEDsumCut( double ee_dsum_cut_new )      
{
	ee_dsum_cut = ee_dsum_cut_new;
}

void MolMechModel::SetEEDampedCut( double ee_damped_cut_new )
{
	ee_damped_cut = ee_damped_cut_new;
}

void MolMechModel::SetSorCoef( double sor_coef_new )
{
	sor_coefficient = sor_coef_new;
}

void MolMechModel::SetTholeExponCoef( double thole_expon_coeff_new )
{
	thole_expon_coeff = thole_expon_coeff_new;
}

void MolMechModel::SetVdwTaper( double vdw_taper_new )
{
	vdw_taper = vdw_taper_new;
}

double MolMechModel::GetGBAtomRad(HaAtom* aptr, int gb_param_type)
{
	double gb_rad = 1.0;
	if( aptr == NULL) return gb_rad;
	int elem = aptr->GetElemNo();
//	std::string ff_symbol = aptr->GetFFSymbol();
//	boost::to_upper(ff_symbol);
//	if( ff_symbol == "DC" || ff_symbol == "DO" || ff_symbol == "DN" || ff_symbol == "DH")
//	{
//		gb_rad = 0.0;
//		return gb_rad;
//	}

	AtomGroup bonded_atoms;
	aptr->GetBondedAtoms(bonded_atoms);

	if( gb_param_type == BONDI_GB_PARAM || gb_param_type == AMBER6_BONDI_GB_PARAM || 
		gb_param_type == MOD_BONDI_GB_PARAM  || gb_param_type == HN_MOD_BONDI_GB_PARAM ) 
	{
         // Bondi or modified Bondi radii 
		switch( elem ) 
		{
		case  1: gb_rad = 1.2; 
			/* make the modifications that hydrogen radii
			depend upon the atoms they are bonded to.  
			gb_param_type = AMBER6_BONDI_GB_PARAM corresponds to Amber 6, JACS 122:2489 (2000);
			gb_param_type = MOD_BONDI_GB_PARAM  adds the update of Biopolymers 56: 275 (2001)   
			*/

			if( bonded_atoms.size() > 0 ) 
			{
				/* For multiply bonded Hydrogen atoms use the first
				* bond for determining modified GB radii.
				* WAT contains multiply bonded Hydrogen atoms 
				* so do not emit a warning.
				*/
				HaAtom* aptr_b = bonded_atoms[0];

				if( gb_param_type == AMBER6_BONDI_GB_PARAM  || gb_param_type == MOD_BONDI_GB_PARAM ) 
				{
					if( aptr_b->GetElemNo() == 6  ) gb_rad = 1.3;
					if( aptr_b->GetElemNo() == 8  ) gb_rad = 0.8;
					if( aptr_b->GetElemNo() == 16 ) gb_rad = 0.8;
					if( aptr_b->GetElemNo() == 7  &&
						gb_param_type == MOD_BONDI_GB_PARAM) gb_rad = 1.3;
				}
				else if( gb_param_type == HN_MOD_BONDI_GB_PARAM ) 
				{ 
					if( aptr_b->GetElemNo() == 7 ) gb_rad = 1.3;
				}
			}
			else 
			{
				std::string at_ref = aptr->GetRef();
				PrintLog("Warning in HaMolMechMod::GetGBAtomRad() \n");
				PrintLog("Hydrogen atom %s is unbonded \n",at_ref.c_str());
				PrintLog("Can not determine unmodified Bondi GB radius \n");
			}
			break;
		case  6: gb_rad = 1.7; break;
		case  7: gb_rad = 1.55; break;
		case  8: gb_rad = 1.5; break;
		case  9: gb_rad = 1.5; break;
		case 14: gb_rad = 2.1; break;
		case 15: gb_rad = 1.85; break;
		case 16: gb_rad = 1.8; break;
		case 17: gb_rad = 1.7; break;
		default: gb_rad = 1.5; break;
		}
    } 
	else if ( gb_param_type == HUO_KOLLMAN_GB_PARAM ) {  // radii from Huo & Kollman 
            switch( elem ) {
                case  1: gb_rad = 1.15; break;
                case  6: gb_rad = 1.85; break;
                case  7: gb_rad = 1.64; break;
                case  8: gb_rad = 1.53; break;
                case  9: gb_rad = 1.53; break;
                case 15: gb_rad = 2.02; break;
                case 16: gb_rad = 2.00; break;
                case 17: gb_rad = 1.97; break;
                case 35: gb_rad = 2.03; break;
                case 53: gb_rad = 2.10; break;
                default: gb_rad = 1.5;  break;  
            }
        }
	return gb_rad;
}

double MolMechModel::GetGBAtomScreening(HaAtom* aptr, int gb_param_type)
{
	double gb_screen = 1.0;
	if( aptr == NULL) return gb_screen;
	int elem = aptr->GetElemNo();
	std::string ff_symbol = aptr->GetFFSymbol();
	boost::to_upper(ff_symbol);
	if( ff_symbol == "DC" || ff_symbol == "DO" || ff_symbol == "DN" || ff_symbol == "DH")
	{
		gb_screen = 0.0;
		return gb_screen;
	}

	if( gb_param_type == BONDI_GB_PARAM || gb_param_type == AMBER6_BONDI_GB_PARAM || 
		gb_param_type == MOD_BONDI_GB_PARAM  || gb_param_type == HUO_KOLLMAN_GB_PARAM ||
		gb_param_type == HN_MOD_BONDI_GB_PARAM )
	{
		/* for now, hardwire the Bondi radii  */
		switch( elem ){
				case 1:  gb_screen = 0.85; break;
				case 6:  gb_screen = 0.72; break;
				case 7:  gb_screen = 0.79; break;
				case 8:  gb_screen = 0.85; break;
				case 9:  gb_screen = 0.88; break;
				case 15: gb_screen = 0.86; break;
				case 16: gb_screen = 0.96; break;
				default: gb_screen = 0.8; break;  
		}
	} 
	else if( gb_param_type == JAYARAM_GB_PARAM )
	{ /* param for Jayaram et al. 'GB' */
		switch( elem ){
				case 1:  gb_screen = 0.8461; break;
				case 6:  gb_screen = 0.9615; break;
				case 7:  gb_screen = 0.9343; break;
				case 8:  gb_screen = 1.0088; break;
				case 11: gb_screen = 1.0000; break; 
				case 12: gb_screen = 1.0000; break; /* set by HG */
				case 15: gb_screen = 1.0700; break; 
				case 16: gb_screen = 1.1733; break; 
				default: gb_screen = 0.8000; break; /* set by HG */
		}
	} 
	else if ( gb_param_type == MOD_JAYARAM_GB_PARAM )
	{  /* param for Jayaram et al. 'MGB' */
		switch( elem ){
				case 1:  gb_screen = 0.8846; break;
				case 6:  gb_screen = 0.9186; break;
				case 7:  gb_screen = 0.8733; break;
				case 8:  gb_screen = 0.8836; break;
				case 11: gb_screen = 1.0000; break;
				case 12: gb_screen = 1.0000; break; /* set by HG */
				case 15: gb_screen = 0.9604; break;
				case 16: gb_screen = 0.9323; break;
				default: gb_screen = 0.8000; break; /* set by HG */
		}
	}
	return gb_screen;
}

int MolMechModel::UpdateDataFromFort()
{
	if( IsAmoebaFF() ) 
	{
		GetIndDipolesFromFort();
		GetMultipolesFromFort();
	}
	return TRUE;
}

int MolMechModel::GetIndDipolesFromFort()
{
	if( !IsAmoebaFF() ) return FALSE;

	int na = this->GetNA();
	if( ind_dip_p.size() != 3*na ) 
	{
		ind_dip_p.resize(3*na);
		ind_dip_p = 0.0;
	}
	
	if( ind_dip_d.size() != 3*na ) 
	{
		ind_dip_d.resize(3*na);
		ind_dip_d = 0.0;
	}
	
	FC_FUNC_MODULE(amoeba_induced_mod,get_ind_dip)( &na, ind_dip_d.v(), ind_dip_p.v() );
	
	return TRUE;
}

int MolMechModel::GetMultipolesFromFort()
{
	if( !IsAmoebaFF() ) return FALSE; 

	int na = this->GetNA();
	if( global_multipole.size() != 10*na ) 
	{
		global_multipole.resize(10*na);
		global_multipole = 0.0;
	}
	FC_FUNC_MODULE(amoeba_multipoles_mod,get_global_multipole)( &na, global_multipole.v() );
	
	return TRUE;
}

void MolMechModel::PrintIndDipoles()
{
	MolSet* pmset = this->GetMolSet();
	int na = pmset->GetNAtoms();

	if( 3*na != ind_dip_p.size() )
	{
		PrintLog(" Error in MolMechModel::PrintIndDipoles() \n");
		PrintLog(" Incorrect size of ind_dip_p array %d \n",ind_dip_p.size() );
		return;
	}

	if( 3*na != ind_dip_d.size() )
	{
		PrintLog(" Error in MolMechModel::PrintIndDipoles() \n");
		PrintLog(" Incorrect size of ind_dip_d array %d \n",ind_dip_d.size() );
		return;
	}

	PrintLog("\n ind_dip_p: \n");
	int i = 0;
	AtomIteratorMolSet aitr(pmset); 
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom() )
	{
		PrintLog(" %s %12.6f %12.6f %12.6f \n", aptr->GetRef().c_str(), ind_dip_p[3*i], ind_dip_p[3*i+1],ind_dip_p[3*i+2]);
		i++;
	}
	PrintLog(" \n");

	i = 0;
	PrintLog("\n ind_dip_d: \n");
	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom() )
	{
		PrintLog(" %s %12.6f %12.6f %12.6f \n", aptr->GetRef().c_str(), ind_dip_d[3*i], ind_dip_d[3*i+1],ind_dip_d[3*i+2]);
		i++;
	}
	PrintLog(" \n");
}

void MolMechModel::PrintMultipoles()
{
	MolSet* pmset = this->GetMolSet();
	int na = pmset->GetNAtoms();

	if( 10*na != global_multipole.size() )
	{
		PrintLog(" Error in MolMechModel::PrintMultipoles() \n");
		PrintLog(" Incorrect size of global_multipole array %d \n",global_multipole.size() );
		return;
	}

	PrintLog("\n charges and dipoles: \n");
	int i = 0;
	AtomIteratorMolSet aitr(pmset); 
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom() )
	{
		PrintLog(" %s  %9.4f     %9.4f %9.4f %9.4f \n", aptr->GetRef().c_str(), 
			 global_multipole[10*i], global_multipole[10*i+1], global_multipole[10*i+2], global_multipole[10*i+3]);
		i++;
	}
	PrintLog(" \n");

	PrintLog("\n quadrupoles: \n");
	i = 0;
	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom() )
	{
		PrintLog(" %s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n", aptr->GetRef().c_str(), 
			 global_multipole[10*i+4], global_multipole[10*i+5], global_multipole[10*i+6], 
			 global_multipole[10*i+7], global_multipole[10*i+8], global_multipole[10*i+9]);
		i++;
	}
	PrintLog(" \n");
}

void MolMechModel::PrintTotMultipoles(const Vec3D* pt_orig )
{
	double x_orig = 0.0; 
	double y_orig = 0.0;
	double z_orig = 0.0;

	if( pt_orig != NULL ) 
	{
		x_orig = pt_orig->GetX();
		y_orig = pt_orig->GetY();
		z_orig = pt_orig->GetZ();
	}

	double ch_tot = 0.0;
	HaVec_double dipole(3); dipole = 0.0;
	HaVec_double ind_dipole(3); ind_dipole = 0.0;
	HaVec_double qpole(6);  qpole  = 0.0;

	MolSet* pmset = this->GetMolSet();
	int na = pmset->GetNAtoms();

	if( 10*na != global_multipole.size() )
	{
		PrintLog(" Error in MolMechModel::PrintMultipolesTot() \n");
		PrintLog(" Incorrect size of global_multipole array %d \n",global_multipole.size() );
		return;
	}

	if( 3*na != ind_dip_p.size() )
	{
		PrintLog(" Error in MolMechModel::PrintMultipolesTot() \n");
		PrintLog(" Incorrect size of ind_dip_p array %d \n",ind_dip_p.size() );
		return;
	}

	if( 3*na != ind_dip_d.size() )
	{
		PrintLog(" Error in MolMechModel::PrintMultipolesTot() \n");
		PrintLog(" Incorrect size of ind_dip_d array %d \n",ind_dip_d.size() );
		return;
	}

	int i = 0;
	AtomIteratorMolSet aitr(pmset); 
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom(), i++ )
	{
		double x = aptr->GetX() - x_orig;
		double y = aptr->GetY() - y_orig;
		double z = aptr->GetZ() - z_orig;

		double ch = global_multipole[10*i];
		ch_tot += ch;
		
		ind_dipole[0] += ind_dip_p[3*i]   + ind_dip_d[3*i];
		ind_dipole[1] += ind_dip_p[3*i+1] + ind_dip_d[3*i+1];
		ind_dipole[2] += ind_dip_p[3*i+2] + ind_dip_d[3*i+2];

		dipole[0] += ch*x + global_multipole[10*i+1] + ind_dip_p[3*i]     + ind_dip_d[3*i];
		dipole[1] += ch*y + global_multipole[10*i+2] + ind_dip_p[3*i+1]   + ind_dip_d[3*i+1];
		dipole[2] += ch*z + global_multipole[10*i+3] + ind_dip_p[3*i+2]   + ind_dip_d[3*i+2];

		qpole[0] += global_multipole[10*i+4];
		qpole[1] += global_multipole[10*i+5];
		qpole[2] += global_multipole[10*i+6];
		qpole[3] += global_multipole[10*i+7];
		qpole[4] += global_multipole[10*i+8];
		qpole[5] += global_multipole[10*i+9];
	}

	PrintLog(" \n");
	PrintLog(" charge = %12.6f", ch_tot );
	PrintLog(" tot dipole:     %12.6f %12.6f %12.6f  debye \n",
		       dipole[0]*EL_ANG_TO_DEBYE,dipole[1]*EL_ANG_TO_DEBYE,dipole[2]*EL_ANG_TO_DEBYE);
	PrintLog(" induced dipole: %12.6f %12.6f %12.6f  debye \n",
		       ind_dipole[0]*EL_ANG_TO_DEBYE,ind_dipole[1]*EL_ANG_TO_DEBYE,ind_dipole[2]*EL_ANG_TO_DEBYE);
	PrintLog(" quadrupoles: %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n",
		       qpole[0],qpole[1],qpole[2],qpole[3],qpole[4],qpole[5]);
	PrintLog(" \n");	
}

Vec3D MolMechModel::GetTotIndDipole(const Vec3D* pt_orig )
{
	Vec3D ind_dipole = GetTotIndDipole1(pt_orig) + GetTotIndDipole2(pt_orig);
	ind_dipole.Scale(0.5);
	
	return ind_dipole;
}

Vec3D MolMechModel::GetTotIndDipole1(const Vec3D* pt_orig )
{
	double x_orig = 0.0; 
	double y_orig = 0.0;
	double z_orig = 0.0;

	if( pt_orig != NULL ) 
	{
		x_orig = pt_orig->GetX();
		y_orig = pt_orig->GetY();
		z_orig = pt_orig->GetZ();
	}

	Vec3D ind_dipole; ind_dipole = 0.0;

	MolSet* pmset = this->GetMolSet();
	int na = pmset->GetNAtoms();

	if( 3*na != ind_dip_d.size() )
	{
		PrintLog(" Error in MolMechModel::GetTotIndDipole1() \n");
		PrintLog(" Incorrect size of ind_dip_d array %d \n",ind_dip_d.size() );
		return ind_dipole;
	}

	int i = 0;
	AtomIteratorMolSet aitr(pmset); 
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom(), i++ )
	{
		double x = aptr->GetX() - x_orig;
		double y = aptr->GetY() - y_orig;
		double z = aptr->GetZ() - z_orig;

		double ch = global_multipole[10*i];
		
		ind_dipole[0] += ind_dip_d[3*i];
		ind_dipole[1] += ind_dip_d[3*i+1];
		ind_dipole[2] += ind_dip_d[3*i+2];
	}
	return ind_dipole;
}

Vec3D MolMechModel::GetTotIndDipole2(const Vec3D* pt_orig )
{
	double x_orig = 0.0; 
	double y_orig = 0.0;
	double z_orig = 0.0;

	if( pt_orig != NULL ) 
	{
		x_orig = pt_orig->GetX();
		y_orig = pt_orig->GetY();
		z_orig = pt_orig->GetZ();
	}

	Vec3D ind_dipole; ind_dipole = 0.0;

	MolSet* pmset = this->GetMolSet();
	int na = pmset->GetNAtoms();

	if( 3*na != ind_dip_p.size() )
	{
		PrintLog(" Error in MolMechModel::GetTotIndDipole2() \n");
		PrintLog(" Incorrect size of ind_dip_p array %d \n",ind_dip_p.size() );
		return ind_dipole;
	}

	int i = 0;
	AtomIteratorMolSet aitr(pmset); 
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom(), i++ )
	{
		double x = aptr->GetX() - x_orig;
		double y = aptr->GetY() - y_orig;
		double z = aptr->GetZ() - z_orig;

		double ch = global_multipole[10*i];
		
		ind_dipole[0] += ind_dip_p[3*i];
		ind_dipole[1] += ind_dip_p[3*i+1];
		ind_dipole[2] += ind_dip_p[3*i+2];
	}
	return ind_dipole;
}


Vec3D MolMechModel::GetTotDipole(const Vec3D* pt_orig )
{
	double x_orig = 0.0; 
	double y_orig = 0.0;
	double z_orig = 0.0;

	if( pt_orig != NULL ) 
	{
		x_orig = pt_orig->GetX();
		y_orig = pt_orig->GetY();
		z_orig = pt_orig->GetZ();
	}

	Vec3D dipole; dipole = 0.0;

	MolSet* pmset = this->GetMolSet();
	int na = pmset->GetNAtoms();

	if( 10*na != global_multipole.size() )
	{
		PrintLog(" Error in MolMechModel::GetTotDipole() \n");
		PrintLog(" Incorrect size of global_multipole array %d \n",global_multipole.size() );
		return dipole;
	}

	//if( 3*na != ind_dip_p.size() )
	//{
	//	PrintLog(" Error in MolMechModel::GetTotDipole() \n");
	//	PrintLog(" Incorrect size of ind_dip_p array %d \n",ind_dip_p.size() );
	//	return dipole;
	//}

	if( 3*na != ind_dip_d.size() )
	{
		PrintLog(" Error in MolMechModel::GetTotDipole() \n");
		PrintLog(" Incorrect size of ind_dip_d array %d \n",ind_dip_d.size() );
		return dipole;
	}

	int i = 0;
	AtomIteratorMolSet aitr(pmset); 
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom(), i++ )
	{
		double x = aptr->GetX() - x_orig;
		double y = aptr->GetY() - y_orig;
		double z = aptr->GetZ() - z_orig;

		double ch = global_multipole[10*i];
		
		dipole[0] += ch*x + global_multipole[10*i+1] + ind_dip_d[3*i];
		dipole[1] += ch*y + global_multipole[10*i+2] + ind_dip_d[3*i+1];
		dipole[2] += ch*z + global_multipole[10*i+3] + ind_dip_d[3*i+2];
	}
	return dipole;
}

HaVec_double MolMechModel::GetTotQpole(const Vec3D* pt_orig )
{
	double x_orig = 0.0; 
	double y_orig = 0.0;
	double z_orig = 0.0;

	if( pt_orig != NULL ) 
	{
		x_orig = pt_orig->GetX();
		y_orig = pt_orig->GetY();
		z_orig = pt_orig->GetZ();
	}

	HaVec_double qpole(6);  qpole  = 0.0;

	MolSet* pmset = this->GetMolSet();
	int na = pmset->GetNAtoms();

	if( 10*na != global_multipole.size() )
	{
		PrintLog(" Error in MolMechModel::PrintMultipolesTot() \n");
		PrintLog(" Incorrect size of global_multipole array %d \n",global_multipole.size() );
		return qpole;
	}

	if( 3*na != ind_dip_p.size() )
	{
		PrintLog(" Error in MolMechModel::PrintMultipolesTot() \n");
		PrintLog(" Incorrect size of ind_dip_p array %d \n",ind_dip_p.size() );
		return qpole;
	}

	if( 3*na != ind_dip_d.size() )
	{
		PrintLog(" Error in MolMechModel::PrintMultipolesTot() \n");
		PrintLog(" Incorrect size of ind_dip_d array %d \n",ind_dip_d.size() );
		return qpole;
	}

	int i = 0;
	AtomIteratorMolSet aitr(pmset); 
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom(), i++ )
	{
		double x = aptr->GetX() - x_orig;
		double y = aptr->GetY() - y_orig;
		double z = aptr->GetZ() - z_orig;

		double ch = global_multipole[10*i];
		HaVec_double dipole(3); 
		HaVec_double ind_dipole(3); 
		
		ind_dipole[0] = ind_dip_p[3*i]   + ind_dip_d[3*i];
		ind_dipole[1] = ind_dip_p[3*i+1] + ind_dip_d[3*i+1];
		ind_dipole[2] = ind_dip_p[3*i+2] + ind_dip_d[3*i+2];

		dipole[0] += ch*x + global_multipole[10*i+1] + ind_dipole[0];
		dipole[1] += ch*y + global_multipole[10*i+2] + ind_dipole[1];
		dipole[2] += ch*z + global_multipole[10*i+3] + ind_dipole[2];

		qpole[0] += global_multipole[10*i+4] + dipole[0]*x;
		qpole[1] += global_multipole[10*i+5] + dipole[1]*y;
		qpole[2] += global_multipole[10*i+6] + dipole[2]*z;
		qpole[3] += global_multipole[10*i+7] + dipole[0]*y + dipole[1]*x;
		qpole[4] += global_multipole[10*i+8] + dipole[0]*z + dipole[2]*x;
		qpole[5] += global_multipole[10*i+9] + dipole[1]*z + dipole[2]*y;
	}
	return qpole;
}

void MolMechModel::SetUseMortLib( int set_par )
{
	setup_params_from_mort_flag = set_par;
}

void  MolMechModel::SetPMECoef( double pme_ew_coeff_par )
{
	pme_ew_coeff = pme_ew_coeff_par;
}
