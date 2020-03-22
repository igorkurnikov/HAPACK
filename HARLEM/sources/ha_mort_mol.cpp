/*! \file ha_mort_mol.cpp

   Classes and Function to define interface bewtween MORT library and HARLEM.

  \author Igor Kurnikov 
  \date 2010-

*/

#include <mpi.h>

#include <boost/format.hpp>

#include "object.hpp"
#include "format.hpp"
#include "pdbent.hpp"
#include "hamolset.h"
#include "mm_force_field.h"
#include "ha_mort_mm.h"

void MMForceField::switch_amoeba_type( mort::molecule_t& mol, const mort::molecule_t& poleff )
//
// atom in AMOEBA has 2 integer force field types atomid (used to set parameters for valence bonds,angles etc)
// and poleid - used to set charges, multipoles and polarizabilities of atoms
//
// when loaded from the residue database : residue template has poleid as TYPE string parameter
// this function finds atomid corresponding to poleid from poleff
// string  parameter TYPE      is set to string representation of atomid
// string  parameter POLTYPE   is set to string representation of poleid
// integer parameter POLTYPEID is set to poleid
{
	int switched = 0;
	if( mol.get_i("amoeba-switched", switched) && switched>0 )
	{
		return; 
	}
	char buf[120];

	mort::atomiter_t ai = mol.atom_begin();
	for( ; ai != mol.atom_end(); ++ai )
	{
		std::string pole_type = ai->get_s( mort::TYPE );
		int pole_type_id = atoi( pole_type.c_str() );

		if( pole_type_id <= 0 )
		{
			std::string str_err = "Error in switch_amoeba_type() \n";
			sprintf(buf," atom %s with ID %d has invalid pole_type_d = %d \n",ai->get_s( mort::NAME ).c_str(),ai->get_s( mort::ID ).c_str() );
			str_err += buf;
			throw std::runtime_error(str_err);
		}

		mort::atom_t pole( poleff, pole_type_id-1 );

		int atom_type_id = pole.get_i(mort::TYPEID);
		sprintf(buf,"%d",atom_type_id);
		std::string atom_type = buf;

		ai->set_s( mort::TYPE, atom_type );
		ai->set_s( mort::POLTYPE, pole_type );
		ai->set_i( mort::POLTYPEID, pole_type_id );
	}

	mol.set_i("amoeba-switched", 1);
}

int MolSet::SetMortMol(mort::molecule_t& mort_mol, const ForceFieldType& ff_type )
{
	MMForceField* p_ff = MMForceField::GetMMForceField(ff_type,TRUE);
	if(p_ff == NULL) return FALSE;

	if( !p_ff->IsMortFFInitiated() ) 
	{
		PrintLog(" Error in MolSet::SetMortMol() \n");
		PrintLog(" MORT Force field structure of the force field %s are not initiated \n",ff_type.label());
		return FALSE;
	}
	char buf[120];
	mort::molecule_t m;
	
	ChainIteratorMolSet chitr(this);
	HaChain* p_chain;

	AtomIntMap at_idx_map;

	int idx_at = 0;

	try
	{
		for(p_chain = chitr.GetFirstChain(); p_chain; p_chain = chitr.GetNextChain())
		{
			ResidueIteratorChain ritr(p_chain);
			HaResidue* pres;

			int first_res = TRUE;

			for( pres = ritr.GetFirstRes();;)
			{
				mort::resd_t resd = mort::resd_t::create( m );

				std::string rtype = MMForceField::GetAmberResName( pres->GetFullName() );

				std::string rname = rtype;
				sprintf(buf,"%d",pres->GetSerNo());
				rname += buf;
				char rchain       = pres->GetHostChain()->ident;

				resd.set_s(mort::NAME, rname );
				resd.set_s(mort::TYPE, rtype );
				resd.set_i(mort::CHAIN, rchain );

				if( first_res )
				{
					resd.set_i( mort::HEAD,  0 );
					resd.set_i( mort::AAPOS, mort::NTERM );
					first_res = FALSE;
				}
				else
				{
					resd.set_i( mort::AAPOS, mort::OTHER );
				}

				AtomIteratorAtomGroup aitr(pres);
				HaAtom* aptr;

				for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
				{
					at_idx_map[aptr] = idx_at;
					idx_at++;

					std::string name = MMForceField::GetAmberAtName(aptr->GetName(),pres->GetFullName());
					int atomseq      = aptr->GetSerNo();
					int element      = aptr->GetElemNo(); 

					double crdx = aptr->GetX();
					double crdy = aptr->GetY();
					double crdz = aptr->GetZ();

					mort::atom_t a = mort::atom_t::create(resd, name);
					a.set_i( mort::ELEMENT,  element );
					a.set_i( mort::SEQUENCE, atomseq );
					a.set_v( mort::POSITION, mort::makevec(crdx,crdy,crdz) );
				}

				pres = ritr.GetNextRes();
				if( pres != NULL) continue;

				int pos = resd.get_i( mort::AAPOS );
				if( pos == mort::NTERM )
				{
					resd.set_i( mort::AAPOS, mort::ALONE );
				}
				else
				{
					resd.set_i( mort::AAPOS, mort::CTERM );
				}

				resd.set_i( mort::TAIL, 0 );
				break;
			}
		}

		BondIteratorMolSet bitr(this);
		HaBond* bptr;

		mort::atom_range atoms_mort = m.atoms();
		
		for( bptr = bitr.GetFirstBond(); bptr; bptr = bitr.GetNextBond() )
		{
			HaAtom* aptr1 = bptr->srcatom;
			HaAtom* aptr2 = bptr->dstatom;

			if( at_idx_map.count(aptr1) == 0 ) continue;
			if( at_idx_map.count(aptr2) == 0 ) continue;

			int idx1 = at_idx_map[aptr1];
			int idx2 = at_idx_map[aptr2];

			mort::atom_t a1 = atoms_mort[idx1];
			mort::atom_t a2 = atoms_mort[idx2];

			mort::bond_t b = mort::bond_t::create( a1, a2 );
			b.set_i( mort::ORDER, 1 );
		}

		mort::database_t* p_mdb = p_ff->p_mort_ff->p_mdb;

//		if( !p_mdb->has("_namemap") )
//		{
//			PrintLog("Warning in MolSet::SetMortMol() \n");
//			PrintLog("cnamemap does not exist when modelize resd \n");
//		}

		map<int,int> idmap;
		mort::molecule_t full;

		mort::resd_t rprev(full, -1);
		mort::resditer_t ri = m.resd_begin();
		for( ; ri != m.resd_end(); ++ri )
		{
			mort::resd_t rcurt = mort::mdlize_resd( *ri, *(p_mdb), full, idmap );
			mort::conect_resd( rprev, rcurt );
			rprev = rcurt;
		}

		m.swap(full);

		if( per_bc->IsSet() )
		{
			double a = per_bc->GetA();
			double b = per_bc->GetB();
			double c = per_bc->GetC();

			double beta = per_bc->GetBeta();

			m.set_v(mort::BOX, mort::makevec(a, b, c, RAD_TO_DEG*beta) );
			m.set_i(mort::SOLUTE, mort::BOX);
		}

		int mol_sol_par = m.get_i( mort::SOLUTE );

		mort_mol.swap(m);

		if( ff_type == ForceFieldType::AMOEBA )
		{
			MMForceField::switch_amoeba_type( mort_mol,*p_ff->p_mort_ff->p_poleff );
		}
	}
	catch(std::exception& ex) 
	{
		PrintLog(" Error in MolSet::SetMortMol() \n");
		PrintLog(" %s \n",ex.what());
		return FALSE;
	}

	return TRUE;
}

int MolSet::SavePDBMort(const char* fname)
{
	mort::molecule_t mol;
	int ires = SetMortMol(mol,ForceFieldType::AMOEBA);
	if( !ires) return FALSE;
	std::ofstream os( fname );
	
	mort::write_pdb(os, mol);
	return TRUE;
}

int MolSet::SaveSDFMort(const char* fname)
{
	mort::molecule_t mol;
	int ires = SetMortMol(mol,ForceFieldType::AMOEBA);
	if( !ires ) return FALSE;
	std::ofstream os( fname );
	
	mort::write_sdf( os, mol );
	return TRUE;
}

int HaResidue::CheckStructMortLib(const ForceFieldType& ff_type)
{
	MMForceField* p_ff = MMForceField::GetMMForceField(ff_type, TRUE );
	if( p_ff == NULL) return FALSE;
	int ires;
	if( !p_ff->IsMortFFInitiated() ) 
	{
		ires = p_ff->InitMortFF();
		if(!ires) return FALSE;
	}
	mort::database_t& mdb = *(p_ff->p_mort_ff->p_mdb);

	std::string res_name_mort = MMForceField::GetAmberResName ( this->GetFullName(), ff_type );
	mort::molecule_ptr p_templ;
	try
	{
		p_templ = mdb.get_mol( res_name_mort );
	}
	catch( std::exception& ex )
	{
		PrintLog(" Error in HaResidue::CheckStructMortLib() \n");
		PrintLog("%s \n",ex.what());
		return FALSE;
	}

	ires = TRUE;

	std::set< std::string, less<std::string> > res_at_names;
	std::set< std::string, less<std::string> > templ_at_names;

	mort::atomiter_t aitr_templ = p_templ->atom_begin();
    for( ; aitr_templ != p_templ->atom_end(); ++aitr_templ )
    {
		std::string at_name = (*aitr_templ).get_s(mort::NAME);
		templ_at_names.insert(at_name);
	}
	
	AtomIteratorAtomGroup aitr_r(this);
	HaAtom* aptr;
	for( aptr = aitr_r.GetFirstAtom(); aptr; aptr = aitr_r.GetNextAtom() )
	{
		std::string at_name = aptr->GetName();
		std::string at_name_mort = MMForceField::GetAmberAtName ( at_name, this->GetFullName(), ff_type );
		if( res_at_names.count(at_name_mort) > 0 )
		{
			PrintLog(" Residue %s  has multiple atoms with name = %s   (mort_lib atom name = %s ) \n",
				      (this->GetRef()).c_str(),at_name.c_str(), at_name_mort.c_str() );
			ires = FALSE;
		}
		else
		{
			res_at_names.insert(at_name_mort);
		}

		if( templ_at_names.count(at_name_mort) == 0 )
        {
             PrintLog(" Template for residue %s does not have atom with name %s \n", 
				        (this->GetRef()).c_str(), at_name_mort.c_str() );
			 ires = FALSE;
        }
	}

	std::set< std::string, less<std::string> >::iterator nm_itr;	
	for( nm_itr = templ_at_names.begin(); nm_itr != templ_at_names.end(); nm_itr++)
	{
		std::string at_name_mort = (*nm_itr);
		if( res_at_names.count( at_name_mort ) == 0 )
		{
			PrintLog(" Template for residue %s has an extra atom with name %s \n",
			         (this->GetRef()).c_str(), at_name_mort.c_str() );
			ires = FALSE;
		}
	}

	if( ires ) // Check if atoms in the same sequence as in the template:
	{
		aitr_templ = p_templ->atom_begin();
		aptr = aitr_r.GetFirstAtom();
		for( ; aptr; aptr = aitr_r.GetNextAtom(),++aitr_templ )
		{
			std::string at_name = aptr->GetName();
			std::string at_name_mort = MMForceField::GetAmberAtName ( at_name, this->GetFullName(), ff_type );
			std::string at_name_templ = (*aitr_templ).get_s(mort::NAME);
			if( at_name_mort != at_name_templ )
			{
				PrintLog(" Residue %s has atom order different from that in the template \n",(this->GetRef()).c_str());
				ires = FALSE;
				break;
			}
		}
	}

	return ires;   
}
