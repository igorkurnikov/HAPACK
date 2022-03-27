/*! \file hamolecule.cpp

   Functions of the class to define molecule object in HARLEM.

   \author Igor Kurnikov 
   \date 1997-2002
*/

#define HAMOLECULE_CPP

#include <mpi.h>

#include <assert.h>
#include <float.h>
#include <math.h>

#include <boost/algorithm/string.hpp>

#include "hamolview.h" 

#include "hamolecule.h"
#include "haatom.h"
#include "habond.h"
#include "command.h"
#include "haqchem.h"
#include "tokens.h"
#include "abstree.h"
#include "haresdb.h"
#include "gaufile.h"

int HaMolecule::SeqFormat=0;

HaMolecule::HaMolecule(MolSet* new_phost_mset, std::string new_name):
Object3D(OBJ3D_MOLECULE)
{
	phost_mset = new_phost_mset;
	name=new_name;
	if( phost_mset != NULL ) phost_mset->AddObject3D(this);

	sec_struct_found= false;
	serno = 0;
	charge = -1000;
	
	structsource=SourceNone;

}

HaMolecule::HaMolecule(const HaMolecule &Mol_ref):
Object3D(Mol_ref)
{
	AddMolCopy((HaMolecule&)Mol_ref);
}

std::string HaMolecule::GetName() const 
{ 
	return mol_name; 
}

int HaMolecule::GetSerNo() const
{
	return serno;
}

std::string HaMolecule::GetRef() const
{
	std::string mol_ref = mol_name;
	const MolSet* pmset = this->GetHostMolSet();
	if (pmset->name_mol_map.count(mol_name) > 2)
	{
		mol_ref += "[" + std::to_string(serno) + "]";
	}
	return mol_ref;
}


bool HaMolecule::AddMolCopy(HaMolecule& Mol_ref, bool create_new_chain, AtomAtomMap* ptr_atom_map)
{
	list<HaChain>::iterator chain_itr;
	std::map<int,HaResidue*, less<int>  >::iterator res_itr;

	HaAtom* faptr;
	MolSet* pmset = GetHostMolSet();
	vector<HaAtom*>::iterator paitr;

	AtomAtomMap at_map;

	HaChain* pch_cur = NULL;
	HaResidue* fres = NULL;

	for( chain_itr= Mol_ref.Chains.begin();chain_itr != Mol_ref.Chains.end(); chain_itr++)
	{
		int ch_new;
		if( pch_cur == NULL || create_new_chain)
		{
			ch_new = 1;
		}
		else
		{
			ch_new = 0;
		}

		ResidueIteratorChain ritr_ch( &(*chain_itr));
		HaResidue* pres;
		for( pres = ritr_ch.GetFirstRes(); pres; pres = ritr_ch.GetNextRes() )
		{
			int res_new=1;
			for(paitr= pres->begin(); paitr != pres->end(); paitr++)
			{
				if(ch_new)
				{
					pch_cur = AddChain((*chain_itr).ident);
					pch_cur->SetParamFrom(*chain_itr);
					ch_new=0;
				}
				if(res_new)
				{
					int new_resno = pres->GetSerNo();
					fres= pch_cur->AddResidue(new_resno);
					if( fres == NULL)
					{
						new_resno = pch_cur->GetUniqResSerNo(false);
						fres = pch_cur->AddResidue(new_resno);
					}
					fres->SetParamFrom(*pres);
					fres->serno = new_resno;
					res_new=0;
				}
				faptr=fres->AddNewAtom();
				if( ptr_atom_map != NULL)
				{
					(*ptr_atom_map)[*paitr] = faptr;
				}
				faptr->SetParamFrom(*(*paitr));
				at_map[*paitr] = faptr;
			}
		}
	}
	
	std::vector<HaBond*>::const_iterator bitr;
	
	const MolSet* pmset_ref = Mol_ref.GetHostMolSet();

	for(bitr = pmset_ref->Bonds.begin(); bitr != pmset_ref->Bonds.end(); bitr++ )
	{
		const HaBond* bptr_ref = *bitr;
		if(bptr_ref->srcatom->GetHostMol() != &Mol_ref || bptr_ref->dstatom->GetHostMol() != &Mol_ref ) 
		{
			continue;
		}
		map<HaAtom*,HaAtom*, less<HaAtom*> >::iterator mitr;
		mitr =at_map.find(bptr_ref->srcatom);
		if(mitr == at_map.end() )
			continue;
		HaAtom* faptr1 = (*mitr).second;
		mitr =at_map.find(bptr_ref->dstatom);
		if(mitr == at_map.end() )
			continue;
		HaAtom* faptr2 = (*mitr).second;
				
		HaBond* fbptr= pmset->AddBond(faptr1,faptr2);
		fbptr->SetParamFrom(*bptr_ref);	
	}
	
//	list<HaHBond>::iterator hbitr;
//	for(hbptr= Mol_ref.HBonds.begin(); hbitr != Mol_ref.HBonds.end(); hbitr++)
//	{		
//			HaAtom* faptr1=at_map[(*hbitr).srcatom];
//			HaAtom* faptr2=at_map[(*hbitr),dstatom];
//			fbptr= AddHBond(faptr1,faptr2,bptr->flag);
//			fbptr->SetParamFrom(*bptr);	
//	}

	return true;
}

int HaMolecule::AttachFragment(HaAtom* catom_host, HaAtom* catom_frag )
{
	if( catom_host == NULL || catom_frag == NULL)
		return FALSE;

    HaMolecule* frag_mol = catom_frag->GetHostMol();
    HaMolecule* host_mol = catom_host->GetHostMol();

	PrintLog(" HaMolecule::AttachFragment() pt 1 \n");

	if( frag_mol == host_mol )
	{
		ErrorInMod(" HaMolecule::AttachFragment() ",
			       " Host and fragment molecules are the same, fragment can't be attached \n");
		return False;
	}

	MolSet* pmset = host_mol->GetHostMolSet();
	MolSet* pmset_frag = frag_mol->GetHostMolSet();

	bool create_new_chain = false;

	AtomAtomMap atom_map;

	PrintLog(" HaMolecule::AttachFragment() pt 2 \n");

	host_mol->AddMolCopy(*frag_mol, create_new_chain, &atom_map );

	PrintLog(" HaMolecule::AttachFragment() pt 3 \n");

	map<HaAtom*, HaAtom*, less<HaAtom*> >::iterator maitr;

	maitr = atom_map.find( catom_frag );
	if( maitr == atom_map.end() )
		return True;

	HaAtom* catom_frag_copy = (*maitr).second;

    AtomGroup bats_host,bats_frag; // the atoms connected to catom_host & catom_frag
	catom_host->GetBondedAtoms(bats_host);
	catom_frag_copy->GetBondedAtoms(bats_frag);

	Vec3D dir_host,dir_frag; // direction vectors to removed hydrogens from atoms they attached to
	                         // difene the direction the fragment will be attached
	Vec3D n1;
	dir_host[0] = 0.0; dir_host[1] = 0.0; dir_host[2] = 0.0;
	dir_frag[0] = 0.0; dir_frag[1] = 0.0; dir_frag[2] = 0.0;

	HaAtom* at2_host = NULL;
	HaAtom* at2_frag = NULL;

	if(!bats_host.empty())
	{
		AtomIteratorAtomGroup aitr_bats_host(&bats_host);
		if(catom_host->IsHydrogen())
		{
			at2_host = aitr_bats_host.GetFirstAtom(); 
			Vec3D::diff(dir_host,*catom_host,*at2_host);
			dir_host.normalize();
		}
		else 
		{ 
			for(at2_host = aitr_bats_host.GetFirstAtom();at2_host; at2_host = aitr_bats_host.GetNextAtom())
			{
				Vec3D::diff(n1,*catom_host,*at2_host);
				n1.normalize();
                Vec3D::sum(dir_host,dir_host,n1);
			}
			dir_host.normalize();
			at2_host = catom_host;
		}
	}


	if(!bats_frag.empty())
	{
		AtomIteratorAtomGroup aitr_bats_frag(&bats_frag);
		if(catom_frag_copy->IsHydrogen())
		{
			at2_frag = aitr_bats_frag.GetFirstAtom();  
			Vec3D::diff(dir_frag,*at2_frag,*catom_frag_copy);
			dir_frag.normalize();
		}
		else
		{
			for(at2_frag = aitr_bats_frag.GetFirstAtom();at2_frag; at2_frag = aitr_bats_frag.GetNextAtom())
			{
				Vec3D::diff(n1,*at2_frag,*catom_frag_copy);
				n1.normalize();
                Vec3D::sum(dir_frag,dir_frag,n1);
			}
			dir_frag.normalize();
			at2_frag = catom_frag_copy;
		}
	}

	double sina = -10.0;
	double cosa = -10.0;

	if( dir_host.length() > 0.1 && dir_frag.length() > 0.1 ) 
	{
		Vec3D::VecProduct(n1,dir_host,dir_frag);
		double sina = n1.length();
		double cosa = Vec3D::DotProduct(dir_host,dir_frag);

		if( sina < DBL_EPSILON)
		{
			Vec3D axx;
			axx[0] = 1.1; axx[1] = 1.2; axx[2] = 1.3;
			Vec3D::VecProduct(n1,axx,dir_host);
		}
	
		n1.normalize();

	// Rotate fragment so the forming bond will be opposite to other bonds in the molecule
		maitr = atom_map.begin();
		for(; maitr != atom_map.end(); maitr++)  
		{
			HaAtom* aptr = (*maitr).second;
			aptr->RotatePt(*catom_frag_copy,n1,cosa,-sina);
		}
	}

	double bdist = 1.5;

	if( at2_host != NULL && at2_frag != NULL && dir_host.length() > 0.1)
	{
		double bdist = HaAtom::StdBondLen(at2_host,at2_frag); 
		Vec3D disp;
		Vec3D::diff(disp,*at2_host,*at2_frag);
		disp = disp + bdist*dir_host;
		maitr = atom_map.begin();
		for(; maitr != atom_map.end(); maitr++)
		{
			HaAtom* aptr = (*maitr).second;
			Vec3D::sum(*aptr,*aptr,disp);
		}
	}

	HaAtom* atb1 = NULL;
	HaAtom* atb2 = NULL;

	if(catom_host->IsHydrogen() )
	{
		pmset->DeleteAtom(catom_host);
		if(at2_host != NULL) atb1 = at2_host;
	}
	else
	{
		atb1 = catom_host;
	}

	if(catom_frag_copy != NULL && catom_frag_copy->IsHydrogen() )
	{
		pmset->DeleteAtom(catom_frag_copy);
		if( at2_frag != NULL) atb2 = at2_frag;
	}
	else
	{
		atb2 = catom_frag_copy;
	}

	HaBond* fbptr;

	if( atb1 != NULL && atb2 != NULL)
	{
		fbptr = pmset->AddBond(atb1, atb2 );
		if( fbptr != NULL ) fbptr->DrawWire();
	}

	if(pmset == pmset_frag)
	{
		pmset->DeleteMol(frag_mol);
	}

	return True; 
}		

int HaMolecule::CombineMolecules(HaMolecule* frag_mol, HaAtom* catom_host, HaAtom* catom_frag )
{
	if( frag_mol == NULL )
		return False;

	MolSet* host_mset = GetHostMolSet();
	MolSet* frag_mset = frag_mol->GetHostMolSet();

	if( host_mset == frag_mset )
	{
		PrintLog(" Molecules %s and %s belong to different Molecular Sets and can't be combined \n",
			       GetObjName(), frag_mol->GetObjName());

	}

	int result = HaMolecule::AttachFragment( catom_host, catom_frag );

	return True;
}


HaMolecule::~HaMolecule(void)
{
	
}

int HaMolecule::SetAtomScreenCoord(HaMolView* pview, int x_sh, int y_sh)
{
	if( pview == NULL )
	{
		return False;
	}

	int ixadd=  pview->pCanv->XRange()/2 + x_sh;
	int iyadd=  pview->pCanv->YRange()/2 + y_sh;
	int izadd=  pview->ZOffset();

	HaAtom* aptr;

//	PrintLog("%8.3f %8.3f %8.3f \n", pview->Orig(1), pview->Orig(2),pview->Orig(3));

	AtomIteratorMolecule aitr(this);

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{	
		double x_tr, y_tr, z_tr;
		
		pview->GetTransfCoord( aptr->GetX(),aptr->GetY(), aptr->GetZ(),x_tr,y_tr,z_tr);
		
		aptr->x = (int)(x_tr*pview->Scale) + ixadd;
		aptr->y = (int)(y_tr*pview->Scale) + iyadd;
		aptr->z = (int)(z_tr*pview->Scale) + izadd;
	}
	return True;
}

int HaMolecule::RotateObjFromWorld( const HaMat_double& rot_mat, const Vec3D& cnt)  
{
	// this matrix rotate the body with respect to the world frame
	//Object3D::RotateObj( rot_mat, cnt);
	
	HaAtom* aptr;
	AtomIteratorMolecule aitr(this);

	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		double x = aptr->GetX();
		double y = aptr->GetY(); 
		double z = aptr->GetZ(); 
		//PrintLog("Xini= %2.3f ", aptr->GetY());
		aptr->SetX(rot_mat(1,1) * x + rot_mat(1,2) * y + rot_mat(1,3) * z);
		aptr->SetY(rot_mat(2,1) * x + rot_mat(2,2) * y + rot_mat(2,3) * z);
		aptr->SetZ(rot_mat(3,1) * x + rot_mat(3,2) * y + rot_mat(3,3) * z);
		//PrintLog("Xafter = %2.3f\n", aptr->GetY());
	}
	return TRUE;
}

int HaMolecule::RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt)  
{
	// this matrix rotate the body with respect to its old attitude to the new attitude (rotation matrix with respect to
	//its own body frame)
	Object3D::RotateObj( rot_mat, cnt);
	
	HaAtom* aptr;
	AtomIteratorMolecule aitr(this);

	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		double x = aptr->GetX() - cnt.GetX(); 
		double y = aptr->GetY() - cnt.GetY();
		double z = aptr->GetZ() - cnt.GetZ();
		//PrintLog("Xini= %2.3f ", aptr->GetY());
		aptr->SetX(cnt.GetX() + rot_mat(1,1) * x + rot_mat(1,2) * y + rot_mat(1,3) * z);
		aptr->SetY(cnt.GetY() + rot_mat(2,1) * x + rot_mat(2,2) * y + rot_mat(2,3) * z);
		aptr->SetZ(cnt.GetZ() + rot_mat(3,1) * x + rot_mat(3,2) * y + rot_mat(3,3) * z);
		//PrintLog("Xafter = %2.3f\n", aptr->GetY());
	}
	return TRUE;
}

int HaMolecule::Translate(const Vec3D& tr_vec )
{	
	bool add_x = ( fabs(tr_vec[0]) > DBL_EPSILON );
	bool add_y = ( fabs(tr_vec[1]) > DBL_EPSILON );
	bool add_z = ( fabs(tr_vec[2]) > DBL_EPSILON );
	
	AtomIteratorMolecule aitr(this);
	HaAtom* aptr;

	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		if( add_x ) aptr->SetX(aptr->GetX() + tr_vec[0]);
		if( add_y ) aptr->SetY(aptr->GetY() + tr_vec[1]);
		if( add_z ) aptr->SetZ(aptr->GetZ() + tr_vec[2]);
	}
	
	return True;
}


int AtomContainer::SetIntCoordFromStr(const char* int_coord_str)
{
	double phi,cost,psi;
	Vec3D trans;
	int num;

	num = sscanf(int_coord_str,"%lf %lf %lf %lf %lf %lf",
		         &trans[0],&trans[1],&trans[2],&phi,&cost,&psi);

	if( num == 6)
	{
		PrintLog(" Setting int coord for mol1 x=%9.4f y=%9.4f z=%9.4f \n phi=%9.4f cos_theta=%9.4f psi=%9.4f \n\n",
			     trans[0],trans[1],trans[2],phi,cost,psi);
		SetPosEulerTrans(phi,cost,psi,trans);
	}
	else
	{
		PrintLog(" Number of internal coordinates for the molecule is less than 6 \n");
		return FALSE;
	}
	return TRUE;
}

int AtomContainer::SaveXYZFile( const char* fout_name, const AtomSaveOptions* p_opt )
{
	ofstream sout(fout_name);
	if(!sout.good())
	{
		PrintLog(" Error in AtomContainer::SaveXYZFile() \n");
		PrintLog(" Error to Open File %s \n",fout_name);
		return FALSE;
	}
	int ires = SaveXYZStream(sout, p_opt );
	return ires;
}


int AtomContainer::SaveXYZStream(std::ostream& sout, const AtomSaveOptions* p_opt_arg )
{
	auto_ptr<AtomSaveOptions> p_opt( p_opt_arg != NULL ? (AtomSaveOptions*) p_opt_arg->clone() : new AtomSaveOptions() ); 

	int na = this->GetNAtoms();
	
	if(!sout.good() )
	{
		PrintLog(" Error in AtomContainer::SaveXYZStream() \n");
		PrintLog(" Can not write to output stream \n");
		return FALSE;
	}
	
	int set_atm_idx_maps = FALSE;
	
	std::auto_ptr<AtomIterator> p_aitr( this->GetAtomIteratorPtr() );
	HaAtom* aptr;

	AtomGroup selected_atoms;
	if(p_opt->save_selected )
	{
		for(aptr = p_aitr->GetFirstAtom(); aptr; aptr = p_aitr->GetNextAtom() )
		{
			selected_atoms.push_back(aptr);
		}
		p_aitr.reset( selected_atoms.GetAtomIteratorPtr() );
		na = selected_atoms.size();
	}
	
	AtomIntMap atm_idx_map;

	if( p_opt->save_connect ) set_atm_idx_maps = TRUE;

	int idx;
	if( set_atm_idx_maps )
	{
		idx = 1;
		for( aptr = p_aitr->GetFirstAtom(); aptr; aptr = p_aitr->GetNextAtom() )
		{
			atm_idx_map[aptr] = idx;
			idx++;
		}
	}
	std::vector< std::vector<int> > bonded_atm_idx;

	if( p_opt->save_connect )
	{
		bonded_atm_idx.resize(na);
		idx = 1;
		AtomGroup bonded_atoms;
		for( aptr = p_aitr->GetFirstAtom(); aptr; aptr = p_aitr->GetNextAtom() )
		{
			aptr->GetBondedAtoms(bonded_atoms);
			AtomIteratorAtomGroup aitr_b( &bonded_atoms);
			HaAtom* aptr_b;
			for( aptr_b = aitr_b.GetFirstAtom(); aptr_b; aptr_b = aitr_b.GetNextAtom())
			{
				AtomIntMap::iterator idx_b_itr = atm_idx_map.find(aptr_b);
				if( idx_b_itr != atm_idx_map.end() )
				{
					int idx_b = (*idx_b_itr).second;
					bonded_atm_idx[idx-1].push_back(idx_b);
				}
			}
			idx++;
		}
	}

	HaMolView* pview = NULL;
	MolSet*  pmset = NULL;
	if( na > 0)
	{
		aptr = p_aitr->GetFirstAtom();
		pmset = aptr->GetHostMolSet();
		pview = pmset->GetActiveMolView();
	}
	
	char buf[256];

	sprintf(buf,"%9d",na);
	sout << buf << endl;
	idx = 1;
	for( aptr = p_aitr->GetFirstAtom(); aptr; aptr = p_aitr->GetNextAtom() )
	{
		int elem_no = aptr->GetElemNo();
	
		double x = aptr->GetX_Ang();
		double y = aptr->GetY_Ang();
		double z = aptr->GetZ_Ang();

		if( p_opt->save_transform && pview != NULL )
		{	
			pview->GetTransfCoord(aptr->GetX_Ang(), aptr->GetY_Ang(), aptr->GetZ_Ang(),x,y,z);
		}

		sprintf(buf," %7d  %4d %16.9f %16.9f %16.9f ",idx,elem_no,x,y,z );
		sout << buf ;
		
		int idx_ff = 0;
		sprintf(buf," %5d ",idx_ff);
		sout << buf;

		int i;
		if( p_opt->save_connect )
		{
			int nb = bonded_atm_idx[idx-1].size();
			for(i = 0; i < nb; i++)
			{
				sprintf(buf," %5d",bonded_atm_idx[idx-1][i]);
				sout << buf;
			}	
		}
		if( p_opt->save_atom_ref )
		{
			std::string at_lbl = aptr->GetRef( p_opt->at_ref_type );
			sprintf(buf,"   # %s",at_lbl.c_str());
			sout << buf;
		}
		sout << endl;
		idx++;
	}
	if( !sout.good() ) return FALSE;

	return TRUE;
}

int AtomContainer::SaveGROFile(const char* fout_name, const AtomSaveOptions* p_opt)
{
	ofstream sout(fout_name);
	if (!sout.good())
	{
		PrintLog(" Error in AtomContainer::SaveGROFile() \n");
		PrintLog(" Error to Open File %s \n", fout_name);
		return FALSE;
	}
	int ires = SaveGROStream(sout, p_opt);
	return ires;
}

int AtomContainer::SaveGROStream(std::ostream& sout, const AtomSaveOptions* p_opt_arg)
{
	auto_ptr<AtomSaveOptions> p_opt(p_opt_arg != NULL ? (AtomSaveOptions*)p_opt_arg->clone() : new AtomSaveOptions());

	int na = this->GetNAtoms();

	if (!sout.good())
	{
		PrintLog(" Error in AtomContainer::SaveGROStream() \n");
		PrintLog(" Can not write to output stream \n");
		return FALSE;
	}

	int set_atm_idx_maps = FALSE;

	std::auto_ptr<AtomIterator> p_aitr(this->GetAtomIteratorPtr());
	HaAtom* aptr;

	AtomGroup selected_atoms;
	if (p_opt->save_selected)
	{
		for (aptr = p_aitr->GetFirstAtom(); aptr; aptr = p_aitr->GetNextAtom())
		{
			selected_atoms.push_back(aptr);
		}
		p_aitr.reset(selected_atoms.GetAtomIteratorPtr());
		na = selected_atoms.size();
	}

	AtomIntMap atm_idx_map;

	if (p_opt->save_connect) set_atm_idx_maps = TRUE;

	int idx;
	if (set_atm_idx_maps)
	{
		idx = 1;
		for (aptr = p_aitr->GetFirstAtom(); aptr; aptr = p_aitr->GetNextAtom())
		{
			atm_idx_map[aptr] = idx;
			idx++;
		}
	}
	std::vector< std::vector<int> > bonded_atm_idx;

	//if (p_opt->save_connect)
	//{
	//	bonded_atm_idx.resize(na);
	//	idx = 1;
	//	AtomGroup bonded_atoms;
	//	for (aptr = p_aitr->GetFirstAtom(); aptr; aptr = p_aitr->GetNextAtom())
	//	{
	//		aptr->GetBondedAtoms(bonded_atoms);
	//		AtomIteratorAtomGroup aitr_b(&bonded_atoms);
	//		HaAtom* aptr_b;
	//		for (aptr_b = aitr_b.GetFirstAtom(); aptr_b; aptr_b = aitr_b.GetNextAtom())
	//		{
	//			AtomIntMap::iterator idx_b_itr = atm_idx_map.find(aptr_b);
	//			if (idx_b_itr != atm_idx_map.end())
	//			{
	//				int idx_b = (*idx_b_itr).second;
	//				bonded_atm_idx[idx - 1].push_back(idx_b);
	//			}
	//		}
	//		idx++;
	//	}
	//}

	HaMolView* pview = NULL;
	MolSet* pmset = NULL;
	if (na > 0)
	{
		aptr = p_aitr->GetFirstAtom();
		pmset = aptr->GetHostMolSet();
		pview = pmset->GetActiveMolView();
	}

	sout << " GROMACS Gro file Generated by HARLEM " << std::endl;
	char buf[256];
	sprintf(buf, "%5d", na);
	sout << buf << endl;
	idx = 1;
	for (aptr = p_aitr->GetFirstAtom(); aptr; aptr = p_aitr->GetNextAtom())
	{
		int elem_no = aptr->GetElemNo();
		HaResidue* pres = aptr->GetHostRes();
		int resno = pres->GetSerNo();
		std::string res_name = pres->GetName();

		double x = aptr->GetX_Ang()*0.1;
		double y = aptr->GetY_Ang()*0.1;
		double z = aptr->GetZ_Ang()*0.1;

		if (p_opt->save_transform && pview != NULL)
		{
			pview->GetTransfCoord(aptr->GetX_Ang(), aptr->GetY_Ang(), aptr->GetZ_Ang(), x, y, z);
		}

		//if( res_name.size() < 4 ) 

		sprintf(buf, "%5d%-4.4s  %4.4s%5d%8.3f%8.3f%8.3f ", resno, res_name.c_str(), aptr->GetName(), idx, x, y, z);
		sout << buf << std::endl;

		idx++;
	}
	double box[3] = { 0.0,0.0,0.0 };
	if (pmset != NULL && pmset->per_bc->IsSet())
	{
		box[0] = pmset->per_bc->GetA() * 0.1;
		box[1] = pmset->per_bc->GetB() * 0.1;
		box[2] = pmset->per_bc->GetC() * 0.1;
	}
	sprintf(buf, "%10.5f%10.5f%10.5f", box[0], box[1], box[2]);
	sout << buf << std::endl;

	if (!sout.good()) return FALSE;

	return TRUE;
}



// function to get from the name of the style( prefix_num) get prefix and num
static int GetUSPrefix(const char* name, std::string& prefix, int& num)
{
	int len = strlen(name);
	prefix =  name;
	num = 0;
	char buf[10];
	buf[0] = 0;
	buf[1] = 0;

	std::string num_str;
	if( !isdigit(name[len-1]) ) return TRUE;
	int i;
	for( i = len ; i > 0; i--)
	{
		buf[0] = name[i-1];
		if( isdigit(buf[0]))
		{
			num_str.insert(0,buf);
			prefix.erase(i-1,1);
			continue;
		}
		if( buf[0] == '_')
		{
			prefix.erase(i-1,1);
		}
		break;
	}
	num = atoi(num_str.c_str());
	return TRUE;
}

bool HaMolecule::SetObjName(const char* new_name)
{
	name = "";
	int len = strlen(new_name);
	for(int i=0; i < len; i++)
	{
		if( !isspace(new_name[i]) ) 
		name+= toupper(new_name[i]);
	}
	MolSet* pmset = GetHostMolSet();
	
// Set Unique name
	int nmol = pmset->GetNMol();
	int imol;
	int unique = TRUE;
	for( imol = 0 ; imol < nmol; imol++)
	{
		HaMolecule* pMol = pmset->HostMolecules[imol];
		if(pMol == this) continue;
		if(!strcmp(name.c_str(),pMol->GetObjName()))
		{
			 unique = FALSE;
			 break;
		}
	}

	if( unique == FALSE )
	{
		std::string prefix1,prefix2;
		int num1, num2;
		GetUSPrefix(name.c_str(), prefix1, num1);

		int maxn = MaxFun(0,num1);

		for( imol = 0 ; imol < nmol; imol++)
		{
			HaMolecule* pMol = pmset->HostMolecules[imol];
			if(pMol == this) continue;

			GetUSPrefix(pMol->GetObjName(), prefix2, num2);

			if(prefix2 == prefix1)
			{
				maxn = MaxFun(maxn,1);
				maxn = MaxFun(maxn,num2);
			}
		}
		char buf[20];
		maxn++;
		if(maxn < 10)
			sprintf(buf,"%1d",maxn);
		else if(maxn < 100 && maxn > 9)
			sprintf(buf,"%2d",maxn);
		else if(maxn < 1000 && maxn > 99)
			sprintf(buf,"%3d",maxn);
		else if(maxn < 10000 && maxn > 999)
			sprintf(buf,"%4d",maxn);
		else if(maxn < 100000 && maxn > 9999)
			sprintf(buf,"%5d",maxn);

		name = prefix1 + "_" + buf;
	}
	
	return true;
}

bool HaMolecule::FillRef(char* buf,int mode) const
{
	std::string mol_ref = this->GetRef();
	sprintf(buf,"$%s$*", mol_ref.c_str() );
	return TRUE;
}

bool HaMolecule::InitAtoms(GauFile& gfile)
// set Atom Info from Gaussian File Mol structure
{
	char at_name[10]; // Atom Symbol
	int iel;
	mol_type gmol;

	int iunit= -GauFile::IO_mol;
	int len = (int)sizeof(gmol) / ((int)sizeof(double));
	int offset = 0;
	void* pdata = &gmol;
	gfile.fileio("read",iunit,len,pdata,offset);

	int old_serno;

	HaResidue* pres    = AddChainAndResidue();
	HaChain*   pch_cur = pres->GetHostChain();

	for(int i=0; i < gmol.natoms; i++)
	{
		iel= gmol.ian[i];
		
		old_serno=pres->serno;
		
		HaAtom* aptr=pres->AddNewAtom();
		aptr->SetElemNo( gmol.ian[i]);
		std::string elem_symbol = HaAtom::GetStdSymbolElem(aptr->GetElemNo());

		if(pres->size() > 50)
		{
			pres = pch_cur->AddResidue(old_serno+1);
			pres->SetName( "MOL");
		}

		int iat = pres->size() + 1;

		if(iat < 10)
			sprintf(at_name,"%s%1d",elem_symbol.c_str(),iat);
		else
			sprintf(at_name,"%s%2d",elem_symbol.c_str(),iat);

		aptr->SetName(at_name);
		aptr->SetElemNo(HaAtom::GetElemNoFromName(at_name, aptr->GetHostRes()) );
		aptr->SetX(gmol.c[i][0]);
		aptr->SetY(gmol.c[i][1]);
		aptr->SetZ(gmol.c[i][2]);
	}	
	return true;
}

HaAtom* HaMolecule::GetAtomByRef(const char* at_ref)
{
	try
	{
		std::string at_ref_str(at_ref);
		int iat_name = at_ref_str.find_first_of('.');
		if(iat_name == -1) throw std::runtime_error(" Atom Reference " + at_ref_str + " doesn't have . separator for atom name ");
	
		std::string AtName= at_ref_str.substr(iat_name+1);
		boost::trim(AtName);
		boost::to_upper(AtName);
		if( AtName.empty() ) throw std::runtime_error(" No Atom name is specified in Atom Reference " + at_ref_str );
	
		HaResidue* pres = GetResByRef(at_ref_str.substr(0,iat_name));
		if(pres == NULL) throw std::runtime_error(" Can not find residue for " + at_ref_str );
	
		return(pres->GetAtomByName(AtName.c_str()));
	}
	catch( const std::exception& ex )
	{
		bool debug = false; 
		if( debug ) PrintLog(" Error in HaMolecule::GetAtomByRef()\n%s\n",ex.what());	
	}
	return NULL;	
}

HaResidue* HaMolecule::GetResByRef(const std::string& res_str_par )
{
	std::string res_str = res_str_par;
	boost::trim( res_str );
	try
	{
		int ich_id= res_str.find_first_of(':');
		int ires_end,i;
		char chain_id;
		if( ich_id == -1)
		{
			chain_id = ' ';
			ires_end = res_str.size()-1;
		}
		else
		{
			ires_end= ich_id -1;
			ich_id++;
			if(ich_id == res_str.size()) // no space between ':' and '.'
			{
				chain_id= ' ';
			}
			else
			{
				chain_id= res_str[ich_id];
			}
		}

		HaChain* chn_ptr= GetChain(chain_id);
		if(chn_ptr == NULL) 
		{
			std::string chain_str;
			chain_str += chain_id;
			throw std::runtime_error( " No Chain " + chain_str  );
		}
	
		int ibeg_res_num = ires_end + 1;
		for( i = ires_end; i >= 0; i-- )
		{
			if(!isdigit(res_str[i])) break;
		}
		ibeg_res_num = i+1;
	
		std::string res_num_str = res_str.substr(ibeg_res_num, ires_end - ibeg_res_num + 1);
		if(res_num_str.size() == 0) throw std::runtime_error( " No Residue Number in Reference " + res_str );
		
		int res_ser_no = 0;
		if( harlem::IsInt(res_num_str) ) res_ser_no = atoi( res_num_str.c_str());
	
		HaResidue* pres= chn_ptr->GetResBySerNo( res_ser_no );
		if(pres == NULL) throw std::runtime_error( " No Residue with number " + harlem::ToString( res_ser_no ) + " in chain " + chain_id );
		if(ibeg_res_num == 0) return pres;

		std::string res_name= res_str.substr(0,ibeg_res_num);

		std::string res_name_found =  pres->GetName();
		boost::trim(res_name_found); 

		if(stricmp_loc(res_name.c_str(), res_name_found.c_str()))
		{
			throw std::runtime_error( " Residue number " + harlem::ToString(res_ser_no ) +  " has a different name from that in " +  res_str );
		}
		return pres;
	}
	catch( const std::exception& ex )
	{
		bool debug = false; 
		if( debug ) PrintLog(" Error in HaMolecule::GetResByRef()\n%s\n",ex.what());	
	}
	return(NULL);

}

HaChain* HaMolecule::GetChain(const char chain_id)
{
	list<HaChain>::iterator ch_itr;
	for(ch_itr = Chains.begin(); ch_itr != Chains.end(); ch_itr++)
	{
		if((*ch_itr).ident == chain_id)
			return &(*ch_itr);
	}
	return NULL;
}


int HaMolecule::GetNAtoms() const
{
	int na = 0;
	AtomIteratorMolecule_const aitr(this);
	const HaAtom* aptr;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		na++;
	}
    return(na);
}


bool
HaMolecule::Print_info(ostream &sout, const int level) const
// Output information about the molecule
{
	sout << "This is the molecule info " << endl;
	sout << "Atoms info: " << endl;

	AtomIteratorMolecule_const aitr(this);
	const HaAtom* aptr;

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		aptr->Print_info(sout,level);
		sout << endl;
	}


	return true;
}

HaBond* MolSet::AddBond(HaAtom* src, HaAtom* dst )
{
	if(src == NULL || dst == NULL)
	{
		PrintLog(" Error in MolSet::AddBond() \n"); 
		PrintLog(" One of the Atom Pointers is NULL \n"); 
		return NULL;
	}

	if(src == dst)
	{
		PrintLog(" Error in MolSet::AddBond() \n"); 
		PrintLog(" Two atom pointers are identical \n"); 
		return NULL;
	}
	
	HaBond* pb = new HaBond(src,dst);
	Bonds.push_back(pb);
	src->p_bonds->push_back(pb); 
	dst->p_bonds->push_back(pb);
	return( pb );
}


bool MolSet::DeleteBond(HaAtom* src, HaAtom* dst)
{
	if(src == NULL || dst == NULL) return false;

	std::vector<HaBond*>& bonds1 = src->GetBonds();
	std::vector<HaBond*>& bonds2 = dst->GetBonds();

	HaBond* bnd_d = NULL;

	std::vector<HaBond*>::iterator bitr = bonds1.begin();
	while( bitr != bonds1.end() )
	{
		HaAtom* aptr1 = (*bitr)->GetFirstAtom();
		HaAtom* aptr2 = (*bitr)->GetSecondAtom();
		
		if( aptr1 != dst && aptr2 != dst ) 
		{
			++bitr;
			continue;
		}		
		bnd_d = *bitr;
		bonds1.erase(bitr);
		break;
	}

	if( bnd_d == NULL )
	{
		PrintLog("Error in MolSet::DeleteBond() \n");
		PrintLog(" Atoms %s  and %s are not bonded \n", src->GetRef().c_str(), dst->GetRef().c_str());
		return false;
	}
	
	for( bitr = bonds2.begin(); bitr != bonds2.end(); ++bitr )
	{
		if( *bitr != bnd_d ) continue; 
		bonds2.erase(bitr);
		break;
	}

	for( bitr = Bonds.begin(); bitr != Bonds.end(); ++bitr )
	{
		if( *bitr != bnd_d ) continue; 
		Bonds.erase(bitr);
		delete bnd_d;
		break;
	}
	return true;
}	

int HaMolecule::GetNRes() const
{
	list<HaChain>::const_iterator citr;
	int icount=0;
	for(citr=Chains.begin(); citr != Chains.end(); citr++)
		icount+= (*citr).GetNRes();
	return icount;
}

HaChain* HaMolecule::AddChain(char ident)
{
	   Chains.push_back(HaChain(this,ident));
	   HaChain* pch_cur = &(Chains.back());

	   return( pch_cur );
}

char HaMolecule::GetChainIdentMax()
{
	ChainIteratorMolecule chitr(this);
	char ident_max = ' ';
	HaChain* chain;
	for( chain = chitr.GetFirstChain(); chain; chain = chitr.GetNextChain() )
	{
		if( chain->ident > ident_max ) ident_max = chain->ident;
	}
	return ident_max;
}

HaChain* HaMolecule::GetFirstChain()
{
	if( Chains.empty() ) return NULL;
	return &(*Chains.begin());
}

SecStructElement* HaMolecule::AddFeature()
{
	Features.push_back(SecStructElement());
	return(&Features.back());
}

void
HaMolecule::DeleteFeatures(const int itype)
{
	list<SecStructElement>::iterator fitr;
	for(fitr=Features.begin(); fitr != Features.end();)
	{
		if( (*fitr).type == itype ) 
		{
			fitr = Features.erase(fitr);
		}
		else
		{
			fitr++;
		}
	}
}

int HaMolecule::GetNumFeatures(const int itype) const
{
	list<SecStructElement>::const_iterator fitr;
	int icount=0;
	for(fitr=Features.begin(); fitr != Features.end(); fitr++)
	{
		if( (*fitr).type == itype)
			icount++;
	}
	return icount;
}





bool HaMolecule::SetUniqueAtomNames()
{
	HaChain* chain;
	ChainIteratorMolecule ch_itr(this);
	for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
	{
		chain->SetUniqueResNo();
	}

	list<HaResidue*> old_groups;
	list<HaResidue*>::iterator ritr;

	ResidueIteratorMolecule res_itr(this);

	HaResidue* pres;
	for(pres = res_itr.GetFirstRes(); pres; pres = res_itr.GetNextRes())
	{
	   old_groups.push_back(pres);
	}

	for(ritr = old_groups.begin(); ritr != old_groups.end(); ritr++)
	{
	   (*ritr)->SplitResidue();
	}

	for(pres = res_itr.GetFirstRes(); pres; pres = res_itr.GetNextRes())
	{
		pres->SetUniqueAtomNames();
	}
    return true;
}



static char ResidueChar[29] = {
        'A', 'G', 'L', 'S', 'V', 'T', 'K', 'D', 'I', 'N',
        'E', 'P', 'R', 'F', 'Q', 'Y', 'H', 'C', 'M', 'W',
        'B', 'Z', '*', 'P',
        'A', 'C', 'G', 'T',
        'U'
    };


void HaMolecule::DescribeSequence()
{
    HaChain  *chain;
    HaResidue  *grp;
    int ichn,count;
    const char *str;

	PrintLog("\n");
   
	ChainIteratorMolecule chitr(this);
    for( chain = chitr.GetFirstChain(); chain; chain= chitr.GetNextChain())
    {
		ichn = (GetNChains() < 2);  count = 0;
		ResidueIteratorChain ritr_ch(chain);
        for( grp = ritr_ch.GetFirstRes(); grp; grp = ritr_ch.GetNextRes() )
            if( (!grp->empty()) && !((grp->front())->flag&HeteroFlag) )
            {
				if( !ichn )
                {
				    PrintLog("Chain %c:\n",chain->ident);
                    ichn = True;
                }

                if( !SeqFormat )
                {
					if( count == 10 )
                    {
						PrintLog("\n");
                        count = 1;
                    }
					else
						count++;

                    str = grp->GetName();

                    PrintLog("%c%c%c%-3d ",str[0],str[1],str[2],grp->serno);
                }
				else
                {
					if( count == 60 )
                    {
						PrintLog("\n");
                        count = 1;
                    }
					else
						count++;

                    if( grp->refno < 29 )
                    {
						PrintLog("%c",ResidueChar[grp->refno]);
                    }
					else
						PrintLog("%c",'*');
                }
            }
        PrintLog("\n");
    }
    PrintLog("\n");
}


SecStructElement::SecStructElement()
{

}

SecStructElement::~SecStructElement()
{

}

bool
SecStructElement::operator==(const SecStructElement& rhs)
{
  if(init != rhs.init) return false;
  if(term != rhs.term) return false;
  if(chain != rhs.chain) return false;
  if(type != rhs.type) return false;

  return true;
}

bool
SecStructElement::operator < (const SecStructElement& rhs)
{
  if(*this == rhs) return false;
  if(type > rhs.type) return false;
  if(chain > rhs.chain) return false;
  if(init > rhs.init) return false;
  if(term > rhs.term) return false;
  return true;
}

AtomIteratorMolecule::AtomIteratorMolecule(HaMolecule* new_pMol)
{
    pMol = new_pMol;
	if(pMol == NULL) throw "Error in: AtomIteratorMolecule::AtomIteratorMolecule() pMol == NULL \n";

	GetFirstAtom();
}

AtomIteratorMolecule::~AtomIteratorMolecule()
{
	
}


HaAtom* 
AtomIteratorMolecule::GetFirstAtom()
{
	if(pMol == NULL) { return NULL; }

	HaChain* chain;
	HaAtom* aptr = NULL;
	for( ch_itr = pMol->Chains.begin(); ch_itr != pMol->Chains.end(); ch_itr++)
	{
		chain = &(*ch_itr);
		aitr.SetForChain(chain); 
		aptr = aitr.GetFirstAtom();
		if( aptr != NULL) break;
	}
	return aptr;

}

HaAtom* 
AtomIteratorMolecule::GetNextAtom()
{
	if(pMol == NULL) { return NULL; }

	HaAtom* aptr = aitr.GetNextAtom();
	
	if( aptr != NULL) return aptr;

	ch_itr++;
	HaChain* chain;
	for(; ch_itr != pMol->Chains.end(); ch_itr++)
	{
		chain = &(*ch_itr);
		aitr. SetForChain(chain); 
		aptr = aitr.GetFirstAtom();
		if( aptr != NULL) break;
	}
	return aptr;
}

AtomIteratorMolecule_const::AtomIteratorMolecule_const(const HaMolecule* new_pMol)
{
    pMol = new_pMol;
	if(pMol == NULL) { return; }
	
	GetFirstAtom();
}

AtomIteratorMolecule_const::~AtomIteratorMolecule_const()
{
	
}

const HaAtom* AtomIteratorMolecule_const::GetFirstAtom()
{
	if(pMol == NULL) { return NULL; }

	const HaChain* chain;
	const HaAtom* aptr = NULL;
	for( ch_itr = pMol->Chains.begin(); ch_itr != pMol->Chains.end(); ch_itr++)
	{
		chain = &(*ch_itr);
		aitr.SetForChain(chain); 
		aptr = aitr.GetFirstAtom();
		if( aptr != NULL) break;
	}
	return aptr;

}

const HaAtom* AtomIteratorMolecule_const::GetNextAtom()
{
	if(pMol == NULL) { return NULL; }

	const HaAtom* aptr = aitr.GetNextAtom();
	
	if( aptr != NULL) return aptr;

	ch_itr++;
	const HaChain* chain;
	for(; ch_itr != pMol->Chains.end(); ch_itr++)
	{
		chain = &(*ch_itr);
		aitr.SetForChain(chain); 
		aptr = aitr.GetFirstAtom();
		if( aptr != NULL) break;
	}
	return aptr;
}

ResidueIteratorMolecule::ResidueIteratorMolecule(HaMolecule* new_pmol)
{
	pmol = new_pmol;
	if(pmol == NULL) { return; }
	
	ch_itr = pmol->Chains.begin(); 
	if( ch_itr != pmol->Chains.end() )
	{
		res_itr = (*ch_itr).res_arr.begin();
	}
}

ResidueIteratorMolecule::ResidueIteratorMolecule(const ResidueIteratorMolecule& ritr_ref)
{
	pmol    = ritr_ref.pmol;
	ch_itr  = ritr_ref.ch_itr;
	res_itr = ritr_ref.res_itr;
}


ResidueIteratorMolecule::~ResidueIteratorMolecule()
{

}

HaResidue* ResidueIteratorMolecule::GetFirstRes()
{
   if(pmol == NULL) return NULL; 
   if(pmol->Chains.empty()) return NULL;
   ch_itr = pmol->Chains.begin();
 
   while( ch_itr != pmol->Chains.end() )
   {
		if( (*ch_itr).res_map.empty() )
		{
			ch_itr++;
		}
		else
		{
            res_itr = (*ch_itr).res_arr.begin();
			return (*res_itr);
		}
   }

   return NULL;
}

HaResidue* ResidueIteratorMolecule::GetNextRes()
{
  if(pmol == NULL) { return NULL; }

  res_itr++;

  if( res_itr != (*ch_itr).res_arr.end())
  {
	   return (*res_itr);
  }

  ch_itr++;
  
  while(ch_itr != pmol->Chains.end())
  {
	  if( (*ch_itr).res_map.empty() ) 
	  { 
		  ch_itr++;
	  }
	  else
	  {
		  res_itr = (*ch_itr).res_arr.begin();
		  return (*res_itr);
	  }
  }
  return NULL;
}


ResidueIteratorMolecule_const::ResidueIteratorMolecule_const(const HaMolecule* new_pmol)
{
	pmol = new_pmol;
	if (pmol == NULL) { return; }

	ch_itr = pmol->Chains.begin();
	if (ch_itr != pmol->Chains.end())
	{
		res_itr = (*ch_itr).res_arr.begin();
	}
}

ResidueIteratorMolecule_const::ResidueIteratorMolecule_const(const ResidueIteratorMolecule_const& ritr_ref)
{
	pmol = ritr_ref.pmol;
	ch_itr = ritr_ref.ch_itr;
	res_itr = ritr_ref.res_itr;
}

ResidueIteratorMolecule_const::~ResidueIteratorMolecule_const()
{

}

const HaResidue* ResidueIteratorMolecule_const::GetFirstRes()
{
	if (pmol == NULL) return NULL;
	if (pmol->Chains.empty()) return NULL;
	ch_itr = pmol->Chains.begin();
	while (ch_itr != pmol->Chains.end())
	{
		if ((*ch_itr).res_map.empty())
		{
			ch_itr++;
		}
		else
		{
			res_itr = (*ch_itr).res_arr.begin();
			return (*res_itr);
		}
	}
	return NULL;
}

const HaResidue* ResidueIteratorMolecule_const::GetNextRes()
{
	if (pmol == NULL) { return NULL; }

	res_itr++;

	if (res_itr != (*ch_itr).res_arr.end())
	{
		return (*res_itr);
	}

	ch_itr++;

	while (ch_itr != pmol->Chains.end())
	{
		if ((*ch_itr).res_map.empty())
		{
			ch_itr++;
		}
		else
		{
			res_itr = (*ch_itr).res_arr.begin();
			return (*res_itr);
		}
	}
	return NULL;
}


ChainIteratorMolecule::ChainIteratorMolecule(HaMolecule* new_pmol)
{
	pmol = new_pmol;
	if(pmol == NULL) { return; }
	ch_itr = pmol->Chains.begin(); 
}

ChainIteratorMolecule::ChainIteratorMolecule(const ChainIteratorMolecule& chitr_ref)
{
	pmol    = chitr_ref.pmol;
	ch_itr  = chitr_ref.ch_itr;
}


ChainIteratorMolecule::~ChainIteratorMolecule()
{

}

HaChain* ChainIteratorMolecule::GetFirstChain()
{
   if(pmol == NULL) return NULL; 

   if(pmol->Chains.empty()) return NULL;
   
   ch_itr = pmol->Chains.begin();
   return &(*ch_itr);
}

HaChain* ChainIteratorMolecule::GetNextChain()
{
  if(pmol == NULL) { return NULL; }

  ch_itr++;

  if(ch_itr == pmol->Chains.end()) return NULL;
  return &(*ch_itr);
}

PointIterator* HaMolecule::GetPointIteratorPtr() 
{ 
	return new AtomIteratorMolecule(this); 
}
PointIterator_const* HaMolecule::GetPointIteratorPtr() const
{ 
	return new AtomIteratorMolecule_const(this); 
}


int HaMolecule::IsMember(const HaAtom* aptr) const
{
	if(aptr == NULL) return FALSE;
	if(aptr->GetHostMol() == this) return TRUE;
	return FALSE;
}
AtomIterator* HaMolecule::GetAtomIteratorPtr()
{
	return new AtomIteratorMolecule(this); 
}

int HaMolecule::GetNBonds()  const  
{
	const MolSet* pmset = GetHostMolSet();
	std::vector<HaBond*>::const_iterator bitr;
	int nb = 0;
	for( bitr = pmset->Bonds.begin(); bitr != pmset->Bonds.end(); bitr++ )
	{
		if( (*bitr)->GetFirstAtom()->GetHostMol() == this && (*bitr)->GetSecondAtom()->GetHostMol() == this )
		{
			nb++;
		}
	}
	return nb;
}

int HaMolecule::GetNHBonds() const
{
	const MolSet* pmset = GetHostMolSet();
	set<HaHBond, less<HaHBond> >::const_iterator bitr;
	int nb = 0;
	for( bitr = pmset->HBonds.begin(); bitr != pmset->HBonds.end(); bitr++ )
	{
		if( (*bitr).src->GetHostMol() == this && (*bitr).dst->GetHostMol() == this )
		{
			nb++;
		}
	}
	return nb;	
}  
	
int HaMolecule::GetNSSBonds() const  
{
	const MolSet* pmset = GetHostMolSet();
	std::vector<HaBond*>::const_iterator bitr;
	int nb = 0;
	for( bitr = pmset->Bonds.begin(); bitr != pmset->Bonds.end(); bitr++ )
	{
		int elem1 = (*bitr)->GetFirstAtom()->GetElemNo();
		int elem2 = (*bitr)->GetSecondAtom()->GetElemNo();
		
		if( elem1 != 16 || elem2 != 16) continue;

		if( (*bitr)->GetFirstAtom()->GetHostMol() == this && (*bitr)->GetSecondAtom()->GetHostMol() == this )
		{
			nb++;
		}
	}
	return nb;	
}

void HaMolecule::DescribeMolecule()
{
 
    PrintLog("\n");

    if( GetObjName() )
    {   
		PrintLog("Molecule name ....... %s\n",GetObjName());
    }

    if( !classification.empty() )
    {   
		PrintLog("Classification ...... %s\n",classification.c_str());
    }

    if( GetNRes() >1 )
    {   
		PrintLog("Secondary Structure . ");
        if( structsource == SourceNone )
        {   
			PrintLog("No Assignment\n");
        } 
		else if( structsource == SourcePDB )
        {   
			PrintLog("PDB Data Records\n");
        } 
		else 
			PrintLog("Calculated\n");
    }

    if( !identcode.empty() )
    {   
		PrintLog("Brookhaven Code ..... %s\n",identcode.c_str());
    }

    if( GetNChains() > 1 )
    {   
        PrintLog("Number of Chains .... %d\n",GetNChains());
    }

    PrintLog("Number of Residues .... %d\n", GetNRes());

    PrintLog("Number of Atoms ..... %ld\n",(long)GetNAtoms());

    if( GetNBonds() )
    {   
		PrintLog("Number of Bonds ..... %ld\n",(long)GetNBonds());
    }
   
    PrintLog("Number of SS Bridges ... %d\n\n", GetNSSBonds() );
    PrintLog("Number of H-Bonds ... %d\n",GetNHBonds());

    if( IsSecStructFound() || structsource == SourceNone)
    {   
        PrintLog("Number of Helices ... %d\n", GetNumFeatures(FeatHelix));
        PrintLog("Number of Strands ... %d\n", GetNumFeatures(FeatSheet));
        PrintLog("Number of Turns ..... %d\n", GetNumFeatures(FeatTurn));
    }
}

HaResidue* HaMolecule::AddChainAndResidue()
{
	HaChain* pch = AddChain(' ');
    HaResidue* pres = pch->AddResidue(1);
	pres->SetName( "RES" );
	return pres;        
}

void HaMolecule::Renumber(int start )
{
    HaChain  *chain;
    HaResidue  *group;
    int resno;

    ChainIteratorMolecule ch_itr(this);
	for(chain = ch_itr.GetFirstChain(); chain; chain = ch_itr.GetNextChain())
    {   
		resno = start;
		ResidueIteratorChain ritr_ch(chain);
		for( group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes() )
		{   
			group->serno = resno++;
		}
    }
}

AtomIntMap HaMolecule::GetAtomSeqNumMap()
{
	AtomIntMap at_seq_num_map;

	HaAtom* aptr = nullptr;
	AtomIteratorMolecule aitr(this);

	int i = 0;

	for (aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		at_seq_num_map[aptr] = i;
		i++;
	}
	return at_seq_num_map;
}

CAtomIntMap HaMolecule::GetAtomSeqNumMap() const
{
	CAtomIntMap at_seq_num_map_loc;

	const HaAtom* aptr = nullptr;
	AtomIteratorMolecule_const aitr(this);
	int i = 0;

	for (aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		at_seq_num_map_loc[aptr] = i;
		i++;
	}
	return at_seq_num_map_loc;
}


