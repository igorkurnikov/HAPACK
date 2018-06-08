/*! \file  halocorb.cpp
 
    Classes 
    to define Local Orbital object in HARLEM.
    implementation
 
    \author Igor Kurnikov  
    \date 1997-2004

*/
#include <mpi.h>

#include "haqchem.h"
#include "halocorb.h"
#include "haintengine.h"
#include "tinyxml.h"


LinCombOrb3D::LinCombOrb3D()
{
	SetStdParams();
}	

LinCombOrb3D::LinCombOrb3D(ArrayOrb3D* new_bas, const double* coef_new, int nv)
{
	bas = new_bas;
	if(bas == NULL) return;
	int nb = bas->GetNBfunc();

	coef.newsize(nb,nv);
	coef.copy(coef_new);
}

LinCombOrb3D::~LinCombOrb3D()
{
	Clear();
}

void LinCombOrb3D::SetStdParams()
{
	internal_basis = FALSE;
	coef.newsize(0,0);
	bas = NULL;
}	

int LinCombOrb3D::GetNOrbs() const
{
	return coef.num_cols();
}

int LinCombOrb3D::CreateEmptyOrbs(int n_orb, ArrayOrb3D* new_bas)
{
	if(new_bas == NULL) return FALSE;

	bas = new_bas;
	int nb = bas->GetNBfunc();
    coef.newsize(nb,n_orb);	
	coef = 0.0;
	return TRUE;
}

int LinCombOrb3D::AddOrbs(LinCombOrb3D& orbs)
{
	if(bas == NULL) 
	{
		bas = orbs.bas;
	}
	if(bas != orbs.bas)
	{
		PrintLog("Error in LinCombOrb3D::AddOrbs() \n");
		PrintLog(" Added orbitals are expanded using a different basis set \n");
		return FALSE;
	}
	if(bas == NULL)
	{
		return FALSE;
	}
	int nadd = orbs.GetNOrbs();

	if(nadd == 0) return TRUE;
	int nold = GetNOrbs();
	int nb = bas->GetNBfunc();
	if(nb == 0) return TRUE;
	
	if(nold == 0) 
	{
		coef = orbs.coef;
	}
	else
	{
		HaMat_double tmp(nb,nold+nadd);
		memcpy(tmp.begin(),coef.begin(), sizeof(double)*nb*nold);
        memcpy(&tmp(1,nold+1),orbs.coef.begin(), sizeof(double)*nb*nadd);
		coef = tmp;
	}
	
	VecPtr hosts_old = at_ptr;

	int norb = nold+nadd;

	at_ptr.resize(norb);
	int i;
	for(i = 0; i < norb; i++)
	{
		if( i < nold )
		{
			if( hosts_old.size() > i)  
			{
				at_ptr[i] = hosts_old[i];
			}
			else
			{
				at_ptr[i] = NULL;
			}
		}
		else
		{
			int idx = i - nold;
			if( orbs.at_ptr.size() > idx )
			{
				at_ptr[i] = orbs.at_ptr[idx];
			}
			else
			{
				at_ptr[i] = NULL;
			}
		}
	}
	return TRUE;
}

void LinCombOrb3D::SetOrbLabel(int idx, const std::string& orb_lbl)
{
	int norb = this->GetNOrbs();
	if( idx < 0 || idx >= norb )
	{
		PrintLog(" Error in LinCombOrb3D::SetOrbLabel() \n");
		PrintLog(" Invalid orbital index \n");
		return;
	}
	int nold = ids.size();
	if( nold != norb )
	{
		ids.resize( norb, "" ); 
	}
	ids[idx] = orb_lbl;
}

int LinCombOrb3D::ProjectToBasis(HaMat_double& coef_new, const ArrayOrb3D* basis_new)
{
	coef_new.resize(0,0);

	if(this->bas == NULL)
	{
		PrintLog("Error in LinCombOrb3D::ProjectToBasis() \n");
		PrintLog("parent basis set is NULL \n");
		return FALSE;
	}


	if(basis_new == NULL)
	{
		PrintLog("Error in LinCombOrb3D::ProjectToBasis() \n");
		PrintLog("new basis is NULL \n");
		return FALSE;
	}

	if( basis_new == this->bas) 
	{	
		coef_new = this->coef;
		return TRUE;
	}

	std::string class_name_new  = basis_new->GetClassName();
	std::string class_name_old  = bas->GetClassName();

	int nb_old = this->bas->GetNBfunc();
	int nb_new = basis_new->GetNBfunc();

	int nf = this->coef.num_cols();

	if( class_name_new == "GauBasisSet" && class_name_old == "GauBasisSet" )
	{
		IntIntMap frag_bas_fun_map;
        ((GauBasisSet*) basis_new)->MatchBasisSet((GauBasisSet*)bas,frag_bas_fun_map);

		if(frag_bas_fun_map.size() != nb_old )
		{
			PrintLog("Warning in LinCombOrb3D::ProjectToBasis() \n");
			PrintLog("Only %d from %d functions of the old basis set are mapped to the new basis set\n",
				      frag_bas_fun_map.size(), nb_old );

		}
		coef_new.newsize(nb_new, nf);
		coef_new = 0.0;

		IntIntMap::iterator itr;
		for(itr = frag_bas_fun_map.begin(); itr != frag_bas_fun_map.end(); itr++ )
		{
			int i1 = (*itr).first;
			int j1 = (*itr).second;
			
			int k;
			for( k = 0; k < nf; k++)
			{
				coef_new.r0(j1,k) = coef.r0(i1,k);
			}
		}
		return TRUE;
	}
	return FALSE;
}

int LinCombOrb3D::CalcOvlpMat(LinCombOrb3D* pbas1, LinCombOrb3D* pbas2, HaMat_double& ovlp_mat)
{
	if( pbas1 == NULL || pbas2 == NULL) return FALSE;

	int nf1 = pbas1->GetNOrbs();
	int nf2 = pbas2->GetNOrbs();

	int nb1 = pbas1->bas->GetNBfunc();
	int nb2 = pbas2->bas->GetNBfunc();	

	HaMat_double bovlp(nb1,nb2);

	int ires = HaBasisSet::CalcOvlpMat(pbas1->bas, pbas2->bas, bovlp);

	if(ires == FALSE) return FALSE;

	ovlp_mat.newsize(nf1,nf2);
	
	HaMat_double tmp;

	matmult_T1(tmp,pbas1->coef,bovlp);
	matmult(ovlp_mat,tmp,pbas2->coef);

	return TRUE;
}

int LinCombOrb3D::Eval1eOp(LinCombOrb3D* pbas1, LinCombOrb3D* pbas2, const HaMat_double& bop_mat, HaMat_double& op_mat)
{
	if( pbas1 == NULL || pbas2 == NULL) return FALSE;

	int nf1 = pbas1->GetNOrbs();
	int nf2 = pbas2->GetNOrbs();

	if( nf1 == 0 || nf2 == 0 )
	{
		op_mat.clear();
		return FALSE;
	}

	int nb1 = pbas1->bas->GetNBfunc();
	int nb2 = pbas2->bas->GetNBfunc();	

	if(nb1 != nb2 || bop_mat.num_rows() != nb1 || bop_mat.num_cols() != nb1)
	{
		PrintLog("LinCombOrb3D::CalcOvlpFromBasOvlp() \n");
	    PrintLog("Invalid size of basis set overlap matrix %d X %d ; nb1 = %d , nb2 = %d \n",
			bop_mat.num_rows(),bop_mat.num_cols(),nb1,nb2 );
		return FALSE;
	}

	op_mat.newsize(nf1,nf2);
	
	HaMat_double tmp;

	matmult_T1(tmp,pbas1->coef,bop_mat);
	matmult(op_mat,tmp,pbas2->coef);

	return TRUE;
}

void LinCombOrb3D::Clear()
{
	if( bas != NULL && internal_basis) delete bas;
	bas = NULL;
	ids.clear();
	at_ptr.clear();
	SetStdParams();
}

bool LinCombOrb3D::IsEmpty()
{
	if( GetNOrbs() > 0) return false;
	return true;
}

int LinCombOrb3D::GetOrbIdxByID(const char* id)
{
	int norb = GetNOrbs();
	int i;
	if(ids.size() != norb )
	{
		PrintLog("Error in: LinCombOrb3D::GetOrbIdxByID() \n");
		PrintLog("Number of orbitals in the array %d not equal to the number of labels %d \n",
			     norb, ids.size());
		return -1;
	}

	for(i = 0; i < norb; i++)
	{
		if( ids[i] == id) return i;
	}

    PrintLog("Error in: LinCombOrb3D::GetOrbIdxByID() \n");
	PrintLog("No orbital with label %s \n",id);

	return -1;
}


std::string LinCombOrb3D::GetLabel(int idx)
{
	if(idx < ids.size() && idx >= 0)
	{
		return ids[idx].c_str();
	}
	return "LIN_COMB_ORB";
}

Vec3D* LinCombOrb3D::GetHostPt(int idx)
{
	if( idx < 0 || idx >= at_ptr.size())
	{
		return NULL;
	}
	return (Vec3D*) at_ptr[idx];
}


int LinCombOrb3D::TransferBetweenAtoms(PtrPtrMap& pt_corr_map)
{
	int ires1 = TRUE;
	int ires2 = TRUE;

	if(bas != NULL) ires1 = bas->TransferBetweenAtoms(pt_corr_map);

	int na = at_ptr.size();
	int i;
	for( i = 0; i < na; i++)
	{
		HaAtom* host_new = (HaAtom*) pt_corr_map.GetVal( at_ptr[i]);
		if( host_new != NULL)
		{
			at_ptr[i] = host_new;
		}
		else
		{
			ires2 = FALSE;
		}
	}
	if( !ires1 || !ires2 )
	{
		PrintLog("Error in LinCombOrb3D::TransferBetweenAtoms() \n");
		PrintLog(" One of the host atoms not found in point correspondence map \n");
		return FALSE;
	}
	return TRUE;
}

TiXmlElement* LinCombOrb3D::AddXml(TiXmlElement* parent_element, const char* name, int option) const
{
	if( parent_element == NULL) return NULL;

	TiXmlElement* lin_comb_element;

	if( strlen(name) > 0 )
	{
		lin_comb_element = new TiXmlElement(name);	
		lin_comb_element->SetAttribute("TYPE","LinCombOrb3D");
	}
	else
	{
		lin_comb_element = new TiXmlElement("LinCombOrb3D");	
	}

	parent_element->LinkEndChild(lin_comb_element);
	
	TiXmlElement* coef_element = coef.AddXml(lin_comb_element, "coef");

	if(bas != NULL) 
	{
		TiXmlElement* bas_element = bas->AddXml(lin_comb_element, "bas");
	}

	return lin_comb_element;
}

int LinCombOrb3D::LoadXml(const TiXmlElement* xml_element, int option )
{
	if( xml_element == NULL) return FALSE;

	const TiXmlElement* coef_elem = xml_element->FirstChildElement("coef");

	if(coef_elem != NULL) coef.LoadXml( coef_elem );

	const TiXmlElement* bas_elem = xml_element->FirstChildElement("bas");
	
	if(bas_elem != NULL )
	{
		if( internal_basis ) delete bas;
		bas = NULL;
		const char* sattr = bas_elem->Attribute("TYPE");
		if( sattr) bas = ArrayOrb3D::CreateObjectWithType(sattr);
		if(bas != NULL) bas->LoadXml( bas_elem );
	}
	return TRUE;
}


const Vec3D* LinCombOrb3D::GetHostPt(int idx) const
{
	if( idx < 0 || idx >= at_ptr.size())
	{
		return NULL;
	}
	return (const Vec3D*) at_ptr[idx];
}

int LinCombOrb3D::TrCoefRot(const HaMat_double& rot_mat)
{
    HaMat_double tmat, coef_mod;
    bas->GetTransfMat(tmat, rot_mat);
	matmult(coef_mod,tmat,coef);
	coef = coef_mod;
	return TRUE;
}


