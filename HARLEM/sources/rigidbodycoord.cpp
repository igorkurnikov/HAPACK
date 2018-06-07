/*! \file rigidbodycoord.cpp

    Classes to define Rigid Body Coordinates manipulations in HARLEM

   \author Igor Kurnikov 
   \date 2009

*/

#include "rapidxml.hpp"

#include "hacoord.h"
#include "haatgroup.h" 
#include "rigidbodycoord.h"

using namespace harlem;

RigidBodyCoord::RigidBodyCoord()
{
	n_obj = 0;
}

RigidBodyCoord::RigidBodyCoord( const RigidBodyCoord& ref )
{
	SetFrom(&ref);
}

RigidBodyCoord::~RigidBodyCoord()
{

}

harlem::Coord* RigidBodyCoord::clone()
{
	RigidBodyCoord* pcrd = new RigidBodyCoord(*this);
	return (harlem::Coord*)pcrd;
}

int RigidBodyCoord::SetFrom(const harlem::Coord* ref_crd)
{
	if( ref_crd->GetClassName() != "RigidBodyCoord") return FALSE;
	RigidBodyCoord* rb_crd = (RigidBodyCoord*) ref_crd;
	n_obj      = rb_crd->n_obj;
	crd_v      = rb_crd->crd_v;
	frozen_idx = rb_crd->frozen_idx; 
	return TRUE;
}


void RigidBodyCoord::SetNumObj( int n_obj_new )
{
	if( n_obj != n_obj_new && n_obj_new >= 0)
	{
		n_obj = n_obj_new;
		crd_v.resize(n_obj*6);
		frozen_idx.resize(n_obj*6);
		frozen_idx = 0;
	}
}

int RigidBodyCoord::GetNumObj() const
{
	return n_obj;
}

int RigidBodyCoord::GetNumCrd() const
{
	return crd_v.size();
}

HaVec_double RigidBodyCoord::AsVecDouble() const
{
	return crd_v;
}
	
int RigidBodyCoord::SetFromVecDouble(const HaVec_double& crd_v_new)
{
	int n_obj = GetNumObj();
	if( crd_v.size() != 6*n_obj )
	{
		PrintLog("Error in RigidBodyCoord::SetFromVecDouble() \n");
		PrintLog("Dimension of Double Vector = %d  doesn't correspond to the number of objects %d \n",n_obj, crd_v.size());
		return FALSE;
	}
	crd_v = crd_v_new;
	return TRUE;
}

double RigidBodyCoord::GetPhi(int iobj) const
{
	if(iobj < 0 || iobj >= n_obj ) return 0.0;
	return crd_v[iobj*6];		
}
    
double RigidBodyCoord::GetCosTheta(int iobj) const
{
	if(iobj < 0 || iobj >= n_obj ) return 0.0;
	return crd_v[iobj*6 + 1];		
}

double RigidBodyCoord::GetPsi(int iobj) const
{
	if(iobj < 0 || iobj >= n_obj ) return 0.0;
	return crd_v[iobj*6 + 2];			
}
    
double RigidBodyCoord::GetTransX(int iobj) const
{
	if(iobj < 0 || iobj >= n_obj ) return 0.0;
	return crd_v[iobj*6 + 3];			
}

double RigidBodyCoord::GetTransY(int iobj) const
{
	if(iobj < 0 || iobj >= n_obj ) return 0.0;
	return crd_v[iobj*6 + 4];	
}
    
double RigidBodyCoord::GetTransZ(int iobj) const
{
	if(iobj < 0 || iobj >= n_obj ) return 0.0;
	return crd_v[iobj*6 + 5];	
}

void RigidBodyCoord::SetPhi(int iobj, double phi_new ) 
{
	if(iobj < 0 || iobj >= n_obj )
	{
		throw "Error in RigidBodyCoord::SetPhi() \n  (iobj < 0 || iobj >= n_obj) \n";
	}
	crd_v[iobj*6] = phi_new;		
}
    
void RigidBodyCoord::SetCosTheta(int iobj, double cos_theta_new) 
{
	if(iobj < 0 || iobj >= n_obj )
	{
		throw "Error in RigidBodyCoord::SetCosTheta() \n  (iobj < 0 || iobj >= n_obj) \n";
	}
	if( cos_theta_new <= -1.0 || cos_theta_new >= 1.0 )
	{
		throw "Error in RigidBodyCoord::SetCosTheta() \n  (cos_theta_new <= -1.0 || cos_theta_new >= 1.0) \n";
	}
	crd_v[iobj*6 + 1] = cos_theta_new;		
}

void RigidBodyCoord::SetPsi(int iobj, double psi_new) 
{
	if(iobj < 0 || iobj >= n_obj ) 
	{
		throw "Error in RigidBodyCoord::SetPsi() \n  (i < 0 || i >= n_obj) \n";
	}
	crd_v[iobj*6 + 2] = psi_new;			
}
    
void RigidBodyCoord::SetTransX(int iobj, double x_new)
{
	if(iobj < 0 || iobj >= n_obj )
	{
		throw "Error in RigidBodyCoord::SetTransX() \n  (iobj < 0 || iobj >= n_obj) \n";
	}
	crd_v[iobj*6 + 3] = x_new;			
}

void RigidBodyCoord::SetTransY(int iobj, double y_new) 
{
	if(iobj < 0 || iobj >= n_obj )
	{
		throw "Error in RigidBodyCoord::SetTransY() \n  (iobj < 0 || iobj >= n_obj) \n";
	}
	crd_v[iobj*6 + 4] = y_new;	
}
    
void RigidBodyCoord::SetTransZ(int iobj, double z_new) 
{
	if(iobj < 0 || iobj >= n_obj )
	{
		throw "Error in RigidBodyCoord::SetTransZ() \n  (iobj < 0 || iobj >= n_obj) \n";
	}
	crd_v[iobj*6 + 5] = z_new;	
}

void RigidBodyCoord::FreezeObject(int iobj)
{
	if( iobj < 0 || iobj >= n_obj) 
	{
		throw "Error in RigidBodyCoord::FreezeObject() \n Invalid iobj value \n";
	}
	int j;
	for(j=0; j < 6; j++)
	{
		frozen_idx[iobj*6 + j] = TRUE;
	}
}

int  RigidBodyCoord::IsObjectFrozen(int iobj) const
{
	if( iobj < 0 || iobj >= n_obj) 
	{
		throw "Error in RigidBodyCoord::IsObjectFrozen() \n Invalid iobj value \n";
	}	
	int j;
	for(j=0; j < 6; j++)
	{
		if( frozen_idx[iobj*6 + j] ) return TRUE;
	}
	return FALSE;
}

int RigidBodyCoord::SetFromCurrAtomCrd(AtomContainer* at_cont, int iobj)
{
	if(at_cont == NULL )
	{
		PrintLog("Error in RigidBodyCoord::SetFromCurrAtomCoord() \n");
		PrintLog("AtomContainer pointer is NULL \n");
		return FALSE;
	}
	if( iobj < 0 || iobj  >= GetNumObj() )
	{
		PrintLog("Error in RigidBodyCoord::SetFromCurrAtomCoord() \n");
		PrintLog("iobj of AtomContainer %d is out of range, Num Objects = %d \n", iobj, GetNumObj() );
		return FALSE;
	}

	Vec3D trans_v;

	at_cont->GetPosEulerTrans( crd_v[iobj*6], crd_v[iobj*6+1], crd_v[iobj*6+2],trans_v);
	crd_v[iobj*6+3] = trans_v[0];
	crd_v[iobj*6+4] = trans_v[1];
	crd_v[iobj*6+5] = trans_v[2];

	return TRUE;
}

int RigidBodyCoord::SetFromCurrAtomCrd(vector<AtomContainer*> vec_at_cont)
{
	int n_obj = GetNumObj();
	if( vec_at_cont.size() != n_obj )
	{
		PrintLog("Error in RigidBodyCoord::SetFromCurrAtomCoord() \n");
		PrintLog("Invalid Number of Atom Containers %d ,Num Objects = %d \n", vec_at_cont.size(), n_obj );
		return FALSE;
	}
	int i;
	for(i = 0; i < n_obj; i++)
	{
		this->SetFromCurrAtomCrd(vec_at_cont[i],i);
	}
	return TRUE;
}

void RigidBodyCoord::FreezeCrd(int idx)
{
	if( idx < 0 || idx >= this->GetNumCrd() ) return;
	frozen_idx[idx] = TRUE;
}

int  RigidBodyCoord::IsCrdFrozen(int idx) const
{
	if( idx < 0 || idx >= this->GetNumCrd() ) return TRUE;
	return frozen_idx[idx];
}

int RigidBodyCoord::LoadFromStream(std::istream& is, const harlem::HashMap* popt_par )
{
	const LoadCrdOptions* popt = dynamic_cast<const LoadCrdOptions*>(popt_par);
	if( popt == NULL) popt = new LoadCrdOptions();

	char buf[256];

	int i_obj;
	int n_obj = this->GetNumObj();

	for( i_obj = 0; i_obj < n_obj; i_obj++ )
	{
		if( IsObjectFrozen(i_obj) && popt->ToLoadNotFrozenCrd() ) continue;

		int junk;
		double x,y,z;
		double phi, cos_theta, psi;

		is >> x >> y >> z;
		is >> phi >> cos_theta >> psi;

		if(is.fail())
		{
			std::string str_err = "Error in RigidBodyCoord::LoadFromStream() \n";
			str_err += " reading line: \n";
			str_err += buf; 
			str_err += "\n";
			throw str_err;
		}
		
		this->SetPhi(i_obj,phi);
		this->SetCosTheta(i_obj,cos_theta);
		this->SetPsi(i_obj,cos_theta);
		this->SetTransX(i_obj,x);
		this->SetTransY(i_obj,y);
		this->SetTransZ(i_obj,z);
	}
	return TRUE;
}

using namespace harlem;

int RigidBodyCoord::SaveToStream(std::ostream& os, const harlem::HashMap* popt_par ) const
{
	const SaveCrdOptions* popt = dynamic_cast<const SaveCrdOptions*>(popt_par);
	if( popt == NULL ) popt = new SaveCrdOptions();

	int i_obj;
	int n_obj = this->GetNumObj();

	for( i_obj = 0; i_obj < n_obj; i_obj++ )
	{
		if( IsObjectFrozen(i_obj) && popt->ToSaveNotFrozenCrd() ) continue;

		double x = this->GetTransX(i_obj);
		double y = this->GetTransY(i_obj);
		double z = this->GetTransZ(i_obj);
		double phi = this->GetPhi(i_obj);
		double cos_theta = this->GetCosTheta(i_obj);
		double psi = this->GetPsi(i_obj);

		os << x;
		os << y;
		os << z;
		os << phi;
		os << cos_theta;
		os << psi;
		os << std::endl;

		if(os.fail())
		{
			throw "Error in RigidBodyCoord::SaveToString \n";
		}
	}
	return TRUE;
}



RigidBodyCoordDiscretized::RigidBodyCoordDiscretized()
{
	
}

RigidBodyCoordDiscretized::~RigidBodyCoordDiscretized()
{
	
}

RigidBodyCoordDiscretized::RigidBodyCoordDiscretized(const RigidBodyCoordDiscretized& ref): RigidBodyCoord(ref)
{
	crd_v_int = ref.crd_v_int;
	dim_crd = ref.dim_crd;
	crd_min = ref.crd_min;
	crd_max = ref.crd_max;
}

harlem::Coord* RigidBodyCoordDiscretized::clone() 
{
	RigidBodyCoordDiscretized* pcrd = new RigidBodyCoordDiscretized(*this);
	return (harlem::Coord*)pcrd;
}

int RigidBodyCoordDiscretized::SetFrom(const harlem::Coord* pcrd)
{
	if( pcrd->GetClassName() != "RigidBodyCoordDiscretized") return FALSE;
	const RigidBodyCoordDiscretized* p_from = (RigidBodyCoordDiscretized*) pcrd;
	RigidBodyCoord* p_to   = (RigidBodyCoord*) this; 
	p_to->SetFrom(p_from);

	crd_v_int = p_from->crd_v_int;
	dim_crd = p_from->dim_crd;
	crd_min = p_from->crd_min;
	crd_max = p_from->crd_max;
	return TRUE;
}

void RigidBodyCoordDiscretized::SetNumObj(int n_obj_new)
{
	if( GetNumObj() == n_obj_new || n_obj_new < 0 ) return;
	
	RigidBodyCoord::SetNumObj(n_obj_new);
	
	crd_v_int.resize(n_obj_new*6);
	dim_crd.resize(n_obj_new*6);
	crd_min.resize(n_obj_new*6);
	crd_min.resize(n_obj_new*6);
	SetStandardLimits();
}

void RigidBodyCoordDiscretized::ConvertDiscrCrdToFloat()
{
	int i;
	int n = this->GetNumCrd();
	for( i = 0; i < n; i++)
	{
		double val = crd_min[i];
		if( dim_crd[i] > 1)
		{
			val = crd_min[i] + (crd_v_int[i]*(crd_max[i] - crd_min[i]))/(dim_crd[i] - 1);
		}
		crd_v[i] = val;
	}
}

void RigidBodyCoordDiscretized::ConvertFloatCrdToDiscr()
{
	int i;
	int n = this->GetNumCrd();
	for( i = 0; i < n; i++)
	{
		int ival = 0;
		if( crd_max[i] > crd_min[i] && dim_crd[i] > 1 ) 
		{
			ival = (int)( (dim_crd[i] - 1)*((crd_v[i] - crd_min[i])/(crd_max[i] - crd_min[i])));
			if(ival < 0) ival = 0; 
			if(ival >= dim_crd[i] ) ival = dim_crd[i] - 1;
		}
		crd_v_int[i] = ival;
	}
}	
	
void RigidBodyCoordDiscretized::SetStandardLimits()
{
	int iobj;
	for(iobj = 0; iobj < n_obj; iobj++)
	{
		SetDiscrNumForCrd(6*iobj   ,10);
		SetDiscrNumForCrd(6*iobj+1 ,10);
		SetDiscrNumForCrd(6*iobj+2 ,10);
		SetDiscrNumForCrd(6*iobj+3  ,10);
		SetDiscrNumForCrd(6*iobj+4  ,10);
		SetDiscrNumForCrd(6*iobj+5  ,10);
		SetLimits(6*iobj , -PI, +PI);
		SetLimits(6*iobj+1, -2.0, +2.0);
		SetLimits(6*iobj+2, -PI, +PI);
		SetLimits(6*iobj+3, 0.0, 10.0);
		SetLimits(6*iobj+4, 0.0, 10.0);
		SetLimits(6*iobj+5, 0.0, 10.0);
	}
}

int RigidBodyCoordDiscretized::SetDiscrNumForCrd(int idx, int npt)
{
	if( idx < 0 || idx >= (6*this->n_obj)) return FALSE;
	dim_crd[idx] = npt;
	return TRUE;
}

int RigidBodyCoordDiscretized::SetLimits(int idx, double amin, double amax)
{
	if( idx < 0 || idx >= (6*this->n_obj)) return FALSE;
	crd_min[idx] = amin;
	crd_max[idx] = amax;
	return TRUE;
}

