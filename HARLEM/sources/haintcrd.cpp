/*! \file haintcrd.cpp

    Classes to define Internal Coordinates representation in HARLEM

   \author Igor Kurnikov 
   \date 2009

*/
#include <sstream>
#include <math.h>

#include <boost/algorithm/string.hpp>

#include "rapidxml.hpp"

#include "hamolecule.h"
#include "hamolset.h"
#include "hacoord.h"
#include "haintcrd.h"

#define HAINTCRD_CPP

ZMatLoadOptions  ZMatCrd::opt_read_default;
ZMatSaveOptions ZMatCrd::opt_save_default;

CrdAssignRule::CrdAssignRule()
{
	priority = 0;
}

int CrdAssignRule::GetPriority() const 
{ 
	return priority; 
}

void CrdAssignRule::SetPriority( int priority_new ) 
{ 
	priority = priority_new; 
}

SingleAtomCrdRule::SingleAtomCrdRule( HaAtom* p_mng_atom_par )
{
	p_mng_atom = p_mng_atom_par;
}

SingleAtomCrdRule::~SingleAtomCrdRule()
{

}

HaAtom* SingleAtomCrdRule::GetManagedAtom() 
{ 
	return p_mng_atom;
}

int SingleAtomCrdRule::SetManagedAtom( HaAtom* aptr )
{
	if( aptr == NULL ) return FALSE;
	p_mng_atom = aptr;
	return TRUE;
}

int SingleAtomCrdRule::SaveXMLToStream(std::ostream& os, const harlem::HashMap* popt ) const
{
	os << "<single_atom_crd_rule atom=\"" << p_mng_atom->GetRef() << "\" />" << std::endl;
	if( os.fail() ) return FALSE;
	return TRUE;
}

SameAtomCrdRule::SameAtomCrdRule(HaAtom* p_mng_atom_new, HaAtom* p_ref_atom_new):
SingleAtomCrdRule(p_mng_atom_new)
{
	p_ref_atom  = p_ref_atom_new;
	priority = 0;
}

SameAtomCrdRule::~SameAtomCrdRule()
{

}

int SameAtomCrdRule::SetManagedAtomCrd()
{
	p_mng_atom->SetCoordFrom(*p_ref_atom);
	return TRUE;
}

int SameAtomCrdRule::ReplaceRefPt( Vec3D* p_ref_old_par, Vec3D* p_ref_new_par)
{
	try
	{
		HaAtom* p_ref_new = dynamic_cast<HaAtom*>(p_ref_new_par);
		HaAtom* p_ref_old = dynamic_cast<HaAtom*>(p_ref_old_par);

		if( p_ref_new == NULL ) throw std::runtime_error(" p_ref_new == NULL ");
		if( p_ref_old == NULL ) throw std::runtime_error(" p_ref_old == NULL ");
		if( p_ref_atom == p_ref_old )
		{
			p_ref_atom = p_ref_new;
			return TRUE;
		}
		throw std::runtime_error(" p_ref_old is not equal to reference point of the rule ");
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in SameAtomCrdRule::ReplaceRefPt() \n");
		PrintLog(" %s\n", ex.what());
		return FALSE;
	}
	return TRUE;
}

void SameAtomCrdRule::SetX( double val )
{
	if( p_ref_atom) p_ref_atom->SetX(val);
}

void SameAtomCrdRule::SetY( double val )
{
	if( p_ref_atom) p_ref_atom->SetY(val);
}

void SameAtomCrdRule::SetZ( double val )
{
	if( p_ref_atom) p_ref_atom->SetZ(val);
}

int SameAtomCrdRule::SaveXMLToStream(std::ostream& os, const harlem::HashMap* popt ) const
{
	os << "<same_atom_crd_rule atom=\"" << p_mng_atom->GetRef() << "\" />" << std::endl;
	if( os.fail() ) return FALSE;
	return TRUE;
}


Pt2CrdRule::Pt2CrdRule(HaAtom* p_mng_atom_new, HaAtom* p_ref_1_new, HaAtom* p_ref_2_new):
SingleAtomCrdRule(p_mng_atom_new)
{
	p_ref_1    = p_ref_1_new;
	p_ref_2    = p_ref_2_new;

	priority = 1;
	bond_len = 0.0;
}

Pt2CrdRule::~Pt2CrdRule()
{

}

int Pt2CrdRule::SetManagedAtomCrd()
{
	if( bond_len < 0.1 )
	{
		PrintLog(" Warning in Pt2CrdRule::SetManagedAtomCrd() \n");
		PrintLog(" bond_len < 0,1 ");
	}
//	PrintLog("Pt2CrdRule  mng_atom = %s ref_1 = %s  ref_2 = %s \n",
//		      (p_mng_atom->GetRef()).c_str(), (p_ref_1->GetRef()).c_str(), (p_ref_2->GetRef()).c_str());
	Vec3D v1_2 = (*p_ref_2) - (*p_ref_1);
	v1_2.normalize();
	v1_2.Scale(bond_len);
	(Vec3D&)(*p_mng_atom) = (*p_ref_1) + v1_2;
	return TRUE;
}

int Pt2CrdRule::ReplaceRefPt( Vec3D* p_ref_old_par, Vec3D* p_ref_new_par)
{
	try
	{
		HaAtom* p_ref_new = dynamic_cast<HaAtom*>(p_ref_new_par);
		HaAtom* p_ref_old = dynamic_cast<HaAtom*>(p_ref_old_par);

		if( p_ref_new == NULL ) throw std::runtime_error(" p_ref_new == NULL ");
		if( p_ref_old == NULL ) throw std::runtime_error(" p_ref_old == NULL ");
		if( p_ref_1 == p_ref_old )
		{
			p_ref_1 = p_ref_new;
			return TRUE;
		}
		if( p_ref_2 == p_ref_old )
		{
			p_ref_2 = p_ref_new;
			return TRUE;
		}
		throw std::runtime_error(" p_ref_old is not equal to any of reference points of the rule ");
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in Pt2CrdRule::ReplaceRefPt() \n");
		PrintLog(" %s\n", ex.what());
		return FALSE;
	}
	return TRUE;
}

void Pt2CrdRule::SetBondLen(double bond_len_new)
{
	bond_len = bond_len_new;
}

int Pt2CrdRule::SaveXMLToStream(std::ostream& os, const harlem::HashMap* popt ) const
{
	char buf[40];
	os << "<pt2_crd_rule atom=\"" << p_mng_atom->GetRef() << "\" ";
	os << " ref1=\"" << p_ref_1->GetRef() << "\"  ref2=\"" << p_ref_2->GetRef() << "\" ";
	sprintf(buf," blen=\"%10.5f\" ",bond_len); os << buf;
	os << "/>" << std::endl;
	if( os.fail() ) return FALSE;
	return TRUE;
}

Pt3CrdRule::Pt3CrdRule(HaAtom* p_mng_atom_new, Vec3D* p_ref_1_new, Vec3D* p_ref_2_new, Vec3D* p_ref_3_new):
SingleAtomCrdRule(p_mng_atom_new)
{
	p_ref_1    = p_ref_1_new;
	p_ref_2    = p_ref_2_new;
	p_ref_3    = p_ref_3_new;
	SetParFromPtPos(p_mng_atom, p_ref_1, p_ref_2, p_ref_3);
}

Pt3CrdRule::~Pt3CrdRule()
{
	
}

Vec3D pt2_std(0.000, 0.0,   1.0);
Vec3D pt3_std(1.0,   0.000, 0.0);

int Pt3CrdRule::SetManagedAtomCrd()
{
	if( p_ref_3 != NULL )
	{
		Vec3D::SetAtomPos(p_mng_atom,p_ref_1,p_ref_2,p_ref_3,bond_len,val_ang,dih_ang);
	}
	else if( p_ref_2 != NULL )
	{
//		PrintLog("Pt3CrdRule  mng_atom = %s ref_1 = %s  ref_2 = %s \n",
//		      (p_mng_atom->GetRef()).c_str(), (((HaAtom*)p_ref_1)->GetRef()).c_str(), (((HaAtom*)p_ref_2)->GetRef()).c_str());
		Vec3D::SetAtomPos(p_mng_atom,p_ref_1,p_ref_2,&pt3_std,bond_len,val_ang,dih_ang);
	}
	else
	{
//		PrintLog("Pt3CrdRule  mng_atom = %s ref_1 = %s  \n",
//		      (p_mng_atom->GetRef()).c_str(), (((HaAtom*)p_ref_1)->GetRef()).c_str() );
		Vec3D::SetAtomPos(p_mng_atom,p_ref_1,&pt2_std,&pt3_std, bond_len, val_ang, dih_ang);
	}
	return TRUE;
}

int Pt3CrdRule::SetParFromPtPos(const Vec3D* pt_mng, const Vec3D* pt_ref_1, const Vec3D* pt_ref_2,const Vec3D* pt_ref_3)
{
	double blen = Vec3D::CalcDistance(pt_mng,pt_ref_1);
	SetBondLen(blen);
	
	if( pt_ref_2 != NULL )
	{
		SetValAng(Vec3D::CalcAngle(pt_mng,pt_ref_1,pt_ref_2));
		if( pt_ref_3 != NULL )
		{
			SetDihAng( Vec3D::CalcTorsion(pt_mng,pt_ref_1,pt_ref_2,pt_ref_3) );
		}
		else
		{
			SetDihAng( 0.0 );
//			SetDihAng( Vec3D::CalcTorsion(pt_mng,pt_ref_1,pt_ref_2,&pt3_std) );
		}
	}
	else
	{
		SetValAng( 0.0 );
		SetDihAng( 0.0 );
//		SetValAng( Vec3D::CalcAngle(pt_mng,pt_ref_1,&pt2_std) );
//		SetDihAng( Vec3D::CalcTorsion(pt_mng,pt_ref_1,&pt2_std,&pt3_std) );
	}
	return TRUE;
}

int Pt3CrdRule::SetParFromCurPos()
{
	return SetParFromPtPos(p_mng_atom, p_ref_1, p_ref_2, p_ref_3);
}

int Pt3CrdRule::SetBondLen(double bond_len_new)
{
	bond_len = bond_len_new;
	return TRUE;
}

int Pt3CrdRule::SetValAng (double val_ang_new)
{
	val_ang = val_ang_new;
	return TRUE;
}

int Pt3CrdRule::SetDihAng(double dih_ang_new)
{
	dih_ang = dih_ang_new;
	return TRUE;
}

double Pt3CrdRule::GetBondLen() const
{
	return bond_len;
}

double Pt3CrdRule::GetValAng() const
{
	return val_ang;
}

double Pt3CrdRule::GetDihAng() const
{
	return dih_ang;
}

Vec3D* Pt3CrdRule::GetRefPt1()
{
	return p_ref_1;
}

Vec3D* Pt3CrdRule::GetRefPt2()
{
	return p_ref_2;
}

Vec3D* Pt3CrdRule::GetRefPt3()
{
	return p_ref_3;
}

int Pt3CrdRule::ReplaceRefPt( Vec3D* p_ref_old, Vec3D* p_ref_new)
{
	try
	{
		if( p_ref_new == NULL ) throw std::runtime_error(" p_ref_new == NULL ");
		if( p_ref_old == NULL ) throw std::runtime_error(" p_ref_old == NULL ");
		if( p_ref_1 == p_ref_old )
		{
			p_ref_1 = p_ref_new;
			return TRUE;
		}
		if( p_ref_2 == p_ref_old )
		{
			p_ref_2 = p_ref_new;
			return TRUE;
		}
		if( p_ref_3 == p_ref_old )
		{
			p_ref_3 = p_ref_new;
			return TRUE;
		}
		throw std::runtime_error(" p_ref_old not equal any of reference points of the rule ");
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in Pt3CrdRule::ReplaceRefPt() \n");
		PrintLog(" %s\n", ex.what());
		return FALSE;
	}
	return TRUE;
}

int Pt3CrdRule::SaveXMLToStream(std::ostream& os, const harlem::HashMap* popt ) const
{
	char buf[40];
	os << "<pt3_crd_rule atom=\"" << p_mng_atom->GetRef() << "\" ";
	HaAtom* aptr1 = dynamic_cast<HaAtom*>(p_ref_1);
	HaAtom* aptr2 = dynamic_cast<HaAtom*>(p_ref_2);
	HaAtom* aptr3 = dynamic_cast<HaAtom*>(p_ref_3);

	os << " ref1=\"";
	if( aptr1 != NULL ) os << aptr1->GetRef();
	os << "\" ";

	os << " ref2=\"";
	if( aptr2 != NULL ) os << aptr2->GetRef();
	os << "\" ";

	os << " ref3=\"";
	if( aptr3 != NULL ) os << aptr3->GetRef();
	os << "\" ";

	sprintf(buf, " blen=\"%10.5f\"  vang=\"%10.5f\"  dih=\"%10.5f\" ", 
		      GetBondLen(), GetValAng()*RAD_TO_DEG, GetDihAng()*RAD_TO_DEG ); os << buf;
	os << "/>" << std::endl;
	if( os.fail() ) return FALSE;
	return TRUE;
}


FixedCrdRule::FixedCrdRule(HaAtom* p_mng_atom_par ):
SingleAtomCrdRule(p_mng_atom_par)
{  
	priority = 0;
	crd.SetX( p_mng_atom->GetX() );
	crd.SetY( p_mng_atom->GetY() );
	crd.SetZ( p_mng_atom->GetZ() );
} 
	
FixedCrdRule::~FixedCrdRule() 
{
	
}

int FixedCrdRule::SetManagedAtomCrd() 
{
	p_mng_atom->SetCoordFrom(crd);
	return TRUE; 
}

void FixedCrdRule::SetCrd( const Vec3D& crd_new)
{
	crd = crd_new;
}

void FixedCrdRule::SetParFromCurPos()
{
	crd.SetCoordFrom(*p_mng_atom);
}

void FixedCrdRule::SetX( double val )
{
	crd.SetX(val);
}
	
void FixedCrdRule::SetY( double val )
{
	crd.SetY(val);
}
	
void FixedCrdRule::SetZ( double val )
{
	crd.SetZ(val);
}

double FixedCrdRule::GetX() const
{
	return crd.GetX();
}

double FixedCrdRule::GetY() const
{
	return crd.GetY();
}

double FixedCrdRule::GetZ() const
{
	return crd.GetZ();
}

int FixedCrdRule::SaveXMLToStream(std::ostream& os, const harlem::HashMap* popt ) const
{
	char buf[80];
	os << "<fixed_crd_rule atom=\"" << p_mng_atom->GetRef() << "\" />";
	sprintf(buf," x=\"%16.9f\" y=\"%16.9f\" z=\"%16.9f\" ", p_mng_atom->GetX(), p_mng_atom->GetY(), p_mng_atom->GetZ());
	os << buf;
	os << " />" << std::endl;
	if( os.fail() ) return FALSE;
	return TRUE;
}

int FixedCrdRule::ReplaceRefPt( Vec3D* p_ref_old, Vec3D* p_ref_new)
{
	return TRUE;
}



DihedralAngleCoord::DihedralAngleCoord( HaAtom* new_aptr1, HaAtom* new_aptr2, HaAtom* new_aptr3, HaAtom* new_aptr4)
{
	aptr1 = new_aptr1;
	aptr2 = new_aptr2;
	aptr3 = new_aptr3;
	aptr4 = new_aptr4;

	FindMovingAtoms();
}

DihedralAngleCoord::~DihedralAngleCoord()
{

}

harlem::Coord* DihedralAngleCoord::clone() 
{
	DihedralAngleCoord* crd_new = new DihedralAngleCoord(aptr1,aptr2,aptr3,aptr4); 
	return crd_new;
}

HaVec_double DihedralAngleCoord::AsVecDouble() const
{
	HaVec_double v(1);
	v[0] = Vec3D::CalcDihedral(aptr1,aptr2,aptr3,aptr4);
	return v;
}

int DihedralAngleCoord::SetFrom(const Coord* pcrd)
{
	if( pcrd->GetClassName() != "DihedralAngleCoord") return FALSE;
	DihedralAngleCoord* p_da_crd = (DihedralAngleCoord*) pcrd;

	aptr1 = p_da_crd->aptr1;
	aptr2 = p_da_crd->aptr2;
	aptr3 = p_da_crd->aptr3;
	aptr4 = p_da_crd->aptr4;

	FindMovingAtoms();
	return TRUE;
}

int DihedralAngleCoord::SetFromVecDouble(const HaVec_double& dbl_vec)
{
	if(dbl_vec.size() < 1) return FALSE;

	// FINISH!!!
	return FALSE;
}

int DihedralAngleCoord::LoadFromStream(std::istream& is, const harlem::HashMap* popt )
{
	// FINISH!!

	return FALSE;
}

int DihedralAngleCoord::SaveToStream(std::ostream&  os, const harlem::HashMap* popt ) const
{
	// FINISH!!

	return FALSE;
}

double DihedralAngleCoord::GetDihVal()
{
	return Vec3D::CalcDihedral(aptr1,aptr2,aptr3,aptr4);
}

int DihedralAngleCoord::SetDihVal(double new_dih_val)
{
	return FALSE;
}

int DihedralAngleCoord::FindMovingAtoms()
{
	if(aptr2 == NULL || aptr3 == NULL)
	{
		return FALSE;
	}
	
	moving_atoms.clear();

	set<HaAtom*> found_atoms;

	HaAtom* aptr;
	AtomGroup bonded_atoms;

	aptr3->GetBondedAtoms(bonded_atoms);
	AtomIteratorAtomGroup aitr_bonded_atoms(&bonded_atoms);

	for( aptr = aitr_bonded_atoms.GetFirstAtom(); aptr != NULL; aptr = aitr_bonded_atoms.GetNextAtom() )
	{
		if(aptr == aptr2) continue;
	}

	return TRUE;
}

ElemCrd::ElemCrd()
{
	type = UNDEF_ELEM_CRD;
	p_rule = NULL;
	tag = "";
	frozen_flag = false;
}

ElemCrd::ElemCrd( ElemCrdType type_par , CrdAssignRule* p_rule_par )
{
	type = type_par;
	p_rule = p_rule_par;
	tag = "";
	frozen_flag = false;
}

ElemCrd::~ElemCrd()
{

}
	
void ElemCrd::SetValue( double val )
{
	try
	{
		if( p_rule == NULL ) throw std::runtime_error(" p_rule == NULL ");
		if( type == X_ELEM_CRD || type == Y_ELEM_CRD || type == Z_ELEM_CRD )
		{
			FixedCrdRule* p_fixed_crd_rule = dynamic_cast< FixedCrdRule* >(p_rule);
			if( p_fixed_crd_rule != NULL )
			{
				if( type == X_ELEM_CRD ) p_fixed_crd_rule->SetX( val );
				if( type == Y_ELEM_CRD ) p_fixed_crd_rule->SetY( val );
				if( type == Z_ELEM_CRD ) p_fixed_crd_rule->SetZ( val );
				return;
			}
			throw std::runtime_error(" elem crd X, Y or Z  however p_rule type is not  FixedCrdRule* ");
		}
		Pt3CrdRule* p_pt3_rule = dynamic_cast< Pt3CrdRule* >(p_rule);
		if( p_pt3_rule != NULL )
		{
			if( type == LEN_ELEM_CRD ) 
			{
				p_pt3_rule->SetBondLen( val );
			}
			else if( type == ANG_ELEM_CRD )
			{
				p_pt3_rule->SetValAng( val );
			}
			else if ( type == DIH_ELEM_CRD )
			{
				p_pt3_rule->SetDihAng( val );
			}
			else
			{
				throw std::runtime_error(" Invalid ElemCrd type for Pt3CrdRule rule ");
			}
			return;
		}
	}
	catch( const std::exception& ex )
	{
		PrintLog("Error in ElemCrd::SetValue() \n");
		PrintLog("%s\n",ex.what() );
	}
}
	
double ElemCrd::GetValue() const
{
	try
	{
		if( p_rule == NULL ) throw std::runtime_error(" p_rule == NULL ");
		FixedCrdRule* p_fixed_crd_rule = dynamic_cast< FixedCrdRule* >(p_rule);
		if( p_fixed_crd_rule != NULL )
		{
			if( type == X_ELEM_CRD ) return p_fixed_crd_rule->GetX();
			if( type == Y_ELEM_CRD ) return p_fixed_crd_rule->GetY();
			if( type == Z_ELEM_CRD ) return p_fixed_crd_rule->GetZ();
			throw std::runtime_error(" incompatible ElemCrd type for FixedCrdRule ");
		}
		Pt3CrdRule* p_pt3_rule = dynamic_cast< Pt3CrdRule* > (p_rule);
		if( p_pt3_rule != NULL )
		{
			if( type == LEN_ELEM_CRD ) return p_pt3_rule->GetBondLen();
			if( type == ANG_ELEM_CRD ) return p_pt3_rule->GetValAng();
			if( type == DIH_ELEM_CRD ) return p_pt3_rule->GetDihAng();
			throw std::runtime_error(" incompatible ElemCrd type for Pt3CrdRule ");
		}
		throw std::runtime_error(" unsupported coordinate assign rule ");
	}
	catch( const std::exception& ex )
	{
		PrintLog("Error in ElemCrd::GetValue() \n");
		PrintLog("%s\n",ex.what() );
	}
	return 0.0;
}

void ElemCrd::SetDisplayValue( double val )
{
	if( type == ANG_ELEM_CRD || type == DIH_ELEM_CRD )
	{
		SetValue( DEG_TO_RAD * val);
		return;
	}
	SetValue( val );
}

double ElemCrd::GetDisplayValue() const 
{
	if( type == ANG_ELEM_CRD || type == DIH_ELEM_CRD ) return RAD_TO_DEG * GetValue();
	return GetValue();
}


HaAtom* ElemCrd::GetManagedAtom()
{
	SingleAtomCrdRule* p_sngl_rule  = dynamic_cast<SingleAtomCrdRule*>(p_rule);
	if( p_sngl_rule != NULL) return p_sngl_rule->GetManagedAtom();
	return NULL;
}

CrdAssignRule* ElemCrd::GetCrdAssignRule()
{
	return p_rule;
}
	
void ElemCrd::SetTag( const std::string& tag_new)
{
	tag = tag_new;
}
	
std::string ElemCrd::GetTag() const
{
	return tag;
}

void ElemCrd::SetFrozen( bool set_flag )
{
	frozen_flag = set_flag;
}

bool ElemCrd::IsFrozen() const
{
	return frozen_flag;
}

ZMatLoadOptions::ZMatLoadOptions()
{
	SetStdOptions();
}

ZMatLoadOptions::ZMatLoadOptions( const ZMatLoadOptions& ref ) 
{ 
	SetStdOptions(); 
	Copy(ref); 
}
	
ZMatLoadOptions::~ZMatLoadOptions()
{

}

void ZMatLoadOptions::SetStdOptions()
{
	SetLoadAtID(false);
	SetLoadAtElem(true);
}

void ZMatLoadOptions::Copy( const harlem::HashMap& ref ) 
{  
	const ZMatLoadOptions* pref = dynamic_cast<const ZMatLoadOptions*>(&ref);
	if( pref != NULL )
	{
		SetLoadAtID( pref->ToLoadAtID() );
		SetLoadAtElem( pref->ToLoadAtElem() );
	}
}

harlem::HashMap* ZMatLoadOptions::clone() const
{ 
	return new ZMatLoadOptions(*this); 
} 

void ZMatLoadOptions::SetLoadAtID(bool set_flag )
{
	load_at_id = set_flag;
}
	
bool ZMatLoadOptions::ToLoadAtID() const
{
	return load_at_id; 
}

void ZMatLoadOptions::SetLoadAtElem(bool set_flag)
{
	load_at_elem = set_flag;
}
	
bool ZMatLoadOptions::ToLoadAtElem() const
{
	return load_at_elem;
}

ZMatSaveOptions::ZMatSaveOptions()
{
	SetStdOptions();
}

ZMatSaveOptions::ZMatSaveOptions( const ZMatSaveOptions& ref ) 
{ 
	SetStdOptions(); 
	Copy(ref); 
}

ZMatSaveOptions::~ZMatSaveOptions()
{

}

void ZMatSaveOptions::SetStdOptions()
{
	SetSaveAtSeqNum(false);
	SetSaveAtElem(false);
	SetSaveAtSymbol(true);
	SetSaveTags(true);
}

void ZMatSaveOptions::Copy( const harlem::HashMap& ref ) 
{  
	const ZMatSaveOptions* pref = dynamic_cast<const ZMatSaveOptions*>(&ref);
	if( pref != NULL )
	{
		SetSaveAtSeqNum( pref->ToSaveAtSeqNum() );
		SetSaveAtElem( pref->ToSaveAtElem() );
		SetSaveAtSymbol( pref->ToSaveAtSymbol() ); 
		SetSaveTags( pref->ToSaveTags() );
	}
}

harlem::HashMap* ZMatSaveOptions::clone() const
{ 
	return new ZMatSaveOptions(*this); 
} 

void ZMatSaveOptions::SetSaveAtSeqNum(bool set_flag )
{
	save_at_seq_num = set_flag;
}
	
bool ZMatSaveOptions::ToSaveAtSeqNum() const
{
	return save_at_seq_num;
}

void ZMatSaveOptions::SetSaveAtElem(bool set_flag )
{
	save_at_elem = set_flag;
	if( set_flag ) save_at_symbol = false;
}
	
bool ZMatSaveOptions::ToSaveAtElem() const
{
	return save_at_elem;
}

void ZMatSaveOptions::SetSaveAtSymbol(bool set_flag )
{
	save_at_symbol = set_flag;
	if( set_flag ) save_at_elem   = false;
}
	
bool ZMatSaveOptions::ToSaveAtSymbol() const
{
	return save_at_symbol;
}

void ZMatSaveOptions::SetSaveTags(bool set_flag )
{
	save_tags = set_flag;	
}

bool ZMatSaveOptions::ToSaveTags() const
{
	return save_tags;
}

ZMatCrd::ZMatCrd(MolSet* pmset_new )
{
	Clear();
	pmset = pmset_new;
}

ZMatCrd::~ZMatCrd()
{
	Clear();
}

void ZMatCrd::Clear()
{
	size_t i,n;
	n = elem_crds.size();
	for( i = 0; i < n; i++)
	{
		delete elem_crds[i];
	}
	elem_crds.clear();

	n = rules.size();
	for( i = 0; i < n; i++)
	{
		delete rules[i];
	}
	rules.clear();

	n = dummy_atoms.size();
	for( i = 0; i < n; i++ )
	{
		if( dummy_atoms[i]->GetHostRes() == NULL ) delete dummy_atoms[i];
	}
	dummy_atoms.clear();
}

bool ZMatCrd::IsEmpty() const
{
	if( rules.empty() && elem_crds.empty() ) return true;
	return false;
}


void ZMatCrd::InitStdZMat()
{
	Clear();

	int na = pmset->GetNAtoms();
	elem_crds.reserve(3*na);
	rules.reserve(na);

	int i = 0;
	HaAtom* aptr;
	AtomIteratorMolSet aitr(pmset);

	HaAtom* aptr1 = NULL;
	HaAtom* aptr2 = NULL;
	HaAtom* aptr3 = NULL;

	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom(), i++ )
	{
		if( i == 0 )
		{
			rules.push_back( new FixedCrdRule(aptr) );
			FixedCrdRule* p_rule = (FixedCrdRule*) rules.back();
			p_rule->SetX(0.0);
			p_rule->SetY(0.0);
			p_rule->SetZ(0.0);
			aptr1 = aptr;
			continue;
		}
		if( i == 1 )
		{
			rules.push_back( new Pt3CrdRule(aptr, aptr1, NULL, NULL) );
			Pt3CrdRule* p_rule = (Pt3CrdRule*) rules.back();
			p_rule->SetValAng(0.0);
			p_rule->SetDihAng(0.0);
			elem_crds.push_back( new ElemCrd( LEN_ELEM_CRD, rules.back() ) );
			aptr2 = aptr1;
			aptr1 = aptr;
			continue;
		}
		if( i == 2 )
		{
			rules.push_back( new Pt3CrdRule(aptr, aptr1, aptr2, NULL) );
			Pt3CrdRule* p_rule = (Pt3CrdRule*) rules.back();
			p_rule->SetDihAng(0.0);
			elem_crds.push_back( new ElemCrd( LEN_ELEM_CRD, rules.back() ) );
			elem_crds.push_back( new ElemCrd( ANG_ELEM_CRD, rules.back() ) );
			aptr3 = aptr2;
			aptr2 = aptr1;
			aptr1 = aptr;
			continue;
		}
		rules.push_back( new Pt3CrdRule(aptr, aptr1, aptr2, aptr3) );
		Pt3CrdRule* p_rule = (Pt3CrdRule*) rules.back();
		elem_crds.push_back( new ElemCrd( LEN_ELEM_CRD, rules.back() ) );
		elem_crds.push_back( new ElemCrd( ANG_ELEM_CRD, rules.back() ) );
		elem_crds.push_back( new ElemCrd( DIH_ELEM_CRD, rules.back() ) );
		aptr3 = aptr2;
		aptr2 = aptr1;
		aptr1 = aptr;
	}
}

void ZMatCrd::InitAllXYZ()
{
	Clear();

	HaAtom* aptr;
	AtomIteratorMolSet aitr(pmset);
	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
	{
		rules.push_back(new FixedCrdRule(aptr) ); 
		elem_crds.push_back(new ElemCrd(X_ELEM_CRD, rules.back()) );
		elem_crds.push_back(new ElemCrd(Y_ELEM_CRD, rules.back()) );
		elem_crds.push_back(new ElemCrd(Z_ELEM_CRD, rules.back()) );
	}
}

int ZMatCrd::OnDelAtoms( AtomContainer& del_atoms_par )
{
	try
	{
		AtomIteratorGen aitr_del(&del_atoms_par );
		HaAtom* aptr;

		std::set<HaAtom*> del_atoms;

		for(aptr = aitr_del.GetFirstAtom(); aptr; aptr = aitr_del.GetNextAtom())
		{
			del_atoms.insert(aptr);
		}

		std::map< HaAtom*,std::set<SingleAtomCrdRule*> >  dep_tree; // dependency tree showing coordinate assign rules in Z matrix that depend on a given atom
		std::map< HaAtom*,std::list<HaAtom*> > dep_tree_rev;    // reverse dependency tree showing atoms that a given atom depends on 
		std::set<SingleAtomCrdRule*> empty_rule_set;
		std::list<HaAtom*> empty_atlist;
		std::set<CrdAssignRule*> del_rules;

		std::vector< SingleAtomCrdRule* >::iterator ritr;
		for( ritr = rules.begin(); ritr != rules.end(); ritr++ )
		{
			SingleAtomCrdRule* pr = (*ritr);
			HaAtom* aptr = pr->GetManagedAtom();

			dep_tree[aptr] = empty_rule_set;
			dep_tree_rev[aptr] = empty_atlist;

			Pt3CrdRule* p3r = dynamic_cast<Pt3CrdRule*>(pr);
			if( p3r != NULL )
			{
				HaAtom* aptr1 = dynamic_cast<HaAtom*>(p3r->GetRefPt1());
				HaAtom* aptr2 = dynamic_cast<HaAtom*>(p3r->GetRefPt2());
				HaAtom* aptr3 = dynamic_cast<HaAtom*>(p3r->GetRefPt3());

				dep_tree_rev[aptr].push_back(aptr1);
				dep_tree_rev[aptr].push_back(aptr2);
				dep_tree_rev[aptr].push_back(aptr3);
				
				if( aptr1 != NULL )  
				{
					if( dep_tree.count(aptr1) == 0 ) dep_tree[aptr1] = empty_rule_set;
					dep_tree[aptr1].insert(pr);
				}
				if( aptr2 != NULL )  
				{
					if( dep_tree.count(aptr2) == 0 ) dep_tree[aptr2] = empty_rule_set;
					dep_tree[aptr2].insert(pr);
				}
				if( aptr3 != NULL )  
				{
					if( dep_tree.count(aptr3) == 0 ) dep_tree[aptr3] = empty_rule_set;
					dep_tree[aptr3].insert(pr);
				}
			}
			Pt2CrdRule* p2r = dynamic_cast<Pt2CrdRule*>(pr);
			if( p2r != NULL )
			{
				HaAtom* aptr1 = p2r->GetRefAtom1();
				HaAtom* aptr2 = p2r->GetRefAtom2();

				dep_tree_rev[aptr].push_back(aptr1);
				dep_tree_rev[aptr].push_back(aptr2);

				if( aptr1 != NULL )  
				{
					if( dep_tree.count(aptr1) == 0 ) dep_tree[aptr1] = empty_rule_set;
					dep_tree[aptr1].insert(pr);
				}
				if( aptr2 != NULL )  
				{
					if( dep_tree.count(aptr2) == 0 ) dep_tree[aptr2] = empty_rule_set;
					dep_tree[aptr2].insert(pr);
				}
			}
			SameAtomCrdRule* psr = dynamic_cast<SameAtomCrdRule*>(pr);
			{
				 if( psr != NULL )
				 {
					 HaAtom* aptr1 = psr->GetRefAtom();

					 dep_tree_rev[aptr].push_back(aptr1);

					 if( dep_tree.count(aptr1) == 0 ) dep_tree[aptr1] = empty_rule_set;
					 dep_tree[aptr1].insert(pr);
				 }
			}
		}
		bool break_cycle = true;
		for(;;) // delete atoms from Z-matrix that are in the list of deleted atoms and no other atoms depend on them (leaves)
		{
			for( ritr = rules.begin(); ritr != rules.end();  )
			{
				SingleAtomCrdRule* pr = (*ritr);
				HaAtom* aptr = pr->GetManagedAtom();
				if( (del_atoms.count(aptr) > 0) && dep_tree[aptr].empty() )
				{
					std::list<HaAtom*>::iterator sitr =  dep_tree_rev[aptr].begin();
					for( ; sitr != dep_tree_rev[aptr].end(); sitr++)
					{
						HaAtom* aptr_dep = (*sitr);
						if( dep_tree.count( aptr_dep) > 0 )
						{
							dep_tree[aptr_dep].erase( pr );
						}
					}

					del_rules.insert(pr);
					delete pr;
					ritr = rules.erase(ritr);
					break_cycle = false;
				}
				else
				{
					ritr++;
				}
			}
			if( break_cycle ) break;
			break_cycle = true;
		}

		for( ritr = rules.begin(); ritr != rules.end(); ritr++ )
		{
			SingleAtomCrdRule* pr = (*ritr);
			HaAtom* aptr = pr->GetManagedAtom();
			if( del_atoms.count(aptr) > 0 )
			{
				HaAtom* aptr_dummy = AddDummyAtom();
				aptr_dummy->SetCoordFrom( *aptr );
				std::set<SingleAtomCrdRule*>::iterator ritr_dep = dep_tree[aptr].begin();
				for( ; ritr_dep != dep_tree[aptr].end(); ritr_dep++ )
				{
					SingleAtomCrdRule* pr_dep = (*ritr_dep);
					pr_dep->ReplaceRefPt( aptr, aptr_dummy);
				}
				pr->SetManagedAtom( aptr_dummy );
			}
		}

		std::vector< ElemCrd* >::iterator el_itr = elem_crds.begin();
		for( ; el_itr != elem_crds.end();  )
		{
			ElemCrd* pcrd = (*el_itr);
			if( del_rules.count( pcrd->GetCrdRule() ) > 0 )
			{
				el_itr =  elem_crds.erase(el_itr);
			}
			else
			{
				el_itr++;
			}
		}
	}
	catch( const std::exception& ex )
	{
		PrintLog(" ZMatCrd::OnDelAtoms() \n");
		PrintLog(" %s\n",ex.what());
		return FALSE;
	}
	return TRUE;

}

int ZMatCrd::LoadFromString( const std::string& str, const harlem::HashMap* popt_par  )
{
	std::istringstream is(str);
	return LoadFromStream(is, popt_par);
}

int ZMatCrd::LoadFromStream(std::istream& is, const harlem::HashMap* popt_par )
{
	this->Clear();

	const ZMatLoadOptions* popt = &opt_read_default;
	if( popt_par != NULL ) popt = dynamic_cast<const ZMatLoadOptions*>(popt_par);

	AtomIteratorMolSet aitr(pmset);
	int i = -1;
	std::string line;

	std::string clbl[3] = {"BOND LEN","VALENCE ANGLE","DIHEDRAL ANGLE"};

	typedef std::multimap<std::string, ElemCrd* >::value_type TagCrdPair;
	std::multimap<std::string, ElemCrd*> tag_crd_map;

	std::vector<HaAtom*> atoms_zmat;

	try
	{
		for(;;)
		{
			std::getline(is,line);
			boost::trim(line);
			if( is.eof() ) break;
			if( i < 0 && line.empty() ) continue;

			if( line.empty() ) break;

			i++;

			std::vector<std::string> str_vec; 
			boost::split(str_vec,line,boost::is_any_of(" ,"),boost::token_compress_on);

			int j = 0;
			HaAtom* aptr_v[4] = { NULL, NULL, NULL, NULL };
		
			double blen = -1000.0; 
			double ang  = -1000.0;
			double dih  = -1000.0;

			std::string blen_tag;
			std::string ang_tag;
			std::string dih_tag;

			if( popt->ToLoadAtID() )
			{
				if( harlem::IsInt(str_vec[j]) ) 
				{
					int num_seq = atoi( str_vec[j].c_str() );
					if( num_seq != (i+1) ) throw std::runtime_error(" Out of sequence number in Atom Sequence number field "); 
					aptr_v[0] = aitr.next();
				}
				else
				{
					std::string at_ref = str_vec[j];
					if( at_ref[0] == 'X' && ( at_ref.size() == 1 || isdigit(at_ref[1]) )  )
					{	
						aptr_v[0] = AddDummyAtom( at_ref );
					}
					else
					{
						aptr_v[0] = pmset->GetAtomByRef( at_ref.c_str() );
						if( aptr_v[0] == NULL ) throw std::runtime_error( " Invalid Atom reference to Atom 0 ");
					}
				}
				j++;
			}

			if( popt->ToLoadAtElem() )
			{
				if( j >= str_vec.size() ) throw std::runtime_error(" No atom element information ");
				std::string elem_str = str_vec[j];
					
				int elemno = 0;
				if( elem_str[0] == 'X' && ( elem_str.size() == 1 || isdigit(elem_str[1]) )  )
				{
					if( aptr_v[0] == NULL ) aptr_v[0] = AddDummyAtom( elem_str );
				}
				else
				{
					if( harlem::IsInt(elem_str) )
					{
						elemno = atoi(elem_str.c_str());
					}
					else
					{
						elemno = HaAtom::GetElemNoFromName(elem_str);
					}
					if( aptr_v[0] == NULL ) aptr_v[0] = aitr.next();
					if( elemno != aptr_v[0]->GetElemNo() ) throw std::runtime_error(" Element name or number is incorrect ");
				}
				j++;
			}

			if( aptr_v[0] == NULL ) aptr_v[0] = aitr.next();

			atoms_zmat.push_back( aptr_v[0] );

			if( i == 0 ) 
			{
				rules.push_back( new FixedCrdRule(aptr_v[0]) );
				FixedCrdRule* p_rule = (FixedCrdRule*) rules.back();
				p_rule->SetX(0.0);
				p_rule->SetY(0.0);
				p_rule->SetZ(0.0);
				continue;
			}
			
			int kmax = 4;
			if( i == 1 ) kmax = 2;
			if( i == 2 ) kmax = 3;

			for(int k = 1; k < kmax; k++)
			{
				if( j >= str_vec.size() ) 
				{
					throw std::runtime_error(" No ID information for reference atom: " + harlem::ToString(k) );
				}
				std::string at_ref = str_vec[j];
				if( harlem::IsInt(at_ref) )
				{
					int seq_no = atoi(at_ref.c_str()) - 1;
					if( seq_no >= atoms_zmat.size() || seq_no < 0 ) throw std::runtime_error(" Invalid stom sequence number " + at_ref );
					aptr_v[k] = atoms_zmat[seq_no];
				}
				else
				{
					aptr_v[k] = pmset->GetAtomByRef(at_ref.c_str());
				}
				if( aptr_v[k] == NULL )  throw std::runtime_error(" Invalid reference atom " + harlem::ToString(k)  );
				j++;

				if( j >= str_vec.size() ) throw std::runtime_error(" No value for " + clbl[k-1]);
				
				std::string tag = str_vec[j];
			
				if( harlem::IsFloat(tag) )
				{
					if( k == 1)      blen = atof( tag.c_str() );
					else if( k == 2) ang  = atof( tag.c_str() )*DEG_TO_RAD;
					else if( k == 3) dih  = atof( tag.c_str() )*DEG_TO_RAD;
				}
				else
				{
					if( k == 1 ) blen_tag     = tag;
					else if( k == 2 ) ang_tag = tag;
					else if( k == 3 ) dih_tag = tag;
				}
				j++;
			}

			rules.push_back( new Pt3CrdRule(aptr_v[0], aptr_v[1], aptr_v[2], aptr_v[3] ) );
			Pt3CrdRule* p_rule = (Pt3CrdRule*) rules.back();

			if( aptr_v[1] != NULL )
			{
				elem_crds.push_back( new ElemCrd( LEN_ELEM_CRD, rules.back() ) );
				ElemCrd* p_el_crd = elem_crds.back();
				if( blen > 0.0 )
				{
					p_el_crd->SetValue(blen);
				}
				else
				{
					p_el_crd->SetTag(blen_tag);
					tag_crd_map.insert( TagCrdPair(blen_tag,p_el_crd) );
				}
			}

			if( aptr_v[2] != NULL )
			{
				elem_crds.push_back( new ElemCrd( ANG_ELEM_CRD, rules.back() ) );
				ElemCrd* p_el_crd = elem_crds.back();
				if( ang > -999.0 )
				{
					p_el_crd->SetValue(ang);
				}
				else
				{
					p_el_crd->SetTag(ang_tag);
					tag_crd_map.insert( TagCrdPair(ang_tag,p_el_crd) );
				}
			}
			if( aptr_v[3] != NULL )
			{
				elem_crds.push_back( new ElemCrd( DIH_ELEM_CRD, rules.back() ) );
				ElemCrd* p_el_crd = elem_crds.back();
				if( dih > -999.0 )
				{
					p_el_crd->SetValue(dih);
				}
				else
				{
					p_el_crd->SetTag(dih_tag);
					tag_crd_map.insert( TagCrdPair(dih_tag,p_el_crd) );
				}
			}
		} // end for(;;) on atoms

		bool frozen_crd = false;
		for(;;) // read tagged elemental coordinate values
		{
			std::getline(is,line);
			boost::trim(line);
			if( is.eof() ) break;
			if( !frozen_crd && line.empty() ) 
			{
				frozen_crd = true;
				continue;
			}
			if( line.empty() ) break;

			std::vector<std::string> str_vec;
			boost::split(str_vec,line,boost::is_any_of(" ,"),boost::token_compress_on);

			if( str_vec.size() != 2) throw std::runtime_error(" Expect two tokens on the line to assign value of the tagged coordinate ");
			std::string tag = str_vec[0];
			if( tag_crd_map.count(tag) == 0 ) throw std::runtime_error(" Unknown coorinate tag " + tag );

			if( !harlem::IsFloat(str_vec[1])) throw std::runtime_error(" Value of the coordinate is not a float number ");

			double val = atof(str_vec[1].c_str());

			std::multimap<std::string, ElemCrd*>::iterator titr1 = tag_crd_map.lower_bound(tag);
			std::multimap<std::string, ElemCrd*>::iterator titr2 = tag_crd_map.upper_bound(tag);
			std::multimap<std::string, ElemCrd*>::iterator titr;
			for( titr = titr1; titr != titr2; titr++ )
			{
				ElemCrd* p_el_crd = (*titr).second;
				p_el_crd->SetDisplayValue( val );
				if( frozen_crd ) p_el_crd->SetFrozen(true);
			}
		}
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in ZMatCrd::LoadFromStream() \n ");
		PrintLog(" Reading line: \n%s\n",line.c_str() );
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return TRUE;
}

int ZMatCrd::SetAtomCrd()
{
	try
	{
		size_t i;
		size_t n = rules.size();
		for( i = 0; i < n; i++ )
		{
			int ires = rules[i]->SetManagedAtomCrd();
			if( !ires ) 
			{
				std::stringstream ss; ss << i; 
				throw std::runtime_error(" Fail to set managed atom coordinates for rule # "  + ss.str() ) ;
			}
		}
		return TRUE;
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in ZMatCrd::SetAtomCrd() \n");
		PrintLog("%s\n",ex.what());
	}
	return FALSE;
}

int ZMatCrd::GetCrdSnapshot( HaVec_double& crd_arr )
{
	SetAtomCrd();
	crd_arr.resize( pmset->GetNAtoms() * 3 );
	HaAtom* aptr;
	int i = 0;
	AtomIteratorMolSet aitr(pmset);
	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom() )
	{
		crd_arr.SetVal_idx0(i,aptr->GetX()); i++;
		crd_arr.SetVal_idx0(i,aptr->GetY()); i++;
		crd_arr.SetVal_idx0(i,aptr->GetZ()); i++;
	}
	return TRUE;
}

int ZMatCrd::SetFromAtomCrd()
{
	try
	{
		if( IsEmpty() ) throw std::runtime_error(" Z-Matrix is not defined ");
		size_t i;
		size_t n = rules.size();
		for( i = 0; i < n; i++ )
		{
			FixedCrdRule* p_fix_rule = dynamic_cast< FixedCrdRule*>(rules[i]);
			Pt3CrdRule*   p_pt3_rule = dynamic_cast< Pt3CrdRule*>  (rules[i]);
			if ( p_fix_rule != NULL )
			{
				p_fix_rule->SetParFromCurPos();
				continue;
			}
			if( p_pt3_rule != NULL )
			{
				p_pt3_rule->SetParFromCurPos();
			}
		}
		return TRUE;
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in ZMatCrd::SetFromAtomCrd() \n");
		PrintLog("%s\n",ex.what());
	}
	return FALSE;
}

int ZMatCrd::SetFromCrdSnapshot( const HaVec_double& crd_arr )
{
	int na = pmset->GetNAtoms(); 
	if( crd_arr.size() != na )
	{
		PrintLog(" Error in ZMatCrd::SetFromCrdSnapshot() \n");
		PrintLog(" size of coordinate array = %d  not equal to 3* natoms = %d \n", crd_arr.size(), na );
		return FALSE;
	}
	int ires = pmset->SetCrdFromArray( crd_arr );
	if( !ires ) return FALSE;
	return SetFromAtomCrd();
}

int ZMatCrd::SaveToStream(std::ostream& os, const harlem::HashMap* popt_par )
{
//	PrintLog(" ZMatCrd::SaveToStream() pt 1 \n");
	const ZMatSaveOptions* popt = dynamic_cast<const ZMatSaveOptions*>(popt_par);
		
	std::auto_ptr<ZMatSaveOptions> popt_auto( popt == NULL ? (ZMatSaveOptions*) opt_save_default.clone(): (ZMatSaveOptions*) popt->clone() );
	popt = popt_auto.get();

	char buf[256];
	try
	{
		AtomIntMap at_sn_map;
		size_t nc = elem_crds.size();
		size_t ic;
		
		std::multimap<std::string, ElemCrd* > tagged_elem_crd;	
		typedef std::multimap<std::string, ElemCrd* >::iterator   TagCrdIterator;
		for(ic = 0; ic < nc; ic++ )
		{
			std::string tag = elem_crds[ic]->GetTag();
			if( !tag.empty() ) tagged_elem_crd.insert( std::pair<std::string,ElemCrd*>(tag,elem_crds[ic]));
		}
		ic = 0;

		size_t nr = rules.size();
		size_t ir;

		for( ir = 0; ir < nr; ir++ )
		{
			HaAtom* aptr = rules[ir]->GetManagedAtom();
			at_sn_map[aptr] = ir;

			if( popt->ToSaveAtSeqNum() ) 
			{
				sprintf(buf," %4d  ",(ir+1)); os << buf;
			}

			if( popt->ToSaveAtSymbol() )
			{
				std::string std_smbl = aptr->GetStdSymbol();
				if( aptr->IsDummy() ) std_smbl = "X";
				sprintf(buf," %3s  ",std_smbl.c_str()); os << buf;
			}
			else if( popt->ToSaveAtElem() )
			{
				if( aptr->IsDummy() )
				{
					os << "  X   ";
				}
				else
				{
					sprintf(buf," %3d  ",aptr->GetElemNo()); os << buf;
				}
			}
			
			FixedCrdRule* p_fix_rule = dynamic_cast< FixedCrdRule*>(rules[ir]);
			Pt3CrdRule*   p_pt3_rule = dynamic_cast< Pt3CrdRule*>  (rules[ir]);
			if ( p_fix_rule != NULL )
			{
				ElemCrd* pc_x = NULL;
				ElemCrd* pc_y = NULL;
				ElemCrd* pc_z = NULL;
				for( ; ic < nc; ic++ )
				{	
					if( elem_crds[ic]->GetManagedAtom() != aptr ) break; 
					if( elem_crds[ic]->GetType() == X_ELEM_CRD ) pc_x = elem_crds[ic];
					if( elem_crds[ic]->GetType() == Y_ELEM_CRD ) pc_y = elem_crds[ic];
					if( elem_crds[ic]->GetType() == Z_ELEM_CRD ) pc_z = elem_crds[ic];
				}

				if( pc_x != NULL && pc_y != NULL && pc_z != NULL )
				{
					if( popt->ToSaveTags() && !pc_x->GetTag().empty())
					{
						sprintf(buf,"%12s  ",pc_x->GetTag().c_str()); os << buf;
					}
					else
					{
						sprintf(buf,"%12.6f  ",pc_x->GetValue()); os << buf;
					}

					if( popt->ToSaveTags() && !pc_y->GetTag().empty())
					{
						sprintf(buf,"%12s  ",pc_y->GetTag().c_str()); os << buf;
					}
					else
					{
						sprintf(buf,"%12.6f  ",pc_y->GetValue()); os << buf;
					}

					if( popt->ToSaveTags() && !pc_z->GetTag().empty())
					{
						sprintf(buf,"%12s  ",pc_z->GetTag().c_str()); os << buf;
					}
					else
					{
						sprintf(buf,"%12.6f  ",pc_z->GetValue()); os << buf;
					}
				}
			}

			if( p_pt3_rule != NULL )
			{
				ElemCrd* pc_len = NULL;
				ElemCrd* pc_ang = NULL;
				ElemCrd* pc_dih = NULL;
				for( ; ic < nc; ic++ )
				{	
					if( elem_crds[ic]->GetManagedAtom() != aptr ) break; 
					if( elem_crds[ic]->GetType() == LEN_ELEM_CRD ) pc_len = elem_crds[ic];
					if( elem_crds[ic]->GetType() == ANG_ELEM_CRD ) pc_ang = elem_crds[ic];
					if( elem_crds[ic]->GetType() == DIH_ELEM_CRD ) pc_dih = elem_crds[ic];
				}

				AtomIntMap::iterator mitr;

				if( pc_len != NULL )
				{
					HaAtom* aptr1 = dynamic_cast<HaAtom*>(p_pt3_rule->GetRefPt1());
					mitr = at_sn_map.find(aptr1);
					if( mitr == at_sn_map.end() ) throw std::runtime_error(" Can not find Atom ref 1 in coordinate assign rule ");
					int ns = (*mitr).second;
					sprintf(buf,"%5i ",(ns + 1)); os << buf;
					if( !popt->ToSaveTags() || pc_len->GetTag().empty() )
					{
						sprintf(buf,"%12.6f   ",pc_len->GetDisplayValue()); os << buf;
					}
					else
					{
						sprintf(buf,"%12s   ",pc_len->GetTag().c_str()); os << buf;
					}
				}
				if( pc_ang != NULL )
				{
					HaAtom* aptr2 = dynamic_cast<HaAtom*>(p_pt3_rule->GetRefPt2());
					mitr = at_sn_map.find(aptr2);
					if( mitr == at_sn_map.end() ) throw std::runtime_error(" Can not find Atom ref 2 in coordinate assign rule ");
					int ns = (*mitr).second;
					sprintf(buf,"%5i ",(ns + 1)); os << buf;
					if( !popt->ToSaveTags() || pc_ang->GetTag().empty() )
					{
						sprintf(buf,"%12.6f   ",pc_ang->GetDisplayValue()); os << buf; 
					}
					else
					{
						sprintf(buf,"%12s   ",pc_ang->GetTag().c_str()); os << buf;
					}
				}
				if( pc_dih != NULL )
				{
					HaAtom* aptr3 = dynamic_cast<HaAtom*>(p_pt3_rule->GetRefPt3());
					mitr = at_sn_map.find(aptr3);
					if( mitr == at_sn_map.end() ) throw std::runtime_error(" Can not find Atom ref 2 in coordinate assign rule ");
					int ns = (*mitr).second;
					sprintf(buf,"%5i ",(ns + 1)); os << buf;
					if( !popt->ToSaveTags() || pc_dih->GetTag().empty() )
					{
						sprintf(buf,"%12.6f   ",pc_dih->GetDisplayValue()); os << buf; 
					}
					else
					{
						sprintf(buf,"%12s   ",pc_dih->GetTag().c_str()); os << buf;
					}
				}				
			}
			os << std::endl;
		}
		if( popt->ToSaveTags() && !tagged_elem_crd.empty() )
		{
			os << "   " << std::endl;
			std::map<std::string,double> frozen_val_map;
			TagCrdIterator titr1 = tagged_elem_crd.begin();
			while( titr1 != tagged_elem_crd.end() )
			{
				std::string tag = (*titr1).first;
				TagCrdIterator titr2 = tagged_elem_crd.upper_bound(tag);
				
				TagCrdIterator titr = titr1;
				int nc = 0;
				double avg_val = 0.0;
				double last_val = 0.0;
				double cur_val  = 0.0;
				bool frozen = false;
				for( ; titr != titr2; ++titr )
				{
					nc++;
					last_val = cur_val;
					cur_val = (*titr).second->GetDisplayValue();
					if( (*titr).second->IsFrozen() ) frozen = true;
					if( nc > 1 ) 
					{
						if( fabs(last_val - cur_val) > 0.000001 ) 
						{
							PrintLog(" some values of coordinates with the tag %s are different\n", tag.c_str() );
							PrintLog(" They will be averaged on output \n");
						}
					}
					avg_val += cur_val;
				}
				avg_val = avg_val/nc;
				if( frozen ) 
				{
					frozen_val_map[tag] = avg_val;
				}
				else
				{
					sprintf(buf,"%12s %12.6f \n",tag.c_str(),avg_val); os << buf;
				}
				titr1 = titr2;
			}
			if( !frozen_val_map.empty() )
			{
				os << "   " << std::endl;
				std::map<std::string,double>::iterator fitr = frozen_val_map.begin();
				for( ; fitr != frozen_val_map.end(); ++fitr )
				{
					sprintf(buf,"%12s %12.6f \n",(*fitr).first.c_str(),(*fitr).second); os << buf;
				}
			}
		}
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in ZMatCrd::SaveToStream() \n");
		PrintLog(" %s\n",ex.what());
		return FALSE;
	}
	return TRUE;
}

std::string ZMatCrd::SaveToString( const harlem::HashMap* popt_par )
{
	std::ostringstream os;
	SaveToStream(os,popt_par);
	return os.str();
}

int ZMatCrd::SaveXMLToStream(std::ostream& os, const harlem::HashMap* popt_arg )
{
	if( os.fail() ) return FALSE;

	const harlem::SaveOptions* popt = dynamic_cast<const harlem::SaveOptions*>(popt_arg);

	std::auto_ptr<harlem::SaveOptions>popt_auto( popt == NULL ? new harlem::SaveOptions() : (harlem::SaveOptions*) popt->clone() );
	popt = popt_auto.get();

	bool save_header = popt->ToSaveHeader();
	bool save_footer = popt->ToSaveFooter();

	if( save_header )
	{
		os << harlem::StdXMLHeader()     << std::endl;
		os << harlem::HarlemDataHeader() << std::endl;
	}

	os << "<zmat mset=\"" << pmset->GetName() << "\" nat=\"" << pmset->GetNAtoms() << "\" nz=\"" << this->rules.size() <<"\" >" << std::endl;
	SaveToStream(os,popt);
	os << "</zmat>" << std::endl;

	if( save_footer )
	{
		os << harlem::HarlemDataFooter() << std::endl;
	}
	return TRUE;
}

int ZMatCrd::LoadXMLNode( rapidxml::xml_node<>* node_zmat, const harlem::HashMap* popt_arg )
{
	using namespace rapidxml;

	try
	{
		std::string tag = node_zmat->name();
		boost::to_lower(tag);
		if( tag != "zmat" ) throw std::runtime_error(" Name of the node is not zmat ");

		const ZMatLoadOptions* popt = dynamic_cast<const ZMatLoadOptions*>(popt_arg);

		std::auto_ptr<ZMatLoadOptions> popt_auto( popt == NULL ? new ZMatLoadOptions() : (ZMatLoadOptions*) popt->clone() );
		popt = popt_auto.get();

		xml_attribute<>* attr;
		for( attr = node_zmat->first_attribute(); attr; attr = attr->next_attribute() )
		{  
			std::string tag_attr = attr->name();
			if( boost::iequals(tag_attr, "mset_name") || boost::iequals(tag_attr, "mset" ) )
			{
				std::string mset_name = attr->value();
				std::string mset_name_curr = pmset->GetName();
				if( !boost::iequals(mset_name, mset_name_curr ) ) 
				{
					PrintLog("Warning in ZMatCrd::LoadXMLNode() \n");
					PrintLog(" Molecular set name of the Z-matrix %s  is not equal to the name of the current molset %s \n",
						     mset_name.c_str(), mset_name_curr.c_str());
						       
				}
			}
			else if( boost::iequals(tag_attr, "nat") )
			{
				int nat = atoi( attr->value() );
				if( pmset->GetNAtoms() != nat ) throw std::runtime_error(" number of atoms in Z-matrix " + harlem::ToString(nat) + "not equal to that of Molecular Set " + harlem::ToString(pmset->GetNAtoms()) );
			}
		}
		
		std::string str_zmat = node_zmat->value();
		int ires = LoadFromString(str_zmat,popt);
		if(!ires ) return FALSE;
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in ZMatCrd::LoadXMLNode() \n");
		PrintLog(" %s\n",ex.what());
		return FALSE;
	}
	return TRUE;
}

int ZMatCrd::GetNZ() const
{
	return rules.size();
}

int ZMatCrd::GetNCrd() const
{
	return elem_crds.size();
}

int ZMatCrd::GetNCrdUnFrozen() const
{
	int nc = elem_crds.size();
	int ic;
	int nc_unf = 0;
	for( ic = 0; ic < nc; ic++ )
	{
		if( !elem_crds[ic]->IsFrozen() ) nc_unf++; 
	}
	return nc_unf;
}

void ZMatCrd::FreezeCrdAll()
{
	int n = elem_crds.size();
	int i;
	for( i = 0; i < n; i++ )
	{
		elem_crds[i]->SetFrozen( true );
	}
}

void ZMatCrd::UnFreezeCrdAll()
{
	int n = elem_crds.size();
	int i;
	for( i = 0; i < n; i++ )
	{
		elem_crds[i]->SetFrozen( false );
	}
}



int ZMatCrd::GetElemCrdVal( HaVec_double& elem_crd_val_arr , bool unfrozen ) const 
{
	int ncrd;
	if( unfrozen )
	{
		ncrd = GetNCrdUnFrozen();
	}
	else
	{
		ncrd = GetNCrd();
	}

	elem_crd_val_arr.resize( ncrd );
	int nc = elem_crds.size();
	int ic;
	int ic_unf = 0;
	for( ic = 0; ic < nc; ic++ )
	{
		if( !unfrozen || !elem_crds[ic]->IsFrozen() ) 
		{
			elem_crd_val_arr[ic_unf] = elem_crds[ic]->GetValue();
			ic_unf++;
		}
	}
	return TRUE;
}

int ZMatCrd::SetFromElemCrdVal( const HaVec_double& elem_crd_val_arr , bool unfrozen )
{
	int ncrd;
	if( unfrozen )
	{
		ncrd = GetNCrdUnFrozen();
	}
	else
	{
		ncrd = GetNCrd();
	}

	if( elem_crd_val_arr.size() != ncrd ) 
	{
		PrintLog(" Error in ZMatCrd::SetFromIntCrdVal() \n");
		PrintLog(" Size of internal(elemental) coordinates values array = %d   is not equal to the number of coordinates = %d \n", 
			       elem_crd_val_arr.size(), ncrd );
		return FALSE;
	}

	int nc = elem_crds.size();
	int ic;
	int ic_unf = 0;
	for( ic = 0; ic < nc; ic++ )
	{
		if( !unfrozen || !elem_crds[ic]->IsFrozen() ) 
		{
			elem_crds[ic]->SetValue( elem_crd_val_arr[ic_unf] );
			ic_unf++;
		}
	}
	return TRUE;
}

extern "C"
{
	void str_gs_ ( int* ic, int* iat0, int* iat1, double* b, int* ib, double* crd );
    void bend_gs_( int* ic, int* iat0, int* iat1, int* iat2, double* b, int* ib, double* crd );
    void tors_gs_( int* ic, int* iat0, int* iat1, int* iat2, int* iat3, double* b, int* ib, double* crd );
	void formg_gs_( int* nt, int* ib  , double* b, double* g);
	void tranf_gs_( int* nc, double* fx , double* f, int* ib, double* b, double* g);
}


int ZMatCrd::TransDerivToIntCrd( const HaVec_double& deriv_cart_crd, HaVec_double& deriv_int_crd, bool unfrozen )
{
	// PrintLog(" ZMatCrd::TransDerivToIntCrd() pt 1 \n");
	int na = pmset->GetNAtoms();
	int ncart = 3*na;
	
	if( deriv_cart_crd.size() != ncart )
	{
		PrintLog(" Error in ZMatCrd::TransDerivToIntCrd() \n" );
		PrintLog(" Dimension of array of cartesian derivatives = %d  is not equalt to the number of cart coords = %d \n",
		          deriv_cart_crd.size(), ncart );
		return FALSE;
	}
	
	int ncrd_tot = GetNCrd();
	int ncrd = ncrd_tot;
	if( unfrozen ) ncrd = GetNCrdUnFrozen();

	int ires = CalcBGMatr();
	if( !ires ) return FALSE;
	
	deriv_int_crd.resize(ncrd);

	HaVec_double deriv_int_loc(ncrd_tot);
	HaVec_double deriv_cart_loc( deriv_cart_crd );
	
	tranf_gs_( &ncrd_tot, deriv_cart_loc.v() , deriv_int_loc.v(), ib_mat.v(), b_mat.v(), g_mat.v() );

	if( ncrd_tot == ncrd ) 
	{
		deriv_int_crd = deriv_int_loc;
	}
	else
	{
		int i;
		int ic_unf = 0;
		for( i = 0; i < ncrd_tot; i++ )
		{
			if( !elem_crds[i]->IsFrozen() )
			{
				deriv_int_crd[ic_unf] = deriv_int_loc[i];
				ic_unf++;
			}
		}
	}
	return TRUE;
}

int ZMatCrd::CalcBGMatr()
{	
//	PrintLog(" ZMatCrd::CalcBGMatr() pt 1 \n");
	int ncrd = GetNCrd();
	int na = pmset->GetNAtoms();

	ib_mat.resize(4,ncrd); ib_mat = 0;
	b_mat.resize(12,ncrd); b_mat = 0.0;

	AtomIntMap at_seqn_map = pmset->GetAtomSeqNumMap();

	SetAtomCrd();

	HaMat_double atcrd(3,na);
	AtomIteratorMolSet aitr(pmset);
	HaAtom* aptr;
	int ia = 0;
	for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		ia++;
		atcrd.SetVal_idx1(1,ia, aptr->GetX() );
		atcrd.SetVal_idx1(2,ia, aptr->GetY() );
		atcrd.SetVal_idx1(3,ia, aptr->GetZ() );
	}
	
	try
	{
		int ic;

		for( ic = 1; ic <= ncrd; ic++ )
		{
			ElemCrd* pcrd = elem_crds[ic-1];
			HaAtom*  aptr_m = pcrd->GetManagedAtom();
			if( aptr_m == NULL ) throw std::runtime_error(" Elemental coordinate is not associated with atom ");
			if( at_seqn_map.count(aptr_m) == 0 ) throw std::runtime_error(" Has not found atom of elemental coordinate in atom - Seq Num map ");
			int iat0 = at_seqn_map[aptr_m] + 1;

			if( pcrd->GetType() ==  X_ELEM_CRD || pcrd->GetType() == Y_ELEM_CRD ||  pcrd->GetType() == Z_ELEM_CRD )
			{
				ib_mat.SetVal_idx1( 1,ic,iat0 );
				int j = 1;
				if( pcrd->GetType() ==  X_ELEM_CRD ) j = 1;
				if( pcrd->GetType() ==  Y_ELEM_CRD ) j = 2;
				if( pcrd->GetType() ==  Z_ELEM_CRD ) j = 3;
				b_mat.SetVal_idx1(j, ic, 1.0 );
				continue;
			}

			if( pcrd->GetType() == LEN_ELEM_CRD || pcrd->GetType() == ANG_ELEM_CRD || pcrd->GetType() == DIH_ELEM_CRD)
			{
				Pt3CrdRule* p_rule = dynamic_cast<Pt3CrdRule*>( pcrd->GetCrdRule() );
				if( p_rule == NULL ) throw std::runtime_error(" Bond length Elemental coordinate is not associated with 3pt coord set rule ");

				HaAtom* aptr1 = dynamic_cast<HaAtom*>(p_rule->GetRefPt1());
				if( at_seqn_map.count(aptr1) == 0 ) throw std::runtime_error(" Can not find Atom ref 1 in coordinate assign rule ");		
				int iat1 = at_seqn_map[aptr1] + 1;

				if( pcrd->GetType() == LEN_ELEM_CRD ) 
				{
        			str_gs_(&ic,&iat0, &iat1, b_mat.v(), ib_mat.v(), atcrd.v() );
					continue;
				}

				HaAtom* aptr2 = dynamic_cast<HaAtom*>(p_rule->GetRefPt2());
				if( at_seqn_map.count(aptr2) == 0 ) throw std::runtime_error(" Can not find Atom ref 2 in coordinate assign rule ");
				int iat2 = at_seqn_map[aptr2] + 1;

				if( pcrd->GetType() == ANG_ELEM_CRD ) 
				{
					bend_gs_(&ic,&iat0,&iat1,&iat2,b_mat.v(), ib_mat.v(), atcrd.v() );
					continue;
				}

				HaAtom* aptr3 = dynamic_cast<HaAtom*>(p_rule->GetRefPt3());
				if( at_seqn_map.count(aptr3) == 0 ) throw std::runtime_error(" Can not find Atom ref 3 in coordinate assign rule ");
				int iat3 = at_seqn_map[aptr3] + 1;

				if( pcrd->GetType() == DIH_ELEM_CRD ) 
				{
					tors_gs_( &ic,&iat0,&iat1,&iat2,&iat3,b_mat.v(),ib_mat.v(),atcrd.v() );
					continue;
				}
			}
		}
		  
		g_mat.resize(ncrd,ncrd); g_mat = 0.0;    
		formg_gs_( &ncrd,ib_mat.v(),b_mat.v(),g_mat.v() );

		int ires = HaMat_double::mat_inverse( g_mat ); 
		if( !ires ) throw runtime_error(" Failed to invert G-matrix");

	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in ZMatCrd::CalcBGMatr() \n");
		PrintLog(" %s\n",ex.what());
		ib_mat.clear();
		b_mat.clear();
		g_mat.clear();
		return FALSE;
	}
	return TRUE;
}

bool ZMatCrd::SetXYZCrd( int iat, double x, double y, double z )
{
	if( iat < 0 || iat >= rules.size() ) return false;
	HaAtom* aptr = rules[iat]->GetManagedAtom();
	if( aptr == NULL ) return false;
	aptr->SetX(x);
	aptr->SetY(y);
	aptr->SetZ(z);
	return true;
}


void ZMatCrd::SetCrdDesc(int      iat, int iat_r,      int iat_ang,      int iat_dih)
{

}

void ZMatCrd::SetCrdDesc(HaAtom* aptr, HaAtom* aptr_r, HaAtom* aptr_ang, HaAtom* aptr_dih)
{

}
	
void ZMatCrd::SetCrdDesc(const std::string& at_nm, const std::string& at_r_nm, const std::string& at_a_nm, 
		                 const std::string& at_d_nm)
{

}

ElemCrd* ZMatCrd::GetRCrdForAtom  ( HaAtom* aptr )
{
	size_t nc = elem_crds.size();
	size_t i;
	for( i = 0; i < nc; i++ )
	{
		ElemCrd* pc = elem_crds[i];
		if( pc->GetManagedAtom() == aptr && pc->GetType() == LEN_ELEM_CRD ) return pc;
	}
	PrintLog("Error in ZMatCrd::GetRCrdForAtom() \n");
	PrintLog(" No R coordinate for atom %s \n", aptr->GetRef().c_str());
	return NULL;
}

ElemCrd* ZMatCrd::GetAngCrdForAtom( HaAtom* aptr )
{
	size_t nc = elem_crds.size();
	size_t i;
	for( i = 0; i < nc; i++ )
	{
		ElemCrd* pc = elem_crds[i];
		if( pc->GetManagedAtom() == aptr && pc->GetType() == ANG_ELEM_CRD ) return pc;
	}
	PrintLog("Error in ZMatCrd::GetAngCrdForAtom() \n");
	PrintLog(" No Valence Angle coordinate for atom %s \n", aptr->GetRef().c_str());
	return NULL;
}
	
ElemCrd* ZMatCrd::GetDihCrdForAtom( HaAtom* aptr )
{
	size_t nc = elem_crds.size();
	size_t i;
	for( i = 0; i < nc; i++ )
	{
		ElemCrd* pc = elem_crds[i];
		if( pc->GetManagedAtom() == aptr && pc->GetType() == DIH_ELEM_CRD ) return pc;
	}
	PrintLog("Error in ZMatCrd::GetDihCrdForAtom() \n");
	PrintLog(" No Dihedral Angle coordinate for atom %s \n", aptr->GetRef().c_str());
	return NULL;
}

ElemCrd* ZMatCrd::GetCrdByTag( const std::string& tag)
{
	int nc = elem_crds.size();
	int i;
	for( i = 0; i < nc; i++ )
	{
		if( boost::iequals(elem_crds[i]->GetTag(),tag) ) return elem_crds[i];
	}
	return NULL;
}

void ZMatCrd::SetTagRCrd  ( HaAtom* aptr, const std::string& tag )
{
	ElemCrd* pel_crd = GetRCrdForAtom(aptr);
	if( pel_crd != NULL ) pel_crd->SetTag(tag);
}

void ZMatCrd::SetTagAngCrd( HaAtom* aptr, const std::string& tag )
{
	ElemCrd* pel_crd = GetAngCrdForAtom(aptr);
	if( pel_crd != NULL ) pel_crd->SetTag(tag);
}

void ZMatCrd::SetTagDihCrd( HaAtom* aptr, const std::string& tag )
{
	ElemCrd* pel_crd = GetDihCrdForAtom(aptr);
	if( pel_crd != NULL ) pel_crd->SetTag(tag);
}

	
void ZMatCrd::SetRVal(HaAtom* aptr, double val, CoordUnits units )
{
	ElemCrd* pcrd = GetRCrdForAtom(aptr);
	if( pcrd == NULL) return;
	if( units == ANGSTROM_U )
	{
		pcrd->SetValue(val);
	}
	else if( units == BOHR_U )
	{
		pcrd->SetValue( val * BOHR_TO_ANG );
	}
}
	
void ZMatCrd::SetRVal(int iat, double val, CoordUnits units )
{
	HaAtom* aptr = pmset->GetAtomBySeqNum( iat );
	SetRVal( aptr, val, units );
}
	
void ZMatCrd::SetRVal(const std::string& at_nm, double val, CoordUnits units )
{
	HaAtom* aptr = pmset->GetAtomByRef( at_nm.c_str() );
	SetRVal( aptr, val, units );
}

void ZMatCrd::SetAngVal(HaAtom* aptr, double val, AngleUnits units )
{
	ElemCrd* pcrd = GetAngCrdForAtom(aptr);
	if( pcrd == NULL) return;
	if( units == DEGREE_U )
	{
		pcrd->SetValue(val * DEG_TO_RAD );
	}
	else if( units == RADIAN_U )
	{
		pcrd->SetValue( val );
	}
}

void ZMatCrd::SetAngVal(int iat, double val, AngleUnits units )
{
	HaAtom* aptr = pmset->GetAtomBySeqNum( iat );
	SetAngVal( aptr, val, units );
}

void ZMatCrd::SetAngVal(const std::string& at_nm, double val, AngleUnits units )
{
	HaAtom* aptr = pmset->GetAtomByRef( at_nm.c_str() );
	SetAngVal( aptr, val, units );
}

void ZMatCrd::SetDihVal(HaAtom* aptr, double val, AngleUnits units )
{
	ElemCrd* pcrd = GetDihCrdForAtom(aptr);
	if( pcrd == NULL) return;
	if( units == DEGREE_U )
	{
		pcrd->SetValue(val * DEG_TO_RAD );
	}
	else if( units == RADIAN_U )
	{
		pcrd->SetValue( val );
	}
}

void ZMatCrd::SetDihVal(int iat, double val, AngleUnits units )
{
	HaAtom* aptr = pmset->GetAtomBySeqNum( iat );
	SetDihVal( aptr, val, units );
}

void ZMatCrd::SetDihVal(const std::string& at_nm, double val, AngleUnits units )
{
	HaAtom* aptr = pmset->GetAtomByRef( at_nm.c_str() );
	SetDihVal( aptr, val, units );
}

HaAtom* ZMatCrd::AddDummyAtom(const std::string& at_name )
{
	HaAtom* aptr = new HaAtom();
	aptr->SetName(at_name);
	dummy_atoms.push_back(aptr);
	return aptr;
}