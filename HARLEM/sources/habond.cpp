/*! \file habond.cpp
 
    Functions of classes 
    to define Chemical Bond object in HARLEM.
   
    \author Igor Kurnikov  
    \date 1997-2002 
  
*/

#define HABOND_CPP

#include <assert.h>
#include <memory>

#include "haatom.h"
#include "habond.h"
#include "hamolecule.h"
#include "hamolset.h"
#include "haresdb.h"


HaBond::HaBond()
{
    srcatom= NULL;
    dstatom= NULL;	

    radius = 0.0;
    irad=0;
    col = 0;
	flag = 0;
	type = SINGLE_BOND;
}

HaBond::HaBond(const HaBond& ref)
{
	srcatom=ref.srcatom;
	dstatom=ref.dstatom;
	this->SetParamFrom(ref);
}


HaBond::HaBond(HaAtom* src, HaAtom* dst )
{
    assign(src,dst);
}

HaBond::~HaBond()
{

}

HaHBond* MolSet::AddHBond(HaAtom* src, HaAtom* dst)
{
	   MolSet* pmset = src->GetHostMolSet();

	   if(pmset != this || dst->GetHostMolSet() != this ) return NULL;
 
	   std::set<HaHBond>::iterator hitr;

	   HaHBond hbond(src,dst);
	
	   hitr = pmset->HBonds.find(hbond);
	   if(hitr != pmset->HBonds.end())
		   return( (HaHBond*) &(*hitr));

	   HaHBond hbond2(dst,src);
	
	   hitr = pmset->HBonds.find(hbond2);
	   if(hitr != pmset->HBonds.end())
		   return( (HaHBond*) &(*hitr));

	   std::pair<std::set<HaHBond>::iterator,bool> pr= pmset->HBonds.insert(hbond);
	   if(pr.second)
	   {
			return( (HaHBond*) &(*pr.first) );
	   }
	   return(NULL);
}

void MolSet::CreateHydrogenBond(HaAtom* src, HaAtom* dst, int energy, int offset )
{
    HaHBond* ptr;
    int flag;
	
	if(src == NULL || dst == NULL)
	{
		PrintLog(" Error In MolSet::CreateHydrogenBond() \n");
		PrintLog(" NULL Atom Pointers \n");
		return;
	}

	HaResidue* src_res = src->GetHostRes();
	HaResidue* dst_res = dst->GetHostRes();

	HaAtom* srcCA = NULL;
	HaAtom* dstCA = NULL;

	if( src_res->IsAmino() )
	{
		srcCA = src_res->GetAtomByName("CA");
		dstCA = dst_res->GetAtomByName("CA");
	}
	else if( src_res->IsNucleo())
	{
		srcCA = src_res->GetAtomByName("P");
		dstCA = dst_res->GetAtomByName("P");
	}

    ptr = AddHBond(src,dst);
    if( !ptr ) 
	{
		ErrorInMod("HaMolset::CreateHydrogenBond()", 
	               "Failed to add new Hydrogen Bond");
		return;
	}
	
    if( (offset>=-128) && (offset<127) )
    {   
		ptr->offset = (char)offset;
    } 
	else 
		ptr->offset = 0;
	
	flag = 1;
    flag = (src->Selected() && dst->Selected());

	if(flag) ptr->Select();
	
    ptr->srcCA = srcCA;
    ptr->dstCA = dstCA;
    ptr->energy = energy;
    ptr->col = 0;
}


bool HaBond::SetParamFrom(const HaBond& ref)
{
	radius = ref.radius;
	irad   = ref.irad;
	col    = ref.col;
	flag   = ref.flag;
	type   = ref.type;
	return true;
}

bool HaBond::SetTypeFrom(const HaBond& ref)
{
	type   = ref.type;
	if( ref.IsSingle() ) this->SetSingle();
	if( ref.IsDouble() ) this->SetDouble();
	if( ref.IsTriple() ) this->SetTriple();
	if( ref.IsAromatic() ) this->SetAromatic();
	return true;
}

std::string HaBond::GetTypeString() const
{
	if( this->IsSingle() ) return "SINGLE";
	if( this->IsDouble() ) return "DOUBLE";
	if( this->IsTriple() ) return "TRIPLE";
	if( this->IsAromatic() ) return "AROMATIC";
	if( this->IsVirtual() )  return "VIRTUAL";
	return "UNKNOWN";
}

HaBond* HaBond::assign(HaAtom* src, HaAtom* dst )
{
	type = SINGLE_BOND;
	flag = NormBondFlag;
	Select();
	// atom with a smaller pointer is always the first 
	if( dst < src)
	{
		srcatom = dst;
		dstatom = src;
	}
	else
	{
		srcatom = src;
		dstatom = dst;
	}

    radius = 0.0;
    col = 0; 
	return(this);
}

bool HaBond::operator==(const HaBond & rhs) const
{
   if( srcatom == rhs.srcatom && dstatom == rhs.dstatom )
		return true;
   return false;
}

bool HaBond::operator < (const HaBond & rhs) const
{
   if( srcatom < rhs.srcatom )
	  return true;
   if( srcatom == rhs.srcatom ) 
   {
	  if( dstatom < rhs.dstatom )
	     return true;
   }
   return false;   
}

// Flag Manipulation Utilities:

void HaBond::Select() { flag |= SelectFlag; }

void HaBond::UnSelect()  { flag &= ~SelectFlag; }

int HaBond::Selected() const 
{ 
	return (flag & SelectFlag); 
}

void HaBond::SetNotDraw()
{
	flag &= ~DrawBondFlag;
}

void HaBond::DrawWire()
{
	SetNotDraw();
	flag |= WireFlag;
}

void HaBond::DrawDashed()
{
	SetNotDraw();
	flag |= DashFlag;
}

void HaBond::DrawCylinder(double new_rad)
{
	SetNotDraw();
	flag |= CylinderFlag;
	radius = new_rad;
//  irad = (int)(Scale*rad);
}

int HaBond::IsToDraw() const 
{ 
	return (flag & DrawBondFlag); 
}

HaHBond::HaHBond()
{           
	src = NULL;
	dst = NULL;
	p_hatom = NULL; 
	            
	srcCA=NULL;            
	dstCA=NULL;             
	energy=0;               
	radius=0.0;               
	irad=0;                 
	offset=0;                 
	flag=0;                 
	col=0;                   
}



HaHBond::HaHBond(HaAtom* new_src, HaAtom* new_dst, HaAtom* new_p_hatom)
{           
	src = new_src;
	dst = new_dst;
	p_hatom = new_p_hatom;
	            
	srcCA=NULL;            
	dstCA=NULL;             
	energy=0;               
	radius=0.0;               
	irad=0;                 
	offset=0;                 
	flag=0;                 
	col=0;                   
}

HaHBond::HaHBond(const HaHBond& ref)
{
	dst = ref.dst;
	src = ref.src;
	p_hatom = ref.p_hatom;
	srcCA=ref.srcCA;            
	dstCA=ref.dstCA;             
	energy= ref.energy;               
	radius= ref.radius;               
	irad= ref.irad;                 
	offset= ref.offset;                 
	flag= ref.flag;                 
	col= ref.col;                   
}


HaHBond::~HaHBond()
{

}

bool HaHBond::operator==(const HaHBond & rhs) const
{
   if( src == rhs.src && dst == rhs.dst )
	  return true;
   return false;
}

bool 
HaHBond::operator < (const HaHBond & rhs) const
{
   if( src < rhs.src )
	  return true;
   if( src ==  rhs.src)
   {
	  if( dst < rhs.dst )
	    return true;
   }
   return false;   
}

// Flag Manipulation Utilities:

void HaHBond::Select() { flag |= SelectFlag; }

void HaHBond::UnSelect()  { flag &= ~SelectFlag; }

int HaHBond::Selected() const 
{ 
	return (flag & SelectFlag); 
}

void HaHBond::SetNotDraw()
{
	flag &= ~DrawBondFlag;
}

void HaHBond::SetDraw()
{
	flag |= DrawBondFlag;
}

void HaHBond::DrawWire()
{
	flag &= ~DashFlag;
	flag &= ~CylinderFlag;
	flag |= WireFlag;
}

void HaHBond::DrawDashed()
{
	flag &= ~WireFlag;
	flag &= ~CylinderFlag;
	flag |= DashFlag;
}

void HaHBond::DrawCylinder(double new_rad)
{
	flag &= ~WireFlag;
	flag &= ~DashFlag;
	flag |= CylinderFlag;
	radius = new_rad;
}

int HaHBond::IsToDrawCylinder() const
{
	return (flag & CylinderFlag);
}

int HaHBond::IsToDrawWire() const
{
	return (flag & WireFlag);
}

int HaHBond::IsToDrawDashed() const
{
	return (flag & DashFlag);
}

int HaHBond::IsToDraw() const 
{ 
	return (flag & DrawBondFlag); 
}

int HaHBond::GetHCoord(Vec3D& pt)
{
	pt[0] = 0.0;
	pt[1] = 0.0;
	pt[2] = 0.0;

	HaAtom* p_don = this->GetDonorAtom();
	HaAtom* p_acc = this->GetAcceptorAtom();
	if(p_don == NULL || p_acc == NULL) return FALSE;
	
	AtomGroup bonded_atoms;

	p_don->GetBondedAtoms(bonded_atoms);
	int nb_at = bonded_atoms.size();
	int i;
	double dist2_min = 9999999.0;
	for(i = 0; i < nb_at; i++)
	{
		HaAtom* aptr_h = bonded_atoms[i];
		if( aptr_h->GetElemNo() != 1) continue;

		double dist2 = Vec3D::CalcDistanceSq(aptr_h,p_acc);
		if(dist2 < dist2_min)
		{
			dist2_min = dist2;
			pt.SetCoordFrom(*aptr_h);
		}
	}
	if( dist2_min > 1000000.0 ) return FALSE;
	return TRUE;
}