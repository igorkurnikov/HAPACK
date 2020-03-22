/*! \file transfor.cpp

   Molecule Display classes
   
   \author Igor Kurnikov
   \date 1998-2002
 
   Derived from:

   transfor.c
   RasMol2 Molecular Graphics
   Roger Sayle, August 1995
   Version 2.6
*/

#include "haconst.h"

#include <stdio.h>
#include <math.h>

#define TRANSFORM

#include <boost/algorithm/string.hpp>

#include "abstree.h"
#include "command.h"
#include "hamolview.h"

#include "hasurface.h"

#include "haatom.h"
#include "habond.h"
#include "hamolecule.h"
#include "harlemapp.h"
#include "math_num.h"
#include "moleditor.h"



void HaMolView::SetAtomScreenRadVal( double rad )
{
    int irad,change;
    HaAtom  *aptr;
		
    irad = (int)(Scale*rad);
    DrawAtoms = False;
    change = False;
	
	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
        if( aptr->Selected() )
        {   
            aptr->SetDrawSphere(true);
			aptr->image_radius = rad;
            aptr->irad = irad;
            change = True;
        } 
		else if( aptr->IsDrawSphere() )
        {   
			DrawAtoms = True;
        }
	}
	
	if( change )
	{   
		DrawAtoms = True;
	}
}

void HaMolView::SetRadiusTemperature()
{
    int irad,change;
	double rad;
    HaAtom  *aptr;

    DrawAtoms = False;
    change = False;

    AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
        if( (aptr->Selected() ) && (aptr->tempf > 0) )
        {   
			rad = (aptr->tempf);
            if( rad > 3.0 ) rad = 3.0;

            irad = (int)(Scale*rad);
            aptr->SetDrawSphere(true);
            aptr->image_radius = rad;
            aptr->irad = irad;
            change = True;
        } 
		else if( aptr->IsDrawSphere() )
        {   
			DrawAtoms = True;
        }
	}

    if( change )
    {   
		DrawAtoms = True;
    }
}


void HaMolView::SetAtomScreenRadVdW()
{
    int change;
    HaAtom *aptr;
 
    DrawAtoms = False;
    change = False;

	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
        if( aptr->Selected() )
        {   
            aptr->image_radius = HaAtom::ElemVDWRadius(aptr->GetElemNo());
            aptr->irad = (int)(Scale*aptr->image_radius);
            change = True;

            aptr->SetDrawSphere(true);
        } 
		else if( aptr->IsDrawSphere() )
        {   
			DrawAtoms = True;
        }
	}

    if( change )
    {   
		DrawAtoms = True;
    }
}


void HaMolView::DisableSpacefill()
{
    if( !DrawAtoms ) return;

    DrawAtoms = False;

	HaAtom  *aptr;

    AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
        if( !(aptr->Selected() ) )
        {   
			if( aptr->IsDrawSphere() )
            {   
                DrawAtoms = True;
            }
        } 
		else if( aptr->IsDrawSphere() )
		{
            aptr->SetDrawSphere(false);
		}
	}
}



void HaMolView::EnableWireframe( int mask, double rad )
{
    HaBond  *bptr;
    int flag, irad;
	
    DrawBonds = False;
    irad = (int)(Scale*rad);

	BondIteratorMolSet bitr(GetMolSet());
    for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond())
    {   
		flag = ZoneBoth ? bptr->dstatom->Selected() && bptr->srcatom->Selected()
			            : bptr->dstatom->Selected() || bptr->srcatom->Selected();
		
        if( flag )
        {   
			DrawBonds = True;
            bptr->SetNotDraw();
            bptr->flag |= mask;
            if( mask == CylinderFlag )
            {   
                bptr->radius = rad;
                bptr->irad = irad;
            }
        } 
		else if( bptr->IsToDraw())
        {    
			DrawBonds = True;
        }
    }
}


void HaMolView::DisableWireframe()
{
    HaBond  *bptr;
    int flag;

    if( !DrawBonds ) return;

    DrawBonds = False;

	BondIteratorMolSet bitr(GetMolSet());
    for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond())
    {   
		flag = ZoneBoth ? bptr->dstatom->Selected() && bptr->srcatom->Selected()
			            : bptr->dstatom->Selected() || bptr->srcatom->Selected();

        if( flag )
        {   
			bptr->SetNotDraw();
        } 
		else if( bptr->IsToDraw() )
        {   
			DrawBonds = True;
        }
    }
}


void HaMolView::EnableBackbone( int mask, double rad )
{
    HaChain  *chain;
    HaBond  *bptr;
    int flag,irad;
		
    irad = (int)(Scale*rad);

	MolSet* pmset = GetMolSet();
	pmset->GetMolEditor()->UpdateBackBone(pmset);
	
	int nb = pmset->BackboneBonds.size();
	int ib;
    for( ib = 0; ib < nb; ib++ )
    {   
		bptr = pmset->BackboneBonds[ib];
		flag = ZoneBoth ? bptr->dstatom->Selected() && bptr->srcatom->Selected()
			            : bptr->dstatom->Selected() || bptr->srcatom->Selected();
		
        if( flag )
        {   
			bptr->SetNotDraw();
            bptr->flag |= mask;
            if( mask == CylinderFlag )
            {   
				bptr->radius = rad;
                bptr->irad = irad;
            }
        } 
    }
}

void HaMolView::DisableBackbone()
{
    HaChain  *chain;
    HaBond  *bptr;

	MolSet* pmset = GetMolSet();
	pmset->GetMolEditor()->UpdateBackBone(pmset);

	int nb = pmset->BackboneBonds.size();
	int ib;

    if( ZoneBoth )
    {   
		for( ib = 0; ib < nb; ib++)
		{
			bptr = pmset->BackboneBonds[ib];
            if( bptr->dstatom->Selected() && bptr->srcatom->Selected() )
                bptr->SetNotDraw();
		}
    } 
	else 
	{
		for( ib = 0; ib < nb; ib++)
		{
			bptr = pmset->BackboneBonds[ib];
			if( bptr->dstatom->Selected() || bptr->srcatom->Selected() )
				bptr->SetNotDraw();
		}
	}
}


void HaMolView::SetHBondStatus( int enable, double rad )
{
    set<HaHBond, less<HaHBond> >::iterator bitr;
    HaAtom  *src;
    HaAtom  *dst;
    int flag, irad;
    HaHBond* phb;

	MolSet* pmset = GetMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);

	if(enable)
	{
		p_mol_editor->CalcHBonds(pmset);
	}
	
	irad = (int)(Scale*rad);
	for( bitr = pmset->HBonds.begin(); bitr != pmset->HBonds.end(); bitr++ )
	{   
		phb = (HaHBond*) &(*bitr);
		src = phb->src;  dst = phb->dst;
		
		flag = ZoneBoth? src->Selected() && dst->Selected() : src->Selected() || dst->Selected() ;
		
		if( flag ) 
		{   
			phb->SetDraw();
			if( enable )
			{   
				if( rad > 0.01)
				{   
					phb->DrawCylinder(rad);
					phb->irad = irad;
				} 
				else 
				{
					phb->DrawWire();
				}
			}
		}
		else
		{
			phb->SetNotDraw();
		}
	}
}

void HaMolView::SetSSBondStatus( int enable, double rad )
{
    HaAtom  *src;
    HaAtom  *dst;
    int flag, irad;
    
	MolSet* pmset = GetMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);

	if( enable && (!pmset->SSBonds_found) ) p_mol_editor->FindDisulphideBridges(pmset);
		
	irad = (int)(Scale*rad);

	BondIteratorMolSet bitr(pmset);
	HaBond* phb;

	for( phb = bitr.GetFirstBond(); phb; phb = bitr.GetNextBond() )
	{   
		src = phb->srcatom;  dst = phb->dstatom;
		
		flag = ZoneBoth? src->Selected() && dst->Selected() : 
		src->Selected() || dst->Selected() ;
		
		if( flag ) 
		{   
			phb->SetNotDraw();
			if( enable )
			{   
				if( rad > 0.01)
				{   
					phb->flag |= CylinderFlag;
					phb->radius = rad;
					phb->irad = irad;
				} 
				else 
					phb->flag |= WireFlag;
			}
		}
	}
}

void HaMolView::SetRibbonStatus( int enable, int flag, double width )
{
    HaChain  *chain;
    HaResidue  *group;
    HaAtom  *ptr;
	vector<HaAtom*>::iterator aitr;

	MolSet* pmset = this->GetMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);

    /* Ribbons already disabled! */
    if( !enable && !DrawRibbon )
        return;

	vector<HaMolecule*>::iterator mol_itr;
	ForEachMol_VIEW
	{
		if( !(*mol_itr)->IsSecStructFound() )
			p_mol_editor->DetermineSecStructure((*mol_itr),False);
	}

    DrawRibbon = False;
	ResidueIteratorMolSet ritr( GetMolSet() );
    for( group = ritr.GetFirstRes(); group; group = ritr.GetNextRes() )
	{
		if( enable )
		{   
			if( group->flag & DrawKnotFlag )
				DrawRibbon = True;
			
			for( aitr=group->begin(); aitr != group->end(); aitr++ )
			{
				ptr=*aitr;
				if( ptr->IsAlphaCarbon() )
				{   
					if( ptr->Selected() )
					{   
						group->flag &= ~DrawKnotFlag;
						group->flag |= flag;
						if( (width < 0.01) )
						{   
							if( group->struc & (HelixFlag|SheetFlag) )
							{      
								group->width = 1.5;
							} 
							else 
								group->width = 0.4;
						} 
						else 
							group->width = width;
						DrawRibbon = True;
					}
					break;
					
				} 
				else if( ptr->IsSugarPhosphate() )
				{   
					if( ptr->Selected() )
					{   
						group->width = (width > 0.01)? width : 2.88 ;
						group->flag &= ~DrawKnotFlag;
						group->flag |= flag;
						DrawRibbon = True;
					}
					break;
				}
			}	
				
		} 
		else  /* Disable Ribbon */
			if( group->flag & DrawKnotFlag )
			{   
				for( aitr=group->begin(); aitr != group->end(); aitr++ )
				{
					ptr=*aitr;
					if( ptr->IsAlphaCarbon() ||
						ptr->IsSugarPhosphate() )
					{   
						if( ptr->Selected() )
							group->flag &= ~DrawKnotFlag;
						break;
					}
                    if( group->flag & DrawKnotFlag ) 
                        DrawRibbon = True;
				}
			}
	}
}


void 
HaMolView::SetRibbonCartoons()
{
    HaChain  *chain;
    HaResidue  *group;
    HaAtom  *ptr;
	vector<HaAtom*>::iterator aitr;

	MoleculesType::iterator mol_itr;
	MolSet* pmset = this->GetMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);

	ForEachMol_VIEW
	{
		if( !(*mol_itr)->IsSecStructFound()  )
			p_mol_editor->DetermineSecStructure((*mol_itr),False);
	}

    DrawRibbon = False;
	ResidueIteratorMolSet ritr( pmset );
    for( group = ritr.GetFirstRes(); group; group = ritr.GetNextRes() )
	{   
		if( group->flag & DrawKnotFlag )
			DrawRibbon = True;
		
		for( aitr=group->begin(); aitr != group->end(); aitr++ )
		{
			ptr=*aitr;
			if( ptr->IsAlphaCarbon() )
			{   
				if( ptr->Selected() )
				{   
					group->flag &= ~DrawKnotFlag;
					if( group->struc & (HelixFlag|SheetFlag) )
					{   
						group->flag |= CartoonFlag;
						group->width = 1.52;
					} 
					else 
					{   
						group->flag |= TraceFlag;
						group->width = 0.4;
					}
					DrawRibbon = True;
				}
				break;
				
			} 
			else if( ptr->IsSugarPhosphate() )
			{   
				if( ptr->Selected() )
				{   
					group->flag &= ~DrawKnotFlag;
					group->flag |= RibbonFlag;
					group->width = 2.88;
					DrawRibbon = True;
				}
				break;
			}
		}
	}
}


void HaMolView::SetTraceTemperature()
{
    HaChain  *chain;
    HaResidue  *group;
    HaAtom  *aptr;
    int init,flag;
    double min,max;
    double coeff;

	MolSet* pmset = this->GetMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);
	
    flag = 0;
    init = False;
	min = 1000000.0;
	max = -1000000.0;

	AtomIteratorMolSet aitr(pmset);
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if( aptr->IsAlphaCarbon() || aptr->IsSugarPhosphate() )
		{   
			if(aptr->Selected()) flag = 1 ;
			if( init )
			{   
				if( aptr->tempf < min ) 
				{   
					min = aptr->tempf;
				} 
				else if( aptr->tempf > max )
				{
					max = aptr->tempf;
				}
			} 
			else
			{   
				min = max = aptr->tempf;
				init = True;
			}
		}
	}

    /* No groups selected! */
    if( !flag )
        return;

	MoleculesType::iterator mol_itr;

	ForEachMol_VIEW
	{
		if( !(*mol_itr)->IsSecStructFound()  )
			p_mol_editor->DetermineSecStructure((*mol_itr),False);
	}

    if( max != min )
    {   
		coeff = 1.5/(max-min);
    } 
	else 
		coeff = 0.0;

    DrawRibbon = False;
	ResidueIteratorMolSet ritr(this->GetMolSet());

	for( group = ritr.GetFirstRes(); group; group = ritr.GetNextRes() )
	{   
		if( group->flag & DrawKnotFlag )
			DrawRibbon = True;

		AtomIteratorAtomGroup aitr(group);
		for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
		{
			if( aptr->IsAlphaCarbon() || aptr->IsSugarPhosphate() )
			{   
				if( aptr->Selected() )
				{   
					group->width = coeff*(aptr->tempf - min) + 0.1;
					group->flag &= ~DrawKnotFlag;
					group->flag |= TraceFlag;
					DrawRibbon = True;
				}
				break;
			}
		}
	}
}


/*===========================*/
/* Atom Selection Functions! */
/*===========================*/
 

void HaMolView::RestrictSelected()
{
    HaBond  *bptr;
    HaChain  *chain;
    HaResidue  *group;
    HaAtom  *aptr;

    DrawAtoms = False;  
    DrawBonds = False;   
    DrawLabels = False;
	
	MolSet* pmset = GetMolSet();
	pmset->GetMolEditor()->UpdateBackBone(pmset);
	
	AtomIteratorMolSet aitr(pmset);
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if( aptr->Selected() )
		{   
			aptr->SetDisplayed(true);
			
			if( aptr->IsDrawSphere() )
			{   
				DrawAtoms = True;
			}
			if( !aptr->label.empty() )
				DrawLabels = True;
			
		}  
		else 
		{   
			aptr->SetDisplayed(false);
			aptr->SetDrawSphere(false);
			if( !aptr->label.empty() )
			{   
				aptr->label = "";
			}
		}
	}
    
	pmset->DisplaySelectCount();
	
	BondIteratorMolSet bitr(pmset);
    for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond()) 
    {  
       
        if( bptr->dstatom->Selected() && bptr->srcatom->Selected() )
        {   
			bptr->Select();
			DrawBonds = True;
		
            if( !bptr->IsToDraw())
            {   
				bptr->flag |= WireFlag;
            } 
        } 
		else 
		{
			bptr->UnSelect();
			bptr->SetNotDraw();
		}
    }
	
	int nb = pmset->BackboneBonds.size();
	int ib;
    for( ib = 0; ib < nb; ib++ )
    {   /* Ignore ZoneBoth setting! */
		bptr = pmset->BackboneBonds[ib];
        if( !(bptr->dstatom->Selected() && bptr->srcatom->Selected()) )
		{
 			bptr->UnSelect();
			bptr->SetNotDraw();           
		}
    }
	
    if( DrawRibbon )
    {   
		DrawRibbon = False;
		ResidueIteratorMolSet ritr(this->GetMolSet());
		for( group= ritr.GetFirstRes(); group; group= ritr.GetNextRes() )
		{
			if( group->flag & DrawKnotFlag )
			{   
				AtomIteratorAtomGroup aitr(group);
				for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
				{
					if( aptr->IsAlphaCarbon() ||
						aptr->IsSugarPhosphate() )
					{   
						if( !(aptr->Selected()) )
							group->flag &= ~DrawKnotFlag;
						break;
					}
					if( group->flag & DrawKnotFlag )
						DrawRibbon = True;
				}
			}
		}
    }
}

int HaMolView::ComputeRevColourMap()
{
	min_color_map.clear();
	rev_color_map.clear();

	uint_4* ptr = (uint_4*)pCanv->FBuffer;
	uint_4* end = (uint_4*)(pCanv->FBuffer + (uint_4) pCanv->XRange() * pCanv->YRange());

	int lut_val;

//	IntIntMap rev_full_color_map;

//	for(i = 0; i< LutSize; i++)
//	{
//		lut_val = pCanv->Lut[i];
//		if( rev_full_color_map.find(lut_val) == rev_full_color_map.end())
//		{
//			rev_full_color_map[lut_val] = i;
//		}
//	}

	IntIntMap::iterator mitr;

	while( ptr < end)
	{
		lut_val = *ptr;
		if( rev_color_map.find(lut_val) == rev_color_map.end())
		{
			min_color_map.push_back(lut_val);
			int idx = min_color_map.size() - 1;
			rev_color_map[lut_val] = idx;
		}
		ptr++;
	}
	return min_color_map.size();
}


void HaMolView::RefreshColors()
//! This function is called in HaMolViewWX::RefreshScreen
//! when RFColour flag in ReDrawFlag is set
//!
//! Management of colors in HARLEM:
//! Array Shade[] - Contain Basic Colors of HARLEM
//! Array Lut[] 
//! as first goes basic color and then with increasing intensity colors obtained by darkening of the basic color
//! ColourDepth - the number of intermediate (darkened) colors to interpolate effects of shaded parts of the object 
//! Shade2Colour(j) - macro return the position of the basic color specified by shade in palitra Lut
//!
{
	Canvas3D::empty_lut_idx = 0; 
	
//	PrintLog(" HaMolView::RefreshColors() \n");

	vector<ColorVal> cvals = HaColor::used_colors;
	
	HaColor::used_colors.clear();
    HaColor::cval_idx_map.clear();

	int n = cvals.size();
	
	int i;
	for( i = 0; i < n ; i++)
	{
		int r = RComp(cvals[i]);
		int g = GComp(cvals[i]);
		int b = BComp(cvals[i]);
		int idx = HaColor::RegisterColor(r,g,b);
//		PrintLog(" color (%d,%d,%d) registered with cidx = %d\n",r,g,b,idx);
	}		
}


int HaMolView::ColorAtomsByProp( const std::string& str_prop_par, DValColorMap* p_col_map )
{
	std::string str_prop = boost::trim_copy(str_prop_par);
	boost::to_upper(str_prop);

	MolSet* pmset = GetMolSet();
	AtomIteratorMolSet aitr(pmset);
	HaAtom* aptr;

	try
	{
		if( p_col_map == NULL ) throw std::runtime_error(" Color map is not specified ");
		if( str_prop == "TEMPERATURE")
		{
			for( aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom() )
			{
				if( !aptr->Selected() ) continue;
				HaColor* p_color = p_col_map->GetColorForVal(aptr->tempf);
				aptr->col = p_color->cidx;
			}
		}
		UpdateThisView(RFColour);
	}
	catch( const std::exception& ex )
	{
		PrintLog("Error in HaMolView::ColorAtomsByProp() \n");
		PrintLog("%s\n", ex.what());
		return FALSE;
	}
	return TRUE;
}

void HaMolView::ColourBondNone()
{
    register HaBond  *bptr;

	BondIteratorMolSet bitr(GetMolSet());
    for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond()) 
	{
		if( bptr->Selected() && bptr->col )
		{   
			bptr->col = 0;
		}
	}
}

void HaMolView::ColourBondAttrib( int r, int g, int b )
{
    HaBond  *bptr;
	BondIteratorMolSet bitr(GetMolSet());
	
	HaColor bcolor(r,g,b);
	
	for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond()) 
	{
		if( bptr->Selected() )
		{   
			bptr->col = bcolor.cidx;
		}
	}
}

void HaMolView::ColourBackNone()
{
    HaChain  *chain;
    HaBond  *bptr;
    int flag;

	MolSet* pmset = GetMolSet();
	pmset->GetMolEditor()->UpdateBackBone(pmset);

	int nb = pmset->BackboneBonds.size();
	int ib;

	for( ib = 0; ib < nb; ib++ )
	{   
		bptr = pmset->BackboneBonds[ib];
		flag = ZoneBoth ? bptr->dstatom->Selected() && bptr->srcatom->Selected()
			            : bptr->dstatom->Selected() || bptr->srcatom->Selected();
		
		if( flag)
		{   
			bptr->Select();
			if( bptr->col )
			{   
				bptr->col = 0;
			}
		} 
		else 
			bptr->UnSelect();
	}
}


void HaMolView::ColourBackAttrib( int r, int g, int b )
{
    HaChain  *chain;
    HaBond  *bptr;

	ColourBackNone();
	HaColor bcolor(r,g,b);
	
	MolSet* pmset = GetMolSet();
	pmset->GetMolEditor()->UpdateBackBone(pmset);

	int nb = pmset->BackboneBonds.size();
	int ib;

	for( ib = 0; ib < nb; ib++ )
	{
		bptr = pmset->BackboneBonds[ib];
		if( bptr->Selected() )
		{   
			bptr->col = bcolor.cidx;
		}
	}
}


void 
HaMolView::ColourHBondNone()
{
    set<HaHBond, less<HaHBond> >  *list_ptr;
    set<HaHBond, less<HaHBond> >::iterator  ptr;
    HaAtom  *src;
    HaAtom  *dst;
    HaHBond* phb;

	MolSet* pmset = GetMolSet();

	if(!pmset->HBonds_found) return;
	list_ptr = &(pmset->HBonds);
	
	if( ZoneBoth )
	{   
		for( ptr=list_ptr->begin(); ptr != list_ptr->end(); ptr++ )
		{   
			phb = (HaHBond*)&(*ptr);
			src = phb->src;  dst = phb->dst;
			
			if( src->Selected() && dst->Selected() )
			{   
				phb->Select();
				if( phb->col )
				{   
					phb->col = 0;
				}
			} 
			else 
				phb->UnSelect();
		}
	} 
	else
		for( ptr=list_ptr->begin(); ptr != list_ptr->end(); ptr++ )
		{   
			phb = (HaHBond*)&(*ptr);
			src = phb->src;  dst = phb->dst;
			
			if( src->Selected() && dst->Selected() )
			{   
				phb->Select();
				if( phb->col )
				{   
					phb->col = 0;
				}
			} 
			else 
				phb->UnSelect();
		}
}

void 
HaMolView::ColourSSBondNone()
{
    HaAtom  *src;
    HaAtom  *dst;

	MolSet* pmset = GetMolSet();

	if(!pmset->SSBonds_found) return;

	BondIteratorMolSet bitr(pmset);
	HaBond* bptr;
	if( ZoneBoth )
	{   
		for( bptr = bitr.GetFirstBond() ; bptr; bptr = bitr.GetNextBond() )
		{   
			if( bptr->srcatom->GetElemNo() != 32 || bptr->dstatom->GetElemNo() != 32 ) continue;

			src = bptr->srcatom;  dst = bptr->dstatom;
			
			if( src->Selected() && dst->Selected() )
			{   
				bptr->Select();
				if( bptr->col )
				{   
					bptr->col = 0;
				}
			} 
			else 
				bptr->UnSelect();
		}
	} 
	else
	{
		for( bptr = bitr.GetFirstBond() ; bptr; bptr = bitr.GetNextBond() )
		{   
			if( bptr->srcatom->GetElemNo() != 32 || bptr->dstatom->GetElemNo() != 32 ) continue;
			src = bptr->srcatom;  dst = bptr->dstatom;
			
			if( src->Selected() && dst->Selected() )
			{   
				bptr->Select();
				if( bptr->col )
				{   
					bptr->col = 0;
				}
			} 
			else 
				bptr->UnSelect();
		}
	}
}


void HaMolView::ColourHBondType()
{
    set<HaHBond, less<HaHBond> >::iterator hptr;
    HaHBond* phb;

	ColourHBondNone();

	MolSet* pmset = GetMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);
	p_mol_editor->CalcHBonds(pmset);

	HaColorMap hbond_cmap;
	
	hbond_cmap.AddColor(255, 255, 255); /* 0  Offset =  2   */
	hbond_cmap.AddColor(255,   0, 255); /* 1  Offset =  3   */
	hbond_cmap.AddColor(255,   0,   0); /* 2  Offset =  4   */
	hbond_cmap.AddColor(255, 165,   0); /* 3  Offset =  5   */
	hbond_cmap.AddColor(  0, 255, 255); /* 4  Offset = -3   */
	hbond_cmap.AddColor( 0, 255,   0 ); /* 5  Offset = -4   */
	hbond_cmap.AddColor(255, 255,  0);  /* 6  Others        */
	
	HaColor* pcol;
	
	for( hptr= pmset->HBonds.begin(); hptr != pmset->HBonds.end() ; hptr++ )
	{
		phb = (HaHBond*) &(*hptr);
		if( (*hptr).Selected() )
		{   
			switch( (*hptr).offset )
			{   
			case(  2 ):  pcol = &hbond_cmap.GetColorByIdx(0);  break;
			case(  3 ):  pcol = &hbond_cmap.GetColorByIdx(1);  break;
			case(  4 ):  pcol = &hbond_cmap.GetColorByIdx(2);  break;
			case(  5 ):  pcol = &hbond_cmap.GetColorByIdx(3);  break;
			case( -3 ):  pcol = &hbond_cmap.GetColorByIdx(4);  break;
			case( -4 ):  pcol = &hbond_cmap.GetColorByIdx(5);  break;
			default:     pcol = &hbond_cmap.GetColorByIdx(6);  break;
			}
			
			phb->col = pcol->cidx;
		}
	}
}


void HaMolView::ColourHBondAttrib( int r, int g, int b )
{
	ColourHBondNone();

	MolSet* pmset = GetMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);
	
	p_mol_editor->CalcHBonds(pmset);
	
	HaColor color(r,g,b);
	
	set<HaHBond, less<HaHBond> >::iterator  bitr;
	for( bitr = pmset->HBonds.begin(); bitr != pmset->HBonds.end(); bitr++)
	{
		HaHBond* phb = (HaHBond*)&(*bitr);
		if( phb->Selected())
		{   
			phb->col = color.cidx;
		}
	}
}

void 
HaMolView::ColourSSBondAttrib( int r, int g, int b )
{
	ColourSSBondNone();

	MolSet* pmset = GetMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);

	if( !pmset->SSBonds_found )
	{   
		p_mol_editor->FindDisulphideBridges(pmset);
	} 
	
	HaColor color(r,g,b);

	BondIteratorMolSet bitr(pmset);
	HaBond* bptr;
	for( bptr = bitr.GetFirstBond(); bptr; bptr = bitr.GetNextBond() )
	{
		if( bptr->srcatom->GetElemNo() != 32 || bptr->dstatom->GetElemNo() != 32 ) continue;
		if( bptr->Selected())
		{   
			bptr->col = color.cidx;
		}
	}
}

void 
HaMolView::ColourRibbonNone( int flag )
{
    register HaResidue  *group;
    register HaAtom  *aptr;

	MoleculesType::iterator mol_itr;

	ForEachMol_VIEW
	{
		if(!(*mol_itr)->IsSecStructFound() ) continue;
        
		AtomIteratorMolecule aitr(*mol_itr);
		for(aptr= aitr.GetFirstAtom(); aptr;aptr= aitr.GetNextAtom())
		{
			if( aptr->Selected() && 
				(aptr->IsAlphaCarbon()||
				aptr->IsSugarPhosphate()) )
			{   
				group=aptr->GetHostRes();
				if( (flag&RibColInside) && group->col1 )
				{   
					group->col1 = 0;
				}
				if( (flag&RibColOutside) && group->col2 )
				{   
					group->col2 = 0;
				}
				break;
			}
		}
	}
}


void 
HaMolView::ColourRibbonAttrib( int flag, int r, int g, int b )
{
    HaResidue  *group;
    HaAtom  *aptr;
	
	ColourRibbonNone( flag );

	MoleculesType::iterator mol_itr;
	MolSet* pmset = this->GetMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);

	ForEachMol_VIEW
	{
		if( !(*mol_itr)->IsSecStructFound()  )
		{
			p_mol_editor->DetermineSecStructure((*mol_itr),False);
		}
	}
	
	HaColor rcolor(r,g,b);
	
	AtomIteratorMolSet aitr(pmset);
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if( aptr->Selected() && 
			(aptr->IsAlphaCarbon()||
			aptr->IsSugarPhosphate()) )
		{  
			group=aptr->GetHostRes();
			if( flag & RibColInside )
			{   
				group->col1 = rcolor.cidx;
			}
			if( flag & RibColOutside )
			{   
				group->col2 = rcolor.cidx;
			}
			break;
		}
	}
}


void 
HaMolView::ColourMonitNone()
{
    int flag;

	list<Monitor>::iterator mitr;

	for( mitr=MonitList.begin(); mitr != MonitList.end(); mitr++ )
		if( (*mitr).col )
		{   
			flag = ZoneBoth ? (*mitr).src->Selected() && (*mitr).dst->Selected()
				            : (*mitr).src->Selected() || (*mitr).dst->Selected();
			if( flag )
			{   
				(*mitr).col = 0;
			}
		}
}

void 
HaMolView::ColourMonitAttrib( int r, int g, int b )
{
    int flag;
    ColourMonitNone();

	HaColor color(r,g,b);

	list<Monitor>::iterator mitr;

    for( mitr=MonitList.begin(); mitr != MonitList.end(); mitr++ )
    {   
		flag = ZoneBoth ? (*mitr).src->Selected() && (*mitr).dst->Selected()
			: (*mitr).src->Selected() || (*mitr).dst->Selected();
        if( flag )
        {   
            (*mitr).col = color.cidx;
        }
    }
}


void 
HaMolView::ColourDotsAttrib( int r, int g, int b )
{
    int i;
		
	HaColor color(r,g,b);
	MolSet* pmset = GetMolSet();

 	list<Object3D*>::iterator oitr;

	for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end();)
	{
		if( (*oitr)->GetObjType() != OBJ3D_DOT_SURFACE )
			continue;

		DotStruct* ptr = (DotStruct*) (*oitr);
		int np = ptr->GetCount();
		for( i=0; i< np; i++ )
		{   
			ptr->dots[i].col = color.cidx;
		}
	}
}



void HaMolView::ColourDotsPotential()
{
    DotStruct  *ptr;
    int i;	
	double result;

	class PotColorMap: public DRangeColorMap
	{
	public:
		PotColorMap()
		{
	       AddColor(255,   0,   0);  // 0  Red     25 < V       
	       AddColor(255, 165,   0);  /* 1  Orange  10 < V <  25 */
           AddColor(255, 255,   0);  /* 2  Yellow   3 < V <  10 */
           AddColor( 0,  255,   0);  /* 3  Green    0 < V <   3 */
	       AddColor( 0, 255, 255 );  /* 4  Cyan    -3 < V <   0 */
           AddColor( 0,   0, 255 );  /* 5  Blue   -10 < V <  -3 */
           AddColor( 160, 32, 240);  /* 6  Purple -25 < V < -10 */
           AddColor( 255, 255, 255); /* 7  White        V < -25 */
		}

		virtual ~PotColorMap()
		{

		}

		virtual HaColor* GetColorForVal(double val)
		{
			int idx = 0;
			if( val >= 0 )
			{   
				if( val > 10.0 )
				{      
					if( val > 24.0 )
					{      
						idx = 0;
					} 
					else 
						idx = 1;
				} 
				else if( val > 3.0 )
				{      
					idx = 2;
				} 
				else 
					idx = 3;
			} 
			else
			{
				if( val > -10.0 )
				{      
					if( val > -3.0 )
					{      
						idx = 4;
					} 
					else 
						idx = 5;
				} 
				else if( val > -24.0 )
				{      
					idx = 6;
				} 
				else 
					idx = 7;
			}
			return &colors[idx];
		}	
	};

	PotColorMap	pot_color_map;

	MolSet* pmset = GetMolSet();

 	list<Object3D*>::iterator oitr;

	for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end();)
	{
// found object of the type dot surface (DotStruct) in the list of 3D objects
		if( (*oitr)->GetObjType() != OBJ3D_DOT_SURFACE ) 
			continue; 

		ptr = (DotStruct*) *oitr;

		int np = ptr->GetCount(); // found number of points in the dot surface(DotStruct) 

		for( i=0; i< np; i++ ) // Cycle over all dots in the DotStruct
		{   
			result = pmset->CalculatePotential( ptr->dots[i].GetX(),
			ptr->dots[i].GetY(),
			ptr->dots[i].GetZ() );
			
			HaColor* pcol = pot_color_map.GetColorForVal(result);	
			if(pcol) ptr->dots[i].col = pcol->cidx;
		}
	}
}

void 
HaMolView::MonoColourAttrib( int r, int g, int b )
{
    HaAtom  *aptr;

	HaColor color(r,g,b);

	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if( aptr->Selected() )
		{   
			aptr->col = color.cidx;
		}
	}
}

void HaMolView::ScaleColourAttrib( int attr )
//! Colour Atoms according to attributes(properties)
{
    int count;
	int num_cols = 0;
	double range, valmin;
	double min,max;

//  range - range of values (difference between minimum and maximum)
//  valmin - minimal value of the property

    HaChain*    chain;
    HaResidue*  group;
	         
    HaAtom*  aptr;
	
    AtomIteratorMolSet aitr(this->GetMolSet());

	if(GetMolSet()->HostMolecules.empty()) return;

	MoleculesType::iterator mol_itr;

    switch( attr )
    {   
	    case(ChainAttr):   
			range= 0.0;
			ForEachMol_VIEW
			{
				range += (*mol_itr)->GetNChains();
				num_cols += (*mol_itr)->GetNChains();
			}
			range -= 1.0;
			if(range < 0.5) range = 1.0;
			valmin = 1.0;
			break;

        case(ResidueAttr): 
			min= 100000.0; max= -100000.0;
			{
				ResidueIteratorMolSet ritr( this->GetMolSet() );
				for( group = ritr.GetFirstRes(); group; group = ritr.GetNextRes() )
				{
					if(group->GetSerNo() < min) min= group->GetSerNo();
					if(group->GetSerNo() > max) max= group->GetSerNo();
				}
			}
			valmin= min;
			range= max - min;
			if(range < 0.5) range = 1.0;
			num_cols = 10;
			break;
		

        case(ChargeAttr):
        case(TempAttr):    
			min= 100000.0; max= -100000.0;
		    for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
			{
				if(aptr->tempf < min) min= aptr->tempf;
				if(aptr->tempf > max) max= aptr->tempf;
			}
			valmin= min;
			range= max - min;
			num_cols = 10;
			if(range < 10E-16) range = 1.0;
			break;

        default:           
			return;
    }

//    if( range < 2 )
//    {   
//		MonoColourAttrib(255,255,255);
//        return;
//    }

//		for( i=0; i<ScaleCount; i++ )
//		{   
//			fract = (int)((1023*i)/(ScaleCount-1));
//			if( fract < 256 )
//			{   
//				r = 0;  g = fract;  b = 255;
//			} 
//			else if( fract < 512 )
//			{   
//				r = 0;  g = 255;  b = 511-fract;
//			} 
//			else if( fract < 768 )
//			{   
//				r = fract-512;  g = 255;  b = 0;
//			} 
//			else /* fract < 1024 */                             
//			{   
//				r = 255;  g = 1023-fract;  b = 0;
//			}
//			ScaleRef[i].r = r;
//			ScaleRef[i].g = g;
//			ScaleRef[i].b = b;
//			ScaleRef[i].shade = 0;
//			ScaleRef[i].col = 0;
//		}
//

	DRangeColorMap attr_map;
	
	HaColor blue_c(0,0,255);
	HaColor red_c(255,0,0);

	attr_map.AddUniformRange(num_cols,blue_c,red_c);
	attr_map.min_val = 	valmin;
	attr_map.max_val =  valmin + range;

    switch( attr )
    {    
	case(ChainAttr):
		count = 0;
		{
			ChainIteratorMolSet chitr( this->GetMolSet() );
			for( chain = chitr.GetFirstChain(); chain; chain = chitr.GetNextChain() )
			{   
				HaColor* pcol = attr_map.GetColorForVal(count);
				if(pcol == NULL) break;

				ResidueIteratorChain ritr_ch(chain);
				for( group = ritr_ch.GetFirstRes(); group; group = ritr_ch.GetNextRes() )
				{
					AtomIteratorAtomGroup aitr(group);
					for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
					{
						if( aptr->Selected() )
						{   
							aptr->col = pcol->cidx;
						}
					}
				}
				count++;
			}
		}
		break;
				
	case(ResidueAttr):
		{
			ResidueIteratorMolSet ritr(this->GetMolSet() );
			for( group = ritr.GetFirstRes(); group; group = ritr.GetNextRes() )
			{   
				HaColor* pcol = attr_map.GetColorForVal(group->serno);
				if(pcol == NULL) break;

				AtomIteratorAtomGroup aitr(group);
				for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
				{
					if( aptr->Selected() )
					{   
						aptr->col = pcol->cidx;
					}
				}
			}
		}
		break;
	
	case(TempAttr):
        for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			if( aptr->Selected() )
			{   
     			HaColor* pcol = attr_map.GetColorForVal(aptr->tempf);
				if(pcol == NULL) break;
				aptr->col = pcol->cidx;
			}
		}
		break;
		
	case(ChargeAttr):
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			if( aptr->Selected() )
			{   
     			HaColor* pcol = attr_map.GetColorForVal(aptr->tempf);
				if(pcol == NULL) break;
				aptr->col = pcol->cidx;
			}
		}
		break;
    }
}


void HaMolView::CPKColourAttrib()
{
    HaAtom  *aptr;

	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
        if( aptr->Selected() ) ColorAtomCPK(aptr);
	}
}

int HaMolView::ColorAtomCPK( HaAtom* aptr)
{
	if( aptr == NULL) return FALSE;
	int idx = ElementArr[aptr->GetElemNo()].cpkcol;
	HaColor& color = cpk_col_map.GetColorByIdx(idx);
    aptr->col = color.cidx;
	return TRUE;
}

void HaMolView::GroupsColourAttrib()
{
    HaAtom*    aptr;
	ChemGroup* gptr;  

	HaColor grey_c(200, 200, 200);

	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
        if( aptr->Selected() )
        {   
            aptr->col = grey_c.cidx;
        }
	}
	
	HaColor blue_c   ( 0,   0, 255);
	HaColor red_c    (255,   0,    0);
    HaColor yellow_c (255, 200,  50);
    HaColor green_c  (0,  255,   0);
	
	HaColor* pcol;

// Colour groups with alternating colours
	int ig=0;
	ChemGroupIterator gitr(GetMolSet());
	for(gptr = gitr.GetFirst(); gptr; gptr = gitr.GetNext())
	{   
		ig++;
		switch(ig%4)
		{
		case 0:
			pcol = &blue_c; // Blue
			break;
		case 1:
			pcol = &red_c; // Red
			break;
		case 2:
			pcol = &yellow_c; // Yellow
			break;
		case 3:
			pcol = &green_c; // Green
			break;
		default:
			pcol = &grey_c;
		}
		
		AtomIteratorAtomGroup aitr(gptr);
		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if(aptr->Selected()) 
			{
				aptr->col = pcol->cidx; 
			}
		}
	}	
}

void HaMolView::RigidClusterColourAttrib()
{
    HaAtom*    aptr;
	AtomGroup* gptr;  

	HaColor grey_c(200, 200, 200);

	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
        if( aptr->Selected() )
        {   
            aptr->col = grey_c.cidx;
        }
	}
	
	HaColor blue_c   ( 0,   0, 255);
	HaColor red_c    (255,   0,    0);
    HaColor yellow_c (255, 200,  50);
    HaColor green_c  (0,  255,   0);
	
	HaColor* pcol;

// Colour groups with alternating colours

	size_t found;
	int ig=0;
	AtomGroupIteratorMolSet gitr(GetMolSet());
	for(gptr = gitr.GetFirst(); gptr; gptr = gitr.GetNext())
	{   
		std::string grp_name = gptr->GetID();
		found = grp_name.find("RIGIDCLST");
		if (found == string::npos) continue;

		ig++;
		switch(ig%4)
		{
		case 0:
			pcol = &blue_c; // Blue
			break;
		case 1:
			pcol = &red_c; // Red
			break;
		case 2:
			pcol = &yellow_c; // Yellow
			break;
		case 3:
			pcol = &green_c; // Green
			break;
		default:
			pcol = &grey_c;
		}
		
		if(gptr->GetNAtoms() < 10) pcol = &grey_c;

		AtomIteratorAtomGroup aitr(gptr);
		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			if(aptr->Selected()) 
			{
				aptr->col = pcol->cidx; 
			}
		}
	}	
}

void  HaMolView::AminoColourAttrib()
{
    HaResidue  *group;
    HaAtom  *aptr;

	StrColorMap AminoColors;

	AminoColors.AddStrColorPair("ASP", 230,  10,  10);
	AminoColors.AddStrColorPair("GLU", 230,  10,  10);
	AminoColors.AddStrColorPair("LYS", 20,   90, 255);
	AminoColors.AddStrColorPair("ARG", 20,   90, 255);
    AminoColors.AddStrColorPair("HIS", 130, 130, 210);
    AminoColors.AddStrColorPair("SER", 250, 150,   0);
    AminoColors.AddStrColorPair("THR", 250, 150,   0);
	AminoColors.AddStrColorPair("ASN", 0,   220, 220);
    AminoColors.AddStrColorPair("GLN", 0,   220, 220);
    AminoColors.AddStrColorPair("CYS", 230, 230,   0);
    AminoColors.AddStrColorPair("MET", 230, 230,   0);
    AminoColors.AddStrColorPair("ALA", 200, 200, 200);
    AminoColors.AddStrColorPair("GLY", 235, 235, 235);
    AminoColors.AddStrColorPair("LEU", 15, 130,   15);
    AminoColors.AddStrColorPair("VAL", 15, 130,   15);
    AminoColors.AddStrColorPair("ILE", 15, 130,   15);
	AminoColors.AddStrColorPair("PHE", 50,  50,  170);
	AminoColors.AddStrColorPair("TYR", 50,  50,  170);
    AminoColors.AddStrColorPair("TRP", 180,  90, 180);
	AminoColors.AddStrColorPair("PRO", 220, 150, 130);
	AminoColors.AddStrColorPair("PCA", 220, 150, 130);
	AminoColors.AddStrColorPair("HYP", 220, 150, 130);
	AminoColors.AddStrColorPair("XXX", 190, 160, 110);

	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
        if( aptr->Selected() )
        {   
			group=aptr->GetHostRes();
			HaColor* pcol = AminoColors.GetColorForStr(group->GetName());
			if(pcol == NULL)
			{
				pcol = AminoColors.GetColorForStr("XXX");
			}

            aptr->col = pcol->cidx;
        }
	}
}


void HaMolView::ShapelyColourAttrib()
{
    HaResidue  *group;
    HaAtom  *aptr;

    HaColorMap shape_col_map;

    shape_col_map.AddColor( 140, 255, 140  );    /* ALA */
    shape_col_map.AddColor(  255, 255, 255 );    /* GLY */
    shape_col_map.AddColor(  69,  94,  69  );    /* LEU */
    shape_col_map.AddColor(  255, 112,  66 );    /* SER */
    shape_col_map.AddColor(  255, 140, 255 );    /* VAL */
    shape_col_map.AddColor(  184,  76,   0 );    /* THR */
    shape_col_map.AddColor(   71,  71, 184 );    /* LYS */
    shape_col_map.AddColor(  160,   0,  66 );    /* ASP */
    shape_col_map.AddColor(    0,  76,   0 );    /* ILE */
    shape_col_map.AddColor(  255, 124, 112 );    /* ASN */
    shape_col_map.AddColor(  102,   0,   0 );    /* GLU */
    shape_col_map.AddColor(   82,  82,  82 );    /* PRO */
    shape_col_map.AddColor(    0,   0, 124 );    /* ARG */
    shape_col_map.AddColor(   83,  76,  66 );    /* PHE */
    shape_col_map.AddColor(  255,  76,  76 );    /* GLN */
    shape_col_map.AddColor(  140, 112,  76 );    /* TYR */
    shape_col_map.AddColor(  112, 112, 255 );    /* HIS */
    shape_col_map.AddColor(  255, 255, 112 );    /* CYS */
    shape_col_map.AddColor(  184, 160,  66 );    /* MET */
    shape_col_map.AddColor(   79,  70,   0 );    /* TRP */

    shape_col_map.AddColor(  255,   0, 255 );    /* ASX */
    shape_col_map.AddColor(  255,   0, 255 );    /* GLX */
    shape_col_map.AddColor(  255,   0, 255 );    /* PCA */
    shape_col_map.AddColor(  255,   0, 255 );    /* HYP */

    shape_col_map.AddColor(  160, 160, 255 );    /*   A */
    shape_col_map.AddColor(  255, 140,  75 );    /*   C */
    shape_col_map.AddColor(  255, 112, 112 );    /*   G */
    shape_col_map.AddColor(  160, 255, 160 );    /*   T */

    shape_col_map.AddColor(  184, 184, 184 );    /* 28 -> BackBone */
    shape_col_map.AddColor(   94,   0,  94 );    /* 29 -> Special  */
    shape_col_map.AddColor(  255,   0, 255 );  /* 30 -> Default  */

	HaColor* pcol;

	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
        if( aptr->Selected() )
        {   
			group=aptr->GetHostRes();
			if( group->IsAminoNucleo() && group->refno < 30)
            {   
				pcol = &shape_col_map.GetColorByIdx(group->refno);
            } 
			else 
				pcol = &shape_col_map.GetColorByIdx(30);

/*  Original Colour Scheme
 *
 *  ref = &(Shapely[26]);
 *  if( IsNucleo(group->refno) )
 *  {   ref = Shapely + group->refno;
 *  } else if( IsShapelyBackbone(aptr->refno) )
 *  {   ref = &(Shapely[24]);
 *  } else if( IsShapelySpecial(aptr->refno) )
 *  {   ref = &(Shapely[25]);
 *  } else if( IsAmino(group->refno) )
 *      ref = Shapely + group->refno;
 */
            aptr->col = pcol->cidx;
        }
	}
}

void HaMolView::StructColourAttrib()
{
    HaResidue  *group;
    HaAtom  *aptr;

	MoleculesType::iterator mol_itr;
	MolSet* pmset = this->GetMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor();

	ForEachMol_VIEW
	{
		if( !(*mol_itr)->IsSecStructFound() )
			p_mol_editor->DetermineSecStructure((*mol_itr),False);
	}
	
	HaColorMap struct_color_map;

    struct_color_map.AddColor(255, 255, 255 );    /* 0  Default     */
    struct_color_map.AddColor(255,   0, 128 );    /* 1  Alpha Helix */
    struct_color_map.AddColor(255, 200,   0 );    /* 2  Beta Sheet  */
    struct_color_map.AddColor(96, 128, 255  );    /* 3  Turn        */
	
    AtomIteratorMolSet aitr(pmset);
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		HaColor* pcol;
        if( aptr->Selected() )
        {   
			group=aptr->GetHostRes();
			if( group->struc & HelixFlag )
            {   
				pcol = &struct_color_map.GetColorByIdx(1);
            } 
			else if( group->struc & SheetFlag )
            {   
				pcol = &struct_color_map.GetColorByIdx(2);
            } 
			else if( group->struc & TurnFlag )
            {   
				pcol = &struct_color_map.GetColorByIdx(3);
            } 
			else 
			{
				pcol = &struct_color_map.GetColorByIdx(0);
			}
			
            aptr->col = pcol->cidx;
        }
	}
}

int HaMolView::IsCPKColour( HaAtom* aptr )
{
	int idx = ElementArr[aptr->GetElemNo()].cpkcol;
	HaColor& color = cpk_col_map.GetColorByIdx(idx);
    return(  aptr->col == color.cidx );
}


void HaMolView::DefaultRepresentation()
{
	ReDrawFlag |= RFRefresh | RFColour;

	MolSet* pmset= this->GetMolSet();

	int nb = pmset->GetNBonds();
	
	if( nb == 0 )
	{   
		EnableBackbone(CylinderFlag, 0.32 );
	} 
	else 
	{
		EnableWireframe(WireFlag,0);
	}
	
	CPKColourAttrib();
}

void
HaMolView::CenterSelected()
{
	MolSet* pmset = GetMolSet();
	if(pmset == NULL) 
		return;

	if(pmset->GetNMol() == 0) return;
	
	AtomIteratorMolSet aitr(pmset);
	HaAtom* aptr;

	double MinX_v, MinY_v, MinZ_v;
	double MaxX_v, MaxY_v, MaxZ_v;

	MinX_v = 10.0E6;
	MinY_v = 10.0E6;
	MinZ_v = 10.0E6;
	MaxX_v = -10.0E6;
	MaxY_v = -10.0E6;
	MaxZ_v = -10.0E6;

	int n_sel_at = 0;

	double x,y,z; 

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if(aptr->Selected())
		{
            x = aptr->GetX();
			y = aptr->GetY();
			z = aptr->GetZ();

//			GetTransfCoord(aptr->GetX(), aptr->GetY(), aptr->GetZ(), x, y, z);
		
			if( x < MinX_v) MinX_v = x;
			if( y < MinY_v) MinY_v = y;
			if( z < MinZ_v) MinZ_v = z;
			if( x > MaxX_v) MaxX_v = x;
			if( y > MaxY_v) MaxY_v = y;
			if( z > MaxZ_v) MaxZ_v = z;

			n_sel_at++;
		}
	}

	if(n_sel_at == 0)
	{
		ErrorInMod(" HaMolView::CenterSelected() ", 
			       "No atoms selected");
		return;
	}

	double dx = MaxX_v - MinX_v; 
	double dy = MaxY_v - MinY_v; 
	double dz = MaxZ_v - MinZ_v; 

// Fix dimensions for the linear and one-atom molecules
	if(dx < 1.0)
	{
		dx = 3.0;
		MinX_v -= dx/2.0;
		MaxX_v += dx/2.0;
	}
	if(dy < 1.0)
	{
		dy = 3.0;
		MinY_v -= dy/2.0;
		MaxY_v += dy/2.0;
	}
	if(dz < 1.0)
	{
		dz = 3.0;
		MinZ_v -= dz/2.0;
		MaxZ_v += dz/2.0;
	}

	Orig(1) = -(dx/2.0 + MinX_v);
	Orig(2) = -(dy/2.0 + MinY_v);
	Orig(3) = -(dz/2.0 + MinZ_v);

	double fmax = MaxFun(dx,dy);
	fmax = MaxFun(fmax,dz);

	if(fmax < 8.0) fmax = 8.0;
//	DScale = 1.0/(1.2*fmax);
    Zoom = 1.0/(DScale*1.2*fmax);
//	if(Zoom < 0.1) Zoom = 0.1;

	CalcRotCenter(TRUE);

//	Zoom = 1.0;

    x = - CenX;
	y = - CenY;
	z = - CenZ;

	Orig(1) += CenX + Rot(1,1)*x + Rot(1,2)*y + Rot(1,3)*z;
	Orig(2) += CenY + Rot(2,1)*x + Rot(2,2)*y + Rot(2,3)*z;
	Orig(3) += CenZ + Rot(3,1)*x + Rot(3,2)*y + Rot(3,3)*z;

	ReDrawFlag |= (RFRefresh | RFApply);
	
	pmset = GetMolSet();
	pmset->RefreshAllViews();
}

void 
HaMolView::InitialTransform()
{
    double fdist,fmax;
    HaAtom  *aptr;
	
	double MinX_v, MinY_v, MinZ_v;
	double MaxX_v, MaxY_v, MaxZ_v;

	MolSet* pmset = GetMolSet();
	if(pmset == NULL) 
		return;

	if(pmset->GetNMol() == 0 && pmset->ViewObjects.size() == 0) return;

	pmset->GetMinMaxCrd(MinX_v,MinY_v,MinZ_v,MaxX_v,MaxY_v,MaxZ_v); 


	list<Object3D*>::iterator oitr;

	for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end(); oitr++)
	{
		if( (*oitr)->GetObjType() ==  OBJ3D_MOLECULE ||  (*oitr)->GetObjType() ==  OBJ3D_SURFACE || (*oitr)->GetObjType() ==  OBJ3D_DOT_SURFACE)
			continue;
		if((*oitr)->IsDisplayed())
		{
			double MinX_v_tmp, MinY_v_tmp, MinZ_v_tmp;
			double MaxX_v_tmp, MaxY_v_tmp, MaxZ_v_tmp;
			if((*oitr)->GetObjectMinMaxCrd(MinX_v_tmp,MinY_v_tmp,MinZ_v_tmp,MaxX_v_tmp,MaxY_v_tmp,MaxZ_v_tmp))
			{
				if( MinX_v_tmp < MinX_v) MinX_v = MinX_v_tmp;
				if( MinY_v_tmp < MinY_v) MinY_v = MinY_v_tmp;
				if( MinZ_v_tmp < MinZ_v) MinZ_v = MinZ_v_tmp;
				if( MaxX_v_tmp > MaxX_v) MaxX_v = MaxX_v_tmp;
				if( MaxY_v_tmp > MaxY_v) MaxY_v = MaxY_v_tmp;
				if( MaxZ_v_tmp > MaxZ_v) MaxZ_v = MaxZ_v_tmp;
			}

		}
	}
	
	double dx = MaxX_v - MinX_v; 
	double dy = MaxY_v - MinY_v; 
	double dz = MaxZ_v - MinZ_v; 

// Fix dimensions for the linear and one-atom molecules
	if(dx < 1.0)
	{
		dx = 3.0;
		MinX_v -= dx/2.0;
		MaxX_v += dx/2.0;
	}
	if(dy < 1.0)
	{
		dy = 3.0;
		MinY_v -= dy/2.0;
		MaxY_v += dy/2.0;
	}
	if(dz < 1.0)
	{
		dz = 3.0;
		MinZ_v -= dz/2.0;
		MaxZ_v += dz/2.0;
	}

	Orig(1) = -(dx/2.0 + MinX_v);
	Orig(2) = -(dy/2.0 + MinY_v);
	Orig(3) = -(dz/2.0 + MinZ_v);
	
	Rot = 0.0;
	Rot(1,1) = Rot(2,2) = Rot(3,3) = 1.0;
	Zoom = 1.0;

    SetZOffset(10000);
		
	fmax = 0.0;

	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{   
		double ax, ay, az;

		GetTransfCoord(aptr->GetX(), aptr->GetY(), aptr->GetZ(), ax, ay, az);
	
		fdist = ax*ax + ay*ay + az*az;
		if( fdist > fmax )
			fmax = fdist;
	}
	for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end(); oitr++)
	{
		if( (*oitr)->GetObjType() ==  OBJ3D_MOLECULE ||  (*oitr)->GetObjType() ==  OBJ3D_SURFACE || (*oitr)->GetObjType() ==  OBJ3D_DOT_SURFACE)
			continue;
		if((*oitr)->IsDisplayed())
		{
			double MinX_v_tmp, MinY_v_tmp, MinZ_v_tmp;
			double MaxX_v_tmp, MaxY_v_tmp, MaxZ_v_tmp;
			if((*oitr)->GetObjectMinMaxCrd(MinX_v_tmp,MinY_v_tmp,MinZ_v_tmp,MaxX_v_tmp,MaxY_v_tmp,MaxZ_v_tmp))
			{
				double ax0, ay0, az0;
				double ax1, ay1, az1;

				GetTransfCoord(MinX_v_tmp, MinY_v_tmp, MinZ_v_tmp, ax0, ay0, az0);
				GetTransfCoord(MaxX_v_tmp, MaxY_v_tmp, MaxZ_v_tmp, ax1, ay1, az1);
	
				fdist = ax0*ax0 + ay0*ay0 + az0*az0;
				if( fdist > fmax )fmax = fdist;
				fdist = ax1*ax1 + ay0*ay0 + az0*az0;
				if( fdist > fmax )fmax = fdist;
				fdist = ax0*ax0 + ay1*ay1 + az0*az0;
				if( fdist > fmax )fmax = fdist;
				fdist = ax1*ax1 + ay1*ay1 + az0*az0;
				if( fdist > fmax )fmax = fdist;
				fdist = ax0*ax0 + ay0*ay0 + az1*az1;
				if( fdist > fmax )fmax = fdist;
				fdist = ax1*ax1 + ay0*ay0 + az1*az1;
				if( fdist > fmax )fmax = fdist;
				fdist = ax0*ax0 + ay1*ay1 + az1*az1;
				if( fdist > fmax )fmax = fdist;
				fdist = ax1*ax1 + ay1*ay1 + az1*az1;
				if( fdist > fmax )fmax = fdist;
			}

		}
	}
 	if(fmax < 6.0) fmax = 6.0;

	CalcRotCenter();
    DScale = 1.0/(2.0*sqrt(fmax));
}


void HaMolView::PrepareTransform()
{
    double theta;
	double cost;  
	double sint;
	double x,y,z;

	MolSet* pmset = GetMolSet();

	list<Object3D*>::iterator oitr;

	double CenX_tr;
	double CenY_tr;
	double CenZ_tr;
	HaMat_double Tr(3,3,0.0), Tr_m(3,3,0.0), scr(3,3);

	GetTransfCoord(CenX, CenY, CenZ, CenX_tr, CenY_tr, CenZ_tr);
    Vec3D cnt;
	cnt[0] = CenX; 
	cnt[1] = CenY; 
	cnt[2] = CenZ; 
	
    if( (ReDrawFlag&RFRotateX) && (CurRX != LastRX) )
    {   
		theta = PI*(CurRX - LastRX);
        LastRX = CurRX;

		cost = cos(theta);  
		sint = sin(theta);
 
		Tr(1,1) = 1.0; Tr(1,2) = 0.0;   Tr(1,3) = 0.0;
        Tr(2,1) = 0.0; Tr(2,2) = cost;  Tr(2,3) = sint;
		Tr(3,1) = 0.0; Tr(3,2) = -sint; Tr(3,3) = cost;

		if(m_screen_transform)
		{						
			matmult(scr, Tr, Rot);
			Rot = scr;
			
			y = Orig(2) - CenY_tr; 
			z = Orig(3) - CenZ_tr;
			Orig(2) = CenY_tr + cost*y + sint*z;
			Orig(3) = CenZ_tr + cost*z - sint*y;
		}
		else
		{
			matmult(scr,Tr,Rot);
			matmult_T1(Tr_m,Rot,scr);

			for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end(); oitr++)
			{
				if( (*oitr)->IsConnected())
				{
					(*oitr)->RotateObj(Tr_m,cnt);
				}
			}
		}
    }
	

    if( (ReDrawFlag&RFRotateY) && (CurRY != LastRY) )
    {   
		theta = PI*(CurRY - LastRY);
        LastRY = CurRY;
		
		cost = cos(theta);  
		sint = sin(theta);
		
		Tr(1,1) = cost;  Tr(1,2) = 0.0; Tr(1,3) = sint;
		Tr(2,1) = 0.0;   Tr(2,2) = 1.0; Tr(2,3) = 0.0;
		Tr(3,1) = -sint; Tr(3,2) = 0.0; Tr(3,3) = cost;
		
		if(m_screen_transform)
		{		
			matmult(scr, Tr, Rot);
			Rot = scr;
			
			x = Orig(1) - CenX_tr; 
			z = Orig(3) - CenZ_tr;
			Orig(1) = CenX_tr + cost*x + sint*z;
			Orig(3) = CenZ_tr + cost*z - sint*x;
		}
		else
		{
			matmult(scr,Tr,Rot);
			matmult_T1(Tr_m,Rot,scr);
			
			for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end(); oitr++)
			{
				if((*oitr)->IsConnected())
				{	
					(*oitr)->RotateObj(Tr_m,cnt);
				}
			}

		}
    }

 	
    if( (ReDrawFlag&RFRotateZ) && (CurRZ != LastRZ) )
    {   
		theta = PI*(CurRZ - LastRZ);
        LastRZ = CurRZ;
		
		cost = cos(theta);  
		sint = sin(theta);
		
		Tr(1,1) = cost; Tr(1,2) = -sint; Tr(1,3) = 0.0;
		Tr(2,1) = sint; Tr(2,2) = cost;  Tr(2,3) = 0.0;
		Tr(3,1) = 0.0;  Tr(3,2) = 0.0;   Tr(3,3) = 1.0;

		if(m_screen_transform)
		{			
			matmult(scr, Tr, Rot);
			Rot = scr;
			
			x = Orig(1) - CenX_tr; 
			y = Orig(2) - CenY_tr;
			Orig(1) = CenX_tr + cost*x - sint*y;
			Orig(2) = CenY_tr + cost*y + sint*x;
		}
		else
		{
			matmult(scr,Tr,Rot);
			matmult_T1(Tr_m,Rot,scr);

			for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end(); oitr++)
			{                    
				if( (*oitr)->IsConnected() )
				{    
					(*oitr)->RotateObj(Tr_m,cnt);
				}
			}
		}
    }
	
    if( (ReDrawFlag & RFTransX) && (CurTX != LastTX) )
    {   
		double deltx = (CurTX - LastTX)*pCanv->XRange()/Scale;
		LastTX= CurTX;
		
		if(m_screen_transform)
		{
			Orig(1) += deltx; 
		}
		else
		{
			Vec3D tr_vec;
			tr_vec[0] = deltx * Rot(1,1);
			tr_vec[1] = deltx * Rot(1,2);
			tr_vec[2] = deltx * Rot(1,3);

			CenX += tr_vec[0];
			CenY += tr_vec[1];
			CenZ += tr_vec[2];
			for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end(); oitr++)
			{
				if( (*oitr)->IsConnected())
				{
					(*oitr)->Translate(tr_vec);				
				}
			}
		}
	}

    if( (ReDrawFlag & RFTransY) && (CurTY != LastTY) )
    {   
		double delty= (CurTY - LastTY)*pCanv->YRange()/Scale;
		LastTY = CurTY;

		if(m_screen_transform)
		{ 
			Orig(2) += delty;
		}
		else
		{
			Vec3D tr_vec;

			tr_vec[0] = delty * Rot(2,1);
			tr_vec[1] = delty * Rot(2,2);
			tr_vec[2] = delty * Rot(2,3);

			CenX += tr_vec[0];
			CenY += tr_vec[1];
			CenZ += tr_vec[2];

			for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end(); oitr++)
			{
				if( (*oitr)->IsConnected() )
				{
					(*oitr)->Translate(tr_vec);				
				}
			}
		}
	}
}

static int apply_count=0;

void HaMolView::ApplyTransform()
{	
    HaChain  *chain;
    set<HaHBond, less<HaHBond> >::iterator  hptr;
    HaBond  *bptr;
    HaAtom  *aptr;
    HaHBond* phb;

    if( ReDrawFlag & RFMagnify )
    {   
		if( Zoom<0.1 ) Zoom=0.1;
		
        Scale = Zoom*DScale*pCanv->Range();
	    pCanv->m_ImageSize= (int)(Scale/DScale);
	    pCanv->m_ImageRadius= pCanv->m_ImageSize/2;
        if( GetImageSize() < 2 )
        {   
	       pCanv->m_ImageSize= 2;
	       pCanv->m_ImageRadius= 1;
        } 
    }

    if( (ReDrawFlag & RFRotate) | (ReDrawFlag & RFTransX) | (ReDrawFlag & RFTransY) )
    {   
		PrepareTransform();
    }

	MolSet* pmset = GetMolSet();

	MoleculesType::iterator mol_itr;

	ForEachMol_VIEW
	{
		(*mol_itr)->SetAtomScreenCoord(this);
	}
	
	
	if( ReDrawFlag & RFMagnify )
	{   
		if( DrawAtoms)
		{
		    AtomIteratorMolSet aitr(pmset);
		    for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
			{
				if( aptr->IsDrawSphere() )
				{   
					aptr->irad = (int)(Scale*aptr->image_radius);
				}					
			}
		}
	    
		BondIteratorMolSet bitr(pmset);
		if( DrawBonds )
		{
			for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond()) 
			{
				if( bptr->flag & CylinderFlag )
				{   
					bptr->irad = (int)(Scale*bptr->radius);
				}
			}
			
			for( hptr= pmset->HBonds.begin(); hptr != pmset->HBonds.end() ; hptr++ )
			{
				phb = (HaHBond*) &(*hptr);
				if( phb->IsToDrawCylinder() )
					phb->irad = (int)(Scale* phb->radius);
			}

			pmset->GetMolEditor()->UpdateBackBone(pmset);
			int nb = pmset->BackboneBonds.size();
			int ib;
			for( ib = 0; ib < nb; ib++ )
			{
				bptr = pmset->BackboneBonds[ib];
				if( bptr->flag&CylinderFlag )
					bptr->irad = (int)(Scale*bptr->radius);
			}
		}
	}
}


void HaMolView::ResetTransform()
{
    LastRX = LastRY = LastRZ = 0.0;
	LastTX=LastTY=0.0;
    CenX = CenY = CenZ = 0.0;
}

void HaMolView::CalcRotCenter(int sel_atoms)
{
	HaAtom* aptr;
	int nn=0;

	CenX=0.0;
	CenY=0.0;
	CenZ=0.0;

	MoleculesType::iterator mol_itr;

	ForEachMol_VIEW
	{
		if( (*mol_itr)->IsConnected() )
		{       
			AtomIteratorMolecule aitr(*mol_itr);
			for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
			{
				if( !sel_atoms || aptr->Selected())
				{
					CenX += aptr->GetX();
					CenY += aptr->GetY();
					CenZ += aptr->GetZ();
					nn++;
				}
			}
		}
	}
	if(nn != 0)
	{
		CenX = CenX/nn;
		CenY = CenY/nn;
		CenZ = CenZ/nn;
	}
}

void HaMolView::GetTransfCoord(double x_abs, double y_abs, double z_abs, double& x_tr, double& y_tr, double& z_tr) const
{						
	x_tr = Orig(1) + x_abs*Rot(1,1)  + y_abs*Rot(1,2)  + z_abs*Rot(1,3);
	y_tr = Orig(2) + x_abs*Rot(2,1)  + y_abs*Rot(2,2)  + z_abs*Rot(2,3);
	z_tr = Orig(3) + x_abs*Rot(3,1)  + y_abs*Rot(3,2)  + z_abs*Rot(3,3);
}
