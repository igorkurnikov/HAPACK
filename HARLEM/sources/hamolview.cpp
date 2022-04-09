/*! \file hamolview.cpp

    Classes to define 3D view of a molecule

    \author Igor Kurnikov 
    \date 1998-2002 

*/

#define HAMOLVIEW_CPP

#include "haconst.h"
#include "hamolview.h"

#include <stdlib.h>
#include <memory>

#include <string>
#include <stdexcept>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#include <chrono>
#include <thread>

#define RENDER

#include "abstree.h"
#include "command.h"

#include "haatom.h"
#include "habond.h"
#include "hamolecule.h"
#include "harlemapp.h"
#include "moleditor.h"
#include "tokens.h"
#include "haflexmod.h"

const double RootSix = 2.44948974278;

/* These define light source position */
const double LightLength = RootSix;
const double LightZComp  = 2.0;

static HaAtom  *Exclude;

/* Identified Atom Info */
static HaAtom* PickHist[4];
static int IdentDist;
static int IdentFound;
static int IdentDepth;
static int PickCount;
static int PickMode;

int HaMolView::UseHourGlass = True;

int HaMolView::MouseMode = MMRasMol;
int HaMolView::FakeSpecular = False;

double HaMolView::Ambient = DefaultAmbient;
int HaMolView::UseBackFade = False;
int HaMolView::SpecPower = 8;
int HaMolView::ZoneBoth  = True;

HaColorMap HaMolView::cpk_col_map;

HaMolView::HaMolView()
{
	host_mol_set = NULL;
	pCanv= new Canvas3D;

	InitializeRenderer();
 	ResetView();

	cpk_col_map.AddColor (200, 200, 200  );       /*  0 Light Grey   */
    cpk_col_map.AddColor ( 143, 143, 255 );       /*  1 Sky Blue     */
    cpk_col_map.AddColor ( 240,   0,   0 );       /*  2 Red          */
    cpk_col_map.AddColor ( 255, 200,  50 );       /*  3 Yellow       */
    cpk_col_map.AddColor ( 255, 255, 255 );       /*  4 White        */
    cpk_col_map.AddColor ( 255, 192, 203 );       /*  5 Pink         */
    cpk_col_map.AddColor ( 218, 165,  32 );       /*  6 Golden Rod   */
    cpk_col_map.AddColor (   0,   0, 255 );       /*  7 Blue         */
    cpk_col_map.AddColor ( 255, 165,   0 );       /*  8 Orange       */
    cpk_col_map.AddColor ( 128, 128, 144 );       /*  9 Dark Grey    */
    cpk_col_map.AddColor ( 165,  42,  42 );       /* 10 Brown        */
    cpk_col_map.AddColor ( 160,  32, 240 );       /* 11 Purple       */
    cpk_col_map.AddColor ( 255,  20, 147 );       /* 12 Deep Pink    */
    cpk_col_map.AddColor (   0, 255,   0 );       /* 13 Green        */
    cpk_col_map.AddColor ( 178,  34,  34 );       /* 14 Fire Brick   */
    cpk_col_map.AddColor ( 34, 139,  34  );        /* 15 Forest Green */

}

HaMolView::~HaMolView()
{
	if(CurMolView == this) CurMolView = NULL;
	ClearImage();
    ResetView();
	if(pCanv != NULL) delete pCanv;
}

void HaMolView::ResetView()
{
	SolventDots = False;
    ProbeRadius = 0.0;
    DrawLabels = False;
	LabelOptFlag = False;

	debug_level = 5;

	Zoom = 1.0;

	DrawMonitDistance = True;
    DrawBetaArrows = True;

	UseTransparent = False;
    UseOutLine = False;

    HetaGroups = True;
    Hydrogens = True;
	DrawBestPath = True;
	DrawContourSurf = True;
	DrawSolidSurfaces = TRUE;
    DrawObj3D = TRUE;
 
	Rot.newsize(3,3);

	Rot(1,1) = 1.0; Rot(2,1) = 0.0; Rot(3,1) = 0.0; 
	Rot(1,2) = 0.0; Rot(2,2) = 1.0; Rot(3,2) = 0.0; 
	Rot(1,3) = 0.0; Rot(2,3) = 0.0; Rot(3,3) = 1.0; 

	Orig.newsize(3);
	Orig(1) = 0.0; Orig(2) = 0.0; Orig(3) = 0.0;

	m_screen_transform = True;

	XOffset=YOffset=0;
	
 	CurRX = CurRY= CurRZ = 0.0;
	CurTX = CurTY = 0.0;
    CurSlabValue = 0.0;
	CurZoom = 0.0;

	UseScreenClip=True;

    ResetTransform();
    ResetRenderer();
   
    DeleteAllMonitors();
	
    DrawLabels = False;

    DrawMonitDistance = True;
    DrawBetaArrows = True;
    CartoonHeight = 0.4;
   
    ClearBuffers();
#ifndef _WIN32
    pCanv->FBClear = False;
#endif
	anim_thread_running = false;
	to_stop_animation = TRUE;
}

MolSet* HaMolView::GetMolSet()
{
	return host_mol_set;
}

unsigned int isqrt(unsigned int val )
{
#ifndef sun386
    register int i,result;
    register unsigned int temp;
    register unsigned int rem;

    i = 16;
    while( !(val&((unsigned int)3<<30)) && i )
    {   val <<= 2;
        i--;
    }

    if( i )
    {   rem = (val>>30)-1;
        val <<= 2;
        result = 1;
        i--;

        while( i )
        {   rem = (rem<<2) | (val>>30);
            result <<= 1;
            val <<= 2;

            temp = result<<1;
            if( rem > temp )
            {   rem -= temp|1;
                result |= 1;
            }
            i--;
        }
        return( result );
    } else return( 0 );
#else
    return( (int)sqrt((double)val) );
#endif
}


/*=============================*/
/*  ClearBuffers Subroutines!  */
/*=============================*/
 
#if defined(_WIN32) | defined(TWIN32)

void 
HaMolView::ClearBuffers()
{	
    if( !pCanv->FBClear )
    {   
		pCanv->FBClear = True;

		uint_4 fill = BackColor.cval;
#if defined(__DECCXX)
		PrintLog("Fill FBuffer with zeros \n");
		fill = 0;
#endif
		uint_4* ptr = (uint_4*)pCanv->FBuffer;
		uint_4* end = (uint_4*)(pCanv->FBuffer + (uint_4) pCanv->XRange() * pCanv->YRange());
		do { 
			*ptr++=fill; *ptr++=fill;
			*ptr++=fill; *ptr++=fill;
		} 
		while( ptr<end );
	}
	
    if( !pCanv->DBClear )
    {   
		pCanv->DBClear = True;
		memset(pCanv->DBuffer,0,(int)pCanv->XRange() * pCanv->YRange() *sizeof(short));
    }
}
#else
void 
HaMolView::ClearBuffers()
{
    uint_4 *ptr;
    uint_4 *end;
    uint_4 fill;
	
    if( !pCanv->FBClear )
    {   
		pCanv->FBClear = True;
		fill = BackColor.cval;

		ptr = (uint_4*)pCanv->FBuffer;
		end = (uint_4*)(pCanv->FBuffer + (uint_4) pCanv->XRange() * pCanv->YRange());
		do { 
			*ptr++=fill; *ptr++=fill;
			*ptr++=fill; *ptr++=fill;
		} 
		while( ptr<end );
    }
	
    if( !pCanv->DBClear )
    {   
		pCanv->DBClear = True;
		ptr = (uint_4*)pCanv->DBuffer;
		end = (uint_4*)(pCanv->DBuffer+(uint_4)pCanv->XRange() * pCanv->YRange());
		do { 
			*ptr++=0; *ptr++=0;
			*ptr++=0; *ptr++=0;
		} 
		while( ptr<end );
    }
}
#endif /* UNIX */



void HaMolView::ReAllocBuffers()
{	
	pCanv->AllocDBuffer();
    pCanv->DBClear=False;
}


void HaMolView::ReSizeScreen()
{	
	if( !pCanv->FBuffer || (FBufX!=pCanv->XRange()) || (FBufY!=pCanv->YRange()) )
	{   
		if( !CreateImage() )
		{ 
            ErrorInMod("HaMolView::ReSizeScreen()",
                       "Error create frame buffer");
            return;
	    }
		FBufX=pCanv->XRange();  FBufY=pCanv->YRange();  
		pCanv->FBClear = False;
		ReAllocBuffers();
		ClearBuffers();
	}
}

int HaMolView::GetImageSize()
{
	return(pCanv->m_ImageSize);
}

int HaMolView::GetImageRadius()
{
	return(pCanv->m_ImageRadius);
}

void HaMolView::BuildHashTable()
{
    int i;
	double xmin, xmax, ymin, ymax, zmin, zmax;
	
    HashTable.SetDimensions(21,21,21);

	MolSet* pmset = GetMolSet();
	
	pmset->GetMinMaxCrd( xmin, ymin, zmin, xmax, ymax, zmax);
	
	xmin -= 0.5;
	ymin -= 0.5;
	zmin -= 0.5;
	xmax += 0.5;
	ymax += 0.5;
	zmax += 0.5;
	
	HashTable.SetBoundaries(xmin, ymin, zmin, xmax, ymax, zmax);	

	HaAtom* aptr;
	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if( aptr->Selected() )
		{   
			i = HashTable.GetPointCellIdx(aptr->GetX(),aptr->GetY(),aptr->GetZ());
			HashTable[i].push_back(aptr);
		}
	}
}


void 
HaMolView::DisplaySpaceFill()
{
    HaAtom  *aptr;
   
   if( UseClipping )
   {   
		AtomIteratorMolSet aitr(this->GetMolSet());
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
		   if( aptr->IsDrawSphere())
			   pCanv->ClipSphere(aptr->x,aptr->y,aptr->z,aptr->irad,aptr->col);
		}
   } 
   else 
   {
		AtomIteratorMolSet aitr(this->GetMolSet());
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
	       if( aptr->IsDrawSphere() )
		      pCanv->DrawSphere(aptr->x,aptr->y,aptr->z,aptr->irad,aptr->col);
		}
   }
}


void HaMolView::DisplayWireframe()
{
    HaBond* bptr;
	HaAtom* aptr;
    HaAtom  *s;
    HaAtom  *d;
    int sc,dc;
	int size;

    if( UseClipping )
    {   
		BondIteratorMolSet bitr(GetMolSet());
		for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond()) // Plot Bonds
		{
			if( bptr->IsToDraw() )
			{   
				s = bptr->srcatom; d = bptr->dstatom;
				if( !bptr->col ) 
				{   
					sc = s->col;  dc = d->col;
				} 
				else
				{
					sc = dc = bptr->col;
				}

//				PrintLog("Draw bond with color idx %d-%d \n",sc, dc);
				
				if( bptr->flag&WireFlag )
				{   
					if (bptr->IsDouble())
					{
						pCanv->ClipTwinVector(s->x - 2, s->y - 2,  s->z - 2, d->x - 2, d->y - 2 , d->z - 2, sc, dc);
						pCanv->ClipTwinVector(s->x + 2, s->y + 2 , s->z + 2, d->x + 2, d->y + 2,  d->z + 2, sc, dc);
					}
					else if ( bptr->IsAromatic() )
					{
						pCanv->ClipTwinVector(s->x - 2, s->y - 2, s->z - 2, d->x - 2, d->y - 2, d->z - 2, sc, dc);
						pCanv->ClipDashVector(s->x + 2, s->y + 2, s->z + 2, d->x + 2, d->y + 2, d->z + 2, sc, dc);
					}
					else
					{
						pCanv->ClipTwinVector(s->x, s->y, s->z, d->x, d->y, d->z, sc, dc);
					}
				} 
				else if( bptr->flag&CylinderFlag )
				{   
					if( bptr->irad>0 )
					{  
						pCanv->ClipCylinder(s->x,s->y,s->z,d->x,d->y,d->z,
							sc,dc,bptr->irad);
					} 
					else 
						pCanv->ClipTwinVector(s->x,s->y,s->z,d->x,d->y,d->z,
						sc,dc);
				} 
				else /* bptr->flag & DashFlag */
				{
					pCanv->ClipDashVector(s->x, s->y, s->z, d->x, d->y, d->z, sc, dc);
				}
			}
		}

        // Plot non-bonded atoms as crosses
		AtomIteratorMolSet aitr(this->GetMolSet());     // Plot non-bonded atoms as crosses
        for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			if( aptr->Selected() && aptr->GetNBonds() == 0 )
			{
				size= (int)(1.0*Scale); 
				pCanv->ClipTwinVector(aptr->x - size, aptr->y, aptr->z,
					aptr->x + size, aptr->y, aptr->z, aptr->col, aptr->col);
				pCanv->ClipTwinVector(aptr->x, aptr->y - size, aptr->z,
					aptr->x, aptr->y + size , aptr->z, aptr->col, aptr->col);
				
			}
		}
		
    } 
	else
	{
		BondIteratorMolSet bitr(GetMolSet());
		for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond()) 
		{
			if( bptr->IsToDraw() )
			{   
				s = bptr->srcatom; d = bptr->dstatom;
				if( !bptr->col )
				{   
					sc = s->col;  dc = d->col;
				} 
				else 
					sc = dc = bptr->col;
				
				if( bptr->flag&WireFlag )
				{      
					pCanv->DrawTwinVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
				} 
				else if( bptr->flag&CylinderFlag )
				{   
					if( bptr->irad>0 )
					{  
						pCanv->DrawCylinder(s->x,s->y,s->z,d->x,d->y,d->z,
							sc,dc,bptr->irad);
					} 
					else 
						pCanv->DrawTwinVector(s->x,s->y,s->z,d->x,d->y,d->z,
						sc,dc);
				} 
				else 
					pCanv->ClipDashVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
			}
		}

		// Plot non-bonded atoms as crosses
		AtomIteratorMolSet aitr(this->GetMolSet());    // Plot non-bonded atoms as crosses
        for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			if( aptr->Selected() && aptr->GetNBonds() )
			{
				size= 10; 
				pCanv->DrawTwinVector(aptr->x - size, aptr->y , aptr->z -size,
					aptr->x + size, aptr->y , aptr->z +size, aptr->col, aptr->col);
				pCanv->DrawTwinVector(aptr->x , aptr->y - size, aptr->z -size,
					aptr->x , aptr->y + size, aptr->z +size, aptr->col, aptr->col);
				
			}
		}
	}

}


/* Used by DisplayDoubleBonds! */
void 
HaMolView::DisplayCylinder(int x1,int y1,int z1,
						   int x2,int y2,int z2,
						   int c1,int c2,int rad)
{
    if( UseClipping )
    {   
		if( rad == 0 )
        {   
			pCanv->ClipTwinVector(x1,y1,z1,x2,y2,z2,c1,c2);
        } 
		else 
			pCanv->ClipCylinder(x1,y1,z1,x2,y2,z2,c1,c2,rad);
    } 
	else
    {   
		if( rad == 0 )
        {   
			pCanv->DrawTwinVector(x1,y1,z1,x2,y2,z2,c1,c2);
        } 
		else 
			pCanv->DrawCylinder(x1,y1,z1,x2,y2,z2,c1,c2,rad);
    }
}


void 
HaMolView::DisplayDoubleBonds()
{
    register HaAtom  *s;
    register HaAtom  *d;
    register HaBond  *bptr;
    register int dx,dy,ix,iy;
    register int ax,ay,sc,dc;
    register int k,flag;
	
	BondIteratorMolSet bitr(GetMolSet());
	for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond()) 
        if( bptr->IsToDraw() )
        {   
			s = bptr->srcatom; d = bptr->dstatom;
            if( !bptr->col ) 
            {   
				sc = s->col;  dc = d->col;
            } 
			else sc = dc = bptr->col;
			
            flag = (bptr->flag&CylinderFlag) && (bptr->irad>4);
            if( !(bptr->flag & DoubBondFlag) || flag )
            {   
				if( bptr->flag&WireFlag )
                {   
					pCanv->ClipTwinVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
                } 
				else if( bptr->flag&CylinderFlag )
                {   
					DisplayCylinder(s->x,s->y,s->z,d->x,d->y,d->z,
						sc,dc,bptr->irad);
                } 
				else /* bptr->flag & DashFlag */
                    pCanv->ClipDashVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
            } 
			
            if( (bptr->flag & (DoubBondFlag|TripBondFlag)) && !flag )
            {   
				if( s->x > d->x )
                {   
					ax = s->x - d->x;  dx = 1;
                } 
				else
                {   
					ax = d->x - s->x;  dx = 0;
                }
				
                if( s->y > d->y )
                {   
					ay = s->y - d->y;  dy = 1;
                } 
				else 
                {  
					ay = d->y - s->y;  dy = 0;
                }
				
                /* Determine Bond Separation */
                if( (bptr->flag&CylinderFlag) && bptr->irad )
                {   
					if( bptr->flag & DoubBondFlag )
                    {   
						k = (3*bptr->irad+1)>>1;
                    } 
					else 
						k = 3*bptr->irad;
                } 
				else 
					k = (bptr->flag&TripBondFlag)?3:2;
				
                if( ax > (ay<<1) )
                {   
					ix = 0;  iy = k;
                } 
				else if( ay > (ax<<1) )
                {   
					iy = 0;  ix = k;
                } 
				else /* diagonal */
                {   
					k = (3*k)>>2;
                    if( dx == dy )
                    {   
						ix = k;  iy = -k;
                    } 
					else
                    {   
						ix = k;  iy = k;
                    }
                }
				
                if( bptr->flag&WireFlag )
                {   
					pCanv->ClipTwinVector(s->x+ix,s->y+iy,s->z,
						d->x+ix,d->y+iy,d->z,sc,dc);
                    pCanv->ClipTwinVector(s->x-ix,s->y-iy,s->z,
						d->x-ix,d->y-iy,d->z,sc,dc);
                } 
				else if( bptr->flag&CylinderFlag )
                {   
					DisplayCylinder(s->x+ix,s->y+iy,s->z,
						d->x+ix,d->y+iy,d->z,
						sc,dc,bptr->irad);
                    DisplayCylinder(s->x-ix,s->y-iy,s->z,
						d->x-ix,d->y-iy,d->z,
						sc,dc,bptr->irad);
                } 
				else /* bptr->flag & DashFlag */
                {  
					pCanv->ClipDashVector(s->x+ix,s->y+iy,s->z,
						d->x+ix,d->y+iy,d->z,sc,dc);
					pCanv->ClipDashVector(s->x+ix,s->y+iy,s->z,
						d->x+ix,d->y+iy,d->z,sc,dc);
                }
            }
        }
}


void HaMolView::DisplayBackbone()
{
    HaBond   *bptr;
    HaAtom   *s;
    HaAtom   *d;
    int sc,dc;

	MolSet* pmset = GetMolSet();

	pmset->GetMolEditor()->UpdateBackBone(pmset);
	int nb = pmset->BackboneBonds.size();
	int ib;
	for( ib = 0; ib < nb; ib++ )
	{
		bptr = pmset->BackboneBonds[ib];
		if( bptr->IsToDraw() )
		{   
			s = bptr->srcatom; d = bptr->dstatom;
			if( !bptr->col ) 
			{   
				sc = s->col;  dc = d->col;
			} 
			else 
				sc = dc = bptr->col;
			
			if( bptr->flag&CylinderFlag )
			{   
				if( bptr->irad>0 )
				{ 
					pCanv->ClipCylinder(s->x,s->y,s->z,d->x,d->y,d->z,
						sc,dc,bptr->irad);
				} 
				else 
					pCanv->ClipTwinVector(s->x,s->y,s->z,d->x,d->y,d->z,
					sc,dc);
			} 
			else if( bptr->flag & WireFlag )
			{      
				pCanv->ClipTwinVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
			} 
			else 
				pCanv->ClipDashVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
		}
	}
}


void HaMolView::DisplayHBonds()
{
	MolSet* pmset = GetMolSet();
	
	set<HaHBond, less<HaHBond> >::iterator  bitr;
	HaAtom  *s;
	HaAtom  *d;
	int sc,dc;
	
	for( bitr= pmset->HBonds.begin(); bitr != pmset->HBonds.end(); bitr++  )
	{
		HaHBond* phbond = (HaHBond*) (&(*bitr));
		if( (*bitr).IsToDraw() )
		{   
//				s = phbond->srcCA; d = phbond->dstCA;
//				if( !s || !d ) continue;
			d = (*bitr).src;
			s = (*bitr).dst;
			
			if( !(*bitr).col )
			{   
				sc = s->col;  dc = d->col;
			} 
			else 
			{
				sc = dc = (*bitr).col;
			}
			
			if( (*bitr).IsToDrawCylinder() )
			{   
				if( (*bitr).irad > 0 )
				{   
					pCanv->ClipCylinder(s->x,s->y,s->z,d->x,d->y,d->z, sc,dc,phbond->irad);
				} 
				else 
				{
					pCanv->ClipTwinVector(s->x,s->y,s->z,d->x,d->y,d->z, sc,dc);
				}
			} 
			else 
			{
				pCanv->ClipDashVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
			}
		}
	}
}

void HaMolView::DisplaySSBonds()
{
	MolSet* pmset = GetMolSet();
		
	HaAtom  *s;
	HaAtom  *d;
	int sc,dc;
	
	BondIteratorMolSet bitr(pmset);
	HaBond* pbond;
	for( pbond = bitr.GetFirstBond(); pbond; pbond = bitr.GetNextBond() )
	{
		if( pbond->srcatom->GetElemNo() != 16 || pbond->dstatom->GetElemNo() != 16 ) continue;

		if( pbond->IsToDraw() )
		{   
			d = pbond->srcatom;
			s = pbond->dstatom;
			
			if( !pbond->col )
			{   
				sc = s->col;  dc = d->col;
			} 
			else 
			{
				sc = dc = pbond->col;
			}
			
			if( pbond->flag & CylinderFlag )
			{   
				if( pbond->irad>0 )
				{   
					pCanv->ClipCylinder(s->x,s->y,s->z,d->x,d->y,d->z,
						sc,dc,pbond->irad);
				} 
				else 
					pCanv->ClipTwinVector(s->x,s->y,s->z,d->x,d->y,d->z,
					sc,dc);
			} 
			else 
				pCanv->ClipDashVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
		}
	}
}

void 
HaMolView::DisplayBoxes()
{
    double lena, lenb, lenc;
    double tmpx, tmpy, tmpz;
    double cosa, cosb, cosg;
    double temp, sing;

    int dxx,dxy,dxz;
    int dyx,dyy,dyz;
    int dzx,dzy,dzz;
    int x, y, z;
	int i;

	MolSet* pmset = GetMolSet();

	if(pmset->HostMolecules.empty()) return;

	int ixadd=  (int)(pCanv->XRange()/2.0);
	int iyadd=  (int)(pCanv->YRange()/2.0);
	int izadd=  (int)(ZOffset());

	MoleculesType::iterator mol_itr;

    if( DrawAxes  || DrawBoundBox )
    {  
		ForEachMol_VIEW
		{
			double minx, miny, minz;
			double maxx, maxy, maxz;
			(*mol_itr)->GetMinMaxCrd(minx, miny, minz, maxx, maxy, maxz);
			GetTransfCoord(maxx, maxy, maxz, maxx, maxy, maxz);
			
			dxx = (int)(maxx*Scale) + ixadd ;
			dxy = iyadd;
			dxz = izadd;
			
			dyx = ixadd;
			dyy = (int)(maxy*Scale) + iyadd;
			dyz = izadd;
			
			dzx = ixadd;
			dzy = iyadd;
			dzz = (int)(maxz*Scale) + izadd;
			
			if( DrawAxes )
			{   /* Line (MinX,0,0) to (MaxX,0,0) */
				x = dxx;  y = dxy;  z = dxz;
				if( ZValid_v(z) ) pCanv->DisplayTextString(x+2,y,z,"X",BoxColor.cidx);
				pCanv->ClipTwinVector(-dxx,-dxy,ZOffset()-dxz,
					x,y,z,BoxColor.cidx,BoxColor.cidx);
				
				/* Line (0,MinY,0) to (0,MaxY,0) */
				x = dyx;  y = dyy;  z = ZOffset()+dyz;
				if( ZValid_v(z) ) pCanv->DisplayTextString(x+2,y,z,"Y",BoxColor.cidx);
				pCanv->ClipTwinVector(-dyx,-dyy,ZOffset()-dyz, 
					x,y,z,BoxColor.cidx,BoxColor.cidx);
				
				
				/* Line (0,0,MinZ) to (0,0,MaxZ) */
				x = -dzx;  y = -dzy;  z = ZOffset()-dzz;
				if( ZValid_v(z) ) pCanv->DisplayTextString(x+2,y,z,"Z",BoxColor.cidx);
				pCanv->ClipTwinVector(+dzx,+dzy,ZOffset()+dzz, 
					x,y,z,BoxColor.cidx,BoxColor.cidx);
				
			}
			
			if( DrawBoundBox )
			{   /* Line (MinX,MinY,MinZ) to (MaxX,MinY,MinZ) */
				x=-dyx-dzx;  y=-dyy-dzy;  z=ZOffset()-dyz-dzz;
				pCanv->ClipTwinVector(x-dxx,y-dxy,z-dxz,x+dxx,y+dxy,z+dxz,BoxColor.cidx,BoxColor.cidx);
				
				/* Line (MaxX,MinY,MinZ) to (MaxX,MaxY,MinZ) */
				x=+dxx-dzx;  y=+dxy-dzy;  z=ZOffset()+dxz-dzz;
				pCanv->ClipTwinVector(x-dyx,y-dyy,z-dyz,x+dyx,y+dyy,z+dyz,BoxColor.cidx,BoxColor.cidx);
				
				/* Line (MaxX,MaxY,MinZ) to (MinX,MaxY,MinZ) */
				x=+dyx-dzx;  y=+dyy-dzy;  z=ZOffset()+dyz-dzz;
				pCanv->ClipTwinVector(x+dxx,y+dxy,z+dxz,x-dxx,y-dxy,z-dxz,BoxColor.cidx,BoxColor.cidx);
				
				/* Line (MinX,MaxY,MinZ) to (MinX,MinY,MinZ) */
				x=-dxx-dzx;  y=-dxy-dzy;  z=ZOffset()-dxz-dzz;
				pCanv->ClipTwinVector(x+dyx,y+dyy,z+dyz,x-dyx,y-dyy,z-dyz,BoxColor.cidx,BoxColor.cidx);
				
				
				/* Line (MinX,MinY,MinZ) to (MinX,MinY,MaxZ) */
				x=-dxx-dyx;  y=-dxy-dyy;  z=ZOffset()-dxz-dyz;
				pCanv->ClipTwinVector(x-dzx,y-dzy,z-dzz,x+dzx,y+dzy,z+dzz,BoxColor.cidx,BoxColor.cidx);
				
				/* Line (MaxX,MinY,MinZ) to (MaxX,MinY,MaxZ) */
				x=+dxx-dyx;  y=+dxy-dyy;  z=ZOffset()+dxz-dyz;
				pCanv->ClipTwinVector(x-dzx,y-dzy,z-dzz,x+dzx,y+dzy,z+dzz,BoxColor.cidx,BoxColor.cidx);
				
				/* Line (MaxX,MaxY,MinZ) to (MaxX,MaxY,MaxZ) */
				x=+dxx+dyx;  y=+dxy+dyy;  z=ZOffset()+dxz+dyz;
				pCanv->ClipTwinVector(x-dzx,y-dzy,z-dzz,x+dzx,y+dzy,z+dzz,BoxColor.cidx,BoxColor.cidx);
				
				/* Line (MinX,MaxY,MinZ) to (MinX,MaxY,MaxZ) */
				x=-dxx+dyx;  y=-dxy+dyy;  z=ZOffset()-dxz+dyz;
				pCanv->ClipTwinVector(x-dzx,y-dzy,z-dzz,x+dzx,y+dzy,z+dzz,BoxColor.cidx,BoxColor.cidx);
				
				
				/* Line (MinX,MinY,MaxZ) to (MaxX,MinY,MaxZ) */
				x=-dyx+dzx;  y=-dyy+dzy;  z=ZOffset()-dyz+dzz;
				pCanv->ClipTwinVector(x-dxx,y-dxy,z-dxz,x+dxx,y+dxy,z+dxz,BoxColor.cidx,BoxColor.cidx);
				
				/* Line (MaxX,MinY,MaxZ) to (MaxX,MaxY,MaxZ) */
				x=+dxx+dzx;  y=+dxy+dzy;  z=ZOffset()+dxz+dzz;
				pCanv->ClipTwinVector(x-dyx,y-dyy,z-dyz,x+dyx,y+dyy,z+dyz,BoxColor.cidx,BoxColor.cidx);
				
				/* Line (MaxX,MaxY,MaxZ) to (MinX,MaxY,MaxZ) */
				x=+dyx+dzx;  y=+dyy+dzy;  z=ZOffset()+dyz+dzz;
				pCanv->ClipTwinVector(x+dxx,y+dxy,z+dxz,x-dxx,y-dxy,z-dxz,BoxColor.cidx,BoxColor.cidx);
				
				/* Line (MinX,MaxY,MaxZ) to (MinX,MinY,MaxZ) */
				x=-dxx+dzx;  y=-dxy+dzy;  z=ZOffset()-dxz+dzz;
				pCanv->ClipTwinVector(x+dyx,y+dyy,z+dyz,x-dyx,y-dyy,z-dyz,BoxColor.cidx,BoxColor.cidx);
			}
		}
    }

    if( DrawUnitCell )
    {   /* Calculate Unit Cell! */
		if( pmset->per_bc->IsSet() )
		{
			lena = pmset->per_bc->GetA();
			lenb = pmset->per_bc->GetB();
			lenc = pmset->per_bc->GetC();
		
			cosa = cos(pmset->per_bc->GetAlpha());
			cosb = cos(pmset->per_bc->GetBeta());
			cosg = cos(pmset->per_bc->GetGamma());  
			sing = sin(pmset->per_bc->GetGamma());
		
			temp = cosa*cosa + cosb*cosb + cosg*cosg - 2.0*cosa*cosb*cosg;
			tmpx = cosb; 
			tmpy = (cosa - cosb*cosg)/sing;
			tmpz = -sqrt(1.0-temp)/sing;
			
			double xv[8],yv[8],zv[8];

			xv[0] = 0.0;    //			xv[0] =  - lena/2.0; 
			yv[0] = 0.0;    //			yv[0] =  - lenb/2.0;
			zv[0] = 0.0;    //			zv[0] =  - lenc/2.0;

			xv[1] = lena;   //			xv[1] =  + lena/2.0;
			yv[1] = 0.0;    //			yv[1] =  - lenb/2.0;
			zv[1] = 0.0;    //			zv[1] =  - lenc/2.0;

			xv[2] = 0.0;    //			xv[2] =  - lena/2.0; 
			yv[2] = lenb;   //			yv[2] =  + lenb/2.0;
			zv[2] = 0.0;    //			zv[2] =  - lenc/2.0;

			xv[3] = lena;  //			xv[3] =  + lena/2.0; 
			yv[3] = lenb;  //			yv[3] =  + lenb/2.0;
			zv[3] = 0.0;   //			zv[3] =  - lenc/2.0;

			xv[4] = 0.0;   //			xv[4] =  - lena/2.0; 
			yv[4] = 0.0;   //			yv[4] =  - lenb/2.0;
			zv[4] = lenc;  //			zv[4] =  + lenc/2.0;

			xv[5] = lena;  //			xv[5] =  + lena/2.0;
			yv[5] = 0.0;   //			yv[5] =  - lenb/2.0;
			zv[5] = lenc;  //			zv[5] =  + lenc/2.0;

			xv[6] = 0.0;   //			xv[6] =  - lena/2.0; 
			yv[6] = lenb;  //			yv[6] =  + lenb/2.0;
			zv[6] = lenc;  //			zv[6] =  + lenc/2.0;

			xv[7] = lena;  //			xv[7] =  + lena/2.0; 
			yv[7] = lenb;  //			yv[7] =  + lenb/2.0;
			zv[7] = lenc;  //			zv[7] =  + lenc/2.0;

			for( i = 0; i < 8; i++)
			{
				double x_tr, y_tr, z_tr;
				GetTransfCoord( xv[i], yv[i], zv[i], x_tr, y_tr, z_tr);
			
				xv[i] = x_tr * Scale + ixadd;
				yv[i] = y_tr * Scale + iyadd;
				zv[i] = z_tr * Scale + izadd;
			}

			/* Draw Unit Cell */
			pCanv->ClipTwinVector((int)xv[0],(int)yv[0],(int)zv[0],(int)xv[1],(int)yv[1],(int)zv[1],BoxColor.cidx,BoxColor.cidx);
			pCanv->ClipTwinVector((int)xv[0],(int)yv[0],(int)zv[0],(int)xv[2],(int)yv[2],(int)zv[2],BoxColor.cidx,BoxColor.cidx);
			pCanv->ClipTwinVector((int)xv[1],(int)yv[1],(int)zv[1],(int)xv[3],(int)yv[3],(int)zv[3],BoxColor.cidx,BoxColor.cidx);
			pCanv->ClipTwinVector((int)xv[2],(int)yv[2],(int)zv[2],(int)xv[3],(int)yv[3],(int)zv[3],BoxColor.cidx,BoxColor.cidx);

			pCanv->ClipTwinVector((int)xv[0],(int)yv[0],(int)zv[0],(int)xv[4],(int)yv[4],(int)zv[4],BoxColor.cidx,BoxColor.cidx);
			pCanv->ClipTwinVector((int)xv[1],(int)yv[1],(int)zv[1],(int)xv[5],(int)yv[5],(int)zv[5],BoxColor.cidx,BoxColor.cidx);
			pCanv->ClipTwinVector((int)xv[2],(int)yv[2],(int)zv[2],(int)xv[6],(int)yv[6],(int)zv[6],BoxColor.cidx,BoxColor.cidx);
			pCanv->ClipTwinVector((int)xv[3],(int)yv[3],(int)zv[3],(int)xv[7],(int)yv[7],(int)zv[7],BoxColor.cidx,BoxColor.cidx);

			pCanv->ClipTwinVector((int)xv[4],(int)yv[4],(int)zv[4],(int)xv[5],(int)yv[5],(int)zv[5],BoxColor.cidx,BoxColor.cidx);
			pCanv->ClipTwinVector((int)xv[4],(int)yv[4],(int)zv[4],(int)xv[6],(int)yv[6],(int)zv[6],BoxColor.cidx,BoxColor.cidx);
			pCanv->ClipTwinVector((int)xv[5],(int)yv[5],(int)zv[5],(int)xv[7],(int)yv[7],(int)zv[7],BoxColor.cidx,BoxColor.cidx);
			pCanv->ClipTwinVector((int)xv[6],(int)yv[6],(int)zv[6],(int)xv[7],(int)yv[7],(int)zv[7],BoxColor.cidx,BoxColor.cidx);
		}
	}
}

void
HaMolView::DisplayPickedAtoms()
{
	MolSet* pmset = GetMolSet();
	AtomIteratorAtomGroup aitr(&pmset->picked_atoms);
	HaAtom* aptr;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		int rad = (int)(0.5*Scale) + 1; 
		HaColor pick_col(255, 255, 0);

		pCanv->ClipSphere( aptr->x,aptr->y,aptr->z, rad, pick_col.cidx );
	}
}

void HaMolView::DisplayOnScreenInfo()
{
	int x = 20;
	int y = pCanv->YRange() - pCanv->m_FontSize - 4;
	int z = pCanv->m_ZOffset;

	int i;
	MolSet* pmset = GetMolSet();
	int nstr = pmset->info_str.size();

	for (i = 0; i < nstr; i++)
	{
		pCanv->DisplayTextString(x,y,z,(pmset->info_str[i]).c_str(),LabelColor.cidx);
		y -= (pCanv->m_FontSize + 4);
	}
}

void HaMolView::RenderFrame()
//! Function to create image by filling View Structure 
//! Enter plotting of new elements here 
{
//	std::cerr << std::endl << " HaMolView::RenderFrame() pt 1 " << std::endl;
	MolSet* pmset = this->GetMolSet();
    HaChain *chain;
	
	MoleculesType::iterator mol_itr;
	
	if( DrawAtoms ) 
		DisplaySpaceFill();
	
	if( !UseSlabPlane() || (SlabMode() != SlabSection) )
	{   
		if( DrawBonds ) 
		{   
			if( DrawDoubleBonds )
			{   
				DisplayDoubleBonds();
			} 
			else 
				DisplayWireframe();
		}
		
		if( DrawRibbon )
		{
			ChainIteratorMolSet chitr( this->GetMolSet() );
			for( chain = chitr.GetFirstChain(); chain; chain = chitr.GetNextChain() )
			{
				if( !chain->res_map.empty() ) DisplayRibbon( chain );
			}
		}
		
		if( DrawDots ) DisplayDotSurfaces();
		if( DrawLabels ) DisplayLabels();
		if( !MonitList.empty() ) DisplayMonitors();
		if( DrawBestPath ) DisplayETBestPath();
		if( DrawContourSurf) DisplayContourSurf();
        if( DrawObj3D) DisplayObj3D();//<mikola 30July06
	    DisplaySSBonds();
    	DisplayHBonds();
		DisplayBackbone();
	}
    DisplayBoxes();
	DisplayOnScreenInfo();
	DisplayPickedAtoms();
	HaFlexMod* p_flex_mod = pmset->GetFlexMod(false);
	if(p_flex_mod) p_flex_mod->Display(this);
}


void HaMolView::DrawFrame()
//! Fill FBuffer and DBuffer that determine the image
//! Use RenderFrame and changing View Structure to plot Stereo 
//! this mechanism may be used to plot multiple molecule on the same screen 
//! for the bulder for example?
{
//	std::cerr << " HaMolView::DrawFrame() pt 1 " << std::endl;
    double temp;
    int wide;
	
    if( GetMolSet()->HostMolecules.empty()  && GetMolSet()->ViewObjects.size() == 0 ) return;
	
    ClearBuffers();
		
    if( UseSlabPlane() )
    {   
		SetSlabValue( (int)(CurSlabValue*GetImageRadius() )+ZOffset() );
		SetSlabInten( (int)(ColourMask*LightZComp/LightLength) );
		SetSliceValue( SlabValue() +16 );
		UseClipping = True;
    } 
	else 
		UseClipping = UseScreenClip;
	
    /* Common View Elements */
    pCanv->View.yskip = pCanv->XRange();
    pCanv->View.ymax =  pCanv->YRange();
	
    if( UseStereo )
    {   
		temp = StereoAngle/180.0;
        wide = pCanv->XRange()>>1;
		
        /* Create 'Left' View structure */
        pCanv->View.fbuf = pCanv->FBuffer;
        pCanv->View.dbuf = pCanv->DBuffer;
        pCanv->View.xmax = wide;
		
        CurRY -= temp;
        ReDrawFlag |= RFRotateY;
        ApplyTransform();
        RenderFrame();        
		
        /* Create 'Right' View structure */
        pCanv->View.fbuf = pCanv->FBuffer+wide;
        pCanv->View.dbuf = pCanv->DBuffer+wide;
        pCanv->View.xmax = wide;
		
        CurRY += temp;
        ReDrawFlag |= RFRotateY;
        ApplyTransform();
        RenderFrame();       
		
    } 
	else /* Mono */
    {   /* Create 'Mono' View structure */
        pCanv->View.fbuf = pCanv->FBuffer;
        pCanv->View.dbuf = pCanv->DBuffer;
        pCanv->View.xmax = pCanv->XRange();
        RenderFrame();
    }
	
    pCanv->DBClear = False;
    pCanv->FBClear = False;

}

void HaMolView::WriteImageFile(const char* name, int type )
{
    if( !type )
        type = PICTTok;
//	type = PPMTok;
	
    switch( type )
    {   
	case(GIFTok):     WriteGIFFile(name);             break;
	case(JPEGTok):     WriteJPEGFile(name);            break;
	case(BMPTok):     WriteBMPFile(name);             break;
	case(PPMTok):     WritePPMFile(name,True);        break;
	case(PICTTok):    WritePICTFile(name);            break;
	case(IRISTok):    WriteIRISFile(name);            break;
	case(VectPSTok):  WriteVectPSFile(name);          break;
		
	case(RasMolTok):
	case(ScriptTok):     WriteScriptFile(name);     break;
	case(MolScriptTok):  WriteMolScriptFile(name);  break;
	case(POVRayTok):     WritePOVRayFile(name);     break;
	case(VRMLTok):       WriteVRMLFile(name);       break;
    }
}


void 
HaMolView::TestAtomProximity( HaAtom* ptr, int xpos, int ypos )
{
    int dist;
    int dx,dy;
	
    if( UseSlabPlane() && (ptr->z > SlabValue()) )
		return;
	
    dx = ptr->x - xpos;
    dy = ptr->y - ypos;
	
    dist = (int)dx*dx + (int)dy*dy;
	
    if( IdentFound )
    {   
		if( dist==IdentDist )
		{   
			if( ptr->z<IdentDepth )
				return;
		} 
		else if( dist>IdentDist ) 
			return;
    }
	
    IdentDepth = ptr->z;
    IdentFound = True;
    IdentDist = dist;
    PkAtom = ptr;
}


void HaMolView::IdentifyAtom( int xpos, int ypos )
{
    int rad, wide, dpth;
    int newf, dx, dy, dz;
    set<HaHBond, less<HaHBond> >::iterator hptr;
	HaChain* chain;
    HaAtom   *aptr;
    HaBond  *bptr;
	
    /* Reset Search */
    PkAtom =  NULL;
    IdentFound = False;

	MoleculesType::iterator mol_itr;

	MolSet* pmset = GetMolSet();

    BondIteratorMolSet bitr(pmset);
    
	if( !UseSlabPlane() || (SlabMode() != SlabSection) )
	{   
		if( DrawBonds )
		{
			AtomIteratorMolSet aitr(pmset);
			for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
			{
				std::vector<HaBond*>& bonds = aptr->GetBonds();
				if( bonds.empty() && aptr->IsDisplayed() )
				{
					TestAtomProximity(aptr, xpos, ypos );
					continue;
				}
				std::vector<HaBond*>::iterator bitr_at = bonds.begin();
				for(; bitr_at != bonds.end(); bitr_at++)
				{
					if( (*bitr_at)->IsToDraw() )
					{
						TestAtomProximity(aptr, xpos, ypos );
						break;
					}
				}
			}

			for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond())
			{
				if( bptr->IsToDraw() )
				{   
					TestAtomProximity(bptr->srcatom,xpos,ypos);
					TestAtomProximity(bptr->dstatom,xpos,ypos);
				}
			}
			
			pmset->GetMolEditor()->UpdateBackBone(pmset);
			int nb = pmset->BackboneBonds.size();
			int ib;
			for( ib = 0; ib < nb; ib++ )
			{
				bptr = pmset->BackboneBonds[ib];
				if( bptr->IsToDraw() )
				{   
					TestAtomProximity(bptr->srcatom,xpos,ypos);
					TestAtomProximity(bptr->dstatom,xpos,ypos);
				}
			}
			
			for( hptr= pmset->HBonds.begin(); hptr != pmset->HBonds.end() ; hptr++ )
			{
				if( (*hptr).IsToDraw() )
				{   
					if( HBondMode )
					{   
						TestAtomProximity((*hptr).srcCA,xpos,ypos);
						TestAtomProximity((*hptr).dstCA,xpos,ypos);
					} 
					else
					{   
						TestAtomProximity((*hptr).src,xpos,ypos);
						TestAtomProximity((*hptr).dst,xpos,ypos);
					}
				}
			}
		}	
	}
	
	AtomIteratorMolSet aitr(this->GetMolSet());

	IdentDist = 3;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{   
		if( aptr->IsDrawSphere() )
		{   
			dy = AbsFun(aptr->y-ypos);
			if( dy>aptr->irad ) continue;
			rad = pCanv->LookUp[aptr->irad][dy];
			dx = AbsFun(aptr->x-xpos);
			if( dx>rad ) continue;
			
			newf = False;
			dpth = aptr->z+ pCanv->LookUp[rad][dx];
			if( UseSlabPlane() && (aptr->z+rad >= SlabValue()) )
			{   
				dz = SlabValue() - aptr->z;
				if( SlabMode() && (dz >= -rad) )
				{   
					wide = pCanv->LookUp[aptr->irad][AbsFun(dz)];
					if( (dy<=wide) && (dx<=(int)(pCanv->LookUp[wide][dy])) )
					{   
						if( SlabMode() == SlabFinal )
						{   
							dpth = SliceValue();
							newf = True;
						} 
						else if( SlabMode() == SlabHollow )
						{   
							dpth = aptr->z- pCanv->LookUp[rad][dx];
							newf = !IdentFound || (dpth>IdentDepth);
						} 
						else if( SlabMode() != SlabHalf )
						{   /* SlabClose, SlabSection */
							dpth = dx*dx+dy*dy+dz*dz+SliceValue();
							if( IdentFound )
							{  
								newf = (IdentDepth<SliceValue()) || (dpth<IdentDepth);
							} 
							else 
								newf=True;
						}
					} 
					else if( (dz>0) && (SlabMode() != SlabSection) )
						newf = !IdentFound || (dpth>IdentDepth);
				}
			} 
			else if( !UseSlabPlane() || (SlabMode() != SlabSection) )
				newf = !IdentFound || IdentDist || (dpth>IdentDepth);
			
			if( newf )
			{   
				IdentFound = True;
				IdentDepth = dpth;
				IdentDist = 0;
				
				PkAtom = aptr;
			}
		} 
	}
	
	
    if( !IdentFound || (IdentDist>=50) )
    {   /* Reset Pick Atom! */
 		PkAtom =      NULL;
    }
}


void HaMolView::SetPickMode( int mode )
{
	if ( mode == PickNone ) PrintLog(" Turn off Atom Pick Mode \n");
	else if ( mode == PickIdent ) PrintLog(" Set Show Atom ID Mode \n");
	else if ( mode == PickDist ) PrintLog(" Set Measure Atom Distance Mode \n");
	else if ( mode == PickAngle ) PrintLog(" Set Measure Valence Angle Mode \n");
	else if ( mode == PickTorsn ) PrintLog(" Set Measure Torsion Angle Mode \n");
	else if ( mode == PickLabel ) PrintLog(" Set Show Atom Label Mode \n");
	else if ( mode == PickMonit ) PrintLog(" Set Add Distance Monitor Mode \n");
	else if ( mode == PickCentr ) PrintLog(" Set Molecular Center Mode \n");
	else if ( mode == PickMolConnect ) PrintLog(" Set Molecule Connect Mode \n");

    PickMode = mode;
    PickCount = 0;
}



void HaMolView::PickAtom( int shift, int xpos, int ypos )
{
    float temp;
    char *str;
    int len;
	std::string label;

    char buffer[256];
	char buf2[256];
	int j;
	
	
    if( PickMode == PickNone )
        return;
    
//	cerr << " HaMolView::PickAtom() " << endl;
//	cerr << " xpos = " << xpos << "  ypos= " << ypos << endl;

		IdentifyAtom(xpos,ypos);
		
		if( PkAtom )
		{
		if( PickMode == PickIdent )
		{
			PrintLog("\n");
			
			BroadcastCurrAtom();
			
			PkAtom->FillRef(buffer);
			PrintMessage(buffer);			
			PrintLog("\tCoord.[A,A,A]=[%f,%f,%f]\n",PkAtom->GetX_Ang(),PkAtom->GetY_Ang(),PkAtom->GetZ_Ang());
		} 
		else if( PickMode == PickLabel )
        {   
			if( PkAtom->label.empty() )
            {   
				if( (PkAtom->GetHostMol())->GetNRes() > 1 )
                {   
					strcpy(buffer,"%n%r");
                    str = buffer+4;
                    if( (PkAtom->GetHostMol())->GetNChains() > 1 )
                    {   
						if( isdigit((PkAtom->GetHostChain())->ident) )
                            *str++ = ':';
                        *str++ = '%';
                        *str++ = 'c';
                    }
                    strcpy(str,".%a");
					
                    len = (str-buffer) + 3;
					std::string tmp_str(buffer,len);
                    label = tmp_str;
                } 
				else 
					label = "%e%i";
				
                PkAtom->label = label;
            } 
			else
            {   
                PkAtom->label = "";
            }
			
            DrawLabels = True;
            ReDrawFlag |= RFRefresh;
			
        } 
		else if( PickMode == PickCentr )
        {   
			CenX = PkAtom->GetX();
            CenY = PkAtom->GetY();
            CenZ = PkAtom->GetZ();
									
			PkAtom->FillRef(buf2);
            PrintLog("\n Rotating about %s ",buf2);
			
        } 
		else if( PickMode == PickMonit )
        {   /* State Machine Implementation */
			
            if( PickCount == 0 )
            {   
				PickHist[0] = PkAtom;
                PickCount = 1;
            } 
			else if( PickCount == 1 )
            {   
				if( !shift )
                {   
					if( PickHist[0] != PkAtom )
                    {   
						AddAtomPairMonitor(PickHist[0],PkAtom);
                        ReDrawFlag |= RFRefresh;
                    }
                    PickCount = 2;
                } 
				else 
					PickHist[0] = PkAtom;
            } 
			else /* PickCount == 2 */
                if( !shift )
                {   
					PickHist[0] = PkAtom;
                    PickCount = 1;
                } 
				else if( PickHist[0] != PkAtom )   
                {   
					AddAtomPairMonitor(PickHist[0],PkAtom);
                    ReDrawFlag |= RFRefresh;
                }
				
        } 
		else if(PickMode == PickMolConnect)  // Pick a molecule to connect
		{
			if(PkAtom != NULL) ConnectObject(PkAtom->GetHostMol());
			PickMode=PickIdent;

		}
		else /* Distance, Angle or Torsion! */
        {   
			if( PickCount )
            {   
				if( shift )
                {   
					PickCount--;
                } 
				else if( PickCount == PickMode )
                    PickCount = 0;
            }
			
            PickHist[PickCount] = PkAtom;
            PickCount++;
			
			PrintLog("\n");
			
		    j= sprintf(buffer,"%s","Atom #");
			j+= sprintf(buffer+j,"%c:",(char)(PickCount+'0'));
		    PkAtom->FillRef(buf2);
			j+= sprintf(buffer+j,"%s",buf2);
            PrintMessage(buffer);
			
            if( PickCount == PickMode )
            {   
				if( PickMode == PickDist )
                {   
					temp = HaAtom::CalcDistance(PickHist[0],PickHist[1], ANGSTROM_U);
					
					j= sprintf(buffer,"%s","Distance ");
                    PickHist[0]->FillRef(buf2);
					j+= sprintf(buffer +j ,"%s - ",buf2);
					PickHist[1]->FillRef(buf2);
					j+= sprintf(buffer +j ,"%s  =",buf2);                
                    j+= sprintf(buffer+j, " %.3f",temp);
					PrintMessage(buffer);
					
                } 
				else if( PickMode == PickAngle )
                {   
					temp = RAD_TO_DEG * HaAtom::CalcAngle(PickHist[0],PickHist[1],PickHist[2]);
					
					j= sprintf(buffer,"%s","Angle ");
                    PickHist[0]->FillRef(buf2);
					j+= sprintf(buffer +j ,"%s - ",buf2);
					PickHist[1]->FillRef(buf2);
					j+= sprintf(buffer +j ,"%s  -",buf2);                
					PickHist[2]->FillRef(buf2);
					j+= sprintf(buffer +j ,"%s  =",buf2);                
                    j+= sprintf(buffer+j, " %.1f",temp);
                    PrintMessage(buffer);
					
                } 
				else /* PickMode == PickTorsn */
                {   
					temp = RAD_TO_DEG * HaAtom::CalcTorsion(PickHist[0],PickHist[1],PickHist[2],PickHist[3]);
					
					j= sprintf(buffer,"%s","Torsion ");
                    PickHist[0]->FillRef(buf2);
					j+= sprintf(buffer +j ,"%s - ",buf2);
					PickHist[1]->FillRef(buf2);
					j+= sprintf(buffer +j ,"%s  -",buf2);                
					PickHist[2]->FillRef(buf2);
					j+= sprintf(buffer +j ,"%s  -",buf2);                
					PickHist[3]->FillRef(buf2);
					j+= sprintf(buffer +j ,"%s  =",buf2);                
                    j+= sprintf(buffer+j, " %.1f",temp);
                    PrintMessage(buffer);
                 }
			}
        }
    }
	this->UpdateThisView(ReDrawFlag);
}

void HaMolView::ClampShiftVal(int  ivar,double  value )
{
    double temp;
	
	double* ptr_val = NULL;

	if( ivar == 3) ptr_val = &CurZoom;
	if( ivar == 4) ptr_val = &CurTX;
	if( ivar == 5) ptr_val = &CurTY;
	if( ivar == 7) ptr_val = &CurSlabValue;

	if(ptr_val == NULL) return;

    temp = *ptr_val + value;

    if( temp > 1.0 )
    {   
		*ptr_val = 1.0;
    } 
	else if( temp < -1.0 )
    {   
		*ptr_val = -1.0;
    } 
	else 
		*ptr_val = temp;
}

void 
HaMolView::SetStereoMode( int enable )
{
    ReDrawFlag |= RFRefresh | RFTransX;
    StereoView = ViewLeft;
    UseStereo = enable;
}


void HaMolView::ResetRenderer()
{
    DrawAtoms = False; 
    DrawBonds = False;  
    DrawRibbon = False; DrawDots = False;

    SetSlabMode(SlabClose);
    SetUseSlabPlane(False);
    UseLabelCol = False;

    UseDepthCue = False;

    SSBondMode = False;
    HBondMode = False;

    DrawDoubleBonds = False;
    DrawBoundBox = False;
    DrawUnitCell = False;
    DrawAxes = False;

    SetStereoMode(False);
    StereoAngle = 6.0;

	UseMolShift = FALSE;
	mol_shift = 8.0;
}


void HaMolView::InitializeTables()
{
    unsigned char  *ptr;
    unsigned int root,root2;
    unsigned int i,rad,arg;

    ptr = pCanv->Array;
    pCanv->LookUp[0] = ptr;  *ptr++ = 0;
    pCanv->LookUp[1] = ptr;  *ptr++ = 1;  *ptr++ = 0;
    
    for( rad=2; rad<MAXRAD; rad++ )
    {   
		pCanv->LookUp[rad] = ptr;
		
        /* i == 0 */
        *ptr++ = (unsigned char)rad;  
		
        root = rad-1;
        root2 = root*root;
		
        arg = rad*rad;
		for( i=1; i<rad; i++ )
        {   /* arg = rad*rad - i*i | correct IK */
            arg -= (i<<1)-1;
			
            /* root = isqrt(arg)   */
            while( arg < root2 )
            {   root2 -= (root<<1)-1;
			root--;
            }
            /* Thanks to James Crook */
            *ptr++ = ((arg-root2)<i)? root : root+1;
        }
		
        /* i == rad */
        *ptr++ = 0;    
    }
}


void 
HaMolView::InitializeRenderer()
{
    int rad,maxval;
		
	Canvas3D::ColConst = Canvas3D::ColConstTable;
	
    InitializeTables();
	
    /* Initialize ColConst! */
    for( rad=0; rad<MAXRAD; rad++ )
    {   
		maxval = (int)(LightLength*rad)+4;
		Canvas3D::ColConst[rad] = ((unsigned int)ColourDepth<<ColBits)/maxval;
    }
	
    PickMode = PickIdent;
	
    ResetRenderer();
    ReSizeScreen();
}

int 
HaMolView::FillCurrAtomRef(char* buf)
{
 	return PkAtom->FillRef(buf);
}

static void CentreZoneExpr(AtomExpr* expr)
{
    register double x, y, z;
    register int count;
	HaAtom* aptr;

	MolSet* pmset = GetCurMolSet();

    if( !pmset) return;

    count = 0;
    x = y = z = 0.0;
    
	AtomIteratorMolSet aitr(pmset);

	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if( expr->EvaluateExprFor(aptr) )
		{   
			x += (double)aptr->GetX();
			y += (double)aptr->GetY();
			z += (double)aptr->GetZ();
			count++;
		}
	}

    if( count )
    {   
		CurMolView->CenX = x/count;
        CurMolView->CenY = y/count;
        CurMolView->CenZ = z/count;
    } 
	else
    {   
        PrintLog("No Atoms to center!\n");
    }
}

int HaMolView::ExecuteCommand(CmdParser& cmd_pr)
{
	int option;
	int temp;
	int done;
	int RVal, GVal, BVal;

	MolSet* pmset = GetMolSet();

	cmd_pr.ResetCursorPosition();
	
    if( !cmd_pr.FetchToken() )  return FALSE;
 
	switch( cmd_pr.CurToken )
    {   
	case(ColourTok):  
		ExecuteColourCommand(cmd_pr);
		break;
	case(SetTok):
		ExecuteSetCommand(cmd_pr);
        break;

	case(WireframeTok):
			cmd_pr.FetchToken();
			if( cmd_pr.CurToken==FalseTok )
			{   
				ReDrawFlag |= RFRefresh;
				DisableWireframe();
			} 
			else if( (cmd_pr.CurToken==TrueTok) || !cmd_pr.CurToken )
			{   
				ReDrawFlag |= RFRefresh;
				EnableWireframe(WireFlag,0.0);
			} 
			else if( cmd_pr.CurToken==DashTok )
			{   
				ReDrawFlag |= RFRefresh;
				EnableWireframe(DashFlag,0.0);
			} 
			else if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken == FloatTok)
			{   				
				if( cmd_pr.TokenValueFloat <= 2.0 )
				{   
					EnableWireframe(CylinderFlag, cmd_pr.TokenValueFloat);
					ReDrawFlag |= RFRefresh;
				} 
				else 
					PrintLog("Parameter value too large");
			} 
			else 
				PrintLog("Invalid command argument\n");
			break;

	case(BackboneTok):
		cmd_pr.FetchToken();
		if( cmd_pr.CurToken==FalseTok )
		{   
			ReDrawFlag |= RFRefresh;
			DisableBackbone();
		} 
		else if( (cmd_pr.CurToken==TrueTok) || !cmd_pr.CurToken )
		{   
			ReDrawFlag |= RFRefresh;
			EnableBackbone(WireFlag,0);
		} 
		else if( cmd_pr.CurToken==DashTok )
        {   
			ReDrawFlag |= RFRefresh;
            EnableBackbone(DashFlag,0);
        } 
		else if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken == FloatTok)
        {   

            if( cmd_pr.TokenValueFloat < 2.0 )
            {   
				EnableBackbone(CylinderFlag, cmd_pr.TokenValueFloat);
                ReDrawFlag |= RFRefresh;
            } 
			else PrintLog("Parameter value too large");
        } 
		else 
			PrintLog("Invalid command argument\n");
        break;
        case(CPKTok):
        case(SpacefillTok):
                          cmd_pr.FetchToken();
                          if( cmd_pr.CurToken==FalseTok )
                          {   
							  ReDrawFlag |= RFRefresh;
                              DisableSpacefill();
                          } 
						  else if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken == FloatTok)
                          {   

                              if( cmd_pr.TokenValueFloat <= 3.0 )
                              {   
								  SetAtomScreenRadVal(MaxFun( cmd_pr.TokenValueFloat, 1.0));
                                  ReDrawFlag |= RFRefresh;
                              } 
							  else PrintLog("Parameter value too large");
                          } 
						  else if( cmd_pr.CurToken==TemperatureTok )
                          {   
							  ReDrawFlag |= RFRefresh;
                              SetRadiusTemperature();
                          } 
						  else if( (cmd_pr.CurToken==TrueTok) || !cmd_pr.CurToken )
                          {   
							  ReDrawFlag |= RFRefresh;
                              SetAtomScreenRadVdW();
                          } 
						  else 
							  PrintLog("Invalid command argument\n");
                          break;

        case(DashTok): 
			cmd_pr.FetchToken();
			if( cmd_pr.CurToken==FalseTok )
			{   
				ReDrawFlag |= RFRefresh;
				DisableWireframe();
			} 
			else if( (cmd_pr.CurToken==TrueTok) || !cmd_pr.CurToken )
			{   
				ReDrawFlag |= RFRefresh;
				EnableWireframe(DashFlag,0);
			} 
			else 
				PrintLog("Invalid command argument\n");
			break;

        case(SSBondTok):  
			cmd_pr.FetchToken();
			if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken == FloatTok)
			{   
				if( cmd_pr.TokenValueFloat <= 2.0 )
				{   
					SetSSBondStatus(True,cmd_pr.TokenValueFloat);
					ReDrawFlag |= RFRefresh;
				} 
				else PrintLog("Parameter value too large");
			} 
			else if( cmd_pr.CurToken==FalseTok )
			{   
				ReDrawFlag |= RFRefresh;
				SetSSBondStatus(False,0);
			} 
			else if( (cmd_pr.CurToken==TrueTok) || !cmd_pr.CurToken )
			{   
				ReDrawFlag |= RFRefresh;
				SetSSBondStatus(True,0);
			} 
			else 
				PrintLog("Invalid command argument\n");
			break;

        case(HBondTok):   
			cmd_pr.FetchToken();
			if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken == FloatTok)
			{   
				if( cmd_pr.TokenValueFloat <= 2.0 )
				{   
					SetHBondStatus(True, cmd_pr.TokenValueFloat);
					ReDrawFlag |= RFRefresh;
				} 
				else 
					PrintLog("Parameter value too large");
			} 
			else if( cmd_pr.CurToken==FalseTok )
			{   
				ReDrawFlag |= RFRefresh;
				SetHBondStatus(False,0);
			} 
			else if( (cmd_pr.CurToken==TrueTok) || !cmd_pr.CurToken )
			{   
				ReDrawFlag |= RFRefresh;
				SetHBondStatus(True,0);
			} 
			else PrintLog("Invalid command argument\n");
			break;

        case(RibbonTok): 
			cmd_pr.FetchToken();
			if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken==FloatTok)
			{   
				if( cmd_pr.TokenValueFloat <= 4.0 )
				{   
					SetRibbonStatus(True,RibbonFlag,cmd_pr.TokenValueFloat);
					ReDrawFlag |= RFRefresh;
				} 
				else 
					PrintLog("Parameter value too large");
			} 
			else if( cmd_pr.CurToken==FalseTok )
			{   
				ReDrawFlag |= RFRefresh;
				SetRibbonStatus(False,RibbonFlag,0);
			} 
			else if( (cmd_pr.CurToken==TrueTok) || !cmd_pr.CurToken )
			{   
				ReDrawFlag |= RFRefresh;
				SetRibbonStatus(True,RibbonFlag,0);
			} 
			else 
				PrintLog("Invalid command argument\n");
			break;

        case(StrandsTok): 
			cmd_pr.FetchToken();
			if( cmd_pr.CurToken == DashTok )
			{   
				option = DashStrandFlag;
				cmd_pr.FetchToken();
			} 
			else 
				option = StrandFlag;
			
			if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken == FloatTok)
			{   
				if( cmd_pr.TokenValueFloat <= 4.0 )
				{   
					SetRibbonStatus(True,option, cmd_pr.TokenValueFloat);
					ReDrawFlag |= RFRefresh;
				} 
				else 
				    PrintLog("Parameter value too large");
			} 
			else if( cmd_pr.CurToken==FalseTok )
			{   
				ReDrawFlag |= RFRefresh;
				SetRibbonStatus(False,option,0);
			} 
			else if( (cmd_pr.CurToken==TrueTok) || !cmd_pr.CurToken )
			{   
				ReDrawFlag |= RFRefresh;
				SetRibbonStatus(True,option,0);
			} 
			else 
				PrintLog("Invalid command argument\n");
			break;

        case(TraceTok):
			cmd_pr.FetchToken();
                          if( cmd_pr.CurToken==FalseTok )
                          {   
							  ReDrawFlag |= RFRefresh;
                              SetRibbonStatus(False,TraceFlag,0.32);
                          } 
						  else if( (cmd_pr.CurToken==TrueTok) || !cmd_pr.CurToken )
                          {   
							  ReDrawFlag |= RFRefresh;
                              SetRibbonStatus(True,TraceFlag,0.32);
                          } 
						  else if( cmd_pr.CurToken==TemperatureTok )
                          {   
							  ReDrawFlag |= RFRefresh;
                              SetTraceTemperature();
                          } 
						  else if( cmd_pr.CurToken==DotsTok )
                          {   
							  ReDrawFlag |= RFRefresh;
                              SetRibbonStatus(True,DotsFlag,0.32);
                          } 
						  else if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken == FloatTok)
                          {   
                              if( cmd_pr.TokenValueFloat <= 2.0 )
                              {   
								  SetRibbonStatus(True,TraceFlag, cmd_pr.TokenValueFloat);
                                  ReDrawFlag |= RFRefresh;
                              } 
							  else 
								  PrintLog("Parameter value too large");
                          } 
						  else 
							  PrintLog("Invalid command argument\n");
                          break;

        case(CartoonTok): 
			cmd_pr.FetchToken();
			if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken == FloatTok)
			{   
				if( cmd_pr.TokenValue <= 4.0 )
				{   
					SetRibbonStatus(True,CartoonFlag, cmd_pr.TokenValueFloat);
					ReDrawFlag |= RFRefresh;
				} 
				else 
					PrintLog("Parameter value too large");
			} 
			else if( cmd_pr.CurToken==FalseTok )
			{   
				ReDrawFlag |= RFRefresh;
				SetRibbonStatus(False,CartoonFlag,0);
			} 
			else if( (cmd_pr.CurToken==TrueTok) || !cmd_pr.CurToken )
			{   
				ReDrawFlag |= RFRefresh;
				SetRibbonStatus(True,CartoonFlag,0);
			} 
			else 
				PrintLog("Invalid command argument\n");
			break;
			
        case(DotsTok):  
			cmd_pr.FetchToken();
			pApp->StartWait();
			if( cmd_pr.CurToken==NumberTok )
			{   
				if( cmd_pr.TokenValue<=1000 )
				{   
					if( cmd_pr.TokenValue )
					{   
						CalculateDotSurface((int)cmd_pr.TokenValue);
					} 
					else 
					{
						CalculateDotSurface(1);
					}
					ReDrawFlag |= RFRefresh;
				} 
				else 
					PrintLog("Parameter value too large");
			} 
			else if( cmd_pr.CurToken==FalseTok )
			{   
				ReDrawFlag |= RFRefresh;
				DeleteDotSurfaces();
			} 
			else if( (cmd_pr.CurToken==TrueTok) || !cmd_pr.CurToken )
			{   
				ReDrawFlag |= RFRefresh;
				CalculateDotSurface(100);
			} 
			else 
				PrintLog("Invalid command argument\n");
			pApp->EndWait();
			break;

        case(MonitorTok): 
			cmd_pr.FetchToken();
			if( cmd_pr.CurToken == NumberTok )
			{   
				temp = cmd_pr.TokenValue;
				cmd_pr.FetchToken();
				if( cmd_pr.CurToken == ',' )
					cmd_pr.FetchToken();
				
				if( cmd_pr.CurToken == NumberTok )
				{   
					CreateMonitor(temp,cmd_pr.TokenValue);
					ReDrawFlag |= RFRefresh;
				} 
				else PrintLog("Integer value expected\n");;
			} 
			else if( cmd_pr.CurToken == StringTok)
			{
				std::string str1,str2;
				str1 = cmd_pr.TokenIdent;
				cmd_pr.FetchToken();
				if( cmd_pr.CurToken == StringTok)
				{
					str2 = cmd_pr.TokenIdent;
					MolSet* pmset = GetMolSet();
					HaAtom* aptr1 = pmset->GetAtomByRef(str1.c_str());
					HaAtom* aptr2 = pmset->GetAtomByRef(str2.c_str());
					if(aptr1 != NULL && aptr2 != NULL)
					{
						AddAtomPairMonitor( aptr1, aptr2 );
					}
				}
				else
				{
					PrintLog("Integer value expected\n");
				}
			}
			else if( cmd_pr.CurToken == FalseTok )
			{   
				ReDrawFlag |= RFRefresh;
				DeleteAllMonitors();
			} 
			else 
				PrintLog("Invalid command argument\n");
			break;

        case(SlabTok): 
			cmd_pr.FetchToken();
			if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken== FloatTok )
			{   
				
				if( cmd_pr.TokenValueFloat <= 100.0 && cmd_pr.TokenValueFloat > 0.0)
				{   
					CurSlabValue = (cmd_pr.TokenValueFloat - 50.0)/50.0;
					ReDrawFlag |= RFSlab;
					SetUseSlabPlane(True);
				} 
				else 
					PrintLog("Parameter value too large");
				
			} 
			else if( cmd_pr.CurToken==FalseTok )
			{   
				if( UseSlabPlane() )
				{   
					ReDrawFlag |= RFRefresh;
					SetUseSlabPlane(False);
				}
			} 
			else if( !cmd_pr.CurToken || (cmd_pr.CurToken==TrueTok) )
			{   
				if( !UseSlabPlane() )
				{   
					ReDrawFlag |= RFRefresh;
					SetUseSlabPlane(True);
				}
			} 
			else 
				PrintLog("Invalid command syntax\n");
			break;

        case(ZoomTok):  
			cmd_pr.FetchToken();
			if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken == FloatTok )
			{   				
				if( cmd_pr.TokenValueFloat > 0.0 )
				{   
					Zoom = cmd_pr.TokenValueFloat/100.0;
					ReDrawFlag |= RFZoom;
				} 
			} 
			else if( cmd_pr.CurToken==TrueTok )
			{   
				ReDrawFlag |= RFZoom;
				Zoom = 1.5;
			} 
			else if( !cmd_pr.CurToken || (cmd_pr.CurToken==FalseTok) )
			{   
				ReDrawFlag |= RFZoom;
				Zoom = 1.0;
			} 
			else 
				PrintLog("Invalid command syntax\n");
			break;

        case(RotateTok):  
			cmd_pr.FetchToken();
			if( cmd_pr.CurToken==XTok )
			{   
				option = 0;
			} 
			else if( cmd_pr.CurToken==YTok )
			{   
				option = 1;
			} 
			else if( cmd_pr.CurToken==ZTok )
			{   
				option = 2;
			} 
			else
			{   
				PrintLog("Invalid command syntax\n");
				break;
			}
			
			cmd_pr.FetchToken();
			if( (done=(cmd_pr.CurToken=='-')) )
				cmd_pr.FetchToken();

			if( option != 1 )
				done = !done;

			if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken== FloatTok )
			{   				
				int redraw_flag = (1<<option);
				if( done ) cmd_pr.TokenValueFloat = -cmd_pr.TokenValueFloat;

				WrapShiftVal(option,cmd_pr.TokenValueFloat/180.0);
				MolSet* pmset = GetMolSet();
				pmset->RefreshAllViews(redraw_flag);	
			} 
			else 
				PrintLog("Integer or Float value is expected\n");;
			break;

        case(TranslateTok):
			cmd_pr.FetchToken();
			if( cmd_pr.CurToken==XTok )
			{   
				option = 4;
			} 
			else if( cmd_pr.CurToken==YTok )
			{   
				option = 5;
			} 
			else if( cmd_pr.CurToken==ZTok )
			{   
				option = 6;
			} 
			else
			{   
			    PrintLog("Invalid command syntax\n");
				break;
			}
			
			cmd_pr.FetchToken();
			if( (done=(cmd_pr.CurToken=='-')) )
				cmd_pr.FetchToken();

			if( option == 5 )
				done = !done;
			
			if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken== FloatTok )
			{   				
				if( cmd_pr.TokenValue <= 100.00 )
				{   
					if( done ) cmd_pr.TokenValueFloat = -cmd_pr.TokenValueFloat;
					
					MoleculesType::iterator mol_itr;
					ForEachMol_VIEW
					{
						HaMolecule* pMol = (*mol_itr);
						if( pMol->IsConnected())
						{
							HaAtom* aptr;
							AtomIteratorMolecule aitr(pMol);
							for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
							{
								if(option == 4)
								{
									double x=aptr->GetX() + cmd_pr.TokenValueFloat;
									aptr->SetX(x);
								}
								else if(option == 5)
								{
									double y=aptr->GetY() + cmd_pr.TokenValueFloat;
									aptr->SetY(y);
								}
								else if(option == 6)
								{
									double z=aptr->GetZ() + cmd_pr.TokenValueFloat;
									aptr->SetZ(z);
								}
							}		
						}
					}
					(GetMolSet())->RefreshAllViews(RFApply | RFRefresh);
				} 
				else PrintLog("Parameter value too large");
			} 
			else PrintLog("Integer value expected\n");;
			break;
			
        case(StereoTok):
			cmd_pr.FetchToken();
			if( !cmd_pr.CurToken || (cmd_pr.CurToken==TrueTok) )
			{   
				SetStereoMode(True);
			} 
			else if( cmd_pr.CurToken==FalseTok )
			{   
				SetStereoMode(False);
			} 
			else if( cmd_pr.CurToken == '-' )
			{   
				if( !cmd_pr.NextIf(NumberTok,"Integer value expected") )
				{   
					StereoAngle = -cmd_pr.TokenValue;
					SetStereoMode(True);
				}
			} 
			else if( cmd_pr.CurToken==NumberTok )
			{   
				StereoAngle = cmd_pr.TokenValue;
				SetStereoMode(True);
			} 
			else 
				PrintLog("Invalid command syntax\n");
			break;
			
        case(CentreTok):  
			cmd_pr.FetchToken();
			if( !cmd_pr.CurToken || (cmd_pr.CurToken==AllTok) )
			{   
				CenX = CenY = CenZ = 0.0;
			} 
			else
			{
				AtomExpr* p_expr;
				if( (p_expr = cmd_pr.ParseExpression(0,pmset)) != NULL )
				{   
					if( !cmd_pr.CurToken )
					{   
						CentreZoneExpr(p_expr);
					} 
					else 
					{
						PrintLog("Invalid command syntax\n");
					}
					delete p_expr;
				}
			}
			break;
				
        case(ResizeTok):  cmd_pr.FetchToken();
                          break;

        case(ResetTok):
			CurRX = CurRY = CurRZ = 0.0;
			CurTX = CurTY = 0.0;
			CurZoom = CurSlabValue = 0.0;
			ReDrawFlag |= RFDials;
			ResetTransform();
			UpdateThisView(0);
			break;

        case(LabelTok):
			cmd_pr.FetchToken();
			if( !cmd_pr.CurToken || (cmd_pr.CurToken==TrueTok) )
			{   
				if( (GetMolSet())->GetNChains() >1 )
				{   
					DefineLabels("%n%r:%c.%a");
				} 
				else if( (GetMolSet())->GetNRes() > 1 )
				{   
					DefineLabels("%n%r.%a");
				} 
				else 
					DefineLabels("%e%i");
			} 
			else if( cmd_pr.CurToken==FalseTok )
			{   
				DeleteLabels();
			} 
			else if( cmd_pr.CurToken!=StringTok )
			{   
				DefineLabels(cmd_pr.GetStartPosSubstr() );
				cmd_pr.CurToken = 0;
			} 
			else 
				DefineLabels(cmd_pr.TokenIdent.c_str());
			ReDrawFlag |= RFRefresh;
			break;

        case(BackgroundTok):
			cmd_pr.FetchToken();
			if( !cmd_pr.CurToken )
			{   
				PrintLog("No colour specified\n");
			} 
			else if( cmd_pr.CurToken == TransparentTok )
			{   
				UseTransparent = True;
			} 
			else if( cmd_pr.CurToken == NormalTok )
			{   
				UseTransparent = False;
			} 
			else if( cmd_pr.ParseColour(RVal,GVal,BVal) )
			{   
				ReDrawFlag |= RFColour;
				BackColor.SetColor(RVal,GVal,BVal);
#ifndef _WIN32
				pCanv->FBClear = False;
#endif
			} 
			else if( cmd_pr.CurToken )
				PrintLog("Unknown or incorrect colour\n");
			break;

        case(ClipboardTok):
                          if( ClipboardImage() )
                          {   
                              PrintLog("\nUnable to copy to clipboard!\n");
                          }
                          break;

        default:          cmd_pr.CommandError("Unrecognised command");
                          break;
    }

	if(ReDrawFlag) UpdateThisView();

    if( cmd_pr.CurToken )
	{
        if( cmd_pr.FetchToken() ) cmd_pr.CommandError("Warning: Ignoring rest of command");
	}
    return( False );
}

void HaMolView::WrapShiftVal(int iaxis, double value )
{
    double temp;

	double* ptr_val = NULL; 

	if(iaxis == 0 ) ptr_val = &CurRX;
	if(iaxis == 1 ) ptr_val = &CurRY;
	if(iaxis == 2 ) ptr_val = &CurRZ;

    temp = *ptr_val + value;
    while( temp < -1.0 )  temp += 2.0;
    while( temp > 1.0 )   temp -= 2.0;
    *ptr_val = temp;
}


int HaMolView::ExecuteSetCommand(CmdParser& cmd_pr)
{
    int option;
	int RVal, GVal, BVal;

	double dtmp;

    switch( cmd_pr.FetchToken() )
    {   
	case(SlabTok):
	case(SlabModeTok):
		option = -1;
		cmd_pr.FetchToken();
		if( cmd_pr.CurToken==RejectTok )
		{   
			option = SlabReject;
		} 
		else if( cmd_pr.CurToken==HalfTok )
		{   
			option = SlabHalf;
		} 
		else if( cmd_pr.CurToken==HollowTok )
		{   
			option = SlabHollow;
		} 
		else if( cmd_pr.CurToken==SolidTok )
		{   
			option = SlabClose;
		} 
		else if( cmd_pr.CurToken==SectionTok )
			option = SlabSection;
		
		if( option != -1 )
		{   
			SetSlabMode(option);
			if( UseSlabPlane() && (SlabMode() != option) )
				UpdateThisView(RFRefresh);
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;
		
	case(SpecularTok):
		cmd_pr.FetchToken();
		if( cmd_pr.CurToken==TrueTok )
		{   
			FakeSpecular = True;
			UpdateThisView(RFColour);
		} 
		else if( cmd_pr.CurToken==FalseTok )
		{   
			FakeSpecular = False;
			UpdateThisView(RFColour);
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;

	case(SpecPowerTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			SpecPower = 8;
			UpdateThisView(RFColour);
		} 
		else if( cmd_pr.CurToken==NumberTok )
		{   
			if( cmd_pr.TokenValue<=100 )
			{   
				SpecPower = (int)cmd_pr.TokenValue;
				UpdateThisView(RFColour);
			} 
			else 
				PrintLog("Parameter value too large");
		} 
		else 
			PrintLog("Integer value expected\n");
		break;
		
	case(AmbientTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			Ambient = DefaultAmbient;
			UpdateThisView(RFColour);
		} 
		else if( cmd_pr.CurToken==NumberTok )
		{   
			if( cmd_pr.TokenValue<=100 )
			{   
				Ambient = cmd_pr.TokenValue/100.0;
				UpdateThisView(RFColour);
			} 
			else
				PrintLog("Parameter value too large"); 
		} 
		else 
			PrintLog("Integer value expected\n");
		break;

	case(HeteroTok):
		cmd_pr.FetchToken();
		if( cmd_pr.CurToken==TrueTok )
		{   
			HetaGroups = True;
			UpdateThisView(RFRefresh);
		} 
		else if( cmd_pr.CurToken==FalseTok )
		{   
			HetaGroups = False;
			UpdateThisView(RFRefresh);
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;
                                  
	case(HydrogenTok):
		cmd_pr.FetchToken();
		if( cmd_pr.CurToken==TrueTok )
		{   
			Hydrogens = True;
			UpdateThisView(RFRefresh);
		} 
		else if( cmd_pr.CurToken==FalseTok )
		{   
			Hydrogens = False;
			UpdateThisView(RFRefresh);
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;
                                  

	case(BackgroundTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			PrintLog("No colour specified\n");
		} 
		else if( cmd_pr.CurToken == TransparentTok )
		{   
			UseTransparent = True;
			UpdateThisView(RFRefresh);
		} 
		else if( cmd_pr.CurToken == NormalTok )
		{   
			UseTransparent = False;
			UpdateThisView(RFRefresh);
		} 
		else if( cmd_pr.ParseColour(RVal,GVal,BVal) )
		{   
			UpdateThisView(RFColour);
#ifndef IBMPC
			pCanv->FBClear = False;
#endif
		} 
		else if( cmd_pr.CurToken )
			PrintLog("Unknown or incorrect colour\n");;
		break;

	case(BondModeTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken || (cmd_pr.CurToken==AndTok) )
		{   
			ZoneBoth = True;
			UpdateThisView(RFRefresh);
		} 
		else if( cmd_pr.CurToken==OrTok )
		{   
			ZoneBoth = False;
			UpdateThisView(RFRefresh);
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;
            
	case(HBondTok):
		cmd_pr.FetchToken();
		if( (cmd_pr.CurToken==BackboneTok) || (cmd_pr.CurToken==MainChainTok) )
		{   
			HBondMode = True;
			UpdateThisView(RFRefresh);
		} 
		else if( !cmd_pr.CurToken || (cmd_pr.CurToken==SidechainTok) )
		{   
			HBondMode = False;
			UpdateThisView(RFRefresh);
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;
		
	case(SSBondTok):
		cmd_pr.FetchToken();
		if( (cmd_pr.CurToken==BackboneTok) || (cmd_pr.CurToken==MainChainTok) )
		{   
			SSBondMode = True;
			UpdateThisView(RFRefresh);
		} 
		else if( !cmd_pr.CurToken || (cmd_pr.CurToken==SidechainTok) )
		{   
			SSBondMode = False;
			UpdateThisView(RFRefresh);
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;
		
	case(HourGlassTok):
		cmd_pr.FetchToken();
		if( cmd_pr.CurToken==TrueTok )
		{   
			UseHourGlass = True;
		} 
		else if( cmd_pr.CurToken==FalseTok )
		{   
			UseHourGlass = False;
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;

	case(StrandsTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			pCanv->m_SplineCount = 5;
			UpdateThisView(RFRefresh);
		} 
		else if( cmd_pr.CurToken==NumberTok )
		{   
			if( (cmd_pr.TokenValue>0) && (cmd_pr.TokenValue<=5) )
			{   
				pCanv->m_SplineCount = (int)cmd_pr.TokenValue;
				UpdateThisView(RFRefresh);
			} 
			else if( cmd_pr.TokenValue==9 )
			{   
				pCanv->m_SplineCount = 9;
				UpdateThisView(RFRefresh);
			} 
			else 
				PrintLog("Invalid parameter setting\n");
		} 
		else 
			PrintLog("Integer value expected\n");
		break;
		
	case(MouseTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken || (cmd_pr.CurToken==RasMolTok) )
		{   
			if( pApp->gui_mode )
				SetMouseMode( MMRasMol );
		} 
		else if( cmd_pr.CurToken==InsightTok )
		{   
			if( pApp->gui_mode )
				SetMouseMode( MMInsight );
		} 
		else if( cmd_pr.CurToken==QuantaTok )
		{   
			if( pApp->gui_mode )
				SetMouseMode( MMQuanta );
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;
				
	case(AxesTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken || (cmd_pr.CurToken==FalseTok) )
		{   
			DrawAxes = False;
			UpdateThisView(RFRefresh);
		} 
		else if( cmd_pr.CurToken == TrueTok )
		{   
			DrawAxes = True;
			UpdateThisView(RFRefresh);
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;
		
	case(BoundBoxTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken || (cmd_pr.CurToken==FalseTok) )
		{   
			DrawBoundBox = False;
			UpdateThisView(RFRefresh);
		} 
		else if( cmd_pr.CurToken == TrueTok )
		{   
			DrawBoundBox = True;
			UpdateThisView(RFRefresh);
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;
		
	case(UnitCellTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken || (cmd_pr.CurToken==FalseTok) )
		{   
			DrawUnitCell = False;
			UpdateThisView(RFRefresh);
		} 
		else if( cmd_pr.CurToken == TrueTok )
		{   
			DrawUnitCell = True;
			UpdateThisView(RFRefresh);
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;
		
	case(VectPSTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken || (cmd_pr.CurToken==FalseTok) )
		{   
			UseOutLine = False;
		} 
		else if( cmd_pr.CurToken == TrueTok )
		{   
			UseOutLine = True;
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;
					
	case(RadiusTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			ProbeRadius = SolventDots? 1.2 : 0.0;
		} 
		else if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken == FloatTok)
		{   
			dtmp = (cmd_pr.CurToken==NumberTok) ? cmd_pr.TokenValue : cmd_pr.TokenValueFloat;

			if( dtmp > 6.0 )
			{   
				PrintLog("Parameter value too large");
			} 
			ProbeRadius = dtmp;
		} 
		else 
			PrintLog("Integer value expected\n");
		break;
		
	case(SolventTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken || (cmd_pr.CurToken==FalseTok) )
		{   
			SolventDots = False;
			ProbeRadius = 0.0;
		} 
		else if( cmd_pr.CurToken == TrueTok )
		{   
			SolventDots = True;
			ProbeRadius = 1.2;
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;
		
	case(FontSizeTok):
		cmd_pr.FetchToken();
		if( cmd_pr.CurToken==NumberTok )
		{   
			if( cmd_pr.TokenValue<=32 )
			{   
				pCanv->SetFontSize((int)cmd_pr.TokenValue);
				if( DrawLabels || (!MonitList.empty() && DrawMonitDistance) )
					UpdateThisView(RFRefresh);
			} 
			else 
				PrintLog("Parameter value too large");
		} 
		else if( !cmd_pr.CurToken )
		{   
			pCanv->SetFontSize(8);
			if( DrawLabels )
				UpdateThisView(RFRefresh);
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;
				
	case(StereoTok):  
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			SetStereoMode(False);
			StereoAngle = 6.0;
			UpdateThisView(RFRefresh);
		} 
		else if( cmd_pr.CurToken==TrueTok )
		{   
			SetStereoMode(True);
			UpdateThisView(RFRefresh);
		} 
		else if( cmd_pr.CurToken==FalseTok )
		{   
			SetStereoMode(False);
			UpdateThisView(RFRefresh);
		} 
		else if( cmd_pr.CurToken == '-' )
		{   
			if( !cmd_pr.NextIf(NumberTok,"Integer value expected") )
			{   
				StereoAngle = -cmd_pr.TokenValue;
				SetStereoMode(True);
				UpdateThisView(RFRefresh);
			}
		} 
		else if( cmd_pr.CurToken == '+' )
		{   
			if( !cmd_pr.NextIf(NumberTok,"Integer value expected") )
			{   
				StereoAngle = cmd_pr.TokenValue;
				SetStereoMode(True);
				UpdateThisView(RFRefresh);
			}
		} 
		else if( cmd_pr.CurToken==NumberTok )
		{   
			StereoAngle = cmd_pr.TokenValue;
			SetStereoMode(True);
			UpdateThisView(RFRefresh);
		} 
		else 
			PrintLog("Invalid command syntax\n");
		break;
		
	case(PickingTok):
		switch( cmd_pr.FetchToken() )
		{   
		case(TrueTok):     
		case(0):
		case(IdentifyTok): 
			SetPickMode(PickIdent); 
			break;
		case(FalseTok):
		case(NoneTok):     
			SetPickMode(PickNone);  
			break;
		case(LabelTok):    
			SetPickMode(PickLabel); 
			break;
		case(DistanceTok): 
			SetPickMode(PickDist);  
			break;
		case(AngleTok):    
			SetPickMode(PickAngle); 
			break;
		case(TorsionTok):  
			SetPickMode(PickTorsn); 
			break;
		case(MonitorTok):  
			SetPickMode(PickMonit);
			break;
		case(CentreTok):   
			SetPickMode(PickCentr); 
			break;
		default:           
			PrintLog("Invalid parameter setting\n");
		}
		break;
			
	case(BondTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken || (cmd_pr.CurToken==FalseTok) )
		{   
			DrawDoubleBonds = False;
			UpdateThisView(RFRefresh);
		} 
		else if( cmd_pr.CurToken == TrueTok )
		{   
			DrawDoubleBonds = True;
			UpdateThisView(RFRefresh);
		} 
		else 
			PrintLog("Invalid parameter setting\n");
		break;
			
        case(MonitorTok):
           cmd_pr.FetchToken();
            if( !cmd_pr.CurToken || (cmd_pr.CurToken==TrueTok) )
            {   
                DrawMonitDistance = True;
				UpdateThisView(RFRefresh);
            } 
			else if( cmd_pr.CurToken == FalseTok )
            {   
                DrawMonitDistance = False;
				UpdateThisView(RFRefresh);
            } 
			else 
				PrintLog("Invalid parameter setting\n");
            break;

        case(CartoonTok):
           cmd_pr.FetchToken();
            if( !cmd_pr.CurToken )
            {   
                DrawBetaArrows = True;
                CartoonHeight = 0.48;
				UpdateThisView(RFRefresh);
            } 
			else if( cmd_pr.CurToken==TrueTok )
            {   
                DrawBetaArrows = True;
				UpdateThisView(RFRefresh);
            } 
			else if( cmd_pr.CurToken==FalseTok )
            {   
                DrawBetaArrows = False;
				UpdateThisView(RFRefresh);
            } 
			else if( cmd_pr.CurToken==NumberTok || cmd_pr.CurToken == FloatTok)
            {   				
				CartoonHeight = cmd_pr.TokenValueFloat;
				UpdateThisView(RFRefresh);

			}
			else 
				PrintLog("Invalid parameter setting\n");
            break;
			
        case(BackFadeTok):
            cmd_pr.FetchToken();
            if( !cmd_pr.CurToken || (cmd_pr.CurToken==FalseTok) )
            {   
                UseBackFade = False;
				UpdateThisView(RFColour);
            } 
			else if( cmd_pr.CurToken == TrueTok )
            {   
                UseBackFade = True;
				UpdateThisView(RFColour);
            } 
			else 
				PrintLog("Invalid parameter setting\n");
            break;

        case(TransparentTok):
            cmd_pr.FetchToken();
            if( !cmd_pr.CurToken || (cmd_pr.CurToken==FalseTok) )
            {   
				UseTransparent = False;
				UpdateThisView(RFColour);
            } 
			else if( cmd_pr.CurToken == TrueTok )
            {   
				UseTransparent = True;
				UpdateThisView(RFColour);
            } 
			else 
				PrintLog("Invalid parameter setting\n");
            break;

        case(DepthCueTok):
            cmd_pr.FetchToken();
            if( !cmd_pr.CurToken || (cmd_pr.CurToken==FalseTok) )
            {   
                UseDepthCue = False;
				UpdateThisView(RFColour);
            } 
			else if( cmd_pr.CurToken == TrueTok )
            {   
                UseDepthCue = True;
				UpdateThisView(RFColour);
            } 
			else 
				PrintLog("Invalid parameter setting\n");
            break;

        case(SequenceTok):
            cmd_pr.FetchToken();
            if( !cmd_pr.CurToken || (cmd_pr.CurToken==FalseTok) )
            {   
				HaMolecule::SeqFormat = False;
            } 
			else if( cmd_pr.CurToken == TrueTok )
            {   
				HaMolecule::SeqFormat = True;
            } 
			else 
				PrintLog("Invalid parameter setting\n");
            break;

        default:
            PrintLog("Invalid parameter name\n");
    }
	return True;
}

int HaMolView::ExecuteColourCommand(CmdParser& cmd_pr)
{
    int flag;
	int RVal,GVal,BVal;
	
    flag = 0;
    switch( cmd_pr.FetchToken() )
    {   
	case(AtomTok):
		cmd_pr.FetchToken();
	default:
		if( !cmd_pr.CurToken )
		{   
			PrintLog("No colour specified\n");
		} 
		else
		{
			switch( cmd_pr.CurToken )
            {   
			case(CPKTok):         
				CPKColourAttrib(); 
				break;
				
			case(AminoTok):       
				AminoColourAttrib();
				break;
				
			case(ShapelyTok):     
				ShapelyColourAttrib();
				break;
                				
			case(GroupTok):       
				GroupsColourAttrib();
				break;

			case(ResidueTok):       
				ScaleColourAttrib(ResidueAttr);
				break;
			
			case(ChainTok):       
				ScaleColourAttrib(ChainAttr);
				break;
				
			case(ChargeTok):      
				ScaleColourAttrib(ChargeAttr);
				break;
				
			case(TemperatureTok): 
				ScaleColourAttrib(TempAttr);
				break;
				
			case(StructureTok):   
				StructColourAttrib();
				break;
				
			default:  
				if( cmd_pr.ParseColour(RVal,GVal,BVal) )
				{   
					MonoColourAttrib(RVal,GVal,BVal);
				} 
				else 
					PrintLog("Unknown or incorrect colour\n");;
            }

		}
		break;
		
	case(BondTok):    
	case(DashTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			PrintLog("No colour specified\n");
		} 
		else if( cmd_pr.CurToken==NoneTok )
		{   
			ColourBondNone();
		} 
		else if( cmd_pr.ParseColour(RVal,GVal,BVal) )
		{   
			ColourBondAttrib(RVal,GVal,BVal);
		} 
		else PrintLog("Unknown or incorrect colour\n");;
		break;
		
	case(BackboneTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			PrintLog("No colour specified\n");
		} 
		else if( cmd_pr.CurToken==NoneTok )
		{   
			ColourBackNone();
		} 
		else if( cmd_pr.ParseColour(RVal,GVal,BVal) )
		{   
			ColourBackAttrib(RVal,GVal,BVal);
		} 
		else PrintLog("Unknown or incorrect colour\n");;
		break;
		
	case(SSBondTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			PrintLog("No colour specified\n");
		} 
		else if( cmd_pr.CurToken==NoneTok )
		{   
			ColourSSBondNone();
		} 
		else if( cmd_pr.ParseColour(RVal,GVal,BVal) )
		{   
			ColourSSBondAttrib(RVal,GVal,BVal);
		} 
		else PrintLog("Unknown or incorrect colour\n");
		break;
		
	case(HBondTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			PrintLog("No colour specified\n");
		} 
		else if( cmd_pr.CurToken==NoneTok )
		{   
			ColourHBondNone();
		} 
		else if( cmd_pr.CurToken==TypeTok )
		{   
			ColourHBondType();
		} 
		else if( cmd_pr.ParseColour(RVal,GVal,BVal) )
		{   
			ColourHBondAttrib(RVal,GVal,BVal);
		} 
		else PrintLog("Unknown or incorrect colour\n");
		break;
		
	case(DotsTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			PrintLog("No colour specified\n");
		} 
		else if( cmd_pr.CurToken== PotentialTok )
		{   
			ColourDotsPotential();
		} 
		else if( cmd_pr.ParseColour(RVal,GVal,BVal) )
		{   
			ColourDotsAttrib(RVal,GVal,BVal);
		} 
		else 
			PrintLog("Unknown or incorrect colour\n");
		break;
		
	case(MonitorTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			PrintLog("No colour specified\n");
		} 
		else if( cmd_pr.CurToken == NoneTok )
		{   
			ColourMonitNone();
		} 
		else if( cmd_pr.ParseColour(RVal,GVal,BVal) )
		{   
			ColourMonitAttrib(RVal,GVal,BVal);
		} 
		else 
			PrintLog("Unknown or incorrect colour\n");
		break;
		
	case(AxesTok):
	case(BoundBoxTok):
	case(UnitCellTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			PrintLog("No colour specified\n");
		} 
		else if( cmd_pr.ParseColour(RVal,GVal,BVal) )
		{   
			BoxColor.SetColor(RVal,GVal,BVal);
		} 
		else 
			PrintLog("Unknown or incorrect colour\n");
		break;
		
	case(LabelTok):
		cmd_pr.FetchToken();
		if( !cmd_pr.CurToken )
		{   
			PrintLog("No colour specified\n");
		} 
		else if( cmd_pr.CurToken==NoneTok )
		{   
			UseLabelCol = False;
		} 
		else if( cmd_pr.ParseColour(RVal,GVal,BVal) )
		{   
			LabelColor.SetColor(RVal,GVal,BVal);
			ReDrawFlag |= RFColour;
			UseLabelCol = True;
		} 
		else 
			PrintLog("Unknown or incorrect colour\n");
		break;
		
	case(TraceTok): 
	case(RibbonTok):
	case(CartoonTok):  flag = RibColBoth;     break;
	case(Ribbon1Tok):  flag = RibColInside;   break;
	case(Ribbon2Tok):  flag = RibColOutside;  break;
    }
	
    if( flag )
    {   
		cmd_pr.FetchToken();
        if( !cmd_pr.CurToken )
        {   
			PrintLog("No colour specified\n");
        } 
		else if( cmd_pr.CurToken==NoneTok )
        {   
            ColourRibbonNone(flag);
        } 
		else if( cmd_pr.ParseColour(RVal,GVal,BVal) )
        {   
            ColourRibbonAttrib(flag,RVal,GVal,BVal);
        } 
		else 
			PrintLog("Unknown or incorrect colour\n");
    }
	ReDrawFlag |= RFRefresh;
	return True;
}

void HaMolView::ConnectObject(Object3D* pObj)
{
	MoleculesType::iterator mol_itr;
	if(pObj != NULL)
	{
		ForEachMol_VIEW
		{
			(*mol_itr)->SetConnected(false);
		}
		pObj->SetConnected(true);
		CalcRotCenter();
		m_screen_transform = false;
		
		char buf[256];
		sprintf(buf,"Transform Dials Connected to the Molecule: %s", pObj->GetObjName() );
		PrintMessage(buf); 
	}
}

void HaMolView::SetMouseMode(int mode )
{
    MouseMode = mode;
}


int HaMolView::CreateImage()
{
	pCanv->AllocImage();
	return((int)pCanv->FBuffer);
}

void run_eigen_vec_animation(HaMolView* p_mol_view, HaVec_double* p_vec, AtomContainer* p_at_cont )
{
	p_mol_view->AnimateEigenVectorInternal(*p_vec, p_at_cont );
}


int HaMolView::AnimateEigenVector( HaVec_double& evec, AtomContainer* p_at_cont )
{
	std::thread anim_t(run_eigen_vec_animation, this, &evec, p_at_cont);
	anim_t.detach();

	return TRUE;
}

int HaMolView::AnimateEigenVectorInternal( HaVec_double& evec, AtomContainer* p_at_coll )
{
	char buf[256];
	HaVec_double ref_crd;
	MolSet* pmset = GetMolSet();
	try
	{
		std::auto_ptr<AtomIterator> paitr(p_at_coll->GetAtomIteratorPtr());
		HaAtom* aptr;
		AtomGroup atgrp;
		for(aptr = paitr->GetFirstAtom(); aptr ; aptr = paitr->GetNextAtom())
		{
			if( aptr->GetHostMolSet() != pmset ) throw std::runtime_error(" Animated Atom Collection and Molecular View associated with different Molecular sets");
			atgrp.push_back(aptr);
		}
		int na = atgrp.GetNAtoms();
		if( na == 0 ) throw std::runtime_error(" Animated Atom Collection is empty ");
		if( evec.size() != 3*na ) 
		{
			sprintf(buf,"Size of Eigenvector = %d is not equal 3*natoms = %d ",evec.size(),3*na);
			throw std::runtime_error( buf );
		}
		ref_crd.resize(3*na);
		int i;
		for(i = 0; i < na; i++)
		{
			ref_crd[3*i  ] = atgrp[i]->GetX_Ang();
			ref_crd[3*i+1] = atgrp[i]->GetY_Ang();
			ref_crd[3*i+2] = atgrp[i]->GetZ_Ang();
		}

		typedef std::chrono::system_clock Time;
		auto update_time = Time::now();

		int npt_cycle = 24;
		int ipt = 0;
		double ampl = 2.0;
		double cycle_time = 2.0;

		std::chrono::duration<long, std::milli> update_interval = std::chrono::milliseconds((long)(1000 * cycle_time / npt_cycle));

		to_stop_animation = FALSE;

		for(;;)
		{
			if( to_stop_animation ) break;
			if( update_time < Time::now() )
			{
				for(i = 0; i < na; i++)
				{
					atgrp[i]->SetX_Ang( ref_crd[3*i   ] + evec[3*i  ]*ampl*sin( 2.0*PI*((double)ipt/npt_cycle) ) );
					atgrp[i]->SetY_Ang( ref_crd[3*i+1 ] + evec[3*i+1]*ampl*sin( 2.0*PI*((double)ipt/npt_cycle) ) );
					atgrp[i]->SetZ_Ang( ref_crd[3*i+2 ] + evec[3*i+2]*ampl*sin( 2.0*PI*((double)ipt/npt_cycle) ) );
				}
				UpdateThisView(RFRefresh | RFApply);
				ipt++;

				update_time = Time::now();
				update_time += update_interval;
			}
		}
		for(i = 0; i < na; i++)
		{
			atgrp[i]->SetX_Ang( ref_crd[3*i   ] );
			atgrp[i]->SetY_Ang( ref_crd[3*i+1 ] );
			atgrp[i]->SetZ_Ang( ref_crd[3*i+2 ] );
		}
	}
	catch( std::exception& ex )
	{
		PrintLog("Error in HaMolView::AnimateEigenVector() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return TRUE;
}

void HaMolView::StopAnimation()
{
	to_stop_animation = TRUE;
}

Monitor::Monitor()
{
  src=dst=NULL;
  dist = 0.0;
  col = 0;
  monitor_show_mode = MONITOR_SHOW_DISTANCE;
}

Monitor::~Monitor()
{

}

bool  Monitor::operator ==( const Monitor& ref) const
{
  return ( (src == ref.src) && (dst == ref.dst));
}

bool Monitor::operator <( const Monitor& ref) const
{
   if( (src < ref.src)) 
         return true;
   else if( src > ref.src )
         return false;
   return ( dst < ref.dst);
}

void Monitor::SetShowValue(bool set )
{
	if (set)
	{
		monitor_show_mode = MONITOR_SHOW_DISTANCE;
	}
	else
	{
		monitor_show_mode = MONITOR_SHOW_NONE;
	}
}

bool Monitor::ToShowValue()
{
	if (monitor_show_mode == MONITOR_SHOW_DISTANCE) return true;
	return false;
}


#if defined(HA_NOGUI)
int 
HaMolView::BroadcastCurrAtom()
{
	return TRUE;
}

void 
HaMolView::UpdateThisView( int lHint)
{
   
}

#endif
