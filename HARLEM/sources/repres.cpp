/*! \file repres.cpp

   Some Plotting function of HaMolView class 

   \author Igor Kurnikov
   \date 1997-2002

   derived from:
 
   repres.c
   RasMol2 Molecular Graphics
   Roger Sayle, August 1995
   Version 2.6
*/
#include "haconst.h"

#include <float.h>
#include <stdlib.h>
#include <string>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#include <boost/algorithm/string.hpp>

#define REPRES
#include "command.h"
#include "abstree.h"
#include "hamolview.h"

#include "haatom.h"
#include "hamolecule.h"
#include "etcoupl.h"
#include "electrostmod.h"
#include "object3d.h"
#include "canvas3d.h"

typedef struct { double dx, dy, dz; } DotVector;
typedef struct {
        DotVector  *probe;
        DotVector  *dots;
        int count;
    } ElemDotStruct;
 

static ElemDotStruct  *ElemDots;
static HaAtom  *Exclude;


/*============================*/
/*  Label Handling Functions  */
/*============================*/
 

int 
HaMolView::DeleteLabels()
{
    HaAtom  *aptr;

	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
        if( aptr->Selected() )
        {   
			if( !aptr->label.empty() )
            {   
                aptr->label.resize(0);
            }
        }
	}
    return( True );
}
 
 
void HaMolView::DefineLabels( const char* label )
{
    HaAtom  *aptr;

	std::string lbl_str = label;
	boost::trim(lbl_str);

    if( lbl_str.size() == 0 ) return;

    DrawLabels = True;

	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
        if( aptr->Selected() )
        {   
			aptr->label = lbl_str;
        }
	}
}
 
void 
HaMolView::DefaultLabels( int enable )
{
    HaChain  *chain;
    HaAtom  *aptr;

	MoleculesType::iterator mol_itr;

	ForEachMol_VIEW
	{
		if( (*mol_itr)->GetNRes() > 1 )
		{   
			AtomIteratorMolecule aitr(*mol_itr);
			for(aptr= aitr.GetFirstAtom(); aptr;aptr= aitr.GetNextAtom())
			{
				if( aptr->IsAlphaCarbon() || aptr->IsSugarPhosphate() )
				{   
					if( aptr->Selected() )
					{   
						chain=aptr->GetHostChain();
						if( enable )
						{   
							if( (*mol_itr)->GetNChains() > 1 )
							{   
								if( isdigit(chain->ident) )
								{   
									aptr->label = "%n%r:%c";
								} 
								else
								{   
									aptr->label = "%n%r%c";
								}
							} 
							else
							{   
								aptr->label = "%n%r";
							}
						} 
						else if( !aptr->label.empty() )
						{   
							aptr->label.resize(0);
						}
						ReDrawFlag |= RFRefresh;
					}
					break;
				}
			}
			
		} 
		else /* Small Molecule! */
		{
			AtomIteratorMolecule aitr(*mol_itr);
			for(aptr= aitr.GetFirstAtom(); aptr;aptr= aitr.GetNextAtom())
			{
				if( (aptr->Selected() ) && (aptr->GetElemNo() != 6)
					&& (aptr->GetElemNo() != 1) )
				{   
					if( enable )
					{   
						aptr->label = "%e";
					} 
					else if( !aptr->label.empty() )
					{   
						aptr->label.resize(0);
					}
					ReDrawFlag |= RFRefresh;
				}
			}
		}
	}
	
    DrawLabels = True;
}


void 
HaMolView::DisplayLabels()
{
    register HaChain  *chain;
    register HaResidue  *group;
    register HaAtom  *aptr;
    register int col,z;
 
    char buffer[256];
  
    if( !UseSlabPlane() )
    {   
		z = GetImageRadius() + ZOffset();
    } 
	else 
		z = SlabValue() - 1;

    AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
        if( !aptr->label.empty() )
        {   
			chain=aptr->GetHostChain();
			group=aptr->GetHostRes();
			/* Peform Label Slabbing! */
            if( !ZValid_v(aptr->z) )
                continue;
 
            FormatLabel(chain,group,aptr,aptr->label.c_str(),buffer);
 
            if( !UseLabelCol )
            {   /* Depth-cue atom labels */
                /* col = aptr->col + (ColorDepth*                  */
                /*       (aptr->z+GetImageRadius()-ZOffset))/GetImageSize(); */
                col = aptr->col + (ColourMask>>1);
            } 
			else 
				col = LabelColor.cidx;
 
            /* (aptr->z+2) + ((aptr->IsDrawSphere())?aptr->irad:0); */
            pCanv->DisplayTextString(aptr->x+4,aptr->y,z,buffer,col);
        }
	}
}

static char *FormatInteger(char* ptr,int value )
{
    char buffer[10];
    register char *tmp;

    if( value<0 )
    {   
		value = -value;
        *ptr++ = '-';
    }

    if( value>9 )
    {   
		tmp = buffer;
        while( value>9 )
        {   
			*tmp++ = (char)(value%10) + '0';
            value /= 10;
        }

        *ptr++ = (char)value + '0';
        do { 
			tmp--; 
            *ptr++ = *tmp;
        } while( tmp != buffer );
    } 
	else 
		*ptr++ = (char)value + '0';
    return( ptr );
}

#define AminoCodeMax  28
static char *AminoCode = "AGLSVTKDINEPRFQYHCMWBZPPACGTX";

void 
HaMolView::FormatLabel(HaChain* chain, HaResidue* group, HaAtom* aptr, const char* label, 
				      char* ptr )
{
    char ch;
    int i,j,len;
	const char* name;

    while( *label )
    {  
	   ch = *label++;
       if( ch=='%' )
       {   
		   ch = *label++;
           if( isupper(ch) )
             ch = tolower(ch);

           switch( ch )
           {   
		       case('a'):  /* Atom Name */
                           name = aptr->GetName();
                           for( j=0; j<4; j++ )
                               if( name[j]!=' ' )
                                   *ptr++ = name[j];
                           break;

               case('b'):  /* Temperature/B-factor */
               case('t'):  j= sprintf(ptr,"%8.3f",aptr->tempf);
				           ptr += j;
                           break;

               case('c'):  /* Chain Identifier */
               case('s'):  *ptr++ = chain->ident;
                           break;
               case('f'):  /* Force Field Atom Symbol name */
						   name = aptr->GetFFSymbol();
				           for( j=0; j<4; j++ )
                               if( name[j]!=' ' )
                                   *ptr++ = name[j];
                           break;
               case('e'):  /* Element Type */
                           i = aptr->GetElemNo();
                           *ptr++ = ElementArr[i].symbol[0];
                           if( ElementArr[i].symbol[1]!=' ' )
                               *ptr++ = ElementArr[i].symbol[1];
                           break;

               case('i'):  /* Atom Number */
                           ptr = FormatInteger(ptr,(int)aptr->GetSerNo());
                           break;

               case('m'):  /* Amino (Nucliec) Acid Code */
                           if( group->refno <= AminoCodeMax )
                           {   
							   *ptr++ = AminoCode[group->refno];
                           } 
						   else 
							   *ptr++ = '?';
                           break;

               case('n'):  /* Residue Name   */
                           name = group->GetName();
						   len = strlen(name);
                           for( j=0; j<len; j++ )
                               if( name[j]!=' ' )
                                   *ptr++ = name[j];
                           break;

               case('r'):  /* Residue Number */
                           ptr = FormatInteger(ptr,group->serno);
                           break;

               case('%'):  *ptr++ = '%';
                           break;
           }
       } else if( (ch>=' ') && (ch<='~') )
           *ptr++ = ch;
    }
    *ptr = '\0';
}



/*==============================*/
/*  Monitor Handling Functions  */
/*==============================*/


 
void HaMolView::DeleteAllMonitors()
{
	MonitList.clear();
}

void HaMolView::DeleteAtomPairMonitor(HaAtom* src, HaAtom* dst)
{
	list<Monitor>::iterator mitr;
	for (mitr = MonitList.begin(); mitr != MonitList.end(); )
	{
		if (((*mitr).src == src) && ((*mitr).dst == dst) ||
			((*mitr).src == dst) && ((*mitr).dst == src))
		{
			mitr = MonitList.erase(mitr);
			continue;
		}
		mitr++;
	}
}
 

void HaMolView::AddAtomPairMonitor( HaAtom* src, HaAtom* dst )
{
    double dx,dy,dz;
    double dist;
	char buf1[256],buf2[256];
 
    /* Delete an already existing monitor! */
	
	src->FillRef(buf1);
	dst->FillRef(buf2);

	PrintLog(" Adding Monitor between atoms %s and %s \n", buf1,buf2);

	list<Monitor>::iterator mitr;

	for(mitr = MonitList.begin(); mitr != MonitList.end(); mitr++)
    {
         if( ((*mitr).src == src) && ((*mitr).dst==dst)  ||
             ((*mitr).src == dst) && ((*mitr).dst==src) )
         {   
			 MonitList.erase(mitr);
             return;
         }
	}
 

	Monitor mtr;
 
    dx = src->GetX() - dst->GetX();
    dy = src->GetY() - dst->GetY();
    dz = src->GetZ() - dst->GetZ();
 
    /* ptr->dist = 100.0*CalcDistance(src,dst) */
    dist = sqrt( dx*dx + dy*dy + dz*dz );
    mtr.dist = dist;
 
    mtr.src = src;
    mtr.dst = dst;
    mtr.col = 0;
	MonitList.push_back(mtr);
}


void  HaMolView::CreateMonitor( int src, int dst )
{
    register HaAtom  *aptr;
    register HaAtom  *sptr;
    register HaAtom  *dptr;
    register int done;
    char buffer[20];
 
    if( src == dst )
    {   
        PrintLog("\n Error: Duplicate atom serial numbers!\n");
        return;
    }
 
    done = False;
    sptr = (HaAtom *)NULL;
    dptr = (HaAtom *)NULL;

    AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{   
		if( aptr->GetSerNo() == src )
		{   
			sptr = aptr;
			if( dptr )
			{   
				done = True;
				break;
			}
		} 
		else if( aptr->GetSerNo() == dst )
		{   
			dptr = aptr;
			if( sptr )
			{   
				done = True;
				break;
			}
		}
	}
 
    if( !done )
    {   
        PrintLog("\n Error: Atom serial number");
        if( sptr )
        {   
			sprintf(buffer," %d",dst);
        } 
		else if( dptr )
        {   
			sprintf(buffer," %d",src);
        } 
		else sprintf(buffer,"s %d and %d",src,dst);
        PrintLog("%s not found!\n",buffer); 
    } 
	else
	{
		AddAtomPairMonitor(sptr, dptr);
	}
}
 
 
void HaMolView::DisplayMonitors()
{
    HaAtom  *s;
    HaAtom  *d;
    Monitor *ptr;
    int x,y,z;
    int sc,dc;
    int col;
 
    double dist;
    char buf[256];
 
    if( !UseSlabPlane() )
    {   
		z = GetImageRadius() + ZOffset();
    } 
	else 
		z = SlabValue() - 1;


	list<Monitor>::iterator mitr;
    for( mitr=MonitList.begin(); mitr != MonitList.end(); mitr++ )
    {   
		ptr = &(*mitr);
		s = ptr->src;
        d = ptr->dst;
 
        if( !ptr->col )
        {   
			sc = s->col;
            dc = d->col;
        } 
		else 
			sc = dc = ptr->col;
 
        pCanv->ClipDashVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
 
		if (DrawMonitDistance)
		{
			if (ZValid_v((s->z + d->z) / 2))
			{
				x = (s->x + d->x) / 2;
				y = (s->y + d->y) / 2;

				if (!UseLabelCol)
				{   /* Use Source atom colour! */
					col = sc + (ColourMask >> 1);
				}
				else
					col = LabelColor.cidx;


				dist = Vec3D::CalcDistance(s, d, ANGSTROM_U);
				sprintf(buf, "%7.3f", dist);

				pCanv->DisplayTextString(x + 4, y, z, buf, col);
			}
		}
    }
}
 
 
/*=========================*/
/*  Dot Surface Functions  */
/*=========================*/

  
void 
HaMolView::DeleteDotSurfaces()
{

	MolSet* pmset = GetMolSet();

 	list<Object3D*>::iterator oitr;
	for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end();)
	{
		if( (*oitr)->GetObjType() == OBJ3D_DOT_SURFACE )
		{
			delete (*oitr);
			oitr= pmset->ViewObjects.erase(oitr);
			continue;
		}
		oitr++;
	}

    DrawDots = False;
}
 
static double max_vdw_rad = 2.5; 

int HaMolView::TestDot( double x, double y, double z , bool solvent_access)
{
    HaAtom  *aptr;
    int lx,ly,lz;
    int ux,uy,uz;
    double dx,dy,dz;
    int ix,iy,iz;
    double dist;
    double rad;
    int i;

	double max_span = max_vdw_rad;
	if( solvent_access ) max_span += ProbeRadius;
 
    lx = (int)((x - max_span - HashTable.xmin)/HashTable.dx);
    if( lx >= HashTable.nx ) return( True );
    ly = (int)((y - max_span - HashTable.ymin)/HashTable.dy);
    if( ly >= HashTable.ny ) return( True );
    lz = (int)((z - max_span - HashTable.zmin)/HashTable.dz);
    if( lz >= HashTable.nz ) return( True );
 
    ux = (int)((x + max_span - HashTable.xmin)/HashTable.dx);
    if( ux < 0 ) return( True );
    uy = (int)((y + max_span - HashTable.ymin)/HashTable.dy);
    if( uy < 0 ) return( True );
    uz = (int)((z + max_span - HashTable.zmin)/HashTable.dz);
    if( uz < 0 ) return( True );
 
    if( lx < 0 ) lx = 0;  if( ux >= HashTable.nx ) ux = HashTable.nx-1;
    if( ly < 0 ) ly = 0;  if( uy >= HashTable.ny ) uy = HashTable.ny-1;
    if( lz < 0 ) lz = 0;  if( uz >= HashTable.nz ) uz = HashTable.nz-1;
 
    for( ix=lx; ix<=ux; ix++ )
       for( iy=ly; iy<=uy; iy++ )
          for( iz=lz; iz<=uz; iz++ )
          {   
			  i = HashTable.GetLinCellIdx(ix,iy,iz);
			  VecPtr::iterator aitr;
              for( aitr = HashTable[i].begin(); aitr != HashTable[i].end(); aitr++)
                  if( (HaAtom*) (*aitr) != Exclude )
                  {   
					  aptr = (HaAtom*)(*aitr);
					  rad = HaAtom::ElemVDWRadius(aptr->GetElemNo());
					  if(solvent_access)
						  rad += ProbeRadius;

                      rad = rad * rad;
 
                      /* Optimized Test! */
                      dx = aptr->GetX() - x;
                      if( (dist=dx*dx) < rad )
                      {   
						  dy = aptr->GetY() - y;
                          if( (dist+=dy*dy) < rad )
                          {   
							  dz = aptr->GetZ() - z;
                              if( (dist+=dz*dz) < rad )
                                  return( False );
                          }
                      }
                  }
          }
    return( True );
}
 
 
void 
HaMolView::InitElemDots()
{
    register int i,size;
 
    size = MAXELEMNO*sizeof(ElemDotStruct);
    ElemDots = (ElemDotStruct *)malloc(size);
    if( !ElemDots ) 
    {
	    ErrorInMod("HaMolView::InitElemDots",
                       "failed to allocate ElemDotStruct");
	    return;
    }
 
    for( i=0; i<MAXELEMNO; i++ )
        ElemDots[i].count = 0;
}
 
 
void HaMolView::AddElemDots( int elem, int density )
{
    DotVector  *ptr;
    DotVector  *probe;
    double x, y, z, p, q, xy;
    int equat,vert,horz;
    int count;
    int i,j,k;
	double temp, rad;
 
    if( SolventDots || fabs(ProbeRadius) < DBL_EPSILON )
    {   
		rad = HaAtom::ElemVDWRadius(elem);
    } 
	else 
		rad = HaAtom::ElemVDWRadius(elem) + ProbeRadius;
 
    count = (int)(density*rad*rad);
    ptr = (DotVector *)malloc(count*sizeof(DotVector));
    if( !ptr ) 
    {
	 ErrorInMod("HaMolView::AddElemDots",
                    "failed to allocate DotVector");
         return;
    }
 
    if( SolventDots )
    {   
		probe = (DotVector *)malloc(count*sizeof(DotVector));
        if( !probe ) 
	{
	    ErrorInMod("HaMolView::AddElemDots",
                       "failed to allocate probe vectors");
            return;
	}
        temp = rad + ProbeRadius;
    } 
	else 
		probe = NULL;
 
    equat = (int)sqrt(PI*count);
    if( !(vert=equat>>1) )
        vert = 1;
 
    i = 0;
    for( j=0; (i<count) && (j<vert); j++ )
    {   
		p = (PI*j)/(double)vert;
        z = cos(p);  xy = sin(p);
        horz = (int)(equat*xy);
        if( !horz ) horz = 1;
 
        for( k=0; (i<count) && (k<horz); k++ )
        {   
			q = (2.0*PI*k)/(double)horz;
            x = xy*sin(q);
            y = xy*cos(q);
 
            ptr[i].dx = rad*x;
            ptr[i].dy = rad*y;
            ptr[i].dz = rad*z;
            if( probe )
            {   
				probe[i].dx = temp*x;
                probe[i].dy = temp*y;
                probe[i].dz = temp*z;
            }
            i++;
        }
    }
    ElemDots[elem].probe = probe;
    ElemDots[elem].dots = ptr;
    ElemDots[elem].count = i;
}
 
 
void 
HaMolView::FreeElemDots()
{
    register int i;
 
    for( i=0; i<MAXELEMNO; i++ )
        if( ElemDots[i].count )
            free( ElemDots[i].dots );
    free( ElemDots );
}
 
 
void HaMolView::CalculateDotSurface( int density )
{
    DotVector  *probe;
    DotVector  *ptr;
    HaAtom  *aptr;
    int i,count;
    int elem;
 
	MolSet* pmset = GetMolSet();

	if(pmset->GetNMol() == 0 )
		return;

	DotStruct* dot_str = new DotStruct;

	if(dot_str != NULL)
		pmset->AddObject3D(dot_str);
	else
		return;
 
    InitElemDots();
    BuildHashTable();

	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		aptr->solv_access_area = 0.0;
        if( aptr->Selected() )
        {   
			elem = aptr->GetElemNo();
            if( !ElemDots[elem].count )
                AddElemDots(elem,density);
 
            Exclude = aptr;
            ptr = ElemDots[elem].dots;
            probe = ElemDots[elem].probe;
            count = ElemDots[elem].count;

            if( SolventDots )
            {   
				for( i=0; i<count; i++ )
                    if( TestDot( aptr->GetX() + probe[i].dx,
						aptr->GetY() + probe[i].dy,
						aptr->GetZ() + probe[i].dz, true ) )

						aptr->solv_access_area += (4*PI)/density;

                        dot_str->AddDot( aptr->GetX() + ptr[i].dx,
						aptr->GetY() + ptr[i].dy,
						aptr->GetZ() + ptr[i].dz,
						aptr->col );
            } 
			else
			{
                for( i=0; i<count; i++ )
				{
                    if( TestDot( aptr->GetX() + ptr[i].dx,
						aptr->GetY() + ptr[i].dy,
						aptr->GetZ() + ptr[i].dz, false ) )

						aptr->solv_access_area += (4*PI)/density;

						dot_str->AddDot( aptr->GetX() + ptr[i].dx,
						aptr->GetY() + ptr[i].dy,
						aptr->GetZ() + ptr[i].dz,
						aptr->col );
				}
			}
		}
	}

	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		aptr->UnSelect();
		if(aptr->solv_access_area > 0.001)
			aptr->Select();
	}

	DrawDots = True;
    FreeElemDots();
}
 
void 
HaMolView::DisplayObj3D()
{
 
  MolSet* pmset = GetMolSet();

  int ixadd=  pCanv->XRange()/2 ;
  int iyadd=  pCanv->YRange()/2 ;
  int izadd=  ZOffset() ;

  list<Object3D*>::iterator oitr;

  for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end(); oitr++)
  {
    if( (*oitr)->GetObjType() ==  OBJ3D_MOLECULE ||  (*oitr)->GetObjType() ==  OBJ3D_SURFACE || (*oitr)->GetObjType() ==  OBJ3D_DOT_SURFACE)
      continue;
		if((*oitr)->IsDisplayed())
			(*oitr)->DrawObj(this);
  }
}

void 
HaMolView::DisplayDotSurfaces()
{
    register int xi,yi,zi;
    register int i;
 
	MolSet* pmset = GetMolSet();

	int ixadd=  pCanv->XRange()/2 ;
	int iyadd=  pCanv->YRange()/2 ;
	int izadd=  ZOffset() ;				

	list<Object3D*>::iterator oitr;

	for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end(); oitr++)
	{
		if( (*oitr)->GetObjType() != OBJ3D_DOT_SURFACE )
			continue;

		DotStruct* ptr = (DotStruct*)(*oitr);

		int np = ptr->GetCount();
        for( i=0; i< np; i++ )
        {   
			double x_tr, y_tr, z_tr; 
			
			GetTransfCoord( ptr->dots[i].GetX(), ptr->dots[i].GetY(),ptr->dots[i].GetZ(),
				            x_tr, y_tr, z_tr);
			

            xi = (int)(x_tr*Scale) + ixadd;
            if( XValid_v(xi) )
            {   
				yi = (int)(y_tr*Scale) + iyadd;
                if( YValid_v(yi) )
                {   
					zi = (int)(z_tr*Scale)+ izadd;
                    if( ZValid_v(zi) )
                        pCanv->PlotDeepPoint(xi,yi,zi,ptr->dots[i].col);
                }
            }
        }

	}
}
 

/*==============================*/
/*  Ribbon & Cartoon Functions  */
/*==============================*/

 
static void CalculateVInten( Knot* ptr )
{
    double inten;
 
    if( !ptr->vsize ) // length of the normal vector 
        ptr->vsize = isqrt( (int)ptr->vnx*ptr->vnx +
                            (int)ptr->vny*ptr->vny +
                            (int)ptr->vnz*ptr->vnz ) + 1;
 
    inten = (double)ptr->vnz/ptr->vsize;  // Cosine of the angle of the normal to the vertical
 
    if( ptr->vnz < 0 ) inten = -inten;
 
    if( inten > 0.0 )
    {   
		ptr->vinten = (char)(ColourMask*inten);
    } 
	else 
		ptr->vinten = 0;
}
 

static void CalculateHInten( Knot* ptr )
{
    register double inten;
 
    /* The intensity of the sides of a protein cartoon
     * may be calculated using ptr->cx,cy,cz and this
     * should save interpolating ptr->hnx,hny,hnz!
     */
 
    if( !ptr->hsize )
        ptr->hsize = isqrt( (int)ptr->hnx*ptr->hnx +
                            (int)ptr->hny*ptr->hny +
                            (int)ptr->hnz*ptr->hnz ) + 1;
 
    inten = (double)ptr->hnz / ptr->hsize;

    if( ptr->hnz < 0 ) inten = -inten;
 
    if( inten > 0.0 )
    {   
		ptr->hinten = (char)(ColourMask*inten);
    } 
	else 
		ptr->hinten = 0;
}
 
 
void HaMolView::DisplayRibbon( HaChain* chain )
{
//	std::cerr << std::endl << " HaMolView::DisplayRibbon() pt 1 " << std::endl;
	HaResidue  *group;
	HaResidue  *gnext;
    HaAtom  *captr;
    HaAtom  *o1ptr;
    HaAtom  *o2ptr;
    HaAtom  *next;
 
    int prev;
	double wide;
    int col1,col2;
    int bx,by,bz;
    int dx,dy,dz;
    int arrow;
    int size;
 
    static Knot mid1, mid2, mid3;
    static Knot knot1, knot2;
 
    prev = False;

	vector<HaResidue*>::iterator gitr, gnext_itr;
	
	if( chain->res_arr.empty() )
	{
		return ;
	}

	gitr = chain->res_arr.begin();
	
    if( (*gitr)->IsProtein() )
    {   
		captr = (*gitr)->GetAtomByName("CA");
    } 
	else 
	{
		captr = (*gitr)->GetAtomByName("P");
	}
 
	gnext_itr = gitr;

	for(;;) // Cycle over residues of the chain
	{   
		gitr = gnext_itr;
		gnext_itr++;
		if( gnext_itr == chain->res_arr.end() )
			break;
		
		gnext = (*gnext_itr);
		group = (*gitr);
		if( gnext->IsProtein() )
		{   
			next =  gnext->GetAtomByName("CA");
			o1ptr = group->GetAtomByName("O");
		} 
		else /* Nucleic Acid */
		{   
			next =  gnext->GetAtomByName("P");
			o1ptr = gnext->GetAtomByName("O5*");
		}
		
		/* When not to have a control point! */
		if( !next || !captr || !o1ptr || (next->flag&BreakFlag) ||
			!( (group->flag|gnext->flag) & DrawKnotFlag) )
		{   
			group = gnext;
			captr = next;
			prev = False;
			continue;
		}
		
		knot2.tx = next->x - captr->x;
		knot2.ty = next->y - captr->y;
		knot2.tz = next->z - captr->z;
		
		if( group->IsProtein() )
		{   
			bx = o1ptr->x - captr->x;
			by = o1ptr->y - captr->y;
			bz = o1ptr->z - captr->z;
		} 
		else if( !group->GetAtomByName("O2") &&
			(o2ptr=group->GetAtomByName("O1P")) )
		{   /* Deoxyribonucleic Acid */
			o2ptr = group->GetAtomByName("O1P");
			bx = (o1ptr->x + o2ptr->x)/2 - captr->x;
			by = (o1ptr->y + o2ptr->y)/2 - captr->y;
			bz = (o1ptr->z + o2ptr->z)/2 - captr->z;
			
		} 
		else /* Ribonucleic Acid */
		{   
			bx = o1ptr->x - captr->x;
			by = o1ptr->y - captr->y;
			bz = o1ptr->z - captr->z;
		}
		
		knot2.px = (captr->x + next->x)/2;
		knot2.py = (captr->y + next->y)/2;
		knot2.pz = (captr->z + next->z)/2;
		
		/* c := a x b */
		knot2.vnx = knot2.ty*bz - knot2.tz*by;
		knot2.vny = knot2.tz*bx - knot2.tx*bz;
		knot2.vnz = knot2.tx*by - knot2.ty*bx;
		
		if( (group->struc & gnext->struc) & HelixFlag )
		{   /* Compensate for narrowing of helices! */
			size = isqrt(knot2.vnx*knot2.vnx +
				knot2.vny*knot2.vny +
				knot2.vnz*knot2.vnz);
			knot2.vsize = size;
			
			if( size )
			{   /* 1.00 Angstrom Displacement */
				wide = Scale;
				knot2.px += (int)(wide*knot2.vnx/size);  
				knot2.py += (int)(wide*knot2.vny/size);  
				knot2.pz += (int)(wide*knot2.vnz/size);  
			}
		} 
		else 
		{
			knot2.vsize = 0;
		}

		if( !(group->flag & gnext->flag & TraceFlag) ) // No trace ? level 1 if 
		{   /* d := c x a */
//			double cf = Scale*0.1;
			double cf = 96.0;
//			double cf = 50.0;
			dx = (int)((knot2.vny*knot2.tz - knot2.vnz*knot2.ty)/cf);
			dy = (int)((knot2.vnz*knot2.tx - knot2.vnx*knot2.tz)/cf);
			dz = (int)((knot2.vnx*knot2.ty - knot2.vny*knot2.tx)/cf);
			
			knot2.hsize = isqrt(dx*dx + dy*dy + dz*dz);
			
			/* Handle Carbonyl Oxygen Flip */
			if( prev && ((knot1.hnx*dx + knot1.hny*dy + knot1.hnz*dz)<0) )
			{   
				knot2.hnx = -dx;   knot2.vnx = -knot2.vnx;
				knot2.hny = -dy;   knot2.vny = -knot2.vny;
				knot2.hnz = -dz;   knot2.vnz = -knot2.vnz;
			} 
			else
			{   
				knot2.hnx = dx;
				knot2.hny = dy;
				knot2.hnz = dz;
			}
			
			arrow = False;
			if( group->flag&CartoonFlag )
			{   
				if( DrawBetaArrows && (group->struc&SheetFlag) &&
					!( gnext->struc&SheetFlag) )
				{   
					wide = (3*group->width)/2.0;
					arrow = True;
				} 
				else 
				{
					wide = group->width;
				}
			} 
			else if( group->flag & WideKnotFlag )
			{   /* Average Ribbon Width */
				if( gnext->flag & WideKnotFlag )
				{   
					wide = (group->width + gnext->width)/2.0;
				} 
				else if( gnext->flag & CartoonFlag )
				{   
					wide = gnext->width;
				} 
				else 
				{
					wide = group->width;
				}
			} 
			else
			{
				wide = gnext->width;
			}
			
			/* Set Ribbon Width */
			wide = wide*Scale;         // Ribbon width in screen coordinates (pixels)
			
			if( knot2.hsize && !arrow )
			{   
				size = knot2.hsize;
				knot2.wx = (int)(wide*knot2.hnx/size);
				knot2.wy = (int)(wide*knot2.hny/size);
				knot2.wz = (int)(wide*knot2.hnz/size);
				knot2.wide = (short)wide;
			} 
			else
			{   
				knot2.wide = 0;
				knot2.wx = 0;
				knot2.wy = 0;
				knot2.wz = 0;
			}
			
			if( group->flag & CartoonFlag ) // level 2 if Cartoon Display 
			{
				if( prev && (knot1.wide!=wide) && knot1.hsize )
				{   
					size = knot1.hsize;
					knot1.wx = (int)(wide*knot1.hnx/size);
					knot1.wy = (int)(wide*knot1.hny/size);
					knot1.wz = (int)(wide*knot1.hnz/size);
				}
			}
				
			if( (group->flag | gnext->flag) & CartoonFlag )
			{   
				CalculateVInten( &knot2 );
				CalculateHInten( &knot2 );
					
				size = knot2.vsize;
				wide = CartoonHeight*Scale;
				knot2.dx = (int)(wide*knot2.vnx/size);
				knot2.dy = (int)(wide*knot2.vny/size);
				knot2.dz = (int)(wide*knot2.vnz/size);
			} 
			else if( (group->flag| gnext->flag) & RibbonFlag )
			{
				CalculateVInten( &knot2 );
			}
		} // end No trace , end level 1 if
		
		if( !(col1 = group->col1) )
			col1 = captr->col;
		
		if( prev ) // if the residue has a previous residue in the chain( is not the first one) level 1 if
		{   
			/* Approximate spline segment with plane! */
			/* SolidRibbon( &knot1, &knot2, col1 );   */
			
			/* Calculate Hermite Spline Points */
			mid1.px = (int)(((int)54*knot1.px + (int)9*knot1.tx +
				(int)10*knot2.px - (int)3*knot2.tx)/64);
			mid1.py = (int)(((int)54*knot1.py + (int)9*knot1.ty +
				(int)10*knot2.py - (int)3*knot2.ty)/64);
			mid1.pz = (int)(((int)54*knot1.pz + (int)9*knot1.tz +
				(int)10*knot2.pz - (int)3*knot2.tz)/64);
			
			mid2.px = (int)(((int)4*knot1.px + knot1.tx +
				(int)4*knot2.px - knot2.tx)/8);
			mid2.py = (int)(((int)4*knot1.py + knot1.ty +
				(int)4*knot2.py - knot2.ty)/8);
			mid2.pz = (int)(((int)4*knot1.pz + knot1.tz +
				(int)4*knot2.pz - knot2.tz)/8);
			
			mid3.px = (int)(((int)10*knot1.px + (int)3*knot1.tx +
				(int)54*knot2.px - (int)9*knot2.tx)/64);
			mid3.py = (int)(((int)10*knot1.py + (int)3*knot1.ty +
				(int)54*knot2.py - (int)9*knot2.ty)/64);
			mid3.pz = (int)(((int)10*knot1.pz + (int)3*knot1.tz +
				(int)54*knot2.pz - (int)9*knot2.tz)/64);
			
			if( group->flag & TraceFlag ) // level 2 if, Trace Display
			{   
				wide = group->width*Scale;
				pCanv->ClipCylinder( knot1.px, knot1.py, knot1.pz, mid1.px, mid1.py, mid1.pz, col1, col1, (int)wide );
				pCanv->ClipCylinder( mid1.px, mid1.py, mid1.pz, mid2.px, mid2.py, mid2.pz, col1, col1, (int)wide );
				pCanv->ClipCylinder( mid2.px, mid2.py, mid2.pz, mid3.px, mid3.py, mid3.pz, col1, col1, (int)wide );
				pCanv->ClipCylinder( mid3.px, mid3.py, mid3.pz, knot2.px, knot2.py, knot2.pz, col1, col1, (int)wide );
			} 
			else if( group->flag & DotsFlag )
			{   
				wide = group->width*Scale;
				pCanv->ClipSphere(knot1.px,knot1.py,knot1.pz, (int)wide,col1);
				pCanv->ClipSphere(mid2.px, mid2.py, mid2.pz,  (int)wide,col1);
			} 
			else  // Strands, Ribbons or cartoon display level 2 if 
			{   /* Calculate Hermite Spline Widths */
				mid1.wx = (27*knot1.wx + 5*knot2.wx)/32;
				mid1.wy = (27*knot1.wy + 5*knot2.wy)/32;
				mid1.wz = (27*knot1.wz + 5*knot2.wz)/32;
				
				mid2.wx = (knot1.wx + knot2.wx)/2;
				mid2.wy = (knot1.wy + knot2.wy)/2;
				mid2.wz = (knot1.wz + knot2.wz)/2;
				
				mid3.wx = (5*knot1.wx + 27*knot2.wx)/32;
				mid3.wy = (5*knot1.wy + 27*knot2.wy)/32;
				mid3.wz = (5*knot1.wz + 27*knot2.wz)/32;
				
				/* Draw the Spline Segments */
				if( group->flag & (StrandFlag|DashStrandFlag) ) // level 3 if strand display
				{   
					if( !(col2 = group->col2) )
						col2 = captr->col;
					if( group->flag & StrandFlag )
					{   
						pCanv->StrandRibbon( &knot1, &mid1,  col1, col2 );
						pCanv->StrandRibbon( &mid1,  &mid2,  col1, col2 );
						pCanv->StrandRibbon( &mid2,  &mid3,  col1, col2 );
						pCanv->StrandRibbon( &mid3,  &knot2, col1, col2 );
					} 
					else /* group->flag & DashStrandFlag */
					{   
						pCanv->DashRibbon( &knot1, &mid1,  col1, col2 );
						pCanv->DashRibbon( &mid1,  &mid2,  col1, col2 );
						pCanv->DashRibbon( &mid2,  &mid3,  col1, col2 );
						pCanv->DashRibbon( &mid3,  &knot2, col1, col2 );
					}
				} 
				else /* Ribbon or Cartoon! */
				{   
					mid1.vsize = 0;
					mid1.vnx = (int)(((int)27*knot1.vnx +
						(int) 5*knot2.vnx)/32);
					mid1.vny = (int)(((int)27*knot1.vny +
						(int) 5*knot2.vny)/32);
					mid1.vnz = (int)(((int)27*knot1.vnz +
						(int) 5*knot2.vnz)/32);
					CalculateVInten( &mid1 );
					
					mid2.vsize = 0;
					mid2.vnx = (knot1.vnx + knot2.vnx)/2;
					mid2.vny = (knot1.vny + knot2.vny)/2;
					mid2.vnz = (knot1.vnz + knot2.vnz)/2;
					CalculateVInten( &mid2 );
					
					mid3.vsize = 0;
					mid3.vnx = (int)(((int) 5*knot1.vnx +
						(int)27*knot2.vnx)/32);
					mid3.vny = (int)(((int) 5*knot1.vny +
						(int)27*knot2.vny)/32);
					mid3.vnz = (int)(((int) 5*knot1.vnz +
						(int)27*knot2.vnz)/32);
					CalculateVInten( &mid3 );
					
					if( group->flag & RibbonFlag ) // level 4 if:  Display Ribbon
					{   
						if( group->struc & HelixFlag )
						{   
							if( !(col2 = group->col2) )
								col2 = captr->col;                          
						} 
						else
						{
							col2 = col1;
						}
						
						if( col1 != col2 )
						{   
							pCanv->SolidRibbon2( &knot1, &mid1,  col1, col2 );
							pCanv->SolidRibbon2( &mid1,  &mid2,  col1, col2 );
							pCanv->SolidRibbon2( &mid2,  &mid3,  col1, col2 );
							pCanv->SolidRibbon2( &mid3,  &knot2, col1, col2 );
						} 
						else
						{   
							pCanv->SolidRibbon( &knot1, &mid1,  col1 );
							pCanv->SolidRibbon( &mid1,  &mid2,  col1 );
							pCanv->SolidRibbon( &mid2,  &mid3,  col1 );
							pCanv->SolidRibbon( &mid3,  &knot2, col1 );
						}
					} 
					else // Display Cartoon 
					{   /* Calculate Spline Heights */
						wide = CartoonHeight*Scale;
						
						size = mid1.vsize;
						mid1.dx = (int)(wide*mid1.vnx/size);
						mid1.dy = (int)(wide*mid1.vny/size);
						mid1.dz = (int)(wide*mid1.vnz/size);
						
						size = mid2.vsize;
						mid2.dx = (int)(wide*mid2.vnx/size);
						mid2.dy = (int)(wide*mid2.vny/size);
						mid2.dz = (int)(wide*mid2.vnz/size);
						
						size = mid3.vsize;
						mid3.dx = (int)(wide*mid3.vnx/size);
						mid3.dy = (int)(wide*mid3.vny/size);
						mid3.dz = (int)(wide*mid3.vnz/size);
						
						/* Calculate Surface Intensity */
						mid1.hsize = 0;
						mid1.hnx = (int)(((int)27*knot1.hnx +
							(int) 5*knot2.hnx)/32);
						mid1.hny = (int)(((int)27*knot1.hny +
							(int) 5*knot2.hny)/32);
						mid1.hnz = (int)(((int)27*knot1.hnz +
							(int) 5*knot2.hnz)/32);
						CalculateHInten( &mid1 );
						
						mid2.hsize = 0;
						mid2.hnx = (knot1.hnx + knot2.hnx)/2;
						mid2.hny = (knot1.hny + knot2.hny)/2;
						mid2.hnz = (knot1.hnz + knot2.hnz)/2;
						CalculateHInten( &mid2 );
						
						mid3.hsize = 0;
						mid3.hnx = (int)(((int) 5*knot1.hnx +
							(int)27*knot2.hnx)/32);
						mid3.hny = (int)(((int) 5*knot1.hny +
							(int)27*knot2.hny)/32);
						mid3.hnz = (int)(((int) 5*knot1.hnz +
							(int)27*knot2.hnz)/32);
						CalculateHInten( &mid3 );
						
						pCanv->RectRibbon( &knot1, &mid1,  col1 );
						pCanv->RectRibbon( &mid1,  &mid2,  col1 );
						pCanv->RectRibbon( &mid2,  &mid3,  col1 );
						pCanv->RectRibbon( &mid3,  &knot2, col1 );
					} // end level 4 if (if ribbon or else if cartoon)
				} // end level 3 (if display group as strand/dashstrand or else as ribbon/cartoon
            } // end level 2 if, Trace Display or else strand, ribbon or cartoon display
        } // end level 1 if , the residue has a previous residue in the chain
		else if( group == chain->GetFirstRes() ) // if the residue a first residue in the chain
        {   
			knot1 = knot2;
            knot1.px = captr->x;
            knot1.py = captr->y;
            knot1.pz = captr->z;
			
            if( group->flag & RibbonFlag )
            {   
				pCanv->SolidRibbon( &knot1, &knot2, col1 );    
            } 
			else if( group->flag & RibbonFlag )
            {   
				pCanv->RectRibbon( &knot1, &knot2, col1 );
            } 
			else if( group->flag & StrandFlag )
            {   
				if( !(col2 = group->col2) )
                    col2 = captr->col;
                pCanv->StrandRibbon( &knot1,  &knot2, col1, col2 );
            } 
			else if( group->flag & DashStrandFlag )
            {   
				if( !(col2 = group->col2) )
                    col2 = captr->col;
                pCanv->DashRibbon( &knot1,  &knot2, col1, col2 );
            } 
			else if( group->flag & TraceFlag )
            {   
				pCanv->ClipCylinder( knot1.px, knot1.py, knot1.pz, knot2.px, knot2.py, knot2.pz, col1, col1, (int)(group->width*Scale) );
            } 
			else if( group->flag & DotsFlag )
            {   
				wide = (int)(group->width*Scale);
                pCanv->ClipSphere(knot1.px,knot1.py,knot1.pz,(int)wide,col1);
            }
            prev = True;
        } 
		else
		{
			prev = True;
		}
        group = gnext;
        captr = next;
		
        knot1 = knot2;
    } // end for cycle over residues of the chain
	
    if( prev ) // Plot the last residue of the chain
    {   
		if( !(col1 = group->col1) )
            col1 = captr->col;
		
        if( group->flag & CartoonFlag )
        {   /* Test for arrow head! */
            if( DrawBetaArrows && (group->struc&SheetFlag) )
            {   
				wide = (3*group->width)/2.0;
                knot2.px = captr->x + (knot2.tx/2);
                knot2.py = captr->y + (knot2.ty/2);
                knot2.pz = captr->z + (knot2.tz/2);
				
                arrow = True;
            } 
			else
            {   
				wide = group->width;
                knot2.px = captr->x;
                knot2.py = captr->y;
                knot2.pz = captr->z;
                arrow = False;
            }
			
            wide = Scale*wide;
            if( (knot1.wide!=wide) && knot1.hsize )
            {   
				size = knot1.hsize;
                knot1.wx = (int)((wide*knot1.hnx)/size);
                knot1.wy = (int)((wide*knot1.hny)/size);
                knot1.wz = (int)((wide*knot1.hnz)/size);
				
                if( !arrow )
                {   
					knot2.wx = knot1.wx;
                    knot2.wy = knot1.wy;
                    knot2.wz = knot1.wz;
                } 
				else
                {   
					knot2.wx = 0;
                    knot2.wy = 0;
                    knot2.wz = 0;
                }
            } 
			else if( arrow )
            {   
				knot2.wx = 0;
                knot2.wy = 0;
                knot2.wz = 0;
            }
			
            pCanv->RectRibbon( &knot1, &knot2, col1 );
        } 
		else /* !Cartoon */
        {   
			knot2.px = captr->x;
            knot2.py = captr->y;
            knot2.pz = captr->z;
			
            if( group->flag & RibbonFlag )
            {   
				pCanv->SolidRibbon( &knot1, &knot2, col1 );    
            } 
			else if( group->flag & StrandFlag )
            {   
				if( !(col2 = group->col2) )
                    col2 = captr->col;
                pCanv->StrandRibbon( &knot1,  &knot2, col1, col2 );
            } 
			else if( group->flag & DashStrandFlag )
            {   
				if( !(col2 = group->col2) )
                    col2 = captr->col;
                pCanv->DashRibbon( &knot1,  &knot2, col1, col2 );
            } 
			else if( group->flag & TraceFlag )
            {   
				pCanv->ClipCylinder( knot1.px, knot1.py, knot1.pz,
					knot2.px, knot2.py, knot2.pz,
					col1, col1, (int)(group->width*Scale) );
            } 
			else if( group->flag & DotsFlag )
            {   
				wide = (int)(group->width*Scale);
                pCanv->ClipSphere(knot1.px,knot1.py,knot1.pz,(int)wide,col1);
                pCanv->ClipSphere(knot2.px,knot2.py,knot2.pz,(int)wide,col1);
            }   
		}
	}  
}

void HaMolView::DisplayETBestPath()
{
	MolSet* pmset = GetMolSet();
	if(pmset == NULL) return;
	ETCouplMod* ptr_et_coupl_mod = pmset->GetETCouplMod(false);
	if(ptr_et_coupl_mod == NULL) return;
	if(ptr_et_coupl_mod->best_path.empty()) return;
	int i;

	HaColor path_col(255,255,0); // yellow

	int sc = path_col.cidx;
	int dc = path_col.cidx;

	int irad = (int) (Scale*0.2); // radius of the cylinder 

	for(i = 0; i < ptr_et_coupl_mod->best_path.size(); i++)
	{
		int isrc=  ptr_et_coupl_mod->best_path[i].source;
		int idest= ptr_et_coupl_mod->best_path[i].destination;
		
		if( isrc < 0 || idest < 0 ) continue; 

		HaAtom* aptr1= ptr_et_coupl_mod->nodes[ isrc ];
		HaAtom* aptr2= ptr_et_coupl_mod->nodes[ idest ];
		if(aptr1 != NULL && aptr2 != NULL)
		{
			if(aptr1->IsBonded(*aptr2))
			{
//				pCanv->ClipCylinder(aptr1->x, aptr1->y, aptr1->z, aptr2->x, aptr2->y, aptr2->z,
//					                aptr1->col,aptr2->col,irad);
				pCanv->ClipCylinder(aptr1->x, aptr1->y, aptr1->z, aptr2->x, aptr2->y, aptr2->z,
					                sc,dc,irad);

			}
			else
			{
				double dist = Vec3D::CalcDistance(aptr1,aptr2,ANGSTROM_U);
				
				int ndiv = (int)(dist/0.3);
				
				if(ndiv < 3) ndiv = 3;
				
				int half_ndiv = ndiv/2;
				ndiv = half_ndiv*2 + 1;
				
				double xstep = ((double)(aptr2->x - aptr1->x))/ndiv;
				double ystep = ((double)(aptr2->y - aptr1->y))/ndiv;
				double zstep = ((double)(aptr2->z - aptr1->z))/ndiv;
				
				int x1 ; double addx = xstep; 
				int y1 ; double addy = ystep; 
				int z1 ; double addz = zstep;
				
				for(int i= 0; i < half_ndiv; i++)
				{
					x1 = aptr1->x + (int)addx;
					y1 = aptr1->y + (int)addy;
					z1 = aptr1->z + (int)addz;

					pCanv->ClipCylinder(x1, y1, z1, 
				          (int)(x1 + xstep/2.0), (int)(y1 + ystep/2.0), (int)(z1 + zstep/2.0), 
						   sc,dc, (int)(irad/1.5));
					
					addx += 2*xstep; addy += 2*ystep; addz += 2*zstep;	
				}				
//				pCanv->ClipDashVector(aptr1->x,aptr1->y,aptr1->z,aptr2->x,aptr2->y,aptr2->z,
//					aptr1->col,aptr2->col);
				
			}
		}
	}

}

void
HaMolView::DisplayContourSurf()
{
	MolSet* pmset= GetMolSet();
	if(pmset == NULL) return;

	list<Object3D*>::iterator oitr;
	for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end(); oitr++)
	{
		if( (*oitr)->GetObjType() != OBJ3D_SURFACE )
			continue;

		if( !(*oitr)->IsDisplayed() )
			continue;

		HaDisplayedSurface* sptr = (HaDisplayedSurface*)(*oitr);
		double transp = sptr->transparency;
		int plot_transp = 0;
		if( transp < -0.0001 || transp > 1.0001) transp = 0.0;
		if( transp > 0.0001) plot_transp = 1;

//		PrintLog(" Transparency = %9.3f \n",transp);

		int i;
		int x1, x2, x3, y1, y2, y3, z1, z2, z3;
		int cl1, cl2, cl3;
		float dx,dy,dz;

		int ixadd=  pCanv->XRange()/2 ;
		int iyadd=  pCanv->YRange()/2 ;
		int izadd=  ZOffset() ;				

		int numverts = sptr->GetNumVerts();
		if(sptr->colors.size() < numverts) continue;
		if(numverts < 3) continue;

		bool subst_vert_1= true;
		static Poly p;

		bool have_norms = (sptr->norms.num_cols() == sptr->GetNumVerts()) && 
				              (sptr->norms.num_rows() == 3 );

		double zn1, zn2, xn3, yn3, zn3;
		
		int ntr = sptr->GetNumTr();
		int j;
		for(i = 1; i <= ntr; i++) 
		{
			for(j = 1; j <= 3; j++)
			{
				int idx_v =  sptr->tr_indx(j,i) + 1;
				dx = sptr->verts(1,idx_v );
				dy = sptr->verts(2,idx_v );
				dz = sptr->verts(3,idx_v );

				double x_tr, y_tr, z_tr;
				GetTransfCoord( dx, dy, dz, x_tr, y_tr, z_tr);
			
				x3 = (int)(x_tr * Scale) + ixadd;
				y3 = (int)(y_tr * Scale) + iyadd;
				z3 = (int)(z_tr * Scale) + izadd;

				dx = sptr->norms(1,idx_v);
				dy = sptr->norms(2,idx_v);
				dz = sptr->norms(3,idx_v);

				xn3 = dx*Rot(1,1)  + dy*Rot(1,2)  + dz*Rot(1,3);
				yn3 = dx*Rot(2,1)  + dy*Rot(2,2)  + dz*Rot(2,3);
				zn3 = dx*Rot(3,1)  + dy*Rot(3,2)  + dz*Rot(3,3);
			
				zn3 = fabs(zn3);
			
				cl3 = sptr->colors(idx_v);
				if( j == 1)
				{
					x1 = x3;
					y1 = y3;
					z1 = z3;
					zn1 = zn3;
					cl1 = cl3;
				}
				if( j == 2)
				{
					x2 = x3;
					y2 = y3;
					z2 = z3;
					zn2 = zn3;
					cl2 = cl3;
				}
			}
			
			if(DrawSolidSurfaces)     // Solid Triangle
			{
//				if( cl3 > 80 && zn3 > 0.1)
//				{
//					x_tr =0.0;
//				}
				
				p.count=3;
				p.v[0].x = x1;  
				p.v[0].y = y1; 
				p.v[0].z = z1;
				p.v[1].x = x2;  
				p.v[1].y = y2;  
				p.v[1].z = z2;
				p.v[2].x = x3;
				p.v[2].y = y3;  
				p.v[2].z = z3;
				if( have_norms) 
				{
					p.v[0].inten = cl1 + (int)(ColourDepth*zn1);
					p.v[1].inten = cl2 + (int)(ColourDepth*zn2); 
					p.v[2].inten = cl3 + (int)(ColourDepth*zn3);
				}
					
				pCanv->ClipPolygon( &p,transp );
			}
			else // WireFrame triangle
			{
				pCanv->ClipTwinVector(x1,y1,z1,x2,y2,z2, cl1,cl2);
				pCanv->ClipTwinVector(x2,y2,z2,x3,y3,z3, cl2,cl3);
				pCanv->ClipTwinVector(x1,y1,z1,x3,y3,z3, cl1,cl3);
			}
		}
	}
}




