/*! \file canvas3d.cpp

    Classes to plot 3D images in a window

   \author Igor Kurnikov
   \date 1998-2002 

   Derived from:

   pixutils.c
   RasMol2 Molecular Graphics
   Roger Sayle, August 1995
   Version 2.6
*/


#define CANVAS3D_CPP
#define PIXUTILS

#if defined(_WIN32) | defined(TWIN32)
#include <malloc.h>
#endif
#ifdef sun386
#include <stdlib.h>
#endif

#include <math.h>

#include "haio.h"
#include "canvas3d.h"
#include "hamolview.h"

#include "abstree.h"
#include "font.h"

unsigned char* Canvas3D::LookUp[MAXRAD]= {0};
unsigned char  Canvas3D::Array[MAXTABLE]={0};
unsigned int  Canvas3D::ColConstTable[MAXRAD]={0};
unsigned int* Canvas3D::ColConst=Canvas3D::ColConstTable;

ColorVal Canvas3D::Lut[]={0};
int Canvas3D::empty_lut_idx = 0;

Canvas3D::Canvas3D()
{
	FBuffer=NULL;
	DBuffer=NULL;
	m_SplineCount = 5;
	SetFontSize(8);
	m_XRange=m_YRange=0;
	m_SplineCount=0;
	m_UseSlabPlane=False;
	m_SlabValue=0;
	m_SlabMode=0;
	m_SlabInten=0;
	m_SliceValue=0;
	m_ImageSize=0;
	m_ImageRadius=0;
	m_ZOffset=0;
}

Canvas3D::~Canvas3D()
{
	DestroyImage();
	DestroyDBuffer();
}

void Canvas3D::resize(int new_XRange, int new_YRange)
{
//	cerr << " Canvas3D::resize() new_XRange & new_YRange = " << 
//		    new_XRange << "  " << new_YRange << endl;

	m_XRange= new_XRange; 
	m_YRange= new_YRange;
	
	/* Ensure int Aligned */
	int dx;
	if( dx = m_XRange%4 )
		m_XRange += 4-dx;
}

ColorVal* Canvas3D::AllocImage()
{
	if(FBuffer) free(FBuffer);
	size_t size= (size_t) XRange() * YRange() *sizeof(ColorVal);
	FBuffer = (ColorVal*)malloc( size+32 );
	return FBuffer;
}

void Canvas3D::DestroyImage()
{
	if(FBuffer) free(FBuffer);
}

short* Canvas3D::AllocDBuffer()
{
	size_t size = (size_t)XRange()*YRange()*sizeof(short)+32;
    if( DBuffer ) free( DBuffer );
    DBuffer = (short*)malloc( size );
    if( !DBuffer ) 
	{
		cerr << " Error allocating depth buffer" << endl;
		exit(1);
	}
	return DBuffer;
}

void Canvas3D::DestroyDBuffer()
{
	if( DBuffer ) free( DBuffer );
}

//#ifdef _WIN32
#	define InvertY(y) ((y))
//#else
//#	define InvertY(y) (-(y))
//#endif

/* Sutherland-Cohen Line Clipping Macros */
#define BitAbove    0x01
#define BitBelow    0x02
#define BitRight    0x04
#define BitLeft     0x08
#define BitFront    0x10

#define Reject(x,y)   ((x)&(y))     // if both non-zero is != 0 - reject pair of poits
#define Accept(x,y)   (!((x)|(y)))  // != 0 if both x and y are zero otherwise = 0 

#define SETPIXEL(dptr,fptr,d,c)    if( (d) > *(dptr) )              \
                                   {   *(dptr) = (d);               \
                                       *(fptr) = (c);               \
                                   }

/* These define light source position */
#define LightDot(x,y,z)  ((x)+(y)+(z)+(z))


/* Note: DrawCylinderCaps currently employs an
 *       extremely crude hack to avoid stripes
 *       appearing aint cylinders.
 */
//#define ARCSIZE  2048
#define ARCSIZE  20000

static ArcEntry  *ArcAcPtr;
static ArcEntry  *ArcDnPtr;

static ArcEntry ArcAc[ARCSIZE];
static ArcEntry ArcDn[ARCSIZE];


static int ClipStatus;


int Canvas3D::OutCode(int x, int y, int z)
// Determine if the point out of view and setup corresponding flags
// in result BitAbove, BitBelow, BitLeft, BitRight, BitFront 
{
    register int result;

    if( y<0 )
    {   
		result = BitAbove;
    } 
	else if( y>=View.ymax )
    {   
		result = BitBelow;
    } 
	else 
		result = 0;

    if( x<0 )
    {   
		result |= BitLeft;
    } 
	else if( x>=View.xmax )
        result |= BitRight;

    if( !ZValid(z) )
        result |= BitFront;
    return result;
}


void Canvas3D::PlotPoint(int x,int y,int z,int col)
{
    ColorVal  *fptr;
    short  *dptr;
    int offset;

    /* SETPIXEL(dptr,fptr,z,Lut[col]); */

    offset = (int)y*View.yskip+x;
    dptr = View.dbuf+offset;
    if( z > *dptr )
    {   
		fptr = View.fbuf+offset;
        *fptr = Lut[col];
        *dptr = z;
    }
}


void Canvas3D::ClipPoint(int x, int y,  int z, int col)
{
    ColorVal  *fptr;
    short  *dptr;
    int offset;

    if( XValid(x) && YValid(y) && ZValid(z) )
    {   /* PlotPoint(x,y,z,col); */
        offset = (int)y*View.yskip+x;
        dptr = View.dbuf+offset;
        if( z > *dptr )
        {   fptr = View.fbuf+offset;
            *fptr = Lut[col];
            *dptr = z;
        }
    }
}


void 
Canvas3D::PlotDeepPoint(int x, int y, int z, int col)
{
    int offset;
    ColorVal  *fptr;
    short  *dptr;
    int inten;

    offset = (int)y*View.yskip+x;
    dptr = View.dbuf+offset;

    if( z > *dptr )
    {  
	   fptr = View.fbuf+offset;
       inten = (ColourDepth*(z+ m_ImageRadius - m_ZOffset ))/m_ImageSize;
       *fptr = Lut[col+inten];
       *dptr = z;
    }
}

void 
Canvas3D::ClipDeepPoint(int x, int y, int z, int col)
{
    int offset;
    ColorVal  *fptr;
    short  *dptr;
    int inten;

    if( XValid(x) && YValid(y) && ZValid(z) )
    {   /* PlotDeepPoint(x,y,z,col); */
        offset = (int)y*View.yskip+x;
        dptr = View.dbuf+offset;

        if( z > *dptr )
        {  fptr = View.fbuf+offset;
           inten = (ColourDepth*(z+ m_ImageRadius - m_ZOffset))/m_ImageSize;
           *fptr = Lut[col+inten];
           *dptr = z;
        }
    }
}


/* Macros for Bresenhams Line Drawing Algorithm */
#define CommonStep(s)  z1 += zrate; SETPIXEL(dptr,fptr,z1,c);     \
                       if( (zerr+=dz)>0 ) { zerr-=(s); z1+=iz; }

#define XStep  { if((err+=dy)>0) { fptr+=ystep; dptr+=ystep; err-=dx; } \
                 fptr+=ix; dptr+=ix; x1+=ix; CommonStep(dx); }

#define YStep  { if((err+=dx)>0) { fptr+=ix; dptr+=ix; err-=dy; } \
                 fptr+=ystep; dptr+=ystep; y1+=iy; CommonStep(dy); }
                     

void 
Canvas3D::DrawTwinLine(int x1,int y1,int z1,int x2,int y2,int z2,
					   int col1,int col2)
{
    int offset;
    ColorVal  *fptr;
    short  *dptr;
    int zrate, zerr;
    int ystep,err;
    int ix,iy,iz;
    int dx,dy,dz;
    int mid;
    ColorVal c;
	
    c = Lut[col1];
	
    offset = (int)y1*View.yskip + x1;
    fptr = View.fbuf+offset;
    dptr = View.dbuf+offset;
	
    SETPIXEL(dptr,fptr,z1,c);
	
    dx = x2-x1;  dy = y2-y1; 
    if( !dx && !dy ) return;
    dz = z2-z1;
	
    if( dy<0 ) 
    {   
		ystep = -View.yskip;
        dy = -dy; 
        iy = -1;
    } 
	else
    {   
		ystep = View.yskip;
        iy = 1;
    }
	
    if( dx<0 ) 
    {   
		dx = -dx;
        ix = -1;
    } 
	else 
		ix = 1;
	
    if( dz<0 ) 
    {   
		dz = -dz;
        iz = -1;
    } 
	else 
		iz = 1;
	
    if( dx>dy )
    {   
		if( dz >= dx )
        {   
			zrate = dz/dx;
            dz -= dx*zrate;
        } 
		else 
			zrate = 0;
        err = zerr = -(dx>>1);
		
        if( col1 != col2 )
        {   
			mid = (x1+x2)>>1;
            while( x1!=mid ) XStep;
            c = Lut[col2];
        }
        while( x1!=x2 ) XStep;
		
    } 
	else
    {   
		if( dz >= dy )
        {   
			zrate = dz/dy;
            dz -= dy*zrate;
        } 
		else 
			zrate = 0;
        err = zerr = -(dy>>1);
		
        if( col1 != col2 )
        {   
			mid = (y1+y2)>>1;
            while( y1!=mid ) YStep;
            c = Lut[col2];
        }
        while( y1!=y2 ) YStep;
    }
}


void 
Canvas3D::ClipLine(int x1,int y1,int z1,int x2,int y2,int z2,int col)
{
    int code1,code2;
    int delta,rest;
    int temp;
	
    while( True )
    {   
		code1 = OutCode(x1,y1,z1);
        code2 = OutCode(x2,y2,z2);
        if( Reject(code1,code2) ) return;
        if( Accept(code1,code2) ) break;
		
        if( !code1 )
        {   
			temp=x1; x1=x2; x2=temp;
            temp=y1; y1=y2; y2=temp;
            temp=z1; z1=z2; z2=temp;
            code1 = code2;
        }
		
        if( code1 & BitAbove )
        {   
			delta = y2-y1;
            x1 += (int)(((int)y1*(x1-x2))/delta);  
            z1 += (int)(((int)y1*(z1-z2))/delta);
            y1 = 0;
        } 
		else if( code1 & BitLeft )
        {   
			delta = x2-x1;
            y1 += (int)(((int)x1*(y1-y2))/delta);
            z1 += (int)(((int)x1*(z1-z2))/delta);
            x1 = 0;
        } 
		else if( code1 & BitRight )
        {   
			delta = x2-x1;
            temp=View.xmax-1; rest=temp-x1;
            y1 += (int)(((int)rest*(y2-y1))/delta);
            z1 += (int)(((int)rest*(z2-z1))/delta);
            x1 = temp;
        } 
		else if( code1 & BitBelow )
        {   
			delta = y2-y1;
            temp=View.ymax-1; rest=temp-y1;
            x1 += (int)(((int)rest*(x2-x1))/delta);
            z1 += (int)(((int)rest*(z2-z1))/delta);
            y1 = temp;
        } 
		else /* SLAB */
        {   
			delta = z2-z1;
            rest = (m_SlabValue-1)-z1;
            x1 += (int)(((int)rest*(x2-x1))/delta);
            y1 += (int)(((int)rest*(y2-y1))/delta);
            z1 =   m_SlabValue-1;
        }
    }
    DrawTwinLine(x1,y1,z1,x2,y2,z2,col,col);
}


void 
Canvas3D::ClipTwinLine(int x1,int y1,int z1,int x2,int y2,int z2,
				  int col1,int col2)
{
    int xmid,ymid,zmid;
    int code1,code2;


    if( col1!=col2 )
    {   code1 = OutCode(x1,y1,z1);
        code2 = OutCode(x2,y2,z2);
        if( !Reject(code1,code2) )
        {   if( !Accept(code1,code2) )
            {  xmid = (x1+x2)/2;
               ymid = (y1+y2)/2;
               zmid = (z1+z2)/2;
               ClipLine(x1,y1,z1,xmid,ymid,zmid,col1);
               ClipLine(xmid,ymid,zmid,x2,y2,z2,col2);
            } else
               DrawTwinLine(x1,y1,z1,x2,y2,z2,col1,col2);
        }
    } else
        ClipLine(x1,y1,z1,x2,y2,z2,col1);
}



/* Macros for 3D Bresenhams Vector Algorithm */
#define CommonVectStep(s)  z1 += zrate;   c1 += crate;                    \
                           SETPIXEL(dptr,fptr,z1,Lut[col+c1]);            \
                           if( (zerr+=dz)>0 ) { zerr -= (s); z1 += iz; }  \
                           if( (cerr+=dc)>0 ) { cerr -= (s); c1 += iz; }

#define XVectStep  { if((err+=dy)>0) { fptr+=ystep; dptr+=ystep; err-=dx; } \
                     fptr+=ix; dptr+=ix; x1+=ix; CommonVectStep(dx); }

#define YVectStep  { if((err+=dx)>0) { fptr+=ix; dptr+=ix; err-=dy; } \
                     fptr+=ystep; dptr+=ystep; y1+=iy; CommonVectStep(dy); }


void Canvas3D::DrawTwinVector(int x1,int y1,int z1,int x2,int y2,int z2, int col1,int col2)
{
    int offset;
    ColorVal  *fptr;
    short  *dptr;
    int dx,dy,dz,dc;
    int crate, cerr;
    int zrate, zerr;
    int ystep,err;
    int ix,iy,iz;
    int col, mid;
    int c1, c2;
	
    c1 = (ColourDepth*(z1+ m_ImageRadius- m_ZOffset))/m_ImageSize;
    c2 = (ColourDepth*(z2+ m_ImageRadius- m_ZOffset))/m_ImageSize;
	
    offset = (int)y1*View.yskip + x1;
    fptr = View.fbuf+offset;
    dptr = View.dbuf+offset;
	
    SETPIXEL(dptr,fptr,z1,Lut[col1+c1]);
	
    dx = x2 - x1;  dy = y2 - y1;
    dz = z2 - z1;  dc = c2 - c1;
    if( !dx && !dy ) return;
	
    if( dy<0 ) 
    {   
		ystep = -View.yskip;
        dy = -dy; 
        iy = -1; 
    } 
	else
    {   
		ystep = View.yskip;
        iy = 1;
    }
	
    if( dx<0 ) 
    {   
		dx = -dx; 
        ix = -1; 
    } 
	else 
		ix = 1;
	
    iz = (dz<0)? -1 : 1;
	
    if( dx>dy )
    {   
		if( dz >= dx )
        {   
			zrate = dz/dx;
            dz -= dx*zrate;
        } 
		else 
			zrate = 0;
		
        if( dc >= dx )
        {   
			crate = dc/dx;
            dc -= dx*crate;
        } 
		else 
			crate = 0;
		
        err = zerr = cerr = -(dx>>1);
        col = col1;
		
        if( dz<0 )
        {   
			dz = -dz;
            dc = -dc;
        }
        
        if( col1 != col2 )
        {   
			mid = (x1+x2)>>1;
            while( x1!=mid ) XVectStep;
            col = col2;
        }
        while( x1!=x2 ) XVectStep;
    } 
	else
    {   
		if( dz >= dy )
        {   
			zrate = dz/dy;
            dz -= dy*zrate;
        } 
		else 
			zrate = 0;
		
        if( dc >= dy )
        {   
			crate = dc/dy;
            dc -= dy*crate;
        } 
		else 
			crate = 0;
		
        err = zerr = cerr = -(dy>>1);
        col = col1;
		
        if( dz<0 )
        {   
			dz = -dz;
            dc = -dc;
        }
		
        if( col1 != col2 )
        {   
			mid = (y1+y2)>>1;
            while( y1!=mid ) YVectStep;
            col=col2;
        }
        while( y1!=y2 ) YVectStep;
    }
}

void Canvas3D::ClipVector(int x1,int y1,int z1, int x2,int y2,int z2,int col)
{
    int code1,code2;
    int delta,rest;
    int temp;

    code1 = OutCode(x1,y1,z1);
    code2 = OutCode(x2,y2,z2);

    while( True )
    {   
		if( Accept(code1,code2) ) break;
        if( Reject(code1,code2) ) return;

        if( !code1 )
        {   
			code1 = code2; code2 = 0;
            temp=x1; x1=x2; x2=temp;
            temp=y1; y1=y2; y2=temp;
            temp=z1; z1=z2; z2=temp;
        }

        if( code1 & BitAbove )
        {   
			delta = y2-y1;
            x1 += (int)(((int)y1*(x1-x2))/delta);  
            z1 += (int)(((int)y1*(z1-z2))/delta);
            y1 = 0;
        } 
		else if( code1 & BitLeft )
        {   
			delta = x2-x1;
            y1 += (int)(((int)x1*(y1-y2))/delta);
            z1 += (int)(((int)x1*(z1-z2))/delta);
            x1 = 0;
        } 
		else if( code1 & BitRight )
        {   
			delta = x2-x1;
            temp=View.xmax-1; rest=temp-x1;
            y1 += (int)(((int)rest*(y2-y1))/delta);
            z1 += (int)(((int)rest*(z2-z1))/delta);
            x1 = temp;
        } 
		else if( code1 & BitBelow )
        {   
			delta = y2-y1;
            temp=View.ymax-1; rest=temp-y1;
            x1 += (int)(((int)rest*(x2-x1))/delta);
            z1 += (int)(((int)rest*(z2-z1))/delta);
            y1 = temp;
        } 
		else /* SLAB */
        {   
			delta = z2-z1;
            rest = (m_SlabValue-1)-z1;
            x1 += (int)(((int)rest*(x2-x1))/delta);
            y1 += (int)(((int)rest*(y2-y1))/delta);
            z1 = m_SlabValue-1;
        }
        code1 = OutCode(x1,y1,z1);
    }
    DrawTwinVector(x1,y1,z1,x2,y2,z2,col,col);
}


void Canvas3D::ClipTwinVector(int x1,int y1,int z1, int x2,int y2,int z2,int col1,int col2)
{
    int xmid,ymid,zmid;
    int code1,code2;
	
    if( col1!=col2 )
    {   
		code1 = OutCode(x1,y1,z1); // Determine if points out of view 
        code2 = OutCode(x2,y2,z2); // OutCode() != 0 if is out of view
        if( !Reject(code1,code2) )
        {   
			if( !Accept(code1,code2) ) // if one of the point is in view recursively calling 
            {                          // the function plot the vector to the border of the view
				xmid = (x1+x2)/2;
				ymid = (y1+y2)/2;
				zmid = (z1+z2)/2;
				ClipVector(x1,y1,z1,xmid,ymid,zmid,col1);
				ClipVector(xmid,ymid,zmid,x2,y2,z2,col2);
            } 
			else
				DrawTwinVector(x1,y1,z1,x2,y2,z2,col1,col2);
        }
    } 
	else
        ClipVector(x1,y1,z1,x2,y2,z2,col1);
}


void Canvas3D::ClipDashVector(int x1,int y1,int z1, int x2,int y2,int z2,int col1,int col2)
{
    int offset;
    ColorVal  *fptr;
    short  *dptr;
    int ix,iy,iz,ic;
    int dx,dy,dz,dc;
    int crate, cerr;
    int zrate, zerr;
    int ystep,err;
    int col, mid;
    int c1, c2;
    int count;

    if( (x1==x2) && (y1==y2) ) return;
    if( Reject(OutCode(x1,y1,z1),OutCode(x2,y2,z2)) )
        return;

    c1 = (ColourDepth*(z1+ m_ImageRadius- m_ZOffset))/m_ImageSize;
    c2 = (ColourDepth*(z2+ m_ImageRadius- m_ZOffset))/m_ImageSize;

    dx = x2 - x1;  dy = y2 - y1;
    dz = z2 - z1;  dc = c2 - c1;

    offset = (int)y1*View.yskip + x1;
    fptr = View.fbuf+offset;
    dptr = View.dbuf+offset;
    count = 0;

    ystep = View.yskip;
    ix = iy = iz = ic = 1;
    if( dy<0 ) { dy = -dy; iy = -1; ystep = -ystep; }
    if( dx<0 ) { dx = -dx; ix = -1; }
    if( dz<0 ) { dz = -dz; iz = -1; }
    if( dc<0 ) { dc = -dc; ic = -1; }


    if( dx>dy )
    {   if( x2<x1 )
        {   mid = col1;
            col1 = col2;
            col2 = mid;
        }
        if( dz >= dx )
        {   zrate = dz/dx;
            dz -= dx*zrate;
        } else zrate = 0;

        if( dc >= dx )
        {   crate = dc/dx;
            dc -= dx*crate;
        } else crate = 0;

        err = zerr = cerr = -(dx>>1);
        mid = (x1+x2)/2;

        while( x1!=x2 )
        {   if( XValid(x1) && YValid(y1) && ZValid(z1) )
            {   if( count<2 )
                {   col = (x1<mid)? col1 : col2;
                    SETPIXEL(dptr,fptr,z1,Lut[col+c1]);
                    count++;
                } else if( count==3 )
                {   count = 0;
                } else count++;
            }

            if( (err+=dy)>0 )
            {   err -= dx;
                fptr+=ystep;
                dptr+=ystep;
                y1+=iy;
            }

            if( (zerr+=dz)>0 )
            {   zerr -= dx;
                z1 += iz;
            }

            if( (cerr+=dc)>0 )
            {   cerr -= dx;
                c1 += ic;
            }

            fptr+=ix; dptr+=ix; x1+=ix;
            z1 += zrate;   c1 += crate;
        }
    } else
    {   if( y1>y2 )
        {   mid = col1;
            col1 = col2;
            col2 = mid;
        }

        if( dz >= dy )
        {   zrate = dz/dy;
            dz -= dy*zrate;
        } else zrate = 0;

        if( dc >= dy )
        {   crate = dc/dy;
            dc -= dy*crate;
        } else crate = 0;

        err = zerr = cerr = -(dy>>1);
        mid = (y1+y2)/2;

        
        while( y1!=y2 )
        {   if( XValid(x1) && YValid(y1) && ZValid(z1) )
            {   if( count<2 )
                {   col = (y1<mid)? col1 : col2;
                    SETPIXEL(dptr,fptr,z1,Lut[col+c1]);
                    count++;
                } else if( count==3 )
                {   count = 0;
                } else count++;
            }

            if( (err+=dx)>0 )
            {   err-=dy;
                fptr+=ix;
                dptr+=ix;
                x1+=ix;
            }

            if( (zerr+=dz)>0 )
            {   zerr -= dy;
                z1 += iz;
            }

            if( (cerr+=dc)>0 )
            {   cerr -= dy;
                c1 += ic;
            }

            fptr+=ystep; dptr+=ystep; y1+=iy;
            z1 += zrate;   c1 += crate;
        }
    }
}


/* m_SplineCount is either 1, 2, 3, 4, 5 or 9! */

void Canvas3D::StrandRibbon( Knot* src, Knot* dst, int col1, int col2 )
{
    int hsx, hsy, hsz;
    int hdx, hdy, hdz;
    int qsx, qsy, qsz;
    int qdx, qdy, qdz;
    int col;

    if( m_SplineCount != 4 )
    {   if( m_SplineCount == 1 ) 
        {   ClipVector( src->px, src->py, src->pz,
                        dst->px, dst->py, dst->pz, col2 );
            return;
        } else if( m_SplineCount != 2 )
            ClipVector( src->px, src->py, src->pz,
                        dst->px, dst->py, dst->pz, col1 );

        ClipVector( src->px+src->wx, src->py+src->wy, src->pz+src->wz,
                    dst->px+dst->wx, dst->py+dst->wy, dst->pz+dst->wz, col2 );
        ClipVector( src->px-src->wx, src->py-src->wy, src->pz-src->wz,
                    dst->px-dst->wx, dst->py-dst->wy, dst->pz-dst->wz, col2 );
        if( m_SplineCount<=3 ) return;

        hsx = src->wx/2;  hsy = src->wy/2;  hsz = src->wz/2;
        hdx = dst->wx/2;  hdy = dst->wy/2;  hdz = dst->wz/2;

        ClipVector( src->px+hsx, src->py+hsy, src->pz+hsz,
                    dst->px+hdx, dst->py+hdy, dst->pz+hdz, col1 );
        ClipVector( src->px-hsx, src->py-hsy, src->pz-hsz,
                    dst->px-hdx, dst->py-hdy, dst->pz-hdz, col1 );
        if( m_SplineCount==5 ) 
            return;
        col = col1;
    } else /* m_SplineCount == 4 */
    {   hsx = src->wx/2;  hsy = src->wy/2;  hsz = src->wz/2;
        hdx = dst->wx/2;  hdy = dst->wy/2;  hdz = dst->wz/2;
        col = col2;
    }

    qsx = hsx/2;  qsy = hsy/2;  qsz = hsz/2;
    qdx = hdx/2;  qdy = hdy/2;  qdz = hdz/2;

    ClipVector( src->px+hsx+qsx, src->py+hsy+qsy, src->pz+hsz+qsz,
                dst->px+hdx+qdx, dst->py+hdy+qdy, dst->pz+hdz+qdz, col );
    ClipVector( src->px+hsx-qsx, src->py+hsy-qsy, src->pz+hsz-qsz,
                dst->px+hdx-qdx, dst->py+hdy-qdy, dst->pz+hdz-qdz, col1 );
    ClipVector( src->px-hsx+qsx, src->py-hsy+qsy, src->pz-hsz+qsz,
                dst->px-hdx+qdx, dst->py-hdy+qdy, dst->pz-hdz+qdz, col1 );
    ClipVector( src->px-hsx-qsx, src->py-hsy-qsy, src->pz-hsz-qsz,
                dst->px-hdx-qdx, dst->py-hdy-qdy, dst->pz-hdz-qdz, col );
}


void Canvas3D::DashRibbon( Knot* src, Knot* dst, int col1, int col2 )
{
    int hsx, hsy, hsz;
    int hdx, hdy, hdz;
    int qsx, qsy, qsz;
    int qdx, qdy, qdz;
    int col;

    if( m_SplineCount != 4 )
    {   if( m_SplineCount == 1 ) 
        {   ClipDashVector( src->px, src->py, src->pz,
                            dst->px, dst->py, dst->pz, col2, col2 );
            return;
        } else if( m_SplineCount != 2 )
            ClipDashVector( src->px, src->py, src->pz,
                            dst->px, dst->py, dst->pz, col1, col1 );

        ClipDashVector(src->px+src->wx,src->py+src->wy,src->pz+src->wz,
                       dst->px+dst->wx,dst->py+dst->wy,dst->pz+dst->wz,
                       col2,col2);
        ClipDashVector(src->px-src->wx,src->py-src->wy,src->pz-src->wz,
                       dst->px-dst->wx,dst->py-dst->wy,dst->pz-dst->wz,
                       col2,col2);
        if( m_SplineCount<=3 ) return;

        hsx = src->wx/2;  hsy = src->wy/2;  hsz = src->wz/2;
        hdx = dst->wx/2;  hdy = dst->wy/2;  hdz = dst->wz/2;

        ClipDashVector( src->px+hsx, src->py+hsy, src->pz+hsz,
                        dst->px+hdx, dst->py+hdy, dst->pz+hdz, col1, col1 );
        ClipDashVector( src->px-hsx, src->py-hsy, src->pz-hsz,
                        dst->px-hdx, dst->py-hdy, dst->pz-hdz, col1, col1 );
        if( m_SplineCount==5 ) 
            return;
        col = col1;
    } else /* m_SplineCount == 4 */
    {   hsx = src->wx/2;  hsy = src->wy/2;  hsz = src->wz/2;
        hdx = dst->wx/2;  hdy = dst->wy/2;  hdz = dst->wz/2;
        col = col2;
    }

    qsx = hsx/2;  qsy = hsy/2;  qsz = hsz/2;
    qdx = hdx/2;  qdy = hdy/2;  qdz = hdz/2;

    ClipDashVector(src->px+hsx+qsx,src->py+hsy+qsy,src->pz+hsz+qsz,
                   dst->px+hdx+qdx,dst->py+hdy+qdy,dst->pz+hdz+qdz,col,col);
    ClipDashVector(src->px+hsx-qsx,src->py+hsy-qsy,src->pz+hsz-qsz,
                   dst->px+hdx-qdx,dst->py+hdy-qdy,dst->pz+hdz-qdz,col1,col1);
    ClipDashVector(src->px-hsx+qsx,src->py-hsy+qsy,src->pz-hsz+qsz,
                   dst->px-hdx+qdx,dst->py-hdy+qdy,dst->pz-hdz+qdz,col1,col1);
    ClipDashVector(src->px-hsx-qsx,src->py-hsy-qsy,src->pz-hsz-qsz,
                   dst->px-hdx-qdx,dst->py-hdy-qdy,dst->pz-hdz-qdz,col,col);
}


#ifndef PIXUTILS  /* Unused Function */
void 
Canvas3D::OutLinePolygon( Poly* p )
{
    register int i;

    for( i=0; i<p->count-1; i++ )
         ClipLine( p->v[i].x, p->v[i].y, p->v[i].z, 
                   p->v[i+1].x, p->v[i+1].y, p->v[i+1].z,
                   p->v[i].inten);
    ClipLine( p->v[i].x, p->v[i].y, p->v[i].z,
              p->v[0].x, p->v[0].y, p->v[0].z,
              p->v[i].inten);
}
#endif


#ifndef PIXUTILS
void 
Canvas3D::DrawPolygon( Poly* p )
{
    static Edge lft, rgt;
    register Edge *pmin, *pmax;
    register ColorVal  *fbase;
    register short  *dbase;
    register short  *dptr;
    register int offset;

    register int dz,di;
    register int z,inten;
    register int ri,li,ry,ly;
    register int xmin,xmax;
    register int dy,ymin;
    register int top,rem;
    register int x,y,i;

    /* Find top vertex */
    top = 0;  
    ymin = p->v[0].y;
    for( i=1; i<p->count; i++ )
       if( p->v[i].y < ymin )
       {   ymin = p->v[i].y;
           top = i;
       }

    rem = p->count;
    ly = ry = y = ymin;
    li = ri = top;

    offset = (int)y*View.yskip;
    fbase = View.fbuf+offset;
    dbase = View.dbuf+offset;

    while( rem )
    {   while( ly<=y && rem )
        {   i = li-1; if( i<0 ) i=p->count-1;
            if( p->v[i].y > y )
            {   dy = p->v[i].y - ly;
                lft.di = (((int)(p->v[i].inten - p->v[li].inten))<<16)/dy;
                lft.dx = (((int)(p->v[i].x - p->v[li].x))<<16)/dy;
                lft.dz = (((int)(p->v[i].z - p->v[li].z))<<16)/dy;

                lft.i = ((int)p->v[li].inten)<<16;
                lft.x = ((int)p->v[li].x)<<16;
                lft.z = ((int)p->v[li].z)<<16;
            }
            ly = p->v[i].y;
            rem--;  li = i;
        }

        while( ry<=y && rem )
        {   i = ri+1; if( i>=p->count ) i = 0;
            if( p->v[i].y > y )
            {   dy = p->v[i].y - ry;
                rgt.di = (((int)(p->v[i].inten - p->v[ri].inten))<<16)/dy;
                rgt.dx = (((int)(p->v[i].x - p->v[ri].x))<<16)/dy;
                rgt.dz = (((int)(p->v[i].z - p->v[ri].z))<<16)/dy;

                rgt.i = ((int)p->v[ri].inten)<<16;
                rgt.x = ((int)p->v[ri].x)<<16;
                rgt.z = ((int)p->v[ri].z)<<16;
            }
            ry = p->v[i].y;
            rem--; ri = i;
        }


        ymin = MinFun(ly,ry);
        
        while( y<ymin )
        {   if( lft.x < rgt.x )
            {   pmin = &lft;
                pmax = &rgt;
            } else
            {   pmin = &rgt;
                pmax = &lft;
            }

            xmax = (int)(pmax->x>>16)+1;
            xmin = (int)(pmin->x>>16);

            di = (int)((pmax->i-pmin->i)/(xmax-xmin));
            dz = (int)((pmax->z-pmin->z)/(xmax-xmin));
            inten = pmin->i;  
            z = pmin->z;

            dptr = dbase+xmin;
            for( x=xmin; x<xmax; x++ )
            {   if( (int)(z>>16) > *dptr )
                {   fbase[x] = Lut[(int)(inten>>16)];
                    *dptr = (int)(z>>16);
                }
                inten += di;
                z += dz;
                dptr++;
            }

            lft.x += lft.dx;  rgt.x += rgt.dx;
            lft.z += lft.dz;  rgt.z += rgt.dz;
            lft.i += lft.di;  rgt.i += rgt.di;
            dbase += View.yskip;
            fbase += View.yskip;
            y++;
        }
    }
}
#endif


void Canvas3D::ClipPolygon( Poly* p, double transp )
//! 
//! To plot a polygon specify number of verticies,
//! coordinates of vertices and color number
//! accounting (for decreasing intensity of the colour with
//! inclination of a surface element
//!
{
    static Edge lft, rgt;
    Edge *pmin, *pmax;
    ColorVal  *fbase;
    short  *dbase;
    short  *dptr;
    int offset;

    register int dz;
	int dr,dg,db; 
	int cr,cg,cb;
    int z;
    int ri,li,ry,ly;
    int xmin,xmax;
    int dy,ymin;
    int top,rem;
    int x,y,i;

	int plot_transp = FALSE;
	double transp_comp = 0.0;
	if(transp > 0.001 && transp < 1.0001)
	{
		plot_transp = TRUE;
		transp_comp = 1.0 - transp;
	}
	

    /* Reject Clip Polygon */
    if( m_UseSlabPlane )
        for( i=0; i<p->count; i++ )
            if( p->v[i].z >= m_SlabValue )
                return;

    /* Find top vertex */
    top = 0;  
    ymin = p->v[0].y;
    for( i=1; i<p->count; i++ )
       if( p->v[i].y < ymin )
       {   ymin = p->v[i].y;
           top = i;
       }

    rem = p->count;
    ly = ry = y = ymin;
    li = ri = top;

    if( y<0 )
    {   rem--;

        while( ly<=0 && rem )
        {   
			i = li-1; 
			if( i<0 ) i=p->count-1;
            if( p->v[i].y > 0 )
            {   
				dy = p->v[i].y - ly;
//              lft.di = (((int)(p->v[i].inten - p->v[li].inten))<<16);
				lft.dr = ( (int)(RComp(Lut[p->v[i].inten]) - RComp(Lut[p->v[li].inten]))<<16)/dy;
				lft.dg = ( (int)(GComp(Lut[p->v[i].inten]) - GComp(Lut[p->v[li].inten]))<<16)/dy;
				lft.db = ( (int)(BComp(Lut[p->v[i].inten]) - BComp(Lut[p->v[li].inten]))<<16)/dy;
				
                lft.dx = (((int)(p->v[i].x - p->v[li].x))<<16)/dy;
                lft.dz = (((int)(p->v[i].z - p->v[li].z))<<16)/dy;

//                lft.i = ((int)p->v[li].inten)<<16;
                lft.r = ( (int) RComp(Lut[p->v[li].inten]))<<16;
                lft.g = ( (int) GComp(Lut[p->v[li].inten]))<<16;
                lft.b = ( (int) BComp(Lut[p->v[li].inten]))<<16;

                lft.x = ((int)p->v[li].x)<<16;
                lft.z = ((int)p->v[li].z)<<16;
            } 
			else 
			{
				rem--;
			}
            ly = p->v[i].y;
            li = i;
        }

        while( ry<=0 && rem )
        {   
			i = ri+1; 
			if( i>=p->count ) i = 0;
            if( p->v[i].y > 0 )
            {   
				dy = p->v[i].y - ry;
//                rgt.di = (((int)(p->v[i].inten - p->v[ri].inten))<<16)/dy;
				rgt.dr = ( (int)(RComp(Lut[p->v[i].inten]) - RComp(Lut[p->v[ri].inten]))<<16)/dy;
				rgt.dg = ( (int)(GComp(Lut[p->v[i].inten]) - GComp(Lut[p->v[ri].inten]))<<16)/dy;
				rgt.db = ( (int)(BComp(Lut[p->v[i].inten]) - BComp(Lut[p->v[ri].inten]))<<16)/dy;

                rgt.dx = (((int)(p->v[i].x - p->v[ri].x))<<16)/dy;
                rgt.dz = (((int)(p->v[i].z - p->v[ri].z))<<16)/dy;

//                rgt.i = ((int)p->v[ri].inten)<<16;
                rgt.r = ((int) RComp(Lut[p->v[ri].inten]))<<16;
                rgt.g = ((int) GComp(Lut[p->v[ri].inten]))<<16;
                rgt.b = ((int) BComp(Lut[p->v[ri].inten]))<<16;

                rgt.x = ((int)p->v[ri].x)<<16;
                rgt.z = ((int)p->v[ri].z)<<16;
            } 
			else
			{
				rem--;
			}
            ry = p->v[i].y;
            ri = i;
        }

        fbase = View.fbuf;
        dbase = View.dbuf;
        y = 0;
    } 
	else /* y >= 0 */
    {   
		offset = (int)y*View.yskip;
        fbase = View.fbuf+offset;
        dbase = View.dbuf+offset;
    }

    while( rem )
    {   
		while( ly<=y && rem )
        {   
			i = li-1; 
			if( i<0 ) i=p->count-1;
            if( p->v[i].y > y )
            {   
				dy = p->v[i].y - ly;
//				lft.di = (((int)(p->v[i].inten - p->v[li].inten))<<16)/dy;
				lft.dr = ( (int)(RComp(Lut[p->v[i].inten]) - RComp(Lut[p->v[li].inten]))<<16)/dy;
				lft.dg = ( (int)(GComp(Lut[p->v[i].inten]) - GComp(Lut[p->v[li].inten]))<<16)/dy;
				lft.db = ( (int)(BComp(Lut[p->v[i].inten]) - BComp(Lut[p->v[li].inten]))<<16)/dy;

				lft.dx = (((int)(p->v[i].x - p->v[li].x))<<16)/dy;
				lft.dz = (((int)(p->v[i].z - p->v[li].z))<<16)/dy;
				
//				lft.i = ((int)p->v[li].inten)<<16;
                lft.r = ((int) RComp(Lut[p->v[li].inten]))<<16;
                lft.g = ((int) GComp(Lut[p->v[li].inten]))<<16;
                lft.b = ((int) BComp(Lut[p->v[li].inten]))<<16;

				lft.x = ((int)p->v[li].x)<<16;
				lft.z = ((int)p->v[li].z)<<16;
            }
            ly = p->v[i].y;
            rem--;  li = i;
        }
		
        while( ry<=y && rem )
        {   
			i = ri+1; if( i>=p->count ) i = 0;
            if( p->v[i].y > y )
            {   
				dy = p->v[i].y - ry;
//                rgt.di = (((int)(p->v[i].inten - p->v[ri].inten))<<16)/dy;
				rgt.dr = ( (int)(RComp(Lut[p->v[i].inten]) - RComp(Lut[p->v[ri].inten]))<<16)/dy;
				rgt.dg = ( (int)(GComp(Lut[p->v[i].inten]) - GComp(Lut[p->v[ri].inten]))<<16)/dy;
				rgt.db = ( (int)(BComp(Lut[p->v[i].inten]) - BComp(Lut[p->v[ri].inten]))<<16)/dy;
				
                rgt.dx = (((int)(p->v[i].x - p->v[ri].x))<<16)/dy;
                rgt.dz = (((int)(p->v[i].z - p->v[ri].z))<<16)/dy;
				
//                rgt.i = ((int)p->v[ri].inten)<<16;
                rgt.r = ((int) RComp(Lut[p->v[ri].inten]))<<16;
                rgt.g = ((int) GComp(Lut[p->v[ri].inten]))<<16;
                rgt.b = ((int) BComp(Lut[p->v[ri].inten]))<<16;

                rgt.x = ((int)p->v[ri].x)<<16;
                rgt.z = ((int)p->v[ri].z)<<16;
            }
            ry = p->v[i].y;
            rem--; ri = i;
        }
		
        ymin = MinFun(ly,ry);
        if( ymin>View.ymax )
        {   
			ymin = View.ymax;
            rem = 0;
        }
        
        while( y<ymin )
        {   
			if( lft.x < rgt.x )
            {   
				pmin = &lft;
                pmax = &rgt;
            } 
			else
            {   
				pmin = &rgt;
                pmax = &lft;
            }

            xmax = (int)(pmax->x>>16)+1;
            xmin = (int)(pmin->x>>16);

            if( (xmin<View.xmax) && (xmax>=0) )
            {   
//				di = (int)((pmax->i-pmin->i)/(xmax-xmin));
				dr = (int)((pmax->r-pmin->r)/(xmax-xmin));
				dg = (int)((pmax->g-pmin->g)/(xmax-xmin));
				db = (int)((pmax->b-pmin->b)/(xmax-xmin));

                dz = (int)((pmax->z-pmin->z)/(xmax-xmin));
                if( xmin<0 )
                {   
//					inten = pmin->i - xmin*di;
					cr = pmin->r - xmin*dr;
					cg = pmin->g - xmin*dg;
					cb = pmin->b - xmin*db;

                    z = pmin->z - xmin*dz;
                    xmin = 0;
                } 
				else /* xmin >= 0 */
                {   
//					inten = pmin->i; 
					cr = pmin->r;
					cg = pmin->g;
					cb = pmin->b;
					
                    z = pmin->z;
                }

                if( xmax>=View.xmax )
                    xmax = View.xmax;

                dptr = dbase+xmin;
                for( x=xmin; x<xmax; x++ )
                {   
//					if( ((int)dptr) %16 != 0 )  // for transparent surfaces
//					{
						if( (int)(z>>16) > *dptr )
						{   							
							uint_4 r=  cr >> 16; 
							uint_4 g = cg >> 16; 
							uint_4 b = cb >> 16;
							if(plot_transp)
							{
								uint_4 r_old = RComp(fbase[x]);
								uint_4 g_old = GComp(fbase[x]);
								uint_4 b_old = BComp(fbase[x]);
								r = (uint_4)(r_old*transp + r*transp_comp);	
                                g = (uint_4)(g_old*transp + g*transp_comp);
								b = (uint_4)(b_old*transp + b*transp_comp);
							}
								
							fbase[x] = (((r<<8)|g)<<8 ) | b;
							
							*dptr = (int)(z>>16);
						}
//						inten += di;
						cr += dr;
						cg += dg;
						cb += db;

						z += dz;
						dptr++;
//					}
                }
            }

            lft.x += lft.dx;  rgt.x += rgt.dx;
            lft.z += lft.dz;  rgt.z += rgt.dz;

//            lft.i += lft.di;  rgt.i += rgt.di;
			  lft.r += lft.dr;  rgt.r += rgt.dr;
			  lft.g += lft.dg;  rgt.g += rgt.dg;
			  lft.b += lft.db;  rgt.b += rgt.db;

            dbase += View.yskip;
            fbase += View.yskip;
            y++;
        }
    }
}


void 
Canvas3D::SolidRibbon( Knot* src, Knot* dst, int col )
{
    static Poly p;

    p.v[0].x = src->px+src->wx;  
    p.v[0].y = src->py+src->wy;  
    p.v[0].z = src->pz+src->wz;
	p.v[0].inten = src->vinten+col;

    p.v[1].x = dst->px+dst->wx;  
    p.v[1].y = dst->py+dst->wy;  
    p.v[1].z = dst->pz+dst->wz;
	p.v[1].inten = dst->vinten+col;

    p.v[2].x = dst->px-dst->wx;
    p.v[2].y = dst->py-dst->wy;  
    p.v[2].z = dst->pz-dst->wz;
	p.v[2].inten = dst->vinten+col;

    p.v[3].x = src->px-src->wx;  
    p.v[3].y = src->py-src->wy;  
    p.v[3].z = src->pz-src->wz;
	p.v[3].inten = src->vinten+col;

    p.count = 4;
    /* OutLinePolygon( &p ); */
    ClipPolygon( &p );
}


void 
Canvas3D::SolidRibbon2( Knot* src, Knot* dst, int col1, int col2 )
{
    static Poly p;

    p.count = 3;
    p.v[0].x = src->px+src->wx;  
    p.v[0].y = src->py+src->wy;  
    p.v[0].z = src->pz+src->wz;
    p.v[1].x = dst->px-dst->wx;  
    p.v[1].y = dst->py-dst->wy;  
    p.v[1].z = dst->pz-dst->wz;


    p.v[2].x = dst->px+dst->wx;
    p.v[2].y = dst->py+dst->wy;  
    p.v[2].z = dst->pz+dst->wz;

	p.v[0].inten = src->vinten+col1;
	p.v[1].inten = dst->vinten+col1;
	p.v[2].inten = dst->vinten+col1;
	
    /* OutLinePolygon( &p ); */
    ClipPolygon( &p );

	p.v[2].x = src->px-src->wx;  
	p.v[2].y = src->py-src->wy;  
	p.v[2].z = src->pz-src->wz;
	
	p.v[0].inten = src->vinten+col2;
	p.v[1].inten = dst->vinten+col2;
	p.v[2].inten = src->vinten+col2;
	
    /* OutLinePolygon( &p ); */
    ClipPolygon( &p );
}


void 
Canvas3D::RectRibbon( Knot* src, Knot* dst, int col )
{
    static Poly p;

    p.count = 4;

	p.v[0].inten = src->vinten+col;
	p.v[1].inten = dst->vinten+col;
	p.v[2].inten = dst->vinten+col;
	p.v[3].inten = src->vinten+col;		
	

    /* Top Surface */
    p.v[0].x = src->px+src->wx+src->dx;  
    p.v[0].y = src->py+src->wy+src->dy;  
    p.v[0].z = src->pz+src->wz+src->dz;

    p.v[1].x = dst->px+dst->wx+dst->dx;  
    p.v[1].y = dst->py+dst->wy+dst->dy;  
    p.v[1].z = dst->pz+dst->wz+dst->dz;

    p.v[2].x = dst->px-dst->wx+dst->dx;
    p.v[2].y = dst->py-dst->wy+dst->dy;  
    p.v[2].z = dst->pz-dst->wz+dst->dz;

    p.v[3].x = src->px-src->wx+src->dx;  
    p.v[3].y = src->py-src->wy+src->dy;  
    p.v[3].z = src->pz-src->wz+src->dz;
    ClipPolygon( &p );

    /* Bottom Surface */
    p.v[0].x = src->px+src->wx-src->dx;  
    p.v[0].y = src->py+src->wy-src->dy;  
    p.v[0].z = src->pz+src->wz-src->dz;

    p.v[1].x = dst->px+dst->wx-dst->dx;  
    p.v[1].y = dst->py+dst->wy-dst->dy;  
    p.v[1].z = dst->pz+dst->wz-dst->dz;

    p.v[2].x = dst->px-dst->wx-dst->dx;
    p.v[2].y = dst->py-dst->wy-dst->dy;  
    p.v[2].z = dst->pz-dst->wz-dst->dz;

    p.v[3].x = src->px-src->wx-src->dx;  
    p.v[3].y = src->py-src->wy-src->dy;  
    p.v[3].z = src->pz-src->wz-src->dz;
    ClipPolygon( &p );

	p.v[0].inten = src->hinten+col;
	p.v[1].inten = dst->hinten+col;
	p.v[2].inten = dst->hinten+col;
	p.v[3].inten = src->hinten+col;

    /* Left Surface */
    p.v[0].x = src->px+src->wx+src->dx;  
    p.v[0].y = src->py+src->wy+src->dy;  
    p.v[0].z = src->pz+src->wz+src->dz;

    p.v[1].x = dst->px+dst->wx+dst->dx;  
    p.v[1].y = dst->py+dst->wy+dst->dy;  
    p.v[1].z = dst->pz+dst->wz+dst->dz;

    p.v[2].x = dst->px+dst->wx-dst->dx;
    p.v[2].y = dst->py+dst->wy-dst->dy;  
    p.v[2].z = dst->pz+dst->wz-dst->dz;

    p.v[3].x = src->px+src->wx-src->dx;  
    p.v[3].y = src->py+src->wy-src->dy;  
    p.v[3].z = src->pz+src->wz-src->dz;
    ClipPolygon( &p );

    /* Right Surface */
    p.v[0].x = src->px-src->wx+src->dx;  
    p.v[0].y = src->py-src->wy+src->dy;  
    p.v[0].z = src->pz-src->wz+src->dz;

    p.v[1].x = dst->px-dst->wx+dst->dx;  
    p.v[1].y = dst->py-dst->wy+dst->dy;  
    p.v[1].z = dst->pz-dst->wz+dst->dz;

    p.v[2].x = dst->px-dst->wx-dst->dx;
    p.v[2].y = dst->py-dst->wy-dst->dy;  
    p.v[2].z = dst->pz-dst->wz-dst->dz;

    p.v[3].x = src->px-src->wx-src->dx;  
    p.v[3].y = src->py-src->wy-src->dy;  
    p.v[3].z = src->pz-src->wz-src->dz;
    ClipPolygon( &p );
}


int 
Canvas3D::TestSphere( register int x, register int y, 
					   register int z, register int rad )
{
    register int temp;

    ClipStatus = 0;

    if( m_UseSlabPlane )
    {   
		if( z-rad>= m_SlabValue )
            return( False );

        if( z+rad>= m_SlabValue )
        {   
			if( m_SlabMode )
            {   
				ClipStatus |= BitFront;
            } 
			else return( False );
        } 
		else if( m_SlabMode==SlabSection )
            return( False );
    }

    temp = x+rad;
    if( temp<0 ) return( False );
    if( temp>=View.xmax ) ClipStatus |= BitRight;

    temp = x-rad;
    if( temp>=View.xmax ) return( False );
    if( temp<0 ) ClipStatus |= BitLeft;

    temp = y+rad;
    if( temp<0 ) return( False );
    if( temp>=View.ymax ) ClipStatus |= BitBelow;

    temp = y-rad;
    if( temp>=View.ymax ) return( False );
    if( temp<0 ) ClipStatus |= BitAbove;

    return True;
}


#define CalcInten(dz)    inten = LightDot(dx,InvertY(dy),(dz))

#define UpdateAcross(dz)    \
        depth = (dz)+z;                    \
        if( depth > *dptr )                \
        {   *dptr = depth;                 \
            fptr = fold+dx;                \
            CalcInten((dz));               \
            if( inten>0 )                  \
            {      inten = (int)((inten*ColConst[rad])>>ColBits); \
                   *fptr = Lut[col+inten]; \
            } else *fptr = Lut[col];       \
        }                                  \
        dptr++;  dx++;


#define UpdateLine  \
        dx = -wide;                   \
        dptr = dold-wide;             \
        tptr = LookUp[wide]+wide;     \
        while( dx<0 ) { UpdateAcross(*tptr); tptr--; }       \
        do { UpdateAcross(*tptr); tptr++; } while(dx<=wide); \
        dold += View.yskip;  fold += View.yskip;             \
        dy++;


void Canvas3D::DrawSphere(int x,int y,int z,int rad,int col)
{
    ColorVal  *fptr,  *fold;
    short  *dptr,  *dold;
    unsigned char  *tptr;

    int offset;
    int depth,wide,inten;
    int dx,dy;

    offset = (int)(y-rad)*View.yskip + x;
    fold=View.fbuf+offset;  
    dold=View.dbuf+offset;

    dy = -rad;
    while( dy<0 ) 
    {   
		wide = LookUp[rad][-dy]; 
        UpdateLine; 
    }

    do { 
        wide = LookUp[rad][dy];  
        UpdateLine; 
    } while( dy<=rad );
}


void Canvas3D::ClipSphere(int x,int y,int z,int rad,int col)
{
    register ColorVal  *fptr,  *fold;
    register short  *dptr,  *dold;
	
    register int lastx,lasty,dx,dy,dz;
    register int depth,wide,inten,side;
    register int crad,cwide,temp;
    register int offset;
	
	
    /* Visibility Tests */
    if( !TestSphere(x,y,z,rad) )
        return;
	
    if( !ClipStatus )
    {   
		DrawSphere(x,y,z,rad,col);
        return;
    }
	
    if( ClipStatus&BitAbove )
    {   
		dy = -y;
        fold = View.fbuf + x;
        dold = View.dbuf + x;
    } 
	else
    {   
		dy = -rad;
        offset = (int)(y+dy)*View.yskip+x;
        fold = View.fbuf + offset;
        dold = View.dbuf + offset;
    }
	
    if( ClipStatus&BitBelow )
    {   
		lasty = (View.ymax-1)-y;
    } 
	else 
		lasty = rad;
	
	
    side = (View.xmax-1)-x;
    /* No Slab Plane Clipping */
    if( !(ClipStatus&BitFront) )
    {   
		while( dy<=lasty )
        {   
			wide = LookUp[rad][AbsFun(dy)];
            lastx = MinFun(wide,side);
            dx = - MinFun(wide,x);
            dptr = dold + dx;
			
            while( dx<=lastx )
            {   
				dz = LookUp[wide][AbsFun(dx)];
                UpdateAcross(dz);
            }
            dold += View.yskip;
            fold += View.yskip;
            dy++;
        }
        return;
    }
	
	
    dz = m_SlabValue-z;
    crad = LookUp[rad][AbsFun(dz)];
	
	
    if( (z > m_SlabValue) || (m_SlabMode==SlabSection) )
    {   
		if( crad<lasty ) lasty = crad;
        if( -crad>dy ) 
        {   
			dy = -crad;
            offset = (int)(y+dy)*View.yskip+x;
            fold = View.fbuf + offset;
            dold = View.dbuf + offset;
        }
    }
	
    while( dy<=lasty )
    {   
		temp = AbsFun(dy);
        wide = LookUp[rad][temp];
        lastx = MinFun(wide,side);
        dx = - MinFun(x,wide);
        dptr = dold + dx;
		
        if( temp<=crad )
        {   
			cwide = LookUp[crad][temp];
		while( dx<=lastx )
		{   
			temp = AbsFun(dx);
			if( temp<=cwide )
			{    /* Slab Plane Clipping Modes */
				switch( m_SlabMode )
				{   
				case( SlabFinal ):
					fold[dx] = Lut[col+ m_SlabInten];
					*dptr = m_SliceValue;
					break;
					
				case( SlabHollow ):
					dz = LookUp[wide][temp];
					depth = z - dz;
					if( depth>*dptr )
					{   
						*dptr = depth;
						inten = LightDot(-dx,-InvertY(dy),dz);
						
						if( inten>0 )
						{   
							inten=(int)( (inten*ColConst[rad])
								>>(ColBits+1));
							fold[dx] = Lut[col+inten];
						} 
						else 
							fold[dx] = Lut[col];
					}
					break;
					
				case( SlabSection ):
				case( SlabClose ):
					dz = m_SlabValue-z;
					depth = dx*dx+dy*dy+dz*dz+ m_SliceValue;
					if( (*dptr< m_SliceValue) || (depth<*dptr) )
					{   
						fold[dx] = Lut[col+ m_SlabInten];
						*dptr = depth;
					}
					break;
				}
				dptr++;  dx++;
			} 
			else if( (z < m_SlabValue) && ( m_SlabMode!=SlabSection) )
			{    
				dz = LookUp[wide][temp];
				UpdateAcross(dz);
			} 
			else
			{   
				dptr++;  dx++;
			}
		}
        } 
		else /* Slabless ScanLine */
            while( dx<=lastx )
            {   
				dz = LookUp[wide][AbsFun(dx)];
                UpdateAcross(dz);
            }
			
			dold += View.yskip;
			fold += View.yskip;
			dy++;
    }
}


/* Function Prototypes */
static void DrawArcDn( short *, ColorVal *, int, int );
static void DrawArcAc( short *, ColorVal *, int, int );
static void ClipArcDn( short *, ColorVal *, int, int, int, int );
static void ClipArcAc( short *, ColorVal *, int, int, int, int );


void Canvas3D::DrawArcAc(short* dbase, ColorVal* fbase, int z, int c)
{
    ArcEntry  *ptr;
    short  *dptr;
    short depth;

    for( ptr=ArcAc; ptr<ArcAcPtr; ptr++ )
    {   
		dptr = dbase+ptr->offset;  depth = ptr->dz+z;
        SETPIXEL(dptr,fbase+ptr->offset,depth,Lut[ptr->inten+c]);
    }
}

void Canvas3D::DrawArcDn(short* dbase, ColorVal* fbase, int z, int c)
{
    ArcEntry  *ptr;
    short  *dptr;
    short depth;

    for( ptr=ArcDn; ptr<ArcDnPtr; ptr++ )
    {   
		dptr = dbase+ptr->offset;  depth = ptr->dz+z;
        SETPIXEL(dptr,fbase+ptr->offset,depth,Lut[ptr->inten+c]);
    }
}


void Canvas3D::DrawCylinderCaps( int x1,int y1,int z1, int x2, int y2, int z2, int c1, int c2, int rad )
{
    // std::cerr << std::endl << " Canvas3D::DrawCylinderCaps() pt 1 " << std::endl;
    short  *dold,  *dptr;
    ColorVal  *fold;
#ifndef PIXUTILS
    int ax,ay,ix,iy;
    int zrate,lz;
#endif
    int offset,temp,end;
    int inten,absx;
    int wide,depth;
    int dx,dy,dz;
    int lx,ly;

    lx = x2-x1;
    ly = y2-y1;

#ifndef PIXUTILS
    lz = z2-z1;
    if( ly>0 ) 
	{ 
		ay = ly; iy = 1; 
	}
    else 
	{ 
		ay = -ly; iy = -1; 
	}
    if( lx>0 ) 
	{ 
		ax = lx; ix = 1; 
	}
    else 
	{ 
		ax = -lx; 
		ix = -1; 
	}
    zrate = lz/MaxFun(ax,ay);
#endif

    end = (int)ly*View.yskip+lx;
    temp = (int)y1*View.yskip+x1;
    fold = View.fbuf+temp;
    dold = View.dbuf+temp;

    ArcAcPtr = ArcAc;
    ArcDnPtr = ArcDn;
	int idx_arc = 0;

    temp = (int)-(rad*View.yskip);
    for( dy= -rad; dy<=rad; dy++ )
    {   
		wide = LookUp[rad][AbsFun(dy)];

        for( dx= -wide; dx<=wide; dx++ )
        {   
			absx = AbsFun(dx);
            dz = LookUp[wide][absx];
            CalcInten(dz);
            if( inten>0 )
            {   
				inten = (int)((inten*ColConst[rad])>>ColBits);
            } 
			else 
				inten = 0;
            offset = temp+dx;

            if( XValid(x1+dx) && YValid(y1+dy) )
            {   
				dptr = dold+offset; depth = z1+dz;
                SETPIXEL(dptr,fold+offset,depth,Lut[c1 + inten]);
            }

            if( XValid(x2+dx) && YValid(y2+dy) )
            {   
				dptr = dold+(offset+end); depth = z2+dz;
                SETPIXEL(dptr,fold+(offset+end),depth,Lut[c2+inten]);
            }

#ifndef PIXUTILS
            k1 = AbsFun(dx+ix); 
            k2 = AbsFun(dx-ix);

            if( ((k1>wide)||(dz>= LookUp[wide][k1]-zrate)) &&
                ((k2>wide)||(dz > LookUp[wide][k2]+zrate)) )
#endif
            {   
				ArcAcPtr->offset = offset; ArcAcPtr->inten = inten;
                ArcAcPtr->dx=dx; ArcAcPtr->dy=dy; ArcAcPtr->dz=dz;
                ArcAcPtr++;
				idx_arc++;
            }

#ifndef PIXUTILS
            k1 = AbsFun(dy+iy);
            k2 = AbsFun(dy-iy);

            high = LookUp[rad][absx];
            if( ((k1>high)||(dz>= LookUp[LookUp[rad][k1]][absx]-zrate)) &&
                ((k2>high)||(dz>  LookUp[LookUp[rad][k2]][absx]+zrate)) )
#endif
            {   
				ArcDnPtr->offset = offset; ArcDnPtr->inten = inten;
                ArcDnPtr->dx=dx; ArcDnPtr->dy=dy; ArcDnPtr->dz=dz;
                ArcDnPtr++;
            }
        }
        temp += View.yskip;
    }
}


void Canvas3D::DrawCylinder( int x1, int y1, int z1, int x2, int y2, int z2, int c1, int c2, int rad )
{
    short  *dbase;
    ColorVal  *fbase;

    int zrate,zerr,ystep,err;
    int ix,iy,ax,ay;
    int lx,ly,lz;
    int mid,tmp;
    int temp;

    /* Trivial Case */
    if( (x1==x2) && (y1==y2) )
    {   
		if( z1>z2 )
        {      
			DrawSphere(x1,y1,z1,rad,c1);
        } 
		else 
			DrawSphere(x2,y2,z2,rad,c2);
        return;
    }

    if( z1<z2 )
    {   
		tmp=x1; x1=x2; x2=tmp;
        tmp=y1; y1=y2; y2=tmp;
        tmp=z1; z1=z2; z2=tmp;
        tmp=c1; c1=c2; c2=tmp;
    }

    DrawCylinderCaps(x1,y1,z1,x2,y2,z2,c1,c2,rad);

    lx = x2-x1;
    ly = y2-y1;
    lz = z2-z1;

    if( ly>0 ) { ystep = View.yskip; ay = ly; iy = 1; }
    else {   ystep = -View.yskip; ay = -ly; iy = -1; }
    if( lx>0 ) { ax = lx; ix = 1; }
    else { ax = -lx; ix = -1; }
    zrate = lz/MaxFun(ax,ay);

    temp = (int)y1*View.yskip+x1;
    fbase = View.fbuf+temp;
    dbase = View.dbuf+temp;

    if( ax>ay )
    {   
		lz -= ax*zrate;
        zerr = err = -(ax>>1);

        if( c1 != c2 )
        {   
			mid = (x1+x2)>>1;
            while( x1!=mid )
            {   
				z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ax; z1--; }
                fbase+=ix; dbase+=ix; x1+=ix;
                if( (err+=ay)>0 )
                {   
					fbase+=ystep; dbase+=ystep; err-=ax;
                       DrawArcDn(dbase,fbase,z1,c1);
                } 
				else 
					DrawArcAc(dbase,fbase,z1,c1);
            }
        }

        while( x1!=x2 )
        {   
			z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ax; z1--; }
            fbase+=ix; dbase+=ix; x1+=ix;
            if( (err+=ay)>0 )
            {   
				fbase+=ystep; dbase+=ystep; err-=ax;
                   DrawArcDn(dbase,fbase,z1,c2);
            } 
			else 
				DrawArcAc(dbase,fbase,z1,c2);
        }
    } 
	else /*ay>=ax*/
    {   
		lz -= ay*zrate;
        zerr = err = -(ay>>1);

        if( c1 != c2 )
        {   
			mid = (y1+y2)>>1;
            while( y1!=mid )
            {   
				z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ay; z1--; }
                fbase+=ystep; dbase+=ystep; y1+=iy;
                if( (err+=ax)>0 )
                {   
					fbase+=ix; dbase+=ix; err-=ay; 
                       DrawArcAc(dbase,fbase,z1,c1);
                } 
				else 
					DrawArcDn(dbase,fbase,z1,c1);
            }
        }

        while( y1!=y2 )
        {   
			z1 += zrate;  
			if( (zerr-=lz)>0 ) 
			{ 
				zerr-=ay; z1--; 
			}
            fbase+=ystep; dbase+=ystep; y1+=iy;
            if( (err+=ax)>0 )
            {   
				fbase+=ix; dbase+=ix; err-=ay; 
                   DrawArcAc(dbase,fbase,z1,c2);
            } 
			else 
				DrawArcDn(dbase,fbase,z1,c2);
        }
    }
}


int Canvas3D::TestCylinder( int x1, int y1, int z1, int x2, int y2, int z2, int rad )
{
    int tmp1, tmp2;

    if( m_UseSlabPlane )
        if( (z1+rad> m_SlabValue) || (z2+rad > m_SlabValue) )
            return(False);

    ClipStatus = False;

    tmp1 = x1+rad;  tmp2 = x2+rad;
    if( (tmp1<0) && (tmp2<0) )
        return( False );
    if( (tmp1>=View.xmax) || (tmp2>=View.xmax) )
        ClipStatus = True;

    tmp1 = x1-rad;  tmp2 = x2-rad;
    if( (tmp1>=View.xmax) && (tmp2>=View.xmax) )
        return( False );
    if( (tmp1<0) || (tmp2<0) )
        ClipStatus = True;

    tmp1 = y1+rad;  tmp2 = y2+rad;
    if( (tmp1<0) && (tmp2<0) )
        return( False );
    if( (tmp1>=View.ymax) || (tmp2>=View.ymax) )
        ClipStatus = True;

    tmp1 = y1-rad;  tmp2 = y2-rad;
    if( (tmp1>=View.ymax) && (tmp2>=View.ymax) )
        return( False );
    if( (tmp1<0) || (tmp2<0) )
        ClipStatus = True;

    return( True );
}



void Canvas3D::ClipArcAc( short* dbase, ColorVal* fbase, int x, int y, int z, int c)
{
    ArcEntry  *ptr;
    short  *dptr;
    short depth;
    int temp;

    ptr = ArcAc;
    while( (temp=y+ptr->dy) < 0 )
        if( ++ptr == ArcAcPtr )
            return;

    while( (temp<View.ymax) && (ptr<ArcAcPtr) )
    {   temp = x+ptr->dx;
        if( XValid(temp) )
        {   dptr = dbase+ptr->offset;  depth = ptr->dz+z;
            SETPIXEL(dptr,fbase+ptr->offset,depth,Lut[ptr->inten+c]);
        }
        ptr++;
        temp = y+ptr->dy;
    }
}

void Canvas3D::ClipArcDn(short* dbase, ColorVal* fbase, int x, int y, int z, int c)
{
    ArcEntry  *ptr;
    short  *dptr;
    short depth;
    int temp;

    ptr = ArcDn;
    while( (temp=y+ptr->dy) < 0 )
        if( ++ptr == ArcDnPtr )
            return;

    while( (temp<View.ymax) && (ptr<ArcDnPtr) )
    {   temp = x+ptr->dx;
        if( XValid(temp) )
        {   dptr = dbase+ptr->offset;  depth = ptr->dz+z;
            SETPIXEL(dptr,fbase+ptr->offset,depth,Lut[ptr->inten+c]);
        }
        ptr++;
        temp = y+ptr->dy;
    }
}


void Canvas3D::ClipCylinder( int x1, int y1, int z1, int x2, int y2, int z2, int c1, int c2, int rad )
{
    short  *dbase;
    ColorVal  *fbase;

    int zrate,zerr,ystep,err;
    int ix,iy,ax,ay;
    int lx,ly,lz;
    int mid,tmp;
    int temp;

    /* Visibility Tests */
    if( !TestCylinder(x1,y1,z1,x2,y2,z2,rad) ) return;

    if( !ClipStatus )
    {   
		DrawCylinder(x1,y1,z1,x2,y2,z2,c1,c2,rad);
        return;
    }
	
    /* Trivial Case */
    if( (x1==x2) && (y1==y2) )
    {   
		if( z1>z2 )
        {      
			ClipSphere(x1,y1,z1,rad,c1);
        } 
		else 
			ClipSphere(x2,y2,z2,rad,c2);
        return;
    }

    if( z1<z2 )
    {   tmp=x1; x1=x2; x2=tmp;
        tmp=y1; y1=y2; y2=tmp;
        tmp=z1; z1=z2; z2=tmp;
        tmp=c1; c1=c2; c2=tmp;
    }

    DrawCylinderCaps(x1,y1,z1,x2,y2,z2,c1,c2,rad);

    lx = x2-x1;
    ly = y2-y1;
    lz = z2-z1;

    if( ly>0 ) { ystep = View.yskip; ay = ly; iy = 1; }
    else {   ystep = -View.yskip; ay = -ly; iy = -1; }
    if( lx>0 ) { ax = lx; ix = 1; }
    else { ax = -lx; ix = -1; }
    zrate = lz/MaxFun(ax,ay);

    temp = (int)y1*View.yskip+x1;
    fbase = View.fbuf+temp;
    dbase = View.dbuf+temp;

    if( ax>ay )
    {   if( x2<x1 )
        {   tmp = c1;
            c1 = c2;
            c2 = tmp;
        }
        lz -= ax*zrate;
        zerr = err = -(ax>>1);
        mid = (x1+x2)/2;

        while( x1!=x2 )
        {   z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ax; z1--; }
            fbase+=ix; dbase+=ix; x1+=ix;
            if( (err+=ay)>0 )
            {   fbase += ystep;  err -= ax;
                dbase += ystep;  y1 += iy;
                   ClipArcDn(dbase,fbase,x1,y1,z1,(x1<mid?c1:c2));
            } else ClipArcAc(dbase,fbase,x1,y1,z1,(x1<mid?c1:c2));
        }
    } else /*ay>=ax*/
    {   if( y2<y1 )
        {   tmp = c1;
            c1 = c2;
            c2 = tmp;
        }
        lz -= ay*zrate;
        zerr = err = -(ay>>1);
        mid = (y1+y2)/2;

        while( y1!=y2 )
        {   z1 += zrate;  if( (zerr-=lz)>0 ) { zerr-=ay; z1--; }
            fbase+=ystep; dbase+=ystep; y1+=iy;
            if( (err+=ax)>0 )
            {   fbase += ix;  err -= ay;
                dbase += ix;  x1 += ix; 
                   ClipArcAc(dbase,fbase,x1,y1,z1,(y1<mid?c1:c2));
            } else ClipArcDn(dbase,fbase,x1,y1,z1,(y1<mid?c1:c2));
        }
    }
}


void Canvas3D::SetFontSize( int size )
{
    int count;
    int i;

    count = 0;
    for( i=0; i<23; i++ )
    {   
		FontDimen[i] = count>>4;
        count += size;
    }
    m_FontSize = size;
}


void Canvas3D::ClipCharacter( int x, int y, int z, int glyph, int col )
{
    register char *ptr;
    register int sx,sy;
    register int ex,ey;

    ptr = VectFont[glyph];
    while( *ptr )
    {   /* Uppercase test */
        if( ptr[0] < 'a' )
        {   
			sx = x + FontDimen[ptr[0]-'A'];
            sy = y + InvertY(FontDimen[ptr[1]-'a']);
            ptr += 2;
        } 
		else
        {   
			sx = ex;
            sy = ey;
        }

        ex = x + FontDimen[ptr[0]-'a'];
        ey = y + InvertY(FontDimen[ptr[1]-'a']);
        if( (ex!=sx) || (ey!=sy) )
        {   
			ClipLine(sx,sy,z,ex,ey,z,col);
        } 
		else 
			ClipPoint(ex,ey,z,col);
        ptr += 2;
    }
}


void 
Canvas3D::DisplayTextString( int x, int y, int z, const char* label, int col )
{
    int clip,high,max;
    char *ptr;
    int sx,sy;
    int ex,ey;
	
    high = (m_FontSize*3)>>1;
	
//    if( (y<0) || ((y-high)>=View.ymax) ) return;
//    clip = (y-high<0) || (y>=View.ymax);
    if( (y + high < 0) || ( y  >=View.ymax) ) return;
    clip = (y<0) || (y +high >=View.ymax);
	
    if( x < 0 )
    {   
		while( *label && (x<=-m_FontSize) )
        {   
			x += m_FontSize;  label++;
        }
		
        if( *label )
        {   
			ClipCharacter(x,y,z,(*label-32),col);
            x += m_FontSize;
            label++;
        } 
		else 
			return;
    }
	
    if( !clip )
    {   
		max = View.xmax-m_FontSize;
        while( *label && (x<max) )
        {  
			if(*label >= 32)
			{
				ptr = VectFont[*label-32];
			}
			else
			{
				ptr = VectFont[0];
			}
			while( *ptr )
			{   /* Uppercase test */
				if( ptr[0] < 'a' )
				{   
					sx = x + FontDimen[ptr[0]-'A'];
					sy = y + InvertY(FontDimen[ptr[1]-'a']);
					ptr += 2;
				} 
				else
				{   
					sx = ex;
					sy = ey;
				}
				
				ex = x + FontDimen[ptr[0]-'a'];
				ey = y + InvertY(FontDimen[ptr[1]-'a']);
				if( (ex!=sx) || (ey!=sy) )
				{   
					DrawTwinVector(sx,sy,z,ex,ey,z,col,col);
				} 
				else 
					PlotPoint(ex,ey,z,col);
				ptr += 2;
			}
			x += m_FontSize;
			label++;
        }
		
        if( *label )
		{
			if(*label >= 32)
				ClipCharacter(x,y,z,(*label-32),col);
			else
				ClipCharacter(x,y,z, 0 ,col);
		}
    } 
	else /* Always Clip! */
	{
        while( *label && (x<View.xmax) )
        {  
			if(*label >= 32)
				ClipCharacter(x,y,z,(*label-32),col);
			else
				ClipCharacter(x,y,z, 0 ,col);
            x += m_FontSize;
            label++;
        }
	}
}

