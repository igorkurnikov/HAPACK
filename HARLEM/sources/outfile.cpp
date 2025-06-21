/*! \file outfile.cpp

   Functions to save a picture into a graphics  file

   \author  Igor Kurnikov
   \date 1998-2002

   derived from:
 
   outfile.c
   RasMol2 Molecular Graphics
   Roger Sayle, August 1995
   Version 2.6

*/
#define OUTFILE

#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#if(wxUSE_GUI != 0)
#include "wx/image.h"
#endif
#include "wx/wfstream.h"

#include "vec3d.h"

#include "haconst.h"

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "command.h"
#include "abstree.h"
#include "hamolview.h"
#include "hasurface.h"

#include "haatom.h"
#include "habond.h"
#include "hamolecule.h"
#include "moleditor.h"

/* Sun rasterfile.h macro defns */
#define RAS_MAGIC       0x59a66a95
#define RAS_RLE         0x80
#define RT_STANDARD     1
#define RT_BYTE_ENCODED 2
#define RMT_NONE        0
#define RMT_EQUAL_RGB   1

/* Standard A4 size page: 8.267x11.811 inches */
/* U.S. Normal size page: 8.500x11.000 inches */
#define PAGEHIGH  (11.811*72.0)
#define PAGEWIDE  (8.267*72.0)
#define BORDER    0.90

/* Compression Ratio   0<x<127 */
#define EPSFCompRatio  32

#define Round(x)       ((int)(x))


#define PSBond      0x00
#define PSHBond     0x01
#define PSSSBond    0x02
#define PSAtom      0x03
#define PSRibbon    0x04
#define PSMonit     0x05


static short ABranch[4096];
static short DBranch[4096];
static short Hash[LutSize];

/* Apple PICT macros */
#define PICTcliprgn         0x0001
#define PICTpicversion      0x0011
#define PICTpackbitsrect    0x0098
#define PICTdirectbitsrect  0x009a
#define PICTendofpict       0x00ff
#define PICTheaderop        0x0c00  
    
static int GIFClrCode; 
static int GIFEOFCode;

static short RLELineSize;
static unsigned int RLEFileSize;
static int RLEEncode;
static int RLEOutput;
static int RLELength;
static int RLEPixel;
static int RLEChar;

static unsigned char Buffer[LutSize];
static int LineLength;
static FILE *OutFile;
static unsigned int BitBuffer;
static int BitBufLen;
static int PacketLen;
static int CodeSize;

static double LineWidth;
static int VectSolid;
static int VectCol;


void 
HaMolView::WriteByte(int  val )
{
    putc( val, OutFile );
}

void 
HaMolView::WriteLSBShort(int  val )
{
    putc( val&0xff, OutFile );
    putc( (val>>8)&0xff, OutFile );
}

void 
HaMolView::WriteMSBShort(int  val )
{
    putc( (val>>8)&0xff, OutFile );
    putc( val&0xff, OutFile );
}

void 
HaMolView::WriteMSBLong(unsigned int val )
{
    putc((int)((val>>24)&0xff),OutFile);
    putc((int)((val>>16)&0xff),OutFile);
    putc((int)((val>>8) &0xff),OutFile);
    putc((int)(val&0xff),OutFile);
}

void 
HaMolView::WritePPMWord(int i )
{
    if( i>99 )
    {   
		putc((i/100)+'0',OutFile); i %= 100;
        putc((i/10) +'0',OutFile); i %= 10;
    } 
	else if( i>9 )
    {   
		putc((i/10)+'0',OutFile);  i %= 10;
    }
    putc(i+'0',OutFile);
}


int 
HaMolView::WritePPMFile( const char* name,int raw )
{
    ColorVal  *ptr;
    int i,col;
    int x,y;

#if defined(_WIN32) 
    OutFile = fopen(name, (raw?"wb":"w") );
#else
    OutFile = fopen(name,"w");
#endif

    if( !OutFile ) 
    {   
		PrintLog("\n Output Error: Unable to create file %s!\n",name);
        return( False );
    }

    fprintf(OutFile,"P%c %d %d 255\n",(raw?'6':'3'), pCanv->XRange(), pCanv->YRange());

    ptr = pCanv->FBuffer;

    if( !raw )
    {   col = 0;
        for( y=pCanv->YRange()-1; y>=0; y-- )
        {
            for( x=0; x < pCanv->XRange(); x++ )
            {   
				i = *ptr++;
                WritePPMWord((int)RComp(i));  WriteByte(' ');
                WritePPMWord((int)GComp(i));  WriteByte(' ');
                WritePPMWord((int)BComp(i));  
                if( ++col == 5 )
                {   
					WriteByte('\n');
                    col = 0;
                } 
				else 
					WriteByte(' ');
            }
        }
    } 
	else
        for( y= pCanv->YRange()-1; y>=0; y-- )
        {
            for( x=0; x < pCanv->XRange(); x++ )
            {   i = *ptr++;
                putc((int)RComp(i),OutFile);
                putc((int)GComp(i),OutFile);
                putc((int)BComp(i),OutFile);
            }
        }

    fclose(OutFile);
    return( True );
}

void 
HaMolView::WriteGIFCode(int code )
{
    register int max;

    max = (code==GIFEOFCode)? 0 : 7;
    BitBuffer |= ((unsigned int)code<<BitBufLen);
    BitBufLen += CodeSize;

    while( BitBufLen > max )
    {    
#ifdef _WIN32
         Buffer[PacketLen++]=(unsigned char)(BitBuffer & 0xff);
#else
         Buffer[PacketLen++]=BitBuffer;
#endif
         BitBuffer >>= 8;
         BitBufLen -= 8;

        if( PacketLen==255 )
        {   
			WriteByte(0xff);
            fwrite((char*)Buffer,1,255,OutFile);
            PacketLen = 0;
        }
    }
}

int 
HaMolView::WriteGIFFile( const char* name )
{
    int i,j,cols;
    int pref,next,last;
    int isize, ilast;
    ColorVal  *ptr;
    short  *prev;
    int x,y,init;

    cols = ComputeRevColourMap();
    if( cols<2 ) return( False );

	if(cols > 256)
	{
		ErrorInMod("HaMolView::WriteGIFFile",
			       " Number of colors in minimal color map larger than 256, image will be broken ");	
	}

    for( isize=2; isize<8; isize++ ) // find minimal n 2^n > ncol in minimal color map 
	{
        if( (1<<isize) >= cols )
		{
			break;
		}
	}
	
	int old_num_cols = cols;

    cols = 1<<isize;
	min_color_map.resize(cols,0);


#if defined(_WIN32) 
    OutFile = fopen(name,"wb");
#else
    OutFile = fopen(name,"w");
#endif
    if( !OutFile ) 
    {    
		 PrintLog("\n Output Error: Unable to create file %s!\n",name);
         return( False );
    }
    fputs("GIF87a",OutFile);
    WriteLSBShort(pCanv->XRange());
    WriteLSBShort(pCanv->YRange());
    WriteByte(0xf0|(isize-1)); 
    WriteByte(0x00); 
    WriteByte(0x00);

    for( j=0; j < cols; j++ )
    {   
		i = min_color_map[j];
        WriteByte(RComp(i));
        WriteByte(GComp(i));
        WriteByte(BComp(i));
    }

    WriteByte(',');
    WriteByte(0x00);  WriteByte(0x00);
    WriteByte(0x00);  WriteByte(0x00);
    WriteLSBShort(pCanv->XRange());
    WriteLSBShort(pCanv->YRange());
    WriteByte(0x00);  WriteByte(isize);

    PacketLen=0;
    BitBuffer=0;
    BitBufLen=0;

    GIFClrCode = (1<<isize);
    GIFEOFCode = GIFClrCode+1;
    ilast = (GIFClrCode<<1)-GIFEOFCode;
    isize++;

    CodeSize = isize;
    last = ilast;
    next = 1;  
   
    WriteGIFCode(GIFClrCode);
    for( i=0; i<cols; i++ )
        Hash[i]=0;

	unsigned char Node[4096];


    /* Avoid Warnings! */
    prev = (short *)0; 
    pref = 0;

// LutInv[x] substitute rev_color_map[x]
// Node[x]   substitute min_color_map[x]

    init = False;
    for( y=pCanv->YRange()-1; y>=0; y-- )
//    for( y= 0; y < pCanv->YRange(); y++ )
    {   
		ptr = pCanv->FBuffer + y * pCanv->XRange();
        for( x=0; x< pCanv->XRange(); x++ )
        {   
			if( !init )
            {   
				pref = rev_color_map[*ptr++];
                prev = Hash+pref;
                init = True;
                continue;
            }

			i = rev_color_map[*ptr++];

            while( *prev && (Node[*prev] != i) )
                prev = ABranch+*prev;

            if( *prev )
            {   
				pref = *prev+GIFEOFCode;
                prev = DBranch+*prev;
            } 
			else
            {   
				WriteGIFCode(pref);
                if( next==last )
                {   
					if( CodeSize==12 )
                    {   
						WriteGIFCode(GIFClrCode);
                        pref = i;  prev = Hash+i;
                        for( i=0; i<cols; i++ )
                            Hash[i] = 0;
                        CodeSize = isize;
                        last = ilast;
                        next = 1; 
                        continue;
                    }
                    last = (last<<1)+GIFEOFCode;
                    CodeSize++;
                }
                *prev = next;
                ABranch[next] = 0;
                DBranch[next] = 0;
                Node[next] = i;
//				min_color_map[next] = i;
                prev = Hash+i;
                pref = i;
                next++;
            }
        }
    }

    WriteGIFCode(pref);
    WriteGIFCode(GIFEOFCode);
    if( PacketLen )
    {   
		WriteByte(PacketLen);
        fwrite((char*)Buffer,1,PacketLen,OutFile);
    }

    WriteByte(0x00);
    WriteByte(';');
    fclose(OutFile);

    return( True );

}

int 
HaMolView::WriteBMPFile( const char* name )
{
#if !defined(HA_NOGUI)
//	::wxInitAllImageHandlers();
	unsigned int x_src = pCanv->XRange();
	unsigned int y_src = pCanv->YRange();
	wxImage img(x_src,y_src);
	int ires = SetWXImage(img);
    if(!ires) return FALSE;
	
    wxFFileOutputStream fstream(name);
	img.SaveFile(fstream, wxBITMAP_TYPE_BMP);
	return ires;
#else
	return FALSE;
#endif
}

int HaMolView::WriteJPEGFile( const char* name )
{
#if !defined(HA_NOGUI)
	::wxInitAllImageHandlers();
	unsigned int x_src = pCanv->XRange();
	unsigned int y_src = pCanv->YRange();
	wxImage img(x_src,y_src);
	int ires = SetWXImage(img);
    if(!ires) return FALSE;
	
    wxFFileOutputStream fstream(name);
	img.SaveFile(fstream, wxBITMAP_TYPE_JPEG);
	return ires;
#else
	return FALSE;
#endif
}

int 
HaMolView::WriteTIFFFile( const char* name )
{
#if !defined(HA_NOGUI)
	::wxInitAllImageHandlers();
	unsigned int x_src = pCanv->XRange();
	unsigned int y_src = pCanv->YRange();
	wxImage img(x_src,y_src);
	int ires = SetWXImage(img);
    if(!ires) return FALSE;
	
    wxFFileOutputStream fstream(name);
	img.SaveFile(fstream, wxBITMAP_TYPE_TIF );
	return ires;
#else
	return FALSE;
#endif
}

int 
HaMolView::WritePNGFile( const char* name )
{
#if !defined(HA_NOGUI)
	::wxInitAllImageHandlers();
	unsigned int x_src = pCanv->XRange();
	unsigned int y_src = pCanv->YRange();
	wxImage img(x_src,y_src);
	int ires = SetWXImage(img);
    if(!ires) return FALSE;
	
    wxFFileOutputStream fstream(name);
	img.SaveFile(fstream, wxBITMAP_TYPE_PNG);
	return ires;
#else
	return FALSE;
#endif
}

int 
HaMolView::WritePCXFile( const char* name )
{
#if !defined(HA_NOGUI)
	::wxInitAllImageHandlers();
	unsigned int x_src = pCanv->XRange();
	unsigned int y_src = pCanv->YRange();
	wxImage img(x_src,y_src);
	int ires = SetWXImage(img);
    if(!ires) return FALSE;
	
    wxFFileOutputStream fstream(name);
	img.SaveFile(fstream, wxBITMAP_TYPE_PCX);
	return ires;
#else
	return FALSE;
#endif
}


int 
HaMolView::FindDepth( PSItemPtr item, int type )
{
    Monitor  *monit;
    HaHBond  *hbond;
    HaAtom   *atom;
    HaBond   *bond;
    int result;

    switch( type )
    {   
	    case(PSAtom):    atom = (HaAtom *)item;
                         return( atom->z );

        case(PSBond):    bond = (HaBond *)item;
                         result = bond->srcatom->z;
                         if( result < bond->dstatom->z )
                             result = bond->dstatom->z;
                         return( result );

        case(PSSSBond):  
        case(PSHBond):   hbond = (HaHBond *)item;
                         if( (type==PSHBond)? HBondMode : SSBondMode )
                         {   
							 result = hbond->srcCA->z;
                             if( result < hbond->dstCA->z )
                                 result = hbond->dstCA->z;
                         } 
						 else
                         {   
							 result = hbond->src->z;
                             if( result < hbond->dst->z )
                                 result = hbond->dst->z;
                         }
                         return( result );

        case(PSMonit):   monit = (Monitor *)item;
                         result = monit->src->z;
                         if( result < monit->dst->z )
                             result = monit->dst->z;
                         return( result );
    }
    return( 0 );
}


void 
HaMolView::DepthSort(PSItemPtr* data, char* type, int count )
{
    register char ttmp;
    register void  *dtmp;
    register int i, j, k;
    register int depth;
    register int temp;

    for( i=1; i<count; i++ )
    {   
		dtmp = data[i];  
        ttmp = type[i];

        j = i-1;
        depth = FindDepth(dtmp,ttmp);
        temp = FindDepth(data[j],type[j]);
        while( (depth<temp) || ((depth==temp)&&(ttmp<type[j])) )
            if( j-- ) 
            {   
				temp = FindDepth(data[j],type[j]);
            } 
			else 
				break;
        j++;

        if( j != i )
        {   
			for( k=i; k>j; k-- )
            {    
				data[k] = data[k-1];
                 type[k] = type[k-1];
            }
            data[j] = dtmp;
            type[j] = ttmp;
        }
    }
}


int 
HaMolView::ClipVectSphere(HaAtom* ptr )
{
    register int rad;

    rad = ptr->irad;

    if( ptr->x + rad < 0 )  return( True );
    if( ptr->y + rad < 0 )  return( True );
    if( ptr->x - rad >= pCanv->XRange() )  return( True );
    if( ptr->y - rad >= pCanv->YRange() )  return( True );
    return( False );
}


int 
HaMolView::ClipVectBond( HaAtom* src, HaAtom* dst )
{
    if( !src || !dst )  return( True );
    if( (src->x<0) && (dst->x<0) )  return( True );
    if( (src->y<0) && (dst->y<0) )  return( True );
    if( (src->x>= pCanv->XRange()) && (dst->x>= pCanv->XRange()) )  return( True );
    if( (src->y>= pCanv->YRange()) && (dst->y>= pCanv->YRange()) )  return( True );
    return( False );
}



void 
HaMolView::WriteVectColour( int col )
{
    if( col != VectCol )
    {   
		fprintf(OutFile,"%g ",(double)RComp(Canvas3D::Lut[col])/255.0);
        fprintf(OutFile,"%g ",(double)GComp(Canvas3D::Lut[col])/255.0);
        fprintf(OutFile,"%g ",(double)BComp(Canvas3D::Lut[col])/255.0);
        fputs("setrgbcolor\n",OutFile);
        VectCol = col;
    }
}


#define MAXSECT 5
typedef struct {
        /* Ellipse */
        double ephi,epsi;
        double etheta;
        double ex,ey;
        double erad;

        /* Sphere */
        double sphi,spsi;
        int sx,sy;
        double srad;
    } SphereSect;


static int VectClipContain( SphereSect *x, SphereSect *y )
{
    if( x->erad != 0.0 )
    {   
		if( y->erad != 0.0 )
            /* Simple segment containment test! */
            return( ((x->sphi+x->spsi)>=(y->sphi+y->spsi)) &&
                    ((x->sphi-x->spsi)<=(y->sphi-y->spsi)) );
    } 
	else if( y->erad == 0.0 )
        return( x->srad >= y->srad );
    return( False );
}


void 
HaMolView::WriteVectSphere( PSItemPtr *data, char *type,int index)
{
    register int ecount, count;
    register HaAtom  *atm;
    register HaAtom  *ptr;
    register int dist2,dist3;
    register int dx, dy, dz;
    register int i,j,k;

    register double b,d,f,g,x;
    register double radf,radb;
    register double phi1,phi2;
    register double temp,psi;
    register double theta;

    register SphereSect *sptr;
    SphereSect sect[MAXSECT];

    ptr = (HaAtom *)data[index];
    radf = ptr->image_radius*Scale;

    count = 0;
    ecount = 0;
    sptr = sect;
    for( i=index-1; i>=0; i-- )
    {   
		if( type[i] != PSAtom )
            continue;

        atm = (HaAtom *)data[i];
        /* Atom can't intersect visibly! */
        if( atm->z + atm->irad < ptr->z )
            continue;

        dx = atm->x - ptr->x; 
        dy = atm->y - ptr->y; 
        dz = atm->z - ptr->z;

        dist2 = (int)dx*dx + (int)dy*dy;
        dist3 = dist2 + dz*dz;

        radb = atm->image_radius*Scale;  
        temp = radf + radb;

        /* Atoms don't intersect! */
        if( dist3 > temp*temp ) continue;


        d = sqrt( (double)dist3 );
        f = (temp*(radf-radb)+dist3)/(2.0*d);
        theta = -dz/d;

        if( f>0 )
        {   
			temp = radf*radf;
            /* Intersection not visible! */
            if( theta*temp > temp-f*f )
                continue;
        } 
		else if( f < -radf )
            return;

        x = sqrt( (radf-f)*(radf+f) );

        if( dx || dy )
        {   
			g = sqrt( (double)dist2 );
            psi = RAD_TO_DEG*atan2((double)dy,(double)dx);
            b = (f*(dz*dz))/(d*g);

            if( AbsFun(b)>x )
                continue;

            phi1 = b + (f*g)/d;
            phi1 = RAD_TO_DEG*acos(phi1/radf);
            if( phi1!=phi1 ) continue;

            phi2 = (d*b)/g;
            if( AbsFun(phi2) < x )
            {   phi2 = RAD_TO_DEG*acos(phi2/x);
                if( phi2!=phi2 ) continue;
                if( phi2 > 90.0 ) 
                    phi2 = 180.0-phi1;
            } 
			else 
				phi2 = 90.0;

            sptr->erad = x;
            sptr->etheta = -theta;
            sptr->ephi = psi;
            sptr->epsi = phi2;

            temp = f/d;
            sptr->ex = ptr->x+temp*dx;
            sptr->ey = ptr->y+temp*dy;

            sptr->srad = radf;
            sptr->sphi = psi;
            sptr->spsi = phi1;
            sptr->sx = ptr->x;
            sptr->sy = ptr->y;

        } 
		else
        {   
			x = sqrt( (radf-g)*(radf+g) );

            sptr->srad = x;
            sptr->erad = 0.0;
            sptr->sx = ptr->x;
            sptr->sy = ptr->y;
            sptr->sphi = 180;
            sptr->spsi = -180;
        }

        /* Optimize Segments */
        j = 0;
        while( j<count )
            if( VectClipContain(sptr,sect+j) )
            {   /* Delete Segment sect[j] */
                for( k=j; k<count; k++ )
                    sect[k] = sect[k+1];
                count--;  sptr--;
            } 
			else if( VectClipContain(sect+j,sptr) )
            {   
				break;  /* Exclude Segment */
            } 
			else 
				j++;
           

        if( j==count )
        {   
			count++;  sptr++;
            if( sptr->erad != 0.0 )
                ecount++;
            if( count==MAXSECT )
                break;
        }
    }

    if( UseOutLine )
    {   
		temp = (ptr->z-ZOffset())/GetImageSize() + 1.0;
        if( temp != LineWidth )
        {   
			fprintf(OutFile,"%g setlinewidth\n",temp);
            LineWidth = temp;
        }
    }

    if( !VectSolid )
    {   
		fputs("[] 0 setdash\n",OutFile);
        VectSolid = True;
    }

    if( count )
    {   
		fputs("gsave\n",OutFile);
        fprintf(OutFile,"%%%% %d %d\n",count,ecount);

        sptr = sect;
        for( i=0; i<count; i++ )
        {   
			if( sptr->erad != 0.0 )
            {   
				fprintf(OutFile,"%g %g %g %g %g %g ClipEllips\n",
                            sptr->erad,sptr->epsi,sptr->etheta,
                            sptr->ephi,sptr->ex,sptr->ey);
            }

            if( (i==count-1) || (sptr->erad==0.0) )
            {   
				fprintf(OutFile,"%g %g %g %d %d ClipSphere\n",sptr->srad,
                                sptr->sphi+sptr->spsi,sptr->sphi-sptr->spsi,
                                sptr->sx, sptr->sy );
            } 
			else 
				fprintf(OutFile,"%g %g %g %d %d ClipBox\n",
                                    sptr->srad+sptr->srad+2,
                                    sptr->srad+1, sptr->ephi,
                                    sptr->sx, sptr->sy );
            sptr++;
        }

        i = ptr->col + ColourMask;
        fprintf(OutFile,"%g ",(double)RComp(Canvas3D::Lut[i])/255.0);
        fprintf(OutFile,"%g ",(double)GComp(Canvas3D::Lut[i])/255.0);
        fprintf(OutFile,"%g ",(double)BComp(Canvas3D::Lut[i])/255.0);
        fprintf(OutFile,"%g Shade\n",radf);
        fputs("grestore\n\n",OutFile);
    } 
	else
    {   
		i = ptr->col + ColourMask;
        fprintf(OutFile,"%g ",(double)RComp(Canvas3D::Lut[i])/255.0);
        fprintf(OutFile,"%g ",(double)GComp(Canvas3D::Lut[i])/255.0);
        fprintf(OutFile,"%g ",(double)BComp(Canvas3D::Lut[i])/255.0);
        fprintf(OutFile,"%g %d %d ",radf,ptr->x,ptr->y);
        fputs("Sphere\n\n",OutFile);
    }
}


void 
HaMolView::WriteVectWire( HaAtom* src, HaAtom* dst, int col, int dash )
{
    HaAtom  *tmp;
    double radius;
    double temp;
    double dist;

    double midx, midy;
    double endx, endy;
    int col1, col2;
    int dx, dy, dz;
    int dist2;
    int inten;


    if( src->z > dst->z )
    {   
		tmp = src;
        src = dst;
        dst = tmp;
    }

    if( !col )
    {   
		col1 = src->col;
        col2 = dst->col;
    } 
	else 
		col1 = col2 = col;

    if( UseBackFade )
    {   
		dz = (src->z+dst->z)>>1;
        inten = (ColourDepth*(dz+GetImageRadius()-ZOffset()))/GetImageSize();
    } 
	else 
		inten = ColourMask;

    dx = dst->x - src->x;  
    dy = dst->y - src->y;
    dist2 = dx*dx + dy*dy;
    dist = sqrt( (double)dist2 );

    if( dst->IsDrawSphere() )
    {   
		radius = dst->image_radius*Scale;
        if( dist <= radius ) return;

        /* Test for second half obscured! */
        if( (col1!=col2) && (0.5*dist < radius) )
            col2 = col1;
    }

    if( src->IsDrawSphere() )
    {   
		radius = src->image_radius*Scale;
        if( dist <= radius ) return;

        /* Test for first half obscured! */
        if( (col1!=col2) && (0.5*dist < radius) )
            col1 = col2;
    }

    WriteVectColour( col1+inten );

    dz = (src->z+dst->z)>>1;
    temp = (double)(dz-ZOffset())/GetImageSize() + 1.0;
    if( temp != LineWidth )
    {   
		fprintf(OutFile,"%g setlinewidth\n",temp);
        LineWidth = temp;
    }

    if( dash )
    {   
		if( VectSolid )
        {   
			fputs("[3 3] 0 setdash\n",OutFile);
            VectSolid = False;
        }
    } 
	else
        if( !VectSolid )
        {   
			fputs("[] 0 setdash\n",OutFile);
            VectSolid = True;
        }


    if( src->IsDrawSphere() )
    {   
		dz = dst->z - src->z;
        dist = sqrt( (double)(dist2 + dz*dz) );
        endx = src->x + (radius*dx)/dist;
        endy = src->y + (radius*dy)/dist;
        fprintf(OutFile,"%g %g ",endx,endy);
    } 
	else
        fprintf(OutFile,"%d %d ",src->x,src->y);

    if( col1 != col2 )
    {   
		midx = 0.5*(src->x + dst->x);
        midy = 0.5*(src->y + dst->y);
        fprintf(OutFile,"%g %g Wire\n",midx,midy);

        WriteVectColour( col2+inten );
        fprintf(OutFile,"%g %g ",midx,midy);
    } 
    fprintf(OutFile,"%d %d Wire\n",dst->x,dst->y);
}


void 
HaMolView::WriteVectStick(HaAtom* src,HaAtom* dst, int col, int rad )
{
    register HaAtom  *tmp;
    register double midx, midy;
    register double relx, rely;
    register double endx, endy;
    register double radius, angle;
    register double dist, dist3;
    register double temp, ratio;
	
    register int dist2;
    register int dx, dy, dz;
    register int col1, col2;
    register int i, inten;
	
    if( !rad )
    {   
		WriteVectWire(src,dst,col,False);
        return;
    }
	
    if( src->z > dst->z )
    {   
		tmp = src;
        src = dst;
        dst = tmp;
    }
	
    if( !col )
    {   
		col1 = src->col;
        col2 = dst->col;
    } 
	else 
		col1 = col2 = col;
	
    dx = dst->x - src->x;  
    dy = dst->y - src->y;
    dz = dst->z - src->z;
    dist2 = dx*dx + dy*dy;
    dist3 = sqrt( (double)(dist2 + dz*dz) );
    dist = sqrt( (double)dist2 );
	
    if( dst->IsDrawSphere() )
    {   
		radius = dst->image_radius*Scale;
        if( dist <= radius ) return;
		
        /* Test for nearest half obscured! */
        if( (col1!=col2) && (0.5*dist < radius) )
            col2 = col1;
    }
	
    if( src->IsDrawSphere() )
    {   
		radius = src->image_radius*Scale;
        if( dist <= radius ) return;
		
        /* Test for furthest half obscured! */
        if( (col1!=col2) && (0.5*dist < radius) )
            col1 = col2;
    }
	
    if( !VectSolid )
    {   
		fputs("[] 0 setdash\n",OutFile);
        VectSolid = True;
    }
	
    temp = ((src->z-ZOffset())+(dst->z-ZOffset()))/GetImageSize() + 1.0;
    if( temp != LineWidth )
    {   
		fprintf(OutFile,"%g setlinewidth\n",temp);
        LineWidth = temp;
    }
	
    radius = rad*Scale;
    angle = RAD_TO_DEG*atan2((double)dy,(double)dx);
    inten = (int)((dist/dist3)*ColourMask);
	
    if( col1 != col2 )
    {   
		midx = 0.5*(src->x + dst->x);
        midy = 0.5*(src->y + dst->y);
        relx = (radius*dx)/dist;
        rely = (radius*dy)/dist;
		
        fprintf(OutFile,"%g %g moveto\n",midx+rely,midy-relx);
        fprintf(OutFile,"%g %g lineto\n",midx-rely,midy+relx);
		
        ratio = dz/dist3;
		
        if( (src->IsDrawSphere()) && (src->image_radius>rad) )
        {   
			temp = (Scale*src->image_radius)/dist3;
            endx = src->x + temp*dx;
            endy = src->y + temp*dy;
			
            fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
            fprintf(OutFile,"%g %g StickEnd\n",endx,endy);
        } 
		else
        {   
			fprintf(OutFile,"%d %d %g ",src->x,src->y,radius);
            fprintf(OutFile,"%g %g arc\n",angle+90,angle-90);
        }
        fputs("closepath ",OutFile);
		
        i = col1 + inten;
        fprintf(OutFile,"%g ",(double)RComp(Canvas3D::Lut[i])/255.0);
        fprintf(OutFile,"%g ",(double)GComp(Canvas3D::Lut[i])/255.0);
        fprintf(OutFile,"%g ",(double)BComp(Canvas3D::Lut[i])/255.0);
        fputs("setrgbcolor fill\n",OutFile);
		
        fprintf(OutFile,"%d %d %g ",dst->x,dst->y,radius);
        fprintf(OutFile,"%g %g arc\n",angle-90,angle+90);
        fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
        fprintf(OutFile,"%g %g StickEnd\n",midx,midy);
        fputs("closepath ",OutFile);
		
        i = col2 + inten;
        fprintf(OutFile,"%g ",(double)RComp(Canvas3D::Lut[i])/255.0);
        fprintf(OutFile,"%g ",(double)GComp(Canvas3D::Lut[i])/255.0);
        fprintf(OutFile,"%g ",(double)BComp(Canvas3D::Lut[i])/255.0);
        fputs("setrgbcolor fill\n",OutFile);
		
        if( UseOutLine )
        {   
			fprintf(OutFile,"%d %d %g ",dst->x,dst->y,radius);
            fprintf(OutFile,"%g %g arc\n",angle-90,angle+90);
            if( (src->IsDrawSphere()) && (src->image_radius>rad) )
            {   
				fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
                fprintf(OutFile,"%g %g StickEnd\n",endx,endy);
            } 
			else
            {   
				fprintf(OutFile,"%d %d %g ",src->x,src->y,radius);
                fprintf(OutFile,"%g %g arc\n",angle+90,angle-90);
            }
            fputs("closepath 0 setgray stroke\n",OutFile);
        }
    } 
	else /* col1 == col2! */
    {   
		fprintf(OutFile,"%d %d %g ",dst->x,dst->y,radius);
        fprintf(OutFile,"%g %g arc\n",angle-90,angle+90);
		
        if( (src->IsDrawSphere() ) && (src->image_radius>rad) )
        {   
			temp = (Scale*src->image_radius)/dist3;
            endx = src->x + temp*dx;
            endy = src->y + temp*dy;
            ratio = dz/dist3;
			
            fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
            fprintf(OutFile,"%g %g StickEnd\n",endx,endy);
        } 
		else
        {   
			fprintf(OutFile,"%d %d %g ",src->x,src->y,radius);
            fprintf(OutFile,"%g %g arc\n",angle+90,angle-90);
        }
		
        i = col1 + inten;
        fprintf(OutFile,"%g ",(double)RComp(Canvas3D::Lut[i])/255.0);
        fprintf(OutFile,"%g ",(double)GComp(Canvas3D::Lut[i])/255.0);
        fprintf(OutFile,"%g ",(double)BComp(Canvas3D::Lut[i])/255.0);
        fputs("Stick\n",OutFile);
    }
    VectCol = 0;
}


void 
HaMolView::WriteVectDots()
{
    double xi,yi;
    int inten;
    int temp;
    int zi;
    int i;

    if( LineWidth != 1.0 )
    {   
		fputs("1 setlinewidth\n",OutFile);
        LineWidth = 1.0;
    }
    temp = SlabValue() - ZOffset();

	MolSet* pmset = GetMolSet();
	list<Object3D*>::iterator oitr;

	int ixadd=  pCanv->XRange()/2;
	int iyadd=  pCanv->YRange()/2;
	int izadd=  ZOffset();

	for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end(); oitr++)
	{
		if( (*oitr)->GetObjType() != OBJ3D_DOT_SURFACE )
			continue;

		DotStruct* ptr = (DotStruct*)(*oitr);

		int np = ptr->GetCount();
        for( i=0; i< np; i++ )
        {   
            xi = ptr->dots[i].GetX()* Scale + ixadd;
            if( (xi<0.0) || (xi>=pCanv->XRange()) ) continue;
            yi = ptr->dots[i].GetY()* Scale + iyadd;
            if( (yi<0.0) || (yi>=pCanv->YRange()) ) continue;
            zi = (int)(ptr->dots[i].GetZ()* Scale) + izadd;
            if( UseSlabPlane() && (zi>=temp) ) continue;

            inten = (ColourDepth*(zi+ GetImageRadius()))/GetImageSize();
            WriteVectColour( ptr->dots[i].col + inten );
            fprintf(OutFile,"%g %g Dot\n",xi,yi);
        }
	}
}

void 
HaMolView::WriteVectLabels()
{
    HaChain  *chain;
    HaResidue  *group;
    HaAtom  *aptr;
    char buffer[80];

    fputs("/Times-Roman",OutFile); /* Courier or Courier-Bold? */
    fprintf(OutFile," findfont %d scalefont setfont\n",pCanv->m_FontSize<<1);

    if( UseLabelCol )
    {   
		if( BackColor.r || BackColor.g || BackColor.b )
        {   
			fprintf(OutFile,"%g %g %g setrgbcolor\n",
                    LabelColor.r/250.0, LabelColor.g/250.0, LabelColor.b/250.0);
        } 
		else 
			fputs("0 setgray\n",OutFile);
    } 
	else 
		VectCol = 0;

	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
        if( !aptr->label.empty() )
        {   
			if( !UseLabelCol && (aptr->col!=VectCol) )
				WriteVectColour( aptr->col );
			
			chain = aptr->GetHostChain();
			group = aptr->GetHostRes();
            FormatLabel(chain,group,aptr,aptr->label.c_str(),buffer);
            fprintf(OutFile,"(%s) %d %d Label\n",buffer,aptr->x,aptr->y);
        }
	}
}


void 
HaMolView::WriteVectMonitors()
{
    register HaAtom  *s;
    register HaAtom  *d;
    register int x,y,col;

    register char *cptr;
    register int dist;
    char buffer[10];
 
    buffer[9] = '\0';
    buffer[6] = '.';

    fputs("/Times-Roman",OutFile); /* Courier or Courier-Bold? */
    fprintf(OutFile," findfont %d scalefont setfont\n",pCanv->m_FontSize<<1);


	list<Monitor>::iterator mtr;
    for( mtr=MonitList.begin(); mtr != MonitList.end(); mtr++ )
    {   
		s = (*mtr).src;
        d = (*mtr).dst;

        if( ZValid_v( (s->z+d->z)/2 ) )
        {   x = (s->x+d->x)/2;
            y = (s->y+d->y)/2;
 
            if( !UseLabelCol )
            {   /* Use Source atom colour! */
                if( (*mtr).col )
                {   
					col = (*mtr).col + (ColourMask>>1);
                } 
				else 
					col = s->col + (ColourMask>>1);
            } 
			else 
				col = LabelColor.cidx;
            WriteVectColour(col);
 
            dist = (int)((*mtr).dist);
            buffer[8] = (dist%10)+'0';  dist /= 10;
            buffer[7] = (dist%10)+'0';
            cptr = &buffer[5];
 
            if( dist > 9 )
            {   
				do {
                    dist /= 10;
                    *cptr-- = (dist%10)+'0';
                } 
				while( dist > 9 );
                cptr++;
            } 
			else 
				*cptr = '0';
 
            fprintf(OutFile,"(%s) %d %d Label\n",cptr,x+4,y);
        }
    }
}


int 
HaMolView::CountPSItems()
{
    HaChain  *chain;
    set<HaHBond, less<HaHBond> >::iterator hptr;
    HaBond  *bptr;
    HaAtom  *aptr;
    int result;

    result = 0;

	MoleculesType::iterator mol_itr;
	MolSet* pmset = GetMolSet();

    if( DrawAtoms )
	{
		AtomIteratorMolSet aitr(this->GetMolSet());
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
            if( aptr->IsDrawSphere() ) 
			{
                if( !UseClipping || !ClipVectSphere(aptr) )
                    result++;
			}
		}
	}

	BondIteratorMolSet bitr( pmset );
    if( DrawBonds )
        for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond()) 
            if( bptr->IsToDraw() && (!UseClipping ||
                !ClipVectBond(bptr->srcatom,bptr->dstatom)) )
                    result++;

	pmset->GetMolEditor()->UpdateBackBone(pmset);
	int nb = pmset->BackboneBonds.size();
	int ib;
	for( ib = 0; ib < nb; ib++ )
	{
		bptr = pmset->BackboneBonds[ib].get();
        if( bptr->IsToDraw() && (!UseClipping ||
            !ClipVectBond(bptr->srcatom,bptr->dstatom)) )
                result++;
	}

	
    for( hptr= pmset->HBonds.begin(); hptr != pmset->HBonds.end() ; hptr++ )
        if( (*hptr).IsToDraw() )
        {   
			if( HBondMode )
            {   
				if( !ClipVectBond((*hptr).srcCA, (*hptr).dstCA) )
                    result++;
            } 
			else if( !ClipVectBond((*hptr).src, (*hptr).dst) )
                result++;
        }

	list<Monitor>::iterator mitr;
    for( mitr=MonitList.begin(); mitr != MonitList.end(); mitr++ )
        if( !UseClipping || !ClipVectBond((*mitr).src,(*mitr).dst) )
            result++;

    return( result );
}


void 
HaMolView::FetchPSItems(PSItemPtr* data, char* type )
{
    HaChain  *chain;
    set<HaHBond, less<HaHBond> >::iterator hptr;
    HaBond  *bptr;
    HaAtom  *aptr;
    int i,flag;

    i = 0;

	MoleculesType::iterator mol_itr;

	MolSet* pmset = GetMolSet();

    if( DrawAtoms )
	{
		AtomIteratorMolSet aitr(pmset);
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
        {
            if( aptr->IsDrawSphere() )
			{
                if( !UseClipping || !ClipVectSphere(aptr) )
                {   
					type[i] = PSAtom; 
                    data[i++] = aptr;
                }
			}
		}
	}

	BondIteratorMolSet bitr(pmset);
    if( DrawBonds )
	{
        for(bptr = bitr.GetFirstBond();bptr;bptr = bitr.GetNextBond()) 
		{
            if( bptr->IsToDraw() && (!UseClipping ||
                !ClipVectBond(bptr->srcatom,bptr->dstatom)) )
            {   
				type[i] = PSBond;
                data[i++] = bptr;
            } 
		}
	}

	pmset->GetMolEditor()->UpdateBackBone(pmset);
	int nb = pmset->BackboneBonds.size();
	int ib;
	for( ib = 0; ib < nb; ib++ )
	{
	   bptr = pmset->BackboneBonds[ib].get();
       if( bptr->IsToDraw() && (!UseClipping || !ClipVectBond(bptr->srcatom,bptr->dstatom)) )
       {   
		   type[i] = PSBond;
           data[i++] = bptr; 
       }
	}
		
	   
	   for( hptr= pmset->HBonds.begin(); hptr != pmset->HBonds.end() ; hptr++ )
	   {
		   if( (*hptr).IsToDraw() )
		   {   
			   if( HBondMode )
			   {   
				   flag = !ClipVectBond((*hptr).srcCA,(*hptr).dstCA);
			   } 
			   else 
				   flag = !ClipVectBond((*hptr).src,(*hptr).dst);
			   
			   if( flag )
			   {   
				   type[i] = PSHBond;
				   data[i++] = (HaHBond*)&(*hptr);
			   }
		   }
	   }  

	list<Monitor>::iterator mitr;
    for( mitr=MonitList.begin(); mitr != MonitList.end(); mitr++ )
        if( !UseClipping || !ClipVectBond((*mitr).src,(*mitr).dst) )
        {   
			type[i] = PSMonit;
            data[i++] = &(*mitr);
        } 
}


void 
HaMolView::WritePSItems(PSItemPtr* data, char* type, int count )
{
    Monitor  *monit;
    HaHBond  *hbond;
    HaBond  *bond;
    HaAtom  *src;
    HaAtom  *dst;
    int i;

    for( i=0; i<count; i++ )
        switch( type[i] )
        {   
			case(PSAtom):   WriteVectSphere(data,type,i);
                            break;

            case(PSBond):   bond = (HaBond *)data[i];
                            src = bond->srcatom;
                            dst = bond->dstatom;

                            if( bond->flag & WireFlag )
                            {   
								WriteVectWire(src,dst,bond->col,False);
                            } 
							else if( bond->flag & CylinderFlag )
                            {   
								WriteVectStick(src,dst,bond->col,(int)(Scale*bond->radius));
                            } 
							else /* bond->flag & DashFlag */
                                WriteVectWire(src,dst,bond->col,True);
                            break;

            case(PSSSBond): 
            case(PSHBond):  hbond = (HaHBond *)data[i];
                            if( (type[i]==PSHBond)? HBondMode : SSBondMode )
                            {   
								src = hbond->srcCA;
                                dst = hbond->dstCA;
                            } 
							else
                            {   
								src = hbond->src;
                                dst = hbond->dst;
                            }

                            if( hbond->IsToDrawWire() )
                            {   
								WriteVectWire(src,dst,hbond->col,True);
                            } 
							else /* bond->flag & CylinderFlag */
                                WriteVectStick(src,dst,hbond->col,
                                                       (int)(Scale*hbond->radius));
                            break;



            case(PSMonit):  monit = (Monitor *)data[i];
                            WriteVectWire(monit->src,monit->dst,
                                          monit->col,True);
                            break;
        }
}



int 
HaMolView::WriteVectPSFile(const char* name )
{
    register double ambi;
    register double temp, inten;
    register int xsize, ysize;
    register int xpos, ypos;
    register int count;
    register int i;

    PSItemPtr  *data;
    char  *type;

    count = CountPSItems();
    if( !count ) return( True );

#ifdef _WIN32
    if( count > 16383 )
    {   
        PrintLog("\n Output Error: Too many PostScript objects!\n");
        return( False );
    }
#endif

    /* Allocate arrays for objects! */
    data = (PSItemPtr *)malloc((size_t)count*sizeof(PSItemPtr));
    type = (char *)malloc((size_t)count*sizeof(char));
    if( !data || !type )
    {   
        PrintLog("\n Output Error: Not enough memory to create PostScript!\n");

        if( data ) free( data );
        if( type ) free( type );
        return( False );
    }

    OutFile = fopen(name,"w");
    if( !OutFile )
    {   
		PrintLog("\n Output Error: Unable to create file %s!\n",name);
        return(False);
    }

    /* Determine the size of the image */
    ysize = (int)(pCanv->YRange()*(BORDER*PAGEWIDE)/pCanv->XRange());
    if( ysize > (int)(BORDER*PAGEHIGH) )
    {   
		xsize = (int)(pCanv->XRange()*(BORDER*PAGEHIGH)/pCanv->YRange());
        ysize = (int)(BORDER*PAGEHIGH);
    } 
	else 
		xsize = (int)(BORDER*PAGEWIDE);

    xpos = (int)(PAGEWIDE-xsize)/2;
    ypos = (int)(PAGEHIGH-ysize)/2;

    fputs("%!PS-Adobe-2.0 EPSF-2.0\n",OutFile);
    fputs("%%Creator: RasMol Version 2.6\n",OutFile);
    fprintf(OutFile,"%%%%Title: %s\n",name);
    fprintf(OutFile,"%%%%BoundingBox: %d %d ",xpos,ypos);
    fprintf(OutFile,"%d %d\n",xpos+xsize,ypos+ysize);

    fputs("%%Pages: 1\n",OutFile);
    fputs("%%EndComments\n",OutFile);
    fputs("%%EndProlog\n",OutFile);
    fputs("%%BeginSetup\n",OutFile);

    fputs("1 setlinecap 1 setlinejoin [] 0 setdash\n",OutFile);
    fputs("1 setlinewidth 0 setgray\n",OutFile);
    fputs("%%EndSetup\n",OutFile);
    fputs("%%Page: 1 1\n",OutFile);

    fputs("gsave\n",OutFile);
    fputs("14 dict begin\n\n",OutFile);
    fputs("/handleerror { showpage } def\n\n",OutFile);
    fputs("/Inten {\n  dup 4 index mul exch\n",OutFile);
    fputs("  dup 4 index mul exch\n",OutFile);
    fputs("  3 index mul setrgbcolor\n} def\n\n",OutFile);

    fputs("/Dot {\n  moveto 0 0 rlineto stroke\n} def\n\n",OutFile);
    fputs("/Wire {\n  moveto lineto stroke\n} def\n\n",OutFile);
    fputs("/Label {\n  moveto 1 -1 scale\n",OutFile);
    fputs("  show mtrx setmatrix\n} def\n\n",OutFile);

    if( UseOutLine )
    {   
		fputs("/Stick {\n  closepath gsave setrgbcolor fill\n",OutFile);
        fputs("  grestore 0 setgray stroke\n} def\n\n",OutFile);
    } 
	else
        fputs("/Stick {\n  closepath setrgbcolor fill\n} def\n\n",OutFile);

    fputs("/StickEnd {\n  matrix currentmatrix 6 1 roll\n",OutFile);
    fputs("  translate rotate 1 scale\n",OutFile);
    fputs("  0 0 3 2 roll 90 -90 arc\n  setmatrix\n} def\n\n",OutFile);

    if( UseOutLine )
    {   
		fputs("/Shade {\n  closepath gsave clip\n",OutFile);
    } 
	else 
		fputs("/Shade {\n  closepath clip\n",OutFile);

    if( Ambient < 0.99 )
    {   
		ambi = 0.5*Ambient;
        fputs("  45 rotate dup -0.81649658092 mul scale\n",OutFile);
        fprintf(OutFile,"  %g Inten fill\n",ambi);
        inten = (1.0-ambi)/31;
        for( i=0; i<31; i++ )
        {   
			temp = (double)(i+1)/32;
            fprintf(OutFile,"  0 %g ",(double)i/32);
            fprintf(OutFile,"%g 0 360 arc ",sqrt(1.0-temp*temp));
            fprintf(OutFile,"%g Inten fill\n",i*inten+ambi);
        }
        if( UseOutLine )
        {   
			fputs("  grestore 0 setgray stroke",OutFile);
        } 
		else 
			fputc(' ',OutFile);
        fputs(" pop pop pop\n} def\n\n",OutFile);

    } 
	else /* Saturated Colours! */
    {   
		fputs("  pop setrgbcolor fill\n",OutFile);
        if( UseOutLine )
            fputs("  grestore 0 setgray stroke\n",OutFile);
        fputs("} def\n\n",OutFile);
    }


    fputs("/ClipSphere {\n  translate 0 0 5 2 roll arc\n} def\n\n",OutFile);
    fputs("/ClipBox {\n  translate rotate\n  dup lineto dup neg ",OutFile);
    fputs("dup\n  0 rlineto 0 exch rlineto 0 rlineto closepath\n",OutFile);
    fputs("  clip newpath mtrx setmatrix\n} def\n\n",OutFile);
    fputs("/ClipEllips {\n  translate rotate 1 scale\n",OutFile);
    fputs("  0 0 4 2 roll dup neg arc\n",OutFile);
    fputs("  reversepath mtrx setmatrix\n} def\n\n",OutFile);

    fputs("/Sphere {\n  gsave\n",OutFile);
    fputs("  translate 0 0 2 index 0 360 arc\n",OutFile);
    if( UseOutLine )
    {   
		fputs("  gsave Shade grestore\n",OutFile);
        fputs("  0 setgray stroke\n",OutFile);
        fputs("  grestore\n} def\n\n",OutFile);
    } 
	else
        fputs("  Shade grestore\n} def\n\n",OutFile);

    fprintf(OutFile,"%d %d translate\n",xpos,ypos+ysize);
    fprintf(OutFile,"%g ",(double)xsize/pCanv->XRange());
    fprintf(OutFile,"%g ",(double)-ysize/pCanv->YRange());

    fputs("scale\n/mtrx matrix currentmatrix def\n\n",OutFile);

    fputs("newpath 0 0 moveto 0 ",OutFile);
    fprintf(OutFile,"%d rlineto %d 0 rlineto 0 %d",pCanv->YRange(),pCanv->XRange(),-pCanv->YRange());
    fputs(" rlineto\nclosepath clip ",OutFile);
    if( BackColor.r || BackColor.g || BackColor.b )
    {   
		fprintf(OutFile,"%g %g %g",BackColor.r/255.0,BackColor.g/255.0,BackColor.b/255.0);
        fputs(" setrgbcolor fill\n\n",OutFile);
    } 
	else 
		fputs("newpath\n\n",OutFile);

    LineWidth = 1.0;
    VectSolid = True;
    VectCol = 0;

    FetchPSItems(data,type);
    if( count>1 )
        DepthSort(data,type,(int)count);

    WritePSItems(data,type,(int)count);
 
    if( !VectSolid )
    {   
		fputs("[] 0 setdash\n",OutFile);
        VectSolid = True;
    }

    if( DrawDots )
        WriteVectDots();
    if( DrawMonitDistance && !MonitList.empty() )
        WriteVectMonitors();
    if( DrawLabels )
        WriteVectLabels();

    fputs("newpath 0 0 moveto 0 ",OutFile);
    fprintf(OutFile,"%d rlineto %d 0 rlineto 0 %d",pCanv->YRange(),pCanv->XRange(),-pCanv->YRange());
    fputs(" rlineto\nclosepath 0 setgray 1 setlinewidth stroke\n",OutFile);
    fputs("end grestore\nshowpage\n",OutFile);
    fputs("%%Trailer\n",OutFile);
    fputs("%%EOF\n",OutFile);

    fclose( OutFile );
    free( data );
    free( type );
    return(True);
}



void 
HaMolView::FlushPICTBuffer()
{
    if( PacketLen )
    {   
		if( RLEOutput )
        {   
			WriteByte(PacketLen-1);
            fwrite((char*)Buffer,1,PacketLen,OutFile);
        } 
		else 
			RLELineSize += PacketLen+1;
        PacketLen = 0;
    }
}


void 
HaMolView::FlushPICTPacket()
{
    register int i;

    if( RLELength>2 )
    {   
		FlushPICTBuffer();

        if( RLEOutput )
        {   
			WriteByte(257-RLELength);
            WriteByte(RLEChar);
        } 
		else 
			RLELineSize += 2;
    } 
	else 
        for( i=0; i<RLELength; i++ )
        {   
			Buffer[PacketLen++] = RLEChar;
            if( PacketLen == 128 ) 
                FlushPICTBuffer();
        }
}


void 
HaMolView::WritePICTCode(int val )
{
    if( !RLELength )
    {   
		RLEChar = val;
        RLELength = 1;
    } 
	else if( (val!=RLEChar) || (RLELength==128) )
    {   
		FlushPICTPacket();
        RLEChar = val;
        RLELength = 1;
    } 
	else 
		RLELength++;
}

void 
HaMolView::WritePICTData()
{
    ColorVal  *ptr;
    ColorVal  *tmp;
    int rowbytes;
    int x,y;

    rowbytes = pCanv->XRange()*3;

    RLEFileSize = 0;

    ptr = pCanv->FBuffer;

//    for( y = 0; y < pCanv->YRange(); y++ )
    for( y= pCanv->YRange()-1; y>=0; y-- )
    {
        RLELineSize = 0;
        RLEOutput = False;
        PacketLen = 0;
        RLELength = 0;

		ptr = pCanv->FBuffer + y* pCanv->XRange();

        tmp = ptr;
        for( x=0; x< pCanv->XRange(); x++ )
            WritePICTCode( (int)RComp(*tmp++) );
        tmp = ptr;
        for( x=0; x< pCanv->XRange(); x++ )
            WritePICTCode( (int)GComp(*tmp++) );
        tmp = ptr;
        for( x=0; x< pCanv->XRange(); x++ )
            WritePICTCode( (int)BComp(*tmp++) );

        FlushPICTPacket();
        FlushPICTBuffer();

        WriteMSBShort(RLELineSize);
        RLEFileSize += (RLELineSize+2);

        RLEOutput = True;
        PacketLen = 0;
        RLELength = 0;

        tmp = ptr;
        for( x=0; x < pCanv->XRange(); x++ )
            WritePICTCode( (int)RComp(*tmp++) );
        tmp = ptr;
        for( x=0; x < pCanv->XRange(); x++ )
            WritePICTCode( (int)GComp(*tmp++) );
		tmp = ptr;
        for( x=0; x < pCanv->XRange(); x++ )
            WritePICTCode( (int)BComp(*tmp++) );

        FlushPICTPacket();
        FlushPICTBuffer();
    }

}


int 
HaMolView::WritePICTFile(const char* name )
{
    int i;

#if defined(_WIN32) 
    OutFile = fopen(name,"wb");
#else
    OutFile = fopen(name,"w");
#endif
    if( !OutFile )
    {    
		 PrintLog("\n Output Error: Unable to create file %s!\n",name);
         return( False );
    }

    /* Write out header */
    for( i=0; i<512; i++ )
        WriteByte( 0x00 );

    WriteMSBShort(0);       /* picSize         */
    WriteMSBShort(0);       /* picFrame.top    */
    WriteMSBShort(0);       /* picFrame.left   */
    WriteMSBShort(pCanv->YRange());  /* picFrame.bottom */
    WriteMSBShort(pCanv->XRange());  /* picFrame.right  */

    WriteMSBShort(PICTpicversion);
    WriteMSBShort(0x02FF);

    WriteMSBShort(PICTheaderop);
    WriteMSBLong((unsigned int)0xffffffff);
    WriteMSBShort(0);      WriteMSBShort(0);
    WriteMSBShort(0);      WriteMSBShort(0);
    WriteMSBShort(pCanv->XRange()); WriteMSBShort(0);
    WriteMSBShort(pCanv->YRange()); WriteMSBShort(0);
    WriteMSBLong(0);

    WriteMSBShort(PICTcliprgn);
    WriteMSBShort(10);      /* rgnSize */
    WriteMSBShort(0);       /* rgnBBox.top    */
    WriteMSBShort(0);       /* rgnBBox.left   */
    WriteMSBShort(pCanv->YRange());  /* rgnBBox.bottom */
    WriteMSBShort(pCanv->XRange());  /* rgnBBox.right  */

    WriteMSBShort(PICTdirectbitsrect);
    WriteMSBShort(0x0000);  /* baseAddr      */
    WriteMSBShort(0x00ff);

    i = (pCanv->XRange()*sizeof(ColorVal)) | 0x8000;
    WriteMSBShort( i );     /* rowBytes      */
    WriteMSBShort(0);       /* bounds.top    */
    WriteMSBShort(0);       /* bounds.left   */
    WriteMSBShort(pCanv->YRange());  /* bounds.bottom */
    WriteMSBShort(pCanv->XRange());  /* bounds.right  */
    WriteMSBShort(0);       /* pmVersion     */

    WriteMSBShort(4);       /* packType      */
    WriteMSBLong(0);        /* packSize      */
    WriteMSBLong(72);       /* hRes          */
    WriteMSBLong(72);       /* vRes          */

    WriteMSBShort(16);      /* pixelType     */
    WriteMSBShort(32);      /* pixelSize     */
    WriteMSBShort(3);       /* cmpCount      */
    WriteMSBShort(8);       /* cmpSize       */

    WriteMSBLong(0);        /* planeBytes    */
    WriteMSBLong(0);        /* pmTable       */
    WriteMSBLong(0);        /* pmReserved    */


    WriteMSBShort(0);       /* srcRect.top    */
    WriteMSBShort(0);       /* srcRect.left   */
    WriteMSBShort(pCanv->YRange());  /* srcRect.bottom */
    WriteMSBShort(pCanv->XRange());  /* srcRect.right  */
    WriteMSBShort(0);       /* dstRect.top    */
    WriteMSBShort(0);       /* dstRect.left   */
    WriteMSBShort(pCanv->YRange());  /* dstRect.bottom */
    WriteMSBShort(pCanv->XRange());  /* dstRect.right  */
    WriteMSBShort(0);       /* mode (srcCopy) */

    WritePICTData();
    if( RLEFileSize & 0x01 ) WriteByte(0x00);
    WriteMSBShort(PICTendofpict);
    fclose(OutFile);
    return( True );
}


void 
HaMolView::FlushIRISBuffer()
{
    if( PacketLen )
    {   
		if( RLEOutput )
        {   
			WriteByte(PacketLen|0x80);
            fwrite((char*)Buffer,1,PacketLen,OutFile);
        } 
		else 
			RLELineSize += PacketLen+1;
        PacketLen = 0;
    }
}


void 
HaMolView::FlushIRISPacket()
{
    register int i;

    if( RLELength>2 )
    {   
		FlushIRISBuffer();

        if( RLEOutput )
        {   
			WriteByte(RLELength);
            WriteByte(RLEChar);
        } 
		else 
			RLELineSize += 2;
    } 
	else
        for( i=0; i<RLELength; i++ )
        {   
			Buffer[PacketLen++] = RLEChar;
            if( PacketLen == 127 )
                FlushIRISBuffer();
        }
    RLELength = 0;
}


void 
HaMolView::WriteIRISCode(int val )
{
    if( !RLELength )
    {   
		RLELength = 1;
        RLEChar = val;
    } 
	else if( (RLEChar!=val) || (RLELength==127) )
    {   
		FlushIRISPacket();
        RLELength = 1;
        RLEChar = val;
    } 
	else 
		RLELength++;
}


void 
HaMolView::DetermineIRISSizes(int* rowstart,short* rowsize, int* min, int* max )
{                    
    ColorVal  *ptr;
    ColorVal  *tmp;
    int i,x,y;

    *max = 0;
    *min = 255;
    RLEFileSize = 512 + 6*pCanv->YRange()*sizeof(int);

    RLEOutput = False;
    PacketLen = 0;
    RLELength = 0;

    for( y=0; y< pCanv->YRange(); y++ )
    {
        ptr = pCanv->FBuffer + (int)((pCanv->YRange()-1)-y)* pCanv->XRange();

        tmp = ptr;
        RLELineSize = 0;
        /* Red Component */
        for( x=0; x < pCanv->XRange(); x++ )
        {   i = RComp(*ptr++);
            if( i<*min ) *min=i;
            if( i>*max ) *max=i;
            WriteIRISCode(i);
        }
        FlushIRISPacket();
        FlushIRISBuffer();

        rowsize[y] = RLELineSize;
        rowstart[y] = RLEFileSize;
        RLEFileSize += RLELineSize;

        ptr = tmp;
        RLELineSize = 0;
        /* Green Component */
        for( x=0; x < pCanv->XRange(); x++ )
        {   
			i = GComp(*ptr++);
            if( i<*min ) *min=i;
            if( i>*max ) *max=i;
            WriteIRISCode(i);
        }
        FlushIRISPacket();
        FlushIRISBuffer();

        i = y+pCanv->YRange();
        rowsize[i] = RLELineSize;
        rowstart[i] = RLEFileSize;
        RLEFileSize += RLELineSize;

        ptr = tmp;
        RLELineSize = 0;
        /* Blue Component */
        for( x=0; x< pCanv->XRange(); x++ )
        {   
			i = BComp(*ptr++);
            if( i<*min ) *min=i;
            if( i>*max ) *max=i;
            WriteIRISCode(i);
        }
        FlushIRISPacket();
        FlushIRISBuffer();

        i = y+(pCanv->YRange()<<1);
        rowsize[i] = RLELineSize;
        rowstart[i] = RLEFileSize;
        RLEFileSize += RLELineSize;
    }
}

                             
void 
HaMolView::WriteIRISHeader(register int* rowstart, register short* rowsize, 
							int min, int max )
{              
    register int i,size;
    
    WriteMSBShort(474);     /* imagic     */
    WriteMSBShort(257);     /* type       */
    WriteMSBShort(3);       /* dim        */
    WriteMSBShort(pCanv->XRange());  /* xsize      */
    WriteMSBShort(pCanv->YRange());  /* ysize      */
    WriteMSBShort(3);       /* zsize      */
    WriteMSBLong(min);      /* min        */
    WriteMSBLong(max);      /* max        */
    WriteMSBLong(0);        /* wastebytes */

                            /* name       */
    fputs("RasMol IRIS RGB format output",OutFile);
    for( i=0; i<51; i++ ) WriteByte(0x00);
    WriteMSBLong(0);        /* colormap   */
                       
    size = 3*pCanv->YRange();
    for( i=0; i<404; i++ )  WriteByte(0x00);
    for( i=0; i<size; i++ ) WriteMSBLong(rowstart[i]);
    for( i=0; i<size; i++ ) WriteMSBLong(rowsize[i]);
}


void 
HaMolView::WriteIRISData()
{
    register ColorVal  *ptr;
    register ColorVal  *tmp;
    register int x,y,i;
	
    RLEOutput = True;
    PacketLen = 0;
    RLELength = 0;
	
	
    for( y=0; y< pCanv->YRange(); y++ )
    {
        ptr = pCanv->FBuffer + (int)((pCanv->YRange()-1)-y)*pCanv->XRange();

		
        tmp = ptr;
        /* Red Component */
        for( x=0; x< pCanv->XRange(); x++ )
        {   
			i = RComp(*ptr++);
            WriteIRISCode(i);
        }
        FlushIRISPacket();
        FlushIRISBuffer();
		
        ptr = tmp;
        /* Green Component */
        for( x=0; x< pCanv->XRange(); x++ )
        {   
			i = GComp(*ptr++);
            WriteIRISCode(i);
        }
        FlushIRISPacket();
        FlushIRISBuffer();
		
        ptr = tmp;
        /* Blue Component */
        for( x=0; x< pCanv->XRange(); x++ )
        {   
			i = BComp(*ptr++);
            WriteIRISCode(i);
        }
        FlushIRISPacket();
        FlushIRISBuffer();
    }
}

                             
int 
HaMolView::WriteIRISFile(const char* name )
{
    register int  *rowstart;
    register short  *rowsize;
    static int min,max;
    register int size;


#if defined(_WIN32) 
    OutFile = fopen(name,"wb");
#else
    OutFile = fopen(name,"w");
#endif
    if( !OutFile )
    {    
		 PrintLog("\n Output Error: Unable to create file %s!\n",name);
         return( False );
    }

    size = 3*pCanv->YRange();
    /* Allocate RLE encoded row length table */
    rowsize = (short *)malloc(size*sizeof(short));
    rowstart = (int *)malloc(size*sizeof(int));

    if( !rowsize || !rowstart )
    {   
        PrintLog("\n Output Error: Unable to allocate memory!\n");
        if( rowstart ) free(rowstart);
        if( rowsize ) free(rowsize);
        fclose(OutFile);
        return( False );
    }
    DetermineIRISSizes(rowstart,rowsize,&min,&max);
    WriteIRISHeader(rowstart,rowsize,min,max);    
    free( rowstart );
    free( rowsize );
                  
    WriteIRISData();
    fclose( OutFile );

    return( True );
}

