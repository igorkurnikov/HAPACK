/*! \file script.cpp
 
    Write Scripts to describe current image in HARLEM
    
	\author Igor Kurnikov
    \date 1998-2008
  
   derived from:
  
   RasMol2 Molecular Graphics
   Roger Sayle, August 1995
   Version 2.6
 */
#define SCRIPT
#include "haconst.h"

#include <stdlib.h> 

#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "command.h"
#include "abstree.h"
#include "hamolview.h"

#include "haatom.h"
#include "habond.h"
#include "hamolecule.h"
#include "moleditor.h"
#include "hasurface.h"
#include "etcoupl.h"

#define InvertY(y) (-(y))

#define Round(x)       ((int)(x))
#define DatWirFlag  (int)0x10000
#define DatDasFlag  (int)0x20000
#define DatCylFlag  (int)0x40000

typedef struct {
        int datum;
        int count;
    } FreqEntry;

#define FREQSIZE  8
static FreqEntry Freq[FREQSIZE];

static HaAtom  *MagePrev;
static char *MageCol;
static FILE *OutFile;
static int SelectAll;

static void IncFreqTable( int );
static int GetBondDatum( HaBond * );
static int GetHBondDatum( HaHBond * );


static void ResetFreqTable()
{
    register int i;

    for( i=0; i<FREQSIZE; i++ )
        Freq[i].count = 0;
}

static void IncFreqTable(int datum )
{
    register int count;
    register int i;

    for( i=0; i<FREQSIZE; i++ )
        if( !Freq[i].count )
        {   Freq[i].datum = datum;
            Freq[i].count = 1;
            return;
        } else if( Freq[i].datum == datum )
        {   count = Freq[i].count+1;
            while( i && (Freq[i-1].count<=count) )
            {   Freq[i] = Freq[i-1];  
                i--;
            }
            Freq[i].datum = datum;
            Freq[i].count = count;
            return;
        }

    /* Replace Singletons! */
    if( Freq[FREQSIZE-1].count == 1 )
        Freq[FREQSIZE-1].datum = datum;
}


static int GetBondDatum( HaBond* bptr )
{
    if( bptr->flag & CylinderFlag )
    {   
		return( DatCylFlag | (int)(bptr->radius*250)  );

    } 
	else if( bptr->flag & WireFlag )
    {   
		return( DatWirFlag );
    } 
	else if( bptr->flag & DashFlag )
    {   
		return( DatDasFlag );
    } 
	else return( (int)0 );
}

static int GetHBondDatum( HaHBond* bptr )
{
    if( bptr->IsToDrawCylinder() )
    {   
		return( DatCylFlag | (int)(bptr->radius*250.0) );
    } 
	else if( bptr->IsToDrawWire() )
    {   
		return( DatWirFlag );
    } 
	else 
		return( (int)0 );
}

 
#ifndef SCRIPT

static void WriteMolScriptAtomSel( HaChain *, HaResidue *, HaAtom * );

static void WriteMolScriptColour( int r, int g, int b )
{
    fprintf(OutFile," rgb %#g %#g %#g",r/255.0,g/255.0,b/255.0);
}

static void WriteMolScriptAtomSel(HaChain* chain, HaResidue* group, 
								  HaAtom* aptr )
{
    register char *ptr;
    register int i;

    fputs("require atom ",OutFile);
    ptr = ElemDesc[aptr->refno];
    for( i=0; i<4; i++ )
        if( ptr[i]=='*' )
        {   fputc('\'',OutFile);
        } else if( ptr[i]!=' ' )
            fputc(ptr[i],OutFile);

    fputs(" and in residue ",OutFile);
    if( chain->ident!=' ' && !isdigit(chain->ident) )
        fputc(chain->ident,OutFile);
    fprintf(OutFile,"%d",group->serno);
}

static void WriteMolScriptAtoms()
{
    HaChain  *chain;
    HaResidue  *group;
    HaAtom  *aptr;
    char *ptr;

	int r,g,b;
    char buffer[80];

    for(aptr= CurMolSet->GetFirstAtom(); aptr; aptr= CurMolSet->GetNextAtom())
    {   if( aptr->IsDrawSphere() )
        {   /* Atom Colour */
            fputs("set atomcolour ",OutFile);
            WriteMolScriptAtomSel(chain,group,aptr);
			HaColor::GetPureRGBbyCIdx(aptr->col,r,g,b);
            WriteMolScriptColour(r,g,b);
            fputs(";\n",OutFile);

            /* CPK Sphere */
            fputs("set atomradius ",OutFile);
            WriteMolScriptAtomSel(chain,group,aptr);
            fprintf(OutFile," %#g;\ncpk ",aptr->image_radius/250.0);
            WriteMolScriptAtomSel(chain,group,aptr);
            fputs(";\n",OutFile);
        }

        if( aptr->label )
        {   /* Atom Label */
            FormatLabel(chain,group,aptr,aptr->label.c_str(),buffer);

            fputs("label ",OutFile);
            WriteMolScriptAtomSel(chain,group,aptr);
            fputs(" \"",OutFile);
            for( ptr=buffer; *ptr; ptr++ )
                if( *ptr!='%' ) fputc(*ptr,OutFile);
            fputs("\";\n",OutFile);
        }
    }
}
#endif

static void MolScriptSegment( char* ptr, int src, int dst, char chain )
{   
    if( (chain!=' ') && !isdigit(chain) ) 
    {   fprintf(OutFile,"  %s from %c%d to %c%d;\n",ptr,chain,src,chain,dst);
    } else fprintf(OutFile,"  %s from %d to %d;\n",ptr,src,dst);
}


int 
HaMolView::WriteMolScriptFile( const char* name )
{
    double temp;
    double psi, phi, theta;
    HaChain  *chain;
    HaResidue  *group;
    HaResidue  *next;
    HaResidue  *prev;
    int flag,len;
    char *ptr;

	HaMolSet* pmset = GetMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);

    OutFile = fopen(name,"w");
    if( !OutFile )
    {   
		PrintLog("Script Error: Unable to create file %s \n",name);
        return(False);
    }
    fprintf(OutFile,"! File: %s\n",name);
    fputs("! Creator: HARLEM Version 1.0\n",OutFile);
    fputs("! Version: MolScript v1.3\n\n",OutFile);

    fputs("plot\n",OutFile);
    if( BackColor.r || BackColor.g || BackColor.b )
    {   
		fputs("  background rgb ",OutFile);
        fprintf(OutFile,"%#g ",  BackColor.r/255.0);
        fprintf(OutFile,"%#g ",  BackColor.g/255.0);
        fprintf(OutFile,"%#g;\n",BackColor.b/255.0);
    }
    temp = 0.004/Scale;
    fprintf(OutFile,"  window %g;\n",temp* pCanv->Range());
    if( UseSlabPlane() )
        fprintf(OutFile,"  slab %g;\n", CurSlabValue); // this is wrong
    fputc('\n',OutFile);

    fprintf(OutFile,"  read mol \"%s\";\n", pmset->GetName());
    fputs("  transform atom *\n",OutFile);
    fputs("    by centre position atom *\n",OutFile);
    fputs("    by rotation x 180.0",OutFile);

    phi = RAD_TO_DEG*asin(Rot(3,1));
    if( (int)phi == 90 )
    {   
		theta = -RAD_TO_DEG*atan2(Rot(1,2),Rot(2,2) );
        psi = 0;
    } 
	else if( (int)phi == -90 )
    {   
		theta = RAD_TO_DEG*atan2(Rot(1,2),Rot(2,2));
        psi = 0;
    } 
	else /* General Case! */
    {   
		theta = RAD_TO_DEG*atan2(Rot(3,2), Rot(3,3));
        psi =  -RAD_TO_DEG*atan2(Rot(2,1), Rot(1,1));
    }

    if( (int)psi )   fprintf(OutFile,"\n    by rotation z %#g",InvertY(psi));
    if( (int)phi )   fprintf(OutFile,"\n    by rotation y %#g",phi);
    if( (int)theta ) fprintf(OutFile,"\n    by rotation x %#g",InvertY(-theta));

    if( UseSlabPlane() )
    {   
		fputs("\n    by translation ",OutFile);
        fprintf(OutFile,"%#g ",-pCanv->WRange()*temp);
        fprintf(OutFile,"%#g ", pCanv->HRange()*temp);
        if( UseSlabPlane() )
        {   
            fprintf(OutFile,"%#g", CurSlabValue); // this is wrong
        } 
		else 
			fputs("0.0",OutFile);
    }
    fputs(";\n\n",OutFile);

    /* fputs("  trace amino-acids;\n",OutFile); */

    MoleculesType::iterator mol_itr;
	mol_itr = pmset->HostMolecules.begin();
	for(; mol_itr != pmset->HostMolecules.end(); mol_itr++)
	{
		if( (*mol_itr)->GetNChains() > 0)
		{   
			if( !(*mol_itr)->IsSecStructFound()   )
				p_mol_editor->DetermineSecStructure( (*mol_itr), False );
		}
	}
	ChainIteratorMolSet chitr(this->GetMolSet() );
    for( chain = chitr.GetFirstChain(); chain; chain = chitr.GetNextChain() )
    {   
		prev = (HaResidue *)0;
		ResidueIteratorChain ritr_ch(chain);
		for( group=ritr_ch.GetFirstRes(); group && group->GetNextResInChain(); group=next )
		{   
			next = group->GetNextResInChain();
			if( next->serno < group->serno )
			{   
				if( prev && prev!=group )
					MolScriptSegment("coil",prev->serno,group->serno,
					chain->ident);
				prev = (HaResidue *)0;
				continue;
			}
			flag = group->struc & next->struc;
			
			if( flag&HelixFlag )
			{   
				flag = HelixFlag;
				ptr = "helix";
			} 
			else if( flag&SheetFlag )
			{   
				flag = SheetFlag;
				ptr = "strand";
			} 
			else 
			{   
				if( flag&TurnFlag )
				{   
					fputs("  turn residue ",OutFile);
					if( chain->ident != ' ' )
						fputc(chain->ident,OutFile);
					fprintf(OutFile,"%d;\n",group->serno);
				}
				if( !prev ) prev = group;
				continue;
			}
			
			len = 2;  /* Determine Structure Length */
			while( next->GetNextResInChain() 
				&& ( (next->GetNextResInChain())->struc&flag )
				&& (next->serno <= (next->GetNextResInChain())->serno) )
			{   
				next = next->GetNextResInChain();
				len++;
			}
			
			if( len>2 )
			{   
				if( prev && prev!=group ) /* MolScript coil or turn? */
					MolScriptSegment("coil",prev->serno,group->serno,
					chain->ident);
				MolScriptSegment(ptr,group->serno,next->serno,
					chain->ident);
				prev = next;
			} 
		}
		
		if( prev && prev!=group )  /* C-terminal coil/turn */
			MolScriptSegment("coil",prev->serno,group->serno,chain->ident);
	}
  

#ifndef SCRIPT
    WriteMolScriptAtoms();
#endif

    fputs("end_plot\n",OutFile);
    fclose(OutFile);
    return( True );
}




static void WriteScriptDatum( char*, int );
static void WriteScriptSelectBond( HaAtom *, HaAtom * );
static void WriteScriptHBonds( char*, list<HaHBond>& lst );


static void WriteScriptAll()
{
    if( !SelectAll )
    {   fputs("select all\n",OutFile);
        SelectAll = True;
    }
}

static void WriteScriptColour( char* ptr, int col )
{
    int r,g,b;
	HaColor::GetPureRGBbyCIdx(col,r,g,b);
    if( col )
    {   
        fprintf(OutFile,"colour %s [%d,%d,%d]\n",ptr,r,g,b);
    } 
	else 
		fprintf(OutFile,"colour %s none\n",ptr);
}


static void WriteScriptBetween( int lo, int hi )
{
    if( lo != hi )
    {   
		fprintf(OutFile,"select (atomno>=%d) and (atomno<=%d)\n",lo,hi);
    } 
	else 
		fprintf(OutFile,"select atomno=%d\n",lo);
    SelectAll = False;
}

static void WriteScriptSelectBond( HaAtom* src, HaAtom* dst )
{
    fprintf(OutFile,"select (atomno=%d) or (atomno==%d)\n",
                    src->GetSerNo(), dst->GetSerNo());
    SelectAll = False;
}


void 
HaMolView::WriteScriptAtoms()
{
    HaAtom  *aptr;
    int first,last;
    int same,init;
    int cpk;
    int col,rad;

    fputs("\n# Atoms\n",OutFile);

    same = True;
    init = False;

	HaMolSet* pmset = GetMolSet();

	AtomIteratorMolSet aitr(pmset);
    for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
        if( !init )
        {   
			first = last = aptr->GetSerNo();
            cpk = IsCPKColour( aptr );
            col = aptr->col;
            init = True;
        } else if( cpk && IsCPKColour(aptr) )
        {   last = aptr->GetSerNo();
            if( aptr->col != col )
                col = 0;
        } else if( aptr->col == col )
        {   last = aptr->GetSerNo();
            cpk = False;
        } else if( aptr->col != col )
        {   WriteScriptBetween( first, last );
            if( !col )
            {   fputs("colour atoms cpk\n",OutFile);
            } else WriteScriptColour("atoms",col);
                
            first = last = aptr->GetSerNo();
            cpk = IsCPKColour( aptr );
            col = aptr->col;
            same = False;
        } else last = aptr->GetSerNo(); 
        
    if( init )
    {   if( !same )
        {   WriteScriptBetween(first,last);
        } else WriteScriptAll();

        if( !col )
        {   fputs("colour atoms cpk\n",OutFile);
        } else WriteScriptColour("atoms",col);
    }

    if( DrawAtoms )
    {   
		same = True;
        init = False;
		AtomIteratorMolSet aitr(pmset);
		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
            if( !init )
            {   
				rad = aptr->IsDrawSphere() ? (int)(Scale*aptr->image_radius) : 0;
                first = last = aptr->GetSerNo();
                init = True;
            } 
			else if( rad == ((aptr->IsDrawSphere())? (int)(Scale*aptr->image_radius) : 0) )
            {   
                last = aptr->GetSerNo();
            } 
			else 
            {   
				WriteScriptBetween(first,last);
                if( rad == -1 )
                {   
					fputs("spacefill on\n",OutFile);
                } 
				else if( rad )
                {   
					fprintf(OutFile,"spacefill %d\n",rad);
                } 
				else 
					fputs("spacefill off\n",OutFile); 
				
                rad = aptr->IsDrawSphere() ? (int)(Scale*aptr->image_radius) : 0;
                first = last = aptr->GetSerNo();
                same = False;
            }
		}
		
        if( !same )
        {   
			WriteScriptBetween(first,last);
        } 
		else 
			WriteScriptAll();

        if( rad == -1 )
        {   
			fputs("spacefill on\n",OutFile);
        } 
		else if( rad )
        {   
			fprintf(OutFile,"spacefill %d\n",rad);
        } 
		else 
			fputs("spacefill off\n",OutFile); 

    } 
	else
    {   
		WriteScriptAll();
        fputs("spacefill off\n",OutFile);
    }
        
}


static void WriteScriptDatum( char* ptr, int datum )
{
    if( datum & DatCylFlag )
    {   
		fprintf(OutFile,"%s %d\n",ptr,(int)(datum-DatCylFlag));
    } 
	else if( datum & DatWirFlag )
    {   
		fprintf(OutFile,"%s on\n",ptr);
    } 
	else if( datum & DatDasFlag )
    {   
		fprintf(OutFile,"%s dash\n",ptr);
    } 
	else 
		fprintf(OutFile,"%s off\n",ptr);
}


void 
HaMolView::WriteScriptBonds()
{
    register HaBond  *bptr;
    register int defdat;
    register int datum;
    register int col;

    fputs("\n# Bonds\n",OutFile);

    ResetFreqTable();

	BondIteratorMolSet bitr(GetMolSet());
    for(bptr= bitr.GetFirstBond();bptr; bptr=bitr.GetNextBond() )
        IncFreqTable(GetBondDatum(bptr));

    WriteScriptAll();
    defdat = Freq[0].datum;
    WriteScriptDatum("wireframe",defdat);

    if( Freq[1].count )
    {   
		for(bptr= bitr.GetFirstBond();bptr; bptr= bitr.GetNextBond() )
        {   
			datum = GetBondDatum(bptr);
            if( datum != defdat )
            {    
				WriteScriptSelectBond(bptr->srcatom,bptr->dstatom);
                 WriteScriptDatum("wireframe",datum);
            }
        }
    } else if( !defdat )
        return;

    ResetFreqTable();
    for(bptr=bitr.GetFirstBond();bptr; bptr=bitr.GetNextBond() )
        IncFreqTable(bptr->col);

    col = (int)Freq[0].datum;
    if( col )
    {   WriteScriptAll();
        WriteScriptColour("bonds",col);
    }

    if( Freq[1].count )
        for(bptr=bitr.GetFirstBond();bptr; bptr= bitr.GetNextBond() )
            if( bptr->col != col )
            {   WriteScriptSelectBond(bptr->srcatom,bptr->dstatom);
                WriteScriptColour("bonds",bptr->col);
            }
}

void HaMolView::WriteScriptBackbone()
{
    HaChain  *chain;
    HaBond  *bptr;

    int defdat;
    int datum;
    int col;

    fputs("\n# Backbone\n",OutFile);

    ResetFreqTable();

	int i;
	HaMolSet* pmset = GetMolSet();
	pmset->GetMolEditor()->UpdateBackBone(pmset);
	int nb = pmset->BackboneBonds.size();
	int ib;
	for( ib = 0; ib < nb; ib++ )
	{
		bptr = pmset->BackboneBonds[ib];
		IncFreqTable(GetBondDatum(bptr));
	}
    WriteScriptAll();
    defdat = Freq[0].datum;
    WriteScriptDatum("backbone",defdat);

	if( Freq[1].count )
	{   
		for( ib = 0; ib < nb; ib++ )
		{
			bptr = pmset->BackboneBonds[ib];
			datum = GetBondDatum(bptr);
			if( datum != defdat )
			{    
				WriteScriptSelectBond(bptr->srcatom,bptr->dstatom);
				WriteScriptDatum("backbone",datum);
			}
		}
	} 
	else if( !defdat )
        return;

    ResetFreqTable();

	for( ib = 0; ib < nb; ib++ )
	{
		bptr = pmset->BackboneBonds[ib];
        IncFreqTable(bptr->col);
	}

    col = (int)Freq[0].datum;
    if( col )
    {   
		WriteScriptAll();
        WriteScriptColour("backbone",col);
    }

	if( Freq[1].count )
	{
		for( ib = 0; ib < nb; ib++ )
		{
			bptr = pmset->BackboneBonds[ib];
			if( bptr->col != col )
			{   
				WriteScriptSelectBond(bptr->srcatom,bptr->dstatom);
				WriteScriptColour("backbone",bptr->col);
			}
		}
	}
}

void HaMolView::WriteScriptRibbons()
{
    HaChain  *chain;
    HaResidue  *group;
    HaAtom  *aptr;
	
    fputs("\n# Ribbons\n",OutFile);
	
    if(DrawRibbon )
    {   
		fprintf(OutFile,"set strands %d\n", pCanv->m_SplineCount);
		fprintf(OutFile,"set cartoon %s\n", DrawBetaArrows?"on":"off");
		fprintf(OutFile,"set cartoon %g\n", CartoonHeight);

		MoleculesType::iterator mol_itr;
		ResidueIteratorMolSet ritr(this->GetMolSet() );
		for( group = ritr.GetFirstRes(); group; group = ritr.GetNextRes() )
		{    
			if( group->IsAmino() )
			{   
				aptr = group->GetAtomByName("CA");
			} 
			else 
				aptr = group->GetAtomByName("P");
			if( !aptr ) continue;

			fprintf(OutFile,"select atomno=%d\n",aptr->GetSerNo());
			SelectAll = False;

			if( group->flag & RibbonFlag )
			{   
				fprintf(OutFile,"ribbons %d\n",group->width);
			} 
			else if( group->flag & CartoonFlag )
			{   
				fprintf(OutFile,"cartoon %d\n",group->width);
			} 
			else if( group->flag & StrandFlag )
			{   
				fprintf(OutFile,"strands %d\n",group->width);
			} 
			else if( group->flag & DashStrandFlag )
			{   
				fprintf(OutFile,"strands dash %d\n",group->width);
			} 
			else if( group->flag & TraceFlag )
			{   
				fprintf(OutFile,"trace %d\n",group->width);
			} 
			else fputs("ribbons off\n",OutFile);

			if( group->col1 != group->col2 )
			{   
				if( group->col1 )
					WriteScriptColour("ribbon1",group->col1);
				if( group->col2 )
					WriteScriptColour("ribbon2",group->col2);
			} 
			else if( group->col1 )
				WriteScriptColour("ribbons",group->col1);
		}
	} 
	else
	{   
		WriteScriptAll();
		fputs("ribbons off\n",OutFile);
	}
}


void 
HaMolView::WriteScriptHBonds( char* obj )
{
    int defdat;
    int datum;
    int col;

	HaMolSet* pmset = GetMolSet();
	HBondIteratorMolSet hbitr(pmset);

	HaHBond* ptr;

    ResetFreqTable();
    for( ptr= hbitr.GetFirstBond(); ptr ; ptr = hbitr.GetNextBond() )
        IncFreqTable(GetHBondDatum(&(*ptr)));

    WriteScriptAll();
    defdat = Freq[0].datum;
    WriteScriptDatum(obj,defdat);

    if( Freq[1].count )
    {   
		for( ptr= hbitr.GetFirstBond(); ptr ; ptr = hbitr.GetNextBond() )
        {   
			datum = GetHBondDatum(&(*ptr));
            if( datum != defdat )
            {    
				WriteScriptSelectBond((*ptr).src,(*ptr).dst);
                WriteScriptDatum(obj,datum);
            }
        }
    } 
	else if( !defdat )
        return;

    ResetFreqTable();
    for( ptr= hbitr.GetFirstBond(); ptr ; ptr = hbitr.GetNextBond() )
        IncFreqTable((*ptr).col);

    col = (int)Freq[0].datum;
    if( col )
    {   
		WriteScriptAll();
        WriteScriptColour(obj,col);
    }

    if( Freq[1].count )
        for( ptr= hbitr.GetFirstBond(); ptr ; ptr = hbitr.GetNextBond()  )
            if( (*ptr).col != col )
            {   
		WriteScriptSelectBond((*ptr).src,(*ptr).dst);
                WriteScriptColour(obj,(*ptr).col);
            }
}

void 
HaMolView::WriteScriptLabels()
{
    HaAtom  *aptr;
    int first,last;

    fputs("\n# Labels\n",OutFile);
    WriteScriptAll();
    fputs("labels off\n",OutFile);
    if( DrawLabels ) return;

    if( UseLabelCol )
    {   
		fprintf(OutFile,"colour labels [%d,%d,%d]\n",
			             LabelColor.r,LabelColor.g,LabelColor.b);
    } 
	else 
		fputs("colour labels none\n",OutFile);
    fprintf(OutFile,"set fontsize %d\n",pCanv->m_FontSize);

	std::string label;
	
	AtomIteratorMolSet aitr(GetMolSet());
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
        if( aptr->label != label )
        {   
			if( !label.empty() )
            {   
				WriteScriptBetween(first,last);
                fprintf(OutFile,"label \"%s\"\n",label.c_str());
            }
            label = aptr->label;
            first = last = aptr->GetSerNo();
        } 
		else 
			last = aptr->GetSerNo();
	}

    if( !label.empty() )
    {   
		WriteScriptBetween(first,last);
        fprintf(OutFile,"label \"%s\"", aptr->label.c_str());
    }
}


void 
HaMolView::WriteScriptMonitors()
{
    register int col;

    fputs("\n# Monitors\n",OutFile);
    if( MonitList.empty() )
    {   fputs("monitors off\n",OutFile);
        return;
    }

    fprintf(OutFile,"set monitors %s\n",DrawMonitDistance?"on":"off");

    ResetFreqTable();

	list<Monitor>::iterator mitr;
    for( mitr= MonitList.begin(); mitr!= MonitList.end(); mitr++ )
    {   

		fprintf(OutFile,"monitor %d %d\n",(*mitr).src->GetSerNo(),

			    (*mitr).dst->GetSerNo());
        IncFreqTable((*mitr).col);
    }

    col = (int)Freq[0].datum;
    if( col )
    {   WriteScriptAll();
        WriteScriptColour("monitors",col);
    }

    if( Freq[1].count )
        for( mitr= MonitList.begin(); mitr != MonitList.end(); mitr++ )
            if( (*mitr).col != col )
            {   
				WriteScriptSelectBond((*mitr).src,(*mitr).dst);
                WriteScriptColour("monitor",(*mitr).col);
            }
}


int 
HaMolView::WriteScriptFile( const char* name )
{
    int theta,phi,psi;
    char *ptr;
    int temp;

    OutFile = fopen(name,"w");
    if( !OutFile )
    {   PrintLog("Script Error: Unable to create file %s \n",name);
        return(False);
    }

    fprintf(OutFile,"#!rasmol -script\n# File: %s\n",name);
    fputs("# Creator: RasMol Version 2.6\n\n",OutFile);
    fputs("zap\n",OutFile);
    fprintf(OutFile,"background [%d,%d,%d]\n", BackColor.r,BackColor.g,BackColor.b);

	HaMolSet* pmset = GetMolSet();

    if( !pmset )
    {   /* No Molecule! */
        fclose(OutFile);
        return(True);
    }

    /* Molecule File Name */
    fprintf(OutFile,"load %s \"%s.pdb\"\n","pdb", pmset->GetName());

    /* Colour Details */
    fprintf(OutFile,"set ambient %d\n", (int)(100*Ambient) );
    fputs("set specular ",OutFile);
    if( FakeSpecular )
    {   
		fprintf(OutFile,"on\nset specpower %d\n", SpecPower);
    } 
	else 
		fputs("off\n",OutFile);
    putc('\n',OutFile);

    /* Transformation */
    fputs("reset\n",OutFile);
    if( UseSlabPlane() )
    {   temp = (int)(50.0* CurSlabValue);
        if( temp )
        {   fprintf(OutFile,"slab %d\n",temp+50);
        } else fputs("slab on\n",OutFile);

        fputs("set slabmode ",OutFile);
        switch( SlabMode() )
        {   default:            
            case(SlabClose):    ptr = "solid";    break;
            case(SlabReject):   ptr = "reject";   break;
            case(SlabHalf):     ptr = "half";     break;
            case(SlabHollow):   ptr = "hollow";   break;
            case(SlabSection):  ptr = "section";
        }
        fputs(ptr,OutFile);
        putc('\n',OutFile);
    } else fputs("slab off\n",OutFile);

    phi = Round(RAD_TO_DEG*asin(Rot(3,1)));
    if( phi == 90 )
    {   
		theta = -Round(RAD_TO_DEG*atan2(Rot(1,2),Rot(2,2)));
        psi = 0;
    } 
	else if( phi == -90 )
    {   
		theta = Round(RAD_TO_DEG*atan2(Rot(1,2),Rot(2,2)));
        psi = 0;
    } 
	else /* General Case! */
    {   
		theta = Round(RAD_TO_DEG*atan2(Rot(3,2),Rot(3,3)));
        psi =  Round(-RAD_TO_DEG*atan2(Rot(2,1),Rot(1,1)));
    }

    if( psi )   fprintf(OutFile,"rotate z %d\n",InvertY(-psi));
    if( phi )   fprintf(OutFile,"rotate y %d\n",phi);
    if( theta ) fprintf(OutFile,"rotate x %d\n",InvertY(-theta));

    temp = (int)(100.0*CurTX);
    if( temp ) fprintf(OutFile,"translate x %d\n",temp);
    temp = (int)(100.0*CurTY);
    if( temp ) fprintf(OutFile,"translate y %d\n",InvertY(-temp));

    if( Zoom != 1.0 )
    {   
        fprintf(OutFile,"zoom %f\n",Zoom);
    }
    putc('\n',OutFile);

    /* Rendering */
    if( DrawAxes || DrawBoundBox || DrawUnitCell )
        fprintf(OutFile,"colour axes [%d,%d,%d]\n", BoxColor.r, BoxColor.g, BoxColor.b);
    if( DrawBonds )
        fprintf(OutFile,"set bonds %s\n", DrawDoubleBonds? "on":"off" );

    fprintf(OutFile,"set axes %s\n", DrawAxes? "on":"off" );
    fprintf(OutFile,"set boundingbox %s\n", DrawBoundBox? "on":"off" );
    fprintf(OutFile,"set unitcell %s\n", DrawUnitCell? "on":"off" );

    if( pmset->HBonds.empty() )
    {   
		fputs("set hbond ",OutFile);
        fputs( HBondMode?"backbone":"sidechain",OutFile);
        fputc('\n',OutFile);
    }

    fputs("set bondmode and\ndots off\n\n",OutFile); 
    fputs("\n# Avoid Colour Problems!\nselect all\n",OutFile);
    fputs("colour bonds none\ncolour backbone none\n",OutFile);
    fputs("colour hbonds none\ncolour ssbonds none\n",OutFile);
    fputs("colour ribbons none\ncolour white\n",OutFile);
    SelectAll = True;

    WriteScriptAtoms();
    if( UseSlabPlane() && ( SlabMode() ==SlabSection) )
    {   /* Section Mode Slabbing! */
        fclose(OutFile);
        return(True);
    }

    WriteScriptBonds();
    WriteScriptRibbons();
    WriteScriptBackbone();
    WriteScriptLabels();
    WriteScriptMonitors();
    fputc('\n',OutFile);
    
    WriteScriptHBonds("ssbonds");
    WriteScriptHBonds("hbonds");
    WriteScriptAll();

    fclose(OutFile);
   return( True );
}


int 
HaMolView::WritePOVRayFile( const char* name )
//!
//! Render POV file, Modifications from RASMOL Code thanks to Lev Gelb
//!
{
    HaAtom  *aptr;
    double x,y,z;

    OutFile = fopen(name,"w");
    if( !OutFile )
    {   
		ErrorInMod("HaMolView::WritePOVRayFile","Can't open POV file for writing");
        return(False);
    }

    fprintf(OutFile,"// File: %s\n",name);
    fputs("// Creator: Harlem Version 0.2\n",OutFile);
    fputs("// Version: POV-Ray Version 3.1\n\n",OutFile);

    fputs("#include \"shapes.inc\"\n",OutFile);
    fputs("#include \"colors.inc\"\n",OutFile);
    fputs("#include \"textures.inc\"\n\n",OutFile);

    fputs("// Camera\ncamera {",OutFile);

	double dd= 1.2*pCanv->Range();
	double dx= pCanv->XRange()/2.0;
	double dy= pCanv->YRange()/2.0;
	double zero = 0.0;

    fprintf(OutFile,"    location < %12.6f ,  %12.6f, %12.6f>\n", dx, dy, -2*dd);

    fprintf(OutFile,"    look_at <%12.6f, %12.6f, %12.6f>\n  angle 40 \n orthographic }\n\n",dx,dy,zero );

    fprintf(OutFile,"// Light\nlight_source {<%12.6f, %12.6f, %12.6f>",dx, dy, -3*dd);
//    fputs(" color red 0.9 green 0.9 blue 0.9 }\n\n",OutFile);
	fputs(" color rgb < 0.9, 0.9, 0.9 > }\n\n",OutFile);

    fprintf(OutFile,"// Light\nlight_source {<%12.6f, %12.6f, %12.6f>",dx, dy - 3*dd, -3*dd);
//    fputs(" color red 0.5 green 0.5 blue 0.5 }\n\n",OutFile);
	fputs(" color rgb < 0.5, 0.5, 0.5 > shadowless }\n\n",OutFile);

    fprintf(OutFile,"// Light\nlight_source {<%12.6f, %12.6f, %12.6f>",dx - 3*dd, dy, -3*dd);
//    fputs(" color red 0.5 green 0.5 blue 0.5 }\n\n",OutFile);
	fputs(" color rgb < 0.5, 0.5, 0.5 > shadowless }\n\n",OutFile);

    fputs("#default { \n",OutFile);
    fputs("texture {  \n",OutFile);
    fputs("    finish \n",OutFile);
    fputs("   {  diffuse 0.6       \n",OutFile);
    fputs("      brilliance 1.4      \n",OutFile);
    fputs("      ambient 0.2  reflection 0   \n",OutFile);
    fputs("      phong 0.2 phong_size 40       \n",OutFile);
    fputs("      crand 0                     \n",OutFile);
    fputs("      specular 0  roughness 0.05  \n",OutFile);
    fputs("    }\n  }\n }\n",OutFile);


	fprintf(OutFile," background { rgb <%12.6f, %12.6f, %12.6f> } \n\n", 
		             (double)BackColor.r/255.0, (double)BackColor.g/255.0, (double)BackColor.b/255.0 );

	fputs("#declare AtomFinish = finish { specular 0.2 }\n\n", OutFile);

	
	double zmin= 1000000.0;

	AtomIteratorMolSet aitr(this->GetMolSet());
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if((double)aptr->z < zmin) zmin = (double)aptr->z; 
	}
	
	fputs("// Objects\n",OutFile);

	if(DrawAtoms)
	{
		fputs("// Draw Atom Spheres \n", OutFile);
		
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			if( aptr->IsDrawSphere() )
			{  
				x = (double)( aptr->x );
				y = (double)( aptr->y );
				z = (double)( aptr->z - zmin);
				
				double srad;
				if( aptr->irad != 0 )
					srad = (double) aptr->irad;
				else
					srad = 1.5*Scale;
				
				fprintf(OutFile,"sphere {<%12.6f, %12.6f, %12.6f> %12.6f\n",
					x,  y, -z, srad);
				
				fputs("  texture {\n",OutFile);
                int r,g,b;
            	HaColor::GetPureRGBbyCIdx(aptr->col,r,g,b);
				
				fprintf(OutFile,"    pigment {color rgb <%g, %g, %g>}\n",
					
                    (double)r/255.0, 
                    (double)g/255.0, 
                    (double)b/255.0 );
				if( FakeSpecular )
					fputs("    finish {phong 1}\n",OutFile);
				else
					fputs("    finish { AtomFinish }\n",OutFile);
				fputs("  }\n}\n",OutFile);
			}
		}
	}

	if(DrawBonds)
	{
		    HaBond* bptr;
			HaAtom  *s;
			HaAtom  *d;
			int sc,dc;
			double rt, gt, bt;
			double x1,y1,z1, x2,y2,z2;
			double rad;

		fprintf(OutFile,"// Draw Bonnds \n");
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

				x1= s->x; y1 = s->y; z1= -(s->z - zmin);
				x2= d->x; y2 = d->y; z2= -(d->z - zmin);

				
				if( bptr->flag&WireFlag )
				{      
					rad = Scale* 0.08; 
				} 
				else if( bptr->flag&CylinderFlag )
				{   
					if( bptr->irad>0 )
					{  
						rad = bptr->irad;
					} 
					else 
					{
						rad = Scale* 0.08; 
					}
				}
				
				int r,g,b;
				HaColor::GetPureRGBbyCIdx(sc,r,g,b);
				rt= r/255.0; gt= g/255.0; bt= b/255.0;

				fprintf(OutFile," cylinder { <%9.3f, %9.3f, %9.3f>, <%9.3f, %9.3f, %9.3f>, %9.3f  \n",
						x1,y1,z1, (x1+x2)/2.0 , (y1+y2)/2.0 , (z1+z2)/2.0, rad  );
				fprintf(OutFile," texture {\n pigment { color rgb<%g, %g, %g> }\n }\n}\n", rt,gt,bt );

				fprintf(OutFile,"sphere {<%12.6f, %12.6f, %12.6f> %12.6f\n", x1,  y1, z1, rad);
				fprintf(OutFile," texture {\n pigment { color rgb<%g, %g, %g> }\n }\n}\n", rt,gt,bt );
				
        
	            HaColor::GetPureRGBbyCIdx(dc,r,g,b);
				rt= r/255.0; gt= g/255.0; bt= b/255.0;

				fprintf(OutFile," cylinder { <%9.3f, %9.3f, %9.3f>, <%9.3f, %9.3f, %9.3f>, %9.3f  \n",
						(x1+x2)/2.0 , (y1+y2)/2.0 , (z1+z2)/2.0, x2,y2,z2, rad  );
				fprintf(OutFile," texture {\n pigment { color rgb<%g, %g, %g> }\n }\n}\n", rt,gt,bt );

				fprintf(OutFile,"sphere {<%12.6f, %12.6f, %12.6f> %12.6f\n", x2,  y2, z2, rad);
				fprintf(OutFile," texture {\n pigment { color rgb<%g, %g, %g> }\n }\n}\n", rt,gt,bt );

			}
		}

// Plot non-bonded atoms as a small sphere
		fprintf(OutFile,"// Draw non-bonded atoms \n");
		AtomIteratorMolSet aitr(this->GetMolSet());
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			if( aptr->GetNBonds() == 0)
			{
				double srad= 0.3*Scale;
				
				x1 = aptr->x;
				y1 = aptr->y;
				z1 = aptr->z - zmin;

				fprintf(OutFile,"sphere {<%12.6f, %12.6f, %12.6f> %12.6f\n",
					x1,  y1, -z1, srad);
				
				fputs("  texture {\n",OutFile);
	            int r,g,b;
				HaColor::GetPureRGBbyCIdx(aptr->col,r,g,b);
				
				fprintf(OutFile,"    pigment {color rgb <%g, %g, %g>}\n",
					              (double)r/255.0, (double)g/255.0, (double)b/255.0 );

				fputs("    finish {phong 1}\n",OutFile);
				fputs("  }\n}\n",OutFile);			
			}
		}
	}
	
	
	HaMolSet* pmset = GetMolSet();

	if(pmset == NULL) return False;
	ETCouplMod* ptr_et_coupl_mod = pmset->GetETCouplMod(false);
	if(ptr_et_coupl_mod != NULL  && !ptr_et_coupl_mod->best_path.empty())
	{
		int i;
		PrintLog("Plot Best Pathway in POV file\n");
		fprintf(OutFile,"// Plot the best PATHWAY \n");

		double rad = Scale* 0.2; // radius of the cylinder 
		
		for(i = 0; i < ptr_et_coupl_mod->best_path.size(); i++)
		{
			int isrc=  ptr_et_coupl_mod->best_path[i].source;
			int idest= ptr_et_coupl_mod->best_path[i].destination;
			
			if( isrc < 0 || idest < 0 ) continue; 
			
			HaAtom* aptr1= ptr_et_coupl_mod->nodes[ isrc ];
			HaAtom* aptr2= ptr_et_coupl_mod->nodes[ idest ];
			if(aptr1 != NULL && aptr2 != NULL)
			{
				double x1= aptr1->x; double y1 = aptr1->y; double z1= -(aptr1->z - zmin);
				double x2= aptr2->x; double y2 = aptr2->y; double z2= -(aptr2->z - zmin);

				if(aptr1->IsBonded(*aptr2))
				{
			    	fprintf(OutFile," cylinder { <%9.3f, %9.3f, %9.3f>, <%9.3f, %9.3f, %9.3f>, %9.3f  \n",
						x1,y1,z1, x2 , y2 , z2, rad  );
				    fprintf(OutFile," texture {\n pigment { color rgb<%g, %g, %g> }\n }\n}\n", 1.0,1.0,0.0 ); //yellow

			    	fprintf(OutFile," sphere {<%12.6f, %12.6f, %12.6f> %12.6f\n", x1,  y1, z1, rad);
				    fprintf(OutFile," texture {\n pigment { color rgb<%g, %g, %g> }\n }\n}\n", 1.0,1.0, 0.0 );

			    	fprintf(OutFile," sphere {<%12.6f, %12.6f, %12.6f> %12.6f\n", x2,  y2, z2, rad);
				    fprintf(OutFile," texture {\n pigment { color rgb<%g, %g, %g> }\n }\n}\n", 1.0,1.0, 0.0 );
				 
				}
				else
				{
   					double dist = Vec3D::CalcDistance(aptr1,aptr2,ANGSTROM_U);
				
					int ndiv = (int)(dist/0.3);
				
					if(ndiv < 3) ndiv = 3;
					
					int half_ndiv = ndiv/2;
					ndiv = half_ndiv*2 + 1;

					double xstep = (x2- x1)/ndiv;
					double ystep = (y2- y1)/ndiv;
					double zstep = (z2- z1)/ndiv;
					
					x1+= xstep; y1+= ystep; z1+= zstep;

					for(int i= 0; i < half_ndiv; i++)
					{
			    		fprintf(OutFile," cylinder { <%9.3f, %9.3f, %9.3f>, <%9.3f, %9.3f, %9.3f>, %9.3f  \n",
							x1,y1,z1, x1+xstep , y1+ystep , z1+zstep, rad  );
				        fprintf(OutFile," texture {\n pigment { color rgb<%g, %g, %g> }\n }\n}\n", 1.0,1.0,0 );

						x1+= 2* xstep; y1+= 2* ystep; z1+= 2*zstep;	
					}
				}
			}
		}	
	}
		
	
	list<Object3D*>::iterator oitr;
	for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end(); oitr++)
	{
		if( (*oitr) == NULL) continue;
		if( (*oitr)->GetObjType() == OBJ3D_SURFACE )
		{
			HaDisplayedSurface* sptr= (HaDisplayedSurface*) (*oitr);
			if(sptr == NULL)
				continue;
			
			int i;
			double x1, x2, x3, y1, y2, y3, z1, z2, z3;
			double xn1, xn2, xn3, yn1, yn2, yn3, zn1, zn2, zn3;
			int cl1, cl2, cl3;
			float dx,dy,dz;
			
			int ixadd=  pCanv->XRange()/2 ;
			int iyadd=  pCanv->YRange()/2 ;
			int izadd=  ZOffset() ;				
			
			int numverts = sptr->GetNumVerts();
			if(sptr->colors.size() < numverts) continue;
			if(numverts < 3) continue;
			
			bool subst_vert_1= true;
			
			bool have_norms = (sptr->norms.num_cols() == sptr->GetNumVerts()) && 
				              (sptr->norms.num_rows() == 3 );

			double x_zero, y_zero, z_zero; // transformed coordinates of the point (0,0,0) for 
			                               // the caclulations of the transformed norm vectors
			GetTransfCoord( 0.0, 0.0, 0.0, x_zero, y_zero, z_zero);


			int ntr = sptr->GetNumTr();
			int j;
			for(i = 1; i <= ntr; i++) 
			{
				for(j = 1; j <= 3; j++)
				{
					int idx_v =  sptr->tr_indx(j,i) +1;
					dx = sptr->verts(1,idx_v );
					dy = sptr->verts(2,idx_v );
					dz = sptr->verts(3,idx_v );
					
					double x_tr, y_tr, z_tr;
					GetTransfCoord( dx, dy, dz, x_tr, y_tr, z_tr);
					
					x3 = x_tr * Scale + ixadd;
					y3 = y_tr * Scale + iyadd;
					z3 = z_tr * Scale + izadd;
					
				    z3= -(z3- zmin);
					
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
						x1 = x3 + 0.001;
						y1 = y3;
						z1 = z3;
						xn1 = xn3;
						yn1 = yn3;
						zn1 = zn3;
						cl1 = cl3;
					}
					if( j == 2)
					{
						x2 = x3;
						y2 = y3 + 0.001;
						z2 = z3;
						xn2 = xn3;
						yn2 = yn3;
						zn2 = zn3;
						cl2 = cl3;
					}
				}
				
				if( (fabs(x1-x2)+ fabs(y1-y2)+ fabs(z1-z2)) > 0.01 && 
					(fabs(x1-x3)+ fabs(y1-y3)+ fabs(z1-z3)) > 0.01 &&
					(fabs(x2-x3)+ fabs(y2-y3)+ fabs(z2-z3)) > 0.01 )
				{
					double rt, gt, bt;
	                int r,g,b;
				    HaColor::GetPureRGBbyCIdx(cl1,r,g,b);
					rt= r/255.0;
					gt= g/255.0;
					bt= b/255.0;
					
					if(have_norms)
					{
						fprintf(OutFile," smooth_triangle { <%9.3f, %9.3f, %9.3f>, <%9.3f, %9.3f, %9.3f>, \n <%9.3f, %9.3f, %9.3f>, <%9.3f, %9.3f, %9.3f>, \n <%9.3f, %9.3f, %9.3f>, <%9.3f, %9.3f, %9.3f> \n",
							x1,y1,z1, xn1,yn1,xn1, x2,y2,z2, xn2,yn2,xn2, x3,y3,z3, xn3,yn3,xn3);	
					}
					else
					{
						fprintf(OutFile," triangle { <%9.3f, %9.3f, %9.3f>, <%9.3f, %9.3f, %9.3f>,<%9.3f, %9.3f, %9.3f> \n",
							x1,y1,z1, x2,y2,z2, x3,y3,z3);
					}
					fprintf(OutFile," texture {\n pigment { color rgb<%g, %g, %g> }\n }\n}\n", rt,gt,bt );
				}
			}
		}
	}
			
    fclose(OutFile);
    return( True );
}


void 
HaMolView::WriteVRMLTriple( double x, double y, double z )
{
#ifdef __STDC__
    fprintf(OutFile,"%.3f %.3f %.3f",x,y,z);
#else
    fprintf(OutFile,"%.3lf %.3lf %.3lf",x,y,z);
#endif
}


void 
HaMolView::WriteVRMLColour( int indent, int col_num )
{
    int i;
	  	    
	int col_idx  = ColourDepth*col_num;
    ColorVal col_val = HaColor::used_colors[col_idx];
    int r = RComp(col_val);
	int g = GComp(col_val);
    int b = BComp(col_val);

    for( i=0; i<indent; i++ )
        fputc(' ',OutFile);
    fputs("Material { diffuseColor ",OutFile);
    WriteVRMLTriple(r/255.0,g/255.0,b/255.0);
    fputs(" }\n",OutFile);
}

void
HaMolView::WriteVRMLAtoms()
{
    HaAtom  *aptr;
    double ox,oy,oz;
    double x,y,z;
    int i,flag;

    for( i=0; i< HaColor::GetNumColors(); i++ )
    {
		int col_idx = ColourDepth*i;
			flag = False;
            ox = 0.0;  oy = 0.0;  oz = 0.0;
					
			AtomIteratorMolSet aitr(GetMolSet());
		    for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
                if( (aptr->IsDrawSphere()) && (col_idx==  aptr->col) )
                {   
					if( !flag )
                    {   
						WriteVRMLColour(2,i);
                        fputs("  Separator {\n",OutFile);
                        flag = True;
                    }

                    x = aptr->GetX();
                    y = aptr->GetY();
                    z = aptr->GetZ();

                    fputs("    Translation { translation ",OutFile);
                    WriteVRMLTriple(x-ox,y-oy,z-oz);
                    ox = x;  oy = y;  oz = z;
                    fputs(" }\n",OutFile);

                    fputs("    Sphere { radius ",OutFile);
                    fprintf(OutFile,"%.3f }\n", aptr->image_radius);
                }

            if( flag )
                fputs("  }\n",OutFile);
        }
}
                            
void 
HaMolView::WriteVRMLLine( int src, int dst, int shade, int* flag )
{
    if( !*flag )
    {   WriteVRMLColour(4,shade);
        fputs("    IndexedLineSet {\n",OutFile);
        fputs("      coordIndex [\n",OutFile);
        *flag = True;
    }
    fprintf(OutFile,"        %5d, %5d, -1,\n",src-1,dst-1);
}

void HaMolView::WriteVRMLWireframe()
{
    HaAtom  *aptr;
    HaBond  *bptr;
    HaAtom  *src;
    HaAtom  *dst;
    double x,y,z;
    int i,j;
    static int flag;

	map<HaAtom*, short, less<HaAtom*> > at_mbox; // Shadow?

	AtomIteratorMolSet aitr(GetMolSet());
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
        at_mbox[aptr] = 0;

    i = 0;
	BondIteratorMolSet bitr(GetMolSet());
    for(bptr=bitr.GetFirstBond();bptr; bptr= bitr.GetNextBond() )
        if( bptr->flag & WireFlag )
        {   src = bptr->srcatom;
            dst = bptr->dstatom;
            if( i == 0 )
            {   fputs("  Separator {\n",OutFile);
                fputs("    Coordinate3 {\n",OutFile);
                fputs("      point [\n",OutFile);
            }

            if( !at_mbox[src])
            {   
				x = src->GetX();
                y = src->GetY();
                z = src->GetZ();

                fputs("        ",OutFile);
                WriteVRMLTriple(x,y,z);
                fputs(",\n",OutFile);
                at_mbox[src] = ++i;
            }
            
            if( !at_mbox[dst] )
            {   x = (double)dst->GetX();
                y = (double)dst->GetY();
                z = (double)dst->GetZ();

                fputs("        ",OutFile);
                WriteVRMLTriple(x,y,z);
                fputs(",\n",OutFile);
                at_mbox[dst] = ++i;
            }

            if( !bptr->col && (src->col!=dst->col) )
            {   
				x = (double)(src->GetX()+dst->GetX())/2.0;
                y = (double)(src->GetY()+dst->GetY())/2.0;
                z = (double)(src->GetZ()+dst->GetZ())/2.0;
                
                fputs("        ",OutFile);
                WriteVRMLTriple(x,y,z);
                fputs(",\n",OutFile);
                i++;
            }
        }

    /* No wireframe! */
    if( !i )  return;

    fputs("      ]\n",OutFile);
    fputs("    }\n",OutFile);
    
	int nc = HaColor::GetNumColors();
    for( j=0; j< nc; j++ )
        {   
			i = 1;
			int col_idx = j*ColourDepth;
            flag = False;
            for(bptr=bitr.GetFirstBond();bptr; bptr= bitr.GetNextBond() )
                if( bptr->flag & WireFlag )
                {   
					src = bptr->srcatom;   
                    dst = bptr->dstatom;

                    if( at_mbox[src] == i ) i++;
                    if( at_mbox[dst] == i ) i++;

                    if( bptr->col )
                    {   
						if( bptr->col == col_idx )
                            WriteVRMLLine(at_mbox[src],at_mbox[dst],j,&flag);
                    } 
					else if( src->col == dst->col )
                    {   
						if( src->col == col_idx )
                            WriteVRMLLine(at_mbox[src],at_mbox[dst],j,&flag);
                    } 
					else /* Two Colour HaBond */
                    {   
						if( src->col == col_idx )
                        {   
							WriteVRMLLine(at_mbox[src],i,j,&flag);
                        } 
						else if( dst->col == col_idx )
                            WriteVRMLLine(at_mbox[dst],i,j,&flag);
                        i++;
                    }
                }

            if( flag )
            {   fputs("      ]\n",OutFile);
                fputs("    }\n",OutFile);
            }
        }
    fputs("  }\n",OutFile);
}


void 
HaMolView::WriteVRMLDots()
{
    vector<int> hist;
    DotStruct  *ptr;
    double x,y,z;
    int count;
    int flag;
    int i,j;
    flag = False;

	int nc = HaColor::GetNumColors();

    for( i=0; i< nc; i++ )
	{
		int col_idx = ColourDepth*i;
		
		count = 0;
		HaMolSet* pmset = GetMolSet();
		
		list<Object3D*>::iterator oitr;
		
		for(oitr = pmset->ViewObjects.begin(); oitr != pmset->ViewObjects.end(); oitr++)
		{
			if( (*oitr)->GetObjType() != OBJ3D_DOT_SURFACE )
				continue;
			
			ptr = (DotStruct*) *oitr;
			int np = ptr->GetCount();
			
			for( j=0; j< np; j++ )
			{
				if( col_idx == i )
				{   
					if( !flag )
					{   
						fputs("  Separator {\n",OutFile);
						fputs("    Coordinate3 {\n",OutFile);
						fputs("      point [\n",OutFile);
						flag = True;
					}
					
					x = (double)ptr->dots[j].GetX_Ang();
					y = (double)ptr->dots[j].GetY_Ang();
					z = (double)ptr->dots[j].GetZ_Ang();
					
					fputs("        ",OutFile);
					WriteVRMLTriple(x,y,z);
					
					fputs(",\n",OutFile);
					count++;					
					hist[i] = count;
				} 
				else 
					hist[i] = 0;
			}

			if( flag )
			{   
				fputs("      ]\n",OutFile);
				fputs("    }\n",OutFile);

				count = 0;
				for( i=0; i< nc; i++ )
				{
					if( hist[i] )
					{   
						WriteVRMLColour(4,i);

						fputs("    PointSet {\n      ",OutFile);
						fprintf(OutFile,"startIndex %d numPoints %d\n",count,hist[i]);
						fputs("    }\n",OutFile);
						count += hist[i];
					}
				}
				fputs("  }\n",OutFile);
			}
		}
	}
}

int 
HaMolView::WriteVRMLFile( const char* name )
{
    OutFile = fopen(name,"w");
    if( !OutFile )
    {   
		PrintLog("Script Error: Unable to create file %s \n",name);
        return(False);
    }

	HaMolSet* pmset = GetMolSet();

    fputs("#VRML V1.0 ascii\n",OutFile);
    fputs("#Created by HARLEM v0.2\n\n",OutFile);

    fputs("DEF Viewer Info { string \"examiner\" }\n",OutFile);
    fprintf(OutFile,"DEF Title Info { string \"%s\" }\n",pmset->GetName());
    fputs("DEF Creator Info { string \"Created by HARLEM v0.1\" }\n",OutFile);
    fputs("DEF Background Info { string \"",OutFile);
    WriteVRMLTriple(BackColor.r/255.0, BackColor.g/255.0, BackColor.b/255.0);
    fputs("\" }\n\n",OutFile);
    
    fputs("Separator {\n",OutFile);

    WriteVRMLAtoms();
    WriteVRMLWireframe();
    if( DrawDots )
        WriteVRMLDots();

    fputs("}\n",OutFile);
    fclose(OutFile);
    return( True );
}


