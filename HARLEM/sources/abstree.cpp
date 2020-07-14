/*! \file abstree.cpp 

   classes for logical expression evaluation
 
   derived from:    
   
   abstree.c    
   RasMol2 Molecular Graphics   
   Roger Sayle, August 1995  
   Version 2.6   

   \author Igor Kurnikov
   \date 1997-2002  
    
*/
#include "haconst.h"

#ifdef _WIN32
#include <malloc.h>
#endif
#ifndef sun386
#include <stdlib.h>
#endif
#include <string>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#include <stdexcept>

#define ABSTREE_CPP

#include "haio.h"
#include "command.h"
#include "hamolset.h"
#include "hamolecule.h"
#include "abstree.h"

HaAtom*    AtomExpr::QAtom  =NULL;
HaResidue* AtomExpr::QGroup =NULL;
HaChain* AtomExpr::QChain =NULL;

//static AtomExpr FalseExpr( (OpConst | OpLftVal | OpRgtVal), 0, 0 );
//static AtomExpr TrueExpr( (OpConst | OpLftVal | OpRgtVal), 1, 1 );


/* Macros for commonly used loops */

#define BitAcidic       0x001
#define BitAliphatic    0x002
#define BitAromatic     0x004
#define BitBasic        0x008
#define BitBuried       0x010
#define BitCyclic       0x020
#define BitHydrophobic  0x040
#define BitMedium       0x080
#define BitNeutral      0x100
#define BitSmall        0x200

#define BitCharged      0x009
#define BitNotLarge     0x280

/* Acyclic = !Cyclic         */
/* Large = !Medium && !Small */
/* Polar = !Hydrophobic      */
/* Surface = !Buried         */


static int AminoProp[] = {
        /*ALA*/  BitAliphatic | BitBuried | BitHydrophobic | BitNeutral |
                 BitSmall,
        /*GLY*/  BitAliphatic | BitHydrophobic | BitNeutral | BitSmall,
        /*LEU*/  BitAliphatic | BitBuried | BitHydrophobic | BitNeutral,
        /*SER*/  BitNeutral | BitSmall,
        /*VAL*/  BitAliphatic | BitBuried | BitHydrophobic | BitMedium |
                 BitNeutral,
        /*THR*/  BitMedium | BitNeutral,
        /*LYS*/  BitBasic,
        /*ASP*/  BitAcidic | BitMedium,
        /*ILE*/  BitAliphatic | BitBuried | BitHydrophobic | BitNeutral,
        /*ASN*/  BitMedium | BitNeutral,
        /*GLU*/  BitAcidic,
        /*PRO*/  BitCyclic | BitHydrophobic | BitMedium | BitNeutral,
        /*ARG*/  BitBasic,
        /*PHE*/  BitAromatic | BitBuried | BitCyclic | BitHydrophobic |
                 BitNeutral,
        /*GLN*/  BitNeutral,
        /*TYR*/  BitAromatic | BitCyclic | BitNeutral,
        /*HIS*/  BitAromatic | BitBasic | BitCyclic | BitNeutral,
        /*CYS*/  BitBuried | BitMedium | BitNeutral,
        /*MET*/  BitBuried | BitHydrophobic | BitNeutral,
        /*TRP*/  BitAromatic | BitBuried | BitCyclic | BitHydrophobic |
                 BitNeutral,

        /*ASX*/  BitMedium | BitNeutral,
        /*GLX*/  BitNeutral,
        /*PCA*/  BitCyclic | BitHydrophobic | BitMedium | BitNeutral,
        /*HYP*/  BitCyclic | BitMedium | BitNeutral
        };


AtomExpr::AtomExpr(int new_type, long rval, long lval)
{
	type=new_type;
	rgt.val=rval;
	lft.val=lval;
}

AtomExpr::AtomExpr()
{
	type = 0;
    rgt.ptr = NULL;
    lft.ptr = NULL;
}

AtomExpr::~AtomExpr()
{
    if( type!=OpWithin )
    {
		if( !( type&(OpLftProp|OpLftVal)) )
            delete lft.ptr;
        if( !( type&(OpRgtProp|OpRgtVal)) )
            delete rgt.ptr;
    }
	else
	{
		delete rgt.set;
	}
}

AtomExpr* AtomExpr::CreateTrueExpr()
{
	return new AtomExpr( (OpConst | OpLftVal | OpRgtVal), 1, 1 );
}

AtomExpr* AtomExpr::CreateFalseExpr()
{
	return new AtomExpr( (OpConst | OpLftVal | OpRgtVal), 0, 0 );
}

bool AtomExpr::IsTrueExpr(AtomExpr* p_expr)
{
	int type_t =  (OpConst | OpLftVal | OpRgtVal);
	if( p_expr->type == type_t && p_expr->lft.val == 1 && p_expr->rgt.val == 1) return true;
	return false;
}
   
bool AtomExpr::IsFalseExpr(AtomExpr* p_expr)
{
	int type_t =  (OpConst | OpLftVal | OpRgtVal);
	if( p_expr->type == type_t && p_expr->lft.val == 0 && p_expr->rgt.val == 0) return true;
	return false;
}

long AtomExpr::EvaluateProperty(long prop )
// Here to include new properties
{
    switch( prop )
    {
	case( PropIdent ):    return( (long) QAtom->GetSerNo() );
	case( PropXCord ):    return( (long) (QAtom->GetX_Ang()*1000.0) );
	case( PropYCord ):    return( (long) (QAtom->GetY_Ang()*1000.0) );
	case( PropZCord ):    return( (long) (QAtom->GetZ_Ang()*1000.0) );
	case( PropTemp ):     return( (long) (QAtom->tempf*1000.0) );
	case( PropName ):     return( QAtom->refno );
	case( PropResId ):    return( QGroup->serno );
	case( PropResName ):  return( QGroup->refno );
	case( PropChain ):    return( (long) QChain->ident );
	case( PropSelect ):   return( QAtom->Selected() );
	case( PropElemNo ):   return( QAtom->GetElemNo() );
	case( PropChemGroup ):
		if( QAtom->GetHostChemGroup() != NULL)
			return( atoi( (QAtom->GetHostChemGroup())->GetID() ) );
	case( PropMolPtr ):   return( (long)QAtom->GetHostMol());
	case( PropRad ):      if( QAtom->IsDrawSphere() )
						  {
							  return( (long)( QAtom->radius*1000.0 ) );
						  }
		else
			return( 0 );
		

        /* Predicates stored in flags */
	case( PredBonded ):       return( !(QAtom->flag&NonBondFlag) );
	case( PredHydrogen ):     return( QAtom->IsHydrogen() );
	case( PredHetero ):       return( QAtom->flag&HeteroFlag );
	case( PredCystine ):      return( QGroup->flag&CystineFlag );
	case( PredHelix ):        return( QGroup->struc&HelixFlag );
	case( PredSheet ):        return( QGroup->struc&SheetFlag );
	case( PredTurn ):         return( QGroup->struc&TurnFlag );
		
        /* Residue type predicates */
	case( PredDNA ):          return( QGroup->IsDNA() );
	case( PredRNA ):          return( QGroup->IsRNA() );
	case( PredNucleic ):      return( QGroup->IsNucleo() );
	case( PredProtein ):      return( QGroup->IsProtein() );
	case( PredAmino ):        return( QGroup->IsAmino() );
	case( PredWater ):        return( QGroup->IsWater() );
	case( PredSolvent ):      return( QGroup->IsSolvent() );
	case( PredIon ):          return( QGroup->IsIon() );
		
        /* General Predicates */
	case( PredAlpha ):        return( QGroup->IsAmino() &&
								  QAtom->IsAlphaCarbon() );
	case( PredMainChain ):    return( (QGroup->IsAmino() &&
								  QAtom->IsAminoBackbone() ) ||
								  (QGroup->IsNucleo() &&
								  QAtom->IsNucleicBackbone() ) );
	case( PredSidechain ):    return( QGroup->IsAmino() &&
								  ! QAtom->IsAminoBackbone() );
	case( PredLigand ):       return( (QAtom->flag&HeteroFlag) &&
								  ! QGroup->IsSolvent() );
		
        /* Nucleic Acid Classifications */
	case( PredAT ):           return( QGroup->IsAdenine() ||
								  QGroup->IsThymine() );
	case( PredCG ):           return( QGroup->IsCytosine() ||
								  QGroup->IsGuanine() );
	case( PredPyrimidine ):   return( QGroup->IsPyrimidine() );
	case( PredPurine ):       return( QGroup->IsPurine() );
		
		
        /* Amino Acid Classifications */
	case( PredAcidic ):       return( QGroup->IsAmino() &&
                                  AminoProp[QGroup->refno]&BitAcidic );
		
	case( PredAcyclic ):      return( QGroup->IsAmino() &&
                                  !(AminoProp[QGroup->refno]&BitCyclic) );
		
	case( PredAliphatic ):    return( QGroup->IsAmino() &&
                                  AminoProp[QGroup->refno]&BitAliphatic );
		
	case( PredAromatic ):     return( QGroup->IsAmino() &&
                                  AminoProp[QGroup->refno]&BitAromatic );
		
	case( PredBasic ):        return( QGroup->IsAmino() &&
                                  AminoProp[QGroup->refno]&BitBasic );
		
	case( PredBuried ):       return( QGroup->IsAmino() &&
                                  AminoProp[QGroup->refno]&BitBuried );
		
	case( PredCharged ):      return( QGroup->IsAmino() &&
                                  AminoProp[QGroup->refno]&BitCharged );
		
	case( PredCyclic ):       return( QGroup->IsAmino() &&
                                  AminoProp[QGroup->refno]&BitCyclic );
		
	case( PredHydrophobic ):  return( QGroup->IsAmino() &&
                                  AminoProp[QGroup->refno]&BitHydrophobic );
		
	case( PredLarge ):        return( QGroup->IsAmino() &&
                                  !(AminoProp[QGroup->refno]&BitNotLarge) );
		
	case( PredMedium ):       return( QGroup->IsAmino() &&
                                  AminoProp[QGroup->refno]&BitMedium );
		
	case( PredNeutral ):      return( QGroup->IsAmino() &&
                                  AminoProp[QGroup->refno]&BitNeutral );
		
	case( PredPolar ):        return( QGroup->IsAmino() &&
                                  !(AminoProp[QGroup->refno]&BitHydrophobic) );
		
	case( PredSmall ):        return( QGroup->IsAmino() &&
                                  AminoProp[QGroup->refno]&BitSmall );
		
	case( PredSurface ):      return( QGroup->IsAmino() &&
                                  !(AminoProp[QGroup->refno]&BitBuried) );
		
    }
    return( True );
}


long AtomExpr::EvaluateExprFor(HaAtom* aptr)
{
	QAtom=aptr;
    QChain=QAtom->GetHostChain();
    QGroup=QAtom->GetHostRes();
	return(EvaluateExpr());
	
}


long AtomExpr::EvaluateExpr()
// Logical, recursive evaluation of the logical expression
{
    long ilft, irgt;

// Within Expression Type

    if( type==OpWithin )
    {
		if( lft.dval > 0.01 )
        {
			return( rgt.set->IsWithinRadius(QAtom, lft.dval) );
        }
		else
			return( rgt.set->IsMember(QAtom) );
    }
	else if( type==OpMember )
        return( rgt.set->IsMember(QAtom) );

//

    if( type & OpLftVal )
    {
		ilft = lft.val;
    }
	else if( type & OpLftProp )
    {
		ilft = EvaluateProperty( lft.val );
    }
	else
		ilft = lft.ptr->EvaluateExpr(); // resursively evaluate left hand side of the expression

    if( OpCode()==OpConst ) return( ilft );
    if( (OpCode()==OpAnd) && !ilft ) return( False );
    if( (OpCode()==OpOr) && ilft ) return( True );
    if( OpCode()==OpNot ) return( !ilft );

// Evaluation of the rigth side of the expression

    if( type & OpRgtVal )
    {
		irgt = rgt.val;
    }
	else if( type & OpRgtProp )
    {
		irgt = EvaluateProperty( rgt.val );
    }
	else
		irgt = rgt.ptr->EvaluateExpr(); // resursively evaluate right hand side  of the expression

    switch( OpCode() )
    {
	case(OpOr):
        case(OpAnd):     return( irgt );
        case(OpLess):    return( ilft < irgt );
        case(OpMore):    return( ilft > irgt );
        case(OpEqual):   return( ilft == irgt );
        case(OpNotEq):   return( ilft != irgt );
        case(OpLessEq):  return( ilft <= irgt );
        case(OpMoreEq):  return( ilft >= irgt );
    }
    return( True );
}

AtomExpr* AtomExpr::LookUpAtGroupExpr(const char* grp_name, MolSet* pmset )
{
    AtomGroup* atl;
	atl = pmset->GetAtomGroupByID(grp_name);
	if(atl == NULL) 
		return NULL;
	
	AtomExpr* expr = new AtomExpr();
    expr->type = OpMember;
    expr->rgt.set = atl;

    return( expr );
}


static int ElemCompare(const char* ident,const char* elem )
{
    while( *elem )
        if( *elem++ != *ident++ )
            return( False );

    /* Handle Plurals */
    if( (ident[0]=='S') && !ident[1] )
        return( (elem[-1]!='S') && (elem[-1]!='Y') );
    return( !*ident );
}


AtomExpr* AtomExpr::LookUpElement(const char* ident )
{
    AtomExpr *expr;
    int elem;

    for( elem=1; elem<MAXELEMNO; elem++ )
        if( ElemCompare(ident,ElementArr[elem].name) )
            break;

    /* Handle Difficult Plurals & US Spelling! */
    if( elem == MAXELEMNO )
    {
		if( *ident=='A' )
        {
			if( ElemCompare(ident,"ALUMINUM") )
            {
				elem = 13;
            }
			else if( !strcmp(ident,"ANTIMONIES") )
                elem = 51;
        }
		else if( *ident=='C' )
        {
			if( ElemCompare(ident,"CESIUM") )
                elem = 55;
        }
		else if( *ident=='M' )
        {
			if( !strcmp(ident,"MERCURIES") )
                elem = 80;
        }
		else if( *ident=='P' )
        {
			if( !strcmp(ident,"PHOSPHORUSES") )
                elem = 8;
        }
		else if( *ident=='S' )
        {
			if( ElemCompare(ident,"SULFUR") )
                elem = 16;
        }
    }

    if( elem<MAXELEMNO )
    {
		expr = new AtomExpr();
        expr->type = OpEqual|OpLftProp|OpRgtVal;
        expr->lft.val = PropElemNo;
        expr->rgt.val = elem;
    }
	else
		expr = (AtomExpr*)0;
    return( expr );
}



static int MatchWildName(const char* src, const char* dst, int size, int len )
{
    int i, left;

    left = size;
    while( *dst==' ' )
    {
		dst++; left--;
    }

    for( i=0; i<len; i++ )
    {
		if( left )
        {
			if( (*dst==*src) || (*src=='?') )
            {
				dst++;  src++;  left--;
            }
			else
				return( False );
        }
		else if( *src++ != '?' )
            return( False );
    }

    while( left )
	{
         if( *dst++!=' ' )
         {
			 return( False );
         }
		 else
			 left--;
	}
    return( True );
}


AtomExpr* AtomExpr::ParsePrimitiveExpr(MolSet* pmset, const char* expr_str, size_t& cr_pos )
//! 
//! Generate a logical expression of subset of atoms of pmset corresponding to a string token 
//! \param pmset    - Molecular Set atom expression is formed on
//! \param expr_str& - string with expression
//! \param cr_pos&   - current cursor position (input and output)
{
	std::string NameBuf;
    AtomExpr *tmp1,*tmp2;
    AtomExpr *wild;
    
    int i, j;
    int neg;
    int ch;

	AtomExpr* p_expr = CreateTrueExpr();

	if(cr_pos >= strlen(expr_str)) return NULL;
 
	try
	{
		ch = expr_str[cr_pos++];	
		if(ch == '$') // Process Molecule Name
		{
			NameBuf.clear();
			while( (ch = expr_str[cr_pos++]) != '$')
			{
				if( ch != 0 )
				{
					NameBuf += toupper(ch);
				}
				else
				{
					throw std::runtime_error("No closing $ in molecule name ");
				}
			}
			if( NameBuf.size() == 0 )
			{
				throw std::runtime_error("Empty molecule name ");
			}

			wild = CreateFalseExpr();
			int nmol= pmset->GetNMol();
			for( i=0; i < nmol; i++)
			{
				HaMolecule* pMol= pmset->HostMolecules[i];
				if( NameBuf == pMol->GetObjName() )
				{
					tmp1 = new AtomExpr();
					tmp1->type = OpEqual | OpLftProp | OpRgtVal;
					tmp1->lft.val = PropMolPtr;
					tmp1->rgt.val = (long)pMol;
				
					tmp2 = new AtomExpr();
					tmp2->type = OpOr;
					tmp2->lft.ptr = tmp1;
					tmp2->rgt.ptr = wild;
					wild = tmp2;	
				}
			}
			p_expr = wild;
			ch = expr_str[cr_pos++];
		} // end of Molecule name parsing

		if( !ch || isspace(ch))
		{
			cr_pos--;
			return p_expr;
		}

		if( ch != ':' && ch != 0) // Parse Residue Name 
		{   
			NameBuf.clear();
			if( ch != '*' )
			{
				if( ch == '[' ) //Parse residue name in [ ] brackets
				{
					while( (ch = expr_str[cr_pos++]) != ']' )
					{
						if( ch != 0 )
						{
							NameBuf += toupper(ch);
						}
						else
						{
							throw std::runtime_error("No closing ] in Residue Name");
						}
					}
				}
				else // Parse residue name without brackets
				{
					for( i=0; i< 20; i++ )
					{
						if( isalpha(ch) )
						{
							NameBuf += toupper(ch);
							ch = expr_str[cr_pos++];
						}
						else if( (ch == '?') || (ch == '%') )
						{
							NameBuf += '?';
							ch = expr_str[cr_pos++];
						}
						else
						{
							break;
						}
					}
				} 
				if( NameBuf.empty() ) throw std::runtime_error(" Empty Residue name ");

				wild = CreateFalseExpr();
				for( j=0; j < HaResidue::ResNames.size(); j++ )
				{
					if( MatchWildName(NameBuf.c_str(),HaResidue::ResNames[j].c_str(), HaResidue::ResNames[j].size(),NameBuf.size() ) )
					{
						tmp1 = new AtomExpr();
						tmp1->type = OpEqual | OpLftProp | OpRgtVal;
						tmp1->lft.val = PropResName;
						tmp1->rgt.val = j;
					
						tmp2 = new AtomExpr();
						tmp2->type = OpOr;
						tmp2->lft.ptr = tmp1;
						tmp2->rgt.ptr = wild;
						wild = tmp2;
					}
				}
			
				if( !IsTrueExpr(p_expr) )
				{
					tmp2 = new AtomExpr();
					tmp2->type = OpAnd;
					tmp2->lft.ptr = wild;
					tmp2->rgt.ptr = p_expr;
					wild = tmp2;
				}
				p_expr = wild;
			}
			else // ch != '*' - for residue name
			{
				ch = expr_str[cr_pos++];
			}

			if( ch != '*' && ch != 0) // Parse Residue Number 
			{
				if( ch == '-' )
				{
					ch = expr_str[cr_pos++];
					neg = True;
				}
				else
				{
					neg = False;
				}

				if( isdigit(ch) )
				{
					i = ch-'0';
					while( isdigit(expr_str[cr_pos]) )
					{
						i = 10*i + (expr_str[cr_pos]-'0');
						cr_pos++;
					}

					tmp1 = new AtomExpr();
					tmp1->type = OpEqual | OpLftProp | OpRgtVal;
					tmp1->rgt.val = neg? -i : i;
					tmp1->lft.val = PropResId;
					if( !IsTrueExpr(p_expr) )
					{
						tmp2 = new AtomExpr();
						tmp2->type = OpAnd;
						tmp2->rgt.ptr = p_expr;
						tmp2->lft.ptr = tmp1;
						p_expr = tmp2;
					}
					else
					{
						p_expr = tmp1;
					}
					ch = expr_str[cr_pos++];
				}
				else if( neg )
				{
					throw std::runtime_error("Only minus sign in place of residue number");
				}
			}
			else if( ch != 0 )  // if (ch == '*') in place of residue number
			{
				ch = expr_str[cr_pos++];
			}
		}	 

		if( ch==':' ) // Parse Chain Ident 
		{
			ch = expr_str[cr_pos++];
		}

		if( isalnum(ch) )
		{
			ch = toupper(ch);

			tmp1 = new AtomExpr();
			tmp1->type = OpEqual | OpLftProp | OpRgtVal;
			tmp1->lft.val = PropChain;
			tmp1->rgt.val = (long) ch;
			if( !IsTrueExpr(p_expr) )
			{
				tmp2 = new AtomExpr();
				tmp2->type = OpAnd;
				tmp2->rgt.ptr = p_expr;
				tmp2->lft.ptr = tmp1;
				p_expr = tmp2;
			}
			else
			{
				p_expr = tmp1;
			}
			ch = expr_str[cr_pos++];
		}
		else if( (ch=='?') || (ch=='%') || (ch=='*') )
		{
			ch = expr_str[cr_pos++];
		}

		if( ch == ':' ) // Parse Model Number 
		{
			ch = expr_str[cr_pos++];
			if( isdigit(ch) )
			{
				i = ch-'0';
				while( isdigit(expr_str[cr_pos]) )
				{
					i = 10*i + (expr_str[cr_pos]-'0');
					cr_pos++;
				}

				tmp1 = new AtomExpr();
				tmp1->type = OpEqual | OpLftProp | OpRgtVal;
				tmp1->lft.val = PropModel;
				tmp1->rgt.val = i;
				if( IsTrueExpr(p_expr) )
				{
					tmp2 = new AtomExpr();
					tmp2->type = OpAnd;
					tmp2->rgt.ptr = p_expr;
					tmp2->lft.ptr = tmp1;
					p_expr = tmp2;
				}
				else
				{
					p_expr = tmp1;
				}
				ch = expr_str[cr_pos++];
			}
			else
			{
				throw std::runtime_error("No model number after : ");
			}
		} // End Parse Model Number

		if( ch == '.' ) // Parse Atom Name 
		{
			NameBuf.clear();
			ch = expr_str[cr_pos++];
			if( ch != '*' )
			{
				for( i=0; i< 20; i++ )
				{
					if( isalnum(ch) || ch=='\'' || ch=='*' )
					{
						NameBuf += toupper(ch);
						ch = expr_str[cr_pos++];
					}
					else if( (ch=='?') || (ch=='%') || (ch=='#') )
					{
						NameBuf += '?';
						ch = expr_str[cr_pos++];
					}
					else
					{
						break;
					}
				}
				if( NameBuf.empty() ) throw std::runtime_error("Empty Atom Name");  
			
//				wild = &FalseExpr;
				wild = CreateFalseExpr();
				for( j=0; j< HaAtom::ElemDesc.size(); j++ )
				{
					if( MatchWildName(NameBuf.c_str(), HaAtom::ElemDesc[j].c_str(), HaAtom::ElemDesc[j].size(),NameBuf.size()) )
					{
						tmp1 = new AtomExpr();
						tmp1->type = OpEqual | OpLftProp | OpRgtVal;
						tmp1->lft.val = PropName;
						tmp1->rgt.val = j;

						tmp2 = new AtomExpr();
						tmp2->type = OpOr;
						tmp2->lft.ptr = tmp1;
						tmp2->rgt.ptr = wild;

						wild = tmp2;
					}
				}

				if( IsTrueExpr(p_expr) || IsFalseExpr(wild) )
				{
					delete p_expr;
					p_expr = wild;
				}
				else
				{
					tmp1 = new AtomExpr();
					tmp1->type = OpAnd;
					tmp1->lft.ptr = p_expr;
					tmp1->rgt.ptr = wild;
					p_expr = tmp1;
				}
			}
			else // if( ch == '*' )
			{
				ch = expr_str[cr_pos++];
			}
		} // if(ch == '.') End Parse Atom Name
	}
	catch( std::exception& ex)
	{
		PrintLog(" Parsing Error %s \n", ex.what());
		cr_pos--;
		delete p_expr;
		return NULL;
	}
    cr_pos--;
    if( !ch || isspace(ch) || ispunct(ch) )
	{
		return p_expr;
	}
	delete p_expr;
	return NULL;
} 

AtomExpr* AtomExpr::ParseExpression(const std::string& expr_str, MolSet* pmset)
{
	CmdParser cmd_parser;
	cmd_parser.SetCmdLine(expr_str);
	cmd_parser.FetchToken();
	AtomExpr* p_expr = cmd_parser.ParseExpression(0,pmset);
	if( cmd_parser.CurToken )
	{   
		PrintLog("Error in AtomExpr::ParseExpression() \n");
		PrintLog("Invalid Expression String syntax \n");
		if(p_expr) delete p_expr;
		p_expr = NULL;
	} 
	return p_expr; 
}
