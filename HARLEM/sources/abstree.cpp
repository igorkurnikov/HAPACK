/*! \file abstree.cpp 

   classes for logical expression evaluation
 
   derived from:    
   
   abstree.c    
   RasMol2 Molecular Graphics   
   Roger Sayle, August 1995  
   Version 2.6   

   \author Igor Kurnikov
   \date 1997-2025  
    
*/
#include "haconst.h"

#ifdef _WIN32
#include <malloc.h>
#endif
#include <stdlib.h>
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


AtomExpr::AtomExpr(int new_type, AtomExprVal rval, AtomExprVal lval)
{
	type=new_type;
	rgt=rval;
	lft=lval;
}

AtomExpr::AtomExpr()
{
	type = (OpConst | OpLftVal | OpRgtVal);
	rgt = false;
	lft = false;
}

AtomExpr::~AtomExpr()
{

}

std::shared_ptr<AtomExpr> AtomExpr::CreateTrueExpr()
{
	std::shared_ptr<AtomExpr> pexpr = std::make_shared<AtomExpr>((OpConst | OpLftVal | OpRgtVal), true, true);
	return pexpr;
}

std::shared_ptr<AtomExpr> AtomExpr::CreateFalseExpr()
{
	std::shared_ptr<AtomExpr> pexpr = std::make_shared<AtomExpr>((OpConst | OpLftVal | OpRgtVal), false, false);
	return pexpr;
}

bool AtomExpr::IsTrueExpr(AtomExpr* p_expr)
{
	int type_t =  (OpConst | OpLftVal | OpRgtVal);
	if( p_expr->type == type_t 
		&& std::holds_alternative<bool>(p_expr->lft) && std::get<bool>(p_expr->lft) == true 
		&& std::holds_alternative<bool>(p_expr->rgt) && std::get<bool>(p_expr->rgt) == true) return true;
	return false;
}
   
bool AtomExpr::IsFalseExpr(AtomExpr* p_expr)
{
	int type_t =  (OpConst | OpLftVal | OpRgtVal);
	if( p_expr->type == type_t  
		&& std::holds_alternative<bool>(p_expr->lft) && std::get<bool>(p_expr->lft) == false
		&& std::holds_alternative<bool>(p_expr->rgt) && std::get<bool>(p_expr->rgt) == false ) return true;
	return false;
}

AtomExprVal AtomExpr::EvaluatePropertyFor(HaAtom* aptr, long prop)
{
	HaResidue* pres = aptr->GetHostRes();
	HaChain* pchain = pres->GetHostChain();

	switch (prop)
	{
	case(PropIdent):    return( aptr->GetSerNo() );
	case(PropXCord):    return( aptr->GetX_Ang() );
	case(PropYCord):    return( aptr->GetY_Ang() );
	case(PropZCord):    return( aptr->GetZ_Ang() );
	case(PropTemp):     return( aptr->tempf );
	case(PropName):     return( aptr->refno);
	case(PropResId):    return( pres->serno);
	case(PropResName):  return( pres->refno);
	case(PropChain):    return( pchain->ident );
	case(PropSelect):   return( aptr->Selected());
	case(PropElemNo):   return( aptr->GetElemNo());
	case(PropChemGroup):
		if (aptr->GetHostChemGroup() != NULL)
		{
			std::string chem_grp_id = aptr->GetHostChemGroup()->GetID();
			return chem_grp_id;
		}
	case(PropMolPtr):   return( aptr->GetHostMol() );
	case(PropRad):      
		if ( aptr->IsDrawSphere() )
		{
			return(( aptr->radius));
		}
		else
			return( (double) 0.0);


		/* Predicates stored in flags */
	case(PredBonded):       return(!(aptr->flag & NonBondFlag));
	case(PredHydrogen):     return( aptr->IsHydrogen());
	case(PredHetero):       return( (bool) (aptr->flag & HeteroFlag) );
	case(PredCystine):      return( (bool) (pres->flag & CystineFlag) );
	case(PredHelix):        return( (bool) (pres->struc & HelixFlag) );
	case(PredSheet):        return( (bool) (pres->struc & SheetFlag) );
	case(PredTurn):         return( (bool) (pres->struc & TurnFlag) );

		/* Residue type predicates */
	case(PredDNA):          return(pres->IsDNA());
	case(PredRNA):          return(pres->IsRNA());
	case(PredNucleic):      return(pres->IsNucleo());
	case(PredProtein):      return(pres->IsProtein());
	case(PredAmino):        return(pres->IsAmino());
	case(PredWater):        return(pres->IsWater());
	case(PredSolvent):      return(pres->IsSolvent());
	case(PredIon):          return(pres->IsIon());

		/* General Predicates */
	case(PredAlpha):        return(pres->IsAmino() && aptr->IsAlphaCarbon());
	case(PredMainChain):    return (pres->IsAmino() && aptr->IsAminoBackbone()) ||
		                           (pres->IsNucleo() && aptr->IsNucleicBackbone());
	case(PredSidechain):    return(pres->IsAmino() && !aptr->IsAminoBackbone());
	case(PredLigand):       return((aptr->flag & HeteroFlag) && !pres->IsSolvent());

		/* Nucleic Acid Classifications */
	case(PredAT):           return(pres->IsAdenine()  || pres->IsThymine());
	case(PredCG):           return(pres->IsCytosine() || pres->IsGuanine());
	case(PredPyrimidine):   return(pres->IsPyrimidine());
	case(PredPurine):       return(pres->IsPurine());

		/* Amino Acid Classifications */
	case(PredAcidic):       return(pres->IsAmino() && AminoProp[pres->refno] & BitAcidic);
	case(PredAcyclic):      return(pres->IsAmino() && !(AminoProp[pres->refno] & BitCyclic));
	case(PredAliphatic):    return(pres->IsAmino() && AminoProp[pres->refno] & BitAliphatic);
	case(PredAromatic):     return(pres->IsAmino() && AminoProp[pres->refno] & BitAromatic);
	case(PredBasic):        return(pres->IsAmino() && AminoProp[pres->refno] & BitBasic);
	case(PredBuried):       return(pres->IsAmino() && AminoProp[pres->refno] & BitBuried);
	case(PredCharged):      return(pres->IsAmino() && AminoProp[pres->refno] & BitCharged);
	case(PredCyclic):       return(pres->IsAmino() && AminoProp[pres->refno] & BitCyclic);
	case(PredHydrophobic):  return(pres->IsAmino() && AminoProp[pres->refno] & BitHydrophobic);
	case(PredLarge):        return(pres->IsAmino() && !(AminoProp[pres->refno] & BitNotLarge));
	case(PredMedium):       return(pres->IsAmino() && AminoProp[pres->refno] & BitMedium);
	case(PredNeutral):      return(pres->IsAmino() && AminoProp[pres->refno] & BitNeutral);
	case(PredPolar):        return(pres->IsAmino() && !(AminoProp[pres->refno] & BitHydrophobic));
	case(PredSmall):        return(pres->IsAmino() && AminoProp[pres->refno] & BitSmall);
	case(PredSurface):      return(pres->IsAmino() && !(AminoProp[pres->refno] & BitBuried));

	}
	return(True);
}

//static AtomExprVal convertBranchToAtomExprVal(const Branch& b)
//{
//	return std::visit([](auto&& val) -> AtomExprVal {
//		using T = std::decay_t<decltype(val)>;
//			if constexpr (std::is_same_v<T, int> ||
//				std::is_same_v<T, long> ||
//				std::is_same_v<T, double>) {
//				return val;
//			}
//			else {
//				return false; 
//			}
//		}, b);
//}


bool AtomExpr::EvaluateExprFor( HaAtom* aptr )
// Logical, recursive evaluation of the logical expression
{
	AtomExprVal val_lft, val_rgt;

// Within Expression Type

    if( type==OpWithin )
    {
		if (!std::holds_alternative<double>(lft) || !std::holds_alternative<std::shared_ptr<AtomSet>>(rgt))
		{
			PrintLog("Wrong parameters of Within(), should be: ( radius , atom expr ) \n");
			return false;
		}

		double radius = std::get<double>(lft);
		std::shared_ptr<AtomSet> pat_set = std::get<std::shared_ptr<AtomSet>>(rgt);

		if( radius > 0.01 )
        {
			AtomGroup at_grp;
			for (HaAtom* aptr : *pat_set.get())
				at_grp.push_back(aptr);

			return( at_grp.IsWithinRadius(aptr, radius) );
        }
		else
			return( pat_set->count(aptr) > 0 );
    }
	else if (type == OpMember)
	{
		if (!std::holds_alternative<std::shared_ptr<AtomSet>>(rgt))
		{
			PrintLog("Wrong parameters of Member(), should be: ( Atom Set ) \n");
			return false;
		}
		std::shared_ptr<AtomSet> pat_set = std::get<std::shared_ptr<AtomSet>>(rgt);
		return(pat_set->count(aptr) > 0);
	}

    if( type & OpLftVal )
    {
		val_lft = lft;
    }
	else if( type & OpLftProp )
    {
		if (std::holds_alternative<long>(lft))
		{
			long prop = std::get<long>(lft);
			val_lft = EvaluatePropertyFor(aptr, prop);
		}
		else
		{
			PrintLog(" OpLftProp operator expects type long Property Code as an argument \n");
			return false;
		}
    }
	else
	{
		if (std::holds_alternative<std::shared_ptr<AtomExpr>>(lft))
		{
			std::shared_ptr<AtomExpr> pexpr = std::get<std::shared_ptr<AtomExpr >>(lft);
			val_lft = pexpr->EvaluateExprFor(aptr); // resursively evaluate left hand side of the expression
		}
	}

	if (std::holds_alternative<bool>(val_lft) )
	{
		bool bres = std::get<bool>(val_lft);
		if ((OpCode() == OpAnd) && !bres) return(false);
		if ((OpCode() == OpOr) && bres) return(true);
		if (OpCode() == OpNot) return(!bres);
	}
	else
	{
		// if (OpCode() == OpConst) return(val_lft);  //IGOR_TMP - not sure how OpConst processed
		//else
		//{
		//	PrintLog("At this point result of the evaluation of the left branch of the expression should be Const Value - return false \n");
		//	return false;
		//}
	}

// Evaluation of the right side of the expression

    if( type & OpRgtVal )
    {
		val_rgt = rgt;
    }
	else if( type & OpRgtProp )
    {
		if (std::holds_alternative<long>(rgt))
		{
			long prop = std::get<long>(rgt);
			val_rgt = EvaluatePropertyFor(aptr, prop);
		}
		else
		{
			PrintLog(" OpRgtProp operator expects type long Property Code as an argument \n");
			return false;
		}
    }
	else
	{
		if (std::holds_alternative<std::shared_ptr<AtomExpr>>(rgt))
		{
			std::shared_ptr<AtomExpr> pexpr = std::get<std::shared_ptr<AtomExpr >>(rgt);
			val_rgt = pexpr->EvaluateExprFor(aptr); // resursively evaluate right hand side of the expression
		}
	}

    switch( OpCode() )
    {
		case(OpOr):
        case(OpAnd):     
		{
			if (std::holds_alternative<bool>(val_rgt))
			{
				bool bres = std::get<bool>(val_rgt);
				return(bres);
			}
			else
			{
				PrintLog("At This point val_rgt should be bool \n");
				return false;
			}
		}

		case(OpLess):    return( val_lft < val_rgt );
        case(OpMore):    return( val_lft > val_rgt );
        case(OpEqual):   return( val_lft == val_rgt );
        case(OpNotEq):   return( val_lft != val_rgt );
        case(OpLessEq):  return( val_lft <= val_rgt );
        case(OpMoreEq):  return( val_lft >= val_rgt );
    }
    return( true );
}

std::shared_ptr<AtomExpr> AtomExpr::LookUpAtGroupExpr( std::string grp_name, MolSet* pmset )
{
	std::shared_ptr<AtomExpr> p_expr = AtomExpr::CreateFalseExpr();

    AtomGroup* atl;
	atl = pmset->GetAtomGroupByID(grp_name);
	if(atl == NULL) return p_expr;
	
	std::shared_ptr<AtomSet> pat_set;
	for( HaAtom* aptr : *atl)
		pat_set->insert(aptr);

    p_expr->type = OpMember;
	p_expr->rgt = pat_set;

    return( p_expr );
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


std::shared_ptr<AtomExpr> AtomExpr::LookUpElement(const char* ident)
{
	std::shared_ptr<AtomExpr> pexpr = AtomExpr::CreateFalseExpr();
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

	if (elem < MAXELEMNO)
	{
		pexpr->type = OpEqual | OpLftProp | OpRgtVal;
		pexpr->lft = (long)PropElemNo;
		pexpr->rgt = (long)elem;
	}

    return( pexpr );
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


std::shared_ptr<AtomExpr> AtomExpr::ParsePrimitiveExpr(MolSet* pmset, std::string expr_str, size_t& cr_pos )
//! 
//! Generate a logical expression of subset of atoms of pmset corresponding to a string token 
//! \param pmset    - Molecular Set atom expression is formed on
//! \param expr_str& - string with expression
//! \param cr_pos&   - current cursor position (input and output)
{
	std::string NameBuf;
    std::shared_ptr<AtomExpr> tmp1,tmp2;
	std::shared_ptr<AtomExpr> wild;
    
    int i, j;
    int neg;
    int ch;

	std::shared_ptr<AtomExpr> p_expr = CreateFalseExpr();
	std::shared_ptr<AtomExpr> p_expr_false = CreateFalseExpr();

	if( cr_pos >= expr_str.size() ) return p_expr_false;
 
	try
	{
		cr_pos++;
		if (cr_pos >= expr_str.size() )
		{
				cr_pos--;
				return p_expr;
		}
		ch = expr_str[cr_pos];

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
					tmp1 = std::make_shared<AtomExpr>();
					tmp1->type = OpEqual | OpLftProp | OpRgtVal;
					tmp1->lft = (long)PropMolPtr;
					tmp1->rgt = pMol;
				
					tmp2 = std::make_shared<AtomExpr>();
					tmp2->type = OpOr;
					tmp2->lft = tmp1;
					tmp2->rgt = wild;
					wild = tmp2;	
				}
			}
			p_expr = wild;
			cr_pos++;
			if (cr_pos >= expr_str.size())
			{
				cr_pos--;
				return p_expr;
			}
			ch = expr_str[cr_pos];
		} // end of Molecule name parsing

		if( isspace(ch))
		{
			cr_pos--;
			return p_expr;
		}

		if( ch != ':' ) // Parse Residue Name 
		{   
			NameBuf.clear();
			if( ch != '*' )
			{
				if( ch == '[' ) //Parse residue name in [ ] brackets
				{
					while( true )
					{
						cr_pos++;
						if ( cr_pos >= expr_str.size() )  throw std::runtime_error("No closing ] in Residue Name");
						ch = expr_str[cr_pos];
						if (ch == ']') break;
						NameBuf += toupper(ch);
					}
				}
				else // Parse residue name without brackets
				{
					for( i=0; i < 20; i++ )
					{
						if( isalpha(ch) )
						{
							NameBuf += toupper(ch);
							cr_pos++;
							if (cr_pos >= expr_str.size()) break;
							ch = expr_str[cr_pos];
						}
						else if( (ch == '?') || (ch == '%') )
						{
							NameBuf += '?';
							cr_pos++;
							if (cr_pos >= expr_str.size()) break;
							ch = expr_str[cr_pos];
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
						tmp1 = std::make_shared<AtomExpr>();
						tmp1->type = OpEqual | OpLftProp | OpRgtVal;
						tmp1->lft = (long)PropResName;
						tmp1->rgt = (long)j;
					
						tmp2 = std::make_shared<AtomExpr>();
						tmp2->type = OpOr;
						tmp2->lft = tmp1;
						tmp2->rgt = wild;
						wild = tmp2;
					}
				}
			
				if( !IsTrueExpr(p_expr.get()) )
				{
					tmp2 = std::make_shared<AtomExpr>();
					tmp2->type = OpAnd;
					tmp2->lft = wild;
					tmp2->rgt = p_expr;
					wild = tmp2;
				}
				p_expr = wild;
			}
			else // ch != '*' - for residue name
			{
				cr_pos++;
				if (cr_pos > expr_str.size()) return p_expr;
				ch = expr_str[cr_pos];
			}

			if( ch != '*') // Parse Residue Number 
			{
				if( ch == '-' )
				{
					cr_pos++;
					if (cr_pos >= expr_str.size()) return p_expr;
					ch = expr_str[cr_pos];
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
						if (cr_pos >= expr_str.size()) break;
					}

					tmp1 = std::make_shared<AtomExpr>();
					tmp1->type = OpEqual | OpLftProp | OpRgtVal;
					tmp1->rgt = neg? (long)-i : (long)i;
					tmp1->lft = (long) PropResId;
					if( !IsTrueExpr(p_expr.get()) )
					{
						tmp2 = std::make_shared<AtomExpr>();
						tmp2->type = OpAnd;
						tmp2->rgt = p_expr;
						tmp2->lft = tmp1;
						p_expr = tmp2;
					}
					else
					{
						p_expr = tmp1;
					}
					if (cr_pos >= expr_str.size()) return p_expr;
					cr_pos++;
					if (cr_pos >= expr_str.size()) return p_expr;
					ch = expr_str[cr_pos];
				}
				else if( neg )
				{
					throw std::runtime_error("Only minus sign in place of residue number");
				}
			}
			else  // if (ch == '*') in place of residue number
			{
				cr_pos++;
				if (cr_pos >= expr_str.size()) return p_expr;
				ch = expr_str[cr_pos];
			}
		}	 

		if( ch==':' ) // Parse Chain Ident 
		{
			cr_pos++;
			if (cr_pos >= expr_str.size()) return p_expr;
			ch = expr_str[cr_pos];
		}

		if( isalnum(ch) )
		{
			ch = toupper(ch);

			tmp1 = std::make_shared<AtomExpr>();
			tmp1->type = OpEqual | OpLftProp | OpRgtVal;
			tmp1->lft = (long) PropChain;
			tmp1->rgt = (char) ch;
			if( !IsTrueExpr(p_expr.get()) )
			{
				tmp2 = std::make_shared<AtomExpr>();
				tmp2->type = OpAnd;
				tmp2->rgt = p_expr;
				tmp2->lft = tmp1;
				p_expr = tmp2;
			}
			else
			{
				p_expr = tmp1;
			}
			cr_pos++;
			if (cr_pos >= expr_str.size()) return p_expr;
			ch = expr_str[cr_pos];
		}
		else if( (ch=='?') || (ch=='%') || (ch=='*') )
		{
			cr_pos++;
			if (cr_pos >= expr_str.size()) return p_expr;
			ch = expr_str[cr_pos];
		}

		if( ch == ':' ) // Parse Model Number 
		{
			cr_pos++;
			if (cr_pos >= expr_str.size()) return p_expr;
			ch = expr_str[cr_pos++];
			if( isdigit(ch) )
			{
				i = ch-'0';
				while( isdigit(expr_str[cr_pos]) )
				{
					i = 10*i + (expr_str[cr_pos]-'0');
					cr_pos++;
					if (cr_pos >= expr_str.size()) break;
				}

				tmp1 = std::make_shared<AtomExpr>();
				tmp1->type = OpEqual | OpLftProp | OpRgtVal;
				tmp1->lft = (long)PropModel;
				tmp1->rgt = (long)i;
				if( IsTrueExpr(p_expr.get()) )
				{
					tmp2 = std::make_shared<AtomExpr>();
					tmp2->type = OpAnd;
					tmp2->rgt = p_expr;
					tmp2->lft = tmp1;
					p_expr = tmp2;
				}
				else
				{
					p_expr = tmp1;
				}
				if (cr_pos >= expr_str.size()) return p_expr;
				cr_pos++;
				if (cr_pos >= expr_str.size()) return p_expr;
				ch = expr_str[cr_pos];
			}
			else
			{
				throw std::runtime_error("No model number after : ");
			}
		} // End Parse Model Number

		if( ch == '.' ) // Parse Atom Name 
		{
			NameBuf.clear();
			cr_pos++;
			if (cr_pos >= expr_str.size()) return p_expr;
			ch = expr_str[cr_pos];
			if( ch != '*' )
			{
				for( i=0; i< 20; i++ )
				{
					if( isalnum(ch) || ch=='\'' || ch=='*' )
					{
						NameBuf += toupper(ch);
						cr_pos++;
						if(cr_pos >= expr_str.size()) break;
						ch = expr_str[cr_pos];
					}
					else if( (ch=='?') || (ch=='%') || (ch=='#') )
					{
						NameBuf += '?';
						cr_pos++;
						if (cr_pos >= expr_str.size()) break;
						ch = expr_str[cr_pos];
					}
					else
					{
						break;
					}
				}
				if( NameBuf.empty() ) throw std::runtime_error("Empty Atom Name");  
			
				wild = CreateFalseExpr();
				for( j=0; j< HaAtom::ElemDesc.size(); j++ )
				{
					if( MatchWildName(NameBuf.c_str(), HaAtom::ElemDesc[j].c_str(), HaAtom::ElemDesc[j].size(),NameBuf.size()) )
					{
						tmp1 = std::make_shared<AtomExpr>();
						tmp1->type = OpEqual | OpLftProp | OpRgtVal;
						tmp1->lft = PropName;
						tmp1->rgt = (long) j;

						tmp2 = std::make_shared<AtomExpr>();
						tmp2->type = OpOr;
						tmp2->lft = tmp1;
						tmp2->rgt = wild;

						wild = tmp2;
					}
				}

				if( IsTrueExpr(p_expr.get()) || IsFalseExpr(wild.get()) )
				{
					p_expr = wild;
				}
				else
				{
					tmp1 = std::make_shared<AtomExpr>();
					tmp1->type = OpAnd;
					tmp1->lft = p_expr;
					tmp1->rgt = wild;
					p_expr = tmp1;
				}
			}
			else // if( ch == '*' )
			{
				cr_pos++;
				if( cr_pos >= expr_str.size())
				ch = expr_str[cr_pos];
			}
		} // if(ch == '.') End Parse Atom Name
	}
	catch( std::exception& ex)
	{
		PrintLog(" Parsing Error %s \n", ex.what());
		cr_pos--;
		return p_expr;
	}
    cr_pos--;
    if( !ch || isspace(ch) || ispunct(ch) )
	{
		return p_expr;
	}
	return p_expr;
} 

std::shared_ptr<AtomExpr> AtomExpr::ParseExpression(const std::string& expr_str, MolSet* pmset)
{
	CmdParser cmd_parser;
	cmd_parser.SetCmdLine(expr_str);
	cmd_parser.FetchToken();
	std::shared_ptr<AtomExpr> p_expr = cmd_parser.ParseExpression(0,pmset);
	if( cmd_parser.CurToken )
	{   
		PrintLog("Error in AtomExpr::ParseExpression() \n");
		PrintLog("Invalid Expression String syntax \n");
		p_expr = std::make_shared<AtomExpr>();
	} 
	return p_expr; 
}
