/*! \file abstree.h

   Logical expression class in HARLEM

   \author Igor Kurnikov
   \date 1997-2010   

   derived from:
   abstree.h
   RasMol2 Molecular Graphics
   Roger Sayle, August 1995
   Version 2.6
*/

#if !defined(ABSTREE_H)
#define ABSTREE_H


/* Operator Types */
const int OpAnd =        0x01;
const int OpOr  =        0x02;
const int OpNot =        0x03;
const int OpEqual =      0x04;
const int OpNotEq =      0x05;
const int OpLess =       0x06;
const int OpMore =       0x07;
const int OpLessEq =     0x08;
const int OpMoreEq =     0x09;
const int OpConst =      0x0a;
const int OpWithin =     0x0b;
const int OpMember =     0xac;

const int OpLftProp =    0x10;
const int OpLftVal  =    0x20;
const int OpRgtProp =    0x40;
const int OpRgtVal  =    0x80;

/* Property fields */
const int PropIdent =       1;
const int PropXCord =       2;
const int PropYCord =       3;
const int PropZCord =       4;
const int PropTemp  =       5;
const int PropRad   =       6;
const int PropResId =       7;
const int PropName  =       8;
const int PropChain =       9;
const int PropResName =     10;
const int PropSelect  =     11;
const int PropElemNo  =     12;
const int PropModel   =     13;
const int PropChemGroup =   14;
const int PropMolPtr  =     15;
const int PropAtGroup =     16;

inline int PredAbsOrd( int x )  { return (x - 20); }
inline int PredAbsChr( int x )  { return (x + 20); }

const int PredAlpha =     20;
const int PredAmino =     21;
const int PredAT    =     22;
const int PredBonded =    23;
const int PredCG     =    24;
const int PredCystine =   25;
const int PredDNA     =   26;
const int PredHelix   =   27;
const int PredHetero  =   28;
const int PredHydrogen =  29;
const int PredIon     =   30;
const int PredLigand  =   31;
const int PredMainChain = 32;
const int PredNucleic =   33;
const int PredProtein =   34;
const int PredPurine  =   35;
const int PredPyrimidine = 36;
const int PredRNA     =    37;
const int PredSelected =   38; /* Unused! */
const int PredSheet    =   39;
const int PredSidechain =  40;
const int PredSolvent   =  41;
const int PredTurn      =  42;
const int PredWater	   =  43;

const int PredAcidic    =  44;
const int PredAcyclic   =  45;
const int PredAliphatic =  46;
const int PredAromatic  =  47;
const int PredBasic     =  48;
const int PredBuried    =  49;
const int PredCharged   =  50;
const int PredCyclic    =  51;
const int PredHydrophobic = 52;
const int PredLarge     =  53;
const int PredMedium    =  54;
const int PredNeutral   =  55;
const int PredPolar     =  56;
const int PredSmall     =  57;
const int PredSurface   =  58;


#include "haatom.h"

class HaMolSet;

class AtomExpr;
class AtomContainer;

typedef union {
	AtomContainer  *set;
	AtomExpr *ptr;
    int limit;
	long val;
	double dval;
	} Branch;  //!< A branch of a logical expression


class HaResidue;
class HaChain;


class AtomExpr
//! Class to define logical or arithmetical expression on atoms
{
public:
   AtomExpr();
   AtomExpr(int new_type, long rval, long lval);
   virtual ~AtomExpr();

   int type;     //!< type of logical operation OpAnd, OpOr .. OpWithin etc
   Branch rgt;   //!< left part of the expression
   Branch lft;   //!< right part of the expression

   long EvaluateExprFor(HaAtom* aptr);  //!< evaluate the expression for the current atom
   long EvaluateProperty(long prop );   //!< Evaluate given Property for a Current Atom

   static AtomExpr* LookUpAtGroupExpr( const char* grp_name, HaMolSet* pmset); //!< Create expression of atom group name
   static AtomExpr* LookUpElement( const char* );    //!< create expression of element name       
   static AtomExpr* ParsePrimitiveExpr( HaMolSet* pmset, const char* expr_str, size_t& cr_pos); //!< Build Logical expression corresponding to a text of atom specifications (possibly wild-carded)

   static AtomExpr* ParseExpression(const std::string& expr_str, HaMolSet* pmset); //!< Build Atom expression from expression string for Molecular set pmset

   static AtomExpr* CreateTrueExpr(); //!< Create TRUE expression
   static AtomExpr* CreateFalseExpr();  //!< Create FALSE Expression
   static bool IsTrueExpr(AtomExpr* p_expr); //!< Check if Expression Is TRUE
   static bool IsFalseExpr(AtomExpr* p_expr); //!< Check if Expression Is FALSE


protected:

   inline int OpCode() { return(type & 0x0f); } //!< return  code of operation
   long EvaluateExpr();                         //!< evaluate expression for current QAtom, QGroup, QChain
	
   static HaAtom*      QAtom;    //!<  Current Atom
   static HaResidue*   QGroup;   //!<  Current Residue
   static HaChain*     QChain;   //!<  Current Chain

};

const double SelectRad = 0.4;

#ifdef ABSTREE_CPP

HaAtom  *PkAtom=NULL;

#else

extern HaAtom  *PkAtom;

#endif // end ifdef ABSTREE_CPP

#endif // end if !defined(ABSTREE_H)
