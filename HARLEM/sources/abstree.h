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


#include <variant>
#include "haatom.h"
#include "haatgroup.h"


class MolSet;

class AtomExpr;
class AtomContainer;

//typedef union {
//	AtomContainer  *set;
//	AtomExpr *ptr;
//   int limit;
//	long val;
//	double dval;
//	} Branch;  //!< A branch of a logical expression

//using Branch = std::variant<
//	std::shared_ptr<AtomContainer>,  // replaces AtomContainer*
//	std::shared_ptr<AtomExpr>,       // replaces AtomExpr*
//	int,
//	long,
//	double,
//	void*
//>;

using AtomExprVal = std::variant<
	bool,
	char,
	int,
	long,
	double,
	std::string,
	HaMolecule*,
	HaResidue*,
	HaChain*,
	std::shared_ptr<AtomExpr>,
	std::shared_ptr<AtomSet>
>;

class HaResidue;
class HaChain;

class AtomExpr
//! Class to define logical or arithmetical expression on atoms
{
public:
   AtomExpr();
   AtomExpr(int new_type, AtomExprVal rval, AtomExprVal lval);
   virtual ~AtomExpr();

   int type;     //!< type of logical operation OpAnd, OpOr .. OpWithin etc
   AtomExprVal rgt;   //!< left part of the expression
   AtomExprVal lft;   //!< right part of the expression

   AtomExprVal EvaluatePropertyFor(HaAtom* aptr, long prop);  //!< evaluate property for the atom, 

   static std::shared_ptr<AtomExpr> LookUpAtGroupExpr( std::string grp_name, MolSet* pmset); //!< Create expression of atom group name
   static std::shared_ptr<AtomExpr> LookUpElement( const char* elem_name );    //!< create expression of element name       
   static std::shared_ptr<AtomExpr> ParsePrimitiveExpr( MolSet* pmset, std::string expr_str, size_t& cr_pos); //!< Build Logical expression corresponding to a text of atom specifications (possibly wild-carded)

   static std::shared_ptr<AtomExpr> ParseExpression(const std::string& expr_str, MolSet* pmset); //!< Build Atom expression from expression string for Molecular Set pmset

   static std::shared_ptr<AtomExpr> CreateTrueExpr(); //!< Create TRUE expression
   static std::shared_ptr<AtomExpr> CreateFalseExpr();  //!< Create FALSE Expression
   static bool IsTrueExpr(AtomExpr* p_expr); //!< Check if Expression Is TRUE
   static bool IsFalseExpr(AtomExpr* p_expr); //!< Check if Expression Is FALSE

   // static AtomExprVal convertBranchToAtomExprVal(const Branch& b);

   inline int OpCode() { return(type & 0x0f); } //!< return  code of operation
   AtomExprVal AtomExpr::EvaluateExprFor(HaAtom* aptr);
};

const double SelectRad = 0.4;

#ifdef ABSTREE_CPP

HaAtom  *PkAtom=NULL;

#else

extern HaAtom  *PkAtom;

#endif // end ifdef ABSTREE_CPP

#endif // end if !defined(ABSTREE_H)
