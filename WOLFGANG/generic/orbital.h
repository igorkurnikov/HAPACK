/*TEX
%
% IPACK - GENERIC
% Quantum Chemistry Project: GENERIC TYPES
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is proprietary software. Unlicensed use and distribution
% of this program or parts thereof are illegal. For details on the 
% license, see the LICENSE section in the file "main.c" distributed 
% with this package or write to: 
%      Wolfgang Wenzel, Theoretical Physics I, 
%      University of Dortmund,
%      D-44221 Dortmund, Germany 
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the LICENSE
% for more details.
%
% $Id: orbital.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $

\section{Types of Orbitals}
The most important criterion to differentiate among the various
orbitals in this calculation is that of orbital type, so a specialized
class is provided here. Orbitals can be either:
\begin{tabular}{ll} 
\name{freeze}   & these orbitals will be frozen \\
\name{soft\_core}  & limited excitations from these orbitals \\
\name{internal_type} & in the sense of an MRCI calculation \\
\name{external} & in the sense of an MRCI calculation \\
\name{bsr\_type}      & these orbitals are treated perturbatively \\
\name{virt}     & these orbitals will be thrown away
\end{tabular} 
These types are defined in an enumeration.

A conversion operator is supplied to ease the use of
\name{OrbitalTypes} as indices of arrays. 
*/
 
#ifndef ORBITAL_DEF
#define ORBITAL_DEF

#include "io_incl.h"
#include <assert.h>
#include <numer.h>
#include "symmetry.h"

#if !defined(GNU) && !defined(_MSC_VER) 
#ifndef BOOLEAN 
#define BOOLEAN
enum Boolean  { false, true  };
#endif
#endif

enum otype {internal_type, external, bsr_type, maxtype}; 

class OrbitalType
{
  int t;
public:
int in_range() const { return (t >= 0 && t < maxtype); } 
                     OrbitalType() { t = internal_type; }   
                     OrbitalType( const otype tt )  { t = tt; }   
OrbitalType&         operator++()         { t++; return *this; } 
int                  operator()() const   { assert(in_range()); return t; }    		      
                     operator int() const { assert(in_range()); return t; } 

int                  operator==(const OrbitalType& tt) const { return t == tt.t;} 
int                  operator!=(const OrbitalType& tt) const { return t != tt.t;} 
int                  operator< (const OrbitalType& tt) const { return t <  tt.t;} 
int                  operator<=(const OrbitalType& tt) const { return t <= tt.t;} 
int                  operator> (const OrbitalType& tt) const { return t >  tt.t;} 

};

ostream & operator<<(ostream & os,const OrbitalType& );
istream & operator>>(istream & is,OrbitalType& type);

/*TEX 
\section{Orbital Class.}  

Each orbital is characterized by a triplet of \name{OrbitalType}, 
\name{Symmetry} and
an integer value that distinguishes orbitals of identical type and 
symmetry. Orbitals are ordered, first by their type, then by their
symmetry and finally by their index. 

\begin{tabular}{p{3cm}p{10cm}}
  s()     & symmetry of orbital   \\
  i()     & index of orbital      \\
  type()  & type of orbital (internal or external) \\
  =       & assignment of two orbitals \\
\name{in\_range} & whether checks whether the orbital is legal. 
\end{tabular} 

The \name{number\_of\_orbitals} function returns the number of
orbitals of a given \name{OrbitalType} or of a given
\name{OrbitalType} and \name{Symmetry}. The \name{number} function
assigns a linear index to all orbital. The numbering agrees with the
ordering of the orbitals.

The static private  array \name{no}  is indexed by \name{OrbitalType}
and \name{Symmetry} (in that order) and records the number of orbitals
of a given type and  symmetry. Similarly \name{start\_no} returns the
index of the first orbital of a given  type and  symmetry.
The increase operator runs through all orbitals in the order defined
above.
*/

#include <numer.h>
#include "vtype.h"
class DataVault;

enum OInitCode  { oinit_generic,oinit_mrci };
 
class Orbital
{
private:
static NUMARRAY2D<int> no;	    
static NUMARRAY2D<int> start_no;     
static NUMARRAY<int>   skip; 
static ARRAY<Orbital>  orb;
static ARRAY<int>      orb_no;
static ARRAY<int>      inactive_id;
static ARRAY<int>      extidx;
static int             no_ext_orb;
static int             inactive_no;
static NUMARRAY<int>   map;

friend class OrbitalIter;
friend class OrbitalSymIter;

int         in_range()              { return v_type.in_range(); }     
Orbital&    operator++(); 

protected:
  int         v_i;
  Symmetry    v_s;
  OrbitalType v_type;

public:
            Orbital()
              : v_i(0), v_type(internal_type)             { v_i = 0; } 
            Orbital(OrbitalType intype,Symmetry sym,int index)
              : v_s(sym), v_i(index)  { v_type = intype; }
            Orbital(OrbitalType intype)
              : v_s(), v_i(0) { v_type = intype; }
            Orbital(Symmetry sym) 
              : v_s(sym), v_i(0), v_type(internal_type)   { }

           ~Orbital()                                { }   

Symmetry    symmetry()    const  { return v_s;    }
int         index()       const  { return v_i;    }
OrbitalType type()        const  { return v_type; }

static int  no_inactive() { return inactive_no; }
static int  number_of_orbitals();
static int  number_of_orbitals(const OrbitalType& type);            
static int  number_of_orbitals( const Symmetry& s );
static int  number_of_orbitals(const OrbitalType& type, const Symmetry& s )
            { return no(type,s); } 
static int  number_of_disk_orbitals() { return no_ext_orb; }
static Orbital& orbital(int number)   { return orb[number]; }
static int      orbno_from_external(int exti) { return orb_no[exti]; }
static Orbital& orb_from_external  (int exti) { return orbital(orb_no[exti]); }
static int      inactive_idx(int ono) { return inactive_id(ono);  }
int         number() const          { return start_no(v_type,v_s) + v_i; }
int         external_index() const  { return extidx[number()];  }

  
static void init(istream& in,DataVault& data,
		 const OInitCode icode = oinit_generic);
static void clear();
static void extidx_set_up();
static void print_table( ostream & os, int verbose = 1 );
static int  set_up_map( const int&, const int& );
static NUMARRAY<int>&  pass_map() { return map; }
                           
};

int         operator< (const Orbital& o1,const Orbital& o2);
int         operator> (const Orbital& o1,const Orbital& o2);
int         operator<=(const Orbital& o1,const Orbital& o2);
int         operator>=(const Orbital& o1,const Orbital& o2);
int         operator==(const Orbital& o1,const Orbital& o2);
int         operator!=(const Orbital& o1,const Orbital& o2);

istream & operator>>(istream & is,Orbital& );
ostream & operator<<(ostream & os,const Orbital& );

/*TEX
\section{OrbitalIter}
Further we provide a class which acts as a iterator over all orbitals
of a given type.
*/

enum OrbitalTypeRange { ALL_TYPE,INT_TYPE,EXT_TYPE,BSR_TYPE,IX_TYPE };  
  
class OrbitalIter
{
  OrbitalTypeRange type;
  Orbital          orb;
  int              ok;
public:  
               OrbitalIter(OrbitalTypeRange target_type);
OrbitalIter&   operator++();
int            in_range()           { return ok;  }  
const Orbital& operator()() const   { return orb; } 
};

class OrbitalSymIter
{
  Symmetry         tsym;
  OrbitalTypeRange type;
  Orbital          orb;
  int              ok;
public:  
                OrbitalSymIter(OrbitalTypeRange target_type,Symmetry s);
OrbitalSymIter& operator++();
int             in_range()           { return ok;  }  
const Orbital&  operator()() const   { return orb; } 
};

#endif

    

