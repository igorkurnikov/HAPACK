/*TEX
%
% IPACK - QCHEM
% Quantum Chemistry Project: Integral Package
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is proprietary software. Unlicensed use and distribution
% of this program or parts thereof are illegal. For details on the 
% license, see the LICENSE section in the file "ipack.c" distributed 
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
% $Id: typesip.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\chapter{Auxiliary Classes\label{auxiliary-classes}}
This section implements two simple classes dealing with 
\name{Location}s and angular momenta, which are useful in the
formulation of the basis sets and recursion relations.
The class \nameindex{FlagList} provides a set of flag variables.
In addition we define a type \nameindex{Direction} which can take the 
values X,Y and Z for the three cartesian directions. Several structures
for passing data are definde at the end of this chapter.
*/
#ifndef TYPES_DEF
#define TYPES_DEF

#include <math.h>
#include <xtime.h>
#include "const.h"
#ifdef MEMCHK_IPACK
#include "memchk.h"
#endif

#if !defined(GNU) && !defined(_MSC_VER) && !defined(__DECCXX)
#ifndef BOOLEAN
#define BOOLEAN
enum Boolean   { false,true};
#endif
#endif 

enum Direction { X,Y,Z };
enum TRANS     { TRANS12 = 1, TRANS34 = 2, TRANSPAIR = 4 };

#include "symop.h"

inline Direction& op_plus(Direction& dir)
{
  dir = (Direction) (dir + 1);
  return dir;
}

inline Direction& operator++(Direction& dir)
{
  dir = (Direction) (dir + 1);
  return dir;
}


// dead: void print(ostream& os, Direction axs);

istream& operator>>(istream& is,Direction& d);
ostream& operator<<(ostream& os,Direction d);

/*TEX
\section{Class FlagList\label{section-flags}}
\index{FlagName}  
\index{FlagList}
\index{flaglist, extern}  
We define a class FlagList which summerizes several switches. They are set
by the user at input time, depending on the input file. This is done in
the BasisSet::parse function. The user can choose which kind of integrals 
he wants to be calculated and he can select whether spherical or 
cartesian functions are generated. The angstrom flag allows distances to be
entered in angstrom although atomic units are always used internally.
The constructor of FlagList is used to predefine default values.
*/
#include <numer.h>
enum FlagName { Overlap,Kinetic,Nuclear,TwoElec,Moment,eField,Angular,SpinOrb,
	     Spherical,Angstrom};

class FlagList
{
 public: 

   IVector list;

   FlagList()   // default settings
   {
     list.reset(10);
     
     list[Overlap]   = true;
     list[Kinetic]   = true;
     list[Nuclear]   = true;
     list[TwoElec]   = true;
     list[Moment]    = false;
     list[eField]    = false;
     list[Angular]   = false;
     list[SpinOrb]   = false;
     list[Spherical] = false;
     list[Angstrom]  = false;
   }
   ~FlagList(){}

   const int& operator()(const FlagName i) const{ return list(i); }
   int& operator[](const FlagName i) { return list[i]; }
};

extern FlagList flaglist; // global flaglist defined in input.c 

/*TEX
\section{Class Location}

This class implements center coordinates, it is so simple that it is
self-explanatory.  
*/
class Location
{
protected:  
  double pos[3];
public:
  Location() { pos[X] = 0; pos[Y] = 0; pos[Z] = 0; }
  Location(const double xxxx,const double yyyy,const double zzzz) 
  { pos[X] = xxxx; pos[Y] = yyyy; pos[Z] = zzzz; }
  Location(const Location& l) 
  { pos[X] = l.pos[X]; pos[Y] = l.pos[Y]; pos[Z] = l.pos[Z]; }
  ~Location() {}
Location& assign(const Location& l) 
  { 
    pos[X] = l.pos[X]; pos[Y] = l.pos[Y]; pos[Z] = l.pos[Z];
    return *this; 
  }
Location& operator=(const Location& l) { return assign(l); } 
double x() const { return pos[X]; }
double y() const { return pos[Y]; }
double z() const { return pos[Z]; }
double operator()(const Direction& d) const;
double coord(const Direction& axs)  { return pos[axs]; }
int    sym_equiv(Direction d) const { return fabs(pos[d]) < delta; }
void   apply(const SymOp& s)  
       { 
	 for(int i=0; i< 3; i++)
	   if (s(i) == reflection) pos[i] = - pos[i];
       }
};

inline int operator==(const Location& a , const Location& b)
{
  return (fabs(a.x()-b.x()) < delta && fabs(a.y()-b.y()) < delta && 
	  fabs(a.z()-b.z()) < delta);
}

inline int operator!=(const Location& a , const Location& b)
{
  return (fabs(a.x()-b.x()) >= delta || fabs(a.y()-b.y()) >= delta || 
	  fabs(a.z()-b.z()) >= delta);
}

ostream& operator<<(ostream& os,const Location& l);
istream& operator>>(istream& os,Location& l);

/*TEX
\section{Angular Momenta}
\index{Momentum}

This section implements the class \name{Momentum} to deal with
the cartesian angular momentum triplets $(n_x,n_y,n_z)$ which index
the angular dependence of all orbitals in this calculation. The total
angular momentum of such a triplet is defined as: $L= n_x + n_y + n_z$.

Upon initialization a \name{Momentum} is initialized to point entirely
in the $x$ direction. The \name{operator++} implements a hirarchy and
can be used in loops. The \name{operator()} can be used to check
\index{Momentum,operators}
whether the current angular momentum is still in range.
This allows loops of the type:

\index{Momentum,loops over}
\begin{verbatim}
     for(Momentum mu(ltot);mu(); ++mu())
     { 
        ....
     }
\end{verbatim}
 
The function \name{decrease} returns a direction $d$ for which $n_d >
0$ if such a direction exists. The function \name{reduce(d)} reduces
the total $L$ and the component in direction $d$ by one.  The function
\name{induce(d)} increases the total $L$ and the component in
direction $d$ by one.

\index{Momentum,decrease}
\index{Momentum,reduce}
\index{Momentum,index}
\index{Momentum,local\_index}

The remaining functions deal with the assignment of linear indices to
a momentum variable. We are cheifly interested in two types of linear
indices: For a given total momentum L, we wish to map the
\name{Momentum} onto a linear index ranging from zero to
\nameindex{number(L)}-1, which designates the total number of such
momenta. This is done by the function \name{local\_index}.
Sometimes, we would also like to know the linear index of a
\name{Momentum} in the set of all mometa of total momemtum {\em less
or equal} to L. This is done by the function \nameindex{index()}. The 
linear index of the first momentum of total momentum L is given by
\nameindex{min}(L).

\subsection{Class Header}

*/

class Momentum
{
       int L;
       int xx,yy,zz;

int         start(const int x) const { return  (L-x)*(L-x+1)/2; } 
int         start() const { return start(xx); }

public:
       Momentum(const int l=0) { L = xx = l; yy = 0;  zz= 0; }
       Momentum(const Momentum& l) { xx = l.xx; yy = l.yy;zz = l.zz; L = l.L;}
       Momentum(const int ax,const int ay,const int az) 
               { xx = ax; yy = ay;zz = az; L = ax+ay+az;}
       Momentum(const char* s);
      ~Momentum() {}
int         x() const { return xx; }
int         y() const { return yy; }
int         z() const { return zz; }
int         l() const { return xx+yy+zz; }
Momentum&   operator=(const Momentum& l)
            { xx = l.xx; yy = l.yy; zz = l.zz; L=l.L; return *this; }
int         operator++();
int         operator()() { return (xx >= 0); }   
int         operator()(const Direction& dir) const;
Direction   decrease() const;
void        reduce(const Direction& dir);
void        induce(const Direction& dir);
  

static int  number(const int ll) { return (ll+1)*(ll+2)/2; }
static int  min(const int l)     { return  l*(l+1)*(l+2)/6; }
int         min() const          { return  min(L); }
int         local_index() const  { return start(xx) + zz; } 
int         index() const  { return min(L) + local_index(); } 
int         operator==(const Momentum& b)
            { return (xx == b.xx && yy == b.yy && zz == b.zz); } 
int         operator<(const Momentum& b)
            {
	      if (xx < b.xx) return 1;
	      if (xx == b.xx && yy  < b.yy) return 1;
	      if (xx == b.xx && yy == b.yy && zz < b.zz) return 1;
	      return 0;
	    } 
int         operator!=(const Momentum& b) { return !(*this == b); } 
             
}; 

ostream& operator<<(ostream& os,const Momentum l);
istream& operator>>(istream& os,Momentum& l);

/*TEX  
\section{Epsilon Operator}
The \nameindex{Epsilon Operator} is needed for evaluating cross
products.
*/
inline double EpsOp(const Direction& i,const Direction& j,const Direction& k) 
{
  switch(i)
  {
  case X:
    switch(j)
    {
    case X:
      return 0;
    case Y:
      if(k == Z)
        return 1;
      else
        return 0;
    case Z:
      if(k == Y)
        return -1;
      else
        return 0;
    }
  case Y:
    switch(j)
    {
    case X:
      if(k == Z)
	return -1;
      else
	return 0;
    case Y:
        return 0;
    case Z:
      if(k == X)
        return 1;
      else
        return 0;
    }
  case Z:
    switch(j)
    {
    case X:
      if(k == Y)
	return 1;
      else
	return 0;
    case Y:
      if(k == X)
        return -1;
      else
        return 0;
    case Z:
        return 0;
    }
  }
  return 0;
}    
      
   
    

/*TEX  
\section{MomentData}
The structure \nameindex{MomentData} is used for passing data needed for 
moment integrals.  
*/

struct MomentData
{
  Location loc;
  int      max_mu;
};

#endif

