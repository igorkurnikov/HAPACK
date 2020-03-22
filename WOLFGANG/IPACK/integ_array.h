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
% $Id: integ_array.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\chapter{Integral Arrays}

\section{Definition}

An \name{IntegArray} is the basic class to store the radial-function
dependent information for a given pair/quadruplet of angular momentum
functions. The definition of an \nameindex{IntegArray} is given in the
discussion of the classes \name{Recursion} and \name{VRR}.

Technically the \name{IntegArray} is simply a matrix/tensor of fourth
rank, implemented through derivation from the \name{Storage} class. 

There are two constructors.  This is done by the \name{Storage} copy
function. Possibly auxiliary indices are required. For secondary, i.e
recursively defined \name{IntegArrays}, we merely allocate storage and
initialize all terms to zero.

In the context of a recursion relation we must consider the following
operations on an \name{IntegArray}:
\begin{itemize}
\item Computation of the integrals of zero angular momentum, the $(ss)$ and
      $(ss|ss)$ integrals which cannot be reduced further. This operation 
      is implemented in the constructors: For single particle
      operators, the constructor takes the \name{Coefficient\_Set},
      the auxiliary indices and an appropriate
      \name{Coefficient\_Type} as  arguments. The primary  integrals
      are then either copied from the \name{Coefficient\_Set} or
      computed in the constructor. For two-particle ERI's the
      constructor takes a \name{Quadruple\_Set} and the auxiliary
      index as arguments.
 
\item The evaluation terms in the single-particle recursion relations. 
      As we discussed in class \name{Recursion}, this can be mapped on 
      the following operation: Given a scalar $f$, a \name{CoefArray} 
      $c_{ij}$ and a source \name{IntegArray} $s_{ij}$, 
      the \name{recursion\_term} member function performs:
      \be
       t_{ij} += f\ c_{ij}\ s_{ij} 
      \ee
      for all pairs $ij$. 
      There are seven member functions implemented which perform the
      computational step in the recursion relations for single 
      particle operators, which take into account that 
      the factor $f$ or $c_{ij}$ may be unity or that $1/c_{ij}$ may be 
      needed, too. These are derived from
      the appropriate member functions of \name{Storage}.

      \index{IntegArray,recursion\_term}   \index{recursion\_term}
\item The evaluation terms in the two-particle recursion relations. 
      As we discussed in class \name{VRR}, this can be mapped on 
      the following operation: Consider a scalar $f$, a \name{CoefArray} 
      $c_{IJ}$ and a source \name{IntegArray} $s_{IJ}$, where 
      $I=(i,j)$ and $J=(kl)$ are the compound indices for the first 
      and second \name{RadialFunctionArray}.
      We consider two types of terms: The \name{quadruplet} member function
      implements the relation: \index{IntegArray,quadruplet} \index{quadruplet}
      \be
       t_{IJ} += f\ c_{IJ}\ s_{IJ}
      \ee
      for all pairs $(I,J)$. Here $c$ depends on all four indices. For terms, 
      where $c$ depends only either the compound index $I$ or the compound 
      index $J$, but not on both, it is useful to introduce the member 
      function: \name{doublet} \index{IntegArray,doublet} \index{doublet}
      which implements either:
      \ba
       t_{IJ}  & += & f\ c_{I} \ s_{IJ} \\
       t_{IJ}  & += & f\ c_{J} \ s_{IJ} 
      \ea
      for all $(I,J)$. These operations are derived from
      the appropriate member functions of \name{Storage}.

\item Contractions for two-index and four-index IntegArrays. 
      This operation is accomplished by the constructors, which take a 
      source \name{IntegArray} and either a \name{Coefficient\_Set} or 
      a \name{Quadruple\_Set} as arguments. \index{contractions}
\item Output of  two-index and four-index IntegArrays.
\end{itemize}

\section{class IntegArray}
*/
#ifndef INTEG_ARRAY_H
#define INTEG_ARRAY_H
#include "storage.h"
#include "typesip.h"
#include "coef_array.h"
#include "integ_file.h"

class Quadruple_Set;
class IntegFile;
class BasisSet;


class IntegArray : public Storage
{
  //   static Allocator ialloc;
   int id;
public:
   IntegArray(Allocator& alloc) : Storage(alloc) 
    { zero(); }

   IntegArray(const IntegArray& source)
     : Storage((Storage&) source) 
      { error("ILLEGAL CONSTRUCTOR FOR INTEGARRAY"); } 
   IntegArray(IntegArray& source,Allocator& alloc) 
     : Storage((Storage&) source,alloc) { }

   IntegArray(IntegArray* source,int m1,int m2,Coefficient_Set& cs,
	      Allocator& alloc);
   IntegArray(IntegArray& source,int m1,int m2,int m3,int m4,
	      Quadruple_Set& qset,Allocator& alloc);

   IntegArray(Coefficient_Set& cset,Coef_Type ct,int aux,
	      Momentum& mu,Location& center,Allocator& alloc); 
   IntegArray(Quadruple_Set& qset,int aux,Allocator& alloc);
  ~IntegArray() { }  
  
void recursion_term(CoefArray& prefactor,IntegArray& term)
{ add_to(1.0,prefactor,term); } 
void recursion_term(double factor,CoefArray& prefactor,
		    IntegArray& term)
{ add_to(factor,prefactor,term); } 
void recursion_term(double factor,IntegArray& term)
{ add_to(factor,term); } 
void recursion_term_inverse(CoefArray& invfactor,IntegArray& term)
{ add_invers_to(1.0,invfactor,term); } 
void recursion_term_inverse(double factor,CoefArray& invfactor,
		    IntegArray& term)
{ add_invers_to(factor,invfactor,term); } 
void quadruplet(double factor,CoefArray& coef,IntegArray& src)
{  add_to(factor,coef,src); }
/*
void* operator new(size_t ts)  
{ StorageUnit su; 
  ialloc.get(su); 
  ((IntegArray*) su.base) -> id = su.id;
  return su.base; 
}
void  operator delete(void* p) 
{ StorageUnit su; su.id = ((IntegArray*) p) -> id; 
  su.base = (double*) p;
  ialloc.hand_back(su); 
}
*/
};

inline ostream& operator<<(ostream& os,const IntegArray& a)
{
  const Storage& s = a;
  return os << s;
}
#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
inline ostream_withassign& operator<<(ostream_withassign& os,
				      const IntegArray& a)
{
  const Storage& s = a;
  (ostream&) os << s; 
  return os;
}
#endif
#endif
