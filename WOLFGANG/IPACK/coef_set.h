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
\section{Coefficient\_Set}
\subsection{Definition}

The \name{Coefficient\_Set} is the set of all \name{CoefArrays}, which
are required for the evaluation of the recursion relations for single
particle operators. The constructor takes the two
\name{RadialfunctionArrays} which describe the primitive orbitals on
the two indices and a list of locations for the computations of single
particle \name{CoeffcientArrays}. Ultimately a \name{Coefficient\_Set}
is nothing but an array of \name{CoefArrays}.  A reference to the
single-particle \name{CoefArray} of \name{Coef\_Type} \name{ct} is
obtained as \name{cs[ct]}

\index{Coefficient Sets, indexing}
\index{Coefficient Sets, creation}
*/

#ifndef COEF_SET
#define COEF_SET

#include <tarray.h>

class SymOp;

enum Integral_Type{OVLP,KIN,NUC,MOM,E_FLD,ANG,SPINORB,TWOEL};
 
class Coefficient_Set : public Array<CoefArray*>
{
  Allocator alloc;
  int sym;
  RadialFunctionArray& rf1;
  RadialFunctionArray& rf2;
  int noloc;
public:
   Coefficient_Set(RadialFunctionArray& f1,int l1,SymOp& s1,
		   RadialFunctionArray& f2,int l2,SymOp& s2,
		   int sflag,Integral_Type I_Type);
  ~Coefficient_Set();

int symmetry()     const { return sym; } 
CoefArray& operator[](int p) { assert(base[p]); return *base[p]; } 
    
Allocator&          allocator() { return alloc; }  
RadialFunctionArray& radials(const int p)
{
  switch(p)
  {
  case 0: return rf1;
  case 1: return rf2;
  default:
    error("Illegal Position in radials. ");
  }
  return rf1; // dummy return
}

};

#endif




