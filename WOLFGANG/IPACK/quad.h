/*TEX
%
% IPACK - QCHEM
% Quantum Chemistry Project: Integral Package
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is propriatary software. Unlicensed use and distribution
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
% $Id: quad.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Quadruple Sets\label{qset}}

\subsection{Definition\label{qset_def}}
\name{Quadruple\_Sets} are the four-index analogue of
\name{Coefficient\_Sets}.  The constructor takes two 
\name{Coefficient\_Sets} and computes all four-index
\name{CoefArrays}, which are required for the evaluation of the ERI
recursion relations. The references to the constituting
\name{Coefficient\_Sets} are stored in the \name{Quadruple\_Set}.

If $\name{qset}$ is a \name{Quadruple\_Set}, we may obtain {\bf both} 
two-index and four-index \name{CoefficientArrays} in the following fashion:
\index{Quadruple\_Set,access}
\index{Quadruple\_Set,obtain}
\index{Quadruple\_Set,initializer}

\begin{itemize}
\item If \name{Coef\_Type} \name{ct} corresponds a to a single particle 
      \name{CoefArray}, then \name{qset.obtain(ct,pos)} returns a reference to 
      the \name{CoefArray} of the required type on the first/second
      index-pair, if \name{pos} is equal to zero or one respectively.
      All other values for \name{pos} are illegal. It is furthermore illegal 
      to call obtain for the four-index \name{Coef\_Type} with the \name{pos}
      parameter.
\item If \name{Coef\_Type} \name{ct} corresponds a to a two particle 
      \name{CoefArray}, then \name{qset.obtain(ct)} returns a reference to 
      the required \name{CoefArray}. It is illegal to call \name{obtain} 
      for a two-index \name{Coef\_Type} without the \name{pos} parameter.
\end{itemize}

The member function \name{coef\_set} allows to retrieve references to the 
constituent \name{Coefficient\_Sets}.

\index{Quadruple\_Set,access to CoefficientSets}

*/
#include "coef_array.h"
#include "coef_set.h"

class Quadruple_Set : public ARRAY<CoefArray*>
{
    Allocator alloc;
    int sym;
    Coefficient_Set& cset1;
    Coefficient_Set& cset2;
public:  
     Quadruple_Set(Coefficient_Set& cs1,Coefficient_Set& cs2,int sflag);
    ~Quadruple_Set();
    
Allocator&       allocator() { return alloc; }  
int symmetry()   const { return sym; } 
CoefArray&       obtain(const Coef_Type& ct,int pos = -1);
CoefArray& operator[](int p) { assert(base[p]); return *base[p]; } 
Coefficient_Set& coef_set(const int index)
{
  switch(index)
  {
  case 0: return cset1;
  case 1: return cset2;
  default:
    error("illegal index in Quadruple_Set::coef_set");
  }
  return cset1; // dummy to suspress compiler warning
}

};

