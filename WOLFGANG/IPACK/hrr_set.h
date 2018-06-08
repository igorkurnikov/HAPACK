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
% $Id: hrr_set.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{HRR\_Set}
\subsection{Definition}

The \name{HRR\_Set} is the set of all \name{CoefArrays}, which are
required for the evaluation of the horizontal recursion relations for
two particle operators. The indices of these \name{CoefArrays} run
over {\bf contracted function indices}, instead of {\bf primitive}
orbital indices.  The constructor takes the two \name{LocationArrays}
which describe the centers of the contracted orbitals on the two
indices.
*/
#ifndef HRR_SET_H
#define HRR_SET_H

class HRR_Set : public ARRAY<CoefArray*>
{
  Allocator alloc;
  int       sym;
public:
           HRR_Set(LocationArray& l1,SymOp& s1,LocationArray& l2,SymOp& s2,int sflag);
          ~HRR_Set();

int        symmetry()                const { return sym; } 
CoefArray& operator[](Coef_Type type)      { return *base[type-ABX]; } 
};

#endif


