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
% $Id: output1.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Output}
\subsection{Single Electron Integrals}

We must distinguish between $b_1 = b_2$ and $b1 \neq b_2$.  If the two
basis sets are identical, each integral with $l_1 \neq l_2$ will be
computed once in the $(l_1,l_2)$ block. Hence all integrals of such
blocks must be printed. If the momenta $\l_1$ and $\l_2$ are the same
we must pay closer attention. In this case the pairs of momenta
$(m_1,m_2)$ appear twice, they have to be printed only once. If
furthermore $m_1 = m_2$ we have \name{integ(i,j) = integ(j,i)} and we
have to output only the upper triangular matrix of \name{IntegArray}.

If the two basis sets are not identical all integrals must be printed
since no symmetries whatever exist.
*/
#include <assert.h>
#include <math.h>
#include "io_incl.h"
#include "const.h"
#include "function.h"
#include "basis.h"
#include "storage.h"
#include "symdesig.h"

void Storage::output(InternalBasis& b1,SymDesignator& s1,
		     InternalBasis& b2,SymDesignator& s2,
		     int l1,int moffset1,int l2,int moffset2,
                     ARRAY<IntegMat>& m,int mu_index,int same,int m_sym)
{
  int i,j;	    
  IVector& start1 = b1.s_index[l1];
  IVector& start2 = b2.s_index[l2];

  if (same && l1 == l2 && moffset1 > moffset2)
    return;

  int both = 1;  
  if (same && l1 == l2 && moffset1 == moffset2)
    both = 0;

  if(!same) both = 0;

  int global1 = b1.sym_start_index(s1) + b1.min_index();
  int global2 = b2.sym_start_index(s2) + b2.min_index();
/*
  cout << "basis offsets   : " << b1.min_index() << " " << b2.min_index() << endl;
  cout << "symm  offsets   : " << b1.sym_start_index(s1) << " " << b2.sym_start_index(s2)
    << endl;
  
  cout << "momentum offsets: " << moffset1 << " " << moffset2 << endl;
  cout << "global offsets  : " << global1 << " " << global2 << endl;
*/  
  int count = 0;  
  for(i= 0 ; i < size1(); i++)
  {
    for(j=0 ; j < size2(); j++)
    {
       int indx1 = start1[i] + moffset1 + global1;
       int indx2 = start2[j] + moffset2 + global2;
       m[mu_index](indx1,indx2) += base[count]; 
       if((indx1 != indx2) && both)
	   {
	      m[mu_index](indx2,indx1) += m_sym * base[count];
	   }
       count++;
    }
  }
}


