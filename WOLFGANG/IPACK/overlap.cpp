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
% $Id: overlap.cpp,v 1.2 2008/05/03 06:23:55 igor Exp $
%
\subsection{Constructor}
\index{OverlapRec, constructor}

The constructor merely assigns the arguments to the local variables
and initializes the f-values according to table \ref{term_flags}. 
Since there are no auxialiary variables we have max\_mu = max\_aux = 0.
*/
#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include "typesip.h"
#include "function.h"
#include "integ_array.h"
#include "coef_array.h"
#include "coef_set.h"
#include "rec.h"
#include "overlap.h"

OverlapRec::OverlapRec(Coefficient_Set& CSet,int momentum1,int momentum2) :
     Recursion(CSet,momentum1,momentum2,0,0,DummyCenter,0)
{
  /* set the f-values */

  f0      = 1;  
  f1      = 0;
  f2      = 0;
  f3      = 1;
  f4      = 0; 

  init_fields(); 
  logic(l1,l2,0,0,0);
}

/*TEX
\subsection{Primary Integrals}
\index{OverlapRec, ss_integ}

Since these values are required in a number of recursion relations, they are 
stored in the \name{Coefficient\_Set} OVERLAP and only copied by calling the
constructor of \name{IntegArray} here. Note that 
aux (and mu) ought be 0. 
*/

void OverlapRec::ss_integ(const int aux,const int mu)
{
  Momentum zero(0);
  if (debug)
    cout << "  SS-Integ: ";
  int idx = compute_index(zero,zero,aux,zero);
  if (fields[idx]) 
  {
    if (debug)
      cout << "  Already done: ";
    return;  // has already been computed
  }
  fields[idx] = new IntegArray(coef_set,OVERLAP,aux,zero,center,alloc_prim);
}








