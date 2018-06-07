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
% $Id: nuclear.cpp,v 1.2 2008/05/03 06:23:55 igor Exp $
%
\subsection{Constructor}
The constructor assigns the arguments to the local variables
and initializes the f-values according to table \ref{term_flags}. 
The charge and the center are initialized. (The later by the 
constructor of the Recursion class.) We set max_mu = 0.

\index{NuclearRec, constructor}
*/
#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include "typesip.h"
#include "function.h"
#include "integ_array.h"
#include "coef_array.h"
#include "nuclear.h"

NuclearRec::NuclearRec(Coefficient_Set& CSet,int momentum1,int momentum2,
		 Location& loc,double chg) :
   Recursion(CSet,momentum1,momentum2,momentum1+momentum2,0,loc,0)
{
  charge = chg;
  
  /* set the f-values */

  f0      =  1;  
  f1      = -1;
  f2      =  0;
  f3      =  1;
  f4      = -1; 

  init_fields();
  logic(l1,l2,0,0,1);
}

/*TEX
\subsection{Primary Integrals}
\index{NuclearRec, ss_integ}  
Calling the constructor of \name{IntegArray} with the flag
NUCLEAR does not simply make a copy of a stored \name{CoefArray},
because the function $F$ in equation \ref{nuc_ss} must be evaluated
for different auxiliary indices $m$. These calculations are done
in the constructor. Finally the result is multipleid with the nuclear charge.
*/

void NuclearRec::ss_integ(const int aux,const int mu)
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
  fields[idx] = new IntegArray(coef_set,NUCLEAR,aux,zero,center,alloc_prim);
  (*fields[idx]) *= charge;
}














