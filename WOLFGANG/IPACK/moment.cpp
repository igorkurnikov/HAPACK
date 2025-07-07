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
\subsection{Constructor}
\index{MomentRec, constructor}
The constructor assigns the arguments to the local variables
and initializes the f-values according to table \ref{term_flags}. 
The location is passed to the base class \name{Recursion}.
Since no auxiliary scalar indices are used we call \name{logic} with the
additional flag \name{use_aux = 0}. This is import for the coefficient of the
f3 term (see section \ref{reduce}).

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

#include "moment.h"

MomentRec::MomentRec(Coefficient_Set& CSet,int momentum1,int momentum2,int mu,
		     Location& loc):
		     Recursion(CSet,momentum1,momentum2,0,mu,loc,0)
{
  /* set the f-values */
 
  f0      = 1;  
  f1      = 0;
  f2      = 1;
  f3      = 1;
  f4      = 0; 

  init_fields();
  logic(l1,l2,0,mu,0);
}

/*TEX
\subsection{Primary Integrals}
\index{MomentRec, ss_integ}
According to equation \ref{mom_ss} in the cas of $\vmu = \zero$ we
just habe to copy the the \name{Coefficient\_Set} OVERLAP by calling
the constructor of \name{IntegArray}.
*/

void MomentRec::ss_integ(const int aux,const int mutot)
{
  // only meaningfull for mutot == 0 !

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
  fields[idx] = new IntegArray(coef_set,OVERLAP ,aux,zero,center,alloc_prim);
}

/*TEX
\subsection{$\vmu$-Logic}
\index{MomentRec, mulogic}
This implements the recursion over $s$ functions for $\mu_{tot} > 0$.
First of all we have to ensure that all integrals with lower
$\mu_{tot}$ have been previously be calculated. 
*/
void MomentRec::mulogic (const int ltot1,const int ltot2,const int aux,const
			 int mutot)
{
  if( mutot > 0 )
  {
    mulogic( ltot1,ltot2,0,mutot-1);
    mureduce(ltot1,ltot2,0,mutot);
  }
  else
    ss_integ(0,0);
}

/*TEX
\subsection{$\vmu$-Reduction}
\index{MomentRec, mureduce}

This sums the terms of equation \ref{mom_rec} using the
\name{CoefArrays} PX to PZ and ZETASUM.

*/
void MomentRec::mureduce(const int l1,const int l2,const int aux,
			 const int mutot)
{
  /* l1,l2 and aux are expected to be zero */

  /* create all targets with total momentum mutot */

  for(Momentum mu(mutot);mu();++mu)
  {
    IntegArray* term;
    Coef_Type   ct;
   

    // find the target IntegArray

    if (debug) cout << "Computing Field: "; 
    IntegArray *target    =  find(l1,l2,aux,mu);
    if (target)
    {
      if (debug) cout << "already done ! " << endl;
      return;  // already done !
    }

//    target = new IntegArray(coef_set[ZETASUM].size1(),
//			    coef_set[ZETASUM].size2());
    target = new IntegArray(alloc_prim);
    
    fields[compute_index(l1,l2,aux,mu)] = target;
    assert(target != 0);
    
    Direction dir  = mu.decrease();
    
    Momentum mubase(mu);    
    mubase.reduce(dir);
    Momentum mubaseminus(mubase);
    mubaseminus.reduce(dir);
    
    // term (p_i-c_i)(0|M(mu-1)|0)

    ct       = (Coef_Type) (PX + dir);
    if (debug)
      cout << "   --> Term Moment 1A : " << ct << endl;
    term     =  find(l1,l2,aux,mubase);
    assert(term != 0);
    target -> recursion_term(coef_set[ct],*term);
    
    double cfac = -center(dir);
    if (debug)
      cout << "   --> Term Moment 1B : Coordinate: " << cfac << endl;
    target -> recursion_term(cfac,*term);
    
    if(mubase(dir) > 0)   
    {
      // term zetasum * N_{dir}(mu) (0|M(mu-2)|0)

      ct = (Coef_Type) (ZETASUM);
      if (debug)
	cout << "   --> Term Moment 2: " << ct << " " << mubase << endl;
      
      term     =  find(l1,l2,aux,mubaseminus);
      assert(term != 0);
      
      target -> recursion_term(mubase(dir),coef_set[ct],*term);
    }
  }
}









