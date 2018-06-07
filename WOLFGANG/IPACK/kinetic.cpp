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
% $id: kinetic.c,v 2.6 1995/10/26 10:07:45 krabbe Exp krabbe $
%
\subsection{Constructor}
\index{KineticRec, constructor}

The constructor assigns the arguments to the local variables
and initializes the f-values according to table \ref{term_flags}. 
Note that the auxiliary operator \name{rec} is
passed to the constructor of Recursion. Since there are no auxialiary 
variables we have $max\_mu = max\_aux = 0$.
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

#include "overlap.h"
#include "kinetic.h"

KineticRec::KineticRec(Coefficient_Set& CSet,int momentum1,int momentum2,
		 Recursion* rec) : Recursion(CSet,momentum1,momentum2,
					     0,0,DummyCenter,rec)
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
\index{KineticRec, ss_integ}

Since these values are required in a number of recursion relations, they are 
stored in the \name{CoefArray} KINETIC and only copied by calling
the constructor of \name{IntegArray} here. Aux (and mu) ought be 0.
*/

void KineticRec::ss_integ(const int aux,const int mu)
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
  fields[idx] = new IntegArray(coef_set,KINETIC,aux,zero,center,alloc_prim);
}
/*TEX
\subsection{Special Logic}
\index{KineticRec, speclogic}
We have to ensure, that the overlap integral $(\va''\,|\,\vb'')$ is
already calculated.
 
*/
void KineticRec::speclogic (const int ltot1,const int ltot2,const int aux,const
			    int mu)
{
  if (ltot1 > ltot2)
  {
    if (ltot1 - 2 > 0)
      auxop -> logic(ltot1-2,ltot2,0,0);
  }
  else 
  {
    if (ltot2 - 2 > 0)
      auxop -> logic(ltot1,ltot2-2,0,0);
  }
}
/*TEX
\subsection{Special Terms}
\index{KineticRec, special}

This finally adds the special term \ref{kinspec} using the
\name{CoefArrays} XI, ZETA1 or ZETA2.

*/
void KineticRec::special(IntegArray* target,const int position,
		      const Direction& dir,
		      const Momentum& a,const Momentum& b,
		      const Momentum& mu,const int aux)
{
  Coef_Type ct = XI;
  if (debug)
    cout << "   --> Special Term 0: " << ct << " ";
  IntegArray* term =  auxop -> find(a,b,aux,mu);
  assert(term != 0);
  target -> recursion_term(coef_set[ct],*term);
  
  ct       = (position) ? ZETA1 : ZETA2;
  
  Momentum a1(a);
  Momentum b1(b);
  
  int fac;
  
  if (position)
  {
    a1.reduce(dir);
    a1.reduce(dir);
    fac = -(a(dir) - 1);
  }
  else 
  {
    b1.reduce(dir);
    b1.reduce(dir);
    fac = -(b(dir) - 1);
  }
  
  if (a1.l() >= 0 && b1.l() >= 0)
  {
    if (debug)
      cout << "   --> Special Term 1: " << ct << " ";
    term     =  auxop -> find(a1,b1,aux,mu);
    assert(term != 0);
    target -> recursion_term(fac,coef_set[ct],*term);
  }  

}



