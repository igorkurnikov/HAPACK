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
% 
\subsection{Constructor}
\index{AngularRec, constructor}

The constructor assigns merely the arguments to the local variables
and initializes the f-values according to table \ref{term_flags}.
Note that the auxiliary operator \name{rec} is passed the \name{Recursion}
base class. \name{use_aux} must be zero even in the case when the aux
parameter will be used to index derivatives with respect to nuclei
coordinates, because it is not used for the recursion!
We also have to restrict the recursion to $\mu_{tot} = 1$ setting
$mumin = mumax = 1$ when calling \name{logic}.
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
#include "angular.h"

AngularRec::AngularRec(Coefficient_Set& CSet,int momentum1,int momentum2,
                 Location& loc,Recursion* rec) : 
                 Recursion(CSet,momentum1,momentum2,0,1,loc,rec)
{
  /* set the f-values */

  f0      = 1;  
  f1      = 0;
  f2      = 0;
  f3      = 1;
  f4      = 0; 

  init_fields();
  logic(l1,l2,0,1,0,1);
}
/*TEX
\subsection{Primary Integrals}
\index{AngularRec, ss_integ}
Note that aux should be 0. This is just the overlap ss function which
is called by mulogic.
*/
void AngularRec::ss_integ(const int aux,const int mu)
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
  if(debug)
  {
    cout << "Contents of fields[" << idx << "] :\n";
    for(int i = 0; i < fields[idx]->size();i++)
    {
      cout << "Element " << i << " : " << fields[idx]->base[i] << endl;
    }
  }
}
/*TEX
\subsection{Special Logic}
\index{AngularRec, speclogic} 
We have to make sure, that the overlap integrals required in
\ref{AngSpec} are provided.
*/

void AngularRec::speclogic (const int ltot1,const int ltot2,const int aux,const
			 int mu)
{
  if(ltot1 > ltot2)
  {
    if(ltot1 - 1 >= 0)
    {
      auxop -> logic(ltot1-1,ltot2,0,0);
      if(ltot2 - 1 >= 0)
	auxop -> logic(ltot1-1,ltot2-1,0,0);
    }
  }
  else
  {  
    if(ltot2 - 1 >= 0)
    {
      auxop -> logic(ltot1,ltot2-1,0,0);
      if(ltot1 - 1 >= 0)
	auxop -> logic(ltot1-1,ltot2-1,0,0);
    }
  }
}
/*TEX
\subsection{Special Terms}
\index{AngularRec, special}
This finally adds the special terms of \ref{AngSpec} using
several specific \name{CoefArrays}.
Here we should mention a point
which is not stated in the paper of \cite{OS}: Since the operator is
defined without the imaginary $i$ we have for the special terms an 
sign change according to the position of reduction. 
 */

void AngularRec::special(IntegArray* target,const int position,
		      const Direction& dir,
		      const Momentum& a,const Momentum& b,
		      const Momentum& mu,const int aux)
{
  if(debug)
  {
    cout << " *** SPECIAL ANGULAR *** a: " 
         << a << " b: " << b << " aux: " << aux 
         << " mu: " << mu << endl;
    cout << "    Starting Target IntegArray " << endl;
    for(int i = 0; i < target->size();i++)
      cout << "     Element " << i << " : " << target->base[i] << endl;
  }

  int i;

  Momentum zero(0);
  Momentum a1(a);
  Momentum b1(b);
  int sign;
  

  if (position)
  {
    a1.reduce(dir);
    sign = 1;
  }
  else 
  {
    b1.reduce(dir);
    sign = -1;
  }

  if( a1.l() >= 0 && b1.l() >= 0 )
  {
    Direction mudir;
    mudir = mu.decrease();
    Coef_Type ct0 = (position) ? (Coef_Type) (ZETA1)
                      : (Coef_Type) (ZETA2);

      IntegArray* term0 = auxop -> find(a1,b1,aux,zero);
    assert(term0 != 0);

    int fac;
    IntegArray* term1;

#ifdef GNU
    for(Direction k = X; k <= Z; op_plus(k))
#else
    for(Direction k = X; k <= Z; ++k)
#endif
    {
      Coef_Type ct1 = (position) ? (Coef_Type) (ZETA1_BX + k)
	                : (Coef_Type) (ZETA2_AX+ k);

      target ->recursion_term(sign*EpsOp(dir,k,mudir),coef_set[ct1],*term0);

      if (debug)
      {
        // reset values for correct output
	if(position)
	{    
	  b1 = b;
	}
	else
	{
	  a1 = a;
	}

	cout << "   --> Term 1A : multiply with " << ct1 << " " << sign
	     << " EpsOp(" << dir << "," << k << "," << mudir<< "): " 
             << EpsOp(dir,k,mudir) << endl;
	cout << "       Coef_Array " << ct1 << endl;
	for(i = 0; i < coef_set[ct1].size();i++)
	  cout << "        Element " << i<< " : " << coef_set[ct1].base[i] 
               << endl;
	cout << "       Source IntegArray " << endl
	     << "        a1: " << a1 << " b1:" << b1 << " aux: " << aux 
             << " mu : " << zero 
             << " INDEX (OVERLAP): "<< compute_index(a1,b1,aux,zero)
             << endl;
	for(i = 0; i < term0->size();i++)
	  cout << "        Element " << i << " : " << term0->base[i] << endl;
	cout << "    Resulting Target IntegArray " << endl;
	for(i = 0; i < target->size();i++)
	  cout << "     Element " << i << " : " << target->base[i] << endl;
      }

      target ->recursion_term(-sign*EpsOp(dir,k,mudir)*center(k),coef_set[ct0],
                               *term0);

      if (debug)
      {
	cout << "   --> Term 1B : multiply with " << ct0 << " " << sign 
	     << " -EpsOp(" << dir << "," << k << "," << mudir<< "): " 
             << -EpsOp(dir,k,mudir) << " C(" << k << ") : " 
             << center(k) << endl;
	cout << "       Coef_Array " << ct0 << endl;
	for(i = 0; i < coef_set[ct0].size();i++)
	  cout << "        Element " << i<< " : " << coef_set[ct0].base[i] 
               << endl;
	cout << "       Source IntegArray " << endl
	     << "        a1: " << a1 << " b1:" << b1 << " aux: " << aux 
             << " mu : " << zero 
             << " INDEX (OVERLAP): "<< compute_index(a1,b1,aux,zero)
             << endl;
	for(i = 0; i < term0->size();i++)
	  cout << "        Element " << i << " : " << term0->base[i] << endl;
	cout << "    Resulting Target IntegArray " << endl;
	for(i = 0; i < target->size();i++)
	  cout << "     Element " << i << " : " << target->base[i] << endl;
      }

      if(position)
      {    
	b1 = b;
        fac = b1(k);
	b1.reduce(k);    
      }
      else
      {
	a1 = a;
        fac = a1(k);
	a1.reduce(k);
      }
    
      if( a1(k) >= 0 && b1(k) >= 0 )
      {
	term1     =  auxop -> find(a1,b1,aux,zero);
        assert(term1 != 0);
	target -> recursion_term(sign*EpsOp(dir,k,mudir)*fac,coef_set[ZETASUM],
				 *term1);

	if (debug)
	{
	  cout << "   --> Term 2 : multiply with" 
	       << " EpsOp(" << dir << "," << k << "," << mudir<< "): " 
	       << EpsOp(dir,k,mudir) 
	       << " N(" << k << "): " << fac << endl;
	  cout << "       Source IntegArray " << endl
	       << "        a1: " << a1 << " b1: " << b1 << " aux: " << aux 
               << " mu: " << zero 
               << " INDEX (OVERLAP) " << compute_index(a1,b1,aux,zero)
               <<  endl;
          for(i = 0; i < term1->size();i++)
	    cout << "        Element " << i << " : " << term1->base[i] << endl;
	  cout << "    Resulting Target IntegArray " << endl;
	  for(i = 0; i < target->size();i++)
	    cout << "     Element " << i << " : " << target->base[i] << endl;
	}
      }
    }
  }  
}

/*TEX
\subsection{Moment Logic}
\index{AngularRec, mulogic}
We Just provide the overlap ss integral and call mureduce to calculate the
appropriate coefficients.
*/
void AngularRec::mulogic (const int ltot1,const int ltot2,const int aux,const
			 int mu)
{
  ss_integ(0,0);
  mureduce(ltot1,ltot2,0,1);
}

/*TEX
\subsection{Moment Reduction}
\index{AngularRec, mureduce}
Now only the coefficients due to \ref{AngMured} have to be collected.
*/
void AngularRec::mureduce(const int l1,const int l2,const int aux,
			 const int mutot)
{
  /* l1,l2 and aux are expected to be zero */

  /* mutot is expected to be one */

  int i;
  Momentum mu(mutot);

  for(mu ;mu();++mu)
  {
    if(debug)
      cout << " *** MUREDUCE ANGULAR *** l1: " << l1 << " l2: " << l2 
           << " aux: " << aux << " mu: " << mu << endl;

    IntegArray* term;
    Coef_Type   ct;
    Momentum zero(0);

    // find the target IntegArray

    if (debug) 
      cout << "Computing Target Field: "<< compute_index(l1,l2,aux,mu)
           << endl;

    IntegArray *target    =  find(l1,l2,aux,mu);
    if (target)
    {
      if (debug) cout << "Target already done ! " << endl;
      return;  // already done !
    }

//    target = new IntegArray(coef_set[ZETASUM].size1(),
//			    coef_set[ZETASUM].size2());
    
    target = new IntegArray(alloc_prim);
    
    fields[compute_index(l1,l2,aux,mu)] = target;
    assert(target != 0);
    
    if(debug)
    {
      cout << "   Contents of fields[" <<compute_index(l1,l2,aux,mu) << "]\n";
      for(i = 0; i < fields[compute_index(l1,l2,aux,mu)]->size();i++)
      {
	cout << "    Element " << i << " : " 
             << fields[compute_index(l1,l2,aux,mu)]->base[i] << endl;
      }
    }

    // term XI [(A-C) X (B-C)]_i (0|0)

    Direction dir, dir1, dir2;
    int sign;
    
    dir  = mu.decrease();

    switch(dir)
    {
    case X:
      dir1 = Y;
      dir2 = Z;
      sign = 1;
      break;
    case Y:
      dir1 = X;
      dir2 = Z;
      sign = -1;
      break;
    case Z:
      dir1 = X;
      dir2 = Y;
      sign = 1;
      break;
    }  


    ct       = (Coef_Type) (XI_AxBX + dir);
    double cfac;
    
    term     =  find(l1,l2,aux,zero);
    assert(term != 0);
    target -> recursion_term(sign,coef_set[ct],*term);

    if (debug)
    {
      cout << "   --> Term 1 :  multiply with " << ct << endl;
      cout << "       Coef_Array " << ct << endl;
      for(i = 0; i < coef_set[ct].size();i++)
	cout << "        Element " << i << " : " << coef_set[ct].base[i] 
             << endl;
      cout << "       Source IntegArray "<< endl
           << "         l1: " << l1 << " l2: " << l2
           << " aux: " << aux << " mu: " << zero 
           << " INDEX: " << compute_index(l1,l2,aux,zero)
           << endl;
      for(i = 0; i < term->size();i++)
	cout << "        Element " << i << " : " << term->base[i] << endl;
      cout << "    Resulting Target IntegArray " << endl;
      for(i = 0; i < target->size();i++)
	cout << "     Element " << i << " : " << target->base[i] << endl;
    }
    
    ct       = (Coef_Type) (XI_AminusBX + dir1);
    cfac     = center(dir2);    
    target -> recursion_term(-sign*cfac,coef_set[ct],*term);

    if (debug)
    {
      cout << "   --> Term 2 : multiply with "<< ct  
           << "  C(" << dir2 << "): " << cfac << endl;
      cout << "       Coef_Array " << ct << endl;
      for(i = 0; i < coef_set[ct].size();i++)
	cout << "        Element " << i << " : " << coef_set[ct].base[i] 
             << endl;
      cout << "       Source IntegArray "<< endl
           << "         l1: " << l1 << " l2: " << l2
           << " aux: " << aux << " mu: " << zero 
           << " INDEX: " << compute_index(l1,l2,aux,zero)
           << endl;
      for(i = 0; i < term->size();i++)
	cout << "        Element " << i << " : " << term->base[i] << endl;
      cout << "    Resulting Target IntegArray " << endl;
      for(i = 0; i < target->size();i++)
	cout << "     Element " << i << " : " << target->base[i] << endl;
    }

    ct       = (Coef_Type) (XI_AminusBX + dir2);
    cfac     = center(dir1);    
    target -> recursion_term(sign*cfac,coef_set[ct],*term);

    if (debug)
    {
      cout << "   --> Term 3 : multiply with "<< ct 
           << " C(" << dir1 << "): " << cfac << endl;
      cout << "       Coef_Array " << ct << endl;
      for(i = 0; i < coef_set[ct].size();i++)
	cout << "        Element " << i << " : " << coef_set[ct].base[i] 
             << endl;
      cout << "       Source IntegArray "<< endl
           << "         l1: " << l1 << " l2: " << l2
           << " aux: " << aux << " mu: " << zero 
           << " INDEX: " << compute_index(l1,l2,aux,zero)
           << endl;
      for(i = 0; i < term->size();i++)
	cout << "        Element " << i << " : " << term->base[i] << endl;
      cout << "    Resulting Target IntegArray " << endl;
      for(i = 0; i < target->size();i++)
	cout << "     Element " << i << " : " << target->base[i] << endl;
    }
  }
  if(debug)
  {
    cout << " *** FINAL RESULT MUREDUCE *** "<< endl
         << "   Contents of fields[" <<compute_index(l1,l2,aux,mu) << "]\n";
    for(i = 0; i < fields[compute_index(l1,l2,aux,mu)]->size();i++)
    {
      cout << "    Element " << i << " : " 
	   << fields[compute_index(l1,l2,aux,mu)]->base[i] << endl;
    }
    cout << " *** MUREDUCE END *** " << endl;
  }
}











