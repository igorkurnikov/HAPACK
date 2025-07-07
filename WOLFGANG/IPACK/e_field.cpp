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
\index{E_FieldRec, constructor}
The constructor assigns the arguments to the local variables
and initializes the f-values according to table \ref{term_flags}. 
The location is passed to the base class \name{Recursion}. The maximum
value of $\mu_{tot}$ is $2$ for the electric field gradient. Of course
calling \name{logic} the \name{use_aux} flag must be $1$.

*/
#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include <errorip.h>
#include "typesip.h"
#include "function.h"
#include "integ_array.h"
#include "coef_array.h"
#include "coef_set.h"

#include "e_field.h"

E_FieldRec::E_FieldRec(Coefficient_Set& CSet,int momentum1,int momentum2,
		       Location& loc):
      Recursion(CSet,momentum1,momentum2,momentum1+momentum2+2,2,loc,0)
{
  /* set the f-values */

  f0      = 1;  
  f1      = -1;
  f2      = 1;
  f3      = 1;
  f4      = -1; 

  init_fields();
  logic(l1,l2,0,2,1);
}

/*TEX
\subsection{Primary Integrals}
\index{E_FieldRec, ss_integ}
As mentioned before we just have to create the nuclear attraction
ss-functions.
*/

void E_FieldRec::ss_integ(const int aux,const int mutot)
{
  // only meaningfull for mutot == 0 !

  Momentum zero(0);
  if (debug)
    cout << "  SS-Integ: aux = " << aux << endl;
  int idx = compute_index(zero,zero,aux,zero);
  if (fields[idx]) 
  {
    if (debug)
      cout << "  Already done: " << endl;
    return;  // has already been computed
  }
  fields[idx] = new IntegArray(coef_set,NUCLEAR,aux,zero,center,alloc_prim);
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
\subsection{$\vmu$-Logic}
\index{E_FieldRec, mulogic}
The equations \ref{efield_rec} and \ref{efieldgrad_rec} can be
combined to two recursive steps.
The values $\vmu =\one_i$ or $\one_i+\one_j$ correspond to the 
electric field or
the electric field gradient, respectively. Only these values are taken into
account here, when we ensure, that all source terms of the mentioned equations
are computed before use.

*/
void E_FieldRec::mulogic(const int ltot1,const int ltot2,const int aux,const
			  int mutot)
{
  if( mutot > 0 )
  {
    mulogic(ltot1,ltot2,aux,mutot-1);
    mulogic(ltot1,ltot2,aux+1,mutot-1);
    mureduce(ltot1,ltot2,aux,mutot);
  }
  else
    ss_integ(aux,0);
}

/*TEX
\subsection{$\vmu$-Reduction}
\index{E_FieldRec, mureduce}
Here we only have to implement the righthand term of \ref{efield_rec}
and for the case $\mu_{tot} = 2$ the first righthand term of 
\ref{efieldgrad_rec}. The second term of \ref{efieldgrad_rec} is
automatically obtained by two recursion step of \ref{efield_rec}.
Several special \name{CoefArrays} are used for these calculations (see the
specific chapter on \name{CoefArrays}).

*/

void E_FieldRec::mureduce(const int l1,const int l2,const int aux,
			  const int mutot)
{
  /* l1,l2 and aux are expected to be zero */

  /* create all targets with total momentum mutot */

  int i;

  for(Momentum mu(mutot);mu();++mu)
  {
    IntegArray* term;
    Coef_Type   ct, inv_ct;
   

    // find the target IntegArray

    if (debug) cout << "Computing Field: (" << l1 << "," << l2 << "," <<
	 aux << "," << mu << ")" << endl; 
    IntegArray *target    =  find(l1,l2,aux,mu);
    if (target)
    {
      if (debug) cout << "target already done ! " << endl;
      return;  // already done !
    }

//    target = new IntegArray(coef_set[ZETASUM].size1(),
//			    coef_set[ZETASUM].size2());
    
    target = new IntegArray(alloc_prim);
    
    assert(target != 0);
    fields[compute_index(l1,l2,aux,mu)] = target;

    if(debug)
    {
      cout << "target has been assigned to fields["	
           << compute_index(l1,l2,aux,mu) << "]\n";
      cout << "Contents of fields[" <<compute_index(l1,l2,aux,mu) << "]\n";
      for(i = 0; i < fields[compute_index(l1,l2,aux,mu)]->size();i++)
      {
	cout << "Element " << i << " : " 
             << fields[compute_index(l1,l2,aux,mu)]->base[i] << endl;
      }
    }
    

    Direction dir1  = mu.decrease();
    
    Momentum mubase(mu);    
    mubase.reduce(dir1);

    Direction dir2  = mubase.decrease();

    Momentum mubaseminus(mubase);
    mubaseminus.reduce(dir2);
    
    // term  1/zetasum * (p_i-c_i)(0|A(mu-1)|0)(aux+1)

    ct       = (Coef_Type) (ZETASUM_PX + dir1);
    
    term     =  find(l1,l2,aux+1,mubase);
    assert(term != 0);

    target -> recursion_term(1.0,coef_set[ct],*term);
    
    if (debug)
    {
      cout << "   --> Term E_Field 1A : multiply with " << ct 
           << " mu : " << mubase << endl;
      cout << "   Coef_Array " << ct << endl;
      for(i = 0; i < coef_set[ct].size();i++)
	cout << "     Element " << i << " : " << coef_set[ct].base[i] << endl;
      cout << "   Source IntegArray " << compute_index(l1,l2,aux+1,mubase)
            << endl;
      for(i = 0; i < term->size();i++)
	cout << "     Element " << i << " : " << term->base[i] << endl;
      cout << "    Resulting Target IntegArray " << endl;
      for(i = 0; i < target->size();i++)
	cout << "     Element " << i << " : " << target->base[i] << endl;
    }
    
    double cfac = -center(dir1);

    target -> recursion_term_inverse(cfac,coef_set[ZETASUM],*term);

    if (debug)
    {
      cout << "   --> Term E_Field 1B : "<< cfac << "/ZETASUM " << endl;

      cout << "   Coef_Array " << ZETASUM << endl;
      for(i = 0; i < coef_set[ZETASUM].size();i++)
	cout << "     Element " << i << " : " << coef_set[ZETASUM].base[i] 
             << endl;
      cout << "    Resulting Target IntegArray " << endl;
      for(i = 0; i < target->size();i++)
	cout << "     Element " << i << " : " << target->base[i] << endl;
    }

    
    if(mubaseminus(dir2) >= 0 && dir1 == dir2)   
    {
      // term - 1/zetasum * (0|M(0)|0)(aux+1)

      term     =  find(l1,l2,aux+1,mubaseminus);
      assert(term != 0);
      
      target -> recursion_term_inverse(-1.0,coef_set[ZETASUM],*term);
      if (debug)
      {
	cout << "   --> Term E_Field 2: divide by" << ct 
             << " mu : " << mubaseminus << endl;
        cout << "   Coef_Array " << ct << endl;
        for(i = 0; i < coef_set[ct].size();i++)
          cout << "     Element " << i << " : " << coef_set[ct].base[i] 
               << endl;
        cout << "   Source IntegArray "<<compute_index(l1,l2,aux+1,mubaseminus)
             << endl;
        for(i = 0; i < term->size();i++)
	  cout << "     Element " << i << " : " << term->base[i] << endl;
        cout << "    Resulting Target IntegArray " << endl;
        for(i = 0; i < target->size();i++)
	  cout << "     Element " << i << " : " << target->base[i] << endl;
      }      
    }
  }
}














