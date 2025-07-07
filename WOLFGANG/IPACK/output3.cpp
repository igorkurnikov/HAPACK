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
#include "basis.h"
#include "sphere.h"

/*TEX
\subsection{Contraction}

The following function contracts the \name{IntegArrays} with auxiliary
index \name{aux} = 0 for all values of $\vmu$. The uncontracted
\name{IntegArrays} are replaced.  
*/

void  Recursion::contract()
{ 
  timer -> start(t_contract);

  for(Momentum m1(l1);m1();++m1)
    for(Momentum m2(l2);m2();++m2)
      for(int l=0; l <= max_mu; l++)
	for(Momentum mu(l); mu(); ++mu)
	{
	  int idx = compute_index(m1,m2,0,mu);
	  if (fields[idx])
	  {
            IntegArray& a = *(fields[idx]);
	    IntegArray* contracted = new IntegArray(fields[idx],l1,l2,
						    coef_set,alloc_contract);
	    assert(contracted !=0);
	    delete fields[idx];
	    fields[idx] = contracted;
	  }
	}
  timer -> stop(t_contract);
}

/*TEX
\subsection{Output} 
We devise a function which allows to combine and output all integrals
of a given auxiliary momentum $\vmu$. The integrals are added to
the matrix in both symmetric positions. The contractions are carried
out in the constructor of the \name{IntegArray} class, the transfer of
the integrals to the matrix is done in the \name{output} member
function of \name{IntegArray}.
\index{recursion relations, output}
\index{operators, output}
*/

void  RecursionStorage::output(InternalBasis& b1,int ll1,SymDesignator& s1,
			       InternalBasis& b2,int ll2,SymDesignator& s2,
			       ARRAY<IntegMat>& m,Momentum& mu,int same,int m_sym)
{
  if(!flaglist(Spherical))
  {  
    for(Momentum m1(ll1);m1();++m1)
      for(Momentum m2(ll2);m2();++m2)
      {
	IntegArray* integ = find(m1,m2,0,mu);      
	integ -> output(b1,s1,b2,s2,
			ll1,m1.local_index(),ll2,m2.local_index(),
                        m,mu.index(),same,m_sym);
      }
  }
  else // spherical gaussians
  {
    Sphere sinteg(*this,false,alloc_prim);
    for(int m1=0; m1 < 2*ll1+1;m1++)
	{
       for(int m2=0; m2 < 2*ll2+1;m2++)
	   {
         IntegArray* integ = sinteg.find(m1,m2,0,mu);      
         assert(integ != 0);
         integ -> output(b1,s1,b2,s2,ll1,m1,ll2,m2,m,mu.index(),same,m_sym);      
	  }
	}
  }  
}



