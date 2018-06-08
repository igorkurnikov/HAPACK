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
% $Id: integ_array.cpp,v 1.2 2008/05/03 06:23:55 igor Exp $
%
\section{Primitive \name{IntegArrays}}

The constructors to IntegArray are used for two purposes: contractions
and the evaluation of (ss) or (ss$|$ss)-type integrals at the bottom of
recursion relatons.

*/
#include <assert.h>
#include <math.h>
#include "parallel.h"
#include "io_incl.h"
#include "vtype.h"
#include "integ_array.h"
#include "coef_array.h"
#include "quad.h"
#include "const.h"
#include "f.h"
#include "function.h"
#include "basis.h"
#include "integ_file.h"
#include "storage.h"

// Allocator IntegArray::ialloc(1,1,1,1,2);

/*TEX
\subsection{Primary Two Electron Integrals}

Given a \name{Quadruple\_Set} of coefficients we are able to construct the 
primitive ERI's of $(ss|ss)$ type. The integrals evaluate as:
\be
(s_as_b|s_cs_d)^{(m)} = 2 \left( \frac{\rho}{\pi} \right)^{1/2} 
(s_a|s_b)  (s_c|s_c) F_m(T)  
\ee
with 
\be
T = \rho (P-Q)^2
\ee

We begin with the computation of $\rho$ and $(P-Q)^2$ and store T in
the \name{IntegArray}. We then compute then $F_m(T)$ for the entire
array. The resulting values are multiplied with the overlap integrals
and the sqrt-prefactor involving $\rho$, which is given in terms of
the $\tilde{\zeta_i}$ as:

\ba
 \rho & = & \frac{\zeta_{12} \zeta_{34}} {\zeta_{12} + \zeta_{34}}   
            \nonumber \\
      & = & \frac {1}{2(\tilde\zeta_{s1} + \tilde\zeta_{s2})}
\ea
*/
IntegArray::IntegArray(Quadruple_Set& qset,int aux,Allocator& alloc) 
    : Storage(alloc)

{
  CoefArray& zeta1 = qset.obtain(ZETASUM,0);
  CoefArray& zeta2 = qset.obtain(ZETASUM,1);

  assert(qset.allocator().size1()  ==  size1());
  assert(qset.allocator().size2()  ==  size2());
  assert(qset.allocator().size3()  ==  size3());
  assert(qset.allocator().size4()  ==  size4());
  assert(qset.allocator().blocks() ==  blocks());
  

  IntegType *px,*py,*pz,*qx,*qy,*qz;
  IntegType *z1,*z2;
 
  z1 = zeta1.get_base();
  z2 = zeta2.get_base();
  px = qset.obtain(PX,0).get_base();
  py = qset.obtain(PY,0).get_base();
  pz = qset.obtain(PZ,0).get_base();

  qx = qset.obtain(PX,1).get_base();
  qy = qset.obtain(PY,1).get_base();
  qz = qset.obtain(PZ,1).get_base();
  
  int mx1 = zeta1.size1();
  int mx2 = zeta1.size2();
  int mx3 = zeta2.size1();  
  int mx4 = zeta2.size2();
  //reset(mx1,mx2,mx3,mx4,1);
  assert(mx1 == size1());
  assert(mx2 == size2());
  assert(mx3 == size3());
  assert(mx4 == size4());
  assert(blocks() == 1);
  
  int i,j;
  int count = 0;
  double dist;
  IntegType* base = get_base();  
  
  for(i=0; i< size12(); i++) 
  for(j=0; j< size34(); j++)
  {
    dist = (px[i]-qx[j])*(px[i]-qx[j])
          +(py[i]-qy[j])*(py[i]-qy[j])
	  +(pz[i]-qz[j])*(pz[i]-qz[j]);
    base[count++] = 0.5 * dist / (z1[i] + z2[j]);
  }

  F(base,base,size(),aux);
  
  IntegType* ovlp1  = qset.obtain(OVERLAP,0).get_base();
  IntegType* ovlp2  = qset.obtain(OVERLAP,1).get_base();
  IntegType  fac     = 2*sqrt(0.5/pi);

  count = 0;
  for(i=0; i< size12(); i++) 
  for(j=0; j< size34(); j++)
    base[count++] *= fac*ovlp1[i]*ovlp2[j] / sqrt(z1[i] + z2[j]);
}
/*TEX
\subsection{Primary Single Particle Integrals}

At the bottom of a recursion relation the $(ss)$-type integrals must be
evaluated.  For the operators OVERLAP and KINETIC the values are
available immediately as \name{CoefArray} and can be copied
element-by-element.

For operators of type NUCLEAR, the values naturally depend on the
nuclear center and the (ss)-integrals are computed from elementary
fields in the coefficient set.

\subsubsection{NUCLEAR primitives}

For NUCLEAR integrals the auxiliary integral of order $m$ at center $C$
is defined as:

\be      
(s_a|V_C|s_b)^{(m)} = 
  \sqrt{\frac{2}{\pi(\zeta_1 + \zeta_2)}} (s_a|s_b)  F_{m}(U)
\ee
with
\be
    U = \frac{(P-C)^2}{2 (\zeta_1 + \zeta_2)}
\ee
We recall that
\be
\zeta_1 + \zeta_2 = \zeta_{s}
\ee

We begin with the computation of the argument to the auxiliary
function F for all elements of the \name{IntegArray}. The function is
then evaluated and in a final step multiplied with the appropriate
prefactors.

*/
IntegArray::IntegArray(Coefficient_Set& cset,Coef_Type ct,
		       int aux,Momentum& mu,Location& center,Allocator& alloc)
     : Storage(alloc)     
{
  assert(cset.allocator().size1()  ==  size1());
  assert(cset.allocator().size2()  ==  size2());
  assert(cset.allocator().size3()  ==  size3());
  assert(cset.allocator().size4()  ==  size4());
  assert(cset.allocator().blocks() ==  blocks());

  switch(ct)
  {
  case NUCLEAR:
  {
    int i;

    int mx1 = cset[ZETASUM].size1();
    int mx2 = cset[ZETASUM].size2();

    assert(mx1 == size1());
    assert(mx2 == size2());
    assert(size3() == 1);
    assert(size4() == 1);
    assert(blocks() == 1);
    
    // first compute (P-C)
    
    IntegType* px  = cset[PX].get_base();
    IntegType* py  = cset[PY].get_base();
    IntegType* pz  = cset[PZ].get_base();
    IntegType* zeta = cset[ZETASUM].get_base();
    IntegType* ovlp = cset[OVERLAP].get_base();
    IntegType* base = get_base();

    for(i=0; i< size(); i++)
    {
      IntegType ddx = px[i] - center(X);
      IntegType ddy = py[i] - center(Y);
      IntegType ddz = pz[i] - center(Z);
      base[i] = 0.5 * (ddx*ddx + ddy*ddy + ddz*ddz) / zeta[i];
    }
      
    F(base,base,size(),aux);

    IntegType fac = sqrt(2.0/pi);
    for(i=0; i< size(); i++)
      base[i] *= fac * ovlp[i] / sqrt(zeta[i]);
    
    //    cout << "NUCLEAR: " << center << endl << *this;
    break;    
  }
    
  default:
  {
    Storage& target = *this;
    Storage& source = cset[ct];
    target = source;
  }
  }    
} 












