/*TEX
%
% IPACK - QCHEM
% Quantum Chemistry Project: Integral Package
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is propriatary software. Unlicensed use and distribution
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
% $Id: hrr2.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%

\section{HRR2} 

We now perform the HRR on the second index pair. During this operation
the first pair of angular indices is inert, we therefore whish to
pack them into a single \name{IntegArray}, labeled by the angular
momentum on the third index$\vc$. The opperation itself is analogous
to that described in the \name{HRR1} class, which should be
consulted first.

For each pair of desired angular momenta $(\va,\vb)$, we obtain the
\name{IntegArray} from the \name{HRR1}. In this \name{IntegArray}, we
obtain the sub-block, which corresponds to angular momentum $\vc$.
These integrals can be thought of as a matrix, which is labeled by the
radial-function dependent compound index-pair (I,J), of dimension
$D_I$ and $D_J$ respectively.

The new \name{IntegArray} will have dimensions $(N_B,D_I,D_J)$
respectively, where $N_B$ is the number of pairs $(\va,\vb)$. Let
$\mu$ be the label for such an angular momentum pair and $n(\mu) =
0,\dots,N_B-1$ the associated linear index.

We now copy the blocks labeled by $\mu$ one after tho other into the new 
\name{IntegArray}. An integral with label: (i,j) and angular momentum doublet 
$\mu$ thus appears in position:

\be
     n(\mu) * N_I N_J + (i * N_J + i) 
\ee
*/

#ifndef HRR2_H
#define HRR2_H
#include "hrr1.h"

class HRR2 : public RecursionStorage
{
  HRR_Set& hset;
  int l1,l2,l3,l4;
  int noblocks;
protected:
  HRR1& hrr;
public:
            HRR2(HRR1& hrr1,HRR_Set& hs);
virtual ~   HRR2() {  }  

void     logic  (const int l1,const int l2,const int aux);
void     reduce (const int l1,const int l2,const int aux);
void     ss_integ(const int aux,const int mu);
void     ss_integ_spherical();
  
/* these routines are pure dummies, they are NEVER called */

virtual  void speclogic (const int l1,const int l2,const int aux,const
			 int mu) {}
virtual  void special   (IntegArray* target,const int position,const 
			 Direction& dir,
			 const Momentum& a,const Momentum& b,
			 const Momentum& mu,const int aux) {}

void     output(InternalBasis& b1,SymDesignator& s1,
		InternalBasis& b2,SymDesignator& s2,
		InternalBasis& b3,SymDesignator& s3,
		InternalBasis& b4,SymDesignator& s4,
		int same12,int same34,int samepair,
		IntegFile& obuf); 
};
#endif
