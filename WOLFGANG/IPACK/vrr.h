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
% $Id: vrr.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Two Electron Vertical Recursion Relation}

\subsection{Concepts}

The {\bf vertical recursion relation} \index{recursion relation,
vertical} \index{VRR} is used in the evaluation of ERI's with angular
momentum multiplet, where one angular momentum index on each side is
zero, i.e. $\vb = \vd = 0$. Hence we have only two {\em active}
angular momentum indices and we may use the framework developed in the
class \name{Recursion} to deal with angular indices. The
\name{IntegArrays} used in the \name{VRR}, however, have two compound
indices, each of which label a pair of \name{Radialfunctions}.

According to OS, the recursion relation is given as:
\ba
    [\,\va+\one_d\, 0\, |\, \vc\, 0\,]^{(m)} & = &
    (P_i-A_i)\ [\,\va\, 0\, |\, \vc\, 0\,]^{(m)} +
    (W_i-P_i)\ [\,\va\, 0\, |\, \vc\, 0\,]^{(m+1)} + \nonumber \\
&&  \frac{N_d(\va)}{2\zeta_{12}} \left( 
    [\,\va-\one_d\, 0\, |\, \vc\, 0\,]^{(m)} -
    \frac{N_d(\va)\rho }{2\zeta_{12}^2} 
    [\,\va-\one_d \,0 \,|\, \vc\, 0\,]^{(m+1)} \right) +
    \nonumber \\    
&&  \frac{N_d(\vc)\rho }{2(\zeta_{12}+\zeta{34})}\  
    \,[\,\va\, 0\, |\, \vc-\one_d\, 0\,]^{(m+1)}
\ea

The evaluation of each term in this recursion relation can again be
viewed as an element by element multiplication of two matrices, one
coefficient-matrix and an \name{IntegArray}, each of which are labeled
by compound indices enumerating  pairs of \name{RadialFunctions}. The
coefficient-matrices  are again represented as \name{CoefArrays} and a
structure to store the required \name{CoefArrays} is provided in
analogy to \name{Coefficient\_Set}. We call this structure a 
\name{Quadruple\_Set} and the reader should consult its definition. 
The function which carries out the multiplication of the
coefficient-matrix and a 4-index \name{CoefArray} is the member
\name{quadruplet} member of \name{IntegArray}. 

In order to save on storage, we have not represented all
coefficient-matrices above in the \name{Quadruple\_Set}. The
coefficients of the first term, e.g. depend only on the compound index
of the first index-pair and not on the second. Furthermore, the 
\name{Coefficient\_Set} stored already the \name{PAX,PAY} and
\name{PAZ} fields for each pair of \name{RadialFunctions}. Hence we
provide a function to multiply the two-index \name{CoefArray}, with a
4-index \name{IntegArray}, the \name{doublet} member of
\name{IntegArray}. The first and third term of the recursion relation
can obviously treated in this fashion, as can the second part of the
$W_i-P_i$ term. Hence we need only four-index \name{CoefArrays} for
the following terms:
\ba
     \rho_1 & = & -\frac{\rho}{2 \zeta_{12}^2} \\
     \rho_2 & = & -\frac{\rho}{2 \zeta_{34}^2} \\
     \rho_s & = & \frac{1}{2 (\zeta_{12} + \zeta_{34})} {\rm\ and} \\ 
     \vec{W}& = & \frac{\zeta_{12} \vec{P}_{12} + \zeta_{34} \vec{P}_{34}}
                       {(\zeta_{12} + \zeta_{34})} 
\ea
These are provided by the quadruple-set.

The \name{logic} of the integral evaluation works pretty much as for
recursion, but completely new functions were provided  for
\name{logic},\name{reduce} and \name{ss\_integ}.  The spirit of these
routines, however, is exactly that of their cousins in \name{Recursion}. 
*/
#include "rec.h"

class InternalBasis;

class VRR : public RecursionStorage
{
  Allocator alloc_contract;
  
  int target_l1,target_l2,target_l3,target_l4;
  
  friend class HRR1;
  friend class HRR2;
protected:
  Quadruple_Set&      qset;

public:
            VRR(Quadruple_Set& qs,int momentum1,int momentum2,
		int momentum3,int momentum4);
           ~VRR(); 

void     logic  (const int l1,const int l2,const int aux);
void     reduce (const int l1,const int l2,const int aux, const int position);
//void     reduce (const int l1,const int l2,const int aux);
void     contract(const int l1, const int l2);

void     ss_integ  (const int aux,const int mu = 0);
void     output(InternalBasis& b1,SymDesignator& s1,
		InternalBasis& b2,SymDesignator& s2,
		InternalBasis& b3,SymDesignator& s3,
		InternalBasis& b4,SymDesignator& s4,
		int same12,int same34,int samepair,
		IntegFile& obuf);

Allocator& allocator() { return alloc_contract; }  
/* these routines are pure dummies, they are NEVER called */

virtual  void speclogic (const int l1,const int l2,const int aux,const
			 int mu) {}
virtual  void special   (IntegArray* target,const int position,const 
			 Direction& dir,
			 const Momentum& a,const Momentum& b,
			 const Momentum& mu,const int aux) {}     
};




