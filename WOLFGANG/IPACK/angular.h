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
\section{Recursion Relation for Angular Momentum Integrals
\label{section-angular-momentum}}
\index{AngluarRec}
\index{angular momentum integrals}

\subsection{Definition}
We now do the recursion relations for orbital angular momentum integrals
at a center $C$. The vector-parameter \name{mu} is used to index the 
three coordinate directions \name{x}, \name{y} and \name{z} of the
angular momentum operator $L_\mu$, i. e.
$mu = (1,0,0), (0,1,0) and (0,0,1)$ respectivlely. This is done in
order to reserve the scalar index \name{aux} as index for derivatives
of the integrals with respect to the coordinates of the nuclei, which
we want to implement in a further developement of the package. 

The following special terms are encoded in \name{speclogic} and
\name{special}
(here given for reduction on a)
\ba
  R & = & \tilde\zeta_1 \{\one_i \times (B-C)\}_{\mu}(\va'|\vb) \nonumber \\
    &   & \tilde\zeta_s \sum_{k = x,y,z} N_k(\vb) \{\one_i \times
          \one_k\}_\mu (\va'|\vb - \one_k) \label{AngSpec} \\
\ea
and the starting integral over $s$ functions 
\be
(s_\alpha|L_\mu|s_\beta) = \tilde\xi\{(A-C)\times(B-C)\}_\mu
(s_\alpha|s_\beta)
\label{AngMured}
\ee
is done by \name{mureduce}, which calls \name{ss_integ} as starting
point. \name{ss_integ} is therefore just the overlap integral 
$(s_\alpha|s_\beta)$. We also need a pointer to a recursion relation
of type \name{OverlapRec} for the special terms, which is passed as
argument to the constructor of \name{AngularRec}.
*/ 

#ifndef ANGULAR_H
#define ANGULAR_H
#include "rec.h"

class AngularRec : public Recursion
{
public:
              AngularRec(Coefficient_Set& CSet,int momentum1,int momentum2,
                      Location& loc,Recursion* rec);
             ~AngularRec() { }

virtual  void ss_integ  (const int aux,const int mu);
virtual  void speclogic (const int l1,const int l2,const int aux,const
			 int mu);
virtual  void special   (IntegArray* target,const int position,
			 const Direction& dir,
			 const Momentum& a,const Momentum& b,
			 const Momentum& mu,const int aux);
virtual  void mulogic   (const int l1,const int l2,const int aux,const
			 int mu);
virtual  void mureduce  (const int l1,const int l2,const int aux,const
                         int mu);
};
#endif



