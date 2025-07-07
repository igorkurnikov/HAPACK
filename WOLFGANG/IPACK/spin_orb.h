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
\section{Recursion Relation for Spin-Orbit Integrals
\label{section-spin-orbit}}
\index{SpinOrbRec}
\index{spin-orbit interaction integrals}

\subsection{Definition}
Similar to the recursion relations for angular momentum integrals
we use the vector-parameter \name{mu} to index the three coordinate
directions of the spin-orbit operator $S_\mu$ ($\mu = x,y,z)$, which
is evaluated at a center $C$. Therefore 
\be
(s_\alpha|S_\mu(m)|s_\beta) = 4\zeta_a \zeta_b \{(A-C)\times(B-C)\}_\mu
(s_\alpha|V_C(\zero,m+1)|s_\beta)
\label{SpinMured}
\ee
has again to be evaluated through \name{mureduce}, which reduces to 
the nuclear attraction ss functions given in \name{ss_integ}.
There are special terms (here given for reduction on a), too:
\ba
  R & = & 2\zeta_b \{\one_i \times (B-C)\}_{\mu}(\va'|V_C(\zero,m+1)|\vb) 
         \nonumber \\
    &   & \sum_{k = x,y,z} N_k(\vb) \{\one_i \times
          \one_k\}_\mu (\va'|V_C(\zero,m+1)|\vb - \one_k). \label{SpinSpec} \\
\ea
As seen from the above formulas, the \name{aux} parameter is also used.
*/ 

#ifndef SPINORB_H
#define SPINORB_H
#include "rec.h"

class SpinOrbRec : public Recursion
{
public:
              SpinOrbRec(Coefficient_Set& CSet,int momentum1,int momentum2,
                      Location& loc,Recursion* rec);
             ~SpinOrbRec() { }

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


