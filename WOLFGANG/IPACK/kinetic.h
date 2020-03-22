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
% $Id: kinetic.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Recursion Relation for Kinetic Energy Integrals
\label{section-kinetic}}
\index{KineticRec}
 \index{kinetic integrals}

 \subsection{Definition}

 In this recursion relation, special terms appear, since auxiliary
 overlap integrals are required. A pointer to a recursion relation of
 type \name{OverlapRec} must therefore be passed as an argument to the
 constructor of \name{KineticRec}

 The \name{speclogic} and \name{special} functions encode the special terms 
 in the recursion relation (for reduction on a)
 \be
 R = \tilde\xi \,(\va\,|\,\vb) - \tilde\zeta_1\ N_i(\va)\ (\va''\,|\,\vb''))
 \label{kinspec}
 \ee 
 These neccessitate four new name{CoefArrays}, named: 
 \ba
 \mbox{KINETIC}  & = & (s_\alpha|K|s_\beta) \nonumber \\
		 & = & \xi\ (3 - 2\xi (A-B)^2)\ (s_\alpha\,|\,s_\beta) \\
 \mbox{\tilde\xi}       & = & 2\xi,                                         \\
 \mbox{\tilde\zeta_1}    & = & \frac{\xi}{\zeta_a} {\rm \ \ \ and} \nonumber \\
 \mbox{\tilde\zeta_2}    & = & \frac{\xi}{\zeta_b} 
 \ea

 \name{mulogic} and \name{mureduce} are not used.
 
  
*/ 

#ifndef KINETIC_H
#define KINETIC_H
#include "rec.h"

class KineticRec : public Recursion
{
public:
              KineticRec(Coefficient_Set& CSet,int momentum1,int momentum2,
		      Recursion* rec);
             ~KineticRec() { }

virtual  void ss_integ  (const int aux,const int mu);
virtual  void speclogic (const int l1,const int l2,const int aux,const
			 int mu);
virtual  void special   (IntegArray* target,const int position,
			 const Direction& dir,
			 const Momentum& a,const Momentum& b,
			 const Momentum& mu,const int aux);
virtual  void mulogic   (const int l1,const int l2,const int aux,const
			 int mu) {}
virtual  void mureduce  (const int l1,const int l2,const int aux,const
                         int mu) {}
};
#endif

