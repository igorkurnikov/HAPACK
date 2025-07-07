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
\section{Recursion Relation for Moment Integrals
\label{section-moment}}
\index{MomentRec}
\index{moment integrals}

\subsection{Definition}

This recursion relation implements all moment integrals at a given
center $C$ up to an arbitrary total $\mu_{tot} = \mu_x + \mu_y + \mu_z\$.

Special terms and auxiliary scalar indices are not present, but
recursion rules for the moment integrals over $s$ functions 
\begin{eqnarray} 
(s_\alpha|M(\vec{\mu})|s_\beta) & = & 
 	(P_i - C_i) (s_i|M(\vec{\mu}-\one_i)|s_j) \nonumber \\
   & &   + \tilde\zeta_s N_i(\vmu) (s_\alpha|M(\vec{\mu}-2\cdot\one_i)|s_\beta)  
\label{mom_rec}
\end{eqnarray}
with  
\be
(s_\alpha|M(\zero)|s_\beta) = (s_\alpha|s_\beta)
\label{mom_ss}
\ee are required.
*/ 

#ifndef MOMENT_H
#define MOMENT_H
#include "rec.h"

class MomentRec : public Recursion
{
public:
              MomentRec(Coefficient_Set& CSet,int Momentum1,int momentum2,
                     int mutot,Location& loc);
             ~MomentRec() { }

virtual  void ss_integ  (const int aux,const int mutot);
virtual  void speclogic (const int l1,const int l2,const int aux,const
			 int mutot) {};
virtual  void special   (IntegArray* target,const int position,
			 const Direction& dir,
			 const Momentum& a,const Momentum& b,
			 const Momentum& mutot,const int aux) {};
virtual  void mulogic   (const int l1,const int l2,const int aux,const
			 int mutot); 
virtual  void mureduce  (const int l1,const int l2,const int aux,const
                         int mutot); 
};
#endif
