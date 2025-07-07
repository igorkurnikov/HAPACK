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
\section{Recursion Relation for Electric Field and Electric Field Gradient
Integrals \label{section-e_field}}
\index{E_FieldRec}
\index{electric field integrals}
\index{electric field gradient integrals}

\subsection{Definition}
This recursion relation implements electric field and electric field
gradient integrals at a given center $C$.
According to table \ref{term_flag} all of the general terms
but no special term are used. Auxiliary scalar indices are used and 
$ 0 <= \mu_{tot} = \mu_x+\mu_y+\mu_z$ <= 2$ is expected. Electric
field integrals are obtained with $\vmu = \one_i$ and
\be
(s_\alpha|V_C(\one_i,m)|s_\beta) =  \frac{1}{\tilde\zeta_s}
 	                  (P_i -C_i)(s_\alpha|V_C(\zero,m+1)|s_\beta) 
\label{efield_rec}
\ee
while electric field gradients integrals 
are obtained with $\vmu = \one_i + \one_j$ and
\ba
(s_\alpha|V_C(\one_i+\one_j,m)|s_\beta) & = &
  - \frac{1}{\tilde\zeta_s}\delta_{ij} (s_\alpha|V_C(\zero,m+1)|s_\beta) 
  \nonumber \\
 & & + \frac{1}{\tilde\zeta_s (P_i -C_i)(P_j -C_j)
           (s_\alpha|V_C(\zero,m+2)|s_\beta) 
\label{efieldgrad_rec}
\ea
Both formulas are combined in the $\cmu$-reduction. 
As the root of the recursion the nuclear attraction
ss-integrals are used.
*/ 

#ifndef EFIELD_H
#define EFIELD_H
#include "rec.h"

class E_FieldRec : public Recursion
{
public:
              E_FieldRec(Coefficient_Set& CSet,int momentum1,int momentum2,
			 Location& loc);
             ~E_FieldRec() { }

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








