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
% $Id: nuclear.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Recursion Relation for Nuclear Attraction Integrals
\label{section-nuclear}}
\index{NuclearRec}
\index{nuclear attraction integrals}

\subsection{Definition}

This recursion relation implements the nuclear attraction integrals for
ONE nuclear center, WITH the factor for the nuclear charge.  

Special terms do not appear, but auxiliary scalar indices are present.
The primitive integrals, which depend on the location of the nuclear
center, can only in part be computed form primary \name{CoefArrays} as:
\be
(s_\alpha|V_C(\zero,m)|s_\beta) = 2 \sqrt{\frac{\zeta}{\pi}} \ 
           (s_\alpha|s_\beta)\ F_m(U)
\label{nuc_ss}
\ee
with
\be
U = \zeta\ (P-C)^2
\ee
where $C$ is the coordinate of the nuclear center.
Recursion on $\mu$ is neccessary.
*/

#ifndef NUCLEAR_H
#define NUCLEAR_H
#include "rec.h"
class NuclearRec : public Recursion
{
  double charge;
public:
              NuclearRec(Coefficient_Set& CSet,int momentum1,int momentum2,
		      Location& loc,double chg);
             ~NuclearRec() { }
virtual  void ss_integ  (const int aux,const int mu);
virtual  void speclogic (const int l1,const int l2,const int aux,const
			 int mu)  { } 
virtual  void special   (IntegArray* target,const int position,
			 const Direction& dir,
			 const Momentum& a,const Momentum& b,
			 const Momentum& mu,const int aux) { } 
virtual  void mulogic   (const int l1,const int l2,const int aux,const
			 int mu) {}
virtual  void mureduce  (const int l1,const int l2,const int aux,const
                         int mu) {}
};
#endif






