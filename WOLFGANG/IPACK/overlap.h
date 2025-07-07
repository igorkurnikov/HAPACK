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
% $Id: overlap.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Recursion Relation for Overlap Integrals\label{section-overlap}}
\index{OverlapRec}
\index{overlap integrals}
  
\subsection{Definition}

We begin with the class for overlaps. In this recursion relation, no
$\mu$ and no special terms appear, hence the \name{mureduce} 
and the \name{special}
functions and their asscociated logic functions are empty.  

*/ 
#include "rec.h"

class OverlapRec : public Recursion
{
public:
              OverlapRec(Coefficient_Set& CSet,int momentum1,int momentum2);
             ~OverlapRec() {}

virtual  void ss_integ  (const int aux,const int mu);
virtual  void speclogic (const int l1,const int l2,const int aux,const
			 int mu) {}
virtual  void special   (IntegArray* target,const int position,const 
			 Direction& dir,
			 const Momentum& a,const Momentum& b,
			 const Momentum& mu,const int aux) {}
virtual  void mulogic   (const int l1,const int l2,const int aux,const
			 int mu) {}
virtual  void mureduce  (const int l1,const int l2,const int aux,const
                         int mu) {}
};








