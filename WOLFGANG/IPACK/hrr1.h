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
% $Id: hrr1.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Two Electron Horizontal Recursion Relation}
\subsection{Concepts}

The horizontal recursion relation is obtained by subtracting the
recursion relations of OS for for: $[\va + \one_d \vb | \vc \vd]^{(m)}$ and 
$[\va \vb + \one_d | \vc \vd]^{(m)}$, resulting in a relation which
allows to increase one angular index  at the expense of another.  This
``trick'' was first noticed by HGP and used to ease the evaluation of
ERI recursion relations.
\be
[\va + \one_d \vb | \vc \vd]^{(m)} - [\va \vb + \one_d | \vc \vd]^{(m)}
  = (B_i - A_i)   [\va  \vb | \vc \vd]^{(m)}
\ee
The following important points are immediately apparent:
\begin{itemize}
\item No auxiliary indices appear.
\item No exponents appear in the recursion relation.
\item The sum of the total angular momenta on all indices is conserved, 
      hence the name {\bf horizontal recursion relation} 
\end{itemize} 

This suggests the following strategy, due to Head-Gordon and Pople,
for the evaluation of ERI's of angular momenta $(l_1,l_2,l_3,l_4)$

\begin{itemize}
\item Use the {\bf vertical recursion relation} to compute the integrals of 
      type $(l_1+l_2,0,l_3+l_4,0)^{(0)}$. 
      All other required integrals appear also in the VRR. 
\item Contract the primitive basis functions immediately after the VRR is 
      complete, since the HRR applies also to contracted basis functions 
      on the same center. 
\item Pack all integrals $(\va,0,\vc,0)$ for a fixed $\va$ into a single
      \name{IntegArray}, labeled by $(\va,\vb = 0)$. This is possible
      since the recursion  relation for the first two indices is
      completely independent of the angular momentum on the last two. 
      Therefore, in order to minimize memory management and maximize
      loop-length, we are able to process all terms for a given
      $(\va,\vb)$ pair at once.
\item Apply the {\bf horizontal recursion relation} to the first two indices,
      increasing the angular momentum on the second index, while reducing 
      the first until the total momentum on the first index is $l_1$ and 
      the total momentum on the second index is $l_2$.
\item Repack all integrals $(\va,\vb|\vc,0)$ for fixed $\vc$ into a single 
      \name{IntegArray}. Label the resulting \name{IntegArray} by 
      $(\vc,\vd =0)$. 
\item Apply the {\bf horizontal recursion relation} on the last two indices,
      increasing the angular momentum on the fourth index, while reducing 
      the third until the total momentum on the first index is $l_3$ and 
      the total momentum on the second index is $l_4$.
\item Output.
\end{itemize}

The first step in this procedure is carried out using the \name{VRR}
class. Here we implement the second step, increasing the $\vb$angular
index to the desired level.

*/
#ifndef HRR1_H
#define HRR1_H
#include "rec.h"
#include "hrr_set.h"

class HRR1 : public RecursionStorage
{
  HRR_Set&     hset;
  friend class HRR2;
  int l3,l4;
  int noblocks;
  
protected:
  VRR& vrr;
public:
            HRR1(VRR& vrr1,HRR_Set& hset,
		 int momentum1,int momentum2,int momentum3,int momentum4);
virtual ~   HRR1() {  }  

void     logic  (const int l1,const int l2,const int aux);
void     reduce (const int l1,const int l2,const int aux);
void     ss_integ(const int aux,const int mu);
void     output(InternalBasis& b1,InternalBasis& b2,IntegFile& obuf,int same);     

/* these routines are pure dummies, they are NEVER called */

virtual  void speclogic (const int l1,const int l2,const int aux,const
			 int mu) {}
virtual  void special   (IntegArray* target,const int position,const 
			 Direction& dir,
			 const Momentum& a,const Momentum& b,
			 const Momentum& mu,const int aux) {}
};
#endif
