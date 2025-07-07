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
% $Id: twoel0.cpp,v 1.2 2008/05/03 06:23:55 igor Exp $
%
\chapter{Computation of Two Electron Repulsion Integrals\label{chapter-twoel}}

\section{Concepts\label{twoel-concepts}}

The evaluation of the recursion relation for two electron repulsion
integrals, ERI's, \index{ERI} proceeds in close analogy to those for
single particle operators. The following chapters dealing with ERI
evaluation require the formulas derived in OS and HGP, the reader should
consult these publications first.

We begin with the formal developments. The goal of the recursion
relation is the evaluation of all ERI's in a shell of angular momenta
$(l_1,l_2,l_3,l_4)$. By this we mean the set of all ERIS, which have 
functions of total angular momentum $l_i$ on the i-th position. The
radial dependence of the orbitals on each index is encoded in a
\name{RadialFunctionArray}. Each function on each index carries
therefore an index which describes its angular momentum $\va,vb,\dots$
and an index describing its radial dependence $i_1,i_2,\dots$. 
As for single-particle operators, we represent the integral-table as
an outer-product of an composite angular index and a composite radial
index, symbolically written as:
\be
      (\va,\vb|\vc,\vd)(i_1,i_2,i_3,i_4)
\ee 
All integrals for with fixed $(\va,\vb|\vc,\vd)$ are again termed an
\name{IntegArray} and a represented as a matrix with compound indices
$I=(i_1,i_2)$ and $J=(i_3,i_4)$. The angular recursion relations
again express relationships between such \name{IntegArrays}.
Individual terms are of the form 
\be
      (\va,\vb|\vc,\vd)(I,J) = f c_{I,J} (\va',\vb'|\vc',\vd')(I,J)
\ee
and the radial-function dependent coefficients are again represented
as \name{CoefArrays}.

Following a trick first introduced by Head-Gordon and Pople, we employ
the following strategy:

\begin{itemize}
\item Beginning with primary ERI's, we compute the
$(l_1+l_2,0|l_3+l_4,0)$ shell, using the {\bf vertical} recursion
relation. This recursion relation is encoded in class \name{VRR}.
\item After the evaluation of this recursion relation, we can  perform
the contractions, since subsequent operations depend only on the
positions of the orbitals and no longer on the exponents.
\item Using the {\bf horizontal} recursion relation on the first two
indices we compute the $(l_1,l_2|l_3+l_4,0)$ shell of ERI's.
The class \name{HRR1} encodes this recursion relation. 
\item Using the {\bf horizontal} recursion relation on the second pair of
indices we compute the $(l_1,l_2|l_3,l_4)$ shell of ERI's.
The class \name{HRR2} encodes this recursion relation. 
\end{itemize}

As is apparent from the definition of the recursion relations in
classes \name{VRR},\name{HRR1} and \name{HRR2}, all steps have only
{\bf two active} angular indices, hence we may use the class
\name{Recursion} to encode the recursion relations.  At this point the
reader should consult the definitions of the individual recursion
classes. Related and required classes are the storage class for the
required coefficient arrays, named \name{Quadruple\_Set} and naturally
the \name{IntegArray} and \name{CoefArray} classes.

\section{Implementation\label{twoel-impl}}

We begin by assembling the \name{Coefficient\_Sets} for all possible
pairs of angular momentum. These arrays are required in the
computation of the \name{Quadruple\_Sets} for the application of the
\name{VRR}. 

We then loop over all unique quadruplets of angular momentum and
evaluate the recursion relations. 

*/
#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include "basis.h"
#include "storage.h"
#include "twoel.h"


extern "C" int getrusage(int ,struct rusage*); 

void twoel(InternalBasisList& blist,IntegFile& obuf)
{
  if (pm -> master())
    cout << " +++ Number of InternalBasisSets: " << blist.size() << endl;
  for(int b1=0; b1 < blist.size(); b1++)
  for(int b2=0; b2 <= b1; b2++)
  for(int b3=0; b3 <= b1; b3++)
  {
    int maxb4 = b3;
    if (b1 == b3)
      maxb4 = b2;
    for(int b4=0; b4 <= maxb4; b4++)
    {
      if (pm -> master())
	cout << " +++ TwoElectronIntegrals called for indices: " 
	     << setw(3) << b1 << " " << setw(3) << b2 << " " 
	     << setw(3) << b3 << " " << setw(3) << b4 << endl;
      TwoElectronIntegrals twoel(*blist[b1],*blist[b2],*blist[b3],*blist[b4],
				 b1 == b2,b3 == b4, b1 == b3 && b2 == b4,obuf);
    }
  }
}

/*
void twoel(InternalBasis& b1,InternalBasis& b2,
	   InternalBasis& b3,InternalBasis& b4,int same12,int same34,
	   int same_pair,IntegFile& obuf)
{  
  int l1,l2,l3,l4;
  struct rusage RUsage;
  SymOp  nosym(1,1,1);
  
  //
  //   Step 1: Compute the coefficient-sets for all pairs of angular momenta 
  //           on the first and on the second pair of indices.
  //
  
  Array2D<Coefficient_Set*> cset12(b1.max_momentum(),b2.max_momentum());
  Array2D<HRR_Set*>         hset12(b1.max_momentum(),b2.max_momentum());
  Array2D<Coefficient_Set*> cset34(b3.max_momentum(),b4.max_momentum());
  Array2D<HRR_Set*>         hset34(b3.max_momentum(),b4.max_momentum());
  cset12.set(0);
  hset12.set(0);
  cset34.set(0);
  hset34.set(0);
  
  Storage dummy;
  dummy.print_stats();
    
  for(l1=0; l1 <  b1.max_momentum(); l1++)
  if (b1.noprim[l1] > 0)
  {
    //
    // the second angular momentum is smaller than the first 
    // if both functions come form the smae basis set.
    //
    int maxl2 = l1 + 1;
    if (!same12)
      maxl2 = b2.max_momentum();

    for(l2=0; l2 < maxl2; l2++) 
    if (b2.noprim[l2] > 0)
    {
      cset12(l1,l2) = new Coefficient_Set(b1.radials[l1],l1,nosym,
					  b2.radials[l2],l2,nosym,TWOEL);
                              // no symmetry: same && (l1 == 0) && (l2 == 0));
      hset12(l1,l2) = new HRR_Set(b1.fcenter[l1],b2.fcenter[l2],false);
                              // no symmetry: same && (l1 == 0) && (l2 == 0));
    }  
  }

  for(l1=0; l1 <  b3.max_momentum(); l1++)
  if (b3.noprim[l1] > 0)
  {
    int maxl2 = l1 + 1;
    if (!same34)
      maxl2 = b4.max_momentum();
    for(l2=0; l2 < maxl2; l2++) // s-functions only 
    if (b4.noprim[l2] > 0)
    {
      cset34(l1,l2) = new Coefficient_Set(b3.radials[l1],l1,nosym,
					  b4.radials[l2],l2,nosym,TWOEL);
                            // no symmetry: same && (l1 == 0) && (l2 == 0));
      hset34(l1,l2) = new HRR_Set(b3.fcenter[l1],b4.fcenter[l2],false);
                            // no symmetry: same && (l1 == 0) && (l2 == 0));  
    }  
  }

  // getrusage(0,&RUsage);	  
  // cout << " Maximum Resident Set Size " << RUsage.ru_maxrss << endl;
  
  // Step 2: Now loop over the non-equvalent pairs of angular momentum
  //           to evaluate the integrals using the recursion relations
 

  for(l1=0; l1 <  b1.max_momentum(); l1++)
  if (b1.noprim[l1] > 0)
  {
    //
    // the second angular momentum is smaller than the first 
    // if both functions come form the smae basis set.
    //
    int maxl2 = l1 + 1;
    if (!same12)
      maxl2 = b2.max_momentum();

    for(l2=0; l2 < maxl2; l2++) 
    if (b2.noprim[l2] > 0)
    {
      //
      // the third angular momentum is smaller than the first
      // if the basis pair 12 equals the basis pair 34
      //
      int maxl3 = b3.max_momentum();
      if (same_pair)
	maxl3 = l1 + 1;
      
      for(l3=0; l3 < maxl3; l3++)   
      if (b3.noprim[l3] > 0)
      {
	// 
	// the last index is doubly constrained:
	// if b3 == b4, the second index l4 < l3
        // if b1 == b3 and b2 == b4 the pair on the left must be larger
        // than that on the right.   
        // 
	int maxl4 = b4.max_momentum();
	if (same34)
	  maxl4 = l3 + 1;
	if (same_pair && l1 == l3)
	  maxl4 = l2 + 1;
	for(l4=0; l4 < maxl4; l4++)   
	if (b2.noprim[l4] > 0)
	{
	}
      }
    }
  }
  cout << "Done. " << endl;
  
  // 
  //   Step 4: Cleanup of Coefficient sets
  //

  for(l1=0; l1 <  b1.max_momentum(); l1++)
  if (b1.noprim[l1] > 0)
  {
    int maxl2 = l1 + 1;
    if (!same12)
      maxl2 = b2.max_momentum();
    for(l2=0; l2 < maxl2; l2++) // s-functions only 
    if (b2.noprim[l2] > 0)
    {
      delete cset12(l1,l2);
      delete hset12(l1,l2);
    }
  }

  for(l1=0; l1 <  b3.max_momentum(); l1++)
  if (b3.noprim[l1] > 0)
  {
    int maxl2 = l1 + 1;
    if (!same34)
      maxl2 = b4.max_momentum();
    for(l2=0; l2 < maxl2; l2++) // s-functions only 
    if (b4.noprim[l2] > 0)
    {
      delete cset34(l1,l2);
      delete hset34(l1,l2); 
    }  
  }


}

*/
