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
% $Id: hrr1.cpp,v 1.2 2008/05/03 06:23:55 igor Exp $
%
\subsection{Constructor}
*/
#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include <errorip.h>
#include "typesip.h"
#include "function.h"
#include "integ_array.h"
#include "quad.h"
#include "vrr.h"
#include "hrr_set.h"
#include "hrr1.h"

HRR1::HRR1(VRR& vrr1,HRR_Set& hs,
	   int momentum1,int momentum2,int momentum3,int momentum4) :
     vrr(vrr1), hset(hs), 
     RecursionStorage(momentum1,momentum2,0,0,
		      vrr1.alloc_contract.size1(),
		      vrr1.alloc_contract.size2(),
		      vrr1.alloc_contract.size3(),
		      vrr1.alloc_contract.size4(),
		      Momentum::min(momentum3+momentum4+1)- 
                      Momentum::min(momentum3))
{
  timer -> start(t_hrr1);
  l3         = momentum3;
  l4         = momentum4;
  debug      = false;
  
  assert(l1+l2 == vrr1.l1);
  assert(l3+l4 == vrr1.l2);

  no_l1      = Momentum::min(l1+l2+1);
  init_fields();
  noblocks   = Momentum::min(l3+l4+1) - Momentum::min(l3);

  logic(l1,l2,0);
  timer -> stop(t_hrr1);
}

/*TEX
\subsection{Packing the IntegralArrays of the HRR}

\index{VRR,indexing} The packing mechanism for the \name{HRR1} class
is provided implicitly in the storage class through the block-index.
The function block-copy is used to transfer one block of contracted
integrals from the \name{VRR} to its appropriate position in the
\name{HRR1} block. 

Example: Consider a VRR with $l_1+l_2 = 3$ and $l_3+l_4 = 2$. 
Suppose there are two functions on the first index  $(a,b)$, 
two on the second $(c,d)$. Three on the third (e,f,g) and two on the fourth 
(h,i). Consider a single \name{IntegArray} in the VRR. The first dimension  
of the \name{IntegArray} is 4, since there are 4 compound indices,
corresponding to 
\be
    I = (ac,ad,bc,bd)
\ee
The second dimension of the \name{IntegArray} is 6, since we have the indices:
\be
    J = (eh,ei,fh,fi,gh,gi)
\ee
The linear index of a four-index pair in the \name{VRR} is therefore given as:
\be
{\rm IDX}_1(i,j) = i D_j + j,
\ee 
where $D_j=6$ is the dimension of compound index $J$.
  
Now consider the angular variables. For a given angular momentum in the
first  index, say $\va =
(1,1,0),\vb = (0,0,0)$, there are 10 \name{IntegArrays}, labeled by:
\be
    \vc = [(0,0,0)(0,0,0)], [(1,0,0)(0,0,0)],[(0,1,0)(0,0,0)],\dots,
    [(0,0,2)(0,0,0)]
\ee
all of which are to be packed into a single new
\name{IntegArray}, which is now labeled by: $\vc \otimes I \otimes
J$.  The new \name{IntegArray} has thus dimension:
$(D_{\vc},D_I,D_J)$.  The linear index of an element of the new 
\name{IntegArray} is:
\ba
{\rm  IDX}_2 & = & k D_I D_J + (i D_j + j)  \nonumber \\
             & = & k D_I D_J + {\rm  IDX}_1
\ea
In other words, the new \name{IntegArray} is constructed by
concatenating the old \name{IntegArrays} in sequence, $D_I D_J$ is
nothing but the size of an individual \name{IntegArray} in the \name{VRR}.
*/

void HRR1::ss_integ(const int,const int)
{
 
  Momentum zero(0);


  for(int ltot1=l1; ltot1 <= vrr.l1;ltot1++)
    for(Momentum m1(ltot1); m1(); ++m1)
    {
      if (debug)
	cout << "HRR1::ss_integ builds " << m1 << endl;
      
      int idx = compute_index(m1,zero,0,zero);
      fields[idx] = 0; 

      int count = 0;
      for(int ltot2=l3; ltot2 <= vrr.l2;ltot2++)
	for(Momentum m2(ltot2); m2(); ++m2)
        {
	  int vrr_idx = vrr.compute_index(m1,m2,0,zero);
	  IntegArray* array = vrr.fields[vrr_idx];
	  assert(array != 0);

          if (fields[idx] == 0)
	  {
	    fields[idx] = new IntegArray(alloc_prim);
	    //  IntegArray(array -> size1(),array -> size2(),
	    //		 array -> size3(),array -> size4(),
	    //		 noblocks);
	    assert(array -> size1() == fields[idx] -> size1());
	    assert(array -> size2() == fields[idx] -> size2());
	    assert(array -> size3() == fields[idx] -> size3());
	    assert(array -> size4() == fields[idx] -> size4());
	    assert (noblocks == fields[idx] -> blocks());
	    assert(fields[idx] != 0);	    
	  }
	  fields[idx] -> block_copy(count++,*array,0);	  

	  // deleting old integ-array:
	  delete vrr.fields[vrr_idx];
	  vrr.fields[vrr_idx] = 0;
	}
      if (debug)
	cout << "done: " << *(fields[idx]) << endl;
    }
}
/*TEX
\subsection{Logic\label{section-hrr1-logic}} 

Consider the following application of the {\bf horizontal
recursion relation}, where the total angular momenta of the first and
second active index are shown for the doublets required to compute the
(0,3) angular momentum pair. We begin with the preparation of the
$(l,0)$ \name{IntegArrays}. If the target index on the second position
is $l_2$, $l$ must run from $l_2$ to $l_1 + l_2$, where $l_1$ is
the target index on the $l_1$ position. 

For example, let $l_1 = 4$ and $l_2 = 3$. The set of required angular
momentum shells are visualized int table~\ref{table-hrr1}
\begin{table}
\caption{Angular Momentum Relations in the HRR}
\label{table-hrr1}
\bc
\begin{tabular}{llll}
  (4,0) \\   
  (5,0) & (4,1)  \\
  (6,0) & (5,1)  & (4,2) \\
  (7,0) & (6,1)  & (5,2) & (4,3)
\end{tabular}
\ec
\end{table}
The left column of this table shows the output shells of the
\name{VRR}, the bottom-right entry the target angular shell. The
application of the \name{HRR} can be visualized as computing one
vertex, say (5,1),  from the vertex directly to the left (6,0) and
the one diagonally upwards to the left, i.e. (5,0). 

This suggests a strategy of computation:
We begin the recursion with the lowest angular momentum on the second 
index, i.e. computing the $(4,1)$ pair. After the computation of this
pair, we no longer need the $(4,0)$ \name{IntegArrays}, which may be
deleted. We continue with the computation of the $(5,1)$ doublet,
deleting the superfluous $(5,0)$ \name{IntegArrays}. Finally we compute
the $(4,1)$ doublet, deleting $(5,0)$ and $(4,0)$ afterwards. We then
continue with the next row to compute the $(4,2)$ pair and delete
$(1,1)$ etc, etc.

Note: The above strategy must still be implemented.
*/

void HRR1::logic(const int l1,const int l2,const int)
{
  ss_integ(0,0);

  Momentum zero;

  for(int cl2 = 1; cl2 <= l2; cl2++)
  {
    for(int cl1 = l1 + l2 - cl2; cl1 >= l1; cl1--)
      reduce(cl1,cl2,0);
  }
  
}
/*TEX
\subsection{Reduction}

Always reduce the first index and increase the second:
\be
  (\va,\vb|XX) = (\va+\one_d,\vb-\one_d|XX) 
                 + (A_d - B_d) (\va,\vb-\one_d|XX)
\ee 
Since the evaluation of the first term requires no computation
whatever, we simply copy the \name{IntegArray} using the
copy-constructor. 

Note: The following is a speculative investigation  whether it is
possible to simply copy the pointer:
For each vertex in diagram~\ref{table-hrr1} the first term corresponds 
to the term diagonally up and left. This term, as discussed in the
section on \name{logic}  (\ref{section-hrr1-logic}) is needed in the
compution of only one angular shell. Here we must still decide,
whether within this angular shell the term is used more then once: We
note that in order to procduce a collision, we must take two target pairs
$(\va,\vb|XX)$ and $(\va',\vb'|XX)$, to which we apply the recursion
relations in two different directions (without loss of generality we
consider X and Y)  to obtain the same source-pair:
\ba
  \va_x + 1 & = & \va_x' \nonumber \\
  \va_y - 1 & = & \va_y' \nonumber \\
  \va_z     & = & \va_z' \nonumber \\
  \vb_x - 1 & = & \vb_x' \nonumber \\
  \vb_y + 1 & = & \vb_y' \nonumber \\
  \vb_z     & = & \vb_z' \nonumber \\
\ea
In this scenario corresponds the first pair is reduces in direction X
the second one in direction Y. Suppose we always reduce on the minimal
non-zero momentum index. (note this is not the present implementation).
....

*/

void HRR1::reduce(const int ltot1,const int ltot2,const int)
{
  // now apply the the recursion formula term by term for all angular momenta
    
  IntegArray* term;
  Coef_Type   ct;
  Momentum    zero(0);

//  if (ltot1 == 1 && ltot2 == 1) debug = 1;
//  else debug = 0;
  
    
  for(Momentum l1(ltot1); l1(); ++l1)
  for(Momentum l2(ltot2); l2(); ++l2)
  {
    // find the target IntegArray

    if (debug) cout << "Computing Field: " << l1 << l2 << endl; 
    IntegArray *target    =  find(l1,l2,0,zero);
    if (target)
    {
      if (debug) cout << "already done ! " << endl;
      return;  // already done !
    }

    Direction dir         =  l2.decrease();

    Momentum l1base(l1);    
    Momentum l2base(l2);
    Momentum l1term(l1);    
    Momentum l2term(l2);

    l1base.induce(dir);
    l2base.reduce(dir);
    l2term.reduce(dir);
    
    IntegArray *source = find(l1base,l2base,0,zero);
    assert(source != 0);
    
    target = new IntegArray(*source,alloc_prim);
    assert(target != 0);
//    cout << "source: " << l1base << " " << l2base << endl << *source
//      << " target: " << endl << *target << endl;
    
    fields[compute_index(l1,l2,0,zero)] = target;

    // term: (a_i-b_i) (a-1b |XX)(m) 
      
    ct       = (Coef_Type) (ABX +  dir);
    if (debug)
      cout << "   --> Term:  " << ct << " Dir: " << dir 
	   << " Factor: " << hset[ct] << endl;
    term     =  find(l1term,l2term,0,zero);
    assert(term != 0);
    target -> doublet(1.0,hset[ct],*term,0);
    if (debug)
      cout << "Result: " << *target << endl;
  }
}
/*TEX
\subsection{Output} 
*/

void  HRR1::output(InternalBasis& b1,InternalBasis& b2,IntegFile& obuf,int same)
{
  error("HRR1::output is no longer supported. use HRR2. ");
  
  Momentum zero(0);

  for(Momentum m1(l1);m1();++m1)
  for(Momentum m2(l2);m2();++m2)
  {
    IntegArray* a = find(m1,m2,0,zero);
    if (a)
    {
      // now run over the momenta which are packed into the IntegArray
      int block = 0;
      for(Momentum m3(l3); m3(); ++m3)
      {
	if (debug)
	  cout  << "HRR1::outblock: " << m3 << " starts at: " << block << endl;
	// a -> output(b1,b2,l1,m1.local_index(),l2,m2.local_index(),
	// 	    l3,m3.local_index(),0,0,obuf,block++,same12,same34,samepair);
      }
    }
  }
}
