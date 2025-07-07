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
% $Id: contract.cpp,v 1.2 2008/05/03 06:23:55 igor Exp $
%
\section{Contractions\label{section-contract}}

Contractions, particularly on four-index arrays, consume a substantial
fraction of the overall computational effort of this package.
(Approximately  30\% for standard applications). The development of an
efficient contraction procedure is therefore one key to an efficient
implementation of this package. In the following we outline the
general strategy we employ, which is then applied to the contraction
of two- and four-index \name{IntegArrays}

\subsection{Strategy}

Consider a four-dimensional, rectangular array $A$ of dimensions 
($n_1,n_2,n_3,n_4$), which is to be contracted on all four indices to
yield a four-dimensional rectangular array of dimension:
($m_1,m_2,m_3,m_4$). In the following we assume a {\em standard
ordering} of the indices, such that the linear index of the 
element $A(i_1,i_2,i_3,i_4)$ is given as: 
$$ ((i_1 n_2 + i_2) n_3 + i_3) n_4 + i_4 $$
We now contract on the first index, i.e. given a contraction matrix:
$c_{ji}$ of dimension $m_1 \times n_1$ we form:
\be
     A'(j_1,i_2,i_3,i_4) = \sum_{i_1} c_(j_1,i_1) A(i_1,i_2,i_3,i_4)
     \label{eq:contract1}
\ee
The resulting array is of dimension: $(m_1,n_2,n_3,n_4)$ We could
implement this operation for all four indices. Instead, we have
chosen the following method:

It is apparent that the indices $i_2,i_3,i_4$ are inert under the
transformation. We can therefore introduce a compound index: 
$$K = (i_2 n_3 + i_3) n_4 + i_4$$ and index the array as:
\be
  A(i_1,K) = A(i_1,i_2,i_3,i_4)
\ee
The transformation then takes the form of a matrix multiplication,
\be
  A(j,K) = \sum_i c_{ji} A(i,K) 
\ee
where $0 \le K < (n_2 n_3 n_4)$. In other words the look length has
increased substantially. Unfortunately this method works only for the
first index of the array, hence we use the following additional trick:
We permute the indices on $A'$, such that:
\be
     A''(i_2,i_3,i_4,j_1) = A'(j_1,i_2,i_3,i_4) \label{eq:contract2}
\ee
In terms of out block indices this simply means:
\be
    A''(K,j_1) = A'(j_1,K)
\ee
i.e. we transpose the index $j_1$ with the compound index $K$. After
this operation, the index $i_2$ has become the `leading' index of the
four-dimensional array and we may repeat the entire procedure to
contract this index following the same scheme. After repeating this
step four times we have contracted all the indices {\bf and} restored their
original order !

The reader may further verify that the same method can be applied to
two-dimensional arrays by simply setting $n_3 = n_4 = 1$.

The following function, the heart of this contraction method,
implements the basic step described in eqns.~\ref{eq:contract1} and
\ref{eq:contract2}. The \name{target} is a pointer to an contracted
array, while \name{source} contains the uncontracted values.
\name{mx[i]}, where $0 \le i < 3$ contains on entry the dimenstions of
\name{source} on exit the dimensions of \name{target}. The
\name{RadialFunctionArray} \name{f} encodes the coefficent-matrix and
$N$ is the number of contracted orbitals encoded in $f$
*/

#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include "const.h"
#include "integ_array.h"
#include "quad.h"
#include "function.h"

inline void contract_index(IntegType* target,IntegType* source,int* mx,
		    RadialFunctionArray& f,int N)
{
  IntegType  coef;
  int i,j,m;
  int M   = mx[1]*mx[2]*mx[3];
  int max = M*N;
  
  int i1,i2,i3,i4;
  for(i=0;  i< f.size(); i++)
    for(int k=0; k < f[i].contractions(); k++)
    {
      coef = f[i].coef(k);
      j    = f[i].target(k);
      IntegType* sbase = source + i*M;
      IntegType* tbase = target + j;
      if (fabs(coef-1.0) < delta)
      {
	for(m = 0; m < max; m +=N)
	  tbase[m] +=  *(sbase++);	
      }
      else
      {
#ifdef RS6000
	int     half   = M/2;
	IntegType* tbase2 = tbase + N*half;
	IntegType* sbase2 = sbase + half;
	int max = half * N;
	
	for(m=0;m < max; m+=N)
	{
	  tbase [m] += coef * *(sbase++);
	  tbase2[m] += coef * *(sbase2++);
	}
	if (2*half != M)
	  tbase2[m] += coef * *(sbase2++);
#else 
	for(m = 0; m < max; m +=N)
	  tbase[m] += coef * *(sbase++);
#endif
      }
    }
  
  /* now permute the max-pointers to reflect the new maximum values. */

  mx[0] = mx[1];
  mx[1] = mx[2];
  mx[2] = mx[3];
  mx[3] = N;

}
/*TEX
\subsection{Two Particle Operators}

The following constructor implements a contraction for two-particle
operators on all indices. The source \name{IntegArray} refers to an
uncontracted set of integrals, the target set will be contracted. The
contraction coefficients are determined from the
\name{RadialFunctionArray}s in the \name{Quadruple\_Set} which is
passed as the argument.

First the source \name{IntegArray} is converted into a rectangular
array using the \name{to\_matrix} function. The primitive function are
multiplied with their exponent dependent normalization factors using
the \name{normalize} function of the \name{Storage} class. 
The indices are then contracted in sequence  using the
\name{contract\_index} function described above. Finally the
resulting, contracted array is converted back into symmetrical
storage using \name{from\_matrix}.

Note: This operation requires temporary storage, since the contraction
is performed one index at a time. We use the static member \name{work1}
and \name{work2} of the \name{Storage} class for temporary storage.
This typically results in the largest demand for temporary storage of
the entire program.
*/

inline void zero(IntegType* buf,int sz)
{
  memset((char*) buf,0,sz*sizeof(IntegType));  
} 

void   compare(IntegType* f1,IntegType* f2,int elem)
{
  for(int i = 0; i < elem; i++)
    if (fabs(f1[i]-f2[i]) > 1E-10)
      abort();

}
 
IntegArray::IntegArray(IntegArray& source,int m1,int m2,int m3,int m4,
		       Quadruple_Set& qset,Allocator& alloc) : Storage(alloc) 
{
  Coefficient_Set& cs1 = qset.coef_set(0);
  Coefficient_Set& cs2 = qset.coef_set(1);

  RadialFunctionArray& f1 = cs1.radials(0);
  RadialFunctionArray& f2 = cs1.radials(1);
  RadialFunctionArray& f3 = cs2.radials(0);
  RadialFunctionArray& f4 = cs2.radials(1);
    
  int mx[4];
  mx[0] = f1.size();
  mx[1] = f2.size();
  mx[2] = f3.size();
  mx[3] = f4.size();
  
  // first normalize

  source.normalize(m1,m2,m3,m4,qset);

  // if there is nothing to contract --> skip the entire step
  /*
  if (mx[0] == 1 && mx[1] == 1 && mx[2] == 1 && mx[3] == 1)
  {
    *this = source; 
    timer -> stop(t_contract);
    return;
  }
  */
  // contract on the first index, source -> work1.
  // make sure that there is enough space in work1.

  IntegType* outbase;

  int N = f1.nofun();
  if (N == 1 && mx[0] == 1)
  {
    mx[0]   = mx[1];
    mx[1]   = mx[2];
    mx[2]   = mx[3];
    mx[3]   = N;

    outbase = source.get_base(); 
  }
  else
  {
    int M = mx[1]*mx[2]*mx[3];
  
    if (work1.size() < N*M)
      work1.reset(N*M);
    ::zero(work1.get_base(),N*M);
    contract_index(work1.get_base(),source.get_base(),mx,f1,N); 
    outbase = work1.get_base();
  }

  // contract on the second index, work1 -> source
  // since contractions always reduce the number of functions
  // the is enough space in source.

  N = f2.nofun();
  if (N == 1  && mx[0] == 1)
  {
    mx[0]   = mx[1];
    mx[1]   = mx[2];
    mx[2]   = mx[3];
    mx[3]   = N;
  }
  else
  {
    int M = mx[1]*mx[2]*mx[3];  
    if (work2.size() < N*M)
      work2.reset(N*M);
    ::zero(work2.get_base(),N*M);
    contract_index(work2.get_base(),outbase,mx,f2,N); 
    outbase = work2.get_base();
  }

  // contract on the third index, work2 -> work1 

  N = f3.nofun();
  if (N == 1  && mx[0] == 1 && 0)
  {
    mx[0]   = mx[1];
    mx[1]   = mx[2];
    mx[2]   = mx[3];
    mx[3]   = N;
  }
  else
  {
    int M = mx[1]*mx[2]*mx[3];  
    if (work3.size() < N*M)
      work3.reset(N*M);
    ::zero(work3.get_base(),N*M);
    contract_index(work3.get_base(),outbase,mx,f3,N); 
    outbase = work3.get_base();
  }

  // contract on the fourth index, work1 -> target
  // first resize target
  
  N = f4.nofun();
  assert(mx[1] == size1());
  assert(mx[2] == size2());
  assert(mx[3] == size3());
  assert(N     == size4());
  assert(1     == blocks());
  if (N == 1  && mx[0] == 1 && 0)
  {
    if (get_base() != outbase)
      memcpy(get_base(),outbase,mx[1]*mx[2]*mx[3]*sizeof(IntegType));
  }
  else
  {
    zero();
    contract_index(get_base(),outbase,mx,f4,N);
  }
}
/*TEX
\subsection{Single Particle Operators}

The following constructor implements a contraction of the uncontracted
integrals in the source \name{IntegArray} to the contracted integrals
in the target array.

The indices in the source refer to an uncontracted, primitive basis
function. The resulting \name{IntegArray} stores the contracted
integrals. The contractions are encoded in the \name{RadialFunctionArrays} of
the \name{CoefficientSet}. The dimension of the resulting set is
determined from the number of functions encoded in the
\name{RadialFunctionArray}. 

First the source \name{IntegArray} is converted into a rectangular
array using the \name{to\_matrix} function. The primitive function are
multiplied with their exponent dependent normalization factors and
then contracted on the first and second index using the
\name{contract\_index} function described above.  Finally the
resulting, contracted array is converted back into symmetrical storage
using \name{from\_matrix}.

Note: This operation requires temporary storage, since the contraction
is performed one index at a time. We use the static member \name{work1}
of the \name{Storage} class for temporary storage.  

*/

IntegArray::IntegArray(IntegArray* source,int m1,int m2,Coefficient_Set& cs,
		       Allocator& alloc) : Storage(alloc)
{
  int mx1,mx2,i,j,maxi;

  RadialFunctionArray& f1 = cs.radials(0);
  RadialFunctionArray& f2 = cs.radials(1);
  
  mx1 = f1.nofun();
  mx2 = f2.nofun();

  if (f1.size() != source -> size1())
    cout << " Error 117A: " << f1.size() << " " << f2.size()
	 << " " << source -> size1() << " " << source -> size2() << endl;
  
  assert(f1.size() == source -> size1());
  assert(f2.size() == source -> size2());

  // cout << "source: " << source << endl;
  
  source -> to_matrix(work1);
  
  // required size of temporary: mx1 * source.size2()

  int hi1 = ::max(source->size1(),mx1);  
  int hi2 = ::max(source->size2(),mx2);  
  if (work2.size() < hi1*hi2)
    work2.reset(hi1*hi2);

  // first use the temporary to store the norms on second index
    
  for(j=0; j < source -> size2();j++)
    work2[j] = f2[j].norm(m2);
  
  for(i=0; i < source -> size1();i++)
  {
    IntegType norm = f1[i].norm(m1);
    for(j=0; j < source -> size2();j++)
    {
      work1[i*source -> size2()+j] *= norm * work2[j];
    } 
  }

  // cout << "normalized: " << work1 << endl;
  
  // contraction on index 1

  int mx[4];
  mx[0] = source -> size1();
  mx[1] = source -> size2();
  mx[2] = 1;
  mx[3] = 1;
  
  work2.zero();
  contract_index(work2.get_base(),work1.get_base(),mx,f1,mx1); 

  if (work1.size() < hi1*hi2)
    work1.reset(hi1*hi2);
  work1.zero();
  contract_index(work1.get_base(),work2.get_base(),mx,f2,mx2); 
  
  from_matrix(work1,mx1,mx2,1,1,1);
  // cout << "contracted: " << *this << endl;
  
}
