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
% $Id: storage.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%

\section{Evaluation of Terms in Recursion Relations}

Here the static counters and work-spaces are defined and intialized.
*/

#include "io_incl.h"
#include "storage.h"
#include "function.h"
#include "quad.h"

NUMARRAY<IntegType> Storage::work1;
NUMARRAY<IntegType> Storage::work2;
NUMARRAY<IntegType> Storage::work3;

long    Storage::c_number  = 0;
long    Storage::c_size    = 0;
long    Storage::c_notot   = 0;
long    Storage::c_totsize = 0;
long    Storage::c_nohigh   = 0;
long    Storage::c_sizehigh = 0;

Storage::Storage(Allocator& al) : alloc(al)
{
  c_stats();
  alloc.get(*this);
}

Storage::Storage(const Storage& scr,Allocator& newalloc) : alloc(newalloc)
{ 
  c_stats(); 
  alloc.get(*this);
  for(int i = 0; i < size(); i++)
    base[i] = scr.base[i];
}
/*
void Storage::IntegTypet(const Storage& coef,const Storage& src,int pos)
{
//  timer -> start(t_doublet);
  int mx12 = size12();
  int mx34 = size34();
  int max  = mx12 * mx34;
  int count = 0,i,j,b;
  
  switch(pos)
  {
  case 0:
    for(b=0; b < alloc.sz[4]; b++)
    {
      for(i=0; i< mx12; i++)
      {
	int max  = mx34;
	double* target1 = base     + count;
	double* source1 = src.base + count;
	double  totfac  = coef.base[i]; 
#ifdef RS6000
	int half = max / 2;
	IntegType* target2 = target1 + half;
	IntegType* source2 = source1 + half;
	for(j= 0; j < half; j++)
	{
	  target1[j] += totfac*source1[j];
	  target2[j] += totfac*source2[j];
	}
	if (2*half != max)
	  target2[j] += totfac*source2[j];
#else
	for(j= 0; j < max; j++)
	  target1[j] += totfac*source1[j];
#endif
	count+=max;
      }
    }
    break;
  case 1:
    for(b=0; b < alloc.sz[4]; b++)
    {
      // int half = mx12/2;
      
      for(j=0; j< mx34; j++)
      {
	IntegType* target1 =     base + j + count;
	IntegType* source1 = src.base + j + count;
	IntegType totfac   = coef.base[j];	  	  
#ifdef RS6000
	IntegType* target2 = target1 + half*mx34;
	IntegType* source2 = source1 + half*mx34;
	max = half*mx34;
	for(i=0; i< max; i+=mx34)
	{
	  target1[i] += totfac*source1[i];
	  target2[i] += totfac*source2[i];
	}
	if (2*half != mx12)
	  target2[i] += totfac*source2[i];
#else
	for(i=0; i< max; i+=mx34)
	  target1[i] += totfac*source1[i];
#endif
      }
      count += mx12*mx34;
    }
    break;
  default:
    error("Illegal position in Storage::doublet");
  }
//  timer -> stop(t_doublet);
}
*/


/*TEX
\subsection{Extracting One Block from A Storage Array. }
 
*/

Storage::Storage(const Storage& src,int block_no,Allocator& al) : alloc(al) 
{
  c_stats();
  assert(src.size1() == al.size1());
  assert(src.size2() == al.size2());
  assert(src.size3() == al.size3());
  assert(src.size4() == al.size4());
  assert(al.blocks() == 1);
 
  al.get(*this);
  
  // reset(src.sz[0],src.sz[1],src.sz[2],src.sz[3],1);
  assert(block_no < src.blocks());
  
  IntegType* start = src.base + block_no * size();
  memcpy(base,start,size()*sizeof(IntegType));
}


/*TEX
\subsection{Single Particle Operators}

We have formulated the recursion relations for SP operators in such a
way that the individual terms can be understood as the
element-by-element product of a set of integrals with a set of
coefficients and possibly a constant scalar factor:
\be
  (a_i|O|b_j) +=  f c_{ij} (a'_i|O|b'_j)
\ee
Here $a,a',b,b'$ designate angular momenta and $i,j$ radial functions.
We provide three functions to carry out this operation, depending which type 
of prefactor is required. 

\subsection{Two Particle Operators I}

The recursion relations for two-particle operators have been
implemented as for single-particle operators. Again the underlying
operation can be understood as an element-by-element product of a
four-dimensional \name{IntegArray} with a four-dimensional
\name{CoefArray} and possibly a constant factor.The result is then
added to the target \name{IntegArray}: 
\be
  (a_ib_j|V|c_kd_l) +=  \ f \ c_{ijkl} (a'_ib'_j|V|c'_kd'_l)
\ee
Here $a,\dots,d,a',\dots,d'$ designate angular momenta and $i,j,k,l$ 
radial functions. This operation, in direct analogy to \name{recursion\_term},
is also implemented by the \name{add\_to} function. For some calculation
we need to multiply with $1/c_{ijkl}$, which is done by \name{add_invers_to}.

*/

void Storage::add_to(IntegType fac,const Storage& src)
{
//  timer -> start(t_addto);
  int mx = size();
  for(int i=0; i< mx; i++)
    base[i] += fac*src.base[i];
//  timer -> stop(t_addto);
}

void Storage::add_to(IntegType fac,const Storage& coef,const Storage& src)
{
//  timer -> start(t_addto);
  int mx = size();
  for(int i=0; i< mx; i++)
    base[i] += fac*coef.base[i]*src.base[i];
//  timer -> stop(t_addto);
}


void Storage::add_invers_to(IntegType fac,const Storage& coef,const Storage& src)
{
//  timer -> start(t_addto);
  int mx = size();
  for(int i=0; i< mx; i++)
    base[i] += fac*src.base[i]/coef.base[i];
//  timer -> stop(t_addto);
}

/*TEX
\subsection{Two Particle Operators II}
For some terms, however the prefactor array depends only on the
indices $ij$ or $kl$ respectively, but not on both. For these terms it
is useful to make use of this simplification in the implementation of the 
recursion term. This reduces the number of four-dimensional 
coefficient arrays which we have to precompute and store.

The function \name{doublet} implements this recursion on two-electron
operators. The index pairs $(ij)$ and $(kl)$ are regarded as compound
indices respectively: $(ij)=I$ and $(kl)=J$. The \name{IntegArray} is
then indexed by $(I,J)$, while the coefficient array is indexed by
either $I$ or $J$ alone. Hence we have in the recursion:

\ba
((ab)_I|V|(cd)_J) &+= &f c_{I} ((a'b')_I|V|(c'd')_J) {\rm \ \ or} \\
((ab)_I|V|(cd)_J) &+= &f c_{J} ((a'b')_I|V|(c'd')_J) 
\ea

The argument \name{position = 0,1} distinguishes between the two cases.

Techinical Note: This is one of the most CPU intensive subroutines of
the package, about 30\% of the total CPU time is spent here. Therefore
great care should be taken to optimize this routine for the particular
architecture. 
*/

void Storage::doublet(IntegType fac,const Storage& coef,const Storage& src,
		      int pos)
{
//  timer -> start(t_doublet);
  int mx12 = size12();
  int mx34 = size34();
  int max  = mx12 * mx34;
  int count = 0,i,j,b;
  int half;
  
  switch(pos)
  {
  case 0:
    for(b=0; b < alloc.sz[4]; b++)
    {
      for(i=0; i< mx12; i++)
      {
	int max  = mx34;
	IntegType* target1 = base     + count;
	IntegType* source1 = src.base + count;
	IntegType  totfac  = fac*coef.base[i]; 
#ifdef RS6000
	int half = max / 2;
	IntegType* target2 = target1 + half;
	IntegType* source2 = source1 + half;
	for(j= 0; j < half; j++)
	{
	  target1[j] += totfac*source1[j];
	  target2[j] += totfac*source2[j];
	}
	if (2*half != max)
	  target2[j] += totfac*source2[j];
#else
	for(j= 0; j < max; j++)
	  target1[j] += totfac*source1[j];
#endif
	count+=max;
      }
    }
    break;
  case 1:
    for(b=0; b < alloc.sz[4]; b++)
    {
      half = mx12/2;
      for(j=0; j< mx34; j++)
      {
	IntegType* target1 =     base + j + count;
	IntegType* source1 = src.base + j + count;
	IntegType totfac   = fac*coef.base[j];	  	  
#ifdef RS6000
	IntegType* target2 = target1 + half*mx34;
	IntegType* source2 = source1 + half*mx34;
	max = half*mx34;
	for(i=0; i< max; i+=mx34)
	{
	  target1[i] += totfac*source1[i];
	  target2[i] += totfac*source2[i];
	}
	if (2*half != mx12)
	  target2[i] += totfac*source2[i];
#else
	for(i=0; i< max; i+=mx34)
	  target1[i] += totfac*source1[i];
#endif
      }
      count += mx12*mx34;
    }
    break;
  default:
    error("Illegal position in Storage::doublet");
  }
//  timer -> stop(t_doublet);
}
/*TEX
\subsection{Two-Particle Operators III}

Same as the previous section, but without the scalar prefactor.
This split is acutalyy not that ciritcal, since the scalar prefactors are
resolved outside the innermost loop.

*/
/*TEX
\subsection{Two Particle Operators IV} 
 
For the integral evaluation in the \name{HRR2}, we need to perform the
following operation: We are given a source array of \name{noblock}
blocks, each of which is indexed by two compound indices $I$ and $J$.
$J$ runs over the length of the coefficient array, $I$ is obtained as
the first dimension of \name{src}, divided by the number of blocks.
For each block, we then perform an operation akin to doublet, 
contracting the first index, while changing the second.

Note: This operation has been subsumed in the doublet function.
*/
void Storage::fold(Storage& coef,Storage& src,int noblocks)
{  
  doublet(1.0,coef,src,1);
  return;
}


/*TEX
\subsection{Block Copy}
When packing the blocks in order to form the staorage classes 
for the horizontal recursion relations, we must be able to tranfer one 
block of data, i.e. a full set of integrals from one storage array to 
another. 
*/

void Storage::block_copy(int target,Storage& source,int src_no)
{
  int blk_size = source.size()/source.alloc.sz[4];
  assert(blk_size*(target+1) <= size());
  memcpy(base + blk_size*target,source.base + blk_size*src_no,
	 blk_size*sizeof(IntegType));
}

/*TEX
\subsection{Converting Storage Arrays to Matrices}

In order to perform the contractions we must be able to unpack the
symmetrized storage-arrays into rectangular objects. This is done by
the \name{to\_matrix} function.  This function is presently not very
efficient. However it is needed only for integrals between
s-functions, since all other integrals are non-symmetric. For these it
is quite efficient using the \name{memcpy} routine. 

*/

void Storage::to_matrix(NUMARRAY<IntegType>& rectangle)
{
  int i1,i2,i3,i4;
  // double* array = base;

  if (rectangle.size() < size())
    rectangle.reset(size());
  memcpy(rectangle.get_base(),base,sizeof(IntegType)*size());
}
/*TEX
\subsection{Converting Storage Arrays to Matrices}

The inverse of the \name{to\_matrix} function.  
*/

void Storage::from_matrix(NUMARRAY<IntegType>& rectangle,int siz1,int siz2,
		  int siz3,int siz4,int mxaux)
{
  // reset(size1,size2,size3,size4,mxaux);
  assert(siz1 == size1());
  assert(siz2 == size2());
  assert(siz3 == size3());
  assert(siz4 == size4());
  assert(blocks() == mxaux);
  
  if (rectangle.size() < size())
    rectangle.reset(size());
  memcpy(base,rectangle.get_base(),sizeof(IntegType)*size());
  return;
}


/*TEX
\subsection{Normalization of Primitive Eri's}

The normalization is discussed in the description of the
\name{CoefArray} class. 
*/ 

void Storage::normalize(int m1,int m2,int m3,int m4,Quadruple_Set& qset)
{	       
  Coefficient_Set& cs1 = qset.coef_set(0);
  Coefficient_Set& cs2 = qset.coef_set(1);

  CoefArray& norm1 = cs1[NORM];
  doublet(norm1,*this,0);
  CoefArray& norm2 = cs2[NORM];
  doublet(norm2,*this,1);
}

/*TEX
\subsection{Statistics}

This function prints some statistics on the nubmer and size of 
\name{IntegArrays} used.
*/ 

void Storage::print_stats()
{
  cout << "  +++++ Statistics on Memory Usage -- IntegArrays: " << endl;
  cout << "                  Present      Total       Maximal " << endl;
  (ostream&) cout << "       Number: " << setw(10) << c_number << " " 
                  << setw(10) << c_notot  << " " 
		  << setw(13) << c_nohigh << endl;
  (ostream&) cout << "       Byte  : " <<  setw(10) << c_size  << " " 
                  << setw(10) << c_totsize  << " " 
                  << setw(13) << c_sizehigh  << endl;
  cout << endl;
}
