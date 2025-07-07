/*TEX
%
% IPACK - QCHEM
% Quantum Chemistry Project: Integral Package
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is propriatary software. Unlicensed use and distribution
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
% $Id: hrr2.cpp,v 1.2 2008/05/03 06:23:55 igor Exp $
%
\subsection{Constructor}

The constructor is a direct analogon to the constructor of \name{HRR1}.
*/
#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include "typesip.h"
#include "function.h"
#include "integ_array.h"
#include "quad.h"
#include "vrr.h"
#include "hrr2.h"
#include "sphere.h"
#include "basis.h"

int nb_func(int l1,int l2)
{
  if (flaglist(Spherical))
    return (2*l1+1)*(2*l2+1);
  else
    return Momentum::number(l1) * Momentum::number(l2);
}

HRR2::HRR2(HRR1& hrr1,HRR_Set& hs) :
    hrr(hrr1), hset(hs), 
    RecursionStorage(hrr1.l3,hrr1.l4,0,0,
		     hrr1.alloc_prim.size1(),hrr1.alloc_prim.size2(),
		     hrr1.alloc_prim.size3(),hrr1.alloc_prim.size4(),
		     nb_func(hrr1.l1,hrr1.l2))
{
  timer -> start(t_hrr2);
  l1         = hrr.l1;
  l2         = hrr.l2;
  l3         = hrr.l3;
  l4         = hrr.l4;
  debug      = false;

/*  
  max_mu     = 0;
  max_aux    = 0;
*/
  Momentum zero(0);
  no_l1      = zero.min(l3+l4+1);
  init_fields();
  
  logic(l3,l4,0);
  timer -> stop(t_hrr2);
}

/*TEX
\subsection{Packing the IntegralArrays of the HRR}
See class \name{HRR1}. The only difference here is that we unpack the
blocks of the \name{HRR1} class and not the \name{VRR1}. The primary
index of this recursion relation is the {\em hidden, inert index}
of the \name{HRR1}. The primary indices of the \name{HRR1} become the
hidden indices, labeling the blocks of \name{IntegArrays} of the \name{HRR2}.
*/

void HRR2::ss_integ(const int,const int)
{  
  Momentum ref1(l1);
  Momentum ref2(l2);
  Momentum zero(0);
  
  noblocks  = zero.number(l1)*zero.number(l2);
  
  IntegArray* ref = hrr.find(ref1,ref2,0,zero); // this one must be there  
  assert(ref != 0);

  Sphere* sphere = 0;  
  if(flaglist(Spherical))
  {
    noblocks = (2*l1+1)*(2*l2+1);
    sphere   = new Sphere(hrr,false,hrr.alloc_prim);
  }
 
  int source_index =  0;
  for(int ltot3=l3; ltot3 <= l3+l4;ltot3++)
    for(Momentum m1(ltot3); m1(); ++m1)
    {
      if (debug)
	cout << " HRR2: ss_integ: builds " << m1 << " Source-Index: " 
	     << source_index << endl; 
      int idx = compute_index(m1,zero,0,zero);
      fields[idx] = new IntegArray(alloc_prim); 
      // IntegArray(ref -> size1(), ref -> size2(), 
      //	    ref -> size3(), ref -> size4(),noblocks);
      assert(ref -> size1() == fields[idx] -> size1());
      assert(ref -> size2() == fields[idx] -> size2());
      assert(ref -> size3() == fields[idx] -> size3());
      assert(ref -> size4() == fields[idx] -> size4());
      if (noblocks != fields[idx] -> blocks())
      {
	(ostream&) cout << " VALUES: " << l1 << " " << l2 << " " << flaglist(Spherical) << endl <<
	  noblocks << " " << fields[idx] -> blocks() << endl;
	abort();
	
      }
      
      int target_index = 0;
      if (!flaglist(Spherical))
      {
	for(Momentum mm1(l1); mm1(); ++mm1)
	for(Momentum mm2(l2); mm2(); ++mm2)
	{
	  IntegArray* array = hrr.find(mm1,mm2,0,zero);
	  assert(array != 0);
	  fields[idx] -> block_copy(target_index++,*array,source_index);
	}      
      }
      else
      {
	for(int mc1 = 0; mc1 < 2*l1+1; ++mc1)
	for(int mc2 = 0; mc2 < 2*l2+1; ++mc2)
	{
	  IntegArray* array = sphere -> find(mc1,mc2,0,zero);
	  assert(array != 0);
	  fields[idx] -> block_copy(target_index++,*array,source_index);
	}      
      }      
      if (debug)
	cout << "done: " << *fields[idx] << endl;
      source_index++;
    }
  if (sphere)
    delete sphere;
}
  
/*TEX
\subsection{Logic} 
See class \name{HRR1}.
*/

void HRR2::logic(const int l1,const int l2,const int)
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
  (XX|\vc\vd) = (XX|\vc+\one_d,\vd-\one_d) 
                 + (C_d - D_d) (XX|\vc,\vd-\one_d)
\ee 
See class \name{HRR1}.
*/

void HRR2::reduce(const int ltot1,const int ltot2,const int)
{
  IntegArray* term;
  Coef_Type   ct;
  Momentum    zero(0);

  for(Momentum l1(ltot1); l1(); ++l1)
  for(Momentum l2(ltot2); l2(); ++l2)
  {
    // find the target IntegArray

    if (debug) 
      cout << "Computing Field: "; 
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
    
    // target = new IntegArray(*source);
    
    target = new IntegArray(*source,alloc_prim);
    
    assert(target != 0);
    fields[compute_index(l1,l2,0,zero)] = target;
        
    // term: (a_i-b_i) (a-1b |XX)(m) 
      
    ct       = (Coef_Type) (ABX + dir);
    if(debug)
      cout << "   --> Term: " << ct << " ";
    term     =  find(l1term,l2term,0,zero);
    assert(term != 0);

    if (debug)
    {
      cout << "Target: " << *target << endl;
      cout << "Source: " << *term << endl;
      cout << "Coef: " << hset[ct] << endl;
    }
    
    target -> fold(hset[ct],*term,noblocks);
    
  }
}
/*TEX
\subsection{Output} 
*/
void  HRR2::output(InternalBasis& b1,SymDesignator& s1,
		   InternalBasis& b2,SymDesignator& s2,
		   InternalBasis& b3,SymDesignator& s3,
		   InternalBasis& b4,SymDesignator& s4,
		   int same12,int same34,int samepair,
		   IntegFile& obuf)
{
  Momentum zero(0);
  cout << " Writing: " << setw(3) << l1 << setw(3) << l2 << setw(3) << l3 
       << setw(3) << l4 << endl;

  if (!flaglist(Spherical))
  {
    for(Momentum m3(l3);m3();++m3)
    for(Momentum m4(l4);m4();++m4)
    {
      IntegArray* a = find(m3,m4,0,zero);
      if (a)
      {
	// now run over the momenta which are packed into the IntegArray
	  int block = 0;
	for(Momentum m1(l1); m1(); ++m1)
	for(Momentum m2(l2); m2(); ++m2)
	{
	  if (debug)
	    cout  << " HRR2::outblock: " << m1 << m2  
	      << " starts at: " << block << endl;
	  a -> output(b1,s1,b2,s2,b3,s3,b4,s4,
		      l1,m1.local_index(),l2,m2.local_index(),
		      l3,m3.local_index(),l4,m4.local_index(),obuf,block++,
		      same12,same34,samepair);
	}
      }
    }
  }
  else // spherical gaussians
  {
    Sphere sphere(*this,false,alloc_prim);
    for(int m3 = 0; m3 < 2*l3+1;++m3)
    for(int m4 = 0; m4 < 2*l4+1;++m4)
    {
      IntegArray* a = sphere.find(m3,m4,0,zero);
      if (a)
      {
	// now run over the momenta which are packed into the IntegArray
	  int block = 0;
	for(int m1 = 0; m1 < 2*l1+1; ++m1)
	for(int m2 = 0; m2 < 2*l2+1; ++m2)
	{
	  if (debug)
	    cout  << " HRR2::outblock: " << m1 << m2  
	      << " starts at: " << block << endl;
	  a -> output(b1,s1,b2,s2,b3,s3,b4,s4,
		      l1,m1,l2,m2,l3,m3,l4,m4,obuf,block++,
		      same12,same34,samepair);
	}
      }
    }
  }  
}
