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
% $Id: vrr.cpp,v 1.2 2008/05/03 06:23:55 igor Exp $
%
\subsection{Constructor}

The constructor is very similar to that of single particle recursion
relations. For the subsequent evaluation of the {\bf horizontal
recursion relation} we require all pairs of angular momentum where the
first index runs from zero to $l_1$ while the second runs from zero to
$l_2$ inclusive. For this reason the function logic is called for all
such pairs. This causes little overhead in the actual computation,
since each term is evaluated at most once. Since all these terms must
be stored, however, the memory requirements are quite large for total
angular momentum.  
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
#include "sphere.h"
#include "basis.h"

/*TEX
\subsection{Constructor} 

The allocator of the VRR holds the contracted (!)
\name{IntegArrays}. The size is determined by the \name{nofun}
function of the \name{RadialFunctionArrays} arrays. 


*/
VRR::VRR(Quadruple_Set& qs,int momentum1,int momentum2,
	 int momentum3,int momentum4) :
   qset(qs) ,
     RecursionStorage(momentum1+momentum2,momentum3+momentum4,l1+l2,0,
		      qs.allocator()),
     alloc_contract(qs.coef_set(0).radials(0).nofun(),
		    qs.coef_set(0).radials(1).nofun(),
		    qs.coef_set(1).radials(0).nofun(),
		    qs.coef_set(1).radials(1).nofun(),1)
{
  target_l1 = momentum1;
  target_l2 = momentum2;
  target_l3 = momentum3;
  target_l4 = momentum4;
  
  timer -> start(t_vrr);

//  max_mu     = 0;
    max_aux    = l1+l2;

  // debug = 1;
  
  init_fields();
  // we need all angular momentum pairs ll1,ll2, where 
  
  logic(l1,l2,0);
  timer -> stop(t_vrr);  
}
/*TEX
\subsection{Primary ERI's}

ERI's of type $(ss|ss)$ are computed by the constructor of the 
\name{IntegArray}. 

Note: This is not as efficient as it might be, since ERIS of many
auxiliary indices need to be computed. This should be done all in one
calculation since a number of pre-computed quantities can be reused
for each value of the auxiliary index.  
*/
void VRR::ss_integ(const int aux,const int mu)
{
  Momentum zero(0);
  if (debug)
    cout << " +++ SS-Integ: ";
  int idx = compute_index(zero,zero,aux,zero);
  if (fields[idx]) 
  {
    if (debug)
      cout << "    +++ Already done: " << endl;
    return;  // has already been computed
  }
  fields[idx] = new IntegArray(qset,aux,alloc_prim);
  // cout << "Primitives: " << aux << endl << *fields[idx];

}
/*TEX
\subsection{Logic}

In spirit this function is a close cousin to the logic function for
single-particle recursion relations. We run through all shells of
total angular momentum $L$ and build the individual terms ($l_1 | l_2$)
such that $|l_1 + l_2| = L.$ Once the shell with total momentum $L$ is
completed, we can free the integrals with nonzero auxiliary indices of
total external momentum $L-2$ since these are no longer needed.

Note: not implemented. we can also immediately contract the integrals
of total momentum $L-2$.
*/    

void VRR::logic(const int ltot1,const int ltot2,const int)
{

  if(debug)
    cout << "logic for: " << ltot1 << "  " << ltot2 <<endl;

  int m, l1, l2, lsum;
  
  for(m = ltot1+ltot2; m >= 0; m--)
  {
    if(ltot1 >= ltot2)
    {
      if(m>=ltot2)
      {
	for(l1 = 0; l1 <= ltot1+ltot2-m; l1++)
        {
	  reduce(l1,0,m,0);
	  if((l1 > 1) && (m < ltot1+ltot2))
	  {
	    clean(l1-2,0,m+1);
	    if(m == 0)
	      if((l1-2 < target_l1) || (0 < target_l3))
		clean(l1-2,0,0);
	      else
		contract(l1-2,0);
	  }
	}
	if(m < ltot1+ltot2)
          clean(ltot1+ltot2-m-1,0,m+1);
	if(m == 0)
	{
          if(ltot1 > 0)
	    if((ltot1-1 < target_l1) || ( 0 < target_l3))
	      clean(ltot1-1,0,0);
	    else
	      contract(ltot1-1,0);
	  contract(ltot1,0);
	}
      }
      else
      {
	for(lsum = 0; lsum <= ltot1+ltot2-m; lsum++)
	{
	  for(l1 = 0; l1 <= lsum; l1++)
	  {
	    l2 = lsum -l1;
	    if(l1 >ltot1) continue;
	    
	    if(l2 > ltot2-m)
	      continue;

            if((l2 > ltot2-1-m) && (l1 < target_l1-m))
	      continue;

            if(l2 == 0)
	      reduce(l1,0,m,0);
	    else
	      reduce(l1,l2,m,1);

            if((l1 > 1) && (l2 < ltot2-m))
            {
              if((l2 != ltot2-m-1) || (l2 == 0) || (l1 > target_l1-m))
		clean(l1-2,l2,m+1);
              if(m == 0)
		if((l1-2 < target_l1) || (l2 < target_l3))
		  clean(l1-2,l2,0);
                else
		  contract(l1-2,l2);
	    }
	    
            if(l1 == ltot1)
	    {
	      if(l2 > 1)
	      {
	        clean(l1,l2-2,m+1);
	        clean(l1-1,l2-2,m+1);
		if(m == 0) 
		{
		  if(l2-2 < target_l3)
		    clean(l1,l2-2,0);
                  else
		    contract(l1,l2-2);
		  
		  if(l2-2 < target_l3)
		    clean(l1-1,l2-2,0);
                  else
		    contract(l1-1,l2-2);
		}
	      }
	      if(l2 == ltot2-m)
	      {
		clean(l1,l2-1,m+1);
                clean(l1-1,l2-1,m+1);
		if(m == 0) 
		{
		  if(l2-1 < target_l3)
		    clean(l1,l2-1,0);
                  else
		    contract(l1,l2-1);
		  if(l2-1 < target_l3)
		    clean(l1-1,l2-1,0);
                  else
		    contract(l1-1,l2-1);
		}
	      }
	    }
	    if((m == 0) && (l2 > ltot2-1) && ( l1 >= target_l1))
		contract(l1,l2);
	  }
	}
      }
    }
    else
    {
      if(m>=ltot1)
      {
	for(l2 = 0; l2 <= ltot1+ltot2-m; l2++)
        {
	  reduce(0,l2,m,1);
	  if((l2 > 1) && (m < ltot1+ltot2))
          {
	    clean(0,l2-2,m+1);
	    if(m == 0)
	      if((l2-2 < target_l3) || ( 0 < target_l1))
		clean(0,l2-2,0);
	      else
		contract(0,l2-2);
	  }
	}
	if(m < ltot1+ltot2) 
          clean(0,ltot1+ltot2-m-1,m+1);
	if(m == 0)
	{
           if(ltot2 > 0)
	     if((ltot2-1 < target_l3) || ( 0 < target_l1))
	       clean(0,ltot2-1,0);
	     else
	       contract(0,ltot2-1);
	   contract(0,ltot2);
	}
      }
      else
      {
	for(lsum = 0; lsum <= ltot1+ltot2-m; lsum++)
	{
	  for(l2 = 0; l2 <= lsum; l2++)
	  {
	    l1 = lsum -l2;
	    if(l2 >ltot2) continue;

	    if(l1 > ltot1-m)
	      continue;

	    if((l1 > ltot1-1-m) && (l2 < target_l3-m))
	      continue;

            if(l1 == 0)
              reduce(0,l2,m,1);
	    else
	      reduce(l1,l2,m,0);

            if((l2 > 1) && (l1 < ltot1-m))
            {
              if((l1 != ltot1-m-1) || (l1 == 0) || (l2 > target_l3-m))
		clean(l1,l2-2,m+1);
              if(m==0)
                if((l1 <target_l1) || (l2-2 <target_l3))
		  clean(l1,l2-2,m);
                else
	          contract(l1,l2-2);
	    }

            if(l2 == ltot2)
	    {
	      if(l1 > 1)
	      {
	        clean(l1-2,l2,m+1);
	        clean(l1-2,l2-1,m+1);
		if(m == 0)
 		{
                  if(l1-2 < target_l1)
		    clean(l1-2,l2,0);
		  else
	            contract(l1-2,l2);
                  if(l1-2 < target_l1)
		    clean(l1-2,l2-1,0);
		  else
	            contract(l1-2,l2-1);
		}
	      }
	      if(l1 == ltot1-m)
	      {
		clean(l1-1,l2,m+1);
		clean(l1-1,l2-1,m+1);
		if(m == 0)
		{
                  if(l1-1 < target_l1)
		    clean(l1-1,l2,0);
		  else
	            contract(l1-1,l2);
                  if(l1-1 < target_l1)
		    clean(l1-1,l2-1,0);
		  else
	            contract(l1-1,l2-1);
		}
	      }
	    }
	    if((m == 0) && (l1 > ltot1-1) && ( l2 >= target_l3))
	      contract(l1,l2);
	  }
	}
      }
    }
  }

}

void VRR::contract(const int l1, const int l2)
{
  // now do the contractions for the remaining arrays
  
  Momentum zero;

  timer -> start(t_contract);
  
  for(Momentum m1(l1); m1(); ++m1)
  {
    for(Momentum m2(l2); m2(); ++m2)
    {	
      if (debug)
	cout << " Contracting: " << m1 << " " << m2 << endl;
      int idx = compute_index(m1,m2,0,zero);
      IntegArray* array = fields[idx];
      assert(array != 0);
      if (debug)
	cout << "Source: " << endl << *array << endl;
	
      IntegArray* contracted = new IntegArray(*array,target_l1,target_l2,
						target_l3,target_l4,qset,
						alloc_contract);
      assert(contracted != 0);
      if (debug)
	cout << "Result: " << endl << *contracted << endl;
      delete array;
      fields[idx] = contracted;
    }
  }
  
  timer -> stop(t_contract);
}
/*TEX
\subsection{Evaluation of Recursion Relations}
  
As in the \name{Recursion} class, the function \name{reduce} is the
workhorse to evaluate the terms of the recursion relation. This
function generates new \name{IntegArrays} for all pairs of angular
momenta $(\va,\vc)$ with total angular momentum $l_1$ on the first and
$l_2$ on the second index.  
*/

void VRR::reduce(const int l1,const int l2,const int aux, const int position)
{
  // now apply the the recursion formula term by term for all angular momenta
    
  IntegArray* term;
  Coef_Type   ct;
  Momentum    zero(0);
     
  if((l1 == 0) && (l2 == 0))
  {
    if(debug) 
      cout << " VRR::reduce -- ss-integ: " << aux << endl;
    ss_integ(aux);
    return;
  }
  

  for(Momentum m1(l1); m1(); ++m1)
  for(Momentum m2(l2); m2(); ++m2)
  {
    
    // find the target IntegArray

    if (debug) 
      cout << "VRR::reduce -- Computing Field: " << m1 << m2 << " aux: " 
	   << aux << endl; 
    IntegArray *target    =  find(m1,m2,aux,zero);
    if (target)
    {
      if (debug) cout << "already done ! " << endl;
      return;  // already done !
    }

    CoefArray& ca = qset.obtain(RHOSUM);
    
    target = new IntegArray(alloc_prim);
    

    fields[compute_index(m1,m2,aux,zero)] = target;
    assert(target != 0);

    Direction dir = (position == 0) ? m1.decrease() : m2.decrease();

    Momentum m1base(m1);    
    Momentum m2base(m2);

    if (position == 0)
      m1base.reduce(dir);
    else
      m2base.reduce(dir);

    Momentum  m1baseminus(m1base);
    m1baseminus.reduce(dir);
    Momentum  m2baseminus(m2base);
    m2baseminus.reduce(dir);

    // term: (p_i-a_i) (a-1b|V|cd)(m) 
      
    ct       = (Coef_Type) (PA1X + dir);
    if (debug)
      cout << "   --> Term  0: " << ct << " ";
    term     =  find(m1base,m2base,aux,zero);
    assert(term != 0);
    target -> doublet(qset.obtain(ct,position),*term,position);

    // term: W_i (a-1b|V|cd)(m+1) 
      
    ct       = (Coef_Type) (WX + dir);
    if (debug)
      cout << "   --> Term 1A: " << ct << " ";
    term     =  find(m1base,m2base,aux+1,zero);
    assert(term != 0);
    //target -> quadruplet(1.0,qset.obtain(ct),*term);
    target -> add_to(1.0,qset.obtain(ct),*term);

    // term: -P_i (a-1b|V|cd)(m+1) 
      
    ct       = (Coef_Type) (PX + dir);
    if (debug)
      cout << "   --> Term 1B: " << ct << " ";
    term     =  find(m1base,m2base,aux+1,zero);
    assert(term != 0);
    target -> doublet(-1,qset.obtain(ct,position),*term,position);
   
    if (position == 0) // reduction on a
    {
      if (m1baseminus(dir) >= 0)
      {
	ct       = ZETASUM;
	if (debug)
	  cout << "   --> Term  2: " << ct << " ";
	term     =  find(m1baseminus,m2base,aux,zero);
	assert(term != 0);
	target -> doublet(m1base(dir),qset.obtain(ct,position),*term,
			  position);

	ct       = RHO1;
	if (debug)
	  cout << "   --> Term  3: " << ct << " ";
	term     =  find(m1baseminus,m2base,aux+1,zero);
	assert(term != 0);
	//target -> quadruplet(l1base(dir),qset.obtain(ct),*term);
	target -> add_to(m1base(dir),qset.obtain(ct),*term);
      }
      
      if (m2baseminus(dir) >= 0)
      {	
	ct       = RHOSUM;
	if (debug)
	  cout << "   --> Term  4: " << ct << " ";
	term     =  find(m1base,m2baseminus,aux+1,zero);
	assert(term != 0);
	// target -> quadruplet(m2base(dir),qset.obtain(ct),*term);
	target -> add_to(m2base(dir),qset.obtain(ct),*term);
      }
      
    }
    else
    {
      if (m2baseminus(dir) >= 0)
      {
	ct       = ZETASUM;
	if (debug)
	  cout << "   --> Term  2: " << ct << " ";
	term     =  find(m1base,m2baseminus,aux,zero);
	assert(term != 0);
	target -> doublet(m2base(dir),qset.obtain(ct,position),*term,position);
	
	ct       = RHO2;
	if (debug)
	  cout << "   --> Term  3: " << ct << " ";
	term     =  find(m1base,m2baseminus,aux+1,zero);
	assert(term != 0);
	// target -> quadruplet(m2base(dir),qset.obtain(ct),*term);
	target -> add_to(m2base(dir),qset.obtain(ct),*term);
      }
      
      if (m1baseminus(dir) >= 0)
      {
	ct       = RHOSUM;
	if (debug)
	  cout << "   --> Term  3: " << ct << " ";
	term     =  find(m1baseminus,m2base,aux+1,zero);
	assert(term != 0);
	// target -> quadruplet(m1base(dir),qset.obtain(ct),*term);
	target -> add_to(m1base(dir),qset.obtain(ct),*term);
      }
      
    }
    if (debug)
      cout << "Result: " << endl << *target << endl;
  }
}

/*TEX
\subsection{Output of ERI's in the VRR} 
*/

void  VRR::output(InternalBasis& b1,SymDesignator& s1,
		  InternalBasis& b2,SymDesignator& s2,
		  InternalBasis& b3,SymDesignator& s3,
		  InternalBasis& b4,SymDesignator& s4,
		  int same12,int same34,int samepair,
		  IntegFile& obuf)
{
  Momentum zero(0);
  if (!flaglist(Spherical))
  {
    for(Momentum m1(l1);m1();++m1)
    for(Momentum m2(l2);m2();++m2)
    {
      IntegArray* primary = find(m1,m2,0,zero);
      if (primary)
      {
	primary -> output(b1,s1,b2,s2,b3,s3,b4,s4,l1,m1.local_index(),0,0,
			  l2,m2.local_index(),0,0,obuf,0,
			  same12,same34,samepair);
      }
    }
  }
  else
  {
    Sphere sinteg(*this,1,alloc_contract);
    for(int m1=0; m1 < 2*l1+1;m1++)
      for(int m2=0; m2 < 2*l2+1;m2++)
      {
	IntegArray* integ = sinteg.find(m1,m2,0,zero);      
	assert(integ != 0);
	integ -> output(b1,s1,b2,s2,b3,s3,b4,s4,l1,m1,0,0,l2,m2,0,0,obuf,0,
			same12,same34,samepair);
      }
  }  
}

VRR::~VRR()
{
  clean();
}
