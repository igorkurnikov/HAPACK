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
% $Id: rec.cpp,v 1.2 2008/05/03 06:23:55 igor Exp $
%
\section{Implementation of Recursion}
*/
#include "parallel.h"
#include "io_incl.h"
// #include <stream.h>
#include <assert.h>
#include <math.h>
#include "typesip.h"
#include "function.h"
#include "integ_array.h"
#include "coef_array.h"
#include "coef_set.h"
#include "rec.h"
#include "basis.h"

Location DummyCenter;

/*TEX
\subsection{Constructor}
\subsection{Constructor}
\index{Recursion, constructor}
  
The constructor merely assigns the arguments to the local variables
and passes the data neccessary for allocating the required
storage to \name{RecursionStorage}.
 Since \name{Recursion} is an abstract class the
initiation of the \name{logic} subroutine and the subsequent evaluation of
the recursion relation is reserved for the constructors of the derived
classes. 

The allocator in \name{Recursion} is dimensioned to the size of the contracted 
arrays. 
*/

Recursion::Recursion(Coefficient_Set& cset,int momentum1,int momentum2,
		     int max_aux,int max_mu,Location cntr,Recursion* rec) : 
     coef_set(cset), RecursionStorage(momentum1,momentum2,max_aux,max_mu,
				      cset.allocator()),
     alloc_contract(cset.radials(0).nofun(),cset.radials(1).nofun(),1,1,1)
{
  auxop      = rec;
  center     = cntr;
  debug      = false;
}

/*TEX 
\subsection{Logic} 
\index{Recursion, logic}
This is the `generating' function for the recursion relations.
Eventually \name{logic} will call the recursion relation to
evaluate the integrals for all pairs of angular momenta which are
designated by the arguments. This means that at the end of a call to
\name{logic} all integrals with total angular momentum $l_1$ on the first
and total angular momentum $l_2$ on the second index, which have the
scalar auxiliary index \name{aux} and ALL angular momentum
auxiliary indices up-to and including total angular momenta
$\mu_{tot}$ will have been computed. (This must be modified for
angular momentum and spin-orbit interaction by setting $mumin = 1$,
which is by default $0$.The \name{minmu} variable is just carried 
through the recursive function calls
and will only be needed in \name{reduce}, because 
angular momentum and spin-orbit integrals do not require the $f2$
term with $O(\vmu -\one_i)$.)

Pointers to the corresponding
\name{IntegArrays}, which hold the values, have been assigned in the
\name{fields} array of the \name{RecursionStorage} class. 

As mentioned before, we must make a difference on whether using the
auxiliary index $m$ (which we now call \name{aux}) or not.
In case \name{aux} is not needed,
logic should be called with \name{aux} = 0 and \name{use_aux} = 0. 
The \name{use_aux} flag is on default 1.

For the ss shell ($l_1 = l_2 = 0$) this function starts the
$\mu-reduction$ or generates the primary integral depending on the value
of $\mu_{tot}$.
For all other values of $l_1$ and $l_2$, the function assures
that all \name{IntegArrays} representing source terms in the recursion
relation, which \name{logic} wants to evaluate, are present.  This is
done by recursively calling \name{logic} and at least \name{speclogic}
for the arguments representing the angular momenta and auxiliary 
parameters of the required terms.

Finally the \name{reduce} function is called to carry out the recursion.

Which terms are required in the recursive calls to \name{logic} is
determined
\begin{itemize}
\item[a)] by the values of the $f_i$ flags and
\item[b)] by the special terms encoded in the \name{speclogic}
function.
\end{itemize}
*/
   
void Recursion::logic(const int ltot1,const int ltot2,
		      const int aux,const int mu,int use_aux,int minmu)
{
  if(debug)
    cout << "logic for: " << ltot1 << "  " << ltot2 << " aux: " << aux <<endl;

  // do we have an (s|O|s) integral ?
    
  if (ltot1 == 0 && ltot2 == 0)
  {
    if( mu > 0 )
      mulogic(ltot1,ltot2,aux,mu);
    else
      ss_integ(aux,mu);

    return;
  }
  
  // always reduce the larger index:
  
  
  int position,l1base,l2base;
  if (ltot1 > ltot2)
  {
    l1base   = ltot1 - 1;
    l2base   = ltot2;
    position = 0;
  }  
  else 
  {
    position = 1;
    l1base   = ltot1;
    l2base   = ltot2-1;
  }
  
  if( use_aux )
  {
    if (f0)
      logic(l1base,l2base,aux,mu,1,minmu);
    if (f1)
      logic(l1base,l2base,aux+1,mu,1,minmu);
    if (f2 && mu > 0 )
      logic(l1base,l2base,aux+1,mu-11,minmu);
    
    if (l1base-1 >= 0) 
    {
      if (f3)
	logic(l1base-1,l2base,aux,mu,1,minmu);
      if (f4)
	logic(l1base-1,l2base,aux+1,mu,1,minmu);
    }
    if (l2base-1 >= 0) 
    {
      if (f3)
	logic(l1base,l2base-1,aux,mu,1,minmu);
      if (f4)
	logic(l1base,l2base-1,aux+1,mu,1,minmu);
    }
  }
  else
  {
    if (f0)
      logic(l1base,l2base,0,mu,0,minmu);
    if (f1)
      logic(l1base,l2base,0,mu,0,minmu);
    if (f2 && mu > 0 )
      logic(l1base,l2base,0,mu-1,0,minmu);
    
    if (l1base-1 >= 0) 
    {
      if (f3)
	logic(l1base-1,l2base,0,mu,0,minmu);
      if (f4)
	logic(l1base-1,l2base,0,mu,0,minmu);
    }
    if (l2base-1 >= 0) 
    {
      if (f3)
	logic(l1base,l2base-1,0,mu,0,minmu);
      if (f4)
	logic(l1base,l2base-1,0,mu,0,minmu);
    }
  }    
    
  speclogic(ltot1,ltot2,aux,mu);
    
  reduce(ltot1,ltot2,aux,mu,use_aux,minmu);
}

/*TEX 
\subsection{Reduction\label{reduce}}
\index{Recursion, reduce}

This function is the workhorse of the \name{Recursion} class, which
applies the recursion relation to generate all integrals with total angular
momentum values (ltot1,ltot2) respecitvely. The scalar auxiliary
index is \name{aux}, while ALL auxiliary momentum indices from
total momentum \name{mumin} up to \name{mumax} are computed. 
The default value of \name{mumin} is $0$.

We first determine the position of reduction, as in logic, we always
reduce the {\em larger} angular momentum. This choice MUST be consistent.

Then follows a loop over all angular momentum indices with the
appropriate angular momenta for $\va,\vb$ and $\vmu$. Within this
loop, we first determine whether the term has already been computed.
If this is the case we quit.

Otherwise we choose a \name{Direction} of reduction and compute the
angular momenta $\va'$ and $\vb'$ as defined in
section\ref{section-general-terms}. 

Next we have to decide whether the scalar auxiliary parameter
\name{aux} is used or not. Default is $use_aux = 1$.

Depending on the values of the
$f_i$ flags, we then compute the appropriate terms.
If, for example, $f_0 = 1$, we need to evaluate a term for a reduction
on the first index (remember $m = aux):  
\be
 (P_i-A_i) (\va'\,|\,O(\vmu,m)|\,\vb')
\ee
The appropriate \name{CoefArray} has the name \name{PA1X,PA1Y} or
\name{PA1Z}, defined as an enumaration to label \name{CoefArrays}.
If $d$ is
the integer value for the the direction of reduction, the appropriate
name is thus: \name{PA1X + d}. (If we were reducing on the second
index, the approrpiate name would be \name{PA2X + d}. 
Next we call \name{find} to obtain the pointer to
the \name{IntegArray} indexed by $(\va',\vb',m,\vmu)$. The function
\name{logic} has assured that this term has previously been computed.
Now we have accumulated the ingredients for the term in the recursion
relation, which is evalued in the function \name{recursion\_term},
which does the actual work.

In a similar fashion all other terms can be treated. The basic
operations to carry out the evaluation of an individual term are
encoded in the \name{IntegArray} class. 
 
If you look at the $f1$-term, you will notice that we first use
the \name{CoefArray} \name{PX + d} and in a second step subtract the
center coordinate. This of course could be combined by introducing
for example the \name{CoefArrays} \name{PCX}, \name{PCY} and
\name{PCZ}, which will make the calculation quicker.
On the other hand we then need a larger \name{CoefficientSet},
because more specific \name{CoefArrays} have to be created. Since the
runtime is mainly determinded by the calculation of the two elctron integrals,
we decided not to change this, yet.
*/

void Recursion::reduce(const int ltot1,const int ltot2,const int aux,
		       const int mumax,int use_aux, int mumin)
{
   // now apply the the recursion formula term by term for all angular momenta
  
  int i;
   
  IntegArray* term;
  Coef_Type   ct;

  int position = (ltot1 > ltot2);

  for(int m=mumin; m <= mumax;m++)
  for(Momentum mu(m);mu();++mu)
  for(Momentum l1(ltot1); l1(); ++l1)
  for(Momentum l2(ltot2); l2(); ++l2)
  {
    if (debug) 
      cout << "  *** REDUCE *** l1 : " << l1 << " l2: " << l2 
           << " aux: " << aux << " mu: " << mu << endl; 

    // find the target IntegArray

    if (debug) 
      cout << "   Computing Target Field: "<< compute_index(l1,l2,aux,mu) 
           << endl;

    IntegArray *target    =  find(l1,l2,aux,mu);
    if (target)
    {
      if (debug) cout << "    Target already done ! " << endl;
      return;  // already done !
    }

    //target = new IntegArray(coef_set[ZETASUM].size1(),
    //	 		     coef_set[ZETASUM].size2());
    target = new IntegArray(alloc_prim);
    assert(target -> size1()  == coef_set[ZETASUM].size1());
    assert(target -> size2()  == coef_set[ZETASUM].size2());
    assert(target -> size3()  == 1);
    assert(target -> size4()  == 1);
    assert(target -> blocks() == 1);

    fields[compute_index(l1,l2,aux,mu)] = target;
    assert(target != 0);

    if(debug)
    {
      cout << "   Contents of fields[" <<compute_index(l1,l2,aux,mu) << "]\n";
      for(i = 0; i < fields[compute_index(l1,l2,aux,mu)]->size();i++)
      {
	cout << "    Element " << i << " : " 
             << fields[compute_index(l1,l2,aux,mu)]->base[i] << endl;
      }
    }

    Direction dir         = (position) ? l1.decrease() : l2.decrease();

    Momentum l1base(l1);    
    Momentum l2base(l2);

    if (position)
      l1base.reduce(dir);
    else
      l2base.reduce(dir);

    Momentum  l1baseminus(l1base);
    l1baseminus.reduce(dir);
    Momentum  l2baseminus(l2base);
    l2baseminus.reduce(dir);

    if( use_aux )
    {
      if (f0)
      {
	// term: (p_i-a_i) (a-1|A(mu)|b)(m) 
      
	 ct       = (position) ? (Coef_Type) (PA1X + dir)
	              : (Coef_Type) (PA2X + dir);

	term     =  find(l1base,l2base,aux,mu);
	assert(term != 0);

	target -> recursion_term(f0,coef_set[ct],*term);

	if (debug)
	{
	  cout << "   --> Term 1 :  multiply with " << ct 
               << " f0: " << f0 << endl;
	  cout << "       Coef_Array " << ct << endl;
	  //	  for(i = 0; i < coef_set[ct].size();i++)
	  //  cout << form("        Element %d : %e ", i, coef_set[ct].base[i])
	  //    << endl;
	  cout << "       Source IntegArray "<< endl
	       << "         l1: " << l1base << " l2: " << l2base
	       << " aux: " << aux << " mu: " << mu 
	       << " INDEX: " << compute_index(l1base,l2base,aux,mu)
	       << endl;
	  // for(i = 0; i < term->size();i++)
	  // cout << form("        Element %d : %e \n", i, term->base[i]);
	  cout << "    Resulting Target IntegArray " << endl;
	  //for(i = 0; i < target->size();i++)
	  // cout << form("        Element %d : %e \n", i, target->base[i]);
	}
      }

      if (f1)
      {
	// term: (p_i-c_i) (a-1|A(mu)|b)(m+1) 
      
	ct       = (Coef_Type) (PX + dir);

	term     =  find(l1base,l2base,aux+1,mu);
	assert(term != 0);
	
	target -> recursion_term(f1,coef_set[ct],*term);

	if (debug)
	{
	  cout << "   --> Term 2A :  multiply with " << ct 
               << " f1: " << f1 << endl;
	  cout << "       Coef_Array " << ct << endl;
	  //for(i = 0; i < coef_set[ct].size();i++)
	    //cout << form("        Element %d : %e ", i, coef_set[ct].base[i])
	    //  << endl;
	  cout << "       Source IntegArray "<< endl
	       << "         l1: " << l1base << " l2: " << l2base
	       << " aux: " << aux+1 << " mu: " << mu 
	       << " INDEX: " << compute_index(l1base,l2base,aux+1,mu)
	       << endl;
	  //for(i = 0; i < term->size();i++)
	  //  cout << form("        Element %d : %e \n", i, term->base[i]);
	  cout << "    Resulting Target IntegArray " << endl;
	  //for(i = 0; i < target->size();i++)
	  //  cout << form("        Element %d : %e \n", i, target->base[i]);
	}

	double cfac = -center(dir);
	
	target -> recursion_term(f1*cfac,*term);

	if (debug)
	{
	  // cout << "   --> Term 2B :  multiply with -C(" << dir << ") : " 
	  // << form("%e",cfac) << " f1: " << f1 << endl;
	  cout << "       Source IntegArray "<< endl
	       << "         l1: " << l1base << " l2: " << l2base
	       << " aux: " << aux+1 << " mu: " << mu 
	       << " INDEX: " << compute_index(l1base,l2base,aux+1,mu)
	       << endl;
	  // for(i = 0; i < term->size();i++)
	  // cout << form("        Element %d : %e \n", i, term->base[i]);
	  cout << "    Resulting Target IntegArray " << endl;
	  // for(i = 0; i < target->size();i++)
	  // cout << form("        Element %d : %e \n", i, target->base[i]);
	}
      }

      if (f2 && mu(dir) > 0 )
      {
	// term:  (a-1|A(mu-1)|b)(m+1) 
	Momentum muminus(mu);
	muminus.reduce(dir);

        term     =  find(l1base,l2base,aux+1,muminus);
        assert(term != 0);

	target -> recursion_term(f2*mu(dir),*term);

	if (debug)
	{
	  cout << "   --> Term 3 :  multiply with mu(" << dir << "): "  
               << " f2: " << f2 << endl;
	  cout << "       Source IntegArray "<< endl
	       << "         l1: " << l1base << " l2: " << l2base
	       << " aux: " << aux+1 << " mu: " << muminus 
	       << " INDEX: " << compute_index(l1base,l2base,aux+1,muminus)
	       << endl;
	  //for(i = 0; i < term->size();i++)
	   // cout << form("        Element %d : %e \n", i, term->base[i]);
	  cout << "    Resulting Target IntegArray " << endl;
	//  for(i = 0; i < target->size();i++)
//	    cout << form("        Element %d : %e \n", i, target->base[i]);
	}
      }

      if(l1base(dir) > 0)
      {
	if (f3)
	{
	  // term: zetasum * N_{dir}(a-1) (a-2|A(mu)|b)(m) 

	    ct       = (Coef_Type) (ZETASUM);
	  
	  term     =  find(l1baseminus,l2base,aux,mu);
	  assert(term != 0);
	  
	  target -> recursion_term(f3*l1base(dir),coef_set[ct],*term);

	  if (debug)
	  {
	    cout << "   --> Term 4/1 :  multiply with " << ct 
                 << " N(" << dir << "): " << l1base(dir) 
                 << " f3: " << f3 << endl;
	    cout << "       Coef_Array " << ct << endl;
	    //	    for(i = 0; i < coef_set[ct].size();i++)
	    //  cout << form("        Element %d : %e ", i, coef_set[ct].base[i])
	    // << endl;
	    cout << "       Source IntegArray "<< endl
	         << "         l1: " << l1baseminus << " l2: " << l2base
		 << " aux: " << aux << " mu: " << mu 
	         << " INDEX: " << compute_index(l1baseminus,l2base,aux,mu)
	         << endl;
	//    for(i = 0; i < term->size();i++)
	 //     cout << form("        Element %d : %e \n", i, term->base[i]);
	    cout << "    Resulting Target IntegArray " << endl;
//	    for(i = 0; i < target->size();i++)
//	      cout << form("        Element %d : %e \n", i, target->base[i]);
	  }
        }
	if (f4)
	{
	  // term: zetasum * N_{dir}(a-1) (a-2|A(mu)|b)(m) 
	  
	  ct       = (Coef_Type) (ZETASUM);

	  term     =  find(l1baseminus,l2base,aux+1,mu);
	  assert(term != 0);

	  target -> recursion_term(f4*l1base(dir),coef_set[ct],*term);

	  if (debug)
	  {
	    cout << "   --> Term 5/1 :  multiply with " << ct 
                 << " N(" << dir << "): " << l1base(dir) 
                 << " f4: " << f4 << endl;
	    cout << "       Coef_Array " << ct << endl;
//	    for(i = 0; i < coef_set[ct].size();i++)
//	      cout << form("        Element %d : %e ", i, coef_set[ct].base[i])
	//	<< endl;
	    cout << "       Source IntegArray "<< endl
	         << "         l1: " << l1baseminus << " l2: " << l2base
		 << " aux: " << aux+1 << " mu: " << mu 
	         << " INDEX: " << compute_index(l1baseminus,l2base,aux+1,mu)
	         << endl;
	  //  for(i = 0; i < term->size();i++)
	   //   cout << form("        Element %d : %e \n", i, term->base[i]);
	    cout << "    Resulting Target IntegArray " << endl;
	   // for(i = 0; i < target->size();i++)
	     // cout << form("        Element %d : %e \n", i, target->base[i]);
	  }
	}
      }

      if(l2base(dir) > 0)
      {
	if(f3)
	{      
	  // term: zeta1 * N_{dir}(a-1) (a-2|A(mu)|b)(m)       
	  ct       = (Coef_Type) (ZETASUM);

	  term     =  find(l1base,l2baseminus,aux,mu);
	  assert(term != 0);

	  target -> recursion_term(f3*l2base(dir),coef_set[ct],*term);

	  if (debug)
	  {
	    cout << "   --> Term 4/2 :  multiply with " << ct 
                 << " N(" << dir << "): " << l2base(dir) 
                 << " f3: " << f3 << endl;
	    cout << "       Coef_Array " << ct << endl;
	    for(i = 0; i < coef_set[ct].size();i++)
	      (ostream&) cout <<  i <<  coef_set[ct].base[i] << endl;
	    cout << "       Source IntegArray "<< endl
	         << "         l1: " << l1base << " l2: " << l2baseminus
		 << " aux: " << aux << " mu: " << mu 
	         << " INDEX: " << compute_index(l1base,l2baseminus,aux,mu)
	         << endl;
	    for(i = 0; i < term->size();i++)
	      (ostream&) cout <<  i <<  term->base[i] << endl;
	    cout << "    Resulting Target IntegArray " << endl;
	    for(i = 0; i < target->size();i++)
	       (ostream&) cout <<  i <<  target->base[i] << endl;
	  }
	}
	if(f4)
	{      
	  // term: zeta1 * N_{dir}(a-1) (a-2|A(mu)|b)(m+1)       
	  ct       = (Coef_Type) (ZETASUM);

	  term     =  find(l1base,l2baseminus,aux+1,mu);
	  assert(term != 0);

	  target -> recursion_term(f4*l2base(dir),coef_set[ct],*term);

	  if (debug)
	  {
	    cout << "   --> Term 5/2 :  multiply with " << ct 
                 << " N(" << dir << "): " << l2base(dir) 
                 << " f4: " << f4 << endl;
	    cout << "       Coef_Array " << ct << endl;
	    for(i = 0; i < coef_set[ct].size();i++)
	       (ostream&) cout <<  i <<  coef_set[ct].base[i] << endl;
	    cout << "       Source IntegArray "<< endl
	         << "         l1: " << l1base << " l2: " << l2baseminus
		 << " aux: " << aux+1 << " mu: " << mu 
	         << " INDEX: " << compute_index(l1base,l2baseminus,aux+1,mu)
	         << endl;
	    for(i = 0; i < term->size();i++)
	       (ostream&) cout << i <<  term->base[i] << endl;
	    cout << "    Resulting Target IntegArray " << endl;
	    for(i = 0; i < target->size();i++)
	       (ostream&) cout <<  i <<  target->base[i] << endl;
	  }
	}
      }
    }
    else
    {
      if (f0)
      {
	// term: (p_i-a_i) (a-1|A(mu)|b)(0) 
      
	ct       = (position) ? (Coef_Type) (PA1X + dir)
                              : (Coef_Type) (PA2X + dir);

	term     =  find(l1base,l2base,0,mu);
	assert(term != 0);

	target -> recursion_term(f0,coef_set[ct],*term);

	if (debug)
	{
	  cout << "   --> Term 1 :  multiply with " << ct 
               << " f0: " << f0 << endl;
	  cout << "       Coef_Array " << ct << endl;
	  for(i = 0; i < coef_set[ct].size();i++)
	    (ostream&) cout <<  i <<  coef_set[ct].base[i] << endl;
	  cout << "       Source IntegArray "<< endl
	       << "         l1: " << l1base << " l2: " << l2base
	       << " aux: " << 0 << " mu: " << mu 
	       << " INDEX: " << compute_index(l1base,l2base,0,mu)
	       << endl;
	  for(i = 0; i < term->size();i++)
	    (ostream&)	    cout << i <<  term->base[i] << endl;
	  cout << "    Resulting Target IntegArray " << endl;
	  for(i = 0; i < target->size();i++)
	    (ostream&) cout << i <<  target->base[i] << endl;
	}
      }

      if (f1)
      {
	// term: (p_i-c_i) (a-1|A(mu)|b)(0) 
      
	ct       = (Coef_Type) (PX + dir);

	term     =  find(l1base,l2base,0,mu);
	assert(term != 0);

	target -> recursion_term(f1,coef_set[ct],*term);

	if (debug)
	{
	  cout << "   --> Term 2A :  multiply with " << ct 
               << " f1: " << f1 << endl;
	  cout << "       Coef_Array " << ct << endl;
	  for(i = 0; i < coef_set[ct].size();i++)
	     (ostream&) cout <<  i << coef_set[ct].base[i] << endl;
	  cout << "       Source IntegArray "<< endl
	       << "         l1: " << l1base << " l2: " << l2base
	       << " aux: " << 0 << " mu: " << mu 
	       << " INDEX: " << compute_index(l1base,l2base,0,mu)
	       << endl;
	  for(i = 0; i < term->size();i++)
	     (ostream&) cout << i <<  term->base[i] << endl;
	  cout << "    Resulting Target IntegArray " << endl;
	  for(i = 0; i < target->size();i++)
	     (ostream&) cout <<  i <<  target->base[i] << endl;
	}

	double cfac = -center(dir);

	target -> recursion_term(f1*cfac,*term);

	if (debug)
	{
	  cout << "   --> Term 2B :  multiply with -C(" << dir << ") : " 
               << cfac << " f1: " << f1 << endl;
	  cout << "       Source IntegArray "<< endl
	       << "         l1: " << l1base << " l2: " << l2base
	       << " aux: " << 0 << " mu: " << mu 
	       << " INDEX: " << compute_index(l1base,l2base,0,mu)
	       << endl;
	  for(i = 0; i < term->size();i++)
	    (ostream&)	    cout << i << term->base[i] << endl;
	  cout << "    Resulting Target IntegArray " << endl;
	  for(i = 0; i < target->size();i++)
	    (ostream&) cout << i << target->base[i] << endl;
	}
      }

      if (f2 && mu(dir) > 0 )
      {
	// term:  zetasum * N_{dir}(mu) (a-1|A(mu-1)|b)(0) 

        ct       = (Coef_Type) (ZETASUM);

	Momentum muminus(mu);
	muminus.reduce(dir);

        term     =  find(l1base,l2base,0,muminus);
	assert(term != 0);

	target -> recursion_term(f2*mu(dir),coef_set[ct],*term);

	if (debug)
	{
	  cout << "   --> Term 3 :  multiply with mu(" << dir << "): "  
               << " f2: " << f2 << endl;
	  cout << "       Source IntegArray "<< endl
	       << "         l1: " << l1base << " l2: " << l2base
	       << " aux: " << 0 << " mu: " << muminus 
	       << " INDEX: " << compute_index(l1base,l2base,0,muminus)
	       << endl;
	  for(i = 0; i < term->size();i++)
	     (ostream&) cout <<  i <<  term->base[i] << endl;
	  cout << "    Resulting Target IntegArray " << endl;
	  for(i = 0; i < target->size();i++)
	     (ostream&) cout << i <<  target->base[i] << endl;
	}
      }

      if(l1base(dir) > 0)
      {
	if (f3)
	{
	  // term: zetasum * N_{dir}(a-1) (a-2|A(mu)|b)(0) 

          ct       = (Coef_Type) (ZETASUM);

	  term     =  find(l1baseminus,l2base,0,mu);
	  assert(term != 0);

	  target -> recursion_term(f3*l1base(dir),coef_set[ct],*term);

	  if (debug)
	  {
	    cout << "   --> Term 4/1 :  multiply with " << ct 
                 << " N(" << dir << "): " << l1base(dir) 
                 << " f3: " << f3 << endl;
	    cout << "       Coef_Array " << ct << endl;
	    for(i = 0; i < coef_set[ct].size();i++)
	       (ostream&) cout <<  i <<  coef_set[ct].base[i] << endl;
	    cout << "       Source IntegArray "<< endl
	         << "         l1: " << l1baseminus << " l2: " << l2base
		 << " aux: " << 0 << " mu: " << mu 
	         << " INDEX: " << compute_index(l1baseminus,l2base,0,mu)
	         << endl;
	    for(i = 0; i < term->size();i++)
	       (ostream&) cout <<  i <<  term->base[i] << endl;
	    cout << "    Resulting Target IntegArray " << endl;
	    for(i = 0; i < target->size();i++)
	      (ostream&)	      cout <<  i <<  target->base[i] << endl;
	  }
	}
	if (f4)
	{
	  // term: zetasum * N_{dir}(a-1) (a-2|A(mu)|b)(0) 

	  ct       = (Coef_Type) (ZETASUM);

	  term     =  find(l1baseminus,l2base,0,mu);
	  assert(term != 0);

	  target -> recursion_term(f4*l1base(dir),coef_set[ct],*term);

	  if (debug)
	  {
	    cout << "   --> Term 5/1 :  multiply with " << ct 
                 << " N(" << dir << "): " << l1base(dir) 
                 << " f4: " << f4 << endl;
	    cout << "       Coef_Array " << ct << endl;
	    for(i = 0; i < coef_set[ct].size();i++)
	      (ostream&)	      cout <<  i << coef_set[ct].base[i] << endl;
	    cout << "       Source IntegArray "<< endl
	         << "         l1: " << l1baseminus << " l2: " << l2base
		 << " aux: " << 0 << " mu: " << mu 
	         << " INDEX: " << compute_index(l1baseminus,l2base,0,mu)
	         << endl;
	    for(i = 0; i < term->size();i++)
	      (ostream&) cout <<  i <<  term->base[i];
	    cout << "    Resulting Target IntegArray " << endl;
	    for(i = 0; i < target->size();i++)
	      (ostream&) cout <<  i << target->base[i] << endl;
	  }
	}
      }

      if(l2base(dir) > 0)
      {
	if(f3)
	{      
	  // term: zeta1 * N_{dir}(a-1) (a-2|A(mu)|b)(0) 
	  ct       = (Coef_Type) (ZETASUM);

	  term     =  find(l1base,l2baseminus,0,mu);
	  assert(term != 0);

	  target -> recursion_term(f3*l2base(dir),coef_set[ct],*term);

	  if (debug)
	  {
	    cout << "   --> Term 4/2 :  multiply with " << ct 
                 << " N(" << dir << "): " << l2base(dir) 
                 << " f3: " << f3 << endl;
	    cout << "       Coef_Array " << ct << endl;
	    for(i = 0; i < coef_set[ct].size();i++)
	      (ostream&) cout <<  i << coef_set[ct].base[i] << endl;
	    cout << "       Source IntegArray "<< endl
	         << "         l1: " << l1base << " l2: " << l2baseminus
		 << " aux: " << 0 << " mu: " << mu 
	         << " INDEX: " << compute_index(l1base,l2baseminus,0,mu)
	         << endl;
	    for(i = 0; i < term->size();i++)
	      (ostream&) cout <<  i << term->base[i] << endl;
	    cout << "    Resulting Target IntegArray " << endl;
	    for(i = 0; i < target->size();i++)
	      (ostream&) cout << i << target->base[i] << endl;
	  }
	}
	if(f4)
	{      
	  // term: zeta1 * N_{dir}(a-1) (a-2|A(mu)|b)(0)       
	  ct       = (Coef_Type) (ZETASUM);

	  term     =  find(l1base,l2baseminus,0,mu);
	  assert(term != 0);

	  target -> recursion_term(f4*l2base(dir),coef_set[ct],*term);

	  if (debug)
	  {
	    cout << "   --> Term 5/2 :  multiply with " << ct 
                 << " N(" << dir << "): " << l2base(dir) 
                 << " f4: " << f4 << endl;
	    cout << "       Coef_Array " << ct << endl;
	    for(i = 0; i < coef_set[ct].size();i++)
	      (ostream&) cout <<  i <<  coef_set[ct].base[i] << endl;
	    cout << "       Source IntegArray "<< endl
	         << "         l1: " << l1base << " l2: " << l2baseminus
		 << " aux: " << 0 << " mu: " << mu 
	         << " INDEX: " << compute_index(l1base,l2baseminus,0,mu)
	         << endl;
	    for(i = 0; i < term->size();i++)
	       (ostream&) cout <<  i <<  term->base[i] << endl;
	    cout << "    Resulting Target IntegArray " << endl;
	    for(i = 0; i < target->size();i++)
	      (ostream&) cout <<  i << target->base[i] << endl;
	  }
	}
      }
    }
    if(debug)
    {
      cout << " *** MEAN RESULT REDUCE *** "<< endl
           << "   Contents of fields[" <<compute_index(l1,l2,aux,mu) << "]\n";
      for(i = 0; i < fields[compute_index(l1,l2,aux,mu)]->size();i++)
      {
	cout << "    Element " << i << " : " 
	     << fields[compute_index(l1,l2,aux,mu)]->base[i] << endl;
      }
      cout << " *** CALL SPECIAL *** " << endl;
    }    
    special(target,position,dir,l1,l2,mu,aux);
    if(debug)
    {
      cout << " *** FINAL RESULT REDUCE *** "<< endl
           << "   Contents of fields[" <<compute_index(l1,l2,aux,mu) << "]\n";
      for(i = 0; i < fields[compute_index(l1,l2,aux,mu)]->size();i++)
      {
	cout << "    Element " << i << " : " 
	     << fields[compute_index(l1,l2,aux,mu)]->base[i] << endl;
      }
      cout << " *** END REDUCE *** " << endl;
    }
  }
}
