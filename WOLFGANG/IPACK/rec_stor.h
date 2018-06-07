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
% $Id: rec_stor.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Storage Module for the Recursion Class \label{section:recstor}}

The following simple class maintains the storage of the pointes to the
\name{IntegArrays}, which are required in all recursion relations.

The routine \name{find} allows to obtain the pointer to a specific
\name{IntegArray} in the storage class, it will fail for illegal
arguments.

*/
#ifndef REC_STOR_H
#define REC_STOR_H

#include <tarray.h>
#include <tsstack.h>

class IntegArray;
class InternalBasis;
class SymDesignator;
class IntegMatrix;
class Allocator;

class RecursionStorage
{
  friend class Sphere;
  friend class HalfSphere;
  
protected:

  int        debug;
  int        l1,l2;

  int        no_l1;   // no of ang. momenta of total momentum less or equal l1 
  int        no_l2;   // no of ang. momenta of total momentum less or equal l2 
  int        max_mu;  // maximal total mu 
  int        no_mu;   // no of ang. momenta of total momentum less
		      // than mumax.
  int        max_aux; // maximal auxiliary index 

  Allocator            alloc_prim;
  ARRAY <IntegArray*>  fields;    // the list of pointers to the
				  // Integarrays in the recursion relations


  // the following function returns the linear index of a specific IntegArray


int compute_index(const Momentum& m1, const Momentum& m2,
		  const int aux, const Momentum& mu) const
{
  int idx = ((m1.index()*no_l2+m2.index())*(max_aux+1)+aux)*no_mu + mu.index();
/*
  if (debug)
  { 
    cout << "COMPUTE INDEX: " << m1 << m2 << " Aux: " << aux << " " << mu 
      << " Index: " << idx << endl;
    cout.flush();
  }
*/
  
  assert(aux     <=  max_aux);
  assert(mu.l()  <=  max_mu);
  assert(idx < fields.size());
  return idx;
}

public:

void     debug_on() { debug = true; }
void     debug_off() { debug = false; }

void     clean();                    // deletes all IntegArrays 
void     clean(int l1,int l2);       // deletes all IntegArrays of the
				     // specificed angular momenta for
				     // all values of the auxiliary indices.
void     clean(int l1,int l2,        // deletes all IntegArrays of the
               int aux);       	     // specificed angular momenta and
				     // auxiliary index aux.
void     clean_aux();                // deletes all IntegArrays with
				     // nonzero auxiliary index
void     clean_aux(int l1,int l2);   // the same for specific angular indices.

         RecursionStorage(int m1,int m2,int mx_aux,int mx_mu,
			  int sz11,int sz12,int sz13,int sz14,int nblock1);
         RecursionStorage(int m1,int m2,int mx_aux,int mx_mu,
			  Allocator& alloc_ref);
virtual ~RecursionStorage() { clean(); }  
 
            // return a pointer to the specified IntegArray

IntegArray* find(const Momentum& l1, const Momentum& l2,
 		 const int aux, const Momentum& mu)
                { 
		  //                  if(debug)
		  //                    cout << "FIND : ";
                  return fields[compute_index(l1,l2,aux,mu)]; 
                }   

void        insert(const Momentum& l1, const Momentum& l2,
		   const int aux, const Momentum& mu,IntegArray* ni)
                { 
		  //                  if(debug)
		  //                    cout << "INSERT : ";
                  fields[compute_index(l1,l2,aux,mu)] = ni; 
                }   

void     init_fields();              // allocate sufficient storage
				     // for the fields array and set
				     // all pointers to NULL. 

void     print (const int l1,const int l2);  
                                     // print all IntegArrays of the specified
			             // angular momenta.  
void     output(InternalBasis& b1,int ll1,SymDesignator& s1,
		InternalBasis& b2,int ll2,SymDesignator& s2,
		ARRAY<IntegMat>& m,Momentum& mu,int same,int m_sym = 1);
Allocator& allocator() { return alloc_prim; }    
};

#endif
