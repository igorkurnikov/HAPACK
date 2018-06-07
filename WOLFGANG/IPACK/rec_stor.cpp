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
% $Id: rec_stor.cpp,v 1.2 2008/05/03 07:00:37 igor Exp $
%
\subsection{Constructor}

We need to store the maximal values for the momenta 
*/
#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include "typesip.h"
#include "storage.h"
#include "rec_stor.h"
#include "integ_array.h"

RecursionStorage::RecursionStorage(int mm1,int mm2,int mx_aux,int mx_mu,
				   int sz1,int sz2,int sz3,int sz4,int nb) 
     : alloc_prim(sz1,sz2,sz3,sz4,nb)
{
  debug      = false;
  l1         = mm1;
  l2         = mm2;
  max_aux    = mx_aux;
  max_mu     = mx_mu;
  
  // compute the number of elements:

  Momentum m1(l1);
  
  no_l1 = m1.min(l1+1);
  no_l2 = m1.min(l2+1);
  no_mu = m1.min(max_mu+1);

  if (debug)
  {
    cout << "RecursionStorage :: Constructor \n";
    cout << "L1 " << l1 << " no_l1 " << no_l1; 
    cout << " L2 " << l2 << " no_l2 " << no_l2<< endl;
    cout << "max_aux " << max_aux ;
    cout << " max_mu " << max_mu << " no_mu " << no_mu << endl;
  }
}


RecursionStorage::RecursionStorage(int mm1,int mm2,int mx_aux,int mx_mu,
				   Allocator& aref) : alloc_prim(aref)
{
  debug      = false;
  l1         = mm1;
  l2         = mm2;
  max_aux    = mx_aux;
  max_mu     = mx_mu;
  
  // compute the number of elements:

  Momentum m1(l1);
  
  no_l1 = m1.min(l1+1);
  no_l2 = m1.min(l2+1);
  no_mu = m1.min(max_mu+1);

  if (debug)
  {
    cout << "RecursionStorage :: Constructor \n";
    cout << "L1 " << l1 << " no_l1 " << no_l1; 
    cout << " L2 " << l2 << " no_l2 " << no_l2<< endl;
    cout << "max_aux " << max_aux ;
    cout << " max_mu " << max_mu << " no_mu " << no_mu << endl;
  }

}


/*TEX
\subsection{InitFields}
  
This routine allocates enough pointers in  the fields array and 
initializes all to zero. It assumes that \name{max\_mu} and \name{max\_aux}
have been set.
*/

void RecursionStorage::init_fields()
{  
  fields.reset(no_l1*no_l2*(max_aux+1)*no_mu);
 
  if (debug)
  {
    cout << "init_fields :: no_l1 " << no_l1 <<  " no_l2 " 
      << no_l2 << " max_aux " << max_aux << " no_mu " << no_mu
      << " total " << no_l1*no_l2*(max_aux+1)*no_mu << endl;
/*
    cout << " ADDRESSES (fields.size() = " << fields.size() << ")\n ";
    for(int i = 0; i < fields.size(); i++)
    { 
      cout << " i = " << i << " address ";
      cout << fields.get_base() + i* sizeof(IntegArray*) << endl;
    }
*/
  }
  
  fields.set(0);
}

/*TEX 
\subsection{Printing} 

For an (ss) type integral this function ss-integ to generate the
values.  For all higher momenta this routine calls itself recursively
according to the setting of the f-flags for all terms required for the
computation of the angular momentum shell (l1,l2,aux).  
*/

void RecursionStorage::print(const int l1,const int l2)
{
  cout << "Printing Integral Arrays for: (l1,l2) = (" << l1 << ","
    << l2 << ")" << endl;
  for(int aux=0; aux <= max_aux; aux++)
  for(int m=0;  m <= max_mu; m++)
  for(Momentum mu(m) ;mu();++mu)
  for(Momentum m1(l1);m1();++m1)
  for(Momentum m2(l2);m2();++m2)
  {
    IntegArray* a = find(m1,m2,aux,mu);
    if (a)
      cout << *a << endl;
  }
}

/*TEX
\subsection{Cleaning Up}

The \name{clean(l1,l2)} routine deletes all \name{IntegArrays} which have 
total angular momentum $l1$ and $l2$ on their respective indices.
*/

void RecursionStorage::clean()
{
  int mx = fields.size();
  for(int i=0;i < mx;i++)
    if(fields[i])
    {
      delete fields[i];
      fields[i] = 0;
    }
  fields.reset(0);
}

void RecursionStorage::clean(int ll1,int ll2)
{
  for(Momentum m1(ll1); m1(); ++m1)
  for(Momentum m2(ll2); m2(); ++m2)
  for(int aux = 0; aux <= max_aux; aux++)
  for(int mut = 0; mut <= max_mu; mut++)
  for(Momentum mu(mut); mu(); ++mu)
  {
    int idx = compute_index(m1,m2,aux,mu);
    if(fields[idx])
    {
      delete fields[idx];
      fields[idx] = 0;
    }
  }
}

void RecursionStorage::clean(int ll1,int ll2, int aux)
{
  for(Momentum m1(ll1); m1(); ++m1)
  for(Momentum m2(ll2); m2(); ++m2)
  for(int mut = 0; mut <= max_mu; mut++)
  for(Momentum mu(mut); mu(); ++mu)
  {
    if(debug)
      cout << "Cleaning m1: " << m1 << " m2: " << m2 << " aux " << aux 
           << " mu " << mu << endl;
    
    int idx = compute_index(m1,m2,aux,mu);
    if(fields[idx])
    {
      delete fields[idx];
      fields[idx] = 0;
    }
  }
}


void RecursionStorage::clean_aux()
{
  for(int ll1=0; ll1 < l1; ll1++)
  for(int ll2=0; ll2 < l2; ll2++)
    clean_aux(ll1,ll2);
}

void RecursionStorage::clean_aux(int ll1,int ll2)
{  
  for(Momentum m1(ll1); m1(); ++m1)
  for(Momentum m2(ll2); m2(); ++m2)
  for(int aux = 1; aux <= max_aux; aux++)
  for(int mut = 0; mut <= max_mu; mut++)
  for(Momentum mu(mut); mu(); ++mu)
  {
    int idx = compute_index(m1,m2,aux,mu);
    if(fields[idx])
    {
      delete fields[idx];
      fields[idx] = 0;
    }
  }
}









