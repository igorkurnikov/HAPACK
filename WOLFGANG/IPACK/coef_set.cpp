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
% $Id: coef_set.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\subsection{Constructor for \name{Coefficient\_Set}}
The constructor for the \name{Coefficient\_Set} calls the initializers of the 
required \name{CoefArrays} for the one electron recursion relations. 
The initializers for 
the individual \name{CoefArrays} are used to fill the fields. 
*/
#include "io_incl.h"
//#include <stream.h>
#include <errorip.h>
#include <math.h>
#include "typesip.h"
#include "function.h"
#include "coef_array.h"
#include "coef_set.h"

Coefficient_Set::Coefficient_Set(RadialFunctionArray& f1,int l1,SymOp& s1,
				 RadialFunctionArray& f2,int l2,SymOp& s2,
				 int sflag,Integral_Type I_Type)  
     : rf1(f1), rf2(f2), alloc(f1.size(),f2.size(),1,1,1)
{
  int debug = false;
  
  if(debug)
  {
    cout << "  Creating Coefficient Set for ";
    cout << "   RadialFunctionArray f1: " << f1 << " l1: " << l1
         << " s1: " << s1 << endl;
    cout << "   RadialFunctionArray f2: " << f2 << " l2: " << l2
         << " s2: " << s2 << endl;
    cout << "   same = " << sflag << " Integral_Type : " << I_Type << endl;
  }

  sym = sflag;

  int i,j;
  reset(ZETAab_AxBZ + 1);  
  set(0);
  
  for(i=0; i<=OVERLAP;i++)
    base[i] = new CoefArray(f1,l1,s1,f2,l2,s2,(Coef_Type) i,sflag,alloc);
  
  if(debug)
  {
    for(i=0; i<=OVERLAP;i++)
    {
      cout << " *** " << (Coef_Type) i << " *** \n";
      for(j=0; j < f1.size()*f2.size();j++)
        (ostream&) cout <<  j << base[i]->base[j] << endl;
    }
  }
  

  switch(I_Type)
  {
  case KIN:
    for(i=ZETA1; i<=KINETIC;i++)
    base[i] = new CoefArray(f1,l1,s1,f2,l2,s2,(Coef_Type) i,sflag,alloc);

    if(debug)
    {
      for(i=ZETA1; i<=KINETIC;i++)
      {
	cout << " *** " << (Coef_Type) i << " *** \n";
	for(j=0; j < f1.size()*f2.size();j++)
           (ostream&) cout <<  j << base[i]->base[j] << endl;
      }
    }
    break;
  case NUC:
  case MOM:
  case TWOEL:
    for(i=PX; i<NUCLEAR;i++)
      base[i] = new CoefArray(f1,l1,s1,f2,l2,s2,(Coef_Type) i,sflag,alloc);
    if(debug)
    {
      for(i=PX; i<NUCLEAR;i++)
      {
	cout << " *** " << (Coef_Type) i << " *** \n";
	for(j=0; j < f1.size()*f2.size();j++)
	  (ostream&) cout << j << base[i]->base[j] << endl;
      }
    }
    break;
  case E_FLD:
    for(i=PX; i<NUCLEAR;i++)
      base[i] = new CoefArray(f1,l1,s1,f2,l2,s2,(Coef_Type) i,sflag,alloc);
    for(i=ZETASUM_PX; i<=ZETASUM_PZ; i++)
      base[i] = new CoefArray(f1,l1,s1,f2,l2,s2,(Coef_Type) i,sflag,alloc);
    if(debug)
    {
      for(i=PX; i<NUCLEAR;i++)
      {
	cout << " *** " << (Coef_Type) i << " *** \n";
	for(j=0; j < f1.size()*f2.size();j++)
           (ostream&) cout <<  j <<  base[i]->base[j] << endl;
	  cout << " Element " << j << " : " << base[i]->base[j] << endl;
      }
      for(i=ZETASUM_PX; i<=ZETASUM_PZ;i++)
      {
	cout << " *** " << (Coef_Type) i << " *** \n";
	for(j=0; j < f1.size()*f2.size();j++)
	  (ostream&)  cout << j << base[i]->base[j] << endl;
      }
    }
    break;
  case ANG:
    base[ZETA1] = new CoefArray(f1,l1,s1,f2,l2,s2,ZETA1,sflag,alloc);    
    base[ZETA2] = new CoefArray(f1,l1,s1,f2,l2,s2,ZETA2,sflag,alloc);    
    for(i=ZETA1_BX; i<=XI_AxBZ; i++)
      base[i] = new CoefArray(f1,l1,s1,f2,l2,s2,(Coef_Type) i,sflag,alloc);
    if(debug)
    {
      cout << " *** ZETA1 *** \n";
      for(j=0; j < f1.size()*f2.size();j++)
         (ostream&) cout <<  j << base[ZETA1]->base[j] << endl;

      cout << " *** ZETA2 *** \n";
      for(j=0; j < f1.size()*f2.size();j++)
         (ostream&) cout <<  j << base[ZETA2]->base[j] << endl;

      for(i=ZETA1_BX; i<=XI_AxBZ;i++)
      {
	cout << " *** " << (Coef_Type) i << " *** \n";
	for(j=0; j < f1.size()*f2.size();j++)
	  (ostream&) cout <<  j << base[i]->base[j] << endl;
      }
    }
    break;
  case SPINORB:
    for(i=PX; i<NUCLEAR;i++)
      base[i] = new CoefArray(f1,l1,s1,f2,l2,s2,(Coef_Type) i,sflag,alloc);
    for(i=ZETAa; i<=ZETAab_AxBZ; i++)
      base[i] = new CoefArray(f1,l1,s1,f2,l2,s2,(Coef_Type) i,sflag,alloc);
    if(debug)
    {
      for(i=PX; i<NUCLEAR;i++)
      {
	cout << " *** " << (Coef_Type) i << " *** \n";
	for(j=0; j < f1.size()*f2.size();j++)
	  (ostream&)          cout <<  j <<  base[i]->base[j] << endl;
      }
      for(i=ZETAa; i<=ZETAab_AxBZ;i++)
      {
	cout << " *** " << (Coef_Type) i << " *** \n";
	for(j=0; j < f1.size()*f2.size();j++)
	  (ostream&)  cout <<  j << base[i]->base[j] << endl;
      }
    }
    break;
  }
}

Coefficient_Set::~Coefficient_Set()
{
  for(int i =0; i < size(); i++)
    if (base[i])
      delete base[i];
}









