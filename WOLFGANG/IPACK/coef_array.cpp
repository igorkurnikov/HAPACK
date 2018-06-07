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
% $Id: coef_array.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\subsection{Field-Names\label{fieldnames}} 
*/ 
#include "io_incl.h"
#include <errorip.h>
#include <math.h>
#include "typesip.h"
#include "function.h"
#include "coef_array.h"
#include "coef_set.h"
#include "const.h"
#include "symop.h"

const char* Coef_Name[] = {
     "ZETASUM","PA1X","PA1Y","PA1Z","PA2X","PA2Y","PA2Z","NORM","OVERLAP",
     "ZETA1","ZETA2","XI","KINETIC",
     "PX","PY","PZ","NUCLEAR",                      // MOMENT
     "ZETASUM_PX","ZETASUM_PY","ZETASUM_PZ",        // E_FIELD
     "ZETA1_BX","ZETA1_BY","ZETA1_BZ",
     "ZETA2_AX","ZETA2_AY","ZETA2_AZ",          
     "XI_AminusBX","XI_AminusBY","XI_AminusBZ",
     "XI_AxBX","XI_AxBY","XI_AxBZ",                 // ANGULAR
     "ZETAa","ZETAa_AX","ZETAa_AY","ZETAa_AZ", 
     "ZETAb","ZETAb_BX","ZETAb_BY","ZETAb_BZ",
     "ZETAab_AminusBX","ZETAab_AminusBY","ZETAab_AminusBZ",
     "ZETAab_AxBX","ZETAab_AxBY","ZETAab_AxBZ",     // SPINORB
     "ABX","ABY","ABZ",
     "RHO1","RHO2","RHOSUM","WX","WY","WZ"          //TWO_EL
   };

ostream& operator<<(ostream& os,Coef_Type ct)
{
  return cout << Coef_Name[ct];
}

/*TEX
\subsection{Initializers \name{CoefArrays}} 

The following function initializes single-particle
\name{CoefficientArrays}.  First the underlying \name{Storage} array
is reset.  Depending of the symmetry-flag \name{same} we then run over
all symmetric or non-symmetric combinations of indices $i,j$ labeling
elements of the \name{RadialFunctionArrays} in the argument. For each
such index-pair we compute the appropriate element of the \name{CoefArray}. 
*/

CoefArray::CoefArray(RadialFunctionArray& f1,int l1,SymOp& s1,
		     RadialFunctionArray& f2,int l2,SymOp& s2,
		     Coef_Type ct,int same,Allocator& alloc) : Storage(alloc) 
{
  int i;
  int j;
  // int cnt = 0;
  
  int mx1 = f1.size();
  int mx2 = f2.size();
  assert(mx1 == size1());
  assert(mx2 == size2());
  assert(size3() == 1);
  assert(size4() == 1);
  assert(blocks() == 1);
  
  //  reset(mx1,mx2,1,1,1);

  double xi,dist,fac;
  Direction d;
  
  // int sym   = symmetry12();
  int count = 0;

  int sym1[3];
  int sym2[3];

  sym1[X] = s1.factor(X);
  sym1[Y] = s1.factor(Y);
  sym1[Z] = s1.factor(Z);

  sym2[X] = s2.factor(X);
  sym2[Y] = s2.factor(Y);
  sym2[Z] = s2.factor(Z);
  
  switch(ct)
  {
  case ZETASUM:
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
      base[count++] =  0.5/(f1[i].exp()+f2[j].exp());
    break;
  case PA1X:
  case PA1Y:
  case PA1Z:
    d = (Direction) (ct - PA1X);
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] =  ( sym1[d]*f1[i].exp()*f1[i](d) 
			  + sym2[d]*f2[j].exp()*f2[j](d))/
			    (f1[i].exp() + f2[j].exp()) 
			      - sym1[d] * f1[i](d);
    break;
  case PA2X:
  case PA2Y:
  case PA2Z:
    d = (Direction) (ct - PA2X);
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] =  (sym1[d]*f1[i].exp()*f1[i](d) +
			  sym2[d]*f2[j].exp()*f2[j](d))/
			    (f1[i].exp() + f2[j].exp()) 
			      - sym2[d]*f2[j](d);
    break;
  case NORM:
    for(i = 0 ; i < mx1; i++)
    { 
      double n1 = f1[i].norm(l1);
      for(j=0 ; j <  mx2; j++)
      	base[count++] = n1 * f2[j].norm(l2) - 1;
    }  
    break;
  case OVERLAP:
    
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
    {
      xi   = (f1[i].exp()*f2[j].exp())/(f1[i].exp()+f2[j].exp());
      dist = ( sym1[X]*f1[i](X) - sym2[X]*f2[j](X) ) 
	   * ( sym1[X]*f1[i](X) - sym2[X]*f2[j](X) ) +
	     ( sym1[Y]*f1[i](Y) - sym2[Y]*f2[j](Y) ) 
	   * ( sym1[Y]*f1[i](Y) - sym2[Y]*f2[j](Y) ) +
	     ( sym1[Z]*f1[i](Z) - sym2[Z]*f2[j](Z) ) 
	   * ( sym1[Z]*f1[i](Z) - sym2[Z]*f2[j](Z) ); 
      fac  = pi / (f1[i].exp()+f2[j].exp());
      base[count++] =  exp(-xi*dist) * fac * sqrt(fac);
    }
    break;
  case ZETA1:
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] =  f2[j].exp()/(f1[i].exp()+f2[j].exp());
    break;
  case ZETA2:
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
      base[count++] =  f1[i].exp()/(f1[i].exp()+f2[j].exp());
    break;
  case XI:
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
      base[count++] = 2*(f1[i].exp()*f2[j].exp())/(f1[i].exp()+f2[j].exp());
    break;
  case KINETIC:
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
      {
        xi   = (f1[i].exp()*f2[j].exp())/(f1[i].exp()+f2[j].exp());
        dist = ( sym1[X]*f1[i](X) - sym2[X]*f2[j](X) ) 
	     * ( sym1[X]*f1[i](X) - sym2[X]*f2[j](X) ) +
	       ( sym1[Y]*f1[i](Y) - sym2[Y]*f2[j](Y) ) 
	     * ( sym1[Y]*f1[i](Y) - sym2[Y]*f2[j](Y) ) +
	       ( sym1[Z]*f1[i](Z) - sym2[Z]*f2[j](Z) ) 
	     * ( sym1[Z]*f1[i](Z) - sym2[Z]*f2[j](Z) );
      
        fac  = pi / (f1[i].exp()+f2[j].exp());
        base[count++] =  xi*(3-2*xi*dist)* exp(-xi*dist) * fac * sqrt(fac);
      }
    break;
  case PX:
  case PY:
  case PZ:
    d = (Direction) (ct - PX);
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
      base[count++] =  (sym1[d]*f1[i].exp()*f1[i](d) +
			sym2[d]*f2[j].exp()*f2[j](d))/
			  (f1[i].exp() + f2[j].exp());
    break;
  case ZETASUM_PX:
  case ZETASUM_PY:
  case ZETASUM_PZ:
    d = (Direction) (ct - ZETASUM_PX);
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
      base[count++] =  2 * (f1[i].exp()+f2[j].exp())
                       *(sym1[d]*f1[i].exp()*f1[i](d) +
			sym2[d]*f2[j].exp()*f2[j](d))/
			  (f1[i].exp() + f2[j].exp());
    break;
  case ZETA1_BX:
  case ZETA1_BY:
  case ZETA1_BZ:
    d = (Direction) (ct - ZETA1_BX);
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] = f2[j].exp()/(f1[i].exp()+f2[j].exp())
                        * sym2[d] * f2[j](d);
    break;
  case ZETA2_AX:
  case ZETA2_AY:
  case ZETA2_AZ:
    d = (Direction) (ct - ZETA2_AX);
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] = f1[i].exp()/(f1[i].exp()+f2[j].exp())
                        * sym1[d] * f1[i](d);
    break;
  case XI_AminusBX:
  case XI_AminusBY:
  case XI_AminusBZ:
    d = (Direction) (ct - XI_AminusBX);
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] = 2*(f1[i].exp()*f2[j].exp())/(f1[i].exp()+f2[j].exp())
                        * (sym1[d] * f1[i](d) - sym2[d] * f2[j](d));
    break;
  case XI_AxBX:
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] = 2*(f1[i].exp()*f2[j].exp())/(f1[i].exp()+f2[j].exp())
                        * (sym1[Y] * f1[i](Y) * sym2[Z] * f2[j](Z)
			   - sym1[Z] * f1[i](Z) * sym2[Y] * f2[j](Y));
    break;
  case XI_AxBY:
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] = 2*(f1[i].exp()*f2[j].exp())/(f1[i].exp()+f2[j].exp())
	                * (sym1[X] * f1[i](X) * sym2[Z] * f2[j](Z)
			   - sym1[Z] * f1[i](Z) * sym2[X] * f2[j](X));
    break;
  case XI_AxBZ:
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] = 2*(f1[i].exp()*f2[j].exp())/(f1[i].exp()+f2[j].exp())
                        * (sym1[X] * f1[i](X) * sym2[Y] * f2[j](Y)
			   - sym1[Y] * f1[i](Y) * sym2[X] * f2[j](X));
    break;
  case ZETAa:
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
      base[count++] =  2 * f1[i].exp();
    break;
  case ZETAa_AX:
  case ZETAa_AY:
  case ZETAa_AZ:
    d = (Direction) (ct - ZETAa_AX);
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] =  2 * f1[i].exp() * sym1[d] * f1[i](d);
    break;
  case ZETAb:
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
      base[count++] =  2 * f2[j].exp();
    break;
  case ZETAb_BX:
  case ZETAb_BY:
  case ZETAb_BZ:
    d = (Direction) (ct - ZETAb_BX);
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] =  2 * f2[j].exp() * sym2[d] * f2[j](d);
    break;
  case ZETAab_AminusBX:
  case ZETAab_AminusBY:
  case ZETAab_AminusBZ:
    d = (Direction) (ct - ZETAab_AminusBX);
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] = 4*(f1[i].exp()*f2[j].exp())
                        * (sym1[d] * f1[i](d) - sym2[d] * f2[j](d));
    break;
  case ZETAab_AxBX:
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] = 4*(f1[i].exp()*f2[j].exp())
                        * (sym1[Y] * f1[i](Y) * sym2[Z] * f2[j](Z)
			   - sym1[Z] * f1[i](Z) * sym2[Y] * f2[j](Y));
    break;
  case ZETAab_AxBY:
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] = 4*(f1[i].exp()*f2[j].exp())
	                * (sym1[X] * f1[i](X) * sym2[Z] * f2[j](Z)
			   - sym1[Z] * f1[i](Z) * sym2[X] * f2[j](X));
    break;
  case ZETAab_AxBZ:
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < mx2; j++)
	base[count++] = 4*(f1[i].exp()*f2[j].exp())
                        * (sym1[X] * f1[i](X) * sym2[Y] * f2[j](Y)
			   - sym1[Y] * f1[i](Y) * sym2[X] * f2[j](X));
    break;    
  default:
    error("Illegal Code in Coef_Array init. ");
  }
}

/*TEX
\subsection{Initializers for contracted function\name{CoefArrays}} 

This section deals with the initialization of the \name{CoefArrays} for 
contracted orbitals. The method is analogous to that for uncontracted orbitals.

*/
CoefArray::CoefArray(LocationArray& f1,SymOp& s1,LocationArray& f2,SymOp& s2,
		     Coef_Type ct,int same,Allocator& alloc) : Storage(alloc)
     
{
  int i;
  int j;
  // int cnt = 0;
  
  int mx1 = f1.size();
  int mx2 = f2.size();

  assert(mx1 == size1());
  assert(mx2 == size2());
  assert(size3() == 1);
  assert(size4() == 1);
  assert(blocks() == 1);
  
  //  reset(mx1,mx2,1,1,1);

//  double fac;
  Direction d;
  
  int count = 0;
  
  switch(ct)
  {
  case ABX:
  case ABY:
  case ABZ:
    d =  (Direction) (ct-ABX);
    for(i = 0 ; i < mx1; i++)
      for(j=0 ; j < ((same) ? (i+1) : mx2); j++)
	base[count++] =  s1.factor(d)*f1[i](d) - s2.factor(d)*f2[j](d);
    break;
  default:
    error("Illegal Code in Coef_Array init. ");
  }
}

/*TEX
\subsection{Initialization of Four-Index \name{CoefArrays}}
A number of four-index arrays are required in the evaluation of the 
vertical recursion relation. Here we show the formulas, how these 
\name{CoefArrays} can be computed from the \name{Coefficient\_Sets} on either
index.

Let $\tilde\zeta_{s1}$ and $\tilde\zeta_{s2}$ label the $\tilde\zeta_s$ \name{CoefArrays} in the first and second \name{Coefficient\_Set} respectively: 
\ba 
\rho   & = & \frac{\zeta_{12} \zeta_{34}}{\zeta_{12} + \zeta_{34}}   \\
\rho_1 & = & \frac{\zeta_{34}}{2\zeta_{12} (\zeta_{12} + \zeta_{34})}  
             \nonumber \\
       & = & \frac{\tilde{\zeta_{s1}}^2 }
                  {\tilde{zeta_{s1}} + \tilde{\zeta_{s2}}}\\
\rho_2 & = & \frac{\zeta_{12}}{2\zeta_{34} (\zeta_{12} + \zeta_{34})}   
             \nonumber \\
       & = & \frac{\tilde{\zeta_{s2}}^2 }
                  { 2(\tilde{\zeta_{s1}} + \tilde{\zeta_{s2}})}\\
\rho_s & = & \frac{1}{2(\zeta_{12}+\zeta_{34})}  \nonumber \\
       & = & \frac{\tilde{\zeta_{1s}}\tilde{\zeta_{2s}}}
                  {\tilde{\zeta_{1s}} + \tilde{\zeta_{2s}}}
\ea

For the remaining \name{CoefArrays} we obtain
\ba
  \vec{W} & = & \frac{ \zeta_{12} \vec{P}_{12} + \zeta_{34} \vec{P}_{34}}
                     { \zeta_{12} + \zeta{34} } \nonumber \\
          & = & \frac{ \tilde\zeta_{34} \vec{P}_{12} 
                     + \tilde\zeta_{12} \vec{P}_{34}}
                     { \tilde\zeta_{12} + \tilde\zeta{34} },
\ea
where $\vec{P}_{12}$ labels the $\vec{P}$ of the first and 
$\vec{P}_{34}$ labels the $\vec{P}$ of the second \name{Coefficient\_Set}.

The indices $i$ and $j$ here refer to the compound indices, each
labeling a pair of \name{RadialFunctions} encoded in the coefficient
arrays from the two \name{Coefficient\_Sets} passed as arguments. The
new \name{CoefArray} inherits the symmetry poroperties of the
underlying single-particle \name{CoefArrays}. 

Note: OS refers to $\vec{P}_{12}$ as $\vec{P}$ and to $\vec{P}_{34}$ 
as $\vec{Q}$.

*/
CoefArray::CoefArray(Coefficient_Set& cs1,Coefficient_Set& cs2,
		     Coef_Type ct,int same,Allocator& alloc) : Storage(alloc) 
{
  int i,j;

  assert(cs1[ZETASUM].size1() == size1());
  assert(cs1[ZETASUM].size2() == size2());
  assert(cs2[ZETASUM].size1() == size3());
  assert(cs2[ZETASUM].size2() == size4());
  assert(blocks() == 1);

  //reset(cs1[ZETASUM].size1(),cs1[ZETASUM].size2(),
  //      cs2[ZETASUM].size1(),cs2[ZETASUM].size2(),1);
    
  IntegType* zeta1 = cs1[ZETASUM].get_base(); // really 0.5 / zeta1
  IntegType* zeta2 = cs2[ZETASUM].get_base();
  
  IntegType* p1[3];
  IntegType* p2[3];
  
  p1[X] = cs1[PX].get_base();
  p1[Y] = cs1[PY].get_base();
  p1[Z] = cs1[PZ].get_base();
  p2[X] = cs2[PX].get_base();
  p2[Y] = cs2[PY].get_base();
  p2[Z] = cs2[PZ].get_base();
  Direction d;
  IntegType *pp1,*pp2;

  int mx1 = size12();
  int mx2 = size34();

  int cnt = 0;
  
  switch(ct)
  {
  case WX:
  case WY:
  case WZ:
    d = (Direction) (ct - WX);
    pp1 = p1[d];
    pp2 = p2[d];
    for(i=0;i< mx1;i++)
    for(j=0;j< mx2;j++)
	base[cnt++] = (zeta2[j]*pp1[i] + zeta1[i]*pp2[j]) /
	  (zeta1[i] + zeta2[j]);
    break;
  case RHO1:
    for(i=0;i < mx1;i++)
    for(j=0;j < mx2;j++)
	base[cnt++] = - zeta1[i] * zeta1[i] / (zeta1[i] + zeta2[j]);
    break;    
  case RHO2:
    for(i=0;i < mx1;i++)
    for(j=0;j < mx2;j++)
	  base[cnt++] = - zeta2[j] * zeta2[j] / (zeta1[i] + zeta2[j]);
    break;    
  case RHOSUM:
    for(i=0;i < mx1;i++)
    for(j=0;j < mx2;j++)
	base[cnt++] = zeta1[i] * zeta2[j] / (zeta1[i] + zeta2[j]);
    break;    
  default:
    error("Illegal Code in Coef_Array init for quadarrays.");
  }
}

