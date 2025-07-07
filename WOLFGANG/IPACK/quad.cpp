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
% $Id: quad.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\subsection{Constructor}
The two coefficient sets contain the pair info for the (AB) and the (CD) 
blocks respectively:

Here we compute six new arrays: \name{WX,WY,WZ} and
\ba     
     \rho_1 & = & -\frac{\rho}{2 \zeta_{12}^2} \\
     \rho_1 & = & -\frac{\rho}{2 \zeta_{34}^2} \\
     \rho_s & = & \frac{1}{2 (\zeta_{12} + \zeta_{34})}
\ea

References to these arrays should be retrieved through the function
\name{obtain}.

*/
#include "io_incl.h"
#include <errorip.h>
#include <math.h>
#include "typesip.h"
#include "function.h"
#include "quad.h"

Quadruple_Set::Quadruple_Set(Coefficient_Set& cs1,Coefficient_Set& cs2, 
			     int sflag) :   cset1(cs1), cset2(cs2),
     alloc(cs1[ZETASUM].size1(),cs1[ZETASUM].size2(),
	   cs2[ZETASUM].size1(),cs2[ZETASUM].size2(),1)

{
  sym = sflag;
  reset(6);
  
  for(int d=0; d<3;d++)
    base[3+d] = new CoefArray(cs1,cs2,(Coef_Type)  (WX+d),sflag,alloc);

  base[0] = new CoefArray(cs1,cs2,  RHO1,sflag,alloc);
  base[1] = new CoefArray(cs1,cs2,  RHO2,sflag,alloc);
  base[2] = new CoefArray(cs1,cs2,RHOSUM,sflag,alloc);
}

Quadruple_Set::~Quadruple_Set()
{
  for(int d=0; d<6;d++)
    delete base[d];
}

/*TEX
\subsection{Retrieving \name{CoefArrays}}

The function \name{obtain} allows th retrieval of both single and
two-particle name{CoefArrays} from a \name{Quadruple\_Set}. See
section~\ref{qset_def} for a description.
*/
CoefArray& Quadruple_Set::obtain(const Coef_Type& ct,int pos)
{
  switch(ct)
  {
  case WX:
  case WY:
  case WZ:
  case RHO1:
  case RHO2:
  case RHOSUM:
    assert(pos == -1);
    return *base[ct-RHO1];
  default:
    switch(pos)
    {
    case 0:
      return cset1[ct];
    case 1:
      return cset2[ct];
    default:
      error("illegal position in obtain.");
    }
  }
  return cset1[ZETA1]; // dummy return statement cant be reached.
}
