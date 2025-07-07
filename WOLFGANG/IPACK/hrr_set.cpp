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
% $Id: hrr_set.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\subsection{Constructor for \name{HRR\_Set}}
The constructor for the \name{HRR\_Set} calls the initializers of the 
three required \name{CoefArrays} for the horizontal recursion relations. 
The initializers for  the individual \name{CoefArrays} are used to fill 
the fields.  
*/
#include "io_incl.h"
#include <errorip.h>
#include <math.h>
#include "typesip.h"
#include "function.h"
#include "coef_array.h"
#include "hrr_set.h"

HRR_Set::HRR_Set(LocationArray& f1,SymOp& s1,LocationArray& f2,SymOp& s2,
		 int sflag) : alloc(f1.size(),f2.size(),1,1,1) 
{
  sym = sflag;
  reset(3);  
  for(int i=0; i<3;i++)
    base[i] = new CoefArray(f1,s1,f2,s2,(Coef_Type) (ABX+i),sflag,alloc);
}

HRR_Set::~HRR_Set()
{
  for(int i=0; i<3;i++)
    delete base[i];
}
