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
% $Id: symdesig.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
\subsection{Integer Indices}
\index{SymDesignator, int operator}

This operator return a integer index to the symmetry adapted basis
functions.
It ss rather tricky, because the number of symmetry adapted
basis functions changes, if symmetry 'none' is set for a direction.
For example:

\begin{tabular}{c|c}
symmetry & index \\ \hline
OON & 0 \\
OEN & 1 \\
EON & 2 \\
EEN & 3 \\ \hline

NON & 0 \\
NEN & 1 \\
\end{tabular}
*/
#include "io_incl.h"
#include <assert.h>
#include "typesip.h"
#include "symdesig.h"
#include "symop.h"
/*
SymDesignator::operator int() const
{ 
  int val = 0;
  int yfac = (SymOp::symmetry(Z)) ? 2      :    1;
  int xfac = (SymOp::symmetry(Y)) ? 2*yfac : yfac;
  if (s[X] == even) val += xfac;
  if (s[Y] == even) val += yfac;
  if (s[Z] == even) val += 1;
  return val;
} 
*/
SymDesignator::operator int() const
{ 
  int val = 0;
  int yfac = (SymOp::symmetry(Z)) ? 2      :    1;
  int xfac = (SymOp::symmetry(Y)) ? 2*yfac : yfac;
  if (s[X] == odd) val += xfac;
  if (s[Y] == odd) val += yfac;
  if (s[Z] == odd) val += 1;
  return val;
} 

/*TEX
\subsection{Output}
*/
ostream& operator<<(ostream& os,const SymType& ss)
{
  switch(ss)
  {
  case none: os << "N"; break;
  case odd : os << "O"; break;
  case even: os << "E"; break;
  }
  return os;
}

ostream& operator<<(ostream& os,const SymDesignator& ss)
{
   return  os << ss(0) << ss(1) << ss(2) << " ";
}





