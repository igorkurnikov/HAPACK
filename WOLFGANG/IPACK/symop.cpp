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
% $Id: symop.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
\subsection{Declaration of symflag}
\index{symflag, declaration}
*/
#include "io_incl.h"
#include "typesip.h"
#include "symop.h"

int SymOp::symflag[3] = { 0,0,0 };
 
/*TEX
\subsection{Combined Symmetry Operations}
\index{SymOp, ^ operator}

The following operator returns the symmetry operations that is
equivalent to performing the connected symmetry operations one after
the other. In other words the result is just the second symmetry
operation applied to the first one, vice versa.
*/
SymOp SymOp::operator^(const SymOp& b) const 
{
  SymOp tmp;
  tmp.legal   = legal && b.legal;
  tmp.same[0] = same[0] && b.same[0];
  tmp.same[1] = same[1] && b.same[1];
  tmp.same[2] = same[2] && b.same[2];

  tmp.sym[0]  = sym[0] ^ b.sym[0];
  tmp.sym[1]  = sym[1] ^ b.sym[1];
  tmp.sym[2]  = sym[2] ^ b.sym[2];
  
  return tmp;
}

/*TEX
\subsection{Output}
\index{SymOp, output}
*/
ostream& operator<<(ostream& os,SymDirOp op)
{
  switch(op)
  {
  case one: os << '1';
    break;
  case reflection: os << 'R';
    break;
  }
  return os;
} 

ostream& operator<<(ostream& os,const SymOp& so)
{
  for(int i=0; i< 3; i++)
    os << so(i);
  return os << " ";  
}
