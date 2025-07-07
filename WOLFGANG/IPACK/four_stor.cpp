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
% $Id: four_stor.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $

\subsection{Implementation}
*/  
#include "io_incl.h"
#include <tarray.h>
#include <assert.h>
#include "typesip.h"
#include "storage.h"
#include "tquadruple.h"
#include "four_stor.h"

void FourStorage::clear()
{
  for(int i=0; i<fields.size(); i++)
    if (fields[i]) 
    {
      delete fields[i];
      fields[i] = 0;
    }
}

void FourStorage::reset(int m1,int m2,int m3,int m4,FOURSTORTYPE ft)
{
  fst = ft;
  Momentum zero;
  
  switch(fst)
  {
  case T0:
    max[0] = m1;
    max[1] = m2;
    max[2] = m3;
    max[3] = m4;
  
    no[0] = zero.number(m1);
    no[1] = zero.number(m2);
    no[2] = zero.number(m3);
    no[3] = zero.number(m4);
    
    break;
  case T1:
    max[0] = 2*m1+1;
    max[1] = m2;
    max[2] = m3;
    max[3] = m4;
  
    no[0] = 2*m1+1;
    no[1] = zero.number(m2);
    no[2] = zero.number(m3);
    no[3] = zero.number(m4);
    
    break;
  case T2:
    max[0] = 2*m1+1;
    max[1] = 2*m2+1;
    max[2] = m3;
    max[3] = m4;
  
    no[0] = 2*m1+1;
    no[1] = 2*m2+1;
    no[2] = zero.number(m3);
    no[3] = zero.number(m4);
    
    break;
  case T3:
    max[0] = 2*m1+1;
    max[1] = 2*m2+1;
    max[2] = 2*m3+1;
    max[3] = m4;
  
    no[0] = 2*m1+1;
    no[1] = 2*m2+1;
    no[2] = 2*m3+1;
    no[3] = zero.number(m4);
    
    break;
  case T4:
    max[0] = 2*m1+1;
    max[1] = 2*m2+1;
    max[2] = 2*m3+1;
    max[3] = 2*m4+1;
  
    no[0] = 2*m1+1;
    no[1] = 2*m2+1;
    no[2] = 2*m3+1;
    no[3] = 2*m4+1;
    
    break;
  }
  
  int nn = no[0]*no[1]*no[2]*no[3];
  
  clear();
  fields.reset(nn);
  fields.set(0);

}

ostream& operator<<(ostream& os,const FourStorage& fs)
{
  return os << "Printing FourStorage Not Implemented. " << endl;
}

