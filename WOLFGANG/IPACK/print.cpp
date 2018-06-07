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
% $Id: print.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\subsection{Printing Storage Arrays}

The following section prints the possibly symmetrized contents of a
storage array. This subroutine defines the canonical order of the
indices of the symmetrized \name{Storage} arrays.
\index{Storage, order}
\index{Storage, canonical index}
*/ 
#include "io_incl.h"
#include "storage.h"

ostream& operator<<(ostream& os,const Storage& s)
{
  char buf[40];

  cout << " +++ Sizes: " << s.size1() << " " << s.size2() << " " 
    << s.size3() << " " << s.size4() << " Aux: " << s.blocks() << endl;
  int count = 0;

  for(int aux = 0; aux < s.blocks(); aux++)
  {
    os << " +++ Block: " << aux << endl;

    if (s.size3() == 1 && s.size4() == 1) // just 2D array
    {      
      for(int i=0; i<s.size1(); i++)
	for(int j=0; j < s.size2(); j++)
	{
	  sprintf(buf," %20.5e ",s(count++));
	  os << buf;
	  if (j % 5 == 4) os << endl;
	}
      os << endl;
    }
    else   // full 4d array
    {
      int i1,i2,i3,i4;      
      for(i1=0;i1 < s.size1(); i1++)
      for(i2=0;i2 < s.size2(); i2++)
      for(i3=0;i3 < s.size3(); i3++)
      for(i4=0;i4 < s.size4(); i4++)
      {
	sprintf(buf," %2i %2i %2i %2i %20.5e \n",i1,i2,i3,i4,s(count++));
	os << buf;
      }
    }
  }
  return os;
}
