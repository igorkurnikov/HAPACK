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
% $Id: transpose.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\subsection{Generic Permutations}
 
The following routine implements a generic permutation of the indices in a 
RECTANGLUAR storage array. The type of transformation is encoded in
the \name{trans} variable, which takes the various constants of type \name{TRANS}.
*/

#include "io_incl.h"
#include <assert.h>
#include "storage.h"

inline void  swap(int& i,int& j)
{
  int tmp = i;
  i = j;
  j = i;
}

Storage::Storage(Storage& source,int trans)
{
  assert(source.sym12 == 0);
  assert(source.sym34 == 0);
  assert(source.sympair == 0);
  c_number++;

  
  Storage target = *this;
  
  int p[4];
  for(int i=0; i < 4; i++)
    p[i] = i;

  if (trans & TRANS12)
    swap(p[0],p[1]);
  if (trans & TRANS34)
    swap(p[2],p[3]);
  if (trans & TRANSPAIR)
  {
    swap(p[0],p[2]);
    swap(p[1],p[3]);
  }
  
  reset(source.sz[p[0]],source.sz[p[1]],false,source.sz[p[2]],source.sz[p[3]],false,
	false,source.sz[4]);

  double *array = base;
  
  for(int block = 0; block < sz[4]; block++)
  {
    double *src = source.base + block*sz[0]*sz[1]*sz[2]*sz[3];
    double *tar =        base + block*sz[0]*sz[1]*sz[2]*sz[3];
    int i[4];
    
    for(i[0]=0; i[0] < source.sz[0]; i[0]++)
    for(i[1]=0; i[1] < source.sz[1]; i[1]++)
    for(i[2]=0; i[2] < source.sz[2]; i[2]++)
    for(i[3]=0; i[3] < source.sz[3]; i[3]++)
    {
      int newpos  = ((i[p[0]]*sz[0] + i[p[1]])*sz[1]+i[p[2]])*sz[2] + i[p[3]];
      tar[newpos] = *src++;
    }
  }
}
