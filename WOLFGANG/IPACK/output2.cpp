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
% $Id: output2.cpp,v 1.2 2008/05/03 06:23:55 igor Exp $
%
\subsection{Two-Electron Integrals}
  
The following function ouputs the two-electron from the \name{IntegArray}.
The $l_i$ are the angular momenta on the four indices, the basis-set
arguments contain the normalization factors and the information encoding 
the ultimate integral-indices.

Suppose there are $m_i$ integrals on the i-th index. The order of the
indices in the integral array is: (l4,l3,l2,l1) because of the fortran
compatible storage mode on matrices. 

A function which carries the index $i$ in the \name{IntegArray}, carries the 
orbital-index: \name{startindex[l][i] + moffset}, where $l$ is the total 
angular momentum of the function, \name{startindex[l][i]} is thus the 
orbital-index of the first function in the angular momentum set associated 
with function $i$ and \name{moffset} is the offset of the orbital induces by
the angular momentum label. 
*/

#include "parallel.h"
#include <assert.h>
#include <math.h>
#include "io_incl.h"
#include "const.h"
#include "function.h"
#include "basis.h"
#include "integ_file.h"
#include "storage.h"
#include "symdesig.h"


void Storage::output(InternalBasis& b1,SymDesignator& sd1,
		     InternalBasis& b2,SymDesignator& sd2,
		     InternalBasis& b3,SymDesignator& sd3,
		     InternalBasis& b4,SymDesignator& sd4,
		     int l1,int moff1,int l2,int moff2,
		     int l3,int moff3,int l4,int moff4,
		     IntegFile& obuf,int block_no,int same12,int same34,int samepair)
{
  // cout << " Output for block: " << block_no << endl << " "  << *this << endl;
  
  if (same12 && l1 == l2 && sd1 == sd2 && moff1 < moff2)
    return;
  if (same34 && l3 == l4 && sd3 == sd4 && moff3 < moff4)
    return;
  if (samepair && l1 == l3 && l2 == l4 && sd1 == sd3 && moff1 < moff3)
    return;
  if (samepair && l1 == l3 && l2 == l4 && sd1 == sd3 && moff1 == moff3 
      && sd2 == sd4 && moff2 < moff4)
    return;

  // the start-indices for the angular momenta

  IVector& s1 = b1.s_index[l1];  
  IVector& s2 = b2.s_index[l2];  
  IVector& s3 = b3.s_index[l3];  
  IVector& s4 = b4.s_index[l4];  

  // the number of functions

  int max1 = b1.nofun[l1];
  int max2 = b2.nofun[l2];
  int max3 = b3.nofun[l3];
  int max4 = b4.nofun[l4];

  int global1 = b1.sym_start_index(sd1);
  int global2 = b2.sym_start_index(sd2);
  int global3 = b3.sym_start_index(sd3);
  int global4 = b4.sym_start_index(sd4);
  
  int i1,i2,i3,i4,count = 0;
  
  int blk_size  = size() / blocks();
  // double *array = base + block_no * blk_size;
  
  for(i1=0;i1 < max1; i1++)
  {
    int    oi1 = s1[i1] + moff1 + global1;
    double n1  = b1.norm[oi1];
    oi1       += b1.min_index();
    
    for(i2=0;i2 < max2; i2++)
    {
      int    oi2 = s2[i2] + moff2 + global2;
      double n2  = n1*b2.norm[oi2];
      oi2       += b2.min_index();
     
      for(i3=0;i3 < max3; i3++)
      {
	int    oi3 = s3[i3] + moff3 + global3;
	double n3  = n2*b3.norm[oi3];
	oi3       += b3.min_index();

	for(i4=0;i4 < size4(); i4++)
	{
	  int    oi4   = s4[i4] + moff4 + global4;
	  double norm  = n3*b4.norm[oi4];
	  oi4         += b4.min_index();

	  double vv = base[block_no * blk_size + count];
	  /*
	    cout << "O:" 
	       <<  setw(4) << oi1 << setw(4) << oi2 << setw(4) << oi3 <<
	    setw(4) << oi4 << " " << *array * norm << endl;
	  */
	  if (same12 && l1 == l2 && moff1 == moff2 && oi2 < oi1) 
	  { count++; continue; } 
	  if (same34 && l3 == l4 && moff3 == moff4 && oi4 < oi3) 
	  { count++; continue;}  
	  if (samepair && l1 == l3 && moff1 == moff3 && l2 == l4 && 
	      moff2 == moff4 && oi3 <  oi1) 
	  { count++; continue; } 
	  if (samepair && l1 == l3 && moff1 == moff3 && l2 == l4 && 
	      moff2 == moff4 && oi1 == oi3 && oi4 < oi2) 
	  { count++; continue; }
	  //cout << "X:" 
	  // <<  setw(4) << oi1 << setw(4) << oi2 << setw(4) << oi3 <<
	  //setw(4) << oi4 << " " << *array * norm << endl;
          obuf.dump(oi1,oi2,oi3,oi4,vv * norm);
	  count++;
	}
      }
    }
  }
}
