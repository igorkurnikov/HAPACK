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
% $Id: symdesig.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
\section{Symmetry Designators}
\index{SymType}
\index{SymDesignator}
\index{same} 

The class \name{SymDesignator} is used to label symmetry adapted basis
functions. For each coordinate 'none' (N), 'odd' (O) or 'even' (E)
reflection symmetry can be specified.
We indicate with \name{same} in which directions 'even and 'odd'
symmetry are not
meaningfull, because the  molecule is not symmetric (\name{symflag})
or the basis function is evaluated on the symmetry-plane.

In order to simply define loops over \name{SymDesignator}, we
introduce an ordering:
\be
OOO < OOE <OEO < OEE <EOO <EOE < EEO < EEE
\ee
This reduces of course if 'none' symmetry is set for a coordinate.
\name{legal} is true if the \name{SymDesignator} is still in an
allowed range.

Additionally we want to label the symmetry adapted basis functions
with integer indeces. This is done by the \name{int} operator.
*/
 
#ifndef SYMDESIG_H
#define SYMDESIG_H

//enum SymType  { none, even, odd };   
enum SymType  { none, odd, even };   
ostream& operator<<(ostream& os,const SymType& ss);
  
class SymDesignator
{
  int     legal;
  int     same[3];
  SymType s[3];
public:
     SymDesignator() { legal = false; } 
     SymDesignator(int sx,int sy,int sz)
       {
	 legal = true;
	 same[X] = sx || !SymOp::symmetry(X); 
	 same[Y] = sy || !SymOp::symmetry(Y); 
	 same[Z] = sz || !SymOp::symmetry(Z); 

	 s[0] = (same[0]) ? none : odd;
	 s[1] = (same[1]) ? none : odd;
	 s[2] = (same[2]) ? none : odd;
       }
    ~SymDesignator() {}
void operator++()
{
  if (s[2] == odd) { s[2] = even; return; }
  if (s[1] == odd) { s[1] = even; if (s[2] != none) s[2] = odd; return; }
  if (s[0] == odd) { s[0] = even; if (s[2] != none) s[2] = odd; 
		     if (s[1] != none) s[1] = odd; return; }
  legal = false;
}

int  in_range() const { return legal; }   
friend ostream& operator<<(ostream& os,const SymDesignator& ss);
int  operator!=(const SymDesignator& b) 
     { return (s[0] != b.s[0] || s[1] != b.s[1] || s[2] != b.s[2]); } 
int  operator==(const SymDesignator& b) 
     { return (s[0] == b.s[0] && s[1] == b.s[1] && s[2] == b.s[2]); } 

int  operator<(const SymDesignator& b) 
     { 
       if (s[0] < b.s[0]) return 1; if (s[0] > b.s[0]) return 0;
       if (s[1] < b.s[1]) return 1; if (s[1] > b.s[1]) return 0;
       if (s[2] < b.s[2]) return 1; return 0;
     }

     operator int() const;
SymType& operator[](int i) { assert(i >=0 && i < 3);  return s[i]; } 
const SymType& operator()(int i) const { assert(i >=0 && i < 3);  return s[i]; } 
}; 


#endif
