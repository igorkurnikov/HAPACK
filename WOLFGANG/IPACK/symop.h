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
% $Id: symop.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
\section{Symmetry Operators}
\index{SymOp}
\index{SymDirOp}  
\index{symflag}
\index{same}
The class \name{SymOp} is created to perform the symmetry operations
on a basis set. We only deal with reflection symmetries in each of the
three cartesian directions. 

Starting \name{ipack} the user can define
the reflection symmetries of the molecule he wants to calculate. They
are stored globally in the static \name{symflag} vector.

With respect to each direction, a basis can be decorated with
either the identity or the reflection operator, provided that the
basis is not invariant under the application of this operator. For a
basis set which is invariant under the application of a given
reflection e. g. when an atom lies on teh reflection-plane,
the only allowed operator for this direction is unity.  
This is also the case, if the molecule has no symmetry in this
direction. Therefore we use \name{same} to indicate whether reflection
is allowed or not. 

For loops over all possible symmetry operations, we define an ordering
like
\be
111\\11R\\1R1\\1RR\\R11\\R1R\\RR1\\RRR
\ee
(if reflection is possible in each direction) and use legal for range
checking.

The combination of two symmetry operations can be achieved by the ^ operator.
*/
#ifndef SYMOP_H
#define SYMOP_H

enum SymDirOp { one, reflection };     

inline SymDirOp operator^(SymDirOp a,SymDirOp b)
{
  if (a == b) return one;
  return reflection;
}

ostream& operator<<(ostream& os,SymDirOp op);

class SymOp
{
static int symflag[3];
  int      legal;  
  int      same[3];
  SymDirOp sym[3];
  
public:
  SymOp() { legal = false; } 
  SymOp(int symx,int symy,int symz) 
       { 
       // set same = 1 if reflection is not allowed
	 same[X] = symx || !symflag[X]; 
	 same[Y] = symy || !symflag[Y]; 
	 same[Z] = symz || !symflag[Z]; 
	 sym[X] = sym[Y] = sym[Z] = one; 
	 legal = true;
       } 
  ~SymOp() {}
static void set_symmetry(Direction d) { symflag[d] = true; }  
static int  symmetry(Direction d) { return symflag[d]; }    

// number of symmetries
static int  no_symmetry() { return symflag[X] + symflag[Y] + symflag[Z]; } 
  
int operator++()
       {
	 if (sym[2] == one
	     && (!same[2])) { sym[2] = reflection;  return 1; } 
	 if (sym[1] == one && !same[1]) { sym[1] = reflection; 
					 sym[2] = one;     return 1; } 
	 if (sym[0] == one && !same[0]) { sym[0] = reflection; 
					 sym[1] = sym[2] = one; return 1; }
	 legal = false;
	 return 0;
       }
int  in_range() const { return legal; } 	 
void reflect(Direction& dir)
     {
       if (!same[dir] && sym[dir] == one) { sym[dir] = reflection; return; } 
       sym[dir] = one;
     }
int  operator!=(const SymOp& b) 
     { return (sym[0] != b.sym[0] || sym[1] != b.sym[1] || sym[2] !=
	       b.sym[2]); } 
int  operator<(const SymOp& b) 
     { if (sym[0] < b.sym[0]) return 1; if (sym[0] > b.sym[0]) return 0;
       if (sym[1] < b.sym[1]) return 1; if (sym[1] > b.sym[1]) return 0;
       if (sym[2] < b.sym[2]) return 1; return 0;
     }
SymDirOp& operator[](int i) { return sym[i]; } 
const SymDirOp& operator()(int i) const { return sym[i]; } 
int   factor(Direction d) const { return (sym[d] == one) ? 1 : -1; } 
      operator int() const { return 4*sym[0]+2*sym[1]+sym[2]; } 
SymOp operator^(const SymOp& b) const; 
};

ostream& operator<<(ostream& os,const SymOp& so);

#endif
