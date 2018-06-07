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
% $Id: tpair.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
\subsection{Pairs of Objects}  
A simple class to store pairs of simple objects and to perform
the symmetry operations typical for single-electron integrals.
*/
#ifndef T_PAIR_H
#define T_PAIR_H

template <class T> class Pair
{
  T p[2];
public:
      Pair() { }   
      Pair(T& e1,T& e2) { p[0] = e1; p[1] = e2; } 
     ~Pair() { }	
void  permute()   { T tmp(p[0]);   p[0] = p[1];  p[1] = tmp;  }

int compare(const Pair<T>& p2) 
                 {  
		   for(int i=0; i< 2; i++)
		     if (p[i] != p2(i))
		       return (p[i] < p2(i)) ? -1 : 1;
		   return 0;
		 }
T&       operator[](int i)       { assert(i >=0 && i < 2);  return p[i]; } 
const T& operator()(int i) const { assert(i >=0 && i < 2);  return p[i]; } 
int      operator< (const Pair<T>& qq) { return (compare(qq) == -1); } 
int      operator!=(const Pair<T>& qq) { return (compare(qq) !=  0); } 
int      operator==(const Pair<T>& qq) { return (compare(qq) ==  0); } 
};

template <class T> ostream& operator<<(ostream& os,const Pair<T> q)
{
  for(int i=0; i<2; i++)
    os << q(i);
  return os;
}

#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
template <class T> ostream_withassign& operator<<(ostream_withassign& os,
						  const Pair<T> q)
{
  for(int i=0; i<2; i++)
    os << q(i);
  return os;
}
#endif


#endif
