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
% $Id: tquadruple.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
\subsection{Quadruples}  
A simple class to store quadruples of simple objects and to perform
the symmetry operations typical for 2-electron integrals.
*/
#ifndef T_QUADRUPLE_H
#define T_QUADRUPLE_H

template <class T> class Quadruple
{
  T p[4];
public:
      Quadruple() { }   
      Quadruple(T& e1,T& e2,T& e3,T& e4) { p[0] = e1; p[1] = e2; 
					   p[2] = e3; p[3] = e4; } 
     ~Quadruple() {}
void  perm12()   { T tmp(p[0]);   p[0] = p[1];  p[1] = tmp;  }
void  perm34()   { T tmp(p[2]);   p[2] = p[3];  p[3] = tmp;  } 
void  permpair() { T tmp1(p[0]);  T tmp2(p[1]); p[0] = p[2];
		   p[2] = tmp1;   p[1] = p[3];  p[3] = tmp2; }
int compare(const Quadruple<T>& p2) 
                 {  
		   for(int i=0; i< 4; i++)
		     if (p[i] != p2(i))
		       return (p[i] < p2(i)) ? -1 : 1;
		   return 0;
		 }
T&       operator[](int i)         { assert(i >=0 && i < 4);  return p[i]; } 
const T& operator()(int i) const   { assert(i >=0 && i < 4);  return p[i]; } 
int      operator <(const Quadruple<T>& qq)  { return (compare(qq) == -1); } 
int      operator!=(const Quadruple<T>& qq)  { return (compare(qq) !=  0); } 
int      operator==(const Quadruple<T>& qq)  { return (compare(qq) ==  0); } 
void     permute  (int permflag)
         {
	   if (permflag & TRANS12  ) perm12(); 
	   if (permflag & TRANS34  ) perm34(); 
	   if (permflag & TRANSPAIR) permpair(); 
	 }
};

template <class T> ostream& operator<<(ostream& os,const Quadruple<T> q)
{
  for(int i=0; i<4; i++)
    os << q(i);
  return os;
}

#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
template <class T> ostream_withassign& operator<<(ostream_withassign& os,
						  const Quadruple<T> q)
{
  for(int i=0; i<4; i++)
    os << q(i);
  return os;
}
#endif

#endif
