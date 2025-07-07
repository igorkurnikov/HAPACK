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
% $Id: function.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\subsection{Implementation}
*/
#include "io_incl.h"
#include "function.h"
 
ostream& operator<<(ostream& os,const RadialFunction& g)
{ 
  char buffer[100];
  sprintf(buffer,"RFUN: EXP = %10.3f at ",g.exp());
  os << buffer << "*";
  const Location& gg = g;
  os << gg;
  os << " --- #Cont. : " << g.contractions() << endl;
  for(int i = 0; i< g.contractions(); i++)
  {
    os << "        +++ contracted to function no:: " << g.target(i);
    sprintf(buffer," with coefficient %10.3f.",g.coef(i));
    os << buffer << endl;
  }
  return os;
}

istream& operator>>(istream& is,RadialFunction& g)
{ 
  is >> g.ee >> g.pos[0] >> g.pos[1] >> g.pos[2];
  return is;
}

void RadialFunctionArray::print(ostream& os)
{
  SStack<RadialFunction,max_radials>& tmp = *this;
  os << tmp;
}
/*TEX
\subsection{Number of Functions}
This function counts the number of functions in the \name{RadialFunctionArray}.
*/
int RadialFunctionArray::nofun() const
{
  int mxind = 0;
  for(int i =0; i< size(); i++)
  {
    for(int j=0; j < base[i].contractions(); j++) 
      if (base[i].target(j) > mxind)
	mxind = base[i].target(j);
  }
  return mxind+1;
}

ostream& operator<<(ostream& os,RadialFunctionArray& ss)
{
  cout << "Radial functions: " << endl;
  
  for(int i=0; i< ss.size();i++)
    os << ss[i] << endl;
  return os;
  
}
/*TEX
\subsection{add}
When adding a function to the \name{RadialFunctionArray}, we first
check whether a function with the same angular momentum, center and 
exponent already exisits. If so, we merge the new function with the
already exisiting one, i.e. we only add the contraction coefficient
and target function index to the existing array. Otherwise a new
function is being enteres into the array.
*/
int  RadialFunctionArray::add(const RadialFunction& rf,int print)
{
  if (print)
    cout << " ++ adding: " << rf << endl;
  for(int i=0; i< size(); i++)
  {
    if (base[i] == rf)
    {
      if (print)
	cout << " ++++ Merged with function: " << i << endl;
      base[i].tr.add((SStack<   int,max_contr>&) rf.tr);
      base[i].cf.add((SStack<double,max_contr>&) rf.cf);
      return i;
    }
  }  
  int nn = SStack<RadialFunction,max_radials>::add(rf);
  if (print)
    cout << " ++++ Function Number is:  " << nn << endl;
  return nn;
}




