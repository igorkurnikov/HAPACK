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
% $Id: function.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%

\section{RadialFunction}

This class implements the representation of gaussian primitives as
they are used in the recursion relations. Each function consists of a
center, see class \name{Location} and and exponent. Furthermore we
store the indices of all functions into which this particular primitive
\name{RadialFunction} contributes in contractions as well as the
contraction coefficients. The \name{exp()} member returns the exponent
and \name{contractions} the number of contracted functions to which
this primitive function contributes. The indices of the individual
contracted functions and the contraction coefficients are returned by
the \name{target} and \name{coef}  member functions. The \name{norm(l)}
member returns the norm of the primitive function, defined as:
\be
{\cal N } = \left(\frac{2 \zeta}{\pi} \right)^{3/4}  (4 \zeta)^{l/2}  
\ee
where $l$ is the angular momentum. Note that the function specific
prefactor:
$$\left[ (2 l_x - 1)!!\ (2 l_y - 1)!!\ (2 l_z - 1)!! \right]{1/2}$$
is missing, it is automatically picked up in the auxiliary
normalization, see section ~\ref{section-normalize}. 

We furthermore provide a comparison operator, which compares exponent
and location of two radialfunctions. Since real number are compared,
we use the criterion that two number are equal, if they differ by less 
than $\delta$, which is defined in section~\ref{section-const}.
*/
#ifndef FUNCTION_DEF
#define FUNCTION_DEF

#include "io_incl.h"
#include <math.h>
#include <tsstack.h>
#include "typesip.h"
#include "const.h"


class RadialFunction : public Location
{
  friend class RadialFunctionArray;
  double  ee;
  SStack<   int,max_contr>  tr;
  SStack<double,max_contr>  cf;
public:
          RadialFunction()  : Location(0,0,0) { ee = -1; }
          RadialFunction(const Location& l,const double e,
			 const int tar,const double cc) :
                   Location(l),ee(e) { tr.add(tar); cf.add(cc); }
         ~RadialFunction() {}
RadialFunction& operator=(const RadialFunction& g)
          { ee = g.ee; tr = g.tr; cf = g.cf; assign(g); return *this; }
double    exp()    const { return ee; }
int       contractions()      const { return tr.size(); } 
int       target(const int i) const { return tr(i); }
double    coef(  const int i) const { return cf(i); }
double    norm(int momentum)  const
          {
	    double x = 2*ee/pi;
	    return (momentum == 0) ?
	      sqrt(sqrt(x*x*x)) :
	      sqrt(sqrt(x*x*x)*pow(4*ee,momentum));
	  }
friend    istream& operator>>(istream& is,RadialFunction& rs);
};

ostream& operator<<(ostream& os,const RadialFunction& rs);

inline int operator==(const RadialFunction& a,const RadialFunction& b)
{
  const Location& ac = a;
  const Location& bc = b;
  
  return (ac == bc && fabs(a.exp()-b.exp()) < delta);
}
/*TEX
\section{RadialFunctionArray} 

A \nameindex{RadialFunctionArray}  contains all
\name{Radialfunction}s of a given angular momentum and is implemented as a
stack of \name{Radialfunction}s.  See template \name{SStack}.
Therefore \name{RadialFunctionArray}s need not be explicitly sized,
extra functions may be added through \name{add}.

Note that add operates somewhat non-trivial: It first checks, whether
there is already a \name{RadialFunction} with the same center present,
in which case the contraction coefficients and targets are added to this
function, otherwise a new \name{Radialfunction} is created and added
to the list.
*/
class RadialFunctionArray : public SStack<RadialFunction,max_radials> 
{
public:
  RadialFunctionArray() {} 
 ~RadialFunctionArray() {}
void print(ostream& os);  
int  nofun()  const;  // return the maximal index of contracted functions
int  add(const RadialFunction& rf,int print = 0);
};

ostream& operator<<(ostream& os,RadialFunctionArray& ss);

#endif
