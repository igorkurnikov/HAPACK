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
% $Id: symlogic.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
\section{Logic for Target Configurations}
  
We are given 4 basis sets $b_i$ with symmetry characteristics $s_i$
and the standard interchange flags $s_{12}$,$s_{34}$ and $s_p$. The
goal of this section is to compute all non-equivalent symmetry
combinations of integrals for the given set of basis-sets. 

Each basis set is characterized by maximally $2^{n(s_i)}$  
symmetry desigators with respect to the symmetry operations $P_x,P_y
and P_z$. Functions from a basis set with coordinates  (1,-1,0), for example
(assuming fill, X,Y and Z) symmetry, will require the construction of
functions which are $(E_xE_y),(E_xO_y),(O_xE_y),(O_xO_y)$. Here $E_x$
indicates an even function with respect to the $x$-reflection, $O_x$
odd functions with respect to $x$ reflections and so on. Note that no
linear combination of functions is required to obtain $z$ symmetry,
since the basis transforms onto itself under $z$ interchange. 

Having established the required symmetry characteristics for the
individual positions we form all combinations of such characteristics
on the 4 indices. Not all such combinations, however amount to unique
sets of integrals. If $b_3 = b_4$ e.g., the sets $[EE|EO]$ and
$[EE|OE]$ are identical. We apply the possible permutations 
indicated by the $s_{12}$,$s_{34}$ and $s_p$ to the allowed
configurations and elmiminate all which appear twice.
*/

#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include <errorip.h>
#include "const.h"
#include "basis.h"
#include "storage.h"
#include "twoel.h"
#include "symdesig.h"
#include "tquadruple.h"
#include "symop.h"

void TwoElectronIntegrals::generate_targets(int l1,int l2,int l3,int l4)
{
  target.reset();
  
  for(SymDesignator s1(b1.sym_equiv(X),b1.sym_equiv(Y),b1.sym_equiv(Z));
      s1.in_range(); ++s1)
  for(SymDesignator s2(b2.sym_equiv(X),b2.sym_equiv(Y),b2.sym_equiv(Z));
      s2.in_range(); ++s2)
  for(SymDesignator s3(b3.sym_equiv(X),b3.sym_equiv(Y),b3.sym_equiv(Z));
      s3.in_range(); ++s3)
  for(SymDesignator s4(b4.sym_equiv(X),b4.sym_equiv(Y),b4.sym_equiv(Z));
      s4.in_range(); ++s4)
  {
    Quadruple<SymDesignator> q(s1,s2,s3,s4);
    // cout << " +++ +++ Considering Target: " << q << endl;
    
    Quadruple<SymDesignator> q12(q);
    Quadruple<SymDesignator> q34(q);
    Quadruple<SymDesignator> q1234(q);

    Quadruple<SymDesignator> pq12(q);
    Quadruple<SymDesignator> pq34(q);
    Quadruple<SymDesignator> pq1234(q);    
    Quadruple<SymDesignator> p(q);

    q12.perm12();     q34.perm34();  q1234.perm12();  q1234.perm34();
    pq12.perm12();   pq34.perm34(); pq1234.perm12(); pq1234.perm34();
    pq12.permpair(); pq34.permpair(); pq1234.permpair(); p.permpair();
    
    if (same12 && l1 == l2 && q12 < q) continue;
    if (same34 && l3 == l4 && q34 < q) continue;
    if (same12 && same34 && l1 == l2 && l3 == l4 && q1234 < q) continue;
    
    int sp_flag = same_pair && l1 == l3 && l2 == l4;
    
    if (sp_flag && p < q) continue;
    if (sp_flag && l1 == l2 && same12 && pq12 < q) continue;
    if (sp_flag && l3 == l4 && same34 && pq34 < q) continue;
    if (sp_flag && l1 == l2 && same12 && l3 == l4 && same34 && pq1234 < q) 
      continue;

    // cout << " +++ +++ +++ Target: " << q << endl;    
    target.add(q);
  }
}

/*TEX
\subsection{Determining Inequivalent Primitive Configurations}

We are given a quadruplet of basis sets with interchange properties
specified by $s_{12},s_{34}$ and $s_p$ and with symmetry
characteristics in the X,Y and Z directions. 

Any primitive configuration is thus charaterized by a quadruplet of
'legal' \name{SymOp's} and our first goal is to determine the minimal
set of such quadruplets to compute the symmetry adapted integrals. To
this end we introduce the canonical ordering of such quadruplets. 
Two primitive configurations produce the same integrals iff one of the
primitive configurations can be transformed into another by allowed
symmetry operations on the quadruplet. Such symmetry operations are
the permutations as indicated by the flags $s_{12},s_{34}$ and $s_p$
and the application of any reflection operator on ALL indices of the
quadruplet. We specify that of each pair of `equivalent' primitive
configurations only the smaller one in the canonical ordering will be 
be generated.

*/

void reflect(Quadruple<SymOp>& p,Direction dir)
{
  for(int i=0; i< 4; i++)
    p[i].reflect(dir);
}

void apply(Quadruple<SymOp>& p,SymOp& op)
{
  for(int i=0; i< 4; i++)
    p[i] = p[i] ^ op;
}

void TwoElectronIntegrals::generate_primitives()
{
  primitive.reset();
  
  for(SymOp s1(b1.sym_equiv(X),b1.sym_equiv(Y),b1.sym_equiv(Z));
      s1.in_range(); ++s1)
  for(SymOp s2(b2.sym_equiv(X),b2.sym_equiv(Y),b2.sym_equiv(Z));
      s2.in_range(); ++s2)
  for(SymOp s3(b3.sym_equiv(X),b3.sym_equiv(Y),b3.sym_equiv(Z));
      s3.in_range(); ++s3)
  for(SymOp s4(b4.sym_equiv(X),b4.sym_equiv(Y),b4.sym_equiv(Z));
      s4.in_range(); ++s4)
  {
    Quadruple<SymOp> q  (s1,s2,s3,s4);
    int prim = true;
    
    for(SymOp oper((b1.sym_equiv(X) || b2.sym_equiv(X) || b3.sym_equiv(X) || b4.sym_equiv(X)),
		   (b1.sym_equiv(Y) || b2.sym_equiv(Y) || b3.sym_equiv(Y) || b4.sym_equiv(Y)),
		   (b1.sym_equiv(Z) || b2.sym_equiv(Z) || b3.sym_equiv(Z) || b4.sym_equiv(Z)));
		   oper.in_range();++oper)
    // SymOp oper(0,0,0);
    {
      Quadruple<SymOp> p(s1,s2,s3,s4);
      apply(p,oper);
      // cout << " Testing Primitive: " << p << " " << oper << endl;  
      if (p < q) 
      {
	prim = false;
	break;
      }
    }
    if (prim)
    {
      // cout << q << " is primitive. " << endl;
      primitive.add(q);
    }
  }  
}
/*TEX
\subsection{Generate Sources}

Given a particular \name{Quadruple} of \name{SymOp}, here we determine all
such permutations and reflections which lead to new symmetry configurations,
such that the integrtals differ at most by a sign from those computed
in the primitive.  
*/

int sort_in(SStack<Quadruple<SymOp>,max_sym >& stack, Quadruple<SymOp>& p)
{
  for(int i=0; i < stack.size(); i++)
    if (stack[i] == p) return false;
  stack.add(p);
  return true;
}


void TwoElectronIntegrals::generate_sources(Quadruple<SymOp>& q)
{
  sources.reset();
  operation.reset();

  for(SymOp oper((b1.sym_equiv(X) || b2.sym_equiv(X) || b3.sym_equiv(X) || b4.sym_equiv(X)),
		   (b1.sym_equiv(Y) || b2.sym_equiv(Y) || b3.sym_equiv(Y) || b4.sym_equiv(Y)),
		   (b1.sym_equiv(Z) || b2.sym_equiv(Z) || b3.sym_equiv(Z) || b4.sym_equiv(Z)));
      oper.in_range();++oper)
  {
    Quadruple<SymOp> p(q);
    apply(p,oper);
    // cout << " Testing Source: " << p << " " << oper << endl;  
    if (sort_in(sources,p))
      operation.add(oper);
  }  

//  cout << " +++ Sources: " << endl;
//  for(int i = 0;i < sources.size(); i++)
//    cout << "     " << sources[i] << " oper: " << operation[i] << endl;
//  cout << endl;
}


/*TEX
\subsection{Compute Factor to contract a Primitive into a Target}

Given a target \name{SymDesignator} and a \name{SymOp} characterizing 
a primitive configuration we must determine the prefactor with which the 
primitive configuration enters to the target configuration. For each direction
we gather:
\begin{itemize}
\item A factor of $1/\sqrt{2}$ for each index decorated with an $E$ or an $O$.
\item A factor of $-1$ for each index, where the source descriptor indicates 
      a reflection {\bf and} the target descriptor indicates $O$.
\item A factor of $-1$ for each index, where the a reflection was required to 
      go from the original source to the present configuration
\end{itemize}
 
*/

double TwoElectronIntegrals::compute_factor(Quadruple<SymDesignator>&  tar,
					    Quadruple<Momentum>     &  m,
					    Quadruple<SymOp>        &  src,
					    SymOp                   & oper)
{
  Quadruple<SymOp> orig(src);
  apply(orig,oper);
  assert(orig == sources[0]);
  
  double root = 1.0/sqrt(2.0);
  double fac  = 1.0;

  Direction dir = X;
  int i;
  for(i = 0 ; i < 4; i++)
  {
    if (tar[i](dir) != none)                                 fac *= root;
    if (tar[i](dir)  == odd && src[i](dir) == reflection)    fac *= -1;
    if (orig[i][dir] == reflection && (m[i](dir) % 2 != 0))  fac *= -1;    
  }

  dir = Y;  
  for(i = 0 ; i < 4; i++)
  {
    if (tar[i](dir) != none) fac *= root;
    if (tar[i](dir) == odd && src[i](dir) == reflection)    fac *= -1;
    if (orig[i][dir] == reflection && (m[i](dir) % 2 != 0)) fac *= -1;
  }

  dir = Z;  
  for(i = 0 ; i < 4; i++)
  {
    if (tar[i](dir) != none) fac *= root;
    if (tar[i](dir) == odd && src[i](dir) == reflection)     fac *= -1;
    if (orig[i][dir] == reflection  && (m[i](dir) % 2 != 0)) fac *= -1;
  }

  return fac;
} 
