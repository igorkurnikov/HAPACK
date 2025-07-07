/*
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
% $Id: one_el.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
\section{Symmetry Adapted Single Electron Integrals}

\subsection{Concepts}
In order to explain the ideas of this chapter, we will use $H_2O$ as an 
example.
Suppose the $H_20$ molecule to be positioned in space
symmetically to the $YZ$-plane and we are given one basis set for the
$H$ atoms located in the positive halfsphere
and another one for the $O$ atom. The goal of this chapter is to create out of 
these informations the symmetry adapted basis functions respectively 
all the information
that is needed to calculate symmetry adapted single electron integrals.
In our example we have two hydrogen atoms, whose basis functions can be 
combined to symmetry adapted basis functions with even or odd symmetry with
respect to reflection on the $YZ$-plane.
We label these functions with the \name{SymDesignators} "ENN" and "ONN". Each
letter corresponds to one of the three coordinate directions. 
"E" stands for even, "O" for odd and "N" for no (known) symmetry. If for an
example the atomic basis for the hydrogen consists only of a $1s$-orbital
we will have to calculate the following symmetry adapted one electron 
integrals:
\be 
\begin{array}{cc}
 \langle ENN~(1s)~|O|~ENN~(1s)\rangle &\langle ENN~(1s)~|O|~ONN~(1s)\rangle \\
 \langle ONN~(1s)~|O|~ENN~(1s)\rangle &\langle ONN~(1s)~|O|~ONN~(1s)\rangle \\
\end{array}
\ee
(Of course the combination of even and odd is zero in this example.)
We call these pairs of \name{SymDesignators} \name{targets} and store
the values of the related integrals in an array called \name{tconf}.

The symmetry adapted \name{target} base functions (respectively
the related one electron integrals) must be constructed via linear
combination of the primitive basis sets for each of the $H$ atoms, but the 
basis set for the $H$ atom in the negative halfsphere does not exist.
Nevertheless it can be represented by a symmetry operation
(\name{SymOp}) acting on the basis set in the positive halfsphere.
In our example it is just the reflection on the $YZ$ plane, which is
expressed by the \name{SymOp} ``R11'', where the ``R'' stands for
reflection in the first and the ``1'' for identity in the second and third
coordinate. 

We therefore enable reflection symmetry for the $X$ direction of the
$H$ basis set and allways build both $H$ atoms from it.
If for example both basis sets passed to the \name{OneElectronIntegrals} class 
are the $H$ basis set, the following pairs of primitives will be constructed:
\be 
\begin{array}{cc}
  [H(111)|H(111)] & [H(R11)|H(111)]\\
  [H(111)|H(R11)] & [H(R11)|H(111)]\\
\end{array}
\ee
We call these pairs of non-symmetry adapted primitive basis functions
\name{sources}. As you can see, not all of the \name{sources} have to
be created, because some of them like
\be
  [H(111)|H(R11)] = SymOp(R11)(~[H(R11)|H(111)]~).
\ee
are equivalent, that is they can be made identically, when acting with
the same allowed \name{SymOp} on both elements of the pair. This
means, that we only have to evaluate the recursions relations for
one electron integrals with a minimal set of inequivalent
\name{sources}, which we label as \name{primitives}.

\subsection{Strategy of implementation\label{OneElec-Strat}}
The resulting strategy is:
\begin{tabbing}
~~~\=~~~\=~~~\=~~~\=~~~ \kill \\
\>Make a list of \name{target} configurations. \\
\>For all pairs of total angular momenta \name{evaluate} integrals, that is \\
\> \>make a list of (inequivalent) \name{primitives}. \\
\> \>For all elements of \name{primitive} \\
\> \> \>\name{generate} non-symmetry adapted integrals via recursion and\\
\> \> \>\name{distribute} these integrals into the target-array \name{tconf},\\
\> \> \> \>thereby taking into account prefactors for the different \\
\> \> \> \>\name{source} configurations. 
\> \>Finally we \name{output} the \name{tconf} field into the output
     matrices (\name{mat_array}).

In order to start specific recursion relations for the different types
of one electron operators, the above strategy is implemented in an
abstract base class, from which we derive subclasses with specific
\name{generate} function.

Thereby the reflection symmetry of the one electron operators must be
recognized, for calculating the correct prefactors for the \name{source}
configurations. This is done by the virtual function \name{operator_fac}.

  

\subsection{Header Of Abstract Base Class}
\index{OneElectronIntegrals, class header}
*/
class Coefficient_Set;
class Matrix;

#include "tpair.h"
#include "symop.h"
#include "symdesig.h"
#include "rec_stor.h"
#include "typesip.h"
#include <tarray.h>

enum RefDir { RX = 1, RY = 2, RZ = 4 }; 

class OneElectronIntegrals
{
 protected:
  int                            debug;
  
  InternalBasis&                 b1,b2;
  int                            same;           // indicates if b1 == b2

  int                            min_mu;
  int                            max_mu;

  ARRAY <IntegMat>&                   mat_array;
  SStack<Pair<SymOp>,max_sym>    primitive;
  SStack<Pair<SymDesignator>,max_sym>   target;

  SStack<Pair<SymOp>,max_sym>    sources;
  SStack<SymOp,max_sym>          operation;      // used to create source

  ARRAY<RecursionStorage*>       tconf;
  
public:
  OneElectronIntegrals(InternalBasis& bb1,InternalBasis& bb2,
		       int s12,ARRAY<IntegMat>& m, int mx_mu = 0, int mn_mu = 0);
 ~OneElectronIntegrals();
  
void              logic();  
void              evaluate(int l1,int l2);
virtual void      output(int l1,int l2);

virtual RecursionStorage* generate(int l1,SymOp& s1,int l2,SymOp& s2) = 0;
  
  
void              generate_targets();
void              generate_sources(Pair<SymOp>& s);
void              generate_primitives();
double            compute_factor(Pair<SymDesignator>& tar,Pair<Momentum>& m,
				 Pair<SymOp>& src,SymOp& oper,Momentum mu);
virtual int       operator_fac(Direction dir, Momentum mu) = 0;
void              distribute(Storage* ss,Pair<Momentum>& m,Pair<SymOp>& source,
			     SymOp& oper,Momentum mu);
  
};

/*TEX
\section{Headers Of Operator Dependend Subclasses}
For each one electron operator a specific recursion relation must be
evaluated. This is achieved by the following subclasses through the 
individual \name{generate} function. The
class headers are summerize here, but the \name{generate} functions
are listed after the implementation of the base class in section 
\ref{Eval-Rec}. All constructors initialize local variables required for
the recursion and start the calculation by calling the base class function
\name{logic}.

\index{OverlapOp, class header}
*/
class OverlapOp : public OneElectronIntegrals
{
public:
  OverlapOp(InternalBasis& bb1,InternalBasis& bb2,int s12,ARRAY<IntegMat>& m)
    : OneElectronIntegrals(bb1,bb2,s12,m) 
  { logic(); } 
  ~OverlapOp() {}
  virtual RecursionStorage* generate(int l1,SymOp& s1,int l2,SymOp& s2);
  virtual int operator_fac(Direction dir, Momentum mu){return 1;}  
};

/*TEX
\index{KineticOp, class header}
*/
class KineticOp : public OneElectronIntegrals
{
public:
  KineticOp(InternalBasis& bb1,InternalBasis& bb2,int s12,ARRAY<IntegMat>& m)
    : OneElectronIntegrals(bb1,bb2,s12,m) 
  { logic(); } 
  ~KineticOp() {}    
  virtual RecursionStorage* generate(int l1,SymOp& s1,int l2,SymOp& s2);
  virtual int operator_fac(Direction dir, Momentum mu){return 1;}
};
/*TEX
\index{MomentOp, class header}
*/
class MomentOp : public OneElectronIntegrals
{
  Location loc;

public:
  MomentOp(InternalBasis& bb1,InternalBasis& bb2,int s12,ARRAY<IntegMat>& m,
           MomentData& data)
    : OneElectronIntegrals(bb1,bb2,s12,m,data.max_mu) 
  { 
    loc = data.loc;
    logic(); 
  } 
 ~MomentOp() {}
 virtual RecursionStorage* generate(int l1,SymOp& s1,int l2,SymOp& s2);
 virtual int operator_fac(Direction dir, Momentum mu)
 { if( (mu(dir) %2 != 0)  && SymOp::symmetry(dir) )
     return -1;
   else   
     return 1;
 }
};


/*TEX
\index{NuclearOp, class header}
*/
class NuclearOp : public OneElectronIntegrals
{
  Location loc;
  double   chg;
  
public:
  NuclearOp(InternalBasis& bb1,InternalBasis& bb2,int s12,ARRAY<IntegMat>& m,
	    Location& l,double cc)
    : OneElectronIntegrals(bb1,bb2,s12,m) 
  { 
    loc = l;
    chg = cc;
    logic(); 
  } 
 ~NuclearOp() {}      
 virtual RecursionStorage* generate(int l1,SymOp& s1,int l2,SymOp& s2);
 virtual int operator_fac(Direction dir, Momentum mu){return 1;}
};

/*TEX
\index{E_FieldOp, class header}
*/
class E_FieldOp : public OneElectronIntegrals
{
  Location loc;
  
public:
  E_FieldOp(InternalBasis& bb1,InternalBasis& bb2,int s12,ARRAY<IntegMat>& m,
	    Location& l)
    : OneElectronIntegrals(bb1,bb2,s12,m,2) 
  { 
    loc = l;
    logic(); 
  } 
 ~E_FieldOp() {}     
 virtual RecursionStorage* generate(int l1,SymOp& s1,int l2,SymOp& s2);
 virtual int operator_fac(Direction dir, Momentum mu)
 { 
   if( (mu(dir) == 1) && SymOp::symmetry(dir) ) 
     return -1;
   else
     return 1;
 }
} ;

/*TEX
\index{AngularOp, class header}
*/
class AngularOp : public OneElectronIntegrals
{
  Location loc;
  
public:
  AngularOp(InternalBasis& bb1,InternalBasis& bb2,int s12,ARRAY<IntegMat>& m,
	    Location& l)
    : OneElectronIntegrals(bb1,bb2,s12,m,1,1) 
  { 
    loc = l;
    logic(); 
  } 
 ~AngularOp() {}  
 virtual void              output(int l1, int l2);
 virtual RecursionStorage* generate(int l1,SymOp& s1,int l2,SymOp& s2);
 virtual int operator_fac(Direction dir, Momentum mu){ return 1; }
};

/*TEX
\index{SpinOrbOp, class header}
*/
class SpinOrbOp : public OneElectronIntegrals
{
  Location loc;
  
public:
  SpinOrbOp(InternalBasis& bb1,InternalBasis& bb2,int s12,ARRAY<IntegMat>& m,
	    Location& l)
    : OneElectronIntegrals(bb1,bb2,s12,m,1,1) 
  { 
    loc = l;
    logic(); 
  } 
 ~SpinOrbOp() {}     
 virtual void              output(int l1, int l2);
 virtual RecursionStorage* generate(int l1,SymOp& s1,int l2,SymOp& s2);
 virtual int operator_fac(Direction dir, Momentum mu){ return 1; }
};


  







