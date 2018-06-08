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
% $Id: one_el.cpp,v 1.2 2008/05/03 07:00:37 igor Exp $
%

\section{Implementation of OneElectronIntegrals}
*/
#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <math.h>
#include "basis.h"
#include "function.h"
#include "integ_array.h"
#include "coef_set.h"
#include <numer.h>
#if !defined(_MSC_VER)
#include <sys/times.h>
#include <sys/resource.h>
#endif
#include "symop.h"
#include "one_el.h"
#include "integ_file.h"
#include "overlap.h"
#include "kinetic.h"
#include "moment.h"
#include "nuclear.h"
#include "e_field.h"
#include "angular.h"
#include "spin_orb.h"

/*TEX
\subsection{Constructor}
\index{OneElectronIntegrals, constructor}
The constructor assigns the arguments to the local variables.
Because \Name{OneElectronIntegrals} is an abstract class the
initiation of the calculation is done by calling \name{logic} in the
constructors of the derived classes.
*/
OneElectronIntegrals::OneElectronIntegrals(
	     InternalBasis& bb1,InternalBasis& bb2,int ss12,
             ARRAY<IntegMat>& m, int mx_mu,int mn_mu) :
     b1(bb1),b2(bb2),same(ss12),mat_array(m)
{ 
  debug  = false; 
  max_mu = mx_mu;
  min_mu = mn_mu;
}

/*TEX
\subsection{Destructor}
\index{OneElectronIntegrals, destructor}
Just throw away the remaining target configurations. 
*/  
OneElectronIntegrals::~OneElectronIntegrals() 
{
  for(int i =0; i < tconf.size(); i++)
    if (tconf[i] != 0)
      delete tconf[i];
}

/*TEX
\section{Initiation of Calculation}
\index{OneElectronIntegrals, logic}  
These function is called by the constructors of the derived subclasses
for the specific one electron operators and starts the calculation of
the symmetry adapted one electron integrals in the above described way
(see section \ref{OneElec-Strat}).
*/

void OneElectronIntegrals::logic()
{
  generate_targets();
  
  int l1,l2;
  
  for(l1=0; l1 < b1.maxl; l1++)
  {
	  if (b1.noprim[l1] > 0)
	  {
		  int maxl2 = l1 + 1;
		  if (!same)
			  maxl2 = b2.max_momentum();
		  for(l2=0; l2 < maxl2; l2++)   // run over all pairs of angular momenta
		  {
			  if (b2.noprim[l2] > 0)
			  {
				  evaluate(l1,l2);
				  output(l1,l2);			  
			  }
		  }
	  }
  }
}


/*TEX
\section{Logic for Target Configurations}
We are given 2 basis sets $b_i$ possibly with reflection symmetry in each 
direction and the \name{same} flag, which indicates whether the basis 
sets are identically or not. The goal of this section
is to compute all non-equivalent symmetry combinations of integrals
for the given two basis-sets.

Each basis set is characterized by maximally $2^3 = 8$ symmetry
designators with respect to the reflection operations $P_x,P_y and
P_z$. Functions from a basis set with coordinates (1,-1,0), for
example (assuming X,Y and Z symmetry), will require the construction of
functions which are $(EEN),(EON),(OEN),(OON)$. 
Here $EON$ indicates an even function with respect to the $x$-reflection and 
odd functions with respect to $y$ reflections and so on. Note that no
linear combination of functions is required to obtain $z$ symmetry,
since the basis transforms onto itself under $z$ interchange.
This indicated by $N$. Which linear combinations are usefull is proved
by \name{sym_equiv}.
  
Having established the possible symmetry characteristics for the
individual direction we form all combinations of such characteristics
on the two indices. Not all such combinations, however amount to
unique sets of integrals. If the two basis sets are the same, 
for example $[EEN|OEN]$ is the same as $[OEN|EEN]$ 
(never mind that both are zero in this case). 
Using the compare operation of the \name{Pair} class such identities are
checked. At least the unique pairs are added to the \name{target} stack and
finally the size of the \name{tconf} array is resetted according to the
number of targets. 

\index{OneElectronIntegrals, generate_targets}  
\index{target}
\index{tconf}
*/
  
void OneElectronIntegrals::generate_targets()
{
  for(SymDesignator s1(b1.sym_equiv(X),b1.sym_equiv(Y),b1.sym_equiv(Z));
      s1.in_range(); ++s1)
  for(SymDesignator s2(b2.sym_equiv(X),b2.sym_equiv(Y),b2.sym_equiv(Z));
      s2.in_range(); ++s2)
  {
    Pair<SymDesignator> q(s1,s2);
    
    if (debug)
      cout << " +++ Checking Target: " << q << endl;
    
    Pair<SymDesignator> q12(q);

    q12.permute();
    if (same && q12 < q) continue;
    
    if (debug) 
      cout << " +++ Target: " << q << endl;    
    target.add(q);
  }

  tconf.reset(target.size());
  tconf.set(0);
  
}

/*TEX
\subsection{Evaluation of Primitives}
  
For our $H_20$ example we are given two basis sets, one for oxigen
located in the $YZ$-plane and one for hydrogen located say in the
positive $X$ halfsphere. Lets call them real basis sets. 
We have to build the symmetry adapted basis function as linear combinations of
given primitive basis functions and therefore we have to 
to create by reflection the primitive virtual basis set for 
hydrogen in the negative $X$ halfsphere. The possible primitive
configurations can be represented by \name{SymOp's} as shown in the
following example.\\
Consider the case of the [ENN|ONN] \name{target}:
\begin{eqnarray}
  [ENN(1s)|ONN(1s)] & = & [\frac{1}{\sqrt{2}}\{(111)(1s) + (R11)(1s)\} | 
                         \frac{1}{\sqrt{2}}\{(111)(1s) - (R11)(1s)\} ] 
                         \nonumber \\
              & = & \frac{1}{2}\left\{[(111)(1s)|(111)(1s)] 
                        - [(111)(1s)|(R11)(1s)] + [(R11)(1s)|(111)(1s)] 
                        - [(R11)(1s)|(R11)(1s)]\right\} \nonumber \\ 
              & = & \frac{1}{2}\left\{[(111)(1s)|(111)(1s)] 
                        - [(111)(1s)|(R11)(1s)] + [(111)(1s)|(R11)(1s)] 
                        - [(111)(1s)|(111)(1s)]\right\} \nonumber \\ 
              & = & 0 \label{source-signs}\\
\end{eqnarray}              
The \name{sources}, that is all possible pairwise combinations of
\name{SymOp's} are not unique. They can be made equivalent through the
application of allowed symmetry operations to the whole pair, as shown
in the last step.
In order to compute the integrals for a choosen pair of angular momenta,
we therefore first construct only these minimal set of
\name{primitives} and start the recursive calculation for them.
For each \name{primitive} we next construct all equivalent
\name{sources}. They are neccessary for getting the correct signs for
each term. At least we \name{distribute} the sources
into the target configurations (keeping an eye on signs and
normalization factors).

This means, that we do not check whether all terms might cancel like
in the above example, before we start calculating the integrals recursivlely.
\index{OneElectronIntegrals, evaluate}
*/

void OneElectronIntegrals::evaluate(int l1,int l2)
{
  /* Configuration of tconf */
  /* only the integrals with aux = 0 are finally relevant */

  if(debug)  
  {
    cout << "EVALUATE \n";
    cout << "tconf.size " << tconf.size() << endl;
  }
  
  for(int i=0; i< tconf.size(); i++)
  {
    if(debug)
      cout << "tconf[" << i << "]\n ";
    RecursionStorage *rec = new RecursionStorage(l1,l2,0,max_mu,
					b1.radial_functions(l1).nofun(),
					b2.radial_functions(l2).nofun(),1,1,1);
    rec -> init_fields();
    assert(rec != 0);
    if (tconf[i] != 0) 
      delete tconf[i];
    tconf[i] = rec;
  } 

  generate_primitives();

  if(debug)  
    cout << "For " << primitive.size() << "primitives do \n";
  
  for(int p=0; p < primitive.size(); ++p)
  {
    /* allocate memory for recursion and calculate recursion by
       calling the specific inherit function of generate */
    RecursionStorage *rec = generate(l1,primitive[p][0],l2,primitive[p][1]);
    assert(rec != 0);

    generate_sources(primitive[p]);

    for(int mutot = min_mu; mutot <= max_mu; mutot++)    
    for(Momentum mu(mutot);mu();++mu)    //for all mu with mutot <= max_mu
    for(Momentum m1(l1); m1();++m1)
    for(Momentum m2(l2); m2();++m2)
    {
      if(debug)
      {
        cout << "p : " << p <<" mutot : " << mutot << " mu : " << mu << endl;
        cout.flush();
      }
      
      Storage* ss = rec -> find(m1,m2,0,mu);
      // cout << "Array: " << *ss << endl;
      
      assert(ss != 0);
      for(int src = 0; src < sources.size(); src++)
      {
	Pair<Momentum> m(m1,m2);
	distribute(ss,m,sources[src],operation[src],mu);	
      }
    }      
    delete rec;
  }
}

/*TEX
\subsection{Determining Inequivalent Primitive Configurations}
\index{OneElectronIntegrals, generate_primitives}

What we are doing here is in a way creating the 'virtual' basis sets
by reflection
of the real basis sets. Therefore we form all allowed symmetry 
combinations of the two basis sets that are required for the symmetry 
adapted targets (see generate_targets for more comment).
Any primitive configuration can thus be charaterized by a pair of
\name{SymOp's} and our first goal is to determine the {\bf minimal}
set of such pairs to compute the symmetry adapted integrals. 

To this end we introduce the canonical ordering of pairs via the
ordering of the \name{SymOp's}. That means the first symmetry
operation of the pair should be 'less' than the second.
 
Two primitive configurations produce the same integrals if one of the
primitive configurations can be transformed into another by allowed
symmetry operations on the pair. Such symmetry operations are
the application of any reflection operator on ALL indices of the
pair. We specify that of each pair of `equivalent' primitive
configurations only the smaller one in the canonical ordering will be 
be generated.
\index{reflect}
*/
void reflect(Pair<SymOp>& p,Direction dir)
{
  for(int i=0; i< 2; i++)
    p[i].reflect(dir);
}

/*TEX
\index{apply}
*/
void apply(Pair<SymOp>& p,SymOp& op)
{
  for(int i=0; i< 2; i++)
    p[i] = p[i] ^ op;
}

void OneElectronIntegrals::generate_primitives()
{
  primitive.reset();
  
  for(SymOp s1(b1.sym_equiv(X),b1.sym_equiv(Y),b1.sym_equiv(Z));
      s1.in_range(); ++s1)
  for(SymOp s2(b2.sym_equiv(X),b2.sym_equiv(Y),b2.sym_equiv(Z));
      s2.in_range(); ++s2)
  {
    Pair<SymOp> q  (s1,s2);
    int prim = true;
    if (debug)
      cout << " Testing Primitive: " << q << " "  << endl;  

    // run over all all operators which perform a coordinate change on BOTH
    // basis sets and prove if a smaller configuration can be generated.

    for(SymOp oper((b1.sym_equiv(X) || b2.sym_equiv(X)),
		   (b1.sym_equiv(Y) || b2.sym_equiv(Y)),
		   (b1.sym_equiv(Z) || b2.sym_equiv(Z)));
		   oper.in_range();++oper)
    {
      Pair<SymOp> p(s1,s2);
      apply(p,oper);
      if (p < q) 
      {
	prim = false;
	break;
      }
    }
    if (prim)
    {
      if (debug)
	cout << q << " is primitive. " << endl;
      primitive.add(q);
    }
  }  
}

/*TEX
\subsection{Generate Sources}
\index{OneElectronIntegrals, generate_sources}
\index{sources}
In the previous section we only created the minimal set of pairs of
\name{SymOp's}, throwing
away equivalent configurations. Here we now generate all 
equivalent configurations (referred to as sources), 
which can be constructed from a
given 'primitive' pair by a symmetry operation acting on it.
\name{sort_in} ensures, that each source is only
stored once. We store the related \name{SymOp}, too.
\index{sort_in}
*/

int sort_in(SStack<Pair<SymOp>,max_sym >& stack, Pair<SymOp>& p)
{
  for(int i=0; i < stack.size(); i++)
    if (stack[i] == p) return false;
  stack.add(p);
  return true;
}

/*TEX
\index{generate_sources}
*/

void OneElectronIntegrals::generate_sources(Pair<SymOp>& q)
{  
  sources.reset();
  operation.reset();
      
  for(SymOp oper((b1.sym_equiv(X) || b2.sym_equiv(X)),
		 (b1.sym_equiv(Y) || b2.sym_equiv(Y)),
		 (b1.sym_equiv(Z) || b2.sym_equiv(Z)));oper.in_range();++oper)
  {
    Pair<SymOp> p(q);
    apply(p,oper);
    if (debug)
      cout << " Testing Source: " << p << " " << oper << endl;  
    if (sort_in(sources,p))
      operation.add(oper);
  }  

  if (debug)
  {
    cout << " +++ +++ Sources for " << q << endl;
    for(int i = 0;i < sources.size(); i++)
      cout << "                    " << sources[i] 
	   << " Oper: " << operation[i] << endl;
    cout << endl;
  }  
}

/*TEX
\subsection{Distribution of Sources into Targets}
\index{OneElectronIntegrals, distribute}

Having determined all allowed source configurations, we have to
collect their contributions to the target integrals. Therefore we use
a copy of the storage that contains the result of the recursion
for the related inequivalent primitive configuration. For each
source configuration we can evaluate with which prefactor and sign
the result for the inequivalent configuration must be added into the
target \name{IntegArray}. The memory for the target integrals in
\name{tconf} is allocated at the first time it is needed.
*/

void OneElectronIntegrals::distribute(Storage* ss,Pair<Momentum>& m,
				   Pair<SymOp>& source,SymOp& oper,Momentum mu)
{
  Storage permute(*ss,0,tconf[0] -> allocator());
  assert(ss -> size1()  == permute.size1());
  assert(ss -> size2()  == permute.size2());
  assert(ss -> size3()  == permute.size3());
  assert(ss -> size4()  == permute.size4());
  assert(ss -> blocks() == permute.blocks());

  if (debug)
    cout << " +++ +++ Distribution for: " << m << " " << source 
         << " Operation: " << oper << endl;
  
  
  for(int itar = 0; itar < tconf.size(); itar++)
  {
    IntegArray *tar = tconf[itar] -> find(m(0),m(1),0,mu);
    if (tar == 0)
    {
      // IntegArray* tt = new IntegArray(permute.size1(),permute.size2(),
      //			         permute.size3(),permute.size3(),
      //				 permute.blocks());
      IntegArray* tt = new IntegArray(tconf[itar] -> allocator());
      assert(tt != 0);
      assert(permute.size1()  == tt -> size1());
      assert(permute.size2()  == tt -> size2());
      assert(permute.size3()  == tt -> size3());
      assert(permute.size4()  == tt -> size4());
      assert(permute.blocks() == tt -> blocks());
      tconf[itar] -> insert(m(0),m(1),0,mu,tt);
      tar = tt;
    }
    double coef = compute_factor(target[itar],m,source,oper,mu);
        
    if (debug)
      cout << " +++ +++ for target: " << target[itar] 
	   << " Factor: " << coef << endl;
    tar -> add_to(coef,permute);  
    if (debug)
      cout << " Target: " << *tar << endl;
    
  }  
}


/*TEX
\subsection{Compute Factor to contract a Source into a Target}
\index{OneElectronIntegrals, compute_factor}

Given a target \name{SymDesignator} and a \name{SymOp} characterizing 
a primitive configuration we must determine the prefactor with which the 
primitive configuration enters into the target configuration. We do not derive
the prefactors here, but simply state the result.
For each direction
we gather:
\begin{itemize}
\item For normalization a factor of $1/\sqrt{2}$ for each index
    	decorated with an $E$ or an $O$.
\item A factor of $-1$ for each index, where the source descriptor indicates 
      a reflection {\bf and} the target descriptor indicates $O$. This is
      obvious from the second right handside term in \ref{source-signs}
\item A factor of $-1$ for each index, where a reflection was required to 
      create the inequivalent primitive. This comes from the reflection
      symmetry of the orbitals. Consider an $p_x$ orbital for an atom located
      in the positive x halsphere. The orbital has odd symmetry with respect 
      to the position. If we create from this atom a 'virtual' atom
      using reflection on the $YZ$ plane, the orbital with respect to the
      virtual position will have an opposite sign.
      This has to be corrected here.
\item At least the symmetry of the one-electron operator may lead to a 
      further sign, if reflection was used to create the source from
      the inequivalent primitive.
\end{itemize}
*/
double OneElectronIntegrals::compute_factor(Pair<SymDesignator>&  tar,
					    Pair<Momentum>     &  m,
					    Pair<SymOp>        &  src,
					    SymOp              &  oper,
                                            Momentum              mu)
{
  Pair<SymOp> orig(src);
  apply(orig,oper);
  assert(orig == sources[0]);
  
  double root = 1.0/sqrt(2.0);
  double fac  = 1.0;
  
  Direction dir = X;
  int i;
  for(i = 0 ; i < 2; i++)
  {
    if (tar[i](dir) != none)                                 fac *= root;
    if (tar[i](dir)  == odd && src[i](dir) == reflection)    fac *= -1;
    if (orig[i][dir] == reflection && (m[i](dir) % 2 != 0))  fac *= -1;
  }
  if(oper(dir) == reflection)              fac *= operator_fac(dir,mu);

  dir = Y;  
  for(i = 0 ; i < 2; i++)
  {
    if (tar[i](dir) != none) fac *= root;
    if (tar[i](dir) == odd && src[i](dir) == reflection)    fac *= -1;
    if (orig[i][dir] == reflection && (m[i](dir) % 2 != 0)) fac *= -1;
  }
  if(oper(dir) == reflection)              fac *= operator_fac(dir,mu);

  dir = Z;  
  for(i = 0 ; i < 2; i++)
  {
    if (tar[i](dir) != none) fac *= root;
    if (tar[i](dir) == odd && src[i](dir) == reflection)     fac *= -1;
    if (orig[i][dir] == reflection  && (m[i](dir) % 2 != 0)) fac *= -1;
  }
  if(oper(dir) == reflection)              fac *= operator_fac(dir,mu);

  return fac;
} 

/*TEX
\subsection{Output}
\index{OneElectronIntegrals, output}
At least we have to copy our results into the array of output matrices.
This is done by calling the output function for the \name{RecursionStorage},
which itself call the output function of the underlying \name{Storage}
class. This function is defined virtual and we give here the default
routine. For angular momentum and spin-orbit interaction special output
functions will be derived, which have to take care of an extra sign.
*/

void OneElectronIntegrals::output(int l1, int l2)
{
	for(int itar = 0; itar < tconf.size(); itar++)
	{
		//    cout << " +++ +++ Output for : " << itar << " "  
		//	 << target[itar] << endl;
		for(int mutot = min_mu; mutot <= max_mu; mutot++)   
		{
			for(Momentum mu(mutot);mu();++mu)    
				//for all mu with mutot <= max_mu
			{
				//	cout << " +++ +++ +++ mu : " << mu << endl;
				tconf[itar] -> output(b1,l1,target[itar][0],
					b2,l2,target[itar][1],
					mat_array,mu,same);
			}
		}
	}
}

/*TEX
\subsection{Evaluation of the Recursion Relation \label{Eval-Rec}}  
The operator dependend recursion relations are started by the
following \name{generate} functions. They have got all the same structure.
First a set of coefficients that are used many times in the recursion
are precomputed. Then the recursion is initiated by calling the
constructor of the related recursion class.
For kinetic energy, angular momentum and spin-orbit
interaction integrals prior recursions have to be done before.
Finally the primitive gaussian basis functions are contracted to orbitals.
We state here the special output function for angular momentum and
spin-orbit interaction functions, too. They are similar to the general
output function but call the output function of the \name{RecursionStorage}
class with the additional sign \name{m_sym} $ 0 -1$.

\index{OverlapOp, generate}
*/
RecursionStorage* OverlapOp::generate(int l1,SymOp& s1,int l2,SymOp& s2)
{
  if(debug)
  {
    cout << " +++ Overlaps for: " << l1 << " " << s1 << " " 
	 << l2 << " " << s2 << endl;
  }
  
  Coefficient_Set cset(b1.radial_functions(l1),l1,s1,
		       b2.radial_functions(l2),l2,s2,0,OVLP);
 
  OverlapRec* rec = new OverlapRec(cset,l1,l2);
  rec -> contract();  
  return rec;
}

/*TEX
\index{KineticOp, generate}
*/
RecursionStorage* KineticOp::generate(int l1,SymOp& s1,int l2,SymOp& s2)
{
  if(debug)
  {
    cout << " +++ Kinetic for: " << l1 << " " << s1 << " " 
	 << l2 << " " << s2 << endl;
  }
  
  Coefficient_Set cset(b1.radial_functions(l1),l1,s1,
		       b2.radial_functions(l2),l2,s2,0,KIN);
 
  OverlapRec ovlp(cset,l1,l2);
  KineticRec* rec = new KineticRec(cset,l1,l2,&ovlp);
  rec -> contract();
  return rec;
}

/*TEX
\index{MomentOp, generate}
*/
RecursionStorage* MomentOp::generate(int l1,SymOp& s1,int l2,SymOp& s2)
{
  if(debug)
  {
    cout << " +++ Moment for: " << l1 << " " << s1 << " " 
	 << l2 << " " << s2 << endl;
  }
  
  Coefficient_Set cset(b1.radial_functions(l1),l1,s1,
		       b2.radial_functions(l2),l2,s2,0,MOM);

  if(debug)
    cout << "NEW MOMENT(CSET,l1=" << l1 << ",l2=" << l2 << "max_mu="
    << max_mu << " loc " << loc << endl;
  
  MomentRec* rec = new MomentRec(cset,l1,l2,max_mu,loc);

  rec -> contract();
  return rec;
}

/*TEX
\index{NuclearOp, generate}
*/
RecursionStorage* NuclearOp::generate(int l1,SymOp& s1,int l2,SymOp& s2)
{
  if(debug)
  {
    cout << " +++ Nuclear for: " << l1 << " " << s1 << " " 
	 << l2 << " " << s2 << endl;
  }
  
  Coefficient_Set cset(b1.radial_functions(l1),l1,s1,
		       b2.radial_functions(l2),l2,s2,0,NUC);

  NuclearRec* rec = new NuclearRec(cset,l1,l2,loc,chg);

  rec -> contract();
  return rec;
}

/*TEX
\index{E_FieldOp, generate}
*/
RecursionStorage* E_FieldOp::generate(int l1,SymOp& s1,int l2,SymOp& s2)
{
  if(debug)
  {
    cout << " +++ E_Field for: " << l1 << " " << s1 << " " 
	 << l2 << " " << s2 << endl;
  }
  
  Coefficient_Set cset(b1.radial_functions(l1),l1,s1,
		       b2.radial_functions(l2),l2,s2,0,E_FLD);

  E_FieldRec* rec = new E_FieldRec(cset,l1,l2,loc);

  rec -> contract();
  return rec;
}

/*TEX
\index{AngularOp, generate}
*/
RecursionStorage* AngularOp::generate(int l1,SymOp& s1,int l2,SymOp& s2)
{
  if(debug)
  {
    cout << " +++ Angular Momentum for: " << l1 << " " << s1 << " " 
	 << l2 << " " << s2 << endl;
    cout << " +++++++++ COEFFICIENT SET +++++++++ " << endl;
  }
  
  Coefficient_Set cset(b1.radial_functions(l1),l1,s1,
		       b2.radial_functions(l2),l2,s2,0,ANG);
  if(debug)
    cout << " ++++++++ OVERLAP RECURSION ++++++++ " << endl;
  OverlapRec ovlp(cset,l1,l2);
  if(debug)
    cout << " ++++++++ ANGULAR RECURSION ++++++++ " << endl;
  AngularRec* rec = new AngularRec(cset,l1,l2,loc,&ovlp);
  rec -> contract();
  return rec;
}

/*TEX
\index{AngularOp, output}
*/

void AngularOp::output(int l1, int l2)
{
  for(int itar = 0; itar < tconf.size(); itar++)
  {
//    cout << " +++ +++ Output for : " << itar << " "  
//	 << target[itar] << endl;
    
    for(int mutot = min_mu; mutot <= max_mu; mutot++)    
      for(Momentum mu(mutot);mu();++mu)    
	//for all mu with mutot <= max_mu
      {
	// cout << " +++ +++ +++ mu : " << mu << endl;
	tconf[itar] -> output(b1,l1,target[itar][0],
			      b2,l2,target[itar][1],
			      mat_array,mu,same,-1);
      }
  }
}
 
/*TEX
\index{SpinOrbOp, generate}
*/
RecursionStorage* SpinOrbOp::generate(int l1,SymOp& s1,int l2,SymOp& s2)
{
  if(debug)
  {
    cout << " +++ Spin_Orbital Interaction  for: " << l1 << " " << s1 << " " 
	 << l2 << " " << s2 << endl;
    cout << " +++++++++ COEFFICIENT SET +++++++++ " << endl;
  }
  
  Coefficient_Set cset(b1.radial_functions(l1),l1,s1,
		       b2.radial_functions(l2),l2,s2,0,SPINORB);
  if(debug)
    cout << " ++++++++ NUCLEAR RECURSION ++++++++ " << endl;
  NuclearRec nuclear(cset,l1,l2,loc,1);
  if(debug)
    cout << " ++++++++ SPINORB RECURSION ++++++++ " << endl;
  SpinOrbRec* rec = new SpinOrbRec(cset,l1,l2,loc,&nuclear);
  rec -> contract();
  return rec;
}

/*TEX
\index{SpinOrbOp, output}
*/

void SpinOrbOp::output(int l1, int l2)
{
  for(int itar = 0; itar < tconf.size(); itar++)
  {
//    cout << " +++ +++ Output for : " << itar << " "  
//	 << target[itar] << endl;
    
    for(int mutot = min_mu; mutot <= max_mu; mutot++)    
      for(Momentum mu(mutot);mu();++mu)    
	//for all mu with mutot <= max_mu
      {
	// cout << " +++ +++ +++ mu : " << mu << endl;
	tconf[itar] -> output(b1,l1,target[itar][0],
			      b2,l2,target[itar][1],
			      mat_array,mu,same,-1);
      }
  }
}



