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
% $Id: rec.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\index{auxiliary indices}

Given a set of primitive gaussians $\phi_(\zeta_\alpha,R_\alpha,\va)$ 
(where $0 \le \alpha \le n$) 
with common angular momentum $\va$ and a second such set of 
primitive gaussians $\phi'_(\zeta_\beta,R_\beta,\vb)$ (where $0 \le \beta \le m$) 
with common angular momentum  $\vb$ we can introduce a $n \times m$
matrix, which holds the matrix elements between all of functions
of (\alpha,\beta) as:
\be 
(\va\,|\,O(\vmu,m)\,|\,\vb) [\alpha\beta] \label{def-integarray}
\ee 
The expression: $(\va\,|\,O(\vmu,m)\,|\,\vb)$ is the name of the
matrix, which we henceforth call an \name{IntegArray}. 
\index{IntegArray, definition}

For angular pairs with $\va=\vb=0$, so called primivtive integrals,
explicit formulas exist for the integrals we wish to compute. 
\index{integrals, primitive}
\index{primitive integrals}

The integrals for other pairs of angular momenta are recursively defined
in terms of angular momentum recursion relations. The terms of
the recursion relation are conveniently in terms of the matrices
defined in equation~\ref{def-integarray}, which therefore constitute a
core concept of this program. All integrals are expressed in a
multi-dimensional array, which separates the angular and the
corresponding radial indices. Suppressing the radial indices, the
recursion relations express the matrix $(\va\,|\,O(\vmu,m)\,|\,\vb)$
in terms of matrices of the type: $(\va'\,|\,O(\vmu',m')\,|\,\vb')$,
where the total angular momentum of the indices $\va'$ and $\vb'$ is
lower than the total momentum of $\va$ and $\vb$. Note that the
auxiliary indices $\vmu'$ and $m'$ can differ from $\vmu$ and $m$
respectively. Furthermore, on occasion, the recursion relations for an
operator $O$ involve the matrix elements of another, auxiliary
operator $O'$, which is again defined by a set of recursion
relations. Most importantly however, we note that the radial indices
involved in {\bf all terms} of the recursion relation are the same, so
that all operations can be understood in terms of operations on the
matrices defined above.

Since there is a great deal of similarity among the recursion
relations for all the operators relevant to this program we develop
one abstract base class which incorporates the general principles for
the evaluation of the recursion relations. From this class we then
specialize through class derivation the specific recursion relations
for individual operators. It may at first seem surprising, but even
the recursion relations for two-particle operators can be expressed in
this class, which therefore may rightfully be called the cornerstone
of this program. The following summarizes our concept for the 
\name{Recursion} class:

\begin{itemize}
\item Each matrix $(\va|O(\vmu,m)|\vb) [\alpha\beta]$ is represented as an 
      \name{IntegArray}.
\item A \name{Recursion} is a vector of (pointers to) such \name{IntegArray}s.
      Each element of this vector corresponds to a specific combination 
      of angular momenta and auxiliary indices. The member function 
      \name{find} can be used to obtain the pointer to a specific 
      \name{IntegArray}. It returns zero if the array in question has 
      not been initialized.
\item To implement the recursion relation for a given \name{IntegArray},
      we need to 
      \begin{itemize}
      \item Compute the required terms, i.e. \name{IntegArray}s,
            of lower angular momentum.
      \item Multiply each such term element-by-element with 
            a scalar angular dependent factor and/or a 
            a radial-function dependent factor $c_{\alpha\beta}$ and 
            add the result to the target \name{IntegArray}.
      \end{itemize}
\item A matrix of radial-function dependent prefactors 
      $c_{\alpha\beta}$ is stored in a class of type \name{CoefArray}. 
      \index{CoefArray}
\item \name{CoefArrays} and \name{IntegArrays} are sufficiently close
      in structure to be derived from a single basic class, 
      named \name{Storage}. This class implements all low-level operations.
\item The set of all such prefactors is maintained in a class of type
      \name{CoefficientSet}, which is passed as an argument to 
      \name{Recursion}. All required \name{CoefArrays} are precomputed
      from the radial-function arrays before \name{Recursion} is initiated.
\item The evaluation of each individual term is therefore a very fast 
      element-by element multiply-and-add operation on the matrices.
\item The evaluation of all terms required to generate a set of 
      \name{IntegArray}s of given total angular momenta $(l_1,l_2)$ 
      is carried out  in the member function \name{reduce}, 
      which is discussed in detail below. This member function assumes that
      all the \name{IntegArray}s required in the individual terms are 
      already in present in the \name{Recursion} class. The function returns
      immediately, when the \name{IntegArray} in question has already been 
      generated. Otherwise a new \name{IntegArray} is created and its values 
      computed. 
\item The availability of the required terms is ensured by the 
      \name{logic} member function which controls the order in which
      the \name{IntegArrays} are evaluated.
      \name{logic} operates on the level of angular momentum shells
      rather than individual angular momenta. An angular momentum shell is the 
      set of all angular momenta of a given total angular momentum.
      \name{logic} may call itself recursively for all required 
      terms of lower angular momenta. On return from these calls 
      it calls \name{reduce} which carries out the reduction.
\item Presently the decision, which route of reduction to follow is simple:
      {\em We always reduce the angular index with the larger total angular 
      momentum}. As a result, a call to \name{logic} with the two target 
      angular momenta will automatically generate all required terms.
      This procedure could be surely optimized.
\end{itemize}
\index{recursion relations, logic}
\index{logic}
\index{reduce}

\section{Recursion Formulas and the \name{Recursion} Class}

\subsection{Generalized Recursion Formulas}
\index{recursion relations, generalized}

Before reading this section the reader should consult the appendix of
\cite{OS}. On inspection of the recursion relations summarized there,we note
a great deal of similarity among them, which we would like to exploit
in our implementation. We therefore implement a {\em abstract class}
\name{Recursion}, which represents the overall structure for a generic
recursion relation. From this {\em abstract class} we then derive
recursion classes for the individual operators we want to implement,
augmenting the basic instrumentarium of \name{Recursion} with the
specifics of the individual operator. Mathematically, the similarities 
in the recursion relations are summarized as follows:

\subsection{General Terms\label{section-general-terms}}

Consider the angular recursion relation of an operator $O$, which
depends on an angular momentum index: $\vmu$ and an auxiliary index
$m$. We list the generic expression for a reduction along direction
$\one_i$, all terms which do not fit this specification are lumped
into the special term $R(\va,\vb,\vmu,m)$:
\ba
  (\va\,|\, O(\vmu,m) \,|\,\vb) & = &  \nonumber \\ 
  & f_0 & (P_i- A_i)\  \ (\va'\,|\, O(       \vmu,  m) \,|\,\vb')  \nonumber \\
+ & f_1 & (P_i- C_i)\ \ (\va'\,|\, O(       \vmu,m+1) \,|\,\vb')  \nonumber \\
+ & f_2 & \{\tilde\zeta_s\} 
           N_i(\vmu) \ \ (\va'\,|\, O(\vmu-\one_i,m+1) \,|\,\vb')  \nonumber \\
+ & f_3 & \tilde\zeta_s 
        \ [\, N_i(\va')\ (\va'-1\,|\, O(\vmu,  m) \,|\,\vb')    \nonumber \\
  &     & \ \ \ \   + N_i \ (\vb') (\va'  \,|\, O(\vmu,  m)
     \,|\,\vb'-1)\ ] \nonumber \\
+ & f_4 & \tilde\zeta_s \ [\  N_i(\va)\  (\va'-1\,|\, O(\vmu,m+1) \,|\,\vb')    \nonumber \\
  &     & \ \ \ \   + N_i(\vb') (\va'  \,|\, O(\vmu,m+1) \,|\,\vb'-1)
     \  ] \nonumber \\
+ & R(\va,\vb,\vmu,m)  & 
\ea
Here we have defined:
\ba
      \va' & = & \left\{ \begin{array}{ll}
                        \va - \one_i & {\rm \ \ for\ reduction\ on\ a} \\
                        \va          & {\rm \ \ for\ reduction\ on\ b} \\ 
                      \end{array} \nonumber \right. \\
      \vb' & = & \left\{ \begin{array}{ll}
                        \vb   & {\rm \ \ for\ reduction\ on\ a}  \\
                        \vb - \one_i& {\rm \ \ for\ reduction\ on\ b}  \\
                      \end{array} \right.
\ea
\index{$\va'$}
\index{$\vb'$}
For two reason we thereby have to make a difference between using the 
auxiliary index $m$ or not. First of all increasing the auxiliary index,
when it is actually not used will not be meaningfull and second the
prefactor in the curly brackets in the term of $f2$ does only occur for
the case, when $m$ is not used (e.g. for the moment integrals). 
\index{auxiliary index, use of}

The integer prefactors $f_k$ depend on the specific operator as follows: 
\index{f prefactors,definition}
\index{f prefactors,table}
\begin{table}[h]
\caption{Flags for The Evaluation of the  General Terms in Single 
Particle Recursion Relations \label{term_flags}}
\bc
\begin{tabular}{|l|rrrrr|}
\hline
           &$f_0$&$f_1$&$f_2$&$f_3$& $f_4$ \\ 
\hline
OVERLAP    & 1 & 0 & 0 & 1 &  0  \\
MOMENT     & 1 & 0 & 1 & 1 &  0  \\ 
KINETIC    & 1 & 0 & 0 & 1 &  0  \\
NUCLEAR    & 1 & -1 & 0 & 1 & -1  \\
ELECTRIC   & 1 & -1 & 1 & 1 & -1  \\
ANGULAR    & 1 & 0 & 0 & 1 &  0  \\
SPINORBIT  & 1 & -1 & 0 & 1 & -1  \\
\hline
\end{tabular}
\ec
\end{table}

\subsection{Special Terms}
\index{special terms}
\index{recursion relation, special terms}

Special terms appear in the recursion relations for KINETIC energy,
ANGULAR momentum and SPINORBIT coupling matrix elements and are
specially coded for each individual case. Each of these special terms
involves an {\em auxialiary operator}, which is passed as a parameter
in a separate auxialiary recursion class. 

Here we define the required angular prefactors as: 
\ba
      \va''& =& \left\{ \begin{array}{ll}
                      \va - 2 \cdot \one_i & {\rm \ \ for\ reduction\ on\ a} \\
                      \va & {\rm \ \ for\ reduction\ on\ b} \\
                      \end{array} \right. \\
      \vb''& = & \left\{ \begin{array}{ll}
                      \vb &  {\rm \ \ for\ reduction\ on\ a}  \\
                      \vb - 2 \cdot \one_i &{\rm \ \ for\ reduction\ on\ b}  \\
                      \end{array} \nonumber \right. 
\ea

\section{Integrals over s functions}
\index{recursion relations, s functions}
In the recursion relations for MOMENT, ELECTRIC FIELD, 
ANGULAR momentum and SPINORBIT coupling matrix elements integrals over
s functions with $\vmu$ greater than $\zero$ occur. In order to reduce
these integrals to the case of $\vmu = \zero$ further specific
recursion rules have to be implemented. We call this $\mu$-reduction.

\section{Strategy of Implementation}
We implement classes to evaluate the recursion relations for each the
various operators, each of which is derived from the abstract class
\name{Recursion}, which implements the generic recursion relations as
described in the previous section. It allows the incorporation of
special terms through the \name{special} and \name{speclogic} member
functions. 
\index{special}
\index{speclogic}
If $\mu$-reduction is neccessary for ss-integrals, this can
be done by \name{mureduce} and \name{mulogic}.
\index{mulogic}
\index{mureduce}

The storage of the \name{IntegArray} pointers is
implemented in a base class \name{RecursionStorage}, see
section~\ref{section:recstor}.

The concept of the recursion class is explained in the sections on the
constructor, the function \name{logic} and \name{reduce}. The basic
concept is as follows: 
\index{logic}
\index{reduce}
We define a set of pointers to the \name{IntegArrays}, which hold the
values of the operators and a function \name{find} to obtain the
pointer to one specific set of angular momentum and auxiliary indices.
Initially all these pointers are set to zero, indicating that the
fields carry no values.

We implement a (possibly recursive) function \name{logic}, which
computes the integrals for the angular momentum shell of its
arguments. (For simplicity we ignore the auxiliary index $m$ for the moment.) 

Logic terminates for the (ss) angular momentum shell
starting the operator-specific $\mu$-reduction if $\vmu \neq \zero$
or calling the primary integral \name{ss_integ} else.
The latter function
must create and fill the \name{IntegArray} for the primary integral
and check the resulting array into the recursion class.

For shells involving nonzero angular momenta $(l_1,l_2)$, the
recursion relation involves angular momentum shells such as
$(l_1-1,l_2)$, $(l_1-2,l_2)$, $(l_1,l_2-1)$ etc. Therefore the
recursion relation can not be applied before the values for these
angular momentum shells have been computed. This is achieved by first
calling the \name{logic} function for the required sets of angular
momentum pairs. Once all the source-terms for the shell (l1,l2) have
been computed we create new \name{IntegArrays} for this angular
momentum shell and compute their values by carrying out the recursion
relation. This is done in the function \name{reduce}.

The table in the previous section lists all possible terms for generic
recursion relations. Depending on the values of the $f_i$ flags,
\name{logic} and \name{reduce} are therefore able to determine which
terms must be computed. To costumize the recursion relation for
specific operators we define two complementary functions
\name{speclogic} and \name{special}. \name{speclogic} implements the
operator-specific calls to \name{logic} for the special terms 
$R(\va,\vb,\vmu,m)$, which are not encoded in the
table in the previous section. This function is always called 
from the \name{logic} routine, {\bf
before} the \name{reduce} function is called.  In analogy to the
\name{logic} routine, \name{speclogic} must ensure the presence
of {\bf all required special terms}. The contributions of the special
terms are eventually computed in the \name{special} recursion
routine. The \name{special} routine is called from within the angular
momentum loop of the reduce routine and evaluates the special terms
for one specific set of angular momenta, not for an entire shell. See
\name{reduce} for details.

The structure of the $\mu$-reduction is similar. 
It is startet out of \name{logic} by calling \name{mulogic} 
with $l_1 = l_2 = 0$ and $\mu_{tot} > 0$. \name{mulogic} ensures that
all  source terms with lower $\mu_{tot}$ are calculated first.
In case of $\mu_{tot} = 0$  again \name{ss_integ} is called.
Finally \name{mulogic} calls \name{mureduce}, which evaluates the 
recursion rules for the whole shell belonging to $\mu_{tot}.

The operator dependent functions \name{ss\_integ},\name{speclogic},
\name{special},\name{mulogic} and \name{mureduce} are defined as 
abstract members in \name{Recursion}.
They {\bf must be defined} in the derived classes implementing
specific operators.

\section{Class Header Recursion}

Using the abstract base class it is possible to instantiate a
number of derived classes which implement specific recursion relations
for a  desired single particle operators. We begin with the
description of the base class:

\subsection{Recursion Class Header}
\index{Recursion, class header}
*/
#ifndef RECURSION_H
#define RECURSION_H

extern Location DummyCenter;

class InternalBasis;

#include "rec_stor.h"
 
class Recursion : public RecursionStorage
{
protected:
  Allocator alloc_contract;
  

  int        f0,f1,f2,f3,f4; // the integer flags indicating the
			     // presense of the individual standart
			     // terms 
  
  Coefficient_Set&     coef_set;  // the coefficient set associated
				  // with the recursion relation
  Recursion*           auxop;     // a pointer to the recursion
				  // relation of an auxialiary
				  // operator if any is present, NULL
				  // otherwise 
  Location             center;    // for recursion relations which
				  // depend on an auxiliary center
				  // the coordinates of the center

public:
         Recursion(Coefficient_Set& CSet,
		   int momentum1,int momentum2,int max_aux,int max_mu,
		   Location cntr,Recursion* auxop);
        ~Recursion() { clean(); } // make sure all IntegArrays are gone }
         void  logic    (const int l1,const int l2,const int aux,const int mu,
                         int use_aux = 1, int minmu = 0);
         void  reduce   (const int l1,const int l2,const int aux,const int mu,
                         int use_aux = 1, int mumin = 0);
virtual  void ss_integ  (const int aux,const int mu) = 0;
virtual  void speclogic (const int l1,const int l2,const int aux,const
			 int mu) = 0;
virtual  void special   (IntegArray* target,const int position,
			 const Direction& cir,
			 const Momentum& a,const Momentum& b,
			 const Momentum& mu,const int aux) = 0;
virtual  void mulogic   (const int l1,const int l2,const int aux,const
			 int mu) = 0;
virtual  void mureduce  (const int l1,const int l2,const int aux,
                         const int mu) = 0;

void     contract();	 // contract the primitive gaussians to
			 // orbitals.
};
#endif

