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
\chapter{BasisSets\label{chapter-basis}}
\section{Introduction}
This program deals with the evaluation of single and two-particle
integrals between contracted gaussian basis functions. A {\bf normalized}
primitive cartesian gaussian basis functions is defined as:
\be
  \phi_{R,\zeta,\vec{n}} (\vec{x}) =  {\cal  N} (\zeta,\vec{n})
  (x-R_{x})^{n_x} (y-R_{y})^{n_y} (z-R_{z})^{n_z} 
  \exp{[-\zeta (\vec{r} - \vec{R})^2}],
\ee
which is characterized by the location $\vec{R}$,the exponent
$\zeta$ and cartesian angular momentum $x^{n_x}  y^{n_y}  z^{n_z}$.
Linear combinations of such functions may be used to construct
spherical gaussian basis functions. The normalization constant is
given as:
\be
  {\cal N } (\zeta,\vec{n}) =
 \left[ (2 n_x - 1)!!\ (2 n_y - 1)!!\ (2 n_z - 1)!! \right]{1/2}
 \left(\frac{2 \zeta}{\pi} \right)^{3/4}  (4 \zeta)^{l/2}  \label{eq-norm}
\ee
 
\index{GTO} \index{degree of contraction} \index{contraction} 
\index{gaussian, primitive}
\index{gaussian, contracted}
\index{gaussian,normalization}
\index{radial dependence}
  
A contracted basis function is a linear combination of primitive
gaussian basis functions of the same angular momentum and
on the same center. A gaussian type orbital, GTO, is thus defined
as:
\be  
   \sum_i^{D} \alpha_i  \phi_{R,\zeta_i,\vec{n}} (\vec{x})
\ee
Note that $\vec{n}$ and $R$ do not vary  with $i$.
$D$ is called the degree of contraction.

\section{Pure Basis Functions}

In the context of a quantum chemistry calculation, the basis set is
approached typically from a different angle. For a given nuclear type,
we define a set of contracted radial basis functions, for the s,p,d,
$\dots$ orbitals. These sets of wavefunctions are then attached as a
{\bf complete set} to the centers of the given type. In other words,
we start with a set of basis functions, which do not yet have a
``center'' attached. In this program, we refer to this type of 
contracted basis function as a \name{PureBasisFunction}. 

The following class stores the relevant information for a single pure
basis functions, which is characterized by a specific total angular
momentum. Radial components are added using the \name{add} function,
which takes the exponent and the contraction coefficient as arguments.
\index{basis function, pure}
*/
#ifndef BASIS_H
#define BASIS_H

#include <numer.h>
#include <tsstack.h>
#include "typesip.h"

class PureBasisFunction 
{
  SStack<double,max_contr>  exp;
  SStack<double,max_contr>  cc;
  int             mm;
public:
        PureBasisFunction(int momentum)  { mm = momentum; }
       ~PureBasisFunction() { }
void    add(double ee,double coef)
        { exp.add(ee); cc.add(coef);  } 
int     contraction()            const { return exp.size(); }
double  exponent   (const int n) const { return exp(n); }
double  coefficient(const int n) const { return cc (n); }
int     momentum()               const { return mm; }
void    print_gaussian_format(ostream& os,int& no_bas,int& no_prim) const;    
};

ostream& operator<<(ostream& os,const PureBasisFunction& bf);
istream& operator>>(istream& os,PureBasisFunction& bf);
/*TEX
\section{Pure Basis Sets}

Having defined pure basis functions, we now proceed to
\name{PureBasisSets}, which are merely a list (stack) of (pointers to)
\name{PureBasisFunctions}. New basis functions are added to the stack
simply by adding them using the function \name{add}. Functions of
different angular momentum and different degree of contraction may be
freely mixed.
\index{basis set, pure}
*/

class PureBasisSetList;

class PureBasisSet : public SStack<PureBasisFunction*,max_pbf>
{
  int printed;
public:
      PureBasisSet() { printed = false; } 
     ~PureBasisSet();
int   max_momentum(); 
int   no_function();    
int   parse(istream& is,PureBasisSetList* list = NULL);
void  print(String& name,ostream& os);    
void  print_gaussian_format(ostream& os,int& no_bas,int& no_prim);    
};

ostream& operator<<(ostream& os,PureBasisSet& bf);

/*TEX
\section{PureBasisList}
This is merely a class which holds a list of all known pure basis sets for 
future use. Each basis set is associated with a unique name, given as a 
string. The introduction of this class-hierarchy is motivated by the
desire to implement a {\em library} of pure basis sets. In most
integral programs I know of the user must copy the information
concerning the basis set into the input file for the program, a
process which is often error prone. Yet the information in the
optimized basis sets never changes. Hence we wish to put all pure 
basis sets for all nuclear types and once and for all into a file. 
The integral program can then parse this file and the user may refer
to the individual basis sets by their symbolic names. 

\index{basis  set, list of pure}
\index{basis  set, library}

*/

class PureBasisSetList 
{
  SStack<       String,max_pbs>   name;
  SStack<PureBasisSet*,max_pbs>     bs;

public:
      PureBasisSetList() { } 
     ~PureBasisSetList();

PureBasisSet* find(String& s)
      {
        for(int i=0; i< name.size(); i++)
	  if(name[i] == s) return bs[i];
	return NULL;
      }
int   add(String& newname,PureBasisSet* set)
      { 
	if(find(newname))
	{
	  cerr << " +++ PureBasisSetList::add PureBasisSet named: " << newname 
	       << " already present. Not Added"  << endl;
	  return -1;
        }
	name.add(newname); return bs.add(set); 
      }
static int   print_flag;  

};

istream& operator>>(istream& is,PureBasisSetList& blist);

/*TEX
\section{Basis Sets}

Having introduced the concept of a \name{PureBasisSets}, 
we are now in a position to define \name{BasisSets}. For standard
calculations, a \name{BasisSet} is the highest organizational
structure dealing with basis functions. A \name{Basisset} is simply a
list (stack) of pairs, each consisting of a \name{PureBasisSet} and a
location. New functions can only be added to the basis set by
attaching such pairs with the add function. There may be several 
\name{PureBasisSets} on the same location.  The function \name{parse}
is provided to read a  \name{BasisSet} from an input stream.
\index{Basis Set}
\index{Basis Set, internal variables}
\index{Orbitals, indexing}

\subsection{Enumeration of BasisFunctions} 

To label the integrals we must attach indices to the orbitals in the
basis set. This program uses an indexing scheme based on an ordering
of \name{PureBasisFunctions} in a given basis set:

\begin{itemize}
\item All basis functions on a given center precede all basis
functions on subsequent centers in the \name{BasisSet}.
\item The orbitals on a given center are ordered by their appearance
in the \name{PureBasisSet}. 
\item All orbitals associated with a given
\name{PureBasisFunction} on a specific center are labeled in sequence
by the \name{starting\_index} of the pure basis function plus the
\name{local\_index} of labeling the specific angular momentum.  For
the indices of the individual angular components see the description
of \name{local\_index} of class \name{Momentum}.
\item The \name{starting\_index} of a \name{PureBasisFunction} is
given as the number of all orbitals in \name{PureBasisFunctions}
of higher precedence. 
\end{itemize}

\subsection{Internal Variables}

To compute matrix elements between basis function, we must apply
angular momentum recursion relations for primitive basis functions. We
therefore must group the primitive basis functions of the same angular
momentum. Since this information is used over and over it is stored
permanently in the class, the information is computed using the
function \name{prepare}. The \name{prepared} flag insures that the 
contents of these fields is kept current.

A set of primitive basis function is encoded in a
\name{RadialFunctionArray}, see the corresponding chapter. A
\name{RadialFunction} is the name for a primitive basis function,
complete with \name{Location} and an exponent. In addition we store
for each primitive function the number of the contracted basis
function to which it contributes and the appropriate contraction coefficient.

\name{maxl} designates the the bound on the angular momentum encoded in
the \name{BasisSet}, \name{nof\_total} the total number of {\bf
contracted} basis functions in the \name{BasisSet}.
\index{maxl} \index{nof\_total}

For each angular momentum $0 \le l <$ \name{maxl} we prepare the
following information:
\begin{itemize}
\item \name{nofun[l]} designates the number of contracted functions of
total angular momentum l. \index{nofun}
\item \name{noprim[l]}  designates the number of primitive functions of
total angular momentum l. \index{noprim}
\item The \name{RadialFunctionArray} \name{radials[l]} encodes
locations, exponents and contraction information for the primitive
basis functions of angular momentum l. \index{radial functions}
\item \name{start\_index[l]} is a vector of length \name{nofun[l]}.
The k-th element of this vector gives the index of the first angular
wavefunction associated with the k-th {\bf contracted orbital} of angular momentum
$l$. This orbital is the one with angular prefactor $x^l$.
\item \name{fcenter[l]} is a vector of length \name{nofun[l]}.
The k-th element of this vector gives location  
of the  k-th {\bf contracted orbital} of angular momentum
$l$. 
\end{itemize}
 
\subsection{Class Header}
*/
#include "vtype.h"
#include "function.h"
#include "orbinfo.h"
#include "element.h"

typedef  ARRAY<Location> LocationArray;
class IntegArray;


class BasisSet  
{
  friend class InternalBasis;
  friend class GauBasisSet;
  
  int debug;

  SStack<     Location,max_cntr>      loc;
  SStack<PureBasisSet*,max_cntr>   pbasis;
  SStack<       String,max_cntr>     name;
  
public:
       BasisSet(const int dbg = false) { debug = dbg; }  
      ~BasisSet() {}  
int    add(PureBasisSet* pset,String& nn,Location& center) 
       { pbasis.add(pset); name.add(nn); return loc.add(center); } 
void   parse(istream& is,Molecule& molec,MomentData& mu_data,
	     Location& e_field_loc,Location& anglr_loc,Location& spin_orb_loc,
	     PureBasisSetList& blist);
};

enum BasisSelectionMode  { ALL, ALL_CENTER } ;
class SymDesignator;  

class InternalBasisList;

class InternalBasis
{
  int debug;
  int global_start_index;
  
  /* auxialiary fields: */
  int maxl;

  ARRAY<int>                 nofun;
  ARRAY<int>                 noprim;
  ARRAY<RadialFunctionArray> radials;
  ARRAY<LocationArray>       fcenter;
  ARRAY<ARRAY<String> >      name;
  ARRAY<IVec>                s_index;  

  OrbInfo&                   oinfo;
  Vec                        norm;
  
public:
  InternalBasis(BasisSet& basis,const Location& target,BasisSelectionMode mode,
		OrbInfo& oinf);
 ~InternalBasis();
friend void compute_overlap(InternalBasis& b1,InternalBasis& b2,
                    ARRAY<IntegMat>& matrix,int same);
friend void compute_kinetic(InternalBasis& b1,InternalBasis& b2,
                    ARRAY<IntegMat>& matrix,int same);     
friend void compute_moment(InternalBasis& b1,InternalBasis& b2,
		    SStack<Location,max_cntr>& centers,
		    SStack<double,max_cntr>& charges,
                    int mu,ARRAY<IntegMat>& matrix,int same);
friend void compute_nuclear(InternalBasis& b1,InternalBasis& b2,
		    SStack<Location,max_cntr>& centers,
		    SStack<  double,max_cntr>& charge,
                    ARRAY<IntegMat>& matrix,int same);
friend void compute_e_field(InternalBasis& b1,InternalBasis& b2,
             Location& loc,
//	            SStack<Location,max_cntr>& centers,
//		    SStack<double,max_cntr>& charge,int mu,
                    ARRAY<IntegMat>& matrix,int same);

friend void normalize(InternalBasis& b1,InternalBasis& b2,
		      ARRAY<IntegMat>& matrix,int max_mu);


friend class TwoElectronIntegrals;  
friend class OneElectronIntegrals;  
friend class IntegArray;  
friend class Storage;  
friend class GauBasisSet;

RadialFunctionArray&  radial_functions(int l) { return radials[l]; }  
int    max_momentum() const { return maxl; } 
int    no_function()      const;
int    no_sym_function()  const;
int    min_index()        const { return global_start_index; }  
int    sym_equiv(Direction dir) const; 
int    sym_start_index(SymDesignator& sd) const;
};

ostream& operator<<(ostream& os,const BasisSet& bf);
ostream& operator<<(ostream& os,const InternalBasis& bf);

class InternalBasisList : public  SStack<InternalBasis*,max_ibs> 
{
public:
  ~InternalBasisList();
};
#endif


