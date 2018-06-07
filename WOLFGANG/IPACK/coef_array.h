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
% $Id: coef_array.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\chapter{Coefficient Arrays \& Coeffcient Sets}
\section{Concepts}

We had seen earlier, see class \name{Recursion}, that the terms in
recursion relations for single particle operators can be 
formulated as element-by-element multiplications of a coefficient
matrix with an integral matrix and possibly a scalar factor. The
indices of the matrix are the \name{RadialFunctions} of the orbitals,
which appear on either side of the operator. 

Only a limited set of such coefficient matrices is required to
implement all recursion relations. These matrices can be computed from
the information contained in the \name{RadialFunctionArrays}
describing the radial dependence of the primitive basis functions.

In this program the coefficient-matrices are precomputed before the
the recursion relations are applied, since most terms appear more than
once. Therefore we have designed a special class for the storage of
such matrices, the class \name{CoefArray}. \index{CoefArray,definition}.
This class in turn is derived from the \name{Storage} class which implements 
all basic operations on arrays of integrals depending on up to four
radial indices. These arrays may be symmetric under the exchange of
various indices as described in the discussion of \name{Storage}.
Three types of \name{CoefArrays} are implemented, for two-index
arrays, for uncontracted and for contracted four-index arrays.

Each \name{CoefArray} is labeled by a code, the \name{Coef\_Type}
\index{Coef\_Type,definition} and is constructed from a pair of
\name{RadialFunctionArrays}. The \name{init}-functions to the \name{CoefArray}
generates the required values. All required \name{CoefArrays} are
stored in a \name{Coefficient\_Set}, which is the passed to the
\name{Recursion} relation to evaluate the integrals.

A number of \name{CoefArrays} are required for the evaluation of a
recursion relation. The set of all required \name{CoefArrays} is
called a \name{Coefficient\_Set}, which is discussed following the
\name{CoefArrays}. 

The evaluation of the two-particle ERI's also requires
\name{CoefArrays}. Here, however each index of the matrix labels a
pair of \name{RadialFunctions} as is discussed in
section~\ref{cfarray2p}. We have also  designed a container class to
hold the set of all such four-index \name{CoefArrays}, which are
required for the evaluation of the two-particle recursion relations.
This container class is called a \name{Quadruple\_Set} and discussed
in section~\ref{qset}.

\subsection{Definition of Single Particle CoefArrays}

\index{Coef\_Type,table of} \index{ZETA1} \index{ZETA2}
\index{Coef\_Array,table of} \index{ZETASUM} \index{XI}
\index{$\tilde \zeta_1$} \index{PA1X} \index{PA1Y} \index{PA1Z}
\index{$\tilde \zeta_2$} \index{PA2X} \index{PA2Y} \index{PA2Z}
\index{$\tilde \zeta_s$} \index{PX} \index{PY} \index{PZ}
\index{$ pa_{1x}$}       \index{ABX} \index{ABY} \index{ABZ}
\index{$ pa_{1y}$}       \index{OVERLAP} \index{KINETIC} 
\index{$ pa_{1z}$}       \index{NORM}
\index{$ pa_{2x}$} \index{$ pa_{2y}$} \index{$ pa_{2z}$}
\index{$ p_{x}$}   \index{$ p_{y}$}   \index{$ p_{z}$}
\index{$ ab_{x}$}  \index{$ ab_{y}$}  \index{$ ab_{z}$}
\index{$\tilde \xi$} \index{$c_{o}$} \index{$c_{k}$} \index{$c_{n}$}
\index{$ pc_{x}$} \index{$ pc_{y}$} \index{$ pc_{z}$}

Here we tabulate the formulas for the various \name{Coef\_Types}.
Consider two \name{RadialFunctions} characterized by exponents
$\zeta_{1}$ and $\zeta_{2}$ respectively, which are centered at
locations, which are centered at locations $\vec{A}_1$ and $\vec{A}_2$
respectively. Suppose, that the \name{RadialFunction} characterized by 
$(\zeta_1,\vec{A}_1)$ is the i-th wavefunction in the first 
\name{RadialFunctionArray}, while the \name{RadialFunction} 
characterized by  $(\zeta_2,\vec{A}_2)$ is the j-th wavefunction 
in the second \name{RadialFunctionArray}. To initialize a given
\name{CoefArray} we run through the standard order of indices as
specified in the storage class. This order depends on the symmetries
of the particular array. \index{Coef\_Array,indexing}

The symbols, names and formulas for the \name{CoefArray} are given in
Table~\ref{coefarray1}. These \name{CoefArrays} are initialized by the
function \name{init}, which takes two \name{RadialFunctionArrays} as
an argument and by the corresponding constructor.  To define the
\name{NORM} array, which is used to normalize the primitive gaussians
before contraction, we introduce:
\be
 {\cal N}'(\zeta) = \left(\frac{2 \zeta}{\pi} \right)^{3/4} (4 \zeta)^{l/2}
\ee
Note that the angular prefactor from equation~\ref{eq-norm} is
missing.

\index{CoefArray,init} 

\begin{table}
\label{coefarray1}
\caption{Definition of Single-Particle Coefficient Arrays}
\bc
\begin{tabular}{ccl}
\hline
Symbol &  \name{Coef\_Type} & Definition \\ 
\hline
$\tilde \zeta_1$ & ZETA1   & $\frac{\zeta_2}{\zeta_1 + \zeta_2}$ \\
$\tilde \zeta_2$ & ZETA2   & $\frac{\zeta_1}{\zeta_1 + \zeta_2}$ \\
$\tilde \zeta_s$ & ZETASUM & $\frac{1}{2(\zeta_1 + \zeta_2)}$    \\
$ p_{x}$         & PX      & $\frac{\zeta_1 A_{1x} + \zeta_2 A_{2x}} 
                             {\zeta_1 + \zeta_2}$ \\
$ p_{y}$         & PY      & $\dots$ \\
$ p_{z}$         & PY      & $\dots$ \\
$ pa_{1x}$       & PA1X    & $p_x - A_{1x}$ \\
$ pa_{1y}$       & PA1Y    & $\dots$ \\
$ pa_{1z}$       & PA1Z    & $\dots$ \\
$ pa_{2x}$       & PA2X    & $p_x - A_{2x}$ \\
$ pa_{2y}$       & PA2Y    & $\dots$ \\
$ pa_{2z}$       & PA1Z    & $\dots$ \\
$ ab_{x}$        & ABX     & $A_{1X} - A_{2X}$ \\
$ ab_{y}$        & ABY     & $\dots$ \\
$ ab_{z}$        & ABZ     & $\dots$ \\
$ {\cal N}$      & NORM    & ${\cal N}'(\zeta_1) {\cal N}'(\zeta_2) - 1$ \\ 
$\tilde \xi$     & XI      & $2 \frac{\zeta_1 \zeta_2}{\zeta_1 + \zeta_2}$ \\
$c_{o}$          & OVERLAP & $(\frac{\pi}{\zeta_1 + \zeta_2})^{3/2} 
                             \exp[\xi (\vec{A}_1 - \vec{A}_2)^2]$ \\    
$c_{k}$          & KINETIC & $\xi [3-2\xi (\vec{A}_1 - \vec{A}_2)^2] \ 
                             c_{o}$ 
\end{tabular}
\ec
\end{table}

\subsection{Definitions for Two Particle CoefArrays\label{cfarray2p}}

Similarly to the evaluation of recursion relations for single particle
operators, the evaluation of recursion relations for the ERI's
requires coefficient arrays. These arrays are indexed (implicitly) by
the indices of the four \name{RadialFunctions}, which make up a single
entry in the \name{CoefficientArray}. In order to keep with the matrix
notation, we combine the first pair of indices as one compound index and
the second pair a second compound index. The names, symbols and
formulas for the required \name{CoefArrays} are listed in
Table~\ref{coefarray3}. 


\index{Coef\_Type,definition}
\index{$w_{x}$}  \index{$w_{y}$} \index{$w_{z}$}
\index{$\rho_1$} \index{$\rho_2$} \index{$\rho_s$}
\index{WX} \index{WY} \index{WZ} \index{RHO1} \index{RHO1} \index{RHOSUM}

Suppose we are given a quadruplet of \name{RadialFunctions}
characterized by exponents $\zeta_i$ and locations $\vec{A}_i$, where
$i = 1, \dots, 4$. We then define:
\ba
   \zeta_{12} & = & \zeta_1 + \zeta_2 \\
   \zeta_{34} & = & \zeta_1 + \zeta_2 {\rm \ and} \\
   \rho       & = & \frac{\zeta_{12}\zeta_{34}}{\zeta_{12} + \zeta_{34}} \\
 P_{\alpha,12}& = & \frac{\zeta_1 A_{\alpha,1} + \zeta_2 A_{\alpha,2}}
                         { \zeta_{12} } \\
 P_{\alpha,34}& = & \frac{\zeta_3 A_{\alpha,3} + \zeta_3 A_{\alpha,4}}
                         { \zeta_{34} },
\ea
where $\alpha = x,y,z$.

\begin{table}
\label{coefarray3}
\caption{Definition of Two-Electron Coefficient Arrays}
\bc
\begin{tabular}{|ccl|}
\hline
Symbol &  \name{Coef\_Type} & Defintion \\ 
\hline
$ w_{x}$  &  WX      &   $\frac{\zeta_{12} P_{x,12} + \zeta_{34} P_{x,34}}
                               {\zeta_{12} + \zeta_{34}} $ \\
$ w_{y}$  &  WY      &   $\dots$ \\
$ w_{z}$  &  WZ      &   $\dots$ \\
$ \rho_1$ &  RHO1    &   $-\frac{\rho}{2\zeta_{12}}$ \\
$ \rho_2$ &  RHO2    &   $-\frac{\rho}{2\zeta_{34}}$ \\
$ \rho_s$ &  RHOSUM  &   $\frac{1}{2(\zeta_{12} + \zeta_{34})}$ \\
\hline
\end{tabular}
\ec
\end{table}

All of these \name{CoefArrays} can be computed from precomputed
single-particle \name{CoefArrays} for pairs of
\name{RadialFunctions}. This is much more efficient than a direct
evaluation from the \name{RadialFunctionArrays}. The two particle,
four-index \name{CoefArrays} are therefore initialized by the
constructor, which takes two \name{Coefficient\_Sets} and the desired
\name{Coef\_Type} as arguments.

\section{\name{CoefArrays}}

\subsection{Definition} As described above a \name{CoefArray} is
nothing but a matrix of radial-function dependent coefficients, which
are required in the evaluation of the recursion relations for single
and two-particle operators. We define an {\bf enumeration} of the
required \name{Coef\_Type}s, which associates each name of a
\name{CoefArray} with an integer index. 
\index{Coef\_Type,definition}
Associated with the {\bf enumeration} is a list of field-names, see
section~\ref{fieldnames} which are appear if a \name{Coef\_Type} is
printed. These names of the {\bf enumeration} and the list of
field-names must be kept current. The {\bf enumeration} is furthermore
used in the constructors to \name{Coefficient\_Set} and
\name{Quadruple\_Set} to initialize the required arrays.

*/
#ifndef COEF_ARRAY_H
#define COEF_ARRAY_H

#include "storage.h"
#include "typesip.h"
enum Coef_Type  { ZETASUM,PA1X,PA1Y,PA1Z,PA2X,PA2Y,PA2Z,NORM,OVERLAP,
                  ZETA1,ZETA2,XI,KINETIC,
                  PX,PY,PZ,NUCLEAR,                      // MOMENT
                  ZETASUM_PX,ZETASUM_PY,ZETASUM_PZ,      // E_FIELD
                  ZETA1_BX,ZETA1_BY,ZETA1_BZ,
		  ZETA2_AX,ZETA2_AY,ZETA2_AZ,          
		  XI_AminusBX,XI_AminusBY,XI_AminusBZ,
                  XI_AxBX,XI_AxBY,XI_AxBZ,               // ANGULAR
                  ZETAa,ZETAa_AX,ZETAa_AY,ZETAa_AZ, 
		  ZETAb,ZETAb_BX,ZETAb_BY,ZETAb_BZ,
                  ZETAab_AminusBX,ZETAab_AminusBY,ZETAab_AminusBZ,
                  ZETAab_AxBX,ZETAab_AxBY,ZETAab_AxBZ,   // SPINORB
                  ABX,ABY,ABZ,
		  RHO1,RHO2,RHOSUM,WX,WY,WZ};            // TWOEL
             

ostream& operator<<(ostream& os,Coef_Type ct);

typedef  ARRAY<Location> LocationArray;
class    Coefficient_Set;
class    RadialFunctionArray;
class    SymOp;

class CoefArray : public Storage
{
  friend class IntegArray;  
public:
~CoefArray() {  }
 CoefArray(LocationArray& l1,SymOp& s1,LocationArray& l2,SymOp& s2,
	   Coef_Type ct,int same,Allocator& alloc);
 CoefArray(RadialFunctionArray& f1,int l1,SymOp& s1,RadialFunctionArray& f2,
	   int l2, SymOp& s2, Coef_Type ct,int same,Allocator& alloc);
 CoefArray(Coefficient_Set& cs1,Coefficient_Set& cs2,Coef_Type ct,int same,
	   Allocator& alloc);
};


inline ostream& operator<<(ostream& os,const CoefArray& a)
{
  const Storage& as = a;
  return os << as; 
}
#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
inline ostream_withassign& operator<<(ostream_withassign& os,
				      const CoefArray& a)
{
  const Storage& as = a;
  return os << as; 
}
#endif
#endif



