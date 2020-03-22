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
% $Id: const.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Constants\label{section-const}}
\index{constants, definition} 
\index{$\pi$}  
\index{$\epsilon$}  
\index{$\delta$}
\index{epsilon}
index{delta}  
This section contains the constansts relevant for this program. 

\begin{table}
\caption{Constants for IPACK}
\bc 
\begin{tabular}{ll}
\hline
pi       & $\pi$ \\
epsilon  & ERI's of magnitude smaller than $\epsilon$ are discarded. \\
delta    & two real numbers are considered equal if the differ by less 
           than $delta$. \\
bohr     & conversion factor from bohr into angstrom\\
max\_contr    & maximal degree of contraction \\
max\_radials  & total maximal number of primitive basis function of a given 
                angular momomentum \\
max\_pbf      & maximal number of pure basis functions in a pure basis set \\
max\_pbs      & maximal number of pure basis sets in a pure basis set list \\
max\_ibs      & maximal number of internal basis sets in an internal 
                basis set list \\
max\_cntr     & maximal number of centetrs in a basis set \\
max\_sym      & maximal number od \name{SymOps} and \name{SymDesignators} \\
max\_sym\_orb & maximal number of orbitals of 1 symmetry \\
\hline
\end{tabular}
\ec
\end{table} 
*/ 
#ifndef CONST_IPACK_H
#define CONST_IPACK_H
const double pi           = 3.1415926535897935;
const double delta        = 1e-13; 
const double bohr         = 0.529177249;
const int    max_momentum =    10;
const int    max_contr    =    40;
const int    max_radials  =  4000;
const int    max_pbf      =   200;
const int    max_pbs      =  2000;
const int    max_ibs      =  2000;
const int    max_cntr     =   200;
const int    max_sym      =  4096;
const int    max_sym_orb  =  2000;
const int    bbuf_size    =  50000;
const int    ibuf_size    = 100000;

// entries for FTABLE
const  int    MAXM =   41;
const  double   DT = 0.05;

/* The following two flags are defined in basis.c */

extern int spherical;    // this flag determines, whether spherical or
                         // cartesian functions are generated.
extern int angstrom;     // this flag determines, whether distances
                         // are entered in angstrom or atomic units
                         // internally ONLY atomic units are used 

enum  timer_names { t_main,t_vrr,t_hrr1,t_hrr2,t_contract,t_addto,
	 	    t_doublet, t_sp, t_f,
 	            t_2p_total = 20,t_2p_output,t_2p_disth,t_2p_distv };

void  init_timer();

#include "vtype.h"
#endif
