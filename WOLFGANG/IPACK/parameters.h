/*TEX
%
% IPACK - ISAM
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
% $Id: parameters.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Parameters \label{Parameters}}
\index{Parameters} 
\index{MAX\_BUCK} 
\index{Index} 

\begin{itemize}
\item \name{LIMIT} is the size below which an integral is considered zero.
\item \name{SHIFT}, \name{SHIFT2}, \name{SHIFT3} are used in the
conversion to and from compound indices.
\item \name{BBUFSIZE} is the length of the vmol compatible buffer.
\item \name{MAX\_BUCK} defines the size of a bucket for a PIP
\item \name{Index} the type of an Index can readily be changed by
changing it definition here.  Note if changing to unsigned or types
with fewer bytes safeties should be implemented.
\end{itemize}
*/

#ifndef PARAMETERS_H
#define PARAMETERS_H

#define LIMIT     1.0e-7

#define SHIFT      8
#define SHIFT2    16
#define SHIFT3    24

#define BBUFSIZE 600

#define MAX_BUCK 100
typedef int Index;

#endif

