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
% $Id: buffer.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Two-Electron Integral Buffer\label{section-Buffers}}
\index{Buffer, two-electron}
\index{Buffer, compound index}
\index{Integrals, output of}
\index{Integrals, input of}

The two electron output file is a binary file consisting of structures
called \name{twoel\_buffers}.  Each buffer holds a maximum
\name{BBUFSIZE} integrals, which are encoded as an array of integer
compound indices and an array of values. The compound indices are
obtained as a bytewise packing of the orbital indices
$(i_1,i_2,i_3,i_4)$ in chemists notation. The entry \name{no} gives
the number of integrals present in the buffer.

The first and last \name{twoel\_buffers} on the file are dummies,
which contain only the following important information: in the first
buffer the first element of \name{compind} reflects the symmetry of
the integrals, and in the last buffer the \name{no}$=-1$

*/
  

#ifndef BUFFER_H
#define BUFFER_H

#include "parameters.h"

typedef struct 
{
  int     no;
  double  value[BBUFSIZE];
  int     compind[BBUFSIZE];
} twoel_buffer;

#endif
