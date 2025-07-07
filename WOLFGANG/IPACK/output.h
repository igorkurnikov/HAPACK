/*TEX
%
% IPACK - 
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
% $Id: output.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Output of Two-Electron Integrals\label{section-IO}}
\index{Input, two-electron} 
\index{Output, two-electron} 
\index{Dump, two-electron} 
\index{File, two-electron}
\index{Buffer, two-electron}
\index{Integrals, output of}
\index{Integrals, input of}

The \name{OutputBuffer} class in IPACK is used to costumize the IO of the
two electron integrals.  Presently, the \name{OutputBuffer} class
implements a VMOL compatible output file, where the two-electron
integrals are stored as blocks consisting of an array of packed
indices followed by an array of corresponding values -- see
\name{twoel\_buffer}.

The entry-point from IPACK is the \name{dump}-function which is called
with a set of four orbital indices in chemists notation and the value
of the corresponding integral, which may be zero.  These integrals are
delivered in sets of angular shells over basis sets as described in
the chapters on two-electron integrals.

\subsection{File Structure} 
The first and last \name{twoel\_buffers} on the file \name{TWOEL} are
dummies, which contain only the following important information: in
the first buffer the first element of \name{compind} reflects the
symmetry of the integrals, and in the last buffer the \name{no}$=-1$

\subsection{Class \name{OutputBuffer}}  

\begin{itemize}
%
\item \name{OutputBuffer()}: the constructor has two strings for arguments -- 
the disc-file name and IO mode ("get" or "put").
%
\item \name{dump()}: converts the four indices of an integral into a
compound index and places that and the integral value in the
two-electron integral buffer.
%
\item \name{write()}: writes the buffer to disc.
%
\item \name{read()}: reads a buffer from disc.
%
\item \name{pass\_data()}: returns a pointer to the two-electron
integral buffer.
%
\item \name{recieve\_data()}: recives and copies a pointer to a two-electron
integral buffer.
%
\item \name{rewind()}: rewinds the file, {\it and} reads the first buffer.
%
\end{itemize} */

#ifndef OUTPUT_H
#define OUTPUT_H

#include "const.h"
#include <math.h>

#include <stdio.h>
#include <stringc.h>
#include "buffer.h"


class  OutputBuffer
{
  int io_mode;
  twoel_buffer buf;
  FILE*  f;
  String name;
public: 
     OutputBuffer( const String& fname, const String& get_or_put );
    ~OutputBuffer();

void dump( unsigned int i1, unsigned int i2, unsigned int i3,
	  unsigned int i4, double val) 
{  
  if ( fabs(val)< 1E-10 ) return; 

  unsigned int i_temp;
  if ( i1 < i2 ) { i_temp = i1; i1 = i2; i2 = i_temp; }
  if ( i3 < i4 ) { i_temp = i3; i3 = i4; i4 = i_temp; }
  if ( i1 < i3 ) { i_temp = i3; i3 = i1; i1 = i_temp;
                   i_temp = i4; i4 = i2; i2 = i_temp; }
  if ( i1 == i3 && i4 > i2 ) { i_temp = i4; i4 = i2; i2 = i_temp; }
  buf.compind[buf.no] = (((((i4 << SHIFT) | i3) << SHIFT) | i2) << SHIFT) | i1;


#ifdef DEBUG
  cout << "O:" 
    <<  setw(4) << i1 << setw(4) << i2 << setw(4) << i3 <<
      setw(4) << i4 << " " << buf.compind[buf.no] << " " << val << endl;
#endif


  buf.value[buf.no++] = val;
  if ( buf.no >= BBUFSIZE ) write();

}


void write(); 
int read();

const twoel_buffer* pass_data() { return &buf; }

void rewind() { fseek( f, 0, 0 ); int no_in=read(); }
 
};

#endif
