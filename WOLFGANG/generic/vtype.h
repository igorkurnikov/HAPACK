/*TEX
%
% QCHEM:GENERIC
% Quantum Chemistry Project: Density Matrix Renormalization Group
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is proprietary software. Unlicensed use and distribution
% of this program or parts thereof are illegal. For details on the 
% license, see the LICENSE section in the file "dmrg.tex" distributed 
% with this package or write to: 
%
%      Wolfgang Wenzel, Theoretical Physics I, 
%      University of Dortmund,
%      D-44221 Dortmund, Germany 
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the LICENSE
% for more details.
%
% $Id: vtype.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $

\subsection{Array Type and Debugging Level Definitions}
*/
#include <numer.h>

class Vector;
class Matrix;
class IVector;

class D_Vector;
class D_Matrix;
class D_IVector;
#include <xtime.h>

/*TEX
The following define determines whether \name{Arrays} or \name{D\_Arrays}
are being used throughout the code.
*/
typedef double        IntegType;
#define MPI_INTEGTYPE MPI_DOUBLE 
typedef float         QCI_INTEGTYPE;
#define MPI_QCI_INTEGTYPE MPI_FLOAT;

#ifdef DEBUG_ARRAYS
typedef D_IVector   IVec;
typedef D_Vector    Vec;
typedef D_Matrix    Mat;
#define ARRAY       D_Array
#define ARRAY2D     D_Array2D
#define NUMARRAY    D_NumArray
#define NUMARRAY2D  D_NumArray2D
#else  
typedef IVector     IVec;
typedef Vector      Vec;
typedef Matrix      Mat;
#define ARRAY       Array
#define ARRAY2D     Array2D
#define NUMARRAY    NumArray
#define NUMARRAY2D  NumArray2D
#endif

#define IntegVec Vec
#define IntegMat Mat

extern Timer* timer;

#ifdef CRAY
#define INT_4 short int
#else
#define INT_4 int
#endif

typedef unsigned char Orbtype;



