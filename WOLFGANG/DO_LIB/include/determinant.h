/*IGN
%  Dortmund C++ Class and Template Library 
%  Copyright (C) 1994 Wolfgang Wenzel and others

%  The code in this file was derived from routines published in 
%  Numerical Recipes, it its the responsibility of the user to insure that
%  he is authorized to use them. Please read the appropriate section entitled
%  LICENSE in the documentation. 

%  This library is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%  Library General Public License for more details.

%  You should have received a copy of the GNU Library General Public
%  License along with this library; if not, write to the Free
%  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%
% $Id: determinant.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Determinant of a Matrix}

Here a function which calculates the determinant of a matrix is implemented. 
We used the algorithms given by numerical recipes. 

\subsection{LU-Decomposition}
First of all the numerical recipes routine for LU-Decomposition is 
implemented in a modified way, in order to suit the Matrix and Vector class. 

\begin{itemize}
\item{\name{void ludcmp(Matrix\& A,int n,double\& d)} is used for the determinat.
   It returns in \name{A} the LU-decomposition of the input matrix.
   \name{n} is the dimension of the matrix.
      \name{d} is on exit $\pm 1$ depending on whether the number of row 
   interchanges was even or odd.
    {\bf Note, that \name{a} must have the offsets equal to 1!}
}
\end{itemize}
*/

#ifndef DOLIB__DETERMINANT
#define DOLIB__DETERMINANT

#include <assert.h>
#include "numer.h"

void ludcmp(Matrix& a,int n,double& d);

/*TEX
\subsubsection{Determinant}
This is the header for the determinat function.
\begin{itemize}
\item{ \name{double determinant(Matrix\& A)} calculates the determinant
   of \name{A} using LU-decomposition. On exit the matrix \name{A} is 
   destroyed.}  
\end{itemize}
*/

double determinant(Matrix& a);

#endif



