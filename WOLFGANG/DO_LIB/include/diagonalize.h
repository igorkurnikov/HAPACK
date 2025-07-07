/*TEX
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

\section{Diagonalize, Orthonormalize and Singular Value Decomposition}

\begin{itemize}
\item{
\name{void diagonalize( Matrix\& mat, Matrix\& eigen, Vector\& eval )}
Diagonalizes a {\bf symmetric} matrix, returning the eigenvectors and
eigenvalues.  The eigenvalues are sorted in ascending order, and the
corresponding eigenvectors are orthonormal.  Machine specific flags
can be used, for example see "SUN", which is the only
non-package diagonalizer presently implemented.  The regular diagonalizer
is from \name{Numerical Recipes}.
}
%
\item{
\name{void orthonormalize( Matrix\& S, Matrix\& u )} 
Determines $u=S^{-1/2}$, by determining the eigenvalues and
eigenvectors of $S$, with \name{diagonalize()}, and then computing $u$.
}
%
\item{
\name{void geneva( NumArray2D<double>\& mat, NumArray2D<double>\& S,
	     NumArray2D<double>\& evec, NumArray<double>\& eval )}
Solves the generalized eigenvalue and vector problem $mat\phi = ES\phi$.
}
%
\item{
\name{void svd( Matrix\& mat, Matrix\& eigen, Vector\& eval )}
Diagonalizes a matrix using a portable (LAPACK based singular value 
decomposition.
}
%
\item{
\name{void LR_svd( Matrix\& mat, Matrix\& U, Matrix\& V, Vector\& eval )}
Performs a singular value decomposion of the matrix \name{mat}, returning
the left and right handed eigenvectors and the eigenvalues. 
}
%
\item{
\name{void straightsort( Vector\& eval, Matrix\& eigen )}
Sorts \name{eval} in ascending order, rearranging \name{eigen} at
the same time.
}
\end{itemize}

These functions diagonalizes or orthonormalize a {\bf symmetric} matrix. */

#ifndef DOLIB__DIAGONALIZE
#define DOLIB__DIAGONALIZE

#include "numer.h"

enum DIAG_ERROR_CODES_ { DIAG_OK, DIAG_DEV, DIAG_TQLI_ITER, DIAG_EMPTY };

int  diagonalize(Matrix& mat, Matrix& eigen, Vector& eval);
void orthonormalize( Matrix& S, Matrix& u );

void geneva( NumArray2D<double>& mat, NumArray2D<double>& S,
	     NumArray2D<double>& evec, NumArray<double>& eval );

void svd( Matrix& mat, Matrix& eigen, Vector& eval );
void LR_svd( const Matrix& mat, Matrix& U, Matrix& V, Vector& eval );
void LR_svd( const NumArray2D<double>& mat, NumArray2D<double>& U,
	     NumArray2D<double>& V, NumArray<double>& eval );

//  For the portable, Numerical-Recipes-based diagonalizer.


void straightsort( Vector& eval, Matrix& eigen );
int  tqli(Vector& d, Vector& e, Matrix& z );
void tred2old( Matrix& a, Vector& d, Vector& e );

#endif

