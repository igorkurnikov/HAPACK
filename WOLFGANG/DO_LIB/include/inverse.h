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
% $Id: inverse.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%

\section{Matrix Inversion}

Here the routines for matrix inversion and multiplications with the inverse 
of the matrix are provided. We used the algorithms given by numerical recipes. 

\subsection{LU-Decomposition}
First of all the numerical recipes routines for LU-Decomposition are 
implemented in a modified way, in order to suit the Matrix and Vector class. 

\begin{itemize}
\item{\name{void ludcmp(Matrix\& A,int n,Array<int>\& indx} is used for 
   calculating the inverse of a matrix.
   It returns in \name{A} the LU-decomposition of the input matrix.
   \name{n} is the dimension of the matrix.
   \name{indx} must be an array of the dimension \name{n} which records the 
   row permutations effected by partial pivoting.
    {\bf Note, that all offsets must be equal to 1!}
}
\item{\name{void ludcmp2(Matrix\& A,int n,Array<int>\& indx} is used for 
   calculating multiply_inverse2.
   It returns in \name{A} the LU-decomposition of the transposed input 
   matrix.
   \name{n} is the dimension of the matrix.
   \name{indx} must be an array of the dimension \name{n} which records the 
   row permutations effected by partial pivoting.
    {\bf Note, that all offsets must be equal to 1!}
}
\item{\name{void lubksb(Matrix\& A,int n,Array<int>\& indx,Array<double>\& b)}
   is used for calculating the inverse of a matrix. It solves the equation
   $ A \dot x = b $, but input for \name{A} must be the LU-decomposition given
   by ludcmp. \name{A} is not effected by this function.
   \name{n} is the dimension of \name{A}.
   \name{indx} is input as given from the above function.
   \name{b} is on input the righthandside vector $b$ of the equation and on 
   output the solution vector $x$.
    {\bf Note, that all offsets must be equal to 1!}
}
\item{\name{void lubksb2(Matrix\& A,int n,Array<int>\& indx,Array<double>\& b)}
   is used for calculating multiply_inverse2. Input for \name{A} must 
   be the LU-decomposition given
   by the ludcmp2. \name{A} is not effected by this function.
   \name{n} is the dimension of \name{A}.
   \name{indx} is input as given from the above function.
   \name{b} is on input the righthandside vector $b$ of the equation and on 
   output the solution vector $x$.
    {\bf Note, that all offsets must be equal to 1!}
}
\end{itemize}   
*/

#ifndef DOLIB__INVERSE
#define DOLIB__INVERSE

#include <assert.h>
#include "numer.h"	  

void ludcmp(Matrix& a,int n,Array<int>& indx);
void ludcmp2(Matrix& a,int n,Array<int>& indx);
void lubksb(Matrix& a,int n,Array<int>& indx,Array<double>& b);
void lubksb2(Matrix& a,int n,Array<int>& indx,Array<double>& b);

/*TEX
\subsection{Inverse}

We now define the functions for the inverse of a matrix and the multiplication
with it.
\begin{itemize}
\item{
\name{void inverse(Matrix\& A,Array<int>\& indxc,Array<int>\& indxr,
Array<int>\& ipiv)} inverts the matrix \name{A} on place using the gauss jordan 
    algorithm. Therefore three index arrays of the dimension of the matrix
    must be provided as working space.\\
    {\bf Note, that all objects must have the offsets equal to 1!}
}
\item{
  \name{inverse(Matrix\& A)} uses the above function but handles the
    index arrays and offsets on its own.
}
\item{
\name{void inverse(Matrix\& A,Matrix\& Y,Array<int>\& indx,Array<double>\& col)}
    returns the invers of \name{A} in \name{Y} using LU-decomposition, which 
    is believed to be more precise than gauss jordan. The operations count
    of both algorithms is about the same.Unfortunately \name{A} is
    destroyed on exit. Two array of the dimension of the matrix must be 
    provided as working space.\\  
    {\bf Note, that all objects must have the offsets equal to 1!}
} 
\item{
\name{void inverse(Matrix\& A,Matrix\& Y)} uses the above function 
    but handles the arrays and offsets on its own.
}
\item{
\name{void multiply\_inverse1(Matrix\& A,Matrix\& B,const Matrix\& C,
Array<int>\& indx,Array<double> col)} returns in \name{A} the product of the 
    inverse of \name{B} multiplicated with \name{C} (from the right).
    It uses LU-decomposition for the inversion, but substitutes
    directly the columns of \name{C} in the backsubstitution.
    This saves a whole matrix multiplication, and is also more precise.
    The matrix \name{B} is destroyed on exit!
    Two array of the dimension of the matrix \name{B} must be 
    provided as working space.\\  
    {\bf Note, that all offsets must be equal to 1!}
}  
\item{
\name{void multiply\_inverse1(Matrix\& A,Matrix\& B, Matrix\& C)} uses the 
     above function, but handles the arrays and offsets on its own.
}
\item{
\name{Matrix multiply\_inverse1(Matrix\& B,Matrix\& C)} is similar to the
    function given before, but returns the result of the inverse of
    \name{B} multiplied with \name{C} as a Matrix with both offsets equal
    to 1.
}
\item{
\name{void multiply\_inverse2(Matrix\& A,const Matrix\& B,Matrix\& C,
Array<int>\& indx,Array<double> col)} returns in \name{A} the product of  
    \name{B} multiplicated with the inverse of \name{C} (from the right).
    It uses LU-decomposition for the inversion, but substitutes
    directly the columns of \name{C} in the backsubstitution.
    This saves a whole matrix multiplication, and is also more precise.
    The matrix \name{C} is destroyed on exit!    
    Two array of the dimension of the matrix \name{B} must be 
    provided as working space.\\  
    {\bf Note, that all offsets must be equal to 1!}
}  
\item{
\name{void multiply\_inverse2(Matrix\& A,Matrix\& B, Matrix\& C)} uses the 
     above function, but handles the arrays and offsets on its own.
}
\item{
\name{Matrix multiply\_inverse2(Matrix\& B,Matrix\& C)} is similar to the
    function given before, but returns the result of the inverse of
    \name{B} multiplied with \name{C} as a Matrix with both offsets equal
    to 1.
}
\item{
\name{void multiply\_inverse(Vector\& a,Matrix\& B,const Vector\& c,
Array<int> indx)}
    returns in \name{a} the product of the 
    inverse of \name{B} multiplicated with \name{c}.
    It uses LU-decomposition for the inversion, but substitutes
    directly the Vector \name{c} in the backsubstitution.
    This saves a complete backsubstitution, and is also more precise.
    Two array of the dimension of the matrix \name{B} must be 
    provided as working space.\\  
    {\bf Note, that all offsets must be equal to 1!}
} 
\item{
\name{void inverse(Vector\& a,Matrix\& B,Vector\& c)} uses the above function 
    but handles the arrays and offsets on its own.
}
\item{
\name{Vector multiply\_inverse(Matrix\& B,Vector\& c)} is similar to the
    function given before, but returns the result of the inverse of
    \name{B} multiplied with \name{C} as a Vector with both offset equal
    to 1.
}
\end{itemize}
*/

void inverse(Matrix& a,Array<int>& indxc,Array<int>& indxr,Array<int>& ipiv); 
void inverse(Matrix& a);
void inverse(Matrix& a,Matrix& y,Array<int>& indx,Array<double>& col);
void inverse(Matrix& a,Matrix& y);
void multiply_inverse1(Matrix& a,Matrix& b,const Matrix& c,
                                         Array<int>& indx,Array<double> col);
void multiply_inverse1(Matrix& a,Matrix& b,Matrix& c);
Matrix multiply_inverse1(Matrix& b,Matrix& c);
void multiply_inverse2(Matrix& a,const Matrix& b,Matrix& c,
                                         Array<int>& indx,Array<double> col);
void multiply_inverse2(Matrix& a,Matrix& b,Matrix& c);
Matrix multiply_inverse2(Matrix& b,Matrix& c);
void multiply_inverse(Vector& a,Matrix& B,const Vector& c,Array<int> indx);
void multiply_inverse(Vector& a,Matrix& B,Vector& c);
Vector multiply_inverse(Matrix& B,Vector& c);

#endif

















