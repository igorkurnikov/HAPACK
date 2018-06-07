/*TEX
%  Dortmund C++ Class and Template Library 
%  Copyright (C) 1994 Wolfgang Wenzel

%  This library is free software; you can redistribute it and/or
%  modify it under the terms of the GNU Library General Public
%  License as published by the Free Software Foundation; either
%  version 2 of the License, or (at your option) any later version.

%  This library is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%  Library General Public License for more details.

%  You should have received a copy of the GNU Library General Public
%  License along with this library; if not, write to the Free
%  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
  
%  $Id: numer.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $

\chapter{Definition of Vectors and Matrices in  DoLib}

\section{Definition of Double Matrices}

This class is derived from the template class NumArray2D
It enables the simple matrix arithmetic defined in NumArray. 
Further given are:
\begin{itemize}
\item{ \name{mat1.multiply(mat2,mat3)} multiplies \name{mat2}
	  with \name{mat3} and adds the resulting matrix to \name{mat1}.
	  The operation fails if the dimensions do
	  not match. Offsets are ignored. An operator version of the
	  operation is defined, its use is discouraged.}
\item{ \name{mat1.transpose()} transposes {\bf SQUARE} matrices
          leaving the offsets unchanged. It fail for non-square matrices.}
\item{ \name{mat1.multiply\_transpose1(mat2,mat3)} multiplies the transposed
          matrix \name{mat2} (not neccessary square)
	  with \name{mat3} and adds the resulting matrix to \name{mat1}.
	  The operation fails if the dimensions do not match. 
          Offsets are ignored.}
\item{ \name{mat1.multiply\_transpose2(mat2,mat3)} multiplies \name{mat2}
	  with the transposed matrix \name{mat3} (not neccessary square) 
          and adds the resulting matrix to \name{mat1}.
	  The operation fails if the dimensions do
	  not match. Offsets are ignored.}
\end{itemize}
*/  

#ifndef DOLIB__NUMER
#define DOLIB__NUMER

#include "tnumarray.h"
#include <math.h>

typedef Array2D<int> IMatrix;

class Matrix : public NumArray2D<double>
{
public:
       Matrix() { } 
       Matrix(const int d1,const int d2,const int o1=0,const int o2=0) :
         NumArray2D<double>(d1,d2,o1,o2) { }
       Matrix(const NumArray2D<double> d) : NumArray2D<double>(d) { }
       Matrix(const Matrix& d) : NumArray2D<double>(d) { }
virtual ~Matrix() { } 

Matrix& operator=(const Matrix& m) { NumArray2D<double>::operator=(m);
				     return *this;}

#ifdef BLAS
void multiply( const Matrix& b, const Matrix& c );
void multiply_transpose1( const Matrix& b, const Matrix& c );
void multiply_transpose2( const Matrix& b, const Matrix& c );
#endif

};

inline Matrix operator+(const Matrix& mat1,const Matrix& mat2)
{
  const NumArray2D<double>& m1 = mat1;  
  const NumArray2D<double>& m2 = mat2;
  NumArray2D<double> result(m1);
  result += m2;
  return Matrix(result);
}


inline Matrix operator-(const Matrix& mat1,const Matrix& mat2)
{
  const NumArray2D<double>& m1 = mat1;  
  const NumArray2D<double>& m2 = mat2;
  NumArray2D<double> result(m1);
  result -= m2;
  return Matrix(result);
}

inline Matrix operator*(const Matrix& mat1,const Matrix& mat2)
{
  const NumArray2D<double>& m1 = mat1;  
  const NumArray2D<double>& m2 = mat2;
  NumArray2D<double> result( m1.size1(), m2.size2() );
  result.set(0);
  result.multiply(m1,m2);
  return Matrix(result);
}


inline ostream& operator<<(ostream& os,const Matrix& m)
{
  const NumArray2D<double>& na = m;
  return os << na;
}

#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
inline ostream_withassign& operator<<(ostream_withassign& os,const Matrix& m)
{
  const NumArray2D<double>& na = m;
  return os << na;
}
#endif
/*TEX
\section{Matrices with Index Checking.} 
*/

class D_Matrix : public Matrix
{
public:
       D_Matrix() { } 
       D_Matrix(const int d1,const int d2,const int o1=0,const int o2=0) :
         Matrix(d1,d2,o1,o2) { }
       D_Matrix(const NumArray2D<double> d) : Matrix(d) { }
       D_Matrix(const Matrix& d) : Matrix(d) { }
virtual ~D_Matrix() { } 

D_Matrix& operator=(const D_Matrix& m) { Matrix::operator=(m);
                                         return *this;}

double& operator()(const int p1,const int p2) 
        { 
	  if(p1 < off1 || p1 >= dim1+off1 || p2 < off2 || p2 >= dim2+off2)
          error("Acces Violation",p1,p2);
	  return pbase[p2][p1];
	}
const double& get(const int p1,const int p2) const
        { 
	  if(p1 < off1 || p1 >= dim1+off1 || p2 < off2 || p2 >= dim2+off2)
          error("Acces Violation",p1,p2);
	  return pbase[p2][p1];
	}
};
	  
inline ostream& operator<<(ostream& os,const D_Matrix& array)
{
  const NumArray2D<double>& a = array; 
  return os << a;
}

#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
inline ostream_withassign& 
operator<<(ostream_withassign& os,const D_Matrix& array)
{
  const NumArray2D<double>& a = array; 
  return os << a;
}
#endif

/*TEX
\section{Definition Integer Vectors}

This class is derived from the template class NumArray. It
enables the simple vector arithmetic and ordering operations as 
defined in NumArray. The scalar product is given, too, but no multiplication
with a matrix.

*/
class IVector : public NumArray<int>
{
public:
         IVector() { }
         IVector(const int new_len, const int new_off = 0) :
	   NumArray<int>(new_len,new_off) { }
         IVector(const NumArray<int>& vec) : NumArray<int>(vec) { }
         IVector(const IVector& vec)  : NumArray<int>(vec) { }
virtual ~IVector () { }

IVector& operator=(const IVector& v) { NumArray<int>::operator=(v);
					return *this; }
}; 

inline IVector operator+(const IVector& vec1,const IVector& vec2)
{
  const NumArray<int>& v1 = vec1;  
  const NumArray<int>& v2 = vec2;
  return IVector(v1+v2);
}

inline IVector operator-(const IVector& vec1,const IVector& vec2)
{
  const NumArray<int>& v1 = vec1;  
  const NumArray<int>& v2 = vec2;
  return IVector(v1-v2);
}

inline ostream& operator<<(ostream& os,const IVector& array)
{
  const NumArray<int>& a = array; 
  return os << a;
}

#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
inline ostream_withassign& 
operator<<(ostream_withassign& os,const IVector& array)
{
  const NumArray<int>& a = array; 
  return os << a;
}
#endif
/*TEX
\section{Definition of Integer Vectors with Index Checking} 
*/

class D_IVector: public IVector
{
public:  
         D_IVector() { }
         D_IVector(const int new_len, const int new_off = 0) :
	   IVector(new_len,new_off) { }
         D_IVector(const NumArray<int>& vec) : IVector(vec)   { } 
         D_IVector(const IVector& vec)  : IVector(vec) { }
virtual ~D_IVector () { }
     
D_IVector& operator=(const D_IVector& iv) { IVector::operator=(iv);
                                         return *this;}

int&     operator[](const int n) 
         {
	   if (n<off || n >= len+off) error("Illegal Index",n);
	   return pbase[n];
	 }
const int& operator()(const int n) const
         {
	   if (n<off || n >= len+off) error("Illegal Index",n);
	   return pbase[n];
	 }

};


inline ostream& operator<<(ostream& os,const D_IVector& array)
{
  const NumArray<int>& a = array; 
  return os << a;
}

#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
inline ostream_withassign& 
operator<<(ostream_withassign& os,const D_IVector& array)
{
  const NumArray<int>& a = array; 
  return os << a;
}
#endif     
     
/*TEX
\section{Definition of Double Vectors} 

This class is derived from the template class NumArray It enables the
simple vector arithmetic and ordering operations as defined in
NumArray. The scalar product and multiplication with a matrix is
given, too.
*/

class Vector : public NumArray<double> {

public:
         Vector() { }
         Vector(const int new_len, const int new_off = 0) :
	   NumArray<double>(new_len,new_off) { }
         Vector(const NumArray<double>& vec)  : NumArray<double>(vec) { }
         Vector(const Vector& vec)  : NumArray<double>(vec) { }
virtual ~Vector () { }

#ifndef BLAS
Vector& operator=(const Vector& v) { NumArray<double>::operator=(v); return *this; }
#else
Vector& operator=(const Vector& v);
Vector& operator+=(const Vector& v);
double  operator*(const Vector& v);
void    multiply(const Matrix& mat, const Vector& vec);
void    set(const double value);
Vector& operator+=(const double value);
#endif

};

inline ostream& operator<<(ostream& os,const Vector& array)
{
  const NumArray<double>& a = array;
  return os << a;
}

#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
inline ostream_withassign& operator<<(ostream_withassign& os,
				      const Vector& array)
{
  const NumArray<double>& a = array;
  return os << a;
}
#endif

#ifndef BLAS

inline Vector operator+(const Vector& vec1,const Vector& vec2) 
{
  const NumArray<double>& v1 = vec1;  
  const NumArray<double>& v2 = vec2;
  NumArray<double> result(v1);
  result += v2;
  return Vector(result);
}

inline Vector operator-(const Vector& vec1,const Vector& vec2) 
{
  const NumArray<double>& v1 = vec1;  
  const NumArray<double>& v2 = vec2;
  NumArray<double> result(v1);
  result -= v2;
  return Vector(result);
}

#endif

inline Vector operator*(const Matrix& mat,const Vector& vec) {
  Vector y( mat.size1() );
  y.set(0.0);
  y.multiply( mat, vec );
  return y;
}

/*

The following function calls the operator NumArray<Base>
operator*(const NumArray<Base>& v,const NumArray2D<Base>& mat) which
is implemented in tnumarray.h but could not be resolved by the GNU
compiler. Therefore it is commented out.
 
The IBM compiler xlC has simple type conversion for templates and does not 
need this function anyway.


inline Vector operator*(const Vector& vec,const Matrix& mat)
{
  const NumArray<double>& v = vec;  
  const NumArray2D<double>& m = mat;
  return Vector(v*m);
} 
*/

/*TEX
\section{Definition of Double Vectors with Index Checking} 
*/

class D_Vector: public Vector
{
public:  
         D_Vector() { }
         D_Vector(const int new_len, const int new_off = 0) :
	   Vector(new_len,new_off) { }
         D_Vector(const NumArray<double>& vec) : Vector(vec)   { } 
         D_Vector(const Vector& vec)  : Vector(vec) { }
virtual ~D_Vector () { }
     
D_Vector& operator=(const D_Vector& v) { Vector::operator=(v);
                                         return *this;}

double&  operator[](const int n) 
         {
	   if (n<off || n >= len+off) error("Illegal Index",n);
	   return pbase[n];
	 }
const double& operator()(const int n) const
         {
	   if (n<off || n >= len+off) error("Illegal Index",n);
	   return pbase[n];
	 }
};

inline ostream& operator<<(ostream& os,const D_Vector& array)
{
  const NumArray<double>& a = array;
  return os << a;
}

#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
inline ostream_withassign& operator<<(ostream_withassign& os,const D_Vector& array)
{
  const NumArray<double>& a = array;
  return os << a;
}
#endif

#ifdef BLAS
#include "blas.h"
#endif

void svdcmp(Matrix& a,Vector& w,Matrix& v);

#endif




