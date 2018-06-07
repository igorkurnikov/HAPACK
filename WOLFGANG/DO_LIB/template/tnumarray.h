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
%  $Id: tnumarray.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%

\chapter{Numeric Arrays}

In the following we define a set of templates for numeric arrays,
which are derived from the standard arrays \name{Array} and
\name{Array2D}. Base types to these templates are required to define
the full set of arithmetic operations and an ordering relation.

\section{One Dimensional Arrays - Numarray1D}

\subsection{Vector Arithmetic}

Arithmetic on vectors and matrices is plagued by the problem of temporaries.
In the expression
$$ A = C + \alpha*B;$$
where $A,B$ and $C$ are vectors, while $\alpha$ is real, two temporaries
are created during evaluation. The slightly less elegant implementation
$$ A = B; A*=\alpha; A+=C; $$
uses none. While we provide operators for both methods of evaluation the
first is strongly discouraged.

We have presently overloaded multiplication, addition, subtraction and
division by a constant. We have implemented a scalar product.

\subsection{Sorting}

Simple sorting routines have been provided. The functions \name{maxpos}
and \name{minpos} return the index of the maximum/minimum value of the 
vector. The functions \name{max} and \name{min} return the corresponding
values. The function \name{sort} implements a simply maximum-sort for
the vector elements and could stand improvement.

\subsection{Description of Functions}

\begin{itemize}
\item{The following vector operations are implemented:
   \begin{itemize}
   \item{Scalar product: \name{vec1*vec2} returns the scalar 
         product of \name{vec1} and \name{vec2}. The operation
         fails if the length of the objects differ. The offsets
         are ignored. }
   \item{Scaling: \name{vec *= factor;} multiplies all elements with
         the scalar \name{factor}.}
   \item{Scaling: \name{vec /= factor;} devides all elements by
         the scalar \name{factor}.}
   \item{\name{vec += value;} adds the scalar \name{value}
         to all elements.} 
    \item{\name{vec -= value;} subtracts the scalar \name{value}
         from all elements.}
    \item{\name{vec += vec2;} adds each element of \name{vec2}
          to the corresponding element of \name{vec}. The operation
         fails if the length of the objects differ. The offsets
         are ignored. }
    \item{\name{vec -= vec2;} subtracts each element of \name{vec2}
         from the corresponding element of \name{vec}. The operation
         fails if the length of the objects differ. The offsets
         are ignored. }
    \item{The element-by-element operations for addition and subtraction
         have been defined. Their use is discouraged. These operations
         fail if the length of the objects differ. The offsets
         are ignored. }       
   \end{itemize}
   Aditionally the following operations are implemented if a matrix
   class of the same type is defined.
   \begin{itemize}
   \item{ \name{vec1.mult(matrix,vec2);} multiplies \name{matrix}
          with \name{vec2} and adds the resulting vector to \name{vec1}.
          The operation fails if the vector and matrix dimensions do
          not match. Offsets are ignored. An operator version of the
          operation is defined, its use is discouraged.}
   \item{ \name{vec1.mult(vec2,matrix);} multiplies \name{vec2}
          with \name{matrix} and adds the resulting vector to \name{vec1}.
          The operation fails if the vector and matrix dimensions do
          not match. Offsets are ignored. An operator version of the
          operation is defined, its use is discouraged.}
   \end{itemize}
   }
\item{The \name{vec.sort(int ascending=1)} function sorts the elements of the
vector in ascending/descending order if the parameter is true/false. The 
default order is ascending. Related functions are:
   \begin{itemize}
   \item{\name{vec.maxpos()} returns the position of the maximum element,
         it fails for an empty vector. }
   \item{\name{vec.minpos()} returns the position of the minimum element,
         it fails for an empty vector. }         
   \item{\name{vec.max()} returns the value of the maximum element,
         it fails for an empty vector. }
   \item{\name{vec.min()} returns the value of the minimum element,
         it fails for an empty vector. }
   \end{itemize}
}
\item{The \name{vec1.gather(vec2,index)} performs the operation:
      $$ vec1[i] = vec2[index[i]]; $$
      for all $i$. The length of \name{index} must exceed the length 
      of $vec1$. The user must assure that all elements in \name{index}
      are valid indices in \name{vec2}. }
\item{The \name{vec1.scatter(vec2,index)} performs the operation:
      $$ vec1[index[i]] = vec2[i]; $$
      for all $i$. The length of \name{index} must exceed the length 
      of $vec2$. The user must assure that all elements in \name{index}
      are valid indices in \name{vec1}, no checking is performed.
      The functions require an index array of type Array<int> e.g.
      IVector. }
\end{itemize}

*/
#ifndef TNUMARRAY_H
#define TNUMARRAY_H
#include <tarray.h>

typedef Array<int> IndexArray;

template <class Base> class NumArray2D;

#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif


template <class Base> class NumArray : public Array<Base>
{
void vecinit() { } 
public:
         NumArray() { vecinit(); }
         NumArray(const int new_len, const int new_off = 0) :
	   Array<Base>(new_len,new_off) { vecinit(); }
         NumArray(const Array<Base>& vec) : Array<Base>(vec) { }
virtual ~NumArray () { }
//compGcc4
using Array<Base>::pbase;
using Array<Base>::base;
using Array<Base>::len;
using Array<Base>::off;
using Array<Base>::error;

void     sort(const int ascending =1);
int      maxpos() const;
int      minpos() const;
Base     max() const { return pbase[maxpos()]; }
Base     min() const { return pbase[minpos()]; }
void     zero() const { memset(base,0,sizeof(Base)*len);} 
void     gather(const NumArray<Base>& source,const IndexArray& index);
void     scatter(NumArray<Base>& source,const IndexArray& index);

NumArray<Base>& operator=(const NumArray<Base>& vec)
                   { Array<Base>::operator=(vec);
                     return *this; }

Base     operator*(const NumArray<Base>& v) const;

NumArray<Base>& operator*=(const Base value);
NumArray<Base>& operator/=(const Base value);
NumArray<Base>& operator+=(const Base value);
NumArray<Base>& operator-=(const Base value);

NumArray<Base>& operator+=(const NumArray<Base>& v);
NumArray<Base>& operator-=(const NumArray<Base>& v);

#if defined(GCC_NEW) || defined(_MSC_VER) || defined(__DECCXX)
friend NumArray<Base> operator+ (const NumArray<Base>& v,
				const NumArray<Base>& vec2);
friend NumArray<Base> operator- (const NumArray<Base>& v,
				const NumArray<Base>& vec2);

friend NumArray<Base> operator* (const NumArray2D<Base>& mat,
				const NumArray<Base>& v);
friend NumArray<Base> operator* (const NumArray<Base>& v,
				const NumArray2D<Base>& mat);
#else
friend NumArray<Base> operator+ <Base>(const NumArray<Base>& v,
				const NumArray<Base>& vec2);
friend NumArray<Base> operator- <Base>(const NumArray<Base>& v,
				const NumArray<Base>& vec2);

friend NumArray<Base> operator* <Base>(const NumArray2D<Base>& mat,
				const NumArray<Base>& v);
friend NumArray<Base> operator* <Base>(const NumArray<Base>& v,
				const NumArray2D<Base>& mat);
#endif
void multiply(const NumArray  <Base>& vec, const NumArray2D<Base>& mat);
void multiply(const NumArray2D<Base>& mat, const NumArray  <Base>& vec);
void scaled_add(Base factor,const NumArray<Base>& vec);

void print(ostream& os,int elem_per_line=10,char* format = 0) const;
};

template <class Base>
ostream& operator<<(ostream& os,const NumArray<Base>& array)
{
  operator<<(os,(const Array<Base>&) array);
  return os;
}

#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
template <class Base>
ostream_withassign& operator<<(ostream_withassign& os,const NumArray<Base>& array)
{
  operator<<(os,(const Array<Base>&) array);
  return os;
}
#endif

/*TEX
\section{Implemetation of the NumArray Template}

\subsubsection{Vector Arithmetic}
*/

template <class Base>
Base  NumArray<Base>::operator*(const NumArray<Base>& v) const 
{
  if (v.len != len)
    error("Opertor* --NumArray Length does not match. ");
  Base sum = 0;
  int i;
  for(i=0; i<len;i++)
    sum += base[i]*v.base[i];
  return sum;
}


template <class Base>
NumArray<Base>& NumArray<Base>::operator*=(const Base value)
{
  int i;
  for(i=0; i<len;i++)
    base[i] *= value;
  return *this;
}

template <class Base>
NumArray<Base>& NumArray<Base>::operator/=(const Base value)
{
  for(int i=0; i<len;i++)
    base[i] /= value;
  return *this;
}

template <class Base>
NumArray<Base>& NumArray<Base>::operator+=(const NumArray<Base>& v)
{
  if (v.len != len)
    error("Opertor* --NumArray Length does not match. ");
  int i;
  for(i=0; i<len;i++)
    base[i] += v.base[i];
  return *this;
}

template <class Base>
NumArray<Base>& NumArray<Base>::operator+=(const Base value)
{
  for(int i=0; i<len;i++)
    base[i] += value;
  return *this;
}

template <class Base>
NumArray<Base>& NumArray<Base>::operator-=(const NumArray<Base>& vec)
{
  if (vec.len != len)
    error("Opertor* --NumArray Length does not match. ");
  int i;
  for(i=0; i<len;i++)
    base[i] -= vec.base[i];
  return *this;
}

template <class Base>
NumArray<Base>& NumArray<Base>::operator-=(const Base value)
{
  int i;
  for(i=0; i<len;i++)
    base[i] -= value;
   return *this;
}

template <class Base>
void NumArray<Base>::multiply(const NumArray2D<Base>& mat,
			      const NumArray<Base>& vec)
{
  int i,j;
  if(len < mat.dim1)
    error("MULTIPLY: Output Vector Length not compatible with Matrix.");
  if (vec.len != mat.dim2)
    error("MULTIPLY: Input Vector Length not compatible with Matrix.");
  
  /*
  for(j=0; j < vec.len;j++)
  {
    Base* mbase = mat.pbase[j+mat.off2] + mat.off1;
    Base  vfac  = vec.base[j];
    for(i=0; i < mat.dim1;i++)
      base[i] += vfac*mbase[i];
  }
  */
  for(i=0; i < mat.dim1;i++)
  {
    Base  sum = 0;
    Base* mbase = mat.base + i;
    int   pos = 0;
    int   skip = mat.dim1;
    Base* vbase = vec.base;
    
    for(j=0; j < vec.len;j++)
    {
      sum += vbase[j]* mbase[pos];
      pos += skip;
    }
    base[i] += sum;
  }
}

template <class Base>
void NumArray<Base>::multiply(const NumArray<Base>& vec,
			      const NumArray2D<Base>& mat)
{
   int i,j;
   if(len < mat.dim2)
     error("MULTIPLY: Vector Length must match 2nd dimension of matrix.");
   if (vec.len != mat.dim1)
     error("MULTIPLY: Vector Length must match 1st dimension of matrix.");

   for(j=0; j < mat.dim2;j++)
     for(i=0; i < vec.len;i++)
       base[j] += mat.pbase[j+mat.off2][i+mat.off1] * vec.base[i];
 }

template <class Base>
void NumArray<Base>::scaled_add(Base factor,const NumArray<Base>& vec)
{
   int i,j;
   if (len < vec.len)
     error("scaled_add: Input Vector Length not compatible.");
   
   for(i=0; i < vec.len;i++)
     base[i] += factor*vec.base[i];
}

template <class Base>
NumArray<Base> operator*(const NumArray<Base>& v,const NumArray2D<Base>& mat)
{
  NumArray<Base> product(mat.size2());
  product.set(0);
  product.multiply(v,mat);
  return product;
}

template <class Base>
NumArray<Base> operator*(const NumArray2D<Base>& mat,const NumArray<Base>& v)
{
  NumArray<Base> product(mat.size1());
  product.set(0); 
  product.multiply(mat,v);
  return product;
}

template <class Base>
NumArray<Base> operator+(const NumArray<Base>& v,const NumArray<Base>& vec2)
{
  NumArray<Base> result(v);
  result += vec2;
  return result;
}

template <class Base>
NumArray<Base> operator-(const NumArray<Base>& v,const NumArray<Base>& vec2)
{
  NumArray<Base> result(v);
  result -= vec2;
  return result;
}

/*TEX
\subsection{Sorting}
*/

template <class Base>
void NumArray<Base>::sort(const int ascending)
{
  if (len <=1) return;
  int pos = 0;
  for(pos = 0; pos < len-1; pos++)
  {
    int p = pos;
    if (ascending)
    {
      for(int i=pos+1; i < len;i++)
	if (base[i] < base[p])
	  p = i;
    }
    else
    {
      for(int i=pos+1; i < len;i++)
	if (base[i] > base[p])
	  p = i;
    }
    Base tmp  = base[pos];
    base[pos] = base[p];
    base[p]   = tmp;
  }
}

template <class Base>
int  NumArray<Base>::maxpos() const
{
  if (len == 0)
    error("maxpos called for vector of length zero. ");
  int p = 0;
  int i;
  for(i=0; i < len;i++)
    if (base[i] > base[p]) p = i;
  return p+off;
}

template <class Base>
int  NumArray<Base>::minpos() const
{
  if (len == 0)
    error("minpos called for vector of length zero. ");
  int p = 0;
  int i;
  for(i=0; i < len;i++)
    if (base[i] < base[p]) p = i;
  return p+off;
}


template <class Base>
void  NumArray<Base>::gather(const NumArray<Base>& source,
			     const Array<int>& idx)
{
  int i;
  for(i=0; i<idx.size(); i++)
    pbase[i] = source.pbase[idx(i)];
}


template <class Base>
void NumArray<Base>::scatter(NumArray<Base>& source,const Array<int>& index)
{
  int i;
  for(i=index.offset(); i<index.size()+index.offset(); i++)
    pbase[index(i)] = source.pbase[i];
}
/*TEX
\subsection{PrettyPrinting}
*/
template <class Base>
void NumArray<Base>::print(ostream& os,int elem_per_line,char* format) const
{

  if (len == 0)
  { 
    os << "Vector is empty. " << endl;
    return;
  }  
  os << "  *** Vector *** Lo: "   << off << " *** Hi: " << len+off-1 << endl;
  int i;
  char buf[100];
  for(i=0; i < len; i++)
  {
    if (i % elem_per_line == 0) 
    { 
      int max = i + elem_per_line - 1; 				   
      if (max > len-1) max = len - 1;  		   
      sprintf(buf,"\n%4i..%4i *",off+i,off+max);
      os << buf;
    }			
    if (format)
    {      
      sprintf(buf,format,base[i]);
      os << buf;   			   
    }
    else 
      os << base[i] << " ";
  }
  os << endl;
}
 
/*TEX
\section{Template for NumArray with Access Checking}
*/

template <class Base> class D_NumArray : public NumArray<Base>
{
public:
         D_NumArray() { }
         D_NumArray(const int new_len, const int new_off = 0) :
	   NumArray<Base>(new_len,new_off) {}
         D_NumArray(const NumArray<Base>& vec) : NumArray<Base>(vec) { }
         D_NumArray(const D_NumArray<Base>& vec)  : NumArray<Base>(vec) { }
virtual ~D_NumArray () { }
//compGcc4
using NumArray<Base>::pbase;
using NumArray<Base>::len;
using NumArray<Base>::off;
using NumArray<Base>::error;

D_NumArray<Base>& operator=(const D_NumArray<Base>& vec)
                   { NumArray<Base>::operator=(vec);
                     return *this; }

Base&    operator[](const int n)
{
  if (n<off || n >= len+off) error("Illegal Index",n);
  return pbase[n];
}
const Base& operator()(const int n) const
{
  if (n<off || n >= len+off) error("Illegal Index",n);
  return pbase[n];
}
};

/*TEX
\section{Two Dimensional Numeric Arrays -- NumArray2D } 

\subsection{Matrix Arithmetic}

The arithmetic operations were implemented in the same way as
for the one dimesional array . We remind to the problem of temporaries.

Presently overloaded are the operators $+=$ and $-=$,
adding/substracting either a scalar value or performing an
element by element addition/subtraction between two matrices.
In the latter case, the dimensions of the matrix must match.
Multiplication and divison $*=$ and $/=$ with a scalar have
been overloaded, too.
Operator versions for these operations have also been
implemented,their use is discouraged.

We implemented matrix-multiplication in the member function \name{multiply},
which adds the result of the product to the final matrix:
$$ A.multiply(B,C)~~\rightarrow~~A += B * C: $$
Offsets are ignored.
The operation fails if the dimensions of the matrices do not match.
An operator version $$ A = B * C $$ of the matrix multiplication 
(ignoring offsets) exists.

Its use is discouraged.

\subsection{Further Operations}

At least matrix-operations like \name{transpose}, \name{inverse},
\name{diagonalize}, \name{determinant} and \name{norms} are planned
to be provided.

Available are yet:
\begin{itemize}
\item{\name{transpose} Is only implemented for square matrices in order
      to avoid a temporary matrix. Therefore only the elements are exchanged
      and the offsets are kept (no reset).
\item{multiply\_transpose} Is similar to A.mulitply(B,C) but transposes
      either B or C.}
\end{itemize}

\subsection{Descriptions of Functions}

The following matrix operations are implemented:
\begin{itemize}
\item{Scaling: \name{mat *= factor;} multiplies all elements with
      the scalar \name{factor}.}
\item{Scaling: \name{mat /= factor;} devides all elements by
         the scalar \name{factor}.}
\item{\name{mat += value;} adds the scalar \name{value}
         to all elements.} 
\item{\name{mat -= value;} subtracts the scalar \name{value}
         from all elements.}
\item{\name{mat1 += mat2;} adds each element of \name{mat2}
          to the corresponding element of \name{mat1}. The operation
         fails if the dimensions of the objects differ. The offsets
         are ignored. }
\item{\name{mat1 -= mat2;} subtracts each element of \name{mat2}
         from the corresponding element of \name{mat1}. The operation
         fails if the dimensions of the objects differ. The offsets
         are ignored. }
\item{The element-by-element operations for the addition and
	 subtraction of matrices
	 have been defined. Their use is discouraged. These operations
	 fail if the dimensions of the objects differ. The offsets
	 are ignored. }
   \end{itemize}

Aditionally the following operations are implemented:
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
	    Offsets are ignored. }

\item{ \name{mat1.multiply\_transpose2(mat2,mat3)} multiplies \name{mat2}
	  with the transposed matrix \name{mat3} (not neccessary square) 
          and adds the resulting matrix to \name{mat1}.
	  The operation fails if the dimensions do
	  not match. Offsets are ignored.} 
\end{itemize}

*/
template <class Base> class NumArray2D : public Array2D<Base>
{
  friend class NumArray<Base>;
public:  
       NumArray2D() { } 
       NumArray2D(const NumArray2D<Base>& src) : Array2D<Base>(src) {}
       NumArray2D(const int d1,const int d2,const int o1=0,const int o2=0) :
         Array2D<Base>(d1,d2,o1,o2) { }
virtual ~NumArray2D() { } 

//compGcc4
using Array2D<Base>::dim1;
using Array2D<Base>::dim2;
using Array2D<Base>::off1;
using Array2D<Base>::off2;
using Array2D<Base>::base;
using Array2D<Base>::pbase;
using Array2D<Base>::error;

#if defined(GCC_NEW) || (_MSC_VER)
friend   NumArray2D<Base> operator* (const NumArray2D<Base>& a,
				    const NumArray2D<Base>& b);
friend   NumArray2D<Base> operator+ (const NumArray2D<Base>& mat1,
				 const NumArray2D<Base>& mat2);
friend   NumArray2D<Base> operator- (const NumArray2D<Base>& mat1,
				 const NumArray2D<Base>& mat2);
#else
friend   NumArray2D<Base> operator* <Base>(const NumArray2D<Base>& a,
				    const NumArray2D<Base>& b);
friend   NumArray2D<Base> operator+ <Base>(const NumArray2D<Base>& mat1,
				 const NumArray2D<Base>& mat2);
friend   NumArray2D<Base> operator- <Base>(const NumArray2D<Base>& mat1,
				 const NumArray2D<Base>& mat2);
#endif
void   multiply(const NumArray2D<Base>& b,const NumArray2D<Base>& c);
void   multiply_transpose1(const NumArray2D<Base>& b,const NumArray2D<Base>& c);
void   multiply_transpose2(const NumArray2D<Base>& b,const NumArray2D<Base>& c);
void   transpose();

NumArray2D<Base>& operator=(const NumArray2D<Base>& mat)
                   { Array2D<Base>::operator=(mat);
                     return *this; }

NumArray2D<Base>& operator+=(const NumArray2D<Base>& mat);
NumArray2D<Base>& operator-=(const NumArray2D<Base>& mat);
NumArray2D<Base>& operator+=(const Base value);
NumArray2D<Base>& operator-=(const Base value);
NumArray2D<Base>& operator*=(const Base value);
NumArray2D<Base>& operator/=(const Base value);
void   print(ostream& os,int elem_per_line = 10,char* format = 0L) const;

};

template <class Base>
ostream& operator<<(ostream& os,const NumArray2D<Base>& array)
{
  operator<<(os,(const Array2D<Base>&) array);
  return os;
}

#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
template <class Base>
ostream_withassign& operator<<(ostream_withassign& os,const NumArray2D<Base>& array)
{
  operator<<(os,(const Array2D<Base>&) array);
  return os;
}
#endif


template <class Base> void NumArray2D<Base>::print(ostream& os,int elem_per_line,
						char* format) const
{							   
  if (dim1*dim2 == 0) 
  { 
    os << "Matrix is empty. " << endl;
    return;
  }
  
  os << "  *** MATRIX *** DIM1: " << off1 << "/" << dim1+off1-1
     << " *** DIM2: " << off2 << "/" << dim2+off2-1 << endl;
     
  int i1,i2;
  char buf[100];
  for(i1=0;i1<dim1;i1++)
  {
    for(i2=0;i2<dim2;i2++)
    {							   
      if (i2 % elem_per_line == 0) 
      { 
        int max = i2 + elem_per_line - 1; 				   
        if (max > dim2-1) max = dim2 - 1;  		   
        sprintf(buf,"\n%4i,%4i..%4i *",off1+i1,off2+i2,off2+max);
        os << buf;
      }			
      if (format)
      {      
	sprintf(buf,format,get(i1+off1,i2+off2));
	os << buf;   			   
      }
      else 
	os << get(i1 + off1,i2 + off2) << " ";
    }
  }
  os << endl;					   
}	

/*TEX
\section{Implimentation of NumArray2D}

\subsection{Matrix Arithmetic}
*/

template <class Base>
NumArray2D<Base> operator+(const NumArray2D<Base>& mat1,
			   const NumArray2D<Base>& mat2)
{
  NumArray2D<Base> result(mat1);
  return result += mat2;    
}

template <class Base>
NumArray2D<Base> operator-(const NumArray2D<Base>& mat1,const NumArray2D<Base>& mat2)
{
  NumArray2D<Base> result(mat1);
  return result -= mat2;   
}

template <class Base>
NumArray2D<Base>& NumArray2D<Base>::operator+=(const NumArray2D<Base>& mat)
{
  if(dim1 != mat.dim1 || dim2 != mat.dim2)
  	error("Operator+= -- Matrix dimensions must match.");
  int i;
  int len = dim1*dim2;
  Base* bb = mat.base;
  for(i=0; i<len;i++) 
  	base[i] += bb[i];   
  return *this;
}

template <class Base>
NumArray2D<Base>& NumArray2D<Base>::operator-=(const NumArray2D<Base>& mat)
{
  if(dim1 != mat.dim1 || dim2 != mat.dim2)
  	error("Operator-= -- Matrix dimensions must match.");
  int i;
  int len = dim1*dim2;
  Base* bb = mat.base;
  for(i=0; i<len;i++) 
  	base[i] -= bb[i];   
  return *this;
}

template <class Base>
NumArray2D<Base>& NumArray2D<Base>::operator+=(const Base value)
{
  int i;
  int len = dim1*dim2;
  for(i=0; i<len;i++) 
  	base[i] += value;   
  return *this;
}
template <class Base>
NumArray2D<Base>& NumArray2D<Base>::operator-=(const Base value)
{
  int i;
  int len = dim1*dim2;
  for(i=0; i<len;i++) 
  	base[i] -= value;   
  return *this;
}

template <class Base>
NumArray2D<Base>& NumArray2D<Base>::operator*=(const Base value)
{
  int i;
  int len = dim1*dim2;
  for(i=0; i<len;i++) 
    base[i] *= value;   
  return *this;
}

template <class Base>
NumArray2D<Base>& NumArray2D<Base>::operator/=(const Base value)
{
  int i;
  int len = dim1*dim2;
  for(i=0; i<len;i++) 
    base[i] /= value;   
  return *this;
}

template <class Base>
NumArray2D<Base> operator*(const NumArray2D<Base>& a,const NumArray2D<Base>& b)
{
  NumArray2D<Base> product(a.dim1,b.dim2);
  product.set(0);
  product.multiply(a,b);
  return product;
}

template <class Base>
void NumArray2D<Base>::multiply(const NumArray2D<Base>& b,const NumArray2D<Base>& c)
{
  if(dim1 != b.dim1)
  	error("Dimension Incompatibilty A with B in MULTIPLY",dim1,b.dim1);
    
  if(b.dim2 != c.dim1)
  	error("Dimension Incompatibilty B with C in MULTIPLY",b.dim2,b.dim1);

  if(dim2 != c.dim2)
  	error("Dimension Incompatibilty A with C in MULTIPLY",dim2,c.dim2);
    
  if (b.dim1 == 0 || b.dim2 == 0 || c.dim1 == 0 || c.dim2 == 0)
    return;
  
  if(base == b.base || base == c.base)
  	error("Overwrite in NumArray2D::multiply()");

  int i,j,k;

  for(i=0; i<dim1;i++)
  for(j=0; j<dim2;j++)
  {
    Base sum =0;

    for(k=0; k<b.dim2; k++)
      sum += b.base[k*b.dim1+i] * c.base[j*c.dim1+k];

    base[j*dim1+i] += sum;
  }
}

/*TEX
\subsection{Further Operations}
*/

template <class Base>
void NumArray2D<Base>::transpose()
{
  int i1,j1,i2,j2;

  if(dim1 != dim2) 
  	error("Transpose implemented only for quadratic matrices. ");

  Base tmp;
  for(i1=off1,i2=off2; i1< dim1+off1; i1++,i2++)
  for(j1=i1+1,j2=i2+1;j1< dim2+off1;j1++,j2++)
  {
     tmp = pbase[j2][i1];
     pbase[j2][i1] = pbase[i2][j1];
     pbase[i2][j1] = tmp;
  }
}

  /* 
     Let *this = a

     a(i,j) += sum b(k,i) c(k,j)
  */

template <class Base>
void NumArray2D<Base>::multiply_transpose1(const NumArray2D<Base>& b,
					   const NumArray2D<Base>& c)
{
  if(dim1 != b.dim2)
    error("Dimension Incompatibilty A with B in MULTIPLY-TRANSPOSE1",
	  dim1,b.dim2);
    
  if(b.dim1 != c.dim1)
    error("Dimension Incompatibilty B with C in MULTIPLY-TRANSPOSE1",
	  b.dim1,c.dim1);

  if(dim2 != c.dim2)
    error("Dimension Incompatibilty A with C in MULTIPLY-TRANSPOSE1",
	  dim2,c.dim2);

  if (b.dim1 == 0 || b.dim2 == 0 || c.dim1 == 0 || c.dim2 == 0)
    return;
  if(base == b.base || base == c.base)
    error("Overwrite in MULTIPLY-TRANSPOSE1");
    
  int i,j,k;

  for(i=0; i<dim1;i++)
  for(j=0; j<dim2;j++)
  {
    Base  sum =0;

    for(k=0; k<b.dim1;k++)
     sum += b.base[i*b.dim1+k] * c.base[j*c.dim1+k];
    base[j*dim1+i] += sum;
  }
}

template <class Base>
void NumArray2D<Base>::multiply_transpose2(const NumArray2D<Base>& b,
					   const NumArray2D<Base>& c)
{

  /* 
     Let *this = a
     a(i,j) += sum b(i,k) c(j,k)
  */

  if(dim1 != b.dim1)
    error("Dimension Incompatibilty A with B in MULTIPLY-TRANSPOSE2",
	  dim1,b.dim1);
    
  if(b.dim2 != c.dim2)
    error("Dimension Incompatibilty B with C in MULTIPLY-TRANSPOSE2",
	  b.dim2,c.dim2);

  if(dim2 != c.dim1)
    error("Dimension Incompatibilty A with C in MULTIPLY-TRANSPOSE2",
	  dim2,c.dim1);

  if (b.dim1 == 0 || b.dim2 == 0 || c.dim1 == 0 || c.dim2 == 0)
    return;
  if(base == b.base || base == c.base)
    error("Overwrite in MULTIPLY-TRANSPOSE2");
    
  int i,j,k;

  for(i=0; i<dim1;i++)
  for(j=0; j<dim2;j++)
  {
    Base  sum =0;
    
    for(k=0; k<b.dim2;k++)
     sum += b.base[k*b.dim1+i] * c.base[k*c.dim1+j];
    base[j*dim1+i] += sum;
  }
}

/*TEX
\section{NumArray2D with Access Checking}
*/
template <class Base> class D_NumArray2D : public NumArray2D<Base>
{
public:
        D_NumArray2D() {}
        D_NumArray2D(const int d1,const int d2,const int o1=0,const int o2=0) :
        NumArray2D<Base>(d1,d2,o1,o2) {}
        D_NumArray2D(const NumArray2D<Base>& mat) : NumArray2D<Base>(mat) { }
        D_NumArray2D(const D_NumArray2D<Base>& mat) : NumArray2D<Base>(mat) { }
virtual ~D_NumArray2D() {} 
//compGcc4
using NumArray2D<Base>::dim1;
using NumArray2D<Base>::dim2;
using NumArray2D<Base>::off1;
using NumArray2D<Base>::off2;
using NumArray2D<Base>::base;
using NumArray2D<Base>::pbase;
using NumArray2D<Base>::error;

D_NumArray2D<Base>& operator=(const D_NumArray2D<Base>& mat)
                   { NumArray2D<Base>::operator=(mat);
                     return *this; }

Base&   operator()(const int p1,const int p2) 
        { 
	  if(p1 < off1 || p1 >= dim1+off1 || p2 < off2 || p2 >= dim2+off2)
          error("Acces Violation",p1,p2);
	  return pbase[p2][p1];
	}
const Base& get(const int p1,const int p2) const
        { 
	  if(p1 < off1 || p1 >= dim1+off1 || p2 < off2 || p2 >= dim2+off2)
          error("Acces Violation",p1,p2);
	  return pbase[p2][p1];
	}
};

#endif
