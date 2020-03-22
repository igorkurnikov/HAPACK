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
%*/

#define BLAS_IMPLEMENT_CPP 

#include "numer.h"

inline void Matrix::multiply( const Matrix& b, const Matrix& c ) {

  if ( size1() == 0 || size2() == 0 ) return;

  if ( base == b.base || base == c.base )
     error("Overwrite in Matrix::multiply");

  int b_dim1 = b.dim1;
  int b_dim2 = b.dim2;
  int c_dim1 = c.dim1;
  int c_dim2 = c.dim2;

  if ( b_dim2 != c_dim1 )
  	error("Dimension Incompatibilty B with C in multiply", b_dim2, c_dim1 );
  if ( dim2 != c_dim2 )
  	error("Dimension Incompatibilty A with C in multiply", dim2, c_dim2 );
  if ( dim1 != b_dim1 )
  	error("Dimension Incompatibilty A with B in multiply", dim1, b_dim1 );

#ifndef SUN

  double d_one = 1.0;
  char N = 'N';
  Mdgemm( &N, &N, &b_dim1, &c_dim2, &b_dim2, 
          &d_one, b.base, &b_dim1, c.base, &c_dim1, &d_one, base, &dim1 );

#else

  Matrix temp = *this;
  Mdmrrrr( &b_dim1, &b_dim2, b.base, &b_dim1, 
           &c_dim1, &c_dim2, c.base, &c_dim1, &dim1, &dim2, base, &dim1 );
  *this += temp;

#endif
    
}

inline void Matrix::multiply_transpose1( const Matrix& b, const Matrix& c ) {

  if ( size1() == 0 || size2() == 0 ) return;

  if( base == b.base || base == c.base )
  	error("Overwrite in Matrix::multiply_transpose1");

  int b_dim1 = b.dim1;
  int b_dim2 = b.dim2;
  int c_dim1 = c.dim1;
  int c_dim2 = c.dim2;

  if ( dim1 != b_dim2 )
  	error("Dimension Incompatibilty A with B in multiply_transpose1", dim1, b_dim2 );
  if ( b_dim1 != c_dim1 )
  	error("Dimension Incompatibilty B with C in multiply_transpose1", b_dim1, c_dim1 );
  if ( dim2 != c_dim2 )
  	error("Dimension Incompatibilty A with C in multiply_transpose1", dim2, c_dim2 );

#ifndef SUN

  double d_one = 1.0;
  char N = 'N';
  char T = 'T';
  Mdgemm( &T, &N, &b_dim2, &c_dim2, &b_dim1, 
          &d_one, b.base, &b_dim1, c.base, &c_dim1, &d_one, base, &dim1 );

#else

  Matrix temp = *this;
  Mdmxtyf( &b_dim1, &b_dim2, b.base, &b_dim1, 
           &c_dim1, &c_dim2, c.base, &c_dim1, &dim1, &dim2, base, &dim1 );
  *this += temp;

#endif

}


inline void Matrix::multiply_transpose2( const Matrix& b, const Matrix& c ) {

  if ( size1() == 0 || size2() == 0 ) return;

  if( base == b.base || base == c.base )
  	error("Overwrite in Matrix::multiply_transpose2");

  int b_dim1 = b.dim1;
  int b_dim2 = b.dim2;
  int c_dim1 = c.dim1;
  int c_dim2 = c.dim2;

  if ( dim1 != b_dim1 )
  	error("Dimension Incompatibilty A with B in multiply_transpose2", dim1, b_dim1 );
  if ( b_dim2 != c_dim2 )
  	error("Dimension Incompatibilty B with C in multiply_transpose2", b_dim2, c_dim2 );
  if ( dim1 != c_dim2 )
  	error("Dimension Incompatibilty A with C in multiply_transpose2", dim1, c_dim2 );

#ifndef SUN

  double d_one = 1.0;
  char N = 'N';
  char T = 'T';
  Mdgemm( &N, &T, &b_dim1, &c_dim1, &b_dim2, 
          &d_one, b.base, &b_dim1, c.base, &c_dim1, &d_one, base, &dim1 );

#else

  Matrix temp = *this;
  Mdmxytf( &b_dim1, &b_dim2, b.base, &b_dim1, 
           &c_dim1, &c_dim2, c.base, &c_dim1, &dim1, &dim2, base, &dim1 );
  *this += temp;

#endif

}


/*TEX 
\subsection{Vector member functions}
*/

Vector& Vector::operator=(const Vector& v) {
  if (v.len != len || v.off != off) reset(v.len,v.off);
  int one = 1;
  Mdcopy( &len, v.base, &one, base, &one);
  return *this;
}

inline Vector& Vector::operator+=(const Vector& v) {
  if (v.len != len)
    error( "Vector::operator+=() Length does not match." );
  int one = 1;
  double d_one = 1.0;
  Mdaxpy( &len, &d_one, v.get_base(), &one, base, &one );
  return *this;
}

inline double Vector::operator*(const Vector& v) {
  if ( v.len != len )
    error( "Vector::Opertor*() -- Length does not match. " );
  int one = 1;
  return Mddot( &len, base, &one, v.base, &one );
}

void Vector::multiply(const Matrix& mat, const Vector& vec) {
  if ( len != mat.size1() )
    error("MULTIPLY: Output Vector Length not compatible with Matrix.");
  if ( vec.size() != mat.size2() )
    error("MULTIPLY: Input Vector Length not compatible with Matrix.");
  char TRANS = 'N';
  int vec_len = vec.size();
  int one = 1;
  double d_one = 1.0;
  Mdgemv( &TRANS, &len, &vec_len, &d_one, mat.get_base(), &len,
          vec.get_base(), &one, &d_one, base, &one );
}

//  The following do not exist in BLAS so we unroll the loop ourselves.

void Vector::set(const double value) {
  int floor = len/4;
  int remainder = len%4;
  int i;
  if ( remainder ) {
    for ( i=0; i<remainder; i++ )
      base[i] = value;
  }
  if ( floor ) {
    for ( i=remainder; i<len; i+=4 ) {
      base[i] = value;
      base[i+1] = value;
      base[i+2] = value;
      base[i+3] = value;
    }
  }
}						

inline Vector& Vector::operator+=(const double value) {
  int floor = len/4;
  int remainder = len%4;
  int i;
  if ( remainder ) {
    for ( i=0; i<remainder; i++ )
      base[i] += value;
  }
  if ( floor ) {
    for ( i=remainder; i<len; i+=4 ) {
      base[i] += value;
      base[i+1] += value;
      base[i+2] += value;
      base[i+3] += value;
    }
  }
  return *this;
}


/*TEX 
\subsection{Vector user functions}
*/

inline Vector operator+(const Vector& vec1,const Vector& vec2) {
  int vec1_len = vec1.size();
  if ( vec1_len != vec2.size() || vec1.offset() != vec2.offset() )
    vec1.error( "Vector::operator+() : lengths or offsets not equal" );
  int one = 1;
  double d_one = 1.0;
  Vector temp( vec1.size(), vec1.offset() );
  Mdcopy( &vec1_len, vec1.get_base(), &one, temp.get_base(), &one );
  Mdaxpy( &vec1_len, &d_one, vec2.get_base(), &one, temp.get_base(), &one );
  return temp;
}

inline Vector operator-(const Vector& vec1,const Vector& vec2) {
  int vec1_len = vec1.size();
  if ( vec1_len != vec2.size() || vec1.offset() != vec2.offset() )
    vec1.error( "Vector::operator-() : lengths or offsets not equal" );
  int one = 1;
  double d_one = 1.0;
  double d_minus_one = -1.0;
  Vector temp( vec1.size(), vec1.offset() );
  Mdcopy( &vec1_len, vec1.get_base(), &one, temp.get_base(), &one );
  Mdaxpy( &vec1_len, &d_minus_one, vec2.get_base(), &one, temp.get_base(), &one );
  return temp;
}

