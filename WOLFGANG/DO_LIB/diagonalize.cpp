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
%  $Id: diagonalize.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\subsection{Implementation}  
*/

#include <diagonalize.h>
#include <errorip.h>


//  Machine dependent diagonalizers.

extern "C" {
#ifdef SUN
  void devcsf_(int*,double*,int*,double*,double*,int*);
  void dgvcsp_(int*,double*,int*,double*,int*,double*,double*,int*);
#define Mdsvdc dsvdc_
#endif
#ifdef XLC
#define Mdsvdc dsvdc
#endif
  void Mdsvdc(double*,int*,int*,int*,double*,double*,double*,int*,
              double*,int*,double*,int*,int*);
}

#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))

int diagonalize( Matrix& mat, Matrix& eigen, Vector& eval) 
{  
  if (mat.size1() == 0 || mat.size2() == 0) return DIAG_EMPTY;
  
  int size = mat.size1();
  if ( size != mat.size2() ) {
    cerr << " diagonalize() : Matrix mat must be quadratic.\n";
    abort();
  }
  if ( eigen.size1() != size ) {
    cerr << " diagonalize() : Matries mat and eigen do not match.\n";
    abort();
  }
  if ( eigen.size2() != size ) {
    cerr << " diagonalize() : Matries mat and eigen do not match.\n";
    abort();
  }
  if ( eval.size() != size ) {
    cerr << " diagonalize() : Matries mat and eval do not match.\n";
    abort();
  }
  
  int i, j, k;
  
#ifdef SUN_IMSL

  // extract pointers to data fields:

  double* matptr = &mat(0,0);
  double* eigptr = &eigen(0,0);
  double*  evptr = &eval[0];

  devcsf_( &size, matptr, &size, evptr, eigptr, &size);

  //  Normalize the eigenvectors

  double norm;
  for( i=0; i<size; i++) {
    norm = 0;
    for( j=0; j<size; j++ )
      norm += eigen(i,j)*eigen(i,j);
    norm = 1 / sqrt(norm);
    for( j=0; j<size; j++)
      eigen(i,j) *= norm;
  }

  return;

#else

  eigen = mat;
 
  Vector tmp( size );  
  Matrix orig( size, size );
  orig = mat;  
  tred2old( eigen, eval, tmp );
  tqli(eval, tmp, eigen );

  straightsort(eval,eigen );

  //  Check the eigensolutions
  
  double s;
  for(i=0; i<size;i++) 
    for(j=0; j<size;j++) 
    {
      s = 0.0;
      for(k=0; k<size;k++)
	s += orig(j,k) * eigen(k,i);
      if (fabs(s - eval[i]*eigen(j,i)) > 1E-7)
	return DIAG_DEV;
    }
  return DIAG_OK;
#endif
}

//################################################################

void orthonormalize( Matrix& S, Matrix& u ) {

  if (S.size1() == 0 || S.size2() == 0) return;

  int i, j, k;
  
  int n = S.size1();
  Vector eval( n );
  Matrix evec( n, n );
  u.reset( n, n );
  u.set(0);
  
  // diagonalize S, and then calculate 1/sqrt(S).

  diagonalize( S, evec, eval);

  for(i=0;i<n;i++)
    eval[i] = 1.0/sqrt(eval[i]);

  double sum;
  for(i=0;i<n;i++)
  for(j=0;j<n;j++) {
    sum = 0;
    for(k=0;k<n;k++)
      sum += evec(j,k)*eval[k]*evec(i,k);
    u(i,j) = sum;
  }

  
}

//################################################################

void geneva( NumArray2D<double>& ham, NumArray2D<double>& S,
	     NumArray2D<double>& evec, NumArray<double>& eval ) {

  //  Check the dimensions.

  int size = ham.size1();
  if ( size != ham.size2() ) {
    cerr << " geneva() : Matrix ham must be quadratic.\n";
    abort();
  }
  if ( S.size2() != size ) {
    cerr << " geneva() : Matrices ham and S do not match.\n";
    abort();
  }
  if ( evec.size2() != size ) {
    cerr << " geneva() : Matrices ham and evec do not match.\n";
    abort();
  }
  if ( eval.size() != size ) {
    cerr << " geneva() : Matrices ham and eval do not match.\n";
    abort();
  }

#ifdef SUN_IMSL
   dgvcsp_( &size, &ham(0,0), &size, &S(0,0), &size,
           &eval[0], &evec(0,0), &size);
#else
  // do with orthonoramlize etc
//  error( " geneva() : not implemented yet" );
#endif

}

//################################################################
/*
void LR_svd( const Matrix& mat, Matrix& U, Matrix& V, Vector& eval ) {

  const NumArray2D<double>& mat1 = mat;  
  NumArray2D<double>& U1 = U;
  NumArray2D<double>& V1 = V;
  NumArray<double>& eval1 = eval;

  LR_svd( mat1, U1, V1, eval1 );

}

//################################################################

void LR_svd( const NumArray2D<double>& mat, NumArray2D<double>& U,
	     NumArray2D<double>& V, NumArray<double>& eval ) {

  //  Check the dimensions.

  int size = mat.size1();
  if ( size != mat.size2() ) {
    cerr << " LR_svd() : Matrix mat must be quadratic.\n";
    abort();
  }
  if ( U.size2() != size ) {
    cerr << " LR_svd() : Matrices mat and U do not match.\n";
    abort();
  }
  if ( V.size2() != size ) {
    cerr << " LR_svd() : Matrices mat and V do not match.\n";
    abort();
  }
  if ( eval.size() != size ) {
    cerr << " LR_svd() : Matrices mat and eval do not match.\n";
    abort();
  }

  //  Create and dimension the workspaces

  Array2D<double> mat_copy = mat;
  Array<double> e( size ), work( size );
  int job = 11;
  int info;
  Mdsvdc( &mat_copy(0,0), &size, &size, &size, &eval[0], &e[0], &U(0,0), &size,
               &V(0,0), &size, &work[0], &job, &info );
  if ( info ) {
    cerr << " LR_svd() : info = " << info <<"\n";
    abort();
  }

}
*/
//################################################################

void straightsort( Vector& eval, Matrix& eigen ) {

  int i;
  int n = eval.size();  

  for(i=0;i<n-1;i++) {
    double v  = eval[i];
    int j;
    int index = i;
    for(j=i+1;j<n;j++)
      if (eval[j] < v) {
	index = j;
	v     = eval[index];
      }
    if (index != i) { // swap index and i 
      double t;
      t           = eval[i];
      eval[i]     = eval[index];
      eval[index] = t;
      for(j=0;j<n;j++) {
        t = eigen(j,i);
	eigen(j,i) = eigen(j,index);
	eigen(j,index) = t;
      }
    }
  }
}

//################################################################

int tqli( Vector& d, Vector& e, Matrix& z ) 
{
  int n = z.size1();
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;
  
  for (i=1;i<n;i++) 
    e[i-1]=e[i];
  e[n-1]=0.0;
  for (l=0;l<n;l++) 
  {
    iter=0;
    do 
    {
      for (m=l;m<n-1;m++) 
      {
	dd=fabs(d[m])+fabs(d[m+1]);
	if (fabs(e[m]) + dd == dd) break;
      }
      
      if (m != l) 
      {
	if (iter++ == 100) 
	  return DIAG_TQLI_ITER;
	
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=sqrt((g*g)+1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) 
	{
	  f=s*e[i];
	  b=c*e[i];
	  if (fabs(f) >= fabs(g)) 
	  {
	    c=g/f;
	    r=sqrt((c*c)+1.0);
	    e[i+1]=f*r;
	    c *= (s=1.0/r);
	  } 
	  else 
	  {
	    s=f/g;
	    r=sqrt((s*s)+1.0);
	    e[i+1]=g*r;
	    s *= (c=1.0/r);
	  }
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  p=s*r;
	  d[i+1]=g+p;
	  g=c*r-b;
	  /* Next loop can be omitted if eigenvectors not wanted */
	  for (k=0;k<n;k++) 
	  {
	    f=z(k,i+1);
	    z(k,i+1)=s*z(k,i)+c*f;
	    z(k,i)=c*z(k,i)-s*f;
	  }
	}
	d[l]=d[l]-p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m != l);
  }
}



void tred2old( Matrix& a, Vector& d, Vector& e ) 
{
  int n = a.size1();
  int l,k,j,i;
  double scale,hh,h,g,f;

  for (i=n;i>=2;i--) 
  {
    l=i-1;
    h=scale=0.0;
    if (l > 1) 
    {
      for (k=1;k<=l;k++)
	scale += fabs(a(i-1,k-1));
      if (scale == 0.0)
	e[i-1]=a(i-1,l-1);
      else {
	for (k=1;k<=l;k++) 
	{
	  a(i-1,k-1) /= scale;
	  h += a(i-1,k-1)*a(i-1,k-1);
	}
	f=a(i-1,l-1);
	g = f>0 ? -sqrt(h) : sqrt(h);
	e[i-1]=scale*g;
	h -= f*g;
	a(i-1,l-1)=f-g;
	f=0.0;
	for (j=1;j<=l;j++) 
	{
	  /* Next statement can be omitted if eigenvectors not wanted */
	  a(j-1,i-1)=a(i-1,j-1)/h;
	  g=0.0;
	  for (k=1;k<=j;k++)
	    g += a(j-1,k-1)*a(i-1,k-1);
	  for (k=j+1;k<=l;k++)
	    g += a(k-1,j-1)*a(i-1,k-1);
	  e[j-1]=g/h;
	  f += e[j-1]*a(i-1,j-1);
	}
	hh=f/(h+h);
	for (j=1;j<=l;j++) {
	  f=a(i-1,j-1);
	  e[j-1]=g=e[j-1]-hh*f;
	  for (k=1;k<=j;k++)
	    a(j-1,k-1) -= (f*e[k-1]+g*a(i-1,k-1));
	}
      }
    } else
      e[i-1]=a(i-1,l-1);
    d[i-1]=h;
  }
  /* Next statement can be omitted if eigenvectors not wanted */
  d[1-1]=0.0;
  e[1-1]=0.0;
  /* Contents of this loop can be omitted if eigenvectors not
     wanted except for statement d[i]=a[i][i]; */
  
  for (i=1;i<=n;i++) 
  {
    l=i-1;
    if (d[i-1]) 
    {
      for (j=1;j<=l;j++) 
      {
	g=0.0;
	for (k=1;k<=l;k++)
	  g += a(i-1,k-1)*a(k-1,j-1);
	for (k=1;k<=l;k++)
	  a(k-1,j-1) -= g*a(k-1,i-1);
      }
    }
    d[i-1]=a(i-1,i-1);
    
    a(i-1,i-1)=1.0;
    for (j=1;j<=l;j++) 
    { 
      a(j-1,i-1)= 0; 
      a(i-1,j-1)=0.0; 
    }
  }    
}

