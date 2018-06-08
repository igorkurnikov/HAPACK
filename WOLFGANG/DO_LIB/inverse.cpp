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
%  $Id: inverse.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
*/

/*TEX
\section{Implementation of Matrix Inversion}

\subsection{LU-Routines}

Now follows the implementation of LU-Routines which are taken from numerical 
recipes and modified to suit the Matrix and Vector class.
Vector class.

*/

#include <inverse.h>

#define TINY 1.0e-20;

void ludcmp(Matrix& a,int n,Array<int>& indx)
{
	int i,imax,j,k;
        double big,dum,sum,temp;
        Array<double> vv(n,1);

	for (i=1;i<=n;i++) {
	      big=0.0;
	      for (j=1;j<=n;j++)
	      	if ((temp=fabs(a(i,j))) > big) big=temp;
              if (big == 0.0) {
                cerr<<"Singular matrix in LU-Decomposition";
                exit(-1);
	      }
	      vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a(i,j);
			for (k=1;k<i;k++) sum -= a(i,k)*a(k,j);
			a(i,j)=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a(i,j);
			for (k=1;k<j;k++)
				sum -= a(i,k)*a(k,j);
			a(i,j)=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a(imax,k);
				a(imax,k)=a(j,k);
				a(j,k)=dum;
			}
	        	vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a(j,j) == 0.0) a(j,j)=TINY;
		if (j != n) {
			dum=1.0/(a(j,j));
			for (i=j+1;i<=n;i++) a(i,j) *= dum;
		}
	}
}

void ludcmp2(Matrix& a,int n,Array<int>& indx)
{
	int i,imax,j,k;
        double big,dum,sum,temp;
        Array<double> vv(n,1);

	for (i=1;i<=n;i++) {
	      big=0.0;
	      for (j=1;j<=n;j++)
	      	if ((temp=fabs(a(j,i))) > big) big=temp;
              if (big == 0.0) {
                cerr<<"Singular matrix in LU-Decomposition";
                exit(-1);
	      }
	      vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a(j,i);
			for (k=1;k<i;k++) sum -= a(k,i)*a(j,k);
			a(j,i)=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a(j,i);
			for (k=1;k<j;k++)
				sum -= a(k,i)*a(j,k);
			a(j,i)=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a(k,imax);
				a(k,imax)=a(k,j);
				a(k,j)=dum;
			}
	        	vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a(j,j) == 0.0) a(j,j)=TINY;
		if (j != n) {
			dum=1.0/(a(j,j));
			for (i=j+1;i<=n;i++) a(j,i) *= dum;
		}
	}
}

#undef TINY

void lubksb(Matrix& a,int n,Array<int>& indx,Array<double>& b)
{
	int i,ii=0,ip,j;
        double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a(i,j)*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a(i,j)*b[j];
                if ( fabs(a(i,i)) <= 1.0e-15 ) {
                   cerr<<"The determinat of the matrix is probably 0";
                   exit(-1);
                }
                else
                   b[i]=sum/a(i,i);
	}
}

void lubksb2(Matrix& a,int n,Array<int>& indx,Array<double>& b)
{
	int i,ii=0,ip,j;
        double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a(j,i)*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a(j,i)*b[j];
                if ( fabs(a(i,i)) <= 1.0e-15 ) {
                   cerr<<"The determinat of the matrix is probably 0";
                   exit(-1);
                }
                else
                   b[i]=sum/a(i,i);
	}
}

/*TEX
\subsection{Inverse}

Now follows the implementation of the routines that calculate the inverse
of a matrix and multiplication with it.
*/

#define SWAP(a,b) {float temp=(a);(a)=(b);(b)=temp;}

void inverse(Matrix& a,Array<int>& indxc,Array<int>& indxr,Array<int>& ipiv) 
{
/* This routine is a modified version of the Numerical Recipes gauss-jordan
   routine gaussj.c . */
   
  assert( a.size1() == a.size2() );  
  assert( indxc.size() == a.size1() );
  assert( indxr.size() == a.size1() );
  assert( ipiv.size() == a.size1() );
  assert( a.offset1() == 1 && a.offset2() == 1);
  assert( indxc.offset() ==1 );  
  assert( indxr.offset() ==1 );  
  assert( ipiv.offset() ==1 );  
 
  int i,icol,irow,j,k,l,ll;
  double big,dum,pivinv;

  for (j=1;j<=a.size1();j++) ipiv[j]=0;
  for (i=1;i<=a.size1();i++) {
      big=0.0;
      for (j=1;j<=a.size1();j++)
	  if (ipiv[j] != 1)
	     for (k=1;k<=a.size1();k++) {
		 if (ipiv[k] == 0) {
		    if (fabs(a(j,k)) >= big) {
	                big=fabs(a(j,k));
	                irow=j;
	                icol=k;
	            }
	         } 
                 else 
                    if (ipiv[k] > 1) { 
                       cerr<<"INVERSE FAILED: Singular Matrix-1";
	            exit(-1);
	            }
	     }
	     ++(ipiv[icol]);
	     if (irow != icol) {
		for (l=1;l<=a.size1();l++) SWAP(a(irow,l),a(icol,l))
	     }
	     indxr[i]=irow;
	     indxc[i]=icol;
	     if (a(icol,icol) == 0.0) { 
                cerr<<"INVERSE FAILED: Singular Matrix-2";
                exit(-1);
	     }	
	     pivinv=1.0/a(icol,icol);
	     a(icol,icol)=1.0;
	     for (l=1;l<=a.size1();l++) a(icol,l) *= pivinv;
	     for (ll=1;ll<=a.size1();ll++)
		 if (ll != icol) {
		    dum=a(ll,icol);
		    a(ll,icol)=0.0;
		    for (l=1;l<=a.size1();l++) a(ll,l) -= a(icol,l)*dum;
		 }
      }
      for (l=a.size1();l>=1;l--) {
	  if (indxr[l] != indxc[l])
	      for (k=1;k<=a.size1();k++)
		  SWAP(a(k,indxr[l]),a(k,indxc[l]));
	}
}

#undef SWAP

void inverse(Matrix& a)
{  
  Array<int>  indxc(a.size1(),1);
  Array<int>  indxr(a.size1(),1);
  Array<int>  ipiv(a.size1(),1);		 

  int o1 =a.offset1();
  int o2 =a.offset2();
  a.rebase(1,1);

  inverse(a,indxc,indxr,ipiv);

  a.rebase(o1,o2);
}	
	 
void inverse(Matrix& a,Matrix& y,Array<int>& indx,Array<double>& col)
{
/* This routine uses the Numerical Recipes LU-decomposition routines 
   ludcmp.c and lubksb.c in a modified way. */

   assert( &a != &y );
   assert( a.size1() == a.size2() );
   assert( y.size1() == y.size2() && y.size1() == a.size1() );
   assert( indx.size() == a.size1() );
   assert( col.size()  == a.size1() );
   assert( a.offset1() == 1 && a.offset2() == 1 );
   assert( y.offset1() == 1 && y.offset2() == 1 );
   assert( indx.offset() ==1);
   assert( col.offset() == 1);

   int i,j;
 
   ludcmp(a,a.size1(),indx);
   for(j=1;j<=a.size1();j++) {  
      for(i=1;i<=a.size1();i++) col[i]=0.0;
      col[j]=1.0;
      lubksb(a,a.size1(),indx,col);
      for(i=1;i<=a.size1();i++) y(i,j)=col[i];
   }
}

void inverse(Matrix& a,Matrix& y)
{
  
  Array<int>    indx(a.size1(),1);
  Array<double> col(a.size1(),1);
  
  int ao1 = a.offset1();
  int ao2 = a.offset2();
  int yo1 = y.offset1();
  int yo2 = y.offset2();
  
  a.rebase(1,1);        
  y.rebase(1,1);
  
  inverse(a,y,indx,col);
  
  a.rebase(ao1,ao2);
  y.rebase(yo1,yo2);
}  

void multiply_inverse1(Matrix& a,Matrix& b,const Matrix& c,
                                         Array<int>& indx,Array<double> col)
{
/* This routine uses the Numerical Recipes LU-decomposition routines 
   ludcmp.c and lubksb.c in a modified way. */

   assert( &a != &b && &a != &c && &b != &c );
   assert( b.size1() == b.size2() );
   assert( a.size1() == b.size1() );
   assert( b.size2() == c.size1() );
   assert( a.size2() == c.size2() );
   assert( indx.size() == b.size1() );
   assert( col.size()  == b.size1() );
   assert( a.offset1() == 1 && a.offset2() == 1 );
   assert( b.offset1() == 1 && b.offset2() == 1 );
   assert( c.offset1() == 1 && c.offset2() == 1 );
   assert( indx.offset() ==1);
   assert( col.offset() == 1);

   int i,j;

   ludcmp(b,b.size1(),indx);
   for(j=1;j<=c.size2();j++) {
      for(i=1;i<=b.size1();i++) col[i]=c.get(i,j);  
      lubksb(b,b.size1(),indx,col);
      for(i=1;i<=b.size1();i++) a(i,j)=col[i];
   }
}

void multiply_inverse1(Matrix& a,Matrix& b,Matrix& c)
{
   int ao1   = a.offset1();
   int ao2   = a.offset2();
   int bo1   = a.offset1();
   int bo2   = a.offset2();
   int co1   = a.offset1();
   int co2   = a.offset2();
   
   a.rebase(1,1);        
   b.rebase(1,1);
   c.rebase(1,1);
   
   Array<int>    indx(b.size1(),1);
   Array<double> col(b.size1(),1);

   multiply_inverse1(a,b,c,indx,col);
   
   a.rebase(ao1,ao2);
   b.rebase(bo1,bo2);
   c.rebase(co1,co2);
}

Matrix multiply_inverse1(Matrix& b,Matrix& c)
{
   int bo1   = b.offset1();
   int bo2   = b.offset2();
   int co1   = c.offset1();
   int co2   = c.offset2();

   Matrix a(b.size1(),c.size2(),1,1);   
   b.rebase(1,1);
   c.rebase(1,1);
   
   Array<int>    indx(b.size1(),1);
   Array<double> col(b.size1(),1);

   multiply_inverse1(a,b,c,indx,col);
   
   b.rebase(bo1,bo2);
   c.rebase(co1,co2);
   
   return a;
}


void multiply_inverse2(Matrix& a,const Matrix& b,Matrix& c,
                                         Array<int>& indx,Array<double> col)
{
/* This routine uses the Numerical Recipes LU-decomposition routines 
   ludcmp.c and lubksb.c in a modified way.
   It is based on the identity B*C=(C^T*B^T)^T  (T=transpose) */

   assert( &a != &b && &a != &c && &b != &c );
   assert( c.size1() == c.size2() );
   assert( a.size1() == b.size1() );
   assert( b.size2() == c.size1() );
   assert( a.size2() == c.size2() );
   assert( indx.size() == c.size1() );
   assert( col.size()  == c.size1() );
   assert( a.offset1() == 1 && a.offset2() == 1 );
   assert( b.offset1() == 1 && b.offset2() == 1 );
   assert( c.offset1() == 1 && c.offset2() == 1 );
   assert( indx.offset() ==1);
   assert( col.offset() == 1);

   int i,j;

   ludcmp2(c,c.size1(),indx);
   for(j=1;j<=b.size1();j++) {
      for(i=1;i<=c.size1();i++) col[i]=b.get(j,i);  
      lubksb2(c,c.size1(),indx,col);
      for(i=1;i<=c.size1();i++) a(j,i)=col[i];
   }
}

void multiply_inverse2(Matrix& a,Matrix& b,Matrix& c)
{
   int ao1   = a.offset1();
   int ao2   = a.offset2();
   int bo1   = a.offset1();
   int bo2   = a.offset2();
   int co1   = a.offset1();
   int co2   = a.offset2();
   
   a.rebase(1,1);        
   b.rebase(1,1);
   c.rebase(1,1);
   
   Array<int>    indx(c.size1(),1);
   Array<double> col(c.size1(),1);

   multiply_inverse2(a,b,c,indx,col);
   
   a.rebase(ao1,ao2);
   b.rebase(bo1,bo2);
   c.rebase(co1,co2);
}

Matrix multiply_inverse2(Matrix& b,Matrix& c)
{
   int bo1   = b.offset1();
   int bo2   = b.offset2();
   int co1   = c.offset1();
   int co2   = c.offset2();

   Matrix a(b.size1(),c.size2(),1,1);   
   b.rebase(1,1);
   c.rebase(1,1);
   
   Array<int>    indx(c.size1(),1);
   Array<double> col(c.size1(),1);

   multiply_inverse2(a,b,c,indx,col);
   
   b.rebase(bo1,bo2);
   c.rebase(co1,co2);
   
   return a;
}

   
void multiply_inverse(Vector& a,Matrix& B,const Vector& c,Array<int> indx)
{
   assert( &a != &c );
   assert( B.size1() == B.size2() );
   assert( a.size() == B.size1() );
   assert( B.size2() == c.size() );
   assert( indx.size() == B.size1() );
   assert( a.offset() == 1 );
   assert( B.offset1() == 1 && B.offset2() == 1 );
   assert( c.offset() == 1 );
   assert( indx.offset() ==1);
   
   for (int i=1;i<=a.size();i++) a[i]=c(i);
   
   ludcmp(B,B.size1(),indx);
   lubksb(B,B.size1(),indx,a);
}


void multiply_inverse(Vector& a,Matrix& B,Vector& c)
{   
   int  ao = a.rebase(1);
   int bo1 = B.offset1();
   int bo2 = B.offset2();
   int  co = c.rebase(1);
   B.rebase(1,1);   

   Array<int> indx(B.size1(),1);

   multiply_inverse(a,B,c,indx);
    
   a.rebase(ao);
   B.rebase(bo1,bo2);
   c.rebase(co);
}

Vector multiply_inverse(Matrix& B,Vector& c)
{   
   int bo1 = B.offset1();
   int bo2 = B.offset2();
   int  co = c.rebase(1);
   B.rebase(1,1);   

   Vector a(B.size1(),1);
   
   Array<int> indx(B.size1(),1);

   multiply_inverse(a,B,c,indx);
    
   B.rebase(bo1,bo2);
   c.rebase(co);
   return a;
}   














