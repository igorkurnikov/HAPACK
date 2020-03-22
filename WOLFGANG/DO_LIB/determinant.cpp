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
% $Id: determinant.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%

\section{Implementation of the Determinant of a Matrix}

\subsection{LU-Routines}

Now follows the implementation of the LU-Routines which are taken from 
numerical recipes and modified to suit the Matrix and Vector class.
Vector class.

*/

#include "io_incl.h"
#include "determinant.h"

#define TINY 1.0e-20;

void ludcmp(Matrix& a,int n,double& d)
{
/* This routine uses the Numerical Recipes LU-decomposition routine 
   ludcmp.c in a modified way. */

	int i,imax,j,k;
        double big,dum,sum,temp;
        Array<double> vv(n,1);

	d=1.0;
	for (i=1;i<=n;i++) {
	      big=0.0;
	      for (j=1;j<=n;j++)
	      	if ((temp=fabs(a(i,j))) > big) big=temp;
              if (big == 0) {
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
			d = -d;
			vv[imax]=vv[j];
		}
		if (a(j,j) == 0.0) a(j,j)=TINY;
		if (j != n) {
			dum=1.0/(a(j,j));
			for (i=j+1;i<=n;i++) a(i,j) *= dum;
		}
	}
}

#undef TINY


/*TEX
\subsection{Determinant}

This is the way how we calculate the determinant.
*/

double determinant(Matrix& a)  
{
  assert( a.size1() == a.size2() );
  
  int i;
  a.rebase(1,1);
  double d;
  
  ludcmp(a,a.size1(),d);

  for(i=1;i<=a.size1();i++) d *= a(i,i);
  return d;
}  











