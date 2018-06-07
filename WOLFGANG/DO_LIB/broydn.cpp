#include <math.h>
#include <errorip.h>
#include <numer.h>
#define MAXITS 200
#define EPS 1.0e-7
#define TOLF 1.0e-4
#define TOLX EPS
#define STPMX 100.0
#define TOLMIN 1.0e-6

int nn;
float *fvec;
void (*nrfuncv)();

class Broydn
{
  int n;
  
  Vector c;
  Vector d;
  Vector fvcold;
  Vector g;
  Vector p;  
  Matrix qt;
  Matrix r;
  Matrix q;
  Vector s;
  Vector w;
  Vector t;
  Vector xold;
  Vector fvec;

double SQR (const double a)                  { return a*a; }
double FMAX(const double a,const double b)   { return (a > b) ? a : b; }
double SIGN(const double a,const double b)   { return (b >= 0) ? a : -a; }

void   rotate(Matrix& r,Matrix& qt,int i,double a,double b);
public:

  Broydn(Vector& x,IVector& check);
  
};


Broydn::Broydn(Vector& x,IVector& check) 
  : c(1,x.size()),d(1,x.size()),fvcold(1,x.size()),g(1,x.size()),
    p(1,x.size()),qt(x.size(),x.size(),1,1),q(x.size(),x.size(),1,1),
    r(x.size(),x.size(),1,1),
    s(1,x.size()),t(1,x.size()),w(1,x.size()),xold(1,x.size()),fvec(1,x.size())
{
  n = x.size();
  
  //double fmin();
  //  void fdjac(),lnsrch(),qrdcmp(),qrupdt(),rsolv();
  
  int i,its,j,k,restrt,sing,skip;
  double den,f,fold,stpmax,sum,temp,test,*c,*d,*fvcold;
  
  nn=n;
  // nrfuncv=vecfunc;
  f=fmin(x);
  
  /* find maximum of fmin */
  test=0.0;
  for (i=1;i<=n;i++)
    if (fabs(fvec[i]) > test)
      test=fabs(fvec[i]);
  
  if (test < 0.01*TOLF) 
  {
    // *check=0;
    return; 
  }
  
  for (sum=0.0,i=1;i<=n;i++) 
    sum += SQR(x[i]);
  stpmax=STPMX*FMAX(sqrt(sum),(double) n);
  restrt=1;
  for (its=1;its<=MAXITS;its++) 
  {
    if (restrt) 
    {
      fdjac(n,x,fvec,r,vecfunc);
      qrdcmp(r,n,c,d,&sing);
      if (sing) 
	error("singular Jacobian in broydn");
      qt.set(0);
      for (i=1;i<=n;i++) 
	qt(i,i) = 1.0;

      for (k=1;k<n;k++)
      {
	if (c[k]) 
	{
	  for (j=1;j<=n;j++) 
	  {
	    sum=0.0;
	    for (i=k;i<=n;i++)
	      sum += r(i,k)*qt(i,j);
	    sum /= c[k];
	    for (i=k;i<=n;i++)
	      qt(i,j) -= sum*r(i,k);
	  }
	}
      }
      for (i=1;i<=n;i++) 
      {
	r(i,i)=d[i];
	for (j=1;j<i;j++) 
	  r(i,j) =0.0;
      }
    } 
    else 
    {
      for (i=1;i<=n;i++) 
	s[i]=x[i]-xold[i];
      for (i=1;i<=n;i++) 
      {
	for (sum=0.0,j=i;j<=n;j++) 
	  sum += r(i,j)*s[j];
	t[i]=sum;
      }
      skip=1;
      for (i=1;i<=n;i++) 
      {
	for (sum=0.0,j=1;j<=n;j++) 
	  sum += qt(j,i)*t[j];
	w[i]=fvec[i]-fvcold[i]-sum;
	if (fabs(w[i]) >= EPS*(fabs(fvec[i])+fabs(fvcold[i]))) 
	  skip=0;
	else 
	  w[i]=0.0;
      }
      if (!skip) 
      {
	for (i=1;i<=n;i++) 
	{
	  for (sum=0.0,j=1;j<=n;j++) 
	    sum += qt(i,j)*w[j];
	  t[i]=sum;
	}
	for (den=0.0,i=1;i<=n;i++) 
	  den += SQR(s[i]);
	for (i=1;i<=n;i++) 
	  s[i] /= den;
	qrupdt(r,qt,n,t,s);
	for (i=1;i<=n;i++) 
	{
	  if (r[i][i] == 0.0) 
	    error("r singular in broydn");
	  d[i]=r[i][i];
	}
      }
    }
    for (i=1;i<=n;i++) 
    {
      for (sum=0.0,j=1;j<=n;j++) 
	sum += qt(i,j)*fvec[j];
      g[i]=sum;
    }
    for (i=n;i>=1;i--) 
    {
      for (sum=0.0,j=1;j<=i;j++) 
	sum += r[j][i]*g[j];
      g[i]=sum;
    }
    for (i=1;i<=n;i++) 
    {
      xold[i]=x[i];
      fvcold[i]=fvec[i];
    }
    fold=f;
    for (i=1;i<=n;i++) 
    {
      for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
      p[i] = -sum;
    }
    rsolv(r,n,d,p);
    lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin);
    test=0.0;
    for (i=1;i<=n;i++)
      if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
    if (test < TOLF) 
    {
      // *check=0;
      return;
    }
    if (*check) 
    {
      if (restrt) 
	return; 
      else 
      {
	test=0.0;
	den=FMAX(f,0.5*n);
	for (i=1;i<=n;i++) 
	{
	  temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
	  if (temp > test) test=temp;
	}
	if (test < TOLMIN) 
	  return; 
	else restrt=1;
      }
    } 
    else 
    {
      restrt=0;
      test=0.0;
      for (i=1;i<=n;i++) 
      {
	temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
	if (temp > test) test=temp;
      }
      if (test < TOLX) 
	return; 
    }
  }
  error("MAXITS exceeded in broydn");
}


void Broydn::rotate(Matrix& r,Matrix& qt,int i,double a,double b)
{
  int j;
  float c,fact,s,w,y;

  if (a == 0.0) 
  {
    c=0.0;
    s=(b >= 0.0 ? 1.0 : -1.0);
  } else if (fabs(a) > fabs(b)) 
  {
    fact=b/a;
    c=SIGN(1.0/sqrt(1.0+(fact*fact)),a);
    s=fact*c;
  } else 
  {
    fact=a/b;
    s=SIGN(1.0/sqrt(1.0+(fact*fact)),b);
    c=fact*s;
  }
  for (j=i;j<=n;j++) 
  {
    y=r[i][j];
    w=r[i+1][j];
    r[i][j]=c*y-s*w;
    r[i+1][j]=s*y+c*w;
  }
  for (j=1;j<=n;j++) 
  {
    y=qt[i][j];
    w=qt[i+1][j];
    qt[i][j]=c*y-s*w;
    qt[i+1][j]=s*y+c*w;
  }
}


void Broydn::fdjac(Vector& x,Vector& fvec,Matrix& df)
{
  const double eps = 1E-4;
  int i,j;
  double h,temp;
  
  Vector f(n,1);
  for (j=1;j<=n;j++) 
  {
    temp=x[j];
    h=EPS*fabs(temp);
    if (h == 0.0) h=EPS;
    x[j]=temp+h;
    h=x[j]-temp;
    (*vecfunc)(n,x,f);
    x[j]=temp;
    for (i=1;i<=n;i++) 
      df[i][j]=(f[i]-fvec[i])/h;
  }
}


#undef MAXITS
#undef EPS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
