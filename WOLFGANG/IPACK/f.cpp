/*TEX
%
% IPACK - QCHEM
% Quantum Chemistry Project: Integral Package
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is proprietary software. Unlicensed use and distribution
% of this program or parts thereof are illegal. For details on the 
% license, see the LICENSE section in the file "ipack.c" distributed 
% with this package or write to: 
%      Wolfgang Wenzel, Theoretical Physics I, 
%      University of Dortmund,
%      D-44221 Dortmund, Germany 
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the LICENSE
% for more details.
%
% $Id: f.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Definitions} 
  
The main function of this chapter is the function \name{F}, which evaluates 
the auxiliary error integral using a sixth-order taylor expansion in
the same spirit as the function \name{FTaylor} below.  
 
Additionally the following methods of comutation have been implemented:
\begin{itemize} 
\item
   \name{numint} -- uses a numerical integration routine to compute 
                F(T,m) based on a trapezoidal rule with convergence
                checks.

   Let (t,t+dt) be an interval, the integral is approximated
   as: 
   \be
     I_{0} = \frac{dt}{2} (f(t) + f(t+dt)).
   \ee 
   Suppose we subdivide by two:
   \be
      \frac{dt}{4} (f(t) + f(t+dt/2)) + \frac{dt}{4} (f(t+dt/2) + f(t+dt))
      = \frac{1}{2} I_{0} + \frac{dt}{2} f(t+dt/2) 
   \ee
 
   The above estimate is recursively evaluated until convergence 
   to machine accuracy is achieved.
\item
   \name{iter}  -- exploits a recursion relation, which can be derived by
               partial integration: 
     
   \be
   F(T,m) = \frac{2T F(T,m+1) + exp(-T)}{2m+1} \label{f-recursion}
   \ee

   Furthermore, we have through direct evaluation of the integral
   for large m:
   $$
     F(T,m) \sim \int_0^1 t^(2m) exp(-T) dt
            = exp(-T) / (2m+1)
   $$
   which yields the starting point for the integration. We begin with
   this approximant for a very large mstart and iterate eqn.~\ref{f-recursion}    until the desired value of m is obtained.

\item \name{fastcomp} -- uses iter, with increasing values of mstart, 
                  until F(T,m) is converged to machine accuracy.
 
\item \name{FTaylor(T,m)} -- computes F by taylor expansion to sixth order
   as described in OS.
   \be
      FTaylor(T,m) = \sum_{k=0}^6 F(TT(T),m+k) (TT(T)-T)^k / k!
   \ee         
   where TT(T) is a set of fixed points spaced DT apart. We have
   chosen $DT = 0.05$ as have OS and Shavitt. The error is then
   bounded by:
   \be
      (TT(T) - T)^7  F(TT(t),m+7) / 7 !
   \ee
   which gives an approximate relative error of $(DT/2)^7/7! = 1.2 \times 
   10^{-15}$, roughly equal to machine accuracy. In order to evaluate
   this function we must tabulate the values of $F$ for the gridpoints
   up to some maximal, $m$ dependent value. These values are stored in the 
   structure \name{ftable}. The contents of this array of vectors is read from 
   the FTABLE file, which must be available in the working directory. It is
   read by the \name{read\_ftable} function. If the file FTABLE is not found 
   the \name{compute\_ftable} function is called to compute the
   values. For each value of $m$ we define a value ${\rm TMAX}(m)$
   which gives the maximal value of $m$ for which the taylor-expansion
   is used to evaluate the function. Otherwise, an approxmate function
   for lrge $T$ is used, whic his described in the next item.
\item   
    \name{highf(T,m)} -- computes the approximate form of f for large
                    T. see OS, formula  (54)
\be
         highf(T,m) = (2m-1) !! / 2 (2T)^m sqrt(pi/T)
\ee
   where $k !! = 1 3 5 \dots k$ for odd k. These values are stored in
   the static array \name{binbin}. This formula is used for $T > TMAX(m)$.
   We choose TMAX(m) as the smallest integer such that the ratio
   of the approximant and the exact value is 1.0 to machine accuracy.
\end{itemize}
   
\section{Tables and Auxiliary Functions}   
*/
#include "io_incl.h"
#include <stdlib.h>
#include <math.h>
#include <numer.h>
#include <errorip.h>
#include "typesip.h"
#include "vtype.h"
#include "qc_utilities.h"  
#include "const.h"
#include "f.h"

ARRAY<Vec*> ftable(MAXM);

void gnufixxxxxxx(){
  ARRAY<Vec*> a(9,9);
  a.print(cout);
}

/* binbin contains the values of k !! for k = 0, ..., 99 */

double* compute_binbin(const int maxm)
{
  int     i;
  double* tmp = new double[maxm];
  double  bb = 1.0;
  double  fac = 1.0;

  for(i=0; i < maxm;i++)
  {
    tmp[i] = bb;
    bb  *= fac;
    fac += 2;
  }
  return tmp;
}

double* binbin = compute_binbin(MAXM);
/*
static int TMAX[] = 
{   28,    35,    35,    41,    44,    44,    48,
    52,    53,    53,    55,    62,    62,    60,
    54,    50,    46,    44,    42,    40,    39,
    38,    37,    36,    35,    35,    34,    34,
    34,    33,    33,    33,    33,    33,    33,
    32,    32,    32,    32,    32,    32 }
*/
static int TMAX[] = 
{   33,    37,    41,    43,    46,    49,    51,
    54,    56,    58,    61,    63,    66,    68,
    70,    72,    74,    76,    78,    80,    82,
    84,    86,    88,    90,    92,    94,    96,
   100,   100,   100,   100,   100,   100,   100,
   100,   100,   100,   100,   100,   100,   100 
};

double highf(double T,int m)
{
  return 0.5*binbin[m] / pow(2*T,m) * sqrt(pi/T);
}

inline double g(double t,double T,int m)
{
  double tt = t*t;
  return pow(tt,m)*exp(-T*tt);
}
/*TEX
\section{Implementation of \name{numint}}
*/  
  
double numint(double T,int m)
{
  double eps   = 1e-4;
  double dt    = 1.0;
  double t;
  double sum = 0.5*dt*( g(0.0,T,m) + g(1.0,T,m)) ;
  double old;
  double tt;  
  do
  {
    double tmp = 0;
    for(t=dt/2;t < 1.0; t += dt)
    {
      tt   = t*t;
      tmp += pow(tt,m)*exp(-T*tt);
    }  
    old = sum;
    sum = 0.5*(sum + dt*tmp);
    dt /= 2;
  } while(dt > eps && fabs(sum-old) > 1e-10);
  if (dt <= eps)
	  cerr << "Convergence error in NUMINT: " << fabs(sum-old) << "\n"; 
  return sum;
}
/*TEX
\section{Implementation of \name{iter and fastcomp}} 
*/  

double iter(double T,int m,int mmax)
{
  int mm = mmax;
  double ff = exp(-T);
  double f = ff/(2*mm+1);
  do
  {
    mm--;
    f = (2*T*f+ff) / (2*mm+1);
  } while(mm > m);
  return f;
}

double fastcomp(IntegType T,int m)
{
  double old;
  int    order = 10;
  double f = iter(T,m,order);
  
  do
  {
    old = f;
    order+=2;
    f = iter(T,m,order);
  } while(fabs(f-old) > 0);
  return f;
}
/*TEX
\subsection{Computation of Lookup Tables}

The following function tabulates the values of $F(T,m)$ in intervals 
from [0,TMAX(m)] in steps of $DT$. The result is written to a file
with name \name{filename}.  
*/  

void compute_ftable(String& filename)
{
  int i;
  int order;

  char *fn = filename();
  fstream fdata(fn,ios::out);
  delete  fn;
  
  
  fdata << "int ftable_maxm   = " << MAXM << ";" <<  endl;
  fdata << "int ftable_size[] = {";

  IVec dim(MAXM);

  for(order = 0; order < MAXM; order++)
  {
    // the number of required elements derives from the maximum of the 7 lower
    // arrays
    int mx = 0;  
    for(i=max(0,order-7); i <= order; i++)
      mx = max(TMAX[order],mx);
    int no = (int) (mx / DT) + 2;
    dim[order] = no;
    fdata << no;
    if (order < MAXM-1)
      fdata << ",";    
    else
      fdata << "};" << endl;
  }
  
  fdata << "double ftable_values[] = {";
  int count = 0;
  char buf[100];
  
  for(order = 0; order < MAXM; order++)
  {
    int no = dim(order);
    (ostream&) cout << "FTABLE::gen: " << order << " " << no << endl;
    for(i=0; i < no; i++)
    {
      sprintf(buf," %20.15e",fastcomp(DT*i,order));
      fdata << buf;
      if (i < no - 1 || order < MAXM-1)
      {
	if (count++ % 8 == 0)
	  fdata << endl;
	fdata << ",";
      }
      else
	fdata << "}\n";
    }
  }
  fdata << "};" << endl;
  
}

/*TEX
\subsection{Reading the ftable from file} 
*/


void read_ftable(String& filename)
{
  char *fn = filename();
  fstream fdata(fn,ios::in);
  if (!fdata)   // file not found --- make ftable
  {
    compute_ftable(filename);
    fdata.open(fn,ios::in);
    if (!fdata)   // death 
    {
      cerr << "Can neither find nor make FTABLE ! " << endl;
      abort();
    }
  }
  
  delete  fn;

  for(int i=0; i< MAXM;i++)
  {
    ftable[i] -> read(fdata);
  }  
}

#include "ftable.h"

void init_ftable()
{
  int order;
  int count = 0;
  for(order = 0; order < ftable_maxm; order++)
  {
    int no = ftable_size[order];
    ftable[order] = new Vec(ftable_size[order]);
    for(int i=0; i < no; i++)
      (*ftable[order])[i] = ftable_values[count++];
  }  
}


/*TEX
\section{Evaluation of F using Lookup Tables\label{section-f-eval}}

The following section implements the evaluation of F using 
either taylor expansion or the large-x expansion for a single 
value.

*/  

const double f2 = 1 /   2.0;
const double f3 = 1 /   6.0;
const double f4 = 1 /  24.0;
const double f5 = 1 / 120.0;
const double f6 = 1 / 720.0;  

Vec   ftmp1;
Vec   ftmp2;
IVec  itmp1;
IVec  itmp2;
Vec   dl;
IVec  idx;

void clean_ftable()
{
  for(int order = 0; order < ftable.size(); order++)
    delete ftable[order];
  ftable.reset(0);
  ftmp1.reset(0);
  ftmp2.reset(0);
  itmp1.reset(0);
  itmp2.reset(0);
  dl.reset(0);
  idx.reset(0);

}

double F(IntegType T,int m)
{
  if (T > TMAX[m])
    return  0.5*binbin[m] / pow(2*T,m) * sqrt(pi/T);    
  else
  {
     int tp = (int) (T / DT + 0.5);
     double dt  = T-tp*DT;
     /* 
      (ostream&) cout << tp << " " << dt << endl;

     (ostream&) cout << "FT: " << ftable[m][tp] << " " << ftable[m+1][tp] << " " << 
       ftable[m+2][tp] << " " <<  ftable[m+3][tp] << " " << endl;
     */
     return (((((  f6* (*ftable[m+6])[tp] *dt 
		 - f5* (*ftable[m+5])[tp])*dt
		 + f4* (*ftable[m+4])[tp])*dt
	         - f3* (*ftable[m+3])[tp])*dt
	         + f2* (*ftable[m+2])[tp])*dt
	         - (*ftable[m+1])[tp])*dt
	         + (*ftable[m])[tp];
     
  }
}


/*inline void gather(Vec& target,Vec& src,IVec& idx,int n)
{
  for(int i=0; i< n; i++)
    target[i] = src[idx[i]];
}
*/
inline void scatter(IntegType* target,IVec& idx,Vec& src,int n)
{
  for(int i=0; i< n; i++)
    target[idx[i]] = src[i];
}

/*inline void gather_multiply(Vec& target,Vec& src,IVec& idx,
		     Vec& factor,double f,int n)
{
  for(int i=0; i< n; i++)
    target[i] = f*target[i]*factor[i] + src[idx[i]];
}
*/
inline void power(Vec& target,Vec& factor,int n)
{
  for(int i=0; i< n; i++)
    target[i] *= factor[i];
}

inline void inv(Vec& target,int n)
{
  for(int i=0; i< n; i++)
    target[i] = 1.0 / target[i];
}

inline void eval_f(double fac,Vec& target,Vec& src,int n)
{
  for(int i=0; i< n; i++)
    target[i] *= fac*sqrt(pi*src[i]);
}

void F(IntegType* T,IntegType* result,int no,int m)
{
  timer -> start(t_f);
  if (ftmp1.size() < no)
  {
    ftmp1.reset(no);
    ftmp2.reset(no);
    itmp1.reset(no);
    itmp2.reset(no);
    dl.reset(no);    
    idx.reset(no);
  }
  
  if (m > MAXM - 7)
    error("F-Evaluation: MAXM too small  -- recompile with larget MAXM");

  int count1 = 0;
  int count2 = 0;
  
  int i;
  int tmax = TMAX[m];
  for(i=0; i < no; i++)
    if(T[i] < tmax)
    {
      ftmp1[count1]    = T[i];
      itmp1[count1++]  = i;
    }
    else
    {
      ftmp2[count2]    = T[i];
      itmp2[count2++]  = i;
    }

  // now evaluate taylor series 
    
  int tmp;
  double dtinv = 1/DT;

  for(i=0; i<count1;i++)
  {  
    idx[i] = tmp = (int) (dtinv*ftmp1[i] + 0.5);
    dl[i] = ftmp1[i] - tmp*DT;
  }
  
  Vec& v0 = *ftable[  m];
  Vec& v1 = *ftable[m+1];
  Vec& v2 = *ftable[m+2];
  Vec& v3 = *ftable[m+3];
  Vec& v4 = *ftable[m+4];
  Vec& v5 = *ftable[m+5];
  Vec& v6 = *ftable[m+6];

  int tp;
  double dt;
  
  for(i = 0; i < count1; i++)
  {
    tp = idx[i];
    dt = dl[i];
    ftmp1[i] = (((((f6*v6[tp] *dt - f5*v5[tp])*dt + f4*v4[tp])*dt
		  - f3*v3[tp])*dt + f2*v2[tp])*dt - v1[tp])*dt + v0[tp];
  }  
  scatter(result,itmp1,ftmp1,count1);

  // now evaluate approximate formula

  double fac = binbin[m] / pow(2.0,m+1);
  
  inv(ftmp2,count2);
  dl = ftmp2;
  if (m == 0)
    ftmp2.set(1);
  else
  {    
    for(i=0; i < m-1;i++)
      power(ftmp2,dl,count2);
  }
  
  eval_f(fac,ftmp2,dl,count2);    
  scatter(result,itmp2,ftmp2,count2);  
  timer -> stop(t_f);
}

/*TEX
\section{Testing f}
*/  

void f_test()
{
  int order;
 
  // compute TMAX(m) 
  /*
  double res;
  double res2;
  double x;
  x = 10.0;  
  for(order=0; order < MAXM; order++)
  {
    for(x = 10; x < 500; x++)
    {
      double ex  = fastcomp(x,order);
      double app = highf(x,order);
      cout << order << SP << x << SP << ex << SP << app << SP << (ex/app)-1.0 << NL; 
      if (app/ex - 1.0 < 0 || ex < 1E-15) 
      {
        fprintf(file," %4.0f, ",x);
	break;
      }
    }
  }
  */

  const int sz = 1000;
  IntegType x[sz];
  IntegType res[sz];
  
  for(int i=0;i<sz;i++)
    x[i] = rand()*10000;
  
  for(order=0; order < MAXM-7; order++)
  {
    cout << "testing order: " << order << endl;
    F(x,res,sz,order);
    for(int i=0;i<sz;i++)
    {      
      double v1 = res[i]; // F(x[i],order);
      double v2 = (x[i] > TMAX[order]) ? 
	highf(x[i],order) : fastcomp(x[i],order);
      
      if (fabs(v1-v2) > 1e-10)
      {
	(ostream&) cout << " ERR: " << i << " " << x[i] << " " << order <<  " " << v1
	  << " " << v2 << endl;
      }
    }
  }
  
  cout << "done. " << endl;
    
}
