/*
   A reallysimple integral program

*/
#include "io_incl.h"
#include <math.h>
#include "const.h"
#include "typesip.h"
#include "function.h"
#include "f.h"

double zeta(RadialFunction& f1,RadialFunction& f2)
{
  return f1.exp()+f2.exp();
}

double rho(RadialFunction& f1,RadialFunction& f2,
	   RadialFunction& f3,RadialFunction& f4)
{
  return zeta(f1,f2)*zeta(f3,f4)/(zeta(f1,f2)+zeta(f3,f4));
}


Location P(RadialFunction& f1,RadialFunction& f2)
{
  double sz = 1.0/zeta(f1,f2);
  return Location(sz*(f1.exp()*f1(X)+f2.exp()*f2(X)),
		  sz*(f1.exp()*f1(Y)+f2.exp()*f2(Y)),
		  sz*(f1.exp()*f1(Z)+f2.exp()*f2(Z)));
}

Location W(RadialFunction& f1,RadialFunction& f2,
	   RadialFunction& f3,RadialFunction& f4)
{
  double z12 = zeta(f1,f2);
  double z34 = zeta(f3,f4);
  Location p(P(f1,f2));
  Location q(P(f3,f4));
  double sz = 1.0/(z12+z34);
  
  return Location(sz*(z12*p(X)+z34*q(X)),
		  sz*(z12*p(Y)+z34*q(Y)),
		  sz*(z12*p(Z)+z34*q(Z)));
}

double distance(Location& f1,Location& f2)
{
  return pow(f1(X)-f2(X),2) + pow(f1(Y)-f2(Y),2) + pow(f1(Z)-f2(Z),2);
}

double ovlp(RadialFunction& f1,RadialFunction& f2)
{
  double zs   = zeta(f1,f2);
  double xi   = f1.exp()*f2.exp()/zs;
  double dist = distance(f1,f2);
  return pow(pi/zs,1.5)*exp(-xi*dist);
}

double ss_integ(RadialFunction& f1,RadialFunction& f2,RadialFunction& f3,
		RadialFunction& f4,int m)
{
  Location p = P(f1,f2);
  Location q = P(f3,f4);
  double rr  = rho(f1,f2,f3,f4);
  double dd  = distance(p,q);
  
  return 2 * sqrt(rr/pi) * ovlp(f1,f2) * ovlp(f3,f4) * fastcomp(rr*dd,m);
}

class Parms
{
public:  
  Momentum        m[4];
  RadialFunction  f[4];
  Parms() {} 
};


double twoel(Parms& p,int m)
{
  
  Parms term1 = p;
  Parms term2 = p;
  Parms term3 = p;
  Parms term4 = p;  
  Parms term5 = p;
  
  // reduce first nonzero index 
  int index = -1;
  for(int i=0; i< 3;i++)
    if (p.m[i].l() > 0) index = i; 
  if (index < 0)
  {
    double res = ss_integ(p.f[0],p.f[1],p.f[2],p.f[3],m);
    cout << "  Twoel: SSSS aux = " << m <<  endl;
    cout << "  Result: " << res << endl << endl; 
    return res;
  }

  
  Direction dir = p.m[index].decrease();
  
  
  term1.m[index].reduce(dir); // term1 is (a   bcd)
  term2.m[index].reduce(dir); // term2 is (a-1 bcd)
  term2.m[index].reduce(dir); 

  int bpos;
  if (index < 2) 
    bpos = 1 - index;
  else
    bpos = 3 - index + 2;   //index=3 then: 3-3+2 = 2 index = 2: 3-2+2 = 3 
   
  term3.m[index].reduce(dir);
  term3.m[bpos] .reduce(dir);
  
  int pos3 = (index + 2) % 4;

  term4.m[index].reduce(dir);
  term4.m[pos3] .reduce(dir);

  int bpos2;
  if (pos3 < 2) 
    bpos2 = 1 - pos3;
  else
    bpos2 = 3 - pos3 + 2;  

  term5.m[index].reduce(dir);
  term5.m[bpos2] .reduce(dir);
  
  // now computing the terms:

  Location pp  = P(p.f[index],p.f[bpos]);
  Location w   = W  (p.f[0],p.f[1],p.f[2],p.f[3]);

  double   z    = zeta(p.f[index],p.f[bpos]);
  double   nu   = zeta(p.f[pos3],p.f[bpos2]);
  double   rho  = z*nu/(z+nu);
  
  double t1 = (pp(dir) - p.f[index](dir)) * twoel(term1,m);
  double t2 = (w(dir) - pp(dir)) * twoel(term1,m+1);

  double t3 = 0;
  double t4 = 0;  
  if (term2.m[index].l() >=0)
  {
    double fac = 0.5*(p.m[index](dir)-1)/z;
    t3 =  fac*twoel(term2,m);
    t4 = -fac*rho*twoel(term2,m+1)/z;
  }

  double t5 = 0;
  double t6 = 0;  
  if (term3.m[bpos].l() >=0)
  {
    double fac = p.m[bpos](dir)/2/z;
    t5 =  fac*twoel(term3,m);
    t6 = -fac*rho*twoel(term3,m+1)/z;
  }

  double t7 = 0;
  if (term4.m[pos3].l() >=0)
    t7 = 0.5*p.m[pos3](dir)*twoel(term4,m+1)/(z+nu);

  double t8 = 0;
  if (term5.m[bpos2].l() >=0)
    t8 = 0.5*p.m[bpos2](dir)*twoel(term5,m+1)/(z+nu);

  double result = t1+t2+t3+t4+t5+t6+t7+t8;

  cout << "Twoel for: " << p.m[0] << p.m[1] << p.m[2] << p.m[3] << " Aux: " 
       << m << endl; 
  cout << "Reducing on index : " << index << " Direction: " << dir << endl;
  cout << "Required Terms: " << endl;
  cout << "  term1: " << term1.m[0] << term1.m[1] << term1.m[2] << term1.m[3] 
    << endl;   
  cout << "  term2: " << term2.m[0] << term2.m[1] << term2.m[2] << term2.m[3] 
    << endl;   
  cout << "  term3: " << term3.m[0] << term3.m[1] << term3.m[2] << term3.m[3] 
    << endl;   
  cout << "  term4: " << term4.m[0] << term4.m[1] << term4.m[2] << term4.m[3] 
    << endl;   
  cout << "  term5: " << term5.m[0] << term5.m[1] << term5.m[2] << term5.m[3] 
    << endl;     
  cout << "RESULT: " << result << endl << endl;
  return result;
}

main()
{
  Parms p;

  for(int i =0 ; i< 4 ;i++)
  {
    cin >> p.f[i] >> p.m[i];
    cout << "RadialFunction : " << i << p.f[i] << " L: " << p.m[i] << endl;
  }
    
  twoel(p,0);
}
