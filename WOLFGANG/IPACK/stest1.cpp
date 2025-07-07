#include "io_incl.h"
#include <string.h>
#include <assert.h>
#include "tquadruple.h"
#include "const.h"
#include "typesip.h"

double i(Quadruple<Momentum>& m,Quadruple<double>& exp,int aux)
{                 
  if (m[0].x() < 0) return 0;
  if (m[1].x() < 0) return 0;
  if (m[2].x() < 0) return 0;
  if (m[3].x() < 0) return 0;

  if (m[0].y() < 0) return 0;
  if (m[1].y() < 0) return 0;
  if (m[2].y() < 0) return 0;
  if (m[3].y() < 0) return 0;

  if (m[0].z() < 0) return 0;
  if (m[1].z() < 0) return 0;
  if (m[2].z() < 0) return 0;
  if (m[3].z() < 0) return 0;
  // cout << " ++++ i: " << m << " " << aux << endl;
  
  double zetasum = (exp[0] + exp[1] + exp[2] + exp[3]);
  double rho     = (exp[0] + exp[1]) * (exp[2] + exp[3]) / zetasum; 
  
  if(m[1].l() > 0)
  {
    // HRR on index 2

    Direction d = m[1].decrease();
    Quadruple<Momentum> m1;
    m1[1].reduce(d);
    m1[0].induce(d);
    return i(m,exp,aux);    
  }
  
  if(m[3].l() > 0)
  {
    // HRR on index 2

    Direction d = m[3].decrease();
    Quadruple<Momentum> m1;
    m1[3].reduce(d);
    m1[2].induce(d);
    return i(m,exp,aux);
  }
  
  if(m[2].l() > 0)
  {
    // VRR on index 3
   
    Direction d = m[2].decrease();
    Quadruple<Momentum> m1(m);
    Quadruple<Momentum> m2(m);
    m1[2].reduce(d);
    m1[2].reduce(d);

    m2[2].reduce(d);
    m2[0].reduce(d);
    
    double zs    = (exp[2]+exp[3]);
    double term1 =  0.5 * (m[2](d) - 1) / zs * i(m1,exp,aux);
    double term2 = -0.5 * (m[2](d) - 1) * rho / zs / zs * i(m1,exp,aux+1);
    double term3 =  0.5 *  m[0](d) / zetasum * i(m2,exp,aux+1);
    return term1+term2+term3;
  }

  if(m[0].l() > 1)
  {
    // VRR on index 1
   
    Direction d = m[0].decrease();
    Quadruple<Momentum> m1(m);
    Quadruple<Momentum> m2(m);
    m1[0].reduce(d);
    m1[0].reduce(d);

    m2[0].reduce(d);
    m2[2].reduce(d);
    
    double zs    = (exp[0]+exp[1]);
    double term1 =  0.5 * (m[0](d) - 1) / zs * i(m1,exp,aux);
    double term2 = -0.5 * (m[0](d) - 1) * rho / zs / zs * i(m1,exp,aux+1);
    double term3 =  0.5 *  m[2](d) / zetasum * i(m2,exp,aux+1);
    return term1+term2+term3;
  }  

  if (m[0].l() == 1)
    return 0;
  
  // (SSSS)
  
  double ovlp1  = pow(pi/(exp[0]+exp[1]),1.5);
  double ovlp2  = pow(pi/(exp[2]+exp[3]),1.5);
  
  double sinteg = 2.0 * sqrt(rho/pi) * ovlp1 * ovlp2 / (2*aux + 1);
  return sinteg;
  
}


double inorm(Quadruple<Momentum>& m,Quadruple<double>& exp,int aux)
{
  double n0 = pow(2*exp[0]/pi,0.75) * pow(4*exp[0],m[0].l()/2); 
  double n1 = pow(2*exp[1]/pi,0.75) * pow(4*exp[1],m[1].l()/2); 
  double n2 = pow(2*exp[2]/pi,0.75) * pow(4*exp[2],m[2].l()/2); 
  double n3 = pow(2*exp[3]/pi,0.75) * pow(4*exp[3],m[3].l()/2); 
  cout << "norm: " << n0 * n1 * n2 * n3 << " "; 
  return i(m,exp,aux) * n0 * n1 * n2 * n3;
}

main()
{
  Quadruple<double> exp;
  exp[0] = 1.0;
  exp[1] = 1.0;
  exp[2] = 1.0;
  exp[3] = 1.0;
  
  Quadruple<Momentum> m;

  int l1 = 8;
  int l2 = 8;
  
  for(Momentum m1(l1); m1(); ++m1)
  for(Momentum m2(l2); m2(); ++m2)
  {
    m[0] = m1;
    m[2] = m2;
    for(int aux = 0; aux < 16-l1-l2+1; aux++)
    {
      double xx =  i(m,exp,aux);
      char buf[100];
      sprintf(buf," %3i  %15.6f ",aux,xx);
      cout << " XXX: " << m[0] << " " << m[2] << " " << buf << endl; 
    }
  }
  
}
