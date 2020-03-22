#include <unistd.h>
#include "tnumarray.h"

main()
{
  const int laenge = 6;  
  int i,offset;

/*
  Array<double> xx(10),yy;
  xx[2] = 2;
  cout << xx;
*/
  NumArray<double> a(laenge),b(laenge,1),c(laenge,10),d;
  for (i=0;i < laenge;i++) 
  {
    a[i]=1.0*i;
    b[i+1]=10.0*i;
    c[i+10]=0.1*i;
  }

  cout<<" ++++++ Testing Vector: \n";
  d=a;
  cout<<"a 0-6  b 0-60  c ''test'' 0-0.6  d=a \n";
  cout<<a;

  IndexArray idx(laenge);
  idx[0] = 2;
  idx[1] = 1;
  idx[2] = 0;
  idx[3] = 3;
  idx[4] = 4;
  idx[5] = 5;
  
  d.gather(a,idx);

  double xval = a*d;

  a *= 3;
  d /= 3;
  a += 5;
  d -= 2;

  a += d;
  a -= d;


  D_NumArray<double> dna(a);

  NumArray2D<double> m1(3,5),m2(3,4),m3(4,5),m3add(4,5),m4(4,4);

  m1  = m2*m3;
  m1 *= 2;
  m3 += m3add;
  m3 += 2;
  m1.multiply(m2,m3);
  m4.transpose();

}
