#include "stringc.h"
#include <numer.h>
#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#define laenge 7

main() 
{
  int i,offset;
  D_Vector a(laenge),b(laenge,1),c(laenge,10),d;
  for (i=0;i < laenge;i++) {
      a[i]=1.0*i;
      b[i+1]=10.0*i;
      c[i+10]=0.1*i;

  }
  cout<<"Testlauf \n";
  d=a;
  cout<<"a 0-6  b 0-60  c ''test'' 0-0.6  d=a \n";
  cout<<a<<b<<c<<d;

  cout<<"a(0): "<<a(0)<<" a(3): "<<a(3)<<" a(laenge-1): "<<a(laenge-1)<<"\n";
  cout<<"a[0]: "<<a[0]<<" a[3]: "<<a[3]<<" a[laenge-1]: "<<a[laenge-1]<<"\n";
  a[0]=4;
  a[3]=5;
  a[laenge-1]=6;
  cout<<"a[0]=4: "<<a[0]<<" a[3]=5: "<<a[3]<<" a[laenge-1]=6: "<<a[laenge-1]<<"\n";
  cin>>i;
  cout<<"Reset und Rebase\n";
  a.reset(laenge+2,10);
  b.reset(laenge);
  offset=c.rebase(15);
  cout<<"a.reset(''reset_a'',laenge+2,10) b.reset(laenge) \n";
  cout<<"laenge: "<<laenge<<" a.size=laenge+2: "<<a.size();
  cout<<" a.offset=10: "<<a.offset()<<"\n";
  cout<<a;
  b.set(11);
  cout<<"b.set(11) \n"<<b;
  cout<<"c.rebase(15)\n"<<c;
  cout<<"previous offset= 10: " <<offset<<"\n";
  cin>>i;
  cout<<"Assignment\n";
  a=c;
  d=c;
  d+=d;
  cout<<"Arithmetik a=c c d=2*c \n";
  cout<<"A:" << a<< "C:" << c<< "D: " << d<<"\n";
  cin >> i;
  a+=c;
  cout<<"A +=c " << a<<"\n";
  a*=2;
  cout<<"a*=2 " << a<<"\n";
  d+=5;
  cout<<"d+=5 " <<d<<"\n";
  cin>>i;
  d-=1;
  cout<<"d-=1 " << d<<"\n";
  d-=d;
  cout<<"d-=d " << d<<"\n";
  a/=2;
  cout<<"a/=2 " << a<<"\n";
  cin >> i;
  b = a;
//  b=a+c;
//  d=a-c;
  cout<<"a c b=a+c d=a-c \n";
  cout<<a<<c<<b<<d;
  cin>>i;
  double var=a*c;
  cout<<"Sorting and so on \n";
  cout<<"a c b=a*c unveraendert b \n"<<a<<c<<b;
  cout<<"a*c=1.82: "<<a*c<<"  a*b=46.2: "<<a*b<<" var=a*c: "<<var<<"\n";
  cin>>i;
  a[16]=50;
  a[17]=-2;
  b.rebase(0);
  b-=b;
  b[0]=2;
  cout<<"a.maxpos()=16: "<<a.maxpos()<<" a.minpos()=17: "<< a.minpos()<<"\n";
  cout<<"a.max()=50: "<<a.max()<<" a.min()=-2: "<< a.min()<<"\n";
  cout<<"b.maxpos()=0: "<<b.maxpos()<<" b.minpos()=1?: "<<b.minpos()<<"\n";
  cout<<"b.max()=2: "<<b.max()<<" b.min()=0: " <<b.min()<<"\n";
  cin>>i;
  a.sort();
  cout<<"a.sort \n"<<a;
  a.sort(0);
  cout<<"c=a.sort(0) \n"<<a;
  cin>>i;

  cout<<"Gather and Scatter \n";
  Vector a1(2*laenge),a2(4);
  a1.set(0);
  a2.set(0);
  IVector index1(laenge,15),index2(4);
  index1[15]=1;
  index1[16]=2;
  index1[17]=4;
  index1[18]=6;
  index1[19]=8;
  index1[20]=9;
  index1[21]=13;
  index2[0]=16;
  index2[1]=15;
  index2[2]=21;
  index2[3]=19;
  cout<<a;
  a2.gather(a,index2);
  cout<<"gather 16 15 21 19 \n"<<a2;
  cout<<a;
  a1.scatter(a,index1);
  cout<<"scatter 1 2 4 6 8 9 13 \n"<<a1;
  cout<<a;
  cin>>i;

  cout<<"Format of output \n";
  cout<<d;

  cin>>i;

  fstream file("data",ios::out);

  d.write(file);
  file.close();
  file.open("data",ios::in);
  a.set(0);
  a.read(file);
  file.close();
  cout<<"fstream file: d.write(file) a.read(file): \n";
  cout<<" d a : \n";
  cout << d;
  cout << a;
  cin>>i;

  cout<<"b(0-6) vorher \n"<<b;
  d.write("datei");
  b.read("datei");
  cout<<"d.write(''datei'') b.read(''datei''):\n";
  cout<<" d(15-21) b(15-21) \n"<<d<<b;
  cin>>i;

  return 0;
}

