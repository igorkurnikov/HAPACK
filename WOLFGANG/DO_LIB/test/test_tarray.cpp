#include <unistd.h>
#include <tarray.h>

main()
{
  const int laenge = 6;  
  int i,offset;
  Array<double> a(laenge),b(laenge,1),c(laenge,10),d;
  for (i=0;i < laenge;i++) 
  {
    a[i]=1.0*i;
    b[i+1]=10.0*i;
    c[i+10]=0.1*i;
  }

  cout<<" ++++++ Testing Vector: \n";
  d=a;
  cout<<"a 0-6  b 0-60  c ''test'' 0-0.6  d=a \n";
  cout<<a<<b<<c<<d;

  cout<<"a(0): "<<a(0)<<" a(3): "<<a(3)<<" a(laenge-1): "<<a(laenge-1)<<"\n";
  cout<<"a[0]: "<<a[0]<<" a[3]: "<<a[3]<<" a[laenge-1]: "<<a[laenge-1]<<"\n";
  a[0]=4;
  a[3]=5;
  a[laenge-1]=6;
  cout << " a[0]=4: "<<a[0]
       << " a[3]=5: "<<a[3]
       <<" a[laenge-1]=6: "<<a[laenge-1]<<"\n";

  cout << "pausing 10 seconds ... " << endl;
  sleep(10);

  cout<<" ++++++ Testing: Reset and Rebase\n";
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
 
  cout << "pausing 10 seconds ... " << endl;
  sleep(10);

  cout<< " +++++ Testing Assignment & Output\n";
  a=c;
  d=c;

  cout<<"a d c " << endl;
  cout<<a<<d<<c;

  cout<<"Format of output \n";
  cout<<d;

  cout<< " +++++ Testing Binary IO\n";

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

  d.write("datei");
  b.read("datei");

  cout<<"d.write(''datei'') b.read(''datei''):\n";
  cout<<" d(15-21) b(15-21) \n"<<d<<b;

  cout<< " +++++ Testing Debugging template, provoking error \n";

  // this is easy inherited class does explicit conversion
  D_Array<double> da(a);  
  cout << da;
  
  //this is trickier: here da is converted into an Arra<double>& 
  // to extract the base-class

  Array<double>  x(da);

  // this is downright wild: da is mare into an Array, which is used 
  // in the constructor of db, to finallyy yield a D_Array again.

  D_Array<double> db(da);

  cout << db;

  // this will kill beacause iit is an illegal index.

  double u = da[-3];

}
