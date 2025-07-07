#include <iostream.h>
#include "bitvector.h"
int main()
{
  Bitvector<2> b,a;
  int i,j;
  while (1)
  {
    cout<<"i's fuer a (ende mit -1) :";
    do
    {
      cin>>i;
      if (i!=-1) a.set(i);
    }
    while (i!=-1);
    cout<<"i's fuer b (ende mit -1) :";
    do
    {
      cin>>i;
      if (i!=-1) b.set(i);
    }
    while (i!=-1);
    cout<<a<<b;
    if (a>b) cout<<"a>b \n";
    if (a<b) cout<<"a<b \n";
    if (a==b) cout<<"a==b \n";
    if (a!=b) cout<<"a!=b \n";
    cout<<"number of set bits in a "<<a.num_set()<<"\n";
    cout<<"bitlist for a \n";
    for(i=0;i<a.set_list().size();i++)
      (ostream &) cout<<a.set_list()(i)<<"\n";
    cout<<"number of set bits in b "<<b.num_set()<<"\n";
    cout<<"bitlist for b \n";
    for(i=0;i<b.set_list().size();i++)
      (ostream &) cout<<b.set_list()(i)<<"\n";
    cout<<"a:number of bits left from 0 ... \n";
    for (i=0;i<a.size();i++) 
      (ostream &)cout<<a.num_set_left_from(i)<<" ";
    cout<<"\n";
    cout<<"b:number of bits left from 0 ... \n";
    for (i=0;i<b.size();i++) 
      (ostream &)cout<<b.num_set_left_from(i)<<" ";
    cout<<"\n";
    
    a.reset(); 
    b.reset();
  }
  
  return 0;
};
  
