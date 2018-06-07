#include <iostream.h>

/*

   Dear Mr. Abel:

   Unfortunately there is a follow-up problem regarding the ostream
   overloading implementation for the templates. If I implement consistently 
   output operators for BOTH ostream& and ostream_withassign& the problem
   shifts.

   Without the definition of the operator<< for ostream_withassign& for the 
   templates the output for the derived non-template type MyT didnt work.
   
   If I introduce the operator<< for ostream_withassign& this output 
   operations is resolved properly, but the output operations for the 
   derived template class bomb.

   In order to further simplify the problem I have constructed a test
   case which has noting to do with ostreams. 

   In addition to my template Base and Derived I define 
  
        a class  A      -- which stands for ostream and 
        a class  AWITH  -- which stands for ostream_withassign

   The I define a function q (simulating operator<<), which takes
   either an A or an AWITH as a first argument and a Base<int>,
   Derived<int> or MyT as a second argument respectively. We are
   completely independent from the fallacies of C++ IO.

   The compiler now fails at: 

   q(aw,d);

   where aw is an AWITH and d is a Derived<int>. The compiler thinks:

           Call matches "q(AWITH&,Derived<T>
           Call matches "q(AWITH&,Base<T>&)"
           Call matches "q(A&,Derived<T>&)".   

   As far as I can tell it shouldnt since there is ONE exactly matching 
   template defined for this call. Both the second and third possiblity are 
   wrong, since "q" might want to manipulate members of AWITH which
   are not members of A and similarly members of Derived<int> which are 
   not members of Base<int>.

   Cheers
   Wolfgang Wenzel

*/

template <class T> class Base 
{
public:
  int f();   
};

template <class T> class Derived : public Base<T> 
{ 
  int u;
  
public:  
  int f();   
};

class MyT : public Derived<int>
{
public:  
    int f();   
};

template <class T> ostream& operator<<(ostream& os,Base<T>& b)    {return os;}
template <class T> ostream& operator<<(ostream& os,Derived<T>& b) {return os;}
                   ostream& operator<<(ostream& os,MyT& m)        {return os;}

template <class T> 
ostream_withassign& operator<<(ostream_withassign& os,Base<T>& b);
template <class T> 
ostream_withassign& operator<<(ostream_withassign& os,Derived<T>& b);
ostream_withassign& operator<<(ostream_withassign& os,MyT& m);
     
ostream&            g(           ostream& os,MyT& b) { return os; }
ostream_withassign& g(ostream_withassign& os,MyT& b) { return os; }

template <class T> ostream& g(ostream& os,Base<T>& b)    { return os; }
template <class T> ostream& g(ostream& os,Derived<T>& b) { return os; }

template <class T> ostream& g(ostream_withassign& os,Base<T>& b);
template <class T> ostream& g(ostream_withassign& os,Derived<T>& b);

ostream& t(MyT& b)                           {  return cout; }
template <class T> ostream& t(Base<T>& b)    {  return cout; }
template <class T> ostream& t(Derived<T>& b) {  return cout; }

class A
{

};

class AWITH : public A
{

};
     
int                    q(A& a,MyT& b);
template <class T> int q(A& a,Base<T>& b) {}
template <class T> int q(A& a,Derived<T>& b);

                   int q(AWITH& a,MyT& b);
template <class T> int q(AWITH& a,Base<T>& b);
template <class T> int q(AWITH& a,Derived<T>& b);


main()
{
  Base<int> b;
  Derived<int> d;
  MyT          m;  

  ostream& os = cout;
  
  cout << b;
  cout << d; // this call fails now ?
  cout << m; // this is the call which used to fail, but it is ok now 

  os << b;
  os << d;
  os << m;  

  g(cout,b);
  g(cout,d);  // this call fails now
  g(cout,m);  // this is the call which used to fail, but it is ok now 

  t(b);
  t(d);
  t(m);      // this call works ?

  A a;
  AWITH aw;
  
  q(a,b);
  q(a,d);
  q(a,m);  

  q(aw,b);
  q(aw,d) ;  // this call fails !
  q(aw,m);    

}
