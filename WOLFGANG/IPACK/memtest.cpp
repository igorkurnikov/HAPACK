#include "io_incl.h"
#include <stringc.h>
#include "memchk.h"


template <class Base> class Array
{
private:
enum Array_Mode { unmarked, marked, mapped, mapped_marked };
void init();
     Array_Mode mode;
  
protected:   

  int    off;
  int    len;			
  
  Base*  base;
  Base*  pbase;
  
public:
         Array();
         Array(const int new_len, const int new_off = 0);
         Array(const Array<Base>& vec);

virtual ~Array () {  delete [] base;  }
void     reset(const int new_len,const int new_off = 0);
int      size()   const { return len; }
int      bound() const { return (len+off); } 
int      offset() const { return off; }
int      rebase(const int newoff);

void     append(const Array<Base>& vec2);
void     steal (Array<Base> & from);
Base&    operator[](const int n) { return pbase[n]; }
const Base& operator()(const int n) const { return pbase[n]; }
Base*    get_base() const { return base; }  
const Base& get(const int n) const { return pbase[n]; }

Array<Base>& operator=(const Array<Base>& rhs);			

void   read(const String& name);
void   write(const String& name) const;
void   read(fstream& f);
void   write(fstream& f) const;
void   read(FILE* fp);
void   write(FILE* fp) const;
void   set  (const Base value); 

void   error(const String& msg) const ; 
void   error(const String& msg, const int n) const; 				
//mms
//template <class Base> friend ostream& operator<<(ostream& os, const Array<Base>& vec);
void   print(ostream& os) const;
};

/*TEX
\section{Implementation of the Array Template}
 
\subsection{Init} 
  
The key initializations are done in the \name{init} member-function.

*/

template <class Base>  void Array<Base>::init()
{
  len     = 0;
  off     = 0;
  base  = NULL;
  pbase = NULL;
  mode  = unmarked;
}

/*TEX
\subsection{Constructors \& Reset} 

This is fairly straigforward. Note, however that the copy constructors
should use the assigment operator of \name{Base} as opposed to memcpy,
since this might have undesired effects. If speed is an issue the
intantiated classes might explicitely overload the copy contructor and
assigments.  

*/ 

template <class Base> Array<Base>::Array() { init();  }

template <class Base> Array<Base>::Array(const Array<Base>& old)
{
  init(); 
  reset(old.len,old.off);
  int i;
  for(i=0; i<len;i++)
    base[i] = old.base[i];
}

template <class Base> Array<Base>::Array(const int new_len, const int new_off)
{ init();  reset(new_len,new_off); }    


template <class Base>  void Array<Base>::set(const Base value)             
{						
  int i;					
  for( i =0; i < len; i++)			
    base[i] = value;				
}						

template <class Base>  void Array<Base>::steal(Array<Base> & from)             
{						
  reset(0);
  off=from.off;
  len=from.len;
  base=from.base;
  pbase=from.pbase;
  mode=from.mode;
  from.off=0;
  from.len=0;
  from.base=NULL;
  from.pbase=NULL;
  from.mode=unmarked;
}						

template <class Base> 
void Array<Base>::reset(const int new_len,const int new_off) 
{							
  if (mode != unmarked)
    error("Attempt to resize marked vector. ");

  if (len > 0) 
    delete [] base;	

  base  = NULL;						
  if (new_len < 0) 
     error("Size must be greater or equal to zero.");	 
  if (new_len > 0)
  {
     base = new Base[new_len];
     if (!base) 
     {
       cerr << "Dimension: " << new_len << "Offset: " << new_off << endl;
       error("Memory Allocation Error.");		 
     }
  }
  else base = 0L;
  len       = new_len;	
  off       = new_off;
  pbase     = base - off;
}

template <class Base> 
int Array<Base>::rebase(const int new_off)
{
  int tmp = off;
  off = new_off;
  pbase = base - off;
  return tmp;
}

template <class Base> 
Array<Base>& Array<Base>::operator=(const Array<Base>& rhs)   
{							   
  if (rhs.len != len || rhs.off != off) 
    reset(rhs.len,rhs.off);
  int i;						   
  for(i=0; i < len; i++)
    base[i] = rhs.base[i];
  return *this;
}

template <class Base>
void Array<Base>::append(const Array<Base>& vec2)
{
  int i, old_len=len;
  Base *tmp= new Base[len];

  for(i=0; i<len; i++)
    tmp[i]=base[i];
  reset(len+vec2.len, off);
  for(i=0; i<old_len; i++)
    base[i]=tmp[i];
  for(i=old_len; i<len; i++)
    base[i]=vec2.base[i-old_len];
  delete [] tmp;
}

/*TEX
\subsection{Binary IO}

We have presently implemented the binary IO operations via fstreams. 
We are not sure this is wise. The issue is the buffering. If anyone
knows how the binary fstream io-operations are buffered or how to
switch possible buffering off, please let us know. 

If you find these routines hiedously slow than it is probably becaus
of buffering problems. If this is a problem we shall revert to \name{FILE*}.

*/

template <class Base> void Array<Base>::read(const String& file_name)
{
  char* tmp = file_name();
  fstream file(tmp,ios::in);  
  read(file);
  delete tmp;
}

template <class Base> void Array<Base>::read(fstream& file)
{
  if (!file)
    error("Input file not found");
  int ll;
  file.read((char*) &ll,sizeof(len));
  if(!file || file.gcount() != sizeof(len))
    error("Read Error for length.");
  int noff;
  file.read((char*) &noff,sizeof(off));
  if (!file || file.gcount() != sizeof(off)) 
    error("Read Error for offset.");
  if (ll != len || off != noff)
    reset(ll,noff);
  if (len > 0)
  {
    file.read((char*) base,sizeof(Base)*len);
    if (!file || file.gcount() != len*sizeof(Base)) 
      error("Read Error for base.");
  }
}

template <class Base> void Array<Base>::read(FILE* fp)
{
  if ( fp == NULL)
    error("Input filepointer not found");
  int ll;
  int code = fread(&ll,sizeof(int),1,fp);
  if (code != 1)
    error("Read Error for length.");
  int noff;
  code = fread(&noff,sizeof(int),1,fp);
  if (code != 1)
    error("Read Error for offset.");
  if (ll != len || off != noff)
    reset(ll,noff);
  if (len > 0)
  {
    code =  fread( base,sizeof(Base),len,fp);
    if (code != len)
      error("Read Error for base.");
  }
}

template <class Base> void Array<Base>::write(const String& file_name) const
{
  char* tmp = file_name();
  fstream file(tmp,ios::out);
  write(file);
  delete tmp;
  
}

template <class Base> void Array<Base>::write(fstream& file) const 
{
  if (!file)
    error("Open error for output file.");
  file.write((char*) &len,sizeof(len));
  if(!file)
    error("Write Error for length.");
  file.write((char*) &off,sizeof(off));
  if(!file)
    error("Write Error for offset.");
  if (len > 0)
    file.write((char*) base,sizeof(Base)*len);
  if (!file)
    error("Write Error for base.");
}

template <class Base> void Array<Base>::write(FILE* fp) const 
{
  if ( fp == NULL )
    error("No filepointer for output");
  if( fwrite( &len,sizeof(int),1,fp) == 0)
    error("Write Error for length.");
  if( fwrite( &off,sizeof(int),1,fp) == 0)
    error("Write Error for offset.");
  
  if(len > 0)
    if (fwrite( base,sizeof(Base),len,fp) == 0)
      error("Write Error for base.");
}

/*TEX
\subsection{Regular IO and Errors}

In my opinion the iomanip gadgets dont help with this, but I am happy
to be proven wrong.

*/

template <class Base> ostream& operator<<(ostream& os, const Array<Base>& vec)
{
  vec.print(os);
  return os;
}

#ifndef GNU
template <class Base> ostream_withassign& 
operator<<(ostream_withassign& os, const Array<Base>& vec)
{
  vec.print(os);
  return os;
}
#endif

template <class Base> void Array<Base>::print(ostream& os) const
{			

  if ( len == 0 ) { 
    os << "Vector is empty.";
  }
  else {
    os << "  *** Vector *** Lo: " << off << " *** Hi: " << len+off-1 << endl;
    int i;
    char buf[100];
    for(i=0; i < len; i++) {
      os << " " << base[i] << endl;
    }
  }
  os << endl;

}	

template <class Base> void  Array<Base>::error(const String& msg) const 
{
  cerr << "Error in Vector  ::: " << msg << endl;
  abort();
}

template <class Base> void  Array<Base>::error(const String& msg,const int n) 
const 
{
  cerr << "Error in Vector  ::: " << msg << " on Element: " << n << endl;
  abort();
}
 

/*TEX
\section{Template for 1D Array with Access Checking}
  
\subsection{Class Header}

The \name{Array} class provides no access checking on the operators []
amd (). Therefore a derived template \nameindex{D\_Array} is defined, which
implements full access checking. This is done by template class inheritance,
hence \name{D\_Array} has all features of \name{Array}
*/
  
template <class Base> class D_Array : public Array<Base>
{
public:
         D_Array() { }
         D_Array(const int new_len, const int new_off = 0) :
	   Array<Base>(new_len,new_off) {}
         D_Array(const Array<Base>& vec) : Array<Base>(vec) { }
         D_Array(const D_Array<Base>& vec)  : Array<Base>(vec) { }
virtual ~D_Array () { }

D_Array<Base>& operator=(const D_Array<Base>& vec)
                   { Array<Base>::operator=(vec);
                     return *this; }

Base&    operator[](const int n)
{
  if (n<off || n >= len+off) error("Illegal Index",n);
  return pbase[n];
}

const Base& operator()(const int n) const
{
  if (n<off || n >= len+off) error("Illegal Index",n);
  return pbase[n];
}

const Base& get(const int n) const
{
  if (n<off || n >= len+off) error("Illegal Index",n);
  return pbase[n];
}

};

class B
{
  double t[3];
};


main(int argc,char** argv)
{
  memory_manager = new Memory_Manager(100000);
  memory_manager -> on();
 
  if (argc > 1)
  {
    int allocno;
    sscanf(argv[1],"%i",&allocno);
    cout << "intercepting: " << allocno << endl;
    memory_manager -> intercept(allocno);
  }

  Array<int> *a = new Array<int>(3);
  Array<B> *b   = new Array<B>(3);
  
  memory_manager -> print(cout);
  
  delete a;
  delete b;

  memory_manager -> print(cout);
}
