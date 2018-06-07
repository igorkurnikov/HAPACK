#ifndef SSTACK_H
#define SSTACK_H

const int USE_DEFAULT = -1;

#include <assert.h>
#include <tarray.h>

template <class Base,int max> class SStack : public Array<Base>
{
protected:  
     int cno;
public:  
     SStack(const int max_elem = -1)  
       : Array<Base>((max_elem == USE_DEFAULT)? max : max_elem)
     { cno = 0; }
    ~SStack() {}
//compGcc4
using Array<Base>::base;
using Array<Base>::error;
int  add(const Base& newel)
     { 
       assert(cno < Array<Base>::size());
       base[cno++] = newel;
       return cno-1;
     }
int  push(const Base& newel)
     { 
       assert(cno < Array<Base>::size());
       base[cno++] = newel;
       return cno-1;
     }
int  add(SStack<Base,max>& stack_1)
     { 
       int r = 0;
       for(int i =  0; i < stack_1.size(); i++)
	 r = add(stack_1[i]);
       return r;
     }
Base remove()
     {
       assert(cno > 0);
       return base[--cno];
     }
Base pop()
     {
       assert(cno > 0);
       return base[--cno];
     }
void reset() { cno = 0; } 
void resize(int newno) { Array<Base>::reset(newno); cno = 0; } 
int  size()   const  { return cno; }      
void write(FILE* file);
void read (FILE* file);
int  max_size() const { return Array<Base>::size(); }  
};

template <class Base,int max> ostream& operator<<(ostream& os,
					  const SStack<Base,max>& sstack)
{
  for(int i=0; i < sstack.size();i++)
    os << sstack(i) << endl;
  return os;
}

template <class Base,int max> void SStack<Base,max>::write(FILE* file)
{
  int code = fwrite(&cno,sizeof(int),1,file);
  if (code != 1)
    error("IOERROR in Stack::write. ");
  Array<Base>::write(file);
}

template <class Base,int max> void SStack<Base,max>::read(FILE* file)
{
  int code = fread(&cno,sizeof(int),1,file);
  if (code != 1)
    error("IOERROR in Stack::read. ");
  Array<Base>::read(file);
}


template <class Base,int max> class D_SStack : public SStack<Base,max>
{
public:  
  D_SStack(const int max_elem = -1)  : SStack<Base,max>(max_elem) {} 
//compGcc4
using SStack<Base,max>::cno;
using SStack<Base,max>::error;
using SStack<Base,max>::pbase;
Base&       operator[](const int n)
            {
	      if (n<0 || n >= cno) error("Illegal Index",n);
	      return pbase[n];
	    }

const Base& operator()(const int n) const
            {
	      if (n<0 || n >= cno) error("Illegal Index",n);
	      return pbase[n];
	    }

const Base& get(const int n) const
            {
	      if (n<0 || n >= cno) error("Illegal Index",n);
	      return pbase[n];
	    }
};

#endif



