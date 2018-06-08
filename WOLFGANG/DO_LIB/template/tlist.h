/*TEX
\section{A Simple List Class for small elements}
*/
#ifndef TSLIST_H
#define TSLIST_H
#include <errorip.h>
#include <stdio.h>

template <class Base> class SListElem
{
       Base elem;
       SListElem<Base>* nxt_ptr;
public:
            SListElem(const Base& e,SListElem<Base>* nxt = 0) : elem(e) 
            { nxt_ptr = nxt; }
SListElem<Base>*  next() const { return nxt_ptr; }
Base&       operator()()       { return elem; }
const Base& get() const { return elem; }
};


template <class Base> class SListIterator;
template <class Base> class SList
{
  friend class SListIterator<Base>;
           int _size;
           SListElem<Base>* start;
	   
public:
            SList ()                  { start = 0; _size = 0; }
            SList (const SList<Base>& src);
           ~SList ()                  { clear();   }
SList<Base>& operator=(const SList<Base>& src);
int         size() const { return _size; }
int         insert(const Base& elem)  
            { start = new SListElem<Base>(elem,start); return ++_size; }  
int         remove()                  
            { 
	      if (!start)
		return 0;
	      SListElem<Base>* tmp = start;
	      start = tmp -> next();
	      delete tmp;
	      return --_size;
	    }
void        clear() 
            {
	      while(start)
		remove();
	    }
Base&       operator[](const int n) 
            { 
	     if (n >= _size)
	       error("SList::No such Item.");
	     SListElem<Base>* tmp = start;
	     for(int i = 0; i < n; i++)
	       tmp = tmp -> next();
	     return (*tmp)();
	    }
const Base& operator()(const int n = 0)  const
            { 
	      if (n >= _size)
		error("SList::No such Item.");
	      SListElem<Base>* tmp = start;
	      for(int i = 0; i < n; i++)
		tmp = tmp -> next();
	     return tmp -> get();
	    }
void        write(FILE* f);
void        read (FILE* f);
};

template <class Base> class SListIterator
{
           SListElem<Base>* current;
public:
	   SListIterator(const SList<Base>& list) { current = list.start; }
Base&      operator()() 
           { 
	     if (!current) error("SListIterator::No Item.");
	     return (*current)();
	   } 
int        operator++() { if (current) current = current -> next(); return 1; }
int        in_range() const { return (current != 0); }
};

template <class Base> SList<Base>::SList(const SList<Base>& src)
{
  start = 0;
  _size = 0;
  *this = src;
}

template <class Base> SList<Base>& SList<Base>::operator=(const SList<Base>& src)
{
  cout << "slist copy " << endl;
  SListIterator<Base> iter(src);
  while(iter.in_range())
  {
    insert(iter());
    ++iter;
  }
  return *this;
}

template <class Base> void SList<Base>::write(FILE* f)
{
  SListIterator<Base> iter(*this);
  int code = fwrite(&(_size),sizeof(int),1,f);
  if (code != 1)
      error("Error in SList::write.");
  
  while(iter.in_range())
  {
    int code = fwrite(&(iter()),sizeof(Base),1,f);
    if (code != 1)
      error("Error in SList::write.");
    ++iter;
  }
}

template <class Base> void SList<Base>::read(FILE* f)
{
  int sz;
  int code = fread(&sz,sizeof(int),1,f);
  if (code != 1)
      error("Error in SList::write.");
  for(int i = 0; i < sz; i++)
  {
    Base b;
    int code = fread(&b,sizeof(Base),1,f);
    if (code != 1)
      error("Error in SList::write.");
    insert(b);
  }
}

template <class Base> ostream& operator<<(ostream& os,const SList<Base>& sl)
{
  SListIterator<Base> iter(sl);
  while(iter.in_range())
  {
    os << iter() << endl;
    ++iter;
  }
  return os;
}

#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
template <class Base> ostream_withassign& operator<<(ostream_withassign& os, 
						     const SList<Base>& sl)
{
  SListIterator<Base> iter(sl);
  while(iter.in_range())
  {
    os << iter() << endl;
    ++iter;
  }
  return os;
}
#endif
#endif
