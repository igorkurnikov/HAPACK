/*TEX
%  Dortmund C++ Class and Template Library 
%  Copyright (C) 1994 Wolfgang Wenzel and others

%  The code in this file was derived from routines published in 
%  Numerical Recipes, it its the responsibility of the user to insure that
%  he is authorized to use them. Please read the appropriate section entitled
%  LICENSE in the documentation. 

%  This library is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%  Library General Public License for more details.

%  You should have received a copy of the GNU Library General Public
%  License along with this library; if not, write to the Free
%  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%
% $Id: tfstack.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Fixed Length Stacks}

The following class implements a fixed-length stack of elements in
{\em static}storage. Such classes are used when a large number of
relatively small arrays consisting of small objects need to be
manipulated.  This class is inefficient, if
\begin{itemize}
\item The base elements are large and difficult to manipulate, copy etc.
      In this case a pointer-based implementation would be more
      appropriate.
\item If many such stacks, with a greatly varying number of elements
      are  required, i.e. when no sensibly upper bound for the number 
      of elements can be found. Again a pointer-based implementation with a
      flexible number of elements is required.
\end{itemize}
The chief advantage is that no dynamic memory management is required in the
creation/destruction of the class, all such operations are carried out
on the stack and are much faster. The class presently assumes a
defined ordering for the base class.

\subsection{Concepts}

 

*/
#ifndef FSTACK_H
#define FSTACK_H

const int SSTACK_MAX_ELEM = 100;

#include <assert.h>

template <class Base,int MaxSize> class FStack
{
protected:  
     int     cno;
     Base    base[MaxSize];
public:  
     FStack() { cno = 0; }
    ~FStack() {}
int  add(const Base& newel)
     { 
       assert(cno < MaxSize);
       base[cno++] = newel;
       return cno-1;
     }
const FStack<Base,MaxSize>& operator=(const FStack<Base,MaxSize>& src)
     {
       cno = src.cno;
       int i;
       for(i = 0; i < cno; i++)
	 base[i] = src.base[i];
       return *this;
     }
int  sort_in(const Base& newel)
     { 
       int p,i_base;
       for(p=0; p < cno; p++)
	 if (base[p] >= newel) break;
       if (p >= MaxSize)
	 return -1;
       if (cno > 0 && p < cno && base[p] == newel) return -1;
       // shift all elements larger than p
       int mm = (cno >= MaxSize) ? MaxSize - 1 : cno;
       for(i_base = mm; i_base > p; i_base--)
	 base[i_base] = base[i_base-1];
       base[p] = newel;
       if (cno < MaxSize) cno++;
       return p;
     }
Base remove()
     {
       assert(cno > 0);
       return base[--cno];
     }
Base remove(const int pos)
     {
       assert(cno > 0);
       assert(pos < cno);
       Base tmp = base[pos];
       int i;
       for(i = pos; i < cno-1;i++)
	 base[i] = base[i+1];
       cno--;
       return tmp;
     }
void reset()        { cno = 0;    }     
int  size()  const  { return cno; }      
const Base&  operator()(const int i) const 
                       { assert( i>=0 && i < cno); return base[i]; }
      Base& operator[](const int i) 
                    { assert( i>=0 && i < cno); return base[i]; }
const Base& get(const int i) const 
                    { assert( i>=0 && i < cno); return base[i]; }
      Base* get_base()    { return base;  }
int   compare(const FStack<Base,MaxSize> b) const
      {
	if (cno != b.cno)
	  return cno - b.cno;
	else
	  return memcmp(base,b.base,cno*sizeof(Base));
      }
};


template <class Base,int MaxSize> 
ostream& operator<<(ostream& os,const FStack<Base,MaxSize>& stack)
{
  int i;
  for(i=0; i < stack.size();i++)
    os << stack(i);
}

template <class Base,int MaxSize> 
istream& operator>>(istream& is,FStack<Base,MaxSize>& stack)
{
  Base tmp;
  char c;
  is >> c;
  assert(c == '(');

  is >> c;
  while(c != ')')
  {
    is.putback(c);
    is >> tmp;
    if(!is) return is;
    stack.add(tmp);
    is >> c;
  }

  return is;
}


#endif



