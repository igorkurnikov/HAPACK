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
% $Id: tarray.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\chapter{Generic Arrays}  
\section{One-Dimensional  Class Template \name{Array}}

In the following we will implement a generalized array class, which is
intended for one-dimensional arrays of generic base classes, which are
required to have a default constructor and overloaded stream io
functions.  

\subsection{Basic Structure}

The vectors we use here are dynamic arrays of fixed length. Upon
initialization, we specify both the \name{length} and the offset,
i.e. the minimal index of the vector. The minimal index is transient,
i.e. may be changed with little effort during the life of the object.
The change of the length requires a \name{reset} operation, during
which the current data is typically lost. As an exception from this
rule, the append function allows the concatenation of two
vectors. Optionally a \name{name} can be given to the vector object,
which also serves as the default file-name for binary io
operations. Both access operators () and [] are available, the former
providing merely read capability, the latter providing both read and
write capability. We implement two templates, one without access
checking and one without. The class additionally provides some primitive
error-handling capabilties, but presently no exxception handling. 

\subsection{Constructors and Memory Managment}

The primary idea behind the class is to hide the memory management
from the user and to provide an easy interface for standard operations
as well as input and output. An empty vector is declared by the default
constructor. An unnamed vectors is created by specifying the length and 
optionally the base-offset of the vector. In addition a name for the vector
can be specified. The copy constructor is also explicitly provided.

\subsection{Use \& Access}

Elements of vectors are accessed using the [] or () operators with a
single argument. The [] allows to both obtain the value of the
element, as well as to change it. The () operator is limited to
retrieve values, but cannot alter the contents of the \name{Vector}.
Length and optionally the name of the vector can be adjusted at any
time by the member function \name{reset}, which can either take a new
length and an optional offset, or a name, new length and optional
offset as arguments.  The previous contents of the vector is then
lost. As the change of the offset requires no memory management, this
field is changed by the member function \name{rebase}, which returns
the old offset. This is particularly useful in numerical recipes
routines, which require an vector offset of 1. Present length and
offset are returned by the functions \name{size} and \name{offset}.

Vectors may be assigned to one another, the target vector is
automatically resized if needed. The member function \name{set} allows
to set all values of a vector to a fixed value. 

\subsection{Input and Output}

The output operator has been overloaded. It prints the name of the vector 
and the elements, 5 per line, using a user-defined format. The
format-string is predefined, but may be altered using the \name{set\_format}
member function, which returns the previous format string. Formatting
is accomplished using the \name{sprintf} function, it can be switched off
by setting the format to NULL.

\subsection{Mapping -- NOT IMPLEMENTED}

We have included the flags to allow the mapping of the allocated space of
one vector to another object of the \name{MappedVector} class. This
allows the user to address the same elements in a variety of different
ways without explicit copying. To insure the safety of the associated 
memory, no reset operations may be performed on any mapped vector and no
vector may presently be mapped more than once at a time.

\subsection{Description of Functions}

In the following we will give a brief description of the public member 
functions of the class:
\begin{itemize}
\item{\name{Array$<$Base$>$ vec;} creates a vector of length 0.}
\item{\name{Array$<$Base$>$ vec(10);} creates a vector of length 10
with indices starting at zero.}
\item{\name{Array$<$Base$>$ vec(10,3);} creates a vector of length 10
with indices starting at 3.}
\item{\name{vec.reset(new\_len,new\_off);} resets the 
vector to the arguments given, \name{new\_off} is optional. 
The content of the previous vector is lost!}
\item{\name{vec.steal(from)} steals contents from vector ''from''. ''from'' 
is empty after this operation.(only pointer operation) }
\item{\name{vec.size()} returns the length of the vector,
      \name{offset()} the offset.}
\item{\name{vec.bound()} returns the bound of the vector} 
\item{\name{vec.rebase(new\_off)} resets the offset to the paramter
       and returns the previous offset.  }
\item{\name{vec.get\_base()} returns the value of the base-pointer to the class.
The use of this function is highly discouraged, since it defeats the purpose 
of creating a contained lcass in the first place. However, it is sometimes 
neccessary, when for instance a FORTRAN subroutine is to be called, for which
no specific header exists, or when sub-blocks of data are to be manipulated 
by functions such as \name{memcpy} etc. }
\item{\name{vec.retype(new\_name)} resets the name and returns the 
    previous name. }
\item{\name{vec.set(value)} sets all elements of the vector 
to \name{value}.  }
\item{\name{vec[i]} returns a reference to the i-th element of 
the vector.
\name{vec(i)} is a \name{const} function, returning the value, 
but not a reference.}

\item{ \name{vec1=vec2; } assigns \name{vec2} to \name{vec1}. The 
target inherits both length and offset of the source, but not the name.
The target is resized if neccessary. }
\item{The \name{vec1.append(vec2)} appends the vector \name{vec2} at the end
 of \name{vec1}. The name and offset of \name{vec1} 
 are kept.}
\item{The \name{read} and \name{write} functions implement binary
  input and output, reading/writing length,offset and elements of the 
    vector. For a \name{String} argument e.g. \name{vec.write("String")},
    the file is first opened, the
operation performed, then the file is closed. If the argument is a
\name{fstream} or a \name{FILE* filepointer} only the io-operations are 
performed. }
\item{ The output operator provides a formatted output capability.
The name, the length and offsets are printed, followed by the elements.
The user can influence the printing of the elements in two ways:
  \begin{itemize}
  \item{ The function \name{vec.set\_format(char *new\_format)}
  specifies a format string like " \%6.4f" (blanks must be included), as used
  in printf to format the individual elements.
  The function returns a pointer to the previous format string. The 
   intial formatis NULL. }
  \item{ The function \name{vec.set\_no\_elem(int n)} sets the number 
   of elements which are printed on a single line. This value defaults to 5. }
  \end{itemize} }
\end{itemize}

\subsection{Class Header}

*/
#ifndef TARRAY__H
#define TARRAY__H

#include <stringc.h>
#include "io_incl.h"
#include <stdlib.h>

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

#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
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
//    char buf[100];
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
//compGcc4
using Array<Base>::off;
using Array<Base>::len;
using Array<Base>::pbase;
using Array<Base>::error;

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

/*TEX
\section{Class Array2D}

\subsection{Basic Structure}
In the following we will define a \name{Array2D} class in the spirit of the 
associated \name{Array} class. In order to enable the use of fortran-routines
the matrix-elements are stored in clusters in the same way as fortran
arrays do. Therefore it suffices to pass a pointer to the first element.
Internally the access is speeded up by an pointer-to-pointer construction (as
possible in C++) according to the offsets.

\subsection{Constructors and Memory Managment}

As for the \name{Array} class the intention is to hide the memory management
from the user and to provide an easy interface for standard operations
as well as input and output. An empty matrix is declared by the default
constructor. An unnamed matrix is created by specifying the dimensions
(number of line/number of rows) and
optionally the base-offsets of the matrix. In addition a name for the matrix
can be specified. The copy constructor is also explicitly provided.

\subsection{Use \& Access}

The element access function is the overloaded operator (), which can
be used both for retrieval and assignment. It needs a double argument:
(n,m). The constant function \name{get(n,m)} can be used to obtain
matrix elements, too.

The functions
\name{reset}, \name{rebase}, \name{size1}, \name{size2}, \name{offset1},
\name{offset2}, \name{set}, \name{set\_format}, \name{set\_no\_elem},
\name{read},
\name{write}, the overloaded output operator and the assignment operator are
analogous to the \name{Array} class. 
Please consult the description of the Array class
for detailed information concerning these functions.

Matrices may be assigned to one another, the target matrix is automatically
resized if needed.

\subsection{Input and Output}

The output operator has been overloaded. It prints the name of the matrix
and the elements, 5 per line, using a user-defined format. The
format-string is predefined, but may be altered using the \name{set\_format}
member function, which returns the previous format string. Formatting
is accomplished using the \name{sprintf} function, it can be switched off
by setting the format to NULL.

\subsection{Description of Functions}

In the following we will give a brief description of the operations and
public member functions of the class:
\begin{itemize}
\item{\name{Array2D$<$Base$>$ mat;} creates an empty matrix 
      \name{mat=NULL}.  }

\item{\name{Array2D$<$Base$>$ mat(n,m);} creates a n$\times$n 
       matrix \name{mat}
with indices starting at zero.}
\item{\name{Array2D$<$Base$>$ mat(n,m,off1,off2);} creates a 
   n$\times$m matrix
  \name{mat} with indices starting at off1 and off2 respectively.}
\item{\name{Array2D$<$Base$>$ mat1=mat2;} creates the vector \name{mat1} and 
   assigns \name{mat2} to it;}
\item{\name{mat.reset(n,m,off1,off2);} resets the
matrix to the arguments given, the offsets
\name{off1}, \name{off2}
are optional. The content of the previous vector is lost!}
\item{\name{mat.size1()} returns the first dimension (number of lines)
of the matrix.}
\item{\name{mat.size2()} returns the second dimension (number of rows)
of the matrix.}
\item{\name{mat.offset1()} and \name{mat.offset2()} return the offsets.}
\item{\name{mat.rebase(off1,off2)} sets the offsets to the paramters.}
\item{\name{mat.get\_base()} returns the value of the base-pointer to the class.
The use of this function is highly discouraged, since it defeats the purpose 
of creating a contained lcass in the first place. However, it is sometimes 
neccessary, when for instance a FORTRAN subroutine is to be called, for which
no specific header exists, or when sub-blocks of data are to be manipulated 
by functions such as \name{memcpy} etc. }
\item{\name{mat.set(value)} sets all elements of the matrix to \name{value}.  }
\item{\name{mat(n,m)} returns a reference to the element (n,m) of the matrix.}
\item{\name{mat.get(n,m)} is a constant function returning the value of
      the element (m,n) of \name{mat}.}
\item{ \name{mat1=mat2; } assigns \name{mat2} to \name{mat1}. The
target inherits both dimensions and offsets of the source, but not the name.
The target is resized if neccessary. }
\item{The \name{read} and \name{write} functions implement binary
input and output, reading/writing dimensions,offsets and elements of the
matrix. For a \name{String} argument e.g. \name{mat.write("String")},
the file is first opened, then
operation performed and finally the file is closed. If the argument is a
\name{fstream} or a \name{FILE* filepointer} only the io-operations are 
performed.  }
\item{ The output operator $<<$ provides a formatted output capability.
The name, the dimensions and offsets are printed, followed by the elements.
The user can influence the printing of the elements in two ways:
  \begin{itemize}
  \item{ The function \name{mat.set\_format(char *new\_format)}
  specifies a format string like " \%6.4f" (blanks must be included), as used
  in printf to format the individual elements.
  The function returns a pointer to the previous format string. The
  intial format is null. }
  \item{ The function \name{mat.set\_no\_elem(int n)} sets the number of elements
	 which are printed on a single line. This value defaults to 5. }
  \end{itemize} }
\end{itemize} 

\subsection{Class Header}
*/

template <class Base> class Array2D;
template <class Base> ostream& operator<<(ostream& os, const Array2D<Base>& vec);
  
template <class Base> class Array2D
{
  enum Array2D_Mode { unmarked, marked, mapped, mapped_marked };  
  Array2D_Mode  mode;
  friend class Array<Base>;
protected:

  int    dim1;
  int    dim2;
  int    off1;
  int    off2;
  
  Base*  base;
  Base** pbase; 

void   init();
public:  
       Array2D();
       Array2D(const int d1,const int d2,const int o1=0,const int o2=0);
       Array2D(const Array2D<Base>& mat);
virtual ~Array2D();
    
void   reset(const int d1,const int d2,const int o1=0,const int o2=0);
void   rebase(const int off1,const int off2);
int    size1() const { return dim1; }
int    size2() const { return dim2; }
int    offset1() const { return off1; }
int    offset2() const { return off2; }
void   set(const Base value);

void   read (const String& filename);
void   write(const String& filename) const;
void   read (fstream& file);
void   write(fstream& file) const;
void   read (FILE* fp);
void   write(FILE* fp) const;
void   error(const String& msg) const; 
void   error(const String& msg, const int n1,const int n2) const;
Array2D<Base>& operator=(const Array2D<Base>& rhs);
  
Base&   operator()(const int p1,const int p2) { return pbase[p2][p1]; }
const Base& get(const int p1,const int p2) const  { return pbase[p2][p1]; }
Base*   get_base() const { return base; }  

void   print(ostream& os,int elem_per_line = 10) const;

};

/*TEX
\section{Implemenation} 
\subsection{Constructors \& Reset} 
*/

template <class Base>
void Array2D<Base>::init()
{  
  dim1 = 0; dim2 =0;
  off1 = 0; off2 =0;
  base  = NULL;
  pbase = NULL;
  mode  = unmarked;
}

template <class Base> Array2D<Base>::Array2D() { init(); }

template <class Base>
Array2D<Base>::Array2D(const int d1,const int d2,const int o1,const int o2)
{  init(); reset(d1,d2,o1,o2);  }


template <class Base>
Array2D<Base>::Array2D(const Array2D<Base>& mat)
{
  init();
  reset(mat.dim1,mat.dim2,mat.off1,mat.off2);
  int i;
  for(i=0; i<dim1*dim2; i++)
    base[i]=mat.base[i];
}

template <class Base> Array2D<Base>::~Array2D()
{
  if (mode != unmarked)
    error("Mapped Matrix destroyed.");
  if (base)
    delete base;  
  if (pbase)
  {
    pbase+=off2;
    delete [] pbase;
  }
}

template <class Base>
void  Array2D<Base>::reset(const int d1,const int d2,const int o1,const int o2)
{
  if (mode != unmarked)
    error("Attempt to resize marked matrix. ");
  
  if (dim1*dim2 > 0) 
  {
    delete [] base;	
    pbase+=off2;
    delete [] pbase;
  }
  
  base  = NULL;	
  off2  = 0;
  pbase = NULL;
  					
  if (d1 < 0 || d2 < 0 ) 
     error("Size must be greater or equal to zero.");	 
  if (d1*d2 > 0)
  {
     base  = new Base[d1*d2];
     pbase = new Base* [d2];
     if (!base || !pbase) 
     {
       cerr << "Dimension: " << d1 << " " << d2 << endl;
       error("Memory Allocation Error.");		 
     }
  }
  dim1 = d1; dim2 = d2;
  rebase(o1,o2);
}

template <class Base>
void Array2D<Base>::rebase(const int o1,const int o2)
{
  if (dim2 > 0 && dim1 > 0)
  {
    pbase+=off2;
    off1=o1; off2=o2;
    for(int j=0; j<dim2;j++)
      pbase[j] = base + j*dim1 -off1;
    pbase-=off2;
  }
  else pbase = NULL;
}

template <class Base>
void Array2D<Base>::set(const Base value)
{
  int len = dim1*dim2;
  int i;
  for(i=0; i<len;i++)
    base[i] = value;
}

template <class Base>
Array2D<Base>& Array2D<Base>::operator=(const Array2D<Base>& rhs)
{
  if ((rhs.dim1 != dim1) || (rhs.dim2 != dim2) ||
      (rhs.off1 != off1) || (rhs.off2 != off2))
    reset(rhs.dim1,rhs.dim2,rhs.off1,rhs.off2);
  int i;
  for(i=0; i<dim1*dim2; i++)
    base[i]=rhs.base[i]; 
  return *this;
}


/*TEX
\subsection{Input/Output}
*/

template <class Base>
void Array2D<Base>::read(const String& file_name)
{
  char* tmp = file_name();
  fstream file(tmp,ios::in);  
  read(file);
  delete tmp;
}

template <class Base>
void  Array2D<Base>::read(fstream& file)
{
  if (!file)
    error("Input file not found");
  int l1;
  file.read((char*) &l1,sizeof(l1));
  if(!file || file.gcount() != sizeof(l1))
    error("Read Error for length.");
  int noff1;
  file.read((char*) &noff1,sizeof(off1));
  if (!file || file.gcount() != sizeof(off1)) 
    error("Read Error for offset.");

  int l2;
  file.read((char*) &l2,sizeof(l2));
  if(!file || file.gcount() != sizeof(l2))
    error("Read Error for length.");
  int noff2;
  file.read((char*) &noff2,sizeof(off2));
  if (!file || file.gcount() != sizeof(off2)) 
    error("Read Error for offset.");
  
  if (l1 != dim1 || l2 != dim2)
    reset(l1,l2,noff1,noff2);
  else
    if(off1 != noff1 || off2 != noff2)
    	rebase(noff1,noff2);
      
  if (dim1*dim2 > 0)
  {
    file.read((char*) base,sizeof(Base)*dim1*dim2);
    if (!file || file.gcount() != dim1*dim2*sizeof(Base)) 
      error("Read Error for base.");
  }
  
}

template <class Base>
void  Array2D<Base>::read(FILE* fp)
{
  if ( fp == NULL )
    error("Input filepointer not found");
  int l1;
  if( fread( &l1,sizeof(int),1,fp ) == 0 )
    error("Read Error for length.");
  int noff1;
  if( fread( &noff1,sizeof(int),1,fp ) == 0 )
    error("Read Error for offset.");

  int l2;
  if( fread( &l2,sizeof(int),1,fp ) == 0 )
    error("Read Error for length.");
  int noff2;
  if( fread( &noff2,sizeof(int),1,fp ) == 0 )
    error("Read Error for offset.");
  
  if (l1 != dim1 || l2 != dim2)
    reset(l1,l2,noff1,noff2);
  else
    if(off1 != noff1 || off2 != noff2)
    	rebase(noff1,noff2);
      
  if (dim1*dim2 > 0)
    if( fread( base,sizeof(Base),dim1*dim2,fp ) == 0 )
      error("Read Error for base.");
}

template <class Base>
void Array2D<Base>::write( const String& file_name ) const 
{
  char* tmp = file_name();
  fstream file(tmp,ios::out);
  write(file);
  delete tmp;  
}

template <class Base>
void Array2D<Base>::write(fstream& file) const 
{
  if (!file)
    error("Open error for output file.");
  file.write((char*) &dim1,sizeof(dim1));
  if(!file)
    error("Write Error for length.");
  file.write((char*) &off1,sizeof(off1));
  if(!file)
    error("Write Error for offset.");

  file.write((char*) &dim2,sizeof(dim2));
  if(!file)
    error("Write Error for length.");
  file.write((char*) &off2,sizeof(off2));
  if(!file)
    error("Write Error for offset.");
  if (dim1*dim2 > 0)
    file.write((char*) base,sizeof(Base)*dim1*dim2);
  if (!file)
    error("Write Error for base.");
}

template <class Base>
void Array2D<Base>::write(FILE* fp) const 
{
  if ( fp == NULL )
    error("No filepointer for output");
  if( fwrite( &dim1,sizeof(int),1,fp ) == 0 )
    error("Write Error for length.");
  if( fwrite( &off1,sizeof(int),1,fp ) == 0 )
    error("Write Error for offset.");

  if( fwrite( &dim2,sizeof(int),1,fp ) == 0 )
    error("Write Error for length.");
  if( fwrite( &off2,sizeof(int),1,fp ) == 0 )
    error("Write Error for offset.");

  if (dim1 * dim2 > 0)
    if( fwrite( base,sizeof(Base),dim1*dim2,fp ) == 0 )
    error("Write Error for base.");
}

template <class Base> 
ostream& operator<<(ostream& os, const Array2D<Base>& mat) 	   
{
  mat.print(os);
  return os;
}

#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
template <class Base> 
ostream_withassign& operator<<(ostream_withassign& os, const Array2D<Base>& mat) 	   
{
  mat.print(os);
  return os;
}
#endif

template <class Base> void Array2D<Base>::print(ostream& os,int elem_per_line) const
{							   
  if (dim1*dim2 == 0) 
  { 
    os << "Matrix is empty. " << endl;
    return;
  }
  
  os << "  *** MATRIX *** DIM1: " << off1 << "/" << dim1+off1-1
     << " *** DIM2: " << off2 << "/" << dim2+off2-1 << endl;
     
  int i1,i2;
  char buf[100];
  for(i1=0;i1<dim1;i1++)
  {
    for(i2=0;i2<dim2;i2++)
    {							   
      if (i2 % elem_per_line == 0) 
      { 
        int max = i2 + elem_per_line - 1; 				   
        if (max > dim2-1) max = dim2 - 1;  		   
        sprintf(buf,"\n%4i,%4i..%4i *",off1+i1,off2+i2,off2+max);
        os << buf;
      }			
      os << " " << get(i1+off1,i2+off2);
    }
  }
  os << endl;					   
}	
			
template <class Base>	
void Array2D<Base>::error(const String& msg)  const 
{  								 
  cerr << "Generic Matrix Error: "          
       << msg << endl;                                             
  abort();
}                                                                
			
template <class Base>					 
void Array2D<Base>::error(const String& msg,const int n1,const int n2) const
{  								 
  cerr << "Generic Matrix "            
       << n1 << " " << n2 << " --- "				 
       << msg << endl;                                             
  abort();
}                                                                

/*TEX
\section{Access Checking for 2D Arrays}
  
As for 1D arrays, an access-checked version of the
template is provided.

*/ 
  

template <class Base> class D_Array2D : public Array2D<Base>
{
public:
        D_Array2D() {}
        D_Array2D(const int d1,const int d2,const int o1=0,const int o2=0) :
        Array2D<Base>(d1,d2,o1,o2) {}
        D_Array2D(const Array2D<Base>& mat) : Array2D<Base>(mat) { }
        D_Array2D(const D_Array2D<Base>& mat)  : Array2D<Base>(mat) { }
virtual ~D_Array2D() {}

//compGcc4
using Array<Base>::off1;
using Array<Base>::dim1;
using Array<Base>::off2;
using Array<Base>::dim2;
using Array<Base>::pbase;
using Array<Base>::error;

D_Array2D<Base>& operator=(const D_Array2D<Base>& mat)
                   { Array2D<Base>::operator=(mat);
                     return *this; }
 
Base&   operator()(const int p1,const int p2) 
        { 
	  if(p1 < off1 || p1 >= dim1+off1 || p2 < off2 || p2 >= dim2+off2)
          error("Acces Violation",p1,p2);
	  return pbase[p2][p1];
	}
const Base&    get(const int p1,const int p2) const
        { 
	  if(p1 < off1 || p1 >= dim1+off1 || p2 < off2 || p2 >= dim2+off2)
          error("Acces Violation",p1,p2);
	  return pbase[p2][p1];
	}
};



#endif
