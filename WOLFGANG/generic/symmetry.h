
#ifndef SYMMETRY_DEF
#define SYMMETRY_DEF

#include "io_incl.h"
#include <tnumarray.h>
#include "vtype.h"

/*TEX 

\section{Symmetry}

\name{obj\_tpye} enumerates whether the object being considered is an
                 atom or a molecule.

The maximum symmetry \name{v\_maxsym} and \name{v\_object} are static
varibles which are set once with the function \name{init()} -- note
there are two \name{init()} functions, so \name{v\_object} does not
have to be set.

Constructor has no arguments and sets symmetry to zero.  With the
reset function the same action can be performed.   The
operators <,>,==,!=,>= and <= are defined. the operator++ is defined
as post and prefix operator.  The function zero tests if the symmetry
is zero. The operator \^ calculates a new symmetry via a xor
operation. The operator () returns a logical value indicating whether
s in the allowed symmetry range or not. The cast operator for
transformation into integer is defined. Maxsym returns the maximum
symmetry */

enum obj_tpye {atom,mol};
class DataVault;

class Symmetry
{
 public:
  friend istream & operator>>(istream & is,Symmetry & s);
  friend ostream & operator<<(ostream & os,const Symmetry & s);
  Symmetry()   { v_sym=0; }
  Symmetry(const Symmetry& sym){v_sym = sym.v_sym;}
  Symmetry& operator=(const Symmetry& sym)
    {
      v_sym=sym.v_sym; 
      return *this;
    }
  void reset() { v_sym=0; } 
  static void init(int sym)  {v_maxsym=sym;} 
  static void init( istream & is );
  static void init( DataVault& vault );
  static void clear();
  int  operator< (const Symmetry & sym) const {return (v_sym<sym.v_sym);}
  int  operator> (const Symmetry & sym) const {return (v_sym>sym.v_sym);}
  int  operator>=(const Symmetry & sym) const {return (v_sym>=sym.v_sym);}
  int  operator<=(const Symmetry & sym) const {return (v_sym<=sym.v_sym);}
  int  operator==(const Symmetry & sym) const {return (v_sym==sym.v_sym);} 
  int  operator!=(const Symmetry & sym)const {return (v_sym!=sym.v_sym);}
  
  Symmetry operator ++()      { v_sym++;return *this;}
  Symmetry operator ++(int)   { Symmetry help=(*this); v_sym++; return help;} 
  int      A() const          { return (v_sym==0); }
  
  Symmetry operator ^(const Symmetry & sym) const
    { Symmetry help; help.v_sym=v_sym^sym.v_sym; return help; }
  int operator () () const { return (v_sym<=v_maxsym); }
  operator int() const { return v_sym; }
  static int size(){ return (v_maxsym+1); }  
  static int object(){ return (v_object); }
  int degeneracy() { return deg[v_sym]; }
protected: 
  unsigned char  v_sym;
  static int v_maxsym; 
  static int v_object; 
  static ARRAY<int> deg;

  friend void init( istream& is );
};

#endif








