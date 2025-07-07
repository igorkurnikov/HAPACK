/*TEX
\subsection{OrbInfo}

This class gathers the Orbital information according to symmetry
classes to compile the external index table. Within the integral
program we enumerate the basis functions according to the hierarchy of
InternalBasis (i.e. Center), total angular momentum, function index
within the angular momentum and momentum index, each distinguished by
a symmetry index. In all further processing the functions are
enumerated first by symmetry class and then by function count. The
\name{OrbInfo} class provides this mapping. In this enumeration the
indices first run over all functions of symmetry A, then B etc.
Unfortunately the map can only be computed AFTER all functions have
been enumerated in the generation of the \name{InternalBasis} sets.
   
*/
#include <numer.h>
#include <tsstack.h>
#include "typesip.h"
#include "symdesig.h"

class DataVault;


class OrbInfo
{
  int   no;
  ARRAY<SStack<int,0> > angular;
  ARRAY<SStack<int,0> > nsym;
  
  SStack<String,0> center;
  SStack<   int,0>    mom;
  SStack<   int,0>  index;
public:
  OrbInfo();
 ~OrbInfo() {}
  
int record(int fno,int l,SymDesignator& s)
{
  angular[l].add(fno);
  no++;
  return nsym[s].add(fno);
}

int  make_orbital_map(ARRAY<unsigned int>& map_vector);
int  no_orbs         () const { return no; } 
void print(fstream& os    ,int atom_flag);
void save (DataVault& data);
};



