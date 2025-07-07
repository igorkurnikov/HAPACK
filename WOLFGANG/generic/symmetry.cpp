#include "parallel.h"
#include "symmetry.h"
#include "vault.h"

#include "io_incl.h"

int Symmetry::v_maxsym;
int Symmetry::v_object;
ARRAY<int> Symmetry::deg;

void Symmetry::clear()
{
  deg.reset(0);
}
 
void Symmetry::init( istream & is ) 
{
  char obj;
  is >> obj;

  switch ( obj ) {
    case 'A' : v_object = atom;
               is >> v_maxsym;
               break;
    case 'M' : v_object = mol;
               is >> v_maxsym;
               break;
    default  : cerr << "main() : Illegal type for obj";
               pm -> abort();
               break;
  }
  //(ostream&) cout << v_maxsym <<endl;
  //cout.flush();
  v_maxsym--;

  deg.reset(v_maxsym+1);
  for ( int i=0; i<=v_maxsym; i++ )
    deg[i] = ( v_object == atom ? 2*i+1 : 1 );

}
 
void Symmetry::init(DataVault& vault) 
{
  String type = vault.retrieve("MOLTYPE",vault_severe);
  char obj = type[0];
  NUMARRAY<int> st;
  vault.retrieve("ORBITALS.TOTAL",st);
  v_maxsym = st.size();
  
  switch ( obj ) 
  {
  case 'A' :
    v_object = atom;
    break;
  case 'M' : 
    v_object = mol;
    break;
  default  : 
    cerr << " Symmetry::init -- Illegal MOLTYPE";
    pm -> abort();
    break;
  }
  //(ostream&) cout << setw(8) << v_maxsym <<endl;
  //cout.flush();
  v_maxsym--;

  deg.reset(v_maxsym+1);
  for ( int i=0; i<=v_maxsym; i++ )
    deg[i] = ( v_object == atom ? 2*i+1 : 1);

}



static char moltext[] = {'A','B','C','D','E','F','G','H'};
static char atomtext[] = {'S','P','D','F','G','H','I','J'};

istream & operator>>(istream & is,Symmetry & s)
{
  char c;
  is>>c;
  char* text_ptr = ( Symmetry::v_object == mol ? moltext : atomtext );
  int i = 0;
  while ( i < 8 && c != text_ptr[i] ) i++;
  if (i>=s.size()) cerr<<"Symmetry out of range: " << c << endl;
  s.v_sym=i;
  return is;
} 

ostream & operator<<(ostream & os,const Symmetry & s)
{
  os << ( Symmetry::v_object == mol ? 
                       moltext[s.v_sym] : atomtext[s.v_sym] );
  return os;
}

  
