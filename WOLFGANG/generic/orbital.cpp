/*TEX
%
% IPACK - BSR
% Quantum Chemistry Project: Basis Set Reduction 
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is proprietary software. Unlicensed use and distribution
% of this program or parts thereof are illegal. For details on the 
% license, see the LICENSE section in the file "main.c" distributed 
% with this package or write to: 
%      Wolfgang Wenzel, Theoretical Physics I, 
%      University of Dortmund,
%      D-44221 Dortmund, Germany 
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the LICENSE
% for more details.
%
% $Id: orbital.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
\section{Implementation of Orbitals}  
\subsection{IO}

We provide overloaded operators for both input and output. The  
convention is to name an orbital by type, symmetry and index in the
order given.  
*/
#include "parallel.h"
#include "io_incl.h"
#include <errorip.h>
#include "vtype.h"
#include "orbital.h"
#include "vault.h"


#ifdef GNU
void gnufix() // never called but necessary for compilation with GNU
{
  Array2D<double> a;
  Array2D<double> b(a);
  Array<Orbital> oa;  
}
#endif

ostream & operator<<(ostream & os,const Orbital& orb)
{
  os << orb.type() << orb.symmetry() << orb.index();
  return os;
}

istream & operator>>(istream & is,Orbital& orb)
{
  Symmetry sym;
  int i;
  OrbitalType t;
  is>>t;
  is>>sym;
  is>>i;
  Orbital o(t,sym,i);
  orb=o;
  return is;
}

/*TEX
\subsection{Order of Orbitals}
Orbitals are ordered first by their type, then by symmetry and finally
by their index.
*/  

int         operator< (const Orbital& o1,const Orbital& o2)
{
  if (o1.type()     < o2.type())     return 1;
  if (o1.type()     > o2.type())     return 0;
  if (o1.symmetry() < o2.symmetry()) return 1;
  if (o1.symmetry() > o2.symmetry()) return 0;
  if (o1.index()    < o2.index())    return 1;
  return 0;
}
int         operator> (const Orbital& o1,const Orbital& o2)
{
  if (o1.type()     > o2.type())     return 1;
  if (o1.type()     < o2.type())     return 0;
  if (o1.symmetry() > o2.symmetry()) return 1;
  if (o1.symmetry() < o2.symmetry()) return 0;
  if (o1.index()    > o2.index())    return 1;
  return 0;
}
int         operator<=(const Orbital& o1,const Orbital& o2)
{
  return !(o1 > o2);
}
int         operator>=(const Orbital& o1,const Orbital& o2)
{
  return !(o1 < o2);
}
int         operator==(const Orbital& o1,const Orbital& o2)
{
  return (o1.type() == o2.type() && o1.symmetry() == o2.symmetry()
	  && o1.index() == o2.index());
}
int         operator!=(const Orbital& o1,const Orbital& o2)
{
  return !(o1 == o2);
}
/*TEX
\subsection{Initialization of Static Variables for the Orbitals} 
  
This function read the number of orbitals of each type and symmetry
*/

NUMARRAY2D<int> Orbital::no;
NUMARRAY  <int> Orbital::skip;
NUMARRAY2D<int> Orbital::start_no;
ARRAY<Orbital>  Orbital::orb;
ARRAY<int>      Orbital::orb_no;

ARRAY<int>      Orbital::extidx;
ARRAY<int>      Orbital::inactive_id;
int             Orbital::inactive_no;
int             Orbital::no_ext_orb;
NUMARRAY<int>   Orbital::map;

void Orbital::clear()
{
  no.reset(0,0);
  skip.reset(0);
  start_no.reset(0,0);
  orb.reset(0);
  orb_no.reset(0);
  extidx.reset(0);
  map.reset(0);
}

void Orbital::init(istream& in,DataVault& vault,const OInitCode code)
{
  //  Set up the orbital arrays

  no      .reset(maxtype,Symmetry::size());
  start_no.reset(maxtype,Symmetry::size());
  skip    .reset(Symmetry::size());
  skip    .set(0);
  
  //  Default set up

  NUMARRAY<int> stable;
  vault.retrieve("ORBITALS.TOTAL",stable);
  no.set(0);
  Symmetry s;
  for(; s(); ++s)
    no(internal_type,s) = stable(s);

  if( vault.check("ORBITALS.INTERNAL"))
  {
    NUMARRAY<int> int_info;
    vault.retrieve("ORBITALS.INTERNAL",int_info);
    if ( int_info.size() != Symmetry::size() )
      error( "Fatal Error in Orbital Init: Wrong DIM of INT_INFO" );
    for( s.reset(); s(); ++s )
    {
      no(external,s) = no(internal_type,s) - int_info(s);
      no(internal_type,s)   = int_info(s);
      if ( no(external,s) < 0 )
	error("init(): number of specified internal orbitals exceeds number of available orbitals" );
    }
  }
  
  NUMARRAY<int> bsr_orbs(Symmetry::size());
  NUMARRAY<int> inactive(Symmetry::size());
  bsr_orbs.set(0);
  inactive.set(0);
  
  if (code == oinit_mrci)
  {
    // check whether orbitals are to be excluded
    char c;
    in >> c;
    if (c != 'I')
    {
      error("MRCI:: expected list of inactive orbitals");
    }
    
    for(s.reset(); s(); ++s)
    {
      in  >> inactive[s];
      if (inactive(s) > Orbital::number_of_orbitals(internal_type,s))
      {
	cout << " +++ Illegal Orbital Number in Segment I: " << s << " " << inactive(s) << endl;
	error("Illegal Orbital Number in Segment I");
      }
    }
    vault.insert("ORBITALS.INACTIVE",inactive);

    in >> c;    
    bsr_orbs.set(0);
    switch(c)
    {
    case 'x':
    case 'X':
      cout << " +++ Exclude orbitals : ";
      for(s.reset(); s(); ++s)
      {
	in  >> bsr_orbs[s];
	cout << setw(4) << bsr_orbs[s];
      }
      cout << endl;
      vault.insert("ORBITALS.BSR",bsr_orbs,vault_ignore);
      break;
    case 'v':
    case 'V':
      {
      NUMARRAY<int> valence(Symmetry::size());
      for(s.reset(); s(); ++s)
      {
	in  >> valence[s];
	bsr_orbs[s] = no(external,s) - valence(s);
	if (bsr_orbs(s) < 0)
	  error("Illegal Orbital Number in Segment.");
      }
      vault.insert("ORBITALS.BSR",bsr_orbs,vault_ignore);
      break;
      }
    case 'n':
    case 'N':
      vault.insert("ORBITALS.BSR",bsr_orbs,vault_ignore);
      break;
    default:
      error("MRCI: Expected token X,V or N to specifiy externals");
    }
  }
  else
  {
    if (vault.check("ORBITALS.BSR"))
      vault.retrieve("ORBITALS.BSR",bsr_orbs);
    if (vault.check("ORBITALS.INACTIVE"))
      vault.retrieve("ORBITALS.INACTIVE",inactive);
  }
  
  // obtain exclude info and adjust the orbital map
  // this is mandatory for the above 2 cases
  
  for(Symmetry ss; ss(); ++ss) 
  {
    no(external,ss) -= bsr_orbs(ss);
    // if ( code == oinit_mrci ) 
    skip[ss] = bsr_orbs(ss);
    //    else
    no(bsr_type,ss) = bsr_orbs(ss);
  }

  //  Print table
  if (pm -> master())
    print_table( cout, 0 );
  
  //  Set up extidx
  extidx_set_up();

  //  Set up start_no and other counters

  int count  = 0;
  for( OrbitalType tt; tt.in_range(); ++tt)
    for(Symmetry ss; ss(); ++ss) 
    {
      start_no(tt,ss) = count;
      count += no(tt,ss);
    }
  
  orb        .reset( count );
  inactive_id.reset( count );
  orb_no     .reset( number_of_disk_orbitals());
  
  int i = 0;
  for(OrbitalIter o(ALL_TYPE); o.in_range(); ++o ) 
  {
    orb   [i++] = o();
    orb_no[o().external_index()] = o().number();
  }

  inactive_id.set(-1);
  count = 0;
  for(s.reset(); s(); ++s)
    for(int i = 0; i < inactive[s];i++)
    {
      Orbital o(internal_type,s,i);
      inactive_id[o.number()] = count++;
    }
  inactive_no = count;
  
}

/*TEX
\subsection{\name{from\_map\_to\_extidx()}} 
This function is called by the init().
*/

void Orbital::extidx_set_up()
{
  extidx.reset(number_of_orbitals());
  int i;

  int xcount = 0;
  int ocount = 0;
  
  Symmetry s;
  for(s.reset(); s(); ++s)
    for( i =0; i < no(internal_type,s); i++)
      extidx[xcount++] = ocount++;

  NUMARRAY<int> bsr_start(Symmetry::size());
  for(s.reset(); s(); ++s)
  {
    for( i = 0; i < no(external,s); i++)
      extidx[xcount++] = ocount++;
    bsr_start[s] = ocount;
    ocount += no(bsr_type,s); //  + skip(s);
  }

  for(s.reset(); s(); ++s)
  {
    ocount = bsr_start(s);
    for( i = 0; i < no(bsr_type,s); i++)
      extidx[xcount++] = ocount++;
  }
  no_ext_orb = ocount;
}

/*TEX
\subsection{\name{set\_up\_map()}}
We define a map which (a) assigns negative numbers to frozen and
virtual orbitals, and (b) reorders the remaining orbitals by 
type.
*/ 

int Orbital::set_up_map( const int& core_done, const int& Back )
{
  int count = 0;
  /*
  int i, j, k;
  int dim, dim_off;
  Orbital orbit;
  OrbitalType t, tt;
  OrbitalType t_freeze;
  ++t_freeze; ++t_freeze; ++t_freeze;    //  Increment t_freeze to freeze
  Symmetry s;
  map.set(-2);                           //  Default is virt

  //  Mark the freeze orbitals

  int offset = 0;
  for( s.reset(); s(); ++s) {
    dim = orbit.number_of_orbitals(freeze,s);
    for ( j=offset; j<offset+dim; j++ ) 
      map[j] = -1;
    offset += orbit.number_of_orbitals(s);
  }

  //  Set up the rest of the map

  for( t=internal_type; t<t_freeze; ++t ) {

    offset = 0;
    for( s.reset(); s(); ++s) {

      dim_off = orbit.number_of_orbitals(freeze,s);
      for ( tt=internal_type; tt<t; ++tt )
        dim_off += orbit.number_of_orbitals(tt,s);
      dim = orbit.number_of_orbitals(t,s);
      for ( j=offset+dim_off; j<offset+dim_off+dim; j++ ) 
        map[j] = count++;
      offset += orbit.number_of_orbitals(s);

    }

  }

  //  If the core has already been frozen remove the negative numbers

  if ( core_done ) {
    NUMARRAY<int> temp(count);
    int i = 0;
    for ( j=0; j<map.size(); j++ )
      if ( map[j] >= 0 ) temp[i++] = map[j];
    map = temp;
  }

  //  If back-transform reverse the map

  if ( Back ) {
    NUMARRAY<int> temp(count);
    for ( j=0; j<map.size(); j++ )
      temp[map[j]] = j;
    map = temp;
  }
  */
  return count;

}

/*TEX
\subsection{Print Orbital Table}
*/
void Orbital::print_table( ostream& os, int verbose ) {

  Symmetry s;

  os << endl << " +++ Orbital Table: \n\n";

  os << "  Type       ";
  for( s.reset(); s(); ++s )
    os  << s << "       ";
  os <<endl<<endl;

  for( OrbitalType tt; tt.in_range(); ++tt) {
    os << "    " << tt << " ";
    for( s.reset(); s(); ++s )
      os << setw(8) << no(tt,s);
    os << endl;
  }
  os << endl;

  if ( verbose ) 
  {
    os << " +++ Orbital Index Table: " << endl;
    for ( int i=0; i<map.size(); i++ ) 
      os << "    Orbital: " << setw(4) << i << " maps to : " 
	 << setw(4) << map[i] << endl;
  }
}

/*TEX
\subsection{Incrementing Orbital}
*/
Orbital&  Orbital::operator++()
{
  if (++v_i < no(v_type,v_s)) return *this;
  v_i = 0;
  ++v_s;
  while(v_type < OrbitalType(maxtype))
  {
    // increase v_s 
    while(v_s < Symmetry::size() && no(v_type,v_s) == 0)
      ++v_s;
    if (v_s < Symmetry::size()) return *this;
    v_s = Symmetry();
    ++v_type;
  }
  return *this; 
}

/*TEX
\subsection{Number of Orbitals}
*/

int Orbital::number_of_orbitals(const OrbitalType& type)
{
  int sum  = 0;
  for(Symmetry s; s(); ++s)
    sum += no(type,s);
  return sum;
}

int Orbital::number_of_orbitals( const Symmetry& s )
{
  int sum  = 0;
  for(OrbitalType t; t.in_range(); ++t)
    sum += no(t,s);
  return sum;
}

int Orbital::number_of_orbitals()
{
  int sum  = 0;
  for(OrbitalType t; t.in_range(); ++t)
    sum += number_of_orbitals(t);
  return sum;
}

/*TEX
\section{Implementation of OrbitalType}
\subsection{IO}
We define a lower case letter for each orbitals type.
*/ 
static char *orbitaltypetext[]={"i","x","b"};

ostream & operator<<(ostream & os,const OrbitalType& type)
{
  os<<orbitaltypetext[type];
  return os;
}

istream & operator>>(istream & is,OrbitalType& type)
{
  char c;
  is>>c;
  switch(c) {
    case 'i': type = internal_type; break;
    case 'x': type = external; break;
    case 'b': type = bsr_type; break;
    default:
      cerr << "Attempt to read Illegal Orbitaltype: Input was" << c << endl;
      abort();
  }
  return is;
}

/*TEX
\section{Implentation of OrbitalIterators} 

\subsection{Constructor}  

The constructor initializes the internal \name{Orbital} variable to
the first orbital in the range, if such an \name{Orbital} exists. In
this case the \name{ok} flag is set to \name{true}, otherwise to
\name{false}. 

*/

OrbitalIter::OrbitalIter(OrbitalTypeRange target_type) : type(target_type)
{
  ok = false;
  switch(type)
  {
  case ALL_TYPE:
  case INT_TYPE:
  case IX_TYPE:
    {  
      for(Symmetry s; s(); ++s)
	if (Orbital::number_of_orbitals(internal_type,s)  > 0)
	{
	  orb = Orbital(internal_type,s,0);
	  ok  = true;
	  break;
	}
      if (ok) break;
      if (type == INT_TYPE) break;
    }
    
  case EXT_TYPE:
  {  
    for(Symmetry s; s(); ++s)
      if (Orbital::number_of_orbitals(external,s) > 0)
      {
	orb = Orbital(external,s,0);
	ok  = true;
	break;
      }
    if (ok) break;
    if (type ==  IX_TYPE) break;
    if (type == EXT_TYPE) break;
   }
  case BSR_TYPE:
    {
      for(Symmetry s; s(); ++s)
	if (Orbital::number_of_orbitals(bsr_type,s) > 0)
	{
	  orb = Orbital(bsr_type,s,0);
	  ok  = true;
	  break;
	}
      if (ok) break; 
    }
  }
} 

OrbitalIter& OrbitalIter::operator++() 
{
  ++orb;

  switch(type)
  {
  case ALL_TYPE:
    if (orb.type()  < OrbitalType(maxtype)) return *this;
  case BSR_TYPE:
    if (orb.type() <= OrbitalType(bsr_type)) return *this;
  case EXT_TYPE:
  case IX_TYPE:
    if (orb.type() <= OrbitalType(external)) return *this;
  case INT_TYPE:
    if (orb.type() == OrbitalType(internal_type)) return *this;
  }
  ok = 0;
  return *this;
}


OrbitalSymIter::OrbitalSymIter(OrbitalTypeRange target_type,Symmetry s) 
     : type(target_type),tsym(s)
{
  ok = false;
  switch(type)
  {
  case ALL_TYPE:
  case INT_TYPE:
  case IX_TYPE:
  {  
    if (Orbital::number_of_orbitals(internal_type,s) > 0)
    {
      orb = Orbital(internal_type,s,0);
      ok  = true;
      break;
    }
    if (ok) break;
    if (type == INT_TYPE) break;
   }
    
  case EXT_TYPE:
  {  
      if (Orbital::number_of_orbitals(external,s) > 0)
      {
	orb = Orbital(external,s,0);
	ok  = true;
	break;
      }
    if (ok) break;
    if (type ==  IX_TYPE) break;
    if (type == EXT_TYPE) break;
   }
  case BSR_TYPE:
  {
    
      if (Orbital::number_of_orbitals(bsr_type,s) > 0)
      {
	orb = Orbital(bsr_type,s,0);
	ok  = true;
	break;
      }
    if (ok) break; 
  }
    
  }
} 

OrbitalSymIter& OrbitalSymIter::operator++() 
{
  OrbitalType tt = orb.type();
  ++orb;

  // this is a little complicated: the question is: if all orbitals of the desired symmetry
  // of a given orbitaltype are exhausted, which is the next allowed orbital type.
  // in searching we must strat with the type of the last allowed orbital, since the ++
  // operation may or may not increment the type
							     
  if (orb.symmetry() != tsym)
  {
    while(++tt < OrbitalType(maxtype) && Orbital::number_of_orbitals(tt,tsym) == 0)
    if (tt == OrbitalType(maxtype)) 
    {
      ok = 0;
      return * this;
    }
    orb = Orbital(tt,tsym,0);
  }
  
  switch(type)
  {
  case ALL_TYPE:
    if (orb.type()  < OrbitalType(maxtype)) return *this;
  case BSR_TYPE:
    if (orb.type() <= OrbitalType(bsr_type)) return *this;
  case EXT_TYPE:
  case IX_TYPE:
    if (orb.type() <= OrbitalType(external)) return *this;
  case INT_TYPE:
    if (orb.type() == OrbitalType(internal_type)) return *this;
  }
  
  ok = 0;
  return *this;
}


