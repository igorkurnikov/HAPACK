#include <errorip.h>
#include "io_incl.h"
#include "orbital.h"
#include "vault.h"

void Orbital::init(istream& is, DataVault& vault, const OInitCode icode)
{
  switch(icode)
  {    
    case oinit_qci:
    {
      
      //  Set up the orbital arrays
      no      .reset(maxtype,Symmetry::size());
      start_no.reset(maxtype,Symmetry::size());
      
      //  Default set up
      
      NUMARRAY<int> stable;
      vault.retrieve("ORBITALS.TOTAL",stable);
      no.set(0);
      Symmetry s;
      for(; s(); ++s)
	no(sd_type,s) = stable(s);
      
      // read the number of I, A, LO, HIGH, VIRT orbitals
      // compute the number of SD orbitals from what is left
      
      OrbitalType t;
      
      for(t = mintype; t < OrbitalType(maxtype); ++t)
      {
	if (t == OrbitalType(sd_type)) continue;
	for( s.reset(); s(); ++s )
	{
	  is >> no(t,s);
	  no(sd_type,s) -= no(t,s);
	}    
      }
      
      // sanity check 
      
      for( s.reset(); s(); ++s )
	if (no(sd_type,s) < 0)
	  error(" Illegal Number of SD orbitals ! ");
      
      
      //  Set up start_no and other counters
      
      int count  = 0;
      for( OrbitalType tt; tt.in_range(); ++tt)
	for(Symmetry ss; ss(); ++ss) 
	{
	  start_no(tt,ss) = count;
	  count += no(tt,ss);
	}
      
      // set up the external index table
      
      extidx.reset(number_of_orbitals());
      int i;
      
      int ocount = 0;
      
      NUMARRAY<int> int_orbs;
      vault.retrieve("ORBITALS.INTERNAL",int_orbs);
      
      // note: the number of internal orbitals in the integral file must 
      // equal the number of inactive plus active orbitals in each symmetry segment
      
      for(s.reset();s();++s)
	if (int_orbs(s) != number_of_orbitals(i_type,s) + number_of_orbitals(a_type,s))
	  error("#internal != #active + #inactive !!!");
      
      
      // fill the indices for the internal/cas_orbitals
      
      ARRAY<int> start(Symmetry::size());
      start[0] = 0;
      s.reset();
      Symmetry ss(s);
      ++ss;
      for(; ss();++ss,++s){
	start[ss] = start[s] + number_of_orbitals(i_type,s) + number_of_orbitals(a_type,s);
      }
      Symmetry s1;
      
      for(s1.reset(); s1(); ++s1)
	for( i = 0; i < no(i_type,s1); i++)
	  extidx[ocount++] = start[s1]++;
      
      for(s1.reset(); s1(); ++s1)
	for( i = 0; i < no(a_type,s1); i++)
	  extidx[ocount++] = start[s1]++;
      
      start[0] = number_of_orbitals(i_type) + number_of_orbitals(a_type);
      s.reset();
      ss = s;
      ++ss;
      for(; ss();++ss,++s)
	start[ss] = start[s] + number_of_orbitals(s) - 
	  number_of_orbitals(i_type,s) - number_of_orbitals(a_type,s);
      
      int x = start[0];
      for(t = lo_type; t.in_range(); ++t)
	for(s1.reset(); s1(); ++s1)
	  for( i = 0; i < no(t,s1); i++)
	    //extidx[ocount++] = x++;
	    extidx[ocount++] = start[s1]++;
      
      no_ext_orb = ocount;
      
      orb.reset( count );
      orb_no.reset( number_of_disk_orbitals());
      i = 0;
      
      for(OrbitalIter o(ALL_TYPE); o.in_range(); ++o ) 
      {
	orb   [i++] = o();
	orb_no[o().external_index()] = o().number();
      }

      sym_from_number.reset(number_of_orbitals());
      for (OrbitalIter oo(NOVIRT_TYPE);oo.in_range();++oo)
	sym_from_number[oo().number()] = oo().symmetry();
      break;

      if (number_of_orbitals(NOVIRT_TYPE) >= 256 && sizeof(Orbtype) == 1){
	cout << "ERROR: #orbitals >= 256. Change Orbtype and recompile!!!" << endl;
	abort();
      }
    }
  default:
    {
      //  Set up the orbital arrays
      no      .reset(maxtype,Symmetry::size());
      start_no.reset(maxtype,Symmetry::size());
      
      //  Default set up
      
      NUMARRAY<int> stable;
      vault.retrieve("ORBITALS.TOTAL",stable);
      no.set(0);
      Symmetry s;
      for(; s(); ++s)
	no(a_type,s) = stable(s);
      
      if( vault.check("ORBITALS.INTERNAL"))
      {
	NUMARRAY<int> int_info;
	vault.retrieve("ORBITALS.INTERNAL",int_info);
	if ( int_info.size() != Symmetry::size() )
	  error( "Fatal Error in Orbital Init: Wrong DIM of INT_INFO" );
	for( s.reset(); s(); ++s )
	{
	  no(sd_type,s) = no(a_type,s) - int_info(s);
	  no( a_type,s) = int_info(s);
	  if (no(sd_type,s) < 0 )
	    error("init(): number of specified internal orbitals exceeds number of available orbitals" );
    }
      }  
  
      // sanity check 
      
      for( s.reset(); s(); ++s )
	if (no(a_type,s) < 0)
	  error(" Illegal Number of SD orbitals ! ");
      
      
      //  Set up start_no and other counters
      
      int count  = 0;
      for( OrbitalType tt; tt.in_range(); ++tt)
	for(Symmetry ss; ss(); ++ss) 
	{
	  start_no(tt,ss) = count;
	  count += no(tt,ss);
	}
      
      // set up the external index table
      
      extidx.reset(number_of_orbitals());
      int i;
      
      int ocount = 0;
      
      // fill the indices for the internal/cas_orbitals
      
      ARRAY<int> start(Symmetry::size());
      start[0] = 0;
      s.reset();
      Symmetry ss(s);
      ++ss;
      for(; ss();++ss,++s){
	start[ss] = start[s] + number_of_orbitals(i_type,s) + number_of_orbitals(a_type,s);
      }
      Symmetry s1;
      
      for(s1.reset(); s1(); ++s1)
	for( i = 0; i < no(i_type,s1); i++)
	  extidx[ocount++] = start[s1]++;
      
      for(s1.reset(); s1(); ++s1)
	for( i = 0; i < no(a_type,s1); i++)
	  extidx[ocount++] = start[s1]++;
      
      start[0] = number_of_orbitals(i_type) + number_of_orbitals(a_type);
      s.reset();
      ss = s;
      ++ss;
      for(; ss();++ss,++s)
	start[ss] = start[s] + number_of_orbitals(s) - 
	  number_of_orbitals(i_type,s) - number_of_orbitals(a_type,s);
      
      int x = start[0];
      OrbitalType t;
      for(t = lo_type; t.in_range(); ++t)
	for(s1.reset(); s1(); ++s1)
	  for( i = 0; i < no(t,s1); i++)
	    extidx[ocount++] = start[s1]++;
      
      no_ext_orb = ocount;
      
      orb.reset( count );
      orb_no.reset( number_of_disk_orbitals());
      i = 0;
      
      for(OrbitalIter o(ALL_TYPE); o.in_range(); ++o ) 
      {
	orb   [i++] = o();
	orb_no[o().external_index()] = o().number();
      }

      sym_from_number.reset(number_of_orbitals());
      for (OrbitalIter oo(NOVIRT_TYPE);oo.in_range();++oo)
	sym_from_number[oo().number()] = oo().symmetry();
    }
  }
  
  if (number_of_orbitals(NOVIRT_TYPE) >= 256 && sizeof(Orbtype) == 1){
    cout << "ERROR: #orbitals >= 256. Change Orbtype and recompile!!!" << endl;
    abort();
  }
}








