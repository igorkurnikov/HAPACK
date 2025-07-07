#include "io_incl.h"
#include "orbinfo.h"
#include "vault.h"

OrbInfo::OrbInfo()
{
  angular.reset(max_momentum);
  nsym   .reset(8);
  int i;
  for(i = 0; i < max_momentum; i++)
    angular[i].resize(max_sym_orb);
  for(i = 0; i < 8; i++)
    nsym[i].resize(max_sym_orb);  
  no = 0;
  
}

int  OrbInfo::make_orbital_map(ARRAY<unsigned int>& map)
{
  int no_orb = 0;
  int s;
  for(s = 0; s < nsym.size(); s++)
    no_orb += nsym[s].size();

  map.reset(no_orb);
  map.set(no_orb+1);
  
  int count = 0;
  for(s = 0; s < nsym.size(); s++)
    for(int f = 0; f < nsym[s].size(); f++)
      map[nsym[s][f]] = count++;

  return no_orb;
}

void OrbInfo::print(fstream& os,int atom_flag)
{
  if (atom_flag)
  {
    os << "A" << endl;
    int maxl = 0;
    for(int l = 0; l < max_momentum; l++)
      if (angular[l].size() > 0) maxl = l;
    os << maxl << endl;

    int s;
    for(s = 0; s < maxl ; s++)
      os << setw(4) << angular[s].size() << " 0 0 0 0 " << endl;
    os << endl;
    
    for(s = 0; s < maxl ; s++)
    {
      for(int f = 0; f < angular[s].size() ; f++)
      {
	if (f % 10 == 0) os << endl;
	os << setw(4) << angular[s][f];
      }	
      os << endl;
    }
  }
  else
  {
    os << "M" << endl;
    int np = SymOp::no_symmetry();
    int ns = (int) pow(2.0,np);
    int ocount = 0;
    
    os << ns << endl;
    
    int s;
    for(s = 0; s < ns ; s++)
    {
      os << setw(4) << nsym[s].size() << " 0 0 0 0 " << endl;
      ocount +=  nsym[s].size();
    }
    os << endl;

    IVector map(ocount);
    ocount = 0;
    for(s = 0; s < ns ; s++)
      for(int f = 0; f < nsym[s].size() ; f++)
      	map[nsym[s][f]] = ocount++;

    os << "M" << endl;
    
    for(int i = 0; i < ocount; i++)
    {
      if (i % 10 == 0) os << endl;      
      os << setw(5) << map[i] << " ";
    }
  }
}


void OrbInfo::save(DataVault& vault)
{
  int np = SymOp::no_symmetry();
  int ns = (int) pow(2.0,np);
  IVector  table(ns);
    
  int ocount = 0;        
  int s;
  for(s = 0; s < ns ; s++)
  {
    table[s] = nsym[s].size();
    ocount +=  nsym[s].size();
  }
  vault.insert("ORBITALS.TOTAL",table);
}

