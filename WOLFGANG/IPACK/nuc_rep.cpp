#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <errorip.h>
#include <math.h>
#include "basis.h"


double nuclear_repulsion(Molecule& molec,int& core_size)
{
  double nuc_e = 0;
  core_size    = 0;
  SStack<Location,max_cntr> eff_center;
  SStack<  double,max_cntr> eff_charge;
    
  for(int i1 = 0; i1 < molec.size(); i1++)
  {
    if (pm -> master())
      cout << " +++ NucRep:: Processing center: " << molec.center(i1) << endl;

    for(SymOp oper1(molec.center(i1).sym_equiv(X),
		    molec.center(i1).sym_equiv(Y),
		    molec.center(i1).sym_equiv(Z));oper1.in_range();++oper1)
    {
      Location l1(molec.center(i1));
      l1.apply(oper1);
      if (fabs((double) molec.element(i1).charge()) > 0)
      {
	eff_center.add(l1);
	core_size += molec.element(i1).core();
	eff_charge.add(molec.element(i1).charge());      
      }
    }
  }
  
  for(int i = 0; i < eff_center.size(); i++)
  for(int j = 0; j < i; j++)
  {
    double dist = 
    (eff_center[i](X) - eff_center[j](X))*(eff_center[i](X) - eff_center[j](X))
   +(eff_center[i](Y) - eff_center[j](Y))*(eff_center[i](Y) - eff_center[j](Y))
   +(eff_center[i](Z) - eff_center[j](Z))*(eff_center[i](Z) - eff_center[j](Z));
    if (fabs(dist) > delta)
    {
      nuc_e += eff_charge[i]*eff_charge[j] / sqrt(dist);
    }
    else
    {
      cout << " +++ Nuclear Repulsion Error: " << eff_center[i] << 
	eff_center[j] << endl;
      error("Divergencye in Nuclear Repulsion. ");
    }
  }
  return nuc_e;
} 
