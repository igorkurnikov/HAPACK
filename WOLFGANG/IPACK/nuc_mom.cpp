#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <errorip.h>
#include <math.h>
#include "basis.h"
#include "vault.h"

void nuclei_moment(Molecule& molec,MomentData mu_data,DataVault &vault)
{
  cout << " +++ NucMom:: Position: " << mu_data.loc << endl;
  SStack<Location,max_cntr> eff_center;
  SStack<  double,max_cntr> eff_charge;

  for(int i1 = 0; i1 < molec.size(); i1++)
  {
//    cout << " +++ NucMom:: Processing center: " << molec.center(i1) << endl;
    for(SymOp oper1(molec.center(i1).sym_equiv(X),
		    molec.center(i1).sym_equiv(Y),
		    molec.center(i1).sym_equiv(Z));oper1.in_range();++oper1)
    {
      Location l1(molec.center(i1));
      l1.apply(oper1);
//      cout << " op: " <<  oper1 << l1 << endl;
      if (fabs((double) molec.element(i1).charge()) > 0)
      {
	eff_center.add(l1);
	eff_charge.add(molec.element(i1).charge());      
      }
    }
  }

  Momentum help;
  Vec Nuc_Moments(help.min(mu_data.max_mu+1)-1);
  for( int m=1; m <= mu_data.max_mu; m++)
  for(Momentum mu(m);mu();++mu)
  {
    Nuc_Moments[mu.index()-1] = 0;
    for(int i = 0; i < eff_center.size(); i++)
    {
      double moment = 
	pow((eff_center[i](X) - mu_data.loc(X)),mu(X))
	*pow((eff_center[i](Y) - mu_data.loc(Y)),mu(Y))
	*pow((eff_center[i](Z) - mu_data.loc(Z)),mu(Z));
  
      Nuc_Moments[mu.index()-1] += eff_charge[i]*moment;
    }
      
  }
  vault.insert("NUC_MOMENTS",Nuc_Moments);
} 





