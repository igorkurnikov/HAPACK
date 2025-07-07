#include "parallel.h"
#include "io_incl.h"
#include <assert.h>
#include <errorip.h>
#include <math.h>
#include "basis.h"
#include "vault.h"
#include "typesip.h"

void nuclei_e_field(Molecule& molec,Location e_field_loc,DataVault &vault)
{
  cout << " +++ NucEField:: Position: " << e_field_loc << endl;
  double nuc_e_field;
  SStack<Location,max_cntr> eff_center;
  SStack<  double,max_cntr> eff_charge;

  for(int i1 = 0; i1 < molec.size(); i1++)
  {
//    cout << " +++ NucEField:: Processing center: " << molec.center(i1) << endl;
    for(SymOp oper1(molec.center(i1).sym_equiv(X),
		    molec.center(i1).sym_equiv(Y),
		    molec.center(i1).sym_equiv(Z));oper1.in_range();++oper1)
    {
      Location l1(molec.center(i1));
      l1.apply(oper1);
      // cout << " op: " <<  oper1 << l1 << endl;
      if (fabs((double) molec.element(i1).charge()) > 0)
      {
	eff_center.add(l1);
	eff_charge.add(molec.element(i1).charge());      
      }
    }
  }

  int index = 0; // used for labeling the output

  for(Direction dir = X; dir <= Z; operator++(dir)) 
  {
    nuc_e_field = 0;
    for(int i = 0; i < eff_center.size(); i++)
    {
     double root = sqrt(
      (eff_center[i](X) - e_field_loc(X))*(eff_center[i](X) - e_field_loc(X))
     +(eff_center[i](Y) - e_field_loc(Y))*(eff_center[i](Y) - e_field_loc(Y))
     +(eff_center[i](Z) - e_field_loc(Z))*(eff_center[i](Z) - e_field_loc(Z)));
     if (root > delta)
     {
       nuc_e_field += eff_charge[i]*(eff_center[i](dir) - e_field_loc(dir)) 
                    / pow(root,3);
     }
     else
     {
       cout << " +++ Nuclear E_Field Error: " << eff_center[i] << 
	  e_field_loc<< endl;
       error("Divergency in Nuclear E_Field. ");
     }
    }
    cout << " Processing Nuclei E_Field in Direction " << dir << endl;
    String name="NUC_E_FIELD";
    ++index;
    cout << name+index << " " << nuc_e_field << endl;
    vault.insert(name+index,nuc_e_field);
  }

  for(Direction dir1 = X; dir1 <= Z ; operator++(dir1))
  for(Direction dir2 = dir1; dir2 <= Z ; operator++(dir2))
  {
    nuc_e_field = 0;
    for(int i = 0; i < eff_center.size(); i++)
    {
     double root = sqrt(
      (eff_center[i](X) - e_field_loc(X))*(eff_center[i](X) - e_field_loc(X))
     +(eff_center[i](Y) - e_field_loc(Y))*(eff_center[i](Y) - e_field_loc(Y))
     +(eff_center[i](Z) - e_field_loc(Z))*(eff_center[i](Z) - e_field_loc(Z)));
     if (root > delta)
     {
      nuc_e_field += 3*eff_charge[i]*(eff_center[i](dir1) - e_field_loc(dir2)) 
	*(eff_center[i](dir2) - e_field_loc(dir1)) / pow(root,5);
      if(dir1 == dir2)
        nuc_e_field -= eff_charge[i] / pow(root,3);
     }
     else
     {
       cout << " +++ Nuclear E_Field_Grad Error: " << eff_center[i] << 
	e_field_loc<< endl;
       error("Divergency in Nuclear E_Field_Grad. ");
     }
    }
    cout << " Processing Nuclei E_Field_Grad (" << dir1 << "," << dir2
         << ")" << endl;
    String name="NUC_E_Field_";
    ++index;
    cout << name+index << " " << nuc_e_field << endl;
    vault.insert(name+index,nuc_e_field);
  }
}






