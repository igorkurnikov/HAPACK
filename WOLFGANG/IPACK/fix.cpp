/*
   This fix function is neccessary for compiling ipack on a linux system
   with the gnu C++ compiler version 2.7.0. There seems to be a bug in
   resolving templates
*/

#include "parallel.h"
#include <tarray.h>
#include "basis.h"
#include "coef_array.h"
#include "symop.h" 
#include "symdesig.h" 
#include "tpair.h"
#include "function.h"
#include "integ_file.h"

class AllocBlock;
class HRR_Set;
class IntegArray;
class RecursionStorage;

void gnu_fix()
{
  IntegCoreFile bf("AA",1);
  
  ARRAY<Vector> av(1,2);

  ARRAY<RadialFunctionArray> rfa1;
  ARRAY<RadialFunctionArray> rfa2(rfa1);

  ARRAY<IntegArray*> irf1;
  ARRAY<IntegArray*> irf2(irf1);

  ARRAY<RecursionStorage*> xirf1;
  ARRAY<RecursionStorage*> xirf2(xirf1);

  ARRAY<RadialFunction> rf1;
  ARRAY<RadialFunction> rf2(rf1);
  rf2 = rf1;

  ARRAY<ARRAY<Location> > l1;
  ARRAY<ARRAY<Location> > l2(l1);

  ARRAY<ARRAY<Momentum> > ml1;
  ARRAY<ARRAY<Momentum> > ml2(ml1);

  
  ARRAY<SymOp>  aaps1;
  ARRAY<SymOp>  aaps2(aaps1);

  ARRAY<Pair<SymOp> > aps1;
  ARRAY<Pair<SymOp> > aps2(aps1);

  ARRAY<Pair<SymDesignator> > apsd1;
  ARRAY<Pair<SymDesignator> > apsd2(apsd1);

  ARRAY<ARRAY<String> > sa1;
  ARRAY<ARRAY<String> > sa2(sa1);

  ARRAY2D<double> x;
  ARRAY2D<double> y(x);

  ARRAY2D<IVector> ix;
  ARRAY2D<IVector> iy(ix);

  ARRAY<IVector> iv1;
  ARRAY<IVector> iiv2(iv1);

  ARRAY<ARRAY<IVector> > ix1;
  ARRAY<ARRAY<IVector> > iy1(ix1);

  ARRAY2D<HRR_Set*> hs1;
  ARRAY2D<HRR_Set*> hs2(hs1);

  ARRAY<InternalBasis *> ib1;
  ARRAY<InternalBasis *> ib2(ib1);
  ib2 = ib1;
  ARRAY<PureBasisFunction *> if1;
  ARRAY<PureBasisFunction *> if2(if1);
  if2 = if1;

  ARRAY<CoefArray* > ic1;
  ARRAY<CoefArray* > ic2(ic1);
  ic2 = ic1;

  ARRAY2D<Coefficient_Set* > ics12d;

  ARRAY<Coefficient_Set* > ics1;
  ARRAY<Coefficient_Set* > ics2(ics1);
  ics2 = ics1;

  ARRAY<AllocBlock* > ab1;
  ARRAY<AllocBlock* > ab2(ab1);
  ab2 = ab1;

  Array<IntegData> idat1;
  Array<IntegData> idat2;
  idat1 = idat2;
}






