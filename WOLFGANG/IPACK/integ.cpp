#include "io_incl.h"
#include <assert.h>
#include "typesip.h"
#include "basis.h"
#include "f.h"

main()
{  
  String ftablename("FTABLE");
  read_ftable(ftablename);
  
  PureBasisFunction pure_f11(0);   // make a pure uncontracted basis function
  PureBasisFunction pure_f12(1);
  PureBasisFunction pure_f13(2);
  pure_f11.add(0.3,1.0);           // of momentum 0 and exponent 0.3
  pure_f12.add(2.2,1.0);           // of momentum 1 and exponent 2.2
  pure_f13.add(3.3,1.0);           // of momentum 2 and exponent 3.3

  PureBasisSet set1;               // make a basis set of these functions
  set1.add(&pure_f11);
  set1.add(&pure_f12);
  set1.add(&pure_f13);

  PureBasisFunction pure_f21(0);   
  PureBasisFunction pure_f22(1);   
  PureBasisFunction pure_f23(2);
  pure_f21.add(4.00,1.0);
  pure_f22.add(0.25,1.0);
  pure_f23.add(5.00,1.0); 
  
  PureBasisSet set2;
  set2.add(&pure_f21);
  set2.add(&pure_f22);
  set2.add(&pure_f23);

  PureBasisFunction pure_f31(0);   
  PureBasisFunction pure_f32(1);   
  PureBasisFunction pure_f33(2);
  pure_f31.add(8.00,1.0);
  pure_f32.add(2.5,1.0);
  pure_f33.add(12.00,1.0); 
  
  PureBasisSet set3;
  set3.add(&pure_f31);
  set3.add(&pure_f32);
  set3.add(&pure_f33);

  PureBasisFunction pure_f41(0);   
  PureBasisFunction pure_f42(1);   
  PureBasisFunction pure_f43(2);
  pure_f41.add(0.20,1.0);
  pure_f42.add(25.2,1.0);
  pure_f43.add(1.00,1.0); 
  
  PureBasisSet set4;
  set4.add(&pure_f41);
  set4.add(&pure_f42);
  set4.add(&pure_f43);

  /* now combine the pure basis sets to make one real basis set */

  Location center1(0, 0, 0);  
  Location center2(1,-1 ,1);
  Location center3(0.3,0.1,0.7);
  Location center4(-2,-3 ,5);
  
  BasisSet basis;
  
  basis.add(&set1,center1);
  basis.add(&set2,center2);  
  basis.add(&set3,center3);
  basis.add(&set4,center4);

  // basis complete

  basis.prepare();
  
  // compute the overlaps

  IntegMat ovlp;
  IntegMat kin;  
  IntegMat nuc;  
  
  basis.overlap(ovlp);
  cout << ovlp;

  basis.kinetic(kin);  
  cout << kin;

  basis.nuclear(nuc);  
  cout << nuc;

  twoel(basis,basis,true);
  
}

