#include <unistd.h>
#include "tarray.h"

main()
{
  const int l1 = 6;  
  const int l2 = 3;  
  int i,j;
  Array2D<double> a(l1,l2);
  for (i=0;i < l1;i++) 
  for (j=0;j < l2;j++) 
    a(i,j) = (3+i)*(2+j);

  cout<<" ++++++ Testing Array2D: \n";
  cout<<a;


}
