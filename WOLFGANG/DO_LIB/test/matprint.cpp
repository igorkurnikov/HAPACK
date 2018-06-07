#include <numer.h>

main()
{
  const int n = 10;
  const int m = 12;
  Matrix a(n,m);
  
  for(int i=0; i< n; i++)
  for(int j=0; j< m; j++)
    a(i,j) = 0.1233234*i + j;
  
  a.print(cout,10," %10.4f ");
  

}
