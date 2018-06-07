#include <numer.h>

main()
{
  Matrix m;
 
  m.reset("hallo",10,10);

  m.set(0.03);
  
  cout << m;
 
}
