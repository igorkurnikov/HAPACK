#include <numer.h>
#include "memchk.h"

main()
{
  memory_manager = new Memory_Manager(100000);
  memory_manager -> on();
 
  if (argc > 1)
  {
    int allocno;
    allocno = atoi(argv[1]);
    cout << "intercepting: " << allocno << endl;
    memory_manager -> intercept(allocno);
  }

  Array<int> *a = new Array<int>(3);
  
  memory_manager -> print(cout);
}
