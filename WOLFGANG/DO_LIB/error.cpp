
#include <errorip.h>

//mms fix for gnu on the SUNs

void error(char* msg)
{  cerr << "Error: " << msg << endl; abort(); }

