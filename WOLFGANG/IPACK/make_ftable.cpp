#include "io_incl.h"
#include <stdlib.h>
#include <math.h>
#include <numer.h>
#include <errorip.h>
#include "typesip.h"
#include "f.h"

// call this routine only once to create the data in ftable.h
main()
{
  compute_ftable("ftable.h");  
}
