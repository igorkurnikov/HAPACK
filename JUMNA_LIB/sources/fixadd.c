//int __imp___fpieee_flt( unsigned long exc_code, struct _EXCEPTION_POINTERS *exc_info, int handler(_FPIEEE_RECORD *) )
#include <stdio.h>
#include <time.h>
#include <string.h>


int _imp___fpieee_flt( unsigned long exc_code, void* p1, void* p2)
{
	return 1;
}

/*
void fdate_( char* f_date )
{ 
  time_t timer;
  char *date;
//  char *fdate;
// fdate = *f_date;
  timer = time((time_t *)NULL);
  date = ctime(&timer);
  strncpy(f_date,date,24);
} 

*/

float dtime_( float* tarray)
{ 
   tarray[0] = (float) (clock()) / CLK_TCK;
   tarray[1] = 0;
   return(tarray[0]+tarray[1]);
} 

//extern int __argc;
//extern char** __argv;
extern int* __p___argc();
extern char*** __p___argv();

/*
int iargc_ () {
    int __argc = *__p___argc();
	return ((int) (__argc - 1)); 
}

void getarg_ (int* num, char* arg, int len) 
{
	int i,lena;
	char** __argv = *__p___argv();
	strncpy ( arg, *(__argv + *num), len );
	lena = strlen(arg);
	for(i = lena; i < len; i++)
	{
        arg[i] = 32;
	}
}

*/