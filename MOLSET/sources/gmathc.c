/*
 *  some basic C functions with gaussian interface 
 */

#include <string.h>
#include <time.h>

#define DllImport   __declspec( dllimport )

/*DllImport int __argc;
DllImport char** __argv;

int iargc_ () {
	return ((int) (__argc - 1)); 
}

void getarg_ (int* num, char* arg, int len) 
{
	int i,lena;
	strncpy ( arg, *(__argv + *num), len );
	lena = strlen(arg);
	for(i = lena; i < len; i++)
	{
        arg[i] = 32;
	}
}*/
extern int* __p___argc();
extern char*** __p___argv();


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


float etime_ (tarray) float *tarray; {
  tarray[0] = (float) (((double) clock()) / ((double) CLK_TCK));
  tarray[1] = 0;
  return(tarray[0]+tarray[1]);
}

//float etime_ (tarray) float *tarray; {
//#ifdef __MAC__
//  unsigned long t;
//  ReadDateTime(&t);
//  tarray[0] = (float) (t - 2940247000);
//  tarray[1] = 0.0;
//#endif
//#else
//  clock_t times(), wall;
//  struct tms tbuf;
//  wall = times(&tbuf);
//  tarray[0] = (float) (((double) tbuf.tms_utime) / ((double) CLK_TCK));
//  tarray[1] = (float) (((double) tbuf.tms_stime) / ((double) CLK_TCK));
//#endif
//  return(tarray[0]+tarray[1]);
//  }

void fdate_( f_date ) char* f_date;
{ 
  time_t timer;
  char *date, *fdate;
/*   fdate = CH_F2C(f_date); */
  fdate = f_date;
  timer = time((time_t *)NULL);
  date = ctime(&timer);
  strncpy(fdate,date,24);
}  
