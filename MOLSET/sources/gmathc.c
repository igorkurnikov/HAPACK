/*
 *  some basic C functions with gaussian interface
 */

#include <string.h>
#include <time.h>

#ifdef _WIN32
/* Windows-specific implementation */
#include <stdlib.h>  /* provides __argc / __argv macros via UCRT */

int iargc_ () {
	return ((int) (__argc - 1));
}

void getarg_ (int* num, char* arg, int len)
{
	int i,lena;
	char** argv_local = __argv;
	strncpy ( arg, *(argv_local + *num), len );
	lena = strlen(arg);
	for(i = lena; i < len; i++)
	{
		arg[i] = 32;
	}
}

#else
/* Linux/Unix - these functions are provided by gfortran runtime */
/* We provide weak symbols that won't conflict with libgfortran */

/* On Linux, iargc_ and getarg_ are provided by libgfortran */
/* If they're not available, provide stub implementations */
__attribute__((weak)) int iargc_() {
	return 0;
}

__attribute__((weak)) void getarg_(int* num, char* arg, int len) {
	int i;
	for(i = 0; i < len; i++) {
		arg[i] = ' ';
	}
}

#endif

#ifndef CLK_TCK
#define CLK_TCK CLOCKS_PER_SEC
#endif

float etime_ (tarray) float *tarray; {
  tarray[0] = (float) (((double) clock()) / ((double) CLK_TCK));
  tarray[1] = 0;
  return(tarray[0]+tarray[1]);
}

void fdate_( f_date ) char* f_date;
{
  time_t timer;
  char *date, *fdate;
  fdate = f_date;
  timer = time((time_t *)NULL);
  date = ctime(&timer);
  strncpy(fdate,date,24);
}
