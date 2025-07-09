#include <stdlib.h>
#if defined(_MSC_VER)
#define NO_RLIMIT_STACK_CTRL
#include <windows.h>
#include <winbase.h>
#include <winnt.h>
#endif

#if !defined(_MSC_VER)
#include <sys/time.h>
#endif

#ifndef NO_RLIMIT_STACK_CTRL
#include <sys/resource.h>
#endif

void    get_bytesize_(void * start, void * end, int * bytes);
void    get_wall_time_(int * sec, int * usec);
void    unlimit_stack_(int * new_limit);

void
get_bytesize_(void * start, void * end, int * bytes)
{
  *bytes = ((char *)end - (char *)start);
}

void get_wall_time_(int * sec, int * usec)
{
#if defined _MSC_VER
	FILETIME file_time;
	UINT64 time_long;

	// may be I need to use QueryPerformanceCounter()
	GetSystemTimeAsFileTime(&file_time);

    time_long = (((UINT64)file_time.dwHighDateTime) << 32) | file_time.dwLowDateTime;

	time_long = time_long/10;
	
	*usec = (int) (time_long % 1000000);

	time_long = time_long/1000000; 
	time_long &= 0x0000000000FFFFFF;
	
	*sec = (int) time_long;
	
//	*sec = time(NULL);
//	*usec = 0;
#else
  struct timeval        wall_time;

  gettimeofday(&wall_time, NULL);

  *sec = wall_time.tv_sec;
  *usec = wall_time.tv_usec;
#endif
  return;
}

/*
 * This routine returns -1 to indicated that the stack has been reset to the
 * maximum.  For machines that define NO_RLIMIT_STACK_CTRL because they don't
 * have the rlimit interface, we are basically lying...  If the routine returns
 * a positive value, that is the actual limit up to 1,000,000,000.  Any limit
 * over 1,000,000,000 is considered "unlimited" for purposes of this code (this
 * actually helps keep 32/64 bit issues from unnecessarily setting off alarms
 * on systems like AIX.
 */

void unlimit_stack_(int * new_limit)
{
#ifdef NO_RLIMIT_STACK_CTRL
  *new_limit = -1;
#else
  struct rlimit         rlim;

  getrlimit(RLIMIT_STACK, &rlim);

  if (rlim.rlim_cur == RLIM_INFINITY)
  {
    /* the stack is already unlimited */
    *new_limit = -1;
  }
  else
  {
    if (rlim.rlim_max == RLIM_INFINITY)
    {
      /* the stack can be set to unlimited by anyone */
      rlim.rlim_cur = rlim.rlim_max;
      setrlimit(RLIMIT_STACK, &rlim);
      *new_limit = -1;
    }
    else
    {
      /* there is a hard limit in effect; bump up to it */

      if (rlim.rlim_max >= 1000000000) /* we consider 1 billion "unlimited" */
        *new_limit = -1;
      else
        *new_limit = (int)rlim.rlim_max;

      if (rlim.rlim_cur < rlim.rlim_max)
      {
        rlim.rlim_cur = rlim.rlim_max;
        setrlimit(RLIMIT_STACK, &rlim);
      }
    }
  }
#endif

  return;
}
