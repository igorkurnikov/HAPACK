#ifndef TIMER_DEF
#define TIMER_DEF
#include <time.h>
#include <assert.h>
#include <string.h>

#ifdef GNU
const long timer_scale = 1;
#else
const long timer_scale = 1000;
#endif

const long timer_fix   = 4294967;

class Timer
{
  int max;
  char         **name;
  long *msec;
  long *count;  
  long *stime;
  int  *code;
  int level;
  int maxdepth;
  
public:
  Timer(const int mx=100,const int mxdep = 20);
 ~Timer();

long get_time(int i) const;
void set_name(int i,char* newname);

void start(const int timer_code );
void stop(const int timer_code );  
void print();
};

#if !defined(XTIME_CPP)
extern Timer* timer;
#else
Timer* timer;
#endif

#endif
