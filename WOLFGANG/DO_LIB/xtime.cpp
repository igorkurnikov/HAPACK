/*
%  Dortmund C++ Class and Template Library 
%  Copyright (C) 1994 Wolfgang Wenzel and others

%  The code in this file was derived from routines published in 
%  Numerical Recipes, it its the responsibility of the user to insure that
%  he is authorized to use them. Please read the appropriate section entitled
%  LICENSE in the documentation. 

%  This library is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%  Library General Public License for more details.

%  You should have received a copy of the GNU Library General Public
%  License along with this library; if not, write to the Free
%  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%
% $Id: xtime.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
*/

#define XTIME_CPP

#include "io_incl.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "xtime.h"


Timer::Timer(const int mx,const int mxdep)
{
  max      = mx;
  maxdepth = mxdep;  

  msec  = new  long[max];
  assert(msec !=0);
  count = new  long[max];
  assert(count !=0);
  name  = new char*[max];
  assert(name != 0);
  stime = new  long[maxdepth];
  assert(stime !=0);
  code = new   int[maxdepth];
  assert(code !=0);
  
  for(int i=0; i<max;i++)
  {
    msec[i]  = 0;
    count[i] = 0;
//    name[i]  = strdup("FUNCTION");
	name[i] = (char*) malloc(strlen("FUNCTION")+1);
	strcpy(name[i],"FUNCTION");
  }
  level = 0;  
}

Timer::~Timer() 
{
  delete msec;
  delete count;
  delete stime;  
  for(int i=0; i < max; i++)
    free(name[i]);
  delete name;
}

long 
Timer::get_time(int i) const 
{ 
	return (long) (1000.0*msec[i]/CLOCKS_PER_SEC*timer_scale); 
} 


void 
Timer::set_name(int i,char* newname)
{
  assert(i >=0 && i < max);
  free(name[i]);
  name[i] = strdup(newname);
}

void 
Timer::start(const int timer_code )  
{ 
  code [++level] = timer_code;
  stime[  level] = clock() / timer_scale;
  if (level >= maxdepth) cerr << "Timer Overflow." << endl;
}

void 
Timer::stop(const int timer_code )   
{
  if (timer_code != code[level])
  {
    cout << "Error in Timer: Trying to Stop: " << timer_code 
      << " expected: " << code[level] << endl;
  }
  
  long etime = (long) clock()/ timer_scale;
  if (etime < stime[level])
  {
    cout << endl << "TIMER OVERFLOW level: " << level
         << " START: " << stime[level] << " END: " << etime << " NEW: ";
    for(int i =0 ; i < level; i++)
    {
      msec[i]  += etime + timer_fix - stime[i];
      stime[i]  = etime;
    }
    etime += timer_fix;
#if !defined(__DECCXX)
        cout <<  etime << " DIFF: " << ( etime - stime[level]) <<  endl;
#endif
  }
  msec [timer_code] += etime - stime[level];  
  count[timer_code]++;
  if (--level<0) 
    cerr << "Timer Underflow" << endl; 
}


void Timer::print()
{
  char buffer[100];
  
  cout << " +++ TIMING INFORMATION: " << endl;
  cout << endl << "               NAME   NO      COUNT     TIME       TIME/CALL " << endl;
  for(int i=0;i<max;i++)
    if (count[i] > 0) 
    {
      double ms = get_time(i);
      sprintf(buffer,"    %15s %4i %10i %10.3f    %10.3f \n",
	      name[i],i,count[i],ms/1E3,1.0*ms/count[i]);
      cout <<  buffer;
    }
  cout << endl;
  
}





