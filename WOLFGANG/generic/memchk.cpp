#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errorip.h>
#include "memchk.h"


void  mallocmap();
/*
void heap_check(char*  mess)
{
  
  char  buf[100];
  sprintf(buf,"MEMCHK START XX%i %s",count++,mess);
  cout << buf << endl;
  cout.flush();
  msg("MEMCHK END ");
}

*/

void operator delete(void*  p)
{
  if(memory_manager) memory_manager -> remove(p);
  if(p) free((char* ) p);
}

void operator delete[](void*  p)
{
  if(memory_manager) memory_manager -> remove(p);
  if(p) free((char* ) p);
}

#if defined(MEMCHK_IPACK)

extern void* operator new(size_t  size)
{
  void*  p;
  p = malloc(size);
  if(memory_manager && memory_manager -> is_trace()) 
    memory_manager -> add(p,size);
  return p;
}

#endif

int  Memory_Manager::trace_flag = 0;
Memory_Manager*  memory_manager = NULL;
const int  STRLEN = 2;

Memory_Manager::Memory_Manager(int  nsz)
{
  tempflag   = trace_flag;
  trace_flag = 0;
  max        = nsz;
  if(max <= 0) error("Memory_Manager:: Illegal Argument.");
  n = 0;

  pointer = new void* [max];
  if(!pointer) error("Memory_Manager::ENOMEM.");
  size = new size_t [max];
  if(!size) error("Memory_Manager::ENOMEM.");

  no_used = new int [max];  
  if(!no_used) error("Memory_Manager::ENOMEM.");
  
  active = new int [max];
  if(!active) error("Memory_Manager::ENOMEM.");
  
  // message = new char [STRLEN];
  // if(!message) error("Memory_Manager::ENOMEM.");
  
  // pmsg = new char* [max];
  // if(!pmsg) error("Memory_Manager::ENOMEM.");
  
  int  i;
  /*
  for(i = 0;i < max;i++)
  {
    pmsg[i] = new char [STRLEN];
    if(!pmsg[i]) error("Memory_Manager::ENOMEM.");
  }
  */
  // strcpy(message,"A");
  hi_water      =  0;  
  tot_mem       =  0;
  check         = -1;
  allocno       =  0;
  tot_mem       =  0;
  hi_water      =  0;
  current_water =  0;
  trace_flag    =  tempflag;
}

Memory_Manager::~Memory_Manager()
{
  tempflag   = trace_flag;
  trace_flag = 0;
  delete  pointer;
  delete  size;
  delete  no_used;
  delete  active;
  
  int  i;
  // for(i = 0;i < max;i++)
  // delete  pmsg[i];
  // delete  pmsg;
  // delete  message;
  trace_flag = tempflag;
}

void Memory_Manager::intercept(int  alloc_no)
{  
  tempflag   = trace_flag; 
  trace_flag = 0;
  check      = alloc_no;
  cout << "    Memory_Manager::Intercept = " << alloc_no << endl;
  trace_flag = tempflag;
}

void Memory_Manager::add(void*  p,size_t  sz)
{
  tempflag = trace_flag;
  trace_flag = 0;
  if(!p) error("Memory_Manager::Out of Memory Error.");
  allocno++;
  current_water += sz;
  if(current_water > hi_water) hi_water = current_water;
  tot_mem += sz;
  
  int  i;  
  for(i = 0;i < n;i++)
    if(!active[i]) 
    {    
      active[i] = 1;
      size[i] = sz;
      no_used[i] = allocno;
      pointer[i] = p;
      // strcpy(pmsg[i],message);
      if(check == allocno) error("Memory_Manager -- Address Intercepted.");
      break;
  }
  
  if(i >= n) 
  {
    if(n < max) 
    {
      active [n] = 1;
      pointer[n] = p;
      no_used[n] = allocno;
      // strcpy(pmsg[n],message);
      if(check == allocno) error("Memory_Manager -- Address Intercepted.");
      size[n] = sz;      
      n++;
    }
    else 
    {
      cerr << "Memory_Manager::Overflow.";
      abort();
    }
  }
  trace_flag = tempflag;
}

void Memory_Manager::remove(void*  p)
{  
  tempflag = trace_flag; 
  trace_flag = 0;
  int  i;
  for(i = 0;i < n;i++)
    if(p == pointer[i]) break;
  
  if(i < n) 
  {
    active[i] = 0;    
    current_water -= size[i];
  }
  trace_flag = tempflag;
}


void Memory_Manager::print(ostream&  os,char*  msg,int  inactive)
{  
  tempflag = trace_flag; 
  trace_flag = 0;
  os << "    +++MEMORY-MAP: ";
  if(msg) cout << msg;
  cout << "\n\n";
  char  buf[100];
  cout << "    Address     Offset   Length  Alloc.No  Allocated in:\n";
  
  int  i;
  for(i = 0;i < n;i++)
    if(active[i]) 
    {
      sprintf(buf," %10x %10x     %4i    %6i \n",
	      pointer[i],(int ) ((char* ) pointer[i] - (char* ) pointer[0]),
	      size[i],no_used[i]);
      cout << buf;
    }
 
  if(inactive) 
  {   
    cout << "    Address     Offset   Length  Used  Allocated in:\n";
    for(i = 0;i < n;i++)
      if(!active[i]) 
      {
	sprintf(buf," %10x %10x     %4i    %2i\n",pointer[i],
		(int) ((char* ) pointer[i] - (char* ) pointer[0]),
		size[i],no_used[i]);
	cout << buf;
      }
  }
  os << "    +++END-MEMORY_MAP\n";
  trace_flag = tempflag;
}

void Memory_Manager::print_memory_usage(ostream&  os)
{  
  os << endl; 
  os << " +++ Memory Currently used: " << current_water << endl;
  os << " +++ Maximum Memory   used: " << hi_water << endl;
  os << " +++ Total Allocations    : " << tot_mem << endl;
  os << endl;
}

void Memory_Manager::mark(char*  msg)
{  
  tempflag = trace_flag; 
  trace_flag = 0;
  if(strlen(msg) >= STRLEN) return;
  // strcpy(message,msg);
  trace_flag = tempflag;
}
