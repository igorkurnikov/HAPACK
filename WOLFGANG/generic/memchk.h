
void  heap_check(char * msg);

class Memory_Manager { 
private: 
static int  trace_flag;
int  tempflag;
int  max;
int  n;
int  allocno;
void**  pointer;
size_t*  size;
int*  no_used;
int*  active;
  //char**  pmsg;
  //char*  message;
int   check;
long  current_water;
long  hi_water;
long  tot_mem;




public: 

Memory_Manager(int  sz = 1000);

~Memory_Manager();

void add(void*  p,size_t  s);

void remove(void*  p);

int  is_trace() {  return trace_flag; }

void on()
{
  //  if TRUE allocations will be recorded,
  //  deletes are always recorded.
  //  records the value of tempflag in the 
  //  routines.
  //  maximum number of pointers to trace.
  //  current number
  //  counter for allocations.
  //  ARRAYS   
  //  pointers, 
  //  size of last allocation to this pointer
  //  how often this adress whas used in allocs.
  //  currently allocated ?
  //  active message at last allocation
  //  currently active message
  //  allocation number 
  //  amount of memory currently allocated
  //  largest amount of memory allocated at one time
  //  total memory allocated
  
  trace_flag = 1;
}


void off()
{
  //  dynamically enable tracing
  
  trace_flag = 0;
}


void print(ostream&  os,char*  msg = NULL,int  inactive = 0);

void print_memory_usage(ostream&  os);

void mark(char*  msg);

void intercept(int  allocno);

//friend void* operator new(size_t  size);

//friend void ::operator delete(void*  p);
};

extern Memory_Manager*  memory_manager;








