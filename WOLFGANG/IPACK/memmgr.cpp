#include "io_incl.h"
#include <errorip.h>
#include "memmgr.h"

ostream& operator<<(ostream& os, const MemEntry& entry)
{
  os << "     " <<  setw(10);
  hex(os);
  os << (int) entry.ptr << setw(10) << entry.sz;

  return os;
}

typedef int (*COMP)(const void*,const void*);

int compare(MemEntry* a, MemEntry* b)
{
  if (a -> ptr >   b -> ptr) return 1;
  if (a -> ptr ==  b -> ptr) return 0;
  return -1;
}

MemoryManager::MemoryManager()
{
  bufsize = 1000000L;
  buffer  = new char[bufsize];
  assert(buffer != 0);
  first_free = 0;
}

char* MemoryManager::alloc(const long size)
{
  if (first_free + size > bufsize)
    error("Memory Overflow. ");

  MemEntry m;
  m.ptr = buffer + first_free;
  m.sz  = size;
  entry.add(m);
  
  first_free += size;  

  return m.ptr;
}

long MemoryManager::find(const char* p)
{
  long mn = 0;
  long mx = entry.size();

  while(mx-mn>1)
  {
    int nxttry = (mn+mx) / 2;
    if (entry[nxttry].ptr < p)
      mn = nxttry;
    else 
      mx = nxttry;
  }
  if (entry[mx].ptr == p) return mx;
  if (entry[mn].ptr == p) return mn;
  
  return -1;
}

void MemoryManager::free (const char* p)
{
  int pos = find(p);
  if (pos == -1)
    error("Illegal Pointer in Free. ");
  entry[pos].ptr = 0;
}

void MemoryManager::garbage_collection()
{
 
  // trow away all empty pointers

  first_free = 0;
  int count = entry.size();
  while(count-- > 0)
  {
    if (entry[count].ptr == 0)
    {
      if (count == entry.size())
	entry.pop();
      else
	entry[count] = entry.pop();
    }
    else
      if (first_free < entry[count].ptr + entry[count].sz - buffer)
	first_free = entry[count].ptr + entry[count].sz - buffer;
  }

  qsort(entry.get_base(),entry.size(),sizeof(MemEntry), (COMP) compare);

}  


void MemoryManager::print(ostream& os)
{
  os << "MemoryManager -- First Free: " << first_free << endl;
  os << "          POINTER     SIZE " << endl;
  os << entry;
}

MemoryManager::~MemoryManager()
{
  delete buffer;
}
