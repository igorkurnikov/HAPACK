/*TEX
\section{Memory Management} 
*/

#ifndef MEMORY_H
#define MEMORY_H

#include "tsstack.h"
typedef struct 
{
  char* ptr;
  int   sz;
} MemEntry;

ostream& operator<<(ostream& os, const MemEntry& entry);

class MemoryManager
{
      SStack<MemEntry,100> entry;
      long  first_free;
      long  bufsize;
      char* buffer;
public:  
      MemoryManager();
     ~MemoryManager();
char* alloc(const long size);
void  free (const char* p);  
long  find (const char* p);  
void  garbage_collection();  
void  print(ostream& os);
      
};

#endif

   
