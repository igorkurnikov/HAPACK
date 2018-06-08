#include "io_incl.h"
#include "storage.h"

         int  AllocBlock::init_flag = 0; 	
         int  AllocBlock::first[256];
unsigned char AllocBlock::one_mask[8];

int     Allocator::all_id   = 0;

void AllocBlock::init()
{
  if (init_flag) return;
  init_flag = 1;
  int i;
  for(i = 0; i < 8; i++)
  { 
    one_mask[i] = 1 << i;
  }
  int lo  = 1;
  int hi  = 2;
  
  for(i=0; i < 256; i++)
  {
    int bit;
    for(bit = 0; bit < 8; bit++)
    {
      int bb = i >> bit;
      int cc = bb & 1;
      if (cc == 0)
	break;
    }
    first[i] = bit;
  }
}

Allocator::Allocator(int s1,int s2,int s3,int s4,int nb) 
{
  AllocBlock::init();      
  sz[0] = s1;
  sz[1] = s2;
  sz[2] = s3;
  sz[3] = s4;
  sz[4] = nb;
  tot_sz = s1*s2*s3*s4*nb;
  idd    = ++all_id;

  // cout << " Creating Allocator (cs): " << idd << " with: " << tot_sz << endl;
  
#ifndef MEM_SAVE
  AllocBlock* ab = new AllocBlock(tot_sz);
  block.add(ab);
#endif
  last_free = 0;
}

Allocator::Allocator(const Allocator& src) 
{
  AllocBlock::init();      
  sz[0] = src.sz[0];
  sz[1] = src.sz[1];
  sz[2] = src.sz[2];
  sz[3] = src.sz[3];
  sz[4] = src.sz[4];
  tot_sz = src.tot_sz;
  idd    = ++all_id;

  // cout << " Creating Allocator (cp): " << idd << " with: " << tot_sz << endl;
#ifndef MEM_SAVE
  AllocBlock* ab = new AllocBlock(tot_sz);
  block.add(ab);
#endif
  last_free = 0;
}

Allocator::~Allocator()
{
  // cout << " Removing Allocator: " << idd << endl;
  idd = -1;
 
  int bz = block.size();
  int i;
  for( i = 0; i < bz; i++)
    if (block[i] -> free() != block[i] -> max_sz())
      cout << "Block not empty: " << i << " " <<  block[i] -> free() << " " 
	   << block[i] -> max_sz() << endl;

  for( i = 0; i < bz; i++)
    delete block.pop();
}



int Allocator::check()
{
  int errflag = 0;
  for(int i = 0; i < block.size(); i++)
    if (block[i] -> error_check(i))
      errflag = 1;
  return errflag;
}













