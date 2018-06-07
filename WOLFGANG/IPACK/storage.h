/*TEX
%
% IPACK - QCHEM
% Quantum Chemistry Project: Integral Package
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is proprietary software. Unlicensed use and distribution
% of this program or parts thereof are illegal. For details on the 
% license, see the LICENSE section in the file "ipack.c" distributed 
% with this package or write to: 
%      Wolfgang Wenzel, Theoretical Physics I, 
%      University of Dortmund,
%      D-44221 Dortmund, Germany 
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the LICENSE
% for more details.
%
% $Id: storage.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\chapter{Storage\label{chapter-storage}}
\subsection{Concepts}  

This is the lowest level class for the handling of both
\name{IntegArrays} and \name{CoefArrays}. The class is based on the
\name{Vector} class of the DOLIB template library. The class can either
implement a four-dimensional, or, as a subset, a two dimensional block
of arrays. The number of blocks is given by the fifth an uppermost
index to the array. 

%The data in class \name{Storage} may additionally
%have symmetries. If the data is symmetric under the interchange of the
%first and the second index, the flag \name{sym12} is true.  If the
%data is symmetric under the interchange of the third and the fourth
%index, the flag \name{sym34} is true.  If the data is symmetric under
%the interchange of the first and the second index-pair, the flag
%\name{sympair} flag is true.

The \name{add\_to},\name{fold} and \name{doublet} functions are used to
evaluate the individual terms which appear in an recursion
relation. The \name{normalize} function normalizes the primitive
gaussian orbitals in an array. Since contractions cannot be performed
on symmetrized arrays, the \name{to\_matrix} and \name{from\_matrix}
routines are needed to pack/unpack the symmetric arrays to/from
rectangular arrays.

Finally, the \name{output} routines are used to output the integrals.

\subsection{Class Header} 

The class maintains certain  resource-statistics in a set of static variables:
\begin{tabular}{ll}
\name{c\_number}    & number of \name{Storage} classes presently allocated. \\
\name{c\_size}      & size   of \name{Storage} classes presently allocated. \\
\name{c\_notot}     & number of \name{Storage} classes allocated so far. \\
\name{c\_totsize}   & size   of \name{Storage} classes allocated so far. \\
\name{c\_sizehigh}  & maximal number of all \name{Storage} classes allocated \\
\name{c\_nohigh}    & maximal size   of all \name{Storage} classes allocated \\
\end{tabular}

The sizes are in units of \name{sizeof(double)}. Since almost the
entire dynamic storage requirements of this program occur in derivates
of \name{Storage} these counts present an accurate balance on the
memory requirements of the program.

The class provided to working arrays, used in various subroutines,
\name{work1} and \name{work2}. Each subroutine using these arrays may
not assume anything about initial their length or contents.

There is a rather obvious set of routines, \name{size} etc, which
returns the present allocation parameters and symmetries. There are two
functions \name{to\_matrix} and \name{from\_matrix} which convert
to and from the symmetric to unsymmetrized storage modes. This is
important for contractions. 
*/

#ifndef STORAGE_H
#define STORAGE_H
#include "typesip.h"
#include <numer.h>
#include <errorip.h>
#include <tsstack.h>
#include "qc_utilities.h"

class InternalBasis;
class IntegFile;
class Quadruple_Set;
class SymDesignator;

typedef struct 
{
  IntegType  *base;
  unsigned      id;  // id = (alloc_no << 16 | block_id << 8 | block_no)
} StorageUnit;

/*
An AllocBloc contains 256 units of storage of size sz;
*/
class AllocBlock
{
         int        free_count;
         int        mode;
         IntegType* base;
unsigned char       mask[32];
	 int        size;
	 
	 int        mxbyte;
static   int        init_flag;	
static   int        first[256];
static   unsigned char    one_mask[8];
	
public:
        AllocBlock(int sz) 
        { 
	  size  = sz;
	  memset(mask,0,32);                   // all blocks free
	  mxbyte = min(16,max(500000/sz/8,1)); // allocate at most 2MB
	  // cout << "Block: " << sz << " " << mxbyte << endl;
	  base = new IntegType[8*mxbyte*sz];
	  if (!base)
	    error("Memory Allocation Error");
	  mode = 0; // sequential handout
	  free_count = 8*mxbyte;  
	}
	AllocBlock(const AllocBlock& src)
        {
	  error("Illegal Constructor for AllocBlock");
	}
       ~AllocBlock()
	{
	  delete base;
	}
const AllocBlock& operator=(const AllocBlock& src)
        {
	  error("Illegal Assignment for AllocBlock");
	  return *this;
	}
int     error_check(int id)
        {
  	  // checking
	  int errflag = 0;  
	  for(int i = 0; i < 32; i++)
	    if (mask[i])
	    {
	      cout << "Error in Alloc: " << id << " byte: " << i 
		   << " mask: " << (int) mask[i] << endl;
	      errflag = 1;
	    }
	  return errflag;
	}
void    get(StorageUnit& su) // get a free pointer
        {
	  if (free_count == 0)
	    error("Not Storage free in AllocBlock. ");
	  
	  free_count--;
	  for(int i = 0; i < mxbyte; i++)
	    if (mask[i] != 0xff)
	    {
	      int bit = first[mask[i]];
	      su.id   = (i << 3) + bit;
	      mask[i] |= one_mask[bit];
	      su.base = base + size*su.id; 
	      return;
	    } 
	  error("AllocBlock::get -- Out of space");
	}
void    hand_back(StorageUnit& su)
        {
	  int byte = (su.id >> 3) & 0x1f;
	  //(ostream&) cout << su.id % 8 << " " <<  byte << " " << (su.id & 0x07)                          << endl;
	  if ((mask[byte] & one_mask[su.id & 0x07]) == 0)
	  {
	    cout << " Error in AllocBlock -- returning free block with id: " << 
	      su.id << endl;
	    abort();
	  }
	  mask[byte] &= ~one_mask[su.id & 0x07];  
	  free_count++;
	}	
int     free() const    { return free_count; }  
int     max_sz() const  { return 8*mxbyte; }  
static  void init();	
};

const int idd_shift =  16;
const int bid_shift =   8;
const int mx_blks   = 255;
class Allocator 
{
  friend class Storage;
  static int all_id;
         SStack<AllocBlock*,mx_blks> block;
         int last_free;
  
	 int idd;
         int sz[5];
         int tot_sz;
public:  
	 Allocator(int size1,int size2,int size3 = 1,int size4 = 1,int mxaux = 1);
	 Allocator(const Allocator& src);
        ~Allocator();

const Allocator& operator=(const Allocator& src)     
{
  error("Illegal Assignment for Allocator");
  return *this;
}

void get(StorageUnit& su)
{
#ifdef MEM_SAVE
  su.base = new IntegType[tot_sz];
  return;
#else
  AllocBlock* ab = block[last_free];
  
  ab -> get(su);
  // cout << "AID: " << idd << " BID: " << last_free << " SID: " << su.id;  
  su.id |= (idd << idd_shift) | (last_free << bid_shift); // build id
  // cout << " ID: " << su.id << endl;
  
  // search for the next free block, make one if neccessary

  if (!(ab -> free()))
  {
    for(int i = last_free + 1; i < block.size(); i++)
      if (block[i] -> free())
      {
	last_free = i; // found a free block
	return;
      }
    
    if (block.size() == block.max_size()) // stack full
    {
      error("need more blocks. got not ids");
      Array<AllocBlock*> tmp(block.size());
      int i;
      for(i =0; i < block.size(); i++)
	tmp[i] = block[i];
      block.resize(block.size() + 1000); // make more blocks
      for(i = 0; i < tmp.size(); i++)
	block.add(tmp[i]);
    }
    
    // make a new block

    ab = new AllocBlock(tot_sz);
    if (!ab) 
      error("Allocation Error for AllocBlock. ");
    last_free = block.add(ab);
  }   
#endif
}
	

void hand_back(StorageUnit& su)
{
#ifdef MEM_SAVE
  delete su.base;
  return;
#else
  int bno = (su.id >> bid_shift) & 0xff;

  // cout << "returning: " << su.id << " ID: " << (su.id >> idd_shift) 
  // << " BLOCK: " << bno << " intern: " << (su.id & 0xff) << endl;
  
  // first check if this block actually belongs to this alloc
  if ((su.id >> idd_shift) != (idd & 0xffff))
  {
    cout << "Error in HandBack -- ID: " << (su.id >> idd_shift) << " expected: " 
         << idd << endl;
    abort();
  }
  
  // determine the block no
    
  if (bno > block.size() || bno < 0)
  {
    cout << "Error in HandBack -- ID: " << (su.id >> idd_shift) 
      << " bno: " << bno 
      << " size is: " << block.size() << endl; 
    abort();
  }
  block[bno] -> hand_back(su);

  if (bno < last_free)
    last_free = bno;
#endif
}

int      check();  
int      size1()    { return sz[0]; }  
int      size2()    { return sz[1]; }  
int      size3()    { return sz[2]; }  
int      size4()    { return sz[3]; }
int      blocks()   { return sz[4]; }  
int      size()     { return tot_sz; } 
int      id() const { return idd; } 	 
};



class Storage  : public StorageUnit
{
static long    c_number;
static long    c_size;
static long    c_notot;
static long    c_totsize;
static long    c_sizehigh;
static long    c_nohigh;

protected:
  Allocator&    alloc;
  
  static NUMARRAY<IntegType> work1;
  static NUMARRAY<IntegType> work2;
  static NUMARRAY<IntegType> work3;
  
public:
  Storage(const Storage& scr) : alloc(scr.alloc)   
  { error("ILLEGAL CONSTRUCTOR FOR STORAGE"); }
  Storage(const Storage& source,int block_no, Allocator& al);
 ~Storage()   { c_number--; c_size -= size(); alloc.hand_back(*this); }  
  Storage(Allocator& al);
  Storage(const Storage& scr,Allocator& newalloc);
  
void add_to        (IntegType fac,const Storage& src);
void add_to        (IntegType fac,const Storage& coef,const Storage& src);
void add_invers_to (IntegType fac,const Storage& coef,const Storage& src);
void doublet(IntegType fac,const Storage& coef,const Storage& src,int pos);
void doublet(const Storage& coef,const Storage& src,int pos)
{
//  timer -> start(t_doublet);
  int mx12 = size12();
  int mx34 = size34();
  int max  = mx12 * mx34;
  int count = 0,i,j,b;
  
  switch(pos)
  {
  case 0:
    for(b=0; b < alloc.sz[4]; b++)
    {
      for(i=0; i< mx12; i++)
      {
	int max  = mx34;
	IntegType* target1 = base     + count;
	IntegType* source1 = src.base + count;
	IntegType  totfac  = coef.base[i]; 
#ifdef RS6000
	int half = max / 2;
	IntegType* target2 = target1 + half;
	IntegType* source2 = source1 + half;
	for(j= 0; j < half; j++)
	{
	  target1[j] += totfac*source1[j];
	  target2[j] += totfac*source2[j];
	}
	if (2*half != max)
	  target2[j] += totfac*source2[j];
#else
	for(j= 0; j < max; j++)
	  target1[j] += totfac*source1[j];
#endif
	count+=max;
      }
    }
    break;
  case 1:
    for(b=0; b < alloc.sz[4]; b++)
    {
      // int half = mx12/2;
      
      for(j=0; j< mx34; j++)
      {
	IntegType* target1 =     base + j + count;
	IntegType* source1 = src.base + j + count;
	IntegType totfac   = coef.base[j];	  	  
#ifdef RS6000
	IntegType* target2 = target1 + half*mx34;
	IntegType* source2 = source1 + half*mx34;
	max = half*mx34;
	for(i=0; i< max; i+=mx34)
	{
	  target1[i] += totfac*source1[i];
	  target2[i] += totfac*source2[i];
	}
	if (2*half != mx12)
	  target2[i] += totfac*source2[i];
#else
	for(i=0; i< max; i+=mx34)
	  target1[i] += totfac*source1[i];
#endif
      }
      count += mx12*mx34;
    }
    break;
  default:
    error("Illegal position in Storage::IntegTypet");
  }
//  timer -> stop(t_doublet);
}

void normalize(int m1,int m2,int m3,int m4,Quadruple_Set& qset);
void fold(Storage& coef,Storage& src,int noblocks);

Storage& operator*=(IntegType fac) 
{
  for(int i =0; i < size(); i++)
    base[i] *= fac;
  return *this;
}

void c_stats()
{ 
  c_notot++; 
  c_number++; 
  if (c_number > c_nohigh) c_nohigh = c_number;
  c_size    += size();
  c_totsize += size();
  if (c_size > c_sizehigh) c_sizehigh = c_size;
}
const Storage& operator=(const Storage& s)
{
  if (size1() != s.size1() || size2()  != s.size2() || size3() != s.size3() || 
      size4() != s.size4() || blocks() != s.blocks())
    error("Illegal Copy for Storage. "); 
  memcpy(base,s.base,sizeof(IntegType)*size());
  return *this;
}

IntegType*  get_base()               { return base;    }
IntegType   operator()(int p) const  { return base[p]; } 
IntegType&  operator[](int p)        { return base[p]; } 
void        zero()     { memset(base,0,sizeof(IntegType)*size()); }      
  
int size ()      const  { return alloc.size(); } 
int size1()      const  { return alloc.sz[0]; } 
int size2()      const  { return alloc.sz[1]; } 
int size3()      const  { return alloc.sz[2]; } 
int size4()      const  { return alloc.sz[3]; } 
int size12()     const  { return alloc.sz[0]*alloc.sz[1];  } 
int size34()     const  { return alloc.sz[2]*alloc.sz[3];  } 
int blocks()     const  { return alloc.sz[4]; } 

void  to_matrix  (NUMARRAY<IntegType>& m); 
void  from_matrix(NUMARRAY<IntegType>& m,int size1,int size2,
		  int size3,int size4,int mxaux);
void  block_copy(int target_block,Storage& source,int source_block);

void output(InternalBasis& b1,SymDesignator& s1,
	    InternalBasis& b2,SymDesignator& s2,
	    InternalBasis& b3,SymDesignator& s3,
	    InternalBasis& b4,SymDesignator& s4,
	    int ll1,int moff1,int ll2,int moff2,
	    int ll3,int moff3,int ll4,int moff4,
	    IntegFile& obuf,int block_no,
	    int same12,int same34,int samepair);
void output(InternalBasis& b1,SymDesignator& s1,
	    InternalBasis& b2,SymDesignator& s2,
	    int l1,int moffset1,int l2,int moffset2,
            ARRAY<IntegMat>& m,int mu_index,int same,int m_sym = 1);
static void print_stats();
static void clear_work() { work1.reset(0); work2.reset(0); work3.reset(0); }
};

ostream& operator<<(ostream& os,const Storage& s);
#if defined(SEPARATE_OSTREAM_WITH_ASSIGN)
ostream_withassign& operator<<(ostream_withassign& os,const Storage& s);
#endif
#endif
