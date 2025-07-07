#ifndef BITVEC_DEF
#define BITVEC_DEF
#include <assert.h>
#include <numer.h>
/*TEX \\Bitvector class. It is a template with an integer value size that
 gives the number of bytes of the bitvector. The number of bits is
8*size. Bit number 0 is the right most.The first byte store 0..7 the
second 8..15 and so on\\        
\begin{tabular}{p{3cm}p{10cm}}
Bitvector()&creates an zero bitvector\\
operator (int i) & returns the i'th bit\\
set(int i) & sets the i'th bit\\
unset (int i) & unsets the i'th bit\\ 
reset() &resets the vector, sets it to zero\\
size() & returns the number of bits\\
char get_byte(int i) & return the ith byte of the bitvector\\
num\_set()&returns the number of set bits\\
num\_set\_left_from(int pos)&returns the number of set bits
                                       left from bit position pos\\
D\_IVectro list\_of\_set\_bits()&returns a list of the set bits\\
operators ==, !=, <,> & are defined. 10000 is bigger than 01000\\
ouput operator << & is defined for bitvectors
\end{tabular}\\ 
*/


template <int Size>
class Bitvector
{
 public:
  Bitvector();
  int operator() (int i) const 
  {
    assert((i>=0)&&(i<Size*8));
    return((vector[v_byte[i]]>>v_bit[i])&1);
  }
  void set (int i)
  { 
    assert((i>=0)&&(i<Size*8));
    vector[v_byte[i]]|=v_setmask[v_bit[i]];
  }
  void unset (int i)
  {
    assert((i>=0)&&(i<Size*8));    
    vector[v_byte[i]]&=v_unsetmask[v_bit[i]];
  }
  void reset();
  int size() const {return Size*8;}
  int num_set() const; 
  D_IVector set_list() const; 
  int num_set_left_from(int pos) const; 
  char get_byte(int i) const {assert(i>=0&&i<Size); return vector[i];}    
  friend int operator==(const Bitvector<Size>& a,const Bitvector<Size>& b);
  friend int operator!=(const Bitvector<Size>& a,const Bitvector<Size>& b); 
  friend int operator<(const Bitvector<Size>& a,const Bitvector<Size>& b); 
  friend int operator>(const Bitvector<Size>& a,const Bitvector<Size>& b); 
 private:
  char vector[Size];
  static char v_byte[8*Size];
  static char v_bit[8*Size];
  static char v_setmask[8];
  static char v_unsetmask[8];
  static char v_num_of_set_bits[256];
  static D_IVector v_list_of_set_bits[256];
  static char v_num_of_ones_left[256][8];
  static char v_static_set;  
 };

template <int Size>
Bitvector<Size>::Bitvector()
{ 
  if (!v_static_set)
  {
    for (int i=0;i<8*Size;i++) 
    {
      v_byte[i]=i/8;
      v_bit[i]=i%8;
    }
    int i0,i1,i2,i3,i4,i5,i6,i7;
    int index,number;
    for (i0=0;i0<2;i0++)
      for (i1=0;i1<2;i1++)
	for (i2=0;i2<2;i2++)
	  for (i3=0;i3<2;i3++)
	    for (i4=0;i4<2;i4++)
	      for (i5=0;i5<2;i5++)
		for (i6=0;i6<2;i6++)
		  for (i7=0;i7<2;i7++)
		  {
		    index=i0+2*i1+4*i2+8*i3+16*i4+32*i5+64*i6+128*i7;
		    number=i0+i1+i2+i3+i4+i5+i6+i7;
		    v_num_of_set_bits[index]=number;
		    v_list_of_set_bits[index]=D_IVector(number);
		    int i=0;
		    if (i7) v_list_of_set_bits[index][i++]=7;
		    if (i6) v_list_of_set_bits[index][i++]=6;
		    if (i5) v_list_of_set_bits[index][i++]=5;
		    if (i4) v_list_of_set_bits[index][i++]=4;
		    if (i3) v_list_of_set_bits[index][i++]=3;
		    if (i2) v_list_of_set_bits[index][i++]=2;
		    if (i1) v_list_of_set_bits[index][i++]=1;
		    if (i0) v_list_of_set_bits[index][i++]=0;
		    v_num_of_ones_left[index][0]=i1+i2+i3+i4+i5+i6+i7;
		    v_num_of_ones_left[index][1]=i2+i3+i4+i5+i6+i7;
    		    v_num_of_ones_left[index][2]=i3+i4+i5+i6+i7;
    		    v_num_of_ones_left[index][3]=i4+i5+i6+i7;
    		    v_num_of_ones_left[index][4]=i5+i6+i7;
    		    v_num_of_ones_left[index][5]=i6+i7;
    		    v_num_of_ones_left[index][6]=i7;
    		    v_num_of_ones_left[index][7]=0;
		  }
    v_static_set=1;
  }
  for (int i=0;i<Size;i++) vector[i]=0;
}

template <int Size>
void Bitvector<Size>::reset()
{ 
  for (int i=0;i<Size;i++) vector[i]=0;
}

template <int Size>
int Bitvector<Size>::num_set() const 
{
  int number=0;
  for (int i=0;i<Size;i++) number+=v_num_of_set_bits[vector[i]];
  return number;
}

template <int Size>
int Bitvector<Size>::num_set_left_from(int pos) const 
{
  int i,number=0;
  for (i=Size-1;i>v_byte[pos];i--) number+=v_num_of_set_bits[vector[i]];
  number+=v_num_of_ones_left[vector[v_byte[pos]]][v_bit[pos]];
  return number;
}

template <int Size>
D_IVector Bitvector<Size>::set_list() const 
{
  int i,j;
  D_IVector list,help;
  for (i=Size-1;i>=0;i--) 
    {
      help=v_list_of_set_bits[vector[i]];
      for (j=0;j<help.size();j++) help[j]+=i*8;
      list.append(help);
    }
  return list;
}

template <int Size>
int operator==(const Bitvector<Size>& a,const Bitvector<Size>& b)
{
  for(int i=0;i<Size;i++) if (a.vector[i]!=b.vector[i]) return 0;
  return 1;
}

template <int Size>
int operator!=(const Bitvector<Size>& a,const Bitvector<Size>& b)
{
  for(int i=0;i<Size;i++) if (a.vector[i]!=b.vector[i]) return 1;
  return 0;
}

template <int Size>
int operator>(const Bitvector<Size>& a,const Bitvector<Size>& b)
{
  for(int i=Size-1;i>=0;i--) 
  {
    if (a.vector[i]>b.vector[i]) return 1;
    if (a.vector[i]<b.vector[i]) return 0;
  }
  return 0;
}

template <int Size>
int operator<(const Bitvector<Size>& a,const Bitvector<Size>& b)
{
  for(int i=Size-1;i>=0;i--) 
  {
    if (a.vector[i]<b.vector[i]) return 1;
    if (a.vector[i]>b.vector[i]) return 0;
  }
  return 0;
}

template <int Size>
ostream & operator <<(ostream & os, const Bitvector<Size> & b)
{
  for (int i=0;i<Size*8;i++)
  {
    if (b(Size*8-i-1)) os <<"1";
    else os<<"0";
    if (((i+1)%70)==0) os<<"\n";
  }
  os<<"\n";
  return os;
}



template <int Size>
char Bitvector<Size>::v_setmask[8]={1,2,4,8,16,32,64,128}; 
 
template <int Size>
char Bitvector<Size>::v_unsetmask[8]={254,253,251,247,239,223,191,127};

template <int Size>
char Bitvector<Size>::v_num_of_set_bits[256];

template <int Size>
D_IVector Bitvector<Size>::v_list_of_set_bits[256];

template <int Size>
char Bitvector<Size>::v_num_of_ones_left[256][8];

template <int Size>
char Bitvector<Size>::v_byte[8*Size];

template <int Size>
char Bitvector<Size>::v_bit[8*Size];


template <int Size>
char Bitvector<Size>::v_static_set=0;





#endif







