/*TEX
%
% IPACK - QCHEM
% Quantum Chemistry Project: Integral Package
% Copyright (C) 1994 : Wolfgang Wenzel, University of Dortmund
%
% This program is propriatary software. Unlicensed use and distribution
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
% $Id: itest.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
*/
#include "io_incl.h"
#include <numer.h>
#include <math.h>
#include <assert.h>

/*
   This program is used to test two-electron integrals agains VMOL.

   Input is: The ASCII   file from VMOL

             INPUT file of integrals to be tested in the same format as ASCII

*/

int MAX = 1000000;  // maximum number of integral in ASCII !


const shift         = 8;
const int BBUFSIZE  = 600;

typedef struct 
{
  int     no;
  double  value[BBUFSIZE];
  int     compind[BBUFSIZE];
} twoel_buffer;

// binary search

typedef struct 
{
  unsigned long  idx;
  double value;
  char   count;
} IVal;


long find(unsigned long idx,IVal* list,long max)
{
  long mn = 0;
  long mx = max;

  while(mx-mn>1)
  {
    int nxttry = (mn+mx) / 2;
    if (list[nxttry].idx < idx)
      mn = nxttry;
    else 
      mx = nxttry;
  }
  if (list[mx].idx == idx) return mx;
  if (list[mn].idx == idx) return mn;
  
  return -1;
}

inline void swap(unsigned long& a,unsigned long& b)
{
  unsigned long tmp = a;
  a = b;
  b = tmp;
}

int compare(IVal* a, IVal* b)
{
  if (a -> idx >   b -> idx) return 1;
  if (a -> idx ==  b -> idx) return 0;
  return -1;
}

typedef int (*COMP)(const void*,const void*);
     
main(int argc,char** argv)
{
  if (argc < 3)
  {
    cout << " +++ ITEST Usage:  itest `asciifile' `map' " << endl;
    exit(1);
  }
  
  
  int i,j,k;
  ifstream in2(argv[2],ios::in);
  if (!in2)
  {
    cout << " Could not find MAP file: " << argv[2] << endl;
    exit(1);
  }
  
  int no_orb;
  in2 >>  no_orb;

  // if           n    is the IPACK orbital number then
  //    translate[n]   is the VMOL  orbial  number

  int* translate = new int[no_orb];
  for(i=0; i<no_orb;i++)
  {
    in2 >> translate[i];
    cout << " Translate: " << i << " to: " << translate[i] << endl;
  }

  IVal* ival = new IVal[MAX];
  fstream in1(argv[1],ios::in);
  if (!in1)
  {
    cout << " Could not find ASCII file: " << argv[1] << endl;
    exit(1);
  }

   // skip single particle operators

   double dummy;
   in1 >> k;

   Matrix snew;
   snew.read("SMATRIX");
   Matrix coreham;
   coreham.read("COREHAM");
 
   // S and H are the matrices from the ASCII file
  
   Matrix s(k,k);
   Matrix h(k,k);

   // read overlap

   for(i=0; i<k;i++)
     for(j=0; j<=i;j++)
     {
       in1 >> dummy;
       s(i,j) = dummy;
       s(j,i) = dummy;
     }

   // norm are the norms of the ASCII file
   
   Vector norm(s.size1());
   for(i=0; i< s.size1();i++)
     norm[i] = 1.0/sqrt(s(i,i));
   norm.print(cout,5," %15.10f");
   cout << endl;
   
   for(i=0; i<k;i++)
     for(j=0; j<=i;j++)
     {
       s(i,j) *= norm[i]*norm[j];
       if (i != j)
	 s(j,i) *= norm[i]*norm[j];
     }
   
   int stest = 0;

  // testing the SMATRIX
  // i , j are indices of the S from SMATRIX
  // ii,jj are indices of the S from ASCII

  // cout << "SMATRIX: " << snew << endl;
  // cout << "ASCII  : " << s    << endl;
  
  for(i=0; i<k;i++)
    for(j=0; j<=i;j++)
    {
      int ii = translate[i];
      int jj = translate[j];
      if (fabs(snew(i,j)/sqrt(snew(i,i)*snew(j,j))-s(ii,jj)) > 1E-7)
      {
	(ostream&) cout << "Error in S: (i,j): " << i << " " << 
	  j << " (it,jt): " << ii << " " << jj << " S: " << s(ii,jj)
	    << " SNew: " << snew(i,j)/sqrt(snew(i,i)*snew(j,j)) 
	      << " " << s(ii,jj) - snew(i,j) << endl;
	stest++;
      }
    }

   cout << " +++ Tested S: " << stest << " Errors detected. " << endl;
   
   // skip coreham
   for(i=0; i<k;i++)
     for(j=0; j<=i;j++)
     {
       in1 >> dummy;
       h(i,j) = dummy*norm[i]*norm[j];
       if (i != j)
	 h(j,i) = dummy*norm[i]*norm[j];
     }
   
   int htest = 0;
   
   for(i=0; i<k;i++)
     for(j=0; j<=i;j++)
     {
       int ii = translate[i];
       int jj = translate[j];
       if (fabs(coreham(i,j)/sqrt(snew(i,i)*snew(j,j))-h(ii,jj)) > 1E-7)
       {
	 (ostream&) cout << "Error in H: (i,j): " << i << " " << 
	   j << " (it,jt): " << ii << " " << jj << " H: " << h(i,j)
	     << " HNEW: " << coreham(ii,jj)/sqrt(snew(i,i)*snew(j,j))
	     << " Deviation: " 
	     << coreham(ii,jj)/sqrt(snew(i,i)*snew(j,j))-h(i,j) << endl;
	 htest++;
       }
     }

   cout << " +++ Tested H: " << htest << " Errors detected. " << endl;

   in1 >> k;
   
   int n = 0;
   unsigned long i1,i2,i3,i4;
   while(in1)
   {  
     assert(n < MAX);
     in1 >> i1 >> i2 >> i3 >> i4 >> ival[n].value;
     // cout << i1 << " " << i2 << " " << i3 << " " << i4 << " "
     //      << ival[n].value << endl;
     
     if (i1 < i2)
       swap(i1,i2);
     if (i3 < i4)
       swap(i3,i4);
     if (i3 > i1 || (i1 == i3  && i4 > i2))
     {
       swap(i1,i3);
       swap(i2,i4);
     }

     ival[n].value *= norm[i1]*norm[i2]*norm[i3]*norm[i4];
     ival[n].count  = 0;   
     ival[n++].idx = (((((i1 << 8) + i2) << 8) + i3) << 8) + i4;
     /*
     cout << i1 << " " << i2 << " " << i3 << " " << i4 
          << " " << ival[n-1].idx   << " " << ival[n-1].value << endl;
     */
   }
   cout << " +++ Read    : " << n << " Integrals. " << endl;
   qsort(ival,n,sizeof(IVal), (COMP) compare);
   
    
   double tval;
   double d2;
   
   int count = 0;
   FILE* infile = fopen("TWOEL","rb");

   int totcount = 0;
   
   twoel_buffer buf;

   /* read header of file */
   int code = fread((char*) &buf,sizeof(buf),1,infile);
   assert(code == 1);

   while(1)
   {
     if (count == 0)
     {
       int code = fread((char*) &buf,sizeof(buf),1,infile);
       assert(code == 1);
       count = buf.no;
       if (count == -1) break;
     }
     if (count == -1) break;
     count--;
     
     
     i1 =  buf.compind[count]        & 0xff;
     i2 = (buf.compind[count] >>  8) & 0xff;
     i3 = (buf.compind[count] >> 16) & 0xff;
     i4 = (buf.compind[count] >> 24) & 0xff;

     double tval = buf.value[count] / 
       sqrt(snew(i1,i1)*snew(i2,i2)*snew(i3,i3)*snew(i4,i4));

     // (ostream&) 
     //  cout << i1 << " " << i2 << " " << i3 << " " << i4 << " " << tval << endl;
     
     // normal ordering:
     i1 = translate[i1];
     i2 = translate[i2];
     i3 = translate[i3];
     i4 = translate[i4];
     
     if (i1 < i2)
       swap(i1,i2);
     if (i3 < i4)
       swap(i3,i4);
     if (i3 > i1 || (i1 == i3  && i4 > i2))
     {
       swap(i1,i3);
       swap(i2,i4);
     }
     int index = (((((i1 << 8) + i2) << 8) + i3) << 8) + i4;
       
     int pos = find(index,ival,n);

     if (fabs(tval) < 1e-10)
     {
       if (pos == -1) 
	 continue; 
       else
	 cout << "Error: " << i1 << " " << i2 << " " << i3 << " " << i4 << " "
	   << " " << tval << " zero expected. " << endl;
     }
     
     if (pos == -1)
     {
       cout << " Error: " << i1 << " " << i2 << " " << i3 << " " << i4 << " "
	 << " " << tval << " Not found. " << endl;
       continue;       
     }

     double v = ival[pos].value;
     ival[pos].count++;
     if (fabs(v-tval) > 1e-5)
     {
       cout << " Error: " << i1 << " " << i2 << " " << i3 << " " << i4 << " "
	 << " " << tval;
       cout << " Error: expected: " << v <<  endl;
     }
     totcount++;
   }
   cout << " +++ Tested  : " << totcount << " Integrals. " << endl;
   
   int errcount = 0;
   for(i=0; i<n;i++)
   {
     if (ival[i].count != 1)
     {
       i1 = ival[i].idx         & 0xff;
       i2 = (ival[i].idx >>  8) & 0xff;
       i3 = (ival[i].idx >> 16) & 0xff;
       i4 = (ival[i].idx >> 24) & 0xff;
       cout << " Count: " << setw(4) << i1 << " " << setw(4) << i2 << " " 
	    << setw(4) << i3 << " " << setw(4) << i4 << " "
	 << " appeared: " << (int) ival[i].count << " times. Value: "
	 << setw(15) << ival[i].value << endl;
       errcount++;
       if (errcount > 100)
       {
	 cout << " ++ Maximum Number of Errors exceeded, quitting. " << endl;
	 break;
       }
     }
   }
   
	   
}







