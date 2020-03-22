/*TEX
%
% IPACK - 
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
% $Id: integ_file.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
\section{Joint Single-Electron and Two-Electron Integral IO\label{section-IO}}
\index{Input, two-electron} 
\index{Output, two-electron} 
\index{Dump, two-electron} 
\index{File, two-electron}
\index{Buffer, two-electron}
*/

#ifndef INTEG_FILE_H
#define INTEG_FILE_H

#include <math.h>
#include <stdio.h>
#include <stringc.h>
#include "flexfile.h"
#include "parallel.h"

#if !defined(GNU) && !defined(_MSC_VER) && !defined(__DECCXX)
#ifndef BOOLEAN 
#define BOOLEAN
enum Boolean  { false, true  };
#endif
#endif
  
typedef struct 
{
  unsigned  short int  idx1;
  unsigned  short int  idx2;
  float                 val;
} IntegData;

ostream& operator<<(ostream& os,IntegData data);
istream& operator>>(istream& is,IntegData data);

int core_file_bucket_size(int no_orbs);

const  int integ_file_twoel_pip_offset = 256;
//const  int integ_file_shift            =  16;
enum   IntegFilePipCodes { IntegInfoC, IntegOvlp, IntegNuclear, IntegKinetic, 
			   IntegAngular    =  4,
			   IntegEField     =  7,
                           IntegSpinOrb    = 10,
			   IntegECP        = 13,
			   IntegCore       = 14,
			   IntegDiagDir    = 18,
			   IntegDiagExc    = 19,
			   IntegRhoUp      = 20,
			   IntegRhoDn      = 25,
			   IntegTrans      = 30,
                           IntegTransTotal = 31,
			   IntegMoments    = 64 };

 
const double epsilon     = 1E-10;
  
class  IntegFile : virtual public BucketFile<IntegData>
{
protected:
  unsigned int              i_temp;
  IntegData                   data;
  ARRAY<unsigned int>  orbital_map;
public: 
        IntegFile(int no_orbs);
       ~IntegFile();

void          set_orbital_map(ARRAY<unsigned int>& map) { orbital_map = map; } 
void          reset_orbital_map();
void          sp_hamiltonian(IntegMat& m);  
void          write_mat(IntegMat& mat,int target_pip);
int           read_mat (IntegMat& mat,int target_pip,int clear = false);
int           no_orbs() const { return orbital_map.size(); }
void          dump(unsigned int ii1, unsigned int ii2, unsigned int ii3,
	           unsigned int ii4, double val)
{  
  if (fabs(val) < epsilon ) return; 

  unsigned int i1 = orbital_map(ii1);
  unsigned int i2 = orbital_map(ii2);
  unsigned int i3 = orbital_map(ii3);
  unsigned int i4 = orbital_map(ii4);

  if ( i1 < i2 ) { i_temp = i1; i1 = i2; i2 = i_temp; }
  if ( i3 < i4 ) { i_temp = i3; i3 = i4; i4 = i_temp; }
  if ( i1 < i3 ) { i_temp = i3; i3 = i1; i1 = i_temp;
                   i_temp = i4; i4 = i2; i2 = i_temp; }
  if ( i1 == i3 && i4 > i2 ) { i_temp = i4; i4 = i2; i2 = i_temp; }

  int pip   = i1*(i1+1)/2 + i2 + integ_file_twoel_pip_offset;
  data.idx1 = i3;
  data.idx2 = i4; 
  data.val = val;

  add(pip,data);

}
virtual void write_file(const char* fname) {};
virtual void read_file (const char* fname) {};  
virtual void clean() { }
virtual void write_all() = 0;
virtual int  local(int key) { return 1; }
virtual int  addr (int key) { return 0; }
};

class IntegDiskFile :  public FlexFile<IntegData>, public IntegFile
{
public:
             IntegDiskFile(const char* fname,int no_orbs,int debug = 0);
             IntegDiskFile(const char* fname);
virtual     ~IntegDiskFile();  
virtual void write_all()
             {
	       FlexFile<IntegData>::finalize();
	     }
};


class IntegCoreFile : public IntegFile,public CoreFile<IntegData>
{
public:
             IntegCoreFile(const char* fname,int no_orbs,int debug = 0);
            ~IntegCoreFile();
virtual void write_file(const char* fname) 
             { CoreFile<IntegData>::write_file(fname,0); };
virtual void read_file (const char* fname) 
             { CoreFile<IntegData>::read_file (fname); };
virtual void clean()  { CoreFile<IntegData>::clean();  }
virtual int  local(int key) { return CoreFile<IntegData>::local(key); } 
virtual int  addr (int key) { return CoreFile<IntegData>::addr(key); } 
virtual void write_all()
             {
	       CoreFile<IntegData>::write_all_basic();
	     }
};

#endif










