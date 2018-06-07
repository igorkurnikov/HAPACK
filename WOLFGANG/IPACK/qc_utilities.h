/*TEX
%
% IPACK - QC_UTILITIES
% Quantum Chemistry Project: Quantum Chemistry Utilities
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
% $Id: qc_utilities.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
*/

#ifndef QC_UTILITIES_H
#define QC_UTILITIES_H

#include <stringc.h>
#include <tarray.h>
#include <errorip.h>

#if !defined(GNU) && !defined(_MSC_VER) && !defined(__DECCXX)
#ifndef BOOLEAN
#define BOOLEAN
enum Boolean   { false, true };
#endif 
#endif 

void  header        (const String& prog_nm,const String& v_numb, 
		     const String& date);
int   input_skip    (istream & is, char* token );
int   cmd_flag      (const char* tokn,int argc,char** argv);
 
const char* fix3(double a);
const char* fix6(double a);
const char* fix15(double a);

const char* sci3(double a);
const char* sci6(double a);
const char* sci15(double a);

void        print_symmetric_mat(const char* name,Mat& m,ostream& os);
void        print_symmetric_mat(const char* name,IntegMat& m,ostream& os);
void        read_symmetric_mat(Mat& m,istream& os);
void        read_symmetric_mat(IntegMat& m,istream& os);

inline int  sign(const int i) { return 1 - 2*(i % 2); }
inline int  min(int a,int b)  { return (a < b) ? a : b; } 
inline int  max(int a,int b)  { return (a < b) ? b : a; } 

int         compare(Mat& a,Mat& b,double eps);
int         compare(Vec& a,Vec& b,double eps);

#endif

