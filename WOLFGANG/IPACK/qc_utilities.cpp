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
% $Id: qc_utilities.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
*/

#include <numer.h>
#include "io_incl.h"
#include "vtype.h"
#include "qc_utilities.h"

void header( const String& prog_nm, const String& v_numb, const String& date ) 
{

  cout << endl << endl;
  cout << "       **********************************************************\n";
  cout << "       ***                                                    ***\n";

  cout << "       ***   Q C H E M   -----";
  int pn_len = prog_nm.length();
  int i;
  for ( i=0; i<(15-pn_len)/2; i++ ) cout << " ";
  cout << prog_nm;
  for ( i=(15-pn_len)/2+pn_len; i<15; i++ ) cout << " ";
  cout << "----  V " << v_numb;
  int vn_len = v_numb.length();
  for ( i=7-vn_len; i<7; i++ ) cout << " ";
  cout << "   ***\n";

  cout << "       ***                                                    ***\n";

  cout << "       ***   " << date;
  int dt_len = date.length();
  for ( i=10-vn_len; i<10; i++ ) cout << " ";
  cout << "                                       ***\n";

  cout << "       ***   Copyright: Wolfgang Wenzel, University Dortmund  ***\n";
  cout << "       ***                                                    ***\n";
  cout << "       ***   Note: You need a license to run this program.    ***\n";
  cout << "       ***         Please consult the documentation on        ***\n";
  cout << "       ***         license and use.                           ***\n";
  cout << "       ***                                                    ***\n";
  cout << "       **********************************************************\n";
  cout << endl << endl << endl;

}  //  End of header()

/*  OLD CODE

int cmd_line_check( Array<String>& cmd_line, char* argv )
{

  int n_cmd_line = cmd_line.size();
  for ( int i=0; i<n_cmd_line; i++ )
    if ( !strcmp( argv, cmd_line[i].contents() ) ) return i;

  error( "cmd_line_check() : a cmd line argument not a possible argument" );

  return -1;

}  //  End of cmd_line_check()

*/

int cmd_flag(const char* token,int argc,char** argv)
{
  for(int i = 0; i < argc; i++)
    if(strcmp(token,argv[i]) == 0)
      return 1;
  return 0;
}

int input_skip(istream & is,char* token )
{
  // skip is until header sequence found 

  char c;
  while (is >> c)
  { 
    if (c == '#')
    {
      String tst;
      is >> tst;
      if ( tst == token ) return 1;
    }
  }
  return 0;
}  //  End of input_skip()

/*TEX
\subsection{Output Routines for Doubles}
*/
const char* sci3(double d)
{
  static char buf[50];
  sprintf(buf,"%10.3e",d);
  return buf;
}

const char* sci6(double d)
{
  static char buf[50];
  sprintf(buf,"%13.6e",d);
  return buf;
}


const char* sci15(double d)
{
  static char buf[50];
  sprintf(buf,"%20.15e",d);
  return buf;
}


const char* fix3(double d)
{
  static char buf[50];
  sprintf(buf,"%10.3f",d);
  return buf;
}

const char* fix6(double d)
{
  static char buf[50];
  sprintf(buf,"%13.6f",d);
  return buf;
}


const char* fix15(double d)
{
  static char buf[50];
  sprintf(buf,"%20.15f",d);
  return buf;
}

void print_symmetric_mat(const char* name,Mat& m,ostream& os)
{
  os << name << " " << setw(4) << m.size1() << endl;
  for(int i = 0; i < m.size1(); i++)
    for(int j = 0; j <= i; j++)
      os << setw(4) << i << " " << setw(4) << j << " " << fix15(m(i,j)) << endl;
  os << endl;
}


void read_symmetric_mat(Mat& m,istream& is)
{
  int n,ii,jj;
  double val;
  is >> n;
  if ((m.size1() != n) || (m.size2() != n))
    m.reset(n,n);

  if (!is)
    error("read_symmetric_mat: could not read matrix size ");

  for(int i = 0; i < m.size1(); i++)
    for(int j = 0; j <= i; j++)
    {
      is >> ii >> jj >> val;
      if (!is)
	error("read_symmetric_mat: input error ");
      if (ii != i) 
	error("read_symmetric_mat: first index mismatch");
      if (jj != j) 
	error("read_symmetric_mat: second index mismatch");
      m(i,j) = val;
      m(j,i) = val;
    }
}
//******************************************************************



/*TEX
\subsection{Comparison for Vectors and Matrices}

These functions are usuful for debugging, returning the number of error encountered.


*/


int    compare(Mat& a,Mat& b,double eps)
{
  if(a.size1() != b.size1())
  {
    cout << " Error in Matrix Compare: size1 does not match: A: " << a.size1() 
      << " B: " << b.size1() << endl;
    return 1;
  }
  if(a.size2() != b.size2())
  {
    cout << " Error in Matrix Compare: size2 does not match: A: " << a.size2() 
      << " B: " << b.size2() << endl;
    return 1;
  }
  
  int errcount = 0;
  for(int i =0 ;i < a.size1(); i++)
  for(int j =0 ;j < a.size2(); j++)
    if (fabs(a(i,j) - b(i,j)) > eps)
    {
      cout << " Error in Matrix Compare on Element: (i,j) " << i << " " << j 
 	   << " A: " << a(i,j) << " B: " << b(i,j) << endl;
      errcount++;
    }
  return errcount;
}

int    compare(Vec& a,Vec& b,double eps)
{
  if(a.size() != b.size())
  {
    cout << " Error in Vector Compare: size does not match: A: " << a.size() 
      << " B: " << b.size() << endl;
    return 1;
  }
  if(a.offset() != b.offset())
  {
    cout << " Error in Vector Compare: Offset does not match: A: " << a.offset() 
      << " B: " << b.offset() << endl;
    return 1;
  }
  int errcount = 0;
  for(int i = a.offset();i < a.bound(); i++)
    if (fabs(a(i) - b(i)) > eps)
    {
      cout << " Error in Vector Compare on Element: i " << i
	<< " A: " << a(i) << " B: " << b(i) << endl;
      errcount++;
    }
  return errcount;
}
