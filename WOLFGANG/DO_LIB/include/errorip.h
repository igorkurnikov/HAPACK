/*
%  Dortmund C++ Class and Template Library 
%  Copyright (C) 1994 Wolfgang Wenzel and others

%  The code in this file was derived from routines published in 
%  Numerical Recipes, it its the responsibility of the user to insure that
%  he is authorized to use them. Please read the appropriate section entitled
%  LICENSE in the documentation. 

%  This library is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%  Library General Public License for more details.

%  You should have received a copy of the GNU Library General Public
%  License along with this library; if not, write to the Free
%  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%
% $Id: errorip.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
*/

#ifndef ERROR_H
#define ERROR_H
#include <stdlib.h>
#include "io_incl.h"

inline void error(const char* msg )
{  cerr << "Error: " << msg << endl; abort(); }

inline void error(const char* msg, int i )
{  cerr << "Error: " << msg <<" : "<< i << endl; abort(); }

inline void error(const char* msg, const char* msg2 )
{  cerr << "Error: " << msg <<" : "<<  msg2 << endl; abort(); }

inline void warning(const char* msg )
{  cerr << "Warning: " << msg << endl; }

inline void warning(const char* msg,int i )
{  cerr << "Warning: " << msg <<" : "<< i << endl; }

#endif
