/*IGN
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
% $Id: stringc.cpp,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
%
*/
/*TEX
\section{Implementation}
*/
#include "io_incl.h"
#include <ctype.h>
#include <assert.h>
#include "stringc.h"
const char EOS   = '\0';

char* my_sdup(const char* str)
{
   int len = strlen(str);
   char* str_new = (char*) malloc(len+1);   
   strcpy(str_new,str);
   return str_new;
}

String::String()
{
	s = 0;
}

String::String(const char* ns)
{
	s = (ns != NULL) ? strdup(ns) : (char*) NULL; 
}

String::String(const String& ns)
{
	s = (ns.s != NULL) ? strdup(ns.s) : (char*) NULL;
}

String::~String() 
{ 
	if (s) delete s; 
}

String& 
String::operator=(char* ns)
{
	  if (s) delete s;
	  s = (ns == (char*) NULL) ? (char*) NULL : strdup(ns);
	  return *this;
}

String& 
String::operator=(const String& ns)
{
	  char *temp = (ns.s) ? strdup(ns.s) : (char*) NULL;
//	  if (s) delete s;  // IGOR CHANGED
      if (s) free(s);
	  s = temp;
	  return *this;
}


/*TEX
\subsection{String::first}
This function returns the first position in the string which matches the character 'c'. The search starts at position \name{pos}, which defaults to zero.
If the character is not found the functions returns -1.
*/
int    String::first(const char c,const int pos) const
{
  char* q = s;
  int i;
  if (q == NULL) return -1;
  for(i = 0; i < pos; i++)
    if (*(q++) == EOS) return -1;
  while(*q != EOS && *q != c) q++;
  if (*q == EOS) return -1;
  return q-s;
}
/*TEX
\subsection{String::firstof}
\name{Firstof} finds the first occurence of any character in the string 
\name{a} in the object to which the function is applied. The function 
returns -1, if none of the characters in \name{a} apprea in \name{*this}. 
*/  
int    String::firstof(const String& a,int pos) const
{
  if (a.s == NULL) return -1;
  if (s == NULL)   return -1;
  int fpos = -1;
  char *q = a.s;
  while(*q != EOS)
  {
    int tpos = first(*q,pos);
    if (tpos >= 0)
       if (fpos < 0) fpos = tpos;
       else fpos = (tpos < fpos) ? tpos : fpos;
    q++;
  }
  return fpos;
}
/*TEX
\subsection{String::substr}

\name{Substr} returns the substring starting at position \name{start}
inclusive up to position \name{pos}. Returns the \name{NULL} string is
start is beyond the length of the object.

*/
String String::substr(const int start,const int end) const
{
   int ss = start;
   char* q = s;
   int i;
   String tmp;
   if (q == NULL) return tmp;
   if (ss < 0) ss = 0;
   for(i = 0; i < ss; i++)
      if (*(q++) == EOS) return tmp;
   char *savepos = q;
 
   for(i=ss; i <= end; i++)
   {
      if (*(q++) == EOS) break;
   }
   char save = *q;
   *q  = EOS;
   tmp = savepos;
   *q  = save;
   return tmp;
}
/*TEX
\subsection{Concatenation}
\name{operator+} returns the concatenation of the two arguments. 
*/

String operator+(const String& a,const String& b)
{
  int newlen = a.length() + b.length()+1;
  char *ns = new char[newlen];
  if (a.s) 
  {
    strcpy(ns,a.s);
    if (b.s) strcat(ns,b.s);
  }
  else
  {
    if (b.s) strcpy(ns,b.s);
    else ns[0] = EOS;
  }
  String temp = ns;
  delete ns;
  return temp;
}


String operator+(const String& str,const int& number)
{
   char *ns;
   ns = new char[strlen(str.s) + sizeof(int) + 1];
   sprintf(ns,"%s%i",str.s,number);
   String temp = ns;
   delete ns;
   return temp;
}
    
/*TEX
\subsection{IO}
Note that the input functions is restricted to strings of maximal length of 
256 characters. Any \name{whitespace} terminates the string.

Note: I am sure fancier things can be done, akin to \name{get} of \name{ios} 
\index{String,input} 
*/

ostream& operator<<(ostream& os,const String& ss)
{
  return os << ((ss.s) ? ss.s : "NULL");
}

static const int maxlen = 256;

istream& operator>>(istream& is,String& ss)
{
  char buf[maxlen];  
  int cnt = 0;
  char c;
  
  while(is && is.get(c) && isspace(c))
  ;
  if (!is) return is;
  is.putback(c);
  
  while (cnt < maxlen-1 && is)
  {
    is.get(c);
    if (isspace(c))
    {
      is.putback(c);
      break;
    }
    buf[cnt++] = c;
  }

  buf[cnt] = '\0';
  if (is) ss = String(buf);
  
  return is;
  
}
/*TEX
\subsection{Binary IO}
*/
void String::write(ostream& file) const
{
  file << length() << ' ' << s;
}
void String::read(istream& file)
{
  int ll;
  file >> ll;
  if (!file) 
  {
    cerr << " Error on String::read. ";
    abort();
  } 
  if (ll < 0)
  {
    cerr << " Illegal length on String::read. " << ll << "\n";
    abort();
  }
  s = new char[ll+1];
  if (!s)
  {
    cerr << " ENOMEM on String::read. " << ll << "\n";
    abort();
  }
  char* p = s;
  file.get(); // skip the blank
  for(int i=0 ; i < ll && file; i++)
   *(p++) = file.get();
  *p = EOS;
  if (!file)
  {
    cerr << " File Error on String::read. " << ll << "\n";
    abort();
  }
}
/*TEX
\subsection{Comparisons}
\index{String,comparison} 
*/
int     operator==(const String& a,const String& b)
{
  if (a.s==NULL && b.s==NULL) return 1;
  if (a.s==NULL || b.s==NULL) return 0;
  if (a.length() != b.length()) return 0;
  int i;
  for(i=0; i< length(a);i++)
    if(a.s[i] != b.s[i]) return 0;
  return 1;
}

int     operator< (const String& a,const String& b)
{
  if (a.s == NULL) return 1;
  if (b.s == NULL) return 0;

  for(int i=0; i < length(a);i++)
   if (a.s[i] >= b.s[i]) return 0;
  return 1;
}

int operator!=(const String& a,const String& b) { return !(a==b);}
int operator<=(const String& a,const String& b) { return !(b<a);}
int operator>=(const String& a,const String& b) { return !(a<b);}
int operator> (const String& a,const String& b) { return (b<a); }

