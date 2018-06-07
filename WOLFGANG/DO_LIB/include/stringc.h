/*IGN
%  Dortmund C++ Class and Template Library 
%  Copyright (C) 1994 Wolfgang Wenzel

%  This library is free software; you can redistribute it and/or
%  modify it under the terms of the GNU Library General Public
%  License as published by the Free Software Foundation; either
%  version 2 of the License, or (at your option) any later version.

%  This library is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%  Library General Public License for more details.

%  You should have received a copy of the GNU Library General Public
%  License along with this library; if not, write to the Free
%  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
  
%  $Id: stringc.h,v 1.1.1.1 2008/04/08 19:44:31 igor Exp $
*/
/*TEX
\chapter{Strings}
\index{String}

This class implements a simple, self-contained string class, which 
simply encapsulates the standard c-implementation of strings as
zero-terminated arrays of \name{char}. The implementation hides
the associated memory management from the user, who is therefore
protected from related common mistakes.

\section{Class Description}

\index{String,description}   
In the following the main features of the class \name{String}
are summarized:

\begin{itemize}
\item{ {\bf Includes: } The file {\bf"stringc.h"} must be included. 
       \index{String,include}}
\item{ {\bf Constructors: \index{String,constructors}} 
  \begin{itemize}
  \item{A \name{String} is constructed with the default constructor
        \name{String s; } Character strings may later be assigned
        to this variable. }
  \item{A \name{String} is constructed from a {\bf char* cstring;} 
        standard null-teriminated c-string as: \name{String s(cstring); }
        The c-string will be copied. }
  \item{A \name{String} is constructed from another \name{String}
        using the copy-constructor. }
  \end{itemize} }
\item{ {\bf Copy \& Assignment: }
  \begin{itemize}
  \item{\name{Strings} can be assigned to each other.}
  \item{c-strings can be assigned to \name{Strings}.}
  \item{\name{Strings} can be converted to c-strings using the 
        \name{operator()} or the member
        function \name{contents}, which return a pointer to a 
        {\bf copy} of the contents of the \name{String}. If the 
        string is empty, zero is returned. Take care with the memory,
        that is allocated for the copy, e.g. delete it after use!}
  \item{\name{Strings} can be concatenated using the \name{+} operator}
  \item{An integer number can appended to a \name{String} by the
        \name{+} operator e.g. \name{String1=String2+integer}.}
  \end{itemize}
\item{ {\bf Access: } 
  \begin{itemize}
  \item{The \name{operator[]} allows retrieval and assignments of 
        individual characters. Full index checking is implemented,
        the program will abort if an index beyond the current
        allocated length is used.}
   \item{The constant \name{operator()(int pos)} allows retrieval 
        of individual characters. Full index checking is implemented,
        the program will abort if an index beyond the current
        allocated length is used.}
  \item{The current allocated length can be retrieved through the
        member function \name{length} or through the function 
        \name{length(const String\& s)}}
  \item{The member function \name{substr(a,b)} returns the substring
        starting from position $a$ to position $b$ inclusive. If $a < 0$ the 
        substring starts at the beginning. If $a$ is beyond the 
        last character of the string or $b < a$, an empty string 
        is returned. If $b$ is positioned beyond the end of the string,
        the substring to the end of the string is returned. }
  \item{The member-function \name{first(c,p)} returns the index of the 
        first occurence of the character $c$, starting the search from
        position p. If the character cannot be found, the function 
        returns -1. The parameter p is optional, if omitted the search 
        starts at the beginning.}
  \item{The member function \name{contains(c)} returns 1 if the character
        is present in the string, 0 otherwise. }
  \item{The member function \name{firstof(stringb,pos)} searches for 
        the first occurence of any character in \name{stringb}, starting
        at position \name{pos}. The parameter \name{pos} is optional.
        The function returns -1, if \name{stringb} is empty of if
        none of the characters in \name{stringb} can be found. 
	}
  \end{itemize} }
\item{ {\bf Input Output: } 
  \begin{itemize}
  \item{Both the input and output stream operators are overloaded. The
        input operator skips all whitespace, and then reads a maximum of
        256 characters until the next white-space is encountered. Whitespace
        is defined by the standard \name{isspace} function. }
  \item{``Binary'' read and write functions also exist, which ignore white
        space entirely. Instead the number of characters is written/read
        preceeding the string. }
  \item{{\bf Not implemented:} The \name{scan(const String\& term)} function
        scans up to 256 characters from the input stream until one of the 
        characters in \name{term} is reached. If succesful it returns the 
        index of the character found, otherwise -1. A similar function
        \name{scan(istream\& is,char t)} implements the scan for a single
        terminating character.}
  \item{ {\bf Not implemented:}
        The \name{skipto(iostream\& is,const String\& s)}
        skip the characters in the the stream until one character in \name{s}
        is found, which is also skipped.  Return values as in \name{scan}.}
  \item{ {\bf Not implemented:}
        The \name{skipto\_word(iostream\& is,const String\& word)}
        skip the characters until \name{word} is found. Returns 0 on 
        success, -1 on failure. }
  \end{itemize} }
\item{ {\bf Comparison: } 
  \index{String,comparison} 
  All lexical comparison operators have been implemented. Strings
  are compared from left to right according to the ASCII code of 
  their characters. If \name{s1} is the initial substring of \name{s2}
  then $s1 < s2$ Two empty strings are defined as equal. } 
  }
  
\end{itemize}


\section{Class Header}
*/
#ifndef __LIB_STRING_H
#define __LIB_STRING_H 

//class ostream;
//class istream;

#include <string.h>
#include "io_incl.h"
#include <stdlib.h>

#undef strdup
#define strdup my_sdup
extern char* my_sdup(const char* str);

class String
{
  char* s;
  public:
        String(); 
        String(const char* ns);         
        String(const String& ns);
       ~String();

String& operator=(char* ns);
String& operator=(const String& ns);
 
char&   operator[](int i)
        {
	  if (i >= 0 && i < strlen(s)) return *(s+i);
	  else 
	  {
	    cerr << "String::Operator[] -- Invalid Index to String: " << *this
	      << " was: " << i << endl;
	    abort();
	  }
	  return *s; // never reached
	}

char    operator()(int i) const
        {
	  if (i >= 0 && i < strlen(s)) return *(s+i);
	  else 
	  {
	    cerr << "String::Operator[] -- Invalid Index to String: " << *this
	      << " was: " << i << endl;
	    abort();
	 }
	  return *s; // never reached
	}
  
String  substr(const int start,const int end) const;
int     first(const char c,const int pos = 0) const;
int     contains(char c) const { return (first(c) >= 0); }
int     firstof(const String& a,int pos = 0) const;

void    write(ostream& file) const;
void    read (istream& file);
friend  istream& operator>>(istream& os,String& ss);
friend  ostream& operator<<(ostream& os,const String& ss);

friend  String  operator+(const String& a,const String& b);
friend  String  operator+(const String& str,const int& number);
 
char*   operator()() const { return (s) ? strdup(s) : (char*) NULL; }
char*   contents()   const { return (s) ? strdup(s) : (char*) NULL; }

int     length() const { return (s) ? strlen(s) : 0; }

friend  int     operator==(const String& a,const String& b);
friend  int     operator< (const String& a,const String& b);
friend  int     operator!=(const String& a,const String& b);
friend  int     operator<=(const String& a,const String& b);
friend  int     operator>=(const String& a,const String& b);
friend  int     operator> (const String& a,const String& b);
};

inline int length(const String& a)              { return a.length(); }

#endif
