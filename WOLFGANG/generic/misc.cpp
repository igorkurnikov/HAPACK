#include <stdio.h>
 
/*TEX
\subsection{Ooutput Routines for Doubles}
*/
const char* sci3(double d)
{
  static char buf[20];
  sprintf(buf,"%10.3e",d);
  return buf;
}

const char* sci6(double d)
{
  static char buf[20];
  sprintf(buf,"%13.6e",d);
  return buf;
}


const char* sci15(double d)
{
  static char buf[20];
  sprintf(buf,"%20.15e",d);
  return buf;
}


const char* fix3(double d)
{
  static char buf[20];
  sprintf(buf,"%10.3f",d);
  return buf;
}

const char* fix6(double d)
{
  static char buf[20];
  sprintf(buf,"%13.6f",d);
  return buf;
}


const char* fix15(double d)
{
  static char buf[20];
  sprintf(buf,"%20.15f",d);
  return buf;
}
