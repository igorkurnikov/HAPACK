#include <iostream.h>
#include <stringc.h>

/*TEX
This function expects that the stream is in the beginning ofthe line 
on entry. It will read a line and scan for the string "#include"
at the very beginning. If this string is found, it returns \name{true}
and the \name{argument} variable contains the argument of the incluide 
directive. The function returns QLOCAL for files included with regular quotes,
QSYS for files with system quotes.
*/

enum ParseCode { QNONE,QLOCAL,QSYS};

ParseCode parse_line(istream& is,String& argument)
{
  char c;
  char token[] = "#include";
  int  match = 0;
  int  tlen  = strlen(token);

  const int maxname = 256; // this must be enough .....
  char name[256];  
  int  ncount = 0;
  int nomatch = 0;

  while(1)
  {
    is >> c;
    if (!is)
      error("Read Error in Parse-Line");

    if (c == '\n')
    {
      if (code != 0)
	error("Parse-Line: found include, but no argument. ");
      return QNONE;
    }

    if (nomatch)
      continue;

    switch(match)
    {
    case tlen: // token found search for " or <
      if (c == '"')
      {
	code = QLOCAL;
	match++;
      }
      if (c == '<')
      {
	code = QSYS; 
	match++;
      }
      break;
    case tlen+1: // find end of argument
      if ((code == QLOCAL && c == '"') ||
	  (code == QSYS   && c == '>'))
      {
	name[ncount] = '\0';
	argument = name;
	return code;
      }
      name[ncount++] = c;
      if (ncount >= maxname)
	error("String Length excceded in parse-line. ");
      break;
    default:  // here we ned when looking for the token
      if (c != token[match])
      {
	match++;
	break;
      }
      // didnt match, continue reading until EOL then quit.
      nomatch = true;
    }
  }
}





