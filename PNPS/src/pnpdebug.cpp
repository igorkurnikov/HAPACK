#include "pnpdebug.h"
#include <string>
#include <stdexcept>
#include <stdio.h>

void pnpThrowErrorRaw(const int line,const char* filename,const char* str,...)
{
	int i;
	char *s=new char[10240];

	i=sprintf(s,"(%s:%d): ",filename,line);

	va_list arg_list;
	va_start(arg_list,str);
	i+=vsprintf (s+i,str, arg_list);
	va_end(arg_list);

	std::runtime_error e(s);
	delete s;

	throw e;
}