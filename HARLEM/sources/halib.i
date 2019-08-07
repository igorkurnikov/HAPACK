%module(package="harlemll") halib

%begin %{
#include <cstddef>  
%}

%include typemaps.i
%include cpointer.i
%include stl.i
%include file.i

%{
#include "haobject.h"
#ifndef HALIB_BASE
//#include "wxgraphs.h"
#endif  
#include "haconsts.h"
#include "tinyxml.h"
#include "haxml.h"
%}

%pointer_class(int,intp);
%pointer_class(double,doublep);

%include "haobject.h"
#ifndef HALIB_BASE
//%include "wxgraphs.h"
#endif
%include "haconsts.h"
%include "tinyxml.h"
%include "haxml.h"