%module halib_base_py
//%import hamolset.i

%include typemaps.i
%include cpointer.i
%include stl.i
%include file.i

%{
#include "haobject.h"
#ifndef HALIB_BASE
#include "wxgraphs.h"
#endif
#include "haconsts.h"
#include "tinyxml.h"
%}

%pointer_class(int,intp);
%pointer_class(double,doublep);

%include "haobject.h"
#ifndef HALIB_BASE
%include "wxgraphs.h"
#endif
%include "haconsts.h"
%include "tinyxml.h"