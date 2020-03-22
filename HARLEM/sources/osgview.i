%module osgview
%import hamolset.i
%import halib.i
%include typemaps.i
%include cpointer.i
%include stl.i
%include file.i
%{
#include "hasurface.h"
#include "mVueFrame.h"
#include "mVueGlCanvas.h"
%}

%pointer_class(int,intp);
%pointer_class(double,doublep);


%include "mVueFrame.h"
%include "mVueGlCanvas.h"

