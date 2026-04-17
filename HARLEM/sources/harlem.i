%module harlemc

%begin %{
//#include <mpi.h>
#include <cstddef>
#ifdef _WIN32
#pragma push_macro("_DEBUG")
#undef _DEBUG
#include <Python.h>
#pragma pop_macro("_DEBUG")
#else
#include <Python.h>
#endif

%}
         
%{
#include <wx/wx.h>
#include "dialogs_wx_2.h"  
void StartHaMainFrameWX();  
%}

%include "dialogs_wx_2.h" 
void StartHaMainFrameWX(); 
