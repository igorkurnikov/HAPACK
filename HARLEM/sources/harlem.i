%module harlemc

%begin %{
#include <mpi.h>  
#include <cstddef>
#include <Python.h>

%}
         
%{
#include <wx/wx.h>
#include "dialogs_wx_2.h"  
void StartHaMainFrameWX();  
void StartHarlemApp();
%}

%include "dialogs_wx_2.h"  
void StartHaMainFrameWX(); 
void StartHarlemApp();
