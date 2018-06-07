# Microsoft Developer Studio Project File - Name="my_vflib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=my_vflib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "my_vflib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "my_vflib.mak" CFG="my_vflib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "my_vflib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "my_vflib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "my_vflib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "c:\testlib"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /include:"Release/" /nologo /warn:nofileopt
# ADD F90 /compile_only /include:"Release/" /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"c:\testlib\vflib.lib"

!ELSEIF  "$(CFG)" == "my_vflib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "c:\testlib"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /include:"Debug/" /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /include:"Debug/" /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"c:\testlib\vflib_deb.lib"

!ENDIF 

# Begin Target

# Name "my_vflib - Win32 Release"
# Name "my_vflib - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=src\argedit.cpp
# End Source File
# Begin Source File

SOURCE=src\argloader.cpp
# End Source File
# Begin Source File

SOURCE=src\argraph.cpp
# End Source File
# Begin Source File

SOURCE=src\error.cpp
# End Source File
# Begin Source File

SOURCE=src\gene.cpp
# End Source File
# Begin Source File

SOURCE=src\gene_mesh.cpp
# End Source File
# Begin Source File

SOURCE=src\match.cpp
# End Source File
# Begin Source File

SOURCE=src\sd_state.cpp
# End Source File
# Begin Source File

SOURCE=src\sortnodes.cpp
# End Source File
# Begin Source File

SOURCE=src\ull_state.cpp
# End Source File
# Begin Source File

SOURCE=src\ull_sub_state.cpp
# End Source File
# Begin Source File

SOURCE=src\vf2_mono_state.cpp
# End Source File
# Begin Source File

SOURCE=src\vf2_state.cpp
# End Source File
# Begin Source File

SOURCE=src\vf2_sub_state.cpp
# End Source File
# Begin Source File

SOURCE=src\vf_mono_state.cpp
# End Source File
# Begin Source File

SOURCE=src\vf_state.cpp
# End Source File
# Begin Source File

SOURCE=src\vf_sub_state.cpp
# End Source File
# Begin Source File

SOURCE=src\xsubgraph.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=src\allocpool.h
# End Source File
# Begin Source File

SOURCE=..\..\MY_VFLIB\src\argedit.h
# End Source File
# Begin Source File

SOURCE=..\..\MY_VFLIB\src\argloader.h
# End Source File
# Begin Source File

SOURCE=..\..\MY_VFLIB\src\argraph.h
# End Source File
# Begin Source File

SOURCE=src\dict.h
# End Source File
# Begin Source File

SOURCE=src\error.h
# End Source File
# Begin Source File

SOURCE=src\gene.h
# End Source File
# Begin Source File

SOURCE=src\gene_mesh.h
# End Source File
# Begin Source File

SOURCE=src\match.h
# End Source File
# Begin Source File

SOURCE=src\sd_state.h
# End Source File
# Begin Source File

SOURCE=src\sortnodes.h
# End Source File
# Begin Source File

SOURCE=src\state.h
# End Source File
# Begin Source File

SOURCE=src\ull_state.h
# End Source File
# Begin Source File

SOURCE=src\ull_sub_state.h
# End Source File
# Begin Source File

SOURCE=src\vf2_mono_state.h
# End Source File
# Begin Source File

SOURCE=src\vf2_state.h
# End Source File
# Begin Source File

SOURCE=src\vf2_sub_state.h
# End Source File
# Begin Source File

SOURCE=src\vf_mono_state.h
# End Source File
# Begin Source File

SOURCE=src\vf_state.h
# End Source File
# Begin Source File

SOURCE=src\vf_sub_state.h
# End Source File
# Begin Source File

SOURCE=src\xsubgraph.h
# End Source File
# End Group
# End Target
# End Project
