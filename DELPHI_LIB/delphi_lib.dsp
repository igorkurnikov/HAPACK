# Microsoft Developer Studio Project File - Name="delphi_lib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=delphi_lib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "delphi_lib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "delphi_lib.mak" CFG="delphi_lib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "delphi_lib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "delphi_lib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "delphi_lib - Win32 Release"

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
# ADD F90 /assume:underscore /compile_only /extend_source:80 /iface:nomixed_str_len_arg /iface:cref /include:"Release/" /libs:dll /names:lowercase /nologo /threads /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "delphi_lib - Win32 Debug"

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
# ADD F90 /assume:underscore /check:bounds /compile_only /debug:full /extend_source:80 /iface:nomixed_str_len_arg /iface:cref /include:"Debug/" /libs:dll /names:lowercase /nologo /threads /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"c:\testlib\delphi_lib_deb.lib"

!ENDIF 

# Begin Target

# Name "delphi_lib - Win32 Release"
# Name "delphi_lib - Win32 Debug"
# Begin Source File

SOURCE=.\source\chrgit4.f
# End Source File
# Begin Source File

SOURCE=.\source\cputime.f
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\source\cputime_pc.f
# End Source File
# Begin Source File

SOURCE=.\source\itit4j.f
DEP_F90_ITIT4=\
	".\source\qdiffpar4.h"\
	
# End Source File
# Begin Source File

SOURCE=.\source\mem4.f
DEP_F90_MEM4_=\
	".\source\qdiffpar4.h"\
	
# End Source File
# Begin Source File

SOURCE=.\source\mkeps4f.f
DEP_F90_MKEPS=\
	".\source\qdiffpar4.h"\
	
# End Source File
# Begin Source File

SOURCE=.\source\non.f
DEP_F90_NON_F=\
	".\source\qdiffpar4.h"\
	
# End Source File
# Begin Source File

SOURCE=.\source\phintp4.f
DEP_F90_PHINT=\
	".\source\qdiffpar4.h"\
	
# End Source File
# Begin Source File

SOURCE=.\source\qdiff4s.f
DEP_F90_QDIFF=\
	".\source\qdiffpar4.h"\
	
# End Source File
# Begin Source File

SOURCE=.\source\qdiffpar4.h
# End Source File
# Begin Source File

SOURCE=.\source\relfac4b.f
DEP_F90_RELFA=\
	".\source\qdiffpar4.h"\
	
# End Source File
# Begin Source File

SOURCE=.\source\scaler4.f
DEP_F90_SCALE=\
	".\source\qdiffpar4.h"\
	
# End Source File
# Begin Source File

SOURCE=.\source\setbc4.f
DEP_F90_SETBC=\
	".\source\qdiffpar4.h"\
	
# End Source File
# Begin Source File

SOURCE=.\source\setcrg4.f
DEP_F90_SETCR=\
	".\source\qdiffpar4.h"\
	
# End Source File
# Begin Source File

SOURCE=.\source\setin4.f
DEP_F90_SETIN=\
	".\source\qdiffpar4.h"\
	
# End Source File
# Begin Source File

SOURCE=.\source\setout4.f
DEP_F90_SETOU=\
	".\source\qdiffpar4.h"\
	
# End Source File
# End Target
# End Project
