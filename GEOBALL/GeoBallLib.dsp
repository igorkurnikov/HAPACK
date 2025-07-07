# Microsoft Developer Studio Project File - Name="GeoBallLib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=GeoBallLib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "GeoBallLib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "GeoBallLib.mak" CFG="GeoBallLib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "GeoBallLib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "GeoBallLib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "GeoBallLib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "GeoBallLib___Win32_Release"
# PROP BASE Intermediate_Dir "GeoBallLib___Win32_Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "GeoBallLib___Win32_Release"
# PROP Intermediate_Dir "GeoBallLib___Win32_Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /include:"GeoBallLib___Win32_Release/" /nologo /warn:nofileopt
# ADD F90 /assume:underscore /compile_only /iface:nomixed_str_len_arg /iface:cref /include:"GeoBallLib___Win32_Release/" /names:lowercase /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"c:\testlib\GeoBallLib.lib"

!ELSEIF  "$(CFG)" == "GeoBallLib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "GeoBallLib___Win32_Debug"
# PROP BASE Intermediate_Dir "GeoBallLib___Win32_Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "G:\testlib"
# PROP Intermediate_Dir "GeoBallLib___Win32_Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /include:"GeoBallLib___Win32_Debug/" /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /assume:underscore /check:bounds /compile_only /debug:full /iface:nomixed_str_len_arg /iface:cref /include:"GeoBallLib___Win32_Debug/" /names:lowercase /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"c:\testlib\GeoBallLib_deb.lib"

!ENDIF 

# Begin Target

# Name "GeoBallLib - Win32 Release"
# Name "GeoBallLib - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\adjust.f
# End Source File
# Begin Source File

SOURCE=.\alf_gmp.c
# End Source File
# Begin Source File

SOURCE=.\alfcx.f
# End Source File
# Begin Source File

SOURCE=.\compute_vol.f
# End Source File
# Begin Source File

SOURCE=.\delcx.f
# End Source File
# Begin Source File

SOURCE=.\GeoBall.c
# End Source File
# Begin Source File

SOURCE=.\prepare_deriv.f
# End Source File
# Begin Source File

SOURCE=.\sos_minor_gmp.c
# End Source File
# Begin Source File

SOURCE=.\spacefill_vol.f
# End Source File
# Begin Source File

SOURCE=.\truncate_real.f
# End Source File
# Begin Source File

SOURCE=.\vector.f
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\defines.h
# End Source File
# Begin Source File

SOURCE=.\geoball.h
# End Source File
# Begin Source File

SOURCE=.\gmp.h
# End Source File
# Begin Source File

SOURCE=.\gmpvar.h
# End Source File
# End Group
# End Target
# End Project
