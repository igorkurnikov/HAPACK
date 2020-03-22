# Microsoft Developer Studio Project File - Name="ipack" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=ipack - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "ipack.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "ipack.mak" CFG="ipack - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "ipack - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "ipack - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "ipack - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "c:\testlib"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /include:"Release/" /nologo /warn:nofileopt
# ADD F90 /compile_only /include:"Release/" /nologo /threads /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /I "..\generic" /I "..\ipack\\" /I "..\DO_LIB\include" /I "..\DO_LIB\template" /I "d:\mpich.nt\sdk\include" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /D "_AFXDLL" /FD /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "ipack - Win32 Debug"

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
# ADD BASE F90 /check:bounds /compile_only /debug:full /include:"Debug/" /nologo /warn:argument_checking /warn:nofileopt
# ADD F90 /browser /check:bounds /compile_only /debug:full /include:"Debug/" /nologo /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /I "..\generic" /I "..\ipack" /I "..\DO_LIB\include" /I "..\DO_LIB\template" /I "d:\mpich.nt\sdk\include" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /FR /YX /FD /GZ /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo /o"c:\testlib\ipack_deb.bsc"
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"c:\testlib\ipack_deb.lib"

!ENDIF 

# Begin Target

# Name "ipack - Win32 Release"
# Name "ipack - Win32 Debug"
# Begin Group "examples"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\integ.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\ipack_main.cpp
# PROP Exclude_From_Build 1
# End Source File
# End Group
# Begin Source File

SOURCE=.\alloc.cpp
# End Source File
# Begin Source File

SOURCE=.\angular.cpp
# End Source File
# Begin Source File

SOURCE=.\angular.h
# End Source File
# Begin Source File

SOURCE=..\generic\avl.h
# End Source File
# Begin Source File

SOURCE=.\basis.cpp
# End Source File
# Begin Source File

SOURCE=.\basis.h
# End Source File
# Begin Source File

SOURCE=.\buffer.h
# End Source File
# Begin Source File

SOURCE=.\coef_array.cpp
# End Source File
# Begin Source File

SOURCE=.\coef_array.h
# End Source File
# Begin Source File

SOURCE=.\coef_set.cpp
# End Source File
# Begin Source File

SOURCE=.\coef_set.h
# End Source File
# Begin Source File

SOURCE=.\const.cpp
# End Source File
# Begin Source File

SOURCE=.\const.h
# End Source File
# Begin Source File

SOURCE=.\contract.cpp
# End Source File
# Begin Source File

SOURCE=..\generic\corefile.h
# End Source File
# Begin Source File

SOURCE=.\e_field.cpp
# End Source File
# Begin Source File

SOURCE=.\e_field.h
# End Source File
# Begin Source File

SOURCE=.\element.cpp
# End Source File
# Begin Source File

SOURCE=.\element.h
# End Source File
# Begin Source File

SOURCE=.\f.cpp
# End Source File
# Begin Source File

SOURCE=.\f.h
# End Source File
# Begin Source File

SOURCE=.\fix.cpp
# End Source File
# Begin Source File

SOURCE=..\generic\flexfile.h
# End Source File
# Begin Source File

SOURCE=.\four_stor.cpp
# End Source File
# Begin Source File

SOURCE=.\four_stor.h
# End Source File
# Begin Source File

SOURCE=.\function.cpp
# End Source File
# Begin Source File

SOURCE=.\function.h
# End Source File
# Begin Source File

SOURCE=..\generic\havl.h
# End Source File
# Begin Source File

SOURCE=.\hrr1.cpp
# End Source File
# Begin Source File

SOURCE=.\hrr1.h
# End Source File
# Begin Source File

SOURCE=.\hrr2.cpp
# End Source File
# Begin Source File

SOURCE=.\hrr2.h
# End Source File
# Begin Source File

SOURCE=.\hrr_set.cpp
# End Source File
# Begin Source File

SOURCE=.\hrr_set.h
# End Source File
# Begin Source File

SOURCE=.\input.cpp
# End Source File
# Begin Source File

SOURCE=.\integ_array.cpp
# End Source File
# Begin Source File

SOURCE=.\integ_array.h
# End Source File
# Begin Source File

SOURCE=..\generic\integ_file.cpp
# End Source File
# Begin Source File

SOURCE=..\generic\integ_file.h
# End Source File
# Begin Source File

SOURCE=.\ipack.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\generic\iphys_file.cpp
# End Source File
# Begin Source File

SOURCE=..\generic\iphys_file.h
# End Source File
# Begin Source File

SOURCE=.\itest.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\kinetic.cpp
# End Source File
# Begin Source File

SOURCE=.\kinetic.h
# End Source File
# Begin Source File

SOURCE=..\generic\memchk.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\generic\memchk.h
# End Source File
# Begin Source File

SOURCE=.\memmgr.cpp

!IF  "$(CFG)" == "ipack - Win32 Release"

!ELSEIF  "$(CFG)" == "ipack - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\memmgr.h
# End Source File
# Begin Source File

SOURCE=..\generic\misc.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\generic\misc.h
# End Source File
# Begin Source File

SOURCE=.\moment.cpp
# End Source File
# Begin Source File

SOURCE=.\moment.h
# End Source File
# Begin Source File

SOURCE=.\nuc_rep.cpp
# End Source File
# Begin Source File

SOURCE=.\nuc_rep.h
# End Source File
# Begin Source File

SOURCE=.\nuclear.cpp
# End Source File
# Begin Source File

SOURCE=.\nuclear.h
# End Source File
# Begin Source File

SOURCE=.\one_el.cpp
# End Source File
# Begin Source File

SOURCE=.\one_el.h
# End Source File
# Begin Source File

SOURCE=.\operators.cpp
# End Source File
# Begin Source File

SOURCE=.\operators.h
# End Source File
# Begin Source File

SOURCE=..\generic\orb_init.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\orbinfo.cpp
# End Source File
# Begin Source File

SOURCE=.\orbinfo.h
# End Source File
# Begin Source File

SOURCE=..\generic\orbital.cpp
# End Source File
# Begin Source File

SOURCE=..\generic\orbital.h
# End Source File
# Begin Source File

SOURCE=.\output.cpp
# End Source File
# Begin Source File

SOURCE=.\output.h
# End Source File
# Begin Source File

SOURCE=.\output1.cpp
# End Source File
# Begin Source File

SOURCE=.\output2.cpp
# End Source File
# Begin Source File

SOURCE=.\output3.cpp
# End Source File
# Begin Source File

SOURCE=.\overlap.cpp
# End Source File
# Begin Source File

SOURCE=.\overlap.h
# End Source File
# Begin Source File

SOURCE=..\generic\parallel.cpp
# End Source File
# Begin Source File

SOURCE=..\generic\parallel.h
# End Source File
# Begin Source File

SOURCE=.\parameters.h
# End Source File
# Begin Source File

SOURCE=.\print.cpp
# End Source File
# Begin Source File

SOURCE=..\generic\qc_utilities.cpp
# End Source File
# Begin Source File

SOURCE=..\generic\qc_utilities.h
# End Source File
# Begin Source File

SOURCE=.\quad.cpp
# End Source File
# Begin Source File

SOURCE=.\quad.h
# End Source File
# Begin Source File

SOURCE=.\rec.cpp
# End Source File
# Begin Source File

SOURCE=.\rec.h
# End Source File
# Begin Source File

SOURCE=.\rec_stor.cpp
# End Source File
# Begin Source File

SOURCE=.\rec_stor.h
# End Source File
# Begin Source File

SOURCE=.\sphere.cpp
# End Source File
# Begin Source File

SOURCE=.\sphere.h
# End Source File
# Begin Source File

SOURCE=.\spin_orb.cpp
# End Source File
# Begin Source File

SOURCE=.\spin_orb.h
# End Source File
# Begin Source File

SOURCE=.\stest.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\stest1.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\storage.cpp
# End Source File
# Begin Source File

SOURCE=.\storage.h
# End Source File
# Begin Source File

SOURCE=.\symdesig.cpp
# End Source File
# Begin Source File

SOURCE=.\symdesig.h
# End Source File
# Begin Source File

SOURCE=.\symlogic.cpp
# End Source File
# Begin Source File

SOURCE=..\generic\symmetry.cpp
# End Source File
# Begin Source File

SOURCE=..\generic\symmetry.h
# End Source File
# Begin Source File

SOURCE=.\symop.cpp
# End Source File
# Begin Source File

SOURCE=.\symop.h
# End Source File
# Begin Source File

SOURCE=.\tpair.h
# End Source File
# Begin Source File

SOURCE=.\tquadruple.h
# End Source File
# Begin Source File

SOURCE=.\transpose.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\twoel.cpp
# End Source File
# Begin Source File

SOURCE=.\twoel.h
# End Source File
# Begin Source File

SOURCE=.\twoel0.cpp
# End Source File
# Begin Source File

SOURCE=.\types.cpp
# End Source File
# Begin Source File

SOURCE=.\typesip.h
# End Source File
# Begin Source File

SOURCE=..\generic\vault.cpp
# End Source File
# Begin Source File

SOURCE=..\generic\vault.h
# End Source File
# Begin Source File

SOURCE=.\vrr.cpp
# End Source File
# Begin Source File

SOURCE=.\vrr.h
# End Source File
# Begin Source File

SOURCE=..\generic\vtype.h
# End Source File
# End Target
# End Project
