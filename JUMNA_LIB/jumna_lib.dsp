# Microsoft Developer Studio Project File - Name="jumna_lib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=jumna_lib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "jumna_lib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "jumna_lib.mak" CFG="jumna_lib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "jumna_lib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "jumna_lib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "jumna_lib - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Yu"stdafx.h" /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /FD /c
# SUBTRACT CPP /YX /Yc /Yu
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "jumna_lib - Win32 Debug"

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
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Yu"stdafx.h" /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /FD /GZ /c
# SUBTRACT CPP /YX /Yc /Yu
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"c:\testlib\jumna_deb_lib.lib"

!ENDIF 

# Begin Target

# Name "jumna_lib - Win32 Release"
# Name "jumna_lib - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\sources\ambout.f
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\sources\ampar91.f
DEP_F90_AMPAR=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\ampar94.f
DEP_F90_AMPAR9=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\ang.f
DEP_F90_ANG_F=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\assemb.f
DEP_F90_ASSEM=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\axeout.f
DEP_F90_AXEOU=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\backbo.f
DEP_F90_BACKB=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\build.f
DEP_F90_BUILD=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\cataxe.f
DEP_F90_CATAX=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\close.f
DEP_F90_CLOSE=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\combi.f
DEP_F90_COMBI=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\cosdir.f
DEP_F90_COSDI=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\deltor.f
DEP_F90_DELTO=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\delval.f
DEP_F90_DELVA=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\dis.f
DEP_F90_DIS_F=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\disth.f
DEP_F90_DISTH=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\dtcabe.f
# End Source File
# Begin Source File

SOURCE=.\sources\ecomp.f
DEP_F90_ECOMP=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\ecomp91.f
DEP_F90_ECOMP9=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\ecomp94.f
DEP_F90_ECOMP94=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\enecyl.f
DEP_F90_ENECY=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\energ91.f
DEP_F90_ENERG=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\energ94.f
DEP_F90_ENERG9=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\energy.f
DEP_F90_ENERGY=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\equim.f
DEP_F90_EQUIM=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\equiv.f
DEP_F90_EQUIV=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\fitsug.f
DEP_F90_FITSU=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\fixadd.c
# End Source File
# Begin Source File

SOURCE=.\sources\flex.f
DEP_F90_FLEX_=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\grads.f
DEP_F90_GRADS=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\gradt.f
DEP_F90_GRADT=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\helix.f
DEP_F90_HELIX=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\helloc.f
DEP_F90_HELLO=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\insout.f
DEP_F90_INSOU=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\joie.f
DEP_F90_JOIE_=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\join.f
DEP_F90_JOIN_=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\jumcall.f
DEP_F90_JUMCA=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\kapa.f
DEP_F90_KAPA_=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\kline.f
DEP_F90_KLINE=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\ligaxe.f
DEP_F90_LIGAX=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\ligput.f
DEP_F90_LIGPU=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\liner.f
DEP_F90_LINER=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\lmoout.f
DEP_F90_LMOOU=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\loops.f
DEP_F90_LOOPS=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\macout.f
DEP_F90_MACOU=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\mc11a.f
DEP_F90_MC11A=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\microb.f
DEP_F90_MICRO=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\minew.f
DEP_F90_MINEW=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\minfor.f
DEP_F90_MINFO=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\minim.f
DEP_F90_MINIM=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\molin.f
DEP_F90_MOLIN=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\move.f
DEP_F90_MOVE_=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\nml.f
DEP_F90_NML_F=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\openb.f
DEP_F90_OPENB=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\pairc.f
DEP_F90_PAIRC=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\pairs.f
DEP_F90_PAIRS=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\pdbout.f
DEP_F90_PDBOU=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\pdbouta.f
DEP_F90_PDBOUT=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\penalty.f
DEP_F90_PENAL=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\present.f
DEP_F90_PRESE=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\putbac.f
DEP_F90_PUTBA=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\renoe.f
DEP_F90_RENOE=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\reord.f
DEP_F90_REORD=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\reset.f
DEP_F90_RESET=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\setbac.f
DEP_F90_SETBA=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\setd.f
DEP_F90_SETD_=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\setgeo.f
DEP_F90_SETGE=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\setvar.f
DEP_F90_SETVA=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\tbasl.f
DEP_F90_TBASL=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\tjunl.f
DEP_F90_TJUNL=\
	".\sources\jumna_data.inc"\
	
# End Source File
# Begin Source File

SOURCE=.\sources\torp.f
DEP_F90_TORP_=\
	".\sources\jumna_data.inc"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# End Group
# Begin Source File

SOURCE=.\Readme.txt
# End Source File
# End Target
# End Project
