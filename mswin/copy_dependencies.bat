REM Copy and link DLLs and other files necessary to run HARLEM
SET script_path=%~dp0

echo "Copying HARLEM dependencies ..."
if not defined CONF (
    echo "Variable CONF is not defined"
    echo "This script should run as post-build event in VS"
    exit 1
)

REM Is it debug
if not x%CONF:Debug=%==x%CONF% (set "IS_DEBUG=Y") else (set "IS_DEBUG=N")
echo "Is it debug: %IS_DEBUG%"

if not defined VCPKG_DLL_PATH (
    echo "Variable VCPKG_DLL_PATH is not defined"
if "%IS_DEBUG%" == "Y" (
    set VCPKG_DLL_PATH="C:\MYPROG\vcpkg\installed\x64-windows\debug\bin"
) else (
	set VCPKG_DLL_PATH="C:\MYPROG\vcpkg\installed\x64-windows\bin"
)
)

echo VCPKG_DLL_PATH set to %VCPKG_DLL_PATH%

if not defined IFORT_DLL_PATH (
    echo "Variable IFORT_DLL_PATH is not defined"
	set IFORT_DLL_PATH="C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\redist\intel64_win\compiler"
)
echo IFORT_DLL_PATH set to %IFORT_DLL_PATH%

if not defined MKL_DLL_PATH (
    echo "Variable MKL_DLL_PATH is not defined"
	set MKL_DLL_PATH="C:\Program Files (x86)\Intel\oneAPI\mkl\latest\redist\intel64"
)
echo MKL_DLL_PATH set to %MKL_DLL_PATH%
    
if not defined WX_DLLS_PATH (
    echo "Variable WX_DLLS_PATH is not defined"
    echo "This script should run as post-build event in VS"
)
echo WX_DLLS_PATH set to %WX_DLLS_PATH%

echo "Configuration: %CONF%"
echo "Script Path: %script_path%"
SET OutputDir="%script_path%%CONF%"
echo "Output Dir: %OutputDir%"

REM ###########################################################################
REM Make directories if not exists
if not exist "%OutputDir%\molset\NUL" (
    mkdir "%OutputDir%\molset"
) 

if not exist "%OutputDir%\pnpsll\NUL" (
    mkdir "%OutputDir%\pnpsll"
) 

REM ###########################################################################
echo "Linking PYTHON Directories"
if not exist "%OutputDir%\DLLs\NUL" (
    mklink /D %OutputDir%\DLLs %PYTHON_HOME_PATH%\DLLs
)

if not exist "%OutputDir%\Lib\NUL" (
    mklink /D "%OutputDir%\Lib" %PYTHON_HOME_PATH%\Lib
) 


if "%IS_DEBUG%" == "Y" (
    echo "Copying Debug Version of Python"
    
    xcopy /y /d %PYTHON_HOME_PATH%\python3_d.dll %OutputDir%
    xcopy /y /d %PYTHON_HOME_PATH%\python3*_d.dll %OutputDir%
    xcopy /y /d %PYTHON_HOME_PATH%\python_d.exe %OutputDir%
	
) else (
    echo "Copying Release Version of Python"
    
    xcopy /y /d %PYTHON_HOME_PATH%\python3.dll %OutputDir%
	xcopy /y /d %PYTHON_HOME_PATH%\python3*.dll %OutputDir%
    xcopy /y /d %PYTHON_HOME_PATH%\python.exe %OutputDir%
)

REM ###########################################################################
REM Copy Things from VCPKG

REM BOOST
echo "Copying BOOST Dlls"
if "%IS_DEBUG%" == "Y" (
    set BOOST_SUFFIX=-vc142-mt-gd-x*-1_79.dll
) else (
    set BOOST_SUFFIX=-vc142-mt-x*-1_79.dll
)

REM set BOOST_LIBS=^
REM boost_atomic      boost_graph      boost_math_tr1          boost_stacktrace_windbg_cached ^
REM boost_chrono      boost_locale     boost_prg_exec_monitor  boost_stacktrace_windbg ^
REM boost_container   boost_log_setup  boost_program_options   boost_system ^
REM boost_context     boost_log        boost_python37          boost_thread ^
REM boost_contract    boost_math_c99f  boost_random            boost_timer ^
REM boost_coroutine   boost_math_c99l  boost_regex             boost_type_erasure ^
REM boost_date_time   boost_math_c99   boost_serialization     boost_unit_test_framework ^
REM boost_fiber       boost_math_tr1f  boost_wserialization    boost_wave ^
REM boost_filesystem  boost_math_tr1l  boost_stacktrace_noop 
set BOOST_LIBS=boost_filesystem boost_system

FOR %%G IN (%BOOST_LIBS%) DO (
    xcopy /y /d %VCPKG_DLL_PATH%\%%G%BOOST_SUFFIX% %OutputDir%\molset
)

REM WXWIDGETS
echo "Copying wxWidgets Dlls"
if "%IS_DEBUG%" == "Y" (
    set WXVER=ud
) else (
    set WXVER=u
)

xcopy /y /d %WX_DLLS_PATH%\wxbase*%WXVER%_vc*.dll  %OutputDir%\molset
xcopy /y /d %WX_DLLS_PATH%\wxmsw*%WXVER%_*core*.dll  %OutputDir%\molset

REM OTHERS
if "%IS_DEBUG%" == "Y" (
    set OTHER_LIBS=mpir.dll jpeg62.dll freetyped.dll libbz2d.dll libpng16d.dll lzma.dll tiffd.dll zlibd1.dll 
) else (
    set OTHER_LIBS=mpir.dll jpeg62.dll freetype.dll libbz2.dll libpng16.dll lzma.dll tiff.dll zlib1.dll
)

FOR %%G IN (%OTHER_LIBS%) DO (
    xcopy /y /d %VCPKG_DLL_PATH%\%%G %OutputDir%\molset
)

REM OTHERS of PNPS
if "%IS_DEBUG%" == "Y" (
    set OTHER_LIBS=zlibd1.dll
) else (
    set OTHER_LIBS=zlib1.dll
)

FOR %%G IN (%OTHER_LIBS%) DO (
    xcopy /y /d %VCPKG_DLL_PATH%\%%G %OutputDir%\pnpsll
)


REM ###########################################################################
echo "Copying IFORT Dlls"
if "%IS_DEBUG%" == "Y" (
    set IFORT_LIBS=libifcoremdd.dll libmmdd.dll libmmd.dll svml_dispmd.dll
) else (	
	set IFORT_LIBS=libifcoremd.dll libmmd.dll svml_dispmd.dll
)
FOR %%G IN (%IFORT_LIBS%) DO (
    xcopy /y /d %IFORT_DLL_PATH%\%%G %OutputDir%\molset
)
REM ###########################################################################
echo "Copying MKL Dlls"
REM set MKL_LIBS=mkl_sequential.dll mkl_core.dll mkl_avx2.dll
set MKL_LIBS=mkl_sequential.2.dll mkl_core.2.dll mkl_avx2.2.dll

FOR %%G IN (%MKL_LIBS%) DO (
    xcopy /y /d %MKL_DLL_PATH%\%%G %OutputDir%\molset
)

REM ###########################################################################
REM Copy MPI

xcopy /y /d C:\Windows\System32\msmpi.dll %OutputDir%\molset

REM ###########################################################################
REM Link DB
if not exist "%OutputDir%\residues_db\NUL" (
    mklink /D %OutputDir%\residues_db %script_path%\..\residues_db 
) 

REM ###########################################################################
echo "Linking molset module python files"

if not exist %OutputDir%\molset\__init__.py (
    mklink %OutputDir%\molset\__init__.py  %script_path%\..\HARLEM\molset\__init__.py
)
if not exist %OutputDir%\molset\harlempy (
    mklink /D %OutputDir%\molset\harlempy     %script_path%\..\HARLEM\molset\harlempy
)
if not exist %OutputDir%\molset\molset_ext (
    mklink /D %OutputDir%\molset\molset_ext   %script_path%\..\HARLEM\molset\molset_ext
)


REM ###########################################################################
REM Link wxextra
if not exist "%OutputDir%\wxextra" (
    mklink /D %OutputDir%\wxextra\ %script_path%\..\PNPS\wxextra
) 



