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

if not defined PYTHON_HOME_PATH (
    echo "Variable PYTHON_HOME_PATH is not defined, defaulting to D:\Python312"
    set PYTHON_HOME_PATH="D:\Python312"
)
echo PYTHON_HOME_PATH set to %PYTHON_HOME_PATH%

echo "Configuration: %CONF%"
echo "Script Path: %script_path%"
SET OutputDir="%script_path%%CONF%"
echo "Output Dir: %OutputDir%"

REM WXWIDGETS
echo "Copying wxWidgets Dlls"
REM Always use RELEASE wxWidgets — HARLEM links release wx even in Debug
REM to stay consistent with wxPython's release wx DLLs (wxTheApp shared instance).
set WXVER=u

xcopy /y /d %WX_DLLS_PATH%\wxbase*%WXVER%_vc*.dll  %OutputDir%\molset
xcopy /y /d %WX_DLLS_PATH%\wxmsw*%WXVER%_*core*.dll  %OutputDir%\molset

REM OTHERS
if "%IS_DEBUG%" == "Y" (
    set OTHER_LIBS=mpir.dll jpeg62.dll freetyped.dll libbz2d.dll libpng16d.dll lzma.dll tiffd.dll zlibd1.dll 
) else (
REM    set OTHER_LIBS=mpir.dll jpeg62.dll freetype.dll libbz2.dll libpng16.dll lzma.dll tiff.dll zlib1.dll
    set OTHER_LIBS=mpir.dll jpeg62.dll tiff.dll zlib1.dll openblas.dll  
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

REM ###########################################################################
REM Copy MPI

xcopy /y /d C:\Windows\System32\msmpi.dll %OutputDir%\molset

REM ###########################################################################
REM Copy Python stdlib (Lib/ and DLLs/) from PYTHON_HOME_PATH.
REM /E = include subdirs (and empty), /I = dest is directory, /Y = no overwrite prompt, /D = copy only newer.
echo "Copying Python stdlib from %PYTHON_HOME_PATH%"
xcopy /E /I /Y /D %PYTHON_HOME_PATH%\Lib  %OutputDir%\Lib
xcopy /E /I /Y /D %PYTHON_HOME_PATH%\DLLs %OutputDir%\DLLs

REM Copy Python interpreter and core DLLs (release variants only — debug Python is not used).
echo "Copying Python interpreter and core DLLs"
xcopy /Y /D %PYTHON_HOME_PATH%\python.exe   %OutputDir%
xcopy /Y /D %PYTHON_HOME_PATH%\pythonw.exe  %OutputDir%
xcopy /Y /D %PYTHON_HOME_PATH%\python312.dll %OutputDir%
xcopy /Y /D %PYTHON_HOME_PATH%\python3.dll  %OutputDir%

REM ###########################################################################
echo "Linking harlem gui python files"

if not exist %OutputDir%\molset\harlempy (
    mklink /D %OutputDir%\molset\harlempy     %script_path%\..\MOLSET\molset\harlempy
)

REM ###########################################################################
REM Link wxextra
if not exist "%OutputDir%\wxextra" (
    mklink /D %OutputDir%\wxextra\ %script_path%\..\PNPS\wxextra
) 

echo "Done library copying and linking"



