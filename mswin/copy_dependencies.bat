REM Copy DLLs and other files necessary to run HARLEM
SET script_path=%~dp0

echo "Copying HARLEM dependencies ..."
if not defined CONF (
    echo "Variable CONF is not defined"
    echo "This script should run as post-build event in VS"
    exit 1
)
if not defined VCPKG_DLL_PATH (
    echo "Variable VCPKG_DLL_PATH is not defined"
    echo "This script should run as post-build event in VS"
    exit 1
)
if not defined PYTHON_DLLS_PATH (
    echo "Variable PYTHON_DLLS_PATH is not defined"
    echo "This script should run as post-build event in VS"
    exit 1
)
if not defined MKL_DLL_PATH (
    echo "Variable MKL_DLL_PATH is not defined"
    echo "This script should run as post-build event in VS"
    exit 1
)
if not defined WX_DLLS_PATH (
    echo "Variable WX_DLLS_PATH is not defined"
    echo "This script should run as post-build event in VS"
    exit 1
)

echo "Configuration: %CONF%"
echo "Script Path: %script_path%"
SET OutputDir="%script_path%%CONF%"
echo "Output Dir: %OutputDir%"

REM Get Python Major Version
if not x%CONF:PY3=%==x%CONF% (set PYTHON_MAJOR_VERSION=3) else (set PYTHON_MAJOR_VERSION=2)
echo "Python Major Version: %PYTHON_MAJOR_VERSION%"

REM Is it debug
if not x%CONF:Debug=%==x%CONF% (set "IS_DEBUG=Y") else (set "IS_DEBUG=N")
echo "Is it debug: %IS_DEBUG%"

REM ###########################################################################
REM Make directories if not exists
if not exist "%OutputDir%\molset\NUL" (
    mkdir "%OutputDir%\molset"
) else (
    echo "%OutputDir%\molset already exists"
)
if not exist "%OutputDir%\pnpsll\NUL" (
    mkdir "%OutputDir%\pnpsll"
) else (
    echo "%OutputDir%\pnpsll already exists"
)

REM ###########################################################################
REM Copy python
echo "Copying PYTHON"
if not exist "%OutputDir%\DLLs\NUL" (
    mkdir "%OutputDir%\DLLs"
) else (
    echo "%OutputDir%\DLLs already exists"
)
if not exist "%OutputDir%\Lib\NUL" (
    mkdir "%OutputDir%\Lib"
) else (
    echo "%OutputDir%\Lib already exists"
)


if "%IS_DEBUG%" == "Y" (
    echo "Copying Debug Version of Python"
    xcopy /y /d %PYTHON_DLLS_PATH%\*_d.pyd %OutputDir%\DLLs
    xcopy /y /d %PYTHON_DLLS_PATH%\*_d.pdb %OutputDir%\DLLs
    
    xcopy /y /s /e /h /d %PYTHON_HOME_PATH%\Lib %OutputDir%\Lib
    
    xcopy /y /d %PYTHON_BIN_PATH%\python%PYTHON_MAJOR_VERSION%?_d.dll %OutputDir%
    xcopy /y /d %PYTHON_BIN_PATH%\python_d.exe %OutputDir%
) else (
    echo "Copying Release Version of Python"
    xcopy /y /d %PYTHON_DLLS_PATH%\*.pyd %OutputDir%\DLLs
    xcopy /y /d %PYTHON_DLLS_PATH%\*.pdb %OutputDir%\DLLs
    
    xcopy /y  /s /e /h /d %PYTHON_HOME_PATH%\Lib %OutputDir%\Lib
    
    xcopy /y  /d %PYTHON_BIN_PATH%\python%PYTHON_MAJOR_VERSION%?.dll %OutputDir%
    xcopy /y /d %PYTHON_BIN_PATH%\python.exe %OutputDir%
)
xcopy /y /d %PYTHON_DLLS_PATH%\*.dll %OutputDir%\DLLs
xcopy /y /d %PYTHON_DLLS_PATH%\*.ico %OutputDir%\DLLs

REM ###########################################################################
REM Copy Things from VCPKG

REM BOOST
echo "Copying BOOST Dlls"
if "%IS_DEBUG%" == "Y" (
    set BOOST_SUFFIX=-vc142-mt-gd-x*-1_70.dll
) else (
    set BOOST_SUFFIX=-vc142-mt-x*-1_70.dll
)
set BOOST_LIBS=^
boost_atomic      boost_graph      boost_math_tr1          boost_stacktrace_windbg_cached ^
boost_chrono      boost_locale     boost_prg_exec_monitor  boost_stacktrace_windbg ^
boost_container   boost_log_setup  boost_program_options   boost_system ^
boost_context     boost_log        boost_python37          boost_thread ^
boost_contract    boost_math_c99f  boost_random            boost_timer ^
boost_coroutine   boost_math_c99l  boost_regex             boost_type_erasure ^
boost_date_time   boost_math_c99   boost_serialization     boost_unit_test_framework ^
boost_fiber       boost_math_tr1f  boost_wserialization    boost_wave ^
boost_filesystem  boost_math_tr1l  boost_stacktrace_noop 

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
xcopy /y /d %WX_DLLS_PATH%\wxbase*%WXVER%_*.dll  %OutputDir%\molset
xcopy /y /d %WX_DLLS_PATH%\wxmsw*%WXVER%_*.dll  %OutputDir%\molset

REM ###########################################################################
REM Copy PLPLOT 
echo "Copying PLPLOT Dlls"
if "%IS_DEBUG%" == "Y" (
    if exist %VCPKG_DLL_PATH%\plplotd.dll (
        set PLPLOT_LIB=plplotd.dll  plplotcxxd.dll  plplotwxwidgetsd.dll csirocsad.dll qsastimed.dll 
    ) else (
        set PLPLOT_LIB=plplot.dll  plplotcxx.dll  plplotwxwidgets.dll csirocsa.dll qsastime.dll 
    )
) else (
    set PLPLOT_LIB=plplot.dll  plplotcxx.dll  plplotwxwidgets.dll csirocsa.dll qsastime.dll 
)
FOR %%G IN (%PLPLOT_LIB%) DO (
    xcopy /y /d %VCPKG_DLL_PATH%\%%G %OutputDir%\molset
)
xcopy /y /d %VCPKG_DLL_PATH%\wxbase*.dll  %OutputDir%\molset
xcopy /y /d %VCPKG_DLL_PATH%\wxmsw*.dll  %OutputDir%\molset

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
REM Copy MKL
echo "Copying MKL Dlls"
set MKL_LIBS=mkl_sequential.dll mkl_core.dll

FOR %%G IN (%MKL_LIBS%) DO (
    xcopy /y /d %MKL_DLL_PATH%\%%G %OutputDir%\molset
)
REM ###########################################################################
REM Copy MPI

REM ###########################################################################
REM Copy DB
if not exist "%OutputDir%\residues_db\NUL" (
    mkdir "%OutputDir%\residues_db"
) else (
    echo "%OutputDir%\residues_db already exists"
)

xcopy /s /e /h /d %script_path%\..\residues_db %OutputDir%\residues_db

REM ###########################################################################
REM Copy harlempy
echo "Copying HARLEMPY"
if not exist "%OutputDir%\harlempy\NUL" (
    mkdir "%OutputDir%\harlempy"
) else (
    echo "%OutputDir%\harlempy already exists"
)

xcopy /y /d %script_path%\..\HARLEM\harlempy\* %OutputDir%\harlempy\
REM ###########################################################################
REM Copy molset(harlemll) module axxiliary python files
echo "Copying molset module axxiliary python files"
if not exist "%OutputDir%\molset\NUL" (
    mkdir "%OutputDir%\molset"
) else (
    echo "%OutputDir%\molset already exists"
)

xcopy /y /d %script_path%\..\HARLEMLL\*.py %OutputDir%\molset\


REM ###########################################################################
REM Copy wxextra
if not exist "%OutputDir%\wxextra\NUL" (
    mkdir "%OutputDir%\wxextra"
) else (
    echo "%OutputDir%\wxextra already exists"
)

xcopy /y /d %script_path%\..\PNPS\wxextra\* %OutputDir%\wxextra\

