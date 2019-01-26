REM Copy DLLs and other files nessesary to run HARLEM
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
REM Copy python
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
    xcopy /d %PYTHON_DLLS_PATH%\*_d.pyd %OutputDir%\DLLs
    xcopy /d %PYTHON_DLLS_PATH%\*_d.pdb %OutputDir%\DLLs
    
    xcopy /s /e /h /d %PYTHON_HOME_PATH%\Lib %OutputDir%\Lib
    
    xcopy /d %PYTHON_BIN_PATH%\python%PYTHON_MAJOR_VERSION%?_d.dll %OutputDir%
    xcopy /d %PYTHON_BIN_PATH%\python_d.exe %OutputDir%
) else (
    echo "Copying Release Version of Python"
    xcopy /d %PYTHON_DLLS_PATH%\*.pyd %OutputDir%\DLLs
    xcopy /d %PYTHON_DLLS_PATH%\*.pdb %OutputDir%\DLLs
    
    xcopy /s /e /h /d %PYTHON_HOME_PATH%\Lib %OutputDir%\Lib
    
    xcopy /d %PYTHON_BIN_PATH%\python%PYTHON_MAJOR_VERSION%?.dll %OutputDir%
    xcopy /d %PYTHON_BIN_PATH%\python.exe %OutputDir%
)

REM ###########################################################################
REM Copy Things from VCPKG

REM BOOST
if "%IS_DEBUG%" == "Y" (
    set BOOST_SUFFIX=-vc141-mt-gd-x32-1_68.dll
) else (
    set BOOST_SUFFIX=-vc141-mt-x32-1_68.dll
)
set BOOST_LIBS=^
boost_atomic      boost_graph      boost_math_tr1          boost_stacktrace_windbg_cached ^
boost_chrono      boost_locale     boost_prg_exec_monitor  boost_stacktrace_windbg ^
boost_container   boost_log_setup  boost_program_options   boost_system ^
boost_context     boost_log        boost_python36          boost_thread ^
boost_contract    boost_math_c99f  boost_random            boost_timer ^
boost_coroutine   boost_math_c99l  boost_regex             boost_type_erasure ^
boost_date_time   boost_math_c99   boost_serialization     boost_unit_test_framework ^
boost_fiber       boost_math_tr1f  boost_signals           boost_wave ^
boost_filesystem  boost_math_tr1l  boost_stacktrace_noop   boost_wserialization

FOR %%G IN (%BOOST_LIBS%) DO (
    xcopy /d %VCPKG_DLL_PATH%\%%G%BOOST_SUFFIX% %OutputDir%
)

REM WXWIDGETS
if "%IS_DEBUG%" == "Y" (
    set WXVER=311ud_
) else (
    set WXVER=311u_
)
set WX_LIBS=^
wxbase%WXVER%net_vc_custom.dll  wxmsw%WXVER%gl_vc_custom.dll        wxmsw%WXVER%richtext_vc_custom.dll ^
wxbase%WXVER%vc_custom.dll      wxmsw%WXVER%html_vc_custom.dll      wxmsw%WXVER%stc_vc_custom.dll ^
wxbase%WXVER%xml_vc_custom.dll  wxmsw%WXVER%media_vc_custom.dll     wxmsw%WXVER%webview_vc_custom.dll ^
wxmsw%WXVER%adv_vc_custom.dll   wxmsw%WXVER%propgrid_vc_custom.dll  wxmsw%WXVER%xrc_vc_custom.dll ^
wxmsw%WXVER%aui_vc_custom.dll   wxmsw%WXVER%qa_vc_custom.dll ^
wxmsw%WXVER%core_vc_custom.dll  wxmsw%WXVER%ribbon_vc_custom.dll ^
plplot.dll  plplotcxx.dll  plplotwxwidgets.dll

FOR %%G IN (%WX_LIBS%) DO (
    xcopy /d %VCPKG_DLL_PATH%\%%G %OutputDir%
)

REM OTHERS
if "%IS_DEBUG%" == "Y" (
    set OTHER_LIBS=csirocsa.dll mpir.dll jpeg62.dll freetyped.dll libbz2d.dll libpng16d.dll lzma.dll qsastime.dll tiffd.dll zlibd1.dll
) else (
    set OTHER_LIBS=csirocsa.dll mpir.dll jpeg62.dll freetyped.dll libbz2.dll libpng16.dll lzma.dll qsastime.dll tiff.dll zlib1.dll
)

FOR %%G IN (%OTHER_LIBS%) DO (
    xcopy /d %VCPKG_DLL_PATH%\%%G %OutputDir%
)

REM ###########################################################################
REM Copy MKL
set MKL_LIBS=mkl_sequential.dll mkl_core.dll

FOR %%G IN (%MKL_LIBS%) DO (
    xcopy /d %MKL_DLL_PATH%\%%G %OutputDir%
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
