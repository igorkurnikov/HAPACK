rem HARLEM conda build for windows

SET script_path=%~dp0

if "%PY_VER%"=="3.6" (
	set PYTHON_LIBRARY=python36.lib
) else if  "%PY_VER%"=="3.7" (
	set PYTHON_LIBRARY=python37.lib
) else if  "%PY_VER%"=="3.8" (
	set PYTHON_LIBRARY=python38.lib
) else if  "%PY_VER%"=="3.9" (
	set PYTHON_LIBRARY=python39.lib
) else (
	echo "Unexpected version of python"
	exit 1
)

REM ###########################################################################
REM Make directories if not exists
if not exist "%PREFIX%\doc\NUL" (
    mkdir "%PREFIX%\doc"
) else (
    echo "%PREFIX%\doc already exists"
)

if not exist "%PREFIX%\doc\advanced_manual_html\NUL" (
    mkdir "%PREFIX%\doc\advanced_manual_html"
) else (
    echo "%PREFIX%\doc\advanced_manual_html already exists"
)

if not exist "%PREFIX%\doc\HARLEM_BeginnerUserManual_files\NUL" (
    mkdir "%PREFIX%\doc\HARLEM_BeginnerUserManual_files"
) else (
    echo "%PREFIX%\doc\HARLEM_BeginnerUserManual_files already exists"
)

if not exist "%PREFIX%\examples\NUL" (
    mkdir "%PREFIX%\examples"
) else (
    echo "%PREFIX%\examples already exists"
)

if not exist "%PREFIX%\residues_db\NUL" (
    mkdir "%PREFIX%\residues_db"
) else (
    echo "%PREFIX%\residues_db already exists"
)

if not exist "%PREFIX%\molset\NUL" (
    mkdir "%PREFIX%\molset"
) else (
    echo "%PREFIX%\molset already exists"
)

if not exist "%PREFIX%\wx\NUL" (
    mkdir "%PREFIX%\wx"
) else (
    echo "%PREFIX%\wx already exists"
)

REM if not exist "%PREFIX%\Lib\site-packages\wx\NUL" (
REM    mkdir "%PREFIX%\Lib\site-packages\wx"
REM ) else (
REM    echo "%PREFIX%\Lib\site-packages\wx already exists"
REM )

rem cd ../mswin/Release_x64/molset
cd c:\MYPROG\HAPACK\conda-recipe
echo %CD%
xcopy /y /s /e /h /d ..\mswin\Release_x64\molset %PREFIX%\molset
xcopy /y /s /e /h /d ..\mswin\Release_x64\harlem.exe %PREFIX%
xcopy /y /s /e /h /d ..\mswin\Release_x64\Lib\site-packages\wx %PREFIX%\wx
xcopy /y /s /e /h /d ..\examples %PREFIX%\examples
xcopy /y /s /e /h /d ..\residues_db %PREFIX%\residues_db
xcopy /y /s /e /h /d ..\doc\advanced_manual_html %PREFIX%\doc\advanced_manual_html
xcopy /y /s /e /h /d ..\doc\HARLEM_BeginnerUserManual.htm %PREFIX%\doc
xcopy /y /s /e /h /d ..\doc\HARLEM_BeginnerUserManual_files %PREFIX%\doc\HARLEM_BeginnerUserManual_files

