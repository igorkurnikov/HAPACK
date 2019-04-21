rem clean all builds but leave python and other dependencies
SET script_path=%~dp0
cd %script_path%

echo off

if exist "obj" (
    echo "Removing obj"
    rmdir  /Q /S "obj"
)

FOR %%D IN (PY3_Debug_Win32 PY3_Release_Win32 PY3_Debug_x64 PY3_Release_x64 Release_Win32 Release_x64) DO (
    if exist "%%D" (
        echo "Cleaning %%D"
        if exist "%%D\scripts" (
            rmdir  /Q /S "%%D\scripts"
        )
        if exist "%%D\harlempy" (
            rmdir  /Q /S "%%D\harlempy
        )
        del  /F /Q "%%D\*.lib"
        del  /F /Q "%%D\*.a"
        FOR %%F IN (*.lib harlem.exe harlem.exp harlem.ilk harlem.lib harlem.pdb) DO (
            if exist "%%D\%%F" (
                del  /F /Q "%%D\%%F"
            )
        )
    )
)

