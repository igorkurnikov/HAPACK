@echo off
SET initial_path=%cd%
SET script_path=%~dp0
cd %script_path%

SET config=PY3_Release_x64

REM Make directories for testing
if exist "%script_path%%config%\pnpstests\" (
    echo "%script_path%%config%\pnpstests already exists, will delete it and create again"
    rmdir "%script_path%%config%\pnpstests" /s /q
    mkdir "%script_path%%config%\pnpstests"
) else (
    mkdir "%script_path%%config%\pnpstests"
)

cd "%script_path%%config%\pnpstests"

%script_path%%config%\python.exe -m pytest C:\MYPROG\HAPACK\PNPS\regtests

cd %initial_path%
