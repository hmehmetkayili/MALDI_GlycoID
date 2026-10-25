@echo off
setlocal
title MALDI-GlycoID v1.0.0 — Local Web Application Launcher
cd /d "%~dp0"

echo ============================================================================
echo   MALDI-GlycoID v1.0.0 — MALDI-TOF MS N-Glycan ^& Glycopeptide Platform
echo ============================================================================
echo.

set "PYTHON_EXE="
if exist "%~dp0.venv\Scripts\python.exe" (
    set "PYTHON_EXE=%~dp0.venv\Scripts\python.exe"
) else if exist "C:\FragPipe\FragPipe-24.0\python\python.exe" (
    set "PYTHON_EXE=C:\FragPipe\FragPipe-24.0\python\python.exe"
) else (
    where python >nul 2>nul
    if %ERRORLEVEL% EQU 0 (
        set "PYTHON_EXE=python"
    )
)

if "%PYTHON_EXE%"=="" (
    echo [ERROR] Could not locate a Python 3.11+ interpreter.
    echo Please install Python 3.11+ or create a virtual environment in .venv
    pause
    exit /b 1
)

echo Using Python interpreter: "%PYTHON_EXE%"
echo Starting MALDI-GlycoID on http://127.0.0.1:8000 and opening browser...
echo (Press Ctrl+C in this window to stop the server)
echo.

"%PYTHON_EXE%" "%~dp0scripts\run_local.py" --host 127.0.0.1 --port 8000 --open-browser

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo [ERROR] MALDI-GlycoID exited with code %ERRORLEVEL%.
    pause
)
endlocal
