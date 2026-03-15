@echo off
REM MetaboDesk v14.0.0 — Quick EXE build (spec dosyasından)
REM Tam build için build_release.bat kullanın.
echo.
echo MetaboDesk PyInstaller Build
echo ════════════════════════════
echo.

if exist "runtime\python\python.exe" (
    echo Runtime Python ile derleniyor...
    runtime\python\python.exe -m PyInstaller --clean --noconfirm metabodesk.spec
) else (
    echo Sistem Python ile derleniyor...
    pyinstaller --clean --noconfirm metabodesk.spec
)

if errorlevel 1 (
    echo.
    echo [HATA] Build basarisiz!
    pause
    exit /b 1
)

echo.
echo [OK] Build tamamlandi: dist\MetaboDesk\MetaboDesk.exe
pause