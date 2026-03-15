@echo off
REM ═══════════════════════════════════════════════════════════════
REM  MetaboDesk v1.0.0 Build Script (Thin Launcher)
REM  ─────────────────────────────────────────────────────────────
REM  Bu script sirasyla:
REM   1) PyInstaller ile ince launcher EXE olusturur (~7 MB)
REM   2) runtime.zip olusturur (GitHub Releases'e yuklenecek)
REM   3) Inno Setup ile installer olusturur
REM
REM  Launcher EXE sadece runtime\python\pythonw.exe'yi bulup
REM  metabodesk.py'yi baslatir. Agir bagimliliklar (cobra, numpy,
REM  PySide6...) runtime.zip icinde dagitilir.
REM ═══════════════════════════════════════════════════════════════
setlocal

echo.
echo ═══════════════════════════════════════════
echo  MetaboDesk v1.0.0 Release Builder
echo  (Thin Launcher Approach)
echo ═══════════════════════════════════════════
echo.

REM ── On kontroller ──
if not exist "launcher.py" (
    echo [HATA] launcher.py bulunamadi! Proje kok dizininde calistirin.
    pause
    exit /b 1
)
if not exist "launcher.spec" (
    echo [HATA] launcher.spec bulunamadi!
    pause
    exit /b 1
)
if not exist "metabodesk.py" (
    echo [HATA] metabodesk.py bulunamadi!
    pause
    exit /b 1
)
if not exist "logo.ico" (
    echo [HATA] logo.ico bulunamadi!
    pause
    exit /b 1
)

echo ═══════════════════════════════════════════
echo  ADIM 1/3: Thin Launcher EXE olustur
echo ═══════════════════════════════════════════
echo.

REM Onceki build ciktisini temizle
if exist "build" rmdir /s /q "build"
if exist "dist\MetaboDesk.exe" del "dist\MetaboDesk.exe"

REM PyInstaller ile ince launcher derle (onefile)
if exist "runtime\python\python.exe" (
    echo Runtime Python kullaniliyor...
    runtime\python\python.exe -m PyInstaller --clean --noconfirm launcher.spec
) else (
    echo Sistem Python kullaniliyor...
    pyinstaller --clean --noconfirm launcher.spec
)

if errorlevel 1 (
    echo.
    echo [HATA] PyInstaller basarisiz!
    echo Olasi nedenler:
    echo   - PyInstaller kurulu degil: pip install pyinstaller
    echo   - logo.ico eksik
    echo   - file_version_info.txt eksik
    pause
    exit /b 1
)

REM EXE'nin olustiguni dogrula
if not exist "dist\MetaboDesk.exe" (
    echo [HATA] MetaboDesk.exe olusturulamadi!
    pause
    exit /b 1
)
echo.
echo [OK] Thin launcher MetaboDesk.exe basariyla olusturuldu!
for %%A in ("dist\MetaboDesk.exe") do echo     Boyut: %%~zA bytes

echo.
echo ═══════════════════════════════════════════
echo  ADIM 2/3: runtime.zip olustur
echo ═══════════════════════════════════════════
echo.

if not exist "runtime" (
    echo [UYARI] runtime klasoru bulunamadi - runtime.zip olusturulmayacak.
    echo         GitHub Releases'e elle runtime.zip yuklemeniz gerekecek.
    goto :step3
)

echo runtime klasoru zipleniyor...
echo (Bu islem birkaç dakika surebilir — ~1.8 GB sikistiriliyor)
if exist "runtime.zip" del "runtime.zip"
powershell -NoProfile -Command "Compress-Archive -Path 'runtime' -DestinationPath 'runtime.zip' -CompressionLevel Optimal"
if errorlevel 1 (
    echo [HATA] runtime.zip olusturulamadi!
    pause
    exit /b 1
)
for %%A in (runtime.zip) do echo [OK] runtime.zip boyutu: %%~zA bytes

:step3
echo.
echo ═══════════════════════════════════════════
echo  ADIM 3/3: Inno Setup ile installer olustur
echo ═══════════════════════════════════════════
echo.

REM Inno Setup yolunu bul
set "ISCC="
if exist "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" (
    set "ISCC=C:\Program Files (x86)\Inno Setup 6\ISCC.exe"
) else if exist "C:\Program Files\Inno Setup 6\ISCC.exe" (
    set "ISCC=C:\Program Files\Inno Setup 6\ISCC.exe"
)

if "%ISCC%"=="" (
    echo [UYARI] Inno Setup 6 bulunamadi.
    echo         MetaboDesk.iss dosyasini Inno Setup ile elle derleyin.
    echo         https://jrsoftware.org/isinfo.php adresinden indirin.
    goto :done
)

echo Inno Setup: %ISCC%
if not exist "dist-installer" mkdir "dist-installer"

"%ISCC%" MetaboDesk.iss
if errorlevel 1 (
    echo.
    echo [HATA] Inno Setup basarisiz!
    pause
    exit /b 1
)
echo.
echo [OK] Installer basariyla olusturuldu!
for %%A in ("dist-installer\MetaboDeskSetup_v1.0.0.exe") do echo     Boyut: %%~zA bytes

:done
echo.
echo ═══════════════════════════════════════════
echo  TAMAMLANDI!
echo ═══════════════════════════════════════════
echo.
echo Cikti dosyalari:
echo   Launcher:  dist\MetaboDesk.exe (~7 MB)
if exist "runtime.zip" echo   Runtime:   runtime.zip (GitHub Releases'e yukleyin)
if exist "dist-installer\MetaboDeskSetup_v1.0.0.exe" echo   Installer: dist-installer\MetaboDeskSetup_v1.0.0.exe
echo.
echo Sonraki adimlar:
echo   1. runtime.zip dosyasini GitHub Releases'e yukleyin
echo      https://github.com/Emiray02/metabodesk/releases/tag/runtime
echo   2. MetaboDeskSetup_v1.0.0.exe dosyasini dagitin
echo.
echo Installer icerik:
echo   - MetaboDesk.exe (ince launcher)
echo   - metabodesk.py + metabodesk_core/*.py (kaynak kod)
echo   - logo.ico, splash.gif, LICENSE.txt
echo   - Kurulum sonrasi runtime.zip otomatik indirilir
echo.
pause
