; MetaboDesk Installer (Inno Setup 6)
; ----------------------------------------------------------------
; Kurulum akisi (Thin Launcher - v2 yaklasimi):
;   1. MetaboDesk.exe (ince launcher, ~7 MB) + kaynak kodu kurar
;   2. Internetten runtime.zip indirir (Python + paketler + tools)
;   3. runtime.zip'i cikartir -> {app}\runtime\
;   4. Masaustu / Baslat Menusu kisayollari olusturur
;
; MetaboDesk.exe sadece bir launcher'dir - runtime\python\pythonw.exe
; uzerinden asil uygulamayi (metabodesk.py) baslatir.
; ----------------------------------------------------------------

#define MyAppName      "MetaboDesk"
#define MyAppVersion   "1.0.0"
#define MyAppPublisher "Emir Ay"
#define MyAppExeName   "MetaboDesk.exe"
#define MyAppURL       "https://github.com/Emiray02/metabodesk"

; -- Runtime indirme URL'si --
; GitHub Releases'e yuklenen runtime.zip dosyasinin adresi
#define RuntimeZipUrl  "https://github.com/Emiray02/metabodesk/releases/download/runtime/runtime.zip"

[Setup]
AppId={{3A4E6B64-3A63-4D8F-A4F2-0C2B9A8E5F9F}}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppVerName={#MyAppName} {#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}
AppUpdatesURL={#MyAppURL}/releases
DefaultDirName={autopf}\{#MyAppName}
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
InfoBeforeFile=beforeinstallation.txt
InfoAfterFile=Afterinstallation.txt
OutputDir=dist-installer
OutputBaseFilename=MetaboDeskSetup_v{#MyAppVersion}
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
SetupIconFile=logo.ico
UninstallDisplayIcon={app}\logo.ico
UninstallDisplayName={#MyAppName} {#MyAppVersion}
VersionInfoVersion={#MyAppVersion}.0
VersionInfoCompany={#MyAppPublisher}
VersionInfoDescription={#MyAppName} - Metabolic Modelling Desktop Application
VersionInfoProductName={#MyAppName}
VersionInfoProductVersion={#MyAppVersion}
PrivilegesRequired=admin
ArchitecturesInstallIn64BitMode=x64compatible
MinVersion=10.0
LicenseFile=LICENSE.txt
CloseApplications=yes
RestartApplications=no
AlwaysRestart=no
RestartIfNeededByRun=no
UninstallRestartComputer=no

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Files]
; -- Thin Launcher EXE (~7 MB) --
Source: "dist\MetaboDesk.exe"; DestDir: "{app}"; Flags: ignoreversion

; -- Uygulama kaynak kodu --
Source: "metabodesk.py";               DestDir: "{app}";                 Flags: ignoreversion
Source: "metabodesk_core\*.py";        DestDir: "{app}\metabodesk_core"; Flags: ignoreversion

; -- Logo ve splash animasyonu --
Source: "logo.ico";    DestDir: "{app}"; Flags: ignoreversion
Source: "splash.gif";  DestDir: "{app}"; Flags: ignoreversion

; -- Lisans --
Source: "LICENSE.txt";  DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{group}\{#MyAppName}";          Filename: "{app}\{#MyAppExeName}"; IconFilename: "{app}\logo.ico"
Name: "{group}\Uninstall {#MyAppName}"; Filename: "{uninstallexe}"
Name: "{commondesktop}\{#MyAppName}";  Filename: "{app}\{#MyAppExeName}"; IconFilename: "{app}\logo.ico"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop icon"; GroupDescription: "Additional icons:"

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "Launch {#MyAppName}"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
; Kurulum sirasinda indirilen runtime klasorunu kaldirma sirasinda temizle
Type: filesandordirs; Name: "{app}\runtime"
; Validation figurleri ve cache
Type: filesandordirs; Name: "{app}\validation\figures"
Type: filesandordirs; Name: "{app}\__pycache__"
Type: filesandordirs; Name: "{app}\metabodesk_core\__pycache__"

[Code]
// =================================================================
//  Inno Setup 6.1+ built-in download page with progress bar.
//  Shows: file name, MB downloaded / total MB, speed, progress %.
// =================================================================

var
  DownloadPage: TDownloadWizardPage;
  RuntimeAlreadyExists: Boolean;

function CloseRunningMetaboDesk(const AppDir: String): Boolean;
var
  ResultCode: Integer;
  Cmd: String;
begin
  Result := True;

  Cmd := '-NoProfile -ExecutionPolicy Bypass -Command ' +
    '"$ErrorActionPreference = ''''SilentlyContinue''''; ' +
    '$app = ''' + AppDir + '''; ' +
    '$targets = Get-CimInstance Win32_Process | Where-Object { ' +
    '($_.Name -eq ''''pythonw.exe'''' -and $_.CommandLine -like (''''*'''' + $app + ''''\\metabodesk.py*'''')) -or ' +
    '($_.Name -eq ''''MetaboDesk.exe'''' -and $_.ExecutablePath -eq (Join-Path $app ''''MetaboDesk.exe'''')) ' +
    '}; ' +
    'foreach ($p in $targets) { Stop-Process -Id $p.ProcessId -Force -ErrorAction SilentlyContinue }; ' +
    'Start-Sleep -Milliseconds 400"';

  if not Exec('powershell.exe', Cmd, '', SW_HIDE, ewWaitUntilTerminated, ResultCode) then
    Result := False
  else
    Result := (ResultCode = 0);
end;

// ----- Download progress callback (shows MB + percent) -----
function OnDownloadProgress(const Url, FileName: string; const Progress, ProgressMax: Int64): Boolean;
var
  ProgressMB, TotalMB: Double;
  Pct: Integer;
  Msg: String;
begin
  if ProgressMax <> 0 then
  begin
    ProgressMB := Progress / 1048576.0;
    TotalMB    := ProgressMax / 1048576.0;
    Pct        := (Progress * 100) div ProgressMax;
    Msg := Format('%.1f MB / %.1f MB   (%d%%)', [ProgressMB, TotalMB, Pct]);
    DownloadPage.SetText('Downloading runtime.zip ...', Msg);
    DownloadPage.SetProgress(Progress, ProgressMax);
  end;
  Result := True;
end;

// ----- Create download page when wizard initializes -----
procedure InitializeWizard;
begin
  DownloadPage := CreateDownloadPage(
    'Downloading Runtime Environment',
    'Please wait while the Python runtime environment is being downloaded...',
    @OnDownloadProgress);
end;

// ----- Start download when "Next" is clicked -----
function NextButtonClick(CurPageID: Integer): Boolean;
var
  PythonExe: String;
  AppDir: String;
begin
  Result := True;

  if CurPageID = wpReady then
  begin
    AppDir := ExpandConstant('{app}');
    CloseRunningMetaboDesk(AppDir);

    // Upgrade: skip download if runtime already exists
    PythonExe := AppDir + '\runtime\python\python.exe';
    if FileExists(PythonExe) then
    begin
      RuntimeAlreadyExists := True;
      Exit;
    end;

    RuntimeAlreadyExists := False;
    DownloadPage.Clear;
    DownloadPage.Add(ExpandConstant('{#RuntimeZipUrl}'), 'runtime.zip', '');
    DownloadPage.Show;
    try
      try
        DownloadPage.Download;
        Result := True;
      except
        if DownloadPage.AbortedByUser then
          Log('Download aborted by user.')
        else
          SuppressibleMsgBox(
            'Runtime download failed.' + #13#10 +
            'Please check your internet connection and try again.' + #13#10#13#10 +
            'Error: ' + AddPeriod(GetExceptionMessage) + #13#10#13#10 +
            'You can also manually download runtime.zip from:' + #13#10 +
            ExpandConstant('{#RuntimeZipUrl}'),
            mbCriticalError, MB_OK, IDOK);
        Result := False;
      end;
    finally
      DownloadPage.Hide;
    end;
  end;
end;

// ----- Extract ZIP (tar + PowerShell fallback) -----
function ExtractZip(const ZipPath, DestDir: String): Boolean;
var
  ResultCode: Integer;
  Cmd: String;
begin
  Result := False;

  // Method 1: tar (built-in on Windows 10+)
  Cmd := '-xf "' + ZipPath + '" -C "' + DestDir + '"';
  if Exec('tar', Cmd, '', SW_HIDE, ewWaitUntilTerminated, ResultCode) then
  begin
    if ResultCode = 0 then
    begin
      Result := True;
      Exit;
    end;
  end;

  // Method 2: PowerShell fallback
  Cmd := '-NoProfile -ExecutionPolicy Bypass -Command "Add-Type -AssemblyName System.IO.Compression.FileSystem; ' +
         '[System.IO.Compression.ZipFile]::ExtractToDirectory(''' + ZipPath + ''', ''' + DestDir + ''')"';
  if Exec('powershell.exe', Cmd, '', SW_HIDE, ewWaitUntilTerminated, ResultCode) then
    Result := (ResultCode = 0);
end;

// ----- Post-install: extract downloaded ZIP -----
procedure CurStepChanged(CurStep: TSetupStep);
var
  AppDir, ZipPath: String;
begin
  if CurStep = ssPostInstall then
  begin
    if RuntimeAlreadyExists then
    begin
      WizardForm.StatusLabel.Caption := 'Existing runtime found - skipping extraction.';
      WizardForm.Refresh();
      Exit;
    end;

    AppDir  := ExpandConstant('{app}');
    ZipPath := ExpandConstant('{tmp}\runtime.zip');

    if not FileExists(ZipPath) then
    begin
      MsgBox('Runtime archive not found. The download may have failed.' + #13#10 +
             'Please re-run the installer.',
             mbError, MB_OK);
      Exit;
    end;

    // Extract ZIP -> {app}\runtime\
    WizardForm.StatusLabel.Caption := 'Extracting runtime files (this may take a few minutes)...';
    WizardForm.Refresh();
    WizardForm.ProgressGauge.Style := npbstMarquee;

    if not ExtractZip(ZipPath, AppDir) then
    begin
      MsgBox('Failed to extract runtime files.' + #13#10 +
             'Please check disk space and try again.',
             mbError, MB_OK);
      Exit;
    end;

    // Delete temp file
    DeleteFile(ZipPath);

    WizardForm.StatusLabel.Caption := 'Installation complete!';
    WizardForm.ProgressGauge.Style := npbstNormal;
    WizardForm.ProgressGauge.Position := WizardForm.ProgressGauge.Max;
  end;
end;
