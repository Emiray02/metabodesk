"""MetaboDesk Launcher — thin EXE wrapper.

This script is compiled to MetaboDesk.exe via PyInstaller (--onefile).
It does NOT bundle any scientific packages — it simply locates the
runtime Python and launches the real application with it.

Flow
----
1.  Find ``runtime/python/pythonw.exe`` next to this EXE.
2.  Launch ``metabodesk.py`` with that interpreter.
3.  Exit immediately (the app runs in its own process).
"""

import os
import sys
import subprocess
import time


def _base_dir() -> str:
    """Return the directory containing this EXE (or .py in dev mode)."""
    if getattr(sys, "frozen", False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))


def _show_error(title: str, msg: str) -> None:
    """Show a native Windows MessageBox (no Qt dependency)."""
    try:
        import ctypes
        ctypes.windll.user32.MessageBoxW(0, msg, title, 0x10)
    except Exception:
        print(f"ERROR: {title}\n{msg}", file=sys.stderr)


def _cleanup_stale_instances(base: str) -> None:
    """Best-effort cleanup of stale headless MetaboDesk pythonw processes.

    Targets only processes started from the same installation path and only
    those without a main window. Also keeps very recent processes (<15s) to
    avoid racing with a currently launching instance.
    """
    if not sys.platform.startswith("win"):
        return

    escaped_base = base.replace("'", "''")
    ps_script = rf"""
$ErrorActionPreference = 'SilentlyContinue'
$base = '{escaped_base}'
$targets = Get-CimInstance Win32_Process | Where-Object {{
    $_.Name -eq 'pythonw.exe' -and
    $_.CommandLine -like "*${{base}}\\metabodesk.py*"
}}
foreach ($p in $targets) {{
    $gp = Get-Process -Id $p.ProcessId -ErrorAction SilentlyContinue
    if ($null -eq $gp) {{ continue }}

    $ageOk = $true
    try {{
        $created = [Management.ManagementDateTimeConverter]::ToDateTime($p.CreationDate)
        if (((Get-Date) - $created).TotalSeconds -lt 15) {{ $ageOk = $false }}
    }} catch {{}}

    if ($ageOk -and $gp.MainWindowHandle -eq 0) {{
        Stop-Process -Id $p.ProcessId -Force -ErrorAction SilentlyContinue
    }}
}}
"""
    try:
        si = subprocess.STARTUPINFO()
        si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        si.wShowWindow = 0  # SW_HIDE
        subprocess.run(
            [
                "powershell.exe",
                "-NoProfile",
                "-ExecutionPolicy",
                "Bypass",
                "-Command",
                ps_script,
            ],
            check=False,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            startupinfo=si,
        )
        time.sleep(0.15)
    except Exception:
        pass


def main() -> None:
    base = _base_dir()

    # Clean stale hidden instances before launching a new one.
    _cleanup_stale_instances(base)

    # ── Locate runtime Python ────────────────────────────────────
    pythonw = os.path.join(base, "runtime", "python", "pythonw.exe")
    python  = os.path.join(base, "runtime", "python", "python.exe")

    # Create a branded copy of pythonw.exe so Task Manager shows
    # "MetaboDesk" instead of "Python" in the process name column.
    branded = os.path.join(base, "runtime", "python", "MetaboDesk.exe")
    if os.path.exists(pythonw) and not os.path.exists(branded):
        try:
            import shutil
            shutil.copy2(pythonw, branded)
        except Exception:
            branded = None

    if branded and os.path.exists(branded):
        interpreter = branded
    elif os.path.exists(pythonw):
        interpreter = pythonw
    elif os.path.exists(python):
        interpreter = python
    else:
        _show_error(
            "MetaboDesk — Runtime Not Found",
            "Python runtime was not found!\n\n"
            f"Expected location:\n{pythonw}\n\n"
            "Please re-run the MetaboDesk installer\n"
            "or manually download and extract runtime.zip.\n\n"
            "https://github.com/Emiray02/metabodesk/releases"
        )
        sys.exit(1)

    # ── Locate the main script ────────────────────────────────
    script = os.path.join(base, "metabodesk.py")
    if not os.path.exists(script):
        _show_error(
            "MetaboDesk — Script Not Found",
            "Application script was not found!\n\n"
            f"Expected location:\n{script}\n\n"
            "Please reinstall MetaboDesk."
        )
        sys.exit(1)

    # ── Launch the real application ──────────────────────────────
    env = os.environ.copy()
    env["METABODESK_BASE"] = base

    # Pass any command-line arguments (e.g. file to open)
    args = [interpreter, script] + sys.argv[1:]

    try:
        # Fully detach the child process so it survives after the
        # launcher EXE exits.  Use only CREATE_NEW_PROCESS_GROUP
        # (not DETACHED_PROCESS which can flash a console on some
        # Windows builds when combined with other flags).
        # IMPORTANT: Do not force SW_HIDE for pythonw/GUI interpreters,
        # otherwise the main Qt window can remain hidden.
        CREATE_NEW_PROCESS_GROUP = 0x00000200
        startupinfo = None
        if os.path.basename(interpreter).lower() == "python.exe":
            si = subprocess.STARTUPINFO()
            si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            si.wShowWindow = 0  # SW_HIDE (console interpreter only)
            startupinfo = si
        subprocess.Popen(
            args,
            cwd=base,
            env=env,
            creationflags=CREATE_NEW_PROCESS_GROUP,
            startupinfo=startupinfo,
            close_fds=True,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except Exception as exc:
        _show_error(
            "MetaboDesk — Launch Failed",
            f"Failed to launch the application:\n\n{exc}\n\n"
            f"Python: {interpreter}\n"
            f"Script: {script}"
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
