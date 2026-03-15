# -*- mode: python ; coding: utf-8 -*-
# MetaboDesk v1.0.0 — Thin Launcher PyInstaller specification
#
# Bu spec SADECE launcher.py'yi derler.
# Sonuç: ~7 MB'lik tek bir MetaboDesk.exe (onefile)
# Ağır bağımlılıklar (cobra, numpy, PySide6…) dahil DEĞİLDİR —
# onlar runtime/python/ içinden çalışır.
#
# Build: pyinstaller --clean --noconfirm launcher.spec

a = Analysis(
    ['launcher.py'],
    pathex=[],
    binaries=[],
    datas=[],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        # Launcher'a hiçbir ağır paket dahil etme
        'numpy', 'scipy', 'pandas', 'matplotlib',
        'PySide6', 'cobra', 'networkx', 'openpyxl',
        'libsbml', 'swiglpk', 'optlang', 'sympy',
        'PIL', 'Pillow', 'docx', 'lxml',
        'tkinter', 'unittest', 'test', 'distutils',
        'setuptools', 'pip', 'wheel',
    ],
    noarchive=False,
    optimize=2,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='MetaboDesk',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['logo.ico'],
    version='file_version_info.txt',
)
