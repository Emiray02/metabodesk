# -*- mode: python ; coding: utf-8 -*-
# MetaboDesk v1.0.0 — PyInstaller build specification
# Build: pyinstaller --clean --noconfirm metabodesk.spec

a = Analysis(
    ['metabodesk.py'],
    pathex=[],
    binaries=[],
    datas=[
        ('logo.ico', '.'),
        ('splash.gif', '.'),
        ('metabodesk_core', 'metabodesk_core'),
    ],
    hiddenimports=[
        # ── MetaboDesk core modules ──
        'metabodesk_core',
        'metabodesk_core.__main__',
        'metabodesk_core.app',
        'metabodesk_core.config',
        'metabodesk_core.constants',
        'metabodesk_core.utils',
        'metabodesk_core.widgets',
        'metabodesk_core.mainwindow',
        # ── Mixin modules ──
        'metabodesk_core.mixin_io',
        'metabodesk_core.mixin_medium',
        'metabodesk_core.mixin_reactions',
        'metabodesk_core.mixin_analysis',
        'metabodesk_core.mixin_advanced',
        'metabodesk_core.mixin_advanced_omics',
        'metabodesk_core.mixin_advanced_design',
        'metabodesk_core.mixin_editor',
        'metabodesk_core.mixin_network',
        'metabodesk_core.mixin_tools',
        'metabodesk_core.mixin_dialogs',
        # ── COBRApy & solvers ──
        'cobra',
        'cobra.io',
        'cobra.io.sbml',
        'cobra.flux_analysis',
        'cobra.flux_analysis.variability',
        'cobra.flux_analysis.deletion',
        'cobra.flux_analysis.loopless',
        'cobra.medium',
        'cobra.util',
        'cobra.util.solver',
        'optlang',
        'optlang.glpk_interface',
        'swiglpk',
        # ── Numerical / Data ──
        'numpy',
        'numpy.core',
        'pandas',
        'scipy',
        'scipy.sparse',
        'scipy.optimize',
        'scipy.stats',
        # ── Plotting ──
        'matplotlib',
        'matplotlib.backends',
        'matplotlib.backends.backend_qtagg',
        'matplotlib.figure',
        # ── PySide6 ──
        'PySide6',
        'PySide6.QtCore',
        'PySide6.QtGui',
        'PySide6.QtWidgets',
        # ── Network ──
        'networkx',
        # ── Optional (lazy imports) ──
        'openpyxl',
        'libsbml',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'tkinter',
        'unittest',
        'test',
        'distutils',
    ],
    noarchive=False,
    optimize=0,
)

pyz = PYZ(a.pure)

splash = Splash(
    'splash.gif',
    binaries=a.binaries,
    datas=a.datas,
    text_pos=None,
    text_size=12,
    minify_script=True,
    always_on_top=True,
)

exe = EXE(
    pyz,
    a.scripts,
    splash,
    [],
    exclude_binaries=True,
    name='MetaboDesk',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['logo.ico'],
    version='file_version_info.txt',
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    splash.binaries,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='MetaboDesk',
)
