"""MetaboDesk - Metabolic Modelling Desktop Application.

This is the entry point.  The actual implementation lives in the
``metabodesk_core`` package (split into ~15 focused modules).

Kept as a single-file entry point for backward-compatibility with
PyInstaller (metabodesk.spec) and direct ``python metabodesk.py`` usage.
"""

from metabodesk_core.app import main

if __name__ == "__main__":
    main()
