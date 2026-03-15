# MetaboDesk

Desktop application for constraint-based metabolic model analysis.

MetaboDesk provides a Windows GUI for loading SBML models, editing constraints, running core COBRA analyses, and exporting reproducible outputs.

## Download

- Latest installer and runtime assets: https://github.com/Emiray02/metabodesk/releases

## System Requirements

- Windows 10/11 (64-bit)
- 4 GB RAM minimum (8 GB recommended)
- ~1.5 GB free disk space
- Internet required only for first runtime download during install

## Quick Start

1. Install `MetaboDeskSetup_v1.0.0.exe`
2. Open MetaboDesk
3. Load an SBML file (`.xml` / `.sbml`)
4. Choose analysis type and run
5. Export results from Results / analysis tabs

## Core Features

### Analysis

- FBA
- pFBA
- FVA
- Single Gene Deletion (SGD)
- Double Gene Deletion (DGD)
- Single Reaction Deletion (SRD)
- Robustness Analysis
- Production Envelope
- Flux Sampling (ACHR / OptGP)

### Modeling Workflow

- SBML load/save
- Medium and reaction bound editing
- Gene knockout / reaction overexpression controls
- Scenario import/export (JSON)
- Recent scenarios menu

### Visualization & Export

- Tabular and chart-based result views
- Network map tab (interactive graph)
- Export to CSV / Excel / JSON / chart image

### Integrated Tools

- Memote integration
- CarveMe integration
- Community-model utility

## Packaging Model

- `dist/MetaboDesk.exe` is a thin launcher (PyInstaller onefile)
- Launcher starts bundled runtime Python at `runtime/python/pythonw.exe`
- Installer script: `MetaboDesk.iss`

## Development Setup

```bash
git clone https://github.com/Emiray02/metabodesk.git
cd metabodesk
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python -m metabodesk_core
```

## Build

### Thin launcher EXE

```bash
runtime\python\python.exe -m PyInstaller --clean --noconfirm launcher.spec
```

Output: `dist/MetaboDesk.exe`

### Installer

- Open `MetaboDesk.iss` in Inno Setup 6
- Compile to generate installer under `dist-installer/`

## Tests

```bash
.venv\Scripts\python.exe -m pytest tests/ -v --tb=short
```

## Project Structure

- `metabodesk_core/` main app package
- `tests/` pytest suite
- `runtime/` bundled Python runtime and wheels
- `launcher.py` thin launcher source
- `launcher.spec` PyInstaller spec for thin EXE
- `MetaboDesk.iss` Inno Setup script

## License

CC BY-NC-ND 4.0 — see `LICENSE`.

## Support

Issues and feature requests:
https://github.com/Emiray02/metabodesk/issues
