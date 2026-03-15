# MetaboDesk Release Notes

This file tracks user-facing changes by version.

## v1.0.0

### Added

- Windows desktop app for SBML-based metabolic analysis
- Core analysis suite: FBA, pFBA, FVA, SGD, DGD, SRD, Robustness, Envelope, Sampling
- Scenario import/export workflow (JSON)
- Results export (CSV/Excel/JSON and chart image export)
- Integrated tools tab (Memote, CarveMe)
- Network map tab and filtering controls

### Packaging

- Thin launcher architecture (`MetaboDesk.exe` + downloaded runtime)
- Installer built with Inno Setup (`MetaboDesk.iss`)

### Quality

- Automated pytest suite included under `tests/`

---

## Upcoming (next patch)

### Improved

- Startup reliability after splash screen
- Silent update-check behavior refined
- Installer behavior adjusted to avoid restart prompt
- Robustness analysis path hardened for solver stability

### Fixed

- Recent startup lock scenarios with hidden background process
- Sensitivity plotting edge cases with infeasible points
- BiGG model import action path issue
- FBA cache invalidation on medium/analysis parameter changes

---

For downloads and binaries:
https://github.com/Emiray02/metabodesk/releases
