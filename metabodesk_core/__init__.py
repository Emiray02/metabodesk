"""MetaboDesk core package — modular metabolic modelling desktop application.

Architecture
------------
The UI is built with PySide6 (Qt 6) and metabolic computations use COBRApy.
Functionality is split across mixin classes that are composed into a
single :class:`~metabodesk_core.mainwindow.MainWindow` via multiple
inheritance (C3-linearised MRO):

- :mod:`mixin_io` — SBML I/O, scenario import/export
- :mod:`mixin_medium` — growth-medium table and presets
- :mod:`mixin_reactions` — reaction/gene tables, knockouts, overexpression
- :mod:`mixin_analysis` — FBA, FVA, pFBA, SGD, DGD, SRD, robustness, etc.
- :mod:`mixin_advanced` — thin facade combining the two sub-mixins below:
    - :mod:`mixin_advanced_omics` — Gene Expression (E-Flux/GIMME/iMAT),
      GECKO enzyme-constrained FBA, Thermodynamic FBA (TMFA), FSEOF
    - :mod:`mixin_advanced_design` — OptKnock, dFBA, Flux Coupling,
      Gap-Filling, PhPP, Model Comparison, Pathway Enrichment, Escher,
      SBML Validation, Jupyter Export
- :mod:`mixin_editor` — patch-based model editing
- :mod:`mixin_network` — interactive network-map visualisation
- :mod:`mixin_tools` — memote, CarveMe, community model builder
- :mod:`mixin_dialogs` — statistics, batch/sensitivity analysis, undo/redo
"""

from metabodesk_core.app import main  # noqa: F401

__all__ = ["main"]
