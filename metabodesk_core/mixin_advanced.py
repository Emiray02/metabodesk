"""Advanced analyses mixin for MetaboDesk -- thin facade.

This module re-exports the `AdvancedMixin` class which combines:

- **OmicsMixin** (`mixin_advanced_omics`): Gene Expression Integration,
  GECKO enzyme-constrained FBA, Thermodynamic FBA (TMFA), FSEOF.
- **DesignMixin** (`mixin_advanced_design`): OptKnock, Dynamic FBA,
  Flux Coupling, Gap-Filling, Model Comparison, Pathway Enrichment,
  Escher Map Viewer, Phenotype Phase Plane, SBML Validation,
  Jupyter Notebook Export.

All downstream code can continue to `from mixin_advanced import AdvancedMixin`
without any changes.
"""

from metabodesk_core.mixin_advanced_omics import OmicsMixin
from metabodesk_core.mixin_advanced_design import DesignMixin


class AdvancedMixin(OmicsMixin, DesignMixin):
    """Composite mixin that combines OmicsMixin and DesignMixin."""
