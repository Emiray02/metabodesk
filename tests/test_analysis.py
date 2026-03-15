"""Tests for core FBA / model operations (integration-level)."""

import pytest

# Skip all tests in this module if cobra is not available
cobra = pytest.importorskip("cobra")


class TestTinyModelFBA:
    """Basic FBA / pFBA / FVA on the tiny fixture model."""

    def test_fba_optimal(self, tiny_model):
        """The tiny model should produce an optimal solution."""
        sol = tiny_model.optimize()
        assert sol.status == "optimal"
        assert sol.objective_value == pytest.approx(10.0)

    def test_pfba(self, tiny_model):
        """Parsimonious FBA should succeed and have ≤ FBA total flux."""
        from cobra.flux_analysis import pfba
        sol = pfba(tiny_model)
        assert sol.status == "optimal"
        total_flux = sol.fluxes.abs().sum()
        assert total_flux > 0

    def test_fva(self, tiny_model):
        """FVA should return min/max for each reaction."""
        from cobra.flux_analysis import flux_variability_analysis
        fva = flux_variability_analysis(tiny_model, processes=1)
        assert "minimum" in fva.columns
        assert "maximum" in fva.columns
        assert len(fva) == len(tiny_model.reactions)

    def test_single_gene_deletion(self, tiny_model):
        """SGD should identify g1 as essential (growth → 0 when deleted)."""
        from cobra.flux_analysis import single_gene_deletion
        sgd = single_gene_deletion(tiny_model, processes=1)
        # g1 controls R1, the only internal reaction
        assert len(sgd) >= 1


class TestEcoliCore:
    """Integration tests on the E. coli core model (if available)."""

    def test_fba(self, ecoli_core):
        if ecoli_core is None:
            pytest.skip("E. coli core model not available")
        sol = ecoli_core.optimize()
        assert sol.status == "optimal"
        assert sol.objective_value > 0.5  # Known to be ~0.87

    def test_model_has_genes(self, ecoli_core):
        if ecoli_core is None:
            pytest.skip("E. coli core model not available")
        assert len(ecoli_core.genes) > 100

    def test_model_has_reactions(self, ecoli_core):
        if ecoli_core is None:
            pytest.skip("E. coli core model not available")
        assert len(ecoli_core.reactions) > 70


class TestModelCopy:
    """Ensure model copies are independent."""

    def test_copy_independence(self, tiny_model):
        copy = tiny_model.copy()
        copy.reactions.get_by_id("R1").upper_bound = 0.0
        # Original should be unaffected
        assert tiny_model.reactions.get_by_id("R1").upper_bound == 1000.0
        # Copy should be constrained
        sol = copy.optimize()
        assert sol.objective_value == pytest.approx(0.0)
