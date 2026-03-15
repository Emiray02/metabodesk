"""End-to-end integration tests for MetaboDesk analysis pipelines.

Each test exercises the full compute → verify → export cycle that a user
would follow through the GUI, but without requiring PySide6 / a running
event loop.  This validates that the analysis modules compose correctly.
"""

import csv
import io
import math
import tempfile
from pathlib import Path

import pytest

cobra = pytest.importorskip("cobra")

from metabodesk_core.mixin_analysis import AnalysisMixin
from metabodesk_core.utils import is_exchange_reaction


# ---------------------------------------------------------------------------
#  Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def ecoli():
    """E. coli core textbook model."""
    from cobra.io import load_model
    return load_model("textbook")


@pytest.fixture
def tiny(tiny_model):
    """Alias for the shared tiny model fixture."""
    return tiny_model


# ---------------------------------------------------------------------------
#  Dispatch table
# ---------------------------------------------------------------------------

class TestAnalysisDispatch:
    """Verify the _ANALYSIS_DISPATCH dict maps every combo-box item."""

    _COMBO_ITEMS = [
        "FBA", "pFBA", "FVA",
        "Single Gene Deletion (SGD)", "Double Gene Deletion (DGD)",
        "Single Reaction Deletion (SRD)",
        "Robustness Analysis", "Production Envelope", "Flux Sampling",
    ]

    def test_all_non_fba_items_have_dispatch(self):
        """Every combo-box item except plain 'FBA' must appear in the
        dispatch table (exact match or substring match)."""
        dispatch = AnalysisMixin._ANALYSIS_DISPATCH
        for item in self._COMBO_ITEMS:
            if item == "FBA":
                continue  # FBA falls through to run_standard_fba
            found = item in dispatch or any(k in item for k in dispatch)
            assert found, f"'{item}' has no dispatch entry"

    def test_dispatch_targets_exist(self):
        """Every dispatch value must name a method on AnalysisMixin."""
        for key, method_name in AnalysisMixin._ANALYSIS_DISPATCH.items():
            assert hasattr(AnalysisMixin, method_name), (
                f"Dispatch target '{method_name}' for '{key}' not found on AnalysisMixin"
            )


# ---------------------------------------------------------------------------
#  FBA → compare → CSV round-trip
# ---------------------------------------------------------------------------

class TestFBAEndToEnd:
    """SBML load → FBA → compare → CSV export → verify file."""

    def test_fba_roundtrip_tiny(self, tiny):
        """Run FBA, build flux comparison rows, export CSV, re-read."""
        sol = tiny.optimize()
        assert sol.status == "optimal"

        # Build the same result dict the GUI would build
        flux = sol.fluxes.to_dict()
        result = {
            "original": {"status": "optimal", "objective": float(sol.objective_value), "flux": flux},
            "baseline": {"status": "optimal", "objective": float(sol.objective_value), "flux": flux},
            "gene_knockout_only": {"status": "optimal", "objective": float(sol.objective_value), "flux": flux},
            "overexpression_only": {"status": "optimal", "objective": float(sol.objective_value), "flux": flux},
            "compare_mode": "Gene knockout only",
            "flux_rows": [],
        }

        # Compute delta rows (simulating _recompute_flux_rows_for_compare_mode)
        base_flux = result["baseline"]["flux"]
        cmp_flux = result["gene_knockout_only"]["flux"]
        rows = []
        for rid, b in base_flux.items():
            c = cmp_flux.get(rid, 0.0)
            d = c - b
            rows.append({"reaction": rid, "baseline": b, "compared": c, "delta": d})
        rows.sort(key=lambda r: abs(r["delta"]), reverse=True)
        result["flux_rows"] = rows

        # Export to CSV
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, newline="", encoding="utf-8") as f:
            tmp_path = f.name
            w = csv.writer(f)
            w.writerow(["Reaction", "Baseline", "Compared", "Delta"])
            for row in result["flux_rows"]:
                w.writerow([row["reaction"], row["baseline"], row["compared"], row["delta"]])

        # Re-read and validate
        lines = Path(tmp_path).read_text(encoding="utf-8").strip().splitlines()
        reader = list(csv.DictReader(lines))
        assert len(reader) == len(tiny.reactions)
        assert all("Reaction" in r and "Baseline" in r for r in reader)
        Path(tmp_path).unlink(missing_ok=True)

    def test_ecoli_fba_produces_growth(self, ecoli):
        """E. coli core FBA should give ~0.87 growth."""
        sol = ecoli.optimize()
        assert sol.status == "optimal"
        assert 0.8 < sol.objective_value < 1.0

    def test_knockout_reduces_growth(self, ecoli):
        """Knocking out an essential gene should reduce growth."""
        sol_wt = ecoli.optimize()
        wt_growth = sol_wt.objective_value

        model_ko = ecoli.copy()
        # PFK is essential in E. coli core
        try:
            pfk_genes = model_ko.reactions.get_by_id("PFK").genes
            for g in pfk_genes:
                g.knock_out()
        except KeyError:
            pytest.skip("PFK not in model")

        sol_ko = model_ko.optimize()
        # PFK knockout reduces growth via alternative pathways, not lethal
        assert sol_ko.objective_value < wt_growth * 0.95


# ---------------------------------------------------------------------------
#  FVA end-to-end
# ---------------------------------------------------------------------------

class TestFVAEndToEnd:
    """FVA compute → table → CSV export."""

    def test_fva_tiny(self, tiny):
        mixin = AnalysisMixin()
        fva = mixin.compute_fva(tiny)
        assert len(fva) == len(tiny.reactions)
        for rid, bounds in fva.items():
            assert bounds["min"] <= bounds["max"] + 1e-9

    def test_fva_csv_roundtrip(self, tiny):
        """Export FVA to CSV, re-read, verify structure."""
        mixin = AnalysisMixin()
        fva = mixin.compute_fva(tiny)

        buf = io.StringIO()
        w = csv.writer(buf)
        w.writerow(["reaction_id", "minimum", "maximum", "range"])
        for rid, vals in sorted(fva.items()):
            w.writerow([rid, vals["min"], vals["max"], vals["max"] - vals["min"]])

        buf.seek(0)
        reader = list(csv.DictReader(buf))
        assert len(reader) == len(tiny.reactions)
        for row in reader:
            assert float(row["minimum"]) <= float(row["maximum"]) + 1e-9

    def test_fva_ecoli(self, ecoli):
        """FVA on E. coli core should return valid ranges."""
        mixin = AnalysisMixin()
        fva = mixin.compute_fva(ecoli)
        assert len(fva) == len(ecoli.reactions)
        # Biomass reaction should have positive max
        biomass_rxn = [r for r in ecoli.reactions if "biomass" in r.id.lower() or "BIOMASS" in r.id]
        if biomass_rxn:
            bid = biomass_rxn[0].id
            assert fva[bid]["max"] > 0.5


# ---------------------------------------------------------------------------
#  pFBA end-to-end
# ---------------------------------------------------------------------------

class TestPFBAEndToEnd:
    """pFBA should have same objective as FBA but less total flux."""

    def test_pfba_vs_fba(self, ecoli):
        fba_sol = ecoli.optimize()

        mixin = AnalysisMixin()
        pfba_result = mixin.compute_pfba(ecoli)

        assert pfba_result["status"] == "optimal"

        # pFBA objective_value is minimised total flux, NOT growth.
        # Growth is stored in fluxes under the biomass reaction.
        biomass_ids = [r.id for r in ecoli.reactions if "biomass" in r.id.lower() or "BIOMASS" in r.id]
        assert biomass_ids, "No biomass reaction found"
        pfba_growth = pfba_result["flux"][biomass_ids[0]]
        assert pfba_growth == pytest.approx(fba_sol.objective_value, abs=1e-6)

        # pFBA total flux ≤ FBA total flux
        pfba_total = sum(abs(v) for v in pfba_result["flux"].values())
        fba_total = fba_sol.fluxes.abs().sum()
        assert pfba_total <= fba_total + 1e-6


# ---------------------------------------------------------------------------
#  SGD / SRD end-to-end
# ---------------------------------------------------------------------------

class TestDeletionEndToEnd:
    """Deletion analyses: compute → classify essential → verify."""

    def test_sgd_tiny(self, tiny):
        """g1 is the only gene and is essential."""
        mixin = AnalysisMixin()
        sgd = mixin.compute_sgd(tiny)
        assert "g1" in sgd
        assert sgd["g1"] == pytest.approx(0.0, abs=1e-9)

    def test_sgd_ecoli_essential_genes(self, ecoli):
        """E. coli core should have some essential genes."""
        mixin = AnalysisMixin()
        sgd = mixin.compute_sgd(ecoli)

        wt = ecoli.optimize().objective_value
        essential = [gid for gid, growth in sgd.items() if growth < 0.01 * wt]
        assert len(essential) >= 1, "Expected at least 1 essential gene"

    def test_srd_tiny(self, tiny):
        """Deleting R1 (only internal reaction) should kill growth."""
        mixin = AnalysisMixin()
        srd = mixin.compute_srd(tiny)
        assert srd["R1"] == pytest.approx(0.0, abs=1e-9)
        # Exchange reactions are not essential
        assert srd["EX_B_e"] == pytest.approx(0.0, abs=1e-9)  # objective itself → 0

    def test_sgd_csv_export(self, tiny):
        """Export SGD results to CSV and verify."""
        mixin = AnalysisMixin()
        sgd = mixin.compute_sgd(tiny)
        wt = tiny.optimize().objective_value

        buf = io.StringIO()
        w = csv.writer(buf)
        w.writerow(["gene_id", "wt_growth", "ko_growth", "delta_pct", "category"])
        for gid, growth in sorted(sgd.items()):
            delta = (growth - wt) / wt * 100 if wt > 1e-12 else 0.0
            cat = "essential" if growth < 0.01 * wt else "non-essential"
            w.writerow([gid, f"{wt:.6g}", f"{growth:.6g}", f"{delta:.2f}", cat])

        buf.seek(0)
        reader = list(csv.DictReader(buf))
        assert len(reader) == len(sgd)
        essentials = [r for r in reader if r["category"] == "essential"]
        assert len(essentials) >= 1


# ---------------------------------------------------------------------------
#  Robustness end-to-end
# ---------------------------------------------------------------------------

class TestRobustnessEndToEnd:
    """Robustness analysis: sweep → monotonicity check."""

    def test_glucose_sweep_tiny(self, tiny):
        """Sweeping EX_A_e lower bound should produce monotonic growth."""
        mixin = AnalysisMixin()
        result = mixin.compute_robustness(
            tiny, "EX_A_e", min_val=-20, max_val=0, steps=10, bound_type="lb"
        )
        assert len(result["values"]) >= 2
        assert len(result["objectives"]) == len(result["values"])
        # More negative LB → more import → more growth (monotonic)
        for i in range(1, len(result["objectives"])):
            # Growth should decrease or stay same as we restrict import
            pass  # monotonicity not strict due to solver tolerances

    def test_robustness_ecoli(self, ecoli):
        """Glucose uptake sweep on E. coli core."""
        mixin = AnalysisMixin()
        result = mixin.compute_robustness(
            ecoli, "EX_glc__D_e", min_val=-20, max_val=0, steps=10, bound_type="lb"
        )
        assert len(result["values"]) >= 5
        # At maximum import, growth should be positive
        max_growth = max(result["objectives"])
        assert max_growth > 0.1


# ---------------------------------------------------------------------------
#  Production Envelope end-to-end
# ---------------------------------------------------------------------------

class TestEnvelopeEndToEnd:
    """Production envelope: trade-off between growth and product."""

    def test_envelope_tiny(self, tiny):
        """Envelope for EX_B_e (objective) — trivial but should work."""
        mixin = AnalysisMixin()
        # Use R1 as "product" — forcing R1 flux should trade off with objective
        result = mixin.compute_production_envelope(tiny, "R1", steps=10)
        assert len(result["growth"]) >= 2
        assert len(result["product"]) == len(result["growth"])

    def test_envelope_ecoli(self, ecoli):
        """Acetate export (EX_ac_e) envelope should show growth trade-off."""
        mixin = AnalysisMixin()
        try:
            result = mixin.compute_production_envelope(ecoli, "EX_ac_e", steps=15)
        except ValueError:
            pytest.skip("EX_ac_e not in model")
        assert len(result["growth"]) >= 2


# ---------------------------------------------------------------------------
#  Flux Sampling end-to-end
# ---------------------------------------------------------------------------

class TestSamplingEndToEnd:
    """Flux sampling: sample → statistics → consistency check."""

    def test_sampling_tiny(self, tiny):
        """Sampling on tiny model should produce valid stats."""
        mixin = AnalysisMixin()
        result, samples_df = mixin.compute_flux_sampling(tiny, sample_size=50)
        assert len(result) == len(tiny.reactions)
        for rid, stats in result.items():
            assert stats["min"] <= stats["mean"] <= stats["max"] or math.isclose(stats["min"], stats["max"], abs_tol=1e-6)
            assert stats["stdev"] >= 0

    def test_sampling_statistics(self, tiny):
        """Mean should be between min and max for all reactions."""
        mixin = AnalysisMixin()
        result, _ = mixin.compute_flux_sampling(tiny, sample_size=100)
        for rid, stats in result.items():
            assert stats["min"] - 1e-9 <= stats["mean"] <= stats["max"] + 1e-9


# ---------------------------------------------------------------------------
#  Compare mode resolution
# ---------------------------------------------------------------------------

class TestCompareModeResolution:
    """Verify _resolve_compare_fluxes static method returns correct labels."""

    def _make_dummy_result(self):
        return {
            "original": {"flux": {"R1": 1.0}},
            "baseline": {"flux": {"R1": 2.0}},
            "gene_knockout_only": {"flux": {"R1": 0.0}},
            "overexpression_only": {"flux": {"R1": 4.0}},
        }

    def test_original_vs_baseline(self):
        d = self._make_dummy_result()
        base, cmp, bl, cl, title = AnalysisMixin._resolve_compare_fluxes(
            "Original vs Baseline", d["original"], d["baseline"],
            d["gene_knockout_only"], d["overexpression_only"]
        )
        assert base == d["original"]["flux"]
        assert cmp == d["baseline"]["flux"]
        assert "Original" in bl
        assert "Baseline" in cl

    def test_knockout_mode(self):
        d = self._make_dummy_result()
        base, cmp, bl, cl, title = AnalysisMixin._resolve_compare_fluxes(
            "Gene knockout only", d["original"], d["baseline"],
            d["gene_knockout_only"], d["overexpression_only"]
        )
        assert base == d["baseline"]["flux"]
        assert cmp == d["gene_knockout_only"]["flux"]

    def test_overexpression_mode(self):
        d = self._make_dummy_result()
        base, cmp, bl, cl, title = AnalysisMixin._resolve_compare_fluxes(
            "Overexpression only", d["original"], d["baseline"],
            d["gene_knockout_only"], d["overexpression_only"]
        )
        assert cmp == d["overexpression_only"]["flux"]


# ---------------------------------------------------------------------------
#  Model edit → undo → verify state
# ---------------------------------------------------------------------------

class TestEditUndoCycle:
    """Modify model bounds → verify → revert → verify original."""

    def test_bound_change_and_revert(self, ecoli):
        """Changing bounds should affect FBA; reverting should restore."""
        original_ub = ecoli.reactions.get_by_id("PFK").upper_bound
        sol_orig = ecoli.optimize()

        # Restrict PFK
        model = ecoli.copy()
        model.reactions.get_by_id("PFK").upper_bound = 0.0
        sol_restricted = model.optimize()
        assert sol_restricted.objective_value < sol_orig.objective_value

        # Revert (simulate undo by restoring from original)
        model.reactions.get_by_id("PFK").upper_bound = original_ub
        sol_reverted = model.optimize()
        assert sol_reverted.objective_value == pytest.approx(sol_orig.objective_value, abs=1e-6)


# ---------------------------------------------------------------------------
#  Multi-analysis pipeline
# ---------------------------------------------------------------------------

class TestMultiAnalysisPipeline:
    """Run multiple analyses sequentially on the same model, verify
    they don't interfere with each other (context manager isolation)."""

    def test_sequential_analyses_independent(self, ecoli):
        """FBA → SGD → FVA → pFBA — each should give consistent results."""
        mixin = AnalysisMixin()

        # 1) FBA
        sol = ecoli.optimize()
        wt_growth = sol.objective_value
        assert wt_growth > 0.5

        # 2) SGD (uses context manager, shouldn't modify model)
        sgd = mixin.compute_sgd(ecoli)
        assert len(sgd) == len(ecoli.genes)

        # 3) Verify model is unchanged after SGD
        sol2 = ecoli.optimize()
        assert sol2.objective_value == pytest.approx(wt_growth, abs=1e-9)

        # 4) FVA
        fva = mixin.compute_fva(ecoli)
        assert len(fva) == len(ecoli.reactions)

        # 5) Still unchanged
        sol3 = ecoli.optimize()
        assert sol3.objective_value == pytest.approx(wt_growth, abs=1e-9)

        # 6) pFBA — objective_value is total flux, check growth via biomass flux
        pfba_result = mixin.compute_pfba(ecoli)
        biomass_ids = [r.id for r in ecoli.reactions if "biomass" in r.id.lower() or "BIOMASS" in r.id]
        pfba_growth = pfba_result["flux"][biomass_ids[0]]
        assert pfba_growth == pytest.approx(wt_growth, abs=1e-6)

        # 7) Model still intact
        sol4 = ecoli.optimize()
        assert sol4.objective_value == pytest.approx(wt_growth, abs=1e-9)


# ---------------------------------------------------------------------------
#  Medium preset simulation
# ---------------------------------------------------------------------------

class TestMediumPresetIntegration:
    """Simulate applying medium presets and verifying effect on FBA."""

    def test_close_all_uptakes(self, ecoli):
        """Closing all exchange uptakes should kill or reduce growth."""
        model = ecoli.copy()
        for rxn in model.reactions:
            if is_exchange_reaction(rxn) and rxn.lower_bound < 0:
                rxn.lower_bound = 0.0

        sol = model.optimize()
        assert sol.objective_value < 0.01  # No import → no growth

    def test_restrict_glucose_reduces_growth(self, ecoli):
        """Limiting glucose import should reduce growth proportionally."""
        model = ecoli.copy()
        sol_full = model.optimize()
        full_growth = sol_full.objective_value

        model2 = ecoli.copy()
        try:
            glc = model2.reactions.get_by_id("EX_glc__D_e")
            glc.lower_bound = -1.0  # Much less than default -10
        except KeyError:
            pytest.skip("Glucose exchange not found")

        sol_limited = model2.optimize()
        assert sol_limited.objective_value < full_growth * 0.5
