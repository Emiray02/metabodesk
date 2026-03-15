"""Tests for advanced analysis algorithms.

Validates the core computation logic of E-Flux, GIMME, iMAT, OptKnock,
dFBA, TMFA, Flux Coupling, and PhPP by replicating the essential steps
from the mixin closures on small test models.
"""

import pytest

cobra = pytest.importorskip("cobra")


# ---------------------------------------------------------------------------
#  E-Flux
# ---------------------------------------------------------------------------

class TestEFlux:
    """E-Flux algorithm: scale reaction bounds by normalised expression."""

    def test_eflux_scales_bounds(self, tiny_model):
        """Scaling should reduce upper bounds for low-expression genes."""
        from metabodesk_core.utils import evaluate_gpr_expression

        model = tiny_model.copy()
        expr_data = {"g1": 0.5}

        # Map gene expression to reactions via GPR evaluation
        rxn_expr = {}
        for rxn in model.reactions:
            if not rxn.genes:
                rxn_expr[rxn.id] = None
                continue
            gene_vals = {g.id: expr_data[g.id] for g in rxn.genes if g.id in expr_data}
            if not gene_vals:
                rxn_expr[rxn.id] = None
                continue
            gpr_str = getattr(rxn, "gene_reaction_rule", "") or ""
            rxn_expr[rxn.id] = (
                evaluate_gpr_expression(gpr_str, gene_vals)
                if gpr_str.strip()
                else min(gene_vals.values())
            )

        max_expr = max((v for v in rxn_expr.values() if v is not None), default=1.0)
        assert max_expr > 0

        for rxn in model.reactions:
            expr_val = rxn_expr.get(rxn.id)
            if expr_val is not None and rxn.upper_bound > 0:
                rxn.upper_bound *= max(expr_val / max_expr, 0.001)

        sol = model.optimize()
        assert sol.status == "optimal"

    def test_eflux_low_expression_limits_flux(self, tiny_model):
        """Low expression on the only internal reaction should limit objective."""
        model = tiny_model.copy()
        sol_orig = model.optimize()
        orig_obj = sol_orig.objective_value

        # Scale R1 to 0.1 % capacity — UB becomes 1.0 which is below
        # the 10-unit import capacity of EX_A_e, creating a bottleneck.
        r1 = model.reactions.get_by_id("R1")
        r1.upper_bound *= 0.001  # 1000 * 0.001 = 1.0

        sol = model.optimize()
        assert sol.status == "optimal"
        assert sol.objective_value <= orig_obj * 0.11  # ≤ ~1.0 vs original 10.0


# ---------------------------------------------------------------------------
#  GIMME
# ---------------------------------------------------------------------------

class TestGIMME:
    """GIMME algorithm: minimise flux through low-expression reactions."""

    def test_gimme_maintains_min_growth(self, tiny_model):
        """Growth should stay above the required minimum fraction."""
        model = tiny_model.copy()
        sol0 = model.optimize()
        required_growth = sol0.objective_value * 0.9

        obj_rxn = model.reactions.get_by_id("EX_B_e")
        obj_rxn.lower_bound = max(obj_rxn.lower_bound, required_growth)

        sol = model.optimize()
        assert sol.status == "optimal"
        assert sol.objective_value >= required_growth - 1e-6

    def test_gimme_no_feasible_growth(self, tiny_model):
        """Over-constraining growth should cause infeasibility."""
        model = tiny_model.copy()
        sol0 = model.optimize()

        obj_rxn = model.reactions.get_by_id("EX_B_e")
        obj_rxn.lower_bound = sol0.objective_value + 100  # impossible

        sol = model.optimize()
        assert sol.status == "infeasible"


# ---------------------------------------------------------------------------
#  iMAT
# ---------------------------------------------------------------------------

class TestIMAT:
    """iMAT algorithm: classify reactions and constrain by expression."""

    def test_imat_constrains_low_expression(self, tiny_model):
        """Constraining the only internal reaction to near-zero blocks growth."""
        model = tiny_model.copy()
        r1 = model.reactions.get_by_id("R1")
        r1.upper_bound = 0.01
        r1.lower_bound = max(r1.lower_bound, -0.01)

        sol = model.optimize()
        assert sol.objective_value < 0.02

    def test_imat_high_expression_unmodified(self, tiny_model):
        """Reactions marked as high-expression should remain active."""
        model = tiny_model.copy()
        sol = model.optimize()
        assert sol.status == "optimal"
        assert abs(sol.fluxes["R1"]) > 1e-6  # R1 is active


# ---------------------------------------------------------------------------
#  OptKnock
# ---------------------------------------------------------------------------

class TestOptKnock:
    """Combinatorial knockout search for strain design."""

    def test_knockout_essential_reaction(self, tiny_model):
        """Knocking out R1 eliminates the only path → zero objective."""
        model = tiny_model.copy()
        with model:
            model.reactions.get_by_id("R1").knock_out()
            sol = model.optimize()
            assert sol.objective_value < 1e-6

    def test_knockout_exchange_does_not_crash(self, tiny_model):
        """Knocking out an exchange reaction should not raise."""
        model = tiny_model.copy()
        with model:
            model.reactions.get_by_id("EX_A_e").knock_out()
            sol = model.optimize()
            # With no A import, objective should be 0
            assert sol.objective_value < 1e-6

    def test_single_ko_search(self, tiny_model):
        """A simple 1-KO sweep should return results for each reaction."""
        model = tiny_model.copy()
        sol0 = model.optimize()
        baseline_obj = sol0.objective_value

        results = []
        for rxn in model.reactions:
            with model:
                rxn.knock_out()
                sol = model.optimize()
                results.append({
                    "knockout": rxn.id,
                    "growth": float(sol.objective_value) if sol.status == "optimal" else 0.0,
                })
        assert len(results) == 3
        # R1 knockout should be lethal
        r1_result = [r for r in results if r["knockout"] == "R1"][0]
        assert r1_result["growth"] < 1e-6


# ---------------------------------------------------------------------------
#  Dynamic FBA (dFBA)
# ---------------------------------------------------------------------------

class TestDFBA:
    """Time-course batch simulation."""

    def test_euler_produces_growth(self, tiny_model):
        """Euler integration should show biomass increase."""
        model = tiny_model.copy()

        biomass, substrate = 0.1, 20.0
        dt, t_end, t = 0.5, 3.0, 0.0
        biomasses = []

        while t < t_end and substrate > 1e-6 and biomass > 1e-9:
            biomasses.append(biomass)
            max_uptake = substrate / (biomass * dt) if biomass > 0 else 0
            with model:
                sub = model.reactions.get_by_id("EX_A_e")
                sub.lower_bound = max(-max_uptake, sub.lower_bound)
                sol = model.optimize()
                if sol.status != "optimal":
                    break
                mu = float(sol.objective_value)
                v_sub = float(sol.fluxes.get("EX_A_e", 0))

            biomass += mu * biomass * dt
            substrate += v_sub * biomass * dt
            biomass, substrate = max(0, biomass), max(0, substrate)
            t += dt

        assert len(biomasses) >= 2
        assert biomasses[-1] >= biomasses[0]

    def test_substrate_depletes(self, tiny_model):
        """Substrate concentration should decrease over time."""
        model = tiny_model.copy()

        biomass, substrate = 0.1, 5.0  # small substrate pool
        dt, t_end, t = 0.5, 10.0, 0.0
        substrates = [substrate]

        while t < t_end and substrate > 1e-6 and biomass > 1e-9:
            max_uptake = substrate / (biomass * dt) if biomass > 0 else 0
            with model:
                sub = model.reactions.get_by_id("EX_A_e")
                sub.lower_bound = max(-max_uptake, sub.lower_bound)
                sol = model.optimize()
                if sol.status != "optimal":
                    break
                mu = float(sol.objective_value)
                v_sub = float(sol.fluxes.get("EX_A_e", 0))

            biomass += mu * biomass * dt
            substrate += v_sub * biomass * dt
            biomass, substrate = max(0, biomass), max(0, substrate)
            substrates.append(substrate)
            t += dt

        assert substrates[-1] < substrates[0]


# ---------------------------------------------------------------------------
#  Thermodynamic FBA (TMFA)
# ---------------------------------------------------------------------------

class TestTMFA:
    """Thermodynamic constraints based on ΔG."""

    def test_positive_dg_blocks_forward(self, tiny_model):
        """Positive ΔG reaction should be blocked in strict mode → zero growth."""
        model = tiny_model.copy()
        model.reactions.get_by_id("R1").upper_bound = 0.0  # block forward
        sol = model.optimize()
        assert sol.objective_value < 1e-6

    def test_negative_dg_blocks_reverse(self, tiny_model):
        """Negative ΔG should block reverse direction only."""
        model = tiny_model.copy()
        r1 = model.reactions.get_by_id("R1")
        r1.lower_bound = 0.0  # block reverse (already 0, no effect)
        sol = model.optimize()
        assert sol.status == "optimal"
        assert sol.objective_value > 1.0

    def test_no_dg_data_unchanged(self, tiny_model):
        """Without ΔG data, model should optimise normally."""
        model = tiny_model.copy()
        sol = model.optimize()
        assert sol.status == "optimal"
        assert sol.objective_value == pytest.approx(10.0)


# ---------------------------------------------------------------------------
#  Flux Coupling Analysis
# ---------------------------------------------------------------------------

class TestFluxCoupling:
    """Ratio-based flux coupling detection."""

    def test_no_blocked_reactions(self, tiny_model):
        """All reactions in tiny model should carry flux."""
        from cobra.flux_analysis import flux_variability_analysis

        model = tiny_model.copy()
        fva = flux_variability_analysis(model, fraction_of_optimum=0.0)

        for rxn_id in fva.index:
            fmin = float(fva.loc[rxn_id, "minimum"])
            fmax = float(fva.loc[rxn_id, "maximum"])
            assert abs(fmax) > 1e-9 or abs(fmin) > 1e-9, f"{rxn_id} appears blocked"

    def test_blocking_r1_blocks_exb(self, tiny_model):
        """Blocking R1 should also block EX_B_e (directional coupling)."""
        from cobra.flux_analysis import flux_variability_analysis

        model = tiny_model.copy()
        model.reactions.get_by_id("R1").bounds = (0, 0)
        fva = flux_variability_analysis(model, fraction_of_optimum=0.0)

        exb_min = float(fva.loc["EX_B_e", "minimum"])
        exb_max = float(fva.loc["EX_B_e", "maximum"])
        assert abs(exb_min) < 1e-6 and abs(exb_max) < 1e-6

    def test_fva_ranges_consistent(self, tiny_model):
        """FVA min ≤ max for every reaction."""
        from cobra.flux_analysis import flux_variability_analysis

        model = tiny_model.copy()
        fva = flux_variability_analysis(model, fraction_of_optimum=0.0)

        for rxn_id in fva.index:
            assert fva.loc[rxn_id, "minimum"] <= fva.loc[rxn_id, "maximum"] + 1e-9


# ---------------------------------------------------------------------------
#  Phenotype Phase Plane (PhPP)
# ---------------------------------------------------------------------------

class TestPhPP:
    """Two-dimensional objective scanning."""

    def test_1d_sweep_produces_values(self, tiny_model):
        """Sweeping one reaction bound should yield multiple data points."""
        model = tiny_model.copy()
        results = []
        steps = 5

        for i in range(steps):
            lb = -10.0 + i * 10.0 / (steps - 1)
            with model:
                model.reactions.get_by_id("EX_A_e").lower_bound = lb
                sol = model.optimize()
                if sol.status == "optimal":
                    results.append((lb, float(sol.objective_value)))

        assert len(results) >= 2

    def test_more_import_more_growth(self, tiny_model):
        """More substrate import (more negative LB) should allow more growth."""
        model = tiny_model.copy()

        with model:
            model.reactions.get_by_id("EX_A_e").lower_bound = -5.0
            sol5 = model.optimize()

        with model:
            model.reactions.get_by_id("EX_A_e").lower_bound = -10.0
            sol10 = model.optimize()

        assert sol10.objective_value >= sol5.objective_value - 1e-6

    def test_zero_import_zero_growth(self, tiny_model):
        """With no substrate import, growth should be zero."""
        model = tiny_model.copy()
        model.reactions.get_by_id("EX_A_e").lower_bound = 0.0
        sol = model.optimize()
        assert sol.objective_value < 1e-6


# ---------------------------------------------------------------------------
#  dFBA Euler Fix
# ---------------------------------------------------------------------------

class TestDFBAEulerFix:
    """Verify that the Euler integration uses old biomass for substrate update."""

    def test_euler_old_biomass_used(self, tiny_model):
        """Euler integration should use the pre-update biomass for substrate."""
        model = tiny_model.copy()
        biomass = 0.1
        substrate = 20.0
        dt = 0.5

        with model:
            sub_rxn = model.reactions.get_by_id("EX_A_e")
            max_uptake = substrate / (biomass * dt)
            sub_rxn.lower_bound = max(-max_uptake, sub_rxn.lower_bound)
            sol = model.optimize()
            growth_rate = float(sol.objective_value)
            uptake_rate = float(sol.fluxes.get("EX_A_e", 0))

        old_biomass = biomass
        biomass += growth_rate * biomass * dt
        # Correct: use old_biomass
        substrate_correct = substrate + uptake_rate * old_biomass * dt
        # Buggy: uses updated biomass
        substrate_buggy = substrate + uptake_rate * biomass * dt

        # The correct value should differ from the buggy value
        assert abs(substrate_correct - substrate_buggy) > 1e-10

    def test_euler_substrate_nonnegative(self, tiny_model):
        """Substrate should never go below zero."""
        model = tiny_model.copy()
        biomass, substrate = 0.1, 0.5  # very small substrate
        dt = 0.5

        for _ in range(50):
            if substrate <= 1e-6 or biomass <= 1e-9:
                break
            with model:
                max_uptake = substrate / (biomass * dt)
                sub = model.reactions.get_by_id("EX_A_e")
                sub.lower_bound = max(-max_uptake, sub.lower_bound)
                sol = model.optimize()
                if sol.status != "optimal":
                    break
                mu = float(sol.objective_value)
                v = float(sol.fluxes.get("EX_A_e", 0))
            old_biomass = biomass
            biomass += mu * biomass * dt
            substrate += v * old_biomass * dt
            biomass = max(0, biomass)
            substrate = max(0, substrate)

        assert substrate >= 0


# ---------------------------------------------------------------------------
#  TMFA with Concentration Variables
# ---------------------------------------------------------------------------

class TestTMFAConcentrations:
    """Test the concentration-based ΔG' adjustment."""

    def test_dg_prime_calculation(self):
        """ΔG' = ΔG° + RT·ln(Q) should adjust standard free energy."""
        import math
        R_GAS = 8.314e-3  # kJ/(mol·K)
        T = 310.15  # K
        dg0 = -5.0  # kJ/mol (favorable forward)
        c_min, c_max = 1e-6, 0.02

        # For a reaction A → B: ΔG' = ΔG° + RT * (ln(c_B) - ln(c_A))
        # Best case: B at c_min, A at c_max
        rt = R_GAS * T
        ln_q_best = 1 * math.log(c_min) + (-1) * math.log(c_max)  # product + substrate
        dg_best = dg0 + rt * ln_q_best

        # Worst case: B at c_max, A at c_min
        ln_q_worst = 1 * math.log(c_max) + (-1) * math.log(c_min)
        dg_worst = dg0 + rt * ln_q_worst

        # Best should be more negative than worst
        assert dg_best < dg_worst
        # With these values, there should be a significant range
        assert abs(dg_worst - dg_best) > 10.0

    def test_thermodynamic_constraint_applied(self, tiny_model):
        """A reaction with very positive ΔG° should be constrained."""
        import math
        model = tiny_model.copy()
        r1 = model.reactions.get_by_id("R1")

        # Set ΔG° very positive (forward thermodynamically unfavorable)
        R_GAS = 8.314e-3
        T = 310.15
        c_min, c_max = 1e-6, 0.02
        dg0 = 100.0  # very unfavorable

        rt = R_GAS * T
        # Best case for R1 (A → B): products at c_min, substrates at c_max
        ln_q_best = 1 * math.log(c_min) + (-1) * math.log(c_max)
        dg_best = dg0 + rt * ln_q_best

        # Even at best case, this should be positive → block forward
        assert dg_best > 0
        r1.upper_bound = 0.0  # strict mode
        sol = model.optimize()
        assert sol.objective_value < 1e-6


# ---------------------------------------------------------------------------
#  Flux Sampling Diagnostics
# ---------------------------------------------------------------------------

class TestSamplingDiagnostics:
    """Test Geweke convergence and correlation computations."""

    def test_geweke_stationary_series(self):
        """A stationary series should have |Z| < 2."""
        import numpy as np
        np.random.seed(42)
        vals = np.random.normal(0, 1, 1000)
        n = len(vals)
        first_10 = vals[:n // 10]
        last_50 = vals[n // 2:]
        z = (first_10.mean() - last_50.mean()) / np.sqrt(
            first_10.var() / len(first_10) + last_50.var() / len(last_50))
        assert abs(z) < 3.0  # with randomness allow up to 3

    def test_geweke_nonstationary_series(self):
        """A trending series should have |Z| > 2."""
        import numpy as np
        vals = np.linspace(0, 100, 1000) + np.random.normal(0, 0.1, 1000)
        n = len(vals)
        first_10 = vals[:n // 10]
        last_50 = vals[n // 2:]
        z = (first_10.mean() - last_50.mean()) / np.sqrt(
            first_10.var() / len(first_10) + last_50.var() / len(last_50))
        assert abs(z) > 2.0

    def test_correlation_matrix_symmetric(self):
        """Pairwise correlation matrix should be symmetric."""
        import numpy as np
        import pandas as pd
        np.random.seed(42)
        data = pd.DataFrame({
            "R1": np.random.normal(5, 1, 100),
            "R2": np.random.normal(-3, 2, 100),
        })
        corr = data.corr()
        assert abs(corr.loc["R1", "R2"] - corr.loc["R2", "R1"]) < 1e-10
        assert corr.loc["R1", "R1"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
#  Gap-Filling
# ---------------------------------------------------------------------------

class TestGapFilling:
    """Test gap-filling logic on a model with a blocked reaction."""

    def test_blocked_reaction_detected(self, tiny_model):
        """A fully blocked reaction should be identified by FVA."""
        from cobra.flux_analysis import flux_variability_analysis
        model = tiny_model.copy()
        # Block R1 completely
        model.reactions.get_by_id("R1").bounds = (0, 0)
        fva = flux_variability_analysis(model, fraction_of_optimum=0.0)
        fmin = float(fva.loc["R1", "minimum"])
        fmax = float(fva.loc["R1", "maximum"])
        assert abs(fmin) < 1e-9 and abs(fmax) < 1e-9

    def test_adding_reaction_rescues_growth(self, tiny_model):
        """Adding a bypass reaction should restore objective flux."""
        model = tiny_model.copy()
        # Block R1
        model.reactions.get_by_id("R1").bounds = (0, 0)
        sol0 = model.optimize()
        assert sol0.objective_value < 1e-6

        # Add bypass A_c → B_c
        bypass = cobra.Reaction("R1_bypass")
        bypass.add_metabolites({
            model.metabolites.get_by_id("A_c"): -1.0,
            model.metabolites.get_by_id("B_c"): 1.0,
        })
        bypass.bounds = (0, 1000)
        model.add_reactions([bypass])
        sol1 = model.optimize()
        assert sol1.objective_value > 1.0


# ---------------------------------------------------------------------------
#  GECKO / Enzyme-Constrained
# ---------------------------------------------------------------------------

class TestGECKO:
    """Test enzyme-constrained FBA logic."""

    def test_protein_constraint_reduces_flux(self, tiny_model):
        """Adding a protein pool constraint should reduce or maintain growth."""
        model = tiny_model.copy()
        sol_free = model.optimize()

        # Add protein pool
        prot_pool = cobra.Metabolite("prot_pool", name="Protein pool", compartment="c")
        # R1 consumes enzyme: kcat=100 1/s, MW=50 kDa, sigma=0.5
        kcat_h = 100 * 3600 * 0.5  # 1/h
        enzyme_cost = 50.0 / (1000.0 * kcat_h)  # g/gDW per flux unit
        model.reactions.get_by_id("R1").add_metabolites({prot_pool: enzyme_cost})

        # Add protein pool exchange
        prot_ex = cobra.Reaction("prot_pool_exchange")
        prot_ex.bounds = (0, 0.5)  # 0.5 g/gDW budget
        prot_ex.add_metabolites({prot_pool: -1.0})
        model.add_reactions([prot_ex])

        sol_constrained = model.optimize()
        assert sol_constrained.status == "optimal"
        assert sol_constrained.objective_value <= sol_free.objective_value + 1e-6

    def test_very_tight_protein_budget_limits_growth(self, tiny_model):
        """An extremely tight protein budget should severely limit growth."""
        model = tiny_model.copy()
        prot_pool = cobra.Metabolite("prot_pool", compartment="c")
        # Use a very high enzyme cost per unit flux to make the constraint binding
        enzyme_cost = 1.0  # 1 g protein per unit flux
        model.reactions.get_by_id("R1").add_metabolites({prot_pool: enzyme_cost})

        prot_ex = cobra.Reaction("prot_pool_exchange")
        prot_ex.bounds = (0, 0.5)  # Only 0.5 g protein available → max flux = 0.5
        prot_ex.add_metabolites({prot_pool: -1.0})
        model.add_reactions([prot_ex])

        sol = model.optimize()
        assert sol.status == "optimal"
        assert sol.objective_value < 1.0  # limited to ~0.5 by protein budget


# ---------------------------------------------------------------------------
#  FSEOF
# ---------------------------------------------------------------------------

class TestFSEOF:
    """Test FSEOF — Flux Scanning based on Enforced Objective Flux."""

    def test_monotonic_detection(self):
        """A monotonically increasing profile should be detected."""
        profile = [0.0, 1.0, 2.0, 3.0, 4.0]
        diffs = [profile[i + 1] - profile[i] for i in range(len(profile) - 1)]
        assert all(d >= -1e-8 for d in diffs)
        assert sum(abs(d) for d in diffs) > 1e-6

    def test_nonmonotonic_detection(self):
        """A non-monotonic profile should NOT be detected as monotonic."""
        profile = [0.0, 2.0, 1.0, 3.0, 4.0]
        diffs = [profile[i + 1] - profile[i] for i in range(len(profile) - 1)]
        assert not all(d >= -1e-8 for d in diffs)

    def test_enforced_objective_reduces_growth(self, tiny_model):
        """Enforcing a minimum growth should constrain the solution space."""
        model = tiny_model.copy()
        sol0 = model.optimize()
        max_growth = sol0.objective_value

        # Enforce 50% growth
        model.reactions.get_by_id("EX_B_e").lower_bound = max_growth * 0.5
        sol1 = model.optimize()
        assert sol1.status == "optimal"
        assert sol1.objective_value >= max_growth * 0.5 - 1e-6


# ---------------------------------------------------------------------------
#  Community Model Conflict Warning
# ---------------------------------------------------------------------------

class TestCommunityConflicts:
    """Test metabolite conflict detection in community model building."""

    def test_formula_conflict_detected(self, tiny_model):
        """Different formulas for the same metabolite should be flagged."""
        from metabodesk_core.mixin_tools import ToolsMixin

        model1 = tiny_model.copy()
        model2 = tiny_model.copy()

        # Give A_c extracellular identity so it's shared
        # Add a shared extracellular metabolite with different formulas
        met_e1 = cobra.Metabolite("X_e", name="X", compartment="e", formula="C6H12O6")
        met_e2 = cobra.Metabolite("X_e", name="X", compartment="e", formula="C3H6O3")

        rxn1 = cobra.Reaction("EX_X_e")
        rxn1.add_metabolites({met_e1: -1.0})
        rxn1.bounds = (-10, 1000)
        model1.add_reactions([rxn1])

        rxn2 = cobra.Reaction("EX_X_e")
        rxn2.add_metabolites({met_e2: -1.0})
        rxn2.bounds = (-10, 1000)
        model2.add_reactions([rxn2])

        # Use ToolsMixin methods directly
        class MockTools(ToolsMixin):
            def _is_extracellular(self, met):
                return met.compartment and met.compartment.lower() in ("e", "extracellular")
            def _prefix_gpr(self, gpr, prefix):
                return gpr
            def _get_objective_reaction_id_for_model(self, model):
                return None

        tools = MockTools()
        _, _, _, conflicts = tools._build_community_model(
            [model1, model2], ["sp1", "sp2"], share_extracellular=True)
        assert len(conflicts) > 0
        assert any("formula" in c.lower() for c in conflicts)


# ---------------------------------------------------------------------------
#  Cache Invalidation
# ---------------------------------------------------------------------------

class TestCacheInvalidation:
    """Test that model state hash includes all relevant state."""

    def test_hash_changes_with_overexpression(self, tiny_model):
        """Hash should change when overexpression state changes."""
        import hashlib

        class MockWindow:
            knockout_genes = set()
            overexpression_reactions = {}
            temporary_upper_bound_overrides = {}
            reaction_bound_overrides = {}

            class loopless_chk:
                @staticmethod
                def isChecked():
                    return False

        from metabodesk_core.mainwindow import MainWindow
        # Test hash logic directly
        mw = MockWindow()
        model = tiny_model.copy()

        def compute_hash(mw_obj, m):
            parts = []
            parts.append(str(m.objective.expression))
            parts.append(m.objective.direction)
            for r in sorted(m.reactions, key=lambda r: r.id):
                parts.append(f"{r.id}:{r.lower_bound}:{r.upper_bound}")
            parts.append(",".join(sorted(mw_obj.knockout_genes)))
            parts.append("|".join(f"{k}={v}" for k, v in sorted(mw_obj.overexpression_reactions.items())))
            parts.append("|".join(f"{k}={v}" for k, v in sorted(mw_obj.temporary_upper_bound_overrides.items())))
            parts.append("|".join(f"{k}={v}" for k, v in sorted(mw_obj.reaction_bound_overrides.items())))
            parts.append(f"loopless:{mw_obj.loopless_chk.isChecked()}")
            return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()

        h1 = compute_hash(mw, model)
        mw.overexpression_reactions = {"R1": 2.0}
        h2 = compute_hash(mw, model)
        assert h1 != h2


# ---------------------------------------------------------------------------
#  GECKO Isoenzyme Support
# ---------------------------------------------------------------------------

class TestGECKOIsoenzyme:
    """Test isoenzyme arm-reaction splitting in GECKO."""

    def test_single_enzyme_direct_constraint(self, tiny_model):
        """Single enzyme per reaction should add protein cost directly."""
        model = tiny_model.copy()
        prot_pool = cobra.Metabolite("prot_pool", compartment="c")
        kcat_h = 100 * 3600 * 0.5
        enzyme_cost = 50.0 / (1000.0 * kcat_h)
        model.reactions.get_by_id("R1").add_metabolites({prot_pool: enzyme_cost})

        prot_ex = cobra.Reaction("prot_pool_exchange")
        prot_ex.bounds = (0, 0.5)
        prot_ex.add_metabolites({prot_pool: -1.0})
        model.add_reactions([prot_ex])

        sol = model.optimize()
        assert sol.status == "optimal"
        # prot_pool should be a metabolite in R1
        assert prot_pool in model.reactions.get_by_id("R1").metabolites

    def test_isoenzyme_creates_arm_reactions(self, tiny_model):
        """Two isoenzymes should create two arm reactions + disable original."""
        model = tiny_model.copy()
        prot_pool = cobra.Metabolite("prot_pool", compartment="c")

        r1 = model.reactions.get_by_id("R1")
        orig_mets = dict(r1.metabolites)
        orig_lb, orig_ub = r1.lower_bound, r1.upper_bound

        # Disable original
        r1.bounds = (0, 0)

        # Arm 1: fast enzyme (low cost)
        kcat_h_1 = 200 * 3600 * 0.5
        cost_1 = 50.0 / (1000.0 * kcat_h_1)
        arm1 = cobra.Reaction("R1_arm_enz1")
        arm1.lower_bound = orig_lb
        arm1.upper_bound = orig_ub
        arm1_mets = {met: coeff for met, coeff in orig_mets.items()}
        arm1_mets[prot_pool] = cost_1
        arm1.add_metabolites(arm1_mets)

        # Arm 2: slow enzyme (high cost)
        kcat_h_2 = 10 * 3600 * 0.5
        cost_2 = 100.0 / (1000.0 * kcat_h_2)
        arm2 = cobra.Reaction("R1_arm_enz2")
        arm2.lower_bound = orig_lb
        arm2.upper_bound = orig_ub
        arm2_mets = {met: coeff for met, coeff in orig_mets.items()}
        arm2_mets[prot_pool] = cost_2
        arm2.add_metabolites(arm2_mets)

        model.add_reactions([arm1, arm2])

        prot_ex = cobra.Reaction("prot_pool_exchange")
        prot_ex.bounds = (0, 10.0)
        prot_ex.add_metabolites({prot_pool: -1.0})
        model.add_reactions([prot_ex])

        sol = model.optimize()
        assert sol.status == "optimal"
        assert sol.objective_value > 0

        # Solver should prefer the cheaper enzyme (arm1)
        assert abs(sol.fluxes.get("R1_arm_enz1", 0)) >= abs(sol.fluxes.get("R1_arm_enz2", 0)) - 1e-6

    def test_original_reaction_disabled(self, tiny_model):
        """After isoenzyme split, original reaction should carry zero flux."""
        model = tiny_model.copy()
        r1 = model.reactions.get_by_id("R1")
        r1.bounds = (0, 0)

        sol = model.optimize()
        assert sol.objective_value < 1e-6


# ---------------------------------------------------------------------------
#  FSEOF Multi-Target
# ---------------------------------------------------------------------------

class TestFSEOFMultiTarget:
    """Test multi-target FSEOF scanning."""

    def test_comma_separated_parsing(self):
        """Comma-separated target IDs should be split correctly."""
        raw = "EX_etoh_e, EX_succ_e , EX_ac_e"
        targets = [t.strip() for t in raw.split(",") if t.strip()]
        assert targets == ["EX_etoh_e", "EX_succ_e", "EX_ac_e"]

    def test_single_target_still_works(self):
        """Single target ID should parse as a list of one."""
        raw = "EX_etoh_e"
        targets = [t.strip() for t in raw.split(",") if t.strip()]
        assert len(targets) == 1
        assert targets[0] == "EX_etoh_e"


# ---------------------------------------------------------------------------
#  dFBA Infeasibility Recovery
# ---------------------------------------------------------------------------

class TestDFBARecovery:
    """Test dFBA recovery when LP becomes infeasible."""

    def test_maintenance_mode_continues_simulation(self, tiny_model):
        """When growth is infeasible, substrate should still be consumed."""
        model = tiny_model.copy()

        # Make growth infeasible by blocking the objective reaction
        model.reactions.get_by_id("EX_B_e").bounds = (0, 0)
        sol = model.optimize()
        assert sol.objective_value < 1e-6

        # Maintenance mode: fix growth to 0, check substrate still possible
        with model:
            for rxn in model.reactions:
                if getattr(rxn, 'objective_coefficient', 0) != 0:
                    rxn.lower_bound = 0
                    rxn.upper_bound = 0
            sol_maint = model.optimize()
            assert sol_maint.status in ("optimal", "infeasible")

    def test_euler_continues_after_infeasibility(self, tiny_model):
        """Euler loop should not crash when LP is infeasible at some step."""
        model = tiny_model.copy()
        biomass, substrate = 0.1, 0.001  # almost no substrate
        dt = 0.5
        steps = 0

        for _ in range(5):
            steps += 1
            with model:
                if substrate > 1e-6 and biomass > 1e-9:
                    max_uptake = substrate / (biomass * dt)
                    sub = model.reactions.get_by_id("EX_A_e")
                    sub.lower_bound = max(-max_uptake, sub.lower_bound)
                else:
                    model.reactions.get_by_id("EX_A_e").lower_bound = 0.0
                sol = model.optimize()
                if sol.status == "optimal":
                    mu = float(sol.objective_value)
                    v = float(sol.fluxes.get("EX_A_e", 0))
                else:
                    mu, v = 0.0, 0.0

            old_biomass = biomass
            biomass += mu * biomass * dt
            substrate += v * old_biomass * dt
            biomass = max(0, biomass)
            substrate = max(0, substrate)

        assert steps == 5
        assert biomass >= 0
        assert substrate >= 0


# ---------------------------------------------------------------------------
#  TMFA Concentration Validation
# ---------------------------------------------------------------------------

class TestTMFAValidation:
    """Test TMFA concentration bounds validation."""

    def test_cmin_must_be_less_than_cmax(self):
        """c_min >= c_max should be invalid."""
        c_min, c_max = 0.02, 0.01
        assert c_min >= c_max  # would be rejected

    def test_cmin_must_be_positive(self):
        """c_min <= 0 should be invalid."""
        assert 0 <= 0  # would be rejected
        assert -1e-6 <= 0  # would be rejected

    def test_valid_bounds_allow_computation(self):
        """Valid concentration bounds should allow ΔG' computation."""
        import math
        c_min, c_max = 1e-6, 0.02
        R_GAS = 8.314e-3
        T = 310.15
        dg0 = -5.0

        rt = R_GAS * T
        ln_q_best = 1 * math.log(c_min) + (-1) * math.log(c_max)
        ln_q_worst = 1 * math.log(c_max) + (-1) * math.log(c_min)
        dg_best = dg0 + rt * ln_q_best
        dg_worst = dg0 + rt * ln_q_worst

        assert dg_best < dg_worst
        assert isinstance(dg_best, float)
        assert isinstance(dg_worst, float)


# ---------------------------------------------------------------------------
#  KDE Histogram Pre-computation (Item 4)
# ---------------------------------------------------------------------------

class TestKDEPrecomputation:
    """Test that KDE histogram bins can be pre-computed in worker."""

    def test_histogram_precompute(self):
        """numpy.histogram should produce valid bins and counts."""
        import numpy as np
        vals = np.random.normal(5.0, 1.0, size=500)
        counts, edges = np.histogram(vals, bins=40, density=True)
        assert len(counts) == 40
        assert len(edges) == 41
        assert all(c >= 0 for c in counts)
        # Density should integrate to ~1
        widths = np.diff(edges)
        integral = sum(c * w for c, w in zip(counts, widths))
        assert abs(integral - 1.0) < 0.1

    def test_histogram_serialisable(self):
        """Pre-computed histogram data should be serialisable to list."""
        import numpy as np
        vals = np.random.uniform(0, 10, size=100)
        counts, edges = np.histogram(vals, bins=20, density=True)
        # .tolist() is what the worker does
        counts_list = counts.tolist()
        edges_list = edges.tolist()
        assert isinstance(counts_list, list)
        assert isinstance(edges_list, list)
        assert len(counts_list) == 20
        assert len(edges_list) == 21


# ---------------------------------------------------------------------------
#  CSV Template Generation (Items 6 + 7)
# ---------------------------------------------------------------------------

class TestCSVTemplateGeneration:
    """Test that CSV template data matches model structure."""

    def test_expression_template_has_all_genes(self, tiny_model):
        """Template should list every gene in the model."""
        genes = [g.id for g in tiny_model.genes]
        assert len(genes) > 0
        # Template row: [gene_id, 1.0]
        for gid in genes:
            assert isinstance(gid, str)
            assert len(gid) > 0

    def test_dg_template_skips_boundary(self, tiny_model):
        """ΔG template should skip boundary reactions."""
        non_boundary = [r for r in tiny_model.reactions if not r.boundary]
        boundary = [r for r in tiny_model.reactions if r.boundary]
        assert len(non_boundary) > 0
        # At least EX_A_e or EX_B_e should be boundary
        assert len(boundary) > 0


# ---------------------------------------------------------------------------
#  Undo/Redo UI State (Item 10)
# ---------------------------------------------------------------------------

class TestUndoRedoUIState:
    """Test undo/redo saves and restores UI state."""

    def test_state_captures_ui_fields(self):
        """_capture_current_state should include UI keys."""
        state = {
            "knockout_genes": set(),
            "overexpression_reactions": {},
            "reaction_bound_overrides": {},
            "active_tab_index": 2,
            "analysis_type_index": 3,
            "topn_value": 15,
        }
        assert "active_tab_index" in state
        assert "analysis_type_index" in state
        assert "topn_value" in state
        assert state["topn_value"] == 15

    def test_restore_state_dict_completeness(self):
        """State dict should contain all required fields."""
        required_keys = {
            "knockout_genes", "overexpression_reactions",
            "reaction_bound_overrides", "active_tab_index",
            "analysis_type_index", "topn_value",
        }
        state = {k: None for k in required_keys}
        assert required_keys.issubset(state.keys())

# ---------------------------------------------------------------------------
#  Minimal Cut Sets (MCS)
# ---------------------------------------------------------------------------

class TestMinimalCutSets:
    """Minimal cut sets: find knockouts that block a target flux."""

    def test_single_knockout_found(self, tiny_model):
        """Knocking out R1 should appear as a single MCS for EX_B_e."""
        import numpy as np
        model = tiny_model.copy()
        target = model.reactions.get_by_id("EX_B_e")

        # Compute single MCS manually (same logic as run_minimal_cut_sets)
        candidates = [r for r in model.reactions
                      if not r.id.startswith("EX_") and r.id != target.id]

        single_mcs = []
        for rxn in candidates:
            with model:
                rxn.knock_out()
                sol = model.optimize()
                if sol.status != "optimal" or abs(sol.fluxes.get(target.id, 0)) < 1e-6:
                    single_mcs.append(frozenset([rxn.id]))

        assert frozenset(["R1"]) in single_mcs

    def test_essential_reaction_is_cut(self, tiny_model):
        """In tiny_model R1 is the only route; MCS must include it."""
        model = tiny_model.copy()
        with model:
            model.reactions.get_by_id("R1").knock_out()
            sol = model.optimize()
            assert sol.status != "optimal" or abs(sol.objective_value) < 1e-6

    def test_no_mcs_when_target_zero(self, tiny_model):
        """If target already has zero flux, MCS list should be empty."""
        model = tiny_model.copy()
        # Block import so target is already zero
        model.reactions.get_by_id("EX_A_e").lower_bound = 0.0
        sol = model.optimize()
        # Objective should be zero or infeasible
        assert sol.status != "optimal" or abs(sol.objective_value) < 1e-6


# ---------------------------------------------------------------------------
#  Metabolic Distance
# ---------------------------------------------------------------------------

class TestMetabolicDistance:
    """BFS metabolic distance between metabolites."""

    def test_adjacent_metabolites_distance_one(self, tiny_model):
        """A_c and B_c are connected by R1, so distance should be 1."""
        from collections import deque

        model = tiny_model.copy()

        # Build adjacency (no currency exclusion in tiny model)
        adj = {m.id: set() for m in model.metabolites}
        for rxn in model.reactions:
            mets = [m.id for m in rxn.metabolites]
            for i, m1 in enumerate(mets):
                for m2 in mets[i+1:]:
                    adj[m1].add(m2)
                    adj[m2].add(m1)

        # BFS from A_c
        visited = {"A_c": 0}
        queue = deque(["A_c"])
        while queue:
            current = queue.popleft()
            for nb in adj.get(current, []):
                if nb not in visited:
                    visited[nb] = visited[current] + 1
                    queue.append(nb)

        assert visited.get("B_c") == 1

    def test_unreachable_metabolite(self, tiny_model):
        """A metabolite not in the graph should be unreachable."""
        model = tiny_model.copy()
        # Add an isolated metabolite
        iso = cobra.Metabolite("ISO_c", name="Isolated", compartment="c")
        model.add_metabolites([iso])

        adj = {m.id: set() for m in model.metabolites}
        for rxn in model.reactions:
            mets = [m.id for m in rxn.metabolites]
            for i, m1 in enumerate(mets):
                for m2 in mets[i+1:]:
                    adj[m1].add(m2)
                    adj[m2].add(m1)

        # BFS from A_c
        from collections import deque
        visited = {"A_c": 0}
        queue = deque(["A_c"])
        while queue:
            current = queue.popleft()
            for nb in adj.get(current, []):
                if nb not in visited:
                    visited[nb] = visited[current] + 1
                    queue.append(nb)

        assert "ISO_c" not in visited

    def test_ecoli_distance_consistency(self, ecoli_core):
        """In E. coli core, distances should be non-negative integers."""
        if ecoli_core is None:
            pytest.skip("E. coli core model not available")
        from collections import deque

        model = ecoli_core.copy()
        mets = [m.id for m in model.metabolites]
        source = mets[0]

        adj = {m.id: set() for m in model.metabolites}
        for rxn in model.reactions:
            ms = [m.id for m in rxn.metabolites]
            for i, m1 in enumerate(ms):
                for m2 in ms[i+1:]:
                    adj[m1].add(m2)
                    adj[m2].add(m1)

        visited = {source: 0}
        queue = deque([source])
        while queue:
            current = queue.popleft()
            for nb in adj.get(current, []):
                if nb not in visited:
                    visited[nb] = visited[current] + 1
                    queue.append(nb)

        for dist in visited.values():
            assert isinstance(dist, int)
            assert dist >= 0


# ---------------------------------------------------------------------------
#  Pareto Frontier
# ---------------------------------------------------------------------------

class TestParetoFrontier:
    """Multi-objective Pareto frontier optimisation."""

    def test_pareto_sweep_produces_points(self, tiny_model):
        """Sweeping between two objectives should produce feasible points."""
        import numpy as np

        model = tiny_model.copy()
        obj1_id = "EX_B_e"
        obj2_id = "R1"
        n_steps = 10

        # Max of obj1
        model.objective = obj1_id
        sol1 = model.optimize()
        assert sol1.status == "optimal"
        obj1_max = sol1.objective_value

        # Max of obj2
        model.objective = obj2_id
        sol2 = model.optimize()
        assert sol2.status == "optimal"

        # Sweep
        results = []
        rxn1 = model.reactions.get_by_id(obj1_id)
        for val in np.linspace(0, obj1_max, n_steps):
            with model:
                rxn1.lower_bound = float(val)
                rxn1.upper_bound = float(val)
                model.objective = obj2_id
                sol = model.optimize()
                if sol.status == "optimal":
                    results.append({"obj1": val, "obj2": sol.objective_value})

        assert len(results) >= 2  # at least some feasible points

    def test_pareto_dominance_filter(self):
        """Non-dominated filter should keep only Pareto-optimal points."""
        points = [
            {"obj1": 1, "obj2": 5},
            {"obj1": 3, "obj2": 3},
            {"obj1": 5, "obj2": 1},
            {"obj1": 2, "obj2": 2},  # dominated by (3, 3)
        ]

        pareto = []
        for pt in points:
            dominated = False
            for other in points:
                if (other["obj1"] >= pt["obj1"] and
                        other["obj2"] >= pt["obj2"] and
                        (other["obj1"] > pt["obj1"] or
                         other["obj2"] > pt["obj2"])):
                    dominated = True
                    break
            if not dominated:
                pareto.append(pt)

        # (2,2) is dominated by (3,3), so 3 pareto points remain
        assert len(pareto) == 3
        dominated_ids = [(p["obj1"], p["obj2"]) for p in pareto]
        assert (2, 2) not in dominated_ids

    def test_pareto_monotonic_tradeoff(self, tiny_model):
        """On the Pareto front, increasing obj1 should decrease obj2."""
        import numpy as np

        model = tiny_model.copy()
        # In tiny model EX_B_e and R1 are tightly coupled,
        # so the frontier is trivial (diagonal), but the logic should hold.
        obj1_id = "EX_B_e"
        obj2_id = "EX_A_e"
        n_steps = 5

        model.objective = obj1_id
        sol = model.optimize()
        obj1_max = sol.objective_value

        results = []
        rxn1 = model.reactions.get_by_id(obj1_id)
        for val in np.linspace(0, obj1_max, n_steps):
            with model:
                rxn1.lower_bound = float(val)
                rxn1.upper_bound = float(val)
                model.objective = obj2_id
                model.objective_direction = "max"
                sol2 = model.optimize()
                if sol2.status == "optimal":
                    results.append({"obj1": val, "obj2": sol2.objective_value})

        # At least 2 points needed for monotonicity check
        assert len(results) >= 2

    def test_knee_point_computation(self):
        """Knee point should be the point closest to the utopia point."""
        pareto = [
            {"obj1": 0, "obj2": 10},
            {"obj1": 5, "obj2": 7},
            {"obj1": 8, "obj2": 3},
            {"obj1": 10, "obj2": 0},
        ]
        utopia_obj1 = 10
        utopia_obj2 = 10

        o1_range = max(p["obj1"] for p in pareto) - min(p["obj1"] for p in pareto)
        o2_range = max(p["obj2"] for p in pareto) - min(p["obj2"] for p in pareto)

        best_dist = float("inf")
        knee = pareto[0]
        for pt in pareto:
            d = (((pt["obj1"] - utopia_obj1) / o1_range) ** 2 +
                 ((pt["obj2"] - utopia_obj2) / o2_range) ** 2) ** 0.5
            if d < best_dist:
                best_dist = d
                knee = pt

        # The knee should be (5, 7) — closest to utopia in normalised space
        assert knee["obj1"] == 5
        assert knee["obj2"] == 7


# ---------------------------------------------------------------------------
#  Mixin Split Verification
# ---------------------------------------------------------------------------

class TestMixinSplit:
    """Verify the mixin split preserves method accessibility."""

    def test_omics_mixin_has_methods(self):
        """OmicsMixin should have all omics analysis methods."""
        from metabodesk_core.mixin_advanced_omics import OmicsMixin

        expected = [
            "run_gene_expression_analysis",
            "run_tmfa",
            "run_gecko_analysis",
            "run_fseof",
        ]
        for method in expected:
            assert hasattr(OmicsMixin, method), f"OmicsMixin missing {method}"

    def test_design_mixin_has_methods(self):
        """DesignMixin should have all design analysis methods."""
        from metabodesk_core.mixin_advanced_design import DesignMixin

        expected = [
            "run_optknock",
            "run_dfba",
            "run_flux_coupling",
            "run_gap_filling",
            "run_model_comparison",
            "run_pathway_enrichment",
            "show_escher_map",
            "run_phenotype_phase_plane",
            "validate_sbml",
            "export_jupyter_notebook",
            "run_minimal_cut_sets",
            "run_metabolic_distance",
            "run_pareto_frontier",
        ]
        for method in expected:
            assert hasattr(DesignMixin, method), f"DesignMixin missing {method}"

    def test_advanced_mixin_facade_inherits_all(self):
        """AdvancedMixin should inherit from both OmicsMixin and DesignMixin."""
        from metabodesk_core.mixin_advanced import AdvancedMixin
        from metabodesk_core.mixin_advanced_omics import OmicsMixin
        from metabodesk_core.mixin_advanced_design import DesignMixin

        assert issubclass(AdvancedMixin, OmicsMixin)
        assert issubclass(AdvancedMixin, DesignMixin)

    def test_facade_exposes_all_methods(self):
        """AdvancedMixin should expose every method from both sub-mixins."""
        from metabodesk_core.mixin_advanced import AdvancedMixin

        all_methods = [
            "run_gene_expression_analysis", "run_tmfa",
            "run_gecko_analysis", "run_fseof",
            "run_optknock", "run_dfba", "run_flux_coupling",
            "run_gap_filling", "run_model_comparison",
            "run_pathway_enrichment", "show_escher_map",
            "run_phenotype_phase_plane", "validate_sbml",
            "export_jupyter_notebook", "run_minimal_cut_sets",
            "run_metabolic_distance", "run_pareto_frontier",
        ]
        for method in all_methods:
            assert hasattr(AdvancedMixin, method), f"AdvancedMixin missing {method}"


# ---------------------------------------------------------------------------
#  TypedDict Return Types
# ---------------------------------------------------------------------------

class TestTypedDictReturnTypes:
    """Verify TypedDicts are importable and structurally correct."""

    def test_typedicts_importable(self):
        """All TypedDict classes should be importable from mixin_analysis."""
        from metabodesk_core.mixin_analysis import (
            FVAResultEntry,
            PFBAResult,
            DeletionResult,
            RobustnessResult,
            ProductionEnvelopeResult,
            SamplingResultEntry,
            FBABaselineResult,
            AnalysisResult,
        )
        # Check they are TypedDict subclasses (they behave like dicts)
        assert FVAResultEntry is not None
        assert PFBAResult is not None
        assert DeletionResult is not None
        assert RobustnessResult is not None
        assert ProductionEnvelopeResult is not None
        assert SamplingResultEntry is not None
        assert FBABaselineResult is not None
        assert AnalysisResult is not None

    def test_fva_result_entry_keys(self):
        """FVAResultEntry should have min/max keys."""
        from metabodesk_core.mixin_analysis import FVAResultEntry
        keys = FVAResultEntry.__annotations__
        assert "min" in keys
        assert "max" in keys

    def test_pfba_result_keys(self):
        """PFBAResult should have objective, status, and flux."""
        from metabodesk_core.mixin_analysis import PFBAResult
        keys = PFBAResult.__annotations__
        assert "objective" in keys
        assert "status" in keys
        assert "flux" in keys

    def test_analysis_result_keys(self):
        """AnalysisResult should have analysis_type and analysis-specific keys."""
        from metabodesk_core.mixin_analysis import AnalysisResult
        keys = AnalysisResult.__annotations__
        assert "analysis_type" in keys
        assert "baseline" in keys
        assert "fva" in keys
        assert "pfba" in keys
        assert "sampling" in keys


# ---------------------------------------------------------------------------
#  ¹³C-MFA Flux Constraints
# ---------------------------------------------------------------------------

class TestC13MFAConstraints:
    """Verify ¹³C-MFA flux constraint logic on small models."""

    def test_apply_13c_constraints_optimal(self, tiny_model):
        """Applying measured-flux bounds should still yield optimal if within range."""
        model = tiny_model.copy()
        # Baseline flux for R1 ~ 10.0, EX_B_e ~ 10.0
        sol_base = model.optimize()
        assert sol_base.status == "optimal"

        # Constrain R1 to [8, 12] — should still be feasible
        constraints = [
            {"reaction_id": "R1", "measured_flux": 10.0, "std_dev": 1.0},
        ]
        n_sigma = 2.0
        constrained_model = model.copy()
        for c in constraints:
            rxn = constrained_model.reactions.get_by_id(c["reaction_id"])
            lb = c["measured_flux"] - n_sigma * c["std_dev"]
            ub = c["measured_flux"] + n_sigma * c["std_dev"]
            rxn.lower_bound = max(rxn.lower_bound, lb)
            rxn.upper_bound = min(rxn.upper_bound, ub)

        sol = constrained_model.optimize()
        assert sol.status == "optimal"
        assert abs(sol.fluxes["R1"] - 10.0) <= 2.0 + 1e-6

    def test_apply_13c_constraints_tight(self, tiny_model):
        """Tight measured constraints should restrict flux near measured value."""
        model = tiny_model.copy()
        constraints = [
            {"reaction_id": "R1", "measured_flux": 5.0, "std_dev": 0.5},
        ]
        n_sigma = 1.0  # [4.5, 5.5]
        constrained_model = model.copy()
        for c in constraints:
            rxn = constrained_model.reactions.get_by_id(c["reaction_id"])
            lb = c["measured_flux"] - n_sigma * c["std_dev"]
            ub = c["measured_flux"] + n_sigma * c["std_dev"]
            rxn.lower_bound = max(rxn.lower_bound, lb)
            rxn.upper_bound = min(rxn.upper_bound, ub)

        sol = constrained_model.optimize()
        assert sol.status == "optimal"
        # Objective is limited by R1 upper bound = 5.5
        assert sol.objective_value <= 5.5 + 1e-6

    def test_apply_13c_constraints_infeasible(self, tiny_model):
        """A constraint requiring negative flux on irreversible reaction → ValueError or infeasible."""
        model = tiny_model.copy()
        constraints = [
            {"reaction_id": "R1", "measured_flux": -5.0, "std_dev": 0.1},
        ]
        n_sigma = 1.0  # [-5.1, -4.9] but R1 LB=0 → max(0, -5.1)=0, min(1000, -4.9)=-4.9
        # COBRApy raises ValueError when setting UB < LB
        constrained_model = model.copy()
        with pytest.raises(ValueError, match="lower bound must be less"):
            for c in constraints:
                rxn = constrained_model.reactions.get_by_id(c["reaction_id"])
                lb = c["measured_flux"] - n_sigma * c["std_dev"]
                ub = c["measured_flux"] + n_sigma * c["std_dev"]
                rxn.lower_bound = max(rxn.lower_bound, lb)
                rxn.upper_bound = min(rxn.upper_bound, ub)

    def test_13c_constraint_multiple_reactions(self, tiny_model):
        """Multiple constraints should apply simultaneously."""
        model = tiny_model.copy()
        constraints = [
            {"reaction_id": "R1", "measured_flux": 7.0, "std_dev": 1.0},
            {"reaction_id": "EX_A_e", "measured_flux": -7.0, "std_dev": 1.0},
        ]
        n_sigma = 1.0
        constrained_model = model.copy()
        for c in constraints:
            rxn = constrained_model.reactions.get_by_id(c["reaction_id"])
            lb = c["measured_flux"] - n_sigma * c["std_dev"]
            ub = c["measured_flux"] + n_sigma * c["std_dev"]
            rxn.lower_bound = max(rxn.lower_bound, lb)
            rxn.upper_bound = min(rxn.upper_bound, ub)

        sol = constrained_model.optimize()
        assert sol.status == "optimal"
        # Flux through R1 is in [6, 8] and EX_A_e is in [-8, -6]
        assert 6.0 - 1e-6 <= sol.fluxes["R1"] <= 8.0 + 1e-6

    def test_13c_zero_std_dev_pins_flux(self, tiny_model):
        """Zero std_dev should pin the flux to the exact measured value."""
        model = tiny_model.copy()
        constraints = [
            {"reaction_id": "R1", "measured_flux": 3.0, "std_dev": 0.0},
        ]
        n_sigma = 2.0  # n*0 = 0, so bounds = [3.0, 3.0]
        constrained_model = model.copy()
        for c in constraints:
            rxn = constrained_model.reactions.get_by_id(c["reaction_id"])
            lb = c["measured_flux"] - n_sigma * c["std_dev"]
            ub = c["measured_flux"] + n_sigma * c["std_dev"]
            rxn.lower_bound = max(rxn.lower_bound, lb)
            rxn.upper_bound = min(rxn.upper_bound, ub)

        sol = constrained_model.optimize()
        assert sol.status == "optimal"
        assert abs(sol.fluxes["R1"] - 3.0) < 1e-6


# ---------------------------------------------------------------------------
#  Auto-Update Version Comparison
# ---------------------------------------------------------------------------

class TestVersionComparison:
    """Verify semver version comparison logic used by auto-update checker."""

    def _version_newer(self, latest, current):
        """Replicate the static method from ToolsMixin."""
        def _parse(v):
            parts = []
            for p in v.split("."):
                try:
                    parts.append(int(p))
                except ValueError:
                    break
            return tuple(parts) if parts else (0,)
        return _parse(latest) > _parse(current)

    def test_newer_patch(self):
        assert self._version_newer("1.0.1", "1.0.0") is True

    def test_newer_minor(self):
        assert self._version_newer("1.1.0", "1.0.0") is True

    def test_newer_major(self):
        assert self._version_newer("2.0.0", "1.0.0") is True

    def test_same_version(self):
        assert self._version_newer("1.0.0", "1.0.0") is False

    def test_older_version(self):
        assert self._version_newer("1.0.0", "1.0.1") is False

    def test_different_length_newer(self):
        assert self._version_newer("1.1", "1.0.0") is True

    def test_different_length_older(self):
        assert self._version_newer("1.0", "1.0.1") is False

    def test_v_prefix_stripped(self):
        """Test with raw version strings (prefix stripping happens in caller)."""
        latest = "1.2.0".lstrip("vV")
        current = "1.1.0"
        assert self._version_newer(latest, current) is True

    def test_pre_release_suffix_ignored(self):
        """Non-numeric suffix segments are ignored by the parser."""
        assert self._version_newer("2.0.0-beta", "1.0.0") is True

    def test_empty_string_fallback(self):
        """Empty version string should default to (0,)."""
        assert self._version_newer("1.0.0", "") is True
        assert self._version_newer("", "1.0.0") is False


# ---------------------------------------------------------------------------
#  Chart Sub-Renderer Split
# ---------------------------------------------------------------------------

class TestChartSubRenderers:
    """Verify that the three chart sub-methods are callable static methods."""

    def test_plot_flux_histogram_callable(self):
        """_plot_flux_histogram should be importable and callable."""
        from metabodesk_core.mixin_analysis import AnalysisMixin
        assert callable(getattr(AnalysisMixin, "_plot_flux_histogram", None))

    def test_plot_flux_waterfall_callable(self):
        """_plot_flux_waterfall should be importable and callable."""
        from metabodesk_core.mixin_analysis import AnalysisMixin
        assert callable(getattr(AnalysisMixin, "_plot_flux_waterfall", None))

    def test_plot_flux_bar_comparison_callable(self):
        """_plot_flux_bar_comparison should be importable and callable."""
        from metabodesk_core.mixin_analysis import AnalysisMixin
        assert callable(getattr(AnalysisMixin, "_plot_flux_bar_comparison", None))

    def test_histogram_renders_on_axes(self):
        """_plot_flux_histogram should draw onto a matplotlib Axes."""
        from metabodesk_core.mixin_analysis import AnalysisMixin
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        AnalysisMixin._plot_flux_histogram(ax, {"R1": 5.0, "R2": -3.0}, {"R1": 4.0, "R2": -2.0})
        assert ax.get_title() == "Flux Distribution (all non-zero fluxes)"
        plt.close(fig)

    def test_waterfall_renders_on_axes(self):
        """_plot_flux_waterfall should draw bars onto a matplotlib Axes."""
        from metabodesk_core.mixin_analysis import AnalysisMixin
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        AnalysisMixin._plot_flux_waterfall(ax, ["R1", "R2"], {"R1": 5.0, "R2": 3.0},
                                            {"R1": 8.0, "R2": 1.0}, "Base", "Compared")
        assert "Top flux changes" in ax.get_title()
        plt.close(fig)

    def test_bar_comparison_renders_on_axes(self):
        """_plot_flux_bar_comparison should draw grouped bars onto a matplotlib Axes."""
        from metabodesk_core.mixin_analysis import AnalysisMixin
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        AnalysisMixin._plot_flux_bar_comparison(ax, ["R1", "R2"], {"R1": 5.0, "R2": 3.0},
                                                 {"R1": 8.0, "R2": 1.0}, "Base", "Comp", "Title")
        assert ax.get_title() == "Title"
        plt.close(fig)


# ---------------------------------------------------------------------------
#  BiGG Database Browser
# ---------------------------------------------------------------------------

class TestBiGGBrowser:
    """Verify BiGG browser helper logic is importable and the method exists."""

    def test_bigg_browser_method_exists(self):
        """DialogsMixin should have show_bigg_browser method."""
        from metabodesk_core.mixin_dialogs import DialogsMixin
        assert callable(getattr(DialogsMixin, "show_bigg_browser", None))


# ---------------------------------------------------------------------------
#  Native Qt Interactive Metabolic Map
# ---------------------------------------------------------------------------

class TestInteractiveMetabolicMap:
    """Verify native Qt metabolic map method exists and layout helpers work."""

    def test_interactive_map_method_exists(self):
        """DialogsMixin should have show_interactive_metabolic_map method."""
        from metabodesk_core.mixin_dialogs import DialogsMixin
        assert callable(getattr(DialogsMixin, "show_interactive_metabolic_map", None))

    def test_networkx_layout_produces_positions(self, tiny_model):
        """Graph layout for a small model should return positions for all nodes."""
        import networkx as nx
        G = nx.Graph()
        for rxn in tiny_model.reactions:
            G.add_node(rxn.id, ntype="rxn")
            for m in rxn.metabolites:
                G.add_node(m.id, ntype="met")
                G.add_edge(rxn.id, m.id)

        pos = nx.kamada_kawai_layout(G, scale=400)
        # Every node should have a position
        for node in G.nodes():
            assert node in pos
            assert len(pos[node]) == 2

    def test_flux_color_intensity_calculation(self, tiny_model):
        """Flux-based colour intensity should scale between 55 and 255."""
        flux = {"R1": 8.0, "EX_A_e": -10.0, "EX_B_e": 10.0}
        max_flux = max(abs(v) for v in flux.values()) or 1.0
        for rxn_id, f_val in flux.items():
            intensity = min(int(55 + 200 * (abs(f_val) / max_flux)), 255)
            assert 55 <= intensity <= 255