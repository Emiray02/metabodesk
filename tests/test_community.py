"""Tests for community model building and GPR prefixing.

Uses a minimal harness that instantiates :class:`ToolsMixin` directly
(no Qt GUI required) to test ``_build_community_model``,
``_prefix_gpr``, and ``_is_extracellular``.
"""

import pytest

cobra = pytest.importorskip("cobra")


# ---------------------------------------------------------------------------
#  Minimal harness – import ToolsMixin at test time so PySide6 is only
#  needed when the test actually runs.
# ---------------------------------------------------------------------------

@pytest.fixture
def tools():
    """Return a bare ToolsMixin instance (no GUI state)."""
    from metabodesk_core.mixin_tools import ToolsMixin

    class _Harness(ToolsMixin):
        """Minimal instance with only the methods under test."""
        pass

    return _Harness()


# ---------------------------------------------------------------------------
#  _prefix_gpr
# ---------------------------------------------------------------------------

class TestPrefixGPR:
    """GPR string prefixing for community model construction."""

    def test_simple_and(self, tools):
        assert tools._prefix_gpr("g1 and g2", "sp__") == "sp__g1 and sp__g2"

    def test_simple_or(self, tools):
        assert tools._prefix_gpr("g1 or g2", "sp__") == "sp__g1 or sp__g2"

    def test_nested(self, tools):
        result = tools._prefix_gpr("(g1 and g2) or g3", "sp__")
        assert "sp__g1" in result
        assert "sp__g2" in result
        assert "sp__g3" in result
        assert " and " in result
        assert " or " in result

    def test_empty_string(self, tools):
        assert tools._prefix_gpr("", "sp__") == ""

    def test_none_passthrough(self, tools):
        assert tools._prefix_gpr(None, "sp__") is None

    def test_word_boundary_or_in_gene_id(self, tools):
        """Gene IDs containing 'or' as a substring must NOT be split."""
        result = tools._prefix_gpr("b_orA1 and g2", "sp__")
        assert "sp__b_orA1" in result
        assert "sp__g2" in result

    def test_word_boundary_and_in_gene_id(self, tools):
        """Gene IDs containing 'and' as a substring must NOT be split."""
        result = tools._prefix_gpr("gandalf or frodo", "sp__")
        assert "sp__gandalf" in result
        assert "sp__frodo" in result

    def test_case_insensitive(self, tools):
        """AND / OR in any case should be recognised as operators."""
        result = tools._prefix_gpr("g1 AND g2 OR g3", "sp__")
        assert "sp__g1" in result
        assert "sp__g2" in result
        assert "sp__g3" in result


# ---------------------------------------------------------------------------
#  _is_extracellular
# ---------------------------------------------------------------------------

class TestIsExtracellular:
    """Compartment detection for extracellular metabolites."""

    def test_compartment_e(self, tools):
        met = cobra.Metabolite("glc_e", compartment="e")
        assert tools._is_extracellular(met) is True

    def test_compartment_c(self, tools):
        met = cobra.Metabolite("glc_c", compartment="c")
        assert tools._is_extracellular(met) is False

    def test_suffix_e(self, tools):
        met = cobra.Metabolite("glc_e", compartment="")
        assert tools._is_extracellular(met) is True

    def test_suffix_bracket(self, tools):
        met = cobra.Metabolite("glc[e]", compartment="")
        assert tools._is_extracellular(met) is True


# ---------------------------------------------------------------------------
#  _build_community_model
# ---------------------------------------------------------------------------

class TestBuildCommunityModel:
    """End-to-end community model construction."""

    def test_two_species_reaction_count(self, tiny_model, tools):
        """Community from two copies should have ~2× reactions."""
        m1 = tiny_model.copy()
        m2 = tiny_model.copy()
        community, obj_map, shared, _ = tools._build_community_model(
            [m1, m2], ["sp1", "sp2"])

        # Each species contributes 3 reactions; some exchange mets shared
        assert len(community.reactions) >= 6

    def test_prefixed_reaction_ids(self, tiny_model, tools):
        """All reaction IDs should be prefixed with the species label."""
        m1 = tiny_model.copy()
        community, obj_map, _, _ = tools._build_community_model(
            [m1], ["org1"])

        for rxn in community.reactions:
            assert rxn.id.startswith("org1__")

    def test_objective_map(self, tiny_model, tools):
        """Each species objective should appear in the objective map."""
        m1 = tiny_model.copy()
        m2 = tiny_model.copy()
        _, obj_map, _, _ = tools._build_community_model(
            [m1, m2], ["sp1", "sp2"])

        assert "sp1" in obj_map
        assert "sp2" in obj_map
        assert obj_map["sp1"].startswith("sp1__")
        assert obj_map["sp2"].startswith("sp2__")

    def test_shared_extracellular_metabolites(self, tiny_model, tools):
        """Extracellular metabolites should be shared, not duplicated."""
        m1 = tiny_model.copy()
        m2 = tiny_model.copy()

        # Add an extracellular metabolite to both models
        ext_met1 = cobra.Metabolite("glc_e", name="Glucose", compartment="e")
        ext_met2 = cobra.Metabolite("glc_e", name="Glucose", compartment="e")
        ex1 = cobra.Reaction("EX_glc_e")
        ex1.add_metabolites({ext_met1: -1.0})
        ex1.lower_bound = -10
        m1.add_reactions([ex1])

        ex2 = cobra.Reaction("EX_glc_e")
        ex2.add_metabolites({ext_met2: -1.0})
        ex2.lower_bound = -10
        m2.add_reactions([ex2])

        community, _, shared, _ = tools._build_community_model(
            [m1, m2], ["sp1", "sp2"])

        assert "glc_e" in shared

    def test_community_optimises(self, tiny_model, tools):
        """The community model should be solvable."""
        m1 = tiny_model.copy()
        m2 = tiny_model.copy()
        community, _, _, _ = tools._build_community_model(
            [m1, m2], ["sp1", "sp2"])

        sol = community.optimize()
        assert sol.status == "optimal"
        assert sol.objective_value > 0

    def test_no_share_extracellular(self, tiny_model, tools):
        """With share_extracellular=False, every met should be prefixed."""
        m1 = tiny_model.copy()
        community, _, shared, _ = tools._build_community_model(
            [m1], ["org"], share_extracellular=False)

        assert len(shared) == 0
        for met in community.metabolites:
            assert met.id.startswith("org__")

    def test_gene_prefixing_in_community(self, tiny_model, tools):
        """Gene reaction rules should have prefixed gene IDs."""
        m1 = tiny_model.copy()
        community, _, _, _ = tools._build_community_model(
            [m1], ["sp1"])

        r1 = community.reactions.get_by_id("sp1__R1")
        assert "sp1__g1" in r1.gene_reaction_rule


# ---------------------------------------------------------------------------
#  SteadyCom Objective Mode
# ---------------------------------------------------------------------------

class TestSteadyComObjective:
    """Test SteadyCom-style equal-growth community objective."""

    @pytest.fixture
    def tools(self):
        from metabodesk_core.mixin_tools import ToolsMixin

        class _Harness(ToolsMixin):
            pass

        return _Harness()

    def test_steadycom_creates_mu_drain(self, tiny_model, tools):
        """SteadyCom mode should create a community_mu_drain reaction."""
        m1 = tiny_model.copy()
        m2 = tiny_model.copy()
        community, obj_map, _, _ = tools._build_community_model(
            [m1, m2], ["sp1", "sp2"], objective_mode="steadycom")

        assert "community_mu_drain" in community.reactions
        assert "community_mu" in community.metabolites

    def test_steadycom_equal_growth(self, tiny_model, tools):
        """Both species should grow at equal rates in SteadyCom mode."""
        m1 = tiny_model.copy()
        m2 = tiny_model.copy()
        community, obj_map, _, _ = tools._build_community_model(
            [m1, m2], ["sp1", "sp2"], objective_mode="steadycom")

        sol = community.optimize()
        assert sol.status == "optimal"

        sp1_growth = sol.fluxes.get(obj_map["sp1"], 0)
        sp2_growth = sol.fluxes.get(obj_map["sp2"], 0)
        assert abs(sp1_growth - sp2_growth) < 1e-6

    def test_sum_mode_still_works(self, tiny_model, tools):
        """Default 'sum' mode should still optimise correctly."""
        m1 = tiny_model.copy()
        community, _, _, _ = tools._build_community_model(
            [m1], ["sp1"], objective_mode="sum")

        sol = community.optimize()
        assert sol.status == "optimal"
        assert sol.objective_value > 0

    def test_steadycom_single_model_falls_back(self, tiny_model, tools):
        """SteadyCom with <2 models should fall back to sum mode."""
        m1 = tiny_model.copy()
        community, _, _, _ = tools._build_community_model(
            [m1], ["sp1"], objective_mode="steadycom")

        # Should not create community_mu_drain for single model
        assert "community_mu_drain" not in community.reactions


# ---------------------------------------------------------------------------
#  Community Scalable Batch Build
# ---------------------------------------------------------------------------

class TestCommunityBatchBuild:
    """Test that batch-add metabolites/reactions works correctly."""

    @pytest.fixture
    def tools(self):
        from metabodesk_core.mixin_tools import ToolsMixin

        class _Harness(ToolsMixin):
            pass

        return _Harness()

    def test_batch_build_produces_valid_model(self, tiny_model, tools):
        """Batch-built community should be a valid solvable model."""
        m1 = tiny_model.copy()
        m2 = tiny_model.copy()
        community, obj_map, shared, conflicts = tools._build_community_model(
            [m1, m2], ["sp1", "sp2"])

        # Should have reactions from both species
        assert "sp1__R1" in community.reactions
        assert "sp2__R1" in community.reactions
        assert "sp1" in obj_map
        assert "sp2" in obj_map

        sol = community.optimize()
        assert sol.status == "optimal"
        assert sol.objective_value > 0

    def test_batch_build_three_species(self, tiny_model, tools):
        """Three-species community should build without error."""
        m1 = tiny_model.copy()
        m2 = tiny_model.copy()
        m3 = tiny_model.copy()
        community, obj_map, shared, conflicts = tools._build_community_model(
            [m1, m2, m3], ["sp1", "sp2", "sp3"])

        assert len(obj_map) == 3
        sol = community.optimize()
        assert sol.status == "optimal"


# ---------------------------------------------------------------------------
#  Community Abundance Constraints (Item 9)
# ---------------------------------------------------------------------------

class TestCommunityAbundance:
    """Test relative abundance constraints in community models."""

    @pytest.fixture
    def tools(self):
        from metabodesk_core.mixin_tools import ToolsMixin

        class _Harness(ToolsMixin):
            pass

        return _Harness()

    def test_abundance_scales_bounds(self, tiny_model, tools):
        """Species with 30% abundance should have 30% of original bounds."""
        m1 = tiny_model.copy()
        m2 = tiny_model.copy()
        abundance = {"sp1": 0.7, "sp2": 0.3}
        community, obj_map, _, _ = tools._build_community_model(
            [m1, m2], ["sp1", "sp2"], abundance=abundance)

        # sp2's R1 should have ~30% of original upper bound (1000 * 0.3 = 300)
        r1_sp2 = community.reactions.get_by_id("sp2__R1")
        assert abs(r1_sp2.upper_bound - 300.0) < 1e-6

        # sp1's R1 should have ~70% of original upper bound (1000 * 0.7 = 700)
        r1_sp1 = community.reactions.get_by_id("sp1__R1")
        assert abs(r1_sp1.upper_bound - 700.0) < 1e-6

    def test_abundance_scales_exchange_bounds(self, tiny_model, tools):
        """Exchange reaction bounds should also be scaled by abundance."""
        m1 = tiny_model.copy()
        m2 = tiny_model.copy()
        abundance = {"sp1": 0.6, "sp2": 0.4}
        community, _, _, _ = tools._build_community_model(
            [m1, m2], ["sp1", "sp2"], abundance=abundance)

        # EX_A_e has LB=-10, so sp2's version should be -10 * 0.4 = -4.0
        ex_a_sp2 = community.reactions.get_by_id("sp2__EX_A_e")
        assert abs(ex_a_sp2.lower_bound - (-4.0)) < 1e-6

    def test_no_abundance_means_full_bounds(self, tiny_model, tools):
        """Without abundance, bounds should be unchanged."""
        m1 = tiny_model.copy()
        m2 = tiny_model.copy()
        community, _, _, _ = tools._build_community_model(
            [m1, m2], ["sp1", "sp2"], abundance=None)

        r1_sp1 = community.reactions.get_by_id("sp1__R1")
        assert abs(r1_sp1.upper_bound - 1000.0) < 1e-6

    def test_abundance_with_steadycom(self, tiny_model, tools):
        """Abundance + SteadyCom should both apply."""
        m1 = tiny_model.copy()
        m2 = tiny_model.copy()
        abundance = {"sp1": 0.5, "sp2": 0.5}
        community, obj_map, _, _ = tools._build_community_model(
            [m1, m2], ["sp1", "sp2"],
            objective_mode="steadycom", abundance=abundance)

        # SteadyCom should still create mu_drain
        assert "community_mu_drain" in community.reactions

        # Bounds should be scaled
        r1_sp1 = community.reactions.get_by_id("sp1__R1")
        assert abs(r1_sp1.upper_bound - 500.0) < 1e-6

        sol = community.optimize()
        assert sol.status == "optimal"
