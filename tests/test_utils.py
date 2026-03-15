"""Tests for metabodesk_core.utils — pure-function helpers."""

import pytest
from metabodesk_core.utils import (
    is_exchange_reaction,
    _safe_float,
    _clean_identifiers_org_id,
    evaluate_gpr_expression,
    parse_reaction_equation,
    _guess_compartment_from_met_id,
    get_kegg_rxn_id,
    get_ec_numbers,
    get_description,
    get_gpr,
)


# ------------------------------------------------------------------ #
#  is_exchange_reaction                                                #
# ------------------------------------------------------------------ #

class TestIsExchangeReaction:
    """Tests for :func:`is_exchange_reaction`."""

    def test_ex_prefix(self, tiny_model):
        rxn = tiny_model.reactions.get_by_id("EX_A_e")
        assert is_exchange_reaction(rxn) is True

    def test_internal_reaction(self, tiny_model):
        rxn = tiny_model.reactions.get_by_id("R1")
        assert is_exchange_reaction(rxn) is False

    def test_demand_reaction_single_metabolite(self, tiny_model):
        rxn = tiny_model.reactions.get_by_id("EX_B_e")
        assert is_exchange_reaction(rxn) is True


# ------------------------------------------------------------------ #
#  _safe_float                                                         #
# ------------------------------------------------------------------ #

class TestSafeFloat:
    """Tests for :func:`_safe_float`."""

    def test_normal_float(self):
        assert _safe_float(3.14) == pytest.approx(3.14)

    def test_string_number(self):
        assert _safe_float("42") == pytest.approx(42.0)

    def test_none_returns_default(self):
        assert _safe_float(None) == 0.0

    def test_nan_returns_default(self):
        assert _safe_float(float("nan"), default=-1.0) == -1.0

    def test_inf_returns_default(self):
        assert _safe_float(float("inf")) == 0.0

    def test_custom_default(self):
        assert _safe_float("invalid", default=99.9) == 99.9


# ------------------------------------------------------------------ #
#  _clean_identifiers_org_id                                           #
# ------------------------------------------------------------------ #

class TestCleanIdentifiersOrgId:
    """Tests for :func:`_clean_identifiers_org_id`."""

    def test_full_url(self):
        assert _clean_identifiers_org_id("https://identifiers.org/kegg.reaction/R00001") == "R00001"

    def test_plain_id(self):
        assert _clean_identifiers_org_id("R00001") == "R00001"

    def test_empty_string(self):
        assert _clean_identifiers_org_id("") == ""

    def test_with_prefix(self):
        assert _clean_identifiers_org_id("kegg.reaction/R00001") == "R00001"


# ------------------------------------------------------------------ #
#  evaluate_gpr_expression                                             #
# ------------------------------------------------------------------ #

class TestEvaluateGPR:
    """Tests for :func:`evaluate_gpr_expression`."""

    def test_single_gene_present(self):
        assert evaluate_gpr_expression("geneA", {"geneA": 5.0}) == 5.0

    def test_single_gene_missing(self):
        assert evaluate_gpr_expression("geneA", {}) is None

    def test_and_logic(self):
        vals = {"g1": 3.0, "g2": 7.0}
        result = evaluate_gpr_expression("g1 and g2", vals)
        assert result == pytest.approx(3.0)  # AND = min

    def test_or_logic(self):
        vals = {"g1": 3.0, "g2": 7.0}
        result = evaluate_gpr_expression("g1 or g2", vals)
        assert result == pytest.approx(7.0)  # OR = max

    def test_nested_expression(self):
        vals = {"g1": 2.0, "g2": 8.0, "g3": 5.0}
        # (g1 and g2) or g3  →  min(2,8) or 5  →  max(2, 5) = 5
        result = evaluate_gpr_expression("(g1 and g2) or g3", vals)
        assert result == pytest.approx(5.0)

    def test_empty_rule(self):
        assert evaluate_gpr_expression("", {"g1": 1.0}) is None

    def test_none_rule(self):
        assert evaluate_gpr_expression(None, {"g1": 1.0}) is None

    def test_partial_data(self):
        # g1 and g2 — g2 missing → should return None (filtered out)
        result = evaluate_gpr_expression("g1 and g2", {"g1": 5.0})
        # Only g1 has data, so AND produces min of [5.0] = 5.0
        assert result == pytest.approx(5.0)


# ------------------------------------------------------------------ #
#  parse_reaction_equation                                             #
# ------------------------------------------------------------------ #

class TestParseReactionEquation:
    """Tests for :func:`parse_reaction_equation`."""

    def test_simple_irreversible(self):
        stoich, rev = parse_reaction_equation("A_c -> B_c")
        assert rev is False
        assert stoich["A_c"] == pytest.approx(-1.0)
        assert stoich["B_c"] == pytest.approx(1.0)

    def test_reversible(self):
        stoich, rev = parse_reaction_equation("A_c <=> B_c")
        assert rev is True

    def test_coefficients(self):
        stoich, rev = parse_reaction_equation("2 A_c + 3 B_c -> C_c")
        assert stoich["A_c"] == pytest.approx(-2.0)
        assert stoich["B_c"] == pytest.approx(-3.0)
        assert stoich["C_c"] == pytest.approx(1.0)

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            parse_reaction_equation("")

    def test_no_arrow_raises(self):
        with pytest.raises(ValueError, match="arrow"):
            parse_reaction_equation("A_c B_c")

    def test_double_arrow(self):
        stoich, rev = parse_reaction_equation("A_c <--> B_c")
        assert rev is True
        assert stoich["A_c"] == pytest.approx(-1.0)
        assert stoich["B_c"] == pytest.approx(1.0)


# ------------------------------------------------------------------ #
#  _guess_compartment_from_met_id                                      #
# ------------------------------------------------------------------ #

class TestGuessCompartment:
    """Tests for :func:`_guess_compartment_from_met_id`."""

    def test_cytoplasm(self):
        assert _guess_compartment_from_met_id("glc__D_c") == "c"

    def test_extracellular(self):
        assert _guess_compartment_from_met_id("glc__D_e") == "e"

    def test_no_suffix(self):
        assert _guess_compartment_from_met_id("glucose") == "c"


# ------------------------------------------------------------------ #
#  Annotation helpers                                                  #
# ------------------------------------------------------------------ #

class TestAnnotationHelpers:
    """Tests for annotation extraction functions."""

    def test_get_description(self, tiny_model):
        rxn = tiny_model.reactions.get_by_id("R1")
        assert get_description(rxn) == "A to B"

    def test_get_gpr(self, tiny_model):
        rxn = tiny_model.reactions.get_by_id("R1")
        assert "g1" in get_gpr(rxn)

    def test_get_kegg_empty(self, tiny_model):
        rxn = tiny_model.reactions.get_by_id("R1")
        assert get_kegg_rxn_id(rxn) == ""

    def test_get_ec_empty(self, tiny_model):
        rxn = tiny_model.reactions.get_by_id("R1")
        assert get_ec_numbers(rxn) == ""
