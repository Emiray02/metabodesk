"""Shared test fixtures and configuration for the MetaboDesk test suite."""

import pytest


@pytest.fixture
def tiny_model():
    """Create a minimal valid COBRApy model for testing.

    The model has:
    - 2 metabolites: A_c (cytoplasm), B_c (cytoplasm)
    - 3 reactions:
        - EX_A_e: A_c ⇌  (exchange, LB=-10, UB=1000, imports up to 10)
        - R1:     A_c → B_c (internal, LB=0, UB=1000)
        - EX_B_e: B_c →  (exchange/demand, LB=0, UB=1000, objective)
    - 1 gene: g1 (associated with R1)
    """
    import cobra

    model = cobra.Model("tiny")

    a = cobra.Metabolite("A_c", name="Metabolite A", compartment="c")
    b = cobra.Metabolite("B_c", name="Metabolite B", compartment="c")

    ex_a = cobra.Reaction("EX_A_e")
    ex_a.name = "A exchange"
    ex_a.lower_bound = -10.0
    ex_a.upper_bound = 1000.0
    ex_a.add_metabolites({a: -1.0})  # standard convention: negative = consumed

    r1 = cobra.Reaction("R1")
    r1.name = "A to B"
    r1.lower_bound = 0.0
    r1.upper_bound = 1000.0
    r1.add_metabolites({a: -1.0, b: 1.0})
    r1.gene_reaction_rule = "g1"

    ex_b = cobra.Reaction("EX_B_e")
    ex_b.name = "B demand"
    ex_b.lower_bound = 0.0
    ex_b.upper_bound = 1000.0
    ex_b.add_metabolites({b: -1.0})

    model.add_reactions([ex_a, r1, ex_b])
    model.objective = "EX_B_e"

    return model


@pytest.fixture
def ecoli_core(tmp_path):
    """Download (or cache) the E. coli core model for integration tests.

    Returns None if the download fails so tests can skip gracefully.
    """
    import cobra
    try:
        # COBRApy ships a bundled test model
        from cobra.io import load_model
        return load_model("textbook")
    except Exception:
        return None
