"""Tests for metabodesk_core.mixin_network — graph construction logic."""

import pytest

cobra = pytest.importorskip("cobra")
nx = pytest.importorskip("networkx")


class TestBuildNetworkGraphData:
    """Tests for the thread-safe graph builder.

    These tests instantiate only the ``NetworkMixin`` (not the full
    ``MainWindow``) and call ``_build_network_graph_data`` directly.
    """

    @staticmethod
    def _default_params(**overrides) -> dict:
        """Return a params dict with sensible defaults."""
        params = {
            "search_text": "",
            "view_mode": "Flux magnitude",
            "depth": 2,
            "threshold": 0.0,
            "max_nodes": 5000,
            "max_edges": 0,
            "min_degree": 0,
            "topn_rxn": 0,
            "layout_choice": "Spring",
            "only_connected": False,
            "hide_orphans": False,
            "show_legend": False,
            "exchange_only": False,
            "objective_only": False,
            "subsystem_filter": "All subsystems",
            "flux_map": {},
            "fva_map": {},
            "sgd_map": {},
            "wt_growth": None,
            "obj_combo_data": None,
        }
        params.update(overrides)
        return params

    def test_tiny_model_full_graph(self, tiny_model):
        """Full graph of tiny model should include all reactions and metabolites."""
        from metabodesk_core.mixin_network import NetworkMixin
        mixin = NetworkMixin()
        params = self._default_params()
        result = mixin._build_network_graph_data(tiny_model, params)

        assert "error" not in result
        G = result["G"]
        assert G.number_of_nodes() > 0
        assert G.number_of_edges() > 0
        assert len(result["rxn_nodes"]) == 3  # EX_A_e, R1, EX_B_e
        assert len(result["met_nodes"]) == 2  # A_c, B_c

    def test_search_filter(self, tiny_model):
        """Searching for 'R1' should produce a focused subgraph."""
        from metabodesk_core.mixin_network import NetworkMixin
        mixin = NetworkMixin()
        params = self._default_params(search_text="r1")
        result = mixin._build_network_graph_data(tiny_model, params)

        assert "error" not in result
        assert "R1" in [n for n in result["rxn_nodes"]]

    def test_no_match_returns_error(self, tiny_model):
        """Searching for a nonexistent ID should return an error dict."""
        from metabodesk_core.mixin_network import NetworkMixin
        mixin = NetworkMixin()
        params = self._default_params(search_text="zzz_nonexistent_zzz")
        result = mixin._build_network_graph_data(tiny_model, params)

        assert "error" in result

    def test_exchange_only_filter(self, tiny_model):
        """exchange_only should exclude internal reactions."""
        from metabodesk_core.mixin_network import NetworkMixin
        mixin = NetworkMixin()
        params = self._default_params(exchange_only=True)
        result = mixin._build_network_graph_data(tiny_model, params)

        assert "error" not in result
        # R1 is internal — should not appear in rxn_nodes
        assert "R1" not in result["rxn_nodes"]

    def test_hide_orphans(self, tiny_model):
        """With hide_orphans, zero-degree nodes should be removed."""
        from metabodesk_core.mixin_network import NetworkMixin
        mixin = NetworkMixin()
        params = self._default_params(hide_orphans=True)
        result = mixin._build_network_graph_data(tiny_model, params)

        assert "error" not in result
        G = result["G"]
        for n in G.nodes():
            assert G.degree(n) > 0

    def test_layout_circular(self, tiny_model):
        """Circular layout should produce valid positions."""
        from metabodesk_core.mixin_network import NetworkMixin
        mixin = NetworkMixin()
        params = self._default_params(layout_choice="Circular")
        result = mixin._build_network_graph_data(tiny_model, params)

        assert "error" not in result
        assert len(result["pos"]) == result["G"].number_of_nodes()

    def test_layout_kamada_kawai(self, tiny_model):
        """Kamada-Kawai layout should produce valid positions."""
        from metabodesk_core.mixin_network import NetworkMixin
        mixin = NetworkMixin()
        params = self._default_params(layout_choice="Kamada-Kawai")
        result = mixin._build_network_graph_data(tiny_model, params)

        assert "error" not in result
        assert len(result["pos"]) == result["G"].number_of_nodes()
