"""Network map visualization mixin for MetaboDesk.

Provides interactive metabolic network graph visualization using NetworkX
for graph construction and layout computation, and matplotlib for rendering.

Key features
------------
- Bipartite graph representation: reaction nodes (circles) and metabolite
  nodes (squares) connected by directed edges based on stoichiometry.
- Multiple layout algorithms: Spring (default), Kamada-Kawai, Circular.
- Colour-mapping by flux magnitude, FVA range, or gene-knockout impact.
- Neighbourhood search with configurable BFS depth.
- Background graph computation via ``QThread`` for large models (≥500
  reactions) so the UI remains responsive during layout computation.
- Interactive click-to-inspect and hover tooltips.
- Graph and image export (CSV, PNG, SVG, PDF).
"""

import csv
import logging

import cobra
# networkx is imported lazily on first use — saves ~0.5 s at startup
_nx = None

def _get_nx():
    global _nx
    if _nx is None:
        import networkx
        _nx = networkx
    return _nx

from pathlib import Path

from PySide6.QtWidgets import QFileDialog, QMessageBox

# matplotlib imports are lazy — loaded only when drawing the map
_cm = None
_Normalize = None


def _get_mpl():
    global _cm, _Normalize
    if _cm is None:
        from matplotlib import cm
        from matplotlib.colors import Normalize
        _cm = cm
        _Normalize = Normalize
    return _cm, _Normalize

from metabodesk_core.utils import (
    is_exchange_reaction, get_kegg_rxn_id, get_ec_numbers,
    get_description, get_gpr,
)

logger = logging.getLogger("MetaboDesk")


class NetworkMixin:
    """Mixin providing network visualization and graph export functionality.

    Builds a ``networkx.DiGraph`` from the loaded metabolic model, applies
    user-selected filters (subsystem, search text, flux threshold, …),
    computes a graph layout, and renders the result onto a matplotlib canvas
    embedded in the Qt UI.

    For models with ≥ 500 reactions the CPU-intensive graph construction
    and layout computation are executed in a background ``QThread`` while
    drawing is always performed on the main (GUI) thread.
    """

    # ================================================================== #
    #  Public entry point                                                  #
    # ================================================================== #

    def _render_network_map(self):
        """Trigger network map rendering.

        For small models (< 500 reactions) computation and drawing happen
        synchronously on the main thread.  For larger models, graph
        construction and layout are offloaded to a ``QThread`` to keep the
        UI responsive.
        """
        if self.is_running:
            return
        if self.base_model is None:
            self._set_map_placeholder("Load model first")
            return

        params = self._collect_map_params()
        model = self.base_model

        if len(model.reactions) < 500:
            # Small model — run synchronously
            self.map_render_btn.setEnabled(False)
            try:
                data = self._build_network_graph_data(model, params)
                if "error" in data:
                    self._set_map_placeholder(data["error"])
                else:
                    self._draw_network_graph(data, params)
            except Exception as e:
                logger.error(f"Map error: {e}")
                self._set_map_placeholder(f"Error: {str(e)[:60]}")
            finally:
                self.map_render_btn.setEnabled(True)
        else:
            # Large model — offload heavy work to QThread
            self.map_render_btn.setEnabled(False)

            def _build(worker=None):
                if worker:
                    worker.report_progress("Building network graph...", 10)
                return self._build_network_graph_data(model, params, worker=worker)

            def _on_done(result):
                try:
                    if "error" in result:
                        self._set_map_placeholder(result["error"])
                    else:
                        self._draw_network_graph(result, params)
                except Exception as e:
                    logger.error(f"Map render error: {e}")
                    self._set_map_placeholder(f"Error: {e}")
                finally:
                    self.map_render_btn.setEnabled(True)
                    self.statusBar().showMessage("Ready.")

            self._launch_worker(_build, _on_done, "Building network graph...")

    # ================================================================== #
    #  Canvas helpers                                                      #
    # ================================================================== #

    def _set_map_placeholder(self, text: str):
        """Display a centred text message on the map canvas."""
        ax = self.map_canvas.ax
        ax.clear()
        ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=11, transform=ax.transAxes)
        ax.axis("off")
        self.map_canvas.draw()

    def _focus_on_selected_map_node(self):
        """Re-render map centred on the last-clicked node."""
        node = getattr(self, "_selected_map_node", None)
        if not node:
            QMessageBox.information(self, "Map", "Click a node on the map first.")
            return
        self.map_search.setText(str(node))
        self._render_network_map()

    # ================================================================== #
    #  Export                                                               #
    # ================================================================== #

    def _export_map_image(self):
        """Save the current map canvas to PNG, SVG, or PDF."""
        if not hasattr(self, "map_canvas"):
            return
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export map image",
            str(Path.home() / "network_map.png"),
            "PNG files (*.png);;SVG files (*.svg);;PDF files (*.pdf)"
        )
        if not file_path:
            return
        try:
            self.map_canvas.figure.savefig(file_path, dpi=300)
            self.statusBar().showMessage(f"Map image exported: {file_path}")
        except Exception as e:
            self._show_error("Export failed", "Failed to export map image.", e)

    def _export_map_csv(self):
        """Export graph nodes and edges as separate CSV files."""
        if not hasattr(self, "_current_graph") or self._current_graph is None:
            QMessageBox.warning(self, "Map", "Render the map first.")
            return
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export graph CSV",
            str(Path.home() / "network_graph.csv"),
            "CSV files (*.csv)"
        )
        if not file_path:
            return
        try:
            base = Path(file_path)
            nodes_path = base.with_name(base.stem + "_nodes.csv")
            edges_path = base.with_name(base.stem + "_edges.csv")

            with open(nodes_path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["id", "type", "name", "weight", "degree"])
                for n, d in self._current_graph.nodes(data=True):
                    w.writerow([
                        n,
                        d.get("type", ""),
                        d.get("name", ""),
                        d.get("weight", ""),
                        self._current_graph.degree(n),
                    ])

            with open(edges_path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["source", "target"])
                for u, v in self._current_graph.edges():
                    w.writerow([u, v])

            self.statusBar().showMessage(f"Graph exported: {nodes_path.name}, {edges_path.name}")
        except Exception as e:
            self._show_error("Export failed", "Failed to export graph CSV.", e)

    # ================================================================== #
    #  Backward-compatible aliases                                         #
    # ================================================================== #

    def _render_network_map_bg(self):
        """Backward-compatible alias — delegates to the refactored entry point."""
        self._render_network_map()

    def _render_network_map_impl(self):
        """Backward-compatible wrapper.

        Previously contained the full rendering pipeline (~450 lines).
        Now delegates to ``_build_network_graph_data`` +
        ``_draw_network_graph``.
        """
        if self.base_model is None:
            self._set_map_placeholder("Load model first")
            return
        try:
            params = self._collect_map_params()
            data = self._build_network_graph_data(self.base_model, params)
            if "error" in data:
                self._set_map_placeholder(data["error"])
            else:
                self._draw_network_graph(data, params)
        except Exception as e:
            logger.error(f"Map error: {e}")
            self._set_map_placeholder(f"Error: {str(e)[:40]}")

    # ================================================================== #
    #  UI state snapshot                                                   #
    # ================================================================== #

    def _collect_map_params(self) -> dict:
        """Snapshot all map-related UI widget values into a plain dict.

        Called on the main thread before handing computation off to a
        background ``QThread``.  The resulting dict is fully serialisable
        and contains no references to Qt widgets.

        Returns
        -------
        dict
            Keys correspond to every user-adjustable map parameter
            (search text, depth, threshold, layout choice, etc.) plus
            the cached flux/FVA/SGD data from the last analysis run.
        """
        flux_map: dict = {}
        fva_map: dict = {}
        sgd_map: dict = {}
        wt_growth = None

        if self.last_run:
            try:
                if "baseline" in self.last_run and "flux" in self.last_run["baseline"]:
                    flux_map = self.last_run["baseline"].get("flux", {}) or {}
            except Exception:
                pass
            try:
                if "fva_result" in self.last_run:
                    fva_map = self.last_run.get("fva_result", {}) or {}
            except Exception:
                pass
            try:
                if "sgd_result" in self.last_run:
                    sgd_map = self.last_run.get("sgd_result", {}) or {}
                    wt_growth = self.last_run.get("wt_growth", None)
            except Exception:
                pass

        obj_id = None
        try:
            obj_id = self.objective_combo.currentData()
        except Exception:
            pass

        return {
            "search_text": (self.map_search.text() or "").strip().lower(),
            "view_mode": self.map_view_combo.currentText(),
            "depth": int(self.map_depth.value()),
            "threshold": float(self.map_threshold.value()),
            "max_nodes": int(self.map_max_nodes.value()) if hasattr(self, "map_max_nodes") else 5000,
            "max_edges": int(self.map_max_edges.value()) if hasattr(self, "map_max_edges") else 0,
            "min_degree": int(self.map_min_degree.value()) if hasattr(self, "map_min_degree") else 0,
            "topn_rxn": int(self.map_topn_rxn.value()) if hasattr(self, "map_topn_rxn") else 0,
            "layout_choice": self.map_layout_combo.currentText() if hasattr(self, "map_layout_combo") else "Spring",
            "only_connected": bool(self.map_only_connected_chk.isChecked()) if hasattr(self, "map_only_connected_chk") else False,
            "hide_orphans": bool(self.map_hide_orphans_chk.isChecked()) if hasattr(self, "map_hide_orphans_chk") else False,
            "show_legend": bool(self.map_show_legend_chk.isChecked()) if hasattr(self, "map_show_legend_chk") else False,
            "exchange_only": bool(self.map_exchange_only_chk.isChecked()) if hasattr(self, "map_exchange_only_chk") else False,
            "objective_only": bool(self.map_objective_only_chk.isChecked()) if hasattr(self, "map_objective_only_chk") else False,
            "subsystem_filter": self.map_subsystem_combo.currentText() if hasattr(self, "map_subsystem_combo") else "All subsystems",
            "flux_map": flux_map,
            "fva_map": fva_map,
            "sgd_map": sgd_map,
            "wt_growth": wt_growth,
            "obj_combo_data": obj_id,
        }

    # ================================================================== #
    #  Graph construction  (thread-safe — no Qt widget access)             #
    # ================================================================== #

    def _build_network_graph_data(self, model: cobra.Model, params: dict,
                                  *, worker=None) -> dict:
        """Build a NetworkX ``DiGraph`` and compute layout positions.

        **Thread-safe**: reads only from *model* (treated as immutable
        during the call) and *params* (a plain ``dict`` captured on the
        main thread).  Never accesses Qt widgets.

        Parameters
        ----------
        model : cobra.Model
            The metabolic model (read-only during the call).
        params : dict
            UI state captured by :meth:`_collect_map_params`.
        worker : AnalysisWorker, optional
            When supplied, progress updates are emitted via
            ``worker.report_progress(msg, pct)``.

        Returns
        -------
        dict
            Either ``{"error": "message"}`` on failure, or a dict with
            keys ``G``, ``pos``, ``rxn_map``, ``met_map``, ``rxn_nodes``,
            ``met_nodes``.
        """
        nx = _get_nx()
        # -- unpack parameters -------------------------------------------------
        search_text = params["search_text"]
        view_mode = params["view_mode"]
        depth = params["depth"]
        threshold = params["threshold"]
        max_nodes = params["max_nodes"]
        max_edges = params["max_edges"]
        min_degree = params["min_degree"]
        topn_rxn = params["topn_rxn"]
        layout_choice = params["layout_choice"]
        only_connected = params["only_connected"]
        hide_orphans = params["hide_orphans"]
        exchange_only = params["exchange_only"]
        objective_only = params["objective_only"]
        subsystem_filter = params["subsystem_filter"]
        flux_map = params["flux_map"]
        fva_map = params["fva_map"]
        sgd_map = params["sgd_map"]
        wt_growth = params["wt_growth"]

        # -- local helper closures (no *self* access) --------------------------
        def rxn_allowed(rxn: cobra.Reaction) -> bool:
            """Check whether *rxn* passes the active UI filters."""
            if exchange_only and not is_exchange_reaction(rxn):
                return False
            if subsystem_filter and subsystem_filter != "All subsystems":
                sub = str(getattr(rxn, "subsystem", "") or "").strip()
                if sub != subsystem_filter:
                    return False
            return True

        def rxn_weight(rid: str) -> float:
            """Return the metric value used for sizing / colouring *rid*."""
            if "FVA" in view_mode:
                mm = fva_map.get(rid, {}) or {}
                return max(float(mm.get("max", 0.0)) - float(mm.get("min", 0.0)), 0.0)
            if "Gene KO" in view_mode:
                if not sgd_map or wt_growth in (None, 0):
                    return 0.0
                try:
                    rxn = model.reactions.get_by_id(rid)
                    gene_vals = [float(sgd_map.get(g.id, None)) for g in rxn.genes if g.id in sgd_map]
                    if not gene_vals:
                        return 0.0
                    min_growth = min(gene_vals)
                    return max((float(wt_growth) - float(min_growth)) / float(wt_growth), 0.0)
                except Exception:
                    return 0.0
            # Default: flux magnitude
            return abs(float(flux_map.get(rid, 0.0)))

        # -- metric availability -----------------------------------------------
        if "FVA" in view_mode:
            has_metric_data = bool(fva_map)
        elif "Gene KO" in view_mode:
            has_metric_data = bool(sgd_map) and wt_growth not in (None, 0)
        else:
            has_metric_data = bool(flux_map)

        effective_threshold = threshold if has_metric_data else 0.0

        filtered_reactions = [r for r in model.reactions if rxn_allowed(r)]

        # -- seed nodes --------------------------------------------------------
        seed_rxns: list[str] = []
        seed_mets: list[str] = []
        if search_text:
            for r in filtered_reactions:
                if search_text in r.id.lower() or search_text in (r.name or "").lower():
                    seed_rxns.append(r.id)
            for m in model.metabolites:
                if search_text in m.id.lower() or search_text in (m.name or "").lower():
                    seed_mets.append(m.id)
        else:
            if not filtered_reactions:
                return {"error": "No reactions found for filters"}

            candidates = filtered_reactions
            if effective_threshold > 0:
                candidates = [r for r in candidates if rxn_weight(r.id) >= effective_threshold]
                if not candidates:
                    candidates = filtered_reactions

            candidates = sorted(candidates, key=lambda r: rxn_weight(r.id), reverse=True)
            if topn_rxn > 0:
                seed_rxns = [r.id for r in candidates[:topn_rxn]]
            else:
                seed_rxns = [r.id for r in candidates[:2000]]

        if objective_only:
            obj_id = params.get("obj_combo_data")
            if not obj_id:
                try:
                    coeffs = model.objective.expression.as_coefficients_dict()
                    obj_rxns = [r.id for r in coeffs.keys()]
                    if obj_rxns:
                        obj_id = obj_rxns[0]
                except Exception:
                    obj_id = None
            if obj_id:
                obj_id = str(obj_id)
                seed_rxns = [obj_id] + [r for r in seed_rxns if r != obj_id]

        if not seed_rxns and not seed_mets:
            return {"error": "No matches found"}

        if worker:
            worker.report_progress("Building graph nodes and edges...", 30)

        # -- build NetworkX graph ----------------------------------------------
        G = nx.DiGraph()
        rxn_map: dict[str, cobra.Reaction] = {}
        met_map: dict[str, cobra.Metabolite] = {}

        if (not search_text and effective_threshold == 0 and topn_rxn == 0
                and not only_connected and not objective_only):
            # Full graph: include all filtered reactions
            for rxn in filtered_reactions:
                w = rxn_weight(rxn.id)
                rxn_map[rxn.id] = rxn
                G.add_node(rxn.id, type="rxn", weight=w, name=(rxn.name or ""))
                for met, coeff in rxn.metabolites.items():
                    met_map[met.id] = met
                    G.add_node(met.id, type="met", name=(met.name or ""))
                    if coeff < 0:
                        G.add_edge(met.id, rxn.id)
                    else:
                        G.add_edge(rxn.id, met.id)
            # Orphan metabolites
            if not hide_orphans:
                for met in model.metabolites:
                    if met.id not in met_map:
                        met_map[met.id] = met
                        G.add_node(met.id, type="met", name=(met.name or ""))
        else:
            # BFS from seed nodes with depth limit
            max_n = max(100, max_nodes)
            queue: list[tuple] = []
            for rid in seed_rxns:
                queue.append(("rxn", rid, 0, True))
            for mid in seed_mets:
                queue.append(("met", mid, 0, True))

            visited: set[tuple] = set()

            while queue and G.number_of_nodes() < max_n:
                ntype, nid, d, is_seed = queue.pop(0)
                key = (ntype, nid)
                if key in visited:
                    continue
                visited.add(key)
                if d > depth:
                    continue

                if ntype == "rxn":
                    try:
                        rxn = model.reactions.get_by_id(nid)
                    except Exception:
                        continue
                    if not rxn_allowed(rxn):
                        continue
                    w = rxn_weight(rxn.id)
                    if not is_seed and effective_threshold > 0 and w < effective_threshold:
                        continue
                    rxn_map[rxn.id] = rxn
                    G.add_node(rxn.id, type="rxn", weight=w, name=(rxn.name or ""))
                    for met, coeff in rxn.metabolites.items():
                        met_map[met.id] = met
                        G.add_node(met.id, type="met", name=(met.name or ""))
                        if coeff < 0:
                            G.add_edge(met.id, rxn.id)
                        else:
                            G.add_edge(rxn.id, met.id)
                        if d + 1 <= depth:
                            queue.append(("met", met.id, d + 1, False))
                else:  # metabolite node
                    try:
                        met = model.metabolites.get_by_id(nid)
                    except Exception:
                        continue
                    met_map[met.id] = met
                    G.add_node(met.id, type="met", name=(met.name or ""))
                    for rxn in met.reactions:
                        if not rxn_allowed(rxn):
                            continue
                        w = rxn_weight(rxn.id)
                        if not is_seed and effective_threshold > 0 and w < effective_threshold:
                            continue
                        rxn_map[rxn.id] = rxn
                        G.add_node(rxn.id, type="rxn", weight=w, name=(rxn.name or ""))
                        coeff = rxn.metabolites.get(met, 0)
                        if coeff < 0:
                            G.add_edge(met.id, rxn.id)
                        else:
                            G.add_edge(rxn.id, met.id)
                        if d + 1 <= depth:
                            queue.append(("rxn", rxn.id, d + 1, False))

        if worker:
            worker.report_progress("Filtering graph...", 50)

        # -- apply post-build filters ------------------------------------------
        if hide_orphans:
            orphans = [n for n in G.nodes() if G.degree(n) == 0]
            if orphans:
                G.remove_nodes_from(orphans)

        if min_degree > 0:
            while True:
                low = [n for n in G.nodes() if G.degree(n) < min_degree]
                if not low:
                    break
                G.remove_nodes_from(low)

        # Enforce max-node limit (drop lowest-weight reactions first)
        if max_nodes and G.number_of_nodes() > max_nodes:
            scored = []
            for n, nd in G.nodes(data=True):
                score = float(nd.get("weight", 0.0)) if nd.get("type") == "rxn" else 0.0
                scored.append((score, n))
            scored.sort(key=lambda x: x[0])
            to_remove = [n for _, n in scored[:max(0, G.number_of_nodes() - max_nodes)]]
            if to_remove:
                G.remove_nodes_from(to_remove)

        # Enforce max-edge limit
        if max_edges and G.number_of_edges() > max_edges:
            edges = list(G.edges())

            def edge_score(e):
                w = 0.0
                for node in e:
                    nd = G.nodes.get(node, {})
                    if nd.get("type") == "rxn":
                        w = max(w, float(nd.get("weight", 0.0)))
                return w

            edges.sort(key=edge_score, reverse=True)
            keep = set(edges[:max_edges])
            remove = [e for e in edges if e not in keep]
            if remove:
                G.remove_edges_from(remove)

        if G.number_of_nodes() == 0:
            return {"error": "Empty network"}

        if worker:
            worker.report_progress("Computing layout...", 65)

        # -- layout computation (most CPU-intensive step) ----------------------
        try:
            if layout_choice == "Circular":
                pos = nx.circular_layout(G)
            elif layout_choice == "Kamada-Kawai":
                pos = nx.kamada_kawai_layout(G)
            else:
                its = 10 if G.number_of_nodes() > 2000 else 25
                pos = nx.spring_layout(G, k=0.7, iterations=its, seed=42)
        except Exception:
            pos = nx.circular_layout(G)

        if worker:
            worker.report_progress("Preparing visualization...", 90)

        rxn_nodes = [n for n, nd in G.nodes(data=True) if nd.get("type") == "rxn"]
        met_nodes = [n for n, nd in G.nodes(data=True) if nd.get("type") == "met"]

        return {
            "G": G,
            "pos": pos,
            "rxn_map": rxn_map,
            "met_map": met_map,
            "rxn_nodes": rxn_nodes,
            "met_nodes": met_nodes,
        }

    # ================================================================== #
    #  Drawing  (main-thread only — accesses Qt canvas)                    #
    # ================================================================== #

    def _draw_network_graph(self, data: dict, params: dict):
        """Render a pre-built graph onto the matplotlib canvas.

        Must be called on the **main (GUI) thread** because it writes to
        the Qt-embedded matplotlib ``FigureCanvas``.

        Parameters
        ----------
        data : dict
            Output of :meth:`_build_network_graph_data` containing ``G``,
            ``pos``, ``rxn_map``, ``met_map``, ``rxn_nodes``, ``met_nodes``.
        params : dict
            UI state from :meth:`_collect_map_params`.
        """
        nx = _get_nx()
        cm, Normalize = _get_mpl()
        G = data["G"]
        pos = data["pos"]
        rxn_map = data["rxn_map"]
        met_map = data["met_map"]
        rxn_nodes = data["rxn_nodes"]
        met_nodes = data["met_nodes"]

        search_text = params["search_text"]
        view_mode = params["view_mode"]
        depth = params["depth"]
        threshold = params["threshold"]
        show_legend = params["show_legend"]

        # -- colour / size computation -----------------------------------------
        weights = [float(G.nodes[n].get("weight", 0.0)) for n in rxn_nodes]
        max_w = max(weights) if weights else 1.0
        norm_w = [(w / max_w) if max_w > 0 else 0 for w in weights]

        rxn_colors = [cm.viridis(0.15 + 0.85 * w) for w in norm_w]
        rxn_sizes = [220 + (w * 380) for w in norm_w]

        ax = self.map_canvas.ax
        ax.clear()

        if getattr(self, "_map_colorbar", None) is not None:
            try:
                self._map_colorbar.remove()
            except Exception:
                pass
            self._map_colorbar = None

        # -- identify highlighted nodes (search matches) -----------------------
        highlighted_rxns: list[str] = []
        highlighted_mets: list[str] = []
        if search_text:
            for n in rxn_nodes:
                if search_text in n.lower() or search_text in (G.nodes[n].get("name", "") or "").lower():
                    highlighted_rxns.append(n)
            for n in met_nodes:
                if search_text in n.lower() or search_text in (G.nodes[n].get("name", "") or "").lower():
                    highlighted_mets.append(n)

        # -- draw edges --------------------------------------------------------
        nx.draw_networkx_edges(
            G, pos, ax=ax,
            arrows=True, arrowstyle="-|>", arrowsize=8,
            alpha=0.25, width=0.6, edge_color="#a0a0a0",
        )

        # -- draw metabolite nodes ---------------------------------------------
        if met_nodes:
            non_hl_mets = [n for n in met_nodes if n not in highlighted_mets]
            if non_hl_mets:
                nx.draw_networkx_nodes(
                    G, pos, nodelist=non_hl_mets, node_shape="s",
                    node_color="#9aa0a6", node_size=90, ax=ax, alpha=0.85,
                    edgecolors="#5f6368", linewidths=0.4,
                )
            if highlighted_mets:
                nx.draw_networkx_nodes(
                    G, pos, nodelist=highlighted_mets, node_shape="s",
                    node_color="#ffeb3b", node_size=140, ax=ax, alpha=1.0,
                    edgecolors="#d32f2f", linewidths=2.0,
                )

        # -- draw reaction nodes -----------------------------------------------
        if rxn_nodes:
            non_hl_rxns = [n for n in rxn_nodes if n not in highlighted_rxns]
            non_hl_colors = [rxn_colors[rxn_nodes.index(n)] for n in non_hl_rxns]
            non_hl_sizes = [rxn_sizes[rxn_nodes.index(n)] for n in non_hl_rxns]
            if non_hl_rxns:
                nx.draw_networkx_nodes(
                    G, pos, nodelist=non_hl_rxns, node_shape="o",
                    node_color=non_hl_colors, node_size=non_hl_sizes,
                    ax=ax, alpha=0.9, edgecolors="#1f2937", linewidths=0.5,
                )
            if highlighted_rxns:
                hl_colors = [rxn_colors[rxn_nodes.index(n)] for n in highlighted_rxns]
                hl_sizes = [rxn_sizes[rxn_nodes.index(n)] + 150 for n in highlighted_rxns]
                nx.draw_networkx_nodes(
                    G, pos, nodelist=highlighted_rxns, node_shape="o",
                    node_color=hl_colors, node_size=hl_sizes, ax=ax, alpha=1.0,
                    edgecolors="#d32f2f", linewidths=2.5,
                )

        # -- labels for small graphs -------------------------------------------
        if G.number_of_nodes() <= 40:
            labels = {n: n for n in rxn_nodes}
            nx.draw_networkx_labels(G, pos, labels, font_size=6, ax=ax)

        # -- colour-bar legend -------------------------------------------------
        if show_legend and rxn_nodes:
            vmax = max_w if max_w > 0 else 1.0
            sm = cm.ScalarMappable(cmap=cm.viridis, norm=Normalize(vmin=0.0, vmax=vmax))
            sm.set_array([])
            try:
                self._map_colorbar = self.map_canvas.figure.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
                self._map_colorbar.set_label(view_mode)
            except Exception:
                self._map_colorbar = None

        ax.set_title(
            f"Network Map ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges) — {view_mode}",
            fontsize=10, fontweight="bold",
        )
        ax.axis("off")
        self.map_canvas.draw()

        # -- store for click / hover handlers ----------------------------------
        self._current_graph = G
        self._current_graph_pos = pos
        self._current_reactions_map = rxn_map
        self._current_metabolites_map = met_map

        self.map_info_lbl.setText(
            f"Reactions: {len(rxn_nodes)} | Metabolites: {len(met_nodes)} "
            f"| Depth: {depth} | Threshold: {threshold:g}"
        )

    # ================================================================== #
    #  Interactive handlers                                                #
    # ================================================================== #

    def _on_map_click(self, event):
        """Handle click on map nodes to show reaction/metabolite details."""
        try:
            if not hasattr(self, "_current_graph") or self._current_graph is None:
                return
            if event.inaxes is None:
                return
            pos = getattr(self, "_current_graph_pos", {})
            if not pos:
                return
            # Find nearest node to clicked coordinates
            x, y = event.xdata, event.ydata
            nearest = None
            best_d = 1e9
            for n, (nx_, ny_) in pos.items():
                d = (nx_ - x) ** 2 + (ny_ - y) ** 2
                if d < best_d:
                    best_d = d
                    nearest = n
            if nearest is None:
                return
            self._selected_map_node = nearest
            # Show details in map_details panel
            G = self._current_graph
            data = G.nodes[nearest]
            ntype = data.get("type", "")
            if ntype == "met":
                met = None
                try:
                    met = self._current_metabolites_map.get(nearest)
                except Exception:
                    met = None
                mname = getattr(met, "name", "") if met else ""
                mform = getattr(met, "formula", "") if met else ""
                mcomp = getattr(met, "compartment", "") if met else ""
                txt = (
                    f"Metabolite: {nearest}\n"
                    f"Name: {mname}\nFormula: {mform}\nCompartment: {mcomp}\n"
                    f"Reactions: {len(getattr(met, 'reactions', [])) if met else 0}"
                )
                self.map_details.setPlainText(txt)
                return

            rxn = None
            try:
                rxn = self._current_reactions_map.get(nearest)
            except Exception:
                rxn = None
            desc = get_description(rxn) if rxn else ""
            eq = getattr(rxn, "reaction", "") if rxn else ""
            gpr = get_gpr(rxn) if rxn else ""
            kegg = get_kegg_rxn_id(rxn) if rxn else ""
            ec = get_ec_numbers(rxn) if rxn else ""
            flux = float(data.get("weight", 0.0))
            txt = (
                f"Reaction: {nearest}\nDescription: {desc}\nEquation: {eq}\n"
                f"Metric: {flux:.6g} ({self.map_view_combo.currentText()})\n"
                f"GPR: {gpr}\nKEGG: {kegg}\nEC: {ec}"
            )
            self.map_details.setPlainText(txt)
        except Exception as e:
            self.map_details.setPlainText(f"Click error: {e}")

    def _on_map_hover(self, event):
        """Handle mouse hover on map to show tooltip near the nearest node."""
        try:
            if not hasattr(self, "_current_graph") or self._current_graph is None:
                return
            if event.inaxes is None:
                # Remove existing tooltip
                if self._map_tooltip:
                    try:
                        self._map_tooltip.remove()
                        self._map_tooltip = None
                        self.map_canvas.draw_idle()
                    except Exception:
                        pass
                return

            pos = getattr(self, "_current_graph_pos", {})
            if not pos:
                return

            x, y = event.xdata, event.ydata
            nearest = None
            best_d = 1e9
            threshold_dist = 0.02  # Proximity threshold

            for n, (nx_, ny_) in pos.items():
                d = (nx_ - x) ** 2 + (ny_ - y) ** 2
                if d < best_d:
                    best_d = d
                    nearest = n

            # Remove old tooltip
            if self._map_tooltip:
                try:
                    self._map_tooltip.remove()
                    self._map_tooltip = None
                except Exception:
                    pass

            if nearest is None or best_d > threshold_dist:
                self.map_canvas.draw_idle()
                return

            # Create tooltip
            G = self._current_graph
            data = G.nodes.get(nearest, {})
            ntype = data.get("type", "")
            name = data.get("name", "") or nearest

            if ntype == "met":
                label = f"[M] {nearest}\n{name}"
            else:
                weight = data.get("weight", 0)
                label = f"[R] {nearest}\n{name}\nFlux: {weight:.4g}"

            px, py = pos[nearest]
            ax = self.map_canvas.ax
            self._map_tooltip = ax.annotate(
                label,
                xy=(px, py),
                xytext=(10, 10),
                textcoords='offset points',
                fontsize=8,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#ffffcc', edgecolor='black', alpha=0.9),
                zorder=1000
            )
            self.map_canvas.draw_idle()
        except Exception as e:
            self.statusBar().showMessage(f"Tooltip error: {e}")

