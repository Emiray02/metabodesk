"""File I/O mixin for MetaboDesk.

Handles SBML model loading and saving (via COBRApy), scenario import/export
(JSON-based snapshots of knockout, overexpression and bound-override state),
and recent-scenario tracking.  Also provides the results-export pipeline
(CSV, Excel, JSON) and the 'Export All' convenience method.
"""

import csv
import json
import logging

import cobra
from datetime import datetime
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QApplication, QFileDialog, QMessageBox, QTableWidgetItem, QCompleter,
)

from metabodesk_core.constants import RECENT_MAX
from metabodesk_core.utils import _recent_file_path, is_exchange_reaction

logger = logging.getLogger("MetaboDesk")

class IOMixin:
    """Mixin providing io functionality."""

    def _load_recent_scenarios(self):
        """Load the recent-scenarios list from the JSON file on disk."""
        self.recent_scenarios = []
        try:
            recent_file = _recent_file_path()
            if recent_file.exists():
                data = json.loads(recent_file.read_text(encoding="utf-8"))
                if isinstance(data, list):
                    cleaned = []
                    for p in data:
                        p = str(p)
                        if Path(p).exists():
                            cleaned.append(p)
                    self.recent_scenarios = cleaned[:RECENT_MAX]
        except Exception:
            self.recent_scenarios = []

    def _save_recent_scenarios(self):
        """Persist the recent-scenarios list to JSON on disk."""
        try:
            recent_file = _recent_file_path()
            recent_file.write_text(json.dumps(self.recent_scenarios[:RECENT_MAX], indent=2), encoding="utf-8")
        except Exception:
            pass

    def _add_recent_scenario(self, file_path: str):
        """Add a file path to the MRU list, de-duplicate, and rebuild the menu."""
        file_path = str(Path(file_path))
        if file_path in self.recent_scenarios:
            self.recent_scenarios.remove(file_path)
        self.recent_scenarios.insert(0, file_path)
        self.recent_scenarios = self.recent_scenarios[:RECENT_MAX]
        self._save_recent_scenarios()
        self._rebuild_recent_menu()

    # ---------------- Busy state ----------------

    def _rebuild_recent_menu(self):
        """Rebuild the File → Recent Scenarios sub-menu from the MRU list."""
        if not hasattr(self, "recent_menu"):
            return
        self.recent_menu.clear()

        if not self.recent_scenarios:
            act_none = QAction("(none)", self)
            act_none.setEnabled(False)
            self.recent_menu.addAction(act_none)
            return

        for p in self.recent_scenarios:
            act = QAction(Path(p).name, self)
            act.setToolTip(p)
            act.triggered.connect(lambda checked=False, path=p: self.import_scenario_from_path(path))
            self.recent_menu.addAction(act)

        self.recent_menu.addSeparator()
        act_clear = QAction("Clear recent list", self)
        act_clear.triggered.connect(self.clear_recent_scenarios)
        self.recent_menu.addAction(act_clear)

    def clear_recent_scenarios(self):
        self.recent_scenarios = []
        self._save_recent_scenarios()
        self._rebuild_recent_menu()

    # ---- Model Statistics Dashboard ----

    def export_scenario_to_path(self, file_path: str):
        scenario = {
            "version": 4,
            "top_n": int(self.topn_spin.value()),
            "objective_reaction_id": self.objective_combo.currentData(),
            "knockout_genes": sorted(self.knockout_genes),
            "overexpression_reactions": {k: float(v) for k, v in self.overexpression_reactions.items()},
            "temporary_upper_bound_overrides": {k: float(v) for k, v in self.temporary_upper_bound_overrides.items()},
            "exchange_bounds": self._capture_exchange_bounds_from_table(),
            "reaction_bounds": {k: {"lb": float(v[0]), "ub": float(v[1])} for k, v in self.reaction_bound_overrides.items()},
        }
        Path(file_path).write_text(json.dumps(scenario, indent=2), encoding="utf-8")

    def export_scenario(self):
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model before exporting a scenario.")
            return
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export scenario", str(Path.home() / "metabodesk_scenario.json"),
            "JSON files (*.json);;All files (*.*)",
        )
        if not file_path:
            return
        try:
            self.export_scenario_to_path(file_path)
        except Exception as e:
            self._show_error("Export failed", "Failed to export scenario.", e)
            return
        self._add_recent_scenario(file_path)
        self.statusBar().showMessage(f"Scenario exported: {file_path}")

    def import_scenario(self):
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model before importing a scenario.")
            return
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Import scenario", str(Path.home()),
            "JSON files (*.json);;All files (*.*)",
        )
        if not file_path:
            return
        self.import_scenario_from_path(file_path)

    def import_scenario_from_path(self, file_path: str):
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model before importing a scenario.")
            return

        try:
            scenario = json.loads(Path(file_path).read_text(encoding="utf-8"))
            if not isinstance(scenario, dict):
                raise ValueError("Invalid scenario format (expected JSON object).")

            self.topn_spin.setValue(int(scenario.get("top_n", self.topn_spin.value())))

            obj_rid = scenario.get("objective_reaction_id", None)
            for i in range(self.objective_combo.count()):
                if self.objective_combo.itemData(i) == obj_rid:
                    self.objective_combo.setCurrentIndex(i)
                    break

            self.knockout_genes = set([str(x) for x in scenario.get("knockout_genes", [])])
            self.update_knockout_label()

            self.overexpression_reactions = {str(k): float(v) for k, v in scenario.get("overexpression_reactions", {}).items()}
            self.temporary_upper_bound_overrides = {str(k): float(v) for k, v in scenario.get("temporary_upper_bound_overrides", {}).items()}
            self.refresh_overexpression_lists()
            self.update_selected_rxn_info()

            exchange_bounds = scenario.get("exchange_bounds", {})
            if isinstance(exchange_bounds, dict):
                self._apply_exchange_bounds_to_table(exchange_bounds)

            rb = scenario.get("reaction_bounds", {})
            if isinstance(rb, dict):
                self.reaction_bound_overrides = {}
                for rid, b in rb.items():
                    self.reaction_bound_overrides[str(rid)] = (float(b.get("lb")), float(b.get("ub")))
                self.populate_reaction_table()

        except Exception as e:
            self._show_error("Import failed", "Failed to import scenario.", e)
            return

        self._add_recent_scenario(file_path)
        self._invalidate_fba_cache()
        self.statusBar().showMessage(f"Scenario imported: {file_path}")

    def _capture_exchange_bounds_from_table(self) -> dict[str, dict[str, float]]:
        """Read LB/UB values for every exchange reaction from the medium table."""
        bounds: dict[str, dict[str, float]] = {}
        for i in range(self.medium_table.rowCount()):
            rxn_id = self.medium_table.item(i, 0).text().strip()
            lb = self._parse_float(self.medium_table.item(i, 2).text())
            ub = self._parse_float(self.medium_table.item(i, 3).text())
            bounds[rxn_id] = {"lb": lb, "ub": ub}
        return bounds

    def _apply_exchange_bounds_to_table(self, bounds: dict):
        """Write saved LB/UB values back into the medium table widget."""
        rxn_to_row = {}
        for i in range(self.medium_table.rowCount()):
            rxn_id = self.medium_table.item(i, 0).text().strip()
            rxn_to_row[rxn_id] = i
        for rxn_id, b in bounds.items():
            if rxn_id not in rxn_to_row:
                continue
            row = rxn_to_row[rxn_id]
            lb = float(b.get("lb"))
            ub = float(b.get("ub"))
            self.medium_table.setItem(row, 2, QTableWidgetItem(str(lb)))
            self.medium_table.setItem(row, 3, QTableWidgetItem(str(ub)))
        self.filter_medium_table()

    # ---------------- Results export / Export all ----------------

    def export_results_to_path(self, file_path: str):
        if not self.last_run:
            raise ValueError("No last run results to export.")
        analysis_type = self.last_run.get("analysis_type", "FBA")
        with open(file_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["# MetaboDesk results export"])
            w.writerow(["# timestamp", self.last_run.get("timestamp", "")])
            w.writerow(["# analysis_type", analysis_type])
            w.writerow([])

            if analysis_type in ("FBA", "Loopless FBA"):
                # Full FBA comparison export
                for section in ("original", "baseline", "gene_knockout_only", "overexpression_only"):
                    data = self.last_run.get(section)
                    if data and isinstance(data, dict):
                        w.writerow([f"{section}_status", data.get("status", "")])
                        w.writerow([f"{section}_objective", data.get("objective", "")])
                w.writerow([])
                w.writerow(["reaction", "baseline_flux", "compared_flux", "delta", "compare_mode"])
                for row in self.last_run.get("flux_rows", []):
                    w.writerow([row.get("reaction", ""), row.get("baseline", ""),
                                row.get("compared", ""), row.get("delta", ""),
                                self.last_run.get("compare_mode", "")])

            elif analysis_type == "FVA":
                baseline = self.last_run.get("baseline", {})
                w.writerow(["baseline_status", baseline.get("status", "")])
                w.writerow(["baseline_objective", baseline.get("objective", "")])
                w.writerow([])
                w.writerow(["reaction", "min_flux", "max_flux", "range"])
                for row in self.last_run.get("flux_rows", []):
                    w.writerow([row.get("reaction", ""), row.get("min_flux", ""),
                                row.get("max_flux", ""), row.get("range", "")])

            elif analysis_type == "pFBA":
                baseline = self.last_run.get("baseline", {})
                pfba_data = self.last_run.get("pfba", {})
                w.writerow(["baseline_objective", baseline.get("objective", "")])
                w.writerow(["pfba_objective", pfba_data.get("objective", "")])
                w.writerow([])
                w.writerow(["reaction", "fba_flux", "pfba_flux", "delta"])
                for row in self.last_run.get("flux_rows", []):
                    w.writerow([row.get("reaction", ""), row.get("fba_flux", ""),
                                row.get("pfba_flux", ""), row.get("delta", "")])

            elif analysis_type == "SGD":
                w.writerow(["wt_growth", self.last_run.get("wt_growth", "")])
                w.writerow([])
                w.writerow(["gene", "growth_after_ko"])
                for gid, gval in self.last_run.get("sgd_result", {}).items():
                    w.writerow([gid, gval])

            elif analysis_type == "DGD":
                w.writerow(["wt_growth", self.last_run.get("wt_growth", "")])
                w.writerow([])
                w.writerow(["gene_pair", "growth_after_ko"])
                for pair, gval in self.last_run.get("dgd_result", {}).items():
                    w.writerow([pair, gval])

            elif analysis_type == "SRD":
                w.writerow(["wt_growth", self.last_run.get("wt_growth", "")])
                w.writerow([])
                w.writerow(["reaction", "growth_after_ko"])
                for rid, gval in self.last_run.get("srd_result", {}).items():
                    w.writerow([rid, gval])

            elif analysis_type == "Robustness":
                w.writerow(["rxn_id", self.last_run.get("rxn_id", "")])
                w.writerow(["bound_type", self.last_run.get("bound_type", "")])
                w.writerow([])
                rob = self.last_run.get("robustness_result", {})
                w.writerow(["bound_value", "objective"])
                for v, o in zip(rob.get("values", []), rob.get("objectives", [])):
                    w.writerow([v, o])

            elif analysis_type == "Production Envelope":
                w.writerow(["product_id", self.last_run.get("product_id", "")])
                w.writerow([])
                env = self.last_run.get("envelope_result", {})
                w.writerow(["growth", "product"])
                for g, p in zip(env.get("growth", []), env.get("product", [])):
                    w.writerow([g, p])

            elif analysis_type == "Sampling":
                w.writerow(["reaction", "mean", "std", "min", "max", "median"])
                for row in self.last_run.get("flux_rows", []):
                    w.writerow([row.get("reaction", ""), row.get("mean", ""),
                                row.get("std", ""), row.get("min", ""),
                                row.get("max", ""), row.get("median", "")])

            else:
                # Generic fallback: export flux_tbl contents
                if hasattr(self, 'flux_tbl') and self.flux_tbl.rowCount() > 0:
                    headers = [self.flux_tbl.horizontalHeaderItem(c).text()
                               for c in range(self.flux_tbl.columnCount())]
                    w.writerow(headers)
                    for row in range(self.flux_tbl.rowCount()):
                        vals = []
                        for col in range(self.flux_tbl.columnCount()):
                            item = self.flux_tbl.item(row, col)
                            vals.append(item.text() if item else "")
                        w.writerow(vals)

    def export_results_csv(self):
        if not self.last_run:
            QMessageBox.warning(self, "No results", "Run FBA at least once before exporting results.")
            return
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export results (CSV)", str(Path.home() / "metabodesk_results.csv"),
            "CSV files (*.csv);;All files (*.*)",
        )
        if not file_path:
            return
        try:
            self.export_results_to_path(file_path)
        except Exception as e:
            self._show_error("Export failed", "Failed to export results.", e)
            return
        self.statusBar().showMessage(f"Results exported: {file_path}")

    def export_all(self):
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return
        if not self.last_run:
            QMessageBox.warning(self, "No results", "Run FBA at least once before Export All.")
            return
        folder = QFileDialog.getExistingDirectory(self, "Select folder for export", str(Path.home()))
        if not folder:
            return
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        folder_path = Path(folder)
        scenario_path = folder_path / f"scenario_{stamp}.json"
        results_path = folder_path / f"results_{stamp}.csv"
        chart_path = folder_path / f"chart_{stamp}.png"
        try:
            self.export_scenario_to_path(str(scenario_path))
            self.export_results_to_path(str(results_path))
            self.canvas.figure.savefig(str(chart_path), dpi=200)
        except Exception as e:
            self._show_error("Export All failed", "Failed to export all results.", e)
            return
        self.statusBar().showMessage(f"Exported: {scenario_path.name}, {results_path.name}, {chart_path.name}")

    # ---------------- Model load ----------------

    def open_sbml(self):
        if self.is_running:
            return

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select SBML model",
            str(Path.home()),
            "SBML files (*.xml *.sbml);;All files (*.*)",
        )
        if not file_path:
            return

        self._load_sbml_model(file_path)

    def open_sbml_from_path(self, file_path: str):
        if self.is_running:
            return
        if not file_path:
            return
        self._load_sbml_model(file_path)

    def _load_sbml_model(self, file_path: str):
        """Common model loading logic used by open_sbml and open_sbml_from_path."""
        # Materialise all lazy tabs before populating them with model data
        self._ensure_all_tabs_built()

        self.set_busy(True, f"Loading: {file_path}")
        QApplication.processEvents()

        try:
            model = cobra.io.read_sbml_model(file_path)
        except Exception as e:
            import traceback
            traceback.print_exc()
            self._show_error("Failed to load SBML", "Could not read the SBML file.", e)
            self.set_busy(False, "Failed to load model.")
            return

        try:
            self.base_model = model
            self.original_model_snapshot = model.copy()
            self.last_run = None
            self._invalidate_fba_cache()
            self.current_sbml_path = Path(file_path)
            self.model_dirty = False
            self.knockout_genes.clear()
            self.overexpression_reactions.clear()
            self.temporary_upper_bound_overrides.clear()
            self.reaction_bound_overrides.clear()
            self.update_knockout_label()
            self.refresh_overexpression_lists()

            self.model_patch = {
                "version": 1,
                "added_metabolites": [],
                "added_genes": [],
                "added_or_updated_reactions": [],
                "disabled_reactions": [],
                "deleted_reactions": [],
                "objective_reaction_id": None,
                "notes": "",
            }
            self._refresh_patch_preview()

            self.exchange_rxns = [r for r in model.reactions if is_exchange_reaction(r)]
            self.exchange_rxns.sort(key=lambda r: r.id)
            self.original_bounds = {r.id: (float(r.lower_bound), float(r.upper_bound)) for r in self.exchange_rxns}

            self.populate_gene_list()
            self.populate_overexpression_reaction_list()
            self.populate_medium_table()
            self.populate_reaction_table()
            self.populate_objective_combo()
            self._populate_reaction_lists()
            self.populate_genes_tab()
            self._setup_map_search_completer()
            self._populate_map_subsystems()

            model_id = getattr(model, "id", "") or "(no id)"
            model_name = getattr(model, "name", "") or "(no name)"
            self.summary_lbl.setText(
                f"Model loaded.\n"
                f"ID: {model_id}\n"
                f"Name: {model_name}\n"
                f"Reactions: {len(model.reactions)} | Metabolites: {len(model.metabolites)} | Genes: {len(model.genes)}\n"
                f"Exchange reactions (heuristic): {len(self.exchange_rxns)}\n"
            )

            self.run_btn.setEnabled(True)
            self.close_uptakes_btn.setEnabled(True)
            self.reset_medium_btn.setEnabled(True)
            self.reset_reactions_btn.setEnabled(True)
            self.medium_apply_preset_btn.setEnabled(True)
            self.add_ko_btn.setEnabled(True)
            self.clear_ko_btn.setEnabled(True)
            self.add_oe_btn.setEnabled(True)
            self.clear_oe_btn.setEnabled(True)
            self.set_temp_ub_btn.setEnabled(True)
            self.clear_temp_ub_btn.setEnabled(True)
            self.clear_rxn_overrides_btn.setEnabled(True)
            self.remove_selected_overrides_btn.setEnabled(True)
            self.objective_combo.setEnabled(True)
            self.analysis_type.setEnabled(True)
            self.solver_combo.setEnabled(True)

            self.memote_btn.setEnabled(True)
            self.carveme_btn.setEnabled(True)

            self._clear_chart()
            self.flux_tbl.setRowCount(0)
            self.results_lbl.setText("Objective: -")

            self.set_busy(False, "Ready.")
            self.tabs.setCurrentWidget(self.tab_medium)
        except Exception as e:
            import traceback
            traceback.print_exc()
            self._show_error("Error after loading", "Model was loaded but UI setup failed.", e)
            self.set_busy(False, "Model loaded with errors.")

    def _populate_reaction_lists(self):
        """Populate reaction ID combos for Robustness and Envelope tabs."""
        try:
            ids = [r.id for r in self.base_model.reactions]
        except Exception:
            ids = []
        ids_sorted = sorted(ids)
        # store for filtering (labels include name)
        items: list[tuple[str, str]] = []
        for rid in ids_sorted:
            name = self._rxn_name_from_id(rid)
            label = f"{rid} — {name}".strip(" —") if name else rid
            items.append((rid, label))
        self._all_reaction_items = items

        # Completers for text inputs (match by id or name)
        labels = [label for _, label in items]
        self._rxn_label_to_id = {label: rid for rid, label in items}
        if hasattr(self, "robustness_rxn"):
            comp = QCompleter(labels, self)
            comp.setCaseSensitivity(Qt.CaseInsensitive)
            comp.setFilterMode(Qt.MatchContains)
            comp.setCompletionMode(QCompleter.PopupCompletion)
            comp.activated.connect(self._on_robustness_completer_selected)
            self.robustness_rxn.setCompleter(comp)

        if hasattr(self, "envelope_product"):
            comp2 = QCompleter(labels, self)
            comp2.setCaseSensitivity(Qt.CaseInsensitive)
            comp2.setFilterMode(Qt.MatchContains)
            comp2.setCompletionMode(QCompleter.PopupCompletion)
            comp2.activated.connect(self._on_envelope_completer_selected)
            self.envelope_product.setCompleter(comp2)

        # Robustness combo
        if hasattr(self, "robustness_combo"):
            self._filter_robustness_combo(self.robustness_search.text() if hasattr(self, "robustness_search") else "")
        # Envelope combo
        if hasattr(self, "envelope_combo"):
            self._filter_envelope_combo(self.envelope_search.text() if hasattr(self, "envelope_search") else "")

    def _filter_robustness_combo(self, text: str):
        """Filter the robustness reaction combo by search text."""
        s = (text or "").strip().lower()
        items = [item for item in getattr(self, "_all_reaction_items", []) if s in item[1].lower()]
        self.robustness_combo.blockSignals(True)
        self.robustness_combo.clear()
        for rid, label in items:
            self.robustness_combo.addItem(label, rid)
        self.robustness_combo.blockSignals(False)

    def _filter_envelope_combo(self, text: str):
        """Filter the envelope reaction combo by search text."""
        s = (text or "").strip().lower()
        items = [item for item in getattr(self, "_all_reaction_items", []) if s in item[1].lower()]
        self.envelope_combo.blockSignals(True)
        self.envelope_combo.clear()
        for rid, label in items:
            self.envelope_combo.addItem(label, rid)
        self.envelope_combo.blockSignals(False)

    def _setup_map_search_completer(self):
        """Attach a QCompleter with all reaction/metabolite IDs to the map search."""
        if self.base_model is None:
            return
        labels = []
        for r in self.base_model.reactions:
            name = getattr(r, "name", "") or ""
            labels.append(f"{r.id} — {name}".strip(" —"))
        for m in self.base_model.metabolites:
            name = getattr(m, "name", "") or ""
            labels.append(f"{m.id} — {name}".strip(" —"))
        comp = QCompleter(sorted(labels), self)
        comp.setCaseSensitivity(Qt.CaseInsensitive)
        comp.setFilterMode(Qt.MatchContains)
        comp.setCompletionMode(QCompleter.PopupCompletion)
        self.map_search.setCompleter(comp)

    def _populate_map_subsystems(self):
        """Fill the subsystem filter combo on the map tab."""
        if self.base_model is None or not hasattr(self, "map_subsystem_combo"):
            return
        subsystems = set()
        for r in self.base_model.reactions:
            sub = str(getattr(r, "subsystem", "") or "").strip()
            if sub:
                subsystems.add(sub)
        items = ["All subsystems"] + sorted(subsystems)
        self.map_subsystem_combo.blockSignals(True)
        self.map_subsystem_combo.clear()
        for s in items:
            self.map_subsystem_combo.addItem(s)
        self.map_subsystem_combo.blockSignals(False)

    def _show_map_help(self):
        """Display a help dialog explaining map controls and interactions."""
        msg = (
            "How to use the map\n\n"
            "Core controls\n"
            "• Search: Type a reaction/metabolite ID or name to focus on matches.\n"
            "• Depth: Neighborhood depth from seeds (1–5). Higher = wider context.\n"
            "• Threshold: Hides low-impact reactions for the selected view.\n"
            "• View: What the color/size represents:\n"
            "   - Flux magnitude: |flux| from last FBA\n"
            "   - Gene KO impact: growth drop from SGD\n"
            "   - FVA bounds: flux range (max–min)\n"
            "\n"
            "Size & layout\n"
            "• Max nodes: Hard cap to keep big graphs responsive.\n"
            "• Max edges: 0 = no limit; otherwise keeps only top edges.\n"
            "• Min degree: Hides nodes with degree below this value.\n"
            "• Top‑N rxns: Keep only the highest‑impact reactions (0 = off).\n"
            "• Layout: Spring / Kamada‑Kawai / Circular.\n"
            "• Show legend: Toggles the color scale.\n"
            "\n"
            "Filters\n"
            "• Only connected to search: Show only the component connected to search seeds.\n"
            "• Hide orphan metabolites: Removes metabolites with no edges.\n"
            "• Exchange only: Keep only exchange reactions.\n"
            "• Subsystem: Filter reactions by subsystem (if available in SBML).\n"
            "• Objective neighborhood: Focus around the objective reaction.\n"
            "\n"
            "Interaction & export\n"
            "• Render Map: Builds the graph with current settings.\n"
            "• Click a node: Details appear below.\n"
            "• Focus selected: Re‑centers the map around the clicked node.\n"
            "• Export map image: Save PNG/SVG/PDF.\n"
            "• Export graph CSV: Save node/edge lists for analysis."
        )
        QMessageBox.information(self, "How to use map?", msg)

    def _on_robustness_combo_changed(self, index: int):
        """Sync robustness text input when the combo selection changes."""
        if index < 0:
            return
        rid = self.robustness_combo.itemData(index)
        if rid:
            self.robustness_rxn.setText(str(rid))

    def _on_envelope_combo_changed(self, index: int):
        """Sync envelope text input when the combo selection changes."""
        if index < 0:
            return
        rid = self.envelope_combo.itemData(index)
        if rid:
            self.envelope_product.setText(str(rid))

    def _on_robustness_completer_selected(self, text: str):
        """Set robustness text input when a completer suggestion is picked."""
        rid = self._rxn_label_to_id.get(text)
        if rid:
            self.robustness_rxn.setText(str(rid))

    def _on_envelope_completer_selected(self, text: str):
        """Set envelope text input when a completer suggestion is picked."""
        rid = self._rxn_label_to_id.get(text)
        if rid:
            self.envelope_product.setText(str(rid))

    # ---------------- Model editor tabs ----------------

    def export_results_json(self):
        """Export FBA results to JSON"""
        if self.last_run is None:
            QMessageBox.warning(self, "Export Error", "No results to export. Run analysis first.")
            return
        
        path, _ = QFileDialog.getSaveFileName(self, "Export Results", "", "JSON Files (*.json)")
        if not path:
            return
        
        try:
            with open(path, 'w') as f:
                json.dump(self.last_run, f, indent=2, default=str)
            
            QMessageBox.information(self, "Export Success", f"Results exported to {path}")
        except Exception as e:
            self._show_error("Export Error", "Failed to export results.", e)

    # -------- 2. SEARCH/FILTER FUNCTIONS --------

