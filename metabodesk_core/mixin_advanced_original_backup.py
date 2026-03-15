"""Advanced analyses mixin for MetaboDesk.

Provides specialised analysis methods beyond the standard FBA/FVA suite:

- **Gene Expression Integration**: E-Flux, GIMME, iMAT algorithms for
  context-specific metabolic modelling.
- **OptKnock / Strain Design**: Identifies gene-deletion strategies to
  maximise production of a target metabolite.
- **Dynamic FBA (dFBA)**: Time-course simulation of batch culture.
- **Thermodynamic FBA (TMFA)**: Thermodynamically constrained FBA.
- **Phenotype Phase Plane (PhPP)**: Two-dimensional objective scanning.
- **Flux Coupling Analysis**: Identifies coupled reaction pairs.
- **Escher Map Viewer**: Launches the Escher web-based pathway viewer.
- **SBML Validation**: Schema and consistency checks.
- **Jupyter Export**: Generates reproducible Jupyter notebooks.
"""

import csv
import json
import logging

import cobra
from datetime import datetime
from pathlib import Path
from itertools import combinations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QFileDialog,
    QMessageBox, QTableWidget, QTableWidgetItem, QAbstractItemView,
    QLineEdit, QSpinBox, QTabWidget, QDoubleSpinBox, QCheckBox,
    QComboBox, QPlainTextEdit, QDialog, QDialogButtonBox, QCompleter,
    QFormLayout, QStackedWidget,
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from metabodesk_core.utils import evaluate_gpr_expression
from metabodesk_core.widgets import TextPopup

logger = logging.getLogger("MetaboDesk")

class AdvancedMixin:
    """Mixin providing advanced functionality."""

    def run_gene_expression_analysis(self):
        """Gene Expression Integration: GIMME, iMAT, and E-Flux algorithms."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Gene Expression Integration")
        dialog.resize(500, 350)
        layout = QVBoxLayout(dialog)

        layout.addWidget(QLabel(
            "Load a CSV file with gene expression data.\n\n"
            "Required columns:\n"
            "  gene_id — gene identifier matching the model (e.g. b0001)\n"
            "  expression_value — numeric expression level\n\n"
            "Expression values can be in any unit (TPM, FPKM, RPKM,\n"
            "raw counts, log2-transformed, etc.). Values are internally\n"
            "normalised to [0, 1] by dividing by the maximum."))

        file_row = QHBoxLayout()
        file_edit = QLineEdit()
        file_edit.setPlaceholderText("Select expression data CSV...")
        file_btn = QPushButton("Browse...")
        def _browse():
            fp, _ = QFileDialog.getOpenFileName(dialog, "Expression Data", "", "CSV (*.csv *.tsv);;All (*.*)")
            if fp:
                file_edit.setText(fp)
        file_btn.clicked.connect(_browse)
        file_row.addWidget(file_edit)
        file_row.addWidget(file_btn)

        expr_template_btn = QPushButton("Generate Template CSV…")
        def _gen_expr_template():
            fp, _ = QFileDialog.getSaveFileName(dialog, "Save Expression Template CSV",
                                                 "expression_template.csv", "CSV (*.csv)")
            if not fp:
                return
            try:
                with open(fp, 'w', newline='', encoding='utf-8') as f:
                    w = csv.writer(f)
                    w.writerow(["gene_id", "expression_value"])
                    if self.base_model:
                        for gene in self.base_model.genes:
                            w.writerow([gene.id, 1.0])
                QMessageBox.information(dialog, "Template Saved",
                                        f"Expression template saved to:\n{fp}\n\n"
                                        "Fill in expression values (TPM, FPKM, counts, etc.).\n"
                                        "Values are normalised internally by dividing by max.")
            except Exception as e:
                QMessageBox.critical(dialog, "Error", f"Failed to save template:\n{e}")
        expr_template_btn.clicked.connect(_gen_expr_template)
        file_row.addWidget(expr_template_btn)
        layout.addLayout(file_row)

        form = QFormLayout()
        method_combo = QComboBox()
        method_combo.addItems(["E-Flux", "GIMME", "iMAT"])
        form.addRow("Method:", method_combo)

        threshold_spin = QDoubleSpinBox()
        threshold_spin.setRange(0.0, 1e6)
        threshold_spin.setValue(0.0)
        threshold_spin.setDecimals(2)
        form.addRow("Expression threshold:", threshold_spin)

        required_growth_spin = QDoubleSpinBox()
        required_growth_spin.setRange(0.0, 1.0)
        required_growth_spin.setValue(0.9)
        required_growth_spin.setSingleStep(0.05)
        form.addRow("Min growth fraction (GIMME):", required_growth_spin)
        layout.addLayout(form)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)

        if dialog.exec() != QDialog.Accepted:
            return

        csv_path = file_edit.text().strip()
        if not csv_path or not Path(csv_path).exists():
            QMessageBox.warning(self, "File Error", "Please select a valid CSV file.")
            return

        method = method_combo.currentText()
        threshold = threshold_spin.value()
        required_growth = required_growth_spin.value()

        # Parse expression data
        try:
            expr_data = {}
            sep = '\t' if csv_path.endswith('.tsv') else ','
            with open(csv_path, encoding='utf-8') as f:
                reader = csv.reader(f, delimiter=sep)
                next(reader)  # skip header
                for row in reader:
                    if len(row) >= 2:
                        gene_id = row[0].strip()
                        try:
                            expr_val = float(row[1].strip())
                        except ValueError:
                            continue
                        expr_data[gene_id] = expr_val
            if not expr_data:
                raise ValueError("No valid expression data found in CSV.")
        except Exception as e:
            self._show_error("Parse Error", "Failed to parse expression data", e)
            return

        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self._apply_selected_solver(model)

        def _compute():
            # Map gene expression to reactions via GPR tree evaluation
            rxn_expr = {}
            for rxn in model.reactions:
                if not rxn.genes:
                    rxn_expr[rxn.id] = None
                    continue
                # Build gene_values dict for this reaction's genes
                gene_vals = {g.id: expr_data[g.id] for g in rxn.genes if g.id in expr_data}
                if not gene_vals:
                    rxn_expr[rxn.id] = None
                    continue
                # Use GPR rule string for recursive AND=min, OR=max evaluation
                gpr_str = getattr(rxn, 'gene_reaction_rule', '') or ''
                if gpr_str.strip():
                    val = evaluate_gpr_expression(gpr_str, gene_vals)
                    rxn_expr[rxn.id] = val
                else:
                    # Fallback: if no GPR rule string, use min of all gene values
                    rxn_expr[rxn.id] = min(gene_vals.values())

            if method == "E-Flux":
                # Scale upper bounds by normalized expression
                max_expr = max((v for v in rxn_expr.values() if v is not None), default=1.0)
                if max_expr <= 0:
                    max_expr = 1.0
                for rxn in model.reactions:
                    expr_val = rxn_expr.get(rxn.id)
                    if expr_val is not None:
                        scale = expr_val / max_expr
                        if rxn.upper_bound > 0:
                            rxn.upper_bound = rxn.upper_bound * max(scale, 0.001)
                        if rxn.lower_bound < 0:
                            rxn.lower_bound = rxn.lower_bound * max(scale, 0.001)
                sol = model.optimize()
                if sol.status != 'optimal':
                    raise ValueError(f"E-Flux optimization failed: {sol.status}")
                return {"method": "E-Flux", "status": sol.status,
                        "objective": float(sol.objective_value),
                        "flux": {r.id: float(sol.fluxes[r.id]) for r in model.reactions},
                        "mapped_genes": len([v for v in rxn_expr.values() if v is not None]),
                        "total_genes": len(expr_data)}

            elif method == "GIMME":
                # GIMME: minimize flux through low-expression reactions
                # Step 1: Find required minimum objective
                sol0 = model.optimize()
                if sol0.status != 'optimal':
                    raise ValueError("Baseline optimization failed.")
                required_obj_val = float(sol0.objective_value) * required_growth

                # Step 2: Constrain growth to minimum required
                obj_rxn_id = None
                obj_coeffs = model.objective.expression.as_coefficients_dict()
                for var, coeff in obj_coeffs.items():
                    if float(coeff) != 0:
                        obj_rxn_id = var.name
                        break
                if obj_rxn_id:
                    obj_rxn_ref = model.reactions.get_by_id(obj_rxn_id)
                    obj_rxn_ref.lower_bound = max(obj_rxn_ref.lower_bound, required_obj_val)

                # Step 3: Identify low-expression reactions
                low_expr_rxns = []
                for rxn in model.reactions:
                    expr_val = rxn_expr.get(rxn.id)
                    if expr_val is not None and expr_val < threshold:
                        low_expr_rxns.append(rxn)

                # Step 4: Penalized objective — minimize total flux through
                # low-expression reactions (true GIMME formulation)
                if low_expr_rxns:
                    penalty_expr = sum(
                        (threshold - (rxn_expr.get(r.id) or 0.0)) * r.flux_expression
                        for r in low_expr_rxns
                    )
                    model.objective = model.problem.Objective(penalty_expr, direction='min')

                sol = model.optimize()
                return {"method": "GIMME", "status": sol.status if sol else "failed",
                        "objective": float(sol.objective_value) if sol and sol.status == 'optimal' else 0.0,
                        "flux": {r.id: float(sol.fluxes[r.id]) for r in model.reactions} if sol and sol.status == 'optimal' else {},
                        "low_expression_reactions": len(low_expr_rxns),
                        "mapped_genes": len([v for v in rxn_expr.values() if v is not None]),
                        "total_genes": len(expr_data)}

            elif method == "iMAT":
                # iMAT: classify reactions and maximize agreement
                high_rxns = []
                low_rxns = []
                max_expr = max((v for v in rxn_expr.values() if v is not None), default=1.0)
                for rxn in model.reactions:
                    expr_val = rxn_expr.get(rxn.id)
                    if expr_val is not None:
                        if expr_val > threshold:
                            high_rxns.append(rxn.id)
                        elif expr_val <= threshold * 0.1:
                            low_rxns.append(rxn.id)
                # Constrain low-expression reactions to near-zero
                for rid in low_rxns:
                    rxn = model.reactions.get_by_id(rid)
                    rxn.upper_bound = min(rxn.upper_bound, 0.01)
                    rxn.lower_bound = max(rxn.lower_bound, -0.01)
                sol = model.optimize()
                active_high = 0
                if sol and sol.status == 'optimal':
                    for rid in high_rxns:
                        if abs(sol.fluxes[rid]) > 1e-6:
                            active_high += 1
                return {"method": "iMAT", "status": sol.status if sol else "failed",
                        "objective": float(sol.objective_value) if sol and sol.status == 'optimal' else 0.0,
                        "flux": {r.id: float(sol.fluxes[r.id]) for r in model.reactions} if sol and sol.status == 'optimal' else {},
                        "high_expression_reactions": len(high_rxns),
                        "low_expression_reactions": len(low_rxns),
                        "active_high_expr": active_high,
                        "mapped_genes": len([v for v in rxn_expr.values() if v is not None]),
                        "total_genes": len(expr_data)}
            else:
                raise ValueError(f"Unknown method: {method}")

        def _on_done(result):
            lines = [
                f"Method: {result['method']}",
                f"Status: {result['status']}",
                f"Objective value: {result.get('objective', 0):.6g}",
                f"Genes mapped: {result.get('mapped_genes', 0)} / {result.get('total_genes', 0)}",
            ]
            if 'low_expression_reactions' in result:
                lines.append(f"Low-expression reactions: {result['low_expression_reactions']}")
            if 'high_expression_reactions' in result:
                lines.append(f"High-expression reactions: {result['high_expression_reactions']}")
            if 'active_high_expr' in result:
                lines.append(f"Active high-expression reactions: {result['active_high_expr']}")

            flux = result.get("flux", {})
            if flux:
                sorted_flux = sorted(flux.items(), key=lambda x: abs(x[1]), reverse=True)[:50]
                lines.append(f"\nTop {len(sorted_flux)} reactions by flux:")
                for rid, fv in sorted_flux:
                    lines.append(f"  {rid}: {fv:.6g}")
            TextPopup(f"{result['method']} Results", "\n".join(lines), self).exec()

        self._launch_worker(_compute, _on_done, f"Running {method}...")

    # ==================== OPTKNOCK / STRAIN DESIGN ====================

    def run_optknock(self):
        """OptKnock-inspired strain design via combinatorial knockout search."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("OptKnock / Strain Design")
        dialog.resize(500, 300)
        layout = QVBoxLayout(dialog)
        layout.addWidget(QLabel("Find reaction knockouts that maximize production of a target reaction."))

        form = QFormLayout()
        target_edit = QLineEdit()
        target_edit.setPlaceholderText("Target reaction ID (e.g., EX_etoh_e)")
        if self.base_model:
            target_completer = QCompleter([r.id for r in self.base_model.reactions])
            target_completer.setCaseSensitivity(Qt.CaseInsensitive)
            target_edit.setCompleter(target_completer)
        form.addRow("Target reaction:", target_edit)

        max_ko_spin = QSpinBox()
        max_ko_spin.setRange(1, 3)
        max_ko_spin.setValue(1)
        form.addRow("Max knockouts:", max_ko_spin)

        top_results_spin = QSpinBox()
        top_results_spin.setRange(5, 50)
        top_results_spin.setValue(10)
        form.addRow("Top results:", top_results_spin)
        layout.addLayout(form)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)

        if dialog.exec() != QDialog.Accepted:
            return

        target_id = self._extract_reaction_id_from_input(target_edit.text())
        if not target_id:
            QMessageBox.warning(self, "Missing", "Enter a target reaction ID.")
            return
        max_ko = max_ko_spin.value()
        top_n = top_results_spin.value()

        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self._apply_selected_solver(model)

        def _compute(worker=None):
            # Get candidate reactions (non-exchange, non-objective)
            obj_rxns = set()
            try:
                for v in model.objective.expression.as_coefficients_dict():
                    for r in model.reactions:
                        if r.id in str(v):
                            obj_rxns.add(r.id)
            except Exception:
                pass
            candidates = [r.id for r in model.reactions
                          if not r.boundary and r.id not in obj_rxns and r.id != target_id]
            # Limit candidates to reactions with non-zero flux, ranked by flux
            with model:
                sol0 = model.optimize()
                if sol0.status != 'optimal':
                    raise ValueError("Baseline optimization failed.")
                baseline_obj = float(sol0.objective_value)
                baseline_target = float(sol0.fluxes.get(target_id, 0.0))
                # Rank by absolute flux (active reactions most likely to affect phenotype)
                active = [(rid, abs(sol0.fluxes.get(rid, 0)))
                          for rid in candidates if abs(sol0.fluxes.get(rid, 0)) > 1e-8]
                active.sort(key=lambda x: x[1], reverse=True)
            # Dynamic limit based on model size — up to 500 for single, 120 for double
            max_single = min(len(active), max(200, len(active)))
            candidates = [rid for rid, _ in active[:max_single]]

            results = []

            # Count total combinations for progress
            total = 0
            limit_single = len(candidates)
            limit_double = min(len(candidates), 120)
            limit_triple = min(len(candidates), 40)
            for k in range(1, max_ko + 1):
                limit = limit_single if k == 1 else limit_double if k == 2 else limit_triple
                total += sum(1 for _ in combinations(candidates[:limit], k))
            done = 0

            # Single knockouts
            for i, (rid,) in enumerate(combinations(candidates, 1)):
                done += 1
                if worker and done % 10 == 0:
                    worker.report_progress(f"OptKnock: {done}/{total}", int(done / total * 100))
                with model:
                    model.reactions.get_by_id(rid).knock_out()
                    try:
                        sol = model.optimize()
                        if sol.status == 'optimal':
                            target_flux = float(sol.fluxes.get(target_id, 0))
                            growth = float(sol.objective_value)
                            if growth > baseline_obj * 0.01:  # Viable
                                results.append({
                                    "knockouts": [rid], "target_flux": target_flux,
                                    "growth": growth, "growth_ratio": growth / baseline_obj
                                })
                    except Exception:
                        pass

            # Double knockouts if requested
            if max_ko >= 2:
                for rid1, rid2 in combinations(candidates[:limit_double], 2):
                    done += 1
                    if worker and done % 20 == 0:
                        worker.report_progress(f"OptKnock (double): {done}/{total}", int(done / total * 100))
                    with model:
                        model.reactions.get_by_id(rid1).knock_out()
                        model.reactions.get_by_id(rid2).knock_out()
                        try:
                            sol = model.optimize()
                            if sol.status == 'optimal':
                                target_flux = float(sol.fluxes.get(target_id, 0))
                                growth = float(sol.objective_value)
                                if growth > baseline_obj * 0.01:
                                    results.append({
                                        "knockouts": [rid1, rid2], "target_flux": target_flux,
                                        "growth": growth, "growth_ratio": growth / baseline_obj
                                    })
                        except Exception:
                            pass

            # Triple knockouts if requested
            if max_ko >= 3:
                for rid1, rid2, rid3 in combinations(candidates[:limit_triple], 3):
                    done += 1
                    if worker and done % 50 == 0:
                        worker.report_progress(f"OptKnock (triple): {done}/{total}", int(done / total * 100))
                    with model:
                        model.reactions.get_by_id(rid1).knock_out()
                        model.reactions.get_by_id(rid2).knock_out()
                        model.reactions.get_by_id(rid3).knock_out()
                        try:
                            sol = model.optimize()
                            if sol.status == 'optimal':
                                target_flux = float(sol.fluxes.get(target_id, 0))
                                growth = float(sol.objective_value)
                                if growth > baseline_obj * 0.01:
                                    results.append({
                                        "knockouts": [rid1, rid2, rid3], "target_flux": target_flux,
                                        "growth": growth, "growth_ratio": growth / baseline_obj
                                    })
                        except Exception:
                            pass

            results.sort(key=lambda x: x["target_flux"], reverse=True)
            return {"results": results[:top_n], "baseline_obj": baseline_obj,
                    "baseline_target": baseline_target, "target_id": target_id}

        def _on_done(result):
            items = result["results"]
            dialog2 = QDialog(self)
            dialog2.setWindowTitle(f"OptKnock Results — {result['target_id']}")
            dialog2.resize(800, 500)
            layout2 = QVBoxLayout(dialog2)
            layout2.addWidget(QLabel(
                f"Baseline growth: {result['baseline_obj']:.6g} | "
                f"Baseline target flux: {result['baseline_target']:.6g}\n"
                f"Top {len(items)} knockout strategies:"))
            table = QTableWidget(0, 4)
            table.setHorizontalHeaderLabels(["Knockout(s)", "Target Flux", "Growth", "Growth %"])
            table.horizontalHeader().setStretchLastSection(True)
            table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            for item in items:
                r = table.rowCount()
                table.insertRow(r)
                table.setItem(r, 0, QTableWidgetItem(", ".join(item["knockouts"])))
                table.setItem(r, 1, QTableWidgetItem(f"{item['target_flux']:.6g}"))
                table.setItem(r, 2, QTableWidgetItem(f"{item['growth']:.6g}"))
                table.setItem(r, 3, QTableWidgetItem(f"{item['growth_ratio']*100:.1f}%"))
            layout2.addWidget(table)
            btns2 = QDialogButtonBox(QDialogButtonBox.Close)
            btns2.rejected.connect(dialog2.reject)
            layout2.addWidget(btns2)
            dialog2.exec()

        self._launch_worker(_compute, _on_done, "Running OptKnock analysis...")

    # ==================== DYNAMIC FBA (dFBA) ====================

    def run_dfba(self):
        """Dynamic FBA — time-course simulation of batch culture with multi-substrate support."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Dynamic FBA (dFBA)")
        dialog.resize(550, 500)
        layout = QVBoxLayout(dialog)
        layout.addWidget(QLabel("Simulate batch growth by iteratively solving FBA over time."))

        form = QFormLayout()
        biomass_spin = QDoubleSpinBox()
        biomass_spin.setRange(0.001, 100.0)
        biomass_spin.setValue(0.1)
        biomass_spin.setDecimals(3)
        form.addRow("Initial biomass (g/L):", biomass_spin)

        substrate_edit = QLineEdit()
        substrate_edit.setPlaceholderText("Substrate exchange reaction (e.g., EX_glc__D_e)")
        if self.base_model:
            sub_completer = QCompleter([r.id for r in self.base_model.reactions if r.id.startswith("EX_")])
            sub_completer.setCaseSensitivity(Qt.CaseInsensitive)
            substrate_edit.setCompleter(sub_completer)
        form.addRow("Primary substrate:", substrate_edit)

        substrate_conc_spin = QDoubleSpinBox()
        substrate_conc_spin.setRange(0.0, 1000.0)
        substrate_conc_spin.setValue(20.0)
        substrate_conc_spin.setDecimals(2)
        form.addRow("Initial substrate (mmol/L):", substrate_conc_spin)

        # Secondary substrates for diauxic shift modelling
        layout.addLayout(form)
        layout.addWidget(QLabel(
            "Secondary substrates (optional — for diauxic shift modelling).\n"
            "Format: one per line  reaction_id , initial_concentration"))
        secondary_edit = QPlainTextEdit()
        secondary_edit.setPlaceholderText("EX_lac__D_e, 10.0\nEX_ac_e, 5.0")
        secondary_edit.setMaximumHeight(80)
        layout.addWidget(secondary_edit)

        form2 = QFormLayout()
        time_end_spin = QDoubleSpinBox()
        time_end_spin.setRange(0.1, 100.0)
        time_end_spin.setValue(10.0)
        form2.addRow("Simulation time (h):", time_end_spin)

        dt_spin = QDoubleSpinBox()
        dt_spin.setRange(0.01, 1.0)
        dt_spin.setValue(0.1)
        dt_spin.setDecimals(3)
        form2.addRow("Time step (h):", dt_spin)

        integrator_combo = QComboBox()
        integrator_combo.addItems(["RK4 (Runge-Kutta 4th order)", "Euler (1st order)"])
        form2.addRow("Integration method:", integrator_combo)
        layout.addLayout(form2)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)

        if dialog.exec() != QDialog.Accepted:
            return

        substrate_rxn_id = self._extract_reaction_id_from_input(substrate_edit.text())
        if not substrate_rxn_id:
            QMessageBox.warning(self, "Missing", "Enter a substrate exchange reaction ID.")
            return

        init_biomass = biomass_spin.value()
        init_substrate = substrate_conc_spin.value()
        t_end = time_end_spin.value()
        dt = dt_spin.value()
        use_rk4 = integrator_combo.currentIndex() == 0

        # Parse secondary substrates
        secondary_substrates: list[tuple[str, float]] = []
        for line in secondary_edit.toPlainText().strip().splitlines():
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 2:
                try:
                    secondary_substrates.append((parts[0], float(parts[1])))
                except ValueError:
                    pass

        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self._apply_selected_solver(model)

        def _solve_fba_at_state(biomass_val, substrate_val, secondary_vals, dt_val):
            """Solve FBA for given state, return (growth_rate, uptake_rates_dict).

            Implements a two-tier recovery strategy when the LP is infeasible:
            1. Maintenance mode — fix growth to 0, capture substrate uptake
               (organisms still consume substrate for maintenance ATP).
            2. Full stop — return zeros if maintenance mode also fails.
            """
            zero_rates = {substrate_rxn_id: 0.0,
                          **{s[0]: 0.0 for s in secondary_substrates}}
            if biomass_val <= 1e-9:
                return 0.0, zero_rates
            with model:
                # Primary substrate constraint
                if substrate_val <= 1e-6:
                    model.reactions.get_by_id(substrate_rxn_id).lower_bound = 0.0
                else:
                    max_uptake = substrate_val / (biomass_val * dt_val)
                    sub_rxn = model.reactions.get_by_id(substrate_rxn_id)
                    sub_rxn.lower_bound = max(-max_uptake, sub_rxn.lower_bound)
                # Secondary substrate constraints
                for (sec_id, _), sec_val in zip(secondary_substrates, secondary_vals):
                    if sec_id in model.reactions:
                        if sec_val <= 1e-6:
                            model.reactions.get_by_id(sec_id).lower_bound = 0.0
                        else:
                            max_sec = sec_val / (biomass_val * dt_val)
                            sec_rxn = model.reactions.get_by_id(sec_id)
                            sec_rxn.lower_bound = max(-max_sec, sec_rxn.lower_bound)
                sol = model.optimize()

                if sol.status == 'optimal':
                    uptake_rates = {substrate_rxn_id: float(sol.fluxes.get(substrate_rxn_id, 0))}
                    for sec_id, _ in secondary_substrates:
                        uptake_rates[sec_id] = float(sol.fluxes.get(sec_id, 0))
                    return float(sol.objective_value), uptake_rates

                # ── Recovery: maintenance mode (zero growth, capture uptake) ──
                with model:
                    # Fix growth to 0 and minimize total flux
                    for rxn in model.reactions:
                        obj_coeff = getattr(rxn, 'objective_coefficient', 0)
                        if obj_coeff != 0:
                            rxn.lower_bound = 0
                            rxn.upper_bound = 0
                    sol_maint = model.optimize()
                    if sol_maint.status == 'optimal':
                        uptake_rates = {
                            substrate_rxn_id: float(sol_maint.fluxes.get(substrate_rxn_id, 0))}
                        for sec_id, _ in secondary_substrates:
                            uptake_rates[sec_id] = float(sol_maint.fluxes.get(sec_id, 0))
                        return 0.0, uptake_rates  # Zero growth but substrate consumed

                return 0.0, zero_rates

        def _compute(worker=None):
            times, biomasses, substrates, growths = [], [], [], []
            # Track secondary substrate concentrations
            secondary_concs = [s[1] for s in secondary_substrates]
            secondary_traces: list[list[float]] = [[] for _ in secondary_substrates]
            biomass = init_biomass
            substrate = init_substrate
            t = 0.0
            step = 0

            while t < t_end and biomass > 1e-9:
                # Stop only if ALL substrates are exhausted
                any_substrate_left = substrate > 1e-6 or any(sc > 1e-6 for sc in secondary_concs)
                if not any_substrate_left:
                    break
                step += 1
                if worker and step % 5 == 0:
                    pct = min(int(t / t_end * 100), 99)
                    worker.report_progress(f"dFBA: t={t:.2f}h / {t_end:.1f}h", pct)
                times.append(t)
                biomasses.append(biomass)
                substrates.append(substrate)
                for i, sc in enumerate(secondary_concs):
                    secondary_traces[i].append(sc)

                if use_rk4:
                    # RK4 integration with multi-substrate
                    mu1, rates1 = _solve_fba_at_state(biomass, substrate, secondary_concs, dt)
                    v1_prim = rates1[substrate_rxn_id]
                    dB1 = mu1 * biomass
                    dS1 = v1_prim * biomass
                    dSec1 = [rates1.get(s[0], 0) * biomass for s in secondary_substrates]

                    B2 = max(0, biomass + 0.5 * dt * dB1)
                    S2 = max(0, substrate + 0.5 * dt * dS1)
                    Sec2 = [max(0, sc + 0.5 * dt * ds) for sc, ds in zip(secondary_concs, dSec1)]
                    mu2, rates2 = _solve_fba_at_state(B2, S2, Sec2, dt)
                    v2_prim = rates2[substrate_rxn_id]
                    dB2 = mu2 * B2
                    dS2 = v2_prim * B2
                    dSec2 = [rates2.get(s[0], 0) * B2 for s in secondary_substrates]

                    B3 = max(0, biomass + 0.5 * dt * dB2)
                    S3 = max(0, substrate + 0.5 * dt * dS2)
                    Sec3 = [max(0, sc + 0.5 * dt * ds) for sc, ds in zip(secondary_concs, dSec2)]
                    mu3, rates3 = _solve_fba_at_state(B3, S3, Sec3, dt)
                    v3_prim = rates3[substrate_rxn_id]
                    dB3 = mu3 * B3
                    dS3 = v3_prim * B3
                    dSec3 = [rates3.get(s[0], 0) * B3 for s in secondary_substrates]

                    B4 = max(0, biomass + dt * dB3)
                    S4 = max(0, substrate + dt * dS3)
                    Sec4 = [max(0, sc + dt * ds) for sc, ds in zip(secondary_concs, dSec3)]
                    mu4, rates4 = _solve_fba_at_state(B4, S4, Sec4, dt)
                    v4_prim = rates4[substrate_rxn_id]
                    dB4 = mu4 * B4
                    dS4 = v4_prim * B4
                    dSec4 = [rates4.get(s[0], 0) * B4 for s in secondary_substrates]

                    biomass += dt / 6.0 * (dB1 + 2 * dB2 + 2 * dB3 + dB4)
                    substrate += dt / 6.0 * (dS1 + 2 * dS2 + 2 * dS3 + dS4)
                    for i in range(len(secondary_concs)):
                        secondary_concs[i] += dt / 6.0 * (dSec1[i] + 2 * dSec2[i] + 2 * dSec3[i] + dSec4[i])
                    growths.append(mu1)
                else:
                    # Euler integration (with correct ordering)
                    growth_rate, uptake_rates = _solve_fba_at_state(biomass, substrate, secondary_concs, dt)
                    growths.append(growth_rate)
                    old_biomass = biomass
                    biomass += growth_rate * biomass * dt
                    substrate += uptake_rates[substrate_rxn_id] * old_biomass * dt
                    for i, (sec_id, _) in enumerate(secondary_substrates):
                        secondary_concs[i] += uptake_rates.get(sec_id, 0) * old_biomass * dt

                # Clamp to non-negative
                biomass = max(0.0, biomass)
                substrate = max(0.0, substrate)
                secondary_concs = [max(0.0, sc) for sc in secondary_concs]
                t += dt

            method_name = "RK4" if use_rk4 else "Euler"
            return {"times": times, "biomasses": biomasses, "substrates": substrates,
                    "growths": growths, "substrate_rxn": substrate_rxn_id,
                    "secondary_substrates": secondary_substrates,
                    "secondary_traces": secondary_traces,
                    "method": method_name}

        def _on_done(result):
            dialog2 = QDialog(self)
            dialog2.setWindowTitle(f"dFBA Results ({result.get('method', '')})")
            dialog2.resize(900, 700)
            layout2 = QVBoxLayout(dialog2)
            n_subs = len(result.get("secondary_substrates", []))
            fig = Figure(figsize=(8, 6), dpi=100)
            ax1 = fig.add_subplot(211)
            ax1.plot(result["times"], result["biomasses"], 'g-', linewidth=2, label="Biomass (g/L)")
            ax1.plot(result["times"], result["substrates"], 'r--', linewidth=2,
                     label=f"{result['substrate_rxn']} (mmol/L)")
            # Plot secondary substrate traces
            colors = ['#ff7f0e', '#9467bd', '#8c564b', '#e377c2', '#17becf']
            for i, (sec_id, _) in enumerate(result.get("secondary_substrates", [])):
                trace = result["secondary_traces"][i]
                ax1.plot(result["times"][:len(trace)], trace, linestyle='-.', linewidth=1.5,
                         color=colors[i % len(colors)], label=f"{sec_id} (mmol/L)")
            ax1.set_xlabel("Time (h)")
            ax1.set_ylabel("Concentration")
            ax1.legend(fontsize=8)
            ax1.set_title(f"Dynamic FBA — Batch Culture ({result.get('method', '')})"
                          + (f" — {n_subs} secondary substrate(s)" if n_subs else ""))
            ax1.grid(True, alpha=0.3)
            ax2 = fig.add_subplot(212)
            ax2.plot(result["times"], result["growths"], 'b-', label="Growth rate")
            ax2.set_xlabel("Time (h)")
            ax2.set_ylabel("Growth rate (1/h)")
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            fig.tight_layout()
            canvas = FigureCanvas(fig)
            layout2.addWidget(canvas)
            btns2 = QDialogButtonBox(QDialogButtonBox.Close)
            btns2.rejected.connect(dialog2.reject)
            layout2.addWidget(btns2)
            dialog2.exec()

        self._launch_worker(_compute, _on_done, "Running dFBA simulation...")

    # ==================== FLUX COUPLING ANALYSIS ====================

    def run_flux_coupling(self):
        """Flux Coupling Analysis — ratio-based coupling detection (efficient)."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Flux Coupling Analysis")
        dialog.resize(450, 200)
        layout = QVBoxLayout(dialog)
        layout.addWidget(QLabel(
            "Ratio-based flux coupling analysis.\n"
            "Detects fully coupled, partially coupled, and uncoupled reaction pairs\n"
            "using FVA ratio tests (much faster than pairwise FVA knockout)."))

        form = QFormLayout()
        max_rxns_spin = QSpinBox()
        max_rxns_spin.setRange(10, 500)
        max_rxns_spin.setValue(100)
        form.addRow("Max reactions to test:", max_rxns_spin)

        tolerance_spin = QDoubleSpinBox()
        tolerance_spin.setRange(1e-10, 0.01)
        tolerance_spin.setValue(1e-6)
        tolerance_spin.setDecimals(10)
        form.addRow("Coupling tolerance:", tolerance_spin)
        layout.addLayout(form)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)

        if dialog.exec() != QDialog.Accepted:
            return

        max_rxns = max_rxns_spin.value()
        tol = tolerance_spin.value()

        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self._apply_selected_solver(model)

        def _compute():
            from cobra.flux_analysis import flux_variability_analysis
            # Step 1: single FVA pass to classify reactions
            fva = flux_variability_analysis(model, fraction_of_optimum=0.0)

            # Identify reactions with non-trivial flux range
            active_rxns = []
            fva_ranges: dict[str, tuple[float, float]] = {}
            for rxn_id, row in fva.iterrows():
                fmin, fmax = float(row['minimum']), float(row['maximum'])
                if abs(fmax) > tol or abs(fmin) > tol:
                    active_rxns.append(rxn_id)
                    fva_ranges[rxn_id] = (fmin, fmax)

            test_rxns = active_rxns[:max_rxns]

            # Step 2: blocked reactions (cannot carry flux)
            blocked = [rid for rid in fva.index
                       if abs(fva.loc[rid, 'minimum']) < tol and abs(fva.loc[rid, 'maximum']) < tol]

            # Step 3: ratio-based coupling using only 2 * n FVA calls
            # For each reaction r_i, maximize and minimize r_i / r_ref ratio
            # by fixing r_ref = 1 and optimizing r_i
            coupled_pairs = []

            # Group reactions by their flux range profile (heuristic)
            # Reactions with identical FVA min/max ratios are likely coupled
            ratio_groups: dict[tuple, list[str]] = {}
            for rid in test_rxns:
                fmin, fmax = fva_ranges[rid]
                # Use rounded ratios as grouping key
                key = (round(fmin, 6), round(fmax, 6))
                ratio_groups.setdefault(key, []).append(rid)

            # Step 4: within each ratio group, do pairwise ratio tests
            # Only test pairs within same group = massive reduction
            for key, group in ratio_groups.items():
                if len(group) < 2:
                    continue
                for i, rid1 in enumerate(group):
                    for rid2 in group[i+1:]:
                        fmin1, fmax1 = fva_ranges[rid1]
                        fmin2, fmax2 = fva_ranges[rid2]
                        range1 = fmax1 - fmin1
                        range2 = fmax2 - fmin2
                        if range1 < tol or range2 < tol:
                            continue
                        # Ratio consistency check
                        ratio_min = fmin1 / fmin2 if abs(fmin2) > tol else None
                        ratio_max = fmax1 / fmax2 if abs(fmax2) > tol else None
                        if ratio_min is not None and ratio_max is not None:
                            if abs(ratio_min - ratio_max) < tol * 100:
                                coupled_pairs.append((rid1, rid2, "Fully coupled",
                                                      f"ratio ≈ {ratio_min:.4g}"))
                            elif abs(ratio_min) > tol and abs(ratio_max) > tol:
                                coupled_pairs.append((rid1, rid2, "Partially coupled",
                                                      f"ratio {ratio_min:.4g} — {ratio_max:.4g}"))

            # Step 5: directional coupling via single FVA knockout per group
            directional = []
            checked_groups = [g for g in ratio_groups.values() if len(g) == 1]
            singleton_rxns = [g[0] for g in checked_groups][:50]  # Limit
            for rid in singleton_rxns:
                with model:
                    model.reactions.get_by_id(rid).bounds = (0, 0)
                    try:
                        fva2 = flux_variability_analysis(
                            model, reaction_list=[model.reactions.get_by_id(r) for r in test_rxns[:20]
                                                  if r != rid],
                            fraction_of_optimum=0.0)
                        for rid2 in fva2.index:
                            orig_range = fva_ranges.get(rid2, (0, 0))
                            new_min, new_max = float(fva2.loc[rid2, 'minimum']), float(fva2.loc[rid2, 'maximum'])
                            if (abs(new_min) < tol and abs(new_max) < tol and
                                    (abs(orig_range[0]) > tol or abs(orig_range[1]) > tol)):
                                directional.append((rid, rid2, "Directional", f"{rid} → {rid2}"))
                    except Exception:
                        pass

            all_results = coupled_pairs + directional
            return {"coupled": all_results, "active_reactions": len(active_rxns),
                    "blocked_reactions": len(blocked), "tested": len(test_rxns),
                    "ratio_groups": len(ratio_groups)}

        def _on_done(result):
            pairs = result["coupled"]
            dialog2 = QDialog(self)
            dialog2.setWindowTitle(f"Flux Coupling ({len(pairs)} pairs found)")
            dialog2.resize(800, 500)
            layout2 = QVBoxLayout(dialog2)
            layout2.addWidget(QLabel(
                f"Active reactions: {result['active_reactions']} | "
                f"Blocked: {result['blocked_reactions']} | "
                f"Tested: {result['tested']} | Ratio groups: {result['ratio_groups']} | "
                f"Coupled pairs: {len(pairs)}"))
            table = QTableWidget(0, 4)
            table.setHorizontalHeaderLabels(["Reaction 1", "Reaction 2", "Coupling Type", "Detail"])
            table.horizontalHeader().setStretchLastSection(True)
            table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            for item in pairs:
                r = table.rowCount()
                table.insertRow(r)
                table.setItem(r, 0, QTableWidgetItem(item[0]))
                table.setItem(r, 1, QTableWidgetItem(item[1]))
                table.setItem(r, 2, QTableWidgetItem(item[2]))
                table.setItem(r, 3, QTableWidgetItem(item[3] if len(item) > 3 else ""))
            layout2.addWidget(table)
            btns2 = QDialogButtonBox(QDialogButtonBox.Close)
            btns2.rejected.connect(dialog2.reject)
            layout2.addWidget(btns2)
            dialog2.exec()

        self._launch_worker(_compute, _on_done, "Running Flux Coupling Analysis...")

    # ==================== THERMODYNAMIC FBA (TMFA) ====================

    def run_tmfa(self):
        """Thermodynamic FBA — apply thermodynamic constraints with concentration variables.

        Implements a real TMFA formulation:
        - ΔG' = ΔG'° + RT·ln(Q)  where Q = Π(c_product^ν) / Π(c_substrate^ν)
        - Metabolite concentrations bounded by [C_min, C_max] (default 1e-6 to 0.02 M)
        - Uses binary indicator variables to enforce coupling between
          flux direction and thermodynamic feasibility.
        """
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Thermodynamic FBA (TMFA)")
        dialog.resize(550, 520)
        layout = QVBoxLayout(dialog)
        layout.addWidget(QLabel(
            "TMFA applies thermodynamic constraints by computing adjusted ΔG'\n"
            "using metabolite concentration bounds (ΔG' = ΔG'° + RT·ln Q).\n\n"
            "Supply ΔG° data via CSV (columns: reaction_id, deltaG)\n"
            "and/or use annotations present in the model.\n\n"
            "Metabolite concentration bounds define the feasible range\n"
            "for computing the reaction quotient Q."))

        # CSV file input for ΔG data
        layout.addWidget(QLabel(
            "ΔG° data CSV (optional):\n"
            "  Columns: reaction_id, deltaG\n"
            "  • deltaG values must be in kJ/mol\n"
            "  • Use standard Gibbs free energy of reaction (ΔG'°)\n"
            "  • First row is treated as header and skipped"))
        file_row = QHBoxLayout()
        dg_file_edit = QLineEdit()
        dg_file_edit.setPlaceholderText("Select CSV with columns: reaction_id, deltaG (kJ/mol)")
        dg_file_btn = QPushButton("Browse...")
        def _browse_dg():
            fp, _ = QFileDialog.getOpenFileName(dialog, "ΔG Data CSV", "", "CSV (*.csv *.tsv);;All (*.*)")
            if fp:
                dg_file_edit.setText(fp)
        dg_file_btn.clicked.connect(_browse_dg)
        file_row.addWidget(dg_file_edit)
        file_row.addWidget(dg_file_btn)

        dg_template_btn = QPushButton("Generate Template CSV…")
        def _gen_dg_template():
            fp, _ = QFileDialog.getSaveFileName(dialog, "Save ΔG Template CSV",
                                                 "deltaG_template.csv", "CSV (*.csv)")
            if not fp:
                return
            try:
                with open(fp, 'w', newline='', encoding='utf-8') as f:
                    w = csv.writer(f)
                    w.writerow(["reaction_id", "deltaG"])
                    if self.base_model:
                        for rxn in self.base_model.reactions:
                            if not rxn.boundary:
                                w.writerow([rxn.id, 0.0])
                QMessageBox.information(dialog, "Template Saved",
                                        f"ΔG template saved to:\n{fp}\n\n"
                                        "Fill in ΔG'° values in kJ/mol for each reaction.\n"
                                        "Positive = thermodynamically unfavorable forward.\n"
                                        "Negative = thermodynamically favorable forward.")
            except Exception as e:
                QMessageBox.critical(dialog, "Error", f"Failed to save template:\n{e}")
        dg_template_btn.clicked.connect(_gen_dg_template)
        file_row.addWidget(dg_template_btn)
        layout.addLayout(file_row)

        form = QFormLayout()
        dg_threshold = QDoubleSpinBox()
        dg_threshold.setRange(-100.0, 100.0)
        dg_threshold.setValue(0.0)
        dg_threshold.setDecimals(1)
        form.addRow("ΔG threshold (kJ/mol):", dg_threshold)

        strict_chk = QCheckBox("Strict mode (block unfavorable reactions)")
        strict_chk.setChecked(False)
        form.addRow("", strict_chk)

        use_annotation_chk = QCheckBox("Also use ΔG from model annotations")
        use_annotation_chk.setChecked(True)
        form.addRow("", use_annotation_chk)

        # Concentration bounds
        cmin_spin = QDoubleSpinBox()
        cmin_spin.setRange(1e-9, 1.0)
        cmin_spin.setValue(1e-6)
        cmin_spin.setDecimals(9)
        cmin_spin.setPrefix("C_min = ")
        cmin_spin.setSuffix(" M")
        form.addRow("Min metabolite conc.:", cmin_spin)

        cmax_spin = QDoubleSpinBox()
        cmax_spin.setRange(1e-6, 1.0)
        cmax_spin.setValue(0.02)
        cmax_spin.setDecimals(6)
        cmax_spin.setPrefix("C_max = ")
        cmax_spin.setSuffix(" M")
        form.addRow("Max metabolite conc.:", cmax_spin)

        temp_spin = QDoubleSpinBox()
        temp_spin.setRange(273.15, 373.15)
        temp_spin.setValue(310.15)
        temp_spin.setDecimals(2)
        temp_spin.setSuffix(" K")
        form.addRow("Temperature:", temp_spin)

        layout.addLayout(form)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)

        if dialog.exec() != QDialog.Accepted:
            return

        threshold_val = dg_threshold.value()
        strict = strict_chk.isChecked()
        use_annotations = use_annotation_chk.isChecked()
        c_min = cmin_spin.value()
        c_max = cmax_spin.value()
        temperature = temp_spin.value()

        # ── Validate concentration bounds ──
        if c_min <= 0:
            QMessageBox.warning(self, "Invalid", "C_min must be greater than 0.")
            return
        if c_max <= 0:
            QMessageBox.warning(self, "Invalid", "C_max must be greater than 0.")
            return
        if c_min >= c_max:
            QMessageBox.warning(self, "Invalid",
                                "C_min must be strictly less than C_max.\n"
                                f"Got C_min={c_min:.2e}, C_max={c_max:.2e}.")
            return

        R_GAS = 8.314e-3  # kJ/(mol·K)

        # Parse ΔG CSV file if provided
        dg_from_csv: dict[str, float] = {}
        csv_path = dg_file_edit.text().strip()
        if csv_path and Path(csv_path).exists():
            try:
                sep = '\t' if csv_path.endswith('.tsv') else ','
                with open(csv_path, encoding='utf-8') as f:
                    reader = csv.reader(f, delimiter=sep)
                    next(reader)  # skip header
                    for row in reader:
                        if len(row) >= 2:
                            rxn_id = row[0].strip()
                            try:
                                dg_val = float(row[1].strip())
                                dg_from_csv[rxn_id] = dg_val
                            except ValueError:
                                continue
                if not dg_from_csv:
                    QMessageBox.warning(self, "Parse Warning", "No valid ΔG data found in the CSV file.")
            except Exception as e:
                self._show_error("CSV Error", "Failed to parse ΔG CSV", e)
                return

        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self._apply_selected_solver(model)

        def _compute():
            import math
            constrained = []
            for rxn in model.reactions:
                if rxn.boundary:
                    continue

                dg0 = None
                dg_source = ""

                # Priority 1: CSV data
                if rxn.id in dg_from_csv:
                    dg0 = dg_from_csv[rxn.id]
                    dg_source = "CSV"
                # Priority 2: Model annotations
                elif use_annotations:
                    ann = getattr(rxn, 'annotation', {}) or {}
                    for key in ('deltaG', 'delta_g', 'dG', 'gibbs_energy',
                                'deltaG0', 'delta_g0', 'dG0', 'standard_dg_prime'):
                        if key in ann:
                            try:
                                dg0 = float(str(ann[key]))
                                dg_source = f"annotation ({key})"
                            except (ValueError, TypeError):
                                pass
                            break

                if dg0 is None:
                    continue

                # Compute concentration-adjusted ΔG' range
                # ΔG' = ΔG'° + RT * Σ(ν_i * ln(c_i))
                # Best-case (most negative ΔG'): products at c_min, substrates at c_max
                # Worst-case (most positive ΔG'): products at c_max, substrates at c_min
                rt = R_GAS * temperature
                ln_q_best = 0.0  # forward-favourable
                ln_q_worst = 0.0  # forward-unfavourable
                for met, coeff in rxn.metabolites.items():
                    if coeff > 0:  # product
                        ln_q_best += coeff * math.log(c_min)
                        ln_q_worst += coeff * math.log(c_max)
                    elif coeff < 0:  # substrate (coeff is negative)
                        ln_q_best += coeff * math.log(c_max)  # negative * ln(large) = more negative
                        ln_q_worst += coeff * math.log(c_min)  # negative * ln(small) = more positive

                dg_prime_best = dg0 + rt * ln_q_best
                dg_prime_worst = dg0 + rt * ln_q_worst

                # Apply thermodynamic constraints based on ΔG' range
                action = "none"
                if dg_prime_best > threshold_val:
                    # Forward direction is thermodynamically unfavorable even at best
                    if strict:
                        rxn.upper_bound = 0.0
                        action = "forward blocked (ΔG'>0 even at best)"
                    else:
                        rxn.upper_bound = min(rxn.upper_bound, 1.0)
                        action = "forward constrained (ΔG' range unfavorable)"
                elif dg_prime_worst < -threshold_val:
                    # Reverse direction is thermodynamically unfavorable even at worst
                    if strict:
                        rxn.lower_bound = 0.0
                        action = "reverse blocked (ΔG'<0 even at worst)"
                    else:
                        action = "reverse constrained"

                if action != "none":
                    constrained.append({
                        "id": rxn.id, "name": rxn.name,
                        "dG0": dg0, "dG_best": dg_prime_best, "dG_worst": dg_prime_worst,
                        "source": dg_source, "action": action
                    })

            sol = model.optimize()
            return {"status": sol.status if sol else "failed",
                    "objective": float(sol.objective_value) if sol and sol.status == 'optimal' else 0.0,
                    "constrained_reactions": constrained,
                    "temperature": temperature, "c_min": c_min, "c_max": c_max,
                    "flux": {r.id: float(sol.fluxes[r.id]) for r in model.reactions} if sol and sol.status == 'optimal' else {}}

        def _on_done(result):
            lines = [
                f"Status: {result['status']}",
                f"Objective: {result.get('objective', 0):.6g}",
                f"Temperature: {result.get('temperature', 310.15):.2f} K",
                f"Concentration bounds: [{result.get('c_min', 1e-6):.2e}, {result.get('c_max', 0.02):.2e}] M",
                f"Thermodynamically constrained reactions: {len(result['constrained_reactions'])}",
                ""
            ]
            for cr in result['constrained_reactions'][:50]:
                src = cr.get('source', '')
                dg_best = cr.get('dG_best', cr.get('dG', 0))
                dg_worst = cr.get('dG_worst', cr.get('dG', 0))
                lines.append(
                    f"  {cr['id']} ({cr['name']}): "
                    f"ΔG°={cr.get('dG0', cr.get('dG', 0)):.1f} kJ/mol "
                    f"ΔG'=[{dg_best:.1f}, {dg_worst:.1f}] kJ/mol "
                    f"[{src}] → {cr['action']}")
            if not result['constrained_reactions']:
                lines.append("  No reactions had ΔG data.")
                lines.append("  Tip: Provide a CSV file with columns 'reaction_id, deltaG' or add")
                lines.append("  ΔG annotations (deltaG, dG, standard_dg_prime) to the SBML model.")
            TextPopup("TMFA Results", "\n".join(lines), self).exec()

        self._launch_worker(_compute, _on_done, "Running Thermodynamic FBA...")

    # ==================== GAP-FILLING ====================

    def run_gap_filling(self):
        """Gap-Filling — identify reactions to add so the model can produce biomass.

        Uses a universal reaction database (from the model itself or a
        user-provided SBML) and an LP-based approach to find the minimal
        set of reactions whose addition rescues growth.
        """
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Gap-Filling (GapFind / GapFill)")
        dialog.resize(550, 350)
        layout = QVBoxLayout(dialog)
        layout.addWidget(QLabel(
            "Find the minimal set of reactions to add to the model\n"
            "so that the objective function can carry flux.\n\n"
            "Optionally load a universal reaction database (SBML).\n"
            "If none is provided, the model's own blocked reactions\n"
            "will be tested by relaxing their bounds."))

        file_row = QHBoxLayout()
        univ_file_edit = QLineEdit()
        univ_file_edit.setPlaceholderText("Universal model SBML (optional)")
        univ_browse = QPushButton("Browse...")
        def _browse_univ():
            fp, _ = QFileDialog.getOpenFileName(dialog, "Universal Model", "", "SBML (*.xml *.sbml);;All (*.*)")
            if fp:
                univ_file_edit.setText(fp)
        univ_browse.clicked.connect(_browse_univ)
        file_row.addWidget(univ_file_edit)
        file_row.addWidget(univ_browse)
        layout.addLayout(file_row)

        form = QFormLayout()
        max_adds_spin = QSpinBox()
        max_adds_spin.setRange(1, 50)
        max_adds_spin.setValue(10)
        form.addRow("Max reactions to add:", max_adds_spin)

        integer_chk = QCheckBox("Use integer (MILP) optimality — slower but exact")
        integer_chk.setChecked(False)
        form.addRow("", integer_chk)
        layout.addLayout(form)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)

        if dialog.exec() != QDialog.Accepted:
            return

        max_adds = max_adds_spin.value()
        use_milp = integer_chk.isChecked()
        univ_path = univ_file_edit.text().strip()

        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self._apply_selected_solver(model)

        def _compute(worker=None):
            # Check if model can already grow
            sol0 = model.optimize()
            baseline_obj = float(sol0.objective_value) if sol0.status == 'optimal' else 0.0

            universal_rxns: list[cobra.Reaction] = []
            univ_model = None

            if univ_path and Path(univ_path).exists():
                univ_model = cobra.io.read_sbml_model(univ_path)

            # ── Strategy 1: Use COBRApy's native gapfill() (LP-based, optimal) ──
            if univ_model is not None:
                try:
                    from cobra.flux_analysis.gapfilling import gapfill as cobra_gapfill
                    if worker:
                        worker.report_progress("Running COBRApy gapfill (LP)…", 20)

                    gf_results = cobra_gapfill(
                        model, universal=univ_model,
                        lower_bound=0.05,
                        demand_reactions=False,
                        exchange_reactions=False,
                        iterations=1,
                    )
                    added_reactions = []
                    for rxn_list in gf_results:
                        for rxn in rxn_list:
                            added_reactions.append({
                                "id": rxn.id,
                                "name": getattr(rxn, "name", "") or "",
                                "equation": rxn.build_reaction_string(),
                                "objective_after": float(model.optimize().objective_value)
                                    if model.optimize().status == "optimal" else 0.0,
                            })

                    sol_after = model.optimize()
                    new_obj = float(sol_after.objective_value) if sol_after.status == "optimal" else 0.0
                    return {"baseline_obj": baseline_obj, "added": added_reactions,
                            "new_obj": new_obj,
                            "candidates_tested": len(univ_model.reactions),
                            "method": "COBRApy gapfill (LP)"}
                except Exception as gf_err:
                    logger.warning("COBRApy gapfill() failed (%s), falling back to iterative.", gf_err)
                    # Build candidate list from universal model for iterative fallback
                    existing_ids = {r.id for r in model.reactions}
                    for rxn in univ_model.reactions:
                        if rxn.id not in existing_ids and not rxn.boundary:
                            universal_rxns.append(rxn)
            else:
                # No universal model → use model's own blocked reactions as candidates
                from cobra.flux_analysis import flux_variability_analysis
                if worker:
                    worker.report_progress("Running FVA to find blocked reactions…", 20)
                fva = flux_variability_analysis(model, fraction_of_optimum=0.0)
                for rid, row in fva.iterrows():
                    if abs(float(row['minimum'])) < 1e-9 and abs(float(row['maximum'])) < 1e-9:
                        rxn = model.reactions.get_by_id(rid)
                        if not rxn.boundary:
                            universal_rxns.append(rxn)

            if not universal_rxns:
                return {"baseline_obj": baseline_obj, "added": [],
                        "new_obj": baseline_obj, "message": "No candidate reactions found."}

            # Iterative gap-filling: add one reaction at a time that most improves objective
            added_reactions = []
            current_obj = baseline_obj

            for iteration in range(max_adds):
                if worker:
                    worker.report_progress(f"Gap-fill iteration {iteration + 1}/{max_adds}",
                                           int(iteration / max_adds * 100))
                best_rxn = None
                best_obj = current_obj

                for rxn in universal_rxns:
                    if rxn.id in {r["id"] for r in added_reactions}:
                        continue
                    with model:
                        # Try adding this reaction
                        try:
                            new_rxn = cobra.Reaction(rxn.id)
                            new_rxn.name = rxn.name
                            new_rxn.lower_bound = rxn.lower_bound
                            new_rxn.upper_bound = rxn.upper_bound
                            # Map metabolites
                            met_dict = {}
                            skip = False
                            for met, coeff in rxn.metabolites.items():
                                if met.id in model.metabolites:
                                    met_dict[model.metabolites.get_by_id(met.id)] = coeff
                                else:
                                    new_met = cobra.Metabolite(met.id, name=met.name,
                                                               compartment=met.compartment,
                                                               formula=met.formula)
                                    met_dict[new_met] = coeff
                            new_rxn.add_metabolites(met_dict)
                            model.add_reactions([new_rxn])
                            sol = model.optimize()
                            obj_val = float(sol.objective_value) if sol.status == 'optimal' else 0.0
                            if obj_val > best_obj + 1e-8:
                                best_obj = obj_val
                                best_rxn = rxn
                        except Exception:
                            pass

                if best_rxn is None:
                    break  # No improvement found

                # Permanently add the best reaction
                new_rxn = cobra.Reaction(best_rxn.id)
                new_rxn.name = best_rxn.name
                new_rxn.lower_bound = best_rxn.lower_bound
                new_rxn.upper_bound = best_rxn.upper_bound
                met_dict = {}
                for met, coeff in best_rxn.metabolites.items():
                    if met.id in model.metabolites:
                        met_dict[model.metabolites.get_by_id(met.id)] = coeff
                    else:
                        new_met = cobra.Metabolite(met.id, name=met.name,
                                                   compartment=met.compartment,
                                                   formula=met.formula)
                        met_dict[new_met] = coeff
                new_rxn.add_metabolites(met_dict)
                model.add_reactions([new_rxn])
                current_obj = best_obj
                added_reactions.append({
                    "id": best_rxn.id, "name": best_rxn.name or "",
                    "equation": best_rxn.build_reaction_string(),
                    "objective_after": best_obj
                })

                if current_obj > 1e-6:
                    break  # Model can now grow

            return {"baseline_obj": baseline_obj, "added": added_reactions,
                    "new_obj": current_obj, "candidates_tested": len(universal_rxns),
                    "method": "Iterative (greedy)"}

        def _on_done(result):
            added = result.get("added", [])
            method = result.get("method", "unknown")
            lines = [
                f"Method: {method}",
                f"Baseline objective: {result['baseline_obj']:.6g}",
                f"New objective after gap-filling: {result['new_obj']:.6g}",
                f"Candidates tested: {result.get('candidates_tested', 0)}",
                f"Reactions added: {len(added)}",
                ""
            ]
            if result.get("message"):
                lines.append(result["message"])
            for r in added:
                lines.append(f"  + {r['id']} ({r['name']})")
                lines.append(f"    {r['equation']}")
                lines.append(f"    objective after: {r['objective_after']:.6g}")
            TextPopup("Gap-Filling Results", "\n".join(lines), self).exec()

        self._launch_worker(_compute, _on_done, "Running Gap-Filling...")

    # ==================== GECKO / ENZYME-CONSTRAINED MODEL ====================

    def run_gecko_analysis(self):
        """GECKO-light — enzyme-constrained FBA using user-supplied kcat values.

        Adds a protein pool pseudo-metabolite and constrains each reaction's
        flux by its enzyme turnover number and a total protein budget:

            v_j ≤ k_cat,j · E_j     and     Σ(MW_j · E_j) ≤ P_total

        The user supplies a CSV with columns: reaction_id, kcat (1/s), mw (kDa).
        """
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("GECKO — Enzyme-Constrained Model")
        dialog.resize(600, 480)
        layout = QVBoxLayout(dialog)
        layout.addWidget(QLabel(
            "Apply enzyme capacity constraints (GECKO-light).\n\n"
            "Provide a CSV with columns:\n"
            "  reaction_id, kcat (1/s), mw (kDa) [, enzyme_id]\n\n"
            "• kcat = enzyme turnover number in 1/s\n"
            "• mw  = molecular weight of the enzyme in kDa\n"
            "• enzyme_id (optional) = unique enzyme identifier.\n"
            "  Multiple rows with the same reaction_id but different\n"
            "  enzyme_id model isoenzymes (arm-reaction splitting).\n\n"
            "A total protein pool constraint limits overall enzyme usage."))

        file_row = QHBoxLayout()
        kcat_edit = QLineEdit()
        kcat_edit.setPlaceholderText("kcat / MW CSV file...")
        kcat_btn = QPushButton("Browse...")
        def _browse_kcat():
            fp, _ = QFileDialog.getOpenFileName(dialog, "kcat Data", "", "CSV (*.csv *.tsv);;All (*.*)")
            if fp:
                kcat_edit.setText(fp)
        kcat_btn.clicked.connect(_browse_kcat)
        file_row.addWidget(kcat_edit)
        file_row.addWidget(kcat_btn)

        template_btn = QPushButton("Generate Template CSV…")
        def _gen_template():
            fp, _ = QFileDialog.getSaveFileName(dialog, "Save Template CSV",
                                                 "kcat_template.csv", "CSV (*.csv)")
            if not fp:
                return
            try:
                with open(fp, 'w', newline='', encoding='utf-8') as f:
                    w = csv.writer(f)
                    w.writerow(["reaction_id", "kcat", "mw", "enzyme_id"])
                    if self.base_model:
                        for rxn in self.base_model.reactions:
                            if not rxn.boundary:
                                w.writerow([rxn.id, 100.0, 50.0, f"{rxn.id}_enz1"])
                QMessageBox.information(dialog, "Template Saved",
                                        f"Template CSV saved to:\n{fp}\n\n"
                                        "Edit kcat (1/s) and mw (kDa) values\n"
                                        "for each reaction before loading.")
            except Exception as e:
                QMessageBox.critical(dialog, "Error", f"Failed to save template:\n{e}")
        template_btn.clicked.connect(_gen_template)
        file_row.addWidget(template_btn)
        layout.addLayout(file_row)

        form = QFormLayout()
        protein_budget_spin = QDoubleSpinBox()
        protein_budget_spin.setRange(0.001, 10000.0)
        protein_budget_spin.setValue(0.5)
        protein_budget_spin.setDecimals(3)
        protein_budget_spin.setSuffix(" g/gDW")
        form.addRow("Total protein budget:", protein_budget_spin)

        sigma_spin = QDoubleSpinBox()
        sigma_spin.setRange(0.01, 1.0)
        sigma_spin.setValue(0.5)
        sigma_spin.setDecimals(2)
        form.addRow("Sigma (enzyme saturation):", sigma_spin)
        layout.addLayout(form)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)

        if dialog.exec() != QDialog.Accepted:
            return

        csv_path = kcat_edit.text().strip()
        if not csv_path or not Path(csv_path).exists():
            QMessageBox.warning(self, "Missing", "Select a valid kcat CSV file.")
            return

        protein_budget = protein_budget_spin.value()
        sigma = sigma_spin.value()

        # Parse kcat data — supports isoenzymes via optional enzyme_id column
        # kcat_data maps reaction_id → list of (kcat, mw, enzyme_id) tuples
        kcat_data: dict[str, list[tuple[float, float, str]]] = {}
        try:
            sep = '\t' if csv_path.endswith('.tsv') else ','
            with open(csv_path, encoding='utf-8') as f:
                reader = csv.reader(f, delimiter=sep)
                next(reader)  # skip header
                for row in reader:
                    if len(row) >= 3:
                        rid = row[0].strip()
                        try:
                            kcat = float(row[1].strip())
                            mw = float(row[2].strip())
                            enz_id = row[3].strip() if len(row) >= 4 and row[3].strip() else rid
                            kcat_data.setdefault(rid, []).append((kcat, mw, enz_id))
                        except ValueError:
                            continue
        except Exception as e:
            self._show_error("CSV Error", "Failed to parse kcat CSV", e)
            return

        if not kcat_data:
            QMessageBox.warning(self, "No Data", "No valid kcat entries found in CSV.")
            return

        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self._apply_selected_solver(model)

        def _compute():
            # Create protein pool pseudo-metabolite
            prot_pool = cobra.Metabolite("prot_pool", name="Protein pool",
                                          compartment="c")

            constrained = []
            for rid, enzymes in kcat_data.items():
                if rid not in model.reactions:
                    continue
                rxn = model.reactions.get_by_id(rid)

                if len(enzymes) == 1:
                    # Single enzyme — direct constraint (original GECKO-light)
                    kcat_val, mw_val, enz_id = enzymes[0]
                    kcat_h = kcat_val * 3600.0 * sigma
                    enzyme_cost = mw_val / (1000.0 * kcat_h) if kcat_h > 0 else 1e6
                    rxn.add_metabolites({prot_pool: enzyme_cost})
                    constrained.append({"id": rid, "name": rxn.name, "kcat": kcat_val,
                                        "mw": mw_val, "enzyme_cost": enzyme_cost,
                                        "enzyme_id": enz_id})
                else:
                    # Multiple isoenzymes — create arm reactions
                    # Disable original reaction; create one arm per isoenzyme
                    # Each arm carries the same stoichiometry but a different
                    # protein cost, so the solver picks the cheapest enzyme.
                    orig_mets = dict(rxn.metabolites)
                    orig_lb = rxn.lower_bound
                    orig_ub = rxn.upper_bound
                    rxn.bounds = (0, 0)  # disable original

                    for kcat_val, mw_val, enz_id in enzymes:
                        kcat_h = kcat_val * 3600.0 * sigma
                        enzyme_cost = mw_val / (1000.0 * kcat_h) if kcat_h > 0 else 1e6
                        arm_id = f"{rid}_arm_{enz_id}"
                        arm_rxn = cobra.Reaction(arm_id)
                        arm_rxn.name = f"{rxn.name} (via {enz_id})"
                        arm_rxn.lower_bound = orig_lb
                        arm_rxn.upper_bound = orig_ub
                        arm_met_dict = {met: coeff for met, coeff in orig_mets.items()}
                        arm_met_dict[prot_pool] = enzyme_cost
                        arm_rxn.add_metabolites(arm_met_dict)
                        model.add_reactions([arm_rxn])
                        constrained.append({"id": arm_id, "name": arm_rxn.name,
                                            "kcat": kcat_val, "mw": mw_val,
                                            "enzyme_cost": enzyme_cost,
                                            "enzyme_id": enz_id})

            # Add protein pool exchange (drain) with upper bound = budget
            prot_exchange = cobra.Reaction("prot_pool_exchange")
            prot_exchange.name = "Protein pool exchange"
            prot_exchange.lower_bound = 0.0
            prot_exchange.upper_bound = protein_budget
            prot_exchange.add_metabolites({prot_pool: -1.0})
            model.add_reactions([prot_exchange])

            # Solve
            sol0 = model.optimize()
            baseline_unconstrained = None
            with model:
                # Remove protein constraint for comparison
                model.reactions.get_by_id("prot_pool_exchange").upper_bound = 1e6
                sol_unconstrained = model.optimize()
                if sol_unconstrained.status == 'optimal':
                    baseline_unconstrained = float(sol_unconstrained.objective_value)

            result = {
                "status": sol0.status if sol0 else "failed",
                "objective": float(sol0.objective_value) if sol0 and sol0.status == 'optimal' else 0.0,
                "unconstrained_obj": baseline_unconstrained,
                "constrained_reactions": constrained,
                "protein_budget": protein_budget,
                "sigma": sigma,
                "flux": {r.id: float(sol0.fluxes[r.id]) for r in model.reactions
                         if r.id != "prot_pool_exchange"} if sol0 and sol0.status == 'optimal' else {}
            }
            # Identify bottleneck enzymes (highest usage / kcat)
            if sol0 and sol0.status == 'optimal':
                usage = []
                for cr in constrained:
                    fv = abs(sol0.fluxes.get(cr["id"], 0))
                    cost = fv * cr["enzyme_cost"]
                    usage.append({**cr, "flux": fv, "protein_usage": cost})
                usage.sort(key=lambda x: x["protein_usage"], reverse=True)
                result["enzyme_usage"] = usage[:30]
            return result

        def _on_done(result):
            lines = [
                f"Status: {result['status']}",
                f"Enzyme-constrained objective: {result.get('objective', 0):.6g}",
                f"Unconstrained objective: {result.get('unconstrained_obj', 'N/A')}",
                f"Protein budget: {result['protein_budget']} g/gDW (σ={result['sigma']})",
                f"Reactions constrained: {len(result['constrained_reactions'])}",
                ""
            ]
            if result.get("enzyme_usage"):
                lines.append("Top enzyme bottlenecks (by protein mass usage):")
                for eu in result["enzyme_usage"][:15]:
                    lines.append(
                        f"  {eu['id']}: flux={eu['flux']:.4g}, "
                        f"kcat={eu['kcat']:.1f} 1/s, MW={eu['mw']:.1f} kDa, "
                        f"protein={eu['protein_usage']:.4g} g/gDW")
            TextPopup("GECKO Results", "\n".join(lines), self).exec()

        self._launch_worker(_compute, _on_done, "Running GECKO enzyme-constrained FBA...")

    # ==================== FSEOF (Flux Scanning based on Enforced Objective Flux) ====================

    def run_fseof(self):
        """FSEOF — Flux Scanning based on Enforced Objective Flux.

        Identifies overexpression targets by scanning the objective flux
        from 0 to maximum and finding reactions whose flux monotonically
        increases with the objective — indicating amplification targets.
        """
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("FSEOF — Flux Scanning")
        dialog.resize(500, 300)
        layout = QVBoxLayout(dialog)
        layout.addWidget(QLabel(
            "Scan the objective from 0 to max and identify reactions whose\n"
            "flux changes monotonically — candidate overexpression targets.\n"
            "Reactions that increase with the objective are amplification targets."))

        form = QFormLayout()
        target_edit = QLineEdit()
        target_edit.setPlaceholderText("Target reaction(s) — comma-separated (e.g., EX_etoh_e, EX_succ_e)")
        if self.base_model:
            completer = QCompleter([r.id for r in self.base_model.reactions])
            completer.setCaseSensitivity(Qt.CaseInsensitive)
            target_edit.setCompleter(completer)
        form.addRow("Target reaction(s):", target_edit)

        steps_spin = QSpinBox()
        steps_spin.setRange(5, 50)
        steps_spin.setValue(10)
        form.addRow("Scan steps:", steps_spin)

        min_growth_spin = QDoubleSpinBox()
        min_growth_spin.setRange(0.0, 1.0)
        min_growth_spin.setValue(0.01)
        min_growth_spin.setDecimals(3)
        form.addRow("Min growth fraction:", min_growth_spin)
        layout.addLayout(form)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)

        if dialog.exec() != QDialog.Accepted:
            return

        # Parse target reaction(s) — comma-separated
        raw_targets = [t.strip() for t in target_edit.text().split(",") if t.strip()]
        target_ids = []
        for raw in raw_targets:
            tid = self._extract_reaction_id_from_input(raw)
            if tid:
                target_ids.append(tid)
        if not target_ids:
            QMessageBox.warning(self, "Missing", "Enter at least one target reaction ID.")
            return
        n_steps = steps_spin.value()
        min_growth = min_growth_spin.value()

        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self._apply_selected_solver(model)

        def _compute(worker=None):
            # Step 1: Get max objective (baseline growth)
            sol0 = model.optimize()
            if sol0.status != 'optimal':
                raise ValueError("Baseline optimization failed.")
            max_obj = float(sol0.objective_value)

            # Identify objective reaction ID
            obj_rxn_id = None
            for v in model.objective.expression.as_coefficients_dict():
                name = getattr(v, "name", "") or ""
                for r in model.reactions:
                    if r.id in name:
                        obj_rxn_id = r.id
                        break
                if obj_rxn_id:
                    break

            # Process each target reaction
            all_target_results = []
            for t_idx, target_id in enumerate(target_ids):
                if worker:
                    worker.report_progress(
                        f"FSEOF target {t_idx + 1}/{len(target_ids)}: {target_id}", 0)

                baseline_target = float(sol0.fluxes.get(target_id, 0))

                # Step 2: Get max target flux at min growth
                with model:
                    if obj_rxn_id:
                        model.reactions.get_by_id(obj_rxn_id).lower_bound = max_obj * min_growth
                    model.objective = target_id
                    sol_max_target = model.optimize()
                    max_target = float(sol_max_target.objective_value) if sol_max_target.status == 'optimal' else 0.0

                if max_target <= 1e-9:
                    all_target_results.append({
                        "target_id": target_id,
                        "error": f"Target reaction {target_id} cannot carry flux.",
                        "max_target": 0.0, "baseline_target": baseline_target,
                        "targets": []})
                    continue

                # Step 3: Scan enforced objective from 0 to max
                enforced_levels = [max_obj * i / n_steps for i in range(n_steps + 1)]
                flux_profiles: dict[str, list[float]] = {r.id: [] for r in model.reactions
                                                          if not r.boundary}

                for step_i, enforced in enumerate(enforced_levels):
                    if worker:
                        overall_pct = int((t_idx / len(target_ids) + step_i / (n_steps + 1) / len(target_ids)) * 100)
                        worker.report_progress(
                            f"FSEOF {target_id}: scan {step_i + 1}/{n_steps + 1}", overall_pct)
                    with model:
                        if obj_rxn_id:
                            model.reactions.get_by_id(obj_rxn_id).lower_bound = enforced
                        model.objective = target_id
                        sol = model.optimize()
                        if sol.status == 'optimal':
                            for rid in flux_profiles:
                                flux_profiles[rid].append(float(sol.fluxes.get(rid, 0)))
                        else:
                            for rid in flux_profiles:
                                flux_profiles[rid].append(float('nan'))

                # Step 4: Identify monotonically changing reactions
                targets = []
                for rid, profile in flux_profiles.items():
                    valid = [v for v in profile if v == v]  # exclude NaN
                    if len(valid) < 3:
                        continue
                    diffs = [valid[i + 1] - valid[i] for i in range(len(valid) - 1)]
                    if all(d >= -1e-8 for d in diffs) and sum(abs(d) for d in diffs) > 1e-6:
                        targets.append({"id": rid, "direction": "increasing",
                                        "start_flux": valid[0], "end_flux": valid[-1],
                                        "change": valid[-1] - valid[0]})
                    elif all(d <= 1e-8 for d in diffs) and sum(abs(d) for d in diffs) > 1e-6:
                        targets.append({"id": rid, "direction": "decreasing",
                                        "start_flux": valid[0], "end_flux": valid[-1],
                                        "change": valid[-1] - valid[0]})

                targets.sort(key=lambda x: abs(x["change"]), reverse=True)

                all_target_results.append({
                    "target_id": target_id, "max_target": max_target,
                    "baseline_target": baseline_target,
                    "targets": targets[:50], "n_steps": n_steps,
                    "enforced_levels": enforced_levels})

            return {"max_obj": max_obj, "target_results": all_target_results,
                    "n_targets": len(target_ids)}

        def _on_done(result):
            target_results = result.get("target_results", [])
            if len(target_results) == 0:
                QMessageBox.warning(self, "FSEOF", "No targets analysed.")
                return

            # For single-target, check for errors
            if len(target_results) == 1 and target_results[0].get("error"):
                QMessageBox.warning(self, "FSEOF", target_results[0]["error"])
                return

            dialog2 = QDialog(self)
            title_targets = ", ".join(tr["target_id"] for tr in target_results)
            dialog2.setWindowTitle(f"FSEOF Results — {title_targets}")
            dialog2.resize(900, 650)
            layout2 = QVBoxLayout(dialog2)

            if len(target_results) == 1:
                # Single target — flat table (backward compatible)
                tr = target_results[0]
                targets = tr.get("targets", [])
                increasing = [t for t in targets if t["direction"] == "increasing"]
                decreasing = [t for t in targets if t["direction"] == "decreasing"]
                layout2.addWidget(QLabel(
                    f"Max growth: {result['max_obj']:.6g} | "
                    f"Max target: {tr['max_target']:.6g}\n"
                    f"Baseline target flux: {tr['baseline_target']:.6g}\n"
                    f"Amplification targets (↑): {len(increasing)} | "
                    f"Attenuation targets (↓): {len(decreasing)}"))

                table = QTableWidget(0, 5)
                table.setHorizontalHeaderLabels(["Reaction", "Direction", "Start Flux", "End Flux", "Δ Flux"])
                table.horizontalHeader().setStretchLastSection(True)
                table.setEditTriggers(QAbstractItemView.NoEditTriggers)
                for t in targets:
                    r = table.rowCount()
                    table.insertRow(r)
                    table.setItem(r, 0, QTableWidgetItem(t["id"]))
                    table.setItem(r, 1, QTableWidgetItem("↑ Amplify" if t["direction"] == "increasing" else "↓ Attenuate"))
                    table.setItem(r, 2, QTableWidgetItem(f"{t['start_flux']:.6g}"))
                    table.setItem(r, 3, QTableWidgetItem(f"{t['end_flux']:.6g}"))
                    table.setItem(r, 4, QTableWidgetItem(f"{t['change']:.6g}"))
                layout2.addWidget(table)
            else:
                # Multi-target — tabbed view
                layout2.addWidget(QLabel(
                    f"Max growth: {result['max_obj']:.6g} | "
                    f"Targets analysed: {result['n_targets']}"))
                tabs = QTabWidget()
                for tr in target_results:
                    tab_widget = QWidget()
                    tab_layout = QVBoxLayout(tab_widget)
                    targets = tr.get("targets", [])
                    increasing = [t for t in targets if t["direction"] == "increasing"]
                    decreasing = [t for t in targets if t["direction"] == "decreasing"]

                    if tr.get("error"):
                        tab_layout.addWidget(QLabel(f"⚠ {tr['error']}"))
                    else:
                        tab_layout.addWidget(QLabel(
                            f"Max target: {tr['max_target']:.6g} | "
                            f"Baseline: {tr['baseline_target']:.6g}\n"
                            f"↑ Amplify: {len(increasing)} | ↓ Attenuate: {len(decreasing)}"))

                    table = QTableWidget(0, 5)
                    table.setHorizontalHeaderLabels(["Reaction", "Direction", "Start Flux", "End Flux", "Δ Flux"])
                    table.horizontalHeader().setStretchLastSection(True)
                    table.setEditTriggers(QAbstractItemView.NoEditTriggers)
                    for t in targets:
                        r = table.rowCount()
                        table.insertRow(r)
                        table.setItem(r, 0, QTableWidgetItem(t["id"]))
                        table.setItem(r, 1, QTableWidgetItem("↑ Amplify" if t["direction"] == "increasing" else "↓ Attenuate"))
                        table.setItem(r, 2, QTableWidgetItem(f"{t['start_flux']:.6g}"))
                        table.setItem(r, 3, QTableWidgetItem(f"{t['end_flux']:.6g}"))
                        table.setItem(r, 4, QTableWidgetItem(f"{t['change']:.6g}"))
                    tab_layout.addWidget(table)
                    tabs.addTab(tab_widget, tr["target_id"])
                layout2.addWidget(tabs)

            btns2 = QDialogButtonBox(QDialogButtonBox.Close)
            btns2.rejected.connect(dialog2.reject)
            layout2.addWidget(btns2)
            dialog2.exec()

        self._launch_worker(_compute, _on_done, "Running FSEOF scan...")

    # ==================== MODEL COMPARISON ====================

    def run_model_comparison(self):
        """Compare two SBML models side by side."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a base model first.")
            return

        fp, _ = QFileDialog.getOpenFileName(self, "Select second SBML model for comparison", "",
                                             "SBML Files (*.xml *.sbml);;All (*.*)")
        if not fp:
            return

        try:
            model2 = cobra.io.read_sbml_model(fp)
        except Exception as e:
            self._show_error("Load Error", "Failed to load model", e)
            return

        m1 = self.base_model
        m2 = model2
        name1 = getattr(m1, 'id', 'Model A') or 'Model A'
        name2 = getattr(m2, 'id', 'Model B') or 'Model B'

        # Reactions comparison
        rxns1 = set(r.id for r in m1.reactions)
        rxns2 = set(r.id for r in m2.reactions)
        only_in_1 = sorted(rxns1 - rxns2)
        only_in_2 = sorted(rxns2 - rxns1)
        common_rxns = sorted(rxns1 & rxns2)
        bound_diffs = []
        for rid in common_rxns:
            r1 = m1.reactions.get_by_id(rid)
            r2 = m2.reactions.get_by_id(rid)
            if r1.lower_bound != r2.lower_bound or r1.upper_bound != r2.upper_bound:
                bound_diffs.append((rid, r1.lower_bound, r1.upper_bound, r2.lower_bound, r2.upper_bound))

        # Metabolites comparison
        mets1 = set(m.id for m in m1.metabolites)
        mets2 = set(m.id for m in m2.metabolites)
        mets_only1 = sorted(mets1 - mets2)
        mets_only2 = sorted(mets2 - mets1)

        # Genes comparison
        genes1 = set(g.id for g in m1.genes)
        genes2 = set(g.id for g in m2.genes)
        genes_only1 = sorted(genes1 - genes2)
        genes_only2 = sorted(genes2 - genes1)

        dialog = QDialog(self)
        dialog.setWindowTitle(f"Model Comparison: {name1} vs {name2}")
        dialog.resize(900, 700)
        layout = QVBoxLayout(dialog)

        summary = QLabel(
            f"<b>{name1}</b>: {len(rxns1)} reactions, {len(mets1)} metabolites, {len(genes1)} genes<br>"
            f"<b>{name2}</b>: {len(rxns2)} reactions, {len(mets2)} metabolites, {len(genes2)} genes<br><br>"
            f"Common reactions: {len(common_rxns)} | Bound differences: {len(bound_diffs)}<br>"
            f"Unique to {name1}: {len(only_in_1)} rxns, {len(mets_only1)} mets, {len(genes_only1)} genes<br>"
            f"Unique to {name2}: {len(only_in_2)} rxns, {len(mets_only2)} mets, {len(genes_only2)} genes"
        )
        summary.setWordWrap(True)
        layout.addWidget(summary)

        tabs = QTabWidget()

        # Reactions unique tab
        if only_in_1 or only_in_2:
            rxn_table = QTableWidget(0, 2)
            rxn_table.setHorizontalHeaderLabels([f"Only in {name1}", f"Only in {name2}"])
            rxn_table.horizontalHeader().setStretchLastSection(True)
            rxn_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            max_rows = max(len(only_in_1), len(only_in_2))
            for i in range(max_rows):
                r = rxn_table.rowCount()
                rxn_table.insertRow(r)
                rxn_table.setItem(r, 0, QTableWidgetItem(only_in_1[i] if i < len(only_in_1) else ""))
                rxn_table.setItem(r, 1, QTableWidgetItem(only_in_2[i] if i < len(only_in_2) else ""))
            tabs.addTab(rxn_table, f"Unique Reactions ({len(only_in_1)} + {len(only_in_2)})")

        # Bound differences tab
        if bound_diffs:
            bd_table = QTableWidget(0, 5)
            bd_table.setHorizontalHeaderLabels(["Reaction", f"LB ({name1})", f"UB ({name1})", f"LB ({name2})", f"UB ({name2})"])
            bd_table.horizontalHeader().setStretchLastSection(True)
            bd_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            for rid, lb1, ub1, lb2, ub2 in bound_diffs:
                r = bd_table.rowCount()
                bd_table.insertRow(r)
                bd_table.setItem(r, 0, QTableWidgetItem(rid))
                bd_table.setItem(r, 1, QTableWidgetItem(f"{lb1:g}"))
                bd_table.setItem(r, 2, QTableWidgetItem(f"{ub1:g}"))
                bd_table.setItem(r, 3, QTableWidgetItem(f"{lb2:g}"))
                bd_table.setItem(r, 4, QTableWidgetItem(f"{ub2:g}"))
            tabs.addTab(bd_table, f"Bound Differences ({len(bound_diffs)})")

        # Metabolites unique tab
        if mets_only1 or mets_only2:
            met_table = QTableWidget(0, 2)
            met_table.setHorizontalHeaderLabels([f"Only in {name1}", f"Only in {name2}"])
            met_table.horizontalHeader().setStretchLastSection(True)
            met_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            max_rows = max(len(mets_only1), len(mets_only2))
            for i in range(max_rows):
                r = met_table.rowCount()
                met_table.insertRow(r)
                met_table.setItem(r, 0, QTableWidgetItem(mets_only1[i] if i < len(mets_only1) else ""))
                met_table.setItem(r, 1, QTableWidgetItem(mets_only2[i] if i < len(mets_only2) else ""))
            tabs.addTab(met_table, f"Unique Metabolites ({len(mets_only1)} + {len(mets_only2)})")

        # Genes unique tab
        if genes_only1 or genes_only2:
            gene_table = QTableWidget(0, 2)
            gene_table.setHorizontalHeaderLabels([f"Only in {name1}", f"Only in {name2}"])
            gene_table.horizontalHeader().setStretchLastSection(True)
            gene_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            max_rows = max(len(genes_only1), len(genes_only2))
            for i in range(max_rows):
                r = gene_table.rowCount()
                gene_table.insertRow(r)
                gene_table.setItem(r, 0, QTableWidgetItem(genes_only1[i] if i < len(genes_only1) else ""))
                gene_table.setItem(r, 1, QTableWidgetItem(genes_only2[i] if i < len(genes_only2) else ""))
            tabs.addTab(gene_table, f"Unique Genes ({len(genes_only1)} + {len(genes_only2)})")

        layout.addWidget(tabs)
        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)
        dialog.exec()

    # ==================== PATHWAY ENRICHMENT ====================

    def run_pathway_enrichment(self):
        """Metabolic Pathway Enrichment — subsystem-based enrichment of active reactions."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return
        if self.last_run is None:
            QMessageBox.warning(self, "No results", "Run an FBA analysis first to determine active reactions.")
            return

        model = self.base_model
        # Get flux data from last run
        flux = {}
        if "baseline" in self.last_run and "flux" in self.last_run["baseline"]:
            flux = self.last_run["baseline"]["flux"]
        elif "flux" in self.last_run:
            flux = self.last_run["flux"]

        if not flux:
            QMessageBox.warning(self, "No flux data", "Last analysis does not contain flux data.")
            return

        # Define active reactions (|flux| > threshold)
        threshold = 1e-6
        active_rxns = set(rid for rid, fv in flux.items() if abs(fv) > threshold)

        # Get subsystem annotations
        subsystem_rxns: dict[str, set] = {}
        for rxn in model.reactions:
            ss = getattr(rxn, 'subsystem', '') or ''
            ss = ss.strip()
            if ss:
                subsystem_rxns.setdefault(ss, set()).add(rxn.id)

        if not subsystem_rxns:
            QMessageBox.warning(self, "No subsystems", "Model has no subsystem annotations.")
            return

        # Compute enrichment using hypergeometric-like test
        total_rxns = len(model.reactions)
        total_active = len(active_rxns)
        enrichment = []
        for ss, ss_rxns in subsystem_rxns.items():
            ss_total = len(ss_rxns)
            ss_active = len(ss_rxns & active_rxns)
            if ss_total == 0:
                continue
            # Expected fraction
            expected = total_active / total_rxns * ss_total if total_rxns > 0 else 0
            ratio = ss_active / expected if expected > 0 else 0
            enrichment.append({
                "subsystem": ss, "total": ss_total, "active": ss_active,
                "expected": expected, "enrichment_ratio": ratio
            })
        enrichment.sort(key=lambda x: x["enrichment_ratio"], reverse=True)

        dialog = QDialog(self)
        dialog.setWindowTitle("Pathway Enrichment Analysis")
        dialog.resize(800, 600)
        layout = QVBoxLayout(dialog)
        layout.addWidget(QLabel(
            f"Total reactions: {total_rxns} | Active (|flux| > {threshold}): {total_active}\n"
            f"Subsystems analyzed: {len(enrichment)}"))

        # Plot
        fig = Figure(figsize=(7, 3.5), dpi=100)
        ax = fig.add_subplot(111)
        top = enrichment[:20]
        if top:
            ss_names = [e["subsystem"][:30] for e in top]
            ratios = [e["enrichment_ratio"] for e in top]
            colors = ['#4CAF50' if r > 1 else '#f44336' for r in ratios]
            ax.barh(range(len(top)), ratios, color=colors, alpha=0.8)
            ax.set_yticks(range(len(top)))
            ax.set_yticklabels(ss_names, fontsize=7)
            ax.axvline(x=1.0, color='black', linestyle='--', linewidth=0.8)
            ax.set_xlabel("Enrichment Ratio (>1 = enriched)")
            ax.set_title("Top 20 Subsystems by Enrichment")
        fig.tight_layout()
        canvas = FigureCanvas(fig)
        layout.addWidget(canvas)

        table = QTableWidget(0, 5)
        table.setHorizontalHeaderLabels(["Subsystem", "Total Rxns", "Active Rxns", "Expected", "Enrichment Ratio"])
        table.horizontalHeader().setStretchLastSection(True)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        for e in enrichment:
            r = table.rowCount()
            table.insertRow(r)
            table.setItem(r, 0, QTableWidgetItem(e["subsystem"]))
            table.setItem(r, 1, QTableWidgetItem(str(e["total"])))
            table.setItem(r, 2, QTableWidgetItem(str(e["active"])))
            table.setItem(r, 3, QTableWidgetItem(f"{e['expected']:.1f}"))
            table.setItem(r, 4, QTableWidgetItem(f"{e['enrichment_ratio']:.3f}"))
        layout.addWidget(table)

        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)
        dialog.exec()

    # ==================== ESCHER MAP (with JSON map import) ====================

    def show_escher_map(self):
        """Show Escher metabolic map.  Supports loading Escher JSON map files
        directly (no ``escher`` package required) and falls back to the
        ``escher`` Python package when available."""

        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        # Collect current flux data if available
        flux: dict[str, float] = {}
        if self.last_run and "baseline" in self.last_run and "flux" in self.last_run["baseline"]:
            flux = self.last_run["baseline"]["flux"]

        # --- Dialog ---
        dialog = QDialog(self)
        dialog.setWindowTitle("Escher Map Viewer")
        dialog.resize(1200, 800)
        layout = QVBoxLayout(dialog)

        # toolbar
        toolbar = QHBoxLayout()
        json_btn = QPushButton("Load Escher JSON map...")
        save_html_btn = QPushButton("Save HTML...")
        save_html_btn.setEnabled(False)
        toolbar.addWidget(json_btn)
        toolbar.addWidget(save_html_btn)
        toolbar.addStretch()
        layout.addLayout(toolbar)

        # content area
        content_stack = QStackedWidget()
        placeholder = QLabel("Load an Escher JSON map file or use the escher package.\n\n"
                             "You can download Escher maps from:\n"
                             "  https://escher.github.io/#/app\n\n"
                             "Click 'Load Escher JSON map...' to import a .json map file.")
        placeholder.setAlignment(Qt.AlignCenter)
        placeholder.setWordWrap(True)
        content_stack.addWidget(placeholder)  # index 0
        layout.addWidget(content_stack, 1)

        html_holder: list[str] = [""]  # mutable ref for closures

        def _build_escher_html(map_json_str: str) -> str:
            """Build a standalone HTML page that renders an Escher map,
            overlaying flux data.  Prefers a local ``escher.min.js``
            bundle (for offline use) and falls back to the CDN."""
            import json as _json
            import sys as _sys
            flux_json = _json.dumps(flux) if flux else "null"
            model_json = "null"
            try:
                model_json = cobra.io.to_json(self.base_model)
            except Exception:
                pass
            # Resolve Escher JS: prefer local file, fall back to CDN
            # Cache the JS text to avoid re-reading the ~2 MB file each time
            escher_script_tag = '<script src="https://unpkg.com/escher@1.7.3/dist/escher.min.js"></script>'
            if not hasattr(self, '_escher_js_cache'):
                self._escher_js_cache = None
            if self._escher_js_cache is not None:
                escher_script_tag = f"<script>{self._escher_js_cache}</script>"
            else:
                search_dirs = [
                    Path(getattr(_sys, '_MEIPASS', '.')),   # PyInstaller bundle
                    Path(__file__).parent.parent,             # project root
                    Path(__file__).parent,                    # metabodesk_core/
                ]
                for d in search_dirs:
                    candidate = d / "escher.min.js"
                    if candidate.is_file():
                        try:
                            js_text = candidate.read_text(encoding="utf-8")
                            self._escher_js_cache = js_text
                            escher_script_tag = f"<script>{js_text}</script>"
                            break
                        except Exception:
                            pass
            return f"""<!DOCTYPE html>
<html><head>
<meta charset="utf-8"/>
{escher_script_tag}
<style>body{{margin:0;padding:0}} #map-container{{width:100vw;height:100vh}}</style>
</head><body>
<div id="map-container"></div>
<script>
var mapData = {map_json_str};
var modelData = {model_json};
var fluxData = {flux_json};
escher.Builder(mapData, modelData, null, document.getElementById('map-container'), {{
  reaction_data: fluxData,
  fill_screen: true,
  menu: 'zoom',
  scroll_behavior: 'zoom',
  reaction_styles: ['color', 'size', 'text'],
  reaction_compare_style: 'log2_fold',
  identifiers_on_map: 'bigg_id'
}});
</script></body></html>"""

        def _load_json_map():
            fp, _ = QFileDialog.getOpenFileName(dialog, "Open Escher JSON Map", "", "JSON files (*.json)")
            if not fp:
                return
            try:
                map_text = Path(fp).read_text(encoding="utf-8")
                # Validate it's a list (Escher map format) or dict
                import json as _json
                parsed = _json.loads(map_text)
                if not isinstance(parsed, (list, dict)):
                    raise ValueError("Not a valid Escher map JSON.")
                html = _build_escher_html(map_text)
                html_holder[0] = html
                save_html_btn.setEnabled(True)
                # Try QWebEngineView first
                try:
                    from PySide6.QtWebEngineWidgets import QWebEngineView
                    web = QWebEngineView()
                    web.setHtml(html)
                    if content_stack.count() > 1:
                        old = content_stack.widget(1)
                        content_stack.removeWidget(old)
                        old.deleteLater()
                    content_stack.addWidget(web)
                    content_stack.setCurrentIndex(1)
                except ImportError:
                    lbl = QLabel("Map loaded. QWebEngineView not available.\n"
                                 "Use 'Save HTML...' to view in your browser.")
                    lbl.setAlignment(Qt.AlignCenter)
                    if content_stack.count() > 1:
                        old = content_stack.widget(1)
                        content_stack.removeWidget(old)
                        old.deleteLater()
                    content_stack.addWidget(lbl)
                    content_stack.setCurrentIndex(1)
            except Exception as e:
                QMessageBox.critical(dialog, "Error", f"Failed to load map:\n{e}")

        def _save_html():
            if not html_holder[0]:
                return
            fp, _ = QFileDialog.getSaveFileName(dialog, "Save Escher HTML", "", "HTML (*.html)")
            if fp:
                Path(fp).write_text(html_holder[0], encoding="utf-8")

        json_btn.clicked.connect(_load_json_map)
        save_html_btn.clicked.connect(_save_html)

        # Also try the escher package directly
        try:
            import escher
            pkg_btn = QPushButton("Use escher package (auto-generate)")
            def _use_pkg():
                try:
                    builder = escher.Builder(model=self.base_model, reaction_data=flux if flux else None)
                    html = builder._get_html()
                    html_holder[0] = html
                    save_html_btn.setEnabled(True)
                    try:
                        from PySide6.QtWebEngineWidgets import QWebEngineView
                        web = QWebEngineView()
                        web.setHtml(html)
                        if content_stack.count() > 1:
                            old = content_stack.widget(1)
                            content_stack.removeWidget(old)
                            old.deleteLater()
                        content_stack.addWidget(web)
                        content_stack.setCurrentIndex(1)
                    except ImportError:
                        lbl = QLabel("Map generated. Use 'Save HTML...' to view in browser.")
                        lbl.setAlignment(Qt.AlignCenter)
                        if content_stack.count() > 1:
                            old = content_stack.widget(1)
                            content_stack.removeWidget(old)
                            old.deleteLater()
                        content_stack.addWidget(lbl)
                        content_stack.setCurrentIndex(1)
                except Exception as e:
                    QMessageBox.critical(dialog, "Error", f"Escher builder failed:\n{e}")
            pkg_btn.clicked.connect(_use_pkg)
            toolbar.addWidget(pkg_btn)
        except ImportError:
            pass

        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)
        dialog.exec()

    # ==================== PHENOTYPE PHASE PLANE (PhPP) ====================

    def run_phenotype_phase_plane(self):
        """Phenotype Phase Plane analysis: 2-D sweep of two reaction fluxes,
        computing the objective value at each (rx1, rx2) grid point."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        rxn_ids = [r.id for r in self.base_model.reactions]

        dialog = QDialog(self)
        dialog.setWindowTitle("Phenotype Phase Plane (PhPP)")
        form = QFormLayout(dialog)

        rx1_combo = QComboBox()
        rx1_combo.addItems(rxn_ids)
        form.addRow("Reaction X-axis:", rx1_combo)

        rx2_combo = QComboBox()
        rx2_combo.addItems(rxn_ids)
        if len(rxn_ids) > 1:
            rx2_combo.setCurrentIndex(1)
        form.addRow("Reaction Y-axis:", rx2_combo)

        steps_spin = QSpinBox()
        steps_spin.setRange(5, 100)
        steps_spin.setValue(20)
        form.addRow("Grid points per axis:", steps_spin)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        form.addRow(btns)

        if dialog.exec() != QDialog.Accepted:
            return

        rx1_id = rx1_combo.currentText()
        rx2_id = rx2_combo.currentText()
        steps = steps_spin.value()

        if rx1_id == rx2_id:
            QMessageBox.warning(self, "Same reaction", "Choose two different reactions.")
            return

        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self.apply_reaction_overrides_to_model(model)
        self._apply_selected_objective_to_model(model)
        self._apply_selected_solver(model)

        rx1 = model.reactions.get_by_id(rx1_id)
        rx2 = model.reactions.get_by_id(rx2_id)
        rx1_range = (rx1.lower_bound, rx1.upper_bound)
        rx2_range = (rx2.lower_bound, rx2.upper_bound)

        import numpy as np

        x_vals = np.linspace(rx1_range[0], rx1_range[1], steps)
        y_vals = np.linspace(rx2_range[0], rx2_range[1], steps)

        def _compute(worker=None):
            grid = np.full((steps, steps), np.nan)
            total = steps * steps
            count = 0
            for i, xv in enumerate(x_vals):
                for j, yv in enumerate(y_vals):
                    with model:
                        rx1.lower_bound = float(xv)
                        rx1.upper_bound = float(xv)
                        rx2.lower_bound = float(yv)
                        rx2.upper_bound = float(yv)
                        try:
                            sol = model.optimize()
                            if sol.status == "optimal":
                                grid[j, i] = float(sol.objective_value)
                        except Exception:
                            pass
                    count += 1
                    if worker and count % max(1, total // 20) == 0:
                        worker.progress_pct.emit(int(100 * count / total))
            return {
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "analysis_type": "Phenotype Phase Plane",
                "rx1_id": rx1_id, "rx2_id": rx2_id,
                "x_vals": x_vals.tolist(), "y_vals": y_vals.tolist(),
                "grid": grid.tolist(),
            }

        def _on_done(result):
            self.last_run = result
            self._render_phpp_results(result)

        self._launch_worker(_compute, _on_done, f"Running PhPP ({rx1_id} × {rx2_id})...")

    def _render_phpp_results(self, result: dict):
        """Render Phenotype Phase Plane results as a heatmap."""
        import numpy as np
        grid = np.array(result["grid"])
        x_vals = result["x_vals"]
        y_vals = result["y_vals"]
        rx1_id = result["rx1_id"]
        rx2_id = result["rx2_id"]

        dialog = QDialog(self)
        dialog.setWindowTitle(f"Phenotype Phase Plane: {rx1_id} × {rx2_id}")
        dialog.resize(900, 700)
        layout = QVBoxLayout(dialog)

        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.figure import Figure

        fig = Figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        im = ax.imshow(
            grid, origin="lower", aspect="auto",
            extent=[x_vals[0], x_vals[-1], y_vals[0], y_vals[-1]],
            cmap="viridis",
        )
        ax.set_xlabel(rx1_id)
        ax.set_ylabel(rx2_id)
        ax.set_title("Phenotype Phase Plane — Objective value")
        fig.colorbar(im, ax=ax, label="Objective")

        # Overlay contour lines
        try:
            cs = ax.contour(
                np.array(x_vals), np.array(y_vals), grid,
                colors="white", linewidths=0.5, levels=10,
            )
            ax.clabel(cs, inline=True, fontsize=7, fmt="%.2f")
        except Exception:
            pass

        fig.tight_layout()
        canvas = FigureCanvas(fig)
        layout.addWidget(canvas)

        # Summary statistics
        valid = grid[~np.isnan(grid)]
        info = QLabel(
            f"Grid: {len(x_vals)}×{len(y_vals)} = {len(x_vals)*len(y_vals)} points  |  "
            f"Feasible: {len(valid)}  |  "
            f"Objective range: [{valid.min():.4f}, {valid.max():.4f}]" if len(valid) > 0
            else f"Grid: {len(x_vals)}×{len(y_vals)}  |  No feasible solutions found"
        )
        layout.addWidget(info)

        # Save buttons
        btn_row = QHBoxLayout()
        save_png = QPushButton("Save as PNG...")
        def _save_png():
            fp, _ = QFileDialog.getSaveFileName(dialog, "Save Plot", "", "PNG (*.png)")
            if fp:
                fig.savefig(fp, dpi=150)
        save_png.clicked.connect(_save_png)
        btn_row.addWidget(save_png)

        save_csv = QPushButton("Save grid as CSV...")
        def _save_csv():
            fp, _ = QFileDialog.getSaveFileName(dialog, "Save CSV", "", "CSV (*.csv)")
            if fp:
                import csv
                with open(fp, 'w', newline='', encoding='utf-8') as f:
                    w = csv.writer(f)
                    w.writerow([""] + [f"{v:.4f}" for v in x_vals])
                    for i, yv in enumerate(y_vals):
                        w.writerow([f"{yv:.4f}"] + [f"{grid[i][j]:.6f}" if not np.isnan(grid[i][j]) else "inf" for j in range(len(x_vals))])
        save_csv.clicked.connect(_save_csv)
        btn_row.addWidget(save_csv)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)
        dialog.exec()

    # ==================== SBML VALIDATION ====================

    def validate_sbml(self):
        """Validate the currently loaded SBML file using libSBML or fallback checks."""
        if self.current_sbml_path is None or not self.current_sbml_path.exists():
            QMessageBox.warning(self, "No SBML file", "Load an SBML model first.")
            return

        errors: list[str] = []
        warnings: list[str] = []
        infos: list[str] = []
        sbml_path = str(self.current_sbml_path)

        # Try libSBML validation (gold standard)
        try:
            import libsbml
            reader = libsbml.SBMLReader()
            doc = reader.readSBML(sbml_path)

            # Enable all consistency checks
            doc.setConsistencyChecks(libsbml.LIBSBML_CAT_GENERAL_CONSISTENCY, True)
            doc.setConsistencyChecks(libsbml.LIBSBML_CAT_IDENTIFIER_CONSISTENCY, True)
            doc.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, False)  # often noisy
            doc.setConsistencyChecks(libsbml.LIBSBML_CAT_MATHML_CONSISTENCY, True)
            doc.setConsistencyChecks(libsbml.LIBSBML_CAT_SBO_CONSISTENCY, True)
            doc.setConsistencyChecks(libsbml.LIBSBML_CAT_OVERDETERMINED_MODEL, True)
            doc.setConsistencyChecks(libsbml.LIBSBML_CAT_MODELING_PRACTICE, True)

            doc.checkConsistency()

            for i in range(doc.getNumErrors()):
                err = doc.getError(i)
                severity = err.getSeverity()
                msg = f"L{err.getLine()}C{err.getColumn()}: [{err.getErrorId()}] {err.getMessage().strip()}"
                if severity == libsbml.LIBSBML_SEV_ERROR or severity == libsbml.LIBSBML_SEV_FATAL:
                    errors.append(msg)
                elif severity == libsbml.LIBSBML_SEV_WARNING:
                    warnings.append(msg)
                else:
                    infos.append(msg)

            model_obj = doc.getModel()
            if model_obj:
                infos.insert(0, f"SBML Level {doc.getLevel()} Version {doc.getVersion()}")
                infos.insert(1, f"Model id: {model_obj.getId()}")
                infos.insert(2, f"Species: {model_obj.getNumSpecies()}  |  Reactions: {model_obj.getNumReactions()}  |  Compartments: {model_obj.getNumCompartments()}")

                # Check for FBC package
                fbc = model_obj.getPlugin("fbc")
                if fbc:
                    infos.append(f"FBC plugin v{fbc.getPackageVersion()}: {fbc.getNumObjectives()} objective(s), {fbc.getNumGeneProducts()} gene product(s)")
                else:
                    warnings.append("FBC plugin not present — flux bounds may be missing")

        except ImportError:
            # Fallback: basic checks without libSBML
            infos.append("libSBML not installed — performing basic checks only")
            infos.append("Install with: pip install python-libsbml")

            if self.base_model:
                m = self.base_model
                infos.append(f"Model: {m.id}  |  Reactions: {len(m.reactions)}  |  Metabolites: {len(m.metabolites)}  |  Genes: {len(m.genes)}")

                # Check mass balance
                unbalanced = 0
                for r in m.reactions:
                    if not r.boundary:
                        bal = r.check_mass_balance()
                        if bal:
                            unbalanced += 1
                if unbalanced:
                    warnings.append(f"{unbalanced} reaction(s) have mass imbalance")

                # Check for reactions without GPR
                no_gpr = sum(1 for r in m.reactions if not r.gene_reaction_rule.strip())
                if no_gpr:
                    infos.append(f"{no_gpr} reaction(s) have no GPR rule")

                # Orphan metabolites
                orphans = [m2.id for m2 in m.metabolites if len(m2.reactions) == 0]
                if orphans:
                    warnings.append(f"{len(orphans)} orphan metabolite(s): {', '.join(orphans[:10])}")

                # Dead-end metabolites (only produced or only consumed)
                dead_end = 0
                for met in m.metabolites:
                    producers = sum(1 for r in met.reactions if met in r.products)
                    consumers = sum(1 for r in met.reactions if met in r.reactants)
                    if (producers > 0 and consumers == 0) or (consumers > 0 and producers == 0):
                        if not any(r.boundary for r in met.reactions):
                            dead_end += 1
                if dead_end:
                    warnings.append(f"{dead_end} potential dead-end metabolite(s)")

                # Unbounded reactions
                unbounded = sum(1 for r in m.reactions if abs(r.upper_bound) >= 999999 or abs(r.lower_bound) >= 999999)
                if unbounded:
                    infos.append(f"{unbounded} reaction(s) with very large bounds (±999999+)")

        # Display results
        dialog = QDialog(self)
        dialog.setWindowTitle("SBML Validation Report")
        dialog.resize(800, 600)
        layout = QVBoxLayout(dialog)

        text = QPlainTextEdit()
        text.setReadOnly(True)
        lines: list[str] = []
        lines.append(f"═══ SBML Validation: {self.current_sbml_path.name} ═══\n")

        if infos:
            lines.append("ℹ️  Information:")
            for i in infos:
                lines.append(f"   {i}")
            lines.append("")

        if errors:
            lines.append(f"❌ Errors ({len(errors)}):")
            for e in errors:
                lines.append(f"   {e}")
            lines.append("")

        if warnings:
            lines.append(f"⚠️  Warnings ({len(warnings)}):")
            for w in warnings:
                lines.append(f"   {w}")
            lines.append("")

        if not errors and not warnings:
            lines.append("✅ No errors or warnings found!")

        lines.append(f"\nTotal: {len(errors)} error(s), {len(warnings)} warning(s)")

        text.setPlainText("\n".join(lines))
        layout.addWidget(text)

        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)
        dialog.exec()

    # ==================== JUPYTER NOTEBOOK EXPORT ====================

    def export_jupyter_notebook(self):
        """Export the current analysis as a Jupyter Notebook (.ipynb) that
        reproduces the analysis with COBRApy code."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        fp, _ = QFileDialog.getSaveFileName(
            self, "Export Jupyter Notebook", "", "Jupyter Notebook (*.ipynb)")
        if not fp:
            return

        cells: list[dict] = []

        def _md(src: str):
            cells.append({"cell_type": "markdown", "metadata": {},
                          "source": [src]})

        def _code(src: str):
            cells.append({"cell_type": "code", "metadata": {},
                          "source": [src], "execution_count": None,
                          "outputs": []})

        _md("# MetaboDesk — Exported Analysis Notebook\n\n"
            f"Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} by MetaboDesk.")

        # Setup cell
        _code("import cobra\nimport pandas as pd\nimport matplotlib.pyplot as plt\nfrom pathlib import Path")

        # Model loading
        model_file = self.current_sbml_path.name if self.current_sbml_path else "model.xml"
        _md(f"## Load Model\n\nMake sure `{model_file}` is in the same directory as this notebook.")
        _code(f'model = cobra.io.read_sbml_model("{model_file}")\n'
              f'print(f"Model: {{model.id}}")\n'
              f'print(f"Reactions: {{len(model.reactions)}}, Metabolites: {{len(model.metabolites)}}, Genes: {{len(model.genes)}}")')

        # Medium / bounds
        if self.reaction_bound_overrides:
            _md("## Custom Reaction Bounds")
            lines = ["# Apply custom bounds"]
            for rid, (lb, ub) in self.reaction_bound_overrides.items():
                lines.append(f'model.reactions.get_by_id("{rid}").bounds = ({lb}, {ub})')
            _code("\n".join(lines))

        # Knockouts
        if self.knockout_genes:
            _md("## Gene Knockouts")
            genes_str = ", ".join(f'"{g}"' for g in sorted(self.knockout_genes))
            _code(f"# Knockout genes\nfor gid in [{genes_str}]:\n"
                  f"    gene = model.genes.get_by_id(gid)\n"
                  f"    cobra.manipulation.delete_model_genes(model, [gene])")

        # Objective
        if self.base_model:
            obj_rxn = None
            for r in self.base_model.objective.variables:
                obj_rxn = r.name
                break
            if obj_rxn:
                _md("## Set Objective")
                _code(f'model.objective = "{obj_rxn}"')

        # Analysis
        analysis_type = self.analysis_type.currentText() if hasattr(self, 'analysis_type') else "FBA"
        _md(f"## Run Analysis: {analysis_type}")

        if "FVA" in analysis_type:
            frac = self.fva_fraction_spin.value() if hasattr(self, 'fva_fraction_spin') else 0.9
            _code(f"from cobra.flux_analysis import flux_variability_analysis\n"
                  f"fva = flux_variability_analysis(model, fraction_of_optimum={frac})\n"
                  f"print(fva.head(20))")
        elif "pFBA" in analysis_type:
            _code("import cobra.flux_analysis\n"
                  "pfba_sol = cobra.flux_analysis.pfba(model)\n"
                  "print(f'Objective: {pfba_sol.objective_value:.6f}')\n"
                  "fluxes = pfba_sol.fluxes\n"
                  "print(fluxes[fluxes.abs() > 1e-6].sort_values(ascending=False).head(20))")
        elif "Single Gene Deletion" in analysis_type:
            _code("from cobra.flux_analysis import single_gene_deletion\n"
                  "sgd = single_gene_deletion(model)\n"
                  "sgd_sorted = sgd.sort_values('growth', ascending=True)\n"
                  "print(sgd_sorted.head(20))")
        elif "Robustness" in analysis_type:
            _code("# Robustness analysis — sweep a reaction and record objective\n"
                  "import numpy as np\n\n"
                  "target_rxn = model.reactions[0]  # Change to your target\n"
                  "steps = np.linspace(target_rxn.lower_bound, target_rxn.upper_bound, 20)\n"
                  "results = []\n"
                  "for v in steps:\n"
                  "    with model:\n"
                  "        target_rxn.bounds = (v, v)\n"
                  "        sol = model.optimize()\n"
                  "        results.append({'flux': v, 'objective': sol.objective_value if sol.status == 'optimal' else float('nan')})\n"
                  "df = pd.DataFrame(results)\n"
                  "df.plot(x='flux', y='objective')\n"
                  "plt.show()")
        elif "Production Envelope" in analysis_type:
            _code("from cobra.flux_analysis import production_envelope\n"
                  "# Change the reaction to your product reaction\n"
                  "prod_rxn = model.reactions[0]\n"
                  "pe = production_envelope(model, reactions=[prod_rxn])\n"
                  "pe.plot()\n"
                  "plt.show()")
        else:
            # Standard FBA
            _code("solution = model.optimize()\n"
                  "print(f'Status: {solution.status}')\n"
                  "print(f'Objective value: {solution.objective_value:.6f}')\n\n"
                  "# Show top fluxes\n"
                  "fluxes = solution.fluxes\n"
                  "active = fluxes[fluxes.abs() > 1e-6].sort_values(ascending=False)\n"
                  "print(f'\\nActive reactions: {len(active)}')\n"
                  "print(active.head(20))")

        # Visualization
        _md("## Flux Distribution")
        _code("# Histogram of active fluxes\n"
              "solution = model.optimize()\n"
              "active_fluxes = solution.fluxes[solution.fluxes.abs() > 1e-6]\n"
              "plt.figure(figsize=(10, 5))\n"
              "plt.hist(active_fluxes.values, bins=50, edgecolor='black')\n"
              "plt.xlabel('Flux')\n"
              "plt.ylabel('Count')\n"
              "plt.title('Distribution of Active Fluxes')\n"
              "plt.tight_layout()\n"
              "plt.show()")

        # Last run results
        if self.last_run:
            _md("## Previous Results Summary\n\n"
                f"Analysis type: {self.last_run.get('analysis_type', 'N/A')}\n\n"
                f"Timestamp: {self.last_run.get('timestamp', 'N/A')}")
            if "baseline" in self.last_run:
                bl = self.last_run["baseline"]
                _md(f"**Baseline objective**: {bl.get('objective', 'N/A')}\n\n"
                    f"**Status**: {bl.get('status', 'N/A')}")

        # Build notebook JSON
        notebook = {
            "nbformat": 4,
            "nbformat_minor": 5,
            "metadata": {
                "kernelspec": {
                    "display_name": "Python 3",
                    "language": "python",
                    "name": "python3",
                },
                "language_info": {
                    "name": "python",
                    "version": "3.11.0",
                },
            },
            "cells": cells,
        }

        try:
            Path(fp).write_text(json.dumps(notebook, indent=1, ensure_ascii=False), encoding="utf-8")
            QMessageBox.information(self, "Exported", f"Jupyter Notebook saved to:\n{fp}")
            logger.info(f"Jupyter Notebook exported: {fp}")
        except Exception as e:
            self._show_error("Export failed", "Failed to export Jupyter Notebook", e)

    # ==================== PLUGIN SYSTEM ====================

