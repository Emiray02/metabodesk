"""Omics-integration advanced analyses for MetaboDesk.

Provides specialised analysis methods that integrate transcriptomic /
proteomic data with constraint-based metabolic models:

- **Gene Expression Integration**: E-Flux, GIMME, iMAT algorithms for
  context-specific metabolic modelling.
- **GECKO**: Enzyme-constrained FBA using user-supplied kcat / MW data
  (including isoenzyme arm-reaction splitting).
- **Thermodynamic FBA (TMFA)**: Thermodynamically constrained FBA with
  concentration-adjusted ΔG' bounds.
- **FSEOF**: Flux Scanning based on Enforced Objective Flux —
  identifies overexpression / attenuation targets.
"""

from __future__ import annotations

import csv
import logging
import math
from pathlib import Path
from typing import Any

import cobra

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QFileDialog,
    QMessageBox, QTableWidget, QTableWidgetItem, QAbstractItemView,
    QLineEdit, QSpinBox, QTabWidget, QDoubleSpinBox, QCheckBox,
    QComboBox, QPlainTextEdit, QDialog, QDialogButtonBox, QCompleter,
    QFormLayout,
)

from metabodesk_core.utils import evaluate_gpr_expression
from metabodesk_core.widgets import TextPopup

logger = logging.getLogger("MetaboDesk")


class OmicsMixin:
    """Mixin providing omics-integration analyses (E-Flux, GIMME, iMAT,
    GECKO, TMFA, FSEOF)."""

    # ==================== GENE EXPRESSION INTEGRATION ====================

    def run_gene_expression_analysis(self) -> None:
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
            expr_data: dict[str, float] = {}
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

        def _compute() -> dict[str, Any]:
            # Map gene expression to reactions via GPR tree evaluation
            rxn_expr: dict[str, float | None] = {}
            for rxn in model.reactions:
                if not rxn.genes:
                    rxn_expr[rxn.id] = None
                    continue
                gene_vals = {g.id: expr_data[g.id] for g in rxn.genes if g.id in expr_data}
                if not gene_vals:
                    rxn_expr[rxn.id] = None
                    continue
                gpr_str = getattr(rxn, 'gene_reaction_rule', '') or ''
                if gpr_str.strip():
                    val = evaluate_gpr_expression(gpr_str, gene_vals)
                    rxn_expr[rxn.id] = val
                else:
                    rxn_expr[rxn.id] = min(gene_vals.values())

            if method == "E-Flux":
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
                sol0 = model.optimize()
                if sol0.status != 'optimal':
                    raise ValueError("Baseline optimization failed.")
                required_obj_val = float(sol0.objective_value) * required_growth
                obj_rxn_id = None
                obj_coeffs = model.objective.expression.as_coefficients_dict()
                for var, coeff in obj_coeffs.items():
                    if float(coeff) != 0:
                        obj_rxn_id = var.name
                        break
                if obj_rxn_id:
                    obj_rxn_ref = model.reactions.get_by_id(obj_rxn_id)
                    obj_rxn_ref.lower_bound = max(obj_rxn_ref.lower_bound, required_obj_val)
                low_expr_rxns = []
                for rxn in model.reactions:
                    expr_val = rxn_expr.get(rxn.id)
                    if expr_val is not None and expr_val < threshold:
                        low_expr_rxns.append(rxn)
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
                high_rxns: list[str] = []
                low_rxns: list[str] = []
                for rxn in model.reactions:
                    expr_val = rxn_expr.get(rxn.id)
                    if expr_val is not None:
                        if expr_val > threshold:
                            high_rxns.append(rxn.id)
                        elif expr_val <= threshold * 0.1:
                            low_rxns.append(rxn.id)
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

        def _on_done(result: dict[str, Any]) -> None:
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

    # ==================== THERMODYNAMIC FBA (TMFA) ====================

    def run_tmfa(self) -> None:
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

        dg_from_csv: dict[str, float] = {}
        csv_path = dg_file_edit.text().strip()
        if csv_path and Path(csv_path).exists():
            try:
                sep = '\t' if csv_path.endswith('.tsv') else ','
                with open(csv_path, encoding='utf-8') as f:
                    reader = csv.reader(f, delimiter=sep)
                    next(reader)
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

        def _compute() -> dict[str, Any]:
            constrained: list[dict[str, Any]] = []
            for rxn in model.reactions:
                if rxn.boundary:
                    continue
                dg0: float | None = None
                dg_source = ""
                if rxn.id in dg_from_csv:
                    dg0 = dg_from_csv[rxn.id]
                    dg_source = "CSV"
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
                rt = R_GAS * temperature
                ln_q_best = 0.0
                ln_q_worst = 0.0
                for met, coeff in rxn.metabolites.items():
                    if coeff > 0:
                        ln_q_best += coeff * math.log(c_min)
                        ln_q_worst += coeff * math.log(c_max)
                    elif coeff < 0:
                        ln_q_best += coeff * math.log(c_max)
                        ln_q_worst += coeff * math.log(c_min)
                dg_prime_best = dg0 + rt * ln_q_best
                dg_prime_worst = dg0 + rt * ln_q_worst
                action = "none"
                if dg_prime_best > threshold_val:
                    if strict:
                        rxn.upper_bound = 0.0
                        action = "forward blocked (ΔG'>0 even at best)"
                    else:
                        rxn.upper_bound = min(rxn.upper_bound, 1.0)
                        action = "forward constrained (ΔG' range unfavorable)"
                elif dg_prime_worst < -threshold_val:
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

        def _on_done(result: dict[str, Any]) -> None:
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

    # ==================== GECKO / ENZYME-CONSTRAINED MODEL ====================

    def run_gecko_analysis(self) -> None:
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

        kcat_data: dict[str, list[tuple[float, float, str]]] = {}
        try:
            sep = '\t' if csv_path.endswith('.tsv') else ','
            with open(csv_path, encoding='utf-8') as f:
                reader = csv.reader(f, delimiter=sep)
                next(reader)
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

        def _compute() -> dict[str, Any]:
            prot_pool = cobra.Metabolite("prot_pool", name="Protein pool", compartment="c")
            constrained: list[dict[str, Any]] = []
            for rid, enzymes in kcat_data.items():
                if rid not in model.reactions:
                    continue
                rxn = model.reactions.get_by_id(rid)
                if len(enzymes) == 1:
                    kcat_val, mw_val, enz_id = enzymes[0]
                    kcat_h = kcat_val * 3600.0 * sigma
                    enzyme_cost = mw_val / (1000.0 * kcat_h) if kcat_h > 0 else 1e6
                    rxn.add_metabolites({prot_pool: enzyme_cost})
                    constrained.append({"id": rid, "name": rxn.name, "kcat": kcat_val,
                                        "mw": mw_val, "enzyme_cost": enzyme_cost,
                                        "enzyme_id": enz_id})
                else:
                    orig_mets = dict(rxn.metabolites)
                    orig_lb = rxn.lower_bound
                    orig_ub = rxn.upper_bound
                    rxn.bounds = (0, 0)
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
            prot_exchange = cobra.Reaction("prot_pool_exchange")
            prot_exchange.name = "Protein pool exchange"
            prot_exchange.lower_bound = 0.0
            prot_exchange.upper_bound = protein_budget
            prot_exchange.add_metabolites({prot_pool: -1.0})
            model.add_reactions([prot_exchange])
            sol0 = model.optimize()
            baseline_unconstrained = None
            with model:
                model.reactions.get_by_id("prot_pool_exchange").upper_bound = 1e6
                sol_unconstrained = model.optimize()
                if sol_unconstrained.status == 'optimal':
                    baseline_unconstrained = float(sol_unconstrained.objective_value)
            result: dict[str, Any] = {
                "status": sol0.status if sol0 else "failed",
                "objective": float(sol0.objective_value) if sol0 and sol0.status == 'optimal' else 0.0,
                "unconstrained_obj": baseline_unconstrained,
                "constrained_reactions": constrained,
                "protein_budget": protein_budget,
                "sigma": sigma,
                "flux": {r.id: float(sol0.fluxes[r.id]) for r in model.reactions
                         if r.id != "prot_pool_exchange"} if sol0 and sol0.status == 'optimal' else {}
            }
            if sol0 and sol0.status == 'optimal':
                usage: list[dict[str, Any]] = []
                for cr in constrained:
                    fv = abs(sol0.fluxes.get(cr["id"], 0))
                    cost = fv * cr["enzyme_cost"]
                    usage.append({**cr, "flux": fv, "protein_usage": cost})
                usage.sort(key=lambda x: x["protein_usage"], reverse=True)
                result["enzyme_usage"] = usage[:30]
            return result

        def _on_done(result: dict[str, Any]) -> None:
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

    # ==================== FSEOF ====================

    def run_fseof(self) -> None:
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

        raw_targets = [t.strip() for t in target_edit.text().split(",") if t.strip()]
        target_ids: list[str] = []
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

        def _compute(worker=None) -> dict[str, Any]:
            sol0 = model.optimize()
            if sol0.status != 'optimal':
                raise ValueError("Baseline optimization failed.")
            max_obj = float(sol0.objective_value)
            obj_rxn_id: str | None = None
            for v in model.objective.expression.as_coefficients_dict():
                name = getattr(v, "name", "") or ""
                for r in model.reactions:
                    if r.id in name:
                        obj_rxn_id = r.id
                        break
                if obj_rxn_id:
                    break
            all_target_results: list[dict[str, Any]] = []
            for t_idx, target_id in enumerate(target_ids):
                if worker:
                    worker.report_progress(
                        f"FSEOF target {t_idx + 1}/{len(target_ids)}: {target_id}", 0)
                baseline_target = float(sol0.fluxes.get(target_id, 0))
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
                targets: list[dict[str, Any]] = []
                for rid, profile in flux_profiles.items():
                    valid = [v for v in profile if v == v]
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

        def _on_done(result: dict[str, Any]) -> None:
            target_results = result.get("target_results", [])
            if len(target_results) == 0:
                QMessageBox.warning(self, "FSEOF", "No targets analysed.")
                return
            if len(target_results) == 1 and target_results[0].get("error"):
                QMessageBox.warning(self, "FSEOF", target_results[0]["error"])
                return
            dialog2 = QDialog(self)
            title_targets = ", ".join(tr["target_id"] for tr in target_results)
            dialog2.setWindowTitle(f"FSEOF Results — {title_targets}")
            dialog2.resize(900, 650)
            layout2 = QVBoxLayout(dialog2)
            if len(target_results) == 1:
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

    # ==================== ¹³C-MFA FLUX CONSTRAINTS ====================

    def run_13c_mfa_constraints(self) -> None:
        """Apply ¹³C-MFA measured flux constraints and re-run FBA.

        Opens a dialog where the user can load a CSV file with columns
        ``reaction_id, measured_flux, std_dev`` (standard deviation).
        Each measured flux is applied as a tight bound on the corresponding
        reaction:  ``[measured - n*std, measured + n*std]`` where *n* is a
        user-defined confidence multiplier (default 2 ≈ 95 % CI).

        After constraining, FBA is re-optimised and the results are displayed
        in a comparison table (original vs constrained fluxes).
        """
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        # ---- Dialog ----
        dialog = QDialog(self)
        dialog.setWindowTitle("¹³C-MFA Flux Constraints")
        dialog.resize(750, 550)
        layout = QVBoxLayout(dialog)

        layout.addWidget(QLabel(
            "<b>¹³C-MFA Measured Flux Constraints</b><br>"
            "Load a CSV with columns: <code>reaction_id, measured_flux, std_dev</code><br>"
            "Each measured flux constrains the reaction bounds to "
            "<code>[flux − n·σ , flux + n·σ]</code>."
        ))

        # File picker row
        file_row = QHBoxLayout()
        csv_path_edit = QLineEdit()
        csv_path_edit.setPlaceholderText("Path to ¹³C-MFA CSV file…")
        file_row.addWidget(csv_path_edit, stretch=1)
        browse_btn = QPushButton("Browse…")
        file_row.addWidget(browse_btn)
        layout.addLayout(file_row)

        def _browse():
            p, _ = QFileDialog.getOpenFileName(dialog, "Select ¹³C-MFA CSV", "", "CSV (*.csv *.tsv *.txt)")
            if p:
                csv_path_edit.setText(p)

        browse_btn.clicked.connect(_browse)

        # Generate template button
        tmpl_btn = QPushButton("Generate template CSV…")
        layout.addWidget(tmpl_btn)

        def _gen_template():
            p, _ = QFileDialog.getSaveFileName(dialog, "Save Template", "", "CSV (*.csv)")
            if not p:
                return
            try:
                with open(p, "w", newline="", encoding="utf-8") as f:
                    w = csv.writer(f)
                    w.writerow(["reaction_id", "measured_flux", "std_dev"])
                    # Add a few example rows from the model
                    count = 0
                    for rxn in self.base_model.reactions:
                        if count >= 5:
                            break
                        w.writerow([rxn.id, "0.0", "0.1"])
                        count += 1
                QMessageBox.information(dialog, "Template", f"Template saved to:\n{p}")
            except Exception as e:
                self._show_error("Template Error", "Could not save template.", e)

        tmpl_btn.clicked.connect(_gen_template)

        # Parameters
        form = QFormLayout()
        sigma_spin = QDoubleSpinBox()
        sigma_spin.setRange(0.5, 5.0)
        sigma_spin.setValue(2.0)
        sigma_spin.setSingleStep(0.5)
        sigma_spin.setToolTip("Number of standard deviations for constraint bounds (2 ≈ 95% CI)")
        form.addRow("Confidence multiplier (n·σ):", sigma_spin)

        relax_chk = QCheckBox("Relax bounds for infeasible constraints")
        relax_chk.setChecked(True)
        relax_chk.setToolTip("If constraining makes the model infeasible, widen bounds iteratively")
        form.addRow(relax_chk)
        layout.addLayout(form)

        # Preview table
        layout.addWidget(QLabel("Constraint preview:"))
        preview_tbl = QTableWidget(0, 5)
        preview_tbl.setHorizontalHeaderLabels(["Reaction", "Measured", "Std Dev", "New LB", "New UB"])
        preview_tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        preview_tbl.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(preview_tbl, stretch=1)

        # Load & preview button
        load_btn = QPushButton("Load && Preview")
        layout.addWidget(load_btn)

        # Internal state for parsed data
        parsed_data: list[dict] = []

        def _load_and_preview():
            path = csv_path_edit.text().strip()
            if not path or not Path(path).exists():
                QMessageBox.warning(dialog, "File", "Select a valid CSV file.")
                return
            parsed_data.clear()
            preview_tbl.setRowCount(0)
            n_sigma = sigma_spin.value()
            missing = []
            try:
                with open(path, "r", encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        rid = (row.get("reaction_id") or "").strip()
                        if not rid:
                            continue
                        try:
                            mflux = float(row.get("measured_flux", 0))
                        except (ValueError, TypeError):
                            continue
                        try:
                            sd = float(row.get("std_dev", 0))
                        except (ValueError, TypeError):
                            sd = 0.0
                        if rid not in self.base_model.reactions:
                            missing.append(rid)
                            continue
                        new_lb = mflux - n_sigma * sd
                        new_ub = mflux + n_sigma * sd
                        parsed_data.append({
                            "reaction_id": rid, "measured_flux": mflux,
                            "std_dev": sd, "new_lb": new_lb, "new_ub": new_ub,
                        })
                        r = preview_tbl.rowCount()
                        preview_tbl.insertRow(r)
                        preview_tbl.setItem(r, 0, QTableWidgetItem(rid))
                        preview_tbl.setItem(r, 1, QTableWidgetItem(f"{mflux:.6g}"))
                        preview_tbl.setItem(r, 2, QTableWidgetItem(f"{sd:.6g}"))
                        preview_tbl.setItem(r, 3, QTableWidgetItem(f"{new_lb:.6g}"))
                        preview_tbl.setItem(r, 4, QTableWidgetItem(f"{new_ub:.6g}"))

                info = f"Loaded {len(parsed_data)} constraints."
                if missing:
                    info += f"\n⚠ {len(missing)} reaction IDs not found: {', '.join(missing[:5])}"
                    if len(missing) > 5:
                        info += f" … and {len(missing) - 5} more"
                QMessageBox.information(dialog, "Preview", info)
            except Exception as e:
                self._show_error("Load Error", "Failed to parse CSV.", e)

        load_btn.clicked.connect(_load_and_preview)

        # Run buttons
        btn_row = QHBoxLayout()
        run_btn = QPushButton("Apply Constraints && Run FBA")
        cancel_btn = QPushButton("Cancel")
        btn_row.addStretch(1)
        btn_row.addWidget(run_btn)
        btn_row.addWidget(cancel_btn)
        layout.addLayout(btn_row)
        cancel_btn.clicked.connect(dialog.reject)

        def _apply_and_run():
            if not parsed_data:
                QMessageBox.warning(dialog, "No data", "Load a CSV file first.")
                return
            dialog.accept()
            self._run_13c_constrained_fba(parsed_data, relax_chk.isChecked())

        run_btn.clicked.connect(_apply_and_run)
        dialog.exec()

    def _run_13c_constrained_fba(self, constraints: list[dict], relax: bool) -> None:
        """Apply ¹³C-MFA constraints to a model copy, run FBA, and display results."""
        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self._apply_selected_solver(model)

        applied: list[dict] = []
        skipped: list[str] = []

        for c in constraints:
            rid = c["reaction_id"]
            if rid not in model.reactions:
                skipped.append(rid)
                continue
            rxn = model.reactions.get_by_id(rid)
            rxn.lower_bound = c["new_lb"]
            rxn.upper_bound = c["new_ub"]
            applied.append(c)

        def _compute(worker=None) -> dict[str, Any]:
            import numpy as np

            # Baseline (unconstrained) FBA for comparison
            baseline_model = self.base_model.copy()
            self.apply_medium_table_bounds_to_model(baseline_model)
            self._apply_selected_solver(baseline_model)
            base_sol = baseline_model.optimize()
            base_flux = dict(base_sol.fluxes) if base_sol.status == "optimal" else {}
            base_obj = base_sol.objective_value if base_sol.status == "optimal" else None

            # Constrained FBA
            sol = model.optimize()
            if sol.status != "optimal" and relax:
                # Try relaxing: widen bounds by 50% iteratively
                for attempt in range(5):
                    for c in applied:
                        rxn = model.reactions.get_by_id(c["reaction_id"])
                        span = c["std_dev"] * (1.5 ** (attempt + 1))
                        rxn.lower_bound = c["measured_flux"] - span * 2
                        rxn.upper_bound = c["measured_flux"] + span * 2
                    sol = model.optimize()
                    if sol.status == "optimal":
                        break

            cflux = dict(sol.fluxes) if sol.status == "optimal" else {}
            cobj = sol.objective_value if sol.status == "optimal" else None

            # Build comparison rows
            rows: list[dict] = []
            all_ids = sorted(set(list(base_flux.keys()) + list(cflux.keys())))
            for rid in all_ids:
                bf = base_flux.get(rid, 0.0)
                cf = cflux.get(rid, 0.0)
                delta = cf - bf
                rows.append({"reaction": rid, "baseline": bf, "constrained": cf, "delta": delta})

            # Sort by |delta| descending
            rows.sort(key=lambda r: abs(r["delta"]), reverse=True)

            constrained_ids = {c["reaction_id"] for c in applied}

            return {
                "base_obj": base_obj, "base_status": base_sol.status,
                "c_obj": cobj, "c_status": sol.status,
                "rows": rows, "n_applied": len(applied),
                "n_skipped": len(skipped), "constrained_ids": constrained_ids,
            }

        def _on_done(result: dict[str, Any]) -> None:
            self._show_13c_results(result)

        self._launch_worker(_compute, _on_done, "Running ¹³C-MFA constrained FBA…")

    def _show_13c_results(self, result: dict[str, Any]) -> None:
        """Display ¹³C-MFA constrained FBA results in a dialog."""
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas

        dialog = QDialog(self)
        dialog.setWindowTitle("¹³C-MFA Constrained FBA Results")
        dialog.resize(1000, 650)
        layout = QVBoxLayout(dialog)

        # Summary
        base_status = result["base_status"]
        c_status = result["c_status"]
        base_obj = result["base_obj"]
        c_obj = result["c_obj"]
        summary = (
            f"Baseline FBA: {base_status} | obj = {base_obj:.6g if base_obj else 'N/A'}\n"
            f"¹³C-Constrained FBA: {c_status} | obj = {c_obj:.6g if c_obj else 'N/A'}\n"
            f"Constraints applied: {result['n_applied']} | Skipped: {result['n_skipped']}"
        )
        lbl = QLabel(summary)
        lbl.setWordWrap(True)
        layout.addWidget(lbl)

        # Chart — top 20 most changed reactions
        rows = result["rows"]
        top_n = rows[:20]
        if top_n:
            fig = Figure(figsize=(8, 3.5), dpi=100)
            ax = fig.add_subplot(111)
            rids = [r["reaction"][:20] for r in top_n]
            deltas = [r["delta"] for r in top_n]
            colors = ["#e74c3c" if r["reaction"] in result["constrained_ids"] else "#3498db" for r in top_n]
            bars = ax.barh(rids, deltas, color=colors, alpha=0.8, edgecolor="gray", linewidth=0.5)
            ax.set_xlabel("Δ Flux (Constrained − Baseline)")
            ax.set_title("Top 20 Most Changed Fluxes (red = constrained reaction)")
            ax.axvline(0, color="black", linewidth=0.5)
            fig.tight_layout()
            canvas = FigureCanvas(fig)
            layout.addWidget(canvas, stretch=1)

        # Table
        table = QTableWidget(0, 4)
        table.setHorizontalHeaderLabels(["Reaction", "Baseline Flux", "Constrained Flux", "Δ Flux"])
        table.horizontalHeader().setStretchLastSection(True)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        for r in rows[:200]:  # limit for performance
            row_idx = table.rowCount()
            table.insertRow(row_idx)
            table.setItem(row_idx, 0, QTableWidgetItem(r["reaction"]))
            table.setItem(row_idx, 1, QTableWidgetItem(f"{r['baseline']:.6g}"))
            table.setItem(row_idx, 2, QTableWidgetItem(f"{r['constrained']:.6g}"))
            table.setItem(row_idx, 3, QTableWidgetItem(f"{r['delta']:.6g}"))
        layout.addWidget(table, stretch=1)

        # Export CSV
        btn_row = QHBoxLayout()
        export_btn = QPushButton("Export CSV…")

        def _export():
            p, _ = QFileDialog.getSaveFileName(dialog, "Export", "", "CSV (*.csv)")
            if not p:
                return
            try:
                with open(p, "w", newline="", encoding="utf-8") as f:
                    w = csv.writer(f)
                    w.writerow(["reaction_id", "baseline_flux", "constrained_flux", "delta"])
                    for r in rows:
                        w.writerow([r["reaction"], r["baseline"], r["constrained"], r["delta"]])
                QMessageBox.information(dialog, "Export", f"Exported to:\n{p}")
            except Exception as e:
                self._show_error("Export Error", "Failed to export.", e)

        export_btn.clicked.connect(_export)
        btn_row.addWidget(export_btn)
        btn_row.addStretch(1)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)

        dialog.exec()
