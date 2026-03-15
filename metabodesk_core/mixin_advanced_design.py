"""Strain-design & model-engineering advanced analyses for MetaboDesk.

Provides analysis methods for rational strain design, model
diagnostics, and data export:

- **OptKnock / Strain Design**: Combinatorial knockout search to
  maximise target-metabolite production.
- **Dynamic FBA (dFBA)**: Time-course batch-culture simulation
  with multi-substrate support (Euler / RK4 integration).
- **Flux Coupling Analysis**: Ratio-based coupling detection.
- **Gap-Filling**: Identifies missing reactions via LP / iterative search.
- **Model Comparison**: Side-by-side diff of two SBML models.
- **Pathway Enrichment**: Subsystem-based enrichment of active fluxes.
- **Escher Map Viewer**: JSON-based or package-based pathway viewer.
- **Phenotype Phase Plane (PhPP)**: 2-D objective scanning.
- **SBML Validation**: libSBML or fallback consistency checks.
- **Jupyter Notebook Export**: Reproducible analysis notebooks.
- **Minimal Cut Sets (MCS)**: Enumerates minimal reaction knockouts
  that eliminate target flux while preserving growth.
- **Metabolic Distance**: BFS shortest-path distance in the bipartite
  metabolite–reaction graph with currency-metabolite exclusion.
- **Pareto Frontier**: Multi-objective optimisation sweep between two
  reactions, producing the Pareto-optimal trade-off curve with
  knee-point detection and CSV export.
"""

from __future__ import annotations

import csv
import json
import logging
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Any

import cobra

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QFileDialog,
    QMessageBox, QTableWidget, QTableWidgetItem, QAbstractItemView,
    QLineEdit, QSpinBox, QTabWidget, QDoubleSpinBox, QCheckBox,
    QComboBox, QPlainTextEdit, QDialog, QDialogButtonBox, QCompleter,
    QFormLayout, QStackedWidget, QWidget,
)

from metabodesk_core.widgets import TextPopup


def _lazy_mpl_canvas():
    """Import matplotlib on demand (avoids ~0.4 s at module load)."""
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    return FigureCanvas, Figure

logger = logging.getLogger("MetaboDesk")


class DesignMixin:
    """Mixin providing strain-design, model diagnostics, and export
    functionality (OptKnock, dFBA, Flux Coupling, Gap-Filling, PhPP,
    Model Comparison, Pathway Enrichment, Escher, SBML Validation,
    Jupyter Export)."""

    # ==================== OPTKNOCK / STRAIN DESIGN ====================

    def run_optknock(self) -> None:
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

        def _compute(worker=None) -> dict[str, Any]:
            obj_rxns: set[str] = set()
            try:
                for v in model.objective.expression.as_coefficients_dict():
                    for r in model.reactions:
                        if r.id in str(v):
                            obj_rxns.add(r.id)
            except Exception:
                pass
            candidates = [r.id for r in model.reactions
                          if not r.boundary and r.id not in obj_rxns and r.id != target_id]
            with model:
                sol0 = model.optimize()
                if sol0.status != 'optimal':
                    raise ValueError("Baseline optimization failed.")
                baseline_obj = float(sol0.objective_value)
                baseline_target = float(sol0.fluxes.get(target_id, 0.0))
                active = [(rid, abs(sol0.fluxes.get(rid, 0)))
                          for rid in candidates if abs(sol0.fluxes.get(rid, 0)) > 1e-8]
                active.sort(key=lambda x: x[1], reverse=True)
            max_single = min(len(active), max(200, len(active)))
            candidates = [rid for rid, _ in active[:max_single]]

            results: list[dict[str, Any]] = []
            total = 0
            limit_single = len(candidates)
            limit_double = min(len(candidates), 120)
            limit_triple = min(len(candidates), 40)
            for k in range(1, max_ko + 1):
                limit = limit_single if k == 1 else limit_double if k == 2 else limit_triple
                total += sum(1 for _ in combinations(candidates[:limit], k))
            done = 0

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
                            if growth > baseline_obj * 0.01:
                                results.append({
                                    "knockouts": [rid], "target_flux": target_flux,
                                    "growth": growth, "growth_ratio": growth / baseline_obj
                                })
                    except Exception:
                        pass

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

        def _on_done(result: dict[str, Any]) -> None:
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

    def run_dfba(self) -> None:
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

        def _solve_fba_at_state(
            biomass_val: float, substrate_val: float,
            secondary_vals: list[float], dt_val: float,
        ) -> tuple[float, dict[str, float]]:
            zero_rates: dict[str, float] = {
                substrate_rxn_id: 0.0,
                **{s[0]: 0.0 for s in secondary_substrates},
            }
            if biomass_val <= 1e-9:
                return 0.0, zero_rates
            with model:
                if substrate_val <= 1e-6:
                    model.reactions.get_by_id(substrate_rxn_id).lower_bound = 0.0
                else:
                    max_uptake = substrate_val / (biomass_val * dt_val)
                    sub_rxn = model.reactions.get_by_id(substrate_rxn_id)
                    sub_rxn.lower_bound = max(-max_uptake, sub_rxn.lower_bound)
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
                with model:
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
                        return 0.0, uptake_rates
                return 0.0, zero_rates

        def _compute(worker=None) -> dict[str, Any]:
            times: list[float] = []
            biomasses: list[float] = []
            substrates: list[float] = []
            growths: list[float] = []
            secondary_concs = [s[1] for s in secondary_substrates]
            secondary_traces: list[list[float]] = [[] for _ in secondary_substrates]
            biomass = init_biomass
            substrate = init_substrate
            t = 0.0
            step = 0

            while t < t_end and biomass > 1e-9:
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
                    growth_rate, uptake_rates = _solve_fba_at_state(biomass, substrate, secondary_concs, dt)
                    growths.append(growth_rate)
                    old_biomass = biomass
                    biomass += growth_rate * biomass * dt
                    substrate += uptake_rates[substrate_rxn_id] * old_biomass * dt
                    for i, (sec_id, _) in enumerate(secondary_substrates):
                        secondary_concs[i] += uptake_rates.get(sec_id, 0) * old_biomass * dt

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

        def _on_done(result: dict[str, Any]) -> None:
            dialog2 = QDialog(self)
            dialog2.setWindowTitle(f"dFBA Results ({result.get('method', '')})")
            dialog2.resize(900, 700)
            layout2 = QVBoxLayout(dialog2)
            n_subs = len(result.get("secondary_substrates", []))
            FigureCanvas, Figure = _lazy_mpl_canvas()
            fig = Figure(figsize=(8, 6), dpi=100)
            ax1 = fig.add_subplot(211)
            ax1.plot(result["times"], result["biomasses"], 'g-', linewidth=2, label="Biomass (g/L)")
            ax1.plot(result["times"], result["substrates"], 'r--', linewidth=2,
                     label=f"{result['substrate_rxn']} (mmol/L)")
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

    def run_flux_coupling(self) -> None:
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

        def _compute() -> dict[str, Any]:
            from cobra.flux_analysis import flux_variability_analysis
            fva = flux_variability_analysis(model, fraction_of_optimum=0.0)
            active_rxns: list[str] = []
            fva_ranges: dict[str, tuple[float, float]] = {}
            for rxn_id, row in fva.iterrows():
                fmin, fmax = float(row['minimum']), float(row['maximum'])
                if abs(fmax) > tol or abs(fmin) > tol:
                    active_rxns.append(rxn_id)
                    fva_ranges[rxn_id] = (fmin, fmax)
            test_rxns = active_rxns[:max_rxns]
            blocked = [rid for rid in fva.index
                       if abs(fva.loc[rid, 'minimum']) < tol and abs(fva.loc[rid, 'maximum']) < tol]
            coupled_pairs: list[tuple[str, str, str, str]] = []
            ratio_groups: dict[tuple[float, float], list[str]] = {}
            for rid in test_rxns:
                fmin, fmax = fva_ranges[rid]
                key = (round(fmin, 6), round(fmax, 6))
                ratio_groups.setdefault(key, []).append(rid)
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
                        ratio_min = fmin1 / fmin2 if abs(fmin2) > tol else None
                        ratio_max = fmax1 / fmax2 if abs(fmax2) > tol else None
                        if ratio_min is not None and ratio_max is not None:
                            if abs(ratio_min - ratio_max) < tol * 100:
                                coupled_pairs.append((rid1, rid2, "Fully coupled",
                                                      f"ratio ≈ {ratio_min:.4g}"))
                            elif abs(ratio_min) > tol and abs(ratio_max) > tol:
                                coupled_pairs.append((rid1, rid2, "Partially coupled",
                                                      f"ratio {ratio_min:.4g} — {ratio_max:.4g}"))
            directional: list[tuple[str, str, str, str]] = []
            checked_groups = [g for g in ratio_groups.values() if len(g) == 1]
            singleton_rxns = [g[0] for g in checked_groups][:50]
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
                            new_min = float(fva2.loc[rid2, 'minimum'])
                            new_max = float(fva2.loc[rid2, 'maximum'])
                            if (abs(new_min) < tol and abs(new_max) < tol and
                                    (abs(orig_range[0]) > tol or abs(orig_range[1]) > tol)):
                                directional.append((rid, rid2, "Directional", f"{rid} → {rid2}"))
                    except Exception:
                        pass
            all_results = coupled_pairs + directional
            return {"coupled": all_results, "active_reactions": len(active_rxns),
                    "blocked_reactions": len(blocked), "tested": len(test_rxns),
                    "ratio_groups": len(ratio_groups)}

        def _on_done(result: dict[str, Any]) -> None:
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

    # ==================== GAP-FILLING ====================

    def run_gap_filling(self) -> None:
        """Gap-Filling — identify reactions to add so the model can produce biomass."""
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

        def _compute(worker=None) -> dict[str, Any]:
            sol0 = model.optimize()
            baseline_obj = float(sol0.objective_value) if sol0.status == 'optimal' else 0.0
            universal_rxns: list[cobra.Reaction] = []
            univ_model = None
            if univ_path and Path(univ_path).exists():
                univ_model = cobra.io.read_sbml_model(univ_path)
            if univ_model is not None:
                try:
                    from cobra.flux_analysis.gapfilling import gapfill as cobra_gapfill
                    if worker:
                        worker.report_progress("Running COBRApy gapfill (LP)…", 20)
                    gf_results = cobra_gapfill(
                        model, universal=univ_model, lower_bound=0.05,
                        demand_reactions=False, exchange_reactions=False, iterations=1)
                    added_reactions: list[dict[str, Any]] = []
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
                    existing_ids = {r.id for r in model.reactions}
                    for rxn in univ_model.reactions:
                        if rxn.id not in existing_ids and not rxn.boundary:
                            universal_rxns.append(rxn)
            else:
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
                        try:
                            new_rxn = cobra.Reaction(rxn.id)
                            new_rxn.name = rxn.name
                            new_rxn.lower_bound = rxn.lower_bound
                            new_rxn.upper_bound = rxn.upper_bound
                            met_dict: dict[cobra.Metabolite, float] = {}
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
                    break
                new_rxn = cobra.Reaction(best_rxn.id)
                new_rxn.name = best_rxn.name
                new_rxn.lower_bound = best_rxn.lower_bound
                new_rxn.upper_bound = best_rxn.upper_bound
                met_dict2: dict[cobra.Metabolite, float] = {}
                for met, coeff in best_rxn.metabolites.items():
                    if met.id in model.metabolites:
                        met_dict2[model.metabolites.get_by_id(met.id)] = coeff
                    else:
                        new_met = cobra.Metabolite(met.id, name=met.name,
                                                   compartment=met.compartment,
                                                   formula=met.formula)
                        met_dict2[new_met] = coeff
                new_rxn.add_metabolites(met_dict2)
                model.add_reactions([new_rxn])
                current_obj = best_obj
                added_reactions.append({
                    "id": best_rxn.id, "name": best_rxn.name or "",
                    "equation": best_rxn.build_reaction_string(),
                    "objective_after": best_obj
                })
                if current_obj > 1e-6:
                    break
            return {"baseline_obj": baseline_obj, "added": added_reactions,
                    "new_obj": current_obj, "candidates_tested": len(universal_rxns),
                    "method": "Iterative (greedy)"}

        def _on_done(result: dict[str, Any]) -> None:
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

    # ==================== MODEL COMPARISON ====================

    def run_model_comparison(self) -> None:
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

        rxns1 = set(r.id for r in m1.reactions)
        rxns2 = set(r.id for r in m2.reactions)
        only_in_1 = sorted(rxns1 - rxns2)
        only_in_2 = sorted(rxns2 - rxns1)
        common_rxns = sorted(rxns1 & rxns2)
        bound_diffs: list[tuple[str, float, float, float, float]] = []
        for rid in common_rxns:
            r1 = m1.reactions.get_by_id(rid)
            r2 = m2.reactions.get_by_id(rid)
            if r1.lower_bound != r2.lower_bound or r1.upper_bound != r2.upper_bound:
                bound_diffs.append((rid, r1.lower_bound, r1.upper_bound, r2.lower_bound, r2.upper_bound))

        mets1 = set(m.id for m in m1.metabolites)
        mets2 = set(m.id for m in m2.metabolites)
        mets_only1 = sorted(mets1 - mets2)
        mets_only2 = sorted(mets2 - mets1)

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
            f"Unique to {name2}: {len(only_in_2)} rxns, {len(mets_only2)} mets, {len(genes_only2)} genes")
        summary.setWordWrap(True)
        layout.addWidget(summary)

        tabs = QTabWidget()
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

    def run_pathway_enrichment(self) -> None:
        """Metabolic Pathway Enrichment — subsystem-based enrichment of active reactions."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return
        if self.last_run is None:
            QMessageBox.warning(self, "No results", "Run an FBA analysis first to determine active reactions.")
            return

        model = self.base_model
        flux: dict[str, float] = {}
        if "baseline" in self.last_run and "flux" in self.last_run["baseline"]:
            flux = self.last_run["baseline"]["flux"]
        elif "flux" in self.last_run:
            flux = self.last_run["flux"]
        if not flux:
            QMessageBox.warning(self, "No flux data", "Last analysis does not contain flux data.")
            return

        threshold = 1e-6
        active_rxns = set(rid for rid, fv in flux.items() if abs(fv) > threshold)

        subsystem_rxns: dict[str, set[str]] = {}
        for rxn in model.reactions:
            ss = getattr(rxn, 'subsystem', '') or ''
            ss = ss.strip()
            if ss:
                subsystem_rxns.setdefault(ss, set()).add(rxn.id)
        if not subsystem_rxns:
            QMessageBox.warning(self, "No subsystems", "Model has no subsystem annotations.")
            return

        total_rxns = len(model.reactions)
        total_active = len(active_rxns)
        enrichment: list[dict[str, Any]] = []
        for ss, ss_rxns in subsystem_rxns.items():
            ss_total = len(ss_rxns)
            ss_active = len(ss_rxns & active_rxns)
            if ss_total == 0:
                continue
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

        FigureCanvas, Figure = _lazy_mpl_canvas()
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

    # ==================== ESCHER MAP ====================

    def show_escher_map(self) -> None:
        """Show Escher metabolic map.  Supports loading Escher JSON map files
        directly (no ``escher`` package required) and falls back to the
        ``escher`` Python package when available."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        flux: dict[str, float] = {}
        if self.last_run and "baseline" in self.last_run and "flux" in self.last_run["baseline"]:
            flux = self.last_run["baseline"]["flux"]

        dialog = QDialog(self)
        dialog.setWindowTitle("Escher Map Viewer")
        dialog.resize(1200, 800)
        layout = QVBoxLayout(dialog)

        toolbar = QHBoxLayout()
        json_btn = QPushButton("Load Escher JSON map...")
        save_html_btn = QPushButton("Save HTML...")
        save_html_btn.setEnabled(False)
        toolbar.addWidget(json_btn)
        toolbar.addWidget(save_html_btn)
        toolbar.addStretch()
        layout.addLayout(toolbar)

        content_stack = QStackedWidget()
        placeholder = QLabel("Load an Escher JSON map file or use the escher package.\n\n"
                             "You can download Escher maps from:\n"
                             "  https://escher.github.io/#/app\n\n"
                             "Click 'Load Escher JSON map...' to import a .json map file.")
        placeholder.setAlignment(Qt.AlignCenter)
        placeholder.setWordWrap(True)
        content_stack.addWidget(placeholder)
        layout.addWidget(content_stack, 1)

        html_holder: list[str] = [""]

        def _build_escher_html(map_json_str: str) -> str:
            import json as _json
            import sys as _sys
            flux_json = _json.dumps(flux) if flux else "null"
            model_json = "null"
            try:
                model_json = cobra.io.to_json(self.base_model)
            except Exception:
                pass
            escher_script_tag = '<script src="https://unpkg.com/escher@1.7.3/dist/escher.min.js"></script>'
            if not hasattr(self, '_escher_js_cache'):
                self._escher_js_cache = None
            if self._escher_js_cache is not None:
                escher_script_tag = f"<script>{self._escher_js_cache}</script>"
            else:
                search_dirs = [
                    Path(getattr(_sys, '_MEIPASS', '.')),
                    Path(__file__).parent.parent,
                    Path(__file__).parent,
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
                import json as _json
                parsed = _json.loads(map_text)
                if not isinstance(parsed, (list, dict)):
                    raise ValueError("Not a valid Escher map JSON.")
                html = _build_escher_html(map_text)
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

    def run_phenotype_phase_plane(self) -> None:
        """Phenotype Phase Plane analysis: 2-D sweep of two reaction fluxes."""
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

        def _compute(worker=None) -> dict[str, Any]:
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

        def _on_done(result: dict[str, Any]) -> None:
            self.last_run = result
            self._render_phpp_results(result)

        self._launch_worker(_compute, _on_done, f"Running PhPP ({rx1_id} × {rx2_id})...")

    def _render_phpp_results(self, result: dict[str, Any]) -> None:
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

        FigureCanvas, Figure = _lazy_mpl_canvas()
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

        valid = grid[~np.isnan(grid)]
        info = QLabel(
            f"Grid: {len(x_vals)}×{len(y_vals)} = {len(x_vals)*len(y_vals)} points  |  "
            f"Feasible: {len(valid)}  |  "
            f"Objective range: [{valid.min():.4f}, {valid.max():.4f}]" if len(valid) > 0
            else f"Grid: {len(x_vals)}×{len(y_vals)}  |  No feasible solutions found"
        )
        layout.addWidget(info)

        btn_row = QHBoxLayout()
        save_png = QPushButton("Save as PNG...")
        def _save_png():
            fp2, _ = QFileDialog.getSaveFileName(dialog, "Save Plot", "", "PNG (*.png)")
            if fp2:
                fig.savefig(fp2, dpi=150)
        save_png.clicked.connect(_save_png)
        btn_row.addWidget(save_png)

        save_csv_btn = QPushButton("Save grid as CSV...")
        def _save_csv():
            fp2, _ = QFileDialog.getSaveFileName(dialog, "Save CSV", "", "CSV (*.csv)")
            if fp2:
                import csv as _csv
                with open(fp2, 'w', newline='', encoding='utf-8') as f:
                    w = _csv.writer(f)
                    w.writerow([""] + [f"{v:.4f}" for v in x_vals])
                    for i2, yv in enumerate(y_vals):
                        w.writerow([f"{yv:.4f}"] + [
                            f"{grid[i2][j2]:.6f}" if not np.isnan(grid[i2][j2]) else ""
                            for j2 in range(len(x_vals))])
        save_csv_btn.clicked.connect(_save_csv)
        btn_row.addWidget(save_csv_btn)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)
        dialog.exec()

    # ==================== SBML VALIDATION ====================

    def validate_sbml(self) -> None:
        """Validate the currently loaded SBML file using libSBML or fallback checks."""
        if self.current_sbml_path is None or not self.current_sbml_path.exists():
            QMessageBox.warning(self, "No SBML file", "Load an SBML model first.")
            return

        errors: list[str] = []
        warnings: list[str] = []
        infos: list[str] = []
        sbml_path = str(self.current_sbml_path)

        try:
            import libsbml
            reader = libsbml.SBMLReader()
            doc = reader.readSBML(sbml_path)
            doc.setConsistencyChecks(libsbml.LIBSBML_CAT_GENERAL_CONSISTENCY, True)
            doc.setConsistencyChecks(libsbml.LIBSBML_CAT_IDENTIFIER_CONSISTENCY, True)
            doc.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, False)
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
                fbc = model_obj.getPlugin("fbc")
                if fbc:
                    infos.append(f"FBC plugin v{fbc.getPackageVersion()}: {fbc.getNumObjectives()} objective(s), {fbc.getNumGeneProducts()} gene product(s)")
                else:
                    warnings.append("FBC plugin not present — flux bounds may be missing")
        except ImportError:
            infos.append("libSBML not installed — performing basic checks only")
            infos.append("Install with: pip install python-libsbml")
            if self.base_model:
                m = self.base_model
                infos.append(f"Model: {m.id}  |  Reactions: {len(m.reactions)}  |  Metabolites: {len(m.metabolites)}  |  Genes: {len(m.genes)}")
                unbalanced = 0
                for r in m.reactions:
                    if not r.boundary:
                        bal = r.check_mass_balance()
                        if bal:
                            unbalanced += 1
                if unbalanced:
                    warnings.append(f"{unbalanced} reaction(s) have mass imbalance")
                no_gpr = sum(1 for r in m.reactions if not r.gene_reaction_rule.strip())
                if no_gpr:
                    infos.append(f"{no_gpr} reaction(s) have no GPR rule")
                orphans = [m2.id for m2 in m.metabolites if len(m2.reactions) == 0]
                if orphans:
                    warnings.append(f"{len(orphans)} orphan metabolite(s): {', '.join(orphans[:10])}")
                dead_end = 0
                for met in m.metabolites:
                    producers = sum(1 for r2 in met.reactions if met in r2.products)
                    consumers = sum(1 for r2 in met.reactions if met in r2.reactants)
                    if (producers > 0 and consumers == 0) or (consumers > 0 and producers == 0):
                        if not any(r2.boundary for r2 in met.reactions):
                            dead_end += 1
                if dead_end:
                    warnings.append(f"{dead_end} potential dead-end metabolite(s)")
                unbounded = sum(1 for r in m.reactions if abs(r.upper_bound) >= 999999 or abs(r.lower_bound) >= 999999)
                if unbounded:
                    infos.append(f"{unbounded} reaction(s) with very large bounds (±999999+)")

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

    def export_jupyter_notebook(self) -> None:
        """Export the current analysis as a Jupyter Notebook (.ipynb)."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        fp, _ = QFileDialog.getSaveFileName(
            self, "Export Jupyter Notebook", "", "Jupyter Notebook (*.ipynb)")
        if not fp:
            return

        cells: list[dict[str, Any]] = []

        def _md(src: str) -> None:
            cells.append({"cell_type": "markdown", "metadata": {}, "source": [src]})

        def _code(src: str) -> None:
            cells.append({"cell_type": "code", "metadata": {},
                          "source": [src], "execution_count": None, "outputs": []})

        _md("# MetaboDesk — Exported Analysis Notebook\n\n"
            f"Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} by MetaboDesk.")
        _code("import cobra\nimport pandas as pd\nimport matplotlib.pyplot as plt\nfrom pathlib import Path")

        model_file = self.current_sbml_path.name if self.current_sbml_path else "model.xml"
        _md(f"## Load Model\n\nMake sure `{model_file}` is in the same directory as this notebook.")
        _code(f'model = cobra.io.read_sbml_model("{model_file}")\n'
              f'print(f"Model: {{model.id}}")\n'
              f'print(f"Reactions: {{len(model.reactions)}}, Metabolites: {{len(model.metabolites)}}, Genes: {{len(model.genes)}}")')

        if self.reaction_bound_overrides:
            _md("## Custom Reaction Bounds")
            lines2 = ["# Apply custom bounds"]
            for rid, (lb, ub) in self.reaction_bound_overrides.items():
                lines2.append(f'model.reactions.get_by_id("{rid}").bounds = ({lb}, {ub})')
            _code("\n".join(lines2))

        if self.knockout_genes:
            _md("## Gene Knockouts")
            genes_str = ", ".join(f'"{g}"' for g in sorted(self.knockout_genes))
            _code(f"# Knockout genes\nfor gid in [{genes_str}]:\n"
                  f"    gene = model.genes.get_by_id(gid)\n"
                  f"    cobra.manipulation.delete_model_genes(model, [gene])")

        if self.base_model:
            obj_rxn = None
            for r in self.base_model.objective.variables:
                obj_rxn = r.name
                break
            if obj_rxn:
                _md("## Set Objective")
                _code(f'model.objective = "{obj_rxn}"')

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
            _code("# Robustness analysis\nimport numpy as np\n\n"
                  "target_rxn = model.reactions[0]  # Change to your target\n"
                  "steps = np.linspace(target_rxn.lower_bound, target_rxn.upper_bound, 20)\n"
                  "results = []\nfor v in steps:\n"
                  "    with model:\n"
                  "        target_rxn.bounds = (v, v)\n"
                  "        sol = model.optimize()\n"
                  "        results.append({'flux': v, 'objective': sol.objective_value if sol.status == 'optimal' else float('nan')})\n"
                  "df = pd.DataFrame(results)\ndf.plot(x='flux', y='objective')\nplt.show()")
        elif "Production Envelope" in analysis_type:
            _code("from cobra.flux_analysis import production_envelope\n"
                  "prod_rxn = model.reactions[0]\n"
                  "pe = production_envelope(model, reactions=[prod_rxn])\n"
                  "pe.plot()\nplt.show()")
        else:
            _code("solution = model.optimize()\n"
                  "print(f'Status: {solution.status}')\n"
                  "print(f'Objective value: {solution.objective_value:.6f}')\n\n"
                  "fluxes = solution.fluxes\n"
                  "active = fluxes[fluxes.abs() > 1e-6].sort_values(ascending=False)\n"
                  "print(f'\\nActive reactions: {len(active)}')\nprint(active.head(20))")

        _md("## Flux Distribution")
        _code("solution = model.optimize()\n"
              "active_fluxes = solution.fluxes[solution.fluxes.abs() > 1e-6]\n"
              "plt.figure(figsize=(10, 5))\n"
              "plt.hist(active_fluxes.values, bins=50, edgecolor='black')\n"
              "plt.xlabel('Flux')\nplt.ylabel('Count')\n"
              "plt.title('Distribution of Active Fluxes')\nplt.tight_layout()\nplt.show()")

        if self.last_run:
            _md("## Previous Results Summary\n\n"
                f"Analysis type: {self.last_run.get('analysis_type', 'N/A')}\n\n"
                f"Timestamp: {self.last_run.get('timestamp', 'N/A')}")
            if "baseline" in self.last_run:
                bl = self.last_run["baseline"]
                _md(f"**Baseline objective**: {bl.get('objective', 'N/A')}\n\n"
                    f"**Status**: {bl.get('status', 'N/A')}")

        notebook = {
            "nbformat": 4, "nbformat_minor": 5,
            "metadata": {
                "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                "language_info": {"name": "python", "version": "3.11.0"},
            },
            "cells": cells,
        }

        try:
            Path(fp).write_text(json.dumps(notebook, indent=1, ensure_ascii=False), encoding="utf-8")
            QMessageBox.information(self, "Exported", f"Jupyter Notebook saved to:\n{fp}")
            logger.info(f"Jupyter Notebook exported: {fp}")
        except Exception as e:
            self._show_error("Export failed", "Failed to export Jupyter Notebook", e)

    # ==================== MINIMAL CUT SETS (MCS) ====================

    def run_minimal_cut_sets(self) -> None:
        """Minimal Cut Sets — find the smallest set of reaction deletions
        that completely block a target reaction.

        Uses an iterative LP-based approach: enumerate single, double, and
        triple knockout combinations that reduce target flux to zero while
        (optionally) preserving a minimum growth rate.
        """
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Minimal Cut Sets (MCS)")
        dialog.resize(500, 320)
        layout = QVBoxLayout(dialog)
        layout.addWidget(QLabel(
            "Find the smallest set of reaction deletions that block\n"
            "flux through a target reaction.\n\n"
            "MCS are useful for identifying metabolic vulnerabilities\n"
            "and essential pathway bottlenecks."))

        form = QFormLayout()
        target_edit = QLineEdit()
        target_edit.setPlaceholderText("Target reaction to block (e.g., EX_etoh_e)")
        if self.base_model:
            completer = QCompleter([r.id for r in self.base_model.reactions])
            completer.setCaseSensitivity(Qt.CaseInsensitive)
            target_edit.setCompleter(completer)
        form.addRow("Target reaction:", target_edit)

        max_size_spin = QSpinBox()
        max_size_spin.setRange(1, 3)
        max_size_spin.setValue(2)
        form.addRow("Max cut-set size:", max_size_spin)

        min_growth_spin = QDoubleSpinBox()
        min_growth_spin.setRange(0.0, 1.0)
        min_growth_spin.setValue(0.01)
        min_growth_spin.setDecimals(3)
        form.addRow("Min growth fraction:", min_growth_spin)

        top_spin = QSpinBox()
        top_spin.setRange(5, 100)
        top_spin.setValue(20)
        form.addRow("Max results:", top_spin)
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
        max_size = max_size_spin.value()
        min_growth = min_growth_spin.value()
        top_n = top_spin.value()

        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self._apply_selected_solver(model)

        def _compute(worker=None) -> dict[str, Any]:
            sol0 = model.optimize()
            if sol0.status != 'optimal':
                raise ValueError("Baseline optimization failed.")
            baseline_obj = float(sol0.objective_value)
            baseline_target = float(sol0.fluxes.get(target_id, 0))

            if abs(baseline_target) < 1e-9:
                return {"cut_sets": [], "baseline_obj": baseline_obj,
                        "baseline_target": baseline_target, "target_id": target_id,
                        "message": "Target reaction already carries no flux."}

            # Get candidate reactions (non-boundary, not the target, not the objective)
            obj_rxns: set[str] = set()
            try:
                for v in model.objective.expression.as_coefficients_dict():
                    for r in model.reactions:
                        if r.id in str(v):
                            obj_rxns.add(r.id)
            except Exception:
                pass

            candidates = [r.id for r in model.reactions
                          if not r.boundary and r.id not in obj_rxns and r.id != target_id]

            # Filter to active reactions only
            active_candidates = [rid for rid in candidates
                                 if abs(sol0.fluxes.get(rid, 0)) > 1e-8]
            active_candidates.sort(key=lambda rid: abs(sol0.fluxes.get(rid, 0)), reverse=True)

            cut_sets: list[dict[str, Any]] = []
            tested = 0
            total_limit = min(len(active_candidates), 200)
            search_pool = active_candidates[:total_limit]

            # Singles
            for rid in search_pool:
                tested += 1
                if worker and tested % 10 == 0:
                    worker.report_progress(f"MCS singles: {tested}", 0)
                with model:
                    model.reactions.get_by_id(rid).knock_out()
                    try:
                        sol = model.optimize()
                        if sol.status == 'optimal':
                            growth = float(sol.objective_value)
                            tgt_flux = float(sol.fluxes.get(target_id, 0))
                            if abs(tgt_flux) < 1e-8 and growth >= baseline_obj * min_growth:
                                cut_sets.append({
                                    "reactions": [rid], "size": 1,
                                    "growth": growth, "growth_pct": growth / baseline_obj * 100,
                                    "target_flux": tgt_flux,
                                })
                    except Exception:
                        pass

            # Doubles
            if max_size >= 2:
                pool2 = search_pool[:min(len(search_pool), 100)]
                total_pairs = len(pool2) * (len(pool2) - 1) // 2
                pair_count = 0
                for i, rid1 in enumerate(pool2):
                    for rid2 in pool2[i+1:]:
                        pair_count += 1
                        if worker and pair_count % 50 == 0:
                            worker.report_progress(
                                f"MCS doubles: {pair_count}/{total_pairs}",
                                int(pair_count / total_pairs * 50))
                        with model:
                            model.reactions.get_by_id(rid1).knock_out()
                            model.reactions.get_by_id(rid2).knock_out()
                            try:
                                sol = model.optimize()
                                if sol.status == 'optimal':
                                    growth = float(sol.objective_value)
                                    tgt_flux = float(sol.fluxes.get(target_id, 0))
                                    if abs(tgt_flux) < 1e-8 and growth >= baseline_obj * min_growth:
                                        # Only add if not a superset of a smaller cut set
                                        is_superset = any(
                                            set(cs["reactions"]).issubset({rid1, rid2})
                                            for cs in cut_sets if cs["size"] < 2
                                        )
                                        if not is_superset:
                                            cut_sets.append({
                                                "reactions": [rid1, rid2], "size": 2,
                                                "growth": growth,
                                                "growth_pct": growth / baseline_obj * 100,
                                                "target_flux": tgt_flux,
                                            })
                            except Exception:
                                pass

            # Triples
            if max_size >= 3:
                pool3 = search_pool[:min(len(search_pool), 40)]
                triple_count = 0
                for i, rid1 in enumerate(pool3):
                    for j, rid2 in enumerate(pool3[i+1:], i+1):
                        for rid3 in pool3[j+1:]:
                            triple_count += 1
                            if worker and triple_count % 100 == 0:
                                worker.report_progress(f"MCS triples: {triple_count}", 75)
                            with model:
                                model.reactions.get_by_id(rid1).knock_out()
                                model.reactions.get_by_id(rid2).knock_out()
                                model.reactions.get_by_id(rid3).knock_out()
                                try:
                                    sol = model.optimize()
                                    if sol.status == 'optimal':
                                        growth = float(sol.objective_value)
                                        tgt_flux = float(sol.fluxes.get(target_id, 0))
                                        if abs(tgt_flux) < 1e-8 and growth >= baseline_obj * min_growth:
                                            trio = {rid1, rid2, rid3}
                                            is_superset = any(
                                                set(cs["reactions"]).issubset(trio)
                                                for cs in cut_sets if cs["size"] < 3
                                            )
                                            if not is_superset:
                                                cut_sets.append({
                                                    "reactions": [rid1, rid2, rid3],
                                                    "size": 3, "growth": growth,
                                                    "growth_pct": growth / baseline_obj * 100,
                                                    "target_flux": tgt_flux,
                                                })
                                except Exception:
                                    pass

            cut_sets.sort(key=lambda x: (x["size"], -x["growth"]))
            return {"cut_sets": cut_sets[:top_n], "baseline_obj": baseline_obj,
                    "baseline_target": baseline_target, "target_id": target_id}

        def _on_done(result: dict[str, Any]) -> None:
            cs_list = result.get("cut_sets", [])
            if result.get("message"):
                QMessageBox.information(self, "MCS", result["message"])
                return
            dialog2 = QDialog(self)
            dialog2.setWindowTitle(f"Minimal Cut Sets — {result['target_id']}")
            dialog2.resize(800, 500)
            layout2 = QVBoxLayout(dialog2)
            layout2.addWidget(QLabel(
                f"Target: {result['target_id']} (baseline flux: {result['baseline_target']:.6g})\n"
                f"Baseline growth: {result['baseline_obj']:.6g}\n"
                f"Cut sets found: {len(cs_list)}"))
            table = QTableWidget(0, 4)
            table.setHorizontalHeaderLabels(["Reactions to Delete", "Size", "Growth", "Growth %"])
            table.horizontalHeader().setStretchLastSection(True)
            table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            for cs in cs_list:
                r = table.rowCount()
                table.insertRow(r)
                table.setItem(r, 0, QTableWidgetItem(", ".join(cs["reactions"])))
                table.setItem(r, 1, QTableWidgetItem(str(cs["size"])))
                table.setItem(r, 2, QTableWidgetItem(f"{cs['growth']:.6g}"))
                table.setItem(r, 3, QTableWidgetItem(f"{cs['growth_pct']:.1f}%"))
            layout2.addWidget(table)
            btns2 = QDialogButtonBox(QDialogButtonBox.Close)
            btns2.rejected.connect(dialog2.reject)
            layout2.addWidget(btns2)
            dialog2.exec()

        self._launch_worker(_compute, _on_done, "Searching Minimal Cut Sets...")

    # ==================== METABOLIC DISTANCE ====================

    def run_metabolic_distance(self) -> None:
        """Metabolic Distance — compute shortest-path distances between
        metabolites through the reaction network.

        Builds a bipartite graph (metabolites ↔ reactions) and computes
        BFS shortest-path distances between a source metabolite and all
        reachable metabolites.  Useful for identifying metabolic proximity,
        pathway connectivity, and biosynthetic potential.
        """
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Metabolic Distance")
        dialog.resize(500, 250)
        layout = QVBoxLayout(dialog)
        layout.addWidget(QLabel(
            "Compute shortest-path metabolic distances from a source\n"
            "metabolite to all reachable metabolites.\n\n"
            "Distance = minimum number of reaction steps between two\n"
            "metabolites in the stoichiometric network."))

        form = QFormLayout()
        source_edit = QLineEdit()
        source_edit.setPlaceholderText("Source metabolite ID (e.g., pyr_c)")
        if self.base_model:
            completer = QCompleter([m.id for m in self.base_model.metabolites])
            completer.setCaseSensitivity(Qt.CaseInsensitive)
            source_edit.setCompleter(completer)
        form.addRow("Source metabolite:", source_edit)

        exclude_currency_chk = QCheckBox("Exclude currency metabolites (H₂O, H⁺, ATP, ADP, …)")
        exclude_currency_chk.setChecked(True)
        form.addRow("", exclude_currency_chk)

        max_dist_spin = QSpinBox()
        max_dist_spin.setRange(1, 50)
        max_dist_spin.setValue(15)
        form.addRow("Max distance:", max_dist_spin)
        layout.addLayout(form)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)

        if dialog.exec() != QDialog.Accepted:
            return

        source_id = source_edit.text().strip()
        if not source_id:
            QMessageBox.warning(self, "Missing", "Enter a source metabolite ID.")
            return
        if source_id not in self.base_model.metabolites:
            QMessageBox.warning(self, "Not found", f"Metabolite '{source_id}' not in model.")
            return

        exclude_currency = exclude_currency_chk.isChecked()
        max_dist = max_dist_spin.value()

        model = self.base_model

        def _compute() -> dict[str, Any]:
            # Currency metabolites to optionally exclude
            currency_ids = set()
            if exclude_currency:
                currency_patterns = {
                    'h_', 'h2o_', 'atp_', 'adp_', 'amp_', 'nad_', 'nadh_',
                    'nadp_', 'nadph_', 'coa_', 'co2_', 'pi_', 'ppi_', 'o2_',
                    'nh4_', 'fad_', 'fadh2_', 'q8_', 'q8h2_', 'accoa_',
                }
                for met in model.metabolites:
                    for pat in currency_patterns:
                        if met.id.startswith(pat) or met.id == pat.rstrip('_'):
                            currency_ids.add(met.id)
                            break

            # Build adjacency list: metabolite → set of reachable metabolites
            # via one reaction step
            adj: dict[str, set[str]] = {m.id: set() for m in model.metabolites
                                        if m.id not in currency_ids}
            for rxn in model.reactions:
                if rxn.boundary:
                    continue
                mets_in_rxn = [m.id for m in rxn.metabolites if m.id not in currency_ids]
                for m1 in mets_in_rxn:
                    for m2 in mets_in_rxn:
                        if m1 != m2:
                            adj.setdefault(m1, set()).add(m2)

            if source_id in currency_ids:
                raise ValueError(
                    f"Source metabolite '{source_id}' is a currency metabolite "
                    "and was excluded. Uncheck 'Exclude currency metabolites' to include it.")

            # BFS from source
            from collections import deque
            distances: dict[str, int] = {source_id: 0}
            queue: deque[str] = deque([source_id])
            while queue:
                current = queue.popleft()
                current_dist = distances[current]
                if current_dist >= max_dist:
                    continue
                for neighbor in adj.get(current, []):
                    if neighbor not in distances:
                        distances[neighbor] = current_dist + 1
                        queue.append(neighbor)

            # Build result sorted by distance
            dist_list = sorted(
                [{"metabolite": mid, "distance": d,
                  "name": model.metabolites.get_by_id(mid).name if mid in model.metabolites else ""}
                 for mid, d in distances.items() if mid != source_id],
                key=lambda x: x["distance"])

            # Distance distribution
            dist_counts: dict[int, int] = {}
            for entry in dist_list:
                d = entry["distance"]
                dist_counts[d] = dist_counts.get(d, 0) + 1

            return {"source": source_id, "distances": dist_list,
                    "dist_counts": dist_counts, "reachable": len(dist_list),
                    "total_metabolites": len(model.metabolites),
                    "excluded_currency": len(currency_ids),
                    "max_dist": max_dist}

        def _on_done(result: dict[str, Any]) -> None:
            dialog2 = QDialog(self)
            dialog2.setWindowTitle(f"Metabolic Distance from {result['source']}")
            dialog2.resize(900, 700)
            layout2 = QVBoxLayout(dialog2)

            layout2.addWidget(QLabel(
                f"Source: {result['source']} | "
                f"Reachable: {result['reachable']} / {result['total_metabolites']} metabolites\n"
                f"Currency metabolites excluded: {result['excluded_currency']}"))

            # Distance distribution chart
            dist_counts = result["dist_counts"]
            FigureCanvas, Figure = _lazy_mpl_canvas()
            if dist_counts:
                fig = Figure(figsize=(7, 3), dpi=100)
                ax = fig.add_subplot(111)
                dists = sorted(dist_counts.keys())
                counts = [dist_counts[d] for d in dists]
                ax.bar(dists, counts, color='#2196F3', alpha=0.8)
                ax.set_xlabel("Distance (reaction steps)")
                ax.set_ylabel("Number of metabolites")
                ax.set_title(f"Metabolic Distance Distribution from {result['source']}")
                ax.set_xticks(dists)
                fig.tight_layout()
                canvas = FigureCanvas(fig)
                layout2.addWidget(canvas)

            # Table
            table = QTableWidget(0, 3)
            table.setHorizontalHeaderLabels(["Metabolite", "Name", "Distance"])
            table.horizontalHeader().setStretchLastSection(True)
            table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            for entry in result["distances"][:200]:
                r = table.rowCount()
                table.insertRow(r)
                table.setItem(r, 0, QTableWidgetItem(entry["metabolite"]))
                table.setItem(r, 1, QTableWidgetItem(entry["name"]))
                table.setItem(r, 2, QTableWidgetItem(str(entry["distance"])))
            layout2.addWidget(table)

            btns2 = QDialogButtonBox(QDialogButtonBox.Close)
            btns2.rejected.connect(dialog2.reject)
            layout2.addWidget(btns2)
            dialog2.exec()

        self._launch_worker(_compute, _on_done, "Computing Metabolic Distances...")
    # ==================== PARETO FRONTIER ANALYSIS ====================

    def run_pareto_frontier(self) -> None:
        """Multi-objective Pareto Frontier between two reactions.

        Sweeps the first objective across its feasible range while
        maximising the second objective at each step, generating
        the Pareto-optimal trade-off curve.  Results are shown as a
        scatter-line chart together with a tabular export.
        """
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        rxn_ids = sorted(r.id for r in self.base_model.reactions)

        dialog = QDialog(self)
        dialog.setWindowTitle("Pareto Frontier — Multi-Objective Optimisation")
        dialog.setMinimumWidth(480)
        layout = QVBoxLayout(dialog)

        layout.addWidget(QLabel(
            "<b>Compute the Pareto-optimal trade-off curve between two "
            "objective reactions.</b><br>"
            "Objective 1 is swept across its feasible range while "
            "Objective 2 is maximised at each step."
        ))

        form = QFormLayout()

        obj1_edit = QLineEdit()
        obj1_completer = QCompleter(rxn_ids, obj1_edit)
        obj1_completer.setCaseSensitivity(Qt.CaseInsensitive)
        obj1_completer.setFilterMode(Qt.MatchContains)
        obj1_edit.setCompleter(obj1_completer)
        # Default to the model's current objective
        try:
            default_obj = str(self.base_model.objective.expression).split("*")[-1].strip()
            obj1_edit.setText(default_obj)
        except Exception:
            pass
        form.addRow("Objective 1 (swept):", obj1_edit)

        obj2_edit = QLineEdit()
        obj2_completer = QCompleter(rxn_ids, obj2_edit)
        obj2_completer.setCaseSensitivity(Qt.CaseInsensitive)
        obj2_completer.setFilterMode(Qt.MatchContains)
        obj2_edit.setCompleter(obj2_completer)
        form.addRow("Objective 2 (maximised):", obj2_edit)

        steps_spin = QSpinBox()
        steps_spin.setRange(5, 500)
        steps_spin.setValue(20)
        form.addRow("Number of steps:", steps_spin)

        fraction_spin = QDoubleSpinBox()
        fraction_spin.setRange(0.0, 1.0)
        fraction_spin.setValue(1.0)
        fraction_spin.setDecimals(2)
        fraction_spin.setSingleStep(0.05)
        form.addRow("Fraction of Obj1 max:", fraction_spin)

        layout.addLayout(form)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)

        if dialog.exec() != QDialog.Accepted:
            return

        obj1_id = obj1_edit.text().strip()
        obj2_id = obj2_edit.text().strip()
        n_steps = steps_spin.value()
        fraction = fraction_spin.value()

        if not obj1_id or not obj2_id:
            QMessageBox.warning(self, "Missing", "Both objectives are required.")
            return
        if obj1_id == obj2_id:
            QMessageBox.warning(self, "Same objective",
                                "Choose two different reactions.")
            return
        if obj1_id not in self.base_model.reactions:
            QMessageBox.warning(self, "Not found",
                                f"Reaction '{obj1_id}' not in model.")
            return
        if obj2_id not in self.base_model.reactions:
            QMessageBox.warning(self, "Not found",
                                f"Reaction '{obj2_id}' not in model.")
            return

        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self.apply_reaction_overrides_to_model(model)
        self._apply_selected_solver(model)

        def _compute() -> dict[str, Any]:
            import numpy as np

            results: list[dict[str, float]] = []

            # --- Phase 1: Find the max of Objective 1 alone ---
            model.objective = obj1_id
            sol1 = model.optimize()
            if sol1.status != "optimal":
                return {"error": f"Cannot maximise '{obj1_id}': {sol1.status}"}
            obj1_max = sol1.objective_value * fraction

            # --- Phase 2: Find the max of Objective 2 alone ---
            model.objective = obj2_id
            sol2 = model.optimize()
            if sol2.status != "optimal":
                return {"error": f"Cannot maximise '{obj2_id}': {sol2.status}"}
            obj2_max = sol2.objective_value

            # --- Phase 3: Find min of Obj1 when Obj2 is maximised ---
            model.objective = obj1_id
            model.objective_direction = "min"
            sol_min = model.optimize()
            obj1_min = sol_min.objective_value if sol_min.status == "optimal" else 0.0
            model.objective_direction = "max"

            # --- Phase 4: Sweep Obj1 and maximise Obj2 at each step ---
            sweep_values = np.linspace(obj1_min, obj1_max, n_steps)
            rxn1 = model.reactions.get_by_id(obj1_id)

            for val in sweep_values:
                with model:
                    # Fix Obj1 to current sweep value
                    rxn1.lower_bound = float(val)
                    rxn1.upper_bound = float(val)

                    model.objective = obj2_id
                    model.objective_direction = "max"
                    sol = model.optimize()

                    if sol.status == "optimal":
                        results.append({
                            "obj1": float(val),
                            "obj2": float(sol.objective_value),
                        })

            # --- Identify the non-dominated (Pareto-optimal) points ---
            pareto: list[dict[str, float]] = []
            for pt in results:
                dominated = False
                for other in results:
                    if (other["obj1"] >= pt["obj1"] and
                            other["obj2"] >= pt["obj2"] and
                            (other["obj1"] > pt["obj1"] or
                             other["obj2"] > pt["obj2"])):
                        dominated = True
                        break
                if not dominated:
                    pareto.append(pt)

            pareto.sort(key=lambda p: p["obj1"])

            return {
                "obj1_id": obj1_id,
                "obj2_id": obj2_id,
                "obj1_max": float(obj1_max),
                "obj2_max": float(obj2_max),
                "all_points": results,
                "pareto_front": pareto,
                "n_steps": n_steps,
            }

        def _on_done(result: dict[str, Any]) -> None:
            if "error" in result:
                self._show_error("Pareto Frontier failed", result["error"])
                return

            pareto = result["pareto_front"]
            all_pts = result["all_points"]
            obj1_id_r = result["obj1_id"]
            obj2_id_r = result["obj2_id"]

            if not pareto:
                QMessageBox.information(
                    self, "Pareto Frontier",
                    "No feasible points found on the Pareto frontier.")
                return

            dialog2 = QDialog(self)
            dialog2.setWindowTitle("Pareto Frontier Results")
            dialog2.setMinimumSize(800, 600)
            layout2 = QVBoxLayout(dialog2)

            tabs = QTabWidget()

            # --- Tab 1: Chart ---
            chart_tab = QWidget()
            chart_layout = QVBoxLayout(chart_tab)

            FigureCanvas, Figure = _lazy_mpl_canvas()
            fig = Figure(figsize=(8, 5))
            canvas = FigureCanvas(fig)
            ax = fig.add_subplot(111)

            # All feasible points
            if all_pts:
                ax.scatter(
                    [p["obj1"] for p in all_pts],
                    [p["obj2"] for p in all_pts],
                    color="#cccccc", s=20, zorder=1,
                    label="Feasible points",
                )

            # Pareto front
            ax.plot(
                [p["obj1"] for p in pareto],
                [p["obj2"] for p in pareto],
                "o-", color="#d32f2f", markersize=6, linewidth=2,
                zorder=2, label="Pareto front",
            )

            # Utopia point
            ax.scatter(
                [result["obj1_max"]], [result["obj2_max"]],
                marker="*", s=200, color="#1976d2", zorder=3,
                label="Utopia point",
            )

            ax.set_xlabel(obj1_id_r, fontsize=11)
            ax.set_ylabel(obj2_id_r, fontsize=11)
            ax.set_title("Pareto Frontier — Multi-Objective Trade-Off", fontsize=13)
            ax.legend(loc="best")
            ax.grid(True, alpha=0.3)
            fig.tight_layout()

            chart_layout.addWidget(canvas)
            tabs.addTab(chart_tab, "📈 Chart")

            # --- Tab 2: Table ---
            table_tab = QWidget()
            table_layout = QVBoxLayout(table_tab)

            table_layout.addWidget(QLabel(
                f"<b>Pareto-optimal points: {len(pareto)}</b>  |  "
                f"Total feasible: {len(all_pts)}  |  "
                f"Steps: {result['n_steps']}"
            ))

            table = QTableWidget(len(pareto), 3)
            table.setHorizontalHeaderLabels([
                obj1_id_r, obj2_id_r, "Trade-off Slope",
            ])
            table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            table.setAlternatingRowColors(True)
            table.horizontalHeader().setStretchLastSection(True)

            for i, pt in enumerate(pareto):
                table.setItem(i, 0, QTableWidgetItem(f"{pt['obj1']:.6g}"))
                table.setItem(i, 1, QTableWidgetItem(f"{pt['obj2']:.6g}"))
                # Compute local trade-off slope
                if i > 0:
                    d_obj1 = pt["obj1"] - pareto[i - 1]["obj1"]
                    d_obj2 = pt["obj2"] - pareto[i - 1]["obj2"]
                    slope = d_obj2 / d_obj1 if abs(d_obj1) > 1e-12 else float("inf")
                    table.setItem(i, 2, QTableWidgetItem(f"{slope:.4g}"))
                else:
                    table.setItem(i, 2, QTableWidgetItem("—"))

            table_layout.addWidget(table)

            # Export button
            export_btn = QPushButton("💾  Export CSV")

            def _export_csv() -> None:
                path, _ = QFileDialog.getSaveFileName(
                    dialog2, "Export Pareto CSV", "", "CSV Files (*.csv)")
                if not path:
                    return
                try:
                    with open(path, "w", newline="", encoding="utf-8") as fh:
                        writer = csv.writer(fh)
                        writer.writerow([obj1_id_r, obj2_id_r, "pareto_optimal"])
                        pareto_set = {(p["obj1"], p["obj2"]) for p in pareto}
                        for pt in all_pts:
                            is_pareto = (pt["obj1"], pt["obj2"]) in pareto_set
                            writer.writerow([pt["obj1"], pt["obj2"],
                                             "yes" if is_pareto else "no"])
                    QMessageBox.information(
                        dialog2, "Exported",
                        f"Pareto frontier exported to:\n{path}")
                except Exception as exc:
                    QMessageBox.critical(dialog2, "Export Error", str(exc))

            export_btn.clicked.connect(_export_csv)
            table_layout.addWidget(export_btn)
            tabs.addTab(table_tab, "📋 Table")

            # --- Tab 3: Summary ---
            summary_tab = QWidget()
            summary_layout = QVBoxLayout(summary_tab)

            best_obj1 = max(pareto, key=lambda p: p["obj1"])
            best_obj2 = max(pareto, key=lambda p: p["obj2"])

            summary_text = (
                f"<h3>Pareto Frontier Summary</h3>"
                f"<table cellpadding='4'>"
                f"<tr><td><b>Objective 1:</b></td><td>{obj1_id_r}</td></tr>"
                f"<tr><td><b>Objective 2:</b></td><td>{obj2_id_r}</td></tr>"
                f"<tr><td><b>Pareto-optimal points:</b></td><td>{len(pareto)}</td></tr>"
                f"<tr><td><b>Total feasible points:</b></td><td>{len(all_pts)}</td></tr>"
                f"<tr><td><b>Best {obj1_id_r}:</b></td>"
                f"<td>{best_obj1['obj1']:.6g} (at {obj2_id_r} = {best_obj1['obj2']:.6g})</td></tr>"
                f"<tr><td><b>Best {obj2_id_r}:</b></td>"
                f"<td>{best_obj2['obj2']:.6g} (at {obj1_id_r} = {best_obj2['obj1']:.6g})</td></tr>"
                f"<tr><td><b>Utopia point:</b></td>"
                f"<td>({result['obj1_max']:.6g}, {result['obj2_max']:.6g})</td></tr>"
                f"</table>"
            )

            # Knee-point detection (point closest to utopia)
            if len(pareto) >= 2:
                # Normalise to [0,1] range
                o1_vals = [p["obj1"] for p in pareto]
                o2_vals = [p["obj2"] for p in pareto]
                o1_range = max(o1_vals) - min(o1_vals) if max(o1_vals) != min(o1_vals) else 1.0
                o2_range = max(o2_vals) - min(o2_vals) if max(o2_vals) != min(o2_vals) else 1.0
                best_dist = float("inf")
                knee = pareto[0]
                for pt in pareto:
                    d = (((pt["obj1"] - result["obj1_max"]) / o1_range) ** 2 +
                         ((pt["obj2"] - result["obj2_max"]) / o2_range) ** 2) ** 0.5
                    if d < best_dist:
                        best_dist = d
                        knee = pt
                summary_text += (
                    f"<br><b>🎯 Knee point (closest to utopia):</b><br>"
                    f"&nbsp;&nbsp;{obj1_id_r} = {knee['obj1']:.6g}, "
                    f"{obj2_id_r} = {knee['obj2']:.6g}"
                )

            summary_label = QLabel(summary_text)
            summary_label.setTextFormat(Qt.RichText)
            summary_label.setWordWrap(True)
            summary_layout.addWidget(summary_label)
            summary_layout.addStretch()
            tabs.addTab(summary_tab, "📊 Summary")

            layout2.addWidget(tabs)

            btns2 = QDialogButtonBox(QDialogButtonBox.Close)
            btns2.rejected.connect(dialog2.reject)
            layout2.addWidget(btns2)
            dialog2.exec()

        self._launch_worker(_compute, _on_done, "Computing Pareto Frontier...")