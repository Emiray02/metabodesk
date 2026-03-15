"""UI dialogs, statistics, and misc mixin for MetaboDesk.

Contains dialogs and functionality that do not fit neatly into one of
the other specialised mixins:

- **Model Statistics Dashboard**: Comprehensive model metrics.
- **Metabolite Connectivity & Annotation Quality**: Structural analysis.
- **Sensitivity Analysis**: Single-parameter sweep of a reaction bound.
- **Batch Analysis**: FBA/FVA/SGD across multiple SBML files.
- **Community Analysis**: Multi-species FBA with shared extracellular pool.
- **Reaction / Metabolite detail dialogs**: Inline model inspection.
- **Undo / Redo**: State-snapshot based undo stack.
- **Keyboard shortcuts, drag-and-drop, theme toggle, plugin loader**.
"""

import csv
import logging
import importlib
import importlib.util

import cobra
from datetime import datetime
from pathlib import Path

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QFileDialog, QMessageBox, QTableWidget, QTableWidgetItem,
    QAbstractItemView, QHeaderView, QLineEdit, QSpinBox, QTabWidget,
    QDoubleSpinBox, QCheckBox, QComboBox, QPlainTextEdit, QDialog,
    QDialogButtonBox, QCompleter, QGroupBox, QFormLayout, QScrollArea,
    QProgressDialog, QApplication,
)

from metabodesk_core.constants import DARK_STYLESHEET
from metabodesk_core.utils import get_kegg_rxn_id, get_ec_numbers, get_gpr

logger = logging.getLogger("MetaboDesk")

class DialogsMixin:
    """Mixin providing dialogs functionality."""

    def show_model_stats_dashboard(self):
        """Comprehensive model statistics dialog.

        Heavy computations (blocked-reaction detection, mass-balance checks)
        run in a background thread so the UI stays responsive.
        """
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        model = self.base_model

        # Show a progress dialog immediately so user knows work is happening
        progress = QProgressDialog("Computing model statistics…", "Cancel", 0, 0, self)
        progress.setWindowTitle("Model Statistics Dashboard")
        progress.setWindowModality(Qt.WindowModal)
        progress.setMinimumDuration(0)
        progress.setValue(0)
        progress.show()
        QApplication.processEvents()

        # ---- Heavy computation in a worker thread ----
        from threading import Thread
        _result_holder: dict = {}

        def _compute():
            n_rxn = len(model.reactions)
            n_met = len(model.metabolites)
            n_gene = len(model.genes)

            compartments = set()
            for m in model.metabolites:
                c = getattr(m, "compartment", None)
                if c:
                    compartments.add(c)
            compartment_names = getattr(model, "compartments", {}) or {}
            comp_lines = []
            for c in sorted(compartments):
                name = compartment_names.get(c, "")
                mets_in = sum(1 for m in model.metabolites if getattr(m, "compartment", "") == c)
                comp_lines.append(f"  {c} ({name}): {mets_in} metabolites" if name else f"  {c}: {mets_in} metabolites")

            subsystems = set()
            for r in model.reactions:
                ss = getattr(r, "subsystem", None) or ""
                if ss.strip():
                    subsystems.add(ss.strip())

            n_reversible = sum(1 for r in model.reactions if r.reversibility)
            pct_rev = (n_reversible / n_rxn * 100) if n_rxn else 0
            exchange_rxns = [r for r in model.reactions if r.boundary]
            n_exchange = len(exchange_rxns)

            blocked: list = []
            try:
                from cobra.flux_analysis import find_blocked_reactions
                blocked = find_blocked_reactions(model.copy(), processes=1)
            except Exception:
                pass
            n_blocked = len(blocked)

            orphan_mets = [m.id for m in model.metabolites if len(m.reactions) <= 1]
            n_orphan = len(orphan_mets)

            dead_ends: list[str] = []
            for m in model.metabolites:
                producers = sum(1 for r in m.reactions if r.metabolites[m] > 0)
                consumers = sum(1 for r in m.reactions if r.metabolites[m] < 0)
                if producers == 0 or consumers == 0:
                    dead_ends.append(m.id)
            n_dead_end = len(dead_ends)

            n_imbalanced = 0
            for rxn in model.reactions:
                if rxn.boundary:
                    continue
                try:
                    bal = rxn.check_mass_balance()
                    if bal:
                        n_imbalanced += 1
                except Exception:
                    pass

            n_gpr = sum(1 for r in model.reactions if getattr(r, "gene_reaction_rule", "").strip())
            pct_gpr = (n_gpr / n_rxn * 100) if n_rxn else 0

            _result_holder.update({
                "n_rxn": n_rxn, "n_met": n_met, "n_gene": n_gene,
                "compartments": compartments, "comp_lines": comp_lines,
                "subsystems": subsystems, "n_reversible": n_reversible,
                "pct_rev": pct_rev, "n_exchange": n_exchange,
                "blocked": blocked, "n_blocked": n_blocked,
                "orphan_mets": orphan_mets, "n_orphan": n_orphan,
                "dead_ends": dead_ends, "n_dead_end": n_dead_end,
                "n_imbalanced": n_imbalanced, "n_gpr": n_gpr,
                "pct_gpr": pct_gpr,
            })

        worker = Thread(target=_compute, daemon=True)
        worker.start()

        # Keep UI responsive while waiting
        while worker.is_alive():
            QApplication.processEvents()
            worker.join(timeout=0.05)

        progress.close()

        if not _result_holder:
            QMessageBox.warning(self, "Error", "Statistics computation failed.")
            return

        R = _result_holder  # shorthand
        model_id = getattr(model, "id", "N/A") or "N/A"
        model_name = getattr(model, "name", "") or ""

        # Build dialog
        dialog = QDialog(self)
        dialog.setWindowTitle("Model Statistics Dashboard")
        dialog.resize(650, 600)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content = QWidget()
        layout = QVBoxLayout(content)

        html = f"""
        <h2 style='color:#2E86C1;'>Model Statistics Dashboard</h2>
        <table style='font-size:11px; border-collapse:collapse; width:100%;'>
        <tr><td style='padding:4px;'><b>Model ID:</b></td><td>{model_id}</td></tr>
        <tr><td style='padding:4px;'><b>Model Name:</b></td><td>{model_name}</td></tr>
        <tr><td colspan='2'><hr></td></tr>
        <tr><td style='padding:4px;'><b>Reactions:</b></td><td>{R['n_rxn']}</td></tr>
        <tr><td style='padding:4px;'><b>Metabolites:</b></td><td>{R['n_met']}</td></tr>
        <tr><td style='padding:4px;'><b>Genes:</b></td><td>{R['n_gene']}</td></tr>
        <tr><td style='padding:4px;'><b>Compartments:</b></td><td>{len(R['compartments'])}</td></tr>
        <tr><td style='padding:4px;'><b>Subsystems:</b></td><td>{len(R['subsystems'])}</td></tr>
        <tr><td colspan='2'><hr></td></tr>
        <tr><td style='padding:4px;'><b>Reversible reactions:</b></td><td>{R['n_reversible']} ({R['pct_rev']:.1f}%)</td></tr>
        <tr><td style='padding:4px;'><b>Exchange/boundary reactions:</b></td><td>{R['n_exchange']}</td></tr>
        <tr><td style='padding:4px;'><b>Blocked reactions:</b></td><td>{R['n_blocked']}</td></tr>
        <tr><td style='padding:4px;'><b>Mass-unbalanced reactions:</b></td><td>{R['n_imbalanced']}</td></tr>
        <tr><td colspan='2'><hr></td></tr>
        <tr><td style='padding:4px;'><b>Orphan metabolites (≤1 rxn):</b></td><td>{R['n_orphan']}</td></tr>
        <tr><td style='padding:4px;'><b>Dead-end metabolites:</b></td><td>{R['n_dead_end']}</td></tr>
        <tr><td colspan='2'><hr></td></tr>
        <tr><td style='padding:4px;'><b>GPR coverage:</b></td><td>{R['n_gpr']}/{R['n_rxn']} ({R['pct_gpr']:.1f}%)</td></tr>
        </table>
        """

        info_label = QLabel(html)
        info_label.setWordWrap(True)
        info_label.setTextFormat(Qt.RichText)
        layout.addWidget(info_label)

        # Compartments detail
        if R['comp_lines']:
            comp_group = QGroupBox("Compartments")
            comp_layout = QVBoxLayout()
            comp_text = QPlainTextEdit("\n".join(R['comp_lines']))
            comp_text.setReadOnly(True)
            comp_text.setMaximumHeight(100)
            comp_layout.addWidget(comp_text)
            comp_group.setLayout(comp_layout)
            layout.addWidget(comp_group)

        # Blocked reactions detail
        if R['blocked']:
            blocked_group = QGroupBox(f"Blocked Reactions ({R['n_blocked']})")
            blocked_layout = QVBoxLayout()
            blocked_text = QPlainTextEdit("\n".join(str(b) for b in R['blocked'][:200]))
            blocked_text.setReadOnly(True)
            blocked_text.setMaximumHeight(120)
            blocked_layout.addWidget(blocked_text)
            blocked_group.setLayout(blocked_layout)
            layout.addWidget(blocked_group)

        # Orphan metabolites detail
        if R['orphan_mets']:
            orphan_group = QGroupBox(f"Orphan Metabolites ({R['n_orphan']})")
            orphan_layout = QVBoxLayout()
            orphan_text = QPlainTextEdit("\n".join(R['orphan_mets'][:200]))
            orphan_text.setReadOnly(True)
            orphan_text.setMaximumHeight(120)
            orphan_layout.addWidget(orphan_text)
            orphan_group.setLayout(orphan_layout)
            layout.addWidget(orphan_group)

        layout.addStretch()
        scroll.setWidget(content)

        dlg_layout = QVBoxLayout(dialog)
        dlg_layout.addWidget(scroll)

        btn_box = QDialogButtonBox(QDialogButtonBox.Close)
        btn_box.rejected.connect(dialog.reject)
        dlg_layout.addWidget(btn_box)

        dialog.exec()

    # ---- Metabolite Connectivity Analysis ----

    def show_metabolite_connectivity(self):
        """Analyse dead-end metabolites and chokepoint reactions."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        model = self.base_model

        # Dead-end metabolites: only produced or only consumed
        dead_ends = []
        for m in model.metabolites:
            producers = 0
            consumers = 0
            for r in m.reactions:
                coeff = r.metabolites[m]
                if coeff > 0:
                    producers += 1
                elif coeff < 0:
                    consumers += 1
                # For reversible reactions, both directions possible
                if r.reversibility:
                    producers += 1
                    consumers += 1
            if producers == 0 or consumers == 0:
                status = "only consumed" if producers == 0 else "only produced"
                dead_ends.append((m.id, m.name or "", getattr(m, "compartment", ""), status))

        # Chokepoint reactions: the only reaction producing or consuming a metabolite
        chokepoints = set()
        for m in model.metabolites:
            producing_rxns = [r for r in m.reactions if r.metabolites[m] > 0]
            consuming_rxns = [r for r in m.reactions if r.metabolites[m] < 0]
            if len(producing_rxns) == 1:
                chokepoints.add(producing_rxns[0].id)
            if len(consuming_rxns) == 1:
                chokepoints.add(consuming_rxns[0].id)

        # Build dialog
        dialog = QDialog(self)
        dialog.setWindowTitle("Metabolite Connectivity Analysis")
        dialog.resize(800, 600)

        tabs = QTabWidget()

        # Tab 1: Dead-end metabolites
        dead_end_tab = QWidget()
        de_layout = QVBoxLayout(dead_end_tab)
        de_layout.addWidget(QLabel(f"<b>Dead-end metabolites: {len(dead_ends)}</b>  "
                                   "(metabolites only produced or only consumed, excluding reversible contributions)"))

        de_filter = QLineEdit()
        de_filter.setPlaceholderText("Filter dead-ends...")
        de_layout.addWidget(de_filter)

        de_table = QTableWidget()
        de_table.setColumnCount(4)
        de_table.setHorizontalHeaderLabels(["ID", "Name", "Compartment", "Status"])
        de_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        de_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        de_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        de_table.setAlternatingRowColors(True)
        de_table.setRowCount(len(dead_ends))

        for i, (mid, mname, mcomp, mstatus) in enumerate(dead_ends):
            de_table.setItem(i, 0, QTableWidgetItem(mid))
            de_table.setItem(i, 1, QTableWidgetItem(mname))
            de_table.setItem(i, 2, QTableWidgetItem(mcomp))
            de_table.setItem(i, 3, QTableWidgetItem(mstatus))

        def de_filter_changed(text):
            t = text.strip().lower()
            for row in range(de_table.rowCount()):
                match = not t
                for col in range(de_table.columnCount()):
                    item = de_table.item(row, col)
                    if item and t in item.text().lower():
                        match = True
                        break
                de_table.setRowHidden(row, not match)

        de_filter.textChanged.connect(de_filter_changed)
        de_layout.addWidget(de_table)
        tabs.addTab(dead_end_tab, f"Dead-end Metabolites ({len(dead_ends)})")

        # Tab 2: Chokepoint reactions
        choke_tab = QWidget()
        ch_layout = QVBoxLayout(choke_tab)
        ch_layout.addWidget(QLabel(f"<b>Chokepoint reactions: {len(chokepoints)}</b>  "
                                   "(only reaction producing/consuming a metabolite — potential drug targets)"))

        ch_filter = QLineEdit()
        ch_filter.setPlaceholderText("Filter chokepoints...")
        ch_layout.addWidget(ch_filter)

        choke_list = sorted(chokepoints)
        ch_table = QTableWidget()
        ch_table.setColumnCount(3)
        ch_table.setHorizontalHeaderLabels(["Reaction ID", "Name", "Subsystem"])
        ch_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        ch_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        ch_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        ch_table.setAlternatingRowColors(True)
        ch_table.setRowCount(len(choke_list))

        for i, rid in enumerate(choke_list):
            rxn = model.reactions.get_by_id(rid)
            ch_table.setItem(i, 0, QTableWidgetItem(rid))
            ch_table.setItem(i, 1, QTableWidgetItem(rxn.name or ""))
            ch_table.setItem(i, 2, QTableWidgetItem(getattr(rxn, "subsystem", "") or ""))

        def ch_filter_changed(text):
            t = text.strip().lower()
            for row in range(ch_table.rowCount()):
                match = not t
                for col in range(ch_table.columnCount()):
                    item = ch_table.item(row, col)
                    if item and t in item.text().lower():
                        match = True
                        break
                ch_table.setRowHidden(row, not match)

        ch_filter.textChanged.connect(ch_filter_changed)
        ch_layout.addWidget(ch_table)
        tabs.addTab(choke_tab, f"Chokepoint Reactions ({len(chokepoints)})")

        dlg_layout = QVBoxLayout(dialog)
        dlg_layout.addWidget(tabs)
        btn_box = QDialogButtonBox(QDialogButtonBox.Close)
        btn_box.rejected.connect(dialog.reject)
        dlg_layout.addWidget(btn_box)

        dialog.exec()

    # ---- Annotation Quality Score ----

    def show_annotation_quality(self):
        """Show annotation coverage (KEGG, EC, GPR, SBO, MetaNetX)."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        model = self.base_model
        n_rxn = len(model.reactions)
        n_met = len(model.metabolites)

        # Reaction annotations
        n_kegg = sum(1 for r in model.reactions if get_kegg_rxn_id(r))
        n_ec = sum(1 for r in model.reactions if get_ec_numbers(r))
        n_gpr = sum(1 for r in model.reactions if get_gpr(r).strip())

        # SBO terms
        n_rxn_sbo = 0
        n_met_sbo = 0
        for r in model.reactions:
            ann = getattr(r, "annotation", {}) or {}
            if ann.get("sbo") or ann.get("sbo_term"):
                n_rxn_sbo += 1
        for m in model.metabolites:
            ann = getattr(m, "annotation", {}) or {}
            if ann.get("sbo") or ann.get("sbo_term"):
                n_met_sbo += 1

        # Metabolite annotations: KEGG compound, ChEBI, MetaNetX, InChI, formula
        n_met_kegg = 0
        n_met_chebi = 0
        n_met_metanetx = 0
        n_met_inchi = 0
        n_met_formula = 0
        for m in model.metabolites:
            ann = getattr(m, "annotation", {}) or {}
            if any(ann.get(k) for k in ("kegg.compound", "kegg_compound", "kegg.metabolite")):
                n_met_kegg += 1
            if any(ann.get(k) for k in ("chebi", "CHEBI")):
                n_met_chebi += 1
            if any(ann.get(k) for k in ("metanetx.chemical", "metanetx", "MetaNetX")):
                n_met_metanetx += 1
            if any(ann.get(k) for k in ("inchi", "InChI", "inchikey", "InChIKey")):
                n_met_inchi += 1
            formula = getattr(m, "formula", None) or ""
            if formula.strip():
                n_met_formula += 1

        def pct(n, total):
            return f"{n}/{total} ({n / total * 100:.1f}%)" if total else "0/0 (0%)"

        def bar_color(n, total):
            if total == 0:
                return "#999"
            ratio = n / total
            if ratio >= 0.8:
                return "#27AE60"
            elif ratio >= 0.5:
                return "#F39C12"
            else:
                return "#E74C3C"

        def bar_html(label, n, total):
            ratio = (n / total * 100) if total else 0
            color = bar_color(n, total)
            return (f"<tr><td style='padding:3px;'><b>{label}</b></td>"
                    f"<td style='padding:3px;'>{pct(n, total)}</td>"
                    f"<td style='padding:3px; width:200px;'>"
                    f"<div style='background:#eee; border-radius:4px; height:16px;'>"
                    f"<div style='background:{color}; border-radius:4px; height:16px; width:{ratio:.0f}%;'></div>"
                    f"</div></td></tr>")

        html = f"""
        <h2 style='color:#2E86C1;'>Annotation Quality Score</h2>
        <h3>Reaction Annotations ({n_rxn} reactions)</h3>
        <table style='font-size:11px; width:100%;'>
        {bar_html("KEGG Reaction", n_kegg, n_rxn)}
        {bar_html("EC Number", n_ec, n_rxn)}
        {bar_html("GPR (gene rules)", n_gpr, n_rxn)}
        {bar_html("SBO Term", n_rxn_sbo, n_rxn)}
        </table>
        <br>
        <h3>Metabolite Annotations ({n_met} metabolites)</h3>
        <table style='font-size:11px; width:100%;'>
        {bar_html("KEGG Compound", n_met_kegg, n_met)}
        {bar_html("ChEBI", n_met_chebi, n_met)}
        {bar_html("MetaNetX", n_met_metanetx, n_met)}
        {bar_html("InChI / InChIKey", n_met_inchi, n_met)}
        {bar_html("Chemical Formula", n_met_formula, n_met)}
        {bar_html("SBO Term", n_met_sbo, n_met)}
        </table>
        """

        # Overall score (weighted average)
        total_checks = n_rxn * 4 + n_met * 6
        total_hits = n_kegg + n_ec + n_gpr + n_rxn_sbo + n_met_kegg + n_met_chebi + n_met_metanetx + n_met_inchi + n_met_formula + n_met_sbo
        overall = (total_hits / total_checks * 100) if total_checks else 0
        ov_color = bar_color(total_hits, total_checks)

        html += f"""
        <br>
        <h3>Overall Annotation Score</h3>
        <div style='text-align:center; font-size:28px; font-weight:bold; color:{ov_color};'>{overall:.1f}%</div>
        <div style='background:#eee; border-radius:6px; height:24px; margin:8px 0;'>
            <div style='background:{ov_color}; border-radius:6px; height:24px; width:{overall:.0f}%;'></div>
        </div>
        """

        dialog = QDialog(self)
        dialog.setWindowTitle("Annotation Quality Score")
        dialog.resize(650, 550)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content = QWidget()
        cl = QVBoxLayout(content)
        lbl = QLabel(html)
        lbl.setWordWrap(True)
        lbl.setTextFormat(Qt.RichText)
        cl.addWidget(lbl)
        cl.addStretch()
        scroll.setWidget(content)

        dlg_layout = QVBoxLayout(dialog)
        dlg_layout.addWidget(scroll)
        btn_box = QDialogButtonBox(QDialogButtonBox.Close)
        btn_box.rejected.connect(dialog.reject)
        dlg_layout.addWidget(btn_box)

        dialog.exec()

    # ---- Add Custom Constraint ----

    def add_custom_constraint(self):
        """Dialog to add an arbitrary linear constraint to the model."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Add Custom Constraint")
        dialog.resize(600, 400)

        layout = QVBoxLayout(dialog)

        info = QLabel(
            "<b>Add a linear constraint to the model.</b><br>"
            "Format: <code>coeff1 * reaction_id1 + coeff2 * reaction_id2 ... [&lt;=, =, &gt;=] bound</code><br>"
            "Example: <code>1.0 * PFK + -0.5 * PYK &lt;= 10.0</code><br><br>"
            "Or simply set reaction bounds below."
        )
        info.setWordWrap(True)
        info.setTextFormat(Qt.RichText)
        layout.addWidget(info)

        # Constraint name
        name_layout = QHBoxLayout()
        name_layout.addWidget(QLabel("Constraint name:"))
        name_edit = QLineEdit()
        name_edit.setPlaceholderText("my_constraint_1")
        name_layout.addWidget(name_edit)
        layout.addLayout(name_layout)

        # Expression mode: simple (reaction bounds) or advanced (linear expression)
        mode_tabs = QTabWidget()

        # Simple mode: pick reaction, set bounds
        simple_tab = QWidget()
        simple_layout = QFormLayout(simple_tab)

        rxn_combo = QComboBox()
        rxn_combo.setEditable(True)
        rxn_ids = sorted(r.id for r in self.base_model.reactions)
        rxn_combo.addItems(rxn_ids)
        rxn_combo.setCompleter(QCompleter(rxn_ids))
        simple_layout.addRow("Reaction:", rxn_combo)

        lb_spin = QDoubleSpinBox()
        lb_spin.setRange(-1e9, 1e9)
        lb_spin.setDecimals(6)
        lb_spin.setValue(0.0)
        simple_layout.addRow("Lower bound:", lb_spin)

        ub_spin = QDoubleSpinBox()
        ub_spin.setRange(-1e9, 1e9)
        ub_spin.setDecimals(6)
        ub_spin.setValue(1000.0)
        simple_layout.addRow("Upper bound:", ub_spin)

        mode_tabs.addTab(simple_tab, "Reaction Bounds")

        # Advanced mode: linear expression
        adv_tab = QWidget()
        adv_layout = QVBoxLayout(adv_tab)
        adv_layout.addWidget(QLabel("Enter linear constraint expression:"))

        expr_edit = QPlainTextEdit()
        expr_edit.setPlaceholderText("1.0 * PFK + -0.5 * PYK <= 10.0")
        expr_edit.setMaximumHeight(80)
        adv_layout.addWidget(expr_edit)

        adv_layout.addWidget(QLabel(
            "<small>Supported operators: <=, >=, =<br>"
            "Use reaction IDs from your model. Coefficients default to 1.0 if omitted.</small>"
        ))
        adv_layout.addStretch()
        mode_tabs.addTab(adv_tab, "Linear Expression")

        layout.addWidget(mode_tabs)

        btn_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btn_box.accepted.connect(dialog.accept)
        btn_box.rejected.connect(dialog.reject)
        layout.addWidget(btn_box)

        if dialog.exec() != QDialog.Accepted:
            return

        try:
            if mode_tabs.currentIndex() == 0:
                # Simple mode: set reaction bounds
                rid = rxn_combo.currentText().strip()
                if rid not in [r.id for r in self.base_model.reactions]:
                    QMessageBox.warning(self, "Invalid reaction", f"Reaction '{rid}' not found in model.")
                    return
                lb = lb_spin.value()
                ub = ub_spin.value()
                if lb > ub:
                    QMessageBox.warning(self, "Invalid bounds", "Lower bound must be ≤ upper bound.")
                    return
                # Store in reaction_bound_overrides
                self.reaction_bound_overrides[rid] = (lb, ub)
                self.model_dirty = True
                QMessageBox.information(self, "Constraint added",
                                        f"Set {rid} bounds to [{lb:.6g}, {ub:.6g}].\n"
                                        f"This will be applied in all subsequent analyses.")
            else:
                # Advanced mode: parse linear expression
                expr_text = expr_edit.toPlainText().strip()
                if not expr_text:
                    QMessageBox.warning(self, "Empty expression", "Please enter a constraint expression.")
                    return
                self._apply_linear_constraint(expr_text, name_edit.text().strip())
        except Exception as e:
            self._show_error("Constraint error", "Failed to apply constraint.", e)

    def _apply_linear_constraint(self, expr_text: str, name: str):
        """Parse and apply a linear constraint like '1.0 * PFK + -0.5 * PYK <= 10.0'."""
        import re as _re

        # Determine operator and split
        if '<=' in expr_text:
            parts = expr_text.split('<=', 1)
            sense = 'L'
        elif '>=' in expr_text:
            parts = expr_text.split('>=', 1)
            sense = 'G'
        elif '=' in expr_text:
            parts = expr_text.split('=', 1)
            sense = 'E'
        else:
            raise ValueError("No comparison operator found. Use <=, >=, or =")

        lhs_str = parts[0].strip()
        rhs_val = float(parts[1].strip())

        # Parse LHS: terms like "1.0 * PFK" or "PFK" or "-0.5*PYK"
        terms = {}
        # Normalize: add spaces around + and -
        lhs_str = _re.sub(r'(?<=[A-Za-z0-9_])\s*\+\s*', ' + ', lhs_str)
        lhs_str = _re.sub(r'(?<=[A-Za-z0-9_])\s*-\s*', ' + -', lhs_str)

        # Split on +
        for token in lhs_str.split('+'):
            token = token.strip()
            if not token:
                continue
            if '*' in token:
                coeff_str, rxn_id = token.split('*', 1)
                coeff = float(coeff_str.strip())
                rxn_id = rxn_id.strip()
            else:
                rxn_id = token.strip()
                coeff = 1.0
                if rxn_id.startswith('-'):
                    rxn_id = rxn_id[1:].strip()
                    coeff = -1.0

            if rxn_id not in [r.id for r in self.base_model.reactions]:
                raise ValueError(f"Reaction '{rxn_id}' not found in model.")
            terms[rxn_id] = terms.get(rxn_id, 0) + coeff

        if not terms:
            raise ValueError("No valid terms found in expression.")

        # Apply constraint using COBRApy
        model = self.base_model
        constraint_name = name or f"custom_{len(model.constraints) + 1}"

        constraint = model.problem.Constraint(
            sum(terms[rid] * model.reactions.get_by_id(rid).flux_expression for rid in terms),
            lb=rhs_val if sense in ('G', 'E') else None,
            ub=rhs_val if sense in ('L', 'E') else None,
            name=constraint_name
        )
        model.add_cons_vars(constraint)
        self.model_dirty = True

        terms_str = " + ".join(f"{c:g}*{r}" for r, c in terms.items())
        op = {"L": "<=", "G": ">=", "E": "="}[sense]
        QMessageBox.information(self, "Constraint added",
                                f"Added constraint '{constraint_name}':\n"
                                f"{terms_str} {op} {rhs_val}\n\n"
                                f"Applied directly to the loaded model.")

    # ---- Generate HTML Report ----

    def generate_html_report(self):
        """Generate a comprehensive HTML report with model info and analysis results."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        path, _ = QFileDialog.getSaveFileName(
            self, "Save HTML Report", "", "HTML files (*.html)"
        )
        if not path:
            return

        try:
            model = self.base_model
            model_id = getattr(model, "id", "N/A") or "N/A"
            n_rxn = len(model.reactions)
            n_met = len(model.metabolites)
            n_gene = len(model.genes)

            # Compartments
            compartments = {}
            for m in model.metabolites:
                c = getattr(m, "compartment", "?")
                compartments[c] = compartments.get(c, 0) + 1

            # Subsystems
            subsystems = {}
            for r in model.reactions:
                ss = getattr(r, "subsystem", "") or "Unassigned"
                subsystems[ss] = subsystems.get(ss, 0) + 1

            # Annotation coverage
            n_kegg = sum(1 for r in model.reactions if get_kegg_rxn_id(r))
            n_ec = sum(1 for r in model.reactions if get_ec_numbers(r))
            n_gpr = sum(1 for r in model.reactions if get_gpr(r).strip())

            # Exchange reactions
            exchange_rxns = [r for r in model.reactions if r.boundary]

            # Build HTML
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            html_parts = [f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>MetaboDesk Report — {model_id}</title>
<style>
  body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 20px 40px; color: #333; line-height: 1.6; }}
  h1 {{ color: #2C3E50; border-bottom: 3px solid #2E86C1; padding-bottom: 8px; }}
  h2 {{ color: #2E86C1; margin-top: 30px; }}
  h3 {{ color: #555; }}
  table {{ border-collapse: collapse; width: 100%; margin: 10px 0 20px 0; }}
  th {{ background: #2E86C1; color: white; padding: 8px 12px; text-align: left; }}
  td {{ padding: 6px 12px; border-bottom: 1px solid #ddd; }}
  tr:nth-child(even) {{ background: #f8f9fa; }}
  .metric {{ display: inline-block; background: #EBF5FB; border-radius: 8px; padding: 15px 25px; margin: 5px; text-align: center; min-width: 120px; }}
  .metric .value {{ font-size: 28px; font-weight: bold; color: #2E86C1; }}
  .metric .label {{ font-size: 12px; color: #666; }}
  .bar {{ background: #eee; border-radius: 4px; height: 16px; }}
  .bar-fill {{ border-radius: 4px; height: 16px; }}
  .green {{ background: #27AE60; }}
  .yellow {{ background: #F39C12; }}
  .red {{ background: #E74C3C; }}
  .footer {{ margin-top: 40px; padding-top: 10px; border-top: 1px solid #ddd; color: #999; font-size: 11px; }}
</style>
</head>
<body>
<h1>MetaboDesk — Model Report</h1>
<p><b>Generated:</b> {timestamp}</p>

<h2>Model Overview</h2>
<div>
  <div class="metric"><div class="value">{model_id}</div><div class="label">Model ID</div></div>
</div>
<div>
  <div class="metric"><div class="value">{n_rxn}</div><div class="label">Reactions</div></div>
  <div class="metric"><div class="value">{n_met}</div><div class="label">Metabolites</div></div>
  <div class="metric"><div class="value">{n_gene}</div><div class="label">Genes</div></div>
  <div class="metric"><div class="value">{len(compartments)}</div><div class="label">Compartments</div></div>
  <div class="metric"><div class="value">{len(exchange_rxns)}</div><div class="label">Exchange Rxns</div></div>
</div>
"""]

            # Compartments table
            html_parts.append("<h2>Compartments</h2><table><tr><th>ID</th><th>Metabolites</th></tr>")
            comp_names = getattr(model, "compartments", {}) or {}
            for c in sorted(compartments.keys()):
                cname = comp_names.get(c, "")
                label = f"{c} ({cname})" if cname else c
                html_parts.append(f"<tr><td>{label}</td><td>{compartments[c]}</td></tr>")
            html_parts.append("</table>")

            # Annotation coverage
            def pct_str(n, total):
                return f"{n}/{total} ({n / total * 100:.1f}%)" if total else "0"

            def bar_cls(n, total):
                if total == 0:
                    return "red"
                r = n / total
                return "green" if r >= 0.8 else ("yellow" if r >= 0.5 else "red")

            def bar_row(label, n, total):
                w = (n / total * 100) if total else 0
                cls = bar_cls(n, total)
                return (f"<tr><td>{label}</td><td>{pct_str(n, total)}</td>"
                        f"<td><div class='bar'><div class='bar-fill {cls}' style='width:{w:.0f}%'></div></div></td></tr>")

            html_parts.append("<h2>Annotation Coverage</h2><table><tr><th>Annotation</th><th>Coverage</th><th>Bar</th></tr>")
            html_parts.append(bar_row("KEGG Reaction", n_kegg, n_rxn))
            html_parts.append(bar_row("EC Number", n_ec, n_rxn))
            html_parts.append(bar_row("GPR Rules", n_gpr, n_rxn))
            html_parts.append("</table>")

            # Subsystem distribution (top 20)
            html_parts.append("<h2>Subsystem Distribution</h2><table><tr><th>Subsystem</th><th>Reactions</th></tr>")
            for ss, cnt in sorted(subsystems.items(), key=lambda x: -x[1])[:30]:
                html_parts.append(f"<tr><td>{ss}</td><td>{cnt}</td></tr>")
            if len(subsystems) > 30:
                html_parts.append(f"<tr><td><i>... and {len(subsystems) - 30} more</i></td><td></td></tr>")
            html_parts.append("</table>")

            # Last analysis results
            if self.last_run:
                html_parts.append("<h2>Last Analysis Results</h2>")
                analysis_type = self.last_run.get("type", "Unknown")
                html_parts.append(f"<p><b>Analysis type:</b> {analysis_type}</p>")

                if "baseline" in self.last_run:
                    bl = self.last_run["baseline"]
                    obj_val = bl.get("objective_value")
                    if obj_val is not None:
                        html_parts.append(f"<p><b>Objective value:</b> {obj_val:.6g}</p>")
                    fluxes = bl.get("fluxes")
                    if fluxes and isinstance(fluxes, dict):
                        html_parts.append("<h3>Flux Distribution (non-zero)</h3>")
                        html_parts.append("<table><tr><th>Reaction</th><th>Flux</th></tr>")
                        nonzero = {k: v for k, v in fluxes.items() if abs(v) > 1e-9}
                        for rid in sorted(nonzero, key=lambda x: -abs(nonzero[x]))[:50]:
                            html_parts.append(f"<tr><td>{rid}</td><td>{nonzero[rid]:.6g}</td></tr>")
                        if len(nonzero) > 50:
                            html_parts.append(f"<tr><td><i>... and {len(nonzero) - 50} more</i></td><td></td></tr>")
                        html_parts.append("</table>")

                if "fva" in self.last_run:
                    fva_df = self.last_run["fva"]
                    if hasattr(fva_df, "shape"):
                        html_parts.append(f"<h3>FVA Results ({fva_df.shape[0]} reactions)</h3>")
                        html_parts.append("<table><tr><th>Reaction</th><th>Minimum</th><th>Maximum</th></tr>")
                        count = 0
                        for idx, row in fva_df.iterrows():
                            if count >= 50:
                                html_parts.append(f"<tr><td><i>... {fva_df.shape[0] - 50} more rows</i></td><td></td><td></td></tr>")
                                break
                            html_parts.append(f"<tr><td>{idx}</td><td>{row.get('minimum', 0):.6g}</td><td>{row.get('maximum', 0):.6g}</td></tr>")
                            count += 1
                        html_parts.append("</table>")

            # Exchange reactions table
            html_parts.append("<h2>Exchange Reactions</h2><table><tr><th>ID</th><th>Name</th><th>Lower Bound</th><th>Upper Bound</th></tr>")
            for r in sorted(exchange_rxns, key=lambda x: x.id)[:50]:
                html_parts.append(f"<tr><td>{r.id}</td><td>{r.name or ''}</td><td>{r.lower_bound:.6g}</td><td>{r.upper_bound:.6g}</td></tr>")
            if len(exchange_rxns) > 50:
                html_parts.append(f"<tr><td><i>... and {len(exchange_rxns) - 50} more</i></td><td></td><td></td><td></td></tr>")
            html_parts.append("</table>")

            # Footer
            html_parts.append(f"""
<div class="footer">
  Generated by MetaboDesk v1.0.0 | {timestamp}
</div>
</body>
</html>""")

            html_content = "\n".join(html_parts)

            with open(path, "w", encoding="utf-8") as f:
                f.write(html_content)

            QMessageBox.information(self, "Report saved",
                                    f"HTML report saved to:\n{path}")

            # Offer to open in browser
            reply = QMessageBox.question(self, "Open report?",
                                         "Open the report in your web browser?",
                                         QMessageBox.Yes | QMessageBox.No)
            if reply == QMessageBox.Yes:
                QDesktopServices.openUrl(QUrl.fromLocalFile(path))

        except Exception as e:
            self._show_error("Report error", "Failed to generate report.", e)

    # ---- About dialog (English) ----

    def show_about(self):
        html = (
            "<div style='font-family:Segoe UI,Ubuntu,Arial;font-size:10px;line-height:1.6;'>"
            "<h2 style='color:#3498db;margin:0 0 4px 0;font-size:14px;'>MetaboDesk</h2>"
            "<p style='margin:0 0 8px 0;font-size:9px;'>"
            "Integrated Metabolic Network Analysis | SBML + COBRApy</p>"
            "<hr style='margin:8px 0;border:none;border-top:1px solid #ccc;'/>"
            "<h4 style='color:#2c3e50;margin:6px 0 4px 0;font-size:10px;'>ANALYSIS TYPES</h4>"
            "<table style='width:100%;border-collapse:collapse;font-size:9px;'>"
            "<tr><td style='padding:4px;font-weight:bold;color:#3498db;width:70px;'>FBA</td>"
            "<td>Flux Balance Analysis. Solves LP to find optimal fluxes. Choose objective, select FBA analysis, press Run.</td></tr>"
            "<tr style='background:#f9fafb;'><td style='padding:4px;font-weight:bold;color:#3498db;'>pFBA</td>"
            "<td>Minimizes total flux after optimization. Most parsimonious solution.</td></tr>"
            "<tr><td style='padding:4px;font-weight:bold;color:#3498db;'>FVA</td>"
            "<td>Flux Variability Analysis. Min/max bounds per reaction at optimal growth.</td></tr>"
            "<tr style='background:#f9fafb;'><td style='padding:4px;font-weight:bold;color:#3498db;'>SGD</td>"
            "<td>Single Gene Deletion. Knocks out each gene, identifies essential genes.</td></tr>"
            "<tr><td style='padding:4px;font-weight:bold;color:#3498db;'>SRD</td>"
            "<td>Single Reaction Deletion. Tests each reaction for essentiality.</td></tr>"
            "<tr style='background:#f9fafb;'><td style='padding:4px;font-weight:bold;color:#3498db;'>Robustness</td>"
            "<td>Parameter sweep. Varies reaction bounds, tracks objective response.</td></tr>"
            "<tr><td style='padding:4px;font-weight:bold;color:#3498db;'>Envelope</td>"
            "<td>Production envelope. Growth vs product trade-off Pareto curve.</td></tr>"
            "<tr style='background:#f9fafb;'><td style='padding:4px;font-weight:bold;color:#3498db;'>Sampling</td>"
            "<td>Monte Carlo flux sampling (ACHR). Distribution statistics per reaction.</td></tr>"
            "</table>"
            "<hr style='margin:8px 0;border:none;border-top:1px solid #ccc;'/>"
            "<h4 style='color:#2c3e50;margin:6px 0 4px 0;font-size:10px;'>QUICK START</h4>"
            "<ol style='margin:0;padding-left:16px;font-size:9px;'>"
            "<li>Open SBML (File → Open or Ctrl+O)</li>"
            "<li>Edit medium bounds (Medium tab, double-click cells)</li>"
            "<li>Select analysis type from dropdown</li>"
            "<li>Configure parameters</li>"
            "<li>Press Run Analysis (Ctrl+R)</li>"
            "<li>View results in tabs</li>"
            "</ol>"
            "<p style='margin:6px 0 0 0;color:#8c92a1;font-size:8px;'>© 2026 Emir Ay</p>"
            "</div>"
        )
        dlg = QMessageBox(self)
        dlg.setWindowTitle("About MetaboDesk")
        dlg.setText(html)
        dlg.setTextFormat(Qt.RichText)
        dlg.setMinimumWidth(600)
        dlg.setMinimumHeight(480)
        dlg.exec()

    # ---------------- Medium ----------------

    def show_metabolite_summary(self):
        """Show detailed summary for a selected metabolite: all producing/consuming reactions and net flux."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        # Metabolite selection dialog
        met_ids = sorted(m.id for m in self.base_model.metabolites)
        dialog = QDialog(self)
        dialog.setWindowTitle("Metabolite Summary")
        dialog.resize(800, 600)
        dlg_layout = QVBoxLayout(dialog)

        # Search + select
        sel_row = QHBoxLayout()
        sel_row.addWidget(QLabel("Metabolite:"))
        met_combo = QComboBox()
        met_combo.setEditable(True)
        met_combo.addItems(met_ids)
        met_combo.setCompleter(QCompleter(met_ids))
        sel_row.addWidget(met_combo, stretch=1)
        show_btn = QPushButton("Show Summary")
        sel_row.addWidget(show_btn)
        dlg_layout.addLayout(sel_row)

        info_label = QLabel("")
        info_label.setWordWrap(True)
        info_label.setTextFormat(Qt.RichText)
        dlg_layout.addWidget(info_label)

        table = QTableWidget()
        table.setColumnCount(5)
        table.setHorizontalHeaderLabels(["Reaction", "Name", "Stoich. Coeff.", "Flux", "Contribution"])
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectRows)
        table.setAlternatingRowColors(True)
        dlg_layout.addWidget(table, stretch=1)

        def _show_met():
            mid = met_combo.currentText().strip()
            if mid not in [m.id for m in self.base_model.metabolites]:
                info_label.setText(f"<b style='color:red;'>Metabolite '{mid}' not found.</b>")
                table.setRowCount(0)
                return

            met = self.base_model.metabolites.get_by_id(mid)

            # Get fluxes from last run if available
            flux_dict = {}
            if self.last_run:
                bl = self.last_run.get("baseline", {})
                flux_dict = bl.get("flux", bl.get("fluxes", {})) or {}

            producing = []
            consuming = []
            for rxn in met.reactions:
                coeff = rxn.metabolites[met]
                flux = flux_dict.get(rxn.id, 0.0)
                contribution = coeff * flux
                entry = (rxn.id, rxn.name or "", coeff, flux, contribution)
                if coeff > 0:
                    producing.append(entry)
                else:
                    consuming.append(entry)

            net_flux = sum(e[4] for e in producing) + sum(e[4] for e in consuming)

            met_name = met.name or ""
            compartment = getattr(met, "compartment", "") or ""
            formula = getattr(met, "formula", "") or ""

            info_html = (
                f"<h3 style='color:#2E86C1;'>{mid}</h3>"
                f"<b>Name:</b> {met_name} | <b>Compartment:</b> {compartment} | "
                f"<b>Formula:</b> {formula}<br>"
                f"<b>Producing reactions:</b> {len(producing)} | "
                f"<b>Consuming reactions:</b> {len(consuming)} | "
                f"<b>Net flux:</b> {net_flux:.6g}"
            )
            info_label.setText(info_html)

            all_entries = sorted(producing + consuming, key=lambda x: -abs(x[4]))
            table.setRowCount(len(all_entries))
            for i, (rid, rname, coeff, flux, contrib) in enumerate(all_entries):
                table.setItem(i, 0, QTableWidgetItem(rid))
                table.setItem(i, 1, QTableWidgetItem(rname))
                table.setItem(i, 2, QTableWidgetItem(f"{coeff:+.4g}"))
                table.setItem(i, 3, QTableWidgetItem(f"{flux:.6g}"))
                item = QTableWidgetItem(f"{contrib:.6g}")
                table.setItem(i, 4, item)

        show_btn.clicked.connect(_show_met)
        met_combo.activated.connect(lambda: _show_met())

        btn_box = QDialogButtonBox(QDialogButtonBox.Close)
        btn_box.rejected.connect(dialog.reject)
        dlg_layout.addWidget(btn_box)

        dialog.exec()

    def show_flux_distribution_comparison(self):
        """Compare flux distributions from multiple analyses side by side.
        
        Runs FBA, pFBA, and Loopless FBA in parallel and shows a comparison table
        so the user can see how each method distributes flux differently.
        """
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self.apply_reaction_overrides_to_model(model)
        self._apply_selected_objective_to_model(model)
        self._apply_selected_solver(model)

        def _compute():
            results = {}
            # Standard FBA
            try:
                sol = model.optimize()
                if sol.status == "optimal":
                    results["FBA"] = {"objective": float(sol.objective_value), "fluxes": sol.fluxes.to_dict()}
            except Exception as e:
                results["FBA"] = {"error": str(e)}

            # pFBA
            try:
                from cobra.flux_analysis import pfba
                psol = pfba(model)
                results["pFBA"] = {"objective": float(psol.objective_value), "fluxes": psol.fluxes.to_dict()}
            except Exception as e:
                results["pFBA"] = {"error": str(e)}

            # Loopless FBA
            try:
                from cobra.flux_analysis import loopless_solution
                lsol = loopless_solution(model)
                results["Loopless"] = {"objective": float(lsol.objective_value), "fluxes": lsol.fluxes.to_dict()}
            except Exception as e:
                results["Loopless"] = {"error": str(e)}

            return results

        def _on_done(results):
            self._show_flux_comparison_dialog(results)

        self._launch_worker(_compute, _on_done, "Running FBA / pFBA / Loopless comparison...")

    def _show_flux_comparison_dialog(self, results: dict):
        """Display side-by-side flux comparison from multiple analysis methods."""
        dialog = QDialog(self)
        dialog.setWindowTitle("Flux Distribution Comparison")
        dialog.resize(1000, 650)

        layout = QVBoxLayout(dialog)

        # Summary header
        methods = [m for m in ["FBA", "pFBA", "Loopless"] if m in results]
        summary_parts = []
        for m in methods:
            r = results[m]
            if "error" in r:
                summary_parts.append(f"<b>{m}:</b> <span style='color:red;'>Error — {r['error']}</span>")
            else:
                summary_parts.append(f"<b>{m}:</b> objective = {r['objective']:.6g}")
        summary_lbl = QLabel("<br>".join(summary_parts))
        summary_lbl.setWordWrap(True)
        summary_lbl.setTextFormat(Qt.RichText)
        layout.addWidget(summary_lbl)

        # Filter
        filter_row = QHBoxLayout()
        filter_row.addWidget(QLabel("Filter:"))
        filter_edit = QLineEdit()
        filter_edit.setPlaceholderText("reaction id contains...")
        filter_row.addWidget(filter_edit, stretch=1)
        nonzero_chk = QCheckBox("Non-zero only")
        nonzero_chk.setChecked(True)
        filter_row.addWidget(nonzero_chk)
        diff_chk = QCheckBox("Differences only")
        diff_chk.setToolTip("Show only reactions where methods disagree")
        filter_row.addWidget(diff_chk)
        layout.addLayout(filter_row)

        # Collect all reaction IDs
        all_rxn_ids = set()
        flux_data = {}
        for m in methods:
            if "fluxes" in results[m]:
                flux_data[m] = results[m]["fluxes"]
                all_rxn_ids.update(flux_data[m].keys())

        # Build table
        ncols = 1 + len(flux_data) + 1  # Reaction, method fluxes..., Max Δ
        col_headers = ["Reaction"] + list(flux_data.keys()) + ["Max Δ"]
        table = QTableWidget()
        table.setColumnCount(ncols)
        table.setHorizontalHeaderLabels(col_headers)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectRows)
        table.setAlternatingRowColors(True)

        # Pre-compute rows
        row_data_list = []
        for rid in sorted(all_rxn_ids):
            fluxes = [flux_data[m].get(rid, 0.0) for m in flux_data]
            max_delta = max(fluxes) - min(fluxes) if fluxes else 0.0
            row_data_list.append((rid, fluxes, max_delta))

        # Sort by max delta descending
        row_data_list.sort(key=lambda x: -abs(x[2]))

        table.setRowCount(len(row_data_list))
        for i, (rid, fluxes, max_delta) in enumerate(row_data_list):
            table.setItem(i, 0, QTableWidgetItem(rid))
            for j, fv in enumerate(fluxes):
                table.setItem(i, 1 + j, QTableWidgetItem(f"{fv:.6g}"))
            table.setItem(i, ncols - 1, QTableWidgetItem(f"{max_delta:.6g}"))

        def _apply_filter():
            text = filter_edit.text().strip().lower()
            only_nonzero = nonzero_chk.isChecked()
            only_diff = diff_chk.isChecked()
            for row in range(table.rowCount()):
                rid, fluxes, max_delta = row_data_list[row]
                visible = True
                if text and text not in rid.lower():
                    visible = False
                if only_nonzero and all(abs(f) < 1e-9 for f in fluxes):
                    visible = False
                if only_diff and abs(max_delta) < 1e-9:
                    visible = False
                table.setRowHidden(row, not visible)

        filter_edit.textChanged.connect(lambda: _apply_filter())
        nonzero_chk.stateChanged.connect(lambda: _apply_filter())
        diff_chk.stateChanged.connect(lambda: _apply_filter())
        _apply_filter()

        layout.addWidget(table, stretch=1)

        # Export buttons
        btn_row = QHBoxLayout()
        export_csv_btn = QPushButton("Export CSV")
        export_latex_btn = QPushButton("Export LaTeX")

        def _export_comparison_csv():
            fp, _ = QFileDialog.getSaveFileName(dialog, "Export CSV", "", "CSV (*.csv)")
            if not fp:
                return
            with open(fp, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(col_headers)
                for rid, fluxes, max_delta in row_data_list:
                    w.writerow([rid] + [f"{fv:.6g}" for fv in fluxes] + [f"{max_delta:.6g}"])
            self.statusBar().showMessage(f"Comparison exported: {fp}")

        def _export_comparison_latex():
            fp, _ = QFileDialog.getSaveFileName(dialog, "Export LaTeX", "", "TeX (*.tex)")
            if not fp:
                return
            col_spec = "l" + "r" * (len(col_headers) - 1)
            lines = [
                "\\begin{table}[htbp]",
                "\\centering",
                "\\caption{Flux Distribution Comparison (FBA vs pFBA vs Loopless)}",
                "\\label{tab:flux_comparison}",
                f"\\begin{{tabular}}{{{col_spec}}}",
                "\\hline",
                " & ".join(f"\\textbf{{{h}}}" for h in col_headers) + " \\\\",
                "\\hline",
            ]
            for rid, fluxes, max_delta in row_data_list[:50]:
                if all(abs(f) < 1e-9 for f in fluxes):
                    continue
                rid_tex = rid.replace("_", "\\_")
                vals = [rid_tex] + [f"{fv:.4g}" for fv in fluxes] + [f"{max_delta:.4g}"]
                lines.append(" & ".join(vals) + " \\\\")
            lines += ["\\hline", "\\end{tabular}", "\\end{table}"]
            with open(fp, "w", encoding="utf-8") as f:
                f.write("\n".join(lines))
            self.statusBar().showMessage(f"LaTeX comparison exported: {fp}")

        export_csv_btn.clicked.connect(_export_comparison_csv)
        export_latex_btn.clicked.connect(_export_comparison_latex)
        btn_row.addWidget(export_csv_btn)
        btn_row.addWidget(export_latex_btn)
        btn_row.addStretch()
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)

        dialog.exec()

    def _show_shadow_prices(self):
        """Show shadow prices from last FBA in a dialog with table and chart."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return
        if not self.last_run:
            QMessageBox.information(self, "No results", "Run FBA first to see shadow prices.")
            return

        # Use cached shadow prices from last FBA run if available
        cached = (self.last_run.get("baseline", {}).get("shadow_prices")
                  or self.last_run.get("original", {}).get("shadow_prices"))
        if cached:
            shadow_prices_dict = cached
        else:
            # Fallback: re-optimize
            try:
                model = self.base_model.copy()
                self.apply_medium_table_bounds_to_model(model)
                self.apply_reaction_overrides_to_model(model)
                self._apply_selected_solver(model)
                sol = model.optimize()
                if sol.status != 'optimal':
                    QMessageBox.warning(self, "Shadow Prices", f"Solution status: {sol.status}.")
                    return
                shadow_prices_dict = sol.shadow_prices.to_dict()
            except Exception as e:
                self._show_error("Error", "Failed to compute shadow prices.", e)
                return

        if not shadow_prices_dict:
            QMessageBox.information(self, "Shadow Prices", "No shadow prices available.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Shadow Prices (Metabolite Dual Values)")
        dialog.resize(600, 500)
        layout = QVBoxLayout(dialog)

        info = QLabel("Shadow prices indicate how much the objective would change per unit increase in metabolite availability.\n"
                       f"Source: {'cached from last FBA run' if cached else 're-optimized model'}")
        info.setWordWrap(True)
        layout.addWidget(info)

        search = QLineEdit()
        search.setPlaceholderText("Search metabolites...")
        layout.addWidget(search)

        table = QTableWidget(0, 3)
        table.setHorizontalHeaderLabels(["Metabolite", "Name", "Shadow Price"])
        table.horizontalHeader().setStretchLastSection(True)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        layout.addWidget(table)

        data = sorted(shadow_prices_dict.items(), key=lambda x: abs(x[1]), reverse=True)
        for met_id, price in data:
            try:
                met = self.base_model.metabolites.get_by_id(met_id)
                met_name = met.name if met else ""
            except Exception:
                met_name = ""
            r = table.rowCount()
            table.insertRow(r)
            table.setItem(r, 0, QTableWidgetItem(met_id))
            table.setItem(r, 1, QTableWidgetItem(met_name))
            table.setItem(r, 2, QTableWidgetItem(f"{price:.6g}"))

        def filter_table():
            txt = search.text().lower()
            for row in range(table.rowCount()):
                mid = (table.item(row, 0).text() or "").lower()
                mname = (table.item(row, 1).text() or "").lower()
                table.setRowHidden(row, txt not in mid and txt not in mname)
        search.textChanged.connect(filter_table)

        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)
        dialog.exec()

    def _show_reduced_costs(self):
        """Show reduced costs from last FBA in a dialog with table and chart."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return
        if not self.last_run:
            QMessageBox.information(self, "No results", "Run FBA first to see reduced costs.")
            return

        # Use cached reduced costs from last FBA run if available
        cached = (self.last_run.get("baseline", {}).get("reduced_costs")
                  or self.last_run.get("original", {}).get("reduced_costs"))
        if cached:
            reduced_costs_dict = cached
        else:
            try:
                model = self.base_model.copy()
                self.apply_medium_table_bounds_to_model(model)
                self.apply_reaction_overrides_to_model(model)
                self._apply_selected_solver(model)
                sol = model.optimize()
                if sol.status != 'optimal':
                    QMessageBox.warning(self, "Reduced Costs", f"Solution status: {sol.status}.")
                    return
                reduced_costs_dict = sol.reduced_costs.to_dict()
            except Exception as e:
                self._show_error("Error", "Failed to compute reduced costs.", e)
                return

        if not reduced_costs_dict:
            QMessageBox.information(self, "Reduced Costs", "No reduced costs available.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Reduced Costs (Reaction Dual Values)")
        dialog.resize(700, 500)
        layout = QVBoxLayout(dialog)

        info = QLabel("Reduced costs indicate how much the objective would change per unit increase in a reaction's flux bound.\n"
                       f"Source: {'cached from last FBA run' if cached else 're-optimized model'}")
        info.setWordWrap(True)
        layout.addWidget(info)

        search = QLineEdit()
        search.setPlaceholderText("Search reactions...")
        layout.addWidget(search)

        table = QTableWidget(0, 3)
        table.setHorizontalHeaderLabels(["Reaction", "Name", "Reduced Cost"])
        table.horizontalHeader().setStretchLastSection(True)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        layout.addWidget(table)

        data = sorted(reduced_costs_dict.items(), key=lambda x: abs(x[1]), reverse=True)
        for rxn_id, cost in data:
            try:
                rxn = self.base_model.reactions.get_by_id(rxn_id)
                rxn_name = rxn.name if rxn else ""
            except Exception:
                rxn_name = ""
            r = table.rowCount()
            table.insertRow(r)
            table.setItem(r, 0, QTableWidgetItem(rxn_id))
            table.setItem(r, 1, QTableWidgetItem(rxn_name))
            table.setItem(r, 2, QTableWidgetItem(f"{cost:.6g}"))

        def filter_table():
            txt = search.text().lower()
            for row in range(table.rowCount()):
                rid = (table.item(row, 0).text() or "").lower()
                rname = (table.item(row, 1).text() or "").lower()
                table.setRowHidden(row, txt not in rid and txt not in rname)
        search.textChanged.connect(filter_table)

        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)
        dialog.exec()

    def _find_essential_reactions(self):
        """Run single-reaction deletion to find essential reactions."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        reply = QMessageBox.question(
            self, "Essential Reactions",
            "This will test each reaction by setting its bounds to zero and checking if growth is still possible.\n\n"
            "This may take a few minutes for large models. Continue?",
            QMessageBox.Yes | QMessageBox.No
        )
        if reply != QMessageBox.Yes:
            return

        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self._apply_selected_solver(model)

        def _compute():
            from cobra.flux_analysis import single_reaction_deletion
            baseline_sol = model.optimize()
            if baseline_sol.status != 'optimal':
                raise ValueError("Baseline solution is not optimal.")
            baseline_growth = baseline_sol.objective_value
            result = single_reaction_deletion(model)
            essential = []
            for idx, row in result.iterrows():
                rxn_ids = list(idx) if hasattr(idx, '__iter__') and not isinstance(idx, str) else [idx]
                growth = row.get('growth', row.get('objective', 0)) or 0
                status = row.get('status', 'unknown')
                if growth < baseline_growth * 0.01 or status != 'optimal':
                    for rxn_id in rxn_ids:
                        rxn = model.reactions.get_by_id(rxn_id)
                        essential.append({
                            'id': rxn_id, 'name': rxn.name,
                            'subsystem': rxn.subsystem or "", 'growth': growth
                        })
            return {"essential": essential, "baseline_growth": baseline_growth}

        def _on_done(result):
            essential = result["essential"]
            baseline_growth = result["baseline_growth"]
            dialog = QDialog(self)
            dialog.setWindowTitle(f"Essential Reactions ({len(essential)} found)")
            dialog.resize(800, 500)
            layout = QVBoxLayout(dialog)
            info = QLabel(f"Found {len(essential)} essential reactions (knockout -> growth < 1% of baseline).\nBaseline growth: {baseline_growth:.6g}")
            info.setWordWrap(True)
            layout.addWidget(info)
            table = QTableWidget(0, 4)
            table.setHorizontalHeaderLabels(["Reaction", "Name", "Subsystem", "Growth after KO"])
            table.horizontalHeader().setStretchLastSection(True)
            table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            layout.addWidget(table)
            for item in essential:
                r = table.rowCount()
                table.insertRow(r)
                table.setItem(r, 0, QTableWidgetItem(item['id']))
                table.setItem(r, 1, QTableWidgetItem(item['name']))
                table.setItem(r, 2, QTableWidgetItem(item['subsystem']))
                table.setItem(r, 3, QTableWidgetItem(f"{item['growth']:.6g}"))
            btns = QDialogButtonBox(QDialogButtonBox.Close)
            btns.rejected.connect(dialog.reject)
            layout.addWidget(btns)
            dialog.exec()

        self._launch_worker(_compute, _on_done, "Finding essential reactions...")

    def _quick_fva(self):
        """Open a quick-FVA options dialog and run FVA in a worker thread."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Quick FVA Options")
        layout = QVBoxLayout(dialog)

        form = QFormLayout()
        fraction_spin = QDoubleSpinBox()
        fraction_spin.setRange(0.0, 1.0)
        fraction_spin.setValue(0.9)
        fraction_spin.setSingleStep(0.05)
        form.addRow("Fraction of optimum:", fraction_spin)

        rxn_filter = QLineEdit()
        rxn_filter.setPlaceholderText("Filter reactions (e.g., 'EX_' for exchange)")
        form.addRow("Reaction filter:", rxn_filter)

        layout.addLayout(form)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        layout.addWidget(btns)

        if dialog.exec() != QDialog.Accepted:
            return

        fraction = fraction_spin.value()
        filter_text = rxn_filter.text().strip().lower()

        model = self.base_model.copy()
        self.apply_medium_table_bounds_to_model(model)
        self._apply_selected_solver(model)

        def _compute():
            from cobra.flux_analysis import flux_variability_analysis
            if filter_text:
                rxns = [r for r in model.reactions if filter_text in r.id.lower()]
            else:
                rxns = model.reactions
            fva_result = flux_variability_analysis(model, reaction_list=rxns, fraction_of_optimum=fraction)
            rows = []
            for rxn_id, row in fva_result.iterrows():
                rows.append({"rxn_id": rxn_id, "min": row['minimum'], "max": row['maximum'],
                             "range": row['maximum'] - row['minimum']})
            return {"fva_rows": rows, "fraction": fraction}

        def _on_done(result):
            fva_rows = result["fva_rows"]
            result_dialog = QDialog(self)
            result_dialog.setWindowTitle(f"FVA Results ({len(fva_rows)} reactions)")
            result_dialog.resize(700, 500)
            layout2 = QVBoxLayout(result_dialog)

            info = QLabel(f"FVA at {result['fraction'] * 100:.0f}% of optimum. Showing min/max flux ranges.")
            layout2.addWidget(info)

            search = QLineEdit()
            search.setPlaceholderText("Search reactions...")
            layout2.addWidget(search)

            table = QTableWidget(0, 4)
            table.setHorizontalHeaderLabels(["Reaction", "Min Flux", "Max Flux", "Range"])
            table.horizontalHeader().setStretchLastSection(True)
            table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            layout2.addWidget(table)

            for frow in fva_rows:
                r = table.rowCount()
                table.insertRow(r)
                table.setItem(r, 0, QTableWidgetItem(frow["rxn_id"]))
                table.setItem(r, 1, QTableWidgetItem(f"{frow['min']:.6g}"))
                table.setItem(r, 2, QTableWidgetItem(f"{frow['max']:.6g}"))
                table.setItem(r, 3, QTableWidgetItem(f"{frow['range']:.6g}"))

            def filter_table():
                txt = search.text().lower()
                for row in range(table.rowCount()):
                    rxn_id = (table.item(row, 0).text() or "").lower()
                    table.setRowHidden(row, txt and txt not in rxn_id)
            search.textChanged.connect(filter_table)

            btns2 = QDialogButtonBox(QDialogButtonBox.Close)
            btns2.rejected.connect(result_dialog.reject)
            layout2.addWidget(btns2)
            result_dialog.exec()

        self._launch_worker(_compute, _on_done, "Running Quick FVA...")

    # ---------------- Run FBA ----------------

    def _apply_selected_solver(self, model: cobra.Model):
        """Apply selected solver to COBRApy model."""
        try:
            s = (self.solver_combo.currentText() or "Auto").lower()
            if s == "auto":
                return
            if s == "glpk":
                model.solver = "glpk"
            elif s == "highs":
                try:
                    import optlang.scipy_interface as iface
                    model.solver = iface
                except ImportError:
                    model.solver = "glpk"
                    self.statusBar().showMessage("HiGHS not available, falling back to GLPK.")
            elif s == "gurobi":
                try:
                    model.solver = "gurobi"
                except Exception:
                    self.statusBar().showMessage("Gurobi not available — install gurobipy with a valid license.")
                    return
            elif s == "cplex":
                try:
                    model.solver = "cplex"
                except Exception:
                    self.statusBar().showMessage("CPLEX not available — install cplex with a valid license.")
                    return
            else:
                return
        except Exception as e:
            # Graceful fallback
            self.statusBar().showMessage(f"Solver set failed ({s}): {e}")

    def _gene_name_from_id(self, gid: str) -> str:
        """Look up a gene's name from its ID in the loaded model."""
        if not gid or self.base_model is None:
            return ""
        try:
            g = self.base_model.genes.get_by_id(str(gid))
            return str(getattr(g, "name", "") or "").strip()
        except Exception:
            return ""

    def _rxn_name_from_id(self, rid: str) -> str:
        """Look up a reaction's name from its ID in the loaded model."""
        if not rid or self.base_model is None:
            return ""
        try:
            r = self.base_model.reactions.get_by_id(str(rid))
            return str(getattr(r, "name", "") or "").strip()
        except Exception:
            return ""

    def _format_gene_label(self, gid: str) -> str:
        """Format a gene ID as ``id — name`` for display."""
        name = self._gene_name_from_id(gid)
        if name:
            return f"{gid} — {name}".strip(" —")
        return str(gid)

    def _format_rxn_label(self, rid: str) -> str:
        """Format a reaction ID as ``id — name`` for display."""
        name = self._rxn_name_from_id(rid)
        if name:
            return f"{rid} — {name}".strip(" —")
        return str(rid)

    def _format_gene_pair_label(self, pair_id: str) -> str:
        """Format a gene-pair ID (e.g. ``g1+g2``) with names for display."""
        text = str(pair_id or "").strip()
        if not text:
            return ""
        if "+" in text:
            parts = [p.strip() for p in text.split("+") if p.strip()]
        elif "-" in text:
            parts = [p.strip() for p in text.split("-") if p.strip()]
        else:
            parts = [text]

        if len(parts) >= 2:
            g1 = self._format_gene_label(parts[0])
            g2 = self._format_gene_label(parts[1])
            return f"{g1} + {g2}"
        return self._format_gene_label(parts[0])

    def _extract_reaction_id_from_input(self, text: str) -> str:
        """Strip the display name from a ``id — name`` string, returning just the ID."""
        t = (text or "").strip()
        if not t:
            return ""
        for sep in (" — ", " - ", " – ", " —", "- "):
            if sep in t:
                return t.split(sep, 1)[0].strip()
        return t

    # ==================== NEW FEATURES ====================

    # -------- 1. EXPORT FUNCTIONS --------

    def show_search_dialog(self):
        """Open search dialog for metabolites and reactions"""
        dialog = QDialog(self)
        dialog.setWindowTitle("Search Metabolites & Reactions")
        dialog.resize(500, 400)
        
        layout = QVBoxLayout()
        
        search_input = QLineEdit()
        search_input.setPlaceholderText("Enter metabolite or reaction ID...")
        layout.addWidget(QLabel("Search:"))
        layout.addWidget(search_input)
        
        results_table = QTableWidget()
        results_table.setColumnCount(3)
        results_table.setHorizontalHeaderLabels(["Type", "ID", "Name"])
        results_table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(results_table)
        
        def perform_search():
            results_table.setRowCount(0)
            if not self.base_model or not search_input.text():
                return
            
            query = search_input.text().lower()
            
            # Search metabolites
            for met in self.base_model.metabolites:
                if query in met.id.lower() or query in met.name.lower():
                    row = results_table.rowCount()
                    results_table.insertRow(row)
                    results_table.setItem(row, 0, QTableWidgetItem("Metabolite"))
                    results_table.setItem(row, 1, QTableWidgetItem(met.id))
                    results_table.setItem(row, 2, QTableWidgetItem(met.name))
            
            # Search reactions
            for rxn in self.base_model.reactions:
                if query in rxn.id.lower() or query in rxn.name.lower():
                    row = results_table.rowCount()
                    results_table.insertRow(row)
                    results_table.setItem(row, 0, QTableWidgetItem("Reaction"))
                    results_table.setItem(row, 1, QTableWidgetItem(rxn.id))
                    results_table.setItem(row, 2, QTableWidgetItem(rxn.name))
        
        search_input.textChanged.connect(perform_search)
        
        buttons = QHBoxLayout()
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        buttons.addStretch()
        buttons.addWidget(close_btn)
        layout.addLayout(buttons)
        
        dialog.setLayout(layout)
        dialog.exec()

    # -------- 3. SENSITIVITY ANALYSIS --------

    def run_sensitivity_analysis(self):
        """Run sensitivity analysis on a chosen reaction parameter."""
        if self.base_model is None:
            QMessageBox.warning(self, "Error", "Load a model first.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Sensitivity Analysis")
        dialog.resize(500, 340)

        layout = QFormLayout()

        # Reaction selector with search
        rxn_combo = QComboBox()
        rxn_combo.setEditable(True)
        rxn_combo.setInsertPolicy(QComboBox.NoInsert)
        obj_id = self._get_objective_reaction_id() or ""
        rxn_ids = sorted(r.id for r in self.base_model.reactions)
        rxn_combo.addItems(rxn_ids)
        if obj_id and obj_id in rxn_ids:
            rxn_combo.setCurrentText(obj_id)
        rxn_combo.setCompleter(QCompleter(rxn_ids, rxn_combo))
        layout.addRow("Reaction:", rxn_combo)

        # Which bound to sweep
        bound_combo = QComboBox()
        bound_combo.addItems(["Lower bound", "Upper bound"])
        layout.addRow("Parameter:", bound_combo)

        min_spin = QDoubleSpinBox()
        min_spin.setRange(-1000, 1000)
        min_spin.setDecimals(4)
        min_spin.setValue(-100)
        layout.addRow("Min value:", min_spin)

        max_spin = QDoubleSpinBox()
        max_spin.setRange(-1000, 1000)
        max_spin.setDecimals(4)
        max_spin.setValue(100)
        layout.addRow("Max value:", max_spin)

        steps_spin = QSpinBox()
        steps_spin.setRange(5, 500)
        steps_spin.setValue(20)
        layout.addRow("Steps:", steps_spin)

        buttons = QHBoxLayout()
        ok_btn = QPushButton("Run")
        cancel_btn = QPushButton("Cancel")

        def run_sens():
            rid = rxn_combo.currentText().strip()
            if rid not in self.base_model.reactions:
                QMessageBox.warning(dialog, "Invalid", f"Reaction '{rid}' not found in model.")
                return
            bound_type = "lb" if bound_combo.currentIndex() == 0 else "ub"
            dialog.close()
            self._perform_sensitivity_analysis(
                rid, bound_type,
                min_spin.value(), max_spin.value(), steps_spin.value(),
            )

        ok_btn.clicked.connect(run_sens)
        cancel_btn.clicked.connect(dialog.close)
        buttons.addStretch()
        buttons.addWidget(ok_btn)
        buttons.addWidget(cancel_btn)

        layout.addRow(buttons)
        dialog.setLayout(layout)
        dialog.exec()

    def _perform_sensitivity_analysis(self, rxn_id: str, bound_type: str,
                                       min_val: float, max_val: float, steps: int):
        """Sweep *bound_type* ('lb' or 'ub') of *rxn_id* and record objective.
        
        Runs the heavy computation in a QThread worker so the UI stays responsive.
        """
        if self.base_model is None or self.base_model.objective is None:
            return

        if rxn_id not in self.base_model.reactions:
            self._show_error("Error", f"Reaction '{rxn_id}' not found.")
            return

        # Prepare a copy so we don't mutate base_model from the worker thread
        model_copy = self.base_model.copy()
        self._apply_selected_solver(model_copy)

        def _compute(worker=None):
            values: list[float] = []
            objectives: list[float | None] = []

            for i in range(steps):
                if worker:
                    worker.report_progress(
                        f"Sensitivity step {i + 1}/{steps}",
                        int(i / max(steps - 1, 1) * 100),
                    )
                param_val = min_val + (max_val - min_val) * (i / max(steps - 1, 1))
                with model_copy:
                    try:
                        rxn = model_copy.reactions.get_by_id(rxn_id)
                        if bound_type == "lb":
                            rxn.lower_bound = param_val
                        else:
                            rxn.upper_bound = param_val
                    except Exception:
                        objectives.append(None)
                        values.append(param_val)
                        continue
                    try:
                        sol = model_copy.optimize()
                        objectives.append(sol.objective_value if sol.status == "optimal" else None)
                    except Exception:
                        objectives.append(None)

                values.append(param_val)

            return {
                "reaction_id": rxn_id,
                "bound_type": bound_type,
                "parameter_values": values,
                "objective_values": objectives,
            }

        def _on_done(result):
            self.sensitivity_results = result
            values = result["parameter_values"]
            objectives = result["objective_values"]
            if values:
                label = f"{rxn_id} {'lower' if bound_type == 'lb' else 'upper'} bound"
                self._show_sensitivity_results(values, objectives, param_label=label)
            else:
                QMessageBox.information(self, "Sensitivity", "No results.")

        self._launch_worker(_compute, _on_done, f"Running sensitivity analysis on {rxn_id}...")

    # -------- 4. BATCH ANALYSIS --------

    def run_batch_analysis(self):
        """Run analysis on multiple SBML files with comprehensive comparison."""
        files, _ = QFileDialog.getOpenFileNames(self, "Select SBML Files", "", "SBML Files (*.xml *.sbml)")
        if not files:
            return
        
        # Options dialog
        opts = QDialog(self)
        opts.setWindowTitle("Batch Analysis Options")
        opts_layout = QFormLayout(opts)

        run_fva_chk = QCheckBox("Include FVA (min/max ranges)")
        run_fva_chk.setToolTip("Compute FVA for each model (slower)")
        opts_layout.addRow(run_fva_chk)

        run_sgd_chk = QCheckBox("Include essential gene count")
        run_sgd_chk.setToolTip("Run single gene deletion to find essential genes (slower)")
        opts_layout.addRow(run_sgd_chk)

        opts_btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        opts_btns.accepted.connect(opts.accept)
        opts_btns.rejected.connect(opts.reject)
        opts_layout.addRow(opts_btns)

        if opts.exec() != QDialog.Accepted:
            return

        do_fva = run_fva_chk.isChecked()
        do_sgd = run_sgd_chk.isChecked()

        def _compute(worker=None):
            results = []
            for i, file_path in enumerate(files):
                if worker:
                    worker.report_progress(
                        f"Analyzing model {i + 1}/{len(files)}: {Path(file_path).name}",
                        int(i / len(files) * 100),
                    )
                try:
                    model = cobra.io.read_sbml_model(file_path)
                    info = {
                        "file": Path(file_path).name,
                        "reactions": len(model.reactions),
                        "metabolites": len(model.metabolites),
                        "genes": len(model.genes),
                        "exchange_rxns": sum(1 for r in model.reactions if len(r.metabolites) == 1),
                    }
                    with model:
                        solution = model.optimize()
                        info["status"] = solution.status if solution else "No solution"
                        info["objective"] = float(solution.objective_value) if solution and solution.status == "optimal" else None
                        info["objective_rxn"] = str(model.objective.expression) if model.objective else ""

                    if do_fva and info.get("status") == "optimal":
                        try:
                            from cobra.flux_analysis import flux_variability_analysis
                            fva = flux_variability_analysis(model, processes=1)
                            ranges = (fva["maximum"] - fva["minimum"]).abs()
                            info["fva_flexible"] = int((ranges > 1e-6).sum())
                            info["fva_fixed"] = int((ranges <= 1e-6).sum())
                        except Exception:
                            info["fva_flexible"] = "Error"
                            info["fva_fixed"] = "Error"

                    if do_sgd and info.get("status") == "optimal":
                        try:
                            from cobra.flux_analysis import single_gene_deletion
                            sgd = single_gene_deletion(model, processes=1)
                            growth_col = "growth" if "growth" in sgd.columns else "growth_rate"
                            wt_growth = float(info["objective"])
                            essential = int((sgd[growth_col].fillna(0) < wt_growth * 0.01).sum())
                            info["essential_genes"] = essential
                        except Exception:
                            info["essential_genes"] = "Error"

                    results.append(info)
                except Exception as e:
                    results.append({
                        "file": Path(file_path).name,
                        "reactions": 0, "metabolites": 0, "genes": 0,
                        "exchange_rxns": 0,
                        "status": "Error",
                        "objective": str(e),
                        "objective_rxn": "",
                    })
            return {"results": results, "do_fva": do_fva, "do_sgd": do_sgd}

        def _on_done(result):
            self._show_batch_results(result["results"], result["do_fva"], result["do_sgd"])

        self._launch_worker(_compute, _on_done, f"Running batch analysis on {len(files)} models...")

    def _show_batch_results(self, results, has_fva=False, has_sgd=False):
        """Display batch analysis results with comprehensive table and export."""
        dialog = QDialog(self)
        dialog.setWindowTitle("Batch Analysis Results")
        dialog.resize(1000, 500)
        
        layout = QVBoxLayout()
        
        # Summary
        ok_count = sum(1 for r in results if r.get("status") == "optimal")
        layout.addWidget(QLabel(f"Models analyzed: {len(results)} | Optimal: {ok_count} | Failed: {len(results) - ok_count}"))

        # Build columns
        cols = ["File", "Reactions", "Metabolites", "Genes", "Exchange Rxns", "Status", "Objective Value"]
        if has_fva:
            cols += ["FVA Flexible", "FVA Fixed"]
        if has_sgd:
            cols += ["Essential Genes"]
        
        table = QTableWidget()
        table.setColumnCount(len(cols))
        table.setHorizontalHeaderLabels(cols)
        table.horizontalHeader().setStretchLastSection(True)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectRows)
        
        for result in results:
            row = table.rowCount()
            table.insertRow(row)
            c = 0
            table.setItem(row, c, QTableWidgetItem(str(result.get("file", "")))); c += 1
            table.setItem(row, c, QTableWidgetItem(str(result.get("reactions", 0)))); c += 1
            table.setItem(row, c, QTableWidgetItem(str(result.get("metabolites", 0)))); c += 1
            table.setItem(row, c, QTableWidgetItem(str(result.get("genes", 0)))); c += 1
            table.setItem(row, c, QTableWidgetItem(str(result.get("exchange_rxns", 0)))); c += 1
            table.setItem(row, c, QTableWidgetItem(str(result.get("status", "")))); c += 1
            obj = result.get("objective")
            table.setItem(row, c, QTableWidgetItem(f"{obj:.6g}" if isinstance(obj, (int, float)) else str(obj or ""))); c += 1
            if has_fva:
                table.setItem(row, c, QTableWidgetItem(str(result.get("fva_flexible", "")))); c += 1
                table.setItem(row, c, QTableWidgetItem(str(result.get("fva_fixed", "")))); c += 1
            if has_sgd:
                table.setItem(row, c, QTableWidgetItem(str(result.get("essential_genes", "")))); c += 1
        
        layout.addWidget(table, stretch=1)
        
        buttons = QHBoxLayout()
        export_csv_btn = QPushButton("Export CSV")
        export_xlsx_btn = QPushButton("Export Excel")

        def _export_batch_csv():
            fp, _ = QFileDialog.getSaveFileName(dialog, "Export Batch (CSV)", str(Path.home() / "batch_results.csv"), "CSV (*.csv)")
            if not fp:
                return
            with open(fp, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(cols)
                for r in results:
                    row_data = [r.get("file", ""), r.get("reactions", 0), r.get("metabolites", 0),
                               r.get("genes", 0), r.get("exchange_rxns", 0), r.get("status", ""),
                               r.get("objective", "")]
                    if has_fva:
                        row_data += [r.get("fva_flexible", ""), r.get("fva_fixed", "")]
                    if has_sgd:
                        row_data += [r.get("essential_genes", "")]
                    w.writerow(row_data)
            self.statusBar().showMessage(f"Batch exported: {fp}")

        def _export_batch_xlsx():
            fp, _ = QFileDialog.getSaveFileName(dialog, "Export Batch (Excel)", str(Path.home() / "batch_results.xlsx"), "Excel (*.xlsx)")
            if not fp:
                return
            try:
                import openpyxl
                wb = openpyxl.Workbook()
                ws = wb.active
                ws.title = "Batch Results"
                ws.append(cols)
                for r in results:
                    row_data = [r.get("file", ""), r.get("reactions", 0), r.get("metabolites", 0),
                               r.get("genes", 0), r.get("exchange_rxns", 0), r.get("status", ""),
                               r.get("objective", "")]
                    if has_fva:
                        row_data += [r.get("fva_flexible", ""), r.get("fva_fixed", "")]
                    if has_sgd:
                        row_data += [r.get("essential_genes", "")]
                    ws.append(row_data)
                wb.save(fp)
                self.statusBar().showMessage(f"Batch exported: {fp}")
            except ImportError:
                QMessageBox.warning(dialog, "Missing", "openpyxl required for Excel export.")

        export_csv_btn.clicked.connect(_export_batch_csv)
        export_xlsx_btn.clicked.connect(_export_batch_xlsx)
        buttons.addWidget(export_csv_btn)
        buttons.addWidget(export_xlsx_btn)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        buttons.addStretch()
        buttons.addWidget(close_btn)
        layout.addLayout(buttons)
        
        dialog.setLayout(layout)
        dialog.exec()

    def _show_sensitivity_results(self, values, objectives, param_label: str = "Parameter value"):
        """Display sensitivity analysis results with chart and table in a dialog."""
        dialog = QDialog(self)
        dialog.setWindowTitle("Sensitivity Analysis Results")
        dialog.resize(800, 500)

        layout = QVBoxLayout()

        # Plot
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.figure import Figure
        fig = Figure(figsize=(6, 3.5), dpi=100)
        ax = fig.add_subplot(111)
        # Filter out None values (infeasible steps) before plotting
        plot_vals = [(v, o) for v, o in zip(values, objectives) if o is not None]
        if plot_vals:
            pv, po = zip(*plot_vals)
            ax.plot(pv, po, marker="o", linewidth=1.5)
        else:
            ax.text(0.5, 0.5, "All steps infeasible", transform=ax.transAxes, ha="center")
        ax.set_xlabel(param_label)
        ax.set_ylabel("Objective value")
        ax.grid(True, alpha=0.3)
        canvas = FigureCanvas(fig)
        layout.addWidget(canvas, stretch=1)

        # Table
        table = QTableWidget()
        table.setColumnCount(2)
        table.setHorizontalHeaderLabels(["Parameter", "Objective"])
        table.horizontalHeader().setStretchLastSection(True)
        for v, obj in zip(values, objectives):
            row = table.rowCount()
            table.insertRow(row)
            table.setItem(row, 0, QTableWidgetItem(f"{v:.6g}"))
            table.setItem(row, 1, QTableWidgetItem("" if obj is None else f"{obj:.6g}"))
        layout.addWidget(table, stretch=1)

        buttons = QHBoxLayout()
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        buttons.addStretch()
        buttons.addWidget(close_btn)
        layout.addLayout(buttons)

        dialog.setLayout(layout)
        dialog.exec()

    # -------- 4b. COMMUNITY ANALYSIS --------

    def run_community_analysis(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Select SBML Files (2+)", "", "SBML Files (*.xml *.sbml)")
        if not files or len(files) < 2:
            QMessageBox.warning(self, "Community", "Select at least two SBML files.")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Community Analysis Options")
        layout = QFormLayout(dialog)

        share_ext_chk = QCheckBox("Share extracellular metabolites")
        share_ext_chk.setChecked(True)
        layout.addRow(QLabel(""), share_ext_chk)

        obj_mode_combo = QComboBox()
        obj_mode_combo.addItems([
            "Sum of objectives (default)",
            "SteadyCom — equal growth rate",
        ])
        layout.addRow(QLabel("Objective mode:"), obj_mode_combo)

        topn_spin = QSpinBox()
        topn_spin.setRange(5, 200)
        topn_spin.setValue(20)
        layout.addRow(QLabel("Top exchange fluxes:"), topn_spin)

        # Abundance constraints
        abundance_edit = QLineEdit()
        abundance_edit.setPlaceholderText("e.g. 0.7, 0.3  (leave empty for equal)")
        abundance_edit.setToolTip(
            "Comma-separated fractional abundances for each species.\n"
            "Values should sum to 1.0 — e.g. for 2 organisms: 0.7, 0.3\n"
            "Leave empty to give all species equal capacity.")
        layout.addRow(QLabel("Relative abundance:"), abundance_edit)

        export_chk = QCheckBox("Export merged community SBML")
        export_chk.setChecked(False)
        layout.addRow(QLabel(""), export_chk)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dialog.accept)
        btns.rejected.connect(dialog.reject)
        layout.addRow(btns)

        if dialog.exec() != QDialog.Accepted:
            return

        labels = [Path(f).stem for f in files]
        share_ext = share_ext_chk.isChecked()
        obj_mode = "steadycom" if obj_mode_combo.currentIndex() == 1 else "sum"
        topn = topn_spin.value()
        do_export = export_chk.isChecked()

        # Parse abundance fractions
        abundance: dict[str, float] | None = None
        raw_abundance = abundance_edit.text().strip()
        if raw_abundance:
            try:
                fracs = [float(x.strip()) for x in raw_abundance.split(",") if x.strip()]
                if len(fracs) != len(files):
                    QMessageBox.warning(self, "Abundance",
                                        f"Expected {len(files)} values but got {len(fracs)}.\n"
                                        "Leave empty for equal abundance or provide one value per species.")
                    return
                total = sum(fracs)
                if abs(total - 1.0) > 0.05:
                    QMessageBox.warning(self, "Abundance",
                                        f"Abundance values should sum to 1.0 (got {total:.3f}).\n"
                                        "Please adjust the values.")
                    return
                abundance = dict(zip(labels, fracs))
            except ValueError:
                QMessageBox.warning(self, "Abundance", "Could not parse abundance values.\n"
                                    "Use comma-separated decimals, e.g. 0.7, 0.3")
                return

        def _compute(worker=None):
            models = []
            for j, f in enumerate(files):
                if worker:
                    worker.report_progress(
                        f"Loading model {j + 1}/{len(files)}: {Path(f).name}",
                        int(j / len(files) * 40),
                    )
                models.append(cobra.io.read_sbml_model(f))

            if worker:
                worker.report_progress("Building community model...", 50)
            community, obj_map, shared_mets, met_conflicts = self._build_community_model(
                models, labels, share_extracellular=share_ext,
                objective_mode=obj_mode, abundance=abundance)

            if worker:
                worker.report_progress("Optimizing community model...", 70)
            solution = community.optimize()
            total_obj = solution.objective_value if solution else None

            if worker:
                worker.report_progress("Analyzing per-species growth...", 85)
            # Per-species growth: solo vs community
            per_species = []
            solo_exchange_maps: dict[str, dict[str, float]] = {}
            for label, model in zip(labels, models):
                solo_val = None
                try:
                    sol = model.optimize()
                    solo_val = sol.objective_value if sol and sol.status == "optimal" else None
                except Exception:
                    sol = None
                    solo_val = None

                comm_val = None
                rxn_id = obj_map.get(label)
                if rxn_id:
                    try:
                        comm_val = solution.fluxes.get(rxn_id, None)
                    except Exception:
                        comm_val = None
                per_species.append((label, rxn_id or "", solo_val, comm_val))

                # Exchange map by metabolite id for solo model
                emap = {}
                try:
                    for rxn in model.reactions:
                        if len(rxn.metabolites) == 1:
                            met = list(rxn.metabolites.keys())[0]
                            met_id = met.id
                            v = sol.fluxes.get(rxn.id, 0.0) if sol else 0.0
                            if abs(v) > 1e-9:
                                emap[met_id] = emap.get(met_id, 0.0) + float(v)
                except Exception:
                    pass
                solo_exchange_maps[label] = emap

            # Community exchange fluxes (shared extracellular)
            comm_exchange = {}
            try:
                for rxn in community.reactions:
                    if len(rxn.metabolites) == 1:
                        met = list(rxn.metabolites.keys())[0]
                        met_id = met.id
                        v = solution.fluxes.get(rxn.id, 0.0)
                        if abs(v) > 1e-9:
                            comm_exchange[met_id] = comm_exchange.get(met_id, 0.0) + float(v)
            except Exception:
                pass

            # Delta exchange: community vs sum of solo
            delta_exch = []
            for met_id, comm_v in comm_exchange.items():
                solo_sum = 0.0
                for emap in solo_exchange_maps.values():
                    solo_sum += emap.get(met_id, 0.0)
                delta = comm_v - solo_sum
                delta_exch.append((met_id, comm_v, solo_sum, delta))
            delta_exch = sorted(delta_exch, key=lambda x: abs(x[3]), reverse=True)[:topn]

            return {"total_obj": total_obj, "per_species": per_species,
                    "delta_exch": delta_exch, "community": community,
                    "do_export": do_export, "met_conflicts": met_conflicts}

        def _on_done(result):
            if result.get("do_export"):
                out_path, _ = QFileDialog.getSaveFileName(
                    self,
                    "Save community SBML",
                    str(Path.home() / "community_model.xml"),
                    "SBML files (*.xml *.sbml);;All files (*.*)",
                )
                if out_path:
                    try:
                        cobra.io.write_sbml_model(result["community"], out_path)
                    except Exception as e:
                        self._show_error("Export failed", "Could not save community model.", e)

            # Show conflict warning if any
            conflicts = result.get("met_conflicts", [])
            if conflicts:
                msg = ("⚠ Metabolite conflicts detected in shared extracellular pool:\n\n"
                       + "\n".join(conflicts[:20]))
                if len(conflicts) > 20:
                    msg += f"\n... and {len(conflicts) - 20} more."
                QMessageBox.warning(self, "Community Model — Metabolite Conflicts", msg)

            self._show_community_results(result["total_obj"], result["per_species"], result["delta_exch"])

        self._launch_worker(_compute, _on_done, "Running Community Analysis...")

    def _show_community_results(self, total_obj, per_species, delta_exch):
        """Display community modeling results (growth chart, tables, exports) in a dialog."""
        dialog = QDialog(self)
        dialog.setWindowTitle("Community Analysis Results")
        dialog.resize(900, 600)

        layout = QVBoxLayout()
        layout.addWidget(QLabel(f"<b>Community Total Objective: {total_obj:.6g if isinstance(total_obj, (int, float)) else total_obj}</b>"))

        # Growth comparison chart
        try:
            from matplotlib.figure import Figure
            from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
            fig = Figure(figsize=(7, 2.5), dpi=100)
            ax = fig.add_subplot(111)
            labels_list = [ps[0] for ps in per_species]
            solo_vals = [ps[2] or 0 for ps in per_species]
            comm_vals = [ps[3] or 0 for ps in per_species]
            x = list(range(len(labels_list)))
            w = 0.35
            ax.bar([i - w/2 for i in x], solo_vals, w, label="Solo", color="#4CAF50", alpha=0.8)
            ax.bar([i + w/2 for i in x], comm_vals, w, label="Community", color="#2196F3", alpha=0.8)
            ax.set_xticks(x)
            ax.set_xticklabels(labels_list, rotation=45, ha="right", fontsize=8)
            ax.set_ylabel("Growth")
            ax.set_title("Solo vs Community Growth")
            ax.legend(fontsize=8)
            fig.tight_layout()
            canvas = FigureCanvas(fig)
            layout.addWidget(canvas)
        except Exception:
            pass

        # Species table
        layout.addWidget(QLabel("<b>Per-species growth comparison:</b>"))
        table = QTableWidget()
        table.setColumnCount(5)
        table.setHorizontalHeaderLabels(["Species", "Objective Rxn", "Solo Growth", "Community Growth", "Δ Growth"])
        table.horizontalHeader().setStretchLastSection(True)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        for label, rxn_id, solo_val, comm_val in per_species:
            row = table.rowCount()
            table.insertRow(row)
            table.setItem(row, 0, QTableWidgetItem(label))
            table.setItem(row, 1, QTableWidgetItem(rxn_id))
            table.setItem(row, 2, QTableWidgetItem("" if solo_val is None else f"{solo_val:.6g}"))
            table.setItem(row, 3, QTableWidgetItem("" if comm_val is None else f"{comm_val:.6g}"))
            delta = ""
            if solo_val is not None and comm_val is not None:
                delta = f"{comm_val - solo_val:.6g}"
            table.setItem(row, 4, QTableWidgetItem(delta))
        layout.addWidget(table)

        if delta_exch:
            layout.addWidget(QLabel("<b>Top exchange flux changes (community vs solo sum):</b>"))
            exch_table = QTableWidget()
            exch_table.setColumnCount(4)
            exch_table.setHorizontalHeaderLabels(["Metabolite", "Community Flux", "Solo Sum", "Delta"])
            exch_table.horizontalHeader().setStretchLastSection(True)
            exch_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            for met_id, comm_v, solo_sum, delta in delta_exch:
                row = exch_table.rowCount()
                exch_table.insertRow(row)
                exch_table.setItem(row, 0, QTableWidgetItem(met_id))
                exch_table.setItem(row, 1, QTableWidgetItem(f"{comm_v:.6g}"))
                exch_table.setItem(row, 2, QTableWidgetItem(f"{solo_sum:.6g}"))
                exch_table.setItem(row, 3, QTableWidgetItem(f"{delta:.6g}"))
            layout.addWidget(exch_table)

        buttons = QHBoxLayout()
        export_csv_btn = QPushButton("Export CSV")

        def _export_community_csv():
            fp, _ = QFileDialog.getSaveFileName(dialog, "Export Community (CSV)", str(Path.home() / "community_results.csv"), "CSV (*.csv)")
            if not fp:
                return
            with open(fp, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["Total Objective", total_obj])
                w.writerow([])
                w.writerow(["Species", "Objective Rxn", "Solo Growth", "Community Growth", "Delta"])
                for label, rxn_id, solo_val, comm_val in per_species:
                    delta = (comm_val - solo_val) if solo_val is not None and comm_val is not None else ""
                    w.writerow([label, rxn_id, solo_val or "", comm_val or "", delta])
                if delta_exch:
                    w.writerow([])
                    w.writerow(["Metabolite", "Community Flux", "Solo Sum", "Delta"])
                    for met_id, comm_v, solo_sum, delta in delta_exch:
                        w.writerow([met_id, comm_v, solo_sum, delta])
            self.statusBar().showMessage(f"Community exported: {fp}")

        export_csv_btn.clicked.connect(_export_community_csv)
        buttons.addWidget(export_csv_btn)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        buttons.addStretch()
        buttons.addWidget(close_btn)
        layout.addLayout(buttons)

        dialog.setLayout(layout)
        dialog.exec()

    # -------- 5. REACTION/METABOLITE DETAILS --------

    def show_reaction_details(self, reaction_id):
        """Show details for a reaction"""
        if not self.base_model or reaction_id not in self.base_model.reactions:
            return
        
        rxn = self.base_model.reactions.get_by_id(reaction_id)
        
        dialog = QDialog(self)
        dialog.setWindowTitle(f"Reaction Details: {reaction_id}")
        dialog.resize(500, 400)
        
        layout = QVBoxLayout()
        
        text = QPlainTextEdit()
        text.setReadOnly(True)
        
        details = f"""
ID: {rxn.id}
Name: {rxn.name}
Bounds: {rxn.lower_bound} to {rxn.upper_bound}
Gene Reaction Rule: {rxn.gene_reaction_rule}
Reversible: {rxn.reversible}

Metabolites:
"""
        for met, coeff in rxn.metabolites.items():
            details += f"  {coeff:+.1f} {met.id} ({met.name})\n"
        
        if hasattr(rxn, 'subsystem') and rxn.subsystem:
            details += f"\nSubsystem: {rxn.subsystem}"
        
        text.setPlainText(details)
        layout.addWidget(text)
        
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        layout.addWidget(close_btn)
        
        dialog.setLayout(layout)
        dialog.exec()

    def show_metabolite_details(self, metabolite_id):
        """Show details for a metabolite"""
        if not self.base_model or metabolite_id not in self.base_model.metabolites:
            return
        
        met = self.base_model.metabolites.get_by_id(metabolite_id)
        
        dialog = QDialog(self)
        dialog.setWindowTitle(f"Metabolite Details: {metabolite_id}")
        dialog.resize(500, 400)
        
        layout = QVBoxLayout()
        
        text = QPlainTextEdit()
        text.setReadOnly(True)
        
        details = f"""
ID: {met.id}
Name: {met.name}
Formula: {met.formula}
Charge: {met.charge}
Compartment: {met.compartment}

Participating Reactions: {len(met.reactions)}
"""
        for rxn in list(met.reactions)[:20]:
            details += f"  - {rxn.id}\n"
        
        text.setPlainText(details)
        layout.addWidget(text)
        
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        layout.addWidget(close_btn)
        
        dialog.setLayout(layout)
        dialog.exec()

    # -------- 6. UNDO/REDO --------

    def save_state_for_undo(self):
        """Save current model + UI state for undo."""
        if self.base_model is None:
            return

        state = {
            "knockout_genes": self.knockout_genes.copy(),
            "overexpression_reactions": self.overexpression_reactions.copy(),
            "reaction_bound_overrides": self.reaction_bound_overrides.copy(),
            "timestamp": datetime.now(),
            # UI state
            "active_tab_index": self.tabs.currentIndex() if hasattr(self, "tabs") else 0,
            "analysis_type_index": self.analysis_type.currentIndex() if hasattr(self, "analysis_type") else 0,
            "topn_value": self.topn_spin.value() if hasattr(self, "topn_spin") else 20,
        }

        self.undo_stack.append(state)
        if len(self.undo_stack) > self.max_undo_steps:
            self.undo_stack.pop(0)

        self.redo_stack.clear()

    def _capture_current_state(self) -> dict:
        """Capture current model + UI state for undo/redo stack."""
        return {
            "knockout_genes": self.knockout_genes.copy(),
            "overexpression_reactions": self.overexpression_reactions.copy(),
            "reaction_bound_overrides": self.reaction_bound_overrides.copy(),
            "active_tab_index": self.tabs.currentIndex() if hasattr(self, "tabs") else 0,
            "analysis_type_index": self.analysis_type.currentIndex() if hasattr(self, "analysis_type") else 0,
            "topn_value": self.topn_spin.value() if hasattr(self, "topn_spin") else 20,
        }

    def _restore_state(self, state: dict):
        """Restore model + UI state from a saved state dict."""
        self.knockout_genes = state["knockout_genes"]
        self.overexpression_reactions = state["overexpression_reactions"]
        self.reaction_bound_overrides = state["reaction_bound_overrides"]
        self.model_dirty = True

        # Restore UI state if available
        if hasattr(self, "tabs") and "active_tab_index" in state:
            idx = state["active_tab_index"]
            if 0 <= idx < self.tabs.count():
                self.tabs.setCurrentIndex(idx)
        if hasattr(self, "analysis_type") and "analysis_type_index" in state:
            idx = state["analysis_type_index"]
            if 0 <= idx < self.analysis_type.count():
                self.analysis_type.setCurrentIndex(idx)
        if hasattr(self, "topn_spin") and "topn_value" in state:
            self.topn_spin.setValue(state["topn_value"])

    def undo(self):
        """Undo last change (restores model + UI state)."""
        if not self.undo_stack:
            QMessageBox.information(self, "Undo", "Nothing to undo.")
            return

        self.redo_stack.append(self._capture_current_state())

        state = self.undo_stack.pop()
        self._restore_state(state)
        self.statusBar().showMessage("Undo completed.")

    def redo(self):
        """Redo last undone change (restores model + UI state)."""
        if not self.redo_stack:
            QMessageBox.information(self, "Redo", "Nothing to redo.")
            return

        self.undo_stack.append(self._capture_current_state())

        state = self.redo_stack.pop()
        self._restore_state(state)
        self.statusBar().showMessage("Redo completed.")

    # -------- 7. KEYBOARD SHORTCUTS --------

    def show_shortcuts(self):
        """Show keyboard shortcuts help"""
        dialog = QDialog(self)
        dialog.setWindowTitle("Keyboard Shortcuts")
        dialog.resize(500, 400)
        
        layout = QVBoxLayout()
        
        text = QPlainTextEdit()
        text.setReadOnly(True)
        
        shortcuts = """
KEYBOARD SHORTCUTS
==================

File Operations:
  Ctrl+O        Open SBML file
  Ctrl+S        Save as (SBML)
  Ctrl+Q        Exit application

Analysis:
  Ctrl+R        Run analysis
  Ctrl+Z        Undo
  Ctrl+Y        Redo
  Ctrl+F        Find/Search

Tools:
  Ctrl+B        Batch analysis
  Ctrl+E        Export results
  Ctrl+H        Show this help
"""
        
        text.setPlainText(shortcuts)
        layout.addWidget(text)
        
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        layout.addWidget(close_btn)
        
        dialog.setLayout(layout)
        dialog.exec()

    # ==================== DRAG & DROP ====================

    def dragEnterEvent(self, event):
        """Accept drag events for SBML files."""
        if event.mimeData().hasUrls():
            for url in event.mimeData().urls():
                fp = url.toLocalFile()
                if fp.endswith(('.xml', '.sbml', '.json')):
                    event.acceptProposedAction()
                    return

    def dropEvent(self, event):
        """Handle dropping SBML files onto the window."""
        for url in event.mimeData().urls():
            fp = url.toLocalFile()
            if fp.endswith(('.xml', '.sbml', '.json')):
                logger.info(f"File dropped: {fp}")
                self.open_sbml_from_path(fp)
                return

    # ==================== THEME TOGGLE ====================

    def toggle_theme(self):
        """Toggle between dark and light theme."""
        self._dark_theme = not self._dark_theme
        if self._dark_theme:
            self.setStyleSheet(DARK_STYLESHEET)
            logger.info("Theme switched to dark")
        else:
            self.setStyleSheet("")
            logger.info("Theme switched to light")
        self.statusBar().showMessage(f"Theme: {'Dark' if self._dark_theme else 'Light (System Default)'}")

    # ==================== GENE EXPRESSION INTEGRATION ====================

    def load_plugins(self):
        """Load plugins from ~/.metabodesk/plugins/ directory."""
        self._plugin_dir.mkdir(parents=True, exist_ok=True)

        plugin_files = list(self._plugin_dir.glob("*.py"))
        if not plugin_files:
            QMessageBox.information(
                self, "Plugins",
                f"No plugins found in:\n{self._plugin_dir}\n\n"
                "To create a plugin, place a .py file with a register(main_window) function\n"
                "in the plugins directory.\n\n"
                "Example plugin:\n"
                "  def register(mw):\n"
                "      from PySide6.QtGui import QAction\n"
                "      act = QAction('My Plugin', mw)\n"
                "      act.triggered.connect(lambda: print('Hello from plugin!'))\n"
                "      mw.menuBar().addAction(act)")
            return

        loaded = 0
        errors = []
        for pf in plugin_files:
            try:
                spec = importlib.util.spec_from_file_location(pf.stem, str(pf))
                if spec and spec.loader:
                    mod = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(mod)
                    if hasattr(mod, 'register'):
                        mod.register(self)
                        self._plugins.append(mod)
                        loaded += 1
                        logger.info(f"Plugin loaded: {pf.name}")
                    else:
                        errors.append(f"{pf.name}: missing register() function")
            except Exception as e:
                errors.append(f"{pf.name}: {e}")
                logger.error(f"Plugin load error: {pf.name}: {e}")

        msg = f"Loaded {loaded} plugin(s) from {self._plugin_dir}"
        if errors:
            msg += "\n\nErrors:\n" + "\n".join(errors)
        QMessageBox.information(self, "Plugins", msg)

    # ==================== INTERACTIVE TUTORIAL WIZARD ====================

    def show_tutorial_wizard(self) -> None:
        """Launch the Interactive Tutorial Wizard.

        A multi-page guided walkthrough teaching MetaboDesk features:
        Welcome → Loading Models → Medium Editor → Running Analyses →
        Advanced Analyses → Visualisation → Export → Finish.
        Each page contains rich-text explanations and contextual tips.
        """
        from PySide6.QtWidgets import QWizard, QWizardPage, QTextBrowser

        # ---------- Tutorial page definitions ----------
        pages: list[tuple[str, str]] = [
            (
                "🎓 Welcome to MetaboDesk",
                "<h2>Welcome!</h2>"
                "<p>This interactive tutorial will walk you through the key "
                "features of <b>MetaboDesk</b> — a desktop application for "
                "constraint-based metabolic modelling.</p>"
                "<p>You will learn how to:</p>"
                "<ol>"
                "<li>Load and inspect SBML metabolic models</li>"
                "<li>Edit the growth medium</li>"
                "<li>Run Flux Balance Analysis (FBA) and related analyses</li>"
                "<li>Use advanced analysis tools</li>"
                "<li>Visualise and export your results</li>"
                "</ol>"
                "<p>Click <b>Next</b> to begin, or <b>Cancel</b> to exit.</p>"
            ),
            (
                "📂 Step 1 — Loading a Model",
                "<h2>Loading an SBML Model</h2>"
                "<p>Go to <b>File → Open SBML...</b> (or press <code>Ctrl+O</code>) "
                "and select an <code>.xml</code> or <code>.sbml</code> file.</p>"
                "<p>Once loaded you will see:</p>"
                "<ul>"
                "<li><b>Reactions tab</b> — all reactions with bounds, GPR rules, "
                "and subsystems</li>"
                "<li><b>Metabolites tab</b> — metabolites with formulas and "
                "compartments</li>"
                "<li><b>Genes tab</b> — gene list (for gene-deletion studies)</li>"
                "</ul>"
                "<p><b>💡 Tip:</b> You can also drag-and-drop an SBML file onto "
                "the main window to open it instantly.</p>"
            ),
            (
                "🧪 Step 2 — Medium Editor",
                "<h2>Editing the Growth Medium</h2>"
                "<p>The <b>Medium</b> tab lists exchange reactions that represent "
                "nutrient uptake.</p>"
                "<p>To change the medium:</p>"
                "<ul>"
                "<li>Double-click a <b>Lower Bound</b> cell to allow/block uptake</li>"
                "<li>Use <b>Medium → Close uptakes</b> to block all imports</li>"
                "<li>Use <b>Medium → Reset medium</b> to restore original bounds</li>"
                "</ul>"
                "<p>You can also <b>export</b> and <b>import</b> medium bounds "
                "as CSV files via the Medium menu.</p>"
                "<p><b>💡 Tip:</b> Negative lower bounds on exchange reactions "
                "mean the model can <i>import</i> that metabolite.</p>"
            ),
            (
                "▶️ Step 3 — Running Analyses",
                "<h2>Flux Balance Analysis &amp; More</h2>"
                "<p>Select an analysis type from the <b>Analysis panel</b> on the "
                "left side, then click <b>Run</b> (or press <code>Ctrl+R</code>).</p>"
                "<p>Available analyses:</p>"
                "<table cellpadding='4'>"
                "<tr><td>📊</td><td><b>FBA</b></td><td>Standard flux balance analysis</td></tr>"
                "<tr><td>🎯</td><td><b>pFBA</b></td><td>Parsimonious FBA (minimum total flux)</td></tr>"
                "<tr><td>📈</td><td><b>FVA</b></td><td>Flux variability analysis</td></tr>"
                "<tr><td>🧬</td><td><b>Gene Deletion</b></td><td>Single/double gene/reaction deletion</td></tr>"
                "<tr><td>📉</td><td><b>Robustness</b></td><td>Objective vs. reaction bound sweep</td></tr>"
                "<tr><td>🔬</td><td><b>Production Envelope</b></td><td>Phenotype phase plane</td></tr>"
                "<tr><td>🎲</td><td><b>Flux Sampling</b></td><td>Monte-Carlo sampling of flux space</td></tr>"
                "</table>"
                "<p><b>💡 Tip:</b> Results appear in the right panel with both a "
                "chart and a data table. You can export them as CSV, Excel, or JSON.</p>"
            ),
            (
                "🔬 Step 4 — Advanced Analyses",
                "<h2>Advanced Analysis Tools</h2>"
                "<p>Under <b>Run → Advanced Analysis</b> you'll find:</p>"
                "<table cellpadding='4'>"
                "<tr><td>🧬</td><td><b>Gene Expression Integration</b></td>"
                "<td>E-Flux / GIMME / iMAT context-specific modelling</td></tr>"
                "<tr><td>⚡</td><td><b>OptKnock</b></td>"
                "<td>Combinatorial knockout search for metabolite production</td></tr>"
                "<tr><td>🔄</td><td><b>Dynamic FBA</b></td>"
                "<td>Time-course batch-culture simulation</td></tr>"
                "<tr><td>🌡️</td><td><b>TMFA</b></td>"
                "<td>Thermodynamically constrained FBA</td></tr>"
                "<tr><td>🔗</td><td><b>Flux Coupling</b></td>"
                "<td>Identify coupled / blocked reactions</td></tr>"
                "<tr><td>🧩</td><td><b>Gap-Filling</b></td>"
                "<td>Find missing reactions to rescue growth</td></tr>"
                "<tr><td>🛡️</td><td><b>GECKO</b></td>"
                "<td>Enzyme-constrained FBA with kcat/MW data</td></tr>"
                "<tr><td>📡</td><td><b>FSEOF</b></td>"
                "<td>Identify overexpression / attenuation targets</td></tr>"
                "<tr><td>✂️</td><td><b>Minimal Cut Sets</b></td>"
                "<td>Find minimal knockouts to block a target flux</td></tr>"
                "<tr><td>📏</td><td><b>Metabolic Distance</b></td>"
                "<td>BFS shortest-path through the metabolic network</td></tr>"
                "<tr><td>⚖️</td><td><b>Pareto Frontier</b></td>"
                "<td>Multi-objective optimisation trade-off curve</td></tr>"
                "</table>"
            ),
            (
                "🗺️ Step 5 — Visualisation & Network",
                "<h2>Visualisation Tools</h2>"
                "<p>MetaboDesk provides multiple visualisation options:</p>"
                "<ul>"
                "<li><b>Results chart</b> — automatic charts for every analysis "
                "(bar, line, scatter, heatmap)</li>"
                "<li><b>Network view</b> — interactive metabolic network graph "
                "with search, filtering, and layout options "
                "(spring / circular / Kamada-Kawai)</li>"
                "<li><b>Phenotype Phase Plane</b> — 2-D heatmap of objective "
                "across two exchange-reaction sweeps</li>"
                "<li><b>Escher Map Viewer</b> — load Escher JSON maps for "
                "pathway-level visualisation</li>"
                "</ul>"
                "<p><b>💡 Tip:</b> Right-click the chart area for a context menu "
                "with save/copy options. You can also zoom and pan.</p>"
            ),
            (
                "💾 Step 6 — Export & Reproducibility",
                "<h2>Exporting Your Work</h2>"
                "<p>MetaboDesk supports many export formats:</p>"
                "<ul>"
                "<li><b>File → Export results (CSV/Excel/JSON)</b></li>"
                "<li><b>File → Export All</b> — scenario + results + chart</li>"
                "<li><b>File → Export Jupyter Notebook</b> — reproducible Python "
                "notebook of your analysis pipeline</li>"
                "<li><b>File → Export/Import scenario</b> — save and reload your "
                "complete medium + knockout + overexpression settings</li>"
                "</ul>"
                "<p><b>Model operations:</b></p>"
                "<ul>"
                "<li><b>File → Save As (SBML)</b> — save modified model</li>"
                "<li><b>Model → Validate SBML</b> — check model consistency</li>"
                "<li><b>Model → Generate HTML Report</b> — full model report</li>"
                "</ul>"
            ),
            (
                "🎉 Tutorial Complete!",
                "<h2>You're All Set!</h2>"
                "<p>You've completed the MetaboDesk interactive tutorial. "
                "Here's a quick-reference of keyboard shortcuts:</p>"
                "<table cellpadding='4'>"
                "<tr><td><code>Ctrl+O</code></td><td>Open SBML model</td></tr>"
                "<tr><td><code>Ctrl+R</code></td><td>Run analysis</td></tr>"
                "<tr><td><code>Ctrl+Z / Ctrl+Y</code></td><td>Undo / Redo</td></tr>"
                "<tr><td><code>Ctrl+F</code></td><td>Search &amp; Filter</td></tr>"
                "<tr><td><code>Ctrl+B</code></td><td>Batch Analysis</td></tr>"
                "<tr><td><code>Ctrl+T</code></td><td>Toggle Theme</td></tr>"
                "<tr><td><code>Ctrl+Q</code></td><td>Quit</td></tr>"
                "</table>"
                "<p>You can re-open this tutorial any time from "
                "<b>Help → 🎓 Interactive Tutorial</b>.</p>"
                "<p><b>Happy modelling! 🧬</b></p>"
            ),
        ]

        wizard = QWizard(self)
        wizard.setWindowTitle("MetaboDesk — Interactive Tutorial")
        wizard.setMinimumSize(700, 520)
        wizard.setWizardStyle(QWizard.ModernStyle)
        wizard.setOption(QWizard.NoBackButtonOnStartPage, True)

        for title, html_content in pages:
            page = QWizardPage()
            page.setTitle(title)
            page_layout = QVBoxLayout(page)

            browser = QTextBrowser()
            browser.setOpenExternalLinks(True)
            browser.setHtml(html_content)
            browser.setStyleSheet(
                "QTextBrowser { border: none; padding: 8px; font-size: 13px; }"
            )
            page_layout.addWidget(browser)

            wizard.addPage(page)

        wizard.exec()

    # ==================== BiGG Models Database Browser ====================

    def show_bigg_browser(self) -> None:
        """Open an interactive BiGG Models database browser dialog.

        Queries the public BiGG API (http://bigg.ucsd.edu/api/v2/) to
        list universal reactions and metabolites.  The user can search,
        inspect details, and — when a model is loaded — add reactions
        from BiGG directly into the current model.
        """
        import urllib.request
        import json as _json

        dialog = QDialog(self)
        dialog.setWindowTitle("BiGG Models Database Browser")
        dialog.resize(1050, 700)
        layout = QVBoxLayout(dialog)

        # ---- top bar ----
        top = QHBoxLayout()
        mode_combo = QComboBox()
        mode_combo.addItems(["Models", "Universal Reactions", "Universal Metabolites"])
        top.addWidget(QLabel("Browse:"))
        top.addWidget(mode_combo)
        search_edit = QLineEdit()
        search_edit.setPlaceholderText("Search by ID or name…")
        top.addWidget(search_edit, stretch=1)
        search_btn = QPushButton("Search")
        top.addWidget(search_btn)
        top.addStretch()
        layout.addLayout(top)

        # ---- results table ----
        table = QTableWidget()
        table.setColumnCount(4)
        table.setHorizontalHeaderLabels(["ID", "Name", "Organism / Formula", "Reaction count"])
        table.horizontalHeader().setStretchLastSection(True)
        table.setSelectionBehavior(QAbstractItemView.SelectRows)
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        layout.addWidget(table, stretch=1)

        # ---- detail area ----
        detail_lbl = QLabel("Select an item above to see details.")
        detail_lbl.setWordWrap(True)
        detail_lbl.setMinimumHeight(80)
        detail_lbl.setStyleSheet("background:#f4f4f4; padding:8px; border:1px solid #ccc;")
        layout.addWidget(detail_lbl)

        # ---- bottom buttons ----
        bottom = QHBoxLayout()
        download_btn = QPushButton("Download Model (SBML)…")
        download_btn.setEnabled(False)
        add_rxn_btn = QPushButton("Add Reaction to Current Model")
        add_rxn_btn.setEnabled(False)
        open_web_btn = QPushButton("Open in Browser")
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        bottom.addWidget(download_btn)
        bottom.addWidget(add_rxn_btn)
        bottom.addWidget(open_web_btn)
        bottom.addStretch()
        bottom.addWidget(close_btn)
        layout.addLayout(bottom)

        # ---- state ----
        _cache: dict[str, object] = {}

        def _api_get(endpoint: str) -> object:
            """Fetch JSON from BiGG API with caching."""
            if endpoint in _cache:
                return _cache[endpoint]
            url = f"http://bigg.ucsd.edu/api/v2/{endpoint}"
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=12) as resp:
                data = _json.loads(resp.read().decode("utf-8"))
            _cache[endpoint] = data
            return data

        def _do_search():
            mode = mode_combo.currentText()
            query = search_edit.text().strip().lower()
            table.setRowCount(0)
            download_btn.setEnabled(False)
            add_rxn_btn.setEnabled(False)
            detail_lbl.setText("Loading…")
            dialog.repaint()

            try:
                if mode == "Models":
                    data = _api_get("models")
                    results = data.get("results", [])
                    table.setHorizontalHeaderLabels(["Model ID", "Organism", "Metabolite count", "Reaction count"])
                    if query:
                        results = [r for r in results
                                   if query in (r.get("bigg_id") or "").lower()
                                   or query in (r.get("organism") or "").lower()]
                    for r in results:
                        row = table.rowCount()
                        table.insertRow(row)
                        table.setItem(row, 0, QTableWidgetItem(r.get("bigg_id") or ""))
                        table.setItem(row, 1, QTableWidgetItem(r.get("organism") or ""))
                        table.setItem(row, 2, QTableWidgetItem(str(r.get("metabolite_count", ""))))
                        table.setItem(row, 3, QTableWidgetItem(str(r.get("reaction_count", ""))))

                elif mode == "Universal Reactions":
                    data = _api_get("universal/reactions")
                    results = data.get("results", [])
                    table.setHorizontalHeaderLabels(["Reaction ID", "Name", "Model count", ""])
                    if query:
                        results = [r for r in results
                                   if query in (r.get("bigg_id") or "").lower()
                                   or query in (r.get("name") or "").lower()]
                    for r in results[:500]:
                        row = table.rowCount()
                        table.insertRow(row)
                        table.setItem(row, 0, QTableWidgetItem(r.get("bigg_id") or ""))
                        table.setItem(row, 1, QTableWidgetItem(r.get("name") or ""))
                        table.setItem(row, 2, QTableWidgetItem(str(r.get("model_bigg_id", r.get("model_count", "")))))
                        table.setItem(row, 3, QTableWidgetItem(""))

                else:  # Universal Metabolites
                    data = _api_get("universal/metabolites")
                    results = data.get("results", [])
                    table.setHorizontalHeaderLabels(["Metabolite ID", "Name", "Model count", ""])
                    if query:
                        results = [r for r in results
                                   if query in (r.get("bigg_id") or "").lower()
                                   or query in (r.get("name") or "").lower()]
                    for r in results[:500]:
                        row = table.rowCount()
                        table.insertRow(row)
                        table.setItem(row, 0, QTableWidgetItem(r.get("bigg_id", "")))
                        table.setItem(row, 1, QTableWidgetItem(r.get("name", "")))
                        table.setItem(row, 2, QTableWidgetItem(str(r.get("model_bigg_id", r.get("model_count", "")))))
                        table.setItem(row, 3, QTableWidgetItem(""))

                detail_lbl.setText(f"{table.rowCount()} results shown.")
            except Exception as e:
                detail_lbl.setText(f"API error: {e}")
                logger.warning("BiGG API error: %s", e)

        def _on_row_selected():
            rows = table.selectionModel().selectedRows()
            if not rows:
                return
            idx = rows[0].row()
            mode = mode_combo.currentText()
            bigg_id = table.item(idx, 0).text()
            download_btn.setEnabled(mode == "Models")
            add_rxn_btn.setEnabled(mode == "Universal Reactions" and self.base_model is not None)

            try:
                if mode == "Models":
                    info = _api_get(f"models/{bigg_id}")
                    detail_lbl.setText(
                        f"<b>{bigg_id}</b> — {info.get('organism', 'N/A')}<br>"
                        f"Genome: {info.get('genome_name', 'N/A')}<br>"
                        f"Reactions: {info.get('reaction_count', '?')} | "
                        f"Metabolites: {info.get('metabolite_count', '?')} | "
                        f"Genes: {info.get('gene_count', '?')}<br>"
                        f"Reference: {info.get('reference_id', 'N/A')} — {info.get('reference_type', '')}"
                    )
                elif mode == "Universal Reactions":
                    info = _api_get(f"universal/reactions/{bigg_id}")
                    db_links = info.get("database_links") or {}
                    link_strs = []
                    for db_name, entries in db_links.items():
                        if not isinstance(entries, list):
                            continue
                        for entry in entries[:2]:
                            if isinstance(entry, dict):
                                link_strs.append(f"{db_name}: {entry.get('id', '?')}")
                            else:
                                link_strs.append(f"{db_name}: {entry}")
                    models_in = [(m.get("bigg_id") or "") if isinstance(m, dict) else str(m)
                                 for m in (info.get("models_containing_reaction") or [])[:8]]
                    detail_lbl.setText(
                        f"<b>{bigg_id}</b> — {info.get('name') or 'N/A'}<br>"
                        f"Pseudoreaction: {info.get('pseudoreaction', False)}<br>"
                        f"Models: {', '.join(models_in)}{'…' if len(info.get('models_containing_reaction') or []) > 8 else ''}<br>"
                        f"DB links: {'; '.join(link_strs[:6])}"
                    )
                else:  # Metabolites
                    info = _api_get(f"universal/metabolites/{bigg_id}")
                    # formulae / charges may be list-of-dicts OR list-of-primitives
                    raw_formulae = info.get("formulae") or []
                    formulae: list[str] = []
                    for f in raw_formulae[:3]:
                        if isinstance(f, dict):
                            formulae.append(str(f.get("formula", "")))
                        else:
                            formulae.append(str(f))
                    raw_charges = info.get("charges") or []
                    charges: list[str] = []
                    for ch in raw_charges[:3]:
                        if isinstance(ch, dict):
                            charges.append(str(ch.get("charge", "")))
                        else:
                            charges.append(str(ch))
                    detail_lbl.setText(
                        f"<b>{bigg_id}</b> — {info.get('name') or 'N/A'}<br>"
                        f"Formula: {', '.join(formulae) or 'N/A'} | "
                        f"Charge: {', '.join(charges) or 'N/A'}<br>"
                        f"Found in {len(info.get('compartments_in_models', []))} model-compartment(s)"
                    )
            except Exception as e:
                detail_lbl.setText(f"Could not fetch details: {e}")

        def _download_model():
            rows = table.selectionModel().selectedRows()
            if not rows:
                return
            bigg_id = table.item(rows[0].row(), 0).text()
            save_path, _ = QFileDialog.getSaveFileName(
                dialog, "Save BiGG Model", f"{bigg_id}.xml", "SBML files (*.xml)")
            if not save_path:
                return
            try:
                url = f"http://bigg.ucsd.edu/static/models/{bigg_id}.xml"
                req = urllib.request.Request(url)
                detail_lbl.setText(f"Downloading {bigg_id}.xml …")
                dialog.repaint()
                with urllib.request.urlopen(req, timeout=60) as resp:
                    Path(save_path).write_bytes(resp.read())
                detail_lbl.setText(f"Downloaded to {save_path}")
                reply = QMessageBox.question(
                    dialog, "Load Model?",
                    f"Model saved to:\n{save_path}\n\nLoad it now?",
                    QMessageBox.Yes | QMessageBox.No)
                if reply == QMessageBox.Yes:
                    dialog.close()
                    self.open_sbml_from_path(save_path)
            except Exception as e:
                detail_lbl.setText(f"Download failed: {e}")

        def _add_reaction_to_model():
            """Add a universal reaction from BiGG to the current model."""
            rows = table.selectionModel().selectedRows()
            if not rows or self.base_model is None:
                return
            bigg_id = table.item(rows[0].row(), 0).text()
            try:
                info = _api_get(f"universal/reactions/{bigg_id}")
                existing_ids = {r.id for r in self.base_model.reactions}
                if bigg_id in existing_ids:
                    QMessageBox.information(dialog, "Exists", f"Reaction {bigg_id} already exists in the model.")
                    return

                rxn = cobra.Reaction(bigg_id)
                rxn.name = info.get("name", bigg_id)
                rxn.lower_bound = -1000.0 if not info.get("pseudoreaction", False) else 0.0
                rxn.upper_bound = 1000.0

                # Try to retrieve stoichiometry from a reference model in BiGG
                model_refs = info.get("models_containing_reaction", [])
                stoich_added = False
                for mref in model_refs[:3]:
                    mid = mref.get("bigg_id", "")
                    if not mid:
                        continue
                    try:
                        rxn_detail = _api_get(f"models/{mid}/reactions/{bigg_id}")
                        metabolites_info = rxn_detail.get("metabolites", [])
                        if not metabolites_info:
                            continue
                        mets: dict[cobra.Metabolite, float] = {}
                        for m in metabolites_info:
                            met_id = m.get("bigg_id", "") + "_" + m.get("compartment_bigg_id", "c")
                            coeff = m.get("stoichiometry", 0.0)
                            try:
                                met = self.base_model.metabolites.get_by_id(met_id)
                            except KeyError:
                                met = cobra.Metabolite(
                                    met_id,
                                    name=m.get("name", met_id),
                                    compartment=m.get("compartment_bigg_id", "c"),
                                )
                            mets[met] = coeff
                        if mets:
                            rxn.add_metabolites(mets)
                            stoich_added = True
                            break
                    except Exception:
                        continue

                self.base_model.add_reactions([rxn])
                self.model_dirty = True
                msg = f"Reaction {bigg_id} added to model"
                if not stoich_added:
                    msg += " (stub — no stoichiometry found, edit manually)"
                detail_lbl.setText(msg)
                QMessageBox.information(dialog, "Added", msg + ".")
            except Exception as e:
                detail_lbl.setText(f"Failed to add reaction: {e}")

        def _open_in_browser():
            rows = table.selectionModel().selectedRows()
            if not rows:
                return
            bigg_id = table.item(rows[0].row(), 0).text()
            mode = mode_combo.currentText()
            if mode == "Models":
                url = f"http://bigg.ucsd.edu/models/{bigg_id}"
            elif mode == "Universal Reactions":
                url = f"http://bigg.ucsd.edu/universal/reactions/{bigg_id}"
            else:
                url = f"http://bigg.ucsd.edu/universal/metabolites/{bigg_id}"
            QDesktopServices.openUrl(QUrl(url))

        search_btn.clicked.connect(_do_search)
        search_edit.returnPressed.connect(_do_search)
        table.selectionModel().selectionChanged.connect(lambda *_: _on_row_selected())
        download_btn.clicked.connect(_download_model)
        add_rxn_btn.clicked.connect(_add_reaction_to_model)
        open_web_btn.clicked.connect(_open_in_browser)

        dialog.exec()

    # ==================== Native Qt Interactive Metabolic Map ====================

    def show_interactive_metabolic_map(self) -> None:
        """Open a native Qt interactive metabolic pathway map using QGraphicsView.

        Provides a zoomable, pannable, interactive metabolic network
        visualisation without requiring Escher or QWebEngine.  Nodes are
        colour-coded by flux magnitude and can be clicked for details.
        """
        from PySide6.QtWidgets import (
            QGraphicsView, QGraphicsScene, QSplitter,
        )
        from PySide6.QtGui import QBrush, QPen, QColor, QFont, QPainter
        from PySide6.QtCore import QRectF, QLineF
        import networkx as nx

        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        model = self.base_model
        flux: dict[str, float] = {}
        if self.last_run and "baseline" in self.last_run and "flux" in self.last_run["baseline"]:
            flux = self.last_run["baseline"]["flux"]

        # ---- Dialog setup ----
        dialog = QDialog(self)
        dialog.setWindowTitle("Interactive Metabolic Map")
        dialog.resize(1300, 850)
        main_layout = QVBoxLayout(dialog)

        # ---- Toolbar ----
        toolbar = QHBoxLayout()
        toolbar.addWidget(QLabel("Search:"))
        search_box = QLineEdit()
        search_box.setPlaceholderText("reaction / metabolite ID…")
        toolbar.addWidget(search_box, stretch=1)

        toolbar.addWidget(QLabel("Subsystem:"))
        sub_combo = QComboBox()
        subsystems = sorted({rxn.subsystem for rxn in model.reactions if rxn.subsystem})
        sub_combo.addItem("All")
        sub_combo.addItems(subsystems)
        toolbar.addWidget(sub_combo)

        toolbar.addWidget(QLabel("Depth:"))
        depth_spin = QSpinBox()
        depth_spin.setRange(1, 4)
        depth_spin.setValue(2)
        toolbar.addWidget(depth_spin)

        render_btn = QPushButton("Render")
        toolbar.addWidget(render_btn)

        fit_btn = QPushButton("Fit View")
        toolbar.addWidget(fit_btn)

        export_btn = QPushButton("Export PNG…")
        toolbar.addWidget(export_btn)

        toolbar.addStretch()
        main_layout.addLayout(toolbar)

        # ---- Splitter: scene + detail panel ----
        splitter = QSplitter(Qt.Horizontal)

        # Custom zoomable graphics view
        class _ZoomableView(QGraphicsView):
            """QGraphicsView subclass with mouse-wheel zoom and hand-drag pan."""

            def __init__(self, s):
                super().__init__(s)
                self.setRenderHints(QPainter.Antialiasing | QPainter.TextAntialiasing)
                self.setDragMode(QGraphicsView.ScrollHandDrag)
                self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
                self._zoom = 1.0

            def wheelEvent(self, event):
                factor = 1.15 if event.angleDelta().y() > 0 else 1 / 1.15
                new_zoom = self._zoom * factor
                if 0.05 < new_zoom < 20:
                    self._zoom = new_zoom
                    self.scale(factor, factor)

        gscene = QGraphicsScene()
        view = _ZoomableView(gscene)
        splitter.addWidget(view)

        # Detail panel
        detail_panel = QLabel(
            "Click a node for details.\n\n"
            "🔵 = Reaction    🟩 = Metabolite\n"
            "Colour intensity = flux magnitude\n\n"
            "Scroll to zoom, drag to pan."
        )
        detail_panel.setWordWrap(True)
        detail_panel.setMinimumWidth(250)
        detail_panel.setMaximumWidth(350)
        detail_panel.setAlignment(Qt.AlignTop)
        detail_panel.setStyleSheet("background:#fafafa; padding:10px; border:1px solid #ddd;")
        splitter.addWidget(detail_panel)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)

        main_layout.addWidget(splitter, stretch=1)

        def _build_graph():
            """Build and render the metabolic graph on the QGraphicsScene."""
            gscene.clear()

            search_text = search_box.text().strip().lower()
            subsystem_filter = sub_combo.currentText()
            depth = depth_spin.value()

            # ---- Filter reactions ----
            filtered_rxns: list[cobra.Reaction] = []
            for rxn in model.reactions:
                if subsystem_filter != "All" and rxn.subsystem != subsystem_filter:
                    continue
                if search_text:
                    if (search_text not in rxn.id.lower()
                            and search_text not in rxn.name.lower()
                            and not any(search_text in m.id.lower() or search_text in m.name.lower()
                                        for m in rxn.metabolites)):
                        continue
                filtered_rxns.append(rxn)

            # BFS expansion
            if search_text and depth > 1:
                rxn_ids_set = {r.id for r in filtered_rxns}
                for _ in range(depth - 1):
                    neighbours: set[str] = set()
                    for rxn in filtered_rxns:
                        for met in rxn.metabolites:
                            for r2 in met.reactions:
                                if r2.id not in rxn_ids_set:
                                    neighbours.add(r2.id)
                    for rxn in model.reactions:
                        if rxn.id in neighbours:
                            filtered_rxns.append(rxn)
                            rxn_ids_set.add(rxn.id)

            # Limit for performance
            if len(filtered_rxns) > 300:
                filtered_rxns = filtered_rxns[:300]

            if not filtered_rxns:
                gscene.addText("No reactions match the filter.").setDefaultTextColor(QColor("#999"))
                return

            # Collect metabolites
            met_set: set[str] = set()
            for rxn in filtered_rxns:
                for m in rxn.metabolites:
                    met_set.add(m.id)

            # Layout via networkx
            G = nx.Graph()
            for rxn in filtered_rxns:
                G.add_node(rxn.id, ntype="rxn")
                for m in rxn.metabolites:
                    G.add_node(m.id, ntype="met")
                    G.add_edge(rxn.id, m.id)

            try:
                pos = nx.kamada_kawai_layout(G, scale=400)
            except Exception:
                pos = nx.spring_layout(G, k=2.5, iterations=60, scale=400)

            max_flux = max((abs(flux.get(r.id, 0)) for r in filtered_rxns), default=1.0) or 1.0

            # ---- Draw edges ----
            edge_pen = QPen(QColor(180, 180, 180, 120), 1.0)
            for rxn in filtered_rxns:
                rx, ry = pos.get(rxn.id, (0, 0))
                for m in rxn.metabolites:
                    mx, my = pos.get(m.id, (0, 0))
                    line = gscene.addLine(QLineF(rx, ry, mx, my), edge_pen)
                    line.setZValue(0)

            # ---- Draw metabolite nodes (squares) ----
            met_font = QFont("Arial", 6)
            for mid in met_set:
                x, y = pos.get(mid, (0, 0))
                size = 12
                rect = gscene.addRect(
                    QRectF(x - size / 2, y - size / 2, size, size),
                    QPen(QColor(34, 139, 34), 1.0),
                    QBrush(QColor(144, 238, 144, 200)),
                )
                rect.setZValue(2)
                rect.setToolTip(mid)
                rect.setData(0, mid)
                rect.setData(1, "met")
                label = gscene.addText(mid[:12], met_font)
                label.setPos(x + 8, y - 6)
                label.setDefaultTextColor(QColor(60, 60, 60))
                label.setZValue(3)

            # ---- Draw reaction nodes (circles) ----
            rxn_font = QFont("Arial", 7, QFont.Bold)
            for rxn in filtered_rxns:
                x, y = pos.get(rxn.id, (0, 0))
                f_val = abs(flux.get(rxn.id, 0.0))
                intensity = min(int(55 + 200 * (f_val / max_flux)), 255)
                colour = QColor(intensity, 80, 255 - intensity, 220)
                radius = 9 + min(f_val / max_flux * 8, 8)
                ellipse = gscene.addEllipse(
                    QRectF(x - radius, y - radius, radius * 2, radius * 2),
                    QPen(QColor(30, 30, 120), 1.5),
                    QBrush(colour),
                )
                ellipse.setZValue(2)
                ellipse.setToolTip(f"{rxn.id}: {rxn.name}\nFlux: {flux.get(rxn.id, 0.0):.4g}")
                ellipse.setData(0, rxn.id)
                ellipse.setData(1, "rxn")
                label = gscene.addText(rxn.id[:14], rxn_font)
                label.setPos(x + radius + 2, y - 8)
                label.setDefaultTextColor(QColor(20, 20, 100))
                label.setZValue(3)

            view.fitInView(gscene.sceneRect().adjusted(-50, -50, 50, 50), Qt.KeepAspectRatio)

        def _on_scene_click(event):
            """Handle click on a node to show details."""
            item = gscene.itemAt(event.scenePos(), view.transform())
            if item is None:
                return
            node_id = item.data(0)
            node_type = item.data(1)
            if node_id is None:
                return

            if node_type == "rxn":
                try:
                    rxn = model.reactions.get_by_id(node_id)
                    f_val = flux.get(node_id, 0.0)
                    gpr = getattr(rxn, "gene_reaction_rule", "") or ""
                    detail_panel.setText(
                        f"<b>Reaction: {node_id}</b><br>"
                        f"Name: {rxn.name}<br>"
                        f"Subsystem: {rxn.subsystem or 'N/A'}<br><br>"
                        f"<b>Equation:</b><br>{rxn.reaction}<br><br>"
                        f"Bounds: [{rxn.lower_bound}, {rxn.upper_bound}]<br>"
                        f"Flux: {f_val:.6g}<br>"
                        f"GPR: {gpr or 'None'}<br>"
                        f"Genes: {', '.join(g.id for g in rxn.genes) or 'None'}"
                    )
                except KeyError:
                    detail_panel.setText(f"Reaction {node_id} not found.")
            elif node_type == "met":
                try:
                    met = model.metabolites.get_by_id(node_id)
                    rxns_using = [r.id for r in met.reactions]
                    detail_panel.setText(
                        f"<b>Metabolite: {node_id}</b><br>"
                        f"Name: {met.name}<br>"
                        f"Formula: {met.formula or 'N/A'}<br>"
                        f"Compartment: {met.compartment}<br>"
                        f"Charge: {met.charge if met.charge is not None else 'N/A'}<br><br>"
                        f"<b>Reactions ({len(rxns_using)}):</b><br>"
                        f"{', '.join(rxns_using[:20])}{'…' if len(rxns_using) > 20 else ''}"
                    )
                except KeyError:
                    detail_panel.setText(f"Metabolite {node_id} not found.")

        gscene.mousePressEvent = _on_scene_click

        def _fit_view():
            view.fitInView(gscene.sceneRect().adjusted(-30, -30, 30, 30), Qt.KeepAspectRatio)

        def _export_png():
            path, _ = QFileDialog.getSaveFileName(dialog, "Export Map", "", "PNG (*.png);;SVG (*.svg)")
            if not path:
                return
            try:
                from PySide6.QtGui import QImage
                rect = gscene.sceneRect()
                img = QImage(int(rect.width() * 2), int(rect.height() * 2), QImage.Format_ARGB32)
                img.fill(QColor(255, 255, 255))
                p = QPainter(img)
                p.setRenderHint(QPainter.Antialiasing)
                gscene.render(p, source=rect)
                p.end()
                img.save(path)
                QMessageBox.information(dialog, "Export", f"Map exported to:\n{path}")
            except Exception as e:
                QMessageBox.critical(dialog, "Export Failed", str(e))

        render_btn.clicked.connect(_build_graph)
        fit_btn.clicked.connect(_fit_view)
        export_btn.clicked.connect(_export_png)

        # Initial render
        _build_graph()

        dialog.exec()
