"""Medium table and presets mixin for MetaboDesk.

Manages the exchange-reaction table that represents the growth medium,
including filtering, inline editing of lower/upper bounds, medium
presets (Minimal M9, Rich LB), and CSV import/export of medium bounds.
"""

import csv
import logging

import cobra
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QFileDialog, QMessageBox, QTableWidgetItem

from metabodesk_core.utils import is_exchange_reaction

logger = logging.getLogger("MetaboDesk")

class MediumMixin:
    """Mixin providing medium functionality."""

    def filter_medium_table(self):
        text = self.medium_search.text().strip().lower()
        for row in range(self.medium_table.rowCount()):
            rxn_id = (self.medium_table.item(row, 0).text() if self.medium_table.item(row, 0) else "").lower()
            name = (self.medium_table.item(row, 1).text() if self.medium_table.item(row, 1) else "").lower()
            hide = bool(text) and (text not in rxn_id) and (text not in name)
            self.medium_table.setRowHidden(row, hide)

    def export_medium_bounds_csv(self):
        if self.medium_table.rowCount() == 0:
            QMessageBox.warning(self, "Medium", "No medium rows to export.")
            return
        file_path, _ = QFileDialog.getSaveFileName(self, "Export Medium bounds (CSV)", str(Path.home() / "medium_bounds.csv"), "CSV files (*.csv);")
        if not file_path:
            return
        try:
            with open(file_path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["reaction", "lb", "ub"]) 
                for i in range(self.medium_table.rowCount()):
                    rid = self.medium_table.item(i, 0).text().strip()
                    lb = self.medium_table.item(i, 2).text().strip()
                    ub = self.medium_table.item(i, 3).text().strip()
                    w.writerow([rid, lb, ub])
            self.statusBar().showMessage(f"Medium bounds exported: {file_path}")
        except Exception as e:
            self._show_error("Export failed", "Failed to export medium bounds.", e)

    def import_medium_bounds_csv(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Import Medium bounds (CSV)", str(Path.home()), "CSV files (*.csv);")
        if not file_path:
            return
        try:
            rxn_to_row = {self.medium_table.item(i, 0).text().strip(): i for i in range(self.medium_table.rowCount())}
            with open(file_path, "r", encoding="utf-8") as f:
                reader = csv.reader(f)
                for rowvals in reader:
                    if not rowvals or rowvals[0] == "reaction":
                        continue
                    rid, lb, ub = rowvals[:3]
                    row = rxn_to_row.get(rid)
                    if row is None:
                        continue
                    self.medium_table.setItem(row, 2, QTableWidgetItem(str(lb)))
                    self.medium_table.setItem(row, 3, QTableWidgetItem(str(ub)))
            self.filter_medium_table()
            self.statusBar().showMessage(f"Medium bounds imported: {file_path}")
        except Exception as e:
            self._show_error("Import failed", "Failed to import medium bounds.", e)

    def populate_medium_table(self):
        self.medium_table.blockSignals(True)
        self.medium_table.setSortingEnabled(False)
        n = len(self.exchange_rxns)
        self.medium_table.setRowCount(n)

        readonly = Qt.ItemIsSelectable | Qt.ItemIsEnabled
        for row, rxn in enumerate(self.exchange_rxns):
            item_id = QTableWidgetItem(rxn.id)
            item_id.setFlags(readonly)

            item_name = QTableWidgetItem(rxn.name or "")
            item_name.setFlags(readonly)

            self.medium_table.setItem(row, 0, item_id)
            self.medium_table.setItem(row, 1, item_name)
            self.medium_table.setItem(row, 2, QTableWidgetItem(str(float(rxn.lower_bound))))
            self.medium_table.setItem(row, 3, QTableWidgetItem(str(float(rxn.upper_bound))))

        self.medium_table.setSortingEnabled(True)
        self.medium_table.blockSignals(False)
        self.filter_medium_table()

    def apply_medium_table_bounds_to_model(self, model: cobra.Model):
        rxn_by_id = {r.id: r for r in model.reactions}
        for i in range(self.medium_table.rowCount()):
            rxn_id = self.medium_table.item(i, 0).text().strip()
            lb_txt = self.medium_table.item(i, 2).text().strip()
            ub_txt = self.medium_table.item(i, 3).text().strip()
            lb = float(lb_txt)
            ub = float(ub_txt)
            if lb > ub:
                raise ValueError(f"Lower bound > upper bound in Medium row {i+1} for reaction {rxn_id}.")
            if rxn_id in rxn_by_id:
                rxn_by_id[rxn_id].lower_bound = lb
                rxn_by_id[rxn_id].upper_bound = ub

    def apply_medium_preset(self):
        if self.base_model is None:
            return
        preset = self.medium_preset_combo.currentText()
        if preset == "Reset to original":
            self.reset_medium()
            return
        # Build a quick heuristic preset
        by_id = {self.medium_table.item(i, 0).text().strip(): i for i in range(self.medium_table.rowCount())}
        def set_bound(rid, lb=None, ub=None):
            row = by_id.get(rid)
            if row is None:
                return
            if lb is not None:
                self.medium_table.setItem(row, 2, QTableWidgetItem(f"{float(lb):g}"))
            if ub is not None:
                self.medium_table.setItem(row, 3, QTableWidgetItem(f"{float(ub):g}"))
        # Close all uptakes by default
        for i in range(self.medium_table.rowCount()):
            lb_item = self.medium_table.item(i, 2)
            try:
                lb = float(lb_item.text())
            except Exception:
                lb = 0.0
            if lb < 0:
                self.medium_table.setItem(i, 2, QTableWidgetItem("0"))
        if preset == "Minimal M9":
            set_bound("EX_glc__D_e", -10, None)
            set_bound("EX_o2_e", -20, None)
            set_bound("EX_nh4_e", -10, None)
            set_bound("EX_pi_e", -10, None)
            set_bound("EX_h2o_e", None, 1000)
            set_bound("EX_h_e", None, 1000)
            set_bound("EX_co2_e", None, 1000)
        elif preset == "Rich LB":
            # Allow common carbon/nitrogen sources generously
            for rid in ("EX_glc__D_e", "EX_o2_e", "EX_nh4_e", "EX_pi_e", "EX_glyc_e", "EX_etoh_e"):
                set_bound(rid, -50, None)
            for rid in ("EX_h2o_e", "EX_h_e", "EX_co2_e"):
                set_bound(rid, None, 1000)
        self.filter_medium_table()
        self.statusBar().showMessage(f"Applied medium preset: {preset}")

    def close_all_uptakes(self):
        if self.base_model is None or self.is_running:
            return

        current = self.tabs.currentWidget()

        if current is self.tab_medium:
            changed = 0
            for rxn in self.exchange_rxns:
                if rxn.lower_bound < 0:
                    rxn.lower_bound = 0.0
                    changed += 1
            self.populate_medium_table()
            self.status_lbl.setText(f"All uptakes closed in Medium (LB<0 → 0). Changed: {changed}")
            self.statusBar().showMessage(f"All uptakes closed in Medium (changed {changed}).")
            if changed:
                self.model_dirty = True
            return

        if current is self.tab_rxns:
            selected_rows = sorted({idx.row() for idx in self.reaction_table.selectionModel().selectedRows()})
            if not selected_rows:
                QMessageBox.information(
                    self,
                    "No selection",
                    "Reactions tab: Select one or more rows to apply 'Close uptakes' as reaction overrides (LB<0 → 0).\n\n"
                    "Tip: Use Medium tab to close all exchange uptakes globally."
                )
                return

            rxn_by_id = {r.id: r for r in self.base_model.reactions}
            changed = 0
            skipped = 0

            for row in selected_rows:
                rxn_id = self._reaction_row_to_id(row)
                if not rxn_id or rxn_id not in rxn_by_id:
                    skipped += 1
                    continue
                rxn = rxn_by_id[rxn_id]

                if not is_exchange_reaction(rxn):
                    skipped += 1
                    continue

                current_lb = float(rxn.lower_bound)
                current_ub = float(rxn.upper_bound)
                if rxn_id in self.reaction_bound_overrides:
                    current_lb, current_ub = self.reaction_bound_overrides[rxn_id]

                if float(current_lb) < 0:
                    self.reaction_bound_overrides[rxn_id] = (0.0, float(current_ub))
                    changed += 1

            self.populate_reaction_table()
            self.statusBar().showMessage(f"Reactions tab: uptake overrides applied. changed={changed}, skipped={skipped}")
            if changed:
                self.model_dirty = True
            return

        QMessageBox.information(self, "Close uptakes", "Go to Medium or Reactions tab to use this action.")

    def reset_medium(self):
        if self.base_model is None or self.is_running:
            return
        for rxn in self.exchange_rxns:
            if rxn.id in self.original_bounds:
                lb, ub = self.original_bounds[rxn.id]
                rxn.lower_bound = lb
                rxn.upper_bound = ub
        self.populate_medium_table()
        self.status_lbl.setText("Medium reset to original.")
        self.statusBar().showMessage("Medium reset to original.")
        self.model_dirty = True

    def reset_reactions_to_original(self):
        if self.base_model is None or self.is_running:
            return
        if self.original_model_snapshot is None:
            QMessageBox.warning(self, "No snapshot", "No original snapshot available.")
            return

        msg = (
            "This will reset ALL reaction bounds in the current model to the original SBML snapshot and clear reaction overrides.\n\n"
            "It will NOT delete your SBML file, but it WILL discard current in-app bound edits/overrides.\n\n"
            "Continue?"
        )
        if QMessageBox.question(self, "Confirm reset reactions", msg, QMessageBox.Yes | QMessageBox.No) != QMessageBox.Yes:
            return

        snap = self.original_model_snapshot
        for rxn in self.base_model.reactions:
            try:
                srxn = snap.reactions.get_by_id(rxn.id)
                rxn.lower_bound = float(srxn.lower_bound)
                rxn.upper_bound = float(srxn.upper_bound)
            except KeyError:
                continue

        self.exchange_rxns = [r for r in self.base_model.reactions if is_exchange_reaction(r)]
        self.exchange_rxns.sort(key=lambda r: r.id)
        self.original_bounds = {r.id: (float(r.lower_bound), float(r.upper_bound)) for r in self.exchange_rxns}

        self.reaction_bound_overrides.clear()

        self.populate_medium_table()
        self.populate_reaction_table()
        self.statusBar().showMessage("Reactions reset to original snapshot (and overrides cleared).")
        self.model_dirty = True

    # ---------------- Genes tab ----------------

