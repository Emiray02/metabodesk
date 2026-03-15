"""Reactions and gene tables mixin for MetaboDesk.

Populates and filters the all-reactions table (with KEGG, EC, GPR metadata),
the gene→reaction mapping table, the gene knockout panel, and the
overexpression (reaction-bound scaling) panel.  Provides inline bound
editing with undo-safe state tracking.
"""

import logging

import cobra
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMessageBox, QTableWidgetItem, QListWidgetItem

from metabodesk_core.utils import (
    get_kegg_rxn_id, get_ec_numbers, get_description, get_gpr,
)

logger = logging.getLogger("MetaboDesk")

class ReactionsMixin:
    """Mixin providing reactions functionality."""

    def populate_genes_tab(self):
        if self.base_model is None:
            self.genes_tab_list.clear()
            self.gene_rxn_table.setRowCount(0)
            self._gene_to_reactions = {}
            return

        self._gene_to_reactions = {}
        for g in self.base_model.genes:
            try:
                rids = sorted([r.id for r in g.reactions])
            except Exception:
                rids = []
            self._gene_to_reactions[str(g.id)] = rids

        self.genes_tab_list.clear()
        for g in sorted(self.base_model.genes, key=lambda gg: gg.id):
            gname = getattr(g, "name", "") or ""
            label = f"{g.id} — {gname}".strip(" —")
            it = QListWidgetItem(label)
            it.setData(Qt.UserRole, str(g.id))
            self.genes_tab_list.addItem(it)

        self.filter_genes_tab_gene_list()
        self.gene_rxn_table.setRowCount(0)
        self.gene_rxn_details.setPlainText("")

    def filter_genes_tab_gene_list(self):
        text = self.gene_global_search.text().strip().lower()
        for i in range(self.genes_tab_list.count()):
            it = self.genes_tab_list.item(i)
            gid = str(it.data(Qt.UserRole) or "")
            gname = it.text().lower()
            it.setHidden(bool(text) and (text not in gid.lower()) and (text not in gname))

    def on_genes_tab_gene_changed(self, current: QListWidgetItem, _prev: QListWidgetItem):
        if self.base_model is None or not current:
            self.gene_rxn_table.setRowCount(0)
            return

        gid = str(current.data(Qt.UserRole) or "").strip()
        if not gid:
            self.gene_rxn_table.setRowCount(0)
            return

        rxn_ids = self._gene_to_reactions.get(gid, [])
        self.populate_gene_rxn_table(rxn_ids)

    def populate_gene_rxn_table(self, rxn_ids: list[str]):
        if self.base_model is None:
            return

        # Collect valid reactions first, then batch-fill the table
        rows_data = []
        readonly = Qt.ItemIsSelectable | Qt.ItemIsEnabled
        for rid in rxn_ids:
            try:
                rxn = self.base_model.reactions.get_by_id(rid)
            except Exception:
                continue
            eq = ""
            try:
                eq = rxn.reaction
            except Exception:
                pass
            rows_data.append((rxn.id, get_description(rxn), eq, get_gpr(rxn)))

        self.gene_rxn_table.blockSignals(True)
        self.gene_rxn_table.setSortingEnabled(False)
        self.gene_rxn_table.setRowCount(len(rows_data))

        for row, (rid, desc, eq, gpr) in enumerate(rows_data):
            item_id = QTableWidgetItem(rid)
            item_desc = QTableWidgetItem(desc)
            item_eq = QTableWidgetItem(eq)
            item_gpr = QTableWidgetItem(gpr)

            for it in (item_id, item_desc, item_eq, item_gpr):
                it.setFlags(readonly)

            self.gene_rxn_table.setItem(row, 0, item_id)
            self.gene_rxn_table.setItem(row, 1, item_desc)
            self.gene_rxn_table.setItem(row, 2, item_eq)
            self.gene_rxn_table.setItem(row, 3, item_gpr)

        self.gene_rxn_table.setSortingEnabled(True)
        self.gene_rxn_table.blockSignals(False)
        self.gene_rxn_details.setPlainText("")

    # ---------------- Reactions table ----------------

    def _parse_metabolite_filter_terms(self) -> list[str]:
        """Parse the metabolite filter text into a list of lowercase search terms."""
        raw = self.met_filter.text().strip()
        if not raw:
            return []
        parts = [p.strip().lower() for p in raw.replace(";", ",").split(",")]
        return [p for p in parts if p]

    def _reaction_row_to_id(self, row: int) -> str | None:
        """Return the reaction ID from the reaction table at the given row index."""
        it = self.reaction_table.item(row, 0)
        return it.text().strip() if it else None

    def _rxn_contains_term(self, rxn: cobra.Reaction, term: str) -> bool:
        """Check if any metabolite in a reaction matches a search term."""
        term = term.strip().lower()
        if not term:
            return True
        for m in rxn.metabolites.keys():
            mid = (m.id or "").lower()
            mn = (m.name or "").lower()
            if term in mid or term in mn:
                return True
        return False

    def _reaction_matches_filters(self, rxn: cobra.Reaction, text: str, met_terms: list[str], mode: str) -> bool:
        """Return True if a reaction matches both text search and metabolite filter terms."""
        text = text.strip().lower()
        if text:
            hay = " ".join([
                rxn.id.lower(),
                (get_description(rxn) or "").lower(),
                (getattr(rxn, "reaction", "") or "").lower(),
                get_kegg_rxn_id(rxn).lower(),
                get_ec_numbers(rxn).lower(),
                get_gpr(rxn).lower(),
            ])
            if text not in hay:
                return False

        if met_terms:
            if mode == "AND":
                for t in met_terms:
                    if not self._rxn_contains_term(rxn, t):
                        return False
            else:
                if not any(self._rxn_contains_term(rxn, t) for t in met_terms):
                    return False

        return True

    def filter_reaction_table(self):
        if self.base_model is None:
            return
        text = self.rxn_global_search.text()
        met_terms = self._parse_metabolite_filter_terms()
        mode = self.met_filter_mode.currentText().strip().upper() or "AND"

        rxn_by_id = {r.id: r for r in self.base_model.reactions}

        for row in range(self.reaction_table.rowCount()):
            rxn_id = self._reaction_row_to_id(row)
            rxn = rxn_by_id.get(rxn_id) if rxn_id else None
            if rxn is None:
                self.reaction_table.setRowHidden(row, True)
                continue
            hide = not self._reaction_matches_filters(rxn, text, met_terms, mode)
            self.reaction_table.setRowHidden(row, hide)

    def populate_reaction_table(self):
        if self.base_model is None:
            return

        self.reaction_table.blockSignals(True)
        self.reaction_table.setSortingEnabled(False)
        sorted_rxns = sorted(self.base_model.reactions, key=lambda r: r.id)
        self.reaction_table.setRowCount(len(sorted_rxns))

        readonly = Qt.ItemIsSelectable | Qt.ItemIsEnabled
        for row, rxn in enumerate(sorted_rxns):
            item_id = QTableWidgetItem(rxn.id)
            item_desc = QTableWidgetItem(get_description(rxn))
            eq = ""
            try:
                eq = rxn.reaction
            except Exception:
                eq = ""
            item_eq = QTableWidgetItem(eq)
            item_kegg = QTableWidgetItem(get_kegg_rxn_id(rxn))
            item_ec = QTableWidgetItem(get_ec_numbers(rxn))
            item_gpr = QTableWidgetItem(get_gpr(rxn))

            for it in (item_id, item_desc, item_eq, item_kegg, item_ec, item_gpr):
                it.setFlags(readonly)

            if rxn.id in self.reaction_bound_overrides:
                lb, ub = self.reaction_bound_overrides[rxn.id]
            else:
                lb, ub = float(rxn.lower_bound), float(rxn.upper_bound)

            item_lb = QTableWidgetItem(str(lb))
            item_ub = QTableWidgetItem(str(ub))

            self.reaction_table.setItem(row, 0, item_id)
            self.reaction_table.setItem(row, 1, item_desc)
            self.reaction_table.setItem(row, 2, item_eq)
            self.reaction_table.setItem(row, 3, item_kegg)
            self.reaction_table.setItem(row, 4, item_ec)
            self.reaction_table.setItem(row, 5, item_gpr)
            self.reaction_table.setItem(row, 6, item_lb)
            self.reaction_table.setItem(row, 7, item_ub)

        self.reaction_table.setSortingEnabled(True)
        self.reaction_table.blockSignals(False)
        self.filter_reaction_table()
        self.clear_rxn_overrides_btn.setEnabled(True)
        self.remove_selected_overrides_btn.setEnabled(True)

    def on_reaction_table_item_changed(self, item: QTableWidgetItem):
        # FIXED: remove invalid usage of rxn_id/lb/ub before definition.
        if self.base_model is None:
            return

        row = item.row()
        col = item.column()
        if col not in (6, 7):
            return

        rxn_id = self._reaction_row_to_id(row)
        if not rxn_id:
            return

        lb_item = self.reaction_table.item(row, 6)
        ub_item = self.reaction_table.item(row, 7)
        if not lb_item or not ub_item:
            return

        try:
            lb = self._parse_float(lb_item.text())
            ub = self._parse_float(ub_item.text())
            if lb > ub:
                raise ValueError("LB > UB")
        except Exception:
            QMessageBox.warning(self, "Invalid bounds", f"Invalid bounds for {rxn_id}. Please enter valid numbers with LB <= UB.")
            self.reaction_table.blockSignals(True)
            if rxn_id in self.reaction_bound_overrides:
                old_lb, old_ub = self.reaction_bound_overrides[rxn_id]
            else:
                rxn = self.base_model.reactions.get_by_id(rxn_id)
                old_lb, old_ub = float(rxn.lower_bound), float(rxn.upper_bound)
            lb_item.setText(str(old_lb))
            ub_item.setText(str(old_ub))
            self.reaction_table.blockSignals(False)
            return

        self.save_state_for_undo()
        self.reaction_bound_overrides[rxn_id] = (lb, ub)
        self.model_dirty = True
        self.statusBar().showMessage(f"Override set: {rxn_id} LB={lb:g}, UB={ub:g}")

    def clear_reaction_overrides(self):
        self.save_state_for_undo()
        self.reaction_bound_overrides.clear()
        self.populate_reaction_table()
        self.model_dirty = True
        self.statusBar().showMessage("Reaction overrides cleared.")

    def remove_selected_reaction_overrides(self):
        if self.base_model is None:
            return

        selected_rows = sorted({idx.row() for idx in self.reaction_table.selectionModel().selectedRows()})
        if not selected_rows:
            QMessageBox.information(self, "No selection", "Select one or more reactions (rows) first.")
            return

        removed = 0
        for row in selected_rows:
            rxn_id = self._reaction_row_to_id(row)
            if not rxn_id:
                continue
            if rxn_id in self.reaction_bound_overrides:
                self.reaction_bound_overrides.pop(rxn_id, None)
                removed += 1

        self.populate_reaction_table()
        if removed:
            self.model_dirty = True
        self.statusBar().showMessage(f"Removed {removed} reaction override(s).")

    def apply_reaction_overrides_to_model(self, model: cobra.Model):
        for rxn_id, (lb, ub) in self.reaction_bound_overrides.items():
            try:
                rxn = model.reactions.get_by_id(rxn_id)
            except KeyError:
                continue
            rxn.lower_bound = float(lb)
            rxn.upper_bound = float(ub)

    # ---------------- KO / OE ----------------

    def filter_gene_list(self):
        text = self.gene_search.text().strip().lower()
        for i in range(self.gene_list.count()):
            item = self.gene_list.item(i)
            item.setHidden(bool(text) and text not in item.text().lower())

    def populate_gene_list(self):
        if self.base_model is None:
            return
        self.gene_list.clear()
        for gid in sorted([g.id for g in self.base_model.genes]):
            self.gene_list.addItem(QListWidgetItem(gid))
        self.filter_gene_list()

    def add_selected_gene_knockout(self):
        item = self.gene_list.currentItem()
        if not item:
            return
        gid = item.text().strip()
        if gid:
            self.save_state_for_undo()
            self.knockout_genes.add(gid)
            self.model_dirty = True
            self.update_knockout_label()

    def clear_knockouts(self):
        self.save_state_for_undo()
        self.knockout_genes.clear()
        self.model_dirty = True
        self.update_knockout_label()

    def update_knockout_label(self):
        if not self.knockout_genes:
            self.ko_lbl.setText("Gene knockouts: (none)")
        else:
            self.ko_lbl.setText(f"Gene knockouts ({len(self.knockout_genes)}): {', '.join(sorted(self.knockout_genes))}")

    def filter_overexpression_reaction_list(self):
        text = self.rxn_search.text().strip().lower()
        for i in range(self.rxn_list.count()):
            item = self.rxn_list.item(i)
            item.setHidden(bool(text) and text not in item.text().lower())

    def populate_overexpression_reaction_list(self):
        if self.base_model is None:
            return
        self.rxn_list.clear()
        for r in sorted(self.base_model.reactions, key=lambda x: x.id):
            desc = get_description(r)
            text = f"{r.id} — {desc}" if desc else r.id
            item = QListWidgetItem(text)
            item.setData(Qt.UserRole, r.id)
            self.rxn_list.addItem(item)
        self.filter_overexpression_reaction_list()

    def refresh_overexpression_lists(self):
        self.oe_active_list.clear()
        for rid in sorted(self.overexpression_reactions):
            factor = self.overexpression_reactions[rid]
            it = QListWidgetItem(f"{rid} ×{factor:g}")
            it.setData(Qt.UserRole, rid)
            self.oe_active_list.addItem(it)

        self.temp_active_list.clear()
        for rid in sorted(self.temporary_upper_bound_overrides):
            ub = self.temporary_upper_bound_overrides[rid]
            it = QListWidgetItem(f"{rid} UB={ub:g}")
            it.setData(Qt.UserRole, rid)
            self.temp_active_list.addItem(it)

    def add_selected_reaction_overexpression(self):
        item = self.rxn_list.currentItem()
        if not item:
            return
        rxn_id = item.data(Qt.UserRole)
        if not rxn_id:
            return
        self.save_state_for_undo()
        self.overexpression_reactions[str(rxn_id)] = float(self.oe_factor.value())
        self.model_dirty = True
        self.refresh_overexpression_lists()
        self.update_selected_rxn_info()

    def remove_selected_reaction_overexpression(self):
        item = self.rxn_list.currentItem()
        if not item:
            return
        rxn_id = item.data(Qt.UserRole)
        if not rxn_id:
            return
        self.save_state_for_undo()
        if str(rxn_id) in self.overexpression_reactions:
            self.model_dirty = True
        self.overexpression_reactions.pop(str(rxn_id), None)
        self.refresh_overexpression_lists()
        self.update_selected_rxn_info()

    def clear_overexpressions(self):
        self.save_state_for_undo()
        if self.overexpression_reactions:
            self.model_dirty = True
        self.overexpression_reactions.clear()
        self.refresh_overexpression_lists()
        self.update_selected_rxn_info()

    def set_temp_upper_bound_for_selected(self):
        item = self.rxn_list.currentItem()
        if not item:
            return
        rxn_id = item.data(Qt.UserRole)
        if not rxn_id:
            return
        self.save_state_for_undo()
        self.temporary_upper_bound_overrides[str(rxn_id)] = float(self.temp_ub_spin.value())
        self.model_dirty = True
        self.refresh_overexpression_lists()
        self.update_selected_rxn_info()

    def remove_temp_upper_bound_for_selected(self):
        item = self.rxn_list.currentItem()
        if not item:
            return
        rxn_id = item.data(Qt.UserRole)
        if not rxn_id:
            return
        self.save_state_for_undo()
        if str(rxn_id) in self.temporary_upper_bound_overrides:
            self.model_dirty = True
        self.temporary_upper_bound_overrides.pop(str(rxn_id), None)
        self.refresh_overexpression_lists()
        self.update_selected_rxn_info()

    def clear_temp_upper_bounds(self):
        self.save_state_for_undo()
        if self.temporary_upper_bound_overrides:
            self.model_dirty = True
        self.temporary_upper_bound_overrides.clear()
        self.refresh_overexpression_lists()
        self.update_selected_rxn_info()

    def jump_to_reaction_from_item(self, item: QListWidgetItem):
        rxn_id = item.data(Qt.UserRole)
        if not rxn_id:
            return
        for i in range(self.rxn_list.count()):
            it = self.rxn_list.item(i)
            if str(it.data(Qt.UserRole)) == str(rxn_id):
                self.rxn_list.setCurrentItem(it)
                self.rxn_list.scrollToItem(it)
                break

    def update_selected_rxn_info(self, *_):
        if self.base_model is None:
            self.rxn_info_lbl.setText("Selected reaction: (none)")
            self.rxn_preview_lbl.setText("Preview: -")
            return
        item = self.rxn_list.currentItem()
        if not item:
            self.rxn_info_lbl.setText("Selected reaction: (none)")
            self.rxn_preview_lbl.setText("Preview: -")
            return
        rxn_id = item.data(Qt.UserRole)
        if not rxn_id:
            return
        rxn = self.base_model.reactions.get_by_id(rxn_id)
        lb = float(rxn.lower_bound)
        ub = float(rxn.upper_bound)
        self.rxn_info_lbl.setText(
            f"Selected: {rxn.id}\n"
            f"Description: {get_description(rxn)}\n"
            f"Current bounds: LB={lb:g}, UB={ub:g}\n"
            f"GPR: {get_gpr(rxn)}"
        )
        factor = float(self.oe_factor.value())
        prev_lb = lb * factor if (self.oe_scale_lb.isChecked() and lb < 0) else lb
        prev_ub = ub * factor
        self.rxn_preview_lbl.setText(f"Preview if overexpression applied: LB={prev_lb:g}, UB={prev_ub:g}")

    def apply_temp_upper_bound_overrides(self, model: cobra.Model):
        for rxn_id, ub in self.temporary_upper_bound_overrides.items():
            try:
                model.reactions.get_by_id(rxn_id).upper_bound = float(ub)
            except KeyError:
                continue

    def apply_overexpression_to_model(self, model: cobra.Model):
        scale_lb = bool(self.oe_scale_lb.isChecked())
        for rxn_id, factor in self.overexpression_reactions.items():
            try:
                rxn = model.reactions.get_by_id(rxn_id)
            except KeyError:
                continue
            if scale_lb and float(rxn.lower_bound) < 0:
                rxn.lower_bound = float(rxn.lower_bound) * float(factor)
            rxn.upper_bound = float(rxn.upper_bound) * float(factor)

    def apply_knockouts_to_model(self, model: cobra.Model):
        for gid in sorted(self.knockout_genes):
            if gid not in model.genes:
                raise ValueError(f"Gene not found in model: {gid}")
            model.genes.get_by_id(gid).knock_out()

    # ---------------- Analysis functions (FVA, PFBA) ----------------

