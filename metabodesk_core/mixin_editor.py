"""Model editor mixin for MetaboDesk.

Implements the patch-based model editing system.  Users can add
metabolites, genes, and reactions, disable or delete existing reactions,
and change the model objective — all stored as a JSON patch that can be
applied to, or exported alongside, the original SBML.  Includes an
equation parser that converts human-readable reaction strings into
COBRApy ``Reaction`` objects.
"""

import json
import logging

import cobra
from pathlib import Path

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QFileDialog, QMessageBox, QTableWidget, QTableWidgetItem,
    QAbstractItemView, QLineEdit, QCheckBox, QPlainTextEdit,
    QGroupBox, QFormLayout,
)

from metabodesk_core.utils import (
    _guess_compartment_from_met_id, parse_reaction_equation,
)

logger = logging.getLogger("MetaboDesk")

class EditorMixin:
    """Mixin providing editor functionality."""

    def _build_editor_metabolite_tab(self):
        """Build the 'Add metabolite' sub-tab in the Editor."""
        tab = QWidget()
        layout = QVBoxLayout(tab)

        layout.addWidget(self._title("Add metabolite"))

        form = QFormLayout()
        self.ed_met_id = QLineEdit()
        self.ed_met_name = QLineEdit()
        self.ed_met_comp = QLineEdit()
        self.ed_met_comp.setPlaceholderText("e.g. c (cytosol), e (extracellular)")
        self.ed_met_formula = QLineEdit()
        self.ed_met_charge = QLineEdit()
        self.ed_met_charge.setPlaceholderText("integer, optional")

        form.addRow("Metabolite ID:", self.ed_met_id)
        form.addRow("Name:", self.ed_met_name)
        form.addRow("Compartment:", self.ed_met_comp)
        form.addRow("Formula (optional):", self.ed_met_formula)
        form.addRow("Charge (optional):", self.ed_met_charge)

        box = QGroupBox("Metabolite fields")
        box.setLayout(form)
        layout.addWidget(box)

        btn_row = QHBoxLayout()
        self.ed_add_met_btn = QPushButton("Add metabolite to model")
        self.ed_add_met_btn.clicked.connect(self.editor_add_metabolite)
        btn_row.addWidget(self.ed_add_met_btn)
        btn_row.addStretch(1)
        layout.addLayout(btn_row)

        self.ed_met_info = QPlainTextEdit()
        self.ed_met_info.setReadOnly(True)
        self.ed_met_info.setMaximumHeight(160)
        layout.addWidget(self.ed_met_info)

        self.editor_tabs.addTab(tab, "Metabolites")

    def _build_editor_gene_tab(self):
        """Build the 'Add/Update gene' sub-tab in the Editor."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.addWidget(self._title("Add/Update gene"))

        form = QFormLayout()
        self.ed_gene_id = QLineEdit()
        self.ed_gene_name = QLineEdit()
        form.addRow("Gene ID:", self.ed_gene_id)
        form.addRow("Name (optional):", self.ed_gene_name)

        box = QGroupBox("Gene fields")
        box.setLayout(form)
        layout.addWidget(box)

        btn_row = QHBoxLayout()
        self.ed_add_gene_btn = QPushButton("Add/Update gene in model")
        self.ed_add_gene_btn.clicked.connect(self.editor_add_or_update_gene)
        btn_row.addWidget(self.ed_add_gene_btn)
        btn_row.addStretch(1)
        layout.addLayout(btn_row)

        self.ed_gene_info = QPlainTextEdit()
        self.ed_gene_info.setReadOnly(True)
        self.ed_gene_info.setMaximumHeight(160)
        layout.addWidget(self.ed_gene_info)

        self.editor_tabs.addTab(tab, "Genes")

    def _build_editor_reaction_tab(self):
        """Build the 'Add/Update reaction' sub-tab with equation parser."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.addWidget(self._title("Add/Update reaction (with equation parser)"))

        form = QFormLayout()
        self.ed_rxn_id = QLineEdit()
        self.ed_rxn_desc = QLineEdit()
        self.ed_rxn_lb = QLineEdit("-1000")
        self.ed_rxn_ub = QLineEdit("1000")
        self.ed_rxn_rev = QCheckBox("Reversible (force)")
        self.ed_rxn_gpr = QLineEdit()
        self.ed_rxn_kegg = QLineEdit()
        self.ed_rxn_ec = QLineEdit()
        self.ed_rxn_autocreate_missing_mets = QCheckBox("Auto-create missing metabolites (recommended)")
        self.ed_rxn_autocreate_missing_mets.setChecked(True)

        form.addRow("Reaction ID:", self.ed_rxn_id)
        form.addRow("Description:", self.ed_rxn_desc)
        form.addRow("LB:", self.ed_rxn_lb)
        form.addRow("UB:", self.ed_rxn_ub)
        form.addRow("", self.ed_rxn_rev)
        form.addRow("GPR:", self.ed_rxn_gpr)
        form.addRow("KEGG (optional):", self.ed_rxn_kegg)
        form.addRow("EC (optional):", self.ed_rxn_ec)
        form.addRow("", self.ed_rxn_autocreate_missing_mets)

        box = QGroupBox("Reaction fields")
        box.setLayout(form)
        layout.addWidget(box)

        layout.addWidget(QLabel("Equation:"))
        self.ed_rxn_eq = QPlainTextEdit()
        self.ed_rxn_eq.setPlaceholderText("e.g. glc__D_c + 2 atp_c -> 2 adp_c + 2 pi_c + h_c")
        self.ed_rxn_eq.setMaximumHeight(110)
        layout.addWidget(self.ed_rxn_eq)

        btn_row = QHBoxLayout()
        self.ed_parse_eq_btn = QPushButton("Parse equation")
        self.ed_parse_eq_btn.clicked.connect(self.editor_parse_equation)
        btn_row.addWidget(self.ed_parse_eq_btn)

        self.ed_apply_rxn_btn = QPushButton("Add/Update reaction in model")
        self.ed_apply_rxn_btn.clicked.connect(self.editor_add_or_update_reaction)
        btn_row.addWidget(self.ed_apply_rxn_btn)
        btn_row.addStretch(1)
        layout.addLayout(btn_row)

        layout.addWidget(self._title("Parsed stoichiometry preview"))
        self.ed_rxn_preview_tbl = QTableWidget(0, 3)
        self.ed_rxn_preview_tbl.setHorizontalHeaderLabels(["Metabolite", "Coefficient", "Side"])
        self.ed_rxn_preview_tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.ed_rxn_preview_tbl.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._configure_table_for_excel_resize(self.ed_rxn_preview_tbl)
        layout.addWidget(self.ed_rxn_preview_tbl, stretch=1)

        self.ed_rxn_info = QPlainTextEdit()
        self.ed_rxn_info.setReadOnly(True)
        self.ed_rxn_info.setMaximumHeight(160)
        layout.addWidget(self.ed_rxn_info)

        self._editor_last_parsed: dict | None = None

        self.editor_tabs.addTab(tab, "Reactions")

    def _build_editor_disable_delete_tab(self):
        """Build the 'Disable/Delete reactions' sub-tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.addWidget(self._title("Disable / Delete reactions"))

        row = QHBoxLayout()
        row.addWidget(QLabel("Reaction ID:"))
        self.ed_disable_rid = QLineEdit()
        self.ed_disable_rid.setPlaceholderText("Type reaction id...")
        row.addWidget(self.ed_disable_rid, stretch=1)

        self.ed_disable_btn = QPushButton("Disable (LB=0,UB=0)")
        self.ed_disable_btn.clicked.connect(self.editor_disable_reaction)
        row.addWidget(self.ed_disable_btn)

        self.ed_delete_btn = QPushButton("Delete reaction")
        self.ed_delete_btn.clicked.connect(self.editor_delete_reaction)
        row.addWidget(self.ed_delete_btn)

        layout.addLayout(row)

        self.ed_disable_info = QPlainTextEdit()
        self.ed_disable_info.setReadOnly(True)
        layout.addWidget(self.ed_disable_info)

        self.editor_tabs.addTab(tab, "Disable/Delete")

    def _build_editor_constraint_tab(self):
        """Build the 'Custom Constraints' sub-tab in the Editor."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.addWidget(self._title("Custom Constraints"))

        info = QLabel(
            "Add linear constraints or simple reaction-bound overrides to the model.\n"
            "Changes are applied instantly and used in all subsequent analyses."
        )
        info.setWordWrap(True)
        layout.addWidget(info)

        # ---- Simple mode: reaction bounds ----
        simple_group = QGroupBox("Reaction Bound Override")
        sg_layout = QFormLayout(simple_group)

        self.ed_con_rxn = QLineEdit()
        self.ed_con_rxn.setPlaceholderText("Reaction ID (e.g. PFK)")
        sg_layout.addRow("Reaction:", self.ed_con_rxn)

        self.ed_con_lb = QLineEdit("-1000")
        sg_layout.addRow("Lower bound:", self.ed_con_lb)

        self.ed_con_ub = QLineEdit("1000")
        sg_layout.addRow("Upper bound:", self.ed_con_ub)

        self.ed_con_simple_btn = QPushButton("Apply Bound Override")
        self.ed_con_simple_btn.clicked.connect(self._editor_apply_bound_constraint)
        sg_layout.addRow(self.ed_con_simple_btn)
        layout.addWidget(simple_group)

        # ---- Advanced mode: linear expression ----
        adv_group = QGroupBox("Linear Constraint (Advanced)")
        ag_layout = QVBoxLayout(adv_group)

        ag_layout.addWidget(QLabel(
            "Format: coeff * rxn_id + coeff * rxn_id ... [<=, =, >=] bound\n"
            "Example: 1.0 * PFK + -0.5 * PYK <= 10.0"
        ))

        name_row = QHBoxLayout()
        name_row.addWidget(QLabel("Name:"))
        self.ed_con_name = QLineEdit()
        self.ed_con_name.setPlaceholderText("my_constraint_1 (optional)")
        name_row.addWidget(self.ed_con_name)
        ag_layout.addLayout(name_row)

        self.ed_con_expr = QPlainTextEdit()
        self.ed_con_expr.setPlaceholderText("1.0 * PFK + -0.5 * PYK <= 10.0")
        self.ed_con_expr.setMaximumHeight(70)
        ag_layout.addWidget(self.ed_con_expr)

        self.ed_con_adv_btn = QPushButton("Apply Linear Constraint")
        self.ed_con_adv_btn.clicked.connect(self._editor_apply_linear_constraint)
        ag_layout.addWidget(self.ed_con_adv_btn)
        layout.addWidget(adv_group)

        # ---- Info ----
        self.ed_con_info = QPlainTextEdit()
        self.ed_con_info.setReadOnly(True)
        self.ed_con_info.setMaximumHeight(100)
        layout.addWidget(self.ed_con_info)

        layout.addStretch()
        self.editor_tabs.addTab(tab, "Constraints")

    def _editor_apply_bound_constraint(self):
        """Apply a simple reaction bound override from the Constraints tab."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        rid = self.ed_con_rxn.text().strip()
        if not rid or rid not in [r.id for r in self.base_model.reactions]:
            QMessageBox.warning(self, "Invalid", f"Reaction '{rid}' not found in model.")
            return

        try:
            lb = float(self.ed_con_lb.text().strip())
            ub = float(self.ed_con_ub.text().strip())
            if lb > ub:
                raise ValueError("LB > UB")
        except ValueError as e:
            QMessageBox.warning(self, "Invalid bounds", f"Check LB/UB values: {e}")
            return

        self.reaction_bound_overrides[rid] = (lb, ub)
        self.model_dirty = True
        self.ed_con_info.setPlainText(
            f"Applied bound override: {rid}\n"
            f"LB = {lb:g}, UB = {ub:g}\n\n"
            f"This will be applied in all subsequent analyses."
        )
        self._refresh_all_model_views()
        if hasattr(self, '_notify_editor_change'):
            self._notify_editor_change(f"Constraint: {rid} bounds [{lb:g}, {ub:g}]")

    def _editor_apply_linear_constraint(self):
        """Apply a linear constraint expression from the Constraints tab."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        expr_text = self.ed_con_expr.toPlainText().strip()
        if not expr_text:
            QMessageBox.warning(self, "Empty", "Enter a constraint expression.")
            return

        name = self.ed_con_name.text().strip()
        try:
            self._apply_linear_constraint(expr_text, name)
            self.ed_con_info.setPlainText(
                f"Applied linear constraint:\n{expr_text}\n"
                f"Name: {name or '(auto)'}\n\n"
                f"Applied directly to the loaded model."
            )
            self._refresh_all_model_views()
            if hasattr(self, '_notify_editor_change'):
                self._notify_editor_change(f"Linear constraint: {name or expr_text[:40]}")
        except Exception as e:
            QMessageBox.warning(self, "Constraint Error", str(e))

    def _build_editor_patch_tab(self):
        """Build the 'Patch Manager' sub-tab for import/export patches."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.addWidget(self._title("Patch Manager"))

        btn_row = QHBoxLayout()
        self.ed_patch_export_btn = QPushButton("Export patch (JSON)...")
        self.ed_patch_export_btn.clicked.connect(self.editor_export_patch)
        btn_row.addWidget(self.ed_patch_export_btn)

        self.ed_patch_import_btn = QPushButton("Import patch (JSON) & apply...")
        self.ed_patch_import_btn.clicked.connect(self.editor_import_patch_and_apply)
        btn_row.addWidget(self.ed_patch_import_btn)

        self.ed_patch_clear_btn = QPushButton("Clear patch (does NOT reset model)")
        self.ed_patch_clear_btn.clicked.connect(self.editor_clear_patch_only)
        btn_row.addWidget(self.ed_patch_clear_btn)

        btn_row.addStretch(1)
        layout.addLayout(btn_row)

        layout.addWidget(QLabel("Patch preview:"))
        self.ed_patch_preview = QPlainTextEdit()
        self.ed_patch_preview.setReadOnly(True)
        layout.addWidget(self.ed_patch_preview, stretch=1)

        self.editor_tabs.addTab(tab, "Patch")

    def _build_editor_export_sbml_tab(self):
        """Build the 'Export SBML' sub-tab for saving the patched model."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.addWidget(self._title("Export patched SBML"))

        info = QLabel(
            "This exports the CURRENT in-app model (including any changes you made in Model Editor) as a new SBML file.\n"
            "Your original SBML is not modified."
        )
        info.setWordWrap(True)
        layout.addWidget(info)

        btn = QPushButton("Export current model as SBML...")
        btn.clicked.connect(self.editor_export_patched_sbml)
        layout.addWidget(btn)

        layout.addStretch(1)
        self.editor_tabs.addTab(tab, "Export SBML")

    def _refresh_patch_preview(self):
        """Refresh the patch JSON preview text widget."""
        try:
            self.ed_patch_preview.setPlainText(json.dumps(self.model_patch, indent=2))
        except Exception:
            pass

    def _rebuild_editor_reaction_dropdowns(self):
        """No-op stub for reaction dropdown rebuild."""
        return

    def _rebuild_editor_gene_dropdowns(self):
        """No-op stub for gene dropdown rebuild."""
        return

    def _rebuild_editor_metabolite_dropdowns(self):
        """No-op stub for metabolite dropdown rebuild."""
        return

    def _refresh_all_model_views(self):
        """Refresh every model-dependent UI list so editor changes are visible everywhere."""
        for fn_name in (
            "populate_reaction_table",
            "populate_medium_table",
            "populate_gene_list",
            "populate_genes_tab",
            "populate_overexpression_reaction_list",
            "populate_objective_combo",
        ):
            fn = getattr(self, fn_name, None)
            if callable(fn):
                try:
                    fn()
                except Exception:
                    pass

    # ---------- editor actions ----------

    def editor_add_metabolite(self):
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        mid = self.ed_met_id.text().strip()
        if not mid:
            QMessageBox.warning(self, "Missing", "Metabolite ID is required.")
            return

        name = self.ed_met_name.text().strip()
        comp = self.ed_met_comp.text().strip() or _guess_compartment_from_met_id(mid)
        formula = self.ed_met_formula.text().strip()
        charge_txt = self.ed_met_charge.text().strip()

        charge = None
        if charge_txt:
            try:
                charge = int(charge_txt)
            except Exception:
                QMessageBox.warning(self, "Invalid charge", "Charge must be an integer (or blank).")
                return

        if mid in self.base_model.metabolites:
            QMessageBox.information(self, "Exists", f"Metabolite already exists: {mid}\n(You can still use it in reactions.)")
            return

        m = cobra.Metabolite(id=mid, name=name or mid, compartment=comp)
        if formula:
            m.formula = formula
        if charge is not None:
            m.charge = charge

        self.base_model.add_metabolites([m])
        self.model_dirty = True

        self.model_patch["added_metabolites"].append({
            "id": mid,
            "name": name,
            "compartment": comp,
            "formula": formula,
            "charge": charge,
        })
        self._refresh_patch_preview()

        self.ed_met_info.setPlainText(f"Added metabolite: {mid}\ncompartment={comp}\nname={name}\nformula={formula}\ncharge={charge}")
        self._refresh_all_model_views()
        if hasattr(self, '_notify_editor_change'):
            self._notify_editor_change(f"Metabolite added: {mid}")

    def editor_add_or_update_gene(self):
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        gid = self.ed_gene_id.text().strip()
        if not gid:
            QMessageBox.warning(self, "Missing", "Gene ID is required.")
            return
        name = self.ed_gene_name.text().strip()

        if gid in self.base_model.genes:
            g = self.base_model.genes.get_by_id(gid)
            if name:
                try:
                    g.name = name
                except Exception:
                    pass
            action = "Updated"
        else:
            gene = cobra.Gene(id=gid)
            try:
                gene.name = name
            except Exception:
                pass
            self.base_model.genes._dict[gid] = gene
            action = "Added"

        self.model_patch["added_genes"].append({"id": gid, "name": name})
        self._refresh_patch_preview()
        self.model_dirty = True

        self.ed_gene_info.setPlainText(f"{action} gene: {gid}\nname={name}\n\nNote: genes become functionally connected via reaction GPRs.")
        self._refresh_all_model_views()
        if hasattr(self, '_notify_editor_change'):
            self._notify_editor_change(f"Gene {action.lower()}: {gid}")

    def editor_parse_equation(self):
        eq = self.ed_rxn_eq.toPlainText().strip()
        try:
            stoich, reversible = parse_reaction_equation(eq)
        except Exception as e:
            QMessageBox.warning(self, "Parse failed", str(e))
            self._editor_last_parsed = None
            self.ed_rxn_preview_tbl.setRowCount(0)
            return

        self._editor_last_parsed = {"stoich": stoich, "reversible": reversible, "equation": eq}

        self.ed_rxn_preview_tbl.setRowCount(0)
        for met_id in sorted(stoich.keys()):
            coeff = float(stoich[met_id])
            side = "Product" if coeff > 0 else "Reactant"
            row = self.ed_rxn_preview_tbl.rowCount()
            self.ed_rxn_preview_tbl.insertRow(row)
            self.ed_rxn_preview_tbl.setItem(row, 0, QTableWidgetItem(met_id))
            self.ed_rxn_preview_tbl.setItem(row, 1, QTableWidgetItem(f"{coeff:g}"))
            self.ed_rxn_preview_tbl.setItem(row, 2, QTableWidgetItem(side))

        self.ed_rxn_info.setPlainText(
            f"Parsed OK.\nReversible={reversible}\nMetabolites={len(stoich)}\n\nTip: If metabolites are missing, Apply will auto-create them (if enabled)."
        )

    def _ensure_metabolites_exist_for_stoich(self, stoich: dict[str, float]) -> list[str]:
        """Auto-create any metabolites needed by *stoich* that are missing from the model."""
        assert self.base_model is not None
        created: list[str] = []

        for mid in stoich.keys():
            if mid in self.base_model.metabolites:
                continue
            comp = _guess_compartment_from_met_id(mid)
            m = cobra.Metabolite(id=mid, name=mid, compartment=comp)
            self.base_model.add_metabolites([m])
            created.append(mid)

            self.model_patch["added_metabolites"].append({
                "id": mid,
                "name": "",
                "compartment": comp,
                "formula": "",
                "charge": None,
                "auto_created": True,
            })

        return created

    def editor_add_or_update_reaction(self):
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        rid = self.ed_rxn_id.text().strip()
        if not rid:
            QMessageBox.warning(self, "Missing", "Reaction ID is required.")
            return

        desc = self.ed_rxn_desc.text().strip()
        gpr = self.ed_rxn_gpr.text().strip()
        kegg = self.ed_rxn_kegg.text().strip()
        ec = self.ed_rxn_ec.text().strip()

        try:
            lb = self._parse_float(self.ed_rxn_lb.text())
            ub = self._parse_float(self.ed_rxn_ub.text())
            if lb > ub:
                raise ValueError("LB > UB")
        except Exception:
            QMessageBox.warning(self, "Invalid bounds", "LB/UB must be numeric and LB <= UB.")
            return

        eq_text = self.ed_rxn_eq.toPlainText().strip()
        parsed = self._editor_last_parsed
        if not parsed or parsed.get("equation") != eq_text:
            try:
                stoich, reversible = parse_reaction_equation(eq_text)
            except Exception as e:
                QMessageBox.warning(self, "Parse failed", f"Equation parse failed:\n{e}")
                return
        else:
            stoich = parsed["stoich"]
            reversible = bool(parsed["reversible"])

        if self.ed_rxn_rev.isChecked():
            reversible = True

        created_mets: list[str] = []
        if self.ed_rxn_autocreate_missing_mets.isChecked():
            created_mets = self._ensure_metabolites_exist_for_stoich(stoich)
            if created_mets:
                QMessageBox.information(
                    self,
                    "Missing metabolites auto-created",
                    "These metabolites were missing and have been auto-created with empty metadata:\n\n"
                    f"{', '.join(created_mets)}\n\n"
                    "Please add proper names/formula/charge later in Model Editor → Metabolites."
                )
        else:
            missing = [m for m in stoich.keys() if m not in self.base_model.metabolites]
            if missing:
                QMessageBox.warning(
                    self,
                    "Missing metabolites",
                    "These metabolites do not exist in the model:\n\n"
                    f"{', '.join(missing)}\n\n"
                    "Enable auto-create or add them first."
                )
                return

        if rid in self.base_model.reactions:
            rxn = self.base_model.reactions.get_by_id(rid)
            action = "Updated"
        else:
            rxn = cobra.Reaction(id=rid)
            self.base_model.add_reactions([rxn])
            action = "Added"

        try:
            rxn.name = desc
        except Exception as e:
            self.statusBar().showMessage(f"Warning: could not set reaction name: {e}")

        try:
            rxn.subtract_metabolites(rxn.metabolites)
        except Exception as e:
            self.statusBar().showMessage(f"Warning: could not clear old stoichiometry: {e}")

        mets = {self.base_model.metabolites.get_by_id(mid): float(coeff) for mid, coeff in stoich.items()}
        rxn.add_metabolites(mets)

        rxn.lower_bound = float(lb)
        rxn.upper_bound = float(ub)
        rxn.reversibility = bool(reversible)

        try:
            rxn.gene_reaction_rule = gpr
        except Exception as e:
            self.statusBar().showMessage(f"Warning: could not set GPR: {e}")

        ann = getattr(rxn, "annotation", {}) or {}
        if kegg:
            ann["kegg.reaction"] = kegg
        if ec:
            ann["ec-code"] = ec
        rxn.annotation = ann

        self.model_patch["added_or_updated_reactions"].append({
            "id": rid,
            "name": desc,
            "equation": eq_text,
            "stoichiometry": stoich,
            "reversible": reversible,
            "lb": lb,
            "ub": ub,
            "gpr": gpr,
            "annotation": {"kegg.reaction": kegg, "ec-code": ec},
        })
        self._refresh_patch_preview()
        self.model_dirty = True

        self.ed_rxn_info.setPlainText(
            f"{action} reaction: {rid}\n"
            f"LB={lb}, UB={ub}, reversible={reversible}\n"
            f"GPR={gpr}\n"
            f"Equation:\n{eq_text}\n"
        )

        self._refresh_all_model_views()
        if hasattr(self, '_notify_editor_change'):
            self._notify_editor_change(f"Reaction {action.lower()}: {rid}")

    def editor_disable_reaction(self):
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        rid = self.ed_disable_rid.text().strip()
        if not rid:
            QMessageBox.warning(self, "Missing", "Reaction ID is required.")
            return
        if rid not in self.base_model.reactions:
            QMessageBox.warning(self, "Not found", f"Reaction not found: {rid}")
            return

        rxn = self.base_model.reactions.get_by_id(rid)
        rxn.lower_bound = 0.0
        rxn.upper_bound = 0.0

        if rid not in self.model_patch["disabled_reactions"]:
            self.model_patch["disabled_reactions"].append(rid)
        self._refresh_patch_preview()
        self.model_dirty = True

        self.ed_disable_info.setPlainText(f"Disabled reaction: {rid} (LB=0, UB=0)")
        self._refresh_all_model_views()
        if hasattr(self, '_notify_editor_change'):
            self._notify_editor_change(f"Reaction disabled: {rid}")

    def editor_delete_reaction(self):
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        rid = self.ed_disable_rid.text().strip()
        if not rid:
            QMessageBox.warning(self, "Missing", "Reaction ID is required.")
            return
        if rid not in self.base_model.reactions:
            QMessageBox.warning(self, "Not found", f"Reaction not found: {rid}")
            return

        msg = (
            "You are about to DELETE a reaction from the in-app model.\n\n"
            f"Reaction: {rid}\n\n"
            "This does NOT modify your original SBML file, but it WILL remove the reaction from the current session and any exported patched SBML.\n\n"
            "Continue?"
        )
        if QMessageBox.question(self, "Confirm delete reaction", msg, QMessageBox.Yes | QMessageBox.No) != QMessageBox.Yes:
            return

        rxn = self.base_model.reactions.get_by_id(rid)
        try:
            self.base_model.remove_reactions([rxn], remove_orphans=False)
        except Exception as e:
            self._show_error("Delete failed", "Failed to delete reaction.", e)
            return

        if rid not in self.model_patch["deleted_reactions"]:
            self.model_patch["deleted_reactions"].append(rid)
        self._refresh_patch_preview()
        self.model_dirty = True

        self.ed_disable_info.setPlainText(f"Deleted reaction: {rid}")
        self._refresh_all_model_views()
        if hasattr(self, '_notify_editor_change'):
            self._notify_editor_change(f"Reaction deleted: {rid}")

    def editor_export_patch(self):
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export patch",
            str(Path.home() / "metabodesk_patch.json"),
            "JSON files (*.json);;All files (*.*)",
        )
        if not file_path:
            return
        try:
            Path(file_path).write_text(json.dumps(self.model_patch, indent=2), encoding="utf-8")
        except Exception as e:
            self._show_error("Export failed", "Failed to export patch.", e)
            return
        self.statusBar().showMessage(f"Patch exported: {file_path}")

    def editor_import_patch_and_apply(self):
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Import patch",
            str(Path.home()),
            "JSON files (*.json);;All files (*.*)",
        )
        if not file_path:
            return
        try:
            patch = json.loads(Path(file_path).read_text(encoding="utf-8"))
            if not isinstance(patch, dict):
                raise ValueError("Patch must be a JSON object.")
        except Exception as e:
            self._show_error("Import failed", "Failed to import patch.", e)
            return

        msg = (
            "This will apply the imported patch to the CURRENT in-app model.\n"
            "Your original SBML file is not modified.\n\n"
            "Continue?"
        )
        if QMessageBox.question(self, "Confirm apply patch", msg, QMessageBox.Yes | QMessageBox.No) != QMessageBox.Yes:
            return

        try:
            self._apply_patch_to_model(patch)
        except Exception as e:
            self._show_error("Apply failed", "Failed to apply patch.", e)
            return

        self.model_patch = patch
        self._refresh_patch_preview()
        self.statusBar().showMessage(f"Patch imported & applied: {file_path}")
        self.model_dirty = True

        self.populate_reaction_table()
        self.populate_medium_table()
        self.populate_overexpression_reaction_list()
        self.populate_objective_combo()
        self.populate_gene_list()
        self.populate_genes_tab()

    def editor_clear_patch_only(self):
        msg = (
            "This will CLEAR the in-app patch record (JSON history), but it will NOT revert the model.\n\n"
            "If you want to revert the model, use 'Reset reactions' (and re-open SBML if needed).\n\n"
            "Continue?"
        )
        if QMessageBox.question(self, "Confirm clear patch", msg, QMessageBox.Yes | QMessageBox.No) != QMessageBox.Yes:
            return

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
        self.statusBar().showMessage("Patch cleared (model not reverted).")

    def editor_export_patched_sbml(self):
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export patched SBML",
            str(Path.home() / "metabodesk_patched_model.xml"),
            "SBML files (*.xml *.sbml);;All files (*.*)",
        )
        if not file_path:
            return
        try:
            cobra.io.write_sbml_model(self.base_model, file_path)
        except Exception as e:
            self._show_error("Export failed", "Failed to export patched SBML.", e)
            return
        self.statusBar().showMessage(f"Patched SBML exported: {file_path}")

    def _apply_patch_to_model(self, patch: dict):
        """Apply a JSON patch dict (metabolites, genes, reactions, deletions) to the model."""
        assert self.base_model is not None

        for m in patch.get("added_metabolites", []) or []:
            mid = str(m.get("id", "")).strip()
            if not mid or mid in self.base_model.metabolites:
                continue
            comp = str(m.get("compartment") or _guess_compartment_from_met_id(mid))
            name = str(m.get("name") or mid)
            met = cobra.Metabolite(id=mid, name=name, compartment=comp)
            self.base_model.add_metabolites([met])

        for g in patch.get("added_genes", []) or []:
            gid = str(g.get("id", "")).strip()
            if not gid:
                continue
            if gid in self.base_model.genes:
                continue
            gene = cobra.Gene(id=gid)
            try:
                gene.name = str(g.get("name") or "")
            except Exception:
                pass
            self.base_model.genes._dict[gid] = gene

        for r in patch.get("added_or_updated_reactions", []) or []:
            rid = str(r.get("id", "")).strip()
            if not rid:
                continue
            eq = str(r.get("equation") or "").strip()
            stoich = r.get("stoichiometry")
            if not isinstance(stoich, dict):
                stoich, _rev = parse_reaction_equation(eq)
            else:
                stoich = {str(k): float(v) for k, v in stoich.items()}

            for mid in list(stoich.keys()):
                if mid not in self.base_model.metabolites:
                    comp = _guess_compartment_from_met_id(mid)
                    self.base_model.add_metabolites([cobra.Metabolite(id=mid, name=mid, compartment=comp)])

            if rid in self.base_model.reactions:
                rxn = self.base_model.reactions.get_by_id(rid)
            else:
                rxn = cobra.Reaction(id=rid)
                self.base_model.add_reactions([rxn])

            try:
                rxn.name = str(r.get("name") or "")
            except Exception as e:
                self.statusBar().showMessage(f"Warning: could not set name for {rid}: {e}")

            try:
                rxn.subtract_metabolites(rxn.metabolites)
            except Exception as e:
                self.statusBar().showMessage(f"Warning: could not clear stoichiometry for {rid}: {e}")
            mets = {self.base_model.metabolites.get_by_id(mid): float(coeff) for mid, coeff in stoich.items()}
            rxn.add_metabolites(mets)

            try:
                rxn.lower_bound = float(r.get("lb"))
                rxn.upper_bound = float(r.get("ub"))
            except Exception as e:
                self.statusBar().showMessage(f"Warning: could not set bounds for {rid}: {e}")

            try:
                rxn.gene_reaction_rule = str(r.get("gpr") or "")
            except Exception as e:
                self.statusBar().showMessage(f"Warning: could not set GPR for {rid}: {e}")

            ann = getattr(rxn, "annotation", {}) or {}
            ann_in = r.get("annotation") or {}
            if isinstance(ann_in, dict):
                ann.update({k: v for k, v in ann_in.items() if v})
            rxn.annotation = ann

        for rid in patch.get("disabled_reactions", []) or []:
            rid = str(rid)
            if rid in self.base_model.reactions:
                rxn = self.base_model.reactions.get_by_id(rid)
                rxn.lower_bound = 0.0
                rxn.upper_bound = 0.0

        for rid in patch.get("deleted_reactions", []) or []:
            rid = str(rid)
            if rid in self.base_model.reactions:
                rxn = self.base_model.reactions.get_by_id(rid)
                self.base_model.remove_reactions([rxn], remove_orphans=False)

    # ---------------- Save As helpers (SBML) ----------------

    def _export_current_model_to_path(self, file_path: Path):
        """Write the current in-app model to an SBML file at *file_path*."""
        if self.base_model is None:
            raise ValueError("No model loaded.")
        cobra.io.write_sbml_model(self.base_model, str(file_path))

    def save_current_model_as(self) -> Path | None:
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "No model loaded.")
            return None

        default_name = "metabodesk_current_model.xml"
        if self.current_sbml_path:
            default_name = f"{self.current_sbml_path.stem}_metabodesk.xml"

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save current model as (SBML)",
            str(Path.home() / default_name),
            "SBML files (*.xml *.sbml);;All files (*.*)",
        )
        if not file_path:
            return None

        out = Path(file_path)
        self._export_current_model_to_path(out)
        self.current_sbml_path = out
        self.model_dirty = False
        self.statusBar().showMessage(f"Saved SBML: {out}")
        return out

    # ---------------- (2) memote helpers ----------------

