import sys
import json
import csv
import re
from datetime import datetime
from pathlib import Path
from threading import Thread

from PySide6.QtCore import Qt, QUrl, QProcess, QSettings, QTimer
from PySide6.QtGui import QAction, QKeySequence, QDesktopServices, QShortcut
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QFileDialog, QMessageBox, QTableWidget,
    QTableWidgetItem, QAbstractItemView, QHeaderView, QListWidget,
    QListWidgetItem, QLineEdit, QSpinBox, QTabWidget, QDoubleSpinBox,
    QCheckBox, QComboBox, QPlainTextEdit, QDialog, QDialogButtonBox,
    QCompleter, QGroupBox, QFormLayout, QScrollArea, QRadioButton,
    QProgressBar, QProgressDialog
)

import cobra
import networkx as nx
try:
    import openpyxl
    HAS_OPENPYXL = True
except ImportError:
    HAS_OPENPYXL = False

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


RECENT_MAX = 5
RECENT_FILE = Path.home() / ".metabodesk_recent.json"


# ---------------- Basic helpers ----------------
def is_exchange_reaction(rxn: cobra.Reaction) -> bool:
    try:
        if rxn.id.startswith("EX_"):
            return True
        return len(rxn.metabolites) == 1
    except Exception:
        return False


def _as_str_list(v) -> list[str]:
    if v is None:
        return []
    if isinstance(v, str):
        s = v.strip()
        return [s] if s else []
    if isinstance(v, (list, tuple, set)):
        out = []
        for x in v:
            if x is None:
                continue
            s = str(x).strip()
            if s:
                out.append(s)
        return out
    s = str(v).strip()
    return [s] if s else []


def _clean_identifiers_org_id(url_or_id: str) -> str:
    s = str(url_or_id).strip()
    if not s:
        return ""
    if "identifiers.org/" in s:
        s = s.split("identifiers.org/")[-1]
    if "/" in s:
        s = s.split("/")[-1]
    return s.strip()


def get_kegg_rxn_id(rxn: cobra.Reaction) -> str:
    ann = getattr(rxn, "annotation", {}) or {}

    keys = [
        "kegg.reaction",
        "kegg_reaction",
        "kegg",
        "kegg.rxn",
        "kegg_rxn",
        "KEGG Reaction",
        "KEGG_REACTION",
        "reaction.kegg",
    ]
    for k in keys:
        if k in ann:
            vals = [_clean_identifiers_org_id(x) for x in _as_str_list(ann.get(k))]
            vals = [v for v in vals if v]
            if vals:
                return ";".join(vals)

    for k in ("identifiers", "identifiers.org"):
        if k in ann:
            vals = _as_str_list(ann.get(k))
            cleaned = []
            for s in vals:
                if "kegg.reaction" in str(s):
                    cleaned.append(_clean_identifiers_org_id(s))
            cleaned = [c for c in cleaned if c]
            if cleaned:
                return ";".join(cleaned)

    return ""


def get_ec_numbers(rxn: cobra.Reaction) -> str:
    ann = getattr(rxn, "annotation", {}) or {}
    keys = [
        "ec-code",
        "ec",
        "ec_number",
        "ec-number",
        "ec-code(s)",
        "EC Number",
        "EC",
        "enzyme.ec",
    ]
    for k in keys:
        if k in ann:
            vals = _as_str_list(ann.get(k))
            vals = [v.strip() for v in vals if str(v).strip()]
            if vals:
                return ";".join(vals)
    return ""


def get_description(rxn: cobra.Reaction) -> str:
    if getattr(rxn, "name", None):
        return str(rxn.name)
    ann = getattr(rxn, "annotation", {}) or {}
    for k in ("description", "desc", "note", "notes"):
        if k in ann:
            vals = _as_str_list(ann.get(k))
            if vals:
                return vals[0]
    return ""


def get_gpr(rxn: cobra.Reaction) -> str:
    try:
        return str(rxn.gene_reaction_rule or "")
    except Exception:
        return ""


# ---------------- Equation parser ----------------
ARROW_PATTERNS = [
    ("<=>", True),
    ("<->", True),
    ("<-->", True),
    ("->", False),
    ("-->", False),
    ("=>", False),
]


def _guess_compartment_from_met_id(met_id: str) -> str:
    """
    Heuristic: metabolite ids like glc__D_c -> compartment 'c'
    If no '_x' suffix, fallback to 'c'.
    """
    m = re.match(r"^(.+)_([A-Za-z])$", met_id.strip())
    if m:
        return m.group(2)
    return "c"


def _parse_coeff_and_met(term: str) -> tuple[float, str]:
    """
    Parse one side term. Accept:
      - "2 atp_c"
      - "0.5 o2_c"
      - "atp_c" (=> coeff=1)
    """
    term = term.strip()
    if not term:
        raise ValueError("Empty term in equation.")
    parts = term.split()
    if len(parts) == 1:
        return 1.0, parts[0]
    try:
        coeff = float(parts[0])
        met = " ".join(parts[1:]).strip()
        if not met:
            raise ValueError("Missing metabolite id after coefficient.")
        return coeff, met
    except ValueError:
        return 1.0, term


def parse_reaction_equation(equation: str) -> tuple[dict[str, float], bool]:
    """
    Return (stoichiometry, reversible)
    stoich: reactants negative, products positive.
    """
    eq = (equation or "").strip()
    if not eq:
        raise ValueError("Equation is empty.")

    arrow = None
    reversible = False
    for a, rev in ARROW_PATTERNS:
        if a in eq:
            arrow = a
            reversible = rev
            break
    if not arrow:
        raise ValueError("No valid arrow found. Use ->, =>, <=>, <->")

    left, right = eq.split(arrow, 1)
    left = left.strip()
    right = right.strip()

    stoich: dict[str, float] = {}

    def add(met_id: str, coeff: float):
        met_id = met_id.strip()
        if not met_id:
            return
        stoich[met_id] = float(stoich.get(met_id, 0.0) + coeff)

    if left:
        for term in [t.strip() for t in left.split("+")]:
            if not term:
                continue
            coeff, met = _parse_coeff_and_met(term)
            add(met, -abs(coeff))

    if right:
        for term in [t.strip() for t in right.split("+")]:
            if not term:
                continue
            coeff, met = _parse_coeff_and_met(term)
            add(met, abs(coeff))

    stoich = {k: v for k, v in stoich.items() if abs(v) > 1e-12}
    if not stoich:
        raise ValueError("Parsed stoichiometry is empty.")

    return stoich, reversible


# ---------------- UI popups ----------------
class TextPopup(QDialog):
    def __init__(self, title: str, text: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        layout = QVBoxLayout(self)

        self.text = QPlainTextEdit()
        self.text.setReadOnly(True)
        self.text.setPlainText(text or "")
        layout.addWidget(self.text)

        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(self.reject)
        btns.accepted.connect(self.accept)
        btns.clicked.connect(self.reject)
        layout.addWidget(btns)

        self.resize(900, 450)


class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=8, height=3.5, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi, tight_layout=True)
        self.ax = fig.add_subplot(111)
        super().__init__(fig)
        self.setParent(parent)
        self.clicked_node = None


# ---------------- Main app ----------------
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        # Speed up startup by deferring repaints during UI construction
        self.setUpdatesEnabled(False)
        self.setWindowTitle("MetaboDesk")

        self.base_model: cobra.Model | None = None
        self.current_sbml_path: Path | None = None
        self.model_dirty: bool = False
        self._tools_proc: QProcess | None = None
        self._memote_proc: QProcess | None = None
        self.original_model_snapshot: cobra.Model | None = None

        self.exchange_rxns: list[cobra.Reaction] = []
        self.original_bounds: dict[str, tuple[float, float]] = {}

        self.knockout_genes: set[str] = set()
        self.overexpression_reactions: dict[str, float] = {}
        self.temporary_upper_bound_overrides: dict[str, float] = {}
        self.reaction_bound_overrides: dict[str, tuple[float, float]] = {}

        self.model_patch: dict = {
            "version": 1,
            "added_metabolites": [],
            "added_genes": [],
            "added_or_updated_reactions": [],
            "disabled_reactions": [],
            "deleted_reactions": [],
            "objective_reaction_id": None,
            "notes": "",
        }

        self._gene_to_reactions: dict[str, list[str]] = {}

        self.is_running = False
        self.last_run: dict | None = None

        self.recent_scenarios: list[str] = []
        self._load_recent_scenarios()

        # Undo/Redo stack
        self.undo_stack: list[dict] = []
        self.redo_stack: list[dict] = []
        self.max_undo_steps = 50

        # Search/Filter
        self.search_results: list[dict] = []

        # Sensitivity analysis
        self.sensitivity_results: dict = {}

        root = QWidget()
        outer = QVBoxLayout(root)

        controls = QHBoxLayout()

        self.open_btn = QPushButton("Open SBML...")
        self.open_btn.clicked.connect(self.open_sbml)

        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.setEnabled(True)
        self.run_btn.setToolTip("Execute the selected analysis type on the loaded model.")
        self.run_btn.clicked.connect(self.run_fba)

        self.analysis_type = QComboBox()
        self.analysis_type.addItems([
            "FBA", "pFBA", "FVA",
            "Single Gene Deletion (SGD)", "Single Reaction Deletion (SRD)",
            "Robustness Analysis", "Production Envelope", "Flux Sampling"
        ])
        self.analysis_type.setEnabled(False)
        self.analysis_type.setToolTip("Çalıştırılacak analiz tipini seçin.")

        self.close_uptakes_btn = QPushButton("Close uptakes")
        self.close_uptakes_btn.setEnabled(False)
        self.close_uptakes_btn.clicked.connect(self.close_all_uptakes)

        self.reset_medium_btn = QPushButton("Reset medium")
        self.reset_medium_btn.setEnabled(False)
        self.reset_medium_btn.clicked.connect(self.reset_medium)

        self.reset_reactions_btn = QPushButton("Reset reactions")
        self.reset_reactions_btn.setEnabled(False)
        self.reset_reactions_btn.clicked.connect(self.reset_reactions_to_original)

        controls.addWidget(self.open_btn)
        controls.addWidget(self.run_btn)
        controls.addWidget(QLabel("Analysis:"))
        controls.addWidget(self.analysis_type)
        controls.addSpacing(15)
        controls.addWidget(self.close_uptakes_btn)
        controls.addWidget(self.reset_medium_btn)
        controls.addWidget(self.reset_reactions_btn)
        controls.addSpacing(15)

        controls.addWidget(QLabel("Top N:"))        
        self.topn_spin = QSpinBox()
        self.topn_spin.setRange(5, 50)
        self.topn_spin.setValue(20)
        self.topn_spin.setToolTip("Grafik/tablolar için en çok gösterilecek öğe sayısı.")
        controls.addWidget(self.topn_spin)

        controls.addSpacing(15)
        controls.addWidget(QLabel("Solver:"))
        self.solver_combo = QComboBox()
        self.solver_combo.setToolTip("Optimizasyon çözücüsü (Auto: mevcutı kullan).")
        self.solver_combo.addItems(["Auto", "glpk", "highs"])  # Free solvers only
        self.solver_combo.setEnabled(False)
        controls.addWidget(self.solver_combo, stretch=1)
        controls.addSpacing(15)
        controls.addWidget(QLabel("Objective:"))
        self.objective_combo = QComboBox()
        self.objective_combo.setEnabled(False)
        self.objective_combo.setMinimumWidth(200)
        self.objective_combo.setToolTip("Hedef reaksiyonu seçin (boş: modelin varsayılan hedefi).")
        self.objective_combo.currentIndexChanged.connect(self.on_objective_changed)
        self.objective_combo.setEditable(True)
        self.objective_combo.setInsertPolicy(QComboBox.NoInsert)
        self.objective_combo.setMaxVisibleItems(25)
        self.objective_combo.lineEdit().setPlaceholderText("Type to search (id/description)...")
        self.objective_combo.lineEdit().textEdited.connect(self._on_objective_text_edited)
        controls.addWidget(self.objective_combo, stretch=2)

        controls.addStretch(1)
        outer.addLayout(controls)

        self.status_lbl = QLabel("No model loaded.")
        self.status_lbl.setWordWrap(True)
        outer.addWidget(self.status_lbl)

        self.summary_lbl = QLabel("")
        self.summary_lbl.setWordWrap(True)
        outer.addWidget(self.summary_lbl)

        self.tabs = QTabWidget()
        outer.addWidget(self.tabs, stretch=1)

        # -------- Tools tab (Installer) --------
        self.tab_tools = QWidget()
        tools_layout = QVBoxLayout(self.tab_tools)
        tools_layout.addWidget(self._title("Tools Installer (CarveMe + memote)"))

        self.tools_status_lbl = QLabel("Tools not checked yet.")
        self.tools_status_lbl.setWordWrap(True)
        tools_layout.addWidget(self.tools_status_lbl)

        btn_row = QHBoxLayout()
        self.tools_check_btn = QPushButton("Check tools")
        self.tools_check_btn.clicked.connect(self.tools_check_status)
        btn_row.addWidget(self.tools_check_btn)

        self.tools_install_btn = QPushButton("Install tools")
        self.tools_install_btn.clicked.connect(self.tools_install_tools)
        btn_row.addWidget(self.tools_install_btn)

        btn_row.addStretch(1)
        tools_layout.addLayout(btn_row)

        # (2) memote run + cancel
        memote_row = QHBoxLayout()
        self.memote_btn = QPushButton("Run memote report (current model)")
        self.memote_btn.setEnabled(False)
        self.memote_btn.clicked.connect(self.run_memote_report_current_model)
        memote_row.addWidget(self.memote_btn)

        self.memote_cancel_btn = QPushButton("Cancel memote")
        self.memote_cancel_btn.setEnabled(False)
        self.memote_cancel_btn.clicked.connect(self.cancel_memote)
        memote_row.addWidget(self.memote_cancel_btn)

        memote_row.addStretch(1)
        tools_layout.addLayout(memote_row)

        self.tools_log = QPlainTextEdit()
        self.tools_log.setReadOnly(True)
        tools_layout.addWidget(self.tools_log, stretch=1)

        self.tabs.addTab(self.tab_tools, "Tools")

        # -------- Medium tab --------
        self.tab_medium = QWidget()
        medium_layout = QVBoxLayout(self.tab_medium)

        medium_layout.addWidget(self._title("Medium / Exchange bounds"))
        self.medium_search = QLineEdit()
        self.medium_search.setPlaceholderText("Search exchange reactions (id or name)...")
        self.medium_search.textChanged.connect(self.filter_medium_table)
        medium_layout.addWidget(self.medium_search)

        self.medium_table = QTableWidget(0, 4)
        self.medium_table.setHorizontalHeaderLabels(["Reaction", "Name", "Lower bound", "Upper bound"])
        self.medium_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.medium_table.setEditTriggers(QAbstractItemView.DoubleClicked | QAbstractItemView.EditKeyPressed)
        self._configure_table_for_excel_resize(self.medium_table)
        medium_layout.addWidget(self.medium_table, stretch=1)

        medium_layout.addWidget(self._title("Cell details (full text)"))
        self.medium_details = QPlainTextEdit()
        self.medium_details.setReadOnly(True)
        self.medium_details.setMaximumHeight(140)
        medium_layout.addWidget(self.medium_details)
        self.medium_table.currentItemChanged.connect(self._update_medium_details_from_current_cell)

        # Presets row
        preset_row = QHBoxLayout()
        preset_row.addWidget(QLabel("Medium preset:"))
        self.medium_preset_combo = QComboBox()
        self.medium_preset_combo.addItems(["Reset to original", "Minimal M9", "Rich LB"])
        preset_row.addWidget(self.medium_preset_combo)
        self.medium_apply_preset_btn = QPushButton("Apply preset")
        self.medium_apply_preset_btn.setEnabled(False)
        self.medium_apply_preset_btn.clicked.connect(self.apply_medium_preset)
        preset_row.addWidget(self.medium_apply_preset_btn)
        preset_row.addStretch(1)
        medium_layout.addLayout(preset_row)

        self.tabs.addTab(self.tab_medium, "Medium")

        # -------- Reactions tab --------
        self.tab_rxns = QWidget()
        rxn_layout = QVBoxLayout(self.tab_rxns)

        rxn_layout.addWidget(self._title("All reactions (bounds & metadata)"))

        rxn_search_row = QHBoxLayout()
        rxn_search_row.addWidget(QLabel("Search (id/equation/GPR/KEGG/EC/description):"))
        self.rxn_global_search = QLineEdit()
        self.rxn_global_search.setPlaceholderText("e.g. biomass / atp / R01070 / 1.1.1.1 / gene ...")
        self.rxn_global_search.textChanged.connect(self.filter_reaction_table)
        rxn_search_row.addWidget(self.rxn_global_search, stretch=2)

        rxn_search_row.addWidget(QLabel("Metabolite filter:"))
        self.met_filter = QLineEdit()
        self.met_filter.setPlaceholderText("e.g. glc__D, atp (comma-separated)")
        self.met_filter.textChanged.connect(self.filter_reaction_table)
        rxn_search_row.addWidget(self.met_filter, stretch=2)

        rxn_search_row.addWidget(QLabel("Mode:"))
        self.met_filter_mode = QComboBox()
        self.met_filter_mode.addItems(["AND", "OR"])
        self.met_filter_mode.currentTextChanged.connect(self.filter_reaction_table)
        rxn_search_row.addWidget(self.met_filter_mode)

        rxn_layout.addLayout(rxn_search_row)

        self.reaction_table = QTableWidget(0, 8)
        self.reaction_table.setHorizontalHeaderLabels(
            ["Reaction", "Description", "Equation", "KEGG", "EC", "GPR", "LB", "UB"]
        )
        self.reaction_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.reaction_table.setEditTriggers(QAbstractItemView.DoubleClicked | QAbstractItemView.EditKeyPressed)
        self._configure_table_for_excel_resize(self.reaction_table)
        self.reaction_table.itemChanged.connect(self.on_reaction_table_item_changed)
        rxn_layout.addWidget(self.reaction_table, stretch=1)

        rxn_btn_row = QHBoxLayout()
        self.clear_rxn_overrides_btn = QPushButton("Clear reaction overrides")
        self.clear_rxn_overrides_btn.setEnabled(False)
        self.clear_rxn_overrides_btn.clicked.connect(self.clear_reaction_overrides)
        rxn_btn_row.addWidget(self.clear_rxn_overrides_btn)

        self.remove_selected_overrides_btn = QPushButton("Remove override (selected)")
        self.remove_selected_overrides_btn.setEnabled(False)
        self.remove_selected_overrides_btn.clicked.connect(self.remove_selected_reaction_overrides)
        rxn_btn_row.addWidget(self.remove_selected_overrides_btn)

        rxn_btn_row.addStretch(1)
        rxn_layout.addLayout(rxn_btn_row)

        rxn_layout.addWidget(self._title("Cell details (full text)"))
        self.reaction_details = QPlainTextEdit()
        self.reaction_details.setReadOnly(True)
        self.reaction_details.setMaximumHeight(180)
        rxn_layout.addWidget(self.reaction_details)
        self.reaction_table.currentItemChanged.connect(self._update_reaction_details_from_current_cell)
        # Double-click popups disabled to keep inline editing simple
        rxn_note = QLabel(
            "Tip: Drag column borders to resize (like Excel).\n"
            "KEGG/EC columns depend on SBML annotations; if the SBML doesn't include them, these will be empty."
        )
        rxn_note.setWordWrap(True)
        rxn_layout.addWidget(rxn_note)

        self.tabs.addTab(self.tab_rxns, "Reactions")

        # -------- Genes tab --------
        self.tab_genes = QWidget()
        genes_layout = QVBoxLayout(self.tab_genes)

        genes_layout.addWidget(self._title("Genes → Related reactions"))

        gene_search_row = QHBoxLayout()
        gene_search_row.addWidget(QLabel("Search gene:"))
        self.gene_global_search = QLineEdit()
        self.gene_global_search.setPlaceholderText("Type gene id or name...")
        self.gene_global_search.textChanged.connect(self.filter_genes_tab_gene_list)
        gene_search_row.addWidget(self.gene_global_search, stretch=1)
        genes_layout.addLayout(gene_search_row)

        genes_split = QHBoxLayout()

        left = QVBoxLayout()
        left.addWidget(QLabel("Genes"))
        self.genes_tab_list = QListWidget()
        self.genes_tab_list.setSelectionMode(QAbstractItemView.SingleSelection)
        self.genes_tab_list.currentItemChanged.connect(self.on_genes_tab_gene_changed)
        left.addWidget(self.genes_tab_list, stretch=1)
        genes_split.addLayout(left, stretch=1)

        right = QVBoxLayout()
        right.addWidget(QLabel("Reactions related to selected gene"))
        self.gene_rxn_table = QTableWidget(0, 4)
        self.gene_rxn_table.setHorizontalHeaderLabels(["Reaction", "Description", "Equation", "GPR"])
        self.gene_rxn_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.gene_rxn_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self._configure_table_for_excel_resize(self.gene_rxn_table)
        right.addWidget(self.gene_rxn_table, stretch=1)

        right.addWidget(self._title("Cell details (full text)"))
        self.gene_rxn_details = QPlainTextEdit()
        self.gene_rxn_details.setReadOnly(True)
        self.gene_rxn_details.setMaximumHeight(180)
        right.addWidget(self.gene_rxn_details)

        self.gene_rxn_table.currentItemChanged.connect(self._update_gene_rxn_details_from_current_cell)
        self.gene_rxn_table.itemDoubleClicked.connect(lambda item: self._open_cell_popup("Gene→Reaction cell", item))

        genes_split.addLayout(right, stretch=3)

        genes_layout.addLayout(genes_split, stretch=1)

        genes_note = QLabel("Note: Reactions are determined from COBRApy gene↔reaction associations (from GPRs).")
        genes_note.setWordWrap(True)
        genes_layout.addWidget(genes_note)

        self.tabs.addTab(self.tab_genes, "Genes")

        # -------- Model Editor tab placeholder --------
        self.tab_editor = QWidget()
        editor_layout = QVBoxLayout(self.tab_editor)
        editor_layout.addWidget(self._title("Model Editor (Patch-based)"))
        self.editor_tabs = QTabWidget()
        editor_layout.addWidget(self.editor_tabs, stretch=1)

        self._build_editor_metabolite_tab()
        self._build_editor_gene_tab()
        self._build_editor_reaction_tab()
        self._build_editor_disable_delete_tab()
        self._build_editor_patch_tab()
        self._build_editor_export_sbml_tab()

        self.tabs.addTab(self.tab_editor, "Model Editor")

        # -------- Gene knockout tab --------
        self.tab_ko = QWidget()
        ko_layout = QVBoxLayout(self.tab_ko)

        ko_layout.addWidget(self._title("Gene knockout"))
        self.gene_search = QLineEdit()
        self.gene_search.setPlaceholderText("Search genes...")
        self.gene_search.textChanged.connect(self.filter_gene_list)
        ko_layout.addWidget(self.gene_search)

        self.gene_list = QListWidget()
        # Allow selecting multiple genes so double deletion can pick exact pairs
        self.gene_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        ko_layout.addWidget(self.gene_list, stretch=1)

        ko_btn_row = QHBoxLayout()
        self.add_ko_btn = QPushButton("Add knockout")
        self.add_ko_btn.setEnabled(False)
        self.add_ko_btn.clicked.connect(self.add_selected_gene_knockout)

        self.clear_ko_btn = QPushButton("Clear knockouts")
        self.clear_ko_btn.setEnabled(False)
        self.clear_ko_btn.clicked.connect(self.clear_knockouts)

        ko_btn_row.addWidget(self.add_ko_btn)
        ko_btn_row.addWidget(self.clear_ko_btn)
        ko_btn_row.addStretch(1)
        ko_layout.addLayout(ko_btn_row)

        self.ko_lbl = QLabel("Gene knockouts: (none)")
        self.ko_lbl.setWordWrap(True)
        ko_layout.addWidget(self.ko_lbl)
        self.tabs.addTab(self.tab_ko, "Gene knockout")

        # -------- Overexpression tab --------
        self.tab_oe = QWidget()
        oe_layout = QVBoxLayout(self.tab_oe)

        oe_layout.addWidget(self._title("Overexpression (reaction bound scaling)"))

        self.rxn_search = QLineEdit()
        self.rxn_search.setPlaceholderText("Search reactions for overexpression (id or description)...")
        self.rxn_search.textChanged.connect(self.filter_overexpression_reaction_list)
        oe_layout.addWidget(self.rxn_search)

        self.rxn_list = QListWidget()
        self.rxn_list.setSelectionMode(QAbstractItemView.SingleSelection)
        self.rxn_list.currentItemChanged.connect(self.update_selected_rxn_info)
        oe_layout.addWidget(self.rxn_list, stretch=2)

        self.rxn_info_lbl = QLabel("Selected reaction: (none)")
        self.rxn_info_lbl.setWordWrap(True)
        oe_layout.addWidget(self.rxn_info_lbl)

        self.rxn_preview_lbl = QLabel("Preview: -")
        self.rxn_preview_lbl.setWordWrap(True)
        oe_layout.addWidget(self.rxn_preview_lbl)

        oe_row = QHBoxLayout()
        oe_row.addWidget(QLabel("Overexpression factor:"))
        self.oe_factor = QDoubleSpinBox()
        self.oe_factor.setRange(1.0, 100.0)
        self.oe_factor.setDecimals(2)
        self.oe_factor.setSingleStep(0.25)
        self.oe_factor.setValue(2.0)
        self.oe_factor.valueChanged.connect(self.update_selected_rxn_info)
        oe_row.addWidget(self.oe_factor)

        self.oe_scale_lb = QCheckBox("Scale LB too (if LB < 0)")
        self.oe_scale_lb.setChecked(True)
        self.oe_scale_lb.stateChanged.connect(self.update_selected_rxn_info)
        oe_row.addWidget(self.oe_scale_lb)

        self.add_oe_btn = QPushButton("Add overexpression")
        self.add_oe_btn.setEnabled(False)
        self.add_oe_btn.clicked.connect(self.add_selected_reaction_overexpression)
        oe_row.addWidget(self.add_oe_btn)

        self.remove_oe_btn = QPushButton("Remove overexpression (selected)")
        self.remove_oe_btn.setEnabled(False)
        self.remove_oe_btn.clicked.connect(self.remove_selected_reaction_overexpression)
        oe_row.addWidget(self.remove_oe_btn)

        self.clear_oe_btn = QPushButton("Clear overexpressions")
        self.clear_oe_btn.setEnabled(False)
        self.clear_oe_btn.clicked.connect(self.clear_overexpressions)
        oe_row.addWidget(self.clear_oe_btn)

        oe_row.addStretch(1)
        oe_layout.addLayout(oe_row)

        temp_row = QHBoxLayout()
        temp_row.addWidget(QLabel("Temporary upper bound (perturbed only):"))
        self.temp_ub_spin = QDoubleSpinBox()
        self.temp_ub_spin.setRange(-1e9, 1e9)
        self.temp_ub_spin.setDecimals(6)
        self.temp_ub_spin.setSingleStep(1.0)
        self.temp_ub_spin.setValue(0.0)
        temp_row.addWidget(self.temp_ub_spin)

        self.set_temp_ub_btn = QPushButton("Set temp upper bound (selected)")
        self.set_temp_ub_btn.setEnabled(False)
        self.set_temp_ub_btn.clicked.connect(self.set_temp_upper_bound_for_selected)
        temp_row.addWidget(self.set_temp_ub_btn)

        self.remove_temp_ub_btn = QPushButton("Remove temp upper bound (selected)")
        self.remove_temp_ub_btn.setEnabled(False)
        self.remove_temp_ub_btn.clicked.connect(self.remove_temp_upper_bound_for_selected)
        temp_row.addWidget(self.remove_temp_ub_btn)

        self.clear_temp_ub_btn = QPushButton("Clear all temp upper bounds")
        self.clear_temp_ub_btn.setEnabled(False)
        self.clear_temp_ub_btn.clicked.connect(self.clear_temp_upper_bounds)
        temp_row.addWidget(self.clear_temp_ub_btn)

        temp_row.addStretch(1)
        oe_layout.addLayout(temp_row)

        oe_layout.addWidget(self._title("Active overexpressions"))
        self.oe_active_list = QListWidget()
        self.oe_active_list.itemClicked.connect(self.jump_to_reaction_from_item)
        oe_layout.addWidget(self.oe_active_list, stretch=1)

        oe_layout.addWidget(self._title("Temp upper bounds"))
        self.temp_active_list = QListWidget()
        self.temp_active_list.itemClicked.connect(self.jump_to_reaction_from_item)
        oe_layout.addWidget(self.temp_active_list, stretch=1)

        self.tabs.addTab(self.tab_oe, "Overexpression")

        # -------- FVA tab --------
        self.tab_fva = QWidget()
        fva_layout = QVBoxLayout(self.tab_fva)

        fva_layout.addWidget(self._title("Flux Variability Analysis (FVA)"))
        fva_help = QLabel(
            "Shows minimum and maximum feasible flux for each reaction while maintaining optimal growth. "
            "Reactions with wide ranges are flexible; those with narrow ranges are constrained."
        )
        fva_help.setStyleSheet("font-size:9px;color:#555;margin:4px;")
        fva_help.setWordWrap(True)
        fva_layout.addWidget(fva_help)

        fva_filter_row = QHBoxLayout()
        fva_filter_row.addWidget(QLabel("Reaction search:"))
        self.fva_flux_search = QLineEdit()
        self.fva_flux_search.setPlaceholderText("reaction id contains...")
        self.fva_flux_search.textChanged.connect(self.filter_fva_table)
        fva_filter_row.addWidget(self.fva_flux_search, stretch=2)
        fva_filter_row.addStretch(1)
        fva_layout.addLayout(fva_filter_row)

        self.fva_canvas = MplCanvas(self, width=10, height=3.8, dpi=100)
        fva_layout.addWidget(self.fva_canvas, stretch=1)

        self.fva_tbl = QTableWidget(0, 4)
        self.fva_tbl.setHorizontalHeaderLabels(["Reaction", "Minimum flux", "Maximum flux", "Range"])
        self.fva_tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.fva_tbl.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._configure_table_for_excel_resize(self.fva_tbl)
        fva_layout.addWidget(self.fva_tbl, stretch=1)

        self.fva_info_lbl = QLabel("FVA results: (none)")
        self.fva_info_lbl.setWordWrap(True)
        fva_layout.addWidget(self.fva_info_lbl)

        fva_btn_row = QHBoxLayout()
        self.fva_export_btn = QPushButton("Export FVA (CSV)")
        self.fva_export_btn.clicked.connect(self.export_fva_csv)
        fva_btn_row.addWidget(self.fva_export_btn)
        fva_btn_row.addStretch(1)
        fva_layout.addLayout(fva_btn_row)

        self.tabs.addTab(self.tab_fva, "FVA")

        # -------- pFBA tab --------
        self.tab_pfba = QWidget()
        pfba_layout = QVBoxLayout(self.tab_pfba)

        pfba_layout.addWidget(self._title("Parsimonious FBA (pFBA)"))
        pfba_help = QLabel(
            "After finding the maximum objective value (like FBA), minimizes total flux. "
            "Results in sparse, efficient flux distributions with fewer active reactions."
        )
        pfba_layout.addWidget(pfba_help)

        pfba_filter_row = QHBoxLayout()
        pfba_filter_row.addWidget(QLabel("Reaction search:"))
        self.pfba_flux_search = QLineEdit()
        self.pfba_flux_search.setPlaceholderText("reaction id contains...")
        self.pfba_flux_search.textChanged.connect(self.filter_pfba_table)
        pfba_filter_row.addWidget(self.pfba_flux_search, stretch=2)
        pfba_filter_row.addStretch(1)
        pfba_layout.addLayout(pfba_filter_row)

        self.pfba_canvas = MplCanvas(self, width=10, height=3.8, dpi=100)
        pfba_layout.addWidget(self.pfba_canvas, stretch=1)

        self.pfba_tbl = QTableWidget(0, 4)
        self.pfba_tbl.setHorizontalHeaderLabels(["Reaction", "FBA flux", "pFBA flux", "Difference"])
        self.pfba_tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.pfba_tbl.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._configure_table_for_excel_resize(self.pfba_tbl)
        pfba_layout.addWidget(self.pfba_tbl, stretch=1)

        self.pfba_info_lbl = QLabel("pFBA results: (none)")
        self.pfba_info_lbl.setWordWrap(True)
        pfba_layout.addWidget(self.pfba_info_lbl)

        pfba_btn_row = QHBoxLayout()
        self.pfba_export_btn = QPushButton("Export pFBA vs FBA (CSV)")
        self.pfba_export_btn.clicked.connect(self.export_pfba_csv)
        pfba_btn_row.addWidget(self.pfba_export_btn)
        pfba_btn_row.addStretch(1)
        pfba_layout.addLayout(pfba_btn_row)

        self.tabs.addTab(self.tab_pfba, "pFBA")

        # -------- SGD/SRD tab --------
        self.tab_deletion = QWidget()
        deletion_layout = QVBoxLayout(self.tab_deletion)

        deletion_layout.addWidget(self._title("Single Gene/Reaction Deletion Analysis"))
        deletion_help = QLabel(
            "Knocks out each gene (SGD) or reaction (SRD) one at a time and measures growth impact. "
            "Essential genes/reactions cause severe growth reduction or lethality when deleted."
        )
        deletion_help.setWordWrap(True)
        deletion_layout.addWidget(deletion_help)

        deletion_search_row = QHBoxLayout()
        deletion_search_row.addWidget(QLabel("Search:"))
        self.deletion_search = QLineEdit()
        self.deletion_search.setPlaceholderText("gene/reaction id contains...")
        self.deletion_search.textChanged.connect(self.filter_deletion_table)
        deletion_search_row.addWidget(self.deletion_search, stretch=2)
        deletion_search_row.addStretch(1)
        deletion_layout.addLayout(deletion_search_row)

        self.deletion_canvas = MplCanvas(self, width=10, height=3.8, dpi=100)
        deletion_layout.addWidget(self.deletion_canvas, stretch=1)

        self.deletion_tbl = QTableWidget(0, 5)
        self.deletion_tbl.setHorizontalHeaderLabels(["Item ID", "WT Growth", "Growth on Deletion", "Δ Growth (%)", "Category"])
        self.deletion_tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.deletion_tbl.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._configure_table_for_excel_resize(self.deletion_tbl)
        try:
            header = self.deletion_tbl.horizontalHeader()
            header.setSectionResizeMode(QHeaderView.Stretch)
            header.setStretchLastSection(True)
        except Exception:
            pass
        deletion_layout.addWidget(self.deletion_tbl, stretch=1)

        self.deletion_info_lbl = QLabel("Deletion analysis: (none)")
        self.deletion_info_lbl.setWordWrap(True)
        deletion_layout.addWidget(self.deletion_info_lbl)

        del_btn_row = QHBoxLayout()
        self.sgd_export_btn = QPushButton("Export SGD (CSV)")
        self.sgd_export_btn.clicked.connect(self.export_sgd_csv)
        del_btn_row.addWidget(self.sgd_export_btn)
        self.srd_export_btn = QPushButton("Export SRD (CSV)")
        self.srd_export_btn.clicked.connect(self.export_srd_csv)
        del_btn_row.addWidget(self.srd_export_btn)
        del_btn_row.addStretch(1)
        deletion_layout.addLayout(del_btn_row)

        self.tabs.addTab(self.tab_deletion, "Deletion Analysis")

        # -------- Robustness tab --------
        self.tab_robustness = QWidget()
        robustness_layout = QVBoxLayout(self.tab_robustness)

        robustness_layout.addWidget(self._title("Robustness Analysis"))
        robustness_help = QLabel(
            "Varies the bounds of a specific reaction (e.g., glucose uptake) over a range of values "
            "and tracks how the objective (growth rate) changes. Shows sensitivity and trade-offs."
        )
        robustness_help.setWordWrap(True)
        robustness_layout.addWidget(robustness_help)

        robustness_param_row = QHBoxLayout()
        robustness_param_row.addWidget(QLabel("Reaction:"))
        self.robustness_rxn = QLineEdit()
        self.robustness_rxn.setPlaceholderText("e.g. EX_glc__D_e")
        robustness_param_row.addWidget(self.robustness_rxn, stretch=1)
        robustness_param_row.addWidget(QLabel("Min:"))
        self.robustness_min = QDoubleSpinBox()
        self.robustness_min.setRange(-1000, 0)
        self.robustness_min.setValue(-20)
        robustness_param_row.addWidget(self.robustness_min)
        robustness_param_row.addWidget(QLabel("Max:"))
        self.robustness_max = QDoubleSpinBox()
        self.robustness_max.setRange(0, 1000)
        self.robustness_max.setValue(20)
        robustness_param_row.addWidget(self.robustness_max)
        robustness_param_row.addWidget(QLabel("Steps:"))
        self.robustness_steps = QSpinBox()
        self.robustness_steps.setRange(5, 100)
        self.robustness_steps.setValue(20)
        robustness_param_row.addWidget(self.robustness_steps)
        robustness_layout.addLayout(robustness_param_row)

        robustness_bound_row = QHBoxLayout()
        robustness_bound_row.addWidget(QLabel("Bound Type:"))
        self.robustness_bound_ub = QCheckBox("Upper Bound (UB)")
        self.robustness_bound_ub.setChecked(True)
        self.robustness_bound_ub.setToolTip("Vary the reaction's upper bound")
        robustness_bound_row.addWidget(self.robustness_bound_ub)
        self.robustness_bound_lb = QCheckBox("Lower Bound (LB)")
        self.robustness_bound_lb.setToolTip("Vary the reaction's lower bound")
        robustness_bound_row.addWidget(self.robustness_bound_lb)
        robustness_bound_row.addStretch(1)
        robustness_layout.addLayout(robustness_bound_row)

        robustness_search_row = QHBoxLayout()
        robustness_search_row.addWidget(QLabel("Search:"))
        self.robustness_search = QLineEdit()
        self.robustness_search.setPlaceholderText("reaction id contains...")
        self.robustness_search.textChanged.connect(self._filter_robustness_combo)
        robustness_search_row.addWidget(self.robustness_search, stretch=2)
        self.robustness_combo = QComboBox()
        self.robustness_combo.currentTextChanged.connect(lambda s: self.robustness_rxn.setText(s))
        robustness_search_row.addWidget(self.robustness_combo, stretch=2)
        robustness_search_row.addStretch(1)
        robustness_layout.addLayout(robustness_search_row)

        self.robustness_canvas = MplCanvas(self, width=10, height=5, dpi=100)
        robustness_layout.addWidget(self.robustness_canvas, stretch=1)

        self.robustness_info_lbl = QLabel("Robustness analysis: (none)")
        self.robustness_info_lbl.setWordWrap(True)
        robustness_layout.addWidget(self.robustness_info_lbl)

        self.tabs.addTab(self.tab_robustness, "Robustness")

        # -------- Production Envelope tab --------
        self.tab_envelope = QWidget()
        envelope_layout = QVBoxLayout(self.tab_envelope)

        envelope_layout.addWidget(self._title("Production Envelope (Growth vs Product)"))
        envelope_help = QLabel(
            "Shows the trade-off between biomass (growth) and production of a target compound. "
            "Useful for metabolic engineering to identify the best balance between growth and yield."
        )
        envelope_help.setWordWrap(True)
        envelope_layout.addWidget(envelope_help)

        envelope_param_row = QHBoxLayout()
        envelope_param_row.addWidget(QLabel("Product Reaction:"))
        self.envelope_product = QLineEdit()
        self.envelope_product.setPlaceholderText("e.g. ACALD (acetaldehyde production)")
        envelope_param_row.addWidget(self.envelope_product, stretch=1)
        envelope_param_row.addWidget(QLabel("Steps:"))
        self.envelope_steps = QSpinBox()
        self.envelope_steps.setRange(10, 100)
        self.envelope_steps.setValue(30)
        envelope_param_row.addWidget(self.envelope_steps)
        envelope_layout.addLayout(envelope_param_row)

        envelope_search_row = QHBoxLayout()
        envelope_search_row.addWidget(QLabel("Search:"))
        self.envelope_search = QLineEdit()
        self.envelope_search.setPlaceholderText("reaction id contains...")
        self.envelope_search.textChanged.connect(self._filter_envelope_combo)
        envelope_search_row.addWidget(self.envelope_search, stretch=2)
        self.envelope_combo = QComboBox()
        self.envelope_combo.currentTextChanged.connect(lambda s: self.envelope_product.setText(s))
        envelope_search_row.addWidget(self.envelope_combo, stretch=2)
        envelope_search_row.addStretch(1)
        envelope_layout.addLayout(envelope_search_row)

        self.envelope_canvas = MplCanvas(self, width=10, height=5, dpi=100)
        envelope_layout.addWidget(self.envelope_canvas, stretch=1)

        self.envelope_info_lbl = QLabel("Production envelope: (none)")
        self.envelope_info_lbl.setWordWrap(True)
        envelope_layout.addWidget(self.envelope_info_lbl)

        self.tabs.addTab(self.tab_envelope, "Envelope")

        # -------- Flux Sampling tab --------
        self.tab_sampling = QWidget()
        sampling_layout = QVBoxLayout(self.tab_sampling)

        sampling_layout.addWidget(self._title("Flux Sampling (ACHR)"))
        sampling_help = QLabel(
            "Monte Carlo sampling of the feasible flux space using the ACHR algorithm. "
            "Returns mean, standard deviation, and min/max ranges for each reaction—useful for probabilistic analysis."
        )
        sampling_help.setWordWrap(True)
        sampling_layout.addWidget(sampling_help)

        sampling_param_row = QHBoxLayout()
        sampling_param_row.addWidget(QLabel("Sample Size:"))
        self.sampling_size = QSpinBox()
        self.sampling_size.setRange(10, 10000)
        self.sampling_size.setValue(500)
        sampling_param_row.addWidget(self.sampling_size)
        sampling_param_row.addStretch(1)
        sampling_layout.addLayout(sampling_param_row)

        self.sampling_canvas = MplCanvas(self, width=10, height=5, dpi=100)
        sampling_layout.addWidget(self.sampling_canvas, stretch=1)

        sampling_search_row = QHBoxLayout()
        sampling_search_row.addWidget(QLabel("Search:"))
        self.sampling_search = QLineEdit()
        self.sampling_search.setPlaceholderText("reaction id contains...")
        self.sampling_search.textChanged.connect(self.filter_sampling_table)
        sampling_search_row.addWidget(self.sampling_search, stretch=2)
        sampling_search_row.addStretch(1)
        sampling_layout.addLayout(sampling_search_row)

        self.sampling_tbl = QTableWidget(0, 4)
        self.sampling_tbl.setHorizontalHeaderLabels(["Reaction", "Mean Flux", "Stdev", "Min/Max"])
        self.sampling_tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.sampling_tbl.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._configure_table_for_excel_resize(self.sampling_tbl)
        sampling_layout.addWidget(self.sampling_tbl, stretch=1)

        self.sampling_info_lbl = QLabel("Flux sampling: (none)")
        self.sampling_info_lbl.setWordWrap(True)
        sampling_layout.addWidget(self.sampling_info_lbl)

        sampling_btn_row = QHBoxLayout()
        self.sampling_export_btn = QPushButton("Export Sampling (CSV)")
        self.sampling_export_btn.clicked.connect(self.export_sampling_csv)
        sampling_btn_row.addWidget(self.sampling_export_btn)
        sampling_btn_row.addStretch(1)
        sampling_layout.addLayout(sampling_btn_row)

        self.tabs.addTab(self.tab_sampling, "Sampling")

        # -------- Network Map tab --------
        self.tab_map = QWidget()
        map_layout = QVBoxLayout(self.tab_map)

        map_layout.addWidget(self._title("Network Map"))
        map_help = QLabel(
            "Visualize metabolic network topology. Nodes: reactions (circles) and metabolites (squares). "
            "Colors reflect flux magnitude or gene knockout impact. Filter by reaction/metabolite search and neighborhood depth."
        )
        map_help.setWordWrap(True)
        map_layout.addWidget(map_help)

        # Search & filter row
        map_filter_row = QHBoxLayout()
        map_filter_row.addWidget(QLabel("Search:"))
        self.map_search = QLineEdit()
        self.map_search.setPlaceholderText("reaction/metabolite id contains...")
        self.map_search.textChanged.connect(self._render_network_map)
        map_filter_row.addWidget(self.map_search, stretch=2)

        map_filter_row.addWidget(QLabel("Depth:"))
        self.map_depth = QSpinBox()
        self.map_depth.setRange(1, 5)
        self.map_depth.setValue(2)
        self.map_depth.setToolTip("Neighborhood depth (1-5)")
        self.map_depth.valueChanged.connect(self._render_network_map)
        map_filter_row.addWidget(self.map_depth)

        map_filter_row.addWidget(QLabel("Flux threshold:"))
        self.map_threshold = QDoubleSpinBox()
        self.map_threshold.setRange(0, 100)
        self.map_threshold.setValue(0.1)
        self.map_threshold.setSingleStep(0.1)
        self.map_threshold.setToolTip("Show only reactions with |flux| > threshold")
        self.map_threshold.valueChanged.connect(self._render_network_map)
        map_filter_row.addWidget(self.map_threshold)

        map_filter_row.addWidget(QLabel("View:"))
        self.map_view_combo = QComboBox()
        self.map_view_combo.addItems(["Flux magnitude", "Gene KO Impact", "FVA Bounds"])
        self.map_view_combo.currentTextChanged.connect(self._render_network_map)
        map_filter_row.addWidget(self.map_view_combo)

        map_filter_row.addStretch(1)
        map_layout.addLayout(map_filter_row)

        self.map_canvas = MplCanvas(self, width=12, height=7, dpi=100)
        map_layout.addWidget(self.map_canvas, stretch=1)
        self.map_canvas.mpl_connect('button_press_event', self._on_map_click)

        # Zoom/Pan toolbar
        self.map_toolbar = NavigationToolbar(self.map_canvas, self.tab_map)
        map_layout.addWidget(self.map_toolbar)

        self.map_info_lbl = QLabel("Network map: (none)")
        self.map_info_lbl.setWordWrap(True)
        map_layout.addWidget(self.map_info_lbl)

        self.map_details = QPlainTextEdit()
        self.map_details.setReadOnly(True)
        self.map_details.setMaximumHeight(180)
        self.map_details.setPlaceholderText("Click on a node in the map to see reaction details...")
        map_layout.addWidget(self.map_details)

        self.tabs.addTab(self.tab_map, "Map")

        # -------- Results tab --------
        self.tab_results = QWidget()
        results_layout = QVBoxLayout(self.tab_results)

        results_layout.addWidget(self._title("Results"))

        compare_row = QHBoxLayout()
        compare_row.addWidget(QLabel("Compare:"))
        self.compare_mode = QComboBox()
        self.compare_mode.addItems(["Gene knockout only", "Overexpression only", "Original vs Baseline"])
        self.compare_mode.currentTextChanged.connect(self._refresh_results_view_from_last_run)
        compare_row.addWidget(self.compare_mode)
        compare_row.addStretch(1)
        results_layout.addLayout(compare_row)

        self.results_lbl = QLabel("Objective: -")
        self.results_lbl.setWordWrap(True)
        results_layout.addWidget(self.results_lbl)

        # (4) Flux table filter controls
        flux_filter_row = QHBoxLayout()
        flux_filter_row.addWidget(QLabel("Flux search:"))
        self.flux_search = QLineEdit()
        self.flux_search.setPlaceholderText("reaction id contains...")
        self.flux_search.textChanged.connect(self.filter_flux_table)
        flux_filter_row.addWidget(self.flux_search, stretch=2)

        flux_filter_row.addWidget(QLabel("|delta| ≥"))
        self.delta_thresh = QDoubleSpinBox()
        self.delta_thresh.setRange(0.0, 1e12)
        self.delta_thresh.setDecimals(9)
        self.delta_thresh.setSingleStep(0.01)
        self.delta_thresh.setValue(0.0)
        self.delta_thresh.valueChanged.connect(self.filter_flux_table)
        flux_filter_row.addWidget(self.delta_thresh)

        flux_filter_row.addStretch(1)
        results_layout.addLayout(flux_filter_row)

        self.canvas = MplCanvas(self, width=10, height=3.8, dpi=100)
        results_layout.addWidget(self.canvas, stretch=1)

        self.flux_tbl = QTableWidget(0, 4)
        self.flux_tbl.setHorizontalHeaderLabels(["Reaction", "Baseline flux", "Compared flux", "Delta"])
        self.flux_tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.flux_tbl.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._configure_table_for_excel_resize(self.flux_tbl)
        results_layout.addWidget(self.flux_tbl, stretch=1)

        self.tabs.addTab(self.tab_results, "Results")

        self.setCentralWidget(root)
        self.resize(1700, 1120)
        self.statusBar().showMessage("Ready.")
        self._clear_chart()
        # Apply dark theme
        self.setStyleSheet(self._dark_stylesheet())
        # Apply persisted UI settings before building menu
        self._init_settings()
        self._build_menu()

        # Re-enable repaints after constructing UI
        self.setUpdatesEnabled(True)

    # ---------------- (1) Unsaved changes prompt on exit ----------------
    def closeEvent(self, event):
        if self.base_model is None or not self.model_dirty:
            event.accept()
            return

        box = QMessageBox(self)
        box.setIcon(QMessageBox.Warning)
        box.setWindowTitle("Unsaved changes")
        box.setText("You have unsaved changes in the current model.\n\nSave as SBML before exiting?")

        btn_save = box.addButton("Save As...", QMessageBox.AcceptRole)
        btn_discard = box.addButton("Don't Save", QMessageBox.DestructiveRole)
        btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
        box.setDefaultButton(btn_save)

        box.exec()
        clicked = box.clickedButton()

        if clicked == btn_cancel:
            event.ignore()
            return
        if clicked == btn_save:
            saved = self.save_current_model_as()
            if saved is None:
                event.ignore()
                return
            event.accept()
            return

        event.accept()

    # ---------------- UI helpers ----------------
    def _title(self, text: str) -> QLabel:
        lbl = QLabel(text)
        lbl.setStyleSheet("font-weight: 600;")
        return lbl

    def _configure_table_for_excel_resize(self, table: QTableWidget):
        hdr = table.horizontalHeader()
        hdr.setSectionResizeMode(QHeaderView.Interactive)
        hdr.setStretchLastSection(False)
        hdr.setMinimumSectionSize(40)
        table.setWordWrap(False)
        table.setHorizontalScrollMode(QAbstractItemView.ScrollPerPixel)
        table.setVerticalScrollMode(QAbstractItemView.ScrollPerPixel)
        table.verticalHeader().setVisible(False)

    # ---------------- Theme & Settings ----------------
    def _init_settings(self):
        self.settings = QSettings("MetaboDesk", "MetaboDesk")
        # Disable custom themes; use default Qt styling
        try:
            self.setStyleSheet("")
        except Exception:
            pass

    def _init_persisted_controls(self):
        """Restore TopN and Solver from settings."""
        try:
            topn = int(self.settings.value("ui/topn", self.topn_spin.value()))
            self.topn_spin.setValue(topn)
        except Exception:
            pass
        try:
            solver = str(self.settings.value("ui/solver", "Auto") or "Auto")
            idx = self.solver_combo.findText(solver)
            if idx >= 0:
                self.solver_combo.setCurrentIndex(idx)
        except Exception:
            pass

    # Placeholder - theme functions removed

    def _light_stylesheet(self) -> str:
        return (
            "QMainWindow{background:#f5f6f8;color:#2c3e50;font-family:'Segoe UI','Ubuntu','Helvetica';font-size:10pt;}\n"
            "QLabel{color:#2c3e50;}\n"
            "QLineEdit,QPlainTextEdit,QComboBox,QSpinBox,QDoubleSpinBox{background:#ffffff;color:#2c3e50;border:1px solid #cbd5e0;border-radius:4px;padding:6px;font-size:10pt;}\n"
            "QLineEdit:focus,QPlainTextEdit:focus,QComboBox:focus,QSpinBox:focus,QDoubleSpinBox:focus{border:2px solid #3498db;background:#fafbfc;}\n"
            "QPushButton{background:#3498db;color:#ffffff;border:none;border-radius:4px;padding:8px 14px;font-weight:600;font-size:10pt;}\n"
            "QPushButton:hover{background:#2980b9;}\n"
            "QPushButton:pressed{background:#1f618d;}\n"
            "QPushButton:disabled{background:#bdc3c7;color:#e0e0e0;}\n"
            "QTableWidget{background:#ffffff;color:#2c3e50;gridline-color:#e8eaed;border:1px solid #cbd5e0;}\n"
            "QTableWidget::item{padding:3px;}\n"
            "QHeaderView::section{background:#e8eaed;color:#2c3e50;padding:8px;border:none;border-right:1px solid #cbd5e0;font-weight:600;}\n"
            "QTabBar::tab{background:#e8eaed;color:#2c3e50;padding:10px 16px;border:none;border-bottom:3px solid #cbd5e0;margin-right:4px;font-weight:500;}\n"
            "QTabBar::tab:selected{background:#ffffff;border-bottom:3px solid #3498db;color:#2c3e50;}\n"
            "QTabBar::tab:hover{background:#dfe4e8;}\n"
            "QStatusBar{background:#e8eaed;color:#2c3e50;border-top:1px solid #cbd5e0;font-size:9pt;}\n"
            "QGroupBox{color:#2c3e50;border:1px solid #cbd5e0;border-radius:4px;padding:10px;margin-top:8px;font-weight:600;}\n"
            "QGroupBox::title{subcontrol-origin:margin;left:10px;padding:0 3px;}\n"
            "QToolTip{background:#2c3e50;color:#ecf0f1;border:1px solid #34495e;border-radius:3px;padding:4px 8px;font-size:9pt;}\n"
            "QScrollBar:vertical{background:#f5f6f8;width:12px;}\n"
            "QScrollBar::handle:vertical{background:#cbd5e0;border-radius:6px;min-height:20px;}\n"
            "QScrollBar::handle:vertical:hover{background:#a8b3bf;}"
        )

    def _dark_stylesheet(self) -> str:
        return (
            "QMainWindow{background:#0d1117;color:#e6edf3;font-family:'Segoe UI','Ubuntu','Helvetica';font-size:10pt;}\n"
            "QLabel{color:#e6edf3;}\n"
            "QLineEdit,QPlainTextEdit,QComboBox,QSpinBox,QDoubleSpinBox{background:#161b22;color:#e6edf3;border:1px solid #30363d;border-radius:4px;padding:6px;font-size:10pt;}\n"
            "QLineEdit:focus,QPlainTextEdit:focus,QComboBox:focus,QSpinBox:focus,QDoubleSpinBox:focus{border:2px solid #58a6ff;background:#0d1117;}\n"
            "QPushButton{background:#238636;color:#ffffff;border:none;border-radius:4px;padding:8px 14px;font-weight:600;font-size:10pt;}\n"
            "QPushButton:hover{background:#2ea043;}\n"
            "QPushButton:pressed{background:#1f6feb;}\n"
            "QPushButton:disabled{background:#30363d;color:#6e7681;}\n"
            "QTableWidget{background:#161b22;color:#e6edf3;gridline-color:#30363d;border:1px solid #30363d;}\n"
            "QTableWidget::item{padding:3px;color:#e6edf3;}\n"
            "QHeaderView::section{background:#0d1117;color:#8b949e;padding:8px;border:none;border-right:1px solid #30363d;font-weight:600;}\n"
            "QTabBar::tab{background:#0d1117;color:#8b949e;padding:10px 16px;border:none;border-bottom:3px solid #30363d;margin-right:4px;font-weight:500;}\n"
            "QTabBar::tab:selected{background:#161b22;border-bottom:3px solid #58a6ff;color:#e6edf3;}\n"
            "QTabBar::tab:hover{background:#21262d;}\n"
            "QStatusBar{background:#0d1117;color:#8b949e;border-top:1px solid #30363d;font-size:9pt;}\n"
            "QGroupBox{color:#e6edf3;border:1px solid #30363d;border-radius:4px;padding:10px;margin-top:8px;font-weight:600;}\n"
            "QGroupBox::title{subcontrol-origin:margin;left:10px;padding:0 3px;color:#e6edf3;}\n"
            "QToolTip{background:#161b22;color:#e6edf3;border:1px solid #30363d;border-radius:3px;padding:4px 8px;font-size:9pt;}\n"
            "QScrollBar:vertical{background:#0d1117;width:12px;}\n"
            "QScrollBar::handle:vertical{background:#30363d;border-radius:6px;min-height:20px;}\n"
            "QScrollBar::handle:vertical:hover{background:#6e7681;}"
        )

    def _open_cell_popup(self, title: str, item: QTableWidgetItem | None):
        if not item:
            return
        dlg = TextPopup(title, item.text(), self)
        dlg.exec()

    def _on_reaction_cell_double_clicked(self, item: QTableWidgetItem | None):
        if not item:
            return
        col = item.column()
        text = (item.text() or "").strip()
        # KEGG column index=3, EC=4
        if col == 3 and text:
            kid = text.split(";")[0].strip()
            if kid:
                url = f"https://www.kegg.jp/dbget-bin/www_bget?rn:{kid}"
                QDesktopServices.openUrl(QUrl(url))
                return
        if col == 4 and text:
            ec = text.split(";")[0].strip()
            if ec:
                url = f"https://www.ebi.ac.uk/intenz/query?cmd=searchEC&ec={ec}"
                QDesktopServices.openUrl(QUrl(url))
                return
        # Fallback: open popup
        self._open_cell_popup("Reaction cell", item)

    def _update_medium_details_from_current_cell(self, current: QTableWidgetItem, _prev: QTableWidgetItem):
        self.medium_details.setPlainText(current.text() if current else "")

    def _update_reaction_details_from_current_cell(self, current: QTableWidgetItem, _prev: QTableWidgetItem):
        self.reaction_details.setPlainText(current.text() if current else "")

    def _update_gene_rxn_details_from_current_cell(self, current: QTableWidgetItem, _prev: QTableWidgetItem):
        self.gene_rxn_details.setPlainText(current.text() if current else "")

    # ---------------- (4) Flux table filtering ----------------
    def filter_flux_table(self):
        text = (self.flux_search.text() or "").strip().lower()
        thr = float(self.delta_thresh.value())
        for row in range(self.flux_tbl.rowCount()):
            rid = (self.flux_tbl.item(row, 0).text() if self.flux_tbl.item(row, 0) else "").lower()
            try:
                d = float((self.flux_tbl.item(row, 3).text() if self.flux_tbl.item(row, 3) else "0").strip())
            except Exception:
                d = 0.0

            hide = False
            if text and text not in rid:
                hide = True
            if abs(d) < thr:
                hide = True
            self.flux_tbl.setRowHidden(row, hide)

    # ---------------- Objective search ----------------
    def _setup_objective_completer(self):
        texts = [self.objective_combo.itemText(i) for i in range(self.objective_combo.count())]
        comp = QCompleter(texts, self)
        comp.setCaseSensitivity(Qt.CaseInsensitive)
        comp.setFilterMode(Qt.MatchContains)
        comp.setCompletionMode(QCompleter.PopupCompletion)
        self.objective_combo.setCompleter(comp)

    def _on_objective_text_edited(self, text: str):
        t = (text or "").strip().lower()
        if not t:
            return
        for i in range(self.objective_combo.count()):
            if t in self.objective_combo.itemText(i).lower():
                self.objective_combo.blockSignals(True)
                self.objective_combo.setCurrentIndex(i)
                self.objective_combo.lineEdit().setText(text)
                self.objective_combo.blockSignals(False)
                return

    def populate_objective_combo(self):
        if self.base_model is None:
            self.objective_combo.clear()
            self.objective_combo.setEnabled(False)
            return

        self.objective_combo.blockSignals(True)
        self.objective_combo.clear()
        self.objective_combo.addItem("Keep model objective", None)

        for rxn in sorted(self.base_model.reactions, key=lambda r: r.id):
            desc = get_description(rxn)
            label = f"{rxn.id} — {desc}" if desc else rxn.id
            self.objective_combo.addItem(label, rxn.id)

        self.objective_combo.setEnabled(True)
        self.objective_combo.blockSignals(False)
        self._setup_objective_completer()

    def _apply_selected_objective_to_model(self, model: cobra.Model):
        rid = self.objective_combo.currentData()
        if not rid:
            return
        model.objective = model.reactions.get_by_id(str(rid))
        model.objective_direction = "max"

    def on_objective_changed(self, *_):
        rid = self.objective_combo.currentData()
        if not rid:
            self.statusBar().showMessage("Objective: Keep model objective")
        else:
            self.statusBar().showMessage(f"Objective set to: {rid} (applies on next Run FBA)")

    # ---------------- Recent scenarios ----------------
    def _load_recent_scenarios(self):
        self.recent_scenarios = []
        try:
            if RECENT_FILE.exists():
                data = json.loads(RECENT_FILE.read_text(encoding="utf-8"))
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
        try:
            RECENT_FILE.write_text(json.dumps(self.recent_scenarios[:RECENT_MAX], indent=2), encoding="utf-8")
        except Exception:
            pass

    def _add_recent_scenario(self, file_path: str):
        file_path = str(Path(file_path))
        if file_path in self.recent_scenarios:
            self.recent_scenarios.remove(file_path)
        self.recent_scenarios.insert(0, file_path)
        self.recent_scenarios = self.recent_scenarios[:RECENT_MAX]
        self._save_recent_scenarios()
        self._rebuild_recent_menu()

    # ---------------- Busy state ----------------
    def set_busy(self, busy: bool, message: str = ""):
        self.is_running = busy
        for w in (
            self.open_btn, self.run_btn, self.close_uptakes_btn,
            self.reset_medium_btn, self.reset_reactions_btn,
            self.analysis_type, self.topn_spin, self.objective_combo,
            self.medium_search, self.medium_table,
            self.rxn_global_search, self.met_filter, self.met_filter_mode, self.reaction_table,
            self.clear_rxn_overrides_btn, self.remove_selected_overrides_btn,
            self.gene_global_search, self.genes_tab_list, self.gene_rxn_table,
            self.gene_search, self.gene_list, self.add_ko_btn, self.clear_ko_btn,
            self.rxn_search, self.rxn_list, self.oe_factor, self.oe_scale_lb,
            self.add_oe_btn, self.remove_oe_btn, self.clear_oe_btn,
            self.temp_ub_spin, self.set_temp_ub_btn, self.remove_temp_ub_btn, self.clear_temp_ub_btn,
            self.oe_active_list, self.temp_active_list,
            # New analysis tabs controls
            self.fva_flux_search, self.pfba_flux_search,
            self.deletion_search,
            self.robustness_rxn, self.robustness_min, self.robustness_max, self.robustness_steps, self.robustness_bound_ub, self.robustness_bound_lb,
            self.envelope_product, self.envelope_steps,
            self.sampling_size, self.sampling_search,
            # Map tab controls
            self.map_search, self.map_depth, self.map_threshold, self.map_view_combo,
            self.compare_mode,
            self.editor_tabs,
            self.tools_check_btn, self.tools_install_btn, self.memote_btn, self.memote_cancel_btn,
        ):
            enabled = not busy
            if w is self.run_btn:
                enabled = enabled and (self.base_model is not None)
            if w is self.analysis_type:
                enabled = enabled and (self.base_model is not None)
            if w is self.memote_btn:
                enabled = enabled and (self.base_model is not None) and (self._memote_proc is None)
            if w is self.memote_cancel_btn:
                enabled = enabled and (self._memote_proc is not None)
            w.setEnabled(enabled)
        if message:
            self.status_lbl.setText(message)
            self.statusBar().showMessage(message)

    # ---------------- Menu ----------------
    def _build_menu(self):
        menubar = self.menuBar()
        self.file_menu = menubar.addMenu("&File")
        run_menu = menubar.addMenu("&Run")
        medium_menu = menubar.addMenu("&Medium")
        model_menu = menubar.addMenu("&Model")
        help_menu = menubar.addMenu("&Help")

        act_open = QAction("Open SBML...", self)
        act_open.setShortcut(QKeySequence.Open)
        act_open.triggered.connect(self.open_sbml)
        self.file_menu.addAction(act_open)

        # (1) Save As (SBML)
        act_save_as = QAction("Save As (SBML)...", self)
        act_save_as.setShortcut(QKeySequence.SaveAs)
        act_save_as.triggered.connect(self.save_current_model_as)
        self.file_menu.addAction(act_save_as)

        self.file_menu.addSeparator()

        act_export = QAction("Export scenario...", self)
        act_export.triggered.connect(self.export_scenario)
        self.file_menu.addAction(act_export)

        act_import = QAction("Import scenario...", self)
        act_import.triggered.connect(self.import_scenario)
        self.file_menu.addAction(act_import)

        self.recent_menu = self.file_menu.addMenu("Recent scenarios")
        self._rebuild_recent_menu()

        self.file_menu.addSeparator()

        act_export_results = QAction("Export results (CSV)...", self)
        act_export_results.triggered.connect(self.export_results_csv)
        self.file_menu.addAction(act_export_results)

        act_export_excel = QAction("Export results (Excel)...", self)
        act_export_excel.triggered.connect(self.export_results_excel)
        self.file_menu.addAction(act_export_excel)

        act_export_json = QAction("Export results (JSON)...", self)
        act_export_json.triggered.connect(self.export_results_json)
        self.file_menu.addAction(act_export_json)

        act_export_all = QAction("Export All (scenario+results+chart)...", self)
        act_export_all.triggered.connect(self.export_all)
        self.file_menu.addAction(act_export_all)

        self.file_menu.addSeparator()

        act_exit = QAction("Exit", self)
        act_exit.setShortcut(QKeySequence.Quit)
        act_exit.triggered.connect(self.close)
        self.file_menu.addAction(act_exit)

        act_run = QAction("Run Analysis", self)
        act_run.setShortcut(QKeySequence("Ctrl+R"))
        act_run.triggered.connect(self.run_fba)
        run_menu.addAction(act_run)

        run_menu.addSeparator()

        act_undo = QAction("Undo", self)
        act_undo.setShortcut(QKeySequence.Undo)
        act_undo.triggered.connect(self.undo)
        run_menu.addAction(act_undo)

        act_redo = QAction("Redo", self)
        act_redo.setShortcut(QKeySequence.Redo)
        act_redo.triggered.connect(self.redo)
        run_menu.addAction(act_redo)

        run_menu.addSeparator()

        act_sensitivity = QAction("Sensitivity Analysis...", self)
        act_sensitivity.setShortcut(QKeySequence("Ctrl+Shift+S"))
        act_sensitivity.triggered.connect(self.run_sensitivity_analysis)
        run_menu.addAction(act_sensitivity)

        act_batch = QAction("Batch Analysis...", self)
        act_batch.setShortcut(QKeySequence("Ctrl+B"))
        act_batch.triggered.connect(self.run_batch_analysis)
        run_menu.addAction(act_batch)

        run_menu.addSeparator()

        act_search = QAction("Search & Filter...", self)
        act_search.setShortcut(QKeySequence.Find)
        act_search.triggered.connect(self.show_search_dialog)
        run_menu.addAction(act_search)

        # Explicit shortcuts to ensure they work on all platforms
        QShortcut(QKeySequence("Ctrl+Q"), self, activated=self.close)
        QShortcut(QKeySequence("Ctrl+Shift+S"), self, activated=self.run_sensitivity_analysis)

        act_close_uptakes = QAction("Close uptakes (LB<0 → 0)", self)
        act_close_uptakes.triggered.connect(self.close_all_uptakes)
        medium_menu.addAction(act_close_uptakes)

        act_reset_medium = QAction("Reset medium", self)
        act_reset_medium.triggered.connect(self.reset_medium)
        medium_menu.addAction(act_reset_medium)

        medium_menu.addSeparator()
        act_export_medium_csv = QAction("Export Medium bounds (CSV)...", self)
        act_export_medium_csv.triggered.connect(self.export_medium_bounds_csv)
        medium_menu.addAction(act_export_medium_csv)
        act_import_medium_csv = QAction("Import Medium bounds (CSV)...", self)
        act_import_medium_csv.triggered.connect(self.import_medium_bounds_csv)
        medium_menu.addAction(act_import_medium_csv)

        # Model menu - validation and checks
        act_validate_model = QAction("Validate model...", self)
        act_validate_model.triggered.connect(self.validate_model)
        model_menu.addAction(act_validate_model)

        act_check_mass_balance = QAction("Check mass balance...", self)
        act_check_mass_balance.triggered.connect(self.check_mass_balance)
        model_menu.addAction(act_check_mass_balance)

        act_check_gpr = QAction("Validate GPR syntax...", self)
        act_check_gpr.triggered.connect(self.check_gpr_syntax)
        model_menu.addAction(act_check_gpr)

        act_shortcuts = QAction("Keyboard shortcuts", self)
        act_shortcuts.triggered.connect(self.show_shortcuts)
        help_menu.addAction(act_shortcuts)

        act_about = QAction("About MetaboDesk", self)
        act_about.triggered.connect(self.show_about_new)
        help_menu.addAction(act_about)

    def _rebuild_recent_menu(self):
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

    def show_shortcuts(self):
        QMessageBox.information(
            self,
            "Keyboard Shortcuts",
            "Ctrl+O            Open SBML\n"
            "Ctrl+R            Run Analysis\n"
            "Ctrl+Shift+S      Save As (SBML)\n"
        )

    def show_about(self):
        html = (
            "<h3>MetaboDesk</h3>"
            "<p>SBML → COBRApy tabanlı metabolik ağ analizi için hepsi bir arada masaüstü uygulaması.</p>"
            "<h4>Analiz Türleri</h4>"
            "<ul>"
            "<li><b>FBA</b>: Akı dengesi analizi. ‘Objective’ ile hedef reaksiyonu seçin, Run FBA ile çalıştırın.</li>"
            "<li><b>pFBA</b>: Parsimonious FBA. Toplam akıyı minimize eder; pFBA sekmesinde sonuçları görüntüler.</li>"
            "<li><b>FVA</b>: Akı değişkenliği analizi. Reaksiyonlar için min–max aralıklarını hesaplar ve çizer.</li>"
            "<li><b>SGD</b>: Tekil gen silme. Her gen için büyüme kaybını raporlar.</li>"
            "<li><b>SRD</b>: Tekil reaksiyon silme. Her reaksiyon için büyüme kaybını raporlar.</li>"
            "<li><b>Robustness</b>: Belirli bir reaksiyon alt/üst sınırı süpürülür; hedef değişimi çizilir.</li>"
            "<li><b>Production Envelope</b>: Ürün üretimi ile büyüme arasındaki değiş tokuşu görselleştirir.</li>"
            "<li><b>Flux Sampling</b>: ACHR örnekleme ile akış dağılımlarını istatistiksel olarak çıkarır.</li>"
            "</ul>"
            "<h4>İpuçları</h4>"
            "<ul>"
            "<li>‘Medium’ ve ‘Reactions’ tablolarında değerleri düzenleyip Run FBA ile etkisini görebilirsiniz.</li>"
            "<li>‘Objective’ kutusuna yazmaya başlayın; otomatik tamamlamadan arayın.</li>"
            "<li>‘View → Theme’ ile Light/Dark tema geçişi yapabilirsiniz.</li>"
            "<li>‘Export All’ menüsü senaryoyu, sonuçları ve grafikleri birlikte kaydeder.</li>"
            "</ul>"
            "<p><i>© Emir Ay, 2026</i></p>"
        )
        QMessageBox.information(self, "About MetaboDesk", html)

    # NEW: Replace show_about with better content
    def show_about_new(self):
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
            QMessageBox.critical(self, "Export failed", str(e))

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
            QMessageBox.critical(self, "Import failed", str(e))

    def populate_medium_table(self):
        self.medium_table.setRowCount(0)
        for rxn in self.exchange_rxns:
            row = self.medium_table.rowCount()
            self.medium_table.insertRow(row)

            item_id = QTableWidgetItem(rxn.id)
            item_id.setFlags(item_id.flags() & ~Qt.ItemIsEditable)

            item_name = QTableWidgetItem(rxn.name or "")
            item_name.setFlags(item_name.flags() & ~Qt.ItemIsEditable)

            self.medium_table.setItem(row, 0, item_id)
            self.medium_table.setItem(row, 1, item_name)
            self.medium_table.setItem(row, 2, QTableWidgetItem(str(float(rxn.lower_bound))))
            self.medium_table.setItem(row, 3, QTableWidgetItem(str(float(rxn.upper_bound))))

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

        self.gene_rxn_table.setRowCount(0)

        for rid in rxn_ids:
            try:
                rxn = self.base_model.reactions.get_by_id(rid)
            except Exception:
                continue

            row = self.gene_rxn_table.rowCount()
            self.gene_rxn_table.insertRow(row)

            item_id = QTableWidgetItem(rxn.id)
            item_desc = QTableWidgetItem(get_description(rxn))
            eq = ""
            try:
                eq = rxn.reaction
            except Exception:
                eq = ""
            item_eq = QTableWidgetItem(eq)
            item_gpr = QTableWidgetItem(get_gpr(rxn))

            for it in (item_id, item_desc, item_eq, item_gpr):
                it.setFlags(it.flags() & ~Qt.ItemIsEditable)

            self.gene_rxn_table.setItem(row, 0, item_id)
            self.gene_rxn_table.setItem(row, 1, item_desc)
            self.gene_rxn_table.setItem(row, 2, item_eq)
            self.gene_rxn_table.setItem(row, 3, item_gpr)

        self.gene_rxn_details.setPlainText("")

    # ---------------- Reactions table ----------------
    def _parse_metabolite_filter_terms(self) -> list[str]:
        raw = self.met_filter.text().strip()
        if not raw:
            return []
        parts = [p.strip().lower() for p in raw.replace(";", ",").split(",")]
        return [p for p in parts if p]

    def _reaction_row_to_id(self, row: int) -> str | None:
        it = self.reaction_table.item(row, 0)
        return it.text().strip() if it else None

    def _rxn_contains_term(self, rxn: cobra.Reaction, term: str) -> bool:
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
        self.reaction_table.setRowCount(0)

        for rxn in sorted(self.base_model.reactions, key=lambda r: r.id):
            row = self.reaction_table.rowCount()
            self.reaction_table.insertRow(row)

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
                it.setFlags(it.flags() & ~Qt.ItemIsEditable)

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

        self.reaction_bound_overrides[rxn_id] = (lb, ub)
        self.model_dirty = True
        self.statusBar().showMessage(f"Override set: {rxn_id} LB={lb:g}, UB={ub:g}")

    def clear_reaction_overrides(self):
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
            self.knockout_genes.add(gid)
            self.model_dirty = True
            self.update_knockout_label()

    def clear_knockouts(self):
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
        if str(rxn_id) in self.overexpression_reactions:
            self.model_dirty = True
        self.overexpression_reactions.pop(str(rxn_id), None)
        self.refresh_overexpression_lists()
        self.update_selected_rxn_info()

    def clear_overexpressions(self):
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
        if str(rxn_id) in self.temporary_upper_bound_overrides:
            self.model_dirty = True
        self.temporary_upper_bound_overrides.pop(str(rxn_id), None)
        self.refresh_overexpression_lists()
        self.update_selected_rxn_info()

    def clear_temp_upper_bounds(self):
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
    def compute_fva(self, model: cobra.Model) -> dict:
        """
        Compute Flux Variability Analysis (FVA).
        Returns a dict with min/max flux for each reaction.
        """
        try:
            from cobra.flux_analysis import flux_variability_analysis
            fva_result = flux_variability_analysis(model)
            result = {}
            for rid in fva_result.index:
                result[rid] = {
                    "min": float(fva_result.loc[rid, "minimum"]),
                    "max": float(fva_result.loc[rid, "maximum"]),
                }
            return result
        except Exception as e:
            raise ValueError(f"FVA computation failed: {e}")

    def compute_pfba(self, model: cobra.Model) -> dict:
        """
        Compute Parsimonious FBA (pFBA).
        First optimizes the objective, then minimizes total flux.
        """
        try:
            from cobra.flux_analysis import pfba
            solution = pfba(model)
            result = {
                "status": str(solution.status),
                "objective": float(solution.objective_value),
                "flux": solution.fluxes.to_dict(),
            }
            return result
        except Exception as e:
            raise ValueError(f"pFBA computation failed: {e}")

    def compute_sgd(self, model: cobra.Model) -> dict:
        """Single Gene Deletion analysis."""
        try:
            from cobra.flux_analysis import single_gene_deletion
            result_df = single_gene_deletion(model)
            result = {}
            for gene_id in result_df.index:
                result[gene_id] = float(result_df.loc[gene_id, "growth"])
            return result
        except Exception as e:
            raise ValueError(f"SGD computation failed: {e}")

    def compute_srd(self, model: cobra.Model) -> dict:
        """Single Reaction Deletion analysis."""
        try:
            from cobra.flux_analysis import single_reaction_deletion
            result_df = single_reaction_deletion(model)
            result = {}
            for rxn_id in result_df.index:
                result[rxn_id] = float(result_df.loc[rxn_id, "growth"])
            return result
        except Exception as e:
            raise ValueError(f"SRD computation failed: {e}")


    def compute_robustness(self, model: cobra.Model, rxn_id: str, min_val: float, max_val: float, steps: int, bound_type: str = "ub") -> dict:
        """Robustness analysis: sweep a reaction's bound and measure objective."""
        result = {"values": [], "objectives": []}
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            original_lb = float(rxn.lower_bound)
            original_ub = float(rxn.upper_bound)
            
            for val in list(dict.fromkeys([min_val + i * (max_val - min_val) / (steps - 1) for i in range(steps)])):
                model_copy = model.copy()
                rxn_copy = model_copy.reactions.get_by_id(rxn_id)
                if bound_type.lower() == "ub":
                    rxn_copy.upper_bound = val
                else:  # lb
                    rxn_copy.lower_bound = val
                
                sol = model_copy.optimize()
                if str(sol.status) == "optimal":
                    result["values"].append(float(val))
                    result["objectives"].append(float(sol.objective_value))
            
            return result
        except Exception as e:
            raise ValueError(f"Robustness computation failed: {e}")

    def compute_production_envelope(self, model: cobra.Model, product_id: str, steps: int) -> dict:
        """Production envelope: measure trade-off between growth and product formation."""
        result = {"growth": [], "product": []}
        try:
            # Find objective (growth)
            growth_rxn_ids = [r.id for r in model.objective.expression.as_coefficients_dict().keys()]
            if not growth_rxn_ids:
                raise ValueError("No objective reaction found")
            
            product_rxn = model.reactions.get_by_id(product_id)
            
            for bound in list(dict.fromkeys([i * 10.0 / (steps - 1) for i in range(steps)])):
                model_copy = model.copy()
                prod_rxn = model_copy.reactions.get_by_id(product_id)
                prod_rxn.lower_bound = 0
                prod_rxn.upper_bound = bound
                
                sol = model_copy.optimize()
                if str(sol.status) == "optimal":
                    result["product"].append(float(bound))
                    result["growth"].append(float(sol.objective_value))
            
            return result
        except Exception as e:
            raise ValueError(f"Production envelope computation failed: {e}")

    def compute_flux_sampling(self, model: cobra.Model, sample_size: int = 500) -> dict:
        """Flux sampling using ACHR method."""
        try:
            from cobra.sampling import ACHRSampler
            sampler = ACHRSampler(model)
            samples = sampler.sample(sample_size)
            
            result = {}
            for rxn_id in samples.columns:
                flux_vals = samples[rxn_id].values
                result[rxn_id] = {
                    "mean": float(flux_vals.mean()),
                    "stdev": float(flux_vals.std()),
                    "min": float(flux_vals.min()),
                    "max": float(flux_vals.max()),
                }
            return result, samples
        except Exception as e:
            raise ValueError(f"Flux sampling failed: {e}")

    # ---------------- Plotting & Results ----------------
    def _clear_chart(self):
        ax = self.canvas.ax
        ax.clear()
        ax.set_title("Run FBA to see flux comparison")
        ax.set_xlabel("Reaction")
        ax.set_ylabel("Flux")
        self.canvas.draw()

    def _plot_flux_comparison(self, rxn_ids, base_flux, compared_flux, base_label: str, compared_label: str, title: str):
        ax = self.canvas.ax
        ax.clear()
        if not rxn_ids:
            ax.set_title("No flux changes above threshold")
            self.canvas.draw()
            return
        bvals = [base_flux.get(r, 0.0) for r in rxn_ids]
        cvals = [compared_flux.get(r, 0.0) for r in rxn_ids]
        x = list(range(len(rxn_ids)))
        width = 0.4
        ax.bar([i - width / 2 for i in x], bvals, width=width, label=base_label)
        ax.bar([i + width / 2 for i in x], cvals, width=width, label=compared_label)
        ax.set_title(title)
        ax.set_xlabel("Reaction")
        ax.set_ylabel("Flux")
        ax.set_xticks(x)
        ax.set_xticklabels(rxn_ids, rotation=90, fontsize=8)
        ax.legend()
        self.canvas.draw()

    def _refresh_results_view_from_last_run(self, *_):
        if not self.last_run:
            return
        self.last_run["compare_mode"] = self.compare_mode.currentText()
        self._recompute_flux_rows_for_compare_mode()
        self._render_results_from_last_run()

    def _set_flux_table_headers_for_mode(self, mode: str):
        if mode == "Original vs Baseline":
            self.flux_tbl.setHorizontalHeaderLabels(["Reaction", "Original flux", "Baseline flux", "Delta (baseline - original)"])
        else:
            self.flux_tbl.setHorizontalHeaderLabels(["Reaction", "Baseline flux", "Compared flux", "Delta (compared - baseline)"])

    def _recompute_flux_rows_for_compare_mode(self):
        assert self.last_run is not None
        mode = self.last_run.get("compare_mode", "Gene knockout only")

        if mode == "Original vs Baseline":
            base_flux = self.last_run["original"]["flux"]
            cmp_flux = self.last_run["baseline"]["flux"]
        elif mode == "Gene knockout only":
            base_flux = self.last_run["baseline"]["flux"]
            cmp_flux = self.last_run["gene_knockout_only"]["flux"]
        else:
            base_flux = self.last_run["baseline"]["flux"]
            cmp_flux = self.last_run["overexpression_only"]["flux"]

        rows = []
        for rxn_id, b in base_flux.items():
            c = cmp_flux.get(rxn_id, 0.0)
            d = c - b
            if abs(d) > 1e-9:
                rows.append((rxn_id, b, c, d))

        rows.sort(key=lambda t: abs(t[3]), reverse=True)
        rows = rows[: int(self.topn_spin.value())]

        self.last_run["flux_rows"] = [
            {"reaction": rid, "baseline": float(b), "compared": float(c), "delta": float(d)}
            for rid, b, c, d in rows
        ]

    def _render_results_from_last_run(self):
        assert self.last_run is not None
        mode = self.last_run.get("compare_mode", "Gene knockout only")

        original = self.last_run.get("original") or self.last_run["baseline"]
        baseline = self.last_run["baseline"]
        gene_knockout_only = self.last_run["gene_knockout_only"]
        overexpression_only = self.last_run["overexpression_only"]

        # Find biomass flux values
        biomass_info = ""
        if self.base_model:
            for rxn in self.base_model.reactions:
                if 'biomass' in rxn.id.lower() or 'biomass' in rxn.name.lower():
                    biomass_id = rxn.id
                    orig_biomass = original['flux'].get(biomass_id, 0.0)
                    base_biomass = baseline['flux'].get(biomass_id, 0.0)
                    ko_biomass = gene_knockout_only['flux'].get(biomass_id, 0.0)
                    oe_biomass = overexpression_only['flux'].get(biomass_id, 0.0)
                    biomass_info = f"\nBiomass ({biomass_id}):\n  Original: {orig_biomass:.6g} | Baseline: {base_biomass:.6g} | Knockout: {ko_biomass:.6g} | Overexpression: {oe_biomass:.6g}"
                    break

        self.results_lbl.setText(
            f"Original (no changes): {original['status']} | obj={original['objective']}\n"
            f"Baseline: {baseline['status']} | obj={baseline['objective']}\n"
            f"Gene knockout only: {gene_knockout_only['status']} | obj={gene_knockout_only['objective']}\n"
            f"Overexpression only: {overexpression_only['status']} | obj={overexpression_only['objective']}"
            f"{biomass_info}\n"
            f"\nCompare mode: {mode} | Top N={int(self.topn_spin.value())}\n"
            f"Reaction overrides: {len(self.reaction_bound_overrides)}"
        )

        self._set_flux_table_headers_for_mode(mode)

        rows = self.last_run.get("flux_rows", [])
        self.flux_tbl.setRowCount(0)
        rxn_ids_for_plot = []
        for row in rows:
            rid = row["reaction"]
            b = row["baseline"]
            c = row["compared"]
            d = row["delta"]
            r = self.flux_tbl.rowCount()
            self.flux_tbl.insertRow(r)
            self.flux_tbl.setItem(r, 0, QTableWidgetItem(rid))
            self.flux_tbl.setItem(r, 1, QTableWidgetItem(f"{b:.6g}"))
            self.flux_tbl.setItem(r, 2, QTableWidgetItem(f"{c:.6g}"))
            self.flux_tbl.setItem(r, 3, QTableWidgetItem(f"{d:.6g}"))
            rxn_ids_for_plot.append(rid)

        # Apply filters after rendering
        self.filter_flux_table()

        if mode == "Original vs Baseline":
            base_flux = original["flux"]
            cmp_flux = baseline["flux"]
            base_label = "Original"
            cmp_label = "Baseline"
            title = "Flux comparison (Original vs Baseline)"
        elif mode == "Gene knockout only":
            base_flux = baseline["flux"]
            cmp_flux = gene_knockout_only["flux"]
            base_label = "Baseline"
            cmp_label = "Gene knockout only"
            title = "Flux comparison (Baseline vs Gene knockout only)"
        else:
            base_flux = baseline["flux"]
            cmp_flux = overexpression_only["flux"]
            base_label = "Baseline"
            cmp_label = "Overexpression only"
            title = "Flux comparison (Baseline vs Overexpression only)"

        self._plot_flux_comparison(rxn_ids_for_plot, base_flux, cmp_flux, base_label, cmp_label, title)
        self.tabs.setCurrentWidget(self.tab_results)

    # ---------------- Run FBA ----------------
    def run_fba(self):
        if self.base_model is None or self.is_running:
            return

        analysis_type = self.analysis_type.currentText()
        
        if "FVA" in analysis_type:
            self.run_fva()
        elif "pFBA" in analysis_type:
            self.run_pfba()
        elif "Single Gene Deletion" in analysis_type:
            self.run_sgd()
        elif "Single Reaction Deletion" in analysis_type:
            self.run_srd()
        elif "Robustness" in analysis_type:
            self.run_robustness()
        elif "Production Envelope" in analysis_type:
            self.run_production_envelope()
        elif "Flux Sampling" in analysis_type:
            self.run_flux_sampling()
        else:
            self.run_standard_fba()

    def run_standard_fba(self):
        """Run standard FBA with baseline, knockout, and overexpression comparisons."""
        self.set_busy(True, "Running FBA (original + baseline + gene knockout only + overexpression only)...")
        QApplication.processEvents()

        try:
            original_model = (self.original_model_snapshot.copy() if self.original_model_snapshot is not None else self.base_model.copy())
            self._apply_selected_objective_to_model(original_model)
            original_sol = original_model.optimize()

            baseline_model = self.base_model.copy()
            self.apply_medium_table_bounds_to_model(baseline_model)
            self.apply_reaction_overrides_to_model(baseline_model)
            self._apply_selected_objective_to_model(baseline_model)
            baseline_sol = baseline_model.optimize()

            knockout_model = self.base_model.copy()
            self.apply_medium_table_bounds_to_model(knockout_model)
            self.apply_reaction_overrides_to_model(knockout_model)
            self.apply_knockouts_to_model(knockout_model)
            self._apply_selected_objective_to_model(knockout_model)
            knockout_sol = knockout_model.optimize()

            overexpression_model = self.base_model.copy()
            self.apply_medium_table_bounds_to_model(overexpression_model)
            self.apply_reaction_overrides_to_model(overexpression_model)
            self.apply_temp_upper_bound_overrides(overexpression_model)
            self.apply_overexpression_to_model(overexpression_model)
            self._apply_selected_objective_to_model(overexpression_model)
            overexpression_sol = overexpression_model.optimize()

            self.last_run = {
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "analysis_type": "FBA",
                "original": {"status": str(original_sol.status), "objective": float(original_sol.objective_value), "flux": original_sol.fluxes.to_dict()},
                "baseline": {"status": str(baseline_sol.status), "objective": float(baseline_sol.objective_value), "flux": baseline_sol.fluxes.to_dict()},
                "gene_knockout_only": {"status": str(knockout_sol.status), "objective": float(knockout_sol.objective_value), "flux": knockout_sol.fluxes.to_dict()},
                "overexpression_only": {"status": str(overexpression_sol.status), "objective": float(overexpression_sol.objective_value), "flux": overexpression_sol.fluxes.to_dict()},
                "compare_mode": self.compare_mode.currentText(),
                "flux_rows": [],
            }
            self._recompute_flux_rows_for_compare_mode()

        except Exception as e:
            QMessageBox.critical(self, "Run failed", str(e))
            self.set_busy(False, "Run failed.")
            return

        self._render_results_from_last_run()
        self.set_busy(False, "Ready.")

    def run_fva(self):
        """Run Flux Variability Analysis (FVA)."""
        self.set_busy(True, "Running FVA (Flux Variability Analysis)...")
        QApplication.processEvents()

        try:
            baseline_model = self.base_model.copy()
            self.apply_medium_table_bounds_to_model(baseline_model)
            self.apply_reaction_overrides_to_model(baseline_model)
            self._apply_selected_objective_to_model(baseline_model)
            
            # First optimize to ensure feasibility
            baseline_sol = baseline_model.optimize()
            if str(baseline_sol.status) != "optimal":
                raise ValueError(f"Model optimization failed: {baseline_sol.status}")
            
            # Compute FVA
            fva_result = self.compute_fva(baseline_model)

            self.last_run = {
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "analysis_type": "FVA",
                "baseline": {"status": str(baseline_sol.status), "objective": float(baseline_sol.objective_value), "flux": baseline_sol.fluxes.to_dict()},
                "fva_result": fva_result,
                "flux_rows": [],
            }
            
            # Prepare flux rows for display (sorted by flux range)
            rows = []
            for rxn_id, minmax in fva_result.items():
                min_flux = minmax["min"]
                max_flux = minmax["max"]
                flux_range = abs(max_flux - min_flux)
                if flux_range > 1e-9:  # Only show reactions with non-zero range
                    rows.append({
                        "reaction": rxn_id,
                        "min_flux": float(min_flux),
                        "max_flux": float(max_flux),
                        "range": float(flux_range),
                    })
            
            rows.sort(key=lambda x: abs(x["range"]), reverse=True)
            rows = rows[: int(self.topn_spin.value())]
            self.last_run["flux_rows"] = rows

        except Exception as e:
            QMessageBox.critical(self, "FVA failed", str(e))
            self.set_busy(False, "FVA failed.")
            return

        self._render_fva_results()
        self.set_busy(False, "Ready.")

    def run_pfba(self):
        """Run Parsimonious FBA (pFBA)."""
        self.set_busy(True, "Running pFBA (Parsimonious FBA)...")
        QApplication.processEvents()

        try:
            baseline_model = self.base_model.copy()
            self.apply_medium_table_bounds_to_model(baseline_model)
            self.apply_reaction_overrides_to_model(baseline_model)
            self._apply_selected_objective_to_model(baseline_model)
            
            # Compute pFBA
            pfba_sol = self.compute_pfba(baseline_model)

            # For comparison, also compute standard FBA
            baseline_model_std = self.base_model.copy()
            self.apply_medium_table_bounds_to_model(baseline_model_std)
            self.apply_reaction_overrides_to_model(baseline_model_std)
            self._apply_selected_objective_to_model(baseline_model_std)
            baseline_sol = baseline_model_std.optimize()

            self.last_run = {
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "analysis_type": "pFBA",
                "baseline": {"status": str(baseline_sol.status), "objective": float(baseline_sol.objective_value), "flux": baseline_sol.fluxes.to_dict()},
                "pfba": pfba_sol,
                "compare_mode": "pFBA vs FBA",
                "flux_rows": [],
            }
            
            # Prepare flux rows for display (flux differences)
            rows = []
            for rxn_id, baseline_flux in baseline_sol.fluxes.items():
                pfba_flux = pfba_sol["flux"].get(rxn_id, 0.0)
                delta = abs(baseline_flux) - abs(pfba_flux)
                if abs(delta) > 1e-9:
                    rows.append({
                        "reaction": rxn_id,
                        "fba_flux": float(baseline_flux),
                        "pfba_flux": float(pfba_flux),
                        "delta": float(delta),
                    })
            
            rows.sort(key=lambda x: abs(x["delta"]), reverse=True)
            rows = rows[: int(self.topn_spin.value())]
            self.last_run["flux_rows"] = rows

        except Exception as e:
            QMessageBox.critical(self, "pFBA failed", str(e))
            self.set_busy(False, "pFBA failed.")
            return

        self._render_pfba_results()
        self.set_busy(False, "Ready.")

    def run_sgd(self):
        """Run Single Gene Deletion analysis."""
        self.set_busy(True, "Running Single Gene Deletion (SGD)...")
        QApplication.processEvents()

        try:
            baseline_model = self.base_model.copy()
            self.apply_medium_table_bounds_to_model(baseline_model)
            self.apply_reaction_overrides_to_model(baseline_model)
            self._apply_selected_objective_to_model(baseline_model)
            
            # Compute baseline growth
            baseline_sol = baseline_model.optimize()
            wt_growth = float(baseline_sol.objective_value)
            
            # Run SGD
            sgd_result = self.compute_sgd(baseline_model)
            
            self.last_run = {
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "analysis_type": "SGD",
                "wt_growth": wt_growth,
                "sgd_result": sgd_result,
            }
            
        except Exception as e:
            QMessageBox.critical(self, "SGD failed", str(e))
            self.set_busy(False, "SGD failed.")
            return

        self._render_sgd_results()
        self.set_busy(False, "Ready.")

    def run_srd(self):
        """Run Single Reaction Deletion analysis."""
        self.set_busy(True, "Running Single Reaction Deletion (SRD)...")
        QApplication.processEvents()

        try:
            baseline_model = self.base_model.copy()
            self.apply_medium_table_bounds_to_model(baseline_model)
            self.apply_reaction_overrides_to_model(baseline_model)
            self._apply_selected_objective_to_model(baseline_model)
            
            # Compute baseline growth
            baseline_sol = baseline_model.optimize()
            wt_growth = float(baseline_sol.objective_value)
            
            # Run SRD
            srd_result = self.compute_srd(baseline_model)
            
            self.last_run = {
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "analysis_type": "SRD",
                "wt_growth": wt_growth,
                "srd_result": srd_result,
            }
            
        except Exception as e:
            QMessageBox.critical(self, "SRD failed", str(e))
            self.set_busy(False, "SRD failed.")
            return

        self._render_srd_results()
        self.set_busy(False, "Ready.")

    def run_robustness(self):
        """Run Robustness Analysis."""
        rxn_id = self.robustness_rxn.text().strip()
        if not rxn_id:
            QMessageBox.warning(self, "Missing", "Enter a reaction ID for robustness analysis.")
            return
        
        sweep_ub = self.robustness_bound_ub.isChecked()
        sweep_lb = self.robustness_bound_lb.isChecked()
        if not sweep_ub and not sweep_lb:
            QMessageBox.warning(self, "Missing", "Select at least one bound type (UB or LB).")
            return
        
        bound_type = "ub" if sweep_ub else "lb"
            
        self.set_busy(True, f"Running Robustness Analysis for {rxn_id} ({bound_type.upper()})...")
        QApplication.processEvents()

        try:
            baseline_model = self.base_model.copy()
            self.apply_medium_table_bounds_to_model(baseline_model)
            self.apply_reaction_overrides_to_model(baseline_model)
            self._apply_selected_objective_to_model(baseline_model)
            
            # Run robustness
            min_val = float(self.robustness_min.value())
            max_val = float(self.robustness_max.value())
            steps = int(self.robustness_steps.value())
            
            robustness_result = self.compute_robustness(baseline_model, rxn_id, min_val, max_val, steps, bound_type)
            
            self.last_run = {
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "analysis_type": "Robustness",
                "rxn_id": rxn_id,
                "bound_type": bound_type,
                "robustness_result": robustness_result,
            }
            
        except Exception as e:
            QMessageBox.critical(self, "Robustness failed", str(e))
            self.set_busy(False, "Robustness failed.")
            return

        self._render_robustness_results()
        self.set_busy(False, "Ready.")

    def run_production_envelope(self):
        """Run Production Envelope Analysis."""
        product_id = self.envelope_product.text().strip()
        if not product_id:
            QMessageBox.warning(self, "Missing", "Enter a product reaction ID for production envelope.")
            return
            
        self.set_busy(True, f"Running Production Envelope for {product_id}...")
        QApplication.processEvents()

        try:
            baseline_model = self.base_model.copy()
            self.apply_medium_table_bounds_to_model(baseline_model)
            self.apply_reaction_overrides_to_model(baseline_model)
            self._apply_selected_objective_to_model(baseline_model)
            
            # Run production envelope
            steps = int(self.envelope_steps.value())
            
            envelope_result = self.compute_production_envelope(baseline_model, product_id, steps)
            
            self.last_run = {
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "analysis_type": "Production Envelope",
                "product_id": product_id,
                "envelope_result": envelope_result,
            }
            
        except Exception as e:
            QMessageBox.critical(self, "Production Envelope failed", str(e))
            self.set_busy(False, "Production Envelope failed.")
            return

        self._render_envelope_results()
        self.set_busy(False, "Ready.")

    def run_flux_sampling(self):
        """Run Flux Sampling (ACHR)."""
        self.set_busy(True, "Running Flux Sampling...")
        QApplication.processEvents()

        try:
            baseline_model = self.base_model.copy()
            self.apply_medium_table_bounds_to_model(baseline_model)
            self.apply_reaction_overrides_to_model(baseline_model)
            self._apply_selected_objective_to_model(baseline_model)
            
            # Run sampling
            sample_size = int(self.sampling_size.value())
            
            sampling_result, samples_df = self.compute_flux_sampling(baseline_model, sample_size)
            
            self.last_run = {
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "analysis_type": "Flux Sampling",
                "sampling_result": sampling_result,
                "samples_df": samples_df,
            }
            
        except Exception as e:
            QMessageBox.critical(self, "Flux Sampling failed", str(e))
            self.set_busy(False, "Flux Sampling failed.")
            return

        self._render_sampling_results()
        self.set_busy(False, "Ready.")

    def filter_fva_table(self):
        """Filter FVA table by reaction ID."""
        text = (self.fva_flux_search.text() or "").strip().lower()
        for row in range(self.fva_tbl.rowCount()):
            rid = (self.fva_tbl.item(row, 0).text() if self.fva_tbl.item(row, 0) else "").lower()
            hide = bool(text) and text not in rid
            self.fva_tbl.setRowHidden(row, hide)

    def _render_fva_results(self):
        """Render FVA results in the FVA tab."""
        assert self.last_run is not None
        baseline = self.last_run.get("baseline", {})
        fva_result = self.last_run.get("fva_result", {})
        
        self.fva_info_lbl.setText(
            f"Baseline: {baseline.get('status', 'N/A')} | Objective: {baseline.get('objective', 'N/A'):.6g}\n"
            f"Top N reactions by flux range: {int(self.topn_spin.value())}"
        )
        
        # Populate flux table
        rows = self.last_run.get("flux_rows", [])
        self.fva_tbl.setRowCount(0)
        rxn_ids_for_plot = []
        
        for row in rows:
            rid = row["reaction"]
            min_f = row["min_flux"]
            max_f = row["max_flux"]
            flux_range = row["range"]
            
            r = self.fva_tbl.rowCount()
            self.fva_tbl.insertRow(r)
            self.fva_tbl.setItem(r, 0, QTableWidgetItem(rid))
            self.fva_tbl.setItem(r, 1, QTableWidgetItem(f"{min_f:.6g}"))
            self.fva_tbl.setItem(r, 2, QTableWidgetItem(f"{max_f:.6g}"))
            self.fva_tbl.setItem(r, 3, QTableWidgetItem(f"{flux_range:.6g}"))
            rxn_ids_for_plot.append(rid)
        
        # Apply filters
        self.filter_fva_table()
        
        # Plot FVA ranges
        ax = self.fva_canvas.ax
        ax.clear()
        if not rxn_ids_for_plot:
            ax.set_title("No reactions with flux variability found")
            self.fva_canvas.draw()
        else:
            # Plot min/max ranges for each reaction
            x = list(range(len(rxn_ids_for_plot)))
            mins = [self.last_run["fva_result"].get(rid, {}).get("min", 0.0) for rid in rxn_ids_for_plot]
            maxs = [self.last_run["fva_result"].get(rid, {}).get("max", 0.0) for rid in rxn_ids_for_plot]
            
            ax.barh(x, [m - n for m, n in zip(maxs, mins)], left=mins, label="Flux range", color="steelblue", alpha=0.7)
            ax.axvline(x=0, color="black", linestyle="-", linewidth=0.5)
            ax.set_title("FVA: Flux Variability Ranges")
            ax.set_ylabel("Reaction")
            ax.set_xlabel("Flux")
            ax.set_yticks(x)
            ax.set_yticklabels(rxn_ids_for_plot, fontsize=8)
            ax.legend()
            self.fva_canvas.draw()
        
        self.tabs.setCurrentWidget(self.tab_fva)

    def export_fva_csv(self):
        if not self.last_run or not self.last_run.get("fva_result"):
            QMessageBox.warning(self, "FVA", "Run FVA first.")
            return
        file_path, _ = QFileDialog.getSaveFileName(self, "Export FVA (CSV)", str(Path.home() / "fva_results.csv"), "CSV files (*.csv);")
        if not file_path:
            return
        try:
            with open(file_path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["reaction", "min", "max", "range"]) 
                for rid, mm in self.last_run["fva_result"].items():
                    rng = float(mm["max"]) - float(mm["min"]) 
                    w.writerow([rid, mm["min"], mm["max"], rng])
            self.statusBar().showMessage(f"FVA exported: {file_path}")
        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e))

    def _render_pfba_results(self):
        """Render pFBA results in the Results tab."""
        assert self.last_run is not None
        baseline = self.last_run.get("baseline", {})
        pfba = self.last_run.get("pfba", {})
        
        self.results_lbl.setText(
            f"pFBA (Parsimonious FBA) Results\n"
            f"Standard FBA: {baseline.get('status', 'N/A')} | obj={baseline.get('objective', 'N/A')}\n"
            f"pFBA: {pfba.get('status', 'N/A')} | obj={pfba.get('objective', 'N/A')}\n"
            f"Top N reactions by flux difference: {int(self.topn_spin.value())}"
        )
        
        # Update table headers for pFBA
        self.flux_tbl.setHorizontalHeaderLabels(["Reaction", "FBA flux", "pFBA flux", "Difference"])
        
        # Populate flux table
        rows = self.last_run.get("flux_rows", [])
        self.flux_tbl.setRowCount(0)
        rxn_ids_for_plot = []
        
        for row in rows:
            rid = row["reaction"]
            fba_f = row["fba_flux"]
            pfba_f = row["pfba_flux"]
            diff = row["delta"]
            
            r = self.flux_tbl.rowCount()
            self.flux_tbl.insertRow(r)
            self.flux_tbl.setItem(r, 0, QTableWidgetItem(rid))
            self.flux_tbl.setItem(r, 1, QTableWidgetItem(f"{fba_f:.6g}"))
            self.flux_tbl.setItem(r, 2, QTableWidgetItem(f"{pfba_f:.6g}"))
            self.flux_tbl.setItem(r, 3, QTableWidgetItem(f"{diff:.6g}"))
            rxn_ids_for_plot.append(rid)
        
        # Apply filters
        self.filter_flux_table()
        
        # Plot pFBA comparison
        ax = self.canvas.ax
        ax.clear()
        if not rxn_ids_for_plot:
            ax.set_title("FBA and pFBA results are identical")
            self.canvas.draw()
        else:
            baseline_flux = self.last_run["baseline"]["flux"]
            pfba_flux = pfba["flux"]
            
            fba_vals = [abs(baseline_flux.get(r, 0.0)) for r in rxn_ids_for_plot]
            pfba_vals = [abs(pfba_flux.get(r, 0.0)) for r in rxn_ids_for_plot]
            
            x = list(range(len(rxn_ids_for_plot)))
            width = 0.4
            ax.bar([i - width / 2 for i in x], fba_vals, width=width, label="FBA", alpha=0.7)
            ax.bar([i + width / 2 for i in x], pfba_vals, width=width, label="pFBA", alpha=0.7)
            
            ax.set_title("pFBA vs Standard FBA: Flux Comparison")
            ax.set_xlabel("Reaction")
            ax.set_ylabel("Absolute Flux")
            ax.set_xticks(x)
            ax.set_xticklabels(rxn_ids_for_plot, rotation=90, fontsize=8)
            ax.legend()
            self.pfba_canvas.draw()
        
        self.tabs.setCurrentWidget(self.tab_pfba)

    def export_pfba_csv(self):
        if not self.last_run or not self.last_run.get("pfba"):
            QMessageBox.warning(self, "pFBA", "Run pFBA first.")
            return
        file_path, _ = QFileDialog.getSaveFileName(self, "Export pFBA (CSV)", str(Path.home() / "pfba_results.csv"), "CSV files (*.csv);")
        if not file_path:
            return
        try:
            rows = self.last_run.get("flux_rows", [])
            with open(file_path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["reaction", "fba_flux", "pfba_flux", "difference"]) 
                for r in rows:
                    w.writerow([r["reaction"], r["fba_flux"], r["pfba_flux"], r["delta"]])
            self.statusBar().showMessage(f"pFBA exported: {file_path}")
        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e))

    def filter_deletion_table(self):
        """Filter deletion analysis table."""
        text = (self.deletion_search.text() or "").strip().lower()
        for row in range(self.deletion_tbl.rowCount()):
            rid = (self.deletion_tbl.item(row, 0).text() if self.deletion_tbl.item(row, 0) else "").lower()
            hide = bool(text) and text not in rid
            self.deletion_tbl.setRowHidden(row, hide)

    def _render_sgd_results(self):
        """Render SGD results."""
        assert self.last_run is not None
        wt_growth = self.last_run.get("wt_growth", 0.0)
        sgd_result = self.last_run.get("sgd_result", {})
        
        self.deletion_info_lbl.setText(
            f"WT Growth: {wt_growth:.4f}\n"
            f"Showing all {len(sgd_result)} genes"
        )
        
        # Populate table
        self.deletion_tbl.setHorizontalHeaderLabels(["Gene ID", "WT Growth", "Growth on Deletion", "Δ Growth (%)", "Category"])
        self.deletion_tbl.setRowCount(0)
        genes_for_plot = []
        growth_vals = []
        
        sorted_genes = sorted(sgd_result.items(), key=lambda x: x[1])
        
        for gene_id, growth in sorted_genes[:int(self.topn_spin.value())]:
            r = self.deletion_tbl.rowCount()
            self.deletion_tbl.insertRow(r)
            self.deletion_tbl.setItem(r, 0, QTableWidgetItem(gene_id))
            self.deletion_tbl.setItem(r, 1, QTableWidgetItem(f"{wt_growth:.4f}"))
            self.deletion_tbl.setItem(r, 2, QTableWidgetItem(f"{growth:.4f}"))
            delta_pct = (0.0 if wt_growth == 0 else ((growth - wt_growth) / wt_growth) * 100.0)
            self.deletion_tbl.setItem(r, 3, QTableWidgetItem(f"{delta_pct:.2f}"))
            if growth < 1e-6:
                cat = "Lethal"
            elif delta_pct < -75:
                cat = "Severe"
            elif delta_pct < -25:
                cat = "Moderate"
            else:
                cat = "Mild"
            self.deletion_tbl.setItem(r, 4, QTableWidgetItem(cat))
            genes_for_plot.append(gene_id)
            growth_vals.append(growth)
        
        self.filter_deletion_table()
        
        # Plot with value labels
        ax = self.deletion_canvas.ax
        ax.clear()
        if genes_for_plot:
            bars = ax.barh(genes_for_plot, growth_vals, color="coral", alpha=0.7, edgecolor="darkred", linewidth=0.8)
            for bar, val in zip(bars, growth_vals):
                ax.text(val, bar.get_y() + bar.get_height()/2, f"{val:.3g}", fontsize=7, va="center", ha="left", style="italic")
            ax.set_xlabel("Growth Rate")
            ax.set_title("Single Gene Deletion - Growth on Deletion")
        self.deletion_canvas.draw()
        
        self.tabs.setCurrentWidget(self.tab_deletion)

    def export_sgd_csv(self):
        if not self.last_run or not self.last_run.get("sgd_result"):
            QMessageBox.warning(self, "SGD", "Run SGD first.")
            return
        file_path, _ = QFileDialog.getSaveFileName(self, "Export SGD (CSV)", str(Path.home() / "sgd_results.csv"), "CSV files (*.csv);")
        if not file_path:
            return
        try:
            wt = float(self.last_run.get("wt_growth", 0.0))
            with open(file_path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["gene_id", "wt_growth", "growth_on_deletion", "delta_pct"]) 
                for gid, growth in sorted(self.last_run["sgd_result"].items(), key=lambda x: x[1]):
                    delta_pct = (0.0 if wt == 0 else ((growth - wt) / wt) * 100.0)
                    w.writerow([gid, wt, growth, delta_pct])
            self.statusBar().showMessage(f"SGD exported: {file_path}")
        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e))

    def _render_srd_results(self):
        """Render SRD results."""
        assert self.last_run is not None
        wt_growth = self.last_run.get("wt_growth", 0.0)
        srd_result = self.last_run.get("srd_result", {})
        
        self.deletion_info_lbl.setText(
            f"WT Growth: {wt_growth:.4f}\n"
            f"Showing all {len(srd_result)} reactions"
        )
        
        # Populate table
        self.deletion_tbl.setHorizontalHeaderLabels(["Reaction ID", "WT Growth", "Growth on Deletion", "Δ Growth (%)", "Category"])
        self.deletion_tbl.setRowCount(0)
        rxns_for_plot = []
        growth_vals = []
        
        sorted_rxns = sorted(srd_result.items(), key=lambda x: x[1])
        
        for rxn_id, growth in sorted_rxns[:int(self.topn_spin.value())]:
            r = self.deletion_tbl.rowCount()
            self.deletion_tbl.insertRow(r)
            self.deletion_tbl.setItem(r, 0, QTableWidgetItem(rxn_id))
            self.deletion_tbl.setItem(r, 1, QTableWidgetItem(f"{wt_growth:.4f}"))
            self.deletion_tbl.setItem(r, 2, QTableWidgetItem(f"{growth:.4f}"))
            delta_pct = (0.0 if wt_growth == 0 else ((growth - wt_growth) / wt_growth) * 100.0)
            self.deletion_tbl.setItem(r, 3, QTableWidgetItem(f"{delta_pct:.2f}"))
            if growth < 1e-6:
                cat = "Lethal"
            elif delta_pct < -75:
                cat = "Severe"
            elif delta_pct < -25:
                cat = "Moderate"
            else:
                cat = "Mild"
            self.deletion_tbl.setItem(r, 4, QTableWidgetItem(cat))
            rxns_for_plot.append(rxn_id)
            growth_vals.append(growth)
        
        self.filter_deletion_table()
        
        # Plot with value labels
        ax = self.deletion_canvas.ax
        ax.clear()
        if rxns_for_plot:
            bars = ax.barh(rxns_for_plot, growth_vals, color="lightgreen", alpha=0.7, edgecolor="darkgreen", linewidth=0.8)
            for bar, val in zip(bars, growth_vals):
                ax.text(val, bar.get_y() + bar.get_height()/2, f"{val:.3g}", fontsize=7, va="center", ha="left", style="italic")
            ax.set_xlabel("Growth Rate")
            ax.set_title("Single Reaction Deletion - Growth on Deletion")
        self.deletion_canvas.draw()
        
        self.tabs.setCurrentWidget(self.tab_deletion)

    def export_srd_csv(self):
        if not self.last_run or not self.last_run.get("srd_result"):
            QMessageBox.warning(self, "SRD", "Run SRD first.")
            return
        file_path, _ = QFileDialog.getSaveFileName(self, "Export SRD (CSV)", str(Path.home() / "srd_results.csv"), "CSV files (*.csv);")
        if not file_path:
            return
        try:
            wt = float(self.last_run.get("wt_growth", 0.0))
            with open(file_path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["reaction_id", "wt_growth", "growth_on_deletion", "delta_pct"]) 
                for rid, growth in sorted(self.last_run["srd_result"].items(), key=lambda x: x[1]):
                    delta_pct = (0.0 if wt == 0 else ((growth - wt) / wt) * 100.0)
                    w.writerow([rid, wt, growth, delta_pct])
            self.statusBar().showMessage(f"SRD exported: {file_path}")
        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e))

    def _render_robustness_results(self):
        """Render Robustness analysis results."""
        assert self.last_run is not None
        rxn_id = self.last_run.get("rxn_id", "")
        robustness_result = self.last_run.get("robustness_result", {})
        
        values = robustness_result.get("values", [])
        objectives = robustness_result.get("objectives", [])
        
        self.robustness_info_lbl.setText(
            f"Robustness Analysis for {rxn_id}\n"
            f"Points computed: {len(values)}"
        )
        
        # Plot with value markers and labels
        ax = self.robustness_canvas.ax
        ax.clear()
        if values and objectives:
            ax.plot(values, objectives, marker="o", linestyle="-", color="steelblue", linewidth=2, markersize=6, markerfacecolor="lightblue", markeredgewidth=1.5)
            for x, y in zip(values, objectives):
                ax.text(x, y, f"{y:.3g}", fontsize=7, ha="center", va="bottom")
            ax.set_xlabel(f"{rxn_id} bound")
            ax.set_ylabel("Objective value")
            ax.set_title(f"Robustness: {rxn_id}")
            ax.grid(True, alpha=0.3)
        self.robustness_canvas.draw()
        
        self.tabs.setCurrentWidget(self.tab_robustness)

    def _render_envelope_results(self):
        """Render Production Envelope results."""
        assert self.last_run is not None
        product_id = self.last_run.get("product_id", "")
        envelope_result = self.last_run.get("envelope_result", {})
        
        product_vals = envelope_result.get("product", [])
        growth_vals = envelope_result.get("growth", [])
        
        self.envelope_info_lbl.setText(
            f"Production Envelope: {product_id} vs Growth\n"
            f"Points computed: {len(product_vals)}"
        )
        
        # Plot
        ax = self.envelope_canvas.ax
        ax.clear()
        if product_vals and growth_vals:
            ax.scatter(product_vals, growth_vals, s=50, alpha=0.6, color="purple")
            ax.fill_between(product_vals, growth_vals, alpha=0.2, color="purple")
            ax.set_xlabel(f"{product_id} production rate")
            ax.set_ylabel("Growth rate")
            ax.set_title(f"Production Envelope: {product_id}")
            ax.grid(True, alpha=0.3)
        self.envelope_canvas.draw()
        
        self.tabs.setCurrentWidget(self.tab_envelope)

    def _render_sampling_results(self):
        """Render Flux Sampling results."""
        assert self.last_run is not None
        sampling_result = self.last_run.get("sampling_result", {})
        samples_df = self.last_run.get("samples_df", None)
        
        self.sampling_info_lbl.setText(
            f"Flux Sampling (ACHR)\n"
            f"Sample size: {len(samples_df) if samples_df is not None else 0}\n"
            f"Reactions analyzed: {len(sampling_result)}"
        )
        
        # Populate table
        self.sampling_tbl.setRowCount(0)
        rxns_for_plot = []
        means = []
        
        sorted_rxns = sorted(sampling_result.items(), key=lambda x: abs(x[1]["mean"]), reverse=True)
        
        for rxn_id, stats in sorted_rxns[:int(self.topn_spin.value())]:
            r = self.sampling_tbl.rowCount()
            self.sampling_tbl.insertRow(r)
            self.sampling_tbl.setItem(r, 0, QTableWidgetItem(rxn_id))
            self.sampling_tbl.setItem(r, 1, QTableWidgetItem(f"{stats['mean']:.6g}"))
            self.sampling_tbl.setItem(r, 2, QTableWidgetItem(f"{stats['stdev']:.6g}"))
            self.sampling_tbl.setItem(r, 3, QTableWidgetItem(f"[{stats['min']:.3g}, {stats['max']:.3g}]"))
            rxns_for_plot.append(rxn_id)
            means.append(stats["mean"])
        
        # Plot
        ax = self.sampling_canvas.ax
        ax.clear()
        if rxns_for_plot:
            colors = ["green" if m > 0 else "red" for m in means]
            ax.barh(rxns_for_plot, means, color=colors, alpha=0.6)
            ax.set_xlabel("Mean flux")
            ax.set_title("Flux Sampling - Mean Flux Distribution")
        self.sampling_canvas.draw()
        
        self.tabs.setCurrentWidget(self.tab_sampling)

    def export_sampling_csv(self):
        if not self.last_run or not self.last_run.get("sampling_result"):
            QMessageBox.warning(self, "Sampling", "Run sampling first.")
            return
        file_path, _ = QFileDialog.getSaveFileName(self, "Export Sampling Stats (CSV)", str(Path.home() / "sampling_stats.csv"), "CSV files (*.csv);")
        if not file_path:
            return
        try:
            sampling_result = self.last_run["sampling_result"]
            with open(file_path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["reaction_id", "mean_flux", "stdev", "min_flux", "max_flux"])
                for rxn_id, stats in sorted(sampling_result.items(), key=lambda x: abs(x[1]["mean"]), reverse=True):
                    w.writerow([rxn_id, stats["mean"], stats["stdev"], stats["min"], stats["max"]])
            self.statusBar().showMessage(f"Sampling exported: {file_path}")
        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e))

    def filter_pfba_table(self):
        """Filter pFBA table by reaction ID."""
        text = (self.pfba_flux_search.text() or "").strip().lower()
        for row in range(self.pfba_tbl.rowCount()):
            rid = (self.pfba_tbl.item(row, 0).text() if self.pfba_tbl.item(row, 0) else "").lower()
            hide = bool(text) and text not in rid
            self.pfba_tbl.setRowHidden(row, hide)
        self.tabs.setCurrentWidget(self.tab_results)

    def filter_sampling_table(self):
        """Filter sampling table."""
        text = (self.sampling_search.text() or "").strip().lower()
        for row in range(self.sampling_tbl.rowCount()):
            rid = (self.sampling_tbl.item(row, 0).text() if self.sampling_tbl.item(row, 0) else "").lower()
            hide = bool(text) and text not in rid
            self.sampling_tbl.setRowHidden(row, hide)

    def _parse_float(self, val, allow_comma: bool = True) -> float:
        """Robust float parser for user inputs. Supports comma decimal.
        Raises ValueError if parsing fails.
        """
        if isinstance(val, (int, float)):
            return float(val)
        if val is None:
            raise ValueError("Missing numeric value")
        s = str(val).strip()
        if allow_comma:
            s = s.replace(",", ".")
        return float(s)
    # ---------------- Scenario export/import ----------------
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
            QMessageBox.critical(self, "Export failed", str(e))
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
            QMessageBox.critical(self, "Import failed", str(e))
            return

        self._add_recent_scenario(file_path)
        self.statusBar().showMessage(f"Scenario imported: {file_path}")

    def _capture_exchange_bounds_from_table(self) -> dict[str, dict[str, float]]:
        bounds: dict[str, dict[str, float]] = {}
        for i in range(self.medium_table.rowCount()):
            rxn_id = self.medium_table.item(i, 0).text().strip()
            lb = self._parse_float(self.medium_table.item(i, 2).text())
            ub = self._parse_float(self.medium_table.item(i, 3).text())
            bounds[rxn_id] = {"lb": lb, "ub": ub}
        return bounds

    def _apply_exchange_bounds_to_table(self, bounds: dict):
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
        with open(file_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["# MetaboDesk results export"])
            w.writerow(["# timestamp", self.last_run.get("timestamp", "")])
            w.writerow(["# flux_filter_search", (self.flux_search.text() if hasattr(self, "flux_search") else "")])
            w.writerow(["# flux_filter_abs_delta_threshold", (float(self.delta_thresh.value()) if hasattr(self, "delta_thresh") else 0.0)])
            w.writerow([])
            w.writerow(["original_status", self.last_run["original"]["status"]])
            w.writerow(["original_objective", self.last_run["original"]["objective"]])
            w.writerow(["baseline_status", self.last_run["baseline"]["status"]])
            w.writerow(["baseline_objective", self.last_run["baseline"]["objective"]])
            w.writerow(["gene_knockout_only_status", self.last_run["gene_knockout_only"]["status"]])
            w.writerow(["gene_knockout_only_objective", self.last_run["gene_knockout_only"]["objective"]])
            w.writerow(["overexpression_only_status", self.last_run["overexpression_only"]["status"]])
            w.writerow(["overexpression_only_objective", self.last_run["overexpression_only"]["objective"]])
            w.writerow([])
            w.writerow(["reaction", "baseline_flux", "compared_flux", "delta", "compare_mode"])
            for row in self.last_run.get("flux_rows", []):
                w.writerow([row["reaction"], row["baseline"], row["compared"], row["delta"], self.last_run.get("compare_mode", "")])

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
            QMessageBox.critical(self, "Export failed", str(e))
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
            QMessageBox.critical(self, "Export All failed", str(e))
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

        self.set_busy(True, f"Loading: {file_path}")
        QApplication.processEvents()

        try:
            model = cobra.io.read_sbml_model(file_path)
        except Exception as e:
            QMessageBox.critical(self, "Failed to load SBML", str(e))
            self.set_busy(False, "Failed to load model.")
            return

        self.base_model = model
        self.original_model_snapshot = model.copy()
        self.last_run = None

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

        # enable memote
        self.memote_btn.setEnabled(True)

        self._clear_chart()
        self.flux_tbl.setRowCount(0)
        self.results_lbl.setText("Objective: -")

        self.set_busy(False, "Ready.")
        self.tabs.setCurrentWidget(self.tab_medium)

    def _populate_reaction_lists(self):
        """Populate reaction ID combos for Robustness and Envelope tabs."""
        try:
            ids = [r.id for r in self.base_model.reactions]
        except Exception:
            ids = []
        ids_sorted = sorted(ids)
        # store for filtering
        self._all_reaction_ids = ids_sorted
        # Robustness combo
        if hasattr(self, "robustness_combo"):
            self.robustness_combo.blockSignals(True)
            self.robustness_combo.clear()
            self.robustness_combo.addItems(ids_sorted)
            self.robustness_combo.blockSignals(False)
        # Envelope combo
        if hasattr(self, "envelope_combo"):
            self.envelope_combo.blockSignals(True)
            self.envelope_combo.clear()
            self.envelope_combo.addItems(ids_sorted)
            self.envelope_combo.blockSignals(False)

    def _filter_robustness_combo(self, text: str):
        s = (text or "").strip().lower()
        items = [x for x in getattr(self, "_all_reaction_ids", []) if s in x.lower()]
        self.robustness_combo.blockSignals(True)
        self.robustness_combo.clear()
        self.robustness_combo.addItems(items)
        self.robustness_combo.blockSignals(False)

    def _filter_envelope_combo(self, text: str):
        s = (text or "").strip().lower()
        items = [x for x in getattr(self, "_all_reaction_ids", []) if s in x.lower()]
        self.envelope_combo.blockSignals(True)
        self.envelope_combo.clear()
        self.envelope_combo.addItems(items)
        self.envelope_combo.blockSignals(False)

    # ---------------- Model editor tabs ----------------
    def _build_editor_metabolite_tab(self):
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

    def _build_editor_patch_tab(self):
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
        try:
            self.ed_patch_preview.setPlainText(json.dumps(self.model_patch, indent=2))
        except Exception:
            pass

    def _rebuild_editor_reaction_dropdowns(self):
        return

    def _rebuild_editor_gene_dropdowns(self):
        return

    def _rebuild_editor_metabolite_dropdowns(self):
        return

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
        self.populate_reaction_table()
        self.populate_medium_table()

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
        self.populate_gene_list()
        self.populate_genes_tab()

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
        except Exception:
            pass

        try:
            rxn.subtract_metabolites(rxn.metabolites)
        except Exception:
            pass

        mets = {self.base_model.metabolites.get_by_id(mid): float(coeff) for mid, coeff in stoich.items()}
        rxn.add_metabolites(mets)

        rxn.lower_bound = float(lb)
        rxn.upper_bound = float(ub)
        rxn.reversibility = bool(reversible)

        try:
            rxn.gene_reaction_rule = gpr
        except Exception:
            pass

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

        self.populate_reaction_table()
        self.populate_medium_table()
        self.populate_overexpression_reaction_list()
        self.populate_objective_combo()
        self.populate_genes_tab()

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
        self.populate_reaction_table()
        self.populate_medium_table()

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
            QMessageBox.critical(self, "Delete failed", str(e))
            return

        if rid not in self.model_patch["deleted_reactions"]:
            self.model_patch["deleted_reactions"].append(rid)
        self._refresh_patch_preview()
        self.model_dirty = True

        self.ed_disable_info.setPlainText(f"Deleted reaction: {rid}")
        self.populate_reaction_table()
        self.populate_medium_table()
        self.populate_overexpression_reaction_list()
        self.populate_objective_combo()
        self.populate_genes_tab()

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
            QMessageBox.critical(self, "Export failed", str(e))
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
            QMessageBox.critical(self, "Import failed", str(e))
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
            QMessageBox.critical(self, "Apply failed", str(e))
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
            QMessageBox.critical(self, "Export failed", str(e))
            return
        self.statusBar().showMessage(f"Patched SBML exported: {file_path}")

    def _apply_patch_to_model(self, patch: dict):
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
            except Exception:
                pass

            try:
                rxn.subtract_metabolites(rxn.metabolites)
            except Exception:
                pass
            mets = {self.base_model.metabolites.get_by_id(mid): float(coeff) for mid, coeff in stoich.items()}
            rxn.add_metabolites(mets)

            try:
                rxn.lower_bound = float(r.get("lb"))
                rxn.upper_bound = float(r.get("ub"))
            except Exception:
                pass

            try:
                rxn.gene_reaction_rule = str(r.get("gpr") or "")
            except Exception:
                pass

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
    def cancel_memote(self):
        if self._memote_proc is None:
            return
        try:
            self.tools_log.appendPlainText("\nCancel requested...\n")
            self._memote_proc.kill()
        except Exception:
            pass

    def _append_proc_output_to_log(self, proc: QProcess):
        data = proc.readAll().data().decode(errors="replace")
        if not data:
            return
        for line in data.splitlines():
            self.tools_log.appendPlainText(line)

    def run_memote_report_current_model(self):
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        sbml_path = self.current_sbml_path
        if self.model_dirty or sbml_path is None or not sbml_path.exists():
            msg = (
                "Your current model has unsaved changes (or no source file path).\n\n"
                "Please save the current model as a NEW SBML file before running memote.\n\n"
                "Continue with Save As?"
            )
            if QMessageBox.question(self, "Save required", msg, QMessageBox.Yes | QMessageBox.No) != QMessageBox.Yes:
                return
            saved = self.save_current_model_as()
            if saved is None:
                return
            sbml_path = saved

        default_html = Path.home() / f"{sbml_path.stem}_memote.html"
        report_path_str, _ = QFileDialog.getSaveFileName(
            self,
            "Save memote report as (HTML)",
            str(default_html),
            "HTML files (*.html);;All files (*.*)",
        )
        if not report_path_str:
            return
        report_path = Path(report_path_str)

        mm = self._micromamba_path()
        env_prefix = self._runtime_dir() / "envs" / "metabodesk"
        if not mm.exists() or not env_prefix.exists():
            QMessageBox.critical(self, "Tools missing", "Tools env not installed. Go to Tools → Install tools first.")
            return

        self.tabs.setCurrentWidget(self.tab_tools)
        self.tools_log.appendPlainText(f"Running memote on: {sbml_path}")
        self.tools_log.appendPlainText(f"Report target: {report_path}\n")

        args = [
            "run",
            "-r", str(self._runtime_dir()),
            "-p", str(env_prefix),
            "memote",
            "report",
            "snapshot",
            str(sbml_path),
            "--filename",
            str(report_path),
        ]

        proc = QProcess(self)
        self._memote_proc = proc
        self.set_busy(True, "Running memote...")
        proc.setProgram(str(mm))
        proc.setArguments(args)
        proc.setProcessChannelMode(QProcess.MergedChannels)

        proc.readyRead.connect(lambda: self._append_proc_output_to_log(proc))

        def on_finished(code, status):
            self._append_proc_output_to_log(proc)
            self.tools_log.appendPlainText(f"\nmemote finished with code={code}, status={status}\n")
            self._memote_proc = None
            self.set_busy(False, "Ready.")

            # If report exists, open even if memote failed (user asked earlier; useful for debug)
            if report_path.exists():
                QDesktopServices.openUrl(QUrl.fromLocalFile(str(report_path)))
                QMessageBox.information(self, "Memote finished", f"Report saved:\n{report_path}")
            else:
                QMessageBox.critical(self, "Memote failed", "memote failed and report was not created. See log in Tools tab.")

        proc.finished.connect(on_finished)
        proc.start()

    # ---------------- Tools installer helpers ----------------
    def _runtime_dir(self) -> Path:
        return Path(__file__).resolve().parent / "runtime"


    def validate_model(self):
        """Comprehensive model validation."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return
        
        issues = []
        try:
            # Check for reactions without genes or bounds
            for rxn in self.base_model.reactions:
                if not rxn.genes and rxn.lower_bound != 0 and rxn.upper_bound != 0:
                    issues.append(f"Reaction {rxn.id}: no genes assigned (orphaned)")
                if rxn.lower_bound > rxn.upper_bound:
                    issues.append(f"Reaction {rxn.id}: lower bound ({rxn.lower_bound}) > upper bound ({rxn.upper_bound})")
                if len(rxn.metabolites) == 0:
                    issues.append(f"Reaction {rxn.id}: no metabolites (empty)")
            
            # Check for orphan metabolites
            consumed = set()
            produced = set()
            for rxn in self.base_model.reactions:
                for met in rxn.metabolites:
                    if rxn.metabolites[met] > 0:
                        produced.add(met.id)
                    else:
                        consumed.add(met.id)
            
            orphan_produced = produced - consumed
            orphan_consumed = consumed - produced
            if orphan_produced:
                issues.append(f"{len(orphan_produced)} metabolites are only produced (no sink reactions)")
            if orphan_consumed:
                issues.append(f"{len(orphan_consumed)} metabolites are only consumed (no source reactions)")
            
            # Check objective
            try:
                obj_val = self.base_model.slim_optimize()
                if obj_val is None or (isinstance(obj_val, float) and obj_val < 1e-6):
                    issues.append("Objective value is zero or infeasible (check bounds and objective)")
            except:
                issues.append("Model optimization failed - check for infeasibility")
            
        except Exception as e:
            issues.append(f"Validation error: {e}")
        
        if not issues:
            QMessageBox.information(self, "Validation", "✓ Model validation passed - no major issues found.")
        else:
            msg = "Issues found:\n\n" + "\n".join(issues[:20])
            if len(issues) > 20:
                msg += f"\n... and {len(issues) - 20} more issues"
            QMessageBox.warning(self, "Model validation", msg)

    def check_mass_balance(self):
        """Check if reactions are mass balanced."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return
        
        try:
            imbalanced = []
            for rxn in self.base_model.reactions:
                if not rxn.metabolites:
                    continue
                # Simple check: sum of stoichiometric coefficients
                # For balanced reactions, should be close to zero (considering atom balance)
                coeff_sum = sum(rxn.metabolites.values())
                if abs(coeff_sum) > 0.01:
                    imbalanced.append(f"{rxn.id}: coeff_sum={coeff_sum:.4f}")
            
            if not imbalanced:
                QMessageBox.information(self, "Mass balance check", "✓ All reactions appear balanced (simple check).")
            else:
                msg = f"Found {len(imbalanced)} potentially imbalanced reactions:\n\n"
                msg += "\n".join(imbalanced[:15])
                if len(imbalanced) > 15:
                    msg += f"\n... and {len(imbalanced) - 15} more"
                QMessageBox.warning(self, "Mass balance check", msg)
        except Exception as e:
            QMessageBox.critical(self, "Check failed", str(e))

    def check_gpr_syntax(self):
        """Validate Gene-Protein-Reaction (GPR) syntax."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return
        
        import re
        issues = []
        try:
            for rxn in self.base_model.reactions:
                gpr = getattr(rxn, 'gene_reaction_rule', '') or ''
                if not gpr:
                    continue
                
                # Basic syntax checks
                open_parens = gpr.count('(')
                close_parens = gpr.count(')')
                if open_parens != close_parens:
                    issues.append(f"{rxn.id}: parentheses mismatch (open={open_parens}, close={close_parens})")
                
                # Check for invalid operators
                if 'OR' in gpr.upper() and 'or' in gpr:
                    issues.append(f"{rxn.id}: mixed case OR/or (should be uppercase)")
                
                # Validate tokens
                tokens = re.split(r'[()and&or| ]', gpr)
                for tok in tokens:
                    if tok and not re.match(r'^[A-Za-z0-9_.-]+$', tok):
                        if tok not in ('and', 'or', 'AND', 'OR'):
                            issues.append(f"{rxn.id}: invalid token '{tok}'")
                            break
            
            if not issues:
                QMessageBox.information(self, "GPR validation", "✓ All GPR expressions appear syntactically valid.")
            else:
                msg = f"Found {len(issues)} GPR syntax issues:\n\n"
                msg += "\n".join(issues[:20])
                if len(issues) > 20:
                    msg += f"\n... and {len(issues) - 20} more"
                QMessageBox.warning(self, "GPR validation", msg)
        except Exception as e:
            QMessageBox.critical(self, "Check failed", str(e))

    def _micromamba_path(self) -> Path:
        # Platform-aware name
        name = "micromamba.exe" if sys.platform.startswith("win") else "micromamba"
        return self._runtime_dir() / name

    def tools_check_status(self):
        mm = self._micromamba_path()
        env_prefix = self._runtime_dir() / "envs" / "metabodesk"
        if mm.exists() and env_prefix.exists():
            self.tools_status_lbl.setText(f"OK: micromamba env exists at {env_prefix}")
            return
        # Fallback: check pip-installed packages in current env
        try:
            import memote  # type: ignore
            import swiglpk  # type: ignore
            self.tools_status_lbl.setText("OK: pip tools available (memote, swiglpk)")
        except Exception:
            self.tools_status_lbl.setText("Tools not installed. Use Install tools (micromamba or pip).")

    def tools_install_tools(self):
        mm = self._micromamba_path()
        root_prefix = self._runtime_dir()
        env_prefix = self._runtime_dir() / "envs" / "metabodesk"

        self.tabs.setCurrentWidget(self.tab_tools)

        if mm.exists():
            msg = (
                "This will install tools via micromamba into:\n\n"
                f"{env_prefix}\n\n"
                "It can take several minutes and needs internet.\n\n"
                "Continue?"
            )
            if QMessageBox.question(self, "Confirm install", msg, QMessageBox.Yes | QMessageBox.No) != QMessageBox.Yes:
                return
            cmd = [
                str(mm), "create", "-y",
                "-r", str(root_prefix),
                "-p", str(env_prefix),
                "-c", "conda-forge",
                "python=3.11", "pip", "memote", "swiglpk",
            ]
            self.tools_log.appendPlainText("Running (micromamba):\n" + " ".join(cmd) + "\n")
            try:
                import subprocess
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)
                for line in p.stdout:
                    self.tools_log.appendPlainText(line.rstrip())
                    QApplication.processEvents()
                code = p.wait()
                if code != 0:
                    QMessageBox.critical(self, "Install failed", f"micromamba exited with code {code}")
                    self.tools_status_lbl.setText("Install failed. See log.")
                    return
            except Exception as e:
                QMessageBox.critical(self, "Install failed", str(e))
                self.tools_status_lbl.setText("Install failed. See error.")
                return
            QMessageBox.information(self, "Install complete", "Tools installed successfully (micromamba).")
            self.tools_check_status()
            return

        # Fallback: pip install into current environment
        msg = (
            "micromamba not found. Will install tools via pip into current environment.\n\n"
            "Packages: memote, swiglpk\n\nContinue?"
        )
        if QMessageBox.question(self, "Confirm pip install", msg, QMessageBox.Yes | QMessageBox.No) != QMessageBox.Yes:
            return
        try:
            import subprocess, sys as _sys
            cmd = [_sys.executable, "-m", "pip", "install", "memote", "swiglpk"]
            self.tools_log.appendPlainText("Running (pip):\n" + " ".join(cmd) + "\n")
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)
            for line in p.stdout:
                self.tools_log.appendPlainText(line.rstrip())
                QApplication.processEvents()
            code = p.wait()
            if code != 0:
                QMessageBox.critical(self, "Install failed", f"pip exited with code {code}")
                self.tools_status_lbl.setText("pip install failed. See log.")
                return
        except Exception as e:
            QMessageBox.critical(self, "Install failed", str(e))
            self.tools_status_lbl.setText("pip install failed. See error.")
            return
        QMessageBox.information(self, "Install complete", "Tools installed successfully (pip).")
        self.tools_check_status()

    # -------- Network Map --------
    def _render_network_map(self):
        """Trigger background map rendering - non-blocking."""
        Thread(target=self._render_network_map_bg, daemon=True).start()
    
    def _render_network_map_bg(self):
        """Background rendering - non-blocking."""
        try:
            if self.base_model is None:
                ax = self.map_canvas.ax
                ax.clear()
                ax.text(0.5, 0.5, "Load model first", ha="center", va="center", fontsize=12, transform=ax.transAxes)
                self.map_canvas.draw()
                return
            
            search_text = (self.map_search.text() or "").strip().lower()
            view_mode = self.map_view_combo.currentText()
            
            reactions = [r for r in self.base_model.reactions 
                        if not search_text or search_text in r.id.lower()]
            
            if not reactions:
                ax = self.map_canvas.ax
                ax.clear()
                ax.text(0.5, 0.5, "No reactions found", ha="center", va="center", fontsize=11, transform=ax.transAxes)
                self.map_canvas.draw()
                return
            
            reactions = reactions[:200]
            
            flux_data = {}
            fva_data = {}
            gene_ko_data = {}
            
            if self.last_run:
                try:
                    if "fba_flux" in self.last_run:
                        flux_data = self.last_run["fba_flux"].get("flux", {})
                except: pass
                try:
                    if "fva_table" in self.last_run:
                        for row in self.last_run["fva_table"]:
                            fva_data[row[0]] = {"min": float(row[1]), "max": float(row[2])}
                except: pass
                try:
                    if "deletion_table" in self.last_run and "Gene KO" in view_mode:
                        for row in self.last_run["deletion_table"]:
                            for rxn in reactions:
                                if any(g.id == row[0] for g in rxn.genes):
                                    gene_ko_data[rxn.id] = True
                except: pass
            
            G = nx.Graph()
            rxn_map = {}
            
            for rxn in reactions:
                flux = abs(float(flux_data.get(rxn.id, 0.0)))
                bounds = fva_data.get(rxn.id, {"min": -1000, "max": 1000})
                is_ko = rxn.id in gene_ko_data
                G.add_node(rxn.id, flux=flux, bounds=bounds, ko_affected=is_ko)
                rxn_map[rxn.id] = rxn
            
            for i, r1 in enumerate(reactions):
                mets1 = set(r1.metabolites.keys())
                for r2 in reactions[i+1:]:
                    if mets1 & set(r2.metabolites.keys()):
                        G.add_edge(r1.id, r2.id)
            
            if G.number_of_nodes() == 0:
                ax = self.map_canvas.ax
                ax.clear()
                ax.text(0.5, 0.5, "Empty network", ha="center", va="center", fontsize=11, transform=ax.transAxes)
                self.map_canvas.draw()
                return
            
            try:
                if G.number_of_nodes() <= 50:
                    pos = nx.spring_layout(G, k=0.5, iterations=10, seed=42)
                else:
                    pos = nx.circular_layout(G)
            except:
                pos = nx.circular_layout(G)
            
            node_colors = []
            node_sizes = []
            
            for node in G.nodes():
                data = G.nodes[node]
                flux = data.get("flux", 0.0)
                bounds = data.get("bounds", {})
                is_ko = data.get("ko_affected", False)
                
                if "Gene KO" in view_mode:
                    node_colors.append(2.0 if is_ko else 0.2)
                    node_sizes.append(300 if is_ko else 150)
                elif "FVA" in view_mode:
                    vmin = bounds.get("min", -1000)
                    vmax = bounds.get("max", 1000)
                    width = max(vmax - vmin, 0)
                    node_colors.append(min(width / 10, 3))
                    node_sizes.append(200)
                else:
                    node_colors.append(min(flux, 2))
                    node_sizes.append(150 + min(flux * 30, 100))
            
            ax = self.map_canvas.ax
            ax.clear()
            
            nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.15, width=0.5, edge_color="gray")
            
            max_color = max(node_colors) if node_colors else 1
            nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes,
                                  cmap="RdYlBu_r", vmin=0, vmax=max_color,
                                  ax=ax, alpha=0.85, edgecolors="black", linewidths=0.5)
            
            self._current_graph = G
            self._current_graph_pos = pos
            self._current_reactions_map = rxn_map
            
            if G.number_of_nodes() <= 30:
                labels = {n: n.split("_")[0][:5] for n in G.nodes()}
                nx.draw_networkx_labels(G, pos, labels, font_size=4, ax=ax)
            
            ax.set_title(f"Network ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges) - {view_mode}", 
                        fontsize=10, fontweight="bold")
            ax.axis("off")
            self.map_canvas.draw()
            
            self.map_info_lbl.setText(
                f"Reactions: {len(reactions)} | Edges: {G.number_of_edges()} | Click node for info | Toolbar: zoom/pan"
            )
        except Exception as e:
            print(f"Map error: {e}")
            ax = self.map_canvas.ax
            ax.clear()
            ax.text(0.5, 0.5, f"Error: {str(e)[:40]}", ha="center", va="center", fontsize=9, transform=ax.transAxes)
            self.map_canvas.draw()
    
    def _on_map_click(self, event):
        """Handle click on map nodes to show reaction details."""
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
            # Show details in map_details panel
            G = self._current_graph
            data = G.nodes[nearest]
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
            bounds = data.get("bounds", {})
            flux = data.get("flux", 0.0)
            is_ko = data.get("ko_affected", False)
            txt = (
                f"Reaction: {nearest}\nDescription: {desc}\nEquation: {eq}\n"
                f"Flux(abs): {float(flux):.6g}\nFVA: [{bounds.get('min','?')}, {bounds.get('max','?')}]\n"
                f"GPR: {gpr}\nKEGG: {kegg}\nEC: {ec}\nKO affected: {is_ko}"
            )
            self.map_details.setPlainText(txt)
        except Exception as e:
            self.map_details.setPlainText(f"Click error: {e}")

    def _apply_selected_solver(self, model: cobra.Model):
        """Apply selected solver to COBRApy model."""
        try:
            s = (self.solver_combo.currentText() or "Auto").lower()
            if s == "auto":
                return
            if s == "glpk":
                import optlang.glpk_interface as iface
            elif s == "highs":
                import optlang.scipy_interface as iface
            else:
                return
            model.solver = iface
        except Exception as e:
            # Graceful fallback
            self.statusBar().showMessage(f"Solver set failed: {e}")

    # ==================== NEW FEATURES ====================

    # -------- 1. EXPORT FUNCTIONS --------
    def export_results_csv(self):
        """Export FBA results to CSV"""
        if self.last_run is None:
            QMessageBox.warning(self, "Export Error", "No results to export. Run analysis first.")
            return
        
        path, _ = QFileDialog.getSaveFileName(self, "Export Results", "", "CSV Files (*.csv)")
        if not path:
            return
        
        try:
            with open(path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(["Reaction ID", "Flux", "Lower Bound", "Upper Bound"])
                
                if "fluxes" in self.last_run:
                    for rxn_id, flux in self.last_run["fluxes"].items():
                        if self.base_model and rxn_id in self.base_model.reactions:
                            rxn = self.base_model.reactions.get_by_id(rxn_id)
                            writer.writerow([rxn_id, flux, rxn.lower_bound, rxn.upper_bound])
            
            QMessageBox.information(self, "Export Success", f"Results exported to {path}")
        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"Failed to export: {str(e)}")

    def export_results_excel(self):
        """Export FBA results to Excel"""
        if not HAS_OPENPYXL:
            QMessageBox.warning(self, "Export Error", "openpyxl not installed. Install with: pip install openpyxl")
            return
        
        if self.last_run is None:
            QMessageBox.warning(self, "Export Error", "No results to export. Run analysis first.")
            return
        
        path, _ = QFileDialog.getSaveFileName(self, "Export Results", "", "Excel Files (*.xlsx)")
        if not path:
            return
        
        try:
            wb = openpyxl.Workbook()
            ws = wb.active
            ws.title = "Results"
            ws.append(["Reaction ID", "Flux", "Lower Bound", "Upper Bound"])
            
            if "fluxes" in self.last_run:
                for rxn_id, flux in self.last_run["fluxes"].items():
                    if self.base_model and rxn_id in self.base_model.reactions:
                        rxn = self.base_model.reactions.get_by_id(rxn_id)
                        ws.append([rxn_id, flux, rxn.lower_bound, rxn.upper_bound])
            
            wb.save(path)
            QMessageBox.information(self, "Export Success", f"Results exported to {path}")
        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"Failed to export: {str(e)}")

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
            QMessageBox.critical(self, "Export Error", f"Failed to export: {str(e)}")

    # -------- 2. SEARCH/FILTER FUNCTIONS --------
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
        """Run sensitivity analysis on objective"""
        if self.base_model is None:
            QMessageBox.warning(self, "Error", "Load a model first.")
            return
        
        dialog = QDialog(self)
        dialog.setWindowTitle("Sensitivity Analysis")
        dialog.resize(400, 250)
        
        layout = QFormLayout()
        
        min_spin = QDoubleSpinBox()
        min_spin.setMinimum(-100)
        min_spin.setMaximum(100)
        min_spin.setValue(-100)
        layout.addRow("Min value:", min_spin)
        
        max_spin = QDoubleSpinBox()
        max_spin.setMinimum(-100)
        max_spin.setMaximum(100)
        max_spin.setValue(100)
        layout.addRow("Max value:", max_spin)
        
        steps_spin = QSpinBox()
        steps_spin.setMinimum(5)
        steps_spin.setMaximum(100)
        steps_spin.setValue(20)
        layout.addRow("Steps:", steps_spin)
        
        buttons = QHBoxLayout()
        ok_btn = QPushButton("Run")
        cancel_btn = QPushButton("Cancel")
        
        def run_sens():
            dialog.close()
            self._perform_sensitivity_analysis(min_spin.value(), max_spin.value(), steps_spin.value())
        
        ok_btn.clicked.connect(run_sens)
        cancel_btn.clicked.connect(dialog.close)
        buttons.addStretch()
        buttons.addWidget(ok_btn)
        buttons.addWidget(cancel_btn)
        
        layout.addRow(buttons)
        dialog.setLayout(layout)
        dialog.exec()

    def _perform_sensitivity_analysis(self, min_val, max_val, steps):
        """Perform the sensitivity analysis"""
        if self.base_model is None or self.base_model.objective is None:
            return
        
        self.is_running = True
        self.run_btn.setEnabled(False)
        
        try:
            progress = QProgressDialog("Running sensitivity analysis...", "Cancel", 0, steps, self)
            progress.setWindowModality(Qt.ApplicationModal)
            
            objective_rxn = list(self.base_model.objective.keys())[0]
            values = []
            objectives = []
            
            for i, val in enumerate(range(steps)):
                if progress.wasCanceled():
                    break
                
                param_val = min_val + (max_val - min_val) * (i / (steps - 1))
                with self.base_model as model:
                    objective_rxn.lower_bound = param_val
                    try:
                        sol = model.optimize()
                        objectives.append(sol.objective_value if sol.status == "optimal" else None)
                    except:
                        objectives.append(None)
                
                values.append(param_val)
                progress.setValue(i + 1)
            
            self.sensitivity_results = {
                "parameter_values": values,
                "objective_values": objectives
            }
            
            QMessageBox.information(self, "Success", "Sensitivity analysis completed!")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Sensitivity analysis failed: {str(e)}")
        finally:
            self.is_running = False
            self.run_btn.setEnabled(True)

    # -------- 4. BATCH ANALYSIS --------
    def run_batch_analysis(self):
        """Run analysis on multiple SBML files"""
        files, _ = QFileDialog.getOpenFileNames(self, "Select SBML Files", "", "SBML Files (*.xml *.sbml)")
        if not files:
            return
        
        progress = QProgressDialog("Running batch analysis...", "Cancel", 0, len(files), self)
        progress.setWindowModality(Qt.ApplicationModal)
        
        results = []
        
        for i, file_path in enumerate(files):
            if progress.wasCanceled():
                break
            
            try:
                model = cobra.io.read_sbml_model(file_path)
                with model:
                    solution = model.optimize()
                    results.append({
                        "file": Path(file_path).name,
                        "status": solution.status if solution else "No solution",
                        "objective": solution.objective_value if solution else None
                    })
            except Exception as e:
                results.append({
                    "file": Path(file_path).name,
                    "status": "Error",
                    "objective": str(e)
                })
            
            progress.setValue(i + 1)
        
        progress.close()
        
        # Show results
        self._show_batch_results(results)

    def _show_batch_results(self, results):
        """Display batch analysis results"""
        dialog = QDialog(self)
        dialog.setWindowTitle("Batch Analysis Results")
        dialog.resize(600, 400)
        
        layout = QVBoxLayout()
        
        table = QTableWidget()
        table.setColumnCount(3)
        table.setHorizontalHeaderLabels(["File", "Status", "Objective Value"])
        table.horizontalHeader().setStretchLastSection(True)
        
        for result in results:
            row = table.rowCount()
            table.insertRow(row)
            table.setItem(row, 0, QTableWidgetItem(result["file"]))
            table.setItem(row, 1, QTableWidgetItem(result["status"]))
            table.setItem(row, 2, QTableWidgetItem(str(result["objective"])))
        
        layout.addWidget(table)
        
        buttons = QHBoxLayout()
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
        """Save current model state for undo"""
        if self.base_model is None:
            return
        
        state = {
            "knockout_genes": self.knockout_genes.copy(),
            "overexpression_reactions": self.overexpression_reactions.copy(),
            "reaction_bound_overrides": self.reaction_bound_overrides.copy(),
            "timestamp": datetime.now()
        }
        
        self.undo_stack.append(state)
        if len(self.undo_stack) > self.max_undo_steps:
            self.undo_stack.pop(0)
        
        self.redo_stack.clear()

    def undo(self):
        """Undo last change"""
        if not self.undo_stack:
            QMessageBox.information(self, "Undo", "Nothing to undo.")
            return
        
        current_state = {
            "knockout_genes": self.knockout_genes.copy(),
            "overexpression_reactions": self.overexpression_reactions.copy(),
            "reaction_bound_overrides": self.reaction_bound_overrides.copy()
        }
        self.redo_stack.append(current_state)
        
        state = self.undo_stack.pop()
        self.knockout_genes = state["knockout_genes"]
        self.overexpression_reactions = state["overexpression_reactions"]
        self.reaction_bound_overrides = state["reaction_bound_overrides"]
        self.model_dirty = True
        self.statusBar().showMessage("Undo completed.")

    def redo(self):
        """Redo last undone change"""
        if not self.redo_stack:
            QMessageBox.information(self, "Redo", "Nothing to redo.")
            return
        
        current_state = {
            "knockout_genes": self.knockout_genes.copy(),
            "overexpression_reactions": self.overexpression_reactions.copy(),
            "reaction_bound_overrides": self.reaction_bound_overrides.copy()
        }
        self.undo_stack.append(current_state)
        
        state = self.redo_stack.pop()
        self.knockout_genes = state["knockout_genes"]
        self.overexpression_reactions = state["overexpression_reactions"]
        self.reaction_bound_overrides = state["reaction_bound_overrides"]
        self.model_dirty = True
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

def main():
    app = QApplication(sys.argv)
    win = MainWindow()
    
    try:
        import pyi_splash
        pyi_splash.close()
    except Exception:
        pass

    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()