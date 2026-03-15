"""MainWindow class for MetaboDesk — assembles all mixin functionality."""

import os
import logging

import cobra
from pathlib import Path

from PySide6.QtCore import Qt, QUrl, QProcess, QTimer
from PySide6.QtGui import QAction, QKeySequence, QDesktopServices, QShortcut, QIcon
from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QMessageBox, QTableWidget, QTableWidgetItem, QAbstractItemView,
    QHeaderView, QListWidget, QLineEdit, QSpinBox, QTabWidget,
    QDoubleSpinBox, QCheckBox, QComboBox, QPlainTextEdit, QCompleter,
    QProgressBar,
)

from metabodesk_core.utils import _resource_path, get_description
from metabodesk_core.widgets import TextPopup, MplCanvas, AnalysisWorker
from metabodesk_core.config import AppConfig

# Defer heavy NavigationToolbar import — loaded on first use via _get_nav_toolbar
_NavigationToolbar = None

def _get_nav_toolbar():
    global _NavigationToolbar
    if _NavigationToolbar is None:
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
        _NavigationToolbar = NavigationToolbar2QT
    return _NavigationToolbar

logger = logging.getLogger("MetaboDesk")


from metabodesk_core.mixin_io import IOMixin
from metabodesk_core.mixin_medium import MediumMixin
from metabodesk_core.mixin_reactions import ReactionsMixin
from metabodesk_core.mixin_analysis import AnalysisMixin
from metabodesk_core.mixin_advanced import AdvancedMixin
from metabodesk_core.mixin_editor import EditorMixin
from metabodesk_core.mixin_network import NetworkMixin
from metabodesk_core.mixin_tools import ToolsMixin
from metabodesk_core.mixin_dialogs import DialogsMixin


class MainWindow(
    IOMixin,
    MediumMixin,
    ReactionsMixin,
    AnalysisMixin,
    AdvancedMixin,
    EditorMixin,
    NetworkMixin,
    ToolsMixin,
    DialogsMixin,
    QMainWindow,
):
    """Main application window — assembles all mixin functionality.

    The class body contains ``__init__`` (UI construction), the menu
    builder, event handlers (close, drag-and-drop), settings persistence,
    and the FBA result cache.  All domain-specific behaviour is provided
    by the nine mixin base classes listed above.
    """
    def __init__(self):
        super().__init__()
        self.setAcceptDrops(True)
        # Speed up startup by deferring repaints during UI construction
        self.setUpdatesEnabled(False)
        self.setWindowTitle("MetaboDesk")
        try:
            icon_path = _resource_path("logo.ico")
            if Path(icon_path).exists():
                self.setWindowIcon(QIcon(icon_path))
        except Exception:
            pass

        self.base_model: cobra.Model | None = None
        self.current_sbml_path: Path | None = None
        self.model_dirty: bool = False
        self._tools_proc: QProcess | None = None
        self._memote_proc: QProcess | None = None
        self._carveme_proc: QProcess | None = None
        self._carveme_last_output: str = ""
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
        self._worker: AnalysisWorker | None = None

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

        # Theme
        self._dark_theme: bool = False

        # Plugin system
        self._plugins: list = []
        self._plugin_dir = Path.home() / ".metabodesk" / "plugins"

        # FBA result cache  (key = model-state hash → result dict)
        self._fba_cache: dict[str, dict] = {}
        self._fba_cache_max = 32

        root = QWidget()
        outer = QVBoxLayout(root)

        self._build_toolbar(outer)

        self.status_lbl = QLabel("No model loaded.")
        self.status_lbl.setWordWrap(True)
        outer.addWidget(self.status_lbl)

        self.summary_lbl = QLabel("")
        self.summary_lbl.setWordWrap(True)
        outer.addWidget(self.summary_lbl)

        self.tabs = QTabWidget()
        outer.addWidget(self.tabs, stretch=1)

        # ── Lazy tab building: only build tabs on first visit ──
        # This avoids creating 15+ tabs and 8 matplotlib canvases at startup.
        # Each (label, builder) pair is registered as a placeholder.  When the
        # user clicks on the tab for the first time, the real builder fires.
        self._lazy_tab_builders: dict[int, callable] = {}
        self._lazy_tab_built: set[int] = set()

        _tab_defs = [
            ("Tools",              self._build_tools_tab),
            ("Medium",             self._build_medium_tab),
            ("Reactions",          self._build_reactions_tab),
            ("Genes",              self._build_genes_tab),
            ("Editor",             self._build_editor_tab),
            ("Knockouts",          self._build_knockout_tab),
            ("Overexpression",     self._build_overexpression_tab),
            ("FVA",                self._build_fva_tab),
            ("pFBA",               self._build_pfba_tab),
            ("Deletion Analysis",  self._build_deletion_tab),
            ("Robustness",         self._build_robustness_tab),
            ("Envelope",           self._build_envelope_tab),
            ("Sampling",           self._build_sampling_tab),
            ("Map",                self._build_network_map_tab),
            ("Results",            self._build_results_tab),
        ]

        for label, builder in _tab_defs:
            placeholder = QWidget()
            idx = self.tabs.addTab(placeholder, label)
            self._lazy_tab_builders[idx] = builder

        self.tabs.currentChanged.connect(self._on_tab_changed)
        # Eagerly build only the first visible tab (Tools)
        self._materialise_tab(0)

        self.setCentralWidget(root)
        self.resize(1700, 1120)
        self.statusBar().showMessage("Ready.")

        # Progress bar in status bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setMaximumWidth(200)
        self.progress_bar.setMaximumHeight(16)
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setVisible(False)
        self.statusBar().addPermanentWidget(self.progress_bar)

        # Use system default theme (no custom stylesheet)
        # Apply persisted UI settings before building menu
        self._init_settings()
        self._build_menu()

        # Apply persisted config from ~/.metabodesk/config.toml
        self._apply_config()

        # Auto-check tools status after UI is fully rendered
        QTimer.singleShot(500, self.tools_check_status)

        # Silent auto-update check 5 seconds after launch (generous on slow PCs)
        QTimer.singleShot(5000, lambda: self.check_for_updates(silent=True))

        # Re-enable repaints after constructing UI
        self.setUpdatesEnabled(True)

    # ---------------- (1) Unsaved changes prompt on exit ----------------

    def closeEvent(self, event):
        # Persist UI settings before exit
        try:
            self._save_config()
        except Exception:
            pass
        try:
            if getattr(self, "_memote_proc", None) is not None:
                self._force_kill_process(self._memote_proc)
                self._memote_proc = None
            if getattr(self, "_carveme_proc", None) is not None:
                self._force_kill_process(self._carveme_proc)
                self._carveme_proc = None
        except Exception:
            pass
        if self.base_model is None or not self.model_dirty:
            event.accept()
            return

        box = QMessageBox(self)
        box.setIcon(QMessageBox.Warning)
        box.setWindowTitle("Unsaved changes")
        box.setText("You have unsaved changes in the current model.\n\nSave as SBML before exiting?")

        btn_save = box.addButton("Save As...", QMessageBox.AcceptRole)
        box.addButton("Don't Save", QMessageBox.DestructiveRole)  # adds button to dialog
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

    # ---------------- Lazy tab system ----------------

    def _on_tab_changed(self, index: int):
        """Build tab content on first visit (lazy initialisation)."""
        self._materialise_tab(index)

    def _materialise_tab(self, index: int):
        """Build a tab's real content, replacing the placeholder widget."""
        if index in self._lazy_tab_built:
            return
        builder = self._lazy_tab_builders.get(index)
        if builder is None:
            return
        self._lazy_tab_built.add(index)
        # The builder creates the real widget and calls self.tabs.addTab().
        # We need to intercept that: each builder sets self.tab_xxx and
        # adds it via addTab.  We temporarily redirect addTab to replaceTab.
        old_label = self.tabs.tabText(index)
        builder()
        # The builder added a NEW tab at the end — swap it into the
        # placeholder's slot and remove the duplicate.
        new_idx = self.tabs.count() - 1
        if new_idx != index:
            new_widget = self.tabs.widget(new_idx)
            self.tabs.removeTab(new_idx)
            old_placeholder = self.tabs.widget(index)
            self.tabs.removeTab(index)
            self.tabs.insertTab(index, new_widget, old_label)
            self.tabs.setCurrentIndex(index)
            old_placeholder.deleteLater()

    def _ensure_tab_built(self, tab_name: str):
        """Ensure a specific tab is materialised by name (for programmatic access)."""
        for idx, builder in self._lazy_tab_builders.items():
            if idx not in self._lazy_tab_built and self.tabs.tabText(idx) == tab_name:
                self._materialise_tab(idx)
                return

    def _ensure_all_tabs_built(self):
        """Materialise all tabs that haven't been built yet."""
        for idx in list(self._lazy_tab_builders.keys()):
            self._materialise_tab(idx)

    # ---------------- UI helpers ----------------

    def _title(self, text: str) -> QLabel:
        """Create a bold QLabel for use as a section heading."""
        lbl = QLabel(text)
        lbl.setStyleSheet("font-weight: 600;")
        return lbl

    def _configure_table_for_excel_resize(self, table: QTableWidget):
        """Configure a QTableWidget with interactive column resize and pixel-level scrolling."""
        hdr = table.horizontalHeader()
        hdr.setSectionResizeMode(QHeaderView.Interactive)
        hdr.setStretchLastSection(False)
        hdr.setMinimumSectionSize(40)
        table.setWordWrap(False)
        table.setHorizontalScrollMode(QAbstractItemView.ScrollPerPixel)
        table.setVerticalScrollMode(QAbstractItemView.ScrollPerPixel)
        table.verticalHeader().setVisible(False)

    # ---- Unified error helper ----

    def _show_error(self, title: str, message: str, exc: Exception | None = None) -> None:
        """Unified error reporting: log + QMessageBox.

        Every mixin should call ``self._show_error(...)`` instead of
        building its own ``QMessageBox.critical`` / ``statusBar`` /
        ``logger.error`` ad-hoc.
        """
        if exc is not None:
            full = f"{message}\n\n{type(exc).__name__}: {exc}"
            logger.error("%s — %s: %s", title, message, exc)
        else:
            full = message
            logger.error("%s — %s", title, message)
        QMessageBox.critical(self, title, full)

    # ---- Extracted tab builders ----

    def _build_toolbar(self, outer: QVBoxLayout) -> None:
        """Build the top-level toolbar / controls row."""
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
            "Single Gene Deletion (SGD)", "Double Gene Deletion (DGD)", "Single Reaction Deletion (SRD)",
            "Robustness Analysis", "Production Envelope", "Flux Sampling"
        ])
        self.analysis_type.setEnabled(False)
        self.analysis_type.setToolTip("Select the analysis type to run.")

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
        self.topn_spin.setToolTip("Maximum number of items to display in charts/tables.")
        controls.addWidget(self.topn_spin)

        controls.addSpacing(15)
        controls.addWidget(QLabel("Solver:"))
        self.solver_combo = QComboBox()
        self.solver_combo.setToolTip("Optimization solver (Auto: use available).\nGurobi/CPLEX require separate licenses.")
        self.solver_combo.addItems(["Auto", "glpk", "highs", "gurobi", "cplex"])
        self.solver_combo.setEnabled(False)
        controls.addWidget(self.solver_combo, stretch=1)

        controls.addSpacing(10)
        self.loopless_chk = QCheckBox("Loopless")
        self.loopless_chk.setToolTip("Thermodynamically feasible FBA (slower, removes internal loops).")
        self.loopless_chk.setChecked(False)
        controls.addWidget(self.loopless_chk)

        controls.addSpacing(15)
        controls.addWidget(QLabel("Objective:"))
        self.objective_combo = QComboBox()
        self.objective_combo.setEnabled(False)
        self.objective_combo.setMinimumWidth(200)
        self.objective_combo.setToolTip("Select objective reaction (empty: model default).")
        self.objective_combo.currentIndexChanged.connect(self.on_objective_changed)
        self.objective_combo.setEditable(True)
        self.objective_combo.setInsertPolicy(QComboBox.NoInsert)
        self.objective_combo.setMaxVisibleItems(25)
        self.objective_combo.lineEdit().setPlaceholderText("Type to search (id/description)...")
        self.objective_combo.lineEdit().textEdited.connect(self._on_objective_text_edited)
        controls.addWidget(self.objective_combo, stretch=2)

        controls.addStretch(1)
        outer.addLayout(controls)

    def _build_tools_tab(self) -> None:
        """Build the Tools (CarveMe + memote) tab."""
        self.tab_tools = QWidget()
        tools_layout = QVBoxLayout(self.tab_tools)
        tools_layout.addWidget(self._title("Tools (CarveMe + memote)"))

        self.tools_status_lbl = QLabel("Checking tools...")
        self.tools_status_lbl.setWordWrap(True)
        tools_layout.addWidget(self.tools_status_lbl)

        btn_row = QHBoxLayout()
        self.tools_check_btn = QPushButton("Check tools")
        self.tools_check_btn.clicked.connect(self.tools_check_status)
        btn_row.addWidget(self.tools_check_btn)

        self.tools_repair_btn = QPushButton("Repair tools (reinstall from bundled wheels)")
        self.tools_repair_btn.clicked.connect(self.tools_repair)
        btn_row.addWidget(self.tools_repair_btn)
        btn_row.addStretch(1)
        tools_layout.addLayout(btn_row)

        memote_row = QHBoxLayout()
        self.memote_btn = QPushButton("Run memote report (current model)")
        self.memote_btn.setMinimumWidth(280)
        self.memote_btn.setEnabled(False)
        self.memote_btn.clicked.connect(self.run_memote_report_current_model)
        memote_row.addWidget(self.memote_btn)

        self.memote_cancel_btn = QPushButton("Cancel memote")
        self.memote_cancel_btn.setMinimumWidth(120)
        self.memote_cancel_btn.setEnabled(False)
        self.memote_cancel_btn.clicked.connect(self.cancel_memote)
        memote_row.addWidget(self.memote_cancel_btn)

        self.memote_time_lbl = QLabel("Elapsed: 00:00")
        self.memote_time_lbl.setMinimumWidth(100)
        memote_row.addWidget(self.memote_time_lbl)
        memote_row.addStretch(1)
        tools_layout.addLayout(memote_row)

        carve_row = QHBoxLayout()
        self.carveme_btn = QPushButton("Run CarveMe (genome → SBML)")
        self.carveme_btn.setMinimumWidth(280)
        self.carveme_btn.setEnabled(False)
        self.carveme_btn.clicked.connect(self.run_carveme_draft)
        carve_row.addWidget(self.carveme_btn)

        self.carveme_cancel_btn = QPushButton("Cancel CarveMe")
        self.carveme_cancel_btn.setMinimumWidth(120)
        self.carveme_cancel_btn.setEnabled(False)
        self.carveme_cancel_btn.clicked.connect(self.cancel_carveme)
        carve_row.addWidget(self.carveme_cancel_btn)

        self.carveme_time_lbl = QLabel("Elapsed: 00:00")
        self.carveme_time_lbl.setMinimumWidth(100)
        carve_row.addWidget(self.carveme_time_lbl)
        carve_row.addStretch(1)
        tools_layout.addLayout(carve_row)

        self.tools_log = QPlainTextEdit()
        self.tools_log.setReadOnly(True)
        tools_layout.addWidget(self.tools_log, stretch=1)

        self.tabs.addTab(self.tab_tools, "Tools")

    def _build_medium_tab(self) -> None:
        """Build the Medium / Exchange bounds tab."""
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

    def _build_reactions_tab(self) -> None:
        """Build the All Reactions (bounds & metadata) tab."""
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
        rxn_note = QLabel(
            "Tip: Drag column borders to resize (like Excel).\n"
            "KEGG/EC columns depend on SBML annotations; if the SBML doesn't include them, these will be empty."
        )
        rxn_note.setWordWrap(True)
        rxn_layout.addWidget(rxn_note)

        self.tabs.addTab(self.tab_rxns, "Reactions")

    def _build_genes_tab(self) -> None:
        """Build the Genes → Related reactions tab."""
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

    def _build_editor_tab(self) -> None:
        """Build the Model Editor tab with live-apply and save controls."""
        self.tab_editor = QWidget()
        editor_layout = QVBoxLayout(self.tab_editor)
        editor_layout.addWidget(self._title("Model Editor"))
        self._editor_change_count = 0

        # Live-apply info banner
        top_bar = QHBoxLayout()
        self._editor_status_lbl = QLabel(
            "ℹ️  Changes are applied instantly to the model and used in all analyses.  "
            "Your SBML file is not modified — use the Save buttons to persist."
        )
        self._editor_status_lbl.setStyleSheet(
            "color:#1565C0; font-weight:bold; padding:6px; "
            "background:#E3F2FD; border:1px solid #90CAF9; border-radius:4px; font-size:11px;"
        )
        self._editor_status_lbl.setWordWrap(True)
        top_bar.addWidget(self._editor_status_lbl, stretch=1)

        # Save button
        self._editor_save_btn = QPushButton("💾  Overwrite Original")
        self._editor_save_btn.setStyleSheet(
            "QPushButton{background:#E65100;color:white;padding:5px 12px;font-weight:bold;border-radius:4px;font-size:11px;}"
            "QPushButton:hover{background:#BF360C;}"
        )
        self._editor_save_btn.setToolTip("Overwrite the loaded SBML file")
        self._editor_save_btn.clicked.connect(lambda: self._editor_save_to_sbml(mode="overwrite"))
        top_bar.addWidget(self._editor_save_btn)

        editor_layout.addLayout(top_bar)

        self.editor_tabs = QTabWidget()
        editor_layout.addWidget(self.editor_tabs, stretch=1)

        self._build_editor_metabolite_tab()
        self._build_editor_gene_tab()
        self._build_editor_reaction_tab()
        self._build_editor_disable_delete_tab()
        self._build_editor_constraint_tab()
        self._build_editor_patch_tab()
        self._build_editor_export_sbml_tab()

        self.tabs.addTab(self.tab_editor, "Model Editor")

    def _notify_editor_change(self, action: str) -> None:
        """Update editor banner after an in-memory model change."""
        self._editor_change_count = getattr(self, "_editor_change_count", 0) + 1
        n = self._editor_change_count
        self._editor_status_lbl.setText(
            f"🔶  {n} unsaved change{'s' if n != 1 else ''}  —  Last: {action}\n"
            f"Changes are already active in analyses. Use the buttons on the right to save to file."
        )
        self._editor_status_lbl.setStyleSheet(
            "color:#E65100; font-weight:bold; padding:6px; "
            "background:#FFF3E0; border:1px solid #FFB74D; border-radius:4px; font-size:11px;"
        )
        self.statusBar().showMessage(f"✔ {action}  (in memory — SBML file unchanged)")

    def _editor_save_to_sbml(self, mode: str = "overwrite"):
        """Save the current in-memory model to SBML.

        *mode* = ``'overwrite'`` → write over the loaded file (with confirm).
        *mode* = ``'saveas'``    → always prompt for a new file path.
        """
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

        target: Path | None = None

        if mode == "overwrite":
            target = self.current_sbml_path
            if target and target.exists():
                reply = QMessageBox.warning(
                    self, "Overwrite Original File",
                    f"This will overwrite the original SBML file:\n\n"
                    f"{target}\n\n"
                    f"This cannot be undone. Continue?",
                    QMessageBox.Yes | QMessageBox.No)
                if reply != QMessageBox.Yes:
                    return
            else:
                # No original file → fall through to Save-As
                mode = "saveas"

        if mode == "saveas" or target is None:
            default_name = "metabodesk_model.xml"
            if self.current_sbml_path:
                stem = self.current_sbml_path.stem
                default_name = f"{stem}_edited.xml"
            path_str, _ = QFileDialog.getSaveFileName(
                self, "Save As New SBML File",
                str(Path.home() / default_name),
                "SBML files (*.xml *.sbml);;All files (*.*)")
            if not path_str:
                return
            target = Path(path_str)

        try:
            cobra.io.write_sbml_model(self.base_model, str(target))
            self.model_dirty = False
            self._editor_change_count = 0
            self._editor_status_lbl.setText(
                f"✅  Saved: {target.name}\n"
                f"Changes are active in analyses. Your SBML file is up to date."
            )
            self._editor_status_lbl.setStyleSheet(
                "color:#2E7D32; font-weight:bold; padding:6px; "
                "background:#E8F5E9; border:1px solid #A5D6A7; border-radius:4px; font-size:11px;"
            )
            self.statusBar().showMessage(f"Model saved: {target}")
        except Exception as e:
            self._show_error("Save Error", "Failed to save SBML file.", e)

    def _build_knockout_tab(self) -> None:
        """Build the Gene knockout tab."""
        self.tab_ko = QWidget()
        ko_layout = QVBoxLayout(self.tab_ko)

        ko_layout.addWidget(self._title("Gene knockout"))
        self.gene_search = QLineEdit()
        self.gene_search.setPlaceholderText("Search genes...")
        self.gene_search.textChanged.connect(self.filter_gene_list)
        ko_layout.addWidget(self.gene_search)

        self.gene_list = QListWidget()
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

    def _build_overexpression_tab(self) -> None:
        """Build the Overexpression (reaction bound scaling) tab."""
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

    def _build_fva_tab(self) -> None:
        """Build the Flux Variability Analysis (FVA) tab."""
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

        fva_params_row = QHBoxLayout()
        fva_params_row.addWidget(QLabel("Fraction of optimum:"))
        self.fva_fraction_spin = QDoubleSpinBox()
        self.fva_fraction_spin.setRange(0.0, 1.0)
        self.fva_fraction_spin.setSingleStep(0.05)
        self.fva_fraction_spin.setDecimals(2)
        self.fva_fraction_spin.setValue(1.0)
        self.fva_fraction_spin.setToolTip("Fraction of optimal objective required (1.0 = optimal only, 0.9 = 90%+ optimal).")
        fva_params_row.addWidget(self.fva_fraction_spin)
        fva_params_row.addSpacing(15)
        fva_params_row.addWidget(QLabel("Processes:"))
        self.fva_processes_spin = QSpinBox()
        self.fva_processes_spin.setRange(1, max(os.cpu_count() or 1, 1))
        self.fva_processes_spin.setValue(1)
        self.fva_processes_spin.setToolTip("Number of parallel processes for FVA.\nHigher = faster but more memory.\nUse 1 if running as frozen .exe.")
        fva_params_row.addWidget(self.fva_processes_spin)
        fva_params_row.addStretch(1)
        fva_layout.addLayout(fva_params_row)

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

    def _build_pfba_tab(self) -> None:
        """Build the Parsimonious FBA (pFBA) tab."""
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

    def _build_deletion_tab(self) -> None:
        """Build the Single Gene/Reaction Deletion Analysis tab."""
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

        self.deletion_filter_optimal = QCheckBox("Optimal only")
        self.deletion_filter_optimal.setToolTip("Show only items with optimal (non-zero) growth after deletion")
        self.deletion_filter_optimal.stateChanged.connect(self.filter_deletion_table)
        deletion_search_row.addWidget(self.deletion_filter_optimal)

        self.deletion_filter_nonoptimal = QCheckBox("Non-optimal only")
        self.deletion_filter_nonoptimal.setToolTip("Show only items with zero or near-zero growth (lethal knockouts)")
        self.deletion_filter_nonoptimal.stateChanged.connect(self.filter_deletion_table)
        deletion_search_row.addWidget(self.deletion_filter_nonoptimal)

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
        self.sgd_export_xlsx_btn = QPushButton("Export SGD (Excel)")
        self.sgd_export_xlsx_btn.clicked.connect(self.export_sgd_excel)
        del_btn_row.addWidget(self.sgd_export_xlsx_btn)
        self.dgd_export_btn = QPushButton("Export DGD (CSV)")
        self.dgd_export_btn.clicked.connect(self.export_dgd_csv)
        del_btn_row.addWidget(self.dgd_export_btn)
        self.dgd_export_xlsx_btn = QPushButton("Export DGD (Excel)")
        self.dgd_export_xlsx_btn.clicked.connect(self.export_dgd_excel)
        del_btn_row.addWidget(self.dgd_export_xlsx_btn)
        self.srd_export_btn = QPushButton("Export SRD (CSV)")
        self.srd_export_btn.clicked.connect(self.export_srd_csv)
        del_btn_row.addWidget(self.srd_export_btn)
        self.srd_export_xlsx_btn = QPushButton("Export SRD (Excel)")
        self.srd_export_xlsx_btn.clicked.connect(self.export_srd_excel)
        del_btn_row.addWidget(self.srd_export_xlsx_btn)
        del_btn_row.addStretch(1)
        deletion_layout.addLayout(del_btn_row)

        self.tabs.addTab(self.tab_deletion, "Deletion Analysis")

    def _build_robustness_tab(self) -> None:
        """Build the Robustness Analysis tab."""
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
        self.robustness_combo.currentIndexChanged.connect(self._on_robustness_combo_changed)
        robustness_search_row.addWidget(self.robustness_combo, stretch=2)
        robustness_search_row.addStretch(1)
        robustness_layout.addLayout(robustness_search_row)

        self.robustness_canvas = MplCanvas(self, width=10, height=4, dpi=100)
        robustness_layout.addWidget(self.robustness_canvas, stretch=1)

        self.robustness_tbl = QTableWidget(0, 2)
        self.robustness_tbl.setHorizontalHeaderLabels(["Bound Value", "Objective Value"])
        self.robustness_tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.robustness_tbl.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._configure_table_for_excel_resize(self.robustness_tbl)
        robustness_layout.addWidget(self.robustness_tbl, stretch=1)

        self.robustness_info_lbl = QLabel("Robustness analysis: (none)")
        self.robustness_info_lbl.setWordWrap(True)
        robustness_layout.addWidget(self.robustness_info_lbl)

        robustness_btn_row = QHBoxLayout()
        self.robustness_export_csv_btn = QPushButton("Export Robustness (CSV)")
        self.robustness_export_csv_btn.clicked.connect(self.export_robustness_csv)
        robustness_btn_row.addWidget(self.robustness_export_csv_btn)
        self.robustness_export_xlsx_btn = QPushButton("Export Robustness (Excel)")
        self.robustness_export_xlsx_btn.clicked.connect(self.export_robustness_excel)
        robustness_btn_row.addWidget(self.robustness_export_xlsx_btn)
        robustness_btn_row.addStretch(1)
        robustness_layout.addLayout(robustness_btn_row)

        self.tabs.addTab(self.tab_robustness, "Robustness")

    def _build_envelope_tab(self) -> None:
        """Build the Production Envelope (Growth vs Product) tab."""
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
        envelope_search_row.addWidget(QLabel("Reaction filter:"))
        self.envelope_search = QLineEdit()
        self.envelope_search.setPlaceholderText("Type to filter reaction list...")
        self.envelope_search.textChanged.connect(self._filter_envelope_combo)
        envelope_search_row.addWidget(self.envelope_search, stretch=2)
        envelope_search_row.addWidget(QLabel("Select reaction:"))
        self.envelope_combo = QComboBox()
        self.envelope_combo.setToolTip("Pick a reaction from the filtered list to set as the product reaction")
        self.envelope_combo.currentIndexChanged.connect(self._on_envelope_combo_changed)
        envelope_search_row.addWidget(self.envelope_combo, stretch=2)
        envelope_search_row.addStretch(1)
        envelope_layout.addLayout(envelope_search_row)

        self.envelope_canvas = MplCanvas(self, width=10, height=4, dpi=100)
        envelope_layout.addWidget(self.envelope_canvas, stretch=1)

        self.envelope_tbl = QTableWidget(0, 2)
        self.envelope_tbl.setHorizontalHeaderLabels(["Product Bound", "Growth Rate"])
        self.envelope_tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.envelope_tbl.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._configure_table_for_excel_resize(self.envelope_tbl)
        envelope_layout.addWidget(self.envelope_tbl, stretch=1)

        self.envelope_info_lbl = QLabel("Production envelope: (none)")
        self.envelope_info_lbl.setWordWrap(True)
        envelope_layout.addWidget(self.envelope_info_lbl)

        envelope_btn_row = QHBoxLayout()
        self.envelope_export_csv_btn = QPushButton("Export Envelope (CSV)")
        self.envelope_export_csv_btn.clicked.connect(self.export_envelope_csv)
        envelope_btn_row.addWidget(self.envelope_export_csv_btn)
        self.envelope_export_xlsx_btn = QPushButton("Export Envelope (Excel)")
        self.envelope_export_xlsx_btn.clicked.connect(self.export_envelope_excel)
        envelope_btn_row.addWidget(self.envelope_export_xlsx_btn)
        envelope_btn_row.addStretch(1)
        envelope_layout.addLayout(envelope_btn_row)

        self.tabs.addTab(self.tab_envelope, "Envelope")

    def _build_sampling_tab(self) -> None:
        """Build the Flux Sampling tab."""
        self.tab_sampling = QWidget()
        sampling_layout = QVBoxLayout(self.tab_sampling)

        sampling_layout.addWidget(self._title("Flux Sampling"))
        sampling_help = QLabel(
            "Monte Carlo sampling of the feasible flux space. "
            "ACHR (Artificially Centered Hit-and-Run) or OptGP (Optima-Guided Parallel, faster for large models). "
            "Returns mean, standard deviation, and min/max ranges for each reaction."
        )
        sampling_help.setWordWrap(True)
        sampling_layout.addWidget(sampling_help)

        sampling_param_row = QHBoxLayout()
        sampling_param_row.addWidget(QLabel("Sampler:"))
        self.sampler_type_combo = QComboBox()
        self.sampler_type_combo.addItems(["ACHR", "OptGP"])
        self.sampler_type_combo.setToolTip("ACHR: standard sampler. OptGP: faster parallel sampler for large models.")
        sampling_param_row.addWidget(self.sampler_type_combo)
        sampling_param_row.addSpacing(15)
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

    def _build_network_map_tab(self) -> None:
        """Build the Network Map visualization tab."""
        self.tab_map = QWidget()
        map_layout = QVBoxLayout(self.tab_map)

        map_layout.addWidget(self._title("Network Map"))
        map_help = QLabel(
            "Visualize metabolic network topology. Nodes: reactions (circles) and metabolites (squares). "
            "Colors reflect flux magnitude or gene knockout impact. Filter by reaction/metabolite search and neighborhood depth."
        )
        map_help.setWordWrap(True)
        map_layout.addWidget(map_help)

        map_filter_row = QHBoxLayout()
        map_filter_row.addWidget(QLabel("Search:"))
        self.map_search = QLineEdit()
        self.map_search.setPlaceholderText("reaction/metabolite id contains...")
        map_filter_row.addWidget(self.map_search, stretch=2)

        map_filter_row.addWidget(QLabel("Depth:"))
        self.map_depth = QSpinBox()
        self.map_depth.setRange(1, 5)
        self.map_depth.setValue(2)
        self.map_depth.setToolTip("Neighborhood depth (1-5)")
        map_filter_row.addWidget(self.map_depth)

        map_filter_row.addWidget(QLabel("Flux threshold:"))
        self.map_threshold = QDoubleSpinBox()
        self.map_threshold.setRange(0, 100)
        self.map_threshold.setValue(0.1)
        self.map_threshold.setSingleStep(0.1)
        self.map_threshold.setToolTip("Show only reactions with |flux| > threshold")
        map_filter_row.addWidget(self.map_threshold)

        map_filter_row.addWidget(QLabel("View:"))
        self.map_view_combo = QComboBox()
        self.map_view_combo.addItems(["Flux magnitude", "Gene KO Impact", "FVA Bounds"])
        map_filter_row.addWidget(self.map_view_combo)

        self.map_render_btn = QPushButton("Render Map")
        self.map_render_btn.clicked.connect(self._render_network_map)
        map_filter_row.addWidget(self.map_render_btn)

        self.map_help_btn = QPushButton("How to use map?")
        self.map_help_btn.clicked.connect(self._show_map_help)
        map_filter_row.addWidget(self.map_help_btn)
        map_filter_row.addStretch(1)
        map_layout.addLayout(map_filter_row)

        map_opts_row1 = QHBoxLayout()
        map_opts_row1.addWidget(QLabel("Max nodes:"))
        self.map_max_nodes = QSpinBox()
        self.map_max_nodes.setRange(100, 50000)
        self.map_max_nodes.setValue(1000)
        self.map_max_nodes.setToolTip("Limit displayed nodes (lower = faster rendering)")
        map_opts_row1.addWidget(self.map_max_nodes)

        map_opts_row1.addWidget(QLabel("Max edges:"))
        self.map_max_edges = QSpinBox()
        self.map_max_edges.setRange(0, 200000)
        self.map_max_edges.setValue(0)
        self.map_max_edges.setToolTip("0 = no edge limit")
        map_opts_row1.addWidget(self.map_max_edges)

        map_opts_row1.addWidget(QLabel("Min degree:"))
        self.map_min_degree = QSpinBox()
        self.map_min_degree.setRange(0, 20)
        self.map_min_degree.setValue(0)
        self.map_min_degree.setToolTip("Hide nodes with degree < min")
        map_opts_row1.addWidget(self.map_min_degree)

        map_opts_row1.addWidget(QLabel("Top-N rxns:"))
        self.map_topn_rxn = QSpinBox()
        self.map_topn_rxn.setRange(0, 5000)
        self.map_topn_rxn.setValue(0)
        self.map_topn_rxn.setToolTip("0 = off")
        map_opts_row1.addWidget(self.map_topn_rxn)

        map_opts_row1.addWidget(QLabel("Layout:"))
        self.map_layout_combo = QComboBox()
        self.map_layout_combo.addItems(["Spring", "Kamada-Kawai", "Circular"])
        map_opts_row1.addWidget(self.map_layout_combo)

        self.map_only_connected_chk = QCheckBox("Only connected to search")
        map_opts_row1.addWidget(self.map_only_connected_chk)

        self.map_hide_orphans_chk = QCheckBox("Hide orphan metabolites")
        map_opts_row1.addWidget(self.map_hide_orphans_chk)

        self.map_show_legend_chk = QCheckBox("Show legend")
        self.map_show_legend_chk.setChecked(True)
        map_opts_row1.addWidget(self.map_show_legend_chk)

        map_opts_row1.addStretch(1)
        map_layout.addLayout(map_opts_row1)

        map_opts_row2 = QHBoxLayout()
        self.map_exchange_only_chk = QCheckBox("Exchange only")
        map_opts_row2.addWidget(self.map_exchange_only_chk)

        self.map_objective_only_chk = QCheckBox("Objective neighborhood only")
        map_opts_row2.addWidget(self.map_objective_only_chk)

        map_opts_row2.addWidget(QLabel("Subsystem:"))
        self.map_subsystem_combo = QComboBox()
        self.map_subsystem_combo.addItem("All subsystems")
        map_opts_row2.addWidget(self.map_subsystem_combo)

        self.map_focus_btn = QPushButton("Focus selected")
        self.map_focus_btn.clicked.connect(self._focus_on_selected_map_node)
        map_opts_row2.addWidget(self.map_focus_btn)

        self.map_export_img_btn = QPushButton("Export map image")
        self.map_export_img_btn.clicked.connect(self._export_map_image)
        map_opts_row2.addWidget(self.map_export_img_btn)

        self.map_export_csv_btn = QPushButton("Export graph CSV")
        self.map_export_csv_btn.clicked.connect(self._export_map_csv)
        map_opts_row2.addWidget(self.map_export_csv_btn)

        map_opts_row2.addStretch(1)
        map_layout.addLayout(map_opts_row2)

        self.map_canvas = MplCanvas(self, width=12, height=7, dpi=100)
        map_layout.addWidget(self.map_canvas, stretch=1)
        self.map_canvas.mpl_connect('button_press_event', self._on_map_click)
        self.map_canvas.mpl_connect('motion_notify_event', self._on_map_hover)
        self._set_map_placeholder("Click 'Render Map' to generate the network.")
        self._map_colorbar = None
        self._map_tooltip = None

        self.map_toolbar = _get_nav_toolbar()(self.map_canvas, self.tab_map)
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

    def _build_results_tab(self) -> None:
        """Build the Results tab (FBA comparison + export)."""
        self.tab_results = QWidget()
        results_layout = QVBoxLayout(self.tab_results)

        results_layout.addWidget(self._title("Results"))

        compare_row = QHBoxLayout()
        compare_row.addWidget(QLabel("Compare:"))
        self.compare_mode = QComboBox()
        self.compare_mode.addItems(["Gene knockout only", "Overexpression only", "Original vs Baseline"])
        self.compare_mode.currentTextChanged.connect(self._refresh_results_view_from_last_run)
        compare_row.addWidget(self.compare_mode)

        compare_row.addWidget(QLabel("Chart:"))
        self.results_chart_type = QComboBox()
        self.results_chart_type.addItems(["Bar comparison", "Flux histogram", "Waterfall (top changes)"])
        self.results_chart_type.currentTextChanged.connect(self._refresh_results_view_from_last_run)
        compare_row.addWidget(self.results_chart_type)

        compare_row.addStretch(1)
        results_layout.addLayout(compare_row)

        self.results_lbl = QLabel("Objective: -")
        self.results_lbl.setWordWrap(True)
        results_layout.addWidget(self.results_lbl)

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

        results_opts_row = QHBoxLayout()
        self.results_nonzero_chk = QCheckBox("Non-zero flux only")
        self.results_nonzero_chk.stateChanged.connect(self.filter_flux_table)
        results_opts_row.addWidget(self.results_nonzero_chk)

        self.results_changed_chk = QCheckBox("Changed only (Δ≠0)")
        self.results_changed_chk.stateChanged.connect(self.filter_flux_table)
        results_opts_row.addWidget(self.results_changed_chk)

        results_opts_row.addStretch(1)

        self.results_export_csv_btn = QPushButton("Export CSV")
        self.results_export_csv_btn.clicked.connect(self._export_results_csv)
        results_opts_row.addWidget(self.results_export_csv_btn)

        self.results_export_excel_btn = QPushButton("Export Excel")
        self.results_export_excel_btn.clicked.connect(self._export_results_excel)
        results_opts_row.addWidget(self.results_export_excel_btn)

        self.results_export_chart_btn = QPushButton("Export Chart")
        self.results_export_chart_btn.clicked.connect(self._export_results_chart)
        results_opts_row.addWidget(self.results_export_chart_btn)

        self.results_export_latex_btn = QPushButton("Export LaTeX")
        self.results_export_latex_btn.setToolTip("Export table as LaTeX for academic papers")
        self.results_export_latex_btn.clicked.connect(self._export_results_latex)
        results_opts_row.addWidget(self.results_export_latex_btn)

        results_layout.addLayout(results_opts_row)

        self.canvas = MplCanvas(self, width=10, height=3.8, dpi=100)
        results_layout.addWidget(self.canvas, stretch=1)

        self.flux_tbl = QTableWidget(0, 4)
        self.flux_tbl.setHorizontalHeaderLabels(["Reaction", "Baseline flux", "Compared flux", "Delta"])
        self.flux_tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.flux_tbl.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._configure_table_for_excel_resize(self.flux_tbl)
        results_layout.addWidget(self.flux_tbl, stretch=1)

        # ── All Fluxes Browser ──────────────────────────────────
        allflux_header = QHBoxLayout()
        allflux_header.addWidget(self._title("All Reaction Fluxes"))
        self.allflux_toggle_btn = QPushButton("▼ Show")
        self.allflux_toggle_btn.setFixedWidth(90)
        self.allflux_toggle_btn.clicked.connect(self._toggle_allflux_panel)
        allflux_header.addWidget(self.allflux_toggle_btn)
        allflux_header.addStretch(1)
        results_layout.addLayout(allflux_header)

        # Container widget for collapsible section
        self.allflux_container = QWidget()
        allflux_inner = QVBoxLayout(self.allflux_container)
        allflux_inner.setContentsMargins(0, 0, 0, 0)

        allflux_filter_row = QHBoxLayout()
        allflux_filter_row.addWidget(QLabel("Search:"))
        self.allflux_search = QLineEdit()
        self.allflux_search.setPlaceholderText("Filter by reaction ID or name...")
        self.allflux_search.textChanged.connect(self._filter_allflux_table)
        allflux_filter_row.addWidget(self.allflux_search, stretch=2)

        allflux_filter_row.addWidget(QLabel("Min |flux|:"))
        self.allflux_min_abs = QDoubleSpinBox()
        self.allflux_min_abs.setRange(0.0, 1e12)
        self.allflux_min_abs.setDecimals(6)
        self.allflux_min_abs.setSingleStep(0.01)
        self.allflux_min_abs.setValue(0.0)
        self.allflux_min_abs.valueChanged.connect(self._filter_allflux_table)
        allflux_filter_row.addWidget(self.allflux_min_abs)

        self.allflux_nonzero_chk = QCheckBox("Non-zero only")
        self.allflux_nonzero_chk.stateChanged.connect(self._filter_allflux_table)
        allflux_filter_row.addWidget(self.allflux_nonzero_chk)

        allflux_filter_row.addWidget(QLabel("Subsystem:"))
        self.allflux_subsystem_combo = QComboBox()
        self.allflux_subsystem_combo.addItem("All")
        self.allflux_subsystem_combo.currentTextChanged.connect(self._filter_allflux_table)
        allflux_filter_row.addWidget(self.allflux_subsystem_combo)

        allflux_filter_row.addStretch(1)
        allflux_inner.addLayout(allflux_filter_row)

        self.allflux_info_lbl = QLabel("Run FBA to see all reaction fluxes.")
        allflux_inner.addWidget(self.allflux_info_lbl)

        self.allflux_tbl = QTableWidget(0, 5)
        self.allflux_tbl.setHorizontalHeaderLabels(
            ["Reaction ID", "Name", "Flux", "Bounds [LB, UB]", "Subsystem"]
        )
        self.allflux_tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.allflux_tbl.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.allflux_tbl.setSortingEnabled(True)
        self._configure_table_for_excel_resize(self.allflux_tbl)
        allflux_inner.addWidget(self.allflux_tbl, stretch=1)

        allflux_export_row = QHBoxLayout()
        self.allflux_export_csv_btn = QPushButton("Export All Fluxes (CSV)")
        self.allflux_export_csv_btn.clicked.connect(self._export_allflux_csv)
        allflux_export_row.addWidget(self.allflux_export_csv_btn)
        self.allflux_export_xlsx_btn = QPushButton("Export All Fluxes (Excel)")
        self.allflux_export_xlsx_btn.clicked.connect(self._export_allflux_excel)
        allflux_export_row.addWidget(self.allflux_export_xlsx_btn)
        allflux_export_row.addStretch(1)
        allflux_inner.addLayout(allflux_export_row)

        self.allflux_container.setVisible(False)
        results_layout.addWidget(self.allflux_container, stretch=1)

        # ── Extra buttons row ──────────────────────────────────
        results_extra_row = QHBoxLayout()
        self.results_shadow_btn = QPushButton("Shadow Prices")
        self.results_shadow_btn.clicked.connect(self._show_shadow_prices)
        results_extra_row.addWidget(self.results_shadow_btn)

        self.results_reduced_btn = QPushButton("Reduced Costs")
        self.results_reduced_btn.clicked.connect(self._show_reduced_costs)
        results_extra_row.addWidget(self.results_reduced_btn)

        self.results_essential_btn = QPushButton("Find Essential Reactions")
        self.results_essential_btn.clicked.connect(self._find_essential_reactions)
        results_extra_row.addWidget(self.results_essential_btn)

        self.results_fva_btn = QPushButton("Quick FVA")
        self.results_fva_btn.clicked.connect(self._quick_fva)
        results_extra_row.addWidget(self.results_fva_btn)

        results_extra_row.addStretch(1)
        results_layout.addLayout(results_extra_row)

        self.tabs.addTab(self.tab_results, "Results")

    # ---------------- Theme & Settings ----------------

    def _init_settings(self):
        """Backward-compat stub — all settings are now in TOML config."""
        pass

    def _apply_config(self):
        """Load ``~/.metabodesk/config.toml`` and apply persisted values to UI widgets.

        Uses ``getattr`` for widgets that may not yet exist (lazy tabs).
        The settings will be re-applied when the tab is materialised via
        _save_config / _apply_config cycle, or the config object retains
        the value for _save_config to read back.
        """
        try:
            self._app_config = AppConfig.load()
            cfg = self._app_config

            # Solver (always exists — in toolbar)
            idx = self.solver_combo.findText(cfg.solver)
            if idx >= 0:
                self.solver_combo.setCurrentIndex(idx)
            self.topn_spin.setValue(cfg.top_n)
            self.loopless_chk.setChecked(cfg.loopless)

            # Widgets that may not exist yet (lazy tabs) ──
            def _set_spin(attr, val):
                w = getattr(self, attr, None)
                if w is not None:
                    w.setValue(val)

            def _set_combo_text(attr, val):
                w = getattr(self, attr, None)
                if w is not None:
                    idx2 = w.findText(val)
                    if idx2 >= 0:
                        w.setCurrentIndex(idx2)

            _set_spin("fva_fraction_spin", cfg.fva_fraction)
            _set_spin("fva_processes_spin", cfg.fva_processes)
            _set_spin("sampling_size", cfg.sampling_size)
            _set_combo_text("sampler_type_combo", cfg.sampler)
            _set_spin("map_max_nodes", cfg.map_max_nodes)
            _set_combo_text("map_layout_combo", cfg.map_layout)
            _set_spin("map_depth", cfg.map_depth)
            _set_spin("map_threshold", cfg.map_threshold)

            # Non-widget settings
            self.max_undo_steps = cfg.undo_max
            self._fba_cache_max = cfg.cache_max

            # Theme
            if cfg.theme == "dark":
                self._dark_theme = True
                from metabodesk_core.constants import DARK_STYLESHEET
                self.setStyleSheet(DARK_STYLESHEET)
            elif cfg.theme == "light":
                self._dark_theme = False
                self.setStyleSheet("")

            logger.info("Config applied from %s", self._app_config._path)
        except Exception as exc:
            logger.warning("Could not apply config: %s", exc)
            self._app_config = AppConfig()

    def _save_config(self):
        """Persist current UI settings to ``~/.metabodesk/config.toml``."""
        try:
            cfg = self._app_config
            cfg.solver = self.solver_combo.currentText()
            cfg.top_n = self.topn_spin.value()
            cfg.loopless = self.loopless_chk.isChecked()

            def _read_spin(attr, default):
                w = getattr(self, attr, None)
                return w.value() if w is not None else default

            def _read_combo(attr, default):
                w = getattr(self, attr, None)
                return w.currentText() if w is not None else default

            cfg.fva_fraction = _read_spin("fva_fraction_spin", cfg.fva_fraction)
            cfg.fva_processes = _read_spin("fva_processes_spin", cfg.fva_processes)
            cfg.sampling_size = _read_spin("sampling_size", cfg.sampling_size)
            cfg.sampler = _read_combo("sampler_type_combo", cfg.sampler)
            cfg.map_max_nodes = _read_spin("map_max_nodes", cfg.map_max_nodes)
            cfg.map_layout = _read_combo("map_layout_combo", cfg.map_layout)
            cfg.map_depth = _read_spin("map_depth", cfg.map_depth)
            cfg.map_threshold = _read_spin("map_threshold", cfg.map_threshold)
            cfg.undo_max = self.max_undo_steps
            cfg.cache_max = self._fba_cache_max
            cfg.theme = "dark" if self._dark_theme else "system"
            cfg.save()
            logger.info("Config saved.")
        except Exception as exc:
            logger.warning("Could not save config: %s", exc)

    def _init_persisted_controls(self):
        """No-op — superseded by ``_apply_config()`` which reads TOML config."""
        pass

    # Placeholder - theme functions removed

    def _light_stylesheet(self) -> str:
        """Return the light-theme Qt stylesheet string."""
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
        """Return the dark-theme Qt stylesheet string."""
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
        """Show a TextPopup dialog displaying the full content of a table cell."""
        if not item:
            return
        dlg = TextPopup(title, item.text(), self)
        dlg.exec()

    def _on_reaction_cell_double_clicked(self, item: QTableWidgetItem | None):
        """Handle double-click on reaction table: open KEGG/EC URL or show text popup."""
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
        """Copy the selected medium-table cell text into the details pane."""
        self.medium_details.setPlainText(current.text() if current else "")

    def _update_reaction_details_from_current_cell(self, current: QTableWidgetItem, _prev: QTableWidgetItem):
        """Copy the selected reaction-table cell text into the details pane."""
        self.reaction_details.setPlainText(current.text() if current else "")

    def _update_gene_rxn_details_from_current_cell(self, current: QTableWidgetItem, _prev: QTableWidgetItem):
        """Copy the selected gene-reaction table cell text into the details pane."""
        self.gene_rxn_details.setPlainText(current.text() if current else "")

    # ---------------- (4) Flux table filtering ----------------

    def filter_flux_table(self):
        text = (self.flux_search.text() or "").strip().lower()
        thr = float(self.delta_thresh.value())
        nonzero_only = self.results_nonzero_chk.isChecked()
        changed_only = self.results_changed_chk.isChecked()

        for row in range(self.flux_tbl.rowCount()):
            rid = (self.flux_tbl.item(row, 0).text() if self.flux_tbl.item(row, 0) else "").lower()
            try:
                baseline = float((self.flux_tbl.item(row, 1).text() if self.flux_tbl.item(row, 1) else "0").strip())
            except Exception:
                baseline = 0.0
            try:
                compared = float((self.flux_tbl.item(row, 2).text() if self.flux_tbl.item(row, 2) else "0").strip())
            except Exception:
                compared = 0.0
            try:
                d = float((self.flux_tbl.item(row, 3).text() if self.flux_tbl.item(row, 3) else "0").strip())
            except Exception:
                d = 0.0

            hide = False
            if text and text not in rid:
                hide = True
            if abs(d) < thr:
                hide = True
            if nonzero_only and baseline == 0 and compared == 0:
                hide = True
            if changed_only and d == 0:
                hide = True
            self.flux_tbl.setRowHidden(row, hide)

    # ---------------- Objective search ----------------

    def _setup_objective_completer(self):
        """Attach a QCompleter to the objective combo box for type-to-search."""
        texts = [self.objective_combo.itemText(i) for i in range(self.objective_combo.count())]
        comp = QCompleter(texts, self)
        comp.setCaseSensitivity(Qt.CaseInsensitive)
        comp.setFilterMode(Qt.MatchContains)
        comp.setCompletionMode(QCompleter.PopupCompletion)
        self.objective_combo.setCompleter(comp)

    def _on_objective_text_edited(self, text: str):
        """Filter objective combo entries as the user types, selecting the first match."""
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
        """Set the model objective to the reaction selected in the objective combo."""
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

    def set_busy(self, busy: bool, message: str = ""):
        self.is_running = busy
        _widget_names = [
            "open_btn", "run_btn", "close_uptakes_btn",
            "reset_medium_btn", "reset_reactions_btn",
            "analysis_type", "topn_spin", "objective_combo",
            "medium_search", "medium_table",
            "rxn_global_search", "met_filter", "met_filter_mode", "reaction_table",
            "clear_rxn_overrides_btn", "remove_selected_overrides_btn",
            "gene_global_search", "genes_tab_list", "gene_rxn_table",
            "gene_search", "gene_list", "add_ko_btn", "clear_ko_btn",
            "rxn_search", "rxn_list", "oe_factor", "oe_scale_lb",
            "add_oe_btn", "remove_oe_btn", "clear_oe_btn",
            "temp_ub_spin", "set_temp_ub_btn", "remove_temp_ub_btn", "clear_temp_ub_btn",
            "oe_active_list", "temp_active_list",
            "fva_flux_search", "pfba_flux_search",
            "fva_fraction_spin", "fva_processes_spin",
            "deletion_search",
            "robustness_rxn", "robustness_min", "robustness_max", "robustness_steps",
            "robustness_bound_ub", "robustness_bound_lb",
            "envelope_product", "envelope_steps",
            "sampling_size", "sampling_search",
            "sampler_type_combo",
            "loopless_chk",
            "map_search", "map_depth", "map_threshold", "map_view_combo",
            "map_max_nodes", "map_max_edges", "map_min_degree", "map_topn_rxn",
            "map_layout_combo", "map_only_connected_chk", "map_hide_orphans_chk",
            "map_show_legend_chk", "map_exchange_only_chk", "map_objective_only_chk",
            "map_subsystem_combo",
            "map_focus_btn", "map_export_img_btn", "map_export_csv_btn",
            "compare_mode",
            "editor_tabs",
            "tools_check_btn", "tools_repair_btn",
            "memote_btn", "memote_cancel_btn", "carveme_btn", "carveme_cancel_btn",
        ]
        for name in _widget_names:
            w = getattr(self, name, None)
            if w is None:
                continue
            if name == "memote_cancel_btn":
                enabled = (self._memote_proc is not None)
            elif name == "carveme_cancel_btn":
                enabled = (self._carveme_proc is not None)
            else:
                enabled = not busy
            if name == "run_btn":
                enabled = enabled and (self.base_model is not None)
            if name == "analysis_type":
                enabled = enabled and (self.base_model is not None)
            if name == "memote_btn":
                enabled = enabled and (self.base_model is not None) and (self._memote_proc is None)
            w.setEnabled(enabled)
        if message:
            self.status_lbl.setText(message)
            self.statusBar().showMessage(message)

    # ---------------- Menu ----------------

    def _build_menu(self):
        """Build the application menu bar (File, Run, Medium, Model, Help)."""
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
        act_export_results.triggered.connect(self._export_results_csv)
        self.file_menu.addAction(act_export_results)

        act_export_excel = QAction("Export results (Excel)...", self)
        act_export_excel.triggered.connect(self._export_results_excel)
        self.file_menu.addAction(act_export_excel)

        act_export_json = QAction("Export results (JSON)...", self)
        act_export_json.triggered.connect(self.export_results_json)
        self.file_menu.addAction(act_export_json)

        act_export_all = QAction("Export All (scenario+results+chart)...", self)
        act_export_all.triggered.connect(self.export_all)
        self.file_menu.addAction(act_export_all)

        act_export_jupyter = QAction("Export Jupyter Notebook...", self)
        act_export_jupyter.triggered.connect(self.export_jupyter_notebook)
        self.file_menu.addAction(act_export_jupyter)

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

        act_community = QAction("Community Analysis...", self)
        act_community.triggered.connect(self.run_community_analysis)
        run_menu.addAction(act_community)

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

        model_menu.addSeparator()

        act_model_stats = QAction("Model Statistics Dashboard...", self)
        act_model_stats.triggered.connect(self.show_model_stats_dashboard)
        model_menu.addAction(act_model_stats)

        act_connectivity = QAction("Metabolite Connectivity Analysis...", self)
        act_connectivity.triggered.connect(self.show_metabolite_connectivity)
        model_menu.addAction(act_connectivity)

        act_annotation = QAction("Annotation Quality Score...", self)
        act_annotation.triggered.connect(self.show_annotation_quality)
        model_menu.addAction(act_annotation)

        act_bigg = QAction("BiGG Database Browser...", self)
        act_bigg.triggered.connect(self.show_bigg_browser)
        model_menu.addAction(act_bigg)

        model_menu.addSeparator()

        act_add_constraint = QAction("Add Custom Constraint...", self)
        act_add_constraint.triggered.connect(self.add_custom_constraint)
        model_menu.addAction(act_add_constraint)

        act_html_report = QAction("Generate HTML Report...", self)
        act_html_report.triggered.connect(self.generate_html_report)
        model_menu.addAction(act_html_report)

        model_menu.addSeparator()

        act_met_summary = QAction("Metabolite Summary...", self)
        act_met_summary.triggered.connect(self.show_metabolite_summary)
        model_menu.addAction(act_met_summary)

        act_flux_compare = QAction("Compare Flux Distributions...", self)
        act_flux_compare.triggered.connect(self.show_flux_distribution_comparison)
        model_menu.addAction(act_flux_compare)

        act_shortcuts = QAction("Keyboard shortcuts", self)
        act_shortcuts.triggered.connect(self.show_shortcuts)
        help_menu.addAction(act_shortcuts)

        act_tutorial = QAction("🎓  Interactive Tutorial...", self)
        act_tutorial.triggered.connect(self.show_tutorial_wizard)
        help_menu.addAction(act_tutorial)

        act_about = QAction("About MetaboDesk", self)
        act_about.triggered.connect(self.show_about)
        help_menu.addAction(act_about)

        act_update = QAction("Check for Updates...", self)
        act_update.triggered.connect(lambda: self.check_for_updates(silent=False))
        help_menu.addAction(act_update)

        # ---- View menu ----
        view_menu = menubar.addMenu("&View")

        act_toggle_theme = QAction("Toggle Dark/Light Theme", self)
        act_toggle_theme.setShortcut(QKeySequence("Ctrl+T"))
        act_toggle_theme.triggered.connect(self.toggle_theme)
        view_menu.addAction(act_toggle_theme)

        # ---- Advanced Analysis submenu ----
        run_menu.addSeparator()
        advanced_menu = run_menu.addMenu("Advanced Analysis")

        act_gene_expr = QAction("Gene Expression Integration...", self)
        act_gene_expr.triggered.connect(self.run_gene_expression_analysis)
        advanced_menu.addAction(act_gene_expr)

        act_optknock = QAction("OptKnock / Strain Design...", self)
        act_optknock.triggered.connect(self.run_optknock)
        advanced_menu.addAction(act_optknock)

        act_dfba = QAction("Dynamic FBA (dFBA)...", self)
        act_dfba.triggered.connect(self.run_dfba)
        advanced_menu.addAction(act_dfba)

        act_flux_coupling = QAction("Flux Coupling Analysis...", self)
        act_flux_coupling.triggered.connect(self.run_flux_coupling)
        advanced_menu.addAction(act_flux_coupling)

        act_tmfa = QAction("Thermodynamic FBA (TMFA)...", self)
        act_tmfa.triggered.connect(self.run_tmfa)
        advanced_menu.addAction(act_tmfa)

        act_gap_fill = QAction("Gap-Filling (GapFind/GapFill)...", self)
        act_gap_fill.triggered.connect(self.run_gap_filling)
        advanced_menu.addAction(act_gap_fill)

        act_gecko = QAction("GECKO — Enzyme-Constrained...", self)
        act_gecko.triggered.connect(self.run_gecko_analysis)
        advanced_menu.addAction(act_gecko)

        act_fseof = QAction("FSEOF — Flux Scanning...", self)
        act_fseof.triggered.connect(self.run_fseof)
        advanced_menu.addAction(act_fseof)

        act_pathway_enrich = QAction("Pathway Enrichment Analysis...", self)
        act_pathway_enrich.triggered.connect(self.run_pathway_enrichment)
        advanced_menu.addAction(act_pathway_enrich)

        act_escher = QAction("Escher Map Viewer...", self)
        act_escher.triggered.connect(self.show_escher_map)
        advanced_menu.addAction(act_escher)

        act_qt_map = QAction("Interactive Metabolic Map (native)...", self)
        act_qt_map.triggered.connect(self.show_interactive_metabolic_map)
        advanced_menu.addAction(act_qt_map)

        act_phpp = QAction("Phenotype Phase Plane (PhPP)...", self)
        act_phpp.triggered.connect(self.run_phenotype_phase_plane)
        advanced_menu.addAction(act_phpp)

        act_mcs = QAction("Minimal Cut Sets (MCS)...", self)
        act_mcs.triggered.connect(self.run_minimal_cut_sets)
        advanced_menu.addAction(act_mcs)

        act_met_dist = QAction("Metabolic Distance...", self)
        act_met_dist.triggered.connect(self.run_metabolic_distance)
        advanced_menu.addAction(act_met_dist)

        act_pareto = QAction("Pareto Frontier...", self)
        act_pareto.triggered.connect(self.run_pareto_frontier)
        advanced_menu.addAction(act_pareto)

        act_13c = QAction("¹³C-MFA Flux Constraints...", self)
        act_13c.triggered.connect(self.run_13c_mfa_constraints)
        advanced_menu.addAction(act_13c)

        # ---- Model comparison ----
        model_menu.addSeparator()

        act_model_compare = QAction("Compare Models...", self)
        act_model_compare.triggered.connect(self.run_model_comparison)
        model_menu.addAction(act_model_compare)

        act_validate_sbml = QAction("Validate SBML...", self)
        act_validate_sbml.triggered.connect(self.validate_sbml)
        model_menu.addAction(act_validate_sbml)

        # ---- Plugins ----
        help_menu.addSeparator()

        act_load_plugins = QAction("Load Plugins...", self)
        act_load_plugins.triggered.connect(self.load_plugins)
        help_menu.addAction(act_load_plugins)

    # ---- FBA result cache helpers ----

    def _model_state_hash(self, model: cobra.Model) -> str:
        """Return a lightweight hash representing the current model state.

        Captures objective, all reaction bounds, gene knockouts,
        overexpression factors, temporary UB overrides, reaction
        bound overrides, medium table state, analysis type,
        objective combo selection, solver, and analysis-specific
        parameters so that identical optimisation requests can be
        served from cache.
        """
        import hashlib
        parts: list[str] = []
        # objective
        parts.append(str(model.objective.expression))
        parts.append(model.objective.direction)
        # reaction bounds (sorted for determinism)
        for r in sorted(model.reactions, key=lambda r: r.id):
            parts.append(f"{r.id}:{r.lower_bound}:{r.upper_bound}")
        # knocked-out genes
        parts.append(",".join(sorted(self.knockout_genes)))
        # overexpression state
        parts.append("|".join(f"{k}={v}" for k, v in sorted(self.overexpression_reactions.items())))
        # temporary UB overrides
        parts.append("|".join(f"{k}={v}" for k, v in sorted(self.temporary_upper_bound_overrides.items())))
        # reaction bound overrides
        parts.append("|".join(f"{k}={v}" for k, v in sorted(self.reaction_bound_overrides.items())))
        # loopless flag
        parts.append(f"loopless:{self.loopless_chk.isChecked()}")
        # ── analysis-affecting inputs ──
        parts.append(f"analysis:{self.analysis_type.currentText()}")
        parts.append(f"solver:{self.solver_combo.currentText()}")
        parts.append(f"obj_combo:{self.objective_combo.currentData()}")
        # medium table state (exchange bounds edited by user)
        medium_parts: list[str] = []
        for row in range(self.medium_table.rowCount()):
            rid_item = self.medium_table.item(row, 0)
            lb_item = self.medium_table.item(row, 2)
            ub_item = self.medium_table.item(row, 3)
            if rid_item:
                medium_parts.append(
                    f"{rid_item.text()}:{lb_item.text() if lb_item else ''}:{ub_item.text() if ub_item else ''}"
                )
        parts.append("med:" + "|".join(medium_parts))
        # analysis-specific parameters (lazy-read from widgets)
        def _sv(attr):
            w = getattr(self, attr, None)
            return str(w.value()) if w is not None else ""
        def _st(attr):
            w = getattr(self, attr, None)
            return w.currentText() if w is not None else ""
        def _se(attr):
            w = getattr(self, attr, None)
            return w.text() if w is not None else ""
        parts.append(f"fva:{_sv('fva_fraction_spin')}:{_sv('fva_processes_spin')}")
        parts.append(f"samp:{_sv('sampling_size')}:{_st('sampler_type_combo')}")
        parts.append(f"rob:{_se('robustness_rxn')}:{_sv('robustness_min')}:{_sv('robustness_max')}:{_sv('robustness_steps')}")
        parts.append(f"env:{_se('envelope_product')}:{_sv('envelope_steps')}")
        blob = "|".join(parts).encode("utf-8")
        return hashlib.sha256(blob).hexdigest()

    def _cache_get(self, key: str) -> dict | None:
        """Look up a cached FBA result by its state-hash key."""
        return self._fba_cache.get(key)

    def _cache_put(self, key: str, result: dict) -> None:
        """Store an FBA result in the LRU cache, evicting the oldest entry if full."""
        if len(self._fba_cache) >= self._fba_cache_max:
            # evict oldest entry
            oldest = next(iter(self._fba_cache))
            del self._fba_cache[oldest]
        self._fba_cache[key] = result

    def _invalidate_fba_cache(self) -> None:
        """Clear the entire FBA result cache."""
        self._fba_cache.clear()

