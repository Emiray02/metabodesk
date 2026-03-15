"""Shared constants, TypedDicts, and stylesheet definitions.

Centralises values that are referenced by multiple modules:

- ``RECENT_MAX`` — maximum number of recent scenarios to track.
- ``BaselineResult`` / ``AnalysisResultDict`` — structured result types.
- ``DARK_STYLESHEET`` — Qt CSS for the dark theme.
- ``ARROW_PATTERNS`` — recognised arrow tokens for reaction equations.
"""

from typing import TypedDict, Any

RECENT_MAX = 5

# --------------- TypedDict definitions ----------------
class BaselineResult(TypedDict, total=False):
    status: str
    objective: float
    flux: dict[str, float]
    shadow_prices: dict[str, float]
    reduced_costs: dict[str, float]


class AnalysisResultDict(TypedDict, total=False):
    timestamp: str
    analysis_type: str
    baseline: BaselineResult
    wt_growth: float
    fva_result: dict[str, dict[str, float]]
    pfba: dict[str, Any]
    sgd_result: dict[str, float]
    dgd_result: dict[str, float]
    srd_result: dict[str, float]
    robustness_result: dict[str, list]
    envelope_result: dict[str, list]
    sampling_result: dict[str, dict]
    flux_rows: list[dict]
    rxn_id: str
    bound_type: str
    product_id: str
    sampler_type: str
    samples_df: Any


# --------------- Dark Theme Stylesheet ----------------
DARK_STYLESHEET = """
QMainWindow, QDialog { background-color: #1e1e1e; color: #d4d4d4; }
QMenuBar { background-color: #2d2d2d; color: #d4d4d4; }
QMenuBar::item:selected { background-color: #3e3e3e; }
QMenu { background-color: #2d2d2d; color: #d4d4d4; border: 1px solid #3e3e3e; }
QMenu::item:selected { background-color: #094771; }
QTabWidget::pane { border: 1px solid #3e3e3e; background: #1e1e1e; }
QTabBar::tab { background: #2d2d2d; color: #d4d4d4; padding: 6px 12px; border: 1px solid #3e3e3e; }
QTabBar::tab:selected { background: #1e1e1e; border-bottom-color: #1e1e1e; }
QTableWidget { background-color: #1e1e1e; color: #d4d4d4; gridline-color: #3e3e3e; selection-background-color: #094771; }
QHeaderView::section { background-color: #2d2d2d; color: #d4d4d4; border: 1px solid #3e3e3e; padding: 4px; }
QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QPlainTextEdit { background-color: #2d2d2d; color: #d4d4d4; border: 1px solid #3e3e3e; padding: 3px; }
QPushButton { background-color: #2d2d2d; color: #d4d4d4; border: 1px solid #3e3e3e; padding: 5px 12px; }
QPushButton:hover { background-color: #3e3e3e; }
QPushButton:pressed { background-color: #094771; }
QLabel { color: #d4d4d4; }
QListWidget { background-color: #1e1e1e; color: #d4d4d4; border: 1px solid #3e3e3e; }
QGroupBox { color: #d4d4d4; border: 1px solid #3e3e3e; margin-top: 6px; padding-top: 10px; }
QGroupBox::title { subcontrol-origin: margin; padding: 0 3px; }
QCheckBox, QRadioButton { color: #d4d4d4; }
QProgressBar { border: 1px solid #3e3e3e; background: #2d2d2d; text-align: center; color: #d4d4d4; }
QProgressBar::chunk { background-color: #094771; }
QScrollArea { background: #1e1e1e; border: none; }
QStatusBar { background: #2d2d2d; color: #d4d4d4; }
"""

# --------------- Equation parser constants ----------------
ARROW_PATTERNS = [
    ("<=>", True),
    ("<->", True),
    ("<-->", True),
    ("->", False),
    ("-->", False),
    ("=>", False),
]
