"""Utility functions shared across modules.

This module contains pure functions that do not depend on any Qt widget
or on the MainWindow instance.
"""

import sys
import re
import math
import logging
from pathlib import Path

from PySide6.QtCore import Qt, QStandardPaths
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import QMessageBox

import cobra

from metabodesk_core.constants import ARROW_PATTERNS

logger = logging.getLogger("MetaboDesk")

# ------------------------------------------------------------------ #
#  Paths                                                               #
# ------------------------------------------------------------------ #

def _recent_file_path() -> Path:
    """Return the Path to the recent-scenarios JSON file in AppData."""
    try:
        app_data = QStandardPaths.writableLocation(QStandardPaths.AppDataLocation)
        if app_data:
            p = Path(app_data)
            p.mkdir(parents=True, exist_ok=True)
            return p / "recent_scenarios.json"
    except Exception:
        pass
    return Path.home() / ".metabodesk_recent.json"


def _resource_path(rel_path: str) -> str:
    """Resolve a relative resource path, aware of PyInstaller frozen bundles."""
    try:
        if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
            return str(Path(sys._MEIPASS) / rel_path)
    except Exception:
        pass
    return str(Path(__file__).resolve().parent.parent / rel_path)


# ------------------------------------------------------------------ #
#  Small helpers                                                       #
# ------------------------------------------------------------------ #

def is_exchange_reaction(rxn: cobra.Reaction) -> bool:
    try:
        if rxn.id.startswith("EX_"):
            return True
        return len(rxn.metabolites) == 1
    except Exception:
        return False


def _as_str_list(v) -> list[str]:
    """Coerce a value (None, str, list, etc.) into a list of non-empty strings."""
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
    """Strip the identifiers.org URL prefix and return the bare identifier."""
    s = str(url_or_id).strip()
    if not s:
        return ""
    if "identifiers.org/" in s:
        s = s.split("identifiers.org/")[-1]
    if "/" in s:
        s = s.split("/")[-1]
    return s.strip()


def _safe_float(val, default: float = 0.0) -> float:
    """Convert *val* to float with NaN/Inf guard, returning *default* on failure."""
    try:
        f = float(val)
        if math.isnan(f) or math.isinf(f):
            return default
        return f
    except Exception:
        return default


# ------------------------------------------------------------------ #
#  Logo / message-box helpers                                          #
# ------------------------------------------------------------------ #

def _logo_pixmap() -> QPixmap | None:
    """Load and return a 64×64 QPixmap of the application logo, or None."""
    try:
        path = _resource_path("logo.ico")
        if Path(path).exists():
            pm = QPixmap(path)
            if not pm.isNull():
                return pm.scaled(64, 64, Qt.KeepAspectRatio, Qt.SmoothTransformation)
    except Exception:
        pass
    return None


def _show_msgbox(parent, title: str, text: str, icon: QMessageBox.Icon,
                 buttons: QMessageBox.StandardButtons = QMessageBox.Ok,
                 default_button: QMessageBox.StandardButton = QMessageBox.NoButton) -> int:
    """Show a branded QMessageBox with the application logo pixmap."""
    box = QMessageBox(parent)
    box.setWindowTitle(title)
    box.setText(text)
    pm = _logo_pixmap()
    if pm:
        box.setIconPixmap(pm)
    else:
        box.setIcon(icon)
    box.setStandardButtons(buttons)
    if default_button != QMessageBox.NoButton:
        box.setDefaultButton(default_button)
    return box.exec()


def _msg_info(parent, title, text, buttons=QMessageBox.Ok, default_button=QMessageBox.NoButton):
    """Shortcut for an Information-level branded message box."""
    return _show_msgbox(parent, title, text, QMessageBox.Information, buttons, default_button)


def _msg_warn(parent, title, text, buttons=QMessageBox.Ok, default_button=QMessageBox.NoButton):
    """Shortcut for a Warning-level branded message box."""
    return _show_msgbox(parent, title, text, QMessageBox.Warning, buttons, default_button)


def _msg_crit(parent, title, text, buttons=QMessageBox.Ok, default_button=QMessageBox.NoButton):
    """Shortcut for a Critical-level branded message box."""
    return _show_msgbox(parent, title, text, QMessageBox.Critical, buttons, default_button)


def _msg_question(parent, title, text, buttons=QMessageBox.Yes | QMessageBox.No, default_button=QMessageBox.NoButton):
    """Shortcut for a Question-level branded message box."""
    return _show_msgbox(parent, title, text, QMessageBox.Question, buttons, default_button)


def apply_custom_message_boxes():
    """Override QMessageBox static helpers with branded versions."""
    QMessageBox.information = _msg_info
    QMessageBox.warning = _msg_warn
    QMessageBox.critical = _msg_crit
    QMessageBox.question = _msg_question


# ------------------------------------------------------------------ #
#  GPR evaluator                                                       #
# ------------------------------------------------------------------ #

def evaluate_gpr_expression(gpr_rule: str, gene_values: dict[str, float]) -> float | None:
    """Recursively evaluate a GPR rule using AND=min, OR=max logic.

    Handles nested boolean expressions like ``(geneA and geneB) or geneC``.
    Returns None if none of the referenced genes have expression data.
    """
    if not gpr_rule or not gpr_rule.strip():
        return None

    rule = gpr_rule.strip()
    # Remove outer parentheses
    while rule.startswith('(') and rule.endswith(')'):
        depth = 0
        matched = True
        for i, ch in enumerate(rule):
            if ch == '(':
                depth += 1
            elif ch == ')':
                depth -= 1
            if depth == 0 and i < len(rule) - 1:
                matched = False
                break
        if matched:
            rule = rule[1:-1].strip()
        else:
            break

    # Split on top-level 'or' first (lowest precedence)
    depth = 0
    or_parts: list[str] = []
    current: list[str] = []
    tokens = rule.split()
    for tok in tokens:
        depth += tok.count('(') - tok.count(')')
        if tok.lower() == 'or' and depth == 0:
            or_parts.append(' '.join(current))
            current = []
        else:
            current.append(tok)
    or_parts.append(' '.join(current))

    if len(or_parts) > 1:
        vals = [evaluate_gpr_expression(p, gene_values) for p in or_parts]
        vals = [v for v in vals if v is not None]
        return max(vals) if vals else None

    # Split on top-level 'and' (higher precedence)
    depth = 0
    and_parts: list[str] = []
    current = []
    tokens = rule.split()
    for tok in tokens:
        depth += tok.count('(') - tok.count(')')
        if tok.lower() == 'and' and depth == 0:
            and_parts.append(' '.join(current))
            current = []
        else:
            current.append(tok)
    and_parts.append(' '.join(current))

    if len(and_parts) > 1:
        vals = [evaluate_gpr_expression(p, gene_values) for p in and_parts]
        vals = [v for v in vals if v is not None]
        return min(vals) if vals else None

    # Leaf node: single gene
    gene_id = rule.strip('() ')
    return gene_values.get(gene_id)


# ------------------------------------------------------------------ #
#  Annotation helpers                                                  #
# ------------------------------------------------------------------ #

def get_kegg_rxn_id(rxn: cobra.Reaction) -> str:
    ann = getattr(rxn, "annotation", {}) or {}

    keys = [
        "kegg.reaction", "kegg_reaction", "kegg", "kegg.rxn",
        "kegg_rxn", "KEGG Reaction", "KEGG_REACTION", "reaction.kegg",
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
        "ec-code", "ec", "ec_number", "ec-number",
        "ec-code(s)", "EC Number", "EC", "enzyme.ec",
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


# ------------------------------------------------------------------ #
#  Equation parser                                                     #
# ------------------------------------------------------------------ #

def _guess_compartment_from_met_id(met_id: str) -> str:
    """Infer compartment letter from a metabolite ID suffix (e.g. ``_e`` → ``e``)."""
    m = re.match(r"^(.+)_([A-Za-z])$", met_id.strip())
    if m:
        return m.group(2)
    return "c"


def _parse_coeff_and_met(term: str) -> tuple[float, str]:
    """Parse a term like ``'2 atp_c'`` into ``(2.0, 'atp_c')``."""
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
    """Return (stoichiometry, reversible).  Reactants negative, products positive."""
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
