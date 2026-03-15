"""External tools mixin for MetaboDesk.

Manages the lifecycle of bundled command-line tools:

- **memote**: Genome-scale metabolic model testing suite.  Runs
  ``memote report snapshot`` on the current model and opens the
  generated HTML report.
- **CarveMe**: Draft model reconstruction from a genome FASTA file.
  Runs ``carve`` with user-selected options and loads the resulting
  SBML model.

Also provides the community-model builder
(:meth:`_build_community_model`) and tool installation/repair from
bundled wheel packages.
"""

import sys
import os
import re
import math
import time
import logging
import socket
import urllib.request

import cobra
from pathlib import Path

from PySide6.QtCore import QUrl, QProcess, QTimer, QThread, Signal, QObject
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QVBoxLayout, QLabel, QFileDialog, QMessageBox, QListWidget,
    QLineEdit, QSpinBox, QCheckBox, QComboBox, QDialog,
    QDialogButtonBox, QFormLayout, QInputDialog,
)

logger = logging.getLogger("MetaboDesk")


class _ToolCheckWorker(QThread):
    """Non-blocking worker that checks tool modules in a background thread."""
    finished = Signal(list)           # emits list of missing module names

    def __init__(self, py_path: str, modules: list[str], parent=None):
        super().__init__(parent)
        self._py = py_path
        self._modules = modules

    def run(self):
        import subprocess
        # On Windows, CREATE_NO_WINDOW (0x08000000) prevents a
        # transient console from flashing while checking modules.
        import sys as _sys
        _cflags = 0x08000000 if _sys.platform.startswith("win") else 0
        missing = []
        for mod in self._modules:
            try:
                result = subprocess.run(
                    [self._py, "-c", f"import {mod}"],
                    capture_output=True, timeout=5,
                    creationflags=_cflags,
                )
                if result.returncode != 0:
                    missing.append(mod)
            except Exception:
                missing.append(mod)
        self.finished.emit(missing)


class _UpdateCheckWorker(QThread):
    """Non-blocking worker for GitHub update check."""
    result_ready = Signal(dict)       # emits {"tag": ..., "url": ..., "body": ...} or {"error": ...}

    def __init__(self, repo: str, parent=None):
        super().__init__(parent)
        self._repo = repo

    def run(self):
        import json as _json
        try:
            url = f"https://api.github.com/repos/{self._repo}/releases/latest"
            req = urllib.request.Request(url, headers={
                "User-Agent": "MetaboDesk-UpdateChecker/1.0",
                "Accept": "application/vnd.github.v3+json",
            })
            with urllib.request.urlopen(req, timeout=8) as resp:
                data = _json.loads(resp.read().decode("utf-8"))
            self.result_ready.emit({
                "tag": str(data.get("tag_name", "")).lstrip("vV").strip(),
                "url": data.get("html_url", f"https://github.com/{self._repo}/releases"),
                "body": (data.get("body") or "")[:500],
            })
        except Exception as e:
            self.result_ready.emit({"error": str(e)})


class ToolsMixin:
    """Mixin providing tools functionality."""

    def cancel_memote(self):
        if self._memote_proc is None:
            return
        try:
            self._memote_cancel_requested = True
            self.tools_log.appendPlainText("\nCancel requested...\n")
            self._force_kill_process(self._memote_proc)
            self._memote_proc = None
            self.set_busy(False, "Ready.")
            self._stop_memote_timer()
        except Exception:
            pass

    def _start_elapsed_timer(self, kind: str):
        """Start a QTimer that updates the elapsed-time label every second."""
        start = time.monotonic()
        timer = QTimer(self)
        timer.setInterval(1000)

        def tick():
            elapsed = int(time.monotonic() - start)
            mm = elapsed // 60
            ss = elapsed % 60
            txt = f"Elapsed: {mm:02d}:{ss:02d}"
            if kind == "memote":
                self.memote_time_lbl.setText(txt)
            else:
                self.carveme_time_lbl.setText(txt)

        timer.timeout.connect(tick)
        timer.start()
        if kind == "memote":
            self._memote_timer = timer
        else:
            self._carveme_timer = timer

    def _stop_memote_timer(self):
        """Stop the memote elapsed-time timer and reset label."""
        try:
            if getattr(self, "_memote_timer", None):
                self._memote_timer.stop()
        except Exception:
            pass
        self.memote_time_lbl.setText("Elapsed: 00:00")

    def _stop_carveme_timer(self):
        """Stop the CarveMe elapsed-time timer and reset label."""
        try:
            if getattr(self, "_carveme_timer", None):
                self._carveme_timer.stop()
        except Exception:
            pass
        self.carveme_time_lbl.setText("Elapsed: 00:00")

    def _get_objective_reaction_id(self, model: cobra.Model | None = None) -> str | None:
        """Return the objective reaction ID for the given model (or self.base_model)."""
        if model is None:
            model = self.base_model
        if model is None or model.objective is None:
            return None
        try:
            obj = model.objective
            if hasattr(obj, "expression"):
                coeffs = obj.expression.as_coefficients_dict()
                for var in coeffs.keys():
                    name = getattr(var, "name", "") or ""
                    rxn_id = name
                    if rxn_id.startswith("R_") or rxn_id.startswith("F_"):
                        rxn_id = rxn_id[2:]
                    if rxn_id in model.reactions:
                        return rxn_id
            if hasattr(obj, "id") and obj.id in model.reactions:
                return obj.id
        except Exception:
            pass
        return None

    # Keep backward-compatible alias used by _build_community_model

    def _get_objective_reaction_id_for_model(self, model: cobra.Model) -> str | None:
        """Backward-compatible alias: delegate to ``_get_objective_reaction_id``."""
        return self._get_objective_reaction_id(model)

    def _is_extracellular(self, met: cobra.Metabolite) -> bool:
        """Check if a metabolite belongs to the extracellular compartment."""
        try:
            if met.compartment and met.compartment.lower() in ("e", "extracellular"):
                return True
            mid = met.id.lower()
            return mid.endswith("_e") or mid.endswith("[e]") or mid.endswith("(e)")
        except Exception:
            return False

    def _prefix_gpr(self, gpr: str, prefix: str) -> str:
        """Prefix all gene identifiers in a GPR rule string.

        Uses word-boundary regex (``\\b``) to avoid splitting gene IDs
        that contain 'and' or 'or' as substrings (e.g. ``b_orA1``).
        """
        if not gpr:
            return gpr
        tokens = re.split(r'(\s+|\(|\)|\band\b|\bor\b)', gpr, flags=re.IGNORECASE)
        out = []
        for t in tokens:
            if t is None:
                continue
            if t.strip() == "":
                out.append(t)
                continue
            if t.lower() in ("and", "or") or t in ("(", ")"):
                out.append(t)
                continue
            out.append(prefix + t)
        return "".join(out)

    def _build_community_model(self, models: list[cobra.Model], labels: list[str],
                               share_extracellular: bool = True,
                               objective_mode: str = "sum",
                               abundance: dict[str, float] | None = None):
        """Build a community model from multiple organism models.

        Parameters
        ----------
        objective_mode : str
            ``"sum"``  — maximise the sum of member objectives (default).
            ``"steadycom"`` — SteadyCom-like equal-growth constraint:
            all species share a common growth rate *μ* which is maximised.
        abundance : dict[str, float] or None
            Relative abundance fractions keyed by species label.
            Values should sum to 1.0 (e.g. ``{"sp1": 0.7, "sp2": 0.3}``).
            When supplied, each species' reaction bounds are scaled by
            its fractional abundance so that a species with 30% abundance
            can use at most 30% of the shared pool capacity.
            ``None`` means equal abundance (no scaling).

        Returns
        -------
        tuple : (community, obj_rxn_map, shared_mets, met_conflicts)
        """
        community = cobra.Model("community")
        shared_mets: dict[str, cobra.Metabolite] = {}
        obj_rxn_map: dict[str, str] = {}
        # Track metabolite conflicts (formula/charge mismatch across organisms)
        met_conflicts: list[str] = []

        # ── Phase 1: Collect all metabolites and reactions per model ──
        # Batch additions for better scalability with large models
        all_new_metabolites: list[cobra.Metabolite] = []
        all_new_reactions: list[cobra.Reaction] = []

        for model, label in zip(models, labels):
            prefix = f"{label}__"
            met_map: dict[str, cobra.Metabolite] = {}

            for met in model.metabolites:
                if share_extracellular and self._is_extracellular(met):
                    if met.id not in shared_mets:
                        new_met = cobra.Metabolite(
                            met.id,
                            name=met.name,
                            compartment=met.compartment,
                            formula=met.formula,
                            charge=met.charge,
                        )
                        all_new_metabolites.append(new_met)
                        shared_mets[met.id] = new_met
                    else:
                        # Check for formula/charge conflicts with existing shared metabolite
                        existing = shared_mets[met.id]
                        if met.formula and existing.formula and met.formula != existing.formula:
                            met_conflicts.append(
                                f"  {met.id}: formula mismatch — "
                                f"'{existing.formula}' vs '{met.formula}' ({label})")
                        if met.charge is not None and existing.charge is not None and met.charge != existing.charge:
                            met_conflicts.append(
                                f"  {met.id}: charge mismatch — "
                                f"{existing.charge} vs {met.charge} ({label})")
                    met_map[met.id] = shared_mets[met.id]
                else:
                    new_id = prefix + met.id
                    new_met = cobra.Metabolite(
                        new_id,
                        name=f"{met.name} ({label})" if met.name else new_id,
                        compartment=met.compartment,
                        formula=met.formula,
                        charge=met.charge,
                    )
                    all_new_metabolites.append(new_met)
                    met_map[met.id] = new_met

            for rxn in model.reactions:
                new_rxn = cobra.Reaction(prefix + rxn.id)
                new_rxn.name = f"{rxn.name} ({label})" if rxn.name else new_rxn.id
                new_rxn.lower_bound = rxn.lower_bound
                new_rxn.upper_bound = rxn.upper_bound
                # Apply abundance scaling: species with 30% abundance can
                # only use at most 30% of its original capacity.
                if abundance and label in abundance:
                    frac = abundance[label]
                    if new_rxn.lower_bound != 0:
                        new_rxn.lower_bound *= frac
                    if new_rxn.upper_bound != 0:
                        new_rxn.upper_bound *= frac
                new_rxn.subsystem = rxn.subsystem
                met_coeffs = {met_map[m.id]: coeff for m, coeff in rxn.metabolites.items()}
                new_rxn.add_metabolites(met_coeffs)
                try:
                    new_rxn.gene_reaction_rule = self._prefix_gpr(getattr(rxn, "gene_reaction_rule", "") or "", prefix)
                except Exception as e:
                    logger.warning(f"GPR not set for {new_rxn.id}: {e}")
                all_new_reactions.append(new_rxn)

            obj_id = self._get_objective_reaction_id_for_model(model)
            if obj_id:
                obj_rxn_map[label] = prefix + obj_id

        # ── Phase 2: Batch-add all metabolites and reactions at once ──
        if all_new_metabolites:
            community.add_metabolites(all_new_metabolites)
        if all_new_reactions:
            community.add_reactions(all_new_reactions)

        # ── Phase 3: Set community objective ──
        if objective_mode == "steadycom" and len(obj_rxn_map) >= 2:
            # SteadyCom-like: all species must grow at the same rate μ.
            # Create a community growth variable and constrain each
            # species' biomass reaction to equal it.
            mu_met = cobra.Metabolite("community_mu", name="Community growth rate",
                                       compartment="c")
            community.add_metabolites([mu_met])

            # Each species' biomass produces community_mu; one drain consumes it
            for label, rxn_id in obj_rxn_map.items():
                if rxn_id in community.reactions:
                    community.reactions.get_by_id(rxn_id).add_metabolites({mu_met: 1.0})

            mu_drain = cobra.Reaction("community_mu_drain")
            mu_drain.name = "Community growth drain"
            mu_drain.lower_bound = 0
            mu_drain.upper_bound = 1000
            mu_drain.add_metabolites({mu_met: -float(len(obj_rxn_map))})
            community.add_reactions([mu_drain])
            community.objective = "community_mu_drain"
        else:
            # Default: sum of member objectives
            for label, rxn_id in obj_rxn_map.items():
                if rxn_id in community.reactions:
                    community.reactions.get_by_id(rxn_id).objective_coefficient = 1.0

        return community, obj_rxn_map, shared_mets, met_conflicts

    def _show_issues_dialog(self, title: str, issues: list[str]):
        """Show a filterable list dialog of model issues/warnings."""
        dialog = QDialog(self)
        dialog.setWindowTitle(title)
        dialog.resize(700, 500)

        layout = QVBoxLayout(dialog)
        search_input = QLineEdit()
        search_input.setPlaceholderText("Filter...")
        layout.addWidget(search_input)

        list_widget = QListWidget()
        layout.addWidget(list_widget, stretch=1)

        def populate(filter_text: str = ""):
            list_widget.clear()
            f = filter_text.strip().lower()
            for item in issues:
                if not f or f in item.lower():
                    list_widget.addItem(item)

        populate()
        search_input.textChanged.connect(populate)

        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(dialog.reject)
        btns.accepted.connect(dialog.accept)
        layout.addWidget(btns)

        dialog.exec()

    def _force_kill_process(self, proc: QProcess):
        """Forcefully terminate a QProcess with taskkill fallback on Windows."""
        try:
            proc.kill()
            if proc.waitForFinished(2000):
                return
        except Exception:
            pass
        try:
            if sys.platform.startswith("win"):
                pid = proc.processId()
                if pid:
                    QProcess.startDetached("taskkill", ["/F", "/T", "/PID", str(pid)])
            else:
                proc.kill()
        except Exception:
            pass

    def cancel_carveme(self):
        if getattr(self, "_carveme_proc", None) is None:
            return
        try:
            self.tools_log.appendPlainText("\nCancel requested (CarveMe)...\n")
            proc = self._carveme_proc
            self._carveme_proc = None
            self._carveme_last_output = ""
            try:
                proc.terminate()
            except Exception:
                pass
            QTimer.singleShot(2000, lambda: proc.kill() if proc.state() != QProcess.NotRunning else None)
            self.set_busy(False, "Ready.")
            self._stop_carveme_timer()
        except Exception:
            pass

    def _start_memote_run(self, sbml_path: Path, report_path: Path, extra_args: list[str] | None = None, title: str = "memote", timeout_seconds: int | None = 900):
        """Launch a memote report subprocess with timeout handling."""
        args = ["report", "snapshot", str(sbml_path), "--filename", str(report_path)]
        if extra_args:
            args += extra_args

        cmd = self._console_script_command("memote", args)
        if not cmd:
            self._show_error("Tools missing", "Tools not found. Use Tools → Repair tools to reinstall.")
            return
        program, proc_args = cmd

        proc = self._run_tool_process(program, proc_args, title)
        try:
            proc.setWorkingDirectory(str(report_path.parent))
        except Exception:
            pass

        self._memote_proc = proc
        self._memote_cancel_requested = False
        self._memote_last_output = time.monotonic()
        self._memote_current_report_path = report_path
        self._memote_current_sbml_path = sbml_path
        self._memote_current_extra_args = extra_args or []
        self.set_busy(True, "Running memote...")
        self._start_elapsed_timer("memote")

        try:
            if getattr(self, "_memote_watchdog", None):
                self._memote_watchdog.stop()
        except Exception:
            pass

        if timeout_seconds is not None:
            self._memote_watchdog = QTimer(self)
            self._memote_watchdog.setInterval(2000)

            def check_memote_hang():
                if self._memote_proc is None:
                    self._memote_watchdog.stop()
                    return
                if time.monotonic() - self._memote_last_output > timeout_seconds:
                    self.tools_log.appendPlainText("\nmemote appears to be stuck. Stopping process...\n")
                    try:
                        self._force_kill_process(self._memote_proc)
                    except Exception:
                        pass
                    self._memote_proc = None
                    self.set_busy(False, "Ready.")
                    self._memote_watchdog.stop()
                    self._stop_memote_timer()

                    if not getattr(self, "_memote_fast_fallback_used", False) and not extra_args:
                        self._memote_fast_fallback_used = True
                        fast_path = report_path.with_name(report_path.stem + "_fast" + report_path.suffix)
                        self.tools_log.appendPlainText("\nRetrying memote with heavy tests skipped...\n")
                        self._start_memote_run(
                            sbml_path,
                            fast_path,
                            extra_args=["--skip", "test_consistency", "--skip", "test_biomass"],
                            title="memote (fast)",
                            timeout_seconds=timeout_seconds,
                        )

            self._memote_watchdog.timeout.connect(check_memote_hang)
            self._memote_watchdog.start()

        def on_finished(code, status):
            self._append_proc_output_to_log(proc)
            self.tools_log.appendPlainText(f"\nmemote finished with code={code}, status={status}\n")
            self._memote_proc = None
            self.set_busy(False, "Ready.")
            self._stop_memote_timer()
            try:
                self._memote_watchdog.stop()
            except Exception:
                pass

            if report_path.exists():
                QDesktopServices.openUrl(QUrl.fromLocalFile(str(report_path)))
                QMessageBox.information(self, "Memote finished", f"Report saved:\n{report_path}")
            else:
                if self._memote_cancel_requested:
                    QMessageBox.information(self, "Memote cancelled", "memote was cancelled.")
                else:
                    if not getattr(self, "_memote_fast_fallback_used", False) and not extra_args:
                        self._memote_fast_fallback_used = True
                        fast_path = report_path.with_name(report_path.stem + "_fast" + report_path.suffix)
                        self.tools_log.appendPlainText("\nmemote crashed. Retrying with heavy tests skipped...\n")
                        self._start_memote_run(
                            sbml_path,
                            fast_path,
                            extra_args=["--skip", "test_consistency", "--skip", "test_biomass"],
                            title="memote (fast)",
                        )
                        return
                    self._show_error("Memote failed", "memote failed and report was not created. See log in Tools tab.")

        proc.finished.connect(on_finished)
        proc.start()

    def _tool_python_path(self) -> str | None:
        """Return bundled Python path, or current Python in dev mode."""
        py = self._bundled_python_path()
        if py:
            return py
        # Dev mode fallback
        if not getattr(sys, "frozen", False):
            exe = str(getattr(sys, "executable", "") or "")
            if exe.lower().endswith("python.exe") or exe.lower().endswith("python"):
                return exe
        return None

    def _app_base_dir(self) -> Path:
        """Return the application base directory (frozen exe dir or source dir)."""
        # Launcher sets METABODESK_BASE to the install directory
        env_base = os.environ.get("METABODESK_BASE")
        if env_base:
            return Path(env_base)
        if getattr(sys, "frozen", False):
            return Path(sys.executable).resolve().parent
        # __file__ is in metabodesk_core/ -> go one level up
        return Path(__file__).resolve().parent.parent

    def _bundled_python_path(self) -> str | None:
        """Return the path to the bundled Python interpreter, if present.

        On Windows, *pythonw.exe* is preferred over *python.exe* because
        it has the GUI subsystem flag — launching it never creates a
        transient console window.
        """
        base = self._runtime_dir()
        cand = [
            base / "python" / "pythonw.exe",   # GUI — no console
            base / "python" / "python.exe",    # fallback
            base / "python" / "bin" / "python",
        ]
        for p in cand:
            if p.exists():
                return str(p)
        return None

    def _download_file(self, url: str, dest: Path, progress_cb=None, cancel_event=None, max_seconds: int = 600, stall_seconds: int = 20):
        """Download a file from *url* with progress callback and timeout/stall guards."""
        dest.parent.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(url, headers={"User-Agent": "MetaboDesk/1.0"})
        start_time = time.monotonic()
        last_progress = start_time
        with urllib.request.urlopen(req, timeout=10) as r, open(dest, "wb") as f:
            total = r.headers.get("Content-Length")
            total = int(total) if total and total.isdigit() else None
            downloaded = 0
            while True:
                if cancel_event and cancel_event.is_set():
                    raise TimeoutError("Download cancelled")
                if time.monotonic() - start_time > max_seconds:
                    raise TimeoutError("Download timed out")
                if time.monotonic() - last_progress > stall_seconds:
                    raise TimeoutError("Download stalled")
                try:
                    chunk = r.read(1024 * 1024)
                except socket.timeout:
                    continue
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                last_progress = time.monotonic()
                if progress_cb:
                    progress_cb(downloaded, total)

    def _get_remote_size(self, url: str) -> int | None:
        """Perform a HEAD request to get Content-Length of a remote file."""
        try:
            req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "MetaboDesk/1.0"})
            with urllib.request.urlopen(req, timeout=10) as r:
                total = r.headers.get("Content-Length")
                return int(total) if total and total.isdigit() else None
        except Exception:
            return None

    def _check_tool_module(self, module_name: str) -> bool:
        """Check if a Python module can be imported in the tools interpreter."""
        py = self._tool_python_path()
        if not py:
            return False
        proc = QProcess(self)
        proc.setProgram(py)
        proc.setArguments(["-c", f"import {module_name}"])
        proc.setProcessChannelMode(QProcess.MergedChannels)
        proc.start()
        proc.waitForFinished(3000)
        return proc.exitCode() == 0

    def _has_module(self, module_name: str) -> bool:
        """Check if a module spec exists in the tools interpreter."""
        py = self._tool_python_path()
        if not py:
            return False
        code = (
            "import importlib.util, sys\n"
            f"spec = importlib.util.find_spec('{module_name}')\n"
            "sys.exit(0 if spec else 1)\n"
        )
        proc = QProcess(self)
        proc.setProgram(py)
        proc.setArguments(["-c", code])
        proc.setProcessChannelMode(QProcess.MergedChannels)
        proc.start()
        proc.waitForFinished(3000)
        return proc.exitCode() == 0

    def _has_console_script(self, script_name: str) -> bool:
        """Check if a named console_scripts entry point exists."""
        py = self._tool_python_path()
        if not py:
            return False
        code = (
            "from importlib.metadata import entry_points\n"
            "eps = entry_points().select(group='console_scripts')\n"
            f"found = any(e.name=='{script_name}' for e in eps)\n"
            "raise SystemExit(0 if found else 1)\n"
        )
        proc = QProcess(self)
        proc.setProgram(py)
        proc.setArguments(["-c", code])
        proc.setProcessChannelMode(QProcess.MergedChannels)
        proc.start()
        proc.waitForFinished(3000)
        return proc.exitCode() == 0

    def _console_script_command(self, script_name: str, argv: list[str]) -> tuple[str, list[str]] | None:
        """Build (program, args) tuple for a console script entry point.

        On Windows we deliberately avoid launching ``.exe`` console scripts
        directly — they have the *console* PE subsystem and would flash a
        command-prompt window.  Instead we prefer ``pythonw -m <pkg>`` or
        the entry-point trampoline via ``pythonw -c ...``.
        """
        py = self._tool_python_path()
        if not py:
            return None

        # ── Prefer module execution (no console window) ──────────
        _mod_map = {"memote": "memote", "carve": "carveme.cli.carve", "carveme": "carveme"}
        mod_name = _mod_map.get(script_name)
        if mod_name and self._has_module(mod_name):
            return py, ["-m", mod_name] + argv

        # ── Non-Windows: try the script executable directly ──────
        py_path = Path(py)
        scripts_dir = py_path.parent / ("Scripts" if sys.platform.startswith("win") else "bin")
        if not sys.platform.startswith("win"):
            exe = scripts_dir / script_name
            if exe.exists():
                return str(exe), argv

        # ── Fallback: entry-point trampoline via interpreter ─────
        code = (
            "from importlib.metadata import entry_points\n"
            "eps = entry_points().select(group='console_scripts')\n"
            f"eps = [e for e in eps if e.name=='{script_name}']\n"
            "if not eps: raise SystemExit(127)\n"
            "func = eps[0].load()\n"
            f"import sys; sys.argv={[script_name] + argv}\n"
            "raise SystemExit(func())\n"
        )
        return py, ["-c", code]

    def _ensure_carveme_shim(self, py_path_str: str):
        """Create a shell shim script for CarveMe if missing."""
        py_path = Path(py_path_str)
        scripts_dir = py_path.parent / ("Scripts" if sys.platform.startswith("win") else "bin")
        scripts_dir.mkdir(parents=True, exist_ok=True)
        if sys.platform.startswith("win"):
            shim = scripts_dir / "carveme.bat"
            if not shim.exists():
                shim.write_text(
                    "@echo off\r\n"
                    "\"%~dp0\\python.exe\" -m carveme %*\r\n",
                    encoding="utf-8",
                )
        else:
            shim = scripts_dir / "carveme"
            if not shim.exists():
                shim.write_text(
                    f"#!/bin/sh\n\"{py_path_str}\" -m carveme \"$@\"\n",
                    encoding="utf-8",
                )
                try:
                    os.chmod(shim, 0o755)
                except Exception:
                    pass

    def _carveme_command(self, genome_path: str | None, out_path: str, extra_args: list[str] | None = None) -> tuple[str, list[str]] | None:
        py = self._tool_python_path()
        if not py:
            return None
        extra_args = extra_args or []
        py_path = Path(py)
        scripts_dir = py_path.parent / ("Scripts" if sys.platform.startswith("win") else "bin")

        # ── Prefer module execution (pythonw, no console) ────────
        if sys.platform.startswith("win"):
            # Try module execution first to avoid console-subsystem .exe
            if self._has_module("carveme.cli.carve"):
                self._carveme_last_mode = "carve"
                args = ["-m", "carveme.cli.carve"] + ([str(genome_path)] if genome_path else []) + ["-o", str(out_path)] + extra_args
                return py, args
            if self._has_module("carveme.__main__"):
                self._carveme_last_mode = "carveme"
                args = ["-m", "carveme", "draft"]
                if genome_path:
                    args += ["-g", str(genome_path)]
                args += ["-o", str(out_path)] + extra_args
                return py, args

        # ── Non-Windows: try console script executables ──────────
        carve_exe = scripts_dir / ("carve.exe" if sys.platform.startswith("win") else "carve")
        if carve_exe.exists() and not sys.platform.startswith("win"):
            self._carveme_last_mode = "carve"
            args = ([str(genome_path)] if genome_path else []) + ["-o", str(out_path)] + extra_args
            return str(carve_exe), args

        carveme_exe = scripts_dir / ("carveme.exe" if sys.platform.startswith("win") else "carveme")
        if carveme_exe.exists() and not sys.platform.startswith("win"):
            self._carveme_last_mode = "carveme"
            args = ["draft"]
            if genome_path:
                args += ["-g", str(genome_path)]
            args += ["-o", str(out_path)] + extra_args
            return str(carveme_exe), args

        if self._has_console_script("carve"):
            self._carveme_last_mode = "carve"
            args = ([str(genome_path)] if genome_path else []) + ["-o", str(out_path)] + extra_args
            return self._console_script_command("carve", args)

        if self._has_console_script("carveme"):
            self._carveme_last_mode = "carveme"
            args = ["draft"]
            if genome_path:
                args += ["-g", str(genome_path)]
            args += ["-o", str(out_path)] + extra_args
            return self._console_script_command("carveme", args)

        # Try module execution (carveme.__main__) if available
        if self._has_module("carveme.__main__"):
            self._carveme_last_mode = "carveme"
            args = ["-m", "carveme", "draft"]
            if genome_path:
                args += ["-g", str(genome_path)]
            args += ["-o", str(out_path)] + extra_args
            return py, args

        # Try module execution for carve
        if self._has_module("carveme.cli.carve"):
            self._carveme_last_mode = "carve"
            args = ["-m", "carveme.cli.carve"] + ([str(genome_path)] if genome_path else []) + ["-o", str(out_path)] + extra_args
            return py, args
        if self._has_module("carveme.carve"):
            self._carveme_last_mode = "carve"
            args = ["-m", "carveme.carve"] + ([str(genome_path)] if genome_path else []) + ["-o", str(out_path)] + extra_args
            return py, args

        self._carveme_last_mode = "carve"
        argv = ["carve"] + ([str(genome_path)] if genome_path else []) + ["-o", str(out_path)] + extra_args
        code = (
            "import sys, importlib, runpy, traceback, pkgutil\n"
            f"sys.argv={argv!r}\n"
            "errors=[]\n"
            "mods = ['carveme.cli.carve', 'carveme.carve', 'carveme', 'carveme.__main__', 'carveme.cli', 'carveme.cli.main', 'carveme.main', 'carveme.carveme']\n"
            "for m in mods:\n"
            "    try:\n"
            "        runpy.run_module(m, run_name='__main__')\n"
            "        raise SystemExit(0)\n"
            "    except Exception as e:\n"
            "        errors.append((m, str(e)))\n"
            "for m in mods:\n"
            "    try:\n"
            "        mod = importlib.import_module(m)\n"
            "    except Exception as e:\n"
            "        errors.append((m, str(e)))\n"
            "        continue\n"
            "    for name in ('main','cli','app','draft','run'):\n"
            "        fn = getattr(mod, name, None)\n"
            "        if callable(fn):\n"
            "            raise SystemExit(fn())\n"
            "try:\n"
            "    import carveme as _carveme_pkg\n"
            "    for finder, name, ispkg in pkgutil.walk_packages(_carveme_pkg.__path__, _carveme_pkg.__name__ + '.'):\n"
            "        if any(x in name for x in ('cli', 'main')):\n"
            "            try:\n"
            "                mod = importlib.import_module(name)\n"
            "            except Exception as e:\n"
            "                errors.append((name, str(e)))\n"
            "                continue\n"
            "            for fn_name in ('main','cli','app','draft','run'):\n"
            "                fn = getattr(mod, fn_name, None)\n"
            "                if callable(fn):\n"
            "                    raise SystemExit(fn())\n"
            "except Exception as e:\n"
            "    errors.append(('carveme package scan', str(e)))\n"
            "from importlib.metadata import entry_points\n"
            "eps = entry_points().select(group='console_scripts')\n"
            "eps = [e for e in eps if e.name=='carveme']\n"
            "if eps:\n"
            "    raise SystemExit(eps[0].load()())\n"
            "print('carveme entrypoint not found')\n"
            "print('errors:')\n"
            "for m, e in errors:\n"
            "    print(' -', m, ':', e)\n"
            "raise SystemExit(1)\n"
        )
        return py, ["-c", code]

    def _run_tool_process(self, program: str, args: list[str], title: str):
        """Create a QProcess with environment, logging, and return it.

        On Windows, if *program* is a console-subsystem executable
        (e.g. ``python.exe``, ``memote.exe``) we silently swap it for
        the windowless ``pythonw.exe`` when it exists side-by-side, so
        that QProcess never flashes a console window.
        """
        # ── Swap console python for pythonw (no console) ─────────
        if sys.platform.startswith("win"):
            prog_lower = os.path.basename(program).lower()
            if prog_lower == "python.exe":
                alt = os.path.join(os.path.dirname(program), "pythonw.exe")
                if os.path.isfile(alt):
                    program = alt
        proc = QProcess(self)
        proc.setProgram(program)
        proc.setArguments(args)
        try:
            env = os.environ.copy()
            tools_bin = self._tools_bin_dir()
            if tools_bin.exists():
                env["PATH"] = str(tools_bin) + os.pathsep + env.get("PATH", "")
            proc.setEnvironment([f"{k}={v}" for k, v in env.items()])
        except Exception:
            pass
        proc.setProcessChannelMode(QProcess.MergedChannels)
        proc.readyRead.connect(lambda: self._append_proc_output_to_log(proc))
        self.tools_log.appendPlainText(f"\n[{title}] {' '.join([program] + args)}\n")
        return proc

    def _reinstall_carveme(self):
        """Reinstall CarveMe via pip (offline wheels or network)."""
        py = self._tool_python_path()
        if not py:
            self._show_error("Tools missing", "Bundled Python not found. Please reinstall the app.")
            return

        wheels_dir = self._wheels_dir()
        if wheels_dir.exists():
            args = [
                "-m", "pip", "install",
                "--no-index", "--find-links", str(wheels_dir),
                "--force-reinstall", "carveme",
            ]
        else:
            args = [
                "-m", "pip", "install", "--force-reinstall", "--no-cache-dir",
                "carveme",
            ]
        proc = self._run_tool_process(py, args, "pip reinstall carveme")
        self.set_busy(True, "Reinstalling CarveMe...")

        def on_finished(code, status):
            self._append_proc_output_to_log(proc)
            self.set_busy(False, "Ready.")
            if code == 0:
                QMessageBox.information(self, "CarveMe", "CarveMe reinstalled. Try again.")
            else:
                self._show_error("CarveMe", f"Reinstall failed (code {code}).")
            self.tools_check_status()

        proc.finished.connect(on_finished)
        proc.start()

    def _append_proc_output_to_log(self, proc: QProcess):
        """Read QProcess output and append it to the tools log widget."""
        data = proc.readAll().data().decode(errors="replace")
        if not data:
            return
        if proc is getattr(self, "_memote_proc", None):
            self._memote_last_output = time.monotonic()
        if proc is getattr(self, "_carveme_proc", None):
            self._carveme_last_output = (getattr(self, "_carveme_last_output", "") + data)[-8000:]
        for line in data.splitlines():
            self.tools_log.appendPlainText(line)

    def _carveme_preflight(self) -> bool:
        """Run a preflight check ensuring CarveMe imports correctly."""
        py = self._tool_python_path()
        if not py:
            return False
        tools_dir = self._tools_bin_dir()
        code = (
            "import sys, traceback, shutil, os\n"
            "try:\n"
            "    import carveme\n"
            "    print('carveme:', getattr(carveme, '__version__', 'unknown'))\n"
            "except Exception as e:\n"
            "    print('carveme import error:', e)\n"
            "    traceback.print_exc()\n"
            "    sys.exit(2)\n"
            f"tools_dir = r'{tools_dir}'\n"
            "if os.path.isdir(tools_dir):\n"
            "    try:\n"
            "        print('tools_dir:', tools_dir)\n"
            "        print('tools_dir_files:', ','.join(os.listdir(tools_dir)))\n"
            "    except Exception:\n"
            "        pass\n"
            "missing=[]\n"
            "for b in ('diamond','prodigal'):\n"
            "    found = shutil.which(b) is not None\n"
            "    if not found and os.name == 'nt':\n"
            "        for ext in (b+'.exe', b+'.windows.exe'):\n"
            "            if os.path.exists(os.path.join(tools_dir, ext)):\n"
            "                found = True\n"
            "                break\n"
            "    if not found:\n"
            "        missing.append(b)\n"
            "if missing:\n"
            "    print('missing binaries:', ','.join(missing))\n"
            "    sys.exit(3)\n"
            "sys.exit(0)\n"
        )
        proc = QProcess(self)
        proc.setProgram(py)
        proc.setArguments(["-c", code])
        try:
            env = os.environ.copy()
            tools_bin = self._tools_bin_dir()
            if tools_bin.exists():
                env["PATH"] = str(tools_bin) + os.pathsep + env.get("PATH", "")
            proc.setEnvironment([f"{k}={v}" for k, v in env.items()])
        except Exception:
            pass
        proc.setProcessChannelMode(QProcess.MergedChannels)
        proc.start()
        proc.waitForFinished(5000)
        self._append_proc_output_to_log(proc)
        if proc.exitCode() == 0:
            return True
        if proc.exitCode() == 3:
            QMessageBox.warning(
                self,
                "CarveMe dependencies",
                "CarveMe requires external binaries (diamond, prodigal).\n"
                "Please bundle them under runtime/tools or add them to PATH.",
            )
            return False
        self._show_error("CarveMe import failed", "CarveMe failed to import. See Tools log for details.")
        return False

    def run_memote_report_current_model(self):
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load an SBML model first.")
            return

        py = self._tool_python_path()
        if not py:
            self._show_error("Tools missing", "Bundled Python not found. Please reinstall the app.")
            return

        if not self._check_tool_module("memote"):
            self._show_error("Tools missing", "memote is not installed. Use Tools → Repair tools to reinstall.")
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

        if report_path.exists():
            if QMessageBox.question(
                self,
                "Overwrite report?",
                "The selected report file already exists. Overwrite it?",
                QMessageBox.Yes | QMessageBox.No,
            ) != QMessageBox.Yes:
                return
            try:
                report_path.unlink()
            except Exception:
                self._show_error("Memote", "Cannot overwrite the report file. Please choose another path.")
                return

        self.tabs.setCurrentWidget(self.tab_tools)
        self.tools_log.appendPlainText(f"Running memote on: {sbml_path}")
        self.tools_log.appendPlainText(f"Report target: {report_path}\n")
        # Options dialog
        dlg = QDialog(self)
        dlg.setWindowTitle("Memote options")
        layout = QVBoxLayout(dlg)
        form = QFormLayout()

        mode_combo = QComboBox()
        mode_combo.addItems(["Full", "Fast (skip heavy tests)"])
        form.addRow(QLabel("Mode:"), mode_combo)

        timeout_spin = QSpinBox()
        timeout_spin.setRange(1, 120)
        timeout_spin.setValue(15)
        form.addRow(QLabel("Timeout (minutes):"), timeout_spin)

        infinity_chk = QCheckBox("Infinity (no timeout)")
        form.addRow(QLabel(""), infinity_chk)

        layout.addLayout(form)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept)
        btns.rejected.connect(dlg.reject)
        layout.addWidget(btns)

        if dlg.exec() != QDialog.Accepted:
            return

        extra_args = []
        if mode_combo.currentIndex() == 1:
            extra_args = ["--skip", "test_consistency", "--skip", "test_biomass"]

        timeout_seconds = None if infinity_chk.isChecked() else int(timeout_spin.value()) * 60

        self._memote_fast_fallback_used = False
        self._start_memote_run(sbml_path, report_path, extra_args=extra_args or None, timeout_seconds=timeout_seconds)

    def run_carveme_draft(self):
        py = self._tool_python_path()
        if not py:
            self._show_error("Tools missing", "Bundled Python not found. Please reinstall the app.")
            return

        self._carveme_no_draft_retry = False

        if not self._check_tool_module("carveme"):
            msg = (
                "CarveMe is not installed. Install bundled tools now?\n\n"
                "This uses offline wheels packaged with the app."
            )
            if QMessageBox.question(self, "Install tools", msg, QMessageBox.Yes | QMessageBox.No) != QMessageBox.Yes:
                return

            def after_install(ok: bool):
                if not ok:
                    self._show_error("Install failed", "CarveMe installation failed. See Tools log.")
                    return
                if not self._check_tool_module("carveme"):
                    self._show_error("Tools missing", "CarveMe still not available after install.")
                    return
                self.run_carveme_draft()

            self._install_tools_offline(confirm=False, on_done=after_install)
            return

        if not self._carveme_preflight():
            return

        # Input type selection
        dlg = QDialog(self)
        dlg.setWindowTitle("CarveMe input options")
        layout = QVBoxLayout(dlg)
        form = QFormLayout()

        input_combo = QComboBox()
        input_combo.addItems(["Protein FASTA", "DNA FASTA", "RefSeq accession"])
        form.addRow(QLabel("Input type:"), input_combo)

        gapfill_edit = QLineEdit()
        gapfill_edit.setPlaceholderText("e.g., M9,LB")
        form.addRow(QLabel("Gapfill media (-g):"), gapfill_edit)

        init_edit = QLineEdit()
        init_edit.setPlaceholderText("e.g., M9")
        form.addRow(QLabel("Initialize media (-i):"), init_edit)

        species_edit = QLineEdit()
        species_edit.setPlaceholderText("Optional")
        form.addRow(QLabel("Species (--species):"), species_edit)

        skip_rebuild_chk = QCheckBox("Skip rebuild")
        form.addRow(QLabel(""), skip_rebuild_chk)
        skip_biolog_chk = QCheckBox("Skip biolog")
        form.addRow(QLabel(""), skip_biolog_chk)
        skip_ess_chk = QCheckBox("Skip essentiality")
        form.addRow(QLabel(""), skip_ess_chk)

        ignore_warn_chk = QCheckBox("Ignore warnings")
        form.addRow(QLabel(""), ignore_warn_chk)

        layout.addLayout(form)
        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept)
        btns.rejected.connect(dlg.reject)
        layout.addWidget(btns)

        if dlg.exec() != QDialog.Accepted:
            return

        extra_args = []
        if ignore_warn_chk.isChecked():
            extra_args.append("--ignore-warnings")

        if gapfill_edit.text().strip():
            extra_args += ["-g", gapfill_edit.text().strip()]
        if init_edit.text().strip():
            extra_args += ["-i", init_edit.text().strip()]
        if species_edit.text().strip():
            extra_args += ["--species", species_edit.text().strip()]
        if skip_rebuild_chk.isChecked():
            extra_args.append("--skip-rebuild")
        if skip_biolog_chk.isChecked():
            extra_args.append("--skip-biolog")
        if skip_ess_chk.isChecked():
            extra_args.append("--skip-essentiality")

        genome_path = None
        if input_combo.currentIndex() == 2:
            refseq, ok = QInputDialog.getText(self, "RefSeq accession", "Enter RefSeq accession (e.g., GCF_000005845.2):")
            if not ok or not refseq.strip():
                return
            extra_args += ["--refseq", refseq.strip()]
            genome_path = None
        else:
            if input_combo.currentIndex() == 1:
                extra_args.append("--dna")

            genome_path, _ = QFileDialog.getOpenFileName(
                self,
                "Select genome file",
                str(Path.home()),
                "FASTA files (*.fa *.fasta *.faa *.fna);;All files (*.*)"
            )
            if not genome_path:
                return

        out_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save SBML model as",
            str(Path.home() / "carveme_model.xml"),
            "SBML files (*.xml *.sbml);;All files (*.*)"
        )
        if not out_path:
            return

        self.tabs.setCurrentWidget(self.tab_tools)
        self.tools_log.appendPlainText(f"Running CarveMe on: {genome_path}")
        self.tools_log.appendPlainText(f"SBML target: {out_path}\n")

        self._carveme_extra_args = extra_args
        cmd = self._carveme_command(genome_path, out_path, extra_args=extra_args)
        if not cmd:
            self._show_error("Tools missing", "Bundled Python not found. Please reinstall the app.")
            return
        program, args = cmd
        proc = self._run_tool_process(program, args, "carveme")
        self._carveme_proc = proc
        self._carveme_last_output = ""
        self.set_busy(True, "Running CarveMe...")
        self._start_elapsed_timer("carveme")

        def on_finished(code, status):
            self._append_proc_output_to_log(proc)
            self.tools_log.appendPlainText(f"\ncarveme finished with code={code}, status={status}\n")
            self._carveme_proc = None
            self.set_busy(False, "Ready.")
            self._stop_carveme_timer()
            out = getattr(self, "_carveme_last_output", "") or ""
            if code != 0 and "unrecognized arguments" in out.lower() and getattr(self, "_carveme_last_mode", "") != "carveme" and not getattr(self, "_carveme_no_draft_retry", False):
                self._carveme_no_draft_retry = True
                self.tools_log.appendPlainText("\nRetrying CarveMe with 'carveme draft' command...\n")
                self._carveme_last_mode = "carveme"
                args2 = ["draft"]
                if genome_path:
                    args2 += ["-g", str(genome_path)]
                args2 += ["-o", str(out_path)] + (self._carveme_extra_args or [])
                cmd2 = self._console_script_command("carveme", args2)
                if not cmd2 and self._tool_python_path():
                    cmd2 = (self._tool_python_path(), ["-m", "carveme", "draft"] + (["-g", str(genome_path)] if genome_path else []) + ["-o", str(out_path)] + (self._carveme_extra_args or []))
                if cmd2:
                    program2, args2 = cmd2
                    proc2 = self._run_tool_process(program2, args2, "carveme")
                    self._carveme_proc = proc2
                    self._carveme_last_output = ""
                    self.set_busy(True, "Running CarveMe...")
                    self._start_elapsed_timer("carveme")

                    def on_finished2(code2, status2):
                        self._append_proc_output_to_log(proc2)
                        self.tools_log.appendPlainText(f"\ncarveme finished with code={code2}, status={status2}\n")
                        self._carveme_proc = None
                        self.set_busy(False, "Ready.")
                        self._stop_carveme_timer()
                        if code2 == 0 and Path(out_path).exists():
                            if QMessageBox.question(self, "CarveMe", "Model created. Open now?", QMessageBox.Yes | QMessageBox.No) == QMessageBox.Yes:
                                try:
                                    self.open_sbml_from_path(out_path)
                                except Exception as e:
                                    self._show_error("Open failed", "Could not open CarveMe model.", e)
                        else:
                            QMessageBox.warning(self, "CarveMe", "CarveMe failed. See Tools log for details.")

                    proc2.finished.connect(on_finished2)
                    proc2.start()
                return
            if code == 0 and Path(out_path).exists():
                if QMessageBox.question(self, "CarveMe", "Model created. Open now?", QMessageBox.Yes | QMessageBox.No) == QMessageBox.Yes:
                    try:
                        self.open_sbml_from_path(out_path)
                    except Exception as e:
                        self._show_error("Open failed", "Could not open CarveMe model.", e)
            else:
                QMessageBox.warning(self, "CarveMe", "CarveMe failed. See Tools log for details.")

        proc.finished.connect(on_finished)
        proc.start()

    # ---------------- Tools installer helpers ----------------

    def _runtime_dir(self) -> Path:
        """Return the ``runtime/`` directory path."""
        return self._app_base_dir() / "runtime"

    def _wheels_dir(self) -> Path:
        """Return the ``runtime/wheels`` directory path."""
        return self._runtime_dir() / "wheels"

    def _tools_bin_dir(self) -> Path:
        """Return the ``runtime/tools`` directory path."""
        return self._runtime_dir() / "tools"

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
                if obj_val is None or (isinstance(obj_val, float) and math.isnan(obj_val)):
                    issues.append("Objective value is NaN or infeasible (check bounds and objective)")
                elif isinstance(obj_val, float) and abs(obj_val) < 1e-9:
                    issues.append("Objective value is zero – model may be over-constrained")
            except Exception:
                issues.append("Model optimization failed - check for infeasibility")
            
        except Exception as e:
            issues.append(f"Validation error: {e}")
        
        if not issues:
            QMessageBox.information(self, "Validation", "✓ Model validation passed - no major issues found.")
        else:
            msg = "Issues found:\n\n" + "\n".join(issues[:20])
            if len(issues) > 20:
                msg += f"\n... and {len(issues) - 20} more issues"
            box = QMessageBox(self)
            box.setWindowTitle("Model validation")
            box.setText(msg)
            view_btn = box.addButton("View all", QMessageBox.ActionRole)
            box.addButton("Close", QMessageBox.RejectRole)
            box.exec()
            if box.clickedButton() == view_btn:
                self._show_issues_dialog("Model validation issues", issues)

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
                # Use COBRApy's element-based mass balance check
                try:
                    balance = rxn.check_mass_balance()
                    if balance:
                        parts = ", ".join(f"{elem}: {val:+.4g}" for elem, val in balance.items())
                        imbalanced.append(f"{rxn.id}: {parts}")
                except Exception:
                    # Fallback: coefficient sum heuristic (only for reactions without formulas)
                    coeff_sum = sum(rxn.metabolites.values())
                    if abs(coeff_sum) > 0.01:
                        imbalanced.append(f"{rxn.id}: coeff_sum={coeff_sum:.4f} (no formula data)")
            
            if not imbalanced:
                QMessageBox.information(self, "Mass balance check", "✓ All reactions appear mass-balanced.")
            else:
                msg = f"Found {len(imbalanced)} potentially imbalanced reactions:\n\n"
                msg += "\n".join(imbalanced[:15])
                if len(imbalanced) > 15:
                    msg += f"\n... and {len(imbalanced) - 15} more"
                box = QMessageBox(self)
                box.setWindowTitle("Mass balance check")
                box.setText(msg)
                view_btn = box.addButton("View all", QMessageBox.ActionRole)
                box.addButton("Close", QMessageBox.RejectRole)
                box.exec()
                if box.clickedButton() == view_btn:
                    self._show_issues_dialog("Mass balance issues", imbalanced)
        except Exception as e:
            self._show_error("Check failed", "Mass balance check failed.", e)

    def check_gpr_syntax(self):
        """Validate Gene-Protein-Reaction (GPR) syntax."""
        if self.base_model is None:
            QMessageBox.warning(self, "No model", "Load a model first.")
            return

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
                
                # Check for mixed case operators (COBRApy expects lowercase 'and'/'or')
                if re.search(r'\bAND\b', gpr) or re.search(r'\bOR\b', gpr):
                    issues.append(f"{rxn.id}: uppercase operator found (COBRApy expects lowercase 'and'/'or')")
                
                # Validate tokens — split on operators and parentheses
                tokens = re.split(r'\s+|\b(?:and|or)\b|[()]', gpr)
                for tok in tokens:
                    tok = tok.strip()
                    if tok and not re.match(r'^[A-Za-z0-9_.\-:]+$', tok):
                        issues.append(f"{rxn.id}: invalid token '{tok}'")
                        break
            
            if not issues:
                QMessageBox.information(self, "GPR validation", "✓ All GPR expressions appear syntactically valid.")
            else:
                msg = f"Found {len(issues)} GPR syntax issues:\n\n"
                msg += "\n".join(issues[:20])
                if len(issues) > 20:
                    msg += f"\n... and {len(issues) - 20} more"
                box = QMessageBox(self)
                box.setWindowTitle("GPR validation")
                box.setText(msg)
                view_btn = box.addButton("View all", QMessageBox.ActionRole)
                box.addButton("Close", QMessageBox.RejectRole)
                box.exec()
                if box.clickedButton() == view_btn:
                    self._show_issues_dialog("GPR syntax issues", issues)
        except Exception as e:
            self._show_error("Check failed", "GPR syntax check failed.", e)

    def tools_check_status(self):
        py = self._tool_python_path()
        if not py:
            self.tools_status_lbl.setText("Bundled Python not found. Please reinstall the app.")
            self.memote_btn.setEnabled(False)
            self.carveme_btn.setEnabled(False)
            return

        # Run module checks in a background thread so the UI stays responsive
        self._tool_check_worker = _ToolCheckWorker(
            py, ["memote", "carveme", "swiglpk"], parent=self
        )
        self._tool_check_worker.finished.connect(self._on_tool_check_done)
        self._tool_check_worker.start()

    def _on_tool_check_done(self, missing: list[str]):
        """Slot called when background tool-check finishes."""
        if not missing:
            self.tools_status_lbl.setText("✔ All tools available (memote, carveme, swiglpk)")
            self.memote_btn.setEnabled(True)
            self.carveme_btn.setEnabled(True)
        else:
            self.tools_status_lbl.setText(
                "Missing: " + ", ".join(missing)
                + ".  Click 'Repair tools' to reinstall from bundled wheels."
            )
            self.memote_btn.setEnabled(False)
            self.carveme_btn.setEnabled(False)

    def tools_repair(self):
        """Reinstall memote / carveme / swiglpk from bundled wheels."""
        self.tabs.setCurrentWidget(self.tab_tools)
        self._install_tools_offline(confirm=True)

    def _install_tools_offline(self, confirm: bool = True, on_done=None):
        """Install tools from bundled wheels (offline pip install)."""
        py = self._bundled_python_path() or self._tool_python_path()
        if not py:
            self._show_error("Tools missing", "Bundled Python not found. Please reinstall the app.")
            return

        wheels_dir = self._wheels_dir()
        if not wheels_dir.exists():
            self._show_error("Tools missing", "Offline wheels not found in runtime/wheels.")
            return

        if confirm:
            msg = (
                "Reinstall tools from bundled wheels?\n\n"
                f"Source: {wheels_dir}\n\n"
                "This does NOT require internet."
            )
            if QMessageBox.question(self, "Repair tools", msg, QMessageBox.Yes | QMessageBox.No) != QMessageBox.Yes:
                return

        args = [
            "-m", "pip", "install",
            "--no-index", "--find-links", str(wheels_dir),
            "--force-reinstall",
            "memote", "carveme", "swiglpk",
        ]
        proc = self._run_tool_process(py, args, "pip install (offline wheels)")
        self.set_busy(True, "Installing tools from bundled wheels...")

        def on_finished(code, status):
            self._append_proc_output_to_log(proc)
            self.set_busy(False, "Ready.")
            if code == 0:
                try:
                    self._ensure_carveme_shim(py)
                except Exception:
                    pass
                if on_done:
                    on_done(True)
                else:
                    QMessageBox.information(self, "Repair complete", "Tools reinstalled successfully.")
            else:
                if on_done:
                    on_done(False)
                else:
                    self._show_error("Repair failed", f"pip exited with code {code}")
            self.tools_check_status()

        proc.finished.connect(on_finished)
        proc.start()

    # -------- Network Map --------

    # ==================== AUTO-UPDATE CHECKER ====================

    _GITHUB_REPO = "Emiray02/metabodesk"
    _CURRENT_VERSION = "1.0.0"

    def check_for_updates(self, silent: bool = False) -> None:
        """Check GitHub Releases API for a newer version of MetaboDesk.

        Parameters
        ----------
        silent : bool
            If True (startup auto-check), do not show modal dialogs.
            If False (manual check), show dialogs for update status.

        The actual HTTP request runs in a background thread so it never
        blocks the UI.
        """
        self._update_check_silent = silent
        worker = _UpdateCheckWorker(self._GITHUB_REPO, parent=self)
        worker.result_ready.connect(self._on_update_check_done)
        self._update_check_worker = worker  # prevent GC
        worker.start()

    def _on_update_check_done(self, data: dict) -> None:
        """Slot called when background update-check finishes."""
        silent = getattr(self, "_update_check_silent", True)

        if "error" in data:
            if not silent:
                self._show_error(
                    "Update Check Failed",
                    f"Could not reach GitHub:\n{data['error']}\n\nCheck your internet connection.",
                )
            logger.warning("Update check failed: %s", data["error"])
            return

        tag = data.get("tag", "")
        if not tag:
            if not silent:
                QMessageBox.information(self, "Update", "Could not determine latest version.")
            return

        html_url = data.get("url", f"https://github.com/{self._GITHUB_REPO}/releases")
        body = data.get("body", "")

        if self._version_newer(tag, self._CURRENT_VERSION):
            if silent:
                try:
                    self.statusBar().showMessage(
                        f"Update available: v{tag} (current v{self._CURRENT_VERSION}). Help → Check for Updates",
                        10000,
                    )
                except Exception:
                    pass
                logger.info("Update available (silent mode): current=%s latest=%s", self._CURRENT_VERSION, tag)
                QTimer.singleShot(
                    1200,
                    lambda: self._show_update_available_dialog(tag, html_url, body, modal=False),
                )
                return
            self._show_update_available_dialog(tag, html_url, body, modal=True)
        else:
            if not silent:
                QMessageBox.information(
                    self, "Up to Date",
                    f"You are running the latest version (v{self._CURRENT_VERSION}).",
                )

    def _show_update_available_dialog(self, tag: str, html_url: str, body: str, modal: bool) -> None:
        """Show update notification dialog.

        modal=True  -> manual check flow (blocking question dialog)
        modal=False -> startup auto-check flow (non-blocking dialog)
        """
        text = (
            "A new version of MetaboDesk is available!\n\n"
            f"Current: v{self._CURRENT_VERSION}\n"
            f"Latest:  v{tag}\n\n"
            f"{body[:300]}{'…' if len(body) > 300 else ''}\n\n"
            "Open the download page?"
        )

        if modal:
            reply = QMessageBox.question(
                self,
                "Update Available 🎉",
                text,
                QMessageBox.Yes | QMessageBox.No,
            )
            if reply == QMessageBox.Yes:
                QDesktopServices.openUrl(QUrl(html_url))
            return

        box = QMessageBox(self)
        box.setIcon(QMessageBox.Information)
        box.setWindowTitle("Update Available 🎉")
        box.setText(text)
        btn_open = box.addButton("Open Download Page", QMessageBox.AcceptRole)
        box.addButton("Later", QMessageBox.RejectRole)

        def _on_clicked(clicked_btn):
            if clicked_btn == btn_open:
                QDesktopServices.openUrl(QUrl(html_url))

        box.buttonClicked.connect(_on_clicked)
        box.open()  # non-blocking; does not stall startup
        self._update_notice_box = box  # keep ref to avoid GC

    @staticmethod
    def _version_newer(latest: str, current: str) -> bool:
        """Return True if *latest* is strictly newer than *current* using semver tuple comparison."""
        def _parse(v: str) -> tuple[int, ...]:
            parts: list[int] = []
            for p in v.split("."):
                try:
                    parts.append(int(p))
                except ValueError:
                    break
            return tuple(parts) if parts else (0,)
        return _parse(latest) > _parse(current)

