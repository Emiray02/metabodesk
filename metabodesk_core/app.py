"""Application entry point.

Creates the ``QApplication``, shows the animated splash screen, constructs
the :class:`~metabodesk_core.mainwindow.MainWindow`, and enters the Qt
event loop.  Called from the thin ``metabodesk.py`` launcher or via
``python -m metabodesk_core``.

IMPORTANT — Only lightweight PySide6 modules are imported at module level.
Heavy libraries (cobra, matplotlib, …) are deferred to *after* the splash
screen is visible so the user gets instant visual feedback.
"""

import os
import sys
import logging
from pathlib import Path

# ── Only lightweight Qt imports — no cobra / matplotlib here ──
from PySide6.QtWidgets import QApplication, QSplashScreen
from PySide6.QtGui import QIcon, QPixmap, QMovie
from PySide6.QtCore import Qt, QTimer, QSize

logger = logging.getLogger("MetaboDesk")


def _startup_trace(msg: str) -> None:
    """Append lightweight startup diagnostics to a writable location."""
    try:
        # Prefer APPDATA\MetaboDesk (always writable, even in Program Files installs)
        appdata = os.environ.get("APPDATA") or os.path.expanduser("~")
        log_dir = Path(appdata) / "MetaboDesk"
        log_dir.mkdir(parents=True, exist_ok=True)
        p = log_dir / "metabodesk_startup.log"
        import datetime
        ts = datetime.datetime.now().strftime("%H:%M:%S.%f")[:-3]
        with p.open("a", encoding="utf-8") as f:
            f.write(f"[{ts}] {msg}\n")
    except Exception:
        pass


def _resolve_resource(rel_path: str) -> Path:
    """Lightweight resource path resolver (no heavy dependencies)."""
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS) / rel_path
    base = os.environ.get("METABODESK_BASE")
    if base:
        return Path(base) / rel_path
    return Path(__file__).resolve().parent.parent / rel_path


class _AnimatedSplash(QSplashScreen):
    """QSplashScreen subclass that plays an animated GIF.

    Uses Qt's native QSplashScreen — avoids the black-rectangle bug
    that occurs with QDialog + WA_TranslucentBackground under pythonw.
    Each GIF frame is painted via ``setPixmap`` so the splash always
    displays correctly regardless of compositor or interpreter.
    """

    def __init__(self, gif_path: Path):
        # Start with the first frame as the initial pixmap
        self._movie = QMovie(str(gif_path))
        self._movie.setScaledSize(QSize(400, 400))
        self._movie.setCacheMode(QMovie.CacheAll)
        self._movie.jumpToFrame(0)
        first_frame = self._movie.currentPixmap()
        if first_frame.isNull():
            first_frame = QPixmap(400, 400)
            first_frame.fill(Qt.transparent)

        super().__init__(first_frame, Qt.WindowStaysOnTopHint)
        self.setFixedSize(400, 400)

        # Each time the movie advances, repaint with the new frame
        self._movie.frameChanged.connect(self._on_frame)
        self._movie.start()

    def _on_frame(self, _frame_number: int):
        pm = self._movie.currentPixmap()
        if not pm.isNull():
            self.setPixmap(pm)

    def close(self):
        self._movie.stop()
        super().close()


def _setup_crash_log():
    """Redirect stderr to a log file so pythonw crashes are not lost."""
    try:
        # Write to APPDATA (always writable, works even from Program Files)
        appdata = os.environ.get("APPDATA") or os.path.expanduser("~")
        log_dir = Path(appdata) / "MetaboDesk"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / "metabodesk_error.log"
        # Append mode so we keep history
        sys.stderr = open(log_file, "a", encoding="utf-8")
    except Exception:
        pass


def main():
    _setup_crash_log()
    _startup_trace("main: start")

    # ── Tell Windows this is "MetaboDesk", not "pythonw.exe" ──────────
    # Without this, the taskbar groups our window under Python's icon
    # and Task-Manager's "Apps" tab shows "Python" instead of our title.
    try:
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(
            "Emiray.MetaboDesk.App.1"   # unique per-app string
        )
    except Exception:
        pass

    app = QApplication(sys.argv)
    app.setApplicationName("MetaboDesk")
    app.setApplicationDisplayName("MetaboDesk")
    app._startup_window_shown = False  # type: ignore[attr-defined]

    # ── Set the application icon FIRST, before any window is created ──
    try:
        icon_path = _resolve_resource("logo.ico")
        if icon_path.exists():
            app.setWindowIcon(QIcon(str(icon_path)))
    except Exception:
        pass

    # ── Show splash IMMEDIATELY (no heavy imports yet) ──
    splash = None
    try:
        splash_path = _resolve_resource("splash.gif")
        if splash_path.exists():
            splash = _AnimatedSplash(splash_path)
            splash.show()
            app.processEvents()
            _startup_trace("main: splash shown")
    except Exception:
        splash = None
        _startup_trace("main: splash unavailable")

    # ── Defer ALL heavy work until the event loop is running ──
    def _build_and_show():
        _startup_trace("build: begin")
        try:
            if splash is not None:
                splash.showMessage("Loading modules…",
                                   Qt.AlignBottom | Qt.AlignHCenter,
                                   Qt.white)
                app.processEvents()

            from metabodesk_core.utils import apply_custom_message_boxes
            apply_custom_message_boxes()

            if splash is not None:
                splash.showMessage("Building interface…",
                                   Qt.AlignBottom | Qt.AlignHCenter,
                                   Qt.white)
                app.processEvents()

            from metabodesk_core.mainwindow import MainWindow
            win = MainWindow()

            try:
                import pyi_splash          # type: ignore[import-not-found]
                pyi_splash.close()
            except Exception:
                pass

            win.show()
            win.raise_()
            win.activateWindow()
            app._startup_window_shown = True  # type: ignore[attr-defined]
            _startup_trace("build: main window shown")
        except Exception as exc:
            logger.exception("Fatal error during startup")
            _startup_trace(f"build: exception: {exc}")
            if splash is not None:
                splash.close()
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.critical(
                None, "MetaboDesk — Startup Error",
                f"Failed to start MetaboDesk:\n\n{exc}\n\n"
                "See %APPDATA%\\MetaboDesk\\metabodesk_error.log for details."
            )
            app.quit()
            return

        if splash is not None:
            splash.close()
            _startup_trace("build: splash closed")

        app._main_window = win         # type: ignore[attr-defined]

    def _startup_watchdog():
        """Fallback if startup gets stuck with no visible top-level window."""
        try:
            if getattr(app, "_startup_window_shown", False):
                return
            for widget in app.topLevelWidgets():
                if widget.isVisible() and widget.objectName() != "":
                    return
        except Exception:
            pass

        _startup_trace("watchdog: no visible window, forcing fallback show")
        try:
            if splash is not None:
                splash.close()
            from metabodesk_core.mainwindow import MainWindow
            win2 = MainWindow()
            win2.show()
            win2.raise_()
            win2.activateWindow()
            app._main_window = win2  # type: ignore[attr-defined]
            app._startup_window_shown = True  # type: ignore[attr-defined]
            _startup_trace("watchdog: fallback main window shown")
        except Exception as exc:
            logger.exception("Startup watchdog fallback failed")
            _startup_trace(f"watchdog: exception: {exc}")

    def _ensure_window_visible():
        """Force the main window to become visible and active.

        Some Windows environments can keep a just-created Qt window hidden,
        minimized, or off-screen when launched via pythonw from Program Files.
        This helper re-normalizes the state and re-activates the window.
        """
        try:
            win = getattr(app, "_main_window", None)
            if win is None:
                _startup_trace("ensure: no main window ref")
                return

            if win.isMinimized():
                win.showNormal()
            else:
                win.show()

            # Guarantee active/non-minimized state
            try:
                win.setWindowState((win.windowState() & ~Qt.WindowMinimized) | Qt.WindowActive)
            except Exception:
                pass

            # If geometry is invalid or tiny, give it a sane default size
            try:
                geom = win.frameGeometry()
                if geom.width() < 200 or geom.height() < 150:
                    win.resize(1700, 1120)
            except Exception:
                pass

            win.raise_()
            win.activateWindow()
            app.processEvents()
            _startup_trace("ensure: main window re-activated")
        except Exception as exc:
            _startup_trace(f"ensure: exception: {exc}")

    _build_and_show()
    QTimer.singleShot(1000, _ensure_window_visible)
    QTimer.singleShot(6000, _ensure_window_visible)
    QTimer.singleShot(12000, _startup_watchdog)
    _startup_trace("main: event loop enter")

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
