"""Reusable Qt widgets for MetaboDesk.

- :class:`TextPopup` — modal dialog showing read-only text (cell content).
- :class:`AnimatedSplash` — frameless GIF splash screen shown at startup.
- :class:`MplCanvas` — matplotlib ``FigureCanvas`` embedded in PySide6.
  Matplotlib is imported lazily on first :class:`MplCanvas` instantiation
  so that startup is not blocked by the ~0.4 s matplotlib import cost.
- :class:`AnalysisWorker` — ``QThread``-based background worker for
  long-running analyses with progress and error signalling.
"""

import inspect
import traceback

from PySide6.QtCore import Qt, Signal, QThread
from PySide6.QtGui import QMovie
from PySide6.QtWidgets import (
    QDialog, QDialogButtonBox, QLabel, QPlainTextEdit, QVBoxLayout, QWidget,
)

# ── Lazy matplotlib imports — loaded on first MplCanvas creation ──────
_FigureCanvas = None
_Figure = None


def _ensure_matplotlib():
    """Import matplotlib on first use, not at module load time."""
    global _FigureCanvas, _Figure
    if _FigureCanvas is None:
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.figure import Figure
        _FigureCanvas = FigureCanvasQTAgg
        _Figure = Figure


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


class AnimatedSplash(QDialog):
    """Frameless animated GIF splash screen.

    The GIF loops continuously (restarts from the beginning) until
    :meth:`close` is called — i.e. until the main window is ready.
    """

    def __init__(self, gif_path: str, parent=None):
        super().__init__(parent)
        self.setWindowFlags(Qt.FramelessWindowHint | Qt.WindowStaysOnTopHint)
        self.setAttribute(Qt.WA_TranslucentBackground)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.label = QLabel()
        self.label.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.label)

        self.movie = QMovie(gif_path)
        # Loop forever (-1) so the animation restarts from the
        # beginning and keeps playing until the app is fully loaded.
        self.movie.setSpeed(100)        # normal speed
        self.movie.setCacheMode(QMovie.CacheAll)
        # Scale to a reasonable size (400×400) so it doesn't fill the screen
        from PySide6.QtCore import QSize
        self.movie.setScaledSize(QSize(400, 400))
        self.label.setMovie(self.movie)
        self.movie.start()
        self.setFixedSize(400, 400)

        # Center on primary screen
        from PySide6.QtWidgets import QApplication
        screen = QApplication.primaryScreen()
        if screen:
            geo = screen.availableGeometry()
            x = (geo.width() - self.width()) // 2 + geo.x()
            y = (geo.height() - self.height()) // 2 + geo.y()
            self.move(x, y)

    def close(self):
        """Stop the animation and close the splash."""
        self.movie.stop()
        super().close()


class MplCanvas(QWidget):
    """Matplotlib canvas that defers the heavy matplotlib import.

    On first creation the class lazily loads ``FigureCanvasQTAgg`` and
    ``Figure``.  A lightweight ``QWidget`` placeholder is used as the
    base class so that the widget can be inserted into layouts at startup
    without triggering an import.  The real canvas is initialised on the
    first call that needs it (paint, draw, or attribute access).
    """

    def __init__(self, parent=None, width=8, height=3.5, dpi=100):
        super().__init__(parent)
        self._mpl_width = width
        self._mpl_height = height
        self._mpl_dpi = dpi
        self._canvas = None  # real FigureCanvasQTAgg — created lazily
        self._ax = None
        self.clicked_node = None
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)

    def _ensure_canvas(self):
        """Materialise the real matplotlib canvas on first use."""
        if self._canvas is not None:
            return
        _ensure_matplotlib()
        fig = _Figure(figsize=(self._mpl_width, self._mpl_height),
                      dpi=self._mpl_dpi, tight_layout=True)
        self._ax = fig.add_subplot(111)
        self._canvas = _FigureCanvas(fig)
        self._canvas.setParent(self)
        self._layout.addWidget(self._canvas)

    @property
    def ax(self):
        """The matplotlib Axes — materialises the canvas on first access."""
        self._ensure_canvas()
        return self._ax

    @ax.setter
    def ax(self, value):
        self._ax = value

    # Forward common matplotlib canvas methods to the real canvas
    def draw(self):
        self._ensure_canvas()
        self._canvas.draw()

    def draw_idle(self):
        self._ensure_canvas()
        self._canvas.draw_idle()

    @property
    def figure(self):
        self._ensure_canvas()
        return self._canvas.figure

    def mpl_connect(self, event, callback):
        self._ensure_canvas()
        return self._canvas.mpl_connect(event, callback)

    def mpl_disconnect(self, cid):
        self._ensure_canvas()
        return self._canvas.mpl_disconnect(cid)

    def get_default_filename(self):
        self._ensure_canvas()
        return self._canvas.get_default_filename()

    def print_figure(self, *args, **kwargs):
        self._ensure_canvas()
        return self._canvas.print_figure(*args, **kwargs)


class AnalysisWorker(QThread):
    """Background worker for long-running analyses.  Prevents UI freezing."""
    finished = Signal(dict)
    error = Signal(str)
    progress = Signal(str)
    progress_pct = Signal(int)

    def __init__(self, func, **kwargs):
        super().__init__()
        self._func = func
        self._kwargs = kwargs

    def report_progress(self, message: str, pct: int = -1):
        """Emit progress from within the worker function."""
        self.progress.emit(message)
        if pct >= 0:
            self.progress_pct.emit(pct)

    def run(self):
        try:
            sig = inspect.signature(self._func)
            if 'worker' in sig.parameters:
                result = self._func(worker=self, **self._kwargs)
            else:
                result = self._func(**self._kwargs)
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(f"{e}\n\n{traceback.format_exc()}")
