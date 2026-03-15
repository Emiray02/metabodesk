"""Persistent application configuration for MetaboDesk.

Manages user preferences via a TOML configuration file stored at
``~/.metabodesk/config.toml``.  Falls back to sensible defaults when
the file does not exist or cannot be parsed.

Usage
-----
::

    from metabodesk_core.config import AppConfig

    cfg = AppConfig.load()          # read (or create) config file
    cfg.solver = "glpk"
    cfg.save()                      # persist to disk

Supported keys
--------------
- ``solver`` — default LP solver (``"Auto"``, ``"glpk"``, ``"highs"``,
  ``"gurobi"``, ``"cplex"``).
- ``theme`` — UI theme (``"system"``, ``"dark"``, ``"light"``).
- ``top_n`` — default Top-N value for charts/tables.
- ``fva_fraction`` — default fraction of optimum for FVA.
- ``fva_processes`` — default number of parallel FVA processes.
- ``loopless`` — whether to use loopless FBA by default.
- ``sampling_size`` — default number of samples for flux sampling.
- ``sampler`` — default sampler type (``"ACHR"`` or ``"OptGP"``).
- ``map_max_nodes`` — max nodes in network map.
- ``map_layout`` — default layout algorithm for network map.
- ``recent_max`` — maximum number of recent scenarios to remember.
- ``undo_max`` — maximum undo stack depth.
- ``cache_max`` — maximum FBA result cache size.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field, fields
from pathlib import Path

logger = logging.getLogger("MetaboDesk")

_CONFIG_DIR = Path.home() / ".metabodesk"
_CONFIG_PATH = _CONFIG_DIR / "config.toml"
_CURRENT_CONFIG_VERSION = 1


@dataclass
class AppConfig:
    """Application-wide configuration with typed defaults.

    Every attribute corresponds to a key in the TOML file.  Unknown keys
    in the file are silently ignored.  When the ``config_version`` in
    the file is older than ``_CURRENT_CONFIG_VERSION``, the
    :meth:`_migrate` hook is called and the file is re-saved.
    """

    # -- Version (for migration) ---------------------------------------------
    config_version: int = _CURRENT_CONFIG_VERSION

    # -- Solver & Analysis defaults ------------------------------------------
    solver: str = "Auto"
    loopless: bool = False
    top_n: int = 20
    fva_fraction: float = 1.0
    fva_processes: int = 1
    sampling_size: int = 500
    sampler: str = "ACHR"

    # -- UI / Theme -----------------------------------------------------------
    theme: str = "system"          # "system", "dark", "light"

    # -- Network map defaults -------------------------------------------------
    map_max_nodes: int = 1000
    map_layout: str = "Spring"     # "Spring", "Kamada-Kawai", "Circular"
    map_depth: int = 2
    map_threshold: float = 0.1

    # -- Capacity limits ------------------------------------------------------
    recent_max: int = 5
    undo_max: int = 50
    cache_max: int = 32

    # -- Internal (not persisted) --------------------------------------------
    _path: Path = field(default=_CONFIG_PATH, repr=False, compare=False)

    # ------------------------------------------------------------------ #
    #  Persistence                                                         #
    # ------------------------------------------------------------------ #

    @classmethod
    def load(cls, path: Path | str | None = None) -> "AppConfig":
        """Load configuration from *path* (default ``~/.metabodesk/config.toml``).

        If the file does not exist or is unparseable, returns an instance
        with default values and writes a fresh config file.

        Parameters
        ----------
        path : Path or str, optional
            Override the default config file location.

        Returns
        -------
        AppConfig
            Populated configuration instance.
        """
        resolved = Path(path) if path else _CONFIG_PATH
        cfg = cls(_path=resolved)
        if resolved.exists():
            try:
                data = _read_toml(resolved)
                valid_names = {f.name for f in fields(cls) if not f.name.startswith("_")}
                for key, value in data.items():
                    if key in valid_names:
                        try:
                            expected = type(getattr(cfg, key))
                            setattr(cfg, key, expected(value))
                        except (TypeError, ValueError) as exc:
                            logger.warning("Config key %r: bad value %r (%s)", key, value, exc)
            except Exception as exc:
                logger.warning("Could not parse config %s: %s", resolved, exc)
        else:
            # Create default config on first run
            cfg.save()
        # Migrate old config versions
        if cfg.config_version < _CURRENT_CONFIG_VERSION:
            cfg._migrate(cfg.config_version)
            cfg.config_version = _CURRENT_CONFIG_VERSION
            cfg.save()
        return cfg

    def _migrate(self, from_version: int) -> None:
        """Apply configuration migrations from *from_version* to current.

        Each future version bump should add a migration step here.
        """
        # Example for future use:
        # if from_version < 2:
        #     self.new_setting = "default"
        logger.info("Config migrated from v%d to v%d", from_version, _CURRENT_CONFIG_VERSION)

    def save(self, path: Path | str | None = None) -> None:
        """Write configuration to disk in TOML format.

        Parameters
        ----------
        path : Path or str, optional
            Override the stored path.
        """
        resolved = Path(path) if path else self._path
        try:
            resolved.parent.mkdir(parents=True, exist_ok=True)
            lines = [
                "# MetaboDesk configuration file",
                "# Edit values below to change application defaults.",
                "# Delete this file to reset all settings.",
                "",
            ]
            for f in fields(self):
                if f.name.startswith("_"):
                    continue
                value = getattr(self, f.name)
                lines.append(_format_toml_line(f.name, value))
            resolved.write_text("\n".join(lines) + "\n", encoding="utf-8")
            logger.info("Config saved to %s", resolved)
        except Exception as exc:
            logger.error("Could not save config to %s: %s", resolved, exc)

    def as_dict(self) -> dict:
        """Return a plain dict of all public settings."""
        return {f.name: getattr(self, f.name) for f in fields(self) if not f.name.startswith("_")}


# ------------------------------------------------------------------ #
#  Minimal TOML helpers (stdlib-only, no 3rd-party dependency)         #
# ------------------------------------------------------------------ #

def _read_toml(path: Path) -> dict:
    """Parse a *simple* TOML file (flat key = value only).

    Handles strings, integers, floats, and booleans.  Does **not**
    support tables, arrays-of-tables, multi-line strings, or inline
    tables — but those are not needed for MetaboDesk's config.
    """
    result: dict = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, _, val = line.partition("=")
        key = key.strip()
        val = val.strip()
        # Remove inline comments
        if "#" in val and not val.startswith('"') and not val.startswith("'"):
            val = val[:val.index("#")].strip()
        result[key] = _parse_toml_value(val)
    return result


def _parse_toml_value(val: str):
    """Coerce a raw TOML value string to a Python type."""
    if val.lower() == "true":
        return True
    if val.lower() == "false":
        return False
    # Quoted string
    if (val.startswith('"') and val.endswith('"')) or (val.startswith("'") and val.endswith("'")):
        return val[1:-1]
    # Integer / Float
    try:
        if "." in val:
            return float(val)
        return int(val)
    except ValueError:
        return val  # return raw string as fallback


def _format_toml_line(key: str, value) -> str:
    """Format a key-value pair as a TOML line."""
    if isinstance(value, bool):
        return f"{key} = {'true' if value else 'false'}"
    if isinstance(value, str):
        return f'{key} = "{value}"'
    if isinstance(value, float):
        return f"{key} = {value}"
    if isinstance(value, int):
        return f"{key} = {value}"
    return f'{key} = "{value}"'
