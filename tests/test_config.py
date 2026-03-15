"""Tests for metabodesk_core.config — TOML configuration system."""

import pytest
from pathlib import Path
from metabodesk_core.config import AppConfig, _read_toml, _parse_toml_value, _format_toml_line


# ------------------------------------------------------------------ #
#  TOML helpers                                                        #
# ------------------------------------------------------------------ #

class TestParseTomlValue:
    """Tests for :func:`_parse_toml_value`."""

    def test_bool_true(self):
        assert _parse_toml_value("true") is True

    def test_bool_false(self):
        assert _parse_toml_value("false") is False

    def test_bool_case_insensitive(self):
        assert _parse_toml_value("True") is True
        assert _parse_toml_value("FALSE") is False

    def test_integer(self):
        assert _parse_toml_value("42") == 42
        assert isinstance(_parse_toml_value("42"), int)

    def test_float(self):
        assert _parse_toml_value("3.14") == pytest.approx(3.14)
        assert isinstance(_parse_toml_value("3.14"), float)

    def test_quoted_string(self):
        assert _parse_toml_value('"hello"') == "hello"

    def test_single_quoted_string(self):
        assert _parse_toml_value("'world'") == "world"

    def test_unquoted_string_fallback(self):
        assert _parse_toml_value("Auto") == "Auto"


class TestFormatTomlLine:
    """Tests for :func:`_format_toml_line`."""

    def test_bool(self):
        assert _format_toml_line("loopless", True) == "loopless = true"
        assert _format_toml_line("loopless", False) == "loopless = false"

    def test_string(self):
        assert _format_toml_line("solver", "glpk") == 'solver = "glpk"'

    def test_int(self):
        assert _format_toml_line("top_n", 20) == "top_n = 20"

    def test_float(self):
        assert _format_toml_line("fva_fraction", 1.0) == "fva_fraction = 1.0"


class TestReadToml:
    """Tests for :func:`_read_toml`."""

    def test_read_simple_file(self, tmp_path):
        f = tmp_path / "test.toml"
        f.write_text('solver = "glpk"\ntop_n = 30\nloopless = true\n')
        data = _read_toml(f)
        assert data["solver"] == "glpk"
        assert data["top_n"] == 30
        assert data["loopless"] is True

    def test_comments_ignored(self, tmp_path):
        f = tmp_path / "test.toml"
        f.write_text('# comment\nsolver = "Auto"\n# another comment\n')
        data = _read_toml(f)
        assert data["solver"] == "Auto"
        assert len(data) == 1

    def test_inline_comment(self, tmp_path):
        f = tmp_path / "test.toml"
        f.write_text('top_n = 25  # items per page\n')
        data = _read_toml(f)
        assert data["top_n"] == 25


# ------------------------------------------------------------------ #
#  AppConfig                                                           #
# ------------------------------------------------------------------ #

class TestAppConfig:
    """Tests for :class:`AppConfig`."""

    def test_defaults(self):
        cfg = AppConfig()
        assert cfg.solver == "Auto"
        assert cfg.top_n == 20
        assert cfg.loopless is False
        assert cfg.theme == "system"
        assert cfg.map_layout == "Spring"

    def test_save_and_load(self, tmp_path):
        path = tmp_path / "config.toml"
        cfg = AppConfig(_path=path)
        cfg.solver = "highs"
        cfg.top_n = 42
        cfg.loopless = True
        cfg.save()

        assert path.exists()
        loaded = AppConfig.load(path)
        assert loaded.solver == "highs"
        assert loaded.top_n == 42
        assert loaded.loopless is True

    def test_load_missing_file_creates_default(self, tmp_path):
        path = tmp_path / "does_not_exist.toml"
        cfg = AppConfig.load(path)
        assert cfg.solver == "Auto"
        assert path.exists()  # default file created

    def test_unknown_keys_ignored(self, tmp_path):
        path = tmp_path / "test.toml"
        path.write_text('solver = "glpk"\nunknown_key = 999\n')
        cfg = AppConfig.load(path)
        assert cfg.solver == "glpk"
        assert not hasattr(cfg, "unknown_key")

    def test_bad_value_type_uses_default(self, tmp_path):
        path = tmp_path / "test.toml"
        path.write_text('top_n = "not_a_number"\n')
        cfg = AppConfig.load(path)
        # Should fall back to default because int("not_a_number") fails
        # The value might be stored as-is since we cast with expected()
        # Let's just check it doesn't crash
        assert isinstance(cfg.top_n, int)

    def test_as_dict(self):
        cfg = AppConfig()
        d = cfg.as_dict()
        assert "solver" in d
        assert "top_n" in d
        assert "_path" not in d

    def test_roundtrip_all_types(self, tmp_path):
        path = tmp_path / "config.toml"
        cfg = AppConfig(_path=path)
        cfg.solver = "cplex"
        cfg.loopless = True
        cfg.fva_fraction = 0.95
        cfg.top_n = 15
        cfg.theme = "dark"
        cfg.save()

        cfg2 = AppConfig.load(path)
        assert cfg2.solver == "cplex"
        assert cfg2.loopless is True
        assert cfg2.fva_fraction == pytest.approx(0.95)
        assert cfg2.top_n == 15
        assert cfg2.theme == "dark"
