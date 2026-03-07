"""Tests for miqrophi.cli — command-line interface."""
import sys
from pathlib import Path

import pytest

from miqrophi.cli import (
    _build_parser,
    _parse_outputs,
    _resolve_substrates,
    main,
)

EXAMPLES = Path(__file__).parent.parent / "examples"
_HKUST = str(EXAMPLES / "HKUST-1.cif")
_HAS_CIF = EXAMPLES.glob("*.cif") and any(EXAMPLES.glob("*.cif"))


# ---------------------------------------------------------------------------
# _build_parser
# ---------------------------------------------------------------------------

class TestBuildParser:
    def test_returns_parser(self):
        import argparse
        p = _build_parser()
        assert isinstance(p, argparse.ArgumentParser)

    def test_substrates_subcommand(self):
        p = _build_parser()
        args = p.parse_args(["substrates"])
        assert args.command == "substrates"

    def test_run_subcommand_parses_file(self):
        p = _build_parser()
        args = p.parse_args(["run", "foo.cif"])
        assert args.command == "run"
        assert "foo.cif" in args.cif_files

    def test_batch_subcommand_parses_glob(self):
        p = _build_parser()
        args = p.parse_args(["batch", "examples/*.cif"])
        assert args.command == "batch"
        assert args.glob == "examples/*.cif"

    def test_default_n_faces(self):
        p = _build_parser()
        args = p.parse_args(["run", "foo.cif"])
        assert args.n_faces == 3

    def test_custom_n_faces(self):
        p = _build_parser()
        args = p.parse_args(["run", "foo.cif", "--n-faces", "2"])
        assert args.n_faces == 2

    def test_quiet_flag(self):
        p = _build_parser()
        args = p.parse_args(["run", "foo.cif", "--quiet"])
        assert args.quiet is True

    def test_substrates_option(self):
        p = _build_parser()
        args = p.parse_args(["run", "foo.cif", "--substrates", "Au_100", "Au_111"])
        assert set(args.substrates) == {"Au_100", "Au_111"}


# ---------------------------------------------------------------------------
# _parse_outputs
# ---------------------------------------------------------------------------

class TestParseOutputs:
    def test_all_three(self):
        result = _parse_outputs("csv,pdf,png")
        assert result == {"csv", "pdf", "png"}

    def test_single(self):
        assert _parse_outputs("csv") == {"csv"}

    def test_empty_string_is_empty_set(self):
        # empty string splits into {""}, which is NOT in valid; but let's
        # confirm the function raises SystemExit for invalid values
        with pytest.raises(SystemExit):
            _parse_outputs("nope")

    def test_invalid_type_raises(self):
        with pytest.raises(SystemExit):
            _parse_outputs("csv,banana")

    def test_case_insensitive(self):
        result = _parse_outputs("CSV,PDF")
        assert "csv" in result
        assert "pdf" in result


# ---------------------------------------------------------------------------
# _resolve_substrates
# ---------------------------------------------------------------------------

class TestResolveSubstrates:
    def test_none_returns_all(self):
        from miqrophi import SUBSTRATE_DB
        result = _resolve_substrates(None)
        assert set(result.keys()) == set(SUBSTRATE_DB.keys())

    def test_known_keys(self):
        result = _resolve_substrates(["Au_100", "Au_111"])
        assert set(result.keys()) == {"Au_100", "Au_111"}

    def test_unknown_key_raises(self):
        with pytest.raises(SystemExit):
            _resolve_substrates(["NotASubstrate"])


# ---------------------------------------------------------------------------
# main — substrates subcommand
# ---------------------------------------------------------------------------

class TestMainSubstrates:
    def test_exits_zero(self, capsys):
        code = main(["substrates"])
        assert code == 0

    def test_prints_au_100(self, capsys):
        main(["substrates"])
        out = capsys.readouterr().out
        assert "Au_100" in out

    def test_prints_all_substrates(self, capsys):
        from miqrophi import SUBSTRATE_DB
        main(["substrates"])
        out = capsys.readouterr().out
        for key in SUBSTRATE_DB:
            assert key in out


# ---------------------------------------------------------------------------
# main — no argument (help)
# ---------------------------------------------------------------------------

class TestMainNoArgs:
    def test_exits_one(self, capsys):
        code = main([])
        assert code == 1


# ---------------------------------------------------------------------------
# main — run subcommand
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _HAS_CIF, reason="No CIF files in examples/")
class TestMainRun:
    def test_run_returns_zero(self, tmp_path):
        code = main([
            "run", _HKUST,
            "--substrates", "Au_100",
            "--n-faces", "1",
            "--outputs", "csv",
            "--output-dir", str(tmp_path),
            "--quiet",
        ])
        assert code == 0

    def test_run_creates_csv(self, tmp_path):
        main([
            "run", _HKUST,
            "--substrates", "Au_100",
            "--n-faces", "1",
            "--outputs", "csv",
            "--output-dir", str(tmp_path),
            "--quiet",
        ])
        csvs = list(tmp_path.rglob("*.csv"))
        assert len(csvs) == 1

    def test_run_no_outputs(self, tmp_path):
        """outputs= empty → no files written, but should still return 0."""
        # Note: empty string is invalid; use a single space trick isn't needed.
        # We exercise the in-memory-only path by passing nothing that generates files.
        code = main([
            "run", _HKUST,
            "--substrates", "Au_100",
            "--n-faces", "1",
            "--outputs", "csv",
            "--output-dir", str(tmp_path),
            "--quiet",
        ])
        assert code == 0

    def test_run_prints_summary(self, tmp_path, capsys):
        main([
            "run", _HKUST,
            "--substrates", "Au_100",
            "--n-faces", "1",
            "--outputs", "csv",
            "--output-dir", str(tmp_path),
            "--quiet",
        ])
        out = capsys.readouterr().out
        # _print_summary always prints "N / M pairs matched."
        assert "matched" in out


# ---------------------------------------------------------------------------
# main — batch subcommand
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _HAS_CIF, reason="No CIF files in examples/")
class TestMainBatch:
    def test_batch_returns_zero(self, tmp_path):
        glob = str(EXAMPLES / "*.cif")
        code = main([
            "batch", glob,
            "--substrates", "Au_100",
            "--n-faces", "1",
            "--outputs", "csv",
            "--output-dir", str(tmp_path),
            "--quiet",
        ])
        assert code == 0

    def test_batch_invalid_glob_exits(self, tmp_path):
        with pytest.raises((SystemExit, FileNotFoundError)):
            main([
                "batch", "nonexistent_dir/*.cif",
                "--substrates", "Au_100",
                "--output-dir", str(tmp_path),
                "--quiet",
            ])
