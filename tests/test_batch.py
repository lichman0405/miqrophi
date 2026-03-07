"""Tests for batch_run pipeline."""
import pytest
from pathlib import Path

import pandas as pd

from miqrophi import batch_run, BatchConfig

EXAMPLES = Path(__file__).parent.parent / "examples"
_CIFS = list(EXAMPLES.glob("*.cif"))


def _cfg(**kw):
    """Minimal no-file, single-substrate, single-face BatchConfig for speed."""
    defaults = dict(outputs=set(), n_faces=1, substrates=["Au_100"], verbose=False)
    defaults.update(kw)
    return BatchConfig(**defaults)


# ── Return type ──────────────────────────────────────────────────────────────

@pytest.mark.skipif(not _CIFS, reason="No CIF files in examples/")
def test_returns_dataframe():
    df = batch_run(_CIFS, config=_cfg())
    assert isinstance(df, pd.DataFrame)


@pytest.mark.skipif(not _CIFS, reason="No CIF files in examples/")
def test_expected_columns():
    df = batch_run(_CIFS, config=_cfg())
    required = {
        "cif_file", "mof_name", "hkl", "d_hkl_A",
        "mof_a", "mof_b", "mof_gamma",
        "substrate", "match_found",
        "theta_deg", "eta", "eps_11", "eps_22", "eps_12",
        "area_A2", "L0_feasible", "png_path", "pdf_path",
    }
    assert required.issubset(df.columns)


@pytest.mark.skipif(not _CIFS, reason="No CIF files in examples/")
def test_row_count_is_n_cifs_times_substrates():
    substrates = ["Au_100", "Au_111"]
    df = batch_run(_CIFS, config=_cfg(substrates=substrates, n_faces=1))
    assert len(df) == len(_CIFS) * len(substrates)


# ── Glob pattern ─────────────────────────────────────────────────────────────

def test_glob_pattern():
    df = batch_run(str(EXAMPLES / "*.cif"), config=_cfg())
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0


def test_glob_no_match_raises():
    with pytest.raises(FileNotFoundError):
        batch_run("nonexistent_dir/*.cif", config=_cfg())


# ── CSV output ───────────────────────────────────────────────────────────────

@pytest.mark.skipif(not _CIFS, reason="No CIF files in examples/")
def test_csv_created(tmp_path):
    df = batch_run(_CIFS, config=_cfg(outputs={"csv"}, output_dir=str(tmp_path)))
    csv_files = list(tmp_path.rglob("summary.csv"))
    assert len(csv_files) == 1
    loaded = pd.read_csv(csv_files[0])
    assert len(loaded) == len(df)


@pytest.mark.skipif(not _CIFS, reason="No CIF files in examples/")
def test_no_csv_when_outputs_empty(tmp_path):
    batch_run(_CIFS, config=_cfg(outputs=set(), output_dir=str(tmp_path)))
    csv_files = list(tmp_path.rglob("*.csv"))
    assert len(csv_files) == 0


# ── Run directory name ───────────────────────────────────────────────────────

@pytest.mark.skipif(not _CIFS, reason="No CIF files in examples/")
def test_run_dir_has_timestamp(tmp_path):
    batch_run(_CIFS, config=_cfg(output_dir=str(tmp_path)))
    subdirs = [d for d in tmp_path.iterdir() if d.is_dir()]
    assert len(subdirs) == 1
    assert subdirs[0].name.startswith("run_")


@pytest.mark.skipif(not _CIFS, reason="No CIF files in examples/")
def test_run_dir_includes_tag(tmp_path):
    batch_run(_CIFS, config=_cfg(output_dir=str(tmp_path), run_tag="mytest"))
    subdirs = [d for d in tmp_path.iterdir() if d.is_dir()]
    assert "mytest" in subdirs[0].name


# ── Parallel mode ────────────────────────────────────────────────────────────

@pytest.mark.skipif(not _CIFS, reason="No CIF files in examples/")
def test_parallel_same_row_count_as_serial(tmp_path):
    substrates = ["Au_100", "Au_111"]
    df_serial = batch_run(_CIFS, config=_cfg(
        substrates=substrates, n_faces=1, n_jobs=1
    ))
    df_parallel = batch_run(_CIFS, config=_cfg(
        substrates=substrates, n_faces=1, n_jobs=2
    ))
    assert len(df_serial) == len(df_parallel)


# ── match_found column ───────────────────────────────────────────────────────

@pytest.mark.skipif(not _CIFS, reason="No CIF files in examples/")
def test_match_found_column_is_bool(tmp_path):
    df = batch_run(_CIFS, config=_cfg())
    assert df["match_found"].dtype == bool or df["match_found"].dtype == object
    # All values must be True or False
    assert set(df["match_found"].unique()).issubset({True, False})
