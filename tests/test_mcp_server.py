"""Smoke tests for the MCP server tool functions (called directly, not via MCP)."""

import math
import os
import tempfile

import pytest

from miqrophi.mcp_server import analyze_epitaxy, list_substrates, parse_cif_surface

EXAMPLES = os.path.join(os.path.dirname(__file__), os.pardir, "examples")
HKUST1_CIF = os.path.join(EXAMPLES, "HKUST-1.cif")


# ---------------------------------------------------------------------------
# list_substrates
# ---------------------------------------------------------------------------

class TestListSubstrates:
    def test_returns_all_entries(self):
        result = list_substrates()
        assert result["count"] >= 10
        assert "Au_111" in result["substrates"]
        au = result["substrates"]["Au_111"]
        assert math.isclose(au["a"], 2.88, abs_tol=0.01)

    def test_each_entry_has_required_keys(self):
        result = list_substrates()
        for key, info in result["substrates"].items():
            assert "a" in info
            assert "b" in info
            assert "gamma_deg" in info
            assert "label" in info


# ---------------------------------------------------------------------------
# parse_cif_surface
# ---------------------------------------------------------------------------

class TestParseCifSurface:
    def test_bfdh_auto(self):
        result = parse_cif_surface(cif_path=HKUST1_CIF, n_faces=3)
        assert "cell_params" in result
        assert len(result["faces"]) == 3
        face = result["faces"][0]
        assert "hkl" in face
        assert "a" in face
        assert "d_hkl" in face

    def test_specific_hkl(self):
        result = parse_cif_surface(cif_path=HKUST1_CIF, hkl=[0, 0, 1])
        assert len(result["faces"]) == 1
        face = result["faces"][0]
        assert face["hkl"] == [0, 0, 1]
        # HKUST-1 cubic: (001) surface should have a ≈ b ≈ 26.3
        assert math.isclose(face["a"], 26.343, abs_tol=0.1)

    def test_cif_content_input(self):
        with open(HKUST1_CIF, "r") as f:
            content = f.read()
        result = parse_cif_surface(cif_content=content, hkl=[1, 1, 1])
        assert len(result["faces"]) == 1


# ---------------------------------------------------------------------------
# analyze_epitaxy
# ---------------------------------------------------------------------------

class TestAnalyzeEpitaxy:
    def test_with_substrate_name_and_lattice_params(self):
        result = analyze_epitaxy(
            substrate_name="Au_100",
            overlayer_a=18.62,
            overlayer_b=18.62,
            overlayer_gamma_deg=90.0,
            eta_tol=0.05,
            top_n=5,
        )
        assert result["level0"]["feasible"] is True
        assert result["level1"]["n_peaks"] >= 5
        assert result["level2"]["n_matches"] >= 1
        best = result["level2"]["matches"][0]
        assert best["eta"] < 0.05
        assert "M" in best
        assert "N" in best

    def test_with_cif_path(self):
        result = analyze_epitaxy(
            substrate_name="Au_111",
            cif_path=HKUST1_CIF,
            eta_tol=0.05,
            top_n=3,
        )
        assert "level0" in result
        assert "level1" in result
        assert "level2" in result
        assert "overlayer" in result

    def test_no_match_returns_empty_list(self):
        # Very tight tolerance — unlikely to find matches
        result = analyze_epitaxy(
            substrate_name="Au_111",
            overlayer_a=18.62,
            overlayer_b=18.62,
            overlayer_gamma_deg=90.0,
            eta_tol=0.0001,
            top_n=5,
        )
        assert result["level2"]["n_matches"] == 0
        assert result["level2"]["matches"] == []

    def test_output_dir_generates_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            result = analyze_epitaxy(
                substrate_name="Au_100",
                overlayer_a=18.62,
                overlayer_b=18.62,
                overlayer_gamma_deg=90.0,
                eta_tol=0.05,
                output_dir=tmpdir,
            )
            assert "files" in result
            assert os.path.isfile(result["files"]["png_path"])
            assert os.path.isfile(result["files"]["pdf_path"])

    def test_custom_substrate_params(self):
        result = analyze_epitaxy(
            substrate_a=4.08,
            substrate_b=4.08,
            substrate_gamma_deg=90.0,
            overlayer_a=18.62,
            overlayer_b=18.62,
            overlayer_gamma_deg=90.0,
            eta_tol=0.05,
        )
        assert result["substrate"]["label"] == "custom-substrate"
        assert result["level0"]["feasible"] is True

    def test_error_on_invalid_substrate(self):
        with pytest.raises(ValueError, match="Unknown substrate"):
            analyze_epitaxy(
                substrate_name="Unobtainium_111",
                overlayer_a=10.0,
                overlayer_b=10.0,
                overlayer_gamma_deg=90.0,
            )

    def test_error_on_missing_input(self):
        with pytest.raises(ValueError):
            analyze_epitaxy(substrate_name="Au_111")
