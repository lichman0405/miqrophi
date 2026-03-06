"""Tests for CIF parsing, systematic absences (centering + screw/glide)."""
import pytest
from pathlib import Path

from miqrocal.cif_parser import (
    _is_allowed,
    _lattice_type,
    _is_zone_allowed,
    bfdh_faces,
    best_surface_lattice,
    read_cell,
)

EXAMPLES = Path(__file__).parent.parent / "examples"
HKUST1   = EXAMPLES / "HKUST-1.cif"
CO_SQL   = next(EXAMPLES.glob("*sql*"), None)


# ── Centering extinction rules ───────────────────────────────────────────────

class TestCenteringP:
    def test_all_hkl_allowed(self):
        for h, k, l in [(1,0,0),(0,1,0),(0,0,1),(1,1,1),(2,3,5)]:
            assert _is_allowed(h, k, l, "P"), f"P should allow {h,k,l}"


class TestCenteringI:
    def test_h_plus_k_plus_l_odd_absent(self):
        assert not _is_allowed(1, 0, 0, "I")  # 1+0+0=1
        assert not _is_allowed(0, 1, 0, "I")  # 0+1+0=1
        assert not _is_allowed(1, 1, 1, "I")  # 1+1+1=3 (odd)
        assert not _is_allowed(2, 1, 0, "I")  # 3 odd

    def test_h_plus_k_plus_l_even_present(self):
        assert _is_allowed(1, 1, 0, "I")  # 2
        assert _is_allowed(2, 0, 0, "I")  # 2
        assert _is_allowed(2, 2, 2, "I")  # 6


class TestCenteringF:
    def test_mixed_parity_absent(self):
        assert not _is_allowed(1, 0, 0, "F")   # mixed
        assert not _is_allowed(0, 1, 0, "F")
        assert not _is_allowed(0, 0, 1, "F")
        assert not _is_allowed(1, 0, 1, "F")   # odd,even,odd → wait... no mixed parity means {odd,even} both present
        assert not _is_allowed(2, 1, 0, "F")   # even,odd,even → mixed

    def test_all_odd_present(self):
        assert _is_allowed(1, 1, 1, "F")
        assert _is_allowed(1, 3, 1, "F")
        assert _is_allowed(3, 3, 1, "F")

    def test_all_even_present(self):
        assert _is_allowed(2, 0, 0, "F")
        assert _is_allowed(2, 2, 0, "F")
        assert _is_allowed(2, 2, 2, "F")
        assert _is_allowed(0, 0, 2, "F")


class TestCenteringC:
    def test_h_plus_k_odd_absent(self):
        assert not _is_allowed(1, 0, 0, "C")
        assert not _is_allowed(0, 1, 0, "C")
        assert not _is_allowed(1, 2, 0, "C")

    def test_h_plus_k_even_present(self):
        assert _is_allowed(1, 1, 0, "C")
        assert _is_allowed(2, 0, 1, "C")


class TestCenteringA:
    def test_k_plus_l_odd_absent(self):
        assert not _is_allowed(0, 1, 0, "A")
        assert not _is_allowed(0, 0, 1, "A")

    def test_k_plus_l_even_present(self):
        assert _is_allowed(1, 1, 1, "A")   # k+l=2
        assert _is_allowed(0, 2, 0, "A")   # k+l=2


class TestCenteringR:
    def test_obverse_absent(self):
        # R obverse: present only when -h+k+l = 3n; else absent
        assert not _is_allowed(1, 1, 2, "R")   # -1+1+2=2 → not 3n → absent
        assert not _is_allowed(1, 2, 0, "R")   # -1+2+0=1 → not 3n → absent
        assert not _is_allowed(0, 0, 1, "R")   # -0+0+1=1 → not 3n → absent

    def test_obverse_present(self):
        assert _is_allowed(0, 0, 3, "R")   # -0+0+3=3  ✓
        assert _is_allowed(0, 0, 6, "R")   # -0+0+6=6  ✓
        assert _is_allowed(2, 1, 1, "R")   # -2+1+1=0  ✓
        assert _is_allowed(0, 1, 2, "R")   # -0+1+2=3  ✓


# ── Lattice type parsing ─────────────────────────────────────────────────────

class TestLatticeTypeParsing:
    def test_F_from_HM_symbol(self):
        cif = "_symmetry_space_group_name_H-M 'Fm-3m'"
        assert _lattice_type(cif) == "F"

    def test_P_from_HM_symbol(self):
        cif = "_symmetry_space_group_name_H-M 'P 21/c'"
        assert _lattice_type(cif) == "P"

    def test_I_from_HM_symbol(self):
        cif = "_symmetry_space_group_name_H-M 'Im-3m'"
        assert _lattice_type(cif) == "I"

    def test_F_from_SG_number_225(self):
        cif = "_symmetry_Int_Tables_number 225"
        assert _lattice_type(cif) == "F"

    def test_I_from_SG_number_229(self):
        cif = "_symmetry_Int_Tables_number 229"
        assert _lattice_type(cif) == "I"

    def test_default_P_when_absent(self):
        assert _lattice_type("_cell_length_a 10.0") == "P"

    def test_C_from_SG_number_15(self):
        cif = "_symmetry_Int_Tables_number 15"
        assert _lattice_type(cif) == "C"

    def test_R_from_SG_number_166(self):
        cif = "_symmetry_Int_Tables_number 166"
        assert _lattice_type(cif) == "R"


# ── Screw/glide zone conditions ──────────────────────────────────────────────

class TestZoneConditions:
    # SG 14: P 2_1/c — most common monoclinic for MOFs
    def test_P21c_0k0_k_odd_absent(self):
        assert not _is_zone_allowed(0, 1, 0, 14)   # k=1 odd → absent
        assert not _is_zone_allowed(0, 3, 0, 14)   # k=3 odd → absent

    def test_P21c_0k0_k_even_present(self):
        assert _is_zone_allowed(0, 2, 0, 14)
        assert _is_zone_allowed(0, 4, 0, 14)

    def test_P21c_h0l_l_odd_absent(self):
        assert not _is_zone_allowed(1, 0, 1, 14)   # l=1 odd → absent
        assert not _is_zone_allowed(2, 0, 1, 14)

    def test_P21c_h0l_l_even_present(self):
        assert _is_zone_allowed(1, 0, 2, 14)
        assert _is_zone_allowed(2, 0, 2, 14)

    def test_P21c_general_hkl_not_filtered(self):
        # Zone conditions only apply in special zones; general hkl always pass
        assert _is_zone_allowed(1, 1, 1, 14)
        assert _is_zone_allowed(2, 3, 1, 14)
        assert _is_zone_allowed(1, 2, 1, 14)

    # SG 19: P 2_1 2_1 2_1
    def test_P212121_h00_odd_absent(self):
        assert not _is_zone_allowed(1, 0, 0, 19)
        assert not _is_zone_allowed(3, 0, 0, 19)

    def test_P212121_0k0_odd_absent(self):
        assert not _is_zone_allowed(0, 1, 0, 19)

    def test_P212121_00l_odd_absent(self):
        assert not _is_zone_allowed(0, 0, 1, 19)

    def test_P212121_even_present(self):
        assert _is_zone_allowed(2, 0, 0, 19)
        assert _is_zone_allowed(0, 2, 0, 19)
        assert _is_zone_allowed(0, 0, 2, 19)

    # SG 62: P n m a
    def test_Pnma_0kl_kplusl_odd_absent(self):
        assert not _is_zone_allowed(0, 1, 0, 62)   # k+l=1 odd
        assert not _is_zone_allowed(0, 0, 1, 62)   # k+l=1 odd
        assert not _is_zone_allowed(0, 2, 1, 62)   # k+l=3 odd

    def test_Pnma_0kl_kplusl_even_present(self):
        assert _is_zone_allowed(0, 1, 1, 62)   # k+l=2 even
        assert _is_zone_allowed(0, 2, 0, 62)   # k+l=2 even

    # SG 227: F d-3m (e.g. diamond, some ZIF-like frameworks)
    def test_Fd3m_h00_h_not_4n_absent(self):
        assert not _is_zone_allowed(2, 0, 0, 227)   # h=2, not 4n → absent
        assert not _is_zone_allowed(1, 0, 0, 227)   # h=1, not 4n → absent

    def test_Fd3m_h00_h_4n_present(self):
        assert _is_zone_allowed(4, 0, 0, 227)
        assert _is_zone_allowed(8, 0, 0, 227)

    def test_Fd3m_0kl_kplusl_not_4n_absent(self):
        assert not _is_zone_allowed(0, 2, 0, 227)   # k+l=2, not 4n → absent
        assert not _is_zone_allowed(0, 1, 1, 227)   # k+l=2, not 4n → absent

    def test_Fd3m_0kl_kplusl_4n_present(self):
        assert _is_zone_allowed(0, 2, 2, 227)    # k+l=4 → present
        assert _is_zone_allowed(0, 0, 4, 227)    # k+l=4 → present

    # SG with no conditions at all (e.g. P-1, Fm-3m)
    def test_P1bar_no_conditions(self):
        # SG 2: P-1, no screw/glide
        assert _is_zone_allowed(1, 0, 0, 2)
        assert _is_zone_allowed(0, 1, 0, 2)
        assert _is_zone_allowed(0, 0, 1, 2)

    def test_Fm3m_no_extra_conditions(self):
        # SG 225: Fm-3m, F centering only, no screw/glide
        assert _is_zone_allowed(2, 0, 0, 225)
        assert _is_zone_allowed(1, 1, 1, 225)

    def test_unknown_sg_no_conditions(self):
        assert _is_zone_allowed(1, 0, 0, 9999)


# ── End-to-end BFDH faces ────────────────────────────────────────────────────

class TestBFDHFaces:
    @pytest.mark.skipif(not HKUST1.exists(), reason="HKUST-1.cif not in examples/")
    def test_hkust1_F_centering_applied(self):
        """All returned faces must have all-odd or all-even hkl (F-centering)."""
        faces = bfdh_faces(HKUST1, n_faces=5)
        assert len(faces) > 0
        for hkl, _ in faces:
            h, k, l = hkl
            parities = {abs(x) % 2 for x in (h, k, l)}
            assert len(parities) == 1, (
                f"Mixed-parity face {hkl} returned — F-centering not applied"
            )

    @pytest.mark.skipif(not HKUST1.exists(), reason="HKUST-1.cif not in examples/")
    def test_hkust1_no_100_face(self):
        """F-centering: (1,0,0) must be absent."""
        faces = bfdh_faces(HKUST1, n_faces=10)
        hkl_list = [hkl for hkl, _ in faces]
        assert (1, 0, 0) not in hkl_list
        assert (0, 1, 0) not in hkl_list
        assert (0, 0, 1) not in hkl_list

    @pytest.mark.skipif(CO_SQL is None, reason="Co-sql CIF not found")
    def test_co_sql_returns_faces(self):
        faces = bfdh_faces(CO_SQL, n_faces=3)
        assert 1 <= len(faces) <= 3

    @pytest.mark.skipif(not HKUST1.exists(), reason="HKUST-1.cif not in examples/")
    def test_hkust1_best_surface_lattice(self):
        results = best_surface_lattice(HKUST1, n_faces=3)
        assert len(results) > 0
        for hkl, d, lat in results:
            assert d > 0
            assert lat.a > 0
            assert lat.b > 0
            # Verify F-centering applied to Miller indices
            parities = {abs(x) % 2 for x in hkl}
            assert len(parities) == 1


# ── read_cell ────────────────────────────────────────────────────────────────

class TestReadCell:
    @pytest.mark.skipif(not HKUST1.exists(), reason="HKUST-1.cif not in examples/")
    def test_hkust1_cubic(self):
        cell = read_cell(HKUST1)
        assert abs(cell["a"] - 26.343) < 0.01
        assert abs(cell["b"] - 26.343) < 0.01
        assert abs(cell["c"] - 26.343) < 0.01
        assert abs(cell["alpha"] - 90.0) < 0.1
        assert abs(cell["gamma"] - 90.0) < 0.1

    @pytest.mark.skipif(not HKUST1.exists(), reason="HKUST-1.cif not in examples/")
    def test_hkust1_has_name(self):
        cell = read_cell(HKUST1)
        assert "name" in cell
        assert len(str(cell["name"])) > 0

    def test_missing_cif_raises(self):
        with pytest.raises(FileNotFoundError):
            read_cell("nonexistent.cif")
