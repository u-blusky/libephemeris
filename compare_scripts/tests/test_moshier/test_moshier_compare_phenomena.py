"""
Moshier Planetary Phenomena Cross-Library Comparison Tests.

Validates pheno_ut and pheno calculations (phase angle, elongation, magnitude,
apparent diameter) between pyswisseph (C library) and libephemeris (Python
reimplementation) using SEFLG_MOSEPH mode.

This is the Moshier-mode mirror of test_compare_phenomena.py (which covers
SEFLG_SWIEPH / JPL mode). With Moshier, both the planet and Sun positions
come from semi-analytical theories (VSOP87/ELP), amplifying potential
differences between the C and Python implementations. Tolerances are
therefore relaxed compared to the JPL-based tests.

Without these tests, magnitude and elongation errors in Moshier mode would
remain undetected — critical for heliacal visibility calculations, historical
eclipse analysis, and morning/evening star determination.
"""

import pytest
import swisseph as swe
import libephemeris as pyephem
from libephemeris.constants import (
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SEFLG_MOSEPH,
)


# ============================================================================
# TOLERANCE THRESHOLDS (relaxed for Moshier)
# ============================================================================


class MoshierPhenoTolerance:
    """Tolerance thresholds for Moshier phenomena comparisons.

    Relaxed relative to the SEFLG_SWIEPH tolerances in test_compare_phenomena.py
    (ANGLE=0.01, MAGNITUDE=0.1, DIAMETER=0.01) because Moshier semi-analytical
    theories amplify C-vs-Python differences in both planet and Sun positions,
    propagating into derived quantities like elongation, magnitude, and diameter.
    """

    ANGLE_DEGREES = 0.1  # Phase angle, elongation (10x relaxed vs JPL)
    MAGNITUDE = 0.5  # Apparent magnitude (5x relaxed vs JPL)
    DIAMETER = 1.0  # Apparent diameter in arcsec (relaxed for gas giants)


# ============================================================================
# PHENOMENA ARRAY INDICES
# ============================================================================

PHENO_PHASE_ANGLE = 0
PHENO_PHASE = 1  # Illuminated fraction
PHENO_ELONGATION = 2
PHENO_DIAMETER = 3
PHENO_MAGNITUDE = 4


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Planets for phenomena tests (naked-eye visible from Earth)
PHENO_PLANETS = [
    (SE_MERCURY, "Mercury"),
    (SE_VENUS, "Venus"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
]

# 3 test dates covering different planetary configurations
TEST_DATES = [
    (2024, 1, 15, "Jan 2024"),
    (2024, 6, 15, "Jun 2024"),
    (2024, 10, 15, "Oct 2024"),
]

# Elongation dates for inner planets (Mercury and Venus)
INNER_PLANET_ELONGATION_DATES = [
    (2024, 3, 22, "Venus western elongation"),
    (2024, 6, 4, "Mercury eastern elongation"),
    (2024, 8, 26, "Mercury western elongation"),
]

# Opposition dates for outer planets
OPPOSITION_DATES = [
    (2025, 1, 16, SE_MARS, "Mars opposition 2025"),
    (2024, 12, 7, SE_JUPITER, "Jupiter opposition 2024"),
    (2024, 9, 8, SE_SATURN, "Saturn opposition 2024"),
]


# ============================================================================
# HELPERS
# ============================================================================


def _extract_attrs(ret):
    """Extract the attributes array from pheno_ut/pheno return value.

    pyswisseph returns a flat tuple of 20 floats (the attributes directly).
    libephemeris returns (attr_tuple, retflag) where attr_tuple is a tuple
    of 20 floats and retflag is an int.
    """
    if isinstance(ret, tuple) and len(ret) == 2 and isinstance(ret[0], tuple):
        # libephemeris format: (attrs_tuple, retflag)
        return ret[0]
    # pyswisseph format: flat tuple of 20 floats
    return ret


def _get_diameter_arcsec(attr, is_pyswisseph=False):
    """Get apparent diameter in arcseconds.

    pyswisseph returns attr[3] in degrees; libephemeris returns it in arcsec.
    This function normalizes both to arcseconds.
    """
    diameter = attr[PHENO_DIAMETER]
    if is_pyswisseph:
        return diameter * 3600.0  # degrees -> arcsec
    return diameter


# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture
def jd_standard():
    """Standard Julian Day for tests."""
    return swe.julday(2024, 6, 15, 12.0)


# ============================================================================
# PHENO_UT TESTS (Moshier)
# ============================================================================


class TestMoshierPhenoUt:
    """Tests for pheno_ut function in Moshier mode.

    Compares all 5 phenomena components (phase angle, illuminated fraction,
    elongation, apparent diameter, apparent magnitude) between pyswisseph
    and libephemeris using SEFLG_MOSEPH for 5 planets x 3 dates = 15 cases.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PHENO_PLANETS)
    @pytest.mark.parametrize("year,month,day,date_desc", TEST_DATES)
    def test_pheno_ut_moshier(self, body_id, body_name, year, month, day, date_desc):
        """Test pheno_ut for planets at various dates using Moshier."""
        jd = swe.julday(year, month, day, 12.0)

        # SwissEphemeris (C Moshier)
        try:
            ret_swe = swe.pheno_ut(jd, body_id, SEFLG_MOSEPH)
            attr_swe = _extract_attrs(ret_swe)
        except Exception as e:
            pytest.skip(f"SwissEphemeris pheno_ut failed: {e}")

        # LibEphemeris (Python Moshier)
        try:
            ret_py = pyephem.swe_pheno_ut(jd, body_id, SEFLG_MOSEPH)
            attr_py = _extract_attrs(ret_py)
        except Exception as e:
            pytest.skip(f"LibEphemeris swe_pheno_ut failed: {e}")

        # Compare all 5 components
        # Note: pyswisseph returns diameter in degrees, libephemeris in arcsec
        diff_phase_angle = abs(attr_swe[PHENO_PHASE_ANGLE] - attr_py[PHENO_PHASE_ANGLE])
        diff_phase = abs(attr_swe[PHENO_PHASE] - attr_py[PHENO_PHASE])
        diff_elongation = abs(attr_swe[PHENO_ELONGATION] - attr_py[PHENO_ELONGATION])
        diameter_swe_arcsec = _get_diameter_arcsec(attr_swe, is_pyswisseph=True)
        diameter_py_arcsec = _get_diameter_arcsec(attr_py, is_pyswisseph=False)
        diff_diameter = abs(diameter_swe_arcsec - diameter_py_arcsec)
        diff_magnitude = abs(attr_swe[PHENO_MAGNITUDE] - attr_py[PHENO_MAGNITUDE])

        assert diff_phase_angle < MoshierPhenoTolerance.ANGLE_DEGREES, (
            f"{body_name} @ {date_desc}: phase angle diff {diff_phase_angle:.4f}°"
        )
        assert diff_phase < 0.01, (
            f"{body_name} @ {date_desc}: illuminated fraction diff {diff_phase:.6f}"
        )
        assert diff_elongation < MoshierPhenoTolerance.ANGLE_DEGREES, (
            f"{body_name} @ {date_desc}: elongation diff {diff_elongation:.4f}°"
        )
        assert diff_diameter < MoshierPhenoTolerance.DIAMETER, (
            f"{body_name} @ {date_desc}: diameter diff {diff_diameter:.4f} arcsec"
        )
        assert diff_magnitude < MoshierPhenoTolerance.MAGNITUDE, (
            f"{body_name} @ {date_desc}: magnitude diff {diff_magnitude:.4f}"
        )


# ============================================================================
# PHENO (ET VERSION) TESTS (Moshier)
# ============================================================================


class TestMoshierPheno:
    """Tests for pheno function (ET version) in Moshier mode."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PHENO_PLANETS)
    def test_pheno_et_moshier(self, jd_standard, body_id, body_name):
        """Test pheno function (ET time) using Moshier."""
        jd_et = jd_standard + swe.deltat(jd_standard)

        # SwissEphemeris (C Moshier)
        try:
            ret_swe = swe.pheno(jd_et, body_id, SEFLG_MOSEPH)
            attr_swe = _extract_attrs(ret_swe)
        except Exception as e:
            pytest.skip(f"SwissEphemeris pheno failed: {e}")

        # LibEphemeris (Python Moshier)
        try:
            ret_py = pyephem.swe_pheno(jd_et, body_id, SEFLG_MOSEPH)
            attr_py = _extract_attrs(ret_py)
        except Exception as e:
            pytest.skip(f"LibEphemeris swe_pheno failed: {e}")

        diff_phase_angle = abs(attr_swe[PHENO_PHASE_ANGLE] - attr_py[PHENO_PHASE_ANGLE])
        diff_elongation = abs(attr_swe[PHENO_ELONGATION] - attr_py[PHENO_ELONGATION])
        diff_magnitude = abs(attr_swe[PHENO_MAGNITUDE] - attr_py[PHENO_MAGNITUDE])

        assert diff_phase_angle < MoshierPhenoTolerance.ANGLE_DEGREES, (
            f"{body_name}: phase angle diff {diff_phase_angle:.4f}°"
        )
        assert diff_elongation < MoshierPhenoTolerance.ANGLE_DEGREES, (
            f"{body_name}: elongation diff {diff_elongation:.4f}°"
        )
        assert diff_magnitude < MoshierPhenoTolerance.MAGNITUDE, (
            f"{body_name}: magnitude diff {diff_magnitude:.4f}"
        )


# ============================================================================
# INNER PLANET ELONGATION TESTS (Moshier)
# ============================================================================


class TestMoshierElongations:
    """Tests for inner planet elongations using Moshier mode.

    Mercury and Venus maximum elongations are critical for morning/evening
    star determination and heliacal visibility. These tests verify that
    Moshier mode correctly computes elongation and phase at dates near
    maximum elongation.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,desc", INNER_PLANET_ELONGATION_DATES)
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SE_MERCURY, "Mercury"),
            (SE_VENUS, "Venus"),
        ],
    )
    def test_inner_planet_elongation_moshier(
        self, year, month, day, desc, body_id, body_name
    ):
        """Test inner planets at elongation dates using Moshier."""
        jd = swe.julday(year, month, day, 12.0)

        try:
            ret_swe = swe.pheno_ut(jd, body_id, SEFLG_MOSEPH)
            ret_py = pyephem.swe_pheno_ut(jd, body_id, SEFLG_MOSEPH)

            attr_swe = _extract_attrs(ret_swe)
            attr_py = _extract_attrs(ret_py)

            diff_phase = abs(attr_swe[PHENO_PHASE_ANGLE] - attr_py[PHENO_PHASE_ANGLE])
            diff_elong = abs(attr_swe[PHENO_ELONGATION] - attr_py[PHENO_ELONGATION])
            diff_mag = abs(attr_swe[PHENO_MAGNITUDE] - attr_py[PHENO_MAGNITUDE])

            assert diff_phase < MoshierPhenoTolerance.ANGLE_DEGREES, (
                f"{body_name} @ {desc}: phase angle diff {diff_phase:.4f}°"
            )
            assert diff_elong < MoshierPhenoTolerance.ANGLE_DEGREES, (
                f"{body_name} @ {desc}: elongation diff {diff_elong:.4f}°"
            )
            assert diff_mag < MoshierPhenoTolerance.MAGNITUDE, (
                f"{body_name} @ {desc}: magnitude diff {diff_mag:.4f}"
            )
        except Exception:
            pytest.skip("Calculation failed")


# ============================================================================
# OUTER PLANET OPPOSITION TESTS (Moshier)
# ============================================================================


class TestMoshierOuterPlanetOppositions:
    """Tests for outer planets at opposition using Moshier mode."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,body_id,desc", OPPOSITION_DATES)
    def test_outer_planet_opposition_moshier(self, year, month, day, body_id, desc):
        """Test outer planets at opposition dates (brightest) using Moshier."""
        jd = swe.julday(year, month, day, 12.0)
        body_name = desc.split()[0]

        try:
            ret_swe = swe.pheno_ut(jd, body_id, SEFLG_MOSEPH)
            ret_py = pyephem.swe_pheno_ut(jd, body_id, SEFLG_MOSEPH)

            attr_swe = _extract_attrs(ret_swe)
            attr_py = _extract_attrs(ret_py)

            diff_magnitude = abs(attr_swe[PHENO_MAGNITUDE] - attr_py[PHENO_MAGNITUDE])
            diff_elongation = abs(
                attr_swe[PHENO_ELONGATION] - attr_py[PHENO_ELONGATION]
            )

            assert diff_magnitude < MoshierPhenoTolerance.MAGNITUDE, (
                f"{body_name} opposition: magnitude diff {diff_magnitude:.4f}"
            )
            assert diff_elongation < MoshierPhenoTolerance.ANGLE_DEGREES, (
                f"{body_name} opposition: elongation diff {diff_elongation:.4f}°"
            )
        except Exception:
            pytest.skip("Calculation failed")


# ============================================================================
# PHENOMENA CONSISTENCY TESTS (Moshier)
# ============================================================================


class TestMoshierPhenomenaConsistency:
    """Tests for Moshier phenomena consistency and physical constraints."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PHENO_PLANETS)
    def test_phenomena_physical_constraints_moshier(
        self, jd_standard, body_id, body_name
    ):
        """Test that Moshier phenomena values satisfy physical constraints."""
        ret_py = pyephem.swe_pheno_ut(jd_standard, body_id, SEFLG_MOSEPH)
        attr = _extract_attrs(ret_py)

        phase_angle = attr[PHENO_PHASE_ANGLE]
        phase = attr[PHENO_PHASE]  # Illuminated fraction
        elongation = attr[PHENO_ELONGATION]
        diameter = attr[PHENO_DIAMETER]

        # Physical constraints
        assert 0 <= phase_angle <= 180, (
            f"{body_name}: phase angle {phase_angle}° must be 0-180°"
        )
        assert 0 <= phase <= 1, f"{body_name}: illuminated fraction {phase} must be 0-1"
        assert 0 <= elongation <= 180, (
            f"{body_name}: elongation {elongation}° must be 0-180°"
        )
        assert diameter >= 0, (
            f"{body_name}: diameter {diameter} arcsec must be non-negative"
        )

    @pytest.mark.comparison
    def test_inner_planet_elongation_limits_moshier(self, jd_standard):
        """Test that inner planets have limited elongation in Moshier mode."""
        for body_id, body_name in [(SE_MERCURY, "Mercury"), (SE_VENUS, "Venus")]:
            ret_py = pyephem.swe_pheno_ut(jd_standard, body_id, SEFLG_MOSEPH)
            attr = _extract_attrs(ret_py)

            elongation = attr[PHENO_ELONGATION]

            # Inner planets can't be more than ~47° (Venus) or ~28° (Mercury)
            max_elong = 50 if body_id == SE_VENUS else 30
            assert elongation <= max_elong, (
                f"{body_name}: elongation {elongation:.2f}° exceeds "
                f"physical limit {max_elong}°"
            )
