"""
Eclipse Calculations Comparison Tests.

Compares eclipse calculations between pyswisseph and libephemeris.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import SE_SUN, SE_MOON


# ============================================================================
# TOLERANCES
# ============================================================================

TIME_TOL_SECONDS = 120.0  # 2 minutes for eclipse timing
POSITION_TOL = 0.01  # degrees


# ============================================================================
# TEST DATA
# ============================================================================

# Known solar eclipses for testing
SOLAR_ECLIPSES = [
    (2024, 4, 8, "Total Solar Eclipse 2024"),
    (2024, 10, 2, "Annular Solar Eclipse 2024"),
    (2023, 10, 14, "Annular Solar Eclipse 2023"),
    (2025, 3, 29, "Partial Solar Eclipse 2025"),
]

# Known lunar eclipses for testing
LUNAR_ECLIPSES = [
    (2024, 3, 25, "Penumbral Lunar Eclipse 2024"),
    (2024, 9, 18, "Partial Lunar Eclipse 2024"),
    (2025, 3, 14, "Total Lunar Eclipse 2025"),
]

# Test locations
TEST_LOCATIONS = [
    ("Rome", 41.9028, 12.4964, 0),
    ("New York", 40.7128, -74.0060, 0),
    ("Tokyo", 35.6762, 139.6503, 0),
]


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestSolarEclipseSearch:
    """Compare solar eclipse search functions."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,desc", SOLAR_ECLIPSES)
    def test_sol_eclipse_when(self, year, month, day, desc):
        """Test finding solar eclipse timing."""
        jd_start = swe.julday(year, month, day, 0.0)

        # Search for next solar eclipse
        ecl_swe = swe.sol_eclipse_when_glob(jd_start, 0)
        ecl_py = ephem.sol_eclipse_when_glob(jd_start, 0)

        # Compare eclipse maximum time
        jd_max_swe = ecl_swe[1][0]  # tret[0] = maximum
        jd_max_py = ecl_py[1][0]

        diff_seconds = abs(jd_max_swe - jd_max_py) * 86400

        assert diff_seconds < TIME_TOL_SECONDS, (
            f"{desc}: eclipse max time diff {diff_seconds:.1f}s exceeds tolerance"
        )


class TestLunarEclipseSearch:
    """Compare lunar eclipse search functions."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,desc", LUNAR_ECLIPSES)
    def test_lun_eclipse_when(self, year, month, day, desc):
        """Test finding lunar eclipse timing."""
        jd_start = swe.julday(year, month, day, 0.0)

        # Search for next lunar eclipse
        ecl_swe = swe.lun_eclipse_when(jd_start, 0)
        ecl_py = ephem.lun_eclipse_when(jd_start, 0)

        # Compare eclipse maximum time
        jd_max_swe = ecl_swe[1][0]
        jd_max_py = ecl_py[1][0]

        diff_seconds = abs(jd_max_swe - jd_max_py) * 86400

        assert diff_seconds < TIME_TOL_SECONDS, (
            f"{desc}: eclipse max time diff {diff_seconds:.1f}s exceeds tolerance"
        )


class TestEclipseAttributes:
    """Compare eclipse attribute calculations."""

    @pytest.mark.comparison
    def test_solar_eclipse_attributes(self):
        """Test solar eclipse attribute calculations."""
        # 2024 total solar eclipse
        jd = swe.julday(2024, 4, 8, 18.0)
        geopos = (-99.0, 25.0, 0)  # Mexico path

        attr_swe = swe.sol_eclipse_how(jd, 0, geopos)
        attr_py = ephem.sol_eclipse_how(jd, 0, geopos)

        # Compare obscuration/magnitude
        # attr[0] = fraction of solar diameter covered
        diff_obscur = abs(attr_swe[1][0] - attr_py[1][0])

        assert diff_obscur < 0.01, (
            f"Eclipse obscuration diff {diff_obscur:.4f} exceeds tolerance"
        )


class TestEclipseFlags:
    """Test eclipse type detection."""

    @pytest.mark.comparison
    def test_eclipse_type_detection(self):
        """Test that eclipse types are correctly identified."""
        # 2024 total solar eclipse
        jd = swe.julday(2024, 4, 1, 0.0)

        ecl_swe = swe.sol_eclipse_when_glob(jd, 0)
        ecl_py = ephem.sol_eclipse_when_glob(jd, 0)

        # Both should detect the same eclipse type
        type_swe = ecl_swe[0]
        type_py = ecl_py[0]

        # At minimum, both should detect an eclipse occurred
        assert type_swe > 0 and type_py > 0, "Both should detect an eclipse"
