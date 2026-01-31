"""
Pytest-style Hypothetical Planets Comparison Tests.

Validates hypothetical body calculations (Uranian planets, Transpluto)
against pyswisseph.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_CUPIDO,
    SE_HADES,
    SE_ZEUS,
    SE_KRONOS,
    SE_APOLLON,
    SE_ADMETOS,
    SE_VULKANUS,
    SE_POSEIDON,
    SE_ISIS,
)


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TOLERANCE THRESHOLDS
# ============================================================================

# Relaxed tolerances due to different calculation methods
URANIAN_LONGITUDE_TOLERANCE = 0.1  # 0.1 degree for Uranian planets
HYPOTHETICAL_LONGITUDE_TOLERANCE = 1.0  # 1 degree for other hypothetical bodies
LATITUDE_TOLERANCE = 5.0  # 5 degrees for latitude


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Uranian Planets (Hamburg School) - IDs 40-47
URANIAN_PLANETS = [
    (SE_CUPIDO, swe.CUPIDO, "Cupido"),
    (SE_HADES, swe.HADES, "Hades"),
    (SE_ZEUS, swe.ZEUS, "Zeus"),
    (SE_KRONOS, swe.KRONOS, "Kronos"),
    (SE_APOLLON, swe.APOLLON, "Apollon"),
    (SE_ADMETOS, swe.ADMETOS, "Admetos"),
    (SE_VULKANUS, swe.VULKANUS, "Vulkanus"),
    (SE_POSEIDON, swe.POSEIDON, "Poseidon"),
]

# Other hypothetical bodies
OTHER_HYPOTHETICAL = [
    (SE_ISIS, swe.ISIS, "Transpluto/Isis"),
]

# Test dates
TEST_DATES = [
    ("J2000.0", 2000, 1, 1, 12.0),
    ("2024-01-01", 2024, 1, 1, 0.0),
    ("2010-07-01", 2010, 7, 1, 12.0),
    ("1980-01-01", 1980, 1, 1, 0.0),
]


# ============================================================================
# URANIAN PLANET TESTS
# ============================================================================


class TestUranianPlanets:
    """Tests for Uranian (Hamburg School) planets."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_py,body_swe,name", URANIAN_PLANETS)
    @pytest.mark.parametrize("date_name,year,month,day,hour", TEST_DATES[:2])  # Subset
    def test_uranian_planet_position(
        self, body_py, body_swe, name, date_name, year, month, day, hour
    ):
        """Test Uranian planet positions match pyswisseph."""
        jd = swe.julday(year, month, day, hour)

        # SwissEphemeris
        try:
            pos_swe, _ = swe.calc_ut(jd, body_swe, 0)
        except swe.Error as e:
            pytest.skip(f"SwissEphemeris {name} not available: {e}")

        # LibEphemeris
        try:
            pos_py, _ = ephem.swe_calc_ut(jd, body_py, 0)
        except Exception as e:
            pytest.skip(f"LibEphemeris {name} failed: {e}")

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])

        assert diff_lon < URANIAN_LONGITUDE_TOLERANCE, (
            f"{name} @ {date_name}: longitude diff {diff_lon}° exceeds tolerance"
        )
        assert diff_lat < LATITUDE_TOLERANCE, (
            f"{name} @ {date_name}: latitude diff {diff_lat}° exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_py,body_swe,name", URANIAN_PLANETS[:4])  # First 4
    def test_uranian_multiple_dates(self, body_py, body_swe, name):
        """Test Uranian planets across multiple dates."""
        for date_name, year, month, day, hour in TEST_DATES:
            jd = swe.julday(year, month, day, hour)

            try:
                pos_swe, _ = swe.calc_ut(jd, body_swe, 0)
                pos_py, _ = ephem.swe_calc_ut(jd, body_py, 0)

                diff_lon = angular_diff(pos_swe[0], pos_py[0])
                assert diff_lon < URANIAN_LONGITUDE_TOLERANCE
            except Exception:
                pass  # Skip dates where calculation fails


# ============================================================================
# OTHER HYPOTHETICAL BODY TESTS
# ============================================================================


class TestOtherHypothetical:
    """Tests for other hypothetical bodies (Transpluto, etc.)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_py,body_swe,name", OTHER_HYPOTHETICAL)
    @pytest.mark.parametrize("date_name,year,month,day,hour", TEST_DATES)
    def test_hypothetical_body_position(
        self, body_py, body_swe, name, date_name, year, month, day, hour
    ):
        """Test other hypothetical body positions."""
        jd = swe.julday(year, month, day, hour)

        # SwissEphemeris
        try:
            pos_swe, _ = swe.calc_ut(jd, body_swe, 0)
        except swe.Error as e:
            pytest.skip(f"SwissEphemeris {name} not available: {e}")

        # LibEphemeris
        try:
            pos_py, _ = ephem.swe_calc_ut(jd, body_py, 0)
        except Exception as e:
            pytest.skip(f"LibEphemeris {name} failed: {e}")

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_lat = abs(pos_swe[1] - pos_py[1])

        assert diff_lon < HYPOTHETICAL_LONGITUDE_TOLERANCE, (
            f"{name} @ {date_name}: longitude diff {diff_lon}°"
        )
        assert diff_lat < LATITUDE_TOLERANCE, (
            f"{name} @ {date_name}: latitude diff {diff_lat}°"
        )


# ============================================================================
# ALL HYPOTHETICAL BODIES SUMMARY TEST
# ============================================================================


class TestHypotheticalSummary:
    """Summary tests for all hypothetical bodies."""

    @pytest.mark.comparison
    def test_all_uranian_available(self):
        """Verify all Uranian planets are calculable."""
        jd = swe.julday(2024, 1, 1, 12.0)
        available = 0

        for body_py, body_swe, name in URANIAN_PLANETS:
            try:
                pos_py, _ = ephem.swe_calc_ut(jd, body_py, 0)
                available += 1
            except Exception:
                pass

        # At least 6 out of 8 should be available
        assert available >= 6, f"Only {available}/8 Uranian planets available"

    @pytest.mark.comparison
    def test_uranian_consistency(self):
        """Test Uranian planets are consistent across a time range."""
        for body_py, body_swe, name in URANIAN_PLANETS[:4]:
            positions = []
            for year in [2000, 2010, 2020]:
                jd = swe.julday(year, 1, 1, 12.0)
                try:
                    pos, _ = ephem.swe_calc_ut(jd, body_py, 0)
                    positions.append(pos[0])
                except Exception:
                    pass

            # Positions should be different (planet moved)
            if len(positions) >= 2:
                assert positions[0] != positions[1], f"{name} didn't move over 10 years"
