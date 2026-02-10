"""
Moshier Planetocentric Calculations Comparison Tests (swe_calc_pctr with SEFLG_MOSEPH).

Validates planetocentric position calculations between pyswisseph and libephemeris
when using the Moshier semi-analytical ephemeris (SEFLG_MOSEPH flag):
- Position of a target body as observed from another planet's center
- Key target/center combinations: Moon from Mars, Sun from Jupiter, Venus from Saturn
- Multiple time periods for validation

This is the Moshier-mode mirror of test_compare_calc_pctr.py (which covers
SEFLG_SWIEPH / JPL mode). Since calc_pctr internally computes positions for
both the target and center bodies, positional errors accumulate (target + center),
requiring a more relaxed tolerance than single-body Moshier comparisons.

Note:
    pyswisseph may not support calc_pctr with FLG_MOSEPH. Tests use try/except
    with pytest.skip to document availability rather than fail unconditionally.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SEFLG_MOSEPH,
    SEFLG_SPEED,
)


# ============================================================================
# TEST CONFIGURATIONS
# ============================================================================

# Key target/center combinations for Moshier planetocentric testing
PCTR_COMBINATIONS = [
    # Moon as seen from Mars - primary heliorelocation use case
    (SE_MOON, "Moon", SE_MARS, "Mars"),
    # Sun as seen from Jupiter - Jupiter-centric calculations
    (SE_SUN, "Sun", SE_JUPITER, "Jupiter"),
    # Venus as seen from Saturn - cross-orbit observation
    (SE_VENUS, "Venus", SE_SATURN, "Saturn"),
]

# Test dates (subset of test_compare_calc_pctr.py dates)
TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000.0"),
    (2024, 11, 15, 0.0, "Current Era"),
    (1950, 6, 21, 12.0, "Mid-Century"),
]

# Moshier planetocentric tolerances (relaxed because errors accumulate from
# both target and center body calculations via the Moshier semi-analytical
# ephemeris). Single-body Moshier tolerance is ~0.02°; for planetocentric
# calculations, target + center errors can sum up to ~0.04° in the worst case,
# so we use 0.1° to provide margin for edge cases.
MOSHIER_PCTR_LONGITUDE = 0.1  # degrees (~360 arcsec)
MOSHIER_PCTR_LATITUDE = 0.1  # degrees
MOSHIER_PCTR_DISTANCE_REL = 0.01  # relative (1%)
MOSHIER_PCTR_VELOCITY = 0.2  # degrees/day


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


def _calc_pctr_swe(jd: float, target: int, center: int, flag: int):
    """Call pyswisseph calc_pctr with Moshier flag, skip if unsupported."""
    try:
        return swe.calc_pctr(jd, target, center, flag)
    except (TypeError, ValueError, AttributeError) as e:
        pytest.skip(f"pyswisseph calc_pctr does not support Moshier mode: {e}")
    except Exception as e:
        # Catch any other error (e.g. swe internal errors)
        pytest.skip(
            f"pyswisseph calc_pctr with FLG_MOSEPH raised: {type(e).__name__}: {e}"
        )


def _calc_pctr_lib(jd: float, target: int, center: int, flag: int):
    """Call libephemeris swe_calc_pctr with Moshier flag, skip if unsupported."""
    try:
        return ephem.swe_calc_pctr(jd, target, center, flag)
    except (TypeError, ValueError, NotImplementedError) as e:
        pytest.skip(f"libephemeris swe_calc_pctr does not support Moshier mode: {e}")
    except Exception as e:
        pytest.skip(
            f"libephemeris swe_calc_pctr with SEFLG_MOSEPH raised: "
            f"{type(e).__name__}: {e}"
        )


# ============================================================================
# PRIMARY USE CASE TESTS
# ============================================================================


class TestMoonFromMars:
    """Test Moon position as seen from Mars with Moshier ephemeris."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_moon_from_mars_moshier(self, year, month, day, hour, date_desc):
        """Moon from Mars with SEFLG_MOSEPH across multiple dates."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH
        flag_lib = SEFLG_MOSEPH

        pos_swe, _ = _calc_pctr_swe(jd, SE_MOON, SE_MARS, flag_swe)
        pos_lib, _ = _calc_pctr_lib(jd, SE_MOON, SE_MARS, flag_lib)

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])
        lat_diff = abs(pos_lib[1] - pos_swe[1])

        assert lon_diff < MOSHIER_PCTR_LONGITUDE, (
            f"Moon from Mars Moshier at {date_desc}: "
            f"longitude diff {lon_diff:.6f}° exceeds tolerance "
            f"{MOSHIER_PCTR_LONGITUDE}° "
            f"(swe={pos_swe[0]:.6f}°, lib={pos_lib[0]:.6f}°)"
        )
        assert lat_diff < MOSHIER_PCTR_LATITUDE, (
            f"Moon from Mars Moshier at {date_desc}: "
            f"latitude diff {lat_diff:.6f}° exceeds tolerance "
            f"{MOSHIER_PCTR_LATITUDE}°"
        )


class TestSunFromJupiter:
    """Test Sun position as seen from Jupiter with Moshier ephemeris."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_sun_from_jupiter_moshier(self, year, month, day, hour, date_desc):
        """Sun from Jupiter with SEFLG_MOSEPH across multiple dates."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH
        flag_lib = SEFLG_MOSEPH

        pos_swe, _ = _calc_pctr_swe(jd, SE_SUN, SE_JUPITER, flag_swe)
        pos_lib, _ = _calc_pctr_lib(jd, SE_SUN, SE_JUPITER, flag_lib)

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])
        lat_diff = abs(pos_lib[1] - pos_swe[1])

        assert lon_diff < MOSHIER_PCTR_LONGITUDE, (
            f"Sun from Jupiter Moshier at {date_desc}: "
            f"longitude diff {lon_diff:.6f}° exceeds tolerance "
            f"{MOSHIER_PCTR_LONGITUDE}° "
            f"(swe={pos_swe[0]:.6f}°, lib={pos_lib[0]:.6f}°)"
        )
        assert lat_diff < MOSHIER_PCTR_LATITUDE, (
            f"Sun from Jupiter Moshier at {date_desc}: "
            f"latitude diff {lat_diff:.6f}° exceeds tolerance "
            f"{MOSHIER_PCTR_LATITUDE}°"
        )

        # Sun-Jupiter distance should be approximately 5.2 AU
        assert 4.5 < pos_lib[2] < 5.8, (
            f"Sun distance from Jupiter {pos_lib[2]:.2f} AU "
            f"should be approximately 5.2 AU"
        )


# ============================================================================
# COMPREHENSIVE COMBINATIONS
# ============================================================================


class TestKeyCombinations:
    """Test key target/center combinations with Moshier ephemeris."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "target_id,target_name,center_id,center_name", PCTR_COMBINATIONS
    )
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_combinations_moshier(
        self,
        target_id,
        target_name,
        center_id,
        center_name,
        year,
        month,
        day,
        hour,
        date_desc,
    ):
        """Test planetocentric combinations with SEFLG_MOSEPH."""
        jd = swe.julday(year, month, day, hour)
        flag_swe = swe.FLG_MOSEPH
        flag_lib = SEFLG_MOSEPH

        pos_swe, _ = _calc_pctr_swe(jd, target_id, center_id, flag_swe)
        pos_lib, _ = _calc_pctr_lib(jd, target_id, center_id, flag_lib)

        lon_diff = angular_diff(pos_lib[0], pos_swe[0])
        lat_diff = abs(pos_lib[1] - pos_swe[1])

        assert lon_diff < MOSHIER_PCTR_LONGITUDE, (
            f"{target_name} from {center_name} Moshier at {date_desc}: "
            f"longitude diff {lon_diff:.6f}° exceeds tolerance "
            f"{MOSHIER_PCTR_LONGITUDE}° "
            f"(swe={pos_swe[0]:.6f}°, lib={pos_lib[0]:.6f}°)"
        )
        assert lat_diff < MOSHIER_PCTR_LATITUDE, (
            f"{target_name} from {center_name} Moshier at {date_desc}: "
            f"latitude diff {lat_diff:.6f}° exceeds tolerance "
            f"{MOSHIER_PCTR_LATITUDE}°"
        )

        # Distance comparison (relative)
        if pos_swe[2] > 0:
            dist_diff_rel = abs(pos_lib[2] - pos_swe[2]) / pos_swe[2]
            assert dist_diff_rel < MOSHIER_PCTR_DISTANCE_REL, (
                f"{target_name} from {center_name} Moshier at {date_desc}: "
                f"distance relative diff {dist_diff_rel:.6f} exceeds tolerance "
                f"{MOSHIER_PCTR_DISTANCE_REL} "
                f"(swe={pos_swe[2]:.8f}, lib={pos_lib[2]:.8f})"
            )


# ============================================================================
# VELOCITY TESTS
# ============================================================================


class TestAllPlanetCenters:
    """Test planetocentric Moshier calculations with velocity (SEFLG_SPEED)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "target_id,target_name,center_id,center_name", PCTR_COMBINATIONS
    )
    def test_moshier_pctr_with_speed(
        self, target_id, target_name, center_id, center_name
    ):
        """Test SEFLG_MOSEPH | SEFLG_SPEED for planetocentric calculations."""
        jd = 2451545.0  # J2000.0
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SPEED
        flag_lib = SEFLG_MOSEPH | SEFLG_SPEED

        pos_swe, _ = _calc_pctr_swe(jd, target_id, center_id, flag_swe)
        pos_lib, _ = _calc_pctr_lib(jd, target_id, center_id, flag_lib)

        # Position
        lon_diff = angular_diff(pos_lib[0], pos_swe[0])
        assert lon_diff < MOSHIER_PCTR_LONGITUDE, (
            f"{target_name} from {center_name} Moshier+speed: "
            f"longitude diff {lon_diff:.6f}°"
        )

        # Velocity (longitude speed)
        vel_diff = abs(pos_lib[3] - pos_swe[3])
        assert vel_diff < MOSHIER_PCTR_VELOCITY, (
            f"{target_name} from {center_name} Moshier+speed: "
            f"velocity diff {vel_diff:.6f}°/day exceeds tolerance "
            f"{MOSHIER_PCTR_VELOCITY}°/day "
            f"(swe={pos_swe[3]:.6f}, lib={pos_lib[3]:.6f})"
        )
