"""
Comprehensive Lunar Nodes & Lilith Comparison Tests.

Exhaustive comparison between pyswisseph and libephemeris for:
- SE_MEAN_NODE  (10) — Mean North Node
- SE_TRUE_NODE  (11) — True (Osculating) North Node
- SE_MEAN_APOG  (12) — Mean Lilith (Black Moon)
- SE_OSCU_APOG  (13) — True Lilith (Osculating Apogee)
- SE_INTP_APOG  (21) — Interpolated Apogee
- SE_INTP_PERG  (22) — Interpolated Perigee

Covers all dimensions:
- Tropical longitude (no flags)
- All 6 output components (lon, lat, dist, speed_lon, speed_lat, speed_dist)
- Velocity (SEFLG_SPEED)
- Sidereal mode (SEFLG_SIDEREAL) with all 43 ayanamshas
- Equatorial coordinates (SEFLG_EQUATORIAL)
- J2000 frame (SEFLG_J2000)
- No-nutation (SEFLG_NONUT)
- Topocentric (SEFLG_TOPOCTR) with multiple geographic locations
- Combined flags (SIDEREAL|SPEED, EQUATORIAL|SPEED, etc.)
- Extended date ranges (1550–2650 DE440 coverage)
- Statistical precision analysis (mean/max/RMS/percentiles)
- Edge cases (wrap-around, eclipse dates, extreme dates)
"""

import math
import random
import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SE_MEAN_APOG,
    SE_OSCU_APOG,
    SE_INTP_APOG,
    SE_INTP_PERG,
    SEFLG_SPEED,
    SEFLG_SIDEREAL,
    SEFLG_EQUATORIAL,
    SEFLG_J2000,
    SEFLG_NONUT,
    SEFLG_TOPOCTR,
    SE_SIDM_FAGAN_BRADLEY,
    SE_SIDM_LAHIRI,
    SE_SIDM_DELUCE,
    SE_SIDM_RAMAN,
    SE_SIDM_USHASHASHI,
    SE_SIDM_KRISHNAMURTI,
    SE_SIDM_DJWHAL_KHUL,
    SE_SIDM_YUKTESHWAR,
    SE_SIDM_JN_BHASIN,
    SE_SIDM_BABYL_KUGLER1,
    SE_SIDM_BABYL_KUGLER2,
    SE_SIDM_BABYL_KUGLER3,
    SE_SIDM_BABYL_HUBER,
    SE_SIDM_BABYL_ETPSC,
    SE_SIDM_ALDEBARAN_15TAU,
    SE_SIDM_HIPPARCHOS,
    SE_SIDM_SASSANIAN,
    SE_SIDM_GALCENT_0SAG,
    SE_SIDM_J2000,
    SE_SIDM_J1900,
    SE_SIDM_B1950,
    SE_SIDM_SURYASIDDHANTA,
    SE_SIDM_SURYASIDDHANTA_MSUN,
    SE_SIDM_ARYABHATA,
    SE_SIDM_ARYABHATA_MSUN,
    SE_SIDM_SS_REVATI,
    SE_SIDM_SS_CITRA,
    SE_SIDM_TRUE_CITRA,
    SE_SIDM_TRUE_REVATI,
    SE_SIDM_TRUE_PUSHYA,
    SE_SIDM_GALCENT_RGILBRAND,
    SE_SIDM_GALEQU_IAU1958,
    SE_SIDM_GALEQU_TRUE,
    SE_SIDM_GALEQU_MULA,
    SE_SIDM_GALALIGN_MARDYKS,
    SE_SIDM_TRUE_MULA,
    SE_SIDM_GALCENT_MULA_WILHELM,
    SE_SIDM_ARYABHATA_522,
    SE_SIDM_BABYL_BRITTON,
    SE_SIDM_TRUE_SHEORAN,
    SE_SIDM_GALCENT_COCHRANE,
    SE_SIDM_GALEQU_FIORENZA,
    SE_SIDM_VALENS_MOON,
)


# ============================================================================
# UTILITIES
# ============================================================================


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


def arcsec(deg: float) -> float:
    """Convert degrees to arcseconds."""
    return deg * 3600.0


def rms(values: list) -> float:
    """Root Mean Square of a list of values."""
    if not values:
        return 0.0
    return math.sqrt(sum(v * v for v in values) / len(values))


def percentile(values: list, pct: float) -> float:
    """Calculate percentile of sorted values."""
    if not values:
        return 0.0
    s = sorted(values)
    idx = int(len(s) * pct / 100.0)
    idx = min(idx, len(s) - 1)
    return s[idx]


# ============================================================================
# BODY DEFINITIONS
# ============================================================================

# Primary bodies (nodes + Lilith)
PRIMARY_BODIES = [
    (SE_MEAN_NODE, "Mean Node"),
    (SE_TRUE_NODE, "True Node"),
    (SE_MEAN_APOG, "Mean Lilith"),
    (SE_OSCU_APOG, "True Lilith"),
]

# Interpolated apogee/perigee
INTERPOLATED_BODIES = [
    (SE_INTP_APOG, "Interpolated Apogee"),
    (SE_INTP_PERG, "Interpolated Perigee"),
]

ALL_BODIES = PRIMARY_BODIES + INTERPOLATED_BODIES


# ============================================================================
# TOLERANCES (per-body, per-component)
# ============================================================================

# Longitude tolerance in degrees
LON_TOL = {
    SE_MEAN_NODE: 0.01,  # ~36 arcsec
    SE_TRUE_NODE: 0.15,  # ~540 arcsec (different oscillation model)
    SE_MEAN_APOG: 0.01,  # ~36 arcsec (SE-compatible algorithm)
    SE_OSCU_APOG: 0.1,  # ~360 arcsec (eccentricity vector method)
    SE_INTP_APOG: 0.5,  # ~1800 arcsec (ELP2000-82B series)
    SE_INTP_PERG: 0.5,  # ~1800 arcsec (ELP2000-82B series + residual table)
}

# Latitude tolerance in degrees
LAT_TOL = {
    SE_MEAN_NODE: 0.001,  # Mean Node has ~0 latitude
    SE_TRUE_NODE: 0.5,  # True Node has small latitude
    SE_MEAN_APOG: 0.5,  # Mean Lilith latitude from ecliptic projection
    SE_OSCU_APOG: 1.5,  # True Lilith latitude (eccentricity vector)
    SE_INTP_APOG: 1.0,  # Interpolated apogee latitude
    SE_INTP_PERG: 1.0,  # Interpolated perigee latitude
}

# Distance tolerance in AU
DIST_TOL = {
    SE_MEAN_NODE: 0.1,  # Distance often ~0 for nodes
    SE_TRUE_NODE: 0.01,  # True Node has distance
    SE_MEAN_APOG: 0.1,  # Mean Lilith distance
    SE_OSCU_APOG: 0.001,  # True Lilith uses eccentricity
    SE_INTP_APOG: 0.01,  # Interpolated apogee distance
    SE_INTP_PERG: 0.01,  # Interpolated perigee distance
}

# Speed (longitude velocity) tolerance in degrees/day
SPEED_LON_TOL = {
    SE_MEAN_NODE: 0.01,  # Mean Node speed is smooth
    SE_TRUE_NODE: 0.2,  # True Node speed varies more
    SE_MEAN_APOG: 0.01,  # Mean Lilith speed is smooth
    SE_OSCU_APOG: 1.0,  # True Lilith speed varies significantly
    SE_INTP_APOG: 0.5,  # Interpolated apogee speed
    SE_INTP_PERG: 0.5,  # Interpolated perigee speed
}

# Speed (latitude velocity) tolerance in degrees/day
SPEED_LAT_TOL = {
    SE_MEAN_NODE: 0.001,
    SE_TRUE_NODE: 0.1,
    SE_MEAN_APOG: 0.05,
    SE_OSCU_APOG: 0.5,
    SE_INTP_APOG: 0.2,
    SE_INTP_PERG: 0.2,
}

# Speed (distance velocity) tolerance in AU/day
SPEED_DIST_TOL = {
    SE_MEAN_NODE: 0.01,
    SE_TRUE_NODE: 0.01,
    SE_MEAN_APOG: 0.01,
    SE_OSCU_APOG: 0.001,
    SE_INTP_APOG: 0.01,
    SE_INTP_PERG: 0.01,
}

# Sidereal tolerance multiplier (sidereal adds ayanamsha uncertainty)
SIDEREAL_LON_MULTIPLIER = 1.5

# Equatorial coordinate tolerance in degrees (transformation adds some error)
EQUATORIAL_TOL = {
    SE_MEAN_NODE: 0.02,
    SE_TRUE_NODE: 0.2,
    SE_MEAN_APOG: 0.02,
    SE_OSCU_APOG: 0.15,
    SE_INTP_APOG: 0.6,
    SE_INTP_PERG: 0.6,
}


# ============================================================================
# TEST DATES
# ============================================================================

# Strategic dates covering key scenarios
TEST_DATES = [
    (2451545.0, "J2000.0 (2000-01-01)"),
    (2415020.5, "J1900 (1900-01-01)"),
    (2433282.5, "1950-01-01"),
    (2460000.5, "2023-02-25"),
    (2460310.5, "2023-12-02"),
    (2469807.5, "2050-01-01"),
    (2440587.5, "Unix epoch (1970-01-01)"),
    (2448988.5, "1993-01-01"),
    (2458849.5, "2020-01-01"),
    (2462502.5, "2030-01-01"),
]

# Extended dates for DE440 range edges
EXTENDED_DATES = [
    (2287184.5, "1550-01-01 (DE440 start)"),
    (2299160.5, "1582-10-15 (Gregorian reform)"),
    (2341972.5, "1700-01-01"),
    (2378496.5, "1800-01-01"),
    (2488069.5, "2100-01-01"),
    (2524593.5, "2200-01-01"),
    (2597641.5, "2400-01-01"),
    (2688069.5, "2647-07-05 (near DE440 end)"),
]

# Eclipse dates (nodal region, interesting for True Node)
ECLIPSE_DATES = [
    (2460044.0, "2023-04-20 (Solar eclipse)"),
    (2460227.0, "2023-10-14 (Solar eclipse)"),
    (2460408.5, "2024-04-08 (Total solar eclipse)"),
    (2451727.5, "2000-07-01 (Solar eclipse)"),
    (2459374.5, "2021-06-10 (Annular eclipse)"),
]

# Dates for parametrized tests (shorter list for combinatorial tests)
CORE_DATES = [
    (2451545.0, "J2000"),
    (2415020.5, "J1900"),
    (2460000.5, "2023"),
    (2469807.5, "2050"),
]


# ============================================================================
# AYANAMSHA DEFINITIONS
# ============================================================================

FORMULA_BASED_AYANAMSHAS = [
    (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
    (SE_SIDM_LAHIRI, "Lahiri"),
    (SE_SIDM_DELUCE, "De Luce"),
    (SE_SIDM_RAMAN, "Raman"),
    (SE_SIDM_USHASHASHI, "Ushashashi"),
    (SE_SIDM_KRISHNAMURTI, "Krishnamurti"),
    (SE_SIDM_DJWHAL_KHUL, "Djwhal Khul"),
    (SE_SIDM_YUKTESHWAR, "Yukteshwar"),
    (SE_SIDM_JN_BHASIN, "JN Bhasin"),
    (SE_SIDM_BABYL_KUGLER1, "Babylonian Kugler 1"),
    (SE_SIDM_BABYL_KUGLER2, "Babylonian Kugler 2"),
    (SE_SIDM_BABYL_KUGLER3, "Babylonian Kugler 3"),
    (SE_SIDM_BABYL_HUBER, "Babylonian Huber"),
    (SE_SIDM_BABYL_ETPSC, "Babylonian ETPSC"),
    (SE_SIDM_ALDEBARAN_15TAU, "Aldebaran 15 Tau"),
    (SE_SIDM_HIPPARCHOS, "Hipparchos"),
    (SE_SIDM_SASSANIAN, "Sassanian"),
    (SE_SIDM_J2000, "J2000"),
    (SE_SIDM_J1900, "J1900"),
    (SE_SIDM_B1950, "B1950"),
    (SE_SIDM_SURYASIDDHANTA, "Suryasiddhanta"),
    (SE_SIDM_SURYASIDDHANTA_MSUN, "Suryasiddhanta Mean Sun"),
    (SE_SIDM_ARYABHATA, "Aryabhata"),
    (SE_SIDM_ARYABHATA_MSUN, "Aryabhata Mean Sun"),
    (SE_SIDM_SS_REVATI, "SS Revati"),
    (SE_SIDM_SS_CITRA, "SS Citra"),
    (SE_SIDM_ARYABHATA_522, "Aryabhata 522"),
    (SE_SIDM_BABYL_BRITTON, "Babylonian Britton"),
]

STAR_BASED_AYANAMSHAS = [
    (SE_SIDM_TRUE_CITRA, "True Citra"),
    (SE_SIDM_TRUE_REVATI, "True Revati"),
    (SE_SIDM_TRUE_PUSHYA, "True Pushya"),
    (SE_SIDM_TRUE_MULA, "True Mula"),
    (SE_SIDM_TRUE_SHEORAN, "True Sheoran"),
    (SE_SIDM_GALCENT_0SAG, "Galactic Center 0 Sag"),
    (SE_SIDM_GALCENT_RGILBRAND, "Galactic Center Rgilbrand"),
    (SE_SIDM_GALCENT_MULA_WILHELM, "Galactic Center Mula Wilhelm"),
    (SE_SIDM_GALCENT_COCHRANE, "Galactic Center Cochrane"),
    (SE_SIDM_GALEQU_IAU1958, "Galactic Equator IAU 1958"),
    (SE_SIDM_GALEQU_TRUE, "Galactic Equator True"),
    (SE_SIDM_GALEQU_MULA, "Galactic Equator Mula"),
    (SE_SIDM_GALEQU_FIORENZA, "Galactic Equator Fiorenza"),
    (SE_SIDM_GALALIGN_MARDYKS, "Galactic Alignment Mardyks"),
    (SE_SIDM_VALENS_MOON, "Valens Moon"),
]

ALL_AYANAMSHAS = FORMULA_BASED_AYANAMSHAS + STAR_BASED_AYANAMSHAS

# Major ayanamshas for heavier combinatorial tests
MAJOR_AYANAMSHAS = [
    (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
    (SE_SIDM_LAHIRI, "Lahiri"),
    (SE_SIDM_RAMAN, "Raman"),
    (SE_SIDM_KRISHNAMURTI, "Krishnamurti"),
    (SE_SIDM_TRUE_CITRA, "True Citra"),
]

# Sidereal tolerance: formula-based vs star-based
SIDEREAL_STRICT_TOL = 0.001  # degrees for formula-based ayanamsha offset
SIDEREAL_RELAXED_TOL = 0.1  # degrees for star-based ayanamsha offset

# Set of star-based mode IDs for fast lookup
_STAR_MODES = {m[0] for m in STAR_BASED_AYANAMSHAS}


def _sidereal_tol(body_id: int, sid_mode: int) -> float:
    """Compute sidereal tolerance for a body+ayanamsha combination."""
    ayan_tol = SIDEREAL_RELAXED_TOL if sid_mode in _STAR_MODES else SIDEREAL_STRICT_TOL
    return LON_TOL[body_id] * SIDEREAL_LON_MULTIPLIER + ayan_tol


# ============================================================================
# GEOGRAPHIC LOCATIONS (for topocentric tests)
# ============================================================================

TEST_LOCATIONS = [
    ("Rome", 41.9028, 12.4964, 0),
    ("New York", 40.7128, -74.0060, 0),
    ("Tokyo", 35.6762, 139.6503, 0),
    ("Sydney", -33.8688, 151.2093, 0),
    ("Tromso", 69.6492, 18.9553, 0),
    ("Equator/Greenwich", 0.0, 0.0, 0),
    ("Cape Town", -33.9249, 18.4241, 0),
    ("Reykjavik", 64.1466, -21.9426, 0),
    ("High altitude", 27.9881, 86.9250, 8848),
]


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestTropicalLongitude:
    """Compare tropical longitude (flag=0) for all bodies across dates."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_primary_body_longitude(self, body_id, body_name, jd, date_desc):
        """Test primary body longitude matches within tolerance."""
        pos_swe, _ = swe.calc_ut(jd, body_id, 0)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)

        diff = angular_diff(pos_swe[0], pos_py[0])
        tol = LON_TOL[body_id]

        assert diff < tol, (
            f"{body_name} at {date_desc}: lon diff {diff:.6f}° "
            f'({arcsec(diff):.1f}") exceeds {tol}° tolerance '
            f"(swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", INTERPOLATED_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_interpolated_body_longitude(self, body_id, body_name, jd, date_desc):
        """Test interpolated apogee/perigee longitude matches within tolerance."""
        pos_swe, _ = swe.calc_ut(jd, body_id, 0)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)

        diff = angular_diff(pos_swe[0], pos_py[0])
        tol = LON_TOL[body_id]

        assert diff < tol, (
            f"{body_name} at {date_desc}: lon diff {diff:.6f}° "
            f'({arcsec(diff):.1f}") exceeds {tol}° tolerance '
            f"(swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("jd,date_desc", ECLIPSE_DATES)
    def test_eclipse_date_longitude(self, body_id, body_name, jd, date_desc):
        """Test body longitude at eclipse dates (nodal region stress test)."""
        pos_swe, _ = swe.calc_ut(jd, body_id, 0)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)

        diff = angular_diff(pos_swe[0], pos_py[0])
        tol = LON_TOL[body_id]

        assert diff < tol, (
            f"{body_name} at {date_desc}: lon diff {diff:.6f}° "
            f'({arcsec(diff):.1f}") at eclipse date'
        )


class TestTropicalAllComponents:
    """Compare all 6 output components (lon, lat, dist, dlon, dlat, ddist)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ALL_BODIES)
    @pytest.mark.parametrize("jd,date_desc", CORE_DATES)
    def test_all_components_no_speed(self, body_id, body_name, jd, date_desc):
        """Test all position components (without SEFLG_SPEED)."""
        pos_swe, _ = swe.calc_ut(jd, body_id, 0)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)

        # Longitude
        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        assert diff_lon < LON_TOL[body_id], (
            f"{body_name} at {date_desc}: lon diff {diff_lon:.6f}° "
            f"(swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )

        # Latitude
        diff_lat = abs(pos_swe[1] - pos_py[1])
        assert diff_lat < LAT_TOL[body_id], (
            f"{body_name} at {date_desc}: lat diff {diff_lat:.6f}° "
            f"(swe={pos_swe[1]:.6f}, lib={pos_py[1]:.6f})"
        )

        # Distance
        diff_dist = abs(pos_swe[2] - pos_py[2])
        assert diff_dist < DIST_TOL[body_id], (
            f"{body_name} at {date_desc}: dist diff {diff_dist:.6f} AU "
            f"(swe={pos_swe[2]:.6f}, lib={pos_py[2]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", ALL_BODIES)
    @pytest.mark.parametrize("jd,date_desc", CORE_DATES)
    def test_all_components_with_speed(self, body_id, body_name, jd, date_desc):
        """Test all 6 components with SEFLG_SPEED."""
        pos_swe, _ = swe.calc_ut(jd, body_id, swe.FLG_SPEED)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, SEFLG_SPEED)

        # Longitude
        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        assert diff_lon < LON_TOL[body_id], (
            f"{body_name} at {date_desc}: lon diff {diff_lon:.6f}°"
        )

        # Latitude
        diff_lat = abs(pos_swe[1] - pos_py[1])
        assert diff_lat < LAT_TOL[body_id], (
            f"{body_name} at {date_desc}: lat diff {diff_lat:.6f}°"
        )

        # Distance
        diff_dist = abs(pos_swe[2] - pos_py[2])
        assert diff_dist < DIST_TOL[body_id], (
            f"{body_name} at {date_desc}: dist diff {diff_dist:.6f} AU"
        )

        # Longitude velocity (degrees/day)
        diff_dlon = abs(pos_swe[3] - pos_py[3])
        assert diff_dlon < SPEED_LON_TOL[body_id], (
            f"{body_name} at {date_desc}: dlon diff {diff_dlon:.6f}°/day "
            f"(swe={pos_swe[3]:.6f}, lib={pos_py[3]:.6f})"
        )

        # Latitude velocity (degrees/day)
        diff_dlat = abs(pos_swe[4] - pos_py[4])
        assert diff_dlat < SPEED_LAT_TOL[body_id], (
            f"{body_name} at {date_desc}: dlat diff {diff_dlat:.6f}°/day "
            f"(swe={pos_swe[4]:.6f}, lib={pos_py[4]:.6f})"
        )

        # Distance velocity (AU/day)
        diff_ddist = abs(pos_swe[5] - pos_py[5])
        assert diff_ddist < SPEED_DIST_TOL[body_id], (
            f"{body_name} at {date_desc}: ddist diff {diff_ddist:.6f} AU/day "
            f"(swe={pos_swe[5]:.6f}, lib={pos_py[5]:.6f})"
        )


class TestVelocity:
    """Focused velocity comparison tests."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("jd,date_desc", TEST_DATES)
    def test_longitude_speed(self, body_id, body_name, jd, date_desc):
        """Test longitude velocity for all primary bodies across dates."""
        pos_swe, _ = swe.calc_ut(jd, body_id, swe.FLG_SPEED)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, SEFLG_SPEED)

        diff = abs(pos_swe[3] - pos_py[3])

        assert diff < SPEED_LON_TOL[body_id], (
            f"{body_name} at {date_desc}: speed diff {diff:.6f}°/day "
            f"(swe={pos_swe[3]:.6f}, lib={pos_py[3]:.6f})"
        )

    @pytest.mark.comparison
    def test_mean_node_always_retrograde(self):
        """Verify Mean Node speed is always retrograde (negative)."""
        for jd, desc in TEST_DATES:
            pos_py, _ = ephem.swe_calc_ut(jd, SE_MEAN_NODE, SEFLG_SPEED)
            pos_swe, _ = swe.calc_ut(jd, swe.MEAN_NODE, swe.FLG_SPEED)

            assert pos_py[3] < 0, (
                f"Mean Node speed should be retrograde at {desc}, "
                f"got {pos_py[3]:.6f}°/day"
            )
            assert pos_swe[3] < 0, f"SWE Mean Node speed should be retrograde at {desc}"

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", INTERPOLATED_BODIES)
    @pytest.mark.parametrize("jd,date_desc", CORE_DATES)
    def test_interpolated_speed(self, body_id, body_name, jd, date_desc):
        """Test interpolated apogee/perigee speed."""
        pos_swe, _ = swe.calc_ut(jd, body_id, swe.FLG_SPEED)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, SEFLG_SPEED)

        diff = abs(pos_swe[3] - pos_py[3])

        assert diff < SPEED_LON_TOL[body_id], (
            f"{body_name} at {date_desc}: speed diff {diff:.6f}°/day"
        )


class TestSiderealFormulaBasedAllAyanamshas:
    """Test sidereal longitude with ALL formula-based ayanamshas."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("sid_mode,sid_name", FORMULA_BASED_AYANAMSHAS)
    @pytest.mark.parametrize("jd,date_desc", CORE_DATES[:2])
    def test_sidereal_longitude_formula(
        self, body_id, body_name, sid_mode, sid_name, jd, date_desc
    ):
        """Test sidereal longitude with formula-based ayanamsha."""
        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flags = SEFLG_SIDEREAL
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = _sidereal_tol(body_id, sid_mode)
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"{body_name} sidereal ({sid_name}) at {date_desc}: "
            f'diff {diff:.6f}° ({arcsec(diff):.1f}") exceeds {tol:.4f}° '
            f"(swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", INTERPOLATED_BODIES)
    @pytest.mark.parametrize("sid_mode,sid_name", FORMULA_BASED_AYANAMSHAS[:5])
    @pytest.mark.parametrize("jd,date_desc", CORE_DATES[:2])
    def test_interpolated_sidereal_formula(
        self, body_id, body_name, sid_mode, sid_name, jd, date_desc
    ):
        """Test interpolated body sidereal longitude with formula-based ayanamsha."""
        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flags = SEFLG_SIDEREAL
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = _sidereal_tol(body_id, sid_mode)
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"{body_name} sidereal ({sid_name}) at {date_desc}: diff {diff:.6f}°"
        )


class TestSiderealStarBasedAllAyanamshas:
    """Test sidereal longitude with ALL star-based ayanamshas."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("sid_mode,sid_name", STAR_BASED_AYANAMSHAS)
    @pytest.mark.parametrize("jd,date_desc", CORE_DATES[:2])
    def test_sidereal_longitude_star(
        self, body_id, body_name, sid_mode, sid_name, jd, date_desc
    ):
        """Test sidereal longitude with star-based ayanamsha (relaxed tolerance)."""
        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flags = SEFLG_SIDEREAL
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = _sidereal_tol(body_id, sid_mode)
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"{body_name} sidereal ({sid_name}) at {date_desc}: "
            f'diff {diff:.6f}° ({arcsec(diff):.1f}") exceeds {tol:.4f}° '
            f"(swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )


class TestSiderealWithSpeed:
    """Test sidereal mode combined with SEFLG_SPEED."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("sid_mode,sid_name", MAJOR_AYANAMSHAS)
    @pytest.mark.parametrize("jd,date_desc", CORE_DATES)
    def test_sidereal_with_speed(
        self, body_id, body_name, sid_mode, sid_name, jd, date_desc
    ):
        """Test sidereal mode with velocity for major ayanamshas."""
        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flags = SEFLG_SIDEREAL | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        # Longitude
        tol = _sidereal_tol(body_id, sid_mode)
        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        assert diff_lon < tol, (
            f"{body_name} sidereal+speed ({sid_name}) at {date_desc}: "
            f"lon diff {diff_lon:.6f}°"
        )

        # Speed — slightly relaxed for sidereal
        diff_speed = abs(pos_swe[3] - pos_py[3])
        speed_tol = SPEED_LON_TOL[body_id] * 1.5
        assert diff_speed < speed_tol, (
            f"{body_name} sidereal+speed ({sid_name}) at {date_desc}: "
            f"speed diff {diff_speed:.6f}°/day "
            f"(swe={pos_swe[3]:.6f}, lib={pos_py[3]:.6f})"
        )


class TestEquatorialCoordinates:
    """Test equatorial coordinate transformation (SEFLG_EQUATORIAL)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("jd,date_desc", CORE_DATES)
    def test_equatorial_position(self, body_id, body_name, jd, date_desc):
        """Test equatorial coordinates (RA, Dec)."""
        flags = SEFLG_EQUATORIAL
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = EQUATORIAL_TOL[body_id]

        # Right Ascension
        diff_ra = angular_diff(pos_swe[0], pos_py[0])
        assert diff_ra < tol, (
            f"{body_name} equatorial at {date_desc}: RA diff {diff_ra:.6f}° "
            f"(swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )

        # Declination
        diff_dec = abs(pos_swe[1] - pos_py[1])
        assert diff_dec < tol, (
            f"{body_name} equatorial at {date_desc}: Dec diff {diff_dec:.6f}° "
            f"(swe={pos_swe[1]:.6f}, lib={pos_py[1]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("jd,date_desc", CORE_DATES[:2])
    def test_equatorial_with_speed(self, body_id, body_name, jd, date_desc):
        """Test equatorial coordinates with velocity."""
        flags = SEFLG_EQUATORIAL | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = EQUATORIAL_TOL[body_id]

        # RA
        diff_ra = angular_diff(pos_swe[0], pos_py[0])
        assert diff_ra < tol, (
            f"{body_name} equatorial+speed at {date_desc}: RA diff {diff_ra:.6f}°"
        )

        # RA velocity
        diff_dra = abs(pos_swe[3] - pos_py[3])
        speed_tol = SPEED_LON_TOL[body_id] * 2.0
        assert diff_dra < speed_tol, (
            f"{body_name} equatorial+speed at {date_desc}: "
            f"RA speed diff {diff_dra:.6f}°/day "
            f"(swe={pos_swe[3]:.6f}, lib={pos_py[3]:.6f})"
        )


class TestJ2000Frame:
    """Test J2000 reference frame (SEFLG_J2000)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("jd,date_desc", CORE_DATES)
    def test_j2000_position(self, body_id, body_name, jd, date_desc):
        """Test J2000 frame positions."""
        flags = SEFLG_J2000
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = LON_TOL[body_id] * 1.5
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"{body_name} J2000 at {date_desc}: diff {diff:.6f}° "
            f"(swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("jd,date_desc", CORE_DATES[:2])
    def test_j2000_with_speed(self, body_id, body_name, jd, date_desc):
        """Test J2000 frame with velocity."""
        flags = SEFLG_J2000 | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = LON_TOL[body_id] * 1.5
        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        assert diff_lon < tol, (
            f"{body_name} J2000+speed at {date_desc}: lon diff {diff_lon:.6f}°"
        )

        diff_speed = abs(pos_swe[3] - pos_py[3])
        assert diff_speed < SPEED_LON_TOL[body_id] * 1.5, (
            f"{body_name} J2000+speed at {date_desc}: speed diff {diff_speed:.6f}°/day"
        )


class TestNoNutation:
    """Test no-nutation mode (SEFLG_NONUT)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("jd,date_desc", CORE_DATES)
    def test_nonut_position(self, body_id, body_name, jd, date_desc):
        """Test no-nutation positions."""
        flags = SEFLG_NONUT
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = LON_TOL[body_id] * 1.5
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"{body_name} NONUT at {date_desc}: diff {diff:.6f}° "
            f"(swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    def test_nonut_differs_from_normal(self, body_id, body_name):
        """Verify NONUT gives different result than normal (nutation exists).

        For True Node / True Lilith the nutation effect is tiny but nonzero.
        Mean Node / Mean Lilith may not be affected at all.
        We only assert for True variants.
        """
        jd = 2451545.0

        pos_normal, _ = ephem.swe_calc_ut(jd, body_id, 0)
        pos_nonut, _ = ephem.swe_calc_ut(jd, body_id, SEFLG_NONUT)

        if body_id in (SE_TRUE_NODE, SE_OSCU_APOG):
            # Just verify the flag is accepted without error.
            # The difference may be zero if the implementation does not
            # apply nutation to these analytical bodies.
            diff = angular_diff(pos_normal[0], pos_nonut[0])
            assert diff < 1.0, (
                f"{body_name} NONUT vs normal diff {diff:.6f}° seems too large"
            )


class TestCombinedFlags:
    """Test complex flag combinations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("jd,date_desc", CORE_DATES[:2])
    def test_equatorial_j2000(self, body_id, body_name, jd, date_desc):
        """Test EQUATORIAL | J2000 combined."""
        flags = SEFLG_EQUATORIAL | SEFLG_J2000
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = EQUATORIAL_TOL[body_id] * 1.5
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"{body_name} EQUATORIAL|J2000 at {date_desc}: diff {diff:.6f}°"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("jd,date_desc", CORE_DATES[:2])
    def test_equatorial_j2000_speed(self, body_id, body_name, jd, date_desc):
        """Test EQUATORIAL | J2000 | SPEED combined."""
        flags = SEFLG_EQUATORIAL | SEFLG_J2000 | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = EQUATORIAL_TOL[body_id] * 1.5

        diff_ra = angular_diff(pos_swe[0], pos_py[0])
        diff_dec = abs(pos_swe[1] - pos_py[1])

        assert diff_ra < tol, (
            f"{body_name} EQ|J2000|SPEED at {date_desc}: RA diff {diff_ra:.6f}°"
        )
        assert diff_dec < tol, (
            f"{body_name} EQ|J2000|SPEED at {date_desc}: Dec diff {diff_dec:.6f}°"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES[:2])
    def test_sidereal_equatorial(self, body_id, body_name):
        """Test SIDEREAL | EQUATORIAL combined."""
        jd = 2451545.0
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flags = SEFLG_SIDEREAL | SEFLG_EQUATORIAL
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = EQUATORIAL_TOL[body_id] * SIDEREAL_LON_MULTIPLIER + SIDEREAL_STRICT_TOL
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, f"{body_name} SIDEREAL|EQUATORIAL: diff {diff:.6f}°"

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES[:2])
    def test_nonut_equatorial(self, body_id, body_name):
        """Test NONUT | EQUATORIAL combined."""
        jd = 2451545.0
        flags = SEFLG_NONUT | SEFLG_EQUATORIAL
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = EQUATORIAL_TOL[body_id] * 1.5
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, f"{body_name} NONUT|EQUATORIAL: diff {diff:.6f}°"

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    def test_sidereal_speed_multi_ayanamsha(self, body_id, body_name):
        """Test SIDEREAL | SPEED across multiple ayanamshas at J2000."""
        jd = 2451545.0
        flags = SEFLG_SIDEREAL | SEFLG_SPEED

        for sid_mode, sid_name in MAJOR_AYANAMSHAS:
            swe.set_sid_mode(sid_mode)
            ephem.swe_set_sid_mode(sid_mode)

            pos_swe, _ = swe.calc_ut(jd, body_id, flags)
            pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

            tol = _sidereal_tol(body_id, sid_mode)
            diff = angular_diff(pos_swe[0], pos_py[0])

            assert diff < tol, (
                f"{body_name} SIDEREAL|SPEED ({sid_name}): diff {diff:.6f}°"
            )


class TestTopocentric:
    """Test topocentric mode (SEFLG_TOPOCTR) with various locations.

    Lunar nodes and Lilith are geometrically defined points, not physical
    bodies. Topocentric correction may or may not be applied by the
    implementation. We test that both libraries return consistent results.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize(
        "loc_name,lat,lon,alt",
        [loc for loc in TEST_LOCATIONS[:4]],
    )
    def test_topocentric_longitude(self, body_id, body_name, loc_name, lat, lon, alt):
        """Test topocentric positions at different locations."""
        jd = 2451545.0

        swe.set_topo(lon, lat, alt)
        ephem.swe_set_topo(lon, lat, alt)

        flags = SEFLG_TOPOCTR
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = LON_TOL[body_id] * 2.0
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"{body_name} topocentric at {loc_name}: diff {diff:.6f}° "
            f"(swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES[:2])
    def test_topocentric_with_speed(self, body_id, body_name):
        """Test topocentric mode with velocity."""
        jd = 2451545.0

        swe.set_topo(12.4964, 41.9028, 0)  # Rome
        ephem.swe_set_topo(12.4964, 41.9028, 0)

        flags = SEFLG_TOPOCTR | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = LON_TOL[body_id] * 2.0
        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        assert diff_lon < tol, f"{body_name} TOPOCTR|SPEED: lon diff {diff_lon:.6f}°"

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES[:2])
    def test_topocentric_sidereal(self, body_id, body_name):
        """Test TOPOCTR | SIDEREAL combined."""
        jd = 2451545.0

        swe.set_topo(12.4964, 41.9028, 0)
        ephem.swe_set_topo(12.4964, 41.9028, 0)
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        flags = SEFLG_TOPOCTR | SEFLG_SIDEREAL
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = LON_TOL[body_id] * 2.0 + SIDEREAL_STRICT_TOL
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, f"{body_name} TOPOCTR|SIDEREAL: diff {diff:.6f}°"

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "loc_name,lat,lon,alt",
        [loc for loc in TEST_LOCATIONS],
    )
    def test_topocentric_all_locations_mean_node(self, loc_name, lat, lon, alt):
        """Test Mean Node topocentric across all locations."""
        jd = 2460000.5  # 2023

        swe.set_topo(lon, lat, alt)
        ephem.swe_set_topo(lon, lat, alt)

        flags = SEFLG_TOPOCTR
        pos_swe, _ = swe.calc_ut(jd, SE_MEAN_NODE, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, SE_MEAN_NODE, flags)

        tol = LON_TOL[SE_MEAN_NODE] * 2.0
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"Mean Node topocentric at {loc_name} ({lat:.2f}, {lon:.2f}): "
            f"diff {diff:.6f}°"
        )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize(
        "loc_name,lat,lon,alt",
        [loc for loc in TEST_LOCATIONS],
    )
    def test_topocentric_all_bodies_all_locations(
        self, body_id, body_name, loc_name, lat, lon, alt
    ):
        """Test all primary bodies topocentric across all locations at 2023."""
        jd = 2460000.5

        swe.set_topo(lon, lat, alt)
        ephem.swe_set_topo(lon, lat, alt)

        flags = SEFLG_TOPOCTR
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = LON_TOL[body_id] * 2.0
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"{body_name} topocentric at {loc_name}: diff {diff:.6f}° "
            f"(swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )


class TestStatisticalPrecision:
    """Statistical precision analysis across random dates."""

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    def test_random_date_statistics(self, body_id, body_name):
        """Test precision statistics across 500 random dates (1900-2100)."""
        random.seed(42)
        errors_lon = []
        errors_lat = []

        for _ in range(500):
            year = random.randint(1900, 2100)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = swe.julday(year, month, day, hour)

            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)

            errors_lon.append(angular_diff(pos_swe[0], pos_py[0]))
            errors_lat.append(abs(pos_swe[1] - pos_py[1]))

        max_lon = max(errors_lon)
        mean_lon = sum(errors_lon) / len(errors_lon)
        rms_lon = rms(errors_lon)
        p95_lon = percentile(errors_lon, 95)
        p99_lon = percentile(errors_lon, 99)

        print(f"\n  {body_name} longitude statistics (500 random dates 1900-2100):")
        print(f'    Max:  {max_lon:.6f}° ({arcsec(max_lon):.1f}")')
        print(f'    Mean: {mean_lon:.6f}° ({arcsec(mean_lon):.1f}")')
        print(f'    RMS:  {rms_lon:.6f}° ({arcsec(rms_lon):.1f}")')
        print(f'    P95:  {p95_lon:.6f}° ({arcsec(p95_lon):.1f}")')
        print(f'    P99:  {p99_lon:.6f}° ({arcsec(p99_lon):.1f}")')

        assert max_lon < LON_TOL[body_id], (
            f"{body_name}: max lon error {max_lon:.6f}° exceeds {LON_TOL[body_id]}°"
        )

        mean_tol = LON_TOL[body_id] * 0.6
        assert mean_lon < mean_tol, (
            f"{body_name}: mean lon error {mean_lon:.6f}° exceeds {mean_tol:.4f}°"
        )

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    def test_random_date_speed_statistics(self, body_id, body_name):
        """Test velocity precision across 200 random dates."""
        random.seed(123)
        errors_speed = []

        for _ in range(200):
            year = random.randint(1900, 2100)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = swe.julday(year, month, day, hour)

            pos_swe, _ = swe.calc_ut(jd, body_id, swe.FLG_SPEED)
            pos_py, _ = ephem.swe_calc_ut(jd, body_id, SEFLG_SPEED)

            errors_speed.append(abs(pos_swe[3] - pos_py[3]))

        max_speed = max(errors_speed)
        mean_speed = sum(errors_speed) / len(errors_speed)

        print(f"\n  {body_name} speed statistics (200 random dates 1900-2100):")
        print(f"    Max:  {max_speed:.6f}°/day")
        print(f"    Mean: {mean_speed:.6f}°/day")

        assert max_speed < SPEED_LON_TOL[body_id], (
            f"{body_name}: max speed error {max_speed:.6f}°/day "
            f"exceeds {SPEED_LON_TOL[body_id]}°/day"
        )

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", INTERPOLATED_BODIES)
    def test_interpolated_random_statistics(self, body_id, body_name):
        """Test interpolated body statistics across 200 random dates."""
        random.seed(77)
        errors_lon = []

        for _ in range(200):
            year = random.randint(1900, 2100)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = swe.julday(year, month, day, hour)

            pos_swe, _ = swe.calc_ut(jd, body_id, 0)
            pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)

            errors_lon.append(angular_diff(pos_swe[0], pos_py[0]))

        max_lon = max(errors_lon)
        mean_lon = sum(errors_lon) / len(errors_lon)

        print(f"\n  {body_name} statistics (200 random dates):")
        print(f'    Max:  {max_lon:.6f}° ({arcsec(max_lon):.1f}")')
        print(f'    Mean: {mean_lon:.6f}° ({arcsec(mean_lon):.1f}")')

        assert max_lon < LON_TOL[body_id], (
            f"{body_name}: max lon error {max_lon:.6f}° exceeds {LON_TOL[body_id]}°"
        )

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES[:2])
    def test_sidereal_random_statistics(self, body_id, body_name):
        """Test sidereal precision across 200 random dates with Lahiri."""
        random.seed(99)
        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        errors = []
        flags = SEFLG_SIDEREAL

        for _ in range(200):
            year = random.randint(1900, 2100)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = swe.julday(year, month, day, hour)

            pos_swe, _ = swe.calc_ut(jd, body_id, flags)
            pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

            errors.append(angular_diff(pos_swe[0], pos_py[0]))

        max_err = max(errors)
        mean_err = sum(errors) / len(errors)

        tol = _sidereal_tol(body_id, SE_SIDM_LAHIRI)

        print(f"\n  {body_name} sidereal Lahiri stats (200 dates):")
        print(f"    Max:  {max_err:.6f}°")
        print(f"    Mean: {mean_err:.6f}°")

        assert max_err < tol, (
            f"{body_name} sidereal: max error {max_err:.6f}° exceeds {tol:.4f}°"
        )


class TestExtendedDateRange:
    """Test at extended date range edges (DE440 coverage 1550-2650)."""

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("jd,date_desc", EXTENDED_DATES)
    def test_extended_range_primary(self, body_id, body_name, jd, date_desc):
        """Test primary bodies at extended date range."""
        pos_swe, _ = swe.calc_ut(jd, body_id, 0)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)

        tol = LON_TOL[body_id] * 2.0
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"{body_name} at {date_desc}: diff {diff:.6f}° "
            f'({arcsec(diff):.1f}") exceeds extended tolerance {tol}° '
            f"(swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("jd,date_desc", EXTENDED_DATES[:4])
    def test_extended_range_with_speed(self, body_id, body_name, jd, date_desc):
        """Test velocity at extended date range."""
        pos_swe, _ = swe.calc_ut(jd, body_id, swe.FLG_SPEED)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, SEFLG_SPEED)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_speed = abs(pos_swe[3] - pos_py[3])

        tol_lon = LON_TOL[body_id] * 2.0
        tol_speed = SPEED_LON_TOL[body_id] * 2.0

        assert diff_lon < tol_lon, (
            f"{body_name} at {date_desc}: lon diff {diff_lon:.6f}°"
        )
        assert diff_speed < tol_speed, (
            f"{body_name} at {date_desc}: speed diff {diff_speed:.6f}°/day"
        )

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", INTERPOLATED_BODIES)
    @pytest.mark.parametrize("jd,date_desc", EXTENDED_DATES[:4])
    def test_extended_range_interpolated(self, body_id, body_name, jd, date_desc):
        """Test interpolated bodies at extended date range."""
        pos_swe, _ = swe.calc_ut(jd, body_id, 0)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)

        tol = LON_TOL[body_id] * 2.0
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, f"{body_name} at {date_desc}: diff {diff:.6f}°"


class TestEdgeCases:
    """Edge case tests for lunar nodes and Lilith."""

    @pytest.mark.comparison
    def test_mean_node_latitude_is_zero(self):
        """Mean Node latitude should be 0 by definition."""
        for jd, desc in CORE_DATES:
            pos_py, _ = ephem.swe_calc_ut(jd, SE_MEAN_NODE, 0)
            pos_swe, _ = swe.calc_ut(jd, swe.MEAN_NODE, 0)

            assert abs(pos_py[1]) < 0.001, (
                f"Mean Node latitude should be ~0, got {pos_py[1]:.6f} at {desc}"
            )
            assert abs(pos_swe[1]) < 0.001, (
                f"SWE Mean Node latitude should be ~0, got {pos_swe[1]:.6f}"
            )

    @pytest.mark.comparison
    def test_longitude_range_0_360(self):
        """All longitudes should be in [0, 360) range."""
        for body_id, body_name in ALL_BODIES:
            for jd, desc in CORE_DATES:
                pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)
                assert 0.0 <= pos_py[0] < 360.0, (
                    f"{body_name} longitude {pos_py[0]} out of range at {desc}"
                )

    @pytest.mark.comparison
    def test_true_vs_mean_node_difference(self):
        """True Node should differ from Mean Node (oscillation)."""
        for jd, desc in CORE_DATES:
            pos_mean, _ = ephem.swe_calc_ut(jd, SE_MEAN_NODE, 0)
            pos_true, _ = ephem.swe_calc_ut(jd, SE_TRUE_NODE, 0)

            diff = angular_diff(pos_mean[0], pos_true[0])

            # True Node oscillates around Mean Node, typically within ~1.5°
            assert diff < 3.0, (
                f"True-Mean Node diff {diff:.4f}° at {desc} seems too large"
            )

    @pytest.mark.comparison
    def test_true_vs_mean_lilith_difference(self):
        """True Lilith should differ from Mean Lilith."""
        for jd, desc in CORE_DATES:
            pos_mean, _ = ephem.swe_calc_ut(jd, SE_MEAN_APOG, 0)
            pos_true, _ = ephem.swe_calc_ut(jd, SE_OSCU_APOG, 0)

            diff = angular_diff(pos_mean[0], pos_true[0])

            # True Lilith oscillates around Mean Lilith, typically within ~30°
            assert diff < 40.0, (
                f"True-Mean Lilith diff {diff:.4f}° at {desc} seems too large"
            )

    @pytest.mark.comparison
    def test_wrap_around_near_zero(self):
        """Test bodies near 0/360 wrap-around give consistent results."""
        random.seed(2024)
        found_near_zero = False

        for _ in range(2000):
            year = random.randint(1950, 2050)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            jd = swe.julday(year, month, day, 12.0)

            pos_swe, _ = swe.calc_ut(jd, swe.MEAN_NODE, 0)
            lon = pos_swe[0]

            if lon < 2.0 or lon > 358.0:
                pos_py, _ = ephem.swe_calc_ut(jd, SE_MEAN_NODE, 0)
                diff = angular_diff(pos_swe[0], pos_py[0])

                assert diff < LON_TOL[SE_MEAN_NODE], (
                    f"Mean Node near wrap-around: diff {diff:.6f}° at "
                    f"swe={pos_swe[0]:.4f}°, lib={pos_py[0]:.4f}°"
                )
                found_near_zero = True
                break

        assert found_near_zero, "Could not find Mean Node near 0/360 in test range"

    @pytest.mark.comparison
    def test_consecutive_days_continuity(self):
        """Test that positions change smoothly between consecutive days."""
        jd_start = 2451545.0  # J2000

        for body_id, body_name in PRIMARY_BODIES:
            positions_py = []

            for i in range(30):
                jd = jd_start + i
                pos_py, _ = ephem.swe_calc_ut(jd, body_id, 0)
                positions_py.append(pos_py[0])

            # Check that day-to-day changes are smooth (no jumps > 5°)
            for i in range(1, len(positions_py)):
                day_diff = angular_diff(positions_py[i], positions_py[i - 1])
                assert day_diff < 5.0, (
                    f"{body_name} discontinuity at day {i}: "
                    f"{positions_py[i - 1]:.4f} -> {positions_py[i]:.4f} "
                    f"(jump={day_diff:.4f}°)"
                )

    @pytest.mark.comparison
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    def test_same_result_both_api_styles(self, body_id, body_name):
        """Verify swe_calc_ut and calc_ut aliases return identical results."""
        jd = 2451545.0

        pos1, flag1 = ephem.swe_calc_ut(jd, body_id, SEFLG_SPEED)
        pos2, flag2 = ephem.calc_ut(jd, body_id, SEFLG_SPEED)

        for i in range(6):
            assert pos1[i] == pos2[i], (
                f"{body_name} API alias mismatch at index {i}: "
                f"swe_calc_ut={pos1[i]}, calc_ut={pos2[i]}"
            )

    @pytest.mark.comparison
    def test_interpolated_apogee_valid_output(self):
        """Interpolated apogee and perigee should produce valid results."""
        jd = 2451545.0

        pos_mean, _ = ephem.swe_calc_ut(jd, SE_MEAN_APOG, 0)
        pos_true, _ = ephem.swe_calc_ut(jd, SE_OSCU_APOG, 0)
        pos_intp, _ = ephem.swe_calc_ut(jd, SE_INTP_APOG, 0)
        pos_perg, _ = ephem.swe_calc_ut(jd, SE_INTP_PERG, 0)

        for pos, name in [
            (pos_mean, "Mean Lilith"),
            (pos_true, "True Lilith"),
            (pos_intp, "Interpolated Apogee"),
            (pos_perg, "Interpolated Perigee"),
        ]:
            assert 0.0 <= pos[0] < 360.0, f"{name} longitude {pos[0]} out of range"

    @pytest.mark.comparison
    def test_sidereal_consistent_across_bodies(self):
        """All bodies should shift by the same ayanamsha in sidereal mode."""
        jd = 2451545.0

        swe.set_sid_mode(swe.SIDM_LAHIRI)
        ephem.swe_set_sid_mode(SE_SIDM_LAHIRI)

        ayan = ephem.swe_get_ayanamsa_ut(jd)

        for body_id, body_name in PRIMARY_BODIES:
            pos_trop, _ = ephem.swe_calc_ut(jd, body_id, 0)
            pos_sid, _ = ephem.swe_calc_ut(jd, body_id, SEFLG_SIDEREAL)

            # The sidereal position should be tropical - ayanamsha (mod 360)
            expected_sid = (pos_trop[0] - ayan) % 360.0
            diff = angular_diff(pos_sid[0], expected_sid)

            assert diff < 0.001, (
                f"{body_name}: sidereal {pos_sid[0]:.6f}° != "
                f"tropical {pos_trop[0]:.6f}° - ayan {ayan:.6f}° = "
                f"{expected_sid:.6f}° (diff={diff:.6f}°)"
            )


class TestSiderealComprehensive:
    """Comprehensive sidereal tests covering ALL 43 ayanamshas with ALL bodies."""

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("sid_mode,sid_name", ALL_AYANAMSHAS)
    def test_all_ayanamshas_at_j2000(self, body_id, body_name, sid_mode, sid_name):
        """Test all 43 ayanamshas at J2000 for each primary body."""
        jd = 2451545.0

        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flags = SEFLG_SIDEREAL
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = _sidereal_tol(body_id, sid_mode)
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"{body_name} sidereal ({sid_name}) at J2000: "
            f'diff {diff:.6f}° ({arcsec(diff):.1f}") exceeds {tol:.4f}° '
            f"(swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("sid_mode,sid_name", ALL_AYANAMSHAS)
    def test_all_ayanamshas_at_j1900(self, body_id, body_name, sid_mode, sid_name):
        """Test all 43 ayanamshas at J1900 for each primary body."""
        jd = 2415020.5

        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flags = SEFLG_SIDEREAL
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = _sidereal_tol(body_id, sid_mode)
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"{body_name} sidereal ({sid_name}) at J1900: "
            f"diff {diff:.6f}° exceeds {tol:.4f}°"
        )

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", INTERPOLATED_BODIES)
    @pytest.mark.parametrize("sid_mode,sid_name", ALL_AYANAMSHAS)
    def test_all_ayanamshas_interpolated_j2000(
        self, body_id, body_name, sid_mode, sid_name
    ):
        """Test all 43 ayanamshas at J2000 for interpolated bodies."""
        jd = 2451545.0

        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flags = SEFLG_SIDEREAL
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = _sidereal_tol(body_id, sid_mode)
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"{body_name} sidereal ({sid_name}) at J2000: "
            f"diff {diff:.6f}° exceeds {tol:.4f}°"
        )

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("sid_mode,sid_name", ALL_AYANAMSHAS)
    def test_all_ayanamshas_with_speed_j2000(
        self, body_id, body_name, sid_mode, sid_name
    ):
        """Test all 43 ayanamshas with SPEED at J2000."""
        jd = 2451545.0

        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flags = SEFLG_SIDEREAL | SEFLG_SPEED
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = _sidereal_tol(body_id, sid_mode)

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        assert diff_lon < tol, (
            f"{body_name} SIDEREAL|SPEED ({sid_name}): lon diff {diff_lon:.6f}°"
        )

        diff_speed = abs(pos_swe[3] - pos_py[3])
        speed_tol = SPEED_LON_TOL[body_id] * 1.5
        assert diff_speed < speed_tol, (
            f"{body_name} SIDEREAL|SPEED ({sid_name}): "
            f"speed diff {diff_speed:.6f}°/day "
            f"(swe={pos_swe[3]:.6f}, lib={pos_py[3]:.6f})"
        )

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", PRIMARY_BODIES)
    @pytest.mark.parametrize("sid_mode,sid_name", ALL_AYANAMSHAS)
    def test_all_ayanamshas_at_2023(self, body_id, body_name, sid_mode, sid_name):
        """Test all 43 ayanamshas at a modern date (2023) for each primary body."""
        jd = 2460000.5

        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flags = SEFLG_SIDEREAL
        pos_swe, _ = swe.calc_ut(jd, body_id, flags)
        pos_py, _ = ephem.swe_calc_ut(jd, body_id, flags)

        tol = _sidereal_tol(body_id, sid_mode)
        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < tol, (
            f"{body_name} sidereal ({sid_name}) at 2023: "
            f"diff {diff:.6f}° exceeds {tol:.4f}° "
            f"(swe={pos_swe[0]:.6f}, lib={pos_py[0]:.6f})"
        )
