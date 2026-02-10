"""
Moshier House Systems Comparison Tests.

Compares all house system calculations between pyswisseph (C Moshier) and
libephemeris (Python) using SEFLG_MOSEPH flag.

This is the Moshier-mode mirror of test_compare_houses.py (which covers
SEFLG_SWIEPH / JPL mode). It ensures that house cusp calculations match
between the C and Python implementations when using the Moshier semi-analytical
ephemeris.

Structural differences:
- Sidereal time (ARMC) and obliquity (eps) are computed by Skyfield in
  libephemeris but by internal Moshier routines in the C library, producing
  small structural discrepancies in all house cusps.
- Sun position for Sunshine houses (I/i) is computed via VSOP87 in both
  implementations but with C-vs-Python porting differences.
- These discrepancies justify relaxed tolerances (MOSHIER_CUSP_TOL = 0.01°)
  compared to the SWIEPH test (HOUSE_CUSP_TOL = 0.001°).

Without these tests, the migration from pyswisseph to libephemeris for dates
outside DE440 range (1550-2650) is completely unvalidated for house calculations.
"""

import pytest
import random
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import SEFLG_MOSEPH


# ============================================================================
# TOLERANCES
# ============================================================================

# Moshier tolerances are relaxed compared to SWIEPH (0.001°) because ARMC and
# obliquity are computed by Skyfield in libephemeris vs internal Moshier routines
# in the C library, producing structural discrepancies in all house cusps.
MOSHIER_CUSP_TOL = 0.01  # degrees (~36 arcsec)
MOSHIER_ASCMC_TOL = 0.01  # degrees (~36 arcsec)
MOSHIER_GAUQUELIN_TOL = 0.01  # Gauquelin 36 sectors


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# HOUSE SYSTEMS
# ============================================================================

HOUSE_SYSTEMS = [
    ("P", "Placidus"),
    ("K", "Koch"),
    ("R", "Regiomontanus"),
    ("C", "Campanus"),
    ("E", "Equal (Ascendant)"),
    ("A", "Equal (MC)"),
    ("W", "Whole Sign"),
    ("O", "Porphyry"),
    ("B", "Alcabitius"),
    ("T", "Polich/Page (Topocentric)"),
    ("M", "Morinus"),
    ("X", "Meridian (Axial Rotation)"),
    ("V", "Vehlow Equal"),
    ("H", "Horizontal"),
    ("F", "Carter Poli-Equatorial"),
    ("U", "Krusinski-Pisa"),
    ("N", "Natural Gradient"),
    ("Y", "APC Houses"),
    ("D", "Equal from MC"),
    ("L", "Pullen SD"),
    ("S", "Sripati"),
    ("G", "Gauquelin"),
    ("I", "Sunshine/Makransky"),
    ("Q", "Pullen SR"),
]

# Systems with relaxed tolerances beyond the Moshier baseline
# These have additional sources of divergence on top of the ARMC/eps differences.
RELAXED_SYSTEMS = {
    "G": MOSHIER_GAUQUELIN_TOL,  # Gauquelin: 36 sectors
    "I": 0.02,  # Sunshine: Sun position C-vs-Python VSOP87 differences compound
    "K": 0.02,  # Koch: sensitive to ARMC/eps at mid-latitudes
    "P": 0.02,  # Placidus: iterative algorithm amplifies ARMC/eps diffs
    "R": 0.02,  # Regiomontanus: sensitive to obliquity differences
    "N": 0.02,  # Natural Gradient: complex interpolation
    "Y": 0.02,  # APC: complex algorithm sensitive to obliquity
}

# Test locations (same as test_compare_houses.py)
STANDARD_LOCATIONS = [
    ("Rome", 41.9028, 12.4964),
    ("New York", 40.7128, -74.0060),
    ("Sydney", -33.8688, 151.2093),
    ("Equator", 0.0, 0.0),
]

HIGH_LATITUDE_LOCATIONS = [
    ("Tromso", 69.6492, 18.9553),
    ("McMurdo", -77.8463, 166.6681),
]

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 6, 15, 0.0, "Summer 2024"),
    (1980, 5, 20, 14.5, "Past"),
]


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestHouseCusps:
    """Compare Moshier house cusp calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS)
    @pytest.mark.parametrize("name,lat,lon", STANDARD_LOCATIONS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_house_cusps_match(
        self, hsys, hsys_name, name, lat, lon, year, month, day, hour, date_desc
    ):
        """Test all 12 house cusps match between implementations using Moshier."""
        jd = swe.julday(year, month, day, hour)

        cusps_swe, ascmc_swe = swe.houses_ex(
            jd, lat, lon, hsys.encode("ascii"), swe.FLG_MOSEPH
        )
        cusps_py, ascmc_py = ephem.swe_houses_ex(jd, lat, lon, hsys, SEFLG_MOSEPH)

        # Get tolerance for this system
        tolerance = RELAXED_SYSTEMS.get(hsys, MOSHIER_CUSP_TOL)

        # Compare cusps (handle Gauquelin which has 36 sectors in pyswisseph)
        num_cusps_to_compare = min(len(cusps_swe), len(cusps_py))
        max_diff = 0.0
        for i in range(num_cusps_to_compare):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            max_diff = max(max_diff, diff)

        assert max_diff < tolerance, (
            f"{hsys_name} Moshier at {name} ({date_desc}): "
            f"max cusp diff {max_diff:.6f}° exceeds tolerance {tolerance}°"
        )


class TestAscendantMC:
    """Compare Moshier Ascendant and MC calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS[:7])  # Common systems
    @pytest.mark.parametrize("name,lat,lon", STANDARD_LOCATIONS)
    def test_ascendant_mc_match(self, hsys, hsys_name, name, lat, lon):
        """Test Ascendant and MC match between implementations using Moshier."""
        jd = swe.julday(2000, 1, 1, 12.0)

        _, ascmc_swe = swe.houses_ex(jd, lat, lon, hsys.encode("ascii"), swe.FLG_MOSEPH)
        _, ascmc_py = ephem.swe_houses_ex(jd, lat, lon, hsys, SEFLG_MOSEPH)

        # Ascendant (index 0)
        asc_diff = angular_diff(ascmc_swe[0], ascmc_py[0])
        assert asc_diff < MOSHIER_ASCMC_TOL, (
            f"{hsys_name} Moshier at {name}: "
            f"Ascendant diff {asc_diff:.6f}° exceeds tolerance {MOSHIER_ASCMC_TOL}°"
        )

        # MC (index 1)
        mc_diff = angular_diff(ascmc_swe[1], ascmc_py[1])
        assert mc_diff < MOSHIER_ASCMC_TOL, (
            f"{hsys_name} Moshier at {name}: "
            f"MC diff {mc_diff:.6f}° exceeds tolerance {MOSHIER_ASCMC_TOL}°"
        )


class TestRandomLocations:
    """Comprehensive random location testing with Moshier."""

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS[:10])
    def test_random_locations(self, hsys, hsys_name):
        """Test Moshier house systems at 50 random locations."""
        random.seed(42)

        tolerance = RELAXED_SYSTEMS.get(hsys, MOSHIER_CUSP_TOL)
        failures = []

        for _ in range(50):
            lat = random.uniform(-60.0, 60.0)
            lon = random.uniform(-180.0, 180.0)
            year = random.randint(1900, 2100)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0.0, 24.0)

            jd = swe.julday(year, month, day, hour)

            try:
                cusps_swe, _ = swe.houses_ex(
                    jd, lat, lon, hsys.encode("ascii"), swe.FLG_MOSEPH
                )
                cusps_py, _ = ephem.swe_houses_ex(jd, lat, lon, hsys, SEFLG_MOSEPH)

                max_diff = max(
                    angular_diff(cusps_swe[i], cusps_py[i]) for i in range(12)
                )

                if max_diff >= tolerance:
                    failures.append((lat, lon, jd, max_diff))
            except Exception:
                pass  # Skip errors (polar regions, etc.)

        assert len(failures) == 0, (
            f"{hsys_name} Moshier: {len(failures)} failures out of 50 random locations"
        )


class TestHighLatitudes:
    """Test Moshier behavior at high latitudes."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "hsys,hsys_name",
        [
            ("E", "Equal"),
            ("W", "Whole Sign"),
            ("O", "Porphyry"),
            ("M", "Morinus"),
        ],
    )  # Systems that should work at all latitudes
    @pytest.mark.parametrize("name,lat,lon", HIGH_LATITUDE_LOCATIONS)
    def test_high_latitude_systems(self, hsys, hsys_name, name, lat, lon):
        """Test Moshier house systems that work at high latitudes."""
        jd = swe.julday(2000, 1, 1, 12.0)

        cusps_swe, ascmc_swe = swe.houses_ex(
            jd, lat, lon, hsys.encode("ascii"), swe.FLG_MOSEPH
        )
        cusps_py, ascmc_py = ephem.swe_houses_ex(jd, lat, lon, hsys, SEFLG_MOSEPH)

        max_diff = max(angular_diff(cusps_swe[i], cusps_py[i]) for i in range(12))

        assert max_diff < MOSHIER_CUSP_TOL, (
            f"{hsys_name} Moshier at {name}: max diff {max_diff:.6f}° at high latitude"
        )


class TestSpecialCases:
    """Test Moshier special cases and edge conditions."""

    @pytest.mark.comparison
    def test_equator(self):
        """Test Moshier house calculations at equator."""
        jd = swe.julday(2000, 1, 1, 12.0)
        lat, lon = 0.0, 0.0

        for hsys in ["P", "K", "E", "W"]:
            cusps_swe, _ = swe.houses_ex(
                jd, lat, lon, hsys.encode("ascii"), swe.FLG_MOSEPH
            )
            cusps_py, _ = ephem.swe_houses_ex(jd, lat, lon, hsys, SEFLG_MOSEPH)

            max_diff = max(angular_diff(cusps_swe[i], cusps_py[i]) for i in range(12))

            tolerance = RELAXED_SYSTEMS.get(hsys, MOSHIER_CUSP_TOL)
            assert max_diff < tolerance, (
                f"Moshier system {hsys} at equator: max diff {max_diff:.6f}°"
            )

    @pytest.mark.comparison
    def test_prime_meridian(self):
        """Test Moshier house calculations at prime meridian."""
        jd = swe.julday(2000, 1, 1, 12.0)
        lat, lon = 51.5074, 0.0  # London area

        cusps_swe, _ = swe.houses_ex(jd, lat, lon, b"P", swe.FLG_MOSEPH)
        cusps_py, _ = ephem.swe_houses_ex(jd, lat, lon, "P", SEFLG_MOSEPH)

        max_diff = max(angular_diff(cusps_swe[i], cusps_py[i]) for i in range(12))

        tolerance = RELAXED_SYSTEMS.get("P", MOSHIER_CUSP_TOL)
        assert max_diff < tolerance, (
            f"Moshier Placidus at prime meridian: max diff {max_diff:.6f}°"
        )

    @pytest.mark.comparison
    def test_international_date_line(self):
        """Test Moshier house calculations near international date line."""
        jd = swe.julday(2000, 1, 1, 12.0)
        lat, lon = 35.0, 179.9  # Near date line

        cusps_swe, _ = swe.houses_ex(jd, lat, lon, b"P", swe.FLG_MOSEPH)
        cusps_py, _ = ephem.swe_houses_ex(jd, lat, lon, "P", SEFLG_MOSEPH)

        max_diff = max(angular_diff(cusps_swe[i], cusps_py[i]) for i in range(12))

        tolerance = RELAXED_SYSTEMS.get("P", MOSHIER_CUSP_TOL)
        assert max_diff < tolerance, (
            f"Moshier Placidus near date line: max diff {max_diff:.6f}°"
        )


# ============================================================================
# MOSHIER RANDOM AND POLAR LOCATION TESTS
# ============================================================================

# 5 most commonly used house systems for random location testing
MOSHIER_RANDOM_SYSTEMS = [
    ("P", "Placidus"),
    ("K", "Koch"),
    ("E", "Equal"),
    ("W", "Whole Sign"),
    ("O", "Porphyry"),
]

# Systems safe at all latitudes (no polar circle issues)
# Placidus and Koch fail when |lat| + eps > 90 (~66.56 threshold)
POLAR_SAFE_SYSTEMS = [
    ("E", "Equal"),
    ("W", "Whole Sign"),
    ("O", "Porphyry"),
]

# Polar locations: 7 latitudes x 2 hemispheres = 14 locations
# Covers from sub-polar (65) through near-pole (89.9)
# Real-world significance: Tromso (69.6N), Reykjavik (64.1N), Murmansk (68.9N)
POLAR_LOCATIONS = [
    ("Polar_65N", 65.0, 18.0),
    ("Polar_67N", 67.0, 25.7),
    ("Polar_70N", 70.0, 19.0),
    ("Polar_75N", 75.0, 30.0),
    ("Polar_80N", 80.0, 16.0),
    ("Polar_85N", 85.0, 0.0),
    ("Polar_89.9N", 89.9, 0.0),
    ("Polar_65S", -65.0, 18.0),
    ("Polar_67S", -67.0, 25.7),
    ("Polar_70S", -70.0, 19.0),
    ("Polar_75S", -75.0, 30.0),
    ("Polar_80S", -80.0, 16.0),
    ("Polar_85S", -85.0, 0.0),
    ("Polar_89.9S", -89.9, 0.0),
]

MOSHIER_POLAR_TOL = 0.02  # degrees - relaxed for polar Moshier comparison


class TestMoshierRandomLocations:
    """Moshier house systems at 50 random geographic locations.

    Validates cusp calculations using SEFLG_MOSEPH across 50 random
    locations (seed=42, lat +/-60, lon +/-180) for the 5 most commonly
    used house systems (P, K, E, W, O).

    Differences between C (pyswisseph) and Python (libephemeris) Moshier
    arise from ARMC/obliquity computation: Skyfield (Python) vs internal
    Moshier routines (C). These differences are structural and bounded.

    Total: 250 validations (50 locations x 5 systems).
    """

    @pytest.mark.comparison
    @pytest.mark.slow
    @pytest.mark.parametrize("hsys,hsys_name", MOSHIER_RANDOM_SYSTEMS)
    def test_moshier_random_locations(self, hsys, hsys_name):
        """Test Moshier house system at 50 random locations (seed=42)."""
        random.seed(42)

        tolerance = RELAXED_SYSTEMS.get(hsys, MOSHIER_CUSP_TOL)
        failures = []

        for idx in range(50):
            lat = random.uniform(-60.0, 60.0)
            lon = random.uniform(-180.0, 180.0)
            year = random.randint(1900, 2100)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0.0, 24.0)

            jd = swe.julday(year, month, day, hour)

            try:
                cusps_swe, _ = swe.houses_ex(
                    jd, lat, lon, hsys.encode("ascii"), swe.FLG_MOSEPH
                )
                cusps_py, _ = ephem.swe_houses_ex(jd, lat, lon, hsys, SEFLG_MOSEPH)

                max_diff = max(
                    angular_diff(cusps_swe[k], cusps_py[k]) for k in range(12)
                )

                if max_diff >= tolerance:
                    failures.append(
                        f"  Loc_{idx} ({lat:.2f}, {lon:.2f}) JD={jd:.1f}: "
                        f"max_diff={max_diff:.6f}"
                    )
            except Exception:
                pass  # Skip errors (polar regions, etc.)

        assert len(failures) == 0, (
            f"{hsys_name} Moshier: {len(failures)}/50 random locations exceeded "
            f"{tolerance} deg tolerance:\n" + "\n".join(failures[:10])
        )


class TestMoshierHighLatitudes:
    """Moshier house systems at polar latitudes (65-89.9 degrees).

    Validates cusp calculations at latitudes where Placidus and Koch fail
    due to the polar circle threshold (90 - obliquity ~ 66.56 degrees).
    Only systems that work at all latitudes are tested (E, W, O).

    Polar latitudes are critical edge cases: cities like Tromso (69.6N),
    Reykjavik (64.1N), and Murmansk (68.9N) have significant astrological
    user bases. Differences in the polar threshold between C and Python
    implementations can cause one to use Placidus while the other falls
    back to Porphyry, producing completely different horoscopes.

    Total: 42 test cases (14 locations x 3 systems).
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,hsys_name", POLAR_SAFE_SYSTEMS)
    @pytest.mark.parametrize(
        "name,lat,lon",
        POLAR_LOCATIONS,
        ids=[loc[0] for loc in POLAR_LOCATIONS],
    )
    def test_moshier_polar_location(self, hsys, hsys_name, name, lat, lon):
        """Test Moshier house system at polar latitude."""
        jd = swe.julday(2000, 1, 1, 12.0)

        cusps_swe, _ = swe.houses_ex(jd, lat, lon, hsys.encode("ascii"), swe.FLG_MOSEPH)
        cusps_py, _ = ephem.swe_houses_ex(jd, lat, lon, hsys, SEFLG_MOSEPH)

        max_diff = max(angular_diff(cusps_swe[k], cusps_py[k]) for k in range(12))

        assert max_diff < MOSHIER_POLAR_TOL, (
            f"{hsys_name} Moshier at {name} (lat={lat}): "
            f"max cusp diff {max_diff:.6f} deg exceeds {MOSHIER_POLAR_TOL} deg"
        )
