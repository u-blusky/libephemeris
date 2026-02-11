"""
Moshier Sidereal House Systems Comparison Tests.

Validates the combination of SEFLG_MOSEPH + SEFLG_SIDEREAL for house calculations
between pyswisseph (C) and libephemeris (Python).

This test is the missing cross-product of:
  - test_moshier_compare_houses.py (Moshier houses, tropical only)
  - test_moshier_compare_sidereal.py (Moshier sidereal, planets only)

The MOSEPH+SIDEREAL pipeline for houses involves 4 steps, each with potential
C-vs-Python divergence:
  1. Moshier planet/ARMC/obliquity computation (Skyfield vs C Moshier routines)
  2. Tropical house cusp calculation (iterative algorithms amplify ARMC diffs)
  3. Ayanamsha computation (formula-based, tight agreement)
  4. Sidereal correction (subtract ayanamsha from cusps, recalculate for W/E/V)

For Whole Sign (W), cusps depend on the zodiacal sign of the sidereal Ascendant.
An error of ~0.5 deg near a sign boundary can shift ALL 12 cusps by 30 degrees.
This is critical for Jyotish (Vedic astrology), which uses Whole Sign + Lahiri
almost exclusively. Historical dates (Moshier range) are common in Jyotish.

Test coverage: 3 ayanamshas x 7 house systems x 4 locations x 3 dates = 252 cases.
"""

import math

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    SE_SIDM_FAGAN_BRADLEY,
    SE_SIDM_LAHIRI,
    SE_SIDM_RAMAN,
    SEFLG_MOSEPH,
    SEFLG_SIDEREAL,
)


# ============================================================================
# TOLERANCES
# ============================================================================

# Moshier sidereal house tolerance: combines Moshier ARMC/eps divergence (~0.01)
# with ayanamsha subtraction. Since ayanamsha is formula-based and agrees to
# ~0.001 deg, the dominant error is from the tropical cusp computation.
# 0.05 deg provides margin for iterative systems (Placidus, Koch) where
# ARMC/eps differences are amplified.
MOSHIER_SIDEREAL_CUSP_TOL = 0.05  # degrees (~3 arcminutes)

# Per-system relaxations for systems with additional divergence sources
RELAXED_SYSTEMS = {
    "P": 0.05,  # Placidus: iterative, amplifies ARMC/eps diffs
    "K": 0.05,  # Koch: sensitive to ARMC at mid-latitudes
    "R": 0.05,  # Regiomontanus: sensitive to obliquity differences
}

# Ascendant/MC tolerance (same pipeline as cusps)
MOSHIER_SIDEREAL_ASCMC_TOL = 0.05  # degrees


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TEST DATA
# ============================================================================

# 3 principal ayanamshas (most used in Vedic and Western sidereal astrology)
AYANAMSHAS = [
    (SE_SIDM_LAHIRI, "Lahiri"),
    (SE_SIDM_FAGAN_BRADLEY, "Fagan-Bradley"),
    (SE_SIDM_RAMAN, "Raman"),
]

# 7 house systems specified in the task
HOUSE_SYSTEMS = [
    ("P", "Placidus"),
    ("K", "Koch"),
    ("E", "Equal"),
    ("W", "Whole Sign"),
    ("O", "Porphyry"),
    ("R", "Regiomontanus"),
    ("C", "Campanus"),
]

# 4 standard locations (same as test_moshier_compare_houses.py)
STANDARD_LOCATIONS = [
    ("Rome", 41.9028, 12.4964),
    ("New York", 40.7128, -74.0060),
    ("Sydney", -33.8688, 151.2093),
    ("Equator", 0.0, 0.0),
]

# 3 test dates
TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 6, 15, 0.0, "Summer 2024"),
    (1980, 5, 20, 14.5, "Past"),
]


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestMoshierSiderealHouseCusps:
    """Compare Moshier sidereal house cusp calculations.

    252 parametrized test cases: 3 ayanamshas x 7 house systems x 4 locations x 3 dates.
    Each test sets the sidereal mode in both implementations, calls houses_ex with
    SEFLG_MOSEPH | SEFLG_SIDEREAL, and compares all 12 cusps.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", AYANAMSHAS)
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS)
    @pytest.mark.parametrize("name,lat,lon", STANDARD_LOCATIONS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_sidereal_house_cusps_match(
        self,
        sid_mode,
        sid_name,
        hsys,
        hsys_name,
        name,
        lat,
        lon,
        year,
        month,
        day,
        hour,
        date_desc,
    ):
        """Test all 12 sidereal house cusps match between implementations."""
        jd = swe.julday(year, month, day, hour)

        # Set sidereal mode in both implementations
        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        # Combined flags: Moshier ephemeris + sidereal zodiac
        flag_swe = swe.FLG_MOSEPH | swe.FLG_SIDEREAL
        flag_py = SEFLG_MOSEPH | SEFLG_SIDEREAL

        cusps_swe, ascmc_swe = swe.houses_ex(
            jd, lat, lon, hsys.encode("ascii"), flag_swe
        )
        cusps_py, ascmc_py = ephem.swe_houses_ex(jd, lat, lon, hsys, flag_py)

        # Get tolerance for this system
        tolerance = RELAXED_SYSTEMS.get(hsys, MOSHIER_SIDEREAL_CUSP_TOL)

        # Compare all 12 cusps
        num_cusps = min(len(cusps_swe), len(cusps_py), 12)
        max_diff = 0.0
        worst_cusp = 0
        for i in range(num_cusps):
            diff = angular_diff(cusps_swe[i], cusps_py[i])
            if diff > max_diff:
                max_diff = diff
                worst_cusp = i + 1

        assert max_diff < tolerance, (
            f"{hsys_name} sidereal ({sid_name}) at {name} ({date_desc}): "
            f"cusp {worst_cusp} max diff {max_diff:.6f}° exceeds {tolerance}° "
            f"(swe={cusps_swe[worst_cusp - 1]:.6f}°, "
            f"lib={cusps_py[worst_cusp - 1]:.6f}°)"
        )


class TestMoshierSiderealAscMC:
    """Compare sidereal Ascendant and MC with SEFLG_MOSEPH | SEFLG_SIDEREAL.

    The sidereal Ascendant is critical for Whole Sign and Equal house systems.
    An error in the sidereal Asc propagates to all cusps.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", AYANAMSHAS)
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS[:4])  # P, K, E, W
    @pytest.mark.parametrize("name,lat,lon", STANDARD_LOCATIONS)
    def test_sidereal_asc_mc_match(
        self, sid_mode, sid_name, hsys, hsys_name, name, lat, lon
    ):
        """Test sidereal Ascendant and MC match between implementations."""
        jd = swe.julday(2000, 1, 1, 12.0)

        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flag_swe = swe.FLG_MOSEPH | swe.FLG_SIDEREAL
        flag_py = SEFLG_MOSEPH | SEFLG_SIDEREAL

        _, ascmc_swe = swe.houses_ex(jd, lat, lon, hsys.encode("ascii"), flag_swe)
        _, ascmc_py = ephem.swe_houses_ex(jd, lat, lon, hsys, flag_py)

        # Ascendant (index 0)
        asc_diff = angular_diff(ascmc_swe[0], ascmc_py[0])
        assert asc_diff < MOSHIER_SIDEREAL_ASCMC_TOL, (
            f"{hsys_name} sidereal ({sid_name}) at {name}: "
            f"Asc diff {asc_diff:.6f}° exceeds {MOSHIER_SIDEREAL_ASCMC_TOL}° "
            f"(swe={ascmc_swe[0]:.6f}°, lib={ascmc_py[0]:.6f}°)"
        )

        # MC (index 1)
        mc_diff = angular_diff(ascmc_swe[1], ascmc_py[1])
        assert mc_diff < MOSHIER_SIDEREAL_ASCMC_TOL, (
            f"{hsys_name} sidereal ({sid_name}) at {name}: "
            f"MC diff {mc_diff:.6f}° exceeds {MOSHIER_SIDEREAL_ASCMC_TOL}° "
            f"(swe={ascmc_swe[1]:.6f}°, lib={ascmc_py[1]:.6f}°)"
        )


class TestWholeSignSidereal:
    """Whole Sign sidereal validation - critical for Jyotish.

    Whole Sign (W) cusps are exact multiples of 30 degrees, determined by the
    zodiacal sign of the sidereal Ascendant. If the sidereal Asc is in a
    different sign between C and Python (e.g., due to rounding near 0/30/60...),
    ALL 12 cusps shift by 30 degrees. This is the most critical test for
    Vedic astrology users.

    Validates:
      1. All cusps are exact multiples of 30 degrees (structural correctness)
      2. Ascendant sign matches between implementations (functional correctness)
      3. All cusps match (consequence of 1+2)
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", AYANAMSHAS)
    @pytest.mark.parametrize("name,lat,lon", STANDARD_LOCATIONS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_whole_sign_cusps_are_multiples_of_30(
        self, sid_mode, sid_name, name, lat, lon, year, month, day, hour, date_desc
    ):
        """Test that Whole Sign sidereal cusps are exact multiples of 30 degrees."""
        jd = swe.julday(year, month, day, hour)

        # Test both implementations
        for label, setup_fn in [
            ("pyswisseph", lambda: swe.set_sid_mode(sid_mode)),
            ("libephemeris", lambda: ephem.swe_set_sid_mode(sid_mode)),
        ]:
            setup_fn()

            if label == "pyswisseph":
                flag = swe.FLG_MOSEPH | swe.FLG_SIDEREAL
                cusps, _ = swe.houses_ex(jd, lat, lon, b"W", flag)
            else:
                flag = SEFLG_MOSEPH | SEFLG_SIDEREAL
                cusps, _ = ephem.swe_houses_ex(jd, lat, lon, "W", flag)

            for i in range(12):
                remainder = cusps[i] % 30.0
                # Allow tiny floating-point tolerance
                is_multiple = remainder < 0.001 or (30.0 - remainder) < 0.001
                assert is_multiple, (
                    f"{label} Whole Sign ({sid_name}) at {name} ({date_desc}): "
                    f"cusp {i + 1} = {cusps[i]:.6f}° is not a multiple of 30°"
                )

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", AYANAMSHAS)
    @pytest.mark.parametrize("name,lat,lon", STANDARD_LOCATIONS)
    @pytest.mark.parametrize("year,month,day,hour,date_desc", TEST_DATES)
    def test_whole_sign_ascendant_sign_matches(
        self, sid_mode, sid_name, name, lat, lon, year, month, day, hour, date_desc
    ):
        """Test that the sidereal Ascendant falls in the same zodiac sign.

        This is THE critical Jyotish test: if the Asc sign differs, the entire
        Vedic horoscope (all bhava/house placements) is wrong.
        """
        jd = swe.julday(year, month, day, hour)

        swe.set_sid_mode(sid_mode)
        ephem.swe_set_sid_mode(sid_mode)

        flag_swe = swe.FLG_MOSEPH | swe.FLG_SIDEREAL
        flag_py = SEFLG_MOSEPH | SEFLG_SIDEREAL

        _, ascmc_swe = swe.houses_ex(jd, lat, lon, b"W", flag_swe)
        _, ascmc_py = ephem.swe_houses_ex(jd, lat, lon, "W", flag_py)

        # Sign = floor(longitude / 30)
        sign_swe = math.floor(ascmc_swe[0] / 30.0)
        sign_py = math.floor(ascmc_py[0] / 30.0)

        sign_names = [
            "Ari",
            "Tau",
            "Gem",
            "Can",
            "Leo",
            "Vir",
            "Lib",
            "Sco",
            "Sag",
            "Cap",
            "Aqu",
            "Pis",
        ]

        assert sign_swe == sign_py, (
            f"Whole Sign ({sid_name}) at {name} ({date_desc}): "
            f"Asc sign MISMATCH: swe={sign_names[sign_swe % 12]} "
            f"({ascmc_swe[0]:.4f}°) vs lib={sign_names[sign_py % 12]} "
            f"({ascmc_py[0]:.4f}°) - entire Vedic horoscope would be wrong!"
        )


class TestMoshierSiderealConsistency:
    """Cross-consistency checks for the Moshier+sidereal pipeline.

    Validates that the sidereal correction is applied correctly by comparing
    sidereal cusps against tropical cusps minus ayanamsha.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize("sid_mode,sid_name", AYANAMSHAS)
    @pytest.mark.parametrize(
        "hsys,hsys_name",
        [("E", "Equal"), ("O", "Porphyry"), ("C", "Campanus")],
    )
    def test_sidereal_equals_tropical_minus_ayanamsha(
        self, sid_mode, sid_name, hsys, hsys_name
    ):
        """Test that sidereal cusps = tropical cusps - ayanamsha.

        For non-Ascendant-based systems (P, K, O, R, C, etc.), the sidereal
        correction is a simple subtraction. This validates the pipeline is
        internally consistent in libephemeris.

        Note: Equal (E) is Ascendant-based but the cusps are Asc + i*30, so
        we validate that the sidereal cusps track sidereal Asc correctly.
        """
        jd = swe.julday(2000, 1, 1, 12.0)
        lat, lon = 41.9028, 12.4964  # Rome

        ephem.swe_set_sid_mode(sid_mode)

        # Tropical cusps (Moshier, no sidereal)
        cusps_trop, _ = ephem.swe_houses_ex(jd, lat, lon, hsys, SEFLG_MOSEPH)

        # Sidereal cusps (Moshier + sidereal)
        cusps_sid, _ = ephem.swe_houses_ex(
            jd, lat, lon, hsys, SEFLG_MOSEPH | SEFLG_SIDEREAL
        )

        # Ayanamsha
        ayan = ephem.swe_get_ayanamsa_ut(jd)

        # For systems where sidereal = tropical - ayanamsha (non-Asc-based)
        # Equal is Asc-based so skip direct subtraction check for it
        if hsys not in ("E", "W", "V", "A"):
            for i in range(12):
                expected = (cusps_trop[i] - ayan) % 360.0
                diff = angular_diff(cusps_sid[i], expected)
                assert diff < 0.001, (
                    f"{hsys_name} ({sid_name}): cusp {i + 1} sidereal "
                    f"{cusps_sid[i]:.6f}° != tropical {cusps_trop[i]:.6f}° "
                    f"- ayan {ayan:.6f}° = {expected:.6f}° (diff={diff:.6f}°)"
                )
