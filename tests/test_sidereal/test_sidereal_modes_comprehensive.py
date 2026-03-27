"""
Comprehensive tests for sidereal ayanamsha modes.

Tests all 43 ayanamsha modes for consistency, value ranges, and
correct application to planetary positions and house cusps.
"""

from __future__ import annotations

import math
import random

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SEFLG_SIDEREAL,
    SEFLG_SPEED,
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


ALL_MODES = [
    (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
    (SE_SIDM_LAHIRI, "Lahiri"),
    (SE_SIDM_DELUCE, "De Luce"),
    (SE_SIDM_RAMAN, "Raman"),
    (SE_SIDM_USHASHASHI, "Ushashashi"),
    (SE_SIDM_KRISHNAMURTI, "Krishnamurti"),
    (SE_SIDM_DJWHAL_KHUL, "Djwhal Khul"),
    (SE_SIDM_YUKTESHWAR, "Yukteshwar"),
    (SE_SIDM_JN_BHASIN, "JN Bhasin"),
    (SE_SIDM_BABYL_KUGLER1, "Babyl Kugler1"),
    (SE_SIDM_BABYL_KUGLER2, "Babyl Kugler2"),
    (SE_SIDM_BABYL_KUGLER3, "Babyl Kugler3"),
    (SE_SIDM_BABYL_HUBER, "Babyl Huber"),
    (SE_SIDM_BABYL_ETPSC, "Babyl ETPSC"),
    (SE_SIDM_ALDEBARAN_15TAU, "Aldebaran 15Tau"),
    (SE_SIDM_HIPPARCHOS, "Hipparchos"),
    (SE_SIDM_SASSANIAN, "Sassanian"),
    (SE_SIDM_GALCENT_0SAG, "Galcent 0Sag"),
    (SE_SIDM_J2000, "J2000"),
    (SE_SIDM_J1900, "J1900"),
    (SE_SIDM_B1950, "B1950"),
    (SE_SIDM_SURYASIDDHANTA, "Surya Siddhanta"),
    (SE_SIDM_SURYASIDDHANTA_MSUN, "Surya MSun"),
    (SE_SIDM_ARYABHATA, "Aryabhata"),
    (SE_SIDM_ARYABHATA_MSUN, "Aryabhata MSun"),
    (SE_SIDM_SS_REVATI, "SS Revati"),
    (SE_SIDM_SS_CITRA, "SS Citra"),
    (SE_SIDM_TRUE_CITRA, "True Citra"),
    (SE_SIDM_TRUE_REVATI, "True Revati"),
    (SE_SIDM_TRUE_PUSHYA, "True Pushya"),
    (SE_SIDM_GALCENT_RGILBRAND, "Galcent Rgilbrand"),
    (SE_SIDM_GALEQU_IAU1958, "Galequ IAU1958"),
    (SE_SIDM_GALEQU_TRUE, "Galequ True"),
    (SE_SIDM_GALEQU_MULA, "Galequ Mula"),
    (SE_SIDM_GALALIGN_MARDYKS, "Galalign Mardyks"),
    (SE_SIDM_TRUE_MULA, "True Mula"),
    (SE_SIDM_GALCENT_MULA_WILHELM, "Galcent Mula Wilhelm"),
    (SE_SIDM_ARYABHATA_522, "Aryabhata 522"),
    (SE_SIDM_BABYL_BRITTON, "Babyl Britton"),
    (SE_SIDM_TRUE_SHEORAN, "True Sheoran"),
    (SE_SIDM_GALCENT_COCHRANE, "Galcent Cochrane"),
    (SE_SIDM_GALEQU_FIORENZA, "Galequ Fiorenza"),
    (SE_SIDM_VALENS_MOON, "Valens Moon"),
]

# Modes known to have unusual ayanamsha values (e.g., mode 40 ~357°)
UNUSUAL_MODES = {40}


class TestAyanamshaValues:
    """Tests for ayanamsha value computation across all modes."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", ALL_MODES)
    def test_ayanamsha_returns_finite(self, mode: int, name: str):
        """Each ayanamsha mode returns a finite value."""
        swe.swe_set_sid_mode(mode)
        jd = 2451545.0
        ayan = swe.swe_get_ayanamsa_ut(jd)
        assert math.isfinite(ayan), f"{name} (mode {mode}): ayanamsha not finite"

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", ALL_MODES)
    def test_ayanamsha_in_range(self, mode: int, name: str):
        """Ayanamsha at J2000 should be in a reasonable range."""
        swe.swe_set_sid_mode(mode)
        jd = 2451545.0
        ayan = swe.swe_get_ayanamsa_ut(jd)
        # Most ayanamshas are 20-30° at J2000; some are very different
        assert 0 <= ayan < 360, (
            f"{name} (mode {mode}): ayanamsha {ayan} out of [0, 360)"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode,name",
        [
            (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
            (SE_SIDM_LAHIRI, "Lahiri"),
            (SE_SIDM_RAMAN, "Raman"),
            (SE_SIDM_KRISHNAMURTI, "Krishnamurti"),
        ],
    )
    def test_major_ayanamsha_approximate_values(self, mode: int, name: str):
        """Major ayanamshas should be approximately 20-25° at J2000."""
        swe.swe_set_sid_mode(mode)
        jd = 2451545.0
        ayan = swe.swe_get_ayanamsa_ut(jd)
        assert 20 < ayan < 30, f"{name}: ayanamsha {ayan}° not in expected 20-30° range"

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", ALL_MODES)
    def test_ayanamsha_increases_over_time(self, mode: int, name: str):
        """Ayanamsha should generally increase with time (precession)."""
        swe.swe_set_sid_mode(mode)
        jd1 = 2451545.0  # J2000
        jd2 = jd1 + 365.25 * 100  # +100 years
        ayan1 = swe.swe_get_ayanamsa_ut(jd1)
        ayan2 = swe.swe_get_ayanamsa_ut(jd2)
        diff = (ayan2 - ayan1) % 360
        if diff > 180:
            diff -= 360
        # Precession is ~1.4°/century. Allow some modes to be zero-like
        assert diff > -5.0, f"{name}: ayanamsha decreased by {diff}° over 100 years"


class TestSiderealPositions:
    """Tests for sidereal planetary positions."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", ALL_MODES)
    def test_sidereal_sun_valid(self, mode: int, name: str):
        """Sun sidereal position is valid for each mode."""
        swe.swe_set_sid_mode(mode)
        jd = 2451545.0
        result, _ = swe.swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)
        assert 0 <= result[0] < 360, (
            f"{name}: Sun sidereal lon {result[0]} out of range"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode,name",
        [
            (SE_SIDM_LAHIRI, "Lahiri"),
            (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
            (SE_SIDM_RAMAN, "Raman"),
            (SE_SIDM_TRUE_CITRA, "True Citra"),
            (SE_SIDM_GALCENT_0SAG, "Galcent 0Sag"),
            (SE_SIDM_J2000, "J2000"),
            (SE_SIDM_BABYL_KUGLER1, "Babyl Kugler1"),
            (SE_SIDM_HIPPARCHOS, "Hipparchos"),
            (SE_SIDM_YUKTESHWAR, "Yukteshwar"),
            (SE_SIDM_SURYASIDDHANTA, "Surya Siddhanta"),
        ],
    )
    def test_sidereal_offset_matches_ayanamsha(self, mode: int, name: str):
        """Tropical - sidereal should approximately equal ayanamsha."""
        swe.swe_set_sid_mode(mode)
        jd = 2451545.0
        r_trop, _ = swe.swe_calc_ut(jd, SE_SUN, 0)
        r_sid, _ = swe.swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)
        ayan = swe.swe_get_ayanamsa_ut(jd)

        diff = (r_trop[0] - r_sid[0]) % 360
        if diff > 350:
            diff -= 360
        ayan_norm = ayan % 360
        if ayan_norm > 350:
            ayan_norm -= 360

        assert abs(diff - ayan_norm) < 0.01, (
            f"{name}: tropical-sidereal diff {diff:.4f} vs ayanamsha {ayan_norm:.4f}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
            (SE_SATURN, "Saturn"),
        ],
    )
    def test_sidereal_lahiri_multiple_bodies(self, body_id: int, body_name: str):
        """Lahiri sidereal positions valid for multiple bodies."""
        swe.swe_set_sid_mode(SE_SIDM_LAHIRI)
        jd = 2451545.0
        result, _ = swe.swe_calc_ut(jd, body_id, SEFLG_SIDEREAL | SEFLG_SPEED)
        assert 0 <= result[0] < 360
        assert math.isfinite(result[3])

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", ALL_MODES)
    def test_sidereal_multiple_dates(self, mode: int, name: str):
        """Sidereal position valid across multiple dates for each mode."""
        swe.swe_set_sid_mode(mode)
        jds = [2415020.0, 2440587.5, 2451545.0, 2460000.0]
        for jd in jds:
            result, _ = swe.swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)
            assert 0 <= result[0] < 360, f"{name}: lon {result[0]} at JD {jd}"
            assert math.isfinite(result[1])

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode,name",
        [
            (SE_SIDM_LAHIRI, "Lahiri"),
            (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
        ],
    )
    def test_sidereal_speed_close_to_tropical(self, mode: int, name: str):
        """Sidereal speed should be close to tropical speed (differ by ~precession)."""
        swe.swe_set_sid_mode(mode)
        jd = 2451545.0
        r_trop, _ = swe.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)
        r_sid, _ = swe.swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL | SEFLG_SPEED)
        # Speed difference should be tiny (~0.004°/day for precession)
        assert abs(r_trop[3] - r_sid[3]) < 0.01, (
            f"{name}: tropical speed {r_trop[3]} vs sidereal {r_sid[3]}"
        )


class TestSiderealHouses:
    """Tests for sidereal house cusps."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode,name",
        [
            (SE_SIDM_LAHIRI, "Lahiri"),
            (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
            (SE_SIDM_RAMAN, "Raman"),
            (SE_SIDM_KRISHNAMURTI, "Krishnamurti"),
            (SE_SIDM_TRUE_CITRA, "True Citra"),
        ],
    )
    def test_sidereal_houses_valid(self, mode: int, name: str):
        """Sidereal houses return 12 valid cusps."""
        swe.swe_set_sid_mode(mode)
        jd = 2451545.0
        cusps, ascmc = swe.swe_houses_ex(jd, 41.9, 12.5, ord("P"), SEFLG_SIDEREAL)
        assert len(cusps) >= 12
        for i, cusp in enumerate(cusps[:12]):
            assert 0 <= cusp < 360, f"{name}: cusp {i + 1} = {cusp} out of range"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "mode,name",
        [
            (SE_SIDM_LAHIRI, "Lahiri"),
            (SE_SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
        ],
    )
    def test_sidereal_cusps_offset_from_tropical(self, mode: int, name: str):
        """Sidereal cusps should differ from tropical by ~ayanamsha."""
        swe.swe_set_sid_mode(mode)
        jd = 2451545.0
        lat, lon = 41.9, 12.5
        cusps_trop, _ = swe.swe_houses(jd, lat, lon, ord("P"))
        cusps_sid, _ = swe.swe_houses_ex(jd, lat, lon, ord("P"), SEFLG_SIDEREAL)
        ayan = swe.swe_get_ayanamsa_ut(jd)

        for i in range(12):
            diff = (cusps_trop[i] - cusps_sid[i]) % 360
            if diff > 350:
                diff -= 360
            ayan_norm = ayan % 360
            if ayan_norm > 350:
                ayan_norm -= 360
            assert abs(diff - ayan_norm) < 0.1, (
                f"{name}: cusp {i + 1} diff {diff:.3f} vs ayanamsha {ayan_norm:.3f}"
            )


class TestAyanamshaExUt:
    """Tests for the extended ayanamsha function."""

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", ALL_MODES)
    def test_get_ayanamsa_ex_ut_returns_correct_retflag(self, mode: int, name: str):
        """get_ayanamsa_ex_ut should return SEFLG_SWIEPH (2) as retflag."""
        swe.swe_set_sid_mode(mode)
        jd = 2451545.0
        retflag, ayan = swe.swe_get_ayanamsa_ex_ut(jd, 0)
        # BUG-004 fix: retflag should be 2 (SEFLG_SWIEPH), not the input flags
        assert retflag == 2, f"{name}: retflag {retflag} != 2 (SEFLG_SWIEPH)"
        assert math.isfinite(ayan)

    @pytest.mark.unit
    @pytest.mark.parametrize("mode,name", ALL_MODES)
    def test_ayanamsa_ex_ut_matches_ayanamsa_ut(self, mode: int, name: str):
        """Extended function should match simple function."""
        swe.swe_set_sid_mode(mode)
        jd = 2451545.0
        ayan_simple = swe.swe_get_ayanamsa_ut(jd)
        _, ayan_ex = swe.swe_get_ayanamsa_ex_ut(jd, 0)
        assert abs(ayan_simple - ayan_ex) < 1e-10, (
            f"{name}: simple {ayan_simple} vs ex {ayan_ex}"
        )
