"""
LEB vs Skyfield Comparison: Distances (Extended Tier).

Validates geocentric and heliocentric distances for all body types
across the extended tier range (-5000 to 5000).  Fills the gap where
no dedicated distance tests existed for the extended tier.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MERCURY,
    SE_VENUS,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_URANUS,
    SE_NEPTUNE,
    SE_PLUTO,
    SE_EARTH,
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SE_MEAN_APOG,
    SE_OSCU_APOG,
    SE_INTP_APOG,
    SE_INTP_PERG,
    SEFLG_SPEED,
    SEFLG_HELCTR,
)

from tests.test_leb.compare.conftest import (
    ICRS_PLANETS,
    ECLIPTIC_BODIES,
    CompareHelper,
)

from .conftest import TOLS_EXT


# ============================================================================
# GEOCENTRIC DISTANCE
# ============================================================================

# Bodies with meaningful geocentric distance (exclude Sun which has ~1 AU)
DISTANCE_BODIES = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_MERCURY, "Mercury"),
    (SE_VENUS, "Venus"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
    (SE_URANUS, "Uranus"),
    (SE_NEPTUNE, "Neptune"),
    (SE_PLUTO, "Pluto"),
]

# Bodies with geocentric distance output
ECLIPTIC_DISTANCE_BODIES = [
    (SE_TRUE_NODE, "TrueNode"),
    (SE_OSCU_APOG, "OscuApog"),
    (SE_INTP_APOG, "IntpApog"),
    (SE_INTP_PERG, "IntpPerg"),
]

# Heliocentric test bodies (exclude Sun and Moon — Sun is center, Moon geocentric)
HELIO_BODIES = [
    (SE_MERCURY, "Mercury"),
    (SE_VENUS, "Venus"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
    (SE_URANUS, "Uranus"),
    (SE_NEPTUNE, "Neptune"),
    (SE_PLUTO, "Pluto"),
]


class TestExtGeocentricDistance:
    """Geocentric distance precision for all ICRS planets."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", DISTANCE_BODIES)
    def test_geocentric_distance(
        self,
        compare: CompareHelper,
        ext_dates_200: list[float],
        body_id: int,
        body_name: str,
    ):
        """Geocentric distance matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_200:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.DISTANCE_AU, (
            f"{body_name}: max geocentric dist error = {max_err:.2e} AU "
            f"at JD {worst_jd:.1f}"
        )


class TestExtEclipticBodyDistance:
    """Distance precision for ecliptic-direct bodies."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ECLIPTIC_DISTANCE_BODIES)
    def test_ecliptic_distance(
        self,
        compare: CompareHelper,
        ext_dates_200: list[float],
        body_id: int,
        body_name: str,
    ):
        """Ecliptic body distance matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_200:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        # Use a generous tolerance: some ecliptic bodies have non-physical
        # or convention-dependent distance values
        dist_tol = 0.01  # AU
        assert max_err < dist_tol, (
            f"{body_name}: max distance error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )


class TestExtDistanceSpeed:
    """Distance velocity (AU/day) precision."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", DISTANCE_BODIES)
    def test_distance_speed(
        self,
        compare: CompareHelper,
        ext_dates_200: list[float],
        body_id: int,
        body_name: str,
    ):
        """Distance velocity matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_200:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[5] - leb[5])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.SPEED_DIST_AU_DAY, (
            f"{body_name}: max dist speed error = {max_err:.2e} AU/day "
            f"at JD {worst_jd:.1f}"
        )


# ============================================================================
# HELIOCENTRIC DISTANCE
# ============================================================================


class TestExtHeliocentricDistance:
    """Heliocentric distance precision (SEFLG_HELCTR)."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", HELIO_BODIES)
    def test_heliocentric_distance(
        self,
        compare: CompareHelper,
        ext_dates_100: list[float],
        body_id: int,
        body_name: str,
    ):
        """Heliocentric distance matches Skyfield within tolerance."""
        flags = SEFLG_SPEED | SEFLG_HELCTR
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_100:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        # HELCTR amplifies Chebyshev error 8-11x for inner planets
        helio_tol = TOLS_EXT.DISTANCE_AU * 15.0
        assert max_err < helio_tol, (
            f"{body_name}: max helio dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )
