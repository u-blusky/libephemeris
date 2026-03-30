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

ECLIPTIC_DISTANCE_BODIES = [
    (SE_TRUE_NODE, "TrueNode"),
    (SE_OSCU_APOG, "OscuApog"),
    (SE_INTP_APOG, "IntpApog"),
    (SE_INTP_PERG, "IntpPerg"),
]

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


class TestExtGeocentricPrecision:
    """Geocentric distance and distance-speed for all ICRS planets."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", DISTANCE_BODIES)
    def test_geocentric_distance_and_speed(
        self,
        compare: CompareHelper,
        ext_dates_200: list[float],
        body_id: int,
        body_name: str,
    ):
        worst_dist = (0.0, 0.0)
        worst_speed = (0.0, 0.0)

        for jd in ext_dates_200:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err_dist = abs(ref[2] - leb[2])
            if err_dist > worst_dist[0]:
                worst_dist = (err_dist, jd)

            err_speed = abs(ref[5] - leb[5])
            if err_speed > worst_speed[0]:
                worst_speed = (err_speed, jd)

        e_dist, jd_dist = worst_dist
        assert e_dist < TOLS_EXT.DISTANCE_AU, (
            f"{body_name}: max geocentric dist error = {e_dist:.2e} AU "
            f"at JD {jd_dist:.1f}"
        )

        e_speed, jd_speed = worst_speed
        assert e_speed < TOLS_EXT.SPEED_DIST_AU_DAY, (
            f"{body_name}: max dist speed error = {e_speed:.2e} AU/day "
            f"at JD {jd_speed:.1f}"
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
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_200:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        dist_tol = 0.01
        assert max_err < dist_tol, (
            f"{body_name}: max distance error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )


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

        helio_tol = TOLS_EXT.DISTANCE_AU * 15.0
        assert max_err < helio_tol, (
            f"{body_name}: max helio dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )
