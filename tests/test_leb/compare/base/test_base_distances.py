"""
LEB vs Skyfield Comparison: Distance Precision (Base Tier).

Dedicated distance validation for geocentric and heliocentric bodies
across the base tier range (1860-2140).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    SEFLG_SPEED,
    SEFLG_HELCTR,
    SE_MARS,
    SE_JUPITER,
    SE_SATURN,
    SE_URANUS,
    SE_NEPTUNE,
)

from tests.test_leb.compare.conftest import (
    ICRS_PLANETS,
    CompareHelper,
)

from .conftest import TOLS_BASE

DISTANCE_BODIES = ICRS_PLANETS + [(15, "Chiron"), (17, "Ceres")]
HELIO_BODIES = [
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
    (SE_SATURN, "Saturn"),
    (SE_URANUS, "Uranus"),
    (SE_NEPTUNE, "Neptune"),
]


class TestBaseGeocentricDistance:
    """Geocentric distance precision for ICRS planets and asteroids."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", DISTANCE_BODIES)
    def test_geocentric_distance(
        self,
        compare: CompareHelper,
        base_dates_150: list[float],
        body_id: int,
        body_name: str,
    ):
        """Geocentric distance matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_150:
            try:
                ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
                leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            except (KeyError, ValueError):
                continue

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_BASE.DISTANCE_AU, (
            f"{body_name}: max dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )


class TestBaseHeliocentricDistance:
    """Heliocentric distance precision for outer planets."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", HELIO_BODIES)
    def test_heliocentric_distance(
        self,
        compare: CompareHelper,
        base_dates_150: list[float],
        body_id: int,
        body_name: str,
    ):
        """Heliocentric distance matches Skyfield within tolerance."""
        flags = SEFLG_SPEED | SEFLG_HELCTR
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_150:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_BASE.DISTANCE_AU, (
            f"{body_name} helio: max dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )


class TestBaseDistanceStatistics:
    """Statistical analysis of distance errors across base tier."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    def test_distance_error_distribution(
        self, compare: CompareHelper, base_dates_150: list[float]
    ):
        """Distance error distribution stays within tolerance for all planets."""
        errors: dict[str, list[float]] = {name: [] for _, name in ICRS_PLANETS}

        for jd in base_dates_150:
            for body_id, body_name in ICRS_PLANETS:
                ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
                leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

                err = abs(ref[2] - leb[2])
                errors[body_name].append(err)

        for name, err_list in errors.items():
            if err_list:
                assert max(err_list) < TOLS_BASE.DISTANCE_AU, (
                    f"{name}: max dist error {max(err_list):.2e} AU exceeds tolerance"
                )
