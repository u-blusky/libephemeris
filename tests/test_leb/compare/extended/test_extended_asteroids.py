"""
LEB vs Skyfield Comparison: Asteroid Positions (Extended Tier).

Validates ecliptic longitude, latitude, distance, and speed for all 5
asteroid bodies: Chiron, Ceres, Pallas, Juno, Vesta.

Asteroid SPK coverage is limited to ~1900-2100 CE.  Test dates are filtered
to this range to avoid catastrophically wrong Keplerian fallback data.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import SEFLG_SPEED

from tests.test_leb.compare.conftest import (
    ASTEROID_BODIES,
    CompareHelper,
    filter_asteroid_dates,
    lon_error_arcsec,
    generate_test_dates,
    year_to_jd,
)

from .conftest import TOLS_EXT


# Generate dates within SPK coverage (~1900-2100)
_SPK_START = year_to_jd(1910)
_SPK_END = year_to_jd(2090)


@pytest.fixture(scope="session")
def asteroid_dates() -> list[float]:
    """100 dates within asteroid SPK coverage (1910-2090)."""
    return generate_test_dates(100, _SPK_START, _SPK_END)


class TestExtAsteroidLongitude:
    """Asteroid longitude precision within SPK coverage."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_longitude(
        self,
        compare: CompareHelper,
        asteroid_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Asteroid ecliptic longitude matches Skyfield within tolerance."""
        dates = filter_asteroid_dates(asteroid_dates, body_id)
        if not dates:
            pytest.skip(f"No valid dates for {body_name}")

        max_err = 0.0
        worst_jd = 0.0

        for jd in dates:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.ASTEROID_ARCSEC, (
            f'{body_name}: max lon error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestExtAsteroidLatitude:
    """Asteroid latitude precision within SPK coverage."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_latitude(
        self,
        compare: CompareHelper,
        asteroid_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Asteroid ecliptic latitude matches Skyfield within tolerance."""
        dates = filter_asteroid_dates(asteroid_dates, body_id)
        if not dates:
            pytest.skip(f"No valid dates for {body_name}")

        max_err = 0.0
        worst_jd = 0.0

        for jd in dates:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[1] - leb[1]) * 3600.0
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.ASTEROID_ARCSEC, (
            f'{body_name}: max lat error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestExtAsteroidDistance:
    """Asteroid distance precision within SPK coverage."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_distance(
        self,
        compare: CompareHelper,
        asteroid_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Asteroid geocentric distance matches Skyfield within tolerance."""
        dates = filter_asteroid_dates(asteroid_dates, body_id)
        if not dates:
            pytest.skip(f"No valid dates for {body_name}")

        max_err = 0.0
        worst_jd = 0.0

        for jd in dates:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.DISTANCE_AU, (
            f"{body_name}: max dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )


class TestExtAsteroidSpeed:
    """Asteroid velocity precision within SPK coverage."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_speed_longitude(
        self,
        compare: CompareHelper,
        asteroid_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Asteroid longitude speed matches Skyfield within tolerance."""
        dates = filter_asteroid_dates(asteroid_dates, body_id)
        if not dates:
            pytest.skip(f"No valid dates for {body_name}")

        max_err = 0.0
        worst_jd = 0.0

        for jd in dates:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[3] - leb[3])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.ASTEROID_SPEED_LON_DEG_DAY, (
            f"{body_name}: max speed error = {max_err:.6f} deg/day at JD {worst_jd:.1f}"
        )

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_speed_latitude(
        self,
        compare: CompareHelper,
        asteroid_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        """Asteroid latitude speed matches Skyfield within tolerance."""
        dates = filter_asteroid_dates(asteroid_dates, body_id)
        if not dates:
            pytest.skip(f"No valid dates for {body_name}")

        max_err = 0.0

        for jd in dates:
            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, SEFLG_SPEED)

            err = abs(ref[4] - leb[4])
            max_err = max(max_err, err)

        assert max_err < TOLS_EXT.ASTEROID_SPEED_LAT_DEG_DAY, (
            f"{body_name}: max lat speed error = {max_err:.6f} deg/day"
        )
