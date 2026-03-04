"""
LEB vs Skyfield Comparison: Sidereal Mode (Base Tier).

Validates sidereal positions for formula-based ayanamsha modes
across the base tier range (1860-2140).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MARS,
    SE_JUPITER,
    SEFLG_SPEED,
    SEFLG_SIDEREAL,
)

from tests.test_leb.compare.conftest import (
    FORMULA_SIDEREAL_MODES,
    CompareHelper,
    generate_test_dates,
    lon_error_arcsec,
    year_to_jd,
)

from .conftest import TOLS_BASE

SIDEREAL_BODIES = [
    (SE_SUN, "Sun"),
    (SE_MOON, "Moon"),
    (SE_MARS, "Mars"),
    (SE_JUPITER, "Jupiter"),
]


@pytest.fixture(scope="module")
def base_sidereal_dates() -> list[float]:
    """50 dates in 1860-2140 range for sidereal tests."""
    return generate_test_dates(50, year_to_jd(1860), year_to_jd(2140))


class TestBaseSiderealPosition:
    """Sidereal position precision for formula-based modes."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", SIDEREAL_BODIES)
    @pytest.mark.parametrize("sid_mode", FORMULA_SIDEREAL_MODES)
    def test_sidereal_position(
        self,
        compare: CompareHelper,
        base_sidereal_dates: list[float],
        body_id: int,
        body_name: str,
        sid_mode: int,
    ):
        """Sidereal longitude matches Skyfield within tolerance."""
        flags = SEFLG_SPEED | SEFLG_SIDEREAL
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_sidereal_dates:
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            err = lon_error_arcsec(ref[0], leb[0])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_BASE.SIDEREAL_ARCSEC, (
            f'{body_name} sidereal mode {sid_mode}: max error = {max_err:.4f}" '
            f"at JD {worst_jd:.1f}"
        )


class TestBaseSiderealSpeed:
    """Sidereal speed precision (includes precession correction)."""

    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", SIDEREAL_BODIES)
    @pytest.mark.parametrize("sid_mode", [0, 1, 2, 3])
    def test_sidereal_speed(
        self,
        compare: CompareHelper,
        base_sidereal_dates: list[float],
        body_id: int,
        body_name: str,
        sid_mode: int,
    ):
        """Sidereal speed matches Skyfield within tolerance."""
        flags = SEFLG_SPEED | SEFLG_SIDEREAL
        max_err = 0.0

        for jd in base_sidereal_dates:
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            ref, _ = compare.skyfield(ephem.swe_calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.swe_calc_ut, jd, body_id, flags)

            err = abs(ref[3] - leb[3])
            max_err = max(max_err, err)

        assert max_err < TOLS_BASE.SPEED_LON_DEG_DAY, (
            f"{body_name} sidereal mode {sid_mode}: max speed error = "
            f"{max_err:.6f} deg/day"
        )
