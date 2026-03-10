"""
LEB vs Skyfield Comparison: Ayanamsha Values.

Validates swe_get_ayanamsa_ut() LEB dispatch directly.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem

from .conftest import (
    TOLS,
    FORMULA_SIDEREAL_MODES,
    STAR_BASED_SIDEREAL_MODES,
    CompareHelper,
    year_to_jd,
    generate_test_dates,
)


@pytest.fixture(scope="module")
def ayanamsha_dates():
    """50 dates in 1560-2640 range for ayanamsha tests."""
    return generate_test_dates(50, year_to_jd(1560), year_to_jd(2640))


class TestAyanamshaValues:
    """Direct ayanamsha value comparison for formula-based modes."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("sid_mode", FORMULA_SIDEREAL_MODES)
    def test_ayanamsha_value(
        self, compare: CompareHelper, ayanamsha_dates: list[float], sid_mode: int
    ):
        """Ayanamsha value matches Skyfield within tolerance."""
        max_err = 0.0
        worst_jd = 0.0

        for jd in ayanamsha_dates:
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            ref = compare.skyfield(ephem.swe_get_ayanamsa_ut, jd)
            leb = compare.leb(ephem.swe_get_ayanamsa_ut, jd)

            err = abs(ref - leb) * 3600.0  # Convert to arcsec
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.AYANAMSHA_ARCSEC, (
            f'Ayanamsha mode {sid_mode}: max error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )


class TestAyanamshaConsistency:
    """Verify ayanamsha consistency with sidereal positions."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("sid_mode", [0, 1, 2, 3])
    def test_ayanamsha_vs_sidereal_offset(
        self, compare: CompareHelper, ayanamsha_dates: list[float], sid_mode: int
    ):
        """Tropical-sidereal difference matches true ayanamsha (with nutation).

        Sidereal positions subtract the *true* ayanamsha (mean + Δψ) so
        that nutation cancels from the ecliptic longitude. We use
        swe_get_ayanamsa_ex_ut() with flags=0 to get the true ayanamsha
        (with nutation included) for comparison.
        """
        from libephemeris.constants import SE_SUN, SEFLG_SPEED, SEFLG_SIDEREAL

        max_err = 0.0

        for jd in ayanamsha_dates[:10]:
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            # Get true ayanamsha (with nutation) via extended API
            _retflag, aya_true = ephem.swe_get_ayanamsa_ex_ut(jd, 0)

            # Get tropical and sidereal positions
            tropical, _ = compare.leb(ephem.swe_calc_ut, jd, SE_SUN, SEFLG_SPEED)
            sidereal, _ = compare.leb(
                ephem.swe_calc_ut, jd, SE_SUN, SEFLG_SPEED | SEFLG_SIDEREAL
            )

            # Difference should match true ayanamsha
            diff = tropical[0] - sidereal[0]
            if diff < 0:
                diff += 360.0

            err = abs(diff - aya_true) * 3600.0
            max_err = max(max_err, err)

        assert max_err < TOLS.AYANAMSHA_ARCSEC, (
            f'Mode {sid_mode}: ayanamsha vs position diff = {max_err:.4f}"'
        )


class TestStarBasedAyanamshaFallback:
    """Star-based modes should produce identical results."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("sid_mode", STAR_BASED_SIDEREAL_MODES[:4])
    def test_star_based_identical(
        self, compare: CompareHelper, ayanamsha_dates: list[float], sid_mode: int
    ):
        """Star-based ayanamsha produces identical results (both fallback)."""
        for jd in ayanamsha_dates[:5]:
            ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

            ref = compare.skyfield(ephem.swe_get_ayanamsa_ut, jd)
            leb = compare.leb(ephem.swe_get_ayanamsa_ut, jd)

            assert ref == pytest.approx(leb, rel=1e-10), (
                f"Star-based mode {sid_mode} differs at JD {jd}"
            )
