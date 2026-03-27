"""Tests for LEB vs Skyfield backend consistency."""

from __future__ import annotations

import math
import pytest
import libephemeris as swe
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
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SE_MEAN_APOG,
    SEFLG_SWIEPH,
    SEFLG_SPEED,
    SEFLG_SIDEREAL,
    SEFLG_HELCTR,
    SE_SIDM_LAHIRI,
)

JD_J2000 = 2451545.0


def reset_to_mode(mode: str):
    """Reset and switch to specified backend mode."""
    swe.swe_close()
    swe.set_calc_mode(mode)


@pytest.mark.unit
class TestLEBvsSkyfield:
    """Compare LEB and Skyfield backends for consistency."""

    @pytest.fixture(autouse=True)
    def _reset_after(self):
        yield
        swe.swe_close()

    @pytest.mark.parametrize(
        "body,name,tol_arcsec",
        [
            (SE_SUN, "Sun", 0.5),
            (SE_MOON, "Moon", 2.5),
            (SE_MERCURY, "Mercury", 1.0),
            (SE_VENUS, "Venus", 1.0),
            (SE_MARS, "Mars", 1.0),
            (SE_JUPITER, "Jupiter", 1.0),
            (SE_SATURN, "Saturn", 1.0),
        ],
    )
    def test_planet_longitude_matches(self, body, name, tol_arcsec):
        """Planet longitude should match between LEB and Skyfield backends."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED

        reset_to_mode("leb")
        leb_result, _ = swe.calc_ut(JD_J2000, body, flags)

        reset_to_mode("skyfield")
        sky_result, _ = swe.calc_ut(JD_J2000, body, flags)

        diff = abs(leb_result[0] - sky_result[0])
        if diff > 180:
            diff = 360 - diff
        diff_arcsec = diff * 3600
        assert diff_arcsec < tol_arcsec, (
            f'{name} lon diff: {diff_arcsec:.2f}" (LEB={leb_result[0]:.6f}, '
            f"Sky={sky_result[0]:.6f})"
        )

    @pytest.mark.parametrize(
        "body,name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
        ],
    )
    def test_planet_latitude_matches(self, body, name):
        """Planet latitude should match between LEB and Skyfield backends."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED

        reset_to_mode("leb")
        leb_result, _ = swe.calc_ut(JD_J2000, body, flags)

        reset_to_mode("skyfield")
        sky_result, _ = swe.calc_ut(JD_J2000, body, flags)

        diff_arcsec = abs(leb_result[1] - sky_result[1]) * 3600
        assert diff_arcsec < 2.5, f'{name} lat diff: {diff_arcsec:.2f}"'

    @pytest.mark.parametrize(
        "body,name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
            (SE_MARS, "Mars"),
        ],
    )
    def test_planet_distance_matches(self, body, name):
        """Planet distance should match between LEB and Skyfield backends."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED

        reset_to_mode("leb")
        leb_result, _ = swe.calc_ut(JD_J2000, body, flags)

        reset_to_mode("skyfield")
        sky_result, _ = swe.calc_ut(JD_J2000, body, flags)

        # Relative distance tolerance
        rel_diff = abs(leb_result[2] - sky_result[2]) / max(sky_result[2], 1e-10)
        assert rel_diff < 1e-4, f"{name} dist rel diff: {rel_diff:.2e}"

    @pytest.mark.parametrize(
        "body,name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
            (SE_MARS, "Mars"),
        ],
    )
    def test_speed_matches(self, body, name):
        """Speed should match between LEB and Skyfield backends."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED

        reset_to_mode("leb")
        leb_result, _ = swe.calc_ut(JD_J2000, body, flags)

        reset_to_mode("skyfield")
        sky_result, _ = swe.calc_ut(JD_J2000, body, flags)

        # Speed tolerance
        tol = 0.001 if body != SE_MOON else 0.01
        assert leb_result[3] == pytest.approx(sky_result[3], abs=tol), (
            f"{name} speed: LEB={leb_result[3]:.6f}, Sky={sky_result[3]:.6f}"
        )

    def test_nodes_match(self):
        """Node positions should match between backends."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED

        for body_id, name in [(SE_MEAN_NODE, "Mean Node"), (SE_TRUE_NODE, "True Node")]:
            reset_to_mode("leb")
            leb_result, _ = swe.calc_ut(JD_J2000, body_id, flags)

            reset_to_mode("skyfield")
            sky_result, _ = swe.calc_ut(JD_J2000, body_id, flags)

            diff = abs(leb_result[0] - sky_result[0])
            if diff > 180:
                diff = 360 - diff
            diff_arcsec = diff * 3600
            assert diff_arcsec < 2.0, f'{name} lon diff: {diff_arcsec:.2f}"'

    def test_sidereal_consistency(self):
        """Sidereal positions should be consistent across backends."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_SIDEREAL

        swe.swe_close()
        swe.set_calc_mode("leb")
        swe.swe_set_sid_mode(SE_SIDM_LAHIRI)
        leb_result, _ = swe.calc_ut(JD_J2000, SE_SUN, flags)

        swe.swe_close()
        swe.set_calc_mode("skyfield")
        swe.swe_set_sid_mode(SE_SIDM_LAHIRI)
        sky_result, _ = swe.calc_ut(JD_J2000, SE_SUN, flags)

        diff = abs(leb_result[0] - sky_result[0])
        if diff > 180:
            diff = 360 - diff
        diff_arcsec = diff * 3600
        assert diff_arcsec < 1.0, f'Sidereal Sun lon diff: {diff_arcsec:.2f}"'


@pytest.mark.unit
class TestLEBMultipleDates:
    """LEB vs Skyfield at multiple dates."""

    @pytest.fixture(autouse=True)
    def _reset_after(self):
        yield
        swe.swe_close()

    @pytest.mark.parametrize(
        "jd_offset",
        [0, 365, 3652, 7305, -3652, -7305],
        ids=["J2000", "+1yr", "+10yr", "+20yr", "-10yr", "-20yr"],
    )
    def test_sun_across_dates(self, jd_offset):
        """Sun position should match between backends at various dates."""
        jd = JD_J2000 + jd_offset
        flags = SEFLG_SWIEPH | SEFLG_SPEED

        reset_to_mode("leb")
        leb_result, _ = swe.calc_ut(jd, SE_SUN, flags)

        reset_to_mode("skyfield")
        sky_result, _ = swe.calc_ut(jd, SE_SUN, flags)

        diff = abs(leb_result[0] - sky_result[0])
        if diff > 180:
            diff = 360 - diff
        diff_arcsec = diff * 3600
        assert diff_arcsec < 0.5, (
            f'Sun lon diff at JD offset {jd_offset}: {diff_arcsec:.2f}"'
        )

    @pytest.mark.parametrize(
        "jd_offset",
        [0, 365, 3652, 7305],
        ids=["J2000", "+1yr", "+10yr", "+20yr"],
    )
    def test_moon_across_dates(self, jd_offset):
        """Moon position should match between backends at various dates."""
        jd = JD_J2000 + jd_offset
        flags = SEFLG_SWIEPH | SEFLG_SPEED

        reset_to_mode("leb")
        leb_result, _ = swe.calc_ut(jd, SE_MOON, flags)

        reset_to_mode("skyfield")
        sky_result, _ = swe.calc_ut(jd, SE_MOON, flags)

        diff = abs(leb_result[0] - sky_result[0])
        if diff > 180:
            diff = 360 - diff
        diff_arcsec = diff * 3600
        assert diff_arcsec < 2.5, (
            f'Moon lon diff at JD offset {jd_offset}: {diff_arcsec:.2f}"'
        )
