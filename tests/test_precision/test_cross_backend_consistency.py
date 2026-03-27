"""
Cross-backend consistency tests.

Verifies that Skyfield and LEB backends produce consistent results
for all core bodies across multiple dates and flag combinations.
"""

from __future__ import annotations

import math
import random

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
    SE_OSCU_APOG,
    SEFLG_SPEED,
    SEFLG_EQUATORIAL,
    SEFLG_SIDEREAL,
    SEFLG_J2000,
    SEFLG_NOABERR,
    SE_SIDM_LAHIRI,
    SE_SIDM_FAGAN_BRADLEY,
)
from libephemeris.state import get_calc_mode


def _random_jds(n: int, seed: int = 42) -> list[float]:
    """Generate n random JDs in 1900-2100."""
    rng = random.Random(seed)
    return [rng.uniform(2415020.0, 2488070.0) for _ in range(n)]


def _angle_diff(a: float, b: float) -> float:
    """Absolute angular difference handling wrap-around."""
    d = abs(a - b) % 360
    return min(d, 360 - d)


# Bodies supported by both Skyfield and LEB
LEB_BODIES = [
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
    (SE_MEAN_NODE, "MeanNode"),
    (SE_TRUE_NODE, "TrueNode"),
    (SE_MEAN_APOG, "MeanApog"),
]


class TestSkyfieldVsLEB:
    """Compare Skyfield backend against LEB backend."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,body_name", LEB_BODIES)
    def test_skyfield_leb_agree_at_j2000(self, body_id: int, body_name: str):
        """Skyfield and LEB agree at J2000 within 1 arcsecond."""
        jd = 2451545.0

        swe.set_calc_mode("skyfield")
        r_sky, _ = swe.swe_calc_ut(jd, body_id, SEFLG_SPEED)

        swe.set_calc_mode("leb")
        r_leb, _ = swe.swe_calc_ut(jd, body_id, SEFLG_SPEED)

        # Longitude agreement
        lon_diff = _angle_diff(r_sky[0], r_leb[0])
        tol = 1 / 3600  # 1 arcsecond
        assert lon_diff < tol, (
            f"{body_name}: Skyfield lon {r_sky[0]:.6f} vs LEB {r_leb[0]:.6f} "
            f'diff {lon_diff * 3600:.3f}"'
        )

        # Latitude agreement
        lat_diff = abs(r_sky[1] - r_leb[1])
        assert lat_diff < tol, (
            f"{body_name}: Skyfield lat {r_sky[1]:.6f} vs LEB {r_leb[1]:.6f} "
            f'diff {lat_diff * 3600:.3f}"'
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
        ],
    )
    def test_skyfield_leb_agree_20_dates(self, body_id: int, body_name: str):
        """Skyfield and LEB agree at 20 random dates."""
        jds = _random_jds(20, seed=body_id * 41)
        tol = 2 / 3600  # 2 arcseconds

        for jd in jds:
            swe.set_calc_mode("skyfield")
            r_sky, _ = swe.swe_calc_ut(jd, body_id, SEFLG_SPEED)

            swe.set_calc_mode("leb")
            r_leb, _ = swe.swe_calc_ut(jd, body_id, SEFLG_SPEED)

            lon_diff = _angle_diff(r_sky[0], r_leb[0])
            assert lon_diff < tol, (
                f'{body_name} @ JD {jd:.1f}: lon diff {lon_diff * 3600:.3f}"'
            )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
            (SE_MARS, "Mars"),
        ],
    )
    def test_skyfield_leb_distance_agree(self, body_id: int, body_name: str):
        """Skyfield and LEB distances agree within 1e-7 AU."""
        jd = 2451545.0

        swe.set_calc_mode("skyfield")
        r_sky, _ = swe.swe_calc_ut(jd, body_id, 0)

        swe.set_calc_mode("leb")
        r_leb, _ = swe.swe_calc_ut(jd, body_id, 0)

        dist_diff = abs(r_sky[2] - r_leb[2])
        assert dist_diff < 1e-7, f"{body_name}: distance diff {dist_diff} AU"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
            (SE_MARS, "Mars"),
        ],
    )
    def test_skyfield_leb_speed_agree(self, body_id: int, body_name: str):
        """Skyfield and LEB speeds agree within 1e-4 °/day."""
        jd = 2451545.0

        swe.set_calc_mode("skyfield")
        r_sky, _ = swe.swe_calc_ut(jd, body_id, SEFLG_SPEED)

        swe.set_calc_mode("leb")
        r_leb, _ = swe.swe_calc_ut(jd, body_id, SEFLG_SPEED)

        speed_diff = abs(r_sky[3] - r_leb[3])
        # Moon speed can differ slightly more between backends (~1e-3)
        tol = 1e-3 if body_id == SE_MOON else 1e-4
        assert speed_diff < tol, f"{body_name}: speed diff {speed_diff}°/day"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags,desc",
        [
            (SEFLG_EQUATORIAL, "equatorial"),
            (SEFLG_J2000, "J2000"),
            (SEFLG_NOABERR, "no aberration"),
        ],
    )
    def test_skyfield_leb_flag_variants(self, flags: int, desc: str):
        """Skyfield and LEB agree with different flag variants."""
        jd = 2451545.0
        tol = 2 / 3600  # 2"

        swe.set_calc_mode("skyfield")
        r_sky, _ = swe.swe_calc_ut(jd, SE_MARS, flags)

        swe.set_calc_mode("leb")
        r_leb, _ = swe.swe_calc_ut(jd, SE_MARS, flags)

        lon_diff = _angle_diff(r_sky[0], r_leb[0])
        assert lon_diff < tol, f'Mars {desc}: lon diff {lon_diff * 3600:.3f}"'

    @pytest.mark.unit
    def test_skyfield_leb_sidereal_agree(self):
        """Skyfield and LEB agree in sidereal mode."""
        swe.swe_set_sid_mode(SE_SIDM_LAHIRI)
        jd = 2451545.0
        tol = 2 / 3600

        swe.set_calc_mode("skyfield")
        r_sky, _ = swe.swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)

        swe.set_calc_mode("leb")
        r_leb, _ = swe.swe_calc_ut(jd, SE_SUN, SEFLG_SIDEREAL)

        lon_diff = _angle_diff(r_sky[0], r_leb[0])
        assert lon_diff < tol, f'Sun sidereal: lon diff {lon_diff * 3600:.3f}"'
