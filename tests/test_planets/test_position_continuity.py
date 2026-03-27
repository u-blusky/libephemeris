"""Tests for position continuity — no discontinuities over time."""

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
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SE_MEAN_APOG,
    SE_CHIRON,
    SEFLG_SWIEPH,
    SEFLG_SPEED,
)

JD_J2000 = 2451545.0


def angle_diff(a: float, b: float) -> float:
    """Minimum angular difference between two angles in degrees."""
    d = (a - b) % 360
    if d > 180:
        d -= 360
    return abs(d)


@pytest.mark.unit
class TestSunContinuity:
    """Sun position should be continuous over time."""

    def test_sun_longitude_smooth_1day(self):
        """Sun longitude changes smoothly day-to-day (~1°/day)."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        prev, _ = swe.calc_ut(JD_J2000, SE_SUN, flags)
        for day in range(1, 365):
            jd = JD_J2000 + day
            curr, _ = swe.calc_ut(jd, SE_SUN, flags)
            diff = angle_diff(curr[0], prev[0])
            assert diff < 1.2, f"Sun jump {diff}° at day {day}"
            prev = curr

    def test_sun_latitude_near_zero(self):
        """Sun latitude should always be near zero (ecliptic)."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        for day in range(0, 365, 10):
            jd = JD_J2000 + day
            result, _ = swe.calc_ut(jd, SE_SUN, flags)
            assert abs(result[1]) < 0.01, f"Sun lat {result[1]}° at day {day}"


@pytest.mark.unit
class TestMoonContinuity:
    """Moon position should be continuous over time."""

    def test_moon_longitude_smooth_1hour(self):
        """Moon longitude changes smoothly hour-to-hour (~0.5°/hr)."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        prev, _ = swe.calc_ut(JD_J2000, SE_MOON, flags)
        for hour in range(1, 48):
            jd = JD_J2000 + hour / 24.0
            curr, _ = swe.calc_ut(jd, SE_MOON, flags)
            diff = angle_diff(curr[0], prev[0])
            assert diff < 0.7, f"Moon jump {diff}° at hour {hour}"
            prev = curr

    def test_moon_latitude_range(self):
        """Moon latitude should stay within ±5.3° (orbital inclination)."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        for day in range(0, 30):
            jd = JD_J2000 + day
            result, _ = swe.calc_ut(jd, SE_MOON, flags)
            assert abs(result[1]) < 5.4, f"Moon lat {result[1]}° at day {day}"

    def test_moon_distance_range(self):
        """Moon distance should stay within 0.0022-0.0028 AU."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        for day in range(0, 30):
            jd = JD_J2000 + day
            result, _ = swe.calc_ut(jd, SE_MOON, flags)
            dist = result[2]
            assert 0.0022 < dist < 0.0028, f"Moon dist {dist} AU at day {day}"


@pytest.mark.unit
class TestPlanetContinuity:
    """Planet positions should be continuous (even through retrograde)."""

    @pytest.mark.parametrize(
        "body,name,max_daily_change",
        [
            (SE_MERCURY, "Mercury", 2.5),
            (SE_VENUS, "Venus", 1.5),
            (SE_MARS, "Mars", 1.0),
            (SE_JUPITER, "Jupiter", 0.3),
            (SE_SATURN, "Saturn", 0.15),
        ],
    )
    def test_planet_longitude_smooth(self, body, name, max_daily_change):
        """Planet longitude changes smoothly day-to-day."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        prev, _ = swe.calc_ut(JD_J2000, body, flags)
        for day in range(1, 90):
            jd = JD_J2000 + day
            curr, _ = swe.calc_ut(jd, body, flags)
            diff = angle_diff(curr[0], prev[0])
            assert diff < max_daily_change, f"{name} jump {diff}° at day {day}"
            prev = curr

    @pytest.mark.parametrize(
        "body,name",
        [
            (SE_MERCURY, "Mercury"),
            (SE_VENUS, "Venus"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
            (SE_SATURN, "Saturn"),
        ],
    )
    def test_planet_latitude_bounded(self, body, name):
        """Planet latitude should be within expected bounds."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        max_lat = 0
        for day in range(0, 365, 5):
            jd = JD_J2000 + day
            result, _ = swe.calc_ut(jd, body, flags)
            max_lat = max(max_lat, abs(result[1]))
        # All planets should have latitude < 10° from ecliptic
        assert max_lat < 10, f"{name} max latitude {max_lat}°"


@pytest.mark.unit
class TestNodeContinuity:
    """Node positions should be continuous."""

    def test_mean_node_smooth(self):
        """Mean node longitude changes smoothly day-to-day."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        prev, _ = swe.calc_ut(JD_J2000, SE_MEAN_NODE, flags)
        for day in range(1, 365):
            jd = JD_J2000 + day
            curr, _ = swe.calc_ut(jd, SE_MEAN_NODE, flags)
            diff = angle_diff(curr[0], prev[0])
            assert diff < 0.1, f"Mean node jump {diff}° at day {day}"
            prev = curr

    def test_true_node_smooth(self):
        """True node longitude changes smoothly day-to-day."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        prev, _ = swe.calc_ut(JD_J2000, SE_TRUE_NODE, flags)
        for day in range(1, 365):
            jd = JD_J2000 + day
            curr, _ = swe.calc_ut(jd, SE_TRUE_NODE, flags)
            diff = angle_diff(curr[0], prev[0])
            # True node has short-period oscillations; allow larger change
            assert diff < 0.5, f"True node jump {diff}° at day {day}"
            prev = curr


@pytest.mark.unit
class TestSpeedMatchesNumericalDerivative:
    """Verify speed values match numerical derivatives."""

    @pytest.mark.parametrize(
        "body,name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
        ],
    )
    def test_speed_matches_numerical(self, body, name):
        """Speed from SEFLG_SPEED should match (pos(t+dt) - pos(t-dt)) / (2*dt)."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        dt = 0.01  # 0.01 days = ~14 minutes
        jd = JD_J2000
        result, _ = swe.calc_ut(jd, body, flags)
        speed_reported = result[3]  # longitude speed

        r_plus, _ = swe.calc_ut(jd + dt, body, SEFLG_SWIEPH)
        r_minus, _ = swe.calc_ut(jd - dt, body, SEFLG_SWIEPH)

        # Handle angle wrapping
        diff = r_plus[0] - r_minus[0]
        if diff > 180:
            diff -= 360
        elif diff < -180:
            diff += 360

        speed_numerical = diff / (2 * dt)
        # Tolerance depends on body — Moon has higher-order terms
        tol = 0.01 if body != SE_MOON else 0.05
        assert speed_reported == pytest.approx(speed_numerical, abs=tol), (
            f"{name} speed {speed_reported} vs numerical {speed_numerical}"
        )

    @pytest.mark.parametrize(
        "body,name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
            (SE_MARS, "Mars"),
        ],
    )
    def test_lat_speed_matches_numerical(self, body, name):
        """Latitude speed should match numerical derivative."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        dt = 0.01
        jd = JD_J2000
        result, _ = swe.calc_ut(jd, body, flags)
        speed_lat = result[4]

        r_plus, _ = swe.calc_ut(jd + dt, body, SEFLG_SWIEPH)
        r_minus, _ = swe.calc_ut(jd - dt, body, SEFLG_SWIEPH)

        speed_numerical = (r_plus[1] - r_minus[1]) / (2 * dt)
        tol = 0.01
        assert speed_lat == pytest.approx(speed_numerical, abs=tol), (
            f"{name} lat speed {speed_lat} vs numerical {speed_numerical}"
        )
