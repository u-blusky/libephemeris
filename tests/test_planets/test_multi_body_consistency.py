"""Tests for multi-body consistency: positions should be physically consistent."""

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
    SE_OSCU_APOG,
    SE_CHIRON,
    SE_EARTH,
    SEFLG_SWIEPH,
    SEFLG_SPEED,
    SEFLG_HELCTR,
    SEFLG_EQUATORIAL,
)

JD_J2000 = 2451545.0


@pytest.mark.unit
class TestDistanceOrdering:
    """Verify planetary distances follow expected ordering."""

    def test_inner_planets_closer_than_outer(self):
        """Inner planets should generally be closer to Earth than outer planets."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        bodies = [SE_MERCURY, SE_VENUS, SE_MARS, SE_JUPITER, SE_SATURN]
        distances = []
        for body in bodies:
            result, _ = swe.calc_ut(JD_J2000, body, flags)
            distances.append(result[2])
        # Jupiter should be farther than Mars (generally true)
        mars_idx = bodies.index(SE_MARS)
        jupiter_idx = bodies.index(SE_JUPITER)
        saturn_idx = bodies.index(SE_SATURN)
        assert distances[jupiter_idx] > distances[mars_idx]
        assert distances[saturn_idx] > distances[jupiter_idx]

    def test_outer_planet_distance_ordering(self):
        """Outer planets should be in approximate distance order."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        bodies = [SE_JUPITER, SE_SATURN, SE_URANUS, SE_NEPTUNE]
        distances = []
        for body in bodies:
            result, _ = swe.calc_ut(JD_J2000, body, flags)
            distances.append(result[2])
        for i in range(1, len(distances)):
            assert distances[i] > distances[i - 1], (
                f"Distance ordering violated at index {i}"
            )

    def test_moon_closest(self):
        """Moon should be the closest body to Earth."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        moon_result, _ = swe.calc_ut(JD_J2000, SE_MOON, flags)
        moon_dist = moon_result[2]
        # Moon ~0.0026 AU
        assert moon_dist < 0.01, f"Moon distance {moon_dist} AU too large"
        # Compare to Mercury
        merc_result, _ = swe.calc_ut(JD_J2000, SE_MERCURY, flags)
        assert merc_result[2] > moon_dist

    def test_sun_distance_1au(self):
        """Sun distance should be ~1 AU."""
        result, _ = swe.calc_ut(JD_J2000, SE_SUN, SEFLG_SWIEPH | SEFLG_SPEED)
        sun_dist = result[2]
        assert 0.98 < sun_dist < 1.02, f"Sun distance {sun_dist} AU not ~1"


@pytest.mark.unit
class TestHeliocentricConsistency:
    """Heliocentric vs geocentric consistency."""

    def test_sun_helio_distance(self):
        """Sun in heliocentric mode should return distance 0 (Sun from itself)."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_HELCTR
        result, _ = swe.calc_ut(JD_J2000, SE_SUN, flags)
        # Sun heliocentric = Sun from Sun = zero distance
        assert result[2] == pytest.approx(0.0, abs=1e-6)

    def test_earth_helio_valid(self):
        """Earth heliocentric returns a valid position."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_HELCTR
        result, _ = swe.calc_ut(JD_J2000, SE_EARTH, flags)
        lon, lat, dist = result[0], result[1], result[2]
        assert 0 <= lon < 360
        assert -90 <= lat <= 90
        assert 0.98 < dist < 1.02  # ~1 AU from Sun

    def test_earth_helio_opposite_sun_geo(self):
        """Earth heliocentric longitude should be ~180° from geocentric Sun."""
        sun_geo, _ = swe.calc_ut(JD_J2000, SE_SUN, SEFLG_SWIEPH | SEFLG_SPEED)
        earth_helio, _ = swe.calc_ut(
            JD_J2000, SE_EARTH, SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_HELCTR
        )
        diff = abs(earth_helio[0] - sun_geo[0])
        if diff > 180:
            diff = 360 - diff
        assert diff == pytest.approx(180.0, abs=1.0)

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
    def test_helio_distance_reasonable(self, body, name):
        """Heliocentric distances should be in known ranges."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_HELCTR
        result, _ = swe.calc_ut(JD_J2000, body, flags)
        dist = result[2]
        expected_ranges = {
            SE_MERCURY: (0.3, 0.5),
            SE_VENUS: (0.7, 0.73),
            SE_MARS: (1.3, 1.7),
            SE_JUPITER: (4.9, 5.5),
            SE_SATURN: (9.0, 10.1),
        }
        lo, hi = expected_ranges[body]
        assert lo < dist < hi, f"{name} helio dist {dist} not in [{lo}, {hi}]"


@pytest.mark.unit
class TestSpeedConsistency:
    """Verify speed values are physically consistent."""

    def test_moon_fastest(self):
        """Moon should have the fastest longitude speed."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        moon, _ = swe.calc_ut(JD_J2000, SE_MOON, flags)
        sun, _ = swe.calc_ut(JD_J2000, SE_SUN, flags)
        assert abs(moon[3]) > abs(sun[3]), "Moon should be faster than Sun"

    def test_outer_planets_slower(self):
        """Outer planets should have slower apparent speed than inner."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        mars, _ = swe.calc_ut(JD_J2000, SE_MARS, flags)
        jupiter, _ = swe.calc_ut(JD_J2000, SE_JUPITER, flags)
        saturn, _ = swe.calc_ut(JD_J2000, SE_SATURN, flags)
        # Mean speeds: Mars ~0.5°/d, Jupiter ~0.08°/d, Saturn ~0.03°/d
        # Use absolute values since bodies can be retrograde
        # Check over many dates to avoid retrograde coincidences
        mars_speeds = []
        jup_speeds = []
        sat_speeds = []
        for offset in range(0, 360, 30):
            jd = JD_J2000 + offset
            m, _ = swe.calc_ut(jd, SE_MARS, flags)
            j, _ = swe.calc_ut(jd, SE_JUPITER, flags)
            s, _ = swe.calc_ut(jd, SE_SATURN, flags)
            mars_speeds.append(abs(m[3]))
            jup_speeds.append(abs(j[3]))
            sat_speeds.append(abs(s[3]))
        avg_mars = sum(mars_speeds) / len(mars_speeds)
        avg_jup = sum(jup_speeds) / len(jup_speeds)
        avg_sat = sum(sat_speeds) / len(sat_speeds)
        assert avg_mars > avg_jup > avg_sat

    def test_sun_speed_approximately_1_deg_per_day(self):
        """Sun's apparent speed should be ~1°/day."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        result, _ = swe.calc_ut(JD_J2000, SE_SUN, flags)
        sun_speed = result[3]
        assert 0.9 < abs(sun_speed) < 1.1, f"Sun speed {sun_speed}°/d not ~1"

    def test_moon_speed_approximately_13_deg_per_day(self):
        """Moon's apparent speed should be ~13°/day."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        result, _ = swe.calc_ut(JD_J2000, SE_MOON, flags)
        moon_speed = result[3]
        assert 11 < abs(moon_speed) < 15, f"Moon speed {moon_speed}°/d not ~13"

    def test_mean_node_always_retrograde(self):
        """Mean node should always have negative speed (retrograde)."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        for offset in range(0, 3650, 365):
            jd = JD_J2000 + offset
            result, _ = swe.calc_ut(jd, SE_MEAN_NODE, flags)
            assert result[3] < 0, (
                f"Mean node speed {result[3]} not retrograde at JD {jd}"
            )


@pytest.mark.unit
class TestEquatorialConsistency:
    """Verify equatorial coordinate transformations are consistent."""

    def test_equatorial_declination_range(self):
        """Declination should be in [-90, 90]."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_EQUATORIAL
        for body in [SE_SUN, SE_MOON, SE_MARS, SE_JUPITER]:
            result, _ = swe.calc_ut(JD_J2000, body, flags)
            dec = result[1]
            assert -90 <= dec <= 90, f"Dec {dec} out of range for body {body}"

    def test_equatorial_ra_range(self):
        """Right ascension should be in [0, 360)."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_EQUATORIAL
        for body in [SE_SUN, SE_MOON, SE_MARS, SE_JUPITER]:
            result, _ = swe.calc_ut(JD_J2000, body, flags)
            ra = result[0]
            assert 0 <= ra < 360, f"RA {ra} out of range for body {body}"

    def test_sun_max_declination(self):
        """Sun's declination should not exceed ~23.44° (obliquity)."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_EQUATORIAL
        max_dec = 0
        for offset in range(0, 365):
            jd = JD_J2000 + offset
            result, _ = swe.calc_ut(jd, SE_SUN, flags)
            max_dec = max(max_dec, abs(result[1]))
        assert 23 < max_dec < 24, f"Sun max dec {max_dec}° not near obliquity"

    def test_equatorial_distance_same_as_ecliptic(self):
        """Distance should be the same regardless of coordinate system."""
        for body in [SE_SUN, SE_MOON, SE_MARS]:
            ecl, _ = swe.calc_ut(JD_J2000, body, SEFLG_SWIEPH | SEFLG_SPEED)
            equ, _ = swe.calc_ut(
                JD_J2000, body, SEFLG_SWIEPH | SEFLG_SPEED | SEFLG_EQUATORIAL
            )
            assert ecl[2] == pytest.approx(equ[2], rel=1e-8), (
                f"Distance mismatch for body {body}"
            )


@pytest.mark.unit
class TestNodeConsistency:
    """Verify node positions are consistent."""

    def test_true_node_near_mean_node(self):
        """True node should be within ~2° of mean node."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        mean, _ = swe.calc_ut(JD_J2000, SE_MEAN_NODE, flags)
        true, _ = swe.calc_ut(JD_J2000, SE_TRUE_NODE, flags)
        diff = abs(mean[0] - true[0])
        if diff > 180:
            diff = 360 - diff
        assert diff < 3.0, f"True-mean node difference {diff}° too large"

    def test_node_latitude_near_zero(self):
        """Node latitude should be ~0 (on the ecliptic by definition)."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        mean, _ = swe.calc_ut(JD_J2000, SE_MEAN_NODE, flags)
        assert abs(mean[1]) < 0.01, f"Mean node lat {mean[1]}° not near zero"

    def test_mean_node_regression_rate(self):
        """Mean node regresses ~19.35° per year."""
        flags = SEFLG_SWIEPH | SEFLG_SPEED
        r1, _ = swe.calc_ut(JD_J2000, SE_MEAN_NODE, flags)
        r2, _ = swe.calc_ut(JD_J2000 + 365.25, SE_MEAN_NODE, flags)
        # Node moves ~19.35° per year retrograde
        diff = r1[0] - r2[0]  # Should be positive (retrograde)
        if diff < 0:
            diff += 360
        assert 18 < diff < 21, f"Node regression rate {diff}°/yr not ~19.35"
