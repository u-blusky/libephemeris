"""
Tests for orbit_max_min_true_distance function.

This function calculates the minimum and maximum geocentric distances
during a planet's orbit, useful for determining perigee and apogee.
"""

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
)


class TestOrbitMaxMinTrueDistanceBasic:
    """Test basic functionality of orbit_max_min_true_distance."""

    @pytest.mark.unit
    def test_returns_tuple_of_three_floats(self):
        """orbit_max_min_true_distance should return a tuple of three floats."""
        jd = 2451545.0  # J2000
        result = ephem.orbit_max_min_true_distance(jd, SE_MARS, 0)

        assert isinstance(result, tuple)
        assert len(result) == 3
        assert isinstance(result[0], float)
        assert isinstance(result[1], float)
        assert isinstance(result[2], float)

    @pytest.mark.unit
    def test_min_less_than_max(self):
        """Minimum distance should be less than maximum distance."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, SE_MARS, 0
        )

        assert min_dist < max_dist

    @pytest.mark.unit
    def test_sun_returns_zeros_for_max_min(self):
        """Sun should return 0.0 for max/min as it has no geocentric orbit."""
        jd = 2451545.0
        result = ephem.orbit_max_min_true_distance(jd, SE_SUN, 0)

        assert result[0] == 0.0  # max
        assert result[1] == 0.0  # min
        # true_dist is the current Earth-Sun distance (~1 AU)
        assert 0.98 < result[2] < 1.02

    @pytest.mark.unit
    def test_earth_returns_zeros_for_max_min(self):
        """Earth should return 0.0 for max/min as it is the observer."""
        jd = 2451545.0
        result = ephem.orbit_max_min_true_distance(jd, SE_EARTH, 0)

        assert result[0] == 0.0  # max
        assert result[1] == 0.0  # min


class TestOrbitMaxMinTrueDistanceMoon:
    """Test orbit_max_min_true_distance for Moon."""

    @pytest.mark.unit
    def test_moon_min_is_perigee(self):
        """Moon minimum distance should be approximately perigee (~0.0024 AU)."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, SE_MOON, 0
        )

        # Moon perigee is about 356,500 km = 0.00238 AU
        assert 0.0020 < min_dist < 0.0026

    @pytest.mark.unit
    def test_moon_max_is_apogee(self):
        """Moon maximum distance should be approximately apogee (~0.0027 AU)."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, SE_MOON, 0
        )

        # Moon apogee is about 406,700 km = 0.00272 AU
        assert 0.0025 < max_dist < 0.0030


class TestOrbitMaxMinTrueDistanceInnerPlanets:
    """Test orbit_max_min_true_distance for inner planets."""

    @pytest.mark.unit
    def test_mercury_distance_range(self):
        """Mercury min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, SE_MERCURY, 0
        )

        # Mercury: min ~0.55 AU (inferior conjunction), max ~1.45 AU (superior conjunction)
        assert 0.4 < min_dist < 0.8
        assert 1.2 < max_dist < 1.6

    @pytest.mark.unit
    def test_venus_distance_range(self):
        """Venus min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, SE_VENUS, 0
        )

        # Venus: min ~0.26 AU (inferior conjunction), max ~1.74 AU (superior conjunction)
        assert 0.2 < min_dist < 0.4
        assert 1.6 < max_dist < 1.9


class TestOrbitMaxMinTrueDistanceOuterPlanets:
    """Test orbit_max_min_true_distance for outer planets."""

    @pytest.mark.unit
    def test_mars_distance_range(self):
        """Mars min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, SE_MARS, 0
        )

        # Mars: min ~0.37 AU (opposition at perihelion), max ~2.68 AU (conjunction at aphelion)
        assert 0.3 < min_dist < 0.6
        assert 2.4 < max_dist < 2.9

    @pytest.mark.unit
    def test_jupiter_distance_range(self):
        """Jupiter min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, SE_JUPITER, 0
        )

        # Jupiter: min ~4.0 AU (opposition), max ~6.5 AU (conjunction)
        assert 3.8 < min_dist < 4.5
        assert 6.0 < max_dist < 6.8

    @pytest.mark.unit
    def test_saturn_distance_range(self):
        """Saturn min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, SE_SATURN, 0
        )

        # Saturn: min ~8.0 AU, max ~11.1 AU
        assert 7.5 < min_dist < 9.0
        assert 10.0 < max_dist < 12.0

    @pytest.mark.unit
    def test_uranus_distance_range(self):
        """Uranus min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, SE_URANUS, 0
        )

        # Uranus: min ~17.3 AU, max ~21.1 AU
        assert 16.5 < min_dist < 19.0
        assert 20.0 < max_dist < 22.0

    @pytest.mark.unit
    def test_neptune_distance_range(self):
        """Neptune min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, SE_NEPTUNE, 0
        )

        # Neptune: min ~28.8 AU, max ~31.3 AU
        assert 28.0 < min_dist < 30.5
        assert 30.0 < max_dist < 32.5

    @pytest.mark.unit
    def test_pluto_distance_range(self):
        """Pluto min/max distances should be reasonable."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, SE_PLUTO, 0
        )

        # Pluto: min ~28.8 AU (at perihelion opposition), max ~50.3 AU (at aphelion conjunction)
        assert 28.0 < min_dist < 32.0
        assert 48.0 < max_dist < 52.0


class TestOrbitMaxMinTrueDistanceAliases:
    """Test that function aliases work."""

    @pytest.mark.unit
    def test_swe_orbit_max_min_true_distance_alias(self):
        """swe_orbit_max_min_true_distance should be available."""
        result = ephem.swe_orbit_max_min_true_distance(2451545.0, SE_MARS, 0)
        assert len(result) == 3

    @pytest.mark.unit
    def test_orbit_max_min_true_distance_alias(self):
        """orbit_max_min_true_distance should be available."""
        result = ephem.orbit_max_min_true_distance(2451545.0, SE_MARS, 0)
        assert len(result) == 3


class TestOrbitMaxMinTrueDistanceAllPlanets:
    """Test orbit_max_min_true_distance for all planets."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet_id,planet_name",
        [
            (SE_MERCURY, "Mercury"),
            (SE_VENUS, "Venus"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
            (SE_SATURN, "Saturn"),
            (SE_URANUS, "Uranus"),
            (SE_NEPTUNE, "Neptune"),
            (SE_PLUTO, "Pluto"),
        ],
    )
    def test_all_planets_return_positive_distances(self, planet_id, planet_name):
        """All planets should return positive min and max distances."""
        jd = 2451545.0
        max_dist, min_dist, true_dist = ephem.orbit_max_min_true_distance(
            jd, planet_id, 0
        )

        assert min_dist > 0, f"{planet_name} min distance should be positive"
        assert max_dist > 0, f"{planet_name} max distance should be positive"
        assert min_dist < max_dist, f"{planet_name} min should be less than max"
        assert true_dist > 0, f"{planet_name} true distance should be positive"
        assert min_dist <= true_dist <= max_dist, (
            f"{planet_name} true distance should be between min and max"
        )
