"""
Tests for get_planet_name function.
"""

import pytest
from libephemeris import (
    get_planet_name,
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
    SE_EARTH,
)


class TestGetPlanetName:
    """Tests for the get_planet_name function."""

    def test_sun(self):
        """Test Sun name lookup."""
        assert get_planet_name(SE_SUN) == "Sun"
        assert get_planet_name(0) == "Sun"

    def test_moon(self):
        """Test Moon name lookup."""
        assert get_planet_name(SE_MOON) == "Moon"
        assert get_planet_name(1) == "Moon"

    def test_planets(self):
        """Test all main planet name lookups."""
        assert get_planet_name(SE_MERCURY) == "Mercury"
        assert get_planet_name(SE_VENUS) == "Venus"
        assert get_planet_name(SE_MARS) == "Mars"
        assert get_planet_name(SE_JUPITER) == "Jupiter"
        assert get_planet_name(SE_SATURN) == "Saturn"
        assert get_planet_name(SE_URANUS) == "Uranus"
        assert get_planet_name(SE_NEPTUNE) == "Neptune"
        assert get_planet_name(SE_PLUTO) == "Pluto"

    def test_earth(self):
        """Test Earth name lookup."""
        assert get_planet_name(SE_EARTH) == "Earth"

    def test_lunar_nodes(self):
        """Test lunar node name lookups."""
        assert get_planet_name(SE_MEAN_NODE) == "Mean Node"
        assert get_planet_name(SE_TRUE_NODE) == "True Node"

    def test_lunar_apogee(self):
        """Test lunar apogee name lookups."""
        assert get_planet_name(SE_MEAN_APOG) == "Mean Apogee"
        assert get_planet_name(SE_OSCU_APOG) == "Osculating Apogee"

    def test_unknown_planet_id(self):
        """Test unknown planet ID returns descriptive string."""
        result = get_planet_name(9999)
        assert result == "Unknown (9999)"

    def test_negative_planet_id(self):
        """Test negative planet ID returns descriptive string."""
        result = get_planet_name(-1)
        assert result == "Unknown (-1)"

    def test_return_type(self):
        """Test that get_planet_name always returns a string."""
        assert isinstance(get_planet_name(SE_SUN), str)
        assert isinstance(get_planet_name(SE_MOON), str)
        assert isinstance(get_planet_name(9999), str)
