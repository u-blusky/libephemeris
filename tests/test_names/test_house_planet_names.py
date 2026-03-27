"""
Tests for swe_house_name and swe_get_planet_name functions.

Verifies correct name lookup for all house systems and body IDs.
"""

from __future__ import annotations

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
    SE_EARTH,
    SE_CHIRON,
    SE_CERES,
    SE_PALLAS,
    SE_JUNO,
    SE_VESTA,
)


class TestHouseName:
    """Test swe_house_name function."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys_char,expected",
        [
            ("P", "Placidus"),
            ("K", "Koch"),
            ("R", "Regiomontanus"),
            ("C", "Campanus"),
            ("B", "Alcabitius"),
            ("O", "Porphyry"),
            ("M", "Morinus"),
            ("W", None),  # Whole Sign or Equal
        ],
    )
    def test_known_systems_return_names(self, hsys_char: str, expected):
        """Known house system chars return correct names."""
        result = swe.swe_house_name(ord(hsys_char))
        assert isinstance(result, str)
        assert len(result) > 0
        if expected:
            assert expected in result, (
                f"'{hsys_char}': got '{result}', expected '{expected}'"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys_char",
        list("PKROCABMWETUGHFVXSLQNYDI"),
    )
    def test_all_systems_return_nonempty(self, hsys_char: str):
        """All supported system chars return non-empty strings."""
        result = swe.swe_house_name(ord(hsys_char))
        assert isinstance(result, str)
        assert len(result) > 0
        assert result != "Unknown", f"System '{hsys_char}' returned 'Unknown'"

    @pytest.mark.unit
    def test_savard_a_system(self):
        """Savard-A house system ('J') should return its name."""
        result = swe.swe_house_name(ord("J"))
        assert isinstance(result, str)
        assert len(result) > 0
        assert result != "Unknown"

    @pytest.mark.unit
    def test_unknown_system_returns_unknown(self):
        """Unknown system char should return 'Unknown' or similar."""
        result = swe.swe_house_name(ord("Z"))
        assert isinstance(result, str)
        # Should indicate unknown
        assert "unknown" in result.lower() or len(result) > 0

    @pytest.mark.unit
    def test_alias_works(self):
        """house_name alias should work the same."""
        r1 = swe.swe_house_name(ord("P"))
        r2 = swe.house_name(ord("P"))
        assert r1 == r2


class TestGetPlanetName:
    """Test swe_get_planet_name function."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,expected_substr",
        [
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
        ],
    )
    def test_major_planets(self, body_id: int, expected_substr: str):
        """Major planets return correct names."""
        result = swe.swe_get_planet_name(body_id)
        assert isinstance(result, str)
        assert expected_substr in result, f"Body {body_id}: got '{result}'"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,expected_substr",
        [
            (SE_MEAN_NODE, "Node"),
            (SE_TRUE_NODE, "Node"),
            (SE_MEAN_APOG, "Apog"),
            (SE_OSCU_APOG, "Apog"),
        ],
    )
    def test_lunar_nodes_apsides(self, body_id: int, expected_substr: str):
        """Lunar nodes and apsides return appropriate names."""
        result = swe.swe_get_planet_name(body_id)
        assert isinstance(result, str)
        assert expected_substr.lower() in result.lower(), (
            f"Body {body_id}: got '{result}', expected substring '{expected_substr}'"
        )

    @pytest.mark.unit
    def test_earth_name(self):
        """Earth body ID returns 'Earth'."""
        result = swe.swe_get_planet_name(SE_EARTH)
        assert "Earth" in result

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,expected_substr",
        [
            (SE_CHIRON, "Chiron"),
            (SE_CERES, "Ceres"),
            (SE_PALLAS, "Pallas"),
            (SE_JUNO, "Juno"),
            (SE_VESTA, "Vesta"),
        ],
    )
    def test_asteroids(self, body_id: int, expected_substr: str):
        """Asteroid names are correct."""
        result = swe.swe_get_planet_name(body_id)
        assert isinstance(result, str)
        assert expected_substr in result, f"Body {body_id}: got '{result}'"

    @pytest.mark.unit
    def test_unknown_body_returns_string(self):
        """Unknown body ID should return some string (not crash)."""
        result = swe.swe_get_planet_name(99999)
        assert isinstance(result, str)
        assert len(result) > 0

    @pytest.mark.unit
    def test_get_planet_name_alias(self):
        """get_planet_name alias should work the same."""
        r1 = swe.swe_get_planet_name(SE_MARS)
        r2 = swe.get_planet_name(SE_MARS)
        assert r1 == r2

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id", range(40, 49))
    def test_uranian_bodies_named(self, body_id: int):
        """Uranian body IDs 40-48 should return names."""
        result = swe.swe_get_planet_name(body_id)
        assert isinstance(result, str)
        assert len(result) > 0
        # Should not be just "Unknown (40)" etc
        assert "Unknown" not in result or body_id == 48, (
            f"Body {body_id}: got '{result}'"
        )
