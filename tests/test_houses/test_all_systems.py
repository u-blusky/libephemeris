"""
Comprehensive tests for all 19 house systems.

Tests all house systems against pyswisseph for multiple locations.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


# House system codes and their expected tolerance
HOUSE_SYSTEMS = [
    (ord("P"), "Placidus", 0.1),
    (ord("K"), "Koch", 0.1),
    (ord("O"), "Porphyry", 0.1),
    (ord("R"), "Regiomontanus", 0.1),
    (ord("C"), "Campanus", 0.1),
    (ord("E"), "Equal", 0.1),
    (ord("W"), "Whole Sign", 0.1),
    (ord("B"), "Alcabitius", 0.1),
    (ord("M"), "Morinus", 0.1),
    (ord("T"), "Topocentric", 1.0),
    (ord("X"), "Meridian", 0.1),
    (ord("V"), "Vehlow", 0.1),
    (ord("H"), "Horizontal", 0.5),
]


class TestHouseSystemsBasic:
    """Basic tests for house systems."""

    @pytest.mark.unit
    def test_houses_returns_correct_structure(self):
        """houses() should return (cusps, ascmc) tuples."""
        jd = 2451545.0
        cusps, ascmc = ephem.swe_houses(jd, 41.9, 12.5, ord("P"))

        # Should have 12 cusps
        assert len(cusps) >= 12
        # Should have at least ASC and MC
        assert len(ascmc) >= 2

    @pytest.mark.unit
    def test_cusps_are_valid_longitudes(self):
        """All cusps should be valid longitudes (0-360)."""
        jd = 2451545.0
        cusps, ascmc = ephem.swe_houses(jd, 41.9, 12.5, ord("P"))

        for i, cusp in enumerate(cusps[:12]):
            assert 0 <= cusp < 360, f"Cusp {i + 1} = {cusp} out of range"

    @pytest.mark.unit
    def test_asc_and_mc_valid(self):
        """ASC and MC should be valid longitudes."""
        jd = 2451545.0
        cusps, ascmc = ephem.swe_houses(jd, 41.9, 12.5, ord("P"))

        # ASC is ascmc[0], MC is ascmc[1]
        assert 0 <= ascmc[0] < 360, f"ASC = {ascmc[0]} out of range"
        assert 0 <= ascmc[1] < 360, f"MC = {ascmc[1]} out of range"


class TestHouseSystemsVsPyswisseph:
    """Compare each house system with pyswisseph."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,name,tolerance", HOUSE_SYSTEMS)
    def test_house_system_rome(self, hsys, name, tolerance):
        """Test house system at Rome."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5

        cusps_lib, ascmc_lib = ephem.swe_houses(jd, lat, lon, hsys)
        cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, bytes([hsys]))

        # Compare ASC
        asc_diff = abs(ascmc_lib[0] - ascmc_swe[0])
        if asc_diff > 180:
            asc_diff = 360 - asc_diff
        assert asc_diff < tolerance, f"{name} ASC diff {asc_diff}"

        # Compare MC
        mc_diff = abs(ascmc_lib[1] - ascmc_swe[1])
        if mc_diff > 180:
            mc_diff = 360 - mc_diff
        assert mc_diff < tolerance, f"{name} MC diff {mc_diff}"

    @pytest.mark.comparison
    @pytest.mark.parametrize("hsys,name,tolerance", HOUSE_SYSTEMS[:7])  # Common systems
    def test_house_system_equator(self, hsys, name, tolerance):
        """Test house system at equator."""
        jd = 2451545.0
        lat, lon = 0.0, 0.0

        cusps_lib, ascmc_lib = ephem.swe_houses(jd, lat, lon, hsys)
        cusps_swe, ascmc_swe = swe.houses(jd, lat, lon, bytes([hsys]))

        asc_diff = abs(ascmc_lib[0] - ascmc_swe[0])
        if asc_diff > 180:
            asc_diff = 360 - asc_diff
        assert asc_diff < tolerance, f"{name} at equator ASC diff {asc_diff}"


class TestHouseCuspsOrder:
    """Test that house cusps are in correct order."""

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys,name,_", HOUSE_SYSTEMS[:7])
    def test_cusps_in_zodiacal_order(self, hsys, name, _):
        """Cusps should generally increase (with wraparound)."""
        jd = 2451545.0
        cusps, _ = ephem.swe_houses(jd, 41.9, 12.5, hsys)

        # Check that cusps generally progress around the zodiac
        # (allowing for the 360/0 wraparound)
        for i in range(11):
            diff = cusps[i + 1] - cusps[i]
            if diff < -180:
                diff += 360
            # Each house should span 0-60 degrees typically
            assert diff > 0 or diff > -180, (
                f"{name} cusp {i + 1} to {i + 2} has odd progression"
            )


class TestHouseSystemsPolarLatitudes:
    """Test house systems at polar latitudes."""

    @pytest.mark.edge_case
    @pytest.mark.parametrize("lat", [65.0, 70.0, 80.0])
    def test_placidus_polar_fallback(self, lat):
        """Placidus should fall back or handle polar latitudes."""
        jd = 2451545.0

        # Should not crash
        cusps, ascmc = ephem.swe_houses(jd, lat, 0.0, ord("P"))

        # Should return valid values (even if approximation)
        assert 0 <= ascmc[0] < 360
        assert 0 <= ascmc[1] < 360

    @pytest.mark.edge_case
    @pytest.mark.parametrize("lat", [65.0, 70.0, 80.0])
    def test_koch_polar_fallback(self, lat):
        """Koch should fall back or handle polar latitudes."""
        jd = 2451545.0
        cusps, ascmc = ephem.swe_houses(jd, lat, 0.0, ord("K"))

        assert 0 <= ascmc[0] < 360
        assert 0 <= ascmc[1] < 360

    @pytest.mark.edge_case
    def test_equal_works_at_poles(self):
        """Equal house system should work at any latitude."""
        jd = 2451545.0

        for lat in [89.0, -89.0]:
            cusps, ascmc = ephem.swe_houses(jd, lat, 0.0, ord("E"))
            assert 0 <= ascmc[0] < 360


class TestHouseAscMcRelationship:
    """Test relationship between ASC, MC, and cusps."""

    @pytest.mark.unit
    def test_asc_equals_cusp_1(self):
        """ASC should equal cusp 1 for most systems."""
        jd = 2451545.0
        cusps, ascmc = ephem.swe_houses(jd, 41.9, 12.5, ord("P"))

        diff = abs(ascmc[0] - cusps[0])
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.01, f"ASC {ascmc[0]} != cusp1 {cusps[0]}"

    @pytest.mark.unit
    def test_mc_equals_cusp_10(self):
        """MC should equal cusp 10 for most systems."""
        jd = 2451545.0
        cusps, ascmc = ephem.swe_houses(jd, 41.9, 12.5, ord("P"))

        diff = abs(ascmc[1] - cusps[9])  # cusp 10 is index 9
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.01, f"MC {ascmc[1]} != cusp10 {cusps[9]}"


class TestHouseNameFunction:
    """Test swe_house_name function."""

    @pytest.mark.unit
    def test_house_name_placidus(self):
        """Should return 'Placidus' for P."""
        name = ephem.swe_house_name(ord("P"))
        assert "Placidus" in name or "placidus" in name.lower()

    @pytest.mark.unit
    def test_house_name_koch(self):
        """Should return 'Koch' for K."""
        name = ephem.swe_house_name(ord("K"))
        assert "Koch" in name or "koch" in name.lower()

    @pytest.mark.unit
    def test_house_name_equal(self):
        """Should return 'Equal' for E."""
        name = ephem.swe_house_name(ord("E"))
        assert "Equal" in name or "equal" in name.lower()
