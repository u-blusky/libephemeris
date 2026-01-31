"""
Tests for house_pos and swe_house_pos functions.

Tests the calculation of house position for celestial bodies.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


class TestHousePosBasic:
    """Basic tests for house_pos function."""

    @pytest.mark.unit
    def test_house_pos_returns_float(self):
        """house_pos() should return a float."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        lon = 15.0

        result = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)

        assert isinstance(result, float)

    @pytest.mark.unit
    def test_house_pos_in_valid_range(self):
        """house_pos() should return value between 1.0 and 13.0."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        lon = 15.0

        result = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)

        assert 1.0 <= result < 13.0

    @pytest.mark.unit
    def test_house_pos_integer_part_is_house_number(self):
        """Integer part of result should be house number (1-12)."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        lon = 15.0

        result = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)
        house_num = int(result)

        assert 1 <= house_num <= 12

    @pytest.mark.unit
    def test_house_pos_alias_matches(self):
        """swe_house_pos should give same result as house_pos."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        lon = 15.0

        result1 = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)
        result2 = ephem.swe_house_pos(armc, lat, eps, ord("P"), lon, 0.0)

        assert result1 == result2


class TestHousePosConsistency:
    """Test house_pos consistency with house cusps."""

    @pytest.mark.unit
    def test_body_at_cusp_start_has_near_zero_fraction(self):
        """Body exactly at cusp start should have near-zero fractional part."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        # Get the house cusps first
        cusps, ascmc = ephem.swe_houses_armc(armc, lat, eps, ord("P"))

        # Place body at the start of house 1 (the ASC)
        lon = cusps[0]

        result = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)
        house_num = int(result)
        fraction = result - house_num

        # Should be in house 1 with near-zero fraction
        assert house_num == 1
        assert fraction < 0.001, f"Expected fraction near 0, got {fraction}"

    @pytest.mark.unit
    def test_body_just_before_next_cusp(self):
        """Body just before next cusp should have fraction close to 1."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        # Get the house cusps
        cusps, ascmc = ephem.swe_houses_armc(armc, lat, eps, ord("P"))

        # Place body just before house 2 cusp (still in house 1)
        # Calculate a position that's 99% through house 1
        cusp1 = cusps[0]
        cusp2 = cusps[1]
        house_size = (cusp2 - cusp1 + 360.0) % 360.0
        lon = (cusp1 + house_size * 0.99) % 360.0

        result = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)
        house_num = int(result)
        fraction = result - house_num

        # Should be in house 1 with fraction close to 1
        assert house_num == 1
        assert fraction > 0.95, f"Expected fraction near 1, got {fraction}"


class TestHousePosMultipleSystems:
    """Test house_pos with different house systems."""

    HOUSE_SYSTEMS = ["P", "K", "R", "C", "E", "W", "O", "B"]

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys", HOUSE_SYSTEMS)
    def test_house_pos_all_systems_return_valid(self, hsys):
        """house_pos should return valid values for all major house systems."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        lon = 45.0  # 15 Taurus

        result = ephem.house_pos(armc, lat, eps, ord(hsys), lon, 0.0)

        assert 1.0 <= result < 13.0, f"Invalid result {result} for system {hsys}"


class TestHousePosComparisonWithSwisseph:
    """Compare house_pos results with Swiss Ephemeris."""

    # Test cases: (armc, lat, eps, lon, lat_body, hsys)
    TEST_CASES = [
        # Rome, Placidus, Sun at 15 Aries
        (292.957, 41.9, 23.4393, 15.0, 0.0, "P"),
        # London, Koch, Sun at 45 degrees
        (180.0, 51.5, 23.4393, 45.0, 0.0, "K"),
        # New York, Regiomontanus, Sun at 120 degrees
        (0.0, 40.7, 23.4393, 120.0, 0.0, "R"),
        # Equator, Placidus, Sun at 270 degrees
        (90.0, 0.0, 23.4393, 270.0, 0.0, "P"),
        # Southern hemisphere, Campanus
        (45.0, -33.9, 23.4393, 180.0, 0.0, "C"),
    ]

    @pytest.mark.unit
    @pytest.mark.parametrize("armc,lat,eps,lon,lat_body,hsys", TEST_CASES)
    def test_house_pos_matches_swisseph(self, armc, lat, eps, lon, lat_body, hsys):
        """house_pos should match Swiss Ephemeris swe_house_pos."""
        # Get result from libephemeris
        result_lib = ephem.house_pos(armc, lat, eps, ord(hsys), lon, lat_body)

        # Get result from Swiss Ephemeris
        # swe.house_pos(armc, geolat, eps, objcoord, hsys)
        result_swe = swe.house_pos(armc, lat, eps, (lon, lat_body), hsys.encode())

        # Compare house numbers
        house_lib = int(result_lib)
        house_swe = int(result_swe)

        assert house_lib == house_swe, (
            f"House mismatch: libephemeris={house_lib}, swisseph={house_swe}"
        )

        # Compare fractions with tolerance
        frac_lib = result_lib - house_lib
        frac_swe = result_swe - house_swe

        # Allow 6% tolerance for fractional position
        # (house calculations can have minor differences between implementations)
        assert abs(frac_lib - frac_swe) < 0.06, (
            f"Fraction mismatch: libephemeris={frac_lib:.4f}, swisseph={frac_swe:.4f}"
        )


class TestHousePosEdgeCases:
    """Test edge cases for house_pos."""

    @pytest.mark.unit
    def test_house_pos_longitude_at_360_normalized(self):
        """Longitude at 360 should be normalized to 0."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        result1 = ephem.house_pos(armc, lat, eps, ord("P"), 0.0, 0.0)
        result2 = ephem.house_pos(armc, lat, eps, ord("P"), 360.0, 0.0)

        assert result1 == result2

    @pytest.mark.unit
    def test_house_pos_negative_longitude(self):
        """Negative longitude should be normalized."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        # -30 should be same as 330
        result1 = ephem.house_pos(armc, lat, eps, ord("P"), -30.0, 0.0)
        result2 = ephem.house_pos(armc, lat, eps, ord("P"), 330.0, 0.0)

        assert result1 == result2

    @pytest.mark.unit
    def test_house_pos_with_ecliptic_latitude(self):
        """house_pos should accept non-zero ecliptic latitude."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393
        lon = 45.0

        # Should not raise exception for non-zero latitude
        result = ephem.house_pos(armc, lat, eps, ord("P"), lon, 5.0)

        assert 1.0 <= result < 13.0


class TestHousePosAllHouses:
    """Test that bodies can be found in all 12 houses."""

    @pytest.mark.unit
    def test_can_find_body_in_each_house(self):
        """Should be able to place a body in each of the 12 houses."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        # Get cusps
        cusps, ascmc = ephem.swe_houses_armc(armc, lat, eps, ord("P"))

        houses_found = set()

        for i in range(12):
            # Place body in middle of each house
            cusp_start = cusps[i]
            cusp_end = cusps[(i + 1) % 12]
            house_size = (cusp_end - cusp_start + 360.0) % 360.0
            lon = (cusp_start + house_size * 0.5) % 360.0

            result = ephem.house_pos(armc, lat, eps, ord("P"), lon, 0.0)
            house_num = int(result)
            houses_found.add(house_num)

        # Should have found all 12 houses
        assert houses_found == set(range(1, 13)), (
            f"Missing houses: {set(range(1, 13)) - houses_found}"
        )
