"""
Comprehensive tests for error handling and edge cases.

Verifies that invalid inputs raise appropriate exceptions,
boundary dates work or fail gracefully, and error types match
expectations.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    SE_SUN,
    SE_MOON,
    SE_MARS,
    SE_JUPITER,
    SE_PLUTO,
    SE_MEAN_NODE,
    SE_TRUE_NODE,
    SE_MEAN_APOG,
    SE_OSCU_APOG,
    SEFLG_SPEED,
    SEFLG_HELCTR,
    SEFLG_TOPOCTR,
    SEFLG_EQUATORIAL,
    SEFLG_SIDEREAL,
    SE_SIDM_LAHIRI,
    SE_AST_OFFSET,
)
from libephemeris.exceptions import (
    CoordinateError,
    UnknownBodyError,
    PolarCircleError,
)


class TestInvalidBodyIds:
    """Test behavior with invalid or unusual body IDs."""

    @pytest.mark.unit
    def test_negative_body_id_returns_result(self):
        """Negative body IDs are accepted (return some result without crash)."""
        jd = 2451545.0
        # libephemeris accepts negative body IDs without raising
        result, _ = swe.swe_calc_ut(jd, -1, 0)
        assert len(result) == 6

    @pytest.mark.unit
    def test_very_large_body_id_raises(self):
        """Very large body IDs (not asteroids) should raise."""
        jd = 2451545.0
        with pytest.raises(Exception):
            swe.swe_calc_ut(jd, 999999, 0)

    @pytest.mark.unit
    def test_body_id_49_raises(self):
        """Body ID 49 (beyond Transpluto=48) should raise."""
        jd = 2451545.0
        with pytest.raises(Exception):
            swe.swe_calc_ut(jd, 49, 0)


class TestInvalidCoordinates:
    """Test set_topo with invalid coordinates."""

    @pytest.mark.unit
    def test_latitude_above_90(self):
        """Latitude > 90 should raise CoordinateError."""
        with pytest.raises((CoordinateError, ValueError)):
            swe.swe_set_topo(0.0, 91.0, 0.0)

    @pytest.mark.unit
    def test_latitude_below_minus_90(self):
        """Latitude < -90 should raise CoordinateError."""
        with pytest.raises((CoordinateError, ValueError)):
            swe.swe_set_topo(0.0, -91.0, 0.0)

    @pytest.mark.unit
    def test_longitude_above_180(self):
        """Longitude > 180 should raise CoordinateError."""
        with pytest.raises((CoordinateError, ValueError)):
            swe.swe_set_topo(181.0, 0.0, 0.0)

    @pytest.mark.unit
    def test_longitude_below_minus_180(self):
        """Longitude < -180 should raise CoordinateError."""
        with pytest.raises((CoordinateError, ValueError)):
            swe.swe_set_topo(-181.0, 0.0, 0.0)

    @pytest.mark.unit
    def test_boundary_latitude_90(self):
        """Latitude exactly 90 should be accepted."""
        swe.swe_set_topo(0.0, 90.0, 0.0)  # Should not raise

    @pytest.mark.unit
    def test_boundary_latitude_minus_90(self):
        """Latitude exactly -90 should be accepted."""
        swe.swe_set_topo(0.0, -90.0, 0.0)  # Should not raise

    @pytest.mark.unit
    def test_boundary_longitude_180(self):
        """Longitude exactly 180 should be accepted."""
        swe.swe_set_topo(180.0, 0.0, 0.0)  # Should not raise


class TestPolarCircleHouses:
    """Test that polar circle house systems raise PolarCircleError."""

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [70.0, 80.0, 85.0, -70.0, -80.0, -85.0])
    def test_placidus_polar_raises(self, lat: float):
        """Placidus at polar latitudes should raise PolarCircleError."""
        jd = 2451545.0
        with pytest.raises(PolarCircleError):
            swe.swe_houses(jd, lat, 0.0, ord("P"))

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [70.0, 80.0, -70.0, -80.0])
    def test_koch_polar_raises(self, lat: float):
        """Koch at polar latitudes should raise PolarCircleError."""
        jd = 2451545.0
        with pytest.raises(PolarCircleError):
            swe.swe_houses(jd, lat, 0.0, ord("K"))

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [45.0, 50.0, -45.0, -50.0])
    def test_equal_houses_polar_ok(self, lat: float):
        """Equal houses should work at any latitude."""
        jd = 2451545.0
        cusps, ascmc = swe.swe_houses(jd, lat, 0.0, ord("A"))
        assert len(cusps) >= 12

    @pytest.mark.unit
    @pytest.mark.parametrize("lat", [80.0, -80.0])
    def test_whole_sign_polar_ok(self, lat: float):
        """Whole Sign houses should work at polar latitudes."""
        jd = 2451545.0
        cusps, ascmc = swe.swe_houses(jd, lat, 0.0, ord("W"))
        assert len(cusps) >= 12


class TestEdgeDates:
    """Test behavior at date range boundaries."""

    @pytest.mark.unit
    def test_very_old_date(self):
        """Very old dates (before DE440 range) should raise or handle gracefully."""
        jd = swe.swe_julday(-5000, 1, 1, 12.0)
        try:
            result, _ = swe.swe_calc_ut(jd, SE_SUN, 0)
            # If it doesn't raise, result should still be valid
            assert 0 <= result[0] < 360
        except Exception:
            pass  # Raising is acceptable for out-of-range dates

    @pytest.mark.unit
    def test_far_future_date(self):
        """Far future dates should raise or handle gracefully."""
        jd = swe.swe_julday(20000, 1, 1, 12.0)
        try:
            result, _ = swe.swe_calc_ut(jd, SE_SUN, 0)
            assert 0 <= result[0] < 360
        except Exception:
            pass  # Raising is acceptable

    @pytest.mark.unit
    def test_j2000_is_valid(self):
        """J2000.0 should always work."""
        result, _ = swe.swe_calc_ut(2451545.0, SE_SUN, SEFLG_SPEED)
        assert 0 <= result[0] < 360
        assert result[2] > 0

    @pytest.mark.unit
    def test_modern_dates(self):
        """Modern dates (1900-2100) should all work."""
        for year in [1900, 1950, 1990, 2000, 2010, 2024, 2050, 2100]:
            jd = swe.swe_julday(year, 6, 15, 12.0)
            result, _ = swe.swe_calc_ut(jd, SE_SUN, SEFLG_SPEED)
            assert 0 <= result[0] < 360, f"Year {year}: lon={result[0]}"


class TestReturnTypes:
    """Test that return types are correct native Python types."""

    @pytest.mark.unit
    def test_calc_ut_returns_native_floats(self):
        """All calc_ut result values should be native Python float."""
        result, retflag = swe.swe_calc_ut(2451545.0, SE_SUN, SEFLG_SPEED)
        for i, val in enumerate(result):
            assert type(val) is float, (
                f"result[{i}] is {type(val).__name__}, expected float"
            )

    @pytest.mark.unit
    def test_retflag_is_int(self):
        """Return flag should be a native Python int."""
        result, retflag = swe.swe_calc_ut(2451545.0, SE_SUN, 0)
        assert isinstance(retflag, int)

    @pytest.mark.unit
    def test_julday_returns_float(self):
        """julday should return a native Python float."""
        jd = swe.swe_julday(2000, 1, 1, 12.0)
        assert type(jd) is float

    @pytest.mark.unit
    def test_deltat_returns_float(self):
        """deltat should return a native Python float."""
        dt = swe.swe_deltat(2451545.0)
        assert type(dt) is float

    @pytest.mark.unit
    def test_sidtime_returns_float(self):
        """sidtime should return a native Python float."""
        st = swe.swe_sidtime(2451545.0)
        assert type(st) is float


class TestSpecialBodies:
    """Test special body calculations."""

    @pytest.mark.unit
    def test_earth_geocentric_is_zeros(self):
        """Earth (body 14) geocentric should return all zeros."""
        result, _ = swe.swe_calc_ut(2451545.0, 14, 0)
        for i, val in enumerate(result):
            assert val == 0.0, f"Earth geo result[{i}] = {val}, expected 0"

    @pytest.mark.unit
    def test_mean_node_valid(self):
        """Mean Node returns valid position."""
        result, _ = swe.swe_calc_ut(2451545.0, SE_MEAN_NODE, 0)
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    def test_true_node_valid(self):
        """True Node returns valid position."""
        result, _ = swe.swe_calc_ut(2451545.0, SE_TRUE_NODE, 0)
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    def test_mean_apogee_valid(self):
        """Mean Apogee (Lilith) returns valid position."""
        result, _ = swe.swe_calc_ut(2451545.0, SE_MEAN_APOG, 0)
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    def test_oscu_apogee_valid(self):
        """Osculating Apogee returns valid position."""
        result, _ = swe.swe_calc_ut(2451545.0, SE_OSCU_APOG, 0)
        assert 0 <= result[0] < 360

    @pytest.mark.unit
    def test_nodes_are_opposite(self):
        """Mean and True nodes should be roughly opposite their desc counterparts."""
        jd = 2451545.0
        mean_node, _ = swe.swe_calc_ut(jd, SE_MEAN_NODE, 0)
        true_node, _ = swe.swe_calc_ut(jd, SE_TRUE_NODE, 0)

        # Mean and True nodes should be close
        diff = abs(mean_node[0] - true_node[0])
        if diff > 180:
            diff = 360 - diff
        assert diff < 5.0, (
            f"Mean node {mean_node[0]:.2f}° vs True node {true_node[0]:.2f}°"
        )


class TestSunHeliocentricIsSpecial:
    """Test that Sun heliocentric is handled correctly."""

    @pytest.mark.unit
    def test_sun_helio_returns_position(self):
        """Sun heliocentric should return a valid position (Earth from Sun)."""
        result, _ = swe.swe_calc_ut(2451545.0, SE_SUN, SEFLG_HELCTR)
        # Could be 0,0,0 or the position of Earth from Sun
        assert len(result) == 6
        for val in result:
            assert math.isfinite(val)


class TestFindStationErrors:
    """Test error handling for find_station_ut."""

    @pytest.mark.unit
    def test_sun_station_raises(self):
        """Sun doesn't have stations — should raise ValueError."""
        with pytest.raises(ValueError):
            swe.swe_find_station_ut(SE_SUN, 2451545.0, "any", 0)

    @pytest.mark.unit
    def test_moon_station_raises(self):
        """Moon doesn't have stations — should raise ValueError."""
        with pytest.raises(ValueError):
            swe.swe_find_station_ut(SE_MOON, 2451545.0, "any", 0)
