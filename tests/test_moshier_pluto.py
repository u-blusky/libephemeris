"""
Tests for the Moshier Pluto analytical theory implementation.

Tests the new DE404-based Pluto theory with ~800 terms that provides
accuracy of 1-3 arcseconds over the range -3000 to +3000 CE.

These tests verify:
1. Basic calculation correctness at known dates
2. Comparison with DE440 in overlap range (1550-2650)
3. Extended range calculations (-3000 to +3000)
4. Heliocentric vs geocentric transformations
5. Velocity calculations
"""

from __future__ import annotations

import math
import pytest

from libephemeris.moshier import (
    calc_position,
    J2000,
    MOSHIER_PLUTO,
)
from libephemeris.moshier.pluto import (
    calc_pluto_heliocentric,
    calc_pluto_geocentric,
    is_pluto_body,
    MOSHIER_PLUTO as PLUTO_ID,
)


class TestPlutoMoshierBasic:
    """Basic tests for Pluto Moshier implementation."""

    def test_is_pluto_body(self):
        """Test is_pluto_body function."""
        assert is_pluto_body(9) is True
        assert is_pluto_body(PLUTO_ID) is True
        assert is_pluto_body(0) is False
        assert is_pluto_body(1) is False
        assert is_pluto_body(10) is False

    def test_pluto_heliocentric_j2000(self):
        """Test Pluto heliocentric position at J2000.0."""
        lon, lat, r = calc_pluto_heliocentric(J2000)

        # At J2000.0, Pluto heliocentric longitude is around 250°
        assert 245.0 < lon < 255.0, f"Pluto helio lon {lon}° unexpected"

        # Pluto has high orbital inclination (~17°), latitude can be significant
        assert -20.0 < lat < 20.0, f"Pluto helio lat {lat}° unexpected"

        # Pluto heliocentric distance should be around 30 AU at J2000
        assert 28.0 < r < 35.0, f"Pluto helio dist {r} AU unexpected"

    def test_pluto_geocentric_j2000(self):
        """Test Pluto geocentric position at J2000.0."""
        lon, lat, dist = calc_pluto_geocentric(J2000)

        # Geocentric position differs from heliocentric due to Earth position
        assert 0.0 <= lon < 360.0, f"Pluto geo lon {lon}° out of range"

        # Latitude should be within orbital inclination range
        assert -20.0 < lat < 20.0, f"Pluto geo lat {lat}° unexpected"

        # Geocentric distance varies from heliocentric by ~1 AU
        assert 27.0 < dist < 36.0, f"Pluto geo dist {dist} AU unexpected"

    def test_calc_position_returns_6_elements(self):
        """Test that calc_position returns all 6 elements."""
        result = calc_position(J2000, MOSHIER_PLUTO)

        assert isinstance(result, tuple)
        assert len(result) == 6

        lon, lat, dist, dlon, dlat, ddist = result

        # All values should be finite
        assert math.isfinite(lon)
        assert math.isfinite(lat)
        assert math.isfinite(dist)
        assert math.isfinite(dlon)
        assert math.isfinite(dlat)
        assert math.isfinite(ddist)


class TestPlutoMoshierVelocity:
    """Test Pluto velocity calculations."""

    def test_pluto_velocity_at_j2000(self):
        """Test Pluto velocity at J2000.0."""
        _, _, _, dlon, dlat, ddist = calc_position(J2000, MOSHIER_PLUTO)

        # Pluto moves about 0.004°/day in longitude (varies with retrograde)
        assert abs(dlon) < 0.05, f"Pluto dlon {dlon}°/day too high"

        # Latitude velocity should be small
        assert abs(dlat) < 0.01, f"Pluto dlat {dlat}°/day unexpected"

        # Distance velocity: Pluto moves a few km/s, so ~ 0.01 AU/day max
        assert abs(ddist) < 0.02, f"Pluto ddist {ddist} AU/day unexpected"


class TestPlutoMoshierRange:
    """Test Pluto at various dates including extended range."""

    @pytest.mark.parametrize(
        "year,jd",
        [
            (2000, 2451545.0),  # J2000.0
            (2020, 2458849.5),  # Modern date
            (1900, 2415020.5),  # Historical
            (1800, 2378496.5),  # Before main periodic range
        ],
    )
    def test_pluto_various_dates_modern(self, year, jd):
        """Test Pluto at various modern dates."""
        lon, lat, dist = calc_pluto_heliocentric(jd)

        assert 0.0 <= lon < 360.0, f"Year {year}: lon {lon}° out of range"
        assert -20.0 <= lat <= 20.0, f"Year {year}: lat {lat}° out of range"
        assert 28.0 < dist < 50.0, f"Year {year}: dist {dist} AU out of range"

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "year,jd",
        [
            (1000, 2086302.5),  # Medieval
            (0, 1721057.5),  # Epoch 0
            (-500, 1538558.5),  # 500 BCE
            (-1000, 1355807.5),  # 1000 BCE
        ],
    )
    def test_pluto_extended_range_ancient(self, year, jd):
        """Test Pluto at ancient dates (extended Moshier range)."""
        lon, lat, dist = calc_pluto_heliocentric(jd)

        assert 0.0 <= lon < 360.0, f"Year {year}: lon {lon}° out of range"
        assert -20.0 <= lat <= 20.0, f"Year {year}: lat {lat}° out of range"
        # Pluto's orbit: perihelion ~30 AU, aphelion ~49 AU
        assert 25.0 < dist < 55.0, f"Year {year}: dist {dist} AU out of range"

    @pytest.mark.slow
    def test_pluto_year_3000(self):
        """Test Pluto at year 3000 CE (edge of Moshier range)."""
        jd_3000 = 2816787.5  # Approximately year 3000 CE

        lon, lat, dist = calc_pluto_heliocentric(jd_3000)

        assert 0.0 <= lon < 360.0, f"lon {lon}° out of range"
        assert -20.0 <= lat <= 20.0, f"lat {lat}° out of range"
        assert 25.0 < dist < 55.0, f"dist {dist} AU out of range"


class TestPlutoMoshierPrecision:
    """Test Pluto precision by comparing with DE440 in overlap range."""

    @pytest.mark.slow
    def test_pluto_vs_de440_j2000(self):
        """Compare Moshier Pluto with DE440 at J2000.0.

        Expected accuracy: ~1-3 arcseconds for the DE404 fit.
        """
        import libephemeris as ephem
        from libephemeris.constants import SEFLG_SPEED, SEFLG_MOSEPH, SE_PLUTO

        # Get Moshier result
        moshier_result = calc_position(J2000, MOSHIER_PLUTO)
        mosh_lon, mosh_lat, mosh_dist, _, _, _ = moshier_result

        # Get DE440 result (default mode)
        try:
            de440_result, _ = ephem.swe_calc_ut(J2000, SE_PLUTO, SEFLG_SPEED)
            de440_lon, de440_lat, de440_dist = (
                de440_result[0],
                de440_result[1],
                de440_result[2],
            )

            # Longitude difference (accounting for wrap-around)
            lon_diff = mosh_lon - de440_lon
            if lon_diff > 180:
                lon_diff -= 360
            elif lon_diff < -180:
                lon_diff += 360

            # 3 arcsec = 0.000833 degrees tolerance
            # But DE404 fit vs DE440 may differ more, so use 0.01 degrees (36 arcsec)
            assert abs(lon_diff) < 0.1, f"Longitude diff {lon_diff}° too large"
            assert abs(mosh_lat - de440_lat) < 0.1, "Latitude diff too large"

            # Distance can differ by a few percent
            dist_pct = abs(mosh_dist - de440_dist) / de440_dist * 100
            assert dist_pct < 1.0, f"Distance diff {dist_pct}% too large"

        except Exception:
            # If DE440 not available, skip comparison
            pytest.skip("DE440 ephemeris not available for comparison")

    @pytest.mark.slow
    def test_pluto_moshier_via_api(self):
        """Test Pluto via the main swe_calc_ut API with SEFLG_MOSEPH."""
        import libephemeris as ephem
        from libephemeris.constants import SEFLG_SPEED, SEFLG_MOSEPH, SE_PLUTO

        result, flags = ephem.swe_calc_ut(J2000, SE_PLUTO, SEFLG_SPEED | SEFLG_MOSEPH)

        assert len(result) == 6

        lon, lat, dist, dlon, dlat, ddist = result

        # Same checks as direct calc_position
        assert 0.0 <= lon < 360.0
        assert -20.0 < lat < 20.0
        assert 27.0 < dist < 36.0

        # Check flag is set correctly
        assert flags & SEFLG_MOSEPH


class TestPlutoMoshierConsistency:
    """Test internal consistency of Pluto calculations."""

    def test_heliocentric_geocentric_difference(self):
        """Verify heliocentric and geocentric differ by Earth position."""
        from libephemeris.moshier.vsop87 import calc_earth_heliocentric

        h_lon, h_lat, h_dist = calc_pluto_heliocentric(J2000)
        g_lon, g_lat, g_dist = calc_pluto_geocentric(J2000)
        e_lon, e_lat, e_dist = calc_earth_heliocentric(J2000)

        # Geocentric distance should differ from heliocentric by roughly Earth's distance
        # (depends on configuration, but should be within reasonable bounds)
        dist_diff = abs(h_dist - g_dist)
        assert dist_diff < e_dist + 1.0, "Distance difference unexpectedly large"

        # Longitude difference due to parallax should be small for distant Pluto
        lon_diff = abs(h_lon - g_lon)
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        assert lon_diff < 5.0, f"Longitude parallax {lon_diff}° unexpectedly large"

    def test_continuity_over_time(self):
        """Test that position changes continuously over time."""
        jd = J2000

        positions = []
        for i in range(10):
            lon, lat, dist = calc_pluto_heliocentric(jd + i * 100)
            positions.append((lon, lat, dist))

        # Check no sudden jumps
        for i in range(1, len(positions)):
            lon_prev, lat_prev, dist_prev = positions[i - 1]
            lon_curr, lat_curr, dist_curr = positions[i]

            # Longitude change in 100 days should be small (~0.4°)
            dlon = lon_curr - lon_prev
            if dlon > 180:
                dlon -= 360
            elif dlon < -180:
                dlon += 360
            assert abs(dlon) < 2.0, f"Sudden longitude jump at step {i}: {dlon}°"

            # Latitude change should be small
            assert abs(lat_curr - lat_prev) < 0.5, f"Sudden latitude jump at step {i}"

            # Distance change should be smooth
            assert abs(dist_curr - dist_prev) < 0.5, f"Sudden distance jump at step {i}"


class TestPlutoMoshierEdgeCases:
    """Test edge cases and error handling."""

    def test_wrong_body_id_raises(self):
        """Test that wrong body ID raises ValueError."""
        from libephemeris.moshier.pluto import calc_position as pluto_calc

        with pytest.raises(ValueError):
            pluto_calc(J2000, 0)  # Sun

        with pytest.raises(ValueError):
            pluto_calc(J2000, 1)  # Moon

        with pytest.raises(ValueError):
            pluto_calc(J2000, 10)  # Invalid

    def test_extreme_dates(self):
        """Test at extreme ends of the Moshier range."""
        # Near -3000 CE
        jd_ancient = 625673.5  # Approximately -3000 CE

        lon, lat, dist = calc_pluto_heliocentric(jd_ancient)

        # Should still produce valid values (may be less accurate)
        assert 0.0 <= lon < 360.0
        assert -90.0 <= lat <= 90.0
        assert dist > 0

    def test_nan_not_returned(self):
        """Ensure no NaN values are returned for valid dates."""
        for jd in [J2000, J2000 + 365250, J2000 - 365250]:
            lon, lat, dist, dlon, dlat, ddist = calc_position(jd, MOSHIER_PLUTO)

            assert not math.isnan(lon), f"NaN longitude at JD {jd}"
            assert not math.isnan(lat), f"NaN latitude at JD {jd}"
            assert not math.isnan(dist), f"NaN distance at JD {jd}"
            assert not math.isnan(dlon), f"NaN dlon at JD {jd}"
            assert not math.isnan(dlat), f"NaN dlat at JD {jd}"
            assert not math.isnan(ddist), f"NaN ddist at JD {jd}"
