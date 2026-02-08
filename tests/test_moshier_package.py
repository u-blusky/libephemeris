"""
Tests for the Moshier semi-analytical ephemeris package.

Tests the basic functionality of the moshier package, including:
- Package imports
- Basic calculations for all supported bodies
- Precession and nutation
- Utility functions
"""

from __future__ import annotations

import math

import pytest

from libephemeris.moshier import (
    # Main functions
    calc_position,
    is_moshier_body,
    # Body IDs
    MOSHIER_SUN,
    MOSHIER_MOON,
    MOSHIER_MERCURY,
    MOSHIER_VENUS,
    MOSHIER_EARTH,
    MOSHIER_MARS,
    MOSHIER_JUPITER,
    MOSHIER_SATURN,
    MOSHIER_URANUS,
    MOSHIER_NEPTUNE,
    MOSHIER_PLUTO,
    MOSHIER_BODIES,
    # Precession
    mean_obliquity,
    true_obliquity,
    nutation_angles,
    # Utilities
    normalize_angle,
    ecliptic_to_equatorial,
    equatorial_to_ecliptic,
    J2000,
)


class TestMoshierPackageImports:
    """Test that the moshier package imports correctly."""

    def test_import_package(self):
        """Test that the package can be imported."""
        import libephemeris.moshier

        assert libephemeris.moshier is not None

    def test_import_submodules(self):
        """Test that all submodules can be imported."""
        from libephemeris.moshier import utils
        from libephemeris.moshier import precession
        from libephemeris.moshier import vsop87
        from libephemeris.moshier import elp82b
        from libephemeris.moshier import pluto

        assert utils is not None
        assert precession is not None
        assert vsop87 is not None
        assert elp82b is not None
        assert pluto is not None

    def test_body_constants(self):
        """Test that body ID constants are correct."""
        assert MOSHIER_SUN == 0
        assert MOSHIER_MOON == 1
        assert MOSHIER_MERCURY == 2
        assert MOSHIER_VENUS == 3
        assert MOSHIER_MARS == 4
        assert MOSHIER_JUPITER == 5
        assert MOSHIER_SATURN == 6
        assert MOSHIER_URANUS == 7
        assert MOSHIER_NEPTUNE == 8
        assert MOSHIER_PLUTO == 9
        assert MOSHIER_EARTH == 14


class TestIsMoshierBody:
    """Test is_moshier_body function."""

    def test_supported_bodies(self):
        """Test that all expected bodies are supported."""
        for body_id in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14]:
            assert is_moshier_body(body_id), f"Body {body_id} should be supported"

    def test_unsupported_bodies(self):
        """Test that unsupported bodies return False."""
        for body_id in [10, 11, 12, 13, 15, 16, 17, 100, -1]:
            assert not is_moshier_body(body_id), (
                f"Body {body_id} should not be supported"
            )


class TestCalcPosition:
    """Test calc_position for all bodies at J2000.0."""

    @pytest.fixture
    def j2000_jd(self):
        """Julian Day of J2000.0 epoch."""
        return J2000  # 2451545.0

    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (MOSHIER_SUN, "Sun"),
            (MOSHIER_MOON, "Moon"),
            (MOSHIER_MERCURY, "Mercury"),
            (MOSHIER_VENUS, "Venus"),
            (MOSHIER_MARS, "Mars"),
            (MOSHIER_JUPITER, "Jupiter"),
            (MOSHIER_SATURN, "Saturn"),
            (MOSHIER_URANUS, "Uranus"),
            (MOSHIER_NEPTUNE, "Neptune"),
            (MOSHIER_PLUTO, "Pluto"),
        ],
    )
    def test_calc_position_returns_tuple(self, j2000_jd, body_id, body_name):
        """Test that calc_position returns a 6-tuple for all bodies."""
        result = calc_position(j2000_jd, body_id)

        assert isinstance(result, tuple), f"{body_name}: Expected tuple"
        assert len(result) == 6, f"{body_name}: Expected 6 elements"

        lon, lat, dist, dlon, dlat, ddist = result

        # All values should be finite
        assert math.isfinite(lon), f"{body_name}: longitude not finite"
        assert math.isfinite(lat), f"{body_name}: latitude not finite"
        assert math.isfinite(dist), f"{body_name}: distance not finite"
        assert math.isfinite(dlon), f"{body_name}: dlon not finite"
        assert math.isfinite(dlat), f"{body_name}: dlat not finite"
        assert math.isfinite(ddist), f"{body_name}: ddist not finite"

    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (MOSHIER_SUN, "Sun"),
            (MOSHIER_MOON, "Moon"),
            (MOSHIER_MERCURY, "Mercury"),
            (MOSHIER_VENUS, "Venus"),
            (MOSHIER_MARS, "Mars"),
            (MOSHIER_JUPITER, "Jupiter"),
            (MOSHIER_SATURN, "Saturn"),
            (MOSHIER_URANUS, "Uranus"),
            (MOSHIER_NEPTUNE, "Neptune"),
            (MOSHIER_PLUTO, "Pluto"),
        ],
    )
    def test_calc_position_reasonable_values(self, j2000_jd, body_id, body_name):
        """Test that calc_position returns physically reasonable values."""
        lon, lat, dist, dlon, dlat, ddist = calc_position(j2000_jd, body_id)

        # Longitude should be in [0, 360)
        assert 0.0 <= lon < 360.0, f"{body_name}: longitude {lon} out of range"

        # Latitude should be in [-90, 90]
        assert -90.0 <= lat <= 90.0, f"{body_name}: latitude {lat} out of range"

        # Distance should be positive
        assert dist > 0, f"{body_name}: distance {dist} should be positive"

        # For planets, distance should be less than 100 AU
        if body_id != MOSHIER_MOON:
            assert dist < 100.0, f"{body_name}: distance {dist} AU too large"

    def test_unsupported_body_raises(self, j2000_jd):
        """Test that unsupported body ID raises ValueError."""
        with pytest.raises(ValueError):
            calc_position(j2000_jd, 100)

    def test_sun_at_j2000(self, j2000_jd):
        """Test Sun position at J2000.0 epoch."""
        lon, lat, dist, _, _, _ = calc_position(j2000_jd, MOSHIER_SUN)

        # At J2000.0, Sun is around 280° ecliptic longitude
        assert 270.0 < lon < 290.0, f"Sun longitude {lon} unexpected at J2000"

        # Sun latitude should be near 0
        assert abs(lat) < 1.0, f"Sun latitude {lat} should be near zero"

        # Sun distance should be about 1 AU
        assert 0.98 < dist < 1.02, f"Sun distance {dist} should be ~1 AU"

    def test_moon_at_j2000(self, j2000_jd):
        """Test Moon position at J2000.0 epoch."""
        lon, lat, dist, dlon, _, _ = calc_position(j2000_jd, MOSHIER_MOON)

        # Moon distance should be about 0.0026 AU (385000 km / 149597870.7 km)
        assert 0.002 < dist < 0.003, f"Moon distance {dist} should be ~0.0026 AU"

        # Moon moves about 13°/day
        assert 10.0 < dlon < 16.0, f"Moon velocity {dlon}°/day unexpected"


class TestMeanObliquity:
    """Test mean obliquity calculation."""

    def test_obliquity_at_j2000(self):
        """Test mean obliquity at J2000.0."""
        eps = mean_obliquity(J2000)

        # IAU 2006 value: 84381.406 arcsec = 23.4392911° (approx)
        assert 23.4 < eps < 23.5, f"Mean obliquity {eps}° unexpected at J2000"

    def test_obliquity_decreases(self):
        """Test that obliquity is decreasing over time (current epoch)."""
        eps_past = mean_obliquity(J2000 - 36525)  # 100 years ago
        eps_now = mean_obliquity(J2000)
        eps_future = mean_obliquity(J2000 + 36525)  # 100 years in future

        # Obliquity is currently decreasing
        assert eps_past > eps_now > eps_future


class TestNutationAngles:
    """Test nutation angle calculations."""

    def test_nutation_at_j2000(self):
        """Test nutation angles at J2000.0."""
        dpsi, deps = nutation_angles(J2000)

        # Nutation in longitude typically < 20 arcsec = 0.0056°
        assert abs(dpsi) < 0.01, f"Nutation dpsi {dpsi}° too large"

        # Nutation in obliquity typically < 10 arcsec = 0.0028°
        assert abs(deps) < 0.005, f"Nutation deps {deps}° too large"

    def test_true_obliquity(self):
        """Test true obliquity includes nutation."""
        eps_mean = mean_obliquity(J2000)
        eps_true = true_obliquity(J2000)
        _, deps = nutation_angles(J2000)

        # True obliquity = mean + nutation in obliquity
        assert abs(eps_true - (eps_mean + deps)) < 1e-10


class TestNormalizeAngle:
    """Test angle normalization function."""

    def test_normalize_positive(self):
        """Test normalization of positive angles."""
        assert normalize_angle(0.0) == 0.0
        assert normalize_angle(180.0) == 180.0
        assert normalize_angle(359.999) == pytest.approx(359.999)
        assert normalize_angle(360.0) == pytest.approx(0.0)
        assert normalize_angle(720.0) == pytest.approx(0.0)
        assert normalize_angle(370.0) == pytest.approx(10.0)

    def test_normalize_negative(self):
        """Test normalization of negative angles."""
        assert normalize_angle(-10.0) == pytest.approx(350.0)
        assert normalize_angle(-360.0) == pytest.approx(0.0)
        assert normalize_angle(-370.0) == pytest.approx(350.0)


class TestCoordinateConversions:
    """Test coordinate conversion functions."""

    def test_ecliptic_equatorial_roundtrip(self):
        """Test ecliptic <-> equatorial roundtrip."""
        lon = 45.0
        lat = 10.0
        obliquity = 23.44

        ra, dec = ecliptic_to_equatorial(lon, lat, obliquity)
        lon2, lat2 = equatorial_to_ecliptic(ra, dec, obliquity)

        assert lon2 == pytest.approx(lon, abs=1e-10)
        assert lat2 == pytest.approx(lat, abs=1e-10)

    def test_ecliptic_equatorial_pole(self):
        """Test ecliptic to equatorial at north ecliptic pole."""
        # North ecliptic pole
        lon = 0.0
        lat = 90.0
        obliquity = 23.44

        ra, dec = ecliptic_to_equatorial(lon, lat, obliquity)

        # At north ecliptic pole, dec = 90 - obliquity
        expected_dec = 90.0 - obliquity
        assert dec == pytest.approx(expected_dec, abs=0.1)


class TestVSOP87:
    """Test VSOP87 planetary calculations."""

    def test_earth_heliocentric(self):
        """Test Earth heliocentric coordinates."""
        from libephemeris.moshier.vsop87 import calc_earth_heliocentric

        lon, lat, r = calc_earth_heliocentric(J2000)

        # Earth longitude at J2000 is around 100°
        assert 90.0 < lon < 110.0, f"Earth longitude {lon}° unexpected"

        # Earth latitude should be very close to 0
        assert abs(lat) < 0.001, f"Earth latitude {lat}° should be ~0"

        # Earth radius should be about 1 AU
        assert 0.98 < r < 1.02, f"Earth distance {r} AU unexpected"

    def test_sun_geocentric(self):
        """Test Sun geocentric coordinates."""
        from libephemeris.moshier.vsop87 import calc_sun_geocentric

        lon, lat, r = calc_sun_geocentric(J2000)

        # Sun is opposite Earth: if Earth is at ~100°, Sun is at ~280°
        # At J2000.0, Sun is around 280° ecliptic longitude
        assert 270.0 < lon < 290.0, f"Sun longitude {lon}° unexpected"


class TestELP82B:
    """Test ELP 2000-82B lunar calculations."""

    def test_moon_position_basic(self):
        """Test basic Moon position calculation."""
        from libephemeris.moshier.elp82b import calc_moon_position

        lon, lat, dist_km = calc_moon_position(J2000)

        # Moon longitude should be in valid range
        assert 0.0 <= lon < 360.0

        # Moon latitude should be within inclination
        assert abs(lat) < 6.0, f"Moon latitude {lat}° too high"

        # Moon distance should be around 385000 km
        assert 356000 < dist_km < 407000, f"Moon distance {dist_km} km unexpected"

    def test_moon_velocity(self):
        """Test Moon velocity is reasonable."""
        _, _, _, dlon, dlat, _ = calc_position(J2000, MOSHIER_MOON)

        # Moon moves about 13.2°/day in longitude
        assert 12.0 < dlon < 15.0, f"Moon dlon {dlon}°/day unexpected"

        # Moon latitude change is small
        assert abs(dlat) < 1.0, f"Moon dlat {dlat}°/day unexpected"


class TestPluto:
    """Test Pluto calculations."""

    def test_pluto_at_j2000(self):
        """Test Pluto position at J2000.0."""
        lon, lat, dist, _, _, _ = calc_position(J2000, MOSHIER_PLUTO)

        # Pluto was around 250° at J2000
        assert 240.0 < lon < 260.0, f"Pluto longitude {lon}° unexpected"

        # Pluto has high inclination
        assert abs(lat) < 20.0, f"Pluto latitude {lat}° too high"

        # Pluto distance 29-49 AU
        assert 28.0 < dist < 50.0, f"Pluto distance {dist} AU unexpected"

    def test_pluto_keplerian_fallback(self):
        """Test Pluto Keplerian fallback for dates outside periodic term range."""
        from libephemeris.moshier.pluto import calc_pluto_heliocentric

        # Year 1800 (before 1885)
        jd_1800 = 2378496.5

        lon, lat, r = calc_pluto_heliocentric(jd_1800)

        # Should still return valid values
        assert 0.0 <= lon < 360.0
        assert -90.0 <= lat <= 90.0
        assert r > 20.0


class TestExtendedDateRange:
    """Test calculations at extreme dates."""

    @pytest.mark.slow
    def test_calculation_1000ce(self):
        """Test calculation at 1000 CE (outside DE440 range)."""
        jd_1000ce = 2086302.5  # Jan 1, 1000 CE

        # Should work without error
        lon, lat, dist, _, _, _ = calc_position(jd_1000ce, MOSHIER_SUN)

        assert 0.0 <= lon < 360.0
        assert 0.98 < dist < 1.02

    @pytest.mark.slow
    def test_calculation_500bce(self):
        """Test calculation at 500 BCE."""
        jd_500bce = 1538558.5  # Approximately 500 BCE

        # Should work without error
        lon, lat, dist, _, _, _ = calc_position(jd_500bce, MOSHIER_MARS)

        assert 0.0 <= lon < 360.0
        assert dist > 0
