"""
Comprehensive tests for coordinate transform functions.

Verifies swe_cotrans (ecliptic <-> equatorial) and related
coordinate transformations work correctly with round-trips,
edge cases, and known values.

Sign convention for swe_cotrans:
  - NEGATIVE obliquity: ecliptic -> equatorial
  - POSITIVE obliquity: equatorial -> ecliptic
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
    SEFLG_SPEED,
    SEFLG_EQUATORIAL,
)


# Mean obliquity at J2000.0 ~23.4393°
J2000_OBLIQUITY = 23.4392911


class TestCotransBasic:
    """Basic swe_cotrans functionality."""

    @pytest.mark.unit
    def test_cotrans_returns_3_tuple(self):
        """swe_cotrans returns a 3-element tuple."""
        result = swe.swe_cotrans((100.0, 20.0, 1.0), -J2000_OBLIQUITY)
        assert len(result) == 3

    @pytest.mark.unit
    def test_cotrans_returns_native_floats(self):
        """All returned values should be native Python floats."""
        result = swe.swe_cotrans((100.0, 20.0, 1.0), -J2000_OBLIQUITY)
        for i, val in enumerate(result):
            assert type(val) is float, (
                f"result[{i}] is {type(val).__name__}, expected float"
            )

    @pytest.mark.unit
    def test_cotrans_all_finite(self):
        """All returned values should be finite."""
        result = swe.swe_cotrans((100.0, 20.0, 1.0), -J2000_OBLIQUITY)
        for i, val in enumerate(result):
            assert math.isfinite(val), f"result[{i}] = {val}"


class TestCotransEclipticToEquatorial:
    """Test ecliptic to equatorial conversion (negative obliquity)."""

    @pytest.mark.unit
    def test_vernal_equinox_point(self):
        """Vernal equinox (0°,0°) maps to RA=0°, Dec=0° in equatorial."""
        result = swe.swe_cotrans((0.0, 0.0, 1.0), -J2000_OBLIQUITY)
        ra, dec, dist = result
        assert abs(ra) < 0.01 or abs(ra - 360) < 0.01, f"RA={ra}"
        assert abs(dec) < 0.01, f"Dec={dec}"

    @pytest.mark.unit
    def test_summer_solstice_point(self):
        """Summer solstice (90°,0°) should have Dec ~+23.44°."""
        result = swe.swe_cotrans((90.0, 0.0, 1.0), -J2000_OBLIQUITY)
        ra, dec, dist = result
        assert abs(ra - 90) < 0.5, f"RA={ra}"
        assert abs(dec - J2000_OBLIQUITY) < 0.5, f"Dec={dec}"

    @pytest.mark.unit
    def test_autumnal_equinox(self):
        """Autumnal equinox (180°,0°) maps to RA=180°, Dec=0°."""
        result = swe.swe_cotrans((180.0, 0.0, 1.0), -J2000_OBLIQUITY)
        ra, dec, dist = result
        assert abs(ra - 180) < 0.5, f"RA={ra}"
        assert abs(dec) < 0.5, f"Dec={dec}"

    @pytest.mark.unit
    def test_winter_solstice(self):
        """Winter solstice (270°,0°) should have Dec ~-23.44°."""
        result = swe.swe_cotrans((270.0, 0.0, 1.0), -J2000_OBLIQUITY)
        ra, dec, dist = result
        assert abs(ra - 270) < 0.5, f"RA={ra}"
        assert abs(dec + J2000_OBLIQUITY) < 0.5, f"Dec={dec}"

    @pytest.mark.unit
    def test_ecliptic_north_pole(self):
        """Ecliptic north pole (lat=90°) maps to Dec ~ 90 - obliquity."""
        result = swe.swe_cotrans((0.0, 90.0, 1.0), -J2000_OBLIQUITY)
        ra, dec, dist = result
        expected_dec = 90.0 - J2000_OBLIQUITY
        assert abs(dec - expected_dec) < 1.0, (
            f"Ecliptic NP: dec={dec}, expected ~{expected_dec}"
        )


class TestCotransRoundTrip:
    """Test ecliptic -> equatorial -> ecliptic round-trip.

    ecl->equ uses negative obliquity, equ->ecl uses positive.
    """

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "lon,lat",
        [
            (0, 0),
            (45, 10),
            (90, 0),
            (135, -15),
            (180, 0),
            (225, 20),
            (270, 0),
            (315, -5),
            (359.99, 0),
            (100, 45),
            (200, -45),
        ],
    )
    def test_round_trip(self, lon: float, lat: float):
        """Converting ecl->equ->ecl should recover original coordinates."""
        dist = 1.5
        # Ecliptic to equatorial: negative obliquity
        equ = swe.swe_cotrans((lon, lat, dist), -J2000_OBLIQUITY)
        # Equatorial to ecliptic: positive obliquity
        ecl = swe.swe_cotrans(equ, J2000_OBLIQUITY)

        lon_diff = abs(ecl[0] - lon)
        if lon_diff > 180:
            lon_diff = 360 - lon_diff
        assert lon_diff < 1e-6, f"Lon: {ecl[0]} != {lon} (diff {lon_diff})"
        assert abs(ecl[1] - lat) < 1e-6, f"Lat: {ecl[1]} != {lat}"
        assert abs(ecl[2] - dist) < 1e-6, f"Dist: {ecl[2]} != {dist}"

    @pytest.mark.unit
    def test_distance_preserved(self):
        """Distance should be preserved through transform."""
        dist = 5.2
        result = swe.swe_cotrans((120.0, 15.0, dist), -J2000_OBLIQUITY)
        assert abs(result[2] - dist) < 1e-6, f"Distance changed: {result[2]} != {dist}"


class TestCotransConsistencyWithCalcUt:
    """Verify cotrans matches calc_ut equatorial output."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (SE_SUN, "Sun"),
            (SE_MOON, "Moon"),
            (SE_MARS, "Mars"),
            (SE_JUPITER, "Jupiter"),
        ],
    )
    def test_cotrans_matches_equatorial_flag(self, body_id: int, name: str):
        """cotrans(ecliptic, -eps) should approximate SEFLG_EQUATORIAL output."""
        jd = 2451545.0
        ecl, _ = swe.swe_calc_ut(jd, body_id, 0)
        equ, _ = swe.swe_calc_ut(jd, body_id, SEFLG_EQUATORIAL)

        # Transform ecliptic to equatorial (negative obliquity)
        transformed = swe.swe_cotrans((ecl[0], ecl[1], ecl[2]), -J2000_OBLIQUITY)

        # Should be close (not exact due to obliquity approximation,
        # nutation, and aberration differences, but within ~1°)
        ra_diff = abs(transformed[0] - equ[0])
        if ra_diff > 180:
            ra_diff = 360 - ra_diff
        dec_diff = abs(transformed[1] - equ[1])

        assert ra_diff < 1.5, (
            f"{name}: RA diff {ra_diff:.4f}° between cotrans and EQUATORIAL"
        )
        assert dec_diff < 1.5, (
            f"{name}: Dec diff {dec_diff:.4f}° between cotrans and EQUATORIAL"
        )


class TestCotransEdgeCases:
    """Edge cases for coordinate transforms."""

    @pytest.mark.unit
    def test_zero_obliquity(self):
        """With zero obliquity, ecliptic = equatorial."""
        result = swe.swe_cotrans((100.0, 20.0, 1.0), 0.0)
        assert abs(result[0] - 100.0) < 1e-10
        assert abs(result[1] - 20.0) < 1e-10

    @pytest.mark.unit
    def test_90_degree_obliquity(self):
        """With 90° obliquity (hypothetical extreme)."""
        result = swe.swe_cotrans((0.0, 0.0, 1.0), 90.0)
        # Should still return finite values
        for val in result:
            assert math.isfinite(val)

    @pytest.mark.unit
    @pytest.mark.parametrize("lon", [0, 90, 180, 270, 359.999])
    def test_ecliptic_plane_dec_bounded(self, lon: float):
        """On ecliptic (lat=0), declination bounded by obliquity."""
        result = swe.swe_cotrans((lon, 0.0, 1.0), -J2000_OBLIQUITY)
        dec = result[1]
        assert abs(dec) <= J2000_OBLIQUITY + 0.1, (
            f"Dec={dec} exceeds obliquity for lon={lon}"
        )

    @pytest.mark.unit
    def test_negative_longitude(self):
        """Negative longitude should be handled (wrapped to 0-360)."""
        r1 = swe.swe_cotrans((-10.0, 0.0, 1.0), -J2000_OBLIQUITY)
        r2 = swe.swe_cotrans((350.0, 0.0, 1.0), -J2000_OBLIQUITY)
        for i in range(3):
            diff = abs(r1[i] - r2[i])
            if i == 0 and diff > 180:
                diff = 360 - diff
            assert diff < 0.01, f"Element {i}: {r1[i]} vs {r2[i]}"

    @pytest.mark.unit
    def test_positive_negative_obliquity_inverse(self):
        """Positive and negative obliquity are inverse operations."""
        coords = (120.0, 30.0, 2.0)
        # Apply -eps then +eps => should get back original
        transformed = swe.swe_cotrans(coords, -J2000_OBLIQUITY)
        recovered = swe.swe_cotrans(transformed, J2000_OBLIQUITY)
        for i in range(3):
            diff = abs(recovered[i] - coords[i])
            if i == 0 and diff > 180:
                diff = 360 - diff
            assert diff < 1e-8, f"Element {i}: {recovered[i]} != {coords[i]}"
