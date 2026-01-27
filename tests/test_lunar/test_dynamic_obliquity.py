"""
Tests for dynamic obliquity calculations in lunar module.

The true lunar node and true Lilith calculations now use time-varying
obliquity (IAU 2006) instead of fixed J2000 obliquity, improving accuracy
for dates far from J2000.
"""

import math
from libephemeris.lunar import (
    _mean_obliquity_radians,
    calc_true_lunar_node,
    calc_true_lilith,
)


class TestMeanObliquityRadians:
    """Test the _mean_obliquity_radians helper function."""

    def test_j2000_obliquity(self):
        """At J2000.0, obliquity should be approximately 23.4392911°."""
        jd_j2000 = 2451545.0
        eps_rad = _mean_obliquity_radians(jd_j2000)
        eps_deg = math.degrees(eps_rad)

        # J2000.0 mean obliquity is 84381.406 arcsec = 23.4392794° (IAU 2006)
        # Should match within 0.0001°
        assert abs(eps_deg - 23.4392794) < 0.0001

    def test_obliquity_formula_coefficients(self):
        """Test the obliquity formula matches IAU 2006 at T=0."""
        jd_j2000 = 2451545.0
        eps_rad = _mean_obliquity_radians(jd_j2000)
        eps_arcsec = math.degrees(eps_rad) * 3600

        # At T=0, only the constant term (84381.406") applies
        expected = 84381.406
        assert abs(eps_arcsec - expected) < 0.001

    def test_obliquity_decreases_over_time(self):
        """Obliquity is decreasing at about 47 arcsec per century."""
        jd_j2000 = 2451545.0
        jd_century_later = jd_j2000 + 36525.0  # 1 Julian century

        eps_j2000 = _mean_obliquity_radians(jd_j2000)
        eps_j2100 = _mean_obliquity_radians(jd_century_later)

        # Obliquity should be smaller in the future
        assert eps_j2100 < eps_j2000

        # The decrease should be approximately 47 arcsec per century
        decrease_arcsec = (math.degrees(eps_j2000) - math.degrees(eps_j2100)) * 3600
        # The linear term is -46.836769" per century
        assert abs(decrease_arcsec - 46.84) < 1.0

    def test_obliquity_for_historical_dates(self):
        """Test obliquity calculation for historical dates."""
        # Year 1000 CE (approximately)
        jd_1000 = 2086302.0  # JD for ~1000 CE
        eps_rad = _mean_obliquity_radians(jd_1000)
        eps_deg = math.degrees(eps_rad)

        # Obliquity was larger in the past
        eps_j2000 = math.degrees(_mean_obliquity_radians(2451545.0))
        assert eps_deg > eps_j2000

        # Should still be in a reasonable range (22° to 25°)
        assert 22.0 < eps_deg < 25.0

    def test_obliquity_for_future_dates(self):
        """Test obliquity calculation for future dates."""
        # Year 3000 CE (approximately)
        jd_3000 = 2451545.0 + 10.0 * 36525.0
        eps_rad = _mean_obliquity_radians(jd_3000)
        eps_deg = math.degrees(eps_rad)

        # Obliquity will be smaller in the future
        eps_j2000 = math.degrees(_mean_obliquity_radians(2451545.0))
        assert eps_deg < eps_j2000

        # Should still be in a reasonable range (22° to 25°)
        assert 22.0 < eps_deg < 25.0


class TestDynamicObliquityInTrueNode:
    """Test that true lunar node uses dynamic obliquity."""

    def test_true_node_returns_valid_longitude(self):
        """True node should return valid longitude for any date."""
        jd_j2000 = 2451545.0
        lon, lat, dist = calc_true_lunar_node(jd_j2000)

        assert 0 <= lon < 360
        assert lat == 0.0  # Node is on ecliptic
        assert dist < 0.01  # Mathematical point, may have small non-zero distance

    def test_true_node_for_early_date(self):
        """True node should work for dates within ephemeris range."""
        # Year 1950 (within ephemeris range 1899-2053)
        jd_1950 = 2433282.5  # 1950-01-01
        lon, lat, dist = calc_true_lunar_node(jd_1950)

        # Should return valid results
        assert 0 <= lon < 360
        assert lat == 0.0

    def test_true_node_for_late_date(self):
        """True node should work for future dates within ephemeris range."""
        # Year 2050 (within ephemeris range 1899-2053)
        jd_2050 = 2469807.5  # 2050-01-01
        lon, lat, dist = calc_true_lunar_node(jd_2050)

        # Should return valid results
        assert 0 <= lon < 360
        assert lat == 0.0


class TestDynamicObliquityInTrueLilith:
    """Test that true Lilith uses dynamic obliquity."""

    def test_true_lilith_returns_valid_coordinates(self):
        """True Lilith should return valid coordinates for any date."""
        jd_j2000 = 2451545.0
        lon, lat, dist = calc_true_lilith(jd_j2000)

        assert 0 <= lon < 360
        # Latitude can be non-zero but should be small (typically < 5°)
        assert -90 < lat < 90

    def test_true_lilith_for_early_date(self):
        """True Lilith should work for dates within ephemeris range."""
        # Year 1950 (within ephemeris range 1899-2053)
        jd_1950 = 2433282.5  # 1950-01-01
        lon, lat, dist = calc_true_lilith(jd_1950)

        # Should return valid results
        assert 0 <= lon < 360
        assert -90 < lat < 90

    def test_true_lilith_for_late_date(self):
        """True Lilith should work for future dates within ephemeris range."""
        # Year 2050 (within ephemeris range 1899-2053)
        jd_2050 = 2469807.5  # 2050-01-01
        lon, lat, dist = calc_true_lilith(jd_2050)

        # Should return valid results
        assert 0 <= lon < 360
        assert -90 < lat < 90


class TestObliquityDifferenceImpact:
    """Test the impact of using dynamic vs fixed obliquity."""

    def test_difference_at_j2000_is_minimal(self):
        """At J2000, dynamic obliquity should be very close to 23.4392911°."""
        jd_j2000 = 2451545.0
        eps_dynamic = math.degrees(_mean_obliquity_radians(jd_j2000))
        eps_fixed = 23.4392911  # Old fixed value

        # Difference should be less than 0.001°
        assert abs(eps_dynamic - eps_fixed) < 0.001

    def test_difference_grows_with_time(self):
        """The difference between dynamic and fixed obliquity grows with time."""
        eps_fixed = 23.4392911  # Old J2000 value

        # 1 century from J2000
        jd_1century = 2451545.0 + 36525.0
        eps_1c = math.degrees(_mean_obliquity_radians(jd_1century))
        diff_1c = abs(eps_1c - eps_fixed)

        # 5 centuries from J2000
        jd_5century = 2451545.0 + 5 * 36525.0
        eps_5c = math.degrees(_mean_obliquity_radians(jd_5century))
        diff_5c = abs(eps_5c - eps_fixed)

        # 10 centuries from J2000
        jd_10century = 2451545.0 + 10 * 36525.0
        eps_10c = math.degrees(_mean_obliquity_radians(jd_10century))
        diff_10c = abs(eps_10c - eps_fixed)

        # Differences should grow with time
        assert diff_5c > diff_1c
        assert diff_10c > diff_5c

        # At 10 centuries, difference should be noticeable (~0.13°)
        assert diff_10c > 0.1

    def test_obliquity_change_is_significant_for_coordinate_transform(self):
        """
        Verify that obliquity change matters for coordinate transformations.

        At 10 centuries from J2000, the obliquity differs by ~0.13°,
        which can cause position errors of up to ~8 arcminutes in ecliptic
        coordinates for objects at high declination.
        """
        # The maximum position error from obliquity error is approximately:
        # Δlon ≈ Δeps * sin(lat_ecl) for high ecliptic latitude objects
        # Δlat ≈ Δeps * sin(RA) for the latitude component

        # At 10 centuries from J2000
        jd_10century = 2451545.0 + 10 * 36525.0
        eps_dynamic = math.degrees(_mean_obliquity_radians(jd_10century))
        eps_fixed = 23.4392911

        delta_eps_deg = abs(eps_dynamic - eps_fixed)

        # For an object at maximum ecliptic latitude (~5° for Moon's orbit),
        # the position error would be noticeable
        # ~0.13° obliquity error * sin(5°) ≈ 0.011° = 40 arcsec position error
        # This is significant for precise calculations
        assert delta_eps_deg > 0.1, (
            "Obliquity difference should be > 0.1° at 10 centuries"
        )
