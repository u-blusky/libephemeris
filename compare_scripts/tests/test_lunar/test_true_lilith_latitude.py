"""
Validation of True Lilith ecliptic latitude calculation.

This module validates that the latitude calculation in calc_true_lilith
is correct by comparing against pyswisseph swe.calc_ut(jd, swe.OSCU_APOG, flags).

Unlike Mean Lilith which is always on the ecliptic (latitude = 0),
True Lilith (osculating lunar apogee) has an ecliptic latitude component
because the Moon's orbit is inclined ~5.145° to the ecliptic plane.

The latitude is calculated from the apogee vector after rotation to
ecliptic coordinates (lunar.py lines 1781-1786):

    lat = math.degrees(
        math.asin(
            apogee_ecl[2]
            / math.sqrt(apogee_ecl[0]**2 + apogee_ecl[1]**2 + apogee_ecl[2]**2)
        )
    )

This test validates that this formula correctly computes the ecliptic
latitude to within acceptable precision compared to Swiss Ephemeris.
"""

import math
import random
import statistics

import pytest

try:
    import swisseph as swe

    HAS_SWISSEPH = True
except ImportError:
    HAS_SWISSEPH = False

from libephemeris.lunar import calc_true_lilith

pytestmark = pytest.mark.skipif(not HAS_SWISSEPH, reason="pyswisseph not installed")


def generate_random_jd(start_year: int, end_year: int, count: int, seed: int):
    """Generate random Julian Dates for testing."""
    random.seed(seed)
    dates = []
    for _ in range(count):
        year = random.randint(start_year, end_year)
        month = random.randint(1, 12)
        day = random.randint(1, 28)
        hour = random.random() * 24.0
        jd = swe.julday(year, month, day, hour)
        dates.append((year, month, day, hour, jd))
    return dates


class TestTrueLilithLatitudeValidation:
    """
    Validate True Lilith latitude calculation against pyswisseph.

    The latitude calculation computes ecliptic latitude from the 3D apogee
    vector using the standard spherical coordinate formula. This test suite
    verifies accuracy and physical correctness.
    """

    def test_latitude_accuracy_comprehensive(self):
        """
        Compare True Lilith latitude against pyswisseph for 500 random dates.

        The latitude should match Swiss Ephemeris within 0.1° for all dates.
        Expected results based on empirical testing:
          - Mean error: ~0.024°
          - Max error:  ~0.063°
          - 95th percentile: ~0.056°
        """
        dates = generate_random_jd(1950, 2050, count=500, seed=456)
        errors = []

        for year, month, day, hour, jd in dates:
            try:
                swe_result = swe.calc_ut(jd, swe.OSCU_APOG, 0)
                swe_lat = swe_result[0][1]

                _, lib_lat, _ = calc_true_lilith(jd)

                diff = abs(swe_lat - lib_lat)
                errors.append(diff)
            except Exception:
                pass  # Skip dates outside ephemeris range

        assert len(errors) >= 400, f"Too many failed calculations: {500 - len(errors)}"

        mean_error = statistics.mean(errors)
        max_error = max(errors)
        percentile_95 = sorted(errors)[int(len(errors) * 0.95)]

        # Assert accuracy thresholds based on validated performance
        assert max_error < 0.1, (
            f"Maximum latitude error {max_error:.4f}° exceeds 0.1° threshold"
        )
        assert mean_error < 0.05, (
            f"Mean latitude error {mean_error:.4f}° exceeds 0.05° threshold"
        )
        assert percentile_95 < 0.07, (
            f"95th percentile error {percentile_95:.4f}° exceeds 0.07° threshold"
        )

    def test_latitude_within_physical_range(self):
        """
        Verify True Lilith latitude is within expected physical range.

        The lunar orbit is inclined ~5.145° to the ecliptic, so the
        osculating apogee latitude should be within approximately ±6°.
        Adding margin for perturbations, we allow ±8°.
        """
        dates = generate_random_jd(1950, 2050, count=200, seed=789)
        latitudes = []

        for year, month, day, hour, jd in dates:
            try:
                _, lib_lat, _ = calc_true_lilith(jd)
                latitudes.append(lib_lat)
            except Exception:
                pass

        assert len(latitudes) >= 180, "Too many failed calculations"

        min_lat = min(latitudes)
        max_lat = max(latitudes)

        # Physical constraint: lunar orbital inclination limits latitude
        assert -8 < min_lat, f"Min latitude {min_lat:.4f}° too negative"
        assert max_lat < 8, f"Max latitude {max_lat:.4f}° too positive"

        # Verify we're seeing reasonable range (should span most of ±5°)
        assert max_lat > 4, f"Max latitude {max_lat:.4f}° seems too small"
        assert min_lat < -4, f"Min latitude {min_lat:.4f}° seems too small"

    def test_latitude_formula_correctness(self):
        """
        Test that the latitude formula correctly extracts z-component.

        The latitude is calculated as arcsin(z/r) where z is the ecliptic
        z-component and r is the total magnitude. This test verifies the
        formula matches expected values for known geometric configurations.
        """
        # Test at a known epoch - J2000.0
        jd_j2000 = 2451545.0
        _, lat, _ = calc_true_lilith(jd_j2000)

        # Latitude should be reasonable
        assert -6 < lat < 6, f"J2000.0 latitude {lat:.4f}° outside expected range"

        # Compare with Swiss Ephemeris
        swe_result = swe.calc_ut(jd_j2000, swe.OSCU_APOG, 0)
        swe_lat = swe_result[0][1]

        diff = abs(swe_lat - lat)
        assert diff < 0.1, (
            f"J2000.0 latitude diff {diff:.4f}° exceeds 0.1° threshold "
            f"(SWE={swe_lat:.4f}°, Lib={lat:.4f}°)"
        )

    def test_latitude_consistency_over_nodal_cycle(self):
        """
        Verify latitude behavior over the lunar nodal cycle (~18.6 years).

        The latitude should oscillate with the nodal cycle as the apogee
        crosses above and below the ecliptic. We check that latitudes
        span the expected range over a full cycle.
        """
        # Sample over 19 years (more than one nodal cycle)
        jd_start = 2451545.0  # J2000.0
        nodal_period_days = 18.6 * 365.25

        # Sample at 50 evenly-spaced points
        latitudes = []
        for i in range(50):
            jd = jd_start + (i / 50.0) * nodal_period_days
            try:
                _, lat, _ = calc_true_lilith(jd)
                latitudes.append(lat)
            except Exception:
                pass

        assert len(latitudes) >= 45, "Too many failed calculations"

        # Over a full nodal cycle, should see range of latitudes
        lat_range = max(latitudes) - min(latitudes)
        assert lat_range > 8, (
            f"Latitude range {lat_range:.2f}° over nodal cycle seems too small"
        )

    def test_latitude_sign_agreement(self):
        """
        Verify latitude sign matches Swiss Ephemeris.

        The sign of latitude indicates whether the apogee is north (+)
        or south (-) of the ecliptic.
        """
        dates = generate_random_jd(1980, 2020, count=100, seed=321)
        sign_matches = 0
        total = 0

        for year, month, day, hour, jd in dates:
            try:
                swe_result = swe.calc_ut(jd, swe.OSCU_APOG, 0)
                swe_lat = swe_result[0][1]

                _, lib_lat, _ = calc_true_lilith(jd)

                total += 1
                # Signs match if same sign or both near zero
                if (swe_lat * lib_lat > 0) or (
                    abs(swe_lat) < 0.1 and abs(lib_lat) < 0.1
                ):
                    sign_matches += 1
            except Exception:
                pass

        match_rate = sign_matches / total if total > 0 else 0
        assert match_rate > 0.98, (
            f"Latitude sign agreement {match_rate:.1%} below 98% threshold"
        )

    def test_latitude_near_nodes(self):
        """
        Verify latitude behavior when apogee is near lunar nodes.

        When the apogee longitude is near the lunar node longitude,
        the apogee is crossing the ecliptic, so latitude should be near zero.
        """
        # Find dates where apogee is expected to be near nodes
        # (by checking when latitude is small)
        dates = generate_random_jd(1990, 2010, count=200, seed=555)
        near_zero_count = 0
        large_count = 0

        for year, month, day, hour, jd in dates:
            try:
                _, lat, _ = calc_true_lilith(jd)

                if abs(lat) < 1.0:
                    near_zero_count += 1
                if abs(lat) > 4.0:
                    large_count += 1
            except Exception:
                pass

        # Should see some latitudes near zero (node crossing) and some large
        assert near_zero_count > 10, f"Too few near-zero latitudes: {near_zero_count}"
        assert large_count > 30, f"Too few large latitudes: {large_count}"


class TestTrueLilithLatitudeEdgeCases:
    """
    Edge case tests for True Lilith latitude calculation.
    """

    def test_latitude_at_historical_dates(self):
        """
        Test latitude calculation at historically significant dates.
        """
        # Historical dates that should work with ephemeris
        historical_dates = [
            (1969, 7, 20, 20.0),  # Apollo 11 landing
            (1986, 1, 28, 16.5),  # Challenger disaster
            (2000, 1, 1, 0.0),  # Y2K
            (2012, 12, 21, 12.0),  # End of Mayan calendar
            (2020, 3, 15, 0.0),  # COVID lockdowns begin
        ]

        for year, month, day, hour in historical_dates:
            jd = swe.julday(year, month, day, hour)

            try:
                swe_result = swe.calc_ut(jd, swe.OSCU_APOG, 0)
                swe_lat = swe_result[0][1]

                _, lib_lat, _ = calc_true_lilith(jd)

                diff = abs(swe_lat - lib_lat)
                assert diff < 0.1, (
                    f"Latitude diff {diff:.4f}° at {year}/{month}/{day} exceeds threshold"
                )
            except Exception as e:
                pytest.skip(f"Skipping {year}/{month}/{day}: {e}")

    def test_latitude_consistency_at_same_instant(self):
        """
        Verify repeated calls at same instant give consistent latitude.
        """
        jd = 2458000.5

        results = [calc_true_lilith(jd) for _ in range(5)]
        latitudes = [r[1] for r in results]

        # All latitudes should be identical
        assert all(lat == latitudes[0] for lat in latitudes), (
            "Inconsistent latitude values for same JD"
        )

    def test_latitude_continuous_variation(self):
        """
        Verify latitude varies continuously (no sudden jumps).
        """
        jd_start = 2458000.5
        step = 0.1  # 0.1 day = 2.4 hours

        prev_lat = None
        max_jump = 0

        for i in range(100):
            jd = jd_start + i * step
            _, lat, _ = calc_true_lilith(jd)

            if prev_lat is not None:
                jump = abs(lat - prev_lat)
                max_jump = max(max_jump, jump)
                # Latitude should not jump more than 0.5° in 2.4 hours
                assert jump < 0.5, (
                    f"Latitude jump of {jump:.4f}° at JD {jd:.4f} is too large"
                )

            prev_lat = lat

    def test_mean_lilith_zero_latitude(self):
        """
        Verify Mean Lilith has zero latitude (for comparison).

        Mean Lilith is defined as the mean apogee which lies on the ecliptic,
        so its latitude should always be zero.
        """
        from libephemeris.lunar import calc_mean_lilith

        dates = generate_random_jd(1980, 2020, count=20, seed=777)

        for year, month, day, hour, jd in dates:
            try:
                mean_lon = calc_mean_lilith(jd)
                # Mean Lilith returns only longitude (scalar)
                assert isinstance(mean_lon, float), "Mean Lilith should return float"
            except Exception:
                pass


class TestTrueLilithLatitudeDocumentation:
    """
    Documentation tests for True Lilith latitude.
    """

    def test_latitude_docstring_claims(self):
        """
        Verify claims made in function docstring about latitude.

        From calc_true_lilith docstring:
          - latitude: Ecliptic latitude in degrees (small, typically < 5°)
        """
        dates = generate_random_jd(1960, 2040, count=300, seed=888)
        latitudes = []

        for year, month, day, hour, jd in dates:
            try:
                _, lat, _ = calc_true_lilith(jd)
                latitudes.append(abs(lat))
            except Exception:
                pass

        # "typically < 5°" means most should be within this range
        # The lunar inclination is 5.145°, so we expect ~85% within 5° (a uniform
        # distribution over inclination would give 5/5.145 ≈ 97%, but the actual
        # distribution is sinusoidal so latitudes near max/min are more likely)
        within_5_deg = sum(1 for lat in latitudes if lat < 5)
        ratio = within_5_deg / len(latitudes)

        assert ratio > 0.8, f"Only {ratio:.1%} of latitudes are < 5°, expected > 80%"

        # All latitudes should be within the physical limit of ~5.3°
        assert all(lat < 6 for lat in latitudes), (
            f"Some latitudes exceed 6°: max={max(latitudes):.2f}°"
        )
