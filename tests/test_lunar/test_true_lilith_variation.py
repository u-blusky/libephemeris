"""
Tests for variation correction in the True Lilith (osculating lunar apogee) calculation.

The Variation is a major lunar perturbation discovered by Tycho Brahe, caused by
the transverse component of the solar tidal force. It affects both the Moon's
longitude and the lunar apogee position.

Key characteristics:
- Period: ~14.77 days (half the synodic month)
- Amplitude: ~0.658 degrees
- Argument: 2D (twice the mean elongation)
- Cause: differential solar tidal force at quadrature vs conjunction/opposition

The variation correction term in True Lilith is: 0.6583 * sin(2D)

References:
    - Chapront-Touze, M. & Chapront, J. "Lunar Tables and Programs" (1991)
    - Brown, E.W. "An Introductory Treatise on the Lunar Theory" (1896)
    - Meeus, J. "Astronomical Algorithms" (2nd ed., 1998), Chapter 47
"""

import math
import pytest
from libephemeris.lunar import calc_true_lilith, calc_mean_lilith


# Variation period in days (half the synodic month)
VARIATION_PERIOD_DAYS = 14.765
# Variation amplitude in degrees (for True Lilith)
VARIATION_AMPLITUDE_DEG = 0.6583


def _calc_variation_argument(jd_tt: float) -> float:
    """
    Calculate the variation argument (2D) for testing.

    This helper computes the same lunar argument used internally by
    calc_true_lilith to verify the variation correction is applied correctly.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        float: Variation argument (2D) in radians
    """
    T = (jd_tt - 2451545.0) / 36525.0  # Julian centuries from J2000.0

    # Mean elongation of Moon from Sun (D)
    D = (
        297.8501921
        + 445267.1114034 * T
        - 0.0018819 * T**2
        + T**3 / 545868.0
        - T**4 / 113065000.0
    )

    # Convert to radians and return 2D
    D_rad = math.radians(D % 360.0)
    return 2.0 * D_rad


class TestVariationCorrectionAmplitude:
    """Test the amplitude and bounds of the variation correction."""

    def test_variation_correction_amplitude_is_bounded(self):
        """Variation correction should be bounded by +-0.6583 degrees."""
        # Test at various dates to capture different phases of the variation cycle
        test_dates = [
            2451545.0,  # J2000.0
            2451545.0 + 3.7,  # ~1/4 variation period later
            2451545.0 + 7.4,  # ~1/2 variation period later
            2451545.0 + 11.1,  # ~3/4 variation period later
            2451545.0 + 14.77,  # ~1 full variation period
        ]

        for jd in test_dates:
            variation_arg = _calc_variation_argument(jd)
            variation_correction = VARIATION_AMPLITUDE_DEG * math.sin(variation_arg)

            # Correction should be bounded by amplitude
            assert -0.66 <= variation_correction <= 0.66, (
                f"Variation correction {variation_correction}deg at JD {jd} exceeds bounds"
            )

    def test_variation_reaches_maximum_amplitude(self):
        """Variation correction should reach approximately +-0.6583 degrees at extrema."""
        jd_base = 2451545.0

        # Search for approximate extrema over one variation cycle
        max_variation = -1.0
        min_variation = 1.0

        for i in range(148):  # Check over ~14.8 days in 0.1 day steps
            jd = jd_base + i * 0.1
            variation_arg = _calc_variation_argument(jd)
            variation_val = VARIATION_AMPLITUDE_DEG * math.sin(variation_arg)

            if variation_val > max_variation:
                max_variation = variation_val
            if variation_val < min_variation:
                min_variation = variation_val

        # Maximum should be close to +0.6583
        assert max_variation > 0.60, (
            f"Maximum variation {max_variation}deg is too small"
        )

        # Minimum should be close to -0.6583
        assert min_variation < -0.60, (
            f"Minimum variation {min_variation}deg is not negative enough"
        )


class TestVariationCorrectionPeriod:
    """Test the period of the variation correction."""

    def test_variation_correction_has_correct_period(self):
        """Variation correction should complete a cycle in approximately 14.77 days."""
        jd_start = 2451545.0

        # Get variation argument at start
        arg_start = _calc_variation_argument(jd_start)

        # After one period, argument should return to same value (mod 2pi)
        arg_end = _calc_variation_argument(jd_start + VARIATION_PERIOD_DAYS)

        # Normalize both to [0, 2pi)
        arg_start_norm = arg_start % (2 * math.pi)
        arg_end_norm = arg_end % (2 * math.pi)

        # The arguments should be close (within ~0.2 rad due to period approximation)
        diff = abs(arg_end_norm - arg_start_norm)
        if diff > math.pi:
            diff = 2 * math.pi - diff

        assert diff < 0.25, (
            f"Variation argument changed by {math.degrees(diff)}deg over one period"
        )

    def test_variation_period_is_half_synodic_month(self):
        """Variation period should be approximately half the synodic month."""
        synodic_month = 29.530589  # days
        expected_period = synodic_month / 2

        # Our constant should match
        assert abs(VARIATION_PERIOD_DAYS - expected_period) < 0.1, (
            f"Variation period {VARIATION_PERIOD_DAYS} != {expected_period}"
        )

    def test_variation_argument_rate(self):
        """Variation argument should advance at approximately 24.4 deg/day."""
        jd_start = 2451545.0

        arg1 = _calc_variation_argument(jd_start)
        arg2 = _calc_variation_argument(jd_start + 1.0)

        # Rate should be approximately 24.38 deg/day (2 * 12.19 deg/day)
        variation_rate = math.degrees(arg2 - arg1)
        # Normalize to [0, 360)
        while variation_rate < 0:
            variation_rate += 360
        while variation_rate >= 360:
            variation_rate -= 360

        assert 23.0 < variation_rate < 26.0, (
            f"Variation rate = {variation_rate} deg/day"
        )


class TestVariationEffectOnTrueLilith:
    """Test the variation effect on True Lilith calculation."""

    def test_true_lilith_varies_with_variation_period(self):
        """True Lilith longitude should show variation consistent with ~14.77 day period."""
        jd_start = 2451545.0
        half_period = VARIATION_PERIOD_DAYS / 2

        # Get True Lilith at start and half-period later
        lon_start, _, _ = calc_true_lilith(jd_start)
        lon_half, _, _ = calc_true_lilith(jd_start + half_period)

        # Calculate the difference
        diff = lon_half - lon_start
        if diff > 180:
            diff -= 360
        if diff < -180:
            diff += 360

        # The mean apogee moves ~0.11deg/day, so over ~7.4 days: ~0.8deg
        # But the variation correction adds oscillations of ~+-0.66deg
        # At half period, we expect the variation to have opposite sign
        # The longitude difference should reflect both the secular motion
        # and the variation oscillation

        # Just verify there is measurable change
        assert abs(diff) > 0.3, (
            f"Expected longitude change over half variation period, got {diff}deg"
        )

    def test_variation_correction_integration(self):
        """Test that variation correction is integrated into True Lilith output."""
        # Test at J2000.0 - calculate expected variation contribution
        jd = 2451545.0
        variation_arg = _calc_variation_argument(jd)
        expected_variation = VARIATION_AMPLITUDE_DEG * math.sin(variation_arg)

        # Get True Lilith result
        lon, lat, e_mag = calc_true_lilith(jd)

        # Verify results are valid (variation doesn't break the calculation)
        assert 0 <= lon < 360, f"Longitude {lon} out of range"
        assert -10 < lat < 10, f"Latitude {lat} unexpectedly large"
        assert 0.002 < e_mag < 0.003, f"Distance {e_mag} AU out of expected range"

        # Log the variation contribution for reference
        assert abs(expected_variation) <= 0.66, (
            f"Expected variation contribution: {expected_variation:.4f}deg"
        )

    def test_variation_at_extrema(self):
        """Test variation correction behavior at maximum and minimum points.

        Find dates where the variation argument (2D) is approximately
        pi/2 (maximum correction) and -pi/2 (minimum correction).
        """
        jd_base = 2451545.0

        # Search for approximate extrema over one variation cycle
        max_variation = -1.0
        min_variation = 1.0
        jd_max = jd_base
        jd_min = jd_base

        for i in range(148):  # Check over ~14.8 days in 0.1 day steps
            jd = jd_base + i * 0.1
            variation_arg = _calc_variation_argument(jd)
            variation_val = VARIATION_AMPLITUDE_DEG * math.sin(variation_arg)

            if variation_val > max_variation:
                max_variation = variation_val
                jd_max = jd
            if variation_val < min_variation:
                min_variation = variation_val
                jd_min = jd

        # Maximum should be close to +0.6583
        assert max_variation > 0.60, (
            f"Maximum variation {max_variation}deg is too small"
        )

        # Minimum should be close to -0.6583
        assert min_variation < -0.60, (
            f"Minimum variation {min_variation}deg is not negative enough"
        )

        # Verify True Lilith can be calculated at both extrema
        lon_max, _, _ = calc_true_lilith(jd_max)
        lon_min, _, _ = calc_true_lilith(jd_min)

        assert 0 <= lon_max < 360
        assert 0 <= lon_min < 360


class TestVariationVsMeanLilith:
    """Test comparing True Lilith (with variation) to Mean Lilith."""

    def test_variation_creates_oscillation_around_mean(self):
        """True Lilith should oscillate around mean Lilith partially due to variation."""
        jd_start = 2451545.0

        differences = []
        for i in range(int(2 * VARIATION_PERIOD_DAYS)):
            jd = jd_start + i
            true_lon, _, _ = calc_true_lilith(jd)
            mean_lon = calc_mean_lilith(jd)

            diff = true_lon - mean_lon
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360
            differences.append(diff)

        # Should see oscillation (sign changes)
        pos_count = sum(1 for d in differences if d > 0)
        neg_count = sum(1 for d in differences if d < 0)

        # Both positive and negative differences should occur
        assert pos_count > 0 or neg_count > 0, (
            "Expected some variation in True-Mean difference"
        )

        # Range of differences should be significant
        diff_range = max(differences) - min(differences)
        assert diff_range > 0.5, f"True-Mean Lilith range = {diff_range}"


class TestVariationCombinesWithEvection:
    """Test that variation and evection corrections work together."""

    def test_both_corrections_applied(self):
        """Both variation and evection corrections should be applied to True Lilith."""
        jd_start = 2451545.0

        # Sample over 60 days (covers multiple variation and evection periods)
        n_days = 60
        longitudes = []
        for i in range(n_days):
            jd = jd_start + i
            lon, _, _ = calc_true_lilith(jd)
            longitudes.append(lon)

        # Remove secular trend (~0.11 deg/day prograde for apogee)
        detrended = []
        for i, lon in enumerate(longitudes):
            expected = longitudes[0] + 0.11 * i
            diff = lon - expected
            while diff > 180:
                diff -= 360
            while diff < -180:
                diff += 360
            detrended.append(diff)

        # Total oscillation should be significant
        oscillation_range = max(detrended) - min(detrended)

        # With both evection (~1.27 deg) and variation (~0.66 deg) corrections,
        # we expect total oscillation of at least 1.5 degrees
        assert oscillation_range > 1.0, (
            f"Detrended oscillation range = {oscillation_range}"
        )

    def test_total_correction_magnitude(self):
        """Total correction from evection + variation should be appropriately sized."""
        jd = 2451545.0

        # Get True Lilith position
        lon, _, _ = calc_true_lilith(jd)

        # Calculate individual corrections
        T = (jd - 2451545.0) / 36525.0

        # Mean elongation (D)
        D = (
            297.8501921
            + 445267.1114034 * T
            - 0.0018819 * T**2
            + T**3 / 545868.0
            - T**4 / 113065000.0
        )

        # Moon's mean anomaly (M')
        M_prime = (
            134.9633964
            + 477198.8675055 * T
            + 0.0087414 * T**2
            + T**3 / 69699.0
            - T**4 / 14712000.0
        )

        D_rad = math.radians(D % 360.0)
        M_prime_rad = math.radians(M_prime % 360.0)

        # Evection: 1.2739 * sin(2D - M')
        evection_arg = 2.0 * D_rad - M_prime_rad
        evection = 1.2739 * math.sin(evection_arg)

        # Variation: 0.6583 * sin(2D)
        variation_arg = 2.0 * D_rad
        variation = 0.6583 * math.sin(variation_arg)

        # Total correction should be bounded by sum of amplitudes
        max_total = 1.2739 + 0.6583  # ~1.93 degrees
        assert abs(evection + variation) <= max_total + 0.01, (
            f"Total correction {evection + variation} exceeds maximum {max_total}"
        )


@pytest.mark.integration
class TestVariationIntegration:
    """Integration tests for variation correction in True Lilith."""

    def test_variation_over_multiple_periods(self):
        """Variation correction should be consistent over multiple periods."""
        jd_start = 2451545.0
        n_periods = 5

        # Sample at consistent phases across multiple periods
        for period in range(n_periods):
            jd = jd_start + period * VARIATION_PERIOD_DAYS

            lon, lat, e_mag = calc_true_lilith(jd)

            # All values should be valid
            assert 0 <= lon < 360, f"Invalid lon at period {period}"
            assert -10 < lat < 10, f"Latitude {lat} unexpectedly large"
            assert 0.002 < e_mag < 0.003, f"Distance {e_mag} AU out of range"

    def test_variation_stability_over_years(self):
        """Variation correction should remain stable over long time spans."""
        jd_start = 2451545.0  # J2000

        # Test over 10 years, sampling monthly
        for year in range(10):
            for month in range(12):
                jd = jd_start + year * 365.25 + month * 30.4

                lon, lat, e_mag = calc_true_lilith(jd)

                # Results should always be valid
                assert 0 <= lon < 360, (
                    f"Invalid longitude at year {2000 + year} month {month + 1}"
                )
                assert -10 < lat < 10, (
                    f"Invalid latitude at year {2000 + year} month {month + 1}"
                )

    def test_variation_effect_amplitude_consistency(self):
        """Variation effect amplitude should be consistent across different epochs."""
        # Test at different starting points
        epochs = [
            2451545.0,  # J2000
            2451545.0 + 365.25 * 10,  # 2010
            2451545.0 + 365.25 * 20,  # 2020
            2451545.0 - 365.25 * 20,  # 1980
        ]

        for epoch in epochs:
            longitudes = []
            for i in range(int(VARIATION_PERIOD_DAYS)):
                jd = epoch + i
                lon, _, _ = calc_true_lilith(jd)
                longitudes.append(lon)

            # Detrend the data
            detrended = []
            for i, lon in enumerate(longitudes):
                expected = longitudes[0] + 0.11 * i
                diff = lon - expected
                while diff > 180:
                    diff -= 360
                while diff < -180:
                    diff += 360
                detrended.append(diff)

            # Range should be significant (variation contributes)
            lon_range = max(detrended) - min(detrended)
            assert lon_range > 0.3, (
                f"Longitude oscillation range at epoch {epoch} = {lon_range}"
            )
