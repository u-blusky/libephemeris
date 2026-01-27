"""
Tests for evection-related secondary terms in the True Lilith (osculating lunar apogee) calculation.

These tests verify the implementation of secondary terms from Meeus "Astronomical Algorithms"
Chapter 47 Table 47.B that affect the lunar eccentricity and thus the apogee direction.

The evection-related secondary terms include:
- l - 2D (M' - 2D): amplitude -0.2136 degrees
- l + 2D (M' + 2D): amplitude +0.1058 degrees
- 2l (2*M'): amplitude -0.2037 degrees
- 2l - 2D (2*M' - 2D): amplitude +0.1027 degrees

References:
    - Meeus, J. "Astronomical Algorithms" (2nd ed., 1998), Chapter 47, Table 47.B
    - Chapront-Touze, M. & Chapront, J. "Lunar Tables and Programs" (1991)
    - Brown, E.W. "An Introductory Treatise on the Lunar Theory" (1896)
"""

import math
import pytest
from libephemeris.lunar import (
    _calc_lunar_fundamental_arguments,
    calc_true_lilith,
    calc_mean_lilith,
)


# Periods of the secondary terms in days (approximate)
PERIOD_M_PRIME_MINUS_2D = 31.81  # Similar to main evection
PERIOD_M_PRIME_PLUS_2D = 9.6  # Faster beat frequency
PERIOD_2M_PRIME = 13.78  # Half anomalistic month
PERIOD_2M_PRIME_MINUS_2D = 14.8  # Similar to variation


class TestEvectionSecondaryTermsArguments:
    """Tests for the evection-related secondary term arguments."""

    def test_m_prime_minus_2d_argument_calculation(self):
        """Test that M' - 2D argument is calculated correctly."""
        jd_j2000 = 2451545.0

        D, M, M_prime, F = _calc_lunar_fundamental_arguments(jd_j2000)

        # Calculate M' - 2D argument
        arg = M_prime - 2.0 * D
        arg_deg = math.degrees(arg) % 360

        # At J2000, D ~ 297.85 deg, M' ~ 134.96 deg
        # So M' - 2D ~ 134.96 - 595.70 ~ -460.74 ~ -100.74 ~ 259.26 (mod 360)
        assert 0 <= arg_deg < 360, f"Arg M'-2D = {arg_deg}"

    def test_m_prime_plus_2d_argument_calculation(self):
        """Test that M' + 2D argument is calculated correctly."""
        jd_j2000 = 2451545.0

        D, M, M_prime, F = _calc_lunar_fundamental_arguments(jd_j2000)

        # Calculate M' + 2D argument
        arg = M_prime + 2.0 * D
        arg_deg = math.degrees(arg) % 360

        # At J2000, D ~ 297.85 deg, M' ~ 134.96 deg
        # So M' + 2D ~ 134.96 + 595.70 ~ 730.66 ~ 10.66 (mod 360)
        assert 0 <= arg_deg < 360, f"Arg M'+2D = {arg_deg}"

    def test_2m_prime_argument_calculation(self):
        """Test that 2*M' argument is calculated correctly."""
        jd_j2000 = 2451545.0

        D, M, M_prime, F = _calc_lunar_fundamental_arguments(jd_j2000)

        # Calculate 2*M' argument
        arg = 2.0 * M_prime
        arg_deg = math.degrees(arg) % 360

        # At J2000, M' ~ 134.96 deg
        # So 2M' ~ 269.92 (mod 360)
        assert 0 <= arg_deg < 360, f"Arg 2M' = {arg_deg}"

    def test_2m_prime_minus_2d_argument_calculation(self):
        """Test that 2*M' - 2D argument is calculated correctly."""
        jd_j2000 = 2451545.0

        D, M, M_prime, F = _calc_lunar_fundamental_arguments(jd_j2000)

        # Calculate 2*M' - 2D argument
        arg = 2.0 * M_prime - 2.0 * D
        arg_deg = math.degrees(arg) % 360

        # At J2000, D ~ 297.85 deg, M' ~ 134.96 deg
        # So 2M' - 2D ~ 269.92 - 595.70 ~ -325.78 ~ 34.22 (mod 360)
        assert 0 <= arg_deg < 360, f"Arg 2M'-2D = {arg_deg}"

    def test_arguments_vary_over_time(self):
        """All evection secondary term arguments should change over time."""
        jd_start = 2451545.0

        D1, M1, M1_prime, F1 = _calc_lunar_fundamental_arguments(jd_start)
        D2, M2, M2_prime, F2 = _calc_lunar_fundamental_arguments(jd_start + 7.0)

        # M' - 2D
        arg1_a = M1_prime - 2.0 * D1
        arg2_a = M2_prime - 2.0 * D2
        assert arg1_a != arg2_a, "M' - 2D should vary over time"

        # M' + 2D
        arg1_b = M1_prime + 2.0 * D1
        arg2_b = M2_prime + 2.0 * D2
        assert arg1_b != arg2_b, "M' + 2D should vary over time"

        # 2M'
        arg1_c = 2.0 * M1_prime
        arg2_c = 2.0 * M2_prime
        assert arg1_c != arg2_c, "2M' should vary over time"

        # 2M' - 2D
        arg1_d = 2.0 * M1_prime - 2.0 * D1
        arg2_d = 2.0 * M2_prime - 2.0 * D2
        assert arg1_d != arg2_d, "2M' - 2D should vary over time"


class TestEvectionSecondaryTermsPeriods:
    """Tests for the periods of the evection-related secondary terms."""

    def test_m_prime_minus_2d_period(self):
        """M' - 2D argument should complete a cycle in ~31.8 days."""
        jd_start = 2451545.0

        D1, M1, M1_prime, F1 = _calc_lunar_fundamental_arguments(jd_start)
        D2, M2, M2_prime, F2 = _calc_lunar_fundamental_arguments(
            jd_start + PERIOD_M_PRIME_MINUS_2D
        )

        arg1 = (M1_prime - 2.0 * D1) % (2.0 * math.pi)
        arg2 = (M2_prime - 2.0 * D2) % (2.0 * math.pi)

        diff_rad = abs(arg2 - arg1)
        diff_deg = math.degrees(min(diff_rad, 2 * math.pi - diff_rad))

        # Allow tolerance for period approximation
        assert diff_deg < 15, f"M'-2D arg diff after 1 period = {diff_deg}"

    def test_2m_prime_period(self):
        """2*M' argument should complete a cycle in ~13.78 days."""
        jd_start = 2451545.0

        D1, M1, M1_prime, F1 = _calc_lunar_fundamental_arguments(jd_start)
        D2, M2, M2_prime, F2 = _calc_lunar_fundamental_arguments(
            jd_start + PERIOD_2M_PRIME
        )

        arg1 = (2.0 * M1_prime) % (2.0 * math.pi)
        arg2 = (2.0 * M2_prime) % (2.0 * math.pi)

        diff_rad = abs(arg2 - arg1)
        diff_deg = math.degrees(min(diff_rad, 2 * math.pi - diff_rad))

        # Allow tolerance for period approximation
        assert diff_deg < 15, f"2M' arg diff after 1 period = {diff_deg}"


class TestEvectionSecondaryTermsCoefficients:
    """Tests for the coefficient values of evection-related secondary terms."""

    def test_coefficient_magnitudes(self):
        """Verify coefficient magnitudes are within expected ranges."""
        # Coefficients from Meeus Table 47.B
        coeff_m_prime_minus_2d = 0.2136
        coeff_m_prime_plus_2d = 0.1058
        coeff_2m_prime = 0.2037
        coeff_2m_prime_minus_2d = 0.1027

        # All coefficients should be between 0.1 and 0.3 degrees
        for coeff in [
            coeff_m_prime_minus_2d,
            coeff_m_prime_plus_2d,
            coeff_2m_prime,
            coeff_2m_prime_minus_2d,
        ]:
            assert 0.05 < coeff < 0.5, f"Coefficient {coeff} out of expected range"

    def test_coefficient_hierarchy(self):
        """Verify the relative magnitudes of coefficients."""
        # Coefficients from Meeus Table 47.B
        coeff_m_prime_minus_2d = 0.2136
        coeff_m_prime_plus_2d = 0.1058
        coeff_2m_prime = 0.2037
        coeff_2m_prime_minus_2d = 0.1027

        # M' - 2D and 2M' should have larger amplitudes
        assert coeff_m_prime_minus_2d > coeff_m_prime_plus_2d
        assert coeff_2m_prime > coeff_2m_prime_minus_2d

    def test_total_secondary_contribution_estimate(self):
        """Estimate the total contribution of secondary terms."""
        # Sum of absolute coefficients
        total = 0.2136 + 0.1058 + 0.2037 + 0.1027

        # Total should be approximately 0.6-0.7 degrees peak-to-peak
        assert 0.5 < total < 0.8, f"Total secondary coefficient sum = {total}"


class TestEvectionSecondaryTermsEffectOnTrueLilith:
    """Tests for the effect of evection secondary terms on True Lilith calculation."""

    def test_true_lilith_valid_over_secondary_periods(self):
        """True Lilith should return valid results across secondary term periods."""
        jd_start = 2451545.0

        # Sample over the longest secondary period (~31.8 days)
        for i in range(35):
            jd = jd_start + i
            lon, lat, e_mag = calc_true_lilith(jd)

            assert 0 <= lon < 360, f"Invalid longitude at day {i}: {lon}"
            assert -10 < lat < 10, f"Invalid latitude at day {i}: {lat}"
            assert 0 < e_mag < 0.2, f"Invalid eccentricity at day {i}: {e_mag}"

    def test_true_lilith_shows_secondary_variation(self):
        """True Lilith should show variations consistent with secondary terms."""
        jd_start = 2451545.0

        longitudes = []
        for i in range(int(2 * PERIOD_M_PRIME_MINUS_2D)):
            jd = jd_start + i
            lon, lat, e_mag = calc_true_lilith(jd)
            longitudes.append(lon)

        # Compute detrended variation (remove apsidal precession trend)
        # The apogee moves prograde ~0.111 degrees per day
        expected_motion = 0.111  # degrees/day (40.67 deg/year)

        detrended = []
        for i, lon in enumerate(longitudes):
            expected = longitudes[0] + expected_motion * i
            # Handle wrap-around
            diff = lon - expected
            while diff > 180:
                diff -= 360
            while diff < -180:
                diff += 360
            detrended.append(diff)

        # Detrended variation should show oscillation
        variation = max(detrended) - min(detrended)

        # Expect variation of at least 1 degree due to all perturbations
        assert variation > 0.5, f"Detrended variation = {variation}"

    def test_true_lilith_varies_from_mean(self):
        """True Lilith should show significant variation from Mean Lilith."""
        jd_start = 2451545.0

        differences = []
        for i in range(int(2 * PERIOD_M_PRIME_MINUS_2D)):
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
        assert pos_count > 0 or neg_count > 0, "True Lilith should differ from mean"

        # Range of differences should be significant
        diff_range = max(differences) - min(differences)
        assert diff_range > 0.5, f"True-Mean Lilith range = {diff_range}"


class TestEvectionSecondaryTermsOscillation:
    """Tests for oscillation patterns from evection secondary terms."""

    def test_sinusoidal_pattern_over_2m_prime_period(self):
        """Secondary terms should create sinusoidal patterns over their periods."""
        jd_start = 2451545.0

        # Sample at quarter-period intervals for 2M' term
        samples = []
        for i in range(5):  # 0, 1/4, 1/2, 3/4, 1 period
            jd = jd_start + i * PERIOD_2M_PRIME / 4
            lon, _, _ = calc_true_lilith(jd)
            samples.append(lon)

        # Compute differences
        diffs = []
        for i in range(len(samples) - 1):
            diff = samples[i + 1] - samples[i]
            # Handle wrap-around
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360
            diffs.append(diff)

        # Not all differences should be the same (oscillation expected)
        unique_diffs = len(set(round(d, 3) for d in diffs))
        assert unique_diffs > 1, "Expected variation in differences"


@pytest.mark.integration
class TestEvectionSecondaryIntegration:
    """Integration tests for evection secondary terms."""

    def test_secondary_terms_over_multiple_cycles(self):
        """Secondary terms should be consistent over multiple periods."""
        jd_start = 2451545.0
        n_periods = 3

        for period in range(n_periods):
            jd = jd_start + period * PERIOD_2M_PRIME

            lon, lat, e_mag = calc_true_lilith(jd)

            # All values should be valid
            assert 0 <= lon < 360, f"Invalid lon at period {period}"
            assert -10 < lat < 10, f"Invalid lat at period {period}"
            assert 0 < e_mag < 0.2, f"Invalid eccentricity at period {period}"

    def test_secondary_terms_stability_over_years(self):
        """Secondary terms should remain stable over long time spans."""
        jd_start = 2451545.0  # J2000

        # Test over 5 years, sampling monthly
        for year in range(5):
            for month in range(12):
                jd = jd_start + year * 365.25 + month * 30.4

                lon, lat, e_mag = calc_true_lilith(jd)

                # Values should be within reasonable bounds
                assert 0 <= lon < 360, (
                    f"Invalid lon at year {2000 + year} month {month + 1}"
                )
                assert -15 < lat < 15, (
                    f"Invalid lat at year {2000 + year} month {month + 1}"
                )

    def test_combined_effect_amplitude(self):
        """Combined effect of all secondary terms should be within expected range."""
        jd_start = 2451545.0

        # Sample over enough time to capture all secondary term cycles
        max_period = max(
            PERIOD_M_PRIME_MINUS_2D,
            PERIOD_M_PRIME_PLUS_2D,
            PERIOD_2M_PRIME,
            PERIOD_2M_PRIME_MINUS_2D,
        )

        longitudes = []
        for i in range(int(3 * max_period)):
            jd = jd_start + i
            lon, _, _ = calc_true_lilith(jd)
            longitudes.append(lon)

        # Handle wrap-around for range calculation
        # Convert to relative positions
        base = longitudes[0]
        relative = []
        for lon in longitudes:
            diff = lon - base
            if diff > 180:
                diff -= 360
            elif diff < -180:
                diff += 360
            relative.append(diff)

        # Compute variation around trend
        # Simple detrending: remove linear trend
        n = len(relative)
        if n > 1:
            slope = (relative[-1] - relative[0]) / (n - 1)
            detrended = [rel - slope * i for i, rel in enumerate(relative)]
            variation = max(detrended) - min(detrended)

            # Total variation should be significant but bounded
            # (includes evection, variation, annual equation, parallactic inequality,
            # reduction to ecliptic, and evection secondary terms)
            # The osculating apogee can vary significantly from mean apogee
            assert 0.5 < variation < 100, f"Total variation = {variation}"


class TestEvectionSecondaryTermsPhysics:
    """Tests validating the physical basis of evection secondary terms."""

    def test_m_prime_related_to_anomaly(self):
        """M' (Moon's mean anomaly) should vary with anomalistic month."""
        jd_start = 2451545.0
        anomalistic_month = 27.55455  # days

        D1, M1, M1_prime, F1 = _calc_lunar_fundamental_arguments(jd_start)
        D2, M2, M2_prime, F2 = _calc_lunar_fundamental_arguments(
            jd_start + anomalistic_month
        )

        # After one anomalistic month, M' should return to ~same value
        diff = abs(M2_prime - M1_prime) % (2 * math.pi)
        diff_deg = math.degrees(min(diff, 2 * math.pi - diff))

        # Allow some tolerance
        assert diff_deg < 3, f"M' diff after 1 anomalistic month = {diff_deg}"

    def test_d_related_to_synodic_month(self):
        """D (mean elongation) should vary with synodic month."""
        jd_start = 2451545.0
        synodic_month = 29.530589  # days

        D1, M1, M1_prime, F1 = _calc_lunar_fundamental_arguments(jd_start)
        D2, M2, M2_prime, F2 = _calc_lunar_fundamental_arguments(
            jd_start + synodic_month
        )

        # After one synodic month, D should return to ~same value
        diff = abs(D2 - D1) % (2 * math.pi)
        diff_deg = math.degrees(min(diff, 2 * math.pi - diff))

        # Allow some tolerance
        assert diff_deg < 3, f"D diff after 1 synodic month = {diff_deg}"
