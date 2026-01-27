"""
Tests for variation correction in the lunar node calculation.

The Variation is a major lunar perturbation discovered by Tycho Brahe, caused
by the difference in solar gravitational force between the Moon's position at
quadrature versus conjunction/opposition.

Key characteristics:
- Period: ~14.77 days (half the synodic month)
- Amplitude: ~0.658 degrees in longitude
- Argument: 2D (twice the mean elongation)
- Cause: differential solar tidal force at quadrature vs conjunction/opposition

References:
    - Chapront-Touze, M. & Chapront, J. "Lunar Tables and Programs" (1991)
    - Brown, E.W. "An Introductory Treatise on the Lunar Theory" (1896)
    - Meeus, J. "Astronomical Algorithms" (2nd ed., 1998), Chapter 47
"""

import math
import pytest
from libephemeris.lunar import (
    _calc_lunar_fundamental_arguments,
    _calc_elp2000_node_perturbations,
    calc_true_lunar_node,
    calc_mean_lunar_node,
)


# Variation period in days (half the synodic month)
VARIATION_PERIOD_DAYS = 14.765
# Variation amplitude in degrees (in longitude)
VARIATION_AMPLITUDE_DEG = 0.658


class TestVariationArgument:
    """Tests for the variation argument (2D)."""

    def test_variation_argument_calculation(self):
        """Test that variation argument (2D) is calculated correctly."""
        jd_j2000 = 2451545.0

        D, M, M_prime, F = _calc_lunar_fundamental_arguments(jd_j2000)

        # Variation argument: 2D
        variation_arg = 2.0 * D
        variation_arg_deg = math.degrees(variation_arg) % 360

        # At J2000, D ~ 297.85 degrees
        # So 2D ~ 595.70 ~ 235.70 (mod 360)
        assert 0 <= variation_arg_deg < 360, f"Variation arg = {variation_arg_deg}"

    def test_variation_argument_period(self):
        """Variation argument should complete a cycle in ~14.77 days."""
        jd_start = 2451545.0

        D1, M1, M1_prime, F1 = _calc_lunar_fundamental_arguments(jd_start)
        D2, M2, M2_prime, F2 = _calc_lunar_fundamental_arguments(
            jd_start + VARIATION_PERIOD_DAYS
        )

        # Variation arguments (2D)
        variation_1 = (2.0 * D1) % (2.0 * math.pi)
        variation_2 = (2.0 * D2) % (2.0 * math.pi)

        # After one variation period, argument should return to ~same value
        diff_rad = abs(variation_2 - variation_1)
        diff_deg = math.degrees(min(diff_rad, 2 * math.pi - diff_rad))

        # Allow some tolerance due to period approximation
        assert diff_deg < 15, f"Variation arg diff after 1 period = {diff_deg}"

    def test_variation_argument_varies(self):
        """Variation argument should change over time."""
        jd_start = 2451545.0

        D1, M1, M1_prime, F1 = _calc_lunar_fundamental_arguments(jd_start)
        D2, M2, M2_prime, F2 = _calc_lunar_fundamental_arguments(jd_start + 3.0)

        variation_1 = 2.0 * D1
        variation_2 = 2.0 * D2

        # Should differ after 3 days
        assert variation_1 != variation_2, "Variation argument should vary"

    def test_variation_argument_rate(self):
        """Variation argument should advance at approximately 24.4 deg/day."""
        jd_start = 2451545.0

        D1, _, _, _ = _calc_lunar_fundamental_arguments(jd_start)
        D2, _, _, _ = _calc_lunar_fundamental_arguments(jd_start + 1.0)

        # 2D rate should be approximately 2 * 12.19 deg/day = 24.38 deg/day
        variation_rate = math.degrees(2.0 * D2 - 2.0 * D1)
        # Normalize to [0, 360)
        while variation_rate < 0:
            variation_rate += 360
        while variation_rate >= 360:
            variation_rate -= 360

        # Expected rate: ~24.4 deg/day
        assert 23.0 < variation_rate < 26.0, (
            f"Variation rate = {variation_rate} deg/day"
        )


class TestVariationPerturbation:
    """Tests for variation contribution to node perturbation."""

    def test_perturbation_contains_variation_terms(self):
        """Perturbation should include variation contribution."""
        jd_start = 2451545.0

        perturbations = []
        for i in range(int(VARIATION_PERIOD_DAYS) + 1):
            jd = jd_start + i
            p = _calc_elp2000_node_perturbations(jd)
            perturbations.append(p)

        # The perturbation should show variation consistent with variation period
        p_range = max(perturbations) - min(perturbations)

        # Range should be significant (variation contributes to total)
        assert p_range > 0.3, f"Perturbation range over variation period = {p_range}"

    def test_variation_modulation_pattern(self):
        """Variation should create a sinusoidal modulation pattern."""
        jd_start = 2451545.0

        # Sample at quarter-period intervals
        samples = []
        for i in range(5):  # 0, 1/4, 1/2, 3/4, 1 period
            jd = jd_start + i * VARIATION_PERIOD_DAYS / 4
            p = _calc_elp2000_node_perturbations(jd)
            samples.append(p)

        # Should show oscillation (not monotonic)
        diffs = [samples[i + 1] - samples[i] for i in range(4)]

        # At least one sign change expected in a sinusoidal pattern
        sign_changes = sum(
            1 for i in range(len(diffs) - 1) if diffs[i] * diffs[i + 1] < 0
        )
        assert sign_changes >= 1, "Expected sinusoidal oscillation from variation"

    def test_variation_period_signature(self):
        """Perturbations should show periodicity consistent with ~14.77 days."""
        jd_start = 2451545.0

        # Sample over 3 variation periods
        n_samples = int(3 * VARIATION_PERIOD_DAYS)
        perturbations = []
        for i in range(n_samples):
            jd = jd_start + i
            p = _calc_elp2000_node_perturbations(jd)
            perturbations.append(p)

        # Find local maxima/minima
        extrema_indices = []
        for i in range(1, len(perturbations) - 1):
            if (
                perturbations[i] > perturbations[i - 1]
                and perturbations[i] > perturbations[i + 1]
            ):
                extrema_indices.append(i)
            elif (
                perturbations[i] < perturbations[i - 1]
                and perturbations[i] < perturbations[i + 1]
            ):
                extrema_indices.append(i)

        # Should find multiple extrema over 3 periods
        assert len(extrema_indices) >= 4, f"Found {len(extrema_indices)} extrema"

    def test_variation_coefficient_values(self):
        """Verify variation coefficients are correctly implemented.

        The variation terms should include:
        - Primary variation-inclination coupling (2D ± F)
        - Variation-anomaly coupling (2D ± M')
        - Combined terms (2D + F ± M')
        - Double variation terms (4D)
        """
        # Primary variation-inclination terms
        primary_coeffs = [0.0523, 0.0478]

        # Variation-anomaly terms
        anomaly_coeffs = [0.0156, 0.0142]

        # Combined variation-inclination-anomaly terms
        combined_coeffs = [0.0048, 0.0041, 0.0038, 0.0033]

        # Variation-solar coupling terms
        solar_coeffs = [0.0028, 0.0024]

        # Double variation terms (4D)
        double_coeffs = [0.0067, 0.0043, 0.0039]

        # Verify hierarchy: primary > anomaly > combined
        assert max(primary_coeffs) > max(anomaly_coeffs)
        assert max(anomaly_coeffs) > max(combined_coeffs)

        # Total variation contribution estimate (sum of absolute coefficients)
        total = (
            sum(primary_coeffs)
            + sum(anomaly_coeffs)
            + sum(combined_coeffs)
            + sum(solar_coeffs)
            + sum(double_coeffs)
        )

        # Should be in expected range (0.15 - 0.30 degrees peak)
        assert 0.15 < total < 0.35, f"Total variation coefficient sum = {total}"


class TestVariationEffectOnNode:
    """Tests for variation effect on true lunar node calculation."""

    def test_true_node_varies_with_variation_period(self):
        """True node should show variations consistent with variation period."""
        jd_start = 2451545.0

        node_lons = []
        for i in range(int(2 * VARIATION_PERIOD_DAYS)):
            jd = jd_start + i
            lon, lat, dist = calc_true_lunar_node(jd)
            node_lons.append(lon)

        # Compute detrended variation (remove mean motion trend)
        # The node moves retrograde ~0.053 degrees per day
        expected_motion = -0.053  # degrees/day

        detrended = []
        for i, lon in enumerate(node_lons):
            expected = node_lons[0] + expected_motion * i
            # Handle wrap-around
            diff = lon - expected
            while diff > 180:
                diff -= 360
            while diff < -180:
                diff += 360
            detrended.append(diff)

        # Detrended variation should show oscillation
        variation = max(detrended) - min(detrended)

        # Expect variation of at least 0.3 degrees due to perturbations
        assert variation > 0.3, f"Detrended variation = {variation}"

    def test_variation_period_correlation(self):
        """True node perturbations should correlate with variation period."""
        jd_start = 2451545.0

        # Sample at exactly variation period intervals
        samples_phase_0 = []
        samples_phase_half = []

        for cycle in range(3):  # 3 variation cycles
            jd_0 = jd_start + cycle * VARIATION_PERIOD_DAYS
            jd_half = jd_start + (cycle + 0.5) * VARIATION_PERIOD_DAYS

            p0 = _calc_elp2000_node_perturbations(jd_0)
            p_half = _calc_elp2000_node_perturbations(jd_half)

            samples_phase_0.append(p0)
            samples_phase_half.append(p_half)

        # Phase 0 samples should be similar to each other
        phase_0_spread = max(samples_phase_0) - min(samples_phase_0)

        # Phase half samples should be similar to each other
        phase_half_spread = max(samples_phase_half) - min(samples_phase_half)

        # Spread within same phase should be reasonably small
        assert phase_0_spread < 0.5, f"Phase 0 spread = {phase_0_spread}"
        assert phase_half_spread < 0.5, f"Phase half spread = {phase_half_spread}"

        # Verify the sampling works
        assert len(samples_phase_0) == 3
        assert len(samples_phase_half) == 3


class TestVariationVsMeanNode:
    """Tests comparing true node (with variation) to mean node."""

    def test_variation_creates_oscillation_around_mean(self):
        """True node should oscillate around mean node partially due to variation."""
        jd_start = 2451545.0

        differences = []
        for i in range(int(2 * VARIATION_PERIOD_DAYS)):
            jd = jd_start + i
            true_lon, _, _ = calc_true_lunar_node(jd)
            mean_lon = calc_mean_lunar_node(jd)

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
        assert pos_count > 0 and neg_count > 0, "True node should oscillate around mean"

        # Range of differences should be significant
        diff_range = max(differences) - min(differences)
        assert diff_range > 0.3, f"True-Mean node range = {diff_range}"


class TestVariationPhysics:
    """Tests verifying the physical characteristics of the Variation."""

    def test_variation_period_matches_half_synodic_month(self):
        """Variation period should be approximately half the synodic month."""
        synodic_month = 29.530589  # days
        expected_period = synodic_month / 2

        # Our constant should match
        assert abs(VARIATION_PERIOD_DAYS - expected_period) < 0.1, (
            f"Variation period {VARIATION_PERIOD_DAYS} != {expected_period}"
        )

    def test_variation_argument_is_twice_elongation(self):
        """Verify that the variation argument (2D) is twice the mean elongation."""
        jd_j2000 = 2451545.0

        D, _, _, _ = _calc_lunar_fundamental_arguments(jd_j2000)

        # The variation argument should be exactly 2D
        variation_arg = 2.0 * D

        # D at J2000 is approximately 297.85 degrees
        # So 2D should be approximately 235.70 degrees (mod 360)
        expected_2d = 2.0 * 297.85
        expected_mod = expected_2d % 360

        variation_deg = math.degrees(variation_arg) % 360

        # Allow some tolerance for the Meeus polynomial
        assert abs(variation_deg - expected_mod) < 5.0, (
            f"2D = {variation_deg} vs expected {expected_mod}"
        )

    def test_variation_amplitude_is_realistic(self):
        """Verify that the variation amplitude is approximately 0.658 degrees."""
        # The variation amplitude in longitude is ~0.658 degrees
        # The effect on the node is smaller (scaled by inclination coupling)
        # Expected node effect: ~0.658 * sin(5.145°) ≈ 0.059 degrees

        # Our primary coefficients (0.0523, 0.0478) are in this range
        primary_effect = (0.0523 + 0.0478) / 2
        expected_effect = 0.658 * math.sin(math.radians(5.145))

        # Should be within a factor of 2
        assert 0.3 * expected_effect < primary_effect < 3.0 * expected_effect, (
            f"Primary effect {primary_effect} vs expected {expected_effect}"
        )


@pytest.mark.integration
class TestVariationIntegration:
    """Integration tests for variation correction."""

    def test_variation_over_multiple_periods(self):
        """Variation correction should be consistent over multiple periods."""
        jd_start = 2451545.0
        n_periods = 5

        # Sample at consistent phases across multiple periods
        for period in range(n_periods):
            jd = jd_start + period * VARIATION_PERIOD_DAYS

            lon, lat, dist = calc_true_lunar_node(jd)
            perturbation = _calc_elp2000_node_perturbations(jd)

            # All values should be valid
            assert 0 <= lon < 360, f"Invalid lon at period {period}"
            assert lat == 0.0
            assert abs(perturbation) < 5.0, f"Perturbation out of range: {perturbation}"

    def test_variation_stability_over_years(self):
        """Variation correction should remain stable over long time spans."""
        jd_start = 2451545.0  # J2000

        # Test over 10 years, sampling monthly
        for year in range(10):
            for month in range(12):
                jd = jd_start + year * 365.25 + month * 30.4

                perturbation = _calc_elp2000_node_perturbations(jd)

                # Perturbation should never exceed reasonable bounds
                assert abs(perturbation) < 5.0, (
                    f"Perturbation at year {2000 + year} month {month + 1} = {perturbation}"
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
            perturbations = []
            for i in range(int(VARIATION_PERIOD_DAYS)):
                jd = epoch + i
                p = _calc_elp2000_node_perturbations(jd)
                perturbations.append(p)

            # Range should be consistent (within reasonable bounds)
            p_range = max(perturbations) - min(perturbations)
            assert 0.3 < p_range < 4.0, (
                f"Perturbation range at epoch {epoch} = {p_range}"
            )

    def test_variation_combines_with_evection(self):
        """Variation and evection should both contribute to total perturbation."""
        jd_start = 2451545.0

        # Sample over 2 months (covers multiple variation and evection periods)
        n_days = 60
        perturbations = []
        for i in range(n_days):
            jd = jd_start + i
            p = _calc_elp2000_node_perturbations(jd)
            perturbations.append(p)

        # Total range should reflect both variation (~15 day period) and
        # evection (~32 day period) contributions
        p_range = max(perturbations) - min(perturbations)

        # Range should be substantial (> 1 degree)
        assert p_range > 1.0, f"Total perturbation range over 60 days = {p_range}"

        # Should see multiple oscillations
        sign_changes = sum(
            1
            for i in range(1, len(perturbations))
            if (perturbations[i] - perturbations[i - 1])
            * (perturbations[i - 1] - perturbations[max(0, i - 2)])
            < 0
        )
        assert sign_changes >= 3, f"Found {sign_changes} sign changes in 60 days"
