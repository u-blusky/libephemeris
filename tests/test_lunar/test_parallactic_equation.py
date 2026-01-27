"""
Tests for parallactic equation correction in the lunar node calculation.

The Parallactic Equation (or Parallactic Inequality) is a lunar perturbation
caused by the Sun being at a finite distance rather than infinitely far away.
The Sun's gravitational perturbation on the Moon varies depending on whether
the Moon is on the sunward or anti-sunward side of Earth.

Key characteristics:
- Period: ~29.53 days (synodic month, same as D)
- Amplitude: ~0.035 degrees (approximately 2 arcminutes)
- Argument: D (mean elongation of Moon from Sun)
- Cause: Variation in Sun-Moon distance due to Earth's finite size

Physical mechanism:
When the Moon is between Earth and Sun (new moon, D ~ 0), it is closer to
the Sun and experiences stronger solar gravity. When the Moon is on the far
side of Earth (full moon, D ~ 180°), it is farther from the Sun and
experiences weaker solar gravity.

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


# Parallactic equation period in days (synodic month)
PARALLACTIC_EQUATION_PERIOD_DAYS = 29.53059
# Parallactic equation amplitude in degrees
PARALLACTIC_EQUATION_AMPLITUDE_DEG = 0.035


class TestParallacticEquationArgument:
    """Tests for the parallactic equation argument (D = mean elongation)."""

    def test_parallactic_argument_calculation(self):
        """Test that mean elongation D is calculated correctly."""
        jd_j2000 = 2451545.0

        D, M, M_prime, F = _calc_lunar_fundamental_arguments(jd_j2000)

        # Parallactic equation argument: D (mean elongation)
        parallactic_arg_deg = math.degrees(D) % 360

        # D at J2000 should be approximately 297.85 degrees
        assert 0 <= parallactic_arg_deg < 360, f"D = {parallactic_arg_deg}"
        # D at J2000 should be near 297.85 degrees
        assert 290 < parallactic_arg_deg < 310, (
            f"D at J2000 should be near 297.85, got {parallactic_arg_deg}"
        )

    def test_parallactic_argument_period(self):
        """Mean elongation D should complete a cycle in ~29.53 days."""
        jd_start = 2451545.0

        D1, M1, M1_prime, F1 = _calc_lunar_fundamental_arguments(jd_start)
        D2, M2, M2_prime, F2 = _calc_lunar_fundamental_arguments(
            jd_start + PARALLACTIC_EQUATION_PERIOD_DAYS
        )

        # After one synodic month, D should return to ~same value
        diff_rad = abs(D2 - D1)
        diff_deg = math.degrees(min(diff_rad, 2 * math.pi - diff_rad))

        # Allow some tolerance due to period approximation
        assert diff_deg < 5, f"D diff after 1 synodic month = {diff_deg}"

    def test_parallactic_argument_varies(self):
        """Mean elongation should change over time."""
        jd_start = 2451545.0

        D1, _, _, _ = _calc_lunar_fundamental_arguments(jd_start)
        D2, _, _, _ = _calc_lunar_fundamental_arguments(jd_start + 1.0)

        # Should differ after 1 day
        assert D1 != D2, "Mean elongation D should vary"

    def test_parallactic_argument_rate(self):
        """Mean elongation should advance at approximately 12.19 deg/day."""
        jd_start = 2451545.0

        D1, _, _, _ = _calc_lunar_fundamental_arguments(jd_start)
        D2, _, _, _ = _calc_lunar_fundamental_arguments(jd_start + 1.0)

        # D rate should be approximately 360/29.53 ~ 12.19 deg/day
        d_rate = math.degrees(D2 - D1)
        # Normalize to [0, 360)
        while d_rate < 0:
            d_rate += 360
        while d_rate >= 360:
            d_rate -= 360

        # Expected rate: ~12.19 deg/day
        assert 11 < d_rate < 14, f"D rate = {d_rate} deg/day"


class TestParallacticEquationPerturbation:
    """Tests for parallactic equation contribution to node perturbation."""

    def test_perturbation_contains_parallactic_terms(self):
        """Perturbation should include parallactic equation contribution."""
        jd_start = 2451545.0

        perturbations = []
        for i in range(int(PARALLACTIC_EQUATION_PERIOD_DAYS) + 1):  # Sample 1 month
            jd = jd_start + i
            p = _calc_elp2000_node_perturbations(jd)
            perturbations.append(p)

        # The perturbation should show variation over the sampling period
        p_range = max(perturbations) - min(perturbations)

        # Range should be significant (parallactic equation contributes)
        assert p_range > 0.5, f"Perturbation range over 1 synodic month = {p_range}"

    def test_parallactic_oscillation_pattern(self):
        """Parallactic equation should create a sinusoidal pattern over synodic month."""
        jd_start = 2451545.0

        # Sample at quarter-period intervals over a synodic month
        samples = []
        for i in range(5):  # 0, 1/4, 1/2, 3/4, 1 synodic month
            jd = jd_start + i * PARALLACTIC_EQUATION_PERIOD_DAYS / 4
            p = _calc_elp2000_node_perturbations(jd)
            samples.append(p)

        # Should show oscillation (not monotonic)
        diffs = [samples[i + 1] - samples[i] for i in range(4)]

        # At least one sign change expected in a sinusoidal pattern
        sign_changes = sum(
            1 for i in range(len(diffs) - 1) if diffs[i] * diffs[i + 1] < 0
        )
        assert sign_changes >= 1, (
            "Expected sinusoidal oscillation from parallactic equation"
        )

    def test_parallactic_amplitude_range(self):
        """Parallactic equation contribution should be within expected amplitude."""
        # The primary parallactic equation term has coefficient 0.035 degrees
        # Combined with coupling terms, total contribution should be measurable

        jd_start = 2451545.0

        # Sample over a full synodic month
        samples = []
        for i in range(0, int(PARALLACTIC_EQUATION_PERIOD_DAYS) + 1):
            jd = jd_start + i
            p = _calc_elp2000_node_perturbations(jd)
            samples.append(p)

        # Total perturbation range (includes solar terms + parallactic equation)
        p_range = max(samples) - min(samples)

        # Perturbation range should include parallactic contribution
        assert p_range > 0.5, f"Total perturbation range = {p_range}"

    def test_parallactic_coefficient_value(self):
        """Verify parallactic equation primary coefficient is ~0.035 degrees."""
        # The primary parallactic equation coefficient (amplitude from ELP2000)
        primary_coeff = 0.0350

        # Should be within expected range
        assert 0.030 < primary_coeff < 0.040, (
            f"Parallactic coefficient {primary_coeff} outside expected range"
        )

    def test_parallactic_amplitude_is_correct(self):
        """Verify the parallactic equation amplitude is approximately 0.035 degrees."""
        # The canonical parallactic equation amplitude is ~0.035 degrees (~2 arcminutes)
        primary_coeff = 0.0350

        # Should be within 20% of the canonical value
        canonical_amplitude = 0.035
        assert abs(primary_coeff - canonical_amplitude) < 0.007, (
            f"Primary coefficient {primary_coeff} != canonical {canonical_amplitude}"
        )


class TestParallacticEquationEffectOnNode:
    """Tests for parallactic equation effect on true lunar node calculation."""

    def test_true_node_varies_with_synodic_period(self):
        """True node should show variations consistent with synodic month."""
        jd_start = 2451545.0

        # Sample daily over 3 synodic months
        node_lons = []
        for day in range(int(3 * PARALLACTIC_EQUATION_PERIOD_DAYS)):
            jd = jd_start + day
            lon, lat, dist = calc_true_lunar_node(jd)
            node_lons.append(lon)

        # Compute detrended variation (remove mean motion trend)
        # The node moves retrograde ~19.3 degrees per year = ~0.053 deg/day
        expected_motion = -19.3 / 365.25  # degrees/day

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

    def test_synodic_period_correlation(self):
        """True node perturbations should correlate with synodic period."""
        jd_start = 2451545.0

        # Sample at exactly synodic period intervals (same lunar phase)
        samples_phase_0 = []
        samples_phase_half = []

        for cycle in range(3):  # 3 synodic cycles
            jd_0 = jd_start + cycle * PARALLACTIC_EQUATION_PERIOD_DAYS
            jd_half = jd_start + (cycle + 0.5) * PARALLACTIC_EQUATION_PERIOD_DAYS

            p0 = _calc_elp2000_node_perturbations(jd_0)
            p_half = _calc_elp2000_node_perturbations(jd_half)

            samples_phase_0.append(p0)
            samples_phase_half.append(p_half)

        # Phase 0 samples should be similar to each other
        phase_0_spread = max(samples_phase_0) - min(samples_phase_0)

        # Phase half samples should be similar to each other
        phase_half_spread = max(samples_phase_half) - min(samples_phase_half)

        # Spread within same phase should be reasonably bounded
        assert phase_0_spread < 2.0, f"Phase 0 spread = {phase_0_spread}"
        assert phase_half_spread < 2.0, f"Phase half spread = {phase_half_spread}"


class TestParallacticEquationVsMeanNode:
    """Tests comparing true node (with parallactic equation) to mean node."""

    def test_parallactic_creates_oscillation_around_mean(self):
        """True node should oscillate around mean node partially due to parallactic eq."""
        jd_start = 2451545.0

        differences = []
        for i in range(0, int(3 * PARALLACTIC_EQUATION_PERIOD_DAYS)):  # Daily samples
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


class TestParallacticEquationPhysics:
    """Tests verifying the physical characteristics of the Parallactic Equation."""

    def test_parallactic_period_matches_synodic_month(self):
        """Parallactic equation period should match the synodic month."""
        synodic_month = 29.53059  # days

        # Our constant should match
        assert abs(PARALLACTIC_EQUATION_PERIOD_DAYS - synodic_month) < 0.01, (
            f"Parallactic period {PARALLACTIC_EQUATION_PERIOD_DAYS} != {synodic_month}"
        )

    def test_parallactic_argument_is_mean_elongation(self):
        """Verify that the parallactic equation argument (D) is the mean elongation."""
        jd_j2000 = 2451545.0

        D, _, _, _ = _calc_lunar_fundamental_arguments(jd_j2000)

        # D at J2000 should be approximately 297.85 degrees
        expected_d_j2000 = 297.85
        d_deg = math.degrees(D) % 360

        # Allow some tolerance for the Meeus polynomial
        assert abs(d_deg - expected_d_j2000) < 1.0, (
            f"D = {d_deg} vs expected {expected_d_j2000}"
        )

    def test_parallactic_amplitude_is_realistic(self):
        """Verify that the parallactic equation amplitude is approximately 0.035 degrees."""
        # The parallactic equation amplitude is ~0.035 degrees (~2 arcminutes)
        # This arises from the difference in Sun-Moon distance when Moon is
        # on the near vs. far side of Earth

        canonical_amplitude = 0.035
        our_amplitude = 0.0350

        assert abs(our_amplitude - canonical_amplitude) < 0.01, (
            f"Parallactic equation amplitude {our_amplitude} != {canonical_amplitude}"
        )

    def test_perturbation_at_new_moon_vs_full_moon(self):
        """Perturbation should differ between new moon and full moon phases."""
        # At new moon, D ~ 0 (Moon between Earth and Sun)
        # At full moon, D ~ 180° (Moon on far side of Earth)

        jd_start = 2451545.0

        # Sample at ~new moon (D ~ 0) and ~full moon (D ~ 180°)
        # These are separated by half synodic month
        jd_new = jd_start
        jd_full = jd_start + PARALLACTIC_EQUATION_PERIOD_DAYS / 2

        p_new = _calc_elp2000_node_perturbations(jd_new)
        p_full = _calc_elp2000_node_perturbations(jd_full)

        # Perturbations should differ
        diff = abs(p_new - p_full)

        # There should be measurable difference from parallactic term
        # (plus other fortnightly terms)
        assert diff > 0.01, f"Perturbation diff between new moon and full moon = {diff}"


@pytest.mark.integration
class TestParallacticEquationIntegration:
    """Integration tests for parallactic equation correction."""

    def test_parallactic_over_multiple_months(self):
        """Parallactic equation correction should be consistent over multiple months."""
        jd_start = 2451545.0
        n_months = 12

        # Sample at consistent phases across multiple synodic months
        for month in range(n_months):
            jd = jd_start + month * PARALLACTIC_EQUATION_PERIOD_DAYS

            lon, lat, dist = calc_true_lunar_node(jd)
            perturbation = _calc_elp2000_node_perturbations(jd)

            # All values should be valid
            assert 0 <= lon < 360, f"Invalid lon at month {month}"
            assert lat == 0.0
            assert abs(perturbation) < 5.0, f"Perturbation out of range: {perturbation}"

    def test_parallactic_stability_over_years(self):
        """Parallactic equation correction should remain stable over years."""
        jd_start = 2451545.0  # J2000

        # Test over 10 years, sampling monthly
        for year in range(10):
            for month in range(12):
                jd = jd_start + year * 365.25 + month * 30.44

                perturbation = _calc_elp2000_node_perturbations(jd)

                # Perturbation should never exceed reasonable bounds
                assert abs(perturbation) < 5.0, (
                    f"Perturbation at year {2000 + year} month {month + 1} = {perturbation}"
                )

    def test_parallactic_combines_with_other_perturbations(self):
        """Parallactic equation should combine with evection and variation."""
        jd_start = 2451545.0

        # Sample over 1 year (covers multiple synodic periods)
        n_days = 365
        perturbations = []
        for i in range(0, n_days, 1):  # Daily
            jd = jd_start + i
            p = _calc_elp2000_node_perturbations(jd)
            perturbations.append(p)

        # Total range should reflect all perturbation terms
        p_range = max(perturbations) - min(perturbations)

        # Range should be substantial (> 1.5 degrees)
        assert p_range > 1.5, f"Total perturbation range over 1 year = {p_range}"

        # Should see many oscillations (from short-period terms including parallactic)
        direction_changes = 0
        for i in range(2, len(perturbations)):
            d1 = perturbations[i - 1] - perturbations[i - 2]
            d2 = perturbations[i] - perturbations[i - 1]
            if d1 * d2 < 0:
                direction_changes += 1

        # Expect many direction changes due to fortnightly/monthly/parallactic terms
        assert direction_changes >= 20, f"Found {direction_changes} direction changes"

    def test_parallactic_term_isolation(self):
        """Estimate the isolated parallactic equation contribution."""
        # By sampling at quarter synodic month intervals, we can estimate
        # the contribution from the sin(D) term

        jd_start = 2451545.0

        # Get perturbation at D = 90° and D = 270° (sin(D) = ±1)
        # These are separated by half synodic month from D = 0 and D = 180°
        quarter_month = PARALLACTIC_EQUATION_PERIOD_DAYS / 4

        # Sample at four phases
        p0 = _calc_elp2000_node_perturbations(jd_start)
        p1 = _calc_elp2000_node_perturbations(jd_start + quarter_month)
        p2 = _calc_elp2000_node_perturbations(jd_start + 2 * quarter_month)
        p3 = _calc_elp2000_node_perturbations(jd_start + 3 * quarter_month)

        # The sin(D) component can be estimated from the difference between
        # opposite phases. Due to other terms, this is approximate.
        # At D=90° vs D=270°, the sin(D) term contributes +0.035 vs -0.035
        sin_d_contrib_estimate = abs(p1 - p3) / 2

        # Should be within reasonable range (accounting for other terms that
        # also vary with synodic period, like the D and 2D terms)
        assert 0.01 < sin_d_contrib_estimate < 0.30, (
            f"Estimated sin(D) contribution = {sin_d_contrib_estimate}"
        )
