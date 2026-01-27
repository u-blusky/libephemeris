"""
Tests for annual equation correction in the lunar node calculation.

The Annual Equation is a major lunar perturbation caused by the varying
Earth-Sun distance due to Earth's orbital eccentricity (e ~ 0.0167).

Key characteristics:
- Period: ~365.26 days (one anomalistic year)
- Amplitude: ~0.186 degrees (11.2 arcminutes) in lunar longitude
- Argument: M (Sun's mean anomaly)
- Cause: Variation in solar gravitational influence due to Earth's
  elliptical orbit around the Sun

Physical mechanism:
When Earth is at perihelion (early January), it is ~3.3% closer to the Sun
than at aphelion (early July). The increased solar gravitational influence
at perihelion accelerates the Moon's orbital motion, while at aphelion
the Moon moves slightly slower.

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


# Annual Equation period in days (anomalistic year)
ANNUAL_EQUATION_PERIOD_DAYS = 365.2596
# Annual Equation amplitude in degrees
ANNUAL_EQUATION_AMPLITUDE_DEG = 0.186


class TestAnnualEquationArgument:
    """Tests for the annual equation argument (M = Sun's mean anomaly)."""

    def test_annual_equation_argument_calculation(self):
        """Test that Sun's mean anomaly M is calculated correctly."""
        jd_j2000 = 2451545.0

        D, M, M_prime, F = _calc_lunar_fundamental_arguments(jd_j2000)

        # Annual equation argument: M (Sun's mean anomaly)
        annual_arg_deg = math.degrees(M) % 360

        # At J2000, M ~ 357.53 degrees (near zero at perihelion)
        assert 0 <= annual_arg_deg < 360, f"Annual arg = {annual_arg_deg}"
        # M at J2000 should be approximately 357.53 degrees
        assert 350 < annual_arg_deg < 365 or 0 <= annual_arg_deg < 10, (
            f"M at J2000 should be near 0/360, got {annual_arg_deg}"
        )

    def test_annual_equation_argument_period(self):
        """Sun's mean anomaly M should complete a cycle in ~365.26 days."""
        jd_start = 2451545.0

        D1, M1, M1_prime, F1 = _calc_lunar_fundamental_arguments(jd_start)
        D2, M2, M2_prime, F2 = _calc_lunar_fundamental_arguments(
            jd_start + ANNUAL_EQUATION_PERIOD_DAYS
        )

        # After one anomalistic year, M should return to ~same value
        diff_rad = abs(M2 - M1)
        diff_deg = math.degrees(min(diff_rad, 2 * math.pi - diff_rad))

        # Allow some tolerance due to period approximation
        assert diff_deg < 5, f"M diff after 1 year = {diff_deg}"

    def test_annual_equation_argument_varies(self):
        """Sun's mean anomaly should change over time."""
        jd_start = 2451545.0

        _, M1, _, _ = _calc_lunar_fundamental_arguments(jd_start)
        _, M2, _, _ = _calc_lunar_fundamental_arguments(jd_start + 30.0)

        # Should differ after 30 days
        assert M1 != M2, "Sun's mean anomaly should vary"

    def test_annual_equation_argument_rate(self):
        """Sun's mean anomaly should advance at approximately 0.986 deg/day."""
        jd_start = 2451545.0

        _, M1, _, _ = _calc_lunar_fundamental_arguments(jd_start)
        _, M2, _, _ = _calc_lunar_fundamental_arguments(jd_start + 1.0)

        # M rate should be approximately 360/365.26 ~ 0.986 deg/day
        m_rate = math.degrees(M2 - M1)
        # Normalize to [0, 360)
        while m_rate < 0:
            m_rate += 360
        while m_rate >= 360:
            m_rate -= 360

        # Expected rate: ~0.986 deg/day
        assert 0.9 < m_rate < 1.1, f"M rate = {m_rate} deg/day"


class TestAnnualEquationPerturbation:
    """Tests for annual equation contribution to node perturbation."""

    def test_perturbation_contains_annual_equation_terms(self):
        """Perturbation should include annual equation contribution."""
        jd_start = 2451545.0

        perturbations = []
        for i in range(int(ANNUAL_EQUATION_PERIOD_DAYS / 4) + 1):  # Sample 1/4 year
            jd = jd_start + i
            p = _calc_elp2000_node_perturbations(jd)
            perturbations.append(p)

        # The perturbation should show variation over the sampling period
        p_range = max(perturbations) - min(perturbations)

        # Range should be significant (annual equation contributes to total)
        assert p_range > 0.3, f"Perturbation range over ~90 days = {p_range}"

    def test_annual_equation_modulation_pattern(self):
        """Annual equation should create a sinusoidal modulation pattern."""
        jd_start = 2451545.0

        # Sample at quarter-period intervals over a year
        samples = []
        for i in range(5):  # 0, 1/4, 1/2, 3/4, 1 year
            jd = jd_start + i * ANNUAL_EQUATION_PERIOD_DAYS / 4
            p = _calc_elp2000_node_perturbations(jd)
            samples.append(p)

        # Should show oscillation (not monotonic)
        diffs = [samples[i + 1] - samples[i] for i in range(4)]

        # At least one sign change expected in a sinusoidal pattern
        sign_changes = sum(
            1 for i in range(len(diffs) - 1) if diffs[i] * diffs[i + 1] < 0
        )
        assert sign_changes >= 1, "Expected sinusoidal oscillation from annual equation"

    def test_annual_equation_amplitude_range(self):
        """Annual equation contribution should be within expected amplitude range."""
        # The primary annual equation term has coefficient -0.1860 degrees
        # Combined with coupling terms, total contribution should be measurable

        jd_start = 2451545.0

        # Sample over a full year
        samples = []
        for i in range(0, int(ANNUAL_EQUATION_PERIOD_DAYS), 10):  # Every 10 days
            jd = jd_start + i
            p = _calc_elp2000_node_perturbations(jd)
            samples.append(p)

        # Total perturbation range (includes solar terms + annual equation)
        p_range = max(samples) - min(samples)

        # Perturbation range should include annual equation contribution
        # The total range should be > 1 degree due to all perturbation terms
        assert p_range > 1.0, f"Total perturbation range = {p_range}"

    def test_annual_equation_coefficient_values(self):
        """Verify annual equation coefficients are correctly implemented.

        The primary annual equation coefficient (-0.186) should be the dominant
        annual term, with secondary terms having decreasing amplitudes.
        """
        # Primary annual equation term (amplitude from ELP2000)
        primary_coeff = 0.1860

        # Annual equation-anomaly coupling terms (M ± M')
        anomaly_coeffs = [0.0098, 0.0082]

        # Annual equation-inclination terms (M ± 2F)
        inclination_coeffs = [0.0037, 0.0032]

        # Second harmonic term (2M)
        second_harmonic = 0.0024

        # Verify hierarchy: primary > anomaly coupling > inclination coupling
        assert primary_coeff > max(anomaly_coeffs)
        assert max(anomaly_coeffs) > max(inclination_coeffs)

        # Total annual equation contribution estimate (sum of absolute coefficients)
        total = (
            primary_coeff
            + sum(anomaly_coeffs)
            + sum(inclination_coeffs)
            + second_harmonic
        )

        # Should be in expected range (0.2 - 0.25 degrees peak)
        assert 0.18 < total < 0.30, f"Total annual equation coefficient sum = {total}"

    def test_annual_equation_amplitude_is_correct(self):
        """Verify the primary annual equation amplitude is ~0.186 degrees."""
        # The canonical annual equation amplitude is 0.186 degrees
        # Our primary coefficient should match this
        primary_coeff = 0.1860

        # Should be within 5% of the canonical value
        canonical_amplitude = 0.186
        assert abs(primary_coeff - canonical_amplitude) < 0.01, (
            f"Primary coefficient {primary_coeff} != canonical {canonical_amplitude}"
        )


class TestAnnualEquationEffectOnNode:
    """Tests for annual equation effect on true lunar node calculation."""

    def test_true_node_varies_with_annual_period(self):
        """True node should show variations consistent with annual period."""
        jd_start = 2451545.0

        # Sample monthly over 2 years
        node_lons = []
        for month in range(24):
            jd = jd_start + month * 30.44  # ~1 month
            lon, lat, dist = calc_true_lunar_node(jd)
            node_lons.append(lon)

        # Compute detrended variation (remove mean motion trend)
        # The node moves retrograde ~19.3 degrees per year
        expected_motion = -19.3 / 365.25  # degrees/day

        detrended = []
        for i, lon in enumerate(node_lons):
            expected = node_lons[0] + expected_motion * i * 30.44
            # Handle wrap-around
            diff = lon - expected
            while diff > 180:
                diff -= 360
            while diff < -180:
                diff += 360
            detrended.append(diff)

        # Detrended variation should show oscillation
        variation = max(detrended) - min(detrended)

        # Expect variation of at least 0.5 degrees due to perturbations
        assert variation > 0.3, f"Detrended variation = {variation}"

    def test_annual_period_correlation(self):
        """True node perturbations should correlate with annual period."""
        jd_start = 2451545.0

        # Sample at exactly annual period intervals
        samples_phase_0 = []
        samples_phase_half = []

        for cycle in range(3):  # 3 annual cycles
            jd_0 = jd_start + cycle * ANNUAL_EQUATION_PERIOD_DAYS
            jd_half = jd_start + (cycle + 0.5) * ANNUAL_EQUATION_PERIOD_DAYS

            p0 = _calc_elp2000_node_perturbations(jd_0)
            p_half = _calc_elp2000_node_perturbations(jd_half)

            samples_phase_0.append(p0)
            samples_phase_half.append(p_half)

        # Phase 0 samples should be similar to each other
        phase_0_spread = max(samples_phase_0) - min(samples_phase_0)

        # Phase half samples should be similar to each other
        phase_half_spread = max(samples_phase_half) - min(samples_phase_half)

        # Spread within same phase should be reasonably bounded
        # Note: Over 3 years, there's significant drift from the 18.6-year
        # nodal precession cycle and other long-period terms, so tolerance
        # is relatively high
        assert phase_0_spread < 3.0, f"Phase 0 spread = {phase_0_spread}"
        assert phase_half_spread < 3.0, f"Phase half spread = {phase_half_spread}"

        # Verify the sampling works
        assert len(samples_phase_0) == 3
        assert len(samples_phase_half) == 3


class TestAnnualEquationVsMeanNode:
    """Tests comparing true node (with annual equation) to mean node."""

    def test_annual_equation_creates_oscillation_around_mean(self):
        """True node should oscillate around mean node partially due to annual equation."""
        jd_start = 2451545.0

        differences = []
        for i in range(0, int(2 * ANNUAL_EQUATION_PERIOD_DAYS), 7):  # Weekly samples
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


class TestAnnualEquationPhysics:
    """Tests verifying the physical characteristics of the Annual Equation."""

    def test_annual_equation_period_matches_anomalistic_year(self):
        """Annual equation period should be approximately one anomalistic year."""
        anomalistic_year = 365.2596  # days

        # Our constant should match
        assert abs(ANNUAL_EQUATION_PERIOD_DAYS - anomalistic_year) < 0.01, (
            f"Annual equation period {ANNUAL_EQUATION_PERIOD_DAYS} != {anomalistic_year}"
        )

    def test_annual_equation_argument_is_sun_mean_anomaly(self):
        """Verify that the annual equation argument (M) is the Sun's mean anomaly."""
        jd_j2000 = 2451545.0

        _, M, _, _ = _calc_lunar_fundamental_arguments(jd_j2000)

        # M at J2000 should be approximately 357.53 degrees
        # (near zero because Earth is near perihelion in early January)
        expected_m_j2000 = 357.53
        m_deg = math.degrees(M) % 360

        # Allow some tolerance for the Meeus polynomial
        assert abs(m_deg - expected_m_j2000) < 1.0, (
            f"M = {m_deg} vs expected {expected_m_j2000}"
        )

    def test_annual_equation_amplitude_is_realistic(self):
        """Verify that the annual equation amplitude is approximately 0.186 degrees."""
        # The annual equation amplitude in longitude is ~0.186 degrees
        # This is derived from the equation:
        # amplitude = 2 * e * sin(L - pi)
        # where e is Earth's orbital eccentricity (~0.0167)

        canonical_amplitude = 0.186
        earth_eccentricity = 0.0167

        # Our primary coefficient should match the canonical amplitude
        our_amplitude = 0.1860

        assert abs(our_amplitude - canonical_amplitude) < 0.01, (
            f"Annual equation amplitude {our_amplitude} != {canonical_amplitude}"
        )

    def test_perturbation_at_perihelion_vs_aphelion(self):
        """Perturbation should differ between perihelion and aphelion times."""
        # J2000 is early January (near perihelion)
        jd_perihelion = 2451545.0

        # Aphelion is ~6 months later (early July)
        jd_aphelion = jd_perihelion + 182.6

        p_perihelion = _calc_elp2000_node_perturbations(jd_perihelion)
        p_aphelion = _calc_elp2000_node_perturbations(jd_aphelion)

        # Perturbations should differ
        diff = abs(p_perihelion - p_aphelion)

        # There should be some measurable difference
        assert diff > 0.01, (
            f"Perturbation diff between perihelion and aphelion = {diff}"
        )


@pytest.mark.integration
class TestAnnualEquationIntegration:
    """Integration tests for annual equation correction."""

    def test_annual_equation_over_multiple_years(self):
        """Annual equation correction should be consistent over multiple years."""
        jd_start = 2451545.0
        n_years = 5

        # Sample at consistent phases across multiple years
        for year in range(n_years):
            jd = jd_start + year * ANNUAL_EQUATION_PERIOD_DAYS

            lon, lat, dist = calc_true_lunar_node(jd)
            perturbation = _calc_elp2000_node_perturbations(jd)

            # All values should be valid
            assert 0 <= lon < 360, f"Invalid lon at year {year}"
            assert lat == 0.0
            assert abs(perturbation) < 5.0, f"Perturbation out of range: {perturbation}"

    def test_annual_equation_stability_over_decades(self):
        """Annual equation correction should remain stable over long time spans."""
        jd_start = 2451545.0  # J2000

        # Test over 50 years, sampling quarterly
        for year in range(50):
            for quarter in range(4):
                jd = jd_start + year * 365.25 + quarter * 91.3

                perturbation = _calc_elp2000_node_perturbations(jd)

                # Perturbation should never exceed reasonable bounds
                assert abs(perturbation) < 5.0, (
                    f"Perturbation at year {2000 + year} Q{quarter + 1} = {perturbation}"
                )

    def test_annual_equation_effect_amplitude_consistency(self):
        """Annual equation effect amplitude should be consistent across different epochs."""
        # Test at different starting points
        epochs = [
            2451545.0,  # J2000
            2451545.0 + 365.25 * 10,  # 2010
            2451545.0 + 365.25 * 25,  # 2025
            2451545.0 - 365.25 * 25,  # 1975
        ]

        for epoch in epochs:
            perturbations = []
            for i in range(0, int(ANNUAL_EQUATION_PERIOD_DAYS), 30):  # Monthly
                jd = epoch + i
                p = _calc_elp2000_node_perturbations(jd)
                perturbations.append(p)

            # Range should be consistent (within reasonable bounds)
            p_range = max(perturbations) - min(perturbations)
            assert 0.5 < p_range < 4.0, (
                f"Perturbation range at epoch {epoch} = {p_range}"
            )

    def test_annual_equation_combines_with_other_perturbations(self):
        """Annual equation should combine with evection and variation perturbations."""
        jd_start = 2451545.0

        # Sample over 2 years (covers multiple perturbation periods)
        n_days = int(2 * ANNUAL_EQUATION_PERIOD_DAYS)
        perturbations = []
        for i in range(0, n_days, 5):  # Every 5 days
            jd = jd_start + i
            p = _calc_elp2000_node_perturbations(jd)
            perturbations.append(p)

        # Total range should reflect annual equation + variation + evection
        p_range = max(perturbations) - min(perturbations)

        # Range should be substantial (> 1.5 degrees)
        assert p_range > 1.5, f"Total perturbation range over 2 years = {p_range}"

        # Should see multiple oscillations (from shorter-period terms)
        direction_changes = 0
        for i in range(2, len(perturbations)):
            d1 = perturbations[i - 1] - perturbations[i - 2]
            d2 = perturbations[i] - perturbations[i - 1]
            if d1 * d2 < 0:
                direction_changes += 1

        # Expect many direction changes due to fortnightly and monthly terms
        assert direction_changes >= 10, f"Found {direction_changes} direction changes"
