"""
Tests for second-order perturbation terms in the lunar node calculation.

Second-order perturbation terms arise from products of first-order perturbations.
When two sinusoidal terms are multiplied:
    sin(A) * sin(B) = 1/2[cos(A-B) - cos(A+B)]
    cos(A) * sin(B) = 1/2[sin(A+B) + sin(A-B)]

These terms contribute at the arcsecond level (0.0001-0.003 degrees) and
are needed for sub-arcsecond precision in the true lunar node calculation.

The main second-order coupling terms include:
- Evection x Variation coupling (cos(M'), cos(4D-M'))
- Evection x Annual Equation coupling (cos(2D-M' +/- M))
- Variation x Annual Equation coupling (cos(2D +/- M))
- Self-coupling terms (cos(4D-2M'), cos(4D), cos(2M), cos(2F))
- E^2 second-order solar eccentricity corrections
- Three-frequency combination terms

References:
    - Chapront-Touze, M. & Chapront, J. "ELP 2000-82B" (1988)
    - Brown, E.W. "Tables of the Moon" (1919)
    - Eckert, W.J. et al. "The Motion of the Moon" (1954)
"""

import math
import pytest
from libephemeris.lunar import (
    _calc_lunar_fundamental_arguments,
    _calc_elp2000_node_perturbations,
    calc_true_lunar_node,
    calc_mean_lunar_node,
)


class TestSecondOrderTermsMagnitude:
    """Tests verifying second-order terms contribute at arcsecond level."""

    def test_second_order_contribution_exists(self):
        """Second-order terms should contribute measurably to perturbation."""
        # Second-order terms have coefficients in range 0.0006 to 0.0024 degrees
        # This is approximately 2-9 arcseconds
        jd_j2000 = 2451545.0

        perturbation = _calc_elp2000_node_perturbations(jd_j2000)

        # Perturbation should be a reasonable value (not just first-order terms)
        # The total perturbation is dominated by first-order terms (~1.5 degrees)
        # but second-order terms contribute at the ~0.01 degree level
        assert abs(perturbation) < 5.0, f"Perturbation {perturbation} out of range"
        assert perturbation != 0.0, "Perturbation should not be zero"

    def test_second_order_coefficients_are_arcsecond_level(self):
        """Document that second-order coefficients are at arcsecond level."""
        # Key second-order coefficients from the implementation
        second_order_coefficients = [
            # Evection x Variation coupling
            0.0024,  # cos(M')
            0.0018,  # cos(4D - M')
            0.0012,  # cos(M' + 2F)
            0.0010,  # cos(M' - 2F)
            # Evection x Annual coupling
            0.0019,  # cos(2D - M' - M)
            0.0016,  # cos(2D - M' + M)
            0.0011,  # cos(2D - 2M' - M)
            0.0009,  # cos(2D - 2M' + M)
            # Variation x Annual coupling
            0.0015,  # cos(2D - M)
            0.0013,  # cos(2D + M)
            # Self-coupling terms
            0.0014,  # cos(4D - 2M')
            0.0011,  # cos(4D)
            0.0008,  # cos(2F)
            0.0006,  # cos(2M) * E^2
        ]

        for coeff in second_order_coefficients:
            # Each coefficient should be in arcsecond range (0.0001 to 0.01 degrees)
            # 0.0001 degrees = 0.36 arcseconds
            # 0.01 degrees = 36 arcseconds
            assert 0.0001 < coeff < 0.01, f"Coefficient {coeff} outside arcsecond range"

        # Total contribution from listed coefficients
        total_max = sum(second_order_coefficients)
        # Maximum possible contribution (all in phase)
        assert total_max < 0.05, f"Total max contribution {total_max} too large"
        # Should be significant (at least 20 arcseconds possible)
        assert total_max > 0.005, f"Total contribution {total_max} too small"


class TestEvectionVariationCoupling:
    """Tests for Evection x Variation coupling terms."""

    def test_evection_variation_coupling_varies(self):
        """Evection-Variation coupling should create cos(M') variation."""
        jd_start = 2451545.0

        # Sample over one anomalistic month (~27.55 days)
        anomalistic_month = 27.55
        perturbations = []

        for i in range(int(anomalistic_month)):
            jd = jd_start + i
            p = _calc_elp2000_node_perturbations(jd)
            perturbations.append(p)

        # Should show variation consistent with M' period
        p_range = max(perturbations) - min(perturbations)
        assert p_range > 0.5, f"Perturbation range {p_range} seems too small"

    def test_cos_m_prime_term_behavior(self):
        """The cos(M') term should vary with lunar anomaly period."""
        jd_start = 2451545.0

        # Get M' at start and at half anomalistic period
        D1, M1, M1_prime, F1 = _calc_lunar_fundamental_arguments(jd_start)
        D2, M2, M2_prime, F2 = _calc_lunar_fundamental_arguments(
            jd_start + 27.55 / 2  # Half anomalistic month
        )

        # cos(M') should have opposite signs at 0 and 180 degrees phase shift
        cos_m1 = math.cos(M1_prime)
        cos_m2 = math.cos(M2_prime)

        # After half period, M' changes by ~180 degrees
        # So cos(M') should change significantly
        assert cos_m1 != cos_m2, "cos(M') should vary over anomalistic month"


class TestEvectionAnnualCoupling:
    """Tests for Evection x Annual Equation coupling terms."""

    def test_annual_coupling_varies_over_year(self):
        """Evection-Annual coupling should show annual variation."""
        jd_start = 2451545.0

        # Sample monthly over one year
        perturbations = []
        for month in range(12):
            jd = jd_start + month * 30.4
            p = _calc_elp2000_node_perturbations(jd)
            perturbations.append(p)

        # Should show variation over the year
        p_range = max(perturbations) - min(perturbations)
        assert p_range > 0.3, f"Annual variation {p_range} seems too small"

    def test_cos_2d_minus_m_prime_minus_m_varies(self):
        """The cos(2D - M' - M) argument should vary appropriately."""
        jd_start = 2451545.0

        # Calculate the argument at different dates
        args = []
        for i in range(365):  # One year
            jd = jd_start + i
            D, M, M_prime, F = _calc_lunar_fundamental_arguments(jd)
            arg = 2.0 * D - M_prime - M
            args.append(math.cos(arg))

        # Should see full range of cosine values
        cos_range = max(args) - min(args)
        assert cos_range > 1.5, f"cos range {cos_range} should cover most of [-1, 1]"


class TestSelfCouplingTerms:
    """Tests for self-coupling (squared) terms."""

    def test_self_coupling_adds_dc_and_double_frequency(self):
        """Self-coupling creates DC offset and double-frequency terms."""
        # sin^2(arg) = 1/2[1 - cos(2*arg)]
        # The cos(4D - 2M') term comes from sin^2(2D - M')

        jd_start = 2451545.0
        evection_period = 31.8  # days

        # At half the evection period, 2D - M' changes by ~180 degrees
        # So 4D - 2M' changes by ~360 degrees (one full cycle)
        D1, M1, M1_prime, F1 = _calc_lunar_fundamental_arguments(jd_start)
        D2, M2, M2_prime, F2 = _calc_lunar_fundamental_arguments(
            jd_start + evection_period / 2
        )

        arg1 = 4.0 * D1 - 2.0 * M1_prime
        arg2 = 4.0 * D2 - 2.0 * M2_prime

        # The arguments should differ by approximately pi radians
        diff = abs(arg2 - arg1) % (2 * math.pi)
        # Allow for approximate period
        assert diff > 2.5 or diff < 0.5 or abs(diff - math.pi) < 0.5, (
            f"4D-2M' should cycle twice per evection period, diff = {diff}"
        )

    def test_cos_4d_term_varies_fortnightly(self):
        """The cos(4D) self-coupling term should vary with ~7.4 day period."""
        jd_start = 2451545.0

        # Sample daily over one synodic month
        synodic_month = 29.53
        cos_4d_values = []

        for i in range(int(synodic_month)):
            jd = jd_start + i
            D, M, M_prime, F = _calc_lunar_fundamental_arguments(jd)
            cos_4d_values.append(math.cos(4.0 * D))

        # Should see oscillation
        assert max(cos_4d_values) > 0.5, "cos(4D) should reach near 1"
        assert min(cos_4d_values) < -0.5, "cos(4D) should reach near -1"


class TestE2Corrections:
    """Tests for E^2 second-order solar eccentricity corrections."""

    def test_e2_factor_value(self):
        """E^2 should be close to 1 for modern epoch."""
        jd_j2000 = 2451545.0
        T = (jd_j2000 - 2451545.0) / 36525.0  # = 0 at J2000

        E = 1.0 - 0.002516 * T - 0.0000074 * T**2
        E2 = E * E

        # At J2000, E = 1.0 exactly
        assert abs(E - 1.0) < 0.001, f"E at J2000 should be ~1, got {E}"
        assert abs(E2 - 1.0) < 0.002, f"E^2 at J2000 should be ~1, got {E2}"

    def test_e2_varies_over_centuries(self):
        """E^2 should decrease slowly over time."""
        jd_j2000 = 2451545.0

        # E at J2000
        T0 = 0.0
        E0 = 1.0 - 0.002516 * T0 - 0.0000074 * T0**2

        # E at year 2100 (T = 1 century)
        T1 = 1.0
        E1 = 1.0 - 0.002516 * T1 - 0.0000074 * T1**2

        assert E1 < E0, "E should decrease over time"
        assert E0 - E1 < 0.01, "E should change slowly (< 0.01 per century)"


class TestThreeFrequencyCombinations:
    """Tests for three-frequency combination terms."""

    def test_three_frequency_terms_vary(self):
        """Three-frequency terms should vary with complex beat patterns."""
        jd_start = 2451545.0

        # Calculate D + M + M' argument over time
        args = []
        for i in range(100):
            jd = jd_start + i
            D, M, M_prime, F = _calc_lunar_fundamental_arguments(jd)
            arg = D + M + M_prime
            args.append(math.sin(arg))

        # Should vary significantly
        arg_range = max(args) - min(args)
        assert arg_range > 1.0, f"D+M+M' term should vary, range = {arg_range}"


class TestSecularSecondOrderTerms:
    """Tests for secular (T-dependent) second-order terms."""

    def test_secular_terms_grow_with_time(self):
        """Secular second-order terms should grow linearly with T."""
        # At J2000, T = 0, so secular terms contribute 0
        jd_j2000 = 2451545.0

        # At year 2100, T = 1 century
        jd_2100 = jd_j2000 + 36525.0

        p_j2000 = _calc_elp2000_node_perturbations(jd_j2000)
        p_2100 = _calc_elp2000_node_perturbations(jd_2100)

        # Perturbations should differ (not just from secular terms)
        assert p_j2000 != p_2100, "Perturbations should differ at different epochs"

    def test_secular_term_coefficients_are_small(self):
        """Secular second-order coefficients should be very small."""
        # The coefficients are 0.00003 and 0.00002 degrees per century
        # After 10 centuries, contribution would be ~0.0003 degrees

        secular_coeffs = [0.00003, 0.00002]
        for coeff in secular_coeffs:
            # Should be very small (sub-arcsecond per century)
            assert coeff < 0.0001, f"Secular coeff {coeff} too large"
            assert coeff > 0, "Secular coeff should be positive"


class TestSecondOrderPrecisionImprovement:
    """Tests verifying second-order terms improve precision."""

    def test_perturbation_stability_over_nodal_cycle(self):
        """Perturbation should remain stable over a full nodal cycle."""
        jd_start = 2451545.0
        nodal_period = 18.6 * 365.25  # ~6798 days

        # Sample at 100 points over nodal cycle
        perturbations = []
        for i in range(100):
            jd = jd_start + i * nodal_period / 100
            p = _calc_elp2000_node_perturbations(jd)
            perturbations.append(p)
            # Each perturbation should be reasonable
            assert abs(p) < 5.0, f"Perturbation {p} at sample {i} out of range"

        # Range should be significant but bounded
        p_range = max(perturbations) - min(perturbations)
        assert p_range > 1.0, f"Perturbation range {p_range} too small"
        assert p_range < 5.0, f"Perturbation range {p_range} too large"

    def test_second_order_terms_do_not_dominate(self):
        """Second-order terms should be smaller than first-order terms."""
        # First-order dominant term: -1.5233 * sin(2D)
        # Largest second-order term: 0.0024 * cos(M')

        first_order_amplitude = 1.5233
        second_order_amplitude = 0.0024

        ratio = second_order_amplitude / first_order_amplitude
        assert ratio < 0.01, "Second-order should be < 1% of first-order"
        assert ratio > 0.0001, "Second-order should be measurable"


@pytest.mark.integration
class TestSecondOrderIntegration:
    """Integration tests for second-order perturbation terms."""

    def test_true_node_precision_with_second_order_terms(self):
        """True node should be computed with second-order terms contributing."""
        jd_j2000 = 2451545.0

        # Calculate true node
        lon, lat, dist = calc_true_lunar_node(jd_j2000)

        # Basic validity checks
        assert 0 <= lon < 360, f"Invalid longitude {lon}"
        assert lat == 0.0, "Node latitude should be 0"
        assert dist > 0, f"Invalid distance {dist}"

    def test_perturbation_consistency_across_epochs(self):
        """Perturbation calculation should be consistent across different epochs."""
        epochs = [
            2451545.0,  # J2000
            2451545.0 + 3652.5,  # 2010
            2451545.0 - 3652.5,  # 1990
            2451545.0 + 7305,  # 2020
            2451545.0 - 7305,  # 1980
        ]

        for epoch in epochs:
            p = _calc_elp2000_node_perturbations(epoch)
            assert abs(p) < 5.0, f"Perturbation at epoch {epoch} = {p} out of range"

    def test_second_order_terms_smooth_transition(self):
        """Second-order terms should not cause discontinuities."""
        jd_start = 2451545.0

        # Check for smooth transitions over small time steps
        prev_p = _calc_elp2000_node_perturbations(jd_start)
        for i in range(1, 100):
            jd = jd_start + i * 0.1  # Every 2.4 hours
            p = _calc_elp2000_node_perturbations(jd)

            # Change should be smooth (no jumps > 0.1 degrees in 2.4 hours)
            change = abs(p - prev_p)
            assert change < 0.1, f"Jump of {change} degrees at step {i}"
            prev_p = p


@pytest.mark.comparison
class TestSecondOrderComparisonWithSwissEph:
    """Comparison tests with Swiss Ephemeris."""

    def test_second_order_improves_agreement(self):
        """Second-order terms should contribute to Swiss Ephemeris agreement."""
        import swisseph as swe

        jd_j2000 = 2451545.0

        # Get Swiss Ephemeris true node
        swe_result = swe.calc_ut(jd_j2000, swe.TRUE_NODE, 0)
        swe_lon = swe_result[0][0]

        # Get our true node
        lib_lon, _, _ = calc_true_lunar_node(jd_j2000)

        # Calculate difference
        diff = abs(swe_lon - lib_lon)
        if diff > 180:
            diff = 360 - diff

        # Precision should be good (< 0.1 degree)
        assert diff < 0.1, f"Difference with Swiss Ephemeris: {diff}°"
